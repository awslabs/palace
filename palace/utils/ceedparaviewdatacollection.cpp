// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/ceedparaviewdatacollection.hpp"

#include <algorithm>
#include <limits>
#include <regex>
#include <type_traits>
#include <utility>
#include <vector>
#include <mfem/mesh/vtk.hpp>

namespace palace
{

namespace
{

void WriteAll(std::ostream &os, const void *data, std::size_t bytes)
{
  const char *ptr = static_cast<const char *>(data);
  constexpr std::size_t max_write = std::size_t{1} << 30;
  while (bytes > 0)
  {
    const std::size_t chunk = std::min(bytes, max_write);
    os.write(ptr, static_cast<std::streamsize>(chunk));
    MFEM_VERIFY(os.good(), "Failed while writing contiguous ParaView point data payload!");
    ptr += chunk;
    bytes -= chunk;
  }
}

}  // namespace

std::map<std::string, CeedParaViewDataCollection::PointField> &
CeedParaViewDataCollection::PointFields(MeshEntityType location)
{
  switch (location)
  {
    case MeshEntityType::Domain:
      return domain_point_fields;
    case MeshEntityType::Boundary:
      return boundary_point_fields;
  }
  MFEM_ABORT("Unknown point field location!");
}

const std::map<std::string, CeedParaViewDataCollection::PointField> &
CeedParaViewDataCollection::PointFields(MeshEntityType location) const
{
  switch (location)
  {
    case MeshEntityType::Domain:
      return domain_point_fields;
    case MeshEntityType::Boundary:
      return boundary_point_fields;
  }
  MFEM_ABORT("Unknown point field location!");
}

int CeedParaViewDataCollection::NumPointFieldEntities(MeshEntityType location) const
{
  return location == MeshEntityType::Domain ? mesh->GetNE() : mesh->GetNBE();
}

mfem::Geometry::Type
CeedParaViewDataCollection::PointFieldBaseGeometry(MeshEntityType location, int i) const
{
  return location == MeshEntityType::Domain ? mesh->GetElementBaseGeometry(i)
                                            : mesh->GetBdrElementBaseGeometry(i);
}

const char *CeedParaViewDataCollection::PointFieldLocationName(MeshEntityType location)
{
  return location == MeshEntityType::Domain ? "Domain" : "Boundary";
}

const char *CeedParaViewDataCollection::PointFieldEntityName(MeshEntityType location)
{
  return location == MeshEntityType::Domain ? "element" : "boundary element";
}

bool CeedParaViewDataCollection::LocationMatchesOutput(MeshEntityType location) const
{
  return location == MeshEntityType::Domain ? !bdr_output : bdr_output;
}

void CeedParaViewDataCollection::RegisterPointField(MeshEntityType location,
                                                    const std::string &field_name,
                                                    const Vector &values,
                                                    const std::vector<int> &bases,
                                                    int num_comp)
{
  MFEM_VERIFY(num_comp > 0, PointFieldLocationName(location)
                                << " point field must have at least one component!");
  PointFields(location)[field_name] =
      PointField{&values, {}, &bases, num_comp, values.Size()};
}

void CeedParaViewDataCollection::RegisterPointEvaluator(
    MeshEntityType location, const std::string &field_name,
    std::function<void(Vector &)> evaluator, const std::vector<int> &bases, int num_comp,
    int buffer_size)
{
  MFEM_VERIFY(evaluator, PointFieldLocationName(location)
                             << " point evaluator must be callable!");
  MFEM_VERIFY(num_comp > 0, PointFieldLocationName(location)
                                << " point field must have at least one component!");
  MFEM_VERIFY(buffer_size >= 0, PointFieldLocationName(location)
                                    << " point evaluator buffer size is invalid!");
  PointFields(location)[field_name] =
      PointField{nullptr, std::move(evaluator), &bases, num_comp, buffer_size};
}

void CeedParaViewDataCollection::DeregisterPointField(MeshEntityType location,
                                                      const std::string &field_name)
{
  PointFields(location).erase(field_name);
}

void CeedParaViewDataCollection::RegisterBoundaryPointField(const std::string &field_name,
                                                            const Vector &values,
                                                            const std::vector<int> &bases,
                                                            int num_comp)
{
  RegisterPointField(MeshEntityType::Boundary, field_name, values, bases, num_comp);
}

void CeedParaViewDataCollection::RegisterBoundaryPointEvaluator(
    const std::string &field_name, std::function<void(Vector &)> evaluator,
    const std::vector<int> &bases, int num_comp, int buffer_size)
{
  RegisterPointEvaluator(MeshEntityType::Boundary, field_name, std::move(evaluator), bases,
                         num_comp, buffer_size);
}

void CeedParaViewDataCollection::RegisterDomainPointEvaluator(
    const std::string &field_name, std::function<void(Vector &)> evaluator,
    const std::vector<int> &bases, int num_comp, int buffer_size)
{
  RegisterPointEvaluator(MeshEntityType::Domain, field_name, std::move(evaluator), bases,
                         num_comp, buffer_size);
}

void CeedParaViewDataCollection::DeregisterBoundaryPointField(const std::string &field_name)
{
  DeregisterPointField(MeshEntityType::Boundary, field_name);
}

void CeedParaViewDataCollection::DeregisterDomainPointField(const std::string &field_name)
{
  DeregisterPointField(MeshEntityType::Domain, field_name);
}

bool CeedParaViewDataCollection::UseAppendedPointFields(MeshEntityType location) const
{
  const auto &fields = PointFields(location);
  return LocationMatchesOutput(location) && !fields.empty() &&
         pv_data_format != mfem::VTKFormat::ASCII && GetCompressionLevel() == 0;
}

int CeedParaViewDataCollection::MaxPointFieldBufferSize(
    const std::map<std::string, PointField> &fields) const
{
  int size = 0;
  for (const auto &[name, field] : fields)
  {
    size = std::max(size, field.buffer_size);
  }
  return size;
}

std::uint64_t
CeedParaViewDataCollection::PointFieldPayloadSize(MeshEntityType location, int ref,
                                                  const PointField &field) const
{
  MFEM_VERIFY((field.values || field.evaluator) && field.bases && field.num_comp > 0,
              "Invalid " << PointFieldLocationName(location)
                         << " point field registration!");
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == NumPointFieldEntities(location),
              PointFieldLocationName(location)
                  << " point field base offsets do not match the mesh!");

  std::uint64_t value_count = 0;
  for (int i = 0; i < NumPointFieldEntities(location); i++)
  {
    const auto *RefG =
        mfem::GlobGeometryRefiner.Refine(PointFieldBaseGeometry(location, i), ref, 1);
    value_count += static_cast<std::uint64_t>(field.num_comp) *
                   static_cast<std::uint64_t>(RefG->RefPts.GetNPoints());
  }
  const std::uint64_t scalar_size =
      pv_data_format == mfem::VTKFormat::BINARY32 ? sizeof(float) : sizeof(double);
  return value_count * scalar_size;
}

void CeedParaViewDataCollection::WritePointFieldValues(MeshEntityType location,
                                                       std::ostream &os, int ref,
                                                       const PointField &field,
                                                       const Vector &values) const
{
  MFEM_VERIFY(values.Size() % field.num_comp == 0,
              PointFieldLocationName(location)
                  << " point field buffer size is not divisible by its component count!");
  const int component_stride = values.Size() / field.num_comp;
  const double *data = values.HostRead();

  const std::uint64_t payload_size = PointFieldPayloadSize(location, ref, field);
  const bool binary32 = pv_data_format == mfem::VTKFormat::BINARY32;
  const std::size_t scalar_size = binary32 ? sizeof(float) : sizeof(double);
  MFEM_VERIFY(payload_size % scalar_size == 0,
              "Point field payload size is not divisible by scalar size!");
  const std::size_t payload_count = static_cast<std::size_t>(payload_size / scalar_size);

  auto PackPayload = [&](auto *payload_data)
  {
    using scalar_t = std::remove_pointer_t<decltype(payload_data)>;
    std::size_t out = 0;
    for (int i = 0; i < NumPointFieldEntities(location); i++)
    {
      const auto *RefG =
          mfem::GlobGeometryRefiner.Refine(PointFieldBaseGeometry(location, i), ref, 1);
      const int base = (*field.bases)[i];
      const int npts = RefG->RefPts.GetNPoints();
      MFEM_VERIFY(base >= 0 && base + npts <= component_stride,
                  PointFieldLocationName(location)
                      << " point field buffer is missing data for "
                      << PointFieldEntityName(location) << " " << i << "!");
      for (int j = 0; j < npts; j++)
      {
        for (int c = 0; c < field.num_comp; c++)
        {
          payload_data[out++] =
              static_cast<scalar_t>(data[base + j + c * component_stride]);
        }
      }
    }
    MFEM_VERIFY(out == payload_count, "Packed point field payload has an invalid size!");
  };

  if (binary32)
  {
    std::vector<float> payload(payload_count);
    PackPayload(payload.data());
    WriteAll(os, payload.data(), payload.size() * sizeof(float));
  }
  else
  {
    std::vector<double> payload(payload_count);
    PackPayload(payload.data());
    WriteAll(os, payload.data(), payload.size() * sizeof(double));
  }
}

void CeedParaViewDataCollection::SavePointFieldVTU(MeshEntityType location,
                                                   std::ostream &os, int ref,
                                                   const std::string &name,
                                                   const PointField &field, Vector *scratch)
{
  MFEM_VERIFY(LocationMatchesOutput(location),
              PointFieldLocationName(location)
                  << " point fields require matching ParaView output!");
  MFEM_VERIFY((field.values || field.evaluator) && field.bases && field.num_comp > 0,
              "Invalid " << PointFieldLocationName(location)
                         << " point field registration!");
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == NumPointFieldEntities(location),
              PointFieldLocationName(location)
                  << " point field base offsets do not match the mesh!");

  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\" "
     << mfem::VTKComponentLabels(field.num_comp) << " "
     << "format=\"" << GetDataFormatString() << "\" >" << '\n';

  Vector values;
  const Vector *value_ptr = field.values;
  if (field.evaluator)
  {
    MFEM_VERIFY(scratch && scratch->Size() >= field.buffer_size,
                PointFieldLocationName(location)
                    << " point evaluator scratch buffer is too small!");
    values.MakeRef(*scratch, 0, field.buffer_size);
    field.evaluator(values);
    value_ptr = &values;
  }
  MFEM_VERIFY(value_ptr, PointFieldLocationName(location)
                             << " point field has no values to write!");

  MFEM_VERIFY(value_ptr->Size() % field.num_comp == 0,
              PointFieldLocationName(location)
                  << " point field buffer size is not divisible by its component count!");
  const int component_stride = value_ptr->Size() / field.num_comp;
  std::vector<char> buf;
  const double *data = value_ptr->HostRead();
  for (int i = 0; i < NumPointFieldEntities(location); i++)
  {
    const auto *RefG =
        mfem::GlobGeometryRefiner.Refine(PointFieldBaseGeometry(location, i), ref, 1);
    const int base = (*field.bases)[i];
    const int npts = RefG->RefPts.GetNPoints();
    MFEM_VERIFY(base >= 0 && base + npts <= component_stride,
                PointFieldLocationName(location)
                    << " point field buffer is missing data for "
                    << PointFieldEntityName(location) << " " << i << "!");
    for (int j = 0; j < npts; j++)
    {
      for (int c = 0; c < field.num_comp; c++)
      {
        mfem::WriteBinaryOrASCII(os, buf, data[base + j + c * component_stride], " ",
                                 pv_data_format);
      }
      if (pv_data_format == mfem::VTKFormat::ASCII)
      {
        os << '\n';
      }
    }
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;
}

void CeedParaViewDataCollection::SavePointFieldVTUAppendedHeader(MeshEntityType location,
                                                                 std::ostream &os, int ref,
                                                                 const std::string &name,
                                                                 const PointField &field,
                                                                 std::uint64_t offset)
{
  MFEM_VERIFY(LocationMatchesOutput(location),
              PointFieldLocationName(location)
                  << " point fields require matching ParaView output!");
  const std::uint64_t payload_size = PointFieldPayloadSize(location, ref, field);
  MFEM_VERIFY(payload_size <= std::numeric_limits<std::uint32_t>::max(),
              PointFieldLocationName(location)
                  << " point field is too large for the current appended VTU writer!");
  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\" "
     << mfem::VTKComponentLabels(field.num_comp) << " "
     << "format=\"appended\" offset=\"" << offset << "\" />" << '\n';
}

void CeedParaViewDataCollection::SavePointFieldVTUAppendedPayload(MeshEntityType location,
                                                                  std::ostream &os, int ref,
                                                                  const PointField &field,
                                                                  Vector *scratch)
{
  const std::uint64_t payload_size64 = PointFieldPayloadSize(location, ref, field);
  MFEM_VERIFY(payload_size64 <= std::numeric_limits<std::uint32_t>::max(),
              PointFieldLocationName(location)
                  << " point field is too large for the current appended VTU writer!");
  const std::uint32_t payload_size = static_cast<std::uint32_t>(payload_size64);
  os.write(reinterpret_cast<const char *>(&payload_size), sizeof(payload_size));

  Vector values;
  const Vector *value_ptr = field.values;
  if (field.evaluator)
  {
    MFEM_VERIFY(scratch && scratch->Size() >= field.buffer_size,
                PointFieldLocationName(location)
                    << " point evaluator scratch buffer is too small!");
    values.MakeRef(*scratch, 0, field.buffer_size);
    field.evaluator(values);
    value_ptr = &values;
  }
  MFEM_VERIFY(value_ptr, PointFieldLocationName(location)
                             << " point field has no values to write!");
  WritePointFieldValues(location, os, ref, field, *value_ptr);
}

void CeedParaViewDataCollection::SaveDataVTU(std::ostream &os, int ref)
{
  os << "<VTKFile type=\"UnstructuredGrid\"";
  if (GetCompressionLevel() != 0)
  {
    os << " compressor=\"vtkZLibDataCompressor\"";
  }
  os << " version=\"2.2\" byte_order=\"" << mfem::VTKByteOrder() << "\">\n";
  os << "<UnstructuredGrid>\n";
  mesh->PrintVTU(os, ref, pv_data_format, high_order_output, GetCompressionLevel(),
                 bdr_output);

  os << "<PointData >\n";
  for (FieldMapIterator it = field_map.begin(); it != field_map.end(); ++it)
  {
    MFEM_VERIFY(!bdr_output,
                "GridFunction output is not supported for ParaViewDataCollection on "
                "domain boundary!");
    SaveGFieldVTU(os, ref, it);
  }
  for (const auto &kv : GetCoeffFieldMap())
  {
    SaveCoeffFieldVTU(os, ref, kv.first, *kv.second);
  }
  for (const auto &kv : GetVCoeffFieldMap())
  {
    SaveVCoeffFieldVTU(os, ref, kv.first, *kv.second);
  }
  const bool appended_boundary_fields = UseAppendedPointFields(MeshEntityType::Boundary);
  const bool appended_domain_fields = UseAppendedPointFields(MeshEntityType::Domain);
  Vector boundary_point_scratch, domain_point_scratch;
  if (!boundary_point_fields.empty())
  {
    boundary_point_scratch.SetSize(MaxPointFieldBufferSize(boundary_point_fields));
    boundary_point_scratch.UseDevice(true);
  }
  if (!domain_point_fields.empty())
  {
    domain_point_scratch.SetSize(MaxPointFieldBufferSize(domain_point_fields));
    domain_point_scratch.UseDevice(true);
  }
  std::uint64_t appended_offset = 0;
  for (const auto &kv : boundary_point_fields)
  {
    if (appended_boundary_fields)
    {
      SavePointFieldVTUAppendedHeader(MeshEntityType::Boundary, os, ref, kv.first,
                                      kv.second, appended_offset);
      appended_offset += sizeof(std::uint32_t) +
                         PointFieldPayloadSize(MeshEntityType::Boundary, ref, kv.second);
    }
    else
    {
      SavePointFieldVTU(MeshEntityType::Boundary, os, ref, kv.first, kv.second,
                        &boundary_point_scratch);
    }
  }
  for (const auto &kv : domain_point_fields)
  {
    if (appended_domain_fields)
    {
      SavePointFieldVTUAppendedHeader(MeshEntityType::Domain, os, ref, kv.first, kv.second,
                                      appended_offset);
      appended_offset += sizeof(std::uint32_t) +
                         PointFieldPayloadSize(MeshEntityType::Domain, ref, kv.second);
    }
    else
    {
      SavePointFieldVTU(MeshEntityType::Domain, os, ref, kv.first, kv.second,
                        &domain_point_scratch);
    }
  }
  os << "</PointData>\n";
  os << "</Piece>\n";
  os << "</UnstructuredGrid>\n";
  if (appended_boundary_fields || appended_domain_fields)
  {
    os << "<AppendedData encoding=\"raw\">\n_";
    for (const auto &kv : boundary_point_fields)
    {
      SavePointFieldVTUAppendedPayload(MeshEntityType::Boundary, os, ref, kv.second,
                                       &boundary_point_scratch);
    }
    for (const auto &kv : domain_point_fields)
    {
      SavePointFieldVTUAppendedPayload(MeshEntityType::Domain, os, ref, kv.second,
                                       &domain_point_scratch);
    }
    os << "\n</AppendedData>\n";
  }
  os << "</VTKFile>" << std::endl;
}

void CeedParaViewDataCollection::Save()
{
  const std::string col_path = GenerateCollectionPath();
  {
    const std::string path = col_path + "/" + GenerateVTUPath();
    const int error_code = create_directory(path, mesh, myid);
    if (error_code)
    {
      error = WRITE_ERROR;
      MFEM_WARNING("Error creating directory: " << path);
      return;
    }
  }

  if (myid == 0 && !pvd_stream.is_open())
  {
    const std::string pvdname = col_path + "/" + GeneratePVDFileName();

    bool write_header = true;
    std::ifstream pvd_in;
    if (restart_mode && (pvd_in.open(pvdname, std::ios::binary), pvd_in.good()))
    {
      std::fstream::pos_type pos_begin = pvd_in.tellg();
      std::fstream::pos_type pos_end = pos_begin;

      std::regex regexp("timestep=\"([^[:space:]]+)\".*file=\"Cycle(\\d+)");
      std::smatch match;

      std::string line;
      while (getline(pvd_in, line))
      {
        if (regex_search(line, match, regexp))
        {
          MFEM_ASSERT(match.size() == 3, "Unable to parse DataSet");
          const double tvalue = std::stod(match[1]);
          if (tvalue >= GetTime())
          {
            break;
          }
          const int cvalue = std::stoi(match[2]);
          MFEM_VERIFY(cvalue < GetCycle(), "Cycle " << GetCycle()
                                                    << " is too small for restart mode: "
                                                       "trying to overwrite existing "
                                                       "data.");
          pos_end = pvd_in.tellg();
        }
      }

      const size_t count = pos_end - pos_begin;
      if (count != 0)
      {
        write_header = false;
        std::vector<char> buf(count);
        pvd_in.clear();
        pvd_in.seekg(pos_begin);
        pvd_in.read(buf.data(), count);
        pvd_in.close();
        pvd_stream.open(pvdname, std::ios::out | std::ios::trunc | std::ios::binary);
        pvd_stream.write(buf.data(), count);
        pvd_stream.close();
        pvd_stream.open(pvdname, std::ios::in | std::ios::out | std::ios::ate);
      }
    }
    if (write_header)
    {
      pvd_stream.open(pvdname, std::ios::out | std::ios::trunc);
      pvd_stream << "<?xml version=\"1.0\"?>\n";
      pvd_stream << "<VTKFile type=\"Collection\" version=\"2.2\"";
      pvd_stream << " byte_order=\"" << mfem::VTKByteOrder() << "\">\n";
      pvd_stream << "<Collection>" << std::endl;
    }
  }

  const std::string vtu_prefix = col_path + "/" + GenerateVTUPath() + "/";
  {
    const std::string os_str = vtu_prefix + GenerateVTUFileName("proc", myid);
    std::ofstream os(os_str, std::ios::binary);
    MFEM_VERIFY(os.is_open(), "Failed to open ofstream " << os_str);
    os.precision(precision);
    SaveDataVTU(os, levels_of_detail);
  }

  for (const auto &qfield : q_field_map)
  {
    MFEM_VERIFY(!bdr_output,
                "QuadratureFunction output is not supported for ParaViewDataCollection "
                "on domain boundary!");
    const std::string &field_name = qfield.first;
    const std::string os_str = vtu_prefix + GenerateVTUFileName(field_name, myid);
    std::ofstream os(os_str);
    MFEM_VERIFY(os.is_open(), "Failed to open ofstream " << os_str);
    qfield.second->SaveVTU(os, pv_data_format, GetCompressionLevel(), field_name);
  }

  if (myid == 0)
  {
    {
      const std::string os_str = vtu_prefix + GeneratePVTUFileName("data");
      std::ofstream pvtu_out(os_str);
      MFEM_VERIFY(pvtu_out.is_open(), "Failed to open ofstream " << os_str);
      WritePVTUHeader(pvtu_out);

      pvtu_out << "<PPointData>\n";
      for (auto &field_it : field_map)
      {
        const int vec_dim = field_it.second->VectorDim();
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << vec_dim << "\" "
                 << mfem::VTKComponentLabels(vec_dim) << " "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
      for (const auto &field_it : GetCoeffFieldMap())
      {
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"1\" "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
      for (const auto &field_it : GetVCoeffFieldMap())
      {
        const int vec_dim = field_it.second->GetVDim();
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << vec_dim << "\" "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
      for (const auto &field_it : boundary_point_fields)
      {
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << field_it.second.num_comp
                 << "\" " << mfem::VTKComponentLabels(field_it.second.num_comp) << " "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
      for (const auto &field_it : domain_point_fields)
      {
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << field_it.second.num_comp
                 << "\" " << mfem::VTKComponentLabels(field_it.second.num_comp) << " "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
      pvtu_out << "</PPointData>\n";

      pvtu_out << "<PCellData>\n";
      pvtu_out << "\t<PDataArray type=\"Int32\" Name=\"attribute\" "
                  "NumberOfComponents=\"1\""
               << " format=\"" << GetDataFormatString() << "\"/>\n";
      pvtu_out << "</PCellData>\n";

      WritePVTUFooter(pvtu_out, "proc");
    }

    pvd_stream << "<DataSet timestep=\"" << GetTime() << "\" group=\"\" part=\"" << 0
               << "\" file=\"" << GeneratePVTUPath() + "/" + GeneratePVTUFileName("data")
               << "\" name=\"mesh\"/>\n";

    for (auto &q_field : q_field_map)
    {
      const std::string &q_field_name = q_field.first;
      const std::string q_fname =
          GeneratePVTUPath() + "/" + GeneratePVTUFileName(q_field_name);
      const std::string os_str = col_path + "/" + q_fname;
      std::ofstream pvtu_out(os_str);
      MFEM_VERIFY(pvtu_out.is_open(), "Failed to open ofstream " << os_str);
      WritePVTUHeader(pvtu_out);
      const int vec_dim = q_field.second->GetVDim();
      pvtu_out << "<PPointData>\n";
      pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
               << q_field_name << "\" NumberOfComponents=\"" << vec_dim << "\" "
               << mfem::VTKComponentLabels(vec_dim) << " "
               << "format=\"" << GetDataFormatString() << "\" />\n";
      pvtu_out << "</PPointData>\n";
      WritePVTUFooter(pvtu_out, q_field_name);

      pvd_stream << "<DataSet timestep=\"" << GetTime() << "\" group=\"\" part=\"" << 0
                 << "\" file=\"" << q_fname << "\" name=\"" << q_field_name << "\"/>\n";
    }
    pvd_stream.flush();
    const std::fstream::pos_type pos = pvd_stream.tellp();
    pvd_stream << "</Collection>\n";
    pvd_stream << "</VTKFile>" << std::endl;
    pvd_stream.seekp(pos);
  }
}

}  // namespace palace
