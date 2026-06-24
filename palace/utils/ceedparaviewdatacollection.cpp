// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/ceedparaviewdatacollection.hpp"

#include <limits>
#include <regex>
#include <utility>
#include <vector>
#include <mfem/mesh/vtk.hpp>

namespace palace
{

void CeedParaViewDataCollection::RegisterBoundaryPointField(
    const std::string &field_name, const Vector &values, const std::vector<int> &bases,
    int num_comp)
{
  MFEM_VERIFY(num_comp > 0, "Boundary point field must have at least one component!");
  boundary_point_fields[field_name] = BoundaryPointField{&values, {}, &bases, num_comp,
                                                        values.Size()};
}

void CeedParaViewDataCollection::RegisterBoundaryPointEvaluator(
    const std::string &field_name, std::function<void(Vector &)> evaluator,
    const std::vector<int> &bases, int num_comp, int buffer_size)
{
  MFEM_VERIFY(evaluator, "Boundary point evaluator must be callable!");
  MFEM_VERIFY(num_comp > 0, "Boundary point field must have at least one component!");
  MFEM_VERIFY(buffer_size >= 0, "Boundary point evaluator buffer size is invalid!");
  boundary_point_fields[field_name] =
      BoundaryPointField{nullptr, std::move(evaluator), &bases, num_comp, buffer_size};
}

void CeedParaViewDataCollection::DeregisterBoundaryPointField(
    const std::string &field_name)
{
  boundary_point_fields.erase(field_name);
}

bool CeedParaViewDataCollection::UseAppendedBoundaryPointFields() const
{
  return bdr_output && !boundary_point_fields.empty() &&
         pv_data_format != mfem::VTKFormat::ASCII && GetCompressionLevel() == 0;
}

std::uint64_t CeedParaViewDataCollection::BoundaryPointFieldPayloadSize(
    int ref, const BoundaryPointField &field) const
{
  MFEM_VERIFY((field.values || field.evaluator) && field.bases && field.num_comp > 0,
              "Invalid boundary point field registration!");
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == mesh->GetNBE(),
              "Boundary point field base offsets do not match the mesh boundary!");

  std::uint64_t value_count = 0;
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(
        mesh->GetBdrElementBaseGeometry(i), ref, 1);
    value_count += static_cast<std::uint64_t>(field.num_comp) *
                   static_cast<std::uint64_t>(RefG->RefPts.GetNPoints());
  }
  const std::uint64_t scalar_size =
      pv_data_format == mfem::VTKFormat::BINARY32 ? sizeof(float) : sizeof(double);
  return value_count * scalar_size;
}

void CeedParaViewDataCollection::WriteBoundaryPointFieldValues(
    std::ostream &os, int ref, const BoundaryPointField &field,
    const Vector &values) const
{
  const double *data = values.HostRead();
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(
        mesh->GetBdrElementBaseGeometry(i), ref, 1);
    const int base = (*field.bases)[i];
    const int npts = RefG->RefPts.GetNPoints();
    MFEM_VERIFY(base >= 0 && base + field.num_comp * npts <= values.Size(),
                "Boundary point field buffer is missing data for boundary element " << i
                                                                                     << "!");
    for (int j = 0; j < npts; j++)
    {
      for (int c = 0; c < field.num_comp; c++)
      {
        const double value = data[base + field.num_comp * j + c];
        if (pv_data_format == mfem::VTKFormat::BINARY32)
        {
          const float value32 = static_cast<float>(value);
          os.write(reinterpret_cast<const char *>(&value32), sizeof(value32));
        }
        else
        {
          os.write(reinterpret_cast<const char *>(&value), sizeof(value));
        }
      }
    }
  }
}

void CeedParaViewDataCollection::SaveBoundaryPointFieldVTU(
    std::ostream &os, int ref, const std::string &name, const BoundaryPointField &field)
{
  MFEM_VERIFY(bdr_output, "Boundary point fields require boundary ParaView output!");
  MFEM_VERIFY((field.values || field.evaluator) && field.bases && field.num_comp > 0,
              "Invalid boundary point field registration!");
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == mesh->GetNBE(),
              "Boundary point field base offsets do not match the mesh boundary!");

  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\""
     << " format=\"" << GetDataFormatString() << "\" >" << '\n';

  Vector values;
  const Vector *value_ptr = field.values;
  if (field.evaluator)
  {
    values.SetSize(field.buffer_size);
    values.UseDevice(true);
    field.evaluator(values);
    value_ptr = &values;
  }
  MFEM_VERIFY(value_ptr, "Boundary point field has no values to write!");

  std::vector<char> buf;
  const double *data = value_ptr->HostRead();
  for (int i = 0; i < mesh->GetNBE(); i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(
        mesh->GetBdrElementBaseGeometry(i), ref, 1);
    const int base = (*field.bases)[i];
    const int npts = RefG->RefPts.GetNPoints();
    MFEM_VERIFY(base >= 0 && base + field.num_comp * npts <= value_ptr->Size(),
                "Boundary point field buffer is missing data for boundary element " << i
                                                                                     << "!");
    for (int j = 0; j < npts; j++)
    {
      for (int c = 0; c < field.num_comp; c++)
      {
        mfem::WriteBinaryOrASCII(os, buf, data[base + field.num_comp * j + c], " ",
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

void CeedParaViewDataCollection::SaveBoundaryPointFieldVTUAppendedHeader(
    std::ostream &os, int ref, const std::string &name, const BoundaryPointField &field,
    std::uint64_t offset)
{
  MFEM_VERIFY(bdr_output, "Boundary point fields require boundary ParaView output!");
  const std::uint64_t payload_size = BoundaryPointFieldPayloadSize(ref, field);
  MFEM_VERIFY(payload_size <= std::numeric_limits<std::uint32_t>::max(),
              "Boundary point field is too large for the current appended VTU writer!");
  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\""
     << " format=\"appended\" offset=\"" << offset << "\" />" << '\n';
}

void CeedParaViewDataCollection::SaveBoundaryPointFieldVTUAppendedPayload(
    std::ostream &os, int ref, const BoundaryPointField &field)
{
  const std::uint64_t payload_size64 = BoundaryPointFieldPayloadSize(ref, field);
  MFEM_VERIFY(payload_size64 <= std::numeric_limits<std::uint32_t>::max(),
              "Boundary point field is too large for the current appended VTU writer!");
  const std::uint32_t payload_size = static_cast<std::uint32_t>(payload_size64);
  os.write(reinterpret_cast<const char *>(&payload_size), sizeof(payload_size));

  Vector values;
  const Vector *value_ptr = field.values;
  if (field.evaluator)
  {
    values.SetSize(field.buffer_size);
    values.UseDevice(true);
    field.evaluator(values);
    value_ptr = &values;
  }
  MFEM_VERIFY(value_ptr, "Boundary point field has no values to write!");
  WriteBoundaryPointFieldValues(os, ref, field, *value_ptr);
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
  const bool appended_boundary_fields = UseAppendedBoundaryPointFields();
  std::uint64_t appended_offset = 0;
  for (const auto &kv : boundary_point_fields)
  {
    if (appended_boundary_fields)
    {
      SaveBoundaryPointFieldVTUAppendedHeader(os, ref, kv.first, kv.second,
                                              appended_offset);
      appended_offset += sizeof(std::uint32_t) +
                         BoundaryPointFieldPayloadSize(ref, kv.second);
    }
    else
    {
      SaveBoundaryPointFieldVTU(os, ref, kv.first, kv.second);
    }
  }
  os << "</PointData>\n";
  os << "</Piece>\n";
  os << "</UnstructuredGrid>\n";
  if (appended_boundary_fields)
  {
    os << "<AppendedData encoding=\"raw\">\n_";
    for (const auto &kv : boundary_point_fields)
    {
      SaveBoundaryPointFieldVTUAppendedPayload(os, ref, kv.second);
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
                 << field_it.first << "\" NumberOfComponents=\""
                 << field_it.second.num_comp << "\" "
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

    pvd_stream << "<DataSet timestep=\"" << GetTime()
               << "\" group=\"\" part=\"" << 0 << "\" file=\""
               << GeneratePVTUPath() + "/" + GeneratePVTUFileName("data")
               << "\" name=\"mesh\"/>\n";

    for (auto &q_field : q_field_map)
    {
      const std::string &q_field_name = q_field.first;
      const std::string q_fname = GeneratePVTUPath() + "/" + GeneratePVTUFileName(q_field_name);
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

      pvd_stream << "<DataSet timestep=\"" << GetTime()
                 << "\" group=\"\" part=\"" << 0 << "\" file=\"" << q_fname
                 << "\" name=\"" << q_field_name << "\"/>\n";
    }
    pvd_stream.flush();
    const std::fstream::pos_type pos = pvd_stream.tellp();
    pvd_stream << "</Collection>\n";
    pvd_stream << "</VTKFile>" << std::endl;
    pvd_stream.seekp(pos);
  }
}

}  // namespace palace
