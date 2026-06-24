// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/ceedparaviewdatacollection.hpp"

#include <cstdint>
#include <cstring>
#include <limits>
#include <regex>
#include <vector>
#include <mfem/mesh/vtk.hpp>

namespace palace
{

namespace
{

template <typename T>
void AppendPod(std::vector<char> &buf, const T &value)
{
  const auto offset = buf.size();
  buf.resize(offset + sizeof(T));
  std::memcpy(buf.data() + offset, &value, sizeof(T));
}

std::size_t VTKRealSize(mfem::VTKFormat format)
{
  return (format == mfem::VTKFormat::BINARY32) ? sizeof(float) : sizeof(double);
}

void ReserveBinary(std::vector<char> &buf, std::size_t n, std::size_t bytes_per_entry)
{
  if (buf.capacity() < n * bytes_per_entry)
  {
    buf.reserve(n * bytes_per_entry);
  }
}

void WriteVTKReal(std::ostream &os, std::vector<char> &buf, double value,
                  const char *suffix, mfem::VTKFormat format)
{
  if (format == mfem::VTKFormat::BINARY32)
  {
    AppendPod(buf, static_cast<float>(value));
  }
  else if (format == mfem::VTKFormat::BINARY)
  {
    AppendPod(buf, value);
  }
  else
  {
    mfem::WriteBinaryOrASCII(os, buf, value, suffix, format);
  }
}

void WriteVTKInt32(std::ostream &os, std::vector<char> &buf, int value,
                   const char *suffix, mfem::VTKFormat format)
{
  if (format == mfem::VTKFormat::ASCII)
  {
    mfem::WriteBinaryOrASCII(os, buf, value, suffix, format);
  }
  else
  {
    AppendPod(buf, static_cast<std::int32_t>(value));
  }
}

void WriteVTKUInt8(std::ostream &os, std::vector<char> &buf, std::uint8_t value,
                   const char *suffix, mfem::VTKFormat format)
{
  if (format == mfem::VTKFormat::ASCII)
  {
    mfem::WriteBinaryOrASCII(os, buf, value, suffix, format);
  }
  else
  {
    AppendPod(buf, value);
  }
}

}  // namespace


void CeedParaViewDataCollection::RegisterDomainPointField(const std::string &field_name,
                                                          const Vector &values,
                                                          const std::vector<int> &bases,
                                                          int num_comp)
{
  MFEM_VERIFY(num_comp > 0, "Domain point field must have at least one component!");
  domain_point_fields[field_name] = PointField{&values, &bases, num_comp};
}

void CeedParaViewDataCollection::DeregisterDomainPointField(const std::string &field_name)
{
  domain_point_fields.erase(field_name);
}

void CeedParaViewDataCollection::RegisterBoundaryPointField(const std::string &field_name,
                                                            const Vector &values,
                                                            const std::vector<int> &bases,
                                                            int num_comp)
{
  MFEM_VERIFY(num_comp > 0, "Boundary point field must have at least one component!");
  boundary_point_fields[field_name] = PointField{&values, &bases, num_comp};
}

void CeedParaViewDataCollection::DeregisterBoundaryPointField(const std::string &field_name)
{
  boundary_point_fields.erase(field_name);
}

void CeedParaViewDataCollection::SavePointFieldVTU(std::ostream &os, int ref,
                                                   const std::string &name,
                                                   const PointField &field, bool boundary)
{
  MFEM_VERIFY(boundary == bdr_output,
              "Precomputed point field registration does not match the ParaView "
              "collection domain/boundary mode!");
  MFEM_VERIFY(field.values && field.bases && field.num_comp > 0,
              "Invalid precomputed point field registration!");

  const int ne = boundary ? mesh->GetNBE() : mesh->GetNE();
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == ne,
              "Point field base offsets do not match the mesh!");

  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\""
     << " format=\"" << GetDataFormatString() << "\" >" << '\n';

  std::vector<char> buf;
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, static_cast<std::size_t>(field.values->Size()),
                  VTKRealSize(pv_data_format));
  }
  const double *data = field.values->HostRead();
  for (int i = 0; i < ne; i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(
        boundary ? mesh->GetBdrElementBaseGeometry(i) : mesh->GetElementBaseGeometry(i),
        ref, 1);
    const int base = (*field.bases)[i];
    const int npts = RefG->RefPts.GetNPoints();
    MFEM_VERIFY(base >= 0 && base + field.num_comp * npts <= field.values->Size(),
                "Point field buffer is missing data for element " << i << "!");
    for (int j = 0; j < npts; j++)
    {
      for (int c = 0; c < field.num_comp; c++)
      {
        WriteVTKReal(os, buf, data[base + field.num_comp * j + c], " ",
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

void CeedParaViewDataCollection::SaveMeshVTU(std::ostream &os, int ref)
{
  // MFEM's high-order VTU mesh writer recomputes the same local VTK connectivity
  // permutation for every element. Palace ParaView output commonly uses one geometry
  // and high-order cells, so cache that per-geometry pattern and pre-reserve the binary
  // buffers before forwarding them to MFEM's base64/zlib encoder.
  if (!high_order_output)
  {
    mesh->PrintVTU(os, ref, pv_data_format, high_order_output, GetCompressionLevel(),
                   bdr_output);
    return;
  }

  const char *fmt_str = (pv_data_format == mfem::VTKFormat::ASCII) ? "ascii" : "binary";
  const char *type_str = (pv_data_format != mfem::VTKFormat::BINARY32) ? "Float64"
                                                                        : "Float32";
  const int ne = bdr_output ? mesh->GetNBE() : mesh->GetNE();

  auto GetGeom = [&](int i)
  {
    return bdr_output ? mesh->GetBdrElementGeometry(i) : mesh->GetElementBaseGeometry(i);
  };

  std::size_t num_points = 0;
  for (int i = 0; i < ne; i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(GetGeom(i), ref, 1);
    num_points += static_cast<std::size_t>(RefG->RefPts.GetNPoints());
  }
  MFEM_VERIFY(num_points <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "VTU output exceeds 32-bit point indexing!");

  os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << ne
     << "\">\n";

  std::vector<char> buf;
  mfem::DenseMatrix pmat;

  os << "<Points>\n";
  os << "<DataArray type=\"" << type_str
     << "\" NumberOfComponents=\"3\" format=\"" << fmt_str << "\">\n";
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, 3 * num_points, VTKRealSize(pv_data_format));
  }
  for (int i = 0; i < ne; i++)
  {
    const auto *RefG = mfem::GlobGeometryRefiner.Refine(GetGeom(i), ref, 1);
    if (bdr_output)
    {
      mesh->GetBdrElementTransformation(i)->Transform(RefG->RefPts, pmat);
    }
    else
    {
      mesh->GetElementTransformation(i)->Transform(RefG->RefPts, pmat);
    }
    for (int j = 0; j < pmat.Width(); j++)
    {
      WriteVTKReal(os, buf, pmat(0, j), " ", pv_data_format);
      WriteVTKReal(os, buf, (pmat.Height() > 1) ? pmat(1, j) : 0.0, " ",
                   pv_data_format);
      WriteVTKReal(os, buf, (pmat.Height() > 2) ? pmat(2, j) : 0.0, "",
                   pv_data_format);
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
  os << "</Points>" << std::endl;

  os << "<Cells>" << std::endl;
  os << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << fmt_str
     << "\">" << std::endl;
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, num_points, sizeof(std::int32_t));
  }
  std::map<mfem::Geometry::Type, mfem::Array<int>> connectivity_cache;
  std::vector<int> offsets;
  offsets.reserve(static_cast<std::size_t>(ne));
  int point_offset = 0;
  for (int i = 0; i < ne; i++)
  {
    const auto geom = GetGeom(i);
    auto &local_connectivity = connectivity_cache[geom];
    if (local_connectivity.Size() == 0)
    {
      mfem::CreateVTKElementConnectivity(local_connectivity, geom, ref);
    }
    for (int j = 0; j < local_connectivity.Size(); j++)
    {
      WriteVTKInt32(os, buf, point_offset + local_connectivity[j], " ",
                    pv_data_format);
    }
    if (pv_data_format == mfem::VTKFormat::ASCII)
    {
      os << '\n';
    }
    point_offset += local_connectivity.Size();
    offsets.push_back(point_offset);
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;

  os << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"" << fmt_str
     << "\">" << std::endl;
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, offsets.size(), sizeof(std::int32_t));
  }
  for (const int offset : offsets)
  {
    WriteVTKInt32(os, buf, offset, "\n", pv_data_format);
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;

  os << "<DataArray type=\"UInt8\" Name=\"types\" format=\"" << fmt_str << "\">"
     << std::endl;
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, static_cast<std::size_t>(ne), sizeof(std::uint8_t));
  }
  const int *vtk_geom_map = mfem::VTKGeometry::HighOrderMap;
  for (int i = 0; i < ne; i++)
  {
    WriteVTKUInt8(os, buf, static_cast<std::uint8_t>(vtk_geom_map[GetGeom(i)]), "\n",
                  pv_data_format);
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;
  os << "</Cells>" << std::endl;

  os << "<CellData Scalars=\"attribute\">" << std::endl;
  os << "<DataArray type=\"Int32\" Name=\"attribute\" format=\"" << fmt_str
     << "\">" << std::endl;
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    ReserveBinary(buf, static_cast<std::size_t>(ne), sizeof(std::int32_t));
  }
  for (int i = 0; i < ne; i++)
  {
    const int attr = bdr_output ? mesh->GetBdrAttribute(i) : mesh->GetAttribute(i);
    WriteVTKInt32(os, buf, attr, "\n", pv_data_format);
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;
  os << "</CellData>" << std::endl;
}

void CeedParaViewDataCollection::SaveGFieldVTU(std::ostream &os, int ref,
                                               const FieldMapIterator &it)
{
  mfem::Vector val;
  mfem::DenseMatrix vval, pmat;
  std::vector<char> buf;
  const int vec_dim = it->second->VectorDim();
  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << it->first
     << "\" NumberOfComponents=\"" << vec_dim << "\" "
     << mfem::VTKComponentLabels(vec_dim) << " "
     << "format=\"" << GetDataFormatString() << "\" >" << '\n';

  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    std::size_t num_values = 0;
    for (int i = 0; i < mesh->GetNE(); i++)
    {
      const auto *RefG = mfem::GlobGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                                          ref, 1);
      num_values += static_cast<std::size_t>(RefG->RefPts.GetNPoints()) * vec_dim;
    }
    ReserveBinary(buf, num_values, VTKRealSize(pv_data_format));
  }

  if (vec_dim == 1)
  {
    for (int i = 0; i < mesh->GetNE(); i++)
    {
      const auto *RefG = mfem::GlobGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                                          ref, 1);
      it->second->GetValues(i, RefG->RefPts, val, pmat);
      for (int j = 0; j < val.Size(); j++)
      {
        WriteVTKReal(os, buf, val(j), "\n", pv_data_format);
      }
    }
  }
  else
  {
    for (int i = 0; i < mesh->GetNE(); i++)
    {
      const auto *RefG = mfem::GlobGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                                          ref, 1);
      it->second->GetVectorValues(i, RefG->RefPts, vval, pmat);
      for (int jj = 0; jj < vval.Width(); jj++)
      {
        for (int ii = 0; ii < vval.Height(); ii++)
        {
          WriteVTKReal(os, buf, vval(ii, jj), " ", pv_data_format);
        }
        if (pv_data_format == mfem::VTKFormat::ASCII)
        {
          os << '\n';
        }
      }
    }
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;
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
  SaveMeshVTU(os, ref);

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
  for (const auto &kv : domain_point_fields)
  {
    SavePointFieldVTU(os, ref, kv.first, kv.second, false);
  }
  for (const auto &kv : boundary_point_fields)
  {
    SavePointFieldVTU(os, ref, kv.first, kv.second, true);
  }
  os << "</PointData>\n";
  os << "</Piece>\n";
  os << "</UnstructuredGrid>\n";
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
      auto WritePointFieldHeader = [&](const auto &field_it)
      {
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << field_it.second.num_comp
                 << "\" "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      };
      for (const auto &field_it : domain_point_fields)
      {
        WritePointFieldHeader(field_it);
      }
      for (const auto &field_it : boundary_point_fields)
      {
        WritePointFieldHeader(field_it);
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
