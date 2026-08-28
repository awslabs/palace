// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/ceedparaviewdatacollection.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>
#include <mfem/mesh/vtk.hpp>
#include "utils/communication.hpp"

namespace palace
{

std::atomic<long long> CeedParaViewDataCollection::individual_evaluator_calls{0};
std::atomic<long long> CeedParaViewDataCollection::provider_calls{0};
std::atomic<long long> CeedParaViewDataCollection::batch_view_fields{0};
std::atomic<long long> CeedParaViewDataCollection::avoided_slice_copies{0};
std::atomic<long long> CeedParaViewDataCollection::avoided_slice_copy_bytes{0};

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
  cached_point_headers_vtu.clear();
  PointFields(location)[field_name] =
      PointField{&values, {}, &bases, num_comp, values.Size(), false};
}

void CeedParaViewDataCollection::RegisterPointEvaluator(
    MeshEntityType location, const std::string &field_name,
    std::function<void(Vector &)> evaluator, const std::vector<int> &bases, int num_comp,
    int buffer_size, bool point_major)
{
  MFEM_VERIFY(evaluator, PointFieldLocationName(location)
                             << " point evaluator must be callable!");
  MFEM_VERIFY(num_comp > 0, PointFieldLocationName(location)
                                << " point field must have at least one component!");
  MFEM_VERIFY(buffer_size >= 0, PointFieldLocationName(location)
                                    << " point evaluator buffer size is invalid!");
  cached_point_headers_vtu.clear();
  PointFields(location)[field_name] =
      PointField{nullptr, std::move(evaluator), &bases, num_comp, buffer_size, point_major};
}

void CeedParaViewDataCollection::DeregisterPointField(MeshEntityType location,
                                                      const std::string &field_name)
{
  cached_point_headers_vtu.clear();
  PointFields(location).erase(field_name);
}

void CeedParaViewDataCollection::ResetPointFieldProviders()
{
  std::set<PointFieldProvider *> providers;
  for (const auto *fields : {&boundary_point_fields, &domain_point_fields})
  {
    for (const auto &[name, field] : *fields)
    {
      if (field.provider && providers.insert(field.provider.get()).second)
      {
        field.provider->cached_values = nullptr;
      }
    }
  }
}

const Vector &CeedParaViewDataCollection::GetPointFieldValues(MeshEntityType location,
                                                              const PointField &field,
                                                              Vector &view,
                                                              Vector *scratch) const
{
  if (field.provider)
  {
    auto &provider = *field.provider;
    if (!provider.cached_values)
    {
      MFEM_VERIFY(provider.evaluator, PointFieldLocationName(location)
                                          << " point field provider must be callable!");
      const Vector &packed = provider.evaluator();
      MFEM_VERIFY(packed.Size() == provider.expected_size,
                  PointFieldLocationName(location)
                      << " point field provider returned an invalid packed buffer size!");
      provider.cached_values = &packed;
      ++provider_calls;
    }
    MFEM_VERIFY(field.provider_offset >= 0 && field.buffer_size >= 0 &&
                    field.provider_offset <= provider.expected_size - field.buffer_size,
                PointFieldLocationName(location)
                    << " point field provider slice is outside its packed buffer!");
    // MakeRef is a non-owning MFEM view. The provider owns the packed Vector for this
    // synchronous save, so no component slice is copied into the normal scratch buffer.
    view.MakeRef(const_cast<Vector &>(*provider.cached_values), field.provider_offset,
                 field.buffer_size);
    ++batch_view_fields;
    ++avoided_slice_copies;
    avoided_slice_copy_bytes += static_cast<long long>(field.buffer_size) * sizeof(double);
    return view;
  }
  if (field.evaluator)
  {
    MFEM_VERIFY(scratch && scratch->Size() >= field.buffer_size,
                PointFieldLocationName(location)
                    << " point evaluator scratch buffer is too small!");
    view.MakeRef(*scratch, 0, field.buffer_size);
    field.evaluator(view);
    ++individual_evaluator_calls;
    return view;
  }
  MFEM_VERIFY(field.values, PointFieldLocationName(location)
                                << " point field has no values to write!");
  return *field.values;
}

void CeedParaViewDataCollection::WriteCellFieldVTU(std::ostream &os,
                                                   const std::string &name,
                                                   const CellField &field) const
{
  MFEM_VERIFY(!bdr_output, "Domain cell fields require domain ParaView output!");
  MFEM_VERIFY(field.values && field.num_comp > 0,
              "Invalid domain cell field registration!");
  MFEM_VERIFY(field.values->Size() == mesh->GetNE() * field.num_comp,
              "Domain cell field size does not match the mesh cell count!");

  os << "<DataArray type=\"" << GetDataTypeString() << "\" Name=\"" << name
     << "\" NumberOfComponents=\"" << field.num_comp << "\" "
     << mfem::VTKComponentLabels(field.num_comp) << " "
     << "format=\"" << GetDataFormatString() << "\" >" << '\n';

  std::vector<char> buf;
  const double *data = field.values->HostRead();
  for (int i = 0; i < mesh->GetNE(); i++)
  {
    for (int c = 0; c < field.num_comp; c++)
    {
      mfem::WriteBinaryOrASCII(os, buf, data[i * field.num_comp + c], " ", pv_data_format);
    }
    if (pv_data_format == mfem::VTKFormat::ASCII)
    {
      os << '\n';
    }
  }
  if (pv_data_format != mfem::VTKFormat::ASCII)
  {
    mfem::WriteBase64WithSizeAndClear(os, buf, GetCompressionLevel());
  }
  os << "</DataArray>" << std::endl;
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

void CeedParaViewDataCollection::RegisterBoundaryPointFieldBatch(
    const std::vector<BoundaryPointFieldDescriptor> &fields, const std::vector<int> &bases,
    std::function<const Vector &()> provider, int expected_packed_size)
{
  MFEM_VERIFY(!fields.empty(), "Boundary point field batch must not be empty!");
  MFEM_VERIFY(provider, "Boundary point field batch provider must be callable!");
  MFEM_VERIFY(expected_packed_size > 0,
              "Boundary point field batch packed buffer size is invalid!");
  MFEM_VERIFY(static_cast<int>(bases.size()) ==
                  NumPointFieldEntities(MeshEntityType::Boundary),
              "Boundary point field batch base offsets do not match the mesh!");

  auto state = std::make_shared<PointFieldProvider>();
  state->evaluator = std::move(provider);
  state->expected_size = expected_packed_size;
  cached_point_headers_vtu.clear();
  for (const auto &field : fields)
  {
    MFEM_VERIFY(!field.name.empty(), "Boundary point field batch name must not be empty!");
    MFEM_VERIFY(field.num_comp > 0,
                "Boundary point field batch field must have at least one component!");
    MFEM_VERIFY(field.size > 0 && field.size % field.num_comp == 0,
                "Boundary point field batch slice has an invalid component size!");
    MFEM_VERIFY(field.offset >= 0 && field.offset <= expected_packed_size - field.size,
                "Boundary point field batch slice is outside its packed buffer!");
    boundary_point_fields[field.name] = PointField{
        nullptr, {}, &bases, field.num_comp, field.size, false, state, field.offset};
  }
}

void CeedParaViewDataCollection::RegisterDomainPointEvaluator(
    const std::string &field_name, std::function<void(Vector &)> evaluator,
    const std::vector<int> &bases, int num_comp, int buffer_size, bool point_major)
{
  RegisterPointEvaluator(MeshEntityType::Domain, field_name, std::move(evaluator), bases,
                         num_comp, buffer_size, point_major);
}

void CeedParaViewDataCollection::DeregisterBoundaryPointField(const std::string &field_name)
{
  DeregisterPointField(MeshEntityType::Boundary, field_name);
}

void CeedParaViewDataCollection::RegisterDomainCellField(const std::string &field_name,
                                                         const Vector &values, int num_comp)
{
  MFEM_VERIFY(num_comp > 0, "Domain cell field must have at least one component!");
  MFEM_VERIFY(values.Size() == mesh->GetNE() * num_comp,
              "Domain cell field size does not match the mesh cell count!");
  cached_point_headers_vtu.clear();
  domain_cell_fields[field_name] = CellField{&values, num_comp};
}

void CeedParaViewDataCollection::DeregisterDomainPointField(const std::string &field_name)
{
  DeregisterPointField(MeshEntityType::Domain, field_name);
}

void CeedParaViewDataCollection::DeregisterDomainCellField(const std::string &field_name)
{
  cached_point_headers_vtu.clear();
  domain_cell_fields.erase(field_name);
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
    // Batch members view their provider's packed storage directly and must not reserve a
    // redundant per-field scratch buffer.
    if (!field.provider)
    {
      size = std::max(size, field.buffer_size);
    }
  }
  return size;
}

std::uint64_t
CeedParaViewDataCollection::PointFieldPayloadSize(MeshEntityType location,
                                                  const PointField &field) const
{
  MFEM_VERIFY((field.values || field.evaluator || field.provider) && field.bases &&
                  field.num_comp > 0,
              "Invalid " << PointFieldLocationName(location)
                         << " point field registration!");
  MFEM_VERIFY(static_cast<int>(field.bases->size()) == NumPointFieldEntities(location),
              PointFieldLocationName(location)
                  << " point field base offsets do not match the mesh!");
  MFEM_VERIFY(field.buffer_size >= 0, PointFieldLocationName(location)
                                          << " point field buffer size is invalid!");

  const std::uint64_t scalar_size =
      pv_data_format == mfem::VTKFormat::BINARY32 ? sizeof(float) : sizeof(double);
  return static_cast<std::uint64_t>(field.buffer_size) * scalar_size;
}

void CeedParaViewDataCollection::WritePointFieldValues(MeshEntityType location,
                                                       std::ostream &os, int ref,
                                                       const std::string &name,
                                                       const PointField &field,
                                                       const Vector &values) const
{
  MFEM_VERIFY(values.Size() % field.num_comp == 0,
              PointFieldLocationName(location)
                  << " point field buffer size is not divisible by its component count!");
  const int component_stride = values.Size() / field.num_comp;
  const double *data = values.HostRead();

  const std::uint64_t payload_size = PointFieldPayloadSize(location, field);
  const bool binary32 = pv_data_format == mfem::VTKFormat::BINARY32;
  const std::size_t scalar_size = binary32 ? sizeof(float) : sizeof(double);
  MFEM_VERIFY(payload_size % scalar_size == 0,
              "Point field payload size is not divisible by scalar size!");
  const std::size_t payload_count = static_cast<std::size_t>(payload_size / scalar_size);

  if (field.point_major)
  {
    int expected_points = 0;
    for (int i = 0; i < NumPointFieldEntities(location); i++)
    {
      expected_points +=
          mfem::GlobGeometryRefiner.Refine(PointFieldBaseGeometry(location, i), ref, 1)
              ->RefPts.GetNPoints();
    }
    MFEM_VERIFY(component_stride == expected_points,
                PointFieldLocationName(location)
                    << " point-major field lattice does not match the VTU mesh!");
    MFEM_VERIFY(static_cast<std::size_t>(values.Size()) == payload_count,
                PointFieldLocationName(location)
                    << " point-major field buffer has an invalid size!");
    if (binary32)
    {
      auto payload = std::make_unique<float[]>(payload_count);
      for (std::size_t i = 0; i < payload_count; i++)
      {
        payload[i] = static_cast<float>(data[i]);
      }
      WriteAll(os, payload.get(), payload_count * sizeof(float));
    }
    else
    {
      WriteAll(os, data, payload_count * sizeof(double));
    }
    return;
  }

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
                      << " point field '" << name << "' buffer is missing data for "
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
  MFEM_VERIFY((field.values || field.evaluator || field.provider) && field.bases &&
                  field.num_comp > 0,
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
  const Vector *value_ptr = &GetPointFieldValues(location, field, values, scratch);

  MFEM_VERIFY(value_ptr->Size() % field.num_comp == 0,
              PointFieldLocationName(location)
                  << " point field buffer size is not divisible by its component count!");
  MFEM_VERIFY(value_ptr->Size() == field.buffer_size,
              PointFieldLocationName(location)
                  << " point field buffer has an invalid size!");
  const int component_stride = value_ptr->Size() / field.num_comp;
  std::vector<char> buf;
  const double *data = value_ptr->HostRead();
  auto WriteTuple = [&](const auto &Value)
  {
    for (int c = 0; c < field.num_comp; c++)
    {
      mfem::WriteBinaryOrASCII(os, buf, Value(c), " ", pv_data_format);
    }
    if (pv_data_format == mfem::VTKFormat::ASCII)
    {
      os << '\n';
    }
  };
  if (field.point_major)
  {
    int expected_points = 0;
    for (int i = 0; i < NumPointFieldEntities(location); i++)
    {
      expected_points +=
          mfem::GlobGeometryRefiner.Refine(PointFieldBaseGeometry(location, i), ref, 1)
              ->RefPts.GetNPoints();
    }
    MFEM_VERIFY(component_stride == expected_points,
                PointFieldLocationName(location)
                    << " point-major field lattice does not match the VTU mesh!");
    for (int p = 0; p < component_stride; p++)
    {
      WriteTuple([&](int c) { return data[p * field.num_comp + c]; });
    }
  }
  else
  {
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
        WriteTuple([&](int c) { return data[base + j + c * component_stride]; });
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
  const std::uint64_t payload_size = PointFieldPayloadSize(location, field);
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
                                                                  const std::string &name,
                                                                  const PointField &field,
                                                                  Vector *scratch)
{
  const std::uint64_t payload_size64 = PointFieldPayloadSize(location, field);
  MFEM_VERIFY(payload_size64 <= std::numeric_limits<std::uint32_t>::max(),
              PointFieldLocationName(location)
                  << " point field is too large for the current appended VTU writer!");
  const std::uint32_t payload_size = static_cast<std::uint32_t>(payload_size64);
  os.write(reinterpret_cast<const char *>(&payload_size), sizeof(payload_size));

  Vector values;
  const Vector *value_ptr = &GetPointFieldValues(location, field, values, scratch);
  WritePointFieldValues(location, os, ref, name, field, *value_ptr);
}

void CeedParaViewDataCollection::SaveDataVTU(std::ostream &os, int ref)
{
  // A provider cache is deliberately scoped to one synchronous VTU save. In particular,
  // solver vectors may retain their address while their values change between saves.
  ResetPointFieldProviders();
  os << "<VTKFile type=\"UnstructuredGrid\"";
  if (GetCompressionLevel() != 0)
  {
    os << " compressor=\"vtkZLibDataCompressor\"";
  }
  os << " version=\"2.2\" byte_order=\"" << mfem::VTKByteOrder() << "\">\n";
  os << "<UnstructuredGrid>\n";
  if (domain_cell_fields.empty())
  {
    if (cached_mesh_vtu.empty() || cached_mesh_vtu_ref != ref ||
        cached_mesh_vtu_bdr_output != bdr_output)
    {
      std::ostringstream mesh_os;
      mesh_os.precision(os.precision());
      mesh->PrintVTU(mesh_os, ref, pv_data_format, high_order_output, GetCompressionLevel(),
                     bdr_output);
      cached_mesh_vtu = mesh_os.str();
      cached_mesh_vtu_ref = ref;
      cached_mesh_vtu_bdr_output = bdr_output;
    }
    os.write(cached_mesh_vtu.data(), static_cast<std::streamsize>(cached_mesh_vtu.size()));
  }
  else
  {
    // MFEM owns mesh/attribute VTU emission and closes the <CellData> block before
    // returning. For lightweight element fields (rank, AMR indicator), splice our arrays
    // into that existing CellData section instead of routing them through a piecewise
    // constant GridFunction, which would expand one cell value to every refined point.
    std::ostringstream mesh_os;
    mesh_os.precision(os.precision());
    mesh->PrintVTU(mesh_os, ref, pv_data_format, high_order_output, GetCompressionLevel(),
                   bdr_output);
    const std::string mesh_vtu = mesh_os.str();
    const std::string cell_data_end = "</CellData>";
    const auto pos = mesh_vtu.find(cell_data_end);
    MFEM_VERIFY(pos != std::string::npos,
                "Unable to locate CellData block in MFEM VTU mesh output!");
    os.write(mesh_vtu.data(), static_cast<std::streamsize>(pos));
    for (const auto &kv : domain_cell_fields)
    {
      WriteCellFieldVTU(os, kv.first, kv.second);
    }
    os.write(mesh_vtu.data() + pos, static_cast<std::streamsize>(mesh_vtu.size() - pos));
  }

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
  const bool cache_appended_headers =
      (appended_boundary_fields || appended_domain_fields) &&
      (boundary_point_fields.empty() || appended_boundary_fields) &&
      (domain_point_fields.empty() || appended_domain_fields);
  auto WriteAppendedPointHeaders = [&](std::ostream &header_os)
  {
    std::uint64_t appended_offset = 0;
    for (const auto &kv : boundary_point_fields)
    {
      SavePointFieldVTUAppendedHeader(MeshEntityType::Boundary, header_os, ref, kv.first,
                                      kv.second, appended_offset);
      appended_offset += sizeof(std::uint32_t) +
                         PointFieldPayloadSize(MeshEntityType::Boundary, kv.second);
    }
    for (const auto &kv : domain_point_fields)
    {
      SavePointFieldVTUAppendedHeader(MeshEntityType::Domain, header_os, ref, kv.first,
                                      kv.second, appended_offset);
      appended_offset +=
          sizeof(std::uint32_t) + PointFieldPayloadSize(MeshEntityType::Domain, kv.second);
    }
  };
  if (cache_appended_headers)
  {
    if (cached_point_headers_vtu.empty() || cached_point_headers_vtu_ref != ref ||
        cached_point_headers_vtu_bdr_output != bdr_output)
    {
      std::ostringstream header_os;
      WriteAppendedPointHeaders(header_os);
      cached_point_headers_vtu = header_os.str();
      cached_point_headers_vtu_ref = ref;
      cached_point_headers_vtu_bdr_output = bdr_output;
    }
    os.write(cached_point_headers_vtu.data(),
             static_cast<std::streamsize>(cached_point_headers_vtu.size()));
  }
  else
  {
    std::uint64_t appended_offset = 0;
    for (const auto &kv : boundary_point_fields)
    {
      if (appended_boundary_fields)
      {
        SavePointFieldVTUAppendedHeader(MeshEntityType::Boundary, os, ref, kv.first,
                                        kv.second, appended_offset);
        appended_offset += sizeof(std::uint32_t) +
                           PointFieldPayloadSize(MeshEntityType::Boundary, kv.second);
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
        SavePointFieldVTUAppendedHeader(MeshEntityType::Domain, os, ref, kv.first,
                                        kv.second, appended_offset);
        appended_offset += sizeof(std::uint32_t) +
                           PointFieldPayloadSize(MeshEntityType::Domain, kv.second);
      }
      else
      {
        SavePointFieldVTU(MeshEntityType::Domain, os, ref, kv.first, kv.second,
                          &domain_point_scratch);
      }
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
      SavePointFieldVTUAppendedPayload(MeshEntityType::Boundary, os, ref, kv.first,
                                       kv.second, &boundary_point_scratch);
    }
    for (const auto &kv : domain_point_fields)
    {
      SavePointFieldVTUAppendedPayload(MeshEntityType::Domain, os, ref, kv.first, kv.second,
                                       &domain_point_scratch);
    }
    os << "\n</AppendedData>\n";
  }
  os << "</VTKFile>" << std::endl;
}

long long CeedParaViewDataCollection::IndividualEvaluatorCallCount()
{
  return individual_evaluator_calls.load();
}

long long CeedParaViewDataCollection::ProviderCallCount()
{
  return provider_calls.load();
}

long long CeedParaViewDataCollection::BatchViewFieldCount()
{
  return batch_view_fields.load();
}

long long CeedParaViewDataCollection::AvoidedSliceCopyCount()
{
  return avoided_slice_copies.load();
}

long long CeedParaViewDataCollection::AvoidedSliceCopyBytes()
{
  return avoided_slice_copy_bytes.load();
}

void CeedParaViewDataCollection::PrintPointFieldProfile() const
{
  if (!std::getenv("PALACE_CEED_PARAVIEW_PROFILE"))
  {
    return;
  }
  Mpi::Print("CeedParaViewDataCollection profile individual_evaluator_calls={} "
             "provider_calls={} batch_view_fields={} avoided_slice_copies={} "
             "avoided_slice_copy_bytes={}\n",
             IndividualEvaluatorCallCount(), ProviderCallCount(), BatchViewFieldCount(),
             AvoidedSliceCopyCount(), AvoidedSliceCopyBytes());
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
  PrintPointFieldProfile();

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
      for (const auto &field_it : domain_cell_fields)
      {
        pvtu_out << "<PDataArray type=\"" << GetDataTypeString() << "\" Name=\""
                 << field_it.first << "\" NumberOfComponents=\"" << field_it.second.num_comp
                 << "\" " << mfem::VTKComponentLabels(field_it.second.num_comp) << " "
                 << "format=\"" << GetDataFormatString() << "\" />\n";
      }
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
