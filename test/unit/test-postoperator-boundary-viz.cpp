// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>
#include <mfem/general/binaryio.hpp>
#ifdef MFEM_USE_ZLIB
#include <zlib.h>
#endif
#include <catch2/catch_test_macros.hpp>
#include "fem/coefficient.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

namespace fs = std::filesystem;

std::unique_ptr<mfem::Mesh> MakeSmallTetInterfaceSerialMesh()
{
  // 2 x 2 x 2 cubes split into tetrahedra: small enough for a unit test, but it has
  // exterior boundaries, an interior material interface, and enough face orientations to
  // exercise boundary-normal conventions through PostOperator's boundary ParaView output.
  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
  for (int e = 0; e < smesh->GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh->GetElementCenter(e, center);
    smesh->SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }

  // Add boundary elements on the interior material interface z = 0.5. The public
  // PostOperator output path should preserve the legacy one-sided exterior and two-sided
  // interior boundary semantics when its internals are replaced by libCEED kernels.
  for (int f = 0; f < smesh->GetNumFaces(); f++)
  {
    int e1, e2;
    smesh->GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh->GetAttribute(e1) != smesh->GetAttribute(e2))
    {
      auto *face_elem = smesh->GetFace(f)->Duplicate(smesh.get());
      face_elem->SetAttribute(7);
      smesh->AddBdrElement(face_elem);
    }
  }
  smesh->FinalizeTopology();
  smesh->Finalize();
  smesh->SetAttributes();
  smesh->EnsureNodes();
  return smesh;
}

std::unique_ptr<mfem::Mesh> MakeSmallTriInterfaceSerialMesh()
{
  // 2D companion to MakeSmallTetInterfaceSerialMesh. It has exterior boundaries and an
  // interior material interface so boundary point fields exercise line-boundary kernels,
  // scalar out-of-plane B, and two-sided averaging through the public PostOperator path.
  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(3, 2, mfem::Element::TRIANGLE, false, 1.0, 1.0));
  for (int e = 0; e < smesh->GetNE(); e++)
  {
    mfem::Vector center(2);
    smesh->GetElementCenter(e, center);
    smesh->SetAttribute(e, (center(1) < 0.5) ? 1 : 2);
  }

  for (int f = 0; f < smesh->GetNumFaces(); f++)
  {
    int e1, e2;
    smesh->GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh->GetAttribute(e1) != smesh->GetAttribute(e2))
    {
      auto *face_elem = smesh->GetFace(f)->Duplicate(smesh.get());
      face_elem->SetAttribute(7);
      smesh->AddBdrElement(face_elem);
    }
  }
  smesh->FinalizeTopology();
  smesh->Finalize();
  smesh->SetAttributes();
  smesh->EnsureNodes();
  return smesh;
}

std::vector<config::MaterialData> MakeTwoMaterials()
{
  config::MaterialData lower, upper;
  lower.attributes = {1};
  upper.attributes = {2};

  // Diagonal anisotropy and different values on each side exercise component ordering and
  // material lookup, not just identity/isotropic paths.
  lower.epsilon_r.s = {1.1, 1.2, 1.3};
  lower.mu_r.s = {0.9, 1.0, 1.1};
  upper.epsilon_r.s = {11.7, 3.1, 2.4};
  upper.mu_r.s = {1.4, 1.8, 2.2};

  return {lower, upper};
}

void ProjectTestFields3D(GridFunction &E, GridFunction &B)
{
  mfem::VectorFunctionCoefficient er(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(1.3 * x(1)) + x(2) * x(2) + 0.2;
                                       v(1) = std::cos(0.7 * x(2)) + 0.4 * x(0);
                                       v(2) = x(0) * x(1) - 0.3 * x(2) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient ei(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) * x(2) - 0.5;
                                       v(1) = std::sin(x(0)) - 0.25 * x(2);
                                       v(2) = std::cos(1.1 * x(1)) + x(0) * x(0);
                                     });
  mfem::VectorFunctionCoefficient br(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) - 0.3 * x(2) + 0.4;
                                       v(1) = std::sin(0.9 * x(2)) + 0.5;
                                       v(2) = std::cos(x(0)) - x(1) * x(2);
                                     });
  mfem::VectorFunctionCoefficient bi(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::cos(x(2)) - 0.2;
                                       v(1) = x(0) * x(2) + 0.1;
                                       v(2) = std::sin(x(1)) - 0.4 * x(0);
                                     });

  E.Real().ProjectCoefficient(er);
  E.Imag().ProjectCoefficient(ei);
  B.Real().ProjectCoefficient(br);
  B.Imag().ProjectCoefficient(bi);

  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();
  B.Imag().ExchangeFaceNbrData();
}

void ProjectTestFields2D(GridFunction &E, GridFunction &B)
{
  mfem::VectorFunctionCoefficient er(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(1.7 * x(1)) + 0.25 * x(0) + 0.1;
                                       v(1) = std::cos(0.8 * x(0)) - 0.35 * x(1) + 0.2;
                                     });
  mfem::VectorFunctionCoefficient ei(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(0) * x(1) - 0.15;
                                       v(1) = std::sin(0.9 * x(0)) + 0.4 * x(1);
                                     });
  mfem::FunctionCoefficient br([](const mfem::Vector &x)
                               { return std::cos(1.1 * x(0)) + 0.3 * x(1) + 0.5; });
  mfem::FunctionCoefficient bi([](const mfem::Vector &x)
                               { return std::sin(1.4 * x(1)) - 0.2 * x(0) + 0.25; });

  E.Real().ProjectCoefficient(er);
  E.Imag().ProjectCoefficient(ei);
  B.Real().ProjectCoefficient(br);
  B.Imag().ProjectCoefficient(bi);

  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();
  B.Imag().ExchangeFaceNbrData();
}

void GetTrueDofs(const GridFunction &gf, ComplexVector &x)
{
  gf.Real().GetTrueDofs(x.Real());
  gf.Imag().GetTrueDofs(x.Imag());
}

void DimensionalizeForPostOperatorOutput(const Units &units, GridFunction &E,
                                         GridFunction &B)
{
  // Reproduce the temporary mesh and Piola scaling performed by WriteParaviewFields so
  // legacy coefficient references and emitted Points are compared in input coordinates.
  auto &mesh = *E.ParFESpace()->GetParMesh();
  const double L = units.GetMeshLengthRelativeScale();
  mesh::DimensionalizeMesh(mesh, L);

  const double e_scale = L * units.GetScaleFactor<Units::ValueType::FIELD_E>();
  E *= e_scale;
  E.Real().FaceNbrData() *= e_scale;
  E.Imag().FaceNbrData() *= e_scale;

  const double b_scale =
      std::pow(L, mesh.Dimension() - 1) * units.GetScaleFactor<Units::ValueType::FIELD_B>();
  B *= b_scale;
  B.Real().FaceNbrData() *= b_scale;
  B.Imag().FaceNbrData() *= b_scale;
}

struct ErrorStats
{
  long long count = 0;
  double max_abs = 0.0;
  double max_scaled = 0.0;
  double sum_sq = 0.0;

  void Add(double val, double ref, double rtol, double atol)
  {
    const double err = std::abs(val - ref);
    const double denom = atol + rtol * std::max({1.0, std::abs(val), std::abs(ref)});
    max_abs = std::max(max_abs, err);
    max_scaled = std::max(max_scaled, err / denom);
    sum_sq += err * err;
    count++;
  }

  double Rms() const
  {
    return count ? std::sqrt(sum_sq / static_cast<double>(count)) : 0.0;
  }
};

void CheckStats(const std::string &name, const ErrorStats &stats)
{
  INFO("field = " << name);
  INFO("count = " << stats.count);
  INFO("max_abs = " << stats.max_abs);
  INFO("max_scaled = " << stats.max_scaled);
  INFO("rms = " << stats.Rms());
  CHECK(stats.count > 0);
  CHECK(stats.max_scaled <= 1.0);
}

std::string ReadFile(const fs::path &path)
{
  std::ifstream in(path, std::ios::binary);
  REQUIRE(in.good());
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

std::string RemoveWhitespace(std::string_view text)
{
  std::string out;
  out.reserve(text.size());
  for (char c : text)
  {
    if (!std::isspace(static_cast<unsigned char>(c)))
    {
      out.push_back(c);
    }
  }
  return out;
}

std::vector<char> DecodeBase64(std::string_view encoded)
{
  std::vector<char> decoded;
  const auto clean = RemoveWhitespace(encoded);
  mfem::bin_io::DecodeBase64(clean.c_str(), clean.size(), decoded);
  return decoded;
}

template <typename T>
T ReadLittleEndian(const char *bytes)
{
  T value;
  std::memcpy(&value, bytes, sizeof(T));
  return value;
}

std::vector<char> DecodeVtkBinaryPayload(std::string_view encoded, bool compressed)
{
  const auto clean = RemoveWhitespace(encoded);
  if (!compressed)
  {
    const auto decoded = DecodeBase64(clean);
    REQUIRE(decoded.size() >= sizeof(std::uint32_t));
    const auto nbytes = ReadLittleEndian<std::uint32_t>(decoded.data());
    REQUIRE(decoded.size() >= sizeof(std::uint32_t) + nbytes);
    return {decoded.begin() + sizeof(std::uint32_t),
            decoded.begin() + sizeof(std::uint32_t) + nbytes};
  }

  // MFEM writes one compressed VTK block with a 4 * uint32 header encoded separately
  // from the zlib payload: [num_blocks=1, uncompressed_size, partial_size=0,
  // compressed_size]. Decode the two base64 chunks independently because the header chunk
  // carries its own padding.
  constexpr std::size_t header_bytes = 4 * sizeof(std::uint32_t);
  const auto header_chars = mfem::bin_io::NumBase64Chars(header_bytes);
  REQUIRE(clean.size() >= header_chars);
  const auto header = DecodeBase64(std::string_view(clean).substr(0, header_chars));
  REQUIRE(header.size() == header_bytes);
  const auto num_blocks = ReadLittleEndian<std::uint32_t>(header.data());
  const auto uncompressed_size = ReadLittleEndian<std::uint32_t>(header.data() + 4);
  const auto partial_size = ReadLittleEndian<std::uint32_t>(header.data() + 8);
  const auto compressed_size = ReadLittleEndian<std::uint32_t>(header.data() + 12);
  REQUIRE(num_blocks == 1);
  REQUIRE(partial_size == 0);

  const auto compressed_payload =
      DecodeBase64(std::string_view(clean).substr(header_chars));
  REQUIRE(compressed_payload.size() >= compressed_size);
  std::vector<char> uncompressed(uncompressed_size);
#ifdef MFEM_USE_ZLIB
  uLongf dest_len = uncompressed_size;
  const int zret = uncompress(reinterpret_cast<Bytef *>(uncompressed.data()), &dest_len,
                              reinterpret_cast<const Bytef *>(compressed_payload.data()),
                              compressed_size);
  REQUIRE(zret == Z_OK);
  REQUIRE(dest_len == uncompressed_size);
#else
  MFEM_ABORT("Cannot decode compressed VTK output without zlib support!");
#endif
  return uncompressed;
}

std::string_view ExtractDataArrayTag(std::string_view xml, const std::string &name)
{
  const std::string needle = "Name=\"" + name + "\"";
  const auto name_pos = xml.find(needle);
  REQUIRE(name_pos != std::string_view::npos);
  const auto tag_begin = xml.rfind("<DataArray", name_pos);
  REQUIRE(tag_begin != std::string_view::npos);
  const auto tag_end = xml.find('>', name_pos);
  REQUIRE(tag_end != std::string_view::npos);
  return xml.substr(tag_begin, tag_end - tag_begin + 1);
}

std::string_view ExtractDataArrayPayload(std::string_view xml, const std::string &name)
{
  const std::string needle = "Name=\"" + name + "\"";
  const auto name_pos = xml.find(needle);
  REQUIRE(name_pos != std::string_view::npos);
  const auto tag_end = xml.find('>', name_pos);
  REQUIRE(tag_end != std::string_view::npos);
  const auto close = xml.find("</DataArray>", tag_end);
  REQUIRE(close != std::string_view::npos);
  return xml.substr(tag_end + 1, close - tag_end - 1);
}

std::uint64_t ExtractUnsignedAttribute(std::string_view tag, std::string_view attr)
{
  const std::string needle = std::string(attr) + "=\"";
  const auto begin = tag.find(needle);
  REQUIRE(begin != std::string_view::npos);
  const auto value_begin = begin + needle.size();
  const auto value_end = tag.find('"', value_begin);
  REQUIRE(value_end != std::string_view::npos);
  return std::stoull(std::string(tag.substr(value_begin, value_end - value_begin)));
}

std::vector<char> ExtractAppendedDataArrayPayload(std::string_view xml,
                                                  std::string_view tag)
{
  REQUIRE(tag.find("format=\"appended\"") != std::string_view::npos);
  const auto offset = ExtractUnsignedAttribute(tag, "offset");
  const auto appended_begin = xml.find("<AppendedData");
  REQUIRE(appended_begin != std::string_view::npos);
  const auto underscore = xml.find('_', appended_begin);
  REQUIRE(underscore != std::string_view::npos);
  const auto data_begin = underscore + 1 + offset;
  REQUIRE(data_begin + sizeof(std::uint32_t) <= xml.size());
  const auto nbytes = ReadLittleEndian<std::uint32_t>(xml.data() + data_begin);
  const auto payload_begin = data_begin + sizeof(std::uint32_t);
  REQUIRE(payload_begin + nbytes <= xml.size());
  return {xml.begin() + static_cast<std::ptrdiff_t>(payload_begin),
          xml.begin() + static_cast<std::ptrdiff_t>(payload_begin + nbytes)};
}

bool VtuUsesCompression(std::string_view xml)
{
  return xml.find("compressor=\"vtkZLibDataCompressor\"") != std::string_view::npos;
}

std::vector<char> ReadDataArrayBytes(std::string_view xml, const std::string &name)
{
  const auto tag = ExtractDataArrayTag(xml, name);
  if (tag.find("format=\"appended\"") != std::string_view::npos)
  {
    return ExtractAppendedDataArrayPayload(xml, tag);
  }
  return DecodeVtkBinaryPayload(ExtractDataArrayPayload(xml, name),
                                VtuUsesCompression(xml));
}

std::vector<double> DecodeFloatDataArray(std::string_view tag,
                                         const std::vector<char> &bytes)
{
  std::vector<double> values;
  if (tag.find("type=\"Float32\"") != std::string_view::npos)
  {
    REQUIRE(bytes.size() % sizeof(float) == 0);
    values.resize(bytes.size() / sizeof(float));
    for (std::size_t i = 0; i < values.size(); i++)
    {
      values[i] = ReadLittleEndian<float>(bytes.data() + i * sizeof(float));
    }
  }
  else
  {
    REQUIRE(tag.find("type=\"Float64\"") != std::string_view::npos);
    REQUIRE(bytes.size() % sizeof(double) == 0);
    values.resize(bytes.size() / sizeof(double));
    for (std::size_t i = 0; i < values.size(); i++)
    {
      values[i] = ReadLittleEndian<double>(bytes.data() + i * sizeof(double));
    }
  }
  return values;
}

std::vector<double> ReadFloatDataArray(std::string_view xml, const std::string &name)
{
  const auto tag = ExtractDataArrayTag(xml, name);
  return DecodeFloatDataArray(tag, ReadDataArrayBytes(xml, name));
}

std::vector<double> ReadPointsDataArray(std::string_view xml)
{
  const auto points_begin = xml.find("<Points>");
  REQUIRE(points_begin != std::string_view::npos);
  const auto tag_begin = xml.find("<DataArray", points_begin);
  REQUIRE(tag_begin != std::string_view::npos);
  const auto tag_end = xml.find('>', tag_begin);
  REQUIRE(tag_end != std::string_view::npos);
  const auto close = xml.find("</DataArray>", tag_end);
  REQUIRE(close != std::string_view::npos);
  const auto tag = xml.substr(tag_begin, tag_end - tag_begin + 1);
  REQUIRE(tag.find("format=\"appended\"") == std::string_view::npos);
  const auto payload = xml.substr(tag_end + 1, close - tag_end - 1);
  return DecodeFloatDataArray(tag,
                              DecodeVtkBinaryPayload(payload, VtuUsesCompression(xml)));
}

std::vector<int> ReadIntDataArray(std::string_view xml, const std::string &name)
{
  const auto bytes = ReadDataArrayBytes(xml, name);
  REQUIRE(bytes.size() % sizeof(std::int32_t) == 0);
  std::vector<int> values(bytes.size() / sizeof(std::int32_t));
  for (std::size_t i = 0; i < values.size(); i++)
  {
    values[i] = ReadLittleEndian<std::int32_t>(bytes.data() + i * sizeof(std::int32_t));
  }
  return values;
}

int CountBoundaryVisualizationPoints(const mfem::ParMesh &pmesh, int lod)
{
  int count = 0;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementBaseGeometry(i), lod, 1);
    count += RefG.RefPts.GetNPoints();
  }
  return count;
}

int CountDomainVisualizationPoints(const mfem::ParMesh &pmesh, int lod)
{
  int count = 0;
  for (int i = 0; i < pmesh.GetNE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetElementBaseGeometry(i), lod, 1);
    count += RefG.RefPts.GetNPoints();
  }
  return count;
}

ErrorStats CompareMeshPoints(const std::vector<double> &values, mfem::ParMesh &pmesh,
                             int lod, bool boundary, double rtol, double atol)
{
  ErrorStats stats;
  int idx = 0;
  const int num_entity = boundary ? pmesh.GetNBE() : pmesh.GetNE();
  mfem::Vector x(pmesh.SpaceDimension());
  for (int i = 0; i < num_entity; i++)
  {
    const auto geom =
        boundary ? pmesh.GetBdrElementBaseGeometry(i) : pmesh.GetElementBaseGeometry(i);
    const auto &RefG = *mfem::GlobGeometryRefiner.Refine(geom, lod, 1);
    auto *T =
        boundary ? pmesh.GetBdrElementTransformation(i) : pmesh.GetElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
    {
      T->Transform(RefG.RefPts.IntPoint(j), x);
      for (int c = 0; c < 3; c++, idx++)
      {
        stats.Add(values[idx], c < x.Size() ? x(c) : 0.0, rtol, atol);
      }
    }
  }
  REQUIRE(idx == static_cast<int>(values.size()));
  return stats;
}

ErrorStats CompareScalarField(const std::vector<double> &values, mfem::ParMesh &pmesh,
                              int lod, mfem::Coefficient &legacy, double rtol, double atol)
{
  ErrorStats stats;
  int idx = 0;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementBaseGeometry(i), lod, 1);
    auto *T = pmesh.GetBdrElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++, idx++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      stats.Add(values[idx], legacy.Eval(*T, ip), rtol, atol);
    }
  }
  REQUIRE(idx == static_cast<int>(values.size()));
  return stats;
}

ErrorStats CompareScalarDomainField(const std::vector<double> &values, mfem::ParMesh &pmesh,
                                    int lod, mfem::Coefficient &legacy, double rtol,
                                    double atol)
{
  ErrorStats stats;
  int idx = 0;
  for (int i = 0; i < pmesh.GetNE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetElementBaseGeometry(i), lod, 1);
    auto *T = pmesh.GetElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++, idx++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      stats.Add(values[idx], legacy.Eval(*T, ip), rtol, atol);
    }
  }
  REQUIRE(idx == static_cast<int>(values.size()));
  return stats;
}

ErrorStats CompareVectorField(const std::vector<double> &values, mfem::ParMesh &pmesh,
                              int lod, mfem::VectorCoefficient &legacy, double rtol,
                              double atol)
{
  ErrorStats stats;
  int idx = 0;
  mfem::Vector ref(legacy.GetVDim());
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementBaseGeometry(i), lod, 1);
    auto *T = pmesh.GetBdrElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      legacy.Eval(ref, *T, ip);
      for (int c = 0; c < legacy.GetVDim(); c++, idx++)
      {
        stats.Add(values[idx], ref(c), rtol, atol);
      }
    }
  }
  REQUIRE(idx == static_cast<int>(values.size()));
  return stats;
}

void RequireInteriorBoundaryWasWritten(std::string_view xml, const mfem::ParMesh &pmesh)
{
  const auto attributes = ReadIntDataArray(xml, "attribute");
  REQUIRE(attributes.size() == static_cast<std::size_t>(pmesh.GetNBE()));
  CHECK(std::find(attributes.begin(), attributes.end(), 7) != attributes.end());
}

ErrorStats CompareVectorDomainField(const std::vector<double> &values, mfem::ParMesh &pmesh,
                                    int lod, mfem::VectorCoefficient &legacy, double rtol,
                                    double atol)
{
  ErrorStats stats;
  int idx = 0;
  mfem::Vector ref(legacy.GetVDim());
  for (int i = 0; i < pmesh.GetNE(); i++)
  {
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetElementBaseGeometry(i), lod, 1);
    auto *T = pmesh.GetElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      legacy.Eval(ref, *T, ip);
      for (int c = 0; c < legacy.GetVDim(); c++, idx++)
      {
        stats.Add(values[idx], ref(c), rtol, atol);
      }
    }
  }
  REQUIRE(idx == static_cast<int>(values.size()));
  return stats;
}

void RequirePointFieldUsesAppendedData(std::string_view xml, const std::string &name)
{
  const auto tag = ExtractDataArrayTag(xml, name);
  CHECK(tag.find("format=\"appended\"") != std::string_view::npos);
  CHECK(xml.find("<AppendedData encoding=\"raw\">") != std::string_view::npos);
}

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir,
                 "PostOperator boundary ParaView fields match legacy coefficients",
                 "[postoperator][boundary-viz][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  REQUIRE(Mpi::Size(comm) == 1);

  constexpr int order = 2;
  // PostOperator intentionally writes ParaView files as Float32, so tolerate single
  // precision roundoff in the file representation while still catching semantic errors.
  constexpr double rtol = 5.0e-5;
  constexpr double atol = 5.0e-7;

  auto serial_mesh = MakeSmallTetInterfaceSerialMesh();
  REQUIRE(serial_mesh->bdr_attributes.Find(7) >= 0);

  IoData iodata(Units(1.0e-4, 2.5e-4));
  iodata.problem.type = ProblemType::EIGENMODE;
  iodata.problem.verbose = 0;
  iodata.problem.output = temp_dir.string();
  iodata.problem.output_formats.paraview = true;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0e-4;
  iodata.model.Lc = 2.5;
  iodata.domains.materials = MakeTwoMaterials();
  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 5, 6};
  iodata.solver.order = order;
  iodata.solver.eigenmode.n_post = 1;
  iodata.CheckConfiguration();
  iodata.NondimensionalizeInputs(serial_mesh);

  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  auto &pmesh = mesh.front()->Get();
  REQUIRE(pmesh.bdr_attributes.Find(7) >= 0);

  SpaceOperator space_op(iodata, mesh);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);

  GridFunction E(space_op.GetNDSpace(), true), B(space_op.GetCurlSpace(), true);
  ProjectTestFields3D(E, B);
  ComplexVector e_true(space_op.GetNDSpace().GetTrueVSize());
  ComplexVector b_true(space_op.GetRTSpace().GetTrueVSize());
  GetTrueDofs(E, e_true);
  GetTrueDofs(B, b_true);

  // Public API under test: this is the same PostOperator entry point used by the solver.
  // On origin/main it writes boundary ParaView fields by evaluating legacy coefficients;
  // this branch replaces selected internals with libCEED buffers while preserving the API
  // and output contract.
  post_op.MeasureAndPrintAll(0, e_true, b_true, {1.0, 0.1}, 0.0, 0.0, 1);

  DimensionalizeForPostOperatorOutput(iodata.units, E, B);

  const auto domain_vtu = ReadFile(fs::path(iodata.problem.output) / "paraview" /
                                   "eigenmode" / "Cycle000001" / "proc000000.vtu");
  const auto boundary_vtu =
      ReadFile(fs::path(iodata.problem.output) / "paraview" / "eigenmode_boundary" /
               "Cycle000001" / "proc000000.vtu");
  RequireInteriorBoundaryWasWritten(boundary_vtu, pmesh);
  RequirePointFieldUsesAppendedData(boundary_vtu, "E_real");
  RequirePointFieldUsesAppendedData(domain_vtu, "U_e");

  const int lod = order;
  const int n_bdr_pts = CountBoundaryVisualizationPoints(pmesh, lod);
  const int n_domain_pts = CountDomainVisualizationPoints(pmesh, lod);
  auto ReadBoundaryChecked = [&](const std::string &name, int components)
  {
    auto values = ReadFloatDataArray(boundary_vtu, name);
    REQUIRE(values.size() == static_cast<std::size_t>(n_bdr_pts * components));
    return values;
  };
  auto ReadDomainChecked = [&](const std::string &name, int components)
  {
    auto values = ReadFloatDataArray(domain_vtu, name);
    REQUIRE(values.size() == static_cast<std::size_t>(n_domain_pts * components));
    return values;
  };

  const auto &mat_op = space_op.GetMaterialOp();
  const double eps_scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) /
                             iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  const double invmu_scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
                               iodata.units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);

  CheckStats("domain mesh points", CompareMeshPoints(ReadPointsDataArray(domain_vtu), pmesh,
                                                     lod, false, rtol, atol));
  CheckStats("boundary mesh points", CompareMeshPoints(ReadPointsDataArray(boundary_vtu),
                                                       pmesh, lod, true, rtol, atol));

  BdrFieldVectorCoefficient E_real_legacy(E.Real()), E_imag_legacy(E.Imag());
  BdrFieldVectorCoefficient B_real_legacy(B.Real()), B_imag_legacy(B.Imag());
  CheckStats("E_real", CompareVectorField(ReadBoundaryChecked("E_real", 3), pmesh, lod,
                                          E_real_legacy, rtol, atol));
  CheckStats("E_imag", CompareVectorField(ReadBoundaryChecked("E_imag", 3), pmesh, lod,
                                          E_imag_legacy, rtol, atol));
  CheckStats("B_real", CompareVectorField(ReadBoundaryChecked("B_real", 3), pmesh, lod,
                                          B_real_legacy, rtol, atol));
  CheckStats("B_imag", CompareVectorField(ReadBoundaryChecked("B_imag", 3), pmesh, lod,
                                          B_imag_legacy, rtol, atol));

  mfem::Vector unused_x0;
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_real_legacy(
      &E.Real(), nullptr, mat_op, true, unused_x0, eps_scaling);
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_imag_legacy(
      &E.Imag(), nullptr, mat_op, true, unused_x0, eps_scaling);
  BdrSurfaceCurrentVectorCoefficient J_real_legacy(B.Real(), mat_op, invmu_scaling);
  BdrSurfaceCurrentVectorCoefficient J_imag_legacy(B.Imag(), mat_op, invmu_scaling);
  CheckStats("Q_s_real", CompareScalarField(ReadBoundaryChecked("Q_s_real", 1), pmesh, lod,
                                            Q_real_legacy, rtol, atol));
  CheckStats("Q_s_imag", CompareScalarField(ReadBoundaryChecked("Q_s_imag", 1), pmesh, lod,
                                            Q_imag_legacy, rtol, atol));
  CheckStats("J_s_real", CompareVectorField(ReadBoundaryChecked("J_s_real", 3), pmesh, lod,
                                            J_real_legacy, rtol, atol));
  CheckStats("J_s_imag", CompareVectorField(ReadBoundaryChecked("J_s_imag", 3), pmesh, lod,
                                            J_imag_legacy, rtol, atol));

  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> Ue_legacy(E, mat_op, eps_scaling);
  EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> Um_legacy(B, mat_op, invmu_scaling);
  PoyntingVectorCoefficient S_legacy(E, B, mat_op, invmu_scaling);
  CheckStats("U_e boundary", CompareScalarField(ReadBoundaryChecked("U_e", 1), pmesh, lod,
                                                Ue_legacy, rtol, atol));
  CheckStats("U_m boundary", CompareScalarField(ReadBoundaryChecked("U_m", 1), pmesh, lod,
                                                Um_legacy, rtol, atol));
  CheckStats("S boundary", CompareVectorField(ReadBoundaryChecked("S", 3), pmesh, lod,
                                              S_legacy, rtol, atol));
  CheckStats("U_e domain", CompareScalarDomainField(ReadDomainChecked("U_e", 1), pmesh, lod,
                                                    Ue_legacy, rtol, atol));
  CheckStats("U_m domain", CompareScalarDomainField(ReadDomainChecked("U_m", 1), pmesh, lod,
                                                    Um_legacy, rtol, atol));
  CheckStats("S domain", CompareVectorDomainField(ReadDomainChecked("S", 3), pmesh, lod,
                                                  S_legacy, rtol, atol));
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "PostOperator 2D boundary ParaView fields match legacy coefficients",
                 "[postoperator][boundary-viz][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  REQUIRE(Mpi::Size(comm) == 1);

  constexpr int order = 2;
  constexpr double rtol = 5.0e-5;
  constexpr double atol = 5.0e-7;

  auto serial_mesh = MakeSmallTriInterfaceSerialMesh();
  REQUIRE(serial_mesh->Dimension() == 2);
  REQUIRE(serial_mesh->SpaceDimension() == 2);
  REQUIRE(serial_mesh->bdr_attributes.Find(7) >= 0);

  IoData iodata(Units(1.0e-4, 2.5e-4));
  iodata.problem.type = ProblemType::EIGENMODE;
  iodata.problem.verbose = 0;
  iodata.problem.output = temp_dir.string();
  iodata.problem.output_formats.paraview = true;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0e-4;
  iodata.model.Lc = 2.5;
  iodata.domains.materials = MakeTwoMaterials();
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};
  iodata.solver.order = order;
  iodata.solver.eigenmode.n_post = 1;
  iodata.CheckConfiguration();
  iodata.NondimensionalizeInputs(serial_mesh);

  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  auto &pmesh = mesh.front()->Get();
  REQUIRE(pmesh.bdr_attributes.Find(7) >= 0);

  SpaceOperator space_op(iodata, mesh);
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);

  GridFunction E(space_op.GetNDSpace(), true), B(space_op.GetCurlSpace(), true);
  REQUIRE(E.VectorDim() == 2);
  REQUIRE(B.VectorDim() == 1);
  ProjectTestFields2D(E, B);
  ComplexVector e_true(space_op.GetNDSpace().GetTrueVSize());
  ComplexVector b_true(space_op.GetCurlSpace().GetTrueVSize());
  GetTrueDofs(E, e_true);
  GetTrueDofs(B, b_true);

  post_op.MeasureAndPrintAll(0, e_true, b_true, {1.0, 0.1}, 0.0, 0.0, 1);

  DimensionalizeForPostOperatorOutput(iodata.units, E, B);

  const auto domain_vtu = ReadFile(fs::path(iodata.problem.output) / "paraview" /
                                   "eigenmode" / "Cycle000001" / "proc000000.vtu");
  const auto boundary_vtu =
      ReadFile(fs::path(iodata.problem.output) / "paraview" / "eigenmode_boundary" /
               "Cycle000001" / "proc000000.vtu");
  RequireInteriorBoundaryWasWritten(boundary_vtu, pmesh);
  RequirePointFieldUsesAppendedData(boundary_vtu, "E_real");
  RequirePointFieldUsesAppendedData(boundary_vtu, "U_m");
  RequirePointFieldUsesAppendedData(domain_vtu, "B_real");

  const int lod = order;
  const int n_bdr_pts = CountBoundaryVisualizationPoints(pmesh, lod);
  const int n_domain_pts = CountDomainVisualizationPoints(pmesh, lod);
  auto ReadBoundaryChecked = [&](const std::string &name, int components)
  {
    auto values = ReadFloatDataArray(boundary_vtu, name);
    REQUIRE(values.size() == static_cast<std::size_t>(n_bdr_pts * components));
    return values;
  };
  auto ReadDomainChecked = [&](const std::string &name, int components)
  {
    auto values = ReadFloatDataArray(domain_vtu, name);
    REQUIRE(values.size() == static_cast<std::size_t>(n_domain_pts * components));
    return values;
  };

  const auto &mat_op = space_op.GetMaterialOp();
  const double eps_scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) /
                             iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  const double invmu_scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
                               iodata.units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);

  CheckStats("domain mesh points 2D", CompareMeshPoints(ReadPointsDataArray(domain_vtu),
                                                        pmesh, lod, false, rtol, atol));
  CheckStats("boundary mesh points 2D", CompareMeshPoints(ReadPointsDataArray(boundary_vtu),
                                                          pmesh, lod, true, rtol, atol));

  BdrFieldVectorCoefficient E_real_legacy(E.Real()), E_imag_legacy(E.Imag());
  CheckStats("E_real boundary 2D",
             CompareVectorField(ReadBoundaryChecked("E_real", 2), pmesh, lod, E_real_legacy,
                                rtol, atol));
  CheckStats("E_imag boundary 2D",
             CompareVectorField(ReadBoundaryChecked("E_imag", 2), pmesh, lod, E_imag_legacy,
                                rtol, atol));

  mfem::Vector unused_x0;
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_real_legacy(
      &E.Real(), nullptr, mat_op, true, unused_x0, eps_scaling);
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_imag_legacy(
      &E.Imag(), nullptr, mat_op, true, unused_x0, eps_scaling);
  CheckStats("Q_s_real boundary 2D",
             CompareScalarField(ReadBoundaryChecked("Q_s_real", 1), pmesh, lod,
                                Q_real_legacy, rtol, atol));
  CheckStats("Q_s_imag boundary 2D",
             CompareScalarField(ReadBoundaryChecked("Q_s_imag", 1), pmesh, lod,
                                Q_imag_legacy, rtol, atol));

  EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> Ue_legacy(E, mat_op, eps_scaling);
  EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> Um_legacy(B, mat_op, invmu_scaling);
  PoyntingVectorCoefficient S_legacy(E, B, mat_op, invmu_scaling);
  CheckStats("U_e boundary 2D", CompareScalarField(ReadBoundaryChecked("U_e", 1), pmesh,
                                                   lod, Ue_legacy, rtol, atol));
  CheckStats("U_m boundary 2D", CompareScalarField(ReadBoundaryChecked("U_m", 1), pmesh,
                                                   lod, Um_legacy, rtol, atol));
  CheckStats("S boundary 2D", CompareVectorField(ReadBoundaryChecked("S", 2), pmesh, lod,
                                                 S_legacy, rtol, atol));

  mfem::VectorGridFunctionCoefficient E_real_domain(&E.Real());
  mfem::VectorGridFunctionCoefficient E_imag_domain(&E.Imag());
  mfem::GridFunctionCoefficient B_real_domain(&B.Real());
  mfem::GridFunctionCoefficient B_imag_domain(&B.Imag());
  CheckStats("E_real domain 2D",
             CompareVectorDomainField(ReadDomainChecked("E_real", 2), pmesh, lod,
                                      E_real_domain, rtol, atol));
  CheckStats("E_imag domain 2D",
             CompareVectorDomainField(ReadDomainChecked("E_imag", 2), pmesh, lod,
                                      E_imag_domain, rtol, atol));
  CheckStats("B_real domain 2D",
             CompareScalarDomainField(ReadDomainChecked("B_real", 1), pmesh, lod,
                                      B_real_domain, rtol, atol));
  CheckStats("B_imag domain 2D",
             CompareScalarDomainField(ReadDomainChecked("B_imag", 1), pmesh, lod,
                                      B_imag_domain, rtol, atol));
  CheckStats("U_e domain 2D", CompareScalarDomainField(ReadDomainChecked("U_e", 1), pmesh,
                                                       lod, Ue_legacy, rtol, atol));
  CheckStats("U_m domain 2D", CompareScalarDomainField(ReadDomainChecked("U_m", 1), pmesh,
                                                       lod, Um_legacy, rtol, atol));
  CheckStats("S domain 2D", CompareVectorDomainField(ReadDomainChecked("S", 2), pmesh, lod,
                                                     S_legacy, rtol, atol));
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "2D magnetostatic ParaView fields match the mesh point lattice",
                 "[postoperator][boundary-viz][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  REQUIRE(Mpi::Size(comm) == 1);

  constexpr int order = 2;
  auto serial_mesh = MakeSmallTriInterfaceSerialMesh();
  IoData iodata(Units(1.0e-4, 2.5e-4));
  iodata.problem.type = ProblemType::MAGNETOSTATIC;
  iodata.problem.verbose = 0;
  iodata.problem.output = temp_dir.string();
  iodata.problem.output_formats.paraview = true;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0e-4;
  iodata.model.Lc = 2.5;
  iodata.domains.materials = MakeTwoMaterials();
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};
  iodata.solver.order = order;
  iodata.solver.magnetostatic.n_post = 1;
  iodata.CheckConfiguration();
  iodata.NondimensionalizeInputs(serial_mesh);

  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  auto &pmesh = mesh.front()->Get();

  CurlCurlOperator curlcurl_op(iodata, mesh);
  PostOperator<ProblemType::MAGNETOSTATIC> post_op(iodata, curlcurl_op);
  const auto &Curl = curlcurl_op.GetCurlMatrix();
  Vector A(Curl.Width()), B(Curl.Height());
  A = 0.0;
  B = 0.0;
  post_op.MeasureAndPrintAll(0, A, B, 0);

  const auto domain_vtu = ReadFile(fs::path(iodata.problem.output) / "paraview" /
                                   "magnetostatic" / "Cycle000000" / "proc000000.vtu");
  const int num_points = CountDomainVisualizationPoints(pmesh, order);
  CHECK(ReadPointsDataArray(domain_vtu).size() == static_cast<std::size_t>(3 * num_points));
  CHECK(ReadFloatDataArray(domain_vtu, "A").size() ==
        static_cast<std::size_t>(2 * num_points));
  CHECK(ReadFloatDataArray(domain_vtu, "B").size() == static_cast<std::size_t>(num_points));
  CHECK(ReadFloatDataArray(domain_vtu, "U_m").size() ==
        static_cast<std::size_t>(num_points));
}

}  // namespace palace
