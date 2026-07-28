// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include "drivers/singularsolver.hpp"
#include "fem/singularfield.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/postoperator.hpp"
#include "models/singularsurfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Tangent frame of the extracted 2D submesh in the parent 3D coordinates. Used by
// Preprocess to rotate material tensors and project path coordinates into the local
// frame.
struct SubmeshFrame
{
  mfem::Vector centroid{0.0, 0.0, 0.0};
  mfem::Vector e1{0.0, 0.0, 0.0};
  mfem::Vector e2{0.0, 0.0, 0.0};
  mfem::Vector normal{0.0, 0.0, 0.0};
};

// Full 3D-boundary → 2D-submesh pipeline on the pre-partitioned serial mesh:
// CreateFromBoundary → attribute remap → add internal-boundary edges → 3D→2D
// projection. Runs on ranks that hold the serial mesh; the returned 2D mesh has its
// own node coordinates and no longer references the parent.
std::unique_ptr<mfem::Mesh>
ExtractBoundary2DSubmesh(mfem::Mesh &parent, const mfem::Array<int> &surface_attrs,
                         const std::vector<int> &internal_bdr_attrs, SubmeshFrame &frame)
{
  auto sub = std::make_unique<mfem::SubMesh>(
      mfem::SubMesh::CreateFromBoundary(parent, surface_attrs));
  // Project to 2D before remapping element attributes: ProjectSubmeshTo2D's call to
  // GetSurfaceNormal walks the mesh's unique-attribute list, which RemapSubMeshAttributes
  // leaves stale relative to the per-element attributes it just rewrote.
  frame.normal = mesh::ProjectSubmeshTo2D(*sub, frame.centroid, frame.e1, frame.e2);
  mesh::RemapSubMeshAttributes(*sub);
  mesh::RemapSubMeshBdrAttributes(*sub, surface_attrs);
  mesh::AddSubMeshInternalBoundaryElements(*sub, surface_attrs, internal_bdr_attrs);
  return sub;
}

// Project each 3D point in `path` onto the 2D local frame (centroid, e1, e2). Leaves
// entries that are already 2D untouched. `path` is the iodata-side vector-of-vectors
// representation; we write back 2-element entries so the downstream PostOperator reads
// them as 2D coordinates.
void ProjectPathTo2D(std::vector<std::vector<double>> &path, const mfem::Vector &centroid,
                     const mfem::Vector &e1, const mfem::Vector &e2)
{
  for (auto &p : path)
  {
    if (p.size() != 3)
    {
      continue;
    }
    mfem::Vector p3d(p.data(), 3);
    mfem::Vector p2d = mesh::Project3Dto2D(p3d, centroid, e1, e2);
    p.assign({p2d(0), p2d(1)});
  }
}

struct SingularModeMeasurement
{
  int mode;
  std::complex<double> propagation_constant;
  std::complex<double> effective_index;
  double backward_error;
  double absolute_error;
  std::complex<double> normalized_power;
  BoundaryModeOperator::SingularCoefficientNorms coefficient_norms;
  BoundaryModeOperator::SingularFieldEnergies field_energies;
};

enum class SingularModeSpace
{
  ND,
  H1
};

struct SingularModeCoefficientMeasurement
{
  int mode;
  SingularModeSpace space;
  HYPRE_BigInt true_dof;
  fem::singular::DofKey key;
  double exponent;
  std::complex<double> coefficient;
};

struct SingularModeTipSlopeMeasurement
{
  int mode;
  fem::singular::H1TipSlopeDiagnostic diagnostic;
};

struct SingularModeSurfaceMeasurement
{
  int mode;
  std::vector<TriangleSingularSurfacePostOperator::Measurement> interfaces;
};

constexpr std::size_t ModeCoefficientIntegerFields = 24;
constexpr std::size_t ModeCoefficientRealFields = 3;
using PackedModeCoefficientIntegers =
    std::array<fem::singular::GlobalVertexId, ModeCoefficientIntegerFields>;
using PackedModeCoefficientReals = std::array<double, ModeCoefficientRealFields>;
constexpr std::size_t ModeTipSlopeIntegerFields = 9;
constexpr std::size_t ModeTipSlopeRealFields = 8;
using PackedModeTipSlopeIntegers =
    std::array<fem::singular::GlobalVertexId, ModeTipSlopeIntegerFields>;
using PackedModeTipSlopeReals = std::array<double, ModeTipSlopeRealFields>;

fem::singular::GlobalVertexId ToDiagnosticInteger(HYPRE_BigInt value)
{
  static_assert(std::numeric_limits<HYPRE_BigInt>::is_signed);
  static_assert(std::numeric_limits<fem::singular::GlobalVertexId>::is_signed);
  static_assert(std::numeric_limits<HYPRE_BigInt>::digits <=
                std::numeric_limits<fem::singular::GlobalVertexId>::digits);
  return static_cast<fem::singular::GlobalVertexId>(value);
}

fem::singular::GlobalVertexId ToDiagnosticInteger(std::size_t value)
{
  if (value >
      static_cast<std::size_t>(std::numeric_limits<fem::singular::GlobalVertexId>::max()))
  {
    throw std::overflow_error("Singular mode diagnostic index is out of range!");
  }
  return static_cast<fem::singular::GlobalVertexId>(value);
}

HYPRE_BigInt ToHypreBigInt(fem::singular::GlobalVertexId value)
{
  if (value < std::numeric_limits<HYPRE_BigInt>::min() ||
      value > std::numeric_limits<HYPRE_BigInt>::max())
  {
    throw std::overflow_error("Singular mode diagnostic true DOF is out of range!");
  }
  return static_cast<HYPRE_BigInt>(value);
}

int ToDiagnosticInt(fem::singular::GlobalVertexId value)
{
  if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
  {
    throw std::overflow_error("Singular mode diagnostic value is out of integer range!");
  }
  return static_cast<int>(value);
}

std::pair<PackedModeCoefficientIntegers, PackedModeCoefficientReals>
PackModeCoefficientMeasurement(const SingularModeCoefficientMeasurement &measurement)
{
  PackedModeCoefficientIntegers integers;
  std::size_t next = 0;
  integers[next++] = measurement.mode;
  integers[next++] = static_cast<int>(measurement.space);
  integers[next++] = ToDiagnosticInteger(measurement.true_dof);
  integers[next++] = static_cast<int>(measurement.key.family);
  integers[next++] = measurement.key.order;
  const auto pack_entity = [&integers, &next](const fem::singular::EntityKey &entity)
  {
    if (entity.size == 0 || entity.size > entity.vertices.size())
    {
      throw std::invalid_argument(
          "Singular mode coefficient has an invalid canonical entity!");
    }
    integers[next++] = ToDiagnosticInteger(entity.size);
    for (auto vertex : entity.vertices)
    {
      integers[next++] = vertex;
    }
  };
  pack_entity(measurement.key.singular_entity);
  pack_entity(measurement.key.support_entity);
  pack_entity(measurement.key.component_entity);
  for (int weight : measurement.key.interpolation_weights)
  {
    integers[next++] = weight;
  }
  if (next != integers.size() || measurement.mode < 1 ||
      !std::isfinite(measurement.exponent) ||
      !std::isfinite(measurement.coefficient.real()) ||
      !std::isfinite(measurement.coefficient.imag()))
  {
    throw std::invalid_argument("Singular mode coefficient contains invalid data!");
  }
  return {integers,
          {measurement.exponent, measurement.coefficient.real(),
           measurement.coefficient.imag()}};
}

SingularModeCoefficientMeasurement
UnpackModeCoefficientMeasurement(const PackedModeCoefficientIntegers &integers,
                                 const PackedModeCoefficientReals &reals)
{
  std::size_t next = 0;
  SingularModeCoefficientMeasurement measurement;
  measurement.mode = ToDiagnosticInt(integers[next++]);
  const int space = ToDiagnosticInt(integers[next++]);
  if (space < static_cast<int>(SingularModeSpace::ND) ||
      space > static_cast<int>(SingularModeSpace::H1))
  {
    throw std::invalid_argument("Singular mode coefficient has an invalid space!");
  }
  measurement.space = static_cast<SingularModeSpace>(space);
  measurement.true_dof = ToHypreBigInt(integers[next++]);
  const int family = ToDiagnosticInt(integers[next++]);
  const bool valid_family =
      family == static_cast<int>(fem::singular::HigherOrderBasisFamily::NODE_GRADIENT) ||
      (measurement.space == SingularModeSpace::ND &&
       family == static_cast<int>(fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL));
  if (!valid_family)
  {
    throw std::invalid_argument(
        "Singular mode coefficient has an invalid triangular basis family!");
  }
  measurement.key.family = static_cast<fem::singular::HigherOrderBasisFamily>(family);
  measurement.key.order = ToDiagnosticInt(integers[next++]);
  const auto unpack_entity = [&integers, &next](fem::singular::EntityKey &entity)
  {
    const auto size = integers[next++];
    if (size <= 0 ||
        size > static_cast<fem::singular::GlobalVertexId>(entity.vertices.size()))
    {
      throw std::invalid_argument(
          "Singular mode coefficient has an invalid canonical entity size!");
    }
    entity.size = static_cast<std::size_t>(size);
    for (auto &vertex : entity.vertices)
    {
      vertex = integers[next++];
    }
  };
  unpack_entity(measurement.key.singular_entity);
  unpack_entity(measurement.key.support_entity);
  unpack_entity(measurement.key.component_entity);
  for (auto &weight : measurement.key.interpolation_weights)
  {
    weight = ToDiagnosticInt(integers[next++]);
  }
  measurement.exponent = reals[0];
  measurement.coefficient = {reals[1], reals[2]};
  if (next != integers.size() || measurement.mode < 1 || measurement.true_dof < 0 ||
      measurement.key.order != 1 || !(measurement.exponent > 0.0) ||
      !(measurement.exponent < 1.0) || !std::isfinite(measurement.coefficient.real()) ||
      !std::isfinite(measurement.coefficient.imag()))
  {
    throw std::invalid_argument("Singular mode coefficient unpacked invalid data!");
  }
  return measurement;
}

std::vector<SingularModeCoefficientMeasurement> GatherModeCoefficientMeasurements(
    MPI_Comm comm, int mode,
    const std::vector<fem::singular::H1CoefficientDiagnostic> &nd_real,
    const std::vector<fem::singular::H1CoefficientDiagnostic> &nd_imag,
    HYPRE_BigInt expected_nd_size,
    const std::vector<fem::singular::H1CoefficientDiagnostic> &h1_real,
    const std::vector<fem::singular::H1CoefficientDiagnostic> &h1_imag,
    HYPRE_BigInt expected_h1_size)
{
  auto pair_diagnostics =
      [mode](SingularModeSpace space,
             const std::vector<fem::singular::H1CoefficientDiagnostic> &real,
             const std::vector<fem::singular::H1CoefficientDiagnostic> &imag,
             std::vector<SingularModeCoefficientMeasurement> &paired)
  {
    if (real.size() != imag.size())
    {
      throw std::runtime_error(
          "Real and imaginary singular mode coefficient sets have different sizes!");
    }
    for (std::size_t i = 0; i < real.size(); i++)
    {
      if (real[i].true_dof != imag[i].true_dof || !(real[i].key == imag[i].key) ||
          real[i].exponent != imag[i].exponent)
      {
        throw std::runtime_error(
            "Real and imaginary singular mode coefficient sets have different keys!");
      }
      paired.push_back({mode,
                        space,
                        real[i].true_dof,
                        real[i].key,
                        real[i].exponent,
                        {real[i].coefficient, imag[i].coefficient}});
    }
  };

  std::vector<SingularModeCoefficientMeasurement> local;
  local.reserve(nd_real.size() + h1_real.size());
  pair_diagnostics(SingularModeSpace::ND, nd_real, nd_imag, local);
  pair_diagnostics(SingularModeSpace::H1, h1_real, h1_imag, local);
  bool valid =
      local.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                          std::max(ModeCoefficientIntegerFields, ModeCoefficientRealFields);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Singular mode coefficient diagnostics exceed MPI integer counts!");
  }

  const int local_records = static_cast<int>(local.size());
  const int ranks = Mpi::Size(comm);
  std::vector<int> record_counts(ranks);
  Mpi::Allgather(1, &local_records, record_counts.data(), comm);
  std::vector<int> integer_counts(ranks), integer_displacements(ranks);
  std::vector<int> real_counts(ranks), real_displacements(ranks);
  std::int64_t total_records = 0;
  for (int rank = 0; rank < ranks; rank++)
  {
    if (record_counts[rank] < 0 ||
        total_records > std::numeric_limits<int>::max() / ModeCoefficientIntegerFields -
                            record_counts[rank])
    {
      throw std::overflow_error(
          "Gathered singular mode coefficient diagnostics exceed MPI integer counts!");
    }
    integer_displacements[rank] =
        static_cast<int>(total_records * ModeCoefficientIntegerFields);
    real_displacements[rank] = static_cast<int>(total_records * ModeCoefficientRealFields);
    integer_counts[rank] =
        static_cast<int>(record_counts[rank] * ModeCoefficientIntegerFields);
    real_counts[rank] = static_cast<int>(record_counts[rank] * ModeCoefficientRealFields);
    total_records += record_counts[rank];
  }

  std::vector<fem::singular::GlobalVertexId> local_integers(local.size() *
                                                            ModeCoefficientIntegerFields);
  std::vector<double> local_reals(local.size() * ModeCoefficientRealFields);
  for (std::size_t i = 0; i < local.size(); i++)
  {
    const auto [integers, reals] = PackModeCoefficientMeasurement(local[i]);
    std::copy(integers.begin(), integers.end(),
              local_integers.begin() + i * ModeCoefficientIntegerFields);
    std::copy(reals.begin(), reals.end(),
              local_reals.begin() + i * ModeCoefficientRealFields);
  }

  std::vector<fem::singular::GlobalVertexId> global_integers(
      static_cast<std::size_t>(total_records) * ModeCoefficientIntegerFields);
  std::vector<double> global_reals(static_cast<std::size_t>(total_records) *
                                   ModeCoefficientRealFields);
  Mpi::Allgatherv(static_cast<int>(local_integers.size()), local_integers.data(),
                  global_integers.data(), integer_counts.data(),
                  integer_displacements.data(), comm);
  Mpi::Allgatherv(static_cast<int>(local_reals.size()), local_reals.data(),
                  global_reals.data(), real_counts.data(), real_displacements.data(), comm);

  std::vector<SingularModeCoefficientMeasurement> measurements;
  measurements.reserve(static_cast<std::size_t>(total_records));
  for (std::size_t i = 0; i < static_cast<std::size_t>(total_records); i++)
  {
    PackedModeCoefficientIntegers integers;
    PackedModeCoefficientReals reals;
    std::copy(global_integers.begin() + i * ModeCoefficientIntegerFields,
              global_integers.begin() + (i + 1) * ModeCoefficientIntegerFields,
              integers.begin());
    std::copy(global_reals.begin() + i * ModeCoefficientRealFields,
              global_reals.begin() + (i + 1) * ModeCoefficientRealFields, reals.begin());
    measurements.push_back(UnpackModeCoefficientMeasurement(integers, reals));
  }
  std::sort(measurements.begin(), measurements.end(),
            [](const auto &left, const auto &right)
            { return std::tie(left.space, left.key) < std::tie(right.space, right.key); });
  if (expected_nd_size < 0 || expected_h1_size < 0 ||
      measurements.size() != static_cast<std::size_t>(expected_nd_size + expected_h1_size))
  {
    throw std::runtime_error(
        "Singular mode coefficients do not cover every global enrichment DOF!");
  }
  for (const auto space : {SingularModeSpace::ND, SingularModeSpace::H1})
  {
    const HYPRE_BigInt expected =
        space == SingularModeSpace::ND ? expected_nd_size : expected_h1_size;
    if (expected > std::numeric_limits<int>::max())
    {
      throw std::overflow_error(
          "Singular mode coefficient space exceeds addressable diagnostics!");
    }
    std::vector<bool> seen(static_cast<std::size_t>(expected), false);
    const SingularModeCoefficientMeasurement *previous = nullptr;
    HYPRE_BigInt count = 0;
    for (const auto &measurement : measurements)
    {
      if (measurement.space == space)
      {
        if (measurement.mode != mode || measurement.true_dof < 0 ||
            measurement.true_dof >= expected ||
            seen[static_cast<std::size_t>(measurement.true_dof)] ||
            (previous && !(previous->key < measurement.key)))
        {
          throw std::runtime_error(
              "Singular mode coefficients have duplicate or missing keys or true "
              "DOFs!");
        }
        seen[static_cast<std::size_t>(measurement.true_dof)] = true;
        previous = &measurement;
        count++;
      }
    }
    if (count != expected || std::find(seen.begin(), seen.end(), false) != seen.end())
    {
      throw std::runtime_error(
          "Singular mode coefficients have an inconsistent space size!");
    }
  }

  // Eigenvectors have an arbitrary global phase. Rotate the first canonical
  // coefficient with appreciable magnitude onto the positive real axis so
  // complex coefficient files are reproducible across MPI partitions.
  double maximum_magnitude = 0.0;
  for (const auto &measurement : measurements)
  {
    maximum_magnitude = std::max(maximum_magnitude, std::abs(measurement.coefficient));
  }
  if (maximum_magnitude > 0.0)
  {
    auto anchor = std::find_if(
        measurements.begin(), measurements.end(),
        [maximum_magnitude](const auto &measurement)
        { return std::abs(measurement.coefficient) > 1.0e-12 * maximum_magnitude; });
    MFEM_VERIFY(anchor != measurements.end(),
                "Unable to select a singular mode phase anchor!");
    const std::complex<double> phase =
        std::conj(anchor->coefficient) / std::abs(anchor->coefficient);
    for (auto &measurement : measurements)
    {
      measurement.coefficient *= phase;
    }
    anchor->coefficient = {std::abs(anchor->coefficient), 0.0};
  }
  return measurements;
}

std::pair<PackedModeTipSlopeIntegers, PackedModeTipSlopeReals>
PackModeTipSlopeMeasurement(const SingularModeTipSlopeMeasurement &measurement)
{
  const auto &diagnostic = measurement.diagnostic;
  PackedModeTipSlopeIntegers integers{measurement.mode,
                                      diagnostic.source_element,
                                      ToDiagnosticInteger(diagnostic.feature),
                                      ToDiagnosticInteger(diagnostic.selected_segment),
                                      diagnostic.canonical_vertices[0],
                                      diagnostic.canonical_vertices[1],
                                      diagnostic.canonical_vertices[2],
                                      diagnostic.sample_count,
                                      diagnostic.valid ? 1 : 0};
  PackedModeTipSlopeReals reals{diagnostic.exponent,
                                diagnostic.expected_slope,
                                diagnostic.fitted_slope,
                                diagnostic.r_squared,
                                diagnostic.minimum_distance,
                                diagnostic.maximum_distance,
                                diagnostic.field_norm_at_minimum_distance,
                                diagnostic.field_norm_at_maximum_distance};
  if (measurement.mode < 1 || diagnostic.source_element < 0 ||
      diagnostic.sample_count < 3 ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }))
  {
    throw std::invalid_argument(
        "Singular BoundaryMode tip slope diagnostic contains invalid data!");
  }
  return {integers, reals};
}

SingularModeTipSlopeMeasurement
UnpackModeTipSlopeMeasurement(const PackedModeTipSlopeIntegers &integers,
                              const PackedModeTipSlopeReals &reals)
{
  SingularModeTipSlopeMeasurement measurement;
  measurement.mode = ToDiagnosticInt(integers[0]);
  auto &diagnostic = measurement.diagnostic;
  diagnostic.source_element = integers[1];
  if (integers[2] < 0 || integers[3] < 0)
  {
    throw std::invalid_argument(
        "Singular BoundaryMode tip slope diagnostic has a negative feature index!");
  }
  diagnostic.feature = static_cast<std::size_t>(integers[2]);
  diagnostic.selected_segment = static_cast<std::size_t>(integers[3]);
  std::copy(integers.begin() + 4, integers.begin() + 7,
            diagnostic.canonical_vertices.begin());
  diagnostic.sample_count = ToDiagnosticInt(integers[7]);
  diagnostic.valid = integers[8] != 0;
  diagnostic.exponent = reals[0];
  diagnostic.expected_slope = reals[1];
  diagnostic.fitted_slope = reals[2];
  diagnostic.r_squared = reals[3];
  diagnostic.minimum_distance = reals[4];
  diagnostic.maximum_distance = reals[5];
  diagnostic.field_norm_at_minimum_distance = reals[6];
  diagnostic.field_norm_at_maximum_distance = reals[7];
  if (measurement.mode < 1 || diagnostic.source_element < 0 ||
      diagnostic.sample_count < 3 || (integers[8] != 0 && integers[8] != 1) ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }) ||
      !(diagnostic.exponent > 0.0) || !(diagnostic.exponent < 1.0) ||
      diagnostic.expected_slope != diagnostic.exponent - 1.0 ||
      (diagnostic.valid && (!(diagnostic.minimum_distance > 0.0) ||
                            !(diagnostic.maximum_distance > diagnostic.minimum_distance) ||
                            !(diagnostic.field_norm_at_minimum_distance > 0.0) ||
                            !(diagnostic.field_norm_at_maximum_distance > 0.0))))
  {
    throw std::invalid_argument(
        "Singular BoundaryMode tip slope diagnostic unpacked invalid data!");
  }
  return measurement;
}

std::vector<SingularModeTipSlopeMeasurement> GatherModeTipSlopeMeasurements(
    MPI_Comm comm, int mode,
    const std::vector<fem::singular::H1TipSlopeDiagnostic> &local_diagnostics)
{
  bool valid = local_diagnostics.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   std::max(ModeTipSlopeIntegerFields, ModeTipSlopeRealFields);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Singular BoundaryMode tip slope diagnostics exceed MPI integer counts!");
  }

  const int local_records = static_cast<int>(local_diagnostics.size());
  const int ranks = Mpi::Size(comm);
  std::vector<int> record_counts(ranks);
  Mpi::Allgather(1, &local_records, record_counts.data(), comm);
  std::vector<int> integer_counts(ranks), integer_displacements(ranks);
  std::vector<int> real_counts(ranks), real_displacements(ranks);
  std::int64_t total_records = 0;
  for (int rank = 0; rank < ranks; rank++)
  {
    if (record_counts[rank] < 0 ||
        total_records > std::numeric_limits<int>::max() / ModeTipSlopeIntegerFields -
                            record_counts[rank])
    {
      throw std::overflow_error(
          "Gathered singular BoundaryMode tip slope diagnostics exceed MPI integer "
          "counts!");
    }
    integer_displacements[rank] =
        static_cast<int>(total_records * ModeTipSlopeIntegerFields);
    real_displacements[rank] = static_cast<int>(total_records * ModeTipSlopeRealFields);
    integer_counts[rank] =
        static_cast<int>(record_counts[rank] * ModeTipSlopeIntegerFields);
    real_counts[rank] = static_cast<int>(record_counts[rank] * ModeTipSlopeRealFields);
    total_records += record_counts[rank];
  }

  std::vector<fem::singular::GlobalVertexId> local_integers(local_diagnostics.size() *
                                                            ModeTipSlopeIntegerFields);
  std::vector<double> local_reals(local_diagnostics.size() * ModeTipSlopeRealFields);
  for (std::size_t i = 0; i < local_diagnostics.size(); i++)
  {
    const auto [integers, reals] =
        PackModeTipSlopeMeasurement({mode, local_diagnostics[i]});
    std::copy(integers.begin(), integers.end(),
              local_integers.begin() + i * ModeTipSlopeIntegerFields);
    std::copy(reals.begin(), reals.end(), local_reals.begin() + i * ModeTipSlopeRealFields);
  }

  std::vector<fem::singular::GlobalVertexId> global_integers(
      static_cast<std::size_t>(total_records) * ModeTipSlopeIntegerFields);
  std::vector<double> global_reals(static_cast<std::size_t>(total_records) *
                                   ModeTipSlopeRealFields);
  Mpi::Allgatherv(static_cast<int>(local_integers.size()), local_integers.data(),
                  global_integers.data(), integer_counts.data(),
                  integer_displacements.data(), comm);
  Mpi::Allgatherv(static_cast<int>(local_reals.size()), local_reals.data(),
                  global_reals.data(), real_counts.data(), real_displacements.data(), comm);

  std::vector<SingularModeTipSlopeMeasurement> measurements;
  measurements.reserve(static_cast<std::size_t>(total_records));
  for (std::size_t i = 0; i < static_cast<std::size_t>(total_records); i++)
  {
    PackedModeTipSlopeIntegers integers;
    PackedModeTipSlopeReals reals;
    std::copy(global_integers.begin() + i * ModeTipSlopeIntegerFields,
              global_integers.begin() + (i + 1) * ModeTipSlopeIntegerFields,
              integers.begin());
    std::copy(global_reals.begin() + i * ModeTipSlopeRealFields,
              global_reals.begin() + (i + 1) * ModeTipSlopeRealFields, reals.begin());
    measurements.push_back(UnpackModeTipSlopeMeasurement(integers, reals));
  }
  std::sort(measurements.begin(), measurements.end(),
            [](const auto &left, const auto &right)
            {
              return std::tie(left.diagnostic.source_element, left.diagnostic.feature) <
                     std::tie(right.diagnostic.source_element, right.diagnostic.feature);
            });
  if (measurements.empty())
  {
    throw std::runtime_error(
        "Singular BoundaryMode tip slope diagnostics found no global triangle sectors!");
  }
  for (std::size_t i = 0; i < measurements.size(); i++)
  {
    if (measurements[i].mode != mode ||
        (i > 0 && std::tie(measurements[i - 1].diagnostic.source_element,
                           measurements[i - 1].diagnostic.feature) ==
                      std::tie(measurements[i].diagnostic.source_element,
                               measurements[i].diagnostic.feature)))
    {
      throw std::runtime_error(
          "Singular BoundaryMode tip slope diagnostics have duplicate triangle "
          "sectors!");
    }
  }
  return measurements;
}

void InsertIntegerDiagnosticColumn(Table &table, std::string name, std::string header)
{
  table.insert(Column(std::move(name), std::move(header), 0, 0, 0, ""));
  table[table.n_cols() - 1].print_as_int = true;
}

void WriteSingularModeMeasurements(const fs::path &post_dir, const IoData &iodata,
                                   const std::vector<SingularModeMeasurement> &measurements,
                                   bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "mode-kn.csv");
  output.table.insert(Column("idx", "m", -1, 0, 2, ""));
  output.table.insert("kn_re", "Re{kn} (1/m)");
  output.table.insert("kn_im", "Im{kn} (1/m)");
  output.table.insert("neff_re", "Re{n_eff}");
  output.table.insert("neff_im", "Im{n_eff}");
  output.table.insert("err_back", "Error (Bkwd.)");
  output.table.insert("err_abs", "Error (Abs.)");
  const double inverse_length =
      1.0 / iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
  for (const auto &measurement : measurements)
  {
    output.table["idx"] << measurement.mode;
    output.table["kn_re"] << inverse_length * measurement.propagation_constant.real();
    output.table["kn_im"] << inverse_length * measurement.propagation_constant.imag();
    output.table["neff_re"] << measurement.effective_index.real();
    output.table["neff_im"] << measurement.effective_index.imag();
    output.table["err_back"] << measurement.backward_error;
    output.table["err_abs"] << measurement.absolute_error;
  }
  output.WriteFullTableTrunc();

  TableWithCSVFile diagnostics(post_dir / "singular-mode-diagnostics.csv");
  diagnostics.table.insert(Column("idx", "m", -1, 0, 2, ""));
  diagnostics.table.insert("power_re", "Re{P} (normalized)");
  diagnostics.table.insert("power_im", "Im{P} (normalized)");
  diagnostics.table.insert("nd_gradient_norm", "ND gradient coefficient norm");
  diagnostics.table.insert("nd_rotational_norm", "ND rotational coefficient norm");
  diagnostics.table.insert("h1_gradient_norm", "H1 gradient coefficient norm");
  diagnostics.table.insert("electric_transverse_energy",
                           "Transverse electric field energy (J)");
  diagnostics.table.insert("electric_normal_energy", "Normal electric field energy (J)");
  diagnostics.table.insert("electric_total_energy", "Total electric field energy (J)");
  diagnostics.table.insert("magnetic_transverse_energy",
                           "Transverse magnetic field energy (J)");
  diagnostics.table.insert("magnetic_normal_energy", "Normal magnetic field energy (J)");
  diagnostics.table.insert("magnetic_total_energy", "Total magnetic field energy (J)");
  const double energy_scale = iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0);
  for (const auto &measurement : measurements)
  {
    diagnostics.table["idx"] << measurement.mode;
    diagnostics.table["power_re"] << measurement.normalized_power.real();
    diagnostics.table["power_im"] << measurement.normalized_power.imag();
    diagnostics.table["nd_gradient_norm"] << measurement.coefficient_norms.nd_gradient;
    diagnostics.table["nd_rotational_norm"] << measurement.coefficient_norms.nd_rotational;
    diagnostics.table["h1_gradient_norm"] << measurement.coefficient_norms.h1_gradient;
    diagnostics.table["electric_transverse_energy"]
        << energy_scale * measurement.field_energies.electric_transverse;
    diagnostics.table["electric_normal_energy"]
        << energy_scale * measurement.field_energies.electric_normal;
    diagnostics.table["electric_total_energy"]
        << energy_scale * (measurement.field_energies.electric_transverse +
                           measurement.field_energies.electric_normal);
    diagnostics.table["magnetic_transverse_energy"]
        << energy_scale * measurement.field_energies.magnetic_transverse;
    diagnostics.table["magnetic_normal_energy"]
        << energy_scale * measurement.field_energies.magnetic_normal;
    diagnostics.table["magnetic_total_energy"]
        << energy_scale * (measurement.field_energies.magnetic_transverse +
                           measurement.field_energies.magnetic_normal);
  }
  diagnostics.WriteFullTableTrunc();
}

void WriteSingularModeCoefficientMeasurements(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularModeCoefficientMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  for (const auto space : {SingularModeSpace::ND, SingularModeSpace::H1})
  {
    const bool nd = space == SingularModeSpace::ND;
    TableWithCSVFile output(post_dir / (nd ? "singular-mode-nd-coefficients.csv"
                                           : "singular-mode-h1-coefficients.csv"));
    InsertIntegerDiagnosticColumn(output.table, "mode", "m");
    InsertIntegerDiagnosticColumn(output.table, "true_dof", "enrichment_true_dof");
    InsertIntegerDiagnosticColumn(output.table, "family",
                                  nd ? "family (0 = node gradient; 1 = node rotational)"
                                     : "family (0 = node gradient)");
    InsertIntegerDiagnosticColumn(output.table, "order", "singular_order");
    output.table.insert("nu", "nu");
    for (const std::string_view prefix : {"singular", "support", "component"})
    {
      InsertIntegerDiagnosticColumn(output.table, fmt::format("{}_size", prefix),
                                    fmt::format("{}_entity_size", prefix));
      for (int vertex = 0; vertex < 4; vertex++)
      {
        InsertIntegerDiagnosticColumn(output.table, fmt::format("{}_v{}", prefix, vertex),
                                      fmt::format("{}_entity_vertex_{}", prefix, vertex));
      }
    }
    for (int weight = 0; weight < 4; weight++)
    {
      InsertIntegerDiagnosticColumn(output.table, fmt::format("weight_{}", weight),
                                    fmt::format("interpolation_weight_{}", weight));
    }
    output.table.insert("coefficient_re",
                        nd ? "Re{basis coefficient} (V)" : "Re{basis coefficient} (V/m)");
    output.table.insert("coefficient_im",
                        nd ? "Im{basis coefficient} (V)" : "Im{basis coefficient} (V/m)");

    const double coefficient_scale =
        nd ? iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(1.0)
           : iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
    for (const auto &measurement : measurements)
    {
      if (measurement.space != space)
      {
        continue;
      }
      output.table["mode"] << measurement.mode;
      output.table["true_dof"] << static_cast<double>(measurement.true_dof);
      output.table["family"] << static_cast<int>(measurement.key.family);
      output.table["order"] << measurement.key.order;
      output.table["nu"] << measurement.exponent;
      const auto write_entity =
          [&output](std::string_view prefix, const fem::singular::EntityKey &entity)
      {
        output.table[fmt::format("{}_size", prefix)] << static_cast<double>(entity.size);
        for (int vertex = 0; vertex < 4; vertex++)
        {
          output.table[fmt::format("{}_v{}", prefix, vertex)]
              << static_cast<double>(entity.vertices[vertex]);
        }
      };
      write_entity("singular", measurement.key.singular_entity);
      write_entity("support", measurement.key.support_entity);
      write_entity("component", measurement.key.component_entity);
      for (int weight = 0; weight < 4; weight++)
      {
        output.table[fmt::format("weight_{}", weight)]
            << measurement.key.interpolation_weights[weight];
      }
      output.table["coefficient_re"] << coefficient_scale * measurement.coefficient.real();
      output.table["coefficient_im"] << coefficient_scale * measurement.coefficient.imag();
    }
    output.WriteFullTableTrunc();
  }
}

void WriteSingularModeTipSlopeMeasurements(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularModeTipSlopeMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "singular-mode-tip-slopes.csv");
  for (const auto &[name, header] :
       std::array<std::pair<std::string_view, std::string_view>, 9>{
           std::pair{"mode", "m"}, std::pair{"element", "source_element"},
           std::pair{"feature", "tip_feature"}, std::pair{"segment", "selected_segment"},
           std::pair{"tip_v", "tip_vertex"}, std::pair{"opposite_v0", "opposite_vertex_0"},
           std::pair{"opposite_v1", "opposite_vertex_1"},
           std::pair{"samples", "sample_count"}, std::pair{"valid", "fit_valid"}})
  {
    InsertIntegerDiagnosticColumn(output.table, std::string(name), std::string(header));
  }
  output.table.insert("nu", "nu");
  output.table.insert("expected", "expected_slope");
  output.table.insert("fitted", "fitted_slope");
  output.table.insert("r_squared", "R_squared");
  output.table.insert("rho_min", "minimum_distance (m)");
  output.table.insert("rho_max", "maximum_distance (m)");
  output.table.insert("field_min", "|Et| at minimum_distance (V/m)");
  output.table.insert("field_max", "|Et| at maximum_distance (V/m)");

  const double length_scale = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
  const double field_scale = iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  for (const auto &measurement : measurements)
  {
    const auto &diagnostic = measurement.diagnostic;
    output.table["mode"] << measurement.mode;
    output.table["element"] << static_cast<double>(diagnostic.source_element);
    output.table["feature"] << static_cast<double>(diagnostic.feature);
    output.table["segment"] << static_cast<double>(diagnostic.selected_segment);
    output.table["tip_v"] << static_cast<double>(diagnostic.canonical_vertices[0]);
    output.table["opposite_v0"] << static_cast<double>(diagnostic.canonical_vertices[1]);
    output.table["opposite_v1"] << static_cast<double>(diagnostic.canonical_vertices[2]);
    output.table["samples"] << diagnostic.sample_count;
    output.table["valid"] << (diagnostic.valid ? 1.0 : 0.0);
    output.table["nu"] << diagnostic.exponent;
    output.table["expected"] << diagnostic.expected_slope;
    output.table["fitted"] << diagnostic.fitted_slope;
    output.table["r_squared"] << diagnostic.r_squared;
    output.table["rho_min"] << length_scale * diagnostic.minimum_distance;
    output.table["rho_max"] << length_scale * diagnostic.maximum_distance;
    output.table["field_min"] << field_scale * diagnostic.field_norm_at_minimum_distance;
    output.table["field_max"] << field_scale * diagnostic.field_norm_at_maximum_distance;
  }
  output.WriteFullTableTrunc();
}

void WriteSingularModeSurfaceMeasurements(
    const fs::path &post_dir,
    const std::vector<SingularModeSurfaceMeasurement> &measurements, bool root)
{
  if (!root || measurements.empty() || measurements.front().interfaces.empty())
  {
    return;
  }
  TableWithCSVFile output(post_dir / "surface-Q.csv");
  output.table.insert(Column("idx", "m", -1, 0, 2, ""));
  for (const auto &interface : measurements.front().interfaces)
  {
    output.table.insert(fmt::format("p_{}", interface.index),
                        fmt::format("p_surf[{}]", interface.index));
    output.table.insert(fmt::format("Q_{}", interface.index),
                        fmt::format("Q_surf[{}]", interface.index));
  }
  for (const auto &measurement : measurements)
  {
    if (measurement.interfaces.size() != measurements.front().interfaces.size())
    {
      throw std::runtime_error(
          "Singular BoundaryMode surface measurements changed interface count!");
    }
    output.table["idx"] << measurement.mode;
    for (std::size_t i = 0; i < measurement.interfaces.size(); i++)
    {
      const auto &interface = measurement.interfaces[i];
      if (interface.index != measurements.front().interfaces[i].index)
      {
        throw std::runtime_error(
            "Singular BoundaryMode surface measurements changed interface ordering!");
      }
      output.table[fmt::format("p_{}", interface.index)] << interface.participation;
      output.table[fmt::format("Q_{}", interface.index)] << interface.quality_factor;
    }
  }
  output.WriteFullTableTrunc();
}

class SingularBoundaryModeParaviewOutput
{
private:
  mfem::ParMesh &mesh;
  mfem::L2_FECollection sample_collection;
  mfem::ParFiniteElementSpace scalar_space;
  mfem::ParFiniteElementSpace vector_space;
  mfem::ParGridFunction et_real, et_imag;
  mfem::ParGridFunction en_real, en_imag;
  mfem::ParGridFunction en_negative_gradient_real, en_negative_gradient_imag;
  mfem::ParGridFunction bz_real, bz_imag;
  mfem::ParGridFunction bt_real, bt_imag;
  mfem::ParaViewDataCollection data_collection;
  const Units &units;

public:
  SingularBoundaryModeParaviewOutput(const fs::path &post_dir, mfem::ParMesh &mesh,
                                     int order, const Units &units)
    : mesh(mesh), sample_collection(order, mesh.Dimension(), mfem::BasisType::GaussLegendre,
                                    mfem::FiniteElement::VALUE),
      scalar_space(&mesh, &sample_collection),
      vector_space(&mesh, &sample_collection, 2, mfem::Ordering::byVDIM),
      et_real(&vector_space), et_imag(&vector_space), en_real(&scalar_space),
      en_imag(&scalar_space), en_negative_gradient_real(&vector_space),
      en_negative_gradient_imag(&vector_space), bz_real(&scalar_space),
      bz_imag(&scalar_space), bt_real(&vector_space), bt_imag(&vector_space),
      data_collection(post_dir / "paraview" / "boundarymode", &mesh), units(units)
  {
    MFEM_VERIFY(order >= 1 && mesh.Dimension() == 2 && mesh.SpaceDimension() == 2,
                "Singular BoundaryMode field visualization requires positive-order "
                "sampling on a two-dimensional mesh!");
    RemovePreviousOutput(post_dir / "paraview", mesh.GetComm());
    data_collection.SetCycle(-1);
    data_collection.SetDataFormat(mfem::VTKFormat::BINARY32);
#if defined(MFEM_USE_ZLIB)
    data_collection.SetCompressionLevel(-1);
#else
    data_collection.SetCompressionLevel(0);
#endif
    data_collection.SetHighOrderOutput(true);
    data_collection.SetLevelsOfDetail(order);
    data_collection.RegisterField("Et_real", &et_real);
    data_collection.RegisterField("Et_imag", &et_imag);
    data_collection.RegisterField("En_real", &en_real);
    data_collection.RegisterField("En_imag", &en_imag);
    data_collection.RegisterField("Bz_real", &bz_real);
    data_collection.RegisterField("Bz_imag", &bz_imag);
    data_collection.RegisterField("Bt_real", &bt_real);
    data_collection.RegisterField("Bt_imag", &bt_imag);
  }

  void Write(fem::singular::TriangleEnrichedNDFieldEvaluator &et_real_evaluator,
             fem::singular::TriangleEnrichedNDFieldEvaluator &et_imag_evaluator,
             fem::singular::TriangleEnrichedH1FieldEvaluator &en_real_evaluator,
             fem::singular::TriangleEnrichedH1FieldEvaluator &en_imag_evaluator,
             std::complex<double> kn, double omega, int step)
  {
    MFEM_VERIFY(std::isfinite(omega) && omega > 0.0 && std::isfinite(kn.real()) &&
                    std::isfinite(kn.imag()),
                "Singular BoundaryMode field output received invalid mode data!");

    et_real_evaluator.ProjectToDiscontinuousGridFunctions(et_real, bz_real);
    et_imag_evaluator.ProjectToDiscontinuousGridFunctions(et_imag, bz_imag);
    en_real_evaluator.ProjectToDiscontinuousGridFunctions(en_real,
                                                          en_negative_gradient_real);
    en_imag_evaluator.ProjectToDiscontinuousGridFunctions(en_imag,
                                                          en_negative_gradient_imag);

    // H1 sampling returns -grad(En). Reconstruct B at the same interior
    // interpolation nodes using the tested BoundaryMode convention helper.
    const int scalar_dofs = vector_space.GetNDofs();
    MFEM_VERIFY(vector_space.GetVSize() == 2 * scalar_dofs,
                "Unexpected discontinuous vector-space ordering!");
    for (int i = 0; i < scalar_dofs; i++)
    {
      const fem::singular::Vector2 etr{et_real[i], et_real[scalar_dofs + i]};
      const fem::singular::Vector2 eti{et_imag[i], et_imag[scalar_dofs + i]};
      const fem::singular::Vector2 grad_enr{-en_negative_gradient_real[i],
                                            -en_negative_gradient_real[scalar_dofs + i]};
      const fem::singular::Vector2 grad_eni{-en_negative_gradient_imag[i],
                                            -en_negative_gradient_imag[scalar_dofs + i]};
      const auto magnetic = fem::singular::ReconstructTriangleBoundaryModeMagneticField(
          etr, eti, bz_real[i], bz_imag[i], grad_enr, grad_eni, kn, omega);
      bz_real[i] = magnetic.normal_real;
      bz_imag[i] = magnetic.normal_imag;
      bt_real[i] = magnetic.transverse_real[0];
      bt_real[scalar_dofs + i] = magnetic.transverse_real[1];
      bt_imag[i] = magnetic.transverse_imag[0];
      bt_imag[scalar_dofs + i] = magnetic.transverse_imag[1];
    }

    const double mesh_scale = units.GetMeshLengthRelativeScale();
    const double electric_scale = units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
    const double magnetic_scale = units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);
    mesh::DimensionalizeMesh(mesh, mesh_scale);
    for (auto *field : {&et_real, &et_imag, &en_real, &en_imag})
    {
      *field *= electric_scale;
    }
    for (auto *field : {&bz_real, &bz_imag, &bt_real, &bt_imag})
    {
      *field *= magnetic_scale;
    }
    data_collection.SetCycle(step);
    data_collection.SetTime(step + 1);
    data_collection.Save();
    for (auto *field : {&et_real, &et_imag, &en_real, &en_imag})
    {
      *field *= 1.0 / electric_scale;
    }
    for (auto *field : {&bz_real, &bz_imag, &bt_real, &bt_imag})
    {
      *field *= 1.0 / magnetic_scale;
    }
    mesh::NondimensionalizeMesh(mesh, mesh_scale);
    Mpi::Barrier(mesh.GetComm());
  }
};

}  // namespace

BoundaryModeSolver::BoundaryModeSolver(const IoData &iodata, bool root, int size,
                                       int num_thread, const char *git_tag)
  : BaseSolver(iodata, root, size, num_thread, git_tag)
{
}

void BoundaryModeSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                                    MPI_Comm comm) const
{
  const auto &bm = iodata.solver.boundary_mode;
  if (bm.attributes.empty())
  {
    // Direct 2D path: caller-supplied mesh is already the solve mesh.
    BaseSolver::Preprocess(iodata, smesh, comm);
  }
  else
  {
    SubmeshFrame frame;

    // Ranks holding the serial mesh run the extraction locally via mfem::SubMesh;
    // non-holders receive the frame via broadcast below.
    if (smesh)
    {
      MFEM_VERIFY(smesh->Dimension() == 3,
                  "BoundaryMode with \"Attributes\" requires a 3D input mesh!");
      Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");

      mfem::Array<int> attr_list;
      attr_list.Append(bm.attributes.data(), bm.attributes.size());

      // Parent boundary face types whose intersection with the cross-section should
      // become 2D boundary edges: BCs that pin the tangential E-field (PEC, auxpec,
      // impedance, conductivity, rational impedance, farfield) plus other-waveports
      // (relabeled to PEC below).
      std::vector<int> internal_bdr_attrs;
      const auto &bdr = iodata.boundaries;
      auto append = [&](const auto &src)
      { internal_bdr_attrs.insert(internal_bdr_attrs.end(), src.begin(), src.end()); };
      append(bdr.pec.attributes);
      append(bdr.auxpec.attributes);
      for (const auto &d : bdr.impedance)
      {
        append(d.attributes);
      }
      for (const auto &d : bdr.conductivity)
      {
        append(d.attributes);
      }
      for (const auto &d : bdr.rational_impedance)
      {
        append(d.attributes);
      }
      append(bdr.farfield.attributes);

      std::unordered_set<int> other_waveport_attrs;
      for (const auto &[idx, data] : bdr.waveport)
      {
        for (auto a : data.attributes)
        {
          if (std::find(bm.attributes.begin(), bm.attributes.end(), a) !=
              bm.attributes.end())
          {
            continue;
          }
          other_waveport_attrs.insert(a);
          internal_bdr_attrs.push_back(a);
        }
      }

      auto extracted =
          ExtractBoundary2DSubmesh(*smesh, attr_list, internal_bdr_attrs, frame);

      // Relabel other-waveport edges as PEC on this cross-section: for the
      // BoundaryMode solve they act as conducting boundaries.
      if (!other_waveport_attrs.empty())
      {
        MFEM_VERIFY(!bdr.pec.attributes.empty(),
                    "BoundaryMode submesh extraction found other-waveport edges on the "
                    "cross-section. Define at least one PEC boundary attribute to "
                    "relabel them to.");
        const int pec_attr =
            *std::min_element(bdr.pec.attributes.begin(), bdr.pec.attributes.end());
        int relabelled = 0;
        for (int sbe = 0; sbe < extracted->GetNBE(); sbe++)
        {
          if (other_waveport_attrs.count(extracted->GetBdrAttribute(sbe)))
          {
            extracted->SetBdrAttribute(sbe, pec_attr);
            relabelled++;
          }
        }
        if (relabelled > 0)
        {
          Mpi::Print(
              " Relabelled {:d} other-waveport boundary edge(s) as PEC (attr {:d})\n",
              relabelled, pec_attr);
        }
      }

      smesh = std::move(extracted);
      Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", frame.normal(0),
                 frame.normal(1), frame.normal(2));
    }

    // Broadcast user-units frame to non-holding ranks.
    Mpi::Broadcast(3, frame.centroid.HostReadWrite(), 0, comm);
    Mpi::Broadcast(3, frame.e1.HostReadWrite(), 0, comm);
    Mpi::Broadcast(3, frame.e2.HostReadWrite(), 0, comm);
    Mpi::Broadcast(3, frame.normal.HostReadWrite(), 0, comm);

    // Bake the tangent frame into iodata (pre-nondim, user units throughout):
    //   - material tensors rotated into the local frame so MaterialOperator constructs
    //     natively-2D downstream;
    //   - impedance / voltage 3D path coordinates projected to 2D local coords so
    //     PostOperator parses them as 2D paths directly.
    RotateMaterialDefinitions(iodata.domains.materials, frame.e1, frame.e2, frame.normal);
    for (auto &[idx, imp] : iodata.boundaries.postpro.impedance)
    {
      ProjectPathTo2D(imp.voltage_path, frame.centroid, frame.e1, frame.e2);
      ProjectPathTo2D(imp.current_path, frame.centroid, frame.e1, frame.e2);
    }
    for (auto &[idx, vol] : iodata.boundaries.postpro.voltage)
    {
      ProjectPathTo2D(vol.voltage_path, frame.centroid, frame.e1, frame.e2);
    }

    // Resolve Lc from the 2D solve mesh and nondimensionalize in one pass.
    BaseSolver::Preprocess(iodata, smesh, comm);
  }

  if (!iodata.solver.singular_elements.Enabled())
  {
    return;
  }

  serial_singular_features = {};
  local_singular_features = {};
  source_vertex_ids.clear();
  source_element_ids.clear();
  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(smesh && smesh->Dimension() == 2 && smesh->SpaceDimension() == 2,
                "BoundaryMode singular enrichment requires a two-dimensional "
                "triangular solve mesh!");
    MFEM_VERIFY(iodata.solver.singular_elements.order == 1,
                "Triangular BoundaryMode singular enrichment currently supports only "
                "Solver.SingularElements.Order = 1!");
    serial_singular_features = fem::singular::ExtractSerialLineFeatures(
        *smesh, iodata.solver.singular_elements.attributes,
        GetSingularTriangleMaterials(iodata));
    ValidateSingularLossTangents(iodata, serial_singular_features);
  }
  fem::singular::BroadcastSerialLineTipFeatures(serial_singular_features, comm);
  MFEM_VERIFY(!serial_singular_features.Empty(),
              "BoundaryMode singular feature extraction produced no PEC line tips or "
              "corners!");
}

bool BoundaryModeSolver::RequiresSourceSerialMeshMetadata() const
{
  return iodata.solver.singular_elements.Enabled();
}

void BoundaryModeSolver::ProcessPartitionedMesh(
    const mfem::ParMesh &parallel_mesh, const mesh::PartitionMetadata &metadata) const
{
  MFEM_VERIFY(iodata.solver.singular_elements.Enabled(),
              "Unexpected singular partition metadata for an unenriched BoundaryMode "
              "solve!");
  MFEM_VERIFY(!serial_singular_features.Empty(),
              "BoundaryMode singular source feature blueprint was not initialized!");
  local_singular_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_singular_features, parallel_mesh, metadata.source_vertex_ids,
      metadata.source_element_ids);
  source_vertex_ids = metadata.source_vertex_ids;
  source_element_ids = metadata.source_element_ids;
}

std::pair<ErrorIndicator, long long int>
BoundaryModeSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &bm = iodata.solver.boundary_mode;
  const int num_modes = bm.n;
  const double freq_GHz = bm.freq;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  BlockTimer bt0(Timer::CONSTRUCT);
  MaterialOperator mat_op(iodata, *mesh.back());

  const bool singular = iodata.solver.singular_elements.Enabled();
  MFEM_VERIFY(!singular || (!local_singular_features.Empty() &&
                            source_vertex_ids.size() ==
                                static_cast<std::size_t>(mesh.back()->Get().GetNV()) &&
                            source_element_ids.size() ==
                                static_cast<std::size_t>(mesh.back()->Get().GetNE())),
              "BoundaryMode singular feature topology was not reconstructed on the "
              "solve mesh!");
  BoundaryModeOperator mode_op(iodata, mesh, mat_op,
                               singular ? &local_singular_features : nullptr,
                               singular ? &source_vertex_ids : nullptr);

  const int nd_size = mode_op.GetNDTrueVSize();
  const auto &dbc_tdof_list = mode_op.GetCombinedDbcTDofList();

  if (const auto n_levels = mode_op.GetNDSpaceHierarchy().GetNumLevels(); n_levels > 1)
  {
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(n_levels));
  }
  const auto which_eig = (bm.target > 0.0) ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                                           : EigenvalueSolver::WhichType::LARGEST_REAL;
  ModeEigenSolver eig(mode_op, dbc_tdof_list, num_modes, bm.max_size, bm.tol, which_eig,
                      iodata.solver.linear, bm.type, iodata.problem.verbose);

  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator;
  std::unique_ptr<PostOperator<ProblemType::BOUNDARYMODE>> post_op;
  if (!singular)
  {
    estimator = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        mode_op.GetMaterialOp(), mode_op.GetNDSpaceHierarchy(),
        mode_op.GetRTSpaceHierarchy(), mode_op.GetCurlSpace(),
        mode_op.GetH1SpaceHierarchy(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
    post_op = std::make_unique<PostOperator<ProblemType::BOUNDARYMODE>>(iodata, mode_op);
  }

  ErrorIndicator indicator;

  // Determine kn_target.
  double kn_target;
  if (bm.target > 0.0)
  {
    kn_target = bm.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", bm.target, kn_target);
  }
  else
  {
    const double mu_eps_max = mode_op.GetMaterialOp().GetMaxMuEpsilon();
    kn_target = omega * std::sqrt(1.1 * mu_eps_max);
    Mpi::Print(" Auto kn_target = {:.6e} (from max(mu_r) * max(epsilon_r) = {:.6e})\n",
               kn_target, mu_eps_max);
  }

  // Solve the eigenvalue problem.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);
  const double sigma = -kn_target * kn_target;
  auto result = eig.Solve(omega, sigma);
  int num_conv = result.num_converged;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", result.sigma);

  for (int i = 0; i < num_conv; i++)
  {
    auto kn = eig.GetPropagationConstant(i);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i, kn.real(),
               kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = eig.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  const int h1_size = mode_op.GetH1TrueVSize();
  std::vector<SingularModeMeasurement> singular_measurements;
  singular_measurements.reserve(std::min(num_conv, num_modes));
  std::vector<SingularModeCoefficientMeasurement> singular_coefficient_measurements;
  std::vector<SingularModeTipSlopeMeasurement> singular_tip_slope_measurements;
  std::vector<SingularModeSurfaceMeasurement> singular_surface_measurements;
  std::unique_ptr<fem::singular::TriangleEnrichedNDFieldEvaluator>
      singular_et_real_evaluator, singular_et_imag_evaluator;
  std::unique_ptr<fem::singular::TriangleEnrichedH1FieldEvaluator>
      singular_en_real_evaluator, singular_en_imag_evaluator;
  std::unique_ptr<SingularBoundaryModeParaviewOutput> singular_paraview;
  std::unique_ptr<TriangleSingularSurfacePostOperator> singular_surface_postoperator;
  if (singular)
  {
    const auto &singular_dofs = mode_op.GetSingularDofTopology();
    const auto &singular_numbering = mode_op.GetSingularDofNumbering();
    singular_et_real_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            singular_dofs, singular_numbering, mode_op.GetNDSpace().Get());
    singular_et_imag_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            singular_dofs, singular_numbering, mode_op.GetNDSpace().Get());
    singular_en_real_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedH1FieldEvaluator>(
            singular_dofs, singular_numbering, mode_op.GetH1Space().Get());
    singular_en_imag_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedH1FieldEvaluator>(
            singular_dofs, singular_numbering, mode_op.GetH1Space().Get());
    singular_surface_postoperator = std::make_unique<TriangleSingularSurfacePostOperator>(
        iodata.boundaries.postpro, mode_op.GetMaterialOp(), mode_op.GetNDSpace().Get());
    if (iodata.problem.output_formats.paraview && bm.n_post > 0)
    {
      singular_paraview = std::make_unique<SingularBoundaryModeParaviewOutput>(
          post_dir, mode_op.GetMesh().Get(), std::max(2, mode_op.GetSolverOrder() + 1),
          iodata.units);
    }
    SaveMetadata(
        "SingularElements",
        nlohmann::json{
            {"Dimension", 2},
            {"Simulation", "BoundaryMode"},
            {"StandardOrder", iodata.solver.order},
            {"SingularOrder", iodata.solver.singular_elements.order},
            {"NDErichmentDegreesOfFreedom",
             mode_op.GetNDGlobalTrueVSize() - mode_op.GetNDSpace().GlobalTrueVSize()},
            {"H1EnrichmentDegreesOfFreedom",
             mode_op.GetH1GlobalTrueVSize() - mode_op.GetH1Space().GlobalTrueVSize()},
            {"ExactCombinedGradient", true},
            {"CombinedPoyntingNormalization", true},
            {"ModePropagationConstantOutput", "mode-kn.csv"},
            {"CoefficientDiagnosticOutput", "singular-mode-diagnostics.csv"},
            {"CoefficientDiagnosticQuantity", "canonical basis coordinate norm"},
            {"CanonicalComplexCoefficientOutput",
             {{"ND", "singular-mode-nd-coefficients.csv"},
              {"H1", "singular-mode-h1-coefficients.csv"},
              {"PhaseConvention",
               "first canonical coefficient above 1e-12 of the mode maximum is "
               "positive real"},
              {"NDUnit", "V"},
              {"H1Unit", "V/m"},
              {"PhysicalAmplitude", false}}},
            {"TransverseFieldTipSlopeOutput",
             {{"File", "singular-mode-tip-slopes.csv"},
              {"Field", "complex magnitude |Et|"},
              {"Ray", "tip to opposite-edge midpoint in each triangle"},
              {"ExpectedSlope", "nu - 1"},
              {"PhysicalAmplitude", false}}},
            {"PhysicalAmplitude", false},
            {"CombinedFieldReconstruction", true},
            {"FieldEnergy",
             {{"TransverseElectric", true},
              {"NormalElectric", true},
              {"NormalMagnetic", true},
              {"TransverseMagnetic", true}}},
            {"FieldPostprocessing",
             {{"DiscontinuousParaViewSampling",
               iodata.problem.output_formats.paraview && bm.n_post > 0},
              {"MFEMGridFunction", false},
              {"VoltageAndImpedance", false},
              {"SurfaceMeasurements", !singular_surface_postoperator->Empty()},
              {"SurfaceQuadrature", "endpoint-weighted Gauss-Jacobi"}}},
            {"ErrorEstimator", false},
            {"SurfaceIntegrability",
             GetSingularSurfaceIntegrabilityMetadata(local_singular_features)},
            {"SurfaceParticipation", GetSingularSurfaceParticipationMetadata(iodata)}});
    Mpi::Warning(
        mode_op.GetComm(),
        "Singular BoundaryMode reconstructs combined ND/H1 fields and reports exact "
        "transverse/normal electric and magnetic energies. Integrable dielectric surface "
        "measurements use the combined field and endpoint-weighted quadrature. Voltage, "
        "impedance, MFEM grid-function output, and the standard flux error estimator "
        "remain disabled. Ideal-sheet surface traces with nu <= 1/2 require an explicit "
        "physical cutoff or response model.\n");
  }

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    // Load eigenvector i, apply the shared VD back-transform so en holds physical En.
    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    eig.GetEigenvector(i, e0);
    ComplexVector et, en;
    const std::complex<double> kn = eig.GetPropagationConstant(i);
    mode_op.ApplyVDBackTransform(e0, kn, et, en);

    // Power-normalize the mode to |P| = 1. PostOperator re-evaluates P on the
    // normalized mode for its impedance postprocessing.
    const std::complex<double> P_initial = mode_op.ComputePoyntingPower(omega, kn, et, en);
    if (std::abs(P_initial) > 0.0)
    {
      const double normalization = 1.0 / std::sqrt(std::abs(P_initial));
      et *= normalization;
      en *= normalization;
    }
    const std::complex<double> P_normalized =
        mode_op.ComputePoyntingPower(omega, kn, et, en);

    const double error_bkwd = eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    const double error_abs = eig.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    if (singular)
    {
      const auto coefficient_norms = mode_op.ComputeSingularCoefficientNorms(et, en);
      const auto field_energies = mode_op.ComputeSingularFieldEnergies(omega, kn, et, en);
      singular_et_real_evaluator->SetFromTrueDofs(et.Real());
      singular_et_imag_evaluator->SetFromTrueDofs(et.Imag());
      singular_en_real_evaluator->SetFromTrueDofs(en.Real());
      singular_en_imag_evaluator->SetFromTrueDofs(en.Imag());
      if (!singular_surface_postoperator->Empty())
      {
        singular_surface_measurements.push_back(
            {i + 1,
             singular_surface_postoperator->Measure(
                 *singular_et_real_evaluator, *singular_et_imag_evaluator,
                 field_energies.electric_transverse + field_energies.electric_normal,
                 {iodata.solver.singular_elements.quadrature_order,
                  iodata.solver.singular_elements.abs_tol,
                  iodata.solver.singular_elements.rel_tol,
                  iodata.solver.singular_elements.max_subdivisions},
                 singular_en_real_evaluator.get(), singular_en_imag_evaluator.get())});
      }
      const auto &numbering = mode_op.GetSingularDofNumbering();
      auto mode_coefficients = GatherModeCoefficientMeasurements(
          mode_op.GetComm(), i + 1,
          singular_et_real_evaluator->GetOwnedCoefficientDiagnostics(),
          singular_et_imag_evaluator->GetOwnedCoefficientDiagnostics(),
          numbering.nd.global_size,
          singular_en_real_evaluator->GetOwnedCoefficientDiagnostics(),
          singular_en_imag_evaluator->GetOwnedCoefficientDiagnostics(),
          numbering.h1.global_size);
      auto mode_tip_slopes = GatherModeTipSlopeMeasurements(
          mode_op.GetComm(), i + 1,
          singular_et_real_evaluator->FitComplexTipSlopes(
              *singular_et_imag_evaluator, local_singular_features, source_vertex_ids,
              source_element_ids));
      if (root)
      {
        singular_coefficient_measurements.insert(
            singular_coefficient_measurements.end(),
            std::make_move_iterator(mode_coefficients.begin()),
            std::make_move_iterator(mode_coefficients.end()));
        singular_tip_slope_measurements.insert(
            singular_tip_slope_measurements.end(),
            std::make_move_iterator(mode_tip_slopes.begin()),
            std::make_move_iterator(mode_tip_slopes.end()));
      }
      if (singular_paraview && i < bm.n_post)
      {
        singular_paraview->Write(*singular_et_real_evaluator, *singular_et_imag_evaluator,
                                 *singular_en_real_evaluator, *singular_en_imag_evaluator,
                                 kn, omega, i);
        Mpi::Print(" Wrote singular mode {:d} to disk (interior-sampled Paraview)\n",
                   i + 1);
      }
      singular_measurements.push_back({i + 1, kn, kn / omega, error_bkwd, error_abs,
                                       P_normalized, coefficient_norms, field_energies});
      continue;
    }

    auto total_domain_energy =
        post_op->MeasureAndPrintAll(i, et, en, kn, omega, error_abs, error_bkwd, n_print);

    if (i < num_modes && ModeEigenSolver::IsPropagating(kn))
    {
      // Bz = curl(Et) / (i·omega) = (Im(curl Et) - i·Re(curl Et)) / omega.
      const int l2_size = mode_op.GetCurlSpace().GetTrueVSize();
      ComplexVector bz(l2_size);
      bz.UseDevice(true);
      const auto &CurlOp =
          mode_op.GetCurlSpace().GetDiscreteInterpolator(mode_op.GetNDSpace());
      Vector curl_etr(l2_size), curl_eti(l2_size);
      curl_etr.UseDevice(true);
      curl_eti.UseDevice(true);
      CurlOp.Mult(et.Real(), curl_etr);
      CurlOp.Mult(et.Imag(), curl_eti);
      bz.Real() = curl_eti;
      bz.Real() *= 1.0 / omega;
      bz.Imag() = curl_etr;
      bz.Imag() *= -1.0 / omega;
      estimator->AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  if (singular)
  {
    WriteSingularModeMeasurements(post_dir, iodata, singular_measurements, root);
    WriteSingularModeCoefficientMeasurements(post_dir, iodata,
                                             singular_coefficient_measurements, root);
    WriteSingularModeTipSlopeMeasurements(post_dir, iodata, singular_tip_slope_measurements,
                                          root);
    WriteSingularModeSurfaceMeasurements(post_dir, singular_surface_measurements, root);
  }
  else
  {
    post_op->MeasureFinalize(indicator);
  }

  return {indicator, mode_op.GetNDGlobalTrueVSize() + mode_op.GetH1GlobalTrueVSize()};
}

}  // namespace palace
