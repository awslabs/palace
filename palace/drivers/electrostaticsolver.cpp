// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/singularsolver.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/laplaceoperator.hpp"
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

struct SingularDomainEnergyMeasurement
{
  int source;
  double total_energy;
  std::map<int, double> domain_energies;
  double algebraic_integral;
  double direct_integral;
  double integration_error_bound;
  double assembly_error_bound;
  double floating_point_allowance;
};

struct SingularCoefficientMeasurement
{
  int source;
  fem::singular::H1CoefficientDiagnostic diagnostic;
};

struct SingularEdgeSlopeMeasurement
{
  int source;
  fem::singular::H1EdgeSlopeDiagnostic diagnostic;
};

struct SingularTipSlopeMeasurement
{
  int source;
  fem::singular::H1TipSlopeDiagnostic diagnostic;
};

struct SingularSurfaceMeasurement
{
  int source;
  std::vector<TriangleSingularSurfacePostOperator::Measurement> interfaces;
};

constexpr std::size_t CoefficientIntegerFields = 23;
constexpr std::size_t CoefficientRealFields = 2;
using PackedCoefficientIntegers =
    std::array<fem::singular::GlobalVertexId, CoefficientIntegerFields>;
using PackedCoefficientReals = std::array<double, CoefficientRealFields>;
constexpr std::size_t EdgeSlopeIntegerFields = 10;
constexpr std::size_t EdgeSlopeRealFields = 8;
using PackedEdgeSlopeIntegers =
    std::array<fem::singular::GlobalVertexId, EdgeSlopeIntegerFields>;
using PackedEdgeSlopeReals = std::array<double, EdgeSlopeRealFields>;
constexpr std::size_t TipSlopeIntegerFields = 9;
constexpr std::size_t TipSlopeRealFields = 8;
using PackedTipSlopeIntegers =
    std::array<fem::singular::GlobalVertexId, TipSlopeIntegerFields>;
using PackedTipSlopeReals = std::array<double, TipSlopeRealFields>;

fem::singular::GlobalVertexId ToDiagnosticInteger(HYPRE_BigInt value)
{
  static_assert(std::numeric_limits<HYPRE_BigInt>::is_signed);
  static_assert(std::numeric_limits<fem::singular::GlobalVertexId>::is_signed);
  static_assert(std::numeric_limits<HYPRE_BigInt>::digits <=
                std::numeric_limits<fem::singular::GlobalVertexId>::digits);
  return static_cast<fem::singular::GlobalVertexId>(value);
}

HYPRE_BigInt ToHypreBigInt(fem::singular::GlobalVertexId value)
{
  if (value < std::numeric_limits<HYPRE_BigInt>::min() ||
      value > std::numeric_limits<HYPRE_BigInt>::max())
  {
    throw std::overflow_error("Singular diagnostic true DOF is out of range!");
  }
  return static_cast<HYPRE_BigInt>(value);
}

fem::singular::GlobalVertexId ToDiagnosticInteger(std::size_t value)
{
  if (value >
      static_cast<std::size_t>(std::numeric_limits<fem::singular::GlobalVertexId>::max()))
  {
    throw std::overflow_error("Singular diagnostic index is out of range!");
  }
  return static_cast<fem::singular::GlobalVertexId>(value);
}

int ToDiagnosticInt(fem::singular::GlobalVertexId value)
{
  if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
  {
    throw std::overflow_error("Singular diagnostic value is out of integer range!");
  }
  return static_cast<int>(value);
}

std::pair<PackedCoefficientIntegers, PackedCoefficientReals>
PackCoefficientMeasurement(const SingularCoefficientMeasurement &measurement)
{
  PackedCoefficientIntegers integers;
  std::size_t next = 0;
  integers[next++] = measurement.source;
  integers[next++] = ToDiagnosticInteger(measurement.diagnostic.true_dof);
  integers[next++] = static_cast<int>(measurement.diagnostic.key.family);
  integers[next++] = measurement.diagnostic.key.order;
  const auto pack_entity = [&integers, &next](const fem::singular::EntityKey &entity)
  {
    if (entity.size == 0 || entity.size > entity.vertices.size())
    {
      throw std::invalid_argument(
          "Singular coefficient diagnostic has an invalid canonical entity!");
    }
    integers[next++] = static_cast<fem::singular::GlobalVertexId>(entity.size);
    for (auto vertex : entity.vertices)
    {
      integers[next++] = vertex;
    }
  };
  pack_entity(measurement.diagnostic.key.singular_entity);
  pack_entity(measurement.diagnostic.key.support_entity);
  pack_entity(measurement.diagnostic.key.component_entity);
  for (int weight : measurement.diagnostic.key.interpolation_weights)
  {
    integers[next++] = weight;
  }
  if (next != integers.size() || !std::isfinite(measurement.diagnostic.exponent) ||
      !std::isfinite(measurement.diagnostic.coefficient))
  {
    throw std::invalid_argument("Singular coefficient diagnostic contains invalid data!");
  }
  return {integers, PackedCoefficientReals{measurement.diagnostic.exponent,
                                           measurement.diagnostic.coefficient}};
}

SingularCoefficientMeasurement
UnpackCoefficientMeasurement(const PackedCoefficientIntegers &integers,
                             const PackedCoefficientReals &reals)
{
  std::size_t next = 0;
  SingularCoefficientMeasurement measurement;
  measurement.source = ToDiagnosticInt(integers[next++]);
  measurement.diagnostic.true_dof = ToHypreBigInt(integers[next++]);
  const int family = ToDiagnosticInt(integers[next++]);
  if (family != static_cast<int>(fem::singular::HigherOrderBasisFamily::NODE_GRADIENT) &&
      family != static_cast<int>(fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT))
  {
    throw std::invalid_argument(
        "Singular coefficient diagnostic has an invalid H1 basis family!");
  }
  measurement.diagnostic.key.family =
      static_cast<fem::singular::HigherOrderBasisFamily>(family);
  measurement.diagnostic.key.order = ToDiagnosticInt(integers[next++]);
  const auto unpack_entity = [&integers, &next](fem::singular::EntityKey &entity)
  {
    const auto size = integers[next++];
    if (size <= 0 ||
        size > static_cast<fem::singular::GlobalVertexId>(entity.vertices.size()))
    {
      throw std::invalid_argument(
          "Singular coefficient diagnostic has an invalid canonical entity size!");
    }
    entity.size = static_cast<std::size_t>(size);
    for (auto &vertex : entity.vertices)
    {
      vertex = integers[next++];
    }
  };
  unpack_entity(measurement.diagnostic.key.singular_entity);
  unpack_entity(measurement.diagnostic.key.support_entity);
  unpack_entity(measurement.diagnostic.key.component_entity);
  for (auto &weight : measurement.diagnostic.key.interpolation_weights)
  {
    weight = ToDiagnosticInt(integers[next++]);
  }
  measurement.diagnostic.exponent = reals[0];
  measurement.diagnostic.coefficient = reals[1];
  if (next != integers.size() || measurement.diagnostic.key.order < 1 ||
      measurement.diagnostic.true_dof < 0 ||
      !std::isfinite(measurement.diagnostic.exponent) ||
      !(measurement.diagnostic.exponent > 0.0) ||
      !(measurement.diagnostic.exponent < 1.0) ||
      !std::isfinite(measurement.diagnostic.coefficient))
  {
    throw std::invalid_argument("Singular coefficient diagnostic unpacked invalid data!");
  }
  return measurement;
}

std::vector<SingularCoefficientMeasurement> GatherCoefficientMeasurements(
    MPI_Comm comm, int source,
    const std::vector<fem::singular::H1CoefficientDiagnostic> &local_diagnostics,
    HYPRE_BigInt expected_global_size)
{
  bool valid = local_diagnostics.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   std::max(CoefficientIntegerFields, CoefficientRealFields);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Singular coefficient diagnostics exceed MPI integer counts!");
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
        total_records > std::numeric_limits<int>::max() / CoefficientIntegerFields -
                            record_counts[rank])
    {
      throw std::overflow_error(
          "Gathered singular coefficient diagnostics exceed MPI integer counts!");
    }
    integer_displacements[rank] =
        static_cast<int>(total_records * CoefficientIntegerFields);
    real_displacements[rank] = static_cast<int>(total_records * CoefficientRealFields);
    integer_counts[rank] = static_cast<int>(record_counts[rank] * CoefficientIntegerFields);
    real_counts[rank] = static_cast<int>(record_counts[rank] * CoefficientRealFields);
    total_records += record_counts[rank];
  }

  std::vector<fem::singular::GlobalVertexId> local_integers(local_diagnostics.size() *
                                                            CoefficientIntegerFields);
  std::vector<double> local_reals(local_diagnostics.size() * CoefficientRealFields);
  for (std::size_t i = 0; i < local_diagnostics.size(); i++)
  {
    const auto [integers, reals] =
        PackCoefficientMeasurement({source, local_diagnostics[i]});
    std::copy(integers.begin(), integers.end(),
              local_integers.begin() + i * CoefficientIntegerFields);
    std::copy(reals.begin(), reals.end(), local_reals.begin() + i * CoefficientRealFields);
  }

  std::vector<fem::singular::GlobalVertexId> global_integers(
      static_cast<std::size_t>(total_records) * CoefficientIntegerFields);
  std::vector<double> global_reals(static_cast<std::size_t>(total_records) *
                                   CoefficientRealFields);
  Mpi::Allgatherv(static_cast<int>(local_integers.size()), local_integers.data(),
                  global_integers.data(), integer_counts.data(),
                  integer_displacements.data(), comm);
  Mpi::Allgatherv(static_cast<int>(local_reals.size()), local_reals.data(),
                  global_reals.data(), real_counts.data(), real_displacements.data(), comm);

  std::vector<SingularCoefficientMeasurement> measurements;
  measurements.reserve(static_cast<std::size_t>(total_records));
  for (std::size_t i = 0; i < static_cast<std::size_t>(total_records); i++)
  {
    PackedCoefficientIntegers integers;
    PackedCoefficientReals reals;
    std::copy(global_integers.begin() + i * CoefficientIntegerFields,
              global_integers.begin() + (i + 1) * CoefficientIntegerFields,
              integers.begin());
    std::copy(global_reals.begin() + i * CoefficientRealFields,
              global_reals.begin() + (i + 1) * CoefficientRealFields, reals.begin());
    measurements.push_back(UnpackCoefficientMeasurement(integers, reals));
  }
  std::sort(measurements.begin(), measurements.end(),
            [](const auto &left, const auto &right)
            { return left.diagnostic.key < right.diagnostic.key; });
  if (expected_global_size < 0 ||
      measurements.size() != static_cast<std::size_t>(expected_global_size))
  {
    throw std::runtime_error(
        "Singular coefficient diagnostics do not cover the global H1 enrichment!");
  }
  std::vector<HYPRE_BigInt> true_dofs;
  true_dofs.reserve(measurements.size());
  for (std::size_t i = 0; i < measurements.size(); i++)
  {
    if (measurements[i].source != source ||
        (i > 0 && !(measurements[i - 1].diagnostic.key < measurements[i].diagnostic.key)))
    {
      throw std::runtime_error(
          "Singular coefficient diagnostics have duplicate or invalid canonical keys!");
    }
    true_dofs.push_back(measurements[i].diagnostic.true_dof);
  }
  std::sort(true_dofs.begin(), true_dofs.end());
  for (std::size_t i = 0; i < true_dofs.size(); i++)
  {
    if (true_dofs[i] != static_cast<HYPRE_BigInt>(i))
    {
      throw std::runtime_error(
          "Singular coefficient diagnostics have duplicate or missing true DOFs!");
    }
  }
  return measurements;
}

std::pair<PackedEdgeSlopeIntegers, PackedEdgeSlopeReals>
PackEdgeSlopeMeasurement(const SingularEdgeSlopeMeasurement &measurement)
{
  const auto &diagnostic = measurement.diagnostic;
  PackedEdgeSlopeIntegers integers{measurement.source,
                                   diagnostic.source_element,
                                   ToDiagnosticInteger(diagnostic.feature),
                                   ToDiagnosticInteger(diagnostic.segment),
                                   diagnostic.canonical_vertices[0],
                                   diagnostic.canonical_vertices[1],
                                   diagnostic.canonical_vertices[2],
                                   diagnostic.canonical_vertices[3],
                                   diagnostic.sample_count,
                                   diagnostic.valid ? 1 : 0};
  PackedEdgeSlopeReals reals{diagnostic.exponent,
                             diagnostic.expected_slope,
                             diagnostic.fitted_slope,
                             diagnostic.r_squared,
                             diagnostic.minimum_distance,
                             diagnostic.maximum_distance,
                             diagnostic.field_norm_at_minimum_distance,
                             diagnostic.field_norm_at_maximum_distance};
  if (diagnostic.source_element < 0 || diagnostic.sample_count < 3 ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }))
  {
    throw std::invalid_argument("Singular edge slope diagnostic contains invalid data!");
  }
  return {integers, reals};
}

SingularEdgeSlopeMeasurement
UnpackEdgeSlopeMeasurement(const PackedEdgeSlopeIntegers &integers,
                           const PackedEdgeSlopeReals &reals)
{
  SingularEdgeSlopeMeasurement measurement;
  measurement.source = ToDiagnosticInt(integers[0]);
  auto &diagnostic = measurement.diagnostic;
  diagnostic.source_element = integers[1];
  if (integers[2] < 0 || integers[3] < 0)
  {
    throw std::invalid_argument(
        "Singular edge slope diagnostic has a negative feature index!");
  }
  diagnostic.feature = static_cast<std::size_t>(integers[2]);
  diagnostic.segment = static_cast<std::size_t>(integers[3]);
  std::copy(integers.begin() + 4, integers.begin() + 8,
            diagnostic.canonical_vertices.begin());
  diagnostic.sample_count = ToDiagnosticInt(integers[8]);
  diagnostic.valid = integers[9] != 0;
  diagnostic.exponent = reals[0];
  diagnostic.expected_slope = reals[1];
  diagnostic.fitted_slope = reals[2];
  diagnostic.r_squared = reals[3];
  diagnostic.minimum_distance = reals[4];
  diagnostic.maximum_distance = reals[5];
  diagnostic.field_norm_at_minimum_distance = reals[6];
  diagnostic.field_norm_at_maximum_distance = reals[7];
  if (diagnostic.source_element < 0 || diagnostic.sample_count < 3 ||
      (integers[9] != 0 && integers[9] != 1) ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }) ||
      !(diagnostic.exponent > 0.0) || !(diagnostic.exponent < 1.0) ||
      diagnostic.expected_slope != diagnostic.exponent - 1.0 ||
      (diagnostic.valid && (!(diagnostic.minimum_distance > 0.0) ||
                            !(diagnostic.maximum_distance > diagnostic.minimum_distance) ||
                            !(diagnostic.field_norm_at_minimum_distance > 0.0) ||
                            !(diagnostic.field_norm_at_maximum_distance > 0.0))))
  {
    throw std::invalid_argument("Singular edge slope diagnostic unpacked invalid data!");
  }
  return measurement;
}

std::vector<SingularEdgeSlopeMeasurement> GatherEdgeSlopeMeasurements(
    MPI_Comm comm, int source,
    const std::vector<fem::singular::H1EdgeSlopeDiagnostic> &local_diagnostics)
{
  bool valid = local_diagnostics.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   std::max(EdgeSlopeIntegerFields, EdgeSlopeRealFields);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error("Singular edge slope diagnostics exceed MPI integer counts!");
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
        total_records >
            std::numeric_limits<int>::max() / EdgeSlopeIntegerFields - record_counts[rank])
    {
      throw std::overflow_error(
          "Gathered singular edge slope diagnostics exceed MPI integer counts!");
    }
    integer_displacements[rank] = static_cast<int>(total_records * EdgeSlopeIntegerFields);
    real_displacements[rank] = static_cast<int>(total_records * EdgeSlopeRealFields);
    integer_counts[rank] = static_cast<int>(record_counts[rank] * EdgeSlopeIntegerFields);
    real_counts[rank] = static_cast<int>(record_counts[rank] * EdgeSlopeRealFields);
    total_records += record_counts[rank];
  }

  std::vector<fem::singular::GlobalVertexId> local_integers(local_diagnostics.size() *
                                                            EdgeSlopeIntegerFields);
  std::vector<double> local_reals(local_diagnostics.size() * EdgeSlopeRealFields);
  for (std::size_t i = 0; i < local_diagnostics.size(); i++)
  {
    const auto [integers, reals] = PackEdgeSlopeMeasurement({source, local_diagnostics[i]});
    std::copy(integers.begin(), integers.end(),
              local_integers.begin() + i * EdgeSlopeIntegerFields);
    std::copy(reals.begin(), reals.end(), local_reals.begin() + i * EdgeSlopeRealFields);
  }

  std::vector<fem::singular::GlobalVertexId> global_integers(
      static_cast<std::size_t>(total_records) * EdgeSlopeIntegerFields);
  std::vector<double> global_reals(static_cast<std::size_t>(total_records) *
                                   EdgeSlopeRealFields);
  Mpi::Allgatherv(static_cast<int>(local_integers.size()), local_integers.data(),
                  global_integers.data(), integer_counts.data(),
                  integer_displacements.data(), comm);
  Mpi::Allgatherv(static_cast<int>(local_reals.size()), local_reals.data(),
                  global_reals.data(), real_counts.data(), real_displacements.data(), comm);

  std::vector<SingularEdgeSlopeMeasurement> measurements;
  measurements.reserve(static_cast<std::size_t>(total_records));
  for (std::size_t i = 0; i < static_cast<std::size_t>(total_records); i++)
  {
    PackedEdgeSlopeIntegers integers;
    PackedEdgeSlopeReals reals;
    std::copy(global_integers.begin() + i * EdgeSlopeIntegerFields,
              global_integers.begin() + (i + 1) * EdgeSlopeIntegerFields, integers.begin());
    std::copy(global_reals.begin() + i * EdgeSlopeRealFields,
              global_reals.begin() + (i + 1) * EdgeSlopeRealFields, reals.begin());
    measurements.push_back(UnpackEdgeSlopeMeasurement(integers, reals));
  }
  std::sort(measurements.begin(), measurements.end(),
            [](const auto &left, const auto &right)
            {
              return std::tie(left.diagnostic.source_element, left.diagnostic.segment,
                              left.diagnostic.feature) <
                     std::tie(right.diagnostic.source_element, right.diagnostic.segment,
                              right.diagnostic.feature);
            });
  if (measurements.empty())
  {
    throw std::runtime_error(
        "Singular edge slope diagnostics found no global edge sectors!");
  }
  for (std::size_t i = 0; i < measurements.size(); i++)
  {
    if (measurements[i].source != source ||
        (i > 0 && std::tie(measurements[i - 1].diagnostic.source_element,
                           measurements[i - 1].diagnostic.segment) ==
                      std::tie(measurements[i].diagnostic.source_element,
                               measurements[i].diagnostic.segment)))
    {
      throw std::runtime_error(
          "Singular edge slope diagnostics have duplicate element sectors!");
    }
  }
  return measurements;
}

std::pair<PackedTipSlopeIntegers, PackedTipSlopeReals>
PackTipSlopeMeasurement(const SingularTipSlopeMeasurement &measurement)
{
  const auto &diagnostic = measurement.diagnostic;
  PackedTipSlopeIntegers integers{measurement.source,
                                  diagnostic.source_element,
                                  ToDiagnosticInteger(diagnostic.feature),
                                  ToDiagnosticInteger(diagnostic.selected_segment),
                                  diagnostic.canonical_vertices[0],
                                  diagnostic.canonical_vertices[1],
                                  diagnostic.canonical_vertices[2],
                                  diagnostic.sample_count,
                                  diagnostic.valid ? 1 : 0};
  PackedTipSlopeReals reals{diagnostic.exponent,
                            diagnostic.expected_slope,
                            diagnostic.fitted_slope,
                            diagnostic.r_squared,
                            diagnostic.minimum_distance,
                            diagnostic.maximum_distance,
                            diagnostic.field_norm_at_minimum_distance,
                            diagnostic.field_norm_at_maximum_distance};
  if (diagnostic.source_element < 0 || diagnostic.sample_count < 3 ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }))
  {
    throw std::invalid_argument("Singular tip slope diagnostic contains invalid data!");
  }
  return {integers, reals};
}

SingularTipSlopeMeasurement
UnpackTipSlopeMeasurement(const PackedTipSlopeIntegers &integers,
                          const PackedTipSlopeReals &reals)
{
  SingularTipSlopeMeasurement measurement;
  measurement.source = ToDiagnosticInt(integers[0]);
  auto &diagnostic = measurement.diagnostic;
  diagnostic.source_element = integers[1];
  if (integers[2] < 0 || integers[3] < 0)
  {
    throw std::invalid_argument(
        "Singular tip slope diagnostic has a negative feature index!");
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
  if (diagnostic.source_element < 0 || diagnostic.sample_count < 3 ||
      (integers[8] != 0 && integers[8] != 1) ||
      !std::all_of(
          reals.begin(), reals.end(), [](double value) { return std::isfinite(value); }) ||
      !(diagnostic.exponent > 0.0) || !(diagnostic.exponent < 1.0) ||
      diagnostic.expected_slope != diagnostic.exponent - 1.0 ||
      (diagnostic.valid && (!(diagnostic.minimum_distance > 0.0) ||
                            !(diagnostic.maximum_distance > diagnostic.minimum_distance) ||
                            !(diagnostic.field_norm_at_minimum_distance > 0.0) ||
                            !(diagnostic.field_norm_at_maximum_distance > 0.0))))
  {
    throw std::invalid_argument("Singular tip slope diagnostic unpacked invalid data!");
  }
  return measurement;
}

std::vector<SingularTipSlopeMeasurement> GatherTipSlopeMeasurements(
    MPI_Comm comm, int source,
    const std::vector<fem::singular::H1TipSlopeDiagnostic> &local_diagnostics)
{
  bool valid = local_diagnostics.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   std::max(TipSlopeIntegerFields, TipSlopeRealFields);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error("Singular tip slope diagnostics exceed MPI integer counts!");
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
        total_records >
            std::numeric_limits<int>::max() / TipSlopeIntegerFields - record_counts[rank])
    {
      throw std::overflow_error(
          "Gathered singular tip slope diagnostics exceed MPI integer counts!");
    }
    integer_displacements[rank] = static_cast<int>(total_records * TipSlopeIntegerFields);
    real_displacements[rank] = static_cast<int>(total_records * TipSlopeRealFields);
    integer_counts[rank] = static_cast<int>(record_counts[rank] * TipSlopeIntegerFields);
    real_counts[rank] = static_cast<int>(record_counts[rank] * TipSlopeRealFields);
    total_records += record_counts[rank];
  }

  std::vector<fem::singular::GlobalVertexId> local_integers(local_diagnostics.size() *
                                                            TipSlopeIntegerFields);
  std::vector<double> local_reals(local_diagnostics.size() * TipSlopeRealFields);
  for (std::size_t i = 0; i < local_diagnostics.size(); i++)
  {
    const auto [integers, reals] = PackTipSlopeMeasurement({source, local_diagnostics[i]});
    std::copy(integers.begin(), integers.end(),
              local_integers.begin() + i * TipSlopeIntegerFields);
    std::copy(reals.begin(), reals.end(), local_reals.begin() + i * TipSlopeRealFields);
  }

  std::vector<fem::singular::GlobalVertexId> global_integers(
      static_cast<std::size_t>(total_records) * TipSlopeIntegerFields);
  std::vector<double> global_reals(static_cast<std::size_t>(total_records) *
                                   TipSlopeRealFields);
  Mpi::Allgatherv(static_cast<int>(local_integers.size()), local_integers.data(),
                  global_integers.data(), integer_counts.data(),
                  integer_displacements.data(), comm);
  Mpi::Allgatherv(static_cast<int>(local_reals.size()), local_reals.data(),
                  global_reals.data(), real_counts.data(), real_displacements.data(), comm);

  std::vector<SingularTipSlopeMeasurement> measurements;
  measurements.reserve(static_cast<std::size_t>(total_records));
  for (std::size_t i = 0; i < static_cast<std::size_t>(total_records); i++)
  {
    PackedTipSlopeIntegers integers;
    PackedTipSlopeReals reals;
    std::copy(global_integers.begin() + i * TipSlopeIntegerFields,
              global_integers.begin() + (i + 1) * TipSlopeIntegerFields, integers.begin());
    std::copy(global_reals.begin() + i * TipSlopeRealFields,
              global_reals.begin() + (i + 1) * TipSlopeRealFields, reals.begin());
    measurements.push_back(UnpackTipSlopeMeasurement(integers, reals));
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
        "Singular tip slope diagnostics found no global triangle sectors!");
  }
  for (std::size_t i = 0; i < measurements.size(); i++)
  {
    if (measurements[i].source != source ||
        (i > 0 && std::tie(measurements[i - 1].diagnostic.source_element,
                           measurements[i - 1].diagnostic.feature) ==
                      std::tie(measurements[i].diagnostic.source_element,
                               measurements[i].diagnostic.feature)))
    {
      throw std::runtime_error(
          "Singular tip slope diagnostics have duplicate triangle sectors!");
    }
  }
  return measurements;
}

class SingularParaviewOutput
{
private:
  mfem::ParMesh &mesh;
  mfem::L2_FECollection sample_collection;
  mfem::ParFiniteElementSpace potential_space;
  mfem::ParFiniteElementSpace electric_space;
  mfem::ParGridFunction potential;
  mfem::ParGridFunction electric_field;
  mfem::ParaViewDataCollection data_collection;
  const Units &units;

public:
  SingularParaviewOutput(const fs::path &post_dir, mfem::ParMesh &mesh, int order,
                         const Units &units)
    : mesh(mesh), sample_collection(order, mesh.Dimension(), mfem::BasisType::GaussLegendre,
                                    mfem::FiniteElement::VALUE),
      potential_space(&mesh, &sample_collection),
      electric_space(&mesh, &sample_collection, mesh.SpaceDimension(),
                     mfem::Ordering::byVDIM),
      potential(&potential_space), electric_field(&electric_space),
      data_collection(post_dir / "paraview" / "electrostatic", &mesh), units(units)
  {
    MFEM_VERIFY(order >= 1,
                "Singular field visualization requires positive sampling order!");
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
    data_collection.RegisterField("V", &potential);
    data_collection.RegisterField("E", &electric_field);
  }

  template <typename Evaluator>
  void Write(Evaluator &evaluator, int step, int source)
  {
    evaluator.ProjectToDiscontinuousGridFunctions(potential, electric_field);

    const double mesh_scale = units.GetMeshLengthRelativeScale();
    const double potential_scale = units.Dimensionalize<Units::ValueType::VOLTAGE>(1.0);
    const double electric_scale = units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
    mesh::DimensionalizeMesh(mesh, mesh_scale);
    potential *= potential_scale;
    electric_field *= electric_scale;
    data_collection.SetCycle(step);
    data_collection.SetTime(source);
    data_collection.Save();
    potential *= 1.0 / potential_scale;
    electric_field *= 1.0 / electric_scale;
    mesh::NondimensionalizeMesh(mesh, mesh_scale);
    Mpi::Barrier(mesh.GetComm());
  }
};

template <typename Evaluator>
SingularDomainEnergyMeasurement
MeasureSingularDomainEnergy(const IoData &iodata, const LaplaceOperator &laplace_op,
                            Evaluator &evaluator, const Vector &potential, int source)
{
  const fem::singular::AdaptiveAssemblyOptions options{
      iodata.solver.singular_elements.quadrature_order,
      iodata.solver.singular_elements.abs_tol, iodata.solver.singular_elements.rel_tol,
      iodata.solver.singular_elements.max_subdivisions};
  const auto &configured_domains = iodata.domains.postpro.energy;
  std::vector<double> reductions(2 + configured_domains.size(), 0.0);

  for (int element = 0; element < laplace_op.GetMesh().GetNE(); element++)
  {
    const int attribute = laplace_op.GetMesh().Get().GetAttribute(element);
    const auto &material = laplace_op.GetMaterialOp();
    MFEM_VERIFY(
        material.IsPermittivityIsotropic(attribute),
        "Singular electrostatic field energy requires isotropic permittivity in every "
        "domain! Domain attribute: "
            << attribute);
    const double epsilon = material.GetPermittivityReal(attribute)(0, 0);
    const auto integral =
        evaluator.IntegrateElementGradientEnergy(element, epsilon, options);
    MFEM_VERIFY(integral.converged && std::isfinite(integral.value) &&
                    std::isfinite(integral.estimated_absolute_error) &&
                    integral.estimated_absolute_error >= 0.0,
                "Combined singular field energy integration failed to converge on "
                "element "
                    << element << "!");
    reductions[0] += integral.value;
    reductions[1] += integral.estimated_absolute_error;

    std::size_t domain = 0;
    for (const auto &[idx, data] : configured_domains)
    {
      if (std::binary_search(data.attributes.begin(), data.attributes.end(), attribute))
      {
        reductions[2 + domain] += integral.value;
      }
      domain++;
    }
  }
  Mpi::GlobalSum(static_cast<int>(reductions.size()), reductions.data(),
                 laplace_op.GetComm());

  const auto &stiffness = laplace_op.GetUnconstrainedStiffnessMatrix();
  Vector stiffness_potential(stiffness.Height());
  stiffness.Mult(potential, stiffness_potential);
  const double algebraic_integral =
      linalg::Dot(laplace_op.GetComm(), potential, stiffness_potential);
  const double assembly_error = laplace_op.GetSingularStiffnessEnergyErrorBound(potential);
  const double floating_point_error =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, static_cast<double>(laplace_op.GlobalTrueVSize())) *
      (std::abs(reductions[0]) + std::abs(algebraic_integral));
  const double discrepancy = std::abs(reductions[0] - algebraic_integral);
  const double error_bound = reductions[1] + assembly_error + floating_point_error;
  MFEM_VERIFY(std::isfinite(algebraic_integral) && algebraic_integral > 0.0 &&
                  std::isfinite(discrepancy) && std::isfinite(error_bound) &&
                  discrepancy <= error_bound,
              "Direct combined-field energy is inconsistent with V^T K V! "
              "Discrepancy: "
                  << discrepancy << ", integration error: " << reductions[1]
                  << ", assembly error: " << assembly_error
                  << ", floating-point allowance: " << floating_point_error);

  SingularDomainEnergyMeasurement measurement{source,
                                              0.5 * reductions[0],
                                              {},
                                              algebraic_integral,
                                              reductions[0],
                                              reductions[1],
                                              assembly_error,
                                              floating_point_error};
  std::size_t domain = 0;
  for (const auto &[idx, data] : configured_domains)
  {
    measurement.domain_energies.emplace(idx, 0.5 * reductions[2 + domain]);
    domain++;
  }
  Mpi::Print(
      " Field energy E = {:.3e} J (combined-field relative discrepancy = {:.3e}, "
      "relative bound = {:.3e})\n",
      iodata.units.Dimensionalize<Units::ValueType::ENERGY>(measurement.total_energy),
      discrepancy / algebraic_integral, error_bound / algebraic_integral);
  return measurement;
}

void WriteSingularDomainEnergy(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularDomainEnergyMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "domain-E.csv");
  output.table.insert(Column("i", "i", 0, 0, 2, ""));
  output.table.insert("Ee", "E_elec (J)");
  output.table.insert("Em", "E_mag (J)");
  output.table.insert("Ec", "E_cap (J)");
  output.table.insert("Ei", "E_ind (J)");
  for (const auto &[idx, data] : iodata.domains.postpro.energy)
  {
    output.table.insert(fmt::format("Ee_{}", idx), fmt::format("E_elec[{}] (J)", idx));
    output.table.insert(fmt::format("pe_{}", idx), fmt::format("p_elec[{}]", idx));
    output.table.insert(fmt::format("Em_{}", idx), fmt::format("E_mag[{}] (J)", idx));
    output.table.insert(fmt::format("pm_{}", idx), fmt::format("p_mag[{}]", idx));
  }

  const double energy_scale = iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0);
  for (const auto &measurement : measurements)
  {
    output.table["i"] << measurement.source;
    output.table["Ee"] << energy_scale * measurement.total_energy;
    output.table["Em"] << 0.0;
    output.table["Ec"] << 0.0;
    output.table["Ei"] << 0.0;
    for (const auto &[idx, data] : iodata.domains.postpro.energy)
    {
      const double energy = measurement.domain_energies.at(idx);
      const double participation =
          measurement.total_energy > 0.0 ? energy / measurement.total_energy : 0.0;
      output.table[fmt::format("Ee_{}", idx)] << energy_scale * energy;
      output.table[fmt::format("pe_{}", idx)] << participation;
      output.table[fmt::format("Em_{}", idx)] << 0.0;
      output.table[fmt::format("pm_{}", idx)] << 0.0;
    }
  }
  output.WriteFullTableTrunc();
}

void WriteSingularSurfaceMeasurements(
    const fs::path &post_dir, const std::vector<SingularSurfaceMeasurement> &measurements,
    bool root)
{
  if (!root || measurements.empty() || measurements.front().interfaces.empty())
  {
    return;
  }
  TableWithCSVFile output(post_dir / "surface-Q.csv");
  output.table.insert(Column("i", "i", 0, 0, 2, ""));
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
          "Singular electrostatic surface measurements changed interface count!");
    }
    output.table["i"] << measurement.source;
    for (std::size_t i = 0; i < measurement.interfaces.size(); i++)
    {
      const auto &interface = measurement.interfaces[i];
      if (interface.index != measurements.front().interfaces[i].index)
      {
        throw std::runtime_error(
            "Singular electrostatic surface measurements changed interface ordering!");
      }
      output.table[fmt::format("p_{}", interface.index)] << interface.participation;
      output.table[fmt::format("Q_{}", interface.index)] << interface.quality_factor;
    }
  }
  output.WriteFullTableTrunc();
}

void InsertIntegerDiagnosticColumn(Table &table, std::string name, std::string header)
{
  table.insert(Column(std::move(name), std::move(header), 0, 0, 0, ""));
  table[table.n_cols() - 1].print_as_int = true;
}

void WriteSingularCoefficientMeasurements(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularCoefficientMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "singular-coefficients.csv");
  InsertIntegerDiagnosticColumn(output.table, "source", "source");
  InsertIntegerDiagnosticColumn(output.table, "true_dof", "true_dof");
  InsertIntegerDiagnosticColumn(output.table, "family", "family (0 = node; 1 = edge)");
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
  output.table.insert("coefficient", "basis_coefficient (V)");

  const double voltage_scale = iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(1.0);
  for (const auto &measurement : measurements)
  {
    const auto &diagnostic = measurement.diagnostic;
    output.table["source"] << measurement.source;
    output.table["true_dof"] << static_cast<double>(diagnostic.true_dof);
    output.table["family"] << ((diagnostic.key.family ==
                                fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
                                   ? 0.0
                                   : 1.0);
    output.table["order"] << diagnostic.key.order;
    output.table["nu"] << diagnostic.exponent;
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
    write_entity("singular", diagnostic.key.singular_entity);
    write_entity("support", diagnostic.key.support_entity);
    write_entity("component", diagnostic.key.component_entity);
    for (int weight = 0; weight < 4; weight++)
    {
      output.table[fmt::format("weight_{}", weight)]
          << diagnostic.key.interpolation_weights[weight];
    }
    output.table["coefficient"] << voltage_scale * diagnostic.coefficient;
  }
  output.WriteFullTableTrunc();
}

void WriteSingularEdgeSlopeMeasurements(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularEdgeSlopeMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "singular-edge-slopes.csv");
  for (const auto &[name, header] :
       std::array<std::pair<std::string_view, std::string_view>, 10>{
           std::pair{"source", "source"}, std::pair{"element", "source_element"},
           std::pair{"feature", "feature"}, std::pair{"segment", "segment"},
           std::pair{"edge_v0", "edge_vertex_0"}, std::pair{"edge_v1", "edge_vertex_1"},
           std::pair{"opposite_v0", "opposite_vertex_0"},
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
  output.table.insert("field_min", "|E| at minimum_distance (V/m)");
  output.table.insert("field_max", "|E| at maximum_distance (V/m)");

  const double length_scale = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
  const double field_scale = iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  for (const auto &measurement : measurements)
  {
    const auto &diagnostic = measurement.diagnostic;
    output.table["source"] << measurement.source;
    output.table["element"] << static_cast<double>(diagnostic.source_element);
    output.table["feature"] << static_cast<double>(diagnostic.feature);
    output.table["segment"] << static_cast<double>(diagnostic.segment);
    output.table["edge_v0"] << static_cast<double>(diagnostic.canonical_vertices[0]);
    output.table["edge_v1"] << static_cast<double>(diagnostic.canonical_vertices[1]);
    output.table["opposite_v0"] << static_cast<double>(diagnostic.canonical_vertices[2]);
    output.table["opposite_v1"] << static_cast<double>(diagnostic.canonical_vertices[3]);
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

void WriteSingularTipSlopeMeasurements(
    const fs::path &post_dir, const IoData &iodata,
    const std::vector<SingularTipSlopeMeasurement> &measurements, bool root)
{
  if (!root)
  {
    return;
  }
  TableWithCSVFile output(post_dir / "singular-tip-slopes.csv");
  for (const auto &[name, header] :
       std::array<std::pair<std::string_view, std::string_view>, 9>{
           std::pair{"source", "source"}, std::pair{"element", "source_element"},
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
  output.table.insert("field_min", "|E| at minimum_distance (V/m)");
  output.table.insert("field_max", "|E| at maximum_distance (V/m)");

  const double length_scale = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
  const double field_scale = iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  for (const auto &measurement : measurements)
  {
    const auto &diagnostic = measurement.diagnostic;
    output.table["source"] << measurement.source;
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

}  // namespace

void ElectrostaticSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                                     MPI_Comm comm) const
{
  BaseSolver::Preprocess(iodata, smesh, comm);
  if (!iodata.solver.singular_elements.Enabled())
  {
    return;
  }

  serial_singular_features = {};
  local_singular_features = {};
  serial_triangle_singular_features = {};
  local_triangle_singular_features = {};
  source_vertex_ids.clear();
  source_element_ids.clear();
  int mesh_dimension = 0;
  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(smesh, "Root rank has no serial mesh for singular feature extraction!");
    mesh_dimension = smesh->Dimension();
    if (mesh_dimension == 3)
    {
      serial_singular_features = fem::singular::ExtractSerialSheetFeatures(
          *smesh, iodata.solver.singular_elements.attributes,
          GetSingularTriangleMaterials(iodata));
    }
    else
    {
      MFEM_VERIFY(smesh->Dimension() == 2 && smesh->SpaceDimension() == 2,
                  "Singular electrostatic enrichment requires a 2D triangular or 3D "
                  "tetrahedral mesh!");
      MFEM_VERIFY(iodata.solver.singular_elements.order == 1,
                  "Triangular singular electrostatic enrichment currently supports only "
                  "Solver.SingularElements.Order = 1!");
      serial_triangle_singular_features = fem::singular::ExtractSerialLineFeatures(
          *smesh, iodata.solver.singular_elements.attributes,
          GetSingularTriangleMaterials(iodata));
    }
  }
  Mpi::Broadcast(1, &mesh_dimension, 0, comm);
  if (mesh_dimension == 2)
  {
    fem::singular::BroadcastSerialLineTipFeatures(serial_triangle_singular_features, comm);
    MFEM_VERIFY(!serial_triangle_singular_features.Empty(),
                "Singular feature extraction produced no PEC line tips or corners!");
  }
  else
  {
    fem::singular::BroadcastSerialSheetFeatures(serial_singular_features, comm);
    MFEM_VERIFY(!serial_singular_features.Empty(),
                "Singular feature extraction produced no sheet-edge features!");
  }
}

bool ElectrostaticSolver::RequiresSourceSerialMeshMetadata() const
{
  return iodata.solver.singular_elements.Enabled();
}

void ElectrostaticSolver::ProcessPartitionedMesh(
    const mfem::ParMesh &parallel_mesh, const mesh::PartitionMetadata &metadata) const
{
  MFEM_VERIFY(iodata.solver.singular_elements.Enabled(),
              "Unexpected singular partition metadata for an unenriched solve!");
  if (parallel_mesh.Dimension() == 2)
  {
    MFEM_VERIFY(!serial_triangle_singular_features.Empty(),
                "Triangular singular source feature blueprint was not initialized!");
    local_triangle_singular_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_triangle_singular_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  else
  {
    MFEM_VERIFY(!serial_singular_features.Empty(),
                "Tetrahedral singular source feature blueprint was not initialized!");
    local_singular_features = fem::singular::DistributeSerialSheetFeatures(
        serial_singular_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  source_vertex_ids = metadata.source_vertex_ids;
  source_element_ids = metadata.source_element_ids;
}

mfem::Array<int>
ElectrostaticSolver::GetRefinementProtection(const mfem::ParMesh &mesh) const
{
  if (!iodata.solver.singular_elements.Enabled())
  {
    return {};
  }
  return mesh.Dimension() == 2
             ? BuildSingularRefinementProtection(mesh, local_triangle_singular_features,
                                                 source_vertex_ids)
             : BuildSingularRefinementProtection(mesh, local_singular_features,
                                                 source_vertex_ids);
}

void ElectrostaticSolver::ProcessRefinedMesh(const mfem::ParMesh &mesh) const
{
  if (!iodata.solver.singular_elements.Enabled())
  {
    return;
  }
  UpdateSingularSourceEntityIds(mesh, source_vertex_ids, source_element_ids);
  ProcessPartitionedMesh(mesh, {source_vertex_ids, source_element_ids});
}

bool ElectrostaticSolver::RebalanceRefinedMesh() const
{
  return !iodata.solver.singular_elements.Enabled();
}

std::pair<ErrorIndicator, long long int>
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  const bool singular = iodata.solver.singular_elements.Enabled();
  const bool triangle_singular = singular && mesh.back()->Dimension() == 2;
  MFEM_VERIFY(!singular ||
                  ((!triangle_singular && !local_singular_features.Empty()) ||
                   (triangle_singular && !local_triangle_singular_features.Empty())) &&
                      source_vertex_ids.size() ==
                          static_cast<std::size_t>(mesh.back()->Get().GetNV()) &&
                      source_element_ids.size() ==
                          static_cast<std::size_t>(mesh.back()->Get().GetNE()),
              "Singular feature topology was not reconstructed on the solve mesh!");
  LaplaceOperator laplace_op(
      iodata, mesh, singular && !triangle_singular ? &local_singular_features : nullptr,
      singular ? &source_vertex_ids : nullptr,
      triangle_singular ? &local_triangle_singular_features : nullptr);
  auto K = laplace_op.GetStiffnessMatrix();
  const auto &Grad = laplace_op.GetGradMatrix();
  SaveMetadata(laplace_op.GetH1Spaces());

  // Set up the linear solver.
  std::unique_ptr<KspSolver> ksp;
  if (singular)
  {
    ksp = MakeSingularPatchKspSolver(iodata, laplace_op.GetH1Spaces(), *K,
                                     laplace_op.GetSingularStandardStiffnessMatrix(),
                                     laplace_op.GetSingularFeaturePatches());
  }
  else
  {
    ksp = std::make_unique<KspSolver>(iodata, laplace_op.GetH1Spaces());
  }
  ksp->SetOperators(*K, *K);

  if (singular)
  {
    const int n_step = static_cast<int>(laplace_op.GetSources().size());
    MFEM_VERIFY(n_step > 0,
                "No terminal boundaries specified for electrostatic simulation!");

    Vector RHS(Grad.Width()), E(Grad.Height());
    std::vector<Vector> V(n_step);
    std::vector<SingularDomainEnergyMeasurement> domain_energy;
    std::vector<SingularCoefficientMeasurement> coefficient_measurements;
    std::vector<SingularEdgeSlopeMeasurement> edge_slope_measurements;
    std::vector<SingularTipSlopeMeasurement> tip_slope_measurements;
    std::vector<SingularSurfaceMeasurement> surface_measurements;
    GradFluxErrorEstimator estimator(
        laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
    ErrorIndicator indicator;
    const fem::singular::EdgeSlopeOptions edge_slope_options;
    domain_energy.reserve(n_step);
    const auto &singular_diagnostics = laplace_op.GetSingularDiagnostics();
    nlohmann::json singular_metadata{
        {"ConventionVersion", singular_diagnostics.convention_version},
        {"StandardOrder", singular_diagnostics.standard_order},
        {"SingularOrder", singular_diagnostics.singular_order},
        {"H1EnrichmentDegreesOfFreedom", singular_diagnostics.h1_enrichment_dofs},
        {"NDEnrichmentDegreesOfFreedom", singular_diagnostics.nd_enrichment_dofs},
        {"Quadrature",
         {{"Rule", "cached feature-aligned and partitioned Duffy reference tables"},
          {"AssemblyScope", "H1 electrostatic diffusion"},
          {"DuffyReferenceOrder", fem::singular::H1DuffyReferenceOrder},
          {"DuffyComparisonOrder", fem::singular::H1DuffyComparisonOrder},
          {"DuffyRadialPower", fem::singular::H1DuffyRadialPower},
          {"DuffyErrorSafetyFactor", fem::singular::H1DuffyErrorSafetyFactor},
          {"DuffyPartitionPower", fem::singular::MultiFeatureDuffyPartitionPower},
          {"DuffyMaximumTableEntriesPerRank",
           singular_diagnostics.duffy_reference_table_maximum_entries},
          {"DuffyTotalCacheHits", singular_diagnostics.duffy_reference_cache_hits},
          {"Order", singular_diagnostics.quadrature_order},
          {"AbsoluteTolerance", singular_diagnostics.quadrature_absolute_tolerance},
          {"RelativeTolerance", singular_diagnostics.quadrature_relative_tolerance},
          {"MaximumSubdivisions", singular_diagnostics.quadrature_maximum_subdivisions},
          {"TotalLeafCount", singular_diagnostics.quadrature_leaf_count},
          {"MaximumDepth", singular_diagnostics.quadrature_maximum_depth}}},
        {"StiffnessDiagonalScaling",
         {{"Metric", "maximum positive diagonal / minimum positive diagonal"},
          {"StandardMinimum", singular_diagnostics.standard_diagonal_minimum},
          {"StandardMaximum", singular_diagnostics.standard_diagonal_maximum},
          {"StandardSpread", singular_diagnostics.standard_diagonal_spread},
          {"EnrichmentMinimum", singular_diagnostics.enrichment_diagonal_minimum},
          {"EnrichmentMaximum", singular_diagnostics.enrichment_diagonal_maximum},
          {"EnrichmentSpread", singular_diagnostics.enrichment_diagonal_spread},
          {"CombinedSpread", singular_diagnostics.combined_diagonal_spread}}},
        {"LinearPreconditioner",
         {{"Type", UsesSingularPatchKspSolver(iodata)
                       ? "symmetric multiplicative standard and additive overlapping "
                         "feature-patch correction"
                       : "configured monolithic fallback"},
          {"StandardSolver",
           UsesSingularPatchKspSolver(iodata) ? "BoomerAMG" : "configured"},
          {"FeaturePatchSolver", UsesSingularPatchKspSolver(iodata) ? "SuperLU" : "none"},
          {"FeaturePatchCount", singular_diagnostics.feature_patch_count},
          {"FeaturePatchSumStandardDegreesOfFreedom",
           singular_diagnostics.feature_patch_sum_standard_dofs},
          {"FeaturePatchSumEnrichmentDegreesOfFreedom",
           singular_diagnostics.feature_patch_sum_enrichment_dofs},
          {"FeaturePatchMaximumStandardDegreesOfFreedom",
           singular_diagnostics.feature_patch_maximum_standard_dofs},
          {"FeaturePatchMaximumEnrichmentDegreesOfFreedom",
           singular_diagnostics.feature_patch_maximum_enrichment_dofs},
          {"FeaturePatchMinimumEnrichmentMultiplicity",
           singular_diagnostics.feature_patch_minimum_enrichment_multiplicity},
          {"FeaturePatchMaximumEnrichmentMultiplicity",
           singular_diagnostics.feature_patch_maximum_enrichment_multiplicity}}},
        {"CoefficientOutput",
         {{"File", "singular-coefficients.csv"},
          {"Quantity", "canonical basis coordinate"},
          {"Ordering", "canonical DofKey"},
          {"TrueDofPartitionDependent", true},
          {"PhysicalAmplitude", false}}},
        {"EdgeSlopeOutput",
         {{"File", "singular-edge-slopes.csv"},
          {"Ray", "edge midpoint to opposite-edge midpoint in each tetrahedron"},
          {"SampleCount", edge_slope_options.sample_count},
          {"MinimumBarycentricRadius", edge_slope_options.minimum_barycentric_radius},
          {"MaximumBarycentricRadius", edge_slope_options.maximum_barycentric_radius},
          {"ExpectedSlope", "nu - 1"},
          {"PhysicalAmplitude", false}}},
        {"SurfaceIntegrability",
         triangle_singular
             ? GetSingularSurfaceIntegrabilityMetadata(local_triangle_singular_features)
             : GetSingularSurfaceIntegrabilityMetadata(local_singular_features)}};
    if (triangle_singular)
    {
      const int high_order =
          std::max(fem::singular::H1DuffyReferenceOrder,
                   2 * iodata.solver.singular_elements.quadrature_order + 15);
      auto &quadrature = singular_metadata["Quadrature"];
      quadrature["Rule"] = "feature-aligned and partitioned triangle Duffy quadrature";
      quadrature["DuffyReferenceOrder"] = high_order;
      quadrature["DuffyComparisonOrder"] = high_order - 8;
      quadrature["DuffyRadialPower"] = fem::singular::TriangleDuffyRadialPower;
      quadrature["DuffyMaximumTableEntriesPerRank"] = 0;
      quadrature["DuffyTotalCacheHits"] = 0;
      quadrature["RecursiveSubdivision"] = false;
      quadrature["TotalDuffyPointEvaluations"] = singular_diagnostics.quadrature_leaf_count;
      quadrature.erase("MaximumSubdivisions");
      quadrature.erase("TotalLeafCount");
      quadrature.erase("MaximumDepth");
      singular_metadata["Dimension"] = 2;
      singular_metadata["SingularFeature"] =
          "internal PEC line tip or one-sided finite-metal corner";
      singular_metadata["Postprocessing"] = {
          {"Capacitance", true},
          {"CombinedDiscreteGradient", true},
          {"CombinedFieldEvaluation", true},
          {"DomainEnergy", true},
          {"SmoothRemainderErrorEstimator", true},
          {"ProtectedSingularPatch", true},
          {"SurfaceMeasurements", true},
          {"SurfaceQuadrature", "basis-aware Gauss-Jacobi power expansion"},
          {"CoefficientDiagnostics", true},
          {"TipSlopeDiagnostics", true},
          {"DiscontinuousParaViewSampling", true}};
      singular_metadata.erase("EdgeSlopeOutput");
      singular_metadata["TipSlopeOutput"] = {
          {"File", "singular-tip-slopes.csv"},
          {"Ray", "tip to opposite-edge midpoint in each triangle"},
          {"SampleCount", edge_slope_options.sample_count},
          {"MinimumBarycentricRadius", edge_slope_options.minimum_barycentric_radius},
          {"MaximumBarycentricRadius", edge_slope_options.maximum_barycentric_radius},
          {"ExpectedSlope", "nu - 1"},
          {"PhysicalAmplitude", false}};
    }
    else
    {
      singular_metadata["Dimension"] = 3;
      singular_metadata["SingularFeature"] =
          "internal PEC sheet edge or one-sided finite-metal wedge edge";
      singular_metadata["Postprocessing"] = {
          {"Capacitance", true},
          {"CombinedDiscreteGradient", true},
          {"CombinedFieldEvaluation", true},
          {"DomainEnergy", true},
          {"SmoothRemainderErrorEstimator", true},
          {"ProtectedSingularPatch", true},
          {"SurfaceMeasurements", true},
          {"SurfaceQuadrature",
           "edge-aligned graded Gauss-Jacobi with logarithmic cutoff map"},
          {"CoefficientDiagnostics", true},
          {"EdgeSlopeDiagnostics", true},
          {"DiscontinuousParaViewSampling", true}};
    }
    SaveMetadata("SingularElements", singular_metadata);
    auto field_evaluator =
        triangle_singular ? nullptr : laplace_op.GetSingularFieldEvaluator();
    auto imaginary_field_evaluator =
        triangle_singular ? nullptr : laplace_op.GetSingularFieldEvaluator();
    auto triangle_field_evaluator =
        triangle_singular ? laplace_op.GetTriangleSingularFieldEvaluator() : nullptr;
    auto triangle_imaginary_field_evaluator =
        triangle_singular ? laplace_op.GetTriangleSingularFieldEvaluator() : nullptr;
    std::unique_ptr<TriangleSingularSurfacePostOperator> singular_surface_postoperator;
    std::unique_ptr<TetrahedronSingularSurfacePostOperator>
        tetrahedron_singular_surface_postoperator;
    if (triangle_singular)
    {
      singular_surface_postoperator = std::make_unique<TriangleSingularSurfacePostOperator>(
          iodata.boundaries.postpro, laplace_op.GetMaterialOp(),
          laplace_op.GetH1Space().Get());
      singular_metadata["Postprocessing"]["SurfaceMeasurements"] =
          !singular_surface_postoperator->Empty();
      singular_metadata["SurfaceParticipation"] =
          GetSingularSurfaceParticipationMetadata(iodata);
    }
    else
    {
      tetrahedron_singular_surface_postoperator =
          std::make_unique<TetrahedronSingularSurfacePostOperator>(
              iodata.boundaries.postpro, laplace_op.GetMaterialOp(),
              laplace_op.GetH1Space().Get());
      singular_metadata["Postprocessing"]["SurfaceMeasurements"] =
          !tetrahedron_singular_surface_postoperator->Empty();
      singular_metadata["Postprocessing"]["SurfaceQuadrature"] =
          "edge-aligned graded Gauss-Jacobi with logarithmic cutoff map";
      singular_metadata["SurfaceParticipation"] =
          GetSingularSurfaceParticipationMetadata(iodata);
    }
    std::unique_ptr<SingularParaviewOutput> field_output;
    if (iodata.problem.output_formats.paraview && iodata.solver.electrostatic.n_post > 0)
    {
      field_output = std::make_unique<SingularParaviewOutput>(
          post_dir, laplace_op.GetMesh().Get(),
          std::max(iodata.solver.order, iodata.solver.singular_elements.order),
          iodata.units);
      Mpi::Print(" Singular field output: discontinuous interior Gauss-Legendre "
                 "sampling of order {:d}\n",
                 std::max(iodata.solver.order, iodata.solver.singular_elements.order));
    }
    Mpi::Warning(
        laplace_op.GetComm(),
        "Singular electrostatic postprocessing is currently limited to capacitance "
        "and domain energy extracted from the combined field plus optional combined-field "
        "ParaView sampling. Integrable two-dimensional dielectric surface measurements "
        "use the combined H1 gradient and basis-aware Gauss-Jacobi quadrature. "
        "Three-dimensional dielectric surface measurements use combined-field "
        "edge-aligned singular quadrature. AMR estimates the standard-space smooth "
        "remainder only and masks every enriched element plus one face layer. MFEM "
        "grid-function output remains disabled. Ideal zero-thickness line-tip or "
        "sheet-edge traces with nu <= 1/2 require explicit physical regularization.\n");

    Mpi::Print("\nComputing singular electrostatic fields for {:d} terminal {}\n", n_step,
               (n_step > 1) ? "boundaries" : "boundary");
    int step = 0;
    auto t0 = Timer::Now();
    for (const auto &[idx, data] : laplace_op.GetSources())
    {
      Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1,
                 n_step, idx, Timer::Duration(Timer::Now() - t0).count());
      Mpi::Print("\n");
      laplace_op.GetExcitationVector(idx, *K, V[step], RHS);
      ksp->Mult(RHS, V[step]);

      BlockTimer bt2(Timer::POSTPRO);
      Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
                 linalg::Norml2(laplace_op.GetComm(), V[step]),
                 linalg::Norml2(laplace_op.GetComm(), RHS));

      // This is the exact discrete gradient in the combined standard-plus-enrichment
      // spaces. Physical field evaluation and output require the Phase 8 evaluator.
      E = 0.0;
      Grad.AddMult(V[step], E, -1.0);
      if (triangle_singular)
      {
        triangle_field_evaluator->SetFromTrueDofs(V[step]);
        Vector zero(V[step].Size());
        zero = 0.0;
        triangle_imaginary_field_evaluator->SetFromTrueDofs(zero);
      }
      else
      {
        field_evaluator->SetFromTrueDofs(V[step]);
        Vector zero(V[step].Size());
        zero = 0.0;
        imaginary_field_evaluator->SetFromTrueDofs(zero);
      }
      auto source_coefficients = GatherCoefficientMeasurements(
          laplace_op.GetComm(), idx,
          triangle_singular ? triangle_field_evaluator->GetOwnedCoefficientDiagnostics()
                            : field_evaluator->GetOwnedCoefficientDiagnostics(),
          singular_diagnostics.h1_enrichment_dofs);
      if (root)
      {
        coefficient_measurements.insert(
            coefficient_measurements.end(),
            std::make_move_iterator(source_coefficients.begin()),
            std::make_move_iterator(source_coefficients.end()));
      }
      if (triangle_singular)
      {
        auto source_tip_slopes = GatherTipSlopeMeasurements(
            laplace_op.GetComm(), idx,
            triangle_field_evaluator->FitTipSlopes(local_triangle_singular_features,
                                                   source_vertex_ids, source_element_ids,
                                                   edge_slope_options));
        if (root)
        {
          tip_slope_measurements.insert(tip_slope_measurements.end(),
                                        std::make_move_iterator(source_tip_slopes.begin()),
                                        std::make_move_iterator(source_tip_slopes.end()));
        }
      }
      else
      {
        auto source_edge_slopes = GatherEdgeSlopeMeasurements(
            laplace_op.GetComm(), idx,
            field_evaluator->FitEdgeSlopes(local_singular_features, source_vertex_ids,
                                           source_element_ids, edge_slope_options));
        if (root)
        {
          edge_slope_measurements.insert(
              edge_slope_measurements.end(),
              std::make_move_iterator(source_edge_slopes.begin()),
              std::make_move_iterator(source_edge_slopes.end()));
        }
      }
      if (triangle_singular)
      {
        domain_energy.push_back(MeasureSingularDomainEnergy(
            iodata, laplace_op, *triangle_field_evaluator, V[step], idx));
        if (!singular_surface_postoperator->Empty())
        {
          surface_measurements.push_back(
              {idx, singular_surface_postoperator->MeasureElectrostatic(
                        *triangle_field_evaluator, *triangle_imaginary_field_evaluator,
                        domain_energy.back().total_energy,
                        {iodata.solver.singular_elements.quadrature_order,
                         iodata.solver.singular_elements.abs_tol,
                         iodata.solver.singular_elements.rel_tol,
                         iodata.solver.singular_elements.max_subdivisions})});
        }
      }
      else
      {
        domain_energy.push_back(MeasureSingularDomainEnergy(
            iodata, laplace_op, *field_evaluator, V[step], idx));
        if (!tetrahedron_singular_surface_postoperator->Empty())
        {
          surface_measurements.push_back(
              {idx, tetrahedron_singular_surface_postoperator->MeasureElectrostatic(
                        *field_evaluator, *imaginary_field_evaluator,
                        domain_energy.back().total_energy,
                        {iodata.solver.singular_elements.quadrature_order,
                         iodata.solver.singular_elements.abs_tol,
                         iodata.solver.singular_elements.rel_tol,
                         iodata.solver.singular_elements.max_subdivisions})});
        }
      }
      Vector standard_e;
      standard_e.MakeRef(E, 0, laplace_op.GetNDSpace().GetTrueVSize());
      estimator.AddErrorIndicator(standard_e, domain_energy.back().total_energy, indicator);
      if (field_output && step < iodata.solver.electrostatic.n_post)
      {
        if (triangle_singular)
        {
          field_output->Write(*triangle_field_evaluator, step, idx);
        }
        else
        {
          field_output->Write(*field_evaluator, step, idx);
        }
        Mpi::Print(" Wrote sampled singular fields to disk (ParaView) for source {:d}\n",
                   idx);
      }
      step++;
    }

    BlockTimer bt1(Timer::POSTPRO);
    SaveMetadata(*ksp);
    singular_metadata["EnergyConsistency"] = nlohmann::json::array();
    for (const auto &measurement : domain_energy)
    {
      const double scale = std::abs(measurement.algebraic_integral);
      MFEM_VERIFY(std::isfinite(scale) && scale > 0.0,
                  "Invalid algebraic energy scale for singular metadata!");
      const double discrepancy =
          std::abs(measurement.direct_integral - measurement.algebraic_integral);
      const double bound = measurement.integration_error_bound +
                           measurement.assembly_error_bound +
                           measurement.floating_point_allowance;
      singular_metadata["EnergyConsistency"].push_back(
          {{"Source", measurement.source},
           {"Normalization", "absolute nondimensional integral divided by |V^T K V|"},
           {"DirectToAlgebraicRelativeDiscrepancy", discrepancy / scale},
           {"CertifiedRelativeErrorBound", bound / scale},
           {"IntegrationRelativeErrorBound", measurement.integration_error_bound / scale},
           {"AssemblyRelativeErrorBound", measurement.assembly_error_bound / scale},
           {"FloatingPointRelativeAllowance",
            measurement.floating_point_allowance / scale}});
    }
    SaveMetadata("SingularElements", singular_metadata);
    WriteSingularDomainEnergy(post_dir, iodata, domain_energy, root);
    WriteSingularSurfaceMeasurements(post_dir, surface_measurements, root);
    WriteSingularCoefficientMeasurements(post_dir, iodata, coefficient_measurements, root);
    if (triangle_singular)
    {
      WriteSingularTipSlopeMeasurements(post_dir, iodata, tip_slope_measurements, root);
    }
    else
    {
      WriteSingularEdgeSlopeMeasurements(post_dir, iodata, edge_slope_measurements, root);
    }
    indicator.ZeroElements(GetRefinementProtection(laplace_op.GetMesh().Get()));
    PostprocessSingularTerminals(laplace_op, laplace_op.GetSources(), V);
    return {std::move(indicator), laplace_op.GlobalTrueVSize()};
  }

  // Terminal indices are the set of boundaries over which to compute the capacitance
  // matrix. Terminal boundaries are aliases for ports.
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0, "No terminal boundaries specified for electrostatic simulation!");

  // Right-hand side term and solution vector storage.
  Vector RHS(Grad.Width()), E(Grad.Height());
  std::vector<Vector> V(n_step);

  // Initialize structures for storing and reducing the results of error estimation.
  GradFluxErrorEstimator estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing electrostatic fields for {:d} terminal {}\n", n_step,
             (n_step > 1) ? "boundaries" : "boundary");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero voltage on the specified
    // terminal.
    Mpi::Print("\n");
    laplace_op.GetExcitationVector(idx, *K, V[step], RHS);
    ksp->Mult(RHS, V[step]);

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), V[step]),
               linalg::Norml2(laplace_op.GetComm(), RHS));

    // Compute E = -∇V on the true dofs.
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, V[step], E, idx);

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, total_domain_energy, indicator);

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(*ksp);
  PostprocessTerminals(post_op, laplace_op.GetSources(), V);
  post_op.MeasureFinalize(indicator);
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void ElectrostaticSolver::PostprocessTerminals(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op,
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &V) const
{
  // Postprocess the Maxwell capacitance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the electric field energy based on a unit voltage
  // excitation for each terminal. Alternatively, we could compute the resulting terminal
  // charges from the prescribed voltage to get C directly as:
  //         Q_i = ∫ ρ dV = ∫ ∇ ⋅ (ε E) dV = ∫ (ε E) ⋅ n dS
  // and C_ij = Q_i/V_j. The energy formulation avoids having to locally integrate E = -∇V.
  mfem::DenseMatrix C(V.size());
  for (int i = 0; i < C.Height(); i++)
  {
    // Diagonal: Cᵢᵢ = 2 Uₑ(Vᵢ) / Vᵢ² = (Vᵢᵀ K Vᵢ) / Vᵢ² (with ∀i, Vᵢ = 1)
    auto &V_gf = post_op.GetVGridFunction().Real();
    auto &D_gf = post_op.GetDomainPostOp().D;
    V_gf.SetFromTrueDofs(V[i]);
    post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    C(i, i) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);

    // Off-diagonals: Cᵢⱼ = Uₑ(Vᵢ + Vⱼ) / (Vᵢ Vⱼ) - 1/2 (Vᵢ/Vⱼ Cᵢᵢ + Vⱼ/Vᵢ Cⱼⱼ)
    //                    = (Vⱼᵀ K Vᵢ) / (Vᵢ Vⱼ)
    for (int j = i + 1; j < C.Width(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      C(i, j) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      C(i, j) = C(j, i);
    }
  }
  WriteTerminalMatrices(terminal_sources, C);
}

void ElectrostaticSolver::PostprocessSingularTerminals(
    const LaplaceOperator &laplace_op,
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &V) const
{
  const auto &K = laplace_op.GetUnconstrainedStiffnessMatrix();
  MFEM_VERIFY(K.Height() == K.Width() && V.size() == terminal_sources.size(),
              "Invalid combined electrostatic vectors for capacitance extraction!");

  mfem::DenseMatrix C(V.size());
  Vector KV(K.Height());
  for (int j = 0; j < C.Width(); j++)
  {
    MFEM_VERIFY(V[j].Size() == K.Width(),
                "Invalid combined electrostatic vector for capacitance extraction!");
    K.Mult(V[j], KV);
    for (int i = 0; i <= j; i++)
    {
      C(i, j) = linalg::Dot<Vector>(laplace_op.GetComm(), V[i], KV);
      C(j, i) = C(i, j);
    }
  }
  WriteTerminalMatrices(terminal_sources, C);
}

void ElectrostaticSolver::WriteTerminalMatrices(
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const mfem::DenseMatrix &C) const
{
  MFEM_VERIFY(C.Height() == C.Width() &&
                  C.Height() == static_cast<int>(terminal_sources.size()),
              "Invalid electrostatic capacitance matrix dimensions!");
  mfem::DenseMatrix Cm(C);
  for (int i = 0; i < Cm.Height(); i++)
  {
    for (int j = 0; j < Cm.Width(); j++)
    {
      if (i != j)
      {
        Cm(i, j) = -C(i, j);
        Cm(i, i) -= Cm(i, j);
      }
    }
  }
  mfem::DenseMatrix Cinv(C);
  Cinv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using VT = Units::ValueType;

  // Write capacitance matrix data.
  auto PrintMatrix = [&terminal_sources, this](const std::string &file,
                                               const std::string &name,
                                               const std::string &unit,
                                               const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    int j = 0;
    for (const auto &[idx2, data2] : terminal_sources)
    {
      output.table.insert(fmt::format("i2{}", idx2),
                          fmt::format("{}[i][{}] {}", name, idx2, unit));
      // Use the fact that iterator over i and j is the same span.
      output.table["i"] << idx2;

      auto &col = output.table[fmt::format("i2{}", idx2)];
      for (std::size_t i = 0; i < terminal_sources.size(); i++)
      {
        col << mat(i, j) * scale;
      }
      j++;
    }
    output.WriteFullTableTrunc();
  };
  const double F = iodata.units.Dimensionalize<VT::CAPACITANCE>(1.0);
  PrintMatrix("terminal-C.csv", "C", "(F)", C, F);
  PrintMatrix("terminal-Cinv.csv", "C⁻¹", "(1/F)", Cinv, 1.0 / F);
  PrintMatrix("terminal-Cm.csv", "C_m", "(F)", Cm, F);

  // Also write out a file with terminal voltage excitations.
  {
    TableWithCSVFile terminal_V(post_dir / "terminal-V.csv");
    terminal_V.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_V.table.insert("Vinc", "V_inc[i] (V)");
    for (const auto &[idx, data] : terminal_sources)
    {
      terminal_V.table["i"] << double(idx);
      terminal_V.table["Vinc"] << iodata.units.Dimensionalize<VT::VOLTAGE>(1.0);
    }
    terminal_V.WriteFullTableTrunc();
  }
}

}  // namespace palace
