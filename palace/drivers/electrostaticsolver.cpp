// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/errorindicator.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>

namespace palace
{

namespace
{

constexpr std::uint64_t response_archive_magic = 0x50414c5253503031ULL;
constexpr std::uint32_t response_archive_version = 1;

enum class ArchivedField : std::uint32_t
{
  POTENTIAL = 1,
  FLUX = 2
};

bool EnvironmentFlag(const char *name)
{
  const char *value = std::getenv(name);
  return value && std::string_view(value) == "1";
}

std::optional<std::filesystem::path> ResponseArchiveDirectory()
{
  const char *value = std::getenv("PALACE_RESPONSE_ARCHIVE_DIR");
  if (!value || std::string_view(value).empty())
  {
    return std::nullopt;
  }
  return std::filesystem::absolute(value).lexically_normal();
}

int ResponseArchiveBlockSize()
{
  const char *value = std::getenv("PALACE_RESPONSE_BLOCK_SIZE");
  if (!value)
  {
    return 3;
  }
  const int size = std::stoi(value);
  MFEM_VERIFY(size > 0, "PALACE_RESPONSE_BLOCK_SIZE must be positive!");
  return size;
}

std::filesystem::path ArchivePath(const std::filesystem::path &directory, int source,
                                  int rank, ArchivedField field)
{
  std::ostringstream name;
  name << "source-" << std::setw(6) << std::setfill('0') << source << "-rank-"
       << std::setw(6) << rank << "-" << (field == ArchivedField::POTENTIAL ? "V" : "D")
       << ".bin";
  return directory / name.str();
}

template <typename T>
void WriteArchiveValue(std::ofstream &stream, const T &value)
{
  stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

template <typename T>
void ReadArchiveValue(std::ifstream &stream, T &value)
{
  stream.read(reinterpret_cast<char *>(&value), sizeof(value));
}

void WriteArchivedVector(const std::filesystem::path &directory, int source,
                         ArchivedField field, const Vector &vector, MPI_Comm comm)
{
  const int rank = Mpi::Rank(comm);
  const int size = Mpi::Size(comm);
  const auto path = ArchivePath(directory, source, rank, field);
  std::ofstream stream(path, std::ios::binary | std::ios::trunc);
  MFEM_VERIFY(stream,
              "Unable to create response archive field \"" << path.string() << "\"!");
  const std::int64_t source_value = source;
  const std::int64_t rank_value = rank;
  const std::int64_t size_value = size;
  const std::int64_t local_size = vector.Size();
  WriteArchiveValue(stream, response_archive_magic);
  WriteArchiveValue(stream, response_archive_version);
  WriteArchiveValue(stream, static_cast<std::uint32_t>(field));
  WriteArchiveValue(stream, source_value);
  WriteArchiveValue(stream, rank_value);
  WriteArchiveValue(stream, size_value);
  WriteArchiveValue(stream, local_size);
  stream.write(reinterpret_cast<const char *>(vector.HostRead()),
               local_size * sizeof(double));
  MFEM_VERIFY(stream, "Failed writing response archive field \"" << path.string() << "\"!");
}

Vector ReadArchivedVector(const std::filesystem::path &directory, int source,
                          ArchivedField field, MPI_Comm comm, int expected_size = -1)
{
  const int rank = Mpi::Rank(comm);
  const int size = Mpi::Size(comm);
  const auto path = ArchivePath(directory, source, rank, field);
  std::ifstream stream(path, std::ios::binary);
  MFEM_VERIFY(stream, "Unable to open response archive field \"" << path.string() << "\"!");
  std::uint64_t magic = 0;
  std::uint32_t version = 0, stored_field = 0;
  std::int64_t stored_source = 0, stored_rank = 0, stored_size = 0, local_size = 0;
  ReadArchiveValue(stream, magic);
  ReadArchiveValue(stream, version);
  ReadArchiveValue(stream, stored_field);
  ReadArchiveValue(stream, stored_source);
  ReadArchiveValue(stream, stored_rank);
  ReadArchiveValue(stream, stored_size);
  ReadArchiveValue(stream, local_size);
  MFEM_VERIFY(stream && magic == response_archive_magic &&
                  version == response_archive_version &&
                  stored_field == static_cast<std::uint32_t>(field) &&
                  stored_source == source && stored_rank == rank && stored_size == size &&
                  local_size >= 0 && local_size <= std::numeric_limits<int>::max() &&
                  (expected_size < 0 || local_size == expected_size),
              "Invalid response archive header in \"" << path.string() << "\"!");
  Vector vector(static_cast<int>(local_size));
  stream.read(reinterpret_cast<char *>(vector.HostWrite()), local_size * sizeof(double));
  MFEM_VERIFY(stream && stream.peek() == std::ifstream::traits_type::eof(),
              "Invalid response archive payload in \"" << path.string() << "\"!");
  if (expected_size >= 0)
  {
    for (int i = 0; i < vector.Size(); i++)
    {
      MFEM_VERIFY(std::isfinite(vector.HostRead()[i]),
                  "Nonfinite archived potential in " << path.string());
    }
  }
  vector.UseDevice(true);
  return vector;
}

}  // namespace

void ValidateArchiveEstimateOptions(const IoData &iodata, MPI_Comm comm, bool check_output)
{
  if (!EnvironmentFlag("PALACE_RESPONSE_ESTIMATE_ONLY"))
  {
    return;
  }
  for (const char *flag :
       {"PALACE_RESPONSE_ARCHIVE_ONLY", "PALACE_RESPONSE_REDUCE_ONLY",
        "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "PALACE_RESPONSE_BLOCK_SIZE"})
  {
    MFEM_VERIFY(!std::getenv(flag), "Archive estimation conflicts with " << flag);
  }
  const auto archive = ResponseArchiveDirectory();
  MFEM_VERIFY(archive && std::filesystem::is_directory(*archive),
              "Archive estimation requires an existing PALACE_RESPONSE_ARCHIVE_DIR!");
  MFEM_VERIFY(std::getenv("PALACE_RESPONSE_ESTIMATE_REQUEST"),
              "Archive estimation requires PALACE_RESPONSE_ESTIMATE_REQUEST!");
  MFEM_VERIFY(iodata.problem.type == ProblemType::ELECTROSTATIC &&
                  !iodata.boundaries.prescribed_potential.empty() &&
                  !iodata.solver.electrostatic.response_correction &&
                  iodata.model.refinement.max_it == 0,
              "Archive estimation requires prescribed electrostatics, no response "
              "correction and no AMR!");
  MFEM_VERIFY(!iodata.model.export_prerefined_mesh,
              "Archive estimation rejects ExportPrerefinedMesh!");
  MFEM_VERIFY(iodata.model.partitioning.empty(),
              "Archive estimation does not support Partitioning!");
  for (const auto &[idx, data] : iodata.boundaries.postpro.dielectric)
  {
    MFEM_VERIFY(!data.edge_refinement,
                "Archive estimation does not support geometry-driven edge refinement!");
    MFEM_VERIFY(data.ownership_data_file.empty(),
                "Archive estimation does not support OwnershipDataFile!");
  }
  const auto output = std::filesystem::weakly_canonical(iodata.problem.output);
  auto Disjoint = [&](const std::filesystem::path &input)
  {
    const auto resolved = std::filesystem::weakly_canonical(input);
    auto Contains = [](const auto &parent, const auto &child)
    {
      auto p = parent.begin(), c = child.begin();
      for (; p != parent.end() && c != child.end() && *p == *c; ++p, ++c)
      {
      }
      return p == parent.end();
    };
    MFEM_VERIFY(!Contains(output, resolved) && !Contains(resolved, output),
                "Archive estimate output overlaps input " << input.string());
  };
  Disjoint(*archive);
  Disjoint(iodata.model.mesh);
  Disjoint(std::getenv("PALACE_RESPONSE_ESTIMATE_REQUEST"));
  for (const auto &[idx, data] : iodata.boundaries.prescribed_potential)
  {
    Disjoint(data.data_file);
  }
  if (check_output)
  {
    MFEM_VERIFY(!std::filesystem::exists(output),
                "Archive estimate output directory must not already exist!");
    // Every rank must finish the read-only preflight before root creates the output.
    Mpi::Barrier(comm);
  }
}

ErrorIndicator
ElectrostaticSolver::EstimateArchivedFields(LaplaceOperator &laplace_op, const Operator &K,
                                            const std::filesystem::path &archive) const
{
  using json = nlohmann::json;
  using VT = Units::ValueType;
  const auto comm = laplace_op.GetComm();
  const auto &Grad = laplace_op.GetGradMatrix();
  std::ifstream request_stream(std::getenv("PALACE_RESPONSE_ESTIMATE_REQUEST"));
  MFEM_VERIFY(request_stream, "Cannot read archive estimate request!");
  const auto request = json::parse(request_stream);
  MFEM_VERIFY(request.at("Version") == 1, "Unsupported archive estimate request version!");
  const auto source_ids = request.at("SourceIds").get<std::vector<int>>();
  std::vector<int> configured_ids;
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    configured_ids.push_back(idx);
  }
  MFEM_VERIFY(source_ids == configured_ids,
              "Estimate request SourceIds must exactly match the ordered configuration!");
  const auto zero_ids = request.at("ZeroTraceIndices").get<std::vector<int>>();
  for (int idx : zero_ids)
  {
    MFEM_VERIFY(std::find(source_ids.begin(), source_ids.end(), idx) != source_ids.end(),
                "Unknown constrained source!");
  }
  const auto &cases = request.at("Excitations");
  MFEM_VERIFY(cases.is_array() && !cases.empty(), "No archive excitations requested!");
  const auto &validation = request.at("Validation");
  const double rtol = validation.at("RelResidualTol");
  const double atol = validation.at("AbsResidualTol");
  const double bc_volts = validation.at("BCAbsTolV");
  const double energy_floor = validation.at("MinEnergyJ");
  const double cancellation_floor = validation.at("MinCancellationRatio");
  for (double value : {rtol, atol, bc_volts, energy_floor, cancellation_floor})
  {
    MFEM_VERIFY(std::isfinite(value) && value >= 0.0, "Invalid diagnostic tolerance!");
  }
  MFEM_VERIFY(energy_floor > 0.0 && cancellation_floor > 0.0 && cancellation_floor < 1.0,
              "MinEnergyJ must be positive and MinCancellationRatio must be in (0,1)!");
  const double voltage_scale = iodata.units.GetScaleFactor<VT::VOLTAGE>();
  const double energy_scale = iodata.units.GetScaleFactor<VT::ENERGY>();
  const double bc_tol = bc_volts / voltage_scale;
  MFEM_VERIFY(std::isfinite(voltage_scale) && voltage_scale > 0.0 &&
                  std::isfinite(energy_scale) && energy_scale > 0.0 &&
                  std::isfinite(bc_tol),
              "Invalid diagnostic unit scales!");
  const bool localize = request.at("WriteElementIndicators").get<bool>();
  std::vector<int> quadrature_extras;
  if (request.contains("SurfaceQuadratureExtras"))
  {
    const auto &extras = request.at("SurfaceQuadratureExtras");
    MFEM_VERIFY(extras.is_array() && !extras.empty() && extras.front() == 0,
                "SurfaceQuadratureExtras must be an array starting at zero!");
    for (const auto &value : extras)
    {
      MFEM_VERIFY(value.is_number_integer() && value >= 0 && value <= 12,
                  "SurfaceQuadratureExtras must contain integer orders in [0,12]!");
      const int extra = value.get<int>();
      MFEM_VERIFY(quadrature_extras.empty() || extra > quadrature_extras.back(),
                  "SurfaceQuadratureExtras must be sorted and unique!");
      quadrature_extras.push_back(extra);
    }
  }
  const auto &linear = iodata.solver.linear;
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op, nullptr,
                                                   &surface_post_geometry);
  GradFluxErrorEstimator<Vector> estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg);
  json report = {
      {"Version", 1},
      {"Status", "incomplete"},
      {"Request", request},
      {"SampleCount", 0},
      {"PDESolves", 0},
      {"Ranks", Mpi::Size(comm)},
      {"Order", iodata.solver.order},
      {"GlobalH1TrueDofs", laplace_op.GlobalTrueVSize()},
      {"GlobalNDTrueDofs", laplace_op.GetNDSpace().GlobalTrueVSize()},
      {"GlobalRTTrueDofs", laplace_op.GetRTSpace().GlobalTrueVSize()},
      {"VoltageScaleV", voltage_scale},
      {"EnergyScaleJ", energy_scale},
      {"BCAbsTolInternal", bc_tol},
      {"ResidualNorm", "Euclidean true-dof norm of eliminated K V - RHS"},
      {"BoundaryMismatch",
       "essential true-dof max versus imposed discrete BC; not ideal P1 error"},
      {"IndicatorDefinition",
       "eta_raw^2 = integral |sqrt(eps) E - invsqrt(eps) D|^2; eta = eta_raw/sqrt(2 U)"},
      {"Interpretation",
       "volume energy-norm heuristic, not a boundary or interface error bound"},
      {"Recovery",
       {{"RelTol", linear.estimator_tol},
        {"AbsTol", std::numeric_limits<double>::epsilon()},
        {"MaxIts", linear.estimator_max_it},
        {"Multigrid", linear.estimator_mg},
        {"ResidualNorm", "preconditioned CG norm"}}},
      {"Checks", json::array()},
      {"Excitations", json::array()}};
  auto Save = [&]()
  {
    for (const auto &value : report.flatten())
    {
      MFEM_VERIFY(!value.is_number_float() || std::isfinite(value.get<double>()),
                  "Nonfinite diagnostic report value!");
    }
    if (root)
    {
      std::ofstream stream(post_dir / "archive-estimates.json");
      stream << report.dump(2) << '\n';
      MFEM_VERIFY(stream, "Cannot write archive diagnostic report!");
    }
  };
  // One small rank-layout record per rank, never a replicated mesh.
  {
    // A rank-local diagnostic fingerprint, not an independent historical archive hash.
    // FNV-1a hashes native-endian geometry bytes and the ordered FE dof maps. The launcher
    // separately records cryptographic hashes of the mesh, executable and these records.
    auto Fingerprint = [&](const FiniteElementSpace &space)
    {
      std::uint64_t hash = 14695981039346656037ULL;
      auto Add = [&](const auto &value)
      {
        const auto *bytes = reinterpret_cast<const unsigned char *>(&value);
        for (std::size_t i = 0; i < sizeof(value); i++)
        {
          hash = (hash ^ bytes[i]) * 1099511628211ULL;
        }
      };
      const auto &mesh = space.GetMesh().Get();
      Add(Mpi::Rank(comm));
      Add(Mpi::Size(comm));
      Add(space.GetVSize());
      Add(space.GetTrueVSize());
      Add(space.GetMaxElementOrder());
      for (int i = 0; i < mesh.GetNV(); i++)
      {
        for (int j = 0; j < mesh.SpaceDimension(); j++)
        {
          Add(mesh.GetVertex(i)[j]);
        }
      }
      if (const auto *nodes = mesh.GetNodes())
      {
        for (int i = 0; i < nodes->Size(); i++)
        {
          Add(nodes->HostRead()[i]);
        }
      }
      mfem::Array<int> dofs, vertices;
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        Add(mesh.GetAttribute(i));
        Add(mesh.GetElementBaseGeometry(i));
        mesh.GetElementVertices(i, vertices);
        for (int v : vertices)
        {
          Add(v);
        }
        space.Get().GetElementDofs(i, dofs);
        for (int dof : dofs)
        {
          Add(dof);
        }
      }
      for (int i = 0; i < space.GetVSize(); i++)
      {
        Add(space.Get().GetLocalTDofNumber(i));
      }
      return fmt::format("{:016x}", hash);
    };
    std::ofstream layout(post_dir /
                         fmt::format("archive-layout-rank-{:06d}.json", Mpi::Rank(comm)));
    layout << json({{"Rank", Mpi::Rank(comm)},
                    {"H1TrueDofs", Grad.Width()},
                    {"H1LayoutFNV1a64", Fingerprint(laplace_op.GetH1Space())},
                    {"NDLayoutFNV1a64", Fingerprint(laplace_op.GetNDSpace())},
                    {"RTLayoutFNV1a64", Fingerprint(laplace_op.GetRTSpace())},
                    {"NDTrueDofs", Grad.Height()},
                    {"RTTrueDofs", laplace_op.GetRTSpace().GetTrueVSize()},
                    {"LocalElements", laplace_op.GetMesh().GetNE()}})
                  .dump(2)
           << '\n';
    MFEM_VERIFY(layout, "Cannot write archive rank layout!");
  }
  auto Finite = [&](const Vector &v)
  {
    int finite = 1;
    const auto *values = v.HostRead();
    for (int i = 0; i < v.Size(); i++)
    {
      finite = finite && std::isfinite(values[i]);
    }
    Mpi::GlobalMin(1, &finite, comm);
    MFEM_VERIFY(finite, "Nonfinite diagnostic vector!");
  };
  // Classify owned true DOFs once. Ground wins at trace/ground intersections, exactly
  // as in the existing prescribed-potential projection; no projection is changed here.
  auto &h1 = laplace_op.GetH1Space().Get();
  const int max_attr = h1.GetParMesh()->bdr_attributes.Max();
  mfem::Array<int> matching_marker(max_attr), ground_marker(max_attr);
  matching_marker = 0;
  ground_marker = 0;
  for (const auto &[idx, data] : iodata.boundaries.prescribed_potential)
  {
    for (int attr : data.attributes)
    {
      if (attr > 0 && attr <= max_attr)
      {
        matching_marker[attr - 1] = 1;
      }
    }
  }
  for (int attr : iodata.boundaries.pec.attributes)
  {
    if (attr > 0 && attr <= max_attr)
    {
      ground_marker[attr - 1] = 1;
    }
  }
  mfem::Array<int> matching_dofs, ground_dofs;
  h1.GetEssentialTrueDofs(matching_marker, matching_dofs);
  h1.GetEssentialTrueDofs(ground_marker, ground_dofs);
  std::vector<int> category(Grad.Width(), 0);
  for (int dof : matching_dofs)
  {
    category[dof] = 1;
  }
  for (int dof : ground_dofs)
  {
    category[dof] = category[dof] == 1 ? 3 : 2;
  }
  long long boundary_counts[3] = {0, 0, 0};
  for (int tag : category)
  {
    if (tag > 0)
    {
      boundary_counts[tag - 1]++;
    }
  }
  Mpi::GlobalSum(3, boundary_counts, comm);
  report["BoundaryCategories"] = {
      {"MatchingOnlyTrueDofs", boundary_counts[0]},
      {"PhysicalGroundTrueDofs", boundary_counts[1] + boundary_counts[2]},
      {"IntersectionTrueDofs", boundary_counts[2]},
      {"Precedence",
       "physical ground includes intersections; matching-only excludes them"}};
  Vector residual(Grad.Width());
  auto Check =
      [&](const Vector &v, const Vector &bc, const Vector &rhs, const std::string &label)
  {
    Finite(v);
    Finite(bc);
    Finite(rhs);
    K.Mult(v, residual);
    residual -= rhs;
    Finite(residual);
    double mismatch = 0.0;
    double boundary_max[3] = {0.0, 0.0, 0.0};
    for (int dof : laplace_op.GetDbcTDofList())
    {
      const double error = std::abs(v.HostRead()[dof] - bc.HostRead()[dof]);
      mismatch = std::max(mismatch, error);
      const int tag = category[dof];
      if (tag > 0)
      {
        boundary_max[tag - 1] = std::max(boundary_max[tag - 1], error);
      }
    }
    Mpi::GlobalMax(1, &mismatch, comm);
    Mpi::GlobalMax(3, boundary_max, comm);
    const double norm = linalg::Norml2(comm, residual);
    const double rhs_norm = linalg::Norml2(comm, rhs);
    const double threshold = atol + rtol * rhs_norm;
    const bool pass = std::isfinite(norm) && std::isfinite(rhs_norm) &&
                      std::isfinite(threshold) && norm <= threshold && mismatch <= bc_tol;
    report["Checks"].push_back(
        {{"Label", label},
         {"ResidualNorm", norm},
         {"RHSNorm", rhs_norm},
         {"ResidualThreshold", threshold},
         {"RelativeResidual", rhs_norm > 0.0 ? json(norm / rhs_norm) : json(nullptr)},
         {"BCMaxInternal", mismatch},
         {"BCMaxV", mismatch * voltage_scale},
         {"MatchingOnlyBCMaxV",
          boundary_counts[0] > 0 ? json(boundary_max[0] * voltage_scale) : json(nullptr)},
         {"PhysicalGroundBCMaxV",
          boundary_counts[1] + boundary_counts[2] > 0
              ? json(std::max(boundary_max[1], boundary_max[2]) * voltage_scale)
              : json(nullptr)},
         {"IntersectionBCMaxV",
          boundary_counts[2] > 0 ? json(boundary_max[2] * voltage_scale) : json(nullptr)},
         {"NearZeroRHS", rhs_norm <= atol},
         {"Passed", pass}});
    Save();
    MFEM_VERIFY(pass, "Archive field validation failed for "
                          << label << ": residual=" << norm << " (limit=" << threshold
                          << "), BC mismatch=" << mismatch * voltage_scale << " V");
  };
  ErrorIndicator aggregate;
  Vector v(Grad.Width()), bc(Grad.Width()), rhs(Grad.Width());
  Vector source_bc, source_rhs, e(Grad.Height()), d(laplace_op.GetRTSpace().GetTrueVSize());
  std::set<std::string> names;
  int case_number = 0;
  for (const auto &excitation : cases)
  {
    const std::string name = excitation.at("Name");
    MFEM_VERIFY(!name.empty() && names.insert(name).second,
                "Duplicate/empty excitation name!");
    const auto coefficients = excitation.at("Coefficients").get<std::vector<double>>();
    MFEM_VERIFY(coefficients.size() == source_ids.size(), "Incomplete coefficient table!");
    v = bc = rhs = 0.0;
    double constituent_norm_sum = 0.0;
    for (std::size_t i = 0; i < source_ids.size(); i++)
    {
      const double c = coefficients[i];
      const int source = source_ids[i];
      MFEM_VERIFY(std::isfinite(c), "Nonfinite excitation coefficient!");
      MFEM_VERIFY(c == 0.0 ||
                      std::find(zero_ids.begin(), zero_ids.end(), source) == zero_ids.end(),
                  "Nonzero coefficient on constrained source " << source);
      if (c == 0.0)
      {
        continue;  // No magnitude cutoff: every nonzero coefficient is retained.
      }
      auto field =
          ReadArchivedVector(archive, source, ArchivedField::POTENTIAL, comm, Grad.Width());
      laplace_op.GetExcitationVector(source, K, source_bc, source_rhs);
      Check(field, source_bc, source_rhs, fmt::format("{}:source-{}", name, source));
      constituent_norm_sum += std::abs(c) * linalg::Norml2(comm, field);
      v.Add(c, field);
      bc.Add(c, source_bc);
      rhs.Add(c, source_rhs);
    }
    Check(v, bc, rhs, name);
    e = 0.0;
    Grad.AddMult(v, e, -1.0);
    Finite(e);
    const auto &projector = estimator.GetFluxProjector();
    const int iterations_before = projector.GetTotalIterations();
    estimator.RecoverFlux(e, d);
    Finite(d);
    const auto energies = post_op.GetElectrostaticEnergies(v, e, &d);
    const double energy_j = energies.domain * energy_scale;
    MFEM_VERIFY(std::isfinite(energy_j) && energy_j >= 0.0, "Invalid domain energy!");
    // Et=0 requests the unnormalized integral from the existing estimator, including for
    // a genuinely zero solution. No division by tiny energy is performed here.
    ErrorIndicator raw;
    estimator.AddErrorIndicator(e, d, 0.0, raw);
    Finite(raw.Local());
    const double raw_norm = raw.Norml2(comm);
    const double raw_j = raw_norm * std::sqrt(energy_scale);
    const bool normalized = energy_j > energy_floor;
    const double solution_norm = linalg::Norml2(comm, v);
    const double cancellation =
        constituent_norm_sum > 0.0 ? solution_norm / constituent_norm_sum : 0.0;
    const double recovery_residual = projector.GetFinalRelativeResidual();
    MFEM_VERIFY(std::isfinite(raw_j) && std::isfinite(cancellation) &&
                    std::isfinite(recovery_residual),
                "Nonfinite diagnostic output!");
    json result = {
        {"Name", name},
        {"SampleCount", 1},
        {"Status", projector.GetConverged() ? "measured" : "recovery_unconverged"},
        {"EnergyJ", energy_j},
        {"EnergyInternal", energies.domain},
        {"EtaRawSqrtJ", raw_j},
        {"NormalizationSqrtJ", std::sqrt(2.0 * energy_j)},
        {"Eta",
         normalized ? json(raw_norm / std::sqrt(2.0 * energies.domain)) : json(nullptr)},
        {"NormalizationStatus",
         normalized ? "defined_above_requested_floor" : "zero_or_small_energy"},
        {"CoefficientCancellationRatio", cancellation},
        {"CancellationStatus",
         constituent_norm_sum == 0.0
             ? "zero_input"
             : (cancellation < cancellation_floor ? "strong_cancellation"
                                                  : "above_requested_floor")},
        {"NormalizationUsable",
         normalized && cancellation >= cancellation_floor && projector.GetConverged()},
        {"RecoveryInitialResidual", projector.GetInitialResidual()},
        {"RecoveryFinalResidual", projector.GetFinalResidual()},
        {"RecoveryConverged", projector.GetConverged()},
        {"RecoveryIterations", projector.GetTotalIterations() - iterations_before},
        {"RecoveryRelativeResidual", recovery_residual},
        {"Interfaces", json::array()},
        {"InterfaceResponses", json::array()},
        {"SurfaceUsesRecoveredFlux", post_op.NeedsRecoveredElectricFlux()}};
    // Reuse exactly the response-matrix integrator, with a single field (no scalar
    // indicator superposition). Supply recovered flux only when the observable uses it.
    const std::vector<Vector> electric_fields{e};
    const std::vector<Vector> recovered_fields = post_op.NeedsRecoveredElectricFlux()
                                                     ? std::vector<Vector>{d}
                                                     : std::vector<Vector>{};
    auto ResponseJSON = [&](const auto &matrices)
    {
      json entries_json = json::array();
      for (const auto &[idx, entries] : matrices)
      {
        for (const auto &entry : entries)
        {
          for (double value :
               {entry.energy_inside(0, 0), entry.energy_inside_normal(0, 0),
                entry.energy_inside_tangential(0, 0), entry.energy_total(0, 0),
                entry.energy_total_normal(0, 0), entry.energy_total_tangential(0, 0)})
          {
            MFEM_VERIFY(std::isfinite(value * energy_scale) && value >= 0.0,
                        "Invalid polarization-resolved diagonal interface energy!");
          }
          entries_json.push_back(
              {{"Index", idx},
               {"DistanceM", iodata.units.Dimensionalize<VT::LENGTH>(entry.distance)},
               {"InsideJ", entry.energy_inside(0, 0) * energy_scale},
               {"InsideNormalJ", entry.energy_inside_normal(0, 0) * energy_scale},
               {"InsideTangentialJ", entry.energy_inside_tangential(0, 0) * energy_scale},
               {"TotalJ", entry.energy_total(0, 0) * energy_scale},
               {"TotalNormalJ", entry.energy_total_normal(0, 0) * energy_scale},
               {"TotalTangentialJ", entry.energy_total_tangential(0, 0) * energy_scale}});
        }
      }
      return entries_json;
    };
    const auto responses = post_op.GetInterfaceElectricFieldEnergyMatrices(
        electric_fields, recovered_fields.empty() ? nullptr : &recovered_fields);
    result["InterfaceResponses"] = ResponseJSON(responses);
    if (!quadrature_extras.empty())
    {
      using Q = fem::DefaultIntegrationOrder;
      const std::array<int, 4> original{Q::p_trial, Q::q_order_jac, Q::q_order_extra_pk,
                                        Q::q_order_extra_qk};
      result["SurfaceQuadratureControls"] = {
          {"TrialOrder", original[0]},
          {"JacobianOrderIncluded", bool(original[1])},
          {"BaseExtraPk", original[2]},
          {"BaseExtraQk", original[3]},
          {"FieldStatus", "same fixed E/D as primary evaluation"}};
      result["SurfaceQuadratureSweeps"] = json::array();
      json local_rules = json::array();
      for (int extra : quadrature_extras)
      {
        std::vector<SurfacePostOperator::InterfaceQuadratureRule> rules;
        const auto matrices = post_op.GetInterfaceElectricFieldEnergyMatrices(
            electric_fields, recovered_fields.empty() ? nullptr : &recovered_fields, extra,
            &rules);
        result["SurfaceQuadratureSweeps"].push_back(
            {{"ExtraOrder", extra},
             {"Status", "measured"},
             {"InterfaceResponses", ResponseJSON(matrices)}});
        for (const auto &r : rules)
        {
          local_rules.push_back({{"ExtraOrder", extra},
                                 {"Interface", r.interface_index},
                                 {"Geometry", r.geometry},
                                 {"RequestedOrder", r.requested_order},
                                 {"ActualRuleOrder", r.actual_order},
                                 {"PointCount", r.points},
                                 {"LocalFaces", r.local_faces},
                                 {"MinimumReferenceWeight", r.minimum_weight},
                                 {"OwnershipRule", r.ownership_rule}});
        }
      }
      MFEM_VERIFY(
          (original == std::array<int, 4>{Q::p_trial, Q::q_order_jac, Q::q_order_extra_pk,
                                          Q::q_order_extra_qk}),
          "Surface-only sweep modified global quadrature defaults!");
      std::ofstream rules_file(
          post_dir / fmt::format("archive-quadrature-case-{:04d}-rank-{:06d}.json",
                                 case_number, Mpi::Rank(comm)));
      rules_file << local_rules.dump(2) << '\n';
      MFEM_VERIFY(rules_file, "Cannot write surface quadrature rule evidence!");
    }
    for (const auto &[idx, surface] : energies.interfaces)
    {
      const double value = surface.energy * energy_scale;
      MFEM_VERIFY(std::isfinite(value), "Nonfinite interface energy!");
      json entry = {{"Index", idx}, {"EnergyJ", value}, {"Windows", json::array()}};
      for (const auto &window : surface.edge_energies)
      {
        const double inside = (surface.energy - window.energy_outside) * energy_scale;
        const double outside = window.energy_outside * energy_scale;
        MFEM_VERIFY(std::isfinite(inside) && std::isfinite(outside),
                    "Nonfinite window energy!");
        entry["Windows"].push_back(
            {{"DistanceM", iodata.units.Dimensionalize<VT::LENGTH>(window.distance)},
             {"InsideJ", inside},
             {"OutsideJ", outside}});
      }
      result["Interfaces"].push_back(entry);
    }
    if (localize)
    {
      const auto filename = fmt::format("archive-elements-case-{:04d}-rank-{:06d}.csv",
                                        case_number, Mpi::Rank(comm));
      std::ofstream elements(post_dir / filename);
      elements
          << "local_element,attribute,center_x_m,center_y_m,center_z_m,eta_raw_squared_J\n"
          << std::setprecision(17);
      auto &local_mesh = laplace_op.GetH1Space().GetMesh().Get();
      mfem::Vector center(local_mesh.SpaceDimension());
      for (int i = 0; i < raw.Local().Size(); i++)
      {
        local_mesh.GetElementCenter(i, center);
        const double eta = raw.Local().HostRead()[i];
        const double eta_squared_j = eta * eta * energy_scale;
        MFEM_VERIFY(std::isfinite(eta_squared_j), "Nonfinite localized indicator!");
        elements << i << ',' << local_mesh.GetAttribute(i);
        for (int j = 0; j < 3; j++)
        {
          const double coordinate =
              j < center.Size() ? iodata.units.Dimensionalize<VT::LENGTH>(center[j]) : 0.0;
          MFEM_VERIFY(std::isfinite(coordinate), "Nonfinite element center!");
          elements << ',' << coordinate;
        }
        elements << ',' << eta_squared_j << '\n';
      }
      MFEM_VERIFY(elements, "Cannot write localized archive indicators!");
    }
    aggregate.AddIndicator(raw.Local());
    report["Excitations"].push_back(result);
    report["SampleCount"] = ++case_number;
    Save();
  }
  report["Status"] = "complete";
  Save();
  Mpi::Print("\nArchive diagnostics: {:d} measured excitations; no PDE solves. "
             "Returned aggregate is unnormalized; use per-excitation JSON, not the AMR "
             "summary, for normalized indicators.\n",
             case_number);
  return aggregate;
}

std::pair<ErrorIndicator, long long int>
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  LaplaceOperator laplace_op(iodata, mesh);
  // The archive reduction solves nothing: it needs the gradient, the electric energy mass
  // operator and the grid functions, not the stiffness matrix nor the preconditioner
  // setup (decision 62(4)).
  const bool archive_reduce_only = EnvironmentFlag("PALACE_RESPONSE_REDUCE_ONLY");
  decltype(laplace_op.GetStiffnessMatrix()) K;
  if (!archive_reduce_only)
  {
    K = laplace_op.GetStiffnessMatrix();
  }
  if (EnvironmentFlag("PALACE_RESPONSE_ESTIMATE_ONLY"))
  {
    MFEM_VERIFY(!archive_reduce_only,
                "PALACE_RESPONSE_ESTIMATE_ONLY and PALACE_RESPONSE_REDUCE_ONLY exclude each "
                "other!");
    ValidateArchiveEstimateOptions(iodata, laplace_op.GetComm(), false);
    SaveMetadata(laplace_op.GetH1Spaces());
    auto indicator = EstimateArchivedFields(laplace_op, *K, *ResponseArchiveDirectory());
    SaveLinearSolverMetadata(laplace_op.GetComm(), 0, 0);
    return {std::move(indicator), laplace_op.GlobalTrueVSize()};
  }
  const auto *response_config = iodata.solver.electrostatic.response_correction
                                    ? &*iodata.solver.electrostatic.response_correction
                                    : nullptr;
  const bool postprocess_response =
      response_config && response_config->IncludesPostprocessing();
  const bool self_consistent_response =
      response_config && response_config->IncludesSelfConsistent();
  MFEM_VERIFY(!archive_reduce_only || !response_config,
              "PALACE_RESPONSE_REDUCE_ONLY reduces archived fields of a configuration "
              "without ResponseCorrection!");
  std::unique_ptr<SurfaceResponseOperator> response_correction;
  std::unique_ptr<SumOperator> corrected_K;
  const Operator *system_K = K.get();
  if (response_config)
  {
    response_correction =
        std::make_unique<SurfaceResponseOperator>(iodata, laplace_op, &response_geometry);
    if (root && !response_correction->GetPatchAssignments().empty())
    {
      TableWithCSVFile patches(post_dir / "surface-response-patches.csv");
      patches.table.insert(Column("model", "model", 0, 0, 2, ""));
      for (const char *coordinate : {"x", "y", "z"})
      {
        patches.table.insert(fmt::format("origin_{}", coordinate),
                             fmt::format("origin {} (m)", coordinate));
      }
      for (const char *axis : {"u", "v", "w"})
      {
        for (const char *coordinate : {"x", "y", "z"})
        {
          patches.table.insert(fmt::format("axis_{}_{}", axis, coordinate),
                               fmt::format("axis {} {}", axis, coordinate));
        }
      }
      patches.table.insert("weight", "weight");
      for (const auto &patch : response_correction->GetPatchAssignments())
      {
        patches.table["model"] << patch.model;
        for (int d = 0; d < 3; d++)
        {
          patches.table[fmt::format("origin_{}", "xyz"[d])]
              << iodata.units.Dimensionalize<Units::ValueType::LENGTH>(patch.origin[d]);
          patches.table[fmt::format("axis_u_{}", "xyz"[d])] << patch.axis_u[d];
          patches.table[fmt::format("axis_v_{}", "xyz"[d])] << patch.axis_v[d];
          patches.table[fmt::format("axis_w_{}", "xyz"[d])] << patch.axis_w[d];
        }
        patches.table["weight"] << patch.weight;
      }
      patches.WriteFullTableTrunc();
    }
    if (self_consistent_response)
    {
      corrected_K = std::make_unique<SumOperator>(*K, *response_correction);
      system_K = corrected_K.get();
    }
  }
  const auto &Grad = laplace_op.GetGradMatrix();
  SaveMetadata(laplace_op.GetH1Spaces());

  // Preserve the historical thin-metal solve and outputs even when response correction is
  // enabled. The corrected solve uses the same assembled thin operator as its
  // preconditioner.
  KspSolver ksp(iodata, laplace_op.GetH1Spaces());
  if (!archive_reduce_only)
  {
    ksp.SetOperators(*K, *K);
  }

  // Source indices are either equipotential terminals or prescribed potential traces.
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op, nullptr,
                                                   &surface_post_geometry);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0,
              "No terminal or prescribed potential boundaries specified for electrostatic "
              "simulation!");
  const auto response_archive = ResponseArchiveDirectory();
  const bool archive_stream_only = EnvironmentFlag("PALACE_RESPONSE_ARCHIVE_ONLY");
  const bool response_source_timing = EnvironmentFlag("PALACE_RESPONSE_SOURCE_TIMING");
  const bool archive_recycle_initial_guess =
      archive_stream_only && EnvironmentFlag("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS");
  MFEM_VERIFY(!archive_stream_only ||
                  (response_archive && !iodata.boundaries.prescribed_potential.empty() &&
                   iodata.solver.electrostatic.response_matrix &&
                   iodata.solver.electrostatic.aggregate_response_matrix),
              "PALACE_RESPONSE_ARCHIVE_ONLY requires PALACE_RESPONSE_ARCHIVE_DIR and an "
              "aggregated prescribed-potential response-matrix configuration!");
  MFEM_VERIFY(!archive_reduce_only ||
                  (response_archive && !iodata.boundaries.prescribed_potential.empty() &&
                   iodata.solver.electrostatic.response_matrix &&
                   iodata.solver.electrostatic.aggregate_response_matrix),
              "PALACE_RESPONSE_REDUCE_ONLY requires PALACE_RESPONSE_ARCHIVE_DIR and an "
              "aggregated prescribed-potential response-matrix configuration!");
  if (response_archive)
  {
    if (root)
    {
      std::filesystem::create_directories(*response_archive);
    }
    Mpi::Barrier(laplace_op.GetComm());
  }

  // Right-hand side term and solution vector storage.
  Vector RHS(Grad.Width()), E(Grad.Height()), D(laplace_op.GetRTSpace().GetTrueVSize());
  D.UseDevice(true);
  std::vector<Vector> V(n_step);
  std::vector<Vector> V_corrected(self_consistent_response ? n_step : 0);
  std::vector<Vector> D_basis(post_op.NeedsRecoveredElectricFlux() &&
                                      iodata.solver.electrostatic.response_matrix &&
                                      !archive_stream_only
                                  ? n_step
                                  : 0);
  using EnergyData = PostOperator<ProblemType::ELECTROSTATIC>::ElectrostaticEnergyData;
  struct CorrectedResult
  {
    int source;
    EnergyData raw;
    EnergyData postprocessed_fixed_trace;
    EnergyData postprocessed_fixed_flux;
    EnergyData corrected;
    std::vector<SurfaceResponseOperator::ModelContribution> raw_model_contributions;
    std::vector<SurfaceResponseOperator::ModelContribution> corrected_model_contributions;
    std::map<int, double> trace_closure_spread;
    double maximum_trace_closure_spread;
    double response_weighted_trace_closure_spread;
    double trace_closure_response_failure_fraction;
    bool has_postprocessed;
    bool has_self_consistent;
    bool confident;
  };
  std::vector<CorrectedResult> corrected_results;
  corrected_results.reserve(response_correction ? n_step : 0);
  std::vector<std::pair<int, std::vector<SurfaceResponseOperator::PatchTrace>>>
      spatial_patch_traces;
  long long int raw_linear_solves = 0;
  long long int raw_linear_iterations = 0;
  long long int corrected_linear_solves = 0;
  long long int corrected_linear_iterations = 0;
  Vector previous_archive_solution;

  // Archive reduction never recovers a field. A streaming worker needs this operator
  // only for requested flux recovery, not AMR estimation. Avoid building the otherwise
  // unused RT hierarchy, which can exceed the memory of the electrostatic solve itself.
  std::unique_ptr<GradFluxErrorEstimator<Vector>> estimator;
  if (!archive_reduce_only &&
      (!archive_stream_only || post_op.NeedsRecoveredElectricFlux()))
  {
    estimator = std::make_unique<GradFluxErrorEstimator<Vector>>(
        laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  ErrorIndicator indicator;

  if (archive_reduce_only)
  {
    Mpi::Print("\nReducing archived electrostatic response fields from {}\n",
               response_archive->string());
    PostprocessArchivedResponseMatrix(post_op, laplace_op, Grad, *response_archive,
                                      ResponseArchiveBlockSize());
    SaveLinearSolverMetadata(laplace_op.GetComm(), 0, 0);
    return {indicator, laplace_op.GlobalTrueVSize()};
  }

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing electrostatic fields for {:d} {}\n", n_step,
             (n_step > 1) ? "sources" : "source");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    const auto source_start = Timer::Now();
    double source_solve_seconds = 0.0;
    long long source_iterations = 0;
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero voltage on the specified
    // terminal.
    Mpi::Print("\n");
    laplace_op.GetExcitationVector(idx, *K, V[step], RHS);
    if (archive_recycle_initial_guess && previous_archive_solution.Size() == V[step].Size())
    {
      Vector essential_values;
      V[step].GetSubVector(laplace_op.GetDbcTDofList(), essential_values);
      V[step] = previous_archive_solution;
      V[step].SetSubVector(laplace_op.GetDbcTDofList(), essential_values);
      ksp.SetInitialGuess(true);
    }
    const double rhs_norm = linalg::Norml2(laplace_op.GetComm(), RHS);
    const bool zero_response = iodata.solver.electrostatic.response_matrix &&
                               rhs_norm <= 100.0 * std::numeric_limits<double>::epsilon();
    if (zero_response)
    {
      Mpi::Print(" Prescribed trace has no active boundary degrees of freedom; "
                 "storing a zero response field\n");
      V[step] = 0.0;
    }
    Vector corrected_rhs;
    if (self_consistent_response)
    {
      V_corrected[step] = V[step];
      corrected_rhs = RHS;
      response_correction->EliminateRHS(V_corrected[step], corrected_rhs);
    }
    if (!zero_response)
    {
      const auto solves_before = ksp.NumTotalMult();
      const auto iterations_before = ksp.NumTotalMultIterations();
      const auto solve_start = Timer::Now();
      ksp.Mult(RHS, V[step]);
      source_solve_seconds = Timer::Duration(Timer::Now() - solve_start).count();
      source_iterations = ksp.NumTotalMultIterations() - iterations_before;
      raw_linear_solves += ksp.NumTotalMult() - solves_before;
      raw_linear_iterations += ksp.NumTotalMultIterations() - iterations_before;
      MFEM_VERIFY(!archive_stream_only || ksp.GetConverged(),
                  "Refusing to archive an unconverged response source " << idx << "!");
    }

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), V[step]), rhs_norm);

    // Compute E = -∇V on the true dofs.
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);

    if (post_op.NeedsRecoveredElectricFlux())
    {
      Mpi::Print(" Recovering electric flux for interface postprocessing\n");
      estimator->RecoverFlux(E, D);
      post_op.SetRecoveredElectricFlux(D);
      if (!D_basis.empty())
      {
        D_basis[step] = D;
      }
    }
    if (response_archive)
    {
      WriteArchivedVector(*response_archive, idx, ArchivedField::POTENTIAL, V[step],
                          laplace_op.GetComm());
      if (post_op.NeedsRecoveredElectricFlux())
      {
        WriteArchivedVector(*response_archive, idx, ArchivedField::FLUX, D,
                            laplace_op.GetComm());
      }
    }
    if (archive_stream_only)
    {
      if (archive_recycle_initial_guess)
      {
        previous_archive_solution = V[step];
      }
      // SetSize(0) retains MFEM's allocated capacity; Destroy releases the completed
      // field so storage cannot grow with the source count.
      V[step].Destroy();
      if (response_source_timing)
      {
        double times[2] = {source_solve_seconds,
                           Timer::Duration(Timer::Now() - source_start).count()};
        Mpi::GlobalMax(2, times, laplace_op.GetComm());
        Mpi::Print("Response source timing: index={}, iterations={}, solve_seconds={:.9e}, "
                   "total_seconds={:.9e}\n",
                   idx, source_iterations, times[0], times[1]);
      }
      step++;
      continue;
    }

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, V[step], E, idx);

    EnergyData raw_energies;
    if (response_correction)
    {
      raw_energies = post_op.GetCachedElectrostaticEnergies();
    }

    if (response_correction)
    {
      BlockTimer response_timer(Timer::POSTPRO_RESPONSE);
      SurfaceResponseOperator::ElectrostaticResponse response;
      if (postprocess_response)
      {
        BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
        response = response_correction->GetElectrostaticResponse(V[step]);
        auto traces = response_correction->GetSpatialPatchTraces(V[step]);
        if (!traces.empty())
        {
          spatial_patch_traces.emplace_back(idx, std::move(traces));
        }
      }
      auto ApplyResponse = [&](EnergyData energies, double domain_correction,
                               const std::map<int, double> &fabricated_surface)
      {
        energies.domain += domain_correction;
        for (const auto &[interface, energy] : fabricated_surface)
        {
          auto it = energies.interfaces.find(interface);
          MFEM_VERIFY(it != energies.interfaces.end(),
                      "Response correction refers to target interface "
                          << interface << " which is not configured for postprocessing!");
          MFEM_VERIFY(
              !it->second.edge_energies.empty(),
              "Response-corrected target interface "
                  << interface << " requires EdgeDistances and EdgeAttributes or AutomaticEdges!");
          // EdgeDistances are sorted. The largest configured radius is the matching
          // distance of the coupon response model.
          it->second.energy = it->second.edge_energies.back().energy_outside + energy;
        }
        MFEM_VERIFY(energies.domain > 0.0,
                    "Response-corrected electrostatic energy is not positive!");
        return energies;
      };
      auto Unavailable = [](EnergyData energies)
      {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        energies.domain = nan;
        for (auto &[interface, data] : energies.interfaces)
        {
          (void)interface;
          data.energy = nan;
          data.edge_energies.clear();
        }
        return energies;
      };
      auto postprocessed_fixed_trace =
          postprocess_response ? ApplyResponse(raw_energies, response.domain_correction,
                                               response.fabricated_surface_energy)
                               : Unavailable(raw_energies);
      auto postprocessed_fixed_flux =
          postprocess_response
              ? ApplyResponse(raw_energies, response.domain_correction_fixed_flux,
                              response.fabricated_surface_energy_fixed_flux)
              : Unavailable(raw_energies);
      if (postprocess_response && !response.confident)
      {
        Mpi::Warning(
            "Electrostatic postprocessing-only surface-response confidence limits were "
            "exceeded: max interface trace-closure spread = {:.3e}, response-weighted "
            "local trace-closure spread = {:.3e}, trace-closure response fraction above "
            "5% = {:.3e}. Corrected values are reported, but the raw thin-metal "
            "field does not determine a closure-independent local response. A "
            "self-consistent result is preferable only when its globally coupled "
            "corrected solve remains well-conditioned and converges.\n",
            response.maximum_trace_closure_spread,
            response.response_weighted_trace_closure_spread,
            response.trace_closure_response_failure_fraction);
      }

      EnergyData corrected_energies = Unavailable(raw_energies);
      std::vector<SurfaceResponseOperator::ModelContribution> corrected_contributions;
      if (self_consistent_response)
      {
        Mpi::Print(" Solving fabrication-response corrected field\n");
        V_corrected[step] = V[step];
        const double solve_tol = ksp.GetRelTol();
        ksp.SetRelTol(response_config->solve_tol);
        ksp.SetOperator(*system_K);
        ksp.SetInitialGuess(true);
        const auto solves_before = ksp.NumTotalMult();
        const auto iterations_before = ksp.NumTotalMultIterations();
        ksp.Mult(corrected_rhs, V_corrected[step]);
        const bool corrected_converged = ksp.GetConverged();
        const double corrected_relative_residual = ksp.GetFinalRelativeResidual();
        corrected_linear_solves += ksp.NumTotalMult() - solves_before;
        corrected_linear_iterations += ksp.NumTotalMultIterations() - iterations_before;
        ksp.SetInitialGuess(iodata.solver.linear.initial_guess);
        ksp.SetOperator(*K);
        ksp.SetRelTol(solve_tol);
        if (!corrected_converged)
        {
          Mpi::Warning(
              "Self-consistent response-corrected solve did not converge (relative "
              "residual = {:.3e}); corrected energies and participations are reported "
              "as unavailable instead of evaluating the unconverged field.\n",
              corrected_relative_residual);
        }
        else
        {
          Vector E_corrected(Grad.Height()), D_corrected;
          E_corrected = 0.0;
          Grad.AddMult(V_corrected[step], E_corrected, -1.0);
          const Vector *D_corrected_ptr = nullptr;
          if (post_op.NeedsRecoveredElectricFlux())
          {
            D_corrected.SetSize(laplace_op.GetRTSpace().GetTrueVSize());
            D_corrected.UseDevice(true);
            estimator->RecoverFlux(E_corrected, D_corrected);
            D_corrected_ptr = &D_corrected;
          }

          const auto target_interfaces = response_correction->GetTargetInterfaces();
          {
            BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
            corrected_energies = post_op.GetElectrostaticEnergies(
                V_corrected[step], E_corrected, D_corrected_ptr, &target_interfaces);
          }
          SurfaceResponseOperator::ElectrostaticResponse corrected_response;
          {
            BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
            corrected_response =
                response_correction->GetElectrostaticResponse(V_corrected[step], false);
          }
          corrected_energies = ApplyResponse(std::move(corrected_energies),
                                             corrected_response.domain_correction,
                                             corrected_response.fabricated_surface_energy);
          corrected_contributions = std::move(corrected_response.model_contributions);
        }
      }

      if (response_correction->HasSurfaceResponse())
      {
        corrected_results.push_back(CorrectedResult{
            idx, std::move(raw_energies), std::move(postprocessed_fixed_trace),
            std::move(postprocessed_fixed_flux), std::move(corrected_energies),
            std::move(response.model_contributions), std::move(corrected_contributions),
            response.trace_closure_spread, response.maximum_trace_closure_spread,
            response.response_weighted_trace_closure_spread,
            response.trace_closure_response_failure_fraction, postprocess_response,
            self_consistent_response, response.confident});
      }
    }

    // Keep AMR driven by the historical raw thin-metal solution. This ensures enabling
    // response correction does not alter the mesh sequence or any ordinary output.
    Mpi::Print(" Updating solution error estimates\n");
    if (post_op.NeedsRecoveredElectricFlux())
    {
      estimator->AddErrorIndicator(E, D, total_domain_energy, indicator);
    }
    else
    {
      estimator->AddErrorIndicator(E, total_domain_energy, indicator);
    }

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix only for equipotential terminal solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveLinearSolverMetadata(laplace_op.GetComm(), raw_linear_solves, raw_linear_iterations);
  if (iodata.boundaries.prescribed_potential.empty())
  {
    PostprocessTerminals(post_op, laplace_op.GetSources(), V);
  }
  else if (iodata.solver.electrostatic.response_matrix && !archive_stream_only)
  {
    PostprocessResponseMatrix(post_op, laplace_op, Grad, V, D_basis);
  }

  if (root && !spatial_patch_traces.empty())
  {
    using VT = Units::ValueType;
    TableWithCSVFile traces(post_dir / "surface-response-traces.csv");
    traces.table.insert(Column("source", "i", 0, 0, 2, ""));
    traces.table.insert(Column("patch", "patch", 0, 0, 2, ""));
    traces.table.insert(Column("model", "model", 0, 0, 2, ""));
    traces.table.insert(Column("coefficient", "coefficient", 0, 0, 2, ""));
    traces.table.insert(Column("conductor_state", "conductor state", 0, 0, 2, ""));
    traces.table.insert("value", "value (V)");
    for (const auto &[source, entries] : spatial_patch_traces)
    {
      for (const auto &entry : entries)
      {
        for (std::size_t i = 0; i < entry.coefficients.size(); i++)
        {
          traces.table["source"] << source;
          traces.table["patch"] << entry.patch + 1;
          traces.table["model"] << entry.model;
          traces.table["coefficient"] << static_cast<int>(i) + 1;
          traces.table["conductor_state"]
              << (static_cast<int>(i) >= entry.contour_size ? 1 : 0);
          traces.table["value"]
              << iodata.units.Dimensionalize<VT::VOLTAGE>(entry.coefficients[i]);
        }
      }
    }
    traces.WriteFullTableTrunc();
  }

  if (root && !corrected_results.empty())
  {
    using VT = Units::ValueType;
    TableWithCSVFile output(post_dir / "surface-Q-corrected.csv");
    output.table.insert(Column("source", "i", 0, 0, 2, ""));
    output.table.insert("domain_raw", "E_elec raw (J)");
    output.table.insert("domain_postprocessed_fixed_trace",
                        "E_elec postprocessed fixed-trace (J)");
    output.table.insert("domain_postprocessed_fixed_flux",
                        "E_elec postprocessed fixed-flux (J)");
    output.table.insert("domain_corrected", "E_elec corrected (J)");
    output.table.insert("maximum_trace_closure_spread", "max trace closure spread");
    output.table.insert("weighted_trace_closure_spread",
                        "response-weighted local trace closure spread");
    output.table.insert("trace_closure_failure_fraction",
                        "trace-closure response fraction above limit");
    output.table.insert("confidence_pass", "confidence pass");
    output.table["confidence_pass"].print_as_int = true;
    const auto &interfaces = corrected_results.front().raw.interfaces;
    for (const auto &[interface, data] : interfaces)
    {
      output.table.insert(fmt::format("energy_raw_{}", interface),
                          fmt::format("E_surf raw[{}] (J)", interface));
      output.table.insert(fmt::format("participation_raw_{}", interface),
                          fmt::format("p_surf raw[{}]", interface));
      output.table.insert(fmt::format("quality_raw_{}", interface),
                          fmt::format("Q_surf raw[{}]", interface));
      output.table.insert(
          fmt::format("energy_postprocessed_fixed_trace_{}", interface),
          fmt::format("E_surf postprocessed fixed-trace[{}] (J)", interface));
      output.table.insert(
          fmt::format("participation_postprocessed_fixed_trace_{}", interface),
          fmt::format("p_surf postprocessed fixed-trace[{}]", interface));
      output.table.insert(fmt::format("quality_postprocessed_fixed_trace_{}", interface),
                          fmt::format("Q_surf postprocessed fixed-trace[{}]", interface));
      output.table.insert(
          fmt::format("energy_postprocessed_fixed_flux_{}", interface),
          fmt::format("E_surf postprocessed fixed-flux[{}] (J)", interface));
      output.table.insert(
          fmt::format("participation_postprocessed_fixed_flux_{}", interface),
          fmt::format("p_surf postprocessed fixed-flux[{}]", interface));
      output.table.insert(fmt::format("quality_postprocessed_fixed_flux_{}", interface),
                          fmt::format("Q_surf postprocessed fixed-flux[{}]", interface));
      output.table.insert(fmt::format("energy_corrected_{}", interface),
                          fmt::format("E_surf corrected[{}] (J)", interface));
      output.table.insert(fmt::format("participation_corrected_{}", interface),
                          fmt::format("p_surf corrected[{}]", interface));
      output.table.insert(fmt::format("quality_corrected_{}", interface),
                          fmt::format("Q_surf corrected[{}]", interface));
      output.table.insert(fmt::format("trace_closure_spread_{}", interface),
                          fmt::format("trace closure spread[{}]", interface));
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto &result : corrected_results)
    {
      output.table["source"] << result.source;
      output.table["domain_raw"]
          << iodata.units.Dimensionalize<VT::ENERGY>(result.raw.domain);
      output.table["domain_postprocessed_fixed_trace"]
          << iodata.units.Dimensionalize<VT::ENERGY>(
                 result.postprocessed_fixed_trace.domain);
      output.table["domain_postprocessed_fixed_flux"]
          << iodata.units.Dimensionalize<VT::ENERGY>(
                 result.postprocessed_fixed_flux.domain);
      output.table["domain_corrected"]
          << iodata.units.Dimensionalize<VT::ENERGY>(result.corrected.domain);
      output.table["maximum_trace_closure_spread"]
          << (result.has_postprocessed ? result.maximum_trace_closure_spread : nan);
      output.table["weighted_trace_closure_spread"]
          << (result.has_postprocessed ? result.response_weighted_trace_closure_spread
                                       : nan);
      output.table["trace_closure_failure_fraction"]
          << (result.has_postprocessed ? result.trace_closure_response_failure_fraction
                                       : nan);
      output.table["confidence_pass"]
          << (result.has_postprocessed ? (result.confident ? 1.0 : 0.0) : nan);
      MFEM_VERIFY(
          result.raw.interfaces.size() == interfaces.size() &&
              result.postprocessed_fixed_trace.interfaces.size() == interfaces.size() &&
              result.postprocessed_fixed_flux.interfaces.size() == interfaces.size() &&
              result.corrected.interfaces.size() == interfaces.size(),
          "Inconsistent corrected surface response entries!");
      for (const auto &[interface, raw] : result.raw.interfaces)
      {
        const auto fixed_trace = result.postprocessed_fixed_trace.interfaces.at(interface);
        const auto fixed_flux = result.postprocessed_fixed_flux.interfaces.at(interface);
        const auto corrected = result.corrected.interfaces.at(interface);
        const double p_raw = raw.energy / result.raw.domain;
        const double p_fixed_trace =
            fixed_trace.energy / result.postprocessed_fixed_trace.domain;
        const double p_fixed_flux =
            fixed_flux.energy / result.postprocessed_fixed_flux.domain;
        const double p_corrected = corrected.energy / result.corrected.domain;
        auto Quality = [](double participation, double loss_tangent)
        {
          return participation == 0.0 || loss_tangent == 0.0
                     ? mfem::infinity()
                     : 1.0 / (participation * loss_tangent);
        };
        output.table[fmt::format("energy_raw_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(raw.energy);
        output.table[fmt::format("participation_raw_{}", interface)] << p_raw;
        output.table[fmt::format("quality_raw_{}", interface)]
            << Quality(p_raw, raw.loss_tangent);
        output.table[fmt::format("energy_postprocessed_fixed_trace_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(fixed_trace.energy);
        output.table[fmt::format("participation_postprocessed_fixed_trace_{}", interface)]
            << p_fixed_trace;
        output.table[fmt::format("quality_postprocessed_fixed_trace_{}", interface)]
            << Quality(p_fixed_trace, fixed_trace.loss_tangent);
        output.table[fmt::format("energy_postprocessed_fixed_flux_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(fixed_flux.energy);
        output.table[fmt::format("participation_postprocessed_fixed_flux_{}", interface)]
            << p_fixed_flux;
        output.table[fmt::format("quality_postprocessed_fixed_flux_{}", interface)]
            << Quality(p_fixed_flux, fixed_flux.loss_tangent);
        output.table[fmt::format("energy_corrected_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(corrected.energy);
        output.table[fmt::format("participation_corrected_{}", interface)] << p_corrected;
        output.table[fmt::format("quality_corrected_{}", interface)]
            << Quality(p_corrected, corrected.loss_tangent);
        const auto closure_spread = result.trace_closure_spread.find(interface);
        output.table[fmt::format("trace_closure_spread_{}", interface)]
            << (!result.has_postprocessed ||
                        closure_spread == result.trace_closure_spread.end()
                    ? nan
                    : closure_spread->second);
      }
    }
    output.WriteFullTableTrunc();

    TableWithCSVFile model_output(post_dir / "surface-response-model-energy.csv");
    model_output.table.insert(Column("source", "source", 0, 0, 2, ""));
    model_output.table.insert(Column("evaluation", "evaluation", 0, 0, 2, ""));
    model_output.table.insert(Column("model", "model", 0, 0, 2, ""));
    model_output.table.insert(Column("patch_count", "patch count", 0, 0, 2, ""));
    model_output.table.insert("patch_weight", "patch weight");
    model_output.table.insert("domain_correction", "domain correction (J)");
    for (const auto &[interface, data] : interfaces)
    {
      (void)data;
      model_output.table.insert(
          fmt::format("interface_{}", interface),
          fmt::format("fabricated surface energy[{}] (J)", interface));
    }
    auto AppendContribution = [&](int source, int evaluation,
                                  const SurfaceResponseOperator::ModelContribution &data,
                                  bool available, bool fixed_flux)
    {
      model_output.table["source"] << source;
      model_output.table["evaluation"] << evaluation;
      model_output.table["model"] << data.model;
      model_output.table["patch_count"] << data.patch_count;
      model_output.table["patch_weight"] << data.patch_weight;
      model_output.table["domain_correction"]
          << (available ? iodata.units.Dimensionalize<VT::ENERGY>(
                              fixed_flux ? data.domain_correction_fixed_flux
                                         : data.domain_correction)
                        : nan);
      for (const auto &[interface, interface_data] : interfaces)
      {
        (void)interface_data;
        const auto &energies = fixed_flux ? data.fabricated_surface_energy_fixed_flux
                                          : data.fabricated_surface_energy;
        const auto energy = energies.find(interface);
        model_output.table[fmt::format("interface_{}", interface)]
            << (available ? iodata.units.Dimensionalize<VT::ENERGY>(
                                energy != energies.end() ? energy->second : 0.0)
                          : nan);
      }
    };
    for (const auto &result : corrected_results)
    {
      for (const auto &contribution : result.raw_model_contributions)
      {
        AppendContribution(result.source, 0, contribution, result.has_postprocessed, false);
        AppendContribution(result.source, 1, contribution, result.has_postprocessed, true);
      }
      for (const auto &contribution : result.corrected_model_contributions)
      {
        AppendContribution(result.source, 2, contribution, result.has_self_consistent,
                           false);
      }
    }
    model_output.WriteFullTableTrunc();
  }
  post_op.MeasureFinalize(indicator);
  if (self_consistent_response)
  {
    SaveSurfaceResponseSolverMetadata(laplace_op.GetComm(), "Electrostatic",
                                      corrected_linear_solves, corrected_linear_iterations);
  }
  if (response_correction)
  {
    SaveMetadata(*response_correction);
  }
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void ElectrostaticSolver::PostprocessArchivedResponseMatrix(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op, const LaplaceOperator &laplace_op,
    const Operator &Grad, const std::filesystem::path &archive, int block_size) const
{
  const auto &sources = laplace_op.GetSources();
  std::vector<int> basis_indices;
  basis_indices.reserve(sources.size());
  for (const auto &[idx, data] : sources)
  {
    (void)data;
    basis_indices.push_back(idx);
  }
  const std::size_t basis_size = basis_indices.size();
  MFEM_VERIFY(basis_size > 0 && block_size > 0,
              "Archived response reduction requires sources and a positive block size!");
  int archive_has_flux = std::filesystem::is_regular_file(ArchivePath(
                             archive, basis_indices.front(),
                             Mpi::Rank(laplace_op.GetComm()), ArchivedField::FLUX))
                             ? 1
                             : 0;
  const int local_archive_has_flux = archive_has_flux;
  Mpi::GlobalMin(1, &archive_has_flux, laplace_op.GetComm());
  int maximum_archive_has_flux = local_archive_has_flux;
  Mpi::GlobalMax(1, &maximum_archive_has_flux, laplace_op.GetComm());
  MFEM_VERIFY(archive_has_flux == maximum_archive_has_flux,
              "Response archive flux fields are inconsistent across MPI ranks!");

  TableWithCSVFile surface_output, domain_output;
  if (root)
  {
    surface_output = TableWithCSVFile(post_dir / "surface-response-matrix.csv");
    surface_output.table.insert(Column("interface", "interface", 0, 0, 2, ""));
    surface_output.table.insert(Column("edge", "edge", 0, 0, 2, ""));
    surface_output.table.insert("distance", "R (m)");
    surface_output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    surface_output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    surface_output.table.insert("Q", "Q_ij (J)");
    surface_output.table.insert("Q_normal", "Q_ij normal (J)");
    surface_output.table.insert("Q_tangential", "Q_ij tangential (J)");
    surface_output.table.insert("Q_total", "Q_total_ij (J)");
    surface_output.table.insert("Q_total_normal", "Q_total_ij normal (J)");
    surface_output.table.insert("Q_total_tangential", "Q_total_ij tangential (J)");

    domain_output = TableWithCSVFile(post_dir / "domain-response-matrix.csv");
    domain_output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    domain_output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    domain_output.table.insert("Q", "Q_ij (J)");
    domain_output.table.reserve(basis_size * (basis_size + 1) / 2, 3);
  }

  using VT = Units::ValueType;
  auto AppendSurface = [&](int interface, double distance, std::size_t i, std::size_t j,
                           const auto &entry, int local_i, int local_j)
  {
    if (!root)
    {
      return;
    }
    surface_output.table["interface"] << interface;
    surface_output.table["edge"] << 1;
    surface_output.table["distance"] << iodata.units.Dimensionalize<VT::LENGTH>(distance);
    surface_output.table["basis_i"] << basis_indices[i];
    surface_output.table["basis_j"] << basis_indices[j];
    surface_output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(
        entry.energy_inside(local_i, local_j));
    surface_output.table["Q_normal"] << iodata.units.Dimensionalize<VT::ENERGY>(
        entry.energy_inside_normal(local_i, local_j));
    surface_output.table["Q_tangential"] << iodata.units.Dimensionalize<VT::ENERGY>(
        entry.energy_inside_tangential(local_i, local_j));
    surface_output.table["Q_total"]
        << iodata.units.Dimensionalize<VT::ENERGY>(entry.energy_total(local_i, local_j));
    surface_output.table["Q_total_normal"] << iodata.units.Dimensionalize<VT::ENERGY>(
        entry.energy_total_normal(local_i, local_j));
    surface_output.table["Q_total_tangential"] << iodata.units.Dimensionalize<VT::ENERGY>(
        entry.energy_total_tangential(local_i, local_j));
  };

  // Streaming one-pass Gram (decision 62(4)): the interface samples are cached once,
  // every archived source is read, differentiated and evaluated ONCE into its amplitude
  // rows (one row per localized interface), the domain Gram takes one mass matvec per
  // source against the resident potentials, and the interface matrices are assembled
  // from the rows at the end. The quadrature rule, order, weights, sample set and
  // ownership selection are those of the batched traversal; `block_size` no longer
  // changes the work or the result (recorded for the log only). Memory: base + N x
  // (sum over interfaces of components x Q_local + L_local) x 8 bytes.
  BlockTimer bt_reduction(Timer::POSTPRO_REDUCTION);
  auto &V_gf = post_op.GetVGridFunction().Real();
  auto &D_gf = post_op.GetDomainPostOp().D;
  const auto samples = post_op.CacheInterfaceResponseSamples();
  std::size_t local_sample_count = 0;
  for (const auto &interface_samples : samples)
  {
    local_sample_count += interface_samples.Count();
  }
  {
    // Per-rank sample counts (the load balance of the interface faces across ranks).
    long long counts[3] = {static_cast<long long>(local_sample_count),
                           static_cast<long long>(local_sample_count),
                           static_cast<long long>(local_sample_count)};
    Mpi::GlobalMin(1, &counts[0], post_op.GetComm());
    Mpi::GlobalMax(1, &counts[1], post_op.GetComm());
    Mpi::GlobalSum(1, &counts[2], post_op.GetComm());
    Mpi::Print(" Archived response reduction: {:d} sources, {:d} interfaces, quadrature "
               "samples per rank min {:d}, max {:d}, total {:d} (streaming; block size {:d} "
               "recorded)\n",
               static_cast<int>(basis_size), static_cast<int>(samples.size()), counts[0],
               counts[1], counts[2], block_size);
    for (const auto &interface_samples : samples)
    {
      long long per_interface[3] = {static_cast<long long>(interface_samples.Count()),
                                    static_cast<long long>(interface_samples.Count()),
                                    static_cast<long long>(interface_samples.Count())};
      Mpi::GlobalMin(1, &per_interface[0], post_op.GetComm());
      Mpi::GlobalMax(1, &per_interface[1], post_op.GetComm());
      Mpi::GlobalSum(1, &per_interface[2], post_op.GetComm());
      Mpi::Print("  interface {:d}: samples per rank min {:d}, max {:d}, total {:d}\n",
                 interface_samples.interface_index, per_interface[0], per_interface[1],
                 per_interface[2]);
    }
  }
  std::vector<std::vector<double>> rows(samples.size());
  for (std::size_t k = 0; k < samples.size(); k++)
  {
    rows[k].assign(basis_size * samples[k].RowSize(), 0.0);
  }
  std::vector<Vector> V_local(basis_size);
  mfem::DenseMatrix local_domain(static_cast<int>(basis_size));
  local_domain = 0.0;
  Vector E_basis(Grad.Height()), D_source;
  for (std::size_t i = 0; i < basis_size; i++)
  {
    const int source = basis_indices[i];
    {
      BlockTimer bt_read(Timer::POSTPRO_REDUCTION_READ);
      Vector V = ReadArchivedVector(archive, source, ArchivedField::POTENTIAL,
                                    laplace_op.GetComm());
      MFEM_VERIFY(V.Size() == Grad.Width(),
                  "Archived potential field has an incompatible local size!");
      if (archive_has_flux)
      {
        D_source =
            ReadArchivedVector(archive, source, ArchivedField::FLUX, laplace_op.GetComm());
      }
      E_basis = 0.0;
      Grad.AddMult(V, E_basis, -1.0);
      V_gf.SetFromTrueDofs(V);
      V_local[i] = V_gf;
    }
    {
      BlockTimer bt_eval(Timer::POSTPRO_REDUCTION_EVAL);
      post_op.SetInterfaceResponseField(E_basis, archive_has_flux ? &D_source : nullptr);
      for (std::size_t k = 0; k < samples.size(); k++)
      {
        post_op.EvaluateInterfaceResponseRow(samples[k], rows[k].data() + i * samples[k].RowSize());
      }
    }
    {
      BlockTimer bt_domain(Timer::POSTPRO_REDUCTION_DOMAIN);
      post_op.GetDomainPostOp().M_elec->Mult(V_local[i], D_gf);
      for (std::size_t j = 0; j <= i; j++)
      {
        local_domain(static_cast<int>(j), static_cast<int>(i)) =
            0.5 * linalg::LocalDot(V_local[j], D_gf);
      }
    }
    Mpi::Print(" Archived response source {:d}/{:d}\n", static_cast<int>(i + 1),
               static_cast<int>(basis_size));
  }
  Mpi::GlobalSum(local_domain.Height() * local_domain.Width(), local_domain.GetData(),
                 post_op.GetComm());

  std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>> surface_matrices;
  {
    BlockTimer bt_gram(Timer::POSTPRO_REDUCTION_GRAM);
    for (std::size_t k = 0; k < samples.size(); k++)
    {
      surface_matrices.emplace(samples[k].interface_index,
                               post_op.AssembleInterfaceResponseMatrices(
                                   samples[k], rows[k].data(), static_cast<int>(basis_size)));
      std::vector<double>().swap(rows[k]);
    }
  }
  for (std::size_t i = 0; i < basis_size; i++)
  {
    for (std::size_t j = i; j < basis_size; j++)
    {
      for (const auto &[interface, entries] : surface_matrices)
      {
        for (const auto &entry : entries)
        {
          AppendSurface(interface, entry.distance, i, j, entry, static_cast<int>(i),
                        static_cast<int>(j));
        }
      }
      if (root)
      {
        domain_output.table["basis_i"] << basis_indices[i];
        domain_output.table["basis_j"] << basis_indices[j];
        domain_output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(
            local_domain(static_cast<int>(i), static_cast<int>(j)));
      }
    }
  }
  Mpi::Print(" Archived response reduction complete: {:d} sources streamed once\n",
             static_cast<int>(basis_size));
  if (root)
  {
    surface_output.WriteFullTableTrunc();
    domain_output.WriteFullTableTrunc();
  }
}

void ElectrostaticSolver::PostprocessResponseMatrix(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op, const LaplaceOperator &laplace_op,
    const Operator &Grad, const std::vector<Vector> &V, const std::vector<Vector> &D) const
{
  using LocalEnergy = SurfacePostOperator::InterfaceLocalEdgeEnergy;
  using EnergyMap = std::map<int, std::vector<LocalEnergy>>;

  const auto &sources = laplace_op.GetSources();
  MFEM_VERIFY(V.size() == sources.size(),
              "Unexpected prescribed-potential basis field count!");
  MFEM_VERIFY(D.empty() || D.size() == V.size(),
              "Unexpected recovered electric flux basis field count!");

  std::vector<int> basis_indices;
  basis_indices.reserve(sources.size());
  for (const auto &[idx, data] : sources)
  {
    basis_indices.push_back(idx);
  }

  Vector E_sum(Grad.Height()), D_sum;
  if (!D.empty())
  {
    D_sum.SetSize(D.front().Size());
    D_sum.UseDevice(true);
  }

  auto Aggregate = [](EnergyMap energies)
  {
    for (auto &[interface, local] : energies)
    {
      (void)interface;
      std::map<double, LocalEnergy> by_distance;
      for (const auto &energy : local)
      {
        auto [it, inserted] = by_distance.emplace(energy.distance, energy);
        auto &aggregate = it->second;
        if (inserted)
        {
          aggregate.edge = 1;
          continue;
        }
        aggregate.energy_inside += energy.energy_inside;
        aggregate.energy_total += energy.energy_total;
        for (int component = 0; component < 2; component++)
        {
          aggregate.energy_inside_polarized[component] +=
              energy.energy_inside_polarized[component];
          aggregate.energy_total_polarized[component] +=
              energy.energy_total_polarized[component];
        }
      }
      local.clear();
      local.reserve(by_distance.size());
      for (auto &[distance, energy] : by_distance)
      {
        (void)distance;
        local.push_back(std::move(energy));
      }
    }
    return energies;
  };

  auto Evaluate = [&](std::size_t i, std::optional<std::size_t> j = std::nullopt)
  {
    E_sum = 0.0;
    Grad.AddMult(V[i], E_sum, -1.0);
    if (j)
    {
      Grad.AddMult(V[*j], E_sum, -1.0);
    }

    const Vector *D_ptr = nullptr;
    if (!D.empty())
    {
      D_sum = D[i];
      if (j)
      {
        D_sum.Add(1.0, D[*j]);
      }
      D_ptr = &D_sum;
    }
    auto energies = post_op.GetInterfaceLocalEdgeElectricFieldEnergies(E_sum, D_ptr);
    if (iodata.solver.electrostatic.aggregate_response_matrix)
    {
      return Aggregate(std::move(energies));
    }
    return energies;
  };

  Mpi::Print("\nAssembling {} interface response matrix for {:d} basis fields\n",
             iodata.solver.electrostatic.aggregate_response_matrix ? "aggregated"
                                                                   : "localized",
             static_cast<int>(V.size()));
  std::vector<EnergyMap> diagonal;
  if (!iodata.solver.electrostatic.aggregate_response_matrix)
  {
    diagonal.reserve(V.size());
    for (std::size_t i = 0; i < V.size(); i++)
    {
      diagonal.push_back(Evaluate(i));
    }
  }

  TableWithCSVFile output;
  if (root)
  {
    output = TableWithCSVFile(post_dir / "surface-response-matrix.csv");
    output.table.insert(Column("interface", "interface", 0, 0, 2, ""));
    output.table.insert(Column("edge", "edge", 0, 0, 2, ""));
    output.table.insert("distance", "R (m)");
    output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    output.table.insert("Q", "Q_ij (J)");
    output.table.insert("Q_normal", "Q_ij normal (J)");
    output.table.insert("Q_tangential", "Q_ij tangential (J)");
    output.table.insert("Q_total", "Q_total_ij (J)");
    output.table.insert("Q_total_normal", "Q_total_ij normal (J)");
    output.table.insert("Q_total_tangential", "Q_total_ij tangential (J)");
    output.table.reserve(V.size() * (V.size() + 1) / 2, 11);
  }

  using VT = Units::ValueType;
  auto Append = [&](int interface, int edge, double distance, std::size_t i, std::size_t j,
                    double q, double q_normal, double q_tangential, double q_total,
                    double q_total_normal, double q_total_tangential)
  {
    if (!root)
    {
      return;
    }
    output.table["interface"] << interface;
    output.table["edge"] << edge;
    output.table["distance"] << iodata.units.Dimensionalize<VT::LENGTH>(distance);
    output.table["basis_i"] << basis_indices[i];
    output.table["basis_j"] << basis_indices[j];
    output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(q);
    output.table["Q_normal"] << iodata.units.Dimensionalize<VT::ENERGY>(q_normal);
    output.table["Q_tangential"] << iodata.units.Dimensionalize<VT::ENERGY>(q_tangential);
    output.table["Q_total"] << iodata.units.Dimensionalize<VT::ENERGY>(q_total);
    output.table["Q_total_normal"]
        << iodata.units.Dimensionalize<VT::ENERGY>(q_total_normal);
    output.table["Q_total_tangential"]
        << iodata.units.Dimensionalize<VT::ENERGY>(q_total_tangential);
  };

  const std::size_t pair_count = V.size() * (V.size() + 1) / 2;
  const std::size_t progress_interval = std::max<std::size_t>(pair_count / 20, 1);
  std::size_t completed_pairs = 0;
  if (iodata.solver.electrostatic.aggregate_response_matrix)
  {
    Mpi::Print(" Using batched surface Gram-matrix assembly\n");
    std::vector<Vector> E_basis(V.size());
    for (std::size_t i = 0; i < V.size(); i++)
    {
      E_basis[i].SetSize(Grad.Height());
      E_basis[i] = 0.0;
      Grad.AddMult(V[i], E_basis[i], -1.0);
    }
    const auto matrices =
        post_op.GetInterfaceElectricFieldEnergyMatrices(E_basis, D.empty() ? nullptr : &D);
    for (std::size_t i = 0; i < V.size(); i++)
    {
      for (std::size_t j = i; j < V.size(); j++)
      {
        for (const auto &[interface, entries] : matrices)
        {
          for (const auto &entry : entries)
          {
            Append(interface, 1, entry.distance, i, j, entry.energy_inside(i, j),
                   entry.energy_inside_normal(i, j), entry.energy_inside_tangential(i, j),
                   entry.energy_total(i, j), entry.energy_total_normal(i, j),
                   entry.energy_total_tangential(i, j));
          }
        }
        completed_pairs++;
        if (completed_pairs % progress_interval == 0 || completed_pairs == pair_count)
        {
          Mpi::Print(" Interface response matrix: {:d}/{:d} basis pairs ({:.0f}%)\n",
                     static_cast<int>(completed_pairs), static_cast<int>(pair_count),
                     100.0 * completed_pairs / pair_count);
        }
      }
    }
  }
  else
  {
    for (std::size_t i = 0; i < V.size(); i++)
    {
      for (std::size_t j = i; j < V.size(); j++)
      {
        std::optional<EnergyMap> combined;
        if (i != j)
        {
          combined = Evaluate(i, j);
        }
        MFEM_VERIFY(diagonal[i].size() == diagonal[j].size() &&
                        (!combined || combined->size() == diagonal[i].size()),
                    "Inconsistent localized interface response data!");

        for (const auto &[interface, energy_i] : diagonal[i])
        {
          const auto it_j = diagonal[j].find(interface);
          const std::vector<LocalEnergy> *energy_sum = nullptr;
          if (combined)
          {
            const auto it_sum = combined->find(interface);
            MFEM_VERIFY(it_sum != combined->end(),
                        "Missing combined-field localized interface response entries!");
            energy_sum = &it_sum->second;
          }
          MFEM_VERIFY(it_j != diagonal[j].end() && energy_i.size() == it_j->second.size() &&
                          (!energy_sum || energy_i.size() == energy_sum->size()),
                      "Inconsistent localized interface response entries!");
          for (std::size_t entry = 0; entry < energy_i.size(); entry++)
          {
            const auto &ei = energy_i[entry];
            const auto &ej = it_j->second[entry];
            MFEM_VERIFY(ei.edge == ej.edge && ei.distance == ej.distance,
                        "Inconsistent localized interface response metadata!");
            if (i == j)
            {
              Append(interface, ei.edge, ei.distance, i, j, ei.energy_inside,
                     ei.energy_inside_polarized[0], ei.energy_inside_polarized[1],
                     ei.energy_total, ei.energy_total_polarized[0],
                     ei.energy_total_polarized[1]);
            }
            else
            {
              const auto &es = (*energy_sum)[entry];
              MFEM_VERIFY(ei.edge == es.edge && ei.distance == es.distance,
                          "Inconsistent localized interface response metadata!");
              Append(interface, ei.edge, ei.distance, i, j,
                     0.5 * (es.energy_inside - ei.energy_inside - ej.energy_inside),
                     0.5 * (es.energy_inside_polarized[0] - ei.energy_inside_polarized[0] -
                            ej.energy_inside_polarized[0]),
                     0.5 * (es.energy_inside_polarized[1] - ei.energy_inside_polarized[1] -
                            ej.energy_inside_polarized[1]),
                     0.5 * (es.energy_total - ei.energy_total - ej.energy_total),
                     0.5 * (es.energy_total_polarized[0] - ei.energy_total_polarized[0] -
                            ej.energy_total_polarized[0]),
                     0.5 * (es.energy_total_polarized[1] - ei.energy_total_polarized[1] -
                            ej.energy_total_polarized[1]));
            }
          }
        }
        completed_pairs++;
        if (completed_pairs % progress_interval == 0 || completed_pairs == pair_count)
        {
          Mpi::Print(" Interface response matrix: {:d}/{:d} basis pairs ({:.0f}%)\n",
                     static_cast<int>(completed_pairs), static_cast<int>(pair_count),
                     100.0 * completed_pairs / pair_count);
        }
      }
    }
  }

  if (root)
  {
    output.WriteFullTableTrunc();
  }

  Mpi::Print("Assembling domain-energy response matrix\n");
  TableWithCSVFile domain_output;
  if (root)
  {
    domain_output = TableWithCSVFile(post_dir / "domain-response-matrix.csv");
    domain_output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    domain_output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    domain_output.table.insert("Q", "Q_ij (J)");
    domain_output.table.reserve(V.size() * (V.size() + 1) / 2, 3);
  }

  auto &V_gf = post_op.GetVGridFunction().Real();
  auto &D_gf = post_op.GetDomainPostOp().D;
  std::vector<Vector> V_local(V.size());
  for (std::size_t i = 0; i < V.size(); i++)
  {
    V_gf.SetFromTrueDofs(V[i]);
    V_local[i] = V_gf;
  }
  mfem::DenseMatrix domain_matrix(static_cast<int>(V.size()));
  domain_matrix = 0.0;
  for (std::size_t i = 0; i < V.size(); i++)
  {
    post_op.GetDomainPostOp().M_elec->Mult(V_local[i], D_gf);
    for (std::size_t j = i; j < V.size(); j++)
    {
      domain_matrix(i, j) = 0.5 * linalg::LocalDot(V_local[j], D_gf);
    }
  }
  Mpi::GlobalSum(domain_matrix.Height() * domain_matrix.Width(), domain_matrix.GetData(),
                 post_op.GetComm());
  if (root)
  {
    for (std::size_t i = 0; i < V.size(); i++)
    {
      for (std::size_t j = i; j < V.size(); j++)
      {
        domain_output.table["basis_i"] << basis_indices[i];
        domain_output.table["basis_j"] << basis_indices[j];
        domain_output.table["Q"]
            << iodata.units.Dimensionalize<VT::ENERGY>(domain_matrix(i, j));
      }
    }
  }
  if (root)
  {
    domain_output.WriteFullTableTrunc();
  }
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
  mfem::DenseMatrix C(V.size()), Cm(V.size());
  for (int i = 0; i < C.Height(); i++)
  {
    // Diagonal: Cᵢᵢ = 2 Uₑ(Vᵢ) / Vᵢ² = (Vᵢᵀ K Vᵢ) / Vᵢ² (with ∀i, Vᵢ = 1)
    auto &V_gf = post_op.GetVGridFunction().Real();
    auto &D_gf = post_op.GetDomainPostOp().D;
    V_gf.SetFromTrueDofs(V[i]);
    post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    C(i, i) = Cm(i, i) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);

    // Off-diagonals: Cᵢⱼ = Uₑ(Vᵢ + Vⱼ) / (Vᵢ Vⱼ) - 1/2 (Vᵢ/Vⱼ Cᵢᵢ + Vⱼ/Vᵢ Cⱼⱼ)
    //                    = (Vⱼᵀ K Vᵢ) / (Vᵢ Vⱼ)
    for (int j = i + 1; j < C.Width(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      C(i, j) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);
      Cm(i, j) = -C(i, j);
      Cm(i, i) -= Cm(i, j);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      C(i, j) = C(j, i);
      Cm(i, j) = Cm(j, i);
      Cm(i, i) -= Cm(i, j);
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
