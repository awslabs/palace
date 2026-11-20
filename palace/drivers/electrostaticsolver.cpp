// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <Eigen/Eigenvalues>
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
#include "utils/collectiveenvironment.hpp"
#include "utils/communication.hpp"
#include "utils/enum_string.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>

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

// Decision 1471 (2) (C10): the PALACE_RESPONSE_* run-mode gates below select code paths
// with collectives (the solve, the archive streams, the source-timing reduction), so each
// is read COLLECTIVELY: the root's value on every rank (utils/collectiveenvironment.hpp).
// Their Solver.Electrostatic.ResponseArchive keys are identical on every rank.
std::filesystem::path ResponseArchivePath(const std::string &value)
{
  return std::filesystem::absolute(value).lexically_normal();
}

void ResponseArchiveConflict(const char *key, const char *variable)
{
  MFEM_ABORT("Solver.Electrostatic.ResponseArchive \""
             << key << "\" and its alias variable " << variable
             << " are both set and differ!");
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

// The extreme Ritz values of the preconditioned operator B⁻¹A seen by a CG solve, from the
// recorded coefficients alpha_k and beta_k / beta_{k-1}: the Lanczos tridiagonal of the m
// iterations has T_kk = 1 / alpha_k + (beta_k / beta_{k-1}) / alpha_{k-1} and
// T_{k,k+1} = sqrt(beta_{k+1} / beta_k) / alpha_k, and its eigenvalues are the Ritz values.
// The tridiagonal is passed to the eigenvalue solve directly (no m x m matrix). NaN when
// the record is empty or when m exceeds corrected_solve_ritz_max_iterations.
std::pair<double, double> CgRitzValueRange(const std::vector<double> &alpha,
                                           const std::vector<double> &beta_ratio)
{
  const int m = static_cast<int>(alpha.size());
  if (m == 0 || m > corrected_solve_ritz_max_iterations)
  {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }
  Eigen::VectorXd diagonal(m), subdiagonal(m - 1);
  for (int k = 0; k < m; k++)
  {
    diagonal(k) = 1.0 / alpha[k] + (k > 0 ? beta_ratio[k] / alpha[k - 1] : 0.0);
    if (k + 1 < m)
    {
      subdiagonal(k) = std::sqrt(std::abs(beta_ratio[k + 1])) / alpha[k];
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensystem;
  eigensystem.computeFromTridiagonal(diagonal, subdiagonal, Eigen::EigenvaluesOnly);
  return {eigensystem.eigenvalues().minCoeff(), eigensystem.eigenvalues().maxCoeff()};
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

CorrectedSolveRecord SolveCorrectedField(KspSolver &ksp, const Operator &K,
                                         const Operator &corrected_K, double solve_tol,
                                         bool initial_guess, const Vector &rhs, Vector &x)
{
  const double raw_tol = ksp.GetRelTol();
  ksp.SetRelTol(solve_tol);
  ksp.SetOperator(corrected_K);
  ksp.SetInitialGuess(true);
  ksp.EnableCgHistory(true);
  ksp.Mult(rhs, x);
  CorrectedSolveRecord record;
  record.converged = ksp.GetConverged();
  record.relative_residual = ksp.GetFinalRelativeResidual();
  record.iterations = static_cast<int>(ksp.GetCgAlphaHistory().size());
  record.negative_curvature_count = ksp.GetCgNegativeCurvatureCount();
  record.first_negative_curvature_iteration = ksp.GetCgFirstNegativeCurvatureIteration();
  std::tie(record.ritz_min, record.ritz_max) =
      CgRitzValueRange(ksp.GetCgAlphaHistory(), ksp.GetCgBetaRatioHistory());
  record.accepted = record.converged && record.negative_curvature_count == 0;
  ksp.EnableCgHistory(false);
  ksp.SetInitialGuess(initial_guess);
  ksp.SetOperator(K);
  ksp.SetRelTol(raw_tol);
  return record;
}

ResponseArchiveSettings
ResolveResponseArchive(const config::ElectrostaticSolverData &electrostatic, MPI_Comm comm)
{
  using Mode = config::ElectrostaticSolverData::ResponseArchiveData::Mode;
  const auto keys = electrostatic.response_archive.value_or(
      config::ElectrostaticSolverData::ResponseArchiveData{});
  ResponseArchiveSettings settings;
  auto SourceOf = [](bool from_key, bool from_variable)
  {
    return from_key ? (from_variable ? ResponseArchiveSource::KEY_AND_VARIABLE
                                     : ResponseArchiveSource::KEY)
                    : (from_variable ? ResponseArchiveSource::VARIABLE
                                     : ResponseArchiveSource::DEFAULT);
  };

  // A boolean alias variable is on when it reads "1". An archive estimate refuses
  // ARCHIVE_ONLY, REDUCE_ONLY, RECYCLE_INITIAL_GUESS and BLOCK_SIZE when they are set at
  // all (in this order), and the BlockSize / RecycleInitialGuess keys when they are given.
  auto Resolve = [&](const char *variable, const char *key,
                     const std::optional<bool> &from_key, bool refused_by_estimate,
                     bool &variable_set)
  {
    const auto value = CollectiveGetenv(comm, variable);
    variable_set = value.has_value();
    if (value && refused_by_estimate)
    {
      settings.estimate_conflicts.emplace_back(variable);
    }
    const bool from_variable = value && *value == "1";
    if (!from_key)
    {
      return from_variable;
    }
    if (value && from_variable != *from_key)
    {
      ResponseArchiveConflict(key, variable);
    }
    return *from_key;
  };
  auto ModeIs = [&](Mode mode)
  { return keys.mode ? std::optional<bool>(*keys.mode == mode) : std::nullopt; };
  bool archive_variable = false, reduce_variable = false, estimate_variable = false;
  settings.archive_only = Resolve("PALACE_RESPONSE_ARCHIVE_ONLY", "Mode",
                                  ModeIs(Mode::ARCHIVE), true, archive_variable);
  settings.reduce_only = Resolve("PALACE_RESPONSE_REDUCE_ONLY", "Mode",
                                 ModeIs(Mode::REDUCE), true, reduce_variable);
  settings.estimate_only = Resolve("PALACE_RESPONSE_ESTIMATE_ONLY", "Mode",
                                   ModeIs(Mode::ESTIMATE), false, estimate_variable);
  settings.mode_source = SourceOf(keys.mode.has_value(),
                                  archive_variable || reduce_variable || estimate_variable);
  bool recycle_variable = false;
  settings.recycle_initial_guess =
      Resolve("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "RecycleInitialGuess",
              keys.recycle_initial_guess, true, recycle_variable);
  settings.recycle_initial_guess_source =
      SourceOf(keys.recycle_initial_guess.has_value(), recycle_variable);
  if (keys.recycle_initial_guess)
  {
    settings.estimate_conflicts.emplace_back(
        "Solver.Electrostatic.ResponseArchive.RecycleInitialGuess");
  }
  bool timing_variable = false;
  settings.source_timing = Resolve("PALACE_RESPONSE_SOURCE_TIMING", "SourceTiming",
                                   keys.source_timing, false, timing_variable);
  settings.source_timing_source = SourceOf(keys.source_timing.has_value(), timing_variable);

  // The reduction's block size: the variable is parsed only where it is used (a reduction)
  // or compared with its key, as before the keys existed.
  bool block_size_variable = false;
  if (const auto value = CollectiveGetenv(comm, "PALACE_RESPONSE_BLOCK_SIZE"))
  {
    settings.estimate_conflicts.emplace_back("PALACE_RESPONSE_BLOCK_SIZE");
    if (settings.reduce_only || keys.block_size)
    {
      settings.block_size = std::stoi(*value);
      MFEM_VERIFY(settings.block_size > 0, "PALACE_RESPONSE_BLOCK_SIZE must be positive!");
      if (keys.block_size && *keys.block_size != settings.block_size)
      {
        ResponseArchiveConflict("BlockSize", "PALACE_RESPONSE_BLOCK_SIZE");
      }
      block_size_variable = true;
    }
  }
  if (keys.block_size)
  {
    settings.estimate_conflicts.emplace_back(
        "Solver.Electrostatic.ResponseArchive.BlockSize");
    settings.block_size = *keys.block_size;
  }
  settings.block_size_source = SourceOf(keys.block_size.has_value(), block_size_variable);

  // A path given by both the key and the variable must name the same path. The comparison
  // is made on the root and broadcast (a relative spelling is taken against the root's
  // working directory), so that every rank refuses or none does.
  auto PathsDiffer = [&](const std::filesystem::path &a, const std::filesystem::path &b)
  {
    bool differ = Mpi::Root(comm) && a != b;
    Mpi::Broadcast(1, &differ, 0, comm);
    return differ;
  };

  // The archive directory: an empty variable is unset.
  if (const auto value = CollectiveGetenv(comm, "PALACE_RESPONSE_ARCHIVE_DIR");
      value && !value->empty())
  {
    settings.directory = ResponseArchivePath(*value);
  }
  settings.directory_source =
      SourceOf(keys.directory.has_value(), settings.directory.has_value());
  if (keys.directory)
  {
    const auto from_key = ResponseArchivePath(*keys.directory);
    if (settings.directory && PathsDiffer(*settings.directory, from_key))
    {
      ResponseArchiveConflict("Directory", "PALACE_RESPONSE_ARCHIVE_DIR");
    }
    settings.directory = from_key;
  }

  // The estimate request file (used as given; an empty variable is set).
  settings.estimate_request = CollectiveGetenv(comm, "PALACE_RESPONSE_ESTIMATE_REQUEST");
  settings.estimate_request_source =
      SourceOf(keys.estimate_request.has_value(), settings.estimate_request.has_value());
  if (keys.estimate_request)
  {
    if (settings.estimate_request &&
        PathsDiffer(ResponseArchivePath(*settings.estimate_request),
                    ResponseArchivePath(*keys.estimate_request)))
    {
      ResponseArchiveConflict("EstimateRequest", "PALACE_RESPONSE_ESTIMATE_REQUEST");
    }
    settings.estimate_request = *keys.estimate_request;
  }
  return settings;
}

std::string ResponseArchiveSummary(const ResponseArchiveSettings &settings)
{
  auto Source = [](ResponseArchiveSource source)
  {
    switch (source)
    {
      case ResponseArchiveSource::KEY:
        return "key";
      case ResponseArchiveSource::VARIABLE:
        return "variable";
      case ResponseArchiveSource::KEY_AND_VARIABLE:
        return "key and variable";
      case ResponseArchiveSource::DEFAULT:
        break;
    }
    return "default";
  };
  std::vector<std::string> modes;
  if (settings.archive_only)
  {
    modes.emplace_back("Archive");
  }
  if (settings.reduce_only)
  {
    modes.emplace_back("Reduce");
  }
  if (settings.estimate_only)
  {
    modes.emplace_back("Estimate");
  }
  std::string mode = modes.empty() ? "None" : modes.front();
  for (std::size_t i = 1; i < modes.size(); i++)
  {
    mode += "+" + modes[i];
  }
  return fmt::format(
      "Response-archive settings: Mode={} ({}), Directory={} ({}), BlockSize={} ({}), "
      "SourceTiming={} ({}), RecycleInitialGuess={} ({}), EstimateRequest={} ({})",
      mode, Source(settings.mode_source),
      settings.directory ? settings.directory->string() : std::string("-"),
      Source(settings.directory_source), settings.block_size,
      Source(settings.block_size_source), settings.source_timing,
      Source(settings.source_timing_source), settings.recycle_initial_guess,
      Source(settings.recycle_initial_guess_source),
      settings.estimate_request ? *settings.estimate_request : std::string("-"),
      Source(settings.estimate_request_source));
}

namespace
{

// The key path only: a Solver.Electrostatic.ResponseArchive key that has no effect in the
// resolved mode is refused, and so is a key-given Archive / Reduce mode without a
// Directory. The estimate refuses its own (RecycleInitialGuess, BlockSize). The alias
// variables keep their silent semantics.
void ValidateResponseArchiveKeys(const config::ElectrostaticSolverData &electrostatic,
                                 const ResponseArchiveSettings &settings)
{
  using Mode = config::ElectrostaticSolverData::ResponseArchiveData::Mode;
  if (!electrostatic.response_archive)
  {
    return;
  }
  const auto &keys = *electrostatic.response_archive;
  auto Inert = [](bool given, bool effective, const char *key, const char *mode)
  {
    MFEM_VERIFY(!given || effective, "Solver.Electrostatic.ResponseArchive \""
                                         << key << "\" has no effect outside Mode \""
                                         << mode << "\"!");
  };
  Inert(keys.source_timing.has_value(), settings.archive_only, "SourceTiming", "Archive");
  Inert(keys.recycle_initial_guess.has_value(),
        settings.archive_only || settings.estimate_only, "RecycleInitialGuess", "Archive");
  Inert(keys.block_size.has_value(), settings.reduce_only || settings.estimate_only,
        "BlockSize", "Reduce");
  Inert(keys.estimate_request.has_value(), settings.estimate_only, "EstimateRequest",
        "Estimate");
  MFEM_VERIFY(!keys.mode || (*keys.mode != Mode::ARCHIVE && *keys.mode != Mode::REDUCE) ||
                  settings.directory,
              "Solver.Electrostatic.ResponseArchive \"Mode\" \""
                  << (*keys.mode == Mode::ARCHIVE ? "Archive" : "Reduce")
                  << "\" requires a \"Directory\" (or PALACE_RESPONSE_ARCHIVE_DIR)!");
}

}  // namespace

void ValidateArchiveEstimateOptions(const IoData &iodata, MPI_Comm comm, bool check_output)
{
  const auto settings = ResolveResponseArchive(iodata.solver.electrostatic, comm);
  ValidateResponseArchiveKeys(iodata.solver.electrostatic, settings);
  if (!settings.estimate_only)
  {
    return;
  }
  MFEM_VERIFY(settings.estimate_conflicts.empty(),
              "Archive estimation conflicts with " << settings.estimate_conflicts.front());
  const auto &archive = settings.directory;
  MFEM_VERIFY(archive && std::filesystem::is_directory(*archive),
              "Archive estimation requires an existing response-archive Directory "
              "(PALACE_RESPONSE_ARCHIVE_DIR)!");
  const auto &estimate_request = settings.estimate_request;
  MFEM_VERIFY(estimate_request, "Archive estimation requires a response-archive "
                                "EstimateRequest (PALACE_RESPONSE_ESTIMATE_REQUEST)!");
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
  Disjoint(*estimate_request);
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
                                            const std::filesystem::path &archive,
                                            const std::string &estimate_request) const
{
  using json = nlohmann::json;
  using VT = Units::ValueType;
  const auto comm = laplace_op.GetComm();
  const auto &Grad = laplace_op.GetGradMatrix();
  std::ifstream request_stream(estimate_request);
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
  // P3-16: every excitation below is evaluated by the batched interface response, so a
  // summed sheet sample is refused before the archive is read (whatever the response-matrix
  // keys of the configuration).
  post_op.VerifyBatchedInterfaceResponse();
  GradFluxErrorEstimator<Vector> estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg,
      iodata.boundaries.GetSheetAttributes());
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
      // P3-21: the recovered D lives on the estimator's flux recovery space (broken across
      // the interior sheet faces when the mesh has any, else the RT space above).
      {"GlobalFluxTrueDofs", estimator.GetFluxSpace().GlobalTrueVSize()},
      {"FluxSpaceBroken", estimator.GetFluxSpace().IsBroken()},
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
                    {"FluxTrueDofs", estimator.GetFluxSpace().GetTrueVSize()},
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
  Vector source_bc, source_rhs, e(Grad.Height()),
      d(estimator.GetFluxSpace().GetTrueVSize());
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
  const auto response_archive_settings =
      ResolveResponseArchive(iodata.solver.electrostatic, laplace_op.GetComm());
  if (response_archive_settings.Active())
  {
    Mpi::Print(laplace_op.GetComm(), "\n{}\n",
               ResponseArchiveSummary(response_archive_settings));
  }
  const bool archive_reduce_only = response_archive_settings.reduce_only;
  decltype(laplace_op.GetStiffnessMatrix()) K;
  if (!archive_reduce_only)
  {
    K = laplace_op.GetStiffnessMatrix();
  }
  if (response_archive_settings.estimate_only)
  {
    MFEM_VERIFY(!archive_reduce_only,
                "The response-archive modes Estimate (PALACE_RESPONSE_ESTIMATE_ONLY) and "
                "Reduce (PALACE_RESPONSE_REDUCE_ONLY) exclude each other!");
    ValidateArchiveEstimateOptions(iodata, laplace_op.GetComm(), false);
    SaveMetadata(laplace_op.GetH1Spaces());
    auto indicator =
        EstimateArchivedFields(laplace_op, *K, *response_archive_settings.directory,
                               *response_archive_settings.estimate_request);
    SaveLinearSolverMetadata(laplace_op.GetComm(), 0, 0);
    return {std::move(indicator), laplace_op.GlobalTrueVSize()};
  }
  const auto *response_config = iodata.solver.electrostatic.response_correction
                                    ? &*iodata.solver.electrostatic.response_correction
                                    : nullptr;
  // Non-const: a whole-run geometry fallback (below) turns every response branch off.
  bool postprocess_response = response_config && response_config->IncludesPostprocessing();
  bool self_consistent_response =
      response_config && response_config->IncludesSelfConsistent();
  // Decision 352 follow-up (1): the per-patch energies (surface-response-patch-energy.csv)
  // are evaluated and written only on request; every other output is unchanged.
  bool patch_energy_response = response_config && response_config->patch_energy;
  MFEM_VERIFY(!archive_reduce_only || !response_config,
              "The response-archive mode Reduce (PALACE_RESPONSE_REDUCE_ONLY) reduces "
              "archived fields of a configuration without ResponseCorrection!");
  std::unique_ptr<SurfaceResponseOperator> response_correction;
  std::unique_ptr<SumOperator> corrected_K;
  const Operator *system_K = K.get();
  // The patch table (surface-response-patches.csv): written at construction and rewritten
  // when the conductor-consistency gate (decision 277) excludes a patch at solve time (its
  // weight becomes 0, like a DomainBoundary cell).
  auto WritePatchTable = [&]()
  {
    if (root && !response_correction->GetPatchAssignments().empty())
    {
      TableWithCSVFile patches(
          response_correction->DiagnosticsPath(post_dir, "surface-response-patches.csv"));
      patches.table.insert(Column::Integer("model", "model"));
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
  };
  // Geometry fallbacks (USER decision 1000 (b) / 1001): a whole-run raw fallback (no model
  // / patch applied: NoMatchedFeature, AllPatchesOutside, NoPerimeter) saves the operator's
  // record + banner and runs the thin solve WITHOUT the operator (every other output = the
  // run without ResponseCorrection); nothing else below sees a response_correction.
  if (response_config)
  {
    if (response_config->impedance_metal_rule ==
        config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceMetalRule::RECORD)
    {
      // The impedance-metal rule (decisions 2006 / 2012): out of stage A here.
      Mpi::Print("\nImpedance-metal rule (RECORD): electrostatic runs are out of stage A: "
                 "every metal is a perfect conductor there (nothing to substitute; no "
                 "candidate, no law reading)\n");
    }
    else if (response_config->impedance_metal_rule ==
             config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceMetalRule::
                 APPLY)
    {
      // The impedance-metal rule's "Apply" (decision 2206 (a); DESIGN-APPLY s8): an
      // electrostatic run solves no impedance boundary, so the rule serves nothing here.
      Mpi::Print("\nImpedance-metal rule (APPLY): electrostatic runs solve no impedance "
                 "boundary: every metal is a perfect conductor there; the rule serves "
                 "nothing; values as \"Off\"\n");
    }
    response_correction = SurfaceResponseOperator::ConstructElectrostatic(
        iodata, laplace_op, &response_geometry);
    // Lane fix-response-memory (decisions 1446 / 1447): the gathered per-patch records live
    // on the root only (their consumers are the root's tables below).
    response_correction->SetPatchRecordResidency(
        SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY);
    // Round 2a R2 (decision 1475 (1)): the root's Records spilled to the output directory.
    response_correction->SpillConductorConsistencyRecords(iodata.problem.output);
    // Round 2b B3 (decision 1484): no AMR construction follows (Refinement.MaxIts 0), so
    // the per-process response-matrix caches serve nothing more.
    if (iodata.model.refinement.max_it == 0)
    {
      ReleaseResponseMatrixCaches();
    }
    if (response_correction->IsRawFallback())
    {
      SaveMetadata(*response_correction);
      response_correction.reset();
      postprocess_response = self_consistent_response = patch_energy_response = false;
    }
    else
    {
      // Decision 2192 (review r2184 MAJOR-1): the census table's path is asked on EVERY
      // rank AHEAD of the root's patch table, so that a failure of the root's directory
      // work ends every rank here with the collective message before the root writes into
      // the directory (PatchEnergy only, the census table's own condition; without it the
      // directory work stays the root writers' own, as before).
      std::filesystem::path resolution_patches_path;
      if (patch_energy_response)
      {
        resolution_patches_path = response_correction->CollectiveDiagnosticsPath(
            post_dir, "surface-response-resolution-patches.csv");
      }
      WritePatchTable();
      if (patch_energy_response)
      {
        // Lane rescensus-product (decision 1817): the census's per-patch readings
        // (PatchEnergy only; collective; its path prepared above).
        response_correction->WriteResolutionPatchTable(resolution_patches_path.string());
      }
      if (self_consistent_response)
      {
        corrected_K = std::make_unique<SumOperator>(*K, *response_correction);
        system_K = corrected_K.get();
      }
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
  SaveMetadata(post_op.GetSurfacePostOp());
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0,
              "No terminal or prescribed potential boundaries specified for electrostatic "
              "simulation!");
  const auto &response_archive = response_archive_settings.directory;
  const bool archive_stream_only = response_archive_settings.archive_only;
  const bool response_source_timing = response_archive_settings.source_timing;
  const bool archive_recycle_initial_guess =
      archive_stream_only && response_archive_settings.recycle_initial_guess;
  MFEM_VERIFY(!archive_stream_only ||
                  (response_archive && !iodata.boundaries.prescribed_potential.empty() &&
                   iodata.solver.electrostatic.response_matrix &&
                   iodata.solver.electrostatic.aggregate_response_matrix),
              "The response-archive mode Archive (PALACE_RESPONSE_ARCHIVE_ONLY) requires a "
              "Directory (PALACE_RESPONSE_ARCHIVE_DIR) and an aggregated "
              "prescribed-potential response-matrix configuration!");
  MFEM_VERIFY(!archive_reduce_only ||
                  (response_archive && !iodata.boundaries.prescribed_potential.empty() &&
                   iodata.solver.electrostatic.response_matrix &&
                   iodata.solver.electrostatic.aggregate_response_matrix),
              "The response-archive mode Reduce (PALACE_RESPONSE_REDUCE_ONLY) requires a "
              "Directory (PALACE_RESPONSE_ARCHIVE_DIR) and an aggregated "
              "prescribed-potential response-matrix configuration!");
  if (response_archive)
  {
    if (root)
    {
      std::filesystem::create_directories(*response_archive);
    }
    Mpi::Barrier(laplace_op.GetComm());
  }

  // Right-hand side term and solution vector storage.
  // D is sized by the flux space of the estimator, below.
  Vector RHS(Grad.Width()), E(Grad.Height()), D;
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
    std::vector<SurfaceResponseOperator::PatchContribution> raw_patch_contributions;
    std::vector<SurfaceResponseOperator::PatchContribution> corrected_patch_contributions;
    // Uncovered requirements (decision 394 F2): the raw within-R interface energy of the
    // unmatched features' portions per target interface, on the raw field (kept in the
    // fixed-trace and fixed-flux energies) and on the corrected field (kept in the
    // self-consistent energy); empty when nothing is uncovered.
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_uncovered;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> corrected_uncovered;
    // DomainBoundary raw claims (F-DB-a, decisions 442 / 454): the raw within-R interface
    // energy of the excluded patches' claims, on the raw field (fixed trace / fixed flux)
    // and on the corrected field (self-consistent); empty when nothing is excluded.
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_domain_boundary;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> corrected_domain_boundary;
    // SpatialClusterPriority raw remainders (OI-SCP-HOLE, decision 922): the raw within-R
    // interface energy of the dropped coupons' windows outside every cluster box, on the
    // raw field and on the corrected field; empty when no coupon is dropped.
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_spatial_cluster_priority;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
        corrected_spatial_cluster_priority;
    // ConductorConsistency raw claims (S7diag R1, decision 758): the raw within-R interface
    // energy of the claims of the patches the conductor-consistency gate has excluded by
    // this excitation (zero entries for an excitation before the first exclusion), on the
    // raw field and on the corrected field; empty when the operator has no probes.
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_conductor_consistency;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> corrected_conductor_consistency;
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
  // The uncovered requirements (decision 394 F2) as the postoperator's perimeter portions.
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> uncovered_portions;
  // The DomainBoundary raw claims (F-DB-a) as the same kind of portions, integrated by the
  // same machinery in a separate call (their own table and share).
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> domain_boundary_portions;
  // Geometry fallbacks (USER decision 1000 (b) / 1001): the converted sites' excluded
  // claims, ALREADY inside domain_boundary_portions (kept raw, F-DB-a); this copy is
  // measured for the record's RawShareByEnergy and surface-response-geometry-fallbacks.csv
  // (written only when a fallback fired).
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> geometry_fallback_portions;
  std::unique_ptr<std::ofstream> geometry_fallback_csv;
  auto MeasureGeometryFallbacks = [&](int source, const EnergyData &energies)
  {
    if (geometry_fallback_portions.empty())
    {
      return;
    }
    BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
    const auto targets = response_correction->GetTargetInterfaces();
    const auto raw_fallback =
        post_op.GetElectrostaticUncoveredEdgeEnergies(targets, geometry_fallback_portions);
    const auto &types = response_correction->GetGeometryFallbackInterfaceTypes();
    std::map<std::string, double> fallback_by_type, within_by_type;
    std::map<std::string, std::map<std::string, double>> fallback_by_class_type;
    for (const auto &[interface, data] : energies.interfaces)
    {
      if (!targets.count(interface) || data.edge_energies.empty())
      {
        continue;
      }
      const auto type = types.find(interface);
      const std::string name = type == types.end() ? "Default" : ToString(type->second);
      within_by_type[name] += data.energy - data.edge_energies.back().energy_outside;
      const auto fallback = raw_fallback.find(interface);
      if (fallback != raw_fallback.end())
      {
        fallback_by_type[name] += fallback->second.energy;
        for (const auto &[cls, energy] : fallback->second.by_type)
        {
          fallback_by_class_type[cls][name] += energy;
        }
      }
    }
    std::map<std::string, double> shares;
    for (const auto &[name, within] : within_by_type)
    {
      shares[name] = within > 0.0 ? fallback_by_type[name] / within : 0.0;
    }
    response_correction->SetGeometryFallbackEnergyShares(shares);
    if (!root)
    {
      return;
    }
    // The entries of a class touching a type (an entry naming no interface counts on every
    // type: the record's rule) and their summed Length in metres (DESIGN 2 (c)).
    std::map<std::string, std::map<std::string, std::pair<int, double>>>
        entries_by_class_type;
    for (const auto &entry : response_correction->GetGeometryFallbacks())
    {
      std::set<std::string> entry_types;
      if (entry.interfaces.empty())
      {
        for (const auto &[name, within] : within_by_type)
        {
          (void)within;
          entry_types.insert(name);
        }
      }
      for (const int interface : entry.interfaces)
      {
        const auto type = types.find(interface);
        entry_types.insert(type == types.end() ? "Default" : ToString(type->second));
      }
      for (const auto &name : entry_types)
      {
        auto &[count, length] = entries_by_class_type[entry.cls][name];
        count++;
        length += entry.length;
      }
    }
    if (!geometry_fallback_csv)
    {
      // OI-ES-EXCLUSION-CSV-SYMLINK-OVERWRITE (fix round P2): a real top-level file on
      // every adaptive cycle (see the exclusion-energy table below).
      std::filesystem::remove(post_dir / "surface-response-geometry-fallbacks.csv");
      geometry_fallback_csv = std::make_unique<std::ofstream>(
          post_dir / "surface-response-geometry-fallbacks.csv");
      *geometry_fallback_csv
          << "source,class,type,count,length (m),raw energy (J),within-R "
             "energy (J),share\n"
          << std::setprecision(17);
    }
    for (const auto &[cls, by_type] : fallback_by_class_type)
    {
      for (const auto &[name, energy] : by_type)
      {
        const double within = within_by_type[name];
        const auto [count, length] = entries_by_class_type[cls][name];
        *geometry_fallback_csv
            << source << "," << cls << "," << name << "," << count << ","
            << iodata.units.Dimensionalize<Units::ValueType::LENGTH>(length) << ","
            << iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy) << ","
            << iodata.units.Dimensionalize<Units::ValueType::ENERGY>(within) << ","
            << (within > 0.0 ? energy / within : 0.0) << "\n";
      }
    }
    geometry_fallback_csv->flush();
  };
  // OI-COVERAGE-EXCLUSION-ENERGY (idfix DESIGN section 8; decisions 1064 (2) / 1074): with
  // Diagnostics.ExclusionEnergy on, the raw within-R energy of the EXCLUDED / unclaimed
  // perimeter portions per class (the identification's classes by name, CornerArmTrim,
  // Unclaimed-other, UnprobedModels, MaxwellPatchExcluded) x interface type x source ->
  // surface-response-exclusion-energy.csv. A diagnostic partition: no energy column of any
  // other table changes; nothing is written when the flag is off.
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> exclusion_energy_portions;
  std::unique_ptr<std::ofstream> exclusion_energy_csv;
  auto MeasureExclusionEnergy = [&](int source, const EnergyData &energies)
  {
    if (!response_correction || !response_correction->ExclusionEnergyRequested())
    {
      return;
    }
    BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
    const auto targets = response_correction->GetTargetInterfaces();
    // The partition against the CLASSIFIED edge tree (decision 1076): every class's energy
    // + "Claimed" (the feature-claimed perimeter's), and the "Total" row = the within-R
    // energy under that tree's distance by ONE integral (the P == S identity: the classes
    // and Claimed sum to it).
    const auto raw_excluded =
        post_op.GetElectrostaticClassifiedEdgeEnergies(targets, exclusion_energy_portions);
    const auto raw_total = post_op.GetElectrostaticClassifiedEdgeEnergies(targets, {});
    if (!root)
    {
      return;
    }
    const auto &types = response_correction->GetGeometryFallbackInterfaceTypes();
    auto TypeName = [&](int interface)
    {
      const auto type = types.find(interface);
      return type == types.end() ? std::string("Default") : ToString(type->second);
    };
    std::map<std::string, double> interface_raw_by_type;
    std::map<std::string, std::map<std::string, double>> excluded_by_class_type;
    for (const auto &[interface, data] : energies.interfaces)
    {
      if (!targets.count(interface))
      {
        continue;
      }
      const std::string name = TypeName(interface);
      interface_raw_by_type[name] += data.energy;
      if (const auto excluded = raw_excluded.find(interface);
          excluded != raw_excluded.end())
      {
        for (const auto &[cls, energy] : excluded->second)
        {
          excluded_by_class_type[cls][name] += energy;
        }
      }
      if (const auto total = raw_total.find(interface); total != raw_total.end())
      {
        for (const auto &[cls, energy] : total->second)
        {
          (void)cls;  // the one remainder group of an empty portion set
          excluded_by_class_type["Total"][name] += energy;
        }
      }
    }
    // The portions of every class (all types) and their summed length in metres; Total =
    // every portion.
    std::map<std::string, std::pair<int, double>> portions_by_class;
    for (const auto &portion : exclusion_energy_portions)
    {
      double length2 = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
      }
      for (const std::string &key : {portion.type, std::string("Total")})
      {
        auto &[count, length] = portions_by_class[key];
        count++;
        length += std::sqrt(length2);
      }
    }
    if (!exclusion_energy_csv)
    {
      // OI-ES-EXCLUSION-CSV-SYMLINK-OVERWRITE (fix round P2): on a later adaptive cycle the
      // path is the symlink SaveIteration left into the previous cycle's archive folder;
      // remove it first so this cycle's table is a real top-level file (the F-18 SaveMesh
      // form) and the archived one keeps its cycle's content.
      std::filesystem::remove(post_dir / "surface-response-exclusion-energy.csv");
      exclusion_energy_csv = std::make_unique<std::ofstream>(
          post_dir / "surface-response-exclusion-energy.csv");
      *exclusion_energy_csv << "source,class,type,class portions,class length (m),raw "
                               "energy (J),interface raw energy (J),share\n"
                            << std::setprecision(17);
    }
    for (const auto &[cls, by_type] : excluded_by_class_type)
    {
      const auto [count, length] = portions_by_class[cls];
      for (const auto &[name, energy] : by_type)
      {
        const double interface_raw = interface_raw_by_type[name];
        *exclusion_energy_csv
            << source << "," << cls << "," << name << "," << count << ","
            << iodata.units.Dimensionalize<Units::ValueType::LENGTH>(length) << ","
            << iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy) << ","
            << iodata.units.Dimensionalize<Units::ValueType::ENERGY>(interface_raw) << ","
            << (interface_raw > 0.0 ? energy / interface_raw : 0.0) << "\n";
      }
    }
    exclusion_energy_csv->flush();
  };
  // The SpatialClusterPriority raw remainders (decision 922) likewise.
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
      spatial_cluster_priority_portions;
  // Ruling 2276: the dropped coupons' windows for the coverage partition (null: none).
  const SurfacePostOperator::CoverageDroppedWindows *dropped_windows = nullptr;
  // PatchEnergy (decision 922 B-1; the 851 (C) export): the applied patches' raw claims by
  // model set (empty without PatchEnergy) and, per source, the raw field's within-R energy
  // partitioned over them.
  SurfacePostOperator::CompactPerimeterPortions model_raw_portions;
  // Ruling 2251 (OI-REGION-TARGET-CORRECTED-NOT-FILTERED): the targets with a Region; their
  // corrected columns count each patch's share inside the region.
  std::set<int> region_targets;
  std::vector<std::pair<int, std::map<int, std::map<std::string, double>>>>
      raw_model_partitions;
  // Lane L1 (decision 1634 (1)): per source, the closed coverage partition of its raw
  // field (SurfacePostOperator::GetInterfaceCoveragePartition), on every corrected run.
  std::vector<
      std::pair<int, std::map<int, SurfacePostOperator::InterfaceCoveragePartition>>>
      coverage_partitions;
  if (response_correction)
  {
    for (const auto &portion : response_correction->GetUncoveredPortions())
    {
      uncovered_portions.push_back(
          {portion.p0, portion.p1, portion.topology, portion.feature});
    }
    for (const auto &portion : response_correction->GetDomainBoundaryPortions())
    {
      domain_boundary_portions.push_back(
          {portion.p0, portion.p1, portion.topology, portion.feature});
    }
    for (const auto &portion : response_correction->GetGeometryFallbackPortions())
    {
      geometry_fallback_portions.push_back(
          {portion.p0, portion.p1, portion.topology, portion.feature});
    }
    for (const auto &portion : response_correction->GetSpatialClusterPriorityPortions())
    {
      spatial_cluster_priority_portions.push_back(
          {portion.p0, portion.p1, portion.topology, portion.feature});
    }
    if (const auto &windows =
            response_correction->GetSpatialClusterPriorityExclusions().coverage_windows;
        !windows.Empty())
    {
      dropped_windows = &windows;
    }
    // Round 2b B1 (decision 1484): moved out of the operator (the same portions in the
    // same order; the operator reads them nowhere after the construction), held once, in
    // the compact exact form (lane L1, ruling 1682).
    model_raw_portions = response_correction->ReleaseModelRawClaimPortionList();
    if (response_correction->HasRegionClaims())
    {
      region_targets =
          post_op.GetRegionInterfaces(response_correction->GetTargetInterfaces());
    }
    // Ruling 1688 (lane L1): the stage-0 Missing LENGTH line (root; log only).
    if (!response_correction->IsRawFallback())
    {
      const auto [missing, method_gap] = SplitCoverageUncoveredByKeyClass(
          uncovered_portions, response_correction->GetCoverageKeyRecords());
      const auto targets = response_correction->GetTargetInterfaces();
      Mpi::Print("\n{}",
                 CoverageStageZeroLine(
                     post_op.GetInterfaceTypes(), targets,
                     post_op.GetSurfacePostOp().GetInterfacePerimeterLengths(
                         targets, {missing, method_gap}),
                     response_correction->GetIdentificationFeatureLengths(), iodata.units));
      const auto not_served = response_correction->GetConfiguredTargetsNotServed();
      if (!not_served.empty())
      {
        Mpi::Warning(
            "{}", CoverageTargetsNotServedWarning(post_op.GetInterfaceTypes(&not_served)));
      }
      // Lane rescensus-product (decision 1817; ruling 1819 F1): the census Types (the
      // coverage interfaces') and the PREFLIGHT line (collective; before the first solve).
      const auto census_interfaces = response_correction->GetCoverageInterfaces();
      response_correction->ConfigureResolutionTypes(
          post_op.GetInterfaceTypes(&census_interfaces));
    }
    if (response_correction->ExclusionEnergyRequested())
    {
      for (const auto &portion : response_correction->GetExclusionEnergyPortions())
      {
        exclusion_energy_portions.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
    }
  }
  // The (ii-a) clipped-correction DIAGNOSTIC of the SpatialClusterPriority-dropped coupons
  // (OI-SCP-BOXGROW, decisions 958 (b) / 960 / 969): per dropped patch, its whole window
  // and its outside remainder as portions labelled "p<patch>" (the raw within-R energy per
  // label = RawWhole / RawOut), paired per source and evaluation with the coupon's
  // fabricated energy on the device trace (Fab, the operator's diagnostic-only patches) in
  // surface-response-spatial-cluster-priority-clipped-correction.csv. DERIVED, never in a
  // total; empty lists when no coupon is dropped (nothing evaluated).
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> clipped_window_portions;
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> clipped_outside_portions;
  struct ClippedCorrectionRow
  {
    int source = 0;
    int evaluation = 0;
    SurfaceResponseOperator::ClippedCorrectionDiagnostic diagnostic;
    std::map<int, double> raw_whole, raw_outside;  // per interface
    bool fixed_flux = false;
  };
  std::vector<ClippedCorrectionRow> clipped_correction_rows;
  if (response_correction)
  {
    for (const auto &entry :
         response_correction->GetSpatialClusterPriorityExclusions().entries)
    {
      for (const auto &portion : entry.window_portions)
      {
        clipped_window_portions.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
      for (const auto &portion : entry.outside_portions)
      {
        clipped_outside_portions.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
    }
  }
  // The per-patch raw energies of the current field (the postoperator's grid functions):
  // RawWhole / RawOut per label, keyed by the dropped patch index.
  auto ClippedRawEnergies = [&](const auto &portions)
  {
    std::map<int, std::map<int, double>> by_patch;  // patch -> interface -> energy
    if (portions.empty())
    {
      return by_patch;
    }
    BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
    const auto energies = post_op.GetElectrostaticUncoveredEdgeEnergies(
        response_correction->GetTargetInterfaces(), portions);
    for (const auto &[interface, data] : energies)
    {
      for (const auto &[label, value] : data.by_type)
      {
        if (label.size() > 1 && label[0] == 'p')
        {
          by_patch[std::stoi(label.substr(1))][interface] = value;
        }
      }
    }
    return by_patch;
  };
  auto AppendClippedRows = [&](int source, const Vector &field, bool include_fixed_flux,
                               const std::map<int, std::map<int, double>> &raw_whole,
                               const std::map<int, std::map<int, double>> &raw_outside)
  {
    const auto diagnostics =
        response_correction->GetSpatialClusterPriorityClippedCorrection(field,
                                                                        include_fixed_flux);
    for (const auto &diagnostic : diagnostics)
    {
      ClippedCorrectionRow row;
      row.source = source;
      row.evaluation = include_fixed_flux ? 0 : 2;
      row.diagnostic = diagnostic;
      const auto whole = raw_whole.find(diagnostic.patch);
      const auto outside = raw_outside.find(diagnostic.patch);
      if (whole != raw_whole.end())
      {
        row.raw_whole = whole->second;
      }
      if (outside != raw_outside.end())
      {
        row.raw_outside = outside->second;
      }
      clipped_correction_rows.push_back(row);
      if (include_fixed_flux)
      {
        row.evaluation = 1;
        row.fixed_flux = true;
        clipped_correction_rows.push_back(row);
      }
    }
  };
  // The corrected capacitance's self-consistent column (PostprocessCorrectedTerminals) is
  // written only when every source's corrected solve was accepted (fail closed).
  bool self_consistent_accepted = self_consistent_response;
  // The spatial coupons' traces and the per-patch energies (PatchEnergy) are written per
  // source as they are gathered (lane fix-response-memory C4, decisions 1439 / 1446: the
  // gathered O(global patches) vectors are not retained across the sources; the files'
  // bytes are the former end-of-run whole-table writes',
  // TableWithCSVFile::AppendRowsAndClear). Root only: the records are gathered to the root
  // (PatchRecordResidency::ROOT_ONLY).
  std::unique_ptr<TableWithCSVFile> traces_output;
  auto WriteSpatialPatchTraces =
      [&](int source, const std::vector<SurfaceResponseOperator::PatchTrace> &entries)
  {
    if (!root || entries.empty())
    {
      return;
    }
    using VT = Units::ValueType;
    if (!traces_output)
    {
      traces_output = std::make_unique<TableWithCSVFile>(
          response_correction->DiagnosticsPath(post_dir, "surface-response-traces.csv"));
      traces_output->table.insert(Column::Integer("source", "i"));
      traces_output->table.insert(Column::Integer("patch", "patch"));
      traces_output->table.insert(Column::Integer("model", "model"));
      traces_output->table.insert(Column::Integer("coefficient", "coefficient"));
      traces_output->table.insert(Column::Integer("conductor_state", "conductor state"));
      traces_output->table.insert("value", "value (V)");
    }
    auto &traces = *traces_output;
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
        traces.table["value"] << iodata.units.Dimensionalize<VT::VOLTAGE>(
            entry.coefficients[i]);
      }
    }
    traces.AppendRowsAndClear();
  };
  std::unique_ptr<TableWithCSVFile> patch_output;
  std::vector<int> patch_output_interfaces;
  std::map<int, const SurfaceResponseOperator::PatchAssignment *> patch_assignments;
  auto WritePatchEnergyRows = [&](CorrectedResult &result)
  {
    // Decision 352 follow-up (1): one row per applied patch, source and evaluation (the
    // codes of surface-response-model-energy.csv); the rows of one (source, evaluation,
    // model) sum to that model's row. `patch` is the 1-based placed patch index of
    // surface-response-traces.csv (the preflight dry run's 0-based `Patch` + 1); the
    // origin and the longitudinal cell (offsets along the patch's AxisW) locate it. The
    // interface columns are the first source's.
    if (!root || !patch_energy_response)
    {
      return;
    }
    using VT = Units::ValueType;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!patch_output)
    {
      for (const auto &assignment : response_correction->GetPatchAssignments())
      {
        patch_assignments.emplace(assignment.global_index, &assignment);
      }
      patch_output =
          std::make_unique<TableWithCSVFile>(response_correction->DiagnosticsPath(
              post_dir, "surface-response-patch-energy.csv"));
      patch_output->table.insert(Column::Integer("source", "source"));
      patch_output->table.insert(Column::Integer("evaluation", "evaluation"));
      patch_output->table.insert(Column::Integer("patch", "patch"));
      patch_output->table.insert(Column::Integer("model", "model"));
      patch_output->table.insert(Column::Integer("feature", "feature"));
      patch_output->table.insert("weight", "weight");
      patch_output->table.insert("cell_begin", "cell begin (m)");
      patch_output->table.insert("cell_end", "cell end (m)");
      for (const char *coordinate : {"x", "y", "z"})
      {
        patch_output->table.insert(fmt::format("origin_{}", coordinate),
                                   fmt::format("origin {} (m)", coordinate));
      }
      patch_output->table.insert("domain_correction", "domain correction (J)");
      for (const auto &[interface, data] : result.raw.interfaces)
      {
        (void)data;
        patch_output_interfaces.push_back(interface);
        patch_output->table.insert(
            fmt::format("interface_{}", interface),
            fmt::format("fabricated surface energy[{}] (J)", interface));
      }
    }
    auto AppendPatch = [&](int source, int evaluation,
                           const SurfaceResponseOperator::PatchContribution &data,
                           bool available, bool fixed_flux)
    {
      const auto assignment = patch_assignments.find(data.patch);
      MFEM_VERIFY(assignment != patch_assignments.end(),
                  "Per-patch response energy of an unknown patch " << data.patch + 1
                                                                   << "!");
      patch_output->table["source"] << source;
      patch_output->table["evaluation"] << evaluation;
      patch_output->table["patch"] << data.patch + 1;
      patch_output->table["model"] << data.model;
      patch_output->table["feature"] << data.feature;
      patch_output->table["weight"] << data.weight;
      patch_output->table["cell_begin"]
          << iodata.units.Dimensionalize<VT::LENGTH>(data.cell[0]);
      patch_output->table["cell_end"]
          << iodata.units.Dimensionalize<VT::LENGTH>(data.cell[1]);
      for (int d = 0; d < 3; d++)
      {
        patch_output->table[fmt::format("origin_{}", "xyz"[d])]
            << iodata.units.Dimensionalize<VT::LENGTH>(assignment->second->origin[d]);
      }
      patch_output->table["domain_correction"]
          << (available ? iodata.units.Dimensionalize<VT::ENERGY>(
                              fixed_flux ? data.domain_correction_fixed_flux
                                         : data.domain_correction)
                        : nan);
      for (const int interface : patch_output_interfaces)
      {
        const auto &energies = fixed_flux ? data.fabricated_surface_energy_fixed_flux
                                          : data.fabricated_surface_energy;
        const auto energy = energies.find(interface);
        patch_output->table[fmt::format("interface_{}", interface)]
            << (available ? iodata.units.Dimensionalize<VT::ENERGY>(
                                energy != energies.end() ? energy->second : 0.0)
                          : nan);
      }
    };
    for (const auto &contribution : result.raw_patch_contributions)
    {
      AppendPatch(result.source, 0, contribution, result.has_postprocessed, false);
      AppendPatch(result.source, 1, contribution, result.has_postprocessed, true);
    }
    for (const auto &contribution : result.corrected_patch_contributions)
    {
      AppendPatch(result.source, 2, contribution, result.has_self_consistent, false);
    }
    patch_output->AppendRowsAndClear();
  };
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
        iodata.solver.linear.estimator_mg, iodata.boundaries.GetSheetAttributes());
  }
  if (estimator)
  {
    // The recovered electric flux lives on the estimator's recovery space (broken across
    // the sheet attributes when the mesh has interior faces there).
    D.SetSize(estimator->GetFluxSpace().GetTrueVSize());
    D.UseDevice(true);
  }
  ErrorIndicator indicator;

  if (archive_reduce_only)
  {
    Mpi::Print("\nReducing archived electrostatic response fields from {}\n",
               response_archive->string());
    PostprocessArchivedResponseMatrix(post_op, laplace_op, Grad, *response_archive,
                                      response_archive_settings.block_size);
    SaveLinearSolverMetadata(laplace_op.GetComm(), 0, 0);
    return {indicator, laplace_op.GlobalTrueVSize()};
  }

  // Lane L1 (decision 1634 (2) / (4); rulings 1641 / 1654 (III)): per source its coverage
  // (root; the cells per Type of its raw field), in source order, and the number of
  // conductor-consistency portions it was measured with (the matrix coverage is re-measured
  // with the FINAL portions: ruling 1654 (III) Q3). The coverage tables are written as each
  // source completes (a killed run keeps the finished sources).
  std::vector<std::pair<int, CoverageMeasurement>> coverage_measurements;
  std::vector<std::size_t> coverage_conductor_consistency_counts;
  // Lane ccpx (Exclusion PerExcitation): per source (the corrected_results order) the
  // patches excluded at its own gate pass and its own portion record: a source is
  // re-measured for the matrix where its exclusion SET differs from the union, and the
  // conductor-consistency table's rows carry their own source's portion columns.
  std::vector<std::vector<int>> coverage_conductor_consistency_patches;
  std::vector<ExcludedPatchPortions> conductor_consistency_source_records;
  const auto coverage_path = post_dir / "surface-response-coverage.csv";
  const auto coverage_by_key_path =
      post_dir / "surface-response-uncovered-energy-by-key.csv";
  // Ruling 2277: the near-edge energy budget beside them (top level).
  const auto budget_path = post_dir / kNearEdgeBudgetTableName;
  // Lane rescensus-product (decision 1817): the census table (top level: small, one row per
  // source x Type x (planes + 1)), written as each source completes.
  const auto resolution_path = post_dir / "surface-response-resolution.csv";
  if (root && response_correction)
  {
    fs::remove(coverage_path);
    fs::remove(coverage_by_key_path);
    fs::remove(budget_path);
    fs::remove(post_dir / "terminal-C-corrected-coverage.csv");
    fs::remove(resolution_path);
  }
  auto MeasureSourceCoverage =
      [&](const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &partition,
          const EnergyData &raw, const std::map<int, double> *fabricated)
  {
    // W per target = the raw interface energy minus its energy outside R (the largest
    // configured radius, the coupon's: the corrected energies' own split); NOT-MEASURED
    // without a postprocessed response (ruling 1654 (III) Q2 (b): a raw response evaluation
    // of a SelfConsistent-only run would move the operator's trace counters).
    // Ruling 1738: the configured targets no loaded model serves are cells' interfaces too
    // (W = raw - outside R, nothing fabricated: uncovered).
    const auto targets = response_correction->GetCoverageInterfaces();
    const auto not_served = response_correction->GetConfiguredTargetsNotServed();
    const auto types = post_op.GetInterfaceTypes(&targets);
    std::map<int, double> within;
    for (const int target : targets)
    {
      const auto interface = raw.interfaces.find(target);
      if (interface == raw.interfaces.end() || interface->second.edge_energies.empty())
      {
        fabricated = nullptr;
        break;
      }
      within[target] =
          interface->second.energy - interface->second.edge_energies.back().energy_outside;
    }
    auto measurement = fabricated
                           ? ComputeCoverageMeasurement(
                                 partition, types, within, *fabricated,
                                 response_correction->GetCoverageKeyRecords(),
                                 response_correction->GetCoverageModelQualifications(),
                                 response_correction->GetCoverageThresholds(),
                                 response_correction->GetCoverageValidatedRegime(),
                                 CoverageExcitationClass::ELECTROSTATIC_TERMINAL)
                           : NotMeasuredCoverageMeasurement(
                                 types, CoverageExcitationClass::ELECTROSTATIC_TERMINAL);
    MarkCoverageTargetsNotServed(measurement, post_op.GetInterfaceTypes(&not_served));
    // Ruling 1694 (1): the postprocessing interfaces no model targets carry their raw
    // energy; their Types' untargeted raw energy and share (NOT-TARGET cells for a Type
    // with untargeted interfaces only).
    std::map<int, std::string> untargeted;
    for (const auto &[interface, type] : post_op.GetInterfaceTypes())
    {
      if (targets.count(interface) == 0)
      {
        untargeted[interface] = type;
      }
    }
    std::map<int, double> raw_energy;
    for (const auto &[interface, data] : raw.interfaces)
    {
      raw_energy[interface] = data.energy;
    }
    AddCoverageUntargetedInterfaces(measurement, untargeted, raw_energy);
    return measurement;
  };
  auto WriteSourceCoverage = [&](int source, const CoverageMeasurement &measurement)
  {
    const auto targets = response_correction->GetCoverageInterfaces();
    const bool first = coverage_measurements.size() == 1;
    for (const auto &path : {coverage_path, coverage_by_key_path, budget_path})
    {
      std::ofstream output(path, first ? std::ios::out | std::ios::trunc
                                       : std::ios::out | std::ios::app);
      MFEM_VERIFY(output, "Unable to open " << path.string() << " for writing!");
      if (path == coverage_path)
      {
        output << (first ? CoverageTableHeader("i") : "");
        WriteCoverageRows(output, std::to_string(source), measurement, iodata.units);
      }
      else if (path == budget_path)
      {
        output << (first ? NearEdgeBudgetTableHeader("i") : "");
        WriteNearEdgeBudgetRows(output, std::to_string(source), measurement, iodata.units);
      }
      else
      {
        output << (first ? CoverageByKeyTableHeader("i", targets) : "");
        WriteCoverageByKeyRows(output, std::to_string(source), measurement, targets,
                               CoveragePortionsByFeature(uncovered_portions), iodata.units);
      }
    }
    response_correction->GetCoverageRunRecord().Add({{"i", static_cast<double>(source)}},
                                                    measurement);
  };
  auto WriteMatrixCoverage = [&](const std::map<int, mfem::Array<int>> &terminal_sources,
                                 const std::vector<CoverageMeasurement> &final_coverage,
                                 const std::vector<bool> &remeasured)
  {
    // terminal-C-corrected-coverage.csv: the matrix shape of terminal-C-corrected.csv
    // (rows i, a column per j); entry (i, j) = the worse coverage status over the Types of
    // sources i and j (CoverageStatusRank) on their FINAL measurement; then, per row,
    // source i's final measurement (its worst status and largest uncovered share over
    // Types, and whether a conductor-consistency exclusion after it made it differ from its
    // own row of surface-response-coverage.csv). A text table, root only.
    MFEM_VERIFY(final_coverage.size() == terminal_sources.size() &&
                    remeasured.size() == terminal_sources.size(),
                "The corrected-matrix coverage needs one measurement per terminal source!");
    auto Worst = [](const CoverageMeasurement &measurement)
    {
      const auto *cell = measurement.WorstCoverageCell();
      return cell ? cell->status : CoverageStatus::NOT_MEASURED;
    };
    auto WorstShare = [](const CoverageMeasurement &measurement)
    {
      double worst = std::numeric_limits<double>::quiet_NaN();
      for (const auto &cell : measurement.cells)
      {
        if (std::isfinite(cell.uncovered_share) && !(cell.uncovered_share <= worst))
        {
          worst = cell.uncovered_share;
        }
      }
      return CsvNumber(worst);
    };
    std::ofstream output(post_dir / "terminal-C-corrected-coverage.csv");
    MFEM_VERIFY(output, "Unable to open terminal-C-corrected-coverage.csv for writing!");
    output << "i";
    for (const auto &[idx2, data2] : terminal_sources)
    {
      (void)data2;
      output << fmt::format(",coverage status[i][{}]", idx2);
    }
    output << ",source final worst coverage status,source final worst uncovered raw share "
              "of W,source re-measured\n";
    std::size_t i = 0;
    for (const auto &[idx, data] : terminal_sources)
    {
      (void)data;
      output << idx;
      for (std::size_t j = 0; j < final_coverage.size(); j++)
      {
        const auto a = Worst(final_coverage[i]), b = Worst(final_coverage[j]);
        output << ','
               << CoverageStatusName(CoverageStatusRank(a) >= CoverageStatusRank(b) ? a
                                                                                    : b);
      }
      output << ',' << CoverageStatusName(Worst(final_coverage[i])) << ','
             << WorstShare(final_coverage[i]) << ',' << (remeasured[i] ? 1 : 0) << '\n';
      i++;
    }
  };

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
      // The initial guess (the essential values; the corrected right-hand side is formed
      // after the conductor-consistency gate below, on the patches it leaves applied).
      V_corrected[step] = V[step];
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
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_uncovered;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_domain_boundary;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_spatial_cluster_priority;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> raw_conductor_consistency;
    // The ConductorConsistency raw claims (S7diag R1, decision 758) as the same kind of
    // portions: the list the gate has built so far, read AFTER this excitation's gate (the
    // exclusion is sticky, so the list only grows; empty on a window with 0 exclusions).
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
        conductor_consistency_portions;
    // OI-ES-COMPARATOR-NC-POSTRAW-STALL (decision 1114): the wall time of the response
    // stages after the raw energies, one line each (a chip-scale run prints nothing between
    // "Field energy" and the corrected solve otherwise); output text only.
    const auto post_raw_started = std::chrono::steady_clock::now();
    auto PostRawStage = [&](const std::string &text)
    {
      Mpi::Print(
          " Response post-raw {} ({:.2f} s)\n", text,
          std::chrono::duration<double>(std::chrono::steady_clock::now() - post_raw_started)
              .count());
    };
    if (response_correction)
    {
      raw_energies = post_op.GetCachedElectrostaticEnergies();
      if (!uncovered_portions.empty())
      {
        // The uncovered requirements' raw within-R energy (decision 394 F2) on the raw
        // field the grid functions hold after MeasureAndPrintAll.
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        raw_uncovered = post_op.GetElectrostaticUncoveredEdgeEnergies(
            response_correction->GetTargetInterfaces(), uncovered_portions);
      }
      if (!domain_boundary_portions.empty())
      {
        // The DomainBoundary raw claims' within-R energy (F-DB-a) on the same raw field.
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        raw_domain_boundary = post_op.GetElectrostaticUncoveredEdgeEnergies(
            response_correction->GetTargetInterfaces(), domain_boundary_portions);
      }
      // Geometry fallbacks (USER decision 1000 (b)): the share of the within-R energy the
      // converted sites left raw, on the same raw field (nothing when none fired).
      PostRawStage("uncovered / domain-boundary raw energies");
      MeasureGeometryFallbacks(idx, raw_energies);
      PostRawStage("geometry fallbacks");
      MeasureExclusionEnergy(idx, raw_energies);
      PostRawStage("exclusion energy");
      if (!spatial_cluster_priority_portions.empty())
      {
        // The SpatialClusterPriority raw remainders' within-R energy (decision 922) on the
        // same raw field.
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        raw_spatial_cluster_priority = post_op.GetElectrostaticUncoveredEdgeEnergies(
            response_correction->GetTargetInterfaces(), spatial_cluster_priority_portions);
      }
      if (!clipped_window_portions.empty())
      {
        // The (ii-a) diagnostic on the raw field (evaluations 0 / 1): the per-patch raw
        // energies, then the dropped coupons' fabricated energies on the raw trace.
        const auto raw_whole = ClippedRawEnergies(clipped_window_portions);
        const auto raw_outside = ClippedRawEnergies(clipped_outside_portions);
        BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
        AppendClippedRows(idx, V[step], true, raw_whole, raw_outside);
      }
      if (patch_energy_response && !model_raw_portions.Empty())
      {
        // PatchEnergy (decision 922 B-1): the raw field's within-R energy per model set
        // (the claims exist on every run since lane L1; the table with PatchEnergy only).
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        raw_model_partitions.emplace_back(
            idx, post_op.GetElectrostaticPartitionedEdgeEnergies(
                     response_correction->GetTargetInterfaces(), model_raw_portions));
      }
      PostRawStage("spatial-cluster-priority / clipped / model raw energies");
    }

    if (response_correction)
    {
      BlockTimer response_timer(Timer::POSTPRO_RESPONSE);
      // The conductor-consistency gate (decision 277): the device potential at the spatial
      // coupons' metal cross-sections on their box faces against the conductor potentials,
      // before any energy of this excitation is formed; an excluded patch leaves this
      // excitation's evaluations (fixed trace, self-consistent) under Exclusion
      // PerExcitation, every later one too under Sticky (lane ccpx); the patch table ends
      // at the union on both rules.
      if (response_correction->HasConductorConsistencyProbes())
      {
        {
          BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
          const auto records =
              response_correction->ApplyConductorConsistencyGate(V[step], idx);
          if (std::any_of(records.begin(), records.end(),
                          [](const auto &record) { return record.excluded; }))
          {
            WritePatchTable();
          }
        }
        // S7diag R1 (decision 758): the excluded patches' raw claims are served RAW. Their
        // within-R energy on this excitation's raw field is integrated by the uncovered
        // requirements' machinery (an empty list yields zero entries without an integral:
        // an excitation before the first exclusion kept the model).
        for (const auto &portion : response_correction->GetConductorConsistencyPortions())
        {
          conductor_consistency_portions.push_back(
              {portion.p0, portion.p1, portion.topology, portion.feature});
        }
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        raw_conductor_consistency = post_op.GetElectrostaticUncoveredEdgeEnergies(
            response_correction->GetTargetInterfaces(), conductor_consistency_portions);
        PostRawStage("conductor-consistency gate + raw energies");
      }
      if (response_correction->HasSurfaceResponse())
      {
        // Lane L1 (decision 1634 (1); rulings 1641 / 1649 Q3): the closed coverage
        // partition of this source's RAW field (the grid functions still hold it), every
        // point in ONE bin by the correction's ownership: the classes the corrected
        // interface energies keep raw (ApplyResponse below) first, whatever claim lies
        // under them (a conductor-consistency-excluded patch's claim included), then the
        // served claims, then the remainder labels. Measured on every corrected run.
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        const auto missing_by_feature =
            SurfacePostOperator::CoverageMissingByFeature(uncovered_portions);
        const std::vector<SurfacePostOperator::CoveragePartitionClass> kept_classes = {
            {"Missing", &missing_by_feature},
            {"DomainBoundary", &domain_boundary_portions},
            {"SpatialClusterPriority", &spatial_cluster_priority_portions},
            {"ConductorConsistency", &conductor_consistency_portions}};
        std::vector<SurfacePostOperator::CoveragePartitionClass> remainder_classes;
        if (response_correction->ExclusionEnergyRequested())
        {
          remainder_classes.push_back({"Excluded", &exclusion_energy_portions});
        }
        // Ruling 2295 (A): a box whose patch this source's gate pass excludes (the portion
        // record in force) holds nothing.
        auto partition = post_op.GetElectrostaticCoveragePartition(
            response_correction->GetCoverageInterfaces(), kept_classes, model_raw_portions,
            remainder_classes, dropped_windows,
            response_correction->GetConductorConsistencyPortionRecord().PatchSet());
        SurfacePostOperator::MoveCoverageServedLabel(partition, "unpatched",
                                                     "UnpatchedPairs");
        coverage_partitions.emplace_back(idx, std::move(partition));
      }
      if (self_consistent_response)
      {
        corrected_rhs = RHS;
        response_correction->EliminateRHS(V_corrected[step], corrected_rhs);
      }
      SurfaceResponseOperator::ElectrostaticResponse response;
      // Lane rescensus-product: this source's census sums (evaluation 0: the postprocessed
      // response of the raw field).
      ResolutionAccumulator resolution_sums =
          response_correction->NewResolutionAccumulator();
      // Ruling 2251: the shares of the targets with a Region, from the footprints of the
      // field the post operator holds (the kept classes of the coverage partition first);
      // every row of surface-response-region-share.csv with PatchEnergy.
      auto RegionShares = [&](int evaluation)
      {
        const auto missing_by_feature =
            SurfacePostOperator::CoverageMissingByFeature(uncovered_portions);
        const std::vector<SurfacePostOperator::CoveragePartitionClass> kept_classes = {
            {"Missing", &missing_by_feature},
            {"DomainBoundary", &domain_boundary_portions},
            {"SpatialClusterPriority", &spatial_cluster_priority_portions},
            {"ConductorConsistency", &conductor_consistency_portions}};
        const auto footprints = post_op.GetElectrostaticRegionPatchFootprints(
            region_targets, kept_classes, response_correction->GetRegionClaimPortions(),
            response_correction->GetRegionClaimPatches(),
            response_correction->GetRegionReferencePoints());
        auto shares = SurfaceResponseOperator::ComputeRegionPatchShares(footprints);
        if (patch_energy_response)
        {
          response_correction->WriteRegionShareRows(post_dir, idx, evaluation, footprints,
                                                    shares, iodata.units);
        }
        return shares;
      };
      // Ruling 2251: the fabricated energies the corrected columns and the coverage read:
      // a Region target's box value in place of its window-wide sum.
      auto RegionFabricated =
          [](std::map<int, double> fabricated,
             const std::map<int, SurfaceResponseOperator::RegionTargetResponse> &targets,
             bool fixed_flux)
      {
        for (const auto &[target, box] : targets)
        {
          const auto it = fabricated.find(target);
          if (it != fabricated.end())
          {
            it->second = fixed_flux ? box.fabricated_surface_energy_fixed_flux
                                    : box.fabricated_surface_energy;
          }
        }
        return fabricated;
      };
      if (postprocess_response)
      {
        BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
        if (region_targets.empty())
        {
          response = response_correction->GetElectrostaticResponse(
              V[step], true, patch_energy_response, &resolution_sums);
        }
        else
        {
          const auto region_shares = RegionShares(0);
          response = response_correction->GetElectrostaticResponse(
              V[step], true, patch_energy_response, &resolution_sums, &region_shares);
          response_correction->RecordRegionResponse(
              idx, "Raw", response.region_targets, response.fabricated_surface_energy,
              response.fabricated_surface_energy_fixed_flux,
              iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0));
        }
        WriteSpatialPatchTraces(idx, response_correction->GetSpatialPatchTraces(V[step]));
      }
      const auto region_fabricated = RegionFabricated(response.fabricated_surface_energy,
                                                      response.region_targets, false);
      const auto region_fabricated_fixed_flux = RegionFabricated(
          response.fabricated_surface_energy_fixed_flux, response.region_targets, true);
      auto ApplyResponse =
          [&](EnergyData energies, double domain_correction,
              const std::map<int, double> &fabricated_surface,
              const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> &uncovered,
              const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
                  &domain_boundary,
              const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
                  &spatial_cluster_priority,
              const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
                  &conductor_consistency)
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
          // distance of the coupon response model. An uncovered requirement (decision 394
          // F2) keeps its raw within-R energy of the same field.
          it->second.energy = it->second.edge_energies.back().energy_outside + energy;
          const auto uncovered_energy = uncovered.find(interface);
          if (uncovered_energy != uncovered.end())
          {
            it->second.energy += uncovered_energy->second.energy;
          }
          // A DomainBoundary-excluded patch's raw claim (F-DB-a, decisions 442 / 454)
          // keeps its raw within-R energy of the same field: never dropped.
          const auto domain_boundary_energy = domain_boundary.find(interface);
          if (domain_boundary_energy != domain_boundary.end())
          {
            it->second.energy += domain_boundary_energy->second.energy;
          }
          // A SpatialClusterPriority-dropped coupon's window outside every cluster box
          // (OI-SCP-HOLE, decision 922) keeps its raw within-R energy of the same field:
          // counted, uncorrected (the inside part is the owner coupon's).
          const auto spatial_cluster_priority_energy =
              spatial_cluster_priority.find(interface);
          if (spatial_cluster_priority_energy != spatial_cluster_priority.end())
          {
            it->second.energy += spatial_cluster_priority_energy->second.energy;
          }
          // A conductor-consistency-excluded patch's raw claim (S7diag R1, decision 758)
          // keeps its raw within-R energy of the same field: never dropped.
          const auto conductor_consistency_energy = conductor_consistency.find(interface);
          if (conductor_consistency_energy != conductor_consistency.end())
          {
            it->second.energy += conductor_consistency_energy->second.energy;
          }
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
          postprocess_response
              ? ApplyResponse(raw_energies, response.domain_correction, region_fabricated,
                              raw_uncovered, raw_domain_boundary,
                              raw_spatial_cluster_priority, raw_conductor_consistency)
              : Unavailable(raw_energies);
      auto postprocessed_fixed_flux =
          postprocess_response
              ? ApplyResponse(raw_energies, response.domain_correction_fixed_flux,
                              region_fabricated_fixed_flux, raw_uncovered,
                              raw_domain_boundary, raw_spatial_cluster_priority,
                              raw_conductor_consistency)
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
      std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> corrected_uncovered;
      std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> corrected_domain_boundary;
      std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
          corrected_spatial_cluster_priority;
      std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
          corrected_conductor_consistency;
      std::vector<SurfaceResponseOperator::ModelContribution> corrected_contributions;
      std::vector<SurfaceResponseOperator::PatchContribution> corrected_patch_contributions;
      if (self_consistent_response)
      {
        Mpi::Print(" Solving fabrication-response corrected field\n");
        V_corrected[step] = V[step];
        const auto solves_before = ksp.NumTotalMult();
        const auto iterations_before = ksp.NumTotalMultIterations();
        const auto corrected_record = SolveCorrectedField(
            ksp, *K, *system_K, response_config->solve_tol,
            iodata.solver.linear.initial_guess, corrected_rhs, V_corrected[step]);
        PostRawStage("corrected solve");
        corrected_linear_solves += ksp.NumTotalMult() - solves_before;
        corrected_linear_iterations += ksp.NumTotalMultIterations() - iterations_before;
        self_consistent_accepted = self_consistent_accepted && corrected_record.accepted;
        if (!corrected_record.converged)
        {
          Mpi::Warning(
              "Self-consistent response-corrected solve did not converge (relative "
              "residual = {:.3e}); corrected energies and participations are reported "
              "as unavailable instead of evaluating the unconverged field.\n",
              corrected_record.relative_residual);
        }
        else if (!corrected_record.accepted)
        {
          const std::string ritz_range =
              corrected_record.iterations > corrected_solve_ritz_max_iterations
                  ? fmt::format("Ritz values of the preconditioned operator not evaluated: "
                                "{:d} iterations exceed the {:d}-iteration bound of the "
                                "Lanczos tridiagonal eigenvalue solve",
                                corrected_record.iterations,
                                corrected_solve_ritz_max_iterations)
                  : fmt::format("Ritz values of the preconditioned operator in [{:.3e}, "
                                "{:.3e}]",
                                corrected_record.ritz_min, corrected_record.ritz_max);
          Mpi::Warning(
              "Self-consistent response-corrected operator is not positive definite: PCG "
              "met "
              "{:d} search direction(s) with (Ap, p) <= 0 (first at iteration {:d} of "
              "{:d}; "
              "{}). The globally coupled correction K + P^T (Q_fab - Q_thin) P has a "
              "negative-energy mode on this mesh (the device's own response to a coupon "
              "trace mode is below the coupon's thin response by more than the fabricated "
              "defect); corrected energies and participations are reported as unavailable "
              "instead of evaluating the meaningless field.\n",
              corrected_record.negative_curvature_count,
              corrected_record.first_negative_curvature_iteration,
              corrected_record.iterations, ritz_range);
        }
        else
        {
          Vector E_corrected(Grad.Height()), D_corrected;
          E_corrected = 0.0;
          Grad.AddMult(V_corrected[step], E_corrected, -1.0);
          const Vector *D_corrected_ptr = nullptr;
          if (post_op.NeedsRecoveredElectricFlux())
          {
            D_corrected.SetSize(estimator->GetFluxSpace().GetTrueVSize());
            D_corrected.UseDevice(true);
            estimator->RecoverFlux(E_corrected, D_corrected);
            D_corrected_ptr = &D_corrected;
          }

          const auto target_interfaces = response_correction->GetTargetInterfaces();
          {
            BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
            corrected_energies = post_op.GetElectrostaticEnergies(
                V_corrected[step], E_corrected, D_corrected_ptr, &target_interfaces);
            if (!uncovered_portions.empty())
            {
              // The uncovered requirements' within-R energy of the corrected field (the
              // field whose outside energy the self-consistent column carries).
              corrected_uncovered = post_op.GetElectrostaticUncoveredEdgeEnergies(
                  target_interfaces, uncovered_portions);
            }
            if (!domain_boundary_portions.empty())
            {
              corrected_domain_boundary = post_op.GetElectrostaticUncoveredEdgeEnergies(
                  target_interfaces, domain_boundary_portions);
            }
            if (!spatial_cluster_priority_portions.empty())
            {
              corrected_spatial_cluster_priority =
                  post_op.GetElectrostaticUncoveredEdgeEnergies(
                      target_interfaces, spatial_cluster_priority_portions);
            }
            if (!clipped_window_portions.empty())
            {
              // The (ii-a) diagnostic on the corrected field (evaluation 2).
              const auto whole = ClippedRawEnergies(clipped_window_portions);
              const auto outside = ClippedRawEnergies(clipped_outside_portions);
              AppendClippedRows(idx, V_corrected[step], false, whole, outside);
            }
            if (!raw_conductor_consistency.empty())
            {
              // The operator has probes: the same list on the corrected field (zero
              // entries while nothing is excluded, like the raw term).
              corrected_conductor_consistency =
                  post_op.GetElectrostaticUncoveredEdgeEnergies(
                      target_interfaces, conductor_consistency_portions);
            }
          }
          SurfaceResponseOperator::ElectrostaticResponse corrected_response;
          {
            BlockTimer coupon_timer(Timer::POSTPRO_RESPONSE_COUPON);
            if (region_targets.empty())
            {
              corrected_response = response_correction->GetElectrostaticResponse(
                  V_corrected[step], false, patch_energy_response);
            }
            else
            {
              // Ruling 2251: the shares of the corrected field (the post operator holds
              // it: the field of this column's kept raw terms).
              const auto region_shares = RegionShares(2);
              corrected_response = response_correction->GetElectrostaticResponse(
                  V_corrected[step], false, patch_energy_response, nullptr, &region_shares);
              response_correction->RecordRegionResponse(
                  idx, "SelfConsistent", corrected_response.region_targets,
                  corrected_response.fabricated_surface_energy,
                  corrected_response.fabricated_surface_energy_fixed_flux,
                  iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0));
            }
          }
          corrected_energies = ApplyResponse(
              std::move(corrected_energies), corrected_response.domain_correction,
              RegionFabricated(corrected_response.fabricated_surface_energy,
                               corrected_response.region_targets, false),
              corrected_uncovered, corrected_domain_boundary,
              corrected_spatial_cluster_priority, corrected_conductor_consistency);
          corrected_contributions = std::move(corrected_response.model_contributions);
          corrected_patch_contributions = std::move(corrected_response.patch_contributions);
        }
      }

      if (response_correction->HasSurfaceResponse())
      {
        corrected_results.push_back(
            CorrectedResult{idx,
                            std::move(raw_energies),
                            std::move(postprocessed_fixed_trace),
                            std::move(postprocessed_fixed_flux),
                            std::move(corrected_energies),
                            std::move(response.model_contributions),
                            std::move(corrected_contributions),
                            std::move(response.patch_contributions),
                            std::move(corrected_patch_contributions),
                            std::move(raw_uncovered),
                            std::move(corrected_uncovered),
                            std::move(raw_domain_boundary),
                            std::move(corrected_domain_boundary),
                            std::move(raw_spatial_cluster_priority),
                            std::move(corrected_spatial_cluster_priority),
                            std::move(raw_conductor_consistency),
                            std::move(corrected_conductor_consistency),
                            response.trace_closure_spread,
                            response.maximum_trace_closure_spread,
                            response.response_weighted_trace_closure_spread,
                            response.trace_closure_response_failure_fraction,
                            postprocess_response,
                            self_consistent_response,
                            response.confident});
        if (response_correction->ConductorConsistencyPerExcitation() &&
            !response_correction->GetConductorConsistencyPassPatches().empty() &&
            corrected_results.back().has_postprocessed)
        {
          // Lane ccpx (decisions 2226 / 2229): the designer's line of this excitation: the
          // raw share of the patches excluded HERE, per target interface and its Type, of
          // the corrected interface energy (evaluation 0; the table's Total rows).
          const auto &result = corrected_results.back();
          const auto targets = response_correction->GetTargetInterfaces();
          std::string shares;
          for (const auto &[interface, type] : post_op.GetInterfaceTypes(&targets))
          {
            double share = std::numeric_limits<double>::quiet_NaN();
            const auto energy = result.raw_conductor_consistency.find(interface);
            const auto corrected =
                result.postprocessed_fixed_trace.interfaces.find(interface);
            if (energy != result.raw_conductor_consistency.end() &&
                corrected != result.postprocessed_fixed_trace.interfaces.end() &&
                corrected->second.energy > 0.0)
            {
              share = energy->second.energy / corrected->second.energy;
            }
            shares += fmt::format("{}{:d} ({}) {:.3e}", shares.empty() ? "" : ", ",
                                  interface, type, share);
          }
          Mpi::Print(" Conductor-consistency raw share at excitation {:d} (the patches "
                     "excluded at this excitation only; of the corrected interface energy, "
                     "evaluation 0): {}\n",
                     idx, shares);
        }
        // The per-patch rows of this source (root; PatchEnergy), then the gathered vectors
        // released: the result keeps its O(interfaces + models) scalars.
        WritePatchEnergyRows(corrected_results.back());
        corrected_results.back().raw_patch_contributions.clear();
        corrected_results.back().raw_patch_contributions.shrink_to_fit();
        corrected_results.back().corrected_patch_contributions.clear();
        corrected_results.back().corrected_patch_contributions.shrink_to_fit();
        // Lane L1: this source's coverage, measured and written before the next source.
        coverage_conductor_consistency_counts.push_back(
            response_correction->GetConductorConsistencyPortions().size());
        if (response_correction->ConductorConsistencyPerExcitation())
        {
          coverage_conductor_consistency_patches.push_back(
              response_correction->GetConductorConsistencyPassPatches());
          conductor_consistency_source_records.push_back(
              response_correction->GetConductorConsistencyPortionRecord());
        }
        if (root)
        {
          coverage_measurements.emplace_back(
              idx, MeasureSourceCoverage(
                       coverage_partitions.back().second, corrected_results.back().raw,
                       postprocess_response ? &region_fabricated : nullptr));
          WriteSourceCoverage(idx, coverage_measurements.back().second);
        }
        // Lane rescensus-product (decision 1817): this source's census cells (collective;
        // the root writes them; a SelfConsistent-only run has no evaluation-0 energy).
        const auto resolution_cells = response_correction->MeasureResolution(
            postprocess_response ? &resolution_sums : nullptr,
            iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0));
        if (root)
        {
          const bool first = coverage_measurements.size() == 1;
          std::ofstream output(resolution_path, first ? std::ios::out | std::ios::trunc
                                                      : std::ios::out | std::ios::app);
          MFEM_VERIFY(output,
                      "Unable to open " << resolution_path.string() << " for writing!");
          output << (first ? ResolutionTableHeader("i") : "");
          WriteResolutionRows(output, std::to_string(idx), resolution_cells);
          for (const auto &cell : resolution_cells)
          {
            response_correction->GetResolutionRunRecord().entries.push_back(
                {{{"i", static_cast<double>(idx)}}, cell});
          }
        }
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

    // Decision 1469 (1) (i): the memory census after the first source (log only; a no-op
    // unless PALACE_RESPONSE_MEMORY_CENSUS=1).
    if (step == 0 && response_correction)
    {
      response_correction->PrintResponseMemoryCensus("after source 1",
                                                     response_geometry.get());
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
    if (response_correction)
    {
      BlockTimer response_timer(Timer::POSTPRO_RESPONSE);
      // Lane L4 (decision 1829; ERRATUM-1829 E-1): the FINAL coverage of every source is
      // measured BEFORE the corrected matrices are written (two independent computations;
      // the order changes no value), so that under Coverage.Policy "MarkNotApplied" entry
      // (i, j) is marked iff ANY cell of source i's or source j's final measurement is.
      // Lane L1 (ruling 1654 (III) Q3): terminal-C-corrected(i, j) is formed from the
      // fields of sources i and j with the FINAL patch weights (after the last
      // conductor-consistency exclusion), so each source's coverage is re-measured with the
      // FINAL conductor-consistency portions where a later exclusion changed them (else its
      // own measurement is the final one); entry (i, j) = the worse of sources i and j.
      // The per-source rows keep their own measurement (what they print).
      // Lane ccpx (ruling 2229 D3 = X1): under PerExcitation the FINAL portions are the
      // union's (every patch excluded at one or more excitations = the Sticky final list)
      // and a source's own measurement is final where its own exclusion SET is the union.
      const auto &final_record =
          response_correction->GetConductorConsistencyUnionPortionRecord();
      const std::size_t final_count = final_record.portions.size();
      const bool per_excitation = response_correction->ConductorConsistencyPerExcitation();
      auto OwnMeasurementIsFinal = [&](std::size_t i)
      {
        return per_excitation
                   ? coverage_conductor_consistency_patches.at(i) ==
                         response_correction->GetConductorConsistencyUnionPatches()
                   : coverage_conductor_consistency_counts.at(i) == final_count;
      };
      std::vector<CoverageMeasurement> final_coverage;
      std::vector<bool> remeasured;
      std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
          final_conductor_consistency;
      for (const auto &portion : final_record.portions)
      {
        final_conductor_consistency.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
      for (std::size_t i = 0; i < coverage_partitions.size() && i < V.size(); i++)
      {
        if (OwnMeasurementIsFinal(i) || V[i].Size() == 0)
        {
          if (root)
          {
            final_coverage.push_back(coverage_measurements.at(i).second);
            remeasured.push_back(false);
          }
          continue;
        }
        BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
        E = 0.0;
        Grad.AddMult(V[i], E, -1.0);
        if (post_op.NeedsRecoveredElectricFlux())
        {
          estimator->RecoverFlux(E, D);
        }
        const auto raw = post_op.GetElectrostaticEnergies(
            V[i], E, post_op.NeedsRecoveredElectricFlux() ? &D : nullptr);
        const auto missing_by_feature =
            SurfacePostOperator::CoverageMissingByFeature(uncovered_portions);
        const std::vector<SurfacePostOperator::CoveragePartitionClass> kept_classes = {
            {"Missing", &missing_by_feature},
            {"DomainBoundary", &domain_boundary_portions},
            {"SpatialClusterPriority", &spatial_cluster_priority_portions},
            {"ConductorConsistency", &final_conductor_consistency}};
        std::vector<SurfacePostOperator::CoveragePartitionClass> remainder_classes;
        if (response_correction->ExclusionEnergyRequested())
        {
          remainder_classes.push_back({"Excluded", &exclusion_energy_portions});
        }
        // Ruling 2295 (A): a box whose patch is in the FINAL set holds nothing.
        auto partition = post_op.GetElectrostaticCoveragePartition(
            response_correction->GetCoverageInterfaces(), kept_classes, model_raw_portions,
            remainder_classes, dropped_windows, final_record.PatchSet());
        SurfacePostOperator::MoveCoverageServedLabel(partition, "unpatched",
                                                     "UnpatchedPairs");
        if (root)
        {
          // The factor of this re-measurement is not printed (the matrix coverage is a
          // status); the source's fabricated energy stands in.
          const auto &own = coverage_measurements.at(i).second;
          std::map<int, double> fabricated;
          for (const auto &[interface, summary] : own.interfaces)
          {
            fabricated[interface] = summary.fabricated;
          }
          final_coverage.push_back(
              own.measured ? MeasureSourceCoverage(partition, raw, &fabricated) : own);
          remeasured.push_back(own.measured);
        }
      }
      if (root && final_coverage.size() == laplace_op.GetSources().size())
      {
        WriteMatrixCoverage(laplace_op.GetSources(), final_coverage, remeasured);
      }
      // E-1: per source, whether its final measurement holds a marked cell (root;
      // fail-closed: a source without a final measurement, or one without cells, is
      // marked).
      std::vector<bool> source_marked;
      if (root && response_correction->CoverageMarkNotApplied())
      {
        for (std::size_t i = 0; i < laplace_op.GetSources().size(); i++)
        {
          source_marked.push_back(final_coverage.size() != laplace_op.GetSources().size() ||
                                  final_coverage[i].MarksWholeMeasurement());
        }
      }
      // Lane ccpx (ruling 2229 D3 = X1): the matrices take ONE operator, construction
      // minus the union (the Sticky final weights; a no-op under Sticky); the corrected
      // fields V_corrected were solved with their own excitation's weights.
      response_correction->ApplyConductorConsistencyUnionWeights();
      PostprocessCorrectedTerminals(
          post_op, laplace_op.GetSources(), V, V_corrected, *response_correction,
          postprocess_response, self_consistent_accepted,
          response_correction->CoverageMarkNotApplied() ? &source_marked : nullptr);
    }
  }
  else if (iodata.solver.electrostatic.response_matrix && !archive_stream_only)
  {
    PostprocessResponseMatrix(post_op, laplace_op, Grad, V, D_basis);
  }

  if (root && !corrected_results.empty())
  {
    using VT = Units::ValueType;
    TableWithCSVFile output(post_dir / "surface-Q-corrected.csv");
    output.table.insert(Column::Integer("source", "i"));
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
    // Lane L1 (DESIGN-v2 s4.3; decision 1636 (1)): the coverage columns APPENDED at the
    // end of the row (every column above keeps its position, name and value): per
    // interface W and its uncovered share, then per Type (first appearance over the
    // ascending target interfaces) the source's cell.
    const auto coverage_targets = response_correction->GetCoverageInterfaces();
    std::vector<std::string> coverage_types;
    for (const auto &[interface, type] : post_op.GetInterfaceTypes(&coverage_targets))
    {
      (void)interface;
      if (std::find(coverage_types.begin(), coverage_types.end(), type) ==
          coverage_types.end())
      {
        coverage_types.push_back(type);
      }
    }
    for (const auto &[interface, data] : interfaces)
    {
      (void)data;
      output.table.insert(fmt::format("coverage_within_{}", interface),
                          fmt::format("W within-R raw[{}] (J)", interface));
      output.table.insert(fmt::format("coverage_uncovered_share_{}", interface),
                          fmt::format("uncovered raw share of W[{}]", interface));
    }
    for (const auto &type : coverage_types)
    {
      output.table.insert(fmt::format("coverage_uncovered_share_{}", type),
                          fmt::format("uncovered raw share of W {}", type));
      output.table.insert(fmt::format("coverage_method_gap_share_{}", type),
                          fmt::format("method-gap raw share of W {}", type));
      output.table.insert(fmt::format("coverage_pending_share_{}", type),
                          fmt::format("pending-model share of W {}", type));
      output.table.insert(fmt::format("coverage_factor_{}", type),
                          fmt::format("covered correction factor {}", type));
      output.table.insert(Column::Text(fmt::format("coverage_status_{}", type),
                                       fmt::format("coverage status {}", type)));
      output.table.insert(Column::Text(fmt::format("regime_status_{}", type),
                                       fmt::format("regime status {}", type)));
    }
    // Lane L4 (decision 1829; DESIGN-L4-v2 s2.3 E1): under Coverage.Policy
    // "MarkNotApplied" the corrected columns print nan where their cell is marked (E_elec:
    // any marked cell of the source; E / p / Q[k]: the Type of interface k), the value in
    // the appended partial column, and the exposure column (the raw within-R energy the
    // marked cells leave uncorrected over E_elec raw) follows; under "Warn" none of these
    // columns exists and every cell prints as before.
    const bool marking = response_correction->CoverageMarkNotApplied();
    std::vector<std::string> marked_source_columns = {"domain_postprocessed_fixed_trace",
                                                      "domain_postprocessed_fixed_flux",
                                                      "domain_corrected"};
    std::vector<std::string> marked_interface_columns;
    for (const auto &[interface, data] : interfaces)
    {
      (void)data;
      for (const char *name :
           {"energy_postprocessed_fixed_trace_{}",
            "participation_postprocessed_fixed_trace_{}",
            "quality_postprocessed_fixed_trace_{}", "energy_postprocessed_fixed_flux_{}",
            "participation_postprocessed_fixed_flux_{}",
            "quality_postprocessed_fixed_flux_{}", "energy_corrected_{}",
            "participation_corrected_{}", "quality_corrected_{}"})
      {
        marked_interface_columns.push_back(fmt::format(fmt::runtime(name), interface));
      }
    }
    if (marking)
    {
      for (const auto *names : {&marked_source_columns, &marked_interface_columns})
      {
        for (const auto &name : *names)
        {
          const std::string header = output.table[name].header_text;
          output.table.insert(name + "_partial", header + kCoveragePartialColumnSuffix);
        }
      }
      output.table.insert("coverage_marked_uncovered_share",
                          "not-applied uncovered raw share of E_elec");
      // Lane L4 round 2 (decision 1914): the explicit marks, last: the source's E_elec
      // cells, then per interface its E / p / Q cells.
      output.table.insert("coverage_mark_source", kCoverageMarkSourceColumn);
      output.table["coverage_mark_source"].print_as_int = true;
      for (const auto &[interface, data] : interfaces)
      {
        (void)data;
        const std::string key = fmt::format("coverage_mark_{}", interface);
        output.table.insert(
            key, fmt::format(fmt::runtime(kCoverageMarkIndexedInterfaceColumn), interface));
        output.table[key].print_as_int = true;
      }
    }
    MFEM_VERIFY(coverage_measurements.size() == corrected_results.size(),
                "The coverage needs one measurement per corrected source!");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::size_t coverage_index = 0;
    for (const auto &result : corrected_results)
    {
      // E1: this source's own measurement (the coverage columns' row).
      const auto &row_coverage = coverage_measurements[coverage_index].second;
      const bool source_marked = marking && row_coverage.MarksWholeMeasurement();
      auto InterfaceMarked = [&](int interface)
      {
        const auto type = row_coverage.types.find(interface);
        if (!marking || type == row_coverage.types.end())
        {
          return false;  // an untargeted interface: not corrected
        }
        return !row_coverage.FindCell(type->second) ||
               row_coverage.IsTypeMarked(type->second);
      };
      auto PrintMarked = [&](const std::string &name, double value, bool marked)
      {
        if (!marking)
        {
          output.table[name] << value;
          return;
        }
        output.table[name] << (marked ? nan : value);
        output.table[name + "_partial"] << (marked ? value : nan);
      };
      output.table["source"] << result.source;
      output.table["domain_raw"]
          << iodata.units.Dimensionalize<VT::ENERGY>(result.raw.domain);
      PrintMarked(
          "domain_postprocessed_fixed_trace",
          iodata.units.Dimensionalize<VT::ENERGY>(result.postprocessed_fixed_trace.domain),
          source_marked);
      PrintMarked(
          "domain_postprocessed_fixed_flux",
          iodata.units.Dimensionalize<VT::ENERGY>(result.postprocessed_fixed_flux.domain),
          source_marked);
      PrintMarked("domain_corrected",
                  iodata.units.Dimensionalize<VT::ENERGY>(result.corrected.domain),
                  source_marked);
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
        const bool interface_marked = InterfaceMarked(interface);
        PrintMarked(fmt::format("energy_postprocessed_fixed_trace_{}", interface),
                    iodata.units.Dimensionalize<VT::ENERGY>(fixed_trace.energy),
                    interface_marked);
        PrintMarked(fmt::format("participation_postprocessed_fixed_trace_{}", interface),
                    p_fixed_trace, interface_marked);
        PrintMarked(fmt::format("quality_postprocessed_fixed_trace_{}", interface),
                    Quality(p_fixed_trace, fixed_trace.loss_tangent), interface_marked);
        PrintMarked(fmt::format("energy_postprocessed_fixed_flux_{}", interface),
                    iodata.units.Dimensionalize<VT::ENERGY>(fixed_flux.energy),
                    interface_marked);
        PrintMarked(fmt::format("participation_postprocessed_fixed_flux_{}", interface),
                    p_fixed_flux, interface_marked);
        PrintMarked(fmt::format("quality_postprocessed_fixed_flux_{}", interface),
                    Quality(p_fixed_flux, fixed_flux.loss_tangent), interface_marked);
        PrintMarked(fmt::format("energy_corrected_{}", interface),
                    iodata.units.Dimensionalize<VT::ENERGY>(corrected.energy),
                    interface_marked);
        PrintMarked(fmt::format("participation_corrected_{}", interface), p_corrected,
                    interface_marked);
        PrintMarked(fmt::format("quality_corrected_{}", interface),
                    Quality(p_corrected, corrected.loss_tangent), interface_marked);
        const auto closure_spread = result.trace_closure_spread.find(interface);
        output.table[fmt::format("trace_closure_spread_{}", interface)]
            << (!result.has_postprocessed ||
                        closure_spread == result.trace_closure_spread.end()
                    ? nan
                    : closure_spread->second);
      }
      const auto &[coverage_source, coverage] = coverage_measurements[coverage_index++];
      MFEM_VERIFY(coverage_source == result.source,
                  "The coverage measurements are out of source order!");
      for (const auto &[interface, data] : interfaces)
      {
        (void)data;
        const auto summary = coverage.interfaces.find(interface);
        // Ruling 1694 (A): an untargeted interface's within-R region is empty (W = 0).
        const double within = summary != coverage.interfaces.end() ? summary->second.within
                              : coverage_targets.count(interface) == 0 ? 0.0
                                                                       : nan;
        output.table[fmt::format("coverage_within_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(within);
        output.table[fmt::format("coverage_uncovered_share_{}", interface)]
            << (summary != coverage.interfaces.end() && within > 0.0
                    ? (within - summary->second.served) / within
                    : nan);
      }
      for (const auto &type : coverage_types)
      {
        const auto *cell = coverage.FindCell(type);
        output.table[fmt::format("coverage_uncovered_share_{}", type)]
            << (cell ? cell->uncovered_share : nan);
        output.table[fmt::format("coverage_method_gap_share_{}", type)]
            << (cell ? cell->method_gap_share : nan);
        output.table[fmt::format("coverage_pending_share_{}", type)]
            << (cell ? cell->pending_share : nan);
        output.table[fmt::format("coverage_factor_{}", type)]
            << (cell ? cell->factor : nan);
        output.table[fmt::format("coverage_status_{}", type)] << std::string(
            CoverageStatusName(cell ? cell->status : CoverageStatus::NOT_MEASURED));
        output.table[fmt::format("regime_status_{}", type)]
            << std::string(RegimeStatusName(cell ? cell->regime : RegimeStatus::NO_MODEL));
      }
      if (marking)
      {
        output.table["coverage_marked_uncovered_share"]
            << row_coverage.MarkedUncoveredEnergy() / result.raw.domain;
        output.table["coverage_mark_source"] << (source_marked ? 1.0 : 0.0);
        for (const auto &[interface, data] : interfaces)
        {
          (void)data;
          output.table[fmt::format("coverage_mark_{}", interface)]
              << (InterfaceMarked(interface) ? 1.0 : 0.0);
        }
      }
    }
    output.WriteFullTableTrunc();

    TableWithCSVFile model_output(post_dir / "surface-response-model-energy.csv");
    model_output.table.insert(Column::Integer("source", "source"));
    model_output.table.insert(Column::Integer("evaluation", "evaluation"));
    model_output.table.insert(Column::Integer("model", "model"));
    model_output.table.insert(Column::Integer("patch_count", "patch count"));
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

    // Decision 394 F2, F-DB-a (decisions 442 / 454) and S7diag R1 (decision 758): the raw
    // within-R energy kept in the corrected interface energies by the uncovered
    // requirements, by the DomainBoundary raw claims and by the conductor-consistency-
    // excluded patches' raw claims, each in its own table, per source, evaluation (the
    // codes of surface-response-model-energy.csv: 0 / 1 the raw field kept in the
    // fixed-trace / fixed-flux energies, 2 the corrected field kept in the self-consistent
    // energy) and feature type (one row per type and a Total row), with the portion count
    // and length of that type, the length removed from that type's portions inside the
    // matched clusters' support boxes at placement (decision 399 MAJOR-1; the uncovered
    // table) or the geometric cell count (the DomainBoundary and ConductorConsistency
    // tables) and, per target interface, the energy and its share of the evaluation's
    // corrected interface energy. Each table is written only when it has portions.
    struct PortionTypeRow
    {
      int portions = 0;
      double length = 0.0;
      double extra = 0.0;  // clipped length (uncovered) or geometric cells (DomainBoundary)
    };
    using PortionEnergies = std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>;
    // `source_portions_by_type` (lane ccpx, Exclusion PerExcitation; the corrected_results
    // order): each source's rows carry its OWN portion columns over the types of
    // `portions_by_type` (0 where the source has none); null: every row the given ones.
    auto WritePortionEnergyTable =
        [&](const std::filesystem::path &path,
            const std::map<std::string, PortionTypeRow> &portions_by_type,
            const std::string &extra_header, const std::string &energy_label,
            const std::function<const PortionEnergies &(const CorrectedResult &, bool)>
                &energies_of,
            const std::vector<std::map<std::string, PortionTypeRow>>
                *source_portions_by_type = nullptr)
    {
      PortionTypeRow portions_total;
      for (const auto &[type, entry] : portions_by_type)
      {
        portions_total.portions += entry.portions;
        portions_total.length += entry.length;
        portions_total.extra += entry.extra;
      }
      // A plain writer: the type column is text (the table columns are numeric). Both
      // conditional tables were removed above, unconditionally, before this conditional
      // write (the previous adaptive iteration's symlink into its archive folder, and any
      // table of an earlier attempt in the same directory, never survive a solve without
      // portions; decision 473 (3)).
      std::ofstream output(path);
      MFEM_VERIFY(output, "Unable to open " << path.filename().string() << " for writing!");
      output << "source,evaluation,type,portions,length (m)," << extra_header;
      for (const auto &[interface, data] : interfaces)
      {
        (void)data;
        output << fmt::format(",{} raw energy[{}] (J),{} share[{}]", energy_label,
                              interface, energy_label, interface);
      }
      output << '\n';
      auto Append = [&](int source, int evaluation, const std::string &type,
                        const PortionTypeRow &portions, const PortionEnergies &energies_map,
                        const EnergyData &energies, bool available, bool extra_is_length)
      {
        output << fmt::format(
            "{},{},{},{},{},", source, evaluation, type, portions.portions,
            CsvNumber(iodata.units.Dimensionalize<VT::LENGTH>(portions.length)));
        if (extra_is_length)
        {
          output << CsvNumber(iodata.units.Dimensionalize<VT::LENGTH>(portions.extra));
        }
        else
        {
          output << fmt::format("{}", static_cast<int>(std::lround(portions.extra)));
        }
        for (const auto &[interface, interface_data] : interfaces)
        {
          (void)interface_data;
          double energy = nan, share = nan;
          const auto it = energies_map.find(interface);
          if (available && it != energies_map.end())
          {
            if (type == "Total")
            {
              energy = it->second.energy;
            }
            else
            {
              const auto by_type = it->second.by_type.find(type);
              energy = by_type != it->second.by_type.end() ? by_type->second : 0.0;
            }
            const auto corrected = energies.interfaces.find(interface);
            share = corrected != energies.interfaces.end() && corrected->second.energy > 0.0
                        ? energy / corrected->second.energy
                        : nan;
            energy = iodata.units.Dimensionalize<VT::ENERGY>(energy);
          }
          output << fmt::format(",{},{}", CsvNumber(energy), CsvNumber(share));
        }
        output << '\n';
      };
      const bool extra_is_length = extra_header.find("(m)") != std::string::npos;
      for (std::size_t k = 0; k < corrected_results.size(); k++)
      {
        const auto &result = corrected_results[k];
        PortionTypeRow source_total = portions_total;
        if (source_portions_by_type)
        {
          source_total = PortionTypeRow();
          for (const auto &[type, entry] : source_portions_by_type->at(k))
          {
            source_total.portions += entry.portions;
            source_total.length += entry.length;
            source_total.extra += entry.extra;
          }
        }
        for (const auto &[type, all_portions] : portions_by_type)
        {
          PortionTypeRow portions = all_portions;
          if (source_portions_by_type)
          {
            const auto own = source_portions_by_type->at(k).find(type);
            portions = own != source_portions_by_type->at(k).end() ? own->second
                                                                   : PortionTypeRow();
          }
          Append(result.source, 0, type, portions, energies_of(result, false),
                 result.postprocessed_fixed_trace, result.has_postprocessed,
                 extra_is_length);
          Append(result.source, 1, type, portions, energies_of(result, false),
                 result.postprocessed_fixed_flux, result.has_postprocessed,
                 extra_is_length);
          Append(result.source, 2, type, portions, energies_of(result, true),
                 result.corrected, result.has_self_consistent, extra_is_length);
        }
        Append(result.source, 0, "Total", source_total, energies_of(result, false),
               result.postprocessed_fixed_trace, result.has_postprocessed, extra_is_length);
        Append(result.source, 1, "Total", source_total, energies_of(result, false),
               result.postprocessed_fixed_flux, result.has_postprocessed, extra_is_length);
        Append(result.source, 2, "Total", source_total, energies_of(result, true),
               result.corrected, result.has_self_consistent, extra_is_length);
      }
    };
    auto PortionLength = [](const auto &portion)
    {
      double length2 = 0.0;
      for (int d = 0; d < 3; d++)
      {
        length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
      }
      return std::sqrt(length2);
    };
    // The four tables are removed UNCONDITIONALLY before the conditional write (decision
    // 473 MAJOR-3): a stale table of an earlier run or attempt in the output directory
    // would otherwise read as a term this solve did not have.
    std::filesystem::remove(post_dir / "surface-response-uncovered-energy.csv");
    std::filesystem::remove(post_dir / "surface-response-domain-boundary-energy.csv");
    std::filesystem::remove(post_dir /
                            "surface-response-spatial-cluster-priority-energy.csv");
    std::filesystem::remove(post_dir / "surface-response-conductor-consistency-energy.csv");
    // Section 8: a stale exclusion-energy table of an earlier run is removed alike; the
    // table of this run (when Diagnostics.ExclusionEnergy is on) was written per source
    // above and is not touched here.
    if (!exclusion_energy_csv)
    {
      std::filesystem::remove(post_dir / "surface-response-exclusion-energy.csv");
    }
    // PatchEnergy (decision 922 B-1): surface-response-model-raw-energy.csv (per source and
    // model-set label the raw within-R energy per postprocessing interface, 0 for a
    // non-target interface) and surface-response-model-raw-groups.csv (per label and
    // model: patches, weight x length), the Maxwell tables' columns with the electrostatic
    // `source` key.
    std::filesystem::remove(post_dir / "surface-response-model-raw-energy.csv");
    std::filesystem::remove(post_dir / "surface-response-model-raw-groups.csv");
    if (!raw_model_partitions.empty())
    {
      std::ofstream energy(post_dir / "surface-response-model-raw-energy.csv");
      MFEM_VERIFY(energy, "Unable to open surface-response-model-raw-energy.csv!");
      energy << "source,evaluation,label";
      for (const auto &[interface, data] : interfaces)
      {
        (void)data;
        energy << fmt::format(",raw energy[{}] (J)", interface);
      }
      energy << '\n';
      for (const auto &[source, partition_by_interface] : raw_model_partitions)
      {
        std::set<std::string> labels;
        for (const auto &[interface, partition] : partition_by_interface)
        {
          (void)interface;
          for (const auto &[label, value] : partition)
          {
            (void)value;
            labels.insert(label);
          }
        }
        for (const auto &label : labels)
        {
          energy << fmt::format("{},{},{}", source, 0, label);
          for (const auto &[interface, data] : interfaces)
          {
            (void)data;
            const auto partition = partition_by_interface.find(interface);
            double value = 0.0;
            if (partition != partition_by_interface.end())
            {
              const auto entry = partition->second.find(label);
              value = entry != partition->second.end() ? entry->second : 0.0;
            }
            energy << "," << CsvNumber(iodata.units.Dimensionalize<VT::ENERGY>(value));
          }
          energy << '\n';
        }
      }
      std::ofstream groups(post_dir / "surface-response-model-raw-groups.csv");
      MFEM_VERIFY(groups, "Unable to open surface-response-model-raw-groups.csv!");
      groups << "label,model,patches,weight length (m)\n";
      for (const auto &member : response_correction->GetModelRawClaimPortions().members)
      {
        groups << fmt::format(
            "{},{},{},{}\n", member.label, member.model, member.patches,
            CsvNumber(iodata.units.Dimensionalize<VT::LENGTH>(member.weight_length)));
      }
    }
    if (!uncovered_portions.empty())
    {
      std::map<std::string, PortionTypeRow> portions_by_type;
      for (const auto &portion : response_correction->GetUncoveredPortions())
      {
        auto &entry = portions_by_type[portion.topology];
        entry.portions++;
        entry.length += PortionLength(portion);
      }
      for (const auto &clip :
           response_correction->GetUncoveredSpatialSupportClipping().clips)
      {
        portions_by_type[clip.topology].extra += clip.length;
      }
      WritePortionEnergyTable(
          post_dir / "surface-response-uncovered-energy.csv", portions_by_type,
          "clipped by spatial support (m)", "uncovered",
          [](const CorrectedResult &result, bool corrected) -> const PortionEnergies &
          { return corrected ? result.corrected_uncovered : result.raw_uncovered; });
    }
    if (!domain_boundary_portions.empty())
    {
      // F-DB-a: the DomainBoundary share per type (type = the excluded patch's model
      // topology; the geometric cell count = translational own-edge intervals, the
      // co-located first-order split patches of one cell counted once).
      std::map<std::string, PortionTypeRow> portions_by_type;
      const auto &record = response_correction->GetDomainBoundaryPortionRecord();
      for (std::size_t i = 0; i < record.portions.size(); i++)
      {
        auto &entry = portions_by_type[record.portions[i].topology];
        entry.portions++;
        entry.length += PortionLength(record.portions[i]);
        entry.extra += record.translational[i] ? 1.0 : 0.0;
      }
      WritePortionEnergyTable(
          post_dir / "surface-response-domain-boundary-energy.csv", portions_by_type,
          "geometric cells", "domainboundary",
          [](const CorrectedResult &result, bool corrected) -> const PortionEnergies &
          {
            return corrected ? result.corrected_domain_boundary
                             : result.raw_domain_boundary;
          });
    }
    if (!spatial_cluster_priority_portions.empty())
    {
      // OI-SCP-HOLE (decision 922): the SpatialClusterPriority share per type (the
      // DomainBoundary table's columns; the extra column = the dropped coupons folded into
      // the portion, 0 for a plain arm window).
      std::map<std::string, PortionTypeRow> portions_by_type;
      const auto &record =
          response_correction->GetSpatialClusterPriorityExclusions().portions;
      for (std::size_t i = 0; i < record.portions.size(); i++)
      {
        auto &entry = portions_by_type[record.portions[i].topology];
        entry.portions++;
        entry.length += PortionLength(record.portions[i]);
        entry.extra += record.translational[i] ? 1.0 : 0.0;
      }
      WritePortionEnergyTable(
          post_dir / "surface-response-spatial-cluster-priority-energy.csv",
          portions_by_type, "geometric cells", "spatialclusterpriority",
          [](const CorrectedResult &result, bool corrected) -> const PortionEnergies &
          {
            return corrected ? result.corrected_spatial_cluster_priority
                             : result.raw_spatial_cluster_priority;
          });
    }
    // The (ii-a) companion table (OI-SCP-BOXGROW, decisions 958 (b) / 960 / 969): one row
    // per dropped coupon x (source, evaluation 0 fixed trace / 1 fixed flux on the raw
    // field / 2 the corrected field); per target interface the DERIVED columns Fab (weight
    // x t^T K t of the dropped coupon on the device trace at its own basis points),
    // RawWhole (the raw within-R energy of its whole window), RawOut (of the window outside
    // every cluster box = the SCP table's add-back, per patch), Ratio = Fab / RawWhole and
    // ClippedCorrection = (Ratio - weight) x RawOut (the proportional clip of the corner's
    // missing correction; the co-located patches of a first-order split sum to the
    // window's). Never in a total.
    std::filesystem::remove(
        post_dir / "surface-response-spatial-cluster-priority-clipped-correction.csv");
    if (!clipped_correction_rows.empty())
    {
      std::ofstream output(
          post_dir / "surface-response-spatial-cluster-priority-clipped-correction.csv");
      MFEM_VERIFY(output,
                  "Unable to open "
                  "surface-response-spatial-cluster-priority-clipped-correction.csv!");
      output
          << "# DERIVED (decisions 958 (b) / 960 / 969; OI-SCP-BOXGROW (ii-a)): parents = "
             "the dropped coupon's response matrices K_type on the device trace sampled "
             "at its own basis points (fab energy), the raw within-R interface energy of "
             "its whole window (raw whole energy) and of the window outside every "
             "cluster box (raw outside energy = the SpatialClusterPriority add-back of "
             "surface-response-spatial-cluster-priority-energy.csv per patch); ratio = "
             "fab / raw whole; clipped correction = (ratio - weight) x raw outside. The "
             "weight column is the dropped coupon's model BLEND weight, the weight it "
             "would carry if applied (1 on a whole-weight coupon), never the ownership "
             "zeroing: the Weight column of "
             "surface-response-spatial-cluster-priority-energy.csv reads 0 on a corner "
             "the vertex rule owns and is not this column. Evaluation 0 = fixed trace on "
             "the raw field, 1 = fixed flux on the raw field, 2 = the self-consistent "
             "corrected field. Diagnostic only: no total carries these columns.\n";
      output << "source,evaluation,patch,model,feature,weight";
      for (const auto &[interface, data] : interfaces)
      {
        (void)data;
        output << fmt::format(",fab energy[{0}] (J),raw whole energy[{0}] (J),raw outside "
                              "energy[{0}] (J),ratio[{0}],clipped correction[{0}] (J)",
                              interface);
      }
      output << '\n';
      const double nan = std::numeric_limits<double>::quiet_NaN();
      for (const auto &row : clipped_correction_rows)
      {
        output << fmt::format("{},{},{},{},{},{}", row.source, row.evaluation,
                              row.diagnostic.patch + 1, row.diagnostic.model,
                              row.diagnostic.feature, CsvNumber(row.diagnostic.weight));
        for (const auto &[interface, data] : interfaces)
        {
          (void)data;
          const auto &fab_map = row.fixed_flux
                                    ? row.diagnostic.fabricated_surface_energy_fixed_flux
                                    : row.diagnostic.fabricated_surface_energy;
          const auto fab_it = fab_map.find(interface);
          const auto whole_it = row.raw_whole.find(interface);
          const auto outside_it = row.raw_outside.find(interface);
          const double fab = fab_it != fab_map.end() ? fab_it->second : nan;
          const double whole = whole_it != row.raw_whole.end() ? whole_it->second : 0.0;
          const double outside =
              outside_it != row.raw_outside.end() ? outside_it->second : 0.0;
          const double ratio = whole > 0.0 && std::isfinite(fab) ? fab / whole : nan;
          const double clipped =
              std::isfinite(ratio) ? (ratio - row.diagnostic.weight) * outside : nan;
          output << fmt::format(
              ",{},{},{},{},{}", CsvNumber(iodata.units.Dimensionalize<VT::ENERGY>(fab)),
              CsvNumber(iodata.units.Dimensionalize<VT::ENERGY>(whole)),
              CsvNumber(iodata.units.Dimensionalize<VT::ENERGY>(outside)), CsvNumber(ratio),
              CsvNumber(iodata.units.Dimensionalize<VT::ENERGY>(clipped)));
        }
        output << '\n';
      }
    }
    if (!response_correction->GetConductorConsistencyUnionPortionRecord().portions.empty())
    {
      // S7diag R1 (decision 758): the ConductorConsistency share per type (the
      // DomainBoundary table's columns). Sticky: the portion rows describe the FINAL list
      // of the solve; an excitation evaluated before the exclusion kept the model and
      // reads 0. PerExcitation (lane ccpx): the rows' types are the union's, each source's
      // portion columns its own exclusions' (an excitation that excluded none reads 0).
      auto ByType = [&PortionLength](const ExcludedPatchPortions &record)
      {
        std::map<std::string, PortionTypeRow> portions_by_type;
        for (std::size_t i = 0; i < record.portions.size(); i++)
        {
          auto &entry = portions_by_type[record.portions[i].topology];
          entry.portions++;
          entry.length += PortionLength(record.portions[i]);
          entry.extra += record.translational[i] ? 1.0 : 0.0;
        }
        return portions_by_type;
      };
      const auto portions_by_type =
          ByType(response_correction->GetConductorConsistencyUnionPortionRecord());
      std::vector<std::map<std::string, PortionTypeRow>> source_portions_by_type;
      for (const auto &record : conductor_consistency_source_records)
      {
        source_portions_by_type.push_back(ByType(record));
      }
      MFEM_VERIFY(!response_correction->ConductorConsistencyPerExcitation() ||
                      source_portions_by_type.size() == corrected_results.size(),
                  "One conductor-consistency portion record per corrected source!");
      WritePortionEnergyTable(
          post_dir / "surface-response-conductor-consistency-energy.csv", portions_by_type,
          "geometric cells", "conductorconsistency",
          [](const CorrectedResult &result, bool corrected) -> const PortionEnergies &
          {
            return corrected ? result.corrected_conductor_consistency
                             : result.raw_conductor_consistency;
          },
          response_correction->ConductorConsistencyPerExcitation()
              ? &source_portions_by_type
              : nullptr);
    }
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
    surface_output.table.insert(Column::Integer("interface", "interface"));
    surface_output.table.insert(Column::Integer("edge", "edge"));
    surface_output.table.insert("distance", "R (m)");
    surface_output.table.insert(Column::Integer("basis_i", "basis_i"));
    surface_output.table.insert(Column::Integer("basis_j", "basis_j"));
    surface_output.table.insert("Q", "Q_ij (J)");
    surface_output.table.insert("Q_normal", "Q_ij normal (J)");
    surface_output.table.insert("Q_tangential", "Q_ij tangential (J)");
    surface_output.table.insert("Q_total", "Q_total_ij (J)");
    surface_output.table.insert("Q_total_normal", "Q_total_ij normal (J)");
    surface_output.table.insert("Q_total_tangential", "Q_total_ij tangential (J)");

    domain_output = TableWithCSVFile(post_dir / "domain-response-matrix.csv");
    domain_output.table.insert(Column::Integer("basis_i", "basis_i"));
    domain_output.table.insert(Column::Integer("basis_j", "basis_j"));
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
    Mpi::Print(
        " Archived response reduction: {:d} sources, {:d} interfaces, quadrature "
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
        // P3-21: an archived recovered flux lives on the flux recovery space of the run
        // that wrote it, this run holds D on its own (the same space for the same mesh and
        // crack form).
        MFEM_VERIFY(D_source.Size() == post_op.GetRecoveredElectricFluxTrueVSize(),
                    "Archived recovered electric flux has an incompatible local size (the "
                    "flux recovery spaces of the archive and of this run differ)!");
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
        post_op.EvaluateInterfaceResponseRow(samples[k],
                                             rows[k].data() + i * samples[k].RowSize());
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
      surface_matrices.emplace(
          samples[k].interface_index,
          post_op.AssembleInterfaceResponseMatrices(samples[k], rows[k].data(),
                                                    static_cast<int>(basis_size)));
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
    output.table.insert(Column::Integer("interface", "interface"));
    output.table.insert(Column::Integer("edge", "edge"));
    output.table.insert("distance", "R (m)");
    output.table.insert(Column::Integer("basis_i", "basis_i"));
    output.table.insert(Column::Integer("basis_j", "basis_j"));
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
    domain_output.table.insert(Column::Integer("basis_i", "basis_i"));
    domain_output.table.insert(Column::Integer("basis_j", "basis_j"));
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

void ElectrostaticSolver::PostprocessCorrectedTerminals(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op,
    const std::map<int, mfem::Array<int>> &terminal_sources, const std::vector<Vector> &V,
    const std::vector<Vector> &V_corrected, const SurfaceResponseOperator &response,
    bool fixed_trace, bool self_consistent, const std::vector<bool> *source_marked) const
{
  // The capacitance of the response-corrected energy. The raw matrix is the bilinear form
  // of the thin-metal operator (PostprocessTerminals); the fixed-trace domain correction
  // adds the bilinear domain defect of the applied coupon patches,
  //         C(i, j) = Vⱼᵀ M Vᵢ + Vⱼᵀ Pᵀ W (Q_fab,dom - Q_thin,dom) P Vᵢ  (∀i, Vᵢ = 1),
  // so that C(i, i) = 2 (E_raw + ΔE_dom,ft) / Vᵢ² (the surface-Q-corrected.csv energies; W
  // are the final patch weights, so the identity holds for the sources solved after the
  // last conductor-consistency exclusion, i.e. for every source when nothing is excluded
  // after the first). The self-consistent matrix is the same fixed-trace form on the
  // corrected fields: the capacitance of the corrected operator K + Pᵀ W D P at the
  // terminal potentials when the translational domain coupling is FixedTrace (under
  // FixedFlux or Disabled the corrected fields solve a different coupling). Fixed flux is
  // not a capacitance (its energy is not that of a fixed-potential ensemble) and is not
  // written.
  const int n = static_cast<int>(V.size());
  MFEM_VERIFY(static_cast<int>(terminal_sources.size()) == n,
              "Terminal source and field counts differ in the corrected capacitance!");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  auto &V_gf = post_op.GetVGridFunction().Real();
  auto &D_gf = post_op.GetDomainPostOp().D;
  Vector defect;
  auto CorrectedMatrix = [&](const std::vector<Vector> &fields)
  {
    mfem::DenseMatrix C(n);
    for (int i = 0; i < n; i++)
    {
      V_gf.SetFromTrueDofs(fields[i]);
      post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
      response.FixedTraceDomainDefectMult(fields[i], defect);
      // Upper triangle, mirrored: the form is symmetric.
      for (int j = i; j < n; j++)
      {
        if (j > i)
        {
          V_gf.SetFromTrueDofs(fields[j]);
        }
        C(i, j) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf) +
                  linalg::Dot<Vector>(post_op.GetComm(), fields[j], defect);
        C(j, i) = C(i, j);
      }
    }
    return C;
  };
  auto Unavailable = [&]()
  {
    mfem::DenseMatrix C(n);
    C = nan;
    return C;
  };
  auto MutualMatrix = [&](const mfem::DenseMatrix &C)
  {
    mfem::DenseMatrix Cm(n);
    for (int i = 0; i < n; i++)
    {
      Cm(i, i) = 0.0;
      for (int j = 0; j < n; j++)
      {
        if (j != i)
        {
          Cm(i, j) = -C(i, j);
          Cm(i, i) += C(i, j);
        }
      }
      Cm(i, i) += C(i, i);
    }
    return Cm;
  };
  auto InverseMatrix = [&](const mfem::DenseMatrix &C, bool available)
  {
    if (!available)
    {
      return Unavailable();
    }
    mfem::DenseMatrix Cinv(C);
    Cinv.Invert();
    return Cinv;
  };
  const mfem::DenseMatrix C_ft = fixed_trace ? CorrectedMatrix(V) : Unavailable();
  const mfem::DenseMatrix C_sc =
      self_consistent ? CorrectedMatrix(V_corrected) : Unavailable();
  const mfem::DenseMatrix Cm_ft = MutualMatrix(C_ft), Cm_sc = MutualMatrix(C_sc);
  const mfem::DenseMatrix Cinv_ft = InverseMatrix(C_ft, fixed_trace),
                          Cinv_sc = InverseMatrix(C_sc, self_consistent);

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using VT = Units::ValueType;
  // Lane L4 (decision 1829; ERRATUM-1829 E-1; DESIGN-L4-v2 s2.3 E2): under Coverage.Policy
  // "MarkNotApplied" (`source_marked` given) C(i, j) is marked iff source i or source j is
  // (any marked cell in its final measurement); C_m: off-diagonal as C(i, j), diagonal iff
  // any entry of row i of C is; C⁻¹: every entry iff any C entry is. A marked entry prints
  // nan, its value in the appended partial column.
  using MarkMatrix = std::vector<std::vector<bool>>;
  MFEM_VERIFY(!source_marked || static_cast<int>(source_marked->size()) == n,
              "The corrected-matrix marks need one entry per terminal source!");
  const auto marks = ComputeCoverageTerminalMatrixMarks(
      source_marked ? *source_marked : std::vector<bool>(n, false));
  const MarkMatrix &marked_C = marks.capacitance, &marked_Cm = marks.mutual,
                   &marked_Cinv = marks.inverse;
  auto PrintMatrices = [&terminal_sources, source_marked,
                        this](const std::string &file, const std::string &name,
                              const std::string &unit,
                              const mfem::DenseMatrix &fixed_trace_mat,
                              const mfem::DenseMatrix &self_consistent_mat, double scale,
                              const MarkMatrix &marks)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    for (const auto &[idx2, data2] : terminal_sources)
    {
      output.table["i"] << idx2;
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto AppendColumns =
        [&](const char *variant, const mfem::DenseMatrix &mat, bool partial)
    {
      int j = 0;
      for (const auto &[idx2, data2] : terminal_sources)
      {
        const std::string key =
            fmt::format("{}_i2{}{}", variant, idx2, partial ? "_partial" : "");
        output.table.insert(key, fmt::format("{} {}[i][{}] {}{}", name, variant, idx2, unit,
                                             partial ? kCoveragePartialColumnSuffix : ""));
        auto &col = output.table[key];
        for (std::size_t i = 0; i < terminal_sources.size(); i++)
        {
          if (!source_marked)
          {
            col << mat(i, j) * scale;
          }
          else
          {
            col << (marks[i][j] != partial ? nan : mat(i, j) * scale);
          }
        }
        j++;
      }
    };
    AppendColumns("fixed-trace", fixed_trace_mat, false);
    AppendColumns("corrected", self_consistent_mat, false);
    if (source_marked)
    {
      AppendColumns("fixed-trace", fixed_trace_mat, true);
      AppendColumns("corrected", self_consistent_mat, true);
      // Lane L4 round 2 (decision 1914): the explicit mark of entry (i, j), last.
      int j = 0;
      for (const auto &[idx2, data2] : terminal_sources)
      {
        const std::string key = fmt::format("coverage_mark_i2{}", idx2);
        output.table.insert(
            key, fmt::format(fmt::runtime(kCoverageMarkMatrixEntryColumn), idx2));
        auto &col = output.table[key];
        col.print_as_int = true;
        for (std::size_t i = 0; i < terminal_sources.size(); i++)
        {
          col << (marks[i][j] ? 1.0 : 0.0);
        }
        j++;
      }
    }
    output.WriteFullTableTrunc();
  };
  const double F = iodata.units.Dimensionalize<VT::CAPACITANCE>(1.0);
  PrintMatrices("terminal-C-corrected.csv", "C", "(F)", C_ft, C_sc, F, marked_C);
  PrintMatrices("terminal-Cinv-corrected.csv", "C⁻¹", "(1/F)", Cinv_ft, Cinv_sc, 1.0 / F,
                marked_Cinv);
  PrintMatrices("terminal-Cm-corrected.csv", "C_m", "(F)", Cm_ft, Cm_sc, F, marked_Cm);

  // The same matrices in palace.json (SurfaceResponse.TerminalCapacitance; an unavailable
  // variant is null).
  auto MatrixRecord = [&](const mfem::DenseMatrix &mat, bool available)
  {
    nlohmann::json record;
    if (!available)
    {
      return record;
    }
    record = nlohmann::json::array();
    for (int i = 0; i < n; i++)
    {
      auto row = nlohmann::json::array();
      for (int j = 0; j < n; j++)
      {
        row.push_back(mat(i, j) * F);
      }
      record.push_back(std::move(row));
    }
    return record;
  };
  nlohmann::json record;
  record["Indices"] = nlohmann::json::array();
  for (const auto &[idx, data] : terminal_sources)
  {
    record["Indices"].push_back(idx);
  }
  record["FixedTrace"] = MatrixRecord(C_ft, fixed_trace);
  record["SelfConsistent"] = MatrixRecord(C_sc, self_consistent);
  if (source_marked)
  {
    // Lane L4 (E2): a marked entry is null; its value in the Partial matrix (null where
    // not marked).
    auto Split = [&](const char *key, const char *partial_key)
    {
      nlohmann::json partial = record[key];
      if (!record[key].is_array())
      {
        record[partial_key] = nullptr;
        return;
      }
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          (marked_C[i][j] ? record[key][i][j] : partial[i][j]) = nullptr;
        }
      }
      record[partial_key] = std::move(partial);
    };
    Split("FixedTrace", "FixedTracePartial");
    Split("SelfConsistent", "SelfConsistentPartial");
  }
  SaveSurfaceResponseMetadata("TerminalCapacitance", record);
}

}  // namespace palace
