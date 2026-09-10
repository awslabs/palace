// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
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

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
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
                          ArchivedField field, MPI_Comm comm)
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
  MFEM_VERIFY(
      stream && magic == response_archive_magic && version == response_archive_version &&
          stored_field == static_cast<std::uint32_t>(field) && stored_source == source &&
          stored_rank == rank && stored_size == size && local_size >= 0,
      "Invalid response archive header in \"" << path.string() << "\"!");
  Vector vector(static_cast<int>(local_size));
  stream.read(reinterpret_cast<char *>(vector.HostWrite()), local_size * sizeof(double));
  MFEM_VERIFY(stream && stream.peek() == std::ifstream::traits_type::eof(),
              "Invalid response archive payload in \"" << path.string() << "\"!");
  vector.UseDevice(true);
  return vector;
}

}  // namespace

std::pair<ErrorIndicator, long long int>
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  LaplaceOperator laplace_op(iodata, mesh);
  auto K = laplace_op.GetStiffnessMatrix();
  const auto *response_config = iodata.solver.electrostatic.response_correction
                                    ? &*iodata.solver.electrostatic.response_correction
                                    : nullptr;
  const bool postprocess_response =
      response_config && response_config->IncludesPostprocessing();
  const bool self_consistent_response =
      response_config && response_config->IncludesSelfConsistent();
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
  ksp.SetOperators(*K, *K);

  // Source indices are either equipotential terminals or prescribed potential traces.
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op, nullptr,
                                                   &surface_post_geometry);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0,
              "No terminal or prescribed potential boundaries specified for electrostatic "
              "simulation!");
  const auto response_archive = ResponseArchiveDirectory();
  const bool archive_reduce_only = EnvironmentFlag("PALACE_RESPONSE_REDUCE_ONLY");
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

  auto &V_gf = post_op.GetVGridFunction().Real();
  auto &D_gf = post_op.GetDomainPostOp().D;
  const std::size_t block = static_cast<std::size_t>(block_size);
  const std::size_t block_count = (basis_size + block - 1) / block;
  std::size_t completed = 0;
  const std::size_t total_blocks = block_count * (block_count + 1) / 2;
  for (std::size_t first_block = 0; first_block < block_count; first_block++)
  {
    const std::size_t first_begin = first_block * block;
    const std::size_t first_end = std::min(first_begin + block, basis_size);
    const std::size_t first_count = first_end - first_begin;
    for (std::size_t second_block = first_block; second_block < block_count; second_block++)
    {
      const std::size_t second_begin = second_block * block;
      const std::size_t second_end = std::min(second_begin + block, basis_size);
      const std::size_t second_count = second_end - second_begin;
      const bool diagonal_block = first_block == second_block;

      std::vector<std::size_t> global_indices;
      global_indices.reserve(first_count + (diagonal_block ? 0 : second_count));
      for (std::size_t i = first_begin; i < first_end; i++)
      {
        global_indices.push_back(i);
      }
      if (!diagonal_block)
      {
        for (std::size_t i = second_begin; i < second_end; i++)
        {
          global_indices.push_back(i);
        }
      }

      std::vector<Vector> V, D, E_basis, V_local;
      V.reserve(global_indices.size());
      if (archive_has_flux)
      {
        D.reserve(global_indices.size());
      }
      E_basis.resize(global_indices.size());
      V_local.resize(global_indices.size());
      for (std::size_t local = 0; local < global_indices.size(); local++)
      {
        const int source = basis_indices[global_indices[local]];
        V.push_back(ReadArchivedVector(archive, source, ArchivedField::POTENTIAL,
                                       laplace_op.GetComm()));
        if (archive_has_flux)
        {
          D.push_back(ReadArchivedVector(archive, source, ArchivedField::FLUX,
                                         laplace_op.GetComm()));
        }
        MFEM_VERIFY(V.back().Size() == Grad.Width(),
                    "Archived potential field has an incompatible local size!");
        E_basis[local].SetSize(Grad.Height());
        E_basis[local] = 0.0;
        Grad.AddMult(V.back(), E_basis[local], -1.0);
        V_gf.SetFromTrueDofs(V.back());
        V_local[local] = V_gf;
      }

      const auto surface_matrices = post_op.GetInterfaceElectricFieldEnergyMatrices(
          E_basis, archive_has_flux ? &D : nullptr);
      mfem::DenseMatrix local_domain(static_cast<int>(first_count),
                                     static_cast<int>(second_count));
      local_domain = 0.0;
      for (std::size_t local_i = 0; local_i < first_count; local_i++)
      {
        post_op.GetDomainPostOp().M_elec->Mult(V_local[local_i], D_gf);
        for (std::size_t local_j = 0; local_j < second_count; local_j++)
        {
          const std::size_t union_j = diagonal_block ? local_j : first_count + local_j;
          if (diagonal_block && local_j < local_i)
          {
            continue;
          }
          local_domain(static_cast<int>(local_i), static_cast<int>(local_j)) =
              0.5 * linalg::LocalDot(V_local[union_j], D_gf);
        }
      }
      Mpi::GlobalSum(local_domain.Height() * local_domain.Width(), local_domain.GetData(),
                     post_op.GetComm());

      for (std::size_t local_i = 0; local_i < first_count; local_i++)
      {
        const std::size_t global_i = first_begin + local_i;
        for (std::size_t local_j = 0; local_j < second_count; local_j++)
        {
          if (diagonal_block && local_j < local_i)
          {
            continue;
          }
          const std::size_t global_j = second_begin + local_j;
          const int union_i = static_cast<int>(local_i);
          const int union_j =
              static_cast<int>(diagonal_block ? local_j : first_count + local_j);
          for (const auto &[interface, entries] : surface_matrices)
          {
            for (const auto &entry : entries)
            {
              AppendSurface(interface, entry.distance, global_i, global_j, entry, union_i,
                            union_j);
            }
          }
          if (root)
          {
            domain_output.table["basis_i"] << basis_indices[global_i];
            domain_output.table["basis_j"] << basis_indices[global_j];
            domain_output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(
                local_domain(static_cast<int>(local_i), static_cast<int>(local_j)));
          }
        }
      }
      completed++;
      Mpi::Print(" Archived response block pair {:d}/{:d}\n", static_cast<int>(completed),
                 static_cast<int>(total_blocks));
    }
  }
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
