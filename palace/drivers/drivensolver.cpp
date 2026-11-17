// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <optional>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
#include "drivers/drivensubstructuredsweep.hpp"
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/condensedenvironment.hpp"
#include "models/drivensubstructure.hpp"
#include "models/floquetportoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
DrivenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Set up the spatial discretization and frequency sweep.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  const auto &port_excitations = space_op.GetPortExcitations();
  SaveMetadata(port_excitations);

  const auto &omega_sample = iodata.solver.driven.sample_f;

  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  if (adaptive && omega_sample.size() <= iodata.solver.driven.prom_indices.size() &&
      !iodata.solver.driven.adaptive_circuit_synthesis)
  {
    Mpi::Warning("Adaptive frequency sweep requires > {} total frequency samples!\n"
                 "Reverting to uniform sweep!\n",
                 iodata.solver.driven.prom_indices.size());
    adaptive = false;
  }
  SaveMetadata(space_op.GetNDSpaces());
  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             port_excitations.FmtLog());

  std::size_t restart = iodata.solver.driven.restart;
  if (restart != 1)
  {
    std::size_t max_iter = omega_sample.size() * space_op.GetPortExcitations().Size();
    MFEM_VERIFY(
        restart - 1 < max_iter,
        fmt::format("\"Restart\" ({}) is greater than the number of total samples ({})!",
                    restart, max_iter));

    Mpi::Print("\nRestarting from solve {}", iodata.solver.driven.restart);
  }

  // Main frequency sweep loop.
  if (iodata.solver.substructuring && CondensedEnvironment::Enabled(iodata.solver))
  {
    // Online driven substructuring by the regular sweep, with the environment condensed
    // onto the interface in SpaceOperator (prototype).
    CondensedEnvironment env(iodata, space_op);
    space_op.SetCondensedEnvironment(&env);
    auto indicator = SweepUniform(space_op);
    space_op.SetCondensedEnvironment(nullptr);
    return {indicator, space_op.GlobalTrueVSize()};
  }
  if (iodata.solver.substructuring)
  {
    return {SweepSubstructured(iodata, space_op, adaptive), space_op.GlobalTrueVSize()};
  }
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

void DrivenSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                              MPI_Comm comm) const
{
  if (const auto &sub = iodata.solver.substructuring;
      sub && sub->mode == SubstructuringMode::ONLINE)
  {
    DrivenSubstructureModel model;
    model.ReadHeader(sub->save_model, comm);
    MFEM_VERIFY(
        !(iodata.model.Lc > 0.0) ||
            std::abs(iodata.model.Lc - model.length_scale) <= 1.0e-12 * model.length_scale,
        "\"Model\": \"Lc\" (" << iodata.model.Lc
                              << ") differs from that of the saved substructuring model ("
                              << model.length_scale << ")!");
    iodata.model.Lc = model.length_scale;
  }
  BaseSolver::Preprocess(iodata, smesh, comm);
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  // With a condensed environment (online driven substructuring), the postprocessing is the
  // region's, and the environment's port voltages come from it.
  auto *env = space_op.GetCondensedEnvironment();
  PostOperator<ProblemType::DRIVEN> post_op(
      iodata.problem, iodata.solver, env ? env->GetDomains() : iodata.domains,
      env ? env->GetBoundaries() : iodata.boundaries, iodata.units, space_op);

  // The fixed system matrices are only needed when the system cannot be applied through
  // the preconditioner (see SpaceOperator::GetSystemOperator), so assemble them on first
  // use. PEC boundaries are handled simply by setting diagonal entries of the system matrix
  // for the corresponding dofs. Because the Dirichlet BC is always homogeneous, no special
  // elimination is required on the RHS.
  std::unique_ptr<ComplexOperator> K, C, M;
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  // With a condensed environment, the preconditioner is its coupled factorization of the
  // region (refactored at each frequency), which needs no preconditioner matrix.
  auto env_pc = env ? env->GetPreconditioner() : nullptr;
  const bool assemble_pc = !env_pc;
  auto ksp =
      env_pc
          ? std::make_unique<ComplexKspSolver>(iodata.solver.linear, iodata.problem.verbose,
                                               space_op.GetComm(), std::move(env_pc))
          : std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                               &space_op.GetH1Spaces());

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation (none
  // with a condensed environment: the field is zero in the environment).
  std::optional<FluxErrorEstimator<ComplexVector>> estimator;
  if (!env)
  {
    estimator.emplace(space_op, iodata);
  }
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Main excitation and frequency loop: the solves by excitation, then frequency, or by
  // frequency, then excitation, when the preconditioner is a sparse direct solve of the
  // system (on a single level, or the coupled factorization of a condensed environment),
  // which is then factored once per frequency for all excitations. Not with field output,
  // which is collected per excitation, nor with a restart, which counts the solves by
  // excitation. PALACE_DRIVEN_EXCITATION_MAJOR=1 in the environment keeps the order by
  // excitation (prototype, for comparison).
  const auto &linear = iodata.solver.linear;
  const char *excitation_major = std::getenv("PALACE_DRIVEN_EXCITATION_MAJOR");
  const bool frequency_major =
      port_excitations.Size() > 1 && iodata.solver.driven.restart == 1 &&
      !post_op.WillWriteFields() &&
      (!assemble_pc ||
       (space_op.GetNDSpaces().GetNumLevels() == 1 &&
        (linear.type == LinearSolver::MUMPS || linear.type == LinearSolver::STRUMPACK ||
         linear.type == LinearSolver::STRUMPACK_MP ||
         linear.type == LinearSolver::SUPERLU || linear.type == LinearSolver::CUDSS))) &&
      !(excitation_major && std::string(excitation_major) == "1");
  struct Solve
  {
    std::size_t excitation_counter;  // 1-based, in the order of the excitations
    int excitation_idx;
    std::size_t omega_i;
  };
  std::vector<Solve> solves;
  {
    std::size_t excitation_counter = 0;
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      excitation_counter++;
      for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
      {
        if ((excitation_counter - 1) * omega_sample.size() + omega_i + 1 >=
            iodata.solver.driven.restart)
        {
          solves.push_back({excitation_counter, excitation_idx, omega_i});
        }
      }
    }
    if (frequency_major)
    {
      std::ranges::stable_sort(solves, {}, &Solve::omega_i);
      Mpi::Print("\nSweeping all excitations at each frequency\n");
    }
  }
  auto t0 = Timer::Now();
  // The operators at the frequency of the last solve. A refers to P and A2, so both must
  // outlive it.
  std::unique_ptr<ComplexOperator> P, A2, A;
  std::size_t omega_A = omega_sample.size();
  std::vector<int> paraview_excitations;
  for (std::size_t s = 0; s < solves.size(); s++)
  {
    const auto [excitation_counter, excitation_idx, omega_i] = solves[s];
    const bool new_excitation = (s == 0 || solves[s - 1].excitation_idx != excitation_idx);
    if (!frequency_major && new_excitation && port_excitations.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 excitation_counter, port_excitations.Size());
    }
    auto omega = omega_sample[omega_i];
    if (omega_i != omega_A)
    {
      // Assemble frequency dependent matrices and initialize operators in linear
      // solver.
      A.reset();
      A2.reset();
      if (assemble_pc)
      {
        P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                              -omega * omega + 0.0i, omega);
        A = space_op.GetSystemOperator(omega, *P);
      }
      if (!A)
      {
        if (!K)
        {
          K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
          C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
          M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
        }
        A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
        A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                     C.get(), M.get(), A2.get());
      }
      ksp->SetOperators(*A, P ? *P : *A);
      omega_A = omega_i;
    }

    Mpi::Print(
        "\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s{})\n",
        omega_i + 1, omega_sample.size(),
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
            (2 * std::numbers::pi),
        Timer::Duration(Timer::Now() - t0).count(),
        (port_excitations.Size() > 1)
            ? fmt::format(", solve {:d}/{:d}{}",
                          1 + (frequency_major
                                   ? omega_i * port_excitations.Size() +
                                         (excitation_counter - 1)
                                   : omega_i + (excitation_counter - 1) * omega_sample.size()),
                          omega_sample.size() * port_excitations.Size(),
                          frequency_major ? fmt::format(", excitation {:d}", excitation_idx)
                                          : "")
            : "");

    // Solve linear system.
    space_op.GetExcitationVector(excitation_idx, omega, RHS);

    Mpi::Print("\n");
    ksp->Mult(RHS, E);

    if (s + 1 == solves.size())
    {
      // That was the last solve: save the solver totals and release the linear solver
      // and system operators, so they do not add to the memory peak of the field output
      // and error estimation which follow.
      SaveMetadata(*ksp);
      ksp.reset();
      A.reset();
      A2.reset();
      P.reset();
      K.reset();
      C.reset();
      M.reset();
    }

    // Start Post-processing.
    BlockTimer bt0(Timer::POSTPRO);
    Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(space_op.GetComm(), E),
               linalg::Norml2(space_op.GetComm(), RHS));

    // Compute B = -1/(iω) ∇ x E on the true dofs.
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }

    if (std::ranges::find(paraview_excitations, excitation_idx) ==
        paraview_excitations.end())
    {
      // Switch ParaView subfolders once per excitation, after the first solve. Delay the
      // potentially large visualization operators until the first solution is ready, so
      // they do not increase the solve/preconditioner peak.
      post_op.InitializeParaviewDataCollection(excitation_idx);
      paraview_excitations.push_back(excitation_idx);
    }
    if (env)
    {
      post_op.SetLumpedPortVoltages(env->PortVoltages(excitation_idx, omega, E));
    }
    auto total_domain_energy =
        post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

    // Calculate and record the error indicators.
    if (estimator)
    {
      Mpi::Print(" Updating solution error estimates\n");
      estimator->AddErrorIndicator(E, B, total_domain_energy, indicator);
    }

    // Final postprocessing & printing of an excitation. After the last solve, the solver
    // totals were already saved when the solver was released.
    if (ksp && !frequency_major && solves[s + 1].excitation_idx != excitation_idx)
    {
      SaveMetadata(*ksp);
    }
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  std::size_t convergence_memory = iodata.solver.driven.adaptive_memory;
  std::size_t max_size_per_excitation = iodata.solver.driven.adaptive_max_size;
  std::size_t nprom_indices = iodata.solver.driven.prom_indices.size();
  MFEM_VERIFY(max_size_per_excitation <= 0 || max_size_per_excitation >= nprom_indices,
              "Adaptive frequency sweep must sample at least " << nprom_indices
                                                               << " frequency points!");

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector E(Curl.Width()), Eh(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  Eh.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  Eh = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  const FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Configure the PROM operator which performs the parameter space sampling and basis
  // construction during the offline phase as well as the PROM solution during the online
  // phase.
  auto t0 = Timer::Now();
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                          (2 * std::numbers::pi);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             omega_sample.size(), omega_sample.front() * unit_GHz,
             omega_sample.back() * unit_GHz);
  RomOperator prom_op(iodata, space_op, max_size_per_excitation);
  space_op.GetWavePortOp().SetSuppressOutput(true);
  space_op.GetWavePortOp().ConfigureReducedModelTraining(
      max_size_per_excitation, port_excitations.Size(),
      iodata.solver.driven.adaptive_circuit_synthesis);

  // Add ports to PROM if we do synthesis.
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.AddLumpedPortModesForSynthesis();
    if (space_op.GetWavePortOp().Size() > 0)
    {
      // Use the band center as the reference frequency for seeding wave-port modes.
      // The choice rescales the basis vector but does not change correctness.
      const double omega_ref = 0.5 * (omega_sample.front() + omega_sample.back());
      prom_op.AddWavePortModesForSynthesis(omega_ref);
    }
  }

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](int excitation_idx, double omega, std::size_t sample_idx)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(E, fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx));
    prom_op.UpdateMRI(excitation_idx, omega, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for energy postprocessing and error estimation.
    BlockTimer bt0(Timer::POSTPRO);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }

    // Measure domain energies for the error indicator only. Don't exchange face_nbr_data,
    // unless printing paraview fields.
    auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(E, B);
    estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
  };

  // Loop excitations to add to PROM.
  //
  // Restart should not really be used for adaptive sweeps, but must work. Construct PROM in
  // the same way same regardless of restart for consistency. Don't shift excitation start.
  int excitation_counter = 0;
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nAdding excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++excitation_counter, port_excitations.Size());
    }
    prom_op.SetExcitationIndex(excitation_idx);  // Pre-compute RHS1

    // Initialize PROM with explicit HDM samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    std::size_t counter_rom_sample = 0;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      auto omega = omega_sample[i];
      prom_op.SolveHDM(excitation_idx, omega, E);
      prom_op.SolvePROM(excitation_idx, omega, Eh);
      linalg::AXPY(-1.0, E, Eh);
      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      UpdatePROM(excitation_idx, omega, counter_rom_sample);
      counter_rom_sample++;
    }
    // The estimates associated to the end points are assumed inaccurate.
    max_errors[0] = std::numeric_limits<double>::infinity();
    max_errors[1] = std::numeric_limits<double>::infinity();
    auto memory = std::distance(max_errors.rbegin(),
                                std::find_if(max_errors.rbegin(), max_errors.rend(),
                                             [=](auto x) { return x > offline_tol; }));
    memory = std::max(0L, memory);  // Ensure memory >= 0 as it should be.

    // Greedy procedure for basis construction (offline phase). Basis is initialized with
    // solutions at frequency sweep endpoints and explicit sample frequencies.
    std::size_t it = max_errors.size();
    for (std::size_t it0 = it; it < max_size_per_excitation && memory < convergence_memory;
         it++)
    {
      // Compute the location of the maximum error in parameter domain (bounded by the
      // previous samples).
      double omega_star = prom_op.FindMaxError(excitation_idx)[0];

      // Sample HDM and add solution to basis.
      prom_op.SolveHDM(excitation_idx, omega_star, E);
      prom_op.SolvePROM(excitation_idx, omega_star, Eh);
      linalg::AXPY(-1.0, E, Eh);

      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      memory = max_errors.back() < offline_tol ? memory + 1 : 0;

      Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
                 "{:.3e}, memory = {:d}/{:d}\n",
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(), memory, convergence_memory);
      UpdatePROM(excitation_idx, omega_star, counter_rom_sample);
      counter_rom_sample++;
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size_per_excitation) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
  }

  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root

  // Circuit synthesis samples the exact cross-section EVP at a tightened tolerance so the
  // fitted pencil does not depend on the port-mode accuracy floor; keep those solves exact
  // (they also enrich the per-port reduced basis).
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.PrintPROMMatrices(iodata.units, iodata.problem.output);
  }

  // Exact wave-port modes computed by the HDM samples (and any synthesis samples) have
  // trained a separate per-port reduced eigenspace. Enable guarded Rayleigh-Ritz evaluation
  // for the online output sweep; failed residual checks transparently fall back to the
  // exact port eigensolver and enrich the basis.
  space_op.GetWavePortOp().EnableReducedModel(iodata.solver.driven.adaptive_tol);
  prom_op.PrepareOnlineExcitations();
  post_op.ConfigureReducedPostprocessing(prom_op);

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  auto solve_online_point = [&](int excitation_idx, std::size_t omega_i)
  {
    const auto omega = omega_sample[omega_i];
    // Refresh every port, including inactive observation-only ports, before either PROM
    // assembly or postprocessing. Repeated calls at the same frequency are cached.
    space_op.GetWavePortOp().PrepareFrequency(omega);
    Mpi::Print("\nIt {:d}/{:d}, excitation {:d}: ω/2π = {:.3e} GHz "
               "(total elapsed time = {:.2e} s)\n",
               omega_i + 1, omega_sample.size(), excitation_idx,
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                   (2 * std::numbers::pi),
               Timer::Duration(Timer::Now() - t0).count());

    // Assemble and solve the PROM linear system.
    prom_op.SolvePROM(excitation_idx, omega, E);
    Mpi::Print("\n");
    if (omega_i == 0 && post_op.WillWriteFields())
    {
      // Switch ParaView subfolders once per excitation. Delay visualization setup until
      // the first online solution is ready.
      post_op.InitializeParaviewDataCollection(excitation_idx);
    }

    if (post_op.HasReducedPostprocessing() &&
        !post_op.WillWriteFields(static_cast<int>(omega_i)))
    {
      post_op.MeasureAndPrintReduced(excitation_idx, int(omega_i), E, omega,
                                     prom_op.GetReducedSolution());
      return;
    }

    // Start full post-processing.
    BlockTimer bt0(Timer::POSTPRO);
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));

    // Compute B = -1/(iω) ∇ x E on the true dofs.
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }
    post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);
  };

  if (!post_op.WillWriteFields())
  {
    // Port modes depend on frequency, not excitation. Frequency-major traversal keeps the
    // one-frequency modal cache hot for every excitation.
    for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
    {
      for (const auto &[excitation_idx, excitation_spec] : port_excitations)
      {
        solve_online_point(excitation_idx, omega_i);
      }
    }
  }
  else
  {
    // Preserve excitation-major ordering for excitation-specific field collections.
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      if (port_excitations.Size() > 1)
      {
        Mpi::Print("\nSweeping excitation index {:d}:\n", excitation_idx);
      }
      for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
      {
        solve_online_point(excitation_idx, omega_i);
      }
    }
  }

  // Final postprocessing & printing: no change to indicator since these are in PROM.
  {
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(prom_op.GetLinearSolver());
  }
  space_op.GetWavePortOp().PrintReducedModelStats();
  post_op.MeasureFinalize(indicator);
  return indicator;
}

}  // namespace palace
