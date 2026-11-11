// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <numeric>
#include <optional>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
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
#include "utils/perfstats.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The fluxes of nearby frequencies are close, so each flux recovery of the error estimator
// starts from the closest combination of the last ones (which, like the PROM basis, span
// the solutions of a sweep well).
constexpr int estimator_num_prev = 20;

// Residual of the eigenvector of a pole of the PROM in the HDM equations, relative to the
// size of their terms, below which the pole is a resonance of the HDM for the error
// indicator of "AdaptiveIndicatorMaxQ": the eigenvectors of the spurious poles of the
// projection leave residuals of order one (above 1e-2 in every case measured), those of the
// resonances which the sampling resolves 1e-6 to 3e-3. A fixed value: the classification
// does not depend on "AdaptiveTol".
constexpr double resonance_residual_tol = 1.0e-2;

}  // namespace

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
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // When the preconditioner's finest level is the sparse part of the system operator, the
  // sweep applies the system through it (see SpaceOperator::GetSystemOperator) and never
  // needs K, C or M. Otherwise, assemble them on first use. PEC boundaries are handled by
  // setting the system-matrix diagonal. The Dirichlet condition is homogeneous, so no
  // special RHS elimination is required.
  std::unique_ptr<ComplexOperator> K, C, M;
  bool fixed_system_operators_assembled = false;
  auto AssembleFixedSystemOperators = [&]()
  {
    if (!fixed_system_operators_assembled)
    {
      K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
      C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      fixed_system_operators_assembled = true;
    }
  };
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  CoarseReuse coarse_reuse(iodata.solver.driven.coarse_reuse_tol);

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  std::optional<FluxErrorEstimator<ComplexVector>> estimator;
  if (ComputeErrorEstimation())
  {
    estimator.emplace(space_op, iodata, estimator_num_prev);
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

  // Main excitation and frequency loop. Per-frequency assembly is reported as a tally of
  // the Operator Construction row.
  const int assembly_tally =
      perf::RegisterTally("Driven/FrequencyAssembly", Timer::CONSTRUCT);
  auto t0 = Timer::Now();
  std::size_t excitation_counter = 0;
  const std::size_t excitation_restart_counter =
      ((iodata.solver.driven.restart - 1) / omega_sample.size()) + 1;
  const std::size_t freq_restart_idx =
      (iodata.solver.driven.restart - 1) % omega_sample.size();
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (++excitation_counter < excitation_restart_counter)
    {
      continue;
    }
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 excitation_counter, port_excitations.Size());
    }
    // Frequency loop. Delay the potentially large visualization operators until the
    // first solution is ready, so they do not increase the solve/preconditioner peak.
    const std::size_t omega_start =
        (excitation_counter == excitation_restart_counter) ? freq_restart_idx : 0;
    for (std::size_t omega_i = omega_start; omega_i < omega_sample.size(); omega_i++)
    {
      auto omega = omega_sample[omega_i];
      // Assemble frequency dependent matrices and initialize operators in linear
      // solver. A refers to P, and to A2 when assembled, so both must outlive it.
      std::unique_ptr<ComplexOperator> A2, A, P;
      {
        perf::ScopedTally st(assembly_tally);
        P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                              -omega * omega + 0.0i, omega);
        A = space_op.GetSystemOperator(omega, *P);
        if (!A)
        {
          AssembleFixedSystemOperators();
          A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
          A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                       K.get(), C.get(), M.get(), A2.get());
        }
      }

      Mpi::Print(
          "\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s{})\n",
          omega_i + 1, omega_sample.size(),
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
              (2 * std::numbers::pi),
          Timer::Duration(Timer::Now() - t0).count(),
          (port_excitations.Size() > 1)
              ? fmt::format(", solve {:d}/{:d}",
                            1 + omega_i + (excitation_counter - 1) * omega_sample.size(),
                            omega_sample.size() * port_excitations.Size())
              : "");

      // Solve linear system.
      space_op.GetExcitationVector(excitation_idx, omega, RHS);

      Mpi::Print("\n");
      coarse_reuse.Solve(ksp, *A, *P, omega, RHS, E, iodata.solver.linear.initial_guess);

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
            E, B,
            space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
      }

      if (omega_i == omega_start)
      {
        // Switch ParaView subfolders once per excitation, after the first solve.
        post_op.InitializeParaviewDataCollection(excitation_idx);
      }
      auto total_domain_energy =
          post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

      // Calculate and record the error indicators.
      if (estimator)
      {
        Mpi::Print(" Updating solution error estimates\n");
        estimator->AddErrorIndicator(E, B, total_domain_energy, indicator);
      }
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
  const double pole_width = iodata.solver.driven.adaptive_pole_exclusion;
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
  std::optional<FluxErrorEstimator<ComplexVector>> estimator;
  if (ComputeErrorEstimation())
  {
    estimator.emplace(space_op, iodata, estimator_num_prev);
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
  // The error indicator for adaptive mesh refinement accumulates over the HDM samples, or,
  // with "AdaptiveIndicatorSamples", over the PROM solutions at fixed frequencies of the
  // sweep after the sampling (within AdaptiveTol of the HDM): the marking then does not
  // depend on the sample set chosen by the greedy search, which follows round-off.
  const auto &indicator_indices = iodata.solver.driven.indicator_indices;
  auto AddErrorIndicator = [&](double omega, double weight)
  {
    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for the domain energy normalizing the error indicator.
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
    if (weight == 1.0)
    {
      estimator->AddErrorIndicator(E, B, total_domain_energy, indicator);
    }
    else
    {
      ErrorIndicator field_indicator;
      estimator->AddErrorIndicator(E, B, total_domain_energy, field_indicator);
      indicator.AddIndicator(field_indicator.Local(), weight);
    }
  };
  auto UpdatePROM = [&](int excitation_idx, double omega, std::size_t sample_idx)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(E, fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx));
    prom_op.UpdateMRI(excitation_idx, omega, E);
    if (estimator && indicator_indices.empty())
    {
      AddErrorIndicator(omega, 1.0);
    }
  };

  // With initial guesses enabled, each HDM sample solve starts from the PROM solution (zero
  // while the reduced basis is empty).
  const bool initial_guess = iodata.solver.linear.initial_guess;

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

    // The error of the PROM within a resonance much narrower than the sampling of the sweep
    // cannot be reduced below the tolerance (the HDM samples do not determine the PROM
    // pole to a small fraction of the linewidth), and the search for the maximum error
    // returns to such resonances, which dominate the indicator |Q(z)|⁻¹ even after the PROM
    // has resolved them. Exclude the poles narrower than the relative width
    // "AdaptivePoleExclusion" from the search, except for the sweep frequencies, once
    // converged to the accuracy which this width and the tolerance require.
    std::vector<std::complex<double>> poles;
    std::vector<std::pair<double, double>> excluded;
    std::vector<double> excluded_sweep;
    auto UpdateExcluded = [&]()
    {
      if (pole_width <= 0.0)
      {
        return;
      }
      auto new_poles = prom_op.ComputePoles();
      excluded = NarrowPoleIntervals(new_poles, poles, omega_sample.front(),
                                     omega_sample.back(), pole_width, offline_tol);
      poles = std::move(new_poles);
      excluded_sweep.clear();
      for (auto omega : omega_sample)
      {
        if (std::ranges::any_of(
                excluded, [=](const auto &interval)
                { return interval.first <= omega && omega <= interval.second; }))
        {
          excluded_sweep.push_back(omega);
        }
      }
    };

    // Initialize PROM with explicit HDM samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    std::size_t counter_rom_sample = 0;
    std::vector<int> hdm_its;
    auto SolveHDM = [&](double omega)
    {
      const int its0 = prom_op.GetLinearSolver().NumTotalMultIterations();
      prom_op.SolveHDM(excitation_idx, omega, E);
      hdm_its.push_back(prom_op.GetLinearSolver().NumTotalMultIterations() - its0);
      perf::SetMax(perf::RegisterCounter("PROM/HDMMaxIts", perf::Reduce::MAX,
                                         perf::Accumulate::MAX),
                   hdm_its.back());
    };

    // Solve the HDM at a sample frequency and return the relative error of the PROM there.
    auto SolveHDMSample = [&](double omega)
    {
      if (initial_guess)
      {
        prom_op.SolvePROM(excitation_idx, omega, Eh);
        E = Eh;
        SolveHDM(omega);
      }
      else
      {
        SolveHDM(omega);
        prom_op.SolvePROM(excitation_idx, omega, Eh);
      }
      linalg::AXPY(-1.0, E, Eh);
      return linalg::Norml2(space_op.GetComm(), Eh) / linalg::Norml2(space_op.GetComm(), E);
    };
    for (auto i : iodata.solver.driven.prom_indices)
    {
      auto omega = omega_sample[i];
      max_errors.push_back(SolveHDMSample(omega));
      UpdatePROM(excitation_idx, omega, counter_rom_sample);
      counter_rom_sample++;
      UpdateExcluded();
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
      double omega_star =
          excluded.empty()
              ? prom_op.FindMaxError(excitation_idx)[0]
              : prom_op.FindMaxError(excitation_idx, 1, excluded, excluded_sweep)[0];

      // Sample HDM and add solution to basis.
      max_errors.push_back(SolveHDMSample(omega_star));
      memory = max_errors.back() < offline_tol ? memory + 1 : 0;

      Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
                 "{:.3e}, memory = {:d}/{:d}\n",
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(), memory, convergence_memory);
      UpdatePROM(excitation_idx, omega_star, counter_rom_sample);
      counter_rom_sample++;
      UpdateExcluded();
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size_per_excitation) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
    if (!excluded.empty())
    {
      std::vector<double> centers;
      for (const auto &[lo, hi] : excluded)
      {
        centers.push_back(0.5 * (lo + hi));
      }
      utils::PrettyPrint(
          centers, unit_GHz,
          fmt::format(" Excluded resonances (GHz, relative width {:.1e}):", pole_width));
    }

    const bool max_reached = (it == max_size_per_excitation);
    perf::Add("PROM/Excitations", 1);
    perf::Add("PROM/Samples", it);
    perf::Add("PROM/GreedyIts", it - iodata.solver.driven.prom_indices.size());
    perf::Add(max_reached ? "PROM/MaxSamplesReached" : "PROM/Converged", 1);
    perf::Add("PROM/HDMSolves", hdm_its.size());
    perf::Add("PROM/HDMIts", std::accumulate(hdm_its.begin(), hdm_its.end(), 0));
    perf::Add("PROM/ExcludedResonances", excluded.size());
    if (std::isfinite(max_errors.back()))
    {
      perf::SetMax(perf::RegisterCounter("PROM/FinalMaxError", perf::Reduce::MAX,
                                         perf::Accumulate::MAX),
                   max_errors.back());
    }
    // Error estimates of the samples after the two end points, in sampling order.
    const auto prefix = fmt::format("PROM/Excitation{:d}/", excitation_idx);
    perf::SetInfo(prefix + "SampleMaxError",
                  std::vector<double>(max_errors.begin() + 2, max_errors.end()));
    perf::SetInfo(prefix + "HDMIts", hdm_its);
  }
  perf::Set("PROM/ReducedDimension", prom_op.GetReducedDimension(), perf::Reduce::MAX,
            perf::Accumulate::MAX);

  if (estimator && !indicator_indices.empty())
  {
    // Error indicator from the final PROM at the fixed frequencies, and at its resonances
    // in the band with a quality factor below "AdaptiveIndicatorMaxQ": the fixed
    // frequencies only weigh a resonance narrower than their spacing when one of them is
    // within a few linewidths. The resonances are the poles whose eigenvectors satisfy the
    // HDM equations to resonance_residual_tol: the spurious poles of the projection change
    // with the samples.
    std::vector<double> omega_ind;
    for (auto i : indicator_indices)
    {
      omega_ind.push_back(omega_sample[i]);
    }
    utils::PrettyPrint(omega_ind, unit_GHz,
                       " Error indicator from the PROM at the frequencies (GHz):");
    const double max_q = iodata.solver.driven.adaptive_indicator_max_q;
    std::optional<std::vector<std::pair<std::complex<double>, double>>> poles;
    std::vector<double> omega_res;
    if (max_q > 0.0)
    {
      poles =
          prom_op.ComputePoleResiduals(omega_sample.front(), omega_sample.back(), max_q);
      if (!poles)
      {
        Mpi::Warning(
            "\"AdaptiveIndicatorMaxQ\" adds no resonances: the poles of the PROM are "
            "not computed with frequency-dependent boundary terms (wave or Floquet "
            "ports, second-order farfield, surface conductivity or rational surface "
            "impedance boundaries)!\n");
      }
    }
    if (poles)
    {
      Mpi::Print(" and at the resonances with Q < {:.3e}:\n", max_q);
      for (const auto &[p, r] : *poles)
      {
        Mpi::Print("   f = {:.6e} GHz, Q = {:.3e}, relative residual = {:.3e}{}\n",
                   p.real() * unit_GHz, std::abs(p) / (2.0 * std::abs(p.imag())), r,
                   (r < resonance_residual_tol) ? "" : " (spurious: not used)");
        if (r < resonance_residual_tol)
        {
          omega_res.push_back(p.real());
        }
      }
      perf::Add("PROM/IndicatorResonances", omega_res.size());
    }
    // Each resonance weighs as much as the fixed frequencies together: one field among
    // theirs would not drive the refinement (it would weigh 1 / (N + 1) of the indicator).
    const double res_weight = static_cast<double>(omega_ind.size());
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      for (auto omega : omega_ind)
      {
        prom_op.SolvePROM(excitation_idx, omega, E);
        AddErrorIndicator(omega, 1.0);
      }
      for (auto omega : omega_res)
      {
        prom_op.SolvePROM(excitation_idx, omega, E);
        AddErrorIndicator(omega, res_weight);
      }
    }
    perf::Add("PROM/IndicatorSamples",
              (omega_ind.size() + omega_res.size()) * port_excitations.Size());
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

  // No change to indicator since these are in PROM.
  space_op.GetWavePortOp().PrintReducedModelStats();
  post_op.MeasureFinalize(indicator);
  return indicator;
}

}  // namespace palace
