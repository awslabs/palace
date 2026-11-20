// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
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
#include "models/surfaceresponseoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

internal::DrivenSystemOperator internal::BuildDrivenSystemOperator(
    SpaceOperator &space_op, double omega, const ComplexOperator &P,
    std::unique_ptr<ComplexOperator> &K, std::unique_ptr<ComplexOperator> &C,
    std::unique_ptr<ComplexOperator> &M, const ComplexOperator *response_mass,
    const ComplexOperator *response_stiffness)
{
  MFEM_VERIFY(!response_stiffness || response_mass,
              "A response-corrected driven stiffness requires a response mass!");
  DrivenSystemOperator system_op;
  system_op.A = space_op.GetSystemOperator(omega, P);
  if (!system_op.A)
  {
    // The fixed system matrices are only needed when the system cannot be applied through
    // the preconditioner, so assemble them on first use.
    system_op.assembled = true;
    if (!K)
    {
      K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
      C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    }
    system_op.A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    system_op.A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                           K.get(), C.get(), M.get(), system_op.A2.get());
  }
  // The response terms are added to whichever system operator A is (applied through the
  // finest level of P, or assembled), so the corrected solve sees the same system.
  if (response_mass)
  {
    system_op.corrected = std::make_unique<SumComplexOperator>(*system_op.A);
    system_op.corrected->AddOperator(*response_mass, -omega * omega + 0.0i);
    if (response_stiffness)
    {
      system_op.corrected->AddOperator(*response_stiffness, 1.0 + 0.0i);
    }
  }
  return system_op;
}

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
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op, &response_geometry,
                                            &surface_post_geometry);

  // The fixed system matrices are only needed when the system cannot be applied through
  // the preconditioner (see SpaceOperator::GetSystemOperator), so assemble them on first
  // use. PEC boundaries are handled simply by setting diagonal entries of the system matrix
  // for the corresponding dofs. Because the Dirichlet BC is always homogeneous, no special
  // elimination is required on the RHS.
  std::unique_ptr<ComplexOperator> K, C, M;
  SurfaceResponseOperator *response_correction = post_op.GetSurfaceResponseOperator();
  if (response_correction && iodata.model.refinement.max_it == 0)
  {
    // Round 2b B3 (decision 1484): no AMR construction follows (Refinement.MaxIts 0).
    ReleaseResponseMatrixCaches();
  }
  const bool self_consistent_response =
      response_correction && !response_correction->IsRawFallback() &&
      iodata.solver.surface_response_correction->IncludesSelfConsistent();
  std::unique_ptr<ComplexWrapperOperator> response_mass;
  // The magnetic (inductive) defect dK (mx-sc-magnetic DESIGN-v2 section 2.1, E1): the
  // K-side twin of the response mass, present only when a loaded model carries a magnetic
  // leaf (a library without them: the operator of record, bitwise).
  std::unique_ptr<SurfaceResponseStiffnessOperator> response_stiffness_real;
  std::unique_ptr<ComplexWrapperOperator> response_stiffness;
  if (self_consistent_response)
  {
    response_mass = std::make_unique<ComplexWrapperOperator>(response_correction, nullptr);
    if (response_correction->HasMagneticModels())
    {
      response_stiffness_real =
          std::make_unique<SurfaceResponseStiffnessOperator>(*response_correction);
      response_stiffness =
          std::make_unique<ComplexWrapperOperator>(response_stiffness_real.get(), nullptr);
    }
  }
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  std::unique_ptr<ComplexKspSolver> corrected_ksp;
  if (self_consistent_response)
  {
    corrected_ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                       &space_op.GetH1Spaces());
    corrected_ksp->SetInitialGuess(true);
    corrected_ksp->SetRelTol(iodata.solver.surface_response_correction->solve_tol);
  }

  // Initialize structures for storing and reducing the results of error estimation.
  const FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), E_corrected, B(Curl.Height()), D,
      D_corrected;
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  if (post_op.NeedsRecoveredElectricFlux())
  {
    // The recovered flux is allocated only when interface postprocessing uses it.
    D.SetSize(estimator.GetElectricFluxSpace().GetTrueVSize());
    D.UseDevice(true);
  }
  E = 0.0;
  B = 0.0;
  if (self_consistent_response)
  {
    E_corrected.SetSize(Curl.Width());
    E_corrected.UseDevice(true);
    E_corrected = 0.0;
    if (post_op.NeedsRecoveredElectricFlux())
    {
      D_corrected.SetSize(estimator.GetElectricFluxSpace().GetTrueVSize());
      D_corrected.UseDevice(true);
    }
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

  // Main excitation and frequency loop.
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
    // Start each excitation from zero rather than the previous excitation's solution.
    E = 0.0;

    // Frequency loop. Delay the potentially large visualization operators until the
    // first solution is ready, so they do not increase the solve/preconditioner peak.
    const std::size_t omega_start =
        (excitation_counter == excitation_restart_counter) ? freq_restart_idx : 0;
    for (std::size_t omega_i = omega_start; omega_i < omega_sample.size(); omega_i++)
    {
      auto omega = omega_sample[omega_i];
      // Assemble frequency dependent matrices and initialize operators in linear
      // solver. The system operator refers to P, so P must outlive it.
      auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
          1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
      auto system_op = internal::BuildDrivenSystemOperator(
          space_op, omega, *P, K, C, M, response_mass.get(), response_stiffness.get());
      ksp->SetOperators(*system_op.A, *P);
      if (corrected_ksp)
      {
        corrected_ksp->SetOperators(*system_op.corrected, *P);
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
      if (corrected_ksp)
      {
        // The path the response-corrected system operator wraps (log only).
        Mpi::Print(" Response-corrected system operator: {}\n",
                   system_op.assembled ? "assembled" : "matrix-free");
      }

      // Solve linear system.
      space_op.GetExcitationVector(excitation_idx, omega, RHS);

      Mpi::Print("\n");
      ksp->Mult(RHS, E);
      if (response_correction && !response_correction->IsRawFallback() &&
          response_correction->HasConductorConsistencyProbes())
      {
        // The conductor-consistency gate on the Maxwell path (idfix DESIGN B.3; decision
        // 1064 (1)): the raw field of this frequency, BEFORE the self-consistent solve (the
        // weights are read at Mult time, so the corrected field sees the exclusion) and
        // before any energy of this frequency is formed; an excluded patch is excluded at
        // this (excitation, frequency) only under Exclusion PerExcitation, at every later
        // one of this operator under Sticky (lane ccpx).
        BlockTimer gate_timer(Timer::POSTPRO_RESPONSE_COUPON);
        response_correction->ApplyConductorConsistencyGate(
            E, excitation_idx,
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI),
            false);
        if (response_correction->HasMagneticModels())
        {
          // The dual R0-a gate on the magnetic term (DESIGN-v2 section 3): the same raw
          // field, before the corrected solve reads the exclusion.
          response_correction->ApplyMagneticConductorConsistencyGate(
              E, excitation_idx,
              iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI));
        }
      }
      if (corrected_ksp)
      {
        Mpi::Print(" Solving self-consistent fabrication-response corrected field\n");
        E_corrected = E;
        corrected_ksp->Mult(RHS, E_corrected);
      }

      if (omega_i + 1 == omega_sample.size() &&
          excitation_counter == port_excitations.Size())
      {
        // That was the last solve (raw and, when self-consistent, corrected): save the
        // solver totals and release the linear solvers and system operators, so they do not
        // add to the memory peak of the field output and error estimation which follow.
        SaveMetadata(*ksp);
        ksp.reset();
        if (corrected_ksp)
        {
          SaveSurfaceResponseSolverMetadata(space_op.GetComm(), "Driven",
                                            corrected_ksp->NumTotalMult(),
                                            corrected_ksp->NumTotalMultIterations());
          corrected_ksp.reset();
        }
        system_op.corrected.reset();
        system_op.A.reset();
        system_op.A2.reset();
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
            E, B,
            space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
      }

      if (omega_i == omega_start)
      {
        // Switch ParaView subfolders once per excitation, after the first solve.
        post_op.InitializeParaviewDataCollection(excitation_idx);
      }
      const ComplexVector *recovered_flux = nullptr;
      if (post_op.NeedsRecoveredElectricFlux())
      {
        Mpi::Print(" Recovering electric flux for interface postprocessing\n");
        estimator.RecoverElectricFlux(E, D);
        post_op.SetRecoveredElectricFlux(D);
        recovered_flux = &D;
      }
      if (self_consistent_response)
      {
        // Not corrected_ksp: it is released after the last solve, before this output.
        const ComplexVector *corrected_flux = nullptr;
        if (post_op.NeedsRecoveredElectricFlux())
        {
          Mpi::Print(" Recovering electric flux for corrected interface postprocessing\n");
          estimator.RecoverElectricFlux(E_corrected, D_corrected);
          corrected_flux = &D_corrected;
        }
        post_op.SetSurfaceResponseCorrectedField(E_corrected, corrected_flux);
      }
      auto total_domain_energy =
          post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

      // Calculate and record the error indicators.
      Mpi::Print(" Updating solution error estimates\n");
      if (recovered_flux)
      {
        estimator.AddErrorIndicator(E, B, *recovered_flux, total_domain_energy, indicator);
      }
      else
      {
        estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
      }
    }

    // Final postprocessing & printing. After the last excitation, the solver totals were
    // already saved when the solver was released.
    BlockTimer bt0(Timer::POSTPRO);
    if (ksp)
    {
      SaveMetadata(*ksp);
    }
  }
  post_op.MeasureFinalize(indicator);
  if (corrected_ksp)
  {
    // Not released (no last solve ran); otherwise its totals were saved at its release.
    SaveSurfaceResponseSolverMetadata(space_op.GetComm(), "Driven",
                                      corrected_ksp->NumTotalMult(),
                                      corrected_ksp->NumTotalMultIterations());
  }
  if (response_correction)
  {
    SaveMetadata(*response_correction);
  }
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op, &response_geometry,
                                            &surface_post_geometry);
  SurfaceResponseOperator *response_correction = post_op.GetSurfaceResponseOperator();
  if (response_correction && iodata.model.refinement.max_it == 0)
  {
    // Round 2b B3 (decision 1484): no AMR construction follows (Refinement.MaxIts 0).
    ReleaseResponseMatrixCaches();
  }
  const bool self_consistent_response =
      response_correction && !response_correction->IsRawFallback() &&
      iodata.solver.surface_response_correction->IncludesSelfConsistent();

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  std::size_t convergence_memory = iodata.solver.driven.adaptive_memory;
  std::size_t max_size_per_excitation = iodata.solver.driven.adaptive_max_size;
  std::size_t nprom_indices = iodata.solver.driven.prom_indices.size();
  MFEM_VERIFY(max_size_per_excitation <= 0 || max_size_per_excitation >= nprom_indices,
              "Adaptive frequency sweep must sample at least " << nprom_indices
                                                               << " frequency points!");

  // Initialize structures for storing and reducing the results of error estimation.
  const FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector E(Curl.Width()), E_corrected, B(Curl.Height()), D, D_corrected;
  E.UseDevice(true);
  B.UseDevice(true);
  if (post_op.NeedsRecoveredElectricFlux())
  {
    // The recovered flux is allocated only when interface postprocessing uses it.
    D.SetSize(estimator.GetElectricFluxSpace().GetTrueVSize());
    D.UseDevice(true);
  }
  E = 0.0;
  B = 0.0;
  if (self_consistent_response)
  {
    E_corrected.SetSize(Curl.Width());
    E_corrected.UseDevice(true);
    E_corrected = 0.0;
    if (post_op.NeedsRecoveredElectricFlux())
    {
      D_corrected.SetSize(estimator.GetElectricFluxSpace().GetTrueVSize());
      D_corrected.UseDevice(true);
    }
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
  std::unique_ptr<RomOperator> corrected_prom_op;
  if (self_consistent_response)
  {
    // E4 (USER 1221 D8): the adaptive PROM carries dK as a projected, omega-independent
    // reduced matrix beside the projected response mass.
    std::unique_ptr<ComplexOperator> prom_stiffness;
    if (response_correction->HasMagneticModels())
    {
      prom_stiffness = std::make_unique<ComplexWrapperOperator>(
          std::make_unique<SurfaceResponseStiffnessOperator>(*response_correction),
          std::unique_ptr<Operator>());
    }
    corrected_prom_op = std::make_unique<RomOperator>(
        iodata, space_op, max_size_per_excitation,
        std::make_unique<ComplexWrapperOperator>(response_correction, nullptr),
        std::move(prom_stiffness));
  }
  space_op.GetWavePortOp().SetSuppressOutput(true);
  space_op.GetWavePortOp().ConfigureReducedModelTraining(
      max_size_per_excitation, port_excitations.Size(),
      iodata.solver.driven.adaptive_circuit_synthesis);

  // Add ports to PROM if we do synthesis.
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.AddLumpedPortModesForSynthesis();
    if (corrected_prom_op)
    {
      corrected_prom_op->AddLumpedPortModesForSynthesis();
    }
    if (space_op.GetWavePortOp().Size() > 0)
    {
      // Use the band center as the reference frequency for seeding wave-port modes.
      // The choice rescales the basis vector but does not change correctness.
      const double omega_ref = 0.5 * (omega_sample.front() + omega_sample.back());
      prom_op.AddWavePortModesForSynthesis(omega_ref);
      if (corrected_prom_op)
      {
        corrected_prom_op->AddWavePortModesForSynthesis(omega_ref);
      }
    }
  }

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator. Returns the relative error of the PROM solution at the sample
  // before the update, from the coefficients of the HDM solution in the updated basis.
  auto UpdatePROM = [&](RomOperator &op, ComplexVector &field, int excitation_idx,
                        double omega, std::size_t sample_idx, bool update_indicator)
  {
    // Add the HDM solution to the PROM reduced basis.
    op.SolvePROM(omega);
    Eigen::VectorXcd d = -op.GetReducedSolution(excitation_idx);
    const auto c =
        op.UpdatePROM(field, fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx));
    op.UpdateMRI(excitation_idx, omega, c);
    d.conservativeResizeLike(Eigen::VectorXcd::Zero(c.size()));
    // The error decides whether to sample again, so it must be the same on every process,
    // while the replicated dense computations can differ in round-off between them.
    double error = op.BasisNorm(c + d) / op.BasisNorm(c);
    Mpi::GlobalMax(1, &error, space_op.GetComm());

    if (!update_indicator)
    {
      return error;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for energy postprocessing and error estimation.
    BlockTimer bt0(Timer::POSTPRO);
    Curl.Mult(field.Real(), B.Real());
    Curl.Mult(field.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k₀, so scale = 1.
      floquet_corr->AddMult(
          field, B,
          space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }

    // Measure domain energies for the error indicator only. Don't exchange face_nbr_data,
    // unless printing paraview fields.
    auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(field, B);
    estimator.AddErrorIndicator(field, B, total_domain_energy, indicator);
    return error;
  };

  auto BuildPROM = [&](RomOperator &op, ComplexVector &field, std::string_view description,
                       bool update_indicator)
  {
    // Restart should not really be used for adaptive sweeps, but must work. Construct each
    // PROM in the same way regardless of restart for consistency.
    int excitation_counter = 0;
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      if (port_excitations.Size() > 1)
      {
        Mpi::Print("\nAdding {}excitation index {:d} ({:d}/{:d}):\n", description,
                   excitation_idx, ++excitation_counter, port_excitations.Size());
      }
      // Initialize PROM with explicit HDM samples, record the estimate but do not act on
      // it.
      std::vector<double> max_errors;
      std::size_t counter_rom_sample = 0;
      for (auto i : iodata.solver.driven.prom_indices)
      {
        auto omega = omega_sample[i];
        op.SolveHDM(excitation_idx, omega, field);
        max_errors.push_back(UpdatePROM(op, field, excitation_idx, omega,
                                        counter_rom_sample++, update_indicator));
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
      for (std::size_t it0 = it;
           it < max_size_per_excitation && memory < convergence_memory; it++)
      {
        // Compute the location of the maximum error in parameter domain (bounded by the
        // previous samples).
        double omega_star = op.FindMaxError(excitation_idx);

        // Sample HDM and add solution to basis.
        op.SolveHDM(excitation_idx, omega_star, field);
        const auto n = op.GetReducedDimension();
        max_errors.push_back(UpdatePROM(op, field, excitation_idx, omega_star,
                                        counter_rom_sample++, update_indicator));
        memory = max_errors.back() < offline_tol ? memory + 1 : 0;

        Mpi::Print(
            "\n{}Greedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
            "{:.3e}, memory = {:d}/{:d}\n",
            description, it - it0 + 1, n, omega_star * unit_GHz, omega_star,
            max_errors.back(), memory, convergence_memory);
      }
      Mpi::Print("\nAdaptive {}sampling{} {:d} frequency samples:\n"
                 " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
                 description,
                 (it == max_size_per_excitation) ? " reached maximum" : " converged with",
                 it, op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
                 convergence_memory);
      utils::PrettyPrint(op.GetSamplePoints(excitation_idx), unit_GHz,
                         fmt::format(" {}sampled frequencies (GHz):", description));
      utils::PrettyPrint(max_errors, 1.0, fmt::format(" {}sample errors:", description));
    }
  };

  BuildPROM(prom_op, E, "", true);
  if (corrected_prom_op)
  {
    Mpi::Print("\nConstructing independent fabrication-response corrected PROM\n");
    BuildPROM(*corrected_prom_op, E_corrected, "fabrication-response corrected ", false);
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
  post_op.ConfigureReducedPostprocessing(prom_op, iodata.solver.driven.adaptive_tol);

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  auto postprocess_online_point = [&](int excitation_idx, std::size_t omega_i)
  {
    const auto omega = omega_sample[omega_i];
    Mpi::Print("\nIt {:d}/{:d}, excitation {:d}: ω/2π = {:.3e} GHz "
               "(total elapsed time = {:.2e} s)\n",
               omega_i + 1, omega_sample.size(), excitation_idx,
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                   (2 * std::numbers::pi),
               Timer::Duration(Timer::Now() - t0).count());

    // Reduced postprocessing evaluates the outputs from the reduced solution, mostly
    // without expanding it. Never with a surface response correction: its gate and
    // corrected energies need the full field postprocessing below.
    const Eigen::VectorXcd y = prom_op.GetReducedSolution(excitation_idx);
    const bool reduced = !response_correction && post_op.HasReducedPostprocessing() &&
                         !post_op.WillWriteFields(static_cast<int>(omega_i));
    const bool expand = !reduced || post_op.ReducedPostprocessingNeedsField(excitation_idx);
    if (expand)
    {
      prom_op.ExpandPROMSolution(excitation_idx, E);
    }
    if (response_correction && !response_correction->IsRawFallback() &&
        response_correction->HasConductorConsistencyProbes())
    {
      // The conductor-consistency gate (idfix DESIGN B.3; decision 1064 (1)) on the
      // reduced-order field of this frequency, before its energies. The corrected PROM
      // was built on the pre-gate operator, so its field is marked SelfConsistentPreGate
      // when an exclusion fires here (the fixed-trace letters are exact).
      BlockTimer gate_timer(Timer::POSTPRO_RESPONSE_COUPON);
      response_correction->ApplyConductorConsistencyGate(
          E, excitation_idx,
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI),
          corrected_prom_op != nullptr);
      if (response_correction->HasMagneticModels())
      {
        response_correction->ApplyMagneticConductorConsistencyGate(
            E, excitation_idx,
            iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI));
      }
    }
    if (corrected_prom_op)
    {
      corrected_prom_op->ExpandPROMSolution(excitation_idx, E_corrected);
    }
    Mpi::Print("\n");
    if (omega_i == 0 && post_op.WillWriteFields())
    {
      // Switch ParaView subfolders once per excitation. Delay visualization setup until
      // the first online solution is ready.
      post_op.InitializeParaviewDataCollection(excitation_idx);
    }

    if (reduced)
    {
      post_op.MeasureAndPrintReduced(excitation_idx, int(omega_i), expand ? &E : nullptr,
                                     omega, y, prom_op);
      return;
    }

    // Start full post-processing.
    BlockTimer bt0(Timer::POSTPRO);
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));
    if (corrected_prom_op)
    {
      Mpi::Print(" Corrected sol. ||E|| = {:.6e}\n",
                 linalg::Norml2(space_op.GetComm(), E_corrected));
    }

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
    if (post_op.NeedsRecoveredElectricFlux())
    {
      Mpi::Print(" Recovering electric flux for interface postprocessing\n");
      estimator.RecoverElectricFlux(E, D);
      post_op.SetRecoveredElectricFlux(D);
    }
    if (corrected_prom_op)
    {
      const ComplexVector *corrected_flux = nullptr;
      if (post_op.NeedsRecoveredElectricFlux())
      {
        Mpi::Print(" Recovering electric flux for corrected interface postprocessing\n");
        estimator.RecoverElectricFlux(E_corrected, D_corrected);
        corrected_flux = &D_corrected;
      }
      post_op.SetSurfaceResponseCorrectedField(E_corrected, corrected_flux);
    }
    post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);
  };

  // The PROM (and the port modes) depend on frequency, not excitation: assemble, factor,
  // and solve it for all excitations once per frequency.
  if (!post_op.WillWriteFields() && !response_correction)
  {
    // Not with a surface response correction: its conductor-consistency gate exclusions
    // are sticky (a later solve sees every earlier exclusion), so its sweep stays
    // excitation-major as below.
    for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
    {
      prom_op.SolvePROM(omega_sample[omega_i]);
      for (const auto &[excitation_idx, excitation_spec] : port_excitations)
      {
        postprocess_online_point(excitation_idx, omega_i);
      }
    }
  }
  else
  {
    // Excitation-major ordering for the excitation-specific field collections and for a
    // surface response correction.
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      if (port_excitations.Size() > 1)
      {
        Mpi::Print("\nSweeping excitation index {:d}:\n", excitation_idx);
      }
      for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
      {
        prom_op.SolvePROM(omega_sample[omega_i]);
        if (corrected_prom_op)
        {
          corrected_prom_op->SolvePROM(omega_sample[omega_i]);
        }
        postprocess_online_point(excitation_idx, omega_i);
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
  if (self_consistent_response)
  {
    MFEM_ASSERT(corrected_prom_op, "Missing corrected adaptive driven solver!");
    const auto &corrected_solver = corrected_prom_op->GetLinearSolver();
    SaveSurfaceResponseSolverMetadata(space_op.GetComm(), "AdaptiveDrivenHDM",
                                      corrected_solver.NumTotalMult(),
                                      corrected_solver.NumTotalMultIterations());
  }
  if (response_correction)
  {
    SaveMetadata(*response_correction);
  }
  return indicator;
}

}  // namespace palace
