// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <iostream>
#include <map>
#include <numbers>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
#include "drivers/drivensubstructuredsweep.hpp"
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/initialguess.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/drivensubstructure.hpp"
#include "models/environmentoperator.hpp"
#include "models/floquetportoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/regionsolver.hpp"
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

namespace
{

// The sweep of a substructured region, for its messages: online, or offline adaptive.
const char *RegionSweep(const IoData &iodata)
{
  return iodata.solver.substructuring->Online() ? "Online" : "Adaptive";
}

// The postprocessing of a sweep: that of the configuration or, for a substructured region,
// whose field is known in the region and on Γ only, that of the region (see
// RegionPostprocessing).
std::pair<config::DomainData, config::BoundaryData>
SweepPostprocessing(const IoData &iodata, const SpaceOperator &space_op)
{
  config::DomainData domains = iodata.domains;
  config::BoundaryData boundaries = iodata.boundaries;
  if (const auto *env = space_op.GetEnvironment())
  {
    RegionPostprocessing(env->RegionAttributes(), space_op.GetMesh().Get(),
                         env->RegionMesh() ? &env->RegionBoundaryAttributes() : nullptr,
                         !env->Ports().empty(), RegionSweep(iodata), domains, boundaries);
  }
  return {domains, boundaries};
}

// For a substructured region, the values of the environment's ports from the model at a
// solution E of an excitation, or its values on Γ on this process, for the next
// measurements.
template <typename T>
void SetEnvironmentPorts(SpaceOperator &space_op,
                         PostOperator<ProblemType::DRIVEN> &post_op, int excitation_idx,
                         double omega, const T &E)
{
  if (auto *env = space_op.GetEnvironment())
  {
    std::map<int, std::complex<double>> V, S;
    env->PortValues(excitation_idx, omega, E, V, S);
    post_op.SetLumpedPortVoltages(std::move(V));
    post_op.SetWavePortOverlaps(std::move(S));
  }
}

void WarnRegionFields(const IoData &iodata, const SpaceOperator &space_op,
                      const PostOperator<ProblemType::DRIVEN> &post_op)
{
  if (space_op.GetEnvironment() && post_op.WillWriteFields())
  {
    Mpi::Warning("{} driven substructuring writes the fields in the region and on the "
                 "interface only (zero in the environment)!\n",
                 RegionSweep(iodata));
  }
}

}  // namespace

std::pair<ErrorIndicator, long long int>
DrivenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Set up the spatial discretization and frequency sweep.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto full_op = std::make_unique<SpaceOperator>(iodata, mesh);
  SpaceOperator &space_op = *full_op;
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

  // Main frequency sweep loop. A substructured region is solved against its environment,
  // with the operators of the region (see SpaceOperator): the environment of the saved
  // model online, or sampled offline for a rational model by an adaptive sweep (on the
  // operators of the full problem, which are then released). An offline exact sweep solves
  // the whole device frequency by frequency.
  if (const auto &sub = iodata.solver.substructuring)
  {
    if (!sub->Online() && !adaptive)
    {
      return {SweepSubstructured(iodata, space_op), space_op.GlobalTrueVSize()};
    }
    std::unique_ptr<EnvironmentOperator> env;
    if (sub->Online())
    {
      env = std::make_unique<EnvironmentOperator>(iodata, space_op);
      if (adaptive && env->NumRecords() == 0)
      {
        Mpi::Warning("An adaptive online sweep needs a model saved by an adaptive sweep: "
                     "every frequency is solved!\n");
        adaptive = false;
      }
    }
    else
    {
      env = SampleEnvironmentModel(iodata, space_op);
    }
    full_op.reset();
    SpaceOperator region_op(iodata, mesh, env.get());
    return {adaptive ? SweepAdaptive(region_op) : SweepUniform(region_op),
            region_op.GlobalTrueVSize()};
  }
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  const auto *env = space_op.GetEnvironment();
  const auto [domains, boundaries] = SweepPostprocessing(iodata, space_op);
  PostOperator<ProblemType::DRIVEN> post_op(iodata.problem, iodata.solver, domains,
                                            boundaries, iodata.units, space_op);
  WarnRegionFields(iodata, space_op, post_op);

  // The fixed system matrices are only needed when the system cannot be applied through
  // the preconditioner (see SpaceOperator::GetSystemOperator), so assemble them on first
  // use. PEC boundaries are handled simply by setting diagonal entries of the system matrix
  // for the corresponding dofs. Because the Dirichlet BC is always homogeneous, no special
  // elimination is required on the RHS.
  std::unique_ptr<ComplexOperator> K, C, M;
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  auto ksp = MakeDrivenKspSolver(iodata, space_op);

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation (none
  // for a substructured region, whose environment field is not known).
  std::unique_ptr<FluxErrorEstimator<ComplexVector>> estimator;
  if (!env)
  {
    estimator = std::make_unique<FluxErrorEstimator<ComplexVector>>(space_op, iodata);
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

  // Main excitation and frequency loop, over the solves in their order, which "Restart"
  // counts: all frequencies of an excitation before the next one or, for a substructured
  // region, all excitations at a frequency before the next one, so that the region is
  // factored once per frequency.
  auto t0 = Timer::Now();
  std::vector<int> excitations;
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    excitations.push_back(excitation_idx);
  }
  const std::size_t n_ex = excitations.size(), n_omega = omega_sample.size();
  const bool by_frequency = (env != nullptr);
  auto excitation_of = [&](std::size_t s) { return by_frequency ? s % n_ex : s / n_omega; };
  auto frequency_of = [&](std::size_t s) { return by_frequency ? s / n_ex : s % n_omega; };
  const std::size_t first = iodata.solver.driven.restart - 1, n_solves = n_ex * n_omega;

  // With an initial guess, each solve of an excitation after the first starts from the
  // combination of its solutions at up to num_guess previous frequencies which minimizes
  // ||RHS - A E|| (for the second, the previous solution), which the iterative solver then
  // scales in the norm of its stopping test. Errors of the previous solutions which a
  // left-preconditioned stopping test does not see are carried over, amplified by up to
  // Σⱼ |cⱼ|; in tests, the outputs stayed at the level of the linear solver tolerance.
  // The images A E_prev only exist during the least squares solve. E_prev keeps fewer
  // solutions than the solves of an excitation in this run, and none if it could keep only
  // one. A substructured region is factored exactly, without initial guesses.
  constexpr std::size_t num_guess = 4;
  const std::size_t max_solves =
      (excitation_of(first) + 1 == n_ex) ? n_omega - frequency_of(first) : n_omega;
  const std::size_t capacity = std::min(num_guess, max_solves - 1);
  linalg::SolutionHistory<ComplexOperator> E_prev(
      (iodata.solver.linear.initial_guess && capacity > 1 && !by_frequency) ? capacity : 0,
      E.Size());

  // The operators of the current frequency. A refers to P and A2, so both must outlive it.
  std::unique_ptr<ComplexOperator> P, A2, A;
  for (std::size_t s = first; s < n_solves; s++)
  {
    const std::size_t k = excitation_of(s), omega_i = frequency_of(s);
    const int excitation_idx = excitations[k];
    const double omega = omega_sample[omega_i];
    const bool excitation_starts = (s == first || excitation_of(s - 1) != k);
    if (!by_frequency && excitation_starts)
    {
      if (n_ex > 1)
      {
        Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx, k + 1,
                   n_ex);
      }
      // Start each excitation from zero rather than the previous excitation's solution.
      E = 0.0;
      E_prev.Clear();
    }

    // Assemble frequency dependent matrices and initialize operators in linear solver.
    if (s == first || frequency_of(s - 1) != omega_i)
    {
      A.reset();
      A2.reset();
      P.reset();
      P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                            -omega * omega + 0.0i, omega);
      A = space_op.GetSystemOperator(omega, *P);
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
      ksp->SetOperators(*A, *P);
    }

    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s{})\n",
               omega_i + 1, n_omega,
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                   (2 * std::numbers::pi),
               Timer::Duration(Timer::Now() - t0).count(),
               (n_ex > 1) ? fmt::format(", solve {:d}/{:d}", s + 1, n_solves) : "");

    // Solve linear system.
    space_op.GetExcitationVector(excitation_idx, omega, RHS);
    if (by_frequency)
    {
      E = 0.0;
    }
    if (E_prev.Size() > 1)
    {
      BlockTimer bt(Timer::KSP_INITIAL_GUESS);
      E_prev.MinimalResidualCombination(space_op.GetComm(), *A, RHS, E);
    }

    Mpi::Print("\n");
    ksp->Mult(RHS, E);
    E_prev.Push(E);

    if (s + 1 == n_solves)
    {
      // That was the last solve: save the solver totals and release the linear solver and
      // system operators, so they do not add to the memory peak of the field output and
      // error estimation which follow.
      SaveMetadata(*ksp);
      ksp.reset();
      E_prev.Release();
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

    if (!by_frequency || n_ex == 1)
    {
      if (excitation_starts)
      {
        // Switch ParaView subfolders once per excitation, after the first solve.
        post_op.InitializeParaviewDataCollection(excitation_idx);
      }
    }
    else if (post_op.WillWriteFields(static_cast<int>(omega_i)))
    {
      // The fields of each excitation go to its own ParaView collection, continued at each
      // frequency.
      post_op.InitializeParaviewDataCollection(excitation_idx, true);
    }
    SetEnvironmentPorts(space_op, post_op, excitation_idx, omega, E);
    auto total_domain_energy =
        post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

    // Calculate and record the error indicators.
    if (estimator)
    {
      Mpi::Print(" Updating solution error estimates\n");
      estimator->AddErrorIndicator(E, B, total_domain_energy, omega, indicator);
    }

    // Final postprocessing & printing of an excitation. After the last solve, the solver
    // totals were already saved when the solver was released.
    if (ksp && !by_frequency && excitation_of(s + 1) != k)
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
  const auto *env = space_op.GetEnvironment();
  const auto [domains, boundaries] = SweepPostprocessing(iodata, space_op);
  PostOperator<ProblemType::DRIVEN> post_op(iodata.problem, iodata.solver, domains,
                                            boundaries, iodata.units, space_op);
  WarnRegionFields(iodata, space_op, post_op);

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
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation (none
  // for a substructured region, whose environment field is not known).
  std::unique_ptr<FluxErrorEstimator<ComplexVector>> estimator;
  if (!env)
  {
    estimator = std::make_unique<FluxErrorEstimator<ComplexVector>>(space_op, iodata);
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

  const auto num_port_modes = prom_op.GetReducedDimension();

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator. Returns the relative error of the PROM solution at the sample
  // before the update, from the coefficients of the HDM solution in the updated basis.
  std::size_t num_refined_samples = 0, num_unreproduced = 0;
  auto UpdatePROM = [&](int excitation_idx, double omega, std::size_t sample_idx)
  {
    // Add the HDM solution to the PROM reduced basis, refined if the PROM does not
    // reproduce it (see RomOperator::AddSamplePROM). A certified solution, started from the
    // PROM solution, can be in the span of the basis to round-off (the linear solver barely
    // changed the PROM solution, or took no iteration): such a part adds no basis vector.
    const auto update = prom_op.AddSamplePROM(
        excitation_idx, omega, E,
        fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx), offline_tol);
    if (update.num_refined > 0)
    {
      Mpi::Print(" Refined the sample at ω = {:.3e} GHz {:d} time{}: the PROM did not "
                 "reproduce it to the tolerance{}\n",
                 omega * unit_GHz, update.num_refined, (update.num_refined == 1) ? "" : "s",
                 update.reproduced ? "" : " (and still does not)");
      num_refined_samples++;
    }
    num_unreproduced += !update.reproduced;
    prom_op.UpdateMRI(excitation_idx, omega, update.c);

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
    if (estimator)
    {
      auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(E, B, omega);
      estimator->AddErrorIndicator(E, B, total_domain_energy, omega, indicator);
    }
    return update.error;
  };

  // Try to certify the PROM at a sample frequency with an HDM solve started from the PROM
  // solution: the PROM is certified if the solve converges (so E is an acceptable HDM
  // solution) with a correction below the tolerance. The solve uses the tighter of the
  // linear solver tolerance and the tolerance of the sweep, so that a PROM solution which
  // only meets a looser linear solver tolerance is not accepted without a correction. The
  // solve is stopped as soon as the correction exceeds the tolerance, since the sample then
  // has to be solved again from zero anyway. Returns the relative correction if certified,
  // and a negative value if not. Eh holds the PROM solution, only with initial guesses.
  const bool initial_guess = iodata.solver.linear.initial_guess;
  const double certification_tol =
      std::min(prom_op.GetLinearSolver().GetRelTol(), offline_tol);
  ComplexVector Eh(initial_guess ? E.Size() : 0);
  Eh.UseDevice(true);
  auto CertifyPROM = [&](int excitation_idx, double omega)
  {
    prom_op.SolvePROM(omega);
    prom_op.ExpandPROMSolution(excitation_idx, Eh);
    if (linalg::Norml2(space_op.GetComm(), Eh) == 0.0)
    {
      return -1.0;
    }
    E = Eh;
    const int num_it = prom_op.GetLinearSolver().NumTotalMultIterations();
    const bool converged = prom_op.SolveHDM(excitation_idx, omega, E, true, offline_tol);
    const int cert_it = prom_op.GetLinearSolver().NumTotalMultIterations() - num_it;
    linalg::AXPY(-1.0, E, Eh);
    const double correction =
        linalg::Norml2(space_op.GetComm(), Eh) / linalg::Norml2(space_op.GetComm(), E);
    const bool certified = converged && correction < offline_tol;
    Mpi::Print("\n Certification from the PROM solution ({:d} iteration{}, tol = {:.3e}): "
               "correction = {:.3e}{}\n",
               cert_it, (cert_it == 1) ? "" : "s", certification_tol, correction,
               certified   ? ", PROM certified"
               : converged ? ", not certified"
                           : ", not certified (solve stopped)");
    return certified ? correction : -1.0;
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
    // Initialize PROM with explicit HDM samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    std::size_t counter_rom_sample = 0;
    num_refined_samples = num_unreproduced = 0;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      auto omega = omega_sample[i];
      prom_op.SolveHDM(excitation_idx, omega, E);
      max_errors.push_back(UpdatePROM(excitation_idx, omega, counter_rom_sample++));
    }
    MFEM_VERIFY(prom_op.GetReducedDimension() > num_port_modes,
                "No HDM solution at the initial samples of the adaptive frequency sweep "
                "was added to the PROM basis: they are zero, or in the span of the port "
                "modes to the orthogonalization tolerance!");
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
    bool certified_last = false;
    for (std::size_t it0 = it; it < max_size_per_excitation && memory < convergence_memory;
         it++)
    {
      // Compute the location of the maximum error in parameter domain (bounded by the
      // previous samples).
      double omega_star = prom_op.FindMaxError(excitation_idx);

      // Sample HDM and add solution to basis. Basis vectors are solved from zero: a
      // solution started from the PROM solution keeps the error of the PROM in the
      // directions which barely affect the residual of the preconditioned system (weakly
      // excited high-Q modes, for example), and would pollute the basis. With initial
      // guesses, the sample which would complete the convergence memory is first solved
      // from the PROM solution, to certify the PROM there: if the solve converges, to the
      // tighter of the linear solver tolerance and the tolerance of the sweep, with a
      // correction below the tolerance, the PROM solution is within the tolerance of an
      // acceptable HDM solution and the sweep has converged. The certified solution is
      // still added to the basis (unless it is in its span to round-off), as the last
      // sample is without initial guesses: the PROM is then already accurate to the
      // tolerance there, so that the solution only differs from it by a correction below
      // the tolerance. Otherwise, the solve is stopped as soon as the correction exceeds
      // the tolerance, and the sample is solved again from zero (reusing the assembled
      // operators and preconditioner) and added to the basis, as without initial guesses.
      // The error of a certified sample is this correction.
      double error = (initial_guess && memory + 1 == convergence_memory)
                         ? CertifyPROM(excitation_idx, omega_star)
                         : -1.0;
      const bool certified = (error >= 0.0);
      if (!certified)
      {
        prom_op.SolveHDM(excitation_idx, omega_star, E);
      }
      const auto n = prom_op.GetReducedDimension();
      const double update_error =
          UpdatePROM(excitation_idx, omega_star, counter_rom_sample++);
      max_errors.push_back(certified ? error : update_error);
      memory = max_errors.back() < offline_tol ? memory + 1 : 0;

      Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
                 "{:.3e}, memory = {:d}/{:d}{}\n",
                 it - it0 + 1, n, omega_star * unit_GHz, omega_star, max_errors.back(),
                 memory, convergence_memory, certified ? " (certified)" : "");
      certified_last = certified;
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size_per_excitation) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    if (num_refined_samples > 0)
    {
      Mpi::Print(" Samples refined: {:d}, not reproduced after refinement: {:d}\n",
                 num_refined_samples, num_unreproduced);
    }
    if (memory < convergence_memory)
    {
      Mpi::Warning("Adaptive sampling reached \"AdaptiveMaxSamples\" = {:d} before "
                   "converging to \"AdaptiveTol\" = {:.3e} (last error estimate {:.3e}): "
                   "the outputs may not be accurate to the tolerance. Increase "
                   "\"AdaptiveMaxSamples\" or loosen \"AdaptiveTol\"!\n",
                   max_size_per_excitation, offline_tol, max_errors.back());
    }
    if (num_unreproduced > 0)
    {
      Mpi::Warning(
          "The PROM did not reproduce {:d} sample{} to \"AdaptiveTol\" = {:.3e} "
          "even after their refinement: the HDM solutions may not reach this "
          "accuracy, and the outputs near these frequencies may not be accurate to "
          "the tolerance. Loosen \"AdaptiveTol\"!\n",
          num_unreproduced, (num_unreproduced == 1) ? "" : "s", offline_tol);
    }
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
    if (certified_last)
    {
      Mpi::Print(" Last sample certified: its error is the correction of the PROM solution "
                 "by the linear solver, converged to {:.3e} (zero if no iteration was "
                 "needed)\n",
                 certification_tol);
    }
  }

  prom_op.ReleaseHDMOperators();
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
  prom_op.PrepareOnlinePhase();
  // For a substructured region, the reduced postprocessing takes the environment's port
  // values from the values on Γ of the reduced solution, those of the basis vectors
  // combined.
  post_op.ConfigureReducedPostprocessing(prom_op, iodata.solver.driven.adaptive_tol);
  std::vector<std::vector<double>> basis_gamma;
  if (env && post_op.HasReducedPostprocessing())
  {
    for (const auto &v : prom_op.GetBasis())
    {
      basis_gamma.push_back(env->LocalInterfaceValues(v));
    }
  }

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
    // without expanding it.
    const Eigen::VectorXcd y = prom_op.GetReducedSolution(excitation_idx);
    const bool reduced = post_op.HasReducedPostprocessing() &&
                         !post_op.WillWriteFields(static_cast<int>(omega_i));
    const bool expand = !reduced || post_op.ReducedPostprocessingNeedsField(excitation_idx);
    if (expand)
    {
      prom_op.ExpandPROMSolution(excitation_idx, E);
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
      if (env)
      {
        BlockTimer bt0(Timer::POSTPRO);
        std::vector<std::complex<double>> e_gamma(
            basis_gamma.empty() ? 0 : basis_gamma.front().size(), 0.0);
        for (std::size_t i = 0; i < basis_gamma.size(); i++)
        {
          for (std::size_t a = 0; a < e_gamma.size(); a++)
          {
            e_gamma[a] += y(static_cast<Eigen::Index>(i)) * basis_gamma[i][a];
          }
        }
        SetEnvironmentPorts(space_op, post_op, excitation_idx, omega, e_gamma);
      }
      post_op.MeasureAndPrintReduced(excitation_idx, int(omega_i), expand ? &E : nullptr,
                                     omega, y, prom_op);
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
    SetEnvironmentPorts(space_op, post_op, excitation_idx, omega, E);
    post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);
  };

  // The PROM (and the port modes) depend on frequency, not excitation: assemble, factor,
  // and solve it for all excitations once per frequency.
  if (!post_op.WillWriteFields())
  {
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
    // Excitation-major ordering for the excitation-specific field collections.
    for (const auto &[excitation_idx, excitation_spec] : port_excitations)
    {
      if (port_excitations.Size() > 1)
      {
        Mpi::Print("\nSweeping excitation index {:d}:\n", excitation_idx);
      }
      for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
      {
        prom_op.SolvePROM(omega_sample[omega_i]);
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
  return indicator;
}

}  // namespace palace
