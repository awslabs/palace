// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
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
  auto excitation_helper = space_op.GetPortExcitations();
  SaveMetadata(excitation_helper);

  // Frequencies will be sampled uniformly in the frequency domain, in addition to explicit
  // samples. Additionally remove any duplicates between the two.
  int n_step = config::GetNumSteps(iodata.solver.driven.min_f, iodata.solver.driven.max_f, iodata.solver.driven.delta_f);
  std::vector<double> omega(n_step);
  std::iota(omega.begin(), omega.end(), 0);
  std::for_each(omega.begin(), omega.end(), [=](double &x)
                { x = iodata.solver.driven.min_f + x * iodata.solver.driven.delta_f; });
  omega.insert(omega.end(), iodata.solver.driven.sample_f.begin(),
               iodata.solver.driven.sample_f.end());
  std::sort(omega.begin(), omega.end());
  constexpr double delta_eps = 1.0e-9; // Precision in frequency comparisons
  omega.erase(std::unique(omega.begin(), omega.end(),
                          [=](auto x, auto y) { return std::abs(x - y) < delta_eps; }),
              omega.end());
  if (iodata.solver.driven.rst > 0)
  {
    MFEM_VERIFY(iodata.solver.driven.rst <= omega.size(),
                "\"Restart\": (" << iodata.solver.driven.rst
                                 << ") is greater than the number of samples ("
                                 << omega.size() << ")!");
    omega.erase(omega.begin(), omega.begin() + iodata.solver.driven.rst - 1);
  }

  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  int nreq_samples = iodata.solver.driven.explicit_prom_sample ? 2 + iodata.solver.driven.sample_f.size() : 2;
  if (adaptive && omega.size() <= nreq_samples)
  {
    Mpi::Warning("Adaptive frequency sweep requires > {} total frequency samples!\n"
                 "Reverting to uniform sweep!\n", nreq_samples);
    adaptive = false;
  }
  SaveMetadata(space_op.GetNDSpaces());
  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             excitation_helper.FmtLog());

  // Main frequency sweep loop.
  return {adaptive ? SweepAdaptive(space_op, omega) : SweepUniform(space_op, omega),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op,
                                          const std::vector<double> &omega_sample) const
{
  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed if multiple excitations.
  PostOperator<config::ProblemData::Type::DRIVEN> post_op(iodata, space_op);
  auto excitation_helper = space_op.GetPortExcitations();

  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver::SetOperators). Compute everything at the first frequency step.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Main frequency sweep loop.
  auto t0 = Timer::Now();
  auto print_counter_excitation = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    if (excitation_helper.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++print_counter_excitation, excitation_helper.Size());
    }

    // Switch paraview subfolders: one for each excitation, if nr_excitations > 1.
    post_op.InitializeParaviewDataCollection(excitation_idx);

    // Frequency loop.
    int step0 = std::max(iodata.solver.driven.rst - 1, 0), step = step0;
    for (auto omega : omega_sample)
    {
      // Assemble frequency dependent matrices and initialize operators in linear solver.
      auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
      auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                        std::complex<double>(-omega * omega, 0.0), K.get(),
                                        C.get(), M.get(), A2.get());
      auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega,
                                                                 omega);
      ksp.SetOperators(*A, *P);

      Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
                 step + 1, omega_sample.size() + step0,
                 iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega),
                 Timer::Duration(Timer::Now() - t0).count());

      // Solve linear system.
      space_op.GetExcitationVector(excitation_idx, omega, RHS);
      Mpi::Print("\n");
      ksp.Mult(RHS, E);

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
        // Calculate B field correction for Floquet BCs.
        // B = -1/(iω) ∇ x E + 1/ω kp x E
        floquet_corr->AddMult(E, B, 1.0 / omega);
      }

      auto total_domain_energy =
          post_op.MeasureAndPrintAll(excitation_idx, step++, E, B, omega);

      // Calculate and record the error indicators.
      Mpi::Print(" Updating solution error estimates\n");
      estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
    }
    // Final postprocessing & printing.
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(ksp);
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op,
                                           const std::vector<double> &omega_sample) const
{
  // Initialize postprocessing for measurement and printers.
  PostOperator<config::ProblemData::Type::DRIVEN> post_op(iodata, space_op);
  auto excitation_helper = space_op.GetPortExcitations();

  // Configure HDM sample array first with end points then any explicit samples.
  std::vector<double> omega_hdm = {iodata.solver.driven.min_f, iodata.solver.driven.max_f};
  if (iodata.solver.driven.explicit_prom_sample)
  {
    omega_hdm.insert(omega_hdm.end(), iodata.solver.driven.sample_f.begin(),
                    iodata.solver.driven.sample_f.end());
  }

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int convergence_memory = iodata.solver.driven.adaptive_memory;
  int max_size_per_excitation = iodata.solver.driven.adaptive_max_size;
  MFEM_VERIFY(max_size_per_excitation <= 0 ||
                  max_size_per_excitation >= static_cast<int>(omega_hdm.size()),
              "Adaptive frequency sweep must sample at least " << omega_hdm.size()
                                                               << " frequency points!");
  // Maximum size — no more than nr steps needed.
  max_size_per_excitation =
      std::min(max_size_per_excitation, static_cast<int>(omega_sample.size()));

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
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
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
  const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             omega_sample.size(), omega_hdm[0] * unit_GHz, omega_hdm[1] * unit_GHz);
  RomOperator prom_op(iodata, space_op, max_size_per_excitation);
  space_op.GetWavePortOp().SetSuppressOutput(true);

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](int excitation_idx, double omega)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(E);
    prom_op.UpdateMRI(excitation_idx, omega, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for energy postprocessing and error estimation.
    BlockTimer bt0(Timer::POSTPRO);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }

    // Measure domain energies for the error indicator only. Don't exchange face_nbr_data,
    // unless printing paraview fields.
    auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(E, B);
    estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
  };

  // Loop excitations to add to PROM.
  auto print_counter_excitation_prom = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    if (excitation_helper.Size() > 1)
    {
      Mpi::Print("\nAdding excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++print_counter_excitation_prom, excitation_helper.Size());
    }
    prom_op.SetExcitationIndex(excitation_idx);  // Pre-compute RHS1

    // Initialize PROM with explicit hdm samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    for (auto omega : omega_hdm)
    {
      prom_op.SolveHDM(excitation_idx, omega, E);
      prom_op.SolvePROM(excitation_idx, omega, Eh);
      linalg::AXPY(-1.0, E, Eh);
      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      UpdatePROM(excitation_idx, omega);
    }
    // The estimates associated to the end points are assumed inaccurate.
    max_errors[0] = std::numeric_limits<double>::infinity();
    max_errors[1] = std::numeric_limits<double>::infinity();
    int memory = std::distance(max_errors.rbegin(),
                               std::find_if(max_errors.rbegin(), max_errors.rend(),
                                            [=](auto x) { return x > offline_tol; }));

    // Greedy procedure for basis construction (offline phase). Basis is initialized with
    // solutions at frequency sweep endpoints and explicit sample frequencies.
    int it = static_cast<int>(max_errors.size());
    for (int it0 = it; it < max_size_per_excitation && memory < convergence_memory; it++)
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
                 "{:.3e}{}\n",
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(),
                 (memory == 0)
                     ? ""
                     : fmt::format(", memory = {:d}/{:d}", memory, convergence_memory));
      UpdatePROM(excitation_idx, omega_star);
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

  // XX TODO: Add output of eigenvalue estimates from the PROM system (and nonlinear EVP
  // in the general case with wave ports, etc.?)

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  auto print_counter_excitation_online = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    if (excitation_helper.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++print_counter_excitation_online, excitation_helper.Size());
    }
    // Switch paraview subfolders: one for each excitation, if nr_excitations > 1.
    post_op.InitializeParaviewDataCollection(excitation_idx);

    // Frequency loop.
    int step0 = std::max(iodata.solver.driven.rst - 1, 0), step = step0;
    for (auto omega : omega_sample)
    {
      const double freq = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega);
      Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
                 omega_sample.size() + step0, freq, Timer::Duration(Timer::Now() - t0).count());

      // Assemble and solve the PROM linear system.
      prom_op.SolvePROM(excitation_idx, omega, E);
      Mpi::Print("\n");

      // Start Post-processing.
      BlockTimer bt0(Timer::POSTPRO);
      Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));

      // Compute B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E.Real(), B.Real());
      Curl.Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (space_op.GetMaterialOp().HasWaveVector())
      {
        // Calculate B field correction for Floquet BCs.
        // B = -1/(iω) ∇ x E + 1/ω kp x E
        floquet_corr->AddMult(E, B, 1.0 / omega);
      }
      post_op.MeasureAndPrintAll(excitation_idx, step++, E, B, omega);
    }
    // Final postprocessing & printing: no change to indicator since these are in PROM.
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(prom_op.GetLinearSolver());
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

}  // namespace palace
