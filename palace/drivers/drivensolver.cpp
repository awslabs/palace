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
  int n_step = GetNumSteps(iodata.solver.driven.min_f, iodata.solver.driven.max_f,
                           iodata.solver.driven.delta_f);
  int step0 = (iodata.solver.driven.rst > 0) ? iodata.solver.driven.rst - 1 : 0;
  double delta_omega = iodata.solver.driven.delta_f;
  double omega0 = iodata.solver.driven.min_f + step0 * delta_omega;
  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  if (adaptive && n_step <= 2)
  {
    Mpi::Warning("Adaptive frequency sweep requires > 2 total frequency samples!\n"
                 "Reverting to uniform sweep!\n");
    adaptive = false;
  }
  SaveMetadata(space_op.GetNDSpaces());

  // Frequencies will be sampled uniformly in the frequency domain. Index sets are for
  // computing things like S-parameters in postprocessing.
  PostOperator post_op(iodata, space_op, "driven");
  {
    Mpi::Print("\nComputing {}frequency response for:\n", adaptive ? "adaptive fast " : "");
    bool first = true;
    for (const auto &[idx, data] : space_op.GetLumpedPortOp())
    {
      if (data.excitation)
      {
        if (first)
        {
          Mpi::Print(" Lumped port excitation specified on port{}",
                     (space_op.GetLumpedPortOp().Size() > 1) ? "s" : "");
          first = false;
        }
        Mpi::Print(" {:d}", idx);
      }
    }
    int excitations = first;
    first = true;
    for (const auto &[idx, data] : space_op.GetWavePortOp())
    {
      if (data.excitation)
      {
        if (first)
        {
          Mpi::Print(" Wave port excitation specified on port{}",
                     (space_op.GetWavePortOp().Size() > 1) ? "s" : "");
          first = false;
        }
        Mpi::Print(" {:d}", idx);
      }
    }
    excitations += first;
    first = true;
    for (const auto &[idx, data] : space_op.GetSurfaceCurrentOp())
    {
      if (first)
      {
        Mpi::Print(" Surface current excitation specified on port{}",
                   (space_op.GetSurfaceCurrentOp().Size() > 1) ? "s" : "");
        first = false;
      }
      Mpi::Print(" {:d}", idx);
    }
    excitations += first;
    MFEM_VERIFY(excitations > 0, "No excitation specified for driven simulation!");
  }
  Mpi::Print("\n");

  // Main frequency sweep loop.
  return {adaptive ? SweepAdaptive(space_op, post_op, n_step, step0, omega0, delta_omega)
                   : SweepUniform(space_op, post_op, n_step, step0, omega0, delta_omega),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op, PostOperator &post_op,
                                          int n_step, int step0, double omega0,
                                          double delta_omega) const
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver::SetOperators). Compute everything at the first frequency step.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega0, Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver and set operators for the first frequency step. The
  // preconditioner for the complex linear system is constructed from a real approximation
  // to the complex system matrix.
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega0,
                                    std::complex<double>(-omega0 * omega0, 0.0), K.get(),
                                    C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, omega0, -omega0 * omega0,
                                                             omega0);
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ksp.SetOperators(*A, *P);

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
  int step = step0;
  double omega = omega0;
  auto t0 = Timer::Now();
  while (step < n_step)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               n_step, freq, Timer::Duration(Timer::Now() - t0).count());

    // Assemble and solve the linear system.
    if (step > step0)
    {
      // Update frequency-dependent excitation and operators.
      A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
      A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                   std::complex<double>(-omega * omega, 0.0), K.get(),
                                   C.get(), M.get(), A2.get());
      P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega,
                                                            omega);
      ksp.SetOperators(*A, *P);
    }
    space_op.GetExcitationVector(omega, RHS);
    Mpi::Print("\n");
    ksp.Mult(RHS, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
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
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.UpdatePorts(space_op.GetLumpedPortOp(), space_op.GetWavePortOp(), omega);
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(space_op.GetComm(), E),
               linalg::Norml2(space_op.GetComm(), RHS));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
                 E_mag * J, (E_elec + E_mag) * J);
    }

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);

    // Postprocess S-parameters and optionally write solution to disk.
    Postprocess(post_op, space_op.GetLumpedPortOp(), space_op.GetWavePortOp(),
                space_op.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                (step == n_step - 1) ? &indicator : nullptr);

    // Increment frequency.
    step++;
    omega += delta_omega;
  }
  BlockTimer bt0(Timer::POSTPRO);
  SaveMetadata(ksp);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op, PostOperator &post_op,
                                           int n_step, int step0, double omega0,
                                           double delta_omega) const
{
  // Configure default parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int max_size = iodata.solver.driven.adaptive_max_size;
  MFEM_VERIFY(max_size <= 0 || max_size > 2,
              "Adaptive frequency sweep must sample at least two frequency points!");
  if (max_size <= 0)
  {
    max_size = 20;  // Default value
  }
  max_size = std::min(max_size, n_step - step0);  // Maximum size dictated by sweep
  int convergence_memory = iodata.solver.driven.adaptive_memory;

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
  const double f0 = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, 1.0);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             n_step - step0, omega0 * f0,
             (omega0 + (n_step - step0 - 1) * delta_omega) * f0);
  RomOperator prom_op(iodata, space_op, max_size);
  space_op.GetWavePortOp().SetSuppressOutput(
      true);  // Suppress wave port output for offline

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](double omega)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(omega, E);

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
    post_op.SetEGridFunction(E, false);
    post_op.SetBGridFunction(B, false);
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);
  };
  prom_op.SolveHDM(omega0, E);
  UpdatePROM(omega0);
  prom_op.SolveHDM(omega0 + (n_step - step0 - 1) * delta_omega, E);
  UpdatePROM(omega0 + (n_step - step0 - 1) * delta_omega);

  // Greedy procedure for basis construction (offline phase). Basis is initialized with
  // solutions at frequency sweep endpoints.
  int it = 2, it0 = it, memory = 0;
  std::vector<double> max_errors = {0.0, 0.0};
  while (true)
  {
    // Compute the location of the maximum error in parameter domain (bounded by the
    // previous samples).
    double omega_star = prom_op.FindMaxError()[0];

    // Compute the actual solution error at the given parameter point.
    prom_op.SolveHDM(omega_star, E);
    prom_op.SolvePROM(omega_star, Eh);
    linalg::AXPY(-1.0, E, Eh);
    max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                         linalg::Norml2(space_op.GetComm(), E));
    if (max_errors.back() < offline_tol)
    {
      if (++memory == convergence_memory)
      {
        break;
      }
    }
    else
    {
      memory = 0;
    }
    if (it == max_size)
    {
      break;
    }

    // Sample HDM and add solution to basis.
    Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
               "{:.3e}{}\n",
               it - it0 + 1, prom_op.GetReducedDimension(), omega_star * f0, omega_star,
               max_errors.back(),
               (memory == 0)
                   ? ""
                   : fmt::format(", memory = {:d}/{:d}", memory, convergence_memory));
    UpdatePROM(omega_star);
    it++;
  }
  Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
             " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
             (it == max_size) ? " reached maximum" : " converged with", it,
             prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
             convergence_memory);
  utils::PrettyPrint(prom_op.GetSamplePoints(), f0, " Sampled frequencies (GHz):");
  utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root

  // XX TODO: Add output of eigenvalue estimates from the PROM system (and nonlinear EVP in
  //          the general case with wave ports, etc.?)

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  int step = step0;
  double omega = omega0;
  while (step < n_step)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               n_step, freq, Timer::Duration(Timer::Now() - t0).count());

    // Assemble and solve the PROM linear system.
    prom_op.SolvePROM(omega, E);
    Mpi::Print("\n");

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
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
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.UpdatePorts(space_op.GetLumpedPortOp(), space_op.GetWavePortOp(), omega);
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
                 E_mag * J, (E_elec + E_mag) * J);
    }

    // Postprocess S-parameters and optionally write solution to disk.
    Postprocess(post_op, space_op.GetLumpedPortOp(), space_op.GetWavePortOp(),
                space_op.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                (step == n_step - 1) ? &indicator : nullptr);

    // Increment frequency.
    step++;
    omega += delta_omega;
  }
  BlockTimer bt0(Timer::POSTPRO);
  SaveMetadata(prom_op.GetLinearSolver());
  return indicator;
}

int DrivenSolver::GetNumSteps(double start, double end, double delta) const
{
  if (end < start)
  {
    return 1;
  }
  MFEM_VERIFY(delta != 0.0, "Zero frequency step is not allowed!");
  constexpr double delta_eps = 1.0e-9;  // 9 digits of precision comparing endpoint
  double dnfreq = std::abs(end - start) / std::abs(delta);
  int n_step = 1 + static_cast<int>(dnfreq);
  double dfinal = start + n_step * delta;
  return n_step + ((delta < 0.0 && dfinal - end > -delta_eps * end) ||
                   (delta > 0.0 && dfinal - end < delta_eps * end));
}

void DrivenSolver::Postprocess(const PostOperator &post_op,
                               const LumpedPortOperator &lumped_port_op,
                               const WavePortOperator &wave_port_op,
                               const SurfaceCurrentOperator &surf_j_op, int step,
                               double omega, double E_elec, double E_mag,
                               const ErrorIndicator *indicator) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main frequency sweep loop.
  const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
  const double E_cap = post_op.GetLumpedCapacitorEnergy(lumped_port_op);
  const double E_ind = post_op.GetLumpedInductorEnergy(lumped_port_op);
  PostprocessCurrents(post_op, surf_j_op, step, omega);
  PostprocessPorts(post_op, lumped_port_op, step, omega);
  if (surf_j_op.Size() == 0)
  {
    PostprocessSParameters(post_op, lumped_port_op, wave_port_op, step, omega);
  }
  PostprocessDomains(post_op, "f (GHz)", step, freq, E_elec, E_mag, E_cap, E_ind);
  PostprocessSurfaces(post_op, "f (GHz)", step, freq, E_elec + E_cap, E_mag + E_ind);
  PostprocessProbes(post_op, "f (GHz)", step, freq);
  if (iodata.solver.driven.delta_post > 0 && step % iodata.solver.driven.delta_post == 0)
  {
    Mpi::Print("\n");
    PostprocessFields(post_op, step / iodata.solver.driven.delta_post, freq);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(post_op, *indicator, iodata.solver.driven.delta_post > 0);
  }
}

namespace
{

struct CurrentData
{
  const int idx;       // Current source index
  const double I_inc;  // Excitation current
};

struct PortVIData
{
  const int idx;                        // Lumped port index
  const bool excitation;                // Flag for excited ports
  const double V_inc, I_inc;            // Incident voltage, current
  const std::complex<double> V_i, I_i;  // Port voltage, current
};

struct PortSData
{
  const int idx;                    // Port index
  const std::complex<double> S_ij;  // Scattering parameter
};

}  // namespace

void DrivenSolver::PostprocessCurrents(const PostOperator &post_op,
                                       const SurfaceCurrentOperator &surf_j_op, int step,
                                       double omega) const
{
  // Postprocess the frequency domain surface current excitations.
  if (post_dir.length() == 0)
  {
    return;
  }
  std::vector<CurrentData> j_data;
  j_data.reserve(surf_j_op.Size());
  for (const auto &[idx, data] : surf_j_op)
  {
    const double I_inc = data.GetExcitationCurrent();
    j_data.push_back({idx, iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_inc)});
  }
  if (root && !j_data.empty())
  {
    std::string path = post_dir + "surface-I.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", "f (GHz)", table.w1);
      for (const auto &data : j_data)
      {
        // clang-format off
        output.print("{:>{}s}{}",
                     "I_inc[" + std::to_string(data.idx) + "] (A)", table.w,
                     (data.idx == j_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    // clang-format off
    output.print("{:{}.{}e},",
                 iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega),
                 table.w1, table.p1);
    // clang-format on
    for (const auto &data : j_data)
    {
      // clang-format off
      output.print("{:+{}.{}e}{}",
                   data.I_inc, table.w, table.p,
                   (data.idx == j_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void DrivenSolver::PostprocessPorts(const PostOperator &post_op,
                                    const LumpedPortOperator &lumped_port_op, int step,
                                    double omega) const
{
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  if (post_dir.length() == 0)
  {
    return;
  }
  std::vector<PortVIData> port_data;
  port_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    const double V_inc = data.GetExcitationVoltage();
    const double I_inc = (std::abs(V_inc) > 0.0) ? data.GetExcitationPower() / V_inc : 0.0;
    const std::complex<double> V_i = post_op.GetPortVoltage(lumped_port_op, idx);
    const std::complex<double> I_i = post_op.GetPortCurrent(lumped_port_op, idx);
    port_data.push_back({idx, data.excitation,
                         iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, V_inc),
                         iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_inc),
                         iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, V_i),
                         iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_i)});
  }
  if (root && !port_data.empty())
  {
    // Write the port voltages.
    {
      std::string path = post_dir + "port-V.csv";
      auto output = OutputFile(path, (step > 0));
      if (step == 0)
      {
        output.print("{:>{}s},", "f (GHz)", table.w1);
        for (const auto &data : port_data)
        {
          if (data.excitation)
          {
            // clang-format off
            output.print("{:>{}s},",
                         "V_inc[" + std::to_string(data.idx) + "] (V)", table.w);
            // clang-format on
          }
        }
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{V[" + std::to_string(data.idx) + "]} (V)", table.w,
                       "Im{V[" + std::to_string(data.idx) + "]} (V)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega),
                   table.w1, table.p1);
      // clang-format on
      for (const auto &data : port_data)
      {
        if (data.excitation)
        {
          // clang-format off
          output.print("{:+{}.{}e},",
                       data.V_inc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.V_i.real(), table.w, table.p,
                     data.V_i.imag(), table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }

    // Write the port currents.
    {
      std::string path = post_dir + "port-I.csv";
      auto output = OutputFile(path, (step > 0));
      if (step == 0)
      {
        output.print("{:>{}s},", "f (GHz)", table.w1);
        for (const auto &data : port_data)
        {
          if (data.excitation)
          {
            // clang-format off
            output.print("{:>{}s},",
                         "I_inc[" + std::to_string(data.idx) + "] (A)", table.w);
            // clang-format on
          }
        }
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{I[" + std::to_string(data.idx) + "]} (A)", table.w,
                       "Im{I[" + std::to_string(data.idx) + "]} (A)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega),
                   table.w1, table.p1);
      // clang-format on
      for (const auto &data : port_data)
      {
        if (data.excitation)
        {
          // clang-format off
          output.print("{:+{}.{}e},",
                       data.I_inc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.I_i.real(), table.w, table.p,
                     data.I_i.imag(), table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
  }
}

void DrivenSolver::PostprocessSParameters(const PostOperator &post_op,
                                          const LumpedPortOperator &lumped_port_op,
                                          const WavePortOperator &wave_port_op, int step,
                                          double omega) const
{
  // Postprocess S-parameters. This computes a column of the S matrix corresponding to the
  // excited port index specified in the configuration file, storing |S_ij| and arg
  // (S_ij) in dB and degrees, respectively. S-parameter output is only available for a
  // single lumped or wave port excitation.
  bool src_lumped_port = false;
  bool src_wave_port = false;
  int source_idx = -1;
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      if (src_lumped_port || src_wave_port)
      {
        return;
      }
      src_lumped_port = true;
      source_idx = idx;
    }
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    if (data.excitation)
    {
      if (src_lumped_port || src_wave_port)
      {
        return;
      }
      src_wave_port = true;
      source_idx = idx;
    }
  }
  if (!src_lumped_port && !src_wave_port)
  {
    return;
  }
  std::vector<PortSData> port_data;
  port_data.reserve(src_lumped_port ? lumped_port_op.Size() : wave_port_op.Size());
  if (src_lumped_port)
  {
    // Compute lumped port S-parameters.
    for (const auto &[idx, data] : lumped_port_op)
    {
      const std::complex<double> S_ij =
          post_op.GetSParameter(lumped_port_op, idx, source_idx);
      port_data.push_back({idx, S_ij});
    }
  }
  else  // src_wave_port
  {
    // Compute wave port S-parameters.
    for (const auto &[idx, data] : wave_port_op)
    {
      const std::complex<double> S_ij =
          post_op.GetSParameter(wave_port_op, idx, source_idx);
      port_data.push_back({idx, S_ij});
    }
  }

  // Print table to stdout.
  for (const auto &data : port_data)
  {
    std::string str =
        "S[" + std::to_string(data.idx) + "][" + std::to_string(source_idx) + "]";
    // clang-format off
    Mpi::Print(" {} = {:+.3e}{:+.3e}i, |{}| = {:+.3e}, arg({}) = {:+.3e}\n",
               str, data.S_ij.real(), data.S_ij.imag(),
               str, 20.0 * std::log10(std::abs(data.S_ij)),
               str, std::arg(data.S_ij) * 180.0 / M_PI);
    // clang-format on
  }

  // Print table to file.
  if (root && post_dir.length() > 0)
  {
    std::string path = post_dir + "port-S.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", "f (GHz)", table.w1);
      for (const auto &data : port_data)
      {
        std::string str =
            "S[" + std::to_string(data.idx) + "][" + std::to_string(source_idx) + "]";
        // clang-format off
        output.print("{:>{}s},{:>{}s}{}",
                     "|" + str + "| (dB)", table.w,
                     "arg(" + str + ") (deg.)", table.w,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    // clang-format off
    output.print("{:{}.{}e},",
                 iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega),
                 table.w1, table.p1);
    for (const auto &data : port_data)
    {
      // clang-format off
      output.print("{:>+{}.{}e},{:>+{}.{}e}{}",
                   20.0 * std::log10(std::abs(data.S_ij)), table.w, table.p,
                   std::arg(data.S_ij) * 180.0 / M_PI, table.w, table.p,
                   (data.idx == port_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

}  // namespace palace
