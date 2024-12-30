// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
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
#include "utils/tablecsv.hpp"
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
  PostprocessPrintResults post_results(root, post_dir, post_op, space_op, n_step,
                                       iodata.solver.driven.delta_post);

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
  return {adaptive ? SweepAdaptive(space_op, post_op, post_results, n_step, step0, omega0,
                                   delta_omega)
                   : SweepUniform(space_op, post_op, post_results, n_step, step0, omega0,
                                  delta_omega),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op, PostOperator &post_op,
                                          PostprocessPrintResults &post_results, int n_step,
                                          int step0, double omega0,
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

  // Main frequency sweep loop.
  double omega = omega0;
  auto t0 = Timer::Now();
  for (int step = step0; step < n_step; step++, omega += delta_omega)
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
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.SetFrequency(omega);
    post_op.MeasureAll();

    Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(space_op.GetComm(), E),
               linalg::Norml2(space_op.GetComm(), RHS));

    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();

    const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
    Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
               E_mag * J, (E_elec + E_mag) * J);

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);

    post_results.PostprocessStep(iodata, post_op, space_op, step);
  }
  // Final postprocessing & printing
  BlockTimer bt0(Timer::POSTPRO);
  SaveMetadata(ksp);
  post_results.PostprocessFinal(post_op, indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op, PostOperator &post_op,
                                           PostprocessPrintResults &post_results,
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
  double omega = omega0;
  for (int step = step0; step < n_step; step++, omega += delta_omega)
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
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.SetFrequency(omega);
    post_op.MeasureAll();

    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));

    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
    Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
               E_mag * J, (E_elec + E_mag) * J);

    post_results.PostprocessStep(iodata, post_op, space_op, step);
  }
  // Final postprocessing & printing
  BlockTimer bt0(Timer::POSTPRO);
  SaveMetadata(prom_op.GetLinearSolver());
  post_results.PostprocessFinal(post_op, indicator);
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

// -----------------
// Measurements / Postprocessing

DrivenSolver::CurrentsPostPrinter::CurrentsPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir,
    const SurfaceCurrentOperator &surf_j_op, int n_expected_rows)
  : root_{root}, do_measurement_{
                     do_measurement             //
                     && (surf_j_op.Size() > 0)  // Needs surface currents
                 }
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  surface_I = TableWithCSVFile(post_dir / "surface-I.csv");
  surface_I.table.reserve(n_expected_rows, surf_j_op.Size());
  surface_I.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));
  for (const auto &[idx, data] : surf_j_op)
  {
    surface_I.table.insert_column(fmt::format("I_{}", idx),
                                  fmt::format("I_inc[{}] (A)", idx));
  }
  surface_I.AppendHeader();
}

void DrivenSolver::CurrentsPostPrinter::AddMeasurement(
    double freq, const SurfaceCurrentOperator &surf_j_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  surface_I.table["idx"] << freq;
  for (const auto &[idx, data] : surf_j_op)
  {
    auto I_inc = data.GetExcitationCurrent();
    surface_I.table[format("I_{}", idx)] << iodata.DimensionalizeValue(VT::CURRENT, I_inc);
  }
  surface_I.AppendRow();
}

DrivenSolver::PortsPostPrinter::PortsPostPrinter(bool do_measurement, bool root,
                                                 const fs::path &post_dir,
                                                 const LumpedPortOperator &lumped_port_op,
                                                 int n_expected_rows)
  : root_{root}, do_measurement_{
                     do_measurement                  //
                     && (lumped_port_op.Size() > 0)  // Only works for lumped ports
                 }
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;
  port_V = TableWithCSVFile(post_dir / "port-V.csv");
  port_V.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_V.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  port_I = TableWithCSVFile(post_dir / "port-I.csv");
  port_I.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_I.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      port_V.table.insert_column(format("inc{}", idx), format("V_inc[{}] (V)", idx));
      port_I.table.insert_column(format("inc{}", idx), format("I_inc[{}] (A)", idx));
    }

    port_V.table.insert_column(format("re{}", idx), format("Re{{V[{}]}} (V)", idx));
    port_V.table.insert_column(format("im{}", idx), format("Im{{V[{}]}} (V)", idx));

    port_I.table.insert_column(format("re{}", idx), format("Re{{I[{}]}} (A)", idx));
    port_I.table.insert_column(format("im{}", idx), format("Im{{I[{}]}} (A)", idx));
  }
  port_V.AppendHeader();
  port_I.AppendHeader();
}

void DrivenSolver::PortsPostPrinter::AddMeasurement(
    double freq, const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
    const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;

  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  port_V.table["idx"] << freq;
  port_I.table["idx"] << freq;

  auto unit_V = iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
  auto unit_A = iodata.DimensionalizeValue(VT::CURRENT, 1.0);

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      double V_inc = data.GetExcitationVoltage();
      double I_inc = (std::abs(V_inc) > 0.0) ? data.GetExcitationPower() / V_inc : 0.0;

      port_V.table[fmt::format("inc{}", idx)] << V_inc * unit_V;
      port_I.table[fmt::format("inc{}", idx)] << I_inc * unit_A;
    }

    std::complex<double> V_i = post_op.GetPortVoltage(idx);
    std::complex<double> I_i = post_op.GetPortCurrent(idx);

    port_V.table[fmt::format("re{}", idx)] << V_i.real() * unit_V;
    port_V.table[fmt::format("im{}", idx)] << V_i.imag() * unit_V;

    port_I.table[fmt::format("re{}", idx)] << I_i.real() * unit_A;
    port_I.table[fmt::format("im{}", idx)] << I_i.imag() * unit_A;
  }
  port_V.AppendRow();
  port_I.AppendRow();
}

DrivenSolver::SParametersPostPrinter::SParametersPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir,
    const LumpedPortOperator &lumped_port_op, const WavePortOperator &wave_port_op,
    int n_expected_rows)
  : root_{root},
    do_measurement_{
        do_measurement  //
        && ((lumped_port_op.Size() > 0) xor
            (wave_port_op.Size() > 0))  // either lumped or wave but not both

    },
    src_lumped_port{lumped_port_op.Size() > 0}
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  // Get excitation index as is currently done: if -1 then no excitation
  // Already ensured that one of lumped or wave ports are empty
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      source_idx = idx;
    }
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    if (data.excitation)
    {
      source_idx = idx;
    }
  }

  do_measurement_ = do_measurement_ && (source_idx > 0);

  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;
  port_S = TableWithCSVFile(post_dir / "port-S.csv");
  port_S.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_S.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  // Already ensured that one of lumped or wave ports are empty
  for (const auto &[o_idx, data] : lumped_port_op)
  {
    port_S.table.insert_column(format("abs_{}_{}", o_idx, source_idx),
                               format("|S[{}][{}]| (dB)", o_idx, source_idx));
    port_S.table.insert_column(format("arg_{}_{}", o_idx, source_idx),
                               format("arg(S[{}][{}]) (deg.)", o_idx, source_idx));
  }
  for (const auto &[o_idx, data] : wave_port_op)
  {
    port_S.table.insert_column(format("abs_{}_{}", o_idx, source_idx),
                               format("|S[{}][{}]| (dB)", o_idx, source_idx));
    port_S.table.insert_column(format("arg_{}_{}", o_idx, source_idx),
                               format("arg(S[{}][{}]) (deg.)", o_idx, source_idx));
  }
  port_S.AppendHeader();
}

void DrivenSolver::SParametersPostPrinter::AddMeasurement(
    double freq, const PostOperator &post_op, const LumpedPortOperator &lumped_port_op,
    const WavePortOperator &wave_port_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Add frequencies
  port_S.table["idx"] << freq;

  std::vector<int> all_port_indices;
  for (const auto &[idx, data] : lumped_port_op)
  {
    all_port_indices.emplace_back(idx);
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    all_port_indices.emplace_back(idx);
  }

  for (const auto o_idx : all_port_indices)
  {
    std::complex<double> S_ij = post_op.GetSParameter(src_lumped_port, o_idx, source_idx);

    auto abs_S_ij = 20.0 * std::log10(std::abs(S_ij));
    auto arg_S_ij = std::arg(S_ij) * 180.0 / M_PI;

    port_S.table[format("abs_{}_{}", o_idx, source_idx)] << abs_S_ij;
    port_S.table[format("arg_{}_{}", o_idx, source_idx)] << arg_S_ij;

    Mpi::Print(" {sij} = {:+.3e}{:+.3e}i, |{sij}| = {:+.3e}, arg({sij}) = {:+.3e}\n",
               S_ij.real(), S_ij.imag(), abs_S_ij, arg_S_ij,
               fmt::arg("sij", format("S[{}][{}]", o_idx, source_idx)));
  }
  // Regenerate from scratch each time since not row-wise (TODO: improve)
  port_S.WriteFullTableTrunc();
}

DrivenSolver::PostprocessPrintResults::PostprocessPrintResults(
    bool root, const fs::path &post_dir, const PostOperator &post_op,
    const SpaceOperator &space_op, int n_expected_rows, int delta_post_)
  : delta_post{delta_post_}, write_paraview_fields{delta_post_ > 0},
    domains{true, root, post_dir, post_op, "f (GHz)", n_expected_rows},
    surfaces{true, root, post_dir, post_op, "f (GHz)", n_expected_rows},
    currents{true, root, post_dir, space_op.GetSurfaceCurrentOp(), n_expected_rows},
    probes{true, root, post_dir, post_op, "f (GHz)", n_expected_rows},
    ports{true, root, post_dir, space_op.GetLumpedPortOp(), n_expected_rows},
    s_parameters{true,
                 root,
                 post_dir,
                 space_op.GetLumpedPortOp(),
                 space_op.GetWavePortOp(),
                 n_expected_rows},
    error_indicator{true, root, post_dir}
{
}

void DrivenSolver::PostprocessPrintResults::PostprocessStep(const IoData &iodata,
                                                            const PostOperator &post_op,
                                                            const SpaceOperator &space_op,
                                                            int step)
{
  double omega = post_op.GetFrequency().real();
  auto freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);

  domains.AddMeasurement(freq, post_op, iodata);
  surfaces.AddMeasurement(freq, post_op, iodata);
  currents.AddMeasurement(freq, space_op.GetSurfaceCurrentOp(), iodata);
  probes.AddMeasurement(freq, post_op, iodata);
  ports.AddMeasurement(freq, post_op, space_op.GetLumpedPortOp(), iodata);
  s_parameters.AddMeasurement(freq, post_op, space_op.GetLumpedPortOp(),
                              space_op.GetWavePortOp(), iodata);
  // The internal GridFunctions in PostOperator have already been set:
  if (write_paraview_fields && (step % delta_post == 0))
  {
    Mpi::Print("\n");
    post_op.WriteFields(step / delta_post, freq);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
}

void DrivenSolver::PostprocessPrintResults::PostprocessFinal(
    const PostOperator &post_op, const ErrorIndicator &indicator)
{
  BlockTimer bt0(Timer::POSTPRO);
  auto indicator_stats = indicator.GetSummaryStatistics(post_op.GetComm());
  error_indicator.PrintIndicatorStatistics(post_op, indicator_stats);
  if (write_paraview_fields)
  {
    post_op.WriteFieldsFinal(&indicator);
  }
}

}  // namespace palace
