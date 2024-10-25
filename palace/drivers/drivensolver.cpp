// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitationhelper.hpp"
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
  auto excitation_helper = space_op.BuildPortExcitationHelper();
  MFEM_VERIFY(!excitation_helper.Empty(), "No excitation specified for driven simulation!");

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
  PostprocessResults post_results{root, post_dir, space_op, excitation_helper, n_step};

  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             excitation_helper.FmtLog());

  // Main frequency sweep loop.
  return {adaptive ? SweepAdaptive(space_op, post_op, post_results, n_step, step0, omega0,
                                   delta_omega)
                   : SweepUniform(space_op, post_op, post_results, n_step, step0, omega0,
                                  delta_omega),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op, PostOperator &post_op,
                                          PostprocessResults &post_results, int n_step,
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

  auto excitation_helper = space_op.BuildPortExcitationHelper();

  // Main frequency sweep loop.
  auto t0 = Timer::Now();
  auto print_counter_excitation = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    print_counter_excitation++;
    Mpi::Print("\nSweeping Excitation Index {:d} ({:d}/{:d}):\n", excitation_idx,
               print_counter_excitation, excitation_helper.Size());

    double omega = omega0;
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
      space_op.GetExcitationVector(excitation_idx, omega, RHS);
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

      post_results.PostprocessStep(iodata, post_op, space_op, excitation_helper, step,
                                   omega, excitation_idx, E_elec, E_mag,
                                   (step == n_step - 1) ? &indicator : nullptr);

      // Postprocess S-parameters and optionally write solution to disk.
      Postprocess(post_op, space_op.GetLumpedPortOp(), space_op.GetWavePortOp(),
                  space_op.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                  (step == n_step - 1) ? &indicator : nullptr);
    }
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(ksp);
  }
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op, PostOperator &post_op,
                                           PostprocessResults &post_results, int n_step,
                                           int step0, double omega0,
                                           double delta_omega) const
{
  // Configure default parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int max_size = iodata.solver.driven.adaptive_max_size;
  MFEM_VERIFY(max_size <= 0 || max_size > 2,
              "Adaptive frequency sweep must sample at least two frequency points!");
  // Fix this max size to include ports
  if (max_size <= 0)
  {
    max_size = 50;  // Default value
  }
  int nr_ports = space_op.GetLumpedPortOp().Size();
  max_size =
      std::min(max_size, n_step - step0 + nr_ports);  // Maximum size dictated by sweep
  // need to add ports (external + virtual)

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

  auto &lumped_port_op = space_op.GetLumpedPortOp();
  for (auto &[port_idx, port_data] : lumped_port_op)
  {
    ComplexVector port_excitation_E;
    space_op.GetLumpedPortExcitationVector(port_idx, port_excitation_E, true);
    prom_op.UpdatePROM(false, 0.0, port_excitation_E);

    // Mpi::Barrier(post_op.GetComm());
    // BlockTimer bt0(Timer::POSTPRO);
    // post_op.SetEGridFunction(port_excitation_E, false);
    // post_op.SetBGridFunction(port_excitation_E, false);
    // post_op.WriteFields(100 + port_idx, 100.0 + port_idx);
    // Mpi::Barrier(post_op.GetComm());

    // // Debugging information:
    // Mpi::Print("GetExcitationVoltage {}\n", port_data.GetExcitationVoltage());

    // port_data.InitializeLinearForms(*(post_op.E)->ParFESpace());
    // auto &voltage_form = *(port_data.v);

    // Vector port_excitation_voltage_form;
    // port_excitation_voltage_form.SetSize(space_op.GetNDSpace().GetTrueVSize());
    // port_excitation_voltage_form.UseDevice(true);
    // port_excitation_voltage_form = 0.0;
    // space_op.GetNDSpace().GetProlongationMatrix()->AddMultTranspose(
    //     voltage_form, port_excitation_voltage_form);

    // // Real GF Only
    // auto GF = std::make_unique<GridFunction>(space_op.GetNDSpace(), false);
    // GF->Real().SetFromTrueDofs(port_excitation_voltage_form);
    // GF->Real().ExchangeFaceNbrData();
    // auto v_self = port_data.GetVoltage(*GF);
    // Mpi::Print("GetVoltage from port_excitation_voltage_form {} {}\n", v_self.real(),
    //            v_self.imag());
    // Mpi::Barrier(post_op.GetComm());

    // GF->Real().SetFromTrueDofs(port_excitation_E.Real());
    // GF->Real().ExchangeFaceNbrData();
    // v_self = port_data.GetVoltage(*GF);
    // Mpi::Print("GetVoltage from port_excitation_E (zero_metal true) {} {}\n",
    // v_self.real(),
    //            v_self.imag());
    // Mpi::Barrier(post_op.GetComm());

    // ComplexVector port_excitation_E_with_metal;
    // port_excitation_E_with_metal.UseDevice(true);
    // space_op.GetLumpedPortExcitationVector(port_idx, port_excitation_E_with_metal,
    // false); prom_op.UpdatePROM(false, 0.0, port_excitation_E_with_metal);
    // GF->Real().SetFromTrueDofs(port_excitation_E_with_metal.Real());
    // GF->Real().ExchangeFaceNbrData();
    // v_self = port_data.GetVoltage(*GF);
    // Mpi::Print("GetVoltage from port_excitation_E_with_metal (zero_metal fasle) {} {}\n",
    //            v_self.real(), v_self.imag());
    // Mpi::Barrier(post_op.GetComm());

    // {
    //   ComplexVector tmp(port_excitation_E);
    //   tmp.UseDevice(true);
    //   tmp.Real() -= port_excitation_voltage_form;

    //   post_op.SetEGridFunction(tmp, false);
    //   post_op.SetBGridFunction(tmp, false);
    //   post_op.WriteFields(200 + port_idx, 200.0 + port_idx);
    //   Mpi::Barrier(post_op.GetComm());
    // }
    // {
    //   ComplexVector tmp(port_excitation_E);
    //   tmp.UseDevice(true);
    //   tmp.Real() -= port_excitation_E_with_metal.Real();
    //   tmp.Imag() -= port_excitation_E_with_metal.Imag();

    //   post_op.SetEGridFunction(tmp, false);
    //   post_op.SetBGridFunction(tmp, false);
    //   post_op.WriteFields(300 + port_idx, 300.0 + port_idx);
    //   Mpi::Barrier(post_op.GetComm());
    // }
  }

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](double omega)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(true, omega, E);

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
  int it = 2;
  int it0 = it;
  int memory = 0;
  std::vector<double> max_errors = {0.0, 0.0};
  for (; it < max_size; it++)
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
      ++memory;
    }
    else
    {
      memory = 0;
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
    if (memory == convergence_memory)  // converged
    {
      break;
    }
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

    auto excitation_helper = space_op.BuildPortExcitationHelper();
    post_results.PostprocessStep(iodata, post_op, space_op, excitation_helper, step, omega,
                                 excitation_helper.begin()->first, E_elec, E_mag,
                                 (step == n_step - 1) ? &indicator : nullptr);

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

  // Print out PROM
  if (root)
  {
    Eigen::IOFormat eigenio(Eigen::FullPrecision, 0, "; ");

    auto output_kr = OutputFile(post_dir + "prom-Kr.csv", false);
    output_kr.print("{}", fmt::streamed(prom_op.Kr.format(eigenio)));

    auto output_mr = OutputFile(post_dir + "prom-Mr.csv", false);
    output_mr.print("{}", fmt::streamed(prom_op.Mr.format(eigenio)));

    auto output_cr = OutputFile(post_dir + "prom-Cr.csv", false);
    output_cr.print("{}", fmt::streamed(prom_op.Cr.format(eigenio)));

    auto output_h = OutputFile(post_dir + "prom-voltage_norm_H.csv", false);
    output_h.print("{}", fmt::streamed(prom_op.voltage_norm_H.format(eigenio)));

    // Print pivot-points
    // Scale Lc, tc / scaling
  }

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

// -----------------
// Measurements / Postprocessing

PostprocessCurrentsHelper::PostprocessCurrentsHelper(
    bool do_measurement, bool root, const std::string &post_dir,
    const SurfaceCurrentOperator &surf_j_op, int n_expected_rows)
  : do_measurement_{do_measurement}, root_{root}
{
  do_measurement_ = do_measurement_             //
                    && post_dir.length() > 0    // Valid output dir
                    && (surf_j_op.Size() > 0);  // Needs surface currents

  if (!do_measurement_ || !root_)
  {
    return;
  }
  surface_I = TableWithCSVFile(post_dir + "surface-I.csv");
  surface_I.table.reserve(n_expected_rows, surf_j_op.Size());
  surface_I.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));
  for (const auto &[idx, data] : surf_j_op)
  {
    surface_I.table.insert_column(fmt::format("I_{}", idx),
                                  fmt::format("I_inc[{}] (A)", idx));
  }
  surface_I.AppendHeader();
}

void PostprocessCurrentsHelper::AddMeasurement(double omega,
                                               const SurfaceCurrentOperator &surf_j_op,
                                               const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  surface_I.table["idx"] << iodata.DimensionalizeValue(VT::FREQUENCY, omega);
  for (const auto &[idx, data] : surf_j_op)
  {
    auto I_inc = data.GetExcitationCurrent();
    surface_I.table[format("I_{}", idx)] << iodata.DimensionalizeValue(VT::CURRENT, I_inc);
  }
  surface_I.AppendRow();
}

PostprocessPortsHelper::PostprocessPortsHelper(bool do_measurement, bool root,
                                               const std::string &post_dir,
                                               const LumpedPortOperator &lumped_port_op,
                                               int n_expected_rows)
  : do_measurement_{do_measurement}, root_{root}
{
  do_measurement_ = do_measurement_                  //
                    && post_dir.length() > 0         // Valid output dir
                    && (lumped_port_op.Size() > 0);  // Only works for lumped ports

  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;
  port_V = TableWithCSVFile(post_dir + "port-V.csv");
  port_V.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_V.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  port_I = TableWithCSVFile(post_dir + "port-I.csv");
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

void PostprocessPortsHelper::AddMeasurement(double omega, const PostOperator &post_op,
                                            const LumpedPortOperator &lumped_port_op,
                                            const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;

  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  auto freq = iodata.DimensionalizeValue(VT::FREQUENCY, omega);
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

    std::complex<double> V_i = post_op.GetPortVoltage(lumped_port_op, idx);
    std::complex<double> I_i = post_op.GetPortCurrent(lumped_port_op, idx);

    port_V.table[fmt::format("re{}", idx)] << V_i.real() * unit_V;
    port_V.table[fmt::format("im{}", idx)] << V_i.imag() * unit_V;

    port_I.table[fmt::format("re{}", idx)] << I_i.real() * unit_A;
    port_I.table[fmt::format("im{}", idx)] << I_i.imag() * unit_A;
  }
  port_V.AppendRow();
  port_I.AppendRow();
}

PostprocessSParametersHelper::PostprocessSParametersHelper(
    bool do_measurement, bool root, const std::string &post_dir,
    const PortExcitationHelper &port_excitations, const LumpedPortOperator &lumped_port_op,
    const WavePortOperator &wave_port_op, int n_expected_rows)
  : do_measurement_{do_measurement}, root_{root}
{
  do_measurement_ = do_measurement_             //
                    && (post_dir.length() > 0)  // valid output dir
                    &&
                    (port_excitations.Size() == 1);  // handle only one excitation for now

  auto ex = port_excitations.begin()->second;
  do_measurement_ = do_measurement_                              //
                    && (ex.current_port.empty())                 // no current ports
                    && (ex.flatten_port_indices().size() == 1);  // only one port excited

  if (!do_measurement_ || !root_)
  {
    return;
  }

  using fmt::format;
  port_S = TableWithCSVFile(post_dir + "port-S.csv");
  port_S.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_S.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  // TODO: Put this in excitation helper
  std::vector<int> all_port_indices;
  for (const auto &[idx, data] : lumped_port_op)
  {
    all_port_indices.emplace_back(idx);
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    all_port_indices.emplace_back(idx);
  }

  for (const auto &[i_idx, exitation] : port_excitations)
  {
    for (const auto &o_idx : all_port_indices)
    {
      port_S.table.insert_column(format("abs_{}_{}", o_idx, i_idx),
                                 format("|S[{}][{}]| (dB)", o_idx, i_idx));
      port_S.table.insert_column(format("arg_{}_{}", o_idx, i_idx),
                                 format("arg(S[{}][{}]) (deg.)", o_idx, i_idx));
    }
  }
  port_S.AppendHeader();
}

void PostprocessSParametersHelper::AddMeasurement(
    double omega, int excitation_idx, const PostOperator &post_op,
    const PortExcitationHelper &port_excitations, const LumpedPortOperator &lumped_port_op,
    const WavePortOperator &wave_port_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Assume that we are adding same frequencies for excitations and in same order
  // column-wise. Add frequencies on first excitation pass.
  if (port_excitations.begin()->first == excitation_idx)
  {
    port_S.table["idx"] << iodata.DimensionalizeValue(VT::FREQUENCY, omega);
  }

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
    auto S_ij = post_op.GetSParameter(lumped_port_op, wave_port_op, port_excitations, o_idx,
                                      excitation_idx);
    auto abs_S_ij = 20.0 * std::log10(std::abs(S_ij));
    auto arg_S_ij = std::arg(S_ij) * 180.8 / M_PI;

    port_S.table[format("abs_{}_{}", o_idx, excitation_idx)] << abs_S_ij;
    port_S.table[format("arg_{}_{}", o_idx, excitation_idx)] << arg_S_ij;

    Mpi::Print(" {sij} = {:+.3e}{:+.3e}i, |{sij}| = {:+.3e}, arg({sij}) = {:+.3e}\n",
               S_ij.real(), S_ij.imag(), abs_S_ij, arg_S_ij,
               fmt::arg("sij", format("S[{}][{}]", o_idx, excitation_idx)));
  }
  // Regenerate from scratch each time since not row-wise (TODO: improve)
  port_S.WriteFullTableTrunc();
}

DrivenSolver::PostprocessResults::PostprocessResults(
    bool root, const std::string &post_dir, const SpaceOperator &space_op,
    const PortExcitationHelper &port_excitations, int n_expected_rows)
  : currents{true, root, post_dir, space_op.GetSurfaceCurrentOp(), n_expected_rows},
    ports{true, root, post_dir, space_op.GetLumpedPortOp(), n_expected_rows},
    s_parameters{true,
                 root,
                 post_dir,
                 port_excitations,
                 space_op.GetLumpedPortOp(),
                 space_op.GetWavePortOp(),
                 n_expected_rows}
{
}

void DrivenSolver::PostprocessResults::PostprocessStep(
    const IoData &iodata, const PostOperator &post_op, const SpaceOperator &space_op,
    const PortExcitationHelper &port_excitations, int step, double omega,
    int excitation_idx, double E_elec, double E_mag, const ErrorIndicator *indicator)
{
  currents.AddMeasurement(omega, space_op.GetSurfaceCurrentOp(), iodata);
  ports.AddMeasurement(omega, post_op, space_op.GetLumpedPortOp(), iodata);
  s_parameters.AddMeasurement(omega, excitation_idx, post_op, port_excitations,
                              space_op.GetLumpedPortOp(), space_op.GetWavePortOp(), iodata);
}

}  // namespace palace
