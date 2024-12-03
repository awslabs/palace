// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "transientsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/timeoperator.hpp"
#include "utils/communication.hpp"
#include "utils/excitations.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
TransientSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Set up the spatial discretization and time integrators for the E and B fields.
  BlockTimer bt0(Timer::CONSTRUCT);
  std::function<double(double)> J_coef = GetTimeExcitation(false);
  std::function<double(double)> dJdt_coef = GetTimeExcitation(true);
  SpaceOperator space_op(iodata, mesh);
  TimeOperator time_op(iodata, space_op, dJdt_coef);

  double delta_t = iodata.solver.transient.delta_t;
  int n_step = GetNumSteps(0.0, iodata.solver.transient.max_t, delta_t);
  SaveMetadata(space_op.GetNDSpaces());

  // Time stepping is uniform in the time domain. Index sets are for computing things like
  // port voltages and currents in postprocessing.
  PostOperator post_op(iodata, space_op, "transient");
  PostprocessPrintResults post_results(root, post_dir, post_op, space_op, n_step,
                                       iodata.solver.transient.delta_post);

  {
    Mpi::Print("\nComputing transient response for:\n");
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
    MFEM_VERIFY(excitations > 0, "No excitation specified for transient simulation!");
  }
  Mpi::Print("\n");

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<Vector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main time integration loop.
  int step = 0;
  double t = -delta_t;
  auto t0 = Timer::Now();
  while (step < n_step)
  {
    const double ts = iodata.DimensionalizeValue(IoData::ValueType::TIME, t + delta_t);
    Mpi::Print("\nIt {:d}/{:d}: t = {:e} ns (elapsed time = {:.2e} s)\n", step, n_step - 1,
               ts, Timer::Duration(Timer::Now() - t0).count());

    // Single time step t -> t + dt.
    BlockTimer bt1(Timer::TS);
    if (step == 0)
    {
      Mpi::Print("\n");
      t += delta_t;
      time_op.Init();  // Initial conditions
    }
    else
    {
      time_op.Step(t, delta_t);  // Advances t internally
    }

    // Postprocess for the time step.
    BlockTimer bt2(Timer::POSTPRO);
    const Vector &E = time_op.GetE();
    const Vector &B = time_op.GetB();
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.UpdatePorts(space_op.GetLumpedPortOp());
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    Mpi::Print(" Sol. ||E|| = {:.6e}, ||B|| = {:.6e}\n",
               linalg::Norml2(space_op.GetComm(), E),
               linalg::Norml2(space_op.GetComm(), B));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
                 E_mag * J, (E_elec + E_mag) * J);
    }
    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);

    post_results.PostprocessStep(iodata, post_op, space_op, step, t, J_coef(t), E_elec,
                                 E_mag);
                                 
    // Increment time step.
    step++;
  }
  // Final postprocessing & printing
  BlockTimer bt1(Timer::POSTPRO);
  time_op.PrintStats();
  SaveMetadata(time_op.GetLinearSolver());
  post_results.PostprocessFinal(post_op, indicator);
  return {indicator, space_op.GlobalTrueVSize()};
}

std::function<double(double)> TransientSolver::GetTimeExcitation(bool dot) const
{
  using namespace excitations;
  using F = std::function<double(double)>;
  const config::TransientSolverData &data = iodata.solver.transient;
  const config::TransientSolverData::ExcitationType &type = data.excitation;
  if (type == config::TransientSolverData::ExcitationType::SINUSOIDAL ||
      type == config::TransientSolverData::ExcitationType::MOD_GAUSSIAN)
  {
    MFEM_VERIFY(data.pulse_f > 0.0,
                "Excitation frequency is missing for transient simulation!");
  }
  if (type == config::TransientSolverData::ExcitationType::GAUSSIAN ||
      type == config::TransientSolverData::ExcitationType::DIFF_GAUSSIAN ||
      type == config::TransientSolverData::ExcitationType::MOD_GAUSSIAN ||
      type == config::TransientSolverData::ExcitationType::SMOOTH_STEP)
  {
    MFEM_VERIFY(data.pulse_tau > 0.0,
                "Excitation width is missing for transient simulation!");
  }
  const double delay =
      (type == config::TransientSolverData::ExcitationType::GAUSSIAN ||
       type == config::TransientSolverData::ExcitationType::DIFF_GAUSSIAN ||
       type == config::TransientSolverData::ExcitationType::MOD_GAUSSIAN)
          ? 4.5 * data.pulse_tau
          : 0.0;
  switch (type)
  {
    case config::TransientSolverData::ExcitationType::SINUSOIDAL:
      if (dot)
      {
        return F{[=](double t) { return dpulse_sinusoidal(t, data.pulse_f, delay); }};
      }
      else
      {
        return F{[=](double t) { return pulse_sinusoidal(t, data.pulse_f, delay); }};
      }
      break;
    case config::TransientSolverData::ExcitationType::GAUSSIAN:
      if (dot)
      {
        return F{[=](double t) { return dpulse_gaussian(t, data.pulse_tau, delay); }};
      }
      else
      {
        return F{[=](double t) { return pulse_gaussian(t, data.pulse_tau, delay); }};
      }
      break;
    case config::TransientSolverData::ExcitationType::DIFF_GAUSSIAN:
      if (dot)
      {
        return F{[=](double t) { return dpulse_gaussian_diff(t, data.pulse_tau, delay); }};
      }
      else
      {
        return F{[=](double t) { return pulse_gaussian_diff(t, data.pulse_tau, delay); }};
      }
      break;
    case config::TransientSolverData::ExcitationType::MOD_GAUSSIAN:
      if (dot)
      {
        return F{[=](double t)
                 { return dpulse_gaussian_mod(t, data.pulse_f, data.pulse_tau, delay); }};
      }
      else
      {
        return F{[=](double t)
                 { return pulse_gaussian_mod(t, data.pulse_f, data.pulse_tau, delay); }};
      }
      break;
    case config::TransientSolverData::ExcitationType::RAMP_STEP:
      if (dot)
      {
        return F{[=](double t) { return dpulse_ramp(t, data.pulse_tau, delay); }};
      }
      else
      {
        return F{[=](double t) { return pulse_ramp(t, data.pulse_tau, delay); }};
      }
      break;
    case config::TransientSolverData::ExcitationType::SMOOTH_STEP:
      if (dot)
      {
        return F{[=](double t) { return dpulse_smootherstep(t, data.pulse_tau, delay); }};
      }
      else
      {
        return F{[=](double t) { return pulse_smootherstep(t, data.pulse_tau, delay); }};
      }
      break;
  }
  return F{};
}

int TransientSolver::GetNumSteps(double start, double end, double delta) const
{
  if (end < start)
  {
    return 1;
  }
  MFEM_VERIFY(delta > 0.0, "Zero time step is not allowed!");
  constexpr double delta_eps = 1.0e-9;  // 9 digits of precision comparing endpoint
  double dnfreq = std::abs(end - start) / std::abs(delta);
  int n_step = 1 + static_cast<int>(dnfreq);
  double dfinal = start + n_step * delta;
  return n_step + ((delta < 0.0 && dfinal - end > -delta_eps * end) ||
                   (delta > 0.0 && dfinal - end < delta_eps * end));
}

// -----------------
// Measurements / Postprocessing

TransientSolver::CurrentsPostPrinter::CurrentsPostPrinter(
    bool do_measurement, bool root, const std::string &post_dir,
    const SurfaceCurrentOperator &surf_j_op, int n_expected_rows)
  : root_{root},                               //
    do_measurement_(do_measurement             //
                    && post_dir.length() > 0   // Valid output dir
                    && (surf_j_op.Size() > 0)  // Needs surface currents
    )
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;

  surface_I = TableWithCSVFile(post_dir + "surface-I.csv");
  surface_I.table.reserve(n_expected_rows, surf_j_op.Size());
  surface_I.table.insert_column(Column("idx", "t (ns)", 0, {}, {}, ""));
  for (const auto &[idx, data] : surf_j_op)
  {
    surface_I.table.insert_column(format("I_{}", idx), format("I_inc[{}] (A)", idx));
  }
  surface_I.AppendHeader();
}

void TransientSolver::CurrentsPostPrinter::AddMeasurement(
    double t, double J_coef, const SurfaceCurrentOperator &surf_j_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  surface_I.table["idx"] << iodata.DimensionalizeValue(VT::TIME, t);
  for (const auto &[idx, data] : surf_j_op)
  {
    auto I_inc = data.GetExcitationCurrent() * J_coef;  // I_inc(t) = g(t) I_inc
    surface_I.table[format("I_{}", idx)] << iodata.DimensionalizeValue(VT::CURRENT, I_inc);
  }
  surface_I.AppendRow();
}

TransientSolver::PortsPostPrinter::PortsPostPrinter(
    bool do_measurement, bool root, const std::string &post_dir,
    const LumpedPortOperator &lumped_port_op, int n_expected_rows)
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
  port_V.table.insert_column(Column("idx", "t (ns)", 0, {}, {}, ""));

  port_I = TableWithCSVFile(post_dir + "port-I.csv");
  port_I.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_I.table.insert_column(Column("idx", "t (ns)", 0, {}, {}, ""));

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      port_V.table.insert_column(format("inc{}", idx), format("V_inc[{}] (V)", idx));
      port_I.table.insert_column(format("inc{}", idx), format("I_inc[{}] (A)", idx));
    }
    port_V.table.insert_column(format("re{}", idx), format("V[{}] (V)", idx));
    port_I.table.insert_column(format("re{}", idx), format("I[{}] (A)", idx));
  }
  port_V.AppendHeader();
  port_I.AppendHeader();
}

void TransientSolver::PortsPostPrinter::AddMeasurement(
    double t, double J_coef, const PostOperator &post_op,
    const LumpedPortOperator &lumped_port_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;
  using VT = IoData::ValueType;

  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  auto time = iodata.DimensionalizeValue(VT::TIME, t);
  port_V.table["idx"] << time;
  port_I.table["idx"] << time;

  auto unit_V = iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
  auto unit_A = iodata.DimensionalizeValue(VT::CURRENT, 1.0);

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation)
    {
      double V_inc = data.GetExcitationVoltage() * J_coef;  // V_inc(t) = g(t) V_inc
      double I_inc = (std::abs(V_inc) > 0.0)
                         ? data.GetExcitationPower() * J_coef * J_coef / V_inc
                         : 0.0;

      port_V.table[format("inc{}", idx)] << V_inc * unit_V;
      port_I.table[format("inc{}", idx)] << I_inc * unit_A;
    }

    std::complex<double> V_i = post_op.GetPortVoltage(lumped_port_op, idx);
    std::complex<double> I_i = post_op.GetPortCurrent(lumped_port_op, idx);

    port_V.table[format("re{}", idx)] << V_i.real() * unit_V;
    port_I.table[format("re{}", idx)] << I_i.real() * unit_A;
  }
  port_V.AppendRow();
  port_I.AppendRow();
}

TransientSolver::PostprocessPrintResults::PostprocessPrintResults(
    bool root, const std::string &post_dir, const PostOperator &post_op,
    const SpaceOperator &space_op, int n_expected_rows, int delta_post_)
  : delta_post{delta_post_},
    domains{true, root, post_dir, post_op.GetDomainPostOp(), "t (ns)", n_expected_rows},
    surfaces{true, root, post_dir, post_op, "t (ns)", n_expected_rows},
    currents{true, root, post_dir, space_op.GetSurfaceCurrentOp(), n_expected_rows},
    probes{true, root, post_dir, post_op, "t (ns)", n_expected_rows},
    ports{true, root, post_dir, space_op.GetLumpedPortOp(), n_expected_rows},
    error_indicator{true, root, post_dir}
{
  // If to print paraview fields
  if (delta_post > 0)
  {
    if (post_dir.length() == 0)
    {
      Mpi::Warning(post_op.GetComm(),
                   "No file specified under [\"Problem\"][\"Output\"]!\nSkipping saving of "
                   "fields to disk in solve!\n");
    }
    else
    {
      write_paraview_fields = true;
    }
  }
}

void TransientSolver::PostprocessPrintResults::PostprocessStep(
    const IoData &iodata, const PostOperator &post_op, const SpaceOperator &space_op,
    int step, double t, double J_coef, double E_elec, double E_mag)
{
  auto time = iodata.DimensionalizeValue(IoData::ValueType::TIME, t);
  auto E_cap = post_op.GetLumpedCapacitorEnergy(space_op.GetLumpedPortOp());
  auto E_ind = post_op.GetLumpedInductorEnergy(space_op.GetLumpedPortOp());

  domains.AddMeasurement(time, post_op, E_elec, E_mag, E_cap, E_ind, iodata);
  surfaces.AddMeasurement(time, post_op, E_elec + E_cap, E_mag + E_ind, iodata);
  currents.AddMeasurement(t, J_coef, space_op.GetSurfaceCurrentOp(), iodata);
  probes.AddMeasurement(time, post_op, iodata);
  ports.AddMeasurement(t, J_coef, post_op, space_op.GetLumpedPortOp(), iodata);
  // The internal GridFunctions in PostOperator have already been set:
  if (write_paraview_fields && (step % delta_post == 0))
  {
    Mpi::Print("\n");
    post_op.WriteFields(step / delta_post, time);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
}

void TransientSolver::PostprocessPrintResults::PostprocessFinal(
    const PostOperator &post_op, const ErrorIndicator &indicator)
{
  BlockTimer bt0(Timer::POSTPRO);
  error_indicator.PrintIndicatorStatistics(post_op, indicator);
  if (write_paraview_fields)
  {
    post_op.WriteFieldsFinal(&indicator);
  }
}

}  // namespace palace
