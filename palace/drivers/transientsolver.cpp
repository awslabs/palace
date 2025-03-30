// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "transientsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitationhelper.hpp"
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
  auto excitation_helper = space_op.GetPortExcitationHelper();

  MFEM_VERIFY(!excitation_helper.Empty(),
              "No excitation specified for transient simulation!");
  MFEM_VERIFY(
      excitation_helper.Size() == 1,
      "Transient solver currently only supports a single excitation per simulation!");

  TimeOperator time_op(iodata, space_op, dJdt_coef);

  double delta_t = iodata.solver.transient.delta_t;
  int n_step = GetNumSteps(0.0, iodata.solver.transient.max_t, delta_t);
  SaveMetadata(space_op.GetNDSpaces());

  // Time stepping is uniform in the time domain. Index sets are for computing things like
  // port voltages and currents in postprocessing.
  PostOperator<config::ProblemData::Type::TRANSIENT> post_op(iodata, space_op);

  Mpi::Print("\nComputing transient response for:\n{}", excitation_helper.FmtLog());

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<Vector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main time integration loop.
  double t = -delta_t;
  auto t0 = Timer::Now();
  for (int step = 0; step < n_step; step++)
  {
    const double ts = iodata.units.Dimensionalize<Units::ValueType::TIME>(t + delta_t);
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
    Mpi::Print(" Sol. ||E|| = {:.6e}, ||B|| = {:.6e}\n",
               linalg::Norml2(space_op.GetComm(), E),
               linalg::Norml2(space_op.GetComm(), B));

    auto total_domain_energy = post_op.MeasureAndPrintAll(step, E, B, t, J_coef(t));

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
  }
  // Final postprocessing & printing
  BlockTimer bt1(Timer::POSTPRO);
  time_op.PrintStats();
  SaveMetadata(time_op.GetLinearSolver());
  post_op.MeasureFinalize(indicator);
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

}  // namespace palace
