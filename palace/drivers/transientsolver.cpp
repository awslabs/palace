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

    // Postprocess port voltages/currents and optionally write solution to disk.
    Postprocess(post_op, space_op.GetLumpedPortOp(), space_op.GetSurfaceCurrentOp(), step,
                t, J_coef(t), E_elec, E_mag, (step == n_step - 1) ? &indicator : nullptr);

    // Increment time step.
    step++;
  }
  BlockTimer bt1(Timer::POSTPRO);
  time_op.PrintStats();
  SaveMetadata(time_op.GetLinearSolver());
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

void TransientSolver::Postprocess(const PostOperator &post_op,
                                  const LumpedPortOperator &lumped_port_op,
                                  const SurfaceCurrentOperator &surf_j_op, int step,
                                  double t, double J_coef, double E_elec, double E_mag,
                                  const ErrorIndicator *indicator) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main time integration loop.
  const double ts = iodata.DimensionalizeValue(IoData::ValueType::TIME, t);
  const double E_cap = post_op.GetLumpedCapacitorEnergy(lumped_port_op);
  const double E_ind = post_op.GetLumpedInductorEnergy(lumped_port_op);
  PostprocessCurrents(post_op, surf_j_op, step, t, J_coef);
  PostprocessPorts(post_op, lumped_port_op, step, t, J_coef);
  PostprocessDomains(post_op, "t (ns)", step, ts, E_elec, E_mag, E_cap, E_ind);
  PostprocessSurfaces(post_op, "t (ns)", step, ts, E_elec + E_cap, E_mag + E_ind);
  PostprocessProbes(post_op, "t (ns)", step, ts);
  if (iodata.solver.transient.delta_post > 0 &&
      step % iodata.solver.transient.delta_post == 0)
  {
    Mpi::Print("\n");
    PostprocessFields(post_op, step / iodata.solver.transient.delta_post, ts);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(post_op, *indicator, iodata.solver.transient.delta_post > 0);
  }
}

namespace
{

struct CurrentData
{
  const int idx;       // Current source index
  const double I_inc;  // Excitation current
};

struct PortData
{
  const int idx;              // Port index
  const bool excitation;      // Flag for excited ports
  const double V_inc, I_inc;  // Incident voltage, current
  const double V_i, I_i;      // Port voltage, current
};

}  // namespace

void TransientSolver::PostprocessCurrents(const PostOperator &post_op,
                                          const SurfaceCurrentOperator &surf_j_op, int step,
                                          double t, double J_coef) const
{
  // Postprocess the time domain surface current excitations.
  if (post_dir.length() == 0)
  {
    return;
  }
  std::vector<CurrentData> j_data;
  j_data.reserve(surf_j_op.Size());
  for (const auto &[idx, data] : surf_j_op)
  {
    const double I_inc = data.GetExcitationCurrent() * J_coef;  // I_inc(t) = g(t) I_inc
    j_data.push_back({idx, iodata.DimensionalizeValue(IoData::ValueType::CURRENT, I_inc)});
  }
  if (root && !j_data.empty())
  {
    std::string path = post_dir + "surface-I.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", "t (ns)", table.w1);
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
                 iodata.DimensionalizeValue(IoData::ValueType::TIME, t),
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

void TransientSolver::PostprocessPorts(const PostOperator &post_op,
                                       const LumpedPortOperator &lumped_port_op, int step,
                                       double t, double J_coef) const
{
  // Postprocess the time domain lumped port voltages and currents, which can then be used
  // to compute S- or Z-parameters.
  if (post_dir.length() == 0)
  {
    return;
  }
  std::vector<PortData> port_data;
  port_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    const double V_inc = data.GetExcitationVoltage() * J_coef;  // V_inc(t) = g(t) V_inc
    const double I_inc =
        (std::abs(V_inc) > 0.0) ? data.GetExcitationPower() * J_coef * J_coef / V_inc : 0.0;
    const double V_i = post_op.GetPortVoltage(lumped_port_op, idx).real();
    const double I_i = post_op.GetPortCurrent(lumped_port_op, idx).real();
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
        output.print("{:>{}s},", "t (ns)", table.w1);
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
          output.print("{:>{}s}{}",
                       "V[" + std::to_string(data.idx) + "] (V)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   iodata.DimensionalizeValue(IoData::ValueType::TIME, t),
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
        output.print("{:+{}.{}e}{}",
                     data.V_i, table.w, table.p,
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
        output.print("{:>{}s},", "t (ns)", table.w1);
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
          output.print("{:>{}s}{}",
                       "I[" + std::to_string(data.idx) + "] (A)", table.w,
                       (data.idx == port_data.back().idx) ? "" : ",");
          // clang-format on
        }
        output.print("\n");
      }
      // clang-format off
      output.print("{:{}.{}e},",
                   iodata.DimensionalizeValue(IoData::ValueType::TIME, t),
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
        output.print("{:+{}.{}e}{}",
                     data.I_i, table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
  }
}

}  // namespace palace
