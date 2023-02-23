// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "transientsolver.hpp"

#include <mfem.hpp>
#include "fem/lumpedportoperator.hpp"
#include "fem/postoperator.hpp"
#include "fem/spaceoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"
#include "fem/timeoperator.hpp"
#include "utils/communication.hpp"
#include "utils/excitations.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

BaseSolver::SolveOutput
TransientSolver::Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh, Timer &timer,
                       int iter) const
{
  // Set up the spatial discretization and time integrators for the E and B fields.
  timer.Lap();
  std::function<double(double)> J_coef = GetTimeExcitation(false);
  std::function<double(double)> dJdt_coef = GetTimeExcitation(true);
  SpaceOperator spaceop(iodata, mesh);
  TimeOperator timeop(iodata, spaceop, dJdt_coef);
  double delta_t = iodata.solver.transient.delta_t;
  if (timeop.isExplicit())
  {
    // Stability limited time step.
    const double dt_max = timeop.GetMaxTimeStep();
    const double dts_max = iodata.DimensionalizeValue(IoData::ValueType::TIME, dt_max);
    Mpi::Print(" Maximum stable time step: {:.6e} ns\n", dts_max);
    delta_t = std::min(delta_t, 0.95 * dt_max);
  }
  int nstep = GetNumSteps(0.0, iodata.solver.transient.max_t, delta_t);
  SaveMetadata(spaceop.GetNDSpace());

  // Time stepping is uniform in the time domain. Index sets are for computing things like
  // port voltages and currents in postprocessing.
  PostOperator postop(iodata, spaceop, "transient");
  {
    Mpi::Print("\nComputing transient response for:\n");
    bool first = true;
    for (const auto &[idx, data] : spaceop.GetLumpedPortOp())
    {
      if (data.IsExcited())
      {
        if (first)
        {
          Mpi::Print(" Lumped port excitation specified on port{}",
                     (spaceop.GetLumpedPortOp().Size() > 1) ? "s" : "");
          first = false;
        }
        Mpi::Print(" {:d}", idx);
      }
    }
    int excitations = first;
    first = true;
    for (const auto &[idx, data] : spaceop.GetSurfaceCurrentOp())
    {
      if (first)
      {
        Mpi::Print(" Surface current excitation specified on port{}",
                   (spaceop.GetSurfaceCurrentOp().Size() > 1) ? "s" : "");
        first = false;
      }
      Mpi::Print(" {:d}", idx);
    }
    excitations += first;
    MFEM_VERIFY(excitations > 0, "No excitation specified for transient simulation!");
  }
  Mpi::Print("\n");
  timer.construct_time += timer.Lap();

  // Main time integration loop.
  int step = 0;
  double t = -delta_t;
  auto t0 = timer.Now();
  while (step < nstep)
  {
    const double ts = iodata.DimensionalizeValue(IoData::ValueType::TIME, t + delta_t);
    Mpi::Print("\nIt {:d}/{:d}: t = {:e} ns (elapsed time = {:.2e} s)\n", step, nstep - 1,
               ts, Timer::Duration(timer.Now() - t0).count());

    // Single time step t => t + dt.
    if (step == 0)
    {
      t += delta_t;
      timeop.Init();  // Initial conditions
    }
    else
    {
      timeop.Step(t, delta_t);  // Advances t internally
    }
    timer.solve_time += timer.Lap();

    double E_elec = 0.0, E_mag = 0.0;
    const mfem::Vector &E = timeop.GetE();
    const mfem::Vector &B = timeop.GetB();
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp());
    // E.Print();
    Mpi::Print(" Sol. ||E|| = {:.6e}, ||B|| = {:.6e}\n",
               std::sqrt(mfem::InnerProduct(mesh.back()->GetComm(), E, E)),
               std::sqrt(mfem::InnerProduct(mesh.back()->GetComm(), B, B)));
    if (!iodata.solver.transient.only_port_post)
    {
      E_elec = postop.GetEFieldEnergy();
      E_mag = postop.GetHFieldEnergy();
      Mpi::Print(" Field energy E ({:.3e}) + H ({:.3e}) = {:.3e}\n", E_elec, E_mag,
                 E_elec + E_mag);
    }

    // Postprocess port voltages/currents and optionally write solution to disk.
    const auto io_time_prev = timer.io_time;
    Postprocess(post_dir_, postop, spaceop.GetLumpedPortOp(), spaceop.GetSurfaceCurrentOp(),
                step, t, J_coef(t), E_elec, E_mag, !iodata.solver.transient.only_port_post,
                timer);
    timer.postpro_time += timer.Lap() - (timer.io_time - io_time_prev);

    // Increment time step.
    step++;
  }
  SaveMetadata(timeop.GetTotalKspMult(), timeop.GetTotalKspIter());

  return BaseSolver::SolveOutput();
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
    default:
      MFEM_ABORT("Unsupported source excitation type!");
  }
  return F{};
}

int TransientSolver::GetNumSteps(double start, double end, double delta) const
{
  MFEM_VERIFY(delta > 0.0, "Zero time step is not allowed!");
  constexpr double delta_eps = 1.0e-9;  // 9 digits of precision comparing endpoint
  double dnfreq = std::abs(end - start) / std::abs(delta);
  int nstep = 1 + static_cast<int>(dnfreq);
  double dfinal = start + nstep * delta;
  return nstep + ((delta < 0.0 && dfinal - end > -delta_eps * end) ||
                  (delta > 0.0 && dfinal - end < delta_eps * end));
}

void TransientSolver::Postprocess(const std::string &post_dir, const PostOperator &postop,
                                  const LumpedPortOperator &lumped_port_op,
                                  const SurfaceCurrentOperator &surf_j_op, int step,
                                  double t, double J_coef, double E_elec, double E_mag,
                                  bool full, Timer &timer) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main time integration loop.
  const double ts = iodata.DimensionalizeValue(IoData::ValueType::TIME, t);
  PostprocessCurrents(post_dir, postop, surf_j_op, step, t, J_coef);
  PostprocessPorts(post_dir, postop, lumped_port_op, step, t, J_coef);
  if (full)
  {
    double E_cap = postop.GetLumpedCapacitorEnergy(lumped_port_op);
    double E_ind = postop.GetLumpedInductorEnergy(lumped_port_op);
    PostprocessDomains(post_dir, postop, "t (ns)", step, ts, E_elec, E_mag, E_cap, E_ind);
    PostprocessSurfaces(post_dir, postop, "t (ns)", step, ts, E_elec + E_cap, E_mag + E_ind,
                        1.0, 1.0);
    PostprocessProbes(post_dir, postop, "t (ns)", step, ts);
  }
  if (iodata.solver.transient.delta_post > 0 &&
      step % iodata.solver.transient.delta_post == 0)
  {
    auto t0 = timer.Now();
    Mpi::Print("\n");
    PostprocessFields(post_dir, postop, step / iodata.solver.transient.delta_post, ts);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step);
    timer.io_time += timer.Now() - t0;
  }
}

namespace
{

struct CurrentData
{
  const int idx;      // Current source index
  const double Iinc;  // Excitation current
};

struct PortData
{
  const int idx;            // Port index
  const bool excitation;    // Flag for excited ports
  const double Vinc, Iinc;  // Incident voltage, current
  const double Vi, Ii;      // Port voltage, current
};

}  // namespace

void TransientSolver::PostprocessCurrents(const std::string &post_dir,
                                          const PostOperator &postop,
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
    const double Iinc = data.GetExcitationCurrent() * J_coef;  // Iinc(t) = g(t) Iinc
    j_data.push_back({idx, Iinc});
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
                     "Iinc[" + std::to_string(data.idx) + "]", table.w,
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
                   data.Iinc, table.w, table.p,
                   (data.idx == j_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void TransientSolver::PostprocessPorts(const std::string &post_dir,
                                       const PostOperator &postop,
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
    const double Vinc = data.GetExcitationVoltage() * J_coef;  // Vinc(t) = g(t) Vinc
    const double Iinc =
        (std::abs(Vinc) > 0.0) ? data.GetExcitationPower() * J_coef * J_coef / Vinc : 0.0;
    const double Vi = postop.GetPortVoltage(lumped_port_op, idx).real();
    const double Ii = postop.GetPortCurrent(lumped_port_op, idx).real();
    port_data.push_back({idx, data.IsExcited(), Vinc, Iinc, Vi, Ii});
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
                         "V_inc[" + std::to_string(data.idx) + "]", table.w);
            // clang-format on
          }
        }
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s}{}",
                       "V[" + std::to_string(data.idx) + "]", table.w,
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
                       data.Vinc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     data.Vi, table.w, table.p,
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
                         "I_inc[" + std::to_string(data.idx) + "]", table.w);
            // clang-format on
          }
        }
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s}{}",
                       "I[" + std::to_string(data.idx) + "]", table.w,
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
                       data.Iinc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     data.Ii, table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
  }
}

}  // namespace palace
