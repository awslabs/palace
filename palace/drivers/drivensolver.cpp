// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <mfem.hpp>
#include "fem/lumpedportoperator.hpp"
#include "fem/postoperator.hpp"
#include "fem/romoperator.hpp"
#include "fem/spaceoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"
#include "fem/waveportoperator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

BaseSolver::ErrorIndicators
DrivenSolver::Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh, Timer &timer,
                    int iter) const
{
  // Set up the spatial discretization and frequency sweep.
  timer.Lap();
  SpaceOperator spaceop(iodata, mesh);
  int nstep = GetNumSteps(iodata.solver.driven.min_f, iodata.solver.driven.max_f,
                          iodata.solver.driven.delta_f);
  int step0 = (iodata.solver.driven.rst > 0) ? iodata.solver.driven.rst - 1 : 0;
  double delta_omega = iodata.solver.driven.delta_f;
  double omega0 = iodata.solver.driven.min_f + step0 * delta_omega;
  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  if (adaptive && nstep <= 2)
  {
    Mpi::Warning("Adaptive frequency sweep requires > 2 total frequency samples!\n"
                 "Reverting to uniform sweep!\n");
    adaptive = false;
  }
  SaveMetadata(spaceop.GetNDSpace());

  // Frequencies will be sampled uniformly in the frequency domain. Index sets are for
  // computing things like S-parameters in postprocessing.
  PostOperator postop(iodata, spaceop, "driven");
  {
    Mpi::Print("\nComputing {}frequency response for:\n", adaptive ? "adaptive fast " : "");
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
    for (const auto &[idx, data] : spaceop.GetWavePortOp())
    {
      if (data.IsExcited())
      {
        if (first)
        {
          Mpi::Print(" Wave port excitation specified on port{}",
                     (spaceop.GetWavePortOp().Size() > 1) ? "s" : "");
          first = false;
        }
        Mpi::Print(" {:d}", idx);
      }
    }
    excitations += first;
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
    MFEM_VERIFY(excitations > 0, "No excitation specified for driven simulation!");
  }
  Mpi::Print("\n");

  // Main frequency sweep loop.
  return adaptive ? SweepAdaptive(spaceop, postop, nstep, step0, omega0, delta_omega, timer)
                  : SweepUniform(spaceop, postop, nstep, step0, omega0, delta_omega, timer);
}

BaseSolver::ErrorIndicators
DrivenSolver::SweepUniform(SpaceOperator &spaceop, PostOperator &postop, int nstep,
                           int step0, double omega0, double delta_omega, Timer &timer) const
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogenous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver:: SetOperators). Compute everything at the first frequency step.
  std::unique_ptr<petsc::PetscParMatrix> A = spaceop.GetSystemMatrixPetsc(
      SpaceOperator::OperatorType::COMPLETE, omega0, mfem::Operator::DIAG_ONE);
  std::unique_ptr<petsc::PetscParMatrix> NegCurl = spaceop.GetNegCurlMatrixPetsc();

  // Set up the linear solver and set operators for the first frequency step. The
  // preconditioner for the complex linear system is constructed from a real approximation
  // to the complex system matrix.
  std::vector<std::unique_ptr<mfem::Operator>> P, AuxP;
  spaceop.GetPreconditionerMatrix(omega0, P, AuxP);

  KspPreconditioner pc(iodata, spaceop.GetDbcMarker(), spaceop.GetNDSpaces(),
                       &spaceop.GetH1Spaces());
  pc.SetOperator(P, &AuxP);

  KspSolver ksp(A->GetComm(), iodata, "ksp_");
  ksp.SetPreconditioner(pc);
  ksp.SetOperator(*A);

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  petsc::PetscParVector RHS(*NegCurl), E(*NegCurl), B(*NegCurl, true);
  E.SetZero();
  B.SetZero();
  timer.construct_time += timer.Lap();

  // Main frequency sweep loop.
  double omega = omega0;
  auto t0 = timer.Now();
  for (int step = step0; step < nstep; ++step)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               nstep, freq, Timer::Duration(timer.Now() - t0).count());

    // Assemble the linear system and solve.
    if (step > step0)
    {
      // Update frequency-dependent excitation and operators.
      A = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::COMPLETE, omega,
                                       mfem::Operator::DIAG_ONE, false);
      spaceop.GetPreconditionerMatrix(omega, P, AuxP, false);
      pc.SetOperator(P, &AuxP);
      ksp.SetOperator(*A);
    }
    spaceop.GetFreqDomainExcitationVector(omega, RHS);
    timer.construct_time += timer.Lap();

    Mpi::Print("\n");
    ksp.Mult(RHS, E);
    timer.solve_time += timer.Lap();

    double E_elec = 0.0, E_mag = 0.0;
    PostOperator::GetBField(omega, *NegCurl, E, B);
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(), omega);
    // E.Print();
    Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n", E.Norml2(), RHS.Norml2());
    if (!iodata.solver.driven.only_port_post)
    {
      E_elec = postop.GetEFieldEnergy();
      E_mag = postop.GetHFieldEnergy();
      Mpi::Print(" Field energy E ({:.3e}) + H ({:.3e}) = {:.3e}\n", E_elec, E_mag,
                 E_elec + E_mag);
    }

    // Postprocess S-parameters and optionally write solution to disk.
    const auto io_time_prev = timer.io_time;
    Postprocess(post_dir_, postop, spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(),
                spaceop.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                !iodata.solver.driven.only_port_post, timer);
    timer.postpro_time += timer.Lap() - (timer.io_time - io_time_prev);

    // Increment frequency.
    omega += delta_omega;
  }
  SaveMetadata(ksp.GetTotalNumMult(), ksp.GetTotalNumIter());

  return BaseSolver::ErrorIndicators();
}

BaseSolver::ErrorIndicators DrivenSolver::SweepAdaptive(SpaceOperator &spaceop,
                                                        PostOperator &postop, int nstep,
                                                        int step0, double omega0,
                                                        double delta_omega,
                                                        Timer &timer) const
{
  // Configure default parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int nmax = iodata.solver.driven.adaptive_nmax;
  int ncand = iodata.solver.driven.adaptive_ncand;
  MFEM_VERIFY(nmax <= 0 || nmax > 2,
              "Adaptive frequency sweep must sample at least two frequency points!");
  if (nmax <= 0)
  {
    nmax = 20;  // Default value
  }
  nmax = std::min(nmax, nstep - step0);  // Maximum number sample points dictated by sweep
  if (ncand > 0)
  {
    if (ncand > nstep - step0)
    {
      Mpi::Warning("Requested candidate points {:d} > number of total frequency sweep "
                   "samples {:d}!\n"
                   "Resetting to the smaller value!\n",
                   ncand, nstep - step0);
      ncand = nstep - step0;
    }
  }
  else
  {
    constexpr int inc = 5;
    ncand = (nstep - step0 + inc - 1) / inc;  // Default value, always >= 1
  }

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  std::unique_ptr<petsc::PetscParMatrix> NegCurl = spaceop.GetNegCurlMatrixPetsc();
  petsc::PetscParVector E(*NegCurl), B(*NegCurl, true);
  E.SetZero();
  B.SetZero();

  // Configure the PROM operator which performs the parameter space sampling and basis
  // construction during the offline phase as well as the PROM solution during the online
  // phase. Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S.
  timer.construct_time += timer.Lap();
  Timer local_timer;
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             nstep - step0,
             iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega0),
             iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY,
                                        omega0 + (nstep - step0 - 1) * delta_omega));
  spaceop.GetWavePortOp().SetSuppressOutput(true);  // Suppress wave port stuff for offline
  RomOperator prom(iodata, spaceop, nmax);
  prom.Initialize(nstep - step0, omega0, delta_omega);
  local_timer.construct_time += local_timer.Lap();

  prom.SolveHDM(omega0, E, true);  // Print matrix stats at first HDM solve
  prom.SolveHDM(omega0 + (nstep - step0 - 1) * delta_omega, E, false);
  local_timer.solve_time += local_timer.Lap();

  // Greedy procedure for basis construction (offline phase). Basis is initialized with
  // solutions at frequency sweep endpoints.
  int iter = static_cast<int>(prom.GetSampleFrequencies().size()), iter0 = iter;
  double max_error = 1.0;
  while (true)
  {
    // Compute maximum error in parameter domain with current PROM.
    double omega_star;
    max_error = prom.ComputeMaxError(ncand, omega_star);
    local_timer.construct_time += local_timer.Lap();
    if (max_error < offline_tol || iter == nmax)
    {
      break;
    }

    // Sample HDM and add solution to basis.
    Mpi::Print(
        "\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = {:.3e}\n",
        iter - iter0 + 1, prom.GetReducedDimension(),
        iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega_star), omega_star,
        max_error);
    prom.SolveHDM(omega_star, E);
    local_timer.solve_time += local_timer.Lap();
    iter++;
  }
  {
    std::vector<double> samples(prom.GetSampleFrequencies());
    // samples.Sort();
    for (auto &sample : samples)
    {
      sample = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, sample);
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}\n",
               (iter == nmax) ? " reached maximum" : " converged with", iter,
               prom.GetReducedDimension(), max_error, offline_tol);
    utils::PrettyPrint(samples, " Sampled frequencies (GHz):");
  }
  SaveMetadata(prom.GetTotalKspMult(), prom.GetTotalKspIter());
  const auto local_construction_time = timer.Lap();
  timer.construct_time += local_construction_time;
  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n"
             " Parameter space sampling: {:.2e} s, HDM solves: {:.2e} s\n",
             Timer::Duration(local_construction_time).count(),
             Timer::Duration(local_timer.construct_time).count(),
             Timer::Duration(local_timer.solve_time).count());  // Timings on rank 0

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  spaceop.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  int step = step0;
  double omega = omega0;
  auto t0 = timer.Now();
  while (step < nstep)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               nstep, freq, Timer::Duration(timer.Now() - t0).count());

    // Assemble the linear system and solve.
    prom.AssemblePROM(omega);
    timer.construct_time += timer.Lap();

    Mpi::Print("\n");
    prom.SolvePROM(E);
    timer.solve_time += timer.Lap();

    double E_elec = 0.0, E_mag = 0.0;
    PostOperator::GetBField(omega, *NegCurl, E, B);
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(), omega);
    // E.Print();
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", E.Norml2());
    if (!iodata.solver.driven.only_port_post)
    {
      E_elec = postop.GetEFieldEnergy();
      E_mag = postop.GetHFieldEnergy();
      Mpi::Print(" Field energy E ({:.3e}) + H ({:.3e}) = {:.3e}\n", E_elec, E_mag,
                 E_elec + E_mag);
    }

    // Postprocess S-parameters and optionally write solution to disk.
    const auto io_time_prev = timer.io_time;
    Postprocess(post_dir_, postop, spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(),
                spaceop.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                !iodata.solver.driven.only_port_post, timer);
    timer.postpro_time += timer.Lap() - (timer.io_time - io_time_prev);

    // Increment frequency.
    step++;
    omega += delta_omega;
  }

  return BaseSolver::ErrorIndicators();
}

int DrivenSolver::GetNumSteps(double start, double end, double delta) const
{
  MFEM_VERIFY(delta != 0.0, "Zero frequency step is not allowed!");
  constexpr double delta_eps = 1.0e-9;  // 9 digits of precision comparing endpoint
  double dnfreq = std::abs(end - start) / std::abs(delta);
  int nstep = 1 + static_cast<int>(dnfreq);
  double dfinal = start + nstep * delta;
  return nstep + ((delta < 0.0 && dfinal - end > -delta_eps * end) ||
                  (delta > 0.0 && dfinal - end < delta_eps * end));
}

void DrivenSolver::Postprocess(const std::string &post_dir, const PostOperator &postop,
                               const LumpedPortOperator &lumped_port_op,
                               const WavePortOperator &wave_port_op,
                               const SurfaceCurrentOperator &surf_j_op, int step,
                               double omega, double E_elec, double E_mag, bool full,
                               Timer &timer) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main frequency sweep loop.
  double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
  PostprocessCurrents(post_dir, postop, surf_j_op, step, omega);
  PostprocessPorts(post_dir, postop, lumped_port_op, step, omega);
  if (surf_j_op.Size() == 0)
  {
    PostprocessSParameters(post_dir, postop, lumped_port_op, wave_port_op, step, omega);
  }
  if (full)
  {
    double E_cap = postop.GetLumpedCapacitorEnergy(lumped_port_op);
    double E_ind = postop.GetLumpedInductorEnergy(lumped_port_op);
    PostprocessDomains(post_dir, postop, "f (GHz)", step, freq, E_elec, E_mag, E_cap,
                       E_ind);
    PostprocessSurfaces(post_dir, postop, "f (GHz)", step, freq, E_elec + E_cap,
                        E_mag + E_ind, 1.0, 1.0);
    PostprocessProbes(post_dir, postop, "f (GHz)", step, freq);
  }
  if (iodata.solver.driven.delta_post > 0 && step % iodata.solver.driven.delta_post == 0)
  {
    auto t0 = timer.Now();
    Mpi::Print("\n");
    PostprocessFields(post_dir, postop, step / iodata.solver.driven.delta_post, freq);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
    timer.io_time += timer.Now() - t0;
  }
}

namespace
{

struct CurrentData
{
  int idx;      // Current source index
  double Iinc;  // Excitation current
};

struct PortVIData
{
  int idx;                      // Port index
  bool excitation;              // Flag for excited ports
  double Vinc, Iinc;            // Incident voltage, current
  std::complex<double> Vi, Ii;  // Port voltage, current
};

struct PortSData
{
  int idx;                   // Port index
  std::complex<double> Sij;  // Scattering parameter
};

}  // namespace

void DrivenSolver::PostprocessCurrents(const std::string &post_dir,
                                       const PostOperator &postop,
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
    const double Iinc = data.GetExcitationCurrent();
    j_data.push_back({idx, Iinc});
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
                     "Iinc[" + std::to_string(data.idx) + "]", table.w,
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
                   data.Iinc, table.w, table.p,
                   (data.idx == j_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void DrivenSolver::PostprocessPorts(const std::string &post_dir, const PostOperator &postop,
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
    const double Vinc = data.GetExcitationVoltage();
    const double Iinc = (std::abs(Vinc) > 0.0) ? data.GetExcitationPower() / Vinc : 0.0;
    const std::complex<double> Vi = postop.GetPortVoltage(lumped_port_op, idx);
    const std::complex<double> Ii = postop.GetPortCurrent(lumped_port_op, idx);
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
        output.print("{:>{}s},", "f (GHz)", table.w1);
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
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{V[" + std::to_string(data.idx) + "]}", table.w,
                       "Im{V[" + std::to_string(data.idx) + "]}", table.w,
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
                       data.Vinc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.Vi.real(), table.w, table.p,
                     data.Vi.imag(), table.w, table.p,
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
                         "I_inc[" + std::to_string(data.idx) + "]", table.w);
            // clang-format on
          }
        }
        for (const auto &data : port_data)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{I[" + std::to_string(data.idx) + "]}", table.w,
                       "Im{I[" + std::to_string(data.idx) + "]}", table.w,
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
                       data.Iinc, table.w, table.p);
          // clang-format on
        }
      }
      for (const auto &data : port_data)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.Ii.real(), table.w, table.p,
                     data.Ii.imag(), table.w, table.p,
                     (data.idx == port_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
  }
}

void DrivenSolver::PostprocessSParameters(const std::string &post_dir,
                                          const PostOperator &postop,
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
    if (data.IsExcited())
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
    if (data.IsExcited())
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
      const std::complex<double> Sij =
          postop.GetSParameter(lumped_port_op, idx, source_idx);
      port_data.push_back({idx, Sij});
    }
  }
  else  // src_wave_port
  {
    // Compute wave port S-parameters.
    for (const auto &[idx, data] : wave_port_op)
    {
      const std::complex<double> Sij = postop.GetSParameter(wave_port_op, idx, source_idx);
      port_data.push_back({idx, Sij});
    }
  }

  // Print table to stdout.
  for (const auto &data : port_data)
  {
    std::string str =
        "S[" + std::to_string(data.idx) + "][" + std::to_string(source_idx) + "]";
    // clang-format off
    Mpi::Print(" {} = {:+.3e}{:+.3e}i, |{}| = {:+.3e}, arg({}) = {:+.3e}\n",
               str, data.Sij.real(), data.Sij.imag(),
               str, 20.0 * std::log10(std::abs(data.Sij)),
               str, std::arg(data.Sij) * 180.0 / M_PI);
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
      output.print("{:+{}.{}e},{:+{}.{}e}{}",
                   20.0 * std::log10(std::abs(data.Sij)), table.w, table.p,
                   std::arg(data.Sij) * 180.0 / M_PI, table.w, table.p,
                   (data.idx == port_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

}  // namespace palace
