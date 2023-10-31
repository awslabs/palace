// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
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
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
DrivenSolver::Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const
{
  // Set up the spatial discretization and frequency sweep.
  BlockTimer bt0(Timer::CONSTRUCT);
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
  SaveMetadata(spaceop.GetNDSpaces());

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
  return {adaptive ? SweepAdaptive(spaceop, postop, nstep, step0, omega0, delta_omega)
                   : SweepUniform(spaceop, postop, nstep, step0, omega0, delta_omega),
          spaceop.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &spaceop, PostOperator &postop,
                                          int nstep, int step0, double omega0,
                                          double delta_omega) const
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogenous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver::SetOperators). Compute everything at the first frequency step.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto K = spaceop.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = spaceop.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = spaceop.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega0, Operator::DIAG_ZERO);
  const auto &Curl = spaceop.GetCurlMatrix();

  // Set up the linear solver and set operators for the first frequency step. The
  // preconditioner for the complex linear system is constructed from a real approximation
  // to the complex system matrix.
  auto A = spaceop.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega0,
                                   std::complex<double>(-omega0 * omega0, 0.0), K.get(),
                                   C.get(), M.get(), A2.get());
  auto P = spaceop.GetPreconditionerMatrix<ComplexOperator>(1.0, omega0, -omega0 * omega0,
                                                            omega0);

  ComplexKspSolver ksp(iodata, spaceop.GetNDSpaces(), &spaceop.GetH1Spaces());
  ksp.SetOperators(*A, *P);

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  CurlFluxErrorEstimator<ComplexVector> estimator(
      spaceop.GetMaterialOp(), spaceop.GetNDSpace(), iodata.solver.linear.estimator_tol,
      iodata.solver.linear.estimator_max_it, 0);
  ErrorIndicator indicator;

  // Main frequency sweep loop.
  int step = step0;
  double omega = omega0;
  auto t0 = Timer::Now();
  while (step < nstep)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               nstep, freq, Timer::Duration(Timer::Now() - t0).count());

    // Assemble the linear system.
    if (step > step0)
    {
      // Update frequency-dependent excitation and operators.
      A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
      A = spaceop.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                  std::complex<double>(-omega * omega, 0.0), K.get(),
                                  C.get(), M.get(), A2.get());
      P = spaceop.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega,
                                                           omega);
      ksp.SetOperators(*A, *P);
    }
    spaceop.GetExcitationVector(omega, RHS);

    // Solve the linear system.
    BlockTimer bt1(Timer::SOLVE);
    Mpi::Print("\n");
    ksp.Mult(RHS, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    BlockTimer bt2(Timer::POSTPRO);
    double E_elec = 0.0, E_mag = 0.0;
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(), omega);
    Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(spaceop.GetComm(), E),
               linalg::Norml2(spaceop.GetComm(), RHS));
    if (!iodata.solver.driven.only_port_post)
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      E_elec = postop.GetEFieldEnergy();
      E_mag = postop.GetHFieldEnergy();
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
                 E_mag * J, (E_elec + E_mag) * J);
    }

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, indicator);

    // Postprocess S-parameters and optionally write solution to disk.
    Postprocess(postop, spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(),
                spaceop.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                !iodata.solver.driven.only_port_post,
                (step == nstep - 1) ? &indicator : nullptr);

    // Increment frequency.
    step++;
    omega += delta_omega;
  }
  SaveMetadata(ksp);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &spaceop, PostOperator &postop,
                                           int nstep, int step0, double omega0,
                                           double delta_omega) const
{
  // Configure default parameters if not specified.
  BlockTimer bt0(Timer::CONSTRUCT);
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int max_size = iodata.solver.driven.adaptive_max_size;
  MFEM_VERIFY(max_size <= 0 || max_size > 2,
              "Adaptive frequency sweep must sample at least two frequency points!");
  if (max_size <= 0)
  {
    max_size = 20;  // Default value
  }
  max_size = std::min(max_size, nstep - step0);  // Maximum size dictated by sweep
  int convergence_memory = iodata.solver.driven.adaptive_memory;

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = spaceop.GetCurlMatrix();
  ComplexVector E(Curl.Width()), E2(Curl.Width()), B(Curl.Height());
  E = 0.0;
  E2 = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  CurlFluxErrorEstimator<ComplexVector> estimator(
      spaceop.GetMaterialOp(), spaceop.GetNDSpace(), iodata.solver.linear.estimator_tol,
      iodata.solver.linear.estimator_max_it, 0);
  ErrorIndicator indicator;

  // Configure the PROM operator which performs the parameter space sampling and basis
  // construction during the offline phase as well as the PROM solution during the online
  // phase.
  auto t0 = Timer::Now();
  const double f0 = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, 1.0);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             nstep - step0, omega0 * f0, (omega0 + (nstep - step0 - 1) * delta_omega) * f0);
  RomOperator promop(iodata, spaceop, max_size);
  spaceop.GetWavePortOp().SetSuppressOutput(true);  // Suppress wave port output for offline

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  BlockTimer bt1(Timer::CONSTRUCTPROM);
  promop.SolveHDM(omega0, E);
  promop.UpdatePROM(omega0, E);
  estimator.AddErrorIndicator(E, indicator);
  promop.SolveHDM(omega0 + (nstep - step0 - 1) * delta_omega, E);
  promop.UpdatePROM(omega0 + (nstep - step0 - 1) * delta_omega, E);
  estimator.AddErrorIndicator(E, indicator);

  // Greedy procedure for basis construction (offline phase). Basis is initialized with
  // solutions at frequency sweep endpoints.
  int it = 2, it0 = it, memory = 0;
  double max_error;
  while (true)
  {
    // Compute the location of the maximum error in parameter domain.
    double omega_star = promop.FindMaxError(omega0, delta_omega, nstep - step0);

    // Compute the actual solution error at the given parameter point.
    promop.SolveHDM(omega_star, E);
    promop.SolvePROM(omega_star, E2);
    linalg::AXPY(-1.0, E, E2);
    max_error =
        linalg::Norml2(spaceop.GetComm(), E2) / linalg::Norml2(spaceop.GetComm(), E);
    if (max_error < offline_tol)
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
    Mpi::Print(
        "\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
        "{:.3e}{}\n",
        it - it0 + 1, promop.GetReducedDimension(), omega_star * f0, omega_star, max_error,
        (memory == 0) ? ""
                      : fmt::format(", memory = {:d}/{:d}", memory, convergence_memory));
    promop.UpdatePROM(omega_star, E);
    estimator.AddErrorIndicator(E, indicator);
    it++;
  }
  Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
             " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
             (it == max_size) ? " reached maximum" : " converged with", it,
             promop.GetReducedDimension(), max_error, offline_tol, memory,
             convergence_memory);
  utils::PrettyPrint(promop.GetSamplePoints(), f0, " Sampled frequencies (GHz):");
  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root
  SaveMetadata(promop.GetLinearSolver());

  // XX TODO: Add output of eigenvalue estimates from the PROM system (and nonlinear EVP in
  //          the general case with wave ports, etc.?)
  const auto eigs = promop.ComputeEigenvalueEstimates(omega0);
  if (Mpi::Root(spaceop.GetComm()))
  {
    std::cout << "Eigenvalues (nev = " << eigs.size() << "):\n";
    for (auto omega : eigs)
    {
      const std::complex<double> f = {
          iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.real()),
          iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.imag())};
      const double Q =
          (f.imag() == 0.0) ? mfem::infinity() : 0.5 * std::abs(f) / std::abs(f.imag());
      std::cout << f << ", " << Q << "\n";
    }
  }
  return indicator;

  // Main fast frequency sweep loop (online phase).
  BlockTimer bt2(Timer::CONSTRUCT);
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  spaceop.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  int step = step0;
  double omega = omega0;
  while (step < nstep)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
               nstep, freq, Timer::Duration(Timer::Now() - t0).count());

    // Assemble and solve the PROM linear system.
    BlockTimer bt3(Timer::SOLVEPROM);
    promop.SolvePROM(omega, E);
    Mpi::Print("\n");  // Newline after PROM assembly

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    BlockTimer bt4(Timer::POSTPRO);
    double E_elec = 0.0, E_mag = 0.0;
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(), omega);
    Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(spaceop.GetComm(), E));
    if (!iodata.solver.driven.only_port_post)
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      E_elec = postop.GetEFieldEnergy();
      E_mag = postop.GetHFieldEnergy();
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", E_elec * J,
                 E_mag * J, (E_elec + E_mag) * J);
    }

    // Postprocess S-parameters and optionally write solution to disk.
    Postprocess(postop, spaceop.GetLumpedPortOp(), spaceop.GetWavePortOp(),
                spaceop.GetSurfaceCurrentOp(), step, omega, E_elec, E_mag,
                !iodata.solver.driven.only_port_post,
                (step == nstep - 1) ? &indicator : nullptr);

    // Increment frequency.
    step++;
    omega += delta_omega;
  }
  return indicator;
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

void DrivenSolver::Postprocess(const PostOperator &postop,
                               const LumpedPortOperator &lumped_port_op,
                               const WavePortOperator &wave_port_op,
                               const SurfaceCurrentOperator &surf_j_op, int step,
                               double omega, double E_elec, double E_mag, bool full,
                               const ErrorIndicator *indicator) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main frequency sweep loop.
  double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
  PostprocessCurrents(postop, surf_j_op, step, omega);
  PostprocessPorts(postop, lumped_port_op, step, omega);
  if (surf_j_op.Size() == 0)
  {
    PostprocessSParameters(postop, lumped_port_op, wave_port_op, step, omega);
  }
  if (full)
  {
    double E_cap = postop.GetLumpedCapacitorEnergy(lumped_port_op);
    double E_ind = postop.GetLumpedInductorEnergy(lumped_port_op);
    PostprocessDomains(postop, "f (GHz)", step, freq, E_elec, E_mag, E_cap, E_ind);
    PostprocessSurfaces(postop, "f (GHz)", step, freq, E_elec + E_cap, E_mag + E_ind, 1.0,
                        1.0);
    PostprocessProbes(postop, "f (GHz)", step, freq);
  }
  if (iodata.solver.driven.delta_post > 0 && step % iodata.solver.driven.delta_post == 0)
  {
    Mpi::Print("\n");
    PostprocessFields(postop, step / iodata.solver.driven.delta_post, freq, indicator);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(postop, *indicator);
  }
}

namespace
{

struct CurrentData
{
  const int idx;      // Current source index
  const double Iinc;  // Excitation current
};

struct PortVIData
{
  const int idx;                      // Lumped port index
  const bool excitation;              // Flag for excited ports
  const double Vinc, Iinc;            // Incident voltage, current
  const std::complex<double> Vi, Ii;  // Port voltage, current
};

struct PortSData
{
  const int idx;                   // Port index
  const std::complex<double> Sij;  // Scattering parameter
};

}  // namespace

void DrivenSolver::PostprocessCurrents(const PostOperator &postop,
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
    j_data.push_back({idx, iodata.DimensionalizeValue(IoData::ValueType::CURRENT, Iinc)});
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
                     "Iinc[" + std::to_string(data.idx) + "] (A)", table.w,
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

void DrivenSolver::PostprocessPorts(const PostOperator &postop,
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
    port_data.push_back({idx, data.IsExcited(),
                         iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, Vinc),
                         iodata.DimensionalizeValue(IoData::ValueType::CURRENT, Iinc),
                         iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, Vi),
                         iodata.DimensionalizeValue(IoData::ValueType::CURRENT, Ii)});
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

void DrivenSolver::PostprocessSParameters(const PostOperator &postop,
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
