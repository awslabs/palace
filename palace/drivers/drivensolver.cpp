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
  SaveMetadata(excitation_helper);

  // Frequencies will be sampled uniformly in the frequency domain.
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

  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             excitation_helper.FmtLog());

  // Main frequency sweep loop.
  return {
      adaptive
          ? SweepAdaptive(space_op, excitation_helper, n_step, step0, omega0, delta_omega)
          : SweepUniform(space_op, excitation_helper, n_step, step0, omega0, delta_omega),
      space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op,
                                          const PortExcitationHelper &excitation_helper,
                                          int n_step, int step0, double omega0,
                                          double delta_omega) const
{
  // Initialize postprocessing operator and printers.
  // Initialize write directory with default path; will be changed if multiple excitations.
  PostOperator post_op(iodata, space_op, "driven");
  PostprocessPrintResults post_results(root, post_dir, post_op, space_op, excitation_helper,
                                       n_step, iodata.solver.driven.delta_post);

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
  ComplexVector RHS(Curl.Width());
  ComplexVector E(Curl.Width());
  ComplexVector B(Curl.Height());
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

  auto t0 = Timer::Now();
  auto print_counter_excitation = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    print_counter_excitation++;
    Mpi::Print("\nSweeping Excitation Index {:d} ({:d}/{:d}):\n", excitation_idx,
               print_counter_excitation, excitation_helper.Size());

    if (excitation_helper.Size() > 1)
    {
      // Switch to multiple paraview subfolders: one for each excitation
      post_op.SetNewParaviewOutput(
          ParaviewPath(iodata, excitation_idx, excitation_helper.Size()));
    }

    // Frequency loop
    double omega = omega0;
    for (int step = step0; step < n_step; step++, omega += delta_omega)
    {
      const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
      Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
                 step + 1, n_step, freq, Timer::Duration(Timer::Now() - t0).count());

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

      post_results.PostprocessStep(iodata, post_op, space_op, step, excitation_idx);
    }
    // Final postprocessing & printing
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(ksp);
  }
  post_results.PostprocessFinal(post_op, indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op,
                                           const PortExcitationHelper &excitation_helper,
                                           int n_step, int step0, double omega0,
                                           double delta_omega) const
{
  // Initialize postprocessing operator and printers.
  // Initialize write directory with default path; will be changed if multiple excitations.
  PostOperator post_op(iodata, space_op, "driven");
  PostprocessPrintResults post_results(root, post_dir, post_op, space_op, excitation_helper,
                                       n_step, iodata.solver.driven.delta_post);

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  int max_size = iodata.solver.driven.adaptive_max_size;
  MFEM_VERIFY(max_size <= 0 || max_size > 2,
              "Adaptive frequency sweep must sample at least two frequency points!");
  if (max_size <= 0)
  {
    max_size = 20 * excitation_helper.Size();  // Default value
  }
  max_size = std::min(max_size,
                      (n_step - step0) *
                          int(excitation_helper.Size()));  // Maximum size dictated by sweep
  int convergence_memory = iodata.solver.driven.adaptive_memory;

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector E(Curl.Width());
  ComplexVector Eh(Curl.Width());
  ComplexVector B(Curl.Height());
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
  const double unit_GHz = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, 1.0);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             n_step - step0, omega0 * unit_GHz,
             (omega0 + (n_step - step0 - 1) * delta_omega) * unit_GHz);
  RomOperator prom_op(iodata, space_op, max_size);
  space_op.GetWavePortOp().SetSuppressOutput(
      true);  // Suppress wave port output for offline

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  int paraview_step = 0;
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
    post_op.SetEGridFunction(E, false);
    post_op.SetBGridFunction(B, false);
    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();
    estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);
  };

  // Loop excitations to add to PROM
  auto print_counter_excitation_prom = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    print_counter_excitation_prom++;
    Mpi::Print("\nAdding Excitation Index {:d} ({:d}/{:d}):\n", excitation_idx,
               print_counter_excitation_prom, excitation_helper.Size());

    prom_op.SetExcitationIndex(excitation_idx);  // Pre-compute RHS1

    prom_op.SolveHDM(excitation_idx, omega0, E);
    UpdatePROM(excitation_idx, omega0);
    prom_op.SolveHDM(excitation_idx, omega0 + (n_step - step0 - 1) * delta_omega, E);
    UpdatePROM(excitation_idx, omega0 + (n_step - step0 - 1) * delta_omega);

    // Greedy procedure for basis construction (offline phase). Basis is initialized with
    // solutions at frequency sweep endpoints.
    int it = 2, it0 = it, memory = 0;
    std::vector<double> max_errors = {0.0, 0.0};
    while (true)
    {
      // Compute the location of the maximum error in parameter domain (bounded by the
      // previous samples).
      double omega_star = prom_op.FindMaxError(excitation_idx)[0];

      // Compute the actual solution error at the given parameter point.
      prom_op.SolveHDM(excitation_idx, omega_star, E);
      prom_op.SolvePROM(excitation_idx, omega_star, Eh);
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
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(),
                 (memory == 0)
                     ? ""
                     : fmt::format(", memory = {:d}/{:d}", memory, convergence_memory));
      UpdatePROM(excitation_idx, omega_star);
      it++;
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
  }

  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root

  // XX TODO: Add output of eigenvalue estimates from the PROM system (and nonlinear EVP in
  //          the general case with wave ports, etc.?)

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  auto print_counter_excitation_online = 0;  // 1 based indexing; will increment at start
  for (const auto &[excitation_idx, spec] : excitation_helper.excitations)
  {
    print_counter_excitation_online++;
    Mpi::Print("\nSweeping Excitation Index {:d} ({:d}/{:d}):\n", excitation_idx,
               print_counter_excitation_online, excitation_helper.Size());

    if (excitation_helper.Size() > 1)
    {
      // Switch to multiple paraview subfolders: one for each excitation
      post_op.SetNewParaviewOutput(
          ParaviewPath(iodata, excitation_idx, excitation_helper.Size()));
    }

    // Frequency loop
    double omega = omega0;
    for (int step = step0; step < n_step; step++, omega += delta_omega)
    {
      const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
      Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (elapsed time = {:.2e} s)\n", step + 1,
                 n_step, freq, Timer::Duration(Timer::Now() - t0).count());

      // Assemble and solve the PROM linear system.
      prom_op.SolvePROM(excitation_idx, omega, E);
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

      post_results.PostprocessStep(iodata, post_op, space_op, step, excitation_idx);
    }
    // Final postprocessing & printing: no change to indicator since prom
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(prom_op.GetLinearSolver());
  }
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

// Frequency checking when adding to a new excitation block: When adding data to data_col
// with freqeucny freq, checks that freq matches what is already written in the
// corresponding row of freq_col. Adds a new frequency row to freq_col if needed.
void set_validate_freq_col_alignment(double freq, Column &freq_col, const Column &data_col)
{
  auto data_curser = data_col.n_rows();
  auto freq_curser = freq_col.n_rows();
  if (data_curser == freq_curser)
  {
    freq_col << freq;
  }
  else
  {
    double current_freq = freq_col.data.at(data_curser);
    MFEM_VERIFY(
        freq == current_freq,
        fmt::format("Writing data table at incorrect frequency. Data has frequency {} "
                    "while table is at {}",
                    freq, current_freq));
  }
}

DrivenSolver::DomainsPostPrinter::DomainsPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir, const PostOperator &post_op,
    const PortExcitationHelper &excitation_helper, int n_expected_rows)
  : do_measurement_{do_measurement}, root_{root}
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;

  domain_E = TableWithCSVFile(post_dir / "domain-E.csv");
  domain_E.table.reserve(n_expected_rows, 4 + post_op.GetDomainPostOp().M_i.size());
  domain_E.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  bool single = (excitation_helper.Size() == 1);
  for (const auto &[ex_idx, data] : excitation_helper.excitations)
  {
    std::string multi_ex = single ? "" : format("[{}]", ex_idx);

    domain_E.table.insert_column(format("Ee_{}", ex_idx), format("E_elec{} (J)", multi_ex));
    domain_E.table.insert_column(format("Em_{}", ex_idx), format("E_mag{} (J)", multi_ex));
    domain_E.table.insert_column(format("Ec_{}", ex_idx), format("E_cap{} (J)", multi_ex));
    domain_E.table.insert_column(format("Ei_{}", ex_idx), format("E_ind{} (J)", multi_ex));

    for (const auto &[idx, data] : post_op.GetDomainPostOp().M_i)
    {
      domain_E.table.insert_column(format("Ee_{}_{}", idx, ex_idx),
                                   format("E_elec[{}]{} (J)", idx, multi_ex));
      domain_E.table.insert_column(format("pe_{}_{}", idx, ex_idx),
                                   format("p_elec[{}]{}", idx, multi_ex));
      domain_E.table.insert_column(format("Em_{}_{}", idx, ex_idx),
                                   format("E_mag[{}]{} (J)", idx, multi_ex));
      domain_E.table.insert_column(format("pm_{}_{}", idx, ex_idx),
                                   format("p_mag[{}]{}", idx, multi_ex));
    }
  }
  domain_E.WriteFullTableTrunc();
}

void DrivenSolver::DomainsPostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                      int excitation_idx,
                                                      const PostOperator &post_op,
                                                      const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  double oneJ = iodata.DimensionalizeValue(VT::ENERGY, 1.0);

  set_validate_freq_col_alignment(idx_value_dimensionful, domain_E.table["idx"],
                                  domain_E.table[format("Ee_{}", excitation_idx)]);

  double E_elec = post_op.GetEFieldEnergy();
  double E_mag = post_op.GetHFieldEnergy();

  domain_E.table[format("Ee_{}", excitation_idx)] << E_elec * oneJ;
  domain_E.table[format("Em_{}", excitation_idx)] << E_mag * oneJ;
  domain_E.table[format("Ec_{}", excitation_idx)]
      << post_op.GetLumpedCapacitorEnergy() * oneJ;
  domain_E.table[format("Ei_{}", excitation_idx)]
      << post_op.GetLumpedInductorEnergy() * oneJ;

  // Write the field and lumped element energies.
  for (const auto &[idx, data] : post_op.GetDomainPostOp().M_i)
  {
    double E_e = post_op.GetEFieldEnergy(idx);
    double E_m = post_op.GetHFieldEnergy(idx);
    domain_E.table[format("Ee_{}_{}", idx, excitation_idx)] << E_e * oneJ;
    domain_E.table[format("pe_{}_{}", idx, excitation_idx)]
        << ((std::abs(E_elec) > 0.0) ? (E_e / E_elec) : 0.0);
    domain_E.table[format("Em_{}_{}", idx, excitation_idx)] << E_m * oneJ;
    domain_E.table[format("pm_{}_{}", idx, excitation_idx)]
        << ((std::abs(E_mag) > 0.0) ? E_m / E_mag : 0.0);
  }

  domain_E.WriteFullTableTrunc();
}

DrivenSolver::SurfacesPostPrinter::SurfacesPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir, const PostOperator &post_op,
    const PortExcitationHelper &excitation_helper, int n_expected_rows)
  : root_{root},
    do_measurement_flux_(do_measurement                                       //
                         && post_op.GetSurfacePostOp().flux_surfs.size() > 0  // Has flux
                         ),
    do_measurement_eps_(do_measurement                                      //
                        && post_op.GetSurfacePostOp().eps_surfs.size() > 0  // Has eps
    )
{
  if (!root_)
  {
    return;
  }
  using fmt::format;

  bool has_imaginary = post_op.HasImag();
  bool single = (excitation_helper.Size() == 1);

  if (do_measurement_flux_)
  {
    surface_F = TableWithCSVFile(post_dir / "surface-F.csv");
    surface_F.table.reserve(n_expected_rows,
                            2 * post_op.GetSurfacePostOp().flux_surfs.size() + 1);
    surface_F.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

    for (const auto &[ex_idx, data] : excitation_helper.excitations)
    {
      std::string multi_ex = single ? "" : format("[{}]", ex_idx);
      for (const auto &[idx, data] : post_op.GetSurfacePostOp().flux_surfs)
      {
        switch (data.type)
        {
          case SurfaceFluxType::ELECTRIC:
            if (has_imaginary)
            {
              surface_F.table.insert_column(
                  format("F_{}_{}_re", idx, ex_idx),
                  format("Re{{Φ_elec[{}]{}}} (C)", idx, multi_ex));
              surface_F.table.insert_column(
                  format("F_{}_{}_im", idx, ex_idx),
                  format("Im{{Φ_elec[{}]{}}} (C)", idx, multi_ex));
            }
            else
            {
              surface_F.table.insert_column(format("F_{}_{}_re", idx, ex_idx),
                                            format("Φ_elec[{}]{} (C)", idx, multi_ex));
            }
            break;
          case SurfaceFluxType::MAGNETIC:
            if (has_imaginary)
            {
              surface_F.table.insert_column(
                  format("F_{}_{}_re", idx, ex_idx),
                  format("Re{{Φ_mag[{}]{}}} (Wb)", idx, multi_ex));
              surface_F.table.insert_column(
                  format("F_{}_{}_im", idx, ex_idx),
                  format("Im{{Φ_mag[{}]{}}} (Wb)", idx, multi_ex));
            }
            else
            {
              surface_F.table.insert_column(format("F_{}_{}_re", idx, ex_idx),
                                            format("Φ_mag[{}]{} (Wb)", idx, multi_ex));
            }
            break;
          case SurfaceFluxType::POWER:
            surface_F.table.insert_column(format("F_{}_{}_re", idx, ex_idx),
                                          format("Φ_pow[{}]{} (W)", idx, multi_ex));
            break;
        }
      }
      surface_F.WriteFullTableTrunc();
    }
  }
  if (do_measurement_eps_)
  {
    surface_Q = TableWithCSVFile(post_dir / "surface-Q.csv");
    surface_Q.table.reserve(n_expected_rows,
                            2 * post_op.GetSurfacePostOp().eps_surfs.size() + 1);
    surface_Q.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

    for (const auto &[ex_idx, data] : excitation_helper.excitations)
    {
      std::string multi_ex = single ? "" : format("[{}]", ex_idx);
      for (const auto &[idx, data] : post_op.GetSurfacePostOp().eps_surfs)
      {
        surface_Q.table.insert_column(format("p_{}_{}", idx, ex_idx),
                                      format("p_surf[{}]{}", idx, multi_ex));
        surface_Q.table.insert_column(format("Q_{}_{}", idx, ex_idx),
                                      format("Q_surf[{}]{}", idx, multi_ex));
      }
    }
  }
}

void DrivenSolver::SurfacesPostPrinter::AddMeasurementFlux(double idx_value_dimensionful,
                                                           int excitation_idx,
                                                           const PostOperator &post_op,
                                                           const IoData &iodata)
{
  if (!do_measurement_flux_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  const bool has_imaginary = post_op.HasImag();
  auto flux_data_vec = post_op.GetSurfaceFluxAll();
  auto dimensionlize_flux = [&iodata](auto Phi, SurfaceFluxType flux_type)
  {
    switch (flux_type)
    {
      case SurfaceFluxType::ELECTRIC:
        Phi *= iodata.DimensionalizeValue(VT::CAPACITANCE, 1.0);
        Phi *= iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
        break;
      case SurfaceFluxType::MAGNETIC:
        Phi *= iodata.DimensionalizeValue(VT::INDUCTANCE, 1.0);
        Phi *= iodata.DimensionalizeValue(VT::CURRENT, 1.0);
        break;
      case SurfaceFluxType::POWER:
        Phi *= iodata.DimensionalizeValue(VT::POWER, 1.0);
        break;
    }
    return Phi;
  };

  set_validate_freq_col_alignment(
      idx_value_dimensionful, surface_F.table["idx"],
      surface_F.table[format("F_{}_{}_re", flux_data_vec.begin()->idx, excitation_idx)]);

  for (const auto &flux_data : flux_data_vec)
  {
    auto Phi_unitful = dimensionlize_flux(flux_data.Phi, flux_data.type);
    surface_F.table[format("F_{}_{}_re", flux_data.idx, excitation_idx)]
        << Phi_unitful.real();
    if (has_imaginary && (flux_data.type == SurfaceFluxType::ELECTRIC ||
                          flux_data.type == SurfaceFluxType::MAGNETIC))
    {
      surface_F.table[format("F_{}_{}_im", flux_data.idx, excitation_idx)]
          << Phi_unitful.imag();
    }
  }
  surface_F.WriteFullTableTrunc();
}

void DrivenSolver::SurfacesPostPrinter::AddMeasurementEps(double idx_value_dimensionful,
                                                          int excitation_idx,
                                                          const PostOperator &post_op,
                                                          const IoData &iodata)
{
  if (!do_measurement_eps_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Interface Participation adds energy contriutions E_elec + E_cap
  // E_cap returns zero if the solver does not supprot lumped ports.
  double E_elec = post_op.GetEFieldEnergy() + post_op.GetLumpedCapacitorEnergy();
  auto eps_data_vec = post_op.GetInterfaceEFieldEnergyAll();

  set_validate_freq_col_alignment(
      idx_value_dimensionful, surface_Q.table["idx"],
      surface_Q.table[format("Q_{}_{}", eps_data_vec.begin()->idx, excitation_idx)]);

  for (const auto &eps_data : eps_data_vec)
  {
    double p = post_op.GetInterfaceParticipation(eps_data.idx, E_elec);
    double tandelta = eps_data.tandelta;
    double Q = (p == 0.0 || tandelta == 0.0) ? mfem::infinity() : 1.0 / (tandelta * p);
    surface_Q.table[format("p_{}_{}", eps_data.idx, excitation_idx)] << p;
    surface_Q.table[format("Q_{}_{}", eps_data.idx, excitation_idx)] << Q;
  }
  surface_Q.WriteFullTableTrunc();
}

void DrivenSolver::SurfacesPostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                       int excitation_idx,
                                                       const PostOperator &post_op,
                                                       const IoData &iodata)
{
  // If surfaces have been specified for postprocessing, compute the corresponding values
  // and write out to disk. The passed in E_elec is the sum of the E-field and lumped
  // capacitor energies, and E_mag is the same for the B-field and lumped inductors.
  AddMeasurementFlux(idx_value_dimensionful, excitation_idx, post_op, iodata);
  AddMeasurementEps(idx_value_dimensionful, excitation_idx, post_op, iodata);
}

DrivenSolver::ProbePostPrinter::ProbePostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir, const PostOperator &post_op,
    const PortExcitationHelper &excitation_helper, int n_expected_rows)
  : root_{root}, do_measurement_E_{do_measurement}, do_measurement_B_{do_measurement},
    has_imag{post_op.HasImag()}, v_dim{post_op.GetInterpolationOpVDim()}
{
#if defined(MFEM_USE_GSLIB)
  do_measurement_E_ = do_measurement_E_                    //
                      && (post_op.GetProbes().size() > 0)  // Has probes defined
                      && post_op.HasE();                   // Has E fields

  do_measurement_B_ = do_measurement_B_                    //
                      && (post_op.GetProbes().size() > 0)  // Has probes defined
                      && post_op.HasB();                   // Has B fields

  if (!root_ || (!do_measurement_E_ && !do_measurement_B_))
  {
    return;
  }
  using fmt::format;
  int scale_col = (has_imag ? 2 : 1) * v_dim;
  auto dim_labeler = [](int i) -> std::string
  {
    switch (i)
    {
      // Note: Zero-based indexing here
      case 0:
        return "x";
      case 1:
        return "y";
      case 2:
        return "z";
      default:
        return format("d{}", i);
    }
  };

  bool single = (excitation_helper.Size() == 1);
  if (do_measurement_E_)
  {
    probe_E = TableWithCSVFile(post_dir / "probe-E.csv");
    probe_E.table.reserve(n_expected_rows, scale_col * post_op.GetProbes().size());
    probe_E.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

    for (const auto &[ex_idx, data] : excitation_helper.excitations)
    {
      std::string multi_ex = single ? "" : format("[{}]", ex_idx);
      for (const auto &idx : post_op.GetProbes())
      {
        for (int i_dim = 0; i_dim < v_dim; i_dim++)
        {
          if (has_imag)
          {
            probe_E.table.insert_column(
                format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                format("Re{{E_{}[{}]{}}} (V/m)", dim_labeler(i_dim), idx, multi_ex));
            probe_E.table.insert_column(
                format("E{}_{}_{}_im", i_dim, idx, ex_idx),
                format("Im{{E_{}[{}]{}}} (V/m)", dim_labeler(i_dim), idx, multi_ex));
          }
          else
          {
            probe_E.table.insert_column(
                format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                format("E_{}[{}]{} (V/m)", dim_labeler(i_dim), idx, multi_ex));
          }
        }
      }
      probe_E.WriteFullTableTrunc();
    }
  }

  if (do_measurement_B_)
  {
    probe_B = TableWithCSVFile(post_dir / "probe-B.csv");
    probe_B.table.reserve(n_expected_rows, scale_col * post_op.GetProbes().size());
    probe_B.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

    for (const auto &[ex_idx, data] : excitation_helper.excitations)
    {
      std::string multi_ex = single ? "" : format("[{}]", ex_idx);

      for (const auto &idx : post_op.GetProbes())
      {
        for (int i_dim = 0; i_dim < v_dim; i_dim++)
        {
          if (has_imag)
          {
            probe_B.table.insert_column(
                format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                format("Re{{B_{}[{}]{}}} (Wb/m²)", dim_labeler(i_dim), idx, multi_ex));
            probe_B.table.insert_column(
                format("B{}_{}_{}_im", i_dim, idx, ex_idx),
                format("Im{{B_{}[{}]{}}} (Wb/m²)", dim_labeler(i_dim), idx, multi_ex));
          }
          else
          {
            probe_B.table.insert_column(
                format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                format("B_{}[{}]{} (Wb/m²)", dim_labeler(i_dim), idx, multi_ex));
          }
        }
      }
    }
    probe_B.WriteFullTableTrunc();
  }
#endif
}

void DrivenSolver::ProbePostPrinter::AddMeasurementE(double idx_value_dimensionful,
                                                     int excitation_idx,
                                                     const PostOperator &post_op,
                                                     const IoData &iodata)
{
  if (!do_measurement_E_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  auto probe_field = post_op.ProbeEField();
  MFEM_VERIFY(probe_field.size() == v_dim * post_op.GetProbes().size(),
              format("Size mismatch: expect vector field to ahve size {} * {} = {}; got {}",
                     v_dim, post_op.GetProbes().size(), v_dim * post_op.GetProbes().size(),
                     probe_field.size()))

  auto ref_col = format("E{}_{}_{}_re", 0, *post_op.GetProbes().begin(), excitation_idx);
  set_validate_freq_col_alignment(idx_value_dimensionful, probe_E.table["idx"],
                                  probe_E.table[ref_col]);

  size_t i = 0;
  for (const auto &idx : post_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = iodata.DimensionalizeValue(VT::FIELD_E, probe_field[i * v_dim + i_dim]);
      probe_E.table[format("E{}_{}_{}_re", i_dim, idx, excitation_idx)] << val.real();
      if (has_imag)
      {
        probe_E.table[format("E{}_{}_{}_im", i_dim, idx, excitation_idx)] << val.imag();
      }
    }
    i++;
  }
  probe_E.WriteFullTableTrunc();
}

void DrivenSolver::ProbePostPrinter::AddMeasurementB(double idx_value_dimensionful,
                                                     int excitation_idx,
                                                     const PostOperator &post_op,
                                                     const IoData &iodata)
{
  if (!do_measurement_B_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  auto probe_field = post_op.ProbeBField();
  MFEM_VERIFY(probe_field.size() == v_dim * post_op.GetProbes().size(),
              format("Size mismatch: expect vector field to ahve size {} * {} = {}; got {}",
                     v_dim, post_op.GetProbes().size(), v_dim * post_op.GetProbes().size(),
                     probe_field.size()))

  auto ref_col = format("B{}_{}_{}_re", 0, *post_op.GetProbes().begin(), excitation_idx);
  set_validate_freq_col_alignment(idx_value_dimensionful, probe_B.table["idx"],
                                  probe_B.table[ref_col]);

  size_t i = 0;
  for (const auto &idx : post_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = iodata.DimensionalizeValue(VT::FIELD_B, probe_field[i * v_dim + i_dim]);
      probe_B.table[format("B{}_{}_{}_re", i_dim, idx, excitation_idx)] << val.real();
      if (has_imag)
      {
        probe_B.table[format("B{}_{}_{}_im", i_dim, idx, excitation_idx)] << val.imag();
      }
    }
    i++;
  }
  probe_B.WriteFullTableTrunc();
}

void DrivenSolver::ProbePostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                    int excitation_idx,
                                                    const PostOperator &post_op,
                                                    const IoData &iodata)
{
#if defined(MFEM_USE_GSLIB)
  AddMeasurementE(idx_value_dimensionful, excitation_idx, post_op, iodata);
  AddMeasurementB(idx_value_dimensionful, excitation_idx, post_op, iodata);
#endif
}

DrivenSolver::CurrentsPostPrinter::CurrentsPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir,
    const SurfaceCurrentOperator &surf_j_op, const PortExcitationHelper &excitation_helper,
    int n_expected_rows)
  : root_{root}, do_measurement_{
                     do_measurement             //
                     && (surf_j_op.Size() > 0)  // Needs surface currents
                 }
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;

  surface_I = TableWithCSVFile(post_dir / "surface-I.csv");
  surface_I.table.reserve(n_expected_rows, surf_j_op.Size() + 1);
  surface_I.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  bool single = (excitation_helper.Size() == 1);
  for (const auto &[ex_idx, data] : excitation_helper.excitations)
  {
    std::string multi_ex = single ? "" : format("[{}]", ex_idx);

    for (const auto &[idx, data] : surf_j_op)
    {
      surface_I.table.insert_column(format("I_{}_{}", idx, ex_idx),
                                    format("I_inc[{}]{} (A)", idx, multi_ex));
    }
  }
  surface_I.WriteFullTableTrunc();
}

void DrivenSolver::CurrentsPostPrinter::AddMeasurement(
    double freq, int excitation_idx, const SurfaceCurrentOperator &surf_j_op,
    const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  set_validate_freq_col_alignment(
      freq, surface_I.table["idx"],
      surface_I.table[format("I_{}_{}", surf_j_op.begin()->first, excitation_idx)]);

  for (const auto &[idx, data] : surf_j_op)
  {
    auto I_inc = data.GetExcitationCurrent();
    surface_I.table[format("I_{}_{}", idx, excitation_idx)]
        << iodata.DimensionalizeValue(VT::CURRENT, I_inc);
  }
  surface_I.WriteFullTableTrunc();
}

DrivenSolver::PortsPostPrinter::PortsPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir,
    const LumpedPortOperator &lumped_port_op, const PortExcitationHelper &excitation_helper,
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

  bool single = (excitation_helper.Size() == 1);
  for (const auto &[ex_idx, data] : excitation_helper.excitations)
  {
    std::string multi_ex = single ? "" : format("[{}]", ex_idx);
    for (const auto &idx : data.lumped_port)
    {
      port_V.table.insert_column(format("inc{}_{}", idx, ex_idx),
                                 format("V_inc[{}]{} (V)", idx, multi_ex));
      port_I.table.insert_column(format("inc{}_{}", idx, ex_idx),
                                 format("I_inc[{}]{} (A)", idx, multi_ex));
    }
    for (const auto &[idx, data] : lumped_port_op)
    {
      port_V.table.insert_column(format("re{}_{}", idx, ex_idx),
                                 format("Re{{V[{}]{}}} (V)", idx, multi_ex));
      port_V.table.insert_column(format("im{}_{}", idx, ex_idx),
                                 format("Im{{V[{}]{}}} (V)", idx, multi_ex));
      port_I.table.insert_column(format("re{}_{}", idx, ex_idx),
                                 format("Re{{I[{}]{}}} (A)", idx, multi_ex));
      port_I.table.insert_column(format("im{}_{}", idx, ex_idx),
                                 format("Im{{I[{}]{}}} (A)", idx, multi_ex));
    }
    port_V.WriteFullTableTrunc();
    port_I.WriteFullTableTrunc();
  }
}

void DrivenSolver::PortsPostPrinter::AddMeasurement(
    double freq, int excitation_idx, const PostOperator &post_op,
    const LumpedPortOperator &lumped_port_op, const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).
  auto ref_col_str = format("re{}_{}", lumped_port_op.begin()->first, excitation_idx);
  set_validate_freq_col_alignment(freq, port_V.table["idx"], port_V.table[ref_col_str]);
  set_validate_freq_col_alignment(freq, port_I.table["idx"], port_I.table[ref_col_str]);

  auto unit_V = iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
  auto unit_A = iodata.DimensionalizeValue(VT::CURRENT, 1.0);

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (data.excitation == excitation_idx)
    {
      double V_inc = data.GetExcitationVoltage();
      double I_inc = (std::abs(V_inc) > 0.0) ? data.GetExcitationPower() / V_inc : 0.0;

      port_V.table[format("inc{}_{}", idx, excitation_idx)] << V_inc * unit_V;
      port_I.table[format("inc{}_{}", idx, excitation_idx)] << I_inc * unit_A;
    }

    std::complex<double> V_i = post_op.GetPortVoltage(idx);
    std::complex<double> I_i = post_op.GetPortCurrent(idx);

    port_V.table[format("re{}_{}", idx, excitation_idx)] << V_i.real() * unit_V;
    port_V.table[format("im{}_{}", idx, excitation_idx)] << V_i.imag() * unit_V;

    port_I.table[format("re{}_{}", idx, excitation_idx)] << I_i.real() * unit_A;
    port_I.table[format("im{}_{}", idx, excitation_idx)] << I_i.imag() * unit_A;
  }
  port_V.WriteFullTableTrunc();
  port_I.WriteFullTableTrunc();
}

DrivenSolver::SParametersPostPrinter::SParametersPostPrinter(
    bool do_measurement, bool root, const fs::path &post_dir,
    const LumpedPortOperator &lumped_port_op, const WavePortOperator &wave_port_op,
    const PortExcitationHelper &excitation_helper, int n_expected_rows)
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
  using fmt::format;

  port_S = TableWithCSVFile(post_dir / "port-S.csv");
  port_S.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_S.table.insert_column(Column("idx", "f (GHz)", 0, {}, {}, ""));

  for (const auto &[ex_idx, data] : excitation_helper.excitations)
  {
    // Already ensured that one of lumped or wave ports are empty
    for (const auto &[o_idx, data] : lumped_port_op)
    {
      port_S.table.insert_column(format("abs_{}_{}", o_idx, ex_idx),
                                 format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      port_S.table.insert_column(format("arg_{}_{}", o_idx, ex_idx),
                                 format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
    for (const auto &[o_idx, data] : wave_port_op)
    {
      port_S.table.insert_column(format("abs_{}_{}", o_idx, ex_idx),
                                 format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      port_S.table.insert_column(format("arg_{}_{}", o_idx, ex_idx),
                                 format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
    port_S.WriteFullTableTrunc();
  }
}

void DrivenSolver::SParametersPostPrinter::AddMeasurement(
    double freq, int excitation_idx, const PostOperator &post_op,
    const LumpedPortOperator &lumped_port_op, const WavePortOperator &wave_port_op,
    const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  std::vector<int> all_port_indices;
  for (const auto &[idx, data] : lumped_port_op)
  {
    all_port_indices.emplace_back(idx);
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    all_port_indices.emplace_back(idx);
  }

  set_validate_freq_col_alignment(
      freq, port_S.table["idx"],
      port_S.table[format("abs_{}_{}", *all_port_indices.begin(), excitation_idx)]);

  for (const auto o_idx : all_port_indices)
  {
    std::complex<double> S_ij =
        post_op.GetSParameter(src_lumped_port, o_idx, excitation_idx);

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

DrivenSolver::PostprocessPrintResults::PostprocessPrintResults(
    bool root, const fs::path &post_dir, const PostOperator &post_op,
    const SpaceOperator &space_op, const PortExcitationHelper &excitation_helper,
    int n_expected_rows, int delta_post_)
  : delta_post{delta_post_}, write_paraview_fields{delta_post_ > 0},
    domains{true, root, post_dir, post_op, excitation_helper, n_expected_rows},
    surfaces{true, root, post_dir, post_op, excitation_helper, n_expected_rows},
    currents{
        true,           root, post_dir, space_op.GetSurfaceCurrentOp(), excitation_helper,
        n_expected_rows},
    probes{true, root, post_dir, post_op, excitation_helper, n_expected_rows},
    ports{true,           root, post_dir, space_op.GetLumpedPortOp(), excitation_helper,
          n_expected_rows},
    s_parameters{true,
                 root,
                 post_dir,
                 space_op.GetLumpedPortOp(),
                 space_op.GetWavePortOp(),
                 excitation_helper,
                 n_expected_rows},
    error_indicator{true, root, post_dir}
{
}

void DrivenSolver::PostprocessPrintResults::PostprocessStep(const IoData &iodata,
                                                            const PostOperator &post_op,
                                                            const SpaceOperator &space_op,
                                                            int step, int excitation_idx)
{
  double omega = post_op.GetFrequency().real();
  auto freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);

  domains.AddMeasurement(freq, excitation_idx, post_op, iodata);
  surfaces.AddMeasurement(freq, excitation_idx, post_op, iodata);
  currents.AddMeasurement(freq, excitation_idx, space_op.GetSurfaceCurrentOp(), iodata);
  probes.AddMeasurement(freq, excitation_idx, post_op, iodata);
  ports.AddMeasurement(freq, excitation_idx, post_op, space_op.GetLumpedPortOp(), iodata);
  s_parameters.AddMeasurement(freq, excitation_idx, post_op, space_op.GetLumpedPortOp(),
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
