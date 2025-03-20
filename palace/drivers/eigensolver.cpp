// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator post_op(iodata, space_op, "eigenmode");
  PostprocessPrintResults post_results(root, post_dir, post_op, space_op,
                                       iodata.solver.eigenmode.n_post);
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  config::EigenSolverData::Type type = iodata.solver.eigenmode.type;
#if defined(PALACE_WITH_ARPACK) && defined(PALACE_WITH_SLEPC)
  if (type == config::EigenSolverData::Type::DEFAULT)
  {
    type = config::EigenSolverData::Type::SLEPC;
  }
#elif defined(PALACE_WITH_ARPACK)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::SLEPC)
  {
    Mpi::Warning("SLEPc eigensolver not available, using ARPACK!\n");
  }
  type = config::EigenSolverData::Type::ARPACK;
#elif defined(PALACE_WITH_SLEPC)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::ARPACK)
  {
    Mpi::Warning("ARPACK eigensolver not available, using SLEPc!\n");
  }
  type = config::EigenSolverData::Type::SLEPC;
#else
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == config::EigenSolverData::Type::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
#endif
  }
  else  // config::EigenSolverData::Type::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (C)
    {
      if (!iodata.solver.eigenmode.pep_linear)
      {
        slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                              iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
    }
    else
    {
      slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    }
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetOrthogonalization(
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::MGS,
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::CGS2);
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (C)
  {
    eigen->SetOperators(*K, *C, *M, scale);
  }
  else
  {
    eigen->SetOperators(*K, *M, scale);
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  eigen->SetTol(iodata.solver.eigenmode.tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  std::unique_ptr<Operator> KM;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  if (iodata.solver.linear.divfree_max_it > 0)
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      space_op.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  const double target = iodata.solver.eigenmode.target;
  {
    const double f_target =
        iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, target);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (C)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen->SetShiftInvert(1i * target);
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. The eigenvalue
      // 1 / (λ - σ) will be a large-magnitude negative imaginary number for an eigenvalue
      // λ with frequency close to but not below the target σ.
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_IMAGINARY);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    eigen->SetShiftInvert(target * target);
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * target,
                                    std::complex<double>(-target * target, 0.0), K.get(),
                                    C.get(), M.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, target, -target * target,
                                                             target);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }
  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Update printer now we know num_conv
  post_results.eigen.stdout_int_print_width = 1 + static_cast<int>(std::log10(num_conv));

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  if (!KM)
  {
    // Normalize the finalized eigenvectors with respect to mass matrix (unit electric field
    // energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_conv);
  }
  Mpi::Print("\n");

  post_results.eigen.PrintStdoutHeader();  // Print headerline for logfile mode
  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (!C)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    eigen->GetEigenvector(i, E);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    post_op.SetEGridFunction(E);
    post_op.SetBGridFunction(B);
    post_op.SetFrequency(omega);
    post_op.MeasureAll();

    const double E_elec = post_op.GetEFieldEnergy();
    const double E_mag = post_op.GetHFieldEnergy();

    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      estimator.AddErrorIndicator(E, B, E_elec + E_mag, indicator);
    }

    // Postprocess state and write fields to file
    post_results.PostprocessStep(iodata, post_op, space_op, i, error_abs, error_bkwd);

    // Final write: Different condition than end of loop (i = num_conv - 1)
    if (i == iodata.solver.eigenmode.n - 1)
    {
      post_results.PostprocessFinal(post_op, indicator);
    }
  }
  return {indicator, space_op.GlobalTrueVSize()};
}

void EigenSolver::EigenPostPrinter::PrintStdoutHeader()
{
  if (!root_)
  {
    return;
  }
  auto save_defaults(eig.table.col_options);

  eig.table.col_options.float_precision = 6;
  eig.table.col_options.min_left_padding = 6;

  // Separate printing due to printing lead as integer not float

  fmt::memory_buffer buf;
  auto to = [&buf](auto f, auto &&...a)
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  to("{}", eig.table[0].format_header(stdout_int_print_width));
  for (int i = 1; i < eig.table.n_cols(); i++)
  {
    if (i > 0)
    {
      to("{:s}", eig.table.print_col_separator);
    }
    to("{:s}", eig.table[i].format_header());
  }
  to("{:s}", eig.table.print_row_separator);

  Mpi::Print("{}{}\n", std::string{buf.data(), buf.size()},
             std::string(stdout_int_print_width + 4 * eig.table[1].col_width(), '='));
  eig.table.col_options = save_defaults;
}

void EigenSolver::EigenPostPrinter::PrintStdoutRow(size_t j)
{
  if (!root_)
  {
    return;
  }
  auto save_defaults(eig.table.col_options);
  eig.table.col_options.float_precision = 6;
  eig.table.col_options.min_left_padding = 6;

  // Separate printing due to integer lead
  fmt::memory_buffer buf;
  auto to = [&buf](auto f, auto &&...a)
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  to("{:{}d}", int(eig.table[0].data[j]), stdout_int_print_width);
  for (int i = 1; i < eig.table.n_cols(); i++)
  {
    if (i > 0)
    {
      to("{:s}", eig.table.print_col_separator);
    }
    to("{:s}", eig.table[i].format_row(j));
  }
  to("{:s}", eig.table.print_row_separator);
  Mpi::Print("{}", fmt::to_string(buf));
  eig.table.col_options = save_defaults;
}

EigenSolver::EigenPostPrinter::EigenPostPrinter(bool do_measurement, bool root,
                                                const fs::path &post_dir, int n_post)
  : root_{root}, do_measurement_(do_measurement),
    stdout_int_print_width(1 + static_cast<int>(std::log10(n_post)))
{
  // Note: we switch to n_eig rather than n_conv for padding since we don't know n_conv
  // until solve
  if (!do_measurement_ || !root_)
  {
    return;
  }
  eig = TableWithCSVFile(post_dir / "eig.csv");
  eig.table.reserve(n_post, 6);
  eig.table.insert_column(Column("idx", "m", 0, {}, {}, ""));
  eig.table.insert_column("f_re", "Re{f} (GHz)");
  eig.table.insert_column("f_im", "Im{f} (GHz)");
  eig.table.insert_column("q", "Q");
  eig.table.insert_column("err_back", "Error (Bkwd.)");
  eig.table.insert_column("err_abs", "Error (Abs.)");
  eig.AppendHeader();
}

void EigenSolver::EigenPostPrinter::AddMeasurement(int eigen_print_idx,
                                                   const PostOperator &post_op,
                                                   double error_bkwd, double error_abs,
                                                   const IoData &iodata)
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  std::complex<double> f =
      iodata.DimensionalizeValue(VT::FREQUENCY, post_op.GetFrequency());
  double Q = (f.imag() == 0.0) ? mfem::infinity() : 0.5 * std::abs(f) / std::abs(f.imag());

  eig.table["idx"] << eigen_print_idx;
  eig.table["f_re"] << f.real();
  eig.table["f_im"] << f.imag();
  eig.table["q"] << Q;
  eig.table["err_back"] << error_bkwd;
  eig.table["err_abs"] << error_abs;

  eig.AppendRow();
  PrintStdoutRow(eig.table.n_rows() - 1);
}

EigenSolver::PortsPostPrinter::PortsPostPrinter(bool do_measurement, bool root,
                                                const fs::path &post_dir,
                                                const LumpedPortOperator &lumped_port_op,
                                                int n_expected_rows)
  : root_{root}, do_measurement_{
                     do_measurement                  //
                     && (lumped_port_op.Size() > 0)  //
                 }
{
  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;
  port_V = TableWithCSVFile(post_dir / "port-V.csv");
  port_V.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_V.table.insert_column(Column("idx", "m", 0, {}, {}, ""));

  port_I = TableWithCSVFile(post_dir / "port-I.csv");
  port_I.table.reserve(n_expected_rows, lumped_port_op.Size());
  port_I.table.insert_column(Column("idx", "m", 0, {}, {}, ""));

  for (const auto &[idx, data] : lumped_port_op)
  {
    port_V.table.insert_column(format("re{}", idx), format("Re{{V[{}]}} (V)", idx));
    port_V.table.insert_column(format("im{}", idx), format("Im{{V[{}]}} (V)", idx));

    port_I.table.insert_column(format("re{}", idx), format("Re{{I[{}]}} (A)", idx));
    port_I.table.insert_column(format("im{}", idx), format("Im{{I[{}]}} (A)", idx));
  }
  port_V.AppendHeader();
  port_I.AppendHeader();
}

void EigenSolver::PortsPostPrinter::AddMeasurement(int eigen_print_idx,
                                                   const PostOperator &post_op,
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
  port_V.table["idx"] << eigen_print_idx;
  port_I.table["idx"] << eigen_print_idx;

  auto unit_V = iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
  auto unit_A = iodata.DimensionalizeValue(VT::CURRENT, 1.0);

  for (const auto &[idx, data] : lumped_port_op)
  {
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

EigenSolver::EPRPostPrinter::EPRPostPrinter(bool do_measurement, bool root,
                                            const fs::path &post_dir,
                                            const LumpedPortOperator &lumped_port_op,
                                            int n_expected_rows)
  : root_{root}, do_measurement_EPR_(do_measurement  //
                                     && lumped_port_op.Size() > 0),
    do_measurement_Q_(do_measurement_EPR_)
{
  // Mode EPR for lumped inductor elements:

  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.L) > 0.)
    {
      ports_with_L.push_back(idx);
    }
    if (std::abs(data.R) > 0.)
    {
      ports_with_R.push_back(idx);
    }
  }

  do_measurement_EPR_ = do_measurement_EPR_ && !ports_with_L.empty();
  do_measurement_Q_ = do_measurement_Q_ && !ports_with_R.empty();

  if (!root_ || (!do_measurement_EPR_ && !do_measurement_Q_))
  {
    return;
  }
  using fmt::format;

  if (do_measurement_EPR_)
  {
    port_EPR = TableWithCSVFile(post_dir / "port-EPR.csv");
    port_EPR.table.reserve(n_expected_rows, 1 + ports_with_L.size());
    port_EPR.table.insert_column(Column("idx", "m", 0, {}, {}, ""));
    for (const auto idx : ports_with_L)
    {
      port_EPR.table.insert_column(format("p_{}", idx), format("p[{}]", idx));
    }
    port_EPR.AppendHeader();
  }
  if (do_measurement_Q_)
  {
    port_Q = TableWithCSVFile(post_dir / "port-Q.csv");
    port_Q.table.reserve(n_expected_rows, 1 + ports_with_R.size());
    port_Q.table.insert_column(Column("idx", "m", 0, {}, {}, ""));
    for (const auto idx : ports_with_R)
    {
      port_Q.table.insert_column(format("Ql_{}", idx), format("Q_ext[{}]", idx));
      port_Q.table.insert_column(format("Kl_{}", idx), format("κ_ext[{}] (GHz)", idx));
    }
    port_Q.AppendHeader();
  }
}

void EigenSolver::EPRPostPrinter::AddMeasurementEPR(
    double eigen_print_idx, const PostOperator &post_op,
    const LumpedPortOperator &lumped_port_op, const IoData &iodata)
{
  if (!do_measurement_EPR_ || !root_)
  {
    return;
  }
  using fmt::format;

  double E_elec = post_op.GetEFieldEnergy();
  double E_cap = post_op.GetLumpedCapacitorEnergy();
  double E_m = E_elec + E_cap;

  port_EPR.table["idx"] << eigen_print_idx;
  for (const auto idx : ports_with_L)
  {
    port_EPR.table[format("p_{}", idx)] << post_op.GetInductorParticipation(idx, E_m);
  }
  port_EPR.AppendRow();
}

void EigenSolver::EPRPostPrinter::AddMeasurementQ(double eigen_print_idx,
                                                  const PostOperator &post_op,
                                                  const LumpedPortOperator &lumped_port_op,
                                                  const IoData &iodata)
{
  if (!do_measurement_Q_ || !root_)
  {
    return;
  }
  using fmt::format;
  using VT = IoData::ValueType;

  auto omega = post_op.GetFrequency();
  double E_elec = post_op.GetEFieldEnergy();
  double E_cap = post_op.GetLumpedCapacitorEnergy();
  double E_m = E_elec + E_cap;

  port_Q.table["idx"] << eigen_print_idx;
  for (const auto idx : ports_with_R)
  {
    double Kl = post_op.GetExternalKappa(idx, E_m);
    double Ql = (Kl == 0.0) ? mfem::infinity() : omega.real() / std::abs(Kl);

    port_Q.table[format("Ql_{}", idx)] << Ql;
    port_Q.table[format("Kl_{}", idx)] << iodata.DimensionalizeValue(VT::FREQUENCY, Kl);
  }
  port_Q.AppendRow();
}

void EigenSolver::EPRPostPrinter::AddMeasurement(double eigen_print_idx,
                                                 const PostOperator &post_op,
                                                 const LumpedPortOperator &lumped_port_op,
                                                 const IoData &iodata)
{
  AddMeasurementEPR(eigen_print_idx, post_op, lumped_port_op, iodata);
  AddMeasurementQ(eigen_print_idx, post_op, lumped_port_op, iodata);
}

EigenSolver::PostprocessPrintResults::PostprocessPrintResults(bool root,
                                                              const fs::path &post_dir,
                                                              const PostOperator &post_op,
                                                              const SpaceOperator &space_op,
                                                              int n_post_)
  : n_post(n_post_), write_paraview_fields(n_post_ > 0),
    domains{true, root, post_dir, post_op, "m", n_post},
    surfaces{true, root, post_dir, post_op, "m", n_post},
    probes{true, root, post_dir, post_op, "m", n_post}, eigen{true, root, post_dir, n_post},
    epr{true, root, post_dir, space_op.GetLumpedPortOp(), n_post},
    error_indicator{true, root, post_dir}
{
}

void EigenSolver::PostprocessPrintResults::PostprocessStep(const IoData &iodata,
                                                           const PostOperator &post_op,
                                                           const SpaceOperator &space_op,
                                                           int step, double error_abs,
                                                           double error_bkward)
{
  int eigen_print_idx = step + 1;

  domains.AddMeasurement(eigen_print_idx, post_op, iodata);
  surfaces.AddMeasurement(eigen_print_idx, post_op, iodata);
  probes.AddMeasurement(eigen_print_idx, post_op, iodata);
  eigen.AddMeasurement(eigen_print_idx, post_op, error_bkward, error_abs, iodata);
  epr.AddMeasurement(eigen_print_idx, post_op, space_op.GetLumpedPortOp(), iodata);
  // The internal GridFunctions in PostOperator have already been set:
  if (write_paraview_fields && step < n_post)
  {
    Mpi::Print("\n");
    post_op.WriteFields(step, eigen_print_idx);
    Mpi::Print(" Wrote mode {:d} to disk\n", eigen_print_idx);
  }
}

void EigenSolver::PostprocessPrintResults::PostprocessFinal(const PostOperator &post_op,
                                                            const ErrorIndicator &indicator)
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
