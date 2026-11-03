// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <map>
#include <numbers>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "fem/shapederivativeform.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/sensitivityexport.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/enum_string.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// The properties of a material which the eigenmodes depend on.
std::vector<MaterialProperty> GetMaterialProperties(const config::MaterialData &material)
{
  auto IsNonzero = [](const auto &tensor)
  { return std::ranges::any_of(tensor.s, [](double v) { return v != 0.0; }); };
  std::vector<MaterialProperty> properties = {MaterialProperty::PERMITTIVITY};
  if (IsNonzero(material.tandelta))
  {
    properties.push_back(MaterialProperty::LOSS_TANGENT);
  }
  if (IsNonzero(material.sigma))
  {
    properties.push_back(MaterialProperty::CONDUCTIVITY);
  }
  properties.push_back(MaterialProperty::PERMEABILITY);
  if (material.lambda_L > 0.0)
  {
    properties.push_back(MaterialProperty::LONDON_DEPTH);
  }
  return properties;
}

// The unconjugated product xᵀ A x.
std::complex<double> TransposeProduct(MPI_Comm comm, const ComplexOperator &A,
                                      const ComplexVector &x, ComplexVector &y)
{
  A.Mult(x, y);
  std::complex<double> dot = {
      linalg::LocalDot(x.Real(), y.Real()) - linalg::LocalDot(x.Imag(), y.Imag()),
      linalg::LocalDot(x.Real(), y.Imag()) + linalg::LocalDot(x.Imag(), y.Real())};
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

}  // namespace

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
  const bool has_C = (C != nullptr);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  // Check if there are nonlinear terms and, if so, setup interpolation operator.
  auto funcA2 = [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
  };
  // As funcA2, plus the matrix-free wave-port modal correction W. Used by the SLP and
  // Quasi-Newton solves (the polynomial seed keeps funcA2, since W cannot enter
  // BuildParSumOperator).
  auto funcA2_full =
      [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  };
  auto funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                           std::complex<double> a2,
                           std::complex<double> a3) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, a3); };
  const double target = iodata.solver.eigenmode.target;
  auto A2 = funcA2(1i * target);
  bool has_A2 = (A2 != nullptr);

  // Freeze the wave-port modal reference at the target so funcA2_full's complex-ω
  // correction W can extrapolate k_n(ω) around it (used by both the SLP and Quasi-Newton
  // solves).
  if (space_op.GetWavePortOp().Size() > 0)
  {
    space_op.GetWavePortOp().InitializeModalReference(target);
  }

  // Extend K, C, M operators with interpolated A2 operator.
  // K' = K + A2_0, C' = C + A2_1, M' = M + A2_2
  std::unique_ptr<ComplexOperator> Kp, Cp, Mp;
  std::unique_ptr<Interpolation> interp_op;
  std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2;
  NonlinearEigenSolver nonlinear_type = iodata.solver.eigenmode.nonlinear_type;
  // SLP is only realized through SLEPc's NEP path. Resolve it against the selected backend
  // here, before the HYBRID interpolation below: an ARPACK backend (or a build without
  // SLEPc) has no SLP solver, and switching later would leave Kp/Cp/Mp unbuilt and crash.
  bool slp_available = true;
#if !defined(PALACE_WITH_SLEPC)
  slp_available = false;
#endif
  if (iodata.solver.eigenmode.type == EigenSolverBackend::ARPACK)
  {
    slp_available = false;
  }
  if (nonlinear_type == NonlinearEigenSolver::SLP && !slp_available)
  {
    Mpi::Warning("SLP nonlinear eigensolver requires the SLEPc backend, using Hybrid!\n");
    nonlinear_type = NonlinearEigenSolver::HYBRID;
  }
  if (nonlinear_type == NonlinearEigenSolver::SLP && !has_A2)
  {
    // The SLP nonlinear eigensolver requires a frequency-dependent term A2(λ) (wave ports
    // or absorbing boundaries). Without one the NEP function/Jacobian shells dereference a
    // null A2; fall back to the standard linear/quadratic eigensolver instead.
    Mpi::Warning(
        "SLP nonlinear eigensolver requires a nonlinear system term (wave ports or "
        "absorbing boundaries), none found; using a linear eigensolver!\n");
    nonlinear_type = NonlinearEigenSolver::HYBRID;
  }
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    const double target_max = iodata.solver.eigenmode.target_upper;
    auto interp = std::make_unique<NewtonInterpolationOperator>(funcA2, A2->Width());
    interp->Interpolate(1i * target, 1i * target_max);
    // Frozen-ABC seed for the NLEPS HYBRID polynomial seed pencil. The 2nd-order farfield
    // ABC contributes a pole term f(λ)·M_ff to A2(λ), with f(λ) = -0.5/λ (λ = iω; the
    // negative sign of β is inside M_ff, not in f), that a polynomial (K' + λC' + λ²M')
    // seed cannot fit accurately. AddFrozenPole removes the pole's interpolated
    // contribution and re-adds it frozen at the target into the K-block. The freeze is
    // deliberately unconditional (no DetermineFrozen): the fit's pointwise window error is
    // small (the pole at λ = 0 sits outside the window), but its curvature term injects a
    // fictitious rank-deficient λ² contribution into the seed's M-block that displaces
    // spurious roots into the unphysical half-plane — the failure mode is structural, not a
    // pointwise approximation error, so the pointwise metric must not be allowed to choose
    // the fit here.
    auto M_ff =
        space_op.GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    if (M_ff)
    {
      interp->AddFrozenPole(
          std::move(M_ff), [](std::complex<double> lambda) { return -0.5 / lambda; },
          1i * target);
    }
    // Fit-or-freeze seed for rational impedance boundaries. Each boundary contributes
    // g(λ)·M_b to A2(λ) with g(λ) = λ·D(λ)/N(λ) = P(λ) + R(λ)/N(λ). The polynomial part
    // P is exactly representable in the seed pencil when deg(P) <= 2, so only the strictly
    // proper pole part R/N is a candidate for freezing. We freeze it at the target when
    // that approximates g over the interpolation window better than the polynomial
    // interpolant does; otherwise keep the fit.
    const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
    for (int idx = 0; idx < surf_rz_op.GetNumBoundaries(); idx++)
    {
      const bool proper_split = (surf_rz_op.GetRobinQuotientDegree(idx) <= 2);
      auto f_full = [&surf_rz_op, idx](std::complex<double> lambda)
      { return surf_rz_op.EvalRobinCoefficient(idx, lambda); };
      auto f_frozen = [&surf_rz_op, idx, proper_split](std::complex<double> lambda)
      {
        return proper_split ? surf_rz_op.EvalRobinRemainder(idx, lambda)
                            : surf_rz_op.EvalRobinCoefficient(idx, lambda);
      };
      double fit_err, freeze_err;
      if (interp->DetermineFrozen(f_full, f_frozen, 1i * target, fit_err, freeze_err))
      {
        Mpi::Print(" Freezing rational impedance boundary (attribute {}) pole part in the "
                   "NLEPS seed (fit error {:.2e} > freeze error {:.2e})\n",
                   fmt::join(surf_rz_op.GetAttrList(idx), ", "), fit_err, freeze_err);
        auto M_b = space_op.GetRationalImpedanceBoundaryMassMatrix<ComplexOperator>(
            idx, Operator::DIAG_ZERO);
        interp->AddFrozenPole(std::move(M_b), f_frozen, 1i * target);
      }
    }
    A2_0 = interp->GetInterpolationOperator(0);
    A2_1 = interp->GetInterpolationOperator(1);
    A2_2 = interp->GetInterpolationOperator(2);
    interp_op = std::move(interp);  // retain: A2_0/A2_1/A2_2 reference its operator DAG
    Kp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {K.get(), A2_0.get()});
    Cp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {C.get(), A2_1.get()});
    Mp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {M.get(), A2_2.get()});
  }

  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  ComplexVector B(Curl.Height());
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = iodata.solver.eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C || has_A2)
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
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      slepc = std::make_unique<slepc::SlepcNEPSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
    }
    else
    {
      if (C || has_A2)
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
    }
    slepc->SetOrthogonalization(iodata.solver.linear.gs_orthog == Orthogonalization::MGS,
                                iodata.solver.linear.gs_orthog == Orthogonalization::CGS2);
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    if (C)
    {
      eigen->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      eigen->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    eigen->SetExtraSystemMatrix(funcA2_full);
    eigen->SetPreconditionerUpdate(funcP);
  }
  else
  {
    if (has_A2)
    {
      eigen->SetOperators(*Kp, *Cp, *Mp, scale);
    }
    else if (C)
    {
      eigen->SetOperators(*K, *C, *M, scale);
    }
    else
    {
      eigen->SetOperators(*K, *M, scale);
    }
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  const double tol = (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
                         ? iodata.solver.eigenmode.linear_tol
                         : iodata.solver.eigenmode.tol;
  eigen->SetTol(tol);
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
  if (iodata.solver.linear.divfree_max_it > 0 &&
      !space_op.GetMaterialOp().HasWaveVector() &&
      !space_op.GetMaterialOp().HasLondonDepth())
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
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
  {
    const double f_target =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) /
        (2 * std::numbers::pi);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen->SetShiftInvert(1i * target);
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
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
    if (type == EigenSolverBackend::ARPACK)
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
  std::unique_ptr<ComplexOperator> A;
  std::unique_ptr<ComplexOperator> A2_shift;
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    // Invert the same W-aware polynomial used by the seed eigensolver.
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, Kp.get(),
                                 Cp.get(), Mp.get());
  }
  else
  {
    // SLP applies the exact nonlinear operator at the target; linear problems have no A2.
    A2_shift = (nonlinear_type == NonlinearEigenSolver::SLP) ? funcA2_full(1i * target)
                                                             : std::move(A2);
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, K.get(),
                                 C.get(), M.get(), A2_shift.get());
  }
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  const FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);
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

  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    Mpi::Print("\n Refining eigenvalues with Quasi-Newton solver\n");
    auto qn = std::make_unique<QuasiNewtonSolver>(space_op.GetComm(), std::move(eigen),
                                                  num_conv, iodata.problem.verbose,
                                                  iodata.solver.eigenmode.refine_nonlinear);
    qn->SetTol(iodata.solver.eigenmode.tol);
    qn->SetMaxIter(iodata.solver.eigenmode.nonlinear_max_it);
    if (C)
    {
      qn->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      qn->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    qn->SetExtraSystemMatrix(funcA2_full);
    qn->SetPreconditionerUpdate(funcP);
    qn->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
    qn->SetPreconditionerLag(iodata.solver.eigenmode.preconditioner_lag,
                             iodata.solver.eigenmode.preconditioner_lag_tol);
    qn->SetMaxRestart(iodata.solver.eigenmode.max_restart);
    qn->SetLinearSolver(*ksp);
    qn->SetShiftInvert(1i * target);
    eigen = std::move(qn);

    // Suppress wave port output during nonlinear eigensolver iterations.
    space_op.GetWavePortOp().SetSuppressOutput(true);
    num_conv = eigen->Solve();
    space_op.GetWavePortOp().SetSuppressOutput(false);
  }

  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  const int num_requested = iodata.solver.eigenmode.n;
  const int num_report = num_conv;
  MFEM_VERIFY(num_report >= num_requested, "Eigenmode solve only found "
                                               << num_report << " modes when "
                                               << num_requested << " were requested!");
  if (!KM)
  {
    // Normalize every converged eigenvector with respect to the mass matrix (unit electric
    // field energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_report);
  }
  std::vector<std::complex<double>> eigenvalues(num_report);
  std::vector<double> errors_bkwd(num_report), errors_abs(num_report);
  const bool stage_eigenvectors = post_op.WillWriteFields();
  // Field output can construct large libCEED visualization operators. Only in that case,
  // stage the converged eigenvectors in host memory so the eigensolver and matrices can
  // release their device allocations before postprocessing.
  std::vector<std::vector<std::complex<double>>> eigenvectors;
  if (stage_eigenvectors)
  {
    eigenvectors.assign(num_report, std::vector<std::complex<double>>(Curl.Width()));
  }
  ComplexVector E(Curl.Width());
  E.UseDevice(true);
  for (int i = 0; i < num_report; i++)
  {
    eigenvalues[i] = eigen->GetEigenvalue(i);
    errors_bkwd[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    errors_abs[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (stage_eigenvectors)
    {
      eigen->GetEigenvector(i, E);
      E.Get(eigenvectors[i].data(), E.Size(), false);
    }
  }
  auto GetOmega = [&](int i)
  {
    // The linear EVP has eigenvalue μ = -λ² = ω², the quadratic one λ = iω.
    return (!has_C && !has_A2) ? std::sqrt(eigenvalues[i]) : eigenvalues[i] / 1i;
  };
  if (!iodata.sensitivity.empty())
  {
    std::vector<std::complex<double>> omega(num_report);
    for (int i = 0; i < num_report; i++)
    {
      omega[i] = GetOmega(i);
    }
    PostprocessSensitivities(space_op, *eigen, C.get(), *M, has_A2, divfree.get(), omega);
  }

  if (stage_eigenvectors)
  {
    eigen.reset();
    KM.reset();
    ksp.reset();
    P.reset();
    A.reset();
    divfree.reset();
    Kp.reset();
    Cp.reset();
    Mp.reset();
    A2_0.reset();
    A2_1.reset();
    A2_2.reset();
    A2.reset();
    interp_op.reset();
    K.reset();
    C.reset();
    M.reset();
  }
  Mpi::Print("\n");

  for (int i = 0; i < num_report; i++)
  {
    // Get the eigenvalue and relative error.
    const std::complex<double> omega = GetOmega(i);
    double error_bkwd = errors_bkwd[i];
    double error_abs = errors_abs[i];

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    if (stage_eigenvectors)
    {
      E.Set(eigenvectors[i].data(), E.Size(), false);
    }
    else
    {
      eigen->GetEigenvector(i, E);
    }
    linalg::NormalizePhase(space_op.GetComm(), E);

    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E.
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_report);
    if (i < num_requested)
    {
      estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
    }
  }
  post_op.MeasureFinalize(indicator);
  return {indicator, space_op.GlobalTrueVSize()};
}

void EigenSolver::PostprocessSensitivities(
    SpaceOperator &space_op, const EigenvalueSolver &eigen, const ComplexOperator *C,
    const ComplexOperator &M, bool has_A2, const DivFreeSolver<ComplexVector> *divfree,
    const std::vector<std::complex<double>> &omega) const
{
  // T(λ) = K + λ C + λ² M + A2(λ), with λ = iω, is complex symmetric, so its left
  // eigenvectors are its right ones: a simple eigenvalue changes by
  // dλ = -xᵀ ∂T x / xᵀ T'(λ) x, and no adjoint solve is needed.
  const auto &sensitivity = iodata.sensitivity;
  Mpi::Print("\nComputing the sensitivities of {:d} objective{}\n",
             sensitivity.objectives.size(), (sensitivity.objectives.size() > 1) ? "s" : "");
  SensitivityExport output(iodata, space_op.GetMesh());
  const auto comm = space_op.GetComm();
  const int size = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector x(size), y(size);
  x.UseDevice(true);
  y.UseDevice(true);

  // The derivatives of ω of a mode: the real and imaginary parts of its gradient, and its
  // derivatives with respect to relative changes of the properties of each material and of
  // the elements (R, L, C) of each lumped port.
  struct ModeDerivatives
  {
    std::array<Vector, 2> shape;
    std::vector<std::vector<std::pair<MaterialProperty, std::complex<double>>>> materials;
    std::map<int, std::array<std::optional<std::complex<double>>, 3>> ports;
  };
  std::map<int, ModeDerivatives> derivatives;
  auto GetDerivatives = [&](int m) -> const ModeDerivatives &
  {
    if (auto it = derivatives.find(m); it != derivatives.end())
    {
      return it->second;
    }
    auto &dm = derivatives[m];
    for (std::size_t k = 0; k < omega.size(); k++)
    {
      if (static_cast<int>(k) != m &&
          std::abs(omega[k] - omega[m]) < 1.0e-6 * std::abs(omega[m]))
      {
        Mpi::Warning("Mode {:d} is degenerate, so its derivatives are not defined: the "
                     "result is for its computed eigenvector!\n",
                     m + 1);
        break;
      }
    }
    const std::complex<double> lambda = 1i * omega[m], lambda2 = lambda * lambda;
    eigen.GetEigenvector(m, x);
    if (divfree)
    {
      // The projection of the eigenvalue solve makes x divergence-free with respect to the
      // permittivity, while the eigenvector satisfies Gᵀ (T(λ) - K) x = 0 with the losses:
      // restore its irrotational part, to which the derivatives of the forms are sensitive.
      std::vector<std::complex<double>> a = {1.0};
      std::vector<const ComplexOperator *> ops = {&M};
      if (C)
      {
        a.push_back(1.0 / lambda);
        ops.push_back(C);
      }
      std::unique_ptr<ComplexOperator> A2;
      if (has_A2)
      {
        A2 = space_op.GetExtraSystemMatrix(omega[m], Operator::DIAG_ZERO);
        if (A2)
        {
          a.push_back(1.0 / lambda2);
          ops.push_back(A2.get());
        }
      }
      divfree->Mult(*BuildParSumOperator(a, ops), x);
    }
    std::complex<double> d = 2.0 * lambda * TransposeProduct(comm, M, x, y);
    if (C)
    {
      d += TransposeProduct(comm, *C, x, y);
    }
    if (has_A2)
    {
      // Fourth-order central difference in λ of xᵀ A2(λ) x.
      const double h = 1.0e-3 * std::abs(lambda);
      auto A2 = [&](double t)
      {
        auto A = space_op.GetExtraSystemMatrix((lambda + t) / 1i, Operator::DIAG_ZERO);
        return A ? TransposeProduct(comm, *A, x, y) : 0.0;
      };
      d += (8.0 * (A2(h) - A2(-h)) - (A2(2.0 * h) - A2(-2.0 * h))) / (12.0 * h);
    }

    // dω = -i dλ = c xᵀ ∂T x, with c = i / xᵀ T'(λ) x.
    const std::complex<double> c = 1i / d;
    auto Derivative = [&](const ComplexOperator *dT) -> std::complex<double>
    { return dT ? c * TransposeProduct(comm, *dT, x, y) : 0.0; };
    if (output.HasShape())
    {
      const auto dT = space_op.GetSystemShapeDerivative(1.0, lambda, lambda2, omega[m]);
      Vector gr, gi;
      GetComplexShapeGradient(*dT[0], *dT[1], x, gr, gi);
      dm.shape[0] = gr;
      dm.shape[0] *= c.real();
      dm.shape[0].Add(-c.imag(), gi);
      dm.shape[1] = gr;
      dm.shape[1] *= c.imag();
      dm.shape[1].Add(c.real(), gi);
    }
    if (sensitivity.material_properties)
    {
      for (const auto &material : iodata.domains.materials)
      {
        auto &dmat = dm.materials.emplace_back();
        for (const auto property : GetMaterialProperties(material))
        {
          const auto dT = space_op.GetSystemMaterialDerivative(
              property, material.attributes, 1.0, lambda, lambda2, omega[m]);
          dmat.emplace_back(property, Derivative(dT.get()));
        }
      }
    }
    if (sensitivity.lumped_ports)
    {
      for (const auto &[idx, data] : space_op.GetLumpedPortOp())
      {
        const auto dT = space_op.GetLumpedPortDerivatives(idx, 1.0, lambda, lambda2);
        auto &dport = dm.ports[idx];
        for (int e = 0; e < 3; e++)
        {
          if (dT[e])
          {
            dport[e] = Derivative(dT[e].get());
          }
        }
      }
    }
    return dm;
  };

  const double GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                     (2.0 * std::numbers::pi);
  constexpr std::array<const char *, 3> port_elements = {"R", "L", "C"};
  for (const auto &[idx, data] : sensitivity.objectives)
  {
    // The derivative of f = Re{ω} / 2π or Q = |ω| / 2 |Im{ω}| is ar Re{dω} + ai Im{dω}.
    const int m = data.mode - 1;
    const auto w = omega[m];
    nlohmann::json entry = {
        {"Index", idx}, {"Type", ToString(data.type)}, {"Mode", data.mode}};
    double ar, ai, scale;
    if (data.type == SensitivityObjective::FREQUENCY)
    {
      entry["Value"] = w.real() * GHz;
      entry["Unit"] = "GHz";
      ar = 1.0;
      ai = 0.0;
      scale = GHz;
    }
    else
    {
      if (w.imag() == 0.0)
      {
        Mpi::Warning("Mode {:d} has no losses, so the quality factor of sensitivity "
                     "objective {:d} is not defined!\n",
                     data.mode, idx);
        continue;
      }
      const double Q = 0.5 * std::abs(w) / std::abs(w.imag());
      entry["Value"] = Q;
      entry["Unit"] = "";
      ar = Q * w.real() / std::norm(w);
      ai = Q * (w.imag() / std::norm(w) - 1.0 / w.imag());
      scale = 1.0;
    }
    auto Objective = [&](std::complex<double> dw)
    { return (ar * dw.real() + ai * dw.imag()) * scale; };
    const auto &dm = GetDerivatives(m);
    if (sensitivity.material_properties)
    {
      auto materials = nlohmann::json::array();
      for (std::size_t k = 0; k < dm.materials.size(); k++)
      {
        nlohmann::json material = {{"Index", k + 1},
                                   {"Attributes", iodata.domains.materials[k].attributes}};
        for (const auto &[property, dw] : dm.materials[k])
        {
          material[ToString(property)] = Objective(dw);
        }
        materials.push_back(std::move(material));
      }
      entry["MaterialProperties"] = std::move(materials);
    }
    if (sensitivity.lumped_ports)
    {
      auto ports = nlohmann::json::array();
      for (const auto &[port, dport] : dm.ports)
      {
        nlohmann::json p = {{"Index", port}};
        for (int e = 0; e < 3; e++)
        {
          if (dport[e])
          {
            p[port_elements[e]] = Objective(*dport[e]);
          }
        }
        ports.push_back(std::move(p));
      }
      entry["LumpedPorts"] = std::move(ports);
    }
    Vector g;
    if (output.HasShape())
    {
      g = dm.shape[0];
      g *= ar;
      g.Add(ai, dm.shape[1]);
    }
    output.AddObjective(std::move(entry), output.HasShape() ? &g : nullptr, scale);
  }
  output.Write(post_dir / "sensitivity");
}

}  // namespace palace
