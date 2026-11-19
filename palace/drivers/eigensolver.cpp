// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <complex>
#include <fstream>
#include <map>
#include <numbers>
#include <set>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mfem.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/gridfunction.hpp"
#include "fem/lumpedelement.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
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
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// DIAGNOSTIC: quadratic forms of a computed eigenvector x (normalized in the mass inner
// product) of the quadratic eigenvalue problem P(λ) x = (K + λ C + λ² M) x = 0, λ = iω,
// and the split of the damping form xᴴ C x by boundary.
struct EigenvectorDiagnostics
{
  // One-sided (xᴴ A x) and two-sided (xᵀ A x) forms of K, C, and M.
  std::complex<double> kH = 0.0, cH = 0.0, mH = 0.0, kT = 0.0, cT = 0.0, mT = 0.0;

  // xᵀ P(λ) x for the Ritz value λ, xᵀ P'(λ) x = xᵀ (C + 2 λ M) x, ||P(λ) x||₂, and ||x||₂.
  std::complex<double> xTr = 0.0, xTdP = 0.0;
  double res = 0.0, xnorm = 0.0;

  // xᴴ C_j x for each boundary damping group.
  std::vector<double> cj;

  // Lumped port element voltages V_e = (1/w_e) ∫ E ⋅ l̂_e dS.
  std::vector<std::complex<double>> Ve;
};

// Root of m μ² + c μ + k = 0 closest to the given value.
std::complex<double> QuadraticRoot(std::complex<double> m, std::complex<double> c,
                                   std::complex<double> k, std::complex<double> guess)
{
  const auto d = std::sqrt(c * c - 4.0 * m * k);
  const auto r1 = (-c + d) / (2.0 * m), r2 = (-c - d) / (2.0 * m);
  return (std::abs(r1 - guess) < std::abs(r2 - guess)) ? r1 : r2;
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
    // ABC contributes a pole term f(λ)·M_ff to A2(λ), with f(λ) = -0.5/λ, that a
    // polynomial (K' + λC' + λ²M') seed cannot fit accurately. AddFrozenPole removes the
    // pole's interpolated contribution and re-adds it frozen at the target into the
    // K-block. The freeze is deliberately unconditional (no DetermineFrozen): the fit's
    // pointwise window error is small (the pole at λ = 0 sits outside the window), but its
    // curvature term injects a fictitious rank-deficient λ² contribution into the seed's
    // M-block that displaces spurious roots into the unphysical half-plane — the failure
    // mode is structural, not a pointwise approximation error, so the pointwise metric
    // must not be allowed to choose the fit here.
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

  // DIAGNOSTIC: boundary damping groups (each lumped port and lumped port element, each
  // impedance boundary with nonzero resistance, and each conductor loss postprocessing
  // index), and the quadratic forms of each converged eigenvector.
  std::vector<std::string> group_names;
  std::vector<std::unique_ptr<ComplexOperator>> group_ops;
  if (C && !has_A2)
  {
    auto AddGroup = [&](const std::string &name, const std::vector<int> &attrs)
    {
      auto Cj = space_op.GetBdrDampingMatrixRestricted(attrs, Operator::DIAG_ZERO);
      if (Cj)
      {
        group_names.push_back(name);
        group_ops.push_back(std::move(Cj));
      }
    };
    for (const auto &[idx, data] : iodata.boundaries.lumpedport)
    {
      std::vector<int> port_attrs;
      for (std::size_t e = 0; e < data.elements.size(); e++)
      {
        const auto &attrs = data.elements[e].attributes;
        port_attrs.insert(port_attrs.end(), attrs.begin(), attrs.end());
        if (data.elements.size() > 1)
        {
          AddGroup(fmt::format("port{:d}_e{:d}", idx, e + 1), attrs);
        }
      }
      AddGroup(fmt::format("port{:d}", idx), port_attrs);
    }
    for (const auto &data : iodata.boundaries.impedance)
    {
      if (data.Rs != 0.0)
      {
        AddGroup(fmt::format("imp[{}]", fmt::join(data.attributes, "+")), data.attributes);
      }
    }
    for (const auto &[idx, data] : iodata.boundaries.postpro.conductor)
    {
      AddGroup(fmt::format("cond{:d}", idx), data.attributes);
    }
  }
  // Lumped port element voltage functional, binned by element: each element contributes
  // V_e = (1/w_e) ∫ E ⋅ l̂_e dS with resistance R_e = n R (elements add in parallel).
  struct PortElement
  {
    int port, elem;
    double R_e;
  };
  std::vector<PortElement> port_elems;
  std::unique_ptr<SurfaceFunctional> elem_vfunc;
  mfem::Array<int> elem_attr_to_bin;
  if (C && !has_A2)
  {
    const auto &pmesh = space_op.GetNDSpace().GetParMesh();
    const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    elem_attr_to_bin.SetSize(bdr_attr_max);
    elem_attr_to_bin = -1;
    std::vector<SurfaceModeCoefficient> modes;
    bool ok = true;
    for (const auto &[idx, data] : space_op.GetLumpedPortOp())
    {
      if (std::abs(data.R) == 0.0)
      {
        continue;
      }
      for (std::size_t k = 0; k < data.elems.size(); k++)
      {
        const auto *uniform = dynamic_cast<const UniformElementData *>(data.elems[k].get());
        if (!uniform)
        {
          ok = false;
          continue;
        }
        SurfaceModeCoefficient mode;
        mode.type = SurfaceModeCoefficient::Type::UNIFORM;
        mode.attr_list = uniform->GetAttrList();
        mode.scale = 1.0 / uniform->GetGeometryWidth();
        for (int d = 0; d < 3; d++)
        {
          mode.direction[d] =
              (d < uniform->GetDirection().Size()) ? uniform->GetDirection()(d) : 0.0;
        }
        for (int attr : mode.attr_list)
        {
          if (attr <= 0 || attr > bdr_attr_max || elem_attr_to_bin[attr - 1] >= 0)
          {
            ok = false;
          }
          else
          {
            elem_attr_to_bin[attr - 1] = static_cast<int>(port_elems.size());
          }
        }
        port_elems.push_back({idx, static_cast<int>(k) + 1,
                              data.R * static_cast<double>(data.elems.size())});
        modes.push_back(std::move(mode));
      }
    }
    if (ok && !modes.empty())
    {
      elem_vfunc = std::make_unique<SurfaceFunctional>(space_op.GetMaterialOp().GetMesh(),
                                                       space_op.GetNDSpace().Get(), modes);
      if (!elem_vfunc->IsValid())
      {
        elem_vfunc.reset();
      }
    }
  }
  std::vector<EigenvectorDiagnostics> diags(num_report);
  if (C && !has_A2)
  {
    const auto comm = space_op.GetComm();
    GridFunction Egf(space_op.GetNDSpace(), true);
    ComplexVector x(Curl.Width()), Kx(Curl.Width()), Cx(Curl.Width()), Mx(Curl.Width()),
        r(Curl.Width()), y(Curl.Width());
    for (auto *v : {&x, &Kx, &Cx, &Mx, &r, &y})
    {
      v->UseDevice(true);
    }
    for (int i = 0; i < num_report; i++)
    {
      const std::complex<double> lambda = eigen->GetEigenvalue(i);
      eigen->GetEigenvector(i, x);
      K->Mult(x, Kx);
      C->Mult(x, Cx);
      M->Mult(x, Mx);
      r = Kx;
      r.AXPY(lambda, Cx);
      r.AXPY(lambda * lambda, Mx);
      y = Cx;
      y.AXPY(2.0 * lambda, Mx);
      std::complex<double> loc[10] = {Kx.Dot(x),          Cx.Dot(x),
                                      Mx.Dot(x),          x.TransposeDot(Kx),
                                      x.TransposeDot(Cx), x.TransposeDot(Mx),
                                      x.TransposeDot(r),  x.TransposeDot(y),
                                      r.Dot(r),           x.Dot(x)};
      Mpi::GlobalSum(10, loc, comm);
      auto &d = diags[i];
      d.kH = loc[0];
      d.cH = loc[1];
      d.mH = loc[2];
      d.kT = loc[3];
      d.cT = loc[4];
      d.mT = loc[5];
      d.xTr = loc[6];
      d.xTdP = loc[7];
      d.res = std::sqrt(std::abs(loc[8]));
      d.xnorm = std::sqrt(std::abs(loc[9]));
      d.cj.resize(group_ops.size());
      for (std::size_t j = 0; j < group_ops.size(); j++)
      {
        group_ops[j]->Mult(x, y);
        std::complex<double> cj = y.Dot(x);
        Mpi::GlobalSum(1, &cj, comm);
        d.cj[j] = cj.real();
      }
      if (elem_vfunc)
      {
        Egf.Real().SetFromTrueDofs(x.Real());
        Egf.Imag().SetFromTrueDofs(x.Imag());
        Egf.Real().ExchangeFaceNbrData();
        Egf.Imag().ExchangeFaceNbrData();
        d.Ve = elem_vfunc->EvalModeOverlapByAttribute(Egf, elem_attr_to_bin,
                                                      static_cast<int>(port_elems.size()));
      }
    }
  }
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
    std::complex<double> omega = eigenvalues[i];
    double error_bkwd = errors_bkwd[i];
    double error_abs = errors_abs[i];
    if (!has_C && !has_A2)
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

  // DIAGNOSTIC: write the eigenvalue consistency checks and the damping split by boundary.
  if (has_C && !has_A2 && Mpi::Root(space_op.GetComm()))
  {
    auto ToFreq = [&](std::complex<double> lambda)
    {
      // λ = iω → f = ω / 2π [GHz].
      const std::complex<double> omega = lambda / 1i;
      return iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
             (2 * std::numbers::pi);
    };
    auto QFromLambda = [](std::complex<double> lambda)
    {
      const std::complex<double> omega = lambda / 1i;
      return 0.5 * std::abs(omega) / std::abs(omega.imag());
    };
    std::ofstream fo(std::string(iodata.problem.output) + "/eig-diag.csv");
    fo << "m,Re{f} (GHz),Im{f} (GHz),Q,Re{f_TS} (GHz),Im{f_TS} (GHz),Q_TS,"
          "Re{f_EB} (GHz),Im{f_EB} (GHz),Q_EB,Re{df_Newton} (GHz),Im{df_Newton} (GHz),"
          "|df| bound (GHz),||r||/||x||,||x||^2,Im{x^T C x}/Re,Im{x^H M x}";
    for (const auto &name : group_names)
    {
      fo << ",Q_damp[" << name << "]";
    }
    fo << ",Q_damp[sum of ports],Q_damp[all C]";
    for (const auto &pe : port_elems)
    {
      fo << ",Re{V[" << pe.port << "_e" << pe.elem << "]},Im{V[" << pe.port << "_e"
         << pe.elem << "]},Q_V[" << pe.port << "_e" << pe.elem << "]";
    }
    {
      std::set<int> ports_seen;
      for (const auto &pe : port_elems)
      {
        if (ports_seen.insert(pe.port).second)
        {
          fo << ",Q_Velems[" << pe.port << "],Q_Vmean[" << pe.port << "]";
        }
      }
    }
    fo << "\n";
    Mpi::Print("\nEigenvector diagnostics (TS = two-sided Rayleigh quotient root, EB = "
               "one-sided root = energy balance, Newton = -x^T r / x^T P' x):\n");
    for (int i = 0; i < num_report; i++)
    {
      const auto &d = diags[i];
      const std::complex<double> lambda = eigenvalues[i];
      const auto lTS = QuadraticRoot(d.mT, d.cT, d.kT, lambda);
      const auto lEB = QuadraticRoot(d.mH.real(), d.cH.real(), d.kH.real(), lambda);
      const auto dl = -d.xTr / d.xTdP;
      const double dl_bound = d.res * d.xnorm / std::abs(d.xTdP);
      const auto f = ToFreq(lambda), fTS = ToFreq(lTS), fEB = ToFreq(lEB);
      const auto df = ToFreq(lambda + dl) - f;
      const double df_bound =
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(dl_bound) /
          (2 * std::numbers::pi);
      // Damping quality factor of each group, consistent with the energy balance:
      // Q_j = Re{ω_EB} xᴴ M x / xᴴ C_j x, such that 1/Q_EB = ∑_j 1/Q_j for a partition.
      const double w_EB = (lEB / 1i).real();
      auto Qj = [&](double c)
      { return (c > 0.0) ? w_EB * d.mH.real() / c : mfem::infinity(); };
      double c_ports = 0.0;
      for (std::size_t j = 0; j < group_names.size(); j++)
      {
        if (group_names[j].rfind("port", 0) == 0 &&
            group_names[j].find("_e") == std::string::npos)
        {
          c_ports += d.cj[j];
        }
      }
      fo << fmt::format(
          "{:d},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},"
          "{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e},{:.15e}",
          i + 1, f.real(), f.imag(), QFromLambda(lambda), fTS.real(), fTS.imag(),
          QFromLambda(lTS), fEB.real(), fEB.imag(), QFromLambda(lEB), df.real(), df.imag(),
          df_bound, d.res / d.xnorm, d.xnorm * d.xnorm, d.cH.imag() / d.cH.real(),
          d.mH.imag());
      for (auto c : d.cj)
      {
        fo << fmt::format(",{:.15e}", Qj(c));
      }
      fo << fmt::format(",{:.15e},{:.15e}", Qj(c_ports), Qj(d.cH.real()));
      // V-based port damping: per element |V_e|² / R_e, summed over the elements of a port
      // (uniform field in each element), and |mean(V_e)|² / R (the circuit postprocessing
      // of port-Q.csv).
      std::map<int, double> c_Velems;
      std::map<int, std::complex<double>> V_sum;
      std::map<int, int> n_elems;
      for (std::size_t k = 0; k < port_elems.size(); k++)
      {
        const auto &pe = port_elems[k];
        const auto V = (k < d.Ve.size()) ? d.Ve[k] : 0.0;
        const double c = std::norm(V) / pe.R_e;
        fo << fmt::format(",{:.15e},{:.15e},{:.15e}", V.real(), V.imag(), Qj(c));
        c_Velems[pe.port] += c;
        V_sum[pe.port] += V;
        n_elems[pe.port]++;
      }
      {
        std::set<int> ports_seen;
        for (const auto &pe : port_elems)
        {
          if (ports_seen.insert(pe.port).second)
          {
            const int n = n_elems[pe.port];
            const double R = pe.R_e / n;
            const auto Vmean = V_sum[pe.port] / static_cast<double>(n);
            fo << fmt::format(",{:.15e},{:.15e}", Qj(c_Velems[pe.port]),
                              Qj(std::norm(Vmean) / R));
          }
        }
      }
      fo << "\n";
      Mpi::Print(
          " m = {:d}: f = {:.12e}{:+.6e}i GHz, Q = {:.6e} | Q_TS = {:.6e}, Q_EB = "
          "{:.6e} | df_Newton = {:+.3e}{:+.3e}i GHz (bound {:.3e}) | Q_damp(ports) = "
          "{:.6e}, Q_damp(all C) = {:.6e}\n",
          i + 1, f.real(), f.imag(), QFromLambda(lambda), QFromLambda(lTS),
          QFromLambda(lEB), df.real(), df.imag(), df_bound, Qj(c_ports), Qj(d.cH.real()));
    }
  }
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
