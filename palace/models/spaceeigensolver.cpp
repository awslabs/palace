// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceeigensolver.hpp"

#include <numbers>
#include <mfem.hpp>
#include "linalg/divfree.hpp"
#include "linalg/nleps.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

SpaceEigenSolver::SpaceEigenSolver(const IoData &iodata, SpaceOperator &space_op,
                                   const ComplexOperator &K, const ComplexOperator *C,
                                   const ComplexOperator &M, ComplexKspSolver &ksp)
  : iodata(iodata), space_op(space_op), K(K), C(C), M(M), ksp(ksp)
{
  const auto &eigenmode = iodata.solver.eigenmode;
  const double target = eigenmode.target;

  // Check if there are nonlinear terms and, if so, setup interpolation operator.
  funcA2 = [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
  };
  // As funcA2, plus the matrix-free wave-port modal correction W. Used by the SLP and
  // Quasi-Newton solves (the polynomial seed keeps funcA2, since W cannot enter
  // BuildParSumOperator).
  funcA2_full = [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  };
  funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                      std::complex<double> a2,
                      std::complex<double> a3) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, a3); };
  A2 = funcA2(1i * target);
  has_A2 = (A2 != nullptr);

  // Freeze the wave-port modal reference at the target so funcA2_full's complex-ω
  // correction W can extrapolate k_n(ω) around it (used by both the SLP and Quasi-Newton
  // solves).
  if (space_op.GetWavePortOp().Size() > 0)
  {
    space_op.GetWavePortOp().InitializeModalReference(target);
  }

  // Extend K, C, M operators with interpolated A2 operator.
  // K' = K + A2_0, C' = C + A2_1, M' = M + A2_2
  nonlinear_type = eigenmode.nonlinear_type;
  // SLP is only realized through SLEPc's NEP path. Resolve it against the selected backend
  // here, before the HYBRID interpolation below: an ARPACK backend (or a build without
  // SLEPc) has no SLP solver, and switching later would leave Kp/Cp/Mp unbuilt and crash.
  bool slp_available = true;
#if !defined(PALACE_WITH_SLEPC)
  slp_available = false;
#endif
  if (eigenmode.type == EigenSolverBackend::ARPACK)
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
    const double target_max = eigenmode.target_upper;
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
    Kp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {&K, A2_0.get()});
    Cp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {C, A2_1.get()});
    Mp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {&M, A2_2.get()});
  }

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  // The SLP fallback for an unavailable backend (or no nonlinear term) is resolved above,
  // before the HYBRID interpolation.
  const EigenSolverBackend type = eigenmode.type;
  quadratic = (C || has_A2);
  nonlinear_slp = (nonlinear_type == NonlinearEigenSolver::SLP);
  eigen =
      BuildEigenvalueSolver(space_op.GetComm(), eigenmode, iodata.solver.linear.gs_orthog,
                            iodata.problem.verbose, quadratic, nonlinear_slp);
  EigenvalueSolver::ScaleType scale = eigenmode.scale ? EigenvalueSolver::ScaleType::NORM_2
                                                      : EigenvalueSolver::ScaleType::NONE;
  if (nonlinear_slp)
  {
    if (C)
    {
      eigen->SetOperators(K, *C, M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      eigen->SetOperators(K, M, EigenvalueSolver::ScaleType::NONE);
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
      eigen->SetOperators(K, *C, M, scale);
    }
    else
    {
      eigen->SetOperators(K, M, scale);
    }
  }
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    // The linearized seed solve for the hybrid approach uses a looser tolerance than the
    // final (Quasi-Newton refined) eigenpairs.
    eigen->SetTol(eigenmode.linear_tol);
  }
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  if (eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, &M);
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, &K, &M);
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
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

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (eigenmode.init_v0)
  {
    ComplexVector v0;
    if (eigenmode.init_v0_const)
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
  SetEigenSolverShiftInvert(*eigen, type, target, quadratic, nonlinear_slp);

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    // Invert the same W-aware polynomial used by the seed eigensolver.
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, Kp.get(),
                                 Cp.get(), Mp.get());
  }
  else
  {
    // SLP applies the exact nonlinear operator at the target; linear problems have no A2.
    A2_shift = nonlinear_slp ? funcA2_full(1i * target) : std::move(A2);
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, &K, C,
                                 &M, A2_shift.get());
  }
  P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  ksp.SetOperators(*A, *P);
  eigen->SetLinearSolver(ksp);
}

SpaceEigenSolver::~SpaceEigenSolver() = default;

int SpaceEigenSolver::Solve()
{
  // Eigenvalue problem solve.
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
    const auto &eigenmode = iodata.solver.eigenmode;
    Mpi::Print("\n Refining eigenvalues with Quasi-Newton solver\n");
    auto qn = std::make_unique<QuasiNewtonSolver>(space_op.GetComm(), std::move(eigen),
                                                  num_conv, iodata.problem.verbose,
                                                  eigenmode.refine_nonlinear);
    qn->SetTol(eigenmode.tol);
    qn->SetMaxIter(eigenmode.nonlinear_max_it);
    if (C)
    {
      qn->SetOperators(K, *C, M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      qn->SetOperators(K, M, EigenvalueSolver::ScaleType::NONE);
    }
    qn->SetExtraSystemMatrix(funcA2_full);
    qn->SetPreconditionerUpdate(funcP);
    qn->SetNumModes(eigenmode.n, eigenmode.max_size);
    qn->SetPreconditionerLag(eigenmode.preconditioner_lag,
                             eigenmode.preconditioner_lag_tol);
    qn->SetMaxRestart(eigenmode.max_restart);
    qn->SetLinearSolver(ksp);
    qn->SetShiftInvert(1i * eigenmode.target);
    eigen = std::move(qn);

    // Suppress wave port output during nonlinear eigensolver iterations.
    space_op.GetWavePortOp().SetSuppressOutput(true);
    num_conv = eigen->Solve();
    space_op.GetWavePortOp().SetSuppressOutput(false);
  }
  return num_conv;
}

void SpaceEigenSolver::NormalizeEigenvectors(int num_eig)
{
  if (!KM)
  {
    // Normalize every converged eigenvector with respect to the mass matrix (unit electric
    // field energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, &M);
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_eig);
  }
}

}  // namespace palace
