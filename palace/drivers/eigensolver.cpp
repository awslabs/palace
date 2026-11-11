// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fieldtransfer.hpp"
#include "fem/mesh.hpp"
#include "fem/phasereference.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/gmg.hpp"
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

// Starting vector from the eigenpairs (λᵢ, xᵢ) of the previous adaptation iteration, with
// xᵢ given as the real and imaginary parts of the stored fields 2i and 2i + 1: v = Σᵢ wᵢ xᵢ
// and v_lambda = Σᵢ wᵢ λᵢ xᵢ, with weights of unit norm (for each xᵢ) and deterministic
// phases, so that no eigenvector cancels from the combination, and a component of the
// default starting vector r of relative weight eta in v for eigenvectors not found in the
// previous iteration. Returns the number of eigenvectors used.
int GetWarmStartVector(const VectorFieldTransfer &prev_eigenvectors,
                       const std::vector<std::complex<double>> &prev_eigenvalues,
                       SpaceOperator &space_op, const ComplexVector &r, double eta,
                       ComplexVector &v, ComplexVector &v_lambda)
{
  auto &nd_fespace = space_op.GetNDSpace().Get();
  const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists().back();
  const int n = nd_fespace.GetTrueVSize();
  ComplexVector x(n);
  x.UseDevice(true);
  v.SetSize(n);
  v.UseDevice(true);
  v_lambda.SetSize(n);
  v_lambda.UseDevice(true);
  v = 0.0;
  v_lambda = 0.0;
  int num_used = 0;
  for (std::size_t i = 0; i < prev_eigenvalues.size(); i++)
  {
    prev_eigenvectors.GetField(2 * i, nd_fespace, x.Real());
    prev_eigenvectors.GetField(2 * i + 1, nd_fespace, x.Imag());
    linalg::SetSubVector(x, dbc_tdof_list, 0.0);
    const double norm = linalg::Norml2(space_op.GetComm(), x);
    if (norm > 0.0)
    {
      // Phases from the golden ratio sequence.
      const double phase =
          2.0 * std::numbers::pi * std::fmod(0.5 * (std::sqrt(5.0) - 1.0) * i, 1.0);
      const std::complex<double> w = std::polar(1.0 / norm, phase);
      v.AXPY(w, x);
      v_lambda.AXPY(w * prev_eigenvalues[i], x);
      num_used++;
    }
  }
  if (eta > 0.0 && num_used > 0)
  {
    const double norm_r = linalg::Norml2(space_op.GetComm(), r);
    if (norm_r > 0.0)
    {
      v.AXPY(eta * linalg::Norml2(space_op.GetComm(), v) / norm_r, r);
    }
  }
  return num_used;
}

}  // namespace

EigenSolver::EigenSolver(const IoData &iodata, bool root, int size, int num_thread,
                         const char *git_tag)
  : BaseSolver(iodata, root, size, num_thread, git_tag)
{
}

EigenSolver::~EigenSolver() = default;

void EigenSolver::UpdateAfterMeshChange(Mesh &mesh) const
{
  if (!prev_eigenvectors)
  {
    return;
  }
  if (prev_eigenvectors->GetParMesh() != &mesh.Get())
  {
    // The mesh object was replaced.
    prev_eigenvectors.reset();
    prev_eigenvalues.clear();
    return;
  }
  prev_eigenvectors->Update();
}

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
  // Adaptation iterations before the last one, whose eigenvectors only drive the mesh
  // refinement, can use looser tolerances. The last iteration is known in advance from
  // Model.Refinement.MaxIts and MaxSize. When the error indicator stops the adaptation
  // instead, SolveEstimateMarkRefine repeats the solve with the final tolerances.
  const auto &eigenmode = iodata.solver.eigenmode;
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = (refinement.max_it > 0);
  const bool final_solve =
      !use_amr || final_amr_solve_ ||
      (refinement.max_size > 0 && space_op.GlobalTrueVSize() > refinement.max_size);
  intermediate_solve =
      !final_solve && (eigenmode.intermediate_tol > eigenmode.tol ||
                       eigenmode.intermediate_linear_tol > iodata.solver.linear.tol);
  const double eig_tol = intermediate_solve ? eigenmode.intermediate_tol : eigenmode.tol;
  const double lin_tol =
      intermediate_solve ? eigenmode.intermediate_linear_tol : iodata.solver.linear.tol;
  if (intermediate_solve)
  {
    Mpi::Print(" Intermediate adaptation iteration: eigenvalue solver tolerance {:.1e}, "
               "linear solver tolerance {:.1e}\n",
               eig_tol, lin_tol);
  }
  eigen->SetNumModes(eigenmode.n, eigenmode.max_size);
  const double tol = (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
                         ? std::max(eigenmode.linear_tol, eig_tol)
                         : eig_tol;
  eigen->SetTol(tol);
  eigen->SetMaxIter(eigenmode.max_it);
  // Each step of the shift-and-invert Krylov iteration costs a linear solve: test
  // convergence after every step instead of only at the end of each cycle, so that a solve
  // converging in its first cycle no longer costs MaxSize linear solves (SLEPc
  // Krylov-Schur, patched). Each check solves the dense projected problem, whose cost grows
  // as the cube of the subspace size: large subspaces keep the end-of-cycle test.
  constexpr int max_check_subspace_size = 128;
  eigen->SetConvergenceCheckInterval((eigenmode.max_size <= max_check_subspace_size) ? 1
                                                                                     : 0);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // Restrict the search to the quality factors of physical modes, and relax the tolerance
  // of the linear solves as the wanted eigenpairs converge (SLEPc Krylov-Schur solver for
  // linear or linearized quadratic problems).
  const bool slepc_eps =
      (type == EigenSolverBackend::SLEPC && nonlinear_type != NonlinearEigenSolver::SLP &&
       (!(C || has_A2) || eigenmode.pep_linear));
  const bool relax_linear_solves = slepc_eps && eigenmode.linear_relaxation > 0.0;
  if (type == EigenSolverBackend::SLEPC && nonlinear_type != NonlinearEigenSolver::SLP &&
      (C || has_A2) && eigenmode.min_q > 0.5)
  {
    Mpi::Print(" Eigenvalue search excludes modes with Q < {:.3g}\n", eigenmode.min_q);
    eigen->SetMinQualityFactor(eigenmode.min_q);
  }
  if (eigenmode.linear_deflation > 0)
  {
    Mpi::Print(" Linear solver deflation with {:d} harmonic Ritz vectors\n",
               eigenmode.linear_deflation);
  }
  else if (eigenmode.linear_deflation < 0)
  {
    Mpi::Print(" Linear solver deflation: automatic\n");
  }
  if (relax_linear_solves)
  {
    Mpi::Print(" Linear solver tolerance relaxation with c = {:.1e}\n",
               eigenmode.linear_relaxation);
    eigen->SetLinearSolverRelaxation(eigenmode.linear_relaxation);
  }

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
  // projected appropriately. With adaptive mesh refinement, it combines the eigenvectors of
  // the previous iteration, with a small component of the constant or random vector.
  if (eigenmode.init_v0)
  {
    ComplexVector v0;
    if (eigenmode.init_v0_const)
    {
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      space_op.GetRandomInitialVector(v0);
    }
    // A warm start cannot save restarts when the previous solve converged within its first
    // Krylov cycle, whose dimension does not depend on the start vector: keep the default
    // start vector then (the Ritz values of a warm-started solve converging in its first
    // cycle are dominated by the error of its first linear solves).
    int num_warm = 0;
    const bool warm_start = prev_eigenvectors && prev_eigenvectors->IsCurrent(
                                                     space_op.GetNDSpace().GetParMesh());
    if (warm_start && prev_num_its == 1)
    {
      Mpi::Print(
          " Previous adaptation iteration converged within one Krylov cycle: no warm "
          "start\n");
    }
    else if (warm_start)
    {
      constexpr double eta = 1.0e-3;
      ComplexVector v, v_lambda;
      num_warm = GetWarmStartVector(*prev_eigenvectors, prev_eigenvalues, space_op, v0, eta,
                                    v, v_lambda);
      if (num_warm > 0)
      {
        Mpi::Print(
            " Using starting vector from {:d} eigenvectors of the previous adaptation "
            "iteration\n",
            num_warm);
        if (divfree)
        {
          divfree->Mult(v);
          divfree->Mult(v_lambda);
        }
        eigen->SetInitialSpaceBlocks(v, v_lambda);  // Copies the vectors

        // The eigenpairs converging in the first Krylov cycle have Ritz values dominated by
        // the error of the first linear solves, applied to the start vector: solve these
        // with a 100 times smaller tolerance (not below 1e-12, where that error is
        // negligible).
        eigen->SetFirstCycleLinearSolverTol(
            std::max(1.0e-2 * lin_tol, std::min(lin_tol, 1.0e-12)));
        eigen->SetWarmStart(true);
      }
    }
    if (num_warm == 0)
    {
      Mpi::Print(" Using {} starting vector\n",
                 eigenmode.init_v0_const ? "constant" : "random");
      if (divfree)
      {
        divfree->Mult(v0);
      }
      eigen->SetInitialSpace(v0);  // Copies the vector
    }

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }
  prev_eigenvectors.reset();
  prev_eigenvalues.clear();

  // Shift of the spectral transformation: the target, or with adaptive mesh refinement
  // (Solver.Eigenmode.AdaptiveShift) a shift inside the band of the requested eigenvalues
  // of the previous adaptation iteration (see GetEigenvalueShift), with the eigenvalue
  // search (sorting criterion and region) still relative to the target. A solve which
  // converged within two Krylov restarts with the shift at the target, as few as a
  // warm-started one needs, keeps it there.
  constexpr double shift_fraction = 0.45;
  double shift = target;
  if (eigenmode.adaptive_shift && eigen->SupportsTarget() && slepc_eps && !has_A2 &&
      static_cast<int>(prev_shift_eigenvalues.size()) >= eigenmode.n &&
      (prev_shifted || prev_shift_its > 2))
  {
    shift = GetEigenvalueShift(target, prev_shift_eigenvalues, eigenmode.n, shift_fraction);
  }
  if (!prev_shift_eigenvalues.empty())
  {
    // The eigenvalues of the previous adaptation iteration estimate the wanted ones.
    std::vector<std::complex<double>> estimates;
    for (const auto &l : prev_shift_eigenvalues)
    {
      estimates.push_back((!has_C && !has_A2) ? -l * l : l);
    }
    eigen->SetEigenvalueEstimates(estimates);
  }
  prev_shift_eigenvalues.clear();

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ, and set up the linear solver required for solving
  // systems involving the shifted operator (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during
  // the eigenvalue solve. The preconditioner for complex linear systems is constructed
  // from a real approximation to the complex system matrix.
  std::unique_ptr<ComplexOperator> A, A2_shift, P;
  std::unique_ptr<ComplexKspSolver> ksp;

  // The choice of the coarse solve of the multigrid preconditioner, made at any setup of
  // the linear solver, is kept for the later shifts and adaptation iterations.
  auto KeepCoarseChoice = [&]()
  {
    const auto *gmg = ksp ? dynamic_cast<const GeometricMultigridSolver<ComplexOperator> *>(
                                ksp->GetPreconditioner())
                          : nullptr;
    if (gmg && gmg->GetCoarseSolveChoice())
    {
      coarse_choice = std::make_shared<CoarseSolveChoice>(*gmg->GetCoarseSolveChoice());
    }
  };
  auto ConfigureShiftInvert = [&](double s)
  {
    const double f_s = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(s) /
                       (2 * std::numbers::pi);
    if (s != target)
    {
      const double f_target =
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) /
          (2 * std::numbers::pi);
      Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e}), target {:.3e} GHz ({:.3e})\n",
                 f_s, s, f_target, target);
    }
    else
    {
      Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_s, s);
    }
    if (C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP)
    {
      // Search for eigenvalues closest to λ = iσ.
      eigen->SetShiftInvert(1i * s);
      if (s != target)
      {
        eigen->SetTarget(1i * target);
      }
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
      eigen->SetShiftInvert(s * s);
      if (s != target)
      {
        eigen->SetTarget(target * target);
      }
      if (type == EigenSolverBackend::ARPACK)
      {
        // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
        // will be a large-magnitude positive real number for an eigenvalue μ with
        // frequency close to but below the target σ².
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
      }
      else
      {
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
      }
    }
    KeepCoarseChoice();
    ksp.reset();
    P.reset();
    A.reset();
    if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
    {
      // Invert the same W-aware polynomial used by the seed eigensolver.
      A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i,
                                   Kp.get(), Cp.get(), Mp.get());
    }
    else
    {
      // SLP applies the exact nonlinear operator at the target; linear problems have no
      // A2.
      if (!A2_shift)
      {
        A2_shift = (nonlinear_type == NonlinearEigenSolver::SLP) ? funcA2_full(1i * target)
                                                                 : std::move(A2);
      }
      A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * s, -s * s + 0.0i, K.get(), C.get(),
                                   M.get(), A2_shift.get());
    }
    P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * s, -s * s + 0.0i,
                                                          s + 0.0i);
    ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                             &space_op.GetH1Spaces());
    auto *gmg =
        dynamic_cast<GeometricMultigridSolver<ComplexOperator> *>(ksp->GetPreconditioner());
    if (gmg && coarse_choice)
    {
      gmg->SetCoarseSolveChoice(*coarse_choice);
    }
    ksp->SetOperators(*A, *P);
    KeepCoarseChoice();
    ksp->SetRelTol(lin_tol);
    ksp->SetDeflation(eigenmode.linear_deflation);
    eigen->SetLinearSolver(*ksp);
  };
  ConfigureShiftInvert(shift);

  // Structures for storing and reducing the results of error estimation (the estimator
  // itself is constructed after the eigenvalue solve, once the solver objects have been
  // released).
  std::optional<FluxErrorEstimator<ComplexVector>> estimator;
  ErrorIndicator indicator;

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  bool shift_repeated = false;
  while (shift != target)
  {
    // The modes nearer to the target than the farthest accepted one were at least as
    // dominant in the transformed operator as the accepted ones if the shift is at most in
    // the middle of the target and that eigenvalue. Otherwise (modes entering the requested
    // set lowered it), repeat the solve once with the shift from the eigenvalues just found
    // if it is lower, else (or with fewer converged eigenpairs than requested) with the
    // shift at the target.
    std::vector<std::complex<double>> lambda;
    for (int i = 0; i < num_conv; i++)
    {
      const auto l = eigen->GetEigenvalue(i);
      lambda.push_back((!has_C && !has_A2) ? 1i * std::sqrt(l) : l);
    }
    if (CheckEigenvalueShift(target, shift, lambda, eigenmode.n))
    {
      break;
    }
    double new_shift = target;
    if (num_conv >= eigenmode.n && !shift_repeated)
    {
      new_shift = GetEigenvalueShift(target, lambda, eigenmode.n, shift_fraction);
    }
    shift = (new_shift < shift) ? new_shift : target;
    shift_repeated = true;
    Mpi::Print("\n Repeating the eigenvalue solve with σ = {:.3e} GHz ({})\n",
               iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(shift) /
                   (2 * std::numbers::pi),
               (num_conv < eigenmode.n)
                   ? fmt::format("{:d} of {:d} eigenpairs converged", num_conv, eigenmode.n)
                   : "the requested modes were not the dominant ones for the shift");
    {
      BlockTimer bt_construct(Timer::CONSTRUCT);
      ConfigureShiftInvert(shift);
    }
    eigen->ResetInitialSpace();
    Mpi::Print("\n");
    num_conv = eigen->Solve();
  }
  const int num_its = eigen->GetNumIterations();
  prev_shift_its = num_its;
  prev_shifted = (shift != target);
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
    qn->SetTol(eig_tol);
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

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\n{}\n", ComputeErrorEstimation()
                           ? "Computing solution error estimates and performing "
                             "postprocessing"
                           : "Performing postprocessing");
  const int num_requested = eigenmode.n;
  MFEM_VERIFY(num_conv >= num_requested, "Eigenmode solve only found "
                                             << num_conv << " modes when " << num_requested
                                             << " were requested!");
  // The relaxation of the linear solver tolerance only targets the requested eigenpairs, so
  // further converged ones can be less accurate and are not reported then (the quasi-Newton
  // solver of the hybrid approach refines every eigenpair).
  const bool report_requested =
      relax_linear_solves && !(has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID);
  const int num_report = report_requested ? std::min(num_requested, num_conv) : num_conv;

  // Keep all converged eigenpairs for the starting vector of the next adaptation iteration.
  // They follow the refinement and rebalancing of nonconforming meshes, which are modified
  // in place.
  const bool store_eigenpairs = !final_solve && eigenmode.warm_start && eigenmode.init_v0 &&
                                space_op.GetNDSpace().GetParMesh().Nonconforming();
  const int num_eig = store_eigenpairs ? num_conv : num_report;
  if (!KM)
  {
    // Normalize every converged eigenvector with respect to the mass matrix (unit electric
    // field energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_eig);
  }
  std::vector<std::complex<double>> eigenvalues(num_eig);
  std::vector<double> errors_bkwd(num_report), errors_abs(num_report);
  // Stage the converged eigenvectors in host memory so that the eigensolver, linear solver,
  // preconditioner and matrices release their memory before postprocessing and error
  // estimation, which then do not add to the peak memory of the solve.
  std::vector<std::vector<std::complex<double>>> eigenvectors(
      num_eig, std::vector<std::complex<double>>(Curl.Width()));
  ComplexVector E(Curl.Width());
  E.UseDevice(true);
  for (int i = 0; i < num_eig; i++)
  {
    eigenvalues[i] = eigen->GetEigenvalue(i);
    prev_shift_eigenvalues.push_back((!has_C && !has_A2) ? 1i * std::sqrt(eigenvalues[i])
                                                         : eigenvalues[i]);
    if (i < num_report)
    {
      errors_bkwd[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
      errors_abs[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    }
    eigen->GetEigenvector(i, E);
    E.Get(eigenvectors[i].data(), E.Size(), false);
  }

  eigen.reset();
  KM.reset();
  KeepCoarseChoice();
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
  Mpi::Print("\n");

  if (store_eigenpairs)
  {
    prev_eigenvectors = std::make_unique<VectorFieldTransfer>(
        space_op.GetNDSpace().GetParMesh(), iodata.solver.order);
    prev_eigenvalues = eigenvalues;
    prev_num_its = num_its;
  }
  auto StoreEigenvector = [&](const ComplexVector &x)
  {
    if (store_eigenpairs)
    {
      prev_eigenvectors->AddField(space_op.GetNDSpace().Get(), x.Real());
      prev_eigenvectors->AddField(space_op.GetNDSpace().Get(), x.Imag());
    }
  };

  // Error estimator, constructed once the solver objects have been released.
  if (ComputeErrorEstimation())
  {
    estimator.emplace(space_op, iodata);
  }

  // Eigenvector phase independent of the basis, the element order and the partition, and
  // the number of eigenvectors which took it from each reference (the fallbacks warn).
  const PhaseReference phase_ref(space_op.GetNDSpace());
  std::array<int, PhaseReference::num_fields + 1> num_ref = {};

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
    E.Set(eigenvectors[i].data(), E.Size(), false);
    if (const int ref = phase_ref.Normalize(E); ref >= 0)
    {
      num_ref[ref]++;
    }
    StoreEigenvector(E);

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
    if (estimator && i < num_requested)
    {
      estimator->AddErrorIndicator(E, B, total_domain_energy, indicator);
    }
  }
  for (int i = num_report; i < num_eig; i++)
  {
    E.Set(eigenvectors[i].data(), E.Size(), false);
    StoreEigenvector(E);
  }
  if (const int num_fallback = num_ref[1] + num_ref[PhaseReference::num_fields];
      num_fallback > 0)
  {
    Mpi::Warning(space_op.GetComm(),
                 "The projection on the first phase reference field is negligible for {:d} "
                 "of {:d} eigenvectors: {:d} took the phase of the second reference field "
                 "and {:d} that of the mean of the degrees of freedom, which depends on "
                 "the numbering of the mesh and on the number of processes.",
                 num_fallback, num_report, num_ref[1], num_ref[PhaseReference::num_fields]);
  }
  post_op.MeasureFinalize(indicator);
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
