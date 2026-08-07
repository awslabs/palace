// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <algorithm>
#include <complex>
#include <map>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "fem/singularfield.hpp"
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
#include "models/hierarchicalmaxwellestimator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/singularsurfacepostoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

void EigenSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                             MPI_Comm comm) const
{
  BaseSolver::Preprocess(iodata, smesh, comm);
  singular_features.Preprocess(iodata, smesh, comm);
}

bool EigenSolver::RequiresSourceSerialMeshMetadata() const
{
  return iodata.solver.singular_elements.Enabled();
}

void EigenSolver::ProcessPartitionedMesh(const mfem::ParMesh &parallel_mesh,
                                         const mesh::PartitionMetadata &metadata) const
{
  singular_features.ProcessPartitionedMesh(iodata, parallel_mesh, metadata);
}

mesh::PartitionMetadata EigenSolver::GetSourceEntityMetadata() const
{
  return singular_features.GetSourceEntityMetadata();
}

mfem::Array<int> EigenSolver::GetRefinementProtection(const mfem::ParMesh &mesh,
                                                      bool *conforming,
                                                      mfem::Array<int> *repair) const
{
  if (!iodata.solver.singular_elements.Enabled())
  {
    if (conforming)
    {
      *conforming = true;
    }
    if (repair)
    {
      repair->SetSize(0);
    }
    return {};
  }
  return singular_features.GetRefinementProtection(mesh, conforming, repair);
}

mfem::Array<int> EigenSolver::GetEnrichedElements(const mfem::ParMesh &mesh) const
{
  if (!iodata.solver.singular_elements.Enabled())
  {
    return {};
  }
  return singular_features.GetEnrichedElements(mesh);
}

void EigenSolver::ReportTraceComponents(const mfem::ParMesh &mesh,
                                        const mfem::Array<int> &primary_marks) const
{
  if (iodata.solver.singular_elements.Enabled())
  {
    singular_features.ReportTraceComponents(iodata, mesh, primary_marks);
  }
}

void EigenSolver::ObserveRefinementAncestry(const mfem::ParMesh &mesh) const
{
  if (iodata.solver.singular_elements.Enabled())
  {
    singular_features.ObserveRefinementAncestry(mesh);
  }
}

void EigenSolver::ProcessRefinedMesh(const mfem::ParMesh &mesh) const
{
  if (iodata.solver.singular_elements.Enabled())
  {
    singular_features.ProcessRefinedMesh(iodata, mesh);
  }
}

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh, singular_features.GetSheetFeatures(),
                         singular_features.GetLineFeatures(),
                         singular_features.GetSourceVertexIds());
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  if (space_op.HasSingularEnrichment())
  {
    return SolveSingular(space_op, std::move(K), std::move(C), std::move(M));
  }

  // Check if there are nonlinear terms and, if so, setup interpolation operator.
  auto funcA2 = [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
  };
  auto funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                           std::complex<double> a2,
                           std::complex<double> a3) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, a3); };
  const double target = iodata.solver.eigenmode.target;
  auto A2 = funcA2(1i * target);
  bool has_A2 = (A2 != nullptr);

  // Extend K, C, M operators with interpolated A2 operator.
  // K' = K + A2_0, C' = C + A2_1, M' = M + A2_2
  std::unique_ptr<ComplexOperator> Kp, Cp, Mp;
  std::unique_ptr<Interpolation> interp_op;
  std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2;
  NonlinearEigenSolver nonlinear_type = iodata.solver.eigenmode.nonlinear_type;
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
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = iodata.solver.eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
#if !defined(PALACE_WITH_SLEPC)
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    Mpi::Warning("SLP nonlinear eigensolver not available without SLEPc, using Hybrid!\n");
  }
  nonlinear_type = NonlinearEigenSolver::HYBRID;
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
    eigen->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    eigen->SetExtraSystemMatrix(funcA2);
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
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) / (2 * M_PI);
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
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i,
                                    K.get(), C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
    {
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    }
    else
    {
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
    }
  };
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
    qn->SetMaxIter(iodata.solver.eigenmode.max_it);
    if (C)
    {
      qn->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      qn->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    qn->SetExtraSystemMatrix(funcA2);
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
  if (!KM)
  {
    // Normalize the finalized eigenvectors with respect to mass matrix (unit electric field
    // energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_conv);
  }
  Mpi::Print("\n");

  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (!C && !has_A2)
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
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_conv);

    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      AddEstimate(E, B, total_domain_energy, indicator);
    }

    // Final write: Different condition than end of loop (i = num_conv - 1).
    if (i == iodata.solver.eigenmode.n - 1)
    {
      post_op.MeasureFinalize(indicator);
    }
  }
  MFEM_VERIFY(num_conv >= iodata.solver.eigenmode.n, "Eigenmode solve only found "
                                                         << num_conv << " modes when "
                                                         << iodata.solver.eigenmode.n
                                                         << " were requested!");
  return {indicator, space_op.GlobalTrueVSize()};
}

std::pair<ErrorIndicator, long long int>
EigenSolver::SolveSingular(SpaceOperator &space_op, std::unique_ptr<ComplexOperator> K,
                           std::unique_ptr<ComplexOperator> C,
                           std::unique_ptr<ComplexOperator> M) const
{
  MFEM_VERIFY(
      K && M && K->Real() && !K->Imag() && (!C || (C->Real() && !C->Imag())) && M->Real(),
      "Eigenmode singular simulations require real stiffness and damping operators and a "
      "complex electric mass operator!");
  const bool damped = (C != nullptr);
  auto K_energy = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto K_bulk = space_op.GetBulkStiffnessMatrix(Operator::DIAG_ZERO);
  auto M_energy = space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  auto M_bulk = space_op.GetBulkMassMatrix(Operator::DIAG_ZERO);
  const double target = iodata.solver.eigenmode.target;

  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = iodata.solver.eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK enriched eigenvalue solver:\n");
    if (damped)
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
  else
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc enriched eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (damped && !iodata.solver.eigenmode.pep_linear)
    {
      slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
    }
    else if (damped)
    {
      slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                            iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    }
    else
    {
      slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    }
    slepc->SetProblemType(damped || M->Imag()
                              ? slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN
                              : slepc::SlepcEigenvalueSolver::ProblemType::GEN_HERMITIAN);
    slepc->SetOrthogonalization(iodata.solver.linear.gs_orthog == Orthogonalization::MGS,
                                iodata.solver.linear.gs_orthog == Orthogonalization::CGS2);
    eigen = std::move(slepc);
#endif
  }
  MFEM_VERIFY(eigen, "Requested eigenvalue solver backend is unavailable!");

  const auto scale = iodata.solver.eigenmode.scale ? EigenvalueSolver::ScaleType::NORM_2
                                                   : EigenvalueSolver::ScaleType::NONE;
  if (damped)
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

  auto mass_inner_product = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    eigen->SetBMat(*mass_inner_product);
  }

  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  const bool has_unconstrained_boundary_mass = [&]()
  {
    std::set<int> constrained_attributes;
    const auto collect_attributes =
        [&constrained_attributes](const std::map<int, double> &coefficients)
    {
      for (const auto &[attribute, coefficient] : coefficients)
      {
        if (coefficient != 0.0)
        {
          constrained_attributes.insert(attribute);
        }
      }
    };
    collect_attributes(space_op.GetLumpedPortOp().GetStiffnessBdrCoefficientMap());
    collect_attributes(space_op.GetLumpedPortOp().GetDampingBdrCoefficientMap());
    collect_attributes(space_op.GetSurfaceImpedanceOp().GetStiffnessBdrCoefficientMap());
    collect_attributes(space_op.GetSurfaceImpedanceOp().GetDampingBdrCoefficientMap());

    const auto has_unconstrained_attribute =
        [&constrained_attributes](const std::map<int, double> &coefficients)
    {
      return std::any_of(coefficients.begin(), coefficients.end(),
                         [&constrained_attributes](const auto &entry)
                         {
                           return entry.second != 0.0 &&
                                  constrained_attributes.find(entry.first) ==
                                      constrained_attributes.end();
                         });
    };
    return has_unconstrained_attribute(
               space_op.GetLumpedPortOp().GetMassBdrCoefficientMap()) ||
           has_unconstrained_attribute(
               space_op.GetSurfaceImpedanceOp().GetMassBdrCoefficientMap());
  }();
  if (iodata.solver.linear.divfree_max_it > 0 && !has_unconstrained_boundary_mass)
  {
    Mpi::Print(" Configuring enriched divergence-free projection\n");
    auto scalar_diffusion = space_op.GetBulkScalarDiffusionMatrix();
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        iodata, space_op.GetComm(), *M_bulk, space_op.GetGradMatrix(),
        std::move(scalar_diffusion), space_op.GetCombinedH1AuxBdrTDofList());
    eigen->SetDivFreeProjector(*divfree);
  }
  else if (iodata.solver.linear.divfree_max_it > 0)
  {
    Mpi::Print(" Skipping enriched divergence-free projection because a purely capacitive "
               "boundary requires a bulk-plus-boundary mass metric\n");
  }

  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector initial;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      space_op.GetConstantInitialVector(initial);
    }
    else
    {
      space_op.GetRandomInitialVector(initial);
    }
    if (divfree)
    {
      divfree->Mult(initial);
    }
    eigen->SetInitialSpace(initial);
  }

  if (damped)
  {
    eigen->SetShiftInvert(1i * target);
    eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
  }
  else
  {
    eigen->SetShiftInvert(target * target);
    eigen->SetWhichEigenpairs(type == EigenSolverBackend::ARPACK
                                  ? EigenvalueSolver::WhichType::LARGEST_REAL
                                  : EigenvalueSolver::WhichType::TARGET_REAL);
  }
  const std::complex<double> damping_shift = damped ? 1i * target : 0.0;
  auto shifted = space_op.GetSystemMatrix(
      1.0 + 0.0i, damping_shift, -target * target + 0.0i, K.get(), C.get(), M.get());
  const auto combined_prolongations = space_op.GetCombinedNDProlongationOperators();
  const auto combined_gradients = space_op.GetCombinedGradientOperators();
  auto ksp = MakeSingularComplexKspSolver(
      iodata, space_op.GetNDSpaces(), space_op.GetH1Spaces(), combined_prolongations,
      combined_gradients, space_op.GetCombinedNDDbcTDofLists(),
      space_op.GetCombinedH1DbcTDofLists());
  auto preconditioner = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, damping_shift, -target * target + 0.0i, target);
  ksp->SetOperators(*shifted, *preconditioner);
  eigen->SetLinearSolver(*ksp);

  {
    const double frequency =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) / (2 * M_PI);
    Mpi::Print(" Shift-and-invert sigma = {:.3e} GHz ({:.3e}{})\n", frequency, target,
               damped ? "i" : "^2");
  }
  int num_converged;
  {
    BlockTimer bt(Timer::EPS);
    Mpi::Print("\n");
    num_converged = eigen->Solve();
    MFEM_VERIFY(ksp->NumFailedMult() == 0,
                "Singular eigenmode shift-and-invert failed to converge "
                    << ksp->NumFailedMult() << " of " << ksp->NumTotalMult()
                    << " inner linear solves!");
    Mpi::Print(" Found {:d} converged enriched eigenvalue{}\n", num_converged,
               num_converged == 1 ? "" : "s");
  }
  if (!iodata.solver.eigenmode.mass_orthog)
  {
    // Match the standard eigenmode path: use Euclidean orthogonalization during
    // the solve, but mass-normalize finalized eigenvectors for energy output.
    eigen->SetBMat(*mass_inner_product);
  }
  eigen->RescaleEigenvectors(num_converged);
  SaveMetadata(*ksp);
  BlockTimer postprocess_timer(Timer::POSTPRO);

  const fem::singular::AdaptiveAssemblyOptions surface_options{
      iodata.solver.singular_elements.quadrature_order,
      iodata.solver.singular_elements.abs_tol,
      iodata.solver.singular_elements.rel_tol,
      iodata.solver.singular_elements.max_subdivisions,
      iodata.solver.singular_elements.UsesFixedSubdivision(),
      iodata.solver.singular_elements.subdivisions};
  std::unique_ptr<fem::singular::TriangleEnrichedNDFieldEvaluator>
      triangular_real_evaluator, triangular_imaginary_evaluator;
  std::unique_ptr<TetrahedronSingularSurfacePostOperator> tetrahedral_surface_postoperator;
  std::unique_ptr<TriangleSingularSurfacePostOperator> triangular_surface_postoperator;
  if (space_op.GetMesh().Dimension() == 3)
  {
    const auto *topology = space_op.GetSingularDofTopology();
    MFEM_VERIFY(topology,
                "Three-dimensional eigenmode singular surface postprocessing requires "
                "tetrahedral singular DOF topology!");
    tetrahedral_surface_postoperator =
        std::make_unique<TetrahedronSingularSurfacePostOperator>(
            iodata.boundaries.postpro, space_op.GetMaterialOp(),
            space_op.GetNDSpace().Get());
  }
  else
  {
    const auto *topology = space_op.GetTriangleSingularDofTopology();
    MFEM_VERIFY(topology,
                "Two-dimensional eigenmode singular surface postprocessing requires "
                "triangular singular DOF topology!");
    triangular_real_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            *topology, space_op.GetSingularParallelNumbering(),
            space_op.GetNDSpace().Get());
    triangular_imaginary_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            *topology, space_op.GetSingularParallelNumbering(),
            space_op.GetNDSpace().Get());
    triangular_surface_postoperator = std::make_unique<TriangleSingularSurfacePostOperator>(
        iodata.boundaries.postpro, space_op.GetMaterialOp(), space_op.GetNDSpace().Get());
  }
  {
    const auto start = Timer::Now();
    space_op.CacheSingularLumpedPortFunctionals(false);
    double elapsed = Timer::Duration(Timer::Now() - start).count();
    Mpi::GlobalMax(1, &elapsed, space_op.GetComm());
    Mpi::Print(" Singular postprocessing setup, cached lumped-port functionals (s): "
               "{:.3f}\n",
               elapsed);
  }

  TableWithCSVFile eig_output, energy_output, surface_output;
  if (root)
  {
    eig_output = TableWithCSVFile(post_dir / "eig.csv");
    eig_output.table.insert("idx", "m", -1, 0, 0, "");
    eig_output.table.insert("f_re", "Re{f} (GHz)");
    eig_output.table.insert("f_im", "Im{f} (GHz)");
    eig_output.table.insert("q", "Q");
    eig_output.table.insert("err_back", "Error (Bkwd.)");
    eig_output.table.insert("err_abs", "Error (Abs.)");
    eig_output.table[0].print_as_int = true;
    eig_output.WriteFullTableTrunc();

    energy_output = TableWithCSVFile(post_dir / "singular-eigenmode.csv");
    energy_output.table.insert("idx", "m", -1, 0, 0, "");
    energy_output.table.insert("electric_energy", "Electric field energy (J)");
    energy_output.table.insert("magnetic_energy", "Magnetic field energy (J)");
    energy_output.table.insert("energy_mismatch", "Relative energy mismatch");
    energy_output.table.insert("weak_divergence", "Relative weak divergence");
    for (const auto &[port_index, port] : space_op.GetLumpedPortOp())
    {
      energy_output.table.insert(fmt::format("voltage_real_{}", port_index),
                                 fmt::format("Re{{V[{}]}} (V)", port_index));
      energy_output.table.insert(fmt::format("voltage_imag_{}", port_index),
                                 fmt::format("Im{{V[{}]}} (V)", port_index));
      energy_output.table.insert(fmt::format("current_real_{}", port_index),
                                 fmt::format("Re{{I[{}]}} (A)", port_index));
      energy_output.table.insert(fmt::format("current_imag_{}", port_index),
                                 fmt::format("Im{{I[{}]}} (A)", port_index));
      energy_output.table.insert(fmt::format("inductive_participation_{}", port_index),
                                 fmt::format("p_ind[{}]", port_index));
    }
    energy_output.table[0].print_as_int = true;
    energy_output.WriteFullTableTrunc();

    const bool has_surface = tetrahedral_surface_postoperator
                                 ? !tetrahedral_surface_postoperator->Empty()
                                 : !triangular_surface_postoperator->Empty();
    if (has_surface)
    {
      surface_output = TableWithCSVFile(post_dir / "surface-Q.csv");
      surface_output.table.insert("idx", "m", -1, 0, 0, "");
      surface_output.table[0].print_as_int = true;
      for (const auto &[index, data] : iodata.boundaries.postpro.dielectric)
      {
        surface_output.table.insert(fmt::format("p_{}", index),
                                    fmt::format("p_surf[{}]", index));
        surface_output.table.insert(fmt::format("Q_{}", index),
                                    fmt::format("Q_surf[{}]", index));
      }
      surface_output.WriteFullTableTrunc();
    }
  }

  ComplexVector electric_field(space_op.GetNDTrueVSize()),
      mass_electric(space_op.GetNDTrueVSize());
  ComplexVector weak_divergence(space_op.GetH1TrueVSize());
  electric_field.UseDevice(true);
  mass_electric.UseDevice(true);
  weak_divergence.UseDevice(true);
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  ErrorIndicator indicator;
  // Hierarchical p+1 residual estimation for three-dimensional singular AMR, opt-in while
  // the element-patch lifting is qualified. Measured on the wedge eigenmode qualification
  // problem, this lifting marks too little of the singular region: it reaches mode-1
  // frequency 1.0359e-1 GHz at 1413 elements after three theta = 0.5 passes, while the
  // sliced flux recovery reaches 1.0318e-1 at 4433 elements against a uniform-refinement
  // trend of about 1.0303e-1. The certified edge/face/interior patch engine must replace
  // the element-patch shapes in parallel before this becomes the default estimator.
  std::unique_ptr<HierarchicalMaxwellDomainData> hierarchy;
  if (!is_2d && iodata.model.refinement.max_it > 0 &&
      iodata.solver.singular_elements.hierarchical_estimator &&
      space_op.GetMesh().Get().Conforming())
  {
    hierarchy = std::make_unique<HierarchicalMaxwellDomainData>(space_op);
  }
  const auto &gradient = space_op.GetGradMatrix();
  std::vector<std::unique_ptr<fem::singular::EnrichedNDFieldEvaluator>>
      tetrahedral_real_evaluators, tetrahedral_imaginary_evaluators;
  std::vector<double> tetrahedral_participation_energies;
  std::vector<std::vector<TetrahedronSingularSurfacePostOperator::Measurement>>
      tetrahedral_surface_measurements;
  if (tetrahedral_surface_postoperator && !tetrahedral_surface_postoperator->Empty())
  {
    const auto *topology = space_op.GetSingularDofTopology();
    MFEM_VERIFY(topology,
                "Three-dimensional eigenmode singular surface postprocessing requires "
                "tetrahedral singular DOF topology!");
    tetrahedral_real_evaluators.reserve(num_converged);
    tetrahedral_imaginary_evaluators.reserve(num_converged);
    tetrahedral_participation_energies.reserve(num_converged);
    std::vector<TetrahedronSingularSurfacePostOperator::NDFieldEvaluatorPair>
        field_evaluators;
    field_evaluators.reserve(num_converged);
    const auto setup_start = Timer::Now();
    for (int mode = 0; mode < num_converged; mode++)
    {
      eigen->GetEigenvector(mode, electric_field);
      linalg::NormalizePhase(space_op.GetComm(), electric_field);
      const double electric_energy =
          0.5 * std::max(0.0, SingularComplexQuadraticForm(space_op.GetComm(), *M_bulk,
                                                           electric_field));
      double capacitor_energy = 0.0;
      for (const auto &[port_index, port] : space_op.GetLumpedPortOp())
      {
        const auto voltage =
            space_op.GetSingularLumpedPortVoltage(port_index, electric_field);
        capacitor_energy +=
            0.5 * std::abs(port.C) * std::real(voltage * std::conj(voltage));
      }
      const double participation_energy = electric_energy + capacitor_energy;
      MFEM_VERIFY(participation_energy > 0.0 && std::isfinite(participation_energy),
                  "Eigenmode singular participation requires positive finite bulk "
                  "electric plus lumped-capacitor energy!");
      tetrahedral_participation_energies.push_back(participation_energy);

      tetrahedral_real_evaluators.push_back(
          std::make_unique<fem::singular::EnrichedNDFieldEvaluator>(
              *topology, space_op.GetSingularParallelNumbering(),
              space_op.GetNDSpace().Get()));
      tetrahedral_imaginary_evaluators.push_back(
          std::make_unique<fem::singular::EnrichedNDFieldEvaluator>(
              *topology, space_op.GetSingularParallelNumbering(),
              space_op.GetNDSpace().Get()));
      auto &real_evaluator = *tetrahedral_real_evaluators.back();
      auto &imaginary_evaluator = *tetrahedral_imaginary_evaluators.back();
      real_evaluator.SetFromTrueDofs(electric_field.Real());
      imaginary_evaluator.SetFromTrueDofs(electric_field.Imag());
      field_evaluators.emplace_back(&real_evaluator, &imaginary_evaluator);
    }
    double setup_elapsed = Timer::Duration(Timer::Now() - setup_start).count();
    Mpi::GlobalMax(1, &setup_elapsed, space_op.GetComm());
    Mpi::Print(" Singular batched surface evaluator setup for {:d} mode{} (s): {:.3f}\n",
               num_converged, num_converged == 1 ? "" : "s", setup_elapsed);

    const auto integration_start = Timer::Now();
    tetrahedral_surface_measurements = tetrahedral_surface_postoperator->Measure(
        field_evaluators, tetrahedral_participation_energies, surface_options);
    double integration_elapsed = Timer::Duration(Timer::Now() - integration_start).count();
    Mpi::GlobalMax(1, &integration_elapsed, space_op.GetComm());
    Mpi::Print(" Singular batched surface integration for {:d} mode{} (s): {:.3f}\n",
               num_converged, num_converged == 1 ? "" : "s", integration_elapsed);
  }
  for (int mode = 0; mode < num_converged; mode++)
  {
    const auto mode_start = Timer::Now();
    auto stage_start = mode_start;
    const auto report_stage = [&](std::string_view stage)
    {
      const auto now = Timer::Now();
      double elapsed = Timer::Duration(now - stage_start).count();
      double minimum = elapsed;
      double maximum = elapsed;
      Mpi::GlobalMin(1, &minimum, space_op.GetComm());
      Mpi::GlobalMax(1, &maximum, space_op.GetComm());
      Mpi::GlobalSum(1, &elapsed, space_op.GetComm());
      Mpi::Print(" Singular mode {:d} timing, {} (s): min. {:.3f}, max. {:.3f}, "
                 "avg. {:.3f}\n",
                 mode + 1, stage, minimum, maximum,
                 elapsed / Mpi::Size(space_op.GetComm()));
      stage_start = Timer::Now();
    };
    const std::complex<double> eigenvalue = eigen->GetEigenvalue(mode);
    const std::complex<double> omega = damped ? eigenvalue / 1i : std::sqrt(eigenvalue);
    const double error_backward =
        eigen->GetError(mode, EigenvalueSolver::ErrorType::BACKWARD);
    const double error_absolute =
        eigen->GetError(mode, EigenvalueSolver::ErrorType::ABSOLUTE);
    eigen->GetEigenvector(mode, electric_field);
    linalg::NormalizePhase(space_op.GetComm(), electric_field);
    report_stage("eigenvector extraction");

    const double augmented_mass_form =
        SingularComplexQuadraticForm(space_op.GetComm(), *M_energy, electric_field);
    const double augmented_stiffness_form =
        SingularComplexQuadraticForm(space_op.GetComm(), *K_energy, electric_field) /
        std::norm(omega);
    const double bulk_mass_form =
        SingularComplexQuadraticForm(space_op.GetComm(), *M_bulk, electric_field);
    const double bulk_stiffness_form =
        SingularComplexQuadraticForm(space_op.GetComm(), *K_bulk, electric_field) /
        std::norm(omega);
    Mpi::Print(" Singular mode {:d} quadratic forms: M_aug = {:.16e}, "
               "K_aug/|omega|^2 = {:.16e}, M_bulk = {:.16e}, "
               "K_bulk/|omega|^2 = {:.16e}\n",
               mode + 1, augmented_mass_form, augmented_stiffness_form, bulk_mass_form,
               bulk_stiffness_form);

    const auto augmented_energy = MeasureSingularFullWaveEnergy(
        space_op.GetComm(), *M_energy, *K_energy, electric_field, omega);
    const auto field_energy = MeasureSingularFullWaveEnergy(space_op.GetComm(), *M_bulk,
                                                            *K_bulk, electric_field, omega);
    const double energy_mismatch =
        std::abs(augmented_energy.electric - augmented_energy.magnetic) /
        std::max({augmented_energy.electric, augmented_energy.magnetic,
                  std::numeric_limits<double>::min()});
    report_stage("quadratic energy products");
    struct LumpedPortMeasurement
    {
      int index;
      std::complex<double> voltage;
      std::complex<double> current;
      double inductive_participation;
    };
    std::vector<LumpedPortMeasurement> lumped_port_measurements;
    double capacitor_energy = 0.0;
    for (const auto &[port_index, port] : space_op.GetLumpedPortOp())
    {
      const auto voltage =
          space_op.GetSingularLumpedPortVoltage(port_index, electric_field);
      const auto current = voltage / port.GetCharacteristicImpedance(omega.real());
      const std::complex<double> inductive_current =
          std::abs(port.L) > 0.0
              ? voltage /
                    port.GetCharacteristicImpedance(omega.real(), LumpedPortData::Branch::L)
              : 0.0;
      const double inductive_energy =
          0.5 * std::abs(port.L) *
          std::real(inductive_current * std::conj(inductive_current));
      capacitor_energy += 0.5 * std::abs(port.C) * std::real(voltage * std::conj(voltage));
      lumped_port_measurements.push_back(
          {port_index, voltage, current,
           std::copysign(inductive_energy, inductive_current.real())});
    }
    const double participation_energy = field_energy.electric + capacitor_energy;
    MFEM_VERIFY(participation_energy > 0.0 && std::isfinite(participation_energy),
                "Eigenmode singular participation requires positive finite bulk electric "
                "plus lumped-capacitor energy!");
    for (auto &measurement : lumped_port_measurements)
    {
      measurement.inductive_participation /= participation_energy;
    }
    if (!tetrahedral_participation_energies.empty())
    {
      const double cached_participation_energy =
          tetrahedral_participation_energies.at(mode);
      const double scale =
          std::max({std::abs(participation_energy), std::abs(cached_participation_energy),
                    std::numeric_limits<double>::min()});
      MFEM_VERIFY(std::abs(participation_energy - cached_participation_energy) <=
                      1.0e-11 * scale,
                  "Batched eigenmode singular surface participation normalization "
                  "changed during output processing!");
    }
    report_stage("lumped-port evaluation");

    std::vector<TriangleSingularSurfacePostOperator::Measurement> surface_measurements;
    if (triangular_surface_postoperator && !triangular_surface_postoperator->Empty())
    {
      triangular_real_evaluator->SetFromTrueDofs(electric_field.Real());
      triangular_imaginary_evaluator->SetFromTrueDofs(electric_field.Imag());
    }
    report_stage("surface evaluator update");
    if (tetrahedral_surface_postoperator && !tetrahedral_surface_postoperator->Empty())
    {
      surface_measurements = tetrahedral_surface_measurements.at(mode);
    }
    else if (triangular_surface_postoperator && !triangular_surface_postoperator->Empty())
    {
      surface_measurements = triangular_surface_postoperator->Measure(
          *triangular_real_evaluator, *triangular_imaginary_evaluator, participation_energy,
          surface_options);
    }
    report_stage("surface integration");
    const double electric_energy =
        iodata.units.Dimensionalize<Units::ValueType::ENERGY>(field_energy.electric);
    const double magnetic_energy =
        iodata.units.Dimensionalize<Units::ValueType::ENERGY>(field_energy.magnetic);

    M_bulk->Mult(electric_field.Real(), mass_electric.Real());
    M_bulk->Mult(electric_field.Imag(), mass_electric.Imag());
    gradient.MultTranspose(mass_electric.Real(), weak_divergence.Real());
    gradient.MultTranspose(mass_electric.Imag(), weak_divergence.Imag());
    linalg::SetSubVector(weak_divergence, space_op.GetCombinedH1DbcTDofList(), 0.0);
    const double relative_weak_divergence =
        linalg::Norml2(space_op.GetComm(), weak_divergence) /
        std::max(linalg::Norml2(space_op.GetComm(), mass_electric),
                 std::numeric_limits<double>::min());

    const std::complex<double> frequency =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI);
    const double quality =
        omega == 0.0 ? mfem::infinity() : 0.5 * std::abs(omega) / std::abs(omega.imag());
    Mpi::Print(" Mode {:d}: f = {:.6e}{:+.6e}i GHz, Q = {:.6e}, "
               "energy mismatch = {:.3e}, weak divergence = {:.3e}\n",
               mode + 1, frequency.real(), frequency.imag(), quality, energy_mismatch,
               relative_weak_divergence);

    if (mode < iodata.solver.eigenmode.n)
    {
      const int standard_nd_size = space_op.GetNDSpace().GetTrueVSize();
      Vector electric_real, electric_imaginary;
      electric_real.MakeRef(electric_field.Real(), 0, standard_nd_size);
      electric_imaginary.MakeRef(electric_field.Imag(), 0, standard_nd_size);
      ComplexVector standard_electric(electric_real, electric_imaginary);
      ComplexVector standard_magnetic(space_op.GetCurlMatrix().Height());
      space_op.GetCurlMatrix().Mult(standard_electric.Real(), standard_magnetic.Real());
      space_op.GetCurlMatrix().Mult(standard_electric.Imag(), standard_magnetic.Imag());
      standard_magnetic *= std::complex<double>(0.0, 1.0) / omega;
      const double total_field_energy = field_energy.electric + field_energy.magnetic;
      if (is_2d)
      {
        estimator_2d->AddErrorIndicator(standard_electric, standard_magnetic,
                                        total_field_energy, indicator);
      }
      else if (hierarchy)
      {
        // The certified edge/face/interior lifting is used whenever available (currently
        // one rank); the element-patch shape is the MPI-capable interim.
        BlockTimer bt_estimate(Timer::ESTIMATION);
        const auto patch_shape = hierarchy->EntityPatchesAvailable()
                                     ? HierarchicalMaxwellDomainData::PatchShape::ENTITY
                                     : HierarchicalMaxwellDomainData::PatchShape::ELEMENT;
        const auto estimate =
            hierarchy->EstimatePolynomialEigenResidual(omega, electric_field, patch_shape);
        Vector local_indicator(estimate.indicator_energy.Size());
        const double scale = total_field_energy > 0.0 ? 0.5 / total_field_energy : 1.0;
        for (int element = 0; element < local_indicator.Size(); element++)
        {
          local_indicator(element) =
              std::sqrt(scale * std::max(estimate.indicator_energy(element), 0.0));
        }
        indicator.AddIndicator(local_indicator);
        ErrorIndicator sliced_diagnostic;
        estimator_3d->AddErrorIndicator(standard_electric, standard_magnetic,
                                        total_field_energy, sliced_diagnostic);
        Mpi::Print(" Mode {:d} hierarchical residual: energy = {:.3e}, indicator = "
                   "{:.3e} (sliced recovery {:.3e})\n",
                   mode + 1, estimate.total_energy,
                   std::sqrt(scale * std::max(estimate.total_energy, 0.0)),
                   sliced_diagnostic.Norml2(space_op.GetComm()));
      }
      else
      {
        estimator_3d->AddErrorIndicator(standard_electric, standard_magnetic,
                                        total_field_energy, indicator);
      }
    }
    report_stage("diagnostics and estimation");
    if (root)
    {
      eig_output.table["idx"] << mode + 1;
      eig_output.table["f_re"] << frequency.real();
      eig_output.table["f_im"] << frequency.imag();
      eig_output.table["q"] << quality;
      eig_output.table["err_back"] << error_backward;
      eig_output.table["err_abs"] << error_absolute;
      eig_output.WriteFullTableTrunc();

      energy_output.table["idx"] << mode + 1;
      energy_output.table["electric_energy"] << electric_energy;
      energy_output.table["magnetic_energy"] << magnetic_energy;
      energy_output.table["energy_mismatch"] << energy_mismatch;
      energy_output.table["weak_divergence"] << relative_weak_divergence;
      for (const auto &measurement : lumped_port_measurements)
      {
        const auto voltage =
            iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(measurement.voltage);
        const auto current =
            iodata.units.Dimensionalize<Units::ValueType::CURRENT>(measurement.current);
        energy_output.table[fmt::format("voltage_real_{}", measurement.index)]
            << voltage.real();
        energy_output.table[fmt::format("voltage_imag_{}", measurement.index)]
            << voltage.imag();
        energy_output.table[fmt::format("current_real_{}", measurement.index)]
            << current.real();
        energy_output.table[fmt::format("current_imag_{}", measurement.index)]
            << current.imag();
        energy_output.table[fmt::format("inductive_participation_{}", measurement.index)]
            << measurement.inductive_participation;
      }
      energy_output.WriteFullTableTrunc();
      if (!surface_measurements.empty())
      {
        surface_output.table["idx"] << mode + 1;
        for (const auto &measurement : surface_measurements)
        {
          surface_output.table[fmt::format("p_{}", measurement.index)]
              << measurement.participation;
          surface_output.table[fmt::format("Q_{}", measurement.index)]
              << measurement.quality_factor;
        }
        surface_output.WriteFullTableTrunc();
      }
    }
    report_stage("output");
    double mode_total = Timer::Duration(Timer::Now() - mode_start).count();
    Mpi::GlobalMax(1, &mode_total, space_op.GetComm());
    Mpi::Print(" Singular mode {:d} timing, total (s): {:.3f}\n", mode + 1, mode_total);
  }

  SaveMetadata(
      "SingularFullWave",
      nlohmann::json{
          {"Enabled", true},
          {"Dimension", space_op.GetMesh().Dimension()},
          {"EigenvalueOutput", "eig.csv"},
          {"EnergyOutput", "singular-eigenmode.csv"},
          {"ElectricEnergy", "0.5 E^H M_bulk E"},
          {"MagneticEnergy", "0.5 E^H K_bulk E / |omega|^2"},
          {"EnergyMismatch", "complete reactive-boundary-augmented eigenproblem"},
          {"DampedQuadraticEigenproblem", damped},
          {"SurfaceParticipationDenominator", "E_electric_bulk + E_capacitor"},
          {"BulkDielectricLoss", space_op.GetMaterialOp().HasLossTangent()},
          {"LumpedPorts", space_op.GetLumpedPortOp().Size()},
          {"DivergenceProjection", divfree != nullptr},
          {"DivergenceProjectionShiftInvert", divfree != nullptr},
          {"FieldGridOutput", false},
          {"ErrorEstimator", hierarchy
                                 ? "hierarchical p+1 polynomial residual lifting"
                                 : "standard-space smooth remainder on the complete mesh"},
          {"SurfaceIntegrability", space_op.GetMesh().Dimension() == 3
                                       ? GetSingularSurfaceIntegrabilityMetadata(
                                             *singular_features.GetSheetFeatures())
                                       : GetSingularSurfaceIntegrabilityMetadata(
                                             *singular_features.GetLineFeatures())},
          {"SurfaceParticipation", GetSingularSurfaceParticipationMetadata(iodata)}});
  MFEM_VERIFY(num_converged >= iodata.solver.eigenmode.n,
              "Singular eigenmode solve found " << num_converged << " modes when "
                                                << iodata.solver.eigenmode.n
                                                << " were requested!");
  return {std::move(indicator), space_op.GlobalTrueVSize()};
}

}  // namespace palace
