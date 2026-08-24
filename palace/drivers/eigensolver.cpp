// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
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

internal::ResponseCorrectedMass
internal::BuildResponseCorrectedMass(const ComplexOperator &mass, const Operator &response)
{
  MFEM_VERIFY(mass.Real(), "Response-corrected eigenmode mass requires a real part!");
  ResponseCorrectedMass corrected;
  auto corrected_real = std::make_unique<SumOperator>(*mass.Real());
  corrected_real->AddOperator(response);
  corrected.real = ApplyEssentialDiagonal(std::move(corrected_real), *mass.Real());
  corrected.op =
      std::make_unique<ComplexWrapperOperator>(corrected.real.get(), mass.Imag());
  return corrected;
}

internal::EigenvalueTarget
internal::GetEigenvalueTarget(double target, bool lambda_eigenproblem,
                              EigenSolverBackend backend,
                              NonlinearEigenSolver nonlinear_type)
{
  if (lambda_eigenproblem)
  {
    EigenvalueSolver::WhichType which;
    if (backend == EigenSolverBackend::ARPACK)
    {
      // ARPACK searches the transformed spectrum. An eigenvalue just above iσ maps to a
      // large-magnitude negative imaginary value under shift-and-invert.
      which = EigenvalueSolver::WhichType::SMALLEST_IMAGINARY;
    }
    else if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      which = EigenvalueSolver::WhichType::TARGET_MAGNITUDE;
    }
    else
    {
      which = EigenvalueSolver::WhichType::TARGET_IMAGINARY;
    }
    return {1i * target, which};
  }

  // For μ = ω², an eigenvalue just below σ² maps to a large positive real value under
  // shift-and-invert.
  const auto which = backend == EigenSolverBackend::ARPACK
                         ? EigenvalueSolver::WhichType::LARGEST_REAL
                         : EigenvalueSolver::WhichType::TARGET_REAL;
  return {target * target, which};
}

std::complex<double> internal::EigenvalueToAngularFrequency(std::complex<double> value,
                                                            bool lambda_eigenproblem)
{
  return lambda_eigenproblem ? value / 1i : std::sqrt(value);
}

std::vector<int>
internal::MaximumWeightModeAssignment(const std::vector<std::vector<double>> &weights)
{
  if (weights.empty())
  {
    return {};
  }
  const int rows = static_cast<int>(weights.size());
  const int columns = static_cast<int>(weights.front().size());
  MFEM_VERIFY(columns >= rows,
              "Mode assignment requires at least as many columns as rows!");
  double maximum_weight = 0.0;
  for (const auto &row : weights)
  {
    MFEM_VERIFY(static_cast<int>(row.size()) == columns,
                "Mode-assignment weights must form a rectangular matrix!");
    for (const double weight : row)
    {
      MFEM_VERIFY(std::isfinite(weight) && weight >= 0.0,
                  "Mode-assignment weights must be finite and nonnegative!");
      maximum_weight = std::max(maximum_weight, weight);
    }
  }

  // Hungarian algorithm for a rectangular minimization problem. Subtracting each
  // overlap from the common maximum converts maximum-weight mode matching to minimum
  // cost without changing the optimal assignment.
  std::vector<double> row_potential(rows + 1), column_potential(columns + 1);
  std::vector<int> column_row(columns + 1), predecessor(columns + 1);
  for (int row = 1; row <= rows; row++)
  {
    column_row[0] = row;
    int current_column = 0;
    std::vector<double> minimum_cost(columns + 1, std::numeric_limits<double>::infinity());
    std::vector<bool> used(columns + 1, false);
    do
    {
      used[current_column] = true;
      const int current_row = column_row[current_column];
      double delta = std::numeric_limits<double>::infinity();
      int next_column = 0;
      for (int column = 1; column <= columns; column++)
      {
        if (used[column])
        {
          continue;
        }
        const double cost = maximum_weight - weights[current_row - 1][column - 1] -
                            row_potential[current_row] - column_potential[column];
        if (cost < minimum_cost[column])
        {
          minimum_cost[column] = cost;
          predecessor[column] = current_column;
        }
        if (minimum_cost[column] < delta)
        {
          delta = minimum_cost[column];
          next_column = column;
        }
      }
      for (int column = 0; column <= columns; column++)
      {
        if (used[column])
        {
          row_potential[column_row[column]] += delta;
          column_potential[column] -= delta;
        }
        else
        {
          minimum_cost[column] -= delta;
        }
      }
      current_column = next_column;
    } while (column_row[current_column] != 0);

    do
    {
      const int next_column = predecessor[current_column];
      column_row[current_column] = column_row[next_column];
      current_column = next_column;
    } while (current_column != 0);
  }

  std::vector<int> assignment(rows, -1);
  for (int column = 1; column <= columns; column++)
  {
    if (column_row[column] > 0)
    {
      assignment[column_row[column] - 1] = column - 1;
    }
  }
  MFEM_ASSERT(std::all_of(assignment.begin(), assignment.end(),
                          [](int column) { return column >= 0; }),
              "Incomplete maximum-weight mode assignment!");
  return assignment;
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
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

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
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op, &response_geometry,
                                               &surface_post_geometry);
  ComplexVector E(Curl.Width()), E_corrected, B(Curl.Height()),
      D(space_op.GetRTSpace().GetTrueVSize()), D_corrected;
  E.UseDevice(true);
  B.UseDevice(true);
  D.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
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
  auto MakeEigenvalueSolver = [&]() -> std::unique_ptr<EigenvalueSolver>
  {
    std::unique_ptr<EigenvalueSolver> result;
    if (type == EigenSolverBackend::ARPACK)
    {
#if defined(PALACE_WITH_ARPACK)
      if (C || has_A2)
      {
        result = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                           iodata.problem.verbose);
      }
      else
      {
        result = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                           iodata.problem.verbose);
      }
#endif
    }
    else  // EigenSolverBackend::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
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
                                  iodata.solver.linear.gs_orthog ==
                                      Orthogonalization::CGS2);
      result = std::move(slepc);
#endif
    }
    MFEM_VERIFY(result, "Requested eigenvalue-solver backend is unavailable!");
    return result;
  };

  Mpi::Print("\nConfiguring {} eigenvalue solver:\n",
             type == EigenSolverBackend::ARPACK ? "ARPACK" : "SLEPc");
  auto eigen = MakeEigenvalueSolver();
  SurfaceResponseOperator *response_correction = post_op.GetSurfaceResponseOperator();
  const bool self_consistent_response =
      response_correction &&
      iodata.solver.surface_response_correction->IncludesSelfConsistent();
  std::unique_ptr<ComplexWrapperOperator> response_mass;
  internal::ResponseCorrectedMass corrected_M, corrected_Mp;
  std::unique_ptr<EigenvalueSolver> corrected_eigen;
  if (self_consistent_response)
  {
    Mpi::Print("Configuring self-consistent response-corrected eigenvalue solver:\n");
    response_mass = std::make_unique<ComplexWrapperOperator>(response_correction, nullptr);
    corrected_M = internal::BuildResponseCorrectedMass(*M, *response_correction);
    corrected_eigen = MakeEigenvalueSolver();
    E_corrected.SetSize(Curl.Width());
    E_corrected.UseDevice(true);
    if (post_op.NeedsRecoveredElectricFlux())
    {
      D_corrected.SetSize(space_op.GetRTSpace().GetTrueVSize());
      D_corrected.UseDevice(true);
    }
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
  if (corrected_eigen)
  {
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      if (C)
      {
        corrected_eigen->SetOperators(*K, *C, *corrected_M.op,
                                      EigenvalueSolver::ScaleType::NONE);
      }
      else
      {
        corrected_eigen->SetOperators(*K, *corrected_M.op,
                                      EigenvalueSolver::ScaleType::NONE);
      }
      corrected_eigen->SetExtraSystemMatrix(funcA2);
      corrected_eigen->SetPreconditionerUpdate(funcP);
    }
    else if (has_A2)
    {
      corrected_Mp = internal::BuildResponseCorrectedMass(*Mp, *response_correction);
      corrected_eigen->SetOperators(*Kp, *Cp, *corrected_Mp.op, scale);
    }
    else if (C)
    {
      corrected_eigen->SetOperators(*K, *C, *corrected_M.op, scale);
    }
    else
    {
      corrected_eigen->SetOperators(*K, *corrected_M.op, scale);
    }
  }
  const double tol = (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
                         ? iodata.solver.eigenmode.linear_tol
                         : iodata.solver.eigenmode.tol;
  auto ConfigureCommon = [&](EigenvalueSolver &solver)
  {
    solver.SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
    solver.SetTol(tol);
    solver.SetMaxIter(iodata.solver.eigenmode.max_it);
  };
  ConfigureCommon(*eigen);
  if (corrected_eigen)
  {
    ConfigureCommon(*corrected_eigen);
  }
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());
  if (corrected_eigen)
  {
    Mpi::Print(" Corrected scaling γ = {:.3e}, δ = {:.3e}\n",
               corrected_eigen->GetScalingGamma(), corrected_eigen->GetScalingDelta());
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
    if (corrected_eigen)
    {
      // The raw positive mass form gives both mode families a common overlap and
      // normalization metric. The response-corrected mass can contain a small
      // indefinite defect and is not available as an assembled ParOperator.
      corrected_eigen->SetBMat(*KM);
    }

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree, corrected_divfree;
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
    if (corrected_eigen)
    {
      Mpi::Print(" Configuring response-corrected divergence-free projection\n");
      corrected_divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
          space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
          space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
          iodata.solver.linear.divfree_max_it, divfree_verbose, response_correction);
      corrected_eigen->SetDivFreeProjector(*corrected_divfree);
    }
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
    ComplexVector corrected_v0;
    if (corrected_divfree)
    {
      corrected_v0 = v0;
      corrected_divfree->Mult(corrected_v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector
    if (corrected_eigen)
    {
      corrected_eigen->SetInitialSpace(corrected_divfree ? corrected_v0 : v0);
    }

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
  const bool lambda_eigenproblem =
      C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP;
  auto ConfigureSpectralTarget = [&](EigenvalueSolver &solver)
  {
    const auto config =
        internal::GetEigenvalueTarget(target, lambda_eigenproblem, type, nonlinear_type);
    solver.SetShiftInvert(config.shift);
    solver.SetWhichEigenpairs(config.which);
  };
  ConfigureSpectralTarget(*eigen);
  if (corrected_eigen)
  {
    ConfigureSpectralTarget(*corrected_eigen);
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
  std::unique_ptr<ComplexOperator> corrected_A;
  std::unique_ptr<ComplexKspSolver> corrected_ksp;
  if (corrected_eigen)
  {
    corrected_A = BuildComplexSumOperator(
        {{1.0 + 0.0i, A.get()}, {-target * target + 0.0i, response_mass.get()}}, *K);
    corrected_ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                       &space_op.GetH1Spaces());
    corrected_ksp->SetOperators(*corrected_A, *P);
    corrected_eigen->SetLinearSolver(*corrected_ksp);
  }

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
  auto AddEstimate = [&](const ComplexVector &E, const ComplexVector &B,
                         const ComplexVector *D, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
    {
      if (D)
        estimator_2d->AddErrorIndicator(E, B, *D, Et, ind);
      else
        estimator_2d->AddErrorIndicator(E, B, Et, ind);
    }
    else
    {
      if (D)
        estimator_3d->AddErrorIndicator(E, B, *D, Et, ind);
      else
        estimator_3d->AddErrorIndicator(E, B, Et, ind);
    }
  };
  auto RecoverElectricFlux = [&](const ComplexVector &E, ComplexVector &D)
  {
    if (is_2d)
      estimator_2d->RecoverElectricFlux(E, D);
    else
      estimator_3d->RecoverElectricFlux(E, D);
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
  int corrected_num_conv = 0;
  if (corrected_eigen)
  {
    Mpi::Print("\nSolving self-consistent fabrication-response corrected eigenproblem\n");
    corrected_num_conv = corrected_eigen->Solve();
    std::complex<double> lambda =
        corrected_num_conv > 0 ? corrected_eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged corrected eigenvalue{}{}\n", corrected_num_conv,
               corrected_num_conv > 1 ? "s" : "",
               corrected_num_conv > 0
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID)
  {
    auto RefineEigenproblem =
        [&](std::unique_ptr<EigenvalueSolver> linear_eigen, int linear_num_conv,
            const ComplexOperator &mass, ComplexKspSolver &linear_solver,
            bool corrected) -> std::pair<std::unique_ptr<EigenvalueSolver>, int>
    {
      Mpi::Print("\n Refining {}eigenvalues with Quasi-Newton solver\n",
                 corrected ? "response-corrected " : "");
      auto qn = std::make_unique<QuasiNewtonSolver>(
          space_op.GetComm(), std::move(linear_eigen), linear_num_conv,
          iodata.problem.verbose, iodata.solver.eigenmode.refine_nonlinear);
      qn->SetTol(iodata.solver.eigenmode.tol);
      qn->SetMaxIter(iodata.solver.eigenmode.max_it);
      if (C)
      {
        qn->SetOperators(*K, *C, mass, EigenvalueSolver::ScaleType::NONE);
      }
      else
      {
        qn->SetOperators(*K, mass, EigenvalueSolver::ScaleType::NONE);
      }
      qn->SetExtraSystemMatrix(funcA2);
      qn->SetPreconditionerUpdate(funcP);
      qn->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
      qn->SetPreconditionerLag(iodata.solver.eigenmode.preconditioner_lag,
                               iodata.solver.eigenmode.preconditioner_lag_tol);
      qn->SetMaxRestart(iodata.solver.eigenmode.max_restart);
      qn->SetLinearSolver(linear_solver);
      qn->SetShiftInvert(1i * target);

      // Suppress wave port output during nonlinear eigensolver iterations.
      space_op.GetWavePortOp().SetSuppressOutput(true);
      const int refined_num_conv = qn->Solve();
      space_op.GetWavePortOp().SetSuppressOutput(false);
      return {std::move(qn), refined_num_conv};
    };

    std::tie(eigen, num_conv) =
        RefineEigenproblem(std::move(eigen), num_conv, *M, *ksp, false);
    if (corrected_eigen)
    {
      std::tie(corrected_eigen, corrected_num_conv) =
          RefineEigenproblem(std::move(corrected_eigen), corrected_num_conv,
                             *corrected_M.op, *corrected_ksp, true);
    }
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
    if (corrected_eigen)
    {
      corrected_eigen->SetBMat(*KM);
      corrected_eigen->RescaleEigenvectors(corrected_num_conv);
    }
  }

  std::vector<int> corrected_mode(num_conv, -1);
  std::vector<double> corrected_mode_overlap(num_conv,
                                             std::numeric_limits<double>::quiet_NaN());
  std::vector<std::complex<double>> corrected_omega(corrected_num_conv);
  auto GetOmega = [&](const EigenvalueSolver &solver, int i)
  {
    return internal::EigenvalueToAngularFrequency(solver.GetEigenvalue(i),
                                                  lambda_eigenproblem);
  };
  if (corrected_eigen)
  {
    MFEM_VERIFY(corrected_num_conv >= iodata.solver.eigenmode.n,
                "Response-corrected eigenmode solve found "
                    << corrected_num_conv << " modes when " << iodata.solver.eigenmode.n
                    << " were requested!");
    for (int j = 0; j < corrected_num_conv; j++)
    {
      corrected_omega[j] = GetOmega(*corrected_eigen, j);
    }

    const int pair_count = std::min(num_conv, corrected_num_conv);
    std::vector<std::vector<double>> overlap(pair_count,
                                             std::vector<double>(corrected_num_conv));
    ComplexVector raw_mode(Curl.Width()), response_mode(Curl.Width());
    raw_mode.UseDevice(true);
    response_mode.UseDevice(true);
    std::vector<double> raw_norm(pair_count), corrected_norm(corrected_num_conv);
    for (int i = 0; i < pair_count; i++)
    {
      eigen->GetEigenvector(i, raw_mode);
      raw_norm[i] =
          std::sqrt(std::abs(linalg::Dot(space_op.GetComm(), raw_mode, *KM, raw_mode)));
    }
    for (int j = 0; j < corrected_num_conv; j++)
    {
      corrected_eigen->GetEigenvector(j, response_mode);
      corrected_norm[j] = std::sqrt(
          std::abs(linalg::Dot(space_op.GetComm(), response_mode, *KM, response_mode)));
    }
    for (int i = 0; i < pair_count; i++)
    {
      eigen->GetEigenvector(i, raw_mode);
      for (int j = 0; j < corrected_num_conv; j++)
      {
        corrected_eigen->GetEigenvector(j, response_mode);
        const double scale = raw_norm[i] * corrected_norm[j];
        MFEM_VERIFY(scale > 0.0, "Cannot match a zero-norm eigenmode!");
        overlap[i][j] = std::min(
            1.0, std::abs(linalg::Dot(space_op.GetComm(), raw_mode, *KM, response_mode)) /
                     scale);
      }
    }
    const auto assignment = internal::MaximumWeightModeAssignment(overlap);
    for (int i = 0; i < pair_count; i++)
    {
      corrected_mode[i] = assignment[i];
      corrected_mode_overlap[i] = overlap[i][assignment[i]];
    }

    if (Mpi::Root(space_op.GetComm()))
    {
      Table table;
      table.col_options = {6, 6};
      table.insert(Column("raw", "Raw", 4, {}, {}, ""));
      table.insert(Column("corrected", "Corrected", 9, {}, {}, ""));
      table.insert(Column("overlap", "M-overlap"));
      table.insert(Column("f_raw", "f raw (GHz)"));
      table.insert(Column("f_corrected", "f corrected (GHz)"));
      table["raw"].print_as_int = true;
      table["corrected"].print_as_int = true;
      for (int i = 0; i < pair_count; i++)
      {
        const auto raw_omega = GetOmega(*eigen, i);
        const auto response_omega = corrected_omega[corrected_mode[i]];
        table["raw"] << i + 1;
        table["corrected"] << corrected_mode[i] + 1;
        table["overlap"] << corrected_mode_overlap[i];
        table["f_raw"] << iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(
                              raw_omega.real()) /
                              (2 * M_PI);
        table["f_corrected"] << iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(
                                    response_omega.real()) /
                                    (2 * M_PI);
      }
      Mpi::Print("\nSelf-consistent response-corrected mode matching:\n{}",
                 table.format_table());
    }
  }
  Mpi::Print("\n");

  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = GetOmega(*eigen, i);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

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

    const ComplexVector *recovered_flux = nullptr;
    if (post_op.NeedsRecoveredElectricFlux())
    {
      Mpi::Print(" Recovering electric flux for interface postprocessing\n");
      RecoverElectricFlux(E, D);
      post_op.SetRecoveredElectricFlux(D);
      recovered_flux = &D;
    }
    if (corrected_mode[i] >= 0)
    {
      corrected_eigen->GetEigenvector(corrected_mode[i], E_corrected);
      linalg::NormalizePhase(space_op.GetComm(), E_corrected);
      const ComplexVector *corrected_flux = nullptr;
      if (post_op.NeedsRecoveredElectricFlux())
      {
        Mpi::Print(" Recovering electric flux for self-consistent corrected mode\n");
        RecoverElectricFlux(E_corrected, D_corrected);
        corrected_flux = &D_corrected;
      }
      post_op.SetSurfaceResponseCorrectedField(E_corrected, corrected_flux,
                                               corrected_omega[corrected_mode[i]],
                                               corrected_mode_overlap[i]);
    }
    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_conv);

    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      AddEstimate(E, B, recovered_flux, total_domain_energy, indicator);
    }

    // Final write: Different condition than end of loop (i = num_conv - 1).
    if (i == iodata.solver.eigenmode.n - 1)
    {
      post_op.MeasureFinalize(indicator);
    }
  }
  if (self_consistent_response)
  {
    MFEM_ASSERT(corrected_ksp, "Missing corrected eigenmode linear solver!");
    SaveSurfaceResponseSolverMetadata(space_op.GetComm(), "Eigenmode",
                                      corrected_ksp->NumTotalMult(),
                                      corrected_ksp->NumTotalMultIterations());
  }
  if (response_correction)
  {
    SaveMetadata(*response_correction);
  }
  MFEM_VERIFY(num_conv >= iodata.solver.eigenmode.n, "Eigenmode solve only found "
                                                         << num_conv << " modes when "
                                                         << iodata.solver.eigenmode.n
                                                         << " were requested!");
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
