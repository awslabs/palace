// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP

#include <filesystem>
#include <map>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

// Called before creating any output; the experimental reader must never reuse output.
void ValidateArchiveEstimateOptions(const IoData &iodata, MPI_Comm comm,
                                    bool check_output = true);

class ErrorIndicator;
class LaplaceOperator;
class Mesh;
class SurfacePostGeometry;
class SurfaceResponseGeometry;
template <ProblemType>
class PostOperator;
template <typename OperType>
class BaseKspSolver;
using KspSolver = BaseKspSolver<Operator>;

// The self-consistent response-corrected solve: PCG on K + Pᵀ (Q_fab - Q_thin) P with the
// thin-metal solver's preconditioner, at the correction's own tolerance, from the raw field
// as the initial guess. The corrected field is ACCEPTED only when PCG converged AND met no
// search direction of non-positive curvature (Ap, p) <= 0: a direction of negative curvature
// proves the corrected operator indefinite, and CG then "converges" to a meaningless field
// (sc-closure diagnostics 2026-09-29: the 120-degree corner model at p5 gave a corrected
// domain correction 300x the fixed-trace one with five such directions). The record carries
// the CG diagnostics; ritz_min / ritz_max are the extreme Ritz values of the preconditioned
// corrected operator from the CG coefficients (a negative ritz_min is the same proof); they
// are the eigenvalues of the CG Lanczos tridiagonal of the m iterations (typically 3-30) and
// are NaN (not evaluated) when m exceeds corrected_solve_ritz_max_iterations.
// The curvature record exists for CgSolver only: with another Krylov solver (GMRES, FGMRES)
// the record stays empty (iterations 0, count 0, Ritz values NaN), the curvature check is a
// no-op and acceptance reduces to convergence, as before this check.
constexpr int corrected_solve_ritz_max_iterations = 1000;
struct CorrectedSolveRecord
{
  bool converged = false;
  int iterations = 0;
  double relative_residual = 0.0;
  int negative_curvature_count = 0;
  int first_negative_curvature_iteration = -1;
  double ritz_min = 0.0;
  double ritz_max = 0.0;
  bool accepted = false;
};
CorrectedSolveRecord SolveCorrectedField(KspSolver &ksp, const Operator &K,
                                         const Operator &corrected_K, double solve_tol,
                                         bool initial_guess, const Vector &rhs, Vector &x);

//
// Driver class for electrostatic simulations.
//
class ElectrostaticSolver : public BaseSolver
{
private:
  mutable std::shared_ptr<const SurfacePostGeometry> surface_post_geometry;
  mutable std::shared_ptr<const SurfaceResponseGeometry> response_geometry;

  void PostprocessTerminals(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                            const std::map<int, mfem::Array<int>> &terminal_sources,
                            const std::vector<Vector> &V) const;
  void PostprocessResponseMatrix(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                                 const LaplaceOperator &laplace_op, const Operator &Grad,
                                 const std::vector<Vector> &V,
                                 const std::vector<Vector> &D) const;
  void PostprocessArchivedResponseMatrix(PostOperator<ProblemType::ELECTROSTATIC> &post_op,
                                         const LaplaceOperator &laplace_op,
                                         const Operator &Grad,
                                         const std::filesystem::path &archive,
                                         int block_size) const;

  ErrorIndicator EstimateArchivedFields(LaplaceOperator &laplace_op, const Operator &K,
                                        const std::filesystem::path &archive) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ELECTROSTATIC_SOLVER_HPP
