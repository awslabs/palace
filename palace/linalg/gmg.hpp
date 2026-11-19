// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
#define PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP

#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

//
// Geometric multigrid preconditioner using a given coarse solver for the provided
// hierarchy of finite element spaces. Optionally can be configured to use auxiliary space
// smoothing at each level, and one-sided V-cycles.
//
template <OperatorType OperType>
class GeometricMultigridSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of V-cycles per preconditioner application.
  const int pc_it;

  // Prolongation operators (not owned).
  std::vector<const Operator *> P;

  // System matrices at each multigrid level (not owned).
  std::vector<const OperType *> A;
  std::vector<const mfem::Array<int> *> dbc_tdof_lists;

  // Smoothers for each level. Coarse-level solver is B[0].
  mutable std::vector<std::unique_ptr<Solver<OperType>>> B;

  // EXPERIMENT (PR 1233 review): post-smoothers of another polynomial order than the
  // pre-smoothers (null where they are the same smoother).
  mutable std::vector<std::unique_ptr<Solver<OperType>>> B_post;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<VecType> X, Y, R;

  // Enable timer contribution for Timer::KSP_COARSE_SOLVE.
  bool use_timer;

  // One-sided V-cycles, which pre-smooth but do not post-smooth (Phillips and Fischer,
  // 2022): the preconditioner is not symmetric, so only for GMRES and FGMRES.
  const bool one_sided;

  // EXPERIMENT (PR 1233 review): polynomial order of the post-smoother (0 for none, as in
  // one-sided V-cycles, which then apply the auxiliary space correction after the
  // coarse-grid correction).
  int post_order;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G, int cycle_it,
                           int smooth_it, int cheby_order, double cheby_sf_max,
                           double cheby_sf_min, bool cheby_4th_kind, bool one_sided = false,
                           int cheby_post_order = -1, int cheby_coarse_order = -1);

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void EnableTimer() { use_timer = true; }
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
