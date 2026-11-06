// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
#define PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP

#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/iodata.hpp"

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
// smoothing at each level.
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

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<VecType> X, Y, R;

  // Enable timer contribution for Timer::KSP_COARSE_SOLVE, unless the coarse solver times
  // itself (setup vs. apply).
  bool use_timer, coarse_timer;

  // Performance tally ids per level: pre-/post-smoothing, residual, restriction,
  // prolongation (index 0 holds the coarse solve).
  struct LevelTallies
  {
    int pre = -1, post = -1, res = -1, restr = -1, prol = -1;
  };
  std::vector<LevelTallies> perf;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G, int cycle_it,
                           int smooth_it, int cheby_order, double cheby_sf_max,
                           double cheby_sf_min, bool cheby_4th_kind);
  GeometricMultigridSolver(const IoData &iodata, MPI_Comm comm,
                           std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G = nullptr)
    : GeometricMultigridSolver(
          comm, std::move(coarse_solver), P, G, iodata.solver.linear.mg_cycle_it,
          iodata.solver.linear.mg_smooth_it, iodata.solver.linear.mg_smooth_order,
          iodata.solver.linear.mg_smooth_sf_max, iodata.solver.linear.mg_smooth_sf_min,
          iodata.solver.linear.mg_smooth_cheby_4th)
  {
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void SetPerfTallies(const std::string &prefix, Timer::Index parent) override;

  void EnableTimer();
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
