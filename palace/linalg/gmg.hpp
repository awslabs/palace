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
class ParFiniteElementSpaceHierarchy;

}  // namespace mfem

namespace palace
{

//
// Geometric multigrid preconditioner using a given coarse solver for the provided
// hierarchy of finite element spaces. Optionally can be configured to use auxiliary space
// smoothing at each level.
//
template <typename OperType>
class GeometricMultigridSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of V-cycles per preconditioner application.
  const int pc_it;

  // System matrices at each multigrid level and prolongation operators (not owned).
  std::vector<const OperType *> A;
  std::vector<const Operator *> P;
  std::vector<const mfem::Array<int> *> dbc_tdof_lists;

  // Smoothers for each level. Coarse level solver is B[0].
  mutable std::vector<std::unique_ptr<Solver<OperType>>> B;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<VecType> X, Y, R;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           const mfem::ParFiniteElementSpaceHierarchy *aux_fespaces,
                           int cycle_it, int smooth_it, int cheby_order,
                           double cheby_sf_max, double cheby_sf_min, bool cheby_4th_kind,
                           int pa_order_threshold, bool pa_discrete_interp);
  GeometricMultigridSolver(const IoData &iodata,
                           std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           const mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
    : GeometricMultigridSolver(
          std::move(coarse_solver), fespaces, aux_fespaces,
          iodata.solver.linear.mg_cycle_it, iodata.solver.linear.mg_smooth_it,
          iodata.solver.linear.mg_smooth_order, iodata.solver.linear.mg_smooth_sf_max,
          iodata.solver.linear.mg_smooth_sf_min, iodata.solver.linear.mg_smooth_cheby_4th,
          iodata.solver.pa_order_threshold, iodata.solver.pa_discrete_interp)
  {
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
