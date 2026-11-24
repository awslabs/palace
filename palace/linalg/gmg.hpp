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

class FiniteElementSpaceHierarchy;

// The smoother which the multigrid preconditioner of a finite element space hierarchy
// uses for the given options. This is where it is decided whether the levels use blocks,
// for the solver (linalg/ksp.cpp) and for the assembly of its operators
// (models/spaceoperator.cpp), before any smoother is constructed: the block Jacobi
// relaxation, and with it the smoother around the coarse solve, only if a mesh of the
// hierarchy is a 3D mesh with flat elements (see linalg/vertexaggregates.hpp), and the
// point Jacobi relaxation otherwise. A hierarchy of a single level whose coarse solve is
// not smoothed has no smoother, and its mesh is not examined. Collective.
config::SmootherData ResolveSmoother(const config::SmootherData &smoother,
                                     const FiniteElementSpaceHierarchy &fespaces);

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

  // Optional smoother applied before and after the coarse-level solver.
  mutable std::unique_ptr<Solver<OperType>> B0;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<VecType> X, Y, R;
  mutable VecType Z0;

  // Enable timer contribution for Timer::KSP_COARSE_SOLVE.
  bool use_timer;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  // The smoothers of the levels are created from the given options: Chebyshev smoothers
  // of the relaxation method, or with the discrete gradients G (one per level), the
  // distributive relaxation smoothers which also smooth in the auxiliary spaces. With the
  // coarse option, a smoother is also applied before and after the coarse solver (and G
  // must then include the discrete gradient of the coarsest level). With print, a block
  // Jacobi relaxation says what became of its blocks.
  GeometricMultigridSolver(MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G, int cycle_it,
                           const config::SmootherData &smoother, bool print = false);
  // Chebyshev smoothers of the point Jacobi method on the levels.
  GeometricMultigridSolver(MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G, int cycle_it,
                           int smooth_it, int cheby_order, double cheby_sf_max,
                           double cheby_sf_min, bool cheby_4th_kind)
    : GeometricMultigridSolver(comm, std::move(coarse_solver), P, G, cycle_it,
                               config::SmootherData{.cheby_order = cheby_order,
                                                    .cheby_4th_kind = cheby_4th_kind,
                                                    .cheby_sf_max = cheby_sf_max,
                                                    .cheby_sf_min = cheby_sf_min,
                                                    .smooth_it = smooth_it})
  {
  }
  GeometricMultigridSolver(const IoData &iodata, MPI_Comm comm,
                           std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G = nullptr)
    : GeometricMultigridSolver(comm, std::move(coarse_solver), P, G,
                               iodata.solver.linear.mg_cycle_it,
                               iodata.solver.linear.mg_smoother)
  {
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void EnableTimer() { use_timer = true; }
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
