// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AMS_HPP
#define PALACE_LINALG_AMS_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

class FiniteElementSpace;

//
// A wrapper for Hypre's AMS solver.
//
class HypreAmsSolver : public mfem::HypreSolver
{
private:
  // The Hypre solver object.
  HYPRE_Solver ams;

  // Parameters used for preconditioner construction.
  const int cycle_type, space_dim, ams_it, ams_smooth_it;
  const bool ams_singular, agg_coarsen;

  // Target number of rows per process for the agglomeration of the coarse levels of the
  // auxiliary space AMG hierarchies during the solve (0: no agglomeration). Ignored
  // without the hypre patch for it in extern/patch/hypre.
  [[maybe_unused]] const int agglomeration_target;

  // Direct (distance-one) interpolation in the auxiliary space AMG hierarchies instead of
  // extended+i (distance two), see LinearSolverData::ams_aux_direct_interp.
  const bool direct_interp;

  // Control print level for debugging.
  const int print;

  // Discrete gradient matrix (not owned).
  const mfem::HypreParMatrix *G;

  // Nedelec interpolation matrix and its components, or, for p = 1, the mesh vertex
  // coordinates.
  std::unique_ptr<mfem::HypreParMatrix> Pi, Pix, Piy, Piz;
  std::unique_ptr<mfem::HypreParVector> x, y, z;

  // Helper function to set up the auxiliary objects required by the AMS solver.
  void ConstructAuxiliaryMatrices(FiniteElementSpace &nd_fespace,
                                  FiniteElementSpace &h1_fespace);

  // Helper function to construct and configure the AMS solver.
  void InitializeSolver();

public:
  // Constructor requires the ND space, but will construct the H1 and (H1)ᵈ spaces
  // internally as needed.
  HypreAmsSolver(FiniteElementSpace &nd_fespace, FiniteElementSpace &h1_fespace,
                 int cycle_it, int smooth_it, bool vector_interp, bool singular_op,
                 bool agg_coarsen, int agglomeration_target, int print,
                 bool direct_interp = false);
  HypreAmsSolver(const IoData &iodata, bool coarse_solver, FiniteElementSpace &nd_fespace,
                 FiniteElementSpace &h1_fespace, int print)
    : HypreAmsSolver(
          nd_fespace, h1_fespace,
          coarse_solver ? iodata.solver.linear.ams_max_it
                        : iodata.solver.linear.mg_cycle_it,
          iodata.solver.linear.mg_smooth_it, iodata.solver.linear.ams_vector_interp,
          iodata.solver.linear.ams_singular_op, iodata.solver.linear.amg_agg_coarsen,
          iodata.solver.linear.ams_agglomeration_target, print,
          iodata.solver.linear.ams_aux_direct_interp)
  {
  }
  ~HypreAmsSolver() override;

  void SetOperator(const Operator &op) override;

  // Flags a zero initial guess so that hypre can skip work on it.
  using mfem::HypreSolver::Mult;
  void Mult(const mfem::HypreParVector &b, mfem::HypreParVector &x) const override;

  operator HYPRE_Solver() const override { return ams; }

  HYPRE_PtrToParSolverFcn SetupFcn() const override
  {
    return (HYPRE_PtrToParSolverFcn)HYPRE_AMSSetup;
  }

  HYPRE_PtrToParSolverFcn SolveFcn() const override
  {
    return (HYPRE_PtrToParSolverFcn)HYPRE_AMSSolve;
  }
};

// Number of smoothing sweeps on the Nedelec space in each AMS iteration, applied before and
// after the auxiliary space corrections, when AMS is the coarse solver of geometric
// multigrid for a singular operator. Takes the global numbers of true dofs of the
// multigrid levels, coarsest first, and the number of smoothing sweeps of the multigrid
// levels (MGSmoothIts), which is a lower bound. For a singular curl-curl operator, the
// number of outer iterations depends strongly on this smoothing, which is cheap next to
// the auxiliary space corrections. More sweeps pay off as long as the coarse solve costs
// little next to the work on the finer levels, so their number grows with the ratio of the
// sizes of the finer levels to the size of the coarsest one, between 8 and 32.
int GetAmsSingularCoarseSmoothIts(const std::vector<HYPRE_BigInt> &level_sizes,
                                  int mg_smooth_it);

}  // namespace palace

#endif  // PALACE_LINALG_AMS_HPP
