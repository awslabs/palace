// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AMS_HPP
#define PALACE_LINALG_AMS_HPP

#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for Hypre's AMS solver.
//
class HypreAmsSolver : public mfem::HypreSolver
{
private:
  // The Hypre solver object.
  HYPRE_Solver ams;

  // Discrete gradient matrix.
  std::unique_ptr<mfem::HypreParMatrix> G;

  // Nedelec interpolation matrix and its components (used even for p = 1).
  std::unique_ptr<mfem::HypreParMatrix> Pi, Pix, Piy, Piz;

  // Parameters used for preconditioner construction.
  const int cycle_type, sdim, ams_it, ams_smooth_it, agg_levels;
  const bool ams_singular;

  // Control print level for debugging.
  const int print;

  // Helper functions to construct the AMS solver and required auxiliary space matrices.
  void Initialize();
  void ConstructAuxiliaryMatrices(mfem::ParFiniteElementSpace &nd_fespace,
                                  mfem::ParFiniteElementSpace *h1_fespace = nullptr);

public:
  // Constructor requires the ND space, but will construct the H1 and (H1)ᵈ spaces
  // internally as needed.
  HypreAmsSolver(mfem::ParFiniteElementSpace &nd_fespace,
                 mfem::ParFiniteElementSpace *h1_fespace, int cycle_it, int smooth_it,
                 int agg_coarsen, bool vector_interp, bool op_singular, int print_lvl);
  HypreAmsSolver(const IoData &iodata, mfem::ParFiniteElementSpace &nd_fespace,
                 mfem::ParFiniteElementSpace *h1_fespace, int print_lvl)
    : HypreAmsSolver(nd_fespace, h1_fespace,
                     iodata.solver.linear.mat_gmg ? 1 : iodata.solver.linear.mg_cycle_it,
                     iodata.solver.linear.mg_smooth_it,
                     (iodata.problem.type == config::ProblemData::Type::TRANSIENT ||
                      iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC)
                         ? 1
                         : 0,
                     iodata.solver.linear.ams_vector,
                     (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC),
                     print_lvl)
  {
  }
  ~HypreAmsSolver() override;

  // Sets matrix associated with the AMS solver.
  void SetOperator(const mfem::Operator &op) override;

  // The typecast to HYPRE_Solver returns the internal ams object.
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

}  // namespace palace

#endif  // PALACE_LINALG_AMS_HPP
