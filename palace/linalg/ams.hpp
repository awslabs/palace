// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AMS_HPP
#define PALACE_LINALG_AMS_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/operator.hpp"
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

  // Parameters used for preconditioner construction.
  const int cycle_type, space_dim, ams_it, ams_smooth_it, amg_agg_levels;
  const bool ams_singular;

  // Control print level for debugging.
  const int print;

  // Discrete gradient matrix.
  std::unique_ptr<ParOperator> G;

  // Nedelec interpolation matrix and its components, or, for p = 1, the mesh vertex
  // coordinates.
  std::unique_ptr<mfem::ParFiniteElementSpace> h1d_fespace;
  std::unique_ptr<ParOperator> Pi;
  std::unique_ptr<mfem::HypreParMatrix> Pix, Piy, Piz;
  std::unique_ptr<mfem::HypreParVector> x, y, z;

  // Helper function to set up the auxiliary objects required by the AMS solver.
  void ConstructAuxiliaryMatrices(mfem::ParFiniteElementSpace &nd_fespace,
                                  mfem::ParFiniteElementSpace &h1_fespace);

  // Helper function to construct and configure the AMS solver.
  void InitializeSolver();

public:
  // Constructor requires the ND space, but will construct the H1 and (H1)ᵈ spaces
  // internally as needed.
  HypreAmsSolver(mfem::ParFiniteElementSpace &nd_fespace,
                 mfem::ParFiniteElementSpace &h1_fespace, int cycle_it, int smooth_it,
                 int agg_coarsen, bool vector_interp, bool op_singular, int print);
  HypreAmsSolver(const IoData &iodata, mfem::ParFiniteElementSpace &nd_fespace,
                 mfem::ParFiniteElementSpace &h1_fespace, int print)
    : HypreAmsSolver(nd_fespace, h1_fespace,
                     iodata.solver.linear.pc_mg ? 1 : iodata.solver.linear.mg_cycle_it,
                     iodata.solver.linear.mg_smooth_it,
                     (iodata.problem.type == config::ProblemData::Type::TRANSIENT ||
                      iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC)
                         ? 1
                         : 0,
                     iodata.solver.linear.ams_vector,
                     (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC),
                     print)
  {
  }
  ~HypreAmsSolver() override;

  void SetOperator(const Operator &op) override;

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
