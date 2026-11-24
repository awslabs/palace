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
                 bool agg_coarsen, int print);
  HypreAmsSolver(const IoData &iodata, bool coarse_solver, FiniteElementSpace &nd_fespace,
                 FiniteElementSpace &h1_fespace, int print)
    : HypreAmsSolver(
          nd_fespace, h1_fespace,
          coarse_solver ? iodata.solver.linear.ams_max_it
                        : iodata.solver.linear.mg_cycle_it,
          iodata.solver.linear.mg_smooth_it, iodata.solver.linear.ams_vector_interp,
          iodata.solver.linear.ams_singular_op, iodata.solver.linear.amg_agg_coarsen, print)
  {
  }
  ~HypreAmsSolver() override;

  void SetOperator(const Operator &op) override;

  // STUDY ONLY (ams-robust): AMSR_AMS_ALPHA=1 passes this vector Laplacian (H1(1)^3,
  // byVDIM) to HYPRE_AMSSetAlphaPoissonMatrix (used by hypre for cycle types < 10).
  inline static mfem::HypreParMatrix *study_alpha = nullptr;
  // STUDY ONLY (ams-robust): AMSR_AMS_PIOP = rap | rapdiv | hx: a user Pi (hypre's
  // vertex-coordinate Pi, built here, byVDIM columns) and a user Pi-space operator:
  // Pi^T A Pi, Pi^T A Pi + study_divdiv (alpha-weighted nodal div-div), or study_alpha.
  inline static mfem::HypreParMatrix *study_divdiv = nullptr;

private:
  std::unique_ptr<mfem::HypreParMatrix> study_Pi, study_APi;
  std::vector<std::unique_ptr<mfem::SparseMatrix>> study_Ed;
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> study_keep;
  mfem::Array<HYPRE_BigInt> study_rows, study_cols;

public:

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
