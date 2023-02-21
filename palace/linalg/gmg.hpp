// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_GEOMETRIC_MULTIGRID_HPP
#define PALACE_GEOMETRIC_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// Geometric multigrid preconditioner using a given coarse solver for the provided
// hierarchy of finite element spaces. Optionally can be configured to use auxiliary space
// smoothing at each level.
//
class GeometricMultigridSolver : public mfem::Solver
{
private:
  // Reference to the underlying finite element space hierarchy used to construct the
  // multilevel preconditioner.
  const mfem::ParFiniteElementSpaceHierarchy &fespaces_;

  // System matrices at each multigrid level (not owned).
  std::vector<mfem::Operator *> A_;

  // Smoothers for each level. Coarse level solver is B_[0].
  std::vector<std::unique_ptr<mfem::Solver>> B_;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<mfem::Vector> x_, y_, r_;
  mutable std::vector<mfem::Array<mfem::Vector *>> X_, Y_, R_;
  mutable std::vector<std::vector<mfem::Vector>> xrefs_, yrefs_, rrefs_;

  // Number of V-cycles per preconditioner application.
  const int pc_it;

  // Returns prolongation operator at given level.
  const mfem::Operator &GetProlongationAtLevel(int l) const
  {
    return *fespaces_.GetProlongationAtLevel(l);
  }

  // Returns the number of levels.
  int GetNumLevels() const { return fespaces_.GetNumLevels(); }

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(std::unique_ptr<mfem::Solver> &&coarse_solver,
                           const mfem::Array<int> &dbc_marker,
                           mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           mfem::ParFiniteElementSpaceHierarchy *aux_fespaces, int cycle_it,
                           int smooth_it, int cheby_order);
  GeometricMultigridSolver(const IoData &iodata,
                           std::unique_ptr<mfem::Solver> &&coarse_solver,
                           const mfem::Array<int> &dbc_marker,
                           mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
    : GeometricMultigridSolver(std::move(coarse_solver), dbc_marker, fespaces, aux_fespaces,
                               iodata.solver.linear.mg_cycle_it,
                               iodata.solver.linear.mg_smooth_it,
                               iodata.solver.linear.mg_smooth_order)
  {
  }

  // Sets the matrices from which to contruct a multilevel preconditioner.
  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "GeometricMultigridSolver, use the other signature instead!");
  }
  void SetOperator(const std::vector<std::unique_ptr<mfem::Operator>> &ops,
                   const std::vector<std::unique_ptr<mfem::Operator>> *aux_ops = nullptr);

  // Application of the solver.
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    mfem::Array<const mfem::Vector *> X(1);
    mfem::Array<mfem::Vector *> Y(1);
    X[0] = &x;
    Y[0] = &y;
    ArrayMult(X, Y);
  }

  void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                 mfem::Array<mfem::Vector *> &Y) const override
  {
    // Initialize.
    const int m = GetNumLevels(), nrhs = X.Size();
    MFEM_VERIFY(!iterative_mode, "Geometric multigrid solver does not use iterative_mode!");

    MFEM_VERIFY(m > 1 || pc_it == 1,
                "Single-level geometric multigrid will not work with multiple iterations!");
    if (nrhs * height != x_[m - 1].Size())
    {
      for (int l = 0; l < m; l++)
      {
        MFEM_VERIFY(A_[l], "Missing operator for geometric multigrid level " << l << "!");
        x_[l].SetSize(nrhs * A_[l]->Height());
        y_[l].SetSize(nrhs * A_[l]->Height());
        r_[l].SetSize(nrhs * A_[l]->Height());
      }
    }
    for (int l = 0; l < m; l++)
    {
      xrefs_[l].resize(nrhs);
      yrefs_[l].resize(nrhs);
      rrefs_[l].resize(nrhs);
      X_[l].SetSize(nrhs);
      Y_[l].SetSize(nrhs);
      R_[l].SetSize(nrhs);
      for (int j = 0; j < nrhs; j++)
      {
        xrefs_[l][j].MakeRef(x_[l], j * A_[l]->Height(), A_[l]->Height());
        yrefs_[l][j].MakeRef(y_[l], j * A_[l]->Height(), A_[l]->Height());
        rrefs_[l][j].MakeRef(r_[l], j * A_[l]->Height(), A_[l]->Height());
        X_[l][j] = &xrefs_[l][j];
        Y_[l][j] = &yrefs_[l][j];
        R_[l][j] = &rrefs_[l][j];
      }
    }

    // Apply V-cycle.
    for (int j = 0; j < nrhs; j++)
    {
      *X_[m - 1][j] = *X[j];
    }
    for (int it = 0; it < pc_it; it++)
    {
      VCycle(m - 1, (it > 0));
    }
    for (int j = 0; j < nrhs; j++)
    {
      *Y[j] = *Y_[m - 1][j];
    }
  }
};

}  // namespace palace

#endif  // PALACE_GEOMETRIC_MULTIGRID_HPP
