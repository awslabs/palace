// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
#define PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
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
  // Number of V-cycles per preconditioner application.
  const int pc_it;

  // System matrices at each multigrid level and prolongation operators (not owned).
  std::vector<const Operator *> A_, P_;

  // Essential Dirichlet boundary conditions at each level (not owned).
  std::vector<const mfem::Array<int> *> dbc_tdof_lists_;

  // Smoothers for each level. Coarse level solver is B_[0].
  std::vector<std::unique_ptr<mfem::Solver>> B_;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<Vector> x_, y_, r_;
  mutable std::vector<std::vector<Vector>> xrefs_, yrefs_, rrefs_;
  mutable std::vector<mfem::Array<Vector *>> X_, Y_, R_;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(std::unique_ptr<mfem::Solver> &&coarse_solver,
                           mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           mfem::ParFiniteElementSpaceHierarchy *aux_fespaces, int cycle_it,
                           int smooth_it, int cheby_order);
  GeometricMultigridSolver(const IoData &iodata,
                           std::unique_ptr<mfem::Solver> &&coarse_solver,
                           mfem::ParFiniteElementSpaceHierarchy &fespaces,
                           mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
    : GeometricMultigridSolver(std::move(coarse_solver), fespaces, aux_fespaces,
                               iodata.solver.linear.mg_cycle_it,
                               iodata.solver.linear.mg_smooth_it,
                               iodata.solver.linear.mg_smooth_order)
  {
  }

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "GeometricMultigridSolver, use the other signature instead!");
  }
  void SetOperator(const std::vector<std::unique_ptr<Operator>> &ops,
                   const std::vector<std::unique_ptr<Operator>> *aux_ops = nullptr);

  void Mult(const Vector &x, Vector &y) const override
  {
    mfem::Array<const Vector *> X(1);
    mfem::Array<Vector *> Y(1);
    X[0] = &x;
    Y[0] = &y;
    ArrayMult(X, Y);
  }

  void ArrayMult(const mfem::Array<const Vector *> &X,
                 mfem::Array<Vector *> &Y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
