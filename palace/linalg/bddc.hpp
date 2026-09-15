// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BDDC_HPP
#define PALACE_LINALG_BDDC_HPP

#include <mfem.hpp>

#if defined(PALACE_WITH_SLEPC)

#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

//
// Two-level balancing domain decomposition by constraints (BDDC) preconditioner using
// PETSc's PCBDDC. Unlike one-level RAS, BDDC adds a coarse correction built from primal
// constraints on the subdomain interface, which keeps the iteration count nearly
// independent of the number of subdomains.
//
// Palace's parallel operators are A = Pᵀ A_loc P, while PETSc's MATIS format needs
// A = Σ Rᵀ Ã R with R a 0/1 injection. On a conforming mesh each local dof row of P has a
// single ±1 entry, so P = S R with S = diag(±1) and Ã = S A_loc S. Nonconforming (hanging)
// dof rows couple several true dofs and are not supported.
//
class BddcSolver : public Solver<Operator>
{
private:
  // Finite element space defining the dof layout of the operator (not owned).
  const FiniteElementSpace &fespace;

  // Auxiliary H1 space used to build the discrete gradient for H(curl) problems, which
  // BDDC needs to handle the curl-curl kernel. May be nullptr for H1 problems.
  const FiniteElementSpace *aux_fespace;

  // Print level.
  int print;

  // Opaque PETSc state, so that PETSc headers stay out of this header.
  struct PetscData;
  std::unique_ptr<PetscData> petsc;

public:
  BddcSolver(const FiniteElementSpace &fespace, const FiniteElementSpace *aux_fespace,
             int print);
  ~BddcSolver() override;

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_BDDC_HPP
