// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_KSP_PRECONDITIONER_HPP
#define PALACE_KSP_PRECONDITIONER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/petsc.hpp"

namespace palace
{

class IoData;

// Global method for preconditioner configuration and construction.
std::unique_ptr<mfem::Solver>
ConfigurePreconditioner(const IoData &iodata, const mfem::Array<int> &dbc_marker,
                        mfem::ParFiniteElementSpaceHierarchy &fespaces,
                        mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr);

//
// Class for preconditioning with interfaces to PETSc linear solvers or those from
// MFEM/Hypre.
//
class KspPreconditioner
{
private:
  // The actual preconditioner solver.
  std::unique_ptr<mfem::Solver> pc_;

  // Temporary vectors for preconditioner application.
  mutable mfem::Vector x_, y_;

  // Helper function for setup.
  void Init(int n);

public:
  KspPreconditioner(std::unique_ptr<mfem::Solver> &&pc) : pc_(std::move(pc)) {}
  KspPreconditioner(const IoData &iodata, const mfem::Array<int> &dbc_marker,
                    mfem::ParFiniteElementSpaceHierarchy &fespaces,
                    mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr)
    : pc_(ConfigurePreconditioner(iodata, dbc_marker, fespaces, aux_fespaces))
  {
    if (pc_->Height())
    {
      Init(pc_->Height());
    }
  }

  // Sets the matrix from which to contruct a preconditioner.
  void SetOperator(const mfem::Operator &op);
  void SetOperator(const std::vector<std::unique_ptr<mfem::Operator>> &ops,
                   const std::vector<std::unique_ptr<mfem::Operator>> *aux_ops = nullptr);

  // Application of the preconditioner.
  void Mult(const petsc::PetscParVector &x, petsc::PetscParVector &y) const;

  // Wrapper functions for PETSc PCSHELL.
  static PetscErrorCode PCSetUp(PC pc);
  static PetscErrorCode PCApply(PC pc, Vec x, Vec y);
  static PetscErrorCode PCDestroy(PC pc);
};

}  // namespace palace

#endif  // PALACE_KSP_PRECONDITIONER_HPP
