// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bddc.hpp"

//#if defined(PALACE_WITH_PETSC)

namespace palace
{

BDDCSolver::BDDCSolver(MPI_Comm comm, int print)
  : mfem::Solver(), comm(comm)
{

  // stuff
  // MFEM's two constructors for PetscBDDCSolver
/*
PetscBDDCSolver::PetscBDDCSolver(PetscParMatrix &A,
                                 const PetscBDDCSolverParams &opts,
                                 const std::string &prefix)
PetscBDDCSolver::PetscBDDCSolver(MPI_Comm comm, Operator &op,
                                 const PetscBDDCSolverParams &opts,
                                 const std::string &prefix)
*/
}

void BDDCSolver::SetOperator(const Operator &op)
{
  //const auto *pA = dynamic_cast<const mfem::PetscParMatrix *>(&op);
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op); // op -> HypreParMatrix
  std::cout << "bddc.cpp L30\n";
  //auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATIS); // op -> PetscParMatrix
  auto pA = std::make_unique<mfem::PetscParMatrix>(hA, Operator::PETSC_MATIS); // HypreParMatrix -> PetscParMatrix
  //auto *pA = const_cast<mfem::PetscParMatrix*>(dynamic_cast<const mfem::PetscParMatrix*>(&op));
  std::cout << "bddc.cpp L32\n";
  mfem::PetscBDDCSolverParams opts;
  std::cout << "bddc.cpp L34\n";
  solver = std::make_unique<mfem::PetscBDDCSolver>(*pA, opts, "prec_bddc_");
  std::cout << "bddc.cpp L36\n";
}

}  // namespace palace

//#endif