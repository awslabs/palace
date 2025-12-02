// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BDDC_HPP
#define PALACE_LINALG_BDDC_HPP

#include <mfem.hpp>

//#if defined(PALACE_WITH_PETSC)

#include "linalg/petsc.hpp"

//#if !defined(PETSC_USE_COMPLEX)
//#error "SLEPc interface requires PETSc compiled with complex scalars!"
//#endif

#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for PETSc's BDDC solver.
//
class BDDCSolver : public mfem::Solver //PetscBDDCSolver
{
private:
  MPI_Comm comm;
  std::unique_ptr<mfem::PetscBDDCSolver> solver;

public:
  BDDCSolver(MPI_Comm comm, int print = 0); // figure out what we need to construct...
  BDDCSolver(const IoData &iodata, MPI_Comm comm, int print)
    : BDDCSolver(comm, print)
  {
  }

  void SetOperator(const Operator &op) override;

  void Mult(const Vector &x, Vector &y) const override { std::cout << "bddc.hpp before mult \n"; solver->Mult(x, y);  std::cout << "bddc.hpp after mult \n";}
};

}  // namespace palace

//#endif

#endif  // PALACE_LINALG_BDDC_HPP
