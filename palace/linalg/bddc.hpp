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

class FiniteElementSpace;

//
// A wrapper for PETSc's BDDC solver.
//
class BDDCSolver : public mfem::Solver
{
private:
  MPI_Comm comm;
  std::unique_ptr<mfem::PetscBDDCSolver> solver;
  std::unique_ptr<mfem::PetscLinearSolver> solver2;
  mfem::Array<int> ess_tdof_list, block_ess_tdof_list, nat_tdof_list;
  FiniteElementSpace &fespace;
  bool complex_coarse;

public:
  BDDCSolver(const IoData &iodata, FiniteElementSpace &fespace, int print = 0);

  void SetOperator(const Operator &op) override;

  //void Mult(const Vector &x, Vector &y) const override { solver2->Mult(x, y);}
  void Mult(const Vector &x, Vector &y) const override {
    //std::cout << "BDDC Mult x.Size(): " << x.Size() << " solver->Width(): " << solver->Width() << "\n";
    solver->Mult(x, y);
  }

  //int Height() const {return solver->Height();}
  //int Width() const {return solver->Width();}
};

}  // namespace palace

//#endif

#endif  // PALACE_LINALG_BDDC_HPP
