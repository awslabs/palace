// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SUPERLU_HPP
#define PALACE_SUPERLU_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_SUPERLU)

#include <memory>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for the SuperLU_DIST direct solver package.
//
class SuperLUSolver : public mfem::Solver
{
private:
  std::unique_ptr<mfem::SuperLURowLocMatrix> Aint;
  mfem::SuperLUSolver solver;

public:
  SuperLUSolver(MPI_Comm comm, int sym_fact_type, bool use_3d, int print_lvl);
  SuperLUSolver(MPI_Comm comm, const IoData &iodata, int print_lvl)
    : SuperLUSolver(comm,
                    (iodata.solver.linear.sym_fact_type ==
                     config::LinearSolverData::SymFactType::PARMETIS)
                        ? 2
                        : ((iodata.solver.linear.sym_fact_type ==
                            config::LinearSolverData::SymFactType::METIS)
                               ? 1
                               : 0),
                    iodata.solver.linear.superlu_3d, print_lvl)
  {
  }

  // Sets matrix associated with the SuperLU solver.
  void SetOperator(const mfem::Operator &op) override;

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { solver.Mult(x, y); }
};

}  // namespace palace

#endif

#endif  // PALACE_SUPERLU_HPP
