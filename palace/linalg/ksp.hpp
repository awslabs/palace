// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_SOLVER_HPP
#define PALACE_LINALG_KSP_SOLVER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class ComplexParOperator;
class ComplexVector;
class IoData;

class KspSolver : public mfem::Solver
{
protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<mfem::IterativeSolver> ksp_;
  std::unique_ptr<mfem::Solver> pc_;

  // Counters for number of calls to Mult method for linear solves, and cumulative number
  // of iterations.
  mutable int ksp_mult, ksp_mult_it;

protected:
  KspSolver() : ksp_(nullptr), pc_(nullptr), ksp_mult(0), ksp_mult_it(0) {}

  void SetOperatorFinalize(const Operator &op);

public:
  KspSolver(const IoData &iodata, mfem::ParFiniteElementSpaceHierarchy &fespaces,
            mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  KspSolver(std::unique_ptr<mfem::IterativeSolver> &&ksp,
            std::unique_ptr<mfem::Solver> &&pc);

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIter() const { return ksp_mult_it; }

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for KspSolver, you "
               "must specify the preconditioner operator as well!");
  }

  virtual void SetOperator(const Operator &op, const Operator &pc_op);
  virtual void
  SetOperator(const Operator &op, const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
              const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr);

  void Mult(const Vector &x, Vector &y) const override;
};

class ComplexKspSolver : public KspSolver
{
public:
  ComplexKspSolver(const IoData &iodata, mfem::ParFiniteElementSpaceHierarchy &fespaces,
                   mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  ComplexKspSolver(std::unique_ptr<mfem::IterativeSolver> &&ksp,
                   std::unique_ptr<mfem::Solver> &&pc);

  using KspSolver::SetOperator;
  void SetOperator(const Operator &op, const Operator &pc_op) override
  {
    MFEM_ABORT("SetOperator with a real-valued operator is not implemented for "
               "ComplexKspSolver, use the complex-valued signature instead!");
  }
  void SetOperator(
      const Operator &op, const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
      const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr) override
  {
    MFEM_ABORT("SetOperator with a real-valued operator is not implemented for "
               "ComplexKspSolver, use the complex-valued signature instead!");
  }

  void SetOperator(const ComplexOperator &op, const Operator &pc_op);
  void SetOperator(const ComplexOperator &op,
                   const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
                   const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr);

  void Mult(const Vector &x, Vector &y) const override
  {
    MFEM_ABORT("Mult with a real-valued vector is not implemented for "
               "ComplexKspSolver, use the complex-valued signature instead!");
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_KSP_SOLVER_HPP
