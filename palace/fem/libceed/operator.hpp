// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_OPERATOR_HPP
#define PALACE_LIBCEED_OPERATOR_HPP

#include <memory>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

namespace hypre
{

class HypreCSRMatrix;

}  // namespace hypre

namespace ceed
{

//
// Wrapper class for libCEED's CeedOperator, supporting composite operator construction and
// application with multiple threads.
//
class Operator : public palace::Operator
{
protected:
  std::vector<CeedOperator> op, op_t;
  std::vector<CeedVector> u, v;
  Vector dof_multiplicity;
  mutable Vector temp;

public:
  Operator(int h, int w);
  ~Operator() override;

  CeedOperator operator[](std::size_t i) const { return op[i]; }

  auto Size() const { return op.size(); }

  void AddSubOperator(CeedOperator sub_op, CeedOperator sub_op_t = nullptr);

  void Finalize();

  void DestroyAssemblyData() const;

  void SetDofMultiplicity(Vector &&mult) { dof_multiplicity = std::move(mult); }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// A symmetric ceed::Operator replaces *MultTranspose with *Mult (by default, libCEED
// operators do not have a transpose operation).
class SymmetricOperator : public Operator
{
public:
  using Operator::Operator;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
  void AddMultTranspose(const Vector &x, Vector &y, double a = 1.0) const override
  {
    AddMult(x, y, a);
  }
};

// Assemble a ceed::Operator as a CSR matrix.
std::unique_ptr<hypre::HypreCSRMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                                bool skip_zeros, bool set);

// Construct a coarse-level ceed::Operator, reusing the quadrature data and quadrature
// function from the fine-level operator. Only available for square, symmetric operators
// (same input and output spaces).
std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_OPERATOR_HPP
