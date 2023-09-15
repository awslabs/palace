// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_OPERATOR_HPP
#define PALACE_LIBCEED_OPERATOR_HPP

#include <memory>
#include <vector>
#include <ceed.h>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace::ceed
{

// Wrapper class for libCEED's CeedOperator.
class Operator : public palace::Operator
{
protected:
  std::vector<CeedOperator> ops, ops_t;
  std::vector<CeedVector> u, v;
  Vector dof_multiplicity;
  mutable Vector temp_x;

public:
  Operator() = default;
  ~Operator() override;

  CeedOperator operator[](std::size_t i) const { return ops[i]; }
  auto Size() const { return ops.size(); }

  void AddOper(CeedOperator op, CeedOperator op_t = nullptr);

  void SetDofMultiplicity(Vector &&mult) { dof_multiplicity = mult; }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    AddMult(x, y);
  }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void MultTranspose(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    AddMultTranspose(x, y);
  }

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// Assemble a ceed::Operator as an mfem::SparseMatrix.
std::unique_ptr<mfem::SparseMatrix>
CeedOperatorFullAssemble(const Operator &op, bool skip_zeros = false, bool set = false);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_OPERATOR_HPP
