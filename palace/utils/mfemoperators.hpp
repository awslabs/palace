// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MFEM_OPERATOR_HPP
#define PALACE_MFEM_OPERATOR_HPP

#include <functional>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

//
// Derived operator classes extending those already in MFEM.
//

// Wraps a transfer operator between finite element spaces to account for eliminated
// essential BC.
class ZeroWrapTransferOperator : public mfem::Operator
{
private:
  const mfem::TrueTransferOperator P;
  mfem::Array<int> coarse_dbc_tdof_list;

public:
  ZeroWrapTransferOperator(mfem::ParFiniteElementSpace &coarse_fespace,
                           mfem::ParFiniteElementSpace &fine_fespace,
                           const mfem::Array<int> &dbc_marker)
    : P(coarse_fespace, fine_fespace)
  {
    height = P.Height();
    width = P.Width();
    coarse_fespace.GetEssentialTrueDofs(dbc_marker, coarse_dbc_tdof_list);
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { P.Mult(x, y); }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    P.MultTranspose(x, y);
    y.SetSubVector(coarse_dbc_tdof_list, 0.0);
  }
};

// Wraps a reference to an existing, externally owned operator.
class ReferenceOperator : public mfem::Operator
{
private:
  const mfem::Operator &op;

public:
  ReferenceOperator(const mfem::Operator &oper)
    : mfem::Operator(oper.Height(), oper.Width()), op(oper)
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { op.Mult(x, y); }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    op.MultTranspose(x, y);
  }
};

// Wrap a sequence of operators of the same dimensions and optional coefficients.
class SumOperator : public mfem::Operator
{
private:
  std::vector<std::pair<std::reference_wrapper<const mfem::Operator>, double>> op;

public:
  SumOperator(int s) : mfem::Operator(s) {}
  SumOperator(int h, int w) : mfem::Operator(h, w) {}

  void AddOperator(const mfem::Operator &oper, double c = 1.0)
  {
    MFEM_VERIFY(oper.Height() == height && oper.Width() == width,
                "Invalid Operator dimensions for SumOperator!");
    op.emplace_back(std::cref(oper), c);
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    y = 0.0;
    for (const auto &[oper, c] : op)
    {
      oper.get().AddMult(x, y, c);
    }
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    y = 0.0;
    for (const auto &[oper, c] : op)
    {
      oper.get().AddMultTranspose(x, y, c);
    }
  }
};

}  // namespace palace

#endif  // PALACE_MFEM_OPERATOR_HPP
