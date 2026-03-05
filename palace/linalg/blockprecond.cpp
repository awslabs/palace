// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "blockprecond.hpp"

#include <mfem.hpp>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

void CopySubVector(const Vector &src, Vector &dst, int src_offset, int dst_offset, int size)
{
  const bool use_dev = dst.UseDevice();
  const auto *sd = src.Read(use_dev);
  auto *dd = dst.Write(use_dev);
  mfem::forall_switch(use_dev, size, [=] MFEM_HOST_DEVICE(int i)
                      { dd[dst_offset + i] = sd[src_offset + i]; });
}

}  // namespace

template <>
void BlockDiagonalPreconditioner<Operator>::ExtractSubVector(const Vector &src, Vector &dst,
                                                             int offset, int size)
{
  CopySubVector(src, dst, offset, 0, size);
}

template <>
void BlockDiagonalPreconditioner<Operator>::InsertSubVector(const Vector &src, Vector &dst,
                                                            int offset, int size)
{
  CopySubVector(src, dst, 0, offset, size);
}

template <>
void BlockDiagonalPreconditioner<ComplexOperator>::ExtractSubVector(
    const ComplexVector &src, ComplexVector &dst, int offset, int size)
{
  CopySubVector(src.Real(), dst.Real(), offset, 0, size);
  CopySubVector(src.Imag(), dst.Imag(), offset, 0, size);
}

template <>
void BlockDiagonalPreconditioner<ComplexOperator>::InsertSubVector(const ComplexVector &src,
                                                                   ComplexVector &dst,
                                                                   int offset, int size)
{
  CopySubVector(src.Real(), dst.Real(), 0, offset, size);
  CopySubVector(src.Imag(), dst.Imag(), 0, offset, size);
}

template <typename OperType>
BlockDiagonalPreconditioner<OperType>::BlockDiagonalPreconditioner(
    int block0_size, std::unique_ptr<Solver<OperType>> &&pc0,
    std::unique_ptr<Solver<OperType>> &&pc1)
  : Solver<OperType>(), block0_size(block0_size), pc0(std::move(pc0)), pc1(std::move(pc1))
{
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::SetBlockOperators(const OperType &op0,
                                                              const OperType &op1)
{
  block0_size = op0.Height();
  pc0->SetOperator(op0);
  pc1->SetOperator(op1);
  this->height = op0.Height() + op1.Height();
  this->width = this->height;
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::SetOperator(const OperType &op)
{
  // No-op: use SetBlockOperators to set each block's operator independently.
  this->height = op.Height();
  this->width = op.Width();
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::Mult(const VecType &x, VecType &y) const
{
  const int n0 = block0_size;
  const int n1 = x.Size() - n0;

  // Lazy allocation of temporary vectors.
  if (x0.Size() != n0)
  {
    x0.SetSize(n0);
    y0.SetSize(n0);
    x0.UseDevice(true);
    y0.UseDevice(true);
  }
  if (x1.Size() != n1)
  {
    x1.SetSize(n1);
    y1.SetSize(n1);
    x1.UseDevice(true);
    y1.UseDevice(true);
    if (L10)
    {
      t1.SetSize(n1);
      t1.UseDevice(true);
    }
  }

  // Block 0: z0 = P0^{-1} r0
  ExtractSubVector(x, x0, 0, n0);
  pc0->Mult(x0, y0);
  InsertSubVector(y0, y, 0, n0);

  // Block 1: z1 = P1^{-1} (r1 - L10 z0)  [or z1 = P1^{-1} r1 if L10 is null]
  ExtractSubVector(x, x1, n0, n1);
  if (L10)
  {
    L10->Mult(y0, t1);
    linalg::AXPBY(1.0, x1, -1.0, t1);
    pc1->Mult(t1, y1);
  }
  else
  {
    pc1->Mult(x1, y1);
  }
  InsertSubVector(y1, y, n0, n1);
}

template class BlockDiagonalPreconditioner<Operator>;
template class BlockDiagonalPreconditioner<ComplexOperator>;

}  // namespace palace
