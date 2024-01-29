// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jacobi.hpp"

#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

template <bool Transpose = false>
inline void Apply(const Vector &dinv, const Vector &x, Vector &y)
{
  const int N = dinv.Size();
  const auto *DI = dinv.Read();
  const auto *X = x.Read();
  auto *Y = y.Write();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = DI[i] * X[i]; });
}

template <bool Transpose = false>
inline void Apply(const ComplexVector &dinv, const ComplexVector &x, ComplexVector &y)
{
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read();
  const auto *DII = dinv.Imag().Read();
  const auto *XR = x.Real().Read();
  const auto *XI = x.Imag().Read();
  auto *YR = y.Real().Write();
  auto *YI = y.Imag().Write();
  if constexpr (!Transpose)
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] = DIR[i] * XR[i] - DII[i] * XI[i];
                   YI[i] = DII[i] * XR[i] + DIR[i] * XI[i];
                 });
  }
  else
  {
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   YR[i] = DIR[i] * XR[i] + DII[i] * XI[i];
                   YI[i] = -DII[i] * XR[i] + DIR[i] * XI[i];
                 });
  }
}

}  // namespace

template <typename OperType>
void JacobiSmoother<OperType>::SetOperator(const OperType &op)
{
  dinv.SetSize(op.Height());
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  this->height = op.Height();
  this->width = op.Width();
}

template <typename OperType>
void JacobiSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  MFEM_ASSERT(!this->initial_guess, "JacobiSmoother does not use initial guess!");
  Apply(dinv, x, y);
}

template class JacobiSmoother<Operator>;
template class JacobiSmoother<ComplexOperator>;

}  // namespace palace
