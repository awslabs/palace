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
  const bool use_dev = dinv.UseDevice() || x.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *X = x.Read(use_dev);
  auto *Y = y.Write(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = DI[i] * X[i]; });
}

template <bool Transpose = false>
inline void Apply(const ComplexVector &dinv, const ComplexVector &x, ComplexVector &y)
{
  const bool use_dev = dinv.UseDevice() || x.UseDevice() || y.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *XR = x.Real().Read(use_dev);
  const auto *XI = x.Imag().Read(use_dev);
  auto *YR = y.Real().Write(use_dev);
  auto *YI = y.Imag().Write(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          YR[i] = DIR[i] * XR[i] - DII[i] * XI[i];
                          YI[i] = DII[i] * XR[i] + DIR[i] * XI[i];
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
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
  dinv.UseDevice(true);
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
