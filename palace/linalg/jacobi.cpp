// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "jacobi.hpp"

#include <mfem/general/forall.hpp>
#include "linalg/rap.hpp"

namespace palace
{

namespace
{

void GetInverseDiagonal(const ParOperator &A, Vector &dinv)
{
  dinv.SetSize(A.Height());
  A.AssembleDiagonal(dinv);
  dinv.Reciprocal();
}

void GetInverseDiagonal(const ComplexParOperator &A, ComplexVector &dinv)
{
  MFEM_VERIFY(A.HasReal() || A.HasImag(),
              "Invalid zero ComplexOperator for JacobiSmoother!");
  dinv.SetSize(A.Height());
  dinv.SetSize(A.Height());
  dinv = 0.0;
  if (A.HasReal())
  {
    A.Real()->AssembleDiagonal(dinv.Real());
  }
  if (A.HasImag())
  {
    A.Imag()->AssembleDiagonal(dinv.Imag());
  }
  dinv.Reciprocal();
}

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
  auto *YR = y.Real().ReadWrite();
  auto *YI = y.Imag().ReadWrite();
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
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  this->height = op.Height();
  this->width = op.Width();

  const auto *PtAP = dynamic_cast<const ParOperType *>(&op);
  MFEM_VERIFY(PtAP,
              "JacobiSmoother requires a ParOperator or ComplexParOperator operator!");
  GetInverseDiagonal(*PtAP, dinv);
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
