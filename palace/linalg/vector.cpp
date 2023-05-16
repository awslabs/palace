// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vector.hpp"

#include <general/forall.hpp>
#include "linalg/complex.hpp"
#include "linalg/operator.hpp"
#include "utils/communication.hpp"

namespace palace::linalg
{

HYPRE_BigInt GlobalSize(MPI_Comm comm, const Vector &x)
{
  HYPRE_BigInt N = x.Size();
  Mpi::GlobalSum(1, &N, comm);
  return N;
}

void SetRandom(MPI_Comm comm, Vector &x, int seed)
{
  seed *= Mpi::Rank(comm) + 1;
  x.Randomize(seed);
}

void SetRandomSign(MPI_Comm comm, Vector &x, int seed)
{
  SetRandom(comm, x, seed);
  const int N = x.Size();
  auto *X = x.ReadWrite();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i)
               { X[i] = (X[i] > 0.0) ? 1.0 : ((X[i] < 0.0) ? -1.0 : 0.0); });
}

double Norml2(MPI_Comm comm, const Vector &x)
{
  return std::sqrt(mfem::InnerProduct(comm, x, x));
}

double Normlinf(MPI_Comm comm, const Vector &x)
{
  double norm = x.Normlinf();
  Mpi::GlobalMax(1, &norm, comm);
  return norm;
}

double Norml1(MPI_Comm comm, const Vector &x)
{
  double norm = x.Norml1();
  Mpi::GlobalSum(1, &norm, comm);
  return norm;
}

double Normalize(MPI_Comm comm, Vector &x)
{
  double norm = Norml2(comm, x);
  MFEM_ASSERT(norm > 0.0, "Zero vector norm in normalization!");
  x *= 1.0 / norm;
  return norm;
}

double Normalize(MPI_Comm comm, Vector &x, const Operator &B, Vector &Bx)
{
  B.Mult(x, Bx);
  double norm = std::sqrt(mfem::InnerProduct(comm, x, Bx));
  MFEM_ASSERT(norm > 0.0, "Zero vector norm in normalization!");
  x *= 1.0 / norm;
  return norm;
}

double Normalize(MPI_Comm comm, ComplexVector &x, const Operator &B, ComplexVector &Bx)
{
  // For SPD B, xᴴ B x is real.
  B.Mult(x.Real(), Bx.Real());
  B.Mult(x.Imag(), Bx.Imag());
  double norm = std::sqrt(mfem::InnerProduct(comm, x, Bx));
  MFEM_ASSERT(norm > 0.0, "Zero vector norm in normalization!");
  x *= 1.0 / norm;
  return norm;
}

}  // namespace palace::linalg
