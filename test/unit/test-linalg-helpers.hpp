// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_LINALG_HELPERS_HPP
#define PALACE_TEST_LINALG_HELPERS_HPP

// Shared operators and vectors of the linear solver and initial guess tests.

#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace::test
{

// The symmetric stencil [s, d_i, s] on the local rows of each process, with
// d_i = d (4 + i % m).
class TridiagonalTestOperator : public Operator
{
  const double s, d;
  const int m;

public:
  TridiagonalTestOperator(int n, double s = -1.0, double d = 1.0, int m = 7)
    : Operator(n), s(s), d(d), m(m)
  {
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    auto *D = diag.HostWrite();
    for (int i = 0; i < height; i++)
    {
      D[i] = d * (4.0 + i % m);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y.SetSize(height);
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = d * (4.0 + i % m) * X[i] +
             s * (((i > 0) ? X[i - 1] : 0.0) + ((i < height - 1) ? X[i + 1] : 0.0));
    }
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// A random vector with entries (real and imaginary parts) in [0, 1), different on each
// process.
inline void SetRandomTestVector(Vector &x, int seed)
{
  linalg::SetRandom(Mpi::World(), x, 1 + seed + 1000 * Mpi::Rank(Mpi::World()));
}

inline void SetRandomTestVector(ComplexVector &x, int seed)
{
  SetRandomTestVector(x.Real(), seed);
  SetRandomTestVector(x.Imag(), seed + 500);
}

// Flag the vectors to use the device, as those of the solvers.
template <typename... VecTypes>
void UseDevice(VecTypes &...v)
{
  (v.UseDevice(true), ...);
}

}  // namespace palace::test

#endif  // PALACE_TEST_LINALG_HELPERS_HPP
