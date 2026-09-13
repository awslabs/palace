// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/densematrix.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace
{

// Rotation matrix for a rotation by the given angle (radians) about the given axis
// (normalized internally), via the Rodrigues formula.
mfem::DenseMatrix RotationMatrix(double ax, double ay, double az, double angle)
{
  const double norm = std::sqrt(ax * ax + ay * ay + az * az);
  const double u[3] = {ax / norm, ay / norm, az / norm};
  const double c = std::cos(angle), s = std::sin(angle);
  mfem::DenseMatrix R(3);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      R(i, j) = (1.0 - c) * u[i] * u[j] + ((i == j) ? c : 0.0);
    }
  }
  R(0, 1) -= s * u[2];
  R(0, 2) += s * u[1];
  R(1, 0) += s * u[2];
  R(1, 2) -= s * u[0];
  R(2, 0) -= s * u[1];
  R(2, 1) += s * u[0];
  return R;
}

// Construct R * diag(d0, d1, d2) * Rᵀ, exactly symmetric by construction.
mfem::DenseMatrix RotatedDiagonal(const mfem::DenseMatrix &R, double d0, double d1,
                                  double d2)
{
  const double d[3] = {d0, d1, d2};
  mfem::DenseMatrix M(3);
  for (int i = 0; i < 3; i++)
  {
    for (int j = i; j < 3; j++)
    {
      double v = 0.0;
      for (int k = 0; k < 3; k++)
      {
        v += R(i, k) * d[k] * R(j, k);
      }
      M(i, j) = M(j, i) = v;
    }
  }
  return M;
}

// Relative Frobenius-norm residual ‖A − B‖ / ‖B‖.
double RelDiff(const mfem::DenseMatrix &A, const mfem::DenseMatrix &B)
{
  mfem::DenseMatrix D(A);
  D -= B;
  return D.FNorm() / B.FNorm();
}

}  // namespace

TEST_CASE("MatrixSqrt matches the closed-form eigendecomposition", "[densematrix][Serial]")
{
  // We compare against exact ground truth instead: RotatedDiagonal(R, d0, d1,
  // d2) is R · diag(d0, d1, d2) · Rᵀ, so its eigenvalues are exactly {d0, d1,
  // d2} and its eigenvectors are the columns of R by construction. Hence the
  // principal square root is exactly R · diag(√d0, √d1, √d2) · Rᵀ, which we
  // build with the same helper. Matching it element-wise directly validates
  // that the eigenvalues (magnitude and sign) and eigenvectors used inside
  // MatrixFunction are correct.
  constexpr double tol = 1.0e-12;
  const auto R = RotationMatrix(1.0, 2.0, 3.0, 0.646);

  SECTION("2x2")
  {
    // 2x2 symmetric M = Q·diag(4, 9)·Qᵀ for a 2D rotation Q(θ); its principal square root
    // is Q·diag(2, 3)·Qᵀ. Exercises the N == 2 path (mfem CalcEigenvalues<2>).
    const double theta = 0.7, c = std::cos(theta), s = std::sin(theta);
    const auto rotated_2x2 = [&](double d0, double d1)
    {
      mfem::DenseMatrix A(2);
      A(0, 0) = c * c * d0 + s * s * d1;
      A(1, 1) = s * s * d0 + c * c * d1;
      A(0, 1) = A(1, 0) = c * s * (d0 - d1);
      return A;
    };
    REQUIRE_THAT(RelDiff(linalg::MatrixSqrt(rotated_2x2(4.0, 9.0)), rotated_2x2(2.0, 3.0)),
                 WithinAbs(0.0, tol));
  }

  SECTION("distinct eigenvalues")
  {
    const auto M = RotatedDiagonal(R, 8.9, 9.6, 11.3);
    const auto expected =
        RotatedDiagonal(R, std::sqrt(8.9), std::sqrt(9.6), std::sqrt(11.3));
    REQUIRE_THAT(RelDiff(linalg::MatrixSqrt(M), expected), WithinAbs(0.0, tol));
  }

  SECTION("repeated eigenvalue")
  {
    // Degenerate spectrum: the eigenvectors of the repeated pair are not unique, but the
    // principal square root R · diag(√d) · Rᵀ still is, so the exact comparison holds.
    const auto M = RotatedDiagonal(R, 9.27, 9.27, 11.34);
    const auto expected =
        RotatedDiagonal(R, std::sqrt(9.27), std::sqrt(9.27), std::sqrt(11.34));
    REQUIRE_THAT(RelDiff(linalg::MatrixSqrt(M), expected), WithinAbs(0.0, tol));
  }
}

TEST_CASE("MatrixPow of symmetric matrices", "[densematrix][Serial]")
{
  constexpr double tol = 1.0e-12;

  // M^(-1/2) · M · M^(-1/2) should recover the identity. Check both a fully coupled SPD
  // input with distinct eigenvalues and a degenerate (rotated uniaxial) one, since the
  // general branch is what the fix targets.
  const auto R = RotationMatrix(1.0, 2.0, 3.0, 0.646);
  mfem::DenseMatrix I(3);
  I = 0.0;
  for (int i = 0; i < 3; i++)
  {
    I(i, i) = 1.0;
  }

  const auto check_inverse_sqrt = [&](const mfem::DenseMatrix &M)
  {
    const mfem::DenseMatrix P = linalg::MatrixPow(M, -0.5);
    mfem::DenseMatrix PM(3), PMP(3);
    mfem::Mult(P, M, PM);
    mfem::Mult(PM, P, PMP);
    REQUIRE_THAT(RelDiff(PMP, I), WithinAbs(0.0, tol));
  };

  SECTION("distinct eigenvalues")
  {
    check_inverse_sqrt(RotatedDiagonal(R, 8.9, 9.6, 11.3));
  }

  SECTION("repeated eigenvalue")
  {
    check_inverse_sqrt(RotatedDiagonal(R, 9.27, 9.27, 11.34));
  }
}
