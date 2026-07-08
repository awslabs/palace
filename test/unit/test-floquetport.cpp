// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Unit tests for Floquet port computational functions.
// Tests call the actual static methods in FloquetPortData.

#include <cmath>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "models/floquetportoperator.hpp"

namespace palace
{

using namespace Catch::Matchers;

// =============================================================================
// Reciprocal lattice: FloquetPortData::ComputeReciprocalLattice
// =============================================================================
TEST_CASE("Floquet Reciprocal Lattice", "[floquetport][Serial]")
{
  SECTION("Square cell: a₁ = (1,0,0), a₂ = (0,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 1.0;
    a2 = 0.0;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
    CHECK_THAT(b1.Norml2(), WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(b2.Norml2(), WithinRel(2.0 * M_PI, 1e-10));
  }

  SECTION("Rectangular cell: a₁ = (4,0,0), a₂ = (0,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 4.0;
    a2 = 0.0;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
    // |b₁| = 2π/4 = π/2, |b₂| = 2π
    CHECK_THAT(b1.Norml2(), WithinRel(M_PI / 2.0, 1e-10));
    CHECK_THAT(b2.Norml2(), WithinRel(2.0 * M_PI, 1e-10));
    // b₁ along x, b₂ along y
    CHECK_THAT(b1(0), WithinRel(M_PI / 2.0, 1e-10));
    CHECK_THAT(std::abs(b1(1)), WithinAbs(0.0, 1e-10));
    CHECK_THAT(b2(1), WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(std::abs(b2(0)), WithinAbs(0.0, 1e-10));
  }

  SECTION("Non-axis-aligned cell: a₁ = (1,0,0), a₂ = (0.5,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 1.0;
    a2 = 0.0;
    a2(0) = 0.5;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * M_PI, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
  }
}

// =============================================================================
// BZ wrapping offset: FloquetPortData::ComputeBZOffset
// =============================================================================
TEST_CASE("Floquet BZ Wrapping Offset", "[floquetport][Serial]")
{
  // Period a = 4 in x → b₁ = 2π/4 = π/2 in x. BZ: [-π/4, π/4].
  mfem::Vector a1(3), a2(3), b1(3), b2(3);
  a1 = 0.0;
  a1(0) = 4.0;
  a2 = 0.0;
  a2(1) = 1.0;
  FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);
  double b1_sq = b1 * b1;

  SECTION("k_F inside BZ: no wrapping")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = 0.5;
    kF_w = kF;  // No wrapping needed

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == 0);
  }

  SECTION("k_F outside BZ: wrapping by +1 b₁")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = 1.05;  // Outside π/4 ≈ 0.785
    kF_w = 0.0;
    kF_w(0) = std::remainder(1.05, 2.0 * M_PI / 4.0);  // Wrapped

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == 1);
    // Verify roundtrip: kF_w + bz_m * b₁ ≈ kF
    CHECK_THAT(kF_w(0) + bz_m * b1(0), WithinRel(kF(0), 1e-12));
  }

  SECTION("Negative k_F outside BZ: wrapping by -1 b₁")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = -1.05;
    kF_w = 0.0;
    kF_w(0) = std::remainder(-1.05, 2.0 * M_PI / 4.0);

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == -1);
    CHECK_THAT(kF_w(0) + bz_m * b1(0), WithinRel(kF(0), 1e-12));
  }
}

}  // namespace palace
