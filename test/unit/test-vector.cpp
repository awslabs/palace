// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <complex>

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("Vector Sum - Real", "[vector][Serial][Parallel][GPU]")
{
  // Create a Vector (on the device, when available), populate it with MPI-dependent known
  // values, and check that the resulting Sum is what we expect.

  int rank = Mpi::Rank(Mpi::World());
  int size = Mpi::Size(Mpi::World());

  Vector v(3);
  v.UseDevice(true);
  auto d_v = v.Write();

  // Rank 0: [1, 2, 3], Rank 1: [4, 5, 6], etc.
  mfem::forall(v.Size(), [=] MFEM_HOST_DEVICE(int i) { d_v[i] = rank * 3 + i + 1; });

  double sum = linalg::Sum(Mpi::World(), v);

  // Expected sum: sum over all ranks of (rank * 3 + 1) + (rank * 3 + 2) +
  // (rank * 3 + 3) = sum over ranks of (3 * rank * 3 + 6) = 9 * sum(rank) + 6
  // * size sum(rank) from 0 to size-1 = size * (size-1) / 2.
  double expected = 9.0 * size * (size - 1) / 2.0 + 6.0 * size;

  CHECK_THAT(sum, WithinRel(expected));
}

TEST_CASE("Vector Sum - Real - Asymmetric Sizes", "[vector][Parallel][GPU]")
{
  // Test that Vector Sum works as expected when different MPI ranks hold different amount
  // of data.
  int rank = Mpi::Rank(Mpi::World());
  int size = Mpi::Size(Mpi::World());
  int local_size = (rank == 0) ? 4 : 2;  // Rank 0 has more elements
  Vector v(local_size);

  v.UseDevice(true);
  auto d_v = v.Write();

  // Fill with rank-dependent values.
  mfem::forall(local_size,
               [=] MFEM_HOST_DEVICE(int i)
               {
                 d_v[i] =
                     (rank + 1) * 10 + i;  // Rank 0: [10,11,12,13], Rank 1: [20,21], etc.
               });

  double sum = linalg::Sum(Mpi::World(), v);

  // Calculate expected sum manually.
  double expected = 0.0;
  for (int r = 0; r < size; r++)
  {
    int r_size = (r == 0) ? 4 : 2;
    for (int i = 0; i < r_size; i++)
      expected += (r + 1) * 10 + i;
  }

  CHECK_THAT(sum, WithinRel(expected));
}

TEST_CASE("Vector Sum - Complex", "[vector][Serial][Parallel][GPU]")
{
  // Create a ComplexVector with rank-dependent data and verify that Sum for both the real
  // and imaginary components is correct.
  int rank = Mpi::Rank(Mpi::World());
  int size = Mpi::Size(Mpi::World());

  ComplexVector cv(2);
  cv.UseDevice(true);
  auto d_real = cv.Real().Write();
  auto d_imag = cv.Imag().Write();

  mfem::forall(cv.Size(),
               [=] MFEM_HOST_DEVICE(int i)
               {
                 // Rank 0: [(0, 1), (0, 2), ...], Rank 1: [(1, 2), (1, 3), ...], etc.
                 d_real[i] = rank;
                 d_imag[i] = rank + i;
               });

  std::complex<double> sum = linalg::Sum(Mpi::World(), cv);

  // Real part: each rank contributes rank * cv.Size().
  // Total = sum over all ranks of (rank * cv.Size()) = cv.Size() * sum(ranks).
  double rank_sum = size * (size - 1) / 2.0;
  double expected_real = cv.Size() * rank_sum;

  // Imaginary part: each rank contributes rank * cv.Size() + sum(0 to cv.Size()-1).
  // sum(0 to cv.Size()-1) = cv.Size() * (cv.Size() - 1) / 2.
  double index_sum = cv.Size() * (cv.Size() - 1) / 2.0;
  double expected_imag = cv.Size() * rank_sum + size * index_sum;

  CHECK_THAT(sum.real(), WithinRel(expected_real));
  CHECK_THAT(sum.imag(), WithinRel(expected_imag));
}

TEST_CASE("ComplexVector Set", "[vector][Serial][Parallel][GPU]")
{
  ComplexVector cv(2);

  // Set requires explicitely specifying on_dev depending on the device.
  bool on_dev = mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  cv.UseDevice(on_dev);

  int rank = Mpi::Rank(Mpi::World());
  std::complex<double> vals[2];

  vals[0] = std::complex<double>(1.0 * rank, 10.0 * rank);
  vals[1] = std::complex<double>(2.0 * rank, 20.0 * rank);

  cv.Set(vals, 2, on_dev);

  CHECK_THAT(cv.Real()[0], WithinRel(1.0 * rank) || WithinAbs(0.0, 1e-14));
  CHECK_THAT(cv.Real()[1], WithinRel(2.0 * rank) || WithinAbs(0.0, 1e-14));

  CHECK_THAT(cv.Imag()[0], WithinRel(10.0 * rank) || WithinAbs(0.0, 1e-14));
  CHECK_THAT(cv.Imag()[1], WithinRel(20.0 * rank) || WithinAbs(0.0, 1e-14));
}

TEST_CASE("StaticVectorConstruction", "[Vector][Serial]")
{
  StaticVector<3> vec;
  CHECK(vec.Size() == 3);
  CHECK(vec.GetData() != nullptr);
}

TEST_CASE("StaticVectorElementAccess", "[Vector][Serial]")
{
  StaticVector<3> vec;
  vec[0] = 1.0;
  vec[1] = 2.0;
  vec[2] = 3.0;

  CHECK_THAT(vec[0], WithinRel(1.0));
  CHECK_THAT(vec[1], WithinRel(2.0));
  CHECK_THAT(vec[2], WithinRel(3.0));
}

TEST_CASE("StaticVectorInterface", "[Vector][Serial]")
{
  // Test Inheritance from mfem::Vector

  StaticVector<4> vec;
  vec = 5.0;  // Set all elements

  CHECK_THAT(vec[0], WithinRel(5.0));
  CHECK_THAT(vec[1], WithinRel(5.0));
  CHECK_THAT(vec[2], WithinRel(5.0));
  CHECK_THAT(vec[3], WithinRel(5.0));
}

TEST_CASE("StaticVectorSizes", "[Vector][Serial]")
{
  StaticVector<1> vec1;
  StaticVector<10> vec10;
  StaticVector<100> vec100;

  CHECK(vec1.Size() == 1);
  CHECK(vec10.Size() == 10);
  CHECK(vec100.Size() == 100);
}

TEST_CASE("StaticVectorPolymorphism", "[Vector][Serial]")
{
  StaticVector<3> static_vec;
  Vector &vec_ref = static_vec;  // Polymorphic usage

  vec_ref[0] = 42.0;
  vec_ref[1] = -3.14;
  vec_ref[2] = 0.0;

  CHECK_THAT(static_vec[0], WithinRel(42.0));
  CHECK_THAT(static_vec[1], WithinRel(-3.14));
  CHECK_THAT(static_vec[2], WithinRel(0.0));
  CHECK(vec_ref.Size() == 3);
}

TEST_CASE("StaticVectorOperations", "[Vector][Serial]")
{
  StaticVector<3> vec1, vec2;

  vec1[0] = 1.0;
  vec1[1] = 2.0;
  vec1[2] = 3.0;
  vec2[0] = 4.0;
  vec2[1] = 5.0;
  vec2[2] = 6.0;

  // Test dot product
  double dot = vec1 * vec2;
  CHECK_THAT(dot, WithinRel(32.0));  // 1*4 + 2*5 + 3*6 = 32

  // Test norm
  double norm1 = vec1.Norml2();
  CHECK_THAT(norm1, WithinRel(std::sqrt(14.0)));  // sqrt(1^2 + 2^2 + 3^2)
}

TEST_CASE("StaticVector", "[Vector][Serial]")
{
  SECTION("Basic StaticVector operations")
  {
    StaticVector<3> vec;
    vec(0) = 1.0;
    vec(1) = 2.0;
    vec(2) = 3.0;

    CHECK_THAT(vec(0), WithinRel(1.0));
    CHECK_THAT(vec(1), WithinRel(2.0));
    CHECK_THAT(vec(2), WithinRel(3.0));
    CHECK(vec.Size() == 3);
  }
}

TEST_CASE("Cross3 function", "[Vector][Serial]")
{
  SECTION("MFEM Vector cross product")
  {
    mfem::Vector A(3), B(3), C(3);
    A(0) = 1.0;
    A(1) = 0.0;
    A(2) = 0.0;
    B(0) = 0.0;
    B(1) = 1.0;
    B(2) = 0.0;

    palace::linalg::Cross3(A, B, C);

    CHECK_THAT(C(0), WithinRel(0.0));
    CHECK_THAT(C(1), WithinRel(0.0));
    CHECK_THAT(C(2), WithinRel(1.0));
  }

  SECTION("MFEM Vector general cross product")
  {
    mfem::Vector A(3), B(3), C(3);
    A(0) = 2.0;
    A(1) = 3.0;
    A(2) = 4.0;
    B(0) = 5.0;
    B(1) = 6.0;
    B(2) = 7.0;

    palace::linalg::Cross3(A, B, C);

    // Expected: A × B = (3*7 - 4*6, 4*5 - 2*7, 2*6 - 3*5) = (-3, 6, -3)
    CHECK_THAT(C(0), WithinRel(-3.0));
    CHECK_THAT(C(1), WithinRel(6.0));
    CHECK_THAT(C(2), WithinRel(-3.0));
  }

  SECTION("MFEM Vector with add=true")
  {
    mfem::Vector A(3), B(3), C(3);
    A(0) = 1.0;
    A(1) = 0.0;
    A(2) = 0.0;
    B(0) = 0.0;
    B(1) = 1.0;
    B(2) = 0.0;
    C(0) = 1.0;
    C(1) = 2.0;
    C(2) = 3.0;

    palace::linalg::Cross3(A, B, C, true);

    CHECK_THAT(C(0), WithinRel(1.0));  // 1.0 + 0.0
    CHECK_THAT(C(1), WithinRel(2.0));  // 2.0 + 0.0
    CHECK_THAT(C(2), WithinRel(4.0));  // 3.0 + 1.0
  }

  SECTION("std::vector cross product")
  {
    std::vector<double> A = {1.0, 2.0, 3.0};
    std::vector<double> B = {4.0, 5.0, 6.0};
    std::vector<double> C(3);

    palace::linalg::Cross3(A, B, C);

    // Expected: A × B = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
    CHECK_THAT(C[0], WithinRel(-3.0));
    CHECK_THAT(C[1], WithinRel(6.0));
    CHECK_THAT(C[2], WithinRel(-3.0));
  }

  SECTION("std::array cross product")
  {
    std::array<double, 3> A = {0.0, 0.0, 1.0};
    std::array<double, 3> B = {1.0, 0.0, 0.0};
    std::array<double, 3> C;

    palace::linalg::Cross3(A, B, C);

    // Expected: A × B = (0*0 - 1*0, 1*1 - 0*0, 0*0 - 0*1) = (0, 1, 0)
    CHECK_THAT(C[0], WithinRel(0.0));
    CHECK_THAT(C[1], WithinRel(1.0));
    CHECK_THAT(C[2], WithinRel(0.0));
  }
}

TEST_CASE("Sqrt function", "[Vector][Serial]")
{
  SECTION("Basic square root")
  {
    Vector vec(4);
    vec(0) = 4.0;
    vec(1) = 9.0;
    vec(2) = 16.0;
    vec(3) = 25.0;

    palace::linalg::Sqrt(vec);

    CHECK_THAT(vec(0), WithinRel(2.0));
    CHECK_THAT(vec(1), WithinRel(3.0));
    CHECK_THAT(vec(2), WithinRel(4.0));
    CHECK_THAT(vec(3), WithinRel(5.0));
  }

  SECTION("Square root with scaling")
  {
    Vector vec(3);
    vec(0) = 1.0;
    vec(1) = 4.0;
    vec(2) = 9.0;

    palace::linalg::Sqrt(vec, 4.0);  // sqrt(4 * x)

    CHECK_THAT(vec(0), WithinRel(2.0));  // sqrt(4 * 1) = 2
    CHECK_THAT(vec(1), WithinRel(4.0));  // sqrt(4 * 4) = 4
    CHECK_THAT(vec(2), WithinRel(6.0));  // sqrt(4 * 9) = 6
  }
}
}  // namespace palace
