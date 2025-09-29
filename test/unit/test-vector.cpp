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

  REQUIRE_THAT(sum, Catch::Matchers::WithinRel(expected));
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

  REQUIRE_THAT(sum, Catch::Matchers::WithinRel(expected));
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

  REQUIRE_THAT(sum.real(), Catch::Matchers::WithinRel(expected_real));
  REQUIRE_THAT(sum.imag(), Catch::Matchers::WithinRel(expected_imag));
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

  REQUIRE_THAT(cv.Real()[0], Catch::Matchers::WithinRel(1.0 * rank));
  REQUIRE_THAT(cv.Real()[1], Catch::Matchers::WithinRel(2.0 * rank));

  REQUIRE_THAT(cv.Imag()[0], Catch::Matchers::WithinRel(10.0 * rank));
  REQUIRE_THAT(cv.Imag()[1], Catch::Matchers::WithinRel(20.0 * rank));
}

}  // namespace palace
