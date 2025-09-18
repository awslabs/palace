#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <complex>

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch;

TEST_CASE("Vector Sum - Real", "[Vector][Serial][Parallel][GPU]")
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

  REQUIRE(sum == Approx(expected));
}

TEST_CASE("Vector Sum - Real - Asymmetric Sizes", "[Vector][Parallel][GPU]")
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

  REQUIRE(sum == Approx(expected));
}

TEST_CASE("Vector Sum - Complex", "[Vector][Serial][Parallel][GPU]")
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

  REQUIRE(sum.real() == Approx(expected_real));
  REQUIRE(sum.imag() == Approx(expected_imag));
}
}  // namespace palace
