#include <iterator>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/romoperator.hpp"

using namespace palace;

// auto complex_circle_sample_points(int nr_sample_points, double radius = 5.5)
// {
//   double end_point_linscale = double(nr_sample_points - 1) / nr_sample_points;
//   Eigen::ArrayXcd zj_sample =
//       Eigen::ArrayXcd::LinSpaced(nr_sample_points, 0, 2 * M_PI * end_point_linscale);
//   zj_sample =
//       zj_sample.unaryExpr([radius](std::complex<double> z)
//                           { return radius * std::exp(std::complex<double>(0., 1.) * z);
//                           });
//   return zj_sample;
// }

TEST_CASE("MinimalRationalInterpolation", "[romoperator]")
{
  auto *comm = Mpi::World();

  auto fn_tan_shift = [](double z)
  { return std::tan(0.5 * M_PI * (z - std::complex<double>(1., 1.))); };

  // Test scalar case: 2 sample points for 2 x 2 vector
  MinimalRationalInterpolation mri_1(2);

  CHECK(mri_1.GetSamplePoints() == std::vector<double>{});
  CHECK_THROWS(mri_1.FindMaxError(1));

  for (double x_sample : {-1.0, 1.0})
  {
    ComplexVector c_vec(1);
    c_vec = fn_tan_shift(x_sample) / double(Mpi::Size(comm));

    mri_1.AddSolutionSample(x_sample, c_vec, comm, Orthogonalization::MGS);
  }

  CHECK(mri_1.GetSamplePoints().size() == 2);
  CHECK(mri_1.GetSamplePoints() == std::vector<double>{-1.0, 1.0});
  // By symmetry of poles max erro should be at zero.
  auto max_err_1 = mri_1.FindMaxError(1);
  REQUIRE(max_err_1.size() == 1);
  CHECK_THAT(max_err_1[0], Catch::Matchers::WithinAbsMatcher(0., 1e-6));
}