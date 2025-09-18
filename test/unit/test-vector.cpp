#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <cmath>

#include "linalg/vector.hpp"

namespace palace
{
using namespace Catch;

TEST_CASE("StaticVectorConstruction", "[vector]")
{
  StaticVector<3> vec;
  REQUIRE(vec.Size() == 3);
  REQUIRE(vec.GetData() != nullptr);
}

TEST_CASE("StaticVectorElementAccess", "[vector]")
{
  StaticVector<3> vec;
  vec[0] = 1.0;
  vec[1] = 2.0;
  vec[2] = 3.0;

  REQUIRE(vec[0] == Approx(1.0));
  REQUIRE(vec[1] == Approx(2.0));
  REQUIRE(vec[2] == Approx(3.0));
}

TEST_CASE("StaticVectorInterface", "[vector]")
{
  // Test Inheritance from mfem::Vector

  StaticVector<4> vec;
  vec = 5.0;  // Set all elements

  REQUIRE(vec[0] == Approx(5.0));
  REQUIRE(vec[1] == Approx(5.0));
  REQUIRE(vec[2] == Approx(5.0));
  REQUIRE(vec[3] == Approx(5.0));
}

TEST_CASE("StaticVectorSizes", "[vector]")
{
  StaticVector<1> vec1;
  StaticVector<10> vec10;
  StaticVector<100> vec100;

  REQUIRE(vec1.Size() == 1);
  REQUIRE(vec10.Size() == 10);
  REQUIRE(vec100.Size() == 100);
}

TEST_CASE("StaticVectorPolymorphism", "[vector]")
{
  StaticVector<3> static_vec;
  Vector &vec_ref = static_vec;  // Polymorphic usage

  vec_ref[0] = 42.0;
  vec_ref[1] = -3.14;
  vec_ref[2] = 0.0;

  REQUIRE(static_vec[0] == Approx(42.0));
  REQUIRE(static_vec[1] == Approx(-3.14));
  REQUIRE(static_vec[2] == Approx(0.0));
  REQUIRE(vec_ref.Size() == 3);
}

TEST_CASE("StaticVectorOperations", "[vector]")
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
  REQUIRE(dot == Approx(32.0));  // 1*4 + 2*5 + 3*6 = 32

  // Test norm
  double norm1 = vec1.Norml2();
  REQUIRE(norm1 == Approx(std::sqrt(14.0)));  // sqrt(1^2 + 2^2 + 3^2)
}

TEST_CASE("StaticVector", "[vector]")
{
  SECTION("Basic StaticVector operations")
  {
    StaticVector<3> vec;
    vec(0) = 1.0;
    vec(1) = 2.0;
    vec(2) = 3.0;

    REQUIRE(vec(0) == Approx(1.0));
    REQUIRE(vec(1) == Approx(2.0));
    REQUIRE(vec(2) == Approx(3.0));
    REQUIRE(vec.Size() == 3);
  }

  SECTION("StaticVector SetSize should abort")
  {
    StaticVector<3> vec;

    // This should abort/throw
    REQUIRE_THROWS([&]() { vec.SetSize(5); }());
  }
}

TEST_CASE("Cross3 function", "[vector]")
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

    REQUIRE(C(0) == Approx(0.0));
    REQUIRE(C(1) == Approx(0.0));
    REQUIRE(C(2) == Approx(1.0));
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
    REQUIRE(C(0) == Approx(-3.0));
    REQUIRE(C(1) == Approx(6.0));
    REQUIRE(C(2) == Approx(-3.0));
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

    REQUIRE(C(0) == Approx(1.0));  // 1.0 + 0.0
    REQUIRE(C(1) == Approx(2.0));  // 2.0 + 0.0
    REQUIRE(C(2) == Approx(4.0));  // 3.0 + 1.0
  }

  SECTION("std::vector cross product")
  {
    std::vector<double> A = {1.0, 2.0, 3.0};
    std::vector<double> B = {4.0, 5.0, 6.0};
    std::vector<double> C(3);

    palace::linalg::Cross3(A, B, C);

    // Expected: A × B = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3)
    REQUIRE(C[0] == Approx(-3.0));
    REQUIRE(C[1] == Approx(6.0));
    REQUIRE(C[2] == Approx(-3.0));
  }

  SECTION("std::array cross product")
  {
    std::array<double, 3> A = {0.0, 0.0, 1.0};
    std::array<double, 3> B = {1.0, 0.0, 0.0};
    std::array<double, 3> C;

    palace::linalg::Cross3(A, B, C);

    // Expected: A × B = (0*0 - 1*0, 1*1 - 0*0, 0*0 - 0*1) = (0, 1, 0)
    REQUIRE(C[0] == Approx(0.0));
    REQUIRE(C[1] == Approx(1.0));
    REQUIRE(C[2] == Approx(0.0));
  }
}

TEST_CASE("Sqrt function", "[vector]")
{
  SECTION("Basic square root")
  {
    Vector vec(4);
    vec(0) = 4.0;
    vec(1) = 9.0;
    vec(2) = 16.0;
    vec(3) = 25.0;

    palace::linalg::Sqrt(vec);

    REQUIRE(vec(0) == Approx(2.0));
    REQUIRE(vec(1) == Approx(3.0));
    REQUIRE(vec(2) == Approx(4.0));
    REQUIRE(vec(3) == Approx(5.0));
  }

  SECTION("Square root with scaling")
  {
    Vector vec(3);
    vec(0) = 1.0;
    vec(1) = 4.0;
    vec(2) = 9.0;

    palace::linalg::Sqrt(vec, 4.0);  // sqrt(4 * x)

    REQUIRE(vec(0) == Approx(2.0));  // sqrt(4 * 1) = 2
    REQUIRE(vec(1) == Approx(4.0));  // sqrt(4 * 4) = 4
    REQUIRE(vec(2) == Approx(6.0));  // sqrt(4 * 9) = 6
  }
}
}  // namespace palace
