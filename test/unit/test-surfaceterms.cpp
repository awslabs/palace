// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceterms.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// The frequency-dependent boundary mass terms of a unit cube (boundary attributes 1 to 6):
// a finite conductivity boundary (attributes 1, 2), a series RL rational impedance boundary
// (attributes 3, 4, with the impedance split between the faces of attribute 4), and a
// parallel RC one (attribute 5), nondimensional.
class CubeSurfaceTerms
{
private:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<MaterialOperator> mat_op;
  std::unique_ptr<SurfaceTerms> terms;
  Units units;

public:
  CubeSurfaceTerms(const std::vector<config::RationalImpedanceData> &rational_impedance)
    : units(1.0, 1.0)
  {
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
    REQUIRE(Mpi::Size(Mpi::World()) <= serial_mesh->GetNE());
    mesh =
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh));
    config::MaterialData material;
    material.attributes = {1};
    config::PeriodicBoundaryData periodic;
    mat_op = std::make_unique<MaterialOperator>(std::vector<config::MaterialData>{material},
                                                periodic, ProblemType::DRIVEN, *mesh);
    config::BoundaryData boundaries;
    auto &cond = boundaries.conductivity.emplace_back();
    cond.sigma = 1.0e3;
    cond.h = 0.1;
    cond.attributes = {1, 2};
    boundaries.rational_impedance = rational_impedance;
    boundaries.split_impedance_attributes = {4};
    terms = std::make_unique<SurfaceTerms>(boundaries, ProblemType::DRIVEN, units, *mat_op,
                                           mesh->Get());
  }

  SurfaceTerms &Get() { return *terms; }
  MaterialPropertyCoefficient Coefficient() const
  {
    return MaterialPropertyCoefficient(mat_op->MaxCeedBdrAttribute());
  }

  // Whether boundary attribute attr has an entry in the coefficient; the attributes are not
  // on all processes.
  bool Has(const MaterialPropertyCoefficient &fb, int attr) const
  {
    bool has = false;
    auto ceed_attrs = mesh->GetCeedBdrAttributes(attr);
    if (ceed_attrs.Size() > 0)
    {
      has = (fb.GetAttributeToMaterial()[ceed_attrs[0] - 1] >= 0);
    }
    Mpi::GlobalOr(1, &has, Mpi::World());
    return has;
  }

  // Coefficient of boundary attribute attr, or 0 when it has none; the attributes are not
  // on all processes.
  double Value(const MaterialPropertyCoefficient &fb, int attr) const
  {
    double value = 0.0;
    auto ceed_attrs = mesh->GetCeedBdrAttributes(attr);
    if (ceed_attrs.Size() > 0)
    {
      const int k = fb.GetAttributeToMaterial()[ceed_attrs[0] - 1];
      value = (k < 0) ? 0.0 : fb.GetMaterialProperties()(0, 0, k);
    }
    double values[2] = {value, -value};
    Mpi::GlobalMax(2, values, Mpi::World());
    return (values[0] != 0.0) ? values[0] : -values[1];
  }
};

config::RationalImpedanceData RationalImpedance(const std::vector<int> &attributes,
                                                const std::vector<double> &num,
                                                const std::vector<double> &den)
{
  config::RationalImpedanceData data;
  data.attributes = attributes;
  data.num = num;
  data.den = den;
  return data;
}

const std::vector<config::RationalImpedanceData> rational_impedance = {
    RationalImpedance({3, 4}, {0.3, 2.0}, {1.0}),
    RationalImpedance({5}, {4.0}, {2.0, 1.0})};

}  // namespace

TEST_CASE("SurfaceTerms boundary mass coefficients are the sum of the terms",
          "[surfaceterms][Serial][Parallel]")
{
  // The frequency-dependent boundary mass coefficients which SurfaceTerms adds at the
  // frequency ω are, on each boundary attribute, the sum over its terms of the unit
  // coefficient times the scale at ω: the conductivity group, then the two rational
  // impedance boundaries, each on its own attributes. Each attribute of a term has an entry
  // in both parts, and the attribute without a term (6) in neither.
  CubeSurfaceTerms cube(rational_impedance);
  auto &terms = cube.Get();
  REQUIRE(terms.GetNumExtraSystemBdrMassTerms() == 3);
  const std::vector<std::vector<int>> term_attrs = {{1, 2}, {3, 4}, {5}};
  for (std::complex<double> omega :
       {std::complex<double>(0.7, 0.0), std::complex<double>(2.5, 0.0),
        std::complex<double>(1.5, 0.4)})
  {
    CAPTURE(omega);
    auto fbr = cube.Coefficient(), fbi = cube.Coefficient();
    if (omega.imag() == 0.0)
    {
      terms.AddExtraSystemBdrMassCoefficients(omega.real(), fbr, fbi);
    }
    else
    {
      auto dfbr = cube.Coefficient(), dfbi = cube.Coefficient(), frbr = cube.Coefficient(),
           frbi = cube.Coefficient();
      terms.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, nullptr,
                                          nullptr);
    }
    for (int i = 0; i < terms.GetNumExtraSystemBdrMassTerms(); i++)
    {
      auto unit = cube.Coefficient();
      terms.AddExtraSystemBdrMassTermCoefficients(i, unit);
      const std::complex<double> scale = terms.GetExtraSystemBdrMassTermScale(i, omega);
      REQUIRE(std::abs(scale.real()) > 0.0);
      REQUIRE(std::abs(scale.imag()) > 0.0);
      for (int attr : term_attrs[i])
      {
        CAPTURE(i, attr);
        REQUIRE(cube.Has(unit, attr));
        const double u = cube.Value(unit, attr);
        REQUIRE(u > 0.0);
        CHECK(cube.Has(fbr, attr));
        CHECK(cube.Has(fbi, attr));
        CHECK_THAT(cube.Value(fbr, attr), WithinRel(scale.real() * u, 1.0e-12));
        CHECK_THAT(cube.Value(fbi, attr), WithinRel(scale.imag() * u, 1.0e-12));
      }
    }
    CHECK_FALSE(cube.Has(fbr, 6));
    CHECK_FALSE(cube.Has(fbi, 6));
  }
}

TEST_CASE("SurfaceTerms exact split of a boundary mass term",
          "[surfaceterms][Serial][Parallel]")
{
  // A rational impedance boundary term has an exact split of its scale c(s), s = iω, into a
  // polynomial P and a strictly proper remainder R, c = P + R, and a conductivity group
  // term none. Here Zs(s) = 0.3 s + 2 (series RL: c(s) = s/(0.3 s + 2), deg P = 0), Zs(s) =
  // 4 / (2 s + 1) (parallel RC: c(s) = 0.5 s² + 0.25 s, no remainder, deg P = 2), and Zs(s)
  // = 1 / (s² + 0.5 s + 2) (c(s) = s³ + 0.5 s² + 2 s, deg P = 3, from an impedance which is
  // not passive).
  CubeSurfaceTerms cube({rational_impedance[0], rational_impedance[1],
                         RationalImpedance({6}, {1.0}, {1.0, 0.5, 2.0})});
  auto &terms = cube.Get();
  REQUIRE(terms.GetNumExtraSystemBdrMassTerms() == 4);
  CHECK_FALSE(terms.GetExtraSystemBdrMassTermSplit(0));
  const std::vector<std::vector<double>> polynomials = {
      {1.0 / 0.3}, {0.5, 0.25, 0.0}, {1.0, 0.5, 2.0, 0.0}};
  const std::vector<int> degrees = {0, 2, 3};
  for (int i = 1; i < 4; i++)
  {
    CAPTURE(i);
    const auto split = terms.GetExtraSystemBdrMassTermSplit(i);
    REQUIRE(split);
    CHECK(split->degree == degrees[i - 1]);
    REQUIRE(split->polynomial.size() == polynomials[i - 1].size());
    for (std::size_t k = 0; k < split->polynomial.size(); k++)
    {
      CAPTURE(k);
      const double p = polynomials[i - 1][k];
      CHECK_THAT(split->polynomial[k], WithinRel(p, 1.0e-14) || WithinAbs(p, 1.0e-15));
    }
    for (std::complex<double> omega :
         {std::complex<double>(0.7, 0.0), std::complex<double>(2.5, 0.0),
          std::complex<double>(1.5, 0.4)})
    {
      CAPTURE(omega);
      const std::complex<double> s = std::complex<double>(0.0, 1.0) * omega;
      const std::complex<double> c = split->coefficient(s);
      CHECK(std::abs(c - terms.GetExtraSystemBdrMassTermScale(i, omega)) <=
            1.0e-14 * std::abs(c));
      std::complex<double> P = 0.0;
      for (double p : split->polynomial)
      {
        P = P * s + p;
      }
      const std::complex<double> R = split->remainder(s);
      CHECK(std::abs(P + R - c) <= 1.0e-13 * std::abs(c));
      if (i == 1)
      {
        CHECK(std::abs(R) > 1.0e-3 * std::abs(c));
      }
      else
      {
        CHECK(std::abs(R) <= 1.0e-14 * std::abs(c));
      }
    }
  }
}

}  // namespace palace
