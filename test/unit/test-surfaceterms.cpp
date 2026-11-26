// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
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

// Whether boundary attribute attr has an entry in the coefficient; the attributes are not
// on all processes.
bool HasCoefficient(const Mesh &mesh, const MaterialPropertyCoefficient &fb, int attr)
{
  bool has = false;
  auto ceed_attrs = mesh.GetCeedBdrAttributes(attr);
  if (ceed_attrs.Size() > 0)
  {
    has = (fb.GetAttributeToMaterial()[ceed_attrs[0] - 1] >= 0);
  }
  Mpi::GlobalOr(1, &has, Mpi::World());
  return has;
}

// Coefficient of boundary attribute attr, or 0 when it has none; the attributes are not on
// all processes.
double CoefficientValue(const Mesh &mesh, const MaterialPropertyCoefficient &fb, int attr)
{
  double value = 0.0;
  auto ceed_attrs = mesh.GetCeedBdrAttributes(attr);
  if (ceed_attrs.Size() > 0)
  {
    const int k = fb.GetAttributeToMaterial()[ceed_attrs[0] - 1];
    value = (k < 0) ? 0.0 : fb.GetMaterialProperties()(0, 0, k);
  }
  double values[2] = {value, -value};
  Mpi::GlobalMax(2, values, Mpi::World());
  return (values[0] != 0.0) ? values[0] : -values[1];
}

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

  bool Has(const MaterialPropertyCoefficient &fb, int attr) const
  {
    return HasCoefficient(*mesh, fb, attr);
  }
  double Value(const MaterialPropertyCoefficient &fb, int attr) const
  {
    return CoefficientValue(*mesh, fb, attr);
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
  // in both parts, and the attribute without a term (6) in neither. This holds for the
  // coefficients of a vector field and, at real ω, of a scalar field (with the unit
  // coefficients of each).
  CubeSurfaceTerms cube(rational_impedance);
  auto &terms = cube.Get();
  REQUIRE(terms.GetNumExtraSystemBdrMassTerms() == 3);
  const std::vector<std::vector<int>> term_attrs = {{1, 2}, {3, 4}, {5}};
  for (std::complex<double> omega :
       {std::complex<double>(0.7, 0.0), std::complex<double>(2.5, 0.0),
        std::complex<double>(1.5, 0.4)})
  {
    CAPTURE(omega);
    auto fbr = cube.Coefficient(), fbi = cube.Coefficient(), fbr_n = cube.Coefficient(),
         fbi_n = cube.Coefficient();
    const bool real = (omega.imag() == 0.0);
    if (real)
    {
      terms.AddExtraSystemBdrMassCoefficients(omega.real(), fbr, fbi);
      terms.AddScalarExtraSystemBdrMassCoefficients(omega.real(), fbr_n, fbi_n);
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
      auto unit = cube.Coefficient(), unit_n = cube.Coefficient();
      terms.AddExtraSystemBdrMassTermCoefficients(i, unit);
      terms.AddScalarExtraSystemBdrMassTermCoefficients(i, unit_n);
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
        if (real)
        {
          REQUIRE(cube.Has(unit_n, attr));
          const double u_n = cube.Value(unit_n, attr);
          REQUIRE(u_n > 0.0);
          CHECK(cube.Has(fbr_n, attr));
          CHECK(cube.Has(fbi_n, attr));
          CHECK_THAT(cube.Value(fbr_n, attr), WithinRel(scale.real() * u_n, 1.0e-12));
          CHECK_THAT(cube.Value(fbi_n, attr), WithinRel(scale.imag() * u_n, 1.0e-12));
        }
      }
    }
    CHECK_FALSE(cube.Has(fbr, 6));
    CHECK_FALSE(cube.Has(fbi, 6));
    CHECK_FALSE(cube.Has(fbr_n, 6));
    CHECK_FALSE(cube.Has(fbi_n, 6));
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

TEST_CASE("SurfaceTerms passivity check of a rational impedance",
          "[surfaceterms][Serial][Parallel]")
{
  // A rational impedance boundary which is not passive at a real frequency (Re{Zs} < 0)
  // warns once: its passivity check returns true the first time and false after, and it
  // is false for a passive boundary. The summed methods check each boundary at a real
  // frequency, but not when only the scale is evaluated or at a complex frequency. Not
  // passive: Zs(s) = 0.3 s - 1 on the attributes 3, 4 (Re{Zs} = -1) and
  // Zs(s) = 4 / (2 s - 1) on the attribute 5 (Re{Zs} = -4 / (1 + 4 ω²)). Passive:
  // Zs(s) = 4 / (2 s + 1) on the attribute 6.
  const std::vector<config::RationalImpedanceData> impedance = {
      RationalImpedance({3, 4}, {0.3, -1.0}, {1.0}),
      RationalImpedance({5}, {4.0}, {2.0, -1.0}),
      RationalImpedance({6}, {4.0}, {2.0, 1.0})};
  const std::vector<int> not_passive = {0, 1};  // Boundary indices.
  const int passive = 2;
  CubeSurfaceTerms cube(impedance);
  auto &terms = cube.Get();
  auto &rz_op = terms.GetRationalImpedanceOp();
  REQUIRE(rz_op.GetNumBoundaries() == 3);
  auto fbr = cube.Coefficient(), fbi = cube.Coefficient(), dfbr = cube.Coefficient(),
       dfbi = cube.Coefficient(), frbr = cube.Coefficient(), frbi = cube.Coefficient();

  SECTION("Once for each boundary")
  {
    for (int b : not_passive)
    {
      CAPTURE(b);
      CHECK(rz_op.CheckPassivity(b, 0.7));
      CHECK_FALSE(rz_op.CheckPassivity(b, 0.7));
      CHECK_FALSE(rz_op.CheckPassivity(b, 1.2));
    }
    CHECK_FALSE(rz_op.CheckPassivity(passive, 0.7));
    CHECK_FALSE(rz_op.CheckPassivity(passive, 1.2));
  }

  SECTION("Checked by the summed methods at a real frequency")
  {
    const auto method = GENERATE(as<std::string>{}, "vector", "scalar", "3D");
    CAPTURE(method);
    if (method == "vector")
    {
      terms.AddExtraSystemBdrMassCoefficients(0.7, fbr, fbi);
    }
    else if (method == "scalar")
    {
      terms.AddScalarExtraSystemBdrMassCoefficients(0.7, fbr, fbi);
    }
    else
    {
      terms.AddExtraSystemBdrCoefficients(0.7, dfbr, dfbi, fbr, fbi, frbr, frbi, nullptr,
                                          nullptr);
    }
    for (int b : not_passive)
    {
      CAPTURE(b);
      CHECK_FALSE(rz_op.CheckPassivity(b, 1.2));
    }
  }

  SECTION("Not checked for the scale or at a complex frequency")
  {
    for (int i = 0; i < terms.GetNumExtraSystemBdrMassTerms(); i++)
    {
      terms.GetExtraSystemBdrMassTermScale(i, 0.7);
    }
    terms.AddExtraSystemBdrCoefficients(std::complex<double>(0.7, 0.1), dfbr, dfbi, fbr,
                                        fbi, frbr, frbi, nullptr, nullptr);
    for (int b : not_passive)
    {
      CAPTURE(b);
      CHECK(rz_op.CheckPassivity(b, 0.7));
    }
  }
}

TEST_CASE("SurfaceTerms absorbing boundary next to a frequency-dependent material",
          "[surfaceterms][Serial][Parallel]")
{
  // A damped Lorentz medium (ε∞ = 2) at x < 1/2 next to vacuum in the unit cube, with
  // absorbing faces at x = 0 (attribute 5, next to the medium) and x = 1 (attribute 3, next
  // to the vacuum). The first-order term next to the medium is the frequency-dependent
  // boundary mass term of the absorbing boundaries: its scale is iω r(ω), with
  // r = √(ε(iω)/ε∞), and its unit coefficient is the inverse impedance √ε∞ of the medium,
  // for both fields, on the attribute 5 only, which the damping terms leave out. The summed
  // methods of the system add the real and imaginary parts of the scale, and those of the
  // preconditioners its magnitude as a loss.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
  REQUIRE(Mpi::Size(Mpi::World()) <= serial_mesh->GetNE());
  for (int e = 0; e < serial_mesh->GetNE(); e++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(e, center);
    serial_mesh->SetAttribute(e, center(0) < 0.5 ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh));
  const double eps_inf = 2.0, delta = 0.7, w0 = 3.0, gamma = 0.4;
  config::MaterialData medium, vacuum;
  medium.attributes = {1};
  medium.epsilon_r.s = {eps_inf, eps_inf, eps_inf};
  medium.permittivity_terms.push_back(config::LorentzData{delta, w0, gamma});
  vacuum.attributes = {2};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({medium, vacuum}, periodic, ProblemType::DRIVEN, mesh);
  config::BoundaryData boundaries;
  boundaries.farfield.attributes = {3, 5};
  boundaries.farfield.order = 1;
  Units units(1.0, 1.0);
  SurfaceTerms terms(boundaries, ProblemType::DRIVEN, units, mat_op, mesh.Get());
  auto Coefficient = [&mat_op]()
  { return MaterialPropertyCoefficient(mat_op.MaxCeedBdrAttribute()); };
  auto Has = [&mesh](const MaterialPropertyCoefficient &fb, int attr)
  { return HasCoefficient(mesh, fb, attr); };
  auto Value = [&mesh](const MaterialPropertyCoefficient &fb, int attr)
  { return CoefficientValue(mesh, fb, attr); };

  REQUIRE(terms.GetNumExtraSystemBdrMassTerms() == 1);
  CHECK(terms.GetExtraSystemBdrMassTermLabel(0) == "farfield_permittivity_0");
  CHECK_FALSE(terms.GetExtraSystemBdrMassTermSplit(0));
  const double u = std::sqrt(eps_inf);
  {
    auto unit = Coefficient(), unit_n = Coefficient(), damp = Coefficient(),
         damp_n = Coefficient();
    terms.AddExtraSystemBdrMassTermCoefficients(0, unit);
    terms.AddScalarExtraSystemBdrMassTermCoefficients(0, unit_n);
    terms.AddDampingBdrCoefficients(1.0, damp);
    terms.AddScalarDampingBdrCoefficients(1.0, damp_n);
    for (const auto *fb : {&unit, &unit_n})
    {
      REQUIRE(Has(*fb, 5));
      CHECK_THAT(Value(*fb, 5), WithinRel(u, 1.0e-14));
      CHECK_FALSE(Has(*fb, 3));
    }
    for (const auto *fb : {&damp, &damp_n})
    {
      REQUIRE(Has(*fb, 3));
      CHECK_THAT(Value(*fb, 3), WithinRel(1.0, 1.0e-14));
      CHECK_FALSE(Has(*fb, 5));
    }
  }

  for (std::complex<double> omega :
       {std::complex<double>(1.3, 0.0), std::complex<double>(1.3, 0.2)})
  {
    CAPTURE(omega);
    const bool real = (omega.imag() == 0.0);
    const std::complex<double> s = std::complex<double>(0.0, 1.0) * omega;
    const std::complex<double> eps =
        eps_inf + delta * w0 * w0 / (s * s + gamma * s + w0 * w0);
    const std::complex<double> scale = s * std::sqrt(eps / eps_inf);
    CHECK(std::abs(terms.GetExtraSystemBdrMassTermScale(0, omega) - scale) <=
          1.0e-14 * std::abs(scale));
    REQUIRE(std::abs(scale.real()) > 0.0);
    REQUIRE(std::abs(scale.imag()) > 0.0);

    // The coefficients of the system, with the real and imaginary parts of the scale, and
    // of the preconditioners, with its magnitude as a loss, on the attribute 5 only.
    auto CheckSystem =
        [&](const MaterialPropertyCoefficient &fbr, const MaterialPropertyCoefficient &fbi)
    {
      CHECK_THAT(Value(fbr, 5), WithinRel(scale.real() * u, 1.0e-12));
      CHECK_THAT(Value(fbi, 5), WithinRel(scale.imag() * u, 1.0e-12));
      CHECK_FALSE(Has(fbr, 3));
      CHECK_FALSE(Has(fbi, 3));
    };
    auto CheckPreconditioner =
        [&](const MaterialPropertyCoefficient &fbr, const MaterialPropertyCoefficient &fbi)
    {
      CHECK_FALSE(Has(fbr, 5));
      CHECK_THAT(Value(fbi, 5), WithinRel(std::abs(scale) * u, 1.0e-12));
      CHECK_FALSE(Has(fbi, 3));
    };
    if (real)
    {
      for (bool scalar : {false, true})
      {
        CAPTURE(scalar);
        auto fbr = Coefficient(), fbi = Coefficient(), pr = Coefficient(),
             pi = Coefficient();
        if (scalar)
        {
          terms.AddScalarExtraSystemBdrMassCoefficients(omega.real(), fbr, fbi);
          terms.AddScalarPreconditionerExtraSystemBdrMassCoefficients(omega.real(), pr, pi);
        }
        else
        {
          terms.AddExtraSystemBdrMassCoefficients(omega.real(), fbr, fbi);
          terms.AddPreconditionerExtraSystemBdrMassCoefficients(omega.real(), pr, pi);
        }
        CheckSystem(fbr, fbi);
        CheckPreconditioner(pr, pi);
      }
    }
    auto dfbr = Coefficient(), dfbi = Coefficient(), fbr = Coefficient(),
         fbi = Coefficient(), frbr = Coefficient(), frbi = Coefficient();
    terms.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, nullptr,
                                        nullptr);
    CheckSystem(fbr, fbi);

    // A real-valued preconditioner of the 3D solve has one coefficient.
    auto dfb = Coefficient(), fb = Coefficient(), frb = Coefficient(),
         empty = Coefficient();
    terms.AddRealPreconditionerExtraSystemBdrCoefficients(omega, dfb, fb, frb);
    CheckPreconditioner(empty, fb);
  }
}

TEST_CASE("SurfaceTerms static limit of the boundary mass terms",
          "[surfaceterms][Serial][Parallel]")
{
  // The static limits of the frequency-dependent boundary mass terms, for the quasistatic
  // anchor solves of circuit synthesis: on each attribute, the unit coefficient of a term
  // times the limit of its scale c(s), s = iω, as s → 0, which is nonzero only for a
  // rational impedance with a simple zero at s = 0. Here an inductive sheet Zs(s) = 0.3 s
  // on the attributes 3, 4 (c(s) = 1/0.3, with the impedance split between the faces of
  // attribute 4), a parallel RL sheet Zs(s) = s / (0.5 s + 2) on the attribute 5 (c(0) =
  // 2), and a series RL sheet Zs(s) = 0.3 s + 2 on the attribute 6 (c(0) = 0), with the
  // finite conductivity boundary of the attributes 1, 2 (c(0) = 0).
  CubeSurfaceTerms cube({RationalImpedance({3, 4}, {0.3, 0.0}, {1.0}),
                         RationalImpedance({5}, {1.0, 0.0}, {0.5, 2.0}),
                         RationalImpedance({6}, {0.3, 2.0}, {1.0})});
  auto &terms = cube.Get();
  REQUIRE(terms.GetNumExtraSystemBdrMassTerms() == 4);
  auto fb = cube.Coefficient();
  terms.AddStaticExtraSystemBdrMassCoefficients(fb);
  const std::vector<std::vector<int>> term_attrs = {{1, 2}, {3, 4}, {5}, {6}};
  const std::vector<double> limits = {0.0, 1.0 / 0.3, 2.0, 0.0};
  for (int i = 0; i < terms.GetNumExtraSystemBdrMassTerms(); i++)
  {
    CAPTURE(i);
    if (i > 0)
    {
      const std::complex<double> c = terms.GetExtraSystemBdrMassTermScale(i, 1.0e-8);
      CHECK(std::abs(c - limits[i]) <= 1.0e-8);
    }
    auto unit = cube.Coefficient();
    terms.AddExtraSystemBdrMassTermCoefficients(i, unit);
    for (int attr : term_attrs[i])
    {
      CAPTURE(attr);
      if (limits[i] == 0.0)
      {
        CHECK_FALSE(cube.Has(fb, attr));
      }
      else
      {
        REQUIRE(cube.Has(fb, attr));
        CHECK_THAT(cube.Value(fb, attr),
                   WithinRel(limits[i] * cube.Value(unit, attr), 1.0e-14));
      }
    }
  }
  CHECK_THAT(cube.Value(fb, 4), WithinRel(0.5 * cube.Value(fb, 3), 1.0e-14));

  // The same parallel RL sheet with a common factor s, Zs(s) = s² / (0.5 s² + 2 s), and a
  // transmission zero of order two, Zs(s) = s², which has no static limit.
  CubeSurfaceTerms common({RationalImpedance({5}, {1.0, 0.0, 0.0}, {0.5, 2.0, 0.0})});
  auto fb_common = common.Coefficient();
  common.Get().AddStaticExtraSystemBdrMassCoefficients(fb_common);
  CHECK_THAT(common.Value(fb_common, 5), WithinRel(cube.Value(fb, 5), 1.0e-14));
  CubeSurfaceTerms singular({RationalImpedance({5}, {1.0, 0.0, 0.0}, {1.0})});
  auto fb_singular = singular.Coefficient();
  CHECK_THROWS(singular.Get().AddStaticExtraSystemBdrMassCoefficients(fb_singular));
}

}  // namespace palace
