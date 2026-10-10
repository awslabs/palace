// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "fem/qfunctions/coeff/pml_qf.h"
#include "models/materialoperator.hpp"
#include "models/pml.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

namespace palace
{
using namespace Catch;
using namespace std::complex_literals;

namespace
{

// Stretch for a PML on the +x face of the physical domain [0, 1]³: the layer is x ∈
// [1, 1.1] with grading order 3.
pml::Stretch MakeXStretch(double sigma_max = 2.0)
{
  pml::Stretch p;
  p.geometry.inner = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  p.geometry.thickness = {0.0, 0.1, 0.0, 0.0, 0.0, 0.0};
  p.sigma_max = {0.0, sigma_max, 0.0, 0.0, 0.0, 0.0};
  p.order = 3;
  p.reference_frequency = 1.5;
  return p;
}

// Box PML of the physical domain [-1, 1]³ with different thicknesses on all faces.
pml::Stretch MakeBoxStretch()
{
  pml::Stretch p;
  p.geometry.inner = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  p.geometry.thickness = {0.2, 0.3, 0.25, 0.2, 0.4, 0.1};
  p.sigma_max = {3.0, 2.0, 4.0, 1.0, 2.5, 5.0};
  p.kappa_max = {2.0, 1.5, 3.0};
  p.alpha_max = {0.1, 0.2, 0.05};
  p.order = 2;
  p.reference_frequency = 2.0;
  return p;
}

// General anisotropic, lossy background (symmetric positive definite tensors,
// column-major).
pml::Background MakeAnisotropicBackground()
{
  pml::Background b;
  b.mu_inv = {0.8, 0.1, 0.0, 0.1, 0.9, 0.05, 0.0, 0.05, 1.1};
  b.epsilon_real = {4.0, 0.3, 0.2, 0.3, 5.0, 0.1, 0.2, 0.1, 6.0};
  b.epsilon_imag = {-4e-3, 0.0, 0.0, 0.0, -5e-3, 0.0, 0.0, 0.0, -6e-3};
  return b;
}

// Independent reference evaluation of the stretch factors with std::complex arithmetic.
std::array<std::complex<double>, 3> ReferenceStretch(const pml::Stretch &p,
                                                     const std::array<double, 3> &x,
                                                     std::complex<double> omega)
{
  if (!p.frequency_dependent)
  {
    omega = p.reference_frequency;
  }
  std::array<std::complex<double>, 3> s;
  for (int a = 0; a < 3; a++)
  {
    s[a] = 1.0;
    const auto &g = p.geometry;
    double r = 0.0, sigma_max = 0.0;
    if (g.thickness[2 * a] > 0.0 && x[a] < g.inner[2 * a])
    {
      r = (g.inner[2 * a] - x[a]) / g.thickness[2 * a];
      sigma_max = p.sigma_max[2 * a];
    }
    else if (g.thickness[2 * a + 1] > 0.0 && x[a] > g.inner[2 * a + 1])
    {
      r = (x[a] - g.inner[2 * a + 1]) / g.thickness[2 * a + 1];
      sigma_max = p.sigma_max[2 * a + 1];
    }
    else
    {
      continue;
    }
    const double shape = std::pow(std::min(r, 1.0), p.order);
    s[a] = 1.0 + (p.kappa_max[a] - 1.0) * shape +
           sigma_max * shape / (p.alpha_max[a] * shape + 1i * omega);
  }
  return s;
}

// Reference tensors μ̃⁻¹ = S μ⁻¹ S / det(S) and ε̃ = det(S) S⁻¹ ε S⁻¹ (column-major).
std::array<std::complex<double>, 9> ReferenceMuInv(const pml::Stretch &p,
                                                   const pml::Background &b,
                                                   const std::array<double, 3> &x,
                                                   std::complex<double> omega)
{
  const auto s = ReferenceStretch(p, x, omega);
  const auto det = s[0] * s[1] * s[2];
  std::array<std::complex<double>, 9> T;
  for (int j = 0; j < 3; j++)
  {
    for (int i = 0; i < 3; i++)
    {
      T[i + 3 * j] = b.mu_inv[i + 3 * j] * s[i] * s[j] / det;
    }
  }
  return T;
}

std::array<std::complex<double>, 9> ReferenceEps(const pml::Stretch &p,
                                                 const pml::Background &b,
                                                 const std::array<double, 3> &x,
                                                 std::complex<double> omega)
{
  const auto s = ReferenceStretch(p, x, omega);
  const auto det = s[0] * s[1] * s[2];
  std::array<std::complex<double>, 9> T;
  for (int j = 0; j < 3; j++)
  {
    for (int i = 0; i < 3; i++)
    {
      T[i + 3 * j] = (b.epsilon_real[i + 3 * j] + 1i * b.epsilon_imag[i + 3 * j]) * det /
                     (s[i] * s[j]);
    }
  }
  return T;
}

// Evaluate the device helpers for libCEED attribute 1 in a PML region with the given
// background material.
bool EvalCoeff(const pml::Stretch &p, const pml::Background &b,
               const pml::ContextHeader &header, bool muinv, const std::array<double, 3> &x,
               std::array<double, 9> &coeff)
{
  const auto ctx = pml::PackContext(header, p, {0}, {b});
  REQUIRE(!ctx.empty());
  const CeedScalar xp[3] = {x[0], x[1], x[2]};
  return muinv ? PMLMuInvCoeff(ctx.data(), 1, xp, coeff.data())
               : PMLEpsCoeff(ctx.data(), 1, xp, coeff.data());
}

void CheckCoeff(const pml::Stretch &p, const pml::Background &b, std::complex<double> c,
                std::complex<double> omega, const std::array<double, 3> &x)
{
  for (bool muinv : {true, false})
  {
    const auto T = muinv ? ReferenceMuInv(p, b, x, omega) : ReferenceEps(p, b, x, omega);
    for (auto part : {pml::TensorPart::REAL, pml::TensorPart::IMAG, pml::TensorPart::ABS})
    {
      pml::ContextHeader header;
      header.part = part;
      (muinv ? header.c_muinv : header.c_eps) = c;
      header.omega = omega;
      std::array<double, 9> coeff;
      REQUIRE(EvalCoeff(p, b, header, muinv, x, coeff));
      for (int k = 0; k < 9; k++)
      {
        const double expected = (part == pml::TensorPart::REAL) ? (c * T[k]).real()
                                : (part == pml::TensorPart::IMAG)
                                    ? (c * T[k]).imag()
                                    : c.real() * std::abs(T[k]);
        CHECK(coeff[k] == Approx(expected).margin(1.0e-14));
      }
    }
  }
}

}  // namespace

TEST_CASE("PML::DetectLayerGeometry", "[pml][Serial]")
{
  SECTION("Single +z layer")
  {
    const auto g = pml::DetectLayerGeometry({{0.0, 0.0, 0.0}}, {{1.0, 2.0, 1.0}},
                                            {{0.0, 0.0, 0.0}}, {{1.0, 2.0, 1.3}});
    for (int f = 0; f < 5; f++)
    {
      CHECK(g.thickness[f] == 0.0);
    }
    CHECK(g.thickness[5] == Approx(0.3));
    CHECK(g.inner[5] == Approx(1.0));
  }

  SECTION("Box layer with asymmetric thicknesses")
  {
    const auto g = pml::DetectLayerGeometry({{-1.0, -1.0, -1.0}}, {{1.0, 1.0, 1.0}},
                                            {{-1.2, -1.3, -1.4}}, {{1.5, 1.6, 1.7}});
    const std::array<double, 6> t = {0.2, 0.5, 0.3, 0.6, 0.4, 0.7};
    const std::array<double, 6> inner = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    for (int f = 0; f < 6; f++)
    {
      CHECK(g.thickness[f] == Approx(t[f]));
      CHECK(g.inner[f] == Approx(inner[f]));
    }
  }

  SECTION("Differences below the tolerance are not PML faces")
  {
    const auto g = pml::DetectLayerGeometry({{0.0, 0.0, 1.0e-9}}, {{1.0, 1.0, 1.0}},
                                            {{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}});
    CHECK(g.thickness == std::array<double, 6>{});
  }
}

TEST_CASE("PML::ConfiguredLayerGeometry", "[pml][Serial]")
{
  config::PMLData data;
  data.directions = {false, true, true, false, false, false};
  data.thickness = {0.5, 0.1, 0.2, 0.0, 0.0, 0.0};  // −x thickness without direction
  const auto g =
      pml::ConfiguredLayerGeometry(data, {{-1.0, -1.0, -1.0}}, {{1.0, 1.0, 1.0}});
  const std::array<double, 6> t = {0.0, 0.1, 0.2, 0.0, 0.0, 0.0};
  for (int f = 0; f < 6; f++)
  {
    CHECK(g.thickness[f] == Approx(t[f]));
  }
  CHECK(g.inner[1] == Approx(0.9));
  CHECK(g.inner[2] == Approx(-0.8));
}

TEST_CASE("PML::BuildStretch", "[pml][Serial]")
{
  config::PMLData data;
  data.order = 3;
  data.reflection_target = 1.0e-6;
  pml::LayerGeometry g;
  g.thickness = {0.0, 0.1, 0.2, 0.4, 0.0, 0.0};

  SECTION("Default σ_max from the reflection target, per face thickness")
  {
    const double n_r = 1.5;
    const auto p = pml::BuildStretch(data, g, n_r);
    for (int f = 0; f < 6; f++)
    {
      const double expected = (g.thickness[f] > 0.0)
                                  ? -4.0 * std::log(1.0e-6) / (2.0 * g.thickness[f] * n_r)
                                  : 0.0;
      CHECK(p.sigma_max[f] == Approx(expected));
    }
  }

  SECTION("Configured σ_max per axis")
  {
    data.sigma_max = {-1.0, 7.0, -1.0};
    const auto p = pml::BuildStretch(data, g, 1.0);
    CHECK(p.sigma_max[2] == Approx(7.0));
    CHECK(p.sigma_max[3] == Approx(7.0));
    CHECK(p.sigma_max[1] == Approx(-4.0 * std::log(1.0e-6) / (2.0 * 0.1)));
  }

  SECTION("Static and frequency-dependent stretch")
  {
    data.order = 2;
    data.kappa_max = {1.0, 2.0, 3.0};
    data.reference_frequency = 5.0;
    auto p = pml::BuildStretch(data, g, 4.0);
    CHECK(p.order == 2);
    CHECK(p.kappa_max == data.kappa_max);
    CHECK(!p.frequency_dependent);
    CHECK(p.reference_frequency == Approx(5.0));

    data.frequency_dependent = true;
    data.alpha_max = {0.0, 0.1, 0.2};
    p = pml::BuildStretch(data, g, 4.0);
    CHECK(p.frequency_dependent);
    CHECK(p.alpha_max == data.alpha_max);
    CHECK(p.reference_frequency == 0.0);
  }
}

TEST_CASE("PML::RefractiveIndex", "[pml][Serial]")
{
  // Smallest principal value of a rotated anisotropic permittivity (eigenvalues 2, 5, 9)
  // and permeability (eigenvalues 1, 1.5, 4).
  const double c = std::cos(0.4), s = std::sin(0.4);
  auto Rotated = [c, s](double a, double b, double d)
  {
    // R diag(a, b, d) Rᵀ for a rotation about z (column-major, symmetric).
    return std::array<double, 9>{c * c * a + s * s * b,
                                 c * s * (a - b),
                                 0.0,
                                 c * s * (a - b),
                                 s * s * a + c * c * b,
                                 0.0,
                                 0.0,
                                 0.0,
                                 d};
  };
  const auto eps = Rotated(5.0, 2.0, 9.0);
  const auto mu_inv = Rotated(1.0, 1.0 / 1.5, 0.25);
  CHECK(pml::RefractiveIndex(mu_inv, eps) == Approx(std::sqrt(2.0 * 1.0)));
}

TEST_CASE("PML::EvaluateStretch", "[pml][Serial]")
{
  auto p = MakeBoxStretch();
  const std::complex<double> omega = 1.3 - 0.2i;
  for (bool frequency_dependent : {false, true})
  {
    p.frequency_dependent = frequency_dependent;
    for (const std::array<double, 3> &x :
         {std::array<double, 3>{0.0, 0.0, 0.0}, std::array<double, 3>{-1.1, 0.5, 1.05},
          std::array<double, 3>{1.2, -1.2, -1.3}, std::array<double, 3>{1.4, 1.3, -1.5}})
    {
      const auto s = pml::EvaluateStretch(p, x, omega);
      const auto s_ref = ReferenceStretch(p, x, omega);
      for (int a = 0; a < 3; a++)
      {
        CHECK(s[a].real() == Approx(s_ref[a].real()).margin(1.0e-14));
        CHECK(s[a].imag() == Approx(s_ref[a].imag()).margin(1.0e-14));
      }
    }
  }
}

TEST_CASE("PML tensors for a single-axis layer", "[pml][Serial]")
{
  // Classic UPML for a +x layer in vacuum: s = 1 - i σ / ω₀, μ̃⁻¹ = diag(s, 1/s, 1/s),
  // ε̃ = diag(1/s, s, s).
  const auto p = MakeXStretch(2.0);
  const pml::Background b;
  const std::array<double, 3> x = {1.05, 0.5, 0.5};
  const double sigma = 2.0 * std::pow(0.5, 3);
  const std::complex<double> s = 1.0 - 1i * sigma / p.reference_frequency;
  const std::array<std::complex<double>, 3> muinv = {s, 1.0 / s, 1.0 / s},
                                            eps = {1.0 / s, s, s};
  for (bool m : {true, false})
  {
    for (auto part : {pml::TensorPart::REAL, pml::TensorPart::IMAG})
    {
      pml::ContextHeader header;
      header.part = part;
      header.c_muinv = header.c_eps = 1.0;
      std::array<double, 9> coeff;
      REQUIRE(EvalCoeff(p, b, header, m, x, coeff));
      for (int i = 0; i < 3; i++)
      {
        const auto T = m ? muinv[i] : eps[i];
        CHECK(coeff[4 * i] ==
              Approx((part == pml::TensorPart::REAL) ? T.real() : T.imag()));
      }
      CHECK(coeff[1] == 0.0);
      CHECK(coeff[3] == 0.0);
    }
  }

  // No stretch in the physical region and at the interface.
  for (const auto &xp : {std::array<double, 3>{0.5, 0.5, 0.5}, {1.0, 0.5, 0.5}})
  {
    pml::ContextHeader header;
    header.c_muinv = 1.0;
    std::array<double, 9> coeff;
    REQUIRE(EvalCoeff(p, b, header, true, xp, coeff));
    CHECK(coeff[0] == Approx(1.0));
    CHECK(coeff[4] == Approx(1.0));
    CHECK(coeff[8] == Approx(1.0));
  }
}

TEST_CASE("PML tensors match the reference evaluation", "[pml][Serial]")
{
  auto p = MakeBoxStretch();
  const auto b = MakeAnisotropicBackground();
  // Points in the physical region, face, edge, and corner regions of the layer.
  const std::array<std::array<double, 3>, 5> points = {{{0.2, -0.3, 0.4},
                                                        {1.2, 0.0, 0.0},
                                                        {-1.15, 1.1, 0.5},
                                                        {1.25, -1.2, -1.35},
                                                        {-1.3, 1.3, 1.2}}};
  SECTION("Static stretch")
  {
    for (const auto &x : points)
    {
      CheckCoeff(p, b, 1.0, 0.0, x);
      CheckCoeff(p, b, {-2.5, 0.7}, 0.0, x);
    }
  }

  SECTION("Frequency-dependent stretch at real and complex frequencies")
  {
    p.frequency_dependent = true;
    p.reference_frequency = 0.0;
    for (const auto &x : points)
    {
      for (std::complex<double> omega :
           {std::complex<double>{1.7, 0.0}, {1.7, 0.3}, {0.9, -0.2}})
      {
        CheckCoeff(p, b, 1.0, omega, x);
        CheckCoeff(p, b, -omega * omega, omega, x);
      }
    }
  }

  SECTION("Frequency-dependent stretch is the analytic continuation in ω")
  {
    // The tensors are holomorphic in ω: check the Cauchy-Riemann equations with central
    // differences for one entry.
    p.frequency_dependent = true;
    const std::array<double, 3> x = {1.25, -1.2, -1.35};
    const std::complex<double> omega = {1.3, 0.2};
    const double h = 1.0e-6;
    auto F = [&](std::complex<double> w)
    {
      pml::ContextHeader header;
      header.c_eps = 1.0;
      header.omega = w;
      std::array<double, 9> re, im;
      header.part = pml::TensorPart::REAL;
      EvalCoeff(p, b, header, false, x, re);
      header.part = pml::TensorPart::IMAG;
      EvalCoeff(p, b, header, false, x, im);
      return std::complex<double>(re[4], im[4]);
    };
    const auto dFdr = (F(omega + h) - F(omega - h)) / (2.0 * h);
    const auto dFdi = (F(omega + 1i * h) - F(omega - 1i * h)) / (2.0 * h);
    CHECK(std::abs(dFdi - 1i * dFdr) < 1.0e-6 * std::abs(dFdr));
  }
}

TEST_CASE("PML::PackContext", "[pml][Serial]")
{
  const auto p = MakeBoxStretch();
  pml::Background b0, b1 = MakeAnisotropicBackground();
  const std::vector<int> attr_background = {1, -1, 0, 1};

  pml::ContextHeader header;
  header.part = pml::TensorPart::IMAG;
  header.c_muinv = {1.0, 2.0};
  header.c_eps = {3.0, 4.0};
  header.omega = {5.0, 6.0};
  header.wave_vector_cross = {0.0, 1.0, -2.0, -1.0, 0.0, 3.0, 2.0, -3.0, 0.0};

  SECTION("Static stretch")
  {
    const auto ctx = pml::PackContext(header, p, attr_background, {b0, b1});
    CHECK(ctx.size() == PALACE_PML_HEADER_SIZE + PALACE_PML_STRETCH_SIZE +
                            attr_background.size() + 2 * PALACE_PML_BACKGROUND_SIZE);
    CHECK(PMLNumAttr(ctx.data()) == 4);
    CHECK(PMLPart(ctx.data()) == PALACE_PML_PART_IM);
    CHECK(ctx[2].second == 1.0);
    CHECK(ctx[5].second == 4.0);
    CHECK(ctx[6].second == p.reference_frequency);  // Stretch at ω₀
    CHECK(ctx[7].second == 0.0);
    CHECK(PMLWaveVectorCross(ctx.data())[5].second == 3.0);
    CHECK(ctx[PALACE_PML_HEADER_SIZE].first == p.order);
    CHECK(ctx[PALACE_PML_HEADER_SIZE + 13 + 5].second == p.sigma_max[5]);
    CHECK(PMLBackgroundData(ctx.data(), 2) == nullptr);
    CHECK(PMLBackgroundData(ctx.data(), 5) == nullptr);  // Out of range
    // Attributes 1 and 4 → background 1, attribute 3 → background 0.
    CHECK(PMLBackgroundData(ctx.data(), 1)[9 + 4].second == b1.epsilon_real[4]);
    CHECK(PMLBackgroundData(ctx.data(), 4)[18].second == b1.epsilon_imag[0]);
    CHECK(PMLBackgroundData(ctx.data(), 3)[9].second == 1.0);
  }

  SECTION("Frequency-dependent stretch")
  {
    auto q = p;
    q.frequency_dependent = true;
    const auto ctx = pml::PackContext(header, q, attr_background, {b0, b1});
    CHECK(ctx[6].second == 5.0);  // Stretch at the solve frequency
    CHECK(ctx[7].second == 6.0);
  }

  SECTION("No attribute in a PML region")
  {
    CHECK(pml::PackContext(header, p, {-1, -1}, {b0}).empty());
  }
}

TEST_CASE("MaterialOperator PML regions", "[pml][materialoperator][Serial]")
{
  // Unit cube with a substrate (y < 0.5) and vacuum (y > 0.5), both of which cross a PML
  // layer on the +z face (z > 0.75) and, optionally, the +x face (x > 0.75).
  auto MakeMesh = [](bool pml_x)
  {
    auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
    for (int i = 0; i < smesh.GetNE(); i++)
    {
      mfem::Vector c;
      smesh.GetElementCenter(i, c);
      const bool sub = c(1) < 0.5, pml_z = c(2) > 0.75, in_x = pml_x && c(0) > 0.75;
      // 1, 2: physical substrate, vacuum. 3, 4: +z PML (incl. corner). 5, 6: +x PML.
      smesh.SetAttribute(i, (pml_z ? 3 : (in_x ? 5 : 1)) + (sub ? 0 : 1));
    }
    smesh.SetAttributes();
    return std::make_unique<Mesh>(Mpi::World(), smesh);
  };
  auto MakeMaterials = [](bool pml_x)
  {
    std::vector<config::MaterialData> materials(2);
    materials[0].attributes = {1, 3};
    materials[1].attributes = {2, 4};
    if (pml_x)
    {
      materials[0].attributes.push_back(5);
      materials[1].attributes.push_back(6);
    }
    materials[0].epsilon_r.s.fill(9.0);
    return materials;
  };
  auto MakePML = [](bool pml_x)
  {
    config::PMLData pml;
    pml.attributes = pml_x ? std::vector<int>{3, 4, 5, 6} : std::vector<int>{3, 4};
    pml.reference_frequency = 2.0;
    return pml;
  };
  config::PeriodicBoundaryData periodic;

  SECTION("Default σ_max for the smallest refractive index of the PML materials")
  {
    auto mesh = MakeMesh(false);
    const auto pml = MakePML(false);
    MaterialOperator mat_op(MakeMaterials(false), periodic, ProblemType::DRIVEN, *mesh,
                            {pml});
    REQUIRE(mat_op.HasPML());
    CHECK(!mat_op.HasFrequencyDependentPML());
    const auto &stretch = mat_op.GetPMLLayers().at(0).GetStretch();
    CHECK(stretch.geometry.thickness[5] == Approx(0.25));
    CHECK(stretch.sigma_max[5] == Approx(-4.0 * std::log(1.0e-6) / (2.0 * 0.25 * 1.0)));
    CHECK(stretch.reference_frequency == Approx(2.0));
    CHECK(mat_op.GetPMLLayers().at(0).GetAttributes() == std::vector<int>{3, 4});
  }

  SECTION("Background material properties and bulk attribute map")
  {
    // The materials are shared by physical and PML attributes. Their properties are not
    // modified, and the PML attributes are only excluded from the bulk attribute map.
    auto mesh = MakeMesh(false);
    const auto pml = MakePML(false);
    MaterialOperator mat_op(MakeMaterials(false), periodic, ProblemType::DRIVEN, *mesh,
                            {pml});
    const auto &loc_attr = mesh->GetCeedAttributes();
    for (int attr = 1; attr <= 4; attr++)
    {
      const int ceed_attr = loc_attr.at(attr);
      const bool in_pml = (attr >= 3);
      CHECK(mat_op.IsPMLCeedAttribute(ceed_attr) == in_pml);
      CHECK(mat_op.GetAttributeToMaterial()[ceed_attr - 1] >= 0);
      CHECK((mat_op.GetBulkAttributeToMaterial()[ceed_attr - 1] < 0) == in_pml);
      CHECK(mat_op.GetPermittivityReal(attr)(0, 0) == Approx((attr % 2) ? 9.0 : 1.0));
    }
    pml::ContextHeader header;
    const auto ctx = mat_op.GetPMLLayers().at(0).PackContext(header);
    REQUIRE(!ctx.empty());
    CHECK(PMLBackgroundData(ctx.data(), loc_attr.at(1)) == nullptr);
    CHECK(PMLBackgroundData(ctx.data(), loc_attr.at(3))[9].second == Approx(9.0));
    CHECK(PMLBackgroundData(ctx.data(), loc_attr.at(4))[9].second == Approx(1.0));
  }

  SECTION("Unsupported materials and problem types")
  {
    auto mesh = MakeMesh(false);
    const auto pml = MakePML(false);
    auto materials = MakeMaterials(false);
    materials[1].sigma.s.fill(1.0);
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh, {pml}));

    // Treated as regular materials for other problem types.
    MaterialOperator mat_op(MakeMaterials(false), periodic, ProblemType::ELECTROSTATIC,
                            *mesh, {pml});
    CHECK(!mat_op.HasPML());
    CHECK(mat_op.GetBulkAttributeToMaterial() == mat_op.GetAttributeToMaterial());
  }

  SECTION("Configured layer geometry")
  {
    auto mesh = MakeMesh(false);
    auto pml = MakePML(false);
    pml.autodetect_geometry = false;
    pml.directions = {false, false, false, false, false, true};
    pml.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25};
    CHECK_NOTHROW(MaterialOperator(MakeMaterials(false), periodic, ProblemType::DRIVEN,
                                   *mesh, {pml}));

    // A thinner layer leaves PML attributes in the physical region.
    auto bad = pml;
    bad.thickness[5] = 0.125;
    CHECK_THROWS(MaterialOperator(MakeMaterials(false), periodic, ProblemType::DRIVEN,
                                  *mesh, {bad}));

    // Vacuum not marked as PML extends into the layer.
    bad = pml;
    bad.attributes = {3};
    CHECK_THROWS(MaterialOperator(MakeMaterials(false), periodic, ProblemType::DRIVEN,
                                  *mesh, {bad}));
  }

  SECTION("PML regions beyond an inactive face")
  {
    // The +x PML regions are outside of the layer of a +z PML.
    auto mesh = MakeMesh(true);
    auto pml = MakePML(true);
    pml.autodetect_geometry = false;
    pml.directions = {false, false, false, false, false, true};
    pml.thickness = {0.0, 0.25, 0.0, 0.0, 0.0, 0.25};
    CHECK_THROWS(
        MaterialOperator(MakeMaterials(true), periodic, ProblemType::DRIVEN, *mesh, {pml}));
    pml.directions[1] = true;
    CHECK_NOTHROW(
        MaterialOperator(MakeMaterials(true), periodic, ProblemType::DRIVEN, *mesh, {pml}));

    // With automatic detection, the layer covers both faces.
    pml = MakePML(true);
    MaterialOperator mat_op(MakeMaterials(true), periodic, ProblemType::DRIVEN, *mesh,
                            {pml});
    CHECK(mat_op.GetPMLLayers().at(0).GetStretch().geometry.thickness[1] == Approx(0.25));
    CHECK(mat_op.GetPMLLayers().at(0).GetStretch().geometry.thickness[5] == Approx(0.25));
  }
}

TEST_CASE("MaterialOperator PML blocks", "[pml][materialoperator][Serial][Parallel]")
{
  config::PeriodicBoundaryData periodic;
  auto MakePML = [](std::vector<int> attributes)
  {
    config::PMLData pml;
    pml.attributes = std::move(attributes);
    pml.reference_frequency = 2.0;
    return pml;
  };
  auto Configure = [](config::PMLData &pml, int face, double thickness)
  {
    pml.autodetect_geometry = false;
    pml.directions.fill(false);
    pml.thickness.fill(0.0);
    pml.directions[face] = true;
    pml.thickness[face] = thickness;
  };

  SECTION("Blocks in contact")
  {
    // Unit cube with a substrate (y < 0.5) and vacuum (y > 0.5), with a PML layer on the +z
    // face (z > 0.75, attributes 3 and 4, including the edge with the +x face) and on the
    // +x face (x > 0.75, attributes 5 and 6), in two blocks. In parallel, the elements of
    // the second block are on the last process, so that the interfaces of the blocks are
    // shared faces.
    const int np = Mpi::Size(Mpi::World());
    auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
    std::vector<int> partitioning(smesh.GetNE());
    for (int i = 0; i < smesh.GetNE(); i++)
    {
      mfem::Vector c;
      smesh.GetElementCenter(i, c);
      const bool sub = c(1) < 0.5;
      const int attr = ((c(2) > 0.75) ? 3 : ((c(0) > 0.75) ? 5 : 1)) + (sub ? 0 : 1);
      smesh.SetAttribute(i, attr);
      partitioning[i] = (np == 1) ? 0 : ((attr >= 5) ? np - 1 : i % (np - 1));
    }
    smesh.SetAttributes();
    Mesh mesh(Mpi::World(), smesh, partitioning.data());
    std::vector<config::MaterialData> materials(2);
    materials[0].attributes = {1, 3, 5};
    materials[1].attributes = {2, 4, 6};
    materials[0].epsilon_r.s.fill(9.0);
    const std::vector<config::PMLData> pml = {MakePML({3, 4}), MakePML({5, 6})};

    // With the same parameters, the stretch of the +x face is the same in both blocks.
    MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, mesh, pml);
    REQUIRE(mat_op.GetPMLLayers().size() == 2);
    CHECK(mat_op.GetPMLAttributes() == std::vector<int>{3, 4, 5, 6});
    const auto &g0 = mat_op.GetPMLLayers()[0].GetStretch().geometry;
    const auto &g1 = mat_op.GetPMLLayers()[1].GetStretch().geometry;
    CHECK(g0.thickness == std::array<double, 6>{0.0, 0.25, 0.0, 0.0, 0.0, 0.25});
    CHECK(g1.thickness == std::array<double, 6>{0.0, 0.25, 0.0, 0.0, 0.0, 0.0});
    CHECK(g0.inner[1] == Approx(0.75));
    CHECK(g1.inner[1] == Approx(0.75));
    const auto &loc_attr = mesh.GetCeedAttributes();
    for (int attr = 1; attr <= 6; attr++)
    {
      if (!loc_attr.contains(attr))
      {
        continue;
      }
      const int ceed_attr = loc_attr.at(attr);
      CHECK(mat_op.IsPMLCeedAttribute(ceed_attr) == (attr >= 3));
      CHECK((mat_op.GetBulkAttributeToMaterial()[ceed_attr - 1] < 0) == (attr >= 3));
      CHECK(mat_op.GetPMLLayers()[0].IsPMLCeedAttribute(ceed_attr) ==
            (attr == 3 || attr == 4));
      CHECK(mat_op.GetPMLLayers()[1].IsPMLCeedAttribute(ceed_attr) == (attr >= 5));
    }

    // A different grading, absorption, or frequency dependence of the +x face of the second
    // block makes the stretch discontinuous across the interface z = 0.75.
    auto bad = pml;
    bad[1].order = 2;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
    bad = pml;
    bad[1].reflection_target = 1.0e-4;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
    bad = pml;
    bad[1].kappa_max = {2.0, 1.0, 1.0};
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
    bad = pml;
    bad[1].reference_frequency = 3.0;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
    bad = pml;
    bad[1].frequency_dependent = true;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));

    // The edge regions of the first block must also be stretched along x.
    bad = pml;
    Configure(bad[0], 5, 0.25);
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));

    // The same holds for frequency-dependent stretches.
    auto fd = pml;
    for (auto &data : fd)
    {
      data.frequency_dependent = true;
      data.alpha_max = {0.1, 0.1, 0.1};
    }
    CHECK_NOTHROW(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, fd));
    fd[1].alpha_max = {0.2, 0.1, 0.1};
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, fd));

    // The PML attributes of the blocks are disjoint.
    bad = pml;
    bad[1].attributes = {4, 5, 6};
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
  }

  SECTION("Separate blocks with different stretches")
  {
    // Unit cube with PML layers on the -z face (z < 0.25, attribute 2) and on the +z face
    // (z > 0.75, attribute 3).
    auto smesh = mfem::Mesh::MakeCartesian3D(2, 2, 4, mfem::Element::HEXAHEDRON);
    for (int i = 0; i < smesh.GetNE(); i++)
    {
      mfem::Vector c;
      smesh.GetElementCenter(i, c);
      smesh.SetAttribute(i, (c(2) < 0.25) ? 2 : ((c(2) > 0.75) ? 3 : 1));
    }
    smesh.SetAttributes();
    Mesh mesh(Mpi::World(), smesh);
    std::vector<config::MaterialData> materials(1);
    materials[0].attributes = {1, 2, 3};
    std::vector<config::PMLData> pml = {MakePML({2}), MakePML({3})};
    pml[0].order = 2;
    pml[0].reflection_target = 1.0e-4;
    pml[1].frequency_dependent = true;
    pml[1].kappa_max = {1.0, 1.0, 2.0};
    MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, mesh, pml);
    REQUIRE(mat_op.GetPMLLayers().size() == 2);
    CHECK(mat_op.HasFrequencyDependentPML());
    const auto &p0 = mat_op.GetPMLLayers()[0].GetStretch();
    const auto &p1 = mat_op.GetPMLLayers()[1].GetStretch();
    CHECK(!p0.frequency_dependent);
    CHECK(p1.frequency_dependent);
    CHECK(p0.geometry.thickness == std::array<double, 6>{0.0, 0.0, 0.0, 0.0, 0.25, 0.0});
    CHECK(p1.geometry.thickness == std::array<double, 6>{0.0, 0.0, 0.0, 0.0, 0.0, 0.25});
    CHECK(p0.sigma_max[4] == Approx(-3.0 * std::log(1.0e-4) / (2.0 * 0.25)));
    CHECK(p1.sigma_max[5] == Approx(-4.0 * std::log(1.0e-6) / (2.0 * 0.25)));

    // A single block for both layers.
    MaterialOperator mat_op_single(materials, periodic, ProblemType::DRIVEN, mesh,
                                   {MakePML({2, 3})});
    REQUIRE(mat_op_single.GetPMLLayers().size() == 1);
    CHECK(mat_op_single.GetPMLLayers()[0].GetStretch().geometry.thickness ==
          std::array<double, 6>{0.0, 0.0, 0.0, 0.0, 0.25, 0.25});
  }

  SECTION("Blocks at different depths of the same face")
  {
    // Two columns x < 0.25 and x > 0.75 on a common base z < 0.25 (the elements with 0.25 <
    // x < 0.75 and z > 0.25 are removed): the first column ends at z = 1, with a PML layer
    // z > 0.75 (attribute 3), and the second one at z = 0.75, with a PML layer z > 0.5
    // (attribute 4).
    auto smesh = mfem::Mesh::MakeCartesian3D(4, 1, 4, mfem::Element::HEXAHEDRON, 1.0, 0.25);
    for (int i = 0; i < smesh.GetNE(); i++)
    {
      mfem::Vector c;
      smesh.GetElementCenter(i, c);
      int attr = 1;
      if (c(0) > 0.25 && c(0) < 0.75 && c(2) > 0.25)
      {
        attr = 9;
      }
      else if (c(0) < 0.25 && c(2) > 0.75)
      {
        attr = 3;
      }
      else if (c(0) > 0.75 && c(2) > 0.75)
      {
        attr = 9;
      }
      else if (c(0) > 0.75 && c(2) > 0.5)
      {
        attr = 4;
      }
      smesh.SetAttribute(i, attr);
    }
    smesh.SetAttributes();
    mfem::Array<int> domain_attr({1, 3, 4});
    auto submesh = mfem::SubMesh::CreateFromDomain(smesh, domain_attr);
    Mesh mesh(Mpi::World(), submesh);
    std::vector<config::MaterialData> materials(1);
    materials[0].attributes = {1, 3, 4};
    std::vector<config::PMLData> pml = {MakePML({3}), MakePML({4})};

    // The layer geometry of the second block is not detected from the bounding box of the
    // physical region.
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, pml));

    // The configured thicknesses are measured from the bounding box of each block. The
    // blocks are not in contact and can have different parameters.
    Configure(pml[0], 5, 0.25);
    Configure(pml[1], 5, 0.25);
    pml[1].order = 2;
    MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, mesh, pml);
    REQUIRE(mat_op.GetPMLLayers().size() == 2);
    CHECK(mat_op.GetPMLLayers()[0].GetStretch().geometry.inner[5] == Approx(0.75));
    CHECK(mat_op.GetPMLLayers()[1].GetStretch().geometry.inner[5] == Approx(0.5));

    // A thicker layer extends into the physical region of the second column.
    auto bad = pml;
    bad[1].thickness[5] = 0.5;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, mesh, bad));
  }
}

}  // namespace palace
