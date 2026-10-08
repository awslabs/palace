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

// Profile for a PML on the +x face of the physical domain [0, 1]³: the layer is
// x ∈ [1, 1.1] with grading order 3.
pml::Profile MakeXProfile(double sigma_max = 2.0)
{
  pml::Profile p;
  p.geometry.inner = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
  p.geometry.thickness = {0.0, 0.1, 0.0, 0.0, 0.0, 0.0};
  p.sigma_max = {0.0, sigma_max, 0.0, 0.0, 0.0, 0.0};
  p.order = 3;
  p.reference_frequency = 1.5;
  return p;
}

// Box PML of the physical domain [-1, 1]³ with different thicknesses on all faces and a
// general anisotropic, lossy background.
pml::Profile MakeBoxProfile()
{
  pml::Profile p;
  p.geometry.inner = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  p.geometry.thickness = {0.2, 0.3, 0.25, 0.2, 0.4, 0.1};
  p.sigma_max = {3.0, 2.0, 4.0, 1.0, 2.5, 5.0};
  p.kappa_max = {2.0, 1.5, 3.0};
  p.alpha_max = {0.1, 0.2, 0.05};
  p.order = 2;
  p.reference_frequency = 2.0;
  // Symmetric positive definite background tensors (column-major).
  p.mu_inv = {0.8, 0.1, 0.0, 0.1, 0.9, 0.05, 0.0, 0.05, 1.1};
  p.epsilon_real = {4.0, 0.3, 0.2, 0.3, 5.0, 0.1, 0.2, 0.1, 6.0};
  p.epsilon_imag = {-4e-3, 0.0, 0.0, 0.0, -5e-3, 0.0, 0.0, 0.0, -6e-3};
  return p;
}

// Independent reference evaluation of the stretch factors with std::complex arithmetic.
std::array<std::complex<double>, 3> ReferenceStretch(const pml::Profile &p,
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
std::array<std::complex<double>, 9> ReferenceMuInv(const pml::Profile &p,
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
      T[i + 3 * j] = p.mu_inv[i + 3 * j] * s[i] * s[j] / det;
    }
  }
  return T;
}

std::array<std::complex<double>, 9> ReferenceEps(const pml::Profile &p,
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
      T[i + 3 * j] = (p.epsilon_real[i + 3 * j] + 1i * p.epsilon_imag[i + 3 * j]) * det /
                     (s[i] * s[j]);
    }
  }
  return T;
}

// Evaluate the device helpers for libCEED attribute 1, mapped to the given profile.
bool EvalCoeff(const pml::Profile &p, const pml::ContextHeader &header, bool muinv,
               const std::array<double, 3> &x, std::array<double, 9> &coeff,
               std::vector<int> attr_to_profile = {0})
{
  const auto ctx = pml::PackContext(header, attr_to_profile, {p}, p.frequency_dependent);
  REQUIRE(!ctx.empty());
  const CeedScalar xp[3] = {x[0], x[1], x[2]};
  return muinv ? PMLMuInvCoeff(ctx.data(), 1, xp, coeff.data())
               : PMLEpsCoeff(ctx.data(), 1, xp, coeff.data());
}

void CheckCoeff(const pml::Profile &p, std::complex<double> c, std::complex<double> omega,
                const std::array<double, 3> &x)
{
  for (bool muinv : {true, false})
  {
    const auto T = muinv ? ReferenceMuInv(p, x, omega) : ReferenceEps(p, x, omega);
    for (auto part : {pml::TensorPart::REAL, pml::TensorPart::IMAG, pml::TensorPart::ABS})
    {
      pml::ContextHeader header;
      header.part = part;
      (muinv ? header.c_muinv : header.c_eps) = c;
      header.omega = omega;
      std::array<double, 9> coeff;
      REQUIRE(EvalCoeff(p, header, muinv, x, coeff));
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
  data.direction_signs = {0, 1, -1, 0, 0, 0};
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

TEST_CASE("PML::ResolveSigmaMax", "[pml][Serial]")
{
  config::PMLData data;
  data.order = 3;
  data.reflection_target = 1.0e-6;
  pml::LayerGeometry g;
  g.thickness = {0.0, 0.1, 0.2, 0.4, 0.0, 0.0};

  SECTION("Default from the reflection target, per face thickness")
  {
    const std::array<double, 6> n_r = {1.0, 2.0, 3.0, 1.5, 1.0, 1.0};
    const auto sigma_max = pml::ResolveSigmaMax(data, g, n_r);
    for (int f = 0; f < 6; f++)
    {
      const double expected =
          (g.thickness[f] > 0.0) ? -4.0 * std::log(1.0e-6) / (2.0 * g.thickness[f] * n_r[f])
                                 : 0.0;
      CHECK(sigma_max[f] == Approx(expected));
    }
  }

  SECTION("Configured values per axis")
  {
    data.sigma_max = {-1.0, 7.0, -1.0};
    std::array<double, 6> n_r;
    n_r.fill(1.0);
    const auto sigma_max = pml::ResolveSigmaMax(data, g, n_r);
    CHECK(sigma_max[2] == Approx(7.0));
    CHECK(sigma_max[3] == Approx(7.0));
    CHECK(sigma_max[1] == Approx(-4.0 * std::log(1.0e-6) / (2.0 * 0.1)));
  }
}

TEST_CASE("PML::BuildProfile", "[pml][Serial]")
{
  config::PMLData data;
  data.order = 2;
  data.reflection_target = 1.0e-4;
  data.kappa_max = {1.0, 2.0, 3.0};
  data.alpha_max = {0.0, 0.1, 0.2};
  data.reference_frequency = 5.0;
  pml::LayerGeometry g;
  g.thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25};
  const std::array<double, 9> mu_inv = {0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5};
  const std::array<double, 9> eps = {8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0};
  // n_r = sqrt(8 · 2) = 4.
  const double n = pml::RefractiveIndex(mu_inv, eps);
  CHECK(n == Approx(4.0));
  std::array<double, 6> n_r;
  n_r.fill(n);
  auto p = pml::BuildProfile(data, g, n_r, mu_inv, eps, {});
  CHECK(p.reference_frequency == Approx(5.0));
  CHECK(p.kappa_max == data.kappa_max);
  CHECK(p.sigma_max[5] == Approx(-3.0 * std::log(1.0e-4) / (2.0 * 0.25 * 4.0)));

  data.frequency_dependent = true;
  p = pml::BuildProfile(data, g, n_r, mu_inv, eps, {});
  CHECK(p.frequency_dependent);
  CHECK(p.reference_frequency == 0.0);
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

TEST_CASE("PML::CheckStretchConsistency", "[pml][Serial]")
{
  auto p = MakeBoxProfile(), q = MakeBoxProfile();
  q.epsilon_real[0] = 11.0;  // Background material properties may differ
  CHECK(pml::CheckStretchConsistency(p, q, 1.0e-9).empty());

  // A face only active in one of the profiles is not shared.
  q.geometry.thickness[0] = 0.0;
  q.sigma_max[0] = 7.0;
  CHECK(pml::CheckStretchConsistency(p, q, 1.0e-9).empty());

  auto Check = [&](auto &&modify, const std::string &param)
  {
    auto r = MakeBoxProfile();
    modify(r);
    const auto msg = pml::CheckStretchConsistency(p, r, 1.0e-9);
    CHECK(msg.find(param) != std::string::npos);
  };
  Check([](auto &r) { r.sigma_max[3] *= 1.01; }, "SigmaMax");
  Check([](auto &r) { r.kappa_max[2] = 1.0; }, "KappaMax");
  Check([](auto &r) { r.alpha_max[0] = 0.0; }, "AlphaMax");
  Check([](auto &r) { r.order = 3; }, "Order");
  Check([](auto &r) { r.geometry.thickness[5] = 0.2; }, "Thickness");
  Check([](auto &r) { r.reference_frequency = 3.0; }, "ReferenceFrequency");
  Check([](auto &r) { r.frequency_dependent = true; }, "FrequencyDependent");
}

TEST_CASE("PML::ComputeDepthFraction", "[pml][Serial]")
{
  const auto p = MakeBoxProfile();
  const auto r = pml::ComputeDepthFraction(p, {{-1.1, 0.5, 1.2}});
  CHECK(r[0] == Approx(0.5));
  CHECK(r[1] == Approx(0.0));
  CHECK(r[2] == Approx(1.0));  // Clamped
}

TEST_CASE("PML tensors for a single-axis layer", "[pml][Serial]")
{
  // Classic UPML for a +x layer in vacuum: s = 1 - i σ / ω₀, μ̃⁻¹ = diag(s, 1/s, 1/s),
  // ε̃ = diag(1/s, s, s).
  const auto p = MakeXProfile(2.0);
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
      REQUIRE(EvalCoeff(p, header, m, x, coeff));
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
    REQUIRE(EvalCoeff(p, header, true, xp, coeff));
    CHECK(coeff[0] == Approx(1.0));
    CHECK(coeff[4] == Approx(1.0));
    CHECK(coeff[8] == Approx(1.0));
  }
}

TEST_CASE("PML tensors match the reference evaluation", "[pml][Serial]")
{
  auto p = MakeBoxProfile();
  // Points in the physical region, face, edge, and corner regions of the layer.
  const std::array<std::array<double, 3>, 5> points = {{{0.2, -0.3, 0.4},
                                                        {1.2, 0.0, 0.0},
                                                        {-1.15, 1.1, 0.5},
                                                        {1.25, -1.2, -1.35},
                                                        {-1.3, 1.3, 1.2}}};
  SECTION("Static profile")
  {
    for (const auto &x : points)
    {
      CheckCoeff(p, 1.0, 0.0, x);
      CheckCoeff(p, {-2.5, 0.7}, 0.0, x);
    }
  }

  SECTION("Frequency-dependent profile at real and complex frequencies")
  {
    p.frequency_dependent = true;
    p.reference_frequency = 0.0;
    for (const auto &x : points)
    {
      for (std::complex<double> omega :
           {std::complex<double>{1.7, 0.0}, {1.7, 0.3}, {0.9, -0.2}})
      {
        CheckCoeff(p, 1.0, omega, x);
        CheckCoeff(p, -omega * omega, omega, x);
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
      EvalCoeff(p, header, false, x, re);
      header.part = pml::TensorPart::IMAG;
      EvalCoeff(p, header, false, x, im);
      return std::complex<double>(re[4], im[4]);
    };
    const auto dFdr = (F(omega + h) - F(omega - h)) / (2.0 * h);
    const auto dFdi = (F(omega + 1i * h) - F(omega - 1i * h)) / (2.0 * h);
    CHECK(std::abs(dFdi - 1i * dFdr) < 1.0e-6 * std::abs(dFdr));
  }
}

TEST_CASE("PML::PackContext", "[pml][Serial]")
{
  auto p0 = MakeXProfile(1.0), p1 = MakeBoxProfile(), p2 = MakeXProfile(3.0);
  p1.frequency_dependent = true;
  const std::vector<pml::Profile> profiles = {p0, p1, p2};
  const std::vector<int> attr_to_profile = {2, -1, 0, 1};

  pml::ContextHeader header;
  header.part = pml::TensorPart::IMAG;
  header.c_muinv = {1.0, 2.0};
  header.c_eps = {3.0, 4.0};
  header.omega = {5.0, 6.0};
  header.wave_vector_cross = {0.0, 1.0, -2.0, -1.0, 0.0, 3.0, 2.0, -3.0, 0.0};

  SECTION("Static profiles")
  {
    const auto ctx = pml::PackContext(header, attr_to_profile, profiles, false);
    CHECK(PMLNumAttr(ctx.data()) == 4);
    CHECK(PMLPart(ctx.data()) == PALACE_PML_PART_IM);
    CHECK(ctx[2].second == 1.0);
    CHECK(ctx[5].second == 4.0);
    CHECK(ctx[7].second == 6.0);
    CHECK(PMLWaveVectorCross(ctx.data())[5].second == 3.0);
    CHECK(ctx.size() ==
          PALACE_PML_HEADER_SIZE + attr_to_profile.size() + 2 * PALACE_PML_PROFILE_SIZE);
    CHECK(PMLProfileData(ctx.data(), 2) == nullptr);
    CHECK(PMLProfileData(ctx.data(), 4) == nullptr);  // Frequency-dependent profile
    CHECK(PMLProfileData(ctx.data(), 5) == nullptr);  // Out of range
    // Attribute 1 → profile 2 (σ_max = 3), attribute 3 → profile 0 (σ_max = 1).
    CHECK(PMLProfileData(ctx.data(), 1)[16].second == 3.0);
    CHECK(PMLProfileData(ctx.data(), 3)[16].second == 1.0);
  }

  SECTION("Frequency-dependent profiles")
  {
    const auto ctx = pml::PackContext(header, attr_to_profile, profiles, true);
    CHECK(ctx.size() ==
          PALACE_PML_HEADER_SIZE + attr_to_profile.size() + PALACE_PML_PROFILE_SIZE);
    CHECK(PMLProfileData(ctx.data(), 1) == nullptr);
    CHECK(PMLProfileData(ctx.data(), 4)[0].first == 1);
    CHECK(PMLProfileData(ctx.data(), 4)[36 + 4].second == p1.epsilon_real[4]);
  }

  SECTION("No matching attribute")
  {
    CHECK(pml::PackContext(header, {-1, 0, 2}, profiles, true).empty());
    CHECK(pml::PackContext(header, {-1, -1}, profiles, false).empty());
  }
}

TEST_CASE("MaterialOperator PML stretch consistency", "[pml][materialoperator][Serial]")
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
    std::vector<config::MaterialData> materials(pml_x ? 6 : 4);
    for (int k = 0; k < static_cast<int>(materials.size()); k++)
    {
      auto &m = materials[k];
      m.attributes = {k + 1};
      m.epsilon_r.s.fill((k % 2 == 0) ? 9.0 : 1.0);
      if (k >= 2)
      {
        m.pml = config::PMLData();
        m.pml->reference_frequency = 2.0;
      }
    }
    return materials;
  };
  config::PeriodicBoundaryData periodic;

  SECTION("Default σ_max for the smallest refractive index of the layer")
  {
    auto mesh = MakeMesh(false);
    MaterialOperator mat_op(MakeMaterials(false), periodic, ProblemType::DRIVEN, *mesh);
    const auto &profiles = mat_op.GetPMLProfiles();
    REQUIRE(profiles.size() == 2);
    const double expected = -4.0 * std::log(1.0e-6) / (2.0 * 0.25 * 1.0);
    CHECK(profiles[0].sigma_max[5] == Approx(expected));
    CHECK(profiles[1].sigma_max[5] == Approx(expected));
    CHECK(profiles[0].epsilon_real[0] == Approx(9.0));
    CHECK(profiles[1].epsilon_real[0] == Approx(1.0));
  }

  SECTION("Inconsistent stretch parameters")
  {
    auto mesh = MakeMesh(false);
    auto materials = MakeMaterials(false);
    materials[3].pml->order = 2;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));
    materials = MakeMaterials(false);
    materials[2].pml->sigma_max = {0.0, 0.0, 10.0};
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));
    materials = MakeMaterials(false);
    materials[2].pml->frequency_dependent = true;
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));
  }

  SECTION("Configured layer geometry")
  {
    auto mesh = MakeMesh(false);
    auto materials = MakeMaterials(false);
    for (int k : {2, 3})
    {
      materials[k].pml->autodetect_geometry = false;
      materials[k].pml->direction_signs = {0, 0, 0, 0, 0, 1};
      materials[k].pml->thickness = {0.0, 0.0, 0.0, 0.0, 0.0, 0.25};
    }
    CHECK_NOTHROW(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));

    // Different thicknesses of the same face.
    auto bad = materials;
    bad[3].pml->thickness[5] = 0.5;
    CHECK_THROWS(MaterialOperator(bad, periodic, ProblemType::DRIVEN, *mesh));

    // A thinner layer leaves PML attributes in the physical region.
    bad = materials;
    for (int k : {2, 3})
    {
      bad[k].pml->thickness[5] = 0.125;
    }
    CHECK_THROWS(MaterialOperator(bad, periodic, ProblemType::DRIVEN, *mesh));

    // Vacuum not marked as PML extends into the layer.
    bad = materials;
    bad[3].pml.reset();
    CHECK_THROWS(MaterialOperator(bad, periodic, ProblemType::DRIVEN, *mesh));
  }

  SECTION("PML attributes extending into faces of other materials")
  {
    // The +z attributes cover the corner region x > 0.75, z > 0.75, so their materials
    // need both the +x and +z directions.
    auto mesh = MakeMesh(true);
    auto materials = MakeMaterials(true);
    for (int k = 2; k < 6; k++)
    {
      materials[k].pml->autodetect_geometry = false;
      materials[k].pml->direction_signs = {0, (k >= 4) ? 1 : 0, 0, 0, 0, (k < 4) ? 1 : 0};
      materials[k].pml->thickness = {0.0, 0.25, 0.0, 0.0, 0.0, 0.25};
    }
    CHECK_THROWS(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));
    for (int k : {2, 3})
    {
      materials[k].pml->direction_signs[1] = 1;
    }
    CHECK_NOTHROW(MaterialOperator(materials, periodic, ProblemType::DRIVEN, *mesh));

    // With automatic detection, all PML materials use both faces.
    materials = MakeMaterials(true);
    MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, *mesh);
    for (const auto &p : mat_op.GetPMLProfiles())
    {
      CHECK(p.geometry.thickness[1] == Approx(0.25));
      CHECK(p.geometry.thickness[5] == Approx(0.25));
    }
  }
}

}  // namespace palace
