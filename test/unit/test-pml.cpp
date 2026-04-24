// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "models/pml.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{
using namespace Catch;

namespace
{

// Construct a minimal Profile for a 1D PML on the +x face (absorption in +x direction
// only). Physical domain is x ∈ [0, 1]; PML is x ∈ [1, 1.1]; d = 0.1, n = 3,
// σ_max = σ, κ_max = 1 (pure UPML), α_max = 0.
pml::Profile MakeSimpleXProfile(double sigma_max = 1.0,
                                PMLStretchFormulation f = PMLStretchFormulation::FIXED,
                                double omega0 = 1.0)
{
  pml::Profile p;
  p.mu_r = 1.0;
  p.epsilon_r = 1.0;
  p.order = 3;
  p.sigma_max = {sigma_max, 0.0, 0.0};
  p.kappa_max = {1.0, 1.0, 1.0};
  p.alpha_max = {0.0, 0.0, 0.0};
  p.thickness = {0.1, 0.0, 0.0};
  p.interface_coord[0][0] = 0.0;  // inner −x interface (unused, −x side inactive)
  p.interface_coord[0][1] = 1.0;  // inner +x interface (x = 1)
  p.interface_coord[1] = {0.0, 0.0};
  p.interface_coord[2] = {0.0, 0.0};
  p.direction_active = {0, 1, 0, 0, 0, 0};  // +x only
  p.formulation = f;
  p.reference_frequency = omega0;
  return p;
}

}  // namespace

TEST_CASE("PML::ComputeDepth axis-aligned", "[pml][Serial]")
{
  const auto p = MakeSimpleXProfile();

  SECTION("Inside physical domain ⇒ depth zero")
  {
    auto d = pml::ComputeDepth(p, {0.5, 0.5, 0.5});
    REQUIRE(d[0] == Approx(0.0));
    REQUIRE(d[1] == Approx(0.0));
    REQUIRE(d[2] == Approx(0.0));
  }

  SECTION("At the inner PML interface ⇒ depth zero")
  {
    auto d = pml::ComputeDepth(p, {1.0, 0.0, 0.0});
    REQUIRE(d[0] == Approx(0.0));
  }

  SECTION("Inside the PML ⇒ depth = x − interface")
  {
    auto d = pml::ComputeDepth(p, {1.04, 0.0, 0.0});
    REQUIRE(d[0] == Approx(0.04));
  }

  SECTION("At the outer PML boundary ⇒ depth = thickness")
  {
    auto d = pml::ComputeDepth(p, {1.1, 0.0, 0.0});
    REQUIRE(d[0] == Approx(0.1));
  }

  SECTION("Negative-x side with only +x active ⇒ no absorption")
  {
    auto d = pml::ComputeDepth(p, {-0.2, 0.0, 0.0});
    REQUIRE(d[0] == Approx(0.0));
  }
}

TEST_CASE("PML::ComputeLocalStretchParams polynomial grading", "[pml][Serial]")
{
  const auto p = MakeSimpleXProfile(/*sigma_max=*/2.0);

  SECTION("At inner interface ⇒ σ = 0, κ = 1, α = 0")
  {
    auto lp = pml::ComputeLocalStretchParams(p, {1.0, 0.0, 0.0});
    REQUIRE(lp.sigma[0] == Approx(0.0));
    REQUIRE(lp.kappa[0] == Approx(1.0));
    REQUIRE(lp.alpha[0] == Approx(0.0));
  }

  SECTION("At outer boundary ⇒ σ = σ_max")
  {
    auto lp = pml::ComputeLocalStretchParams(p, {1.1, 0.0, 0.0});
    REQUIRE(lp.sigma[0] == Approx(2.0));
  }

  SECTION("Midway through: σ = σ_max · (d/D)^n with n=3")
  {
    // Depth = 0.05, D = 0.1 ⇒ t = 0.5, t^3 = 0.125 ⇒ σ = 0.25.
    auto lp = pml::ComputeLocalStretchParams(p, {1.05, 0.0, 0.0});
    REQUIRE(lp.sigma[0] == Approx(2.0 * 0.125));
  }

  SECTION("Inactive axis ⇒ σ = 0, κ = 1")
  {
    auto lp = pml::ComputeLocalStretchParams(p, {1.05, 0.0, 0.0});
    REQUIRE(lp.sigma[1] == Approx(0.0));
    REQUIRE(lp.kappa[1] == Approx(1.0));
    REQUIRE(lp.sigma[2] == Approx(0.0));
    REQUIRE(lp.kappa[2] == Approx(1.0));
  }
}

TEST_CASE("PML::ComputeStretch formulations", "[pml][Serial]")
{
  SECTION("FIXED: s_x = 1 − i σ / ω₀ (Palace e^{+iωt} convention)")
  {
    pml::LocalStretchParams lp;
    lp.sigma = {0.5, 0.0, 0.0};
    lp.kappa = {1.0, 1.0, 1.0};
    lp.alpha = {0.0, 0.0, 0.0};
    auto s = pml::ComputeStretch(lp, /*omega=*/2.0, PMLStretchFormulation::FIXED);
    REQUIRE(s[0].real() == Approx(1.0));
    REQUIRE(s[0].imag() == Approx(-0.25));  // −σ/ω = −0.5/2
    REQUIRE(s[1] == std::complex<double>(1.0, 0.0));
    REQUIRE(s[2] == std::complex<double>(1.0, 0.0));
  }

  SECTION("FREQUENCY_DEPENDENT behaves like FIXED with live ω")
  {
    pml::LocalStretchParams lp;
    lp.sigma = {0.5, 0.0, 0.0};
    lp.kappa = {1.0, 1.0, 1.0};
    auto s_fd =
        pml::ComputeStretch(lp, /*omega=*/4.0, PMLStretchFormulation::FREQUENCY_DEPENDENT);
    auto s_fixed = pml::ComputeStretch(lp, /*omega=*/4.0, PMLStretchFormulation::FIXED);
    REQUIRE(s_fd[0] == s_fixed[0]);
  }

  SECTION("CFS: s = κ + σ/(α + iω) (Palace e^{+iωt} convention)")
  {
    pml::LocalStretchParams lp;
    lp.sigma = {1.0, 0.0, 0.0};
    lp.kappa = {1.5, 1.0, 1.0};
    lp.alpha = {0.2, 0.0, 0.0};
    const double omega = 3.0;
    auto s = pml::ComputeStretch(lp, omega, PMLStretchFormulation::CFS);
    // σ / (α + i ω) = σ (α − i ω) / (α² + ω²)
    const double denom = 0.2 * 0.2 + omega * omega;
    const double re_expected = 1.5 + 1.0 * 0.2 / denom;
    const double im_expected = -1.0 * omega / denom;
    REQUIRE(s[0].real() == Approx(re_expected));
    REQUIRE(s[0].imag() == Approx(im_expected));
  }

  SECTION("CFS with α = 0 and κ = 1 reduces to FIXED-equivalent")
  {
    pml::LocalStretchParams lp;
    lp.sigma = {0.8, 0.0, 0.0};
    lp.kappa = {1.0, 1.0, 1.0};
    lp.alpha = {0.0, 0.0, 0.0};
    const double omega = 2.5;
    auto s_cfs = pml::ComputeStretch(lp, omega, PMLStretchFormulation::CFS);
    auto s_fixed = pml::ComputeStretch(lp, omega, PMLStretchFormulation::FIXED);
    REQUIRE(s_cfs[0].real() == Approx(s_fixed[0].real()));
    REQUIRE(s_cfs[0].imag() == Approx(s_fixed[0].imag()));
  }
}

TEST_CASE("PML::ComputeStretchTensors produces diagonal UPML Jacobian", "[pml][Serial]")
{
  auto p = MakeSimpleXProfile(/*sigma_max=*/1.0);

  SECTION("At inner interface ⇒ stretch is identity ⇒ Λ = I")
  {
    auto t = pml::ComputeStretchTensors(p, {1.0, 0.0, 0.0}, /*omega=*/1.0);
    for (int i = 0; i < 3; i++)
    {
      REQUIRE(t.mu_inv_re[i] == Approx(1.0));
      REQUIRE(t.mu_inv_im[i] == Approx(0.0));
      REQUIRE(t.eps_re[i] == Approx(1.0));
      REQUIRE(t.eps_im[i] == Approx(0.0));
    }
  }

  SECTION("Inside 1D x-PML: μ̃⁻¹_xx = s_x (only x is stretched)")
  {
    // 1D PML along x: s_y = s_z = 1, so the effective tensors are:
    //   μ̃⁻¹_xx = s_x/μ_r (longitudinal: weak absorption)
    //   μ̃⁻¹_yy = μ̃⁻¹_zz = 1/(μ_r·s_x) (transverse: strong absorption)
    //   ε̃_xx   = ε_r/s_x
    //   ε̃_yy   = ε̃_zz = ε_r·s_x
    // (μ_r = ε_r = 1 here.)
    const std::array<double, 3> x = {1.05, 0.0, 0.0};  // midway
    auto lp = pml::ComputeLocalStretchParams(p, x);
    auto s = pml::ComputeStretch(lp, 1.0, p.formulation);
    const std::complex<double> sx = s[0];
    const std::complex<double> inv_sx = 1.0 / sx;

    auto t = pml::ComputeStretchTensors(p, x, /*omega=*/1.0);
    // Longitudinal (xx): μ̃⁻¹_xx = s_x, ε̃_xx = 1/s_x
    REQUIRE(t.mu_inv_re[0] == Approx(sx.real()));
    REQUIRE(t.mu_inv_im[0] == Approx(sx.imag()));
    REQUIRE(t.eps_re[0] == Approx(inv_sx.real()));
    REQUIRE(t.eps_im[0] == Approx(inv_sx.imag()));
    // Transverse (yy, zz): μ̃⁻¹ = 1/s_x, ε̃ = s_x
    REQUIRE(t.mu_inv_re[1] == Approx(inv_sx.real()));
    REQUIRE(t.mu_inv_im[1] == Approx(inv_sx.imag()));
    REQUIRE(t.eps_re[1] == Approx(sx.real()));
    REQUIRE(t.eps_im[1] == Approx(sx.imag()));
  }
}

TEST_CASE("PML::ResolveSigmaMaxDefaults matches textbook formula", "[pml][Serial]")
{
  config::PMLData data;
  data.order = 3;
  data.reflection_target = 1.0e-6;
  data.thickness = {0.0, 0.0, 0.0, 0.1, 0.0, 0.0};  // +y face, d = 0.1
  data.sigma_max = {-1.0, -1.0, -1.0};              // all auto

  std::array<double, 3> sigma_max{};
  pml::ResolveSigmaMaxDefaults(data, /*mu_r=*/1.0, /*epsilon_r=*/1.0, sigma_max);

  // σ_max_y = −(n+1) ln(R) / (2 d √(εμ)) = −4 · ln(1e−6) / (2 · 0.1 · 1)
  //        = −4 · (−13.8155...) / 0.2 = 276.31...
  const double expected = -4.0 * std::log(1.0e-6) / (2.0 * 0.1 * 1.0);
  REQUIRE(sigma_max[1] == Approx(expected));

  // Inactive axes (thickness == 0) ⇒ σ_max = 0.
  REQUIRE(sigma_max[0] == Approx(0.0));
  REQUIRE(sigma_max[2] == Approx(0.0));
}

TEST_CASE("PML::BuildProfile wires config into Profile", "[pml][Serial]")
{
  config::PMLData data;
  data.order = 3;
  data.direction_signs = {0, 1, 0, 0, 0, 0};
  data.thickness = {0.0, 0.1, 0.0, 0.0, 0.0, 0.0};
  data.reflection_target = 1.0e-6;
  data.sigma_max = {-1.0, -1.0, -1.0};
  data.kappa_max = {1.0, 1.0, 1.0};
  data.alpha_max = {0.0, 0.0, 0.0};
  data.formulation = PMLStretchFormulation::FIXED;
  data.reference_frequency = 5.0;
  data.allow_refinement = false;

  auto p = pml::BuildProfile(data, /*mu_r=*/1.0, /*epsilon_r=*/1.0);

  REQUIRE(p.direction_active[1] == 1);  // +x active
  REQUIRE(p.direction_active[0] == 0);
  REQUIRE(p.thickness[0] == Approx(0.1));
  REQUIRE(p.reference_frequency == Approx(5.0));  // FIXED picks up the config value

  // Auto σ_max is resolved.
  const double expected = -4.0 * std::log(1.0e-6) / (2.0 * 0.1 * 1.0);
  REQUIRE(p.sigma_max[0] == Approx(expected));

  SECTION("FrequencyDependent leaves reference_frequency at zero")
  {
    data.formulation = PMLStretchFormulation::FREQUENCY_DEPENDENT;
    data.reference_frequency = 5.0;
    auto p2 = pml::BuildProfile(data, 1.0, 1.0);
    REQUIRE(p2.reference_frequency == Approx(0.0));
  }
}

TEST_CASE("PML::DetectSlabGeometry", "[pml][Serial]")
{
  // Physical region spans [−1, 1]³. PML attribute bboxes test each case.
  const std::array<double, 3> g_min{{-1.0, -1.0, -1.0}};
  const std::array<double, 3> g_max{{+1.0, +1.0, +1.0}};

  SECTION("+x face slab")
  {
    const std::array<double, 3> a_min{{+1.0, -1.0, -1.0}};
    const std::array<double, 3> a_max{{+1.2, +1.0, +1.0}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs == std::array<int, 6>{{0, +1, 0, 0, 0, 0}});
    REQUIRE(g.thickness[1] == Approx(0.2));  // +x thickness
    REQUIRE(g.thickness[0] == Approx(0.0));
  }

  SECTION("−y face slab")
  {
    const std::array<double, 3> a_min{{-1.0, -1.3, -1.0}};
    const std::array<double, 3> a_max{{+1.0, -1.0, +1.0}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs == std::array<int, 6>{{0, 0, -1, 0, 0, 0}});
    REQUIRE(g.thickness[2] == Approx(0.3));
  }

  SECTION("+x +y edge column (two axes active)")
  {
    const std::array<double, 3> a_min{{+1.0, +1.0, -1.0}};
    const std::array<double, 3> a_max{{+1.2, +1.2, +1.0}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs == std::array<int, 6>{{0, +1, 0, +1, 0, 0}});
    REQUIRE(g.thickness[1] == Approx(0.2));
    REQUIRE(g.thickness[3] == Approx(0.2));
    REQUIRE(g.thickness[5] == Approx(0.0));
  }

  SECTION("+x +y +z corner (three axes active)")
  {
    const std::array<double, 3> a_min{{+1.0, +1.0, +1.0}};
    const std::array<double, 3> a_max{{+1.2, +1.2, +1.2}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs == std::array<int, 6>{{0, +1, 0, +1, 0, +1}});
    REQUIRE(g.thickness[1] == Approx(0.2));
    REQUIRE(g.thickness[3] == Approx(0.2));
    REQUIRE(g.thickness[5] == Approx(0.2));
  }

  SECTION("Interior region (no faces touch bbox) ⇒ no absorption")
  {
    const std::array<double, 3> a_min{{-0.5, -0.5, -0.5}};
    const std::array<double, 3> a_max{{+0.5, +0.5, +0.5}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs == std::array<int, 6>{{0, 0, 0, 0, 0, 0}});
  }

  SECTION("Attribute spanning full x extent is inactive in x")
  {
    const std::array<double, 3> a_min{{-1.0, +1.0, -1.0}};
    const std::array<double, 3> a_max{{+1.0, +1.2, +1.0}};
    auto g = pml::DetectSlabGeometry(a_min, a_max, g_min, g_max);
    REQUIRE(g.direction_signs[0] == 0);   // not −x
    REQUIRE(g.direction_signs[1] == 0);   // not +x
    REQUIRE(g.direction_signs[3] == +1);  // +y active
  }
}

TEST_CASE("PML::ComputeSlabCentroid", "[pml][Serial]")
{
  config::PMLData data;
  const std::array<double, 3> bbmin{{-1.0, -1.0, -1.0}};
  const std::array<double, 3> bbmax{{+1.0, +1.0, +1.0}};

  SECTION("+x slab centroid lies inside the +x half of the box")
  {
    data.direction_signs = {0, 1, 0, 0, 0, 0};
    data.thickness = {0.0, 0.2, 0.0, 0.0, 0.0, 0.0};
    auto c = pml::ComputeSlabCentroid(data, bbmin, bbmax);
    REQUIRE(c[0] == Approx(0.9));  // bbmax − 0.5·thickness = 1.0 − 0.1
    REQUIRE(c[1] == Approx(0.0));  // inactive → midpoint
    REQUIRE(c[2] == Approx(0.0));
  }

  SECTION("−y slab centroid lies inside the −y half")
  {
    data.direction_signs = {0, 0, -1, 0, 0, 0};
    data.thickness = {0.0, 0.0, 0.3, 0.0, 0.0, 0.0};
    auto c = pml::ComputeSlabCentroid(data, bbmin, bbmax);
    REQUIRE(c[1] == Approx(-0.85));  // bbmin + 0.5·thickness = −1 + 0.15
  }

  SECTION("Corner active on three axes")
  {
    data.direction_signs = {0, 1, 0, 1, 0, 1};
    data.thickness = {0.0, 0.2, 0.0, 0.2, 0.0, 0.2};
    auto c = pml::ComputeSlabCentroid(data, bbmin, bbmax);
    REQUIRE(c[0] == Approx(0.9));
    REQUIRE(c[1] == Approx(0.9));
    REQUIRE(c[2] == Approx(0.9));
  }
}

}  // namespace palace
