// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "pml.hpp"

#include <cmath>
#include <limits>
#include "fem/libceed/ceed.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"

namespace palace::pml
{

std::array<double, 3> ComputeDepth(const Profile &profile, const std::array<double, 3> &x)
{
  std::array<double, 3> d{{0.0, 0.0, 0.0}};
  for (int axis = 0; axis < 3; axis++)
  {
    const double inner_neg = profile.interface_coord[axis][0];
    const double inner_pos = profile.interface_coord[axis][1];
    const bool neg_active = profile.direction_active[2 * axis + 0] != 0;
    const bool pos_active = profile.direction_active[2 * axis + 1] != 0;
    if (neg_active && x[axis] < inner_neg)
    {
      d[axis] = inner_neg - x[axis];
    }
    else if (pos_active && x[axis] > inner_pos)
    {
      d[axis] = x[axis] - inner_pos;
    }
    // else: element is on the inner side of the interface for this axis ⇒ depth 0.
  }
  return d;
}

LocalStretchParams ComputeLocalStretchParams(const Profile &profile,
                                             const std::array<double, 3> &x)
{
  LocalStretchParams lp;
  const auto depth = ComputeDepth(profile, x);
  for (int axis = 0; axis < 3; axis++)
  {
    const double d = profile.thickness[axis];
    if (d <= 0.0 || depth[axis] <= 0.0)
    {
      // Not in PML along this axis ⇒ s = 1 (σ = α = 0, κ = 1).
      lp.sigma[axis] = 0.0;
      lp.kappa[axis] = 1.0;
      lp.alpha[axis] = 0.0;
      continue;
    }
    const double t = depth[axis] / d;  // 0 at inner interface, 1 at outer.
    const double shape = std::pow(t, profile.order);
    lp.sigma[axis] = profile.sigma_max[axis] * shape;
    lp.alpha[axis] = profile.alpha_max[axis] * shape;
    // κ transitions from 1 (inner) to κ_max (outer) using the same polynomial grading.
    lp.kappa[axis] = 1.0 + (profile.kappa_max[axis] - 1.0) * shape;
  }
  return lp;
}

std::array<std::complex<double>, 3> ComputeStretch(const LocalStretchParams &local,
                                                   double omega,
                                                   PMLStretchFormulation formulation)
{
  using cdouble = std::complex<double>;
  std::array<cdouble, 3> s{{cdouble{1.0, 0.0}, cdouble{1.0, 0.0}, cdouble{1.0, 0.0}}};

  // Guard against ω = 0: the FIXED and FREQUENCY_DEPENDENT formulas divide by ω. At ω=0
  // the stretch is degenerate; callers are expected to avoid this, but we fall back to
  // s=1 (inactive PML) rather than producing NaNs.
  const double safe_omega =
      (std::abs(omega) > 0.0) ? omega : std::numeric_limits<double>::infinity();

  for (int axis = 0; axis < 3; axis++)
  {
    const double sigma = local.sigma[axis];
    const double kappa = local.kappa[axis];
    const double alpha = local.alpha[axis];
    // Palace uses the e^{+iωt} time convention. An outgoing spherical wave is
    // ∝ e^{i(ωt − kr)}, and analytic continuation r → r̃ = r − (i/ω)∫₀^r σ dr′ yields
    // decay e^{−kΣ/ω}. The resulting Jacobian is s = 1 − iσ/ω (FIXED/FREQUENCY_DEPENDENT)
    // or s = κ + σ/(α + iωε₀) (CFS, which reduces to FIXED when α → 0, κ → 1). The
    // difference between FIXED and FREQUENCY_DEPENDENT is only which ω is supplied by
    // the caller: a fixed ω₀ at setup versus the live solve frequency. The literature's
    // e^{−iωt} convention has the opposite sign (s = 1 + iσ/ω); that flip is load-bearing.
    if (formulation == PMLStretchFormulation::CFS)
    {
      const double denom = alpha * alpha + safe_omega * safe_omega;
      if (denom > 0.0)
      {
        const double re = kappa + sigma * alpha / denom;
        const double im = -sigma * safe_omega / denom;
        s[axis] = cdouble{re, im};
      }
      else
      {
        s[axis] = cdouble{kappa, 0.0};
      }
    }
    else
    {
      s[axis] = cdouble{1.0, -sigma / safe_omega};
    }
  }
  return s;
}

StretchTensors ComputeStretchTensors(const Profile &profile, const std::array<double, 3> &x,
                                     double omega)
{
  using cdouble = std::complex<double>;
  const auto local = ComputeLocalStretchParams(profile, x);
  const auto s = ComputeStretch(local, omega, profile.formulation);

  // Chew–Weedon UPML with isotropic background: ε̃ = ε_r · Λ, μ̃ = μ_r · Λ, where
  // Λ_ii = s_j · s_k / s_i (i, j, k cyclic). The curl-curl bilinear form uses μ̃⁻¹,
  // so store μ̃⁻¹_ii = (1/μ_r) · s_i / (s_j · s_k).
  StretchTensors t;
  for (int i = 0; i < 3; i++)
  {
    const int j = (i + 1) % 3;
    const int k = (i + 2) % 3;
    const cdouble eps = profile.epsilon_r * (s[j] * s[k]) / s[i];
    const cdouble mu_inv = s[i] / (profile.mu_r * s[j] * s[k]);
    t.mu_inv_re[i] = mu_inv.real();
    t.mu_inv_im[i] = mu_inv.imag();
    t.eps_re[i] = eps.real();
    t.eps_im[i] = eps.imag();
  }
  return t;
}

void ResolveSigmaMaxDefaults(const config::PMLData &data, double mu_r, double epsilon_r,
                             std::array<double, 3> &sigma_max_out)
{
  // σ_max_i = −(n + 1) · ln(R) / (2 · d_i · √(εᵣ · μᵣ)).
  // Operates in whatever unit system the Profile is in; callers pass nondimensional
  // thicknesses (data.thickness is already nondim), so the output σ_max is nondim.
  const double n_plus_1 = static_cast<double>(data.order + 1);
  const double ln_R = std::log(data.reflection_target);  // negative
  const double mat = std::sqrt(mu_r * epsilon_r);
  for (int axis = 0; axis < 3; axis++)
  {
    sigma_max_out[axis] = data.sigma_max[axis];
    if (sigma_max_out[axis] >= 0.0)
    {
      continue;  // user-specified, leave alone
    }
    // Use the maximum thickness between the negative and positive faces for this axis.
    const double d = std::max(data.thickness[2 * axis + 0], data.thickness[2 * axis + 1]);
    if (d > 0.0)
    {
      sigma_max_out[axis] = -n_plus_1 * ln_R / (2.0 * d * mat);
    }
    else
    {
      sigma_max_out[axis] = 0.0;  // inactive axis
    }
  }
}

Profile BuildProfile(const config::PMLData &data, double mu_r, double epsilon_r)
{
  Profile p;
  p.mu_r = mu_r;
  p.epsilon_r = epsilon_r;
  p.order = data.order;
  p.kappa_max = data.kappa_max;
  p.alpha_max = data.alpha_max;
  p.formulation = data.formulation;
  for (std::size_t i = 0; i < data.direction_signs.size(); i++)
  {
    p.direction_active[i] = (data.direction_signs[i] != 0) ? 1 : 0;
  }

  // Per-axis thickness: use the max of the negative and positive faces (most PML setups
  // use symmetric layers; asymmetric layers are supported but the profile uses the
  // larger d for the grading denominator).
  for (int axis = 0; axis < 3; axis++)
  {
    p.thickness[axis] =
        std::max(data.thickness[2 * axis + 0], data.thickness[2 * axis + 1]);
  }

  ResolveSigmaMaxDefaults(data, mu_r, epsilon_r, p.sigma_max);

  // Reference frequency: FIXED and CFS require a positive config value (enforced at
  // config parse time); FREQUENCY_DEPENDENT is set per solve and stays zero here.
  p.reference_frequency = (data.formulation == PMLStretchFormulation::FREQUENCY_DEPENDENT)
                              ? 0.0
                              : data.reference_frequency;

  // interface_coord is populated by the caller (pml region setup in MaterialOperator)
  // because it depends on the mesh bounding box.
  p.allow_refinement = data.allow_refinement;
  return p;
}

SlabGeometry DetectSlabGeometry(const std::array<double, 3> &attr_min,
                                const std::array<double, 3> &attr_max,
                                const std::array<double, 3> &global_min,
                                const std::array<double, 3> &global_max, double rel_tol)
{
  // A PML face on axis a is active on the + side if the attribute's bbox touches the
  // global +a boundary but not the −a boundary, and vice versa. Thickness = attr bbox
  // extent along that axis. This handles:
  //   - Outer slab touching one face: standard.
  //   - Corner/edge slabs touching multiple faces: multiple axes activated.
  //   - Attribute spanning the full axis extent (i.e., inner slab in a middle position):
  //     ambiguous, no PML active on that axis.
  // A previous version gated this on the centroid of the attribute's bbox being on one
  // side of the global centroid. That broke down when the PML occupied more than half
  // the domain on an axis (common for short guides): the PML bbox then straddled the
  // global centroid even though the axis was clearly active. Using "touches global
  // boundary" is more robust and matches users' expectation that PML attributes always
  // live on the outer boundary of the mesh.
  SlabGeometry g;
  for (int axis = 0; axis < 3; axis++)
  {
    const double extent = global_max[axis] - global_min[axis];
    const double eps = rel_tol * std::max(extent, 1.0);
    const double attr_thickness = std::max(attr_max[axis] - attr_min[axis], 0.0);
    if (attr_thickness <= eps)
    {
      continue;
    }
    const bool touches_neg = (attr_min[axis] <= global_min[axis] + eps);
    const bool touches_pos = (attr_max[axis] >= global_max[axis] - eps);
    if (touches_neg && !touches_pos)
    {
      g.direction_signs[2 * axis + 0] = -1;
      g.thickness[2 * axis + 0] = attr_thickness;
    }
    else if (touches_pos && !touches_neg)
    {
      g.direction_signs[2 * axis + 1] = +1;
      g.thickness[2 * axis + 1] = attr_thickness;
    }
    // Spans (touches both faces) or is interior (touches neither) ⇒ no PML on this axis.
  }
  return g;
}

std::array<double, 3> ComputeSlabCentroid(const config::PMLData &data,
                                          const std::array<double, 3> &bbmin,
                                          const std::array<double, 3> &bbmax)
{
  // Representative sample point for the PML attribute. For each axis with a PML face
  // active, sample the midpoint of the PML slab along that axis; otherwise sample the
  // middle of the physical extent (where σ = 0 anyway).
  //
  // With both +face and −face active on the same axis (a non-contiguous PML on one
  // attribute), this returns the midpoint of the +slab. That case is unusual; users
  // with distinct ± slabs should mesh them as separate attributes so each gets its own
  // centroid and profile.
  std::array<double, 3> c{};
  for (int axis = 0; axis < 3; axis++)
  {
    const bool neg_active = data.direction_signs[2 * axis + 0] < 0;
    const bool pos_active = data.direction_signs[2 * axis + 1] > 0;
    const double t_neg = data.thickness[2 * axis + 0];
    const double t_pos = data.thickness[2 * axis + 1];
    if (pos_active && t_pos > 0.0)
    {
      c[axis] = bbmax[axis] - 0.5 * t_pos;
    }
    else if (neg_active && t_neg > 0.0)
    {
      c[axis] = bbmin[axis] + 0.5 * t_neg;
    }
    else
    {
      c[axis] = 0.5 * (bbmin[axis] + bbmax[axis]);
    }
  }
  return c;
}

void PackProfileContext(const Profile &profile, CeedIntScalar *out)
{
  // Layout must match pml_qf.h's PMLEvalStretchTensors unpacker.
  out[0].first = static_cast<int>(profile.formulation);
  out[1].first = profile.order;
  out[2].second = profile.mu_r;
  out[3].second = profile.epsilon_r;
  out[4].second = profile.reference_frequency;
  out[5].second = profile.sigma_max[0];
  out[6].second = profile.sigma_max[1];
  out[7].second = profile.sigma_max[2];
  out[8].second = profile.kappa_max[0];
  out[9].second = profile.kappa_max[1];
  out[10].second = profile.kappa_max[2];
  out[11].second = profile.alpha_max[0];
  out[12].second = profile.alpha_max[1];
  out[13].second = profile.alpha_max[2];
  out[14].second = profile.thickness[0];
  out[15].second = profile.thickness[1];
  out[16].second = profile.thickness[2];
  // Interface coordinates are interleaved neg/pos per axis to match the unpacker.
  out[17].second = profile.interface_coord[0][0];  // -x
  out[18].second = profile.interface_coord[0][1];  // +x
  out[19].second = profile.interface_coord[1][0];  // -y
  out[20].second = profile.interface_coord[1][1];  // +y
  out[21].second = profile.interface_coord[2][0];  // -z
  out[22].second = profile.interface_coord[2][1];  // +z
  out[23].first = profile.direction_active[0];     // -x
  out[24].first = profile.direction_active[1];     // +x
  out[25].first = profile.direction_active[2];     // -y
  out[26].first = profile.direction_active[3];     // +y
  out[27].first = profile.direction_active[4];     // -z
  out[28].first = profile.direction_active[5];     // +z
}

std::vector<CeedIntScalar> PackProfileContextAll(const std::vector<int> &attr_to_profile,
                                                 const std::vector<Profile> &profiles,
                                                 double scale)
{
  const std::size_t num_attr = attr_to_profile.size();
  const std::size_t num_profiles = profiles.size();
  std::vector<CeedIntScalar> ctx(2 + num_attr + kPMLRegionStride * num_profiles);
  ctx[0].second = scale;
  ctx[1].first = static_cast<int>(num_attr);
  for (std::size_t i = 0; i < num_attr; i++)
  {
    ctx[2 + i].first = attr_to_profile[i];
  }
  for (std::size_t p = 0; p < num_profiles; p++)
  {
    PackProfileContext(profiles[p], ctx.data() + 2 + num_attr + kPMLRegionStride * p);
  }
  return ctx;
}

void SetPMLContextScale(CeedIntScalar *ctx, double scale)
{
  ctx[0].second = scale;
}

void RefreshPMLContextFrequency(CeedIntScalar *ctx, int num_attr, int num_profiles,
                                double omega)
{
  for (int p = 0; p < num_profiles; p++)
  {
    CeedIntScalar *region = ctx + 2 + num_attr + kPMLRegionStride * p;
    // region[0] is formulation, region[4] is omega (see PackProfileContext layout).
    const int formulation = region[0].first;
    if (formulation == static_cast<int>(PMLStretchFormulation::FREQUENCY_DEPENDENT))
    {
      region[4].second = omega;
    }
  }
}

}  // namespace palace::pml
