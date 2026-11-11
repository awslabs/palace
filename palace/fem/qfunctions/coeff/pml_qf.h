// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_PML_QF_H
#define PALACE_LIBCEED_PML_QF_H

#ifndef CEED_RUNNING_JIT_PASS
#include <math.h>
#endif
#include "coeff_qf.h"

// Device-callable evaluation of the Cartesian PML material tensors at a quadrature point,
// from the physical coordinate x and a packed QFunction context.
//
// With the diagonal coordinate stretch S = diag(s_x, s_y, s_z), where (Palace e^{+iωt}
// time convention)
//                       s_a(x, ω) = κ_a(x) + σ_a(x) / (α_a(x) + iω) ,
// Maxwell's equations in the stretched coordinates are equivalent to Maxwell's equations
// in the real coordinates with the transformed (complex, symmetric) material tensors
//                 μ̃⁻¹ = S μ⁻¹ S / det(S) ,     ε̃ = det(S) S⁻¹ ε S⁻¹ ,
// that is, μ̃⁻¹_ij = μ⁻¹_ij s_i s_j / (s_x s_y s_z) and ε̃_ij = ε_ij (s_x s_y s_z) / (s_i
// s_j). These hold for general (anisotropic, lossy) background tensors μ⁻¹ and ε = ε' + i
// ε''. For an isotropic background they reduce to the classic uniaxial PML tensors μ̃⁻¹_ii =
// μ⁻¹ s_i / (s_j s_k), ε̃_ii = ε s_j s_k / s_i. The frequency ω is the complex solve
// frequency for a frequency-dependent stretch (so that the stretch is the analytic
// continuation for complex eigenfrequencies), or a fixed real reference frequency ω₀ for a
// static stretch.
//
// Context layout (array of CeedIntScalar):
//
//   Header (PALACE_PML_HEADER_SIZE entries):
//     [0]       number of (1-based) libCEED attributes, num_attr
//     [1]       output part of the scaled tensor c T: 0 → Re{c T}, 1 → Im{c T},
//               3 → Re{c T̂}, 4 → Im{c T̂} for the tensor T̂ with the stretch factors
//               replaced by their magnitudes (exact for a unit stretch),
//               2 → Re{c} |T| (entrywise magnitude, for real-valued approximations)
//     [2, 3]    (Re, Im) of the prefactor c for μ̃⁻¹ terms (curl-curl, Floquet)
//     [4, 5]    (Re, Im) of the prefactor c for ε̃ terms (mass, diffusion)
//     [6, 7]    (Re, Im) of the frequency ω of the stretch: the solve frequency for a
//               frequency-dependent stretch, or the reference frequency ω₀ otherwise
//     [8..16]   Floquet wave vector cross product matrix [k ×] (3 x 3, column-major)
//   Stretch (PALACE_PML_STRETCH_SIZE entries) of the PML block of the integrator:
//     [0]       polynomial grading order n
//     [1..6]    inner interface coordinate of each face {-x, +x, -y, +y, -z, +z}
//     [7..12]   layer thickness of each face (≤ 0 for inactive faces)
//     [13..18]  σ_max of each face
//     [19..21]  κ_max per axis
//     [22..24]  α_max per axis
//   Attribute map (num_attr entries): background material index for each attribute, or -1
//   for attributes outside of the PML regions.
//   Background materials (PALACE_PML_BACKGROUND_SIZE entries each):
//     [0..8]    μ⁻¹ (3 x 3, column-major)
//     [9..17]   Re{ε}
//     [18..26]  Im{ε}
//
// The profile parameters are graded with depth d into the layer of thickness t on each face
// as σ = σ_max (d/t)ⁿ, κ = 1 + (κ_max - 1)(d/t)ⁿ, α = α_max (d/t)ⁿ.
//
// For a PML extruded from a boundary of the mesh, the stretch maps the point x at the depth
// ξ into the layer, along the unit vector d from its closest point p on a core segment, to
//                 x̃ = x + (ξ̃(ξ) - ξ) d ,     ξ̃(ξ) = ∫₀^ξ s(ξ') dξ' ,
// with the profile s(ξ) above for the depth ξ into the layer. Its Jacobian is
//                 J̃ = ∂x̃/∂x = I + (s - 1) d ∇ξᵀ + (ξ̃ - ξ) ∇d ,
// with ∇d = (I - d dᵀ) / |x - p| for p at an end of the core segment and (I - d dᵀ - e eᵀ)
// / |x - p| inside of it, for its unit direction e. The material tensors are those above
// for the general stretch J̃ in place of S,
//                 μ̃⁻¹ = J̃ᵀ μ⁻¹ J̃ / det(J̃) ,     ε̃ = det(J̃) J̃⁻¹ ε J̃⁻ᵀ ,
// which are complex symmetric, and J̃ = I at the interface with the physical region (ξ =
// 0). The depth ξ and its gradient are interpolated from the depth of the mesh nodes. On a
// sphere around a point core, J̃ = s d dᵀ + (r̃ / r) (I - d dᵀ): the spherical PML. With
// d ∇ξᵀ = n nᵀ and ∇d = 0, it is the Cartesian stretch of a face with normal n. The tensors
// with the magnitudes of the stretch factors use J̃ with s and |x - p| + ξ̃ - ξ replaced by
// their magnitudes. The extruded stretch section of the context is:
//     [25]      1 for an extruded stretch (the Cartesian entries above are then zero)
//     [26..31]  end points of the core segment (a point if they are equal)
//     [32..35]  layer thickness t, σ_max, κ_max, and α_max

#define PALACE_PML_HEADER_SIZE 17
#define PALACE_PML_STRETCH_SIZE 36
#define PALACE_PML_BACKGROUND_SIZE 27

enum
{
  PALACE_PML_PART_RE = 0,
  PALACE_PML_PART_IM = 1,
  PALACE_PML_PART_ABS = 2,
  PALACE_PML_PART_MAG_RE = 3,
  PALACE_PML_PART_MAG_IM = 4
};

CEED_QFUNCTION_HELPER CeedInt PMLNumAttr(const CeedIntScalar *ctx)
{
  return ctx[0].first;
}

CEED_QFUNCTION_HELPER CeedInt PMLPart(const CeedIntScalar *ctx)
{
  return ctx[1].first;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *PMLWaveVectorCross(const CeedIntScalar *ctx)
{
  return ctx + 8;
}

// Return the background material data for the given (1-based) libCEED attribute, or a null
// pointer if the attribute is not in a PML region.
CEED_QFUNCTION_HELPER const CeedIntScalar *PMLBackgroundData(const CeedIntScalar *ctx,
                                                             CeedInt attr)
{
  const CeedInt num_attr = PMLNumAttr(ctx);
  if (attr < 1 || attr > num_attr)
  {
    return 0;
  }
  const CeedIntScalar *attr_map = ctx + PALACE_PML_HEADER_SIZE + PALACE_PML_STRETCH_SIZE;
  const CeedInt k = attr_map[attr - 1].first;
  return (k < 0) ? 0 : attr_map + num_attr + PALACE_PML_BACKGROUND_SIZE * k;
}

CEED_QFUNCTION_HELPER CeedScalar PMLIntPow(CeedScalar x, CeedInt n)
{
  CeedScalar r = 1.0;
  for (CeedInt i = 0; i < n; i++)
  {
    r *= x;
  }
  return r;
}

CEED_QFUNCTION_HELPER void PMLComplexMult(CeedScalar ar, CeedScalar ai, CeedScalar br,
                                          CeedScalar bi, CeedScalar *cr, CeedScalar *ci)
{
  *cr = ar * br - ai * bi;
  *ci = ar * bi + ai * br;
}

CEED_QFUNCTION_HELPER void PMLComplexDiv(CeedScalar ar, CeedScalar ai, CeedScalar br,
                                         CeedScalar bi, CeedScalar *cr, CeedScalar *ci)
{
  const CeedScalar den = br * br + bi * bi;
  *cr = (ar * br + ai * bi) / den;
  *ci = (ai * br - ar * bi) / den;
}

// Stretch factor s = κ + σ / (α + iω) at the normalized depth r ∈ [0, 1] into the layer,
// for the frequency ω of the context.
CEED_QFUNCTION_HELPER void PMLStretchFactor(const CeedIntScalar *ctx, CeedScalar r,
                                            CeedScalar sigma_max, CeedScalar kappa_max,
                                            CeedScalar alpha_max, CeedScalar *s_re,
                                            CeedScalar *s_im)
{
  const CeedInt order = ctx[PALACE_PML_HEADER_SIZE].first;
  const CeedScalar omega_re = ctx[6].second, omega_im = ctx[7].second;
  const CeedScalar shape = PMLIntPow((r < 1.0) ? r : 1.0, order);
  const CeedScalar sigma = sigma_max * shape;
  const CeedScalar kappa = 1.0 + (kappa_max - 1.0) * shape;
  const CeedScalar alpha = alpha_max * shape;

  // σ / (α + iω) with α + iω = (α - Im{ω}) + i Re{ω}. The degenerate α = ω = 0 case falls
  // back to the real coordinate scaling κ.
  const CeedScalar den_re = alpha - omega_im, den_im = omega_re;
  *s_re = kappa;
  *s_im = 0.0;
  if (den_re != 0.0 || den_im != 0.0)
  {
    CeedScalar q_re, q_im;
    PMLComplexDiv(sigma, 0.0, den_re, den_im, &q_re, &q_im);
    *s_re += q_re;
    *s_im = q_im;
  }
}

// Compute the complex stretch factors s_a for each axis at the physical point x.
CEED_QFUNCTION_HELPER void PMLStretch(const CeedIntScalar *ctx, const CeedScalar x[3],
                                      CeedScalar s_re[3], CeedScalar s_im[3])
{
  const CeedIntScalar *p = ctx + PALACE_PML_HEADER_SIZE;
  for (CeedInt a = 0; a < 3; a++)
  {
    s_re[a] = 1.0;
    s_im[a] = 0.0;
    const CeedScalar t_neg = p[7 + 2 * a].second, t_pos = p[8 + 2 * a].second;
    const CeedScalar d_neg = p[1 + 2 * a].second - x[a], d_pos = x[a] - p[2 + 2 * a].second;
    CeedScalar r, sigma_max;
    if (t_neg > 0.0 && d_neg > 0.0)
    {
      r = d_neg / t_neg;
      sigma_max = p[13 + 2 * a].second;
    }
    else if (t_pos > 0.0 && d_pos > 0.0)
    {
      r = d_pos / t_pos;
      sigma_max = p[14 + 2 * a].second;
    }
    else
    {
      continue;
    }
    PMLStretchFactor(ctx, r, sigma_max, p[19 + a].second, p[22 + a].second, &s_re[a],
                     &s_im[a]);
  }
}

// Requested part of the complex prefactor c times a complex tensor entry t: Re{c t},
// Im{c t}, or the real-valued approximation Re{c} |t|.
CEED_QFUNCTION_HELPER CeedScalar PMLScaledPart(CeedInt part, CeedScalar c_re,
                                               CeedScalar c_im, CeedScalar t_re,
                                               CeedScalar t_im)
{
  switch (part)
  {
    case PALACE_PML_PART_IM:
    case PALACE_PML_PART_MAG_IM:
      return c_re * t_im + c_im * t_re;
    case PALACE_PML_PART_ABS:
      return c_re * sqrt(t_re * t_re + t_im * t_im);
    default:
      return c_re * t_re - c_im * t_im;
  }
}

// Stretch factor f = f_re + i f_im of a tensor entry, replaced by its magnitude for the
// parts with the magnitudes of the stretch factors.
CEED_QFUNCTION_HELPER void PMLStretchFactorPart(CeedInt part, CeedScalar *f_re,
                                                CeedScalar *f_im)
{
  if (part == PALACE_PML_PART_MAG_RE || part == PALACE_PML_PART_MAG_IM)
  {
    *f_re = sqrt(*f_re * *f_re + *f_im * *f_im);
    *f_im = 0.0;
  }
}

// Stretch factors s_a and their product det(S) = s_x s_y s_z at x.
CEED_QFUNCTION_HELPER void PMLStretchDet(const CeedIntScalar *ctx, const CeedScalar x[3],
                                         CeedScalar s_re[3], CeedScalar s_im[3],
                                         CeedScalar *det_re, CeedScalar *det_im)
{
  PMLStretch(ctx, x, s_re, s_im);
  PMLComplexMult(s_re[0], s_im[0], s_re[1], s_im[1], det_re, det_im);
  PMLComplexMult(*det_re, *det_im, s_re[2], s_im[2], det_re, det_im);
}

// Requested part of c μ̃⁻¹ (3 x 3, column-major), with μ̃⁻¹_ij = μ⁻¹_ij s_i s_j / det(S) for
// the background material b and c the μ̃⁻¹ prefactor in the context header.
CEED_QFUNCTION_HELPER void PMLMuInvCoeffStretch(const CeedIntScalar *ctx,
                                                const CeedIntScalar *b,
                                                const CeedScalar s_re[3],
                                                const CeedScalar s_im[3], CeedScalar det_re,
                                                CeedScalar det_im, CeedScalar coeff[9])
{
  const CeedInt part = PMLPart(ctx);
  const CeedScalar c_re = ctx[2].second, c_im = ctx[3].second;
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      const CeedScalar mu_inv = b[i + 3 * j].second;
      if (mu_inv == 0.0)
      {
        coeff[i + 3 * j] = 0.0;
        continue;
      }
      CeedScalar f_re, f_im;
      PMLComplexMult(s_re[i], s_im[i], s_re[j], s_im[j], &f_re, &f_im);
      PMLComplexDiv(f_re, f_im, det_re, det_im, &f_re, &f_im);
      PMLStretchFactorPart(part, &f_re, &f_im);
      coeff[i + 3 * j] = PMLScaledPart(part, c_re, c_im, mu_inv * f_re, mu_inv * f_im);
    }
  }
}

// Requested part of c ε̃ (3 x 3, column-major), with ε̃_ij = ε_ij det(S) / (s_i s_j) for the
// complex permittivity ε_ij = ε'_ij + i ε''_ij of the background material b and c the ε̃
// prefactor in the context header.
CEED_QFUNCTION_HELPER void PMLEpsCoeffStretch(const CeedIntScalar *ctx,
                                              const CeedIntScalar *b,
                                              const CeedScalar s_re[3],
                                              const CeedScalar s_im[3], CeedScalar det_re,
                                              CeedScalar det_im, CeedScalar coeff[9])
{
  const CeedInt part = PMLPart(ctx);
  const CeedScalar c_re = ctx[4].second, c_im = ctx[5].second;
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      const CeedScalar eps_re = b[9 + i + 3 * j].second, eps_im = b[18 + i + 3 * j].second;
      if (eps_re == 0.0 && eps_im == 0.0)
      {
        coeff[i + 3 * j] = 0.0;
        continue;
      }
      CeedScalar f_re, f_im, t_re, t_im;
      PMLComplexMult(s_re[i], s_im[i], s_re[j], s_im[j], &f_re, &f_im);
      PMLComplexDiv(det_re, det_im, f_re, f_im, &f_re, &f_im);
      PMLStretchFactorPart(part, &f_re, &f_im);
      PMLComplexMult(eps_re, eps_im, f_re, f_im, &t_re, &t_im);
      coeff[i + 3 * j] = PMLScaledPart(part, c_re, c_im, t_re, t_im);
    }
  }
}

// Requested part of c μ̃⁻¹ at x. Returns false (and leaves coeff unset) for attributes
// outside of the PML regions.
CEED_QFUNCTION_HELPER bool PMLMuInvCoeff(const CeedIntScalar *ctx, CeedInt attr,
                                         const CeedScalar x[3], CeedScalar coeff[9])
{
  const CeedIntScalar *b = PMLBackgroundData(ctx, attr);
  if (!b)
  {
    return false;
  }
  CeedScalar s_re[3], s_im[3], det_re, det_im;
  PMLStretchDet(ctx, x, s_re, s_im, &det_re, &det_im);
  PMLMuInvCoeffStretch(ctx, b, s_re, s_im, det_re, det_im, coeff);
  return true;
}

// Requested part of c ε̃ at x. Returns false (and leaves coeff unset) for attributes outside
// of the PML regions.
CEED_QFUNCTION_HELPER bool PMLEpsCoeff(const CeedIntScalar *ctx, CeedInt attr,
                                       const CeedScalar x[3], CeedScalar coeff[9])
{
  const CeedIntScalar *b = PMLBackgroundData(ctx, attr);
  if (!b)
  {
    return false;
  }
  CeedScalar s_re[3], s_im[3], det_re, det_im;
  PMLStretchDet(ctx, x, s_re, s_im, &det_re, &det_im);
  PMLEpsCoeffStretch(ctx, b, s_re, s_im, det_re, det_im, coeff);
  return true;
}

// Both of the above, sharing the stretch evaluation.
CEED_QFUNCTION_HELPER bool PMLMuInvEpsCoeff(const CeedIntScalar *ctx, CeedInt attr,
                                            const CeedScalar x[3], CeedScalar mu_coeff[9],
                                            CeedScalar eps_coeff[9])
{
  const CeedIntScalar *b = PMLBackgroundData(ctx, attr);
  if (!b)
  {
    return false;
  }
  CeedScalar s_re[3], s_im[3], det_re, det_im;
  PMLStretchDet(ctx, x, s_re, s_im, &det_re, &det_im);
  PMLMuInvCoeffStretch(ctx, b, s_re, s_im, det_re, det_im, mu_coeff);
  PMLEpsCoeffStretch(ctx, b, s_re, s_im, det_re, det_im, eps_coeff);
  return true;
}

// Whether the stretch of the context is that of an extruded PML.
CEED_QFUNCTION_HELPER bool PMLIsExtruded(const CeedIntScalar *ctx)
{
  return ctx[PALACE_PML_HEADER_SIZE + 25].first == 1;
}

// ξ̃ - ξ = ∫₀^ξ (s - 1) dξ' for the extruded stretch: in closed form without the complex
// frequency shift, and otherwise with 8-point Gauss-Legendre quadrature. Depths beyond the
// layer (from rounding) have the stretch of its outer edge.
CEED_QFUNCTION_HELPER void PMLStretchIntegral(const CeedIntScalar *ctx, CeedScalar xi,
                                              CeedScalar *re, CeedScalar *im)
{
  const CeedIntScalar *p = ctx + PALACE_PML_HEADER_SIZE;
  const CeedInt order = p[0].first;
  const CeedScalar t = p[32].second, sigma_max = p[33].second, kappa_max = p[34].second,
                   alpha_max = p[35].second;
  const CeedScalar omega_re = ctx[6].second, omega_im = ctx[7].second;
  const CeedScalar r = (xi <= 0.0) ? 0.0 : (xi < t) ? xi / t : 1.0;
  const CeedScalar beyond = (xi > t) ? xi - t : 0.0;
  const CeedScalar m = t * PMLIntPow(r, order + 1) / (order + 1) + beyond;  // ∫ (ξ'/t)ⁿ
  *re = (kappa_max - 1.0) * m;
  *im = 0.0;
  CeedScalar q_re = 0.0, q_im = 0.0;
  if (alpha_max == 0.0)
  {
    if (omega_re != 0.0 || omega_im != 0.0)
    {
      PMLComplexDiv(sigma_max * m, 0.0, -omega_im, omega_re, &q_re, &q_im);
    }
  }
  else
  {
    const CeedScalar gl_x[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                0.9602898564975363};
    const CeedScalar gl_w[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                                0.1012285362903763};
    for (CeedInt k = 0; k < 8; k++)
    {
      const CeedScalar u = 0.5 * r * (1.0 + ((k < 4) ? -gl_x[k] : gl_x[k - 4]));
      const CeedScalar w = 0.5 * r * t * gl_w[(k < 4) ? k : k - 4];
      CeedScalar s_re, s_im;
      PMLStretchFactor(ctx, u, sigma_max, 1.0, alpha_max, &s_re, &s_im);
      q_re += w * (s_re - 1.0);
      q_im += w * s_im;
    }
    CeedScalar s_re, s_im;
    PMLStretchFactor(ctx, 1.0, sigma_max, 1.0, alpha_max, &s_re, &s_im);
    q_re += beyond * (s_re - 1.0);
    q_im += beyond * s_im;
  }
  *re += q_re;
  *im += q_im;
}

// Unit vector d from the closest point p of x on the core segment of the extruded stretch,
// the distance |x - p|, and ∇d (3 x 3, column-major, (∇d)_ij = ∂d_i / ∂x_j).
CEED_QFUNCTION_HELPER void PMLRay(const CeedIntScalar *ctx, const CeedScalar x[3],
                                  CeedScalar d[3], CeedScalar *dist, CeedScalar grad_d[9])
{
  const CeedIntScalar *a = ctx + PALACE_PML_HEADER_SIZE + 26, *b = a + 3;
  CeedScalar e[3], ee = 0.0, t = 0.0;
  for (CeedInt c = 0; c < 3; c++)
  {
    e[c] = b[c].second - a[c].second;
    ee += e[c] * e[c];
    t += (x[c] - a[c].second) * e[c];
  }
  t = (ee > 0.0) ? t / ee : 0.0;
  const bool inside = (t > 0.0 && t < 1.0);
  t = (t < 0.0) ? 0.0 : (t > 1.0) ? 1.0 : t;
  for (CeedInt c = 0; c < 3; c++)
  {
    d[c] = x[c] - a[c].second - t * e[c];
  }
  *dist = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  for (CeedInt c = 0; c < 3; c++)
  {
    d[c] /= *dist;
  }
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      grad_d[i + 3 * j] =
          ((i == j) - d[i] * d[j] - (inside ? e[i] * e[j] / ee : 0.0)) / *dist;
    }
  }
}

// Jacobian J̃ (complex 3 x 3, column-major) of the extruded stretch at x, for the depth ξ
// and its gradient, or the one with the magnitudes of the stretch factors.
CEED_QFUNCTION_HELPER void PMLExtrudedJacobian(const CeedIntScalar *ctx,
                                               const CeedScalar x[3], CeedScalar xi,
                                               const CeedScalar grad_xi[3], bool magnitude,
                                               CeedScalar J_re[9], CeedScalar J_im[9])
{
  const CeedIntScalar *p = ctx + PALACE_PML_HEADER_SIZE;
  const CeedScalar t = p[32].second;
  CeedScalar s_re, s_im, D_re, D_im, d[3], dist, grad_d[9];
  PMLStretchFactor(ctx, (xi > 0.0) ? xi / t : 0.0, p[33].second, p[34].second, p[35].second,
                   &s_re, &s_im);
  PMLStretchIntegral(ctx, xi, &D_re, &D_im);
  PMLRay(ctx, x, d, &dist, grad_d);
  if (magnitude)
  {
    s_re = sqrt(s_re * s_re + s_im * s_im);
    s_im = 0.0;
    D_re = sqrt((dist + D_re) * (dist + D_re) + D_im * D_im) - dist;
    D_im = 0.0;
  }
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      J_re[i + 3 * j] =
          (i == j) + (s_re - 1.0) * d[i] * grad_xi[j] + D_re * grad_d[i + 3 * j];
      J_im[i + 3 * j] = s_im * d[i] * grad_xi[j] + D_im * grad_d[i + 3 * j];
    }
  }
}

// Adjugate B = adj(A) and determinant of a complex 3 x 3 matrix A (column-major).
CEED_QFUNCTION_HELPER void PMLComplexAdj33(const CeedScalar A_re[9],
                                           const CeedScalar A_im[9], CeedScalar B_re[9],
                                           CeedScalar B_im[9], CeedScalar *det_re,
                                           CeedScalar *det_im)
{
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      // adj(A)_ij = A_(j+1)(i+1) A_(j+2)(i+2) - A_(j+1)(i+2) A_(j+2)(i+1).
      const CeedInt j1 = (j + 1) % 3, j2 = (j + 2) % 3, i1 = (i + 1) % 3, i2 = (i + 2) % 3;
      CeedScalar p_re, p_im, q_re, q_im;
      PMLComplexMult(A_re[j1 + 3 * i1], A_im[j1 + 3 * i1], A_re[j2 + 3 * i2],
                     A_im[j2 + 3 * i2], &p_re, &p_im);
      PMLComplexMult(A_re[j1 + 3 * i2], A_im[j1 + 3 * i2], A_re[j2 + 3 * i1],
                     A_im[j2 + 3 * i1], &q_re, &q_im);
      B_re[i + 3 * j] = p_re - q_re;
      B_im[i + 3 * j] = p_im - q_im;
    }
  }
  *det_re = *det_im = 0.0;
  for (CeedInt k = 0; k < 3; k++)
  {
    CeedScalar p_re, p_im;
    PMLComplexMult(A_re[3 * k], A_im[3 * k], B_re[k], B_im[k], &p_re, &p_im);
    *det_re += p_re;
    *det_im += p_im;
  }
}

// Requested part of c μ̃⁻¹ (3 x 3, column-major), with μ̃⁻¹ = J̃ᵀ μ⁻¹ J̃ / det(J̃) for the
// background material b and c the μ̃⁻¹ prefactor in the context header.
CEED_QFUNCTION_HELPER void PMLMuInvCoeffJacobian(const CeedIntScalar *ctx,
                                                 const CeedIntScalar *b,
                                                 const CeedScalar J_re[9],
                                                 const CeedScalar J_im[9],
                                                 CeedScalar coeff[9])
{
  const CeedInt part = PMLPart(ctx);
  const CeedScalar c_re = ctx[2].second, c_im = ctx[3].second;
  CeedScalar adj_re[9], adj_im[9], det_re, det_im, M_re[9], M_im[9];
  PMLComplexAdj33(J_re, J_im, adj_re, adj_im, &det_re, &det_im);
  for (CeedInt j = 0; j < 3; j++)  // M = μ⁻¹ J̃
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      M_re[i + 3 * j] = M_im[i + 3 * j] = 0.0;
      for (CeedInt k = 0; k < 3; k++)
      {
        M_re[i + 3 * j] += b[i + 3 * k].second * J_re[k + 3 * j];
        M_im[i + 3 * j] += b[i + 3 * k].second * J_im[k + 3 * j];
      }
    }
  }
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      CeedScalar t_re = 0.0, t_im = 0.0;
      for (CeedInt k = 0; k < 3; k++)
      {
        CeedScalar p_re, p_im;
        PMLComplexMult(J_re[k + 3 * i], J_im[k + 3 * i], M_re[k + 3 * j], M_im[k + 3 * j],
                       &p_re, &p_im);
        t_re += p_re;
        t_im += p_im;
      }
      PMLComplexDiv(t_re, t_im, det_re, det_im, &t_re, &t_im);
      coeff[i + 3 * j] = PMLScaledPart(part, c_re, c_im, t_re, t_im);
    }
  }
}

// Requested part of c ε̃ (3 x 3, column-major), with ε̃ = det(J̃) J̃⁻¹ ε J̃⁻ᵀ = adj(J̃) ε
// adj(J̃)ᵀ / det(J̃) for the complex permittivity of the background material b and c the ε̃
// prefactor in the context header.
CEED_QFUNCTION_HELPER void
PMLEpsCoeffJacobian(const CeedIntScalar *ctx, const CeedIntScalar *b,
                    const CeedScalar J_re[9], const CeedScalar J_im[9], CeedScalar coeff[9])
{
  const CeedInt part = PMLPart(ctx);
  const CeedScalar c_re = ctx[4].second, c_im = ctx[5].second;
  CeedScalar adj_re[9], adj_im[9], det_re, det_im, M_re[9], M_im[9];
  PMLComplexAdj33(J_re, J_im, adj_re, adj_im, &det_re, &det_im);
  for (CeedInt j = 0; j < 3; j++)  // M = ε adj(J̃)ᵀ
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      M_re[i + 3 * j] = M_im[i + 3 * j] = 0.0;
      for (CeedInt k = 0; k < 3; k++)
      {
        CeedScalar p_re, p_im;
        PMLComplexMult(b[9 + i + 3 * k].second, b[18 + i + 3 * k].second, adj_re[j + 3 * k],
                       adj_im[j + 3 * k], &p_re, &p_im);
        M_re[i + 3 * j] += p_re;
        M_im[i + 3 * j] += p_im;
      }
    }
  }
  for (CeedInt j = 0; j < 3; j++)
  {
    for (CeedInt i = 0; i < 3; i++)
    {
      CeedScalar t_re = 0.0, t_im = 0.0;
      for (CeedInt k = 0; k < 3; k++)
      {
        CeedScalar p_re, p_im;
        PMLComplexMult(adj_re[i + 3 * k], adj_im[i + 3 * k], M_re[k + 3 * j],
                       M_im[k + 3 * j], &p_re, &p_im);
        t_re += p_re;
        t_im += p_im;
      }
      PMLComplexDiv(t_re, t_im, det_re, det_im, &t_re, &t_im);
      coeff[i + 3 * j] = PMLScaledPart(part, c_re, c_im, t_re, t_im);
    }
  }
}

// Requested parts of c μ̃⁻¹ and/or c ε̃ (either may be null) at a quadrature point with the
// coordinates x, for the Cartesian stretch, or for the extruded one at the depth ξ with the
// gradient grad_xi. Returns false (and leaves the coefficients unset) for attributes
// outside of the PML regions.
CEED_QFUNCTION_HELPER bool PMLCoeffAt(const CeedIntScalar *ctx, CeedInt attr,
                                      const CeedScalar x[3], CeedScalar xi,
                                      const CeedScalar grad_xi[3], CeedScalar mu_coeff[9],
                                      CeedScalar eps_coeff[9])
{
  if (!PMLIsExtruded(ctx))
  {
    return (mu_coeff && eps_coeff) ? PMLMuInvEpsCoeff(ctx, attr, x, mu_coeff, eps_coeff)
           : mu_coeff              ? PMLMuInvCoeff(ctx, attr, x, mu_coeff)
                                   : PMLEpsCoeff(ctx, attr, x, eps_coeff);
  }
  const CeedIntScalar *b = PMLBackgroundData(ctx, attr);
  if (!b)
  {
    return false;
  }
  const CeedInt part = PMLPart(ctx);
  CeedScalar J_re[9], J_im[9];
  PMLExtrudedJacobian(ctx, x, xi, grad_xi,
                      part == PALACE_PML_PART_MAG_RE || part == PALACE_PML_PART_MAG_IM,
                      J_re, J_im);
  if (mu_coeff)
  {
    PMLMuInvCoeffJacobian(ctx, b, J_re, J_im, mu_coeff);
  }
  if (eps_coeff)
  {
    PMLEpsCoeffJacobian(ctx, b, J_re, J_im, eps_coeff);
  }
  return true;
}

#endif  // PALACE_LIBCEED_PML_QF_H
