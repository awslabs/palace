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

#define PALACE_PML_HEADER_SIZE 17
#define PALACE_PML_STRETCH_SIZE 25
#define PALACE_PML_BACKGROUND_SIZE 27

enum
{
  PALACE_PML_PART_RE = 0,
  PALACE_PML_PART_IM = 1,
  PALACE_PML_PART_ABS = 2
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

// Compute the complex stretch factors s_a for each axis at the physical point x.
CEED_QFUNCTION_HELPER void PMLStretch(const CeedIntScalar *ctx, const CeedScalar x[3],
                                      CeedScalar s_re[3], CeedScalar s_im[3])
{
  const CeedIntScalar *p = ctx + PALACE_PML_HEADER_SIZE;
  const CeedInt order = p[0].first;
  const CeedScalar omega_re = ctx[6].second, omega_im = ctx[7].second;
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
    const CeedScalar shape = PMLIntPow((r < 1.0) ? r : 1.0, order);
    const CeedScalar sigma = sigma_max * shape;
    const CeedScalar kappa = 1.0 + (p[19 + a].second - 1.0) * shape;
    const CeedScalar alpha = p[22 + a].second * shape;

    // σ / (α + iω) with α + iω = (α - Im{ω}) + i Re{ω}. The degenerate α = ω = 0 case falls
    // back to the real coordinate scaling κ.
    const CeedScalar den_re = alpha - omega_im, den_im = omega_re;
    s_re[a] = kappa;
    if (den_re != 0.0 || den_im != 0.0)
    {
      CeedScalar q_re, q_im;
      PMLComplexDiv(sigma, 0.0, den_re, den_im, &q_re, &q_im);
      s_re[a] += q_re;
      s_im[a] = q_im;
    }
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
      return c_re * t_im + c_im * t_re;
    case PALACE_PML_PART_ABS:
      return c_re * sqrt(t_re * t_re + t_im * t_im);
    default:
      return c_re * t_re - c_im * t_im;
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

#endif  // PALACE_LIBCEED_PML_QF_H
