// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_PML_QF_H
#define PALACE_LIBCEED_PML_QF_H

#ifndef CEED_RUNNING_JIT_PASS
#include <math.h>
#endif
#include "coeff_qf.h"

// Device-callable evaluation of the PML material tensors at a quadrature point, from the
// physical coordinate x and a packed QFunction context.
//
// For a coordinate stretch with the Jacobian J, Maxwell's equations in the stretched
// coordinates are equivalent to Maxwell's equations in the real coordinates with the
// transformed (complex, symmetric) material tensors
//                 μ̃⁻¹ = Jᵀ μ⁻¹ J / det(J) ,     ε̃ = det(J) J⁻¹ ε J⁻ᵀ ,
// for general (anisotropic, lossy) background tensors μ⁻¹ and ε = ε' + i ε''. The Cartesian
// stretch has J = S = diag(s_x, s_y, s_z), where (Palace e^{+iωt} time convention)
//                       s_a(x, ω) = κ_a(x) + σ_a(x) / (α_a(x) + iω) ,
// that is, μ̃⁻¹_ij = μ⁻¹_ij s_i s_j / (s_x s_y s_z) and ε̃_ij = ε_ij (s_x s_y s_z) / (s_i
// s_j): for an isotropic background, the classic uniaxial PML tensors μ̃⁻¹_ii = μ⁻¹ s_i /
// (s_j s_k), ε̃_ii = ε s_j s_k / s_i. The frequency ω is the complex solve
// frequency for a frequency-dependent stretch (so that the stretch is the analytic
// continuation for complex eigenfrequencies), or a fixed real reference frequency ω₀ for a
// static stretch.
//
// Context layout (array of CeedIntScalar):
//
//   Header (PALACE_PML_HEADER_SIZE entries):
//     [0]       number of (1-based) libCEED attributes, num_attr
//     [1]       output part of the scaled tensor c T: 0 → Re{c T}, 1 → Im{c T},
//               2 → Re{c T̂}, 3 → Im{c T̂} for the tensor T̂ with the stretch factors
//               replaced by their magnitudes (exact for a unit stretch)
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
// as σ = σ_max (d/t)ⁿ, κ = 1 + (κ_max - 1)(d/t)ⁿ, α = α_max (d/t)ⁿ. The tensors T̂ of the
// magnitudes of the stretch factors use J with the stretch factors replaced by their
// magnitudes.
//
// For a PML extruded from a boundary of the mesh, the stretch maps the point x at the depth
// ξ into the layer, along the unit direction d of the extrusion, to
//                 x̃ = x + (ξ̃(ξ) - ξ) d ,     ξ̃(ξ) = ∫₀^ξ s(ξ') dξ' ,
// with the profile s(ξ) above for the depth ξ into the layer. Its Jacobian is
//                 J̃ = ∂x̃/∂x = I + (s - 1) d ∇ξᵀ + (ξ̃ - ξ) ∇d .
// The material tensors are those above for J = J̃, and J̃ = I at the interface with the
// physical region (ξ = 0). The depth is a closed form of x, for one of two shapes of the
// surface the layer is extruded from (the layer of a box extruded per component is the
// Cartesian PML above):
//   - Normal extrusion of a flat surface on the plane n · x = h: d = n, ξ = n · x - h, and
//     J̃ = I + (s - 1) n nᵀ, the Cartesian stretch along n.
//   - Radial extrusion of a surface at the distance r₀ from its core (a sphere around a
//     point, the wall of a cylinder around a segment): d = (x - p) / |x - p| for the
//     closest point p of x on the core, an axis-aligned box (a point, a segment, or a
//     rectangle where its corners share coordinates), ξ = |x - p| - r₀, and ∇ξ = d. For
//     the diagonal P = ∂(x - p)/∂x, one along the axes on which x is outside of the core
//     and zero along the others, ∇d = (P - d dᵀ) / |x - p|.
// Negative depths, of points of a sampled curved surface inside of r₀, are zero. On a
// sphere around a point core, J̃ = s d dᵀ + (r̃ / r) (I - d dᵀ): the spherical PML. For the
// tensors T̂, s and |x - p| + ξ̃ - ξ in J̃ are replaced by their magnitudes. The extruded
// stretch section of the context is:
//     [25]      form of the depth: PALACE_PML_CARTESIAN (then the entries below are zero,
//               and otherwise the Cartesian entries [1..24] above), PALACE_PML_NORMAL, or
//               PALACE_PML_DISTANCE
//     [26..31]  lower and upper corners of the core box
//     [32..35]  layer thickness t, σ_max, κ_max, and α_max
//     [36]      r₀ for PALACE_PML_DISTANCE
//     [37..40]  the plane (n, h) of the surface for PALACE_PML_NORMAL

#define PALACE_PML_HEADER_SIZE 17
#define PALACE_PML_STRETCH_SIZE 41
#define PALACE_PML_BACKGROUND_SIZE 27

enum
{
  PALACE_PML_PART_RE = 0,
  PALACE_PML_PART_IM = 1,
  PALACE_PML_PART_MAG_RE = 2,
  PALACE_PML_PART_MAG_IM = 3
};

enum
{
  PALACE_PML_CARTESIAN = 0,
  PALACE_PML_NORMAL = 1,
  PALACE_PML_DISTANCE = 2
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

// Requested part of the complex prefactor c times a complex tensor entry t: Re{c t} or
// Im{c t} (also for the real tensor entries with the magnitudes of the stretch factors).
CEED_QFUNCTION_HELPER CeedScalar PMLScaledPart(CeedInt part, CeedScalar c_re,
                                               CeedScalar c_im, CeedScalar t_re,
                                               CeedScalar t_im)
{
  return (part == PALACE_PML_PART_IM || part == PALACE_PML_PART_MAG_IM)
             ? c_re * t_im + c_im * t_re
             : c_re * t_re - c_im * t_im;
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

// Requested part of c μ̃⁻¹ (3 x 3, column-major), with μ̃⁻¹ = Jᵀ μ⁻¹ J / det(J) for the
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
  for (CeedInt j = 0; j < 3; j++)  // M = μ⁻¹ J
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

// Requested part of c ε̃ (3 x 3, column-major), with ε̃ = det(J) J⁻¹ ε J⁻ᵀ = adj(J) ε
// adj(J)ᵀ / det(J) for the complex permittivity of the background material b and c the ε̃
// prefactor in the context header.
CEED_QFUNCTION_HELPER void
PMLEpsCoeffJacobian(const CeedIntScalar *ctx, const CeedIntScalar *b,
                    const CeedScalar J_re[9], const CeedScalar J_im[9], CeedScalar coeff[9])
{
  const CeedInt part = PMLPart(ctx);
  const CeedScalar c_re = ctx[4].second, c_im = ctx[5].second;
  CeedScalar adj_re[9], adj_im[9], det_re, det_im, M_re[9], M_im[9];
  PMLComplexAdj33(J_re, J_im, adj_re, adj_im, &det_re, &det_im);
  for (CeedInt j = 0; j < 3; j++)  // M = ε adj(J)ᵀ
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

// Jacobian J = diag(s_x, s_y, s_z) (complex 3 x 3, column-major) of the Cartesian stretch
// at x, or the one with the magnitudes of the stretch factors.
CEED_QFUNCTION_HELPER void PMLCartesianJacobian(const CeedIntScalar *ctx,
                                                const CeedScalar x[3], bool magnitude,
                                                CeedScalar J_re[9], CeedScalar J_im[9])
{
  CeedScalar s_re[3], s_im[3];
  PMLStretch(ctx, x, s_re, s_im);
  for (CeedInt k = 0; k < 9; k++)
  {
    J_re[k] = J_im[k] = 0.0;
  }
  for (CeedInt a = 0; a < 3; a++)
  {
    J_re[4 * a] = magnitude ? sqrt(s_re[a] * s_re[a] + s_im[a] * s_im[a]) : s_re[a];
    J_im[4 * a] = magnitude ? 0.0 : s_im[a];
  }
}

// Whether the stretch of the context is that of an extruded PML.
CEED_QFUNCTION_HELPER bool PMLIsExtruded(const CeedIntScalar *ctx)
{
  return ctx[PALACE_PML_HEADER_SIZE + 25].first != PALACE_PML_CARTESIAN;
}

// ξ̃ - ξ = ∫₀^ξ (s - 1) dξ' for the extruded stretch, to roundoff for any shift and order.
// Depths beyond the layer (from rounding) have the stretch of its outer edge.
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
    // σ_max t F for F = ∫₀^r uⁿ / (a uⁿ + b) du, a = α_max, b = iω, with p = a rⁿ / b: for
    // |p| < 1/4 the series F = (r rⁿ / b) Σ_j (-p)ʲ / (n (j + 1) + 1) to |p|ʲ < 1e-17 (at
    // most 25 terms), and otherwise partial fractions over the roots z_k of zⁿ = -b / a,
    //   F = (1 / a) (r + (1 / n) Σ_k z_k log(1 - r / z_k)) ,
    // which loses a factor of about n / |p| to cancellation. On the path [0, r] the
    // principal logarithm crosses no branch cut unless s is infinite on it. b = 0 is
    // s - 1 = σ / α.
    const CeedScalar b_re = -omega_im, b_im = omega_re,
                     b_abs = sqrt(b_re * b_re + b_im * b_im);
    const CeedScalar rn = PMLIntPow(r, order);
    if (alpha_max * rn < 0.25 * b_abs)
    {
      CeedScalar p_re, p_im, S_re = 0.0, S_im = 0.0, T_re = 1.0, T_im = 0.0;
      PMLComplexDiv(-alpha_max * rn, 0.0, b_re, b_im, &p_re, &p_im);
      for (CeedInt j = 0; j < 25 && T_re * T_re + T_im * T_im > 1.0e-34; j++)
      {
        const CeedScalar c = 1.0 / (order * (j + 1) + 1);
        S_re += c * T_re;
        S_im += c * T_im;
        PMLComplexMult(T_re, T_im, p_re, p_im, &T_re, &T_im);
      }
      PMLComplexDiv(r * rn, 0.0, b_re, b_im, &q_re, &q_im);
      PMLComplexMult(q_re, q_im, S_re, S_im, &q_re, &q_im);
    }
    else if (b_abs == 0.0)
    {
      q_re = r / alpha_max;
    }
    else
    {
      // z_k = |z| e^{iφ_k}, φ_k = (arg(-b) + 2πk) / n, and 1 - r / z_k = 1 - ρ e^{-iφ_k}
      // for ρ = r / |z|.
      const CeedScalar z = pow(b_abs / alpha_max, 1.0 / order), rho = r / z,
                       theta = atan2(-b_im, -b_re);
      CeedScalar L_re = order * rho, L_im = 0.0;
      for (CeedInt k = 0; k < order; k++)
      {
        const CeedScalar phi = (theta + 6.283185307179586 * k) / order, c = cos(phi),
                         s = sin(phi), w_re = 1.0 - rho * c, w_im = rho * s;
        const CeedScalar l_re = 0.5 * log(w_re * w_re + w_im * w_im),
                         l_im = atan2(w_im, w_re);
        L_re += c * l_re - s * l_im;
        L_im += c * l_im + s * l_re;
      }
      q_re = z * L_re / (order * alpha_max);
      q_im = z * L_im / (order * alpha_max);
    }
    q_re *= sigma_max * t;
    q_im *= sigma_max * t;
    CeedScalar s_re, s_im;
    PMLStretchFactor(ctx, 1.0, sigma_max, 1.0, alpha_max, &s_re, &s_im);
    q_re += beyond * (s_re - 1.0);
    q_im += beyond * s_im;
  }
  *re += q_re;
  *im += q_im;
}

// Unit vector d from the closest point p of x on the core box of the extruded stretch to x,
// the distance |x - p|, and ∇d (3 x 3, column-major, (∇d)_ij = ∂d_i / ∂x_j), with the
// diagonal Jacobian dp = ∂p/∂x.
CEED_QFUNCTION_HELPER void PMLRay(const CeedIntScalar *ctx, const CeedScalar x[3],
                                  CeedScalar d[3], CeedScalar *dist, CeedScalar grad_d[9])
{
  const CeedIntScalar *lo = ctx + PALACE_PML_HEADER_SIZE + 26, *hi = lo + 3;
  CeedScalar p[3], dp[3];
  for (CeedInt c = 0; c < 3; c++)
  {
    const bool inside = (x[c] > lo[c].second && x[c] < hi[c].second);
    p[c] = inside ? x[c] : (x[c] <= lo[c].second) ? lo[c].second : hi[c].second;
    dp[c] = inside ? 1.0 : 0.0;
    d[c] = x[c] - p[c];
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
      grad_d[i + 3 * j] = ((i == j) * (1.0 - dp[i]) - d[i] * d[j]) / *dist;
    }
  }
}

// Depth ξ of x into the extruded layer and its gradient, the unit direction d of the
// extrusion at x and ∇d (3 x 3, column-major), and the distance |x - p| of x from the core
// (1 for normal extrusion, where ∇d = 0).
CEED_QFUNCTION_HELPER CeedScalar PMLDepth(const CeedIntScalar *ctx, const CeedScalar x[3],
                                          CeedScalar grad_xi[3], CeedScalar d[3],
                                          CeedScalar *dist, CeedScalar grad_d[9])
{
  const CeedIntScalar *q = ctx + PALACE_PML_HEADER_SIZE, *plane = q + 37;
  const CeedInt form = q[25].first;
  CeedScalar xi;
  if (form == PALACE_PML_NORMAL)
  {
    xi = -plane[3].second;
    for (CeedInt c = 0; c < 3; c++)
    {
      d[c] = grad_xi[c] = plane[c].second;
      xi += d[c] * x[c];
    }
    for (CeedInt k = 0; k < 9; k++)
    {
      grad_d[k] = 0.0;
    }
    *dist = 1.0;
  }
  else
  {
    PMLRay(ctx, x, d, dist, grad_d);
    for (CeedInt c = 0; c < 3; c++)
    {
      grad_xi[c] = d[c];
    }
    xi = *dist - q[36].second;
  }
  return (xi > 0.0) ? xi : 0.0;
}

// Jacobian J̃ (complex 3 x 3, column-major) of the extruded stretch at x, or the one with
// the magnitudes of the stretch factors.
CEED_QFUNCTION_HELPER void PMLExtrudedJacobian(const CeedIntScalar *ctx,
                                               const CeedScalar x[3], bool magnitude,
                                               CeedScalar J_re[9], CeedScalar J_im[9])
{
  const CeedIntScalar *p = ctx + PALACE_PML_HEADER_SIZE;
  const CeedScalar t = p[32].second;
  CeedScalar s_re, s_im, D_re, D_im, grad_xi[3], d[3], dist, grad_d[9];
  const CeedScalar xi = PMLDepth(ctx, x, grad_xi, d, &dist, grad_d);
  PMLStretchFactor(ctx, xi / t, p[33].second, p[34].second, p[35].second, &s_re, &s_im);
  PMLStretchIntegral(ctx, xi, &D_re, &D_im);
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

// Requested parts of c μ̃⁻¹ and/or c ε̃ (either may be null) at a quadrature point with the
// coordinates x, from the Jacobian of the Cartesian or the extruded stretch. Returns false
// (and leaves the coefficients unset) for attributes outside of the PML regions.
CEED_QFUNCTION_HELPER bool PMLCoeffAt(const CeedIntScalar *ctx, CeedInt attr,
                                      const CeedScalar x[3], CeedScalar mu_coeff[9],
                                      CeedScalar eps_coeff[9])
{
  const CeedIntScalar *b = PMLBackgroundData(ctx, attr);
  if (!b)
  {
    return false;
  }
  const CeedInt part = PMLPart(ctx);
  const bool magnitude = (part == PALACE_PML_PART_MAG_RE || part == PALACE_PML_PART_MAG_IM);
  CeedScalar J_re[9], J_im[9];
  if (PMLIsExtruded(ctx))
  {
    PMLExtrudedJacobian(ctx, x, magnitude, J_re, J_im);
  }
  else
  {
    PMLCartesianJacobian(ctx, x, magnitude, J_re, J_im);
  }
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
