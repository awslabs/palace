// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_PML_QF_H
#define PALACE_LIBCEED_PML_QF_H

#include "coeff_qf.h"

// Device-callable PML stretch evaluator. Evaluates the diagonal Cartesian UPML stretch
// tensors ε̃ and μ̃⁻¹ at a single quadrature point, given the physical coordinate x and
// a packed profile context.
//
// Context layout: the first entry is the number of libCEED attributes; the next
// `num_attr` entries are the attr→profile map (-1 for non-PML attributes); the
// remaining entries hold the per-profile data in blocks of PALACE_PML_REGION_STRIDE
// entries each.
//
// Per-region layout (29 entries, indexed 0..28):
//   [0]         formulation (int; 0=FIXED, 1=CFS, 2=FREQUENCY_DEPENDENT)
//   [1]         polynomial grading order (int)
//   [2]         mu_r
//   [3]         eps_r
//   [4]         omega (reference frequency for FIXED/CFS, live ω for FD)
//   [5..7]      sigma_max[axis] for axis ∈ {x, y, z}
//   [8..10]     kappa_max[axis]
//   [11..13]    alpha_max[axis]
//   [14..16]    thickness[axis] (per-axis max of ±face thickness)
//   [17..22]    interface_coord interleaved: {neg_x, pos_x, neg_y, pos_y, neg_z, pos_z}
//   [23..28]    direction_active interleaved: {neg_x, pos_x, neg_y, pos_y, neg_z, pos_z}
//               (int; nonzero ⇒ face active)
//
// PML formulation tags must match PMLStretchFormulation enum in labels.hpp.

#define PALACE_PML_REGION_STRIDE 29

CEED_QFUNCTION_HELPER CeedInt PMLNumProfiles(const CeedIntScalar *ctx)
{
  return ctx[0].first;
}

CEED_QFUNCTION_HELPER CeedInt PMLAttrToProfile(const CeedIntScalar *ctx, CeedInt attr)
{
  // attr is 1-based; ctx[1 + attr - 1] holds the profile index (-1 for non-PML).
  return ctx[1 + attr - 1].first;
}

CEED_QFUNCTION_HELPER const CeedIntScalar *PMLRegion(const CeedIntScalar *ctx,
                                                     CeedInt num_attr, CeedInt pidx)
{
  return ctx + 1 + num_attr + PALACE_PML_REGION_STRIDE * pidx;
}

// Integer pow for non-negative exponent.
CEED_QFUNCTION_HELPER CeedScalar PMLIntPow(CeedScalar x, CeedInt n)
{
  CeedScalar r = 1.0;
  for (CeedInt i = 0; i < n; i++)
  {
    r *= x;
  }
  return r;
}

// Compute the diagonal Cartesian PML stretch tensors μ̃⁻¹ and ε̃ at physical point x,
// for a single PML profile. Writes 12 scalars: [mu_inv_re[3], mu_inv_im[3], eps_re[3],
// eps_im[3]].
//
// Palace's e^{+iωt} convention: s_axis = 1 − iσ/ω (FIXED/FREQUENCY_DEPENDENT) or
// κ + σ/(α + iω) (CFS). See palace/models/pml.cpp for the host-side derivation.
CEED_QFUNCTION_HELPER void PMLEvalStretchTensors(const CeedIntScalar *region,
                                                 const CeedScalar x[3],
                                                 CeedScalar mu_inv_re[3],
                                                 CeedScalar mu_inv_im[3],
                                                 CeedScalar eps_re[3], CeedScalar eps_im[3])
{
  // Unpack region fields. Layout (27 entries):
  //   [0]     formulation (int)
  //   [1]     order (int)
  //   [2]     mu_r
  //   [3]     eps_r
  //   [4]     omega
  //   [5..7]  sigma_max[axis]
  //   [8..10] kappa_max[axis]
  //   [11..13] alpha_max[axis]
  //   [14..16] thickness[axis]
  //   [17, 19, 21] iface_neg[axis]     // per axis: neg, pos interleaved
  //   [18, 20, 22] iface_pos[axis]
  //   [23, 25, 27-2] active_neg[axis]
  //   [24, 26, 27-1] active_pos[axis]
  const CeedInt formulation = region[0].first;
  const CeedInt order = region[1].first;
  const CeedScalar mu_r = region[2].second;
  const CeedScalar eps_r = region[3].second;
  const CeedScalar omega = region[4].second;
  const CeedScalar sigma_max[3] = {region[5].second, region[6].second, region[7].second};
  const CeedScalar kappa_max[3] = {region[8].second, region[9].second, region[10].second};
  const CeedScalar alpha_max[3] = {region[11].second, region[12].second, region[13].second};
  const CeedScalar thickness[3] = {region[14].second, region[15].second, region[16].second};
  const CeedScalar iface_neg[3] = {region[17].second, region[19].second, region[21].second};
  const CeedScalar iface_pos[3] = {region[18].second, region[20].second, region[22].second};

  // Compute per-axis depth, then σ, κ, α via polynomial grading.
  CeedScalar sigma_loc[3] = {0.0, 0.0, 0.0};
  CeedScalar kappa_loc[3] = {1.0, 1.0, 1.0};
  CeedScalar alpha_loc[3] = {0.0, 0.0, 0.0};
  for (CeedInt axis = 0; axis < 3; axis++)
  {
    const CeedInt a_neg = region[23 + 2 * axis].first;
    const CeedInt a_pos = region[24 + 2 * axis].first;
    const CeedScalar d = thickness[axis];
    if (d <= 0.0)
    {
      continue;
    }
    CeedScalar depth = 0.0;
    if (a_neg != 0 && x[axis] < iface_neg[axis])
    {
      depth = iface_neg[axis] - x[axis];
    }
    else if (a_pos != 0 && x[axis] > iface_pos[axis])
    {
      depth = x[axis] - iface_pos[axis];
    }
    if (depth <= 0.0)
    {
      continue;
    }
    const CeedScalar t = depth / d;
    const CeedScalar shape = PMLIntPow(t, order);
    sigma_loc[axis] = sigma_max[axis] * shape;
    alpha_loc[axis] = alpha_max[axis] * shape;
    kappa_loc[axis] = 1.0 + (kappa_max[axis] - 1.0) * shape;
  }

  // Compute complex stretch s_axis per axis. Store as (re, im) pairs.
  CeedScalar s_re[3] = {1.0, 1.0, 1.0};
  CeedScalar s_im[3] = {0.0, 0.0, 0.0};
  const CeedScalar safe_omega = (omega > 0.0 || omega < 0.0) ? omega : 1.0e300;
  for (CeedInt axis = 0; axis < 3; axis++)
  {
    const CeedScalar sg = sigma_loc[axis];
    const CeedScalar kp = kappa_loc[axis];
    const CeedScalar al = alpha_loc[axis];
    if (formulation == 1)  // CFS
    {
      const CeedScalar denom = al * al + safe_omega * safe_omega;
      if (denom > 0.0)
      {
        s_re[axis] = kp + sg * al / denom;
        s_im[axis] = -sg * safe_omega / denom;
      }
      else
      {
        s_re[axis] = kp;
        s_im[axis] = 0.0;
      }
    }
    else  // FIXED or FREQUENCY_DEPENDENT
    {
      s_re[axis] = 1.0;
      s_im[axis] = -sg / safe_omega;
    }
  }

  // UPML tensors: μ̃ = μ_r · Λ, ε̃ = ε_r · Λ, Λ_ii = s_j s_k / s_i.
  // We store μ̃⁻¹, so μ̃⁻¹_ii = (1/μ_r) / Λ_ii = (1/μ_r) · s_i / (s_j · s_k).
  for (CeedInt i = 0; i < 3; i++)
  {
    const CeedInt j = (i + 1) % 3;
    const CeedInt k = (i + 2) % 3;
    // eps = eps_r * (s[j] * s[k]) / s[i]
    // num = s[j] * s[k]
    const CeedScalar num_re = s_re[j] * s_re[k] - s_im[j] * s_im[k];
    const CeedScalar num_im = s_re[j] * s_im[k] + s_im[j] * s_re[k];
    // eps = eps_r * num / s[i]
    const CeedScalar si_mag2 = s_re[i] * s_re[i] + s_im[i] * s_im[i];
    const CeedScalar eps_axis_re = eps_r * (num_re * s_re[i] + num_im * s_im[i]) / si_mag2;
    const CeedScalar eps_axis_im = eps_r * (num_im * s_re[i] - num_re * s_im[i]) / si_mag2;
    eps_re[i] = eps_axis_re;
    eps_im[i] = eps_axis_im;
    // mu_inv = (1/mu_r) * s[i] / (s[j] * s[k])
    const CeedScalar num2_re = s_re[j] * s_re[k] - s_im[j] * s_im[k];
    const CeedScalar num2_im = s_re[j] * s_im[k] + s_im[j] * s_re[k];
    const CeedScalar num2_mag2 = num2_re * num2_re + num2_im * num2_im;
    const CeedScalar mu_inv_axis_re =
        (s_re[i] * num2_re + s_im[i] * num2_im) / (mu_r * num2_mag2);
    const CeedScalar mu_inv_axis_im =
        (s_im[i] * num2_re - s_re[i] * num2_im) / (mu_r * num2_mag2);
    mu_inv_re[i] = mu_inv_axis_re;
    mu_inv_im[i] = mu_inv_axis_im;
  }
}

#endif  // PALACE_LIBCEED_PML_QF_H
