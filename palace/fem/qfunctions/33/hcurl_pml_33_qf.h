// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// QFunctions building the quadrature data of the PML terms in 3D, with the material tensors
// evaluated at each quadrature point from the physical coordinate (see coeff/pml_qf.h for
// the context layout and the tensor definitions), for the generic f_apply_33 and
// f_apply_33_cross QFunctions. Each QFunction computes one part (real, imaginary, or
// magnitude) of a complex prefactor times the PML tensor:
//
//   curl          : c μ̃⁻¹ (curl u, curl v), the curl lives in H(div) and is pulled back
//                   with J / |J| (mirrors f_apply_hdiv_33).
//   mass          : c ε̃ (u, v) for u in H(curl), pulled back with adj(J)ᵀ / |J| (mirrors
//                   f_apply_hcurl_33). With the gradient evaluation mode this is also the
//                   H1 diffusion form c ε̃ (∇u, ∇v), since gradients pull back the same way.
//   curlmass      : the sum of the two above, with separate prefactors.
//   floquet_mass  : c [k ×]ᵀ μ̃⁻¹ [k ×] (u, v) (or the corresponding H1 diffusion form).
//   floquet_cross : c ([k ×]ᵀ μ̃⁻¹ curl u, v) - c (μ̃⁻¹ [k ×] u, curl v).
//
// geom_data layout: {attr, w|J|, adj(J)ᵀ / |J|, x}, (2 + 9 + 3) components.

CEED_QFUNCTION_HELPER void PMLGeomUnpack33(const CeedScalar *geom, CeedInt Q, CeedInt i,
                                           CeedScalar adjJt[9], CeedScalar x[3])
{
  MatUnpack33(geom + 2 * Q + i, Q, adjJt);
  x[0] = geom[11 * Q + i];
  x[1] = geom[12 * Q + i];
  x[2] = geom[13 * Q + i];
}

// Coefficient [k ×]ᵀ (c μ̃⁻¹) [k ×] for the Floquet mass term.
CEED_QFUNCTION_HELPER bool PMLFloquetMassCoeff(const CeedIntScalar *ctx, CeedInt attr,
                                               const CeedScalar x[3], CeedScalar coeff[9])
{
  CeedScalar mu_inv[9], K[9];
  if (!PMLCoeffAt(ctx, attr, x, mu_inv, 0))
  {
    return false;
  }
  const CeedIntScalar *kx = PMLWaveVectorCross(ctx);
  for (CeedInt k = 0; k < 9; k++)
  {
    K[k] = kx[k].second;
  }
  MultAtBC33(K, mu_inv, K, coeff);
  return true;
}

template <int FORM>
CEED_QFUNCTION_HELPER int f_build_hcurl_pml_33_impl(void *__restrict__ ctx, CeedInt Q,
                                                    const CeedScalar *const *in,
                                                    CeedScalar *const *out)
{
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd = out[0];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], coeff[9],
        qd_loc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    const CeedInt attr_i = (CeedInt)attr[i];
    const bool active = (FORM == 0)   ? PMLCoeffAt(pml_ctx, attr_i, x_loc, coeff, 0)
                        : (FORM == 1) ? PMLCoeffAt(pml_ctx, attr_i, x_loc, 0, coeff)
                                      : PMLFloquetMassCoeff(pml_ctx, attr_i, x_loc, coeff);
    if (active)
    {
      if (FORM == 0)
      {
        CeedScalar J_loc[9];
        AdjJt33(adjJt_loc, J_loc);
        MultAtBA33(J_loc, coeff, qd_loc);
      }
      else
      {
        MultAtBA33(adjJt_loc, coeff, qd_loc);
      }
    }
    for (CeedInt k = 0; k < 9; k++)
    {
      qd[i + Q * k] = wdetJ[i] * qd_loc[k];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_hcurl_pml_curl_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<0>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<1>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_floquet_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                                  const CeedScalar *const *in,
                                                  CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<2>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_curlmass_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // Quadrature data for f_apply_33: the mass block (interp) followed by the curl block.
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 9 * Q;
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], mu_coeff[9], eps_coeff[9],
        qd1_loc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        qd2_loc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    if (PMLCoeffAt(pml_ctx, (CeedInt)attr[i], x_loc, mu_coeff, eps_coeff))
    {
      MultAtBA33(adjJt_loc, eps_coeff, qd1_loc);

      CeedScalar J_loc[9];
      AdjJt33(adjJt_loc, J_loc);
      MultAtBA33(J_loc, mu_coeff, qd2_loc);
    }
    for (CeedInt k = 0; k < 9; k++)
    {
      qd1[i + Q * k] = wdetJ[i] * qd1_loc[k];
      qd2[i + Q * k] = wdetJ[i] * qd2_loc[k];
    }
  }
  return 0;
}

CEED_QFUNCTION(f_build_hcurl_pml_floquet_cross_33)(void *__restrict__ ctx, CeedInt Q,
                                                   const CeedScalar *const *in,
                                                   CeedScalar *const *out)
{
  // Quadrature data for f_apply_33_cross: the curl-to-field block ([k ×]ᵀ μ̃⁻¹ curl u, v)
  // followed by the field-to-curl block -(μ̃⁻¹ [k ×] u, curl v).
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q;
  CeedScalar *qd1 = out[0], *qd2 = out[0] + 9 * Q;
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], coeff[9],
        qd1_loc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        qd2_loc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    if (PMLCoeffAt(pml_ctx, (CeedInt)attr[i], x_loc, coeff, 0))
    {
      CeedScalar J_loc[9], K[9], KA[9];
      AdjJt33(adjJt_loc, J_loc);
      const CeedIntScalar *kx = PMLWaveVectorCross(pml_ctx);
      for (CeedInt k = 0; k < 9; k++)
      {
        K[k] = kx[k].second;
      }

      // With A = adj(J)ᵀ / |J| and B = J / |J|: (K A)ᵀ μ̃⁻¹ B maps curl u to v, and
      // -Bᵀ μ̃⁻¹ (K A) maps u to curl v.
      MultBA33(adjJt_loc, K, KA);
      MultAtBC33(KA, coeff, J_loc, qd1_loc);
      MultAtBC33(J_loc, coeff, KA, qd2_loc);
    }
    for (CeedInt k = 0; k < 9; k++)
    {
      qd1[i + Q * k] = wdetJ[i] * qd1_loc[k];
      qd2[i + Q * k] = -wdetJ[i] * qd2_loc[k];
    }
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
