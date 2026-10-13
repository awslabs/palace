// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_PML_33_QF_H
#define PALACE_LIBCEED_HCURL_PML_33_QF_H

#include "../coeff/pml_qf.h"
#include "utils_33_qf.h"

// QFunctions for the Cartesian PML terms in 3D, with the material tensors evaluated at each
// quadrature point from the physical coordinate (see coeff/pml_qf.h for the context layout
// and the tensor definitions). Each QFunction computes one part (real, imaginary, or
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
// The build variants assemble the quadrature data for the generic f_apply_3 / f_apply_33
// QFunctions so that the PML tensors are only evaluated once per quadrature point.
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

// Compute y = Aᵀ x for a 3x3 matrix A stored column-major.
CEED_QFUNCTION_HELPER void PMLMultAtx33(const CeedScalar A[9], const CeedScalar x[3],
                                        CeedScalar y[3])
{
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

// Coefficient [k ×]ᵀ (c μ̃⁻¹) [k ×] for the Floquet mass term.
CEED_QFUNCTION_HELPER bool PMLFloquetMassCoeff(const CeedIntScalar *ctx, CeedInt attr,
                                               const CeedScalar x[3], CeedScalar coeff[9])
{
  CeedScalar mu_inv[9], K[9];
  if (!PMLMuInvCoeff(ctx, attr, x, mu_inv))
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
CEED_QFUNCTION_HELPER int f_apply_hcurl_pml_33_impl(void *__restrict__ ctx, CeedInt Q,
                                                    const CeedScalar *const *in,
                                                    CeedScalar *const *out)
{
  // FORM: 0 = curl (J), 1 = mass (adj(J)ᵀ), 2 = Floquet mass (adj(J)ᵀ).
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *u = in[1];
  CeedScalar *v = out[0];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], coeff[9], v_loc[3] = {0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    const bool active = (FORM == 0) ? PMLMuInvCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff)
                        : (FORM == 1)
                            ? PMLEpsCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff)
                            : PMLFloquetMassCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff);
    if (active)
    {
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      if (FORM == 0)
      {
        CeedScalar J_loc[9];
        AdjJt33(adjJt_loc, J_loc);
        MultAtBCx33(J_loc, coeff, J_loc, u_loc, v_loc);
      }
      else
      {
        MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);
      }
    }
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
  }
  return 0;
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
    const bool active = (FORM == 0) ? PMLMuInvCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff)
                        : (FORM == 1)
                            ? PMLEpsCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff)
                            : PMLFloquetMassCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff);
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

CEED_QFUNCTION(f_apply_hcurl_pml_curl_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<0>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_curl_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<0>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<1>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                          const CeedScalar *const *in,
                                          CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<1>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                                  const CeedScalar *const *in,
                                                  CeedScalar *const *out)
{
  return f_apply_hcurl_pml_33_impl<2>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_build_hcurl_pml_floquet_mass_33)(void *__restrict__ ctx, CeedInt Q,
                                                  const CeedScalar *const *in,
                                                  CeedScalar *const *out)
{
  return f_build_hcurl_pml_33_impl<2>(ctx, Q, in, out);
}

CEED_QFUNCTION(f_apply_hcurl_pml_curlmass_33)(void *__restrict__ ctx, CeedInt Q,
                                              const CeedScalar *const *in,
                                              CeedScalar *const *out)
{
  // Active inputs/outputs are ordered as (interp, curl).
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *u = in[1], *curl_u = in[2];
  CeedScalar *v = out[0], *curl_v = out[1];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], mu_coeff[9], eps_coeff[9],
        v_loc[3] = {0.0, 0.0, 0.0}, curl_v_loc[3] = {0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    if (PMLMuInvEpsCoeff(pml_ctx, (CeedInt)attr[i], x_loc, mu_coeff, eps_coeff))
    {
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      MultAtBCx33(adjJt_loc, eps_coeff, adjJt_loc, u_loc, v_loc);

      CeedScalar J_loc[9];
      const CeedScalar curl_u_loc[3] = {curl_u[i + Q * 0], curl_u[i + Q * 1],
                                        curl_u[i + Q * 2]};
      AdjJt33(adjJt_loc, J_loc);
      MultAtBCx33(J_loc, mu_coeff, J_loc, curl_u_loc, curl_v_loc);
    }
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
    curl_v[i + Q * 0] = wdetJ[i] * curl_v_loc[0];
    curl_v[i + Q * 1] = wdetJ[i] * curl_v_loc[1];
    curl_v[i + Q * 2] = wdetJ[i] * curl_v_loc[2];
  }
  return 0;
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
    if (PMLMuInvEpsCoeff(pml_ctx, (CeedInt)attr[i], x_loc, mu_coeff, eps_coeff))
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

CEED_QFUNCTION(f_apply_hcurl_pml_floquet_cross_33)(void *__restrict__ ctx, CeedInt Q,
                                                   const CeedScalar *const *in,
                                                   CeedScalar *const *out)
{
  // Active inputs/outputs are ordered as (interp, curl). The two halves match the
  // transposed MixedVectorCurlIntegrator and the MixedVectorWeakCurlIntegrator used for the
  // Floquet terms outside of the PML.
  const CeedScalar *attr = in[0], *wdetJ = in[0] + Q, *u = in[1], *curl_u = in[2];
  CeedScalar *v = out[0], *curl_v = out[1];
  const CeedIntScalar *pml_ctx = (const CeedIntScalar *)ctx;

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    CeedScalar adjJt_loc[9], x_loc[3], coeff[9], v_loc[3] = {0.0, 0.0, 0.0},
                                                 curl_v_loc[3] = {0.0, 0.0, 0.0};
    PMLGeomUnpack33(in[0], Q, i, adjJt_loc, x_loc);
    if (PMLMuInvCoeff(pml_ctx, (CeedInt)attr[i], x_loc, coeff))
    {
      CeedScalar J_loc[9], K[9], t1[3], t2[3], t3[3];
      AdjJt33(adjJt_loc, J_loc);
      const CeedIntScalar *kx = PMLWaveVectorCross(pml_ctx);
      for (CeedInt k = 0; k < 9; k++)
      {
        K[k] = kx[k].second;
      }

      // Curl-to-field half: ([k ×]ᵀ μ̃⁻¹ curl u, v).
      const CeedScalar curl_u_loc[3] = {curl_u[i + Q * 0], curl_u[i + Q * 1],
                                        curl_u[i + Q * 2]};
      MultAx33(J_loc, curl_u_loc, t1);
      MultAx33(coeff, t1, t2);
      PMLMultAtx33(K, t2, t3);
      PMLMultAtx33(adjJt_loc, t3, v_loc);

      // Field-to-curl half: -(μ̃⁻¹ [k ×] u, curl v).
      const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
      MultAx33(adjJt_loc, u_loc, t1);
      MultAx33(K, t1, t2);
      MultAx33(coeff, t2, t3);
      PMLMultAtx33(J_loc, t3, curl_v_loc);
      curl_v_loc[0] = -curl_v_loc[0];
      curl_v_loc[1] = -curl_v_loc[1];
      curl_v_loc[2] = -curl_v_loc[2];
    }
    v[i + Q * 0] = wdetJ[i] * v_loc[0];
    v[i + Q * 1] = wdetJ[i] * v_loc[1];
    v[i + Q * 2] = wdetJ[i] * v_loc[2];
    curl_v[i + Q * 0] = wdetJ[i] * curl_v_loc[0];
    curl_v[i + Q * 1] = wdetJ[i] * curl_v_loc[1];
    curl_v[i + Q * 2] = wdetJ[i] * curl_v_loc[2];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_PML_33_QF_H
