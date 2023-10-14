// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_VECFEMASS_QF_H
#define PALACE_LIBCEED_VECFEMASS_QF_H

struct VectorFEMassContext
{
  CeedInt dim, space_dim;
  bool sym;
  CeedScalar coeff;
};

// libCEED QFunction for applying a symmetric or nonsymmetric vector FE mass operator.
CEED_QFUNCTION(f_apply_vecfemass)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                                  CeedScalar *const *out)
{
  // in[0], out[0] have shape [dim, ncomp=1, Q]
  VectorFEMassContext *bc = (VectorFEMassContext *)ctx;
  const CeedScalar *u = in[0], *qd = in[1];
  CeedScalar *v = out[0];
  switch (bc->dim)
  {
    case 1:
      CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
      {
        v[i] = qd[i] * u[i];
      }
      break;
    case 2:
      if (bc->sym)
      {
        CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
        {
          const CeedScalar u0 = u[i + Q * 0];
          const CeedScalar u1 = u[i + Q * 1];
          v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1;
          v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 2] * u1;
        }
      }
      else
      {
        CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
        {
          const CeedScalar u0 = u[i + Q * 0];
          const CeedScalar u1 = u[i + Q * 1];
          v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 2] * u1;
          v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 3] * u1;
        }
      }
      break;
    case 3:
      if (bc->sym)
      {
        CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
        {
          const CeedScalar u0 = u[i + Q * 0];
          const CeedScalar u1 = u[i + Q * 1];
          const CeedScalar u2 = u[i + Q * 2];
          v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 1] * u1 + qd[i + Q * 2] * u2;
          v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 4] * u2;
          v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 5] * u2;
        }
      }
      else
      {
        CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
        {
          const CeedScalar u0 = u[i + Q * 0];
          const CeedScalar u1 = u[i + Q * 1];
          const CeedScalar u2 = u[i + Q * 2];
          v[i + Q * 0] = qd[i + Q * 0] * u0 + qd[i + Q * 3] * u1 + qd[i + Q * 6] * u2;
          v[i + Q * 1] = qd[i + Q * 1] * u0 + qd[i + Q * 4] * u1 + qd[i + Q * 7] * u2;
          v[i + Q * 2] = qd[i + Q * 2] * u0 + qd[i + Q * 5] * u1 + qd[i + Q * 8] * u2;
        }
      }
      break;
  }
  return 0;
}

#endif  // PALACE_LIBCEED_VECFEMASS_QF_H
