// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <algorithm>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "models/materialoperator.hpp"

#include "fem/qfunctions/types_qf.h"

namespace palace::ceed
{

namespace internal
{

inline void MakeDiagonalCoefficient(MatCoeffContext1 &ctx, CeedInt k = 0,
                                    CeedScalar a = 1.0)
{
  ctx.mat_coeff[k][0] = a;
}

inline void MakeDiagonalCoefficient(MatCoeffContext2 &ctx, CeedInt k = 0,
                                    CeedScalar a = 1.0)
{
  ctx.mat_coeff[k][0] = a;
  ctx.mat_coeff[k][1] = 0.0;
  ctx.mat_coeff[k][2] = a;
}

inline void MakeDiagonalCoefficient(MatCoeffContext3 &ctx, CeedInt k = 0,
                                    CeedScalar a = 1.0)
{
  ctx.mat_coeff[k][0] = a;
  ctx.mat_coeff[k][1] = 0.0;
  ctx.mat_coeff[k][2] = 0.0;
  ctx.mat_coeff[k][3] = a;
  ctx.mat_coeff[k][4] = 0.0;
  ctx.mat_coeff[k][5] = a;
}

template <int NUM_COEFF_COMP>
inline MatCoeffContextN<NUM_COEFF_COMP>
PopulateCoefficientContext(const MaterialPropertyCoefficient *Q)
{
  MatCoeffContextN<NUM_COEFF_COMP> ctx = {{0}};
  if (!Q)
  {
    // All attributes map to identity coefficient.
    MakeDiagonalCoefficient(ctx);
  }
  else
  {
    const auto &attr_mat = Q->GetAttributeToMaterial();
    const auto &mat_coeff = Q->GetMaterialProperties();
    MFEM_VERIFY(attr_mat.Size() > 0 && attr_mat.Size() <= ctx.MaxAttr(),
                "Overflow in number of attributes for MaterialPropertyCoefficient ("
                    << attr_mat.Size() << " vs. " << ctx.MaxAttr() << ")!");
    MFEM_VERIFY(
        mat_coeff.SizeK() > 0 && mat_coeff.SizeK() < ctx.MaxNumMat(),
        "Overflow in number of material properties for MaterialPropertyCoefficient ("
            << mat_coeff.SizeK() << " vs. " << ctx.MaxNumMat() - 1 << ")!");

    // Map unassigned attributes to zero material property coefficient (the last material
    // property is reserved for zero).
    MFEM_VERIFY(attr_mat.Max() < mat_coeff.SizeK(),
                "Invalid attribute material property for MaterialPropertyCoefficient ("
                    << attr_mat.Max() << " vs. " << mat_coeff.SizeK() << ")!");
    const int zero_mat = mat_coeff.SizeK();
    std::fill(ctx.attr_mat, ctx.attr_mat + ctx.MaxAttr(), zero_mat);
    std::transform(attr_mat.begin(), attr_mat.end(), ctx.attr_mat,
                   [zero_mat](CeedInt v) { return (v < 0) ? zero_mat : v; });

    // Copy material properties: Matrix-valued material properties are always assumed to be
    // symmetric and we store only the lower triangular part.
    MFEM_VERIFY(
        mat_coeff.SizeI() == mat_coeff.SizeJ() &&
            (mat_coeff.SizeI() == 1 ||
             mat_coeff.SizeI() * (mat_coeff.SizeI() + 1) / 2 == NUM_COEFF_COMP),
        "Dimension mismatch for MaterialPropertyCoefficient and ceed::MatCoeffContext!");
    const int dim = mat_coeff.SizeI();
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      if (dim == 1)
      {
        // Copy as diagonal matrix coefficient.
        MakeDiagonalCoefficient(ctx, k, mat_coeff(0, 0, k));
      }
      else
      {
        for (int dj = 0; dj < dim; ++dj)
        {
          for (int di = dj; di < dim; ++di)
          {
            const int idx = (dj * dim) - (((dj - 1) * dj) / 2) + di - dj;
            ctx.mat_coeff[k][idx] = mat_coeff(di, dj, k);  // Column-major
          }
        }
      }
    }
  }
  return ctx;
}

}  // namespace internal

inline MatCoeffContext1
PopulateCoefficientContext1(const MaterialPropertyCoefficient *Q = nullptr)
{
  // Scalar coefficients.
  return internal::PopulateCoefficientContext<1>(Q);
}

inline MatCoeffContext2
PopulateCoefficientContext2(const MaterialPropertyCoefficient *Q = nullptr)
{
  // 2D matrix-valued coefficients.
  return internal::PopulateCoefficientContext<3>(Q);
}

inline MatCoeffContext3
PopulateCoefficientContext3(const MaterialPropertyCoefficient *Q = nullptr)
{
  // 3D matrix-valued coefficients.
  return internal::PopulateCoefficientContext<6>(Q);
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
