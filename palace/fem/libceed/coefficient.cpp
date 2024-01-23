// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "coefficient.hpp"

#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "models/materialoperator.hpp"

#include "fem/qfunctions/coeff/coeff_qf.h"

namespace palace::ceed
{

namespace
{

inline constexpr auto DefaultNumAttr()
{
  return 64;
}

inline auto CoeffDim(int dim)
{
  return dim * (dim + 1) / 2;
}

auto InitDefaultCoefficient(int dim)
{
  // All entries are value-initialized to zero, including the material property coefficient.
  std::vector<CeedIntScalar> ctx(2 + DefaultNumAttr() + CoeffDim(dim), {0});
  ctx[0].first = DefaultNumAttr();
  ctx[1 + DefaultNumAttr()].first = 1;
  return ctx;
}

void MakeDiagonalCoefficient(int dim, CeedIntScalar *mat_coeff, CeedScalar a, CeedInt k)
{
  const int coeff_dim = CoeffDim(dim);
  for (int i = 0; i < coeff_dim; i++)
  {
    mat_coeff[coeff_dim * k + i].second = 0.0;
  }
  for (int di = 0; di < dim; ++di)
  {
    const int idx = (di * dim) - (((di - 1) * di) / 2);
    mat_coeff[coeff_dim * k + idx].second = a;
  }
}

inline auto *AttrMat(CeedIntScalar *ctx)
{
  return ctx + 1;
}

inline auto *MatCoeff(CeedIntScalar *ctx)
{
  const CeedInt num_attr = ctx[0].first;
  return ctx + 2 + num_attr;
}

}  // namespace

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim, const MaterialPropertyCoefficient *Q, double a)
{
  if (!Q)
  {
    // All attributes map to identity coefficient.
    auto ctx = InitDefaultCoefficient(dim);
    MakeDiagonalCoefficient(dim, MatCoeff(ctx.data()), a, 0);
    return ctx;
  }

  // Material property coefficients might be empty if all attributes map to zero
  // coefficient.
  const auto &attr_mat = Q->GetAttributeToMaterial();
  const auto &mat_coeff = Q->GetMaterialProperties();
  MFEM_VERIFY(attr_mat.Size() > 0, "Empty attributes for MaterialPropertyCoefficient!");
  MFEM_VERIFY(attr_mat.Max() < mat_coeff.SizeK(),
              "Invalid attribute material property for MaterialPropertyCoefficient ("
                  << attr_mat.Max() << " vs. " << mat_coeff.SizeK() << ")!");
  MFEM_VERIFY(mat_coeff.SizeI() == mat_coeff.SizeJ() &&
                  (mat_coeff.SizeI() == 1 || mat_coeff.SizeI() == dim),
              "Dimension mismatch for MaterialPropertyCoefficient and libCEED integrator!");

  // Map unassigned attributes to zero material property coefficient (the last material
  // property is reserved for zero).
  const int coeff_dim = CoeffDim(dim);
  std::vector<CeedIntScalar> ctx(2 + attr_mat.Size() + coeff_dim * (mat_coeff.SizeK() + 1));
  ctx[0].first = attr_mat.Size();
  const int zero_mat = mat_coeff.SizeK();
  for (int i = 0; i < attr_mat.Size(); i++)
  {
    const int k = attr_mat[i];
    AttrMat(ctx.data())[i].first = (k < 0) ? zero_mat : k;
  }

  // Copy material properties: Matrix-valued material properties are always assumed to be
  // symmetric and we store only the lower triangular part.
  ctx[1 + attr_mat.Size()].first = mat_coeff.SizeK() + 1;
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    if (mat_coeff.SizeI() == 1)
    {
      // Copy as diagonal matrix coefficient.
      MakeDiagonalCoefficient(dim, MatCoeff(ctx.data()), a * mat_coeff(0, 0, k), k);
    }
    else
    {
      for (int dj = 0; dj < dim; ++dj)
      {
        for (int di = dj; di < dim; ++di)
        {
          // Column-major ordering.
          const int idx = (dj * dim) - (((dj - 1) * dj) / 2) + di - dj;
          MatCoeff(ctx.data())[coeff_dim * k + idx].second = a * mat_coeff(di, dj, k);
        }
      }
    }
  }
  for (int d = 0; d < coeff_dim; d++)
  {
    MatCoeff(ctx.data())[coeff_dim * zero_mat + d].second = 0.0;
  }

  return ctx;
}

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim, int dim_mass, const MaterialPropertyCoefficient *Q,
                           const MaterialPropertyCoefficient *Q_mass, double a,
                           double a_mass)
{
  // Mass coefficient comes first, then the other one for the QFunction.
  auto ctx = PopulateCoefficientContext(dim, Q, a);
  auto ctx_mass = PopulateCoefficientContext(dim_mass, Q_mass, a_mass);
  ctx_mass.insert(ctx_mass.end(), ctx.begin(), ctx.end());
  return ctx_mass;
}

}  // namespace palace::ceed
