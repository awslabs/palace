// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "coefficient.hpp"

#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "models/materialoperator.hpp"

#include "fem/qfunctions/coeff_qf.h"

namespace palace::ceed
{

namespace
{

inline constexpr auto DefaultNumAttr()
{
  return 64;
}

template <int DIM>
inline constexpr auto CoeffDim()
{
  return DIM * (DIM + 1) / 2;
}

template <int DIM>
auto InitDefaultCoefficient()
{
  // All entries are value-initialized to zero, including the material property coefficient.
  std::vector<CeedIntScalar> ctx(2 + DefaultNumAttr() + CoeffDim<DIM>(), {0});
  ctx[0].first = DefaultNumAttr();
  ctx[1 + DefaultNumAttr()].first = 1;
  return ctx;
}

template <int DIM>
void MakeDiagonalCoefficient(CeedIntScalar *mat_coeff, CeedScalar a, CeedInt k);

template <>
void MakeDiagonalCoefficient<1>(CeedIntScalar *mat_coeff, CeedScalar a, CeedInt k)
{
  mat_coeff[k].second = a;
}

template <>
void MakeDiagonalCoefficient<2>(CeedIntScalar *mat_coeff, CeedScalar a, CeedInt k)
{
  mat_coeff[3 * k + 0].second = a;
  mat_coeff[3 * k + 1].second = 0.0;
  mat_coeff[3 * k + 2].second = a;
}

template <>
void MakeDiagonalCoefficient<3>(CeedIntScalar *mat_coeff, CeedScalar a, CeedInt k)
{
  mat_coeff[6 * k + 0].second = a;
  mat_coeff[6 * k + 1].second = 0.0;
  mat_coeff[6 * k + 2].second = 0.0;
  mat_coeff[6 * k + 3].second = a;
  mat_coeff[6 * k + 4].second = 0.0;
  mat_coeff[6 * k + 5].second = a;
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

template <int DIM>
std::vector<CeedIntScalar> PopulateCoefficientContext(const MaterialPropertyCoefficient *Q,
                                                      double a)
{
  if (!Q)
  {
    // All attributes map to identity coefficient.
    auto ctx = InitDefaultCoefficient<DIM>();
    MakeDiagonalCoefficient<DIM>(MatCoeff(ctx.data()), a, 0);
    return ctx;
  }

  const auto &attr_mat = Q->GetAttributeToMaterial();
  const auto &mat_coeff = Q->GetMaterialProperties();
  MFEM_VERIFY(attr_mat.Size() > 0, "Empty attributes for MaterialPropertyCoefficient!");
  MFEM_VERIFY(mat_coeff.SizeK() > 0,
              "Empty material properties for MaterialPropertyCoefficient!");
  MFEM_VERIFY(attr_mat.Max() < mat_coeff.SizeK(),
              "Invalid attribute material property for MaterialPropertyCoefficient ("
                  << attr_mat.Max() << " vs. " << mat_coeff.SizeK() << ")!");
  MFEM_VERIFY(mat_coeff.SizeI() == mat_coeff.SizeJ() &&
                  (mat_coeff.SizeI() == 1 || mat_coeff.SizeI() == DIM),
              "Dimension mismatch for MaterialPropertyCoefficient and libCEED integrator!");

  // Map unassigned attributes to zero material property coefficient (the last material
  // property is reserved for zero).
  std::vector<CeedIntScalar> ctx(2 + attr_mat.Size() +
                                 CoeffDim<DIM>() * (mat_coeff.SizeK() + 1));
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
  const int dim = mat_coeff.SizeI();
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    if (dim == 1)
    {
      // Copy as diagonal matrix coefficient.
      MakeDiagonalCoefficient<DIM>(MatCoeff(ctx.data()), a * mat_coeff(0, 0, k), k);
    }
    else
    {
      for (int dj = 0; dj < dim; ++dj)
      {
        for (int di = dj; di < dim; ++di)
        {
          const int idx = (dj * dim) - (((dj - 1) * dj) / 2) + di - dj;
          MatCoeff(ctx.data())[CoeffDim<DIM>() * k + idx].second =
              a * mat_coeff(di, dj, k);  // Column-major
        }
      }
    }
  }

  return ctx;
}

template <int DIM, int DIM_MASS>
std::vector<CeedIntScalar>
PopulateCoefficientContext(const MaterialPropertyCoefficient *Q,
                           const MaterialPropertyCoefficient *Q_mass, double a,
                           double a_mass)
{
  auto ctx = PopulateCoefficientContext<DIM>(Q, a);
  auto ctx_mass = PopulateCoefficientContext<DIM_MASS>(Q_mass, a_mass);
  ctx.insert(ctx.end(), ctx_mass.begin(), ctx_mass.end());
  return ctx;
}

template std::vector<CeedIntScalar>
PopulateCoefficientContext<1>(const MaterialPropertyCoefficient *, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<2>(const MaterialPropertyCoefficient *, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<3>(const MaterialPropertyCoefficient *, double);

template std::vector<CeedIntScalar>
PopulateCoefficientContext<2, 1>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<3, 1>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<1, 2>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<1, 3>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<2, 2>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);
template std::vector<CeedIntScalar>
PopulateCoefficientContext<3, 3>(const MaterialPropertyCoefficient *,
                                 const MaterialPropertyCoefficient *, double, double);

}  // namespace palace::ceed
