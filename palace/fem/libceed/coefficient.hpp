// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <vector>
#include <mfem.hpp>
#include <mfem/linalg/dtensor.hpp>

namespace palace::ceed
{

struct QuadratureCoefficient
{
  int ncomp;
  mfem::Vector data;
};

inline void InitCoefficient(mfem::Coefficient &Q, mfem::Mesh &mesh,
                            const mfem::IntegrationRule &ir,
                            const std::vector<int> &indices, bool use_bdr,
                            QuadratureCoefficient &coeff)
{
  const auto ne = indices.size();
  const auto nqpts = ir.GetNPoints();
  coeff.ncomp = 1;
  coeff.data.SetSize(ne * nqpts);
  auto C = mfem::Reshape(coeff.data.HostWrite(), nqpts, ne);
  mfem::IsoparametricTransformation T;
  for (std::size_t i = 0; i < ne; ++i)
  {
    const auto e = indices[i];
    if (use_bdr)
    {
      mesh.GetBdrElementTransformation(e, &T);
    }
    else
    {
      mesh.GetElementTransformation(e, &T);
    }
    for (int q = 0; q < nqpts; ++q)
    {
      const mfem::IntegrationPoint &ip = ir.IntPoint(q);
      T.SetIntPoint(&ip);
      C(q, i) = Q.Eval(T, ip);
    }
  }
}

inline void InitCoefficient(mfem::VectorCoefficient &VQ, mfem::Mesh &mesh,
                            const mfem::IntegrationRule &ir,
                            const std::vector<int> &indices, bool use_bdr,
                            QuadratureCoefficient &coeff)
{
  const auto ne = indices.size();
  const auto vdim = VQ.GetVDim();
  const auto nqpts = ir.GetNPoints();
  coeff.ncomp = vdim;
  coeff.data.SetSize(ne * nqpts * vdim);
  auto C = mfem::Reshape(coeff.data.HostWrite(), vdim, nqpts, ne);
  mfem::IsoparametricTransformation T;
  mfem::DenseMatrix Q_ip;
  for (std::size_t i = 0; i < ne; ++i)
  {
    const auto e = indices[i];
    if (use_bdr)
    {
      mesh.GetBdrElementTransformation(e, &T);
    }
    else
    {
      mesh.GetElementTransformation(e, &T);
    }
    VQ.Eval(Q_ip, T, ir);
    for (int q = 0; q < nqpts; ++q)
    {
      for (int d = 0; d < vdim; ++d)
      {
        C(d, q, i) = Q_ip(d, q);
      }
    }
  }
}

inline void InitCoefficient(mfem::MatrixCoefficient &MQ, mfem::Mesh &mesh,
                            const mfem::IntegrationRule &ir,
                            const std::vector<int> &indices, bool use_bdr,
                            QuadratureCoefficient &coeff)
{
  // Assumes matrix coefficient is symmetric.
  const auto ne = indices.size();
  const auto vdim = MQ.GetVDim();
  const auto ncomp = (vdim * (vdim + 1)) / 2;
  const auto nqpts = ir.GetNPoints();
  coeff.ncomp = ncomp;
  coeff.data.SetSize(ne * nqpts * ncomp);
  auto C = mfem::Reshape(coeff.data.HostWrite(), ncomp, nqpts, ne);
  mfem::IsoparametricTransformation T;
  mfem::DenseMatrix Q_ip;
  for (std::size_t i = 0; i < ne; ++i)
  {
    const auto e = indices[i];
    if (use_bdr)
    {
      mesh.GetBdrElementTransformation(e, &T);
    }
    else
    {
      mesh.GetElementTransformation(e, &T);
    }
    for (int q = 0; q < nqpts; ++q)
    {
      const mfem::IntegrationPoint &ip = ir.IntPoint(q);
      T.SetIntPoint(&ip);
      MQ.Eval(Q_ip, T, ip);
      for (int dj = 0; dj < vdim; ++dj)
      {
        for (int di = dj; di < vdim; ++di)
        {
          const int idx = (dj * vdim) - (((dj - 1) * dj) / 2) + di - dj;
          C(idx, q, i) = Q_ip(di, dj);  // Column-major
        }
      }
    }
  }
}

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
