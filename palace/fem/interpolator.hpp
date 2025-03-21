// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTERPOLATOR_HPP
#define PALACE_FEM_INTERPOLATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class GridFunction;
class IoData;
class FiniteElementSpace;

//
// A class which wraps MFEM's GSLIB interface for high-order field interpolation.
//
class InterpolationOperator
{
private:
#if defined(MFEM_USE_GSLIB)
  mfem::FindPointsGSLIB op;
#endif
  std::vector<int> op_idx;

  int v_dim_fes;  // dimension of interpolated vector from NDSpace

  std::vector<double> ProbeField(const mfem::ParGridFunction &U);

public:
  InterpolationOperator(const IoData &iodata, FiniteElementSpace &nd_space);

  auto GetVDim() const { return v_dim_fes; }
  const auto &GetProbes() const { return op_idx; }

  std::vector<std::complex<double>> ProbeField(const GridFunction &U);
};

namespace fem
{

// Interpolate a function on a serial or parallel mesh to a different mesh, using GSLIB.
// Similar to MFEM's field-interp miniapp.
void InterpolateFunction(const mfem::GridFunction &U, mfem::GridFunction &V);

// Interpolate a function at a specific list of points, specified using the provided
// ordering. The output vector values are always arranged byVDIM.
void InterpolateFunction(const mfem::Vector &xyz, const mfem::GridFunction &U,
                         mfem::Vector &V,
                         mfem::Ordering::Type ordering = mfem::Ordering::byNODES);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_INTERPOLATOR_HPP
