// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTERPOLATOR_HPP
#define PALACE_FEM_INTERPOLATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace mfem
{

// Forward declaration so the GSLIB point-locator interpolation helpers below can be
// declared unconditionally (matching the other functions in this header, which gate only
// their definitions on MFEM_USE_GSLIB). The type is fully defined by <mfem.hpp> when GSLIB
// is enabled; without it, only this declaration exists and the functions abort at runtime.
class FindPointsGSLIB;

}  // namespace mfem

namespace palace
{

class CeedProbeEvaluator;
class GridFunction;
class IoData;
class Units;
class FiniteElementSpace;

//
// A class which wraps MFEM's GSLIB interface for high-order field interpolation.
//
class InterpolationOperator
{
private:
#if defined(MFEM_USE_GSLIB)
  mfem::FindPointsGSLIB op;

  // libCEED evaluators at the located probe points, avoiding per-sample GSLIB
  // interpolation), constructed lazily per source finite element space.
  mutable std::map<const mfem::FiniteElementSpace *, std::unique_ptr<CeedProbeEvaluator>>
      ceed_probes;
#endif
  std::vector<int> op_idx;

  int v_dim_fes;  // dimension of interpolated vector from NDSpace

  std::vector<double> ProbeField(const mfem::ParGridFunction &U);

public:
  InterpolationOperator(const std::map<int, config::ProbeData> &probe, const Units &units,
                        FiniteElementSpace &nd_space);
  InterpolationOperator(const IoData &iodata, FiniteElementSpace &nd_space);
  ~InterpolationOperator();

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

// Compute the line integral V = ∫ F · dl along a straight line from p1 to p2 using GSLIB
// interpolation. The field F is a real-valued vector ParGridFunction. Uses Gauss-Legendre
// quadrature of the specified order on [0,1]. Returns the scalar integral value.
double ComputeLineIntegral(const mfem::Vector &p1, const mfem::Vector &p2,
                           const mfem::ParGridFunction &field, int quad_order);

// Setup a reusable GSLIB point locator on a fixed mesh. The geometric Setup is the
// expensive step (builds a spatial hash over all elements) and depends only on the mesh,
// so callers that interpolate repeatedly on it (e.g. wave-port voltage-path line integrals
// over a frequency sweep) should Setup once and reuse the operator via the overloads below.
// Requires MFEM_USE_GSLIB.
void SetupInterpolator(mfem::FindPointsGSLIB &op, mfem::Mesh &mesh);

// As above, but using a pre-Setup point locator (op.Setup already called on U's mesh).
void InterpolateFunction(mfem::FindPointsGSLIB &op, const mfem::Vector &xyz,
                         const mfem::GridFunction &U, mfem::Vector &V,
                         mfem::Ordering::Type ordering = mfem::Ordering::byNODES);

// As ComputeLineIntegral, but reusing a pre-Setup point locator.
double ComputeLineIntegral(mfem::FindPointsGSLIB &op, const mfem::Vector &p1,
                           const mfem::Vector &p2, const mfem::ParGridFunction &field,
                           int quad_order);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_INTERPOLATOR_HPP
