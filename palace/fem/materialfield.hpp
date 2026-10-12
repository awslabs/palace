// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MATERIALFIELD_HPP
#define PALACE_FEM_MATERIALFIELD_HPP

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace palace
{

class Mesh;

//
// Quadrature-point samples of a scalar material property over a group of elements of one
// geometry type, stored like the mesh geometry factor data so that it can be attached to a
// libCEED operator as an extra inactive (CEED_EVAL_NONE) input.
//
struct CeedMaterialFieldData
{
  CeedVector data = nullptr;
  CeedElemRestriction data_restr = nullptr;
};

//
// A spatially-varying scalar material coefficient defined by a grid function on the solve
// mesh. The grid function is sampled at the quadrature points used by libCEED operators
// (which are fixed by the mesh nodal finite element and the solver order, and so are the
// same at every level of a p-multigrid hierarchy), per element geometry. The resulting
// quadrature data is cached and shared by all operators on the mesh.
//
class MaterialField
{
public:
  // Construct from a (real) scalar grid function on the solve mesh. `scaling` multiplies
  // every sampled value (e.g. use 1.0 to supply the coefficient itself).
  MaterialField(const mfem::ParGridFunction &gf, double scaling = 1.0);
  ~MaterialField();

  MaterialField(const MaterialField &) = delete;
  MaterialField &operator=(const MaterialField &) = delete;

  // Return (building on first use) the quadrature data for the given libCEED context and
  // element geometry. `indices` must be the element ordering used by CeedGeomFactorData for
  // the same geometry. Callers should warm the cache before entering a parallel region
  // (see MaterialField::WarmCache).
  const CeedMaterialFieldData &GetData(Ceed ceed, mfem::Geometry::Type geom,
                                       const std::vector<int> &indices) const;

  // Build quadrature data for all element geometries of the mesh on the given libCEED
  // context, so that subsequent GetData calls are read-only (thread-safe).
  void WarmCache(Ceed ceed, const Mesh &mesh) const;

  // Maximum grid function coefficient over this process (reduce across ranks for a global
  // bound).
  double GetLocalMax() const { return local_max; }

private:
  const mfem::ParGridFunction *gf;
  double scaling;
  double local_max;
  mutable ceed::CeedObjectMap<CeedMaterialFieldData> data;

  CeedMaterialFieldData BuildData(Ceed ceed, mfem::Geometry::Type geom,
                                  const std::vector<int> &indices) const;
};

}  // namespace palace

#endif  // PALACE_FEM_MATERIALFIELD_HPP
