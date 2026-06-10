// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SURFACE_FUNCTIONAL_HPP
#define PALACE_FEM_SURFACE_FUNCTIONAL_HPP

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class Mesh;

//
// Class to compute output functionals (integrals of functions of solution fields) over
// surfaces (sets of boundary elements) of a 3D mesh using libCEED, supporting full
// (non-trace) evaluation of volume fields at boundary element quadrature points. This
// enables postprocessing measurements (interface dielectric energy participation,
// surface fluxes, port powers, etc.) to execute on the device, in contrast to the
// legacy mfem::Coefficient-based paths which are host-only.
//
// The key construction: for each boundary element, the field is evaluated from the
// attached volume element (element 1, following the conventions of
// BdrGridFunctionCoefficient). Boundary elements are grouped by the mapped positions of
// the face quadrature points in the volume element's reference space, so that each
// group shares a single tabulated basis (the volume element basis evaluated at the
// mapped face quadrature points) and element restriction (the volume element dofs,
// reusing the standard volume restriction machinery including H(curl)/H(div) dof
// orientations and transformations).
//
class SurfaceFunctional
{
public:
  enum class Kind
  {
    AREA,        // ∫ dS (no field input, for validation)
    HCURL_NORM2  // ∫ |u|² dS for an H(curl) field u
  };

private:
  // Computation kind.
  Kind kind;

  // Field finite element space (not owned, may be nullptr for field-less functionals).
  const FiniteElementSpace *fespace;

  // MPI communicator from the mesh.
  MPI_Comm comm;

  // Per-group assembled libCEED operators. Each operator integrates over one group of
  // boundary elements and accumulates per-element integrals into the local output
  // vector (CeedOperatorApplyAdd with all field inputs passive).
  struct Group
  {
    Ceed ceed;
    CeedOperator op;
  };
  std::vector<Group> groups;

  // Staging vector used to initialize the field input CeedVector at construction. The
  // field CeedVector is re-pointed at the caller's data on each Eval() call.
  mutable Vector field_staging;

  // Local output vector with one slot per marked boundary element on this process.
  mutable Vector local_out;

  void Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);

public:
  // Construct a functional over the boundary elements with marked attributes (marker
  // over global mfem boundary attributes). For field-less functionals (AREA), fespace
  // may be nullptr but the mesh is still required.
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const FiniteElementSpace *fespace = nullptr);
  ~SurfaceFunctional();

  SurfaceFunctional(const SurfaceFunctional &) = delete;
  SurfaceFunctional &operator=(const SurfaceFunctional &) = delete;

  // Evaluate the functional for the given field (L-vector, e.g. the local vector of a
  // GridFunction on the field space). Collective on the mesh communicator. For
  // field-less functionals, u is ignored and may be nullptr.
  double Eval(const Vector *u = nullptr) const;
};

}  // namespace palace

#endif  // PALACE_FEM_SURFACE_FUNCTIONAL_HPP
