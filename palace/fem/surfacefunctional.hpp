// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SURFACE_FUNCTIONAL_HPP
#define PALACE_FEM_SURFACE_FUNCTIONAL_HPP

#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FiniteElementSpace;
class GridFunction;
class MaterialOperator;
class Mesh;

//
// Class to compute output functionals (integrals of functions of solution fields) over
// surfaces (sets of boundary elements) of a 3D mesh using libCEED, supporting full
// (non-trace) evaluation of volume fields at boundary element quadrature points. This
// enables postprocessing measurements (interface dielectric energy participation,
// surface fluxes, port powers, etc.) to execute on the device, in contrast to the
// legacy mfem::Coefficient-based paths which are host-only.
//
// The key construction: for each boundary element, the field is evaluated from an
// attached volume element (or both, with averaging, for interior boundaries, following
// the conventions of BdrGridFunctionCoefficient and its derived legacy coefficients).
// Boundary elements are grouped by the mapped positions of the face quadrature points
// in the volume element reference space(s), so that each group shares tabulated bases
// (the volume element basis evaluated at the mapped face quadrature points) and element
// restrictions (the volume element dofs, reusing the standard volume restriction
// machinery including H(curl)/H(div) dof orientations and transformations).
//
class SurfaceFunctional
{
public:
  enum class Kind
  {
    AREA,          // ∫ dS (no field input, for validation)
    HCURL_NORM2,   // ∫ |u|² dS for an H(curl) field u (single-sided, for validation)
    INTERFACE_EPR  // Interface dielectric energy following InterfaceDielectricCoefficient
  };

private:
  // Computation kind and interface dielectric parameters (INTERFACE_EPR only).
  Kind kind;
  InterfaceDielectric epr_type = InterfaceDielectric::DEFAULT;
  double epr_t = 0.0, epr_epsilon = 0.0;

  // Field finite element space (not owned, may be nullptr for field-less functionals)
  // and material operator (not owned, required for INTERFACE_EPR).
  const FiniteElementSpace *fespace;
  const MaterialOperator *mat_op;

  // MPI communicator from the mesh.
  MPI_Comm comm;

  // Per-group assembled libCEED operators. Each operator integrates over one group of
  // boundary elements and accumulates per-element integrals into the local output
  // vector (CeedOperatorApplyAdd with all field inputs passive). Groups may have one
  // (u_1) or two (u_1, u_2) field inputs which are re-pointed at the caller's data on
  // each evaluation.
  struct GroupOp
  {
    Ceed ceed;
    CeedOperator op;
    int num_fields;
  };
  std::vector<GroupOp> groups;

  // Staging vector used to initialize the field input CeedVectors at construction. The
  // field CeedVectors are re-pointed at the caller's data on each Eval() call.
  mutable Vector field_staging;

  // Local output vector with one slot per marked boundary element on this process.
  mutable Vector local_out;

  void Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);

  // Apply all group operators with the field inputs pointed at the given field vector,
  // accumulating into the local output vector.
  void ApplyAdd(const Vector *u) const;

public:
  // Construct a functional over the boundary elements with marked attributes (marker
  // over global mfem boundary attributes). For field-less functionals (AREA), fespace
  // may be nullptr but the mesh is still required.
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const FiniteElementSpace *fespace = nullptr);

  // Construct an interface dielectric energy participation functional with the given
  // interface type, thickness, and permittivity (see InterfaceDielectricCoefficient).
  SurfaceFunctional(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const FiniteElementSpace &nd_fespace, const MaterialOperator &mat_op,
                    InterfaceDielectric type, double t_i, double epsilon_i);

  ~SurfaceFunctional();

  SurfaceFunctional(const SurfaceFunctional &) = delete;
  SurfaceFunctional &operator=(const SurfaceFunctional &) = delete;

  // Evaluate the functional for the given field (L-vector, e.g. the local vector of a
  // GridFunction on the field space). Collective on the mesh communicator. For
  // field-less functionals, u is ignored and may be nullptr.
  double Eval(const Vector *u = nullptr) const;

  // Evaluate the functional for the given (possibly complex-valued) grid function. For
  // complex fields, the real and imaginary part contributions add (all implemented
  // integrands are quadratic in the field). Collective on the mesh communicator.
  double Eval(const GridFunction &u) const;
};

}  // namespace palace

#endif  // PALACE_FEM_SURFACE_FUNCTIONAL_HPP
