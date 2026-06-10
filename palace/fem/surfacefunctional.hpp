// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SURFACE_FUNCTIONAL_HPP
#define PALACE_FEM_SURFACE_FUNCTIONAL_HPP

#include <array>
#include <complex>
#include <string>
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
// attached volume element (or both, with averaging or differencing, for interior
// boundaries, following the conventions of BdrGridFunctionCoefficient and its derived
// legacy coefficients). Boundary elements are grouped by the mapped positions of the
// face quadrature points in the volume element reference space(s), so that each group
// shares tabulated bases (the volume element basis evaluated at the mapped face
// quadrature points) and element restrictions (the volume element dofs, reusing the
// standard volume restriction machinery including H(curl)/H(div) dof orientations and
// transformations).
//
class SurfaceFunctional
{
public:
  enum class Kind
  {
    AREA,           // ∫ dS (no field input, for validation)
    HCURL_NORM2,    // ∫ |u|² dS for an H(curl) field u (single-sided, for validation)
    INTERFACE_EPR,  // Interface dielectric energy following InterfaceDielectricCoefficient
    SURFACE_FLUX    // Surface flux following BdrSurfaceFluxCoefficient
  };

private:
  // Computation kind and integrand parameters.
  Kind kind;
  InterfaceDielectric epr_type = InterfaceDielectric::DEFAULT;
  double epr_t = 0.0, epr_epsilon = 0.0;
  SurfaceFlux flux_type = SurfaceFlux::ELECTRIC;
  bool flux_two_sided = false;
  mfem::Vector flux_x0;

  // Field finite element spaces (not owned): fespace_e for H(curl) fields (source index
  // 0), fespace_b for H(div) fields (source index 1). Either may be nullptr depending
  // on the functional kind. Material operator (not owned) for material property lookups
  // and side selection.
  const FiniteElementSpace *fespace_e;
  const FiniteElementSpace *fespace_b;
  const MaterialOperator *mat_op;

  // MPI communicator from the mesh.
  MPI_Comm comm;

  // Per-group assembled libCEED operators. Each operator integrates over one group of
  // boundary elements and accumulates per-element integrals into the local output
  // vector (CeedOperatorApplyAdd with all field inputs passive). The field inputs
  // (QFunction input name, source vector index) are re-pointed at the caller's data on
  // each evaluation.
  struct GroupOp
  {
    Ceed ceed;
    CeedOperator op;
    std::vector<std::pair<std::string, int>> field_sources;
  };
  std::vector<GroupOp> groups;

  // Staging vector used to initialize the field input CeedVectors at construction. The
  // field CeedVectors are re-pointed at the caller's data on each Eval() call.
  mutable Vector field_staging;

  // Local output vector with one slot per marked boundary element on this process.
  mutable Vector local_out;

  void Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);

  // Apply all group operators with the field inputs pointed at the given source
  // vectors, accumulating into the local output vector.
  void ApplyAdd(const std::array<const Vector *, 2> &srcs) const;

  // Zero the local output vector, apply, and return the local sum (no MPI reduction).
  double EvalLocal(const std::array<const Vector *, 2> &srcs) const;

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

  // Construct a surface flux functional (see BdrSurfaceFluxCoefficient). The required
  // finite element spaces depend on the flux type: ELECTRIC requires nd_fespace,
  // MAGNETIC requires rt_fespace, POWER requires both.
  SurfaceFunctional(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const FiniteElementSpace *nd_fespace,
                    const FiniteElementSpace *rt_fespace, const MaterialOperator &mat_op,
                    SurfaceFlux type, bool two_sided, const mfem::Vector &x0);

  ~SurfaceFunctional();

  SurfaceFunctional(const SurfaceFunctional &) = delete;
  SurfaceFunctional &operator=(const SurfaceFunctional &) = delete;

  // Evaluate the functional for the given field (L-vector, e.g. the local vector of a
  // GridFunction on the field space). Collective on the mesh communicator. For
  // field-less functionals, u is ignored and may be nullptr.
  double Eval(const Vector *u = nullptr) const;

  // Evaluate the functional for the given (possibly complex-valued) grid function. For
  // complex fields, the real and imaginary part contributions add (the implemented
  // integrands are quadratic in the field). Collective on the mesh communicator.
  double Eval(const GridFunction &u) const;

  // Evaluate a surface flux functional for the given fields (either of which may be
  // nullptr if not required by the flux type). For complex-valued fields, returns the
  // real and imaginary parts of the flux (ELECTRIC, MAGNETIC), or the stationary real
  // power (POWER). Collective on the mesh communicator.
  std::complex<double> EvalFlux(const GridFunction *E, const GridFunction *B) const;
};

}  // namespace palace

#endif  // PALACE_FEM_SURFACE_FUNCTIONAL_HPP
