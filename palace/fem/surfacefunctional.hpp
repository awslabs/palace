// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SURFACE_FUNCTIONAL_HPP
#define PALACE_FEM_SURFACE_FUNCTIONAL_HPP

#include <array>
#include <complex>
#include <deque>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class GridFunction;
class MaterialOperator;
class Mesh;

namespace fem
{

// An assembled libCEED operator over one group of (boundary) elements, with the named
// passive field inputs re-pointed at caller data (by source vector index) on each
// evaluation.
struct CeedGroupOperator
{
  Ceed ceed;
  CeedOperator op;
  std::vector<std::pair<std::string, int>> field_sources;
};

}  // namespace fem

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
    SURFACE_FLUX,   // Surface flux following BdrSurfaceFluxCoefficient
    FARFIELD,       // Stratton-Chu far-field following AddStrattonChuIntegrandAtElement
    BDR_FIELD_E,    // H(curl) field values at boundary visualization points
    BDR_FIELD_B,    // H(div) field values at boundary visualization points
    BDR_FLUX_Q,     // Surface charge (eps E) . n at boundary visualization points
    BDR_CURRENT_J,  // Surface current n x (mu^-1 B) at boundary visualization points
    BDR_ENERGY_E,   // Electric energy density at boundary visualization points
    BDR_ENERGY_M    // Magnetic energy density at boundary visualization points
  };

  // Whether the kind fills a per-point visualization buffer (vs. computing integrals).
  static bool IsBufferKind(Kind kind)
  {
    return kind == Kind::BDR_FIELD_E || kind == Kind::BDR_FIELD_B ||
           kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_CURRENT_J ||
           kind == Kind::BDR_ENERGY_E || kind == Kind::BDR_ENERGY_M;
  }

  // Number of components per visualization point for buffer kinds.
  static int BufferNumComp(Kind kind)
  {
    return (kind == Kind::BDR_FLUX_Q || kind == Kind::BDR_ENERGY_E ||
            kind == Kind::BDR_ENERGY_M)
               ? 1
               : 3;
  }

  // Total buffer size (all boundary elements, lattice points, components) and
  // per-element base offsets for the boundary visualization field kinds.
  int BufferSize() const { return buffer_size; }
  const std::vector<int> &BufferBases() const { return buffer_bases; }

private:
  // Computation kind and integrand parameters.
  Kind kind;
  InterfaceDielectric epr_type = InterfaceDielectric::DEFAULT;
  double epr_t = 0.0, epr_epsilon = 0.0;
  SurfaceFlux flux_type = SurfaceFlux::ELECTRIC;
  bool flux_two_sided = false;
  mfem::Vector flux_x0;
  std::vector<std::array<double, 3>> farfield_dirs;
  double farfield_omega_re = 0.0, farfield_omega_im = 0.0;

  // Mesh and marker for reassembly when the far-field frequency changes (the frequency
  // enters the QFunction context).
  const Mesh *farfield_mesh = nullptr;
  mfem::Array<int> farfield_marker;

  // Boundary visualization field kinds: lattice refinement level, output scaling,
  // total output buffer size, and per-boundary-element base offsets into the buffer.
  int viz_lod = 0;
  double viz_scaling = 1.0;
  int buffer_size = 0;
  std::vector<int> buffer_bases;

  // Field finite element spaces (not owned): fespace_e for H(curl) fields (source index
  // 0), fespace_b for H(div) fields (source index 1). Either may be nullptr depending
  // on the functional kind. Material operator (not owned) for material property lookups
  // and side selection.
  const mfem::ParFiniteElementSpace *fespace_e;
  const mfem::ParFiniteElementSpace *fespace_b;
  const MaterialOperator *mat_op;

  // Whether the functional could be assembled (false when the configuration is not yet
  // supported, e.g. non-3D meshes or two-sided evaluation on process-boundary interior
  // surfaces, in which case callers should fall back to the legacy evaluation paths).
  bool valid = true;

  // MPI communicator from the mesh.
  MPI_Comm comm;

  // Per-group assembled libCEED operators, accumulating per-element integrals into the
  // local output vector. The element attribute vectors are operator inputs and must
  // outlive the operators.
  std::vector<fem::CeedGroupOperator> groups;
  std::deque<Vector> elem_attrs;

  // Staging vector used to initialize the field input CeedVectors at construction. The
  // field CeedVectors are re-pointed at the caller's data on each Eval() call.
  mutable Vector field_staging;

  // Local output vector with one slot per marked boundary element on this process.
  mutable Vector local_out;

  void Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);
  void AssembleLocal(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);

  // Apply all group operators with the field inputs pointed at the given source
  // vectors, accumulating into the local output vector.
  void ApplyAdd(const std::array<const Vector *, 4> &srcs) const;

  // Zero the local output vector, apply, and return the local sum (no MPI reduction).
  double EvalLocal(const std::array<const Vector *, 4> &srcs) const;

public:
  // Returns false when libCEED surface functionals have been globally disabled via the
  // PALACE_LEGACY_SURFACE_POSTPRO environment variable (legacy mfem::Coefficient paths
  // are used instead, for debugging and benchmarking).
  static bool Enabled();

  // Construct a functional over the boundary elements with marked attributes (marker
  // over global mfem boundary attributes). For field-less functionals (AREA), fespace
  // may be nullptr but the mesh is still required.
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace *fespace = nullptr);

  // Construct a boundary visualization field evaluator (BDR_FIELD_E or BDR_FIELD_B),
  // evaluating at the order-lod lattice points of each boundary element (the ParaView
  // output sampling, see mfem::RefinedGeometry).
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &fespace, int lod);

  // Construct a boundary visualization field evaluator with material properties and
  // output scaling (BDR_FLUX_Q, BDR_CURRENT_J, BDR_ENERGY_E, BDR_ENERGY_M).
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &fespace,
                    const MaterialOperator &mat_op, int lod, double scaling);

  // Construct an interface dielectric energy participation functional with the given
  // interface type, thickness, and permittivity (see InterfaceDielectricCoefficient).
  SurfaceFunctional(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &nd_fespace,
                    const MaterialOperator &mat_op, InterfaceDielectric type, double t_i,
                    double epsilon_i);

  // Construct a surface flux functional (see BdrSurfaceFluxCoefficient). The required
  // finite element spaces depend on the flux type: ELECTRIC requires nd_fespace,
  // MAGNETIC requires rt_fespace, POWER requires both.
  SurfaceFunctional(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace *nd_fespace,
                    const mfem::ParFiniteElementSpace *rt_fespace,
                    const MaterialOperator &mat_op, SurfaceFlux type, bool two_sided,
                    const mfem::Vector &x0);

  // Construct a Stratton-Chu far-field functional for the given observation directions
  // (see AddStrattonChuIntegrandAtElement; external boundaries only).
  SurfaceFunctional(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &nd_fespace,
                    const mfem::ParFiniteElementSpace &rt_fespace,
                    const MaterialOperator &mat_op,
                    const std::vector<std::array<double, 3>> &r_naughts);

  ~SurfaceFunctional();

  SurfaceFunctional(const SurfaceFunctional &) = delete;
  SurfaceFunctional &operator=(const SurfaceFunctional &) = delete;

  // Whether the functional was successfully assembled. When false, evaluation is not
  // possible and callers should use the legacy evaluation paths.
  bool IsValid() const { return valid; }

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

  // Evaluate the complex power P = ∫ E ⋅ (n x H) dS following the conventions of
  // LumpedPortData::GetPower (two-sided POWER flux functionals only). Collective on the
  // mesh communicator.
  std::complex<double> EvalComplexPower(const GridFunction &E, const GridFunction &B) const;

  // Evaluate the far-field rE integrals for all observation directions at the given
  // (complex) frequency, following SurfacePostOperator::GetFarFieldrE. Reassembles when
  // the frequency changes. Collective on the mesh communicator.
  std::vector<std::array<std::complex<double>, 3>> EvalFarField(const GridFunction &E,
                                                                const GridFunction &B,
                                                                double omega_re,
                                                                double omega_im);

  // Fill the boundary visualization buffer with the pointwise field values (local
  // operation, buffer kinds only). The grid function overload accumulates the real and
  // imaginary part contributions (energy density kinds).
  void EvalBuffer(const Vector &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &u, Vector &buffer) const;
};

//
// Class to evaluate derived field quantities (energy densities, Poynting vector) at the
// nodal points of an interpolatory output space using libCEED, filling a grid function
// for visualization output without per-point host coefficient evaluation. Follows the
// conventions of EnergyDensityCoefficient and PoyntingVectorCoefficient.
//
class DomainFieldEvaluator
{
public:
  enum class Kind
  {
    ENERGY_E,  // 1/2 (ε E)ᴴ E, scalar output
    ENERGY_M,  // 1/2 (μ⁻¹ B)ᴴ B, scalar output
    POYNTING   // Re{E x (μ⁻¹ B)⥁}, vector output
  };

private:
  Kind kind;

  // Source field finite element spaces (not owned): fespace_e for H(curl), fespace_b
  // for H(div); either may be nullptr depending on the kind.
  const mfem::ParFiniteElementSpace *fespace_e;
  const mfem::ParFiniteElementSpace *fespace_b;

  // Whether the evaluator could be assembled (3D meshes only).
  bool valid = true;

  // Per-geometry assembled libCEED operators, evaluating at the target space nodal
  // points and scattering directly into the output grid function. The element attribute
  // vectors are operator inputs and must outlive the operators.
  std::vector<fem::CeedGroupOperator> groups;
  std::deque<Vector> elem_attrs;

  // Staging vector used to initialize the field input CeedVectors at construction.
  mutable Vector field_staging;

  void Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                const mfem::ParFiniteElementSpace &target_fespace, double scaling);

public:
  // Construct an evaluator filling grid functions on target_fespace (an interpolatory
  // L2 space). The scaling multiplies the output as for the legacy coefficients.
  DomainFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                       const mfem::ParFiniteElementSpace *nd_fespace,
                       const mfem::ParFiniteElementSpace *rt_fespace,
                       const mfem::ParFiniteElementSpace &target_fespace, double scaling);
  ~DomainFieldEvaluator();

  DomainFieldEvaluator(const DomainFieldEvaluator &) = delete;
  DomainFieldEvaluator &operator=(const DomainFieldEvaluator &) = delete;

  // Whether the evaluator was successfully assembled.
  bool IsValid() const { return valid; }

  // Fill the output vector (L-vector of the target space, e.g. a GridFunction) with
  // the pointwise quantity. Real and imaginary field contributions add. Local
  // operation (no MPI communication).
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;
};

}  // namespace palace

#endif  // PALACE_FEM_SURFACE_FUNCTIONAL_HPP
