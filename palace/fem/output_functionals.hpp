// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_OUTPUT_FUNCTIONALS_HPP
#define PALACE_FEM_OUTPUT_FUNCTIONALS_HPP

#include <array>
#include <complex>
#include <deque>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/ceed_group_operator.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

union CeedIntScalar;

namespace palace::ceed
{

struct CeedQFunctionInfo;

}  // namespace palace::ceed

namespace palace
{

class BoundaryPhysicalTraceCache;
class FaceSamplingPlan;
class GridFunction;
class MaterialOperator;
class Mesh;
class PointFieldEvaluator;
class FaceNbrFieldExchange;

// Description of a real vector-valued mode coefficient on a marked boundary surface.
// UNIFORM represents scale * direction. COAXIAL represents scale * (x - origin) /
// |x - origin|^2, matching CoaxialElementData::GetModeCoefficient.
struct SurfaceModeCoefficient
{
  enum class Type
  {
    UNIFORM,
    COAXIAL
  };

  Type type = Type::UNIFORM;
  mfem::Array<int> attr_list;
  double scale = 1.0;
  std::array<double, 3> direction = {0.0, 0.0, 0.0};
  std::array<double, 3> origin = {0.0, 0.0, 0.0};
};

//
// Class to compute reducing output functionals (integrals of functions of solution
// fields) over boundary element sets using libCEED, supporting full (non-trace)
// evaluation of volume fields at boundary element quadrature points. Non-reducing
// visualization point fields are exposed through PointFieldEvaluator, which uses private
// boundary point-field assembly hooks here. This enables postprocessing measurements
// (interface dielectric energy participation, surface fluxes, port powers, etc.) to
// execute on the device instead of using host-only mfem::Coefficient evaluation for
// supported paths.
//
// The key construction: for each boundary element, the field is evaluated from an
// attached volume element (or both, with averaging or differencing, for interior
// boundaries, following the conventions of BdrGridFunctionCoefficient and its derived
// coefficient semantics). Groups use fixed mapped integration rules identified by the
// exact local face/subface transformation matrix, rather than rounded coordinates or
// face identity. Processor-boundary ghost values use the same canonical keys through
// FaceNbrFieldExchange.
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
    FARFIELD        // Stratton-Chu far-field following AddStrattonChuIntegrandAtElement
  };

private:
  friend class BoundaryPhysicalTraceCache;
  friend class PointFieldEvaluator;

  // Internal backend selector. Public SurfaceFunctional::Kind is reduction-only;
  // non-reducing boundary visualization entries are reachable only through
  // PointFieldEvaluator's private backend hooks.
  enum class KernelKind
  {
    AREA,
    HCURL_NORM2,
    INTERFACE_EPR,
    SURFACE_FLUX,
    FARFIELD,
    MODE_OVERLAP,
    BDR_FIELD_E,
    BDR_FIELD_B,
    BDR_FIELD_H1,
    BDR_FLUX_Q,
    BDR_CURRENT_J,
    BDR_ENERGY_E,
    BDR_ENERGY_M,
    BDR_POYNTING
  };

  static KernelKind ToKernelKind(Kind kind);
  static KernelKind ToKernelKind(PointFieldKind kind);
  static const char *KindName(KernelKind kind);

  // Whether the backend kind fills a per-point visualization buffer (vs. computing
  // reductions).
  static bool IsBufferKind(KernelKind kind)
  {
    return kind == KernelKind::BDR_FIELD_E || kind == KernelKind::BDR_FIELD_B ||
           kind == KernelKind::BDR_FIELD_H1 || kind == KernelKind::BDR_FLUX_Q ||
           kind == KernelKind::BDR_CURRENT_J || kind == KernelKind::BDR_ENERGY_E ||
           kind == KernelKind::BDR_ENERGY_M || kind == KernelKind::BDR_POYNTING;
  }

  // Number of components per visualization point for buffer kinds.
  static int BufferNumComp(KernelKind kind, int sdim)
  {
    return (kind == KernelKind::BDR_FIELD_H1 || kind == KernelKind::BDR_FLUX_Q ||
            kind == KernelKind::BDR_ENERGY_E || kind == KernelKind::BDR_ENERGY_M)
               ? 1
               : sdim;
  }

  // Total buffer size (all boundary elements, lattice points, components) and
  // per-element point-base offsets for the boundary visualization field kinds. Vector
  // buffers are component-major: x[points], y[points], (z[points] in 3D).
  int BufferSize() const { return buffer_size; }
  int BufferNumComp() const { return buffer_num_comp; }
  const std::vector<int> &BufferBases() const { return buffer_bases; }

  // Computation kind and integrand parameters.
  KernelKind kind;
  InterfaceDielectric epr_type = InterfaceDielectric::DEFAULT;
  double epr_t = 0.0, epr_epsilon = 0.0;
  SurfaceFlux flux_type = SurfaceFlux::ELECTRIC;
  bool flux_two_sided = false;
  mfem::Vector flux_x0;
  std::vector<std::array<double, 3>> farfield_dirs;
  std::complex<double> farfield_omega = 0.0;
  std::map<int, SurfaceModeCoefficient> mode_coeff_by_attr;

  // Boundary visualization field kinds: lattice refinement level, output scaling,
  // total output buffer size, and per-boundary-element point-base offsets into the
  // component-major buffer.
  int viz_lod = 0;
  double viz_scaling = 1.0;
  int buffer_size = 0;
  int buffer_num_comp = 0;
  std::vector<int> buffer_bases;

  // Field finite element spaces (not owned): nd_fespace for H(curl)/H1 fields (source
  // index 0), rt_fespace for H(div) fields (source index 1). Either may be nullptr
  // depending on the functional kind. Material operator (not owned) for material property
  // lookups and side selection.
  const mfem::ParFiniteElementSpace *nd_fespace;
  const mfem::ParFiniteElementSpace *rt_fespace;
  const MaterialOperator *mat_op;

  // Whether the functional could be assembled. False means the configuration is outside
  // the current support matrix; supported model-level paths treat invalid assembly as an
  // error rather than silently selecting a different implementation.
  bool valid = true;

  // MPI communicator from the mesh.
  MPI_Comm comm;

  // Per-group assembled libCEED operators, accumulating per-element integrals into the
  // local output vector. The element attribute vectors are operator inputs and must
  // outlive the operators.
  std::vector<fem::CeedGroupOperator> groups;
  std::deque<Vector> elem_attrs;

  // Face neighbor field exchange for two-sided interior boundaries crossing parallel
  // interfaces: the owning process pulls the neighbor (ghost) volume field values at the
  // face quadrature points, fed to the ghost side of the two-sided operators (nullptr
  // when no marked boundary element has a ghost neighbor). Refilled before each apply.
  std::unique_ptr<FaceNbrFieldExchange> face_nbr_exchange;

  // Shared immutable boundary sampling setup for private point-field evaluators only.
  // Public reduction functional constructors intentionally do not accept this plan.
  std::shared_ptr<const FaceSamplingPlan> sampling_plan;

  // Boundary point evaluators supplied by BoundaryPhysicalTraceCache consume the
  // already Piola-mapped side traces through EVAL_NONE restrictions. Trace producers
  // set trace_side (0/1) and keep their private buffers in canonical point order.
  std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache;
  int trace_side = -1;

  // Staging vector used to initialize the field input CeedVectors at construction. The
  // field CeedVectors are re-pointed at the caller's data on each Eval() call.
  mutable Vector field_staging;

  // Local output vector with one slot per marked boundary element on this process.
  // Integral functionals also keep the originating boundary attribute for each slot so
  // batched model-level callers can recover several independent reductions from one
  // assembled operator without changing the per-element libCEED kernels.
  mutable Vector local_out;
  std::vector<int> local_out_attrs;

  void Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);
  std::vector<CeedIntScalar> BuildBaseContext(int dim, bool is_2d) const;
  void ConfigureGroupQFunction(bool is_2d, bool has_b, double side_scale,
                               double normal_scale, int bdr_attr,
                               const std::vector<CeedIntScalar> &base_ctx,
                               std::vector<CeedIntScalar> &ctx,
                               ceed::CeedQFunctionInfo &info) const;
  void AssembleLocal(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker);
  void WarmUpBufferOperators() const;

  // Apply all group operators with the field inputs pointed at the given source
  // vectors, accumulating into the local output vector.
  void ApplyAdd(const std::array<const Vector *, 4> &srcs) const;

  // Zero the local output vector, apply, and return the local sum (no MPI reduction).
  double EvalLocal(const std::array<const Vector *, 4> &srcs) const;

  // Add the current local output element slots into per-attribute bins.
  void BinLocalOutByAttribute(const mfem::Array<int> &attr_to_bin, int num_bins,
                              std::vector<double> &bins, double scale) const;

  // Construct boundary point-field evaluators. These are intentionally private to keep
  // SurfaceFunctional reduction-oriented at call sites; PointFieldEvaluator owns the
  // non-reducing visualization API.
  SurfaceFunctional(PointFieldKind kind, const Mesh &mesh,
                    const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &fespace, int lod,
                    std::shared_ptr<const FaceSamplingPlan> sampling_plan,
                    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr,
                    int trace_side = -1);
  SurfaceFunctional(PointFieldKind kind, const Mesh &mesh,
                    const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &fespace,
                    const MaterialOperator &mat_op, int lod, double scaling,
                    std::shared_ptr<const FaceSamplingPlan> sampling_plan,
                    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr);
  SurfaceFunctional(PointFieldKind kind, const Mesh &mesh,
                    const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace &nd_fespace,
                    const mfem::ParFiniteElementSpace &rt_fespace,
                    const MaterialOperator &mat_op, int lod, double scaling,
                    std::shared_ptr<const FaceSamplingPlan> sampling_plan,
                    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr);

  // Fill boundary visualization buffers. Friend-only; non-reducing callers use
  // PointFieldEvaluator.
  void EvalBuffer(const Vector &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &E, const GridFunction &B, Vector &buffer) const;

  // The trace overload receives canonical, component-major physical side buffers.
  // It is private because PointFieldEvaluator controls source identity and cache
  // generation.
  void EvalBufferTrace(const std::array<const Vector *, 4> &traces, Vector &buffer,
                       bool zero = true) const;
  std::size_t TraceGroupCount() const { return groups.size(); }

public:
  // Construct a functional over the boundary elements with marked attributes (marker
  // over global mfem boundary attributes). For field-less functionals (AREA), fespace
  // may be nullptr but the mesh is still required.
  SurfaceFunctional(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                    const mfem::ParFiniteElementSpace *fespace = nullptr);

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

  // Construct a linear H(curl) mode-overlap functional, ∫ E · f_mode dS, for lumped-port
  // voltage and S-parameter extraction. The mode coefficient data is real-valued; complex
  // fields are evaluated by applying the same functional to the real and imaginary parts.
  SurfaceFunctional(const Mesh &mesh, const mfem::ParFiniteElementSpace &nd_fespace,
                    const std::vector<SurfaceModeCoefficient> &mode_coeffs);

  ~SurfaceFunctional();

  SurfaceFunctional(const SurfaceFunctional &) = delete;
  SurfaceFunctional &operator=(const SurfaceFunctional &) = delete;

  // Whether the functional was successfully assembled. When false, evaluation is not
  // possible; supported model-level paths should error rather than silently selecting a
  // different implementation.
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

  // Evaluate the complex power P = ∫ (E x H*) ⋅ n dS = -∫ E ⋅ (n x H*) dS following
  // the conventions of LumpedPortData::GetPower (two-sided POWER flux functionals only).
  // Collective on the mesh communicator.
  std::complex<double> EvalComplexPower(const GridFunction &E, const GridFunction &B) const;

  // Same complex-power convention as EvalComplexPower, but return one result per
  // boundary-attribute bin. attr_to_bin is indexed by boundary attribute - 1 and uses -1
  // for attributes not assigned to any output bin. This enables safe batching of many
  // disjoint port surfaces while retaining per-port scalar outputs. Collective on the
  // mesh communicator.
  std::vector<std::complex<double>>
  EvalComplexPowerByAttribute(const GridFunction &E, const GridFunction &B,
                              const mfem::Array<int> &attr_to_bin, int num_bins) const;

  // Evaluate the far-field rE integrals for all observation directions at the given
  // (complex) frequency, following SurfacePostOperator::GetFarFieldrE. Reassembles when
  // the frequency changes. Collective on the mesh communicator.
  std::vector<std::array<std::complex<double>, 3>>
  EvalFarField(const GridFunction &E, const GridFunction &B, std::complex<double> omega);

  // Evaluate the linear mode-overlap functional, returning the complex integral when the
  // supplied field has real and imaginary parts.
  std::complex<double> EvalModeOverlap(const GridFunction &E) const;

  // Same mode-overlap convention as EvalModeOverlap, but return one result per
  // boundary-attribute bin. attr_to_bin is indexed by boundary attribute - 1 and uses -1
  // for attributes not assigned to an output bin. Collective on the mesh communicator.
  std::vector<std::complex<double>>
  EvalModeOverlapByAttribute(const GridFunction &E, const mfem::Array<int> &attr_to_bin,
                             int num_bins) const;
};

}  // namespace palace

#endif  // PALACE_FEM_OUTPUT_FUNCTIONALS_HPP
