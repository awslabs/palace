// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP

#include <array>
#include <complex>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/surfaceresponsemirror.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class FiniteElementSpace;
class GridFunction;
class IoData;
class Mesh;
class BoundaryModeOperator;
class DistributedPointLocator;
class LaplaceOperator;
class MaterialOperator;
class SpaceOperator;
struct IdentifiedFeature;
struct IdentificationResult;

// The placement's clip of the uncovered requirements (decision 394 F2) by the matched
// spatial clusters' support boxes (decision 399 MAJOR-1): the part of an uncovered portion
// strictly inside a matched cluster's box (the same bounds and strict-interior test the
// continuation ownership clips the translational cells with, CellInsideBox) is removed —
// the coupon models its whole box, so that raw within-R energy would be counted twice.
// One record per (portion, box); lengths in patch units.
struct UncoveredSpatialSupportClip
{
  int feature = -1;
  std::string topology;
  int segment = -1;
  std::size_t spatial_patch = 0;
  double length = 0.0;  // removed from the portion by this box
  // The removed piece's ends (mesh units): the F-DB-a footprint of the box's coupon when
  // that coupon is a DomainBoundary exclusion (DESIGN 2.1).
  std::array<double, 3> p0{}, p1{};
};
struct UncoveredSpatialSupportClipping
{
  std::vector<UncoveredSpatialSupportClip> clips;  // in portion order, then box order
  int clipped_portions = 0;  // portions that lost a part (wholly removed ones included)
  int removed_portions = 0;  // portions wholly inside the boxes (nothing kept)
  int split_portions = 0;    // portions whose kept part is two pieces (a box crossed)
  double removed_length = 0.0;
  std::map<int, double> removed_by_feature;  // feature id -> removed length
};

// F-DB-a (decisions 442 / 454, DESIGN 2.1): the raw claims of the DomainBoundary-excluded
// patches as perimeter portions (see CollectDomainBoundaryPortions).
struct DomainBoundaryPortions
{
  std::vector<config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
      portions;
  // Per portion (parallel): the excluded patches it came from (ascending; two for a
  // deduplicated first-order split cell).
  std::vector<std::vector<std::size_t>> patches;
  // Per portion (parallel): a translational own-edge interval (a geometric cell).
  std::vector<bool> translational;
  int geometric_cells = 0;    // translational own-edge intervals, split patches once
  int duplicate_patches = 0;  // excluded patches folded into an existing interval
  double length = 0.0;        // mesh units, over all portions
};

// Mesh-independent automatic coupon layout. A solver retains this across AMR iterations;
// finite-element point interpolation is still rebuilt for every refined mesh.
class SurfaceResponseGeometry
{
private:
  struct Impl;
  std::shared_ptr<const Impl> impl;

  explicit SurfaceResponseGeometry(std::shared_ptr<const Impl> impl_)
    : impl(std::move(impl_))
  {
  }

  friend class SurfaceResponseOperator;
};

//
// Low-rank Schur-complement correction which replaces the local response of an ideal
// thin-metal edge by a fabrication-resolved coupon response. The correction has the form
//
//                         C = Pᵀ (S_fabricated - S_thin) P,
//
// where P evaluates the global potential relative to the metal potential at coupon
// contour knots. The assembled thin-metal Laplace operator remains the preconditioner.
//
class SurfaceResponseOperator : public Operator
{
public:
  struct PatchAssignment
  {
    int model = 0;
    std::array<double, 3> origin{};
    std::array<double, 3> axis_u{};
    std::array<double, 3> axis_v{};
    std::array<double, 3> axis_w{};
    double weight = 0.0;
    // The patch's global (placed, 0-based) index: the key of the solve-time records.
    int global_index = -1;
  };

  // Conductor-consistency gate (decision 277 (A)): the record of one spatial patch at one
  // excitation. A spatial coupon holds its metal cross-sections on the box faces at the
  // conductor potentials, while the surface mortar imposes the DEVICE potential on the free
  // knots around them; where the coupon's metal is not the device's (a straight 3R
  // continuation past a device corner, a mis-keyed or misplaced coupon) the device
  // potential at the coupon's metal differs from the conductor's and the coupon's energies
  // under the device trace are wrong (the S1p 19-edge MS surplus, decisions 271-277). The
  // probe: the device potential at every conductor vertex of the trace mesh ON THE PROCESS
  // PLANE (w = 0, the coupon's metal bottom = the device's metal sheet, a Dirichlet surface
  // for real metal) minus the potential at that conductor's reference point, relative to
  // the patch's trace amplitude max(|c_i|, |V_state|) floored at
  // kConductorConsistencyAmplitudeFloor x the excitation's largest potential. Real metal
  // reads ~0 (the knot lies on the device conductor), fictitious metal reads the device gap
  // potential (S1p: 0.20 to 0.43). Information only (recorded, not gated): the conductor
  // vertices OFF the plane
  // (the coupon's metal top, 0.1 um into the device's gap vacuum for a thin device: reads
  // the normal field x the thickness, S1p 0.013-0.028) and the ADJACENT free knots (sharing
  // a trace-triangle edge with a plane conductor vertex in the same column, the 50-nm
  // trench-floor row): their mortar coefficient minus the conductor value reads the
  // near-edge field x 50 nm, up to 0.27 of the amplitude on the transmon's real metal, so
  // it cannot discriminate (the calibration of conductor-consistency-20261003).
  struct ConductorConsistencyRecord
  {
    int patch = -1;  // global 0-based index
    int model = 0;   // runtime model index (ModelCatalog)
    int source = 0;  // the excitation index
    int plane_knots = 0;
    int off_plane_knots = 0;
    int adjacent_knots = 0;
    double amplitude = 0.0;  // max |trace coefficient| incl. the conductor states
    double state = 0.0;      // max |conductor state| (0 for a single-conductor coupon)
    // The excitation's largest potential, max |x| over the device (the largest terminal
    // potential: the maximum principle), and the ratio's denominator max(amplitude,
    // kConductorConsistencyAmplitudeFloor x excitation_potential).
    double excitation_potential = 0.0;
    double normalization = 0.0;
    double max_deviation = 0.0;  // max |V_device(knot) - V_conductor| over the plane knots
    double max_ratio = 0.0;      // max_deviation / normalization
    int worst_vertex = -1;       // trace-mesh vertex (0-based) of max_deviation
    int worst_conductor = 0;     // its conductor (1-based)
    // Nondimensional mesh coordinates; the record and the log multiply by the mesh
    // coordinate scale (mesh-file units = the device coordinates).
    std::array<double, 3> worst_point{};
    double off_plane_max_ratio = 0.0;
    double adjacent_max_ratio = 0.0;
    double claim_length =
        0.0;  // the cluster's claimed portions (nondimensional; 0 otherwise)
    double cell_length =
        0.0;  // the longitudinal cell (nondimensional; 0 for a spatial patch)
    bool excluded = false;
  };
  // The dimensionless tolerance of the gate: max_ratio above it excludes the patch (weight
  // 0, like a DomainBoundary cell). Calibrated (conductor-consistency-20261003/gate): the
  // plane knots of real metal read <= 7e-7 on the S1p c0 field (31 knots) and exactly the
  // Dirichlet value in the operator; the two S1p fictitious blocks read 0.20-0.43 of the
  // amplitude (0.22-0.48 of the state); 0.02 is x10 below the weakest flagged knot and
  // ~3e4 above the measured real-metal maximum, and of the order of a metal-thickness
  // (0.1 um) offset normal to the plane. A 5 % potential mismatch would already shift the
  // 15 %-residual MS closure by O(several points), so the tolerance is not larger.
  static constexpr double kConductorConsistencyTolerance = 0.02;
  // The normalization floor of the ratio (decision 279, MINOR-3): the denominator is
  // max(amplitude, kConductorConsistencyAmplitudeFloor x the excitation's largest
  // potential), so a patch whose trace amplitude is a near-zero fraction of the excitation
  // (noise over noise) cannot be excluded spuriously; the record notes FloorApplied. A
  // patch below the floor carries at most floor^2 = 1e-6 of a unit-amplitude patch's
  // energy (energy ~ amplitude^2), and a defective coupon is still excluded at the
  // excitation that fields it (the exclusion is sticky). The S1p flagged coupons read
  // amplitude ~1.1 x the 1-V state, far above the floor.
  static constexpr double kConductorConsistencyAmplitudeFloor = 1.0e-3;

private:
  struct MaxwellLine
  {
    int point_offset = 0;
    int point_count = 0;
  };

  struct MaxwellLineGeometry
  {
    std::array<double, 3> begin{};
    std::array<double, 3> end{};
  };

  struct MaxwellContourPath
  {
    int anchor_line = -1;
    int contour_line_offset = 0;
    int contour_line_count = 0;
    int end_line = -1;
    int start_conductor_trace = -1;
    int end_conductor_trace = -1;
    bool closed = true;
    std::vector<int> trace_indices;
  };

  struct MaxwellConductorPath
  {
    int contour_path = -1;
    int parent_trace_offset = -1;
    int trace_offset = 0;
    int integral_sign = 0;
  };

  struct OpenContourPath
  {
    std::vector<int> indices;
    int start_conductor = 0;
    int end_conductor = 0;
  };

  enum class DomainCorrectionMode : char
  {
    DISABLED,
    FIXED_TRACE,
    FIXED_FLUX
  };

  struct MortarSegment
  {
    int begin = 0;
    int end = 0;
    double length = 0.0;
    int subdivisions = 1;
  };

  // A vertex of a spatial model's trace triangulation: a basis knot (basis >= 0, weight 1),
  // a conductor vertex (basis < 0, conductor > 0) or a SLAVE vertex (corner-family trace
  // basis rule): a geometric vertex of the box (a box corner that is no knot) whose trace
  // is the linear interpolation between two knots, `weight` on `basis` and 1 - weight on
  // `second_basis`; its hat contributions are those of its parents.
  struct MortarVertex
  {
    std::array<double, 3> point{};
    int basis = -1;
    int conductor = 0;
    int second_basis = -1;
    double weight = 1.0;

    // The (basis, hat value) pairs the vertex carries.
    template <typename F>
    void ForEachBasis(F &&f) const
    {
      if (basis >= 0)
      {
        f(basis, weight);
      }
      if (second_basis >= 0)
      {
        f(second_basis, 1.0 - weight);
      }
    }
  };

  struct MortarTriangle
  {
    std::array<int, 3> vertices{};
    double area = 0.0;
    double maximum_edge_length = 0.0;
  };

  struct ResponseModel
  {
    int idx = 0;
    std::string name;
    std::string topology;
    int contour_size = 0;
    int basis_size = 0;
    int conductor_state_count = 0;
    mfem::DenseMatrix fabricated_domain;
    mfem::DenseMatrix thin_domain;
    mfem::DenseMatrix domain_defect;
    mfem::DenseMatrix fixed_flux_transform;
    mfem::DenseMatrix fixed_flux_domain_defect;
    DomainCorrectionMode domain_correction_mode = DomainCorrectionMode::FIXED_TRACE;
    std::map<int, mfem::DenseMatrix> fabricated_surfaces;
    std::map<int, mfem::DenseMatrix> surface_defects;
    bool spatial_basis = false;
    std::vector<int> contour_groups;
    std::vector<int> zero_trace_indices;
    std::vector<OpenContourPath> open_contour_paths;
    // Trailing basis points in the interior of the matching-box caps (on no contour).
    int interior_trace_count = 0;
    bool surface_mortar = false;
    bool spatial_mortar = false;
    // Translational surface mortar: the hat segments between consecutive contour vertices.
    // A segment end >= contour_size is a constrained metal-band vertex
    // consistent_mortar_vertices[end - contour_size] (decision 404 D1): no basis
    // coefficient, no load, so the adjacent free hat ramps to zero there; the band between
    // two such vertices carries no segment.
    std::vector<MortarSegment> mortar_segments;
    std::vector<std::array<double, 3>> consistent_mortar_vertices;
    std::string consistent_mortar_rule;
    // A translational mortar whose hat basis has constrained vertices (inserted band
    // vertices or band knots listed as ZeroTraceIndices): the free hats do not sum to one
    // beside them, so the projection removes the reference through the hats' integrals.
    bool constrained_translational_mortar = false;
    std::vector<MortarVertex> mortar_vertices;
    std::vector<MortarTriangle> mortar_triangles;
    mfem::DenseMatrix mortar_mass_inverse;
    mfem::Vector mortar_constant_load;
    std::vector<mfem::Vector> mortar_conductor_loads;
    // Conductor-consistency gate (decision 277): the conductor vertices of the trace mesh
    // on the process plane (the probes of the gate), those off it (the metal top rows,
    // recorded) and the free knots adjacent to a plane conductor vertex in its column
    // (free vertex, plane conductor vertex; recorded).
    std::vector<int> plane_conductor_vertices;
    std::vector<int> off_plane_conductor_vertices;
    std::vector<std::pair<int, int>> adjacent_free_vertices;
  };

  struct Patch
  {
    int global_index = -1;
    int model = -1;
    int point_offset = 0;
    int trace_offset = 0;
    int point_count = 0;
    // Translational surface mortar: the longitudinal strip [begin, end] (offsets along the
    // patch's AxisW, mesh units; the patch's cell of its feature portion) is sampled in
    // mortar_longitudinal_subdivisions cross-sections of the local mortar resolution; a
    // {0, 0} strip is one cross-section at the origin (2D, spatial, explicit patches).
    std::array<double, 2> mortar_longitudinal_strip = {0.0, 0.0};
    int mortar_longitudinal_subdivisions = 1;
    double mortar_resolution = 0.0;
    double weight = 1.0;
    // Conductor-consistency probes (decision 277): the device potential at the conductor
    // vertices of a spatial surface-mortar trace mesh, sampled as the last
    // probe_point_count points of the patch (after the quadrature points and the
    // conductor references; the trace walks stop before them). Their mesh coordinates
    // and the patch's provenance for the record.
    int probe_point_count = 0;
    std::vector<std::array<double, 3>> probe_points;
    int feature = -1;
    double claim_length = 0.0;
    double cell_length = 0.0;
  };

  const FiniteElementSpace &fespace;
  mfem::Array<int> dbc_tdof_list;
  std::vector<ResponseModel> models;
  std::vector<Patch> patches;
  int basis_size;
  int global_basis_size = 0;
  int global_patch_count = 0;
  std::vector<PatchAssignment> patch_assignments;
  int point_query_count = 0;
  std::vector<int> point_send_counts;
  std::vector<int> point_send_offsets;
  std::vector<int> point_receive_counts;
  std::vector<int> point_receive_offsets;
  std::vector<int> point_send_counts_pair;
  std::vector<int> point_send_offsets_pair;
  std::vector<int> point_receive_counts_pair;
  std::vector<int> point_receive_offsets_pair;
  std::vector<int> point_send_indices;
  std::vector<int> point_dof_offsets;
  std::vector<int> point_dofs;
  std::vector<double> point_weights;
  bool maxwell = false;

  // Maxwell postprocessing uses local coupon contour line integrals instead of H1 point
  // values. The line quadrature representation supplies both the trace action and its
  // transpose, so the same map is used for postprocessing and self-consistent correction.
  std::vector<std::vector<mfem::Vector>> maxwell_contours;
  std::vector<std::vector<mfem::Vector>> maxwell_conductor_anchors;
  std::vector<MaxwellLine> maxwell_lines;
  std::vector<MaxwellContourPath> maxwell_paths;
  std::vector<MaxwellConductorPath> maxwell_conductor_paths;
  std::vector<std::pair<int, int>> maxwell_patch_paths;
  int maxwell_quadrature_order = 0;

  double matching_radius = 0.0;
  double mesh_coordinate_scale = 1.0;
  double minimum_wave_speed = mfem::infinity();
  double matched_length_fraction = 1.0;
  double corner_neighborhood_fraction = 0.0;
  std::map<int, double> matched_length_fraction_by_interface;
  std::map<int, double> corner_neighborhood_fraction_by_interface;
  double maximum_curvature_ratio = 0.0;
  double maximum_library_distance = 0.0;
  bool boundary_law_verified = true;

  mutable Vector x_free, local_x, local_x_imag, local_y, trace, response, correction;
  mutable Vector mortar_load, mortar_coefficients;
  mutable std::vector<double> point_owned_values, point_packed_values;
  mutable std::vector<double> point_owned_values_pair, point_packed_values_pair;
  mutable Vector maxwell_point_values, maxwell_conductor_adjoint, maxwell_path_adjoint;

  // Mesh-independent matching statistics are cached with the automatic geometry. The
  // remaining counters describe this mesh-specific interpolation and runtime operator.
  nlohmann::json automatic_statistics;
  // The ownership records of the translational stretches inside the spatial supports
  // (decision 236), reported under Diagnostics with the statistics.
  nlohmann::json ownership_diagnostics;
  // The uncovered requirements (decision 394 F2): the portions of the unmatched features
  // as perimeter sub-segments in mesh coordinates, whose raw within-R surface energy the
  // electrostatic driver keeps in the corrected interface energies.
  std::vector<config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
      uncovered_portions;
  // The placement's clip of those portions by the matched clusters' support boxes
  // (decision 399 MAJOR-1); uncovered_portions holds the clipped portions.
  UncoveredSpatialSupportClipping uncovered_spatial_support_clipping;
  // The raw claims of the DomainBoundary-excluded patches (F-DB-a; DESIGN 2.1).
  DomainBoundaryPortions domain_boundary_portions;
  // The Natural mirror planes and band (mesh units) of the mirror-point evaluation
  // (F-DB-c; DESIGN 2.2.3); empty / 0 without the mirror.
  std::vector<config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData>
      mirror_planes;
  double mirror_band = 0.0;
  double mirror_tolerance = 0.0;
  long long int mirrored_point_count = 0;
  long long int candidate_query_count = 0;
  long long int fallback_query_count = 0;
  long long int point_send_peer_count = 0;
  long long int point_receive_peer_count = 0;
  long long int point_send_item_count = 0;
  long long int point_receive_item_count = 0;
  long long int stencil_nonzero_count = 0;
  long long int contour_line_count = 0;
  mutable long long int operator_mult_count = 0;
  mutable long long int eliminate_rhs_count = 0;
  mutable long long int trace_forward_count = 0;
  mutable long long int trace_transpose_count = 0;

  void ConfigurePointCommunication(
      DistributedPointLocator &locator, const mfem::Vector &xyz, int dimension,
      const std::vector<std::array<double, 3>> *weighted_tangents = nullptr);
  void ConfigureMaxwellLines(const std::vector<MaxwellLineGeometry> &line_geometry);
  void EvaluatePointValues(const Vector &x, Vector &values) const;
  void EvaluatePointValues(const Vector &xr, const Vector &xi, Vector &vr,
                           Vector &vi) const;
  void AddPointValuesTranspose(const Vector &values, Vector &y) const;
  void EvaluatePoints(const Vector &x, Vector &values) const;
  void EvaluateMaxwellLines(const Vector &x, Vector &values) const;
  void EvaluateMaxwellLines(const Vector &xr, const Vector &xi, Vector &vr,
                            Vector &vi) const;
  void AddMaxwellLinesTranspose(const Vector &values, Vector &y) const;
  void BuildMaxwellTrace(const Vector &line_values, Vector &values) const;
  void BuildMaxwellTraceTranspose(const Vector &values, Vector &line_values) const;
  void ApplyTrace(const Vector &x, Vector &values) const;
  void ApplyTraceTranspose(const Vector &values, Vector &y) const;
  // The conductor-consistency probes of a spatial surface-mortar model (decision 277):
  // the conductor vertices on the process plane, off it, and the adjacent free knots; fails
  // closed when a conductor of the trace mesh has no vertex on the plane. A trace mesh
  // without any conductor vertex gets no probes (the constructor records the model under
  // Diagnostics.ConductorConsistency.UnprobedModels and warns).
  static void ConfigureConductorConsistencyProbes(ResponseModel &model);
  // The Diagnostics.ConductorConsistency object with its defaults (tolerance, floor, the
  // counters, the empty lists and the rule), created on first use.
  nlohmann::json &ConductorConsistencyDiagnostics();
  void ApplyUneliminated(const Vector &x, Vector &y) const;
  // y = Pᵀ W D P x on the full potential x (no essential-dof masking): the weighted
  // domain defect of every patch, D = the model's fixed-trace defect when fixed_trace is
  // set, else the defect of its domain-coupling mode (ApplyUneliminated).
  void ApplyDomainDefect(const Vector &x, Vector &y, bool fixed_trace) const;
  void ConfigureMaxwellResponse(
      const IoData &iodata, const MaterialOperator &mat_op,
      const mfem::Array<int> &dbc_tdof_list,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry);

public:
  struct EnergyCorrection
  {
    double domain = 0.0;
    std::map<int, double> interfaces;
  };

  struct PatchTrace
  {
    int patch = 0;
    int model = 0;
    int contour_size = 0;
    std::vector<double> coefficients;
  };

  struct ModelContribution
  {
    int model = 0;
    double patch_count = 0.0;
    double patch_weight = 0.0;
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;
  };

  // The energies of one applied patch (decision 352 follow-up (1)): the increments the
  // patch adds to its model's ModelContribution (the gate-updated weight included), with
  // its provenance (the placed patch index, feature, longitudinal cell in mesh units).
  // Ordered by patch index and replicated on all ranks when requested.
  struct PatchContribution
  {
    int patch = 0;
    int model = 0;
    int feature = -1;
    double weight = 0.0;
    std::array<double, 2> cell = {0.0, 0.0};
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;
  };

  struct ElectrostaticResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;
    std::map<int, double> trace_closure_spread;
    std::vector<ModelContribution> model_contributions;
    std::vector<PatchContribution> patch_contributions;
    double maximum_trace_closure_spread = 0.0;
    double response_weighted_trace_closure_spread = 0.0;
    double trace_closure_response_failure_fraction = 0.0;
    bool confident = true;
  };

  struct MaxwellResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;

    double kR = 0.0;
    double loop_residual = 0.0;
    double response_weighted_loop_residual = 0.0;
    double loop_response_failure_fraction = 0.0;
    double matched_length_fraction = 1.0;
    double corner_neighborhood_fraction = 0.0;
    std::map<int, double> matched_length_fraction_by_interface;
    std::map<int, double> corner_neighborhood_fraction_by_interface;
    double maximum_curvature_ratio = 0.0;
    double maximum_library_distance = 0.0;
    bool boundary_law_verified = true;
    double maximum_trace_closure_spread = 0.0;
    double response_weighted_trace_closure_spread = 0.0;
    double trace_closure_response_failure_fraction = 0.0;
    bool closure_independent_confident = true;
    bool confident = true;
  };

  SurfaceResponseOperator(
      const IoData &iodata, const LaplaceOperator &laplace_op,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry = nullptr);
  SurfaceResponseOperator(
      const IoData &iodata, const SpaceOperator &space_op,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry = nullptr);
  SurfaceResponseOperator(
      const IoData &iodata, const BoundaryModeOperator &mode_op,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry = nullptr);

  void Mult(const Vector &x, Vector &y) const override;
  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }

  // Add the contribution from prescribed essential values to an already assembled
  // thin-metal right-hand side.
  void EliminateRHS(const Vector &x, Vector &rhs) const;

  // Evaluate the nondimensional domain- and surface-energy defects for a global field.
  EnergyCorrection GetEnergyCorrection(const Vector &x) const;

  // The fixed-trace domain defect as a bilinear form on full potentials: y = Pᵀ W
  // (Q_fab,dom - Q_thin,dom) P x with the gate-updated patch weights W and no essential-dof
  // masking (the conductor part of the trace stays, as in GetElectrostaticResponse), every
  // model in its fixed-trace form irrespective of its domain-coupling mode. Its quadratic
  // form is twice the fixed-trace domain correction: 1/2 xᵀ y =
  // GetElectrostaticResponse(x).domain_correction. Collective; x and y are true-dof
  // vectors.
  void FixedTraceDomainDefectMult(const Vector &x, Vector &y) const;

  // Evaluate the complete fabricated-coupon surface energy for every mapped target
  // interface. Corrected participation replaces the measured global core with this data.
  std::map<int, double> GetFabricatedSurfaceEnergy(const Vector &x) const;

  // Evaluate coupon responses on an electrostatic potential. With fixed flux enabled,
  // return both complete postprocessing-only closure ensembles. Otherwise use the active
  // per-model domain-coupling policy for corrected-domain accounting while retaining the
  // fabricated fixed-trace surface evaluation. With include_patches, also return every
  // applied patch's own increments (PatchContribution; their sums per model are the
  // ModelContribution values to roundoff), gathered on all ranks.
  ElectrostaticResponse GetElectrostaticResponse(const Vector &x,
                                                 bool include_fixed_flux = true,
                                                 bool include_patches = false) const;

  // Collect the actual local contour and conductor-state coefficients for every
  // three-dimensional spatial response patch. The result is ordered by global patch
  // index and replicated on all ranks.
  std::vector<PatchTrace> GetSpatialPatchTraces(const Vector &x) const;
  // The same for every applied patch (translational patches included; their coefficients
  // are the surface-mortar projection or the collocated values of the trace coupling).
  std::vector<PatchTrace> GetPatchTraces(const Vector &x) const;

private:
  std::vector<PatchTrace> GatherPatchTraces(const Vector &x, bool spatial_only) const;

public:
  // Conductor-consistency gate (decision 277 (A)), solve time only: measure every applied
  // spatial surface-mortar patch on the potential x of excitation `source`
  // (ConductorConsistencyRecord), exclude the patches whose max_ratio exceeds
  // kConductorConsistencyTolerance (weight 0 from now on, like a DomainBoundary cell: the
  // fixed-trace energies, the self-consistent operator and the ModelCatalog weights of
  // this and every later excitation omit them), append the records to the Diagnostics
  // ("ConductorConsistency": every tested patch of every excitation, the excluded ones
  // with their claim / cell lengths next to the DomainBoundary record) and return the
  // records of this excitation (ordered by global patch index, replicated on all ranks).
  // Collective. The preflight (dry run) cannot evaluate it: it needs the device trace.
  std::vector<ConductorConsistencyRecord> ApplyConductorConsistencyGate(const Vector &x,
                                                                        int source);
  // Whether any patch carries conductor-consistency probes (a spatial surface-mortar
  // model with conductor vertices on the process plane).
  bool HasConductorConsistencyProbes() const;

  // Evaluate a postprocessing-only response for a complex Nedelec Maxwell field. Coupon
  // voltages are reconstructed from transverse contour integrals and applied through
  // Hermitian response-matrix quadratic forms.
  MaxwellResponse GetMaxwellResponse(const GridFunction &E,
                                     std::complex<double> omega) const;

  bool HasSurfaceResponse() const;
  std::set<int> GetTargetInterfaces() const;
  // The uncovered requirements (decision 394 F2; empty when every feature is matched).
  const std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData> &
  GetUncoveredPortions() const
  {
    return uncovered_portions;
  }
  // The parts of the uncovered requirements removed inside the matched clusters' support
  // boxes at placement (decision 399 MAJOR-1; GetUncoveredPortions holds the remainder).
  const UncoveredSpatialSupportClipping &GetUncoveredSpatialSupportClipping() const
  {
    return uncovered_spatial_support_clipping;
  }
  // The raw claims of the DomainBoundary-excluded patches (F-DB-a, decisions 442 / 454;
  // empty when nothing is excluded): kept in the corrected interface energies like the
  // uncovered requirements and reported as the DomainBoundary share.
  const std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData> &
  GetDomainBoundaryPortions() const
  {
    return domain_boundary_portions.portions;
  }
  const DomainBoundaryPortions &GetDomainBoundaryPortionRecord() const
  {
    return domain_boundary_portions;
  }

  int GetBasisSize() const { return global_basis_size; }
  int GetPatchCount() const { return global_patch_count; }
  const std::vector<PatchAssignment> &GetPatchAssignments() const
  {
    return patch_assignments;
  }
  int GetEdgeCount() const { return GetPatchCount(); }
  // Runtime model index (ModelContribution::model) -> model name.
  std::map<int, std::string> GetModelNames() const
  {
    std::map<int, std::string> names;
    for (const auto &model : models)
    {
      names[model.idx] = model.name;
    }
    return names;
  }
  double GetPatchWeight() const;
  double GetMatchingRadius() const { return matching_radius; }

  // Collect deterministic setup/work-distribution counters and runtime application
  // counts. This method is collective over the response operator communicator.
  nlohmann::json GetStatistics() const;

private:
  // Gather the local per-patch response records of GetElectrostaticResponse on every rank
  // and decode them, ordered by patch index.
  std::vector<PatchContribution>
  GatherPatchContributions(const std::vector<double> &local_records) const;
};

// Classify the configured automatic surface-response neighborhoods without assembling a
// finite-element operator or solving a field problem. The deterministic JSON manifest
// reports exact, interpolated, and missing process-library requirements.
void WriteSurfaceResponseRequirements(const IoData &iodata, const Mesh &mesh,
                                      const std::string &path);

// The response-geometry cache (PALACE_RESPONSE_GEOMETRY_CACHE): the automatically
// constructed models and UNPLACED patches with the provenance the continuation ownership
// needs (feature, mesh segment, chain stretch, own-edge offset, the cluster patches'
// claims) and the matching radius; version 6, older caches refused. The reader returns the
// request with its library, models and patches replaced by the cached ones.
// Legacy-contract alias of a library model (USER decision 283): an explicit mapping from a
// contract-3 key (a SpatialEdgeCluster signature carrying Box + Context) to a legacy model
// built under the decision-236 straight-continuation contract, listed by the library under
// the model's "LegacyContractAliases" with the key's context digest. Resolved by the
// matching pass ONLY for the listed key, never as a fallback for any other key.
struct LegacyContractAlias
{
  std::string key;             // the v3 feature hash (64 hex)
  std::string context_digest;  // SpatialSupportContextDigest of the v3 signature (64 hex)
  std::string reason;
  nlohmann::json context;  // the recorded Box + Context (informative)
};

// The alias record of a feature whose hash is the alias key: fails closed (MFEM_ABORT) when
// the feature's context digest differs from the alias's or when the feature's claims-only
// signature is not the legacy model's Signature (the alias names another geometry).
nlohmann::json ResolveLegacyContractAlias(const std::string &model_name,
                                          const nlohmann::json &model_signature,
                                          const LegacyContractAlias &alias,
                                          const IdentifiedFeature &feature);

// Distance from a point to a straight segment a-b.
double SegmentDistance(const std::array<double, 3> &q, const std::array<double, 3> &a,
                       const std::array<double, 3> &b);

// Distance from a point to the arc of the fitted circle (centre, radius rho) that the chord
// a-b subtends (block (b) DESIGN section 2 (a)): the point's in-plane projection inside the
// chord's angular interval reads the distance to the circle (radial residual + out-of-plane
// offset), outside it the distance to the nearer chord end; a chord collinear with the
// centre falls back to the straight distance.
double ArcChordDistance(const std::array<double, 3> &q, const std::array<double, 3> &a,
                        const std::array<double, 3> &b, const std::array<double, 3> &center,
                        double rho);

// Distance from a point to the device perimeter of an identification result: the segment
// keys as straight chords, a segment on a fitted arc (Segments[].Arc) read on its circle's
// arc. The A10 check extended to the context (decision 282 section 3) reads it for every
// placed context piece end (tolerance kSignatureParameterToleranceOverRadius x R).
double DevicePerimeterDistance(const IdentificationResult &identification,
                               const std::array<double, 3> &q);

void WriteResponseGeometryCache(
    const std::filesystem::path &path,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config);
config::ElectrostaticSolverData::ResponseCorrectionData ReadResponseGeometryCache(
    const std::filesystem::path &path,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request);

// Ownership record of the translational patches against the spatial supports (decisions
// 224 / 236): a translational STRETCH — the maximal contiguous stretch of one feature side
// along its chain (patch provenance feature / stretch), the identification's portion unit —
// whose every longitudinal cell lies strictly inside one spatial support's box is RECORDED
// (manifest / metadata Diagnostics.TranslationalStretchesInsideSpatialSupport and a
// warning), never an abort: the coupon is calibrated on the cluster's claims and their
// straight continuations to the box face, so a stretch that continues one of the cluster's
// claims through its claim cut (class Continuation) is corrected by both the coupon and its
// own patches, while any other stretch (class Foreign) is absent from the coupon's twins —
// a model mismatch, not a double count. A stretch continues a claim only through the
// side's OWN edge (decision 252): when one of its cells lies on the claim's mesh segment
// (the claim boundary cut that segment) or when it runs parallel to the claim (within the
// signature angle tolerance), one of its cell ends ON ITS OWN EDGE (the cell ends shifted
// by the provenance edge_offset along AxisU: a pair's cells sit on the midline, a stack's
// on the first side) abuts a claim end along the chain within continuation_tolerance and
// within continuation_tolerance of it transversely (a claim boundary snapped onto a mesh
// vertex of the same edge), AND it extends beyond that claim end (every cell end on the
// outward side of the abutting end, away from the claim's other end, within
// continuation_tolerance; a parallel stretch alongside the claim over the claim's own
// range is Foreign). The side of a pair or stack whose own edge is not claimed never
// continues the claim, whatever its cells' proximity. Judged per stretch, never per cell:
// the box extends 3R past every claim-cut end along its edge (2R continuation + R padding;
// 2R transversely), so the first cells of every stack portion adjacent to a cluster lie
// inside its box legitimately. Cell ends are the strip ends along AxisW from the origin
// (patch units); boxes are (spatial patch index, min, max, claims) in patch units;
// continuation_tolerance is the abutment tolerance along the chain (the signature parameter
// tolerance 1e-3 R). Records in (feature, stretch, spatial patch) order.
struct TranslationalOwnershipRecord
{
  int feature = -1;
  int stretch = -1;
  std::size_t first_patch = 0;
  std::size_t patch_count = 0;
  std::size_t spatial_patch = 0;
  double length = 0.0;               // sum of the stretch's cell lengths
  std::array<double, 3> lo{}, hi{};  // extent of the stretch's cell ends
  bool continuation = false;
};
struct SpatialSupportBounds
{
  std::size_t patch = 0;
  std::array<double, 3> min{}, max{};
  std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionPatchData::Provenance::Claim>
      claims;
  // A contract-3 model's support box and chain pieces in the patch's local frame (units of
  // the matching radius; rule B4), copied from the patch provenance; empty for a legacy
  // model.
  bool has_support_box = false;
  std::array<double, 4> support_box{};
  std::vector<std::array<double, 4>> chain;
  // The bounds come from the Signature's box (a placeholder without basis points).
  bool from_signature_box = false;
};
std::vector<TranslationalOwnershipRecord> FindTranslationalStretchInsideSpatialSupport(
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const std::vector<SpatialSupportBounds> &supports, int dimension,
    double continuation_tolerance);

// The spatial supports of a response configuration: the bounding box of every spatial
// model's basis points (mesh units) placed by its patch frame, with the cluster's claims;
// models without basis points (preflight placeholders) are reported by name in skipped.
std::vector<SpatialSupportBounds> CollectSpatialSupports(
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const std::function<const std::vector<std::array<double, 3>> *(int model_idx)>
        &basis_points,
    double coordinate_scale, int dimension, std::vector<std::string> *skipped = nullptr);

// The Diagnostics entry of the ownership records (lengths and coordinates in mesh units).
struct ContinuationOwnership;
nlohmann::json DescribeTranslationalOwnershipRecords(
    const std::vector<TranslationalOwnershipRecord> &records,
    const std::vector<SpatialSupportBounds> &supports,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale, const std::vector<std::string> &skipped,
    const ContinuationOwnership *ownership = nullptr);

// The warning logged for a non-empty Diagnostics entry (one line per record).
std::string DescribeTranslationalOwnershipWarning(const nlohmann::json &diagnostics);

// Continuation ownership at placement (decisions 236 (2) / 244): a translational cell of a
// stretch that continues a claim of a spatial support (the Continuation criterion of the
// stretch record above, judged on the whole stretch, inside the box or not) is owned by
// that coupon inside the box — the coupon's twins carry the straight continuation of its
// claims to the box face, so the cell's own patch there is the double count. The cell is
// CLIPPED exactly at the box face: its kept part is the part outside EVERY box whose claims
// its stretch continues (the complement intersection: symmetric for a cell on the
// continuations of two coupons, exact, idempotent); the patch keeps weight x kept / cell
// (the weight is linear in the cell length), its origin moves to the kept interval's
// midpoint with the cell re-expressed symmetric about it, its provenance quadrature weight
// scales by the same fraction (the portion [s0, s1) stays, so the dry run's weight formula
// holds and the quadrature weights of a portion sum to 1 - owned / portion length), and
// the Maxwell conductor anchors move with the origin. A cell wholly inside keeps weight 0
// (the operator skips it). Cells of foreign stretches, cells of a continuing stretch
// outside the box and the stack-end cells of a stretch that does not continue a claim are
// untouched. The owned length of a cell shared by two coupons is attributed per coupon by
// the midpoint between the two coupons' inside intervals along the cell (the per-coupon
// lengths sum to the removed length). Fails closed when a box's inside interval lies
// strictly inside a cell (two kept pieces: a cell longer than a coupon box along its own
// direction is no translational cell of a coupon library). Curved cells never continue a
// claim (the Continuation criterion is parallel within the signature angle tolerance), so
// an arc continuing an arc claim keeps its patches: a residual double count the record
// lengths show. Lengths in patch units.
struct ContinuationOwnership
{
  struct Cell
  {
    std::size_t patch = 0;
    int feature = -1;
    int stretch = -1;
    double cell_length = 0.0;
    double owned_length = 0.0;  // removed from the cell (once, whatever the owner count)
    std::vector<std::size_t> owners;  // spatial patches, ascending
    std::vector<double> attributed;   // owned length per owner (midpoint rule)
    // The owned part attributed to each owner as an own-edge sub-segment (the cell line
    // shifted by the provenance edge offset; mesh units, parallel to `owners`): the F-DB-a
    // footprint of an owner that is a DomainBoundary exclusion (DESIGN 2.1).
    std::vector<std::array<std::array<double, 3>, 2>> attributed_intervals;
  };
  std::vector<Cell> cells;  // in patch order
  // Owned length per (feature, stretch, spatial patch) and per spatial patch.
  std::map<std::tuple<int, int, std::size_t>, double> owned_by_stretch;
  std::map<std::size_t, double> owned_by_support;
  double owned_length = 0.0;
  double shared_length = 0.0;  // of cells with two or more owners
  int shared_cells = 0;
  int wholly_owned_cells = 0;
  int clipped_cells = 0;
  // Vertex ownership (decision 282 rule B4, decision 285 (4)): a vertex feature's patch
  // (corner / junction / endpoint coupon: coupon_depth 0, no claims) whose vertex lies on a
  // chain piece END of a contract-3 coupon's continuation chain, inside that coupon's box,
  // is owned by the coupon — the device-plan coupon contains the real corner, so the corner
  // coupon's band on the chain arms would be counted twice. The patch keeps weight 0 (once,
  // whatever the owner count; a vertex inside two boxes lists both owners, review MINOR-1
  // (b)); its face distance is recorded and a vertex closer than R to a face (an arm partly
  // outside the box, MINOR-1 (a)) is flagged. A vertex on another cluster's CLAIMS is never
  // a chain vertex (the chain stops at those claims) and is never owned here.
  struct Vertex
  {
    std::size_t patch = 0;
    int feature = -1;
    std::vector<std::size_t> owners;    // spatial patches, ascending
    double face_distance_over_r = 0.0;  // from the nearest face of the first owner's box
    double chain_end_distance_over_r = 0.0;
    bool arm_outside_box = false;
    // max(0, 1 - face_distance_over_r): the part of the corner's R window beyond the first
    // owner's box (R1 final review MINOR-2; 0 unless arm_outside_box).
    double lost_arm_length_over_r = 0.0;
  };
  std::vector<Vertex> vertices;  // in patch order
  int shared_vertices = 0;
};
// `continuation_tolerance` (patch units) is the cell / claim-end tolerance of the
// translational ownership; `matching_radius` (R, patch units) scales the vertex ownership's
// local frames (passed explicitly: R1 final review MINOR-7).
ContinuationOwnership ApplyContinuationOwnership(
    std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData> &patches,
    const std::vector<SpatialSupportBounds> &supports, int dimension,
    double continuation_tolerance, double matching_radius);

// Clip the uncovered requirements' portions by the matched spatial clusters' support boxes
// (decision 399 MAJOR-1): for every support with claims (a matched SpatialEdgeCluster; a
// vertex coupon's support has none) the part of every portion strictly inside the
// support's bounds — the interval the continuation ownership's CellInsideBox finds on a
// translational cell, the same bounds and the same 1e-12 x max(1, length) tolerance — is
// removed; the kept part (nothing, one piece, or two pieces when the portion crosses the
// box) replaces the portion in place, in portion order. Without such a support the
// portions are left untouched (bitwise). Called at placement on a copy of the
// configuration's portions (the boxes need the models' basis points), never on the cache.
UncoveredSpatialSupportClipping ClipUncoveredPortionsBySpatialSupport(
    std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
        &portions,
    const std::vector<SpatialSupportBounds> &supports, int dimension);

// Coupon-vs-coupon margin overlap (decision 244): two spatial cluster supports whose boxes
// overlap in their interiors are recorded, not aborted, when the overlap is MARGINS ONLY —
// no claim SEGMENT of either enters the other's claims hull (the bounding box of its
// claims, patch units; a coplanar direction widened to the other's box) by a positive
// length beyond continuation_tolerance — tested on the segment, so a claim crossing the
// hull with both ends outside counts (decision 252). Each coupon's twins continue every
// claim CUT end (a claim end no other claim of the same cluster shares within
// continuation_tolerance) straight to its own box face; where that continuation lies on a
// claim of the other coupon (margin-vs-claim) or on a continuation of the other coupon
// (margin-vs-margin) both coupons correct the same edge: a double count the placement
// cannot remove (the dense coupon operator is not clippable), quantified here per pair
// (the continuations are tested parallel within the signature angle tolerance and within
// continuation_tolerance transversely). A claim inside the other's claims hull is a true
// overlap of two coupons' domains: claim_in_hull = true, and the placement aborts.
struct SpatialSupportMarginOverlap
{
  std::size_t first_patch = 0, second_patch = 0;
  std::array<double, 3> overlap_min{}, overlap_max{};
  double first_margin_over_second_claims = 0.0;
  double second_margin_over_first_claims = 0.0;
  double margin_over_margin = 0.0;
  bool claim_in_hull = false;
};
std::vector<SpatialSupportMarginOverlap>
FindSpatialSupportMarginOverlaps(const std::vector<SpatialSupportBounds> &supports,
                                 int dimension, double continuation_tolerance);

// The Diagnostics entries of the continuation ownership and of the margin overlaps
// (lengths and coordinates in mesh units).
nlohmann::json DescribeContinuationOwnership(
    const ContinuationOwnership &ownership,
    const std::vector<SpatialSupportBounds> &supports,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);
nlohmann::json DescribeSpatialSupportMarginOverlaps(
    const std::vector<SpatialSupportMarginOverlap> &overlaps,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);

// The warning logged for a non-empty margin-overlap entry and the summary printed for a
// non-empty continuation-ownership entry.
std::string DescribeSpatialSupportMarginOverlapWarning(const nlohmann::json &diagnostics);
std::string DescribeContinuationOwnershipSummary(const nlohmann::json &diagnostics);

// The Diagnostics entry of the consistent translational mortar (decision 404 D1): per
// model the constrained metal-band vertices of its surface-mortar hat basis and the rule
// that produced them (canonical coupon frame, the library's length unit).
nlohmann::json DescribeConsistentMortar(
    const config::ElectrostaticSolverData::ResponseCorrectionData &config);
// The Diagnostics entries of the corner-arm trim (decision 394 F1) and of the uncovered
// requirements (decision 394 F2) of a response configuration (lengths and coordinates in
// mesh units); both empty-but-complete when nothing was trimmed / nothing is uncovered.
nlohmann::json DescribeCornerArmTrims(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::CornerArmTrimData> &trims,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);
// The trimmed corners whose vertex coupon the placement does not apply (decision 399
// MINOR-7): excluded_reason names the exclusion of a patch index (nullopt: applied). Their
// second arm's [R, s) is then modelled by nothing (a recorded KNOWN LIMIT).
nlohmann::json DescribeCornerArmTrimExcludedCoupons(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::CornerArmTrimData> &trims,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const std::function<std::optional<std::string>(std::size_t patch_idx)> &excluded_reason,
    double coordinate_scale);
// The uncovered entry describes the placed (clipped) portions and, under
// ClippedBySpatialSupport, the parts removed inside the matched clusters' boxes (decision
// 399 MAJOR-1; nullptr when the portions were not placed against boxes, e.g. a 2D request).
nlohmann::json DescribeUncoveredPortions(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
        &portions,
    const UncoveredSpatialSupportClipping *clipping,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);

// Domain-boundary exclusion (decision 258): a coupon's trace coupling is undefined beyond
// the device domain, which a placed coupon reaches wherever a metal edge meets an
// artificial domain cut (a window cut, a chip outline) within ~R |sin theta| of it (theta
// the angle between the edge and the cut's normal: the cross-section at the edge end sticks
// out of the cut by R |sin theta| on one side). A library-placed 3D patch any of whose
// placed coupon points — its model's basis (contour) points and its conductor references,
// at the patch origin cross-section and, for a translational patch with a longitudinal
// cell, at both cell ends moved kSignatureParameterToleranceOverRadius x R inward along
// AxisW (a cut lead's first cell ends exactly on the cut; a strip ending on a domain face
// keeps its sample slices inside) — lies outside the mesh is NOT applied: its weight
// becomes 0 (the operator skips it, the dry run writes it like a wholly owned cell) and it
// is recorded here (the preflight manifest's
// Identification.Diagnostics.DomainBoundaryExclusions and the operator's Diagnostics; the
// inventory total next to Missing). The containment test is the operator's own element
// point locator on every rank over every tested point, the found flags OR-reduced over the
// communicator: the decision is the same for any rank count and locator path (the
// operator's later point location runs on the applied patches only, and any point it cannot
// locate still fails closed naming the patch). The record's nearest outside point is the
// excluded patch's outside point with the smallest distance to an element bounding box
// (exact for an axis-aligned cut; a lower bound of the distance to the mesh otherwise).
// Fail closed (decision 260): a candidate whose first conductor reference at the origin
// section (the metal-edge point, on a mesh face for every correctly placed patch) is not
// located, or none of whose tested points is, is a misplaced or mis-scaled coupon and
// aborts naming the patch (0-based), model and point; when there were candidates at least
// one applied patch must remain. Every patch of `skipped` (already not applied) is left
// alone; `points` are the local model basis points in patch units (origin + points /
// scale); `model_name` names a model for the abort; coordinates and lengths of the result
// in mesh units. Collective over the mesh's communicator.
struct DomainBoundaryExclusion
{
  std::size_t patch = 0;
  int tested_points = 0;
  int outside_points = 0;
  std::array<double, 3> nearest_outside_point{};
  double nearest_distance = 0.0;
};
// A patch classified Mirrored (boundary-cut DESIGN 2.2.3; decisions 442 / 454): every one
// of its placed points outside the mesh reflects through the Natural mirror planes it lies
// beyond (canonical order, within the band) into the mesh, so its trace is taken by even
// extension (mirror-point evaluation in ConfigurePointCommunication) and the patch stays
// applied with its weight.
struct DomainBoundaryMirrored
{
  std::size_t patch = 0;
  int tested_points = 0;
  int outside_points = 0;  // all reflected into the mesh
};
struct DomainBoundaryExclusions
{
  std::vector<DomainBoundaryExclusion> patches;  // ascending patch index
  std::vector<DomainBoundaryMirrored> mirrored;  // ascending patch index
  long long int tested_patches = 0;
  long long int tested_points = 0;
  long long int reflected_points = 0;  // outside points reflected into the mesh
  double wall_time = 0.0;              // of the containment test, seconds
};
// `mirror_planes` (empty: no mirror, every cut-crossing patch is DomainBoundary) and the
// band (mesh units) classify the outside points: a point beyond a Natural plane within the
// band is reflected and located again (the found flags of the reflections OR-reduced like
// the originals); a patch whose every outside point is located after reflection is
// Mirrored, any other patch with an outside point is DomainBoundary (weight 0). A
// metal-edge reference outside the mesh whose reflection is located is a cut-adjacent
// coupon (decision 314's follow-up), never the misplaced-coupon abort.
DomainBoundaryExclusions FindDomainBoundaryExclusions(
    mfem::ParMesh &mesh,
    std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData> &patches,
    const std::function<const std::vector<std::array<double, 3>> *(int model_idx)>
        &basis_points,
    const std::function<bool(int model_idx)> &spatial_basis,
    const std::function<std::string(int model_idx)> &model_name, double coordinate_scale,
    double matching_radius, const std::set<std::size_t> &skipped,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData>
        &mirror_planes = {},
    double mirror_band = 0.0);

// The mirror planes of a response configuration as the mirror module's planes.
std::vector<MirrorPlane> MirrorPlanesOf(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData> &planes);

// Reflect every point outside the mesh's Natural mirror planes into the domain (the even
// extension of the trace; boundary-cut DESIGN 2.2.3): returns the number of reflected
// points; a point beyond a plane that does not mirror or beyond the band is left in place
// (the locator then fails closed naming the patch).
int ReflectPointsIntoDomain(
    mfem::Vector &xyz, int dimension,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData>
        &mirror_planes,
    double mirror_band, double tolerance);

// The Diagnostics entry of the exclusions (per patch: feature, model, cell and portion
// with their lengths, outside-point count, the nearest outside point; totals: the CELL
// length left uncorrected and the portion sum; lengths and coordinates in mesh units) and
// the summary printed for a non-empty entry.
nlohmann::json DescribeDomainBoundaryExclusions(
    const DomainBoundaryExclusions &exclusions,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);
std::string DescribeDomainBoundaryExclusionSummary(const nlohmann::json &diagnostics);

// F-DB-a (decisions 442 / 454, DESIGN 2.1): NEVER DROP. The raw claim of every
// DomainBoundary-excluded patch as perimeter portions whose within-R raw energy the
// electrostatic driver keeps in the corrected interface energies exactly as the uncovered
// requirements' (decision 394 F2), reported per type as the DomainBoundary share. The claim
// of a translational cell is its OWN-EDGE interval (the clipped cell shifted by the
// provenance edge offset along AxisU: a pair's two sides and a stack's n sides are n
// intervals); the co-located first-order split patches of one cell (model weights summing
// to 1) map to ONE interval, deduplicated within the signature tolerance and counted once.
// The claim of a vertex coupon is its provenance raw_claims (the arms' R windows and a
// trimmed second arm's [R, s)); the claim of a spatial cluster coupon is its claims plus
// the parts its box removed from others: the continuation-owned cell parts attributed to it
// and the uncovered portions clipped by its box. Overlaps are excluded by construction (a
// DB cell is a kept part outside every matched box; a vertex claim ends where the arm cells
// start). Portions in mesh units; `types` names each portion's type (the model topology).
DomainBoundaryPortions CollectDomainBoundaryPortions(
    const DomainBoundaryExclusions &exclusions,
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const ContinuationOwnership &ownership,
    const UncoveredSpatialSupportClipping &uncovered_clipping, double matching_radius);

// The Diagnostics entry of the raw portions (count, geometric cells, duplicates, length,
// per type and per feature, the portion list; mesh units scaled by coordinate_scale).
nlohmann::json DescribeDomainBoundaryPortions(const DomainBoundaryPortions &portions,
                                              double coordinate_scale);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
