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
#include <unordered_map>
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

// The raw claims of EXCLUDED patches as perimeter portions, whose within-R raw energy the
// drivers keep in the corrected interface energies: the DomainBoundary-excluded patches of
// the placement (F-DB-a, decisions 442 / 454, DESIGN 2.1; CollectDomainBoundaryPortions)
// and the conductor-consistency-excluded patches of the solve (S7diag R1, decision 758;
// SurfaceResponseOperator::GetConductorConsistencyPortionRecord). One record shape, one
// per-patch claim rule (AppendExcludedPatchClaim), one share per consumer.
struct ExcludedPatchPortions
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

// OI-SCP-HOLE (decisions 915 / 918 / 922): the corner / vertex coupons the
// SpatialClusterPriority rule drops whole (a spatial-basis non-cluster patch whose basis-
// point box interior-overlaps an exact spatial cluster's), each with its owner, and the
// parts of their raw claims (the arms' R windows) OUTSIDE every matched cluster's box as
// perimeter portions the drivers keep RAW in the corrected interface energies (the part
// inside the box is modelled by the owner coupon). The record is the 915 MINOR-7
// Diagnostics entry; it is written only when at least one coupon is dropped.
struct SpatialClusterPriorityExclusions
{
  struct Entry
  {
    std::size_t patch = 0;        // the dropped (subordinate) patch
    std::size_t owner = 0;        // the cluster patch owning the overlap
    double window_length = 0.0;   // the dropped patch's raw claims, mesh units
    double outside_length = 0.0;  // the part outside every matched cluster's box
  };
  std::vector<Entry> entries;
  ExcludedPatchPortions portions;
};

// PatchEnergy (decision 922 B-1; the 851 (C) diagnosis export formalised, scratch
// 3232e29548): the raw claims of every APPLIED patch (weight > 0, not
// DomainBoundary-excluded), built by the DomainBoundary rule (AppendExcludedPatchClaim),
// labelled by the SET of models whose patches cover the claim ("m13" alone, "m13+m14" for a
// co-located first-order split cell). The surface postoperator partitions the within-R raw
// energy over these claims (GetInterfacePartitionedEdgeEnergies): per label, the raw energy
// on exactly the cells those models correct, i.e. r_raw of decision 850 (d) "on the model's
// own patches".
struct ModelRawClaimPortions
{
  // The claims; `topology` carries the model-set label.
  std::vector<config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
      portions;
  // Per label and model of that label: the patches of the model on the label's claims and
  // the sum over those claims of (patch weight x claim length, mesh units), the share the
  // diagnosis attributes a shared label's raw energy with.
  struct GroupMember
  {
    std::string label;
    int model = 0;
    int patches = 0;
    double weight_length = 0.0;
  };
  std::vector<GroupMember> members;  // ascending (label, model)
  double length = 0.0;               // mesh units, over all portions
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
    // The patch's conductor probes read on the device chord polyline (decision 611 (2),
    // ResponseCorrectionPatchData::Provenance::ProbeSnap): the snapped knots in the order
    // of the placed patch's record; max_probe_snap = the largest displacement
    // (nondimensional; the record multiplies by the mesh coordinate scale).
    struct ProbeSnap
    {
      int vertex = -1;
      int segment = -1;
      int arc = -1;
      double displacement = 0.0;
    };
    std::vector<ProbeSnap> probe_snaps;
    double max_probe_snap = 0.0;
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

  // One chain of line integrals through consecutive trace coefficients (trace_indices; a
  // PEC knot is in no chain: its trace is the zero initialisation). The first value is the
  // start trace (a conductor state, or the anchored contour knot a cap-interior hat hangs
  // from; 0 when there is none) minus the anchor line (from the patch's conductor anchor
  // or from the start knot); every next value subtracts the lines of the step to the next
  // knot. A step between two consecutive ring knots follows the trace mesh's ring polyline
  // (decision 818 R3): one straight line per leg, through the slave vertices (box corners)
  // of that edge, never a chord across a box corner; step_line_counts holds the number of
  // legs of every step between consecutive trace knots, and the contour lines after the
  // last trace step close the loop (a closed ring's step back to its start knot; a free
  // arc's step to the next PEC knot). A closed ring has no anchor start trace; an open path
  // ends with end_line to the end conductor's anchor; a free arc of a ZeroTraceIndices
  // ring ends with return_line from its PEC knot back to the anchor (a true loop for the
  // residual, decision 818 R1 / R4); a hat path has one knot and hat_check_line, the
  // diagnostic line from the conductor anchor straight to the hat (the hat residual, R2).
  struct MaxwellContourPath
  {
    int anchor_line = -1;
    int contour_line_offset = 0;
    int contour_line_count = 0;
    int end_line = -1;
    int return_line = -1;
    int hat_check_line = -1;
    int start_trace = -1;
    int end_conductor_trace = -1;
    bool closed = true;
    std::vector<int> trace_indices;
    std::vector<int> step_line_counts;
  };

  struct MaxwellConductorPath
  {
    int contour_path = -1;
    int parent_trace_offset = -1;
    int trace_offset = 0;
    int integral_sign = 0;
  };

  // OI-MX-LOOPGATE (decisions 950 / 952 / 954 / 955): the per-patch weighted least-squares
  // lift of the Maxwell trace. The chain above sets every knot by ONE route, so on a patch
  // whose field carries a circulation (an inductive patch) the trace is route-dependent.
  // The fit treats every line the paths evaluate as an EQUATION phi(to) - phi(from) = -L
  // between two NODES (a trace coefficient, or the ground super-node -1 = the reference
  // anchor and every PEC knot), weighted by 1 / length (the polyline L2 projection of the
  // tangential field), and solves the normal equations in the CYCLE SPACE: the chain is the
  // exact solution of the TREE equations (the lines it uses), the CLOSING equations (a
  // ring's last step, a free arc's step onto its PEC knot and its return line, a non-tree
  // open path's end line; the hat check lines are information only, decision 818 R2) carry
  // the loop misclosures r; the loop currents I = G^-1 r with G = D_C + K D_T K^T (K = the
  // fundamental cycles' +-1 incidence on the tree steps, D = the step lengths) give the
  // correction u_j = -l_j sum_c K_cj I_c on every tree step, accumulated along the tree.
  // When every loop closes (r = 0) the fit is the chain to roundoff. The W-norm of the
  // circulating remainder is r^T G^-1 r; its share of the patch's route norm sum_j L_j^2 /
  // l_j is the patch's circulating share s_p (GetMaxwellResponse).
  struct MaxwellLoopFit
  {
    // A tree step in topological order (parent set before child): the chain's route into
    // `child` (a trace index; parent -1 = ground), its legs [line_begin, line_begin +
    // line_count) and their total length (0 for an identification by a coincident line).
    struct TreeStep
    {
      int parent = -1;
      int child = -1;
      int line_begin = 0;
      int line_count = 0;
      double length = 0.0;
    };
    // A closing equation phi(to) - phi(from) = -L (from / to as the legs run; -1 = ground)
    // with its fundamental cycle: the tree steps on the tree path from `from` to `to` with
    // the sign +1 where the path runs parent -> child.
    struct ClosingStep
    {
      int from = -1;
      int to = -1;
      int line_begin = 0;
      int line_count = 0;
      double length = 0.0;
      std::vector<std::pair<int, double>> cycle;
    };
    std::vector<TreeStep> tree;
    std::vector<ClosingStep> closing;
    mfem::DenseMatrix cycle_inverse;  // G^-1, closing.size() square
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
    // The rim knot each cap-interior hat hangs from on the Maxwell path (decision 818 R2):
    // the nearest contour knot of a path (never a PEC knot) in the hat's own cap plane,
    // chosen from the model data on every rank (a hat without one is refused at
    // configuration). "Own cap plane" = the same local w within
    // kSignatureParameterToleranceOverRadius x the library matching radius, which is 0 for
    // an explicit (non-library) model: there the knot's w must equal the hat's exactly.
    std::vector<int> hat_rim_knots;
    // The trace mesh's slave vertices (the corner family's box corners, between two
    // consecutive ring knots): the Maxwell contour line of that ring edge runs knot ->
    // slave -> knot instead of the chord (decision 818 R3). Parents are 0-based basis
    // indices; the point is in the model's local coupon coordinates.
    struct RingSlaveVertex
    {
      int parent_a = 0;
      int parent_b = 0;
      double weight_a = 0.0;
      std::array<double, 3> point{};
    };
    std::vector<RingSlaveVertex> ring_slave_vertices;
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
    // Factor on the domain correction only (ResponseCorrectionPatchData::domain_weight): 0
    // for a ForeignExcludedEdge cell (decision 553), whose surface matrices keep `weight`.
    double domain_weight = 1.0;
    // Conductor-consistency probes (decision 277): the device potential at the conductor
    // vertices of a spatial surface-mortar trace mesh, sampled as the last
    // probe_point_count points of the patch (after the quadrature points and the
    // conductor references; the trace walks stop before them). Their mesh coordinates
    // and the patch's provenance for the record.
    int probe_point_count = 0;
    std::vector<std::array<double, 3>> probe_points;
    // The probes read on the device arc chords (decision 611 (2)): the placed patch's
    // ProbeSnaps applied to probe_points, for the record.
    std::vector<ConductorConsistencyRecord::ProbeSnap> probe_snaps;
    double max_probe_snap = 0.0;
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
  // OI-MX-LOOPGATE: per local patch the least-squares structure, the line lengths (mesh
  // units) and, keyed by the placed patch index on every rank, every applied patch's raw
  // claim (the DomainBoundary rule, AppendExcludedPatchClaim) for the remainder energy.
  std::vector<MaxwellLoopFit> maxwell_loop_fits;
  std::vector<double> maxwell_line_lengths;
  std::map<int, ExcludedPatchPortions> maxwell_patch_raw_claims;
  long long int maxwell_loop_fit_cycle_count = 0;
  long long int maxwell_loop_fit_unconstrained_count = 0;

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
  mutable Vector maxwell_point_values, maxwell_trace_adjoint, maxwell_path_adjoint;
  mutable Vector maxwell_fit_adjoint, maxwell_fit_line_adjoint;
  mutable std::vector<double> maxwell_loop_workspace;

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
  ExcludedPatchPortions domain_boundary_portions;
  // OI-SCP-HOLE (decisions 915 / 918 / 922): the corner / vertex coupons dropped whole by
  // the SpatialClusterPriority rule, with the parts of their raw claims OUTSIDE every
  // matched cluster's box as the drivers' fourth raw class (counted, uncorrected).
  SpatialClusterPriorityExclusions spatial_cluster_priority;
  // S7diag R1 (decision 758): the raw claim of every APPLIED patch the
  // conductor-consistency gate can probe (its model has conductor vertices on the process
  // plane), built at construction by the DomainBoundary rule (AppendExcludedPatchClaim) and
  // keyed by the global patch index on every rank; the gate appends an EXCLUDED patch's
  // claim to conductor_consistency_portions (replicated: the records are gathered first),
  // which the electrostatic driver keeps raw in E_ft / E_ff / E_sc like the uncovered
  // requirements and the DomainBoundary claims. The list only grows (an exclusion is
  // sticky).
  std::map<int, ExcludedPatchPortions> conductor_consistency_claims;
  ExcludedPatchPortions conductor_consistency_portions;
  // PatchEnergy (decision 922 B-1): the applied patches' raw claims by model set; empty
  // unless the configuration requests PatchEnergy.
  ModelRawClaimPortions model_raw_claim_portions;
  // kSignatureParameterToleranceOverRadius x R (mesh units): the claim portions' length
  // floor and deduplication tolerance (the DomainBoundary rule's).
  double conductor_consistency_claim_tolerance = 0.0;
  // The domain-boundary containment result was reused from the first construction on the
  // cached automatic geometry rather than recomputed on this mesh (decision 734 F2).
  bool domain_boundary_containment_reused = false;
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
      const std::vector<std::array<double, 3>> *weighted_tangents = nullptr,
      const std::vector<double> *point_band = nullptr);
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
  // The chain lift (one route per knot; the loop-residual detector's trace, decision 955
  // D3) and its transpose.
  void BuildMaxwellChainTrace(const Vector &line_values, Vector &values) const;
  void BuildMaxwellChainTraceTranspose(const Vector &values, Vector &line_values) const;
  // The least-squares lift (the chain plus the cycle-space correction of every patch) and
  // its transpose: the trace of the energies and of the self-consistent Mult.
  void BuildMaxwellTrace(const Vector &line_values, Vector &values) const;
  void BuildMaxwellTraceTranspose(const Vector &values, Vector &line_values) const;
  // One patch's cycle-space correction added to its chain values in `values` (a real
  // line-value vector); returns the circulating W-norm r^T G^-1 r of that vector.
  double AddMaxwellLoopFitCorrection(std::size_t patch_index, const Vector &line_values,
                                     Vector &values) const;
  // One patch's route norm sum_j L_j^2 / l_j over its equations with a positive length.
  double MaxwellLoopFitRouteNorm(std::size_t patch_index, const Vector &line_values) const;
  // Builds maxwell_loop_fits from the paths (ConfigureMaxwellResponse, after the paths).
  void ConfigureMaxwellLoopFit();
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
  // The placement-ownership stage shared by the electrostatic and Maxwell constructors
  // (decision 818 R5): the configured patches with the continuation ownership applied
  // (placed_patches) and the set of patch indices that are not applied (owned by a
  // spatial cluster, wholly trimmed, vertex-owned or DomainBoundary-excluded). Records the
  // uncovered / DomainBoundary / conductor-consistency portions and the ownership
  // Diagnostics as members. Collective (the domain-boundary containment test).
  struct PlacementOwnership
  {
    std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        placed_patches;
    std::set<std::size_t> owned_patches;
  };
  PlacementOwnership ApplyPlacementOwnership(
      const config::ElectrostaticSolverData::ResponseCorrectionData &config,
      const std::unordered_map<int, int> &model_indices,
      const std::vector<std::vector<std::array<double, 3>>> &basis_points,
      double coordinate_scale, int dimension,
      const std::shared_ptr<const SurfaceResponseGeometry> &cached_geometry,
      bool setup_timing);

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
    // Cap-interior hats (decision 818 R2): the maximum and the response-weighted RMS of
    // |V_hat(rim route) - V_hat(anchor route)| / scale over the hat paths; information,
    // not a gate (decision 814 (1)).
    double hat_residual = 0.0;
    double response_weighted_hat_residual = 0.0;
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
    // Per model (decision 845 O-C, M-5): the same per-patch increments as the totals above,
    // summed over the model's applied patches (the electrostatic driver's aggregation of
    // surface-response-model-energy.csv), one entry per runtime model in model-index order,
    // gathered by the one batched reduction and replicated on all ranks.
    std::vector<ModelContribution> model_contributions;
    // PatchEnergy on the Maxwell path (decisions 913 / 920 / 922 B; include_patches): every
    // applied patch's own increments (the electrostatic PatchContribution form: the
    // per-model sums are the entries above to roundoff; `cell` = the mortar longitudinal
    // strip), the spatial patches' complex trace coefficients (the electrostatic traces
    // file's form, real and imaginary parts), and the per-path loop residuals of the
    // line-integral trace lift (920: one record per contour loop (kind 0) or cap-interior
    // hat loop (kind 1) with a non-zero scale, `residual` = |Σ terms| / Σ |terms|, `scale`
    // = Σ |terms| in the field's units; the maximum over the records is `loop_residual`,
    // the quantity behind the 0.05 gate). All three ordered by global patch index and
    // replicated on all ranks; empty without include_patches.
    std::vector<PatchContribution> patch_contributions;
    struct SpatialTrace
    {
      int patch = 0;
      int model = 0;
      int contour_size = 0;
      std::vector<std::complex<double>> coefficients;
    };
    std::vector<SpatialTrace> spatial_traces;
    struct PathLoopResidual
    {
      int patch = 0;
      int model = 0;
      int path = 0;
      int kind = 0;
      double residual = 0.0;
      double scale = 0.0;
    };
    std::vector<PathLoopResidual> path_loop_residuals;
    // OI-MX-LOOPGATE (decision 955): per patch the circulating share s_p = ||r||_W^2 /
    // ||L||_W^2 of its line integrals (0 for a curl-free field), the detector weight (the
    // patch's response energy on the CHAIN trace, the loop-RMS weight) and the fallback
    // flag; every patch when include_patches, ordered by global patch index, replicated.
    struct PatchLoopShare
    {
      int patch = 0;
      double share = 0.0;
      double detector_energy = 0.0;
      int fallback = 0;
    };
    std::vector<PatchLoopShare> patch_loop_shares;
    // The patches whose circulating remainder is counted RAW: rho_p = share above the
    // floor kLoopFitShareFloor (fallback 0), or rho_p = 1 above kLoopFitFallbackShare
    // (fallback 1: the patch is served raw, its fabricated energies and domain defects are
    // not added). Always gathered, ordered by global patch index, replicated.
    struct LoopRemainder
    {
      int patch = 0;
      int model = 0;
      double share = 0.0;
      double circulation = 0.0;  // r^T G^-1 r, real + imaginary parts (V^2 / mesh unit)
      double total = 0.0;        // sum_j |L_j|^2 / l_j over the patch's equations
      int fallback = 0;
    };
    std::vector<LoopRemainder> loop_remainders;
    double maximum_loop_fit_share = 0.0;
    double response_weighted_loop_fit_share = 0.0;
    int loop_remainder_patch_count = 0;
    int loop_fallback_patch_count = 0;
  };
  // OI-MX-LOOPGATE (decision 955 D5 / D6, pre-registered): a patch whose circulating share
  // is at most the floor contributes no remainder (the uncounted part is at most the floor
  // times its within-R energy); above the fallback share its potential is not worth
  // correcting and the whole patch is served raw.
  static constexpr double kLoopFitShareFloor = 1.0e-4;
  static constexpr double kLoopFitFallbackShare = 0.25;
  // The raw claims of the patches of `remainders` as perimeter portions labelled
  // "m<model>" with the per-portion weight rho_p (the surface postoperator's weighted
  // partitioned integral), and per label the patch count and the rho-weighted claim length
  // (mesh units).
  struct LoopRemainderPortions
  {
    std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData>
        portions;
    std::vector<double> weights;
    struct Row
    {
      int patches = 0;
      double weighted_length = 0.0;
    };
    std::map<std::string, Row> rows;
  };
  LoopRemainderPortions GetLoopRemainderPortions(
      const std::vector<MaxwellResponse::LoopRemainder> &remainders) const;

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
  // are the surface-mortar projection or the collocated values of the trace coupling). On
  // the Maxwell Collocated path the coefficients are the trace ApplyTrace builds from the
  // real vector x (the contour line integrals of that one field); no production Maxwell
  // caller, a unit-test read (decision 821).
  std::vector<PatchTrace> GetPatchTraces(const Vector &x) const;

private:
  std::vector<PatchTrace> GatherPatchTraces(const Vector &x, bool spatial_only) const;
  // Every rank's local records concatenated in rank order (Allgatherv), the per-patch
  // record gathers' shared step.
  std::vector<double> AllgatherRecords(const std::vector<double> &local_records) const;

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
  MaxwellResponse GetMaxwellResponse(const GridFunction &E, std::complex<double> omega,
                                     bool include_patches = false) const;

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
  const ExcludedPatchPortions &GetDomainBoundaryPortionRecord() const
  {
    return domain_boundary_portions;
  }
  // The parts of the SpatialClusterPriority-dropped coupons' raw claims outside every
  // matched cluster's box (OI-SCP-HOLE, decision 922; empty when no coupon is dropped or
  // every dropped window lies inside its owner's box): kept in the corrected interface
  // energies like the uncovered requirements and reported as the SpatialClusterPriority
  // share. The inside part is the owner coupon's; the corner's own correction outside the
  // box stays missing (a box growth of the cluster coupon is a library item).
  const std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData> &
  GetSpatialClusterPriorityPortions() const
  {
    return spatial_cluster_priority.portions.portions;
  }
  const SpatialClusterPriorityExclusions &GetSpatialClusterPriorityExclusions() const
  {
    return spatial_cluster_priority;
  }
  // The raw claims of the patches the conductor-consistency gate has excluded so far
  // (S7diag R1, decision 758; empty before the first exclusion and on every window with 0
  // exclusions): kept in the corrected interface energies like the uncovered requirements
  // and reported as the ConductorConsistency share. The gate runs on the electrostatic path
  // only (ApplyConductorConsistencyGate refuses a Maxwell operator), so a Maxwell
  // operator's list is empty by construction: the Maxwell path consumes no such term today
  // (decision 763); when it gains the gate, it consumes this list as the electrostatic
  // driver does.
  const std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData> &
  GetConductorConsistencyPortions() const
  {
    return conductor_consistency_portions.portions;
  }
  const ExcludedPatchPortions &GetConductorConsistencyPortionRecord() const
  {
    return conductor_consistency_portions;
  }
  // PatchEnergy (decision 922 B-1): the applied patches' raw claims labelled by model set
  // (empty without PatchEnergy).
  const ModelRawClaimPortions &GetModelRawClaimPortions() const
  {
    return model_raw_claim_portions;
  }
  // Whether this construction reused the first construction's domain-boundary containment
  // result through the cached automatic geometry (decision 734 F2; false on the first
  // construction, on a curved mesh, under PALACE_RESPONSE_CONTAINMENT_RECHECK, or when the
  // tested patch set changed).
  bool DomainBoundaryContainmentReused() const
  {
    return domain_boundary_containment_reused;
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
  // Runtime model index -> (basis size, number of zero-trace (PEC) knots): the constructed
  // basis of an interpolated corner is checked against its rule by the unit tests.
  std::map<int, std::pair<int, int>> GetModelBasisSizes() const
  {
    std::map<int, std::pair<int, int>> sizes;
    for (const auto &model : models)
    {
      sizes[model.idx] = {model.basis_size,
                          static_cast<int>(model.zero_trace_indices.size())};
    }
    return sizes;
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

// The role-insensitive Box + Context digest (decision 899 B2, OI-ROLEKEY): SHA-256 of
// {Box, Context} with the Conductor label removed from every Chain: false piece (a piece
// not connected to the claims inside the box, whose label depends on the conductor roles
// outside it: a claw shorted to ground reads the claims' label where the cut-short claw's
// ground reads a new one). Empty without a Box. The matching pass resolves an otherwise
// unmatched SpatialEdgeCluster through a claims-only legacy model's alias whose stored
// Context has the feature's role-insensitive digest and whose model key is the feature's
// claims-only key (ResolveRoleInsensitiveLegacyContractAlias); the exact alias path is
// untouched and tried first.
std::string RoleInsensitiveContextDigest(const nlohmann::json &signature);
nlohmann::json ResolveRoleInsensitiveLegacyContractAlias(
    const std::string &model_name, const nlohmann::json &model_signature,
    const LegacyContractAlias &alias, const IdentifiedFeature &feature);

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

// The same distance on the MIRROR-COMPLETE plan (decisions 732 / 733, F-M1): a point beyond
// an acting Natural mirror plane (within `band` of it; outside the domain, where the plan
// is the reflected real plan) is reflected back into the domain through the planes it lies
// beyond (ReflectIntoDomain) and read there as well; the smaller of the two distances. The
// A10 check extended to the context reads it with the support box span as the band: a
// reflected continuation piece ends on a box face beyond the mirror band, on no image
// segment. Without planes, or for a point inside the domain, exactly
// DevicePerimeterDistance.
double MirrorCompletePerimeterDistance(const IdentificationResult &identification,
                                       const std::array<double, 3> &q,
                                       const std::vector<MirrorPlane> &planes, double band,
                                       double tolerance);

// The conductor-consistency probe snap (decision 611 (2); the R1 conductor diagnosis
// section 6): the nearest REAL device perimeter segment to q (every segment read as the
// straight chord of its key, the image segments of a mirror band excluded; a segment of
// ANY conductor of the device - the identified segments carry no conductor, and a segment
// of another conductor within the bound below would need a device gap smaller than the
// arc's sagitta, so the patch's own conductor is the one found in practice). When that
// segment is a chord of a fitted arc (Segments[].Arc) and q's in-plane distance to the
// chord - in the arc's plane, the out-of-plane offset set aside - lies within the arc's
// recorded largest chord sagitta (MaxChordSagittaOverR x R) plus the arc-fit tolerance
// (kArcFitToleranceOverRadius x R, the band the device vertices lie on the fitted circle
// within), the result is the displacement from q to its in-plane foot on the chord (the
// out-of-plane offset unchanged) with the segment and arc; otherwise (a straight segment
// nearest, an arc chord farther than its sagitta - coupon metal the device polyline does
// not have - or a probe already on the chord at the signature's resolution,
// kSignatureLengthQuantumOverRadius x R) nullopt - the probe is not moved. A coupon knot
// built ON the fitted circle lies within that bound of the device's chord by construction,
// a knot claiming metal farther away is read where it is.
struct ProbeChordSnap
{
  int segment = -1;
  int arc = -1;
  std::array<double, 3> displacement{};
  double distance = 0.0;
};
std::optional<ProbeChordSnap>
SnapProbeToDeviceChord(const IdentificationResult &identification,
                       const std::array<double, 3> &q);

// The probe snap under the mirror rule (decision 774; the GROUP B 79c3 diagnosis of
// decision 771): the solve evaluates a probe beyond a Natural truncation plane at its
// reflection into the domain (ReflectIntoDomain: the even extension, the planes in index
// order, within `band` = the global mirror band widened to a mirror-formed coupon's reach),
// so the device chord such a probe must be read on is the one at its REFLECTED position -
// the mapped position of an image-side knot lies beyond the plane, where no real segment
// is. The snap is computed at the reflected position and its displacement carried back
// through the same planes in reverse order (ReflectDirection, the linear part of each
// reflection, its own inverse), so that the stored displacement is applied to the MAPPED
// probe exactly as a real-side snap's is and the solve's reflection of (probe +
// displacement) lands on the chord. A probe inside every plane, or one the solve would not
// reflect (beyond the band, a plane that does not mirror, no planes, band 0), is read by
// SnapProbeToDeviceChord at q, bitwise as before.
std::optional<ProbeChordSnap> SnapProbeToDeviceChordAcrossMirror(
    const IdentificationResult &identification, const std::array<double, 3> &q,
    const std::vector<MirrorPlane> &planes, double band, double tolerance);

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
// own patches. Any other stretch is Foreign to the claims; under spatial-support contract
// v3 (Box + Context, decisions 281 / 282 / 285 / 286) the coupon's twins CONTAIN the
// context metal inside the box, so a Foreign stretch lying on a context piece of the coupon
// is a double count too (decision 553, the rerun-2 diagnosis review): on a Chain: false
// piece (the library's ForeignEdges = the qualification's EdgeExcludeSegments, rule B5) the
// coupon's within-R SURFACES leave that edge to the cells but its DOMAIN matrices integrate
// the whole box (sub-class ForeignExcludedEdge: the cells keep their surface matrices and
// drop their domain correction); on a Chain: true piece, or any context piece of a coupon
// without ForeignEdges, surface AND domain are both in the coupon (sub-class
// ForeignContextEdge: the cells are owned by the coupon, weight 0, as a continuation). A
// Foreign stretch on no context piece of a legacy (claims-only) coupon is absent from the
// coupon's twins — a model mismatch, not a double count (the decision-236 premise, true for
// legacy coupons only). A stretch continues a claim only through the
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
// The foreign-context class of a Foreign stretch inside a contract-3 support (the B5 test
// at placement: a CELL lies on a context piece when its OWN-edge ends and midpoint lie
// within continuation_tolerance of one piece of that class, the EdgeExcludeSegments form;
// decision 559 (1): the treatment is PER CELL by piece class, the stretch class summarises
// its cells): NONE (not Foreign, or no cell on any context piece: the legacy class, every
// cell untouched), EXCLUDED_EDGE (every cell on a Chain: false piece: every cell keeps its
// surface matrices and drops its domain correction), CONTEXT_EDGE (every cell on a Chain:
// true piece: every cell owned by the coupon), MIXED_EDGES (cells of different classes:
// each treated by its own class, a cell on no piece untouched; the three counts recorded).
enum class ForeignContextClass
{
  NONE,
  EXCLUDED_EDGE,
  CONTEXT_EDGE,
  MIXED_EDGES
};
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
  ForeignContextClass foreign_context = ForeignContextClass::NONE;
  // Cells of the stretch whose own-edge ends and midpoint lie on a Chain: false / Chain:
  // true context piece of the support (the per-cell classes the treatment applies; a cell
  // on no piece counts in neither: patch_count - both).
  std::size_t cells_on_excluded_edges = 0;
  std::size_t cells_on_chain_edges = 0;
};
struct SpatialSupportBounds
{
  std::size_t patch = 0;
  // The mesh-axis-aligned bounding box of the placed basis points (patch units): read by
  // the margin-overlap diagnostics and the stretch records' Box entry. NOT an ownership
  // test: for a rotated patch frame it circumscribes the coupon's box (decision 715).
  std::array<double, 3> min{}, max{};
  // The support box in the cluster's PATCH FRAME (the canonical frame of contract v3):
  // `axes[k]` = the spatial patch's axis_u / axis_v / axis_w in mesh coordinates, and
  // `frame_min[k]` / `frame_max[k]` the extremes over the placed basis points P of the
  // projection dot(P, axes[k]) (patch units; NOT relative to the frame origin, so that a
  // mesh-aligned frame — every axis a signed unit vector — reproduces `min` / `max` exactly
  // and the ownership arithmetic stays bitwise). LineInsideBox / StretchStrictlyInside test
  // against THIS box. Identity axes by default: a hand-built support sets both boxes at
  // once with SetMeshAlignedBox.
  std::array<std::array<double, 3>, 3> axes = {
      {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  std::array<double, 3> frame_min{}, frame_max{};
  void SetMeshAlignedBox(const std::array<double, 3> &lo, const std::array<double, 3> &hi)
  {
    min = frame_min = lo;
    max = frame_max = hi;
    axes = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  }
  std::vector<
      config::ElectrostaticSolverData::ResponseCorrectionPatchData::Provenance::Claim>
      claims;
  // A contract-3 model's support box, chain pieces (Chain: true) and foreign pieces
  // (Chain: false = ForeignEdges, rule B5) in the patch's local frame (units of the
  // matching radius; rules B4 / decision 553), copied from the patch provenance; empty for
  // a legacy model.
  bool has_support_box = false;
  std::array<double, 4> support_box{};
  std::vector<std::array<double, 4>> chain;
  std::vector<std::array<double, 4>> foreign;
  // The matching radius R in patch units (the scale of `support_box` / `chain` / `foreign`;
  // CollectSpatialSupports copies the configuration's); 0 for a legacy model. REQUIRED
  // (> 0) for a contract-3 support that carries context pieces: the foreign-context
  // classification fails closed (MFEM_VERIFY) without it, so a hand-built support of a unit
  // fixture with `chain` / `foreign` pieces must set it (the B4 vertex ownership scales by
  // the explicit `matching_radius` argument instead).
  double matching_radius = 0.0;
  // The bounds come from the Signature's box (a placeholder without basis points).
  bool from_signature_box = false;
};
std::vector<TranslationalOwnershipRecord> FindTranslationalStretchInsideSpatialSupport(
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const std::vector<SpatialSupportBounds> &supports, int dimension,
    double continuation_tolerance);

// The Diagnostics "Class" of a record: Continuation, Foreign, ForeignExcludedEdge or
// ForeignContextEdge.
const char *TranslationalOwnershipClassName(const TranslationalOwnershipRecord &record);

// The spatial supports of a response configuration: for every spatial model's basis points
// (mesh units) placed by its patch frame, the mesh AABB (`min` / `max`) and the box in the
// patch frame (`axes`, `frame_min` / `frame_max`: the coupon's own box, decision 715), with
// the cluster's claims; a contract-3 placeholder without basis points gets the Signature's
// support box placed by the patch frame (its corners `support_box` x R, built in mesh units
// like the basis points, so the bounds come out in the same patch units as every other
// support); models without basis points and without a support box (legacy preflight
// placeholders) are reported by name in skipped.
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
// that coupon inside the box — the support's box IN THE CLUSTER'S PATCH FRAME (`axes`,
// `frame_min` / `frame_max`; decision 715: the mesh AABB of a rotated frame over-owned the
// cells between the coupon's face and the AABB's) — the coupon's twins carry the straight
// continuation of its claims to the box face, so the cell's own patch there is the double
// count. The cell is
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
// lengths show. FOREIGN-CONTEXT OWNERSHIP (decisions 553 / 559, the same pass): the cells
// of a Foreign stretch wholly inside a contract-3 support (the record above) are classified
// ONCE on the intact patches, before the continuation pass (the same classification the
// stretch records carry), and treated PER CELL by piece class after it — a cell on a
// Chain: true piece is owned by the coupon (weight 0, quadrature weight 0, cell {0, 0},
// listed under `cells` like a wholly owned continuation cell so the portion audit
// reconciles); a cell on a Chain: false piece (an excluded ForeignEdge) keeps its weight
// and sets `domain_weight` 0 (listed under `foreign_excluded_cells`); a cell on no piece is
// untouched. A cell inside two supports takes the stronger class (owned over
// domain-dropped) and lists every support of that class as an owner (ascending spatial
// patch; the owned length split equally): deterministic in the support order. A cell
// already at weight 0 (owned by the continuation pass) is skipped. Lengths in patch units.
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
    // The cell is owned as a ForeignContextEdge cell (decision 553 sub-class (ii)), not
    // as the continuation of a claim.
    bool foreign_context = false;
  };
  std::vector<Cell> cells;  // in patch order
  // Cells on an excluded ForeignEdge (decision 553 sub-class (i)): weight kept,
  // domain_weight 0.
  struct ForeignExcludedCell
  {
    std::size_t patch = 0;
    int feature = -1;
    int stretch = -1;
    // The spatial patches whose ForeignEdge the cell lies on (ascending; every support that
    // classes the cell so, decision 559 MINOR-6).
    std::vector<std::size_t> owners;
    double cell_length = 0.0;
  };
  std::vector<ForeignExcludedCell> foreign_excluded_cells;  // in patch order
  double foreign_excluded_length = 0.0;
  double foreign_context_length = 0.0;  // the ForeignContextEdge part of owned_length
  int foreign_context_cells = 0;
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
// The Diagnostics entry of the corner-arm extension (decision 511 O2 (i); fillet-basis
// design section 7.1): the matched corners whose claim ends beyond the square exit on an
// arm (rounded corners), the cells extended back to the exit and any unhosted stretch;
// empty-but-complete when nothing was extended. The one-line summary and the warning on an
// unhosted stretch read that entry.
nlohmann::json DescribeCornerArmExtensions(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::CornerArmExtensionData>
        &extensions,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);
std::string DescribeCornerArmExtensionSummary(const nlohmann::json &diagnostics);
std::string DescribeCornerArmExtensionUnhostedWarning(const nlohmann::json &diagnostics);
// The Diagnostics entry of the mirror arm trim (boundary-cut DESIGN 2.2.3): the virtual
// corners placed at weight 1 / 2 and the real-arm cells beginning at s_half.
nlohmann::json DescribeMirrorArmTrims(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorArmTrimData> &trims,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale);
// The per-Edge weights of a model placed on an unmerged mirror-formed configuration by its
// REQUIREMENT CONTRACT (decision 557 (4); impl-B5 CONTRACT.md sections 2-4): 1.0 for an
// Edge on one of the contract's RealPortions (a real edge), 0.0 for every other Edge (an
// image edge, never applied to a cell). `edge_portions` are the model Edges' Signature
// portion indices (VerifySpatialEdgesInSignatureFrame), `edge_weights` the ENTRY's per-Edge
// Weight when stamped (nullopt otherwise), `entry` the model's MirrorFormed ENTRY record
// (null for an ordinary coupon of the same key, e.g. a real coupon of a symmetric key).
// Fails closed by name (MFEM_ABORT) on an inconsistent record: an Edge on no portion, a
// stamped Weight disagreeing with the contract, ENTRY RealPortions differing from the
// contract's, an ENTRY RealLengthFraction off the contract's by more than the signature
// tolerance. Returns the weights (parallel to the Edges) and the real length fraction.
struct MirrorFormedEdgeWeights
{
  std::vector<double> weights;
  double real_length_fraction = 1.0;
};
MirrorFormedEdgeWeights MirrorFormedEdgeWeightsOf(
    const nlohmann::json &contract, const std::vector<int> &edge_portions,
    const std::vector<std::optional<double>> &edge_weights, const nlohmann::json &entry,
    const std::string &model_name, int feature_id);

// The Diagnostics entry of the mirror-formed contract placement (decision 557 (4); impl-B5
// CONTRACT.md section 4): every spatial patch placed on an unmerged mirror-formed
// configuration by its contract - the model, the weight (the real length fraction), the
// per-Edge weights (1 real / 0 image), the real claims, the even-extension reach.
nlohmann::json DescribeMirrorFormedPlacements(
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
  // Why the patch is DomainBoundary rather than Mirrored (boundary-cut DESIGN 2.2.5;
  // Diagnostics.DomainBoundary.Reasons): MirrorOff (no mirror: the band is 0 - Mirror
  // "Off" or a Maxwell build), UnmergedTopology (its own footprint overlaps a real portion
  // of a mirror-formed configuration without a mirror placement, read as a Missing feature:
  // decision 481), NonMirroringPlane (an outside point lies beyond an Unsupported /
  // NonPlanar plane), BeyondBand (beyond a Natural plane by more
  // than the band), ReflectionNotLocated (reflected into the domain but not located: a
  // non-convex domain, a hole), NoPlane (an outside point beyond no truncation plane).
  std::string reason;
  // Reason UnmergedTopology (decision 481): the unmerged configuration whose REAL portion
  // the patch's own footprint overlaps (index into MirrorBand.UnmergedFeatures), that
  // portion (mesh units) and its distance from the patch origin (0 for an own cell).
  int unmerged_topology = -1;
  std::array<double, 3> unmerged_p0{}, unmerged_p1{};
  bool unmerged_image = false;
  double unmerged_distance = 0.0;
  bool unmerged_own = false;  // the patch's own footprint overlaps a real portion
};
// Information only (decision 481): a patch whose coupon support reaches
// a portion, real or image, of an unmerged mirror-formed configuration, with the nearest
// such portion's distance from the patch origin (mesh units). The patch keeps its
// classification.
struct DomainBoundaryUnmergedReach
{
  std::size_t patch = 0;
  int topology = -1;
  double distance = 0.0;
  bool image = false;
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
  std::vector<DomainBoundaryExclusion> patches;             // ascending patch index
  std::vector<DomainBoundaryMirrored> mirrored;             // ascending patch index
  std::vector<DomainBoundaryUnmergedReach> unmerged_reach;  // ascending patch index
  long long int tested_patches = 0;
  long long int tested_points = 0;
  long long int reflected_points = 0;  // outside points reflected into the mesh
  double wall_time = 0.0;              // of the containment test, seconds
  // Cumulative phase stamps of the test (seconds since its start: points, locator, find,
  // reduce, reflect, distances), printed with PALACE_RESPONSE_SETUP_TIMING (decision 538).
  std::map<std::string, double> phase_times;
};
// `mirror_planes` (empty: no mirror, every cut-crossing patch is DomainBoundary) and the
// band (mesh units) classify the outside points: a point beyond a Natural plane within the
// band is reflected and located again (the found flags of the reflections OR-reduced like
// the originals); a patch whose every outside point is located after reflection is
// Mirrored, any other patch with an outside point is DomainBoundary (weight 0). A
// metal-edge reference outside the mesh whose reflection is located is a cut-adjacent
// coupon (decision 314's follow-up), never the misplaced-coupon abort. A patch whose own
// footprint overlaps a REAL portion of an unmerged mirror-formed configuration
// (`unmerged_portions`) is DomainBoundary whatever its points' containment (decision 481:
// the configuration is read as a Missing feature); the coupon supports' reach into the
// configurations (real or image portions) is recorded as information.
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
    double mirror_band = 0.0,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::UnmergedPortionData>
        &unmerged_portions = {});

// The per-process cache of the DERIVED model response matrices (decision 734 F1; the
// response-setup PROFILE 7.1): BuildDomainResponseMatrices (the unit-scaled fabricated /
// thin matrices, the defect, the fixed-flux transform and defect: an Eigen eigensolve and
// five n x n products per model) and BuildSurfaceResponseMatrices (the per-interface sums)
// are recomputed at every AMR cycle on inputs that do not change across cycles - the
// model's files, its ZeroTraceIndices and the units. Their outputs are cached per process,
// keyed on the source files' ResponseMatrixFileKey (path / size / time / expected size /
// within-R), the blend weights, the zero-trace indices, the interface mappings and the two
// unit scale factors; a hit returns the first computation's bits, so every downstream
// quantity is byte-identical by construction. PALACE_RESPONSE_DERIVED_CACHE=0 recomputes
// on every hit and fails closed unless the recomputation equals the cached matrices bit
// for bit (the validation switch). The counters are process-wide (every construction,
// electrostatic and Maxwell).
struct DerivedResponseMatrixCacheCounters
{
  long long int hits = 0;    // cached derived matrices returned (recomputed-and-verified
                             // ones included under PALACE_RESPONSE_DERIVED_CACHE=0)
  long long int misses = 0;  // derived matrices computed and stored
  long long int entries = 0;
};
DerivedResponseMatrixCacheCounters GetDerivedResponseMatrixCacheCounters();

// The mirror planes of a response configuration as the mirror module's planes.
std::vector<MirrorPlane> MirrorPlanesOf(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData> &planes);

// Reflect every point outside the mesh's Natural mirror planes into the domain (the even
// extension of the trace; boundary-cut DESIGN 2.2.3): returns the number of reflected
// points; a point beyond a plane that does not mirror or beyond the band is left in place
// (the locator then fails closed naming the patch). `point_band` (optional, one entry per
// point, mesh units) widens the band for the points of a mirror-formed cluster coupon
// placed by its contract (decision 557 (4)): the coupon straddles its planes and its image
// half, the exact reflection of its real half, reaches the coupon's own extent (the
// patch's Provenance::mirror_reach); 0 keeps the global band.
int ReflectPointsIntoDomain(
    mfem::Vector &xyz, int dimension,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::MirrorPlaneData>
        &mirror_planes,
    double mirror_band, double tolerance, const std::vector<double> *point_band = nullptr);

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

// The real arm of a virtual (mirror-formed) corner for the mirror arm trim: the direction
// from the vertex along the feature's first real straight portion (its arc chords, when the
// corner is rounded, are skipped; the first real portion without a straight one; empty
// without a real portion). Segments at or beyond `real_segments` are images.
std::optional<std::array<double, 3>>
MirrorArmDirection(const IdentificationResult &identification,
                   const IdentifiedFeature &feature, std::size_t real_segments);

// Clip the longitudinal cell of a translational patch to the kept offsets [kept_lo,
// kept_hi] of its cell (the corner / mirror arm trims and the continuation ownership): the
// origin moves to the kept midpoint, the cell becomes symmetric about it, the weight and
// the provenance quadrature weight scale by kept / cell, the Maxwell anchors move with the
// origin, the recorded own-segment pre-image (Provenance::own_cell, decision 537) maps
// linearly onto the kept part; an empty kept interval leaves weight 0 and cell {0, 0}.
void ClipLongitudinalCell(
    config::ElectrostaticSolverData::ResponseCorrectionPatchData &patch, double kept_lo,
    double kept_hi);

// The own-edge point of a translational patch at the cell offset c: the recorded pre-image
// on the own segment (Provenance::own_cell) mapped linearly from the cell, the frame
// reconstruction origin + EdgeOffset AxisU + c AxisW without the record (decisions 537 /
// 545: the raw claims of a DomainBoundary cell and the parts of a cell attributed to a
// DomainBoundary spatial coupon are built from it).
std::array<double, 3>
OwnEdgePointAt(const config::ElectrostaticSolverData::ResponseCorrectionPatchData &patch,
               double c);

// The SpatialClusterPriority pass of ApplyPlacementOwnership on a set of spatial supports
// (the preflight's reading of the operator's rule): every pair of supports not already
// owned and not of one interpolation group whose boxes interior-overlap with exactly one
// of them a spatial edge cluster yields (subordinate patch, owner patch); a placeholder
// support (from_signature_box: a Missing key, no coupon) and cluster / cluster pairs (the
// margin-overlap record) are skipped. Deterministic in the supports' order.
std::vector<std::pair<std::size_t, std::size_t>> FindSpatialClusterPriorityExclusions(
    const std::vector<SpatialSupportBounds> &supports,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const std::set<std::size_t> &owned_patches, int dimension);

// The dropped coupons' raw claims (provenance raw_claims: the arms' R windows) clipped by
// ClipUncoveredPortionsBySpatialSupport against every matched cluster's box; the kept
// (outside) pieces are appended by the DomainBoundary rule (AppendExcludedPatchPortion:
// the tolerance floor, the deduplication).
SpatialClusterPriorityExclusions CollectSpatialClusterPriorityPortions(
    const std::vector<std::pair<std::size_t, std::size_t>> &exclusions,
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const std::vector<SpatialSupportBounds> &supports, int dimension,
    double matching_radius);

// The Diagnostics entry "SpatialClusterPriority" (decision 915 MINOR-7 / 922): Count, the
// dropped patches (patch, feature, model, topology, origin, weight, box, owner, window /
// outside lengths), the raw portions (DescribeExcludedPatchPortions) and the rule.
nlohmann::json DescribeSpatialClusterPriorityExclusions(
    const SpatialClusterPriorityExclusions &exclusions,
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const std::vector<SpatialSupportBounds> &supports, double coordinate_scale);

ExcludedPatchPortions CollectDomainBoundaryPortions(
    const DomainBoundaryExclusions &exclusions,
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const ContinuationOwnership &ownership,
    const UncoveredSpatialSupportClipping &uncovered_clipping, double matching_radius);

// The raw claim of ONE excluded patch appended to `result` by the rule above (the per-patch
// step of CollectDomainBoundaryPortions, shared with the conductor-consistency consumer of
// S7diag R1, decision 758): `translational` marks the own-edge interval of a cell; a
// portion shorter than `tolerance` (kSignatureParameterToleranceOverRadius x R, mesh units)
// is dropped; a portion of the same feature already in `result` within `tolerance` at both
// ends folds the patch into that entry (DuplicatePatches).
void AppendExcludedPatchClaim(
    ExcludedPatchPortions &result, std::size_t patch_idx,
    const std::vector<config::ElectrostaticSolverData::ResponseCorrectionPatchData>
        &patches,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    const ContinuationOwnership &ownership,
    const UncoveredSpatialSupportClipping &uncovered_clipping, double tolerance);
void AppendExcludedPatchPortion(
    ExcludedPatchPortions &result, std::size_t patch_idx,
    const config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData
        &portion,
    bool translational, double tolerance);

// The Diagnostics entry of the raw portions (count, geometric cells, duplicates, length,
// per type and per feature, the portion list; mesh units scaled by coordinate_scale) with
// the consumer's rule text; DescribeDomainBoundaryPortions is the F-DB-a entry.
nlohmann::json DescribeExcludedPatchPortions(const ExcludedPatchPortions &portions,
                                             double coordinate_scale,
                                             const std::string &rule);
nlohmann::json DescribeDomainBoundaryPortions(const ExcludedPatchPortions &portions,
                                              double coordinate_scale);

// The integration intervals of one Maxwell response contour line (ConfigureMaxwellLines)
// from the exact element intersection intervals gathered over the ranks (line parameter
// in [0, 1]): the ends are snapped to 0 / 1 and the breakpoints uniqued at the parameter
// tolerance 1e-10, every piece between breakpoints shorter than that is skipped, and a
// remaining piece covered by no element interval is either a FLOATING-POINT CRACK or a
// true gap (decision 745 (1), lane contourcover). Two tets sharing the face the line
// crosses clip it independently from their own vertex arithmetic; the barycentric
// coordinates carry a rounding error of O(eps x coordinate magnitude / altitude), so the
// two clips of one face differ by a few tens of ULP of the coordinate magnitude in PHYSICAL
// length (measured 20 / 68 x eps x coordinate_scale on the two cracks of record), a gap
// or an overlap at random; in parameter space that is ~2e-10 on a 158-nm line, above the
// 1e-10 constant, and the Maxwell driven runs aborted at AMR depth. A crack is closed
// (dropped from the integration like the < 1e-10 pieces, counted and its length kept) iff
// its physical length (piece x line_length) <= kMaxwellLineCrackUlps x eps x
// max(coordinate_scale, line_length), i.e. 2^12 ULP of the coordinate magnitude: >= 60x
// the measured cracks and a near-tangential face crossing's 1 / sin(theta) down to ~1 deg,
// while 3 orders below a mesher's geometric tolerance (1e-8 x Lc), so a line that leaves
// the mesh over any geometric length still fails closed (MFEM_VERIFY naming the line, the
// parameter range, the gap and the tolerance). Where no uncovered piece occurs the
// intervals are the pre-existing ones to the bit.
//
// A piece that is a true gap for the exact intervals is still integrated when the line
// lies IN a mesh face no element's exact intersection covers (decisions 829 / 831, lane
// inplane): Palace cracks the internal PEC boundaries and moves the two sides' duplicated
// face vertices apart by CrackDisplacementFactor x h_min (geodata.cpp), so a line lying in
// such a sheet plane (a cluster's anchor / hat line on the metal) sits inside the open
// crack, outside both adjacent tets by delta / h ~ 1e-14 .. 1e-12 in barycentric terms,
// beyond the 64-eps clip, while a line in an uncracked interior face is kept by both tets.
// `in_face_intervals` are the same candidates' intersections with the face the line is
// coplanar with (within the point locator's 1e-9 barycentric acceptance over the line's
// portion inside the element's bounding box: a property of the line near the element, not
// of its far ends, so that the two grazing faces through a mesh edge the line follows are
// skipped as well and the interval is the clip by the transversal faces only, decision 862
// / OI-LOC-HANG: the C4 iteration-11 abort of the Maxwell-at-scale pilot, where a 2.4-um
// hat line following the trace edge over 20-nm elements of a nonconforming fan was cut by
// those faces at ill-conditioned crossings and a 0.09-nm stretch stayed uncovered) taken as
// not clipping: such a piece is split at their breakpoints (element-aligned quadrature)
// and integrated once, each quadrature point then located by the rank / element rule of
// DistributedPointLocator (the lowest candidate rank, its lowest element). The one-sided
// assignment is exact where the two sides agree: across an uncracked (conforming) face
// the tangential field is continuous in H(curl), and on a cracked PEC sheet both sides'
// tangential traces are essential (zero), so the integral is the same from either side
// (zero on the sheet). On a cracked face that is not PEC the two sides carry independent
// tangential dofs and this rule integrates the owner's side only; within the scope of
// record (clusters on PEC metal; finite-impedance Maxwell patches with ZeroTraceIndices
// refused) no such line occurs. Pieces covered by exact intervals and closed cracks are
// untouched by this pass: a line without a true gap is bitwise the pre-existing
// resolution. A piece neither pass covers still fails closed.
constexpr double kMaxwellLineCrackUlps = 4096.0;

// The bounding-box tolerance of the Maxwell line locator (ConfigureMaxwellLines): the
// candidate search and the box portion of the coplanarity test above.
double MaxwellLineBoxTolerance(double coordinate_scale);

// The exact and in-face intersection intervals (line parameter in [0, 1]) of one segment
// with the elements of the LOCAL mesh (ElementPointLocator::FindSegmentIntersections with
// the box tolerance of `coordinate_scale`): the per-rank input of the resolution below,
// exposed for the unit tests of the locator (decision 862).
struct MaxwellSegmentIntersections
{
  std::vector<std::pair<double, double>> exact;
  std::vector<std::pair<double, double>> in_face;
};
MaxwellSegmentIntersections
FindMaxwellSegmentIntersections(mfem::ParMesh &mesh, int dimension,
                                const std::array<double, 3> &begin,
                                const std::array<double, 3> &end, double coordinate_scale);

struct MaxwellLineCoverage
{
  std::vector<std::pair<double, double>> intervals;
  int closed_crack_count = 0;
  double max_closed_crack = 0.0;  // Mesh units.
  // Pieces integrated from the in-face intervals (no exact interval covers them).
  int in_face_count = 0;
  double in_face_length = 0.0;  // Mesh units.
};
MaxwellLineCoverage
ResolveMaxwellLineCoverage(std::vector<std::pair<double, double>> element_intervals,
                           std::vector<std::pair<double, double>> in_face_intervals,
                           double line_length, double coordinate_scale, int line_index,
                           const std::array<double, 3> &begin,
                           const std::array<double, 3> &end);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
