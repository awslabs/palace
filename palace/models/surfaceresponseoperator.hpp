// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP

#include <array>
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class FiniteElementSpace;
class GridFunction;
class IoData;
class Mesh;
class BoundaryModeOperator;
class LaplaceOperator;
class MaterialOperator;
class SpaceOperator;

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
  };

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
    std::vector<MortarSegment> mortar_segments;
    std::vector<MortarVertex> mortar_vertices;
    std::vector<MortarTriangle> mortar_triangles;
    mfem::DenseMatrix mortar_mass_inverse;
    mfem::Vector mortar_constant_load;
    std::vector<mfem::Vector> mortar_conductor_loads;
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
      const mfem::Vector &xyz, int dimension,
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
  void ApplyUneliminated(const Vector &x, Vector &y) const;
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

  struct ElectrostaticResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;
    std::map<int, double> trace_closure_spread;
    std::vector<ModelContribution> model_contributions;
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

  // Evaluate the complete fabricated-coupon surface energy for every mapped target
  // interface. Corrected participation replaces the measured global core with this data.
  std::map<int, double> GetFabricatedSurfaceEnergy(const Vector &x) const;

  // Evaluate coupon responses on an electrostatic potential. With fixed flux enabled,
  // return both complete postprocessing-only closure ensembles. Otherwise use the active
  // per-model domain-coupling policy for corrected-domain accounting while retaining the
  // fabricated fixed-trace surface evaluation.
  ElectrostaticResponse GetElectrostaticResponse(const Vector &x,
                                                 bool include_fixed_flux = true) const;

  // Collect the actual local contour and conductor-state coefficients for every
  // three-dimensional spatial response patch. The result is ordered by global patch
  // index and replicated on all ranks.
  std::vector<PatchTrace> GetSpatialPatchTraces(const Vector &x) const;

  // Evaluate a postprocessing-only response for a complex Nedelec Maxwell field. Coupon
  // voltages are reconstructed from transverse contour integrals and applied through
  // Hermitian response-matrix quadratic forms.
  MaxwellResponse GetMaxwellResponse(const GridFunction &E,
                                     std::complex<double> omega) const;

  bool HasSurfaceResponse() const;
  std::set<int> GetTargetInterfaces() const;

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
};

// Classify the configured automatic surface-response neighborhoods without assembling a
// finite-element operator or solving a field problem. The deterministic JSON manifest
// reports exact, interpolated, and missing process-library requirements.
void WriteSurfaceResponseRequirements(const IoData &iodata, const Mesh &mesh,
                                      const std::string &path);

// Ownership record of the translational patches against the spatial supports (decisions
// 224 / 236): a translational STRETCH — the maximal contiguous stretch of one feature side
// along its chain (patch provenance feature / stretch), the identification's portion unit —
// whose every longitudinal cell lies strictly inside one spatial support's box is RECORDED
// (manifest / metadata Diagnostics.TranslationalStretchesInsideSpatialSupport and a
// warning), never an abort: the coupon is calibrated on the cluster's claims and their
// straight continuations to the box face, so a stretch that continues one of the cluster's
// claims through its claim cut (class Continuation) is corrected by both the coupon and its
// own patches, while any other stretch (class Foreign) is absent from the coupon's twins —
// a model mismatch, not a double count. A stretch continues a claim when one of its cells
// lies on the claim's mesh segment (the claim boundary cut that segment) or when it runs
// parallel to the claim (within the signature angle tolerance) and one of its ends abuts a
// claim end along the chain within continuation_tolerance, within R of it transversely (a
// claim boundary snapped onto a mesh vertex); the cells of a pair sit on the pair's
// midline, so collinearity with the edge is not tested. Judged per stretch, never per cell:
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
nlohmann::json DescribeTranslationalOwnershipRecords(
    const std::vector<TranslationalOwnershipRecord> &records,
    const std::vector<SpatialSupportBounds> &supports,
    const config::ElectrostaticSolverData::ResponseCorrectionData &config,
    double coordinate_scale, const std::vector<std::string> &skipped);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
