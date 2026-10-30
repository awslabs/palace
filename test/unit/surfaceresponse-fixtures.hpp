// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_SURFACERESPONSE_FIXTURES_HPP
#define PALACE_TEST_SURFACERESPONSE_FIXTURES_HPP

#include "fixtures.hpp"
#include "utils/tablecsv.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

namespace palace
{

class IoData;

namespace config
{

struct BoundaryData;

}  // namespace config

}  // namespace palace

namespace palace::test
{

using json = nlohmann::json;

// Metal edge segments of one interface of a mesh (ExtractMetalEdgeGeometry): the segment
// indices, their total length and the counts of the physical vertex types they touch.
struct InterfaceEdgeSummary
{
  std::vector<std::size_t> segments;
  double length = 0.0;
  int corners = 0;
  int junctions = 0;
};

// Scope of the response-geometry cache environment of a test: the constructor sets
// PALACE_RESPONSE_GEOMETRY_CACHE to the cache path (and
// PALACE_RESPONSE_GEOMETRY_CACHE_WRITE when write is true), DisableWrite unsets the WRITE
// variable before a reload, and the destructor unsets both, so a failing assertion inside
// the scope cannot leak the cache variables into the following test cases of the same
// process.
class GeometryCacheEnvGuard
{
public:
  GeometryCacheEnvGuard(const std::string &cache_path, bool write);
  GeometryCacheEnvGuard(const GeometryCacheEnvGuard &) = delete;
  GeometryCacheEnvGuard &operator=(const GeometryCacheEnvGuard &) = delete;
  ~GeometryCacheEnvGuard();
  void DisableWrite();
};

// Shared fixture of the SurfaceResponseOperator unit tests (TEST_CASE_METHOD): the
// temporary directory with every synthetic response matrix, basis-point file and
// fabrication-process library the cases read (written by the root rank in the constructor;
// the writes cost nothing measurable), the configuration builders the cases derive their
// IoData from, the synthetic meshes and the check helpers shared by several cases.
struct SurfaceResponseFiles
{
  SharedTempDir temp;
  const fs::path points_path = temp.temp_dir / "basis-points.csv";
  const fs::path fabricated_path = temp.temp_dir / "fabricated.csv";
  const fs::path thin_path = temp.temp_dir / "thin.csv";
  const fs::path fabricated_surface_path = temp.temp_dir / "fabricated-surface.csv";
  const fs::path thin_surface_path = temp.temp_dir / "thin-surface.csv";
  const fs::path compact_fabricated_surface_path =
      temp.temp_dir / "fabricated-surface-compact.csv";
  const fs::path compact_thin_surface_path = temp.temp_dir / "thin-surface-compact.csv";
  const fs::path library_path = temp.temp_dir / "fabrication-process.json";
  const fs::path legacy_library_path = temp.temp_dir / "fabrication-process-legacy.json";
  const fs::path missing_layer_library_path =
      temp.temp_dir / "fabrication-process-missing-layer.json";
  const fs::path invalid_library_path =
      temp.temp_dir / "fabrication-process-invalid-depth.json";
  const fs::path impedance_library_path =
      temp.temp_dir / "fabrication-process-impedance.json";
  const fs::path legacy_impedance_library_path =
      temp.temp_dir / "fabrication-process-impedance-legacy.json";
  const fs::path conductivity_library_path =
      temp.temp_dir / "fabrication-process-conductivity.json";
  const fs::path rational_impedance_library_path =
      temp.temp_dir / "fabrication-process-rational-impedance.json";
  const fs::path invalid_boundary_law_library_path =
      temp.temp_dir / "fabrication-process-invalid-boundary-law.json";
  const fs::path library_3d_path = temp.temp_dir / "fabrication-process-3d.json";
  const fs::path impedance_library_3d_path =
      temp.temp_dir / "fabrication-process-impedance-3d.json";
  const fs::path exact_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-exact-pair-2d.json";
  const fs::path different_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-different-pair-2d.json";
  const fs::path reference_strip_library_2d_path =
      temp.temp_dir / "fabrication-process-reference-strip-2d.json";
  const fs::path reference_gap_library_2d_path =
      temp.temp_dir / "fabrication-process-reference-gap-2d.json";
  const fs::path interpolated_pair_library_2d_path =
      temp.temp_dir / "fabrication-process-interpolated-pair-2d.json";
  const fs::path parallel_cluster_library_2d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-2d.json";
  const fs::path impedance_parallel_cluster_library_2d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-impedance-2d.json";
  const fs::path parallel_cluster_points_2d_path =
      temp.temp_dir / "parallel-cluster-basis-points-2d.csv";
  const fs::path coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-3d.json";
  const fs::path missing_pair_library_3d_path =
      temp.temp_dir / "fabrication-process-missing-pair-3d.json";
  const fs::path interpolated_coupled_library_3d_path =
      temp.temp_dir / "fabrication-process-coupled-interpolated-3d.json";
  const fs::path parallel_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-3d.json";
  const fs::path parallel_cluster_only_library_3d_path =
      temp.temp_dir / "fabrication-process-parallel-cluster-only-3d.json";
  const fs::path disconnected_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-disconnected-cluster-3d.json";
  const fs::path coupled_fabricated_path = temp.temp_dir / "coupled-fabricated.csv";
  const fs::path coupled_thin_path = temp.temp_dir / "coupled-thin.csv";
  const fs::path coupled_fabricated_surface_path =
      temp.temp_dir / "coupled-fabricated-surface.csv";
  const fs::path coupled_thin_surface_path = temp.temp_dir / "coupled-thin-surface.csv";
  const fs::path cluster_fabricated_path = temp.temp_dir / "cluster-fabricated.csv";
  const fs::path cluster_thin_path = temp.temp_dir / "cluster-thin.csv";
  const fs::path cluster_fabricated_surface_path =
      temp.temp_dir / "cluster-fabricated-surface.csv";
  const fs::path cluster_thin_surface_path = temp.temp_dir / "cluster-thin-surface.csv";
  const fs::path zero_trace_path = temp.temp_dir / "zero-trace.csv";
  const fs::path shared_boundary_trace_path = temp.temp_dir / "shared-boundary-trace.csv";
  const fs::path corner_points_path = temp.temp_dir / "corner-basis-points.csv";
  const fs::path corner_fabricated_path = temp.temp_dir / "corner-fabricated.csv";
  const fs::path corner_thin_path = temp.temp_dir / "corner-thin.csv";
  const fs::path corner_constrained_perturbed_fabricated_path =
      temp.temp_dir / "corner-constrained-perturbed-fabricated.csv";
  const fs::path corner_constrained_perturbed_thin_path =
      temp.temp_dir / "corner-constrained-perturbed-thin.csv";
  const fs::path corner_fabricated_surface_path =
      temp.temp_dir / "corner-fabricated-surface.csv";
  const fs::path corner_thin_surface_path = temp.temp_dir / "corner-thin-surface.csv";
  const fs::path convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-3d.json";
  // Corner (3D box) surface files of the within-R regression: whole-box Q_total inflated
  // 3x with the within-R Q_ij unchanged; both scaled 3x; the legacy compact format without
  // the within-R column; the within-R rows at another radius.
  const fs::path inflated_box_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-3d-inflated-box.json";
  const fs::path scaled_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-3d-scaled.json";
  const fs::path legacy_compact_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-3d-legacy-compact.json";
  const fs::path other_radius_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-3d-other-radius.json";
  const fs::path finite_impedance_convex_library_3d_path =
      temp.temp_dir / "fabrication-process-convex-finite-impedance-3d.json";
  const fs::path concave_library_3d_path =
      temp.temp_dir / "fabrication-process-concave-3d.json";
  const fs::path strip_library_3d_path =
      temp.temp_dir / "fabrication-process-strip-3d.json";
  const fs::path rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-3d.json";
  // The rounded corner with PEC knots on BOTH metal rings (the sheet, z = 0, and the slab
  // top, z = +0.02): the top ring's one free knot is reached by a free arc that starts at a
  // slab-top PEC knot on main (the E1 defect) and at the conductor anchor from decision
  // 818 R1 on.
  const fs::path slab_top_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-slab-top-3d.json";
  // The slab-top rounded corner with an explicit TraceMesh whose ring edges carry slave
  // vertices (the corner family's trace basis form: a vertex between two consecutive ring
  // knots, no coefficient of its own): refused on the Maxwell path before decision 818 R3,
  // routed knot -> slave -> knot from then on (two slaves on the free bottom ring, one on
  // the top ring's free-to-PEC step, one on its PEC-to-free step that the anchor line
  // replaces).
  const fs::path slab_top_rounded_slave_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-slab-top-slaves-3d.json";
  const fs::path slab_top_rounded_slave_trace_vertices_path =
      temp.temp_dir / "rounded-slab-top-slaves-trace-vertices.csv";
  const fs::path slab_top_rounded_slave_trace_triangles_path =
      temp.temp_dir / "rounded-slab-top-slaves-trace-triangles.csv";
  const fs::path constrained_perturbed_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-constrained-perturbed-3d.json";
  const fs::path finite_impedance_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-finite-impedance-3d.json";
  const fs::path rounded_concave_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-concave-3d.json";
  const fs::path interpolated_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-interpolated-3d.json";
  const fs::path unqualified_interpolated_rounded_library_3d_path =
      temp.temp_dir / "fabrication-process-rounded-interpolated-unqualified-3d.json";
  const fs::path endpoint_library_3d_path =
      temp.temp_dir / "fabrication-process-endpoint-3d.json";
  const fs::path junction_library_3d_path =
      temp.temp_dir / "fabrication-process-junction-3d.json";
  const fs::path spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-3d.json";
  const fs::path cross_layer_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cross-layer-3d.json";
  const fs::path incomplete_cross_layer_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cross-layer-incomplete-3d.json";
  const fs::path spatial_cluster_position_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-position-mismatch-3d.json";
  const fs::path spatial_cluster_orientation_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-orientation-mismatch-3d.json";
  const fs::path spatial_cluster_interval_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-interval-mismatch-3d.json";
  const fs::path spatial_cluster_extra_edge_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-extra-edge-3d.json";
  const fs::path spatial_cluster_impedance_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-impedance-mismatch-3d.json";
  const fs::path spatial_cluster_mixed_impedance_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-mixed-impedance-3d.json";
  const fs::path spatial_cluster_parameter_mismatch_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-parameter-mismatch-3d.json";
  const fs::path spatial_cluster_points_path =
      temp.temp_dir / "spatial-cluster-basis-points.csv";
  // Cap-interior hats (InteriorTraceCount, decision 112(b)): the spatial cluster with two
  // trailing basis points on no contour.
  const fs::path cap_hat_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cap-hats-3d.json";
  const fs::path cap_hat_loaded_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cap-hats-loaded-3d.json";
  const fs::path cap_hat_full_partition_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cap-hats-full-partition-3d.json";
  const fs::path cap_hat_without_trace_mesh_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cap-hats-no-trace-mesh-3d.json";
  const fs::path cap_hat_translational_library_3d_path =
      temp.temp_dir / "fabrication-process-cap-hats-translational-3d.json";
  // The cap-hat cluster whose first hat sits at a plane no ring knot shares (local w = 0):
  // the Maxwell path refuses it (decision 821 (2)); the electrostatic path loads it.
  const fs::path cap_hat_off_plane_spatial_cluster_library_3d_path =
      temp.temp_dir / "fabrication-process-spatial-cluster-cap-hats-off-plane-3d.json";
  // The cap-hat cluster (open paths between its two conductors) declaring ZeroTraceIndices,
  // resp. ContourGroups, as well: a closed-contour key on an open-path model, refused at
  // the library load on both paths.
  const fs::path cap_hat_zero_trace_open_paths_library_3d_path =
      temp.temp_dir /
      "fabrication-process-spatial-cluster-cap-hats-zero-trace-open-paths-3d.json";
  const fs::path cap_hat_contour_groups_open_paths_library_3d_path =
      temp.temp_dir /
      "fabrication-process-spatial-cluster-cap-hats-contour-groups-open-paths-3d.json";
  const fs::path cap_hat_off_plane_points_path =
      temp.temp_dir / "cap-hat-off-plane-basis-points.csv";
  const fs::path cap_hat_off_plane_trace_vertices_path =
      temp.temp_dir / "cap-hat-off-plane-trace-vertices.csv";
  const fs::path cap_hat_points_path = temp.temp_dir / "cap-hat-basis-points.csv";
  const fs::path cap_hat_trace_vertices_path = temp.temp_dir / "cap-hat-trace-vertices.csv";
  const fs::path cap_hat_trace_triangles_path =
      temp.temp_dir / "cap-hat-trace-triangles.csv";
  const fs::path cap_hat_fabricated_path = temp.temp_dir / "cap-hat-fabricated.csv";
  const fs::path cap_hat_thin_path = temp.temp_dir / "cap-hat-thin.csv";
  const fs::path cap_hat_fabricated_surface_path =
      temp.temp_dir / "cap-hat-fabricated-surface.csv";
  const fs::path cap_hat_thin_surface_path = temp.temp_dir / "cap-hat-thin-surface.csv";
  const fs::path cap_hat_loaded_fabricated_path =
      temp.temp_dir / "cap-hat-loaded-fabricated.csv";
  const fs::path cap_hat_loaded_thin_path = temp.temp_dir / "cap-hat-loaded-thin.csv";
  const fs::path cap_hat_loaded_fabricated_surface_path =
      temp.temp_dir / "cap-hat-loaded-fabricated-surface.csv";
  const fs::path cap_hat_loaded_thin_surface_path =
      temp.temp_dir / "cap-hat-loaded-thin-surface.csv";
  const fs::path cross_layer_fabricated_surface_path =
      temp.temp_dir / "cross-layer-fabricated-surface.csv";
  const fs::path cross_layer_thin_surface_path =
      temp.temp_dir / "cross-layer-thin-surface.csv";
  const json impedance_law = {{"Type", "Impedance"}, {"Ls", 1.0e-13}};
  const json second_impedance_law = {{"Type", "Impedance"}, {"Ls", 2.0e-13}};
  const json conductivity_law = {{"Type", "Conductivity"},
                                 {"Conductivity", 5.8e7},
                                 {"Permeability", 1.2},
                                 {"Thickness", 1.0e-7},
                                 {"External", true}};
  const json rational_impedance_law = {{"Type", "RationalImpedance"},
                                       {"Numerator", {5.0e-8, 0.0}},
                                       {"Denominator", {5.0e-20, 1.0e-9, 50.0}}};

  SurfaceResponseFiles();

  // Configurations (Problem.Output = the temporary directory).
  // The 2D automatic library problem: an 8 x 4 triangle mesh with the cracked metal
  // attributes 9 / 10 on y = 0.5 (MakeAutomatic2DMesh), one SA interface at R = 0.1.
  json AutomaticConfig2D() const;
  // The same problem as a BoundaryMode solve (PEC metal, Freq 5).
  json BoundaryModeConfig2D() const;
  // The automatic problem with the exact same-conductor pair library at R = 0.3.
  json ExactPairConfig2D() const;
  // The cpw3d_surface example (examples/cpw3d_surface/cpw3d_surface_validation_thin.json)
  // on the test mesh cpw3d-surface-nc.msh with the legacy 3D library and PatchConstruction
  // "Legacy".
  json Cpw3dLegacyConfig() const;
  // A closed rectangular PEC island (attribute 9) on the plane y = 0.5 of a box with the
  // concave 3D library and PatchConstruction "Legacy" (MakeIslandMesh).
  json IslandConfig() const;
  json ConvexIslandConfig() const;
  json ConvexMaxwellIslandConfig() const;
  json RoundedIslandConfig() const;
  // Two islands (attributes 9 / 10, two terminals) with the spatial cluster library and
  // UnmatchedPolicy "Warn" (MakeIslandMesh(..., neighboring_island = true)).
  json HighOrderSpatialConfig() const;

  // An IoData whose operators are built on the mesh another IoData read with mesh::ReadMesh
  // (the reading IoData records the cracked boundary attributes; both must describe the
  // same boundaries).
  static void ShareCrackedAttributes(IoData &iodata, const IoData &reader);

  // Meshes.
  static mfem::Mesh MakeAutomatic2DMesh();
  // normal_elements: the element count across the plane normal y (4: planes at 0.25,
  // 0.5, 0.75; 16: a plane every 0.0625, so a trace with kinks at multiples of 0.0625
  // above and below the metal is a Q1 function). MakeSerialIslandMesh is the serial mesh
  // MakeIslandMesh distributes (written to a file by the driver cases).
  static mfem::Mesh MakeSerialIslandMesh(bool rounded = false, bool tetrahedral = false,
                                         bool aperture = false,
                                         bool neighboring_island = false,
                                         bool second_layer = false,
                                         bool high_order_rounded = false,
                                         int normal_elements = 4);
  static std::unique_ptr<mfem::ParMesh>
  MakeIslandMesh(bool rounded = false, bool tetrahedral = false, bool aperture = false,
                 bool neighboring_island = false, bool second_layer = false,
                 bool high_order_rounded = false, int normal_elements = 4);
  static std::unique_ptr<mfem::ParMesh> MakeTouchingIslandMesh();
  // The plain island whose metal boundary elements are split by the island's diagonal
  // x + z = 1 into attributes 9 (below) and 10 (above): one conductor, two attributes, so
  // a dielectric interface per attribute gives the two diagonal corners (0.25, 0.75) and
  // (0.75, 0.25) one arm per target map (the OI-MX-R2A-COMPARATOR-FORM of readiness
  // hardening-2 (iii): a listed gap next to a differently listed one at a one-slot corner
  // coupon).
  static std::unique_ptr<mfem::ParMesh> MakeDiagonalIslandMesh();
  // The neighboring-island mesh (two 0.5 x 0.5 islands, attributes 9 / 10, on the plane
  // y = 0.5 of a 2 x 1 x 2 box; MakeIslandMesh(neighboring_island = true)) rotated by
  // `angle_degrees` about the plane normal y through the box centre (1, y, 1): the islands'
  // corner frames lie off the mesh axes while the geometry (the corners 0.5 apart along
  // each island's edges, the facing pair 0.25 apart) is the unrotated one's. The
  // operator-level B2 fixture of readiness hardening-2 (iv): the MESH-ALIGNED box of a
  // coupon cube widens with the rotation (half-width h sqrt 2 at 45 deg), its ORIENTED box
  // does not.
  static std::unique_ptr<mfem::ParMesh> MakeRotatedIslandMesh(double angle_degrees);
  // Two disconnected rectangular islands (attributes 9 / 10) on the plane y = 0.5 of a
  // 2 x 1 x 2 box whose facing corners (1.0, 0.75) and (1.125, 0.875) in (x, z) form the
  // offset-corner-pair spatial cluster; MakeSerialOffsetCornerPairMesh is the serial mesh
  // MakeOffsetCornerPairMesh distributes (written to a file by the driver cases). With
  // `straddling_island` the box is 2 x 1 x 3 and a third 0.5 x 0.5 island (attribute 11)
  // x in [0.125, 0.625], z in [1.75, 2.25] sits across the cluster's contract-box face
  // z = 1.93 (R = 0.2): its lower corners' +z arm windows exit the box by 0.02 (the T5
  // excluded-corner containment fixture, OI-SCP-BOXGROW); 0.5 from island 2's top-left
  // corner and 0.5 wide (beyond the 2R interaction: a free island, never a claim).
  static mfem::Mesh MakeSerialOffsetCornerPairMesh(bool straddling_island = false);
  static std::unique_ptr<mfem::ParMesh>
  MakeOffsetCornerPairMesh(bool straddling_island = false);
  // The SA edge segments of interface `interface` (with their length and vertex types).
  static InterfaceEdgeSummary
  SummarizeInterfaceEdges(const mfem::ParMesh &mesh, const config::BoundaryData &boundaries,
                          int interface);

  // The constant Maxwell probe field (0.7, -0.4, 0.2) of the 3D cases.
  static mfem::VectorConstantCoefficient ConstantFieldCoefficient();

  // Check helpers of the signature (version-2) library contracts.
  static std::vector<std::string> SignatureInterfaceSets(const json &signature);
  static json SignatureInterfaces(const json &signature);
  static json ChordedSignatureEdges(const json &signature, double radius);
  static std::vector<std::vector<std::string>> ReadPatchRows(const fs::path &path);
  static std::size_t
  CheckPlacedModelEdges(const json &feature, const json &identification,
                        const std::vector<std::vector<std::string>> &patch_rows,
                        const json &model_edges, double tolerance);
  // The geometry-fallback invariant (readiness hardening DESIGN section 4, U-INV; USER
  // decision 1000 (b), supervisor 1001 / 1008): the preflight manifest's
  // Identification.Diagnostics.GeometryFallbacks.Entries (P) and the operator's
  // GetGeometryFallbackRecord() Entries (S; null without a record) agree on (Class,
  // Features, Patches, Models, Length to 1e-12: one unit on both records) for every
  // geometry class; the solve-only NoMatchedFeature (A10.1) has no entry in P and is
  // predicted by Summary.Counts (Exact + Interpolated + LegacyContract + QuantumNearMatched
  // == 0), which is asserted when S carries it. Returns the number of entries compared.
  static std::size_t CheckGeometryFallbackInvariant(const json &manifest,
                                                    const json &solve_record);
};

// Shared helpers of the electrostatic surface-response tests: a whole electrostatic run of
// a config (output folder created), the raw text of a file, a numeric CSV table and one of
// its columns by header, and the plain text table of the per-type energy CSVs
// (surface-response-uncovered-energy.csv, surface-response-domain-boundary-energy.csv).
void RunElectrostatic(json config);
std::string ReadFile(const std::filesystem::path &path);
Table LoadCsv(const std::filesystem::path &path);
const Column &ColumnByHeader(const Table &table, const std::string &header);
struct TextCsv
{
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;
  std::size_t Column(const std::string &name) const;
};
TextCsv ReadTextCsv(const std::filesystem::path &path);

}  // namespace palace::test

#endif  // PALACE_TEST_SURFACERESPONSE_FIXTURES_HPP
