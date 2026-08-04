// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_SINGULAR_SOLVER_HPP
#define PALACE_DRIVERS_SINGULAR_SOLVER_HPP

#include <complex>
#include <memory>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json_fwd.hpp>

#include "fem/singularfeatures.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;

namespace mesh
{

struct PartitionMetadata;

}  // namespace mesh

std::vector<fem::singular::TriangleMaterial>
GetSingularTriangleMaterials(const IoData &iodata);

// A real singular exponent remains exact for complex permittivities only when
// every material sector of a finite transmission wedge has the same loss
// tangent. Thin-sheet endpoints and edges are geometry-driven and are exempt.
void ValidateSingularLossTangents(const IoData &iodata, const mfem::Mesh &mesh,
                                  const fem::singular::FeatureTopology &features);
void ValidateSingularLossTangents(const IoData &iodata,
                                  const fem::singular::TriangleFeatureTopology &features);

// Finite Robin coefficients are lower-order edge terms. Free singular
// tangential traces are retained when square integrable (nu > 1/2). At the
// thin-sheet threshold they are excluded from the impedance energy space and
// constrained like PEC enrichment traces, while standard impedance trace DOFs
// remain free. Mixed PEC/impedance feature junctions still require a distinct
// trace space and are rejected.
void ValidateSingularImpedanceFeatures(const IoData &iodata,
                                       const fem::singular::FeatureTopology &features);
void ValidateSingularImpedanceFeatures(
    const IoData &iodata, const fem::singular::TriangleFeatureTopology &features);
std::set<int>
GetConstrainedSingularImpedanceAttributes(const IoData &iodata,
                                          const fem::singular::FeatureTopology &features);
std::set<int> GetConstrainedSingularImpedanceAttributes(
    const IoData &iodata, const fem::singular::TriangleFeatureTopology &features);

nlohmann::json GetSingularSurfaceParticipationMetadata(const IoData &iodata);
nlohmann::json
GetSingularSurfaceIntegrabilityMetadata(const fem::singular::FeatureTopology &features);
nlohmann::json GetSingularSurfaceIntegrabilityMetadata(
    const fem::singular::TriangleFeatureTopology &features);

// Exact decomposition-independent identities for vertices of a nonconforming
// refinement tree. The root rank retains a descriptor-to-ID dictionary so existing
// vertices keep their identities while newly refined vertices receive fresh IDs.
class NonconformingVertexIdentity
{
private:
  class Impl;
  std::unique_ptr<Impl> impl;

public:
  NonconformingVertexIdentity();
  ~NonconformingVertexIdentity();
  NonconformingVertexIdentity(const NonconformingVertexIdentity &) = delete;
  NonconformingVertexIdentity &operator=(const NonconformingVertexIdentity &) = delete;

  void Clear();
  void Observe(const mfem::ParMesh &mesh,
               const std::vector<fem::singular::GlobalVertexId> &vertex_ids);
  void Update(const mfem::ParMesh &mesh,
              std::vector<fem::singular::GlobalVertexId> &vertex_ids);
};

// Exact topological ancestry for vertices created by conforming refinement. MFEM's
// coarse-to-fine reference embeddings identify the two half-weight barycentric parent
// vertices of every appended edge midpoint. Combining those parent vertices' persistent
// identities yields a descriptor independent of child ordering, MPI decomposition,
// transient global numbering, and physical coordinates.
// Topologically distinct coincident edges stay distinct because their endpoint
// identities differ. Repeated refinement composes: each pass keys on the identities its
// parents already received.
class ConformingVertexAncestry
{
private:
  class Impl;
  std::unique_ptr<Impl> impl;

public:
  ConformingVertexAncestry();
  ~ConformingVertexAncestry();
  ConformingVertexAncestry(const ConformingVertexAncestry &) = delete;
  ConformingVertexAncestry &operator=(const ConformingVertexAncestry &) = delete;

  void Clear();

  // Record the coarse edge-to-endpoint-identity map immediately before a conforming
  // refinement pass. Must be called while mesh is still the unrefined mesh.
  void Observe(const mfem::ParMesh &mesh,
               const std::vector<fem::singular::GlobalVertexId> &vertex_ids);

  // Assign identities to vertices appended by the refinement pass which followed the
  // matching Observe call. Returns false when no usable snapshot is available.
  bool Assign(const mfem::ParMesh &mesh,
              std::vector<fem::singular::GlobalVertexId> &vertex_ids);
};

// Rebuild exact rank-local source identities after in-place refinement. Conforming
// refinement preserves existing vertex indices and keys appended vertices by the
// persistent identities of their parent edge, using the snapshot in conforming_ancestry.
// Nonconforming refinement uses exact refinement-tree descriptors managed by
// vertex_identity. Element IDs are refreshed for diagnostics; singular incidence is
// reconstructed from vertex and edge identities rather than element numbering.
void UpdateSingularSourceEntityIds(
    const mfem::ParMesh &mesh,
    std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    std::vector<fem::singular::GlobalVertexId> &source_element_ids,
    NonconformingVertexIdentity &vertex_identity,
    ConformingVertexAncestry *conforming_ancestry = nullptr);

// Rebuild the immutable global feature blueprint after h-refinement. Existing physical
// singular features retain their exponent and straight-feature identity; selected
// boundary entities and any child feature segments are reconstructed from stable mesh
// vertex IDs. The solve is rebuilt from scratch, so no coefficient transfer is involved.
void RebuildRefinedSingularFeatures(
    const mfem::ParMesh &mesh, const std::vector<int> &boundary_attributes,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    fem::singular::FeatureTopology &features);
void RebuildRefinedSingularFeatures(
    const mfem::ParMesh &mesh, const std::vector<int> &boundary_attributes,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    fem::singular::TriangleFeatureTopology &features);

// Return a marker for every enriched element and its complete face-neighbor layer. On a
// nonconforming mesh this is the refinement closure used to keep all custom enrichment
// away from hanging interfaces. If supplied, "conforming" reports collectively whether
// every enriched element currently has only conforming or physical-boundary faces, and
// "repair" identifies local coarse elements for diagnostics and focused tests. Production
// nonconforming AMR protects the complete closure and fails closed if it becomes
// nonconforming.
mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    bool *conforming = nullptr, mfem::Array<int> *repair = nullptr);
mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::TriangleFeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    bool *conforming = nullptr, mfem::Array<int> *repair = nullptr);

// Immutable-source feature extraction and partition transport shared by the
// driven and eigenmode singular-element drivers.
class FullWaveSingularFeatures
{
private:
  fem::singular::FeatureTopology serial_sheet_features;
  fem::singular::FeatureTopology local_sheet_features;
  fem::singular::TriangleFeatureTopology serial_line_features;
  fem::singular::TriangleFeatureTopology local_line_features;
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids;
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  NonconformingVertexIdentity vertex_identity;
  ConformingVertexAncestry conforming_ancestry;
  int dimension = 0;

public:
  void Preprocess(const IoData &iodata, const std::unique_ptr<mfem::Mesh> &serial_mesh,
                  MPI_Comm comm);
  void ProcessPartitionedMesh(const IoData &iodata, const mfem::ParMesh &parallel_mesh,
                              const mesh::PartitionMetadata &metadata);
  void ProcessRefinedMesh(const IoData &iodata, const mfem::ParMesh &parallel_mesh);

  // Snapshot exact conforming-refinement ancestry immediately before the mesh is
  // refined. Nonconforming meshes ignore this because they carry an NCMesh tree.
  void ObserveRefinementAncestry(const mfem::ParMesh &parallel_mesh);

  // Report trace-connected component structure on the current mesh, for the three trace
  // policies, given the primary marks about to be refined. Diagnostic only.
  void ReportTraceComponents(const IoData &iodata, const mfem::ParMesh &parallel_mesh,
                             const mfem::Array<int> &primary_marks) const;

  const fem::singular::FeatureTopology *GetSheetFeatures() const;
  const fem::singular::TriangleFeatureTopology *GetLineFeatures() const;
  const std::vector<fem::singular::GlobalVertexId> *GetSourceVertexIds() const;
  mesh::PartitionMetadata GetSourceEntityMetadata() const;
  mfem::Array<int> GetRefinementProtection(const mfem::ParMesh &mesh,
                                           bool *conforming = nullptr,
                                           mfem::Array<int> *repair = nullptr) const;

  // Per-element marker of the elements carrying active singular basis functions.
  mfem::Array<int> GetEnrichedElements(const mfem::ParMesh &mesh) const;
};

// Bit mask over the four local faces of a tetrahedron (bit f set == local face f carries a
// nonzero singular trace). A tetrahedral local face f is the face opposite local vertex f,
// matching MFEM's tetrahedron face ordering convention only after permutation, so callers
// must map through the mesh face list rather than assuming an ordering.
//
// A face is INACTIVE only when every active singular H1 and ND basis of every feature in
// the element has identically zero relevant trace there: zero scalar potential for H1, and
// zero tangential component for ND. This follows from the basis algebra:
//
//   - node feature at local node n: the H1 potential carries the factor
//     (1 - rho^(nu-1)) with rho = 1 - lambda[n], which vanishes where lambda[n] = 0, and
//     the surviving ND value is parallel to grad(lambda[n]), which is normal to that face.
//     Hence the single face opposite n is inactive.
//   - edge feature on local edge (a, b): every H1 potential carries lambda[a]*lambda[b],
//     and the paired ND tangential trace vanishes likewise, so the two faces which do not
//     contain both a and b are inactive.
//   - the edge-rotational family is an H(curl) bubble whose tangential trace vanishes on
//     all four faces, so it never activates a face on its own.
//
// The reduction over features is a logical OR: any feature which activates a face wins.
// Essential-boundary exclusions (PEC, and impedance only where the exponent is
// constrained) are deliberately NOT applied here; that is a further refinement and its
// omission is conservative.
constexpr int SingularFaceMaskAllActive = 0xF;

int GetTetrahedronSingularFaceMask(const fem::singular::ElementFeatureIncidence &incidence);

// Which singular spaces a trace-conformity policy must protect. Maxwell solves need only
// the ND tangential trace; electrostatics needs only the H1 scalar trace; SHARED is the
// conservative union used when one closure serves both.
enum class SingularTracePolicy
{
  ND_ONLY,
  H1_ONLY,
  SHARED
};

// Diagnostic report on the trace-connected structure of the enriched patch, used to decide
// between up-front component closure, conforming refinement, and hanging constraints. It
// answers one question: how large is the transitive component of elements joined through
// faces that actually carry a singular trace?
struct SingularTraceComponentReport
{
  // Feature multiplicity per enriched element (index = count, value = elements).
  std::vector<long long> node_feature_histogram;
  std::vector<long long> edge_feature_histogram;

  // Active-face counts 0..4 over enriched elements, for the raw feature-incidence mask and
  // for the mask restricted to unconstrained DOFs under the requested policy.
  std::array<long long, 5> raw_active_face_histogram{};
  std::array<long long, 5> unconstrained_active_face_histogram{};

  long long enriched_elements = 0;
  // Faces where raw and unconstrained masks disagree, i.e. bases deactivated by PEC or by a
  // constrained impedance exponent.
  long long mask_differences = 0;
  // Active faces whose other side is not enriched. These bound the singular patch and are
  // worth inspecting: a trace should not normally be active there.
  long long active_faces_adjoining_non_enriched = 0;

  // Trace-connected components over enriched elements, largest first.
  std::vector<long long> component_sizes;
  long long largest_component = 0;
  // Components touched by the supplied primary marks, and their total size.
  long long marked_components = 0;
  long long marked_component_elements = 0;
};

// Compute the report on the current mesh without modifying it. primary_marks may be empty.
SingularTraceComponentReport MeasureSingularTraceComponents(
    const mfem::ParMesh &mesh, const fem::singular::FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &vertex_ids,
    const mfem::Array<int> &primary_marks, const std::set<int> &constrained_attributes,
    SingularTracePolicy policy, int singular_order);

void PrintSingularTraceComponentReport(const SingularTraceComponentReport &report,
                                       MPI_Comm comm, SingularTracePolicy policy);

struct SingularFullWaveEnergy
{
  double electric = 0.0;
  double magnetic = 0.0;
};

double SingularComplexQuadraticForm(MPI_Comm comm, const mfem::Operator &op,
                                    const ComplexVector &x);
SingularFullWaveEnergy MeasureSingularFullWaveEnergy(MPI_Comm comm,
                                                     const mfem::Operator &mass,
                                                     const mfem::Operator &stiffness,
                                                     const ComplexVector &electric_field,
                                                     std::complex<double> omega);

}  // namespace palace

#endif  // PALACE_DRIVERS_SINGULAR_SOLVER_HPP
