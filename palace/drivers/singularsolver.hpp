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

// Rebuild exact rank-local source identities after in-place refinement. Conforming
// refinement preserves existing vertex indices and assigns decomposition-independent IDs
// to appended vertices in canonical physical-coordinate order. Distinct coincident new
// vertices are rejected because they require explicit refinement ancestry. Nonconforming
// refinement uses exact refinement-tree descriptors managed by vertex_identity. Element
// IDs are refreshed for diagnostics; singular incidence is reconstructed from vertex and
// edge identities rather than element numbering.
void UpdateSingularSourceEntityIds(
    const mfem::ParMesh &mesh,
    std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    std::vector<fem::singular::GlobalVertexId> &source_element_ids,
    NonconformingVertexIdentity &vertex_identity);

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
// "repair" marks the local coarse elements which must be refined to restore that property.
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
  int dimension = 0;

public:
  void Preprocess(const IoData &iodata, const std::unique_ptr<mfem::Mesh> &serial_mesh,
                  MPI_Comm comm);
  void ProcessPartitionedMesh(const IoData &iodata, const mfem::ParMesh &parallel_mesh,
                              const mesh::PartitionMetadata &metadata);
  void ProcessRefinedMesh(const IoData &iodata, const mfem::ParMesh &parallel_mesh);

  const fem::singular::FeatureTopology *GetSheetFeatures() const;
  const fem::singular::TriangleFeatureTopology *GetLineFeatures() const;
  const std::vector<fem::singular::GlobalVertexId> *GetSourceVertexIds() const;
  mesh::PartitionMetadata GetSourceEntityMetadata() const;
  mfem::Array<int> GetRefinementProtection(const mfem::ParMesh &mesh,
                                           bool *conforming = nullptr,
                                           mfem::Array<int> *repair = nullptr) const;
};

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
