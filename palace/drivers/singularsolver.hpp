// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_SINGULAR_SOLVER_HPP
#define PALACE_DRIVERS_SINGULAR_SOLVER_HPP

#include <complex>
#include <memory>
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

nlohmann::json GetSingularSurfaceParticipationMetadata(const IoData &iodata);
nlohmann::json
GetSingularSurfaceIntegrabilityMetadata(const fem::singular::FeatureTopology &features);
nlohmann::json GetSingularSurfaceIntegrabilityMetadata(
    const fem::singular::TriangleFeatureTopology &features);

// Rebuild exact rank-local source identities after in-place refinement. Conforming
// refinement preserves existing vertex indices and assigns globally consistent IDs to
// appended vertices. Nonconforming refinement uses persistent NCMesh node IDs directly.
// Element IDs are refreshed for diagnostics; singular incidence is reconstructed from
// vertex and edge identities rather than element numbering.
void UpdateSingularSourceEntityIds(
    const mfem::ParMesh &mesh,
    std::vector<fem::singular::GlobalVertexId> &source_vertex_ids,
    std::vector<fem::singular::GlobalVertexId> &source_element_ids);

// Return a marker for every enriched element and its complete face-neighbor layer. The
// protected layer keeps custom enrichment away from hanging interfaces until a certified
// parent/child enrichment constraint operator is available.
mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::FeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids);
mfem::Array<int> BuildSingularRefinementProtection(
    const mfem::ParMesh &mesh, const fem::singular::TriangleFeatureTopology &features,
    const std::vector<fem::singular::GlobalVertexId> &source_vertex_ids);

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
  mfem::Array<int> GetRefinementProtection(const mfem::ParMesh &mesh) const;
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
