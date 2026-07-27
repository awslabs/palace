// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARFIELD_HPP
#define PALACE_FEM_SINGULARFIELD_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>

#include "fem/singularassembly.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

struct H1FieldValue
{
  double potential = 0.0;
  Vector3 gradient{};
};

// One uniquely owned coefficient in the canonical singular H1 basis. The
// coefficient is a basis coordinate, not a mesh-independent physical edge
// amplitude.
struct H1CoefficientDiagnostic
{
  HYPRE_BigInt true_dof;
  DofKey key;
  double exponent;
  double coefficient;
};

struct EdgeSlopeOptions
{
  int sample_count = 9;
  double minimum_barycentric_radius = 1.0 / 262144.0;
  double maximum_barycentric_radius = 1.0 / 256.0;
};

// A log-log fit of |E| along one interior element-sector ray approaching a
// straight singular edge. The canonical vertices are ordered as singular-edge
// endpoints followed by the opposite-edge endpoints.
struct H1EdgeSlopeDiagnostic
{
  GlobalVertexId source_element;
  std::size_t feature;
  std::size_t segment;
  std::array<GlobalVertexId, 4> canonical_vertices;
  int sample_count;
  double exponent;
  double expected_slope;
  double fitted_slope;
  double r_squared;
  double minimum_distance;
  double maximum_distance;
  double field_norm_at_minimum_distance;
  double field_norm_at_maximum_distance;
  bool valid;
};

// Evaluate one element's scalar singular potential and physical gradient from
// rank-local canonical enrichment coefficients.
H1FieldValue EvaluateElementH1Enrichment(const ElementDofMap &element_dofs,
                                         const mfem::Vector &local_coefficients,
                                         const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda);

// Reconstruct and evaluate a combined standard-plus-singular H1 field. Input
// true vectors use the process-local ordering [standard, owned enrichment]
// produced by BuildParallelEnrichedOperator.
class EnrichedH1FieldEvaluator
{
private:
  const DofTopology &topology;
  const ParallelDofNumbering &numbering;
  mfem::ParFiniteElementSpace &fespace;
  mfem::ParGridFunction standard_field;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_prolongation;
  mfem::Vector local_enrichment;
  std::vector<double> h1_exponents;
  bool initialized;

public:
  EnrichedH1FieldEvaluator(const DofTopology &topology,
                           const ParallelDofNumbering &numbering,
                           mfem::ParFiniteElementSpace &fespace);

  void SetFromTrueDofs(const mfem::Vector &combined_true_dofs);

  // The integration point must be strictly inside the reference tetrahedron.
  H1FieldValue Evaluate(int element, const mfem::IntegrationPoint &point);

  // Integrate epsilon * |grad(V)|^2 over one physical element. The factor 1/2
  // used for stored electrostatic energy is deliberately not included.
  AdaptiveQuadratureResult
  IntegrateElementGradientEnergy(int element, double electric_coefficient,
                                 const AdaptiveAssemblyOptions &options);

  // Sample the combined potential and E = -grad(V) into discontinuous scalar
  // and three-component vector spaces on the same mesh. Gauss-Legendre L2
  // spaces keep every interpolation point strictly inside its tetrahedron.
  void ProjectToDiscontinuousGridFunctions(mfem::ParGridFunction &potential,
                                           mfem::ParGridFunction &electric_field);

  // Return exactly one record for every enrichment true DOF owned by this MPI
  // rank. Coefficients use the nondimensional voltage units of the solve.
  std::vector<H1CoefficientDiagnostic> GetOwnedCoefficientDiagnostics() const;

  std::vector<H1EdgeSlopeDiagnostic>
  FitEdgeSlopes(const FeatureTopology &features,
                const std::vector<GlobalVertexId> &source_vertex_ids,
                const std::vector<GlobalVertexId> &source_element_ids,
                const EdgeSlopeOptions &options = {});

  const mfem::Vector &GetLocalEnrichmentDofs() const { return local_enrichment; }
};

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARFIELD_HPP
