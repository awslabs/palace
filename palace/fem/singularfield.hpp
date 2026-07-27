// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARFIELD_HPP
#define PALACE_FEM_SINGULARFIELD_HPP

#include <complex>
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

struct TriangleH1FieldValue
{
  double potential = 0.0;
  Vector2 gradient{};
};

struct TriangleNDFieldValue
{
  Vector2 value{};
  double curl = 0.0;
};

struct TriangleBoundaryModeMagneticFieldValue
{
  Vector2 transverse_real{};
  Vector2 transverse_imag{};
  double normal_real = 0.0;
  double normal_imag = 0.0;
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

// A log-log fit of |E| along one interior triangle-sector ray approaching a
// singular PEC line tip. Canonical vertices are ordered as the tip followed by
// the two opposite vertices.
struct H1TipSlopeDiagnostic
{
  GlobalVertexId source_element;
  std::size_t feature;
  std::size_t selected_segment;
  std::array<GlobalVertexId, 3> canonical_vertices;
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

TriangleH1FieldValue
EvaluateElementTriangleH1Enrichment(const TriangleElementDofMap &element_dofs,
                                    const mfem::Vector &local_coefficients,
                                    const TriangleBarycentricPoint &lambda,
                                    const TriangleBarycentricGradients &grad_lambda);

TriangleNDFieldValue
EvaluateElementTriangleNDEnrichment(const TriangleElementDofMap &element_dofs,
                                    const mfem::Vector &local_coefficients,
                                    const TriangleBarycentricPoint &lambda,
                                    const TriangleBarycentricGradients &grad_lambda);

// Reconstruct B from the physical BoundaryMode electric field convention:
//   Bz = curl(Et) / (i omega),
//   Bt = -(kn / omega) (z x Et) + (grad(En) x z) / (i omega).
TriangleBoundaryModeMagneticFieldValue ReconstructTriangleBoundaryModeMagneticField(
    const Vector2 &et_real, const Vector2 &et_imag, double curl_et_real,
    double curl_et_imag, const Vector2 &grad_en_real, const Vector2 &grad_en_imag,
    std::complex<double> kn, double omega);

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

// Two-dimensional counterpart of EnrichedH1FieldEvaluator for the additive
// triangular basis. Combined vectors use the same process-local
// [standard, owned enrichment] ordering.
class TriangleEnrichedH1FieldEvaluator
{
private:
  const TriangleDofTopology &topology;
  const ParallelDofNumbering &numbering;
  mfem::ParFiniteElementSpace &fespace;
  mfem::ParGridFunction standard_field;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_prolongation;
  mfem::Vector local_enrichment;
  std::vector<double> h1_exponents;
  bool initialized;

  TriangleH1FieldValue EvaluateBarycentric(int element, const mfem::IntegrationPoint &point,
                                           const TriangleBarycentricPoint &lambda);

public:
  TriangleEnrichedH1FieldEvaluator(const TriangleDofTopology &topology,
                                   const ParallelDofNumbering &numbering,
                                   mfem::ParFiniteElementSpace &fespace);

  void SetFromTrueDofs(const mfem::Vector &combined_true_dofs);

  // The integration point must be strictly inside the reference triangle.
  TriangleH1FieldValue Evaluate(int element, const mfem::IntegrationPoint &point);

  // Integrate epsilon * |grad(V)|^2 over one physical triangle. A partitioned
  // node-Duffy rule is used when one triangle contains multiple singular tips.
  AdaptiveQuadratureResult
  IntegrateElementGradientEnergy(int element, double electric_coefficient,
                                 const AdaptiveAssemblyOptions &options);

  void ProjectToDiscontinuousGridFunctions(mfem::ParGridFunction &potential,
                                           mfem::ParGridFunction &electric_field);

  std::vector<H1CoefficientDiagnostic> GetOwnedCoefficientDiagnostics() const;

  std::vector<H1TipSlopeDiagnostic>
  FitTipSlopes(const TriangleFeatureTopology &features,
               const std::vector<GlobalVertexId> &source_vertex_ids,
               const std::vector<GlobalVertexId> &source_element_ids,
               const EdgeSlopeOptions &options = {});

  const mfem::Vector &GetLocalEnrichmentDofs() const { return local_enrichment; }
};

// Reconstruct a real standard-plus-singular triangular H(curl) field and its
// scalar curl. The caller uses one evaluator for each part of a complex field.
class TriangleEnrichedNDFieldEvaluator
{
private:
  const TriangleDofTopology &topology;
  const ParallelDofNumbering &numbering;
  mfem::ParFiniteElementSpace &fespace;
  mfem::ParGridFunction standard_field;
  std::unique_ptr<mfem::HypreParMatrix> enrichment_prolongation;
  mfem::Vector local_enrichment;
  std::vector<double> nd_exponents;
  bool initialized;

  TriangleNDFieldValue EvaluateBarycentric(int element, const mfem::IntegrationPoint &point,
                                           const TriangleBarycentricPoint &lambda);

public:
  TriangleEnrichedNDFieldEvaluator(const TriangleDofTopology &topology,
                                   const ParallelDofNumbering &numbering,
                                   mfem::ParFiniteElementSpace &fespace);

  void SetFromTrueDofs(const mfem::Vector &combined_true_dofs);

  // The integration point must be strictly inside the reference triangle.
  TriangleNDFieldValue Evaluate(int element, const mfem::IntegrationPoint &point);

  // Integrate coefficient * |field|^2 over one physical triangle.
  AdaptiveQuadratureResult
  IntegrateElementFieldEnergy(int element, double coefficient,
                              const AdaptiveAssemblyOptions &options);

  // Sample the vector field and scalar curl into discontinuous Gauss-Legendre
  // value spaces on the solve mesh.
  void ProjectToDiscontinuousGridFunctions(mfem::ParGridFunction &field,
                                           mfem::ParGridFunction &curl);

  std::vector<H1CoefficientDiagnostic> GetOwnedCoefficientDiagnostics() const;

  // Fit the magnitude of a complex transverse field represented by this
  // evaluator and `imaginary` along rays approaching every singular line tip.
  std::vector<H1TipSlopeDiagnostic>
  FitComplexTipSlopes(TriangleEnrichedNDFieldEvaluator &imaginary,
                      const TriangleFeatureTopology &features,
                      const std::vector<GlobalVertexId> &source_vertex_ids,
                      const std::vector<GlobalVertexId> &source_element_ids,
                      const EdgeSlopeOptions &options = {});

  const mfem::Vector &GetLocalEnrichmentDofs() const { return local_enrichment; }
};

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARFIELD_HPP
