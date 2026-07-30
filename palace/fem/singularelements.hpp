// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARELEMENTS_HPP
#define PALACE_FEM_SINGULARELEMENTS_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <variant>
#include <vector>

namespace palace
{

namespace fem
{

namespace singular
{

using Vector3 = std::array<double, 3>;
using BarycentricPoint = std::array<double, 4>;
using BarycentricGradients = std::array<Vector3, 4>;
using BarycentricPermutation = std::array<int, 4>;
using ReferenceCoefficientTensor = std::array<std::array<double, 3>, 3>;
using InterpolationIndices = std::array<int, 4>;

using Vector2 = std::array<double, 2>;
using TriangleBarycentricPoint = std::array<double, 3>;
using TriangleBarycentricGradients = std::array<Vector2, 3>;

struct WeightedSegmentQuadraturePoint
{
  double coordinate;
  double weight;
};

struct VectorBasisValue
{
  Vector3 value;
  Vector3 curl;
};

struct TriangleVectorBasisValue
{
  Vector2 value;
  double curl;
};

struct ScalarPolynomialValue
{
  double value;
  double derivative;
};

enum class FirstOrderBasisFamily
{
  STANDARD_H1_GRADIENT,
  STANDARD_NEDELEC,
  NODE_GRADIENT,
  NODE_ROTATIONAL,
  EDGE_GRADIENT,
  EDGE_ROTATIONAL
};

enum class HigherOrderBasisFamily
{
  NODE_GRADIENT,
  NODE_ROTATIONAL,
  EDGE_GRADIENT,
  EDGE_ROTATIONAL
};

enum class ReferenceQuadratureRule
{
  RECURSIVE_TETRAHEDRON,
  ADAPTIVE_TETRAHEDRON,
  NODE_DUFFY,
  EDGE_DUFFY,
  PARTITIONED_DUFFY
};

inline constexpr double MultiFeatureDuffyPartitionPower = 2.0;

struct FirstOrderBasis
{
  FirstOrderBasisFamily family;
  std::array<int, 4> indices;
  double nu;
};

struct HigherOrderBasis
{
  HigherOrderBasisFamily family;
  std::array<int, 4> nodes;
  InterpolationIndices interpolation_indices;
  int order;
  double nu;
};

using ReferenceBasis = std::variant<FirstOrderBasis, HigherOrderBasis>;

struct CanonicalFirstOrderBasisPair
{
  FirstOrderBasis row_basis;
  FirstOrderBasis column_basis;
  BarycentricPermutation input_to_canonical;
};

struct CanonicalReferenceBasisPair
{
  ReferenceBasis row_basis;
  ReferenceBasis column_basis;
  BarycentricPermutation input_to_canonical;
};

struct FirstOrderReferenceIntegral
{
  static constexpr int ConventionVersion = 1;
  static constexpr int StandardOrder = 1;
  static constexpr int SingularOrder = 1;

  FirstOrderBasis row_basis;
  FirstOrderBasis column_basis;
  BarycentricPermutation input_to_canonical;
  ReferenceQuadratureRule quadrature_rule;
  int quadrature_order;
  int subdivisions;
  double radial_power;
  ReferenceCoefficientTensor mass;
  ReferenceCoefficientTensor curl_curl;
};

struct ReferenceIntegral
{
  static constexpr int ConventionVersion = 1;

  ReferenceBasis row_basis;
  ReferenceBasis column_basis;
  BarycentricPermutation input_to_canonical;
  ReferenceQuadratureRule quadrature_rule;
  int quadrature_order;
  int subdivisions;
  double radial_power;
  ReferenceCoefficientTensor mass;
  ReferenceCoefficientTensor curl_curl;
};

struct AdaptiveQuadratureResult
{
  double value;
  double estimated_absolute_error;
  std::size_t leaf_count;
  int maximum_subdivision_depth;
  bool converged;
};

struct AdaptiveReferenceIntegral
{
  ReferenceIntegral integral;
  double absolute_tolerance;
  double relative_tolerance;
  ReferenceCoefficientTensor mass_estimated_absolute_error;
  ReferenceCoefficientTensor curl_curl_estimated_absolute_error;
  std::size_t total_leaf_count;
  int maximum_subdivision_depth;
  bool converged;
};

struct FirstOrderAdaptiveReferenceIntegral
{
  FirstOrderReferenceIntegral integral;
  double absolute_tolerance;
  double relative_tolerance;
  ReferenceCoefficientTensor mass_estimated_absolute_error;
  ReferenceCoefficientTensor curl_curl_estimated_absolute_error;
  std::size_t total_leaf_count;
  int maximum_subdivision_depth;
  bool converged;
};

// MFEM reference tetrahedron:
//   v0 = (0,0,0), v1 = (1,0,0), v2 = (0,1,0), v3 = (0,0,1).
BarycentricPoint ReferenceBarycentricPoint(const Vector3 &point);
const BarycentricGradients &ReferenceBarycentricGradients();

// MFEM reference triangle:
//   v0 = (0,0), v1 = (1,0), v2 = (0,1).
TriangleBarycentricPoint ReferenceTriangleBarycentricPoint(const Vector2 &point);
const TriangleBarycentricGradients &ReferenceTriangleBarycentricGradients();

// Silvester-Lagrange factors from equations (25) and (27), where
// grid_denominator is the superscript q. The derivative is with respect to the
// scalar coordinate.
ScalarPolynomialValue EvaluateSilvesterLagrange(int grid_denominator, int index,
                                                double coordinate);
ScalarPolynomialValue EvaluateShiftedSilvesterLagrange(int grid_denominator, int index,
                                                       double coordinate);

// Stable radial coordinates for a singular node or edge. These use
// complementary barycentric sums instead of subtracting nearly unit
// coordinates, which is essential at strongly graded quadrature points.
double NodeRadialCoordinate(const BarycentricPoint &lambda, int i);
double EdgeRadialCoordinate(const BarycentricPoint &lambda, int i, int j);
double TriangleNodeRadialCoordinate(const TriangleBarycentricPoint &lambda, int i);

// Gauss-Jacobi rule on [0, 1] for
//
//   integral_0^1 t^alpha (1-t)^beta f(t) dt.
//
// Both endpoint exponents must be greater than -1. The returned weights
// include the Jacobi weight and therefore sum to B(alpha + 1, beta + 1).
std::vector<WeightedSegmentQuadraturePoint>
BuildWeightedSegmentQuadrature(int order, double alpha, double beta);

// Lowest-order additive triangular bases from Graglia-Lombardi (2004),
// expressed in the barycentric and normalization conventions of Elkin et al.
// The scalar curl is d(value_y)/dx - d(value_x)/dy.
TriangleVectorBasisValue
EvaluateTriangleStandardEdge(const TriangleBarycentricPoint &lambda,
                             const TriangleBarycentricGradients &grad_lambda, int i, int j);
double EvaluateTriangleNodeGradientPotential(const TriangleBarycentricPoint &lambda, int i,
                                             int j, double nu);
TriangleVectorBasisValue
EvaluateTriangleNodeGradient(const TriangleBarycentricPoint &lambda,
                             const TriangleBarycentricGradients &grad_lambda, int i, int j,
                             double nu);
TriangleVectorBasisValue
EvaluateTriangleNodeRotational(const TriangleBarycentricPoint &lambda,
                               const TriangleBarycentricGradients &grad_lambda, int i,
                               int j, int k, double nu);

using TriangleQuadraturePointVisitor =
    std::function<void(const TriangleBarycentricPoint &lambda, double weight)>;
using TriangleReferenceIntegrand =
    std::function<double(const TriangleBarycentricPoint &lambda)>;

// Tensor-product Gauss rule under the singular-node Duffy map
//
//   lambda_i = 1-r, lambda_j = r(1-t), lambda_k = r t,
//   r = u^radial_power.
//
// The weights integrate over the MFEM reference triangle of area 1/2.
void ForEachReferenceTriangleNodeDuffyQuadraturePoint(
    int order, int singular_node, double radial_power,
    const TriangleQuadraturePointVisitor &visitor);
double IntegrateReferenceTriangleNodeDuffy(int order, int singular_node,
                                           double radial_power,
                                           const TriangleReferenceIntegrand &integrand);
std::size_t ReferenceTriangleNodeDuffyQuadraturePointCount(int order);

// Apply a simultaneous node relabeling where permutation[a] is the new index of
// old index a. Reference tables use these routines to map every feature tuple to
// a canonical tuple before quadrature.
BarycentricPoint ApplyBarycentricPermutation(const BarycentricPoint &lambda,
                                             const BarycentricPermutation &permutation);
BarycentricGradients ApplyBarycentricPermutation(const BarycentricGradients &grad_lambda,
                                                 const BarycentricPermutation &permutation);
FirstOrderBasis ApplyBarycentricPermutation(const FirstOrderBasis &basis,
                                            const BarycentricPermutation &permutation);
HigherOrderBasis ApplyBarycentricPermutation(const HigherOrderBasis &basis,
                                             const BarycentricPermutation &permutation);
ReferenceBasis ApplyBarycentricPermutation(const ReferenceBasis &basis,
                                           const BarycentricPermutation &permutation);
CanonicalFirstOrderBasisPair
CanonicalizeFirstOrderBasisPair(const FirstOrderBasis &row_basis,
                                const FirstOrderBasis &column_basis);
CanonicalReferenceBasisPair
CanonicalizeReferenceBasisPair(const ReferenceBasis &row_basis,
                               const ReferenceBasis &column_basis);

// Equation (5), including the paper's orientation convention:
//   N_ij = lambda_j grad(lambda_i) - lambda_i grad(lambda_j).
VectorBasisValue EvaluateStandardEdge(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i,
                                      int j);

// Scalar potential and gradient-singular vector basis from equation (6).
double EvaluateNodeGradientPotential(const BarycentricPoint &lambda, int i, int j,
                                     double nu);
VectorBasisValue EvaluateNodeGradient(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i, int j,
                                      double nu);

// Higher-order node-gradient scalar potential and vector basis from equation
// (28). The ordered nodes are i,j,k,l and the interpolation tuple is a,b,c,d.
double EvaluateHigherOrderNodeGradientPotential(
    const BarycentricPoint &lambda, int i, int j, int k, int l,
    const InterpolationIndices &interpolation_indices, int order, double nu);
VectorBasisValue EvaluateHigherOrderNodeGradient(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu);

// Node-singular rotational basis from equation (11).
VectorBasisValue EvaluateNodeRotational(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda, int i,
                                        int j, int k, double nu);
// Higher-order node-rotational basis from equation (30).
VectorBasisValue EvaluateHigherOrderNodeRotational(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu);

// Scalar potential and edge-gradient basis from equation (21).
double EvaluateEdgeGradientPotential(const BarycentricPoint &lambda, int i, int j, int k,
                                     double nu);
VectorBasisValue EvaluateEdgeGradient(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i, int j,
                                      int k, double nu);
// Higher-order edge-gradient scalar potential and vector basis from equation
// (32).
double EvaluateHigherOrderEdgeGradientPotential(
    const BarycentricPoint &lambda, int i, int j, int k, int l,
    const InterpolationIndices &interpolation_indices, int order, double nu);
VectorBasisValue EvaluateHigherOrderEdgeGradient(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu);

// Edge-singular rotational basis from equation (23). The four indices must be a
// permutation of the tetrahedron nodes, with ij the singular edge.
VectorBasisValue EvaluateEdgeRotational(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda, int i,
                                        int j, int k, int l, double nu);
// Higher-order edge-rotational basis from equation (33).
VectorBasisValue EvaluateHigherOrderEdgeRotational(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu);

// Enumerate the retained higher-order families after the dependency
// eliminations described below equations (28), (30), and (32). For node
// families, canonical_nodes[0] is the singular node and the remaining entries
// give a stable order for its neighbors. For edge families, canonical_nodes[0:2]
// is the oriented singular edge and the remaining entries give a stable order
// for the opposite nodes.
std::vector<HigherOrderBasis>
EnumerateHigherOrderNodeGradientBases(const std::array<int, 4> &canonical_nodes, int order,
                                      double nu);
std::vector<HigherOrderBasis>
EnumerateHigherOrderNodeRotationalBases(const std::array<int, 4> &canonical_nodes,
                                        int order, double nu);
std::vector<HigherOrderBasis>
EnumerateHigherOrderEdgeGradientBases(const std::array<int, 4> &canonical_nodes, int order,
                                      double nu);
std::vector<HigherOrderBasis>
EnumerateHigherOrderEdgeRotationalBases(const std::array<int, 4> &canonical_nodes,
                                        int order, double nu);
VectorBasisValue EvaluateHigherOrderBasis(const BarycentricPoint &lambda,
                                          const BarycentricGradients &grad_lambda,
                                          const HigherOrderBasis &basis);
// Value-only counterpart used by electric surface observables. This avoids
// evaluating curls for rotational bases when they are not consumed.
Vector3 EvaluateHigherOrderBasisValue(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda,
                                      const HigherOrderBasis &basis);
// Evaluate the scalar potential whose physical gradient is a retained
// higher-order gradient-family basis. Rotational families have no scalar H1
// potential and are rejected.
double EvaluateHigherOrderGradientPotential(const BarycentricPoint &lambda,
                                            const HigherOrderBasis &basis);
VectorBasisValue EvaluateReferenceBasis(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda,
                                        const ReferenceBasis &basis);

FirstOrderBasis MakeStandardH1Gradient(int i);
FirstOrderBasis MakeStandardNedelec(int i, int j);
FirstOrderBasis MakeNodeGradient(int i, int j, double nu);
FirstOrderBasis MakeNodeRotational(int i, int j, int k, double nu);
FirstOrderBasis MakeEdgeGradient(int i, int j, int k, double nu);
FirstOrderBasis MakeEdgeRotational(int i, int j, int k, int l, double nu);
VectorBasisValue EvaluateFirstOrderBasis(const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda,
                                         const FirstOrderBasis &basis);

// Stream an order-q tetrahedral Gaussian rule over a uniform recursive one-to-eight
// subdivision of the reference tetrahedron. All returned points are strictly interior,
// including for leaves touching a singular feature.
using QuadraturePointVisitor =
    std::function<void(const BarycentricPoint &lambda, double weight)>;
void ForEachReferenceTetrahedronQuadraturePoint(int order, int subdivisions,
                                                const QuadraturePointVisitor &visitor);
std::size_t ReferenceTetrahedronQuadraturePointCount(int order, int subdivisions);

// Independent tensor-product Gauss rules under feature-aligned Duffy maps with
// radial grading r = s^radial_power. A radial power of six polynomializes the
// supported exponents nu = 1/2 and 2/3. The order and radial power must leave
// every graded point representably interior in double precision. These are
// validation rules for integrals with one known singular node or edge, not
// replacements for the general recursive rule.
void ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
    int order, int singular_node, double radial_power,
    const QuadraturePointVisitor &visitor);
void ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
    int order, int singular_node_i, int singular_node_j, double radial_power,
    const QuadraturePointVisitor &visitor);

// Integrate with compensated long-double accumulation over the streamed rule.
// The adaptive rule compares each tetrahedron with its eight children, assigns
// absolute tolerance in proportion to reference volume, and uses local
// relative tolerances to select refinements. Accepted differences are summed
// into a global estimate; if that estimate misses the requested global
// absolute/relative tolerance, the local budgets are tightened and the
// integration is repeated. A result is converged only when the global estimate
// meets the original tolerance.
using ReferenceIntegrand = std::function<double(const BarycentricPoint &lambda)>;
double IntegrateReferenceTetrahedron(int order, int subdivisions,
                                     const ReferenceIntegrand &integrand);
AdaptiveQuadratureResult
IntegrateReferenceTetrahedronAdaptive(int order, double absolute_tolerance,
                                      double relative_tolerance, int max_subdivisions,
                                      const ReferenceIntegrand &integrand);
double IntegrateReferenceTetrahedronNodeDuffy(int order, int singular_node,
                                              double radial_power,
                                              const ReferenceIntegrand &integrand);
double IntegrateReferenceTetrahedronEdgeDuffy(int order, int singular_node_i,
                                              int singular_node_j, double radial_power,
                                              const ReferenceIntegrand &integrand);

// Equations (37)-(41), using grad(lambda_1..lambda_3) as the independent
// barycentric gradients and
//   e_1 = grad(lambda_2) x grad(lambda_3),
//   e_2 = grad(lambda_3) x grad(lambda_1),
//   e_3 = grad(lambda_1) x grad(lambda_2).
// The mass tensor also gives H1 diffusion entries whenever both descriptors are
// gradient families.
FirstOrderReferenceIntegral
ComputeFirstOrderReferenceIntegral(const FirstOrderBasis &row_basis,
                                   const FirstOrderBasis &column_basis, int order,
                                   int subdivisions);
FirstOrderAdaptiveReferenceIntegral ComputeFirstOrderAdaptiveReferenceIntegral(
    const FirstOrderBasis &row_basis, const FirstOrderBasis &column_basis, int order,
    double absolute_tolerance, double relative_tolerance, int max_subdivisions);
FirstOrderReferenceIntegral
ComputeFirstOrderNodeDuffyReferenceIntegral(const FirstOrderBasis &row_basis,
                                            const FirstOrderBasis &column_basis, int order,
                                            double radial_power);
FirstOrderReferenceIntegral
ComputeFirstOrderEdgeDuffyReferenceIntegral(const FirstOrderBasis &row_basis,
                                            const FirstOrderBasis &column_basis, int order,
                                            double radial_power);

// Generic reference tensors support first-order, higher-order, and mixed basis
// pairs. Every pair is canonically relabeled before quadrature so the
// deterministic recursive subdivision does not make equivalent element-local
// node numberings produce different tables.
ReferenceIntegral ComputeReferenceIntegral(const ReferenceBasis &row_basis,
                                           const ReferenceBasis &column_basis, int order,
                                           int subdivisions);
AdaptiveReferenceIntegral ComputeAdaptiveReferenceIntegral(
    const ReferenceBasis &row_basis, const ReferenceBasis &column_basis, int order,
    double absolute_tolerance, double relative_tolerance, int max_subdivisions);
ReferenceIntegral ComputeNodeDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                                    const ReferenceBasis &column_basis,
                                                    int order, double radial_power);
ReferenceIntegral ComputeEdgeDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                                    const ReferenceBasis &column_basis,
                                                    int order, double radial_power);
// Integrate a pair associated with two distinct singular features using two
// full-tetrahedron Duffy charts and the exact partition of unity
//
//   w_a = rho_b^q / (rho_a^q + rho_b^q),  w_b = 1 - w_a.
//
// Chart a is aligned with feature a and w_a cancels feature b; chart b does the
// converse. The fixed q = MultiFeatureDuffyPartitionPower is part of the
// quadrature convention. The routine rejects pairs without exactly two
// distinct node/edge singular features.
ReferenceIntegral
ComputePartitionedDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                         const ReferenceBasis &column_basis, int order,
                                         double radial_power);

// Contract a reference tensor with affine physical-element geometry. The
// Jacobian determinant maps the reference tetrahedron, whose volume is 1/6, to
// the physical tetrahedron. Scalar material coefficients are applied by the
// caller.
double ContractFirstOrderMass(const FirstOrderReferenceIntegral &integral,
                              const BarycentricGradients &grad_lambda,
                              double jacobian_determinant);
double ContractFirstOrderCurlCurl(const FirstOrderReferenceIntegral &integral,
                                  const BarycentricGradients &grad_lambda,
                                  double jacobian_determinant);
double ContractMass(const ReferenceIntegral &integral,
                    const BarycentricGradients &grad_lambda, double jacobian_determinant);
double ContractCurlCurl(const ReferenceIntegral &integral,
                        const BarycentricGradients &grad_lambda,
                        double jacobian_determinant);

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARELEMENTS_HPP
