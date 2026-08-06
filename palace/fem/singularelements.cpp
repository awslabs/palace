// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularelements.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <Eigen/Eigenvalues>
#include <mfem.hpp>

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

constexpr double BarycentricTolerance = 1.0e-12;
constexpr int MaximumQuadratureOrder = 21;
constexpr int MaximumDuffyQuadratureOrder = 127;
constexpr int MaximumSubdivisionDepth = 8;
constexpr int MaximumAdaptiveSubdivisionDepth = 12;

using AffineTetrahedron = std::array<BarycentricPoint, 4>;

constexpr AffineTetrahedron ReferenceTetrahedron{
    BarycentricPoint{1.0, 0.0, 0.0, 0.0}, BarycentricPoint{0.0, 1.0, 0.0, 0.0},
    BarycentricPoint{0.0, 0.0, 1.0, 0.0}, BarycentricPoint{0.0, 0.0, 0.0, 1.0}};

struct CompensatedAccumulator
{
  long double sum = 0.0L;
  long double correction = 0.0L;

  void Add(long double term)
  {
    const long double updated = sum + term;
    if (std::abs(sum) >= std::abs(term))
    {
      correction += (sum - updated) + term;
    }
    else
    {
      correction += (term - updated) + sum;
    }
    sum = updated;
  }

  double Value() const { return static_cast<double>(sum + correction); }
};

void ValidateIndex(int i)
{
  if (i < 0 || i >= 4)
  {
    throw std::invalid_argument("Singular-element barycentric index must be in [0, 3]!");
  }
}

void ValidateTriangleIndex(int i)
{
  if (i < 0 || i >= 3)
  {
    throw std::invalid_argument(
        "Singular-element triangle barycentric index must be in [0, 2]!");
  }
}

template <typename... T>
void ValidateDistinctTriangleIndices(T... indices)
{
  const std::array<int, sizeof...(T)> values{indices...};
  for (std::size_t a = 0; a < values.size(); a++)
  {
    ValidateTriangleIndex(values[a]);
    for (std::size_t b = a + 1; b < values.size(); b++)
    {
      if (values[a] == values[b])
      {
        throw std::invalid_argument(
            "Singular-element triangle barycentric indices must be distinct!");
      }
    }
  }
}

template <typename... T>
void ValidateDistinctIndices(T... indices)
{
  const std::array<int, sizeof...(T)> values{indices...};
  for (std::size_t a = 0; a < values.size(); a++)
  {
    ValidateIndex(values[a]);
    for (std::size_t b = a + 1; b < values.size(); b++)
    {
      if (values[a] == values[b])
      {
        throw std::invalid_argument(
            "Singular-element barycentric indices must be distinct!");
      }
    }
  }
}

void ValidateExponent(double nu)
{
  if (!std::isfinite(nu) || nu <= 0.0 || nu >= 1.0)
  {
    throw std::invalid_argument(
        "Singular-element exponent must be finite and satisfy 0 < nu < 1!");
  }
}

void ValidateInterpolationPolynomial(int grid_denominator, int index, double coordinate,
                                     bool shifted)
{
  if (grid_denominator < 1)
  {
    throw std::invalid_argument(
        "Singular-element interpolation grid denominator must be positive!");
  }
  if (index < (shifted ? 1 : 0) || index > grid_denominator)
  {
    throw std::invalid_argument(
        "Singular-element interpolation index is incompatible with the grid!");
  }
  if (!std::isfinite(coordinate) || coordinate < -BarycentricTolerance ||
      coordinate > 1.0 + BarycentricTolerance)
  {
    throw std::invalid_argument(
        "Singular-element interpolation coordinate must be finite and in [0, 1]!");
  }
}

void ValidateSingularOrder(int order, int grid_offset)
{
  if (order < 1 || order > std::numeric_limits<int>::max() - grid_offset)
  {
    throw std::invalid_argument(
        "Singular-element interpolation order must be positive and representable!");
  }
}

void ValidateHigherOrderNodeGradientTuple(const InterpolationIndices &indices, int order)
{
  ValidateSingularOrder(order, 1);
  if (indices[0] < 1 || indices[0] > order || indices[1] < 1 || indices[1] > order ||
      indices[2] < 0 || indices[2] >= order || indices[3] < 0 || indices[3] >= order)
  {
    throw std::invalid_argument(
        "Equation-(28) interpolation indices must satisfy a,b in [1,s] and "
        "c,d in [0,s-1]!");
  }
  const long long sum =
      static_cast<long long>(indices[0]) + indices[1] + indices[2] + indices[3];
  if (sum != static_cast<long long>(order) + 1)
  {
    throw std::invalid_argument("Equation-(28) interpolation indices must sum to s+1!");
  }
}

void ValidateHigherOrderFaceTuple(const InterpolationIndices &indices, int order)
{
  ValidateSingularOrder(order, 2);
  if (indices[0] < 1 || indices[0] > order || indices[1] < 1 || indices[1] > order ||
      indices[2] < 1 || indices[2] > order || indices[3] < 0 || indices[3] >= order)
  {
    throw std::invalid_argument(
        "Equation-(30)/(32) interpolation indices must satisfy a,b,c in [1,s] "
        "and d in [0,s-1]!");
  }
  const long long sum =
      static_cast<long long>(indices[0]) + indices[1] + indices[2] + indices[3];
  if (sum != static_cast<long long>(order) + 2)
  {
    throw std::invalid_argument(
        "Equation-(30)/(32) interpolation indices must sum to s+2!");
  }
}

void ValidateHigherOrderVolumeTuple(const InterpolationIndices &indices, int order)
{
  ValidateSingularOrder(order, 3);
  if (indices[0] < 1 || indices[0] > order || indices[1] < 1 || indices[1] > order ||
      indices[2] < 1 || indices[2] > order || indices[3] < 1 || indices[3] > order)
  {
    throw std::invalid_argument(
        "Equation-(33) interpolation indices must satisfy a,b,c,d in [1,s]!");
  }
  const long long sum =
      static_cast<long long>(indices[0]) + indices[1] + indices[2] + indices[3];
  if (sum != static_cast<long long>(order) + 3)
  {
    throw std::invalid_argument("Equation-(33) interpolation indices must sum to s+3!");
  }
}

void ValidateHigherOrderBasis(const HigherOrderBasis &basis)
{
  ValidateDistinctIndices(basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3]);
  ValidateExponent(basis.nu);
  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
      ValidateHigherOrderNodeGradientTuple(basis.interpolation_indices, basis.order);
      return;
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
    case HigherOrderBasisFamily::EDGE_GRADIENT:
      ValidateHigherOrderFaceTuple(basis.interpolation_indices, basis.order);
      return;
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      ValidateHigherOrderVolumeTuple(basis.interpolation_indices, basis.order);
      return;
  }
  throw std::invalid_argument("Unknown higher-order singular-element basis family!");
}

int BasisIndexCount(FirstOrderBasisFamily family)
{
  switch (family)
  {
    case FirstOrderBasisFamily::STANDARD_H1_GRADIENT:
      return 1;
    case FirstOrderBasisFamily::STANDARD_NEDELEC:
    case FirstOrderBasisFamily::NODE_GRADIENT:
      return 2;
    case FirstOrderBasisFamily::NODE_ROTATIONAL:
    case FirstOrderBasisFamily::EDGE_GRADIENT:
      return 3;
    case FirstOrderBasisFamily::EDGE_ROTATIONAL:
      return 4;
  }
  throw std::invalid_argument("Unknown first-order singular-element basis family!");
}

void ValidateFirstOrderBasis(const FirstOrderBasis &basis)
{
  const int index_count = BasisIndexCount(basis.family);
  switch (index_count)
  {
    case 1:
      ValidateIndex(basis.indices[0]);
      break;
    case 2:
      ValidateDistinctIndices(basis.indices[0], basis.indices[1]);
      break;
    case 3:
      ValidateDistinctIndices(basis.indices[0], basis.indices[1], basis.indices[2]);
      break;
    case 4:
      ValidateDistinctIndices(basis.indices[0], basis.indices[1], basis.indices[2],
                              basis.indices[3]);
      break;
  }
  for (int i = index_count; i < 4; i++)
  {
    if (basis.indices[i] != -1)
    {
      throw std::invalid_argument(
          "Unused first-order singular-element basis indices must equal -1!");
    }
  }

  const bool standard = basis.family == FirstOrderBasisFamily::STANDARD_H1_GRADIENT ||
                        basis.family == FirstOrderBasisFamily::STANDARD_NEDELEC;
  if (standard)
  {
    if (basis.nu != 0.0)
    {
      throw std::invalid_argument(
          "Standard first-order basis descriptors must have zero exponent!");
    }
  }
  else
  {
    ValidateExponent(basis.nu);
  }
}

void ValidateReferenceBasis(const ReferenceBasis &basis)
{
  std::visit(
      [](const auto &entry)
      {
        using Basis = std::decay_t<decltype(entry)>;
        if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
        {
          ValidateFirstOrderBasis(entry);
        }
        else
        {
          ValidateHigherOrderBasis(entry);
        }
      },
      basis);
}

FirstOrderBasis ApplyBasisPermutation(const FirstOrderBasis &basis,
                                      const BarycentricPermutation &permutation)
{
  FirstOrderBasis result = basis;
  const int index_count = BasisIndexCount(basis.family);
  for (int i = 0; i < index_count; i++)
  {
    result.indices[i] = permutation[basis.indices[i]];
  }
  return result;
}

auto BasisCanonicalKey(const FirstOrderBasis &basis)
{
  return std::tuple{static_cast<int>(basis.family), basis.indices, basis.nu};
}

auto BasisPairCanonicalKey(const FirstOrderBasis &row_basis,
                           const FirstOrderBasis &column_basis,
                           const BarycentricPermutation &permutation)
{
  auto row_key = BasisCanonicalKey(row_basis);
  auto column_key = BasisCanonicalKey(column_basis);
  if (column_key < row_key)
  {
    std::swap(row_key, column_key);
  }
  return std::tuple{row_key, column_key, permutation};
}

auto ReferenceBasisCanonicalKey(const ReferenceBasis &basis)
{
  return std::visit(
      [](const auto &entry)
      {
        using Basis = std::decay_t<decltype(entry)>;
        if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
        {
          return std::tuple{0,
                            static_cast<int>(entry.family),
                            entry.indices,
                            InterpolationIndices{-1, -1, -1, -1},
                            1,
                            entry.nu};
        }
        else
        {
          return std::tuple{1,           static_cast<int>(entry.family),
                            entry.nodes, entry.interpolation_indices,
                            entry.order, entry.nu};
        }
      },
      basis);
}

auto ReferenceBasisPairCanonicalKey(const ReferenceBasis &row_basis,
                                    const ReferenceBasis &column_basis,
                                    const BarycentricPermutation &permutation)
{
  auto row_key = ReferenceBasisCanonicalKey(row_basis);
  auto column_key = ReferenceBasisCanonicalKey(column_basis);
  if (column_key < row_key)
  {
    std::swap(row_key, column_key);
  }
  return std::tuple{row_key, column_key, permutation};
}

bool IsNodeSingularFamily(FirstOrderBasisFamily family)
{
  return family == FirstOrderBasisFamily::NODE_GRADIENT ||
         family == FirstOrderBasisFamily::NODE_ROTATIONAL;
}

bool IsEdgeSingularFamily(FirstOrderBasisFamily family)
{
  return family == FirstOrderBasisFamily::EDGE_GRADIENT ||
         family == FirstOrderBasisFamily::EDGE_ROTATIONAL;
}

bool IsNodeSingularFamily(HigherOrderBasisFamily family)
{
  return family == HigherOrderBasisFamily::NODE_GRADIENT ||
         family == HigherOrderBasisFamily::NODE_ROTATIONAL;
}

bool IsEdgeSingularFamily(HigherOrderBasisFamily family)
{
  return family == HigherOrderBasisFamily::EDGE_GRADIENT ||
         family == HigherOrderBasisFamily::EDGE_ROTATIONAL;
}

enum class SingularFeatureKind
{
  NONE,
  NODE,
  EDGE
};

struct SingularFeature
{
  SingularFeatureKind kind = SingularFeatureKind::NONE;
  std::array<int, 2> nodes{-1, -1};
};

SingularFeature GetSingularFeature(const ReferenceBasis &basis)
{
  ValidateReferenceBasis(basis);
  return std::visit(
      [](const auto &entry)
      {
        using Basis = std::decay_t<decltype(entry)>;
        const auto node = [&]()
        {
          if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
          {
            return entry.indices[0];
          }
          else
          {
            return entry.nodes[0];
          }
        };
        const auto edge = [&]()
        {
          std::array<int, 2> result;
          if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
          {
            result = {entry.indices[0], entry.indices[1]};
          }
          else
          {
            result = {entry.nodes[0], entry.nodes[1]};
          }
          std::sort(result.begin(), result.end());
          return result;
        };
        if (IsNodeSingularFamily(entry.family))
        {
          return SingularFeature{SingularFeatureKind::NODE, {node(), -1}};
        }
        if (IsEdgeSingularFamily(entry.family))
        {
          return SingularFeature{SingularFeatureKind::EDGE, edge()};
        }
        return SingularFeature{};
      },
      basis);
}

bool SameFeature(const SingularFeature &a, const SingularFeature &b)
{
  return a.kind == b.kind && a.nodes == b.nodes;
}

double FeatureRadialCoordinate(const BarycentricPoint &lambda,
                               const SingularFeature &feature)
{
  switch (feature.kind)
  {
    case SingularFeatureKind::NODE:
      return NodeRadialCoordinate(lambda, feature.nodes[0]);
    case SingularFeatureKind::EDGE:
      return EdgeRadialCoordinate(lambda, feature.nodes[0], feature.nodes[1]);
    case SingularFeatureKind::NONE:
      break;
  }
  throw std::invalid_argument("Duffy quadrature requires a singular feature!");
}

double DuffyPartitionWeight(const BarycentricPoint &lambda,
                            const SingularFeature &aligned_feature,
                            const SingularFeature &other_feature)
{
  const double aligned_rho = FeatureRadialCoordinate(lambda, aligned_feature);
  const double other_rho = FeatureRadialCoordinate(lambda, other_feature);
  if (!(aligned_rho > 0.0) || !(other_rho > 0.0))
  {
    throw std::domain_error(
        "Partitioned Duffy quadrature requires strictly interior points!");
  }

  // Evaluate rho_other^q / (rho_aligned^q + rho_other^q) through a ratio so
  // strongly graded coordinates cannot underflow simultaneously.
  if (aligned_rho <= other_rho)
  {
    const double ratio = aligned_rho / other_rho;
    return 1.0 / (1.0 + std::pow(ratio, MultiFeatureDuffyPartitionPower));
  }
  const double ratio = other_rho / aligned_rho;
  const double power = std::pow(ratio, MultiFeatureDuffyPartitionPower);
  return power / (1.0 + power);
}

void ForEachFeatureDuffyQuadraturePoint(const SingularFeature &feature, int order,
                                        double radial_power,
                                        const QuadraturePointVisitor &visitor)
{
  switch (feature.kind)
  {
    case SingularFeatureKind::NODE:
      ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(order, feature.nodes[0],
                                                          radial_power, visitor);
      return;
    case SingularFeatureKind::EDGE:
      ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
          order, feature.nodes[0], feature.nodes[1], radial_power, visitor);
      return;
    case SingularFeatureKind::NONE:
      break;
  }
  throw std::invalid_argument("Duffy quadrature requires a singular feature!");
}

int GetCommonNodeFeature(const FirstOrderBasis &row_basis,
                         const FirstOrderBasis &column_basis)
{
  int singular_node = -1;
  for (const auto *basis : {&row_basis, &column_basis})
  {
    ValidateFirstOrderBasis(*basis);
    if (IsEdgeSingularFamily(basis->family))
    {
      throw std::invalid_argument(
          "Node-aligned Duffy quadrature does not support edge-singular basis "
          "functions!");
    }
    if (IsNodeSingularFamily(basis->family))
    {
      if (singular_node >= 0 && singular_node != basis->indices[0])
      {
        throw std::invalid_argument(
            "Node-aligned Duffy quadrature requires one common singular node!");
      }
      singular_node = basis->indices[0];
    }
  }
  if (singular_node < 0)
  {
    throw std::invalid_argument(
        "Node-aligned Duffy quadrature requires a node-singular basis function!");
  }
  return singular_node;
}

std::array<int, 2> GetCommonEdgeFeature(const FirstOrderBasis &row_basis,
                                        const FirstOrderBasis &column_basis)
{
  std::array<int, 2> singular_edge{-1, -1};
  for (const auto *basis : {&row_basis, &column_basis})
  {
    ValidateFirstOrderBasis(*basis);
    if (IsNodeSingularFamily(basis->family))
    {
      throw std::invalid_argument(
          "Edge-aligned Duffy quadrature does not support node-singular basis "
          "functions!");
    }
    if (IsEdgeSingularFamily(basis->family))
    {
      auto edge = std::array{basis->indices[0], basis->indices[1]};
      std::sort(edge.begin(), edge.end());
      if (singular_edge[0] >= 0 && singular_edge != edge)
      {
        throw std::invalid_argument(
            "Edge-aligned Duffy quadrature requires one common singular edge!");
      }
      singular_edge = edge;
    }
  }
  if (singular_edge[0] < 0)
  {
    throw std::invalid_argument(
        "Edge-aligned Duffy quadrature requires an edge-singular basis function!");
  }
  return singular_edge;
}

int GetCommonNodeFeature(const ReferenceBasis &row_basis,
                         const ReferenceBasis &column_basis)
{
  int singular_node = -1;
  for (const auto *basis : {&row_basis, &column_basis})
  {
    ValidateReferenceBasis(*basis);
    std::visit(
        [&](const auto &entry)
        {
          if (IsEdgeSingularFamily(entry.family))
          {
            throw std::invalid_argument(
                "Node-aligned Duffy quadrature does not support edge-singular basis "
                "functions!");
          }
          if (IsNodeSingularFamily(entry.family))
          {
            const int node = [&]()
            {
              using Basis = std::decay_t<decltype(entry)>;
              if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
              {
                return entry.indices[0];
              }
              else
              {
                return entry.nodes[0];
              }
            }();
            if (singular_node >= 0 && singular_node != node)
            {
              throw std::invalid_argument(
                  "Node-aligned Duffy quadrature requires one common singular node!");
            }
            singular_node = node;
          }
        },
        *basis);
  }
  if (singular_node < 0)
  {
    throw std::invalid_argument(
        "Node-aligned Duffy quadrature requires a node-singular basis function!");
  }
  return singular_node;
}

std::array<int, 2> GetCommonEdgeFeature(const ReferenceBasis &row_basis,
                                        const ReferenceBasis &column_basis)
{
  std::array<int, 2> singular_edge{-1, -1};
  for (const auto *basis : {&row_basis, &column_basis})
  {
    ValidateReferenceBasis(*basis);
    std::visit(
        [&](const auto &entry)
        {
          if (IsNodeSingularFamily(entry.family))
          {
            throw std::invalid_argument(
                "Edge-aligned Duffy quadrature does not support node-singular basis "
                "functions!");
          }
          if (IsEdgeSingularFamily(entry.family))
          {
            auto edge = [&]()
            {
              using Basis = std::decay_t<decltype(entry)>;
              if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
              {
                return std::array{entry.indices[0], entry.indices[1]};
              }
              else
              {
                return std::array{entry.nodes[0], entry.nodes[1]};
              }
            }();
            std::sort(edge.begin(), edge.end());
            if (singular_edge[0] >= 0 && singular_edge != edge)
            {
              throw std::invalid_argument(
                  "Edge-aligned Duffy quadrature requires one common singular edge!");
            }
            singular_edge = edge;
          }
        },
        *basis);
  }
  if (singular_edge[0] < 0)
  {
    throw std::invalid_argument(
        "Edge-aligned Duffy quadrature requires an edge-singular basis function!");
  }
  return singular_edge;
}

void ValidateBarycentricPoint(const BarycentricPoint &lambda)
{
  double sum = 0.0;
  for (double value : lambda)
  {
    if (!std::isfinite(value) || value < -BarycentricTolerance ||
        value > 1.0 + BarycentricTolerance)
    {
      throw std::invalid_argument(
          "Singular-element barycentric coordinates must describe a tetrahedron point!");
    }
    sum += value;
  }
  if (std::abs(sum - 1.0) > BarycentricTolerance)
  {
    throw std::invalid_argument(
        "Singular-element barycentric coordinates must sum to one!");
  }
}

void ValidateBarycentricGradients(const BarycentricGradients &grad_lambda)
{
  Vector3 sum{0.0, 0.0, 0.0};
  double scale = 1.0;
  for (const auto &gradient : grad_lambda)
  {
    for (int d = 0; d < 3; d++)
    {
      if (!std::isfinite(gradient[d]))
      {
        throw std::invalid_argument(
            "Singular-element barycentric gradients must be finite!");
      }
      sum[d] += gradient[d];
      scale = std::max(scale, std::abs(gradient[d]));
    }
  }
  for (double value : sum)
  {
    if (std::abs(value) > BarycentricTolerance * scale)
    {
      throw std::invalid_argument(
          "Singular-element barycentric gradients must sum to zero!");
    }
  }
}

void ValidateTriangleBarycentricPoint(const TriangleBarycentricPoint &lambda)
{
  double sum = 0.0;
  for (double value : lambda)
  {
    if (!std::isfinite(value) || value < -BarycentricTolerance ||
        value > 1.0 + BarycentricTolerance)
    {
      throw std::invalid_argument(
          "Singular-element barycentric coordinates must describe a triangle point!");
    }
    sum += value;
  }
  if (std::abs(sum - 1.0) > BarycentricTolerance)
  {
    throw std::invalid_argument(
        "Singular-element triangle barycentric coordinates must sum to one!");
  }
}

void ValidateTriangleBarycentricGradients(const TriangleBarycentricGradients &grad_lambda)
{
  Vector2 sum{0.0, 0.0};
  double scale = 1.0;
  for (const auto &gradient : grad_lambda)
  {
    for (int d = 0; d < 2; d++)
    {
      if (!std::isfinite(gradient[d]))
      {
        throw std::invalid_argument(
            "Singular-element triangle barycentric gradients must be finite!");
      }
      sum[d] += gradient[d];
      scale = std::max(scale, std::abs(gradient[d]));
    }
  }
  for (double value : sum)
  {
    if (std::abs(value) > BarycentricTolerance * scale)
    {
      throw std::invalid_argument(
          "Singular-element triangle barycentric gradients must sum to zero!");
    }
  }
}

void ValidateTriangleEvaluation(const TriangleBarycentricPoint &lambda,
                                const TriangleBarycentricGradients &grad_lambda, double nu)
{
  ValidateTriangleBarycentricPoint(lambda);
  ValidateTriangleBarycentricGradients(grad_lambda);
  ValidateExponent(nu);
}

void ValidateEvaluation(const BarycentricPoint &lambda,
                        const BarycentricGradients &grad_lambda, double nu)
{
  ValidateBarycentricPoint(lambda);
  ValidateBarycentricGradients(grad_lambda);
  ValidateExponent(nu);
}

double PositivePower(double base, double exponent)
{
  return std::exp(exponent * std::log(base));
}

double OneMinusPositivePower(double base, double exponent)
{
  return -std::expm1(exponent * std::log(base));
}

double PositivePowerMinusOne(double base, double exponent)
{
  return std::expm1(exponent * std::log(base));
}

void ValidateRadialProxy(double rho)
{
  if (!(rho > 0.0))
  {
    throw std::domain_error(
        "Singular-element basis cannot be evaluated on its singular feature!");
  }
}

Vector3 Scale(double scale, const Vector3 &x)
{
  return {scale * x[0], scale * x[1], scale * x[2]};
}

void Add(double scale, const Vector3 &x, Vector3 &y)
{
  for (int d = 0; d < 3; d++)
  {
    y[d] += scale * x[d];
  }
}

Vector3 Cross(const Vector3 &x, const Vector3 &y)
{
  return {x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
}

Vector2 Scale(double scale, const Vector2 &x)
{
  return {scale * x[0], scale * x[1]};
}

void Add(double scale, const Vector2 &x, Vector2 &y)
{
  for (int d = 0; d < 2; d++)
  {
    y[d] += scale * x[d];
  }
}

double Cross(const Vector2 &x, const Vector2 &y)
{
  return x[0] * y[1] - x[1] * y[0];
}

std::array<ScalarPolynomialValue, 4>
EvaluateHigherOrderNodeGradientPolynomials(const BarycentricPoint &lambda,
                                           const std::array<int, 4> &nodes,
                                           const InterpolationIndices &indices, int order)
{
  ValidateHigherOrderNodeGradientTuple(indices, order);
  const int grid_denominator = order + 1;
  return {EvaluateShiftedSilvesterLagrange(grid_denominator, indices[0], lambda[nodes[0]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[1], lambda[nodes[1]]),
          EvaluateSilvesterLagrange(grid_denominator, indices[2], lambda[nodes[2]]),
          EvaluateSilvesterLagrange(grid_denominator, indices[3], lambda[nodes[3]])};
}

std::array<ScalarPolynomialValue, 4>
EvaluateHigherOrderFacePolynomials(const BarycentricPoint &lambda,
                                   const std::array<int, 4> &nodes,
                                   const InterpolationIndices &indices, int order)
{
  ValidateHigherOrderFaceTuple(indices, order);
  const int grid_denominator = order + 2;
  return {EvaluateShiftedSilvesterLagrange(grid_denominator, indices[0], lambda[nodes[0]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[1], lambda[nodes[1]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[2], lambda[nodes[2]]),
          EvaluateSilvesterLagrange(grid_denominator, indices[3], lambda[nodes[3]])};
}

std::array<ScalarPolynomialValue, 4>
EvaluateHigherOrderVolumePolynomials(const BarycentricPoint &lambda,
                                     const std::array<int, 4> &nodes,
                                     const InterpolationIndices &indices, int order)
{
  ValidateHigherOrderVolumeTuple(indices, order);
  const int grid_denominator = order + 3;
  return {EvaluateShiftedSilvesterLagrange(grid_denominator, indices[0], lambda[nodes[0]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[1], lambda[nodes[1]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[2], lambda[nodes[2]]),
          EvaluateShiftedSilvesterLagrange(grid_denominator, indices[3], lambda[nodes[3]])};
}

double PolynomialProduct(const std::array<ScalarPolynomialValue, 4> &polynomials)
{
  double product = 1.0;
  for (const auto &polynomial : polynomials)
  {
    product *= polynomial.value;
  }
  return product;
}

Vector3 PolynomialProductGradient(const std::array<ScalarPolynomialValue, 4> &polynomials,
                                  const BarycentricGradients &grad_lambda,
                                  const std::array<int, 4> &nodes)
{
  Vector3 gradient{};
  for (int a = 0; a < 4; a++)
  {
    double coefficient = polynomials[a].derivative;
    for (int b = 0; b < 4; b++)
    {
      if (b != a)
      {
        coefficient *= polynomials[b].value;
      }
    }
    Add(coefficient, grad_lambda[nodes[a]], gradient);
  }
  return gradient;
}

double Dot(const Vector3 &x, const Vector3 &y)
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

VectorBasisValue StandardEdge(const BarycentricPoint &lambda,
                              const BarycentricGradients &grad_lambda, int i, int j)
{
  VectorBasisValue result{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Add(lambda[j], grad_lambda[i], result.value);
  Add(-lambda[i], grad_lambda[j], result.value);
  result.curl = Scale(-2.0, Cross(grad_lambda[i], grad_lambda[j]));
  return result;
}

void ValidateQuadratureParameters(int order, int subdivisions)
{
  if (order < 1 || order > MaximumQuadratureOrder)
  {
    throw std::invalid_argument("Singular-element tetrahedron quadrature order must be "
                                "in [1, 21]!");
  }
  if (subdivisions < 0 || subdivisions > MaximumSubdivisionDepth)
  {
    throw std::invalid_argument(
        "Singular-element tetrahedron subdivision depth must be in [0, 8]!");
  }
}

void ValidateQuadratureVisitor(const QuadraturePointVisitor &visitor)
{
  if (!visitor)
  {
    throw std::invalid_argument(
        "Singular-element tetrahedron quadrature visitor must be callable!");
  }
}

void ValidateAdaptiveQuadratureParameters(double absolute_tolerance,
                                          double relative_tolerance, int max_subdivisions)
{
  if (!std::isfinite(absolute_tolerance) || absolute_tolerance < 0.0 ||
      !std::isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
      !(absolute_tolerance > 0.0 || relative_tolerance > 0.0))
  {
    throw std::invalid_argument(
        "Singular-element adaptive quadrature tolerances must be finite and "
        "nonnegative, with at least one positive!");
  }
  if (max_subdivisions < 1 || max_subdivisions > MaximumAdaptiveSubdivisionDepth)
  {
    throw std::invalid_argument(
        "Singular-element adaptive quadrature maximum subdivision depth must be "
        "in [1, 12]!");
  }
}

void ValidateDuffyQuadratureParameters(int order, double radial_power)
{
  if (order < 1 || order > MaximumDuffyQuadratureOrder)
  {
    throw std::invalid_argument(
        "Singular-element Duffy quadrature order must be in [1, 127]!");
  }
  if (!std::isfinite(radial_power) || !(radial_power > 0.0))
  {
    throw std::invalid_argument(
        "Singular-element Duffy radial power must be finite and positive!");
  }
}

const mfem::IntegrationRule &GetDuffyQuadratureRule(int order, double radial_power)
{
  ValidateDuffyQuadratureParameters(order, radial_power);
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, order);
  for (int q = 0; q < rule.GetNPoints(); q++)
  {
    const double parameter = rule.IntPoint(q).x;
    const double radial_coordinate = std::pow(parameter, radial_power);
    if (!(radial_coordinate > 0.0) || !(radial_coordinate < 1.0))
    {
      throw std::invalid_argument(
          "Singular-element Duffy order and radial power produce a graded point "
          "which is not representable in double precision!");
    }
  }
  return rule;
}

void ValidateNodeDuffyQuadratureRule(const mfem::IntegrationRule &rule, double radial_power)
{
  for (int q = 0; q < rule.GetNPoints(); q++)
  {
    const double radial_coordinate = std::pow(rule.IntPoint(q).x, radial_power);
    if (!(1.0 - radial_coordinate < 1.0))
    {
      throw std::invalid_argument(
          "Singular-element node Duffy order and radial power produce a graded "
          "point which is not representably interior in double precision!");
    }
  }
}

void ValidateFiniteBasisValue(const VectorBasisValue &basis)
{
  for (double value : basis.value)
  {
    if (!std::isfinite(value))
    {
      throw std::domain_error("Singular-element basis returned a non-finite vector value!");
    }
  }
  for (double value : basis.curl)
  {
    if (!std::isfinite(value))
    {
      throw std::domain_error("Singular-element basis returned a non-finite curl!");
    }
  }
}

void ValidateAffineGeometry(const BarycentricGradients &grad_lambda,
                            double jacobian_determinant)
{
  ValidateBarycentricGradients(grad_lambda);
  if (!std::isfinite(jacobian_determinant) || !(jacobian_determinant > 0.0))
  {
    throw std::invalid_argument(
        "Singular-element affine Jacobian determinant must be finite and positive!");
  }

  const double inverse_jacobian_determinant =
      std::abs(Dot(grad_lambda[1], Cross(grad_lambda[2], grad_lambda[3])));
  if (!(inverse_jacobian_determinant > 0.0))
  {
    throw std::invalid_argument(
        "Singular-element barycentric gradients must be linearly independent!");
  }
  const double determinant_product = jacobian_determinant * inverse_jacobian_determinant;
  if (std::abs(determinant_product - 1.0) > 1.0e-10 * std::max(1.0, determinant_product))
  {
    throw std::invalid_argument(
        "Singular-element barycentric gradients and Jacobian determinant are "
        "inconsistent!");
  }
}

BarycentricPoint Midpoint(const BarycentricPoint &a, const BarycentricPoint &b)
{
  BarycentricPoint midpoint;
  for (int i = 0; i < 4; i++)
  {
    midpoint[i] = 0.5 * (a[i] + b[i]);
  }
  return midpoint;
}

std::array<AffineTetrahedron, 8> Subdivide(const AffineTetrahedron &tetrahedron)
{
  const auto m01 = Midpoint(tetrahedron[0], tetrahedron[1]);
  const auto m02 = Midpoint(tetrahedron[0], tetrahedron[2]);
  const auto m03 = Midpoint(tetrahedron[0], tetrahedron[3]);
  const auto m12 = Midpoint(tetrahedron[1], tetrahedron[2]);
  const auto m13 = Midpoint(tetrahedron[1], tetrahedron[3]);
  const auto m23 = Midpoint(tetrahedron[2], tetrahedron[3]);

  // The four central tetrahedra use the m01-m23 diagonal of the midpoint
  // octahedron. This is a deterministic reference-table convention.
  return {AffineTetrahedron{tetrahedron[0], m01, m02, m03},
          AffineTetrahedron{tetrahedron[1], m01, m12, m13},
          AffineTetrahedron{tetrahedron[2], m02, m12, m23},
          AffineTetrahedron{tetrahedron[3], m03, m13, m23},
          AffineTetrahedron{m01, m02, m03, m23},
          AffineTetrahedron{m01, m02, m12, m23},
          AffineTetrahedron{m01, m03, m13, m23},
          AffineTetrahedron{m01, m12, m13, m23}};
}

double Determinant(const AffineTetrahedron &tetrahedron)
{
  Vector3 edge_1, edge_2, edge_3;
  for (int d = 0; d < 3; d++)
  {
    edge_1[d] = tetrahedron[1][d + 1] - tetrahedron[0][d + 1];
    edge_2[d] = tetrahedron[2][d + 1] - tetrahedron[0][d + 1];
    edge_3[d] = tetrahedron[3][d + 1] - tetrahedron[0][d + 1];
  }
  return Dot(Cross(edge_1, edge_2), edge_3);
}

BarycentricPoint MapToTetrahedron(const AffineTetrahedron &tetrahedron,
                                  const mfem::IntegrationPoint &point)
{
  const BarycentricPoint reference{1.0 - point.x - point.y - point.z, point.x, point.y,
                                   point.z};
  BarycentricPoint mapped{};
  for (int vertex = 0; vertex < 4; vertex++)
  {
    for (int i = 0; i < 4; i++)
    {
      mapped[i] += reference[vertex] * tetrahedron[vertex][i];
    }
  }
  return mapped;
}

void VisitQuadraturePoints(const AffineTetrahedron &tetrahedron, int subdivisions,
                           const mfem::IntegrationRule &rule,
                           const QuadraturePointVisitor &visitor)
{
  if (subdivisions > 0)
  {
    for (const auto &child : Subdivide(tetrahedron))
    {
      VisitQuadraturePoints(child, subdivisions - 1, rule, visitor);
    }
    return;
  }

  const double determinant = std::abs(Determinant(tetrahedron));
  if (!(determinant > 0.0))
  {
    throw std::runtime_error(
        "Singular-element tetrahedron subdivision produced a degenerate leaf!");
  }
  for (int q = 0; q < rule.GetNPoints(); q++)
  {
    const auto &point = rule.IntPoint(q);
    const auto lambda = MapToTetrahedron(tetrahedron, point);
    ValidateBarycentricPoint(lambda);
    for (double value : lambda)
    {
      if (!(value > 0.0))
      {
        throw std::runtime_error(
            "Singular-element volume quadrature point is not strictly interior!");
      }
    }
    visitor(lambda, determinant * point.weight);
  }
}

template <typename Q>
double IntegrateReferenceQuadrature(Q &&quadrature, const ReferenceIntegrand &integrand)
{
  if (!integrand)
  {
    throw std::invalid_argument("Singular-element reference integrand must be callable!");
  }

  CompensatedAccumulator accumulator;
  quadrature(
      [&](const BarycentricPoint &lambda, double weight)
      {
        const double value = integrand(lambda);
        if (!std::isfinite(value))
        {
          throw std::domain_error(
              "Singular-element reference integrand returned a non-finite value!");
        }
        accumulator.Add(static_cast<long double>(weight) * static_cast<long double>(value));
      });
  return accumulator.Value();
}

double IntegrateAffineTetrahedron(const AffineTetrahedron &tetrahedron,
                                  const mfem::IntegrationRule &rule,
                                  const ReferenceIntegrand &integrand)
{
  CompensatedAccumulator accumulator;
  VisitQuadraturePoints(
      tetrahedron, 0, rule,
      [&](const BarycentricPoint &lambda, double weight)
      {
        const double value = integrand(lambda);
        if (!std::isfinite(value))
        {
          throw std::domain_error(
              "Singular-element reference integrand returned a non-finite value!");
        }
        accumulator.Add(static_cast<long double>(weight) * static_cast<long double>(value));
      });
  return accumulator.Value();
}

std::vector<double> IntegrateAffineTetrahedron(const AffineTetrahedron &tetrahedron,
                                               const mfem::IntegrationRule &rule,
                                               std::size_t number_components,
                                               const ReferenceVectorIntegrand &integrand)
{
  std::vector<CompensatedAccumulator> accumulators(number_components);
  std::vector<double> values(number_components);
  VisitQuadraturePoints(
      tetrahedron, 0, rule,
      [&](const BarycentricPoint &lambda, double weight)
      {
        std::fill(values.begin(), values.end(), std::numeric_limits<double>::quiet_NaN());
        integrand(lambda, values);
        if (values.size() != number_components)
        {
          throw std::domain_error(
              "Singular-element vector integrand changed its output size!");
        }
        for (std::size_t component = 0; component < number_components; component++)
        {
          if (!std::isfinite(values[component]))
          {
            throw std::domain_error(
                "Singular-element vector integrand returned a non-finite value!");
          }
          accumulators[component].Add(static_cast<long double>(weight) *
                                      static_cast<long double>(values[component]));
        }
      });

  std::vector<double> result(number_components);
  for (std::size_t component = 0; component < number_components; component++)
  {
    result[component] = accumulators[component].Value();
  }
  return result;
}

AdaptiveQuadratureResult IntegrateAffineTetrahedronAdaptive(
    const AffineTetrahedron &tetrahedron, const mfem::IntegrationRule &rule,
    const ReferenceIntegrand &integrand, double coarse_value, double absolute_tolerance,
    double relative_tolerance, int depth, int max_subdivisions)
{
  const auto children = Subdivide(tetrahedron);
  std::array<double, 8> child_values;
  CompensatedAccumulator fine_accumulator;
  for (std::size_t child = 0; child < children.size(); child++)
  {
    child_values[child] = IntegrateAffineTetrahedron(children[child], rule, integrand);
    fine_accumulator.Add(static_cast<long double>(child_values[child]));
  }
  const double fine_value = fine_accumulator.Value();
  const double estimated_error = std::abs(fine_value - coarse_value);
  const double volume_fraction = std::abs(Determinant(tetrahedron));
  const double local_tolerance =
      absolute_tolerance * volume_fraction + relative_tolerance * std::abs(fine_value);
  if (estimated_error <= local_tolerance || depth + 1 == max_subdivisions)
  {
    return {fine_value, estimated_error, children.size(), depth + 1,
            estimated_error <= local_tolerance};
  }

  CompensatedAccumulator value_accumulator;
  CompensatedAccumulator error_accumulator;
  std::size_t leaf_count = 0;
  int maximum_depth = depth;
  bool converged = true;
  for (std::size_t child = 0; child < children.size(); child++)
  {
    const auto result = IntegrateAffineTetrahedronAdaptive(
        children[child], rule, integrand, child_values[child], absolute_tolerance,
        relative_tolerance, depth + 1, max_subdivisions);
    value_accumulator.Add(static_cast<long double>(result.value));
    error_accumulator.Add(static_cast<long double>(result.estimated_absolute_error));
    if (leaf_count > std::numeric_limits<std::size_t>::max() - result.leaf_count)
    {
      throw std::overflow_error(
          "Singular-element adaptive quadrature leaf count overflow!");
    }
    leaf_count += result.leaf_count;
    maximum_depth = std::max(maximum_depth, result.maximum_subdivision_depth);
    converged = converged && result.converged;
  }
  return {value_accumulator.Value(), error_accumulator.Value(), leaf_count, maximum_depth,
          converged};
}

AdaptiveVectorQuadratureResult IntegrateAffineTetrahedronAdaptive(
    const AffineTetrahedron &tetrahedron, const mfem::IntegrationRule &rule,
    std::size_t number_components, const ReferenceVectorIntegrand &integrand,
    const std::vector<double> &coarse_value, double absolute_tolerance,
    double relative_tolerance, int depth, int max_subdivisions)
{
  MFEM_ASSERT(coarse_value.size() == number_components,
              "Invalid coarse vector quadrature value!");
  const auto children = Subdivide(tetrahedron);
  std::array<std::vector<double>, 8> child_values;
  std::vector<CompensatedAccumulator> fine_accumulators(number_components);
  for (std::size_t child = 0; child < children.size(); child++)
  {
    child_values[child] =
        IntegrateAffineTetrahedron(children[child], rule, number_components, integrand);
    for (std::size_t component = 0; component < number_components; component++)
    {
      fine_accumulators[component].Add(
          static_cast<long double>(child_values[child][component]));
    }
  }

  std::vector<double> fine_value(number_components);
  std::vector<double> estimated_error(number_components);
  const double volume_fraction = std::abs(Determinant(tetrahedron));
  bool locally_converged = true;
  for (std::size_t component = 0; component < number_components; component++)
  {
    fine_value[component] = fine_accumulators[component].Value();
    estimated_error[component] = std::abs(fine_value[component] - coarse_value[component]);
    const double local_tolerance = absolute_tolerance * volume_fraction +
                                   relative_tolerance * std::abs(fine_value[component]);
    locally_converged = locally_converged && estimated_error[component] <= local_tolerance;
  }
  if (locally_converged || depth + 1 == max_subdivisions)
  {
    return {std::move(fine_value), std::move(estimated_error), children.size(), depth + 1,
            locally_converged};
  }

  std::vector<CompensatedAccumulator> value_accumulators(number_components);
  std::vector<CompensatedAccumulator> error_accumulators(number_components);
  std::size_t leaf_count = 0;
  int maximum_depth = depth;
  bool converged = true;
  for (std::size_t child = 0; child < children.size(); child++)
  {
    const auto result = IntegrateAffineTetrahedronAdaptive(
        children[child], rule, number_components, integrand, child_values[child],
        absolute_tolerance, relative_tolerance, depth + 1, max_subdivisions);
    for (std::size_t component = 0; component < number_components; component++)
    {
      value_accumulators[component].Add(static_cast<long double>(result.value[component]));
      error_accumulators[component].Add(
          static_cast<long double>(result.estimated_absolute_error[component]));
    }
    if (leaf_count > std::numeric_limits<std::size_t>::max() - result.leaf_count)
    {
      throw std::overflow_error(
          "Singular-element adaptive vector quadrature leaf count overflow!");
    }
    leaf_count += result.leaf_count;
    maximum_depth = std::max(maximum_depth, result.maximum_subdivision_depth);
    converged = converged && result.converged;
  }

  std::vector<double> value(number_components), error(number_components);
  for (std::size_t component = 0; component < number_components; component++)
  {
    value[component] = value_accumulators[component].Value();
    error[component] = error_accumulators[component].Value();
  }
  return {std::move(value), std::move(error), leaf_count, maximum_depth, converged};
}

template <typename Q>
ReferenceIntegral ComputeCanonicalReferenceIntegral(
    const CanonicalReferenceBasisPair &canonical, ReferenceQuadratureRule quadrature_rule,
    int quadrature_order, int subdivisions, double radial_power, Q &&quadrature)
{
  ReferenceIntegral result{canonical.row_basis,
                           canonical.column_basis,
                           canonical.input_to_canonical,
                           quadrature_rule,
                           quadrature_order,
                           subdivisions,
                           radial_power,
                           {},
                           {}};
  std::array<std::array<CompensatedAccumulator, 3>, 3> mass;
  std::array<std::array<CompensatedAccumulator, 3>, 3> curl_curl;
  const auto &grad_lambda = ReferenceBarycentricGradients();

  quadrature(
      [&](const BarycentricPoint &lambda, double weight)
      {
        const auto row = EvaluateReferenceBasis(lambda, grad_lambda, canonical.row_basis);
        const auto column =
            EvaluateReferenceBasis(lambda, grad_lambda, canonical.column_basis);
        ValidateFiniteBasisValue(row);
        ValidateFiniteBasisValue(column);
        for (int u = 0; u < 3; u++)
        {
          for (int v = 0; v < 3; v++)
          {
            mass[u][v].Add(static_cast<long double>(weight) *
                           static_cast<long double>(row.value[u]) *
                           static_cast<long double>(column.value[v]));
            curl_curl[u][v].Add(static_cast<long double>(weight) *
                                static_cast<long double>(row.curl[u]) *
                                static_cast<long double>(column.curl[v]));
          }
        }
      });

  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      result.mass[u][v] = mass[u][v].Value();
      result.curl_curl[u][v] = curl_curl[u][v].Value();
    }
  }
  return result;
}

FirstOrderReferenceIntegral
ConvertToFirstOrderReferenceIntegral(const ReferenceIntegral &integral)
{
  return {std::get<FirstOrderBasis>(integral.row_basis),
          std::get<FirstOrderBasis>(integral.column_basis),
          integral.input_to_canonical,
          integral.quadrature_rule,
          integral.quadrature_order,
          integral.subdivisions,
          integral.radial_power,
          integral.mass,
          integral.curl_curl};
}

FirstOrderAdaptiveReferenceIntegral
ConvertToFirstOrderReferenceIntegral(const AdaptiveReferenceIntegral &integral)
{
  return {ConvertToFirstOrderReferenceIntegral(integral.integral),
          integral.absolute_tolerance,
          integral.relative_tolerance,
          integral.mass_estimated_absolute_error,
          integral.curl_curl_estimated_absolute_error,
          integral.total_leaf_count,
          integral.maximum_subdivision_depth,
          integral.converged};
}

}  // namespace

BarycentricPoint ReferenceBarycentricPoint(const Vector3 &point)
{
  BarycentricPoint lambda{1.0 - point[0] - point[1] - point[2], point[0], point[1],
                          point[2]};
  ValidateBarycentricPoint(lambda);
  return lambda;
}

const BarycentricGradients &ReferenceBarycentricGradients()
{
  static constexpr BarycentricGradients GradLambda{
      Vector3{-1.0, -1.0, -1.0}, Vector3{1.0, 0.0, 0.0}, Vector3{0.0, 1.0, 0.0},
      Vector3{0.0, 0.0, 1.0}};
  return GradLambda;
}

TriangleBarycentricPoint ReferenceTriangleBarycentricPoint(const Vector2 &point)
{
  TriangleBarycentricPoint lambda{1.0 - point[0] - point[1], point[0], point[1]};
  ValidateTriangleBarycentricPoint(lambda);
  return lambda;
}

const TriangleBarycentricGradients &ReferenceTriangleBarycentricGradients()
{
  static constexpr TriangleBarycentricGradients GradLambda{
      Vector2{-1.0, -1.0}, Vector2{1.0, 0.0}, Vector2{0.0, 1.0}};
  return GradLambda;
}

ScalarPolynomialValue EvaluateSilvesterLagrange(int grid_denominator, int index,
                                                double coordinate)
{
  ValidateInterpolationPolynomial(grid_denominator, index, coordinate, false);
  ScalarPolynomialValue result{1.0, 0.0};
  for (int m = 0; m < index; m++)
  {
    const double denominator = static_cast<double>(m + 1);
    const double factor =
        (static_cast<double>(grid_denominator) * coordinate - m) / denominator;
    const double factor_derivative = static_cast<double>(grid_denominator) / denominator;
    result.derivative = result.derivative * factor + result.value * factor_derivative;
    result.value *= factor;
  }
  return result;
}

ScalarPolynomialValue EvaluateShiftedSilvesterLagrange(int grid_denominator, int index,
                                                       double coordinate)
{
  ValidateInterpolationPolynomial(grid_denominator, index, coordinate, true);
  ScalarPolynomialValue result{1.0, 0.0};
  for (int m = 1; m < index; m++)
  {
    const double denominator = static_cast<double>(m);
    const double factor =
        (static_cast<double>(grid_denominator) * coordinate - m) / denominator;
    const double factor_derivative = static_cast<double>(grid_denominator) / denominator;
    result.derivative = result.derivative * factor + result.value * factor_derivative;
    result.value *= factor;
  }
  return result;
}

double NodeRadialCoordinate(const BarycentricPoint &lambda, int i)
{
  ValidateBarycentricPoint(lambda);
  ValidateIndex(i);
  double rho = 0.0;
  for (int a = 0; a < 4; a++)
  {
    if (a != i)
    {
      rho += lambda[a];
    }
  }
  return rho;
}

double EdgeRadialCoordinate(const BarycentricPoint &lambda, int i, int j)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(i, j);
  double rho = 0.0;
  for (int a = 0; a < 4; a++)
  {
    if (a != i && a != j)
    {
      rho += lambda[a];
    }
  }
  return rho;
}

double TriangleNodeRadialCoordinate(const TriangleBarycentricPoint &lambda, int i)
{
  ValidateTriangleBarycentricPoint(lambda);
  ValidateTriangleIndex(i);
  double rho = 0.0;
  for (int a = 0; a < 3; a++)
  {
    if (a != i)
    {
      rho += lambda[a];
    }
  }
  return rho;
}

std::vector<WeightedSegmentQuadraturePoint>
BuildWeightedSegmentQuadrature(int order, double alpha, double beta)
{
  if (order < 1 || !std::isfinite(alpha) || !std::isfinite(beta) || !(alpha > -1.0) ||
      !(beta > -1.0))
  {
    throw std::invalid_argument(
        "Weighted segment quadrature requires positive order and finite endpoint "
        "exponents greater than -1!");
  }

  // Golub-Welsch for Jacobi weight (1-x)^beta (1+x)^alpha on [-1,1].
  // Mapping x = 2t-1 gives t^alpha (1-t)^beta on [0,1].
  const double jacobi_alpha = beta;
  const double jacobi_beta = alpha;
  const double exponent_sum = jacobi_alpha + jacobi_beta;
  Eigen::MatrixXd jacobi = Eigen::MatrixXd::Zero(order, order);
  for (int n = 0; n < order; n++)
  {
    const double two_n_sum = 2.0 * n + exponent_sum;
    jacobi(n, n) =
        (n == 0 && std::abs(exponent_sum) <= 32.0 * std::numeric_limits<double>::epsilon())
            ? (jacobi_beta - jacobi_alpha) / (exponent_sum + 2.0)
            : (jacobi_beta * jacobi_beta - jacobi_alpha * jacobi_alpha) /
                  (two_n_sum * (two_n_sum + 2.0));
    if (n == 0)
    {
      continue;
    }
    const double numerator =
        4.0 * n * (n + jacobi_alpha) * (n + jacobi_beta) * (n + exponent_sum);
    const double denominator =
        two_n_sum * two_n_sum * (two_n_sum - 1.0) * (two_n_sum + 1.0);
    const double recurrence = numerator / denominator;
    if (!std::isfinite(recurrence) || !(recurrence > 0.0))
    {
      throw std::runtime_error(
          "Weighted segment quadrature produced an invalid Jacobi recurrence!");
    }
    const double off_diagonal = std::sqrt(recurrence);
    jacobi(n - 1, n) = off_diagonal;
    jacobi(n, n - 1) = off_diagonal;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(jacobi);
  if (eigensolver.info() != Eigen::Success)
  {
    throw std::runtime_error(
        "Weighted segment quadrature eigendecomposition did not converge!");
  }
  const double total_weight = std::exp(std::lgamma(alpha + 1.0) + std::lgamma(beta + 1.0) -
                                       std::lgamma(alpha + beta + 2.0));
  if (!std::isfinite(total_weight) || !(total_weight > 0.0))
  {
    throw std::runtime_error(
        "Weighted segment quadrature produced an invalid total weight!");
  }

  std::vector<WeightedSegmentQuadraturePoint> rule;
  rule.reserve(order);
  for (int q = 0; q < order; q++)
  {
    const double coordinate = 0.5 * (eigensolver.eigenvalues()[q] + 1.0);
    const double leading_component = eigensolver.eigenvectors()(0, q);
    const double weight = total_weight * leading_component * leading_component;
    if (!std::isfinite(coordinate) || !(coordinate > 0.0) || !(coordinate < 1.0) ||
        !std::isfinite(weight) || !(weight > 0.0))
    {
      throw std::runtime_error(
          "Weighted segment quadrature produced an invalid point or weight!");
    }
    rule.push_back({coordinate, weight});
  }
  return rule;
}

TriangleVectorBasisValue
EvaluateTriangleStandardEdge(const TriangleBarycentricPoint &lambda,
                             const TriangleBarycentricGradients &grad_lambda, int i, int j)
{
  ValidateTriangleBarycentricPoint(lambda);
  ValidateTriangleBarycentricGradients(grad_lambda);
  ValidateDistinctTriangleIndices(i, j);

  TriangleVectorBasisValue result{{0.0, 0.0}, -2.0 * Cross(grad_lambda[i], grad_lambda[j])};
  Add(lambda[j], grad_lambda[i], result.value);
  Add(-lambda[i], grad_lambda[j], result.value);
  return result;
}

double EvaluateTriangleNodeGradientPotential(const TriangleBarycentricPoint &lambda, int i,
                                             int j, double nu)
{
  ValidateTriangleBarycentricPoint(lambda);
  ValidateDistinctTriangleIndices(i, j);
  ValidateExponent(nu);

  const double rho = TriangleNodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  return lambda[j] * OneMinusPositivePower(rho, nu - 1.0);
}

TriangleVectorBasisValue
EvaluateTriangleNodeGradient(const TriangleBarycentricPoint &lambda,
                             const TriangleBarycentricGradients &grad_lambda, int i, int j,
                             double nu)
{
  ValidateTriangleEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctTriangleIndices(i, j);

  const double rho = TriangleNodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  const double factor = OneMinusPositivePower(rho, nu - 1.0);
  const double derivative = lambda[j] * (nu - 1.0) * PositivePower(rho, nu - 2.0);

  TriangleVectorBasisValue result{{0.0, 0.0}, 0.0};
  Add(factor, grad_lambda[j], result.value);
  Add(derivative, grad_lambda[i], result.value);
  return result;
}

TriangleVectorBasisValue
EvaluateTriangleNodeRotational(const TriangleBarycentricPoint &lambda,
                               const TriangleBarycentricGradients &grad_lambda, int i,
                               int j, int k, double nu)
{
  ValidateTriangleEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctTriangleIndices(i, j, k);

  const double rho = TriangleNodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  const auto edge = EvaluateTriangleStandardEdge(lambda, grad_lambda, j, k);
  const double factor = PositivePowerMinusOne(rho, nu);
  const auto grad_factor = Scale(-nu * PositivePower(rho, nu - 1.0), grad_lambda[i]);

  return {Scale(factor, edge.value), factor * edge.curl + Cross(grad_factor, edge.value)};
}

BarycentricPoint ApplyBarycentricPermutation(const BarycentricPoint &lambda,
                                             const BarycentricPermutation &permutation)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(permutation[0], permutation[1], permutation[2], permutation[3]);
  BarycentricPoint result;
  for (int i = 0; i < 4; i++)
  {
    result[permutation[i]] = lambda[i];
  }
  return result;
}

BarycentricGradients ApplyBarycentricPermutation(const BarycentricGradients &grad_lambda,
                                                 const BarycentricPermutation &permutation)
{
  ValidateBarycentricGradients(grad_lambda);
  ValidateDistinctIndices(permutation[0], permutation[1], permutation[2], permutation[3]);
  BarycentricGradients result;
  for (int i = 0; i < 4; i++)
  {
    result[permutation[i]] = grad_lambda[i];
  }
  return result;
}

FirstOrderBasis ApplyBarycentricPermutation(const FirstOrderBasis &basis,
                                            const BarycentricPermutation &permutation)
{
  ValidateFirstOrderBasis(basis);
  ValidateDistinctIndices(permutation[0], permutation[1], permutation[2], permutation[3]);
  return ApplyBasisPermutation(basis, permutation);
}

HigherOrderBasis ApplyBarycentricPermutation(const HigherOrderBasis &basis,
                                             const BarycentricPermutation &permutation)
{
  ValidateHigherOrderBasis(basis);
  ValidateDistinctIndices(permutation[0], permutation[1], permutation[2], permutation[3]);
  HigherOrderBasis result = basis;
  for (int i = 0; i < 4; i++)
  {
    result.nodes[i] = permutation[basis.nodes[i]];
  }
  return result;
}

ReferenceBasis ApplyBarycentricPermutation(const ReferenceBasis &basis,
                                           const BarycentricPermutation &permutation)
{
  ValidateReferenceBasis(basis);
  return std::visit([&](const auto &entry) -> ReferenceBasis
                    { return ApplyBarycentricPermutation(entry, permutation); }, basis);
}

CanonicalFirstOrderBasisPair
CanonicalizeFirstOrderBasisPair(const FirstOrderBasis &row_basis,
                                const FirstOrderBasis &column_basis)
{
  ValidateFirstOrderBasis(row_basis);
  ValidateFirstOrderBasis(column_basis);

  BarycentricPermutation permutation{0, 1, 2, 3};
  CanonicalFirstOrderBasisPair result{row_basis, column_basis, permutation};
  auto result_key = BasisPairCanonicalKey(result.row_basis, result.column_basis,
                                          result.input_to_canonical);
  do
  {
    const auto row = ApplyBasisPermutation(row_basis, permutation);
    const auto column = ApplyBasisPermutation(column_basis, permutation);
    const auto key = BasisPairCanonicalKey(row, column, permutation);
    if (key < result_key)
    {
      result = {row, column, permutation};
      result_key = key;
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return result;
}

CanonicalReferenceBasisPair
CanonicalizeReferenceBasisPair(const ReferenceBasis &row_basis,
                               const ReferenceBasis &column_basis)
{
  ValidateReferenceBasis(row_basis);
  ValidateReferenceBasis(column_basis);

  BarycentricPermutation permutation{0, 1, 2, 3};
  CanonicalReferenceBasisPair result{row_basis, column_basis, permutation};
  auto result_key = ReferenceBasisPairCanonicalKey(result.row_basis, result.column_basis,
                                                   result.input_to_canonical);
  do
  {
    const auto row = ApplyBarycentricPermutation(row_basis, permutation);
    const auto column = ApplyBarycentricPermutation(column_basis, permutation);
    const auto key = ReferenceBasisPairCanonicalKey(row, column, permutation);
    if (key < result_key)
    {
      result = {row, column, permutation};
      result_key = key;
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return result;
}

VectorBasisValue EvaluateStandardEdge(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i, int j)
{
  ValidateBarycentricPoint(lambda);
  ValidateBarycentricGradients(grad_lambda);
  ValidateDistinctIndices(i, j);
  return StandardEdge(lambda, grad_lambda, i, j);
}

double EvaluateNodeGradientPotential(const BarycentricPoint &lambda, int i, int j,
                                     double nu)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(i, j);
  ValidateExponent(nu);

  const double rho = NodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  return lambda[j] * OneMinusPositivePower(rho, nu - 1.0);
}

VectorBasisValue EvaluateNodeGradient(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i, int j,
                                      double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j);

  const double rho = NodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  const double factor = OneMinusPositivePower(rho, nu - 1.0);
  const double derivative = lambda[j] * (nu - 1.0) * PositivePower(rho, nu - 2.0);

  VectorBasisValue result{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Add(factor, grad_lambda[j], result.value);
  Add(derivative, grad_lambda[i], result.value);
  return result;
}

double EvaluateHigherOrderNodeGradientPotential(
    const BarycentricPoint &lambda, int i, int j, int k, int l,
    const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(i, j, k, l);
  ValidateExponent(nu);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials = EvaluateHigherOrderNodeGradientPolynomials(
      lambda, nodes, interpolation_indices, order);
  return PolynomialProduct(polynomials) * EvaluateNodeGradientPotential(lambda, i, j, nu);
}

VectorBasisValue EvaluateHigherOrderNodeGradient(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k, l);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials = EvaluateHigherOrderNodeGradientPolynomials(
      lambda, nodes, interpolation_indices, order);
  const double polynomial = PolynomialProduct(polynomials);
  const auto grad_polynomial = PolynomialProductGradient(polynomials, grad_lambda, nodes);
  const double first_order_potential = EvaluateNodeGradientPotential(lambda, i, j, nu);
  const auto first_order = EvaluateNodeGradient(lambda, grad_lambda, i, j, nu);

  VectorBasisValue result{Scale(polynomial, first_order.value), {0.0, 0.0, 0.0}};
  Add(first_order_potential, grad_polynomial, result.value);
  return result;
}

VectorBasisValue EvaluateNodeRotational(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda, int i,
                                        int j, int k, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k);

  const double rho = NodeRadialCoordinate(lambda, i);
  ValidateRadialProxy(rho);
  const auto edge = StandardEdge(lambda, grad_lambda, j, k);
  const double factor = PositivePowerMinusOne(rho, nu);
  const auto grad_factor = Scale(-nu * PositivePower(rho, nu - 1.0), grad_lambda[i]);

  VectorBasisValue result{Scale(factor, edge.value), Scale(factor, edge.curl)};
  Add(1.0, Cross(grad_factor, edge.value), result.curl);
  return result;
}

VectorBasisValue EvaluateHigherOrderNodeRotational(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k, l);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials =
      EvaluateHigherOrderFacePolynomials(lambda, nodes, interpolation_indices, order);
  const double polynomial = PolynomialProduct(polynomials);
  const auto grad_polynomial = PolynomialProductGradient(polynomials, grad_lambda, nodes);
  const auto first_order = EvaluateNodeRotational(lambda, grad_lambda, i, j, k, nu);

  VectorBasisValue result{Scale(polynomial, first_order.value),
                          Scale(polynomial, first_order.curl)};
  Add(1.0, Cross(grad_polynomial, first_order.value), result.curl);
  return result;
}

double EvaluateEdgeGradientPotential(const BarycentricPoint &lambda, int i, int j, int k,
                                     double nu)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(i, j, k);
  ValidateExponent(nu);

  const double rho = EdgeRadialCoordinate(lambda, i, j);
  ValidateRadialProxy(rho);
  return lambda[i] * lambda[j] * lambda[k] * OneMinusPositivePower(rho, nu - 1.0);
}

VectorBasisValue EvaluateEdgeGradient(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda, int i, int j,
                                      int k, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k);

  const double rho = EdgeRadialCoordinate(lambda, i, j);
  ValidateRadialProxy(rho);
  const double polynomial = lambda[i] * lambda[j] * lambda[k];
  Vector3 grad_polynomial{0.0, 0.0, 0.0};
  Add(lambda[j] * lambda[k], grad_lambda[i], grad_polynomial);
  Add(lambda[i] * lambda[k], grad_lambda[j], grad_polynomial);
  Add(lambda[i] * lambda[j], grad_lambda[k], grad_polynomial);
  Vector3 grad_rho{0.0, 0.0, 0.0};
  Add(-1.0, grad_lambda[i], grad_rho);
  Add(-1.0, grad_lambda[j], grad_rho);

  VectorBasisValue result{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Add(OneMinusPositivePower(rho, nu - 1.0), grad_polynomial, result.value);
  Add(polynomial * (1.0 - nu) * PositivePower(rho, nu - 2.0), grad_rho, result.value);
  return result;
}

double EvaluateHigherOrderEdgeGradientPotential(
    const BarycentricPoint &lambda, int i, int j, int k, int l,
    const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateBarycentricPoint(lambda);
  ValidateDistinctIndices(i, j, k, l);
  ValidateExponent(nu);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials =
      EvaluateHigherOrderFacePolynomials(lambda, nodes, interpolation_indices, order);
  return PolynomialProduct(polynomials) *
         EvaluateEdgeGradientPotential(lambda, i, j, k, nu);
}

VectorBasisValue EvaluateHigherOrderEdgeGradient(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k, l);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials =
      EvaluateHigherOrderFacePolynomials(lambda, nodes, interpolation_indices, order);
  const double polynomial = PolynomialProduct(polynomials);
  const auto grad_polynomial = PolynomialProductGradient(polynomials, grad_lambda, nodes);
  const double first_order_potential = EvaluateEdgeGradientPotential(lambda, i, j, k, nu);
  const auto first_order = EvaluateEdgeGradient(lambda, grad_lambda, i, j, k, nu);

  VectorBasisValue result{Scale(polynomial, first_order.value), {0.0, 0.0, 0.0}};
  Add(first_order_potential, grad_polynomial, result.value);
  return result;
}

VectorBasisValue EvaluateEdgeRotational(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda, int i,
                                        int j, int k, int l, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k, l);

  const double rho = EdgeRadialCoordinate(lambda, i, j);
  ValidateRadialProxy(rho);
  const auto edge = StandardEdge(lambda, grad_lambda, k, l);
  const double radial_factor = PositivePowerMinusOne(rho, nu);
  const double polynomial = lambda[i] * lambda[j];
  const double factor = radial_factor * polynomial;

  Vector3 grad_polynomial{0.0, 0.0, 0.0};
  Add(lambda[j], grad_lambda[i], grad_polynomial);
  Add(lambda[i], grad_lambda[j], grad_polynomial);
  Vector3 grad_rho{0.0, 0.0, 0.0};
  Add(-1.0, grad_lambda[i], grad_rho);
  Add(-1.0, grad_lambda[j], grad_rho);
  Vector3 grad_factor = Scale(radial_factor, grad_polynomial);
  Add(polynomial * nu * PositivePower(rho, nu - 1.0), grad_rho, grad_factor);

  VectorBasisValue result{Scale(factor, edge.value), Scale(factor, edge.curl)};
  Add(1.0, Cross(grad_factor, edge.value), result.curl);
  return result;
}

VectorBasisValue EvaluateHigherOrderEdgeRotational(
    const BarycentricPoint &lambda, const BarycentricGradients &grad_lambda, int i, int j,
    int k, int l, const InterpolationIndices &interpolation_indices, int order, double nu)
{
  ValidateEvaluation(lambda, grad_lambda, nu);
  ValidateDistinctIndices(i, j, k, l);
  const std::array<int, 4> nodes{i, j, k, l};
  const auto polynomials =
      EvaluateHigherOrderVolumePolynomials(lambda, nodes, interpolation_indices, order);
  const double polynomial = PolynomialProduct(polynomials);
  const auto grad_polynomial = PolynomialProductGradient(polynomials, grad_lambda, nodes);
  const auto first_order = EvaluateEdgeRotational(lambda, grad_lambda, i, j, k, l, nu);

  VectorBasisValue result{Scale(polynomial, first_order.value),
                          Scale(polynomial, first_order.curl)};
  Add(1.0, Cross(grad_polynomial, first_order.value), result.curl);
  return result;
}

std::vector<HigherOrderBasis>
EnumerateHigherOrderNodeGradientBases(const std::array<int, 4> &canonical_nodes, int order,
                                      double nu)
{
  ValidateDistinctIndices(canonical_nodes[0], canonical_nodes[1], canonical_nodes[2],
                          canonical_nodes[3]);
  ValidateSingularOrder(order, 1);
  ValidateExponent(nu);
  const int grid_denominator = order + 1;
  std::vector<HigherOrderBasis> basis;
  for (int w0 = 1; w0 <= order; w0++)
  {
    for (int w1 = 0; w1 <= grid_denominator - w0; w1++)
    {
      for (int w2 = 0; w2 <= grid_denominator - w0 - w1; w2++)
      {
        const std::array<int, 4> weights{w0, w1, w2, grid_denominator - w0 - w1 - w2};
        int edge_position = -1;
        for (int position = 1; position < 4; position++)
        {
          if (weights[position] > 0)
          {
            edge_position = position;
            break;
          }
        }
        if (edge_position < 0)
        {
          continue;
        }

        std::array<int, 4> positions{0, edge_position, -1, -1};
        int opposite = 2;
        for (int position = 1; position < 4; position++)
        {
          if (position != edge_position)
          {
            positions[opposite++] = position;
          }
        }
        HigherOrderBasis entry{
            HigherOrderBasisFamily::NODE_GRADIENT,
            {canonical_nodes[positions[0]], canonical_nodes[positions[1]],
             canonical_nodes[positions[2]], canonical_nodes[positions[3]]},
            {weights[positions[0]], weights[positions[1]], weights[positions[2]],
             weights[positions[3]]},
            order,
            nu};
        ValidateHigherOrderBasis(entry);
        basis.push_back(entry);
      }
    }
  }
  return basis;
}

std::vector<HigherOrderBasis>
EnumerateHigherOrderNodeRotationalBases(const std::array<int, 4> &canonical_nodes,
                                        int order, double nu)
{
  ValidateDistinctIndices(canonical_nodes[0], canonical_nodes[1], canonical_nodes[2],
                          canonical_nodes[3]);
  ValidateSingularOrder(order, 2);
  ValidateExponent(nu);
  const int grid_denominator = order + 2;
  std::vector<HigherOrderBasis> basis;
  for (int w0 = 1; w0 <= order; w0++)
  {
    for (int w1 = 0; w1 <= grid_denominator - w0; w1++)
    {
      for (int w2 = 0; w2 <= grid_denominator - w0 - w1; w2++)
      {
        const std::array<int, 4> weights{w0, w1, w2, grid_denominator - w0 - w1 - w2};
        std::array<std::array<int, 2>, 3> face_positions;
        int face_count = 0;
        for (int first = 1; first < 4; first++)
        {
          for (int second = first + 1; second < 4; second++)
          {
            if (weights[first] > 0 && weights[second] > 0)
            {
              face_positions[face_count++] = {first, second};
            }
          }
        }
        const int retained_faces = std::min(face_count, 2);
        for (int face = 0; face < retained_faces; face++)
        {
          const int first = face_positions[face][0];
          const int second = face_positions[face][1];
          int opposite = -1;
          for (int position = 1; position < 4; position++)
          {
            if (position != first && position != second)
            {
              opposite = position;
              break;
            }
          }
          HigherOrderBasis entry{
              HigherOrderBasisFamily::NODE_ROTATIONAL,
              {canonical_nodes[0], canonical_nodes[first], canonical_nodes[second],
               canonical_nodes[opposite]},
              {weights[0], weights[first], weights[second], weights[opposite]},
              order,
              nu};
          ValidateHigherOrderBasis(entry);
          basis.push_back(entry);
        }
      }
    }
  }
  return basis;
}

std::vector<HigherOrderBasis>
EnumerateHigherOrderEdgeGradientBases(const std::array<int, 4> &canonical_nodes, int order,
                                      double nu)
{
  ValidateDistinctIndices(canonical_nodes[0], canonical_nodes[1], canonical_nodes[2],
                          canonical_nodes[3]);
  ValidateSingularOrder(order, 2);
  ValidateExponent(nu);
  const int grid_denominator = order + 2;
  std::vector<HigherOrderBasis> basis;
  for (int w0 = 1; w0 <= order; w0++)
  {
    for (int w1 = 1; w1 <= order; w1++)
    {
      const int remaining = grid_denominator - w0 - w1;
      if (remaining <= 0)
      {
        continue;
      }
      for (int w2 = 0; w2 <= remaining; w2++)
      {
        const std::array<int, 4> weights{w0, w1, w2, remaining - w2};
        const int face_position = weights[2] > 0 ? 2 : 3;
        const int opposite_position = face_position == 2 ? 3 : 2;
        HigherOrderBasis entry{
            HigherOrderBasisFamily::EDGE_GRADIENT,
            {canonical_nodes[0], canonical_nodes[1], canonical_nodes[face_position],
             canonical_nodes[opposite_position]},
            {weights[0], weights[1], weights[face_position], weights[opposite_position]},
            order,
            nu};
        ValidateHigherOrderBasis(entry);
        basis.push_back(entry);
      }
    }
  }
  return basis;
}

std::vector<HigherOrderBasis>
EnumerateHigherOrderEdgeRotationalBases(const std::array<int, 4> &canonical_nodes,
                                        int order, double nu)
{
  ValidateDistinctIndices(canonical_nodes[0], canonical_nodes[1], canonical_nodes[2],
                          canonical_nodes[3]);
  ValidateSingularOrder(order, 3);
  ValidateExponent(nu);
  const int grid_denominator = order + 3;
  std::vector<HigherOrderBasis> basis;
  for (int w0 = 1; w0 <= order; w0++)
  {
    for (int w1 = 1; w1 <= order; w1++)
    {
      for (int w2 = 1; w2 <= order; w2++)
      {
        const int w3 = grid_denominator - w0 - w1 - w2;
        if (w3 < 1 || w3 > order)
        {
          continue;
        }
        HigherOrderBasis entry{HigherOrderBasisFamily::EDGE_ROTATIONAL,
                               canonical_nodes,
                               {w0, w1, w2, w3},
                               order,
                               nu};
        ValidateHigherOrderBasis(entry);
        basis.push_back(entry);
      }
    }
  }
  return basis;
}

VectorBasisValue EvaluateHigherOrderBasis(const BarycentricPoint &lambda,
                                          const BarycentricGradients &grad_lambda,
                                          const HigherOrderBasis &basis)
{
  ValidateHigherOrderBasis(basis);
  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
      return EvaluateHigherOrderNodeGradient(
          lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
          basis.nodes[3], basis.interpolation_indices, basis.order, basis.nu);
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
      return EvaluateHigherOrderNodeRotational(
          lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
          basis.nodes[3], basis.interpolation_indices, basis.order, basis.nu);
    case HigherOrderBasisFamily::EDGE_GRADIENT:
      return EvaluateHigherOrderEdgeGradient(
          lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
          basis.nodes[3], basis.interpolation_indices, basis.order, basis.nu);
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      return EvaluateHigherOrderEdgeRotational(
          lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
          basis.nodes[3], basis.interpolation_indices, basis.order, basis.nu);
  }
  throw std::invalid_argument("Unknown higher-order singular-element basis family!");
}

Vector3 EvaluateHigherOrderBasisValue(const BarycentricPoint &lambda,
                                      const BarycentricGradients &grad_lambda,
                                      const HigherOrderBasis &basis)
{
  ValidateHigherOrderBasis(basis);
  if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT ||
      basis.family == HigherOrderBasisFamily::EDGE_GRADIENT)
  {
    return EvaluateHigherOrderBasis(lambda, grad_lambda, basis).value;
  }

  const auto &nodes = basis.nodes;
  if (basis.family == HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    ValidateEvaluation(lambda, grad_lambda, basis.nu);
    const auto polynomials = EvaluateHigherOrderFacePolynomials(
        lambda, nodes, basis.interpolation_indices, basis.order);
    const double rho = NodeRadialCoordinate(lambda, nodes[0]);
    constexpr double endpoint_tolerance = 128.0 * std::numeric_limits<double>::epsilon();
    if (rho <= endpoint_tolerance)
    {
      // The rotational field vanishes continuously at its singular node even
      // though its curl is singular there.
      return {0.0, 0.0, 0.0};
    }
    const auto edge = StandardEdge(lambda, grad_lambda, nodes[1], nodes[2]);
    const double factor =
        PolynomialProduct(polynomials) * PositivePowerMinusOne(rho, basis.nu);
    return Scale(factor, edge.value);
  }
  if (basis.family == HigherOrderBasisFamily::EDGE_ROTATIONAL)
  {
    ValidateEvaluation(lambda, grad_lambda, basis.nu);
    const auto polynomials = EvaluateHigherOrderVolumePolynomials(
        lambda, nodes, basis.interpolation_indices, basis.order);
    const double rho = EdgeRadialCoordinate(lambda, nodes[0], nodes[1]);
    constexpr double endpoint_tolerance = 128.0 * std::numeric_limits<double>::epsilon();
    if (rho <= endpoint_tolerance)
    {
      // The damping and transverse Nedelec factor give a continuous zero
      // value on the singular edge. Its curl remains singular.
      return {0.0, 0.0, 0.0};
    }
    const auto edge = StandardEdge(lambda, grad_lambda, nodes[2], nodes[3]);
    const double factor = PolynomialProduct(polynomials) * lambda[nodes[0]] *
                          lambda[nodes[1]] * PositivePowerMinusOne(rho, basis.nu);
    return Scale(factor, edge.value);
  }
  throw std::invalid_argument("Unknown higher-order singular-element basis family!");
}

double EvaluateHigherOrderGradientPotential(const BarycentricPoint &lambda,
                                            const HigherOrderBasis &basis)
{
  ValidateHigherOrderBasis(basis);
  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
      return EvaluateHigherOrderNodeGradientPotential(
          lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
          basis.interpolation_indices, basis.order, basis.nu);
    case HigherOrderBasisFamily::EDGE_GRADIENT:
      return EvaluateHigherOrderEdgeGradientPotential(
          lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
          basis.interpolation_indices, basis.order, basis.nu);
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      throw std::invalid_argument(
          "A rotational singular basis has no scalar H1 potential!");
  }
  throw std::invalid_argument("Unknown higher-order singular-element basis family!");
}

VectorBasisValue EvaluateReferenceBasis(const BarycentricPoint &lambda,
                                        const BarycentricGradients &grad_lambda,
                                        const ReferenceBasis &basis)
{
  ValidateReferenceBasis(basis);
  return std::visit(
      [&](const auto &entry)
      {
        using Basis = std::decay_t<decltype(entry)>;
        if constexpr (std::is_same_v<Basis, FirstOrderBasis>)
        {
          return EvaluateFirstOrderBasis(lambda, grad_lambda, entry);
        }
        else
        {
          return EvaluateHigherOrderBasis(lambda, grad_lambda, entry);
        }
      },
      basis);
}

FirstOrderBasis MakeStandardH1Gradient(int i)
{
  ValidateIndex(i);
  return {FirstOrderBasisFamily::STANDARD_H1_GRADIENT, {i, -1, -1, -1}, 0.0};
}

FirstOrderBasis MakeStandardNedelec(int i, int j)
{
  ValidateDistinctIndices(i, j);
  return {FirstOrderBasisFamily::STANDARD_NEDELEC, {i, j, -1, -1}, 0.0};
}

FirstOrderBasis MakeNodeGradient(int i, int j, double nu)
{
  ValidateDistinctIndices(i, j);
  ValidateExponent(nu);
  return {FirstOrderBasisFamily::NODE_GRADIENT, {i, j, -1, -1}, nu};
}

FirstOrderBasis MakeNodeRotational(int i, int j, int k, double nu)
{
  ValidateDistinctIndices(i, j, k);
  ValidateExponent(nu);
  return {FirstOrderBasisFamily::NODE_ROTATIONAL, {i, j, k, -1}, nu};
}

FirstOrderBasis MakeEdgeGradient(int i, int j, int k, double nu)
{
  ValidateDistinctIndices(i, j, k);
  ValidateExponent(nu);
  return {FirstOrderBasisFamily::EDGE_GRADIENT, {i, j, k, -1}, nu};
}

FirstOrderBasis MakeEdgeRotational(int i, int j, int k, int l, double nu)
{
  ValidateDistinctIndices(i, j, k, l);
  ValidateExponent(nu);
  return {FirstOrderBasisFamily::EDGE_ROTATIONAL, {i, j, k, l}, nu};
}

VectorBasisValue EvaluateFirstOrderBasis(const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda,
                                         const FirstOrderBasis &basis)
{
  switch (basis.family)
  {
    case FirstOrderBasisFamily::STANDARD_H1_GRADIENT:
      ValidateBarycentricPoint(lambda);
      ValidateBarycentricGradients(grad_lambda);
      ValidateIndex(basis.indices[0]);
      return {grad_lambda[basis.indices[0]], {0.0, 0.0, 0.0}};
    case FirstOrderBasisFamily::STANDARD_NEDELEC:
      return EvaluateStandardEdge(lambda, grad_lambda, basis.indices[0], basis.indices[1]);
    case FirstOrderBasisFamily::NODE_GRADIENT:
      return EvaluateNodeGradient(lambda, grad_lambda, basis.indices[0], basis.indices[1],
                                  basis.nu);
    case FirstOrderBasisFamily::NODE_ROTATIONAL:
      return EvaluateNodeRotational(lambda, grad_lambda, basis.indices[0], basis.indices[1],
                                    basis.indices[2], basis.nu);
    case FirstOrderBasisFamily::EDGE_GRADIENT:
      return EvaluateEdgeGradient(lambda, grad_lambda, basis.indices[0], basis.indices[1],
                                  basis.indices[2], basis.nu);
    case FirstOrderBasisFamily::EDGE_ROTATIONAL:
      return EvaluateEdgeRotational(lambda, grad_lambda, basis.indices[0], basis.indices[1],
                                    basis.indices[2], basis.indices[3], basis.nu);
  }
  throw std::invalid_argument("Unknown first-order singular-element basis family!");
}

void ForEachReferenceTetrahedronQuadraturePoint(int order, int subdivisions,
                                                const QuadraturePointVisitor &visitor)
{
  ValidateQuadratureParameters(order, subdivisions);
  ValidateQuadratureVisitor(visitor);

  const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, order);
  VisitQuadraturePoints(ReferenceTetrahedron, subdivisions, rule, visitor);
}

void ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
    int order, int singular_node, double radial_power,
    const QuadraturePointVisitor &visitor)
{
  ValidateIndex(singular_node);
  ValidateQuadratureVisitor(visitor);

  std::array<int, 3> other_nodes;
  int other = 0;
  for (int i = 0; i < 4; i++)
  {
    if (i != singular_node)
    {
      other_nodes[other++] = i;
    }
  }

  const auto &rule = GetDuffyQuadratureRule(order, radial_power);
  ValidateNodeDuffyQuadratureRule(rule, radial_power);
  for (int q_r = 0; q_r < rule.GetNPoints(); q_r++)
  {
    const auto &point_r = rule.IntPoint(q_r);
    const double parameter = point_r.x;
    const double r = std::pow(parameter, radial_power);
    const double derivative = radial_power * r / parameter;
    for (int q_u = 0; q_u < rule.GetNPoints(); q_u++)
    {
      const auto &point_u = rule.IntPoint(q_u);
      const double u = point_u.x;
      for (int q_v = 0; q_v < rule.GetNPoints(); q_v++)
      {
        const auto &point_v = rule.IntPoint(q_v);
        const double v = point_v.x;
        BarycentricPoint lambda;
        lambda[singular_node] = 1.0 - r;
        lambda[other_nodes[0]] = r * (1.0 - u);
        lambda[other_nodes[1]] = r * u * (1.0 - v);
        lambda[other_nodes[2]] = r * u * v;
        const double weight =
            point_r.weight * point_u.weight * point_v.weight * derivative * r * r * u;
        ValidateBarycentricPoint(lambda);
        visitor(lambda, weight);
      }
    }
  }
}

void ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
    int order, int singular_node_i, int singular_node_j, double radial_power,
    const QuadraturePointVisitor &visitor)
{
  ValidateDistinctIndices(singular_node_i, singular_node_j);
  ValidateQuadratureVisitor(visitor);

  std::array<int, 2> other_nodes;
  int other = 0;
  for (int i = 0; i < 4; i++)
  {
    if (i != singular_node_i && i != singular_node_j)
    {
      other_nodes[other++] = i;
    }
  }

  const auto &rule = GetDuffyQuadratureRule(order, radial_power);
  for (int q_rho = 0; q_rho < rule.GetNPoints(); q_rho++)
  {
    const auto &point_rho = rule.IntPoint(q_rho);
    const double parameter = point_rho.x;
    const double rho = std::pow(parameter, radial_power);
    const double derivative = radial_power * rho / parameter;
    for (int q_t = 0; q_t < rule.GetNPoints(); q_t++)
    {
      const auto &point_t = rule.IntPoint(q_t);
      const double t = point_t.x;
      for (int q_u = 0; q_u < rule.GetNPoints(); q_u++)
      {
        const auto &point_u = rule.IntPoint(q_u);
        const double u = point_u.x;
        BarycentricPoint lambda;
        lambda[singular_node_i] = (1.0 - rho) * (1.0 - t);
        lambda[singular_node_j] = (1.0 - rho) * t;
        lambda[other_nodes[0]] = rho * (1.0 - u);
        lambda[other_nodes[1]] = rho * u;
        const double weight = point_rho.weight * point_t.weight * point_u.weight *
                              derivative * rho * (1.0 - rho);
        ValidateBarycentricPoint(lambda);
        visitor(lambda, weight);
      }
    }
  }
}

void ForEachReferenceTriangleNodeDuffyQuadraturePoint(
    int order, int singular_node, double radial_power,
    const TriangleQuadraturePointVisitor &visitor)
{
  ValidateTriangleIndex(singular_node);
  if (!visitor)
  {
    throw std::invalid_argument(
        "Singular-element triangle quadrature visitor must be callable!");
  }

  std::array<int, 2> other_nodes;
  int other = 0;
  for (int i = 0; i < 3; i++)
  {
    if (i != singular_node)
    {
      other_nodes[other++] = i;
    }
  }

  const auto &rule = GetDuffyQuadratureRule(order, radial_power);
  for (int q_r = 0; q_r < rule.GetNPoints(); q_r++)
  {
    const auto &point_r = rule.IntPoint(q_r);
    const double parameter = point_r.x;
    const double r = std::pow(parameter, radial_power);
    const double derivative = radial_power * r / parameter;
    for (int q_t = 0; q_t < rule.GetNPoints(); q_t++)
    {
      const auto &point_t = rule.IntPoint(q_t);
      const double t = point_t.x;
      TriangleBarycentricPoint lambda;
      lambda[singular_node] = 1.0 - r;
      lambda[other_nodes[0]] = r * (1.0 - t);
      lambda[other_nodes[1]] = r * t;
      const double weight = point_r.weight * point_t.weight * derivative * r;
      ValidateTriangleBarycentricPoint(lambda);
      visitor(lambda, weight);
    }
  }
}

double IntegrateReferenceTriangleNodeDuffy(int order, int singular_node,
                                           double radial_power,
                                           const TriangleReferenceIntegrand &integrand)
{
  if (!integrand)
  {
    throw std::invalid_argument(
        "Singular-element triangle quadrature integrand must be callable!");
  }
  CompensatedAccumulator result;
  ForEachReferenceTriangleNodeDuffyQuadraturePoint(
      order, singular_node, radial_power,
      [&result, &integrand](const TriangleBarycentricPoint &lambda, double weight)
      {
        const double value = integrand(lambda);
        if (!std::isfinite(value))
        {
          throw std::domain_error(
              "Singular-element triangle integrand returned a non-finite value!");
        }
        result.Add(static_cast<long double>(weight) * value);
      });
  return result.Value();
}

std::size_t ReferenceTriangleNodeDuffyQuadraturePointCount(int order)
{
  if (order < 1 || order > MaximumDuffyQuadratureOrder)
  {
    throw std::invalid_argument(
        "Singular-element Duffy quadrature order must be in [1, 127]!");
  }
  const std::size_t points = static_cast<std::size_t>(
      mfem::IntRules.Get(mfem::Geometry::SEGMENT, order).GetNPoints());
  if (points == 0 || points > std::numeric_limits<std::size_t>::max() / points)
  {
    throw std::overflow_error("Singular-element triangle quadrature point count overflow!");
  }
  return points * points;
}

std::size_t ReferenceTetrahedronQuadraturePointCount(int order, int subdivisions)
{
  ValidateQuadratureParameters(order, subdivisions);
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, order);
  std::size_t count = static_cast<std::size_t>(rule.GetNPoints());
  for (int level = 0; level < subdivisions; level++)
  {
    if (count > std::numeric_limits<std::size_t>::max() / 8)
    {
      throw std::overflow_error("Singular-element quadrature point count overflow!");
    }
    count *= 8;
  }
  return count;
}

double IntegrateReferenceTetrahedron(int order, int subdivisions,
                                     const ReferenceIntegrand &integrand)
{
  return IntegrateReferenceQuadrature(
      [&](const QuadraturePointVisitor &visitor)
      { ForEachReferenceTetrahedronQuadraturePoint(order, subdivisions, visitor); },
      integrand);
}

std::vector<double> IntegrateReferenceTetrahedron(int order, int subdivisions,
                                                  std::size_t number_components,
                                                  const ReferenceVectorIntegrand &integrand)
{
  ValidateQuadratureParameters(order, subdivisions);
  if (!integrand || number_components == 0)
  {
    throw std::invalid_argument(
        "Singular-element reference vector integrand must be callable and nonempty!");
  }
  std::vector<CompensatedAccumulator> accumulators(number_components);
  std::vector<double> values(number_components);
  ForEachReferenceTetrahedronQuadraturePoint(
      order, subdivisions,
      [&](const BarycentricPoint &lambda, double weight)
      {
        std::fill(values.begin(), values.end(), std::numeric_limits<double>::quiet_NaN());
        integrand(lambda, values);
        if (values.size() != number_components ||
            std::any_of(values.begin(), values.end(),
                        [](double value) { return !std::isfinite(value); }))
        {
          throw std::domain_error(
              "Singular-element reference vector integrand returned invalid values!");
        }
        for (std::size_t component = 0; component < number_components; component++)
        {
          accumulators[component].Add(static_cast<long double>(weight) * values[component]);
        }
      });
  std::vector<double> result(number_components);
  for (std::size_t component = 0; component < number_components; component++)
  {
    result[component] = accumulators[component].Value();
  }
  return result;
}

AdaptiveQuadratureResult
IntegrateReferenceTetrahedronAdaptive(int order, double absolute_tolerance,
                                      double relative_tolerance, int max_subdivisions,
                                      const ReferenceIntegrand &integrand)
{
  ValidateQuadratureParameters(order, 0);
  ValidateAdaptiveQuadratureParameters(absolute_tolerance, relative_tolerance,
                                       max_subdivisions);
  if (!integrand)
  {
    throw std::invalid_argument("Singular-element reference integrand must be callable!");
  }
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, order);
  const double coarse_value =
      IntegrateAffineTetrahedron(ReferenceTetrahedron, rule, integrand);
  double local_absolute_tolerance = absolute_tolerance;
  double local_relative_tolerance = relative_tolerance;
  AdaptiveQuadratureResult result{};
  for (int attempt = 0; attempt <= max_subdivisions; attempt++)
  {
    result = IntegrateAffineTetrahedronAdaptive(
        ReferenceTetrahedron, rule, integrand, coarse_value, local_absolute_tolerance,
        local_relative_tolerance, 0, max_subdivisions);
    const double requested_tolerance =
        absolute_tolerance + relative_tolerance * std::abs(result.value);
    result.converged = result.estimated_absolute_error <= requested_tolerance;
    if (result.converged || result.maximum_subdivision_depth == max_subdivisions)
    {
      return result;
    }
    const double tightening =
        std::min(0.9, 0.8 * requested_tolerance / result.estimated_absolute_error);
    local_absolute_tolerance *= tightening;
    local_relative_tolerance *= tightening;
  }
  return result;
}

AdaptiveVectorQuadratureResult IntegrateReferenceTetrahedronAdaptive(
    int order, double absolute_tolerance, double relative_tolerance, int max_subdivisions,
    std::size_t number_components, const ReferenceVectorIntegrand &integrand)
{
  ValidateQuadratureParameters(order, 0);
  ValidateAdaptiveQuadratureParameters(absolute_tolerance, relative_tolerance,
                                       max_subdivisions);
  if (!integrand || number_components == 0)
  {
    throw std::invalid_argument(
        "Singular-element reference vector integrand must be callable and nonempty!");
  }

  const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, order);
  const auto coarse_value =
      IntegrateAffineTetrahedron(ReferenceTetrahedron, rule, number_components, integrand);
  double local_absolute_tolerance = absolute_tolerance;
  double local_relative_tolerance = relative_tolerance;
  AdaptiveVectorQuadratureResult result;
  for (int attempt = 0; attempt <= max_subdivisions; attempt++)
  {
    result = IntegrateAffineTetrahedronAdaptive(
        ReferenceTetrahedron, rule, number_components, integrand, coarse_value,
        local_absolute_tolerance, local_relative_tolerance, 0, max_subdivisions);
    result.converged = true;
    double tightening = 0.9;
    for (std::size_t component = 0; component < number_components; component++)
    {
      const double requested_tolerance =
          absolute_tolerance + relative_tolerance * std::abs(result.value[component]);
      const bool component_converged =
          result.estimated_absolute_error[component] <= requested_tolerance;
      result.converged = result.converged && component_converged;
      if (!component_converged)
      {
        tightening = std::min(tightening, 0.8 * requested_tolerance /
                                              result.estimated_absolute_error[component]);
      }
    }
    if (result.converged || result.maximum_subdivision_depth == max_subdivisions)
    {
      return result;
    }
    local_absolute_tolerance *= tightening;
    local_relative_tolerance *= tightening;
  }
  return result;
}

double IntegrateReferenceTetrahedronNodeDuffy(int order, int singular_node,
                                              double radial_power,
                                              const ReferenceIntegrand &integrand)
{
  return IntegrateReferenceQuadrature(
      [&](const QuadraturePointVisitor &visitor)
      {
        ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(order, singular_node,
                                                            radial_power, visitor);
      },
      integrand);
}

double IntegrateReferenceTetrahedronEdgeDuffy(int order, int singular_node_i,
                                              int singular_node_j, double radial_power,
                                              const ReferenceIntegrand &integrand)
{
  return IntegrateReferenceQuadrature(
      [&](const QuadraturePointVisitor &visitor)
      {
        ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
            order, singular_node_i, singular_node_j, radial_power, visitor);
      },
      integrand);
}

ReferenceIntegral ComputeReferenceIntegral(const ReferenceBasis &row_basis,
                                           const ReferenceBasis &column_basis, int order,
                                           int subdivisions)
{
  const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
  return ComputeCanonicalReferenceIntegral(
      canonical, ReferenceQuadratureRule::RECURSIVE_TETRAHEDRON, order, subdivisions, 0.0,
      [&](const QuadraturePointVisitor &visitor)
      { ForEachReferenceTetrahedronQuadraturePoint(order, subdivisions, visitor); });
}

AdaptiveReferenceIntegral ComputeAdaptiveReferenceIntegral(
    const ReferenceBasis &row_basis, const ReferenceBasis &column_basis, int order,
    double absolute_tolerance, double relative_tolerance, int max_subdivisions)
{
  const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
  AdaptiveReferenceIntegral result{{canonical.row_basis,
                                    canonical.column_basis,
                                    canonical.input_to_canonical,
                                    ReferenceQuadratureRule::ADAPTIVE_TETRAHEDRON,
                                    order,
                                    max_subdivisions,
                                    0.0,
                                    {},
                                    {}},
                                   absolute_tolerance,
                                   relative_tolerance,
                                   {},
                                   {},
                                   0,
                                   0,
                                   true};
  const auto &grad_lambda = ReferenceBarycentricGradients();
  const auto record_result = [&](const AdaptiveQuadratureResult &entry,
                                 ReferenceCoefficientTensor &tensor,
                                 ReferenceCoefficientTensor &error, int u, int v)
  {
    tensor[u][v] = entry.value;
    error[u][v] = entry.estimated_absolute_error;
    if (result.total_leaf_count >
        std::numeric_limits<std::size_t>::max() - entry.leaf_count)
    {
      throw std::overflow_error(
          "Singular-element adaptive tensor quadrature leaf count overflow!");
    }
    result.total_leaf_count += entry.leaf_count;
    result.maximum_subdivision_depth =
        std::max(result.maximum_subdivision_depth, entry.maximum_subdivision_depth);
    result.converged = result.converged && entry.converged;
  };

  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      const auto mass = IntegrateReferenceTetrahedronAdaptive(
          order, absolute_tolerance, relative_tolerance, max_subdivisions,
          [&](const BarycentricPoint &lambda)
          {
            const auto row =
                EvaluateReferenceBasis(lambda, grad_lambda, canonical.row_basis);
            const auto column =
                EvaluateReferenceBasis(lambda, grad_lambda, canonical.column_basis);
            ValidateFiniteBasisValue(row);
            ValidateFiniteBasisValue(column);
            return row.value[u] * column.value[v];
          });
      record_result(mass, result.integral.mass, result.mass_estimated_absolute_error, u, v);

      const auto curl_curl = IntegrateReferenceTetrahedronAdaptive(
          order, absolute_tolerance, relative_tolerance, max_subdivisions,
          [&](const BarycentricPoint &lambda)
          {
            const auto row =
                EvaluateReferenceBasis(lambda, grad_lambda, canonical.row_basis);
            const auto column =
                EvaluateReferenceBasis(lambda, grad_lambda, canonical.column_basis);
            ValidateFiniteBasisValue(row);
            ValidateFiniteBasisValue(column);
            return row.curl[u] * column.curl[v];
          });
      record_result(curl_curl, result.integral.curl_curl,
                    result.curl_curl_estimated_absolute_error, u, v);
    }
  }
  return result;
}

ReferenceIntegral ComputeNodeDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                                    const ReferenceBasis &column_basis,
                                                    int order, double radial_power)
{
  GetCommonNodeFeature(row_basis, column_basis);
  const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
  const int singular_node =
      GetCommonNodeFeature(canonical.row_basis, canonical.column_basis);
  return ComputeCanonicalReferenceIntegral(
      canonical, ReferenceQuadratureRule::NODE_DUFFY, order, 0, radial_power,
      [&](const QuadraturePointVisitor &visitor)
      {
        ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(order, singular_node,
                                                            radial_power, visitor);
      });
}

ReferenceIntegral ComputeEdgeDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                                    const ReferenceBasis &column_basis,
                                                    int order, double radial_power)
{
  GetCommonEdgeFeature(row_basis, column_basis);
  const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
  const auto singular_edge =
      GetCommonEdgeFeature(canonical.row_basis, canonical.column_basis);
  return ComputeCanonicalReferenceIntegral(
      canonical, ReferenceQuadratureRule::EDGE_DUFFY, order, 0, radial_power,
      [&](const QuadraturePointVisitor &visitor)
      {
        ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
            order, singular_edge[0], singular_edge[1], radial_power, visitor);
      });
}

ReferenceIntegral
ComputePartitionedDuffyReferenceIntegral(const ReferenceBasis &row_basis,
                                         const ReferenceBasis &column_basis, int order,
                                         double radial_power)
{
  const auto row_feature = GetSingularFeature(row_basis);
  const auto column_feature = GetSingularFeature(column_basis);
  if (row_feature.kind == SingularFeatureKind::NONE ||
      column_feature.kind == SingularFeatureKind::NONE ||
      SameFeature(row_feature, column_feature))
  {
    throw std::invalid_argument(
        "Partitioned Duffy quadrature requires two distinct singular features!");
  }

  const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
  const auto canonical_row_feature = GetSingularFeature(canonical.row_basis);
  const auto canonical_column_feature = GetSingularFeature(canonical.column_basis);
  return ComputeCanonicalReferenceIntegral(
      canonical, ReferenceQuadratureRule::PARTITIONED_DUFFY, order, 0, radial_power,
      [&](const QuadraturePointVisitor &visitor)
      {
        const auto emit_chart = [&](const SingularFeature &aligned_feature,
                                    const SingularFeature &other_feature)
        {
          ForEachFeatureDuffyQuadraturePoint(
              aligned_feature, order, radial_power,
              [&](const BarycentricPoint &lambda, double weight)
              {
                visitor(lambda, weight * DuffyPartitionWeight(lambda, aligned_feature,
                                                              other_feature));
              });
        };
        emit_chart(canonical_row_feature, canonical_column_feature);
        emit_chart(canonical_column_feature, canonical_row_feature);
      });
}

FirstOrderReferenceIntegral
ComputeFirstOrderReferenceIntegral(const FirstOrderBasis &row_basis,
                                   const FirstOrderBasis &column_basis, int order,
                                   int subdivisions)
{
  return ConvertToFirstOrderReferenceIntegral(
      ComputeReferenceIntegral(row_basis, column_basis, order, subdivisions));
}

FirstOrderAdaptiveReferenceIntegral ComputeFirstOrderAdaptiveReferenceIntegral(
    const FirstOrderBasis &row_basis, const FirstOrderBasis &column_basis, int order,
    double absolute_tolerance, double relative_tolerance, int max_subdivisions)
{
  return ConvertToFirstOrderReferenceIntegral(
      ComputeAdaptiveReferenceIntegral(row_basis, column_basis, order, absolute_tolerance,
                                       relative_tolerance, max_subdivisions));
}

FirstOrderReferenceIntegral
ComputeFirstOrderNodeDuffyReferenceIntegral(const FirstOrderBasis &row_basis,
                                            const FirstOrderBasis &column_basis, int order,
                                            double radial_power)
{
  return ConvertToFirstOrderReferenceIntegral(
      ComputeNodeDuffyReferenceIntegral(row_basis, column_basis, order, radial_power));
}

FirstOrderReferenceIntegral
ComputeFirstOrderEdgeDuffyReferenceIntegral(const FirstOrderBasis &row_basis,
                                            const FirstOrderBasis &column_basis, int order,
                                            double radial_power)
{
  return ConvertToFirstOrderReferenceIntegral(
      ComputeEdgeDuffyReferenceIntegral(row_basis, column_basis, order, radial_power));
}

double ContractMass(const ReferenceIntegral &integral,
                    const BarycentricGradients &grad_lambda, double jacobian_determinant)
{
  ValidateAffineGeometry(grad_lambda, jacobian_determinant);
  const auto canonical_grad_lambda =
      ApplyBarycentricPermutation(grad_lambda, integral.input_to_canonical);
  CompensatedAccumulator accumulator;
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      accumulator.Add(static_cast<long double>(jacobian_determinant) *
                      static_cast<long double>(integral.mass[u][v]) *
                      static_cast<long double>(
                          Dot(canonical_grad_lambda[u + 1], canonical_grad_lambda[v + 1])));
    }
  }
  return accumulator.Value();
}

double ContractCurlCurl(const ReferenceIntegral &integral,
                        const BarycentricGradients &grad_lambda,
                        double jacobian_determinant)
{
  ValidateAffineGeometry(grad_lambda, jacobian_determinant);
  const auto canonical_grad_lambda =
      ApplyBarycentricPermutation(grad_lambda, integral.input_to_canonical);
  const std::array<Vector3, 3> curl_geometry{
      Cross(canonical_grad_lambda[2], canonical_grad_lambda[3]),
      Cross(canonical_grad_lambda[3], canonical_grad_lambda[1]),
      Cross(canonical_grad_lambda[1], canonical_grad_lambda[2])};
  CompensatedAccumulator accumulator;
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      accumulator.Add(static_cast<long double>(jacobian_determinant) *
                      static_cast<long double>(integral.curl_curl[u][v]) *
                      static_cast<long double>(Dot(curl_geometry[u], curl_geometry[v])));
    }
  }
  return accumulator.Value();
}

double ContractFirstOrderMass(const FirstOrderReferenceIntegral &integral,
                              const BarycentricGradients &grad_lambda,
                              double jacobian_determinant)
{
  return ContractMass({integral.row_basis, integral.column_basis,
                       integral.input_to_canonical, integral.quadrature_rule,
                       integral.quadrature_order, integral.subdivisions,
                       integral.radial_power, integral.mass, integral.curl_curl},
                      grad_lambda, jacobian_determinant);
}

double ContractFirstOrderCurlCurl(const FirstOrderReferenceIntegral &integral,
                                  const BarycentricGradients &grad_lambda,
                                  double jacobian_determinant)
{
  return ContractCurlCurl({integral.row_basis, integral.column_basis,
                           integral.input_to_canonical, integral.quadrature_rule,
                           integral.quadrature_order, integral.subdivisions,
                           integral.radial_power, integral.mass, integral.curl_curl},
                          grad_lambda, jacobian_determinant);
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
