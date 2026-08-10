// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularassembly.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>
#include <unistd.h>
#include <fmt/format.h>

#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

double Dot(const Vector3 &a, const Vector3 &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Dot(const Vector2 &a, const Vector2 &b)
{
  return a[0] * b[0] + a[1] * b[1];
}

Vector3 Cross(const Vector3 &a, const Vector3 &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

bool SameBasis(const HigherOrderBasis &a, const HigherOrderBasis &b)
{
  return a.family == b.family && a.nodes == b.nodes &&
         a.interpolation_indices == b.interpolation_indices && a.order == b.order &&
         a.nu == b.nu;
}

bool SameBasis(const TriangleBasis &a, const TriangleBasis &b)
{
  return a.family == b.family && a.nodes == b.nodes && a.order == b.order && a.nu == b.nu;
}

bool SameBasis(const FirstOrderBasis &a, const FirstOrderBasis &b)
{
  return a.family == b.family && a.indices == b.indices && a.nu == b.nu;
}

bool SameBasis(const ReferenceBasis &a, const ReferenceBasis &b)
{
  if (a.index() != b.index())
  {
    return false;
  }
  return std::visit(
      [](const auto &left, const auto &right)
      {
        using Left = std::decay_t<decltype(left)>;
        using Right = std::decay_t<decltype(right)>;
        if constexpr (std::is_same_v<Left, Right>)
        {
          return SameBasis(left, right);
        }
        else
        {
          return false;
        }
      },
      a, b);
}

enum class DuffyFeatureKind
{
  NONE,
  NODE,
  EDGE
};

struct DuffyFeature
{
  DuffyFeatureKind kind = DuffyFeatureKind::NONE;
  std::array<int, 2> nodes{-1, -1};
};

DuffyFeature GetDuffyFeature(const FirstOrderBasis &basis)
{
  switch (basis.family)
  {
    case FirstOrderBasisFamily::NODE_GRADIENT:
    case FirstOrderBasisFamily::NODE_ROTATIONAL:
      return {DuffyFeatureKind::NODE, {basis.indices[0], -1}};
    case FirstOrderBasisFamily::EDGE_GRADIENT:
    case FirstOrderBasisFamily::EDGE_ROTATIONAL:
      {
        std::array<int, 2> edge{basis.indices[0], basis.indices[1]};
        std::sort(edge.begin(), edge.end());
        return {DuffyFeatureKind::EDGE, edge};
      }
    case FirstOrderBasisFamily::STANDARD_H1_GRADIENT:
    case FirstOrderBasisFamily::STANDARD_NEDELEC:
      return {};
  }
  throw std::invalid_argument("Unknown first-order singular-element basis family!");
}

DuffyFeature GetDuffyFeature(const HigherOrderBasis &basis)
{
  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
      return {DuffyFeatureKind::NODE, {basis.nodes[0], -1}};
    case HigherOrderBasisFamily::EDGE_GRADIENT:
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      {
        std::array<int, 2> edge{basis.nodes[0], basis.nodes[1]};
        std::sort(edge.begin(), edge.end());
        return {DuffyFeatureKind::EDGE, edge};
      }
  }
  throw std::invalid_argument("Unknown higher-order singular-element basis family!");
}

DuffyFeature GetDuffyFeature(const ReferenceBasis &basis)
{
  return std::visit([](const auto &entry) { return GetDuffyFeature(entry); }, basis);
}

std::optional<DuffyFeature> GetCommonDuffyFeature(const ReferenceBasis &row_basis,
                                                  const ReferenceBasis &column_basis)
{
  DuffyFeature common;
  for (const auto *basis : {&row_basis, &column_basis})
  {
    const auto feature = GetDuffyFeature(*basis);
    if (feature.kind == DuffyFeatureKind::NONE)
    {
      continue;
    }
    if (common.kind == DuffyFeatureKind::NONE)
    {
      common = feature;
    }
    else if (common.kind != feature.kind || common.nodes != feature.nodes)
    {
      return std::nullopt;
    }
  }
  return (common.kind == DuffyFeatureKind::NONE) ? std::nullopt
                                                 : std::optional<DuffyFeature>(common);
}

struct DuffyH1ReferenceIntegral
{
  ReferenceIntegral integral;
  ReferenceCoefficientTensor estimated_absolute_error;
};

class DuffyH1ReferenceTable
{
private:
  struct Entry
  {
    ReferenceBasis row_basis;
    ReferenceBasis column_basis;
    ReferenceIntegral integral;
    ReferenceCoefficientTensor estimated_absolute_error;
  };

  std::vector<Entry> entries;
  std::size_t hits = 0;

public:
  std::optional<DuffyH1ReferenceIntegral> Get(const ReferenceBasis &row_basis,
                                              const ReferenceBasis &column_basis)
  {
    const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
    const auto feature = GetCommonDuffyFeature(canonical.row_basis, canonical.column_basis);
    const auto row_feature = GetDuffyFeature(canonical.row_basis);
    const auto column_feature = GetDuffyFeature(canonical.column_basis);
    const bool partitioned = !feature && row_feature.kind != DuffyFeatureKind::NONE &&
                             column_feature.kind != DuffyFeatureKind::NONE;
    if (!feature && !partitioned)
    {
      return std::nullopt;
    }

    for (const auto &entry : entries)
    {
      if (SameBasis(entry.row_basis, canonical.row_basis) &&
          SameBasis(entry.column_basis, canonical.column_basis))
      {
        hits++;
        auto integral = entry.integral;
        integral.input_to_canonical = canonical.input_to_canonical;
        return DuffyH1ReferenceIntegral{std::move(integral),
                                        entry.estimated_absolute_error};
      }
    }

    const auto integrate = [&](int order)
    {
      if (feature && feature->kind == DuffyFeatureKind::NODE)
      {
        return ComputeNodeDuffyReferenceIntegral(
            canonical.row_basis, canonical.column_basis, order, H1DuffyRadialPower);
      }
      if (feature)
      {
        return ComputeEdgeDuffyReferenceIntegral(
            canonical.row_basis, canonical.column_basis, order, H1DuffyRadialPower);
      }
      return ComputePartitionedDuffyReferenceIntegral(
          canonical.row_basis, canonical.column_basis, order, H1DuffyRadialPower);
    };
    auto integral = integrate(H1DuffyReferenceOrder);
    const auto comparison = integrate(H1DuffyComparisonOrder);
    if (!SameBasis(integral.row_basis, canonical.row_basis) ||
        !SameBasis(integral.column_basis, canonical.column_basis) ||
        !SameBasis(comparison.row_basis, canonical.row_basis) ||
        !SameBasis(comparison.column_basis, canonical.column_basis))
    {
      throw std::logic_error(
          "Singular H1 Duffy reference table canonicalization is inconsistent!");
    }

    ReferenceCoefficientTensor estimated_absolute_error{};
    for (int u = 0; u < 3; u++)
    {
      for (int v = 0; v < 3; v++)
      {
        const double value = integral.mass[u][v];
        const double comparison_value = comparison.mass[u][v];
        if (!std::isfinite(value) || !std::isfinite(comparison_value))
        {
          throw std::runtime_error(
              "Singular H1 Duffy reference table produced a nonfinite entry!");
        }
        const double scale = std::max({1.0, std::abs(value), std::abs(comparison_value)});
        estimated_absolute_error[u][v] =
            H1DuffyErrorSafetyFactor * std::abs(value - comparison_value) +
            64.0 * std::numeric_limits<double>::epsilon() * scale;
      }
    }

    entries.push_back(
        {canonical.row_basis, canonical.column_basis, integral, estimated_absolute_error});
    integral.input_to_canonical = canonical.input_to_canonical;
    return DuffyH1ReferenceIntegral{std::move(integral), estimated_absolute_error};
  }

  std::size_t Size() const { return entries.size(); }
  std::size_t Hits() const { return hits; }
};

std::size_t FixedSubdivisionLeafCount(int subdivisions)
{
  std::size_t leaves = 1;
  for (int level = 0; level < subdivisions; level++)
  {
    leaves *= 8;
  }
  return leaves;
}

AdaptiveQuadratureResult IntegrateConfigured(int order,
                                             const AdaptiveAssemblyOptions &options,
                                             const ReferenceIntegrand &integrand)
{
  if (!options.fixed_subdivision)
  {
    return IntegrateReferenceTetrahedronAdaptive(order, options.absolute_tolerance,
                                                 options.relative_tolerance,
                                                 options.maximum_subdivisions, integrand);
  }
  return {IntegrateReferenceTetrahedron(order, options.subdivisions, integrand), 0.0,
          FixedSubdivisionLeafCount(options.subdivisions), options.subdivisions, true};
}

AdaptiveVectorQuadratureResult
IntegrateConfigured(int order, const AdaptiveAssemblyOptions &options,
                    std::size_t number_components,
                    const ReferenceVectorIntegrand &integrand)
{
  if (!options.fixed_subdivision)
  {
    return IntegrateReferenceTetrahedronAdaptive(
        order, options.absolute_tolerance, options.relative_tolerance,
        options.maximum_subdivisions, number_components, integrand);
  }
  return {IntegrateReferenceTetrahedron(order, options.subdivisions, number_components,
                                        integrand),
          std::vector<double>(number_components),
          FixedSubdivisionLeafCount(options.subdivisions), options.subdivisions, true};
}

AdaptiveReferenceIntegral
ComputeReferenceIntegralConfigured(const ReferenceBasis &row_basis,
                                   const ReferenceBasis &column_basis,
                                   const AdaptiveAssemblyOptions &options)
{
  if (!options.fixed_subdivision)
  {
    return ComputeAdaptiveReferenceIntegral(
        row_basis, column_basis, options.quadrature_order, options.absolute_tolerance,
        options.relative_tolerance, options.maximum_subdivisions);
  }
  constexpr std::size_t tensor_entries = 18;
  return {ComputeReferenceIntegral(row_basis, column_basis, options.quadrature_order,
                                   options.subdivisions),
          0.0,
          0.0,
          {},
          {},
          tensor_entries * FixedSubdivisionLeafCount(options.subdivisions),
          options.subdivisions,
          true};
}

class AdaptiveReferenceTable
{
private:
  struct Entry
  {
    ReferenceBasis row_basis;
    ReferenceBasis column_basis;
    AdaptiveReferenceIntegral reference;
  };

  AdaptiveAssemblyOptions options;
  std::vector<Entry> entries;
  std::size_t hits = 0;

public:
  explicit AdaptiveReferenceTable(const AdaptiveAssemblyOptions &options_in)
    : options(options_in)
  {
  }

  AdaptiveReferenceIntegral Get(const ReferenceBasis &row_basis,
                                const ReferenceBasis &column_basis)
  {
    const auto canonical = CanonicalizeReferenceBasisPair(row_basis, column_basis);
    for (const auto &entry : entries)
    {
      if (SameBasis(entry.row_basis, canonical.row_basis) &&
          SameBasis(entry.column_basis, canonical.column_basis))
      {
        hits++;
        auto reference = entry.reference;
        reference.integral.input_to_canonical = canonical.input_to_canonical;
        return reference;
      }
    }

    auto reference = ComputeReferenceIntegralConfigured(canonical.row_basis,
                                                        canonical.column_basis, options);
    if (!SameBasis(reference.integral.row_basis, canonical.row_basis) ||
        !SameBasis(reference.integral.column_basis, canonical.column_basis))
    {
      throw std::logic_error(
          "Singular adaptive reference table canonicalization is inconsistent!");
    }
    entries.push_back({canonical.row_basis, canonical.column_basis, reference});
    reference.integral.input_to_canonical = canonical.input_to_canonical;
    return reference;
  }

  std::size_t Size() const { return entries.size(); }
  std::size_t Hits() const { return hits; }
};

void ValidateAdaptiveAssemblyOptions(const AdaptiveAssemblyOptions &options)
{
  const bool valid_fixed = options.fixed_subdivision && options.quadrature_order >= 1 &&
                           options.quadrature_order <= 21 && options.subdivisions >= 0 &&
                           options.subdivisions <= 8;
  const bool valid_adaptive =
      !options.fixed_subdivision && options.quadrature_order >= 1 &&
      std::isfinite(options.absolute_tolerance) && options.absolute_tolerance >= 0.0 &&
      std::isfinite(options.relative_tolerance) && options.relative_tolerance >= 0.0 &&
      (options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) &&
      options.maximum_subdivisions >= 1;
  if (!valid_fixed && !valid_adaptive)
  {
    throw std::invalid_argument(
        "Singular element assembly has invalid quadrature options!");
  }
}

TriangleVectorBasisValue
EvaluateTriangleBasis(const TriangleBarycentricPoint &lambda,
                      const TriangleBarycentricGradients &grad_lambda,
                      const TriangleBasis &basis)
{
  if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
  {
    return EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0], basis.nodes[1],
                                        basis.nu);
  }
  if (basis.family == HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    return EvaluateTriangleNodeRotational(lambda, grad_lambda, basis.nodes[0],
                                          basis.nodes[1], basis.nodes[2], basis.nu);
  }
  throw std::invalid_argument("Unsupported triangular singular basis family!");
}

Vector2 EvaluateTriangleBasisValue(const TriangleBarycentricPoint &lambda,
                                   const TriangleBarycentricGradients &grad_lambda,
                                   const TriangleBasis &basis)
{
  if (basis.family != HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    return EvaluateTriangleBasis(lambda, grad_lambda, basis).value;
  }
  constexpr double endpoint_tolerance = 128.0 * std::numeric_limits<double>::epsilon();
  if (TriangleNodeRadialCoordinate(lambda, basis.nodes[0]) <= endpoint_tolerance)
  {
    // The rotational field vanishes continuously at the singular node. Do not
    // evaluate its singular curl when constructing a value-only interpolant.
    return {0.0, 0.0};
  }
  return EvaluateTriangleNodeRotational(lambda, grad_lambda, basis.nodes[0], basis.nodes[1],
                                        basis.nodes[2], basis.nu)
      .value;
}

int TriangleSingularNode(const TriangleBasis &basis)
{
  if (basis.family != HigherOrderBasisFamily::NODE_GRADIENT &&
      basis.family != HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    throw std::invalid_argument("Unsupported triangular singular basis family!");
  }
  return basis.nodes[0];
}

void ValidateTriangleInputs(const TriangleElementDofMap &element_dofs,
                            const TriangleBarycentricGradients &grad_lambda,
                            double jacobian_determinant,
                            const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  if (!std::isfinite(jacobian_determinant) || !(jacobian_determinant > 0.0))
  {
    throw std::invalid_argument(
        "Triangular singular assembly requires a positive affine area Jacobian!");
  }
  Vector2 gradient_sum{};
  double gradient_scale = 1.0;
  for (const auto &gradient : grad_lambda)
  {
    for (int d = 0; d < 2; d++)
    {
      if (!std::isfinite(gradient[d]))
      {
        throw std::invalid_argument(
            "Triangular singular assembly received nonfinite barycentric gradients!");
      }
      gradient_sum[d] += gradient[d];
      gradient_scale = std::max(gradient_scale, std::abs(gradient[d]));
    }
  }
  const double inverse_jacobian_determinant = std::abs(
      grad_lambda[1][0] * grad_lambda[2][1] - grad_lambda[1][1] * grad_lambda[2][0]);
  const double determinant_product = jacobian_determinant * inverse_jacobian_determinant;
  for (double value : gradient_sum)
  {
    if (std::abs(value) > 1.0e-12 * gradient_scale)
    {
      throw std::invalid_argument(
          "Triangular singular barycentric gradients do not sum to zero!");
    }
  }
  if (!(inverse_jacobian_determinant > 0.0) ||
      std::abs(determinant_product - 1.0) > 1.0e-10 * std::max(1.0, determinant_product))
  {
    throw std::invalid_argument(
        "Triangular singular gradients and Jacobian are inconsistent!");
  }

  const auto validate_basis = [](const TriangleBasis &basis)
  {
    auto sorted = basis.nodes;
    std::sort(sorted.begin(), sorted.end());
    if (basis.order != 1 || sorted != std::array<int, 3>{0, 1, 2} ||
        !std::isfinite(basis.nu) || basis.nu <= 0.0 || basis.nu >= 1.0)
    {
      throw std::invalid_argument(
          "Triangular singular assembly received an invalid basis descriptor!");
    }
    (void)TriangleSingularNode(basis);
  };
  std::set<std::size_t> h1_dofs, nd_dofs;
  for (const auto &dof : element_dofs.h1)
  {
    validate_basis(dof.basis);
    if (!h1_dofs.insert(dof.dof).second ||
        dof.basis.family != HigherOrderBasisFamily::NODE_GRADIENT)
    {
      throw std::invalid_argument("Triangular singular H1 basis map is inconsistent!");
    }
  }
  for (const auto &dof : element_dofs.nd)
  {
    validate_basis(dof.basis);
    if (!nd_dofs.insert(dof.dof).second)
    {
      throw std::invalid_argument(
          "Triangular singular ND basis map contains duplicate DOFs!");
    }
  }
}

std::vector<int> BuildTriangleH1ToNDMap(const TriangleElementDofMap &element_dofs)
{
  std::vector<int> h1_to_nd(element_dofs.h1.size(), -1);
  for (std::size_t h1 = 0; h1 < element_dofs.h1.size(); h1++)
  {
    for (std::size_t nd = 0; nd < element_dofs.nd.size(); nd++)
    {
      if (SameBasis(element_dofs.h1[h1].basis, element_dofs.nd[nd].basis))
      {
        if (h1_to_nd[h1] >= 0)
        {
          throw std::invalid_argument(
              "Triangular singular H1 basis has multiple matching ND bases!");
        }
        h1_to_nd[h1] = static_cast<int>(nd);
      }
    }
    if (h1_to_nd[h1] < 0)
    {
      throw std::invalid_argument(
          "Triangular singular H1 basis is absent from the ND basis map!");
    }
  }
  return h1_to_nd;
}

struct TriangleQuadratureResult
{
  double value;
  double estimated_absolute_error;
  std::size_t point_count;
};

double TriangleDuffyPartitionWeight(const TriangleBarycentricPoint &lambda,
                                    int aligned_node, int other_node)
{
  const double aligned_rho = TriangleNodeRadialCoordinate(lambda, aligned_node);
  const double other_rho = TriangleNodeRadialCoordinate(lambda, other_node);
  if (!(aligned_rho > 0.0) || !(other_rho > 0.0))
  {
    throw std::domain_error(
        "Partitioned triangle Duffy quadrature requires interior points!");
  }
  if (aligned_rho <= other_rho)
  {
    const double ratio = aligned_rho / other_rho;
    return 1.0 / (1.0 + std::pow(ratio, MultiFeatureDuffyPartitionPower));
  }
  const double ratio = other_rho / aligned_rho;
  const double power = std::pow(ratio, MultiFeatureDuffyPartitionPower);
  return power / (1.0 + power);
}

template <typename Integrand>
TriangleQuadratureResult
IntegrateTriangleDuffy(int row_node, int column_node, const char *quantity, int matrix_row,
                       int matrix_column, const AdaptiveAssemblyOptions &options,
                       Integrand &&integrand)
{
  if (!quantity)
  {
    throw std::invalid_argument(
        "Triangle Duffy quadrature requires a matrix quantity description!");
  }
  const int initial_high_order =
      std::max(H1DuffyReferenceOrder, 2 * options.quadrature_order + 15);
  constexpr int order_increment = 8;
  const bool partitioned = row_node >= 0 && column_node >= 0 && row_node != column_node;
  const int aligned_node = row_node >= 0 ? row_node : column_node;
  if (aligned_node < 0 || aligned_node >= 3)
  {
    throw std::invalid_argument(
        "Triangle Duffy quadrature requires at least one singular basis!");
  }

  const auto integrate = [&](int order)
  {
    if (!partitioned)
    {
      return IntegrateReferenceTriangleNodeDuffy(order, aligned_node,
                                                 TriangleDuffyRadialPower, integrand);
    }
    const auto chart = [&](int node, int other)
    {
      return IntegrateReferenceTriangleNodeDuffy(
          order, node, TriangleDuffyRadialPower,
          [&](const TriangleBarycentricPoint &lambda)
          {
            return TriangleDuffyPartitionWeight(lambda, node, other) * integrand(lambda);
          });
    };
    return chart(row_node, column_node) + chart(column_node, row_node);
  };

  const std::size_t charts = partitioned ? 2 : 1;
  int comparison_order = initial_high_order - order_increment;
  double comparison = integrate(comparison_order);
  std::size_t point_count =
      charts * ReferenceTriangleNodeDuffyQuadraturePointCount(comparison_order);
  double value = std::numeric_limits<double>::quiet_NaN();
  double error = std::numeric_limits<double>::infinity();
  double tolerance = 0.0;
  int high_order = initial_high_order;
  for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
  {
    value = integrate(high_order);
    point_count += charts * ReferenceTriangleNodeDuffyQuadraturePointCount(high_order);
    const double scale = std::max({1.0, std::abs(value), std::abs(comparison)});
    error = H1DuffyErrorSafetyFactor * std::abs(value - comparison) +
            64.0 * std::numeric_limits<double>::epsilon() * scale;
    tolerance = options.absolute_tolerance + options.relative_tolerance * std::abs(value);
    if (std::isfinite(value) && std::isfinite(comparison) && std::isfinite(error) &&
        error >= 0.0 && error <= tolerance)
    {
      return {value, error, point_count};
    }
    if (refinement == options.maximum_subdivisions)
    {
      break;
    }
    comparison = value;
    comparison_order = high_order;
    high_order += order_increment;
  }
  throw std::runtime_error(fmt::format(
      "Triangular singular Duffy integral did not meet tolerance: value = {:.17g}, "
      "comparison = {:.17g}, estimated absolute error = {:.17g}, tolerance = "
      "{:.17g}, quantity = {}, row/column = {}/{}, singular nodes = {}/{}, partitioned "
      "= {}, orders = {}/{}, order refinements = {}!",
      value, comparison, error, tolerance, quantity, matrix_row, matrix_column, row_node,
      column_node, partitioned, high_order, comparison_order,
      options.maximum_subdivisions));
}

void ValidateInputs(const ElementDofMap &element_dofs,
                    const BarycentricGradients &grad_lambda, double jacobian_determinant,
                    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  if (!std::isfinite(jacobian_determinant) || !(jacobian_determinant > 0.0))
  {
    throw std::invalid_argument(
        "Singular element assembly requires a positive affine Jacobian determinant!");
  }

  Vector3 gradient_sum{};
  double gradient_scale = 1.0;
  for (const auto &gradient : grad_lambda)
  {
    for (int d = 0; d < 3; d++)
    {
      if (!std::isfinite(gradient[d]))
      {
        throw std::invalid_argument(
            "Singular element assembly received nonfinite barycentric gradients!");
      }
      gradient_sum[d] += gradient[d];
      gradient_scale = std::max(gradient_scale, std::abs(gradient[d]));
    }
  }
  const double inverse_jacobian_determinant =
      std::abs(Dot(grad_lambda[1], Cross(grad_lambda[2], grad_lambda[3])));
  const double determinant_product = jacobian_determinant * inverse_jacobian_determinant;
  for (double value : gradient_sum)
  {
    if (std::abs(value) > 1.0e-12 * gradient_scale)
    {
      throw std::invalid_argument(
          "Singular element barycentric gradients do not sum to zero!");
    }
  }
  if (!(inverse_jacobian_determinant > 0.0) ||
      std::abs(determinant_product - 1.0) > 1.0e-10 * std::max(1.0, determinant_product))
  {
    throw std::invalid_argument(
        "Singular element gradients and Jacobian determinant are inconsistent!");
  }

  std::set<std::size_t> h1_dofs, nd_dofs;
  for (const auto &dof : element_dofs.h1)
  {
    if (!h1_dofs.insert(dof.dof).second || !IsGradientFamily(dof.basis.family))
    {
      throw std::invalid_argument("Singular element H1 basis map is inconsistent!");
    }
  }
  for (const auto &dof : element_dofs.nd)
  {
    if (!nd_dofs.insert(dof.dof).second)
    {
      throw std::invalid_argument("Singular element ND basis map contains duplicate DOFs!");
    }
  }
}

std::vector<int> BuildH1ToNDMap(const ElementDofMap &element_dofs)
{
  std::vector<int> h1_to_nd(element_dofs.h1.size(), -1);
  for (std::size_t h1 = 0; h1 < element_dofs.h1.size(); h1++)
  {
    for (std::size_t nd = 0; nd < element_dofs.nd.size(); nd++)
    {
      if (SameBasis(element_dofs.h1[h1].basis, element_dofs.nd[nd].basis))
      {
        if (h1_to_nd[h1] >= 0)
        {
          throw std::invalid_argument(
              "Singular element H1 basis has multiple matching ND bases!");
        }
        h1_to_nd[h1] = static_cast<int>(nd);
      }
    }
    if (h1_to_nd[h1] < 0)
    {
      throw std::invalid_argument(
          "Singular element H1 basis is absent from the ND basis map!");
    }
  }
  return h1_to_nd;
}

double ContractMassError(const ReferenceIntegral &integral,
                         const ReferenceCoefficientTensor &estimated_absolute_error,
                         const BarycentricGradients &grad_lambda,
                         double jacobian_determinant)
{
  const auto canonical_grad_lambda =
      ApplyBarycentricPermutation(grad_lambda, integral.input_to_canonical);
  long double error = 0.0L;
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      error += static_cast<long double>(jacobian_determinant) *
               static_cast<long double>(estimated_absolute_error[u][v]) *
               std::abs(static_cast<long double>(
                   Dot(canonical_grad_lambda[u + 1], canonical_grad_lambda[v + 1])));
    }
  }
  return static_cast<double>(error);
}

double ContractMassError(const AdaptiveReferenceIntegral &reference,
                         const BarycentricGradients &grad_lambda,
                         double jacobian_determinant)
{
  return ContractMassError(reference.integral, reference.mass_estimated_absolute_error,
                           grad_lambda, jacobian_determinant);
}

double ContractCurlCurlError(const AdaptiveReferenceIntegral &reference,
                             const BarycentricGradients &grad_lambda,
                             double jacobian_determinant)
{
  const auto canonical_grad_lambda =
      ApplyBarycentricPermutation(grad_lambda, reference.integral.input_to_canonical);
  const std::array<Vector3, 3> curl_geometry{
      Cross(canonical_grad_lambda[2], canonical_grad_lambda[3]),
      Cross(canonical_grad_lambda[3], canonical_grad_lambda[1]),
      Cross(canonical_grad_lambda[1], canonical_grad_lambda[2])};
  long double error = 0.0L;
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      error +=
          static_cast<long double>(jacobian_determinant) *
          static_cast<long double>(reference.curl_curl_estimated_absolute_error[u][v]) *
          std::abs(static_cast<long double>(Dot(curl_geometry[u], curl_geometry[v])));
    }
  }
  return static_cast<double>(error);
}

void RecordQuadratureStatistics(const AdaptiveReferenceIntegral &reference,
                                ElementEnrichmentMatrices &result)
{
  if (result.total_quadrature_leaf_count >
      std::numeric_limits<std::size_t>::max() - reference.total_leaf_count)
  {
    throw std::overflow_error("Singular element assembly quadrature leaf count overflow!");
  }
  result.total_quadrature_leaf_count += reference.total_leaf_count;
  result.maximum_subdivision_depth =
      std::max(result.maximum_subdivision_depth, reference.maximum_subdivision_depth);
}

template <typename T>
void RecordQuadratureStatistics(const AdaptiveQuadratureResult &integral, T &result)
{
  if (result.total_quadrature_leaf_count >
      std::numeric_limits<std::size_t>::max() - integral.leaf_count)
  {
    throw std::overflow_error("Singular element assembly quadrature leaf count overflow!");
  }
  result.total_quadrature_leaf_count += integral.leaf_count;
  result.maximum_subdivision_depth =
      std::max(result.maximum_subdivision_depth, integral.maximum_subdivision_depth);
}

template <typename T>
void AccumulateQuadratureStatistics(const T &element, std::size_t &total_leaf_count,
                                    int &maximum_subdivision_depth)
{
  if (total_leaf_count >
      std::numeric_limits<std::size_t>::max() - element.total_quadrature_leaf_count)
  {
    throw std::overflow_error("Singular element assembly quadrature leaf count overflow!");
  }
  total_leaf_count += element.total_quadrature_leaf_count;
  maximum_subdivision_depth =
      std::max(maximum_subdivision_depth, element.maximum_subdivision_depth);
}

void ValidateStandardRowTransformation(const mfem::DofTransformation &dof_transformation,
                                       const mfem::DenseMatrix &standard_enrichment,
                                       const mfem::DenseMatrix &enrichment_standard,
                                       const mfem::DenseMatrix &estimated_absolute_error)
{
  const int standard_size = standard_enrichment.Height();
  const int enrichment_size = standard_enrichment.Width();
  if (enrichment_standard.Height() != enrichment_size ||
      enrichment_standard.Width() != standard_size ||
      estimated_absolute_error.Height() != standard_size ||
      estimated_absolute_error.Width() != enrichment_size)
  {
    throw std::invalid_argument(
        "Singular element coupling blocks have inconsistent dimensions!");
  }
  if (const auto *transformation = dof_transformation.GetDofTransformation();
      transformation && transformation->Size() != standard_size)
  {
    throw std::invalid_argument(
        "Singular element coupling block does not match its MFEM DOF transformation!");
  }
}

void ApplyStandardRowTransformation(const mfem::DofTransformation &dof_transformation,
                                    mfem::DenseMatrix &standard_enrichment,
                                    mfem::DenseMatrix &enrichment_standard,
                                    mfem::DenseMatrix &estimated_absolute_error)
{
  const int standard_size = standard_enrichment.Height();
  const int enrichment_size = standard_enrichment.Width();
  if (!dof_transformation.IsIdentity())
  {
    mfem::DenseMatrix dual_transformation(standard_size);
    dual_transformation = 0.0;
    for (int i = 0; i < standard_size; i++)
    {
      dual_transformation(i, i) = 1.0;
    }
    dof_transformation.TransformDualCols(dual_transformation);
    dof_transformation.TransformDualCols(standard_enrichment);

    mfem::DenseMatrix transformed_error(standard_size, enrichment_size);
    for (int i = 0; i < standard_size; i++)
    {
      for (int j = 0; j < enrichment_size; j++)
      {
        long double error = 0.0L;
        for (int k = 0; k < standard_size; k++)
        {
          error += std::abs(static_cast<long double>(dual_transformation(i, k))) *
                   static_cast<long double>(estimated_absolute_error(k, j));
        }
        transformed_error(i, j) = static_cast<double>(error);
      }
    }
    estimated_absolute_error.Swap(transformed_error);
  }
  enrichment_standard.Transpose(standard_enrichment);
}

void ValidateMaterialCoefficients(const IsotropicMaterialCoefficients &coefficients)
{
  if (!std::isfinite(coefficients.electric) ||
      !std::isfinite(coefficients.inverse_magnetic) ||
      !(coefficients.inverse_magnetic > 0.0))
  {
    throw std::invalid_argument(
        "Singular element electric coefficient must be finite and its "
        "inverse-magnetic coefficient must be finite and positive!");
  }
}

void ValidateMatrixForScaling(const mfem::DenseMatrix &matrix, double coefficient,
                              bool estimated_absolute_error)
{
  const double magnitude = std::abs(coefficient);
  const double maximum_unscaled = magnitude > 0.0
                                      ? std::numeric_limits<double>::max() / magnitude
                                      : std::numeric_limits<double>::infinity();
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      const double value = matrix(i, j);
      if (!std::isfinite(value) || (estimated_absolute_error && value < 0.0))
      {
        throw std::invalid_argument(
            "Singular element material scaling received an invalid matrix entry!");
      }
      if (std::abs(value) > maximum_unscaled)
      {
        throw std::overflow_error("Singular element material scaling would overflow!");
      }
    }
  }
}

void ValidateSymmetricMatrixForScaling(const mfem::DenseMatrix &matrix,
                                       const mfem::DenseMatrix &estimated_absolute_error,
                                       double coefficient)
{
  if (matrix.Height() != matrix.Width() ||
      estimated_absolute_error.Height() != matrix.Height() ||
      estimated_absolute_error.Width() != matrix.Width())
  {
    throw std::invalid_argument(
        "Singular element material scaling requires matching square matrices!");
  }
  ValidateMatrixForScaling(matrix, coefficient, false);
  ValidateMatrixForScaling(estimated_absolute_error, coefficient, true);
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = i + 1; j < matrix.Width(); j++)
    {
      if (matrix(i, j) != matrix(j, i) ||
          estimated_absolute_error(i, j) != estimated_absolute_error(j, i))
      {
        throw std::invalid_argument(
            "Singular element material scaling requires symmetric matrices!");
      }
    }
  }
}

void ValidateCouplingMatricesForScaling(const mfem::DenseMatrix &standard_enrichment,
                                        const mfem::DenseMatrix &enrichment_standard,
                                        const mfem::DenseMatrix &estimated_absolute_error,
                                        double coefficient)
{
  if (enrichment_standard.Height() != standard_enrichment.Width() ||
      enrichment_standard.Width() != standard_enrichment.Height() ||
      estimated_absolute_error.Height() != standard_enrichment.Height() ||
      estimated_absolute_error.Width() != standard_enrichment.Width())
  {
    throw std::invalid_argument(
        "Singular element material scaling received inconsistent coupling blocks!");
  }
  ValidateMatrixForScaling(standard_enrichment, coefficient, false);
  ValidateMatrixForScaling(enrichment_standard, coefficient, false);
  ValidateMatrixForScaling(estimated_absolute_error, coefficient, true);
  for (int i = 0; i < standard_enrichment.Height(); i++)
  {
    for (int j = 0; j < standard_enrichment.Width(); j++)
    {
      if (standard_enrichment(i, j) != enrichment_standard(j, i))
      {
        throw std::invalid_argument(
            "Singular element material scaling requires exact transpose blocks!");
      }
    }
  }
}

void ScaleMatrix(mfem::DenseMatrix &matrix, double coefficient)
{
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      matrix(i, j) *= coefficient;
    }
  }
}

void ScaleSymmetricMatrices(mfem::DenseMatrix &matrix,
                            mfem::DenseMatrix &estimated_absolute_error, double coefficient)
{
  ScaleMatrix(matrix, coefficient);
  ScaleMatrix(estimated_absolute_error, std::abs(coefficient));
}

void ScaleCouplingMatrices(mfem::DenseMatrix &standard_enrichment,
                           mfem::DenseMatrix &enrichment_standard,
                           mfem::DenseMatrix &estimated_absolute_error, double coefficient)
{
  ScaleMatrix(standard_enrichment, coefficient);
  enrichment_standard.Transpose(standard_enrichment);
  ScaleMatrix(estimated_absolute_error, std::abs(coefficient));
}

template <typename ElementDofType>
mfem::Array<int> GetElementEnrichmentDofs(const std::vector<ElementDofType> &element_dofs,
                                          std::size_t local_size)
{
  if (element_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::overflow_error("Singular element local DOF count exceeds integer range!");
  }
  mfem::Array<int> result(static_cast<int>(element_dofs.size()));
  std::set<std::size_t> unique_dofs;
  for (int i = 0; i < result.Size(); i++)
  {
    const std::size_t dof = element_dofs[i].dof;
    if (dof >= local_size ||
        dof > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        !unique_dofs.insert(dof).second)
    {
      throw std::invalid_argument(
          "Singular element sparse assembly received an invalid local DOF map!");
    }
    result[i] = static_cast<int>(dof);
  }
  return result;
}

mfem::Array<int> UnsignedDofs(const mfem::Array<int> &dofs)
{
  mfem::Array<int> result(dofs.Size());
  for (int i = 0; i < dofs.Size(); i++)
  {
    if (dofs[i] == std::numeric_limits<int>::min())
    {
      throw std::invalid_argument("Singular sparse assembly received an invalid VDof!");
    }
    result[i] = (dofs[i] >= 0) ? dofs[i] : -1 - dofs[i];
  }
  return result;
}

void InitializeLocalSparseBlock(LocalSparseOperatorBlocks &block, int standard_size,
                                int enrichment_size)
{
  block.enrichment_enrichment =
      std::make_unique<mfem::SparseMatrix>(enrichment_size, enrichment_size);
  block.standard_enrichment =
      std::make_unique<mfem::SparseMatrix>(standard_size, enrichment_size);
  block.enrichment_enrichment_estimated_absolute_error =
      std::make_unique<mfem::SparseMatrix>(enrichment_size, enrichment_size);
  block.standard_enrichment_estimated_absolute_error =
      std::make_unique<mfem::SparseMatrix>(standard_size, enrichment_size);
}

void FinalizeLocalSparseBlock(LocalSparseOperatorBlocks &block)
{
  block.enrichment_enrichment->Finalize();
  block.standard_enrichment->Finalize();
  block.enrichment_enrichment_estimated_absolute_error->Finalize();
  block.standard_enrichment_estimated_absolute_error->Finalize();
  block.enrichment_standard.reset(mfem::Transpose(*block.standard_enrichment));
}

std::vector<HYPRE_BigInt> BuildPartition(MPI_Comm comm, HYPRE_BigInt local_offset,
                                         HYPRE_BigInt local_size, HYPRE_BigInt global_size)
{
  const int ranks = Mpi::Size(comm);
  std::vector<HYPRE_BigInt> offsets(ranks), sizes(ranks), global_sizes(ranks);
  Mpi::Allgather(1, &local_offset, offsets.data(), comm);
  Mpi::Allgather(1, &local_size, sizes.data(), comm);
  Mpi::Allgather(1, &global_size, global_sizes.data(), comm);
  HYPRE_BigInt expected_offset = 0;
  bool valid = global_size >= 0 && std::all_of(global_sizes.begin(), global_sizes.end(),
                                               [global_size](HYPRE_BigInt size)
                                               { return size == global_size; });
  for (int rank = 0; rank < ranks; rank++)
  {
    if (offsets[rank] != expected_offset || sizes[rank] < 0 ||
        sizes[rank] > std::numeric_limits<HYPRE_BigInt>::max() - expected_offset)
    {
      valid = false;
      break;
    }
    expected_offset += sizes[rank];
  }
  if (!valid || expected_offset != global_size)
  {
    throw std::invalid_argument("Singular parallel DOF partition is inconsistent!");
  }
  if (HYPRE_AssumedPartitionCheck())
  {
    return {local_offset, local_offset + local_size};
  }

  std::vector<HYPRE_BigInt> partition(ranks + 1);
  std::copy(offsets.begin(), offsets.end(), partition.begin());
  partition.back() = global_size;
  return partition;
}

void ValidateTrueDofMap(MPI_Comm comm, const TrueDofMap &dofs)
{
  bool locally_valid =
      dofs.local_size == static_cast<HYPRE_BigInt>(dofs.local_to_true.size()) &&
      dofs.owner.size() == dofs.local_to_true.size() && dofs.local_size >= 0 &&
      dofs.local_size <= static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      dofs.global_size >= 0 && dofs.owned_offset >= 0 && dofs.owned_size >= 0 &&
      dofs.owned_offset <= dofs.global_size &&
      dofs.owned_size <= dofs.global_size - dofs.owned_offset;
  if (locally_valid)
  {
    for (std::size_t i = 0; i < dofs.local_to_true.size(); i++)
    {
      const auto column = dofs.local_to_true[i];
      const bool locally_owned =
          dofs.owned_offset <= column && column < dofs.owned_offset + dofs.owned_size;
      if (column < 0 || column >= dofs.global_size || dofs.owner[i] < 0 ||
          dofs.owner[i] >= Mpi::Size(comm) ||
          (dofs.owner[i] == Mpi::Rank(comm)) != locally_owned)
      {
        locally_valid = false;
        break;
      }
    }
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument(
        "Singular enrichment prolongation received an invalid local DOF map!");
  }
}

std::unique_ptr<mfem::HypreParMatrix> BuildEnrichmentProlongation(MPI_Comm comm,
                                                                  const TrueDofMap &dofs)
{
  ValidateTrueDofMap(comm, dofs);
  const int local_size = static_cast<int>(dofs.local_size);
  std::vector<int> rows(local_size + 1);
  std::vector<HYPRE_BigInt> columns(local_size);
  std::vector<double> values(local_size, 1.0);
  for (int i = 0; i < local_size; i++)
  {
    rows[i] = i;
    columns[i] = dofs.local_to_true[i];
  }
  rows[local_size] = local_size;
  const auto local_partition =
      BuildPartition(comm, dofs.local_offset, dofs.local_size, dofs.global_local_size);
  const auto true_partition =
      BuildPartition(comm, dofs.owned_offset, dofs.owned_size, dofs.global_size);
  return std::make_unique<mfem::HypreParMatrix>(
      comm, local_size, dofs.global_local_size, dofs.global_size, rows.data(),
      columns.data(), values.data(), local_partition.data(), true_partition.data());
}

std::unique_ptr<mfem::HypreParMatrix>
BuildAbsoluteProlongation(const mfem::HypreParMatrix &prolongation)
{
  auto *matrix = static_cast<hypre_ParCSRMatrix *>(prolongation);
  if (hypre_ParCSRMatrixMemoryLocation(matrix) != HYPRE_MEMORY_HOST)
  {
    throw std::invalid_argument(
        "Singular sparse error assembly currently requires CPU prolongation matrices!");
  }
  auto *absolute = hypre_ParCSRMatrixClone(matrix, 1);
  for (auto *block : {hypre_ParCSRMatrixDiag(absolute), hypre_ParCSRMatrixOffd(absolute)})
  {
    auto *data = hypre_CSRMatrixData(block);
    for (HYPRE_Int i = 0; i < hypre_CSRMatrixNumNonzeros(block); i++)
    {
      data[i] = std::abs(data[i]);
    }
  }
  return std::make_unique<mfem::HypreParMatrix>(absolute, true);
}

std::unique_ptr<mfem::HypreParMatrix>
AssembleParallelRAP(MPI_Comm comm, const mfem::SparseMatrix &local,
                    const mfem::HypreParMatrix &test_prolongation,
                    const mfem::HypreParMatrix &trial_prolongation)
{
  if (local.Height() != test_prolongation.Height() ||
      local.Width() != trial_prolongation.Height())
  {
    throw std::invalid_argument(
        "Singular local sparse block does not match its parallel prolongations!");
  }
  mfem::HypreParMatrix local_parallel(
      comm, test_prolongation.GetGlobalNumRows(), trial_prolongation.GetGlobalNumRows(),
      test_prolongation.GetRowStarts(), trial_prolongation.GetRowStarts(),
      const_cast<mfem::SparseMatrix *>(&local));
  return std::unique_ptr<mfem::HypreParMatrix>(
      mfem::RAP(&test_prolongation, &local_parallel, &trial_prolongation));
}

void AssembleParallelSparseBlock(MPI_Comm comm, const LocalSparseOperatorBlocks &local,
                                 const mfem::HypreParMatrix &standard_prolongation,
                                 const mfem::HypreParMatrix &absolute_standard_prolongation,
                                 const mfem::HypreParMatrix &enrichment_prolongation,
                                 ParallelSparseOperatorBlocks &parallel)
{
  if (!local.enrichment_enrichment || !local.standard_enrichment ||
      !local.enrichment_enrichment_estimated_absolute_error ||
      !local.standard_enrichment_estimated_absolute_error)
  {
    throw std::invalid_argument("Singular local sparse block is incomplete!");
  }
  parallel.enrichment_enrichment = AssembleParallelRAP(
      comm, *local.enrichment_enrichment, enrichment_prolongation, enrichment_prolongation);
  parallel.standard_enrichment = AssembleParallelRAP(
      comm, *local.standard_enrichment, standard_prolongation, enrichment_prolongation);
  parallel.enrichment_standard.reset(parallel.standard_enrichment->Transpose());
  parallel.enrichment_enrichment_estimated_absolute_error =
      AssembleParallelRAP(comm, *local.enrichment_enrichment_estimated_absolute_error,
                          enrichment_prolongation, enrichment_prolongation);
  parallel.standard_enrichment_estimated_absolute_error =
      AssembleParallelRAP(comm, *local.standard_enrichment_estimated_absolute_error,
                          absolute_standard_prolongation, enrichment_prolongation);
  if (local.transformed_enrichment_diagonal)
  {
    parallel.transformed_enrichment_diagonal =
        std::make_unique<mfem::Vector>(enrichment_prolongation.Width());
    enrichment_prolongation.MultTranspose(*local.transformed_enrichment_diagonal,
                                          *parallel.transformed_enrichment_diagonal);
  }
}

void ValidateStandardFiniteElements(const mfem::FiniteElement &h1_fe,
                                    const mfem::FiniteElement &nd_fe,
                                    const mfem::ElementTransformation &transformation)
{
  if (h1_fe.GetDim() != 3 || h1_fe.GetGeomType() != mfem::Geometry::TETRAHEDRON ||
      h1_fe.GetRangeType() != mfem::FiniteElement::SCALAR ||
      h1_fe.GetMapType() != mfem::FiniteElement::VALUE ||
      h1_fe.GetDerivType() != mfem::FiniteElement::GRAD || h1_fe.GetOrder() < 1 ||
      transformation.GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
      transformation.GetDimension() != 3 || transformation.GetSpaceDim() != 3)
  {
    throw std::invalid_argument(
        "Singular element coupling requires a three-dimensional tetrahedral H1 finite "
        "element of positive order!");
  }
  if (nd_fe.GetDim() != 3 || nd_fe.GetGeomType() != mfem::Geometry::TETRAHEDRON ||
      nd_fe.GetRangeType() != mfem::FiniteElement::VECTOR ||
      nd_fe.GetMapType() != mfem::FiniteElement::H_CURL ||
      nd_fe.GetDerivType() != mfem::FiniteElement::CURL || nd_fe.GetCurlDim() != 3)
  {
    throw std::invalid_argument(
        "Singular element coupling requires three-dimensional tetrahedral "
        "H1 and H(curl) finite elements!");
  }
  if (h1_fe.GetOrder() != nd_fe.GetOrder())
  {
    throw std::invalid_argument(
        "Singular element coupling requires matching positive H1 and ND orders!");
  }
}

void ValidateStandardTriangleFiniteElements(
    const mfem::FiniteElement &h1_fe, const mfem::FiniteElement &nd_fe,
    const mfem::ElementTransformation &transformation)
{
  if (h1_fe.GetDim() != 2 || h1_fe.GetGeomType() != mfem::Geometry::TRIANGLE ||
      h1_fe.GetRangeType() != mfem::FiniteElement::SCALAR ||
      h1_fe.GetMapType() != mfem::FiniteElement::VALUE ||
      h1_fe.GetDerivType() != mfem::FiniteElement::GRAD || h1_fe.GetOrder() < 1 ||
      transformation.GetGeometryType() != mfem::Geometry::TRIANGLE ||
      transformation.GetDimension() != 2 || transformation.GetSpaceDim() != 2)
  {
    throw std::invalid_argument(
        "Triangular singular coupling requires a two-dimensional H1 finite element of "
        "positive order!");
  }
  if (nd_fe.GetDim() != 2 || nd_fe.GetGeomType() != mfem::Geometry::TRIANGLE ||
      nd_fe.GetRangeType() != mfem::FiniteElement::VECTOR ||
      nd_fe.GetMapType() != mfem::FiniteElement::H_CURL ||
      nd_fe.GetDerivType() != mfem::FiniteElement::CURL || nd_fe.GetCurlDim() != 1)
  {
    throw std::invalid_argument(
        "Triangular singular coupling requires two-dimensional H1 and H(curl) finite "
        "elements!");
  }
  if (h1_fe.GetOrder() != nd_fe.GetOrder())
  {
    throw std::invalid_argument(
        "Triangular singular coupling requires matching H1 and ND orders!");
  }
}

Vector3 DenseMatrixRow(const mfem::DenseMatrix &matrix, int row)
{
  return {matrix(row, 0), matrix(row, 1), matrix(row, 2)};
}

Vector2 DenseMatrixRow2(const mfem::DenseMatrix &matrix, int row)
{
  return {matrix(row, 0), matrix(row, 1)};
}

double
NormalizeTriangleBarycentricGradients(const TriangleBarycentricGradients &grad_lambda,
                                      TriangleBarycentricGradients &normalized_grad_lambda)
{
  double scale = 1.0;
  for (const auto &gradient : grad_lambda)
  {
    for (double value : gradient)
    {
      scale = std::max(scale, std::abs(value));
    }
  }
  if (!std::isfinite(scale))
  {
    throw std::invalid_argument(
        "Triangular singular affine normalization requires finite gradients!");
  }
  normalized_grad_lambda = grad_lambda;
  for (auto &gradient : normalized_grad_lambda)
  {
    for (double &value : gradient)
    {
      value /= scale;
    }
  }
  return scale;
}

double ScaleTriangleJacobianDeterminant(double jacobian_determinant, double gradient_scale,
                                        int gradient_power)
{
  long double scaled = jacobian_determinant;
  for (int i = 0; i < gradient_power; i++)
  {
    scaled *= gradient_scale;
  }
  const double result = static_cast<double>(scaled);
  if (!std::isfinite(result) || !(result > 0.0))
  {
    throw std::overflow_error(
        "Triangular singular affine geometry scaling produced an invalid prefactor!");
  }
  return result;
}

TriangleVectorBasisValue EvaluateAffineNormalizedTriangleBasis(
    const TriangleBarycentricPoint &lambda,
    const TriangleBarycentricGradients &normalized_grad_lambda,
    double normalized_jacobian_determinant, const TriangleBasis &basis)
{
  const auto reference =
      EvaluateTriangleBasis(lambda, ReferenceTriangleBarycentricGradients(), basis);
  TriangleVectorBasisValue result{{0.0, 0.0},
                                  reference.curl / normalized_jacobian_determinant};
  for (int d = 0; d < 2; d++)
  {
    result.value[d] = reference.value[0] * normalized_grad_lambda[1][d] +
                      reference.value[1] * normalized_grad_lambda[2][d];
  }
  return result;
}

void CalcAffinePhysVShape(const mfem::FiniteElement &finite_element,
                          const mfem::IntegrationPoint &point,
                          const mfem::DenseMatrix &inverse_jacobian,
                          mfem::DenseMatrix &reference_shape,
                          mfem::DenseMatrix &physical_shape)
{
  finite_element.CalcVShape(point, reference_shape);
  mfem::Mult(reference_shape, inverse_jacobian, physical_shape);
}

void CalcAffinePhysCurlShape(const mfem::FiniteElement &finite_element,
                             const mfem::IntegrationPoint &point,
                             const mfem::DenseMatrix &jacobian, double jacobian_determinant,
                             mfem::DenseMatrix &reference_curl,
                             mfem::DenseMatrix &physical_curl)
{
  finite_element.CalcCurlShape(point, reference_curl);
  if (finite_element.GetDim() == 2)
  {
    physical_curl = reference_curl;
  }
  else if (finite_element.GetDim() == 3)
  {
    mfem::MultABt(reference_curl, jacobian, physical_curl);
  }
  else
  {
    throw std::invalid_argument(
        "Singular affine H(curl) transformation requires dimension two or three!");
  }
  physical_curl *= 1.0 / jacobian_determinant;
}

void CalcAffinePhysDShape(const mfem::FiniteElement &finite_element,
                          const mfem::IntegrationPoint &point,
                          const mfem::DenseMatrix &inverse_jacobian,
                          mfem::DenseMatrix &reference_gradient,
                          mfem::DenseMatrix &physical_gradient)
{
  finite_element.CalcDShape(point, reference_gradient);
  mfem::Mult(reference_gradient, inverse_jacobian, physical_gradient);
}

void ValidateAdaptiveEntry(const AdaptiveQuadratureResult &entry, std::string_view quantity,
                           int standard, int enrichment, const HigherOrderBasis &basis)
{
  if (!entry.converged)
  {
    throw std::runtime_error(fmt::format(
        "Singular element adaptive {} coupling integral did not converge: standard ND "
        "DOF = {}, enrichment ND DOF = {}, family = {}, order = {}, nu = {:.17g}, "
        "nodes = [{}, {}, {}, {}], interpolation = [{}, {}, {}, {}], value = {:.17g}, "
        "estimated absolute error = {:.17g}, leaves = {}, maximum depth = {}!",
        quantity, standard, enrichment, static_cast<int>(basis.family), basis.order,
        basis.nu, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
        basis.interpolation_indices[0], basis.interpolation_indices[1],
        basis.interpolation_indices[2], basis.interpolation_indices[3], entry.value,
        entry.estimated_absolute_error, entry.leaf_count, entry.maximum_subdivision_depth));
  }
  if (!std::isfinite(entry.value) || !std::isfinite(entry.estimated_absolute_error) ||
      entry.estimated_absolute_error < 0.0)
  {
    throw std::runtime_error(fmt::format(
        "Singular element adaptive {} coupling produced an invalid entry: standard ND "
        "DOF = {}, enrichment ND DOF = {}, family = {}, order = {}, nu = {:.17g}, "
        "value = {:.17g}, estimated absolute error = {:.17g}!",
        quantity, standard, enrichment, static_cast<int>(basis.family), basis.order,
        basis.nu, entry.value, entry.estimated_absolute_error));
  }
}

void ValidateAdaptiveH1Entry(const AdaptiveQuadratureResult &entry,
                             std::string_view quantity, int row, int column,
                             const HigherOrderBasis &basis)
{
  if (!entry.converged)
  {
    throw std::runtime_error(fmt::format(
        "Singular element adaptive H1 {} integral did not converge: row DOF = {}, "
        "column DOF = {}, family = {}, order = {}, nu = {:.17g}, nodes = [{}, {}, {}, "
        "{}], interpolation = [{}, {}, {}, {}], value = {:.17g}, estimated absolute "
        "error = {:.17g}, leaves = {}, maximum depth = {}!",
        quantity, row, column, static_cast<int>(basis.family), basis.order, basis.nu,
        basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
        basis.interpolation_indices[0], basis.interpolation_indices[1],
        basis.interpolation_indices[2], basis.interpolation_indices[3], entry.value,
        entry.estimated_absolute_error, entry.leaf_count, entry.maximum_subdivision_depth));
  }
  if (!std::isfinite(entry.value) || !std::isfinite(entry.estimated_absolute_error) ||
      entry.estimated_absolute_error < 0.0)
  {
    throw std::runtime_error(fmt::format(
        "Singular element adaptive H1 {} integral produced an invalid entry: row DOF = "
        "{}, column DOF = {}, family = {}, order = {}, nu = {:.17g}, value = {:.17g}, "
        "estimated absolute error = {:.17g}!",
        quantity, row, column, static_cast<int>(basis.family), basis.order, basis.nu,
        entry.value, entry.estimated_absolute_error));
  }
}

void ValidateDuffyH1Entry(double value, double estimated_absolute_error,
                          std::string_view quantity, int row, int column,
                          const HigherOrderBasis &basis,
                          const AdaptiveAssemblyOptions &options)
{
  const double tolerance =
      options.absolute_tolerance + options.relative_tolerance * std::abs(value);
  if (!std::isfinite(value) || !std::isfinite(estimated_absolute_error) ||
      estimated_absolute_error < 0.0)
  {
    throw std::runtime_error(fmt::format(
        "Singular H1 Duffy {} integral produced an invalid entry: row DOF = {}, "
        "column DOF = {}, family = {}, order = {}, nu = {:.17g}, value = {:.17g}, "
        "estimated absolute error = {:.17g}!",
        quantity, row, column, static_cast<int>(basis.family), basis.order, basis.nu, value,
        estimated_absolute_error));
  }
  if (estimated_absolute_error > tolerance)
  {
    throw std::runtime_error(fmt::format(
        "Singular H1 Duffy {} integral did not converge: row DOF = {}, column DOF = "
        "{}, family = {}, order = {}, nu = {:.17g}, nodes = [{}, {}, {}, {}], "
        "interpolation = [{}, {}, {}, {}], value = {:.17g}, estimated absolute error "
        "= {:.17g}, tolerance = {:.17g}, reference orders = {}/{}!",
        quantity, row, column, static_cast<int>(basis.family), basis.order, basis.nu,
        basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
        basis.interpolation_indices[0], basis.interpolation_indices[1],
        basis.interpolation_indices[2], basis.interpolation_indices[3], value,
        estimated_absolute_error, tolerance, H1DuffyComparisonOrder,
        H1DuffyReferenceOrder));
  }
}

bool BasisLess(const HigherOrderBasis &a, const HigherOrderBasis &b)
{
  return std::tie(a.family, a.nodes, a.interpolation_indices, a.order, a.nu) <
         std::tie(b.family, b.nodes, b.interpolation_indices, b.order, b.nu);
}

struct AdaptiveCoefficientTensor
{
  ReferenceCoefficientTensor value{};
  ReferenceCoefficientTensor estimated_absolute_error{};
  std::size_t total_leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

struct AffineGeometryTensors
{
  ReferenceCoefficientTensor mass{};
  ReferenceCoefficientTensor curl_curl{};
};

AffineGeometryTensors BuildAffineGeometryTensors(const BarycentricGradients &grad_lambda)
{
  AffineGeometryTensors result;
  const std::array<Vector3, 3> curl_geometry{Cross(grad_lambda[2], grad_lambda[3]),
                                             Cross(grad_lambda[3], grad_lambda[1]),
                                             Cross(grad_lambda[1], grad_lambda[2])};
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      result.mass[u][v] = Dot(grad_lambda[u + 1], grad_lambda[v + 1]);
      result.curl_curl[u][v] = Dot(curl_geometry[u], curl_geometry[v]);
    }
  }
  return result;
}

double ContractAffineTensor(const ReferenceCoefficientTensor &tensor,
                            const ReferenceCoefficientTensor &geometry,
                            double jacobian_determinant, bool error)
{
  long double result = 0.0L;
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      result += static_cast<long double>(jacobian_determinant) *
                static_cast<long double>(tensor[u][v]) *
                (error ? std::abs(static_cast<long double>(geometry[u][v]))
                       : static_cast<long double>(geometry[u][v]));
    }
  }
  return static_cast<double>(result);
}

class AffineStandardReferenceTable
{
public:
  struct Entry
  {
    AdaptiveCoefficientTensor mass;
    AdaptiveCoefficientTensor curl_curl;
  };

private:
  struct Key
  {
    const mfem::FiniteElement *finite_element;
    HigherOrderBasis basis;
  };

  struct KeyLess
  {
    bool operator()(const Key &a, const Key &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return BasisLess(a.basis, b.basis);
    }
  };

  struct PatternKey
  {
    const mfem::FiniteElement *finite_element;
    std::vector<HigherOrderBasis> bases;
  };

  struct PatternKeyLess
  {
    bool operator()(const PatternKey &a, const PatternKey &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return std::lexicographical_compare(a.bases.begin(), a.bases.end(), b.bases.begin(),
                                          b.bases.end(), BasisLess);
    }
  };

  AdaptiveAssemblyOptions options;
  std::map<Key, std::vector<Entry>, KeyLess> entries;
  std::map<PatternKey, std::vector<const std::vector<Entry> *>, PatternKeyLess> patterns;
  std::size_t hits = 0;
  double generation_time = 0.0;

  const std::vector<Entry> &GetEntries(const mfem::FiniteElement &finite_element,
                                       const HigherOrderBasis &basis)
  {
    const int standard_size = finite_element.GetDof();
    const Key key{&finite_element, basis};
    const auto cached = entries.find(key);
    if (cached != entries.end())
    {
      hits += static_cast<std::size_t>(standard_size);
      return cached->second;
    }

    const auto generation_start = std::chrono::steady_clock::now();
    mfem::IntegrationPoint point;
    mfem::DenseMatrix standard_value(standard_size, 3);
    mfem::DenseMatrix standard_curl(standard_size, 3);
    const bool has_curl = !IsGradientFamily(basis.family);
    constexpr int tensor_size = 9;
    const std::size_t mass_size = static_cast<std::size_t>(standard_size) * tensor_size;
    const std::size_t number_components = has_curl ? 2 * mass_size : mass_size;
    const auto component = [standard_size](int standard_dof, int u, int v, bool curl)
    {
      constexpr int tensor_size = 9;
      return static_cast<std::size_t>(curl ? standard_size * tensor_size : 0) +
             static_cast<std::size_t>(standard_dof * tensor_size + 3 * u + v);
    };
    const auto integral = IntegrateConfigured(
        options.quadrature_order, options, number_components,
        [&](const BarycentricPoint &lambda, std::vector<double> &value)
        {
          point.Set3(lambda[1], lambda[2], lambda[3]);
          finite_element.CalcVShape(point, standard_value);
          if (has_curl)
          {
            finite_element.CalcCurlShape(point, standard_curl);
          }
          const auto singular =
              EvaluateHigherOrderBasis(lambda, ReferenceBarycentricGradients(), basis);
          for (int standard_dof = 0; standard_dof < standard_size; standard_dof++)
          {
            for (int u = 0; u < 3; u++)
            {
              for (int v = 0; v < 3; v++)
              {
                value[component(standard_dof, u, v, false)] =
                    standard_value(standard_dof, u) * singular.value[v];
                if (has_curl)
                {
                  value[component(standard_dof, u, v, true)] =
                      standard_curl(standard_dof, u) * singular.curl[v];
                }
              }
            }
          }
        });
    if (!integral.converged)
    {
      std::size_t failed_component = 0;
      for (; failed_component < number_components; failed_component++)
      {
        const double tolerance =
            options.absolute_tolerance +
            options.relative_tolerance * std::abs(integral.value[failed_component]);
        if (integral.estimated_absolute_error[failed_component] > tolerance)
        {
          break;
        }
      }
      const bool failed_curl = failed_component >= mass_size;
      const std::size_t tensor_component = failed_component - (failed_curl ? mass_size : 0);
      const int failed_standard = static_cast<int>(tensor_component / tensor_size);
      const int failed_u = static_cast<int>((tensor_component % tensor_size) / 3);
      const int failed_v = static_cast<int>(tensor_component % 3);
      throw std::runtime_error(fmt::format(
          "Cached affine singular standard-enrichment tensors did not converge: "
          "standard ND DOF = {}, quantity = {}, component = [{}, {}], "
          "family = {}, order = {}, nu = {:.17g}, nodes = [{}, {}, {}, {}], "
          "interpolation = [{}, {}, {}, {}], components = {}, leaves = {}, maximum "
          "depth = {}!",
          failed_standard, failed_curl ? "curl-curl" : "mass", failed_u, failed_v,
          static_cast<int>(basis.family), basis.order, basis.nu, basis.nodes[0],
          basis.nodes[1], basis.nodes[2], basis.nodes[3], basis.interpolation_indices[0],
          basis.interpolation_indices[1], basis.interpolation_indices[2],
          basis.interpolation_indices[3], number_components, integral.leaf_count,
          integral.maximum_subdivision_depth));
    }

    std::vector<Entry> basis_entries(standard_size);
    for (int standard_dof = 0; standard_dof < standard_size; standard_dof++)
    {
      auto &entry = basis_entries[standard_dof];
      entry.mass.total_leaf_count = integral.leaf_count;
      entry.mass.maximum_subdivision_depth = integral.maximum_subdivision_depth;
      entry.curl_curl.total_leaf_count = has_curl ? integral.leaf_count : 0;
      entry.curl_curl.maximum_subdivision_depth =
          has_curl ? integral.maximum_subdivision_depth : 0;
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          const auto mass_index = component(standard_dof, u, v, false);
          entry.mass.value[u][v] = integral.value[mass_index];
          entry.mass.estimated_absolute_error[u][v] =
              integral.estimated_absolute_error[mass_index];
          if (has_curl)
          {
            const auto curl_index = component(standard_dof, u, v, true);
            entry.curl_curl.value[u][v] = integral.value[curl_index];
            entry.curl_curl.estimated_absolute_error[u][v] =
                integral.estimated_absolute_error[curl_index];
          }
        }
      }
    }
    hits += static_cast<std::size_t>(standard_size - 1);
    const auto &result = entries.emplace(key, std::move(basis_entries)).first->second;
    generation_time +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
            .count();
    return result;
  }

public:
  explicit AffineStandardReferenceTable(const AdaptiveAssemblyOptions &options_in)
    : options(options_in)
  {
  }

  const std::vector<const std::vector<Entry> *> &
  GetPattern(const mfem::FiniteElement &finite_element,
             const std::vector<ElementDof> &element_dofs)
  {
    const int standard_size = finite_element.GetDof();
    PatternKey key{&finite_element, {}};
    key.bases.reserve(element_dofs.size());
    for (const auto &element_dof : element_dofs)
    {
      key.bases.push_back(element_dof.basis);
    }
    const auto cached = patterns.find(key);
    if (cached != patterns.end())
    {
      hits += static_cast<std::size_t>(standard_size) * element_dofs.size();
      return cached->second;
    }

    std::vector<const std::vector<Entry> *> pattern(element_dofs.size());
    for (std::size_t enrichment = 0; enrichment < element_dofs.size(); enrichment++)
    {
      const auto &basis_entries =
          GetEntries(finite_element, element_dofs[enrichment].basis);
      pattern[enrichment] = &basis_entries;
    }
    return patterns.emplace(std::move(key), std::move(pattern)).first->second;
  }

  std::size_t Size() const
  {
    std::size_t size = 0;
    for (const auto &[key, basis_entries] : entries)
    {
      (void)key;
      size += basis_entries.size();
    }
    return size;
  }
  std::size_t Hits() const { return hits; }
  double GenerationTime() const { return generation_time; }
};

template <typename T>
void AppendReferenceCacheKey(std::vector<std::uint8_t> &key, const T &value)
{
  static_assert(std::is_trivially_copyable_v<T>);
  const auto *data = reinterpret_cast<const std::uint8_t *>(&value);
  key.insert(key.end(), data, data + sizeof(T));
}

void UpdateReferenceCacheHash(std::uint64_t &hash, const void *data, std::size_t size)
{
  const auto *bytes = static_cast<const std::uint8_t *>(data);
  for (std::size_t i = 0; i < size; i++)
  {
    hash ^= bytes[i];
    hash *= 1099511628211ULL;
  }
}

std::uint64_t HashReferenceCacheKey(const std::vector<std::uint8_t> &key)
{
  std::uint64_t hash = 1469598103934665603ULL;
  UpdateReferenceCacheHash(hash, key.data(), key.size());
  return hash;
}

template <typename T>
bool ReadReferenceCacheValue(std::istream &stream, T &value)
{
  static_assert(std::is_trivially_copyable_v<T>);
  stream.read(reinterpret_cast<char *>(&value), sizeof(T));
  return static_cast<bool>(stream);
}

template <typename T>
void WriteReferenceCacheValue(std::ostream &stream, const T &value)
{
  static_assert(std::is_trivially_copyable_v<T>);
  stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
  if (!stream)
  {
    throw std::runtime_error("Failed to write a singular reference cache value!");
  }
}

struct AffineCommonNDMassBlocks
{
  mfem::DenseMatrix standard_enrichment;
  mfem::DenseMatrix enrichment_standard;
  mfem::DenseMatrix enrichment_enrichment;
  mfem::DenseMatrix standard_enrichment_estimated_absolute_error;
  mfem::DenseMatrix enrichment_enrichment_estimated_absolute_error;
  std::size_t leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

// All standard-enrichment and enrichment-enrichment mass entries in one
// element use the same positive adaptive quadrature rule. For affine
// tetrahedra, a base rule exact for the standard-standard polynomial products
// makes the completed matrix a Gram matrix up to floating-point roundoff.
class AffineCommonNDMassReferenceTable
{
private:
  struct Key
  {
    const mfem::FiniteElement *finite_element;
    std::vector<HigherOrderBasis> bases;
  };

  struct KeyLess
  {
    bool operator()(const Key &a, const Key &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return std::lexicographical_compare(a.bases.begin(), a.bases.end(), b.bases.begin(),
                                          b.bases.end(), BasisLess);
    }
  };

  struct Entry
  {
    int standard_size = 0;
    int enrichment_size = 0;
    std::vector<AdaptiveCoefficientTensor> standard_enrichment;
    std::vector<AdaptiveCoefficientTensor> enrichment_enrichment;
    std::size_t leaf_count = 0;
    int maximum_subdivision_depth = 0;
  };

  AdaptiveAssemblyOptions options;
  std::map<Key, Entry, KeyLess> entries;
  std::size_t hits = 0;
  std::size_t pattern_hits = 0;
  std::size_t persistent_hits = 0;
  std::size_t persistent_writes = 0;
  std::size_t generated_leaf_count = 0;
  double generation_time = 0.0;

  static std::size_t EnrichmentPairIndex(int row, int column, int size)
  {
    if (row > column)
    {
      std::swap(row, column);
    }
    return static_cast<std::size_t>(row * size - row * (row - 1) / 2 + column - row);
  }

  const Entry &GetEntry(const mfem::FiniteElement &finite_element,
                        const std::vector<ElementDof> &element_dofs)
  {
    Key key{&finite_element, {}};
    key.bases.reserve(element_dofs.size());
    for (const auto &element_dof : element_dofs)
    {
      key.bases.push_back(element_dof.basis);
    }
    const auto cached = entries.find(key);
    if (cached != entries.end())
    {
      hits += cached->second.standard_enrichment.size() +
              cached->second.enrichment_enrichment.size();
      pattern_hits++;
      return cached->second;
    }

    const auto persistent_key = BuildPersistentKey(finite_element, element_dofs);
    const auto persistent_path = GetPersistentPath(persistent_key);
    const auto generation_start = std::chrono::steady_clock::now();
    Entry entry;
    entry.standard_size = finite_element.GetDof();
    entry.enrichment_size = static_cast<int>(element_dofs.size());
    const std::size_t standard_enrichment_size =
        static_cast<std::size_t>(entry.standard_size) * entry.enrichment_size;
    const std::size_t enrichment_enrichment_size =
        static_cast<std::size_t>(entry.enrichment_size) * (entry.enrichment_size + 1) / 2;
    constexpr std::size_t tensor_size = 9;
    const std::size_t number_tensors =
        standard_enrichment_size + enrichment_enrichment_size;
    if (number_tensors == 0 ||
        number_tensors > std::numeric_limits<std::size_t>::max() / tensor_size)
    {
      throw std::overflow_error(
          "Common affine ND mass quadrature has invalid matrix dimensions!");
    }
    const bool persistent_file_existed =
        !persistent_path.empty() && std::filesystem::is_regular_file(persistent_path);
    if (LoadPersistentEntry(persistent_path, persistent_key, standard_enrichment_size,
                            enrichment_enrichment_size, entry))
    {
      return entries.emplace(std::move(key), std::move(entry)).first->second;
    }
    if (persistent_file_existed)
    {
      std::error_code error;
      std::filesystem::remove(persistent_path, error);
    }

    const int quadrature_order =
        std::max(options.quadrature_order, 2 * finite_element.GetOrder());
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, quadrature_order);
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      if (!std::isfinite(rule.IntPoint(q).weight) || !(rule.IntPoint(q).weight > 0.0))
      {
        throw std::runtime_error(
            "Common affine ND mass quadrature requires a positive base rule!");
      }
    }

    mfem::IntegrationPoint point;
    mfem::DenseMatrix standard_value(entry.standard_size, 3);
    std::vector<Vector3> enrichment_value(static_cast<std::size_t>(entry.enrichment_size));
    const auto component = [](std::size_t tensor, int u, int v)
    { return tensor_size * tensor + static_cast<std::size_t>(3 * u + v); };
    const auto integral = IntegrateConfigured(
        quadrature_order, options, tensor_size * number_tensors,
        [&](const BarycentricPoint &lambda, std::vector<double> &value)
        {
          point.Set3(lambda[1], lambda[2], lambda[3]);
          finite_element.CalcVShape(point, standard_value);
          for (int enrichment = 0; enrichment < entry.enrichment_size; enrichment++)
          {
            enrichment_value[enrichment] =
                EvaluateHigherOrderBasis(lambda, ReferenceBarycentricGradients(),
                                         element_dofs[enrichment].basis)
                    .value;
          }

          std::size_t tensor = 0;
          for (int standard = 0; standard < entry.standard_size; standard++)
          {
            for (int enrichment = 0; enrichment < entry.enrichment_size; enrichment++)
            {
              for (int u = 0; u < 3; u++)
              {
                for (int v = 0; v < 3; v++)
                {
                  value[component(tensor, u, v)] =
                      standard_value(standard, u) * enrichment_value[enrichment][v];
                }
              }
              tensor++;
            }
          }
          for (int row = 0; row < entry.enrichment_size; row++)
          {
            for (int column = row; column < entry.enrichment_size; column++)
            {
              for (int u = 0; u < 3; u++)
              {
                for (int v = 0; v < 3; v++)
                {
                  value[component(tensor, u, v)] =
                      enrichment_value[row][u] * enrichment_value[column][v];
                }
              }
              tensor++;
            }
          }
          if (tensor != number_tensors)
          {
            throw std::logic_error(
                "Common affine ND mass quadrature filled inconsistent dimensions!");
          }
        });
    if (!integral.converged)
    {
      throw std::runtime_error(fmt::format(
          "Common affine ND mass quadrature did not converge: components = {}, "
          "leaves = {}, maximum depth = {}!",
          integral.value.size(), integral.leaf_count, integral.maximum_subdivision_depth));
    }

    entry.standard_enrichment.resize(standard_enrichment_size);
    entry.enrichment_enrichment.resize(enrichment_enrichment_size);
    for (std::size_t tensor = 0; tensor < number_tensors; tensor++)
    {
      auto &destination =
          tensor < standard_enrichment_size
              ? entry.standard_enrichment[tensor]
              : entry.enrichment_enrichment[tensor - standard_enrichment_size];
      destination.total_leaf_count = integral.leaf_count;
      destination.maximum_subdivision_depth = integral.maximum_subdivision_depth;
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          const std::size_t index = component(tensor, u, v);
          destination.value[u][v] = integral.value[index];
          destination.estimated_absolute_error[u][v] =
              integral.estimated_absolute_error[index];
        }
      }
    }
    entry.leaf_count = integral.leaf_count;
    entry.maximum_subdivision_depth = integral.maximum_subdivision_depth;
    if (generated_leaf_count >
        std::numeric_limits<std::size_t>::max() - integral.leaf_count)
    {
      throw std::overflow_error("Common affine ND mass generated leaf count overflow!");
    }
    generated_leaf_count += integral.leaf_count;
    if (SavePersistentEntry(persistent_path, persistent_key, entry))
    {
      persistent_writes++;
    }
    hits += number_tensors - 1;
    const auto &result = entries.emplace(std::move(key), std::move(entry)).first->second;
    generation_time +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
            .count();
    return result;
  }

  std::vector<std::uint8_t>
  BuildPersistentKey(const mfem::FiniteElement &finite_element,
                     const std::vector<ElementDof> &element_dofs) const
  {
    std::vector<std::uint8_t> key;
    constexpr std::uint32_t format_version = 2;
    AppendReferenceCacheKey(key, format_version);
    AppendReferenceCacheKey(key, static_cast<std::uint32_t>(MFEM_VERSION));
    AppendReferenceCacheKey(
        key, static_cast<std::uint32_t>(ReferenceIntegral::ConventionVersion));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(finite_element.GetGeomType()));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(finite_element.GetDim()));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(finite_element.GetOrder()));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(finite_element.GetDof()));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(options.quadrature_order));
    AppendReferenceCacheKey(key, static_cast<std::int32_t>(options.subdivisions));
    AppendReferenceCacheKey(key, static_cast<std::uint64_t>(element_dofs.size()));
    for (const auto &element_dof : element_dofs)
    {
      const auto &basis = element_dof.basis;
      AppendReferenceCacheKey(key, static_cast<std::int32_t>(basis.family));
      for (int node : basis.nodes)
      {
        AppendReferenceCacheKey(key, static_cast<std::int32_t>(node));
      }
      for (int index : basis.interpolation_indices)
      {
        AppendReferenceCacheKey(key, static_cast<std::int32_t>(index));
      }
      AppendReferenceCacheKey(key, static_cast<std::int32_t>(basis.order));
      std::uint64_t exponent;
      static_assert(sizeof(exponent) == sizeof(basis.nu));
      std::memcpy(&exponent, &basis.nu, sizeof(exponent));
      AppendReferenceCacheKey(key, exponent);
    }
    return key;
  }

  std::uint64_t ComputePersistentChecksum(const std::vector<std::uint8_t> &key,
                                          const Entry &entry) const
  {
    std::uint64_t checksum = 1469598103934665603ULL;
    UpdateReferenceCacheHash(checksum, key.data(), key.size());
    const std::uint64_t standard_size = entry.standard_enrichment.size();
    const std::uint64_t enrichment_size = entry.enrichment_enrichment.size();
    const std::uint64_t leaf_count = entry.leaf_count;
    const std::int32_t maximum_depth = entry.maximum_subdivision_depth;
    UpdateReferenceCacheHash(checksum, &standard_size, sizeof(standard_size));
    UpdateReferenceCacheHash(checksum, &enrichment_size, sizeof(enrichment_size));
    UpdateReferenceCacheHash(checksum, &leaf_count, sizeof(leaf_count));
    UpdateReferenceCacheHash(checksum, &maximum_depth, sizeof(maximum_depth));
    for (const auto &tensor : entry.standard_enrichment)
    {
      UpdateReferenceCacheHash(checksum, &tensor.value, sizeof(tensor.value));
    }
    for (const auto &tensor : entry.enrichment_enrichment)
    {
      UpdateReferenceCacheHash(checksum, &tensor.value, sizeof(tensor.value));
    }
    return checksum;
  }

  std::filesystem::path GetPersistentPath(const std::vector<std::uint8_t> &key) const
  {
    if (!options.fixed_subdivision || options.reference_cache.empty())
    {
      return {};
    }
    std::ostringstream name;
    name << "affine-h1-mass-v2-" << std::hex << std::setw(16) << std::setfill('0')
         << HashReferenceCacheKey(key) << ".bin";
    return std::filesystem::path(options.reference_cache) / name.str();
  }

  bool LoadPersistentEntry(const std::filesystem::path &path,
                           const std::vector<std::uint8_t> &key,
                           std::size_t standard_enrichment_size,
                           std::size_t enrichment_enrichment_size, Entry &entry)
  {
    if (path.empty() || !std::filesystem::is_regular_file(path))
    {
      return false;
    }
    std::ifstream stream(path, std::ios::binary);
    constexpr std::uint64_t magic = 0x50414c4852454631ULL;
    constexpr std::uint32_t format_version = 2;
    std::uint64_t file_magic, checksum, key_size, standard_size, enrichment_size,
        leaf_count;
    std::uint32_t file_version;
    std::int32_t maximum_depth;
    if (!ReadReferenceCacheValue(stream, file_magic) ||
        !ReadReferenceCacheValue(stream, file_version) ||
        !ReadReferenceCacheValue(stream, checksum) ||
        !ReadReferenceCacheValue(stream, key_size) || file_magic != magic ||
        file_version != format_version || key_size != key.size() ||
        !ReadReferenceCacheValue(stream, standard_size) ||
        !ReadReferenceCacheValue(stream, enrichment_size) ||
        !ReadReferenceCacheValue(stream, leaf_count) ||
        !ReadReferenceCacheValue(stream, maximum_depth) ||
        standard_size != standard_enrichment_size ||
        enrichment_size != enrichment_enrichment_size ||
        leaf_count > std::numeric_limits<std::size_t>::max() || maximum_depth < 0)
    {
      return false;
    }
    std::vector<std::uint8_t> file_key(key.size());
    stream.read(reinterpret_cast<char *>(file_key.data()),
                static_cast<std::streamsize>(file_key.size()));
    if (!stream || file_key != key)
    {
      return false;
    }

    entry.standard_enrichment.resize(standard_enrichment_size);
    entry.enrichment_enrichment.resize(enrichment_enrichment_size);
    const auto read_tensors = [&](std::vector<AdaptiveCoefficientTensor> &tensors)
    {
      for (auto &tensor : tensors)
      {
        if (!ReadReferenceCacheValue(stream, tensor.value))
        {
          return false;
        }
        for (const auto &row : tensor.value)
        {
          if (!std::all_of(row.begin(), row.end(),
                           [](double value) { return std::isfinite(value); }))
          {
            return false;
          }
        }
        tensor.total_leaf_count = static_cast<std::size_t>(leaf_count);
        tensor.maximum_subdivision_depth = maximum_depth;
      }
      return true;
    };
    if (!read_tensors(entry.standard_enrichment) ||
        !read_tensors(entry.enrichment_enrichment))
    {
      return false;
    }
    entry.leaf_count = static_cast<std::size_t>(leaf_count);
    entry.maximum_subdivision_depth = maximum_depth;
    if (ComputePersistentChecksum(key, entry) != checksum ||
        stream.peek() != std::char_traits<char>::eof())
    {
      return false;
    }
    persistent_hits++;
    return true;
  }

  bool SavePersistentEntry(const std::filesystem::path &path,
                           const std::vector<std::uint8_t> &key, const Entry &entry) const
  {
    if (path.empty())
    {
      return false;
    }
    std::error_code error;
    std::filesystem::create_directories(path.parent_path(), error);
    if (error)
    {
      throw std::runtime_error(
          fmt::format("Failed to create singular reference cache directory '{}': {}",
                      path.parent_path().string(), error.message()));
    }
    if (std::filesystem::is_regular_file(path))
    {
      return false;
    }
    std::string temporary_template = path.string() + ".tmp.XXXXXX";
    std::vector<char> temporary_name(temporary_template.begin(), temporary_template.end());
    temporary_name.push_back('\0');
    const int temporary_descriptor = mkstemp(temporary_name.data());
    if (temporary_descriptor < 0)
    {
      throw std::runtime_error(
          fmt::format("Failed to create a temporary singular reference cache file in "
                      "'{}'!",
                      path.parent_path().string()));
    }
    close(temporary_descriptor);
    const std::filesystem::path temporary(temporary_name.data());
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream)
    {
      std::filesystem::remove(temporary);
      throw std::runtime_error(
          fmt::format("Failed to open singular reference cache file '{}' for writing!",
                      temporary.string()));
    }
    constexpr std::uint64_t magic = 0x50414c4852454631ULL;
    constexpr std::uint32_t format_version = 2;
    WriteReferenceCacheValue(stream, magic);
    WriteReferenceCacheValue(stream, format_version);
    WriteReferenceCacheValue(stream, ComputePersistentChecksum(key, entry));
    WriteReferenceCacheValue(stream, static_cast<std::uint64_t>(key.size()));
    WriteReferenceCacheValue(stream,
                             static_cast<std::uint64_t>(entry.standard_enrichment.size()));
    WriteReferenceCacheValue(
        stream, static_cast<std::uint64_t>(entry.enrichment_enrichment.size()));
    WriteReferenceCacheValue(stream, static_cast<std::uint64_t>(entry.leaf_count));
    WriteReferenceCacheValue(stream,
                             static_cast<std::int32_t>(entry.maximum_subdivision_depth));
    stream.write(reinterpret_cast<const char *>(key.data()),
                 static_cast<std::streamsize>(key.size()));
    for (const auto &tensor : entry.standard_enrichment)
    {
      WriteReferenceCacheValue(stream, tensor.value);
    }
    for (const auto &tensor : entry.enrichment_enrichment)
    {
      WriteReferenceCacheValue(stream, tensor.value);
    }
    stream.close();
    if (!stream)
    {
      std::filesystem::remove(temporary);
      throw std::runtime_error(fmt::format(
          "Failed to finalize singular reference cache file '{}'!", temporary.string()));
    }
    std::filesystem::create_hard_link(temporary, path, error);
    std::filesystem::remove(temporary);
    if (error)
    {
      if (error == std::errc::file_exists && std::filesystem::is_regular_file(path))
      {
        return false;
      }
      throw std::runtime_error(
          fmt::format("Failed to install singular reference cache file '{}': {}",
                      path.string(), error.message()));
    }
    return true;
  }

public:
  explicit AffineCommonNDMassReferenceTable(const AdaptiveAssemblyOptions &options_in)
    : options(options_in)
  {
  }

  AffineCommonNDMassBlocks Get(const mfem::FiniteElement &finite_element,
                               const std::vector<ElementDof> &element_dofs,
                               const BarycentricGradients &grad_lambda,
                               double jacobian_determinant)
  {
    const auto &entry = GetEntry(finite_element, element_dofs);
    const auto geometry = BuildAffineGeometryTensors(grad_lambda);
    AffineCommonNDMassBlocks result;
    result.standard_enrichment.SetSize(entry.standard_size, entry.enrichment_size);
    result.enrichment_standard.SetSize(entry.enrichment_size, entry.standard_size);
    result.enrichment_enrichment.SetSize(entry.enrichment_size);
    result.standard_enrichment_estimated_absolute_error.SetSize(entry.standard_size,
                                                                entry.enrichment_size);
    result.enrichment_enrichment_estimated_absolute_error.SetSize(entry.enrichment_size);
    for (int standard = 0; standard < entry.standard_size; standard++)
    {
      for (int enrichment = 0; enrichment < entry.enrichment_size; enrichment++)
      {
        const auto &reference =
            entry.standard_enrichment[static_cast<std::size_t>(standard) *
                                          entry.enrichment_size +
                                      enrichment];
        const double value = ContractAffineTensor(reference.value, geometry.mass,
                                                  jacobian_determinant, false);
        const double error = ContractAffineTensor(
            reference.estimated_absolute_error, geometry.mass, jacobian_determinant, true);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND mass contraction produced an invalid coupling entry!");
        }
        result.standard_enrichment(standard, enrichment) = value;
        result.enrichment_standard(enrichment, standard) = value;
        result.standard_enrichment_estimated_absolute_error(standard, enrichment) = error;
      }
    }
    for (int row = 0; row < entry.enrichment_size; row++)
    {
      for (int column = row; column < entry.enrichment_size; column++)
      {
        const auto &reference = entry.enrichment_enrichment[EnrichmentPairIndex(
            row, column, entry.enrichment_size)];
        const double value = ContractAffineTensor(reference.value, geometry.mass,
                                                  jacobian_determinant, false);
        const double error = ContractAffineTensor(
            reference.estimated_absolute_error, geometry.mass, jacobian_determinant, true);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND mass contraction produced an invalid enrichment entry!");
        }
        result.enrichment_enrichment(row, column) =
            result.enrichment_enrichment(column, row) = value;
        result.enrichment_enrichment_estimated_absolute_error(row, column) =
            result.enrichment_enrichment_estimated_absolute_error(column, row) = error;
      }
    }
    result.leaf_count = entry.leaf_count;
    result.maximum_subdivision_depth = entry.maximum_subdivision_depth;
    return result;
  }

  std::size_t Size() const
  {
    std::size_t size = 0;
    for (const auto &[key, entry] : entries)
    {
      (void)key;
      size += entry.standard_enrichment.size() + entry.enrichment_enrichment.size();
    }
    return size;
  }
  std::size_t Hits() const { return hits; }
  std::size_t PatternCount() const { return entries.size(); }
  std::size_t PatternHits() const { return pattern_hits; }
  std::size_t PersistentHits() const { return persistent_hits; }
  std::size_t PersistentWrites() const { return persistent_writes; }
  std::size_t GeneratedLeafCount() const { return generated_leaf_count; }
  double GenerationTime() const { return generation_time; }
};

struct AffineCommonNDCurlBlocks
{
  mfem::DenseMatrix standard_enrichment;
  mfem::DenseMatrix enrichment_standard;
  mfem::DenseMatrix enrichment_enrichment;
  mfem::DenseMatrix standard_enrichment_estimated_absolute_error;
  mfem::DenseMatrix enrichment_enrichment_estimated_absolute_error;
  std::size_t leaf_count = 0;
  int maximum_subdivision_depth = 0;
};

// The completed curl-curl matrix is another Gram matrix and must use one
// positive quadrature rule just like the ND mass matrix. Pairwise adaptive
// integration can make nearly dependent standard and singular curls
// materially indefinite even when every entry meets its scalar tolerance.
class AffineCommonNDCurlReferenceTable
{
private:
  struct Key
  {
    const mfem::FiniteElement *finite_element;
    std::vector<HigherOrderBasis> bases;
  };

  struct KeyLess
  {
    bool operator()(const Key &a, const Key &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return std::lexicographical_compare(a.bases.begin(), a.bases.end(), b.bases.begin(),
                                          b.bases.end(), BasisLess);
    }
  };

  struct Entry
  {
    int standard_size = 0;
    int enrichment_size = 0;
    std::vector<int> rotational_indices;
    std::vector<AdaptiveCoefficientTensor> standard_enrichment;
    std::vector<AdaptiveCoefficientTensor> enrichment_enrichment;
    std::size_t leaf_count = 0;
    int maximum_subdivision_depth = 0;
  };

  AdaptiveAssemblyOptions options;
  std::map<Key, Entry, KeyLess> entries;
  std::size_t hits = 0;
  double generation_time = 0.0;

  static std::size_t EnrichmentPairIndex(int row, int column, int size)
  {
    if (row > column)
    {
      std::swap(row, column);
    }
    return static_cast<std::size_t>(row * size - row * (row - 1) / 2 + column - row);
  }

  const Entry &GetEntry(const mfem::FiniteElement &finite_element,
                        const std::vector<ElementDof> &element_dofs)
  {
    Key key{&finite_element, {}};
    key.bases.reserve(element_dofs.size());
    for (const auto &element_dof : element_dofs)
    {
      key.bases.push_back(element_dof.basis);
    }
    const auto cached = entries.find(key);
    if (cached != entries.end())
    {
      hits += cached->second.standard_enrichment.size() +
              cached->second.enrichment_enrichment.size();
      return cached->second;
    }

    const auto generation_start = std::chrono::steady_clock::now();
    Entry entry;
    entry.standard_size = finite_element.GetDof();
    entry.enrichment_size = static_cast<int>(element_dofs.size());
    for (int enrichment = 0; enrichment < entry.enrichment_size; enrichment++)
    {
      if (!IsGradientFamily(element_dofs[enrichment].basis.family))
      {
        entry.rotational_indices.push_back(enrichment);
      }
    }
    const int rotational_size = static_cast<int>(entry.rotational_indices.size());
    const std::size_t standard_enrichment_size =
        static_cast<std::size_t>(entry.standard_size) * rotational_size;
    const std::size_t enrichment_enrichment_size =
        static_cast<std::size_t>(rotational_size) * (rotational_size + 1) / 2;
    const std::size_t number_tensors =
        standard_enrichment_size + enrichment_enrichment_size;
    if (number_tensors == 0)
    {
      const auto &result = entries.emplace(std::move(key), std::move(entry)).first->second;
      generation_time +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
              .count();
      return result;
    }

    constexpr std::size_t tensor_size = 9;
    if (number_tensors > std::numeric_limits<std::size_t>::max() / tensor_size)
    {
      throw std::overflow_error(
          "Common affine ND curl quadrature has invalid matrix dimensions!");
    }
    const int quadrature_order =
        std::max(options.quadrature_order, 2 * finite_element.GetOrder());
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, quadrature_order);
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      if (!std::isfinite(rule.IntPoint(q).weight) || !(rule.IntPoint(q).weight > 0.0))
      {
        throw std::runtime_error(
            "Common affine ND curl quadrature requires a positive base rule!");
      }
    }

    mfem::IntegrationPoint point;
    mfem::DenseMatrix standard_curl(entry.standard_size, 3);
    std::vector<Vector3> enrichment_curl(static_cast<std::size_t>(rotational_size));
    const auto component = [](std::size_t tensor, int u, int v)
    { return tensor_size * tensor + static_cast<std::size_t>(3 * u + v); };
    const auto integral = IntegrateConfigured(
        quadrature_order, options, tensor_size * number_tensors,
        [&](const BarycentricPoint &lambda, std::vector<double> &value)
        {
          point.Set3(lambda[1], lambda[2], lambda[3]);
          finite_element.CalcCurlShape(point, standard_curl);
          for (int rotational = 0; rotational < rotational_size; rotational++)
          {
            enrichment_curl[rotational] =
                EvaluateHigherOrderBasis(
                    lambda, ReferenceBarycentricGradients(),
                    element_dofs[entry.rotational_indices[rotational]].basis)
                    .curl;
          }

          std::size_t tensor = 0;
          for (int standard = 0; standard < entry.standard_size; standard++)
          {
            for (int rotational = 0; rotational < rotational_size; rotational++)
            {
              for (int u = 0; u < 3; u++)
              {
                for (int v = 0; v < 3; v++)
                {
                  value[component(tensor, u, v)] =
                      standard_curl(standard, u) * enrichment_curl[rotational][v];
                }
              }
              tensor++;
            }
          }
          for (int row = 0; row < rotational_size; row++)
          {
            for (int column = row; column < rotational_size; column++)
            {
              for (int u = 0; u < 3; u++)
              {
                for (int v = 0; v < 3; v++)
                {
                  value[component(tensor, u, v)] =
                      enrichment_curl[row][u] * enrichment_curl[column][v];
                }
              }
              tensor++;
            }
          }
          if (tensor != number_tensors)
          {
            throw std::logic_error(
                "Common affine ND curl quadrature filled inconsistent dimensions!");
          }
        });
    if (!integral.converged)
    {
      throw std::runtime_error(fmt::format(
          "Common affine ND curl quadrature did not converge: components = {}, "
          "leaves = {}, maximum depth = {}!",
          integral.value.size(), integral.leaf_count, integral.maximum_subdivision_depth));
    }

    entry.standard_enrichment.resize(standard_enrichment_size);
    entry.enrichment_enrichment.resize(enrichment_enrichment_size);
    for (std::size_t tensor = 0; tensor < number_tensors; tensor++)
    {
      auto &destination =
          tensor < standard_enrichment_size
              ? entry.standard_enrichment[tensor]
              : entry.enrichment_enrichment[tensor - standard_enrichment_size];
      destination.total_leaf_count = integral.leaf_count;
      destination.maximum_subdivision_depth = integral.maximum_subdivision_depth;
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          const std::size_t index = component(tensor, u, v);
          destination.value[u][v] = integral.value[index];
          destination.estimated_absolute_error[u][v] =
              integral.estimated_absolute_error[index];
        }
      }
    }
    entry.leaf_count = integral.leaf_count;
    entry.maximum_subdivision_depth = integral.maximum_subdivision_depth;
    hits += number_tensors - 1;
    const auto &result = entries.emplace(std::move(key), std::move(entry)).first->second;
    generation_time +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
            .count();
    return result;
  }

public:
  explicit AffineCommonNDCurlReferenceTable(const AdaptiveAssemblyOptions &options_in)
    : options(options_in)
  {
  }

  AffineCommonNDCurlBlocks Get(const mfem::FiniteElement &finite_element,
                               const std::vector<ElementDof> &element_dofs,
                               const BarycentricGradients &grad_lambda,
                               double jacobian_determinant)
  {
    const auto &entry = GetEntry(finite_element, element_dofs);
    const auto geometry = BuildAffineGeometryTensors(grad_lambda);
    AffineCommonNDCurlBlocks result;
    result.standard_enrichment.SetSize(entry.standard_size, entry.enrichment_size);
    result.enrichment_standard.SetSize(entry.enrichment_size, entry.standard_size);
    result.enrichment_enrichment.SetSize(entry.enrichment_size);
    result.standard_enrichment_estimated_absolute_error.SetSize(entry.standard_size,
                                                                entry.enrichment_size);
    result.enrichment_enrichment_estimated_absolute_error.SetSize(entry.enrichment_size);
    result.standard_enrichment = 0.0;
    result.enrichment_standard = 0.0;
    result.enrichment_enrichment = 0.0;
    result.standard_enrichment_estimated_absolute_error = 0.0;
    result.enrichment_enrichment_estimated_absolute_error = 0.0;

    const int rotational_size = static_cast<int>(entry.rotational_indices.size());
    for (int standard = 0; standard < entry.standard_size; standard++)
    {
      for (int rotational = 0; rotational < rotational_size; rotational++)
      {
        const int enrichment = entry.rotational_indices[rotational];
        const auto &reference =
            entry.standard_enrichment[static_cast<std::size_t>(standard) * rotational_size +
                                      rotational];
        const double value = ContractAffineTensor(reference.value, geometry.curl_curl,
                                                  jacobian_determinant, false);
        const double error =
            ContractAffineTensor(reference.estimated_absolute_error, geometry.curl_curl,
                                 jacobian_determinant, true);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND curl contraction produced an invalid coupling entry!");
        }
        result.standard_enrichment(standard, enrichment) = value;
        result.enrichment_standard(enrichment, standard) = value;
        result.standard_enrichment_estimated_absolute_error(standard, enrichment) = error;
      }
    }
    for (int row = 0; row < rotational_size; row++)
    {
      for (int column = row; column < rotational_size; column++)
      {
        const int enrichment_row = entry.rotational_indices[row];
        const int enrichment_column = entry.rotational_indices[column];
        const auto &reference =
            entry.enrichment_enrichment[EnrichmentPairIndex(row, column, rotational_size)];
        const double value = ContractAffineTensor(reference.value, geometry.curl_curl,
                                                  jacobian_determinant, false);
        const double error =
            ContractAffineTensor(reference.estimated_absolute_error, geometry.curl_curl,
                                 jacobian_determinant, true);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND curl contraction produced an invalid enrichment entry!");
        }
        result.enrichment_enrichment(enrichment_row, enrichment_column) =
            result.enrichment_enrichment(enrichment_column, enrichment_row) = value;
        result.enrichment_enrichment_estimated_absolute_error(enrichment_row,
                                                              enrichment_column) =
            result.enrichment_enrichment_estimated_absolute_error(enrichment_column,
                                                                  enrichment_row) = error;
      }
    }
    result.leaf_count = entry.leaf_count;
    result.maximum_subdivision_depth = entry.maximum_subdivision_depth;
    return result;
  }

  std::size_t Size() const
  {
    std::size_t size = 0;
    for (const auto &[key, entry] : entries)
    {
      (void)key;
      size += entry.standard_enrichment.size() + entry.enrichment_enrichment.size();
    }
    return size;
  }
  std::size_t Hits() const { return hits; }
  double GenerationTime() const { return generation_time; }
};

void ApplyCommonAffineNDMassBlocks(const ElementDofMap &element_dofs,
                                   const mfem::FiniteElement &h1_fe,
                                   const mfem::FiniteElement &nd_fe,
                                   mfem::ElementTransformation &transformation,
                                   const AffineCommonNDMassBlocks &common,
                                   ElementEnrichmentMatrices &enrichment,
                                   ElementStandardEnrichmentMatrices &coupling)
{
  const int standard_h1_size = h1_fe.GetDof();
  const int standard_nd_size = nd_fe.GetDof();
  const int enrichment_h1_size = static_cast<int>(element_dofs.h1.size());
  const int enrichment_nd_size = static_cast<int>(element_dofs.nd.size());
  if (common.standard_enrichment.Height() != standard_nd_size ||
      common.standard_enrichment.Width() != enrichment_nd_size ||
      common.enrichment_standard.Height() != enrichment_nd_size ||
      common.enrichment_standard.Width() != standard_nd_size ||
      common.enrichment_enrichment.Height() != enrichment_nd_size ||
      common.enrichment_enrichment.Width() != enrichment_nd_size ||
      common.standard_enrichment_estimated_absolute_error.Height() != standard_nd_size ||
      common.standard_enrichment_estimated_absolute_error.Width() != enrichment_nd_size ||
      common.enrichment_enrichment_estimated_absolute_error.Height() !=
          enrichment_nd_size ||
      common.enrichment_enrichment_estimated_absolute_error.Width() != enrichment_nd_size)
  {
    throw std::invalid_argument(
        "Common affine ND mass blocks have inconsistent element dimensions!");
  }

  enrichment.nd_mass = common.enrichment_enrichment;
  enrichment.nd_mass_estimated_absolute_error =
      common.enrichment_enrichment_estimated_absolute_error;
  coupling.nd_mass_standard_enrichment = common.standard_enrichment;
  coupling.nd_mass_enrichment_standard = common.enrichment_standard;
  coupling.nd_mass_estimated_absolute_error =
      common.standard_enrichment_estimated_absolute_error;

  const auto h1_to_nd = BuildH1ToNDMap(element_dofs);
  for (int row = 0; row < enrichment_h1_size; row++)
  {
    for (int column = 0; column < enrichment_h1_size; column++)
    {
      enrichment.h1_diffusion(row, column) =
          enrichment.nd_mass(h1_to_nd[row], h1_to_nd[column]);
      enrichment.h1_diffusion_estimated_absolute_error(row, column) =
          enrichment.nd_mass_estimated_absolute_error(h1_to_nd[row], h1_to_nd[column]);
    }
  }

  mfem::DenseMatrix discrete_gradient(standard_nd_size, standard_h1_size);
  nd_fe.ProjectGrad(h1_fe, transformation, discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      long double value = 0.0L;
      long double error = 0.0L;
      const int enrichment_nd = h1_to_nd[enrichment_h1];
      for (int standard_nd = 0; standard_nd < standard_nd_size; standard_nd++)
      {
        const double coefficient = discrete_gradient(standard_nd, standard_h1);
        if (!std::isfinite(coefficient))
        {
          throw std::runtime_error(
              "Common affine ND mass discrete gradient contains a nonfinite entry!");
        }
        value += static_cast<long double>(coefficient) *
                 coupling.nd_mass_standard_enrichment(standard_nd, enrichment_nd);
        error += std::abs(static_cast<long double>(coefficient)) *
                 coupling.nd_mass_estimated_absolute_error(standard_nd, enrichment_nd);
      }
      const double entry = static_cast<double>(value);
      const double entry_error = static_cast<double>(error);
      if (!std::isfinite(entry) || !std::isfinite(entry_error) || entry_error < 0.0)
      {
        throw std::runtime_error(
            "Common affine ND mass produced an invalid H1 coupling entry!");
      }
      coupling.h1_standard_enrichment(standard_h1, enrichment_h1) = entry;
      coupling.h1_enrichment_standard(enrichment_h1, standard_h1) = entry;
      coupling.h1_estimated_absolute_error(standard_h1, enrichment_h1) = entry_error;
    }
  }

  if (enrichment.total_quadrature_leaf_count >
      std::numeric_limits<std::size_t>::max() - common.leaf_count)
  {
    throw std::overflow_error("Common affine ND mass quadrature leaf count overflow!");
  }
  enrichment.total_quadrature_leaf_count += common.leaf_count;
  enrichment.maximum_subdivision_depth =
      std::max(enrichment.maximum_subdivision_depth, common.maximum_subdivision_depth);
}

class AffineScalarReferenceTable
{
private:
  struct PairKey
  {
    HigherOrderBasis row;
    HigherOrderBasis column;
  };

  struct PairKeyLess
  {
    bool operator()(const PairKey &a, const PairKey &b) const
    {
      return BasisLess(a.row, b.row) ||
             (!BasisLess(b.row, a.row) && BasisLess(a.column, b.column));
    }
  };

  struct StandardKey
  {
    const mfem::FiniteElement *finite_element;
    HigherOrderBasis basis;
  };

  struct StandardKeyLess
  {
    bool operator()(const StandardKey &a, const StandardKey &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return BasisLess(a.basis, b.basis);
    }
  };

  struct StandardPatternKey
  {
    const mfem::FiniteElement *finite_element;
    std::vector<HigherOrderBasis> bases;
  };

  struct StandardPatternKeyLess
  {
    bool operator()(const StandardPatternKey &a, const StandardPatternKey &b) const
    {
      const auto pointer_less = std::less<const mfem::FiniteElement *>{};
      if (pointer_less(a.finite_element, b.finite_element))
      {
        return true;
      }
      if (pointer_less(b.finite_element, a.finite_element))
      {
        return false;
      }
      return std::lexicographical_compare(a.bases.begin(), a.bases.end(), b.bases.begin(),
                                          b.bases.end(), BasisLess);
    }
  };

  AdaptiveAssemblyOptions options;
  std::map<PairKey, AdaptiveQuadratureResult, PairKeyLess> pair_entries;
  std::map<StandardKey, std::vector<AdaptiveQuadratureResult>, StandardKeyLess>
      standard_entries;
  std::map<StandardPatternKey, std::vector<const AdaptiveQuadratureResult *>,
           StandardPatternKeyLess>
      standard_patterns;
  std::size_t hits = 0;
  double generation_time = 0.0;

  static void Validate(const AdaptiveQuadratureResult &entry, std::string_view quantity)
  {
    if (!entry.converged || !std::isfinite(entry.value) ||
        !std::isfinite(entry.estimated_absolute_error) ||
        entry.estimated_absolute_error < 0.0)
    {
      throw std::runtime_error(
          fmt::format("Cached affine singular {} integral did not converge: value = "
                      "{:.17g}, estimated absolute error = {:.17g}, leaves = {}, "
                      "maximum depth = {}!",
                      quantity, entry.value, entry.estimated_absolute_error,
                      entry.leaf_count, entry.maximum_subdivision_depth));
    }
  }

public:
  explicit AffineScalarReferenceTable(const AdaptiveAssemblyOptions &options_in)
    : options(options_in)
  {
  }

  AdaptiveQuadratureResult GetPair(const HigherOrderBasis &row_basis,
                                   const HigherOrderBasis &column_basis)
  {
    const auto canonical = CanonicalizeReferenceBasisPair(ReferenceBasis{row_basis},
                                                          ReferenceBasis{column_basis});
    const auto &row = std::get<HigherOrderBasis>(canonical.row_basis);
    const auto &column = std::get<HigherOrderBasis>(canonical.column_basis);
    const PairKey key{row, column};
    const auto cached = pair_entries.find(key);
    if (cached != pair_entries.end())
    {
      hits++;
      return cached->second;
    }
    const auto generation_start = std::chrono::steady_clock::now();
    const auto entry =
        IntegrateConfigured(options.quadrature_order, options,
                            [&](const BarycentricPoint &lambda)
                            {
                              return EvaluateHigherOrderGradientPotential(lambda, row) *
                                     EvaluateHigherOrderGradientPotential(lambda, column);
                            });
    Validate(entry, "H1 enrichment mass");
    const auto result = pair_entries.emplace(key, entry).first->second;
    generation_time +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
            .count();
    return result;
  }

private:
  const std::vector<AdaptiveQuadratureResult> &
  GetStandardEntries(const mfem::FiniteElement &finite_element,
                     const HigherOrderBasis &basis)
  {
    const int standard_size = finite_element.GetDof();
    const StandardKey key{&finite_element, basis};
    const auto cached = standard_entries.find(key);
    if (cached != standard_entries.end())
    {
      hits += static_cast<std::size_t>(standard_size);
      return cached->second;
    }
    const auto generation_start = std::chrono::steady_clock::now();
    mfem::IntegrationPoint point;
    mfem::Vector standard_shape(standard_size);
    const auto integral = IntegrateConfigured(
        options.quadrature_order, options, static_cast<std::size_t>(standard_size),
        [&](const BarycentricPoint &lambda, std::vector<double> &value)
        {
          point.Set3(lambda[1], lambda[2], lambda[3]);
          finite_element.CalcShape(point, standard_shape);
          const double singular = EvaluateHigherOrderGradientPotential(lambda, basis);
          for (int standard_dof = 0; standard_dof < standard_size; standard_dof++)
          {
            value[standard_dof] = standard_shape[standard_dof] * singular;
          }
        });
    if (!integral.converged)
    {
      throw std::runtime_error(fmt::format(
          "Cached affine standard-enrichment H1 mass integrals did not converge: "
          "family = {}, order = {}, nu = {:.17g}, nodes = [{}, {}, {}, {}], "
          "interpolation = [{}, {}, {}, {}], components = {}, leaves = {}, maximum "
          "depth = {}!",
          static_cast<int>(basis.family), basis.order, basis.nu, basis.nodes[0],
          basis.nodes[1], basis.nodes[2], basis.nodes[3], basis.interpolation_indices[0],
          basis.interpolation_indices[1], basis.interpolation_indices[2],
          basis.interpolation_indices[3], standard_size, integral.leaf_count,
          integral.maximum_subdivision_depth));
    }
    std::vector<AdaptiveQuadratureResult> basis_entries(standard_size);
    for (int standard_dof = 0; standard_dof < standard_size; standard_dof++)
    {
      basis_entries[standard_dof] = {
          integral.value[standard_dof], integral.estimated_absolute_error[standard_dof],
          integral.leaf_count, integral.maximum_subdivision_depth, true};
      Validate(basis_entries[standard_dof], "standard-enrichment H1 mass");
    }
    hits += static_cast<std::size_t>(standard_size - 1);
    const auto &result =
        standard_entries.emplace(key, std::move(basis_entries)).first->second;
    generation_time +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start)
            .count();
    return result;
  }

public:
  const std::vector<const AdaptiveQuadratureResult *> &
  GetStandardPattern(const mfem::FiniteElement &finite_element,
                     const std::vector<ElementDof> &element_dofs)
  {
    const int standard_size = finite_element.GetDof();
    StandardPatternKey key{&finite_element, {}};
    key.bases.reserve(element_dofs.size());
    for (const auto &element_dof : element_dofs)
    {
      key.bases.push_back(element_dof.basis);
    }
    const auto cached = standard_patterns.find(key);
    if (cached != standard_patterns.end())
    {
      hits += static_cast<std::size_t>(standard_size) * element_dofs.size();
      return cached->second;
    }

    std::vector<const AdaptiveQuadratureResult *> pattern(
        static_cast<std::size_t>(standard_size) * element_dofs.size());
    for (std::size_t enrichment = 0; enrichment < element_dofs.size(); enrichment++)
    {
      const auto &basis_entries =
          GetStandardEntries(finite_element, element_dofs[enrichment].basis);
      for (int standard = 0; standard < standard_size; standard++)
      {
        pattern[static_cast<std::size_t>(standard) * element_dofs.size() + enrichment] =
            &basis_entries[standard];
      }
    }
    return standard_patterns.emplace(std::move(key), std::move(pattern)).first->second;
  }

  std::size_t Size() const
  {
    std::size_t size = pair_entries.size();
    for (const auto &[key, basis_entries] : standard_entries)
    {
      (void)key;
      size += basis_entries.size();
    }
    return size;
  }
  std::size_t Hits() const { return hits; }
  double GenerationTime() const { return generation_time; }
};

}  // namespace

BarycentricGradients
GetAffineBarycentricGradients(mfem::ElementTransformation &transformation,
                              double &jacobian_determinant)
{
  if (transformation.GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
      transformation.GetDimension() != 3 || transformation.GetSpaceDim() != 3 ||
      !IsAffineElementTransformation(transformation))
  {
    throw std::invalid_argument(
        "Singular field evaluation requires an affine tetrahedron in three dimensions!");
  }

  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  return GetBarycentricGradients(transformation, center, jacobian_determinant);
}

TriangleBarycentricGradients
GetAffineTriangleBarycentricGradients(mfem::ElementTransformation &transformation,
                                      double &jacobian_determinant)
{
  if (transformation.GetGeometryType() != mfem::Geometry::TRIANGLE ||
      transformation.GetDimension() != 2 || transformation.GetSpaceDim() != 2 ||
      !IsAffineElementTransformation(transformation))
  {
    throw std::invalid_argument(
        "Singular field evaluation requires an affine triangle in two dimensions!");
  }

  mfem::IntegrationPoint center;
  center.Set2(1.0 / 3.0, 1.0 / 3.0);
  return GetTriangleBarycentricGradients(transformation, center, jacobian_determinant);
}

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichmentProlongation(MPI_Comm comm, const TrueDofMap &dofs)
{
  return BuildEnrichmentProlongation(comm, dofs);
}

std::unique_ptr<mfem::Vector> AssembleParallelEnrichmentVector(MPI_Comm comm,
                                                               const TrueDofMap &dofs,
                                                               const mfem::Vector &local)
{
  auto prolongation = BuildParallelEnrichmentProlongation(comm, dofs);
  if (local.Size() != prolongation->Height())
  {
    throw std::invalid_argument(
        "Cannot parallel-assemble an inconsistent local enrichment vector!");
  }
  auto assembled = std::make_unique<mfem::Vector>(prolongation->Width());
  prolongation->MultTranspose(local, *assembled);
  return assembled;
}

std::unique_ptr<mfem::SparseMatrix>
BuildLocalEnrichmentInterpolant(const mfem::HypreParMatrix &true_interpolant,
                                const mfem::ParFiniteElementSpace &standard_fespace,
                                const TrueDofMap &enrichment_numbering)
{
  const MPI_Comm comm = standard_fespace.GetComm();
  const auto *standard_prolongation = standard_fespace.Dof_TrueDof_Matrix();
  if (!standard_prolongation ||
      true_interpolant.Height() != standard_prolongation->Width() ||
      true_interpolant.GetGlobalNumRows() != standard_prolongation->GetGlobalNumCols() ||
      true_interpolant.GetGlobalNumCols() != enrichment_numbering.global_size)
  {
    throw std::invalid_argument(
        "Cannot localize an inconsistent singular enrichment interpolant!");
  }

  auto enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, enrichment_numbering);
  std::unique_ptr<mfem::HypreParMatrix> enrichment_restriction(
      enrichment_prolongation->Transpose());
  std::unique_ptr<mfem::HypreParMatrix> local_to_true(
      mfem::ParMult(standard_prolongation, &true_interpolant, true));
  std::unique_ptr<mfem::HypreParMatrix> local(
      mfem::ParMult(local_to_true.get(), enrichment_restriction.get(), true));
  if (!local || local->Height() != standard_fespace.GetVSize() ||
      local->Width() != static_cast<int>(enrichment_numbering.local_size))
  {
    throw std::runtime_error("Failed to localize a singular enrichment interpolant!");
  }

  mfem::SparseMatrix diagonal, off_diagonal;
  HYPRE_BigInt *off_diagonal_columns = nullptr;
  local->GetDiag(diagonal);
  local->GetOffd(off_diagonal, off_diagonal_columns);
  // P_e^T has one column for every rank-local copy of a shared true DOF.
  // Consequently P_s D P_e^T repeats the same interpolant in off-rank local
  // columns. Each rank assembles only its own elements and canonical local
  // enrichment columns; the final P_e^T vector reduction sums those element
  // contributions into the shared true DOF.
  return std::make_unique<mfem::SparseMatrix>(diagonal);
}

namespace
{

mfem::DenseMatrix
ExtractElementInterpolant(const mfem::SparseMatrix &local_interpolant,
                          const mfem::Array<int> &signed_standard_dofs,
                          const mfem::DofTransformation &standard_dof_transformation,
                          const mfem::Array<int> &enrichment_dofs)
{
  mfem::DenseMatrix element(signed_standard_dofs.Size(), enrichment_dofs.Size());
  element = 0.0;
  std::map<int, int> enrichment_column;
  for (int column = 0; column < enrichment_dofs.Size(); column++)
  {
    const auto [entry, inserted] =
        enrichment_column.emplace(enrichment_dofs[column], column);
    if (!inserted || entry->first < 0 || entry->first >= local_interpolant.Width())
    {
      throw std::invalid_argument(
          "Element singular interpolant has invalid enrichment DOFs!");
    }
  }

  const auto *offsets = local_interpolant.HostReadI();
  const auto *columns = local_interpolant.HostReadJ();
  const auto *values = local_interpolant.HostReadData();
  for (int row = 0; row < signed_standard_dofs.Size(); row++)
  {
    const int signed_dof = signed_standard_dofs[row];
    if (signed_dof == std::numeric_limits<int>::min())
    {
      throw std::invalid_argument(
          "Element singular interpolant has an invalid signed standard DOF!");
    }
    const int standard_dof = signed_dof >= 0 ? signed_dof : -1 - signed_dof;
    const double sign = signed_dof >= 0 ? 1.0 : -1.0;
    if (standard_dof < 0 || standard_dof >= local_interpolant.Height())
    {
      throw std::invalid_argument(
          "Element singular interpolant has an out-of-range standard DOF!");
    }
    for (int entry = offsets[standard_dof]; entry < offsets[standard_dof + 1]; entry++)
    {
      const auto column = enrichment_column.find(columns[entry]);
      if (column != enrichment_column.end())
      {
        element(row, column->second) = sign * values[entry];
      }
    }
  }

  for (int column = 0; column < element.Width(); column++)
  {
    mfem::Vector coefficients(element.GetColumn(column), element.Height());
    standard_dof_transformation.InvTransformPrimal(coefficients);
  }
  return element;
}

template <typename Topology>
std::vector<LocalInterpolatedNDDiagonals> AssembleLocalInterpolatedNDDomainDiagonalsImpl(
    const Topology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches)
{
  auto *mesh = nd_fespace.GetMesh();
  if (!mesh || topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      local_interpolant.Height() != nd_fespace.GetVSize() ||
      local_interpolant.Width() != static_cast<int>(topology.nd_dofs.size()) ||
      material_batches.empty())
  {
    throw std::invalid_argument(
        "Interpolated singular ND diagonal assembly received inconsistent inputs!");
  }
  for (const auto &materials : material_batches)
  {
    if (materials.size() != static_cast<std::size_t>(mesh->GetNE()))
    {
      throw std::invalid_argument(
          "Interpolated singular ND diagonal materials do not cover the mesh!");
    }
    for (const auto &material : materials)
    {
      ValidateMaterialCoefficients(material);
    }
  }

  std::vector<LocalInterpolatedNDDiagonals> result(material_batches.size());
  for (auto &diagonal : result)
  {
    diagonal.mass.SetSize(local_interpolant.Width());
    diagonal.curl_curl.SetSize(local_interpolant.Width());
    diagonal.mass = 0.0;
    diagonal.curl_curl = 0.0;
  }

  mfem::VectorFEMassIntegrator mass_integrator;
  mfem::CurlCurlIntegrator curl_integrator;
  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element].nd;
    if (element_dofs.empty())
    {
      continue;
    }
    const auto enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs, topology.nd_dofs.size());
    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    nd_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    const auto interpolant = ExtractElementInterpolant(local_interpolant, standard_dofs,
                                                       dof_transformation, enrichment_dofs);

    auto &transformation = *mesh->GetElementTransformation(element);
    mfem::DenseMatrix mass, curl_curl;
    mass_integrator.AssembleElementMatrix(*nd_fespace.GetFE(element), transformation, mass);
    curl_integrator.AssembleElementMatrix(*nd_fespace.GetFE(element), transformation,
                                          curl_curl);
    mfem::Vector work(interpolant.Height());
    for (int column = 0; column < interpolant.Width(); column++)
    {
      mfem::Vector coefficients(const_cast<double *>(interpolant.GetColumn(column)),
                                interpolant.Height());
      mass.Mult(coefficients, work);
      const double mass_energy = coefficients * work;
      curl_curl.Mult(coefficients, work);
      const double curl_energy = coefficients * work;
      const int enrichment = enrichment_dofs[column];
      for (std::size_t batch = 0; batch < material_batches.size(); batch++)
      {
        const auto &material = material_batches[batch][element];
        result[batch].mass[enrichment] += material.electric * mass_energy;
        result[batch].curl_curl[enrichment] += material.inverse_magnetic * curl_energy;
      }
    }
  }
  return result;
}

int SetBoundaryIntegrationPoint(mfem::Mesh &mesh, int boundary,
                                const mfem::IntegrationPoint &boundary_point,
                                mfem::FaceElementTransformations &face,
                                mfem::IsoparametricTransformation &element1,
                                mfem::IsoparametricTransformation &element2);

template <typename Topology>
mfem::Vector AssembleLocalInterpolatedNDBoundaryDiagonalImpl(
    const Topology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::map<int, double> &boundary_coefficients)
{
  auto *mesh = nd_fespace.GetMesh();
  const int dimension = mesh ? mesh->Dimension() : 0;
  if (!mesh || mesh->SpaceDimension() != dimension || (dimension != 2 && dimension != 3) ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      local_interpolant.Height() != nd_fespace.GetVSize() ||
      local_interpolant.Width() != static_cast<int>(topology.nd_dofs.size()))
  {
    throw std::invalid_argument(
        "Interpolated singular ND boundary diagonal received inconsistent inputs!");
  }
  for (const auto &[attribute, coefficient] : boundary_coefficients)
  {
    if (attribute <= 0 || !std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Interpolated singular ND boundary diagonal requires positive attributes "
          "and finite coefficients!");
    }
  }

  mfem::Vector result(local_interpolant.Width());
  result = 0.0;
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const auto coefficient = boundary_coefficients.find(mesh->GetBdrAttribute(boundary));
    if (coefficient == boundary_coefficients.end() || coefficient->second == 0.0)
    {
      continue;
    }
    mfem::FaceElementTransformations face;
    mfem::IsoparametricTransformation element1, element2, boundary_transformation;
    mesh->GetBdrElementTransformation(boundary, &boundary_transformation);
    const auto geometry = boundary_transformation.GetGeometryType();
    const auto &center = mfem::Geometries.GetCenter(geometry);
    const int element =
        SetBoundaryIntegrationPoint(*mesh, boundary, center, face, element1, element2);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Interpolated singular ND boundary diagonal has no adjacent element!");
    }
    const auto &element_dofs = topology.elements[element].nd;
    if (element_dofs.empty())
    {
      continue;
    }

    const auto enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs, topology.nd_dofs.size());
    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    nd_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    const auto interpolant = ExtractElementInterpolant(local_interpolant, standard_dofs,
                                                       dof_transformation, enrichment_dofs);
    const auto *finite_element = nd_fespace.GetFE(element);
    if (!finite_element || !face.Elem1 || geometry == mfem::Geometry::INVALID)
    {
      throw std::runtime_error(
          "Interpolated singular ND boundary diagonal has invalid transformations!");
    }

    const int quadrature_order = fem::DefaultIntegrationOrder::Get(boundary_transformation);
    const auto &rule = mfem::IntRules.Get(geometry, quadrature_order);
    mfem::DenseMatrix standard_value(finite_element->GetDof(), dimension);
    mfem::Vector normal_or_tangent(dimension), value(dimension);
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &point = rule.IntPoint(q);
      const int current =
          SetBoundaryIntegrationPoint(*mesh, boundary, point, face, element1, element2);
      if (current != element)
      {
        throw std::runtime_error(
            "Interpolated singular ND boundary diagonal changed adjacent element!");
      }
      finite_element->CalcPhysVShape(*face.Elem1, standard_value);
      boundary_transformation.SetIntPoint(&point);
      const auto &jacobian = boundary_transformation.Jacobian();
      const double weight =
          coefficient->second * point.weight * boundary_transformation.Weight();
      if (!std::isfinite(weight) || weight < 0.0)
      {
        throw std::runtime_error(
            "Interpolated singular ND boundary diagonal found invalid quadrature!");
      }

      if (dimension == 3)
      {
        mfem::CalcOrtho(jacobian, normal_or_tangent);
      }
      else
      {
        if (jacobian.Height() != 2 || jacobian.Width() != 1)
        {
          throw std::runtime_error(
              "Interpolated triangular ND boundary diagonal has an invalid tangent!");
        }
        normal_or_tangent[0] = jacobian(0, 0);
        normal_or_tangent[1] = jacobian(1, 0);
      }
      normal_or_tangent /= normal_or_tangent.Norml2();

      for (int enrichment = 0; enrichment < interpolant.Width(); enrichment++)
      {
        mfem::Vector coefficients(const_cast<double *>(interpolant.GetColumn(enrichment)),
                                  interpolant.Height());
        standard_value.MultTranspose(coefficients, value);
        const double trace = value * normal_or_tangent;
        const double tangential_norm =
            dimension == 3 ? value * value - trace * trace : trace * trace;
        result[enrichment_dofs[enrichment]] += weight * tangential_norm;
      }
    }
  }
  return result;
}

double SparseRowInnerProduct(const mfem::SparseMatrix &first, int first_row,
                             const mfem::SparseMatrix &second, int second_row)
{
  const auto *first_offsets = first.HostReadI();
  const auto *first_columns = first.HostReadJ();
  const auto *first_values = first.HostReadData();
  const auto *second_offsets = second.HostReadI();
  const auto *second_columns = second.HostReadJ();
  const auto *second_values = second.HostReadData();
  int i = first_offsets[first_row];
  int j = second_offsets[second_row];
  const int i_end = first_offsets[first_row + 1];
  const int j_end = second_offsets[second_row + 1];
  double product = 0.0;
  while (i < i_end && j < j_end)
  {
    if (first_columns[i] < second_columns[j])
    {
      i++;
    }
    else if (second_columns[j] < first_columns[i])
    {
      j++;
    }
    else
    {
      product += first_values[i] * second_values[j];
      i++;
      j++;
    }
  }
  return product;
}

ElementEnrichmentMatrices AssembleAffineElementEnrichmentMatrices(
    const ElementDofMap &element_dofs, const BarycentricGradients &grad_lambda,
    double jacobian_determinant, const AdaptiveAssemblyOptions &options,
    AdaptiveReferenceTable *reference_table, AffineScalarReferenceTable *scalar_table,
    const AffineCommonNDMassBlocks *common_mass,
    const AffineCommonNDCurlBlocks *common_curl)
{
  ValidateInputs(element_dofs, grad_lambda, jacobian_determinant, options);
  const auto h1_to_nd = BuildH1ToNDMap(element_dofs);

  ElementEnrichmentMatrices result;
  const int h1_size = static_cast<int>(element_dofs.h1.size());
  const int nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_diffusion.SetSize(h1_size);
  result.h1_diffusion_estimated_absolute_error.SetSize(h1_size);
  result.h1_mass.SetSize(h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(h1_size);
  result.nd_mass.SetSize(nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(nd_size);
  result.nd_curl_curl.SetSize(nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(nd_size);

  for (int row = 0; row < nd_size; row++)
  {
    for (int column = row; column < nd_size; column++)
    {
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      const double mass = common_mass ? common_mass->enrichment_enrichment(row, column)
                                      : std::numeric_limits<double>::quiet_NaN();
      const double mass_error =
          common_mass
              ? common_mass->enrichment_enrichment_estimated_absolute_error(row, column)
              : std::numeric_limits<double>::quiet_NaN();
      const bool curl_free =
          IsGradientFamily(row_basis.family) || IsGradientFamily(column_basis.family);
      double curl_curl =
          common_curl ? common_curl->enrichment_enrichment(row, column) : 0.0;
      double curl_curl_error =
          common_curl
              ? common_curl->enrichment_enrichment_estimated_absolute_error(row, column)
              : 0.0;
      double final_mass = mass;
      double final_mass_error = mass_error;
      if (!common_mass || (!common_curl && !curl_free))
      {
        const auto reference =
            reference_table
                ? reference_table->Get(ReferenceBasis{row_basis},
                                       ReferenceBasis{column_basis})
                : ComputeReferenceIntegralConfigured(ReferenceBasis{row_basis},
                                                     ReferenceBasis{column_basis}, options);
        if (!reference.converged)
        {
          throw std::runtime_error(
              "Singular element adaptive reference tensor did not converge!");
        }
        RecordQuadratureStatistics(reference, result);
        if (!common_mass)
        {
          final_mass = ContractMass(reference.integral, grad_lambda, jacobian_determinant);
          final_mass_error =
              ContractMassError(reference, grad_lambda, jacobian_determinant);
        }
        if (!common_curl && !curl_free)
        {
          curl_curl =
              ContractCurlCurl(reference.integral, grad_lambda, jacobian_determinant);
          curl_curl_error =
              ContractCurlCurlError(reference, grad_lambda, jacobian_determinant);
        }
      }
      if (!std::isfinite(final_mass) || !std::isfinite(final_mass_error) ||
          !std::isfinite(curl_curl) || !std::isfinite(curl_curl_error) ||
          final_mass_error < 0.0 || curl_curl_error < 0.0)
      {
        throw std::runtime_error(
            "Singular element assembly produced a nonfinite matrix entry!");
      }

      result.nd_mass(row, column) = result.nd_mass(column, row) = final_mass;
      result.nd_mass_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(column, row) = final_mass_error;
      result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = curl_curl;
      result.nd_curl_curl_estimated_absolute_error(row, column) =
          result.nd_curl_curl_estimated_absolute_error(column, row) = curl_curl_error;
    }
  }

  for (int row = 0; row < h1_size; row++)
  {
    for (int column = 0; column < h1_size; column++)
    {
      result.h1_diffusion(row, column) = result.nd_mass(h1_to_nd[row], h1_to_nd[column]);
      result.h1_diffusion_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(h1_to_nd[row], h1_to_nd[column]);
    }
  }
  for (int row = 0; row < h1_size; row++)
  {
    for (int column = row; column < h1_size; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      const auto integral = [&]()
      {
        if (scalar_table)
        {
          auto reference = scalar_table->GetPair(row_basis, column_basis);
          reference.value *= jacobian_determinant;
          reference.estimated_absolute_error *= jacobian_determinant;
          return reference;
        }
        return IntegrateConfigured(
            options.quadrature_order, options,
            [&](const BarycentricPoint &lambda)
            {
              return jacobian_determinant *
                     EvaluateHigherOrderGradientPotential(lambda, row_basis) *
                     EvaluateHigherOrderGradientPotential(lambda, column_basis);
            });
      }();
      if (!integral.converged || !std::isfinite(integral.value) ||
          !std::isfinite(integral.estimated_absolute_error) ||
          integral.estimated_absolute_error < 0.0)
      {
        throw std::runtime_error(
            "Singular H1 mass integration did not meet its requested tolerance!");
      }
      result.h1_mass(row, column) = result.h1_mass(column, row) = integral.value;
      result.h1_mass_estimated_absolute_error(row, column) =
          result.h1_mass_estimated_absolute_error(column, row) =
              integral.estimated_absolute_error;
      result.total_quadrature_leaf_count += integral.leaf_count;
      result.maximum_subdivision_depth =
          std::max(result.maximum_subdivision_depth, integral.maximum_subdivision_depth);
    }
  }
  return result;
}

}  // namespace

ElementEnrichmentMatrices AssembleElementEnrichmentMatrices(
    const ElementDofMap &element_dofs, const BarycentricGradients &grad_lambda,
    double jacobian_determinant, const AdaptiveAssemblyOptions &options)
{
  AdaptiveReferenceTable reference_table(options);
  AffineScalarReferenceTable scalar_table(options);
  return AssembleAffineElementEnrichmentMatrices(
      element_dofs, grad_lambda, jacobian_determinant, options, &reference_table,
      &scalar_table, nullptr, nullptr);
}

ElementEnrichmentMatrices
AssembleElementEnrichmentMatrices(const ElementDofMap &element_dofs,
                                  mfem::ElementTransformation &transformation,
                                  const AdaptiveAssemblyOptions &options)
{
  if (IsAffineElementTransformation(transformation))
  {
    double jacobian_determinant;
    const auto grad_lambda =
        GetAffineBarycentricGradients(transformation, jacobian_determinant);
    return AssembleElementEnrichmentMatrices(element_dofs, grad_lambda,
                                             jacobian_determinant, options);
  }

  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetBarycentricGradients(transformation, center, center_jacobian_determinant);
  ValidateInputs(element_dofs, center_grad_lambda, center_jacobian_determinant, options);
  const auto h1_to_nd = BuildH1ToNDMap(element_dofs);

  ElementEnrichmentMatrices result;
  const int h1_size = static_cast<int>(element_dofs.h1.size());
  const int nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_diffusion.SetSize(h1_size);
  result.h1_diffusion_estimated_absolute_error.SetSize(h1_size);
  result.h1_mass.SetSize(h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(h1_size);
  result.nd_mass.SetSize(nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(nd_size);
  result.nd_curl_curl.SetSize(nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(nd_size);

  mfem::IntegrationPoint point;
  for (int row = 0; row < nd_size; row++)
  {
    for (int column = row; column < nd_size; column++)
    {
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      const auto mass = IntegrateConfigured(
          options.quadrature_order, options,
          [&](const BarycentricPoint &lambda)
          {
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant;
            const auto grad_lambda =
                GetBarycentricGradients(transformation, point, jacobian_determinant);
            const auto row_value = EvaluateHigherOrderBasis(lambda, grad_lambda, row_basis);
            const auto column_value =
                EvaluateHigherOrderBasis(lambda, grad_lambda, column_basis);
            return jacobian_determinant * Dot(row_value.value, column_value.value);
          });
      ValidateAdaptiveH1Entry(mass, "curved ND mass", row, column, column_basis);
      RecordQuadratureStatistics(mass, result);
      result.nd_mass(row, column) = result.nd_mass(column, row) = mass.value;
      result.nd_mass_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;

      const bool curl_free =
          IsGradientFamily(row_basis.family) || IsGradientFamily(column_basis.family);
      if (curl_free)
      {
        result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(row, column) =
            result.nd_curl_curl_estimated_absolute_error(column, row) = 0.0;
        continue;
      }
      const auto curl_curl = IntegrateConfigured(
          options.quadrature_order, options,
          [&](const BarycentricPoint &lambda)
          {
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant;
            const auto grad_lambda =
                GetBarycentricGradients(transformation, point, jacobian_determinant);
            const auto row_value = EvaluateHigherOrderBasis(lambda, grad_lambda, row_basis);
            const auto column_value =
                EvaluateHigherOrderBasis(lambda, grad_lambda, column_basis);
            return jacobian_determinant * Dot(row_value.curl, column_value.curl);
          });
      ValidateAdaptiveH1Entry(curl_curl, "curved ND curl-curl", row, column, column_basis);
      RecordQuadratureStatistics(curl_curl, result);
      result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(row, column) =
          result.nd_curl_curl_estimated_absolute_error(column, row) =
              curl_curl.estimated_absolute_error;
    }
  }

  for (int row = 0; row < h1_size; row++)
  {
    for (int column = 0; column < h1_size; column++)
    {
      result.h1_diffusion(row, column) = result.nd_mass(h1_to_nd[row], h1_to_nd[column]);
      result.h1_diffusion_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(h1_to_nd[row], h1_to_nd[column]);
    }
  }
  for (int row = 0; row < h1_size; row++)
  {
    for (int column = row; column < h1_size; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      const auto mass = IntegrateConfigured(
          options.quadrature_order, options,
          [&](const BarycentricPoint &lambda)
          {
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant;
            (void)GetBarycentricGradients(transformation, point, jacobian_determinant);
            return jacobian_determinant *
                   EvaluateHigherOrderGradientPotential(lambda, row_basis) *
                   EvaluateHigherOrderGradientPotential(lambda, column_basis);
          });
      ValidateAdaptiveH1Entry(mass, "curved H1 mass", row, column, column_basis);
      RecordQuadratureStatistics(mass, result);
      result.h1_mass(row, column) = result.h1_mass(column, row) = mass.value;
      result.h1_mass_estimated_absolute_error(row, column) =
          result.h1_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;
    }
  }
  return result;
}

namespace
{

ElementEnrichmentMatrices AssembleElementEnrichmentMatricesCached(
    const ElementDofMap &element_dofs, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options, AdaptiveReferenceTable &reference_table,
    AffineScalarReferenceTable &scalar_table, const AffineCommonNDMassBlocks *common_mass,
    const AffineCommonNDCurlBlocks *common_curl)
{
  if (!IsAffineElementTransformation(transformation))
  {
    return AssembleElementEnrichmentMatrices(element_dofs, transformation, options);
  }
  double jacobian_determinant;
  const auto grad_lambda =
      GetAffineBarycentricGradients(transformation, jacobian_determinant);
  return AssembleAffineElementEnrichmentMatrices(
      element_dofs, grad_lambda, jacobian_determinant, options, &reference_table,
      &scalar_table, common_mass, common_curl);
}

}  // namespace

ElementEnrichmentMatrices
AssembleTriangleElementEnrichmentMatrices(const TriangleElementDofMap &element_dofs,
                                          const TriangleBarycentricGradients &grad_lambda,
                                          double jacobian_determinant,
                                          const AdaptiveAssemblyOptions &options)
{
  ValidateTriangleInputs(element_dofs, grad_lambda, jacobian_determinant, options);
  const auto h1_to_nd = BuildTriangleH1ToNDMap(element_dofs);
  TriangleBarycentricGradients normalized_grad_lambda;
  const double gradient_scale =
      NormalizeTriangleBarycentricGradients(grad_lambda, normalized_grad_lambda);
  const double mass_prefactor =
      ScaleTriangleJacobianDeterminant(jacobian_determinant, gradient_scale, 2);
  const double curl_curl_prefactor =
      ScaleTriangleJacobianDeterminant(jacobian_determinant, gradient_scale, 4);

  ElementEnrichmentMatrices result;
  const int h1_size = static_cast<int>(element_dofs.h1.size());
  const int nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_diffusion.SetSize(h1_size);
  result.h1_diffusion_estimated_absolute_error.SetSize(h1_size);
  result.h1_mass.SetSize(h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(h1_size);
  result.nd_mass.SetSize(nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(nd_size);
  result.nd_curl_curl.SetSize(nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(nd_size);

  for (int row = 0; row < nd_size; row++)
  {
    for (int column = row; column < nd_size; column++)
    {
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      const int row_node = TriangleSingularNode(row_basis);
      const int column_node = TriangleSingularNode(column_basis);
      const auto mass = IntegrateTriangleDuffy(
          row_node, column_node, "enrichment ND mass", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            const auto row_value = EvaluateAffineNormalizedTriangleBasis(
                lambda, normalized_grad_lambda, mass_prefactor, row_basis);
            const auto column_value = EvaluateAffineNormalizedTriangleBasis(
                lambda, normalized_grad_lambda, mass_prefactor, column_basis);
            return mass_prefactor * Dot(row_value.value, column_value.value);
          });
      result.nd_mass(row, column) = result.nd_mass(column, row) = mass.value;
      result.nd_mass_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;

      const bool curl_free =
          IsGradientFamily(row_basis.family) || IsGradientFamily(column_basis.family);
      if (curl_free)
      {
        result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(row, column) =
            result.nd_curl_curl_estimated_absolute_error(column, row) = 0.0;
        continue;
      }
      const auto curl_curl = IntegrateTriangleDuffy(
          row_node, column_node, "enrichment ND curl-curl", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            const auto row_value = EvaluateAffineNormalizedTriangleBasis(
                lambda, normalized_grad_lambda, mass_prefactor, row_basis);
            const auto column_value = EvaluateAffineNormalizedTriangleBasis(
                lambda, normalized_grad_lambda, mass_prefactor, column_basis);
            return curl_curl_prefactor * row_value.curl * column_value.curl;
          });
      result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(row, column) =
          result.nd_curl_curl_estimated_absolute_error(column, row) =
              curl_curl.estimated_absolute_error;
      result.total_quadrature_leaf_count += curl_curl.point_count;
    }
  }

  for (int row = 0; row < h1_size; row++)
  {
    for (int column = 0; column < h1_size; column++)
    {
      result.h1_diffusion(row, column) = result.nd_mass(h1_to_nd[row], h1_to_nd[column]);
      result.h1_diffusion_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(h1_to_nd[row], h1_to_nd[column]);
    }
  }
  for (int row = 0; row < h1_size; row++)
  {
    for (int column = row; column < h1_size; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      const auto mass = IntegrateTriangleDuffy(
          TriangleSingularNode(row_basis), TriangleSingularNode(column_basis),
          "enrichment H1 mass", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            return jacobian_determinant *
                   EvaluateTriangleNodeGradientPotential(lambda, row_basis.nodes[0],
                                                         row_basis.nodes[1], row_basis.nu) *
                   EvaluateTriangleNodeGradientPotential(lambda, column_basis.nodes[0],
                                                         column_basis.nodes[1],
                                                         column_basis.nu);
          });
      result.h1_mass(row, column) = result.h1_mass(column, row) = mass.value;
      result.h1_mass_estimated_absolute_error(row, column) =
          result.h1_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;
    }
  }
  return result;
}

ElementEnrichmentMatrices
AssembleTriangleElementEnrichmentMatrices(const TriangleElementDofMap &element_dofs,
                                          mfem::ElementTransformation &transformation,
                                          const AdaptiveAssemblyOptions &options)
{
  if (IsAffineElementTransformation(transformation))
  {
    double jacobian_determinant;
    const auto grad_lambda =
        GetAffineTriangleBarycentricGradients(transformation, jacobian_determinant);
    return AssembleTriangleElementEnrichmentMatrices(element_dofs, grad_lambda,
                                                     jacobian_determinant, options);
  }

  mfem::IntegrationPoint center;
  center.Set2(1.0 / 3.0, 1.0 / 3.0);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetTriangleBarycentricGradients(transformation, center, center_jacobian_determinant);
  ValidateTriangleInputs(element_dofs, center_grad_lambda, center_jacobian_determinant,
                         options);
  const auto h1_to_nd = BuildTriangleH1ToNDMap(element_dofs);

  ElementEnrichmentMatrices result;
  const int h1_size = static_cast<int>(element_dofs.h1.size());
  const int nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_diffusion.SetSize(h1_size);
  result.h1_diffusion_estimated_absolute_error.SetSize(h1_size);
  result.h1_mass.SetSize(h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(h1_size);
  result.nd_mass.SetSize(nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(nd_size);
  result.nd_curl_curl.SetSize(nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(nd_size);

  mfem::IntegrationPoint point;
  for (int row = 0; row < nd_size; row++)
  {
    for (int column = row; column < nd_size; column++)
    {
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      const int row_node = TriangleSingularNode(row_basis);
      const int column_node = TriangleSingularNode(column_basis);
      const auto mass = IntegrateTriangleDuffy(
          row_node, column_node, "curved enrichment ND mass", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant;
            const auto grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                                     jacobian_determinant);
            const auto row_value = EvaluateTriangleBasis(lambda, grad_lambda, row_basis);
            const auto column_value =
                EvaluateTriangleBasis(lambda, grad_lambda, column_basis);
            return jacobian_determinant * Dot(row_value.value, column_value.value);
          });
      result.nd_mass(row, column) = result.nd_mass(column, row) = mass.value;
      result.nd_mass_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;

      const bool curl_free =
          IsGradientFamily(row_basis.family) || IsGradientFamily(column_basis.family);
      if (curl_free)
      {
        result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(row, column) =
            result.nd_curl_curl_estimated_absolute_error(column, row) = 0.0;
        continue;
      }
      const auto curl_curl = IntegrateTriangleDuffy(
          row_node, column_node, "curved enrichment ND curl-curl", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant;
            const auto grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                                     jacobian_determinant);
            const auto row_value = EvaluateTriangleBasis(lambda, grad_lambda, row_basis);
            const auto column_value =
                EvaluateTriangleBasis(lambda, grad_lambda, column_basis);
            return jacobian_determinant * row_value.curl * column_value.curl;
          });
      result.nd_curl_curl(row, column) = result.nd_curl_curl(column, row) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(row, column) =
          result.nd_curl_curl_estimated_absolute_error(column, row) =
              curl_curl.estimated_absolute_error;
      result.total_quadrature_leaf_count += curl_curl.point_count;
    }
  }

  for (int row = 0; row < h1_size; row++)
  {
    for (int column = 0; column < h1_size; column++)
    {
      result.h1_diffusion(row, column) = result.nd_mass(h1_to_nd[row], h1_to_nd[column]);
      result.h1_diffusion_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(h1_to_nd[row], h1_to_nd[column]);
    }
  }
  for (int row = 0; row < h1_size; row++)
  {
    for (int column = row; column < h1_size; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      const auto mass = IntegrateTriangleDuffy(
          TriangleSingularNode(row_basis), TriangleSingularNode(column_basis),
          "curved enrichment H1 mass", row, column, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant;
            (void)GetTriangleBarycentricGradients(transformation, point,
                                                  jacobian_determinant);
            return jacobian_determinant *
                   EvaluateTriangleNodeGradientPotential(lambda, row_basis.nodes[0],
                                                         row_basis.nodes[1], row_basis.nu) *
                   EvaluateTriangleNodeGradientPotential(lambda, column_basis.nodes[0],
                                                         column_basis.nodes[1],
                                                         column_basis.nu);
          });
      result.h1_mass(row, column) = result.h1_mass(column, row) = mass.value;
      result.h1_mass_estimated_absolute_error(row, column) =
          result.h1_mass_estimated_absolute_error(column, row) =
              mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;
    }
  }
  return result;
}

namespace
{

ElementStandardEnrichmentMatrices AssembleElementStandardEnrichmentMatricesImpl(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options,
    AffineStandardReferenceTable *standard_reference_table,
    AffineScalarReferenceTable *scalar_reference_table,
    const AffineCommonNDMassBlocks *common_mass,
    const AffineCommonNDCurlBlocks *common_curl)
{
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); };
  const auto setup_start = Clock::now();
  ValidateStandardFiniteElements(h1_fe, nd_fe, transformation);
  const bool affine_geometry = IsAffineElementTransformation(transformation);
  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetBarycentricGradients(transformation, center, center_jacobian_determinant);
  const auto affine_geometry_tensors = BuildAffineGeometryTensors(center_grad_lambda);
  const mfem::DenseMatrix center_jacobian(transformation.Jacobian());
  const mfem::DenseMatrix center_inverse_jacobian(transformation.InverseJacobian());
  ValidateInputs(element_dofs, center_grad_lambda, center_jacobian_determinant, options);
  const auto h1_to_nd = BuildH1ToNDMap(element_dofs);

  ElementStandardEnrichmentMatrices result;
  const int standard_h1_size = h1_fe.GetDof();
  const int standard_nd_size = nd_fe.GetDof();
  const int enrichment_h1_size = static_cast<int>(element_dofs.h1.size());
  const int enrichment_nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_standard_enrichment.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_enrichment_standard.SetSize(enrichment_h1_size, standard_h1_size);
  result.h1_estimated_absolute_error.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_mass_standard_enrichment.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_mass_enrichment_standard.SetSize(enrichment_h1_size, standard_h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(standard_h1_size, enrichment_h1_size);
  result.nd_mass_standard_enrichment.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_mass_enrichment_standard.SetSize(enrichment_nd_size, standard_nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_curl_curl_standard_enrichment.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_curl_curl_enrichment_standard.SetSize(enrichment_nd_size, standard_nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(standard_nd_size,
                                                       enrichment_nd_size);

  mfem::IntegrationPoint point;
  mfem::DenseMatrix standard_value(standard_nd_size, 3);
  mfem::DenseMatrix reference_standard_value(standard_nd_size, 3);
  mfem::DenseMatrix standard_curl(standard_nd_size, 3);
  mfem::DenseMatrix reference_standard_curl(standard_nd_size, 3);
  const std::vector<const std::vector<AffineStandardReferenceTable::Entry> *>
      *affine_nd_pattern =
          affine_geometry && standard_reference_table && (!common_mass || !common_curl)
              ? &standard_reference_table->GetPattern(nd_fe, element_dofs.nd)
              : nullptr;
  result.setup_time = elapsed(setup_start);
  const auto nd_coupling_start = Clock::now();
  for (int enrichment = 0; enrichment < enrichment_nd_size; enrichment++)
  {
    const auto &basis = element_dofs.nd[enrichment].basis;
    std::optional<AdaptiveVectorQuadratureResult> affine_mass_reintegration;
    std::optional<AdaptiveVectorQuadratureResult> affine_curl_reintegration;
    for (int standard = 0; standard < standard_nd_size; standard++)
    {
      const AffineStandardReferenceTable::Entry *affine_reference =
          affine_nd_pattern ? &(*(*affine_nd_pattern)[enrichment])[standard] : nullptr;
      if (common_mass)
      {
        const double value = common_mass->standard_enrichment(standard, enrichment);
        const double error =
            common_mass->standard_enrichment_estimated_absolute_error(standard, enrichment);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND mass supplied an invalid coupling entry!");
        }
        result.nd_mass_standard_enrichment(standard, enrichment) = value;
        result.nd_mass_enrichment_standard(enrichment, standard) = value;
        result.nd_mass_estimated_absolute_error(standard, enrichment) = error;
        result.affine_nd_mass_contraction_count++;
      }
      else
      {
        const auto integrate_mass = [&]()
        {
          if (affine_reference)
          {
            const double value = ContractAffineTensor(affine_reference->mass.value,
                                                      affine_geometry_tensors.mass,
                                                      center_jacobian_determinant, false);
            const double error = ContractAffineTensor(
                affine_reference->mass.estimated_absolute_error,
                affine_geometry_tensors.mass, center_jacobian_determinant, true);
            if (error <=
                options.absolute_tolerance + options.relative_tolerance * std::abs(value))
            {
              result.affine_nd_mass_contraction_count++;
              return AdaptiveQuadratureResult{
                  value, error, affine_reference->mass.total_leaf_count,
                  affine_reference->mass.maximum_subdivision_depth, true};
            }
            result.affine_nd_mass_reintegration_count++;
            if (!affine_mass_reintegration)
            {
              result.affine_nd_mass_reintegration_batch_count++;
              affine_mass_reintegration = IntegrateConfigured(
                  options.quadrature_order, options,
                  static_cast<std::size_t>(standard_nd_size),
                  [&](const BarycentricPoint &lambda, std::vector<double> &value)
                  {
                    point.Set3(lambda[1], lambda[2], lambda[3]);
                    CalcAffinePhysVShape(nd_fe, point, center_inverse_jacobian,
                                         reference_standard_value, standard_value);
                    const auto singular =
                        EvaluateHigherOrderBasis(lambda, center_grad_lambda, basis);
                    for (int standard_dof = 0; standard_dof < standard_nd_size;
                         standard_dof++)
                    {
                      value[standard_dof] =
                          center_jacobian_determinant *
                          Dot(DenseMatrixRow(standard_value, standard_dof), singular.value);
                    }
                  });
            }
            const double reintegrated_value = affine_mass_reintegration->value[standard];
            const double reintegrated_error =
                affine_mass_reintegration->estimated_absolute_error[standard];
            return AdaptiveQuadratureResult{
                reintegrated_value, reintegrated_error,
                affine_mass_reintegration->leaf_count,
                affine_mass_reintegration->maximum_subdivision_depth,
                std::isfinite(reintegrated_value) && std::isfinite(reintegrated_error) &&
                    reintegrated_error >= 0.0 &&
                    reintegrated_error <=
                        options.absolute_tolerance +
                            options.relative_tolerance * std::abs(reintegrated_value)};
          }
          return IntegrateConfigured(
              options.quadrature_order, options,
              [&](const BarycentricPoint &lambda)
              {
                point.Set3(lambda[1], lambda[2], lambda[3]);
                double jacobian_determinant = center_jacobian_determinant;
                auto grad_lambda = center_grad_lambda;
                if (affine_geometry)
                {
                  CalcAffinePhysVShape(nd_fe, point, center_inverse_jacobian,
                                       reference_standard_value, standard_value);
                }
                else
                {
                  grad_lambda =
                      GetBarycentricGradients(transformation, point, jacobian_determinant);
                  nd_fe.CalcPhysVShape(transformation, standard_value);
                }
                const auto singular = EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
                return jacobian_determinant *
                       Dot(DenseMatrixRow(standard_value, standard), singular.value);
              });
        };
        const auto mass = integrate_mass();
        ValidateAdaptiveEntry(mass, "mass", standard, enrichment, basis);
        RecordQuadratureStatistics(mass, result);
        result.nd_mass_standard_enrichment(standard, enrichment) = mass.value;
        result.nd_mass_enrichment_standard(enrichment, standard) = mass.value;
        result.nd_mass_estimated_absolute_error(standard, enrichment) =
            mass.estimated_absolute_error;
      }

      if (IsGradientFamily(basis.family))
      {
        result.nd_curl_curl_standard_enrichment(standard, enrichment) = 0.0;
        result.nd_curl_curl_enrichment_standard(enrichment, standard) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(standard, enrichment) = 0.0;
        continue;
      }

      if (common_curl)
      {
        const double value = common_curl->standard_enrichment(standard, enrichment);
        const double error =
            common_curl->standard_enrichment_estimated_absolute_error(standard, enrichment);
        if (!std::isfinite(value) || !std::isfinite(error) || error < 0.0)
        {
          throw std::runtime_error(
              "Common affine ND curl supplied an invalid coupling entry!");
        }
        result.nd_curl_curl_standard_enrichment(standard, enrichment) = value;
        result.nd_curl_curl_enrichment_standard(enrichment, standard) = value;
        result.nd_curl_curl_estimated_absolute_error(standard, enrichment) = error;
        result.affine_nd_curl_contraction_count++;
        continue;
      }

      const auto integrate_curl_curl = [&]()
      {
        if (affine_reference)
        {
          const double value = ContractAffineTensor(affine_reference->curl_curl.value,
                                                    affine_geometry_tensors.curl_curl,
                                                    center_jacobian_determinant, false);
          const double error = ContractAffineTensor(
              affine_reference->curl_curl.estimated_absolute_error,
              affine_geometry_tensors.curl_curl, center_jacobian_determinant, true);
          if (error <=
              options.absolute_tolerance + options.relative_tolerance * std::abs(value))
          {
            result.affine_nd_curl_contraction_count++;
            return AdaptiveQuadratureResult{
                value, error, affine_reference->curl_curl.total_leaf_count,
                affine_reference->curl_curl.maximum_subdivision_depth, true};
          }
          result.affine_nd_curl_reintegration_count++;
          if (!affine_curl_reintegration)
          {
            result.affine_nd_curl_reintegration_batch_count++;
            affine_curl_reintegration = IntegrateConfigured(
                options.quadrature_order, options,
                static_cast<std::size_t>(standard_nd_size),
                [&](const BarycentricPoint &lambda, std::vector<double> &value)
                {
                  point.Set3(lambda[1], lambda[2], lambda[3]);
                  CalcAffinePhysCurlShape(nd_fe, point, center_jacobian,
                                          center_jacobian_determinant,
                                          reference_standard_curl, standard_curl);
                  const auto singular =
                      EvaluateHigherOrderBasis(lambda, center_grad_lambda, basis);
                  for (int standard_dof = 0; standard_dof < standard_nd_size;
                       standard_dof++)
                  {
                    value[standard_dof] =
                        center_jacobian_determinant *
                        Dot(DenseMatrixRow(standard_curl, standard_dof), singular.curl);
                  }
                });
          }
          const double reintegrated_value = affine_curl_reintegration->value[standard];
          const double reintegrated_error =
              affine_curl_reintegration->estimated_absolute_error[standard];
          return AdaptiveQuadratureResult{
              reintegrated_value, reintegrated_error, affine_curl_reintegration->leaf_count,
              affine_curl_reintegration->maximum_subdivision_depth,
              std::isfinite(reintegrated_value) && std::isfinite(reintegrated_error) &&
                  reintegrated_error >= 0.0 &&
                  reintegrated_error <=
                      options.absolute_tolerance +
                          options.relative_tolerance * std::abs(reintegrated_value)};
        }
        return IntegrateConfigured(
            options.quadrature_order, options,
            [&](const BarycentricPoint &lambda)
            {
              point.Set3(lambda[1], lambda[2], lambda[3]);
              double jacobian_determinant = center_jacobian_determinant;
              auto grad_lambda = center_grad_lambda;
              if (affine_geometry)
              {
                CalcAffinePhysCurlShape(nd_fe, point, center_jacobian,
                                        center_jacobian_determinant,
                                        reference_standard_curl, standard_curl);
              }
              else
              {
                grad_lambda =
                    GetBarycentricGradients(transformation, point, jacobian_determinant);
                nd_fe.CalcPhysCurlShape(transformation, standard_curl);
              }
              const auto singular = EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
              return jacobian_determinant *
                     Dot(DenseMatrixRow(standard_curl, standard), singular.curl);
            });
      };
      const auto curl_curl = integrate_curl_curl();
      ValidateAdaptiveEntry(curl_curl, "curl-curl", standard, enrichment, basis);
      RecordQuadratureStatistics(curl_curl, result);
      result.nd_curl_curl_standard_enrichment(standard, enrichment) = curl_curl.value;
      result.nd_curl_curl_enrichment_standard(enrichment, standard) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(standard, enrichment) =
          curl_curl.estimated_absolute_error;
    }
  }

  result.nd_coupling_time = elapsed(nd_coupling_start);
  const auto h1_gradient_coupling_start = Clock::now();
  mfem::DenseMatrix discrete_gradient(standard_nd_size, standard_h1_size);
  nd_fe.ProjectGrad(h1_fe, transformation, discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      long double value = 0.0L;
      long double error = 0.0L;
      const int enrichment_nd = h1_to_nd[enrichment_h1];
      for (int standard_nd = 0; standard_nd < standard_nd_size; standard_nd++)
      {
        const double coefficient = discrete_gradient(standard_nd, standard_h1);
        if (!std::isfinite(coefficient))
        {
          throw std::runtime_error(
              "Singular element discrete gradient contains a nonfinite entry!");
        }
        value += static_cast<long double>(coefficient) *
                 static_cast<long double>(
                     result.nd_mass_standard_enrichment(standard_nd, enrichment_nd));
        error += std::abs(static_cast<long double>(coefficient)) *
                 static_cast<long double>(
                     result.nd_mass_estimated_absolute_error(standard_nd, enrichment_nd));
      }
      const double entry = static_cast<double>(value);
      const double entry_error = static_cast<double>(error);
      if (!std::isfinite(entry) || !std::isfinite(entry_error) || entry_error < 0.0)
      {
        throw std::runtime_error(
            "Singular element assembly produced a nonfinite H1 coupling entry!");
      }
      result.h1_standard_enrichment(standard_h1, enrichment_h1) = entry;
      result.h1_enrichment_standard(enrichment_h1, standard_h1) = entry;
      result.h1_estimated_absolute_error(standard_h1, enrichment_h1) = entry_error;
    }
  }

  result.h1_gradient_coupling_time = elapsed(h1_gradient_coupling_start);
  const auto h1_mass_coupling_start = Clock::now();
  mfem::Vector standard_shape(standard_h1_size);
  const std::vector<const AdaptiveQuadratureResult *> *affine_h1_pattern =
      affine_geometry && scalar_reference_table
          ? &scalar_reference_table->GetStandardPattern(h1_fe, element_dofs.h1)
          : nullptr;
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      const auto &basis = element_dofs.h1[enrichment_h1].basis;
      const auto integrate_h1_mass = [&]()
      {
        if (affine_h1_pattern)
        {
          auto reference = *(*affine_h1_pattern)[static_cast<std::size_t>(standard_h1) *
                                                     element_dofs.h1.size() +
                                                 enrichment_h1];
          reference.value *= center_jacobian_determinant;
          reference.estimated_absolute_error *= center_jacobian_determinant;
          if (reference.estimated_absolute_error <=
              options.absolute_tolerance +
                  options.relative_tolerance * std::abs(reference.value))
          {
            return reference;
          }
        }
        return IntegrateConfigured(
            options.quadrature_order, options,
            [&](const BarycentricPoint &lambda)
            {
              point.Set3(lambda[1], lambda[2], lambda[3]);
              double jacobian_determinant = center_jacobian_determinant;
              if (!affine_geometry)
              {
                (void)GetBarycentricGradients(transformation, point, jacobian_determinant);
              }
              h1_fe.CalcShape(point, standard_shape);
              return jacobian_determinant * standard_shape[standard_h1] *
                     EvaluateHigherOrderGradientPotential(lambda, basis);
            });
      };
      const auto mass = integrate_h1_mass();
      ValidateAdaptiveEntry(mass, "H1 mass", standard_h1, enrichment_h1, basis);
      RecordQuadratureStatistics(mass, result);
      result.h1_mass_standard_enrichment(standard_h1, enrichment_h1) = mass.value;
      result.h1_mass_enrichment_standard(enrichment_h1, standard_h1) = mass.value;
      result.h1_mass_estimated_absolute_error(standard_h1, enrichment_h1) =
          mass.estimated_absolute_error;
    }
  }
  result.h1_mass_coupling_time = elapsed(h1_mass_coupling_start);
  return result;
}

}  // namespace

ElementStandardEnrichmentMatrices AssembleElementStandardEnrichmentMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options)
{
  AffineStandardReferenceTable standard_reference_table(options);
  AffineScalarReferenceTable scalar_reference_table(options);
  return AssembleElementStandardEnrichmentMatricesImpl(
      element_dofs, h1_fe, nd_fe, transformation, options, &standard_reference_table,
      &scalar_reference_table, nullptr, nullptr);
}

namespace
{

ElementStandardEnrichmentMatrices AssembleElementStandardEnrichmentMatricesCached(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options,
    AffineStandardReferenceTable &standard_reference_table,
    AffineScalarReferenceTable &scalar_reference_table,
    const AffineCommonNDMassBlocks *common_mass,
    const AffineCommonNDCurlBlocks *common_curl)
{
  return AssembleElementStandardEnrichmentMatricesImpl(
      element_dofs, h1_fe, nd_fe, transformation, options, &standard_reference_table,
      &scalar_reference_table, common_mass, common_curl);
}

}  // namespace

ElementStandardEnrichmentMatrices AssembleTriangleElementStandardEnrichmentMatrices(
    const TriangleElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options)
{
  ValidateStandardTriangleFiniteElements(h1_fe, nd_fe, transformation);
  const bool affine_geometry = IsAffineElementTransformation(transformation);
  mfem::IntegrationPoint center;
  center.Set2(1.0 / 3.0, 1.0 / 3.0);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetTriangleBarycentricGradients(transformation, center, center_jacobian_determinant);
  const mfem::DenseMatrix center_jacobian(transformation.Jacobian());
  const mfem::DenseMatrix center_inverse_jacobian(transformation.InverseJacobian());
  ValidateTriangleInputs(element_dofs, center_grad_lambda, center_jacobian_determinant,
                         options);
  const auto h1_to_nd = BuildTriangleH1ToNDMap(element_dofs);
  TriangleBarycentricGradients normalized_center_grad_lambda;
  const double gradient_scale =
      affine_geometry ? NormalizeTriangleBarycentricGradients(center_grad_lambda,
                                                              normalized_center_grad_lambda)
                      : 1.0;
  if (!affine_geometry)
  {
    normalized_center_grad_lambda = center_grad_lambda;
  }
  const double affine_mass_prefactor =
      ScaleTriangleJacobianDeterminant(center_jacobian_determinant, gradient_scale, 2);
  const double affine_curl_curl_prefactor =
      ScaleTriangleJacobianDeterminant(center_jacobian_determinant, gradient_scale, 4);
  mfem::DenseMatrix normalized_center_inverse_jacobian(center_inverse_jacobian);
  normalized_center_inverse_jacobian *= 1.0 / gradient_scale;

  ElementStandardEnrichmentMatrices result;
  const int standard_h1_size = h1_fe.GetDof();
  const int standard_nd_size = nd_fe.GetDof();
  const int enrichment_h1_size = static_cast<int>(element_dofs.h1.size());
  const int enrichment_nd_size = static_cast<int>(element_dofs.nd.size());
  result.h1_standard_enrichment.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_enrichment_standard.SetSize(enrichment_h1_size, standard_h1_size);
  result.h1_estimated_absolute_error.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_mass_standard_enrichment.SetSize(standard_h1_size, enrichment_h1_size);
  result.h1_mass_enrichment_standard.SetSize(enrichment_h1_size, standard_h1_size);
  result.h1_mass_estimated_absolute_error.SetSize(standard_h1_size, enrichment_h1_size);
  result.nd_mass_standard_enrichment.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_mass_enrichment_standard.SetSize(enrichment_nd_size, standard_nd_size);
  result.nd_mass_estimated_absolute_error.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_curl_curl_standard_enrichment.SetSize(standard_nd_size, enrichment_nd_size);
  result.nd_curl_curl_enrichment_standard.SetSize(enrichment_nd_size, standard_nd_size);
  result.nd_curl_curl_estimated_absolute_error.SetSize(standard_nd_size,
                                                       enrichment_nd_size);

  mfem::IntegrationPoint point;
  mfem::DenseMatrix standard_value(standard_nd_size, 2);
  mfem::DenseMatrix reference_standard_value(standard_nd_size, 2);
  mfem::DenseMatrix standard_curl(standard_nd_size, 1);
  mfem::DenseMatrix reference_standard_curl(standard_nd_size, 1);
  for (int standard = 0; standard < standard_nd_size; standard++)
  {
    for (int enrichment = 0; enrichment < enrichment_nd_size; enrichment++)
    {
      const auto &basis = element_dofs.nd[enrichment].basis;
      const int singular_node = TriangleSingularNode(basis);
      const auto mass = IntegrateTriangleDuffy(
          singular_node, singular_node, "standard-enrichment ND mass", standard, enrichment,
          options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double integration_prefactor = center_jacobian_determinant;
            auto grad_lambda = normalized_center_grad_lambda;
            if (affine_geometry)
            {
              integration_prefactor = affine_mass_prefactor;
              CalcAffinePhysVShape(nd_fe, point, normalized_center_inverse_jacobian,
                                   reference_standard_value, standard_value);
            }
            else
            {
              grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                            integration_prefactor);
              nd_fe.CalcPhysVShape(transformation, standard_value);
            }
            const auto singular = affine_geometry
                                      ? EvaluateAffineNormalizedTriangleBasis(
                                            lambda, normalized_center_grad_lambda,
                                            affine_mass_prefactor, basis)
                                      : EvaluateTriangleBasis(lambda, grad_lambda, basis);
            return integration_prefactor *
                   Dot(DenseMatrixRow2(standard_value, standard), singular.value);
          });
      result.nd_mass_standard_enrichment(standard, enrichment) = mass.value;
      result.nd_mass_enrichment_standard(enrichment, standard) = mass.value;
      result.nd_mass_estimated_absolute_error(standard, enrichment) =
          mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;

      if (IsGradientFamily(basis.family))
      {
        result.nd_curl_curl_standard_enrichment(standard, enrichment) = 0.0;
        result.nd_curl_curl_enrichment_standard(enrichment, standard) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(standard, enrichment) = 0.0;
        continue;
      }
      const auto curl_curl = IntegrateTriangleDuffy(
          singular_node, singular_node, "standard-enrichment ND curl-curl", standard,
          enrichment, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double integration_prefactor = center_jacobian_determinant;
            auto grad_lambda = normalized_center_grad_lambda;
            if (affine_geometry)
            {
              integration_prefactor = affine_curl_curl_prefactor;
              CalcAffinePhysCurlShape(nd_fe, point, center_jacobian, affine_mass_prefactor,
                                      reference_standard_curl, standard_curl);
            }
            else
            {
              grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                            integration_prefactor);
              nd_fe.CalcPhysCurlShape(transformation, standard_curl);
            }
            const auto singular = affine_geometry
                                      ? EvaluateAffineNormalizedTriangleBasis(
                                            lambda, normalized_center_grad_lambda,
                                            affine_mass_prefactor, basis)
                                      : EvaluateTriangleBasis(lambda, grad_lambda, basis);
            return integration_prefactor * standard_curl(standard, 0) * singular.curl;
          });
      result.nd_curl_curl_standard_enrichment(standard, enrichment) = curl_curl.value;
      result.nd_curl_curl_enrichment_standard(enrichment, standard) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(standard, enrichment) =
          curl_curl.estimated_absolute_error;
      result.total_quadrature_leaf_count += curl_curl.point_count;
    }
  }

  mfem::DenseMatrix discrete_gradient(standard_nd_size, standard_h1_size);
  nd_fe.ProjectGrad(h1_fe, transformation, discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      long double value = 0.0L;
      long double error = 0.0L;
      const int enrichment_nd = h1_to_nd[enrichment_h1];
      for (int standard_nd = 0; standard_nd < standard_nd_size; standard_nd++)
      {
        const double coefficient = discrete_gradient(standard_nd, standard_h1);
        if (!std::isfinite(coefficient))
        {
          throw std::runtime_error(
              "Triangular singular discrete gradient contains a nonfinite entry!");
        }
        value += static_cast<long double>(coefficient) *
                 result.nd_mass_standard_enrichment(standard_nd, enrichment_nd);
        error += std::abs(static_cast<long double>(coefficient)) *
                 result.nd_mass_estimated_absolute_error(standard_nd, enrichment_nd);
      }
      const double entry = static_cast<double>(value);
      const double entry_error = static_cast<double>(error);
      if (!std::isfinite(entry) || !std::isfinite(entry_error) || entry_error < 0.0)
      {
        throw std::runtime_error(
            "Triangular singular assembly produced an invalid H1 coupling entry!");
      }
      result.h1_standard_enrichment(standard_h1, enrichment_h1) = entry;
      result.h1_enrichment_standard(enrichment_h1, standard_h1) = entry;
      result.h1_estimated_absolute_error(standard_h1, enrichment_h1) = entry_error;
    }
  }

  mfem::Vector standard_shape(standard_h1_size);
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      const auto &basis = element_dofs.h1[enrichment_h1].basis;
      const int singular_node = TriangleSingularNode(basis);
      const auto mass = IntegrateTriangleDuffy(
          singular_node, singular_node, "standard-enrichment H1 mass", standard_h1,
          enrichment_h1, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant = center_jacobian_determinant;
            if (!affine_geometry)
            {
              (void)GetTriangleBarycentricGradients(transformation, point,
                                                    jacobian_determinant);
            }
            h1_fe.CalcShape(point, standard_shape);
            return jacobian_determinant * standard_shape[standard_h1] *
                   EvaluateTriangleNodeGradientPotential(lambda, basis.nodes[0],
                                                         basis.nodes[1], basis.nu);
          });
      result.h1_mass_standard_enrichment(standard_h1, enrichment_h1) = mass.value;
      result.h1_mass_enrichment_standard(enrichment_h1, standard_h1) = mass.value;
      result.h1_mass_estimated_absolute_error(standard_h1, enrichment_h1) =
          mass.estimated_absolute_error;
      result.total_quadrature_leaf_count += mass.point_count;
    }
  }
  return result;
}

namespace
{

ElementH1EnrichmentMatrices AssembleElementH1EnrichmentMatricesImpl(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    mfem::ElementTransformation &transformation, const AdaptiveAssemblyOptions &options,
    DuffyH1ReferenceTable *duffy_table,
    AffineCommonNDMassReferenceTable *affine_reference_table,
    const mfem::FiniteElement *nd_fe)
{
  if (h1_fe.GetDim() != 3 || h1_fe.GetGeomType() != mfem::Geometry::TETRAHEDRON ||
      h1_fe.GetRangeType() != mfem::FiniteElement::SCALAR ||
      h1_fe.GetMapType() != mfem::FiniteElement::VALUE ||
      h1_fe.GetDerivType() != mfem::FiniteElement::GRAD || h1_fe.GetOrder() < 1 ||
      transformation.GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
      transformation.GetDimension() != 3 || transformation.GetSpaceDim() != 3)
  {
    throw std::invalid_argument(
        "Singular H1 coupling requires a three-dimensional tetrahedral finite element "
        "of positive order!");
  }

  const bool affine_geometry = IsAffineElementTransformation(transformation);
  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetBarycentricGradients(transformation, center, center_jacobian_determinant);
  const mfem::DenseMatrix center_inverse_jacobian(transformation.InverseJacobian());
  ValidateInputs(element_dofs, center_grad_lambda, center_jacobian_determinant, options);

  ElementH1EnrichmentMatrices result;
  const int standard_size = h1_fe.GetDof();
  const int enrichment_size = static_cast<int>(element_dofs.h1.size());
  result.enrichment_enrichment.SetSize(enrichment_size);
  result.enrichment_enrichment_estimated_absolute_error.SetSize(enrichment_size);
  result.standard_enrichment.SetSize(standard_size, enrichment_size);
  result.enrichment_standard.SetSize(enrichment_size, standard_size);
  result.standard_enrichment_estimated_absolute_error.SetSize(standard_size,
                                                              enrichment_size);

  if (affine_geometry && affine_reference_table && nd_fe)
  {
    const auto common = affine_reference_table->Get(
        *nd_fe, element_dofs.h1, center_grad_lambda, center_jacobian_determinant);
    const int standard_nd_size = nd_fe->GetDof();
    if (common.standard_enrichment.Height() != standard_nd_size ||
        common.standard_enrichment.Width() != enrichment_size ||
        common.enrichment_enrichment.Height() != enrichment_size ||
        common.enrichment_enrichment.Width() != enrichment_size ||
        common.standard_enrichment_estimated_absolute_error.Height() != standard_nd_size ||
        common.standard_enrichment_estimated_absolute_error.Width() != enrichment_size ||
        common.enrichment_enrichment_estimated_absolute_error.Height() != enrichment_size ||
        common.enrichment_enrichment_estimated_absolute_error.Width() != enrichment_size)
    {
      throw std::runtime_error(
          "Cached affine singular H1 reference blocks have inconsistent dimensions!");
    }

    result.enrichment_enrichment = common.enrichment_enrichment;
    result.enrichment_enrichment_estimated_absolute_error =
        common.enrichment_enrichment_estimated_absolute_error;
    mfem::DenseMatrix discrete_gradient(standard_nd_size, standard_size);
    nd_fe->ProjectGrad(h1_fe, transformation, discrete_gradient);
    for (int standard = 0; standard < standard_size; standard++)
    {
      for (int enrichment = 0; enrichment < enrichment_size; enrichment++)
      {
        long double value = 0.0L;
        long double error = 0.0L;
        for (int standard_nd = 0; standard_nd < standard_nd_size; standard_nd++)
        {
          const double coefficient = discrete_gradient(standard_nd, standard);
          if (!std::isfinite(coefficient))
          {
            throw std::runtime_error(
                "Cached affine singular H1 discrete gradient contains a nonfinite "
                "entry!");
          }
          value += static_cast<long double>(coefficient) *
                   common.standard_enrichment(standard_nd, enrichment);
          error +=
              std::abs(static_cast<long double>(coefficient)) *
              common.standard_enrichment_estimated_absolute_error(standard_nd, enrichment);
        }
        const double entry = static_cast<double>(value);
        const double entry_error = static_cast<double>(error);
        if (!std::isfinite(entry) || !std::isfinite(entry_error) || entry_error < 0.0)
        {
          throw std::runtime_error(
              "Cached affine singular H1 contraction produced an invalid entry!");
        }
        result.standard_enrichment(standard, enrichment) = entry;
        result.enrichment_standard(enrichment, standard) = entry;
        result.standard_enrichment_estimated_absolute_error(standard, enrichment) =
            entry_error;
      }
    }
    result.total_quadrature_leaf_count = common.leaf_count;
    result.maximum_subdivision_depth = common.maximum_subdivision_depth;
    return result;
  }

  for (int row = 0; row < enrichment_size; row++)
  {
    for (int column = row; column < enrichment_size; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      if (duffy_table && affine_geometry)
      {
        const auto reference =
            duffy_table->Get(ReferenceBasis{row_basis}, ReferenceBasis{column_basis});
        if (reference)
        {
          const double value = ContractMass(reference->integral, center_grad_lambda,
                                            center_jacobian_determinant);
          const double error =
              ContractMassError(reference->integral, reference->estimated_absolute_error,
                                center_grad_lambda, center_jacobian_determinant);
          ValidateDuffyH1Entry(value, error, "enrichment-enrichment", row, column,
                               column_basis, options);
          result.enrichment_enrichment(row, column) =
              result.enrichment_enrichment(column, row) = value;
          result.enrichment_enrichment_estimated_absolute_error(row, column) =
              result.enrichment_enrichment_estimated_absolute_error(column, row) = error;
          continue;
        }
      }
      const auto integral = IntegrateConfigured(
          options.quadrature_order, options,
          [&](const BarycentricPoint &lambda)
          {
            mfem::IntegrationPoint point;
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant = center_jacobian_determinant;
            auto grad_lambda = center_grad_lambda;
            if (!affine_geometry)
            {
              grad_lambda =
                  GetBarycentricGradients(transformation, point, jacobian_determinant);
            }
            const auto row_value = EvaluateHigherOrderBasis(lambda, grad_lambda, row_basis);
            const auto column_value =
                EvaluateHigherOrderBasis(lambda, grad_lambda, column_basis);
            return jacobian_determinant * Dot(row_value.value, column_value.value);
          });
      ValidateAdaptiveH1Entry(integral, "enrichment-enrichment", row, column, column_basis);
      RecordQuadratureStatistics(integral, result);
      result.enrichment_enrichment(row, column) =
          result.enrichment_enrichment(column, row) = integral.value;
      result.enrichment_enrichment_estimated_absolute_error(row, column) =
          result.enrichment_enrichment_estimated_absolute_error(column, row) =
              integral.estimated_absolute_error;
    }
  }

  mfem::IntegrationPoint point;
  mfem::DenseMatrix standard_gradient(standard_size, 3);
  mfem::DenseMatrix reference_standard_gradient(standard_size, 3);
  for (int standard = 0; standard < standard_size; standard++)
  {
    for (int enrichment = 0; enrichment < enrichment_size; enrichment++)
    {
      const auto &basis = element_dofs.h1[enrichment].basis;
      if (duffy_table && affine_geometry)
      {
        if (standard_size != 4)
        {
          throw std::logic_error(
              "First-order singular H1 reference coupling requires four standard DOFs!");
        }
        const auto reference = duffy_table->Get(
            ReferenceBasis{MakeStandardH1Gradient(standard)}, ReferenceBasis{basis});
        if (!reference)
        {
          throw std::logic_error(
              "A standard H1 gradient and one singular gradient must share a Duffy "
              "feature!");
        }
        const double value = ContractMass(reference->integral, center_grad_lambda,
                                          center_jacobian_determinant);
        const double error =
            ContractMassError(reference->integral, reference->estimated_absolute_error,
                              center_grad_lambda, center_jacobian_determinant);
        ValidateDuffyH1Entry(value, error, "standard-enrichment", standard, enrichment,
                             basis, options);
        result.standard_enrichment(standard, enrichment) = value;
        result.enrichment_standard(enrichment, standard) = value;
        result.standard_enrichment_estimated_absolute_error(standard, enrichment) = error;
        continue;
      }
      const auto integral = IntegrateConfigured(
          options.quadrature_order, options,
          [&](const BarycentricPoint &lambda)
          {
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant = center_jacobian_determinant;
            auto grad_lambda = center_grad_lambda;
            if (affine_geometry)
            {
              CalcAffinePhysDShape(h1_fe, point, center_inverse_jacobian,
                                   reference_standard_gradient, standard_gradient);
            }
            else
            {
              grad_lambda =
                  GetBarycentricGradients(transformation, point, jacobian_determinant);
              h1_fe.CalcPhysDShape(transformation, standard_gradient);
            }
            const auto singular = EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
            return jacobian_determinant *
                   Dot(DenseMatrixRow(standard_gradient, standard), singular.value);
          });
      ValidateAdaptiveH1Entry(integral, "standard-enrichment", standard, enrichment, basis);
      RecordQuadratureStatistics(integral, result);
      result.standard_enrichment(standard, enrichment) = integral.value;
      result.enrichment_standard(enrichment, standard) = integral.value;
      result.standard_enrichment_estimated_absolute_error(standard, enrichment) =
          integral.estimated_absolute_error;
    }
  }
  return result;
}

}  // namespace

ElementH1EnrichmentMatrices AssembleElementH1EnrichmentMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    mfem::ElementTransformation &transformation, const AdaptiveAssemblyOptions &options)
{
  DuffyH1ReferenceTable duffy_table;
  return AssembleElementH1EnrichmentMatricesImpl(
      element_dofs, h1_fe, transformation, options,
      (!options.fixed_subdivision && h1_fe.GetOrder() == 1) ? &duffy_table : nullptr,
      nullptr, nullptr);
}

void ApplyStandardDofTransformations(const mfem::DofTransformation &h1_dof_transformation,
                                     const mfem::DofTransformation &nd_dof_transformation,
                                     ElementStandardEnrichmentMatrices &matrices)
{
  ValidateStandardRowTransformation(h1_dof_transformation, matrices.h1_standard_enrichment,
                                    matrices.h1_enrichment_standard,
                                    matrices.h1_estimated_absolute_error);
  ValidateStandardRowTransformation(
      nd_dof_transformation, matrices.nd_mass_standard_enrichment,
      matrices.nd_mass_enrichment_standard, matrices.nd_mass_estimated_absolute_error);
  ValidateStandardRowTransformation(
      h1_dof_transformation, matrices.h1_mass_standard_enrichment,
      matrices.h1_mass_enrichment_standard, matrices.h1_mass_estimated_absolute_error);
  ValidateStandardRowTransformation(nd_dof_transformation,
                                    matrices.nd_curl_curl_standard_enrichment,
                                    matrices.nd_curl_curl_enrichment_standard,
                                    matrices.nd_curl_curl_estimated_absolute_error);

  ApplyStandardRowTransformation(h1_dof_transformation, matrices.h1_standard_enrichment,
                                 matrices.h1_enrichment_standard,
                                 matrices.h1_estimated_absolute_error);
  ApplyStandardRowTransformation(
      h1_dof_transformation, matrices.h1_mass_standard_enrichment,
      matrices.h1_mass_enrichment_standard, matrices.h1_mass_estimated_absolute_error);
  ApplyStandardRowTransformation(
      nd_dof_transformation, matrices.nd_mass_standard_enrichment,
      matrices.nd_mass_enrichment_standard, matrices.nd_mass_estimated_absolute_error);
  ApplyStandardRowTransformation(nd_dof_transformation,
                                 matrices.nd_curl_curl_standard_enrichment,
                                 matrices.nd_curl_curl_enrichment_standard,
                                 matrices.nd_curl_curl_estimated_absolute_error);
}

void ApplyIsotropicMaterialCoefficients(const IsotropicMaterialCoefficients &coefficients,
                                        ElementEnrichmentMatrices &matrices)
{
  ValidateMaterialCoefficients(coefficients);
  ValidateSymmetricMatrixForScaling(matrices.h1_diffusion,
                                    matrices.h1_diffusion_estimated_absolute_error,
                                    coefficients.electric);
  ValidateSymmetricMatrixForScaling(
      matrices.h1_mass, matrices.h1_mass_estimated_absolute_error, coefficients.electric);
  ValidateSymmetricMatrixForScaling(
      matrices.nd_mass, matrices.nd_mass_estimated_absolute_error, coefficients.electric);
  ValidateSymmetricMatrixForScaling(matrices.nd_curl_curl,
                                    matrices.nd_curl_curl_estimated_absolute_error,
                                    coefficients.inverse_magnetic);

  ScaleSymmetricMatrices(matrices.h1_diffusion,
                         matrices.h1_diffusion_estimated_absolute_error,
                         coefficients.electric);
  ScaleSymmetricMatrices(matrices.h1_mass, matrices.h1_mass_estimated_absolute_error,
                         coefficients.electric);
  ScaleSymmetricMatrices(matrices.nd_mass, matrices.nd_mass_estimated_absolute_error,
                         coefficients.electric);
  ScaleSymmetricMatrices(matrices.nd_curl_curl,
                         matrices.nd_curl_curl_estimated_absolute_error,
                         coefficients.inverse_magnetic);
}

void ApplyIsotropicMaterialCoefficients(const IsotropicMaterialCoefficients &coefficients,
                                        ElementStandardEnrichmentMatrices &matrices)
{
  ValidateMaterialCoefficients(coefficients);
  ValidateCouplingMatricesForScaling(
      matrices.h1_standard_enrichment, matrices.h1_enrichment_standard,
      matrices.h1_estimated_absolute_error, coefficients.electric);
  ValidateCouplingMatricesForScaling(
      matrices.h1_mass_standard_enrichment, matrices.h1_mass_enrichment_standard,
      matrices.h1_mass_estimated_absolute_error, coefficients.electric);
  ValidateCouplingMatricesForScaling(
      matrices.nd_mass_standard_enrichment, matrices.nd_mass_enrichment_standard,
      matrices.nd_mass_estimated_absolute_error, coefficients.electric);
  ValidateCouplingMatricesForScaling(
      matrices.nd_curl_curl_standard_enrichment, matrices.nd_curl_curl_enrichment_standard,
      matrices.nd_curl_curl_estimated_absolute_error, coefficients.inverse_magnetic);

  ScaleCouplingMatrices(matrices.h1_standard_enrichment, matrices.h1_enrichment_standard,
                        matrices.h1_estimated_absolute_error, coefficients.electric);
  ScaleCouplingMatrices(matrices.h1_mass_standard_enrichment,
                        matrices.h1_mass_enrichment_standard,
                        matrices.h1_mass_estimated_absolute_error, coefficients.electric);
  ScaleCouplingMatrices(matrices.nd_mass_standard_enrichment,
                        matrices.nd_mass_enrichment_standard,
                        matrices.nd_mass_estimated_absolute_error, coefficients.electric);
  ScaleCouplingMatrices(
      matrices.nd_curl_curl_standard_enrichment, matrices.nd_curl_curl_enrichment_standard,
      matrices.nd_curl_curl_estimated_absolute_error, coefficients.inverse_magnetic);
}

namespace
{

struct ElementNDBoundaryMassMatrices
{
  mfem::DenseMatrix standard_enrichment;
  mfem::DenseMatrix enrichment_standard;
  mfem::DenseMatrix enrichment_enrichment;
  mfem::DenseMatrix standard_enrichment_estimated_absolute_error;
  mfem::DenseMatrix enrichment_enrichment_estimated_absolute_error;
};

struct WeightedSegmentIntegral
{
  std::vector<double> value;
  std::vector<double> estimated_absolute_error;
};

struct TetrahedronFaceTracePowers
{
  std::array<double, 3> node{};
  std::array<std::array<double, 3>, 3> edge{};
};

struct TetrahedronBoundaryTraceValues
{
  std::vector<Vector3> standard;
  std::vector<Vector3> enrichment;
  Vector3 normal;
  double scale;
};

int SetBoundaryIntegrationPoint(mfem::Mesh &mesh, int boundary,
                                const mfem::IntegrationPoint &boundary_point,
                                mfem::FaceElementTransformations &face,
                                mfem::IsoparametricTransformation &element1,
                                mfem::IsoparametricTransformation &element2)
{
  int face_index = -1, orientation = -1;
  mesh.GetBdrElementFace(boundary, &face_index, &orientation);
  if (auto *parallel_mesh = dynamic_cast<mfem::ParMesh *>(&mesh))
  {
    BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
        boundary, *parallel_mesh, face, element1, element2);
  }
  else
  {
    mesh.GetFaceElementTransformations(face_index, face, element1, element2);
  }
  if (!face.Elem1 || face.Elem1No < 0 || face.GetGeometryType() == mfem::Geometry::INVALID)
  {
    throw std::domain_error(
        "Singular lumped-port traces require a conforming boundary element with a "
        "locally owned first adjacent element!");
  }
  const auto face_point = mfem::Mesh::TransformBdrElementToFace(
      face.GetGeometryType(), orientation, boundary_point);
  face.SetAllIntPoints(&face_point);
  return face.Elem1No;
}

int GetTetrahedronVertex(const mfem::IntegrationPoint &point)
{
  const std::array<double, 4> lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                     point.z};
  const auto maximum = std::max_element(lambda.begin(), lambda.end());
  if (*maximum < 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::runtime_error(
        "Singular boundary trace could not identify a tetrahedron vertex!");
  }
  return static_cast<int>(std::distance(lambda.begin(), maximum));
}

std::array<int, 3> GetTetrahedronBoundaryNodes(mfem::Mesh &mesh, int boundary,
                                               mfem::FaceElementTransformations &face,
                                               mfem::IsoparametricTransformation &element1,
                                               mfem::IsoparametricTransformation &element2,
                                               int &element)
{
  std::array<mfem::IntegrationPoint, 3> points;
  points[0].Set2(0.0, 0.0);
  points[1].Set2(1.0, 0.0);
  points[2].Set2(0.0, 1.0);
  std::array<int, 3> nodes{-1, -1, -1};
  element = -1;
  for (int vertex = 0; vertex < 3; vertex++)
  {
    const int current = SetBoundaryIntegrationPoint(mesh, boundary, points[vertex], face,
                                                    element1, element2);
    if (element >= 0 && current != element)
    {
      throw std::runtime_error(
          "Singular boundary trace changed adjacent element across one face!");
    }
    element = current;
    nodes[vertex] = GetTetrahedronVertex(face.Elem1->GetIntPoint());
  }
  auto sorted = nodes;
  std::sort(sorted.begin(), sorted.end());
  if (sorted[0] < 0 || sorted[2] > 3 ||
      std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
  {
    throw std::runtime_error(
        "Singular boundary trace could not map a complete tetrahedron face!");
  }
  return nodes;
}

bool HasTetrahedronBoundaryTrace(const HigherOrderBasis &basis,
                                 const std::array<int, 3> &nodes)
{
  for (int i = 0; i < 4; i++)
  {
    if (basis.interpolation_indices[i] > 0 &&
        std::find(nodes.begin(), nodes.end(), basis.nodes[i]) == nodes.end())
    {
      return false;
    }
  }
  return true;
}

TetrahedronFaceTracePowers GetTetrahedronBoundaryPowers(const HigherOrderBasis &basis,
                                                        const std::array<int, 3> &nodes)
{
  TetrahedronFaceTracePowers powers;
  if (!IsGradientFamily(basis.family))
  {
    return powers;
  }
  const auto face_node = [&nodes](int element_node)
  {
    const auto node = std::find(nodes.begin(), nodes.end(), element_node);
    if (node == nodes.end())
    {
      throw std::invalid_argument(
          "A singular gradient trace does not contain its singular entity!");
    }
    return static_cast<int>(std::distance(nodes.begin(), node));
  };
  const double power = basis.nu - 1.0;
  if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
  {
    powers.node[face_node(basis.nodes[0])] = power;
  }
  else
  {
    const int first = face_node(basis.nodes[0]);
    const int second = face_node(basis.nodes[1]);
    powers.edge[first][second] = powers.edge[second][first] = power;
  }
  return powers;
}

TetrahedronFaceTracePowers AddTetrahedronBoundaryPowers(const TetrahedronFaceTracePowers &a,
                                                        const TetrahedronFaceTracePowers &b)
{
  TetrahedronFaceTracePowers result;
  for (int i = 0; i < 3; i++)
  {
    result.node[i] = a.node[i] + b.node[i];
    for (int j = 0; j < 3; j++)
    {
      result.edge[i][j] = a.edge[i][j] + b.edge[i][j];
    }
  }
  return result;
}

template <typename Integrand>
WeightedSegmentIntegral
IntegrateTetrahedronBoundaryTrace(const TetrahedronFaceTracePowers &powers, int value_size,
                                  const AdaptiveAssemblyOptions &options,
                                  const Integrand &integrand)
{
  if (value_size < 1)
  {
    throw std::invalid_argument(
        "Tetrahedral singular boundary trace received invalid integration dimensions!");
  }
  for (int i = 0; i < 3; i++)
  {
    if (!std::isfinite(powers.node[i]))
    {
      throw std::invalid_argument(
          "Tetrahedral singular boundary trace received a nonfinite node power!");
    }
    for (int j = 0; j < 3; j++)
    {
      const double scale =
          std::max({1.0, std::abs(powers.edge[i][j]), std::abs(powers.edge[j][i])});
      if (!std::isfinite(powers.edge[i][j]) ||
          std::abs(powers.edge[i][j] - powers.edge[j][i]) >
              64.0 * std::numeric_limits<double>::epsilon() * scale)
      {
        throw std::invalid_argument(
            "Tetrahedral singular boundary trace received invalid edge powers!");
      }
    }
  }

  const auto integrate = [&](int order)
  {
    std::vector<long double> accumulation(static_cast<std::size_t>(value_size));
    constexpr double radial_power = 3.0;
    constexpr double chart_determinant = 1.0 / 6.0;
    for (int vertex = 0; vertex < 3; vertex++)
    {
      for (int neighbor = 0; neighbor < 3; neighbor++)
      {
        if (neighbor == vertex)
        {
          continue;
        }
        int opposite = 0;
        while (opposite == vertex || opposite == neighbor)
        {
          opposite++;
        }
        const double outer_edge_power = powers.edge[vertex][neighbor];
        const double radial_singular_power = powers.node[vertex] +
                                             powers.edge[vertex][neighbor] +
                                             powers.edge[vertex][opposite];
        const double effective_radial_power =
            2.0 * radial_power - 1.0 + radial_power * radial_singular_power;
        if (!(outer_edge_power > -1.0) || !(effective_radial_power > -1.0))
        {
          throw std::domain_error(
              "A boundary face contains a nonintegrable singular edge or point. Split "
              "the boundary operator from the enriched feature or exclude enrichment "
              "near that boundary!");
        }
        const auto radial_rule =
            BuildWeightedSegmentQuadrature(order, effective_radial_power, 0.0);
        const double effective_tangential_power =
            radial_power * (outer_edge_power + 1.0) - 1.0;
        const auto tangential_rule =
            BuildWeightedSegmentQuadrature(order, effective_tangential_power, 0.0);
        for (const auto &radial_quadrature : radial_rule)
        {
          const double parameter = radial_quadrature.coordinate;
          const double radial = std::pow(parameter, radial_power);
          const double radial_weight =
              std::pow(parameter, radial_power * radial_singular_power);
          for (const auto &tangential_quadrature : tangential_rule)
          {
            const double tangential_parameter = tangential_quadrature.coordinate;
            const double tangential = std::pow(tangential_parameter, radial_power);
            const double tangential_weight =
                std::pow(tangential_parameter, radial_power * outer_edge_power);
            if (!std::isfinite(radial_weight) || !(radial_weight > 0.0) ||
                !std::isfinite(tangential_weight) || !(tangential_weight > 0.0))
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary trace produced an invalid chart "
                  "weight!");
            }

            std::array<double, 3> lambda{};
            for (int i = 0; i < 3; i++)
            {
              const double a = i == vertex ? 1.0 : 0.0;
              const double b = (i == vertex || i == neighbor) ? 0.5 : 0.0;
              constexpr double c = 1.0 / 3.0;
              lambda[i] =
                  (1.0 - radial) * a + radial * ((1.0 - tangential) * b + tangential * c);
            }
            mfem::IntegrationPoint point;
            point.Set2(lambda[1], lambda[2]);
            const auto values = integrand(point);
            if (values.size() != accumulation.size())
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary trace integrand changed dimensions!");
            }
            const double chart_scale =
                radial_quadrature.weight * tangential_quadrature.weight * radial_power *
                radial_power * chart_determinant / (radial_weight * tangential_weight);
            if (!std::isfinite(chart_scale) || !(chart_scale > 0.0))
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary trace produced an invalid chart "
                  "Jacobian!");
            }
            for (std::size_t i = 0; i < values.size(); i++)
            {
              const double contribution = chart_scale * values[i];
              if (!std::isfinite(contribution))
              {
                throw std::runtime_error(
                    "Tetrahedral singular boundary trace produced a nonfinite "
                    "contribution!");
              }
              accumulation[i] += contribution;
            }
          }
        }
      }
    }
    std::vector<double> result(accumulation.size());
    std::transform(accumulation.begin(), accumulation.end(), result.begin(),
                   [](long double value) { return static_cast<double>(value); });
    return result;
  };

  if (options.fixed_subdivision)
  {
    return {integrate(std::max(4, options.quadrature_order)),
            std::vector<double>(static_cast<std::size_t>(value_size))};
  }

  constexpr int order_increment = 2;
  int comparison_order = std::max(4, options.quadrature_order);
  auto comparison = integrate(comparison_order);
  std::vector<double> value;
  std::vector<double> estimated_absolute_error(static_cast<std::size_t>(value_size));
  int order = comparison_order + order_increment;
  for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
  {
    value = integrate(order);
    bool converged = true;
    for (std::size_t i = 0; i < value.size(); i++)
    {
      const double scale = std::max({1.0, std::abs(value[i]), std::abs(comparison[i])});
      const double error = 8.0 * std::abs(value[i] - comparison[i]) +
                           256.0 * std::numeric_limits<double>::epsilon() * scale;
      estimated_absolute_error[i] = error;
      const double tolerance =
          options.absolute_tolerance + options.relative_tolerance * std::abs(value[i]);
      converged = converged && std::isfinite(error) && error <= tolerance;
    }
    if (converged)
    {
      return {std::move(value), std::move(estimated_absolute_error)};
    }
    if (refinement < options.maximum_subdivisions)
    {
      comparison = value;
      comparison_order = order;
      order += order_increment;
    }
  }
  throw std::runtime_error(
      fmt::format("Tetrahedral singular boundary trace quadrature did not meet tolerance "
                  "(orders = {}/{}, order refinements = {})!",
                  order, comparison_order, options.maximum_subdivisions));
}

double TangentialDot(const Vector3 &a, const Vector3 &b, const Vector3 &normal)
{
  return Dot(a, b) - Dot(a, normal) * Dot(b, normal);
}

int GetTriangleVertex(const mfem::IntegrationPoint &point)
{
  const std::array<double, 3> lambda{1.0 - point.x - point.y, point.x, point.y};
  const auto maximum = std::max_element(lambda.begin(), lambda.end());
  if (*maximum < 1.0 - 256.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::runtime_error(
        "Singular boundary trace could not identify a triangle vertex!");
  }
  return static_cast<int>(std::distance(lambda.begin(), maximum));
}

std::array<int, 2> GetTriangleBoundaryNodes(mfem::Mesh &mesh, int boundary,
                                            mfem::FaceElementTransformations &face,
                                            mfem::IsoparametricTransformation &element1,
                                            mfem::IsoparametricTransformation &element2,
                                            int &element)
{
  std::array<int, 2> nodes{-1, -1};
  element = -1;
  for (int endpoint = 0; endpoint < 2; endpoint++)
  {
    mfem::IntegrationPoint point;
    point.x = static_cast<double>(endpoint);
    const int current =
        SetBoundaryIntegrationPoint(mesh, boundary, point, face, element1, element2);
    if (element >= 0 && current != element)
    {
      throw std::runtime_error(
          "Singular boundary trace changed adjacent element across one segment!");
    }
    element = current;
    nodes[endpoint] = GetTriangleVertex(face.Elem1->GetIntPoint());
  }
  if (nodes[0] < 0 || nodes[1] < 0 || nodes[0] == nodes[1])
  {
    throw std::runtime_error(
        "Singular boundary trace could not map a complete triangle edge!");
  }
  return nodes;
}

std::array<double, 2> GetTriangleBoundaryPowers(const TriangleBasis &basis,
                                                const std::array<int, 2> &nodes)
{
  if (basis.family != HigherOrderBasisFamily::NODE_GRADIENT)
  {
    throw std::invalid_argument(
        "Only triangular node-gradient functions have a nonzero boundary trace!");
  }
  std::array<double, 2> powers{};
  const int singular_node = TriangleSingularNode(basis);
  for (int endpoint = 0; endpoint < 2; endpoint++)
  {
    if (singular_node == nodes[endpoint])
    {
      powers[endpoint] = basis.nu - 1.0;
    }
  }
  return powers;
}

std::array<double, 2> GetTriangleBoundaryMonomialPowers(const TriangleBasis &basis,
                                                        const std::array<int, 2> &nodes,
                                                        double power)
{
  if (basis.family != HigherOrderBasisFamily::NODE_GRADIENT)
  {
    throw std::invalid_argument(
        "Only triangular node-gradient functions have a boundary monomial trace!");
  }
  std::array<double, 2> powers{};
  const int singular_node = TriangleSingularNode(basis);
  for (int endpoint = 0; endpoint < 2; endpoint++)
  {
    if (singular_node == nodes[endpoint])
    {
      powers[endpoint] = power;
    }
  }
  return powers;
}

bool HasTriangleBoundaryTrace(const TriangleBasis &basis, const std::array<int, 2> &nodes)
{
  if (basis.family == HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    return false;
  }
  if (basis.family != HigherOrderBasisFamily::NODE_GRADIENT)
  {
    throw std::invalid_argument(
        "Triangular boundary trace received an unsupported basis family!");
  }
  return std::find(nodes.begin(), nodes.end(), basis.nodes[0]) != nodes.end() &&
         std::find(nodes.begin(), nodes.end(), basis.nodes[1]) != nodes.end();
}

void AddScaledIntegral(WeightedSegmentIntegral &sum, double coefficient,
                       const WeightedSegmentIntegral &term)
{
  if (sum.value.empty())
  {
    sum.value.resize(term.value.size());
    sum.estimated_absolute_error.resize(term.estimated_absolute_error.size());
  }
  if (sum.value.size() != term.value.size() ||
      sum.estimated_absolute_error.size() != term.estimated_absolute_error.size())
  {
    throw std::invalid_argument(
        "Cannot combine triangular boundary integrals with different dimensions!");
  }
  for (std::size_t i = 0; i < sum.value.size(); i++)
  {
    sum.value[i] += coefficient * term.value[i];
    sum.estimated_absolute_error[i] +=
        std::abs(coefficient) * term.estimated_absolute_error[i];
  }
}

template <typename Integrand>
WeightedSegmentIntegral
IntegrateTriangleBoundaryTrace(const std::array<double, 2> &powers, int value_size,
                               const AdaptiveAssemblyOptions &options,
                               const Integrand &integrand)
{
  if (value_size < 1 || !std::isfinite(powers[0]) || !std::isfinite(powers[1]))
  {
    throw std::invalid_argument(
        "Triangular singular boundary trace received invalid integration dimensions!");
  }
  if (!(powers[0] > -1.0) || !(powers[1] > -1.0))
  {
    throw std::domain_error(
        "A boundary segment contains a nonintegrable singular endpoint. Split the "
        "boundary operator from the enriched feature or exclude enrichment near that "
        "boundary!");
  }

  const auto integrate = [&](int order)
  {
    std::vector<long double> accumulation(static_cast<std::size_t>(value_size));
    for (const auto &quadrature :
         BuildWeightedSegmentQuadrature(order, powers[0], powers[1]))
    {
      const auto values = integrand(quadrature.coordinate);
      if (values.size() != accumulation.size())
      {
        throw std::runtime_error(
            "Triangular singular boundary trace integrand changed dimensions!");
      }
      for (std::size_t i = 0; i < values.size(); i++)
      {
        const double contribution = quadrature.weight * values[i];
        if (!std::isfinite(contribution))
        {
          throw std::runtime_error(
              "Triangular singular boundary trace produced a nonfinite contribution!");
        }
        accumulation[i] += contribution;
      }
    }
    std::vector<double> result(accumulation.size());
    std::transform(accumulation.begin(), accumulation.end(), result.begin(),
                   [](long double value) { return static_cast<double>(value); });
    return result;
  };

  constexpr int order_increment = 4;
  int comparison_order = std::max(4, 2 * options.quadrature_order);
  auto comparison = integrate(comparison_order);
  std::vector<double> value;
  std::vector<double> estimated_absolute_error(static_cast<std::size_t>(value_size));
  double maximum_ratio = 0.0, maximum_error = 0.0, maximum_tolerance = 0.0;
  double maximum_value = 0.0, maximum_comparison = 0.0;
  std::size_t maximum_index = 0;
  int order = comparison_order + order_increment;
  for (int refinement = 0; refinement <= options.maximum_subdivisions; refinement++)
  {
    value = integrate(order);
    bool converged = true;
    for (std::size_t i = 0; i < value.size(); i++)
    {
      const double scale = std::max({1.0, std::abs(value[i]), std::abs(comparison[i])});
      const double error = 8.0 * std::abs(value[i] - comparison[i]) +
                           256.0 * std::numeric_limits<double>::epsilon() * scale;
      estimated_absolute_error[i] = error;
      const double tolerance =
          options.absolute_tolerance + options.relative_tolerance * std::abs(value[i]);
      const double ratio = error / tolerance;
      if (ratio > maximum_ratio)
      {
        maximum_ratio = ratio;
        maximum_error = error;
        maximum_tolerance = tolerance;
        maximum_value = value[i];
        maximum_comparison = comparison[i];
        maximum_index = i;
      }
      converged = converged && std::isfinite(error) && error <= tolerance;
    }
    if (converged)
    {
      return {std::move(value), std::move(estimated_absolute_error)};
    }
    if (refinement < options.maximum_subdivisions)
    {
      comparison = value;
      comparison_order = order;
      order += order_increment;
    }
  }
  throw std::runtime_error(
      fmt::format("Triangular singular boundary trace quadrature did not meet tolerance "
                  "(orders = {}/{}, order refinements = {}, worst value = {:.16e}, "
                  "comparison = {:.16e}, error = {:.3e}, tolerance = {:.3e}, "
                  "endpoint powers = [{:.12g}, {:.12g}], value size = {}, index = {})!",
                  order, comparison_order, options.maximum_subdivisions, maximum_value,
                  maximum_comparison, maximum_error, maximum_tolerance, powers[0],
                  powers[1], value_size, maximum_index));
}

ElementNDBoundaryMassMatrices AssembleTetrahedronElementNDBoundaryMassMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &nd_fe, mfem::Mesh &mesh,
    int boundary, double coefficient, const AdaptiveAssemblyOptions &options)
{
  if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3 ||
      nd_fe.GetGeomType() != mfem::Geometry::TETRAHEDRON ||
      nd_fe.GetRangeType() != mfem::FiniteElement::VECTOR ||
      nd_fe.GetMapType() != mfem::FiniteElement::H_CURL || !std::isfinite(coefficient))
  {
    throw std::invalid_argument(
        "Tetrahedral singular boundary mass assembly received invalid input!");
  }
  ValidateAdaptiveAssemblyOptions(options);
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
  int element = -1;
  const auto boundary_nodes =
      GetTetrahedronBoundaryNodes(mesh, boundary, face, element1, element2, element);

  const int standard_size = nd_fe.GetDof();
  const int enrichment_size = static_cast<int>(element_dofs.nd.size());
  auto *boundary_transformation = mesh.GetBdrElementTransformation(boundary);
  if (!boundary_transformation ||
      boundary_transformation->GetGeometryType() != mfem::Geometry::TRIANGLE)
  {
    throw std::runtime_error(
        "Tetrahedral singular boundary mass assembly requires triangle boundaries!");
  }

  mfem::DenseMatrix standard_value(standard_size, 3);
  const auto evaluate_traces = [&](const mfem::IntegrationPoint &boundary_point)
  {
    const int current = SetBoundaryIntegrationPoint(mesh, boundary, boundary_point, face,
                                                    element1, element2);
    if (current != element)
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary mass changed adjacent element!");
    }
    const auto &element_point = face.Elem1->GetIntPoint();
    double jacobian_determinant = 0.0;
    const auto grad_lambda =
        GetBarycentricGradients(*face.Elem1, element_point, jacobian_determinant);
    nd_fe.CalcPhysVShape(*face.Elem1, standard_value);

    boundary_transformation->SetIntPoint(&boundary_point);
    const auto &jacobian = boundary_transformation->Jacobian();
    if (jacobian.Height() != 3 || jacobian.Width() != 2)
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary mass found an invalid surface Jacobian!");
    }
    const Vector3 first{jacobian(0, 0), jacobian(1, 0), jacobian(2, 0)};
    const Vector3 second{jacobian(0, 1), jacobian(1, 1), jacobian(2, 1)};
    Vector3 normal = Cross(first, second);
    const double boundary_weight = std::sqrt(Dot(normal, normal));
    if (!std::isfinite(boundary_weight) || !(boundary_weight > 0.0))
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary mass found a degenerate triangle!");
    }
    for (double &value : normal)
    {
      value /= boundary_weight;
    }

    const BarycentricPoint lambda{1.0 - element_point.x - element_point.y - element_point.z,
                                  element_point.x, element_point.y, element_point.z};
    TetrahedronBoundaryTraceValues result;
    result.standard.resize(static_cast<std::size_t>(standard_size));
    for (int i = 0; i < standard_size; i++)
    {
      result.standard[i] = {standard_value(i, 0), standard_value(i, 1),
                            standard_value(i, 2)};
    }
    result.enrichment.reserve(static_cast<std::size_t>(enrichment_size));
    for (const auto &dof : element_dofs.nd)
    {
      result.enrichment.push_back(
          EvaluateHigherOrderBasis(lambda, grad_lambda, dof.basis).value);
    }
    result.normal = normal;
    result.scale = coefficient * boundary_weight;
    return result;
  };

  ElementNDBoundaryMassMatrices result;
  result.standard_enrichment.SetSize(standard_size, enrichment_size);
  result.enrichment_standard.SetSize(enrichment_size, standard_size);
  result.enrichment_enrichment.SetSize(enrichment_size);
  result.standard_enrichment_estimated_absolute_error.SetSize(standard_size,
                                                              enrichment_size);
  result.enrichment_enrichment_estimated_absolute_error.SetSize(enrichment_size);
  for (int enrichment = 0; enrichment < enrichment_size; enrichment++)
  {
    const auto &basis = element_dofs.nd[enrichment].basis;
    if (!HasTetrahedronBoundaryTrace(basis, boundary_nodes))
    {
      continue;
    }
    const auto powers = GetTetrahedronBoundaryPowers(basis, boundary_nodes);
    const auto integral = IntegrateTetrahedronBoundaryTrace(
        powers, standard_size, options,
        [&](const mfem::IntegrationPoint &point)
        {
          const auto values = evaluate_traces(point);
          std::vector<double> result(static_cast<std::size_t>(standard_size));
          for (int standard = 0; standard < standard_size; standard++)
          {
            result[standard] =
                values.scale * TangentialDot(values.standard[standard],
                                             values.enrichment[enrichment], values.normal);
          }
          return result;
        });
    for (int standard = 0; standard < standard_size; standard++)
    {
      result.standard_enrichment(standard, enrichment) = integral.value[standard];
      result.standard_enrichment_estimated_absolute_error(standard, enrichment) =
          integral.estimated_absolute_error[standard];
    }
  }
  result.enrichment_standard.Transpose(result.standard_enrichment);
  for (int row = 0; row < enrichment_size; row++)
  {
    const auto &row_basis = element_dofs.nd[row].basis;
    if (!HasTetrahedronBoundaryTrace(row_basis, boundary_nodes))
    {
      continue;
    }
    const auto row_powers = GetTetrahedronBoundaryPowers(row_basis, boundary_nodes);
    for (int column = row; column < enrichment_size; column++)
    {
      const auto &column_basis = element_dofs.nd[column].basis;
      if (!HasTetrahedronBoundaryTrace(column_basis, boundary_nodes))
      {
        continue;
      }
      const auto powers = AddTetrahedronBoundaryPowers(
          row_powers, GetTetrahedronBoundaryPowers(column_basis, boundary_nodes));
      const auto integral = IntegrateTetrahedronBoundaryTrace(
          powers, 1, options,
          [&](const mfem::IntegrationPoint &point)
          {
            const auto values = evaluate_traces(point);
            return std::vector<double>{
                values.scale * TangentialDot(values.enrichment[row],
                                             values.enrichment[column], values.normal)};
          });
      result.enrichment_enrichment(row, column) =
          result.enrichment_enrichment(column, row) = integral.value[0];
      result.enrichment_enrichment_estimated_absolute_error(row, column) =
          result.enrichment_enrichment_estimated_absolute_error(column, row) =
              integral.estimated_absolute_error[0];
    }
  }
  return result;
}

ElementNDBoundaryMassMatrices AssembleTriangleElementNDBoundaryMassMatrices(
    const TriangleElementDofMap &element_dofs, const mfem::FiniteElement &nd_fe,
    mfem::Mesh &mesh, int boundary, double coefficient,
    const AdaptiveAssemblyOptions &options)
{
  if (mesh.Dimension() != 2 || mesh.SpaceDimension() != 2 ||
      nd_fe.GetGeomType() != mfem::Geometry::TRIANGLE ||
      nd_fe.GetRangeType() != mfem::FiniteElement::VECTOR ||
      nd_fe.GetMapType() != mfem::FiniteElement::H_CURL || !std::isfinite(coefficient))
  {
    throw std::invalid_argument(
        "Triangular singular boundary mass assembly received invalid input!");
  }
  ValidateAdaptiveAssemblyOptions(options);
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
  int element = -1;
  const auto boundary_nodes =
      GetTriangleBoundaryNodes(mesh, boundary, face, element1, element2, element);

  const int standard_size = nd_fe.GetDof();
  const int enrichment_size = static_cast<int>(element_dofs.nd.size());
  mfem::IsoparametricTransformation boundary_transformation;
  mesh.GetBdrElementTransformation(boundary, &boundary_transformation);
  if (boundary_transformation.GetGeometryType() != mfem::Geometry::SEGMENT)
  {
    throw std::runtime_error(
        "Triangular singular boundary mass assembly requires segment boundaries!");
  }
  const auto reference_vertex = [](int node)
  { return Vector2{node == 1 ? 1.0 : 0.0, node == 2 ? 1.0 : 0.0}; };
  const auto reference_start = reference_vertex(boundary_nodes[0]);
  const auto reference_end = reference_vertex(boundary_nodes[1]);
  const Vector2 reference_tangent{reference_end[0] - reference_start[0],
                                  reference_end[1] - reference_start[1]};

  mfem::DenseMatrix standard_value(standard_size, 2);
  const auto evaluate_traces = [&](double coordinate)
  {
    mfem::IntegrationPoint boundary_point;
    boundary_point.x = coordinate;
    const int current = SetBoundaryIntegrationPoint(mesh, boundary, boundary_point, face,
                                                    element1, element2);
    if (current != element)
    {
      throw std::runtime_error(
          "Triangular singular boundary mass changed adjacent element!");
    }
    const auto &element_point = face.Elem1->GetIntPoint();
    nd_fe.CalcVShape(element_point, standard_value);

    boundary_transformation.SetIntPoint(&boundary_point);
    const double boundary_weight = boundary_transformation.Weight();
    if (!std::isfinite(boundary_weight) || !(boundary_weight > 0.0))
    {
      throw std::runtime_error(
          "Triangular singular boundary mass found a degenerate segment!");
    }

    std::vector<double> standard_trace(static_cast<std::size_t>(standard_size));
    for (int i = 0; i < standard_size; i++)
    {
      standard_trace[i] = standard_value(i, 0) * reference_tangent[0] +
                          standard_value(i, 1) * reference_tangent[1];
    }
    std::vector<double> enrichment_derivative(static_cast<std::size_t>(enrichment_size));
    for (int i = 0; i < enrichment_size; i++)
    {
      const auto &basis = element_dofs.nd[i].basis;
      if (!HasTriangleBoundaryTrace(basis, boundary_nodes))
      {
        continue;
      }
      const int singular_node = TriangleSingularNode(basis);
      enrichment_derivative[i] = singular_node == boundary_nodes[0] ? 1.0 : -1.0;
    }
    return std::make_tuple(std::move(standard_trace), std::move(enrichment_derivative),
                           coefficient / boundary_weight);
  };

  ElementNDBoundaryMassMatrices result;
  result.standard_enrichment.SetSize(standard_size, enrichment_size);
  result.enrichment_standard.SetSize(enrichment_size, standard_size);
  result.enrichment_enrichment.SetSize(enrichment_size);
  result.standard_enrichment_estimated_absolute_error.SetSize(standard_size,
                                                              enrichment_size);
  result.enrichment_enrichment_estimated_absolute_error.SetSize(enrichment_size);
  for (int enrichment = 0; enrichment < enrichment_size; enrichment++)
  {
    if (!HasTriangleBoundaryTrace(element_dofs.nd[enrichment].basis, boundary_nodes))
    {
      continue;
    }
    const auto &basis = element_dofs.nd[enrichment].basis;
    const auto integrate_term = [&](const std::array<double, 2> &powers)
    {
      return IntegrateTriangleBoundaryTrace(powers, standard_size, options,
                                            [&](double coordinate)
                                            {
                                              auto [standard_trace, enrichment_derivative,
                                                    scale] = evaluate_traces(coordinate);
                                              for (double &value : standard_trace)
                                              {
                                                value *= scale *
                                                         enrichment_derivative[enrichment];
                                              }
                                              return standard_trace;
                                            });
    };
    WeightedSegmentIntegral integral;
    AddScaledIntegral(integral, 1.0, integrate_term({0.0, 0.0}));
    AddScaledIntegral(integral, -basis.nu,
                      integrate_term(GetTriangleBoundaryPowers(basis, boundary_nodes)));
    for (int standard = 0; standard < standard_size; standard++)
    {
      result.standard_enrichment(standard, enrichment) = integral.value[standard];
      result.standard_enrichment_estimated_absolute_error(standard, enrichment) =
          integral.estimated_absolute_error[standard];
    }
  }
  result.enrichment_standard.Transpose(result.standard_enrichment);
  for (int row = 0; row < enrichment_size; row++)
  {
    if (!HasTriangleBoundaryTrace(element_dofs.nd[row].basis, boundary_nodes))
    {
      continue;
    }
    const auto row_powers =
        GetTriangleBoundaryPowers(element_dofs.nd[row].basis, boundary_nodes);
    for (int column = row; column < enrichment_size; column++)
    {
      if (!HasTriangleBoundaryTrace(element_dofs.nd[column].basis, boundary_nodes))
      {
        continue;
      }
      const auto column_powers =
          GetTriangleBoundaryPowers(element_dofs.nd[column].basis, boundary_nodes);
      const auto integrate_term = [&](const std::array<double, 2> &powers)
      {
        return IntegrateTriangleBoundaryTrace(
            powers, 1, options,
            [&](double coordinate)
            {
              auto [standard_trace, enrichment_derivative, scale] =
                  evaluate_traces(coordinate);
              return std::vector<double>{scale * enrichment_derivative[row] *
                                         enrichment_derivative[column]};
            });
      };
      const std::array<double, 2> combined_powers{row_powers[0] + column_powers[0],
                                                  row_powers[1] + column_powers[1]};
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      WeightedSegmentIntegral integral;
      AddScaledIntegral(integral, 1.0, integrate_term({0.0, 0.0}));
      AddScaledIntegral(integral, -row_basis.nu, integrate_term(row_powers));
      AddScaledIntegral(integral, -column_basis.nu, integrate_term(column_powers));
      AddScaledIntegral(integral, row_basis.nu * column_basis.nu,
                        integrate_term(combined_powers));
      result.enrichment_enrichment(row, column) =
          result.enrichment_enrichment(column, row) = integral.value[0];
      result.enrichment_enrichment_estimated_absolute_error(row, column) =
          result.enrichment_enrichment_estimated_absolute_error(column, row) =
              integral.estimated_absolute_error[0];
    }
  }
  return result;
}

ElementNDBoundaryMassMatrices AssembleTriangleElementH1BoundaryMassMatrices(
    const TriangleElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    mfem::Mesh &mesh, int boundary, double coefficient,
    const AdaptiveAssemblyOptions &options)
{
  if (mesh.Dimension() != 2 || mesh.SpaceDimension() != 2 ||
      h1_fe.GetGeomType() != mfem::Geometry::TRIANGLE ||
      h1_fe.GetRangeType() != mfem::FiniteElement::SCALAR ||
      h1_fe.GetMapType() != mfem::FiniteElement::VALUE || !std::isfinite(coefficient))
  {
    throw std::invalid_argument(
        "Triangular singular H1 boundary mass assembly received invalid input!");
  }
  ValidateAdaptiveAssemblyOptions(options);
  mfem::FaceElementTransformations face;
  mfem::IsoparametricTransformation element1, element2;
  int element = -1;
  const auto boundary_nodes =
      GetTriangleBoundaryNodes(mesh, boundary, face, element1, element2, element);

  const int standard_size = h1_fe.GetDof();
  const int enrichment_size = static_cast<int>(element_dofs.h1.size());
  mfem::IsoparametricTransformation boundary_transformation;
  mesh.GetBdrElementTransformation(boundary, &boundary_transformation);
  if (boundary_transformation.GetGeometryType() != mfem::Geometry::SEGMENT)
  {
    throw std::runtime_error(
        "Triangular singular H1 boundary mass assembly requires segment boundaries!");
  }

  mfem::Vector standard_value(standard_size);
  const auto evaluate_traces = [&](double coordinate)
  {
    mfem::IntegrationPoint boundary_point;
    boundary_point.x = coordinate;
    const int current = SetBoundaryIntegrationPoint(mesh, boundary, boundary_point, face,
                                                    element1, element2);
    if (current != element)
    {
      throw std::runtime_error(
          "Triangular singular H1 boundary mass changed adjacent element!");
    }
    const auto &element_point = face.Elem1->GetIntPoint();
    h1_fe.CalcShape(element_point, standard_value);

    boundary_transformation.SetIntPoint(&boundary_point);
    const double boundary_weight = boundary_transformation.Weight();
    if (!std::isfinite(boundary_weight) || !(boundary_weight > 0.0))
    {
      throw std::runtime_error(
          "Triangular singular H1 boundary mass found a degenerate segment!");
    }

    std::vector<double> standard_trace(static_cast<std::size_t>(standard_size));
    for (int i = 0; i < standard_size; i++)
    {
      standard_trace[i] = standard_value[i];
    }
    return std::make_pair(std::move(standard_trace), coefficient * boundary_weight);
  };

  ElementNDBoundaryMassMatrices result;
  result.standard_enrichment.SetSize(standard_size, enrichment_size);
  result.enrichment_standard.SetSize(enrichment_size, standard_size);
  result.enrichment_enrichment.SetSize(enrichment_size);
  result.standard_enrichment_estimated_absolute_error.SetSize(standard_size,
                                                              enrichment_size);
  result.enrichment_enrichment_estimated_absolute_error.SetSize(enrichment_size);
  for (int enrichment = 0; enrichment < enrichment_size; enrichment++)
  {
    const auto &basis = element_dofs.h1[enrichment].basis;
    if (!HasTriangleBoundaryTrace(basis, boundary_nodes))
    {
      continue;
    }
    const auto integrate_term = [&](const std::array<double, 2> &powers)
    {
      return IntegrateTriangleBoundaryTrace(powers, standard_size, options,
                                            [&](double coordinate)
                                            {
                                              auto [standard_trace, scale] =
                                                  evaluate_traces(coordinate);
                                              for (double &value : standard_trace)
                                              {
                                                value *= scale;
                                              }
                                              return standard_trace;
                                            });
    };
    WeightedSegmentIntegral integral;
    AddScaledIntegral(
        integral, 1.0,
        integrate_term(GetTriangleBoundaryMonomialPowers(basis, boundary_nodes, 1.0)));
    AddScaledIntegral(
        integral, -1.0,
        integrate_term(GetTriangleBoundaryMonomialPowers(basis, boundary_nodes, basis.nu)));
    for (int standard = 0; standard < standard_size; standard++)
    {
      result.standard_enrichment(standard, enrichment) = integral.value[standard];
      result.standard_enrichment_estimated_absolute_error(standard, enrichment) =
          integral.estimated_absolute_error[standard];
    }
  }
  result.enrichment_standard.Transpose(result.standard_enrichment);
  for (int row = 0; row < enrichment_size; row++)
  {
    const auto &row_basis = element_dofs.h1[row].basis;
    if (!HasTriangleBoundaryTrace(row_basis, boundary_nodes))
    {
      continue;
    }
    const auto row_linear_powers =
        GetTriangleBoundaryMonomialPowers(row_basis, boundary_nodes, 1.0);
    const auto row_singular_powers =
        GetTriangleBoundaryMonomialPowers(row_basis, boundary_nodes, row_basis.nu);
    for (int column = row; column < enrichment_size; column++)
    {
      const auto &column_basis = element_dofs.h1[column].basis;
      if (!HasTriangleBoundaryTrace(column_basis, boundary_nodes))
      {
        continue;
      }
      const auto column_linear_powers =
          GetTriangleBoundaryMonomialPowers(column_basis, boundary_nodes, 1.0);
      const auto column_singular_powers =
          GetTriangleBoundaryMonomialPowers(column_basis, boundary_nodes, column_basis.nu);
      const auto add_powers =
          [](const std::array<double, 2> &a, const std::array<double, 2> &b)
      { return std::array<double, 2>{a[0] + b[0], a[1] + b[1]}; };
      const auto integrate_term = [&](const std::array<double, 2> &powers)
      {
        return IntegrateTriangleBoundaryTrace(powers, 1, options,
                                              [&](double coordinate)
                                              {
                                                auto [standard_trace, scale] =
                                                    evaluate_traces(coordinate);
                                                return std::vector<double>{scale};
                                              });
      };
      WeightedSegmentIntegral integral;
      AddScaledIntegral(
          integral, 1.0,
          integrate_term(add_powers(row_linear_powers, column_linear_powers)));
      AddScaledIntegral(
          integral, -1.0,
          integrate_term(add_powers(row_linear_powers, column_singular_powers)));
      AddScaledIntegral(
          integral, -1.0,
          integrate_term(add_powers(row_singular_powers, column_linear_powers)));
      AddScaledIntegral(
          integral, 1.0,
          integrate_term(add_powers(row_singular_powers, column_singular_powers)));
      result.enrichment_enrichment(row, column) =
          result.enrichment_enrichment(column, row) = integral.value[0];
      result.enrichment_enrichment_estimated_absolute_error(row, column) =
          result.enrichment_enrichment_estimated_absolute_error(column, row) =
              integral.estimated_absolute_error[0];
    }
  }
  return result;
}

}  // namespace

LocalSparseOperatorBlocks AssembleLocalSparseNDBoundaryMassMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options,
    std::vector<LocalNDBoundaryPatchMatrices> *retained_patches,
    const std::set<int> &excluded_singular_attributes)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = nd_fespace.GetMesh();
  if (!mesh || mesh->Dimension() != 3 || mesh->SpaceDimension() != 3 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Tetrahedral singular boundary mass topology and space must share one "
        "three-dimensional mesh!");
  }
  for (const auto &[attribute, coefficient] : boundary_coefficients)
  {
    if (attribute <= 0 || !std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Tetrahedral singular boundary mass requires positive attributes and finite "
          "coefficients!");
    }
  }

  LocalSparseOperatorBlocks result;
  InitializeLocalSparseBlock(result, nd_fespace.GetVSize(),
                             static_cast<int>(topology.nd_dofs.size()));
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const auto coefficient = boundary_coefficients.find(mesh->GetBdrAttribute(boundary));
    if (coefficient == boundary_coefficients.end() || coefficient->second == 0.0)
    {
      continue;
    }

    int element = -1, local_face = -1;
    mesh->GetBdrElementAdjacentElement(boundary, element, local_face);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary mass could not find its adjacent element!");
    }
    const auto &element_dofs = topology.elements[element];
    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    nd_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    const bool exclude_singular =
        excluded_singular_attributes.find(mesh->GetBdrAttribute(boundary)) !=
        excluded_singular_attributes.end();
    mfem::Array<int> enrichment_dofs;
    if (!exclude_singular)
    {
      enrichment_dofs = GetElementEnrichmentDofs(element_dofs.nd, topology.nd_dofs.size());
    }

    // Retain a complete standard-plus-enrichment boundary facet in final unsigned local
    // coordinates. Standard trace assembly uses the boundary finite element and its own
    // face orientation transformation; singular coupling uses the adjacent volume element
    // transformation below.
    const auto retain_patch = [&](const ElementNDBoundaryMassMatrices *singular_matrices)
    {
      if (!retained_patches)
      {
        return;
      }
      LocalNDBoundaryPatchMatrices patch;
      patch.boundary = boundary;
      patch.element = element;
      for (int dof : standard_dofs)
      {
        patch.dofs.push_back(dof >= 0 ? dof : -1 - dof);
      }
      for (int dof : enrichment_dofs)
      {
        patch.dofs.push_back(nd_fespace.GetVSize() + dof);
      }
      patch.mass.SetSize(static_cast<int>(patch.dofs.size()));
      patch.mass = 0.0;

      mfem::Array<int> boundary_dofs;
      mfem::DofTransformation boundary_transformation;
      nd_fespace.GetBdrElementVDofs(boundary, boundary_dofs, boundary_transformation);
      mfem::ConstantCoefficient standard_coefficient(coefficient->second);
      mfem::VectorFEMassIntegrator standard_integrator(standard_coefficient);
      mfem::DenseMatrix standard_mass;
      standard_integrator.AssembleElementMatrix(
          *nd_fespace.GetBE(boundary), *mesh->GetBdrElementTransformation(boundary),
          standard_mass);
      boundary_transformation.TransformDual(standard_mass);
      for (int i = 0; i < boundary_dofs.Size(); i++)
      {
        const int global_i =
            boundary_dofs[i] >= 0 ? boundary_dofs[i] : -1 - boundary_dofs[i];
        const auto found_i = std::find(patch.dofs.begin(),
                                       patch.dofs.begin() + standard_dofs.Size(), global_i);
        MFEM_VERIFY(found_i != patch.dofs.begin() + standard_dofs.Size(),
                    "Boundary trace DOF is absent from its adjacent ND element!");
        const int row = static_cast<int>(found_i - patch.dofs.begin());
        const double row_sign = boundary_dofs[i] >= 0 ? 1.0 : -1.0;
        for (int j = 0; j < boundary_dofs.Size(); j++)
        {
          const int global_j =
              boundary_dofs[j] >= 0 ? boundary_dofs[j] : -1 - boundary_dofs[j];
          const auto found_j = std::find(
              patch.dofs.begin(), patch.dofs.begin() + standard_dofs.Size(), global_j);
          MFEM_VERIFY(found_j != patch.dofs.begin() + standard_dofs.Size(),
                      "Boundary trace DOF is absent from its adjacent ND element!");
          const int column = static_cast<int>(found_j - patch.dofs.begin());
          const double column_sign = boundary_dofs[j] >= 0 ? 1.0 : -1.0;
          patch.mass(row, column) += row_sign * column_sign * standard_mass(i, j);
        }
      }
      if (singular_matrices)
      {
        const int standard_size = standard_dofs.Size();
        for (int row = 0; row < standard_size; row++)
        {
          const double sign = standard_dofs[row] >= 0 ? 1.0 : -1.0;
          for (int column = 0; column < enrichment_dofs.Size(); column++)
          {
            patch.mass(row, standard_size + column) =
                sign * singular_matrices->standard_enrichment(row, column);
            patch.mass(standard_size + column, row) =
                sign * singular_matrices->enrichment_standard(column, row);
          }
        }
        for (int row = 0; row < enrichment_dofs.Size(); row++)
        {
          for (int column = 0; column < enrichment_dofs.Size(); column++)
          {
            patch.mass(standard_size + row, standard_size + column) =
                singular_matrices->enrichment_enrichment(row, column);
          }
        }
      }
      retained_patches->push_back(std::move(patch));
    };

    if (exclude_singular || element_dofs.nd.empty())
    {
      retain_patch(nullptr);
      continue;
    }

    ElementNDBoundaryMassMatrices matrices;
    try
    {
      matrices = AssembleTetrahedronElementNDBoundaryMassMatrices(
          element_dofs, *nd_fespace.GetFE(element), *mesh, boundary, coefficient->second,
          options);
    }
    catch (const std::exception &error)
    {
      throw std::runtime_error(fmt::format(
          "Tetrahedral singular boundary mass assembly failed on local boundary element "
          "{} (attribute {}, adjacent element {}): {}",
          boundary, mesh->GetBdrAttribute(boundary), element, error.what()));
    }

    ValidateStandardRowTransformation(
        dof_transformation, matrices.standard_enrichment, matrices.enrichment_standard,
        matrices.standard_enrichment_estimated_absolute_error);
    ApplyStandardRowTransformation(dof_transformation, matrices.standard_enrichment,
                                   matrices.enrichment_standard,
                                   matrices.standard_enrichment_estimated_absolute_error);
    const auto unsigned_standard_dofs = UnsignedDofs(standard_dofs);

    result.enrichment_enrichment->AddSubMatrix(enrichment_dofs, enrichment_dofs,
                                               matrices.enrichment_enrichment);
    result.standard_enrichment->AddSubMatrix(standard_dofs, enrichment_dofs,
                                             matrices.standard_enrichment);
    result.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        enrichment_dofs, enrichment_dofs,
        matrices.enrichment_enrichment_estimated_absolute_error);
    result.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        unsigned_standard_dofs, enrichment_dofs,
        matrices.standard_enrichment_estimated_absolute_error);
    retain_patch(&matrices);
  }
  FinalizeLocalSparseBlock(result);
  return result;
}

LocalSparseOperatorBlocks AssembleLocalSparseNDBoundaryMassMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = nd_fespace.GetMesh();
  if (!mesh || mesh->Dimension() != 2 || mesh->SpaceDimension() != 2 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Triangular singular boundary mass topology and space must share one "
        "two-dimensional mesh!");
  }
  for (const auto &[attribute, coefficient] : boundary_coefficients)
  {
    if (attribute <= 0 || !std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Triangular singular boundary mass requires positive attributes and finite "
          "coefficients!");
    }
  }

  LocalSparseOperatorBlocks result;
  InitializeLocalSparseBlock(result, nd_fespace.GetVSize(),
                             static_cast<int>(topology.nd_dofs.size()));
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const auto coefficient = boundary_coefficients.find(mesh->GetBdrAttribute(boundary));
    if (coefficient == boundary_coefficients.end() || coefficient->second == 0.0)
    {
      continue;
    }

    int element = -1, local_face = -1;
    mesh->GetBdrElementAdjacentElement(boundary, element, local_face);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Triangular singular boundary mass could not find its adjacent element!");
    }
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.nd.empty())
    {
      continue;
    }
    const auto enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.nd, topology.nd_dofs.size());

    ElementNDBoundaryMassMatrices matrices;
    try
    {
      matrices = AssembleTriangleElementNDBoundaryMassMatrices(
          element_dofs, *nd_fespace.GetFE(element), *mesh, boundary, coefficient->second,
          options);
    }
    catch (const std::exception &error)
    {
      throw std::runtime_error(fmt::format(
          "Triangular singular boundary mass assembly failed on local boundary element "
          "{} (attribute {}, adjacent element {}): {}",
          boundary, mesh->GetBdrAttribute(boundary), element, error.what()));
    }

    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    nd_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    ValidateStandardRowTransformation(
        dof_transformation, matrices.standard_enrichment, matrices.enrichment_standard,
        matrices.standard_enrichment_estimated_absolute_error);
    ApplyStandardRowTransformation(dof_transformation, matrices.standard_enrichment,
                                   matrices.enrichment_standard,
                                   matrices.standard_enrichment_estimated_absolute_error);
    const auto unsigned_standard_dofs = UnsignedDofs(standard_dofs);

    result.enrichment_enrichment->AddSubMatrix(enrichment_dofs, enrichment_dofs,
                                               matrices.enrichment_enrichment);
    result.standard_enrichment->AddSubMatrix(standard_dofs, enrichment_dofs,
                                             matrices.standard_enrichment);
    result.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        enrichment_dofs, enrichment_dofs,
        matrices.enrichment_enrichment_estimated_absolute_error);
    result.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        unsigned_standard_dofs, enrichment_dofs,
        matrices.standard_enrichment_estimated_absolute_error);
  }
  FinalizeLocalSparseBlock(result);
  return result;
}

LocalSparseOperatorBlocks AssembleLocalSparseH1BoundaryMassMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::map<int, double> &boundary_coefficients,
    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || mesh->Dimension() != 2 || mesh->SpaceDimension() != 2 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Triangular singular H1 boundary mass topology and space must share one "
        "two-dimensional mesh!");
  }
  for (const auto &[attribute, coefficient] : boundary_coefficients)
  {
    if (attribute <= 0 || !std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Triangular singular H1 boundary mass requires positive attributes and finite "
          "coefficients!");
    }
  }

  LocalSparseOperatorBlocks result;
  InitializeLocalSparseBlock(result, h1_fespace.GetVSize(),
                             static_cast<int>(topology.h1_dofs.size()));
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const auto coefficient = boundary_coefficients.find(mesh->GetBdrAttribute(boundary));
    if (coefficient == boundary_coefficients.end() || coefficient->second == 0.0)
    {
      continue;
    }

    int element = -1, local_face = -1;
    mesh->GetBdrElementAdjacentElement(boundary, element, local_face);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Triangular singular H1 boundary mass could not find its adjacent element!");
    }
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.h1.empty())
    {
      continue;
    }
    const auto enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.h1, topology.h1_dofs.size());

    ElementNDBoundaryMassMatrices matrices;
    try
    {
      matrices = AssembleTriangleElementH1BoundaryMassMatrices(
          element_dofs, *h1_fespace.GetFE(element), *mesh, boundary, coefficient->second,
          options);
    }
    catch (const std::exception &error)
    {
      throw std::runtime_error(fmt::format(
          "Triangular singular H1 boundary mass assembly failed on local boundary "
          "element {} (attribute {}, adjacent element {}): {}",
          boundary, mesh->GetBdrAttribute(boundary), element, error.what()));
    }

    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    h1_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    ValidateStandardRowTransformation(
        dof_transformation, matrices.standard_enrichment, matrices.enrichment_standard,
        matrices.standard_enrichment_estimated_absolute_error);
    ApplyStandardRowTransformation(dof_transformation, matrices.standard_enrichment,
                                   matrices.enrichment_standard,
                                   matrices.standard_enrichment_estimated_absolute_error);
    const auto unsigned_standard_dofs = UnsignedDofs(standard_dofs);

    result.enrichment_enrichment->AddSubMatrix(enrichment_dofs, enrichment_dofs,
                                               matrices.enrichment_enrichment);
    result.standard_enrichment->AddSubMatrix(standard_dofs, enrichment_dofs,
                                             matrices.standard_enrichment);
    result.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        enrichment_dofs, enrichment_dofs,
        matrices.enrichment_enrichment_estimated_absolute_error);
    result.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        unsigned_standard_dofs, enrichment_dofs,
        matrices.standard_enrichment_estimated_absolute_error);
  }
  FinalizeLocalSparseBlock(result);
  return result;
}

ParallelSparseOperatorBlocks
AssembleParallelSparseNDBoundaryMassMatrices(const LocalSparseOperatorBlocks &local,
                                             const ParallelDofNumbering &parallel_numbering,
                                             const mfem::ParFiniteElementSpace &nd_fespace)
{
  const MPI_Comm comm = nd_fespace.GetComm();
  auto enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, parallel_numbering.nd);
  const auto *standard_prolongation = nd_fespace.Dof_TrueDof_Matrix();
  if (!standard_prolongation)
  {
    throw std::runtime_error(
        "MFEM did not provide an ND finite-element prolongation matrix!");
  }
  auto absolute_prolongation = BuildAbsoluteProlongation(*standard_prolongation);

  ParallelSparseOperatorBlocks result;
  AssembleParallelSparseBlock(comm, local, *standard_prolongation, *absolute_prolongation,
                              *enrichment_prolongation, result);
  return result;
}

ParallelSparseOperatorBlocks
AssembleParallelSparseH1BoundaryMassMatrices(const LocalSparseOperatorBlocks &local,
                                             const ParallelDofNumbering &parallel_numbering,
                                             const mfem::ParFiniteElementSpace &h1_fespace)
{
  const MPI_Comm comm = h1_fespace.GetComm();
  auto enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, parallel_numbering.h1);
  const auto *standard_prolongation = h1_fespace.Dof_TrueDof_Matrix();
  if (!standard_prolongation)
  {
    throw std::runtime_error(
        "MFEM did not provide an H1 finite-element prolongation matrix!");
  }
  auto absolute_prolongation = BuildAbsoluteProlongation(*standard_prolongation);

  ParallelSparseOperatorBlocks result;
  AssembleParallelSparseBlock(comm, local, *standard_prolongation, *absolute_prolongation,
                              *enrichment_prolongation, result);
  return result;
}

ParallelSparseOperatorBlocks
RestrictParallelSparseOperatorBlocks(const ParallelSparseOperatorBlocks &fine,
                                     const mfem::HypreParMatrix &standard_prolongation)
{
  const bool empty = !fine.enrichment_enrichment && !fine.standard_enrichment &&
                     !fine.enrichment_standard &&
                     !fine.enrichment_enrichment_estimated_absolute_error &&
                     !fine.standard_enrichment_estimated_absolute_error;
  if (empty)
  {
    return {};
  }
  if (!fine.enrichment_enrichment || !fine.standard_enrichment ||
      !fine.enrichment_standard || !fine.enrichment_enrichment_estimated_absolute_error ||
      !fine.standard_enrichment_estimated_absolute_error ||
      standard_prolongation.Height() != fine.standard_enrichment->Height())
  {
    throw std::invalid_argument(
        "Cannot restrict an incomplete or incompatible singular sparse block!");
  }

  auto restriction =
      std::unique_ptr<mfem::HypreParMatrix>(standard_prolongation.Transpose());
  auto absolute_prolongation = BuildAbsoluteProlongation(standard_prolongation);
  auto absolute_restriction =
      std::unique_ptr<mfem::HypreParMatrix>(absolute_prolongation->Transpose());

  ParallelSparseOperatorBlocks coarse;
  coarse.enrichment_enrichment =
      std::make_unique<mfem::HypreParMatrix>(*fine.enrichment_enrichment);
  coarse.standard_enrichment.reset(
      mfem::ParMult(restriction.get(), fine.standard_enrichment.get()));
  coarse.enrichment_standard.reset(coarse.standard_enrichment->Transpose());
  coarse.enrichment_enrichment_estimated_absolute_error =
      std::make_unique<mfem::HypreParMatrix>(
          *fine.enrichment_enrichment_estimated_absolute_error);
  coarse.standard_enrichment_estimated_absolute_error.reset(mfem::ParMult(
      absolute_restriction.get(), fine.standard_enrichment_estimated_absolute_error.get()));
  if (!coarse.standard_enrichment || !coarse.standard_enrichment_estimated_absolute_error)
  {
    throw std::runtime_error("Failed to restrict a singular sparse block!");
  }
  return coarse;
}

ParallelSparseEnrichmentMatrices RestrictParallelSparseEnrichmentMatrices(
    const ParallelSparseEnrichmentMatrices &fine,
    const mfem::HypreParMatrix &h1_standard_prolongation,
    const mfem::HypreParMatrix &nd_standard_prolongation)
{
  ParallelSparseEnrichmentMatrices coarse;
  coarse.h1_diffusion =
      RestrictParallelSparseOperatorBlocks(fine.h1_diffusion, h1_standard_prolongation);
  coarse.h1_mass =
      RestrictParallelSparseOperatorBlocks(fine.h1_mass, h1_standard_prolongation);
  coarse.nd_mass =
      RestrictParallelSparseOperatorBlocks(fine.nd_mass, nd_standard_prolongation);
  coarse.nd_curl_curl =
      RestrictParallelSparseOperatorBlocks(fine.nd_curl_curl, nd_standard_prolongation);
  return coarse;
}

ParallelSparseOperatorBlocks
ProjectParallelSparseOperatorBlocksToH1(const ParallelSparseOperatorBlocks &nd,
                                        const mfem::HypreParMatrix &standard_gradient,
                                        const mfem::HypreParMatrix &enrichment_gradient)
{
  if (!nd.standard_enrichment || !nd.enrichment_standard || !nd.enrichment_enrichment)
  {
    throw std::invalid_argument(
        "Cannot project incomplete singular ND operator blocks into H1!");
  }
  if (standard_gradient.Height() != nd.standard_enrichment->Height() ||
      enrichment_gradient.Height() != nd.enrichment_enrichment->Height() ||
      nd.standard_enrichment->Width() != enrichment_gradient.Height() ||
      nd.enrichment_standard->Width() != standard_gradient.Height())
  {
    throw std::invalid_argument(
        "Cannot project incompatible singular ND operator blocks into H1!");
  }

  auto standard_restriction =
      std::unique_ptr<mfem::HypreParMatrix>(standard_gradient.Transpose());
  auto enrichment_restriction =
      std::unique_ptr<mfem::HypreParMatrix>(enrichment_gradient.Transpose());
  auto standard_enrichment_gradient = std::unique_ptr<mfem::HypreParMatrix>(
      mfem::ParMult(nd.standard_enrichment.get(), &enrichment_gradient, true));
  auto enrichment_enrichment_gradient = std::unique_ptr<mfem::HypreParMatrix>(
      mfem::ParMult(nd.enrichment_enrichment.get(), &enrichment_gradient, true));
  if (!standard_restriction || !enrichment_restriction || !standard_enrichment_gradient ||
      !enrichment_enrichment_gradient)
  {
    throw std::runtime_error("Failed to form singular ND-to-H1 gradient products!");
  }

  ParallelSparseOperatorBlocks h1;
  h1.standard_enrichment.reset(
      mfem::ParMult(standard_restriction.get(), standard_enrichment_gradient.get()));
  h1.enrichment_enrichment.reset(
      mfem::ParMult(enrichment_restriction.get(), enrichment_enrichment_gradient.get()));
  if (!h1.standard_enrichment || !h1.enrichment_enrichment)
  {
    throw std::runtime_error("Failed to project singular ND operator blocks into H1!");
  }
  h1.enrichment_standard.reset(h1.standard_enrichment->Transpose());
  return h1;
}

void AssembleParallelNDBoundaryLinearForm(const DofTopology &topology,
                                          const ParallelDofNumbering &parallel_numbering,
                                          mfem::ParFiniteElementSpace &nd_fespace,
                                          mfem::VectorCoefficient &coefficient,
                                          const AdaptiveAssemblyOptions &options,
                                          mfem::Vector &enrichment_true_dofs,
                                          const mfem::Array<int> *attribute_marker)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = nd_fespace.GetParMesh();
  if (!mesh || mesh->Dimension() != 3 || mesh->SpaceDimension() != 3 ||
      coefficient.GetVDim() != 3 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.nd_dofs.size() != parallel_numbering.nd.local_to_true.size() ||
      parallel_numbering.nd.owned_size >
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Tetrahedral singular boundary linear form received inconsistent topology, "
        "numbering, coefficient, or space!");
  }

  mfem::Vector local(static_cast<int>(topology.nd_dofs.size()));
  local = 0.0;
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const int attribute = mesh->GetBdrAttribute(boundary);
    if (attribute_marker && (attribute <= 0 || attribute > attribute_marker->Size() ||
                             !(*attribute_marker)[attribute - 1]))
    {
      continue;
    }
    int element = -1, local_face = -1;
    mesh->GetBdrElementAdjacentElement(boundary, element, local_face);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary linear form could not find its adjacent "
          "element!");
    }
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.nd.empty())
    {
      continue;
    }

    mfem::FaceElementTransformations face;
    mfem::IsoparametricTransformation element1, element2;
    int mapped_element = -1;
    const auto boundary_nodes = GetTetrahedronBoundaryNodes(*mesh, boundary, face, element1,
                                                            element2, mapped_element);
    if (mapped_element != element)
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary linear form found inconsistent adjacent "
          "elements!");
    }
    auto *boundary_transformation = mesh->GetBdrElementTransformation(boundary);
    if (!boundary_transformation ||
        boundary_transformation->GetGeometryType() != mfem::Geometry::TRIANGLE)
    {
      throw std::runtime_error(
          "Tetrahedral singular boundary linear form requires triangle boundaries!");
    }

    const int enrichment_size = static_cast<int>(element_dofs.nd.size());
    for (int i = 0; i < enrichment_size; i++)
    {
      const auto &descriptor = element_dofs.nd[i].basis;
      if (!HasTetrahedronBoundaryTrace(descriptor, boundary_nodes))
      {
        continue;
      }
      const auto powers = GetTetrahedronBoundaryPowers(descriptor, boundary_nodes);
      const auto integral = IntegrateTetrahedronBoundaryTrace(
          powers, 1, options,
          [&](const mfem::IntegrationPoint &boundary_point)
          {
            const int current = SetBoundaryIntegrationPoint(*mesh, boundary, boundary_point,
                                                            face, element1, element2);
            if (current != element)
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary linear form changed adjacent element!");
            }
            const auto &element_point = face.Elem1->GetIntPoint();
            double jacobian_determinant = 0.0;
            const auto grad_lambda =
                GetBarycentricGradients(*face.Elem1, element_point, jacobian_determinant);
            const BarycentricPoint lambda{
                1.0 - element_point.x - element_point.y - element_point.z, element_point.x,
                element_point.y, element_point.z};

            boundary_transformation->SetIntPoint(&boundary_point);
            const double boundary_weight = boundary_transformation->Weight();
            if (!std::isfinite(boundary_weight) || !(boundary_weight > 0.0))
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary linear form found a degenerate "
                  "triangle!");
            }
            mfem::Vector coefficient_value(3);
            coefficient.Eval(coefficient_value, *boundary_transformation, boundary_point);
            if (!std::isfinite(coefficient_value[0]) ||
                !std::isfinite(coefficient_value[1]) ||
                !std::isfinite(coefficient_value[2]))
            {
              throw std::runtime_error(
                  "Tetrahedral singular boundary linear form coefficient is nonfinite!");
            }

            const Vector3 value{coefficient_value[0], coefficient_value[1],
                                coefficient_value[2]};
            const auto basis = EvaluateHigherOrderBasis(lambda, grad_lambda, descriptor);
            return std::vector<double>{boundary_weight * Dot(basis.value, value)};
          });
      const std::size_t dof = element_dofs.nd[i].dof;
      if (dof >= topology.nd_dofs.size())
      {
        throw std::invalid_argument(
            "Tetrahedral singular boundary linear form has an invalid local DOF!");
      }
      local[static_cast<int>(dof)] += integral.value[0];
    }
  }

  auto prolongation =
      BuildParallelEnrichmentProlongation(nd_fespace.GetComm(), parallel_numbering.nd);
  enrichment_true_dofs.SetSize(static_cast<int>(parallel_numbering.nd.owned_size));
  enrichment_true_dofs = 0.0;
  prolongation->MultTranspose(local, enrichment_true_dofs);
}

void AssembleParallelNDBoundaryLinearForm(const TriangleDofTopology &topology,
                                          const ParallelDofNumbering &parallel_numbering,
                                          mfem::ParFiniteElementSpace &nd_fespace,
                                          mfem::VectorCoefficient &coefficient,
                                          const AdaptiveAssemblyOptions &options,
                                          mfem::Vector &enrichment_true_dofs,
                                          const mfem::Array<int> *attribute_marker)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = nd_fespace.GetParMesh();
  if (!mesh || mesh->Dimension() != 2 || mesh->SpaceDimension() != 2 ||
      coefficient.GetVDim() != 2 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.nd_dofs.size() != parallel_numbering.nd.local_to_true.size() ||
      parallel_numbering.nd.owned_size >
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Triangular singular boundary linear form received inconsistent topology, "
        "numbering, coefficient, or space!");
  }

  mfem::Vector local(static_cast<int>(topology.nd_dofs.size()));
  local = 0.0;
  for (int boundary = 0; boundary < mesh->GetNBE(); boundary++)
  {
    const int attribute = mesh->GetBdrAttribute(boundary);
    if (attribute_marker && (attribute <= 0 || attribute > attribute_marker->Size() ||
                             !(*attribute_marker)[attribute - 1]))
    {
      continue;
    }
    int element = -1, local_face = -1;
    mesh->GetBdrElementAdjacentElement(boundary, element, local_face);
    if (element < 0 || element >= mesh->GetNE())
    {
      throw std::runtime_error(
          "Triangular singular boundary linear form could not find its adjacent element!");
    }
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.nd.empty())
    {
      continue;
    }

    mfem::FaceElementTransformations face;
    mfem::IsoparametricTransformation element1, element2;
    int mapped_element = -1;
    const auto boundary_nodes =
        GetTriangleBoundaryNodes(*mesh, boundary, face, element1, element2, mapped_element);
    if (mapped_element != element)
    {
      throw std::runtime_error(
          "Triangular singular boundary linear form found inconsistent adjacent "
          "elements!");
    }
    auto *boundary_transformation = mesh->GetBdrElementTransformation(boundary);
    if (!boundary_transformation ||
        boundary_transformation->GetGeometryType() != mfem::Geometry::SEGMENT)
    {
      throw std::runtime_error(
          "Triangular singular boundary linear form requires segment boundaries!");
    }

    const int enrichment_size = static_cast<int>(element_dofs.nd.size());
    for (int i = 0; i < enrichment_size; i++)
    {
      const auto &basis = element_dofs.nd[i].basis;
      if (!HasTriangleBoundaryTrace(basis, boundary_nodes))
      {
        continue;
      }
      const int singular_node = TriangleSingularNode(basis);
      const double orientation = singular_node == boundary_nodes[0] ? 1.0 : -1.0;
      const auto integrate_term = [&](const std::array<double, 2> &powers)
      {
        return IntegrateTriangleBoundaryTrace(
            powers, 1, options,
            [&](double coordinate)
            {
              mfem::IntegrationPoint boundary_point;
              boundary_point.x = coordinate;
              const int current = SetBoundaryIntegrationPoint(
                  *mesh, boundary, boundary_point, face, element1, element2);
              if (current != element)
              {
                throw std::runtime_error(
                    "Triangular singular boundary linear form changed adjacent element!");
              }

              boundary_transformation->SetIntPoint(&boundary_point);
              const auto &jacobian = boundary_transformation->Jacobian();
              Vector2 tangent{jacobian(0, 0), jacobian(1, 0)};
              const double boundary_weight = std::sqrt(Dot(tangent, tangent));
              if (!std::isfinite(boundary_weight) || !(boundary_weight > 0.0))
              {
                throw std::runtime_error(
                    "Triangular singular boundary linear form found a degenerate "
                    "segment!");
              }
              tangent[0] /= boundary_weight;
              tangent[1] /= boundary_weight;

              mfem::Vector coefficient_value(2);
              coefficient.Eval(coefficient_value, *boundary_transformation, boundary_point);
              if (!std::isfinite(coefficient_value[0]) ||
                  !std::isfinite(coefficient_value[1]))
              {
                throw std::runtime_error(
                    "Triangular singular boundary linear form coefficient is nonfinite!");
              }
              const Vector2 value{coefficient_value[0], coefficient_value[1]};
              return std::vector<double>{orientation * Dot(tangent, value)};
            });
      };
      WeightedSegmentIntegral integral;
      AddScaledIntegral(integral, 1.0, integrate_term({0.0, 0.0}));
      AddScaledIntegral(integral, -basis.nu,
                        integrate_term(GetTriangleBoundaryPowers(basis, boundary_nodes)));
      const std::size_t dof = element_dofs.nd[i].dof;
      if (dof >= topology.nd_dofs.size())
      {
        throw std::invalid_argument(
            "Triangular singular boundary linear form has an invalid local DOF!");
      }
      local[static_cast<int>(dof)] += integral.value[0];
    }
  }

  auto prolongation =
      BuildParallelEnrichmentProlongation(nd_fespace.GetComm(), parallel_numbering.nd);
  enrichment_true_dofs.SetSize(static_cast<int>(parallel_numbering.nd.owned_size));
  enrichment_true_dofs = 0.0;
  prolongation->MultTranspose(local, enrichment_true_dofs);
}

LocalSparseEnrichmentMatrices AssembleLocalSparseEnrichmentMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options)
{
  auto batches = AssembleLocalSparseEnrichmentMatricesBatch(
      topology, h1_fespace, nd_fespace, {materials}, options);
  MFEM_ASSERT(batches.size() == 1, "Invalid single-material singular assembly result!");
  return std::move(batches.front());
}

LocalSparseEnrichmentMatrices AssembleLocalSparseEnrichmentMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || nd_fespace.GetMesh() != mesh || mesh->Dimension() != 2 ||
      mesh->SpaceDimension() != 2 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      materials.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_to_nd.size() != topology.h1_dofs.size() ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Triangular singular sparse topology, spaces, and materials must share one "
        "two-dimensional mesh!");
  }
  for (std::size_t h1 = 0; h1 < topology.h1_to_nd.size(); h1++)
  {
    const std::size_t nd = topology.h1_to_nd[h1];
    if (nd >= topology.nd_dofs.size() || !(topology.h1_dofs[h1] == topology.nd_dofs[nd]))
    {
      throw std::invalid_argument(
          "Triangular singular sparse assembly requires a consistent H1-to-ND "
          "topology!");
    }
  }
  for (const auto &material : materials)
  {
    ValidateMaterialCoefficients(material);
  }

  const int h1_enrichment_size = static_cast<int>(topology.h1_dofs.size());
  const int nd_enrichment_size = static_cast<int>(topology.nd_dofs.size());
  LocalSparseEnrichmentMatrices result;
  InitializeLocalSparseBlock(result.h1_diffusion, h1_fespace.GetVSize(),
                             h1_enrichment_size);
  InitializeLocalSparseBlock(result.h1_mass, h1_fespace.GetVSize(), h1_enrichment_size);
  InitializeLocalSparseBlock(result.nd_mass, nd_fespace.GetVSize(), nd_enrichment_size);
  InitializeLocalSparseBlock(result.nd_curl_curl, nd_fespace.GetVSize(),
                             nd_enrichment_size);

  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.h1.empty() && element_dofs.nd.empty())
    {
      continue;
    }
    const auto h1_enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.h1, topology.h1_dofs.size());
    const auto nd_enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.nd, topology.nd_dofs.size());

    auto &transformation = *mesh->GetElementTransformation(element);
    ElementEnrichmentMatrices enrichment;
    ElementStandardEnrichmentMatrices coupling;
    try
    {
      enrichment =
          AssembleTriangleElementEnrichmentMatrices(element_dofs, transformation, options);
      coupling = AssembleTriangleElementStandardEnrichmentMatrices(
          element_dofs, *h1_fespace.GetFE(element), *nd_fespace.GetFE(element),
          transformation, options);
    }
    catch (const std::exception &error)
    {
      const auto &basis = element_dofs.nd.front().basis;
      throw std::runtime_error(fmt::format(
          "Triangular singular sparse assembly failed on local mesh element {} "
          "(relative Jacobian variation = {:.3e}, first singular basis family = {}, "
          "nodes = [{}, {}, {}], nu = {:.17g}): {}",
          element, GetElementTransformationRelativeJacobianVariation(transformation),
          static_cast<int>(basis.family), basis.nodes[0], basis.nodes[1], basis.nodes[2],
          basis.nu, error.what()));
    }
    ApplyIsotropicMaterialCoefficients(materials[element], enrichment);
    ApplyIsotropicMaterialCoefficients(materials[element], coupling);

    mfem::Array<int> h1_standard_dofs, nd_standard_dofs;
    mfem::DofTransformation h1_dof_transformation, nd_dof_transformation;
    h1_fespace.GetElementVDofs(element, h1_standard_dofs, h1_dof_transformation);
    nd_fespace.GetElementVDofs(element, nd_standard_dofs, nd_dof_transformation);
    ApplyStandardDofTransformations(h1_dof_transformation, nd_dof_transformation, coupling);
    const auto h1_unsigned_dofs = UnsignedDofs(h1_standard_dofs);
    const auto nd_unsigned_dofs = UnsignedDofs(nd_standard_dofs);

    result.h1_diffusion.enrichment_enrichment->AddSubMatrix(
        h1_enrichment_dofs, h1_enrichment_dofs, enrichment.h1_diffusion);
    result.h1_diffusion.standard_enrichment->AddSubMatrix(
        h1_standard_dofs, h1_enrichment_dofs, coupling.h1_standard_enrichment);
    result.h1_diffusion.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        h1_enrichment_dofs, h1_enrichment_dofs,
        enrichment.h1_diffusion_estimated_absolute_error);
    result.h1_diffusion.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        h1_unsigned_dofs, h1_enrichment_dofs, coupling.h1_estimated_absolute_error);

    result.h1_mass.enrichment_enrichment->AddSubMatrix(
        h1_enrichment_dofs, h1_enrichment_dofs, enrichment.h1_mass);
    result.h1_mass.standard_enrichment->AddSubMatrix(h1_standard_dofs, h1_enrichment_dofs,
                                                     coupling.h1_mass_standard_enrichment);
    result.h1_mass.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        h1_enrichment_dofs, h1_enrichment_dofs,
        enrichment.h1_mass_estimated_absolute_error);
    result.h1_mass.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        h1_unsigned_dofs, h1_enrichment_dofs, coupling.h1_mass_estimated_absolute_error);

    result.nd_mass.enrichment_enrichment->AddSubMatrix(
        nd_enrichment_dofs, nd_enrichment_dofs, enrichment.nd_mass);
    result.nd_mass.standard_enrichment->AddSubMatrix(nd_standard_dofs, nd_enrichment_dofs,
                                                     coupling.nd_mass_standard_enrichment);
    result.nd_mass.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        nd_enrichment_dofs, nd_enrichment_dofs,
        enrichment.nd_mass_estimated_absolute_error);
    result.nd_mass.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        nd_unsigned_dofs, nd_enrichment_dofs, coupling.nd_mass_estimated_absolute_error);

    result.nd_curl_curl.enrichment_enrichment->AddSubMatrix(
        nd_enrichment_dofs, nd_enrichment_dofs, enrichment.nd_curl_curl);
    result.nd_curl_curl.standard_enrichment->AddSubMatrix(
        nd_standard_dofs, nd_enrichment_dofs, coupling.nd_curl_curl_standard_enrichment);
    result.nd_curl_curl.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        nd_enrichment_dofs, nd_enrichment_dofs,
        enrichment.nd_curl_curl_estimated_absolute_error);
    result.nd_curl_curl.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        nd_unsigned_dofs, nd_enrichment_dofs,
        coupling.nd_curl_curl_estimated_absolute_error);

    AccumulateQuadratureStatistics(enrichment, result.total_quadrature_leaf_count,
                                   result.maximum_subdivision_depth);
    AccumulateQuadratureStatistics(coupling, result.total_quadrature_leaf_count,
                                   result.maximum_subdivision_depth);
  }

  FinalizeLocalSparseBlock(result.h1_diffusion);
  FinalizeLocalSparseBlock(result.h1_mass);
  FinalizeLocalSparseBlock(result.nd_mass);
  FinalizeLocalSparseBlock(result.nd_curl_curl);
  return result;
}

namespace
{

mfem::DenseMatrix
BuildCombinedNDElementMatrix(const mfem::DenseMatrix &standard_standard,
                             double standard_coefficient,
                             const mfem::DenseMatrix &standard_enrichment,
                             const mfem::DenseMatrix &enrichment_standard,
                             const mfem::DenseMatrix &enrichment_enrichment)
{
  const int standard_size = standard_standard.Height();
  const int enrichment_size = enrichment_enrichment.Height();
  if (standard_standard.Width() != standard_size ||
      standard_enrichment.Height() != standard_size ||
      standard_enrichment.Width() != enrichment_size ||
      enrichment_standard.Height() != enrichment_size ||
      enrichment_standard.Width() != standard_size ||
      enrichment_enrichment.Width() != enrichment_size)
  {
    throw std::invalid_argument(
        "Coupled singular element patch received inconsistent matrix dimensions: ss=" +
        std::to_string(standard_size) +
        ", se=" + std::to_string(standard_enrichment.Height()) + "x" +
        std::to_string(standard_enrichment.Width()) +
        ", es=" + std::to_string(enrichment_standard.Height()) + "x" +
        std::to_string(enrichment_standard.Width()) +
        ", ee=" + std::to_string(enrichment_enrichment.Height()) + "x" +
        std::to_string(enrichment_enrichment.Width()) + "!");
  }

  mfem::DenseMatrix combined(standard_size + enrichment_size);
  combined = 0.0;
  for (int i = 0; i < standard_size; i++)
  {
    for (int j = 0; j < standard_size; j++)
    {
      combined(i, j) = standard_coefficient * standard_standard(i, j);
    }
    for (int j = 0; j < enrichment_size; j++)
    {
      combined(i, standard_size + j) = standard_enrichment(i, j);
      combined(standard_size + j, i) = enrichment_standard(j, i);
    }
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    for (int j = 0; j < enrichment_size; j++)
    {
      combined(standard_size + i, standard_size + j) = enrichment_enrichment(i, j);
    }
  }
  return combined;
}

mfem::DenseMatrix
BuildNDElementEnrichmentColumns(const mfem::DenseMatrix &standard_enrichment,
                                const mfem::DenseMatrix &enrichment_enrichment)
{
  const int standard_size = standard_enrichment.Height();
  const int enrichment_size = enrichment_enrichment.Height();
  if (standard_enrichment.Width() != enrichment_size ||
      enrichment_enrichment.Width() != enrichment_size)
  {
    throw std::invalid_argument(
        "Singular element enrichment strip received inconsistent dimensions!");
  }
  mfem::DenseMatrix columns(standard_size + enrichment_size, enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    for (int j = 0; j < enrichment_size; j++)
    {
      columns(i, j) = standard_enrichment(i, j);
    }
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    for (int j = 0; j < enrichment_size; j++)
    {
      columns(standard_size + i, j) = enrichment_enrichment(i, j);
    }
  }
  return columns;
}

mfem::DenseMatrix
BuildCombinedNDElementErrorMatrix(int standard_size,
                                  const mfem::DenseMatrix &standard_enrichment,
                                  const mfem::DenseMatrix &enrichment_enrichment)
{
  const int enrichment_size = enrichment_enrichment.Height();
  if (standard_size < 0 || standard_enrichment.Height() != standard_size ||
      standard_enrichment.Width() != enrichment_size ||
      enrichment_enrichment.Width() != enrichment_size)
  {
    throw std::invalid_argument(
        "Coupled singular element patch received inconsistent error matrix dimensions!");
  }

  mfem::DenseMatrix combined(standard_size + enrichment_size);
  combined = 0.0;
  for (int i = 0; i < standard_size; i++)
  {
    for (int j = 0; j < enrichment_size; j++)
    {
      const double error = standard_enrichment(i, j);
      combined(i, standard_size + j) = error;
      combined(standard_size + j, i) = error;
    }
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    for (int j = 0; j < enrichment_size; j++)
    {
      combined(standard_size + i, standard_size + j) = enrichment_enrichment(i, j);
    }
  }
  return combined;
}

template <typename Topology, typename ElementAssembler>
std::vector<LocalSparseEnrichmentMatrices> AssembleLocalSparseEnrichmentMatricesBatchImpl(
    const Topology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches,
    int retained_patch_batch, ElementAssembler &&assemble_element)
{
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); };
  auto *mesh = h1_fespace.GetMesh();
  if (material_batches.empty())
  {
    throw std::invalid_argument(
        "Singular sparse batch assembly requires at least one material field!");
  }
  if (retained_patch_batch < RetainNDElementPatchStripsOnly ||
      retained_patch_batch >= static_cast<int>(material_batches.size()))
  {
    throw std::invalid_argument(
        "Singular sparse batch assembly received an invalid retained patch batch!");
  }
  for (const auto &materials : material_batches)
  {
    if (materials.size() != static_cast<std::size_t>(mesh->GetNE()))
    {
      throw std::invalid_argument(
          "Singular sparse batch assembly material fields must cover every element!");
    }
    for (const auto &material : materials)
    {
      ValidateMaterialCoefficients(material);
    }
  }
  for (std::size_t h1 = 0; h1 < topology.h1_to_nd.size(); h1++)
  {
    const std::size_t nd = topology.h1_to_nd[h1];
    if (nd >= topology.nd_dofs.size() || !(topology.h1_dofs[h1] == topology.nd_dofs[nd]))
    {
      throw std::invalid_argument(
          "Singular sparse batch assembly requires a consistent H1-to-ND topology!");
    }
  }

  const int h1_enrichment_size = static_cast<int>(topology.h1_dofs.size());
  const int nd_enrichment_size = static_cast<int>(topology.nd_dofs.size());
  std::vector<LocalSparseEnrichmentMatrices> results(material_batches.size());
  for (auto &result : results)
  {
    InitializeLocalSparseBlock(result.h1_diffusion, h1_fespace.GetVSize(),
                               h1_enrichment_size);
    InitializeLocalSparseBlock(result.h1_mass, h1_fespace.GetVSize(), h1_enrichment_size);
    InitializeLocalSparseBlock(result.nd_mass, nd_fespace.GetVSize(), nd_enrichment_size);
    InitializeLocalSparseBlock(result.nd_curl_curl, nd_fespace.GetVSize(),
                               nd_enrichment_size);
  }

  double material_transformation_time = 0.0;
  double sparse_insertion_time = 0.0;
  mfem::VectorFEMassIntegrator standard_mass_integrator;
  mfem::CurlCurlIntegrator standard_curl_integrator;
  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.h1.empty() && element_dofs.nd.empty())
    {
      continue;
    }
    const auto h1_enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.h1, topology.h1_dofs.size());
    const auto nd_enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.nd, topology.nd_dofs.size());
    auto &transformation = *mesh->GetElementTransformation(element);
    const auto [unweighted_enrichment, unweighted_coupling] =
        assemble_element(element, element_dofs, transformation);

    mfem::Array<int> h1_standard_dofs, nd_standard_dofs;
    mfem::DofTransformation h1_dof_transformation, nd_dof_transformation;
    h1_fespace.GetElementVDofs(element, h1_standard_dofs, h1_dof_transformation);
    nd_fespace.GetElementVDofs(element, nd_standard_dofs, nd_dof_transformation);
    const auto h1_unsigned_dofs = UnsignedDofs(h1_standard_dofs);
    const auto nd_unsigned_dofs = UnsignedDofs(nd_standard_dofs);

    mfem::DenseMatrix standard_mass, standard_curl_curl;
    if (retained_patch_batch == RetainAllNDElementPatchBatches || retained_patch_batch >= 0)
    {
      standard_mass_integrator.AssembleElementMatrix(*nd_fespace.GetFE(element),
                                                     transformation, standard_mass);
      standard_curl_integrator.AssembleElementMatrix(*nd_fespace.GetFE(element),
                                                     transformation, standard_curl_curl);
      nd_dof_transformation.TransformDual(standard_mass);
      nd_dof_transformation.TransformDual(standard_curl_curl);
    }

    for (std::size_t batch = 0; batch < material_batches.size(); batch++)
    {
      const auto material_start = Clock::now();
      auto enrichment = unweighted_enrichment;
      auto coupling = unweighted_coupling;
      ApplyIsotropicMaterialCoefficients(material_batches[batch][element], enrichment);
      ApplyIsotropicMaterialCoefficients(material_batches[batch][element], coupling);
      ApplyStandardDofTransformations(h1_dof_transformation, nd_dof_transformation,
                                      coupling);
      material_transformation_time += elapsed(material_start);
      auto &result = results[batch];

      if (retained_patch_batch == RetainNDElementPatchStripsOnly ||
          retained_patch_batch == RetainAllNDElementPatchBatches ||
          static_cast<int>(batch) == retained_patch_batch)
      {
        const auto &material = material_batches[batch][element];
        LocalNDElementPatchMatrices patch;
        patch.element = element;
        patch.standard_dofs = nd_standard_dofs;
        patch.enrichment_dofs = nd_enrichment_dofs;
        if (retained_patch_batch == RetainNDElementPatchStripsOnly)
        {
          patch.mass = BuildNDElementEnrichmentColumns(coupling.nd_mass_standard_enrichment,
                                                       enrichment.nd_mass);
          patch.curl_curl = BuildNDElementEnrichmentColumns(
              coupling.nd_curl_curl_standard_enrichment, enrichment.nd_curl_curl);
          result.nd_element_patches.push_back(std::move(patch));
          continue;
        }
        const int combined_nd_size = standard_mass.Height() + nd_enrichment_dofs.Size();
        if (material.electric == 0.0)
        {
          patch.mass.SetSize(combined_nd_size);
          patch.mass = 0.0;
          patch.mass_estimated_absolute_error.SetSize(combined_nd_size);
          patch.mass_estimated_absolute_error = 0.0;
        }
        else
        {
          patch.mass = BuildCombinedNDElementMatrix(
              standard_mass, material.electric, coupling.nd_mass_standard_enrichment,
              coupling.nd_mass_enrichment_standard, enrichment.nd_mass);
          patch.mass_estimated_absolute_error = BuildCombinedNDElementErrorMatrix(
              standard_mass.Height(), coupling.nd_mass_estimated_absolute_error,
              enrichment.nd_mass_estimated_absolute_error);
        }
        if (material.inverse_magnetic == 0.0)
        {
          patch.curl_curl.SetSize(combined_nd_size);
          patch.curl_curl = 0.0;
          patch.curl_curl_estimated_absolute_error.SetSize(combined_nd_size);
          patch.curl_curl_estimated_absolute_error = 0.0;
        }
        else
        {
          patch.curl_curl = BuildCombinedNDElementMatrix(
              standard_curl_curl, material.inverse_magnetic,
              coupling.nd_curl_curl_standard_enrichment,
              coupling.nd_curl_curl_enrichment_standard, enrichment.nd_curl_curl);
          patch.curl_curl_estimated_absolute_error = BuildCombinedNDElementErrorMatrix(
              standard_curl_curl.Height(), coupling.nd_curl_curl_estimated_absolute_error,
              enrichment.nd_curl_curl_estimated_absolute_error);
        }
        result.nd_element_patches.push_back(std::move(patch));
      }

      const auto insertion_start = Clock::now();
      result.h1_diffusion.enrichment_enrichment->AddSubMatrix(
          h1_enrichment_dofs, h1_enrichment_dofs, enrichment.h1_diffusion);
      result.h1_diffusion.standard_enrichment->AddSubMatrix(
          h1_standard_dofs, h1_enrichment_dofs, coupling.h1_standard_enrichment);
      result.h1_diffusion.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
          h1_enrichment_dofs, h1_enrichment_dofs,
          enrichment.h1_diffusion_estimated_absolute_error);
      result.h1_diffusion.standard_enrichment_estimated_absolute_error->AddSubMatrix(
          h1_unsigned_dofs, h1_enrichment_dofs, coupling.h1_estimated_absolute_error);

      result.h1_mass.enrichment_enrichment->AddSubMatrix(
          h1_enrichment_dofs, h1_enrichment_dofs, enrichment.h1_mass);
      result.h1_mass.standard_enrichment->AddSubMatrix(
          h1_standard_dofs, h1_enrichment_dofs, coupling.h1_mass_standard_enrichment);
      result.h1_mass.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
          h1_enrichment_dofs, h1_enrichment_dofs,
          enrichment.h1_mass_estimated_absolute_error);
      result.h1_mass.standard_enrichment_estimated_absolute_error->AddSubMatrix(
          h1_unsigned_dofs, h1_enrichment_dofs, coupling.h1_mass_estimated_absolute_error);

      result.nd_mass.enrichment_enrichment->AddSubMatrix(
          nd_enrichment_dofs, nd_enrichment_dofs, enrichment.nd_mass);
      result.nd_mass.standard_enrichment->AddSubMatrix(
          nd_standard_dofs, nd_enrichment_dofs, coupling.nd_mass_standard_enrichment);
      result.nd_mass.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
          nd_enrichment_dofs, nd_enrichment_dofs,
          enrichment.nd_mass_estimated_absolute_error);
      result.nd_mass.standard_enrichment_estimated_absolute_error->AddSubMatrix(
          nd_unsigned_dofs, nd_enrichment_dofs, coupling.nd_mass_estimated_absolute_error);

      result.nd_curl_curl.enrichment_enrichment->AddSubMatrix(
          nd_enrichment_dofs, nd_enrichment_dofs, enrichment.nd_curl_curl);
      result.nd_curl_curl.standard_enrichment->AddSubMatrix(
          nd_standard_dofs, nd_enrichment_dofs, coupling.nd_curl_curl_standard_enrichment);
      result.nd_curl_curl.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
          nd_enrichment_dofs, nd_enrichment_dofs,
          enrichment.nd_curl_curl_estimated_absolute_error);
      result.nd_curl_curl.standard_enrichment_estimated_absolute_error->AddSubMatrix(
          nd_unsigned_dofs, nd_enrichment_dofs,
          coupling.nd_curl_curl_estimated_absolute_error);
      sparse_insertion_time += elapsed(insertion_start);

      AccumulateQuadratureStatistics(enrichment, result.total_quadrature_leaf_count,
                                     result.maximum_subdivision_depth);
      AccumulateQuadratureStatistics(coupling, result.total_quadrature_leaf_count,
                                     result.maximum_subdivision_depth);
    }
  }

  const auto finalization_start = Clock::now();
  for (auto &result : results)
  {
    FinalizeLocalSparseBlock(result.h1_diffusion);
    FinalizeLocalSparseBlock(result.h1_mass);
    FinalizeLocalSparseBlock(result.nd_mass);
    FinalizeLocalSparseBlock(result.nd_curl_curl);
  }
  const double sparse_finalization_time = elapsed(finalization_start);
  for (auto &result : results)
  {
    result.material_transformation_time = material_transformation_time;
    result.sparse_insertion_time = sparse_insertion_time;
    result.sparse_finalization_time = sparse_finalization_time;
  }
  return results;
}

}  // namespace

std::vector<LocalSparseEnrichmentMatrices> AssembleLocalSparseEnrichmentMatricesBatch(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches,
    const AdaptiveAssemblyOptions &options, int retained_patch_batch)
{
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); };
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || nd_fespace.GetMesh() != mesh ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_to_nd.size() != topology.h1_dofs.size() ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Singular sparse batch topology and spaces must share one mesh!");
  }
  AdaptiveReferenceTable reference_table(options);
  AffineStandardReferenceTable standard_reference_table(options);
  AffineCommonNDMassReferenceTable common_mass_reference_table(options);
  AffineCommonNDCurlReferenceTable common_curl_reference_table(options);
  AffineScalarReferenceTable scalar_reference_table(options);
  double enrichment_evaluation_time = 0.0;
  double standard_enrichment_evaluation_time = 0.0;
  double standard_reference_generation_time = 0.0;
  double standard_enrichment_setup_time = 0.0;
  double nd_coupling_time = 0.0;
  double h1_gradient_coupling_time = 0.0;
  double h1_mass_coupling_time = 0.0;
  std::size_t affine_nd_mass_contraction_count = 0;
  std::size_t affine_nd_mass_reintegration_count = 0;
  std::size_t affine_nd_mass_reintegration_batch_count = 0;
  std::size_t affine_nd_curl_contraction_count = 0;
  std::size_t affine_nd_curl_reintegration_count = 0;
  std::size_t affine_nd_curl_reintegration_batch_count = 0;
  auto results = AssembleLocalSparseEnrichmentMatricesBatchImpl(
      topology, h1_fespace, nd_fespace, material_batches, retained_patch_batch,
      [&](int element, const ElementDofMap &element_dofs,
          mfem::ElementTransformation &transformation)
      {
        try
        {
          const double generation_start = standard_reference_table.GenerationTime() +
                                          common_mass_reference_table.GenerationTime() +
                                          common_curl_reference_table.GenerationTime() +
                                          scalar_reference_table.GenerationTime();
          std::optional<AffineCommonNDMassBlocks> common_mass;
          std::optional<AffineCommonNDCurlBlocks> common_curl;
          if (IsAffineElementTransformation(transformation))
          {
            double jacobian_determinant;
            const auto grad_lambda =
                GetAffineBarycentricGradients(transformation, jacobian_determinant);
            common_mass =
                common_mass_reference_table.Get(*nd_fespace.GetFE(element), element_dofs.nd,
                                                grad_lambda, jacobian_determinant);
            common_curl =
                common_curl_reference_table.Get(*nd_fespace.GetFE(element), element_dofs.nd,
                                                grad_lambda, jacobian_determinant);
          }
          const auto enrichment_start = Clock::now();
          auto enrichment = AssembleElementEnrichmentMatricesCached(
              element_dofs, transformation, options, reference_table,
              scalar_reference_table, common_mass ? &*common_mass : nullptr,
              common_curl ? &*common_curl : nullptr);
          enrichment_evaluation_time += elapsed(enrichment_start);
          const auto coupling_start = Clock::now();
          auto coupling = AssembleElementStandardEnrichmentMatricesCached(
              element_dofs, *h1_fespace.GetFE(element), *nd_fespace.GetFE(element),
              transformation, options, standard_reference_table, scalar_reference_table,
              common_mass ? &*common_mass : nullptr, common_curl ? &*common_curl : nullptr);
          if (common_mass)
          {
            ApplyCommonAffineNDMassBlocks(element_dofs, *h1_fespace.GetFE(element),
                                          *nd_fespace.GetFE(element), transformation,
                                          *common_mass, enrichment, coupling);
          }
          standard_enrichment_evaluation_time += elapsed(coupling_start);
          standard_reference_generation_time +=
              standard_reference_table.GenerationTime() +
              common_mass_reference_table.GenerationTime() +
              common_curl_reference_table.GenerationTime() +
              scalar_reference_table.GenerationTime() - generation_start;
          standard_enrichment_setup_time += coupling.setup_time;
          nd_coupling_time += coupling.nd_coupling_time;
          h1_gradient_coupling_time += coupling.h1_gradient_coupling_time;
          h1_mass_coupling_time += coupling.h1_mass_coupling_time;
          affine_nd_mass_contraction_count += coupling.affine_nd_mass_contraction_count;
          affine_nd_mass_reintegration_count += coupling.affine_nd_mass_reintegration_count;
          affine_nd_mass_reintegration_batch_count +=
              coupling.affine_nd_mass_reintegration_batch_count;
          affine_nd_curl_contraction_count += coupling.affine_nd_curl_contraction_count;
          affine_nd_curl_reintegration_count += coupling.affine_nd_curl_reintegration_count;
          affine_nd_curl_reintegration_batch_count +=
              coupling.affine_nd_curl_reintegration_batch_count;
          return std::make_pair(std::move(enrichment), std::move(coupling));
        }
        catch (const std::exception &error)
        {
          throw std::runtime_error(
              fmt::format("Singular sparse batch assembly failed on local mesh element "
                          "{}: {}",
                          element, error.what()));
        }
      });
  for (auto &result : results)
  {
    result.affine_reference_table_entries =
        reference_table.Size() + standard_reference_table.Size() +
        common_mass_reference_table.Size() + common_curl_reference_table.Size() +
        scalar_reference_table.Size();
    result.affine_reference_cache_hits =
        reference_table.Hits() + standard_reference_table.Hits() +
        common_mass_reference_table.Hits() + common_curl_reference_table.Hits() +
        scalar_reference_table.Hits();
    result.enrichment_evaluation_time = enrichment_evaluation_time;
    result.standard_enrichment_evaluation_time = standard_enrichment_evaluation_time;
    result.standard_reference_generation_time = standard_reference_generation_time;
    result.standard_enrichment_setup_time = standard_enrichment_setup_time;
    result.nd_coupling_time = nd_coupling_time;
    result.h1_gradient_coupling_time = h1_gradient_coupling_time;
    result.h1_mass_coupling_time = h1_mass_coupling_time;
    result.affine_nd_mass_contraction_count = affine_nd_mass_contraction_count;
    result.affine_nd_mass_reintegration_count = affine_nd_mass_reintegration_count;
    result.affine_nd_mass_reintegration_batch_count =
        affine_nd_mass_reintegration_batch_count;
    result.affine_nd_curl_contraction_count = affine_nd_curl_contraction_count;
    result.affine_nd_curl_reintegration_count = affine_nd_curl_reintegration_count;
    result.affine_nd_curl_reintegration_batch_count =
        affine_nd_curl_reintegration_batch_count;
  }
  return results;
}

std::vector<LocalSparseEnrichmentMatrices> AssembleLocalSparseEnrichmentMatricesBatch(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches,
    const AdaptiveAssemblyOptions &options, int retained_patch_batch)
{
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](const Clock::time_point &start)
  { return std::chrono::duration<double>(Clock::now() - start).count(); };
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || nd_fespace.GetMesh() != mesh || mesh->Dimension() != 2 ||
      mesh->SpaceDimension() != 2 ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_to_nd.size() != topology.h1_dofs.size() ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Triangular singular sparse batch topology and spaces must share one "
        "two-dimensional mesh!");
  }
  double enrichment_evaluation_time = 0.0;
  double standard_enrichment_evaluation_time = 0.0;
  auto results = AssembleLocalSparseEnrichmentMatricesBatchImpl(
      topology, h1_fespace, nd_fespace, material_batches, retained_patch_batch,
      [&](int element, const TriangleElementDofMap &element_dofs,
          mfem::ElementTransformation &transformation)
      {
        try
        {
          const auto enrichment_start = Clock::now();
          auto enrichment = AssembleTriangleElementEnrichmentMatrices(
              element_dofs, transformation, options);
          enrichment_evaluation_time += elapsed(enrichment_start);
          const auto coupling_start = Clock::now();
          auto coupling = AssembleTriangleElementStandardEnrichmentMatrices(
              element_dofs, *h1_fespace.GetFE(element), *nd_fespace.GetFE(element),
              transformation, options);
          standard_enrichment_evaluation_time += elapsed(coupling_start);
          return std::make_pair(std::move(enrichment), std::move(coupling));
        }
        catch (const std::exception &error)
        {
          const auto &basis = element_dofs.nd.front().basis;
          throw std::runtime_error(fmt::format(
              "Triangular singular sparse batch assembly failed on local mesh element {} "
              "(relative Jacobian variation = {:.3e}, first singular basis family = {}, "
              "nodes = [{}, {}, {}], nu = {:.17g}): {}",
              element, GetElementTransformationRelativeJacobianVariation(transformation),
              static_cast<int>(basis.family), basis.nodes[0], basis.nodes[1],
              basis.nodes[2], basis.nu, error.what()));
        }
      });
  for (auto &result : results)
  {
    result.enrichment_evaluation_time = enrichment_evaluation_time;
    result.standard_enrichment_evaluation_time = standard_enrichment_evaluation_time;
  }
  return results;
}

std::vector<LocalInterpolatedNDDiagonals> AssembleLocalInterpolatedNDDomainDiagonals(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches)
{
  return AssembleLocalInterpolatedNDDomainDiagonalsImpl(
      topology, nd_fespace, local_interpolant, material_batches);
}

mfem::Vector AssembleLocalInterpolatedNDBoundaryDiagonal(
    const DofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::map<int, double> &boundary_coefficients)
{
  return AssembleLocalInterpolatedNDBoundaryDiagonalImpl(
      topology, nd_fespace, local_interpolant, boundary_coefficients);
}

mfem::Vector AssembleLocalInterpolatedNDBoundaryDiagonal(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::map<int, double> &boundary_coefficients)
{
  return AssembleLocalInterpolatedNDBoundaryDiagonalImpl(
      topology, nd_fespace, local_interpolant, boundary_coefficients);
}

std::vector<LocalInterpolatedNDDiagonals> AssembleLocalInterpolatedNDDomainDiagonals(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &nd_fespace,
    const mfem::SparseMatrix &local_interpolant,
    const std::vector<std::vector<IsotropicMaterialCoefficients>> &material_batches)
{
  return AssembleLocalInterpolatedNDDomainDiagonalsImpl(
      topology, nd_fespace, local_interpolant, material_batches);
}

void SetLocalTransformedEnrichmentDiagonal(
    LocalSparseOperatorBlocks &blocks, const mfem::SparseMatrix &local_interpolant,
    const mfem::Vector &interpolated_standard_diagonal)
{
  if (!blocks.enrichment_enrichment || !blocks.standard_enrichment ||
      !blocks.enrichment_standard ||
      local_interpolant.Height() != blocks.standard_enrichment->Height() ||
      local_interpolant.Width() != blocks.standard_enrichment->Width() ||
      interpolated_standard_diagonal.Size() != local_interpolant.Width())
  {
    throw std::invalid_argument(
        "Cannot assemble an inconsistent transformed enrichment diagonal!");
  }
  std::unique_ptr<mfem::SparseMatrix> interpolant_transpose(
      mfem::Transpose(local_interpolant));
  std::unique_ptr<mfem::SparseMatrix> coupling_transpose(
      mfem::Transpose(*blocks.standard_enrichment));
  if (!interpolant_transpose || !coupling_transpose)
  {
    throw std::runtime_error("Failed to transpose local transformed-diagonal factors!");
  }

  auto diagonal = std::make_unique<mfem::Vector>(interpolated_standard_diagonal);
  mfem::Vector enrichment_diagonal(local_interpolant.Width());
  blocks.enrichment_enrichment->GetDiag(enrichment_diagonal);
  *diagonal += enrichment_diagonal;
  for (int enrichment = 0; enrichment < diagonal->Size(); enrichment++)
  {
    (*diagonal)[enrichment] -=
        2.0 * SparseRowInnerProduct(*interpolant_transpose, enrichment, *coupling_transpose,
                                    enrichment);
  }
  blocks.transformed_enrichment_diagonal = std::move(diagonal);
}

LocalSparseH1EnrichmentMatrices AssembleLocalSparseH1EnrichmentMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      materials.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_to_nd.size() != topology.h1_dofs.size() ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Singular H1 sparse assembly topology, space, and materials must share one mesh!");
  }
  for (std::size_t h1 = 0; h1 < topology.h1_to_nd.size(); h1++)
  {
    const std::size_t nd = topology.h1_to_nd[h1];
    if (nd >= topology.nd_dofs.size() || !(topology.h1_dofs[h1] == topology.nd_dofs[nd]))
    {
      throw std::invalid_argument(
          "Singular H1 sparse assembly requires a consistent H1-to-ND topology!");
    }
  }
  for (const auto &material : materials)
  {
    ValidateMaterialCoefficients(material);
  }

  const int enrichment_size = static_cast<int>(topology.h1_dofs.size());
  LocalSparseH1EnrichmentMatrices result;
  InitializeLocalSparseBlock(result.diffusion, h1_fespace.GetVSize(), enrichment_size);
  DuffyH1ReferenceTable duffy_table;
  std::unique_ptr<mfem::ND_FECollection> nd_collection;
  std::unique_ptr<AffineCommonNDMassReferenceTable> affine_reference_table;
  if (options.fixed_subdivision)
  {
    nd_collection = std::make_unique<mfem::ND_FECollection>(h1_fespace.GetMaxElementOrder(),
                                                            mesh->Dimension());
    affine_reference_table = std::make_unique<AffineCommonNDMassReferenceTable>(options);
  }

  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element];
    if (element_dofs.h1.empty())
    {
      continue;
    }
    const auto enrichment_dofs =
        GetElementEnrichmentDofs(element_dofs.h1, topology.h1_dofs.size());

    auto &transformation = *mesh->GetElementTransformation(element);
    ElementH1EnrichmentMatrices matrices;
    try
    {
      const auto &h1_fe = *h1_fespace.GetFE(element);
      const auto *nd_fe = nd_collection
                              ? nd_collection->FiniteElementForGeometry(h1_fe.GetGeomType())
                              : nullptr;
      if (nd_fe && nd_fe->GetOrder() != h1_fe.GetOrder())
      {
        nd_fe = nullptr;
      }
      matrices = AssembleElementH1EnrichmentMatricesImpl(
          element_dofs, h1_fe, transformation, options,
          (!options.fixed_subdivision && h1_fe.GetOrder() == 1) ? &duffy_table : nullptr,
          affine_reference_table.get(), nd_fe);
    }
    catch (const std::exception &error)
    {
      throw std::runtime_error(
          fmt::format("Singular H1 sparse assembly failed on local mesh element {}: {}",
                      element, error.what()));
    }

    const double electric = materials[element].electric;
    ValidateSymmetricMatrixForScaling(
        matrices.enrichment_enrichment,
        matrices.enrichment_enrichment_estimated_absolute_error, electric);
    ValidateCouplingMatricesForScaling(
        matrices.standard_enrichment, matrices.enrichment_standard,
        matrices.standard_enrichment_estimated_absolute_error, electric);
    ScaleSymmetricMatrices(matrices.enrichment_enrichment,
                           matrices.enrichment_enrichment_estimated_absolute_error,
                           electric);
    ScaleCouplingMatrices(matrices.standard_enrichment, matrices.enrichment_standard,
                          matrices.standard_enrichment_estimated_absolute_error, electric);

    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    h1_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    ApplyStandardRowTransformation(dof_transformation, matrices.standard_enrichment,
                                   matrices.enrichment_standard,
                                   matrices.standard_enrichment_estimated_absolute_error);
    const auto unsigned_standard_dofs = UnsignedDofs(standard_dofs);

    result.diffusion.enrichment_enrichment->AddSubMatrix(enrichment_dofs, enrichment_dofs,
                                                         matrices.enrichment_enrichment);
    result.diffusion.standard_enrichment->AddSubMatrix(standard_dofs, enrichment_dofs,
                                                       matrices.standard_enrichment);
    result.diffusion.enrichment_enrichment_estimated_absolute_error->AddSubMatrix(
        enrichment_dofs, enrichment_dofs,
        matrices.enrichment_enrichment_estimated_absolute_error);
    result.diffusion.standard_enrichment_estimated_absolute_error->AddSubMatrix(
        unsigned_standard_dofs, enrichment_dofs,
        matrices.standard_enrichment_estimated_absolute_error);

    AccumulateQuadratureStatistics(matrices, result.total_quadrature_leaf_count,
                                   result.maximum_subdivision_depth);
  }

  FinalizeLocalSparseBlock(result.diffusion);
  result.duffy_reference_table_entries = duffy_table.Size();
  result.duffy_reference_cache_hits = duffy_table.Hits();
  if (affine_reference_table)
  {
    result.affine_reference_table_entries = affine_reference_table->Size();
    result.affine_reference_cache_hits = affine_reference_table->Hits();
    result.affine_reference_pattern_count = affine_reference_table->PatternCount();
    result.affine_reference_pattern_hits = affine_reference_table->PatternHits();
    result.affine_reference_persistent_hits = affine_reference_table->PersistentHits();
    result.affine_reference_persistent_writes = affine_reference_table->PersistentWrites();
    result.affine_reference_generated_leaf_count =
        affine_reference_table->GeneratedLeafCount();
    result.affine_reference_generation_time = affine_reference_table->GenerationTime();
  }
  return result;
}

LocalSparseH1EnrichmentMatrices AssembleLocalSparseH1EnrichmentMatrices(
    const TriangleDofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options)
{
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh)
  {
    throw std::invalid_argument("Triangular singular H1 sparse assembly requires a mesh!");
  }

  // The full triangular topology also contains rotational ND functions. H1 diffusion
  // needs only the gradient ND functions paired one-to-one with the H1 enrichment.
  TriangleDofTopology gradient_topology;
  gradient_topology.h1_dofs = topology.h1_dofs;
  gradient_topology.nd_dofs = topology.h1_dofs;
  gradient_topology.h1_to_nd.resize(topology.h1_dofs.size());
  std::iota(gradient_topology.h1_to_nd.begin(), gradient_topology.h1_to_nd.end(), 0);
  gradient_topology.elements.resize(topology.elements.size());
  for (std::size_t element = 0; element < topology.elements.size(); element++)
  {
    gradient_topology.elements[element].h1 = topology.elements[element].h1;
    gradient_topology.elements[element].nd.reserve(topology.elements[element].h1.size());
    for (const auto &dof : topology.elements[element].h1)
    {
      gradient_topology.elements[element].nd.push_back({dof.dof, dof.basis});
    }
  }

  mfem::ND_FECollection nd_collection(h1_fespace.GetMaxElementOrder(), mesh->Dimension());
  mfem::FiniteElementSpace nd_fespace(mesh, &nd_collection);
  auto full = AssembleLocalSparseEnrichmentMatrices(gradient_topology, h1_fespace,
                                                    nd_fespace, materials, options);
  LocalSparseH1EnrichmentMatrices result;
  result.diffusion = std::move(full.h1_diffusion);
  result.total_quadrature_leaf_count = full.total_quadrature_leaf_count;
  result.maximum_subdivision_depth = full.maximum_subdivision_depth;
  return result;
}

namespace
{

double EvaluateH1InterpolationPotential(const BarycentricPoint &lambda,
                                        const HigherOrderBasis &basis)
{
  double rho;
  switch (basis.family)
  {
    case HigherOrderBasisFamily::NODE_GRADIENT:
      rho = 1.0 - lambda[basis.nodes[0]];
      break;
    case HigherOrderBasisFamily::EDGE_GRADIENT:
      rho = 1.0 - lambda[basis.nodes[0]] - lambda[basis.nodes[1]];
      break;
    case HigherOrderBasisFamily::NODE_ROTATIONAL:
    case HigherOrderBasisFamily::EDGE_ROTATIONAL:
      throw std::invalid_argument(
          "An H1 singular interpolant requires a gradient-family basis!");
    default:
      throw std::invalid_argument("Unknown singular H1 basis family!");
  }
  constexpr double endpoint_tolerance = 128.0 * std::numeric_limits<double>::epsilon();
  if (rho <= endpoint_tolerance)
  {
    // Both node- and edge-gradient scalar potentials have the continuous
    // limiting behavior O(rho^nu) on their singular feature.
    return 0.0;
  }
  return EvaluateHigherOrderGradientPotential(lambda, basis);
}

double EvaluateH1InterpolationPotential(const TriangleBarycentricPoint &lambda,
                                        const TriangleBasis &basis)
{
  if (basis.family != HigherOrderBasisFamily::NODE_GRADIENT)
  {
    throw std::invalid_argument(
        "A triangular H1 singular interpolant requires a node-gradient basis!");
  }
  const double rho = 1.0 - lambda[basis.nodes[0]];
  constexpr double endpoint_tolerance = 128.0 * std::numeric_limits<double>::epsilon();
  if (rho <= endpoint_tolerance)
  {
    return 0.0;
  }
  return EvaluateTriangleNodeGradientPotential(lambda, basis.nodes[0], basis.nodes[1],
                                               basis.nu);
}

template <typename Topology, typename BarycentricPointBuilder>
mfem::SparseMatrix
BuildLocalH1EnrichmentInterpolant(const Topology &topology,
                                  mfem::FiniteElementSpace &h1_fespace,
                                  BarycentricPointBuilder &&make_barycentric_point)
{
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Singular H1 interpolation topology and finite-element space are inconsistent!");
  }

  using Entry = std::pair<int, int>;
  std::map<Entry, double> entries;
  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element].h1;
    if (element_dofs.empty())
    {
      continue;
    }
    const auto &finite_element = *h1_fespace.GetFE(element);
    const auto &nodes = finite_element.GetNodes();
    if (nodes.GetNPoints() != finite_element.GetDof())
    {
      throw std::runtime_error(
          "Singular H1 interpolation requires a nodal standard finite element!");
    }

    mfem::DenseMatrix values(finite_element.GetDof(),
                             static_cast<int>(element_dofs.size()));
    for (int row = 0; row < finite_element.GetDof(); row++)
    {
      const auto lambda = make_barycentric_point(nodes.IntPoint(row));
      for (int column = 0; column < static_cast<int>(element_dofs.size()); column++)
      {
        values(row, column) =
            EvaluateH1InterpolationPotential(lambda, element_dofs[column].basis);
      }
    }

    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    h1_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    if (standard_dofs.Size() != values.Height())
    {
      throw std::runtime_error(
          "Singular H1 interpolation received inconsistent standard element DOFs!");
    }
    dof_transformation.TransformPrimalCols(values);
    for (int row = 0; row < standard_dofs.Size(); row++)
    {
      const int signed_dof = standard_dofs[row];
      if (signed_dof == std::numeric_limits<int>::min())
      {
        throw std::runtime_error(
            "Singular H1 interpolation received an invalid signed standard DOF!");
      }
      const int standard_dof = signed_dof >= 0 ? signed_dof : -1 - signed_dof;
      const double sign = signed_dof >= 0 ? 1.0 : -1.0;
      for (int column = 0; column < static_cast<int>(element_dofs.size()); column++)
      {
        const std::size_t enrichment = element_dofs[column].dof;
        if (enrichment >= topology.h1_dofs.size())
        {
          throw std::runtime_error(
              "Singular H1 interpolation received an invalid enrichment DOF!");
        }
        const double value = sign * values(row, column);
        const auto [entry, inserted] =
            entries.emplace(Entry{standard_dof, static_cast<int>(enrichment)}, value);
        if (!inserted)
        {
          const double scale = std::max({1.0, std::abs(value), std::abs(entry->second)});
          if (std::abs(value - entry->second) >
              4096.0 * std::numeric_limits<double>::epsilon() * scale)
          {
            throw std::runtime_error(fmt::format(
                "Singular H1 interpolant is not conforming on local standard DOF {} "
                "for enrichment DOF {}: values {:.17g} and {:.17g} disagree!",
                standard_dof, enrichment, entry->second, value));
          }
        }
      }
    }
  }

  mfem::SparseMatrix interpolation(h1_fespace.GetVSize(),
                                   static_cast<int>(topology.h1_dofs.size()));
  for (const auto &[indices, value] : entries)
  {
    if (value != 0.0)
    {
      interpolation.Set(indices.first, indices.second, value);
    }
  }
  interpolation.Finalize();
  return interpolation;
}

class TetrahedronRotationalBasisCoefficient : public mfem::VectorCoefficient
{
private:
  HigherOrderBasis basis;

public:
  explicit TetrahedronRotationalBasisCoefficient(const HigherOrderBasis &basis)
    : mfem::VectorCoefficient(3), basis(basis)
  {
    if (IsGradientFamily(basis.family))
    {
      throw std::invalid_argument(
          "A rotational Nedelec interpolant received a gradient basis!");
    }
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    double jacobian_determinant = 0.0;
    const auto grad_lambda =
        GetBarycentricGradients(transformation, point, jacobian_determinant);
    const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                  point.z};
    const auto evaluated = EvaluateHigherOrderBasisValue(lambda, grad_lambda, basis);
    value.SetSize(3);
    for (int d = 0; d < 3; d++)
    {
      value[d] = evaluated[d];
    }
  }
};

class TriangleRotationalBasisCoefficient : public mfem::VectorCoefficient
{
private:
  TriangleBasis basis;

public:
  explicit TriangleRotationalBasisCoefficient(const TriangleBasis &basis)
    : mfem::VectorCoefficient(2), basis(basis)
  {
    if (IsGradientFamily(basis.family))
    {
      throw std::invalid_argument(
          "A rotational triangular Nedelec interpolant received a gradient basis!");
    }
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    double jacobian_determinant = 0.0;
    const auto grad_lambda =
        GetTriangleBarycentricGradients(transformation, point, jacobian_determinant);
    const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
    const auto evaluated = EvaluateTriangleBasisValue(lambda, grad_lambda, basis);
    value.SetSize(2);
    for (int d = 0; d < 2; d++)
    {
      value[d] = evaluated[d];
    }
  }
};

template <typename Topology, typename ProjectBasis>
mfem::SparseMatrix
BuildLocalNDRotationalEnrichmentInterpolant(const Topology &topology,
                                            mfem::FiniteElementSpace &nd_fespace,
                                            ProjectBasis &&project_basis)
{
  auto *mesh = nd_fespace.GetMesh();
  if (!mesh || topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Singular rotational interpolation topology and Nedelec space are "
        "inconsistent!");
  }

  using Entry = std::pair<int, int>;
  std::map<Entry, double> entries;
  for (int element = 0; element < mesh->GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element].nd;
    if (element_dofs.empty())
    {
      continue;
    }
    const auto &finite_element = *nd_fespace.GetFE(element);
    if (finite_element.GetRangeType() != mfem::FiniteElement::VECTOR ||
        finite_element.GetMapType() != mfem::FiniteElement::H_CURL)
    {
      throw std::runtime_error(
          "Singular rotational interpolation requires a Nedelec finite element!");
    }
    auto *transformation = mesh->GetElementTransformation(element);
    if (!transformation)
    {
      throw std::runtime_error(
          "Singular rotational interpolation could not obtain an element "
          "transformation!");
    }

    mfem::DenseMatrix values(finite_element.GetDof(),
                             static_cast<int>(element_dofs.size()));
    values = 0.0;
    for (int column = 0; column < static_cast<int>(element_dofs.size()); column++)
    {
      const auto &basis = element_dofs[column].basis;
      if (IsGradientFamily(basis.family))
      {
        continue;
      }
      mfem::Vector projected(values.GetColumn(column), values.Height());
      project_basis(basis, finite_element, *transformation, projected);
    }

    mfem::Array<int> standard_dofs;
    mfem::DofTransformation dof_transformation;
    nd_fespace.GetElementVDofs(element, standard_dofs, dof_transformation);
    if (standard_dofs.Size() != values.Height())
    {
      throw std::runtime_error(
          "Singular rotational interpolation received inconsistent standard element "
          "DOFs!");
    }
    dof_transformation.TransformPrimalCols(values);
    for (int row = 0; row < standard_dofs.Size(); row++)
    {
      const int signed_dof = standard_dofs[row];
      if (signed_dof == std::numeric_limits<int>::min())
      {
        throw std::runtime_error(
            "Singular rotational interpolation received an invalid signed standard "
            "DOF!");
      }
      const int standard_dof = signed_dof >= 0 ? signed_dof : -1 - signed_dof;
      const double sign = signed_dof >= 0 ? 1.0 : -1.0;
      for (int column = 0; column < static_cast<int>(element_dofs.size()); column++)
      {
        const std::size_t enrichment = element_dofs[column].dof;
        if (enrichment >= topology.nd_dofs.size())
        {
          throw std::runtime_error(
              "Singular rotational interpolation received an invalid enrichment DOF!");
        }
        const double value = sign * values(row, column);
        const auto [entry, inserted] =
            entries.emplace(Entry{standard_dof, static_cast<int>(enrichment)}, value);
        if (!inserted)
        {
          const double scale = std::max({1.0, std::abs(value), std::abs(entry->second)});
          if (std::abs(value - entry->second) > 2.0e-11 * scale)
          {
            throw std::runtime_error(fmt::format(
                "Singular rotational interpolant is not conforming on local standard "
                "DOF {} for enrichment DOF {}: values {:.17g} and {:.17g} disagree!",
                standard_dof, enrichment, entry->second, value));
          }
        }
      }
    }
  }

  mfem::SparseMatrix interpolation(nd_fespace.GetVSize(),
                                   static_cast<int>(topology.nd_dofs.size()));
  for (const auto &[indices, value] : entries)
  {
    if (value != 0.0)
    {
      interpolation.Set(indices.first, indices.second, value);
    }
  }
  interpolation.Finalize();
  return interpolation;
}

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichmentInterpolantImpl(const mfem::SparseMatrix &local,
                                       const TrueDofMap &enrichment_numbering,
                                       const mfem::ParFiniteElementSpace &standard_fespace)
{
  const MPI_Comm comm = standard_fespace.GetComm();
  const auto *standard_prolongation = standard_fespace.Dof_TrueDof_Matrix();
  if (!standard_prolongation || local.Height() != standard_prolongation->Height() ||
      local.Width() != static_cast<int>(enrichment_numbering.local_size))
  {
    throw std::invalid_argument(
        "Cannot parallelize an inconsistent singular enrichment interpolant!");
  }
  auto enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, enrichment_numbering);
  auto interpolation =
      AssembleParallelRAP(comm, local, *standard_prolongation, *enrichment_prolongation);

  Vector local_ones(standard_prolongation->Height());
  Vector true_multiplicity(standard_prolongation->Width());
  local_ones = 1.0;
  standard_prolongation->AbsMultTranspose(local_ones, true_multiplicity);
  for (int i = 0; i < true_multiplicity.Size(); i++)
  {
    if (!(true_multiplicity[i] > 0.0) || !std::isfinite(true_multiplicity[i]))
    {
      throw std::runtime_error(
          "Singular enrichment interpolation found an invalid standard true-DOF "
          "multiplicity!");
    }
  }
  interpolation->InvScaleRows(true_multiplicity);
  return interpolation;
}

}  // namespace

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelH1EnrichmentInterpolant(const DofTopology &topology,
                                     const ParallelDofNumbering &parallel_numbering,
                                     const mfem::ParFiniteElementSpace &h1_fespace)
{
  auto local = BuildLocalH1EnrichmentInterpolant(
      topology, const_cast<mfem::ParFiniteElementSpace &>(h1_fespace),
      [](const mfem::IntegrationPoint &point)
      {
        return BarycentricPoint{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
      });
  return BuildParallelEnrichmentInterpolantImpl(local, parallel_numbering.h1, h1_fespace);
}

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelH1EnrichmentInterpolant(const TriangleDofTopology &topology,
                                     const ParallelDofNumbering &parallel_numbering,
                                     const mfem::ParFiniteElementSpace &h1_fespace)
{
  auto local = BuildLocalH1EnrichmentInterpolant(
      topology, const_cast<mfem::ParFiniteElementSpace &>(h1_fespace),
      [](const mfem::IntegrationPoint &point)
      { return TriangleBarycentricPoint{1.0 - point.x - point.y, point.x, point.y}; });
  return BuildParallelEnrichmentInterpolantImpl(local, parallel_numbering.h1, h1_fespace);
}

std::unique_ptr<mfem::HypreParMatrix> BuildParallelNDRotationalEnrichmentInterpolant(
    const DofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    const mfem::ParFiniteElementSpace &nd_fespace)
{
  auto local = BuildLocalNDRotationalEnrichmentInterpolant(
      topology, const_cast<mfem::ParFiniteElementSpace &>(nd_fespace),
      [](const HigherOrderBasis &basis, const mfem::FiniteElement &finite_element,
         mfem::ElementTransformation &transformation, mfem::Vector &projected)
      {
        TetrahedronRotationalBasisCoefficient coefficient(basis);
        finite_element.Project(coefficient, transformation, projected);
      });
  return BuildParallelEnrichmentInterpolantImpl(local, parallel_numbering.nd, nd_fespace);
}

std::unique_ptr<mfem::HypreParMatrix> BuildParallelNDRotationalEnrichmentInterpolant(
    const TriangleDofTopology &topology, const ParallelDofNumbering &parallel_numbering,
    const mfem::ParFiniteElementSpace &nd_fespace)
{
  auto local = BuildLocalNDRotationalEnrichmentInterpolant(
      topology, const_cast<mfem::ParFiniteElementSpace &>(nd_fespace),
      [](const TriangleBasis &basis, const mfem::FiniteElement &finite_element,
         mfem::ElementTransformation &transformation, mfem::Vector &projected)
      {
        TriangleRotationalBasisCoefficient coefficient(basis);
        finite_element.Project(coefficient, transformation, projected);
      });
  return BuildParallelEnrichmentInterpolantImpl(local, parallel_numbering.nd, nd_fespace);
}

ParallelSparseEnrichmentMatrices
AssembleParallelSparseEnrichmentMatrices(const LocalSparseEnrichmentMatrices &local,
                                         const ParallelDofNumbering &parallel_numbering,
                                         const mfem::ParFiniteElementSpace &h1_fespace,
                                         const mfem::ParFiniteElementSpace &nd_fespace)
{
  if (h1_fespace.GetComm() != nd_fespace.GetComm() ||
      h1_fespace.GetParMesh() != nd_fespace.GetParMesh())
  {
    throw std::invalid_argument(
        "Singular parallel sparse assembly requires spaces on one communicator and mesh!");
  }
  const MPI_Comm comm = h1_fespace.GetComm();
  auto h1_enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, parallel_numbering.h1);
  auto nd_enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, parallel_numbering.nd);
  const auto *h1_standard_prolongation = h1_fespace.Dof_TrueDof_Matrix();
  const auto *nd_standard_prolongation = nd_fespace.Dof_TrueDof_Matrix();
  if (!h1_standard_prolongation || !nd_standard_prolongation)
  {
    throw std::runtime_error(
        "MFEM did not provide standard finite-element prolongation matrices!");
  }
  auto h1_absolute_prolongation = BuildAbsoluteProlongation(*h1_standard_prolongation);
  auto nd_absolute_prolongation = BuildAbsoluteProlongation(*nd_standard_prolongation);

  ParallelSparseEnrichmentMatrices result;
  AssembleParallelSparseBlock(comm, local.h1_diffusion, *h1_standard_prolongation,
                              *h1_absolute_prolongation, *h1_enrichment_prolongation,
                              result.h1_diffusion);
  AssembleParallelSparseBlock(comm, local.h1_mass, *h1_standard_prolongation,
                              *h1_absolute_prolongation, *h1_enrichment_prolongation,
                              result.h1_mass);
  AssembleParallelSparseBlock(comm, local.nd_mass, *nd_standard_prolongation,
                              *nd_absolute_prolongation, *nd_enrichment_prolongation,
                              result.nd_mass);
  AssembleParallelSparseBlock(comm, local.nd_curl_curl, *nd_standard_prolongation,
                              *nd_absolute_prolongation, *nd_enrichment_prolongation,
                              result.nd_curl_curl);
  return result;
}

ParallelSparseOperatorBlocks
AssembleParallelSparseH1EnrichmentMatrices(const LocalSparseH1EnrichmentMatrices &local,
                                           const ParallelDofNumbering &parallel_numbering,
                                           const mfem::ParFiniteElementSpace &h1_fespace)
{
  const MPI_Comm comm = h1_fespace.GetComm();
  auto enrichment_prolongation =
      BuildParallelEnrichmentProlongation(comm, parallel_numbering.h1);
  const auto *standard_prolongation = h1_fespace.Dof_TrueDof_Matrix();
  if (!standard_prolongation)
  {
    throw std::runtime_error(
        "MFEM did not provide an H1 finite-element prolongation matrix!");
  }
  auto absolute_prolongation = BuildAbsoluteProlongation(*standard_prolongation);

  ParallelSparseOperatorBlocks result;
  AssembleParallelSparseBlock(comm, local.diffusion, *standard_prolongation,
                              *absolute_prolongation, *enrichment_prolongation, result);
  return result;
}

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichmentGradient(MPI_Comm comm,
                                const ParallelDofNumbering &parallel_numbering)
{
  ValidateTrueDofMap(comm, parallel_numbering.h1);
  ValidateTrueDofMap(comm, parallel_numbering.nd);
  const auto &h1 = parallel_numbering.h1;
  const auto &nd = parallel_numbering.nd;
  bool locally_valid =
      parallel_numbering.h1_to_nd_true.size() == h1.local_to_true.size() &&
      nd.owned_size <= static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max());
  std::vector<HYPRE_BigInt> row_columns(static_cast<std::size_t>(nd.owned_size), -1);
  if (locally_valid)
  {
    for (std::size_t local_h1 = 0; local_h1 < h1.local_to_true.size(); local_h1++)
    {
      if (h1.owner[local_h1] != Mpi::Rank(comm))
      {
        continue;
      }
      const HYPRE_BigInt nd_true = parallel_numbering.h1_to_nd_true[local_h1];
      if (nd_true < nd.owned_offset || nd_true >= nd.owned_offset + nd.owned_size)
      {
        locally_valid = false;
        break;
      }
      auto &column = row_columns[nd_true - nd.owned_offset];
      if (column >= 0)
      {
        locally_valid = false;
        break;
      }
      column = h1.local_to_true[local_h1];
    }
  }
  Mpi::GlobalAnd(1, &locally_valid, comm);
  if (!locally_valid)
  {
    throw std::invalid_argument("Parallel singular H1-to-ND gradient map is inconsistent!");
  }

  const int local_nd_size = static_cast<int>(nd.owned_size);
  std::vector<int> rows(local_nd_size + 1);
  std::vector<HYPRE_BigInt> columns;
  std::vector<double> values;
  columns.reserve(h1.owned_size);
  values.reserve(h1.owned_size);
  for (int row = 0; row < local_nd_size; row++)
  {
    rows[row] = static_cast<int>(columns.size());
    if (row_columns[row] >= 0)
    {
      columns.push_back(row_columns[row]);
      values.push_back(1.0);
    }
  }
  rows[local_nd_size] = static_cast<int>(columns.size());
  if (columns.size() != static_cast<std::size_t>(h1.owned_size))
  {
    throw std::logic_error("Parallel singular gradient omitted an owned H1 true DOF!");
  }

  const auto nd_partition =
      BuildPartition(comm, nd.owned_offset, nd.owned_size, nd.global_size);
  const auto h1_partition =
      BuildPartition(comm, h1.owned_offset, h1.owned_size, h1.global_size);
  return std::make_unique<mfem::HypreParMatrix>(
      comm, local_nd_size, nd.global_size, h1.global_size, rows.data(), columns.data(),
      values.data(), nd_partition.data(), h1_partition.data());
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
