// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularassembly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <vector>
#include <fmt/format.h>

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

void ValidateAdaptiveAssemblyOptions(const AdaptiveAssemblyOptions &options)
{
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Singular element assembly has invalid adaptive quadrature options!");
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
TriangleQuadratureResult IntegrateTriangleDuffy(int row_node, int column_node,
                                                const AdaptiveAssemblyOptions &options,
                                                Integrand &&integrand)
{
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
      "{:.17g}, orders = {}/{}, order refinements = {}!",
      value, comparison, error, tolerance, high_order, comparison_order,
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

ElementEnrichmentMatrices AssembleElementEnrichmentMatrices(
    const ElementDofMap &element_dofs, const BarycentricGradients &grad_lambda,
    double jacobian_determinant, const AdaptiveAssemblyOptions &options)
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
      const auto reference = ComputeAdaptiveReferenceIntegral(
          ReferenceBasis{row_basis}, ReferenceBasis{column_basis}, options.quadrature_order,
          options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions);
      if (!reference.converged)
      {
        throw std::runtime_error(
            "Singular element adaptive reference tensor did not converge!");
      }
      RecordQuadratureStatistics(reference, result);

      const double mass =
          ContractMass(reference.integral, grad_lambda, jacobian_determinant);
      const double mass_error =
          ContractMassError(reference, grad_lambda, jacobian_determinant);
      const bool curl_free =
          IsGradientFamily(row_basis.family) || IsGradientFamily(column_basis.family);
      const double curl_curl = curl_free ? 0.0
                                         : ContractCurlCurl(reference.integral, grad_lambda,
                                                            jacobian_determinant);
      const double curl_curl_error =
          curl_free ? 0.0
                    : ContractCurlCurlError(reference, grad_lambda, jacobian_determinant);
      if (!std::isfinite(mass) || !std::isfinite(mass_error) || !std::isfinite(curl_curl) ||
          !std::isfinite(curl_curl_error) || mass_error < 0.0 || curl_curl_error < 0.0)
      {
        throw std::runtime_error(
            "Singular element assembly produced a nonfinite matrix entry!");
      }

      result.nd_mass(row, column) = result.nd_mass(column, row) = mass;
      result.nd_mass_estimated_absolute_error(row, column) =
          result.nd_mass_estimated_absolute_error(column, row) = mass_error;
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
      const auto integral = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
          [&](const BarycentricPoint &lambda)
          {
            return jacobian_determinant *
                   EvaluateHigherOrderGradientPotential(lambda, row_basis) *
                   EvaluateHigherOrderGradientPotential(lambda, column_basis);
          });
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
      const auto mass = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      const auto curl_curl = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      const auto mass = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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

ElementEnrichmentMatrices
AssembleTriangleElementEnrichmentMatrices(const TriangleElementDofMap &element_dofs,
                                          const TriangleBarycentricGradients &grad_lambda,
                                          double jacobian_determinant,
                                          const AdaptiveAssemblyOptions &options)
{
  ValidateTriangleInputs(element_dofs, grad_lambda, jacobian_determinant, options);
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

  for (int row = 0; row < nd_size; row++)
  {
    for (int column = row; column < nd_size; column++)
    {
      const auto &row_basis = element_dofs.nd[row].basis;
      const auto &column_basis = element_dofs.nd[column].basis;
      const int row_node = TriangleSingularNode(row_basis);
      const int column_node = TriangleSingularNode(column_basis);
      const auto mass = IntegrateTriangleDuffy(
          row_node, column_node, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
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
          row_node, column_node, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
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
          TriangleSingularNode(row_basis), TriangleSingularNode(column_basis), options,
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
          row_node, column_node, options,
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
          row_node, column_node, options,
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
          TriangleSingularNode(row_basis), TriangleSingularNode(column_basis), options,
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

ElementStandardEnrichmentMatrices AssembleElementStandardEnrichmentMatrices(
    const ElementDofMap &element_dofs, const mfem::FiniteElement &h1_fe,
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    const AdaptiveAssemblyOptions &options)
{
  ValidateStandardFiniteElements(h1_fe, nd_fe, transformation);
  const bool affine_geometry = IsAffineElementTransformation(transformation);
  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  double center_jacobian_determinant;
  const auto center_grad_lambda =
      GetBarycentricGradients(transformation, center, center_jacobian_determinant);
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
  for (int standard = 0; standard < standard_nd_size; standard++)
  {
    for (int enrichment = 0; enrichment < enrichment_nd_size; enrichment++)
    {
      const auto &basis = element_dofs.nd[enrichment].basis;
      const auto mass = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      ValidateAdaptiveEntry(mass, "mass", standard, enrichment, basis);
      RecordQuadratureStatistics(mass, result);
      result.nd_mass_standard_enrichment(standard, enrichment) = mass.value;
      result.nd_mass_enrichment_standard(enrichment, standard) = mass.value;
      result.nd_mass_estimated_absolute_error(standard, enrichment) =
          mass.estimated_absolute_error;

      if (IsGradientFamily(basis.family))
      {
        result.nd_curl_curl_standard_enrichment(standard, enrichment) = 0.0;
        result.nd_curl_curl_enrichment_standard(enrichment, standard) = 0.0;
        result.nd_curl_curl_estimated_absolute_error(standard, enrichment) = 0.0;
        continue;
      }

      const auto curl_curl = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
          [&](const BarycentricPoint &lambda)
          {
            point.Set3(lambda[1], lambda[2], lambda[3]);
            double jacobian_determinant = center_jacobian_determinant;
            auto grad_lambda = center_grad_lambda;
            if (affine_geometry)
            {
              CalcAffinePhysCurlShape(nd_fe, point, center_jacobian,
                                      center_jacobian_determinant, reference_standard_curl,
                                      standard_curl);
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
      ValidateAdaptiveEntry(curl_curl, "curl-curl", standard, enrichment, basis);
      RecordQuadratureStatistics(curl_curl, result);
      result.nd_curl_curl_standard_enrichment(standard, enrichment) = curl_curl.value;
      result.nd_curl_curl_enrichment_standard(enrichment, standard) = curl_curl.value;
      result.nd_curl_curl_estimated_absolute_error(standard, enrichment) =
          curl_curl.estimated_absolute_error;
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

  mfem::Vector standard_shape(standard_h1_size);
  for (int standard_h1 = 0; standard_h1 < standard_h1_size; standard_h1++)
  {
    for (int enrichment_h1 = 0; enrichment_h1 < enrichment_h1_size; enrichment_h1++)
    {
      const auto &basis = element_dofs.h1[enrichment_h1].basis;
      const auto mass = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      ValidateAdaptiveEntry(mass, "H1 mass", standard_h1, enrichment_h1, basis);
      RecordQuadratureStatistics(mass, result);
      result.h1_mass_standard_enrichment(standard_h1, enrichment_h1) = mass.value;
      result.h1_mass_enrichment_standard(enrichment_h1, standard_h1) = mass.value;
      result.h1_mass_estimated_absolute_error(standard_h1, enrichment_h1) =
          mass.estimated_absolute_error;
    }
  }
  return result;
}

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
          singular_node, singular_node, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant = center_jacobian_determinant;
            auto grad_lambda = center_grad_lambda;
            if (affine_geometry)
            {
              CalcAffinePhysVShape(nd_fe, point, center_inverse_jacobian,
                                   reference_standard_value, standard_value);
            }
            else
            {
              grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                            jacobian_determinant);
              nd_fe.CalcPhysVShape(transformation, standard_value);
            }
            const auto singular = EvaluateTriangleBasis(lambda, grad_lambda, basis);
            return jacobian_determinant *
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
          singular_node, singular_node, options,
          [&](const TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant = center_jacobian_determinant;
            auto grad_lambda = center_grad_lambda;
            if (affine_geometry)
            {
              CalcAffinePhysCurlShape(nd_fe, point, center_jacobian,
                                      center_jacobian_determinant, reference_standard_curl,
                                      standard_curl);
            }
            else
            {
              grad_lambda = GetTriangleBarycentricGradients(transformation, point,
                                                            jacobian_determinant);
              nd_fe.CalcPhysCurlShape(transformation, standard_curl);
            }
            const auto singular = EvaluateTriangleBasis(lambda, grad_lambda, basis);
            return jacobian_determinant * standard_curl(standard, 0) * singular.curl;
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
          singular_node, singular_node, options,
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
    DuffyH1ReferenceTable *duffy_table)
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
      const auto integral = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      const auto integral = IntegrateReferenceTetrahedronAdaptive(
          options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
          options.maximum_subdivisions,
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
      (h1_fe.GetOrder() == 1) ? &duffy_table : nullptr);
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

LocalSparseEnrichmentMatrices AssembleLocalSparseEnrichmentMatrices(
    const DofTopology &topology, mfem::FiniteElementSpace &h1_fespace,
    mfem::FiniteElementSpace &nd_fespace,
    const std::vector<IsotropicMaterialCoefficients> &materials,
    const AdaptiveAssemblyOptions &options)
{
  ValidateAdaptiveAssemblyOptions(options);
  auto *mesh = h1_fespace.GetMesh();
  if (!mesh || nd_fespace.GetMesh() != mesh ||
      topology.elements.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      materials.size() != static_cast<std::size_t>(mesh->GetNE()) ||
      topology.h1_to_nd.size() != topology.h1_dofs.size() ||
      topology.h1_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      topology.nd_dofs.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    throw std::invalid_argument(
        "Singular sparse assembly topology, spaces, and materials must share one mesh!");
  }
  for (std::size_t h1 = 0; h1 < topology.h1_to_nd.size(); h1++)
  {
    const std::size_t nd = topology.h1_to_nd[h1];
    if (nd >= topology.nd_dofs.size() || !(topology.h1_dofs[h1] == topology.nd_dofs[nd]))
    {
      throw std::invalid_argument(
          "Singular sparse assembly requires a consistent H1-to-ND topology!");
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
      enrichment = AssembleElementEnrichmentMatrices(element_dofs, transformation, options);
      coupling = AssembleElementStandardEnrichmentMatrices(
          element_dofs, *h1_fespace.GetFE(element), *nd_fespace.GetFE(element),
          transformation, options);
    }
    catch (const std::exception &error)
    {
      throw std::runtime_error(
          fmt::format("Singular sparse assembly failed on local mesh element {}: {}",
                      element, error.what()));
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
      matrices = AssembleElementH1EnrichmentMatricesImpl(
          element_dofs, h1_fe, transformation, options,
          (h1_fe.GetOrder() == 1) ? &duffy_table : nullptr);
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
  mfem::ND_FECollection nd_collection(h1_fespace.GetMaxElementOrder(), mesh->Dimension());
  mfem::FiniteElementSpace nd_fespace(mesh, &nd_collection);
  auto full = AssembleLocalSparseEnrichmentMatrices(topology, h1_fespace, nd_fespace,
                                                    materials, options);
  LocalSparseH1EnrichmentMatrices result;
  result.diffusion = std::move(full.h1_diffusion);
  result.total_quadrature_leaf_count = full.total_quadrature_leaf_count;
  result.maximum_subdivision_depth = full.maximum_subdivision_depth;
  return result;
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
