// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/singularassembly.hpp"
#include "fem/singularelements.hpp"

namespace palace
{

using namespace Catch::Matchers;
using fem::singular::BarycentricGradients;
using fem::singular::BarycentricPoint;
using fem::singular::TriangleBarycentricGradients;
using fem::singular::TriangleBarycentricPoint;
using fem::singular::TriangleVectorBasisValue;
using fem::singular::Vector2;
using fem::singular::Vector3;
using fem::singular::VectorBasisValue;

namespace
{

double Dot(const Vector3 &x, const Vector3 &y)
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

double Norm(const Vector3 &x)
{
  return std::sqrt(Dot(x, x));
}

double Dot(const Vector2 &x, const Vector2 &y)
{
  return x[0] * y[0] + x[1] * y[1];
}

double Norm(const Vector2 &x)
{
  return std::sqrt(Dot(x, x));
}

Vector3 Cross(const Vector3 &x, const Vector3 &y)
{
  return {x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
}

void Add(double scale, const Vector3 &x, Vector3 &y)
{
  for (int d = 0; d < 3; d++)
  {
    y[d] += scale * x[d];
  }
}

void CheckVector(const Vector3 &actual, const Vector3 &expected, double tolerance)
{
  for (int d = 0; d < 3; d++)
  {
    CAPTURE(d, actual[d], expected[d]);
    CHECK_THAT(actual[d], WithinAbs(expected[d], tolerance));
  }
}

void CheckVector(const Vector2 &actual, const Vector2 &expected, double tolerance)
{
  for (int d = 0; d < 2; d++)
  {
    CAPTURE(d, actual[d], expected[d]);
    CHECK_THAT(actual[d], WithinAbs(expected[d], tolerance));
  }
}

template <typename F>
Vector3 NumericalGradient(const Vector3 &point, F &&function, double step = 1.0e-6)
{
  Vector3 gradient{};
  for (int d = 0; d < 3; d++)
  {
    auto plus = point;
    auto minus = point;
    plus[d] += step;
    minus[d] -= step;
    gradient[d] = (function(plus) - function(minus)) / (2.0 * step);
  }
  return gradient;
}

template <typename F>
Vector2 NumericalGradient(const Vector2 &point, F &&function, double step = 1.0e-6)
{
  Vector2 gradient{};
  for (int d = 0; d < 2; d++)
  {
    auto plus = point;
    auto minus = point;
    plus[d] += step;
    minus[d] -= step;
    gradient[d] = (function(plus) - function(minus)) / (2.0 * step);
  }
  return gradient;
}

template <typename F>
Vector3 NumericalCurl(const Vector3 &point, F &&function, double step = 1.0e-6)
{
  std::array<Vector3, 3> derivative{};
  for (int d = 0; d < 3; d++)
  {
    auto plus = point;
    auto minus = point;
    plus[d] += step;
    minus[d] -= step;
    const auto value_plus = function(plus);
    const auto value_minus = function(minus);
    for (int c = 0; c < 3; c++)
    {
      derivative[d][c] = (value_plus[c] - value_minus[c]) / (2.0 * step);
    }
  }
  return {derivative[1][2] - derivative[2][1], derivative[2][0] - derivative[0][2],
          derivative[0][1] - derivative[1][0]};
}

template <typename F>
double NumericalCurl(const Vector2 &point, F &&function, double step = 1.0e-6)
{
  auto plus_x = point;
  auto minus_x = point;
  auto plus_y = point;
  auto minus_y = point;
  plus_x[0] += step;
  minus_x[0] -= step;
  plus_y[1] += step;
  minus_y[1] -= step;
  const double derivative_x = (function(plus_x)[1] - function(minus_x)[1]) / (2.0 * step);
  const double derivative_y = (function(plus_y)[0] - function(minus_y)[0]) / (2.0 * step);
  return derivative_x - derivative_y;
}

Vector3 TangentialPart(const Vector3 &value, const Vector3 &normal)
{
  Vector3 tangential = value;
  Add(-Dot(value, normal) / Dot(normal, normal), normal, tangential);
  return tangential;
}

BarycentricPoint FacePoint(int opposite)
{
  BarycentricPoint lambda{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
  lambda[opposite] = 0.0;
  return lambda;
}

double LogSlope(double rho_1, double value_1, double rho_2, double value_2)
{
  return std::log(value_2 / value_1) / std::log(rho_2 / rho_1);
}

VectorBasisValue EvaluatePolynomialWeightedEdge(const BarycentricPoint &lambda,
                                                const BarycentricGradients &grad_lambda,
                                                int i, int j, int k, int l)
{
  const auto edge = fem::singular::EvaluateStandardEdge(lambda, grad_lambda, k, l);
  const double polynomial = lambda[i] * lambda[j];
  Vector3 grad_polynomial{0.0, 0.0, 0.0};
  Add(lambda[j], grad_lambda[i], grad_polynomial);
  Add(lambda[i], grad_lambda[j], grad_polynomial);

  VectorBasisValue result{
      {polynomial * edge.value[0], polynomial * edge.value[1], polynomial * edge.value[2]},
      {polynomial * edge.curl[0], polynomial * edge.curl[1], polynomial * edge.curl[2]}};
  Add(1.0, Cross(grad_polynomial, edge.value), result.curl);
  return result;
}

double IntegrateRadialPowerOnLogGrid(double power, double lower, int intervals)
{
  const double u_min = std::log(lower);
  const double step = -u_min / intervals;
  double integral = 0.0;
  for (int n = 0; n < intervals; n++)
  {
    const double u = u_min + (n + 0.5) * step;
    integral += std::exp((power + 1.0) * u) * step;
  }
  return integral;
}

double ExactBarycentricMonomialIntegral(const std::array<int, 4> &exponents)
{
  int total = 0;
  double numerator = 1.0;
  for (int exponent : exponents)
  {
    total += exponent;
    numerator *= std::tgamma(exponent + 1.0);
  }
  return numerator / std::tgamma(total + 4.0);
}

double TensorDifferenceMaxNorm(const fem::singular::ReferenceCoefficientTensor &a,
                               const fem::singular::ReferenceCoefficientTensor &b)
{
  double norm = 0.0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      norm = std::max(norm, std::abs(a[i][j] - b[i][j]));
    }
  }
  return norm;
}

struct ReferenceTensorPair
{
  fem::singular::ReferenceCoefficientTensor mass;
  fem::singular::ReferenceCoefficientTensor curl_curl;
};

template <typename Q>
ReferenceTensorPair
DirectReferencePairIntegral(const fem::singular::FirstOrderBasis &row_basis,
                            const fem::singular::FirstOrderBasis &column_basis,
                            Q &&quadrature)
{
  const auto canonical =
      fem::singular::CanonicalizeFirstOrderBasisPair(row_basis, column_basis);
  std::array<std::array<long double, 3>, 3> mass{};
  std::array<std::array<long double, 3>, 3> curl_curl{};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  quadrature(
      [&](const BarycentricPoint &lambda, double weight)
      {
        const auto row = fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda,
                                                                canonical.row_basis);
        const auto column = fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda,
                                                                   canonical.column_basis);
        for (int u = 0; u < 3; u++)
        {
          for (int v = 0; v < 3; v++)
          {
            mass[u][v] += static_cast<long double>(weight) * row.value[u] * column.value[v];
            curl_curl[u][v] +=
                static_cast<long double>(weight) * row.curl[u] * column.curl[v];
          }
        }
      });

  ReferenceTensorPair result{};
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      result.mass[u][v] = static_cast<double>(mass[u][v]);
      result.curl_curl[u][v] = static_cast<double>(curl_curl[u][v]);
    }
  }
  return result;
}

template <typename Q>
ReferenceTensorPair
DirectReferencePairIntegral(const fem::singular::ReferenceBasis &row_basis,
                            const fem::singular::ReferenceBasis &column_basis,
                            Q &&quadrature)
{
  const auto canonical =
      fem::singular::CanonicalizeReferenceBasisPair(row_basis, column_basis);
  std::array<std::array<long double, 3>, 3> mass{};
  std::array<std::array<long double, 3>, 3> curl_curl{};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  quadrature(
      [&](const BarycentricPoint &lambda, double weight)
      {
        const auto row =
            fem::singular::EvaluateReferenceBasis(lambda, grad_lambda, canonical.row_basis);
        const auto column = fem::singular::EvaluateReferenceBasis(lambda, grad_lambda,
                                                                  canonical.column_basis);
        for (int u = 0; u < 3; u++)
        {
          for (int v = 0; v < 3; v++)
          {
            mass[u][v] += static_cast<long double>(weight) * row.value[u] * column.value[v];
            curl_curl[u][v] +=
                static_cast<long double>(weight) * row.curl[u] * column.curl[v];
          }
        }
      });

  ReferenceTensorPair result{};
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      result.mass[u][v] = static_cast<double>(mass[u][v]);
      result.curl_curl[u][v] = static_cast<double>(curl_curl[u][v]);
    }
  }
  return result;
}

template <typename Integral>
void CheckEntrywiseReferenceConvergence(const Integral &coarse, const Integral &fine,
                                        const ReferenceTensorPair &reference)
{
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      CAPTURE(u, v);
      const double coarse_mass_error = std::abs(coarse.mass[u][v] - reference.mass[u][v]);
      const double fine_mass_error = std::abs(fine.mass[u][v] - reference.mass[u][v]);
      const double coarse_curl_error =
          std::abs(coarse.curl_curl[u][v] - reference.curl_curl[u][v]);
      const double fine_curl_error =
          std::abs(fine.curl_curl[u][v] - reference.curl_curl[u][v]);
      CHECK(fine_mass_error <= std::max(coarse_mass_error, 5.0e-11));
      CHECK(fine_curl_error <= std::max(coarse_curl_error, 5.0e-11));
    }
  }
}

BarycentricGradients AffineBarycentricGradients(const Vector3 &edge_1,
                                                const Vector3 &edge_2,
                                                const Vector3 &edge_3,
                                                double &jacobian_determinant)
{
  const double determinant = Dot(edge_1, Cross(edge_2, edge_3));
  jacobian_determinant = std::abs(determinant);
  BarycentricGradients gradients;
  const auto cross_23 = Cross(edge_2, edge_3);
  const auto cross_31 = Cross(edge_3, edge_1);
  const auto cross_12 = Cross(edge_1, edge_2);
  for (int d = 0; d < 3; d++)
  {
    gradients[1][d] = cross_23[d] / determinant;
    gradients[2][d] = cross_31[d] / determinant;
    gradients[3][d] = cross_12[d] / determinant;
    gradients[0][d] = -gradients[1][d] - gradients[2][d] - gradients[3][d];
  }
  return gradients;
}

std::array<double, 2>
DirectAffinePairIntegral(const fem::singular::FirstOrderBasis &row_basis,
                         const fem::singular::FirstOrderBasis &column_basis,
                         const BarycentricGradients &grad_lambda,
                         double jacobian_determinant, int order, int subdivisions)
{
  const auto canonical =
      fem::singular::CanonicalizeFirstOrderBasisPair(row_basis, column_basis);
  std::array<long double, 2> values{};
  fem::singular::ForEachReferenceTetrahedronQuadraturePoint(
      order, subdivisions,
      [&](const BarycentricPoint &canonical_lambda, double weight)
      {
        BarycentricPoint lambda;
        for (int i = 0; i < 4; i++)
        {
          lambda[i] = canonical_lambda[canonical.input_to_canonical[i]];
        }
        const auto row =
            fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda, row_basis);
        const auto column =
            fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda, column_basis);
        values[0] += static_cast<long double>(jacobian_determinant) * weight *
                     Dot(row.value, column.value);
        values[1] += static_cast<long double>(jacobian_determinant) * weight *
                     Dot(row.curl, column.curl);
      });
  return {static_cast<double>(values[0]), static_cast<double>(values[1])};
}

std::array<double, 2>
DirectAffinePairIntegral(const fem::singular::ReferenceBasis &row_basis,
                         const fem::singular::ReferenceBasis &column_basis,
                         const BarycentricGradients &grad_lambda,
                         double jacobian_determinant, int order, int subdivisions)
{
  const auto canonical =
      fem::singular::CanonicalizeReferenceBasisPair(row_basis, column_basis);
  std::array<long double, 2> values{};
  fem::singular::ForEachReferenceTetrahedronQuadraturePoint(
      order, subdivisions,
      [&](const BarycentricPoint &canonical_lambda, double weight)
      {
        BarycentricPoint lambda;
        for (int i = 0; i < 4; i++)
        {
          lambda[i] = canonical_lambda[canonical.input_to_canonical[i]];
        }
        const auto row =
            fem::singular::EvaluateReferenceBasis(lambda, grad_lambda, row_basis);
        const auto column =
            fem::singular::EvaluateReferenceBasis(lambda, grad_lambda, column_basis);
        values[0] += static_cast<long double>(jacobian_determinant) * weight *
                     Dot(row.value, column.value);
        values[1] += static_cast<long double>(jacobian_determinant) * weight *
                     Dot(row.curl, column.curl);
      });
  return {static_cast<double>(values[0]), static_cast<double>(values[1])};
}

std::vector<double> SymmetricEigenvalues(mfem::DenseMatrix matrix)
{
  const int size = matrix.Height();
  const double tolerance = 1.0e-14 * std::max(1.0, matrix.MaxMaxNorm());
  for (int iteration = 0; iteration < 100 * size * size; iteration++)
  {
    int p = 0;
    int q = 0;
    double maximum = 0.0;
    for (int i = 0; i < size; i++)
    {
      for (int j = i + 1; j < size; j++)
      {
        if (std::abs(matrix(i, j)) > maximum)
        {
          maximum = std::abs(matrix(i, j));
          p = i;
          q = j;
        }
      }
    }
    if (maximum <= tolerance)
    {
      break;
    }

    const double off_diagonal = matrix(p, q);
    const double tau = (matrix(q, q) - matrix(p, p)) / (2.0 * off_diagonal);
    const double tangent =
        std::copysign(1.0 / (std::abs(tau) + std::sqrt(1.0 + tau * tau)), tau);
    const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
    const double sine = tangent * cosine;
    for (int k = 0; k < size; k++)
    {
      if (k == p || k == q)
      {
        continue;
      }
      const double value_kp = matrix(k, p);
      const double value_kq = matrix(k, q);
      matrix(k, p) = matrix(p, k) = cosine * value_kp - sine * value_kq;
      matrix(k, q) = matrix(q, k) = sine * value_kp + cosine * value_kq;
    }
    matrix(p, p) -= tangent * off_diagonal;
    matrix(q, q) += tangent * off_diagonal;
    matrix(p, q) = matrix(q, p) = 0.0;
  }

  std::vector<double> eigenvalues(size);
  for (int i = 0; i < size; i++)
  {
    eigenvalues[i] = matrix(i, i);
  }
  return eigenvalues;
}

void CheckBasisValue(const VectorBasisValue &actual, const VectorBasisValue &expected,
                     double tolerance)
{
  CheckVector(actual.value, expected.value, tolerance);
  CheckVector(actual.curl, expected.curl, tolerance);
}

double RadicalInverse(int index, int base)
{
  double value = 0.0;
  double scale = 1.0 / base;
  while (index > 0)
  {
    value += scale * (index % base);
    index /= base;
    scale /= base;
  }
  return value;
}

std::vector<BarycentricPoint> InteriorSamplePoints(int count)
{
  constexpr std::array<int, 4> bases{2, 3, 5, 7};
  std::vector<BarycentricPoint> points;
  points.reserve(count);
  for (int sample = 1; sample <= count; sample++)
  {
    BarycentricPoint lambda;
    double sum = 0.0;
    for (int i = 0; i < 4; i++)
    {
      lambda[i] = -std::log(RadicalInverse(sample, bases[i]));
      sum += lambda[i];
    }
    for (double &value : lambda)
    {
      value /= sum;
    }
    points.push_back(lambda);
  }
  return points;
}

std::vector<double>
SampledBasisSingularValues(const std::vector<fem::singular::HigherOrderBasis> &basis)
{
  const auto points = InteriorSamplePoints(32);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  Eigen::MatrixXd samples(6 * points.size(), basis.size());
  for (int column = 0; column < samples.cols(); column++)
  {
    for (std::size_t sample = 0; sample < points.size(); sample++)
    {
      const auto value = fem::singular::EvaluateHigherOrderBasis(
          points[sample], grad_lambda, basis[column]);
      for (int d = 0; d < 3; d++)
      {
        samples(6 * sample + d, column) = value.value[d];
        samples(6 * sample + 3 + d, column) = value.curl[d];
      }
    }

    const double norm = samples.col(column).norm();
    REQUIRE(norm > 0.0);
    samples.col(column) /= norm;
  }

  const Eigen::JacobiSVD<Eigen::MatrixXd> svd(samples);
  const auto &singular_values = svd.singularValues();
  return {singular_values.data(), singular_values.data() + singular_values.size()};
}

std::vector<fem::singular::HigherOrderBasis>
EnumerateUnreducedHigherOrderBases(fem::singular::HigherOrderBasisFamily family, int order,
                                   double nu)
{
  using fem::singular::HigherOrderBasis;
  using fem::singular::HigherOrderBasisFamily;
  const std::array<int, 4> nodes{0, 1, 2, 3};
  std::vector<HigherOrderBasis> basis;

  if (family == HigherOrderBasisFamily::NODE_GRADIENT)
  {
    const int denominator = order + 1;
    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 0; w1 <= denominator - w0; w1++)
      {
        for (int w2 = 0; w2 <= denominator - w0 - w1; w2++)
        {
          const std::array<int, 4> weights{w0, w1, w2, denominator - w0 - w1 - w2};
          for (int edge = 1; edge < 4; edge++)
          {
            if (weights[edge] == 0)
            {
              continue;
            }
            std::array<int, 4> positions{0, edge, -1, -1};
            int next = 2;
            for (int position = 1; position < 4; position++)
            {
              if (position != edge)
              {
                positions[next++] = position;
              }
            }
            basis.push_back({family,
                             {nodes[positions[0]], nodes[positions[1]], nodes[positions[2]],
                              nodes[positions[3]]},
                             {weights[positions[0]], weights[positions[1]],
                              weights[positions[2]], weights[positions[3]]},
                             order,
                             nu});
          }
        }
      }
    }
    return basis;
  }

  if (family == HigherOrderBasisFamily::NODE_ROTATIONAL)
  {
    const int denominator = order + 2;
    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 0; w1 <= denominator - w0; w1++)
      {
        for (int w2 = 0; w2 <= denominator - w0 - w1; w2++)
        {
          const std::array<int, 4> weights{w0, w1, w2, denominator - w0 - w1 - w2};
          for (int first = 1; first < 4; first++)
          {
            for (int second = first + 1; second < 4; second++)
            {
              if (weights[first] == 0 || weights[second] == 0)
              {
                continue;
              }
              int opposite = 1;
              while (opposite == first || opposite == second)
              {
                opposite++;
              }
              basis.push_back(
                  {family,
                   {nodes[0], nodes[first], nodes[second], nodes[opposite]},
                   {weights[0], weights[first], weights[second], weights[opposite]},
                   order,
                   nu});
            }
          }
        }
      }
    }
    return basis;
  }

  if (family == HigherOrderBasisFamily::EDGE_GRADIENT)
  {
    const int denominator = order + 2;
    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 1; w1 <= order; w1++)
      {
        const int remaining = denominator - w0 - w1;
        for (int w2 = 0; w2 <= remaining; w2++)
        {
          const std::array<int, 4> weights{w0, w1, w2, remaining - w2};
          for (int face = 2; face < 4; face++)
          {
            if (weights[face] == 0)
            {
              continue;
            }
            const int opposite = face == 2 ? 3 : 2;
            basis.push_back({family,
                             {nodes[0], nodes[1], nodes[face], nodes[opposite]},
                             {weights[0], weights[1], weights[face], weights[opposite]},
                             order,
                             nu});
          }
        }
      }
    }
    return basis;
  }

  REQUIRE(family == HigherOrderBasisFamily::EDGE_ROTATIONAL);
  const int denominator = order + 3;
  for (int w0 = 1; w0 <= order; w0++)
  {
    for (int w1 = 1; w1 <= order; w1++)
    {
      for (int w2 = 1; w2 <= order; w2++)
      {
        const int w3 = denominator - w0 - w1 - w2;
        if (w3 >= 1 && w3 <= order)
        {
          basis.push_back({family, nodes, {w0, w1, w2, w3}, order, nu});
        }
      }
    }
  }
  return basis;
}

void CheckScaledBasisIdentity(double scale_a, const VectorBasisValue &a, double scale_b,
                              const VectorBasisValue &b, double tolerance)
{
  for (int d = 0; d < 3; d++)
  {
    CHECK_THAT(scale_a * a.value[d] - scale_b * b.value[d], WithinAbs(0.0, tolerance));
    CHECK_THAT(scale_a * a.curl[d] - scale_b * b.curl[d], WithinAbs(0.0, tolerance));
  }
}

void CheckBasisDescriptor(const fem::singular::FirstOrderBasis &actual,
                          const fem::singular::FirstOrderBasis &expected)
{
  CHECK(actual.family == expected.family);
  CHECK(actual.indices == expected.indices);
  CHECK(actual.nu == expected.nu);
}

int CanonicalNodeFeature(const fem::singular::CanonicalFirstOrderBasisPair &canonical)
{
  for (const auto *basis : {&canonical.row_basis, &canonical.column_basis})
  {
    if (basis->family == fem::singular::FirstOrderBasisFamily::NODE_GRADIENT ||
        basis->family == fem::singular::FirstOrderBasisFamily::NODE_ROTATIONAL)
    {
      return basis->indices[0];
    }
  }
  throw std::invalid_argument("Canonical basis pair has no node-singular function!");
}

std::array<int, 2>
CanonicalEdgeFeature(const fem::singular::CanonicalFirstOrderBasisPair &canonical)
{
  for (const auto *basis : {&canonical.row_basis, &canonical.column_basis})
  {
    if (basis->family == fem::singular::FirstOrderBasisFamily::EDGE_GRADIENT ||
        basis->family == fem::singular::FirstOrderBasisFamily::EDGE_ROTATIONAL)
    {
      return {basis->indices[0], basis->indices[1]};
    }
  }
  throw std::invalid_argument("Canonical basis pair has no edge-singular function!");
}

int CanonicalNodeFeature(const fem::singular::CanonicalReferenceBasisPair &canonical)
{
  for (const auto *basis : {&canonical.row_basis, &canonical.column_basis})
  {
    const int node = std::visit(
        [](const auto &entry)
        {
          using Basis = std::decay_t<decltype(entry)>;
          if constexpr (std::is_same_v<Basis, fem::singular::FirstOrderBasis>)
          {
            return entry.family == fem::singular::FirstOrderBasisFamily::NODE_GRADIENT ||
                           entry.family ==
                               fem::singular::FirstOrderBasisFamily::NODE_ROTATIONAL
                       ? entry.indices[0]
                       : -1;
          }
          else
          {
            return entry.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT ||
                           entry.family ==
                               fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL
                       ? entry.nodes[0]
                       : -1;
          }
        },
        *basis);
    if (node >= 0)
    {
      return node;
    }
  }
  throw std::invalid_argument("Canonical basis pair has no node-singular function!");
}

std::array<int, 2>
CanonicalEdgeFeature(const fem::singular::CanonicalReferenceBasisPair &canonical)
{
  for (const auto *basis : {&canonical.row_basis, &canonical.column_basis})
  {
    const auto edge = std::visit(
        [](const auto &entry)
        {
          using Basis = std::decay_t<decltype(entry)>;
          if constexpr (std::is_same_v<Basis, fem::singular::FirstOrderBasis>)
          {
            return entry.family == fem::singular::FirstOrderBasisFamily::EDGE_GRADIENT ||
                           entry.family ==
                               fem::singular::FirstOrderBasisFamily::EDGE_ROTATIONAL
                       ? std::array{entry.indices[0], entry.indices[1]}
                       : std::array{-1, -1};
          }
          else
          {
            return entry.family == fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT ||
                           entry.family ==
                               fem::singular::HigherOrderBasisFamily::EDGE_ROTATIONAL
                       ? std::array{entry.nodes[0], entry.nodes[1]}
                       : std::array{-1, -1};
          }
        },
        *basis);
    if (edge[0] >= 0)
    {
      return edge;
    }
  }
  throw std::invalid_argument("Canonical basis pair has no edge-singular function!");
}

}  // namespace

TEST_CASE("Singular elements reference tetrahedron and edge convention",
          "[singularelements][Serial]")
{
  const Vector3 point{0.2, 0.3, 0.1};
  const auto lambda = fem::singular::ReferenceBarycentricPoint(point);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  CheckVector({lambda[0], lambda[1], lambda[2]}, {0.4, 0.2, 0.3}, 1.0e-15);
  CHECK_THAT(lambda[3], WithinAbs(0.1, 1.0e-15));

  Vector3 gradient_sum{};
  for (const auto &gradient : grad_lambda)
  {
    Add(1.0, gradient, gradient_sum);
  }
  CheckVector(gradient_sum, {0.0, 0.0, 0.0}, 0.0);

  const auto edge_12 = fem::singular::EvaluateStandardEdge(lambda, grad_lambda, 1, 2);
  CheckVector(edge_12.value, {0.3, -0.2, 0.0}, 1.0e-15);
  CheckVector(edge_12.curl, {0.0, 0.0, -2.0}, 1.0e-15);

  const auto edge_21 = fem::singular::EvaluateStandardEdge(lambda, grad_lambda, 2, 1);
  CheckVector(edge_21.value, {-0.3, 0.2, 0.0}, 1.0e-15);
  CheckVector(edge_21.curl, {0.0, 0.0, 2.0}, 1.0e-15);
}

TEST_CASE("Singular element Silvester-Lagrange polynomial kernels",
          "[singularelements][Serial]")
{
  for (int grid_denominator = 2; grid_denominator <= 7; grid_denominator++)
  {
    for (int index = 0; index <= grid_denominator; index++)
    {
      for (int grid_index = 0; grid_index <= index; grid_index++)
      {
        CAPTURE(grid_denominator, index, grid_index);
        const double coordinate = static_cast<double>(grid_index) / grid_denominator;
        const auto polynomial =
            fem::singular::EvaluateSilvesterLagrange(grid_denominator, index, coordinate);
        CHECK_THAT(polynomial.value, WithinAbs(grid_index == index ? 1.0 : 0.0, 2.0e-14));
        CHECK(std::isfinite(polynomial.derivative));
      }
    }

    for (int index = 1; index <= grid_denominator; index++)
    {
      for (int grid_index = 1; grid_index <= index; grid_index++)
      {
        CAPTURE(grid_denominator, index, grid_index);
        const double coordinate = static_cast<double>(grid_index) / grid_denominator;
        const auto polynomial = fem::singular::EvaluateShiftedSilvesterLagrange(
            grid_denominator, index, coordinate);
        CHECK_THAT(polynomial.value, WithinAbs(grid_index == index ? 1.0 : 0.0, 2.0e-14));
        CHECK(std::isfinite(polynomial.derivative));
      }
    }
  }

  constexpr double coordinate = 0.371;
  constexpr double step = 1.0e-7;
  for (int grid_denominator : {3, 5, 8})
  {
    for (int index = 0; index <= grid_denominator; index++)
    {
      CAPTURE(grid_denominator, index);
      const auto polynomial =
          fem::singular::EvaluateSilvesterLagrange(grid_denominator, index, coordinate);
      const double numerical_derivative = (fem::singular::EvaluateSilvesterLagrange(
                                               grid_denominator, index, coordinate + step)
                                               .value -
                                           fem::singular::EvaluateSilvesterLagrange(
                                               grid_denominator, index, coordinate - step)
                                               .value) /
                                          (2.0 * step);
      CHECK(std::abs(polynomial.derivative - numerical_derivative) <=
            2.0e-7 * std::max(1.0, std::abs(numerical_derivative)));
    }
    for (int index = 1; index <= grid_denominator; index++)
    {
      CAPTURE(grid_denominator, index);
      const auto polynomial = fem::singular::EvaluateShiftedSilvesterLagrange(
          grid_denominator, index, coordinate);
      const double numerical_derivative = (fem::singular::EvaluateShiftedSilvesterLagrange(
                                               grid_denominator, index, coordinate + step)
                                               .value -
                                           fem::singular::EvaluateShiftedSilvesterLagrange(
                                               grid_denominator, index, coordinate - step)
                                               .value) /
                                          (2.0 * step);
      CHECK(std::abs(polynomial.derivative - numerical_derivative) <=
            2.0e-7 * std::max(1.0, std::abs(numerical_derivative)));
    }
  }

  for (int grid_denominator : {2, 3, 4})
  {
    const auto shifted =
        fem::singular::EvaluateShiftedSilvesterLagrange(grid_denominator, 1, coordinate);
    const auto unshifted =
        fem::singular::EvaluateSilvesterLagrange(grid_denominator, 0, coordinate);
    CHECK(shifted.value == 1.0);
    CHECK(shifted.derivative == 0.0);
    CHECK(unshifted.value == 1.0);
    CHECK(unshifted.derivative == 0.0);
  }

  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(0, 0, coordinate),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(2, -1, coordinate),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(2, 3, coordinate),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateShiftedSilvesterLagrange(2, 0, coordinate),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateShiftedSilvesterLagrange(2, 3, coordinate),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(
                      2, 1, std::numeric_limits<double>::infinity()),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(2, 1, -0.1),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateSilvesterLagrange(2, 1, 1.1),
                  std::invalid_argument);
}

TEST_CASE("Singular element analytic gradients and curls", "[singularelements][Serial]")
{
  const Vector3 point{0.2, 0.25, 0.15};
  const auto lambda = fem::singular::ReferenceBarycentricPoint(point);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  for (double nu : {0.5, 2.0 / 3.0, 0.731})
  {
    CAPTURE(nu);

    const auto node_gradient =
        fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu);
    const auto node_potential = [=](const Vector3 &x)
    {
      return fem::singular::EvaluateNodeGradientPotential(
          fem::singular::ReferenceBarycentricPoint(x), 0, 1, nu);
    };
    const auto node_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateNodeGradient(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, nu)
          .value;
    };
    CheckVector(node_gradient.value, NumericalGradient(point, node_potential), 2.0e-9);
    CheckVector(node_gradient.curl, {0.0, 0.0, 0.0}, 0.0);
    CheckVector(NumericalCurl(point, node_vector), {0.0, 0.0, 0.0}, 2.0e-9);

    const auto edge_gradient =
        fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu);
    const auto edge_potential = [=](const Vector3 &x)
    {
      return fem::singular::EvaluateEdgeGradientPotential(
          fem::singular::ReferenceBarycentricPoint(x), 0, 1, 2, nu);
    };
    const auto edge_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateEdgeGradient(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, nu)
          .value;
    };
    CheckVector(edge_gradient.value, NumericalGradient(point, edge_potential), 2.0e-9);
    CheckVector(edge_gradient.curl, {0.0, 0.0, 0.0}, 0.0);
    CheckVector(NumericalCurl(point, edge_vector), {0.0, 0.0, 0.0}, 2.0e-9);

    const auto node_rotational =
        fem::singular::EvaluateNodeRotational(lambda, grad_lambda, 0, 1, 2, nu);
    const auto node_rotational_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateNodeRotational(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, nu)
          .value;
    };
    CheckVector(node_rotational.curl, NumericalCurl(point, node_rotational_vector), 2.0e-9);

    const auto edge_rotational =
        fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3, nu);
    const auto edge_rotational_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateEdgeRotational(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, 3, nu)
          .value;
    };
    CheckVector(edge_rotational.curl, NumericalCurl(point, edge_rotational_vector), 2.0e-9);
  }
}

TEST_CASE("Singular element higher-order node-gradient basis from equation 28",
          "[singularelements][Serial]")
{
  constexpr double nu = 2.0 / 3.0;
  const Vector3 point{0.2, 0.25, 0.15};
  const auto lambda = fem::singular::ReferenceBarycentricPoint(point);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  const fem::singular::InterpolationIndices first_order_indices{1, 1, 0, 0};
  CHECK_THAT(
      fem::singular::EvaluateHigherOrderNodeGradientPotential(lambda, 0, 1, 2, 3,
                                                              first_order_indices, 1, nu),
      WithinAbs(fem::singular::EvaluateNodeGradientPotential(lambda, 0, 1, nu), 0.0));
  CheckBasisValue(fem::singular::EvaluateHigherOrderNodeGradient(
                      lambda, grad_lambda, 0, 1, 2, 3, first_order_indices, 1, nu),
                  fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu), 0.0);

  const std::array higher_order_cases{
      std::pair{2, fem::singular::InterpolationIndices{1, 1, 1, 0}},
      std::pair{2, fem::singular::InterpolationIndices{2, 1, 0, 0}},
      std::pair{3, fem::singular::InterpolationIndices{1, 1, 1, 1}},
      std::pair{4, fem::singular::InterpolationIndices{2, 1, 1, 1}}};
  for (const auto &higher_order_case : higher_order_cases)
  {
    const int order = higher_order_case.first;
    const auto indices = higher_order_case.second;
    CAPTURE(order, indices);
    const auto basis = fem::singular::EvaluateHigherOrderNodeGradient(
        lambda, grad_lambda, 0, 1, 2, 3, indices, order, nu);
    const auto potential = [=](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderNodeGradientPotential(
          fem::singular::ReferenceBarycentricPoint(x), 0, 1, 2, 3, indices, order, nu);
    };
    const auto vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderNodeGradient(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, 3,
                 indices, order, nu)
          .value;
    };
    CheckVector(basis.value, NumericalGradient(point, potential), 2.0e-8);
    CheckVector(basis.curl, {0.0, 0.0, 0.0}, 0.0);
    CheckVector(NumericalCurl(point, vector), {0.0, 0.0, 0.0}, 2.0e-8);
  }

  const fem::singular::InterpolationIndices trace_indices{1, 1, 1, 1};
  for (int opposite : {0, 1})
  {
    const auto basis = fem::singular::EvaluateHigherOrderNodeGradient(
        FacePoint(opposite), grad_lambda, 0, 1, 2, 3, trace_indices, 3, nu);
    CHECK_THAT(Norm(TangentialPart(basis.value, grad_lambda[opposite])),
               WithinAbs(0.0, 2.0e-14));
  }

  const fem::singular::InterpolationIndices permutation_indices{2, 1, 1, 0};
  const auto reference = fem::singular::EvaluateHigherOrderNodeGradient(
      lambda, grad_lambda, 0, 1, 2, 3, permutation_indices, 3, nu);
  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    const auto permuted_lambda =
        fem::singular::ApplyBarycentricPermutation(lambda, permutation);
    const auto permuted_gradients =
        fem::singular::ApplyBarycentricPermutation(grad_lambda, permutation);
    CheckBasisValue(fem::singular::EvaluateHigherOrderNodeGradient(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], permutation[3], permutation_indices, 3, nu),
                    reference, 2.0e-14);
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderNodeGradient(lambda, grad_lambda, 0, 1,
                                                                 2, 3, {1, 1, 0, 0}, 0, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderNodeGradient(lambda, grad_lambda, 0, 1,
                                                                 2, 3, {0, 2, 0, 0}, 1, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderNodeGradient(lambda, grad_lambda, 0, 1,
                                                                 2, 3, {1, 1, 1, 0}, 1, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderNodeGradient(lambda, grad_lambda, 0, 1,
                                                                 2, 2, {1, 1, 0, 0}, 1, nu),
                  std::invalid_argument);
}

TEST_CASE("Singular element higher-order bases from equations 30 32 and 33",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.61;
  const Vector3 point{0.2, 0.25, 0.15};
  const auto lambda = fem::singular::ReferenceBarycentricPoint(point);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const fem::singular::InterpolationIndices first_order_face{1, 1, 1, 0};
  const fem::singular::InterpolationIndices first_order_volume{1, 1, 1, 1};

  CheckBasisValue(fem::singular::EvaluateHigherOrderNodeRotational(
                      lambda, grad_lambda, 0, 1, 2, 3, first_order_face, 1, nu),
                  fem::singular::EvaluateNodeRotational(lambda, grad_lambda, 0, 1, 2, nu),
                  0.0);
  CHECK_THAT(
      fem::singular::EvaluateHigherOrderEdgeGradientPotential(lambda, 0, 1, 2, 3,
                                                              first_order_face, 1, nu),
      WithinAbs(fem::singular::EvaluateEdgeGradientPotential(lambda, 0, 1, 2, nu), 0.0));
  CheckBasisValue(fem::singular::EvaluateHigherOrderEdgeGradient(
                      lambda, grad_lambda, 0, 1, 2, 3, first_order_face, 1, nu),
                  fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu),
                  0.0);
  CheckBasisValue(
      fem::singular::EvaluateHigherOrderEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3,
                                                       first_order_volume, 1, nu),
      fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3, nu), 0.0);

  const std::array face_cases{
      std::pair{2, fem::singular::InterpolationIndices{1, 1, 1, 1}},
      std::pair{2, fem::singular::InterpolationIndices{2, 1, 1, 0}},
      std::pair{3, fem::singular::InterpolationIndices{2, 1, 1, 1}}};
  for (const auto &face_case : face_cases)
  {
    const int order = face_case.first;
    const auto indices = face_case.second;
    CAPTURE(order, indices);
    const auto node_rotational = fem::singular::EvaluateHigherOrderNodeRotational(
        lambda, grad_lambda, 0, 1, 2, 3, indices, order, nu);
    const auto node_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderNodeRotational(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, 3,
                 indices, order, nu)
          .value;
    };
    CheckVector(node_rotational.curl, NumericalCurl(point, node_vector), 3.0e-8);

    const auto edge_gradient = fem::singular::EvaluateHigherOrderEdgeGradient(
        lambda, grad_lambda, 0, 1, 2, 3, indices, order, nu);
    const auto edge_potential = [=](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderEdgeGradientPotential(
          fem::singular::ReferenceBarycentricPoint(x), 0, 1, 2, 3, indices, order, nu);
    };
    const auto edge_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderEdgeGradient(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, 3,
                 indices, order, nu)
          .value;
    };
    CheckVector(edge_gradient.value, NumericalGradient(point, edge_potential), 3.0e-8);
    CheckVector(edge_gradient.curl, {0.0, 0.0, 0.0}, 0.0);
    CheckVector(NumericalCurl(point, edge_vector), {0.0, 0.0, 0.0}, 3.0e-8);
  }

  const std::array volume_cases{
      std::pair{2, fem::singular::InterpolationIndices{2, 1, 1, 1}},
      std::pair{3, fem::singular::InterpolationIndices{2, 2, 1, 1}}};
  for (const auto &volume_case : volume_cases)
  {
    const int order = volume_case.first;
    const auto indices = volume_case.second;
    CAPTURE(order, indices);
    const auto edge_rotational = fem::singular::EvaluateHigherOrderEdgeRotational(
        lambda, grad_lambda, 0, 1, 2, 3, indices, order, nu);
    const auto edge_vector = [=, &grad_lambda](const Vector3 &x)
    {
      return fem::singular::EvaluateHigherOrderEdgeRotational(
                 fem::singular::ReferenceBarycentricPoint(x), grad_lambda, 0, 1, 2, 3,
                 indices, order, nu)
          .value;
    };
    CheckVector(edge_rotational.curl, NumericalCurl(point, edge_vector), 3.0e-8);
  }

  const fem::singular::InterpolationIndices face_trace_indices{1, 1, 1, 1};
  for (int opposite : {0, 1, 2})
  {
    const auto node_rotational = fem::singular::EvaluateHigherOrderNodeRotational(
        FacePoint(opposite), grad_lambda, 0, 1, 2, 3, face_trace_indices, 2, nu);
    const auto edge_gradient = fem::singular::EvaluateHigherOrderEdgeGradient(
        FacePoint(opposite), grad_lambda, 0, 1, 2, 3, face_trace_indices, 2, nu);
    CHECK_THAT(Norm(TangentialPart(node_rotational.value, grad_lambda[opposite])),
               WithinAbs(0.0, 2.0e-14));
    CHECK_THAT(Norm(TangentialPart(edge_gradient.value, grad_lambda[opposite])),
               WithinAbs(0.0, 2.0e-14));
  }
  const fem::singular::InterpolationIndices volume_trace_indices{2, 1, 1, 1};
  for (int opposite = 0; opposite < 4; opposite++)
  {
    const auto edge_rotational = fem::singular::EvaluateHigherOrderEdgeRotational(
        FacePoint(opposite), grad_lambda, 0, 1, 2, 3, volume_trace_indices, 2, nu);
    CHECK_THAT(Norm(TangentialPart(edge_rotational.value, grad_lambda[opposite])),
               WithinAbs(0.0, 2.0e-14));
  }

  const fem::singular::InterpolationIndices permutation_face_indices{2, 1, 1, 0};
  const fem::singular::InterpolationIndices permutation_volume_indices{2, 1, 1, 1};
  const auto reference_node = fem::singular::EvaluateHigherOrderNodeRotational(
      lambda, grad_lambda, 0, 1, 2, 3, permutation_face_indices, 2, nu);
  const auto reference_gradient = fem::singular::EvaluateHigherOrderEdgeGradient(
      lambda, grad_lambda, 0, 1, 2, 3, permutation_face_indices, 2, nu);
  const auto reference_rotational = fem::singular::EvaluateHigherOrderEdgeRotational(
      lambda, grad_lambda, 0, 1, 2, 3, permutation_volume_indices, 2, nu);
  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    const auto permuted_lambda =
        fem::singular::ApplyBarycentricPermutation(lambda, permutation);
    const auto permuted_gradients =
        fem::singular::ApplyBarycentricPermutation(grad_lambda, permutation);
    CheckBasisValue(fem::singular::EvaluateHigherOrderNodeRotational(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], permutation[3], permutation_face_indices, 2, nu),
                    reference_node, 2.0e-14);
    CheckBasisValue(fem::singular::EvaluateHigherOrderEdgeGradient(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], permutation[3], permutation_face_indices, 2, nu),
                    reference_gradient, 2.0e-14);
    CheckBasisValue(fem::singular::EvaluateHigherOrderEdgeRotational(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], permutation[3], permutation_volume_indices, 2, nu),
                    reference_rotational, 2.0e-14);
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderNodeRotational(
                      lambda, grad_lambda, 0, 1, 2, 3, {1, 1, 1, 0}, 0, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderEdgeGradient(lambda, grad_lambda, 0, 1,
                                                                 2, 3, {1, 1, 1, 2}, 2, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderEdgeRotational(
                      lambda, grad_lambda, 0, 1, 2, 3, {1, 1, 1, 0}, 1, nu),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderEdgeRotational(
                      lambda, grad_lambda, 0, 1, 2, 3, {2, 2, 2, 2}, 2, nu),
                  std::invalid_argument);
}

TEST_CASE("Singular element generic H1 potential evaluation", "[singularelements][Serial]")
{
  constexpr int order = 2;
  constexpr double nu = 2.0 / 3.0;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const fem::singular::BarycentricPoint lambda{0.31, 0.23, 0.19, 0.27};

  const auto node_gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, order, nu);
  const auto edge_gradients =
      fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, order, nu);
  const auto node_rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, order, nu);
  const auto edge_rotations =
      fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, order, nu);
  REQUIRE_FALSE(node_gradients.empty());
  REQUIRE_FALSE(edge_gradients.empty());
  REQUIRE_FALSE(node_rotations.empty());
  REQUIRE_FALSE(edge_rotations.empty());

  const auto &node = node_gradients.front();
  CHECK_THAT(fem::singular::EvaluateHigherOrderGradientPotential(lambda, node),
             WithinAbs(fem::singular::EvaluateHigherOrderNodeGradientPotential(
                           lambda, node.nodes[0], node.nodes[1], node.nodes[2],
                           node.nodes[3], node.interpolation_indices, node.order, node.nu),
                       0.0));
  const auto &edge = edge_gradients.front();
  CHECK_THAT(fem::singular::EvaluateHigherOrderGradientPotential(lambda, edge),
             WithinAbs(fem::singular::EvaluateHigherOrderEdgeGradientPotential(
                           lambda, edge.nodes[0], edge.nodes[1], edge.nodes[2],
                           edge.nodes[3], edge.interpolation_indices, edge.order, edge.nu),
                       0.0));
  CHECK_THROWS_AS(
      fem::singular::EvaluateHigherOrderGradientPotential(lambda, node_rotations.front()),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::EvaluateHigherOrderGradientPotential(lambda, edge_rotations.front()),
      std::invalid_argument);
}

TEST_CASE("Singular element higher-order retained sets match Table I",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const auto check_unique = [](const std::vector<fem::singular::HigherOrderBasis> &basis)
  {
    std::set<std::tuple<int, std::array<int, 4>, fem::singular::InterpolationIndices, int>>
        keys;
    for (const auto &entry : basis)
    {
      CHECK(keys.emplace(static_cast<int>(entry.family), entry.nodes,
                         entry.interpolation_indices, entry.order)
                .second);
    }
  };

  for (int order = 1; order <= 4; order++)
  {
    CAPTURE(order);
    const auto node_gradient =
        fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, order, nu);
    const auto node_rotational =
        fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, order, nu);
    const auto edge_gradient =
        fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, order, nu);
    const auto edge_rotational =
        fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, order, nu);

    const std::size_t expected_node_gradient =
        3 * order + 3 * order * (order - 1) / 2 + order * (order - 1) * (order - 2) / 6;
    const std::size_t expected_node_rotational =
        3 * (order + 1) * order / 2 + (order + 1) * order * (order - 1) / 3;
    const std::size_t expected_edge_gradient =
        (order + 1) * order + (order + 1) * order * (order - 1) / 6;
    const std::size_t expected_edge_rotational = (order + 2) * (order + 1) * order / 6;
    CHECK(node_gradient.size() == expected_node_gradient);
    CHECK(node_rotational.size() == expected_node_rotational);
    CHECK(edge_gradient.size() == expected_edge_gradient);
    CHECK(edge_rotational.size() == expected_edge_rotational);
    check_unique(node_gradient);
    check_unique(node_rotational);
    check_unique(edge_gradient);
    check_unique(edge_rotational);
  }

  const BarycentricPoint lambda{0.37, 0.21, 0.18, 0.24};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  for (const auto &basis :
       fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu))
  {
    const auto evaluated =
        fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
    CHECK(fem::singular::EvaluateHigherOrderBasisValue(lambda, grad_lambda, basis) ==
          evaluated.value);
    CheckBasisValue(evaluated,
                    fem::singular::EvaluateNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                                        basis.nodes[1], nu),
                    0.0);
  }
  for (const auto &basis :
       fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu))
  {
    const auto evaluated =
        fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
    CHECK(fem::singular::EvaluateHigherOrderBasisValue(lambda, grad_lambda, basis) ==
          evaluated.value);
    CheckBasisValue(evaluated,
                    fem::singular::EvaluateNodeRotational(lambda, grad_lambda,
                                                          basis.nodes[0], basis.nodes[1],
                                                          basis.nodes[2], nu),
                    0.0);
  }
  for (const auto &basis :
       fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 1, nu))
  {
    const auto evaluated =
        fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
    CHECK(fem::singular::EvaluateHigherOrderBasisValue(lambda, grad_lambda, basis) ==
          evaluated.value);
    CheckBasisValue(evaluated,
                    fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, basis.nodes[0],
                                                        basis.nodes[1], basis.nodes[2], nu),
                    0.0);
  }
  for (const auto &basis :
       fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, 1, nu))
  {
    const auto evaluated =
        fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
    CHECK(fem::singular::EvaluateHigherOrderBasisValue(lambda, grad_lambda, basis) ==
          evaluated.value);
    CheckBasisValue(evaluated,
                    fem::singular::EvaluateEdgeRotational(
                        lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
                        basis.nodes[3], nu),
                    0.0);
  }

  const auto node_rotational =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu)
          .front();
  BarycentricPoint singular_node{0.0, 0.0, 0.0, 0.0};
  singular_node[node_rotational.nodes[0]] = 1.0;
  CHECK(fem::singular::EvaluateHigherOrderBasisValue(
            singular_node, grad_lambda, node_rotational) == Vector3{0.0, 0.0, 0.0});
  CHECK_THROWS_AS(
      fem::singular::EvaluateHigherOrderBasis(singular_node, grad_lambda, node_rotational),
      std::domain_error);

  const auto edge_rotational =
      fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, 1, nu)
          .front();
  BarycentricPoint singular_edge{0.0, 0.0, 0.0, 0.0};
  singular_edge[edge_rotational.nodes[0]] = 0.5;
  singular_edge[edge_rotational.nodes[1]] = 0.5;
  CHECK(fem::singular::EvaluateHigherOrderBasisValue(
            singular_edge, grad_lambda, edge_rotational) == Vector3{0.0, 0.0, 0.0});
  CHECK_THROWS_AS(
      fem::singular::EvaluateHigherOrderBasis(singular_edge, grad_lambda, edge_rotational),
      std::domain_error);

  const fem::singular::BarycentricPermutation permutation{2, 0, 3, 1};
  std::array<int, 4> permuted_nodes;
  for (int i = 0; i < 4; i++)
  {
    permuted_nodes[i] = permutation[canonical_nodes[i]];
  }
  const auto check_permuted =
      [&](const std::vector<fem::singular::HigherOrderBasis> &reference,
          const std::vector<fem::singular::HigherOrderBasis> &permuted)
  {
    REQUIRE(permuted.size() == reference.size());
    for (std::size_t i = 0; i < reference.size(); i++)
    {
      const auto expected =
          fem::singular::ApplyBarycentricPermutation(reference[i], permutation);
      CHECK(permuted[i].family == expected.family);
      CHECK(permuted[i].nodes == expected.nodes);
      CHECK(permuted[i].interpolation_indices == expected.interpolation_indices);
      CHECK(permuted[i].order == expected.order);
      CHECK(permuted[i].nu == expected.nu);
    }
  };
  check_permuted(
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 4, nu),
      fem::singular::EnumerateHigherOrderNodeGradientBases(permuted_nodes, 4, nu));
  check_permuted(
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 4, nu),
      fem::singular::EnumerateHigherOrderNodeRotationalBases(permuted_nodes, 4, nu));
  check_permuted(
      fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 4, nu),
      fem::singular::EnumerateHigherOrderEdgeGradientBases(permuted_nodes, 4, nu));
  check_permuted(
      fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, 4, nu),
      fem::singular::EnumerateHigherOrderEdgeRotationalBases(permuted_nodes, 4, nu));

  CHECK_THROWS_AS(
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 0, nu),
      std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EnumerateHigherOrderEdgeGradientBases({0, 1, 2, 2}, 1, nu),
                  std::invalid_argument);
  const fem::singular::HigherOrderBasis unknown{
      static_cast<fem::singular::HigherOrderBasisFamily>(99),
      canonical_nodes,
      {1, 1, 0, 0},
      1,
      nu};
  CHECK_THROWS_AS(fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, unknown),
                  std::invalid_argument);
}

TEST_CASE("Singular element higher-order retained sets have the required rank",
          "[singularelements][Serial]")
{
  using fem::singular::HigherOrderBasisFamily;
  constexpr double nu = 0.5;
  constexpr std::array<int, 4> nodes{0, 1, 2, 3};
  const std::array families{
      HigherOrderBasisFamily::NODE_GRADIENT, HigherOrderBasisFamily::NODE_ROTATIONAL,
      HigherOrderBasisFamily::EDGE_GRADIENT, HigherOrderBasisFamily::EDGE_ROTATIONAL};

  for (int order = 1; order <= 4; order++)
  {
    for (const auto family : families)
    {
      CAPTURE(order, family);
      std::vector<fem::singular::HigherOrderBasis> retained;
      switch (family)
      {
        case HigherOrderBasisFamily::NODE_GRADIENT:
          retained = fem::singular::EnumerateHigherOrderNodeGradientBases(nodes, order, nu);
          break;
        case HigherOrderBasisFamily::NODE_ROTATIONAL:
          retained =
              fem::singular::EnumerateHigherOrderNodeRotationalBases(nodes, order, nu);
          break;
        case HigherOrderBasisFamily::EDGE_GRADIENT:
          retained = fem::singular::EnumerateHigherOrderEdgeGradientBases(nodes, order, nu);
          break;
        case HigherOrderBasisFamily::EDGE_ROTATIONAL:
          retained =
              fem::singular::EnumerateHigherOrderEdgeRotationalBases(nodes, order, nu);
          break;
      }
      const auto unreduced = EnumerateUnreducedHigherOrderBases(family, order, nu);
      const auto retained_singular_values = SampledBasisSingularValues(retained);
      const auto unreduced_singular_values = SampledBasisSingularValues(unreduced);
      REQUIRE(retained_singular_values.size() == retained.size());
      REQUIRE(unreduced_singular_values.size() == unreduced.size());

      const double retained_threshold = 1.0e-10 * retained_singular_values.front();
      const double unreduced_threshold = 1.0e-10 * unreduced_singular_values.front();
      CHECK(std::count_if(retained_singular_values.begin(), retained_singular_values.end(),
                          [=](double value) { return value > retained_threshold; }) ==
            static_cast<std::ptrdiff_t>(retained.size()));
      CHECK(std::count_if(unreduced_singular_values.begin(),
                          unreduced_singular_values.end(),
                          [=](double value) { return value > unreduced_threshold; }) ==
            static_cast<std::ptrdiff_t>(retained.size()));
    }
  }
}

TEST_CASE("Singular element discarded higher-order functions satisfy dependencies",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  const auto points = InteriorSamplePoints(8);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  for (int order = 2; order <= 4; order++)
  {
    const int node_denominator = order + 1;
    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 0; w1 <= node_denominator - w0; w1++)
      {
        for (int w2 = 0; w2 <= node_denominator - w0 - w1; w2++)
        {
          const std::array<int, 4> weights{w0, w1, w2, node_denominator - w0 - w1 - w2};
          for (int first = 1; first < 4; first++)
          {
            for (int second = first + 1; second < 4; second++)
            {
              if (weights[first] == 0 || weights[second] == 0)
              {
                continue;
              }
              std::array<int, 2> others;
              int next = 0;
              for (int position = 1; position < 4; position++)
              {
                if (position != first)
                {
                  others[next++] = position;
                }
              }
              const fem::singular::HigherOrderBasis first_basis{
                  fem::singular::HigherOrderBasisFamily::NODE_GRADIENT,
                  {0, first, others[0], others[1]},
                  {weights[0], weights[first], weights[others[0]], weights[others[1]]},
                  order,
                  nu};
              next = 0;
              for (int position = 1; position < 4; position++)
              {
                if (position != second)
                {
                  others[next++] = position;
                }
              }
              const fem::singular::HigherOrderBasis second_basis{
                  fem::singular::HigherOrderBasisFamily::NODE_GRADIENT,
                  {0, second, others[0], others[1]},
                  {weights[0], weights[second], weights[others[0]], weights[others[1]]},
                  order,
                  nu};
              for (const auto &lambda : points)
              {
                CAPTURE(order, weights, first, second, lambda);
                CheckScaledBasisIdentity(weights[second],
                                         fem::singular::EvaluateHigherOrderBasis(
                                             lambda, grad_lambda, first_basis),
                                         weights[first],
                                         fem::singular::EvaluateHigherOrderBasis(
                                             lambda, grad_lambda, second_basis),
                                         2.0e-12);
              }
            }
          }
        }
      }
    }

    const int face_denominator = order + 2;
    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 1; w1 <= order; w1++)
      {
        for (int w2 = 1; w2 <= order; w2++)
        {
          const int w3 = face_denominator - w0 - w1 - w2;
          if (w3 < 1 || w3 > order - 1)
          {
            continue;
          }
          const std::array<int, 4> weights{w0, w1, w2, w3};
          const fem::singular::HigherOrderBasis face_12{
              fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
              {0, 1, 2, 3},
              weights,
              order,
              nu};
          const fem::singular::HigherOrderBasis face_13{
              fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
              {0, 1, 3, 2},
              {w0, w1, w3, w2},
              order,
              nu};
          const fem::singular::HigherOrderBasis face_23{
              fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
              {0, 2, 3, 1},
              {w0, w2, w3, w1},
              order,
              nu};
          for (const auto &lambda : points)
          {
            CAPTURE(order, weights, lambda);
            const auto value_12 =
                fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, face_12);
            const auto value_13 =
                fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, face_13);
            const auto value_23 =
                fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, face_23);
            for (int d = 0; d < 3; d++)
            {
              CHECK_THAT(w3 * value_12.value[d] - w2 * value_13.value[d] +
                             w1 * value_23.value[d],
                         WithinAbs(0.0, 2.0e-12));
              CHECK_THAT(w3 * value_12.curl[d] - w2 * value_13.curl[d] +
                             w1 * value_23.curl[d],
                         WithinAbs(0.0, 2.0e-12));
            }
          }
        }
      }
    }

    for (int w0 = 1; w0 <= order; w0++)
    {
      for (int w1 = 1; w1 <= order; w1++)
      {
        for (int w2 = 1; w2 <= order; w2++)
        {
          const int w3 = face_denominator - w0 - w1 - w2;
          if (w3 < 1 || w3 > order - 1)
          {
            continue;
          }
          const fem::singular::HigherOrderBasis face_2{
              fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT,
              {0, 1, 2, 3},
              {w0, w1, w2, w3},
              order,
              nu};
          const fem::singular::HigherOrderBasis face_3{
              fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT,
              {0, 1, 3, 2},
              {w0, w1, w3, w2},
              order,
              nu};
          for (const auto &lambda : points)
          {
            CAPTURE(order, w0, w1, w2, w3, lambda);
            CheckScaledBasisIdentity(
                w3, fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, face_2),
                w2, fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, face_3),
                2.0e-12);
          }
        }
      }
    }
  }
}

TEST_CASE("Singular element tangential traces", "[singularelements][Serial]")
{
  constexpr double nu = 2.0 / 3.0;
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  for (int opposite : {0, 1})
  {
    const auto basis =
        fem::singular::EvaluateNodeGradient(FacePoint(opposite), grad_lambda, 0, 1, nu);
    CHECK_THAT(Norm(TangentialPart(basis.value, grad_lambda[opposite])),
               WithinAbs(0.0, 1.0e-14));
  }

  for (int opposite : {0, 1, 2})
  {
    const auto node_basis = fem::singular::EvaluateNodeRotational(FacePoint(opposite),
                                                                  grad_lambda, 0, 1, 2, nu);
    CHECK_THAT(Norm(TangentialPart(node_basis.value, grad_lambda[opposite])),
               WithinAbs(0.0, 1.0e-14));

    const auto edge_basis =
        fem::singular::EvaluateEdgeGradient(FacePoint(opposite), grad_lambda, 0, 1, 2, nu);
    CHECK_THAT(Norm(TangentialPart(edge_basis.value, grad_lambda[opposite])),
               WithinAbs(0.0, 1.0e-14));
  }

  for (int opposite = 0; opposite < 4; opposite++)
  {
    const auto basis = fem::singular::EvaluateEdgeRotational(FacePoint(opposite),
                                                             grad_lambda, 0, 1, 2, 3, nu);
    CHECK_THAT(Norm(TangentialPart(basis.value, grad_lambda[opposite])),
               WithinAbs(0.0, 1.0e-14));
  }
}

TEST_CASE("Singular element asymptotic exponents", "[singularelements][Serial]")
{
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  constexpr double rho_1 = 1.0e-10;
  constexpr double rho_2 = 1.0e-8;

  for (double nu : {0.5, 2.0 / 3.0, 0.731})
  {
    CAPTURE(nu);

    BarycentricPoint node_1{1.0 - rho_1, rho_1, 0.0, 0.0};
    BarycentricPoint node_2{1.0 - rho_2, rho_2, 0.0, 0.0};
    const double node_norm_1 =
        Norm(fem::singular::EvaluateNodeGradient(node_1, grad_lambda, 0, 1, nu).value);
    const double node_norm_2 =
        Norm(fem::singular::EvaluateNodeGradient(node_2, grad_lambda, 0, 1, nu).value);
    CHECK_THAT(LogSlope(rho_1, node_norm_1, rho_2, node_norm_2),
               WithinAbs(nu - 1.0, 2.0e-3));

    BarycentricPoint edge_1{0.5 * (1.0 - rho_1), 0.5 * (1.0 - rho_1), rho_1, 0.0};
    BarycentricPoint edge_2{0.5 * (1.0 - rho_2), 0.5 * (1.0 - rho_2), rho_2, 0.0};
    const double edge_norm_1 =
        Norm(fem::singular::EvaluateEdgeGradient(edge_1, grad_lambda, 0, 1, 2, nu).value);
    const double edge_norm_2 =
        Norm(fem::singular::EvaluateEdgeGradient(edge_2, grad_lambda, 0, 1, 2, nu).value);
    CHECK_THAT(LogSlope(rho_1, edge_norm_1, rho_2, edge_norm_2),
               WithinAbs(nu - 1.0, 2.0e-3));

    auto corrected_node_curl = [&](double rho)
    {
      const BarycentricPoint lambda{1.0 - rho, 0.2 * rho, 0.3 * rho, 0.5 * rho};
      auto curl =
          fem::singular::EvaluateNodeRotational(lambda, grad_lambda, 0, 1, 2, nu).curl;
      Add(1.0, fem::singular::EvaluateStandardEdge(lambda, grad_lambda, 1, 2).curl, curl);
      return Norm(curl);
    };
    CHECK_THAT(
        LogSlope(rho_1, corrected_node_curl(rho_1), rho_2, corrected_node_curl(rho_2)),
        WithinAbs(nu, 2.0e-8));

    auto corrected_edge_curl = [&](double rho)
    {
      const BarycentricPoint lambda{0.5 * (1.0 - rho), 0.5 * (1.0 - rho), 0.3 * rho,
                                    0.7 * rho};
      auto curl =
          fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3, nu).curl;
      Add(1.0, EvaluatePolynomialWeightedEdge(lambda, grad_lambda, 0, 1, 2, 3).curl, curl);
      return Norm(curl);
    };
    CHECK_THAT(
        LogSlope(rho_1, corrected_edge_curl(rho_1), rho_2, corrected_edge_curl(rho_2)),
        WithinAbs(nu, 2.0e-3));
  }
}

TEST_CASE("Singular element tetrahedron permutation covariance",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.61;
  const BarycentricPoint lambda{0.37, 0.21, 0.18, 0.24};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const auto standard = fem::singular::EvaluateStandardEdge(lambda, grad_lambda, 0, 1);
  const auto node_gradient =
      fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu);
  const auto node_rotational =
      fem::singular::EvaluateNodeRotational(lambda, grad_lambda, 0, 1, 2, nu);
  const auto edge_gradient =
      fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu);
  const auto edge_rotational =
      fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3, nu);
  const double node_potential =
      fem::singular::EvaluateNodeGradientPotential(lambda, 0, 1, nu);
  const double edge_potential =
      fem::singular::EvaluateEdgeGradientPotential(lambda, 0, 1, 2, nu);

  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    const auto permuted_lambda =
        fem::singular::ApplyBarycentricPermutation(lambda, permutation);
    const auto permuted_gradients =
        fem::singular::ApplyBarycentricPermutation(grad_lambda, permutation);

    CheckBasisValue(fem::singular::EvaluateStandardEdge(permuted_lambda, permuted_gradients,
                                                        permutation[0], permutation[1]),
                    standard, 1.0e-14);
    CheckBasisValue(fem::singular::EvaluateNodeGradient(permuted_lambda, permuted_gradients,
                                                        permutation[0], permutation[1], nu),
                    node_gradient, 1.0e-14);
    CheckBasisValue(fem::singular::EvaluateNodeRotational(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], nu),
                    node_rotational, 1.0e-14);
    CheckBasisValue(fem::singular::EvaluateEdgeGradient(permuted_lambda, permuted_gradients,
                                                        permutation[0], permutation[1],
                                                        permutation[2], nu),
                    edge_gradient, 1.0e-14);
    CheckBasisValue(fem::singular::EvaluateEdgeRotational(
                        permuted_lambda, permuted_gradients, permutation[0], permutation[1],
                        permutation[2], permutation[3], nu),
                    edge_rotational, 1.0e-14);
    CHECK_THAT(fem::singular::EvaluateNodeGradientPotential(permuted_lambda, permutation[0],
                                                            permutation[1], nu),
               WithinAbs(node_potential, 1.0e-14));
    CHECK_THAT(fem::singular::EvaluateEdgeGradientPotential(
                   permuted_lambda, permutation[0], permutation[1], permutation[2], nu),
               WithinAbs(edge_potential, 1.0e-14));
  } while (std::next_permutation(permutation.begin(), permutation.end()));
}

TEST_CASE("Singular element reference-table canonicalization", "[singularelements][Serial]")
{
  constexpr double nu = 2.0 / 3.0;
  const auto row = fem::singular::MakeNodeRotational(2, 0, 3, nu);
  const auto column = fem::singular::MakeStandardNedelec(1, 2);
  const auto canonical = fem::singular::CanonicalizeFirstOrderBasisPair(row, column);
  const auto swapped = fem::singular::CanonicalizeFirstOrderBasisPair(column, row);
  CHECK(swapped.input_to_canonical == canonical.input_to_canonical);
  CheckBasisDescriptor(swapped.row_basis, canonical.column_basis);
  CheckBasisDescriptor(swapped.column_basis, canonical.row_basis);

  constexpr int order = 8;
  constexpr int subdivisions = 2;
  const auto reference =
      fem::singular::ComputeFirstOrderReferenceIntegral(row, column, order, subdivisions);
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const double reference_mass =
      fem::singular::ContractFirstOrderMass(reference, grad_lambda, 1.0);
  const double reference_curl =
      fem::singular::ContractFirstOrderCurlCurl(reference, grad_lambda, 1.0);

  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    CAPTURE(permutation);
    const auto permuted_row = fem::singular::ApplyBarycentricPermutation(row, permutation);
    const auto permuted_column =
        fem::singular::ApplyBarycentricPermutation(column, permutation);
    const auto permuted_canonical =
        fem::singular::CanonicalizeFirstOrderBasisPair(permuted_row, permuted_column);
    CheckBasisDescriptor(permuted_canonical.row_basis, canonical.row_basis);
    CheckBasisDescriptor(permuted_canonical.column_basis, canonical.column_basis);

    const auto integral = fem::singular::ComputeFirstOrderReferenceIntegral(
        permuted_row, permuted_column, order, subdivisions);
    CHECK(integral.mass == reference.mass);
    CHECK(integral.curl_curl == reference.curl_curl);
    const auto permuted_grad_lambda =
        fem::singular::ApplyBarycentricPermutation(grad_lambda, permutation);
    CHECK_THAT(fem::singular::ContractFirstOrderMass(integral, permuted_grad_lambda, 1.0),
               WithinAbs(reference_mass, 2.0e-14));
    CHECK_THAT(
        fem::singular::ContractFirstOrderCurlCurl(integral, permuted_grad_lambda, 1.0),
        WithinAbs(reference_curl, 2.0e-14));
  } while (std::next_permutation(permutation.begin(), permutation.end()));
}

TEST_CASE("Singular element Duffy reference tensors are permutation covariant",
          "[singularelements][Serial]")
{
  constexpr double nu = 2.0 / 3.0;
  constexpr int order = 45;
  constexpr double radial_power = 6.0;
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  const auto check_pair = [&](const fem::singular::FirstOrderBasis &row,
                              const fem::singular::FirstOrderBasis &column, auto &&compute)
  {
    const auto reference = compute(row, column);
    const auto transpose = compute(column, row);
    const double reference_mass =
        fem::singular::ContractFirstOrderMass(reference, grad_lambda, 1.0);
    const double reference_curl =
        fem::singular::ContractFirstOrderCurlCurl(reference, grad_lambda, 1.0);
    for (int u = 0; u < 3; u++)
    {
      for (int v = 0; v < 3; v++)
      {
        CAPTURE(u, v);
        CHECK_THAT(transpose.mass[u][v], WithinAbs(reference.mass[v][u], 2.0e-15));
        CHECK_THAT(transpose.curl_curl[u][v],
                   WithinAbs(reference.curl_curl[v][u], 2.0e-15));
      }
    }

    std::array<int, 4> permutation{0, 1, 2, 3};
    do
    {
      CAPTURE(permutation);
      const auto permuted_row =
          fem::singular::ApplyBarycentricPermutation(row, permutation);
      const auto permuted_column =
          fem::singular::ApplyBarycentricPermutation(column, permutation);
      const auto integral = compute(permuted_row, permuted_column);
      CHECK(integral.mass == reference.mass);
      CHECK(integral.curl_curl == reference.curl_curl);

      const auto permuted_grad_lambda =
          fem::singular::ApplyBarycentricPermutation(grad_lambda, permutation);
      CHECK_THAT(fem::singular::ContractFirstOrderMass(integral, permuted_grad_lambda, 1.0),
                 WithinAbs(reference_mass, 2.0e-14));
      CHECK_THAT(
          fem::singular::ContractFirstOrderCurlCurl(integral, permuted_grad_lambda, 1.0),
          WithinAbs(reference_curl, 2.0e-14));
    } while (std::next_permutation(permutation.begin(), permutation.end()));
  };

  check_pair(fem::singular::MakeNodeRotational(2, 0, 3, nu),
             fem::singular::MakeStandardNedelec(1, 2),
             [&](const fem::singular::FirstOrderBasis &row,
                 const fem::singular::FirstOrderBasis &column)
             {
               return fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(
                   row, column, order, radial_power);
             });
  check_pair(fem::singular::MakeEdgeRotational(2, 0, 3, 1, nu),
             fem::singular::MakeStandardNedelec(1, 3),
             [&](const fem::singular::FirstOrderBasis &row,
                 const fem::singular::FirstOrderBasis &column)
             {
               return fem::singular::ComputeFirstOrderEdgeDuffyReferenceIntegral(
                   row, column, order, radial_power);
             });
}

TEST_CASE("Singular element input validation", "[singularelements][Serial]")
{
  const BarycentricPoint lambda{0.4, 0.2, 0.25, 0.15};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  CHECK_THROWS_AS(fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, 0.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, 1.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 0, 2, 0.5),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 2, 0.5),
      std::invalid_argument);

  auto bad_lambda = lambda;
  bad_lambda[0] += 0.1;
  CHECK_THROWS_AS(fem::singular::EvaluateNodeGradient(bad_lambda, grad_lambda, 0, 1, 0.5),
                  std::invalid_argument);

  auto bad_gradients = grad_lambda;
  bad_gradients[0][0] += 0.1;
  CHECK_THROWS_AS(fem::singular::EvaluateNodeGradient(lambda, bad_gradients, 0, 1, 0.5),
                  std::invalid_argument);

  const BarycentricPoint singular_node{1.0, 0.0, 0.0, 0.0};
  CHECK_THROWS_AS(
      fem::singular::EvaluateNodeGradient(singular_node, grad_lambda, 0, 1, 0.5),
      std::domain_error);

  const BarycentricPoint singular_edge{0.4, 0.6, 0.0, 0.0};
  CHECK_THROWS_AS(
      fem::singular::EvaluateEdgeGradient(singular_edge, grad_lambda, 0, 1, 2, 0.5),
      std::domain_error);
}

TEST_CASE("Singular element reference tetrahedron quadrature", "[singularelements][Serial]")
{
  constexpr int order = 8;
  const auto base_count = fem::singular::ReferenceTetrahedronQuadraturePointCount(order, 0);
  CHECK(base_count == 46);

  for (int subdivisions = 0; subdivisions <= 3; subdivisions++)
  {
    CAPTURE(subdivisions);
    std::size_t count = 0;
    double weight_sum = 0.0;
    fem::singular::ForEachReferenceTetrahedronQuadraturePoint(
        order, subdivisions,
        [&](const BarycentricPoint &lambda, double weight)
        {
          count++;
          weight_sum += weight;
          CHECK(weight > 0.0);
          for (double value : lambda)
          {
            CHECK(value > 0.0);
            CHECK(value < 1.0);
          }
        });
    CHECK(count ==
          fem::singular::ReferenceTetrahedronQuadraturePointCount(order, subdivisions));
    CHECK_THAT(weight_sum, WithinAbs(1.0 / 6.0, 2.0e-13));
    CHECK_THAT(fem::singular::IntegrateReferenceTetrahedron(
                   order, subdivisions, [](const BarycentricPoint &) { return 1.0; }),
               WithinAbs(1.0 / 6.0, 2.0e-15));
    const auto vector_integral = fem::singular::IntegrateReferenceTetrahedron(
        order, subdivisions, 2,
        [](const BarycentricPoint &lambda, std::vector<double> &value)
        {
          value[0] = 1.0;
          value[1] = lambda[0];
        });
    REQUIRE(vector_integral.size() == 2);
    CHECK_THAT(vector_integral[0], WithinAbs(1.0 / 6.0, 2.0e-15));
    CHECK_THAT(vector_integral[1], WithinAbs(1.0 / 24.0, 2.0e-15));
  }

  std::array<double, 8> leaf_weights{};
  std::size_t point = 0;
  fem::singular::ForEachReferenceTetrahedronQuadraturePoint(
      order, 1,
      [&](const BarycentricPoint &, double weight)
      {
        leaf_weights[point / base_count] += weight;
        point++;
      });
  for (double weight : leaf_weights)
  {
    CHECK_THAT(weight, WithinAbs(1.0 / 48.0, 2.0e-15));
  }
}

TEST_CASE("Singular element quadrature polynomial exactness", "[singularelements][Serial]")
{
  constexpr int order = 8;
  for (int subdivisions : {0, 2})
  {
    for (int a = 0; a <= order; a++)
    {
      for (int b = 0; b <= order - a; b++)
      {
        for (int c = 0; c <= order - a - b; c++)
        {
          for (int d = 0; d <= order - a - b - c; d++)
          {
            CAPTURE(subdivisions, a, b, c, d);
            const std::array<int, 4> exponents{a, b, c, d};
            const double integral = fem::singular::IntegrateReferenceTetrahedron(
                order, subdivisions,
                [&](const BarycentricPoint &lambda)
                {
                  double value = 1.0;
                  for (int i = 0; i < 4; i++)
                  {
                    value *= std::pow(lambda[i], exponents[i]);
                  }
                  return value;
                });
            CHECK_THAT(integral,
                       WithinAbs(ExactBarycentricMonomialIntegral(exponents), 2.0e-14));
          }
        }
      }
    }
  }
}

TEST_CASE("Singular element quadrature resolves integrable powers",
          "[singularelements][Serial]")
{
  constexpr int order = 8;
  constexpr double nu = 0.5;
  constexpr double power = 2.0 * nu - 2.0;
  constexpr double exact_node = 1.0 / (2.0 * (2.0 * nu + 1.0));
  constexpr double exact_edge = 1.0 / (2.0 * nu * (2.0 * nu + 1.0));

  std::array<double, 7> node_error{};
  std::array<double, 7> edge_error{};
  for (int subdivisions = 0; subdivisions <= 6; subdivisions++)
  {
    const double node = fem::singular::IntegrateReferenceTetrahedron(
        order, subdivisions,
        [](const BarycentricPoint &lambda) { return std::pow(1.0 - lambda[0], power); });
    const double edge = fem::singular::IntegrateReferenceTetrahedron(
        order, subdivisions, [](const BarycentricPoint &lambda)
        { return std::pow(lambda[2] + lambda[3], power); });
    node_error[subdivisions] = std::abs(node - exact_node);
    edge_error[subdivisions] = std::abs(edge - exact_edge);
  }
  CAPTURE(node_error, edge_error);
  CHECK(node_error[6] < node_error[4]);
  CHECK(node_error[4] < node_error[2]);
  CHECK(node_error[2] < node_error[0]);
  CHECK(edge_error[6] < edge_error[4]);
  CHECK(edge_error[4] < edge_error[2]);
  CHECK(edge_error[2] < edge_error[0]);
  CHECK(node_error[6] < 2.0e-7);
  CHECK(edge_error[6] < 3.0e-4);

  std::array<double, 4> node_integrals{};
  for (int i = 0; i < 4; i++)
  {
    node_integrals[i] = fem::singular::IntegrateReferenceTetrahedron(
        order, 4,
        [=](const BarycentricPoint &lambda) { return std::pow(1.0 - lambda[i], power); });
    CHECK_THAT(node_integrals[i], WithinAbs(exact_node, 2.0e-6));
  }
  const auto [node_min, node_max] =
      std::minmax_element(node_integrals.begin(), node_integrals.end());
  CHECK_THAT(*node_max - *node_min, WithinAbs(0.0, 1.0e-10));

  std::array<double, 6> edge_integrals{};
  int edge_index = 0;
  for (int i = 0; i < 4; i++)
  {
    for (int j = i + 1; j < 4; j++)
    {
      edge_integrals[edge_index++] = fem::singular::IntegrateReferenceTetrahedron(
          order, 4, [=](const BarycentricPoint &lambda)
          { return std::pow(1.0 - lambda[i] - lambda[j], power); });
    }
  }
  for (double integral : edge_integrals)
  {
    CHECK_THAT(integral, WithinAbs(exact_edge, 1.0e-3));
  }

  // Equivalent feature orientations must map to one canonical reference
  // integral. Generating a separate table for each raw label would retain the
  // fixed midpoint-octahedron diagonal's finite-depth quadrature error.
  const double canonical_edge = edge_integrals[0];
  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    const double permuted_edge = fem::singular::IntegrateReferenceTetrahedron(
        order, 4,
        [&](const BarycentricPoint &canonical_lambda)
        {
          const auto lambda =
              fem::singular::ApplyBarycentricPermutation(canonical_lambda, permutation);
          return std::pow(1.0 - lambda[permutation[0]] - lambda[permutation[1]], power);
        });
    CHECK_THAT(permuted_edge, WithinAbs(canonical_edge, 1.0e-14));
  } while (std::next_permutation(permutation.begin(), permutation.end()));
}

TEST_CASE("Singular element feature-aligned Duffy quadrature", "[singularelements][Serial]")
{
  constexpr int order = 45;
  constexpr double radial_power = 6.0;
  const auto check_rule = [](auto &&visit)
  {
    double weight_sum = 0.0;
    std::size_t count = 0;
    visit(
        [&](const BarycentricPoint &lambda, double weight)
        {
          count++;
          weight_sum += weight;
          CHECK(weight > 0.0);
          for (double value : lambda)
          {
            CHECK(value > 0.0);
            CHECK(value < 1.0);
          }
        });
    CHECK(count > 0);
    CHECK_THAT(weight_sum, WithinAbs(1.0 / 6.0, 2.0e-15));
  };
  check_rule(
      [](const fem::singular::QuadraturePointVisitor &visitor)
      {
        fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
            order, 0, radial_power, visitor);
      });
  check_rule(
      [](const fem::singular::QuadraturePointVisitor &visitor)
      {
        fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
            order, 0, 1, radial_power, visitor);
      });

  for (double nu : {0.5, 2.0 / 3.0, 0.731})
  {
    CAPTURE(nu);
    const double power = 2.0 * nu - 2.0;
    const double exact_node = 1.0 / (2.0 * (2.0 * nu + 1.0));
    const double exact_edge = 1.0 / (2.0 * nu * (2.0 * nu + 1.0));
    const double node = fem::singular::IntegrateReferenceTetrahedronNodeDuffy(
        order, 0, radial_power, [=](const BarycentricPoint &lambda)
        { return std::pow(fem::singular::NodeRadialCoordinate(lambda, 0), power); });
    const double edge = fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
        order, 0, 1, radial_power, [=](const BarycentricPoint &lambda)
        { return std::pow(fem::singular::EdgeRadialCoordinate(lambda, 0, 1), power); });
    CHECK_THAT(node, WithinAbs(exact_node, 2.0e-8));
    CHECK_THAT(edge, WithinAbs(exact_edge, 2.0e-8));
  }

  constexpr double nu = 0.5;
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const auto node_energy = [&](const BarycentricPoint &lambda)
  {
    const auto value =
        fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu).value;
    return Dot(value, value);
  };
  const auto edge_energy = [&](const BarycentricPoint &lambda)
  {
    const auto value =
        fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu).value;
    return Dot(value, value);
  };
  const double node_reference = fem::singular::IntegrateReferenceTetrahedronNodeDuffy(
      order, 0, radial_power, node_energy);
  const double edge_reference = fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
      order, 0, 1, radial_power, edge_energy);
  CHECK_THAT(fem::singular::IntegrateReferenceTetrahedronNodeDuffy(43, 0, radial_power,
                                                                   node_energy),
             WithinAbs(node_reference, 2.0e-10));
  CHECK_THAT(fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(43, 0, 1, radial_power,
                                                                   edge_energy),
             WithinAbs(edge_reference, 2.0e-10));

  std::array<double, 6> node_error;
  std::array<double, 6> edge_error;
  for (int subdivisions = 0; subdivisions <= 5; subdivisions++)
  {
    node_error[subdivisions] = std::abs(
        fem::singular::IntegrateReferenceTetrahedron(8, subdivisions, node_energy) -
        node_reference);
    edge_error[subdivisions] = std::abs(
        fem::singular::IntegrateReferenceTetrahedron(8, subdivisions, edge_energy) -
        edge_reference);
  }
  CAPTURE(node_reference, edge_reference, node_error, edge_error);
  CHECK(node_error[5] < node_error[3]);
  CHECK(node_error[3] < node_error[1]);
  CHECK(edge_error[5] < edge_error[3]);
  CHECK(edge_error[3] < edge_error[1]);
}

TEST_CASE("Singular element adaptive tetrahedron quadrature", "[singularelements][Serial]")
{
  const auto polynomial = fem::singular::IntegrateReferenceTetrahedronAdaptive(
      8, 1.0e-13, 1.0e-13, 4,
      [](const BarycentricPoint &lambda) { return lambda[0] * lambda[1] * lambda[2]; });
  CHECK(polynomial.converged);
  CHECK(polynomial.leaf_count == 8);
  CHECK(polynomial.maximum_subdivision_depth == 1);
  CHECK_THAT(polynomial.value,
             WithinAbs(ExactBarycentricMonomialIntegral({1, 1, 1, 0}), 2.0e-15));
  CHECK(polynomial.estimated_absolute_error < 2.0e-15);

  constexpr double nu = 0.5;
  constexpr double radial_power = 6.0;
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const auto node_energy = [&](const BarycentricPoint &lambda)
  {
    const auto value =
        fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu).value;
    return Dot(value, value);
  };
  const auto edge_energy = [&](const BarycentricPoint &lambda)
  {
    const auto value =
        fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu).value;
    return Dot(value, value);
  };
  const double node_reference = fem::singular::IntegrateReferenceTetrahedronNodeDuffy(
      45, 0, radial_power, node_energy);
  const double edge_reference = fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
      45, 0, 1, radial_power, edge_energy);
  const auto node = fem::singular::IntegrateReferenceTetrahedronAdaptive(8, 1.0e-6, 1.0e-6,
                                                                         9, node_energy);
  const auto edge = fem::singular::IntegrateReferenceTetrahedronAdaptive(8, 1.0e-6, 1.0e-6,
                                                                         9, edge_energy);
  CAPTURE(node_reference, edge_reference, node.value, edge.value,
          node.estimated_absolute_error, edge.estimated_absolute_error, node.leaf_count,
          edge.leaf_count, node.maximum_subdivision_depth, edge.maximum_subdivision_depth);
  CHECK(node.converged);
  CHECK(edge.converged);
  CHECK(std::abs(node.value - node_reference) <=
        2.0 * node.estimated_absolute_error + 1.0e-12);
  CHECK(std::abs(edge.value - edge_reference) <=
        2.0 * edge.estimated_absolute_error + 1.0e-12);

  const auto vector = fem::singular::IntegrateReferenceTetrahedronAdaptive(
      8, 1.0e-6, 1.0e-6, 9, 3,
      [&](const BarycentricPoint &lambda, std::vector<double> &value)
      {
        value[0] = lambda[0] * lambda[1] * lambda[2];
        value[1] = node_energy(lambda);
        value[2] = edge_energy(lambda);
      });
  REQUIRE(vector.value.size() == 3);
  REQUIRE(vector.estimated_absolute_error.size() == 3);
  CHECK(vector.converged);
  CHECK(vector.maximum_subdivision_depth >=
        std::max(node.maximum_subdivision_depth, edge.maximum_subdivision_depth));
  CHECK_THAT(vector.value[0],
             WithinAbs(ExactBarycentricMonomialIntegral({1, 1, 1, 0}), 2.0e-15));
  CHECK(std::abs(vector.value[1] - node_reference) <=
        2.0 * vector.estimated_absolute_error[1] + 1.0e-12);
  CHECK(std::abs(vector.value[2] - edge_reference) <=
        2.0 * vector.estimated_absolute_error[2] + 1.0e-12);

  const auto depth_limited =
      fem::singular::IntegrateReferenceTetrahedronAdaptive(8, 1.0e-14, 0.0, 1, edge_energy);
  CHECK_FALSE(depth_limited.converged);
  CHECK(depth_limited.leaf_count == 8);
  CHECK(depth_limited.maximum_subdivision_depth == 1);
  CHECK(depth_limited.estimated_absolute_error > 1.0e-14);
}

TEST_CASE("Singular element adaptive quadrature resolves mixed features",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  constexpr double absolute_tolerance = 2.0e-6;
  constexpr double relative_tolerance = 2.0e-6;
  constexpr int max_subdivisions = 9;
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const std::array mixed_pairs{std::pair{fem::singular::MakeNodeGradient(0, 1, nu),
                                         fem::singular::MakeNodeGradient(1, 0, nu)},
                               std::pair{fem::singular::MakeNodeGradient(0, 1, nu),
                                         fem::singular::MakeEdgeGradient(0, 1, 2, nu)},
                               std::pair{fem::singular::MakeEdgeGradient(0, 1, 2, nu),
                                         fem::singular::MakeEdgeGradient(2, 3, 0, nu)}};

  for (std::size_t pair = 0; pair < mixed_pairs.size(); pair++)
  {
    CAPTURE(pair);
    const auto row_basis = mixed_pairs[pair].first;
    const auto column_basis = mixed_pairs[pair].second;
    const auto integrand = [&](const BarycentricPoint &lambda)
    {
      const auto row =
          fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda, row_basis);
      const auto column =
          fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda, column_basis);
      return Dot(row.value, column.value);
    };
    const auto order_8 = fem::singular::IntegrateReferenceTetrahedronAdaptive(
        8, absolute_tolerance, relative_tolerance, max_subdivisions, integrand);
    const auto order_10 = fem::singular::IntegrateReferenceTetrahedronAdaptive(
        10, absolute_tolerance, relative_tolerance, max_subdivisions, integrand);
    CAPTURE(order_8.value, order_10.value, order_8.estimated_absolute_error,
            order_10.estimated_absolute_error, order_8.leaf_count, order_10.leaf_count,
            order_8.maximum_subdivision_depth, order_10.maximum_subdivision_depth);
    CHECK(order_8.converged);
    CHECK(order_10.converged);
    CHECK(std::abs(order_8.value - order_10.value) <=
          2.0 * (order_8.estimated_absolute_error + order_10.estimated_absolute_error) +
              1.0e-12);
  }

  const auto node_rotational = fem::singular::MakeNodeRotational(0, 1, 3, nu);
  const auto edge_rotational = fem::singular::MakeEdgeRotational(1, 2, 0, 3, nu);
  const auto tensor_order_8 = fem::singular::ComputeFirstOrderAdaptiveReferenceIntegral(
      node_rotational, edge_rotational, 8, absolute_tolerance, relative_tolerance,
      max_subdivisions);
  const auto tensor_order_10 = fem::singular::ComputeFirstOrderAdaptiveReferenceIntegral(
      node_rotational, edge_rotational, 10, absolute_tolerance, relative_tolerance,
      max_subdivisions);
  CHECK(tensor_order_8.integral.quadrature_rule ==
        fem::singular::ReferenceQuadratureRule::ADAPTIVE_TETRAHEDRON);
  CHECK(tensor_order_8.integral.quadrature_order == 8);
  CHECK(tensor_order_8.integral.subdivisions == max_subdivisions);
  CHECK(tensor_order_8.integral.radial_power == 0.0);
  CHECK(tensor_order_8.absolute_tolerance == absolute_tolerance);
  CHECK(tensor_order_8.relative_tolerance == relative_tolerance);
  CHECK(tensor_order_8.total_leaf_count > 0);
  CHECK(tensor_order_8.maximum_subdivision_depth <= max_subdivisions);
  CHECK(tensor_order_8.converged);
  CHECK(tensor_order_10.converged);
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      CAPTURE(u, v, tensor_order_8.integral.mass[u][v], tensor_order_10.integral.mass[u][v],
              tensor_order_8.mass_estimated_absolute_error[u][v],
              tensor_order_10.mass_estimated_absolute_error[u][v],
              tensor_order_8.integral.curl_curl[u][v],
              tensor_order_10.integral.curl_curl[u][v],
              tensor_order_8.curl_curl_estimated_absolute_error[u][v],
              tensor_order_10.curl_curl_estimated_absolute_error[u][v]);
      CHECK(std::abs(tensor_order_8.integral.mass[u][v] -
                     tensor_order_10.integral.mass[u][v]) <=
            2.0 * (tensor_order_8.mass_estimated_absolute_error[u][v] +
                   tensor_order_10.mass_estimated_absolute_error[u][v]) +
                1.0e-12);
      CHECK(std::abs(tensor_order_8.integral.curl_curl[u][v] -
                     tensor_order_10.integral.curl_curl[u][v]) <=
            2.0 * (tensor_order_8.curl_curl_estimated_absolute_error[u][v] +
                   tensor_order_10.curl_curl_estimated_absolute_error[u][v]) +
                1.0e-12);
    }
  }
}

TEST_CASE("Singular element quadrature input validation", "[singularelements][Serial]")
{
  CHECK_THROWS_AS(fem::singular::ReferenceTetrahedronQuadraturePointCount(0, 0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ReferenceTetrahedronQuadraturePointCount(8, -1),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ReferenceTetrahedronQuadraturePointCount(8, 9),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ApplyBarycentricPermutation(
                      BarycentricPoint{0.25, 0.25, 0.25, 0.25},
                      fem::singular::BarycentricPermutation{0, 1, 1, 3}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronQuadraturePoint(
                      8, 0, fem::singular::QuadraturePointVisitor{}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                      8, 4, 1.0, [](const BarycentricPoint &, double) {}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
                      8, 0, 0, 1.0, [](const BarycentricPoint &, double) {}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                      8, 0, 1.0, fem::singular::QuadraturePointVisitor{}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                      8, 0, 0.0, [](const BarycentricPoint &, double) {}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                      63, 0, 6.0, [](const BarycentricPoint &, double) {}),
                  std::invalid_argument);
  CHECK_NOTHROW(fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
      63, 0, 1, 6.0, [](const BarycentricPoint &, double) {}));
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedron(
                      8, 0, fem::singular::ReferenceIntegrand{}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
                      8, 0, 1, 1.0, fem::singular::ReferenceIntegrand{}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, 0.0, 0.0, 8, [](const BarycentricPoint &) { return 1.0; }),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, -1.0, 1.0e-6, 8, [](const BarycentricPoint &) { return 1.0; }),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, 1.0e-6, 1.0e-6, 0, [](const BarycentricPoint &) { return 1.0; }),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, 1.0e-6, 1.0e-6, 13, [](const BarycentricPoint &) { return 1.0; }),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, 1.0e-6, 1.0e-6, 8, fem::singular::ReferenceIntegrand{}),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::IntegrateReferenceTetrahedronAdaptive(
          8, 1.0e-6, 1.0e-6, 8, 0, [](const BarycentricPoint &, std::vector<double> &) {}),
      std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedronAdaptive(
                      8, 1.0e-6, 1.0e-6, 8, 2, fem::singular::ReferenceVectorIntegrand{}),
                  std::invalid_argument);
  const auto standard = fem::singular::MakeStandardNedelec(0, 1);
  const auto node_0 = fem::singular::MakeNodeGradient(0, 1, 0.5);
  const auto node_1 = fem::singular::MakeNodeGradient(1, 0, 0.5);
  const auto edge_01 = fem::singular::MakeEdgeGradient(0, 1, 2, 0.5);
  const auto edge_02 = fem::singular::MakeEdgeGradient(0, 2, 1, 0.5);
  CHECK_THROWS_AS(fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(
                      standard, standard, 21, 6.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(node_0, node_1, 21, 6.0),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(node_0, edge_01, 21, 6.0),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::ComputeFirstOrderEdgeDuffyReferenceIntegral(edge_01, edge_02, 21, 6.0),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::ComputeFirstOrderEdgeDuffyReferenceIntegral(edge_01, node_0, 21, 6.0),
      std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTetrahedron(
                      8, 0, [](const BarycentricPoint &)
                      { return std::numeric_limits<double>::infinity(); }),
                  std::domain_error);
}

TEST_CASE("Singular element generic reference tensors support higher-order bases",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  constexpr int order = 8;
  constexpr int subdivisions = 2;
  const fem::singular::HigherOrderBasis node_gradient{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT,
      {0, 1, 2, 3},
      {2, 1, 1, 0},
      3,
      nu};
  const fem::singular::HigherOrderBasis node_rotational{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
      {0, 1, 2, 3},
      {1, 1, 2, 1},
      3,
      nu};
  const fem::singular::ReferenceBasis standard = fem::singular::MakeStandardNedelec(1, 2);
  const fem::singular::ReferenceBasis singular = node_gradient;

  const auto mixed =
      fem::singular::ComputeReferenceIntegral(standard, singular, order, subdivisions);
  const auto direct = DirectReferencePairIntegral(
      standard, singular,
      [&](const fem::singular::QuadraturePointVisitor &visitor)
      {
        fem::singular::ForEachReferenceTetrahedronQuadraturePoint(order, subdivisions,
                                                                  visitor);
      });
  CHECK(mixed.quadrature_rule ==
        fem::singular::ReferenceQuadratureRule::RECURSIVE_TETRAHEDRON);
  CHECK(mixed.quadrature_order == order);
  CHECK(mixed.subdivisions == subdivisions);
  CHECK(mixed.radial_power == 0.0);
  CHECK_THAT(TensorDifferenceMaxNorm(mixed.mass, direct.mass), WithinAbs(0.0, 5.0e-15));
  CHECK_THAT(TensorDifferenceMaxNorm(mixed.curl_curl, direct.curl_curl),
             WithinAbs(0.0, 5.0e-15));

  const fem::singular::HigherOrderBasis first_order_node{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT,
      {0, 1, 2, 3},
      {1, 1, 0, 0},
      1,
      nu};
  const auto generic_first_order = fem::singular::ComputeReferenceIntegral(
      first_order_node, first_order_node, order, subdivisions);
  const auto legacy_first_order = fem::singular::ComputeFirstOrderReferenceIntegral(
      fem::singular::MakeNodeGradient(0, 1, nu), fem::singular::MakeNodeGradient(0, 1, nu),
      order, subdivisions);
  CHECK_THAT(TensorDifferenceMaxNorm(generic_first_order.mass, legacy_first_order.mass),
             WithinAbs(0.0, 0.0));
  CHECK_THAT(
      TensorDifferenceMaxNorm(generic_first_order.curl_curl, legacy_first_order.curl_curl),
      WithinAbs(0.0, 0.0));

  const auto reference = fem::singular::ComputeReferenceIntegral(
      node_gradient, node_rotational, order, subdivisions);
  std::array<int, 4> permutation{0, 1, 2, 3};
  do
  {
    const auto permuted_row = fem::singular::ApplyBarycentricPermutation(
        fem::singular::ReferenceBasis{node_gradient}, permutation);
    const auto permuted_column = fem::singular::ApplyBarycentricPermutation(
        fem::singular::ReferenceBasis{node_rotational}, permutation);
    const auto permuted = fem::singular::ComputeReferenceIntegral(
        permuted_row, permuted_column, order, subdivisions);
    CAPTURE(permutation);
    CHECK_THAT(TensorDifferenceMaxNorm(permuted.mass, reference.mass), WithinAbs(0.0, 0.0));
    CHECK_THAT(TensorDifferenceMaxNorm(permuted.curl_curl, reference.curl_curl),
               WithinAbs(0.0, 0.0));
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  double jacobian_determinant;
  const auto physical_gradients = AffineBarycentricGradients(
      {1.7, 0.2, -0.1}, {0.3, 1.4, 0.25}, {-0.2, 0.1, 0.9}, jacobian_determinant);
  const auto physical_direct = DirectAffinePairIntegral(
      fem::singular::ReferenceBasis{node_gradient},
      fem::singular::ReferenceBasis{node_rotational}, physical_gradients,
      jacobian_determinant, order, subdivisions);
  CHECK_THAT(
      fem::singular::ContractMass(reference, physical_gradients, jacobian_determinant),
      WithinAbs(physical_direct[0], 2.0e-13));
  CHECK_THAT(
      fem::singular::ContractCurlCurl(reference, physical_gradients, jacobian_determinant),
      WithinAbs(physical_direct[1], 2.0e-13));
}

TEST_CASE("Singular element higher-order reference tensors pass Duffy and adaptive checks",
          "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  constexpr int duffy_order = 45;
  constexpr double radial_power = 6.0;
  const fem::singular::HigherOrderBasis node_rotational{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
      {0, 1, 2, 3},
      {1, 1, 2, 1},
      3,
      nu};
  const fem::singular::HigherOrderBasis edge_rotational{
      fem::singular::HigherOrderBasisFamily::EDGE_ROTATIONAL,
      {0, 1, 2, 3},
      {2, 2, 1, 1},
      3,
      nu};

  const auto check_pair = [&](const fem::singular::ReferenceBasis &basis, bool node_aligned)
  {
    const auto canonical = fem::singular::CanonicalizeReferenceBasisPair(basis, basis);
    const auto direct = DirectReferencePairIntegral(
        basis, basis,
        [&](const fem::singular::QuadraturePointVisitor &visitor)
        {
          if (node_aligned)
          {
            fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                duffy_order, CanonicalNodeFeature(canonical), radial_power, visitor);
          }
          else
          {
            const auto edge = CanonicalEdgeFeature(canonical);
            fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
                duffy_order, edge[0], edge[1], radial_power, visitor);
          }
        });
    const auto duffy = node_aligned ? fem::singular::ComputeNodeDuffyReferenceIntegral(
                                          basis, basis, duffy_order, radial_power)
                                    : fem::singular::ComputeEdgeDuffyReferenceIntegral(
                                          basis, basis, duffy_order, radial_power);
    const auto lower_order_duffy = node_aligned
                                       ? fem::singular::ComputeNodeDuffyReferenceIntegral(
                                             basis, basis, duffy_order - 2, radial_power)
                                       : fem::singular::ComputeEdgeDuffyReferenceIntegral(
                                             basis, basis, duffy_order - 2, radial_power);
    CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, direct.mass), WithinAbs(0.0, 5.0e-14));
    CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, direct.curl_curl),
               WithinAbs(0.0, 5.0e-14));
    CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, lower_order_duffy.mass),
               WithinAbs(0.0, 2.0e-9));
    CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, lower_order_duffy.curl_curl),
               WithinAbs(0.0, 2.0e-9));

    const auto coarse = fem::singular::ComputeReferenceIntegral(basis, basis, 8, 2);
    const auto fine = fem::singular::ComputeReferenceIntegral(basis, basis, 8, 4);
    CheckEntrywiseReferenceConvergence(coarse, fine, direct);

    const auto adaptive =
        fem::singular::ComputeAdaptiveReferenceIntegral(basis, basis, 8, 2.0e-6, 2.0e-6, 9);
    CHECK(adaptive.converged);
    for (int u = 0; u < 3; u++)
    {
      for (int v = 0; v < 3; v++)
      {
        CAPTURE(node_aligned, u, v, adaptive.integral.mass[u][v], duffy.mass[u][v],
                adaptive.mass_estimated_absolute_error[u][v],
                adaptive.integral.curl_curl[u][v], duffy.curl_curl[u][v],
                adaptive.curl_curl_estimated_absolute_error[u][v]);
        CHECK(std::abs(adaptive.integral.mass[u][v] - duffy.mass[u][v]) <=
              2.0 * adaptive.mass_estimated_absolute_error[u][v] + 2.0e-11);
        CHECK(std::abs(adaptive.integral.curl_curl[u][v] - duffy.curl_curl[u][v]) <=
              2.0 * adaptive.curl_curl_estimated_absolute_error[u][v] + 2.0e-11);
      }
    }
  };

  check_pair(node_rotational, true);
  check_pair(edge_rotational, false);
  CHECK_THROWS_AS(fem::singular::ComputeNodeDuffyReferenceIntegral(
                      node_rotational, edge_rotational, duffy_order, radial_power),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ComputeEdgeDuffyReferenceIntegral(
                      edge_rotational, node_rotational, duffy_order, radial_power),
                  std::invalid_argument);
}

TEST_CASE("Partitioned Duffy tensors resolve distinct singular features",
          "[singularelements][Serial]")
{
  constexpr int comparison_order = 39;
  constexpr int reference_order = 47;
  constexpr int higher_order = 49;
  constexpr double radial_power = 6.0;
  constexpr double adaptive_tolerance = 2.0e-6;
  constexpr int adaptive_max_subdivisions = 9;

  for (double nu : {0.5, 2.0 / 3.0, 0.87165534177406845})
  {
    CAPTURE(nu);
    const std::array pairs{
        std::pair{fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(0, 1, nu)},
                  fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(1, 0, nu)}},
        std::pair{
            fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(0, 2, nu)},
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(0, 1, 2, nu)}},
        std::pair{
            fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(0, 1, nu)},
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(1, 2, 3, nu)}},
        std::pair{
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(0, 1, 2, nu)},
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(0, 2, 1, nu)}},
        std::pair{
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(0, 1, 2, nu)},
            fem::singular::ReferenceBasis{fem::singular::MakeEdgeGradient(2, 3, 0, nu)}}};

    for (std::size_t pair = 0; pair < pairs.size(); pair++)
    {
      CAPTURE(pair);
      const auto &[row, column] = pairs[pair];
      const auto comparison = fem::singular::ComputePartitionedDuffyReferenceIntegral(
          row, column, comparison_order, radial_power);
      const auto reference = fem::singular::ComputePartitionedDuffyReferenceIntegral(
          row, column, reference_order, radial_power);
      const auto higher = fem::singular::ComputePartitionedDuffyReferenceIntegral(
          row, column, higher_order, radial_power);
      const auto adaptive = fem::singular::ComputeAdaptiveReferenceIntegral(
          row, column, 10, adaptive_tolerance, adaptive_tolerance,
          adaptive_max_subdivisions);
      CHECK(reference.quadrature_rule ==
            fem::singular::ReferenceQuadratureRule::PARTITIONED_DUFFY);
      CHECK(reference.quadrature_order == reference_order);
      CHECK(reference.subdivisions == 0);
      CHECK(reference.radial_power == radial_power);
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          const double mass_estimate =
              fem::singular::H1DuffyErrorSafetyFactor *
                  std::abs(reference.mass[u][v] - comparison.mass[u][v]) +
              64.0 * std::numeric_limits<double>::epsilon() *
                  std::max({1.0, std::abs(reference.mass[u][v]),
                            std::abs(comparison.mass[u][v])});
          const double curl_estimate =
              fem::singular::H1DuffyErrorSafetyFactor *
                  std::abs(reference.curl_curl[u][v] - comparison.curl_curl[u][v]) +
              64.0 * std::numeric_limits<double>::epsilon() *
                  std::max({1.0, std::abs(reference.curl_curl[u][v]),
                            std::abs(comparison.curl_curl[u][v])});
          CAPTURE(u, v, comparison.mass[u][v], reference.mass[u][v], higher.mass[u][v],
                  adaptive.integral.mass[u][v],
                  adaptive.mass_estimated_absolute_error[u][v], mass_estimate,
                  comparison.curl_curl[u][v], reference.curl_curl[u][v],
                  higher.curl_curl[u][v], adaptive.integral.curl_curl[u][v],
                  adaptive.curl_curl_estimated_absolute_error[u][v], curl_estimate);
          CHECK(std::abs(reference.mass[u][v] - higher.mass[u][v]) <= mass_estimate);
          CHECK(std::abs(reference.curl_curl[u][v] - higher.curl_curl[u][v]) <=
                curl_estimate);
          CHECK(std::abs(reference.mass[u][v] - adaptive.integral.mass[u][v]) <=
                mass_estimate + 2.0 * adaptive.mass_estimated_absolute_error[u][v] +
                    1.0e-12);
          CHECK(std::abs(reference.curl_curl[u][v] - adaptive.integral.curl_curl[u][v]) <=
                curl_estimate + 2.0 * adaptive.curl_curl_estimated_absolute_error[u][v] +
                    1.0e-12);
        }
      }
    }
  }

  const auto node =
      fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(0, 1, 0.5)};
  const auto same_node =
      fem::singular::ReferenceBasis{fem::singular::MakeNodeGradient(0, 2, 0.5)};
  const auto standard =
      fem::singular::ReferenceBasis{fem::singular::MakeStandardH1Gradient(1)};
  CHECK_THROWS_AS(fem::singular::ComputePartitionedDuffyReferenceIntegral(
                      node, same_node, reference_order, radial_power),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ComputePartitionedDuffyReferenceIntegral(
                      standard, node, reference_order, radial_power),
                  std::invalid_argument);
}

TEST_CASE("First-order H1 Duffy table estimates bound an exhaustive higher-order audit",
          "[singularelements][Serial]")
{
  constexpr int comparison_order = 39;
  constexpr int reference_order = 47;
  constexpr int audit_order = 49;
  constexpr double radial_power = 6.0;

  const auto check_entry = [](double comparison, double reference, double audit)
  {
    const double estimate =
        fem::singular::H1DuffyErrorSafetyFactor * std::abs(reference - comparison) +
        64.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(reference), std::abs(comparison)});
    CAPTURE(comparison, reference, audit, estimate);
    CHECK(std::abs(reference - audit) <= estimate);
  };
  const auto same_feature = [](const fem::singular::HigherOrderBasis &row,
                               const fem::singular::HigherOrderBasis &column)
  {
    if (row.family != column.family)
    {
      return false;
    }
    if (row.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
    {
      return row.nodes[0] == column.nodes[0];
    }
    REQUIRE(row.family == fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT);
    std::array<int, 2> row_edge{row.nodes[0], row.nodes[1]};
    std::array<int, 2> column_edge{column.nodes[0], column.nodes[1]};
    std::sort(row_edge.begin(), row_edge.end());
    std::sort(column_edge.begin(), column_edge.end());
    return row_edge == column_edge;
  };

  for (double nu : {0.5, 2.0 / 3.0})
  {
    CAPTURE(nu);
    std::vector<fem::singular::HigherOrderBasis> bases;
    for (int node = 0; node < 4; node++)
    {
      std::array<int, 4> nodes{node, -1, -1, -1};
      int next = 1;
      for (int other = 0; other < 4; other++)
      {
        if (other != node)
        {
          nodes[next++] = other;
        }
      }
      const auto entries =
          fem::singular::EnumerateHigherOrderNodeGradientBases(nodes, 1, nu);
      bases.insert(bases.end(), entries.begin(), entries.end());
    }
    for (int first = 0; first < 4; first++)
    {
      for (int second = first + 1; second < 4; second++)
      {
        std::array<int, 4> nodes{first, second, -1, -1};
        int next = 2;
        for (int other = 0; other < 4; other++)
        {
          if (other != first && other != second)
          {
            nodes[next++] = other;
          }
        }
        const auto entries =
            fem::singular::EnumerateHigherOrderEdgeGradientBases(nodes, 1, nu);
        bases.insert(bases.end(), entries.begin(), entries.end());
      }
    }
    REQUIRE(bases.size() == 24);

    for (std::size_t row = 0; row < bases.size(); row++)
    {
      for (std::size_t column = row; column < bases.size(); column++)
      {
        CAPTURE(row, column, static_cast<int>(bases[row].family),
                static_cast<int>(bases[column].family), bases[row].nodes,
                bases[column].nodes);
        const auto integrate = [&](int order)
        {
          const auto row_basis = fem::singular::ReferenceBasis{bases[row]};
          const auto column_basis = fem::singular::ReferenceBasis{bases[column]};
          if (!same_feature(bases[row], bases[column]))
          {
            return fem::singular::ComputePartitionedDuffyReferenceIntegral(
                row_basis, column_basis, order, radial_power);
          }
          if (bases[row].family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
          {
            return fem::singular::ComputeNodeDuffyReferenceIntegral(row_basis, column_basis,
                                                                    order, radial_power);
          }
          return fem::singular::ComputeEdgeDuffyReferenceIntegral(row_basis, column_basis,
                                                                  order, radial_power);
        };
        const auto comparison = integrate(comparison_order);
        const auto reference = integrate(reference_order);
        const auto audit = integrate(audit_order);
        for (int u = 0; u < 3; u++)
        {
          for (int v = 0; v < 3; v++)
          {
            CAPTURE(u, v);
            check_entry(comparison.mass[u][v], reference.mass[u][v], audit.mass[u][v]);
          }
        }
      }
    }

    for (int standard = 0; standard < 4; standard++)
    {
      for (std::size_t singular = 0; singular < bases.size(); singular++)
      {
        CAPTURE(standard, singular, static_cast<int>(bases[singular].family),
                bases[singular].nodes);
        const auto integrate = [&](int order)
        {
          const auto standard_basis = fem::singular::ReferenceBasis{
              fem::singular::MakeStandardH1Gradient(standard)};
          const auto singular_basis = fem::singular::ReferenceBasis{bases[singular]};
          if (bases[singular].family ==
              fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
          {
            return fem::singular::ComputeNodeDuffyReferenceIntegral(
                standard_basis, singular_basis, order, radial_power);
          }
          return fem::singular::ComputeEdgeDuffyReferenceIntegral(
              standard_basis, singular_basis, order, radial_power);
        };
        const auto comparison = integrate(comparison_order);
        const auto reference = integrate(reference_order);
        const auto audit = integrate(audit_order);
        for (int u = 0; u < 3; u++)
        {
          for (int v = 0; v < 3; v++)
          {
            CAPTURE(u, v);
            check_entry(comparison.mass[u][v], reference.mass[u][v], audit.mass[u][v]);
          }
        }
      }
    }
  }
}

TEST_CASE("Singular element reference tensors converge to feature-aligned quadrature",
          "[singularelements][Serial]")
{
  constexpr int recursive_order = 8;
  constexpr int coarse_subdivisions = 3;
  constexpr int fine_subdivisions = 5;
  constexpr int duffy_order = 45;
  constexpr double radial_power = 6.0;

  for (double nu : {0.5, 2.0 / 3.0})
  {
    CAPTURE(nu);
    const auto node_gradient = fem::singular::MakeNodeGradient(0, 1, nu);
    const auto node_rotational = fem::singular::MakeNodeRotational(0, 1, 2, nu);
    const std::array node_pairs{
        std::pair{fem::singular::MakeStandardH1Gradient(1), node_gradient},
        std::pair{node_gradient, node_gradient},
        std::pair{fem::singular::MakeStandardNedelec(1, 2), node_rotational},
        std::pair{node_rotational, node_rotational}};
    for (std::size_t pair = 0; pair < node_pairs.size(); pair++)
    {
      CAPTURE(pair);
      const auto &[row, column] = node_pairs[pair];
      const auto canonical = fem::singular::CanonicalizeFirstOrderBasisPair(row, column);
      const int singular_node = CanonicalNodeFeature(canonical);
      const auto reference = DirectReferencePairIntegral(
          row, column,
          [&](const fem::singular::QuadraturePointVisitor &visitor)
          {
            fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
                duffy_order, singular_node, radial_power, visitor);
          });
      const auto duffy = fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(
          row, column, duffy_order, radial_power);
      const auto lower_order_duffy =
          fem::singular::ComputeFirstOrderNodeDuffyReferenceIntegral(
              row, column, duffy_order - 2, radial_power);
      CHECK(duffy.quadrature_rule == fem::singular::ReferenceQuadratureRule::NODE_DUFFY);
      CHECK(duffy.quadrature_order == duffy_order);
      CHECK(duffy.subdivisions == 0);
      CHECK(duffy.radial_power == radial_power);
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, reference.mass),
                 WithinAbs(0.0, 5.0e-15));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, reference.curl_curl),
                 WithinAbs(0.0, 5.0e-15));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, lower_order_duffy.mass),
                 WithinAbs(0.0, 2.0e-10));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, lower_order_duffy.curl_curl),
                 WithinAbs(0.0, 2.0e-10));
      const auto coarse = fem::singular::ComputeFirstOrderReferenceIntegral(
          row, column, recursive_order, coarse_subdivisions);
      const auto fine = fem::singular::ComputeFirstOrderReferenceIntegral(
          row, column, recursive_order, fine_subdivisions);
      CheckEntrywiseReferenceConvergence(coarse, fine, reference);
    }

    const auto edge_gradient = fem::singular::MakeEdgeGradient(0, 1, 2, nu);
    const auto edge_rotational = fem::singular::MakeEdgeRotational(0, 1, 2, 3, nu);
    const std::array edge_pairs{
        std::pair{fem::singular::MakeStandardH1Gradient(2), edge_gradient},
        std::pair{edge_gradient, edge_gradient},
        std::pair{fem::singular::MakeStandardNedelec(2, 3), edge_rotational},
        std::pair{edge_rotational, edge_rotational}};
    for (std::size_t pair = 0; pair < edge_pairs.size(); pair++)
    {
      CAPTURE(pair);
      const auto &[row, column] = edge_pairs[pair];
      const auto canonical = fem::singular::CanonicalizeFirstOrderBasisPair(row, column);
      const auto singular_edge = CanonicalEdgeFeature(canonical);
      const auto reference = DirectReferencePairIntegral(
          row, column,
          [&](const fem::singular::QuadraturePointVisitor &visitor)
          {
            fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
                duffy_order, singular_edge[0], singular_edge[1], radial_power, visitor);
          });
      const auto duffy = fem::singular::ComputeFirstOrderEdgeDuffyReferenceIntegral(
          row, column, duffy_order, radial_power);
      const auto lower_order_duffy =
          fem::singular::ComputeFirstOrderEdgeDuffyReferenceIntegral(
              row, column, duffy_order - 2, radial_power);
      CHECK(duffy.quadrature_rule == fem::singular::ReferenceQuadratureRule::EDGE_DUFFY);
      CHECK(duffy.quadrature_order == duffy_order);
      CHECK(duffy.subdivisions == 0);
      CHECK(duffy.radial_power == radial_power);
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, reference.mass),
                 WithinAbs(0.0, 5.0e-15));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, reference.curl_curl),
                 WithinAbs(0.0, 5.0e-15));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.mass, lower_order_duffy.mass),
                 WithinAbs(0.0, 2.0e-10));
      CHECK_THAT(TensorDifferenceMaxNorm(duffy.curl_curl, lower_order_duffy.curl_curl),
                 WithinAbs(0.0, 2.0e-10));
      const auto coarse = fem::singular::ComputeFirstOrderReferenceIntegral(
          row, column, recursive_order, coarse_subdivisions);
      const auto fine = fem::singular::ComputeFirstOrderReferenceIntegral(
          row, column, recursive_order, fine_subdivisions);
      CheckEntrywiseReferenceConvergence(coarse, fine, reference);
    }
  }
}

TEST_CASE("Singular element first-order basis descriptors", "[singularelements][Serial]")
{
  constexpr double nu = 0.61;
  const BarycentricPoint lambda{0.37, 0.21, 0.18, 0.24};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();

  CheckBasisValue(fem::singular::EvaluateFirstOrderBasis(
                      lambda, grad_lambda, fem::singular::MakeStandardH1Gradient(2)),
                  {grad_lambda[2], {0.0, 0.0, 0.0}}, 0.0);
  CheckBasisValue(fem::singular::EvaluateFirstOrderBasis(
                      lambda, grad_lambda, fem::singular::MakeStandardNedelec(0, 1)),
                  fem::singular::EvaluateStandardEdge(lambda, grad_lambda, 0, 1), 0.0);
  CheckBasisValue(fem::singular::EvaluateFirstOrderBasis(
                      lambda, grad_lambda, fem::singular::MakeNodeGradient(0, 1, nu)),
                  fem::singular::EvaluateNodeGradient(lambda, grad_lambda, 0, 1, nu), 0.0);
  CheckBasisValue(fem::singular::EvaluateFirstOrderBasis(
                      lambda, grad_lambda, fem::singular::MakeNodeRotational(0, 1, 2, nu)),
                  fem::singular::EvaluateNodeRotational(lambda, grad_lambda, 0, 1, 2, nu),
                  0.0);
  CheckBasisValue(fem::singular::EvaluateFirstOrderBasis(
                      lambda, grad_lambda, fem::singular::MakeEdgeGradient(0, 1, 2, nu)),
                  fem::singular::EvaluateEdgeGradient(lambda, grad_lambda, 0, 1, 2, nu),
                  0.0);
  CheckBasisValue(
      fem::singular::EvaluateFirstOrderBasis(
          lambda, grad_lambda, fem::singular::MakeEdgeRotational(0, 1, 2, 3, nu)),
      fem::singular::EvaluateEdgeRotational(lambda, grad_lambda, 0, 1, 2, 3, nu), 0.0);

  CHECK_THROWS_AS(fem::singular::MakeStandardH1Gradient(4), std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::MakeStandardNedelec(1, 1), std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::MakeNodeGradient(0, 1, 1.0), std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::MakeEdgeRotational(0, 1, 2, 2, nu), std::invalid_argument);

  const fem::singular::FirstOrderBasis unknown{
      static_cast<fem::singular::FirstOrderBasisFamily>(99), {0, 1, 2, 3}, nu};
  CHECK_THROWS_AS(fem::singular::EvaluateFirstOrderBasis(lambda, grad_lambda, unknown),
                  std::invalid_argument);
}

TEST_CASE("Singular element first-order reference coefficient integrals",
          "[singularelements][Serial]")
{
  constexpr int order = 8;
  constexpr int subdivisions = 3;
  constexpr double nu = 0.5;

  const auto h1_0 = fem::singular::MakeStandardH1Gradient(0);
  const auto h1_00 =
      fem::singular::ComputeFirstOrderReferenceIntegral(h1_0, h1_0, order, subdivisions);
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      CAPTURE(u, v);
      CHECK_THAT(h1_00.mass[u][v], WithinAbs(1.0 / 6.0, 2.0e-15));
      CHECK_THAT(h1_00.curl_curl[u][v], WithinAbs(0.0, 0.0));
    }
  }

  const auto standard = fem::singular::MakeStandardNedelec(0, 1);
  const auto standard_integral = fem::singular::ComputeFirstOrderReferenceIntegral(
      standard, standard, order, subdivisions);
  const auto standard_value = fem::singular::EvaluateFirstOrderBasis(
      BarycentricPoint{0.25, 0.25, 0.25, 0.25},
      fem::singular::ReferenceBarycentricGradients(), standard);
  for (int u = 0; u < 3; u++)
  {
    for (int v = 0; v < 3; v++)
    {
      CAPTURE(u, v);
      CHECK_THAT(standard_integral.curl_curl[u][v],
                 WithinAbs(standard_value.curl[u] * standard_value.curl[v] / 6.0, 2.0e-15));
    }
  }

  const std::array<fem::singular::FirstOrderBasis, 6> basis{
      h1_0,
      standard,
      fem::singular::MakeNodeGradient(0, 1, nu),
      fem::singular::MakeNodeRotational(0, 1, 2, nu),
      fem::singular::MakeEdgeGradient(0, 1, 2, nu),
      fem::singular::MakeEdgeRotational(0, 1, 2, 3, nu)};
  for (std::size_t i = 0; i < basis.size(); i++)
  {
    for (std::size_t j = i; j < basis.size(); j++)
    {
      CAPTURE(i, j);
      const auto ij = fem::singular::ComputeFirstOrderReferenceIntegral(
          basis[i], basis[j], order, subdivisions);
      const auto ji = fem::singular::ComputeFirstOrderReferenceIntegral(
          basis[j], basis[i], order, subdivisions);
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          CHECK_THAT(ij.mass[u][v], WithinAbs(ji.mass[v][u], 2.0e-15));
          CHECK_THAT(ij.curl_curl[u][v], WithinAbs(ji.curl_curl[v][u], 2.0e-15));
        }
      }
    }
  }
}

TEST_CASE("Singular element reference tensors factor affine geometry",
          "[singularelements][Serial]")
{
  constexpr int order = 8;
  constexpr int subdivisions = 3;
  constexpr double nu = 2.0 / 3.0;
  double jacobian_determinant;
  const auto grad_lambda = AffineBarycentricGradients(
      {1.7, 0.2, -0.1}, {0.3, 1.2, 0.4}, {-0.2, 0.1, 0.9}, jacobian_determinant);
  const std::array<fem::singular::FirstOrderBasis, 6> basis{
      fem::singular::MakeStandardH1Gradient(0),
      fem::singular::MakeStandardNedelec(0, 2),
      fem::singular::MakeNodeGradient(0, 1, nu),
      fem::singular::MakeNodeRotational(0, 1, 3, nu),
      fem::singular::MakeEdgeGradient(0, 1, 2, nu),
      fem::singular::MakeEdgeRotational(0, 1, 2, 3, nu)};

  for (std::size_t i = 0; i < basis.size(); i++)
  {
    for (std::size_t j = 0; j < basis.size(); j++)
    {
      CAPTURE(i, j);
      const auto integral = fem::singular::ComputeFirstOrderReferenceIntegral(
          basis[i], basis[j], order, subdivisions);
      const auto direct = DirectAffinePairIntegral(
          basis[i], basis[j], grad_lambda, jacobian_determinant, order, subdivisions);
      CHECK_THAT(fem::singular::ContractFirstOrderMass(integral, grad_lambda,
                                                       jacobian_determinant),
                 WithinAbs(direct[0], 2.0e-12));
      CHECK_THAT(fem::singular::ContractFirstOrderCurlCurl(integral, grad_lambda,
                                                           jacobian_determinant),
                 WithinAbs(direct[1], 2.0e-12));
    }
  }

  const auto integral =
      fem::singular::ComputeFirstOrderReferenceIntegral(basis[0], basis[0], order, 0);
  CHECK_THROWS_AS(fem::singular::ContractFirstOrderMass(integral, grad_lambda, -1.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ContractFirstOrderMass(integral, grad_lambda,
                                                        1.01 * jacobian_determinant),
                  std::invalid_argument);
}

TEST_CASE("Singular element first-order reference matrices are semidefinite",
          "[singularelements][Serial]")
{
  constexpr int order = 8;
  constexpr int subdivisions = 3;
  constexpr double nu = 0.5;
  const std::array<fem::singular::FirstOrderBasis, 6> basis{
      fem::singular::MakeStandardH1Gradient(1),
      fem::singular::MakeStandardNedelec(0, 1),
      fem::singular::MakeNodeGradient(0, 1, nu),
      fem::singular::MakeNodeRotational(0, 1, 2, nu),
      fem::singular::MakeEdgeGradient(0, 1, 2, nu),
      fem::singular::MakeEdgeRotational(0, 1, 2, 3, nu)};
  mfem::DenseMatrix mass(static_cast<int>(basis.size()));
  mfem::DenseMatrix curl_curl(static_cast<int>(basis.size()));
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  for (std::size_t i = 0; i < basis.size(); i++)
  {
    for (std::size_t j = i; j < basis.size(); j++)
    {
      const auto integral = fem::singular::ComputeFirstOrderReferenceIntegral(
          basis[i], basis[j], order, subdivisions);
      mass(static_cast<int>(i), static_cast<int>(j)) =
          mass(static_cast<int>(j), static_cast<int>(i)) =
              fem::singular::ContractFirstOrderMass(integral, grad_lambda, 1.0);
      curl_curl(static_cast<int>(i), static_cast<int>(j)) =
          curl_curl(static_cast<int>(j), static_cast<int>(i)) =
              fem::singular::ContractFirstOrderCurlCurl(integral, grad_lambda, 1.0);
    }
  }

  const auto mass_eigenvalues = SymmetricEigenvalues(mass);
  const auto curl_curl_eigenvalues = SymmetricEigenvalues(curl_curl);
  const double mass_tolerance = 2.0e-12 * mass.MaxMaxNorm();
  const double curl_curl_tolerance = 2.0e-12 * curl_curl.MaxMaxNorm();
  for (std::size_t i = 0; i < mass_eigenvalues.size(); i++)
  {
    CAPTURE(i, mass_eigenvalues[i]);
    CHECK(mass_eigenvalues[i] >= -mass_tolerance);
  }
  for (std::size_t i = 0; i < curl_curl_eigenvalues.size(); i++)
  {
    CAPTURE(i, curl_curl_eigenvalues[i]);
    CHECK(curl_curl_eigenvalues[i] >= -curl_curl_tolerance);
  }
}

TEST_CASE("Singular element first-order reference tensors converge with subdivision",
          "[singularelements][Serial]")
{
  constexpr int order = 8;
  constexpr double nu = 0.5;
  const auto node = fem::singular::MakeNodeGradient(0, 1, nu);
  const auto edge = fem::singular::MakeEdgeGradient(0, 1, 2, nu);

  std::array<fem::singular::FirstOrderReferenceIntegral, 6> node_integrals;
  std::array<fem::singular::FirstOrderReferenceIntegral, 6> edge_integrals;
  for (int subdivisions = 0; subdivisions <= 5; subdivisions++)
  {
    node_integrals[subdivisions] =
        fem::singular::ComputeFirstOrderReferenceIntegral(node, node, order, subdivisions);
    edge_integrals[subdivisions] =
        fem::singular::ComputeFirstOrderReferenceIntegral(edge, edge, order, subdivisions);
  }

  std::array<double, 5> node_difference;
  std::array<double, 5> edge_difference;
  for (int subdivisions = 0; subdivisions < 5; subdivisions++)
  {
    node_difference[subdivisions] = TensorDifferenceMaxNorm(
        node_integrals[subdivisions].mass, node_integrals[subdivisions + 1].mass);
    edge_difference[subdivisions] = TensorDifferenceMaxNorm(
        edge_integrals[subdivisions].mass, edge_integrals[subdivisions + 1].mass);
  }
  CAPTURE(node_difference, edge_difference);
  CHECK(node_difference[4] < node_difference[2]);
  CHECK(node_difference[2] < node_difference[0]);
  CHECK(edge_difference[4] < edge_difference[2]);
  CHECK(edge_difference[2] < edge_difference[0]);

  // Uniform h-refinement halves the edge-singular quadrature error and reduces
  // the node-singular error by roughly four. This verifies convergence without
  // claiming that depth five is accurate enough for a checked-in table.
  CHECK(node_difference[4] / node_difference[3] > 0.2);
  CHECK(node_difference[4] / node_difference[3] < 0.3);
  CHECK(edge_difference[4] / edge_difference[3] > 0.45);
  CHECK(edge_difference[4] / edge_difference[3] < 0.55);
}

TEST_CASE("Singular element volume and surface integrability", "[singularelements][Serial]")
{
  constexpr double nu = 0.5;
  constexpr int intervals = 4096;

  // The volume radial Jacobian contributes rho, so the integrand is
  // rho^(2*nu - 1), whose limit is finite for nu > 0.
  const double volume_2 = IntegrateRadialPowerOnLogGrid(2.0 * nu - 1.0, 1.0e-2, intervals);
  const double volume_6 = IntegrateRadialPowerOnLogGrid(2.0 * nu - 1.0, 1.0e-6, intervals);
  CHECK_THAT(volume_2, WithinAbs(1.0 - 1.0e-2, 5.0e-7));
  CHECK_THAT(volume_6, WithinAbs(1.0 - 1.0e-6, 2.0e-6));
  CHECK(volume_6 < 1.0);

  // A surface trace has no rho Jacobian. At nu = 1/2 its integrand is 1/rho,
  // and each fixed logarithmic refinement range adds a nonzero amount.
  const double surface_2 = IntegrateRadialPowerOnLogGrid(2.0 * nu - 2.0, 1.0e-2, intervals);
  const double surface_4 = IntegrateRadialPowerOnLogGrid(2.0 * nu - 2.0, 1.0e-4, intervals);
  const double surface_8 = IntegrateRadialPowerOnLogGrid(2.0 * nu - 2.0, 1.0e-8, intervals);
  CHECK_THAT(surface_2, WithinAbs(std::log(1.0e2), 1.0e-12));
  CHECK_THAT(surface_4, WithinAbs(std::log(1.0e4), 1.0e-12));
  CHECK_THAT(surface_8, WithinAbs(std::log(1.0e8), 1.0e-11));
  CHECK_THAT(surface_4 - surface_2, WithinAbs(surface_2, 1.0e-12));
  CHECK_THAT(surface_8 - surface_4, WithinAbs(surface_4, 1.0e-11));
}

TEST_CASE("Singular triangle lowest-order basis identities",
          "[singularelements][triangle][Serial]")
{
  const Vector2 point{0.27, 0.31};
  const auto lambda = fem::singular::ReferenceTriangleBarycentricPoint(point);
  const auto &grad_lambda = fem::singular::ReferenceTriangleBarycentricGradients();

  for (double nu : {0.5, 2.0 / 3.0})
  {
    CAPTURE(nu);
    for (int j : {1, 2})
    {
      const auto gradient =
          fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, j, nu);
      const auto numerical_gradient = NumericalGradient(
          point,
          [j, nu](const Vector2 &x)
          {
            return fem::singular::EvaluateTriangleNodeGradientPotential(
                fem::singular::ReferenceTriangleBarycentricPoint(x), 0, j, nu);
          });
      CheckVector(gradient.value, numerical_gradient, 2.0e-9);
      CHECK(gradient.curl == 0.0);

      const double numerical_curl =
          NumericalCurl(point,
                        [j, nu, &grad_lambda](const Vector2 &x)
                        {
                          return fem::singular::EvaluateTriangleNodeGradient(
                                     fem::singular::ReferenceTriangleBarycentricPoint(x),
                                     grad_lambda, 0, j, nu)
                              .value;
                        });
      CHECK_THAT(numerical_curl, WithinAbs(0.0, 2.0e-9));
    }

    const auto rotational =
        fem::singular::EvaluateTriangleNodeRotational(lambda, grad_lambda, 0, 1, 2, nu);
    const double rho = fem::singular::TriangleNodeRadialCoordinate(lambda, 0);
    CHECK_THAT(rotational.curl, WithinAbs(2.0 - (2.0 + nu) * std::pow(rho, nu), 2.0e-14));
    const double numerical_curl =
        NumericalCurl(point,
                      [nu, &grad_lambda](const Vector2 &x)
                      {
                        return fem::singular::EvaluateTriangleNodeRotational(
                                   fem::singular::ReferenceTriangleBarycentricPoint(x),
                                   grad_lambda, 0, 1, 2, nu)
                            .value;
                      });
    CHECK_THAT(rotational.curl, WithinAbs(numerical_curl, 2.0e-9));
  }
}

TEST_CASE("Singular triangle bases have the published edge traces",
          "[singularelements][triangle][Serial]")
{
  constexpr double nu = 0.5;
  constexpr std::array<Vector2, 3> vertices{Vector2{0.0, 0.0}, Vector2{1.0, 0.0},
                                            Vector2{0.0, 1.0}};
  const auto &grad_lambda = fem::singular::ReferenceTriangleBarycentricGradients();

  for (int opposite = 0; opposite < 3; opposite++)
  {
    std::array<int, 2> edge_vertices;
    int next = 0;
    for (int vertex = 0; vertex < 3; vertex++)
    {
      if (vertex != opposite)
      {
        edge_vertices[next++] = vertex;
      }
    }
    TriangleBarycentricPoint lambda{0.0, 0.0, 0.0};
    lambda[edge_vertices[0]] = 0.37;
    lambda[edge_vertices[1]] = 0.63;
    const Vector2 tangent{vertices[edge_vertices[1]][0] - vertices[edge_vertices[0]][0],
                          vertices[edge_vertices[1]][1] - vertices[edge_vertices[0]][1]};

    for (int j : {1, 2})
    {
      const auto gradient =
          fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, j, nu);
      const int associated_opposite = j == 1 ? 2 : 1;
      if (opposite == associated_opposite)
      {
        CHECK(std::abs(Dot(gradient.value, tangent)) > 1.0e-3);
      }
      else
      {
        CHECK_THAT(Dot(gradient.value, tangent), WithinAbs(0.0, 2.0e-14));
      }
    }

    const auto rotational =
        fem::singular::EvaluateTriangleNodeRotational(lambda, grad_lambda, 0, 1, 2, nu);
    CHECK_THAT(Dot(rotational.value, tangent), WithinAbs(0.0, 2.0e-14));
  }
}

TEST_CASE("Singular triangle basis is covariant under vertex relabeling",
          "[singularelements][triangle][Serial]")
{
  constexpr double nu = 2.0 / 3.0;
  const TriangleBarycentricPoint lambda{0.42, 0.23, 0.35};
  const TriangleBarycentricGradients grad_lambda{Vector2{-1.7, -0.4}, Vector2{0.6, -0.9},
                                                 Vector2{1.1, 1.3}};
  const auto gradient =
      fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, nu);
  const auto rotational =
      fem::singular::EvaluateTriangleNodeRotational(lambda, grad_lambda, 0, 1, 2, nu);

  std::array<int, 3> permutation{0, 1, 2};
  do
  {
    TriangleBarycentricPoint permuted_lambda;
    TriangleBarycentricGradients permuted_gradients;
    for (int i = 0; i < 3; i++)
    {
      permuted_lambda[permutation[i]] = lambda[i];
      permuted_gradients[permutation[i]] = grad_lambda[i];
    }
    const auto permuted_gradient = fem::singular::EvaluateTriangleNodeGradient(
        permuted_lambda, permuted_gradients, permutation[0], permutation[1], nu);
    const auto permuted_rotational = fem::singular::EvaluateTriangleNodeRotational(
        permuted_lambda, permuted_gradients, permutation[0], permutation[1], permutation[2],
        nu);
    CheckVector(permuted_gradient.value, gradient.value, 2.0e-14);
    CHECK(permuted_gradient.curl == gradient.curl);
    CheckVector(permuted_rotational.value, rotational.value, 2.0e-14);
    CHECK_THAT(permuted_rotational.curl, WithinAbs(rotational.curl, 2.0e-14));
  } while (std::next_permutation(permutation.begin(), permutation.end()));
}

TEST_CASE("Singular triangle gradient has the Meixner slope",
          "[singularelements][triangle][Serial]")
{
  const auto &grad_lambda = fem::singular::ReferenceTriangleBarycentricGradients();
  constexpr double angular_coordinate = 0.37;
  constexpr double rho_1 = 0x1p-48;
  constexpr double rho_2 = 0x1p-36;

  for (double nu : {0.5, 2.0 / 3.0})
  {
    const auto magnitude = [nu, &grad_lambda](double rho)
    {
      const TriangleBarycentricPoint lambda{1.0 - rho, rho * (1.0 - angular_coordinate),
                                            rho * angular_coordinate};
      return Norm(
          fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, nu).value);
    };
    const double slope = LogSlope(rho_1, magnitude(rho_1), rho_2, magnitude(rho_2));
    CAPTURE(nu, slope);
    CHECK_THAT(slope, WithinAbs(nu - 1.0, 1.0e-4));
  }
}

TEST_CASE("Singular triangle node Duffy quadrature", "[singularelements][triangle][Serial]")
{
  constexpr int order = 31;
  constexpr double radial_power = 6.0;

  double weight_sum = 0.0;
  std::size_t point_count = 0;
  fem::singular::ForEachReferenceTriangleNodeDuffyQuadraturePoint(
      order, 0, radial_power,
      [&weight_sum, &point_count](const TriangleBarycentricPoint &lambda, double weight)
      {
        for (double value : lambda)
        {
          CHECK(value > 0.0);
          CHECK(value <= 1.0);
        }
        CHECK(fem::singular::TriangleNodeRadialCoordinate(lambda, 0) > 0.0);
        CHECK(fem::singular::TriangleNodeRadialCoordinate(lambda, 0) < 1.0);
        CHECK(weight > 0.0);
        weight_sum += weight;
        point_count++;
      });
  CHECK(point_count > 0);
  CHECK(point_count ==
        fem::singular::ReferenceTriangleNodeDuffyQuadraturePointCount(order));
  CHECK_THAT(weight_sum, WithinAbs(0.5, 2.0e-15));

  const auto factorial = [](int n) { return std::tgamma(static_cast<double>(n) + 1.0); };
  for (int a = 0; a <= 4; a++)
  {
    for (int b = 0; b <= 4 - a; b++)
    {
      for (int c = 0; c <= 4 - a - b; c++)
      {
        const double value = fem::singular::IntegrateReferenceTriangleNodeDuffy(
            order, 0, radial_power,
            [a, b, c](const TriangleBarycentricPoint &lambda)
            {
              return std::pow(lambda[0], a) * std::pow(lambda[1], b) *
                     std::pow(lambda[2], c);
            });
        const double exact =
            factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 2);
        CAPTURE(a, b, c, value, exact);
        CHECK_THAT(value, WithinAbs(exact, 3.0e-15));
      }
    }
  }

  for (double exponent : {-1.0, -2.0 / 3.0, 0.0, 2.0})
  {
    const double value = fem::singular::IntegrateReferenceTriangleNodeDuffy(
        order, 0, radial_power,
        [exponent](const TriangleBarycentricPoint &lambda)
        {
          return std::pow(fem::singular::TriangleNodeRadialCoordinate(lambda, 0), exponent);
        });
    CAPTURE(exponent, value);
    CHECK_THAT(value, WithinAbs(1.0 / (exponent + 2.0), 3.0e-14));
  }

  constexpr double transmission_nu = 0.5255535;
  const double transmission_exponent = 2.0 * transmission_nu - 2.0;
  const double transmission_value = fem::singular::IntegrateReferenceTriangleNodeDuffy(
      95, 0, fem::singular::TriangleDuffyRadialPower,
      [transmission_exponent](const TriangleBarycentricPoint &lambda)
      {
        return std::pow(fem::singular::TriangleNodeRadialCoordinate(lambda, 0),
                        transmission_exponent);
      });
  CHECK_THAT(transmission_value, WithinAbs(1.0 / (transmission_exponent + 2.0), 2.0e-13));

  const auto &grad_lambda = fem::singular::ReferenceTriangleBarycentricGradients();
  const auto gradient_mass = [&grad_lambda](int quadrature_order, double nu)
  {
    return fem::singular::IntegrateReferenceTriangleNodeDuffy(
        quadrature_order, 0, fem::singular::TriangleDuffyRadialPower,
        [&grad_lambda, nu](const TriangleBarycentricPoint &lambda)
        {
          const auto basis =
              fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, nu);
          return Dot(basis.value, basis.value);
        });
  };
  for (double nu : {0.5, 2.0 / 3.0})
  {
    const double mass_23 = gradient_mass(23, nu);
    const double mass_31 = gradient_mass(31, nu);
    CAPTURE(nu, mass_23, mass_31);
    CHECK(std::isfinite(mass_31));
    CHECK(mass_31 > 0.0);
    CHECK_THAT(mass_23, WithinAbs(mass_31, 2.0e-13));
  }
}

TEST_CASE("Singular triangle rejects invalid inputs",
          "[singularelements][triangle][Serial]")
{
  const TriangleBarycentricPoint lambda{0.4, 0.3, 0.3};
  const auto &grad_lambda = fem::singular::ReferenceTriangleBarycentricGradients();
  CHECK_THROWS_AS(
      fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 0, 0.5),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, 0, 1, 0.0),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::EvaluateTriangleNodeGradient({1.0, 0.0, 0.0}, grad_lambda, 0, 1, 0.5),
      std::domain_error);
  CHECK_THROWS_AS(fem::singular::ForEachReferenceTriangleNodeDuffyQuadraturePoint(
                      0, 0, 6.0, [](const auto &, double) {}),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::IntegrateReferenceTriangleNodeDuffy(
                      12, 3, 6.0, [](const auto &) { return 1.0; }),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::ReferenceTriangleNodeDuffyQuadraturePointCount(0),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::IntegrateReferenceTriangleNodeDuffy(
          12, 0, 6.0, [](const auto &) { return std::numeric_limits<double>::infinity(); }),
      std::domain_error);
}

TEST_CASE("Weighted segment quadrature integrates Jacobi moments",
          "[singularelements][quadrature][Serial]")
{
  for (const auto [alpha, beta] : std::array<std::pair<double, double>, 4>{
           std::pair{0.0, 0.0}, std::pair{-0.75, 0.0}, std::pair{0.25, 1.5},
           std::pair{2.0 * 0.525553491856 - 2.0, 2.0 * 0.525553491856 - 2.0}})
  {
    const auto rule = fem::singular::BuildWeightedSegmentQuadrature(8, alpha, beta);
    REQUIRE(rule.size() == 8);
    for (int degree = 0; degree <= 15; degree++)
    {
      long double value = 0.0L;
      for (const auto &point : rule)
      {
        value += point.weight * std::pow(point.coordinate, degree);
      }
      const double exact =
          std::exp(std::lgamma(alpha + degree + 1.0) + std::lgamma(beta + 1.0) -
                   std::lgamma(alpha + beta + degree + 2.0));
      CAPTURE(alpha, beta, degree, value, exact);
      CHECK_THAT(static_cast<double>(value), WithinRel(exact, 2.0e-12));
    }
  }

  CHECK_THROWS_AS(fem::singular::BuildWeightedSegmentQuadrature(0, 0.0, 0.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::BuildWeightedSegmentQuadrature(4, -1.0, 0.0),
                  std::invalid_argument);
  CHECK_THROWS_AS(fem::singular::BuildWeightedSegmentQuadrature(
                      4, std::numeric_limits<double>::quiet_NaN(), 0.0),
                  std::invalid_argument);
}

namespace
{

// Phase 7 support. Interior points of the tetrahedral face opposite local vertex
// "opposite", deliberately nonsymmetric and away from both interpolation nodes and the
// singular loci, so a zero result cannot come from landing on a special point.
std::vector<BarycentricPoint> AsymmetricFacePoints(int opposite)
{
  const std::array<std::array<double, 3>, 4> interior{{{0.5137, 0.2791, 0.2072},
                                                       {0.1723, 0.6031, 0.2246},
                                                       {0.2611, 0.1907, 0.5482},
                                                       {0.3931, 0.3557, 0.2512}}};
  std::vector<BarycentricPoint> points;
  points.reserve(interior.size());
  for (const auto &weights : interior)
  {
    BarycentricPoint lambda{0.0, 0.0, 0.0, 0.0};
    int next = 0;
    for (int node = 0; node < 4; node++)
    {
      if (node != opposite)
      {
        lambda[node] = weights[next++];
      }
    }
    points.push_back(lambda);
  }
  return points;
}

// Largest H1 potential magnitude and largest ND tangential-trace magnitude of one basis
// over the sampled interior points of one face.
std::pair<double, double> MaximumFaceTrace(const fem::singular::HigherOrderBasis &basis,
                                           int opposite)
{
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  const bool gradient_family =
      basis.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT ||
      basis.family == fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT;
  double potential = 0.0, tangential = 0.0;
  for (const auto &lambda : AsymmetricFacePoints(opposite))
  {
    if (gradient_family)
    {
      potential = std::max(
          potential,
          std::abs(fem::singular::EvaluateHigherOrderGradientPotential(lambda, basis)));
    }
    const auto value = fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, basis);
    tangential =
        std::max(tangential, Norm(TangentialPart(value.value, grad_lambda[opposite])));
  }
  return {potential, tangential};
}

// Every retained basis of a node feature at canonical_nodes, and of an edge feature on the
// first two canonical nodes, at one singular order and exponent.
std::vector<fem::singular::HigherOrderBasis>
NodeFeatureBases(const std::array<int, 4> &canonical_nodes, int order, double nu)
{
  auto bases =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, order, nu);
  if (order >= 2)
  {
    const auto rotational =
        fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, order, nu);
    bases.insert(bases.end(), rotational.begin(), rotational.end());
  }
  return bases;
}

std::vector<fem::singular::HigherOrderBasis>
EdgeFeatureBases(const std::array<int, 4> &canonical_nodes, int order, double nu)
{
  std::vector<fem::singular::HigherOrderBasis> bases;
  if (order >= 2)
  {
    bases =
        fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, order, nu);
    const auto rotational =
        fem::singular::EnumerateHigherOrderEdgeRotationalBases(canonical_nodes, order, nu);
    bases.insert(bases.end(), rotational.begin(), rotational.end());
  }
  return bases;
}

}  // namespace

TEST_CASE("Singular face-trace classification matches direct basis evaluation",
          "[singularelements][Serial]")
{
  // Phase 7 of the singular AMR plan. The basis algebra proves which faces carry zero
  // trace; this test ensures the classifier used by refinement agrees with that proof for
  // every retained descriptor. Both directions are asserted: a face the classifier calls
  // inactive must be numerically zero, and a face it calls active must have at least one
  // demonstrably nonzero trace somewhere among its bases.
  const int order = GENERATE(1, 2, 3, 4);
  // Thin sheet, finite-thickness wedge, and a generic material-derived exponent.
  const double nu = GENERATE(0.5, 2.0 / 3.0, 0.731);
  CAPTURE(order, nu);

  constexpr double zero_tolerance = 1.0e-12;
  // All 24 vertex permutations, so the classification cannot depend on node ordering.
  std::array<int, 4> permutation{0, 1, 2, 3};
  std::sort(permutation.begin(), permutation.end());
  int permutations = 0;
  do
  {
    CAPTURE(permutation);
    permutations++;

    // Node feature: the singular node is canonical_nodes[0]. Expect exactly the face
    // opposite that node to be inactive.
    {
      const auto bases = NodeFeatureBases(permutation, order, nu);
      REQUIRE_FALSE(bases.empty());
      const int singular = permutation[0];
      for (int opposite = 0; opposite < 4; opposite++)
      {
        double potential = 0.0, tangential = 0.0;
        for (const auto &basis : bases)
        {
          const auto [p, t] = MaximumFaceTrace(basis, opposite);
          potential = std::max(potential, p);
          tangential = std::max(tangential, t);
        }
        CAPTURE(opposite, singular, potential, tangential);
        if (opposite == singular)
        {
          // Face opposite the singular node: both traces must vanish.
          CHECK(potential <= zero_tolerance);
          CHECK(tangential <= zero_tolerance);
        }
        else
        {
          // Face containing the singular node: something must be nonzero, otherwise the
          // classifier would be needlessly conservative and the test vacuous.
          CHECK(std::max(potential, tangential) > zero_tolerance);
        }
      }
    }

    // Edge feature on canonical nodes 0 and 1. Expect the two faces which do not contain
    // both endpoints to be inactive.
    if (order >= 2)
    {
      const auto bases = EdgeFeatureBases(permutation, order, nu);
      REQUIRE_FALSE(bases.empty());
      const int first = permutation[0], second = permutation[1];
      for (int opposite = 0; opposite < 4; opposite++)
      {
        const bool contains_edge = (opposite != first) && (opposite != second);
        double potential = 0.0, tangential = 0.0;
        for (const auto &basis : bases)
        {
          const auto [p, t] = MaximumFaceTrace(basis, opposite);
          potential = std::max(potential, p);
          tangential = std::max(tangential, t);
        }
        CAPTURE(opposite, first, second, contains_edge, potential, tangential);
        if (!contains_edge)
        {
          CHECK(potential <= zero_tolerance);
          CHECK(tangential <= zero_tolerance);
        }
        else
        {
          CHECK(std::max(potential, tangential) > zero_tolerance);
        }
      }
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  CHECK(permutations == 24);
}

TEST_CASE("Edge-rotational singular bases are H(curl) bubbles on every face",
          "[singularelements][Serial]")
{
  // The edge-rotational family must never activate a face on its own: its tangential trace
  // vanishes on all four faces. This is what lets an edge feature exempt two faces even
  // though its rotational bases are supported throughout the element.
  const int order = GENERATE(2, 3, 4);
  const double nu = GENERATE(0.5, 2.0 / 3.0, 0.731);
  CAPTURE(order, nu);
  constexpr std::array<int, 4> Nodes{0, 1, 2, 3};
  const auto bases =
      fem::singular::EnumerateHigherOrderEdgeRotationalBases(Nodes, order, nu);
  REQUIRE_FALSE(bases.empty());
  for (const auto &basis : bases)
  {
    for (int opposite = 0; opposite < 4; opposite++)
    {
      const auto [potential, tangential] = MaximumFaceTrace(basis, opposite);
      MFEM_CONTRACT_VAR(potential);
      CAPTURE(opposite, tangential);
      CHECK(tangential <= 1.0e-12);
    }
  }
}

}  // namespace palace
