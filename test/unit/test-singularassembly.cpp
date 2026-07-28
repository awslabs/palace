// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/singularassembly.hpp"
#include "fem/singularsystem.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/iterative.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

double Dot(const fem::singular::Vector3 &a, const fem::singular::Vector3 &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Dot(const fem::singular::Vector2 &a, const fem::singular::Vector2 &b)
{
  return a[0] * b[0] + a[1] * b[1];
}

fem::singular::Vector3 Cross(const fem::singular::Vector3 &a,
                             const fem::singular::Vector3 &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

fem::singular::BarycentricGradients
AffineBarycentricGradients(const fem::singular::Vector3 &a, const fem::singular::Vector3 &b,
                           const fem::singular::Vector3 &c, double &jacobian_determinant)
{
  jacobian_determinant = Dot(a, Cross(b, c));
  REQUIRE(jacobian_determinant > 0.0);
  fem::singular::BarycentricGradients gradients;
  const auto cross_bc = Cross(b, c);
  const auto cross_ca = Cross(c, a);
  const auto cross_ab = Cross(a, b);
  for (int d = 0; d < 3; d++)
  {
    gradients[1][d] = cross_bc[d] / jacobian_determinant;
    gradients[2][d] = cross_ca[d] / jacobian_determinant;
    gradients[3][d] = cross_ab[d] / jacobian_determinant;
    gradients[0][d] = -gradients[1][d] - gradients[2][d] - gradients[3][d];
  }
  return gradients;
}

fem::singular::BarycentricGradients
PhysicalBarycentricGradients(mfem::ElementTransformation &transformation,
                             double &jacobian_determinant)
{
  mfem::IntegrationPoint center;
  center.Set3(0.25, 0.25, 0.25);
  transformation.SetIntPoint(&center);
  jacobian_determinant = transformation.Weight();
  const auto &inverse_jacobian = transformation.InverseJacobian();
  REQUIRE(inverse_jacobian.Height() == 3);
  REQUIRE(inverse_jacobian.Width() == 3);
  fem::singular::BarycentricGradients gradients{};
  for (int i = 1; i < 4; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      gradients[i][d] = inverse_jacobian(i - 1, d);
      gradients[0][d] -= gradients[i][d];
    }
  }
  return gradients;
}

mfem::Mesh AffineTetrahedronMesh(const fem::singular::Vector3 &origin,
                                 const fem::singular::Vector3 &a,
                                 const fem::singular::Vector3 &b,
                                 const fem::singular::Vector3 &c)
{
  mfem::Mesh mesh(3, 4, 1, 0, 3);
  mesh.AddVertex(origin.data());
  std::array<fem::singular::Vector3, 3> vertices{a, b, c};
  for (auto &vertex : vertices)
  {
    for (int d = 0; d < 3; d++)
    {
      vertex[d] += origin[d];
    }
    mesh.AddVertex(vertex.data());
  }
  mesh.AddTet(0, 1, 2, 3, 1);
  mesh.FinalizeTopology();
  mesh.Finalize(false, false);
  return mesh;
}

mfem::Mesh AffineTriangleMesh(const fem::singular::Vector2 &origin,
                              const fem::singular::Vector2 &a,
                              const fem::singular::Vector2 &b)
{
  mfem::Mesh mesh(2, 3, 1, 0, 2);
  mesh.AddVertex(origin.data());
  const fem::singular::Vector2 vertex_a{origin[0] + a[0], origin[1] + a[1]};
  const fem::singular::Vector2 vertex_b{origin[0] + b[0], origin[1] + b[1]};
  mesh.AddVertex(vertex_a.data());
  mesh.AddVertex(vertex_b.data());
  mesh.AddTriangle(0, 1, 2, 1);
  mesh.FinalizeTopology();
  mesh.Finalize(false, false);
  return mesh;
}

void SetQuadraticGeometry(mfem::Mesh &mesh, bool curved)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  if (!curved)
  {
    return;
  }

  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        value = x;
        if (dimension == 2)
        {
          value[1] += 0.2 * x[0] * x[1];
        }
        else
        {
          value[2] += 0.15 * x[0] * x[1];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

mfem::Mesh InternalLineTipTriangleMesh()
{
  mfem::Mesh mesh(2, 9, 8, 1, 2);
  mesh.AddVertex(-1.0, -1.0);
  mesh.AddVertex(0.0, -1.0);
  mesh.AddVertex(1.0, -1.0);
  mesh.AddVertex(-1.0, 0.0);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(-1.0, 1.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);
  mesh.AddTriangle(0, 1, 4, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(3, 4, 7, 1);
  mesh.AddTriangle(3, 7, 6, 1);
  mesh.AddTriangle(1, 2, 5, 1);
  mesh.AddTriangle(1, 5, 4, 1);
  mesh.AddTriangle(4, 5, 8, 1);
  mesh.AddTriangle(4, 8, 7, 1);
  mesh.AddBdrSegment(3, 4, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

std::array<double, 2>
DirectPhysicalPairIntegral(const fem::singular::HigherOrderBasis &row,
                           const fem::singular::HigherOrderBasis &column,
                           const fem::singular::BarycentricGradients &grad_lambda,
                           double jacobian_determinant)
{
  long double mass = 0.0L, curl_curl = 0.0L;
  fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
      45, 0, 6.0,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        const auto row_value =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, row);
        const auto column_value =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, column);
        mass += static_cast<long double>(jacobian_determinant) *
                static_cast<long double>(weight) *
                static_cast<long double>(Dot(row_value.value, column_value.value));
        curl_curl += static_cast<long double>(jacobian_determinant) *
                     static_cast<long double>(weight) *
                     static_cast<long double>(Dot(row_value.curl, column_value.curl));
      });
  return {static_cast<double>(mass), static_cast<double>(curl_curl)};
}

double DirectH1StandardEnrichmentIntegral(
    const mfem::FiniteElement &h1_fe, mfem::ElementTransformation &transformation,
    int standard, const fem::singular::HigherOrderBasis &enrichment,
    const fem::singular::BarycentricGradients &grad_lambda, double jacobian_determinant)
{
  long double value = 0.0L;
  mfem::DenseMatrix standard_gradient(h1_fe.GetDof(), 3);
  mfem::IntegrationPoint point;
  fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
      35, 0, 6.0,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        point.Set3(lambda[1], lambda[2], lambda[3]);
        transformation.SetIntPoint(&point);
        h1_fe.CalcPhysDShape(transformation, standard_gradient);
        const auto singular =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, enrichment);
        const fem::singular::Vector3 standard_value{standard_gradient(standard, 0),
                                                    standard_gradient(standard, 1),
                                                    standard_gradient(standard, 2)};
        value += static_cast<long double>(jacobian_determinant) *
                 static_cast<long double>(weight) *
                 static_cast<long double>(Dot(standard_value, singular.value));
      });
  return static_cast<double>(value);
}

void ForEachAlignedDuffyQuadraturePoint(
    const fem::singular::HigherOrderBasis &basis, int order,
    const fem::singular::QuadraturePointVisitor &visitor)
{
  switch (basis.family)
  {
    case fem::singular::HigherOrderBasisFamily::NODE_GRADIENT:
      fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
          order, basis.nodes[0], fem::singular::H1DuffyRadialPower, visitor);
      return;
    case fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT:
      fem::singular::ForEachReferenceTetrahedronEdgeDuffyQuadraturePoint(
          order, basis.nodes[0], basis.nodes[1], fem::singular::H1DuffyRadialPower,
          visitor);
      return;
    case fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL:
    case fem::singular::HigherOrderBasisFamily::EDGE_ROTATIONAL:
      FAIL("H1 Duffy quadrature requires a gradient-family basis");
  }
}

double
DirectH1EnrichmentPairDuffyIntegral(const fem::singular::HigherOrderBasis &row,
                                    const fem::singular::HigherOrderBasis &column,
                                    const fem::singular::BarycentricGradients &grad_lambda,
                                    double jacobian_determinant, int order)
{
  long double value = 0.0L;
  ForEachAlignedDuffyQuadraturePoint(
      row, order,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        const auto row_value =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, row);
        const auto column_value =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, column);
        value += static_cast<long double>(jacobian_determinant) *
                 static_cast<long double>(weight) *
                 static_cast<long double>(Dot(row_value.value, column_value.value));
      });
  return static_cast<double>(value);
}

double DirectH1StandardEnrichmentDuffyIntegral(
    const mfem::FiniteElement &h1_fe, mfem::ElementTransformation &transformation,
    int standard, const fem::singular::HigherOrderBasis &enrichment,
    const fem::singular::BarycentricGradients &grad_lambda, double jacobian_determinant,
    int order)
{
  long double value = 0.0L;
  mfem::DenseMatrix standard_gradient(h1_fe.GetDof(), 3);
  mfem::IntegrationPoint point;
  ForEachAlignedDuffyQuadraturePoint(
      enrichment, order,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        point.Set3(lambda[1], lambda[2], lambda[3]);
        transformation.SetIntPoint(&point);
        h1_fe.CalcPhysDShape(transformation, standard_gradient);
        const auto singular =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, enrichment);
        const fem::singular::Vector3 standard_value{standard_gradient(standard, 0),
                                                    standard_gradient(standard, 1),
                                                    standard_gradient(standard, 2)};
        value += static_cast<long double>(jacobian_determinant) *
                 static_cast<long double>(weight) *
                 static_cast<long double>(Dot(standard_value, singular.value));
      });
  return static_cast<double>(value);
}

std::array<double, 2> DirectNDStandardEnrichmentIntegral(
    const mfem::FiniteElement &nd_fe, mfem::ElementTransformation &transformation,
    int standard, const fem::singular::HigherOrderBasis &enrichment,
    const fem::singular::BarycentricGradients &grad_lambda, double jacobian_determinant)
{
  std::array<long double, 2> value{};
  mfem::DenseMatrix standard_shape(nd_fe.GetDof(), 3);
  mfem::DenseMatrix standard_curl(nd_fe.GetDof(), 3);
  mfem::IntegrationPoint point;
  fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
      35, 0, 6.0,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        point.Set3(lambda[1], lambda[2], lambda[3]);
        transformation.SetIntPoint(&point);
        nd_fe.CalcPhysVShape(transformation, standard_shape);
        nd_fe.CalcPhysCurlShape(transformation, standard_curl);
        const auto singular =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, enrichment);
        const fem::singular::Vector3 standard_value{standard_shape(standard, 0),
                                                    standard_shape(standard, 1),
                                                    standard_shape(standard, 2)};
        const fem::singular::Vector3 standard_curl_value{standard_curl(standard, 0),
                                                         standard_curl(standard, 1),
                                                         standard_curl(standard, 2)};
        value[0] += static_cast<long double>(jacobian_determinant) *
                    static_cast<long double>(weight) *
                    static_cast<long double>(Dot(standard_value, singular.value));
        value[1] += static_cast<long double>(jacobian_determinant) *
                    static_cast<long double>(weight) *
                    static_cast<long double>(Dot(standard_curl_value, singular.curl));
      });
  return {static_cast<double>(value[0]), static_cast<double>(value[1])};
}

Eigen::MatrixXd ToEigen(const mfem::DenseMatrix &matrix)
{
  Eigen::MatrixXd result(matrix.Height(), matrix.Width());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      result(i, j) = matrix(i, j);
    }
  }
  return result;
}

double FrobeniusNorm(const mfem::DenseMatrix &matrix)
{
  double norm_squared = 0.0;
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      norm_squared += matrix(i, j) * matrix(i, j);
    }
  }
  return std::sqrt(norm_squared);
}

void CheckSymmetric(const mfem::DenseMatrix &matrix)
{
  REQUIRE(matrix.Height() == matrix.Width());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      CHECK(matrix(i, j) == matrix(j, i));
    }
  }
}

void CheckExactTranspose(const mfem::DenseMatrix &matrix,
                         const mfem::DenseMatrix &transpose)
{
  REQUIRE(matrix.Height() == transpose.Width());
  REQUIRE(matrix.Width() == transpose.Height());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      CHECK(matrix(i, j) == transpose(j, i));
    }
  }
}

void CheckMatricesEqual(const mfem::DenseMatrix &left, const mfem::DenseMatrix &right)
{
  REQUIRE(left.Height() == right.Height());
  REQUIRE(left.Width() == right.Width());
  for (int i = 0; i < left.Height(); i++)
  {
    for (int j = 0; j < left.Width(); j++)
    {
      CHECK(left(i, j) == right(i, j));
    }
  }
}

void CheckScaledMatrix(const mfem::DenseMatrix &scaled, const mfem::DenseMatrix &reference,
                       double factor)
{
  REQUIRE(scaled.Height() == reference.Height());
  REQUIRE(scaled.Width() == reference.Width());
  for (int i = 0; i < scaled.Height(); i++)
  {
    for (int j = 0; j < scaled.Width(); j++)
    {
      CHECK(scaled(i, j) == factor * reference(i, j));
    }
  }
}

void CheckAffineScaling(const mfem::DenseMatrix &scaled,
                        const mfem::DenseMatrix &scaled_error,
                        const mfem::DenseMatrix &reference,
                        const mfem::DenseMatrix &reference_error, double factor)
{
  REQUIRE(scaled.Height() == reference.Height());
  REQUIRE(scaled.Width() == reference.Width());
  REQUIRE(scaled_error.Height() == scaled.Height());
  REQUIRE(scaled_error.Width() == scaled.Width());
  REQUIRE(reference_error.Height() == reference.Height());
  REQUIRE(reference_error.Width() == reference.Width());
  for (int i = 0; i < scaled.Height(); i++)
  {
    for (int j = 0; j < scaled.Width(); j++)
    {
      const double expected = factor * reference(i, j);
      const double roundoff = 256.0 * std::numeric_limits<double>::epsilon() *
                              std::max({1.0, std::abs(scaled(i, j)), std::abs(expected)});
      CAPTURE(i, j, scaled(i, j), expected, scaled_error(i, j), reference_error(i, j),
              factor);
      CHECK(std::abs(scaled(i, j) - expected) <=
            scaled_error(i, j) + factor * reference_error(i, j) + roundoff);
    }
  }
}

void CheckExactSparseTranspose(const mfem::SparseMatrix &matrix,
                               const mfem::SparseMatrix &transpose)
{
  REQUIRE(matrix.Height() == transpose.Width());
  REQUIRE(matrix.Width() == transpose.Height());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      CHECK(matrix(i, j) == transpose(j, i));
    }
  }
}

void CheckSparseNonnegative(const mfem::SparseMatrix &matrix)
{
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      CHECK(matrix(i, j) >= 0.0);
    }
  }
}

double SparseBilinear(const mfem::Vector &left, const mfem::SparseMatrix &matrix,
                      const mfem::Vector &right)
{
  REQUIRE(left.Size() == matrix.Height());
  REQUIRE(right.Size() == matrix.Width());
  mfem::Vector result(matrix.Height());
  matrix.Mult(right, result);
  return left * result;
}

mfem::Vector AbsoluteValue(const mfem::Vector &vector)
{
  mfem::Vector result(vector.Size());
  for (int i = 0; i < vector.Size(); i++)
  {
    result[i] = std::abs(vector[i]);
  }
  return result;
}

void CheckPositiveSemidefinite(const mfem::DenseMatrix &matrix,
                               const mfem::DenseMatrix &error)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(ToEigen(matrix));
  REQUIRE(solver.info() == Eigen::Success);
  const double error_bound = FrobeniusNorm(error);
  CHECK(solver.eigenvalues().minCoeff() >= -2.0 * error_bound - 1.0e-11);
}

mfem::Mesh SharedFaceTetrahedronMesh(bool permute_vertices)
{
  mfem::Mesh mesh(3, 5, 2, 0, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  mesh.AddVertex(0.0, 0.0, -1.0);
  std::array<std::array<int, 4>, 2> elements{std::array<int, 4>{0, 1, 2, 3},
                                             std::array<int, 4>{0, 2, 1, 4}};
  if (permute_vertices)
  {
    constexpr std::array<int, 4> permutation{1, 2, 0, 3};
    for (auto &element : elements)
    {
      const auto original = element;
      for (int i = 0; i < 4; i++)
      {
        element[i] = original[permutation[i]];
      }
    }
  }
  for (int i = 0; i < static_cast<int>(elements.size()); i++)
  {
    mesh.AddTet(elements[i].data(), i + 1);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(false, false);
  return mesh;
}

std::array<int, 4> CanonicalNodeTuple(const mfem::Element &element, int singular_vertex)
{
  int singular_local = -1;
  std::array<std::pair<int, int>, 3> remaining;
  int next = 0;
  for (int local = 0; local < 4; local++)
  {
    const int vertex = element.GetVertices()[local];
    if (vertex == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {vertex, local};
    }
  }
  REQUIRE(singular_local >= 0);
  REQUIRE(next == 3);
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second, remaining[2].second};
}

fem::singular::ElementDofMap SharedFaceSingularDofs(const mfem::Element &element)
{
  constexpr double nu = 0.5;
  const auto canonical_nodes = CanonicalNodeTuple(element, 0);
  const auto gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu);

  const fem::singular::HigherOrderBasis *shared_gradient = nullptr;
  for (const auto &basis : gradients)
  {
    if (element.GetVertices()[basis.nodes[1]] == 1)
    {
      REQUIRE(shared_gradient == nullptr);
      shared_gradient = &basis;
    }
  }
  const fem::singular::HigherOrderBasis *shared_rotation = nullptr;
  for (const auto &basis : rotations)
  {
    std::array<int, 2> component{element.GetVertices()[basis.nodes[1]],
                                 element.GetVertices()[basis.nodes[2]]};
    std::sort(component.begin(), component.end());
    if (component == std::array<int, 2>{1, 2})
    {
      REQUIRE(shared_rotation == nullptr);
      shared_rotation = &basis;
    }
  }
  REQUIRE(shared_gradient != nullptr);
  REQUIRE(shared_rotation != nullptr);

  fem::singular::ElementDofMap result;
  result.h1.push_back({0, *shared_gradient});
  result.nd.push_back({0, *shared_gradient});
  result.nd.push_back({1, *shared_rotation});
  return result;
}

std::array<int, 4> CanonicalStableNodeTuple(
    const mfem::Element &element,
    const std::vector<fem::singular::GlobalVertexId> &serial_vertex_ids,
    fem::singular::GlobalVertexId singular_vertex)
{
  int singular_local = -1;
  std::array<std::pair<fem::singular::GlobalVertexId, int>, 3> remaining;
  int next = 0;
  for (int local = 0; local < 4; local++)
  {
    const int mesh_vertex = element.GetVertices()[local];
    REQUIRE(mesh_vertex >= 0);
    REQUIRE(mesh_vertex < static_cast<int>(serial_vertex_ids.size()));
    const auto stable_vertex = serial_vertex_ids[mesh_vertex];
    if (stable_vertex == singular_vertex)
    {
      singular_local = local;
    }
    else
    {
      remaining[next++] = {stable_vertex, local};
    }
  }
  REQUIRE(singular_local >= 0);
  REQUIRE(next == 3);
  std::sort(remaining.begin(), remaining.end());
  return {singular_local, remaining[0].second, remaining[1].second, remaining[2].second};
}

fem::singular::ElementDofMap
SharedFaceSingularDofs(const mfem::Element &element,
                       const std::vector<fem::singular::GlobalVertexId> &serial_vertex_ids)
{
  constexpr double nu = 0.5;
  const auto canonical_nodes = CanonicalStableNodeTuple(element, serial_vertex_ids, 0);
  const auto gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu);
  const auto stable_vertex = [&](int local)
  {
    const int mesh_vertex = element.GetVertices()[local];
    REQUIRE(mesh_vertex >= 0);
    REQUIRE(mesh_vertex < static_cast<int>(serial_vertex_ids.size()));
    return serial_vertex_ids[mesh_vertex];
  };

  const fem::singular::HigherOrderBasis *shared_gradient = nullptr;
  for (const auto &basis : gradients)
  {
    if (stable_vertex(basis.nodes[1]) == 1)
    {
      REQUIRE(shared_gradient == nullptr);
      shared_gradient = &basis;
    }
  }
  const fem::singular::HigherOrderBasis *shared_rotation = nullptr;
  for (const auto &basis : rotations)
  {
    std::array<fem::singular::GlobalVertexId, 2> component{stable_vertex(basis.nodes[1]),
                                                           stable_vertex(basis.nodes[2])};
    std::sort(component.begin(), component.end());
    if (component == std::array<fem::singular::GlobalVertexId, 2>{1, 2})
    {
      REQUIRE(shared_rotation == nullptr);
      shared_rotation = &basis;
    }
  }
  REQUIRE(shared_gradient != nullptr);
  REQUIRE(shared_rotation != nullptr);

  fem::singular::ElementDofMap result;
  result.h1.push_back({0, *shared_gradient});
  result.nd.push_back({0, *shared_gradient});
  result.nd.push_back({1, *shared_rotation});
  return result;
}

fem::singular::EntityKey
Entity(std::initializer_list<fem::singular::GlobalVertexId> vertices)
{
  REQUIRE(vertices.size() > 0);
  REQUIRE(vertices.size() <= 4);
  fem::singular::EntityKey result;
  result.size = vertices.size();
  std::copy(vertices.begin(), vertices.end(), result.vertices.begin());
  std::sort(result.vertices.begin(), result.vertices.begin() + result.size);
  return result;
}

fem::singular::DofTopology
SharedFaceDofTopology(const mfem::Mesh &mesh,
                      const std::vector<fem::singular::GlobalVertexId> &serial_vertex_ids)
{
  fem::singular::DofKey gradient;
  gradient.family = fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
  gradient.order = 1;
  gradient.singular_entity = Entity({0});
  gradient.support_entity = Entity({0, 1});
  gradient.component_entity = Entity({0, 1});
  gradient.interpolation_weights = {1, 1, 0, 0};

  fem::singular::DofKey rotation;
  rotation.family = fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
  rotation.order = 1;
  rotation.singular_entity = Entity({0});
  rotation.support_entity = Entity({0, 1, 2});
  rotation.component_entity = Entity({1, 2});
  rotation.interpolation_weights = {1, 1, 1, 0};

  fem::singular::DofTopology result;
  result.h1_dofs = {gradient};
  result.nd_dofs = {gradient, rotation};
  result.h1_to_nd = {0};
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    result.elements.push_back(
        SharedFaceSingularDofs(*mesh.GetElement(element), serial_vertex_ids));
  }
  return result;
}

mfem::DenseMatrix GatherParallelMatrix(const mfem::HypreParMatrix &matrix)
{
  auto *parallel = static_cast<hypre_ParCSRMatrix *>(matrix);
  const HYPRE_BigInt global_rows = matrix.GetGlobalNumRows();
  const HYPRE_BigInt global_columns = matrix.GetGlobalNumCols();
  REQUIRE(global_rows >= 0);
  REQUIRE(global_rows <= std::numeric_limits<int>::max());
  REQUIRE(global_columns >= 0);
  REQUIRE(global_columns <= std::numeric_limits<int>::max());
  const int rows = static_cast<int>(global_rows);
  const int columns = static_cast<int>(global_columns);
  const HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parallel);
  const HYPRE_BigInt first_column = hypre_ParCSRMatrixFirstColDiag(parallel);
  const int local_rows = hypre_CSRMatrixNumRows(hypre_ParCSRMatrixDiag(parallel));
  const int local_columns = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixDiag(parallel));
  REQUIRE(matrix.Height() == local_rows);
  REQUIRE(matrix.Width() == local_columns);

  mfem::DenseMatrix result(rows, columns);
  result = 0.0;
  mfem::Vector x(local_columns), y(local_rows);
  std::vector<double> global_column(rows);
  for (int column = 0; column < columns; column++)
  {
    x = 0.0;
    if (first_column <= column && column < first_column + local_columns)
    {
      x[column - first_column] = 1.0;
    }
    matrix.Mult(x, y);
    std::fill(global_column.begin(), global_column.end(), 0.0);
    for (int row = 0; row < local_rows; row++)
    {
      REQUIRE(first_row + row >= 0);
      REQUIRE(first_row + row < rows);
      global_column[first_row + row] = y[row];
    }
    Mpi::GlobalSum(rows, global_column.data(), matrix.GetComm());
    for (int row = 0; row < rows; row++)
    {
      result(row, column) = global_column[row];
    }
  }
  return result;
}

mfem::DenseMatrix DenseMatrix(const mfem::SparseMatrix &matrix)
{
  mfem::DenseMatrix result(matrix.Height(), matrix.Width());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      result(i, j) = matrix(i, j);
    }
  }
  return result;
}

double ParallelBilinear(const mfem::Vector &left, const mfem::HypreParMatrix &matrix,
                        const mfem::Vector &right)
{
  REQUIRE(left.Size() == matrix.Height());
  REQUIRE(right.Size() == matrix.Width());
  mfem::Vector result(matrix.Height());
  matrix.Mult(right, result);
  double value = left * result;
  Mpi::GlobalSum(1, &value, matrix.GetComm());
  return value;
}

void CheckClose(double value, double reference)
{
  CHECK(std::abs(value - reference) <=
        2.0e-11 * std::max({1.0, std::abs(value), std::abs(reference)}));
}

void CheckClose(const mfem::DenseMatrix &matrix, const mfem::DenseMatrix &reference)
{
  REQUIRE(matrix.Height() == reference.Height());
  REQUIRE(matrix.Width() == reference.Width());
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      CAPTURE(i, j);
      CheckClose(matrix(i, j), reference(i, j));
    }
  }
}

double BilinearErrorBound(const mfem::Vector &standard, const mfem::DenseMatrix &error,
                          const mfem::Vector &enrichment)
{
  REQUIRE(standard.Size() == error.Height());
  REQUIRE(enrichment.Size() == error.Width());
  long double result = 0.0L;
  for (int i = 0; i < error.Height(); i++)
  {
    for (int j = 0; j < error.Width(); j++)
    {
      result += std::abs(static_cast<long double>(standard[i])) *
                static_cast<long double>(error(i, j)) *
                std::abs(static_cast<long double>(enrichment[j]));
    }
  }
  return static_cast<double>(result);
}

struct CouplingEnergies
{
  double h1;
  double nd_mass;
  double nd_curl_curl;
  double h1_error;
  double nd_mass_error;
  double nd_curl_curl_error;
};

CouplingEnergies AssembleSharedFaceCoupling(bool permute_vertices)
{
  auto mesh = SharedFaceTetrahedronMesh(permute_vertices);
  mfem::H1_FECollection h1_collection(2, 3);
  mfem::ND_FECollection nd_collection(2, 3);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);

  mfem::FunctionCoefficient h1_coefficient(
      [](const mfem::Vector &x) { return x[0] + 2.0 * x[1] - 0.5 * x[2] + x[0] * x[1]; });
  mfem::VectorFunctionCoefficient nd_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value[0] = x[1] + 2.0 * x[2];
        value[1] = -x[0] + 0.5 * x[2];
        value[2] = 3.0 * x[0] - x[1];
      });
  mfem::GridFunction h1_field(&h1_space);
  mfem::GridFunction nd_field(&nd_space);
  h1_field.ProjectCoefficient(h1_coefficient);
  nd_field.ProjectCoefficient(nd_coefficient);

  mfem::Vector h1_enrichment(1);
  h1_enrichment[0] = 0.37;
  mfem::Vector nd_enrichment(2);
  nd_enrichment[0] = 0.37;
  nd_enrichment[1] = -0.21;
  mfem::Array<int> h1_enrichment_dofs(1);
  h1_enrichment_dofs[0] = 0;
  mfem::Array<int> nd_enrichment_dofs(2);
  nd_enrichment_dofs[0] = 0;
  nd_enrichment_dofs[1] = 1;

  mfem::SparseMatrix h1_forward(h1_space.GetVSize(), 1);
  mfem::SparseMatrix h1_reverse(1, h1_space.GetVSize());
  mfem::SparseMatrix nd_mass_forward(nd_space.GetVSize(), 2);
  mfem::SparseMatrix nd_mass_reverse(2, nd_space.GetVSize());
  mfem::SparseMatrix nd_curl_curl_forward(nd_space.GetVSize(), 2);
  mfem::SparseMatrix nd_curl_curl_reverse(2, nd_space.GetVSize());
  CouplingEnergies local{};
  bool has_signed_nd_dof = false;
  bool has_nontrivial_face_transform = false;
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto element_dofs = SharedFaceSingularDofs(*mesh.GetElement(element));
    auto matrices = fem::singular::AssembleElementStandardEnrichmentMatrices(
        element_dofs, *h1_space.GetFE(element), *nd_space.GetFE(element),
        *mesh.GetElementTransformation(element), options);

    mfem::Array<int> h1_vdofs, nd_vdofs;
    mfem::DofTransformation h1_dof_transformation, nd_dof_transformation;
    h1_space.GetElementVDofs(element, h1_vdofs, h1_dof_transformation);
    nd_space.GetElementVDofs(element, nd_vdofs, nd_dof_transformation);
    has_signed_nd_dof = has_signed_nd_dof || std::any_of(nd_vdofs.begin(), nd_vdofs.end(),
                                                         [](int dof) { return dof < 0; });

    mfem::Vector h1_local, nd_local;
    h1_field.GetSubVector(h1_vdofs, h1_local);
    nd_field.GetSubVector(nd_vdofs, nd_local);
    h1_dof_transformation.InvTransformPrimal(h1_local);
    nd_dof_transformation.InvTransformPrimal(nd_local);
    mfem::Vector local_result(h1_local.Size());
    matrices.h1_standard_enrichment.Mult(h1_enrichment, local_result);
    local.h1 += h1_local * local_result;
    local_result.SetSize(nd_local.Size());
    matrices.nd_mass_standard_enrichment.Mult(nd_enrichment, local_result);
    local.nd_mass += nd_local * local_result;
    matrices.nd_curl_curl_standard_enrichment.Mult(nd_enrichment, local_result);
    local.nd_curl_curl += nd_local * local_result;
    local.h1_error +=
        BilinearErrorBound(h1_local, matrices.h1_estimated_absolute_error, h1_enrichment);
    local.nd_mass_error += BilinearErrorBound(
        nd_local, matrices.nd_mass_estimated_absolute_error, nd_enrichment);
    local.nd_curl_curl_error += BilinearErrorBound(
        nd_local, matrices.nd_curl_curl_estimated_absolute_error, nd_enrichment);

    const auto untransformed_mass = matrices.nd_mass_standard_enrichment;
    fem::singular::ApplyStandardDofTransformations(h1_dof_transformation,
                                                   nd_dof_transformation, matrices);
    for (int i = 0; i < untransformed_mass.Height(); i++)
    {
      for (int j = 0; j < untransformed_mass.Width(); j++)
      {
        has_nontrivial_face_transform =
            has_nontrivial_face_transform ||
            matrices.nd_mass_standard_enrichment(i, j) != untransformed_mass(i, j);
      }
    }

    h1_forward.AddSubMatrix(h1_vdofs, h1_enrichment_dofs, matrices.h1_standard_enrichment);
    h1_reverse.AddSubMatrix(h1_enrichment_dofs, h1_vdofs, matrices.h1_enrichment_standard);
    nd_mass_forward.AddSubMatrix(nd_vdofs, nd_enrichment_dofs,
                                 matrices.nd_mass_standard_enrichment);
    nd_mass_reverse.AddSubMatrix(nd_enrichment_dofs, nd_vdofs,
                                 matrices.nd_mass_enrichment_standard);
    nd_curl_curl_forward.AddSubMatrix(nd_vdofs, nd_enrichment_dofs,
                                      matrices.nd_curl_curl_standard_enrichment);
    nd_curl_curl_reverse.AddSubMatrix(nd_enrichment_dofs, nd_vdofs,
                                      matrices.nd_curl_curl_enrichment_standard);
  }
  REQUIRE(has_signed_nd_dof);
  REQUIRE(has_nontrivial_face_transform);

  h1_forward.Finalize();
  h1_reverse.Finalize();
  nd_mass_forward.Finalize();
  nd_mass_reverse.Finalize();
  nd_curl_curl_forward.Finalize();
  nd_curl_curl_reverse.Finalize();
  const mfem::SparseMatrix &h1_forward_const = h1_forward;
  const mfem::SparseMatrix &h1_reverse_const = h1_reverse;
  const mfem::SparseMatrix &nd_mass_forward_const = nd_mass_forward;
  const mfem::SparseMatrix &nd_mass_reverse_const = nd_mass_reverse;
  const mfem::SparseMatrix &nd_curl_curl_forward_const = nd_curl_curl_forward;
  const mfem::SparseMatrix &nd_curl_curl_reverse_const = nd_curl_curl_reverse;
  for (int i = 0; i < h1_forward.Height(); i++)
  {
    CHECK(h1_forward_const(i, 0) == h1_reverse_const(0, i));
  }
  for (int i = 0; i < nd_mass_forward.Height(); i++)
  {
    for (int j = 0; j < nd_mass_forward.Width(); j++)
    {
      CHECK(nd_mass_forward_const(i, j) == nd_mass_reverse_const(j, i));
      CHECK(nd_curl_curl_forward_const(i, j) == nd_curl_curl_reverse_const(j, i));
    }
  }

  CouplingEnergies global = local;
  mfem::Vector global_result(h1_space.GetVSize());
  h1_forward.Mult(h1_enrichment, global_result);
  global.h1 = h1_field * global_result;
  global_result.SetSize(nd_space.GetVSize());
  nd_mass_forward.Mult(nd_enrichment, global_result);
  global.nd_mass = nd_field * global_result;
  nd_curl_curl_forward.Mult(nd_enrichment, global_result);
  global.nd_curl_curl = nd_field * global_result;
  CHECK(std::abs(global.h1 - local.h1) <=
        128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(local.h1)));
  CHECK(std::abs(global.nd_mass - local.nd_mass) <=
        128.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(local.nd_mass)));
  CHECK(std::abs(global.nd_curl_curl - local.nd_curl_curl) <=
        128.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(local.nd_curl_curl)));
  return global;
}

}  // namespace

TEST_CASE("Triangular singular assembly preserves the complete exact sequence",
          "[singularelements][singularassembly][triangle][Serial]")
{
  constexpr double nu = 0.5;
  fem::singular::TriangleElementDofMap element_dofs;
  const fem::singular::TriangleBasis gradient_1{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, nu};
  const fem::singular::TriangleBasis gradient_2{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 2, 1}, 1, nu};
  const fem::singular::TriangleBasis rotational{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL, {0, 1, 2}, 1, nu};
  element_dofs.h1 = {{0, gradient_1}, {1, gradient_2}};
  element_dofs.nd = {{0, gradient_1}, {1, gradient_2}, {2, rotational}};

  auto mesh = AffineTriangleMesh({0.3, -0.2}, {1.7, 0.2}, {0.3, 1.4});
  mfem::H1_FECollection h1_collection(1, 2);
  mfem::ND_FECollection nd_collection(1, 2);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  auto &transformation = *mesh.GetElementTransformation(0);
  double jacobian_determinant;
  const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
      transformation, jacobian_determinant);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-10, 2.0e-10, 9};

  const auto enrichment = fem::singular::AssembleTriangleElementEnrichmentMatrices(
      element_dofs, grad_lambda, jacobian_determinant, options);
  const auto coupling = fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
      element_dofs, *h1_space.GetFE(0), *nd_space.GetFE(0), transformation, options);

  REQUIRE(enrichment.h1_diffusion.Height() == 2);
  REQUIRE(enrichment.h1_mass.Height() == 2);
  REQUIRE(enrichment.nd_mass.Height() == 3);
  REQUIRE(enrichment.nd_curl_curl.Height() == 3);
  CheckSymmetric(enrichment.h1_diffusion);
  CheckSymmetric(enrichment.h1_diffusion_estimated_absolute_error);
  CheckSymmetric(enrichment.h1_mass);
  CheckSymmetric(enrichment.h1_mass_estimated_absolute_error);
  CheckSymmetric(enrichment.nd_mass);
  CheckSymmetric(enrichment.nd_mass_estimated_absolute_error);
  CheckSymmetric(enrichment.nd_curl_curl);
  CheckSymmetric(enrichment.nd_curl_curl_estimated_absolute_error);
  CheckExactTranspose(coupling.h1_standard_enrichment, coupling.h1_enrichment_standard);
  CheckExactTranspose(coupling.h1_mass_standard_enrichment,
                      coupling.h1_mass_enrichment_standard);
  CheckExactTranspose(coupling.nd_mass_standard_enrichment,
                      coupling.nd_mass_enrichment_standard);
  CheckExactTranspose(coupling.nd_curl_curl_standard_enrichment,
                      coupling.nd_curl_curl_enrichment_standard);

  for (int row = 0; row < 2; row++)
  {
    for (int column = 0; column < 2; column++)
    {
      CHECK(enrichment.h1_diffusion(row, column) == enrichment.nd_mass(row, column));
      CHECK(enrichment.h1_diffusion_estimated_absolute_error(row, column) ==
            enrichment.nd_mass_estimated_absolute_error(row, column));
    }
    for (int column = 0; column < 3; column++)
    {
      CHECK(enrichment.nd_curl_curl(row, column) == 0.0);
      CHECK(enrichment.nd_curl_curl_estimated_absolute_error(row, column) == 0.0);
    }
  }
  CHECK(enrichment.nd_curl_curl(2, 2) > 0.0);
  CheckPositiveSemidefinite(enrichment.h1_diffusion,
                            enrichment.h1_diffusion_estimated_absolute_error);
  CheckPositiveSemidefinite(enrichment.h1_mass,
                            enrichment.h1_mass_estimated_absolute_error);
  CheckPositiveSemidefinite(enrichment.nd_mass,
                            enrichment.nd_mass_estimated_absolute_error);
  CheckPositiveSemidefinite(enrichment.nd_curl_curl,
                            enrichment.nd_curl_curl_estimated_absolute_error);

  mfem::DenseMatrix standard_nd_mass;
  mfem::VectorFEMassIntegrator nd_mass_integrator;
  nd_mass_integrator.AssembleElementMatrix(*nd_space.GetFE(0), transformation,
                                           standard_nd_mass);
  mfem::DenseMatrix standard_h1_diffusion;
  mfem::DiffusionIntegrator diffusion_integrator;
  diffusion_integrator.AssembleElementMatrix(*h1_space.GetFE(0), transformation,
                                             standard_h1_diffusion);
  mfem::DenseMatrix standard_h1_mass;
  mfem::MassIntegrator h1_mass_integrator;
  h1_mass_integrator.AssembleElementMatrix(*h1_space.GetFE(0), transformation,
                                           standard_h1_mass);
  mfem::DenseMatrix standard_gradient(nd_space.GetFE(0)->GetDof(),
                                      h1_space.GetFE(0)->GetDof());
  nd_space.GetFE(0)->ProjectGrad(*h1_space.GetFE(0), transformation, standard_gradient);

  const int standard_h1 = standard_h1_diffusion.Height();
  const int standard_nd = standard_nd_mass.Height();
  const int enriched_h1 = enrichment.h1_diffusion.Height();
  const int enriched_nd = enrichment.nd_mass.Height();
  Eigen::MatrixXd combined_mass =
      Eigen::MatrixXd::Zero(standard_nd + enriched_nd, standard_nd + enriched_nd);
  combined_mass.topLeftCorner(standard_nd, standard_nd) = ToEigen(standard_nd_mass);
  combined_mass.topRightCorner(standard_nd, enriched_nd) =
      ToEigen(coupling.nd_mass_standard_enrichment);
  combined_mass.bottomLeftCorner(enriched_nd, standard_nd) =
      ToEigen(coupling.nd_mass_enrichment_standard);
  combined_mass.bottomRightCorner(enriched_nd, enriched_nd) = ToEigen(enrichment.nd_mass);

  Eigen::MatrixXd combined_gradient =
      Eigen::MatrixXd::Zero(standard_nd + enriched_nd, standard_h1 + enriched_h1);
  combined_gradient.topLeftCorner(standard_nd, standard_h1) = ToEigen(standard_gradient);
  combined_gradient.block(standard_nd, standard_h1, enriched_h1, enriched_h1).setIdentity();
  Eigen::MatrixXd combined_diffusion =
      Eigen::MatrixXd::Zero(standard_h1 + enriched_h1, standard_h1 + enriched_h1);
  combined_diffusion.topLeftCorner(standard_h1, standard_h1) =
      ToEigen(standard_h1_diffusion);
  combined_diffusion.topRightCorner(standard_h1, enriched_h1) =
      ToEigen(coupling.h1_standard_enrichment);
  combined_diffusion.bottomLeftCorner(enriched_h1, standard_h1) =
      ToEigen(coupling.h1_enrichment_standard);
  combined_diffusion.bottomRightCorner(enriched_h1, enriched_h1) =
      ToEigen(enrichment.h1_diffusion);
  const Eigen::MatrixXd exact_sequence =
      combined_gradient.transpose() * combined_mass * combined_gradient;
  CHECK((combined_diffusion - exact_sequence).cwiseAbs().maxCoeff() < 3.0e-12);

  Eigen::MatrixXd combined_h1_mass =
      Eigen::MatrixXd::Zero(standard_h1 + enriched_h1, standard_h1 + enriched_h1);
  combined_h1_mass.topLeftCorner(standard_h1, standard_h1) = ToEigen(standard_h1_mass);
  combined_h1_mass.topRightCorner(standard_h1, enriched_h1) =
      ToEigen(coupling.h1_mass_standard_enrichment);
  combined_h1_mass.bottomLeftCorner(enriched_h1, standard_h1) =
      ToEigen(coupling.h1_mass_enrichment_standard);
  combined_h1_mass.bottomRightCorner(enriched_h1, enriched_h1) =
      ToEigen(enrichment.h1_mass);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> h1_mass_solver(combined_h1_mass);
  REQUIRE(h1_mass_solver.info() == Eigen::Success);
  CHECK(h1_mass_solver.eigenvalues().minCoeff() > 0.0);

  for (int row = 0; row < enriched_h1; row++)
  {
    for (int column = row; column < enriched_h1; column++)
    {
      const auto &row_basis = element_dofs.h1[row].basis;
      const auto &column_basis = element_dofs.h1[column].basis;
      const double direct = fem::singular::IntegrateReferenceTriangleNodeDuffy(
          fem::singular::H1DuffyReferenceOrder + 4, 0,
          fem::singular::TriangleDuffyRadialPower,
          [&](const fem::singular::TriangleBarycentricPoint &lambda)
          {
            return jacobian_determinant *
                   fem::singular::EvaluateTriangleNodeGradientPotential(
                       lambda, row_basis.nodes[0], row_basis.nodes[1], row_basis.nu) *
                   fem::singular::EvaluateTriangleNodeGradientPotential(
                       lambda, column_basis.nodes[0], column_basis.nodes[1],
                       column_basis.nu);
          });
      CAPTURE(row, column, enrichment.h1_mass(row, column), direct,
              enrichment.h1_mass_estimated_absolute_error(row, column));
      CHECK(std::abs(enrichment.h1_mass(row, column) - direct) <=
            enrichment.h1_mass_estimated_absolute_error(row, column) + 2.0e-13);
    }
  }
  mfem::IntegrationPoint point;
  mfem::Vector standard_shape(standard_h1);
  for (int standard = 0; standard < standard_h1; standard++)
  {
    for (int enrichment_index = 0; enrichment_index < enriched_h1; enrichment_index++)
    {
      const auto &basis = element_dofs.h1[enrichment_index].basis;
      const double direct = fem::singular::IntegrateReferenceTriangleNodeDuffy(
          fem::singular::H1DuffyReferenceOrder + 4, 0,
          fem::singular::TriangleDuffyRadialPower,
          [&](const fem::singular::TriangleBarycentricPoint &lambda)
          {
            point.Set2(lambda[1], lambda[2]);
            h1_space.GetFE(0)->CalcShape(point, standard_shape);
            return jacobian_determinant * standard_shape[standard] *
                   fem::singular::EvaluateTriangleNodeGradientPotential(
                       lambda, basis.nodes[0], basis.nodes[1], basis.nu);
          });
      CAPTURE(standard, enrichment_index,
              coupling.h1_mass_standard_enrichment(standard, enrichment_index), direct,
              coupling.h1_mass_estimated_absolute_error(standard, enrichment_index));
      CHECK(std::abs(coupling.h1_mass_standard_enrichment(standard, enrichment_index) -
                     direct) <=
            coupling.h1_mass_estimated_absolute_error(standard, enrichment_index) +
                2.0e-13);
    }
  }

  auto missing_gradient = element_dofs;
  missing_gradient.nd.erase(missing_gradient.nd.begin());
  CHECK_THROWS_AS(fem::singular::AssembleTriangleElementEnrichmentMatrices(
                      missing_gradient, grad_lambda, jacobian_determinant, options),
                  std::invalid_argument);
  const fem::singular::AdaptiveAssemblyOptions impossible{8, 1.0e-30, 1.0e-30, 1};
  CHECK_THROWS_AS(fem::singular::AssembleTriangleElementEnrichmentMatrices(
                      element_dofs, grad_lambda, jacobian_determinant, impossible),
                  std::runtime_error);
}

TEST_CASE("Quadratic affine simplices retain singular reference assembly",
          "[singularelements][singularassembly][curved][Serial]")
{
  constexpr double nu = 0.5;
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 8};

  fem::singular::TriangleElementDofMap triangle_dofs;
  const fem::singular::TriangleBasis triangle_gradient{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, nu};
  const fem::singular::TriangleBasis triangle_rotation{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL, {0, 1, 2}, 1, nu};
  triangle_dofs.h1 = {{0, triangle_gradient}};
  triangle_dofs.nd = {{0, triangle_gradient}, {1, triangle_rotation}};

  auto linear_triangle = AffineTriangleMesh({0.3, -0.2}, {1.7, 0.2}, {0.3, 1.4});
  auto quadratic_triangle = AffineTriangleMesh({0.3, -0.2}, {1.7, 0.2}, {0.3, 1.4});
  SetQuadraticGeometry(quadratic_triangle, false);
  auto &linear_triangle_transformation = *linear_triangle.GetElementTransformation(0);
  auto &quadratic_triangle_transformation = *quadratic_triangle.GetElementTransformation(0);
  REQUIRE(quadratic_triangle_transformation.OrderJ() > 0);
  CHECK(fem::singular::IsAffineElementTransformation(quadratic_triangle_transformation));
  double quadratic_triangle_jacobian;
  CHECK_NOTHROW(fem::singular::GetAffineTriangleBarycentricGradients(
      quadratic_triangle_transformation, quadratic_triangle_jacobian));
  const auto linear_triangle_matrices =
      fem::singular::AssembleTriangleElementEnrichmentMatrices(
          triangle_dofs, linear_triangle_transformation, options);
  const auto quadratic_triangle_matrices =
      fem::singular::AssembleTriangleElementEnrichmentMatrices(
          triangle_dofs, quadratic_triangle_transformation, options);
  CheckClose(quadratic_triangle_matrices.h1_diffusion,
             linear_triangle_matrices.h1_diffusion);
  CheckClose(quadratic_triangle_matrices.h1_mass, linear_triangle_matrices.h1_mass);
  CheckClose(quadratic_triangle_matrices.nd_mass, linear_triangle_matrices.nd_mass);
  CheckClose(quadratic_triangle_matrices.nd_curl_curl,
             linear_triangle_matrices.nd_curl_curl);
  mfem::H1_FECollection triangle_h1_collection(1, 2);
  mfem::ND_FECollection triangle_nd_collection(1, 2);
  mfem::FiniteElementSpace linear_triangle_h1_space(&linear_triangle,
                                                    &triangle_h1_collection);
  mfem::FiniteElementSpace linear_triangle_nd_space(&linear_triangle,
                                                    &triangle_nd_collection);
  mfem::FiniteElementSpace quadratic_triangle_h1_space(&quadratic_triangle,
                                                       &triangle_h1_collection);
  mfem::FiniteElementSpace quadratic_triangle_nd_space(&quadratic_triangle,
                                                       &triangle_nd_collection);
  const auto linear_triangle_coupling =
      fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
          triangle_dofs, *linear_triangle_h1_space.GetFE(0),
          *linear_triangle_nd_space.GetFE(0), linear_triangle_transformation, options);
  const auto quadratic_triangle_coupling =
      fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
          triangle_dofs, *quadratic_triangle_h1_space.GetFE(0),
          *quadratic_triangle_nd_space.GetFE(0), quadratic_triangle_transformation,
          options);
  CheckClose(quadratic_triangle_coupling.h1_standard_enrichment,
             linear_triangle_coupling.h1_standard_enrichment);
  CheckClose(quadratic_triangle_coupling.h1_mass_standard_enrichment,
             linear_triangle_coupling.h1_mass_standard_enrichment);
  CheckClose(quadratic_triangle_coupling.nd_mass_standard_enrichment,
             linear_triangle_coupling.nd_mass_standard_enrichment);
  CheckClose(quadratic_triangle_coupling.nd_curl_curl_standard_enrichment,
             linear_triangle_coupling.nd_curl_curl_standard_enrichment);

  auto linear_tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                                  {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  auto quadratic_tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                                     {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  SetQuadraticGeometry(quadratic_tetrahedron, false);
  const auto tetrahedron_gradient =
      fem::singular::EnumerateHigherOrderNodeGradientBases({0, 1, 2, 3}, 1, nu).front();
  const auto tetrahedron_rotation =
      fem::singular::EnumerateHigherOrderNodeRotationalBases({0, 1, 2, 3}, 1, nu).front();
  fem::singular::ElementDofMap tetrahedron_dofs;
  tetrahedron_dofs.h1 = {{0, tetrahedron_gradient}};
  tetrahedron_dofs.nd = {{0, tetrahedron_gradient}, {1, tetrahedron_rotation}};
  auto &linear_tetrahedron_transformation = *linear_tetrahedron.GetElementTransformation(0);
  auto &quadratic_tetrahedron_transformation =
      *quadratic_tetrahedron.GetElementTransformation(0);
  REQUIRE(quadratic_tetrahedron_transformation.OrderJ() > 0);
  CHECK(fem::singular::IsAffineElementTransformation(quadratic_tetrahedron_transformation));
  double quadratic_tetrahedron_jacobian;
  CHECK_NOTHROW(fem::singular::GetAffineBarycentricGradients(
      quadratic_tetrahedron_transformation, quadratic_tetrahedron_jacobian));
  const auto linear_tetrahedron_matrices = fem::singular::AssembleElementEnrichmentMatrices(
      tetrahedron_dofs, linear_tetrahedron_transformation, options);
  const auto quadratic_tetrahedron_matrices =
      fem::singular::AssembleElementEnrichmentMatrices(
          tetrahedron_dofs, quadratic_tetrahedron_transformation, options);
  CheckClose(quadratic_tetrahedron_matrices.h1_diffusion,
             linear_tetrahedron_matrices.h1_diffusion);
  CheckClose(quadratic_tetrahedron_matrices.h1_mass, linear_tetrahedron_matrices.h1_mass);
  CheckClose(quadratic_tetrahedron_matrices.nd_mass, linear_tetrahedron_matrices.nd_mass);
  CheckClose(quadratic_tetrahedron_matrices.nd_curl_curl,
             linear_tetrahedron_matrices.nd_curl_curl);
  mfem::H1_FECollection tetrahedron_h1_collection(1, 3);
  mfem::ND_FECollection tetrahedron_nd_collection(1, 3);
  mfem::FiniteElementSpace linear_tetrahedron_h1_space(&linear_tetrahedron,
                                                       &tetrahedron_h1_collection);
  mfem::FiniteElementSpace linear_tetrahedron_nd_space(&linear_tetrahedron,
                                                       &tetrahedron_nd_collection);
  mfem::FiniteElementSpace quadratic_tetrahedron_h1_space(&quadratic_tetrahedron,
                                                          &tetrahedron_h1_collection);
  mfem::FiniteElementSpace quadratic_tetrahedron_nd_space(&quadratic_tetrahedron,
                                                          &tetrahedron_nd_collection);
  const auto linear_tetrahedron_coupling =
      fem::singular::AssembleElementStandardEnrichmentMatrices(
          tetrahedron_dofs, *linear_tetrahedron_h1_space.GetFE(0),
          *linear_tetrahedron_nd_space.GetFE(0), linear_tetrahedron_transformation,
          options);
  const auto quadratic_tetrahedron_coupling =
      fem::singular::AssembleElementStandardEnrichmentMatrices(
          tetrahedron_dofs, *quadratic_tetrahedron_h1_space.GetFE(0),
          *quadratic_tetrahedron_nd_space.GetFE(0), quadratic_tetrahedron_transformation,
          options);
  CheckClose(quadratic_tetrahedron_coupling.h1_standard_enrichment,
             linear_tetrahedron_coupling.h1_standard_enrichment);
  CheckClose(quadratic_tetrahedron_coupling.h1_mass_standard_enrichment,
             linear_tetrahedron_coupling.h1_mass_standard_enrichment);
  CheckClose(quadratic_tetrahedron_coupling.nd_mass_standard_enrichment,
             linear_tetrahedron_coupling.nd_mass_standard_enrichment);
  CheckClose(quadratic_tetrahedron_coupling.nd_curl_curl_standard_enrichment,
             linear_tetrahedron_coupling.nd_curl_curl_standard_enrichment);
}

TEST_CASE("Numerically affine quadratic simplices retain singular reference assembly",
          "[singularelements][singularassembly][curved][Serial]")
{
  auto triangle = AffineTriangleMesh({0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0});
  SetQuadraticGeometry(triangle, false);
  mfem::VectorFunctionCoefficient triangle_map(
      2,
      [](const mfem::Vector &point, mfem::Vector &value)
      {
        value.SetSize(2);
        value = point;
        value[1] += 3.0e-9 * point[0] * point[1];
      });
  triangle.GetNodes()->ProjectCoefficient(triangle_map);
  auto &triangle_transformation = *triangle.GetElementTransformation(0);
  const double triangle_variation =
      fem::singular::GetElementTransformationRelativeJacobianVariation(
          triangle_transformation);
  CAPTURE(triangle_variation);
  CHECK(triangle_variation > 1.0e-9);
  CHECK(triangle_variation < fem::singular::AffineElementRelativeJacobianTolerance);
  CHECK_FALSE(
      fem::singular::IsAffineElementTransformation(triangle_transformation, 1.0e-9));
  CHECK(fem::singular::IsAffineElementTransformation(triangle_transformation));

  auto tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  SetQuadraticGeometry(tetrahedron, false);
  mfem::VectorFunctionCoefficient tetrahedron_map(
      3,
      [](const mfem::Vector &point, mfem::Vector &value)
      {
        value.SetSize(3);
        value = point;
        value[2] += 3.0e-9 * point[0] * point[1];
      });
  tetrahedron.GetNodes()->ProjectCoefficient(tetrahedron_map);
  auto &tetrahedron_transformation = *tetrahedron.GetElementTransformation(0);
  const double tetrahedron_variation =
      fem::singular::GetElementTransformationRelativeJacobianVariation(
          tetrahedron_transformation);
  CAPTURE(tetrahedron_variation);
  CHECK(tetrahedron_variation > 1.0e-9);
  CHECK(tetrahedron_variation < fem::singular::AffineElementRelativeJacobianTolerance);
  CHECK_FALSE(
      fem::singular::IsAffineElementTransformation(tetrahedron_transformation, 1.0e-9));
  CHECK(fem::singular::IsAffineElementTransformation(tetrahedron_transformation));
}

TEST_CASE("Tiny translated affine simplices include coordinate roundoff",
          "[singularelements][singularassembly][curved][Serial]")
{
  constexpr double h = 1.0e-8;
  auto triangle = AffineTriangleMesh({0.5, 0.25}, {0.5 + h, 0.25}, {0.5, 0.25 + h});
  SetQuadraticGeometry(triangle, false);
  auto &triangle_transformation = *triangle.GetElementTransformation(0);
  const double triangle_variation =
      fem::singular::GetElementTransformationRelativeJacobianVariation(
          triangle_transformation);
  const double triangle_allowance =
      fem::singular::GetElementTransformationAffineRoundoffAllowance(
          triangle_transformation);
  CAPTURE(triangle_variation, triangle_allowance);
  CHECK(std::isfinite(triangle_variation));
  CHECK(triangle_allowance > 0.0);
  CHECK(fem::singular::IsAffineElementTransformation(triangle_transformation));

  auto tetrahedron = AffineTetrahedronMesh({0.5, 0.25, 0.125}, {0.5 + h, 0.25, 0.125},
                                           {0.5, 0.25 + h, 0.125}, {0.5, 0.25, 0.125 + h});
  SetQuadraticGeometry(tetrahedron, false);
  auto &tetrahedron_transformation = *tetrahedron.GetElementTransformation(0);
  const double tetrahedron_variation =
      fem::singular::GetElementTransformationRelativeJacobianVariation(
          tetrahedron_transformation);
  const double tetrahedron_allowance =
      fem::singular::GetElementTransformationAffineRoundoffAllowance(
          tetrahedron_transformation);
  CAPTURE(tetrahedron_variation, tetrahedron_allowance);
  CHECK(std::isfinite(tetrahedron_variation));
  CHECK(tetrahedron_allowance > 0.0);
  CHECK(fem::singular::IsAffineElementTransformation(tetrahedron_transformation));
}

TEST_CASE("Singular geometry rejects inverted and folded simplex maps",
          "[singularelements][singularassembly][curved][Serial]")
{
  mfem::IntegrationPoint triangle_center;
  triangle_center.Set2(1.0 / 3.0, 1.0 / 3.0);
  auto inverted_triangle = AffineTriangleMesh({0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0});
  for (int d = 0; d < 2; d++)
  {
    std::swap(inverted_triangle.GetVertex(1)[d], inverted_triangle.GetVertex(2)[d]);
  }
  auto &inverted_triangle_transformation = *inverted_triangle.GetElementTransformation(0);
  inverted_triangle_transformation.SetIntPoint(&triangle_center);
  REQUIRE(inverted_triangle_transformation.Jacobian().Det() < 0.0);
  REQUIRE(inverted_triangle_transformation.Weight() < 0.0);
  double jacobian_determinant;
  CHECK_THROWS_AS(
      fem::singular::GetTriangleBarycentricGradients(inverted_triangle_transformation,
                                                     triangle_center, jacobian_determinant),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::IsAffineElementTransformation(inverted_triangle_transformation),
      std::invalid_argument);

  mfem::IntegrationPoint tetrahedron_center;
  tetrahedron_center.Set3(0.25, 0.25, 0.25);
  auto inverted_tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                                    {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  for (int d = 0; d < 3; d++)
  {
    std::swap(inverted_tetrahedron.GetVertex(1)[d], inverted_tetrahedron.GetVertex(2)[d]);
  }
  auto &inverted_tetrahedron_transformation =
      *inverted_tetrahedron.GetElementTransformation(0);
  inverted_tetrahedron_transformation.SetIntPoint(&tetrahedron_center);
  REQUIRE(inverted_tetrahedron_transformation.Jacobian().Det() < 0.0);
  REQUIRE(inverted_tetrahedron_transformation.Weight() < 0.0);
  CHECK_THROWS_AS(
      fem::singular::GetBarycentricGradients(inverted_tetrahedron_transformation,
                                             tetrahedron_center, jacobian_determinant),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::IsAffineElementTransformation(inverted_tetrahedron_transformation),
      std::invalid_argument);

  auto folded_triangle = AffineTriangleMesh({0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0});
  folded_triangle.SetCurvature(2, false, 2, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient folded_triangle_map(
      2,
      [](const mfem::Vector &point, mfem::Vector &value)
      {
        value.SetSize(2);
        value[0] = point[0];
        value[1] = (1.0 - 2.0 * point[0]) * point[1];
      });
  folded_triangle.GetNodes()->ProjectCoefficient(folded_triangle_map);
  auto &folded_triangle_transformation = *folded_triangle.GetElementTransformation(0);
  mfem::IntegrationPoint folded_triangle_point;
  folded_triangle_point.Set2(0.75, 0.1);
  folded_triangle_transformation.SetIntPoint(&folded_triangle_point);
  REQUIRE(folded_triangle_transformation.Jacobian().Det() < 0.0);
  CHECK_THROWS_AS(
      fem::singular::GetTriangleBarycentricGradients(
          folded_triangle_transformation, folded_triangle_point, jacobian_determinant),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::IsAffineElementTransformation(folded_triangle_transformation),
      std::invalid_argument);
}

TEST_CASE("Curved simplex singular matrices use pointwise physical geometry",
          "[singularelements][singularassembly][curved][Serial]")
{
  constexpr double nu = 0.5;
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};

  auto triangle = AffineTriangleMesh({0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0});
  SetQuadraticGeometry(triangle, true);
  auto &triangle_transformation = *triangle.GetElementTransformation(0);
  REQUIRE_FALSE(fem::singular::IsAffineElementTransformation(triangle_transformation));
  double affine_jacobian;
  CHECK_THROWS_AS(fem::singular::GetAffineTriangleBarycentricGradients(
                      triangle_transformation, affine_jacobian),
                  std::invalid_argument);

  const fem::singular::TriangleBasis triangle_gradient{
      fem::singular::HigherOrderBasisFamily::NODE_GRADIENT, {0, 1, 2}, 1, nu};
  const fem::singular::TriangleBasis triangle_rotation{
      fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL, {0, 1, 2}, 1, nu};
  fem::singular::TriangleElementDofMap triangle_dofs;
  triangle_dofs.h1 = {{0, triangle_gradient}};
  triangle_dofs.nd = {{0, triangle_gradient}, {1, triangle_rotation}};
  const auto triangle_matrices = fem::singular::AssembleTriangleElementEnrichmentMatrices(
      triangle_dofs, triangle_transformation, options);
  mfem::H1_FECollection triangle_h1_collection(1, 2);
  mfem::ND_FECollection triangle_nd_collection(1, 2);
  mfem::FiniteElementSpace triangle_h1_space(&triangle, &triangle_h1_collection);
  mfem::FiniteElementSpace triangle_nd_space(&triangle, &triangle_nd_collection);
  const auto triangle_coupling =
      fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
          triangle_dofs, *triangle_h1_space.GetFE(0), *triangle_nd_space.GetFE(0),
          triangle_transformation, options);

  const auto evaluate_triangle =
      [](const fem::singular::TriangleBarycentricPoint &lambda,
         const fem::singular::TriangleBarycentricGradients &grad_lambda,
         const fem::singular::TriangleBasis &basis)
  {
    return basis.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT
               ? fem::singular::EvaluateTriangleNodeGradient(
                     lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nu)
               : fem::singular::EvaluateTriangleNodeRotational(
                     lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2],
                     basis.nu);
  };
  for (int row = 0; row < 2; row++)
  {
    for (int column = row; column < 2; column++)
    {
      long double mass = 0.0L;
      long double curl_curl = 0.0L;
      mfem::IntegrationPoint point;
      fem::singular::ForEachReferenceTriangleNodeDuffyQuadraturePoint(
          fem::singular::H1DuffyReferenceOrder + 2, 0,
          fem::singular::TriangleDuffyRadialPower,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight)
          {
            point.Set2(lambda[1], lambda[2]);
            double jacobian_determinant;
            const auto grad_lambda = fem::singular::GetTriangleBarycentricGradients(
                triangle_transformation, point, jacobian_determinant);
            const auto row_value =
                evaluate_triangle(lambda, grad_lambda, triangle_dofs.nd[row].basis);
            const auto column_value =
                evaluate_triangle(lambda, grad_lambda, triangle_dofs.nd[column].basis);
            mass += static_cast<long double>(weight * jacobian_determinant) *
                    Dot(row_value.value, column_value.value);
            curl_curl += static_cast<long double>(weight * jacobian_determinant) *
                         row_value.curl * column_value.curl;
          });
      CAPTURE(row, column);
      CHECK(std::abs(triangle_matrices.nd_mass(row, column) - static_cast<double>(mass)) <=
            triangle_matrices.nd_mass_estimated_absolute_error(row, column) + 2.0e-11);
      CHECK(std::abs(triangle_matrices.nd_curl_curl(row, column) -
                     static_cast<double>(curl_curl)) <=
            triangle_matrices.nd_curl_curl_estimated_absolute_error(row, column) + 2.0e-11);
    }
  }
  CHECK(triangle_matrices.h1_diffusion(0, 0) == triangle_matrices.nd_mass(0, 0));
  mfem::DenseMatrix triangle_discrete_gradient(triangle_nd_space.GetFE(0)->GetDof(),
                                               triangle_h1_space.GetFE(0)->GetDof());
  triangle_nd_space.GetFE(0)->ProjectGrad(
      *triangle_h1_space.GetFE(0), triangle_transformation, triangle_discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < triangle_h1_space.GetFE(0)->GetDof();
       standard_h1++)
  {
    double exact_sequence_value = 0.0;
    for (int standard_nd = 0; standard_nd < triangle_nd_space.GetFE(0)->GetDof();
         standard_nd++)
    {
      exact_sequence_value += triangle_discrete_gradient(standard_nd, standard_h1) *
                              triangle_coupling.nd_mass_standard_enrichment(standard_nd, 0);
    }
    CHECK(std::abs(triangle_coupling.h1_standard_enrichment(standard_h1, 0) -
                   exact_sequence_value) <=
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(exact_sequence_value)));
  }

  auto tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  SetQuadraticGeometry(tetrahedron, true);
  auto &tetrahedron_transformation = *tetrahedron.GetElementTransformation(0);
  REQUIRE_FALSE(fem::singular::IsAffineElementTransformation(tetrahedron_transformation));
  CHECK_THROWS_AS(fem::singular::GetAffineBarycentricGradients(tetrahedron_transformation,
                                                               affine_jacobian),
                  std::invalid_argument);

  const auto tetrahedron_gradient =
      fem::singular::EnumerateHigherOrderNodeGradientBases({0, 1, 2, 3}, 1, nu).front();
  const auto tetrahedron_rotation =
      fem::singular::EnumerateHigherOrderNodeRotationalBases({0, 1, 2, 3}, 1, nu).front();
  fem::singular::ElementDofMap tetrahedron_dofs;
  tetrahedron_dofs.h1 = {{0, tetrahedron_gradient}};
  tetrahedron_dofs.nd = {{0, tetrahedron_gradient}, {1, tetrahedron_rotation}};
  const auto tetrahedron_matrices = fem::singular::AssembleElementEnrichmentMatrices(
      tetrahedron_dofs, tetrahedron_transformation, options);
  mfem::H1_FECollection tetrahedron_h1_collection(1, 3);
  mfem::ND_FECollection tetrahedron_nd_collection(1, 3);
  mfem::FiniteElementSpace tetrahedron_h1_space(&tetrahedron, &tetrahedron_h1_collection);
  mfem::FiniteElementSpace tetrahedron_nd_space(&tetrahedron, &tetrahedron_nd_collection);
  const auto tetrahedron_coupling =
      fem::singular::AssembleElementStandardEnrichmentMatrices(
          tetrahedron_dofs, *tetrahedron_h1_space.GetFE(0), *tetrahedron_nd_space.GetFE(0),
          tetrahedron_transformation, options);
  const auto tetrahedron_h1_matrices = fem::singular::AssembleElementH1EnrichmentMatrices(
      tetrahedron_dofs, *tetrahedron_h1_space.GetFE(0), tetrahedron_transformation,
      options);

  std::array<long double, 2> direct_mass{};
  mfem::IntegrationPoint point;
  fem::singular::ForEachReferenceTetrahedronNodeDuffyQuadraturePoint(
      fem::singular::H1DuffyReferenceOrder + 2, 0, fem::singular::H1DuffyRadialPower,
      [&](const fem::singular::BarycentricPoint &lambda, double weight)
      {
        point.Set3(lambda[1], lambda[2], lambda[3]);
        double jacobian_determinant;
        const auto grad_lambda = fem::singular::GetBarycentricGradients(
            tetrahedron_transformation, point, jacobian_determinant);
        const auto gradient = fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda,
                                                                      tetrahedron_gradient);
        const auto rotation = fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda,
                                                                      tetrahedron_rotation);
        direct_mass[0] += static_cast<long double>(weight * jacobian_determinant) *
                          Dot(gradient.value, gradient.value);
        direct_mass[1] += static_cast<long double>(weight * jacobian_determinant) *
                          Dot(gradient.value, rotation.value);
      });
  for (int column = 0; column < 2; column++)
  {
    CAPTURE(column);
    CHECK(std::abs(tetrahedron_matrices.nd_mass(0, column) -
                   static_cast<double>(direct_mass[column])) <=
          tetrahedron_matrices.nd_mass_estimated_absolute_error(0, column) + 2.0e-8);
  }
  CHECK(tetrahedron_matrices.h1_diffusion(0, 0) == tetrahedron_matrices.nd_mass(0, 0));
  CHECK(std::abs(tetrahedron_h1_matrices.enrichment_enrichment(0, 0) -
                 tetrahedron_matrices.h1_diffusion(0, 0)) <=
        tetrahedron_h1_matrices.enrichment_enrichment_estimated_absolute_error(0, 0) +
            tetrahedron_matrices.h1_diffusion_estimated_absolute_error(0, 0) + 2.0e-8);
  CHECK(tetrahedron_matrices.nd_curl_curl(0, 0) == 0.0);
  CHECK(tetrahedron_matrices.nd_curl_curl(0, 1) == 0.0);
  CHECK(tetrahedron_matrices.nd_curl_curl(1, 1) > 0.0);
  mfem::DenseMatrix tetrahedron_discrete_gradient(tetrahedron_nd_space.GetFE(0)->GetDof(),
                                                  tetrahedron_h1_space.GetFE(0)->GetDof());
  tetrahedron_nd_space.GetFE(0)->ProjectGrad(*tetrahedron_h1_space.GetFE(0),
                                             tetrahedron_transformation,
                                             tetrahedron_discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < tetrahedron_h1_space.GetFE(0)->GetDof();
       standard_h1++)
  {
    double exact_sequence_value = 0.0;
    for (int standard_nd = 0; standard_nd < tetrahedron_nd_space.GetFE(0)->GetDof();
         standard_nd++)
    {
      exact_sequence_value +=
          tetrahedron_discrete_gradient(standard_nd, standard_h1) *
          tetrahedron_coupling.nd_mass_standard_enrichment(standard_nd, 0);
    }
    CHECK(std::abs(tetrahedron_coupling.h1_standard_enrichment(standard_h1, 0) -
                   exact_sequence_value) <=
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(exact_sequence_value)));
  }
}

TEST_CASE("Adaptive singular element local matrices preserve exact-sequence structure",
          "[singularelements][singularassembly][Serial]")
{
  constexpr double nu = 0.5;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const auto gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu);
  REQUIRE(gradients.size() >= 2);
  REQUIRE(rotations.size() >= 2);

  fem::singular::ElementDofMap element_dofs;
  for (int i = 0; i < 2; i++)
  {
    element_dofs.h1.push_back({static_cast<std::size_t>(i), gradients[i]});
    element_dofs.nd.push_back({static_cast<std::size_t>(i), gradients[i]});
  }
  for (int i = 0; i < 2; i++)
  {
    element_dofs.nd.push_back({static_cast<std::size_t>(i + 2), rotations[i]});
  }

  double jacobian_determinant;
  const auto grad_lambda = AffineBarycentricGradients(
      {1.7, 0.2, -0.1}, {0.3, 1.4, 0.25}, {-0.2, 0.1, 0.9}, jacobian_determinant);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto matrices = fem::singular::AssembleElementEnrichmentMatrices(
      element_dofs, grad_lambda, jacobian_determinant, options);

  CHECK(matrices.total_quadrature_leaf_count > 0);
  CHECK(matrices.maximum_subdivision_depth <= options.maximum_subdivisions);
  CHECK(matrices.h1_diffusion.Height() == 2);
  CHECK(matrices.nd_mass.Height() == 4);
  CheckSymmetric(matrices.h1_diffusion);
  CheckSymmetric(matrices.h1_diffusion_estimated_absolute_error);
  CheckSymmetric(matrices.nd_mass);
  CheckSymmetric(matrices.nd_mass_estimated_absolute_error);
  CheckSymmetric(matrices.nd_curl_curl);
  CheckSymmetric(matrices.nd_curl_curl_estimated_absolute_error);

  for (int row = 0; row < matrices.nd_mass.Height(); row++)
  {
    for (int column = 0; column < matrices.nd_mass.Width(); column++)
    {
      const auto direct = DirectPhysicalPairIntegral(element_dofs.nd[row].basis,
                                                     element_dofs.nd[column].basis,
                                                     grad_lambda, jacobian_determinant);
      CAPTURE(row, column, matrices.nd_mass(row, column), direct[0],
              matrices.nd_mass_estimated_absolute_error(row, column),
              matrices.nd_curl_curl(row, column), direct[1],
              matrices.nd_curl_curl_estimated_absolute_error(row, column));
      CHECK(std::abs(matrices.nd_mass(row, column) - direct[0]) <=
            2.0 * matrices.nd_mass_estimated_absolute_error(row, column) + 2.0e-10);
      CHECK(std::abs(matrices.nd_curl_curl(row, column) - direct[1]) <=
            2.0 * matrices.nd_curl_curl_estimated_absolute_error(row, column) + 2.0e-10);
    }
  }

  for (int row = 0; row < matrices.h1_diffusion.Height(); row++)
  {
    for (int column = 0; column < matrices.h1_diffusion.Width(); column++)
    {
      CHECK(matrices.h1_diffusion(row, column) == matrices.nd_mass(row, column));
      CHECK(matrices.h1_diffusion_estimated_absolute_error(row, column) ==
            matrices.nd_mass_estimated_absolute_error(row, column));
    }
  }
  for (int gradient = 0; gradient < 2; gradient++)
  {
    for (int column = 0; column < matrices.nd_curl_curl.Width(); column++)
    {
      CHECK(matrices.nd_curl_curl(gradient, column) == 0.0);
      CHECK(matrices.nd_curl_curl_estimated_absolute_error(gradient, column) == 0.0);
    }
  }

  CheckPositiveSemidefinite(matrices.h1_diffusion,
                            matrices.h1_diffusion_estimated_absolute_error);
  CheckPositiveSemidefinite(matrices.nd_mass, matrices.nd_mass_estimated_absolute_error);
  CheckPositiveSemidefinite(matrices.nd_curl_curl,
                            matrices.nd_curl_curl_estimated_absolute_error);

  auto missing_gradient = element_dofs;
  missing_gradient.nd.erase(missing_gradient.nd.begin());
  CHECK_THROWS_AS(fem::singular::AssembleElementEnrichmentMatrices(
                      missing_gradient, grad_lambda, jacobian_determinant, options),
                  std::invalid_argument);

  const fem::singular::AdaptiveAssemblyOptions invalid{8, 0.0, 0.0, 1};
  CHECK_THROWS_AS(fem::singular::AssembleElementEnrichmentMatrices(
                      element_dofs, grad_lambda, jacobian_determinant, invalid),
                  std::invalid_argument);

  const fem::singular::AdaptiveAssemblyOptions impossible{8, 1.0e-30, 1.0e-30, 1};
  CHECK_THROWS_AS(fem::singular::AssembleElementEnrichmentMatrices(
                      element_dofs, grad_lambda, jacobian_determinant, impossible),
                  std::runtime_error);
}

TEST_CASE("MFEM tetrahedral edge bases have the opposite sign from equation 5",
          "[singularelements][singularassembly][Serial]")
{
  mfem::ND_TetrahedronElement nd_fe(1);
  const fem::singular::BarycentricPoint lambda{0.37, 0.21, 0.18, 0.24};
  mfem::IntegrationPoint point;
  point.Set3(lambda[1], lambda[2], lambda[3]);
  mfem::DenseMatrix shape(nd_fe.GetDof(), 3);
  mfem::DenseMatrix curl(nd_fe.GetDof(), 3);
  nd_fe.CalcVShape(point, shape);
  nd_fe.CalcCurlShape(point, curl);
  constexpr std::array<std::array<int, 2>, 6> edges{
      std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2}, std::array<int, 2>{0, 3},
      std::array<int, 2>{1, 2}, std::array<int, 2>{1, 3}, std::array<int, 2>{2, 3}};
  const auto &grad_lambda = fem::singular::ReferenceBarycentricGradients();
  for (int edge = 0; edge < 6; edge++)
  {
    const auto paper = fem::singular::EvaluateStandardEdge(lambda, grad_lambda,
                                                           edges[edge][0], edges[edge][1]);
    for (int d = 0; d < 3; d++)
    {
      CAPTURE(edge, d);
      CHECK(std::abs(shape(edge, d) + paper.value[d]) <= 2.0e-15);
      CHECK(std::abs(curl(edge, d) + paper.curl[d]) <= 2.0e-15);
    }
  }
}

TEST_CASE("First-order H1 assembly uses estimated cached Duffy entries",
          "[singularelements][singularassembly][Serial]")
{
  constexpr double nu = 0.5;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const auto node_gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto edge_gradients =
      fem::singular::EnumerateHigherOrderEdgeGradientBases(canonical_nodes, 1, nu);
  REQUIRE(node_gradients.size() >= 2);
  REQUIRE(edge_gradients.size() >= 2);

  const fem::singular::Vector3 origin{0.4, -0.3, 0.2};
  const fem::singular::Vector3 a{1.7, 0.2, -0.1};
  const fem::singular::Vector3 b{0.3, 1.4, 0.25};
  const fem::singular::Vector3 c{-0.2, 0.1, 0.9};
  auto mesh = AffineTetrahedronMesh(origin, a, b, c);
  auto &transformation = *mesh.GetElementTransformation(0);
  double jacobian_determinant;
  const auto grad_lambda =
      PhysicalBarycentricGradients(transformation, jacobian_determinant);
  mfem::H1_TetrahedronElement h1_fe(1);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};

  const auto check_aligned_family =
      [&](const std::vector<fem::singular::HigherOrderBasis> &bases, bool node_aligned)
  {
    fem::singular::ElementDofMap element_dofs;
    for (int i = 0; i < 2; i++)
    {
      element_dofs.h1.push_back({static_cast<std::size_t>(i), bases[i]});
    }
    const auto matrices = fem::singular::AssembleElementH1EnrichmentMatrices(
        element_dofs, h1_fe, transformation, options);
    CHECK(matrices.total_quadrature_leaf_count == 0);
    CHECK(matrices.maximum_subdivision_depth == 0);

    for (int row = 0; row < 2; row++)
    {
      for (int column = 0; column < 2; column++)
      {
        const double direct = DirectH1EnrichmentPairDuffyIntegral(
            bases[row], bases[column], grad_lambda, jacobian_determinant,
            fem::singular::H1DuffyReferenceOrder);
        const auto lower_reference = node_aligned
                                         ? fem::singular::ComputeNodeDuffyReferenceIntegral(
                                               fem::singular::ReferenceBasis{bases[row]},
                                               fem::singular::ReferenceBasis{bases[column]},
                                               fem::singular::H1DuffyComparisonOrder,
                                               fem::singular::H1DuffyRadialPower)
                                         : fem::singular::ComputeEdgeDuffyReferenceIntegral(
                                               fem::singular::ReferenceBasis{bases[row]},
                                               fem::singular::ReferenceBasis{bases[column]},
                                               fem::singular::H1DuffyComparisonOrder,
                                               fem::singular::H1DuffyRadialPower);
        const double lower =
            fem::singular::ContractMass(lower_reference, grad_lambda, jacobian_determinant);
        const double roundoff = 256.0 * std::numeric_limits<double>::epsilon() *
                                std::max({1.0, std::abs(direct), std::abs(lower)});
        CAPTURE(node_aligned, row, column, matrices.enrichment_enrichment(row, column),
                direct, lower,
                matrices.enrichment_enrichment_estimated_absolute_error(row, column));
        CHECK(std::abs(matrices.enrichment_enrichment(row, column) - direct) <= roundoff);
        CHECK(std::abs(matrices.enrichment_enrichment(row, column) - lower) <=
              matrices.enrichment_enrichment_estimated_absolute_error(row, column) +
                  roundoff);
      }
    }

    for (int standard = 0; standard < h1_fe.GetDof(); standard++)
    {
      for (int enrichment = 0; enrichment < 2; enrichment++)
      {
        const double direct = DirectH1StandardEnrichmentDuffyIntegral(
            h1_fe, transformation, standard, bases[enrichment], grad_lambda,
            jacobian_determinant, fem::singular::H1DuffyReferenceOrder);
        const auto standard_basis =
            fem::singular::ReferenceBasis{fem::singular::MakeStandardH1Gradient(standard)};
        const auto singular_basis = fem::singular::ReferenceBasis{bases[enrichment]};
        const auto lower_reference =
            node_aligned
                ? fem::singular::ComputeNodeDuffyReferenceIntegral(
                      standard_basis, singular_basis, fem::singular::H1DuffyComparisonOrder,
                      fem::singular::H1DuffyRadialPower)
                : fem::singular::ComputeEdgeDuffyReferenceIntegral(
                      standard_basis, singular_basis, fem::singular::H1DuffyComparisonOrder,
                      fem::singular::H1DuffyRadialPower);
        const double lower =
            fem::singular::ContractMass(lower_reference, grad_lambda, jacobian_determinant);
        const double roundoff = 256.0 * std::numeric_limits<double>::epsilon() *
                                std::max({1.0, std::abs(direct), std::abs(lower)});
        CAPTURE(
            node_aligned, standard, enrichment,
            matrices.standard_enrichment(standard, enrichment), direct, lower,
            matrices.standard_enrichment_estimated_absolute_error(standard, enrichment));
        CHECK(std::abs(matrices.standard_enrichment(standard, enrichment) - direct) <=
              roundoff);
        CHECK(std::abs(matrices.standard_enrichment(standard, enrichment) - lower) <=
              matrices.standard_enrichment_estimated_absolute_error(standard, enrichment) +
                  roundoff);
      }
    }
  };
  check_aligned_family(node_gradients, true);
  check_aligned_family(edge_gradients, false);

  fem::singular::ElementDofMap mixed_features;
  mixed_features.h1.push_back({0, node_gradients[0]});
  const auto second_node_gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases({1, 0, 2, 3}, 1, nu);
  REQUIRE_FALSE(second_node_gradients.empty());
  mixed_features.h1.push_back({1, second_node_gradients[0]});
  const auto mixed = fem::singular::AssembleElementH1EnrichmentMatrices(
      mixed_features, h1_fe, transformation, options);
  CHECK(mixed.total_quadrature_leaf_count == 0);
  CHECK(mixed.maximum_subdivision_depth == 0);
  const auto mixed_reference = fem::singular::ComputePartitionedDuffyReferenceIntegral(
      fem::singular::ReferenceBasis{node_gradients[0]},
      fem::singular::ReferenceBasis{second_node_gradients[0]},
      fem::singular::H1DuffyReferenceOrder, fem::singular::H1DuffyRadialPower);
  const auto mixed_comparison = fem::singular::ComputePartitionedDuffyReferenceIntegral(
      fem::singular::ReferenceBasis{node_gradients[0]},
      fem::singular::ReferenceBasis{second_node_gradients[0]},
      fem::singular::H1DuffyComparisonOrder, fem::singular::H1DuffyRadialPower);
  const auto mixed_higher_order = fem::singular::ComputePartitionedDuffyReferenceIntegral(
      fem::singular::ReferenceBasis{node_gradients[0]},
      fem::singular::ReferenceBasis{second_node_gradients[0]},
      fem::singular::H1DuffyReferenceOrder + 2, fem::singular::H1DuffyRadialPower);
  const double mixed_value =
      fem::singular::ContractMass(mixed_reference, grad_lambda, jacobian_determinant);
  const double mixed_comparison_value =
      fem::singular::ContractMass(mixed_comparison, grad_lambda, jacobian_determinant);
  const double mixed_higher_order_value =
      fem::singular::ContractMass(mixed_higher_order, grad_lambda, jacobian_determinant);
  const double mixed_roundoff =
      256.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(mixed_value), std::abs(mixed_comparison_value),
                std::abs(mixed_higher_order_value)});
  CAPTURE(mixed.enrichment_enrichment(0, 1), mixed_value, mixed_comparison_value,
          mixed_higher_order_value,
          mixed.enrichment_enrichment_estimated_absolute_error(0, 1));
  CHECK(std::abs(mixed.enrichment_enrichment(0, 1) - mixed_value) <= mixed_roundoff);
  CHECK(std::abs(mixed.enrichment_enrichment(0, 1) - mixed_comparison_value) <=
        mixed.enrichment_enrichment_estimated_absolute_error(0, 1) + mixed_roundoff);
  CHECK(std::abs(mixed.enrichment_enrichment(0, 1) - mixed_higher_order_value) <=
        mixed.enrichment_enrichment_estimated_absolute_error(0, 1) + mixed_roundoff);

  const fem::singular::AdaptiveAssemblyOptions impossible{8, 1.0e-30, 1.0e-30, 1};
  fem::singular::ElementDofMap one_node;
  one_node.h1.push_back({0, node_gradients[0]});
  CHECK_THROWS_WITH(fem::singular::AssembleElementH1EnrichmentMatrices(
                        one_node, h1_fe, transformation, impossible),
                    Catch::Matchers::ContainsSubstring("Duffy"));

  auto one_tetrahedron = AffineTetrahedronMesh({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                               {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
  const std::vector<fem::singular::GlobalVertexId> one_vertex_ids{0, 1, 2, 3};
  const auto one_topology = SharedFaceDofTopology(one_tetrahedron, one_vertex_ids);
  mfem::H1_FECollection one_collection(1, 3);
  mfem::FiniteElementSpace one_space(&one_tetrahedron, &one_collection);
  const std::vector<fem::singular::IsotropicMaterialCoefficients> one_material{{1.0, 1.0}};
  const auto one_sparse = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
      one_topology, one_space, one_material, options);

  auto repeated_tetrahedra = SharedFaceTetrahedronMesh(false);
  std::vector<fem::singular::GlobalVertexId> repeated_vertex_ids(
      repeated_tetrahedra.GetNV());
  for (int vertex = 0; vertex < repeated_tetrahedra.GetNV(); vertex++)
  {
    repeated_vertex_ids[vertex] = vertex;
  }
  const auto repeated_topology =
      SharedFaceDofTopology(repeated_tetrahedra, repeated_vertex_ids);
  mfem::H1_FECollection repeated_collection(1, 3);
  mfem::FiniteElementSpace repeated_space(&repeated_tetrahedra, &repeated_collection);
  const std::vector<fem::singular::IsotropicMaterialCoefficients> repeated_materials{
      {1.0, 1.0}, {1.0, 1.0}};
  const auto repeated_sparse = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
      repeated_topology, repeated_space, repeated_materials, options);

  CHECK(one_sparse.total_quadrature_leaf_count == 0);
  CHECK(repeated_sparse.total_quadrature_leaf_count == 0);
  CHECK(one_sparse.duffy_reference_table_entries > 0);
  CHECK(repeated_sparse.duffy_reference_table_entries ==
        one_sparse.duffy_reference_table_entries);
  CHECK(repeated_sparse.duffy_reference_cache_hits > one_sparse.duffy_reference_cache_hits);
}

TEST_CASE("MFEM standard-enrichment coupling preserves the exact sequence",
          "[singularelements][singularassembly][Serial]")
{
  constexpr double nu = 0.5;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const auto gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu);
  REQUIRE(gradients.size() >= 2);
  REQUIRE(rotations.size() >= 2);

  fem::singular::ElementDofMap element_dofs;
  for (int i = 0; i < 2; i++)
  {
    element_dofs.h1.push_back({static_cast<std::size_t>(i), gradients[i]});
    element_dofs.nd.push_back({static_cast<std::size_t>(i), gradients[i]});
  }
  for (int i = 0; i < 2; i++)
  {
    element_dofs.nd.push_back({static_cast<std::size_t>(i + 2), rotations[i]});
  }

  const fem::singular::Vector3 origin{0.4, -0.3, 0.2};
  const fem::singular::Vector3 a{1.7, 0.2, -0.1};
  const fem::singular::Vector3 b{0.3, 1.4, 0.25};
  const fem::singular::Vector3 c{-0.2, 0.1, 0.9};
  double jacobian_determinant;
  const auto grad_lambda = AffineBarycentricGradients(a, b, c, jacobian_determinant);
  auto mesh = AffineTetrahedronMesh(origin, a, b, c);
  auto &transformation = *mesh.GetElementTransformation(0);
  mfem::H1_TetrahedronElement h1_fe(2);
  mfem::ND_TetrahedronElement nd_fe(2);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto matrices = fem::singular::AssembleElementStandardEnrichmentMatrices(
      element_dofs, h1_fe, nd_fe, transformation, options);
  const auto h1_matrices = fem::singular::AssembleElementH1EnrichmentMatrices(
      element_dofs, h1_fe, transformation, options);

  CHECK(matrices.total_quadrature_leaf_count > 0);
  CHECK(matrices.maximum_subdivision_depth <= options.maximum_subdivisions);
  CHECK(matrices.h1_standard_enrichment.Height() == h1_fe.GetDof());
  CHECK(matrices.h1_standard_enrichment.Width() == 2);
  CHECK(matrices.nd_mass_standard_enrichment.Height() == nd_fe.GetDof());
  CHECK(matrices.nd_mass_standard_enrichment.Width() == 4);
  CheckExactTranspose(matrices.h1_standard_enrichment, matrices.h1_enrichment_standard);
  CheckExactTranspose(matrices.nd_mass_standard_enrichment,
                      matrices.nd_mass_enrichment_standard);
  CheckExactTranspose(matrices.nd_curl_curl_standard_enrichment,
                      matrices.nd_curl_curl_enrichment_standard);
  CHECK(h1_matrices.total_quadrature_leaf_count > 0);
  CHECK(h1_matrices.total_quadrature_leaf_count < matrices.total_quadrature_leaf_count);
  CHECK(h1_matrices.maximum_subdivision_depth <= options.maximum_subdivisions);
  CheckSymmetric(h1_matrices.enrichment_enrichment);
  CheckSymmetric(h1_matrices.enrichment_enrichment_estimated_absolute_error);
  CheckExactTranspose(h1_matrices.standard_enrichment, h1_matrices.enrichment_standard);

  for (int row = 0; row < h1_matrices.enrichment_enrichment.Height(); row++)
  {
    for (int column = 0; column < h1_matrices.enrichment_enrichment.Width(); column++)
    {
      const double direct = DirectPhysicalPairIntegral(
          element_dofs.h1[row].basis, element_dofs.h1[column].basis, grad_lambda,
          jacobian_determinant)[0];
      CAPTURE(row, column, h1_matrices.enrichment_enrichment(row, column), direct,
              h1_matrices.enrichment_enrichment_estimated_absolute_error(row, column));
      CHECK(std::abs(h1_matrices.enrichment_enrichment(row, column) - direct) <=
            2.0 * h1_matrices.enrichment_enrichment_estimated_absolute_error(row, column) +
                1.0e-12);
    }
  }
  for (int standard = 0; standard < h1_fe.GetDof(); standard++)
  {
    for (int enrichment = 0; enrichment < h1_matrices.standard_enrichment.Width();
         enrichment++)
    {
      const double direct = DirectH1StandardEnrichmentIntegral(
          h1_fe, transformation, standard, element_dofs.h1[enrichment].basis, grad_lambda,
          jacobian_determinant);
      CAPTURE(
          standard, enrichment, h1_matrices.standard_enrichment(standard, enrichment),
          direct,
          h1_matrices.standard_enrichment_estimated_absolute_error(standard, enrichment));
      CHECK(std::abs(h1_matrices.standard_enrichment(standard, enrichment) - direct) <=
            2.0 * h1_matrices.standard_enrichment_estimated_absolute_error(standard,
                                                                           enrichment) +
                1.0e-12);
    }
  }

  const fem::singular::AdaptiveAssemblyOptions impossible_coupling{8, 1.0e-30, 1.0e-30, 1};
  CHECK_THROWS_WITH(fem::singular::AssembleElementStandardEnrichmentMatrices(
                        element_dofs, h1_fe, nd_fe, transformation, impossible_coupling),
                    Catch::Matchers::ContainsSubstring("standard ND DOF"));

  for (int standard = 0; standard < nd_fe.GetDof(); standard++)
  {
    for (int enrichment = 0; enrichment < 4; enrichment++)
    {
      const auto direct = DirectNDStandardEnrichmentIntegral(
          nd_fe, transformation, standard, element_dofs.nd[enrichment].basis, grad_lambda,
          jacobian_determinant);
      CAPTURE(standard, enrichment,
              matrices.nd_mass_standard_enrichment(standard, enrichment), direct[0],
              matrices.nd_mass_estimated_absolute_error(standard, enrichment),
              matrices.nd_curl_curl_standard_enrichment(standard, enrichment), direct[1],
              matrices.nd_curl_curl_estimated_absolute_error(standard, enrichment));
      CHECK(std::abs(matrices.nd_mass_standard_enrichment(standard, enrichment) -
                     direct[0]) <=
            2.0 * matrices.nd_mass_estimated_absolute_error(standard, enrichment) +
                2.0e-10);
      CHECK(std::abs(matrices.nd_curl_curl_standard_enrichment(standard, enrichment) -
                     direct[1]) <=
            2.0 * matrices.nd_curl_curl_estimated_absolute_error(standard, enrichment) +
                2.0e-10);
    }
  }

  for (int standard = 0; standard < h1_fe.GetDof(); standard++)
  {
    for (int enrichment = 0; enrichment < 2; enrichment++)
    {
      const double direct = DirectH1StandardEnrichmentIntegral(
          h1_fe, transformation, standard, element_dofs.h1[enrichment].basis, grad_lambda,
          jacobian_determinant);
      CAPTURE(standard, enrichment, matrices.h1_standard_enrichment(standard, enrichment),
              direct, matrices.h1_estimated_absolute_error(standard, enrichment));
      CHECK(std::abs(matrices.h1_standard_enrichment(standard, enrichment) - direct) <=
            2.0 * matrices.h1_estimated_absolute_error(standard, enrichment) + 2.0e-10);
    }
  }

  mfem::DenseMatrix discrete_gradient(nd_fe.GetDof(), h1_fe.GetDof());
  nd_fe.ProjectGrad(h1_fe, transformation, discrete_gradient);
  for (int standard_h1 = 0; standard_h1 < h1_fe.GetDof(); standard_h1++)
  {
    for (int enrichment = 0; enrichment < 2; enrichment++)
    {
      long double value = 0.0L;
      long double scale = 0.0L;
      for (int standard_nd = 0; standard_nd < nd_fe.GetDof(); standard_nd++)
      {
        const long double term =
            static_cast<long double>(discrete_gradient(standard_nd, standard_h1)) *
            static_cast<long double>(
                matrices.nd_mass_standard_enrichment(standard_nd, enrichment));
        value += term;
        scale += std::abs(term);
      }
      const long double residual =
          std::abs(static_cast<long double>(
                       matrices.h1_standard_enrichment(standard_h1, enrichment)) -
                   value);
      CHECK(residual <= 64.0L * std::numeric_limits<double>::epsilon() * scale);
    }
  }
  for (int standard = 0; standard < nd_fe.GetDof(); standard++)
  {
    for (int gradient = 0; gradient < 2; gradient++)
    {
      CHECK(matrices.nd_curl_curl_standard_enrichment(standard, gradient) == 0.0);
      CHECK(matrices.nd_curl_curl_estimated_absolute_error(standard, gradient) == 0.0);
    }
  }

  const auto rotate = [](const fem::singular::Vector3 &x)
  { return fem::singular::Vector3{x[2], x[0], x[1]}; };
  const fem::singular::Vector3 rotated_origin{-0.8, 0.6, 1.1};
  auto rotated_mesh =
      AffineTetrahedronMesh(rotated_origin, rotate(a), rotate(b), rotate(c));
  const auto rotated = fem::singular::AssembleElementStandardEnrichmentMatrices(
      element_dofs, h1_fe, nd_fe, *rotated_mesh.GetElementTransformation(0), options);
  for (int standard = 0; standard < nd_fe.GetDof(); standard++)
  {
    for (int enrichment = 0; enrichment < 4; enrichment++)
    {
      CHECK(std::abs(rotated.nd_mass_standard_enrichment(standard, enrichment) -
                     matrices.nd_mass_standard_enrichment(standard, enrichment)) <=
            rotated.nd_mass_estimated_absolute_error(standard, enrichment) +
                matrices.nd_mass_estimated_absolute_error(standard, enrichment) + 2.0e-11);
      CHECK(std::abs(rotated.nd_curl_curl_standard_enrichment(standard, enrichment) -
                     matrices.nd_curl_curl_standard_enrichment(standard, enrichment)) <=
            rotated.nd_curl_curl_estimated_absolute_error(standard, enrichment) +
                matrices.nd_curl_curl_estimated_absolute_error(standard, enrichment) +
                2.0e-11);
    }
  }

  mfem::DofTransformation h1_dof_transformation;
  mfem::ND_TetDofTransformation nd_stateless_transformation(2);
  mfem::DofTransformation nd_dof_transformation(nd_stateless_transformation);
  mfem::Array<int> face_orientations(4);
  face_orientations[0] = 1;
  face_orientations[1] = 2;
  face_orientations[2] = 4;
  face_orientations[3] = 5;
  nd_dof_transformation.SetFaceOrientations(face_orientations);

  auto oriented = matrices;
  fem::singular::ApplyStandardDofTransformations(h1_dof_transformation,
                                                 nd_dof_transformation, oriented);
  CheckExactTranspose(oriented.h1_standard_enrichment, oriented.h1_enrichment_standard);
  CheckExactTranspose(oriented.nd_mass_standard_enrichment,
                      oriented.nd_mass_enrichment_standard);
  CheckExactTranspose(oriented.nd_curl_curl_standard_enrichment,
                      oriented.nd_curl_curl_enrichment_standard);

  mfem::DenseMatrix expected_mass = matrices.nd_mass_standard_enrichment;
  mfem::DenseMatrix expected_curl_curl = matrices.nd_curl_curl_standard_enrichment;
  mfem::DenseMatrix dual_transformation(nd_fe.GetDof());
  dual_transformation = 0.0;
  for (int i = 0; i < nd_fe.GetDof(); i++)
  {
    dual_transformation(i, i) = 1.0;
  }
  nd_dof_transformation.TransformDualCols(dual_transformation);
  nd_dof_transformation.TransformDualCols(expected_mass);
  nd_dof_transformation.TransformDualCols(expected_curl_curl);
  for (int standard = 0; standard < nd_fe.GetDof(); standard++)
  {
    for (int enrichment = 0; enrichment < 4; enrichment++)
    {
      CHECK(oriented.nd_mass_standard_enrichment(standard, enrichment) ==
            expected_mass(standard, enrichment));
      CHECK(oriented.nd_curl_curl_standard_enrichment(standard, enrichment) ==
            expected_curl_curl(standard, enrichment));
      long double expected_mass_error = 0.0L;
      long double expected_curl_curl_error = 0.0L;
      for (int k = 0; k < nd_fe.GetDof(); k++)
      {
        const long double weight =
            std::abs(static_cast<long double>(dual_transformation(standard, k)));
        expected_mass_error +=
            weight * static_cast<long double>(
                         matrices.nd_mass_estimated_absolute_error(k, enrichment));
        expected_curl_curl_error +=
            weight * static_cast<long double>(
                         matrices.nd_curl_curl_estimated_absolute_error(k, enrichment));
      }
      CHECK(oriented.nd_mass_estimated_absolute_error(standard, enrichment) ==
            static_cast<double>(expected_mass_error));
      CHECK(oriented.nd_curl_curl_estimated_absolute_error(standard, enrichment) ==
            static_cast<double>(expected_curl_curl_error));
    }
  }

  mfem::DenseMatrix oriented_gradient = discrete_gradient;
  mfem::TransformPrimal(nd_dof_transformation, h1_dof_transformation, oriented_gradient);
  for (int standard_h1 = 0; standard_h1 < h1_fe.GetDof(); standard_h1++)
  {
    for (int enrichment = 0; enrichment < 2; enrichment++)
    {
      long double value = 0.0L;
      long double scale = 0.0L;
      for (int standard_nd = 0; standard_nd < nd_fe.GetDof(); standard_nd++)
      {
        const long double term =
            static_cast<long double>(oriented_gradient(standard_nd, standard_h1)) *
            static_cast<long double>(
                oriented.nd_mass_standard_enrichment(standard_nd, enrichment));
        value += term;
        scale += std::abs(term);
      }
      const long double residual =
          std::abs(static_cast<long double>(
                       oriented.h1_standard_enrichment(standard_h1, enrichment)) -
                   value);
      CHECK(residual <= 64.0L * std::numeric_limits<double>::epsilon() * scale);
    }
  }

  mfem::ND_TetDofTransformation wrong_size_stateless_transformation(3);
  mfem::DofTransformation wrong_size_transformation(wrong_size_stateless_transformation);
  wrong_size_transformation.SetFaceOrientations(face_orientations);
  auto wrong_size_matrices = matrices;
  CHECK_THROWS_AS(fem::singular::ApplyStandardDofTransformations(h1_dof_transformation,
                                                                 wrong_size_transformation,
                                                                 wrong_size_matrices),
                  std::invalid_argument);
  CheckMatricesEqual(wrong_size_matrices.h1_standard_enrichment,
                     matrices.h1_standard_enrichment);
  CheckMatricesEqual(wrong_size_matrices.h1_enrichment_standard,
                     matrices.h1_enrichment_standard);
  CheckMatricesEqual(wrong_size_matrices.h1_estimated_absolute_error,
                     matrices.h1_estimated_absolute_error);
  CheckMatricesEqual(wrong_size_matrices.nd_mass_standard_enrichment,
                     matrices.nd_mass_standard_enrichment);
  CheckMatricesEqual(wrong_size_matrices.nd_mass_enrichment_standard,
                     matrices.nd_mass_enrichment_standard);
  CheckMatricesEqual(wrong_size_matrices.nd_mass_estimated_absolute_error,
                     matrices.nd_mass_estimated_absolute_error);
  CheckMatricesEqual(wrong_size_matrices.nd_curl_curl_standard_enrichment,
                     matrices.nd_curl_curl_standard_enrichment);
  CheckMatricesEqual(wrong_size_matrices.nd_curl_curl_enrichment_standard,
                     matrices.nd_curl_curl_enrichment_standard);
  CheckMatricesEqual(wrong_size_matrices.nd_curl_curl_estimated_absolute_error,
                     matrices.nd_curl_curl_estimated_absolute_error);

  mfem::H1_TetrahedronElement wrong_order_h1(1);
  CHECK_THROWS_AS(fem::singular::AssembleElementStandardEnrichmentMatrices(
                      element_dofs, wrong_order_h1, nd_fe, transformation, options),
                  std::invalid_argument);
  const fem::singular::AdaptiveAssemblyOptions impossible{8, 1.0e-30, 1.0e-30, 1};
  CHECK_THROWS_AS(fem::singular::AssembleElementStandardEnrichmentMatrices(
                      element_dofs, h1_fe, nd_fe, transformation, impossible),
                  std::runtime_error);
}

TEST_CASE("Adaptive coupling uses its full global error budget",
          "[singularelements][singularassembly][Serial]")
{
  constexpr double scale = 1.0 / 640.0;
  const fem::singular::Vector3 origin{-160.0 * scale, -1.5 * scale, 0.0};
  const fem::singular::Vector3 a{0.0, 1.5 * scale, 0.0};
  const fem::singular::Vector3 b{0.0, 1.5 * scale, (10.0 / 3.0) * scale};
  const fem::singular::Vector3 c{0.36596077291304 * scale, 1.5 * scale,
                                 (10.0 / 3.0) * scale};
  auto mesh = AffineTetrahedronMesh(origin, a, b, c);
  auto &transformation = *mesh.GetElementTransformation(0);
  mfem::H1_TetrahedronElement h1_fe(1);
  mfem::ND_TetrahedronElement nd_fe(1);

  fem::singular::ElementDofMap element_dofs;
  element_dofs.nd.push_back({0,
                             {fem::singular::HigherOrderBasisFamily::EDGE_ROTATIONAL,
                              {1, 2, 0, 3},
                              {1, 1, 1, 1},
                              1,
                              0.5}});
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-3, 1.0e-3, 6};
  const auto matrices = fem::singular::AssembleElementStandardEnrichmentMatrices(
      element_dofs, h1_fe, nd_fe, transformation, options);

  CHECK(matrices.maximum_subdivision_depth > 2);
  CHECK(matrices.maximum_subdivision_depth <= options.maximum_subdivisions);
  CHECK(matrices.total_quadrature_leaf_count > 0);
}

TEST_CASE("Isotropic materials and affine dilation scale singular element blocks",
          "[singularelements][singularassembly][Serial]")
{
  constexpr double nu = 0.5;
  const std::array<int, 4> canonical_nodes{0, 1, 2, 3};
  const auto gradients =
      fem::singular::EnumerateHigherOrderNodeGradientBases(canonical_nodes, 1, nu);
  const auto rotations =
      fem::singular::EnumerateHigherOrderNodeRotationalBases(canonical_nodes, 1, nu);
  REQUIRE(gradients.size() >= 2);
  REQUIRE(rotations.size() >= 2);

  fem::singular::ElementDofMap element_dofs;
  for (int i = 0; i < 2; i++)
  {
    element_dofs.h1.push_back({static_cast<std::size_t>(i), gradients[i]});
    element_dofs.nd.push_back({static_cast<std::size_t>(i), gradients[i]});
  }
  for (int i = 0; i < 2; i++)
  {
    element_dofs.nd.push_back({static_cast<std::size_t>(i + 2), rotations[i]});
  }

  const fem::singular::Vector3 origin{0.4, -0.3, 0.2};
  const fem::singular::Vector3 a{1.7, 0.2, -0.1};
  const fem::singular::Vector3 b{0.3, 1.4, 0.25};
  const fem::singular::Vector3 c{-0.2, 0.1, 0.9};
  double jacobian_determinant;
  const auto grad_lambda = AffineBarycentricGradients(a, b, c, jacobian_determinant);
  auto mesh = AffineTetrahedronMesh(origin, a, b, c);
  mfem::H1_TetrahedronElement h1_fe(2);
  mfem::ND_TetrahedronElement nd_fe(2);
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto enrichment = fem::singular::AssembleElementEnrichmentMatrices(
      element_dofs, grad_lambda, jacobian_determinant, options);
  const auto coupling = fem::singular::AssembleElementStandardEnrichmentMatrices(
      element_dofs, h1_fe, nd_fe, *mesh.GetElementTransformation(0), options);

  constexpr fem::singular::IsotropicMaterialCoefficients material{3.7, 0.4};
  auto material_enrichment = enrichment;
  auto material_coupling = coupling;
  fem::singular::ApplyIsotropicMaterialCoefficients(material, material_enrichment);
  fem::singular::ApplyIsotropicMaterialCoefficients(material, material_coupling);

  CHECK(material_enrichment.total_quadrature_leaf_count ==
        enrichment.total_quadrature_leaf_count);
  CHECK(material_enrichment.maximum_subdivision_depth ==
        enrichment.maximum_subdivision_depth);
  CheckScaledMatrix(material_enrichment.h1_diffusion, enrichment.h1_diffusion,
                    material.electric);
  CheckScaledMatrix(material_enrichment.h1_diffusion_estimated_absolute_error,
                    enrichment.h1_diffusion_estimated_absolute_error, material.electric);
  CheckScaledMatrix(material_enrichment.nd_mass, enrichment.nd_mass, material.electric);
  CheckScaledMatrix(material_enrichment.nd_mass_estimated_absolute_error,
                    enrichment.nd_mass_estimated_absolute_error, material.electric);
  CheckScaledMatrix(material_enrichment.nd_curl_curl, enrichment.nd_curl_curl,
                    material.inverse_magnetic);
  CheckScaledMatrix(material_enrichment.nd_curl_curl_estimated_absolute_error,
                    enrichment.nd_curl_curl_estimated_absolute_error,
                    material.inverse_magnetic);

  CheckScaledMatrix(material_coupling.h1_standard_enrichment,
                    coupling.h1_standard_enrichment, material.electric);
  CheckScaledMatrix(material_coupling.h1_enrichment_standard,
                    coupling.h1_enrichment_standard, material.electric);
  CheckScaledMatrix(material_coupling.h1_estimated_absolute_error,
                    coupling.h1_estimated_absolute_error, material.electric);
  CheckScaledMatrix(material_coupling.nd_mass_standard_enrichment,
                    coupling.nd_mass_standard_enrichment, material.electric);
  CheckScaledMatrix(material_coupling.nd_mass_enrichment_standard,
                    coupling.nd_mass_enrichment_standard, material.electric);
  CheckScaledMatrix(material_coupling.nd_mass_estimated_absolute_error,
                    coupling.nd_mass_estimated_absolute_error, material.electric);
  CheckScaledMatrix(material_coupling.nd_curl_curl_standard_enrichment,
                    coupling.nd_curl_curl_standard_enrichment, material.inverse_magnetic);
  CheckScaledMatrix(material_coupling.nd_curl_curl_enrichment_standard,
                    coupling.nd_curl_curl_enrichment_standard, material.inverse_magnetic);
  CheckScaledMatrix(material_coupling.nd_curl_curl_estimated_absolute_error,
                    coupling.nd_curl_curl_estimated_absolute_error,
                    material.inverse_magnetic);
  CheckExactTranspose(material_coupling.h1_standard_enrichment,
                      material_coupling.h1_enrichment_standard);
  CheckExactTranspose(material_coupling.nd_mass_standard_enrichment,
                      material_coupling.nd_mass_enrichment_standard);
  CheckExactTranspose(material_coupling.nd_curl_curl_standard_enrichment,
                      material_coupling.nd_curl_curl_enrichment_standard);
  for (int i = 0; i < material_enrichment.h1_diffusion.Height(); i++)
  {
    for (int j = 0; j < material_enrichment.h1_diffusion.Width(); j++)
    {
      CHECK(material_enrichment.h1_diffusion(i, j) == material_enrichment.nd_mass(i, j));
      CHECK(material_enrichment.h1_diffusion_estimated_absolute_error(i, j) ==
            material_enrichment.nd_mass_estimated_absolute_error(i, j));
    }
  }

  for (const auto invalid : {fem::singular::IsotropicMaterialCoefficients{0.0, 1.0},
                             fem::singular::IsotropicMaterialCoefficients{-1.0, 1.0},
                             fem::singular::IsotropicMaterialCoefficients{
                                 std::numeric_limits<double>::quiet_NaN(), 1.0},
                             fem::singular::IsotropicMaterialCoefficients{
                                 1.0, std::numeric_limits<double>::infinity()}})
  {
    auto rejected = enrichment;
    CHECK_THROWS_AS(fem::singular::ApplyIsotropicMaterialCoefficients(invalid, rejected),
                    std::invalid_argument);
    CheckMatricesEqual(rejected.h1_diffusion, enrichment.h1_diffusion);
    CheckMatricesEqual(rejected.nd_mass, enrichment.nd_mass);
    CheckMatricesEqual(rejected.nd_curl_curl, enrichment.nd_curl_curl);
  }

  auto malformed = coupling;
  malformed.nd_curl_curl_enrichment_standard(0, 0) += 1.0;
  CHECK_THROWS_AS(fem::singular::ApplyIsotropicMaterialCoefficients(material, malformed),
                  std::invalid_argument);
  CheckMatricesEqual(malformed.h1_standard_enrichment, coupling.h1_standard_enrichment);
  CheckMatricesEqual(malformed.nd_mass_standard_enrichment,
                     coupling.nd_mass_standard_enrichment);

  auto overflow = enrichment;
  overflow.nd_curl_curl(2, 2) = std::numeric_limits<double>::max();
  CHECK_THROWS_AS(fem::singular::ApplyIsotropicMaterialCoefficients(
                      fem::singular::IsotropicMaterialCoefficients{1.0, 2.0}, overflow),
                  std::overflow_error);
  CheckMatricesEqual(overflow.h1_diffusion, enrichment.h1_diffusion);
  CHECK(overflow.nd_curl_curl(2, 2) == std::numeric_limits<double>::max());

  constexpr double length_scale = 2.5;
  const auto scale = [](const fem::singular::Vector3 &value)
  {
    fem::singular::Vector3 result = value;
    for (double &entry : result)
    {
      entry *= length_scale;
    }
    return result;
  };
  double scaled_jacobian_determinant;
  const auto scaled_grad_lambda =
      AffineBarycentricGradients(scale(a), scale(b), scale(c), scaled_jacobian_determinant);
  auto scaled_mesh = AffineTetrahedronMesh(origin, scale(a), scale(b), scale(c));
  const auto scaled_enrichment = fem::singular::AssembleElementEnrichmentMatrices(
      element_dofs, scaled_grad_lambda, scaled_jacobian_determinant, options);
  const auto scaled_coupling = fem::singular::AssembleElementStandardEnrichmentMatrices(
      element_dofs, h1_fe, nd_fe, *scaled_mesh.GetElementTransformation(0), options);

  CheckAffineScaling(scaled_enrichment.h1_diffusion,
                     scaled_enrichment.h1_diffusion_estimated_absolute_error,
                     enrichment.h1_diffusion,
                     enrichment.h1_diffusion_estimated_absolute_error, length_scale);
  CheckAffineScaling(scaled_enrichment.nd_mass,
                     scaled_enrichment.nd_mass_estimated_absolute_error, enrichment.nd_mass,
                     enrichment.nd_mass_estimated_absolute_error, length_scale);
  CheckAffineScaling(scaled_enrichment.nd_curl_curl,
                     scaled_enrichment.nd_curl_curl_estimated_absolute_error,
                     enrichment.nd_curl_curl,
                     enrichment.nd_curl_curl_estimated_absolute_error, 1.0 / length_scale);
  CheckAffineScaling(
      scaled_coupling.h1_standard_enrichment, scaled_coupling.h1_estimated_absolute_error,
      coupling.h1_standard_enrichment, coupling.h1_estimated_absolute_error, length_scale);
  CheckAffineScaling(scaled_coupling.nd_mass_standard_enrichment,
                     scaled_coupling.nd_mass_estimated_absolute_error,
                     coupling.nd_mass_standard_enrichment,
                     coupling.nd_mass_estimated_absolute_error, length_scale);
  CheckAffineScaling(scaled_coupling.nd_curl_curl_standard_enrichment,
                     scaled_coupling.nd_curl_curl_estimated_absolute_error,
                     coupling.nd_curl_curl_standard_enrichment,
                     coupling.nd_curl_curl_estimated_absolute_error, 1.0 / length_scale);
}

TEST_CASE("Local sparse singular assembly preserves element bilinear forms",
          "[singularelements][singularassembly][Serial]")
{
  auto mesh = SharedFaceTetrahedronMesh(false);
  mfem::H1_FECollection h1_collection(2, 3);
  mfem::ND_FECollection nd_collection(2, 3);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);

  fem::singular::DofTopology topology;
  topology.h1_dofs.resize(1);
  topology.h1_dofs[0].family = fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
  topology.nd_dofs.resize(2);
  topology.nd_dofs[0] = topology.h1_dofs[0];
  topology.nd_dofs[1].family = fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
  topology.h1_to_nd = {0};
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    topology.elements.push_back(SharedFaceSingularDofs(*mesh.GetElement(element)));
  }

  const std::vector<fem::singular::IsotropicMaterialCoefficients> materials{{2.3, 0.8},
                                                                            {4.1, 1.7}};
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto sparse = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      topology, h1_space, nd_space, materials, options);

  REQUIRE(sparse.h1_diffusion.enrichment_enrichment);
  REQUIRE(sparse.h1_diffusion.standard_enrichment);
  REQUIRE(sparse.h1_diffusion.enrichment_standard);
  REQUIRE(sparse.h1_mass.enrichment_enrichment);
  REQUIRE(sparse.h1_mass.standard_enrichment);
  REQUIRE(sparse.h1_mass.enrichment_standard);
  REQUIRE(sparse.nd_mass.enrichment_enrichment);
  REQUIRE(sparse.nd_mass.standard_enrichment);
  REQUIRE(sparse.nd_mass.enrichment_standard);
  REQUIRE(sparse.nd_curl_curl.enrichment_enrichment);
  REQUIRE(sparse.nd_curl_curl.standard_enrichment);
  REQUIRE(sparse.nd_curl_curl.enrichment_standard);
  CHECK(sparse.total_quadrature_leaf_count > 0);
  CHECK(sparse.maximum_subdivision_depth <= options.maximum_subdivisions);
  CHECK(sparse.h1_diffusion.enrichment_enrichment->Height() == 1);
  CHECK(sparse.h1_diffusion.standard_enrichment->Height() == h1_space.GetVSize());
  CHECK(sparse.h1_mass.enrichment_enrichment->Height() == 1);
  CHECK(sparse.h1_mass.standard_enrichment->Height() == h1_space.GetVSize());
  CHECK(sparse.nd_mass.enrichment_enrichment->Height() == 2);
  CHECK(sparse.nd_mass.standard_enrichment->Height() == nd_space.GetVSize());
  CheckExactSparseTranspose(*sparse.h1_diffusion.standard_enrichment,
                            *sparse.h1_diffusion.enrichment_standard);
  CheckExactSparseTranspose(*sparse.h1_mass.standard_enrichment,
                            *sparse.h1_mass.enrichment_standard);
  CheckExactSparseTranspose(*sparse.nd_mass.standard_enrichment,
                            *sparse.nd_mass.enrichment_standard);
  CheckExactSparseTranspose(*sparse.nd_curl_curl.standard_enrichment,
                            *sparse.nd_curl_curl.enrichment_standard);
  CheckSparseNonnegative(
      *sparse.h1_diffusion.enrichment_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.h1_diffusion.standard_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.h1_mass.enrichment_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.h1_mass.standard_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.nd_mass.enrichment_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.nd_mass.standard_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(
      *sparse.nd_curl_curl.enrichment_enrichment_estimated_absolute_error);
  CheckSparseNonnegative(*sparse.nd_curl_curl.standard_enrichment_estimated_absolute_error);

  mfem::FunctionCoefficient h1_coefficient(
      [](const mfem::Vector &x) { return x[0] + 2.0 * x[1] - 0.5 * x[2] + x[0] * x[1]; });
  mfem::VectorFunctionCoefficient nd_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value[0] = x[1] + 2.0 * x[2];
        value[1] = -x[0] + 0.5 * x[2];
        value[2] = 3.0 * x[0] - x[1];
      });
  mfem::GridFunction h1_field(&h1_space);
  mfem::GridFunction nd_field(&nd_space);
  h1_field.ProjectCoefficient(h1_coefficient);
  nd_field.ProjectCoefficient(nd_coefficient);
  mfem::Vector h1_enrichment(1);
  h1_enrichment[0] = 0.37;
  mfem::Vector nd_enrichment(2);
  nd_enrichment[0] = 0.37;
  nd_enrichment[1] = -0.21;

  struct Energies
  {
    double h1_coupling = 0.0;
    double h1_enrichment = 0.0;
    double nd_mass_coupling = 0.0;
    double nd_mass_enrichment = 0.0;
    double nd_curl_curl_coupling = 0.0;
    double nd_curl_curl_enrichment = 0.0;
    double h1_coupling_error = 0.0;
    double h1_enrichment_error = 0.0;
    double nd_mass_coupling_error = 0.0;
    double nd_mass_enrichment_error = 0.0;
    double nd_curl_curl_coupling_error = 0.0;
    double nd_curl_curl_enrichment_error = 0.0;
  } element_sum;

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &element_dofs = topology.elements[element];
    auto &transformation = *mesh.GetElementTransformation(element);
    double jacobian_determinant;
    const auto grad_lambda =
        PhysicalBarycentricGradients(transformation, jacobian_determinant);
    auto enrichment = fem::singular::AssembleElementEnrichmentMatrices(
        element_dofs, grad_lambda, jacobian_determinant, options);
    auto coupling = fem::singular::AssembleElementStandardEnrichmentMatrices(
        element_dofs, *h1_space.GetFE(element), *nd_space.GetFE(element), transformation,
        options);
    fem::singular::ApplyIsotropicMaterialCoefficients(materials[element], enrichment);
    fem::singular::ApplyIsotropicMaterialCoefficients(materials[element], coupling);

    mfem::Array<int> h1_vdofs, nd_vdofs;
    mfem::DofTransformation h1_dof_transformation, nd_dof_transformation;
    h1_space.GetElementVDofs(element, h1_vdofs, h1_dof_transformation);
    nd_space.GetElementVDofs(element, nd_vdofs, nd_dof_transformation);
    mfem::Vector h1_local, nd_local;
    h1_field.GetSubVector(h1_vdofs, h1_local);
    nd_field.GetSubVector(nd_vdofs, nd_local);
    h1_dof_transformation.InvTransformPrimal(h1_local);
    nd_dof_transformation.InvTransformPrimal(nd_local);

    mfem::Vector result(h1_local.Size());
    coupling.h1_standard_enrichment.Mult(h1_enrichment, result);
    element_sum.h1_coupling += h1_local * result;
    result.SetSize(h1_enrichment.Size());
    enrichment.h1_diffusion.Mult(h1_enrichment, result);
    element_sum.h1_enrichment += h1_enrichment * result;
    element_sum.h1_coupling_error +=
        BilinearErrorBound(h1_local, coupling.h1_estimated_absolute_error, h1_enrichment);
    element_sum.h1_enrichment_error += BilinearErrorBound(
        h1_enrichment, enrichment.h1_diffusion_estimated_absolute_error, h1_enrichment);

    result.SetSize(nd_local.Size());
    coupling.nd_mass_standard_enrichment.Mult(nd_enrichment, result);
    element_sum.nd_mass_coupling += nd_local * result;
    coupling.nd_curl_curl_standard_enrichment.Mult(nd_enrichment, result);
    element_sum.nd_curl_curl_coupling += nd_local * result;
    result.SetSize(nd_enrichment.Size());
    enrichment.nd_mass.Mult(nd_enrichment, result);
    element_sum.nd_mass_enrichment += nd_enrichment * result;
    enrichment.nd_curl_curl.Mult(nd_enrichment, result);
    element_sum.nd_curl_curl_enrichment += nd_enrichment * result;
    element_sum.nd_mass_coupling_error += BilinearErrorBound(
        nd_local, coupling.nd_mass_estimated_absolute_error, nd_enrichment);
    element_sum.nd_mass_enrichment_error += BilinearErrorBound(
        nd_enrichment, enrichment.nd_mass_estimated_absolute_error, nd_enrichment);
    element_sum.nd_curl_curl_coupling_error += BilinearErrorBound(
        nd_local, coupling.nd_curl_curl_estimated_absolute_error, nd_enrichment);
    element_sum.nd_curl_curl_enrichment_error += BilinearErrorBound(
        nd_enrichment, enrichment.nd_curl_curl_estimated_absolute_error, nd_enrichment);
  }

  const auto h1_absolute = AbsoluteValue(h1_field);
  const auto nd_absolute = AbsoluteValue(nd_field);
  const auto h1_enrichment_absolute = AbsoluteValue(h1_enrichment);
  const auto nd_enrichment_absolute = AbsoluteValue(nd_enrichment);
  const auto check_energy = [](double sparse_value, double element_value)
  {
    CHECK(std::abs(sparse_value - element_value) <=
          128.0 * std::numeric_limits<double>::epsilon() *
              std::max({1.0, std::abs(sparse_value), std::abs(element_value)}));
  };
  const auto check_bound = [](double sparse_bound, double element_bound)
  {
    CHECK(sparse_bound + 128.0 * std::numeric_limits<double>::epsilon() *
                             std::max({1.0, sparse_bound, element_bound}) >=
          element_bound);
  };
  check_energy(
      SparseBilinear(h1_field, *sparse.h1_diffusion.standard_enrichment, h1_enrichment),
      element_sum.h1_coupling);
  check_energy(SparseBilinear(h1_enrichment, *sparse.h1_diffusion.enrichment_enrichment,
                              h1_enrichment),
               element_sum.h1_enrichment);
  check_energy(SparseBilinear(nd_field, *sparse.nd_mass.standard_enrichment, nd_enrichment),
               element_sum.nd_mass_coupling);
  check_energy(
      SparseBilinear(nd_enrichment, *sparse.nd_mass.enrichment_enrichment, nd_enrichment),
      element_sum.nd_mass_enrichment);
  check_energy(
      SparseBilinear(nd_field, *sparse.nd_curl_curl.standard_enrichment, nd_enrichment),
      element_sum.nd_curl_curl_coupling);
  check_energy(SparseBilinear(nd_enrichment, *sparse.nd_curl_curl.enrichment_enrichment,
                              nd_enrichment),
               element_sum.nd_curl_curl_enrichment);

  check_bound(
      SparseBilinear(h1_absolute,
                     *sparse.h1_diffusion.standard_enrichment_estimated_absolute_error,
                     h1_enrichment_absolute),
      element_sum.h1_coupling_error);
  check_energy(
      SparseBilinear(h1_enrichment_absolute,
                     *sparse.h1_diffusion.enrichment_enrichment_estimated_absolute_error,
                     h1_enrichment_absolute),
      element_sum.h1_enrichment_error);
  check_bound(SparseBilinear(nd_absolute,
                             *sparse.nd_mass.standard_enrichment_estimated_absolute_error,
                             nd_enrichment_absolute),
              element_sum.nd_mass_coupling_error);
  check_energy(
      SparseBilinear(nd_enrichment_absolute,
                     *sparse.nd_mass.enrichment_enrichment_estimated_absolute_error,
                     nd_enrichment_absolute),
      element_sum.nd_mass_enrichment_error);
  check_bound(
      SparseBilinear(nd_absolute,
                     *sparse.nd_curl_curl.standard_enrichment_estimated_absolute_error,
                     nd_enrichment_absolute),
      element_sum.nd_curl_curl_coupling_error);
  check_energy(
      SparseBilinear(nd_enrichment_absolute,
                     *sparse.nd_curl_curl.enrichment_enrichment_estimated_absolute_error,
                     nd_enrichment_absolute),
      element_sum.nd_curl_curl_enrichment_error);

  auto invalid_topology = topology;
  invalid_topology.elements[0].h1[0].dof = 1;
  CHECK_THROWS_AS(fem::singular::AssembleLocalSparseEnrichmentMatrices(
                      invalid_topology, h1_space, nd_space, materials, options),
                  std::invalid_argument);
  auto missing_material = materials;
  missing_material.pop_back();
  CHECK_THROWS_AS(fem::singular::AssembleLocalSparseEnrichmentMatrices(
                      topology, h1_space, nd_space, missing_material, options),
                  std::invalid_argument);
}

TEST_CASE("MFEM standard orientations preserve assembled singular coupling",
          "[singularelements][singularassembly][Serial]")
{
  const auto canonical = AssembleSharedFaceCoupling(false);
  const auto permuted = AssembleSharedFaceCoupling(true);
  CHECK(std::abs(canonical.h1 - permuted.h1) <=
        canonical.h1_error + permuted.h1_error + 2.0e-12);
  CHECK(std::abs(canonical.nd_mass - permuted.nd_mass) <=
        canonical.nd_mass_error + permuted.nd_mass_error + 2.0e-12);
  CHECK(std::abs(canonical.nd_curl_curl - permuted.nd_curl_curl) <=
        canonical.nd_curl_curl_error + permuted.nd_curl_curl_error + 2.0e-12);
}

TEST_CASE("Additive overlapping patch correction is symmetric positive definite",
          "[singularelements][singularassembly][Serial]")
{
  mfem::DenseMatrix full(3);
  full = 0.0;
  full(0, 0) = 4.0;
  full(0, 1) = full(1, 0) = -1.0;
  full(1, 1) = 5.0;
  full(1, 2) = full(2, 1) = -1.5;
  full(2, 2) = 3.0;

  mfem::DenseMatrix first(2), second(2);
  first(0, 0) = full(0, 0);
  first(0, 1) = first(1, 0) = full(0, 1);
  first(1, 1) = full(1, 1);
  second(0, 0) = full(1, 1);
  second(0, 1) = second(1, 0) = full(1, 2);
  second(1, 1) = full(2, 2);

  mfem::Array<int> first_dofs, second_dofs;
  first_dofs.Append(0);
  first_dofs.Append(1);
  second_dofs.Append(1);
  second_dofs.Append(2);
  std::vector<mfem::Array<int>> patch_dofs;
  patch_dofs.push_back(first_dofs);
  patch_dofs.push_back(second_dofs);
  std::vector<std::unique_ptr<mfem::Solver>> patch_solvers;
  for (int patch = 0; patch < 2; patch++)
  {
    auto solver = std::make_unique<mfem::CGSolver>(Mpi::World());
    solver->SetRelTol(1.0e-14);
    solver->SetAbsTol(1.0e-15);
    solver->SetMaxIter(20);
    solver->SetPrintLevel(0);
    patch_solvers.push_back(std::move(solver));
  }
  AdditivePatchSolver additive(3, std::move(patch_dofs), std::move(patch_solvers));
  additive.SetPatchOperators({&first, &second});
  additive.SetOperator(full);

  mfem::Vector left({0.3, -0.7, 1.1});
  mfem::Vector right({-0.2, 0.9, 0.4});
  mfem::Vector corrected_left, corrected_right;
  additive.Mult(left, corrected_left);
  additive.Mult(right, corrected_right);
  CheckClose(left * corrected_right, right * corrected_left);
  CHECK(left * corrected_left > 0.0);
  CHECK(right * corrected_right > 0.0);

  for (int coordinate = 0; coordinate < 3; coordinate++)
  {
    mfem::Vector basis(3), corrected;
    basis = 0.0;
    basis[coordinate] = 1.0;
    additive.Mult(basis, corrected);
    CAPTURE(coordinate);
    CHECK(basis * corrected > 0.0);
  }
}

TEST_CASE("Parallel triangular sparse assembly preserves global operators",
          "[singularelements][singularassembly][triangle][Parallel]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  auto serial_mesh = InternalLineTipTriangleMesh();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  const auto serial_topology =
      fem::singular::BuildSerialTriangleDofTopology(serial_mesh, serial_features, 1);
  const std::array<int, 8> partition{0, 0, 0, 0, 1, 1, 1, 1};
  mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh, partition.data());
  const auto parallel_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
      serial_mesh, parallel_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids;
  for (int element = 0; element < serial_mesh.GetNE(); element++)
  {
    if (partition[element] == Mpi::Rank(Mpi::World()))
    {
      source_element_ids.push_back(element);
    }
  }
  const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, parallel_mesh, parallel_vertex_ids, source_element_ids);
  const auto local_topology = fem::singular::BuildLocalTriangleDofTopology(
      parallel_mesh, local_features, parallel_vertex_ids, 1);
  const auto numbering =
      fem::singular::BuildParallelDofNumbering(Mpi::World(), local_topology);

  mfem::H1_FECollection serial_h1_collection(1, 2);
  mfem::ND_FECollection serial_nd_collection(1, 2);
  mfem::FiniteElementSpace serial_h1_space(&serial_mesh, &serial_h1_collection);
  mfem::FiniteElementSpace serial_nd_space(&serial_mesh, &serial_nd_collection);
  mfem::H1_FECollection parallel_h1_collection(1, 2);
  mfem::ND_FECollection parallel_nd_collection(1, 2);
  mfem::ParFiniteElementSpace parallel_h1_space(&parallel_mesh, &parallel_h1_collection);
  mfem::ParFiniteElementSpace parallel_nd_space(&parallel_mesh, &parallel_nd_collection);

  const std::vector<fem::singular::IsotropicMaterialCoefficients> serial_materials(
      serial_mesh.GetNE(), {2.3, 0.8});
  const std::vector<fem::singular::IsotropicMaterialCoefficients> local_materials(
      parallel_mesh.GetNE(), {2.3, 0.8});
  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto serial = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      serial_topology, serial_h1_space, serial_nd_space, serial_materials, options);
  const auto local = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      local_topology, parallel_h1_space, parallel_nd_space, local_materials, options);
  const auto parallel = fem::singular::AssembleParallelSparseEnrichmentMatrices(
      local, numbering, parallel_h1_space, parallel_nd_space);
  const auto enrichment_gradient =
      fem::singular::BuildParallelEnrichmentGradient(Mpi::World(), numbering);

  CHECK(numbering.h1.global_size == 6);
  CHECK(numbering.nd.global_size == 12);
  bool has_signed_nd_dof = false;
  for (int element = 0; element < parallel_nd_space.GetNE(); element++)
  {
    mfem::Array<int> vdofs;
    parallel_nd_space.GetElementVDofs(element, vdofs);
    has_signed_nd_dof = has_signed_nd_dof || std::any_of(vdofs.begin(), vdofs.end(),
                                                         [](int dof) { return dof < 0; });
  }
  Mpi::GlobalOr(1, &has_signed_nd_dof, Mpi::World());
  CHECK(has_signed_nd_dof);

  const auto key_coefficient = [](const fem::singular::DofKey &key)
  {
    double value = 0.17 * (1 + static_cast<int>(key.family));
    for (std::size_t i = 0; i < key.support_entity.size; i++)
    {
      value += 0.031 * (i + 1) * (key.support_entity.vertices[i] + 1);
    }
    return value;
  };
  const auto serial_coefficients = [&](const std::vector<fem::singular::DofKey> &keys)
  {
    mfem::Vector result(static_cast<int>(keys.size()));
    for (int i = 0; i < result.Size(); i++)
    {
      result[i] = key_coefficient(keys[static_cast<std::size_t>(i)]);
    }
    return result;
  };
  const auto owned_coefficients = [&](const std::vector<fem::singular::DofKey> &keys,
                                      const fem::singular::TrueDofMap &map)
  {
    REQUIRE(map.owned_size <= std::numeric_limits<int>::max());
    mfem::Vector result(static_cast<int>(map.owned_size));
    result = 0.0;
    for (std::size_t local = 0; local < keys.size(); local++)
    {
      if (map.owner[local] != Mpi::Rank(Mpi::World()))
      {
        continue;
      }
      const HYPRE_BigInt owned = map.local_to_true[local] - map.owned_offset;
      REQUIRE(owned >= 0);
      REQUIRE(owned < map.owned_size);
      result[static_cast<int>(owned)] = key_coefficient(keys[local]);
    }
    return result;
  };
  const auto serial_h1 = serial_coefficients(serial_topology.h1_dofs);
  const auto serial_nd = serial_coefficients(serial_topology.nd_dofs);
  const auto parallel_h1 = owned_coefficients(local_topology.h1_dofs, numbering.h1);
  const auto parallel_nd = owned_coefficients(local_topology.nd_dofs, numbering.nd);
  CheckClose(
      SparseBilinear(serial_h1, *serial.h1_diffusion.enrichment_enrichment, serial_h1),
      ParallelBilinear(parallel_h1, *parallel.h1_diffusion.enrichment_enrichment,
                       parallel_h1));
  CheckClose(
      SparseBilinear(serial_h1, *serial.h1_mass.enrichment_enrichment, serial_h1),
      ParallelBilinear(parallel_h1, *parallel.h1_mass.enrichment_enrichment, parallel_h1));
  CheckClose(
      SparseBilinear(serial_nd, *serial.nd_mass.enrichment_enrichment, serial_nd),
      ParallelBilinear(parallel_nd, *parallel.nd_mass.enrichment_enrichment, parallel_nd));
  CheckClose(
      SparseBilinear(serial_nd, *serial.nd_curl_curl.enrichment_enrichment, serial_nd),
      ParallelBilinear(parallel_nd, *parallel.nd_curl_curl.enrichment_enrichment,
                       parallel_nd));

  mfem::ParDiscreteLinearOperator standard_gradient(&parallel_h1_space, &parallel_nd_space);
  standard_gradient.AddDomainInterpolator(new mfem::GradientInterpolator);
  standard_gradient.Assemble();
  standard_gradient.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_gradient_matrix(
      standard_gradient.ParallelAssemble());
  REQUIRE(standard_gradient_matrix);

  mfem::ConstantCoefficient electric(2.3);
  mfem::ConstantCoefficient inverse_magnetic(0.8);
  mfem::ParBilinearForm standard_h1_form(&parallel_h1_space);
  standard_h1_form.AddDomainIntegrator(new mfem::DiffusionIntegrator(electric));
  standard_h1_form.Assemble();
  standard_h1_form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_h1_matrix(
      standard_h1_form.ParallelAssemble());
  mfem::ParBilinearForm standard_nd_mass_form(&parallel_nd_space);
  standard_nd_mass_form.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(electric));
  standard_nd_mass_form.Assemble();
  standard_nd_mass_form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_nd_mass_matrix(
      standard_nd_mass_form.ParallelAssemble());
  mfem::ParBilinearForm standard_nd_curl_form(&parallel_nd_space);
  standard_nd_curl_form.AddDomainIntegrator(new mfem::CurlCurlIntegrator(inverse_magnetic));
  standard_nd_curl_form.Assemble();
  standard_nd_curl_form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_nd_curl_matrix(
      standard_nd_curl_form.ParallelAssemble());
  REQUIRE(standard_h1_matrix);
  REQUIRE(standard_nd_mass_matrix);
  REQUIRE(standard_nd_curl_matrix);

  const auto combined_h1 = fem::singular::BuildParallelEnrichedOperator(
      *standard_h1_matrix, parallel.h1_diffusion);
  const auto combined_nd_mass = fem::singular::BuildParallelEnrichedOperator(
      *standard_nd_mass_matrix, parallel.nd_mass);
  const auto combined_nd_curl = fem::singular::BuildParallelEnrichedOperator(
      *standard_nd_curl_matrix, parallel.nd_curl_curl);
  const auto combined_gradient = fem::singular::BuildParallelEnrichedGradient(
      *standard_gradient_matrix, *enrichment_gradient);
  REQUIRE(combined_h1);
  REQUIRE(combined_nd_mass);
  REQUIRE(combined_nd_curl);
  REQUIRE(combined_gradient);

  const Eigen::MatrixXd h1 = ToEigen(GatherParallelMatrix(*combined_h1));
  const Eigen::MatrixXd mass = ToEigen(GatherParallelMatrix(*combined_nd_mass));
  const Eigen::MatrixXd curl = ToEigen(GatherParallelMatrix(*combined_nd_curl));
  const Eigen::MatrixXd gradient = ToEigen(GatherParallelMatrix(*combined_gradient));
  const double scale = std::max(1.0, h1.cwiseAbs().maxCoeff());
  CHECK((h1 - gradient.transpose() * mass * gradient).cwiseAbs().maxCoeff() <
        2.0e-11 * scale);
  CHECK((curl * gradient).cwiseAbs().maxCoeff() <
        2.0e-11 * std::max(1.0, curl.cwiseAbs().maxCoeff()));
}

TEST_CASE("Parallel singular sparse assembly matches its serial true-DOF operator",
          "[singularelements][singularassembly][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Parallel sparse assembly is exercised by the [Parallel] test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  auto serial_mesh = SharedFaceTetrahedronMesh(false);
  std::vector<fem::singular::GlobalVertexId> serial_vertex_ids(serial_mesh.GetNV());
  for (int vertex = 0; vertex < serial_mesh.GetNV(); vertex++)
  {
    serial_vertex_ids[vertex] = vertex;
  }
  const auto serial_topology = SharedFaceDofTopology(serial_mesh, serial_vertex_ids);
  const std::array<int, 2> partition{0, 1};
  mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh, partition.data());
  const auto parallel_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
      serial_mesh, parallel_mesh, partition.data());
  const auto local_topology = SharedFaceDofTopology(parallel_mesh, parallel_vertex_ids);
  const auto numbering =
      fem::singular::BuildParallelDofNumbering(Mpi::World(), local_topology);

  mfem::H1_FECollection serial_h1_collection(2, 3);
  mfem::ND_FECollection serial_nd_collection(2, 3);
  mfem::FiniteElementSpace serial_h1_space(&serial_mesh, &serial_h1_collection);
  mfem::FiniteElementSpace serial_nd_space(&serial_mesh, &serial_nd_collection);
  mfem::H1_FECollection parallel_h1_collection(2, 3);
  mfem::ND_FECollection parallel_nd_collection(2, 3);
  mfem::ParFiniteElementSpace parallel_h1_space(&parallel_mesh, &parallel_h1_collection);
  mfem::ParFiniteElementSpace parallel_nd_space(&parallel_mesh, &parallel_nd_collection);

  const std::vector<fem::singular::IsotropicMaterialCoefficients> serial_materials{
      {2.3, 0.8}, {4.1, 1.7}};
  std::vector<fem::singular::IsotropicMaterialCoefficients> local_materials;
  for (int element = 0; element < parallel_mesh.GetNE(); element++)
  {
    bool upper_tetrahedron = false;
    for (int local = 0; local < 4; local++)
    {
      const int mesh_vertex = parallel_mesh.GetElement(element)->GetVertices()[local];
      upper_tetrahedron = upper_tetrahedron || parallel_vertex_ids[mesh_vertex] == 3;
    }
    local_materials.push_back(serial_materials[upper_tetrahedron ? 0 : 1]);
  }

  const fem::singular::AdaptiveAssemblyOptions options{8, 2.0e-6, 2.0e-6, 9};
  const auto serial = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      serial_topology, serial_h1_space, serial_nd_space, serial_materials, options);
  const auto serial_h1 = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
      serial_topology, serial_h1_space, serial_materials, options);
  const auto local = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      local_topology, parallel_h1_space, parallel_nd_space, local_materials, options);
  const auto local_h1 = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(
      local_topology, parallel_h1_space, local_materials, options);
  const auto parallel = fem::singular::AssembleParallelSparseEnrichmentMatrices(
      local, numbering, parallel_h1_space, parallel_nd_space);
  const auto parallel_h1 = fem::singular::AssembleParallelSparseH1EnrichmentMatrices(
      local_h1, numbering, parallel_h1_space);
  const auto enrichment_gradient =
      fem::singular::BuildParallelEnrichmentGradient(Mpi::World(), numbering);

  REQUIRE(numbering.h1.global_local_size == 2);
  REQUIRE(numbering.h1.global_size == 1);
  REQUIRE(numbering.nd.global_local_size == 4);
  REQUIRE(numbering.nd.global_size == 2);
  const auto check_enrichment_block =
      [](const fem::singular::LocalSparseOperatorBlocks &serial_block,
         const fem::singular::ParallelSparseOperatorBlocks &parallel_block)
  {
    REQUIRE(serial_block.enrichment_enrichment);
    REQUIRE(serial_block.enrichment_enrichment_estimated_absolute_error);
    REQUIRE(parallel_block.enrichment_enrichment);
    REQUIRE(parallel_block.enrichment_enrichment_estimated_absolute_error);
    REQUIRE(parallel_block.standard_enrichment);
    REQUIRE(parallel_block.enrichment_standard);
    CheckClose(GatherParallelMatrix(*parallel_block.enrichment_enrichment),
               DenseMatrix(*serial_block.enrichment_enrichment));
    CheckClose(GatherParallelMatrix(
                   *parallel_block.enrichment_enrichment_estimated_absolute_error),
               DenseMatrix(*serial_block.enrichment_enrichment_estimated_absolute_error));

    const auto forward = GatherParallelMatrix(*parallel_block.standard_enrichment);
    const auto reverse = GatherParallelMatrix(*parallel_block.enrichment_standard);
    REQUIRE(forward.Height() == reverse.Width());
    REQUIRE(forward.Width() == reverse.Height());
    for (int i = 0; i < forward.Height(); i++)
    {
      for (int j = 0; j < forward.Width(); j++)
      {
        CHECK(forward(i, j) == reverse(j, i));
      }
    }
  };
  check_enrichment_block(serial.h1_diffusion, parallel.h1_diffusion);
  check_enrichment_block(serial.nd_mass, parallel.nd_mass);
  check_enrichment_block(serial.nd_curl_curl, parallel.nd_curl_curl);
  check_enrichment_block(serial_h1.diffusion, parallel_h1);

  mfem::ParDiscreteLinearOperator standard_gradient(&parallel_h1_space, &parallel_nd_space);
  standard_gradient.AddDomainInterpolator(new mfem::GradientInterpolator);
  standard_gradient.Assemble();
  standard_gradient.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_gradient_matrix(
      standard_gradient.ParallelAssemble());
  REQUIRE(standard_gradient_matrix);

  mfem::Vector electric_coefficients(2);
  electric_coefficients[0] = serial_materials[0].electric;
  electric_coefficients[1] = serial_materials[1].electric;
  mfem::PWConstCoefficient electric_coefficient(electric_coefficients);
  mfem::ParBilinearForm standard_h1_form(&parallel_h1_space);
  standard_h1_form.AddDomainIntegrator(new mfem::DiffusionIntegrator(electric_coefficient));
  standard_h1_form.Assemble();
  standard_h1_form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_h1_matrix(
      standard_h1_form.ParallelAssemble());
  REQUIRE(standard_h1_matrix);

  mfem::ParBilinearForm standard_nd_mass_form(&parallel_nd_space);
  standard_nd_mass_form.AddDomainIntegrator(
      new mfem::VectorFEMassIntegrator(electric_coefficient));
  standard_nd_mass_form.Assemble();
  standard_nd_mass_form.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> standard_nd_mass_matrix(
      standard_nd_mass_form.ParallelAssemble());
  REQUIRE(standard_nd_mass_matrix);

  const auto combined_h1 = fem::singular::BuildParallelEnrichedOperator(
      *standard_h1_matrix, parallel.h1_diffusion);
  const auto combined_nd_mass = fem::singular::BuildParallelEnrichedOperator(
      *standard_nd_mass_matrix, parallel.nd_mass);
  const auto combined_gradient = fem::singular::BuildParallelEnrichedGradient(
      *standard_gradient_matrix, *enrichment_gradient);
  REQUIRE(combined_h1);
  REQUIRE(combined_nd_mass);
  REQUIRE(combined_gradient);
  CHECK(combined_h1->Height() ==
        parallel_h1_space.GetTrueVSize() + numbering.h1.owned_size);
  CHECK(combined_nd_mass->Height() ==
        parallel_nd_space.GetTrueVSize() + numbering.nd.owned_size);
  CHECK(combined_gradient->Height() == combined_nd_mass->Height());
  CHECK(combined_gradient->Width() == combined_h1->Width());

  const Eigen::MatrixXd gradient_enrichment =
      ToEigen(GatherParallelMatrix(*enrichment_gradient));
  const Eigen::MatrixXd gradient_standard =
      ToEigen(GatherParallelMatrix(*standard_gradient_matrix));
  REQUIRE(gradient_enrichment.rows() == 2);
  REQUIRE(gradient_enrichment.cols() == 1);
  CHECK(gradient_enrichment(0, 0) == 1.0);
  CHECK(gradient_enrichment(1, 0) == 0.0);

  const Eigen::MatrixXd h1_ee =
      ToEigen(GatherParallelMatrix(*parallel.h1_diffusion.enrichment_enrichment));
  const Eigen::MatrixXd nd_mass_ee =
      ToEigen(GatherParallelMatrix(*parallel.nd_mass.enrichment_enrichment));
  const Eigen::MatrixXd h1_se =
      ToEigen(GatherParallelMatrix(*parallel.h1_diffusion.standard_enrichment));
  const Eigen::MatrixXd nd_mass_se =
      ToEigen(GatherParallelMatrix(*parallel.nd_mass.standard_enrichment));
  const Eigen::MatrixXd nd_curl_ee =
      ToEigen(GatherParallelMatrix(*parallel.nd_curl_curl.enrichment_enrichment));
  const Eigen::MatrixXd nd_curl_se =
      ToEigen(GatherParallelMatrix(*parallel.nd_curl_curl.standard_enrichment));
  const auto check_eigen_close =
      [](const Eigen::MatrixXd &value, const Eigen::MatrixXd &reference)
  {
    REQUIRE(value.rows() == reference.rows());
    REQUIRE(value.cols() == reference.cols());
    for (Eigen::Index i = 0; i < value.rows(); i++)
    {
      for (Eigen::Index j = 0; j < value.cols(); j++)
      {
        CAPTURE(i, j);
        CheckClose(value(i, j), reference(i, j));
      }
    }
  };
  check_eigen_close(h1_ee,
                    gradient_enrichment.transpose() * nd_mass_ee * gradient_enrichment);
  check_eigen_close(h1_se,
                    gradient_standard.transpose() * nd_mass_se * gradient_enrichment);
  check_eigen_close(gradient_enrichment.transpose() * nd_curl_ee * gradient_enrichment,
                    Eigen::MatrixXd::Zero(h1_ee.rows(), h1_ee.cols()));
  check_eigen_close(nd_curl_se * gradient_enrichment,
                    Eigen::MatrixXd::Zero(nd_curl_se.rows(), gradient_enrichment.cols()));

  const Eigen::MatrixXd combined_h1_dense = ToEigen(GatherParallelMatrix(*combined_h1));
  const Eigen::MatrixXd combined_nd_mass_dense =
      ToEigen(GatherParallelMatrix(*combined_nd_mass));
  const Eigen::MatrixXd combined_gradient_dense =
      ToEigen(GatherParallelMatrix(*combined_gradient));
  check_eigen_close(combined_h1_dense, combined_gradient_dense.transpose() *
                                           combined_nd_mass_dense *
                                           combined_gradient_dense);

  mfem::Array<int> standard_essential, enrichment_essential;
  REQUIRE(parallel_h1_space.GetTrueVSize() > 0);
  standard_essential.Append(0);
  if (numbering.h1.owned_size > 0)
  {
    enrichment_essential.Append(0);
  }
  auto dirichlet = fem::singular::BuildParallelDirichletSystem(
      std::make_unique<mfem::HypreParMatrix>(*combined_h1),
      parallel_h1_space.GetTrueVSize(), standard_essential, enrichment_essential);
  REQUIRE(dirichlet.constrained);
  REQUIRE(dirichlet.eliminated);
  CHECK(dirichlet.essential_true_dofs.Size() ==
        standard_essential.Size() + enrichment_essential.Size());

  const auto constrained_blocks = fem::singular::BuildParallelConstrainedOperatorBlocks(
      *standard_h1_matrix, parallel.h1_diffusion, standard_essential, enrichment_essential);
  mfem::Array2D<const mfem::HypreParMatrix *> block_view(2, 2);
  block_view(0, 0) = constrained_blocks.standard_standard.get();
  block_view(0, 1) = constrained_blocks.standard_enrichment.get();
  block_view(1, 0) = constrained_blocks.enrichment_standard.get();
  block_view(1, 1) = constrained_blocks.enrichment_enrichment.get();
  std::unique_ptr<mfem::HypreParMatrix> blockwise_constrained(
      mfem::HypreParMatrixFromBlocks(block_view));
  REQUIRE(blockwise_constrained);
  CheckClose(GatherParallelMatrix(*blockwise_constrained),
             GatherParallelMatrix(*dirichlet.constrained));

  mfem::Array<int> no_enrichment_essential;
  auto correction_dirichlet = fem::singular::BuildParallelDirichletSystem(
      std::make_unique<mfem::HypreParMatrix>(*combined_h1),
      parallel_h1_space.GetTrueVSize(), standard_essential, no_enrichment_essential);
  const auto correction_blocks = fem::singular::BuildParallelConstrainedOperatorBlocks(
      *standard_h1_matrix, parallel.h1_diffusion, standard_essential,
      no_enrichment_essential);
  const auto coupled_patch = fem::singular::BuildParallelCoupledPatch(
      *correction_dirichlet.constrained, *correction_blocks.standard_enrichment,
      parallel_h1_space.GetTrueVSize());
  REQUIRE(coupled_patch.matrix);
  CHECK(coupled_patch.global_standard_dofs > 0);
  CHECK(coupled_patch.global_enrichment_dofs == numbering.h1.global_size);
  CHECK(coupled_patch.matrix->GetGlobalNumRows() ==
        coupled_patch.global_standard_dofs + coupled_patch.global_enrichment_dofs);

  const std::vector<std::vector<std::size_t>> one_feature_membership(
      local_topology.h1_dofs.size(), std::vector<std::size_t>{0});
  const auto feature_patches = fem::singular::BuildParallelFeaturePatches(
      *correction_dirichlet.constrained, *correction_blocks.standard_enrichment,
      parallel_h1_space.GetTrueVSize(), one_feature_membership, numbering.h1, 1);
  REQUIRE(feature_patches.patches.size() == 1);
  const auto &feature_patch = feature_patches.patches.front();
  CHECK(feature_patch.feature == 0);
  CHECK(feature_patch.true_dofs == coupled_patch.true_dofs);
  CHECK(feature_patch.global_standard_dofs == coupled_patch.global_standard_dofs);
  CHECK(feature_patch.global_enrichment_dofs == coupled_patch.global_enrichment_dofs);
  CHECK(feature_patches.sum_global_standard_dofs == coupled_patch.global_standard_dofs);
  CHECK(feature_patches.sum_global_enrichment_dofs == coupled_patch.global_enrichment_dofs);
  CHECK(feature_patches.maximum_global_standard_dofs == coupled_patch.global_standard_dofs);
  CHECK(feature_patches.maximum_global_enrichment_dofs ==
        coupled_patch.global_enrichment_dofs);
  CHECK(feature_patches.minimum_enrichment_multiplicity == 1);
  CHECK(feature_patches.maximum_enrichment_multiplicity == 1);
  CheckClose(GatherParallelMatrix(*feature_patch.matrix),
             GatherParallelMatrix(*coupled_patch.matrix));

  auto uncovered_membership = one_feature_membership;
  REQUIRE_FALSE(uncovered_membership.empty());
  uncovered_membership.front().clear();
  CHECK_THROWS_AS(
      fem::singular::BuildParallelFeaturePatches(
          *correction_dirichlet.constrained, *correction_blocks.standard_enrichment,
          parallel_h1_space.GetTrueVSize(), uncovered_membership, numbering.h1, 1),
      std::invalid_argument);

  auto overlapping_membership = one_feature_membership;
  for (auto &membership : overlapping_membership)
  {
    membership = {0, 1};
  }
  const auto overlapping_feature_patches = fem::singular::BuildParallelFeaturePatches(
      *correction_dirichlet.constrained, *correction_blocks.standard_enrichment,
      parallel_h1_space.GetTrueVSize(), overlapping_membership, numbering.h1, 2);
  REQUIRE(overlapping_feature_patches.patches.size() == 2);
  CHECK(overlapping_feature_patches.minimum_enrichment_multiplicity == 2);
  CHECK(overlapping_feature_patches.maximum_enrichment_multiplicity == 2);
  CHECK(overlapping_feature_patches.sum_global_enrichment_dofs ==
        2 * numbering.h1.global_size);
  CHECK(overlapping_feature_patches.patches[0].true_dofs == coupled_patch.true_dofs);
  CheckClose(GatherParallelMatrix(*overlapping_feature_patches.patches[0].matrix),
             GatherParallelMatrix(*coupled_patch.matrix));
  CHECK(overlapping_feature_patches.patches[1].global_enrichment_dofs > 0);

  auto inconsistent_membership = one_feature_membership;
  bool changed_shared_membership = false;
  for (std::size_t i = 0; i < inconsistent_membership.size(); i++)
  {
    if (numbering.h1.owner[i] != Mpi::Rank(Mpi::World()))
    {
      inconsistent_membership[i] = {1};
      changed_shared_membership = true;
      break;
    }
  }
  Mpi::GlobalOr(1, &changed_shared_membership, Mpi::World());
  REQUIRE(changed_shared_membership);
  CHECK_THROWS_AS(
      fem::singular::BuildParallelFeaturePatches(
          *correction_dirichlet.constrained, *correction_blocks.standard_enrichment,
          parallel_h1_space.GetTrueVSize(), inconsistent_membership, numbering.h1, 2),
      std::invalid_argument);

  mfem::Vector patch_vector(coupled_patch.true_dofs.Size());
  mfem::Vector restricted_full(correction_dirichlet.constrained->Width());
  restricted_full = 0.0;
  for (int i = 0; i < patch_vector.Size(); i++)
  {
    patch_vector[i] = 0.17 * (1 + Mpi::Rank(Mpi::World())) - 0.09 * i;
    restricted_full[coupled_patch.true_dofs[i]] = patch_vector[i];
  }
  CheckClose(ParallelBilinear(patch_vector, *coupled_patch.matrix, patch_vector),
             ParallelBilinear(restricted_full, *correction_dirichlet.constrained,
                              restricted_full));

  const auto make_exact_cg = []()
  {
    auto solver = std::make_unique<mfem::CGSolver>(Mpi::World());
    solver->SetRelTol(1.0e-14);
    solver->SetAbsTol(1.0e-15);
    solver->SetMaxIter(100);
    solver->SetPrintLevel(0);
    return solver;
  };
  std::vector<mfem::Array<int>> additive_patch_dofs;
  std::vector<std::unique_ptr<mfem::Solver>> additive_patch_solvers;
  std::vector<const Operator *> additive_patch_operators;
  for (const auto &patch : overlapping_feature_patches.patches)
  {
    additive_patch_dofs.push_back(patch.true_dofs);
    additive_patch_solvers.push_back(make_exact_cg());
    additive_patch_operators.push_back(patch.matrix.get());
  }
  auto additive_patch = std::make_unique<AdditivePatchSolver>(
      correction_dirichlet.constrained->Height(), std::move(additive_patch_dofs),
      std::move(additive_patch_solvers));
  additive_patch->SetPatchOperators(additive_patch_operators);
  SymmetricPatchSubspacePreconditioner patch_preconditioner(
      parallel_h1_space.GetTrueVSize(), make_exact_cg(), std::move(additive_patch));
  patch_preconditioner.SetSubspaceOperators(*correction_dirichlet.constrained,
                                            *correction_blocks.standard_standard);

  mfem::Vector test_left(correction_dirichlet.constrained->Width());
  mfem::Vector test_right(correction_dirichlet.constrained->Width());
  for (int i = 0; i < test_left.Size(); i++)
  {
    test_left[i] = 0.11 * (1 + Mpi::Rank(Mpi::World())) + 0.03 * i;
    test_right[i] = -0.07 * (1 + Mpi::Rank(Mpi::World())) + 0.05 * i;
  }
  mfem::Vector corrected_left, corrected_right;
  patch_preconditioner.Mult(test_left, corrected_left);
  patch_preconditioner.Mult(test_right, corrected_right);
  double left_right = test_left * corrected_right;
  double right_left = test_right * corrected_left;
  double left_positive = test_left * corrected_left;
  Mpi::GlobalSum(1, &left_right, Mpi::World());
  Mpi::GlobalSum(1, &right_left, Mpi::World());
  Mpi::GlobalSum(1, &left_positive, Mpi::World());
  CheckClose(left_right, right_left);
  CHECK(left_positive > 0.0);

  mfem::Vector patch_exact(correction_dirichlet.constrained->Width());
  for (int i = 0; i < patch_exact.Size(); i++)
  {
    patch_exact[i] = 0.19 * (1 + Mpi::Rank(Mpi::World())) - 0.04 * i;
  }
  mfem::Vector patch_rhs(correction_dirichlet.constrained->Height());
  correction_dirichlet.constrained->Mult(patch_exact, patch_rhs);
  CgSolver<Operator> patch_solver(Mpi::World(), 0);
  patch_solver.SetOperator(*correction_dirichlet.constrained);
  patch_solver.SetPreconditioner(patch_preconditioner);
  patch_solver.SetRelTol(1.0e-12);
  patch_solver.SetAbsTol(1.0e-14);
  patch_solver.SetMaxIter(100);
  mfem::Vector patch_solution(patch_exact.Size());
  patch_solution = 0.0;
  patch_solver.Mult(patch_rhs, patch_solution);
  CHECK(patch_solver.GetConverged());
  for (int i = 0; i < patch_solution.Size(); i++)
  {
    CAPTURE(i, patch_solver.GetNumIterations(), patch_solver.GetFinalRes());
    CHECK(std::abs(patch_solution[i] - patch_exact[i]) <=
          2.0e-10 * std::max(1.0, std::abs(patch_exact[i])));
  }

  mfem::Vector exact(combined_h1->Width()), prescribed(combined_h1->Width());
  for (int i = 0; i < exact.Size(); i++)
  {
    exact[i] = 0.13 * (1 + Mpi::Rank(Mpi::World())) - 0.07 * i;
  }
  prescribed = 0.0;
  for (const int dof : dirichlet.essential_true_dofs)
  {
    prescribed[dof] = exact[dof];
  }
  mfem::Vector rhs(combined_h1->Height()), constrained_action(combined_h1->Height());
  combined_h1->Mult(exact, rhs);
  dirichlet.EliminateRHS(prescribed, rhs);
  dirichlet.constrained->Mult(exact, constrained_action);
  for (int i = 0; i < rhs.Size(); i++)
  {
    CAPTURE(i);
    CheckClose(rhs[i], constrained_action[i]);
  }

  fem::singular::SymmetricDiagonalScaling diagonal_scaling(*dirichlet.constrained);
  CHECK(diagonal_scaling.GetScaledDiagonalMinimum() ==
        Catch::Approx(1.0).epsilon(64.0 * std::numeric_limits<double>::epsilon()));
  CHECK(diagonal_scaling.GetScaledDiagonalMaximum() ==
        Catch::Approx(1.0).epsilon(64.0 * std::numeric_limits<double>::epsilon()));
  mfem::Vector scaled_exact, scaled_rhs, recovered_exact, scaled_action;
  diagonal_scaling.ScaleInitialGuess(exact, scaled_exact);
  diagonal_scaling.ScaleRHS(rhs, scaled_rhs);
  diagonal_scaling.RecoverSolution(scaled_exact, recovered_exact);
  scaled_action.SetSize(diagonal_scaling.GetMatrix().Height());
  diagonal_scaling.GetMatrix().Mult(scaled_exact, scaled_action);
  for (int i = 0; i < exact.Size(); i++)
  {
    CAPTURE(i);
    CheckClose(recovered_exact[i], exact[i]);
    CheckClose(scaled_action[i], scaled_rhs[i]);
  }
  CheckClose(ParallelBilinear(scaled_exact, diagonal_scaling.GetMatrix(), scaled_exact),
             ParallelBilinear(exact, *dirichlet.constrained, exact));

  mfem::CGSolver solver(Mpi::World());
  solver.SetOperator(*dirichlet.constrained);
  solver.SetRelTol(1.0e-13);
  solver.SetAbsTol(1.0e-15);
  solver.SetMaxIter(200);
  solver.SetPrintLevel(0);
  mfem::Vector solution(exact.Size());
  solution = 0.0;
  solver.Mult(rhs, solution);
  CHECK(solver.GetConverged());
  for (int i = 0; i < solution.Size(); i++)
  {
    CAPTURE(i, solver.GetNumIterations(), solver.GetFinalNorm());
    CHECK(std::abs(solution[i] - exact[i]) <= 2.0e-10 * std::max(1.0, std::abs(exact[i])));
  }

  mfem::Vector prescribed_action(combined_h1->Height());
  dirichlet.constrained->Mult(prescribed, prescribed_action);
  for (int i = 0; i < prescribed_action.Size(); i++)
  {
    const bool essential = std::find(dirichlet.essential_true_dofs.begin(),
                                     dirichlet.essential_true_dofs.end(),
                                     i) != dirichlet.essential_true_dofs.end();
    CAPTURE(i, essential);
    CheckClose(prescribed_action[i], essential ? prescribed[i] : 0.0);
  }

  mfem::Array<int> duplicate_standard_essential;
  duplicate_standard_essential.Append(0);
  duplicate_standard_essential.Append(0);
  CHECK_THROWS_AS(fem::singular::BuildParallelDirichletSystem(
                      std::make_unique<mfem::HypreParMatrix>(*combined_h1),
                      parallel_h1_space.GetTrueVSize(), duplicate_standard_essential,
                      enrichment_essential),
                  std::invalid_argument);

  mfem::FunctionCoefficient h1_coefficient(
      [](const mfem::Vector &x) { return x[0] + 2.0 * x[1] - 0.5 * x[2] + x[0] * x[1]; });
  mfem::VectorFunctionCoefficient nd_coefficient(
      3,
      [](const mfem::Vector &x, mfem::Vector &value)
      {
        value[0] = x[1] + 2.0 * x[2];
        value[1] = -x[0] + 0.5 * x[2];
        value[2] = 3.0 * x[0] - x[1];
      });
  mfem::GridFunction serial_h1_field(&serial_h1_space);
  mfem::GridFunction serial_nd_field(&serial_nd_space);
  mfem::ParGridFunction parallel_h1_field(&parallel_h1_space);
  mfem::ParGridFunction parallel_nd_field(&parallel_nd_space);
  serial_h1_field.ProjectCoefficient(h1_coefficient);
  serial_nd_field.ProjectCoefficient(nd_coefficient);
  parallel_h1_field.ProjectCoefficient(h1_coefficient);
  parallel_nd_field.ProjectCoefficient(nd_coefficient);
  mfem::Vector parallel_h1_true(parallel_h1_space.GetTrueVSize());
  mfem::Vector parallel_nd_true(parallel_nd_space.GetTrueVSize());
  parallel_h1_field.ParallelProject(parallel_h1_true);
  parallel_nd_field.ParallelProject(parallel_nd_true);

  mfem::Vector h1_enrichment(1);
  h1_enrichment[0] = 0.37;
  mfem::Vector nd_enrichment(2);
  nd_enrichment[0] = 0.37;
  nd_enrichment[1] = -0.21;
  mfem::Vector parallel_h1_enrichment(numbering.h1.owned_size);
  mfem::Vector parallel_nd_enrichment(numbering.nd.owned_size);
  parallel_h1_enrichment = 0.0;
  parallel_nd_enrichment = 0.0;
  for (std::size_t local_dof = 0; local_dof < local_topology.h1_dofs.size(); local_dof++)
  {
    if (numbering.h1.owner[local_dof] == Mpi::Rank(Mpi::World()))
    {
      parallel_h1_enrichment[numbering.h1.local_to_true[local_dof] -
                             numbering.h1.owned_offset] = h1_enrichment[0];
    }
  }
  for (std::size_t local_dof = 0; local_dof < local_topology.nd_dofs.size(); local_dof++)
  {
    if (numbering.nd.owner[local_dof] == Mpi::Rank(Mpi::World()))
    {
      const auto family = local_topology.nd_dofs[local_dof].family;
      const double value = (family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
                               ? nd_enrichment[0]
                               : nd_enrichment[1];
      parallel_nd_enrichment[numbering.nd.local_to_true[local_dof] -
                             numbering.nd.owned_offset] = value;
    }
  }

  CheckClose(ParallelBilinear(parallel_h1_true, *parallel.h1_diffusion.standard_enrichment,
                              parallel_h1_enrichment),
             SparseBilinear(serial_h1_field, *serial.h1_diffusion.standard_enrichment,
                            h1_enrichment));
  CheckClose(
      ParallelBilinear(parallel_nd_true, *parallel.nd_mass.standard_enrichment,
                       parallel_nd_enrichment),
      SparseBilinear(serial_nd_field, *serial.nd_mass.standard_enrichment, nd_enrichment));
  CheckClose(ParallelBilinear(parallel_nd_true, *parallel.nd_curl_curl.standard_enrichment,
                              parallel_nd_enrichment),
             SparseBilinear(serial_nd_field, *serial.nd_curl_curl.standard_enrichment,
                            nd_enrichment));

  const auto serial_h1_absolute = AbsoluteValue(serial_h1_field);
  const auto serial_nd_absolute = AbsoluteValue(serial_nd_field);
  const auto h1_enrichment_absolute = AbsoluteValue(h1_enrichment);
  const auto nd_enrichment_absolute = AbsoluteValue(nd_enrichment);
  const auto parallel_h1_absolute = AbsoluteValue(parallel_h1_true);
  const auto parallel_nd_absolute = AbsoluteValue(parallel_nd_true);
  const auto parallel_h1_enrichment_absolute = AbsoluteValue(parallel_h1_enrichment);
  const auto parallel_nd_enrichment_absolute = AbsoluteValue(parallel_nd_enrichment);
  const auto check_conservative_bound = [](double parallel_bound, double serial_bound)
  {
    CHECK(parallel_bound + 256.0 * std::numeric_limits<double>::epsilon() *
                               std::max({1.0, parallel_bound, serial_bound}) >=
          serial_bound);
  };
  check_conservative_bound(
      ParallelBilinear(parallel_h1_absolute,
                       *parallel.h1_diffusion.standard_enrichment_estimated_absolute_error,
                       parallel_h1_enrichment_absolute),
      SparseBilinear(serial_h1_absolute,
                     *serial.h1_diffusion.standard_enrichment_estimated_absolute_error,
                     h1_enrichment_absolute));
  check_conservative_bound(
      ParallelBilinear(parallel_nd_absolute,
                       *parallel.nd_mass.standard_enrichment_estimated_absolute_error,
                       parallel_nd_enrichment_absolute),
      SparseBilinear(serial_nd_absolute,
                     *serial.nd_mass.standard_enrichment_estimated_absolute_error,
                     nd_enrichment_absolute));
  check_conservative_bound(
      ParallelBilinear(parallel_nd_absolute,
                       *parallel.nd_curl_curl.standard_enrichment_estimated_absolute_error,
                       parallel_nd_enrichment_absolute),
      SparseBilinear(serial_nd_absolute,
                     *serial.nd_curl_curl.standard_enrichment_estimated_absolute_error,
                     nd_enrichment_absolute));

  auto inconsistent_numbering = numbering;
  if (Mpi::Rank(Mpi::World()) == 0)
  {
    inconsistent_numbering.h1.local_offset++;
  }
  CHECK_THROWS_AS(fem::singular::AssembleParallelSparseEnrichmentMatrices(
                      local, inconsistent_numbering, parallel_h1_space, parallel_nd_space),
                  std::invalid_argument);
}

}  // namespace palace
