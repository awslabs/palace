// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/hierarchicalerrorestimator.hpp"
#include "fem/singularassembly.hpp"
#include "fem/singulardofs.hpp"
#include "fem/singularelements.hpp"
#include "fem/singularfeatures.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch::Matchers;
namespace
{

constexpr double MaxwellNu = 2.0 / 3.0;
constexpr double MaxwellBeta = 0.35;
constexpr double MaxwellGamma = 1.0;
constexpr int MaxwellReentrantAttribute = 1;
constexpr int MaxwellOuterAttribute = 2;

mfem::Mesh MaxwellLShapeMesh(int n)
{
  struct Key
  {
    int i, j;
    bool operator<(const Key &other) const
    {
      return i != other.i ? i < other.i : j < other.j;
    }
  };
  std::map<Key, int> ids;
  std::vector<std::array<double, 2>> points;
  const auto vertex = [&](int i, int j)
  {
    const Key key{i, j};
    if (const auto found = ids.find(key); found != ids.end())
    {
      return found->second;
    }
    const int id = static_cast<int>(points.size());
    ids.emplace(key, id);
    points.push_back({static_cast<double>(i) / n, static_cast<double>(j) / n});
    return id;
  };
  std::vector<std::array<int, 3>> triangles;
  for (const auto &origin : std::array<std::array<int, 2>, 3>{{{-n, 0}, {0, 0}, {-n, -n}}})
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        const int x = origin[0] + i, y = origin[1] + j;
        const int p00 = vertex(x, y), p10 = vertex(x + 1, y);
        const int p01 = vertex(x, y + 1), p11 = vertex(x + 1, y + 1);
        triangles.push_back({p00, p10, p11});
        triangles.push_back({p00, p11, p01});
      }
    }
  }
  std::vector<std::array<int, 3>> segments;
  for (int i = 0; i < n; i++)
  {
    segments.push_back({vertex(i, 0), vertex(i + 1, 0), MaxwellReentrantAttribute});
    segments.push_back({vertex(0, -i - 1), vertex(0, -i), MaxwellReentrantAttribute});
    segments.push_back({vertex(-n, i), vertex(-n, i + 1), MaxwellOuterAttribute});
    segments.push_back({vertex(-n, -i - 1), vertex(-n, -i), MaxwellOuterAttribute});
    segments.push_back({vertex(n, i), vertex(n, i + 1), MaxwellOuterAttribute});
    segments.push_back({vertex(-n + i, -n), vertex(-n + i + 1, -n), MaxwellOuterAttribute});
    segments.push_back({vertex(-n + i, n), vertex(-n + i + 1, n), MaxwellOuterAttribute});
    segments.push_back({vertex(i, n), vertex(i + 1, n), MaxwellOuterAttribute});
  }
  mfem::Mesh mesh(2, static_cast<int>(points.size()), static_cast<int>(triangles.size()),
                  static_cast<int>(segments.size()), 2);
  for (const auto &point : points)
  {
    mesh.AddVertex(point[0], point[1]);
  }
  for (const auto &triangle : triangles)
  {
    mesh.AddTriangle(triangle.data(), 1);
  }
  for (const auto &segment : segments)
  {
    mesh.AddBdrSegment(segment[0], segment[1], segment[2]);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

double MaxwellS(double x, double y)
{
  const double r = std::hypot(x, y);
  if (r == 0.0)
  {
    return 0.0;
  }
  double theta = std::atan2(y, x);
  if (theta < 0.0)
  {
    theta += 2.0 * M_PI;
  }
  return std::pow(r, MaxwellNu) * std::sin(MaxwellNu * theta);
}

std::array<double, 2> MaxwellGradS(double x, double y)
{
  const double r = std::hypot(x, y);
  double theta = std::atan2(y, x);
  if (theta < 0.0)
  {
    theta += 2.0 * M_PI;
  }
  const double scale = MaxwellNu * std::pow(r, MaxwellNu - 1.0);
  return {scale * std::sin((MaxwellNu - 1.0) * theta),
          scale * std::cos((MaxwellNu - 1.0) * theta)};
}

double MaxwellU(double x, double y)
{
  return (1.0 - x * x) * (1.0 - y * y) * MaxwellS(x, y);
}

std::array<double, 2> MaxwellGradU(double x, double y)
{
  const double B = (1.0 - x * x) * (1.0 - y * y);
  const auto gradient = MaxwellGradS(x, y);
  const double singular = MaxwellS(x, y);
  return {B * gradient[0] - 2.0 * x * (1.0 - y * y) * singular,
          B * gradient[1] - 2.0 * y * (1.0 - x * x) * singular};
}

double MaxwellC(double x, double y)
{
  return x * y * (1.0 - x * x) * (1.0 - y * y);
}

std::array<double, 2> MaxwellField(double x, double y)
{
  const auto gradient = MaxwellGradU(x, y);
  return {gradient[0] + MaxwellBeta * MaxwellC(x, y), gradient[1]};
}

double MaxwellCurl(double x, double y)
{
  return -MaxwellBeta * x * (1.0 - x * x) * (1.0 - 3.0 * y * y);
}

std::array<double, 2> MaxwellSource(double x, double y)
{
  const auto field = MaxwellField(x, y);
  const double curl_curl_x = 6.0 * MaxwellBeta * x * y * (1.0 - x * x);
  const double curl_curl_y = MaxwellBeta * (1.0 - 3.0 * x * x) * (1.0 - 3.0 * y * y);
  return {curl_curl_x + MaxwellGamma * field[0], curl_curl_y + MaxwellGamma * field[1]};
}

int MaxwellCornerNode(const mfem::Mesh &mesh, int element)
{
  mfem::Array<int> vertices;
  mesh.GetElementVertices(element, vertices);
  for (int local = 0; local < vertices.Size(); local++)
  {
    const double *point = mesh.GetVertex(vertices[local]);
    if (std::hypot(point[0], point[1]) < 1.0e-13)
    {
      return local;
    }
  }
  return -1;
}

template <typename Visitor>
void MaxwellQuadrature(mfem::Mesh &mesh, int element, const Visitor &visitor)
{
  if (const int node = MaxwellCornerNode(mesh, element); node >= 0)
  {
    fem::singular::ForEachReferenceTriangleNodeDuffyQuadraturePoint(18, node, 6.0, visitor);
  }
  else
  {
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, 18);
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &ip = rule.IntPoint(q);
      visitor(fem::singular::TriangleBarycentricPoint{1.0 - ip.x - ip.y, ip.x, ip.y},
              ip.weight);
    }
  }
}

fem::singular::TriangleVectorBasisValue
MaxwellSingularBasis(const fem::singular::TriangleBarycentricPoint &lambda,
                     const fem::singular::TriangleBarycentricGradients &grad_lambda,
                     const fem::singular::TriangleBasis &basis)
{
  if (basis.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
  {
    return fem::singular::EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                                       basis.nodes[1], basis.nu);
  }
  return fem::singular::EvaluateTriangleNodeRotational(
      lambda, grad_lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nu);
}

struct MaxwellSolveReport
{
  long long dofs = 0;
  long long enrichment_dofs = 0;
  long long free_enrichment_dofs = 0;
  long long free_gradient_dofs = 0;
  long long free_rotational_dofs = 0;
  double residual = 0.0;
  double graph_error = 0.0;
  double exact_mass_energy = 0.0;
  double exact_curl_energy = 0.0;
  double gradient_roundtrip_error = 0.0;
  double rotational_roundtrip_error = 0.0;
  mfem::Vector standard_coefficients;
  mfem::Vector enrichment_coefficients;
  std::vector<bool> essential_mask;
  std::vector<double> element_error;
  std::shared_ptr<mfem::SparseMatrix> combined_matrix;
  mfem::Vector combined_rhs;
};

MaxwellSolveReport MaxwellSolve(mfem::Mesh &mesh, int order, bool enriched,
                                bool include_rotational = true,
                                bool synthetic_roundtrip = false)
{
  mfem::ND_FECollection nd_collection(order, 2);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  mfem::H1_FECollection h1_collection(order, 2);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  const int standard_size = nd_space.GetVSize();

  fem::singular::TriangleFeatureTopology features;
  fem::singular::TriangleDofTopology topology;
  fem::singular::LocalSparseEnrichmentMatrices singular;
  int enrichment_size = 0;
  if (enriched)
  {
    features = fem::singular::ExtractSerialLineFeatures(mesh, {MaxwellReentrantAttribute},
                                                        {{1, 1.0}});
    REQUIRE_THAT(features.vertices.front().nu, WithinRel(MaxwellNu, 1.0e-12));
    topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
    enrichment_size = static_cast<int>(topology.nd_dofs.size());
    const std::vector<fem::singular::IsotropicMaterialCoefficients> materials(mesh.GetNE(),
                                                                              {1.0, 1.0});
    singular = fem::singular::AssembleLocalSparseEnrichmentMatrices(
        topology, h1_space, nd_space, materials, {8, 1.0e-9, 1.0e-9, 8});
  }
  const int combined_size = standard_size + enrichment_size;
  mfem::SparseMatrix matrix(combined_size, combined_size);
  {
    mfem::BilinearForm graph(&nd_space);
    graph.AddDomainIntegrator(new mfem::CurlCurlIntegrator());
    graph.AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(new mfem::ConstantCoefficient(MaxwellGamma)));
    graph.Assemble();
    graph.Finalize();
    const auto &standard = graph.SpMat();
    for (int row = 0; row < standard.Height(); row++)
    {
      for (int entry = standard.GetI()[row]; entry < standard.GetI()[row + 1]; entry++)
      {
        matrix.Add(row, standard.GetJ()[entry], standard.GetData()[entry]);
      }
    }
  }
  if (enriched)
  {
    const auto add = [&](const mfem::SparseMatrix &block, int row_offset, int column_offset,
                         double scale)
    {
      for (int row = 0; row < block.Height(); row++)
      {
        for (int entry = block.GetI()[row]; entry < block.GetI()[row + 1]; entry++)
        {
          matrix.Add(row_offset + row, column_offset + block.GetJ()[entry],
                     scale * block.GetData()[entry]);
        }
      }
    };
    for (const auto &[blocks, scale] :
         std::array<std::pair<const fem::singular::LocalSparseOperatorBlocks *, double>, 2>{
             std::pair{&singular.nd_curl_curl, 1.0},
             std::pair{&singular.nd_mass, MaxwellGamma}})
    {
      add(*blocks->standard_enrichment, 0, standard_size, scale);
      add(*blocks->enrichment_standard, standard_size, 0, scale);
      add(*blocks->enrichment_enrichment, standard_size, standard_size, scale);
    }
  }
  matrix.Finalize();

  mfem::Vector rhs(combined_size);
  rhs = 0.0;
  mfem::Array<int> standard_dofs;
  mfem::DenseMatrix shape;
  mfem::Vector point(3);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    const auto &fe = *nd_space.GetFE(element);
    mfem::DofTransformation transformation;
    nd_space.GetElementVDofs(element, standard_dofs, transformation);
    mfem::Vector local(fe.GetDof());
    local = 0.0;
    shape.SetSize(fe.GetDof(), 2);
    double determinant;
    const auto grad_lambda =
        fem::singular::GetAffineTriangleBarycentricGradients(T, determinant);
    MaxwellQuadrature(
        mesh, element,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          T.Transform(ip, point);
          fe.CalcPhysVShape(T, shape);
          const auto source = MaxwellSource(point(0), point(1));
          const double weight = weight_ref * T.Weight();
          for (int i = 0; i < fe.GetDof(); i++)
          {
            local(i) += weight * (shape(i, 0) * source[0] + shape(i, 1) * source[1]);
          }
          if (enriched)
          {
            for (const auto &dof : topology.elements[element].nd)
            {
              const auto basis = MaxwellSingularBasis(lambda, grad_lambda, dof.basis);
              rhs(standard_size + static_cast<int>(dof.dof)) +=
                  weight * (basis.value[0] * source[0] + basis.value[1] * source[1]);
            }
          }
        });
    transformation.TransformDual(local);
    rhs.AddElementVector(standard_dofs, local);
  }

  std::vector<bool> essential(combined_size, false);
  mfem::Array<int> standard_essential;
  nd_space.GetBoundaryTrueDofs(standard_essential);
  for (int dof : standard_essential)
  {
    essential[dof] = true;
  }
  if (enriched)
  {
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), topology);
    mfem::Array<int> attributes(2);
    attributes[0] = MaxwellReentrantAttribute;
    attributes[1] = MaxwellOuterAttribute;
    const auto singular_essential = fem::singular::GetEssentialTriangleNDTrueDofs(
        Mpi::World(), features, topology, numbering, attributes);
    for (int dof : singular_essential)
    {
      essential[standard_size + dof] = true;
    }
    if (!include_rotational)
    {
      for (int dof = 0; dof < enrichment_size; dof++)
      {
        if (topology.nd_dofs[dof].family ==
            fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL)
        {
          essential[standard_size + dof] = true;
        }
      }
    }
  }

  mfem::Vector expected(combined_size);
  expected = 0.0;
  int expected_gradient = -1, expected_rotational = -1;
  if (synthetic_roundtrip)
  {
    REQUIRE(enriched);
    REQUIRE(include_rotational);
    for (int dof = 0; dof < enrichment_size; dof++)
    {
      if (essential[standard_size + dof])
      {
        continue;
      }
      if (expected_gradient < 0 && topology.nd_dofs[dof].family ==
                                       fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
      {
        expected_gradient = standard_size + dof;
        expected(expected_gradient) = 0.37;
      }
      if (expected_rotational < 0 &&
          topology.nd_dofs[dof].family ==
              fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL)
      {
        expected_rotational = standard_size + dof;
        expected(expected_rotational) = -0.29;
      }
    }
    REQUIRE(expected_gradient >= 0);
    REQUIRE(expected_rotational >= 0);
    matrix.Mult(expected, rhs);
  }

  std::vector<int> free, combined_to_free(combined_size, -1);
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (!essential[dof])
    {
      combined_to_free[dof] = static_cast<int>(free.size());
      free.push_back(dof);
    }
  }
  mfem::SparseMatrix reduced(static_cast<int>(free.size()), static_cast<int>(free.size()));
  mfem::Vector reduced_rhs(static_cast<int>(free.size()));
  for (int row = 0; row < static_cast<int>(free.size()); row++)
  {
    reduced_rhs(row) = rhs(free[row]);
    for (int entry = matrix.GetI()[free[row]]; entry < matrix.GetI()[free[row] + 1];
         entry++)
    {
      const int column = combined_to_free[matrix.GetJ()[entry]];
      if (column >= 0)
      {
        reduced.Add(row, column, matrix.GetData()[entry]);
      }
    }
  }
  reduced.Finalize();
  mfem::Vector reduced_solution(static_cast<int>(free.size()));
  reduced_solution = 0.0;
  mfem::GSSmoother preconditioner(reduced);
  mfem::CGSolver cg;
  cg.SetOperator(reduced);
  cg.SetPreconditioner(preconditioner);
  cg.SetRelTol(1.0e-13);
  cg.SetAbsTol(1.0e-30);
  cg.SetMaxIter(50000);
  cg.SetPrintLevel(-1);
  cg.Mult(reduced_rhs, reduced_solution);
  REQUIRE(cg.GetConverged());
  mfem::Vector solution(combined_size);
  solution = 0.0;
  for (int row = 0; row < static_cast<int>(free.size()); row++)
  {
    solution(free[row]) = reduced_solution(row);
  }

  MaxwellSolveReport report;
  report.dofs = combined_size;
  report.enrichment_dofs = enrichment_size;
  for (int dof = 0; dof < enrichment_size; dof++)
  {
    const bool is_free = !essential[standard_size + dof];
    report.free_enrichment_dofs += is_free;
    if (is_free && topology.nd_dofs[dof].family ==
                       fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
    {
      report.free_gradient_dofs++;
    }
    if (is_free && topology.nd_dofs[dof].family ==
                       fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL)
    {
      report.free_rotational_dofs++;
    }
  }
  mfem::Vector algebraic(static_cast<int>(free.size()));
  reduced.Mult(reduced_solution, algebraic);
  algebraic -= reduced_rhs;
  report.residual = algebraic.Norml2() / std::max(reduced_rhs.Norml2(), 1.0e-30);
  report.standard_coefficients.SetSize(standard_size);
  for (int i = 0; i < standard_size; i++)
  {
    report.standard_coefficients(i) = solution(i);
  }
  report.enrichment_coefficients.SetSize(enrichment_size);
  for (int i = 0; i < enrichment_size; i++)
  {
    report.enrichment_coefficients(i) = solution(standard_size + i);
  }
  report.essential_mask = essential;
  report.combined_matrix = std::make_shared<mfem::SparseMatrix>(matrix);
  report.combined_rhs = rhs;
  if (synthetic_roundtrip)
  {
    report.gradient_roundtrip_error =
        std::abs(solution(expected_gradient) - expected(expected_gradient));
    report.rotational_roundtrip_error =
        std::abs(solution(expected_rotational) - expected(expected_rotational));
    return report;
  }

  mfem::GridFunction standard_solution(&nd_space);
  for (int i = 0; i < standard_size; i++)
  {
    standard_solution(i) = solution(i);
  }
  mfem::Vector discrete(2), discrete_curl;
  report.element_error.assign(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double determinant;
    const auto grad_lambda =
        fem::singular::GetAffineTriangleBarycentricGradients(T, determinant);
    MaxwellQuadrature(
        mesh, element,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          T.Transform(ip, point);
          standard_solution.GetVectorValue(T, ip, discrete);
          standard_solution.GetCurl(T, discrete_curl);
          double curl = discrete_curl(0);
          if (enriched)
          {
            for (const auto &dof : topology.elements[element].nd)
            {
              const double coefficient =
                  solution(standard_size + static_cast<int>(dof.dof));
              const auto basis = MaxwellSingularBasis(lambda, grad_lambda, dof.basis);
              discrete(0) += coefficient * basis.value[0];
              discrete(1) += coefficient * basis.value[1];
              curl += coefficient * basis.curl;
            }
          }
          const auto exact = MaxwellField(point(0), point(1));
          const double exact_curl = MaxwellCurl(point(0), point(1));
          const double dx = exact[0] - discrete(0), dy = exact[1] - discrete(1);
          const double dc = exact_curl - curl;
          const double weight = weight_ref * T.Weight();
          const double error = weight * (dc * dc + MaxwellGamma * (dx * dx + dy * dy));
          report.element_error[element] += error;
          report.graph_error += error;
          report.exact_mass_energy += weight * (exact[0] * exact[0] + exact[1] * exact[1]);
          report.exact_curl_energy += weight * exact_curl * exact_curl;
        });
  }
  return report;
}

using MaxwellSparseColumn = fem::hierarchical::SparseColumn;
using MaxwellLocalElementData = fem::hierarchical::LocalOperatorContribution;

struct MaxwellLocalSystem
{
  int standard_size = 0;
  int enrichment_size = 0;
  std::vector<bool> essential;
  fem::singular::TriangleDofTopology topology;
  std::vector<MaxwellLocalElementData> elements;
};

int MaxwellUnsignedDof(int dof)
{
  return dof >= 0 ? dof : -1 - dof;
}

double MaxwellDofSign(int dof)
{
  return dof >= 0 ? 1.0 : -1.0;
}

MaxwellLocalSystem MaxwellAssembleLocalSystem(mfem::Mesh &mesh, int order)
{
  MaxwellLocalSystem system;
  mfem::ND_FECollection nd_collection(order, 2);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  mfem::H1_FECollection h1_collection(order, 2);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  system.standard_size = nd_space.GetVSize();
  const auto features = fem::singular::ExtractSerialLineFeatures(
      mesh, {MaxwellReentrantAttribute}, {{1, 1.0}});
  system.topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
  system.enrichment_size = static_cast<int>(system.topology.nd_dofs.size());
  const int combined_size = system.standard_size + system.enrichment_size;
  system.essential.assign(combined_size, false);
  mfem::Array<int> standard_essential;
  nd_space.GetBoundaryTrueDofs(standard_essential);
  for (int dof : standard_essential)
  {
    system.essential[dof] = true;
  }
  const auto numbering =
      fem::singular::BuildParallelDofNumbering(Mpi::World(), system.topology);
  mfem::Array<int> attributes(2);
  attributes[0] = MaxwellReentrantAttribute;
  attributes[1] = MaxwellOuterAttribute;
  const auto singular_essential = fem::singular::GetEssentialTriangleNDTrueDofs(
      Mpi::World(), features, system.topology, numbering, attributes);
  for (int dof : singular_essential)
  {
    system.essential[system.standard_size + dof] = true;
  }

  mfem::CurlCurlIntegrator curl_curl;
  mfem::ConstantCoefficient gamma(MaxwellGamma);
  mfem::VectorFEMassIntegrator mass(gamma);
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-9, 1.0e-9, 8};
  system.elements.resize(mesh.GetNE());
  mfem::Array<int> nd_dofs, h1_dofs;
  mfem::DenseMatrix shape;
  mfem::Vector point(3);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &data = system.elements[element];
    data.support_element = element;
    const auto &singular_dofs = system.topology.elements[element].nd;
    mfem::DofTransformation nd_transformation, h1_transformation;
    nd_space.GetElementVDofs(element, nd_dofs, nd_transformation);
    h1_space.GetElementVDofs(element, h1_dofs, h1_transformation);
    const int standard_count = nd_dofs.Size();
    for (int dof : nd_dofs)
    {
      data.dofs.push_back(MaxwellUnsignedDof(dof));
    }
    for (const auto &dof : singular_dofs)
    {
      data.dofs.push_back(system.standard_size + static_cast<int>(dof.dof));
    }
    const int local_size = static_cast<int>(data.dofs.size());
    data.matrix.SetSize(local_size);
    data.matrix = 0.0;
    data.rhs.SetSize(local_size);
    data.rhs = 0.0;
    const auto &nd_fe = *nd_space.GetFE(element);
    const auto &h1_fe = *h1_space.GetFE(element);
    auto &T = *mesh.GetElementTransformation(element);
    mfem::DenseMatrix local_curl, local_mass;
    curl_curl.AssembleElementMatrix(nd_fe, T, local_curl);
    mass.AssembleElementMatrix(nd_fe, T, local_mass);
    local_curl += local_mass;
    nd_transformation.TransformDual(local_curl);
    for (int row = 0; row < standard_count; row++)
    {
      for (int column = 0; column < standard_count; column++)
      {
        data.matrix(row, column) = MaxwellDofSign(nd_dofs[row]) *
                                   MaxwellDofSign(nd_dofs[column]) *
                                   local_curl(row, column);
      }
    }

    mfem::Vector local_rhs(standard_count);
    local_rhs = 0.0;
    shape.SetSize(standard_count, 2);
    double determinant;
    const auto grad_lambda =
        fem::singular::GetAffineTriangleBarycentricGradients(T, determinant);
    MaxwellQuadrature(
        mesh, element,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          T.Transform(ip, point);
          nd_fe.CalcPhysVShape(T, shape);
          const auto source = MaxwellSource(point(0), point(1));
          const double weight = weight_ref * T.Weight();
          for (int i = 0; i < standard_count; i++)
          {
            local_rhs(i) += weight * (shape(i, 0) * source[0] + shape(i, 1) * source[1]);
          }
          for (int i = 0; i < static_cast<int>(singular_dofs.size()); i++)
          {
            const auto basis =
                MaxwellSingularBasis(lambda, grad_lambda, singular_dofs[i].basis);
            data.rhs(standard_count + i) +=
                weight * (basis.value[0] * source[0] + basis.value[1] * source[1]);
          }
        });
    nd_transformation.TransformDual(local_rhs);
    for (int i = 0; i < standard_count; i++)
    {
      data.rhs(i) = MaxwellDofSign(nd_dofs[i]) * local_rhs(i);
    }

    if (!singular_dofs.empty())
    {
      const auto singular = fem::singular::AssembleTriangleElementEnrichmentMatrices(
          system.topology.elements[element], T, options);
      auto coupling = fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
          system.topology.elements[element], h1_fe, nd_fe, T, options);
      fem::singular::ApplyStandardDofTransformations(h1_transformation, nd_transformation,
                                                     coupling);
      for (int row = 0; row < standard_count; row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(row, standard_count + column) =
              MaxwellDofSign(nd_dofs[row]) *
              (coupling.nd_curl_curl_standard_enrichment(row, column) +
               MaxwellGamma * coupling.nd_mass_standard_enrichment(row, column));
          data.matrix(standard_count + column, row) =
              MaxwellDofSign(nd_dofs[row]) *
              (coupling.nd_curl_curl_enrichment_standard(column, row) +
               MaxwellGamma * coupling.nd_mass_enrichment_standard(column, row));
        }
      }
      for (int row = 0; row < static_cast<int>(singular_dofs.size()); row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(standard_count + row, standard_count + column) =
              singular.nd_curl_curl(row, column) +
              MaxwellGamma * singular.nd_mass(row, column);
        }
      }
    }
  }
  return system;
}

std::vector<MaxwellSparseColumn> MaxwellBuildSparseTransfer(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, double &maximum_consistency_error,
    int *signed_coarse_dofs = nullptr, int *signed_fine_dofs = nullptr,
    int *nonidentity_transformations = nullptr)
{
  std::vector<std::map<int, double>> entries(coarse_space.GetVSize());
  mfem::IsoparametricTransformation identity;
  identity.SetIdentityTransformation(mfem::Geometry::TRIANGLE);
  mfem::DenseMatrix local_transfer;
  mfem::Array<int> coarse_dofs, fine_dofs;
  maximum_consistency_error = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::DofTransformation coarse_transformation, fine_transformation;
    coarse_space.GetElementVDofs(element, coarse_dofs, coarse_transformation);
    fine_space.GetElementVDofs(element, fine_dofs, fine_transformation);
    if (signed_coarse_dofs)
    {
      *signed_coarse_dofs += std::count_if(coarse_dofs.begin(), coarse_dofs.end(),
                                           [](int dof) { return dof < 0; });
    }
    if (signed_fine_dofs)
    {
      *signed_fine_dofs += std::count_if(fine_dofs.begin(), fine_dofs.end(),
                                         [](int dof) { return dof < 0; });
    }
    if (nonidentity_transformations)
    {
      *nonidentity_transformations += !coarse_transformation.IsIdentity();
      *nonidentity_transformations += !fine_transformation.IsIdentity();
    }
    fine_space.GetFE(element)->GetTransferMatrix(*coarse_space.GetFE(element), identity,
                                                 local_transfer);
    std::set<int> element_coarse;
    for (int dof : coarse_dofs)
    {
      element_coarse.insert(MaxwellUnsignedDof(dof));
    }
    for (int global_coarse : element_coarse)
    {
      mfem::Vector local_coarse(coarse_dofs.Size());
      local_coarse = 0.0;
      for (int i = 0; i < coarse_dofs.Size(); i++)
      {
        if (MaxwellUnsignedDof(coarse_dofs[i]) == global_coarse)
        {
          local_coarse(i) = MaxwellDofSign(coarse_dofs[i]);
        }
      }
      coarse_transformation.InvTransformPrimal(local_coarse);
      mfem::Vector local_fine(fine_dofs.Size());
      local_transfer.Mult(local_coarse, local_fine);
      fine_transformation.TransformPrimal(local_fine);
      auto &column = entries[global_coarse];
      for (int i = 0; i < fine_dofs.Size(); i++)
      {
        const double value = MaxwellDofSign(fine_dofs[i]) * local_fine(i);
        if (std::abs(value) <= 1.0e-13)
        {
          continue;
        }
        const int row = MaxwellUnsignedDof(fine_dofs[i]);
        const auto [found, inserted] = column.emplace(row, value);
        if (!inserted)
        {
          maximum_consistency_error =
              std::max(maximum_consistency_error, std::abs(found->second - value));
        }
      }
    }
  }
  std::vector<MaxwellSparseColumn> result(entries.size());
  for (std::size_t column = 0; column < entries.size(); column++)
  {
    for (const auto &[row, value] : entries[column])
    {
      result[column].dofs.push_back(row);
      result[column].values.push_back(value);
    }
  }
  return result;
}

struct MaxwellPatchReport
{
  std::vector<double> indicator;
  double energy = 0.0;
  double work = 0.0;
  double transfer_consistency_error = 0.0;
  double transfer_reference_error = 0.0;
  double maximum_patch_residual = 0.0;
  double maximum_patch_condition = 0.0;
  int edge_patches = 0;
  int interior_patches = 0;
  int owned_modes = 0;
  int maximum_support_elements = 0;
  int maximum_patch_dimension = 0;
  int maximum_element_overlap = 0;
  int covered_gradient_guests = 0;
  int covered_rotational_guests = 0;
  int unique_gradient_guests = 0;
  int unique_rotational_guests = 0;
  int unique_coarse_guests = 0;
  int signed_coarse_dofs = 0;
  int signed_fine_dofs = 0;
  int nonidentity_transformations = 0;
};

MaxwellPatchReport MaxwellBuildLocalPatches(mfem::Mesh &mesh, int order,
                                            const MaxwellSolveReport &coarse,
                                            bool include_rotational_guests = true)
{
  MaxwellPatchReport report;
  auto fine = MaxwellAssembleLocalSystem(mesh, order + 1);
  mfem::ND_FECollection coarse_collection(order, 2), fine_collection(order + 1, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  const int enrichment = fine.enrichment_size;
  const int combined_size = fine_standard + enrichment;
  REQUIRE(coarse.enrichment_dofs == enrichment);
  auto injection = MaxwellBuildSparseTransfer(
      mesh, coarse_space, fine_space, report.transfer_consistency_error,
      &report.signed_coarse_dofs, &report.signed_fine_dofs,
      &report.nonidentity_transformations);

  // Independent orientation/reference check against MFEM's p-transfer.
  mfem::Vector probe(coarse_standard), sparse_image(fine_standard),
      reference(fine_standard);
  for (int i = 0; i < coarse_standard; i++)
  {
    probe(i) = std::sin(0.37 * i + 0.2);
  }
  sparse_image = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection[column].dofs.size(); entry++)
    {
      sparse_image(injection[column].dofs[entry]) +=
          probe(column) * injection[column].values[entry];
    }
  }
  mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
  transfer.Mult(probe, reference);
  sparse_image -= reference;
  report.transfer_reference_error = sparse_image.Normlinf();

  mfem::Vector injected(combined_size);
  injected = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection[column].dofs.size(); entry++)
    {
      injected(injection[column].dofs[entry]) +=
          coarse.standard_coefficients(column) * injection[column].values[entry];
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    injected(fine_standard + dof) = coarse.enrichment_coefficients(dof);
  }
  mfem::Vector residual = fem::hierarchical::AssembleResidual(combined_size, fine.elements,
                                                              injected, fine.essential);
  std::vector<std::set<int>> dof_elements(combined_size);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    for (int dof : fine.elements[element].dofs)
    {
      dof_elements[dof].insert(element);
    }
  }

  const auto entity_complement =
      [&](mfem::Array<int> fine_entity, mfem::Array<int> coarse_entity)
  {
    mfem::FiniteElementSpace::AdjustVDofs(fine_entity);
    mfem::FiniteElementSpace::AdjustVDofs(coarse_entity);
    std::vector<int> rows, columns;
    for (int dof : fine_entity)
    {
      if (!fine.essential[dof])
      {
        rows.push_back(dof);
      }
    }
    for (int dof : coarse_entity)
    {
      if (!coarse.essential_mask[dof])
      {
        columns.push_back(dof);
      }
    }
    std::vector<mfem::Vector> range;
    for (int column : columns)
    {
      mfem::Vector vector(static_cast<int>(rows.size()));
      vector = 0.0;
      for (int i = 0; i < static_cast<int>(rows.size()); i++)
      {
        for (std::size_t entry = 0; entry < injection[column].dofs.size(); entry++)
        {
          if (injection[column].dofs[entry] == rows[i])
          {
            vector(i) = injection[column].values[entry];
          }
        }
      }
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-12)
      {
        vector /= norm;
        range.push_back(vector);
      }
    }
    std::vector<mfem::Vector> complement;
    for (int direction = 0; direction < static_cast<int>(rows.size()); direction++)
    {
      mfem::Vector vector(static_cast<int>(rows.size()));
      vector = 0.0;
      vector(direction) = 1.0;
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
        for (const auto &basis : complement)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-10)
      {
        vector /= norm;
        complement.push_back(vector);
      }
    }
    std::vector<MaxwellSparseColumn> result;
    for (const auto &local : complement)
    {
      MaxwellSparseColumn column;
      for (int i = 0; i < static_cast<int>(rows.size()); i++)
      {
        if (std::abs(local(i)) > 1.0e-13)
        {
          column.dofs.push_back(rows[i]);
          column.values.push_back(local(i));
        }
      }
      result.push_back(std::move(column));
    }
    return result;
  };

  struct Patch
  {
    int owned = 0;
    std::vector<MaxwellSparseColumn> basis;
    std::vector<int> coarse_guests;
    std::vector<int> singular_guests;
    mfem::Vector coefficients;
    mfem::DenseMatrix restricted;
    std::set<int> support;
  };
  std::vector<Patch> patches;
  std::vector<int> overlap(mesh.GetNE(), 0);
  mfem::Array<int> element_dofs;
  const auto add_patch = [&](std::vector<MaxwellSparseColumn> owned,
                             const std::vector<int> &owner_elements, bool edge)
  {
    if (owned.empty())
    {
      return;
    }
    Patch patch;
    patch.owned = static_cast<int>(owned.size());
    patch.basis = std::move(owned);
    report.owned_modes += patch.owned;
    report.edge_patches += edge;
    report.interior_patches += !edge;
    std::set<int> coarse_guests, singular_guests;
    for (int element : owner_elements)
    {
      mfem::DofTransformation transformation;
      coarse_space.GetElementVDofs(element, element_dofs, transformation);
      for (int dof : element_dofs)
      {
        dof = MaxwellUnsignedDof(dof);
        if (!coarse.essential_mask[dof])
        {
          coarse_guests.insert(dof);
        }
      }
      for (const auto &dof : fine.topology.elements[element].nd)
      {
        const bool rotational = fine.topology.nd_dofs[dof.dof].family ==
                                fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
        if (!fine.essential[fine_standard + static_cast<int>(dof.dof)] &&
            (include_rotational_guests || !rotational))
        {
          singular_guests.insert(static_cast<int>(dof.dof));
        }
      }
    }
    for (int dof : coarse_guests)
    {
      patch.coarse_guests.push_back(dof);
      patch.basis.push_back(injection[dof]);
    }
    for (int dof : singular_guests)
    {
      MaxwellSparseColumn column;
      column.dofs.push_back(fine_standard + dof);
      column.values.push_back(1.0);
      patch.singular_guests.push_back(dof);
      patch.basis.push_back(column);
      report.covered_gradient_guests +=
          fine.topology.nd_dofs[dof].family ==
          fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
      report.covered_rotational_guests +=
          fine.topology.nd_dofs[dof].family ==
          fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
    }
    const int dimension = static_cast<int>(patch.basis.size());
    report.maximum_patch_dimension = std::max(report.maximum_patch_dimension, dimension);
    std::map<int, std::vector<std::pair<int, double>>> row_basis;
    bool touches_essential = false;
    for (int basis = 0; basis < dimension; basis++)
    {
      for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
      {
        const int dof = patch.basis[basis].dofs[entry];
        touches_essential = touches_essential || fine.essential[dof];
        row_basis[dof].push_back({basis, patch.basis[basis].values[entry]});
        patch.support.insert(dof_elements[dof].begin(), dof_elements[dof].end());
      }
    }
    REQUIRE_FALSE(touches_essential);
    report.maximum_support_elements =
        std::max(report.maximum_support_elements, static_cast<int>(patch.support.size()));
    for (int element : patch.support)
    {
      overlap[element]++;
    }
    mfem::DenseMatrix restricted;
    mfem::Vector patch_rhs;
    fem::hierarchical::AssembleRestrictedOperator(fine.elements, patch.support, patch.basis,
                                                  residual, restricted, patch_rhs);
    patch.restricted = restricted;
    mfem::DenseMatrixInverse inverse(restricted, true);
    mfem::DenseMatrix inverse_matrix;
    inverse.GetInverseMatrix(inverse_matrix);
    report.maximum_patch_condition = std::max(report.maximum_patch_condition,
                                              restricted.FNorm() * inverse_matrix.FNorm());
    patch.coefficients.SetSize(dimension);
    inverse.Mult(patch_rhs, patch.coefficients);
    mfem::Vector solve_residual(dimension);
    restricted.Mult(patch.coefficients, solve_residual);
    solve_residual -= patch_rhs;
    report.maximum_patch_residual =
        std::max(report.maximum_patch_residual,
                 solve_residual.Norml2() / std::max(patch_rhs.Norml2(), 1.0e-30));
    patches.push_back(std::move(patch));
  };

  std::vector<std::vector<int>> edge_elements(mesh.GetNEdges());
  mfem::Array<int> edges, orientations;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementEdges(element, edges, orientations);
    for (int edge : edges)
    {
      edge_elements[edge].push_back(element);
    }
  }
  for (int edge = 0; edge < mesh.GetNEdges(); edge++)
  {
    mfem::Array<int> coarse_entity, fine_entity;
    coarse_space.GetEdgeDofs(edge, coarse_entity);
    fine_space.GetEdgeDofs(edge, fine_entity);
    add_patch(entity_complement(fine_entity, coarse_entity), edge_elements[edge], true);
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::Array<int> coarse_entity, fine_entity;
    coarse_space.GetElementInteriorDofs(element, coarse_entity);
    fine_space.GetElementInteriorDofs(element, fine_entity);
    add_patch(entity_complement(fine_entity, coarse_entity), {element}, false);
  }
  int coarse_free = 0, fine_free = 0;
  for (int dof = 0; dof < coarse_standard; dof++)
  {
    coarse_free += !coarse.essential_mask[dof];
  }
  for (int dof = 0; dof < fine_standard; dof++)
  {
    fine_free += !fine.essential[dof];
  }
  REQUIRE(report.owned_modes == fine_free - coarse_free);
  report.maximum_element_overlap = *std::max_element(overlap.begin(), overlap.end());

  mfem::Vector raw(combined_size);
  raw = 0.0;
  std::vector<double> coarse_sum(coarse_standard, 0.0), singular_sum(enrichment, 0.0);
  std::vector<int> coarse_count(coarse_standard, 0), singular_count(enrichment, 0);
  for (const auto &patch : patches)
  {
    for (int basis = 0; basis < patch.owned; basis++)
    {
      for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
      {
        raw(patch.basis[basis].dofs[entry]) +=
            patch.coefficients(basis) * patch.basis[basis].values[entry];
      }
    }
    int coefficient = patch.owned;
    for (int dof : patch.coarse_guests)
    {
      coarse_sum[dof] += patch.coefficients(coefficient++);
      coarse_count[dof]++;
    }
    for (int dof : patch.singular_guests)
    {
      singular_sum[dof] += patch.coefficients(coefficient++);
      singular_count[dof]++;
    }
  }
  for (int dof = 0; dof < coarse_standard; dof++)
  {
    if (!coarse.essential_mask[dof])
    {
      REQUIRE(coarse_count[dof] > 0);
      report.unique_coarse_guests++;
    }
    if (coarse_count[dof] > 0)
    {
      const double coefficient = coarse_sum[dof] / coarse_count[dof];
      for (std::size_t entry = 0; entry < injection[dof].dofs.size(); entry++)
      {
        raw(injection[dof].dofs[entry]) += coefficient * injection[dof].values[entry];
      }
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    const bool rotational = fine.topology.nd_dofs[dof].family ==
                            fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
    if (!fine.essential[fine_standard + dof] && (include_rotational_guests || !rotational))
    {
      REQUIRE(singular_count[dof] > 0);
      report.unique_gradient_guests += fine.topology.nd_dofs[dof].family ==
                                       fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
      report.unique_rotational_guests +=
          fine.topology.nd_dofs[dof].family ==
          fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
    }
    if (singular_count[dof] > 0)
    {
      raw(fine_standard + dof) += singular_sum[dof] / singular_count[dof];
    }
  }
  double raw_energy = fem::hierarchical::Energy(fine.elements, raw);
  double raw_work = mfem::InnerProduct(raw, residual);
  const double alpha = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
  raw *= alpha;

  // Additional additive-Schwarz defect-correction sweeps. These remain patch-local and are
  // especially important for p1, where one independent patch response poorly approximates
  // the long-range coarse Schur relaxation.
  for (int sweep = 1; sweep < 4; sweep++)
  {
    mfem::Vector current_residual(residual);
    for (const auto &data : fine.elements)
    {
      mfem::Vector local(static_cast<int>(data.dofs.size()));
      for (int i = 0; i < local.Size(); i++)
      {
        local(i) = raw(data.dofs[i]);
      }
      mfem::Vector action(local.Size());
      data.matrix.Mult(local, action);
      for (int i = 0; i < local.Size(); i++)
      {
        current_residual(data.dofs[i]) -= action(i);
      }
    }
    for (auto &patch : patches)
    {
      mfem::Vector patch_rhs(static_cast<int>(patch.basis.size()));
      for (int basis = 0; basis < static_cast<int>(patch.basis.size()); basis++)
      {
        patch_rhs(basis) = 0.0;
        for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
        {
          patch_rhs(basis) += patch.basis[basis].values[entry] *
                              current_residual(patch.basis[basis].dofs[entry]);
        }
      }
      mfem::DenseMatrixInverse inverse(patch.restricted, true);
      inverse.Mult(patch_rhs, patch.coefficients);
    }
    mfem::Vector delta(combined_size);
    delta = 0.0;
    std::fill(coarse_sum.begin(), coarse_sum.end(), 0.0);
    std::fill(singular_sum.begin(), singular_sum.end(), 0.0);
    std::fill(coarse_count.begin(), coarse_count.end(), 0);
    std::fill(singular_count.begin(), singular_count.end(), 0);
    for (const auto &patch : patches)
    {
      for (int basis = 0; basis < patch.owned; basis++)
      {
        for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
        {
          delta(patch.basis[basis].dofs[entry]) +=
              patch.coefficients(basis) * patch.basis[basis].values[entry];
        }
      }
      int coefficient = patch.owned;
      for (int dof : patch.coarse_guests)
      {
        coarse_sum[dof] += patch.coefficients(coefficient++);
        coarse_count[dof]++;
      }
      for (int dof : patch.singular_guests)
      {
        singular_sum[dof] += patch.coefficients(coefficient++);
        singular_count[dof]++;
      }
    }
    for (int dof = 0; dof < coarse_standard; dof++)
    {
      if (coarse_count[dof] > 0)
      {
        const double coefficient = coarse_sum[dof] / coarse_count[dof];
        for (std::size_t entry = 0; entry < injection[dof].dofs.size(); entry++)
        {
          delta(injection[dof].dofs[entry]) += coefficient * injection[dof].values[entry];
        }
      }
    }
    for (int dof = 0; dof < enrichment; dof++)
    {
      if (singular_count[dof] > 0)
      {
        delta(fine_standard + dof) += singular_sum[dof] / singular_count[dof];
      }
    }
    const double delta_energy = fem::hierarchical::Energy(fine.elements, delta);
    const double delta_work = mfem::InnerProduct(delta, current_residual);
    if (delta_energy > 0.0)
    {
      raw.Add(delta_work / delta_energy, delta);
    }
  }

  // A final scalar projection restores the energy/work identity without changing ranking.
  raw_energy = fem::hierarchical::Energy(fine.elements, raw);
  raw_work = mfem::InnerProduct(raw, residual);
  const double final_scale = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
  raw *= final_scale;
  report.energy = final_scale * final_scale * raw_energy;
  report.work = final_scale * raw_work;
  report.indicator.assign(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &data = fine.elements[element];
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = raw(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    report.indicator[element] = mfem::InnerProduct(local, action);
  }
  return report;
}

std::vector<double> MaxwellGlobalHierarchyIndicator(mfem::Mesh &mesh, int order,
                                                    const MaxwellSolveReport &coarse)
{
  const auto fine_solution = MaxwellSolve(mesh, order + 1, true);
  const auto fine = MaxwellAssembleLocalSystem(mesh, order + 1);
  mfem::ND_FECollection coarse_collection(order, 2), fine_collection(order + 1, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  double consistency;
  const auto injection =
      MaxwellBuildSparseTransfer(mesh, coarse_space, fine_space, consistency);
  REQUIRE(consistency < 1.0e-12);
  mfem::Vector correction(fine.standard_size + fine.enrichment_size);
  correction = 0.0;
  for (int i = 0; i < fine.standard_size; i++)
  {
    correction(i) = fine_solution.standard_coefficients(i);
  }
  for (int column = 0; column < coarse_space.GetVSize(); column++)
  {
    for (std::size_t entry = 0; entry < injection[column].dofs.size(); entry++)
    {
      correction(injection[column].dofs[entry]) -=
          coarse.standard_coefficients(column) * injection[column].values[entry];
    }
  }
  for (int dof = 0; dof < fine.enrichment_size; dof++)
  {
    correction(fine.standard_size + dof) =
        fine_solution.enrichment_coefficients(dof) - coarse.enrichment_coefficients(dof);
  }
  std::vector<double> indicator(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &data = fine.elements[element];
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = correction(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    indicator[element] = mfem::InnerProduct(local, action);
  }
  return indicator;
}

std::vector<std::size_t> MaxwellDorfler(const std::vector<double> &indicator, double theta)
{
  std::vector<std::size_t> indices(indicator.size());
  std::iota(indices.begin(), indices.end(), std::size_t{0});
  std::sort(indices.begin(), indices.end(),
            [&](std::size_t a, std::size_t b) { return indicator[a] > indicator[b]; });
  const double target = theta * std::accumulate(indicator.begin(), indicator.end(), 0.0);
  double sum = 0.0;
  std::size_t count = 0;
  while (count < indices.size() && sum < target)
  {
    sum += indicator[indices[count++]];
  }
  indices.resize(count);
  return indices;
}

double MaxwellSumOn(const std::vector<double> &value,
                    const std::vector<std::size_t> &indices)
{
  double sum = 0.0;
  for (std::size_t index : indices)
  {
    sum += value[index];
  }
  return sum;
}

}  // namespace

TEST_CASE("Manufactured Maxwell wedge source and tangential trace are consistent",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  for (const auto &point : std::vector<std::array<double, 2>>{
           {-1.0, 0.2}, {1.0, 0.4}, {-0.3, 1.0}, {-0.5, -1.0}, {0.0, -0.5}, {0.5, 0.0}})
  {
    const auto field = MaxwellField(point[0], point[1]);
    std::array<double, 2> tangent{};
    if (std::abs(std::abs(point[0]) - 1.0) < 1.0e-12 || std::abs(point[0]) < 1.0e-12)
    {
      tangent = {0.0, 1.0};
    }
    else
    {
      tangent = {1.0, 0.0};
    }
    CHECK_THAT(field[0] * tangent[0] + field[1] * tangent[1], WithinAbs(0.0, 1.0e-12));
  }
  CHECK(std::abs(MaxwellCurl(-0.4, 0.3)) > 1.0e-3);
  constexpr double h = 1.0e-6;
  for (const auto &point :
       std::vector<std::array<double, 2>>{{-0.4, 0.3}, {-0.7, -0.4}, {0.4, 0.6}})
  {
    const auto field = MaxwellField(point[0], point[1]);
    const double dc_dx =
        (MaxwellCurl(point[0] + h, point[1]) - MaxwellCurl(point[0] - h, point[1])) /
        (2.0 * h);
    const double dc_dy =
        (MaxwellCurl(point[0], point[1] + h) - MaxwellCurl(point[0], point[1] - h)) /
        (2.0 * h);
    const auto source = MaxwellSource(point[0], point[1]);
    CHECK_THAT(source[0], WithinAbs(dc_dy + MaxwellGamma * field[0], 1.0e-9));
    CHECK_THAT(source[1], WithinAbs(-dc_dx + MaxwellGamma * field[1], 1.0e-9));
  }
}

TEST_CASE("Manufactured Maxwell wedge gradient and rotational blocks round trip",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellLShapeMesh(2);
  const auto report = MaxwellSolve(mesh, order, true, true, true);
  CAPTURE(order, report.residual, report.gradient_roundtrip_error,
          report.rotational_roundtrip_error);
  CHECK(report.residual < 1.0e-10);
  CHECK(report.gradient_roundtrip_error < 1.0e-10);
  CHECK(report.rotational_roundtrip_error < 1.0e-10);
}

TEST_CASE("Manufactured Maxwell wedge sparse ND patches preserve orientation and families",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellLShapeMesh(2);
  const auto coarse = MaxwellSolve(mesh, order, true);
  const auto local = MaxwellBuildLocalPatches(mesh, order, coarse);
  const auto gradient_guest_local = MaxwellBuildLocalPatches(mesh, order, coarse, false);
  CAPTURE(order, local.energy, local.work, local.transfer_consistency_error,
          local.transfer_reference_error, local.maximum_patch_residual,
          local.maximum_patch_condition, local.edge_patches, local.interior_patches,
          local.owned_modes, local.maximum_support_elements, local.maximum_patch_dimension,
          local.maximum_element_overlap, local.covered_gradient_guests,
          local.covered_rotational_guests, local.unique_gradient_guests,
          local.unique_rotational_guests, local.unique_coarse_guests,
          local.signed_coarse_dofs, local.signed_fine_dofs,
          local.nonidentity_transformations);
  CHECK(local.transfer_consistency_error < 1.0e-12);
  CHECK(local.transfer_reference_error < 1.0e-12);
  CHECK(local.maximum_patch_residual < 1.0e-10);
  CHECK(local.maximum_patch_condition < 1.0e7);
  CHECK(local.edge_patches > 0);
  CHECK(local.interior_patches > 0);
  CHECK(local.owned_modes > 0);
  CHECK(local.maximum_support_elements < 30);
  CHECK(local.maximum_patch_dimension < 30);
  CHECK(local.maximum_element_overlap < 100);
  CHECK(local.covered_gradient_guests > 0);
  CHECK(local.covered_rotational_guests > 0);
  CHECK(local.unique_gradient_guests > 0);
  CHECK(local.unique_rotational_guests > 0);
  CHECK(local.unique_coarse_guests > 0);
  CHECK(local.signed_coarse_dofs > 0);
  CHECK(local.signed_fine_dofs > 0);
  CHECK(local.energy > 0.0);
  CHECK_THAT(local.energy, WithinRel(local.work, 1.0e-10));
  CHECK_THAT(std::accumulate(local.indicator.begin(), local.indicator.end(), 0.0),
             WithinRel(local.energy, 1.0e-12));

  const auto marked = MaxwellDorfler(local.indicator, 0.5);
  const auto oracle = MaxwellDorfler(coarse.element_error, 0.5);
  const auto top_k = [](const std::vector<double> &value, std::size_t count)
  {
    std::vector<std::size_t> indices(value.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    std::partial_sort(indices.begin(), indices.begin() + std::min(count, indices.size()),
                      indices.end(),
                      [&](std::size_t a, std::size_t b) { return value[a] > value[b]; });
    indices.resize(std::min(count, indices.size()));
    return indices;
  };
  const double rank =
      MaxwellSumOn(coarse.element_error, marked) /
      MaxwellSumOn(coarse.element_error, top_k(coarse.element_error, marked.size()));
  const double extend =
      MaxwellSumOn(coarse.element_error, top_k(local.indicator, oracle.size())) /
      MaxwellSumOn(coarse.element_error, oracle);
  const auto gradient_marked = MaxwellDorfler(gradient_guest_local.indicator, 0.5);
  const double gradient_rank =
      MaxwellSumOn(coarse.element_error, gradient_marked) /
      MaxwellSumOn(coarse.element_error,
                   top_k(coarse.element_error, gradient_marked.size()));
  const double gradient_extend =
      MaxwellSumOn(coarse.element_error,
                   top_k(gradient_guest_local.indicator, oracle.size())) /
      MaxwellSumOn(coarse.element_error, oracle);
  const auto global_indicator = MaxwellGlobalHierarchyIndicator(mesh, order, coarse);
  const auto global_marked = MaxwellDorfler(global_indicator, 0.5);
  const double global_rank =
      MaxwellSumOn(coarse.element_error, global_marked) /
      MaxwellSumOn(coarse.element_error, top_k(coarse.element_error, global_marked.size()));
  const double global_extend =
      MaxwellSumOn(coarse.element_error, top_k(global_indicator, oracle.size())) /
      MaxwellSumOn(coarse.element_error, oracle);
  CAPTURE(rank, extend, gradient_rank, gradient_extend, global_rank, global_extend,
          marked.size(), oracle.size());
  if (order == 1)
  {
    // Multiple local defect-correction sweeps close the poor first-iterate p1 response
    // while preserving bounded patch support and orientation-aware transfer.
    CHECK(std::max(rank, gradient_rank) > 0.8);
    CHECK(std::max(extend, gradient_extend) > 0.8);
    CHECK(global_rank - std::max(rank, gradient_rank) < 0.2);
    CHECK(global_extend - std::max(extend, gradient_extend) < 0.2);
  }
  else
  {
    CHECK(rank > 0.85);
    CHECK(extend > 0.8);
    CHECK(rank >= global_rank - 0.15);
    CHECK(extend >= global_extend - 0.15);
  }

  // Independent element-local/global graph-operator oracle, including signed ND rows and
  // both singular families.
  const auto local_system = MaxwellAssembleLocalSystem(mesh, order + 1);
  const auto global_system = MaxwellSolve(mesh, order + 1, true);
  const int combined_size = local_system.standard_size + local_system.enrichment_size;
  mfem::DenseMatrix scattered(combined_size), assembled(combined_size);
  scattered = 0.0;
  assembled = 0.0;
  mfem::Vector scattered_rhs(combined_size);
  scattered_rhs = 0.0;
  for (const auto &data : local_system.elements)
  {
    for (int row = 0; row < static_cast<int>(data.dofs.size()); row++)
    {
      scattered_rhs(data.dofs[row]) += data.rhs(row);
      for (int column = 0; column < static_cast<int>(data.dofs.size()); column++)
      {
        scattered(data.dofs[row], data.dofs[column]) += data.matrix(row, column);
      }
    }
  }
  for (int row = 0; row < combined_size; row++)
  {
    for (int entry = global_system.combined_matrix->GetI()[row];
         entry < global_system.combined_matrix->GetI()[row + 1]; entry++)
    {
      assembled(row, global_system.combined_matrix->GetJ()[entry]) =
          global_system.combined_matrix->GetData()[entry];
    }
  }
  double matrix_difference = 0.0, rhs_difference = 0.0;
  for (int row = 0; row < combined_size; row++)
  {
    REQUIRE(local_system.essential[row] == global_system.essential_mask[row]);
    rhs_difference = std::max(
        rhs_difference, std::abs(scattered_rhs(row) - global_system.combined_rhs(row)));
    for (int column = 0; column < combined_size; column++)
    {
      matrix_difference = std::max(
          matrix_difference, std::abs(scattered(row, column) - assembled(row, column)));
    }
  }
  CAPTURE(matrix_difference, rhs_difference);
  CHECK(matrix_difference < 1.0e-10);
  CHECK(rhs_difference < 1.0e-10);
}

TEST_CASE("Manufactured Maxwell wedge sparse ND support remains local",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  std::vector<int> elements, support, overlap, dimension;
  for (int n : {2, 3})
  {
    auto mesh = MaxwellLShapeMesh(n);
    const auto coarse = MaxwellSolve(mesh, 2, true);
    const auto local = MaxwellBuildLocalPatches(mesh, 2, coarse);
    elements.push_back(mesh.GetNE());
    support.push_back(local.maximum_support_elements);
    overlap.push_back(local.maximum_element_overlap);
    dimension.push_back(local.maximum_patch_dimension);
  }
  CAPTURE(elements, support, overlap, dimension);
  CHECK(*std::max_element(support.begin(), support.end()) < 20);
  CHECK(*std::max_element(overlap.begin(), overlap.end()) < 30);
  CHECK(*std::max_element(dimension.begin(), dimension.end()) < 30);
}

TEST_CASE("Manufactured Maxwell wedge sparse ND AMR beats uniform",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 3;
  auto mesh = MaxwellLShapeMesh(2);
  std::vector<long long> dofs;
  std::vector<double> errors;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto solve = MaxwellSolve(mesh, order, true);
    const auto local = MaxwellBuildLocalPatches(mesh, order, solve);
    dofs.push_back(solve.dofs);
    errors.push_back(solve.graph_error);
    if (pass + 1 < passes)
    {
      mfem::Array<mfem::Refinement> refinements;
      for (std::size_t element : MaxwellDorfler(local.indicator, 0.5))
      {
        refinements.Append(mfem::Refinement(static_cast<int>(element)));
      }
      mesh.GeneralRefinement(refinements, -1, 1);
      REQUIRE(mesh.Conforming());
    }
  }
  std::vector<std::pair<long long, double>> uniform;
  for (int n : {2, 3, 4, 6})
  {
    auto uniform_mesh = MaxwellLShapeMesh(n);
    const auto solve = MaxwellSolve(uniform_mesh, order, true);
    uniform.emplace_back(solve.dofs, solve.graph_error);
  }
  const auto envelope = [&](long long at_dofs)
  {
    for (std::size_t i = 0; i + 1 < uniform.size(); i++)
    {
      if (uniform[i].first <= at_dofs && at_dofs <= uniform[i + 1].first)
      {
        const double t = (std::log(static_cast<double>(at_dofs)) -
                          std::log(static_cast<double>(uniform[i].first))) /
                         (std::log(static_cast<double>(uniform[i + 1].first)) -
                          std::log(static_cast<double>(uniform[i].first)));
        return std::exp(std::log(uniform[i].second) + t * (std::log(uniform[i + 1].second) -
                                                           std::log(uniform[i].second)));
      }
    }
    return std::numeric_limits<double>::quiet_NaN();
  };
  std::vector<double> ratios;
  for (std::size_t pass = 0; pass < errors.size(); pass++)
  {
    const double reference = envelope(dofs[pass]);
    REQUIRE(std::isfinite(reference));
    ratios.push_back(errors[pass] / reference);
    if (pass > 0)
    {
      CHECK(errors[pass] < errors[pass - 1]);
    }
  }
  CAPTURE(dofs, errors, uniform, ratios);
  CHECK(ratios.back() < 1.0);
}

TEST_CASE("Manufactured Maxwell wedge singular enrichment improves the graph error",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  auto mesh = MaxwellLShapeMesh(2);
  const auto standard = MaxwellSolve(mesh, order, false);
  const auto gradient_only = MaxwellSolve(mesh, order, true, false);
  const auto enriched = MaxwellSolve(mesh, order, true);
  CAPTURE(standard.dofs, enriched.dofs, enriched.enrichment_dofs,
          enriched.free_enrichment_dofs, enriched.free_gradient_dofs,
          enriched.free_rotational_dofs, standard.residual, gradient_only.residual,
          enriched.residual, standard.graph_error, gradient_only.graph_error,
          enriched.graph_error, enriched.exact_mass_energy, enriched.exact_curl_energy);
  CHECK(standard.residual < 1.0e-10);
  CHECK(gradient_only.residual < 1.0e-10);
  CHECK(enriched.residual < 1.0e-10);
  CHECK(enriched.enrichment_dofs > 0);
  CHECK(enriched.free_gradient_dofs > 0);
  CHECK(enriched.free_rotational_dofs > 0);
  CHECK(gradient_only.free_rotational_dofs == 0);
  CHECK(enriched.exact_mass_energy > 1.0e-3);
  CHECK(enriched.exact_curl_energy > 1.0e-3);
  // The exact singular component is a gradient. This ablation proves that family supplies
  // the principal improvement while retained rotational modes remain active and can relax
  // the nonzero-curl remainder.
  CHECK(gradient_only.graph_error < 0.9 * standard.graph_error);
  CHECK(enriched.graph_error <= gradient_only.graph_error * (1.0 + 1.0e-10));
}

}  // namespace palace
