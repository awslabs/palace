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
constexpr double MaxwellOmega2 = 0.36;
constexpr double MaxwellBoundaryConductance = 0.75;
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

// Indefinite frequency-domain source: f = curl curl E - omega^2 E for the same
// manufactured field.
std::array<double, 2> MaxwellDrivenSource(double x, double y)
{
  const auto field = MaxwellField(x, y);
  const double curl_curl_x = 6.0 * MaxwellBeta * x * y * (1.0 - x * x);
  const double curl_curl_y = MaxwellBeta * (1.0 - 3.0 * x * x) * (1.0 - 3.0 * y * y);
  return {curl_curl_x - MaxwellOmega2 * field[0], curl_curl_y - MaxwellOmega2 * field[1]};
}

// Manufactured impedance data on the outer boundary: the exact field satisfies
// curl E = -g (E . t) + h with E . t = 0 there, so h is the exact scalar curl trace. It is
// nonzero on the y = +-1 segments, which makes the facet load an active oracle for the
// boundary-contribution orientation.
double MaxwellFacetData(double x, double y)
{
  return MaxwellCurl(x, y);
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

// Measure the graph-norm error of a combined standard-plus-singular solution against the
// manufactured field, with edge-aligned Duffy quadrature at the corner. The solution vector
// uses the [standard, enrichment] combined layout.
void MaxwellMeasureGraphError(mfem::Mesh &mesh, mfem::FiniteElementSpace &nd_space,
                              const fem::singular::TriangleDofTopology *topology,
                              const mfem::Vector &solution, MaxwellSolveReport &report)
{
  const int standard_size = nd_space.GetVSize();
  mfem::GridFunction standard_solution(&nd_space);
  for (int i = 0; i < standard_size; i++)
  {
    standard_solution(i) = solution(i);
  }
  mfem::Vector point(3), discrete(2), discrete_curl;
  report.graph_error = 0.0;
  report.exact_mass_energy = 0.0;
  report.exact_curl_energy = 0.0;
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
          if (topology)
          {
            for (const auto &dof : topology->elements[element].nd)
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
}

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

  MaxwellMeasureGraphError(mesh, nd_space, enriched ? &topology : nullptr, solution,
                           report);
  return report;
}

using MaxwellSparseColumn = fem::hierarchical::SparseColumn;
using MaxwellLocalElementData = fem::hierarchical::LocalOperatorContribution;

// Assembly options for the local uneliminated system. The defaults reproduce the coercive
// manufactured fixture. Driven variants flip the mass sign, activate outer-boundary facet
// contributions, and restrict the essential set to the reentrant conductor; the metric
// variant is the coercive graph operator used for patch solves and element energies.
struct MaxwellSystemOptions
{
  double mass_coefficient = MaxwellGamma;
  double boundary_conductance = 0.0;
  bool boundary_rhs = false;
  bool pec_everywhere = true;
  bool driven_source = false;
  bool enriched = true;
  int boundary_quadrature_bump = 0;
};

MaxwellSystemOptions MaxwellDrivenResidualOptions()
{
  MaxwellSystemOptions options;
  options.mass_coefficient = -MaxwellOmega2;
  options.boundary_conductance = MaxwellBoundaryConductance;
  options.boundary_rhs = true;
  options.pec_everywhere = false;
  options.driven_source = true;
  return options;
}

MaxwellSystemOptions MaxwellDrivenMetricOptions()
{
  MaxwellSystemOptions options;
  options.mass_coefficient = MaxwellOmega2;
  options.boundary_conductance = MaxwellBoundaryConductance;
  options.boundary_rhs = false;
  options.pec_everywhere = false;
  options.driven_source = false;
  return options;
}

struct MaxwellLocalSystem
{
  int standard_size = 0;
  int enrichment_size = 0;
  std::vector<bool> essential;
  fem::singular::TriangleDofTopology topology;
  std::vector<MaxwellLocalElementData> elements;
  std::vector<MaxwellLocalElementData> facets;
};

std::vector<MaxwellLocalElementData>
MaxwellAllContributions(const MaxwellLocalSystem &system)
{
  auto contributions = system.elements;
  contributions.insert(contributions.end(), system.facets.begin(), system.facets.end());
  return contributions;
}

int MaxwellUnsignedDof(int dof)
{
  return dof >= 0 ? dof : -1 - dof;
}

double MaxwellDofSign(int dof)
{
  return dof >= 0 ? 1.0 : -1.0;
}

MaxwellLocalSystem MaxwellAssembleLocalSystem(mfem::Mesh &mesh, int order,
                                              const MaxwellSystemOptions &options = {})
{
  MaxwellLocalSystem system;
  mfem::ND_FECollection nd_collection(order, 2);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  mfem::H1_FECollection h1_collection(order, 2);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  system.standard_size = nd_space.GetVSize();
  fem::singular::TriangleFeatureTopology features;
  if (options.enriched)
  {
    features = fem::singular::ExtractSerialLineFeatures(mesh, {MaxwellReentrantAttribute},
                                                        {{1, 1.0}});
    system.topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
    system.enrichment_size = static_cast<int>(system.topology.nd_dofs.size());
  }
  const int combined_size = system.standard_size + system.enrichment_size;
  system.essential.assign(combined_size, false);
  mfem::Array<int> standard_essential;
  if (options.pec_everywhere)
  {
    nd_space.GetBoundaryTrueDofs(standard_essential);
  }
  else
  {
    mfem::Array<int> marker(mesh.bdr_attributes.Max());
    marker = 0;
    marker[MaxwellReentrantAttribute - 1] = 1;
    nd_space.GetEssentialTrueDofs(marker, standard_essential);
  }
  for (int dof : standard_essential)
  {
    system.essential[dof] = true;
  }
  if (options.enriched)
  {
    const auto numbering =
        fem::singular::BuildParallelDofNumbering(Mpi::World(), system.topology);
    mfem::Array<int> attributes(options.pec_everywhere ? 2 : 1);
    attributes[0] = MaxwellReentrantAttribute;
    if (options.pec_everywhere)
    {
      attributes[1] = MaxwellOuterAttribute;
    }
    const auto singular_essential = fem::singular::GetEssentialTriangleNDTrueDofs(
        Mpi::World(), features, system.topology, numbering, attributes);
    for (int dof : singular_essential)
    {
      system.essential[system.standard_size + dof] = true;
    }
  }
  static const std::vector<fem::singular::TriangleElementDof> no_singular_dofs;
  const auto element_singular_dofs =
      [&](int element) -> const std::vector<fem::singular::TriangleElementDof> &
  { return options.enriched ? system.topology.elements[element].nd : no_singular_dofs; };

  mfem::CurlCurlIntegrator curl_curl;
  mfem::ConstantCoefficient mass_coefficient(options.mass_coefficient);
  mfem::VectorFEMassIntegrator mass(mass_coefficient);
  const fem::singular::AdaptiveAssemblyOptions assembly_options{8, 1.0e-9, 1.0e-9, 8};
  system.elements.resize(mesh.GetNE());
  mfem::Array<int> nd_dofs, h1_dofs;
  mfem::DenseMatrix shape;
  mfem::Vector point(3);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &data = system.elements[element];
    data.support_element = element;
    const auto &singular_dofs = element_singular_dofs(element);
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
          const auto source = options.driven_source
                                  ? MaxwellDrivenSource(point(0), point(1))
                                  : MaxwellSource(point(0), point(1));
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
          system.topology.elements[element], T, assembly_options);
      auto coupling = fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
          system.topology.elements[element], h1_fe, nd_fe, T, assembly_options);
      fem::singular::ApplyStandardDofTransformations(h1_transformation, nd_transformation,
                                                     coupling);
      for (int row = 0; row < standard_count; row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(row, standard_count + column) =
              MaxwellDofSign(nd_dofs[row]) *
              (coupling.nd_curl_curl_standard_enrichment(row, column) +
               options.mass_coefficient *
                   coupling.nd_mass_standard_enrichment(row, column));
          data.matrix(standard_count + column, row) =
              MaxwellDofSign(nd_dofs[row]) *
              (coupling.nd_curl_curl_enrichment_standard(column, row) +
               options.mass_coefficient *
                   coupling.nd_mass_enrichment_standard(column, row));
        }
      }
      for (int row = 0; row < static_cast<int>(singular_dofs.size()); row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(standard_count + row, standard_count + column) =
              singular.nd_curl_curl(row, column) +
              options.mass_coefficient * singular.nd_mass(row, column);
        }
      }
    }
  }

  if (options.boundary_conductance != 0.0 || options.boundary_rhs)
  {
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT,
                                          2 * order + 4 + options.boundary_quadrature_bump);
    mfem::Vector normal(2), tangent(2);
    for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
    {
      if (mesh.GetBdrAttribute(boundary) != MaxwellOuterAttribute)
      {
        continue;
      }
      auto *face = mesh.GetBdrFaceTransformations(boundary);
      MFEM_VERIFY(face && face->Elem1No >= 0,
                  "Manufactured boundary facet has no interior neighbor!");
      const int element = face->Elem1No;
      mfem::DofTransformation nd_transformation;
      nd_space.GetElementVDofs(element, nd_dofs, nd_transformation);
      const int standard_count = nd_dofs.Size();
      const auto &singular_dofs = element_singular_dofs(element);
      const int singular_count = static_cast<int>(singular_dofs.size());
      auto &data = system.facets.emplace_back();
      data.support_element = element;
      for (int dof : nd_dofs)
      {
        data.dofs.push_back(MaxwellUnsignedDof(dof));
      }
      for (const auto &dof : singular_dofs)
      {
        data.dofs.push_back(system.standard_size + static_cast<int>(dof.dof));
      }
      const int local_size = standard_count + singular_count;
      data.matrix.SetSize(local_size);
      data.matrix = 0.0;
      data.rhs.SetSize(local_size);
      data.rhs = 0.0;
      const auto &nd_fe = *nd_space.GetFE(element);
      auto &element_transformation = *mesh.GetElementTransformation(element);
      double determinant;
      const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
          element_transformation, determinant);
      mfem::Vector centroid(2);
      {
        mfem::Array<int> vertices;
        mesh.GetElementVertices(element, vertices);
        centroid = 0.0;
        for (int vertex : vertices)
        {
          centroid(0) += mesh.GetVertex(vertex)[0] / vertices.Size();
          centroid(1) += mesh.GetVertex(vertex)[1] / vertices.Size();
        }
      }
      mfem::DenseMatrix facet_standard(standard_count);
      facet_standard = 0.0;
      mfem::DenseMatrix facet_coupling(standard_count, singular_count);
      facet_coupling = 0.0;
      mfem::Vector facet_rhs(standard_count);
      facet_rhs = 0.0;
      shape.SetSize(standard_count, 2);
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        face->SetAllIntPoints(&ip);
        const auto &eip = face->GetElement1IntPoint();
        mfem::CalcOrtho(face->Jacobian(), normal);
        face->Transform(ip, point);
        // The counterclockwise boundary tangent is the +90 degree rotation of the outward
        // normal. Verify outwardness directly instead of assuming an orientation.
        const bool outward =
            normal(0) * (point(0) - centroid(0)) + normal(1) * (point(1) - centroid(1)) >
            0.0;
        REQUIRE(outward);
        tangent(0) = -normal(1);
        tangent(1) = normal(0);
        tangent /= tangent.Norml2();
        const double weight = ip.weight * face->Weight();
        auto &T1 = *face->Elem1;
        nd_fe.CalcPhysVShape(T1, shape);
        const double facet_data =
            options.boundary_rhs ? MaxwellFacetData(point(0), point(1)) : 0.0;
        std::vector<double> singular_trace(singular_count, 0.0);
        if (singular_count > 0)
        {
          const fem::singular::TriangleBarycentricPoint lambda{1.0 - eip.x - eip.y, eip.x,
                                                               eip.y};
          for (int k = 0; k < singular_count; k++)
          {
            const auto basis =
                MaxwellSingularBasis(lambda, grad_lambda, singular_dofs[k].basis);
            singular_trace[k] = basis.value[0] * tangent(0) + basis.value[1] * tangent(1);
          }
        }
        for (int i = 0; i < standard_count; i++)
        {
          const double trace_i = shape(i, 0) * tangent(0) + shape(i, 1) * tangent(1);
          facet_rhs(i) += weight * facet_data * trace_i;
          for (int j = 0; j < standard_count; j++)
          {
            const double trace_j = shape(j, 0) * tangent(0) + shape(j, 1) * tangent(1);
            facet_standard(i, j) +=
                options.boundary_conductance * weight * trace_i * trace_j;
          }
          for (int k = 0; k < singular_count; k++)
          {
            facet_coupling(i, k) +=
                options.boundary_conductance * weight * trace_i * singular_trace[k];
          }
        }
        for (int k = 0; k < singular_count; k++)
        {
          data.rhs(standard_count + k) += weight * facet_data * singular_trace[k];
          for (int l = 0; l < singular_count; l++)
          {
            data.matrix(standard_count + k, standard_count + l) +=
                options.boundary_conductance * weight * singular_trace[k] *
                singular_trace[l];
          }
        }
      }
      nd_transformation.TransformDual(facet_standard);
      nd_transformation.TransformDual(facet_rhs);
      for (int k = 0; k < singular_count; k++)
      {
        mfem::Vector column;
        facet_coupling.GetColumnReference(k, column);
        nd_transformation.TransformDual(column);
      }
      for (int i = 0; i < standard_count; i++)
      {
        data.rhs(i) = MaxwellDofSign(nd_dofs[i]) * facet_rhs(i);
        for (int j = 0; j < standard_count; j++)
        {
          data.matrix(i, j) = MaxwellDofSign(nd_dofs[i]) * MaxwellDofSign(nd_dofs[j]) *
                              facet_standard(i, j);
        }
        for (int k = 0; k < singular_count; k++)
        {
          const double value = MaxwellDofSign(nd_dofs[i]) * facet_coupling(i, k);
          data.matrix(i, standard_count + k) = value;
          data.matrix(standard_count + k, i) = value;
        }
      }
    }
  }
  return system;
}

// Solve the (possibly indefinite) manufactured system assembled from the element-local
// uneliminated contributions, including outer-boundary facet terms, with a dense LU on the
// reduced free equations. The scattered operator equals the certified global assembly for
// the coercive options, so this doubles as an assembly-path oracle.
MaxwellSolveReport MaxwellDrivenSolve(mfem::Mesh &mesh, int order,
                                      const MaxwellSystemOptions &options)
{
  auto system = MaxwellAssembleLocalSystem(mesh, order, options);
  const int combined_size = system.standard_size + system.enrichment_size;
  const auto contributions = MaxwellAllContributions(system);
  mfem::Vector rhs(combined_size);
  rhs = 0.0;
  std::vector<int> free, combined_to_free(combined_size, -1);
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (!system.essential[dof])
    {
      combined_to_free[dof] = static_cast<int>(free.size());
      free.push_back(dof);
    }
  }
  mfem::DenseMatrix reduced(static_cast<int>(free.size()));
  reduced = 0.0;
  for (const auto &data : contributions)
  {
    for (int row = 0; row < static_cast<int>(data.dofs.size()); row++)
    {
      rhs(data.dofs[row]) += data.rhs(row);
      const int reduced_row = combined_to_free[data.dofs[row]];
      if (reduced_row < 0)
      {
        continue;
      }
      for (int column = 0; column < static_cast<int>(data.dofs.size()); column++)
      {
        const int reduced_column = combined_to_free[data.dofs[column]];
        if (reduced_column >= 0)
        {
          reduced(reduced_row, reduced_column) += data.matrix(row, column);
        }
      }
    }
  }
  mfem::Vector reduced_rhs(static_cast<int>(free.size()));
  for (int row = 0; row < static_cast<int>(free.size()); row++)
  {
    reduced_rhs(row) = rhs(free[row]);
  }
  mfem::Vector reduced_solution(static_cast<int>(free.size()));
  mfem::DenseMatrixInverse inverse(reduced);
  inverse.Mult(reduced_rhs, reduced_solution);
  mfem::Vector algebraic(static_cast<int>(free.size()));
  reduced.Mult(reduced_solution, algebraic);
  algebraic -= reduced_rhs;
  mfem::Vector solution(combined_size);
  solution = 0.0;
  for (int row = 0; row < static_cast<int>(free.size()); row++)
  {
    solution(free[row]) = reduced_solution(row);
  }

  MaxwellSolveReport report;
  report.dofs = combined_size;
  report.enrichment_dofs = system.enrichment_size;
  report.residual = algebraic.Norml2() / std::max(reduced_rhs.Norml2(), 1.0e-30);
  REQUIRE(report.residual < 1.0e-9);
  for (int dof = 0; dof < system.enrichment_size; dof++)
  {
    const bool is_free = !system.essential[system.standard_size + dof];
    report.free_enrichment_dofs += is_free;
    if (is_free && system.topology.nd_dofs[dof].family ==
                       fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
    {
      report.free_gradient_dofs++;
    }
    if (is_free && system.topology.nd_dofs[dof].family ==
                       fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL)
    {
      report.free_rotational_dofs++;
    }
  }
  report.standard_coefficients.SetSize(system.standard_size);
  for (int i = 0; i < system.standard_size; i++)
  {
    report.standard_coefficients(i) = solution(i);
  }
  report.enrichment_coefficients.SetSize(system.enrichment_size);
  for (int i = 0; i < system.enrichment_size; i++)
  {
    report.enrichment_coefficients(i) = solution(system.standard_size + i);
  }
  report.essential_mask = system.essential;
  mfem::ND_FECollection nd_collection(order, 2);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  MaxwellMeasureGraphError(mesh, nd_space, options.enriched ? &system.topology : nullptr,
                           solution, report);
  return report;
}

std::vector<MaxwellSparseColumn> MaxwellBuildSparseTransfer(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, double &maximum_consistency_error,
    int *signed_coarse_dofs = nullptr, int *signed_fine_dofs = nullptr,
    int *nonidentity_transformations = nullptr)
{
  auto injection = fem::hierarchical::BuildSparsePInjection(mesh, coarse_space, fine_space);
  maximum_consistency_error = injection.consistency_error;
  if (signed_coarse_dofs)
  {
    *signed_coarse_dofs += injection.signed_coarse_dofs;
  }
  if (signed_fine_dofs)
  {
    *signed_fine_dofs += injection.signed_fine_dofs;
  }
  if (nonidentity_transformations)
  {
    *nonidentity_transformations += injection.nonidentity_transformations;
  }
  return std::move(injection.columns);
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

MaxwellPatchReport
MaxwellBuildLocalPatches(mfem::Mesh &mesh, int order, const MaxwellSolveReport &coarse,
                         bool include_rotational_guests = true,
                         const MaxwellSystemOptions &residual_options = {},
                         const MaxwellSystemOptions *metric_options = nullptr)
{
  MaxwellPatchReport report;
  auto fine = MaxwellAssembleLocalSystem(mesh, order + 1, residual_options);
  const auto residual_all = MaxwellAllContributions(fine);
  // Patch solves, defect sweeps, and element energies use a coercive graph metric. For the
  // certified coercive fixture the metric is the residual operator itself; the driven
  // variant lifts the indefinite residual in the positive graph metric instead.
  std::vector<MaxwellLocalElementData> metric_all;
  if (metric_options)
  {
    auto metric_system = MaxwellAssembleLocalSystem(mesh, order + 1, *metric_options);
    REQUIRE(metric_system.essential == fine.essential);
    metric_all = MaxwellAllContributions(metric_system);
  }
  else
  {
    metric_all = residual_all;
  }
  mfem::ND_FECollection coarse_collection(order, 2), fine_collection(order + 1, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  const int enrichment = fine.enrichment_size;
  REQUIRE(coarse.enrichment_dofs == enrichment);
  const auto injection =
      fem::hierarchical::BuildSparsePInjection(mesh, coarse_space, fine_space);
  report.transfer_consistency_error = injection.consistency_error;
  report.signed_coarse_dofs = injection.signed_coarse_dofs;
  report.signed_fine_dofs = injection.signed_fine_dofs;
  report.nonidentity_transformations = injection.nonidentity_transformations;

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
    for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
    {
      sparse_image(injection.columns[column].dofs[entry]) +=
          probe(column) * injection.columns[column].values[entry];
    }
  }
  mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
  transfer.Mult(probe, reference);
  sparse_image -= reference;
  report.transfer_reference_error = sparse_image.Normlinf();

  // Combined coarse solution and per-element enrichment guest lists, including the
  // optional rotational-family ablation filter.
  mfem::Vector coarse_combined(coarse_standard + enrichment);
  for (int i = 0; i < coarse_standard; i++)
  {
    coarse_combined(i) = coarse.standard_coefficients(i);
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    coarse_combined(coarse_standard + dof) = coarse.enrichment_coefficients(dof);
  }
  std::vector<std::vector<int>> element_enrichment_guests(mesh.GetNE());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    for (const auto &dof : fine.topology.elements[element].nd)
    {
      const bool rotational = fine.topology.nd_dofs[dof.dof].family ==
                              fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
      if (include_rotational_guests || !rotational)
      {
        element_enrichment_guests[element].push_back(static_cast<int>(dof.dof));
      }
    }
  }

  const auto lifting = fem::hierarchical::EstimateByPatchLifting(
      mesh, coarse_space, fine_space, injection, residual_all, metric_all, fine.essential,
      coarse.essential_mask, coarse_combined, element_enrichment_guests);
  REQUIRE(lifting.face_patches == 0);
  report.indicator = lifting.indicator;
  report.energy = lifting.energy;
  report.work = lifting.work;
  report.maximum_patch_residual = lifting.maximum_patch_residual;
  report.maximum_patch_condition = lifting.maximum_patch_condition;
  report.edge_patches = lifting.edge_patches;
  report.interior_patches = lifting.interior_patches;
  report.owned_modes = lifting.owned_modes;
  report.maximum_support_elements = lifting.maximum_support_elements;
  report.maximum_patch_dimension = lifting.maximum_patch_dimension;
  report.maximum_element_overlap = lifting.maximum_element_overlap;
  for (int dof = 0; dof < coarse_standard; dof++)
  {
    if (!coarse.essential_mask[dof])
    {
      REQUIRE(lifting.coarse_guest_counts[dof] > 0);
      report.unique_coarse_guests++;
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    const bool gradient = fine.topology.nd_dofs[dof].family ==
                          fem::singular::HigherOrderBasisFamily::NODE_GRADIENT;
    const bool rotational = fine.topology.nd_dofs[dof].family ==
                            fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
    report.covered_gradient_guests += gradient ? lifting.enrichment_guest_counts[dof] : 0;
    report.covered_rotational_guests +=
        rotational ? lifting.enrichment_guest_counts[dof] : 0;
    if (!fine.essential[fine_standard + dof] && (include_rotational_guests || !rotational))
    {
      REQUIRE(lifting.enrichment_guest_counts[dof] > 0);
      report.unique_gradient_guests += gradient;
      report.unique_rotational_guests += rotational;
    }
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

TEST_CASE("Manufactured driven Maxwell wedge reduces to the coercive oracle",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellLShapeMesh(2);
  const auto reference = MaxwellSolve(mesh, order, true);
  // Default options reproduce the certified coercive fixture: same operator, source, and
  // essential set, assembled from element-local contributions and solved by dense LU.
  const auto scattered = MaxwellDrivenSolve(mesh, order, MaxwellSystemOptions{});
  REQUIRE(scattered.dofs == reference.dofs);
  double coefficient_difference = 0.0;
  for (int i = 0; i < reference.standard_coefficients.Size(); i++)
  {
    coefficient_difference =
        std::max(coefficient_difference, std::abs(scattered.standard_coefficients(i) -
                                                  reference.standard_coefficients(i)));
  }
  for (int i = 0; i < reference.enrichment_coefficients.Size(); i++)
  {
    coefficient_difference =
        std::max(coefficient_difference, std::abs(scattered.enrichment_coefficients(i) -
                                                  reference.enrichment_coefficients(i)));
  }
  CAPTURE(order, coefficient_difference, scattered.residual, reference.residual,
          scattered.graph_error, reference.graph_error);
  CHECK(coefficient_difference < 1.0e-6);
  CHECK_THAT(scattered.graph_error, WithinRel(reference.graph_error, 1.0e-8));
}

TEST_CASE("Manufactured driven Maxwell wedge boundary facets are consistent",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellLShapeMesh(2);
  const auto options = MaxwellDrivenResidualOptions();
  const auto system = MaxwellAssembleLocalSystem(mesh, order, options);
  int outer_boundary_elements = 0;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    outer_boundary_elements += mesh.GetBdrAttribute(boundary) == MaxwellOuterAttribute;
  }
  REQUIRE(outer_boundary_elements > 0);
  REQUIRE(static_cast<int>(system.facets.size()) == outer_boundary_elements);

  // Facet matrices are symmetric and positive semidefinite; the manufactured facet load is
  // active on the y = +-1 segments.
  double asymmetry = 0.0, maximum_rhs = 0.0, minimum_quadratic_form = 0.0;
  for (const auto &facet : system.facets)
  {
    for (int i = 0; i < facet.matrix.Height(); i++)
    {
      maximum_rhs = std::max(maximum_rhs, std::abs(facet.rhs(i)));
      for (int j = 0; j < facet.matrix.Width(); j++)
      {
        asymmetry = std::max(asymmetry, std::abs(facet.matrix(i, j) - facet.matrix(j, i)));
      }
    }
    mfem::Vector probe(facet.matrix.Height()), action(facet.matrix.Height());
    for (int trial = 0; trial < 3; trial++)
    {
      for (int i = 0; i < probe.Size(); i++)
      {
        probe(i) = std::sin(0.61 * i + 0.17 * trial + 0.4);
      }
      facet.matrix.Mult(probe, action);
      minimum_quadratic_form =
          std::min(minimum_quadratic_form, mfem::InnerProduct(probe, action));
    }
  }
  CAPTURE(order, asymmetry, maximum_rhs, minimum_quadratic_form);
  CHECK(asymmetry < 1.0e-14);
  CHECK(maximum_rhs > 1.0e-3);
  CHECK(minimum_quadratic_form > -1.0e-13);

  // Independent quadrature-refinement oracle for every facet entry.
  auto refined_options = options;
  refined_options.boundary_quadrature_bump = 8;
  const auto refined = MaxwellAssembleLocalSystem(mesh, order, refined_options);
  REQUIRE(refined.facets.size() == system.facets.size());
  double matrix_difference = 0.0, rhs_difference = 0.0;
  for (std::size_t facet = 0; facet < system.facets.size(); facet++)
  {
    const auto &base = system.facets[facet];
    const auto &bumped = refined.facets[facet];
    REQUIRE(base.dofs == bumped.dofs);
    for (int i = 0; i < base.matrix.Height(); i++)
    {
      rhs_difference = std::max(rhs_difference, std::abs(base.rhs(i) - bumped.rhs(i)));
      for (int j = 0; j < base.matrix.Width(); j++)
      {
        matrix_difference =
            std::max(matrix_difference, std::abs(base.matrix(i, j) - bumped.matrix(i, j)));
      }
    }
  }
  CAPTURE(matrix_difference, rhs_difference);
  CHECK(matrix_difference < 1.0e-12);
  CHECK(rhs_difference < 1.0e-12);

  // The driven essential set frees the outer boundary while keeping the reentrant PEC
  // conductor constrained.
  const auto coercive = MaxwellAssembleLocalSystem(mesh, order, MaxwellSystemOptions{});
  const int driven_essential =
      static_cast<int>(std::count(system.essential.begin(), system.essential.end(), true));
  const int coercive_essential = static_cast<int>(
      std::count(coercive.essential.begin(), coercive.essential.end(), true));
  CAPTURE(driven_essential, coercive_essential);
  CHECK(driven_essential > 0);
  CHECK(driven_essential < coercive_essential);
}

TEST_CASE("Manufactured driven Maxwell wedge converges to the manufactured field",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  const auto options = MaxwellDrivenResidualOptions();
  auto coarse_mesh = MaxwellLShapeMesh(2);
  auto fine_mesh = MaxwellLShapeMesh(4);
  const auto coarse = MaxwellDrivenSolve(coarse_mesh, order, options);
  const auto fine = MaxwellDrivenSolve(fine_mesh, order, options);
  auto standard_options = options;
  standard_options.enriched = false;
  const auto standard = MaxwellDrivenSolve(coarse_mesh, order, standard_options);
  CAPTURE(order, coarse.dofs, fine.dofs, coarse.graph_error, fine.graph_error,
          standard.graph_error, coarse.free_gradient_dofs, coarse.free_rotational_dofs);
  CHECK(coarse.free_gradient_dofs > 0);
  CHECK(coarse.free_rotational_dofs > 0);
  CHECK(coarse.graph_error < 0.9 * standard.graph_error);
  CHECK(fine.graph_error < 0.5 * coarse.graph_error);
}

TEST_CASE("Manufactured driven Maxwell wedge local patches rank the true error",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellLShapeMesh(2);
  const auto residual_options = MaxwellDrivenResidualOptions();
  const auto metric_options = MaxwellDrivenMetricOptions();
  const auto coarse = MaxwellDrivenSolve(mesh, order, residual_options);
  const auto local = MaxwellBuildLocalPatches(mesh, order, coarse, true, residual_options,
                                              &metric_options);
  CAPTURE(order, local.energy, local.work, local.transfer_consistency_error,
          local.transfer_reference_error, local.maximum_patch_residual,
          local.maximum_patch_condition, local.edge_patches, local.interior_patches,
          local.owned_modes, local.maximum_support_elements, local.maximum_patch_dimension,
          local.maximum_element_overlap, local.unique_coarse_guests,
          local.unique_gradient_guests, local.unique_rotational_guests);
  CHECK(local.transfer_consistency_error < 1.0e-12);
  CHECK(local.transfer_reference_error < 1.0e-12);
  CHECK(local.maximum_patch_residual < 1.0e-10);
  CHECK(local.edge_patches > 0);
  CHECK(local.interior_patches > 0);
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
  CAPTURE(rank, extend, marked.size(), oracle.size());
  if (order == 1)
  {
    // The p1 first iterate is the known weaker case (see the coercive fixture, which also
    // gates p1 below p2). Measured driven values: rank 0.798, extend 0.837.
    CHECK(rank > 0.75);
    CHECK(extend > 0.8);
  }
  else
  {
    // Measured driven values: rank 0.904, extend 0.874.
    CHECK(rank > 0.85);
    CHECK(extend > 0.8);
  }
}

TEST_CASE("Manufactured driven Maxwell wedge sparse ND AMR beats uniform",
          "[singularmaxwellmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 3;
  const auto residual_options = MaxwellDrivenResidualOptions();
  const auto metric_options = MaxwellDrivenMetricOptions();
  auto mesh = MaxwellLShapeMesh(2);
  std::vector<long long> dofs;
  std::vector<double> errors;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto solve = MaxwellDrivenSolve(mesh, order, residual_options);
    const auto local = MaxwellBuildLocalPatches(mesh, order, solve, true, residual_options,
                                                &metric_options);
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
    const auto solve = MaxwellDrivenSolve(uniform_mesh, order, residual_options);
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

namespace
{

// Extruded three-dimensional L-shaped wedge: the reentrant vertical sheet faces meet at
// the singular edge x = y = 0 with the same 3 pi / 2 opening as the two-dimensional
// fixture. Every prism is split by the same Freudenthal diagonal so shared faces agree.
mfem::Mesh MaxwellExtrudedLShapeMesh(int n, int nz)
{
  REQUIRE(n > 0);
  REQUIRE(nz > 0);
  const int grid_vertices = (2 * n + 1) * (2 * n + 1) * (nz + 1);
  const int active_vertices = grid_vertices - n * n * (nz + 1);
  const int elements = 18 * n * n * nz;
  mfem::Mesh mesh(3, active_vertices, elements, 0, 3);
  std::vector<int> vertex_ids(grid_vertices, -1);
  const auto flat = [n](int i, int j, int k)
  { return (k * (2 * n + 1) + j) * (2 * n + 1) + i; };
  const auto vertex = [&](int i, int j, int k)
  {
    const int index = flat(i, j, k);
    if (vertex_ids[index] >= 0)
    {
      return vertex_ids[index];
    }
    const double point[3]{-1.0 + static_cast<double>(i) / n,
                          -1.0 + static_cast<double>(j) / n,
                          -1.0 + 2.0 * static_cast<double>(k) / nz};
    vertex_ids[index] = mesh.AddVertex(point);
    return vertex_ids[index];
  };
  const auto add_oriented_tet = [&](std::array<int, 4> tet)
  {
    const double *a = mesh.GetVertex(tet[0]);
    const double *b = mesh.GetVertex(tet[1]);
    const double *c = mesh.GetVertex(tet[2]);
    const double *d = mesh.GetVertex(tet[3]);
    const double determinant =
        (b[0] - a[0]) * ((c[1] - a[1]) * (d[2] - a[2]) - (c[2] - a[2]) * (d[1] - a[1])) -
        (b[1] - a[1]) * ((c[0] - a[0]) * (d[2] - a[2]) - (c[2] - a[2]) * (d[0] - a[0])) +
        (b[2] - a[2]) * ((c[0] - a[0]) * (d[1] - a[1]) - (c[1] - a[1]) * (d[0] - a[0]));
    REQUIRE(std::abs(determinant) > 1.0e-14);
    if (determinant < 0.0)
    {
      std::swap(tet[1], tet[2]);
    }
    mesh.AddTet(tet.data(), 1);
  };
  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < 2 * n; j++)
    {
      for (int i = 0; i < 2 * n; i++)
      {
        if (i >= n && j < n)
        {
          continue;
        }
        const int p000 = vertex(i, j, k), p100 = vertex(i + 1, j, k);
        const int p010 = vertex(i, j + 1, k), p110 = vertex(i + 1, j + 1, k);
        const int p001 = vertex(i, j, k + 1), p101 = vertex(i + 1, j, k + 1);
        const int p011 = vertex(i, j + 1, k + 1), p111 = vertex(i + 1, j + 1, k + 1);
        for (const auto &tet :
             std::array<std::array<int, 4>, 6>{std::array<int, 4>{p000, p100, p110, p111},
                                               std::array<int, 4>{p000, p100, p101, p111},
                                               std::array<int, 4>{p000, p010, p110, p111},
                                               std::array<int, 4>{p000, p010, p011, p111},
                                               std::array<int, 4>{p000, p001, p101, p111},
                                               std::array<int, 4>{p000, p001, p011, p111}})
        {
          add_oriented_tet(tet);
        }
      }
    }
  }
  mesh.FinalizeTopology(true);
  mfem::Array<int> face_vertices;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    mesh.GetBdrElementVertices(boundary, face_vertices);
    bool on_x_reentrant = true, on_y_reentrant = true;
    for (int v : face_vertices)
    {
      const double *point = mesh.GetVertex(v);
      on_x_reentrant = on_x_reentrant && std::abs(point[0]) < 1.0e-13 &&
                       point[1] <= 1.0e-13 && point[1] >= -1.0 - 1.0e-13;
      on_y_reentrant = on_y_reentrant && std::abs(point[1]) < 1.0e-13 &&
                       point[0] >= -1.0e-13 && point[0] <= 1.0 + 1.0e-13;
    }
    mesh.GetBdrElement(boundary)->SetAttribute((on_x_reentrant || on_y_reentrant)
                                                   ? MaxwellReentrantAttribute
                                                   : MaxwellOuterAttribute);
  }
  mesh.Finalize(true, false);
  return mesh;
}

struct MaxwellTetSystem
{
  int standard_size = 0;
  int enrichment_size = 0;
  std::vector<bool> essential;
  fem::singular::FeatureTopology features;
  fem::singular::DofTopology topology;
  std::vector<MaxwellLocalElementData> elements;
};

// Assemble the coercive combined graph system curl-curl plus mass on tetrahedra from the
// production retained element patch matrices for enriched elements and standard MFEM
// element assembly elsewhere. This is exactly the data path a driver-level estimator
// consumes from SpaceOperator.
MaxwellTetSystem MaxwellAssembleTetLocalSystem(mfem::Mesh &mesh, int order)
{
  MaxwellTetSystem system;
  mfem::ND_FECollection nd_collection(order, 3);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  mfem::H1_FECollection h1_collection(order, 3);
  mfem::FiniteElementSpace h1_space(&mesh, &h1_collection);
  system.standard_size = nd_space.GetVSize();
  system.features = fem::singular::ExtractSerialSheetFeatures(
      mesh, {MaxwellReentrantAttribute},
      std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  REQUIRE(!system.features.features.empty());
  system.topology = fem::singular::BuildSerialDofTopology(mesh, system.features, 1);
  system.enrichment_size = static_cast<int>(system.topology.nd_dofs.size());
  REQUIRE(system.enrichment_size > 0);
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
  const auto singular_essential = fem::singular::GetEssentialNDTrueDofs(
      Mpi::World(), system.features, system.topology, numbering, attributes);
  for (int dof : singular_essential)
  {
    system.essential[system.standard_size + dof] = true;
  }

  const std::vector<fem::singular::IsotropicMaterialCoefficients> materials(mesh.GetNE(),
                                                                            {1.0, 1.0});
  const fem::singular::AdaptiveAssemblyOptions options{8, 5.0e-7, 1.0e-6, 8};
  const auto batches = fem::singular::AssembleLocalSparseEnrichmentMatricesBatch(
      system.topology, h1_space, nd_space, {materials}, options, 0);
  const auto &retained = batches.front().nd_element_patches;
  REQUIRE(!retained.empty());
  std::vector<int> retained_index(mesh.GetNE(), -1);
  for (std::size_t patch = 0; patch < retained.size(); patch++)
  {
    retained_index[retained[patch].element] = static_cast<int>(patch);
  }

  mfem::CurlCurlIntegrator curl_curl;
  mfem::VectorFEMassIntegrator mass;
  system.elements.resize(mesh.GetNE());
  mfem::Array<int> nd_dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &data = system.elements[element];
    data.support_element = element;
    if (retained_index[element] >= 0)
    {
      const auto &patch = retained[retained_index[element]];
      const int standard_count = patch.standard_dofs.Size();
      const int local_size = standard_count + patch.enrichment_dofs.Size();
      for (int dof : patch.standard_dofs)
      {
        data.dofs.push_back(MaxwellUnsignedDof(dof));
      }
      for (int dof : patch.enrichment_dofs)
      {
        data.dofs.push_back(system.standard_size + dof);
      }
      data.matrix.SetSize(local_size);
      data.rhs.SetSize(local_size);
      data.rhs = 0.0;
      for (int row = 0; row < local_size; row++)
      {
        const double row_sign =
            row < standard_count ? MaxwellDofSign(patch.standard_dofs[row]) : 1.0;
        for (int column = 0; column < local_size; column++)
        {
          const double column_sign =
              column < standard_count ? MaxwellDofSign(patch.standard_dofs[column]) : 1.0;
          data.matrix(row, column) =
              row_sign * column_sign *
              (patch.curl_curl(row, column) + patch.mass(row, column));
        }
      }
    }
    else
    {
      mfem::DofTransformation transformation;
      nd_space.GetElementVDofs(element, nd_dofs, transformation);
      const int standard_count = nd_dofs.Size();
      for (int dof : nd_dofs)
      {
        data.dofs.push_back(MaxwellUnsignedDof(dof));
      }
      mfem::DenseMatrix local_curl, local_mass;
      auto &T = *mesh.GetElementTransformation(element);
      const auto &fe = *nd_space.GetFE(element);
      curl_curl.AssembleElementMatrix(fe, T, local_curl);
      mass.AssembleElementMatrix(fe, T, local_mass);
      local_curl += local_mass;
      transformation.TransformDual(local_curl);
      data.matrix.SetSize(standard_count);
      data.rhs.SetSize(standard_count);
      data.rhs = 0.0;
      for (int row = 0; row < standard_count; row++)
      {
        for (int column = 0; column < standard_count; column++)
        {
          data.matrix(row, column) = MaxwellDofSign(nd_dofs[row]) *
                                     MaxwellDofSign(nd_dofs[column]) *
                                     local_curl(row, column);
        }
      }
    }
  }
  return system;
}

// Conjugate-gradient solve of the reduced coercive combined system scattered from the
// element contributions.
mfem::Vector MaxwellTetReducedSolve(const MaxwellTetSystem &system, const mfem::Vector &rhs)
{
  const int combined_size = system.standard_size + system.enrichment_size;
  std::vector<int> free, combined_to_free(combined_size, -1);
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (!system.essential[dof])
    {
      combined_to_free[dof] = static_cast<int>(free.size());
      free.push_back(dof);
    }
  }
  mfem::SparseMatrix reduced(static_cast<int>(free.size()), static_cast<int>(free.size()));
  for (const auto &data : system.elements)
  {
    for (int row = 0; row < static_cast<int>(data.dofs.size()); row++)
    {
      const int reduced_row = combined_to_free[data.dofs[row]];
      if (reduced_row < 0)
      {
        continue;
      }
      for (int column = 0; column < static_cast<int>(data.dofs.size()); column++)
      {
        const int reduced_column = combined_to_free[data.dofs[column]];
        if (reduced_column >= 0 && data.matrix(row, column) != 0.0)
        {
          reduced.Add(reduced_row, reduced_column, data.matrix(row, column));
        }
      }
    }
  }
  reduced.Finalize();
  mfem::Vector reduced_rhs(static_cast<int>(free.size()));
  for (int row = 0; row < static_cast<int>(free.size()); row++)
  {
    reduced_rhs(row) = rhs(free[row]);
  }
  mfem::Vector reduced_solution(static_cast<int>(free.size()));
  reduced_solution = 0.0;
  mfem::GSSmoother preconditioner(reduced);
  mfem::CGSolver cg;
  cg.SetOperator(reduced);
  cg.SetPreconditioner(preconditioner);
  cg.SetRelTol(1.0e-12);
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
  return solution;
}

}  // namespace

TEST_CASE("Extruded Maxwell wedge sparse ND patches certify the 3D engine",
          "[.][singularmaxwellmms3d][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  auto mesh = MaxwellExtrudedLShapeMesh(1, 2);
  auto coarse = MaxwellAssembleTetLocalSystem(mesh, order);
  auto fine = MaxwellAssembleTetLocalSystem(mesh, order + 1);
  REQUIRE(coarse.enrichment_size == fine.enrichment_size);
  const int enrichment = fine.enrichment_size;
  const int combined_size = fine.standard_size + enrichment;
  for (int dof = 0; dof < enrichment; dof++)
  {
    REQUIRE(coarse.essential[coarse.standard_size + dof] ==
            fine.essential[fine.standard_size + dof]);
  }

  mfem::ND_FECollection coarse_collection(order, 3), fine_collection(order + 1, 3);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  REQUIRE(coarse_standard == coarse.standard_size);
  REQUIRE(fine_standard == fine.standard_size);

  // Tetrahedral signed injection against MFEM's p-transfer; three-dimensional ND spaces
  // must exercise genuine nonidentity DOF transformations.
  const auto injection =
      fem::hierarchical::BuildSparsePInjection(mesh, coarse_space, fine_space);
  mfem::Vector probe(coarse_standard), sparse_image(fine_standard),
      reference(fine_standard);
  for (int i = 0; i < coarse_standard; i++)
  {
    probe(i) = std::sin(0.37 * i + 0.2);
  }
  sparse_image = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
    {
      sparse_image(injection.columns[column].dofs[entry]) +=
          probe(column) * injection.columns[column].values[entry];
    }
  }
  mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
  transfer.Mult(probe, reference);
  sparse_image -= reference;
  CAPTURE(order, injection.consistency_error, injection.signed_coarse_dofs,
          injection.signed_fine_dofs, injection.nonidentity_transformations,
          sparse_image.Normlinf());
  CHECK(injection.consistency_error < 1.0e-12);
  CHECK(sparse_image.Normlinf() < 1.0e-12);
  CHECK(injection.signed_coarse_dofs > 0);
  CHECK(injection.signed_fine_dofs > 0);
  CHECK(injection.nonidentity_transformations > 0);

  // Galerkin consistency of the injected fine operator against the coarse assembly.
  {
    mfem::Vector coarse_probe(coarse_standard + enrichment), injected_probe(combined_size);
    for (int i = 0; i < coarse_probe.Size(); i++)
    {
      coarse_probe(i) = std::cos(0.23 * i + 0.11);
    }
    injected_probe = 0.0;
    for (int column = 0; column < coarse_standard; column++)
    {
      for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
      {
        injected_probe(injection.columns[column].dofs[entry]) +=
            coarse_probe(column) * injection.columns[column].values[entry];
      }
    }
    for (int dof = 0; dof < enrichment; dof++)
    {
      injected_probe(fine_standard + dof) = coarse_probe(coarse_standard + dof);
    }
    const double coarse_form = fem::hierarchical::Energy(coarse.elements, coarse_probe);
    const double fine_form = fem::hierarchical::Energy(fine.elements, injected_probe);
    CAPTURE(coarse_form, fine_form);
    CHECK_THAT(fine_form, WithinRel(coarse_form, 1.0e-5));
  }

  // Synthetic fine truth: a smooth projected standard field plus explicit free gradient
  // and rotational singular coefficients. The manufactured load is its exact fine action.
  mfem::Vector truth(combined_size);
  truth = 0.0;
  {
    mfem::GridFunction smooth(&fine_space);
    mfem::VectorFunctionCoefficient smooth_field(
        3,
        [](const mfem::Vector &x, mfem::Vector &value)
        {
          value(0) = x(1) * x(2) + 0.3 * x(0);
          value(1) = x(0) - 0.5 * x(2) * x(2);
          value(2) = x(0) * x(1) + 0.25 * x(2);
        });
    smooth.ProjectCoefficient(smooth_field);
    for (int i = 0; i < fine_standard; i++)
    {
      truth(i) = smooth(i);
    }
  }
  int free_gradient = -1, free_rotational = -1;
  for (int dof = 0; dof < enrichment; dof++)
  {
    if (fine.essential[fine_standard + dof])
    {
      continue;
    }
    const bool gradient =
        fem::singular::IsGradientFamily(fine.topology.nd_dofs[dof].family);
    if (gradient && free_gradient < 0)
    {
      free_gradient = dof;
    }
    if (!gradient && free_rotational < 0)
    {
      free_rotational = dof;
    }
  }
  CAPTURE(order, enrichment, free_gradient, free_rotational);
  REQUIRE(free_gradient >= 0);
  REQUIRE(free_rotational >= 0);
  truth(fine_standard + free_gradient) = 0.3;
  truth(fine_standard + free_rotational) = -0.2;
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (fine.essential[dof])
    {
      truth(dof) = 0.0;
    }
  }
  mfem::Vector fine_rhs(combined_size);
  fine_rhs = 0.0;
  for (const auto &data : fine.elements)
  {
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = truth(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    for (int i = 0; i < local.Size(); i++)
    {
      fine_rhs(data.dofs[i]) += action(i);
    }
  }
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (fine.essential[dof])
    {
      fine_rhs(dof) = 0.0;
    }
  }

  // Coarse Galerkin problem for the same load, solved on the coarse combined layout.
  mfem::Vector coarse_rhs(coarse_standard + enrichment);
  coarse_rhs = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
    {
      coarse_rhs(column) += injection.columns[column].values[entry] *
                            fine_rhs(injection.columns[column].dofs[entry]);
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    coarse_rhs(coarse_standard + dof) = fine_rhs(fine_standard + dof);
  }
  const auto coarse_solution = MaxwellTetReducedSolve(coarse, coarse_rhs);

  // Attach the manufactured load to the fine contributions so the engine sees b - A x.
  auto fine_loaded = fine.elements;
  {
    std::vector<bool> assigned(combined_size, false);
    for (auto &data : fine_loaded)
    {
      for (int i = 0; i < static_cast<int>(data.dofs.size()); i++)
      {
        if (!assigned[data.dofs[i]])
        {
          data.rhs(i) = fine_rhs(data.dofs[i]);
          assigned[data.dofs[i]] = true;
        }
      }
    }
  }

  // Global hierarchical benchmark: the exact fine correction energy per element.
  mfem::Vector injected_coarse(combined_size);
  injected_coarse = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
    {
      injected_coarse(injection.columns[column].dofs[entry]) +=
          coarse_solution(column) * injection.columns[column].values[entry];
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    injected_coarse(fine_standard + dof) = coarse_solution(coarse_standard + dof);
  }
  const auto residual = fem::hierarchical::AssembleResidual(
      combined_size, fine_loaded, injected_coarse, fine.essential);
  const auto global_correction = MaxwellTetReducedSolve(fine, residual);
  std::vector<double> global_indicator(mesh.GetNE(), 0.0);
  double global_energy = 0.0;
  for (const auto &data : fine.elements)
  {
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = global_correction(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    global_indicator[data.support_element] += mfem::InnerProduct(local, action);
    global_energy += mfem::InnerProduct(local, action);
  }
  REQUIRE(global_energy > 0.0);

  // Element-local sparse lifting through the shared engine.
  std::vector<std::vector<int>> element_enrichment_guests(mesh.GetNE());
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    for (const auto &dof : fine.topology.elements[element].nd)
    {
      element_enrichment_guests[element].push_back(static_cast<int>(dof.dof));
    }
  }
  const auto lifting = fem::hierarchical::EstimateByPatchLifting(
      mesh, coarse_space, fine_space, injection, fine_loaded, fine_loaded, fine.essential,
      coarse.essential, coarse_solution, element_enrichment_guests);
  CAPTURE(lifting.energy, lifting.work, lifting.edge_patches, lifting.face_patches,
          lifting.interior_patches, lifting.owned_modes, lifting.maximum_support_elements,
          lifting.maximum_patch_dimension, lifting.maximum_element_overlap,
          lifting.maximum_patch_residual, lifting.maximum_patch_condition, global_energy);
  CHECK(lifting.energy > 0.0);
  CHECK_THAT(lifting.energy, WithinRel(lifting.work, 1.0e-10));
  CHECK_THAT(std::accumulate(lifting.indicator.begin(), lifting.indicator.end(), 0.0),
             WithinRel(lifting.energy, 1.0e-12));
  CHECK(lifting.edge_patches > 0);
  CHECK(lifting.face_patches > 0);
  if (order >= 2)
  {
    CHECK(lifting.interior_patches > 0);
  }
  CHECK(lifting.maximum_patch_residual < 1.0e-10);
  for (int dof = 0; dof < coarse_standard; dof++)
  {
    if (!coarse.essential[dof])
    {
      REQUIRE(lifting.coarse_guest_counts[dof] > 0);
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    if (!fine.essential[fine_standard + dof])
    {
      REQUIRE(lifting.enrichment_guest_counts[dof] > 0);
    }
  }

  // Marking effectivity against the global hierarchical correction.
  const auto marked = MaxwellDorfler(lifting.indicator, 0.5);
  const auto oracle = MaxwellDorfler(global_indicator, 0.5);
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
      MaxwellSumOn(global_indicator, marked) /
      MaxwellSumOn(global_indicator, top_k(global_indicator, marked.size()));
  const double extend =
      MaxwellSumOn(global_indicator, top_k(lifting.indicator, oracle.size())) /
      MaxwellSumOn(global_indicator, oracle);
  const double capture = lifting.energy / global_energy;
  CAPTURE(rank, extend, capture, marked.size(), oracle.size());
  // Measured: p1 rank 1.0, extend 1.0, capture 0.713; p2 rank 0.9995, extend 0.942,
  // capture 0.929.
  CHECK(rank > 0.95);
  CHECK(extend > 0.9);
  CHECK(capture > 0.6);
  CHECK(capture < 1.5);
}

}  // namespace palace
