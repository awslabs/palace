// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Manufactured-solution certification of the singular H1 space and the error estimator on
// the L-shaped domain with a reentrant 2/3 corner. See
// ~/bedrock-tests/mms-lshape-20260803/SPEC.md for the derivation and the independently
// cross-checked reference energy.
//
//   Omega = (-1,1)^2 \ ([0,1] x [-1,0]),  nu = 2/3,  theta in [0, 3pi/2]
//   s = r^nu sin(nu theta)   (harmonic)
//   B = (1-x^2)(1-y^2)
//   u = B s                  (vanishes on all six faces)
//   f = -s Lap(B) - 2 grad(B).grad(s),  so -Lap(u) = f
//
// Reference: a(u,u) = int f u = 1.710627311944, established by two method-independent
// high-order quadratures agreeing to 2.3e-13.
//
// The certification is deliberately staged, and each stage gates the next: assemble and
// solve, then algebraic checks, then agreement of three independent error expressions,
// then convergence, and only then any statement about the estimator.

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
#include "fem/singularfield.hpp"
#include "utils/communication.hpp"

namespace palace
{

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

constexpr double MmsNu = 2.0 / 3.0;
// Reference energy a(u,u) = int_Omega f u, cross-checked by graded-Cartesian and polar
// quadratures to 2.3e-13. See SPEC.md; this is an empirical cross-method agreement, not a
// rigorous bound, but it sits ~10 orders below any discretization error measured here.
constexpr double MmsReferenceEnergy = 1.710627311944;

double MmsTheta(double x, double y)
{
  const double t = std::atan2(y, x);
  return (t < 0.0) ? t + 2.0 * M_PI : t;
}

double MmsS(double x, double y)
{
  const double r = std::hypot(x, y);
  return (r == 0.0) ? 0.0 : std::pow(r, MmsNu) * std::sin(MmsNu * MmsTheta(x, y));
}

std::array<double, 2> MmsGradS(double x, double y)
{
  const double r = std::hypot(x, y);
  if (r == 0.0)
  {
    return {0.0, 0.0};
  }
  const double t = MmsTheta(x, y);
  const double c = MmsNu * std::pow(r, MmsNu - 1.0);
  return {c * std::sin((MmsNu - 1.0) * t), c * std::cos((MmsNu - 1.0) * t)};
}

double MmsB(double x, double y)
{
  return (1.0 - x * x) * (1.0 - y * y);
}

std::array<double, 2> MmsGradB(double x, double y)
{
  return {-2.0 * x * (1.0 - y * y), -2.0 * y * (1.0 - x * x)};
}

double MmsLapB(double x, double y)
{
  return -2.0 * (1.0 - y * y) - 2.0 * (1.0 - x * x);
}

double MmsU(double x, double y)
{
  return MmsB(x, y) * MmsS(x, y);
}

std::array<double, 2> MmsGradU(double x, double y)
{
  const double B = MmsB(x, y), s = MmsS(x, y);
  const auto gB = MmsGradB(x, y), gs = MmsGradS(x, y);
  return {B * gs[0] + s * gB[0], B * gs[1] + s * gB[1]};
}

double MmsF(double x, double y)
{
  const auto gB = MmsGradB(x, y), gs = MmsGradS(x, y);
  return -MmsS(x, y) * MmsLapB(x, y) - 2.0 * (gB[0] * gs[0] + gB[1] * gs[1]);
}

// Boundary attributes of the L-shape. Both carry homogeneous Dirichlet data; they are
// distinguished only so the singular-feature extractor can select the two reentrant faces
// that bound the 3 pi / 2 wedge.
constexpr int MmsReentrantAttribute = 1;
constexpr int MmsOuterAttribute = 2;

// Structured triangulation of the L-shape: three unit squares, each split n x n into
// right triangles. The reentrant corner sits at the origin and is a mesh vertex.
mfem::Mesh MmsLShapeMesh(int n)
{
  struct Key
  {
    int i, j;
    bool operator<(const Key &other) const
    {
      return (i != other.i) ? i < other.i : j < other.j;
    }
  };
  std::map<Key, int> ids;
  std::vector<std::array<double, 2>> coordinates;
  const double h = 1.0 / n;
  const auto vertex = [&](int i, int j)
  {
    const Key key{i, j};
    const auto found = ids.find(key);
    if (found != ids.end())
    {
      return found->second;
    }
    const int id = static_cast<int>(coordinates.size());
    ids.emplace(key, id);
    coordinates.push_back({i * h, j * h});
    return id;
  };
  std::vector<std::array<int, 3>> triangles;
  // Squares (-1,0)x(0,1), (0,1)x(0,1), (-1,0)x(-1,0) in units of the full square.
  const std::array<std::array<int, 2>, 3> origins{{{-n, 0}, {0, 0}, {-n, -n}}};
  for (const auto &origin : origins)
  {
    for (int a = 0; a < n; a++)
    {
      for (int b = 0; b < n; b++)
      {
        const int i = origin[0] + a, j = origin[1] + b;
        const int p00 = vertex(i, j), p10 = vertex(i + 1, j);
        const int p01 = vertex(i, j + 1), p11 = vertex(i + 1, j + 1);
        triangles.push_back({p00, p10, p11});
        triangles.push_back({p00, p11, p01});
      }
    }
  }

  // Boundary segments. The two reentrant faces (the positive x-axis from the origin and
  // the negative y-axis from the origin) get their own attribute so the singular-feature
  // extractor can identify the 3 pi / 2 wedge between them.
  std::vector<std::array<int, 3>> segments;
  for (int a = 0; a < n; a++)
  {
    // Reentrant faces, walking outward from the corner.
    segments.push_back({vertex(a, 0), vertex(a + 1, 0), MmsReentrantAttribute});
    segments.push_back({vertex(0, -a - 1), vertex(0, -a), MmsReentrantAttribute});
    // Outer faces x = -1, x = +1, y = -1, y = +1.
    segments.push_back({vertex(-n, a), vertex(-n, a + 1), MmsOuterAttribute});
    segments.push_back({vertex(-n, -a - 1), vertex(-n, -a), MmsOuterAttribute});
    segments.push_back({vertex(n, a), vertex(n, a + 1), MmsOuterAttribute});
    segments.push_back({vertex(-n + a, -n), vertex(-n + a + 1, -n), MmsOuterAttribute});
    segments.push_back({vertex(-n + a, n), vertex(-n + a + 1, n), MmsOuterAttribute});
    segments.push_back({vertex(a, n), vertex(a + 1, n), MmsOuterAttribute});
  }

  mfem::Mesh mesh(2, static_cast<int>(coordinates.size()),
                  static_cast<int>(triangles.size()), static_cast<int>(segments.size()), 2);
  for (const auto &c : coordinates)
  {
    mesh.AddVertex(c[0], c[1]);
  }
  for (const auto &t : triangles)
  {
    mesh.AddTriangle(t[0], t[1], t[2], 1);
  }
  for (const auto &s : segments)
  {
    mesh.AddBdrSegment(s[0], s[1], s[2]);
  }
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

// Affine tetrahedral extrusion of the L-shape through z in [-1,1]. Each Cartesian cube is
// split by the same Freudenthal diagonal, so all shared quadrilateral faces agree. Only
// vertices touched by the three retained square columns are inserted.
mfem::Mesh MmsExtrudedLShapeMesh(int n, int nz)
{
  REQUIRE(n > 0);
  REQUIRE(nz > 0);
  const int grid_vertices = (2 * n + 1) * (2 * n + 1) * (nz + 1);
  const int active_vertices = grid_vertices - n * n * (nz + 1);
  const int elements = 18 * n * n * nz;
  mfem::Mesh mesh(3, active_vertices, elements, 0, 3);
  std::vector<int> vertex_ids(grid_vertices, -1);
  const auto flat = [n, nz](int i, int j, int k)
  {
    (void)nz;
    return (k * (2 * n + 1) + j) * (2 * n + 1) + i;
  };
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
        // Remove [0,1] x [-1,0], preserving the same 3 pi / 2 cross-section as the 2D MMS.
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

  // Attribute one selects only the two vertical faces bounding the reentrant wedge.
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
    mesh.GetBdrElement(boundary)->SetAttribute(
        (on_x_reentrant || on_y_reentrant) ? MmsReentrantAttribute : MmsOuterAttribute);
  }
  mesh.Finalize(true, false);
  return mesh;
}

double MmsAxial(double z)
{
  return (1.0 - z * z) * (1.0 + 0.25 * z);
}

double MmsAxialDerivative(double z)
{
  return 0.25 - 2.0 * z - 0.75 * z * z;
}

double MmsU3(double x, double y, double z)
{
  return MmsU(x, y) * MmsAxial(z);
}

std::array<double, 3> MmsGradU3(double x, double y, double z)
{
  const auto gradient = MmsGradU(x, y);
  return {MmsAxial(z) * gradient[0], MmsAxial(z) * gradient[1],
          MmsAxialDerivative(z) * MmsU(x, y)};
}

double MmsF3(double x, double y, double z)
{
  return MmsAxial(z) * MmsF(x, y) + (2.0 + 1.5 * z) * MmsU(x, y);
}

// Coefficient wrapper so MFEM can integrate the manufactured source.
class MmsSourceCoefficient : public mfem::Coefficient
{
public:
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    double data[3];
    mfem::Vector point(data, 3);
    T.Transform(ip, point);
    return MmsF(point(0), point(1));
  }
};

class MmsSource3Coefficient : public mfem::Coefficient
{
public:
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    double data[3];
    mfem::Vector point(data, 3);
    T.Transform(ip, point);
    return MmsF3(point(0), point(1), point(2));
  }
};

struct Mms3DSolveReport
{
  long long dofs = 0;
  long long enrichment_dofs = 0;
  long long constrained_enrichment_dofs = 0;
  long long free_enrichment_dofs = 0;
  double algebraic_residual = 0.0;
  double rhs_quadrature_difference = 0.0;
  double error = 0.0;
  double error_quadrature_difference = 0.0;
  double axial_error = 0.0;
  std::vector<double> element_error;
  mfem::Vector standard_coefficients;
  mfem::Vector enrichment_coefficients;
  std::shared_ptr<mfem::SparseMatrix> combined_matrix;
  mfem::Vector combined_rhs;
  mfem::Vector combined_solution;
  std::vector<bool> essential_mask;
};

Mms3DSolveReport MmsSolve3D(mfem::Mesh &mesh, int order, bool enriched)
{
  mfem::H1_FECollection collection(order, 3);
  mfem::FiniteElementSpace space(&mesh, &collection);
  REQUIRE(space.GetVSize() == space.GetTrueVSize());
  const int standard_size = space.GetVSize();

  fem::singular::FeatureTopology features;
  fem::singular::DofTopology topology;
  fem::singular::LocalSparseH1EnrichmentMatrices enrichment;
  int enrichment_size = 0;
  if (enriched)
  {
    features = fem::singular::ExtractSerialSheetFeatures(
        mesh, {MmsReentrantAttribute},
        std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
    REQUIRE(features.features.size() == 1);
    REQUIRE_THAT(features.features.front().nu, WithinRel(MmsNu, 1.0e-12));
    topology = fem::singular::BuildSerialDofTopology(mesh, features, 1);
    enrichment_size = static_cast<int>(topology.h1_dofs.size());
    REQUIRE(enrichment_size > 0);
    const std::vector<fem::singular::IsotropicMaterialCoefficients> materials(mesh.GetNE(),
                                                                              {1.0, 1.0});
    const fem::singular::AdaptiveAssemblyOptions options{8, 5.0e-7, 1.0e-6, 8};
    enrichment = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(topology, space,
                                                                        materials, options);
  }
  const int combined_size = standard_size + enrichment_size;

  mfem::SparseMatrix combined(combined_size, combined_size);
  {
    mfem::BilinearForm standard(&space);
    standard.AddDomainIntegrator(new mfem::DiffusionIntegrator());
    standard.Assemble();
    standard.Finalize();
    const auto &block = standard.SpMat();
    for (int row = 0; row < block.Height(); row++)
    {
      for (int entry = block.GetI()[row]; entry < block.GetI()[row + 1]; entry++)
      {
        combined.Add(row, block.GetJ()[entry], block.GetData()[entry]);
      }
    }
  }
  if (enriched)
  {
    const auto add_block =
        [&combined](const mfem::SparseMatrix &block, int row_offset, int column_offset)
    {
      for (int row = 0; row < block.Height(); row++)
      {
        for (int entry = block.GetI()[row]; entry < block.GetI()[row + 1]; entry++)
        {
          combined.Add(row_offset + row, column_offset + block.GetJ()[entry],
                       block.GetData()[entry]);
        }
      }
    };
    add_block(*enrichment.diffusion.standard_enrichment, 0, standard_size);
    add_block(*enrichment.diffusion.enrichment_standard, standard_size, 0);
    add_block(*enrichment.diffusion.enrichment_enrichment, standard_size, standard_size);
  }
  combined.Finalize();

  mfem::Vector rhs(combined_size);
  rhs = 0.0;
  double rhs_quadrature_difference = 0.0;
  {
    const auto assemble_standard_rhs = [&](int quadrature_order)
    {
      MmsSource3Coefficient source;
      auto *integrator = new mfem::DomainLFIntegrator(source);
      integrator->SetIntRule(
          &mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, quadrature_order));
      mfem::LinearForm standard(&space);
      standard.AddDomainIntegrator(integrator);
      standard.Assemble();
      return mfem::Vector(standard);
    };
    const auto reference = assemble_standard_rhs(16);
    const auto lower_order = assemble_standard_rhs(12);
    mfem::Vector difference(reference);
    difference -= lower_order;
    rhs_quadrature_difference = difference.Norml2() / std::max(reference.Norml2(), 1.0e-30);
    for (int i = 0; i < standard_size; i++)
    {
      rhs(i) = reference(i);
    }
  }
  if (enriched)
  {
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &element_dofs = topology.elements[element].h1;
      if (element_dofs.empty())
      {
        continue;
      }
      auto &T = *mesh.GetElementTransformation(element);
      double jacobian_determinant;
      fem::singular::GetAffineBarycentricGradients(T, jacobian_determinant);
      const auto local = fem::singular::IntegrateReferenceTetrahedronAdaptive(
          8, 1.0e-8, 1.0e-7, 8, element_dofs.size(),
          [&](const fem::singular::BarycentricPoint &lambda, std::vector<double> &value)
          {
            mfem::IntegrationPoint ip;
            ip.Set3(lambda[1], lambda[2], lambda[3]);
            T.SetIntPoint(&ip);
            double data[3];
            mfem::Vector point(data, 3);
            T.Transform(ip, point);
            const double source = MmsF3(point(0), point(1), point(2));
            for (std::size_t i = 0; i < element_dofs.size(); i++)
            {
              value[i] = jacobian_determinant * source *
                         fem::singular::EvaluateHigherOrderGradientPotential(
                             lambda, element_dofs[i].basis);
            }
          });
      REQUIRE(local.converged);
      for (std::size_t i = 0; i < element_dofs.size(); i++)
      {
        rhs(standard_size + static_cast<int>(element_dofs[i].dof)) += local.value[i];
      }
    }
  }

  std::vector<bool> essential(combined_size, false);
  mfem::Array<int> standard_essential;
  space.GetBoundaryTrueDofs(standard_essential);
  for (int dof : standard_essential)
  {
    essential[dof] = true;
  }
  if (enriched)
  {
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), topology);
    REQUIRE(numbering.h1.owned_offset == 0);
    REQUIRE(numbering.h1.owned_size == enrichment_size);
    std::vector<std::array<fem::singular::GlobalVertexId, 3>> boundary_faces;
    mfem::Array<int> vertices;
    for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
    {
      mesh.GetBdrElementVertices(boundary, vertices);
      REQUIRE(vertices.Size() == 3);
      std::array<fem::singular::GlobalVertexId, 3> face{vertices[0], vertices[1],
                                                        vertices[2]};
      std::sort(face.begin(), face.end());
      boundary_faces.push_back(face);
    }
    const auto enrichment_essential = fem::singular::GetEssentialH1TrueDofsOnFaces(
        Mpi::World(), boundary_faces, topology, numbering);
    for (int dof : enrichment_essential)
    {
      essential[standard_size + dof] = true;
    }
  }

  std::vector<int> free_dofs, combined_to_free(combined_size, -1);
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (!essential[dof])
    {
      combined_to_free[dof] = static_cast<int>(free_dofs.size());
      free_dofs.push_back(dof);
    }
  }
  mfem::SparseMatrix system(static_cast<int>(free_dofs.size()),
                            static_cast<int>(free_dofs.size()));
  mfem::Vector reduced_rhs(static_cast<int>(free_dofs.size()));
  for (int row = 0; row < static_cast<int>(free_dofs.size()); row++)
  {
    const int combined_row = free_dofs[row];
    reduced_rhs(row) = rhs(combined_row);
    for (int entry = combined.GetI()[combined_row];
         entry < combined.GetI()[combined_row + 1]; entry++)
    {
      const int column = combined_to_free[combined.GetJ()[entry]];
      if (column >= 0)
      {
        system.Add(row, column, combined.GetData()[entry]);
      }
    }
  }
  system.Finalize();

  mfem::Vector reduced_solution(static_cast<int>(free_dofs.size()));
  reduced_solution = 0.0;
  mfem::GSSmoother preconditioner(system);
  mfem::CGSolver cg;
  cg.SetOperator(system);
  cg.SetPreconditioner(preconditioner);
  cg.SetRelTol(1.0e-14);
  cg.SetAbsTol(1.0e-30);
  cg.SetMaxIter(50000);
  cg.SetPrintLevel(-1);
  cg.Mult(reduced_rhs, reduced_solution);
  REQUIRE(cg.GetConverged());

  mfem::Vector combined_solution(combined_size);
  combined_solution = 0.0;
  for (int row = 0; row < static_cast<int>(free_dofs.size()); row++)
  {
    combined_solution(free_dofs[row]) = reduced_solution(row);
  }
  Mms3DSolveReport report;
  report.dofs = combined_size;
  report.enrichment_dofs = enrichment_size;
  report.rhs_quadrature_difference = rhs_quadrature_difference;
  for (int dof = 0; dof < enrichment_size; dof++)
  {
    report.constrained_enrichment_dofs += essential[standard_size + dof];
    report.free_enrichment_dofs += !essential[standard_size + dof];
  }
  {
    mfem::Vector algebraic(static_cast<int>(free_dofs.size()));
    system.Mult(reduced_solution, algebraic);
    algebraic -= reduced_rhs;
    report.algebraic_residual =
        algebraic.Norml2() / std::max(reduced_rhs.Norml2(), 1.0e-30);
  }

  mfem::GridFunction standard_solution(&space);
  for (int i = 0; i < standard_size; i++)
  {
    standard_solution(i) = combined_solution(i);
  }
  report.element_error.assign(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double jacobian_determinant;
    const auto grad_lambda =
        fem::singular::GetAffineBarycentricGradients(T, jacobian_determinant);
    const auto error_at = [&](const fem::singular::BarycentricPoint &lambda)
    {
      mfem::IntegrationPoint ip;
      ip.Set3(lambda[1], lambda[2], lambda[3]);
      T.SetIntPoint(&ip);
      double data[3];
      mfem::Vector point(data, 3), gradient(3);
      T.Transform(ip, point);
      standard_solution.GetGradient(T, gradient);
      if (enriched)
      {
        for (const auto &dof : topology.elements[element].h1)
        {
          const double coefficient =
              combined_solution(standard_size + static_cast<int>(dof.dof));
          const auto basis =
              fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, dof.basis);
          for (int d = 0; d < 3; d++)
          {
            gradient(d) += coefficient * basis.value[d];
          }
        }
      }
      const auto exact = MmsGradU3(point(0), point(1), point(2));
      const double dx = exact[0] - gradient(0), dy = exact[1] - gradient(1);
      const double dz = exact[2] - gradient(2);
      return std::array<double, 2>{jacobian_determinant * (dx * dx + dy * dy + dz * dz),
                                   jacobian_determinant * dz * dz};
    };

    // Feature-aligned Duffy quadrature resolves r^{-2/3} exactly enough for the direct
    // energy error. Tetrahedra containing only a feature endpoint use the node map; those
    // containing a feature segment use the edge map. Ordinary elements need no grading.
    mfem::Array<int> element_vertices;
    mesh.GetElementVertices(element, element_vertices);
    std::vector<int> singular_nodes;
    for (int local = 0; local < element_vertices.Size(); local++)
    {
      const double *point = mesh.GetVertex(element_vertices[local]);
      if (std::hypot(point[0], point[1]) < 1.0e-13)
      {
        singular_nodes.push_back(local);
      }
    }
    const auto integrate_component = [&](int component, int quadrature_order)
    {
      const auto integrand = [&](const fem::singular::BarycentricPoint &lambda)
      { return error_at(lambda)[component]; };
      if (singular_nodes.size() >= 2)
      {
        return fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
            quadrature_order, singular_nodes[0], singular_nodes[1], 6.0, integrand);
      }
      if (singular_nodes.size() == 1)
      {
        return fem::singular::IntegrateReferenceTetrahedronNodeDuffy(
            quadrature_order, singular_nodes[0], 6.0, integrand);
      }
      return fem::singular::IntegrateReferenceTetrahedron(quadrature_order, 1, integrand);
    };
    const double error = integrate_component(0, 20);
    const double lower_order_error = integrate_component(0, 18);
    const double axial_error = integrate_component(1, 20);
    report.element_error[element] = error;
    report.error += error;
    report.error_quadrature_difference += std::abs(error - lower_order_error);
    report.axial_error += axial_error;
  }
  report.standard_coefficients.SetSize(standard_size);
  for (int i = 0; i < standard_size; i++)
  {
    report.standard_coefficients(i) = combined_solution(i);
  }
  report.enrichment_coefficients.SetSize(enrichment_size);
  for (int i = 0; i < enrichment_size; i++)
  {
    report.enrichment_coefficients(i) = combined_solution(standard_size + i);
  }
  report.combined_matrix = std::make_shared<mfem::SparseMatrix>(combined);
  report.combined_rhs = rhs;
  report.combined_solution = combined_solution;
  report.essential_mask = essential;
  return report;
}

std::vector<double> Mms3DCorrectionIndicator(mfem::Mesh &mesh, int order,
                                             const mfem::Vector &correction,
                                             int enrichment_size)
{
  mfem::H1_FECollection collection(order, 3);
  mfem::FiniteElementSpace space(&mesh, &collection);
  const int standard_size = space.GetVSize();
  REQUIRE(correction.Size() == standard_size + enrichment_size);
  const auto features = fem::singular::ExtractSerialSheetFeatures(
      mesh, {MmsReentrantAttribute},
      std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  const auto topology = fem::singular::BuildSerialDofTopology(mesh, features, 1);
  REQUIRE(static_cast<int>(topology.h1_dofs.size()) == enrichment_size);
  mfem::GridFunction standard(&space);
  for (int i = 0; i < standard_size; i++)
  {
    standard(i) = correction(i);
  }

  std::vector<double> indicator(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double jacobian_determinant;
    const auto grad_lambda =
        fem::singular::GetAffineBarycentricGradients(T, jacobian_determinant);
    const auto energy = [&](const fem::singular::BarycentricPoint &lambda)
    {
      mfem::IntegrationPoint ip;
      ip.Set3(lambda[1], lambda[2], lambda[3]);
      T.SetIntPoint(&ip);
      mfem::Vector gradient(3);
      standard.GetGradient(T, gradient);
      for (const auto &dof : topology.elements[element].h1)
      {
        const double coefficient = correction(standard_size + static_cast<int>(dof.dof));
        const auto basis =
            fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, dof.basis);
        for (int d = 0; d < 3; d++)
        {
          gradient(d) += coefficient * basis.value[d];
        }
      }
      return jacobian_determinant * mfem::InnerProduct(gradient, gradient);
    };
    mfem::Array<int> vertices;
    mesh.GetElementVertices(element, vertices);
    std::vector<int> singular_nodes;
    for (int local = 0; local < vertices.Size(); local++)
    {
      const double *point = mesh.GetVertex(vertices[local]);
      if (std::hypot(point[0], point[1]) < 1.0e-13)
      {
        singular_nodes.push_back(local);
      }
    }
    if (singular_nodes.size() >= 2)
    {
      indicator[element] = fem::singular::IntegrateReferenceTetrahedronEdgeDuffy(
          18, singular_nodes[0], singular_nodes[1], 6.0, energy);
    }
    else if (singular_nodes.size() == 1)
    {
      indicator[element] = fem::singular::IntegrateReferenceTetrahedronNodeDuffy(
          18, singular_nodes[0], 6.0, energy);
    }
    else
    {
      indicator[element] = fem::singular::IntegrateReferenceTetrahedron(14, 1, energy);
    }
  }
  return indicator;
}

struct Mms3DHierarchicalReport
{
  std::vector<double> element_error;
  std::vector<double> indicator_global;
  std::vector<double> indicator_local;
  long long coarse_dofs = 0;
  double coarse_error = 0.0;
  double global_energy = 0.0;
  double local_energy = 0.0;
  double captured_fraction = 0.0;
  double global_residual = 0.0;
  double maximum_patch_residual = 0.0;
  double maximum_patch_condition = 0.0;
  double overlap_error = 0.0;
  int edge_patches = 0;
  int face_patches = 0;
  int interior_patches = 0;
  int owned_modes = 0;
  int maximum_patch_dimension = 0;
};

Mms3DHierarchicalReport Mms3DHierarchicalBenchmark(mfem::Mesh &mesh, int order)
{
  Mms3DHierarchicalReport report;
  const auto coarse = MmsSolve3D(mesh, order, true);
  const auto fine = MmsSolve3D(mesh, order + 1, true);
  report.element_error = coarse.element_error;
  report.coarse_dofs = coarse.dofs;
  report.coarse_error = coarse.error;
  REQUIRE(coarse.enrichment_dofs == fine.enrichment_dofs);
  const int enrichment = static_cast<int>(coarse.enrichment_dofs);
  mfem::H1_FECollection coarse_collection(order, 3), fine_collection(order + 1, 3);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  const int fine_combined = fine_standard + enrichment;

  mfem::DenseMatrix injection(fine_combined, coarse_standard + enrichment);
  injection = 0.0;
  {
    mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
    mfem::GridFunction unit(&coarse_space), image(&fine_space);
    for (int column = 0; column < coarse_standard; column++)
    {
      unit = 0.0;
      unit(column) = 1.0;
      transfer.Mult(unit, image);
      for (int row = 0; row < fine_standard; row++)
      {
        injection(row, column) = image(row);
      }
    }
    for (int e = 0; e < enrichment; e++)
    {
      injection(fine_standard + e, coarse_standard + e) = 1.0;
    }
  }
  mfem::Vector coarse_combined(coarse_standard + enrichment);
  for (int i = 0; i < coarse_standard; i++)
  {
    coarse_combined(i) = coarse.standard_coefficients(i);
  }
  for (int e = 0; e < enrichment; e++)
  {
    coarse_combined(coarse_standard + e) = coarse.enrichment_coefficients(e);
  }
  mfem::Vector injected(fine_combined);
  injection.Mult(coarse_combined, injected);
  mfem::Vector residual(fine.combined_rhs);
  mfem::Vector applied(fine_combined);
  fine.combined_matrix->Mult(injected, applied);
  residual -= applied;
  mfem::Vector global(fine.combined_solution);
  global -= injected;
  fine.combined_matrix->Mult(global, applied);
  report.global_energy = mfem::InnerProduct(global, applied);
  mfem::Vector global_residual(applied);
  global_residual -= residual;
  for (int i = 0; i < fine_combined; i++)
  {
    if (fine.essential_mask[i])
    {
      global_residual(i) = 0.0;
    }
  }
  double free_residual_norm2 = 0.0;
  for (int i = 0; i < fine_combined; i++)
  {
    if (!fine.essential_mask[i])
    {
      free_residual_norm2 += residual(i) * residual(i);
    }
  }
  report.global_residual =
      global_residual.Norml2() / std::max(std::sqrt(free_residual_norm2), 1.0e-30);
  report.indicator_global = Mms3DCorrectionIndicator(mesh, order + 1, global, enrichment);

  const auto features = fem::singular::ExtractSerialSheetFeatures(
      mesh, {MmsReentrantAttribute},
      std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  const auto topology = fem::singular::BuildSerialDofTopology(mesh, features, 1);
  struct Patch
  {
    int owned = 0;
    std::vector<mfem::Vector> basis;
    std::vector<int> coarse_guests;
    std::vector<int> singular_guests;
    mfem::Vector coefficients;
  };
  std::vector<Patch> patches;

  const auto entity_complement =
      [&](mfem::Array<int> fine_entity, mfem::Array<int> coarse_entity)
  {
    mfem::FiniteElementSpace::AdjustVDofs(fine_entity);
    mfem::FiniteElementSpace::AdjustVDofs(coarse_entity);
    std::vector<int> rows, columns;
    for (int dof : fine_entity)
    {
      if (!fine.essential_mask[dof])
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
      for (int row = 0; row < vector.Size(); row++)
      {
        vector(row) = injection(rows[row], column);
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
    std::vector<mfem::Vector> expanded;
    for (const auto &local : complement)
    {
      mfem::Vector vector(fine_combined);
      vector = 0.0;
      for (int row = 0; row < static_cast<int>(rows.size()); row++)
      {
        vector(rows[row]) = local(row);
      }
      expanded.push_back(vector);
    }
    return expanded;
  };

  const auto add_patch = [&](std::vector<mfem::Vector> owned,
                             const std::vector<int> &elements, int entity_type)
  {
    if (owned.empty())
    {
      return;
    }
    Patch patch;
    patch.owned = static_cast<int>(owned.size());
    patch.basis = std::move(owned);
    report.owned_modes += patch.owned;
    report.edge_patches += entity_type == 0;
    report.face_patches += entity_type == 1;
    report.interior_patches += entity_type == 2;
    std::set<int> coarse_guests, singular_guests;
    mfem::Array<int> element_dofs;
    for (int element : elements)
    {
      coarse_space.GetElementDofs(element, element_dofs);
      mfem::FiniteElementSpace::AdjustVDofs(element_dofs);
      for (int dof : element_dofs)
      {
        if (!coarse.essential_mask[dof])
        {
          coarse_guests.insert(dof);
        }
      }
      for (const auto &dof : topology.elements[element].h1)
      {
        if (!fine.essential_mask[fine_standard + static_cast<int>(dof.dof)])
        {
          singular_guests.insert(static_cast<int>(dof.dof));
        }
      }
    }
    for (int dof : coarse_guests)
    {
      mfem::Vector column(fine_combined);
      for (int row = 0; row < fine_combined; row++)
      {
        column(row) = injection(row, dof);
      }
      patch.coarse_guests.push_back(dof);
      patch.basis.push_back(column);
    }
    for (int dof : singular_guests)
    {
      mfem::Vector column(fine_combined);
      column = 0.0;
      column(fine_standard + dof) = 1.0;
      patch.singular_guests.push_back(dof);
      patch.basis.push_back(column);
    }
    const int size = static_cast<int>(patch.basis.size());
    report.maximum_patch_dimension = std::max(report.maximum_patch_dimension, size);
    mfem::DenseMatrix matrix(size);
    mfem::Vector patch_rhs(size);
    std::vector<mfem::Vector> images(size);
    for (int column = 0; column < size; column++)
    {
      images[column].SetSize(fine_combined);
      fine.combined_matrix->Mult(patch.basis[column], images[column]);
      patch_rhs(column) = mfem::InnerProduct(patch.basis[column], residual);
    }
    for (int row = 0; row < size; row++)
    {
      for (int column = 0; column < size; column++)
      {
        matrix(row, column) = mfem::InnerProduct(patch.basis[row], images[column]);
      }
    }
    mfem::DenseMatrixInverse inverse(matrix, true);
    mfem::DenseMatrix inverse_matrix;
    inverse.GetInverseMatrix(inverse_matrix);
    report.maximum_patch_condition =
        std::max(report.maximum_patch_condition, matrix.FNorm() * inverse_matrix.FNorm());
    patch.coefficients.SetSize(size);
    inverse.Mult(patch_rhs, patch.coefficients);
    mfem::Vector check(size);
    matrix.Mult(patch.coefficients, check);
    check -= patch_rhs;
    report.maximum_patch_residual =
        std::max(report.maximum_patch_residual,
                 check.Norml2() / std::max(patch_rhs.Norml2(), 1.0e-30));
    patches.push_back(std::move(patch));
  };

  std::vector<std::vector<int>> edge_elements(mesh.GetNEdges());
  std::vector<std::vector<int>> face_elements(mesh.GetNumFaces());
  mfem::Array<int> entities, orientations;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementEdges(element, entities, orientations);
    for (int edge : entities)
    {
      edge_elements[edge].push_back(element);
    }
    mesh.GetElementFaces(element, entities, orientations);
    for (int face : entities)
    {
      face_elements[face].push_back(element);
    }
  }
  for (int edge = 0; edge < mesh.GetNEdges(); edge++)
  {
    mfem::Array<int> fine_dofs, coarse_dofs;
    fine_space.GetEdgeInteriorDofs(edge, fine_dofs);
    coarse_space.GetEdgeInteriorDofs(edge, coarse_dofs);
    add_patch(entity_complement(fine_dofs, coarse_dofs), edge_elements[edge], 0);
  }
  for (int face = 0; face < mesh.GetNumFaces(); face++)
  {
    mfem::Array<int> fine_dofs, coarse_dofs;
    fine_space.GetFaceInteriorDofs(face, fine_dofs);
    coarse_space.GetFaceInteriorDofs(face, coarse_dofs);
    add_patch(entity_complement(fine_dofs, coarse_dofs), face_elements[face], 1);
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::Array<int> fine_dofs, coarse_dofs;
    fine_space.GetElementInteriorDofs(element, fine_dofs);
    coarse_space.GetElementInteriorDofs(element, coarse_dofs);
    add_patch(entity_complement(fine_dofs, coarse_dofs), {element}, 2);
  }
  int coarse_free = 0, fine_free = 0;
  for (int i = 0; i < coarse_standard; i++)
  {
    coarse_free += !coarse.essential_mask[i];
  }
  for (int i = 0; i < fine_standard; i++)
  {
    fine_free += !fine.essential_mask[i];
  }
  REQUIRE(report.owned_modes == fine_free - coarse_free);

  mfem::Vector raw(fine_combined);
  raw = 0.0;
  std::vector<double> coarse_sum(coarse_standard, 0.0), singular_sum(enrichment, 0.0);
  std::vector<int> coarse_count(coarse_standard, 0), singular_count(enrichment, 0);
  for (const auto &patch : patches)
  {
    for (int i = 0; i < patch.owned; i++)
    {
      raw.Add(patch.coefficients(i), patch.basis[i]);
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
    }
    if (coarse_count[dof] == 0)
    {
      continue;
    }
    const double coefficient = coarse_sum[dof] / coarse_count[dof];
    for (int row = 0; row < fine_combined; row++)
    {
      raw(row) += coefficient * injection(row, dof);
    }
    double weight_sum = 0.0;
    for (const auto &patch : patches)
    {
      weight_sum += static_cast<double>(std::count(patch.coarse_guests.begin(),
                                                   patch.coarse_guests.end(), dof)) /
                    coarse_count[dof];
    }
    report.overlap_error = std::max(report.overlap_error, std::abs(weight_sum - 1.0));
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    if (!coarse.essential_mask[coarse_standard + dof])
    {
      REQUIRE(singular_count[dof] > 0);
    }
    if (singular_count[dof] > 0)
    {
      raw(fine_standard + dof) += singular_sum[dof] / singular_count[dof];
      double weight_sum = 0.0;
      for (const auto &patch : patches)
      {
        weight_sum += static_cast<double>(std::count(patch.singular_guests.begin(),
                                                     patch.singular_guests.end(), dof)) /
                      singular_count[dof];
      }
      report.overlap_error = std::max(report.overlap_error, std::abs(weight_sum - 1.0));
    }
  }
  mfem::Vector raw_applied(fine_combined);
  fine.combined_matrix->Mult(raw, raw_applied);
  const double raw_energy = mfem::InnerProduct(raw, raw_applied);
  const double alpha =
      raw_energy > 0.0 ? mfem::InnerProduct(raw, residual) / raw_energy : 0.0;
  raw *= alpha;
  fine.combined_matrix->Mult(raw, raw_applied);
  report.local_energy = mfem::InnerProduct(raw, raw_applied);
  report.captured_fraction = report.local_energy / report.global_energy;
  report.indicator_local = Mms3DCorrectionIndicator(mesh, order + 1, raw, enrichment);
  return report;
}

}  // namespace

TEST_CASE("Manufactured L-shape source and solution are self-consistent",
          "[singularmms][singularelements][Serial]")
{
  // Stage 0. Before any assembly, confirm the analytic pieces used by the harness agree
  // with the specification: s is harmonic, u vanishes on every face, and the reference
  // energy is reproduced by an independent in-test quadrature of int f u.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);

  // grad(s) must satisfy |grad s|^2 = nu^2 r^(2 nu - 2) exactly.
  for (const auto &p : std::vector<std::array<double, 2>>{
           {-0.5, 0.5}, {0.4, 0.6}, {-0.7, -0.3}, {0.2, 0.9}, {-0.9, 0.1}})
  {
    const auto gs = MmsGradS(p[0], p[1]);
    const double r = std::hypot(p[0], p[1]);
    CHECK_THAT(gs[0] * gs[0] + gs[1] * gs[1],
               WithinRel(MmsNu * MmsNu * std::pow(r, 2.0 * MmsNu - 2.0), 1.0e-12));
  }

  // u vanishes on all six faces: the four outer ones through B, the two reentrant ones
  // through s itself.
  for (const auto &p : std::vector<std::array<double, 2>>{
           {-1.0, 0.3}, {1.0, 0.5}, {0.3, 1.0}, {-0.4, -1.0}, {0.0, -0.5}, {0.5, 0.0}})
  {
    CHECK_THAT(MmsU(p[0], p[1]), WithinAbs(0.0, 1.0e-14));
  }

  // Independent in-test evaluation of int f u by corner-graded quadrature on the three
  // unit squares. This is a third implementation of the reference, cheaper and lower
  // order than the two in SPEC.md, so it is only asked to agree to 1e-6.
  const auto graded_square = [](double sx, double sy, int order, int levels)
  {
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::SQUARE, order);
    double total = 0.0, a = 1.0;
    for (int level = 0; level < levels; level++)
    {
      const double b = (level + 1 < levels) ? 0.5 * a : 0.0;
      const std::array<std::array<double, 4>, 2> shells{{{b, a, 0.0, b}, {0.0, a, b, a}}};
      for (const auto &shell : shells)
      {
        const double du = shell[1] - shell[0], dv = shell[3] - shell[2];
        if (!(du > 0.0) || !(dv > 0.0))
        {
          continue;
        }
        for (int q = 0; q < rule.GetNPoints(); q++)
        {
          const auto &ip = rule.IntPoint(q);
          const double uu = shell[0] + du * ip.x, vv = shell[2] + dv * ip.y;
          const double x = sx * uu, y = sy * vv;
          total += ip.weight * du * dv * MmsF(x, y) * MmsU(x, y);
        }
      }
      a = b;
    }
    return total;
  };
  const double energy = graded_square(1.0, 1.0, 20, 30) + graded_square(-1.0, 1.0, 20, 30) +
                        graded_square(-1.0, -1.0, 20, 30);
  CAPTURE(energy, MmsReferenceEnergy);
  CHECK_THAT(energy, WithinRel(MmsReferenceEnergy, 1.0e-6));
}

TEST_CASE("Manufactured extruded L-shape is a genuine 3D singular-edge MMS",
          "[singularmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(1, 2);
  constexpr int nz = 2;
  CAPTURE(n, nz);
  auto mesh = MmsExtrudedLShapeMesh(n, nz);
  CHECK(mesh.Dimension() == 3);
  CHECK(mesh.SpaceDimension() == 3);
  CHECK(mesh.Conforming());
  CHECK(mesh.GetNE() == 18 * n * n * nz);
  CHECK_THAT(mesh.GetElementVolume(0), WithinAbs(1.0 / (3.0 * n * n * nz), 1.0e-13));
  double volume = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    CHECK(mesh.GetElementBaseGeometry(element) == mfem::Geometry::TETRAHEDRON);
    CHECK(mesh.GetElementVolume(element) > 0.0);
    volume += mesh.GetElementVolume(element);
  }
  CHECK_THAT(volume, WithinAbs(6.0, 1.0e-12));

  int reentrant_faces = 0, outer_faces = 0;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    const int attribute = mesh.GetBdrAttribute(boundary);
    reentrant_faces += attribute == MmsReentrantAttribute;
    outer_faces += attribute == MmsOuterAttribute;
  }
  CHECK(reentrant_faces == 4 * n * nz);
  CHECK(outer_faces > 0);

  const auto features = fem::singular::ExtractSerialSheetFeatures(
      mesh, {MmsReentrantAttribute},
      std::vector<fem::singular::TriangleMaterial>{{1, 1.0}});
  REQUIRE(features.features.size() == 1);
  const auto &feature = features.features.front();
  CHECK_FALSE(feature.closed);
  CHECK_THAT(feature.nu, WithinRel(MmsNu, 1.0e-12));
  CHECK(feature.segments.size() == static_cast<std::size_t>(nz));
  double feature_length = 0.0;
  for (std::size_t segment_index : feature.segments)
  {
    REQUIRE(segment_index < features.segments.size());
    const auto &segment = features.segments[segment_index];
    CHECK(segment.type == fem::singular::FeatureSegmentType::TRANSMISSION_WEDGE);
    const double *a = mesh.GetVertex(segment.mesh_vertices[0]);
    const double *b = mesh.GetVertex(segment.mesh_vertices[1]);
    CHECK(std::hypot(a[0], a[1]) < 1.0e-13);
    CHECK(std::hypot(b[0], b[1]) < 1.0e-13);
    feature_length += std::abs(b[2] - a[2]);
  }
  CHECK_THAT(feature_length, WithinAbs(2.0, 1.0e-13));
  const auto topology = fem::singular::BuildSerialDofTopology(mesh, features, 1);
  CHECK_FALSE(topology.h1_dofs.empty());

  // The axial factor makes this more than a 2D solve embedded in a 3D mesh: it vanishes on
  // the caps but has a nonzero z-gradient and the exact source includes -Z'' v.
  for (const auto &point : std::vector<std::array<double, 3>>{{-0.4, 0.3, -1.0},
                                                              {-0.4, 0.3, 1.0},
                                                              {0.0, -0.5, 0.2},
                                                              {0.5, 0.0, -0.3},
                                                              {-1.0, 0.2, 0.4},
                                                              {0.2, 1.0, -0.4}})
  {
    CHECK_THAT(MmsU3(point[0], point[1], point[2]), WithinAbs(0.0, 1.0e-13));
  }
  const auto gradient = MmsGradU3(-0.4, 0.3, 0.2);
  CHECK(std::abs(gradient[2]) > 1.0e-4);

  // Independent finite-difference check of -Delta u = f away from the singular edge.
  constexpr double h = 2.0e-4;
  for (const auto &point : std::vector<std::array<double, 3>>{
           {-0.4, 0.3, 0.2}, {-0.7, -0.4, -0.3}, {0.4, 0.6, 0.5}})
  {
    double laplacian = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
      auto plus = point, minus = point;
      plus[axis] += h;
      minus[axis] -= h;
      laplacian +=
          (MmsU3(plus[0], plus[1], plus[2]) - 2.0 * MmsU3(point[0], point[1], point[2]) +
           MmsU3(minus[0], minus[1], minus[2])) /
          (h * h);
    }
    CHECK_THAT(-laplacian, WithinAbs(MmsF3(point[0], point[1], point[2]), 2.0e-6));
  }

  const auto &line_rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 12);
  double axial_mass = 0.0, axial_energy = 0.0;
  for (int q = 0; q < line_rule.GetNPoints(); q++)
  {
    const auto &ip = line_rule.IntPoint(q);
    const double z = 2.0 * ip.x - 1.0;
    axial_mass += 2.0 * ip.weight * MmsAxial(z) * MmsAxial(z);
    axial_energy += 2.0 * ip.weight * MmsAxialDerivative(z) * MmsAxialDerivative(z);
  }
  CHECK_THAT(axial_mass, WithinRel(113.0 / 105.0, 1.0e-13));
  CHECK_THAT(axial_energy, WithinRel(83.0 / 30.0, 1.0e-13));
}

TEST_CASE("Manufactured extruded L-shape singular electrostatics improves the 3D solve",
          "[singularmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  auto mesh = MmsExtrudedLShapeMesh(1, 2);
  const auto standard = MmsSolve3D(mesh, order, false);
  const auto enriched = MmsSolve3D(mesh, order, true);
  CAPTURE(standard.dofs, enriched.dofs, enriched.enrichment_dofs,
          enriched.constrained_enrichment_dofs, enriched.free_enrichment_dofs,
          standard.algebraic_residual, enriched.algebraic_residual,
          standard.rhs_quadrature_difference, enriched.rhs_quadrature_difference,
          standard.error, enriched.error, standard.error_quadrature_difference,
          enriched.error_quadrature_difference, standard.axial_error, enriched.axial_error);
  CHECK(standard.algebraic_residual < 1.0e-10);
  CHECK(enriched.algebraic_residual < 1.0e-10);
  CHECK(enriched.enrichment_dofs > 0);
  CHECK(enriched.constrained_enrichment_dofs > 0);
  CHECK(enriched.free_enrichment_dofs > 0);
  CHECK(enriched.constrained_enrichment_dofs + enriched.free_enrichment_dofs ==
        enriched.enrichment_dofs);
  CHECK(standard.rhs_quadrature_difference < 1.0e-4);
  CHECK(enriched.rhs_quadrature_difference < 1.0e-4);
  // Orders 18 and 20 agree to one percent, while enrichment reduces the error by at
  // least 25 percent. The improvement is therefore far larger than quadrature uncertainty.
  CHECK(standard.error_quadrature_difference < 1.0e-2 * standard.error);
  CHECK(enriched.error_quadrature_difference < 1.0e-2 * enriched.error);
  CHECK(enriched.error < 0.75 * standard.error);
  CHECK(enriched.axial_error > 0.01 * enriched.error);
}

TEST_CASE("Manufactured L-shape mesh resolves the reentrant corner",
          "[singularmms][singularelements][Serial]")
{
  // Stage 1a. The mesh must place a vertex exactly at the reentrant corner, exclude the
  // missing quadrant, and contain no element straddling it.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4, 8);
  CAPTURE(n);
  auto mesh = MmsLShapeMesh(n);

  CHECK(mesh.Dimension() == 2);
  CHECK(mesh.SpaceDimension() == 2);
  CHECK(mesh.GetNE() == 6 * n * n);

  // A vertex sits at the origin.
  int corner = -1;
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    const double *c = mesh.GetVertex(v);
    if (std::hypot(c[0], c[1]) < 1.0e-14)
    {
      corner = v;
    }
  }
  CHECK(corner >= 0);

  // Every element centroid lies inside the L-shape, i.e. never in the removed quadrant.
  mfem::Vector centroid;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementCenter(element, centroid);
    const double x = centroid(0), y = centroid(1);
    CAPTURE(element, x, y);
    CHECK(x >= -1.0 - 1.0e-12);
    CHECK(x <= 1.0 + 1.0e-12);
    CHECK(y >= -1.0 - 1.0e-12);
    CHECK(y <= 1.0 + 1.0e-12);
    CHECK_FALSE((x > 1.0e-12 && y < -1.0e-12));
  }

  // Total area is three unit squares.
  double area = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    area += mesh.GetElementVolume(element);
  }
  CHECK_THAT(area, WithinRel(3.0, 1.0e-12));
}

namespace
{

// Everything the staged certification reports for one solve.
struct MmsSolveReport
{
  int elements = 0;
  long long dofs = 0;
  long long enrichment_dofs = 0;

  // Stage 2 algebra checks.
  double algebraic_residual = 0.0;  // ||K x - b|| / ||b||
  double a_u_uh = 0.0;              // a(u, u_h) by quadrature against the exact gradient
  double l_uh = 0.0;                // l(u_h) = int f u_h
  double a_uh_uh = 0.0;             // x^T K x

  // Stage 3 error expressions, which must agree.
  double error_direct = 0.0;    // sum_K int_K |grad u - grad u_h|^2
  double error_expanded = 0.0;  // a(u,u) - 2 a(u,u_h) + a(u_h,u_h)
  double error_galerkin = 0.0;  // a(u,u) - a(u_h,u_h)

  std::vector<double> element_error;  // e_K^2, for later ranking metrics
  double corner_error = 0.0;          // e_K^2 summed over corner-incident elements only

  // Stage 6 recovered-flux indicators, eta_K^2. Two variants, differing ONLY in whether the
  // singular enrichment contributes to the field handed to the estimator.
  std::vector<double> indicator_full;      // enrichment included
  std::vector<double> indicator_standard;  // enrichment sliced off, as production does

  // Naive "augmented" recovery: recover only the standard field, then re-add the singular
  // flux carried at its EXISTING finite-element coefficient. Computed solely to demonstrate
  // that this is algebraically identical to indicator_standard and therefore useless as a
  // candidate improvement.
  std::vector<double> indicator_naive_augmented;

  // Free-coefficient RT (+) S_flux recovery: the singular amplitude is fitted INDEPENDENTLY
  // of the solution's own enrichment coefficient, which is what makes it differ from the
  // naive add-back.
  std::vector<double> indicator_augmented;

  // rho = |<r, s_perp>|^2 / (||r||^2 ||s_perp||^2), the fraction of the RECOVERY RESIDUAL
  // that one free singular mode can remove. This, not ||s_perp||^2/||s||^2, is the quantity
  // that decides whether augmentation can change what the estimator sees.
  double augmentation_rho = 0.0;
  double residual_energy_before = 0.0;
  double residual_energy_after = 0.0;

  // Retained so a two-level hierarchical benchmark can inject this solution into a
  // higher-order space and form delta u there. The standard block is the H1 grid-function
  // vector; the enrichment block is the canonical singular coefficient vector.
  mfem::Vector standard_coefficients;
  mfem::Vector enrichment_block;
  int standard_order = 0;

  // Full combined system, retained so a hierarchical ablation can form the true residual
  // r_f = b_f - A_f P_c u_c INCLUDING the singular coupling block A_se u_s. Restricting the
  // residual to the standard block alone silently drops that coupling.
  std::shared_ptr<mfem::SparseMatrix> combined_matrix;
  mfem::Vector combined_rhs;
  mfem::Vector combined_solution_vector;
  std::vector<bool> essential_mask;
};

// How the corner-incident elements are integrated. Both schemes resolve the algebraic
// r^(-1/3) gradient, which no plain rule can do at any order, and they are independent of
// one another. Sweeping both is a check on the quadrature itself.
enum class MmsCornerQuadrature
{
  // Geometrically graded self-similar shells, each integrated directly. Mirrors the
  // reference quadrature in SPEC.md.
  GRADED,
  // Palace's production singular rule: the node Duffy map with radial power 6. Never
  // places a point at the singular node, which is what the enriched space requires.
  DUFFY
};

// Which local vertex of an element, if any, sits at the reentrant corner.
int MmsCornerVertex(const mfem::Mesh &mesh, int element)
{
  mfem::Array<int> vertices;
  mesh.GetElementVertices(element, vertices);
  for (int local = 0; local < vertices.Size(); local++)
  {
    const double *c = mesh.GetVertex(vertices[local]);
    if (std::hypot(c[0], c[1]) < 1.0e-13)
    {
      return local;
    }
  }
  return -1;
}

// Visit reference-triangle quadrature points for one element, choosing a corner-resolving
// rule when the element touches the reentrant corner. Weights are reference weights, so
// they still need the geometric Jacobian.
//
// The visitor receives the full barycentric TRIPLE, never just (xi, eta). That matters:
// under the radial-power-6 Duffy map the singular coordinate reaches ~1e-18, so if the
// third coordinate is rebuilt as 1 - xi - eta it is destroyed by cancellation whenever the
// singular node is local vertex 1 or 2 (rho comes out exactly zero, or silently wrong by
// three orders of magnitude). This is the same reason
// fem::singular::TriangleNodeRadialCoordinate sums the complementary coordinates instead of
// subtracting a nearly unit one.
template <typename Visitor>
void MmsForEachQuadraturePoint(const mfem::Mesh &mesh, int element, int quadrature_order,
                               MmsCornerQuadrature scheme, const Visitor &visit)
{
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, quadrature_order);
  const int corner = MmsCornerVertex(mesh, element);
  if (corner < 0)
  {
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &ip = rule.IntPoint(q);
      visit(fem::singular::TriangleBarycentricPoint{1.0 - ip.x - ip.y, ip.x, ip.y},
            ip.weight);
    }
    return;
  }
  if (scheme == MmsCornerQuadrature::DUFFY)
  {
    // Palace's own assembly maps a requested order q to 2q + 15 before building this rule,
    // so the harness applies the same mapping. Duffy weights already integrate over the
    // reference triangle of area 1/2, matching MFEM's convention.
    fem::singular::ForEachReferenceTriangleNodeDuffyQuadraturePoint(
        2 * quadrature_order + 15, corner, fem::singular::TriangleDuffyRadialPower, visit);
    return;
  }

  // Reference-triangle vertices, with the corner one identified.
  const std::array<std::array<double, 2>, 3> reference{
      {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}};
  const auto &c = reference[corner];
  const auto &p1 = reference[(corner + 1) % 3];
  const auto &p2 = reference[(corner + 2) % 3];
  // Point on the edge c->p at parameter t, i.e. the scaled similar triangle's vertex.
  const auto along = [&c](const std::array<double, 2> &p, double t)
  { return std::array<double, 2>{c[0] + t * (p[0] - c[0]), c[1] + t * (p[1] - c[1])}; };
  // Integrate a reference-space triangle (A,B,C) with the plain rule, weighting by its
  // area relative to the reference triangle (whose area is 1/2).
  const auto integrate_triangle = [&](const std::array<double, 2> &A,
                                      const std::array<double, 2> &B,
                                      const std::array<double, 2> &C)
  {
    const double jac = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    const double factor = std::abs(jac);  // reference weights already carry the 1/2
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &ip = rule.IntPoint(q);
      const double l1 = ip.x, l2 = ip.y, l0 = 1.0 - l1 - l2;
      const double xi = l0 * A[0] + l1 * B[0] + l2 * C[0];
      const double eta = l0 * A[1] + l1 * B[1] + l2 * C[1];
      // The graded scheme's innermost shell only reaches 2^-27, where rebuilding the third
      // barycentric coordinate is still accurate; it is written this way for uniformity
      // with the Duffy branch.
      visit(fem::singular::TriangleBarycentricPoint{1.0 - xi - eta, xi, eta},
            ip.weight * factor);
    }
  };
  // Geometric shells toward the corner. Each shell between scales b < a is the
  // quadrilateral (b*p1, a*p1, a*p2, b*p2), split into two triangles and integrated
  // DIRECTLY. Integrating T(a) minus T(b) instead would telescope to the plain rule
  // over the whole element and gain nothing.
  constexpr int levels = 28;
  double a_scale = 1.0;
  for (int level = 0; level < levels; level++)
  {
    const double b_scale = (level + 1 < levels) ? 0.5 * a_scale : 0.0;
    const auto b1 = along(p1, b_scale), a1 = along(p1, a_scale);
    const auto b2 = along(p2, b_scale), a2 = along(p2, a_scale);
    if (b_scale > 0.0)
    {
      integrate_triangle(b1, a1, a2);
      integrate_triangle(b1, a2, b2);
    }
    else
    {
      // Innermost shell degenerates to the small triangle at the corner itself.
      integrate_triangle(c, a1, a2);
    }
    a_scale = b_scale;
  }
}

// Solve the manufactured problem and report every certification quantity. With
// enriched = false this is the pure standard H1 space, which is what certifies the harness
// itself; with enriched = true the same harness additionally carries Palace's singular
// node-gradient enrichment at the reentrant corner, so the two are directly comparable on
// an identical mesh.
// Recovered-flux error indicator, the same principle as Palace's GradFluxErrorEstimator:
//   eta_K = || grad(u_h) - D ||_K,
// where D is a smoother reconstruction of grad(u_h) obtained by global L2 projection into a
// flux space with continuous normal component (RT). Palace projects into RT of the matching
// order; this reproduces that construction locally so the indicator can be compared against
// the exact element error on the same mesh.
//
// `sample_gradient` supplies grad(u_h) at a quadrature point. Passing a version that omits
// the enrichment reproduces exactly what production does today, where the estimator
// receives only the standard block of the field vector.
// `sample_added` supplies a field ADDED to the recovered flux before the error is measured.
// It exists to express "recover the standard part, then re-add the singular flux at its
// existing coefficient", which is the obvious but ineffective augmentation.
template <typename SampleRecovered, typename SampleMeasured, typename SampleAdded>
std::vector<double> MmsRecoveredFluxIndicator(mfem::Mesh &mesh, int order,
                                              int quadrature_order,
                                              MmsCornerQuadrature scheme,
                                              const SampleRecovered &sample_recovered,
                                              const SampleMeasured &sample_measured,
                                              const SampleAdded &sample_added)
{
  mfem::RT_FECollection flux_collection(std::max(order - 1, 0), 2);
  mfem::FiniteElementSpace flux_space(&mesh, &flux_collection);

  // Global L2 projection of grad(u_h) into the RT flux space: (D, w) = (grad(u_h), w).
  mfem::BilinearForm mass(&flux_space);
  mass.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  mass.Assemble();
  mass.Finalize();

  mfem::Vector load(flux_space.GetVSize());
  load = 0.0;
  mfem::Array<int> flux_dofs;
  mfem::DenseMatrix shape;
  mfem::Vector point(3);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &fe = *flux_space.GetFE(element);
    auto &T = *mesh.GetElementTransformation(element);
    mfem::DofTransformation dof_transformation;
    flux_space.GetElementVDofs(element, flux_dofs, dof_transformation);
    mfem::Vector local(fe.GetDof());
    local = 0.0;
    shape.SetSize(fe.GetDof(), 2);
    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, scheme,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          const double weight = weight_ref * T.Weight();
          fe.CalcPhysVShape(T, shape);
          const auto gradient = sample_recovered(element, T, ip, lambda);
          for (int i = 0; i < fe.GetDof(); i++)
          {
            local(i) += weight * (shape(i, 0) * gradient[0] + shape(i, 1) * gradient[1]);
          }
        });
    // A linear-form contribution is a DUAL vector, so it takes the dual DOF transformation,
    // not the primal one. On the 2D RT collections used here the transformation is the
    // identity, so this does not change the numbers, but omitting it would silently break
    // on any collection with a nontrivial transformation.
    dof_transformation.TransformDual(local);
    load.AddElementVector(flux_dofs, local);
  }

  mfem::Vector flux(flux_space.GetVSize());
  flux = 0.0;
  {
    mfem::GSSmoother preconditioner(mass.SpMat());
    mfem::CGSolver cg;
    cg.SetOperator(mass.SpMat());
    cg.SetPreconditioner(preconditioner);
    cg.SetRelTol(1.0e-14);
    cg.SetAbsTol(1.0e-30);
    cg.SetMaxIter(50000);
    cg.SetPrintLevel(-1);
    cg.Mult(load, flux);
    // The recovery is part of the measurement chain, so its solve must be certified too:
    // an unconverged projection would masquerade as indicator error.
    REQUIRE(cg.GetConverged());
    mfem::Vector residual(load.Size());
    mass.SpMat().Mult(flux, residual);
    residual -= load;
    const double norm = load.Norml2();
    const double relative = (norm > 0.0) ? residual.Norml2() / norm : residual.Norml2();
    CAPTURE(cg.GetNumIterations(), relative);
    CHECK(relative < 1.0e-10);
  }
  mfem::GridFunction recovered(&flux_space, flux.GetData());

  // eta_K^2 = int_K |grad(u_h) - D|^2.
  std::vector<double> indicator(mesh.GetNE(), 0.0);
  mfem::Vector recovered_value(2);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double local = 0.0;
    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, scheme,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          const double weight = weight_ref * T.Weight();
          recovered.GetVectorValue(T, ip, recovered_value);
          const auto gradient = sample_measured(element, T, ip, lambda);
          const auto added = sample_added(element, T, ip, lambda);
          const double dx = gradient[0] - (recovered_value(0) + added[0]);
          const double dy = gradient[1] - (recovered_value(1) + added[1]);
          local += weight * (dx * dx + dy * dy);
        });
    indicator[element] = local;
  }
  return indicator;
}

// FEATURE-LEVEL singular flux mode for the recovery space, constructed from FEATURE
// GEOMETRY AND MATERIAL DATA ONLY. The exponent nu is whatever the singular-feature
// extractor computed from the wedge angle and the material fan; the angular factor sin(nu
// theta) is the leading Dirichlet wedge eigenfunction for that geometry. The manufactured
// exact solution is never consulted.
//
//   s_mode = grad(r^nu sin(nu theta))
//
// Properties that make it admissible in an H(div) recovery space:
//   * real-analytic away from the corner, hence normal-trace continuous across every
//   interior
//     face (verified numerically below, not merely asserted);
//   * divergence free, since r^nu sin(nu theta) is harmonic for this nu;
//   * in L2: |s_mode| ~ r^(nu-1) = r^(-1/3), so |s_mode|^2 ~ r^(-2/3), integrable in 2D.
//
// There is ONE mode per singular feature. Deliberately NOT one per local enrichment
// gradient: an unrestricted set of element-local copies could absorb ordinary polynomial
// approximation error and would make the recovery look good for the wrong reason.
std::array<double, 2> MmsSingularFluxMode(double x, double y, double nu)
{
  const double r = std::hypot(x, y);
  if (r == 0.0)
  {
    return {0.0, 0.0};
  }
  const double t = MmsTheta(x, y);
  const double c = nu * std::pow(r, nu - 1.0);
  return {c * std::sin((nu - 1.0) * t), c * std::cos((nu - 1.0) * t)};
}

// Free-coefficient RT (+) S_flux flux recovery, and the residual-alignment fraction rho.
//
// Plain recovery solves  min over q in RT  of ||E - q||.  Augmenting with one free singular
// mode s solves  min over (q, alpha)  of ||E - q - alpha s||, whose optimum removes exactly
//
//   |<r, s_perp>|^2 / ||s_perp||^2,        r = E - P_RT E,   s_perp = s - P_RT s
//
// from the squared recovery error. The USEFUL fraction is therefore
//
//   rho = |<r, s_perp>|^2 / (||r||^2 ||s_perp||^2),
//
// which is what this returns. Note ||s_perp||^2 / ||s||^2 -- how much of the MODE lies
// outside RT -- is a different and much less relevant quantity: a mode holding a tiny share
// of the total field energy can still be strongly aligned with the residual and so remove
// most of it.
//
// The singular amplitude alpha is fitted here, NOT copied from the solution's enrichment
// coefficient; that distinction is what makes this differ from the naive add-back, which is
// provably identical to plain sliced recovery.
struct MmsAugmentedRecovery
{
  std::vector<double> indicator;
  double rho = 0.0;
  double residual_before = 0.0;
  double residual_after = 0.0;
};

template <typename SampleField>
MmsAugmentedRecovery MmsAugmentedFluxIndicator(mfem::Mesh &mesh, int order,
                                               int quadrature_order,
                                               MmsCornerQuadrature scheme, double nu,
                                               const SampleField &sample_field)
{
  mfem::RT_FECollection flux_collection(std::max(order - 1, 0), 2);
  mfem::FiniteElementSpace flux_space(&mesh, &flux_collection);
  const int flux_size = flux_space.GetVSize();

  mfem::BilinearForm mass(&flux_space);
  mass.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  mass.Assemble();
  mass.Finalize();

  // Assemble the RT load for the field, the RT-mode cross term, and the mode's own norms.
  mfem::Vector field_load(flux_size), cross(flux_size);
  field_load = 0.0;
  cross = 0.0;
  double mode_mass = 0.0, mode_field = 0.0;
  {
    mfem::Array<int> flux_dofs;
    mfem::DenseMatrix shape;
    mfem::Vector point(3);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &fe = *flux_space.GetFE(element);
      auto &T = *mesh.GetElementTransformation(element);
      mfem::DofTransformation dof_transformation;
      flux_space.GetElementVDofs(element, flux_dofs, dof_transformation);
      mfem::Vector local_field(fe.GetDof()), local_cross(fe.GetDof());
      local_field = 0.0;
      local_cross = 0.0;
      shape.SetSize(fe.GetDof(), 2);
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, scheme,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            T.Transform(ip, point);
            const double weight = weight_ref * T.Weight();
            fe.CalcPhysVShape(T, shape);
            const auto field = sample_field(element, T, ip, lambda);
            const auto mode = MmsSingularFluxMode(point(0), point(1), nu);
            mode_mass += weight * (mode[0] * mode[0] + mode[1] * mode[1]);
            mode_field += weight * (mode[0] * field[0] + mode[1] * field[1]);
            for (int i = 0; i < fe.GetDof(); i++)
            {
              local_field(i) += weight * (shape(i, 0) * field[0] + shape(i, 1) * field[1]);
              local_cross(i) += weight * (shape(i, 0) * mode[0] + shape(i, 1) * mode[1]);
            }
          });
      dof_transformation.TransformDual(local_field);
      dof_transformation.TransformDual(local_cross);
      field_load.AddElementVector(flux_dofs, local_field);
      cross.AddElementVector(flux_dofs, local_cross);
    }
  }

  const auto solve_mass = [&](const mfem::Vector &rhs)
  {
    mfem::Vector x(rhs.Size());
    x = 0.0;
    mfem::GSSmoother preconditioner(mass.SpMat());
    mfem::CGSolver cg;
    cg.SetOperator(mass.SpMat());
    cg.SetPreconditioner(preconditioner);
    cg.SetRelTol(1.0e-14);
    cg.SetAbsTol(1.0e-30);
    cg.SetMaxIter(50000);
    cg.SetPrintLevel(-1);
    cg.Mult(rhs, x);
    REQUIRE(cg.GetConverged());
    return x;
  };
  const auto plain_rt = solve_mass(field_load);
  const auto mass_inverse_cross = solve_mass(cross);

  // s_perp energy and the residual-mode inner product, both by algebra rather than extra
  // quadrature: <r, s_perp> = <E, s> - <P_RT E, s> and ||s_perp||^2 = M_SS - c^T M^-1 c.
  const double s_perp_energy = mode_mass - mfem::InnerProduct(cross, mass_inverse_cross);
  const double residual_mode = mode_field - mfem::InnerProduct(cross, plain_rt);
  const double alpha = (s_perp_energy > 0.0) ? residual_mode / s_perp_energy : 0.0;

  // Free-coefficient optimum: q = P_RT(E - alpha s).
  mfem::Vector augmented_rt(plain_rt);
  augmented_rt.Add(-alpha, mass_inverse_cross);

  MmsAugmentedRecovery result;
  result.indicator.assign(mesh.GetNE(), 0.0);
  mfem::GridFunction plain(&flux_space, const_cast<double *>(plain_rt.GetData()));
  mfem::GridFunction augmented(&flux_space, augmented_rt.GetData());
  mfem::Vector plain_value(2), augmented_value(2), point(3);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double local_before = 0.0, local_after = 0.0;
    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, scheme,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          T.Transform(ip, point);
          const double weight = weight_ref * T.Weight();
          const auto field = sample_field(element, T, ip, lambda);
          const auto mode = MmsSingularFluxMode(point(0), point(1), nu);
          plain.GetVectorValue(T, ip, plain_value);
          augmented.GetVectorValue(T, ip, augmented_value);
          const double bx = field[0] - plain_value(0), by = field[1] - plain_value(1);
          const double ax = field[0] - (augmented_value(0) + alpha * mode[0]);
          const double ay = field[1] - (augmented_value(1) + alpha * mode[1]);
          local_before += weight * (bx * bx + by * by);
          local_after += weight * (ax * ax + ay * ay);
        });
    result.indicator[element] = local_after;
    result.residual_before += local_before;
    result.residual_after += local_after;
  }
  result.rho =
      (result.residual_before > 0.0 && s_perp_energy > 0.0)
          ? (residual_mode * residual_mode) / (result.residual_before * s_perp_energy)
          : 0.0;
  return result;
}

// Test-local operation counters make candidate independence observable. The sparse local
// builder snapshots them; any future call into this global assembly/solve path will fail
// its zero-delta regression instead of being hidden behind a hard-coded report value.
long long MmsGlobalFineMatrixAssemblies = 0;
long long MmsGlobalFineSolves = 0;

// `frozen_enrichment`, when non-null, holds the singular coefficients FIXED at the supplied
// values instead of solving for them. Used by the hierarchical ablations to suppress the
// singular part of the coarse response dc = -A_cc^{-1} A_cn dn.
MmsSolveReport MmsSolveOnMesh(mfem::Mesh &mesh, int order, int quadrature_order = 20,
                              bool enriched = false,
                              MmsCornerQuadrature scheme = MmsCornerQuadrature::GRADED,
                              double override_nu = 0.0,
                              const mfem::Vector *frozen_enrichment = nullptr,
                              bool solve_system = true)
{
  mfem::H1_FECollection collection(order, 2);
  mfem::FiniteElementSpace space(&mesh, &collection);
  const int standard_size = space.GetVSize();
  REQUIRE(standard_size == space.GetTrueVSize());

  // Singular enrichment at the reentrant corner. The two selected faces bound a
  // homogeneous 3 pi / 2 dielectric wedge, so the extracted exponent must be exactly the
  // nu = 2/3 the manufactured solution was built from.
  fem::singular::TriangleFeatureTopology features;
  fem::singular::TriangleDofTopology topology;
  fem::singular::LocalSparseH1EnrichmentMatrices enrichment;
  int enrichment_size = 0;
  if (enriched)
  {
    features =
        fem::singular::ExtractSerialLineFeatures(mesh, {MmsReentrantAttribute}, {{1, 1.0}});
    REQUIRE(features.vertices.size() == 1);
    REQUIRE(features.vertices[0].type == fem::singular::FeatureVertexType::CORNER);
    // The extractor must recover exactly the exponent the manufactured solution was built
    // from, with no help from the harness.
    REQUIRE_THAT(features.vertices[0].nu, WithinRel(MmsNu, 1.0e-12));
    if (override_nu > 0.0)
    {
      features.vertices[0].nu = override_nu;
    }
    topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
    enrichment_size = static_cast<int>(topology.h1_dofs.size());
    REQUIRE(enrichment_size > 0);

    const std::vector<fem::singular::IsotropicMaterialCoefficients> materials(mesh.GetNE(),
                                                                              {1.0, 1.0});
    const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-10, 1.0e-10, 8};
    enrichment = fem::singular::AssembleLocalSparseH1EnrichmentMatrices(topology, space,
                                                                        materials, options);
  }
  const int combined_size = standard_size + enrichment_size;

  // Combined stiffness. The standard block comes from MFEM; the coupling and
  // enrichment-enrichment blocks come from Palace's singular assembly.
  mfem::SparseMatrix combined(combined_size, combined_size);
  {
    mfem::BilinearForm standard(&space);
    standard.AddDomainIntegrator(new mfem::DiffusionIntegrator());
    standard.Assemble();
    standard.Finalize();
    const auto &block = standard.SpMat();
    for (int row = 0; row < block.Height(); row++)
    {
      for (int entry = block.GetI()[row]; entry < block.GetI()[row + 1]; entry++)
      {
        combined.Add(row, block.GetJ()[entry], block.GetData()[entry]);
      }
    }
  }
  if (enriched)
  {
    const auto add_block =
        [&combined](const mfem::SparseMatrix &block, int row_offset, int column_offset)
    {
      for (int row = 0; row < block.Height(); row++)
      {
        for (int entry = block.GetI()[row]; entry < block.GetI()[row + 1]; entry++)
        {
          combined.Add(row_offset + row, column_offset + block.GetJ()[entry],
                       block.GetData()[entry]);
        }
      }
    };
    add_block(*enrichment.diffusion.standard_enrichment, 0, standard_size);
    add_block(*enrichment.diffusion.enrichment_standard, standard_size, 0);
    add_block(*enrichment.diffusion.enrichment_enrichment, standard_size, standard_size);
  }
  combined.Finalize();
  MmsGlobalFineMatrixAssemblies++;

  // Right-hand side. The standard part is an ordinary MFEM linear form at high order; the
  // enrichment part is integrated against the singular basis with the same corner-resolving
  // rule used everywhere else in the harness.
  mfem::Vector rhs(combined_size);
  rhs = 0.0;
  {
    MmsSourceCoefficient source;
    auto *integrator = new mfem::DomainLFIntegrator(source);
    integrator->SetIntRule(&mfem::IntRules.Get(mfem::Geometry::TRIANGLE, quadrature_order));
    mfem::LinearForm standard(&space);
    standard.AddDomainIntegrator(integrator);
    standard.Assemble();
    const mfem::Vector &assembled = standard;
    for (int i = 0; i < standard_size; i++)
    {
      rhs(i) = assembled(i);
    }
  }
  if (enriched)
  {
    mfem::Vector point(3);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &element_dofs = topology.elements[element];
      if (element_dofs.h1.empty())
      {
        continue;
      }
      auto &T = *mesh.GetElementTransformation(element);
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, scheme,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            T.Transform(ip, point);
            const double weight = weight_ref * T.Weight();
            const double f = MmsF(point(0), point(1));
            for (const auto &dof : element_dofs.h1)
            {
              rhs(standard_size + static_cast<int>(dof.dof)) +=
                  weight * f *
                  fem::singular::EvaluateTriangleNodeGradientPotential(
                      lambda, dof.basis.nodes[0], dof.basis.nodes[1], dof.basis.nu);
            }
          });
    }
  }

  // Essential DOFs. u vanishes on all six faces, so every standard boundary DOF is
  // constrained, together with every enrichment function whose trace lies on a boundary
  // segment. The enrichment classification uses the production API rather than an ad hoc
  // rule, since that classification is itself part of what is being certified.
  std::vector<bool> essential(combined_size, false);
  {
    mfem::Array<int> standard_essential;
    space.GetBoundaryTrueDofs(standard_essential);
    for (int dof : standard_essential)
    {
      essential[dof] = true;
    }
  }
  if (enriched)
  {
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), topology);
    // Serial: the canonical local ordering is already the true-DOF ordering, which is what
    // makes the returned indices usable as offsets into the combined vector.
    REQUIRE(numbering.h1.owned_offset == 0);
    REQUIRE(numbering.h1.owned_size == enrichment_size);
    for (int local = 0; local < enrichment_size; local++)
    {
      REQUIRE(numbering.h1.local_to_true[local] == local);
    }
    std::set<std::array<fem::singular::GlobalVertexId, 2>> boundary_segments;
    mfem::Array<int> segment_vertices;
    for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
    {
      mesh.GetBdrElementVertices(boundary, segment_vertices);
      REQUIRE(segment_vertices.Size() == 2);
      std::array<fem::singular::GlobalVertexId, 2> segment{segment_vertices[0],
                                                           segment_vertices[1]};
      std::sort(segment.begin(), segment.end());
      boundary_segments.insert(segment);
    }
    const auto enrichment_essential =
        fem::singular::GetEssentialTriangleH1TrueDofsOnSegments(
            Mpi::World(),
            std::vector<std::array<fem::singular::GlobalVertexId, 2>>(
                boundary_segments.begin(), boundary_segments.end()),
            topology, numbering);
    for (int dof : enrichment_essential)
    {
      essential[standard_size + dof] = true;
    }
  }

  // Freezing the enrichment: treat every singular DOF as essential at the supplied value.
  // Unlike the homogeneous Dirichlet data this is a NONZERO prescription, so the reduction
  // below must carry a load correction rather than being a plain submatrix extraction.
  mfem::Vector prescribed(combined_size);
  prescribed = 0.0;
  if (enriched && frozen_enrichment)
  {
    REQUIRE(frozen_enrichment->Size() >= enrichment_size);
    for (int i = 0; i < enrichment_size; i++)
    {
      essential[standard_size + i] = true;
      prescribed(standard_size + i) = (*frozen_enrichment)(i);
    }
  }

  // Candidate construction needs only the fine test-space operator, load, and constraints.
  // Returning here makes the absence of a global p+1 solve measurable rather than merely a
  // data-flow observation inside the global benchmark.
  if (!solve_system)
  {
    MmsSolveReport report;
    report.elements = mesh.GetNE();
    report.dofs = combined_size;
    report.enrichment_dofs = enrichment_size;
    report.combined_matrix = std::make_shared<mfem::SparseMatrix>(combined);
    report.combined_rhs = rhs;
    report.essential_mask = essential;
    return report;
  }

  // Reduce to the free DOFs. Essential values are zero except for a frozen enrichment,
  // whose nonzero prescription is moved to the right-hand side.
  std::vector<int> free_dofs, combined_to_free(combined_size, -1);
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (!essential[dof])
    {
      combined_to_free[dof] = static_cast<int>(free_dofs.size());
      free_dofs.push_back(dof);
    }
  }
  const int free_size = static_cast<int>(free_dofs.size());
  mfem::SparseMatrix system(free_size, free_size);
  mfem::Vector reduced_rhs(free_size);
  for (int row = 0; row < free_size; row++)
  {
    const int combined_row = free_dofs[row];
    reduced_rhs(row) = rhs(combined_row);
    for (int entry = combined.GetI()[combined_row];
         entry < combined.GetI()[combined_row + 1]; entry++)
    {
      const int combined_column = combined.GetJ()[entry];
      const int column = combined_to_free[combined_column];
      if (column >= 0)
      {
        system.Add(row, column, combined.GetData()[entry]);
      }
      else
      {
        // Eliminate the prescribed column into the load.
        reduced_rhs(row) -= combined.GetData()[entry] * prescribed(combined_column);
      }
    }
  }
  system.Finalize();

  // Tight iterative solve. A direct solver is not available in this build, so drive CG to
  // near machine precision instead; the algebraic-residual check below verifies that the
  // solver contributes far less than the discretization error.
  mfem::Vector reduced_solution(free_size);
  reduced_solution = 0.0;
  MmsGlobalFineSolves++;
  {
    mfem::GSSmoother preconditioner(system);
    mfem::CGSolver cg;
    cg.SetOperator(system);
    cg.SetPreconditioner(preconditioner);
    cg.SetRelTol(1.0e-14);
    cg.SetAbsTol(1.0e-30);
    cg.SetMaxIter(50000);
    cg.SetPrintLevel(-1);
    cg.Mult(reduced_rhs, reduced_solution);
  }
  mfem::Vector combined_solution(prescribed);
  for (int row = 0; row < free_size; row++)
  {
    combined_solution(free_dofs[row]) = reduced_solution(row);
  }

  MmsSolveReport report;
  report.elements = mesh.GetNE();
  report.dofs = combined_size;
  report.enrichment_dofs = enrichment_size;

  // Stage 2(i): algebraic residual of the reduced system.
  {
    mfem::Vector residual(free_size);
    system.Mult(reduced_solution, residual);
    residual -= reduced_rhs;
    const double norm = reduced_rhs.Norml2();
    report.algebraic_residual = (norm > 0.0) ? residual.Norml2() / norm : residual.Norml2();
  }

  // The standard part of the discrete solution as a grid function, so MFEM's own evaluator
  // supplies its value and gradient.
  mfem::GridFunction solution(&space);
  for (int i = 0; i < standard_size; i++)
  {
    solution(i) = combined_solution(i);
  }
  mfem::Vector enrichment_coefficients(std::max(enrichment_size, 1));
  enrichment_coefficients = 0.0;
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_coefficients(i) = combined_solution(standard_size + i);
  }

  // Element-wise quantities by quadrature against the analytic solution. grad(u) carries
  // r^(-1/3) at the reentrant corner, so a plain rule cannot resolve it on the incident
  // elements no matter how high its order: the integrand is not polynomial. Those elements
  // therefore use a corner-resolving rule.
  report.element_error.assign(mesh.GetNE(), 0.0);
  mfem::Vector discrete_gradient(2), point(3);

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double local_error = 0.0, local_a = 0.0, local_l = 0.0;

    // Physical barycentric gradients, needed to evaluate the singular basis. These are
    // constant on an affine triangle.
    fem::singular::TriangleBarycentricGradients grad_lambda{};
    const bool element_enriched = enriched && !topology.elements[element].h1.empty();
    if (element_enriched)
    {
      double jacobian_determinant;
      grad_lambda =
          fem::singular::GetAffineTriangleBarycentricGradients(T, jacobian_determinant);
    }

    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, scheme,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          T.Transform(ip, point);
          const double weight = weight_ref * T.Weight();
          const auto exact = MmsGradU(point(0), point(1));
          solution.GetGradient(T, discrete_gradient);
          double discrete_value = solution.GetValue(T, ip);
          double gx = discrete_gradient(0), gy = discrete_gradient(1);
          if (element_enriched)
          {
            const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
                topology.elements[element], enrichment_coefficients, lambda, grad_lambda);
            discrete_value += singular.potential;
            gx += singular.gradient[0];
            gy += singular.gradient[1];
          }
          const double dx = exact[0] - gx, dy = exact[1] - gy;
          local_error += weight * (dx * dx + dy * dy);
          local_a += weight * (exact[0] * gx + exact[1] * gy);
          local_l += weight * MmsF(point(0), point(1)) * discrete_value;
        });

    report.element_error[element] = local_error;
    report.error_direct += local_error;
    if (MmsCornerVertex(mesh, element) >= 0)
    {
      report.corner_error += local_error;
    }
    report.a_u_uh += local_a;
    report.l_uh += local_l;
  }

  // a(u_h,u_h) from the assembled operator, i.e. algebraically rather than by quadrature.
  {
    mfem::Vector product(combined_size);
    combined.Mult(combined_solution, product);
    report.a_uh_uh = mfem::InnerProduct(combined_solution, product);
  }

  report.error_expanded = MmsReferenceEnergy - 2.0 * report.a_u_uh + report.a_uh_uh;
  report.error_galerkin = MmsReferenceEnergy - report.a_uh_uh;
  report.standard_order = order;
  report.combined_matrix = std::make_shared<mfem::SparseMatrix>(combined);
  report.combined_rhs = rhs;
  report.combined_solution_vector = combined_solution;
  report.essential_mask = essential;
  report.standard_coefficients.SetSize(standard_size);
  for (int i = 0; i < standard_size; i++)
  {
    report.standard_coefficients(i) = solution(i);
  }
  report.enrichment_block.SetSize(std::max(enrichment_size, 1));
  report.enrichment_block = 0.0;
  for (int i = 0; i < enrichment_size; i++)
  {
    report.enrichment_block(i) = enrichment_coefficients(i);
  }

  // Stage 6 indicators. The two variants differ only in whether the enrichment contributes
  // to the field the estimator sees, which is precisely the difference between the enriched
  // discrete solution and what production hands to GradFluxErrorEstimator today.
  const auto sample = [&](bool include_enrichment)
  {
    return [&, include_enrichment](int element, mfem::ElementTransformation &T,
                                   const mfem::IntegrationPoint &ip,
                                   const fem::singular::TriangleBarycentricPoint &lambda)
    {
      mfem::Vector standard(2);
      solution.GetGradient(T, standard);
      std::array<double, 2> gradient{standard(0), standard(1)};
      if (include_enrichment && enriched && !topology.elements[element].h1.empty())
      {
        double jacobian_determinant;
        const auto grad_lambda =
            fem::singular::GetAffineTriangleBarycentricGradients(T, jacobian_determinant);
        // GetAffineTriangleBarycentricGradients moves the integration point, so restore it.
        T.SetIntPoint(&ip);
        const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
            topology.elements[element], enrichment_coefficients, lambda, grad_lambda);
        gradient[0] += singular.gradient[0];
        gradient[1] += singular.gradient[1];
      }
      return gradient;
    };
  };
  // The singular flux on its own, i.e. the enrichment's contribution to grad(u_h).
  const auto singular_only = [&](int element, mfem::ElementTransformation &T,
                                 const mfem::IntegrationPoint &ip,
                                 const fem::singular::TriangleBarycentricPoint &lambda)
  {
    std::array<double, 2> value{0.0, 0.0};
    if (enriched && !topology.elements[element].h1.empty())
    {
      double jacobian_determinant;
      const auto grad_lambda =
          fem::singular::GetAffineTriangleBarycentricGradients(T, jacobian_determinant);
      T.SetIntPoint(&ip);
      const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
          topology.elements[element], enrichment_coefficients, lambda, grad_lambda);
      value = {singular.gradient[0], singular.gradient[1]};
    }
    return value;
  };
  const auto zero = [](int, mfem::ElementTransformation &, const mfem::IntegrationPoint &,
                       const fem::singular::TriangleBarycentricPoint &)
  { return std::array<double, 2>{0.0, 0.0}; };

  // Recover the full field, measure against the full field.
  report.indicator_full = MmsRecoveredFluxIndicator(mesh, order, quadrature_order, scheme,
                                                    sample(true), sample(true), zero);
  // Production: recover the standard block, measure against the standard block.
  report.indicator_standard = MmsRecoveredFluxIndicator(
      mesh, order, quadrature_order, scheme, sample(false), sample(false), zero);
  // Naive "augmentation": recover only the STANDARD field, then re-add the singular flux at
  // its EXISTING finite-element coefficient before measuring against the full field.
  // Algebraically
  //   (E_std + E_sing) - (R(E_std) + E_sing) = E_std - R(E_std),
  // so this must come out IDENTICAL to indicator_standard. It is computed only to prove
  // that, because it is the obvious candidate improvement and it is a dead end: a genuine
  // augmented recovery has to project the full field into a global RT + S_flux space with
  // INDEPENDENTLY recovered singular coefficients, not reuse the solution's own.
  report.indicator_naive_augmented = MmsRecoveredFluxIndicator(
      mesh, order, quadrature_order, scheme, sample(false), sample(true), singular_only);

  // Free-coefficient RT (+) S_flux recovery on the PRODUCTION field, i.e. the sliced
  // standard block. That is the correct comparison: it changes only the recovery space,
  // holding the input field fixed at what production actually supplies.
  if (enriched)
  {
    const auto augmented = MmsAugmentedFluxIndicator(
        mesh, order, quadrature_order, scheme, features.vertices[0].nu, sample(false));
    report.indicator_augmented = augmented.indicator;
    report.augmentation_rho = augmented.rho;
    report.residual_energy_before = augmented.residual_before;
    report.residual_energy_after = augmented.residual_after;
  }
  return report;
}

// Structured-mesh convenience wrapper.
MmsSolveReport MmsSolve(int n, int order, int quadrature_order = 20, bool enriched = false,
                        MmsCornerQuadrature scheme = MmsCornerQuadrature::GRADED,
                        double override_nu = 0.0)
{
  auto mesh = MmsLShapeMesh(n);
  return MmsSolveOnMesh(mesh, order, quadrature_order, enriched, scheme, override_nu);
}

MmsSolveReport MmsSolveStandard(int n, int order, int quadrature_order = 20)
{
  return MmsSolve(n, order, quadrature_order, false, MmsCornerQuadrature::GRADED);
}

// Solve with the enrichment exponent deliberately overridden to a wrong value. Used as a
// negative control: the enriched space still contains the standard one, so it can never be
// worse, but a wrong exponent should recover far less of the corner error than nu = 2/3.
MmsSolveReport MmsSolveWrongExponent(int n, int order, double nu)
{
  return MmsSolve(n, order, 20, true, MmsCornerQuadrature::GRADED, nu);
}

}  // namespace

TEST_CASE("Manufactured L-shape standard-space harness satisfies its algebra checks",
          "[singularmms][singularelements][Serial]")
{
  // Stages 1-3 on the STANDARD space. This certifies the harness itself: if the algebra
  // identities do not hold without enrichment, no conclusion about the enriched space
  // would be trustworthy. Enrichment is added only after this passes.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4, 8);
  const int order = GENERATE(1, 2);
  CAPTURE(n, order);

  // Quadrature order is swept so the a(u,u_h) vs l(u_h) gap can be attributed: if it
  // shrinks with the rule, it is quadrature error on the r^(2/3) integrand rather than an
  // assembly inconsistency.
  const int quadrature_order = GENERATE(20, 30);
  CAPTURE(quadrature_order);
  const auto report = MmsSolveStandard(n, order, quadrature_order);
  CAPTURE(report.elements, report.dofs, report.algebraic_residual, report.a_u_uh,
          report.l_uh, report.a_uh_uh, report.error_direct, report.error_expanded,
          report.error_galerkin);

  // Stage 2(i): the direct solve must be essentially exact.
  CHECK(report.algebraic_residual < 1.0e-10);

  // Stage 2(ii): a(u,u_h) == l(u_h). Both are quadratures of different integrands, and
  // they coincide only if the source and the discrete solution are mutually consistent,
  // so this is a genuine assembly check rather than an identity by construction.
  // With corner-graded quadrature this agrees to ~1e-11 relative; the gate is set two
  // orders looser to tolerate coarser meshes and higher orders, and is still six orders
  // below any discretization effect.
  const double scale = std::max({1.0, std::abs(report.a_u_uh), std::abs(report.l_uh)});
  CHECK_THAT(report.a_u_uh - report.l_uh, WithinAbs(0.0, 1.0e-9 * scale));

  // Stage 3: the three error expressions must agree well below the error itself.
  REQUIRE(report.error_direct > 0.0);
  // Achieved agreement is ~3e-8 relative, so a 1e-5 relative gate is comfortably met
  // while remaining far below the discretization error being measured.
  const double tolerance = 1.0e-5 * report.error_direct;
  CHECK_THAT(report.error_expanded, WithinAbs(report.error_direct, tolerance));
  CHECK_THAT(report.error_galerkin, WithinAbs(report.error_direct, tolerance));

  // The element errors must sum to the direct total.
  const double summed =
      std::accumulate(report.element_error.begin(), report.element_error.end(), 0.0);
  CHECK_THAT(summed, WithinRel(report.error_direct, 1.0e-12));
}

TEST_CASE("Manufactured L-shape enriched-space harness satisfies its algebra checks",
          "[singularmms][singularelements][Serial]")
{
  // Stage 5a. The SAME stage-2 and stage-3 gates, now with Palace's singular enrichment
  // active. Nothing is compared against the standard space until these pass: an enriched
  // error that failed the algebra checks would be meaningless.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4, 8);
  const int order = GENERATE(1, 2);
  // Both corner rules are swept. They resolve the algebraic corner behaviour by entirely
  // different means, so agreement is evidence about the quadrature and not just about one
  // implementation of it.
  const auto scheme = GENERATE(MmsCornerQuadrature::GRADED, MmsCornerQuadrature::DUFFY);
  CAPTURE(n, order, scheme == MmsCornerQuadrature::DUFFY);

  const auto report = MmsSolve(n, order, 20, true, scheme);
  CAPTURE(report.elements, report.dofs, report.enrichment_dofs, report.algebraic_residual,
          report.a_u_uh, report.l_uh, report.a_uh_uh, report.error_direct,
          report.error_expanded, report.error_galerkin);

  // The enrichment must actually be present, or this test would silently degenerate into
  // the standard-space one.
  CHECK(report.enrichment_dofs > 0);

  // Stage 2(i).
  CHECK(report.algebraic_residual < 1.0e-10);

  // Stage 2(ii): a(u,u_h) == l(u_h) for the enriched u_h. This is the assembly check that
  // matters most here, since it couples Palace's singular element matrices to an
  // independently integrated singular right-hand side. An error in either would break it.
  // Both corner rules achieve ~1e-11 or better, so the gate sits two orders looser.
  const double scale = std::max({1.0, std::abs(report.a_u_uh), std::abs(report.l_uh)});
  CHECK_THAT(report.a_u_uh - report.l_uh, WithinAbs(0.0, 1.0e-9 * scale));

  // Stage 3.
  REQUIRE(report.error_direct > 0.0);
  const double tolerance = 1.0e-5 * report.error_direct;
  CHECK_THAT(report.error_expanded, WithinAbs(report.error_direct, tolerance));
  CHECK_THAT(report.error_galerkin, WithinAbs(report.error_direct, tolerance));

  const double summed =
      std::accumulate(report.element_error.begin(), report.element_error.end(), 0.0);
  CHECK_THAT(summed, WithinRel(report.error_direct, 1.0e-12));
}

namespace
{

// Stage 6 ranking metrics for one indicator against the exact element errors.
struct MmsRankingMetrics
{
  double effectivity = 0.0;           // eta / e, both energy-norm (i.e. sqrt of the sums)
  double correlation = 0.0;           // Spearman rank correlation of eta_K^2 against e_K^2
  double dorfler_true_capture = 0.0;  // fraction of true error in the marked set
  double marked_overlap = 0.0;        // Jaccard overlap with the true-error marked set
  std::size_t marked = 0;
  std::size_t true_marked = 0;
  bool corner_marked = false;  // is any corner-incident element marked?

  // Global Spearman is diluted on this problem: the error concentrates so sharply at the
  // corner that most elements carry almost none of it, and their relative ordering is
  // numerical noise which the rank correlation nonetheless weights equally. The decile
  // figure restricts attention to the elements that actually carry error, and top_decile
  // _overlap is the quantity AMR truly depends on -- whether the indicator's top decile is
  // the true top decile.
  double top_decile_overlap = 0.0;
  double top_decile_true_capture = 0.0;
};

// Dorfler marking: smallest set whose indicator sum reaches theta of the total.
std::vector<std::size_t> MmsDorflerMark(const std::vector<double> &value, double theta)
{
  std::vector<std::size_t> order(value.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(),
            [&value](std::size_t a, std::size_t b) { return value[a] > value[b]; });
  const double total = std::accumulate(value.begin(), value.end(), 0.0);
  std::vector<std::size_t> marked;
  double running = 0.0;
  for (std::size_t element : order)
  {
    if (running >= theta * total)
    {
      break;
    }
    running += value[element];
    marked.push_back(element);
  }
  return marked;
}

struct MmsOracleComparison
{
  double rank = 0.0;
  double extend = 0.0;
  double concentration = 0.0;
  double true_capture = 0.0;
  double correction_capture = 0.0;
  std::size_t marked = 0;
  std::size_t oracle_marked = 0;
};

// Compare one candidate ordering to exact-error Dörfler marking at identical theta and on
// the identical mesh. RANK holds the candidate's own cardinality fixed; EXTEND holds the
// oracle cardinality fixed. correction_capture is the share of an independently supplied
// global hierarchical correction lying in the candidate's marked set.
MmsOracleComparison MmsCompareToOracle(const std::vector<double> &indicator,
                                       const std::vector<double> &exact,
                                       const std::vector<double> &global_correction,
                                       double theta)
{
  REQUIRE(indicator.size() == exact.size());
  REQUIRE(global_correction.size() == exact.size());
  const auto marked = MmsDorflerMark(indicator, theta);
  const auto oracle = MmsDorflerMark(exact, theta);
  REQUIRE(!marked.empty());
  REQUIRE(!oracle.empty());
  const auto top_k = [](const std::vector<double> &value, std::size_t count)
  {
    std::vector<std::size_t> index(value.size());
    std::iota(index.begin(), index.end(), std::size_t{0});
    count = std::min(count, value.size());
    std::partial_sort(index.begin(), index.begin() + count, index.end(),
                      [&value](std::size_t a, std::size_t b)
                      { return value[a] > value[b]; });
    return std::vector<std::size_t>(index.begin(), index.begin() + count);
  };
  const auto sum_on =
      [](const std::vector<double> &value, const std::vector<std::size_t> &indices)
  {
    double sum = 0.0;
    for (std::size_t index : indices)
    {
      sum += value[index];
    }
    return sum;
  };
  const double exact_total = std::accumulate(exact.begin(), exact.end(), 0.0);
  const double correction_total =
      std::accumulate(global_correction.begin(), global_correction.end(), 0.0);
  const double marked_capture = sum_on(exact, marked);
  const double oracle_capture = sum_on(exact, oracle);

  MmsOracleComparison comparison;
  comparison.marked = marked.size();
  comparison.oracle_marked = oracle.size();
  comparison.true_capture = marked_capture / exact_total;
  comparison.rank = marked_capture / sum_on(exact, top_k(exact, marked.size()));
  comparison.extend = sum_on(exact, top_k(indicator, oracle.size())) / oracle_capture;
  comparison.concentration =
      static_cast<double>(marked.size()) / static_cast<double>(oracle.size());
  comparison.correction_capture = sum_on(global_correction, marked) / correction_total;
  return comparison;
}

MmsRankingMetrics MmsComputeRanking(const mfem::Mesh &mesh,
                                    const std::vector<double> &indicator,
                                    const std::vector<double> &exact, double theta = 0.3)
{
  MmsRankingMetrics metrics;
  const double indicator_total = std::accumulate(indicator.begin(), indicator.end(), 0.0);
  const double exact_total = std::accumulate(exact.begin(), exact.end(), 0.0);
  // Effectivity compares energy norms, so both sums are square-rooted. Comparing eta^2
  // against e would be a units error.
  metrics.effectivity = std::sqrt(indicator_total) / std::sqrt(exact_total);

  // Spearman rank correlation, which is what actually matters for marking: only the
  // ordering of the indicator affects which elements get refined.
  const auto ranks = [](const std::vector<double> &value)
  {
    std::vector<std::size_t> order(value.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&value](std::size_t a, std::size_t b) { return value[a] < value[b]; });
    std::vector<double> rank(value.size());
    for (std::size_t i = 0; i < order.size(); i++)
    {
      rank[order[i]] = static_cast<double>(i);
    }
    return rank;
  };
  const auto indicator_rank = ranks(indicator), exact_rank = ranks(exact);
  const double mean = 0.5 * static_cast<double>(indicator.size() - 1);
  double covariance = 0.0, indicator_variance = 0.0, exact_variance = 0.0;
  for (std::size_t i = 0; i < indicator.size(); i++)
  {
    const double a = indicator_rank[i] - mean, b = exact_rank[i] - mean;
    covariance += a * b;
    indicator_variance += a * a;
    exact_variance += b * b;
  }
  metrics.correlation = covariance / std::sqrt(indicator_variance * exact_variance);

  // Dorfler sets from the indicator and from the true error.
  const auto marked = MmsDorflerMark(indicator, theta);
  const auto true_marked = MmsDorflerMark(exact, theta);
  metrics.marked = marked.size();
  metrics.true_marked = true_marked.size();

  // How much of the TRUE error the indicator-marked set actually captures. This is the
  // question that matters for AMR: a badly scaled but correctly ordered indicator still
  // refines the right elements.
  double captured = 0.0;
  for (std::size_t element : marked)
  {
    captured += exact[element];
  }
  metrics.dorfler_true_capture = captured / exact_total;

  const std::set<std::size_t> marked_set(marked.begin(), marked.end());
  const std::set<std::size_t> true_set(true_marked.begin(), true_marked.end());
  std::vector<std::size_t> intersection, union_set;
  std::set_intersection(marked_set.begin(), marked_set.end(), true_set.begin(),
                        true_set.end(), std::back_inserter(intersection));
  std::set_union(marked_set.begin(), marked_set.end(), true_set.begin(), true_set.end(),
                 std::back_inserter(union_set));
  metrics.marked_overlap = union_set.empty() ? 1.0
                                             : static_cast<double>(intersection.size()) /
                                                   static_cast<double>(union_set.size());

  for (std::size_t element : marked)
  {
    if (MmsCornerVertex(mesh, static_cast<int>(element)) >= 0)
    {
      metrics.corner_marked = true;
    }
  }

  // Top-decile agreement. This is the AMR-relevant question, and it is not the same as the
  // global rank correlation.
  {
    const auto top = [](const std::vector<double> &value, std::size_t count)
    {
      std::vector<std::size_t> order(value.size());
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::partial_sort(order.begin(), order.begin() + count, order.end(),
                        [&value](std::size_t a, std::size_t b)
                        { return value[a] > value[b]; });
      return std::set<std::size_t>(order.begin(), order.begin() + count);
    };
    const std::size_t count = std::max<std::size_t>(1, (indicator.size() + 9) / 10);
    const auto indicator_top = top(indicator, count), exact_top = top(exact, count);
    std::vector<std::size_t> intersection;
    std::set_intersection(indicator_top.begin(), indicator_top.end(), exact_top.begin(),
                          exact_top.end(), std::back_inserter(intersection));
    metrics.top_decile_overlap =
        static_cast<double>(intersection.size()) / static_cast<double>(count);
    double decile_captured = 0.0;
    for (std::size_t element : indicator_top)
    {
      decile_captured += exact[element];
    }
    metrics.top_decile_true_capture = decile_captured / exact_total;
  }
  return metrics;
}

}  // namespace

TEST_CASE("Manufactured L-shape indicator is certified against the exact element error",
          "[singularmms][singularelements][Serial]")
{
  // Stage 6. Only now, with the static-mesh algebra and error quantities certified, are
  // effectivity and element-ranking metrics meaningful.
  //
  // The indicator is a recovered-flux estimator built the same way as Palace's
  // GradFluxErrorEstimator, and it is evaluated TWICE per solve: once on the full enriched
  // field, and once on the standard block alone, which is what production actually passes
  // (electrostaticsolver.cpp:1592 slices E down to GetNDSpace().GetTrueVSize() before
  // calling AddErrorIndicator). Comparing the two isolates the effect of that slice.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  const int n = GENERATE(4, 8);
  CAPTURE(order, n);

  auto mesh = MmsLShapeMesh(n);
  const auto report = MmsSolve(n, order, 20, true);
  REQUIRE(report.indicator_full.size() == report.element_error.size());
  REQUIRE(report.indicator_standard.size() == report.element_error.size());

  const auto full = MmsComputeRanking(mesh, report.indicator_full, report.element_error);
  const auto standard =
      MmsComputeRanking(mesh, report.indicator_standard, report.element_error);
  CAPTURE(full.effectivity, full.correlation, full.dorfler_true_capture,
          full.marked_overlap, full.marked, full.true_marked, full.corner_marked,
          full.top_decile_overlap, full.top_decile_true_capture);
  CAPTURE(standard.effectivity, standard.correlation, standard.dorfler_true_capture,
          standard.marked_overlap, standard.marked, standard.corner_marked,
          standard.top_decile_overlap);

  // Effectivity must be within a usable band. A recovered-flux indicator is not guaranteed
  // to bracket the error, but an order-of-magnitude miss would make it unusable for AMR.
  // MEASURED: 1.18 to 2.08, i.e. mildly conservative and drifting up with order, which is
  // the expected behaviour for flux recovery on a corner-dominated problem.
  CHECK(full.effectivity > 0.5);
  CHECK(full.effectivity < 4.0);

  // What AMR actually needs is that the indicator's top decile IS the true top decile.
  // MEASURED across (order, n) = (1,4) (1,8) (2,4) (2,8):
  //   decile overlap  0.50  0.72  0.90  0.80
  //   decile capture  0.23  0.33  0.65  0.91
  //   Dorfler overlap 0.29  0.52  1.00  1.00
  // Both improve markedly with order and refinement, which is the regime that matters: at
  // order 2 the indicator's Dorfler set is EXACTLY the true one and its top decile carries
  // 91 percent of the true error. Gates are set to the weakest configuration, so they hold
  // everywhere rather than only where the indicator looks best.
  CHECK(full.top_decile_overlap > 0.45);
  CHECK(full.top_decile_true_capture > 0.2);
  CHECK(full.marked_overlap > 0.25);
  CHECK(full.dorfler_true_capture > 0.15);

  // At order 2 the agreement is much stronger, and that is worth pinning separately so a
  // regression cannot hide behind the coarse-order gates above.
  if (order == 2)
  {
    CHECK(full.top_decile_overlap > 0.75);
    CHECK(full.top_decile_true_capture > 0.6);
    CHECK(full.marked_overlap > 0.99);
  }

  // The corner is where the error concentrates, so a usable indicator must mark it.
  CHECK(full.corner_marked);

  // NOT asserted: a high global Spearman correlation. MEASURED 0.34, 0.75, 0.40, 0.44 --
  // note it does NOT improve with order while every AMR-relevant metric does. That is the
  // tell that the low values are a property of the metric, not a defect in the indicator:
  // the error concentrates so sharply that most elements carry a negligible share, and
  // their mutual ordering is noise which global Spearman weights equally with the elements
  // that matter. At order 2 / n = 4 the global correlation is 0.40 while the Dorfler marked
  // set is exactly correct. Gating on global rank correlation would reject a sound
  // indicator, so it is reported but not asserted.
}

TEST_CASE("Manufactured L-shape indicator degrades when the enrichment is sliced off",
          "[singularmms][singularelements][Serial]")
{
  // Stage 6b, the diagnostic that matters for AMR. Production feeds the estimator only the
  // standard block of the field, discarding the enrichment coefficients. If the singular
  // part of grad(u_h) is invisible to the estimator, then on precisely the elements where
  // the enrichment is doing its work the indicator is computed from an incomplete field.
  //
  // This test pins the measured consequence so a future change to the production slice has
  // something to compare against.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, n = 8;

  auto mesh = MmsLShapeMesh(n);
  const auto report = MmsSolve(n, order, 20, true);
  const auto full = MmsComputeRanking(mesh, report.indicator_full, report.element_error);
  const auto sliced =
      MmsComputeRanking(mesh, report.indicator_standard, report.element_error);

  // Corner-element indicator values under both variants.
  double full_corner = 0.0, sliced_corner = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    if (MmsCornerVertex(mesh, element) >= 0)
    {
      full_corner += report.indicator_full[element];
      sliced_corner += report.indicator_standard[element];
    }
  }
  CAPTURE(full.effectivity, sliced.effectivity, full.correlation, sliced.correlation,
          full.dorfler_true_capture, sliced.dorfler_true_capture, full.corner_marked,
          sliced.corner_marked, full_corner, sliced_corner);

  // Both variants are computed from the same mesh and the same exact errors, so any
  // difference is attributable to the slice alone.
  CHECK(report.indicator_full.size() == report.indicator_standard.size());

  // The indicators genuinely differ: slicing is not a no-op.
  CHECK(std::abs(full_corner - sliced_corner) > 1.0e-12 * std::max(full_corner, 1.0));

  // MEASURED consequence of the slice, at order 2:
  //   n = 4: top-decile overlap 0.90 with the enrichment, 0.60 without
  //   n = 8: top-decile overlap 0.795 both ways, effectivity 2.08 vs 2.28
  // So on this problem the slice degrades the indicator's static-mesh element ranking but
  // does not blind it to the corner: the sliced indicator still marks the corner element,
  // because the STANDARD part of the discrete solution is itself large and irregular there.
  //
  // IMPORTANT: worse static-mesh ranking does NOT make it worse for AMR. The stage-7
  // experiment measures the opposite -- driving refinement with this sliced indicator is
  // MORE efficient than driving it with the full-field one (final adaptive point 2.178e-3
  // at 345 DOFs versus 4.145e-3 at 353). Ranking against the exact error on a fixed mesh
  // and efficiency over an adaptive sequence are different questions, and here they
  // disagree, so single-mesh overlap is not a reliable predictor of a path-dependent
  // adaptive sequence.
  //
  // Better than full-field is still not GOOD: cost-matched against a dense quasi-uniform
  // envelope the sliced indicator trails uniform refinement by 1.54x, while an exact-error
  // oracle beats it by 0.51x. See the stage-7 test for that comparison.
  CHECK(sliced.corner_marked);
  CHECK(sliced.effectivity > full.effectivity);
  CHECK(sliced.top_decile_overlap <= full.top_decile_overlap + 1.0e-12);
}

template <typename CandidateFactory>
std::vector<double> MmsSelectSingularOnlyIndicator(const MmsSolveReport &report,
                                                   const CandidateFactory &candidate)
{
  // Preserve the production estimator and avoid even constructing hierarchical patch data
  // unless the solved system actually contains singular enrichment DOFs.
  if (report.enrichment_dofs == 0)
  {
    return report.indicator_standard;
  }
  auto indicator = candidate();
  REQUIRE(indicator.size() == report.indicator_standard.size());
  return indicator;
}

TEST_CASE("Manufactured L-shape test selector is lazy on standard systems",
          "[singularmms][singularelements][Serial]")
{
  // This is a harness-level dispatch invariant: candidate construction is lazy and returns
  // the already-computed production indicator when no singular DOFs exist. Production
  // estimator regression is supplied separately by the unchanged ordinary regression suite.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  const int n = GENERATE(2, 4);
  CAPTURE(order, n);
  auto mesh = MmsLShapeMesh(n);
  const auto report = MmsSolveOnMesh(mesh, order, 20, false, MmsCornerQuadrature::GRADED);
  REQUIRE(report.enrichment_dofs == 0);
  bool candidate_constructed = false;
  const auto selected = MmsSelectSingularOnlyIndicator(
      report,
      [&]()
      {
        candidate_constructed = true;
        return std::vector<double>(report.indicator_standard.size(), -1.0);
      });
  CHECK_FALSE(candidate_constructed);
  REQUIRE(selected.size() == report.indicator_standard.size());
  double worst = 0.0, scale = 0.0;
  for (std::size_t element = 0; element < selected.size(); element++)
  {
    worst =
        std::max(worst, std::abs(selected[element] - report.indicator_standard[element]));
    scale = std::max(scale, std::abs(report.indicator_standard[element]));
  }
  CAPTURE(worst, scale);
  CHECK(worst <= 1.0e-12 * std::max(scale, 1.0));

  const auto selected_ranking = MmsComputeRanking(mesh, selected, report.element_error);
  const auto production =
      MmsComputeRanking(mesh, report.indicator_standard, report.element_error);
  CHECK_THAT(selected_ranking.effectivity, WithinAbs(production.effectivity, 1.0e-13));
  CHECK_THAT(selected_ranking.dorfler_true_capture,
             WithinAbs(production.dorfler_true_capture, 1.0e-13));
  CHECK(selected_ranking.marked == production.marked);
}

namespace
{

// One of the three sequences the AMR experiment compares.
enum class MmsAmrVariant
{
  STANDARD,         // no enrichment at all
  ENRICHED_FULL,    // enrichment, indicator from the full enriched field
  ENRICHED_SLICED,  // enrichment, indicator from the standard block only, AS PRODUCTION
                    // DOES
  ENRICHED_ORACLE   // enrichment, marked by the EXACT element errors
};

const char *MmsAmrVariantName(MmsAmrVariant variant)
{
  switch (variant)
  {
    case MmsAmrVariant::STANDARD:
      return "standard";
    case MmsAmrVariant::ENRICHED_ORACLE:
      return "enriched-oracle";
    case MmsAmrVariant::ENRICHED_FULL:
      return "enriched-full";
    default:
      return "enriched-sliced";
  }
}

// Per-pass record, so where refinement is actually spent can be reported rather than
// inferred from whether a single corner element happened to be marked.
struct MmsAmrPass
{
  int elements = 0;
  long long dofs = 0;
  double error = 0.0;

  // Shares carried by a FIXED-RADIUS neighbourhood of the corner, so passes are comparable.
  double corner_true_share = 0.0;
  double corner_indicator_share = 0.0;

  // A SECOND, smaller radius. Reported so the mechanism can be shown not to be an artifact
  // of one particular radius choice.
  double corner_true_share_tight = 0.0;
  double corner_indicator_share_tight = 0.0;

  // The shrinking-region version, kept only to show the difference between the two
  // accountings. Not used for any conclusion.
  double incident_true_share = 0.0;

  // Marking outcome.
  std::size_t marked = 0;
  std::size_t marked_corner = 0;
  double marked_true_capture = 0.0;
};

template <typename IndicatorSelector>
std::vector<MmsAmrPass> MmsRunAmrWithIndicator(bool enriched, int order, int passes,
                                               double theta,
                                               const IndicatorSelector &select_indicator)
{
  auto mesh = MmsLShapeMesh(2);
  std::vector<MmsAmrPass> history;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto report =
        MmsSolveOnMesh(mesh, order, 20, enriched, MmsCornerQuadrature::GRADED);
    // Which indicator drives marking is supplied explicitly so candidate estimators use
    // this exact same conforming AMR loop and algebra certification as the production,
    // full-field, and oracle references.
    const std::vector<double> indicator = select_indicator(mesh, report);
    REQUIRE(indicator.size() == report.element_error.size());

    MmsAmrPass record;
    record.elements = report.elements;
    record.dofs = report.dofs;
    record.error = report.error_direct;

    // Every stage-2 and stage-3 gate must continue to hold on every adapted mesh. An AMR
    // loop that silently broke the algebra would otherwise look like a convergence result.
    CHECK(report.algebraic_residual < 1.0e-10);
    CHECK_THAT(report.error_expanded,
               WithinAbs(report.error_direct, 1.0e-4 * report.error_direct));
    CHECK_THAT(report.error_galerkin,
               WithinAbs(report.error_direct, 1.0e-4 * report.error_direct));

    // Where the error and the indicator mass actually live, measured over a FIXED PHYSICAL
    // REGION rather than over "corner-incident elements". The latter is a shrinking region:
    // after each refinement the elements touching the corner are smaller, so their share of
    // anything necessarily falls, and reading that as "the corner is resolved" would be
    // circular. The fixed radius is a geometric neighbourhood of the corner that does not
    // move, so shares across passes are comparable.
    //
    // The radius is the initial corner-element size, so on pass zero the region is exactly
    // the initial corner patch and thereafter it holds that patch's descendants.
    constexpr double corner_radius = 0.5, corner_radius_tight = 0.2;
    double error_total = 0.0, indicator_total = 0.0;
    mfem::Vector centroid;
    for (int element = 0; element < report.elements; element++)
    {
      error_total += report.element_error[element];
      indicator_total += indicator[element];
      mesh.GetElementCenter(element, centroid);
      const double distance = std::hypot(centroid(0), centroid(1));
      if (distance < corner_radius)
      {
        record.corner_true_share += report.element_error[element];
        record.corner_indicator_share += indicator[element];
      }
      if (distance < corner_radius_tight)
      {
        record.corner_true_share_tight += report.element_error[element];
        record.corner_indicator_share_tight += indicator[element];
      }
      if (MmsCornerVertex(mesh, element) >= 0)
      {
        record.incident_true_share += report.element_error[element];
      }
    }
    record.corner_true_share /= error_total;
    record.corner_indicator_share /= indicator_total;
    record.corner_true_share_tight /= error_total;
    record.corner_indicator_share_tight /= indicator_total;
    record.incident_true_share /= error_total;

    const auto marked = MmsDorflerMark(indicator, theta);
    REQUIRE(!marked.empty());
    record.marked = marked.size();
    double captured = 0.0;
    for (std::size_t element : marked)
    {
      captured += report.element_error[element];
      if (MmsCornerVertex(mesh, static_cast<int>(element)) >= 0)
      {
        record.marked_corner++;
      }
    }
    record.marked_true_capture = captured / error_total;
    history.push_back(record);

    if (pass + 1 == passes)
    {
      break;
    }
    mfem::Array<mfem::Refinement> refinements;
    for (std::size_t element : marked)
    {
      refinements.Append(mfem::Refinement(static_cast<int>(element)));
    }
    mesh.GeneralRefinement(refinements);
    REQUIRE(mesh.Conforming());
  }
  return history;
}

std::vector<MmsAmrPass> MmsRunAmr(MmsAmrVariant variant, int order, int passes,
                                  double theta)
{
  const bool enriched = variant != MmsAmrVariant::STANDARD;
  return MmsRunAmrWithIndicator(enriched, order, passes, theta,
                                [variant](mfem::Mesh &, const MmsSolveReport &report)
                                {
                                  // The sliced variant reproduces production exactly. The
                                  // oracle is the exact-error upper reference, and the full
                                  // variant measures the complete recovered field.
                                  if (variant == MmsAmrVariant::ENRICHED_ORACLE)
                                  {
                                    return report.element_error;
                                  }
                                  if (variant == MmsAmrVariant::ENRICHED_FULL)
                                  {
                                    return report.indicator_full;
                                  }
                                  return report.indicator_standard;
                                });
}

// Does `candidate` Pareto-dominate `reference`, i.e. no more DOFs AND less error?
bool MmsDominates(const std::pair<long long, double> &candidate,
                  const std::pair<long long, double> &reference)
{
  return candidate.first <= reference.first && candidate.second < reference.second;
}

}  // namespace

TEST_CASE("Manufactured L-shape naive singular-flux augmentation is identical to slicing",
          "[singularmms][singularelements][Serial]")
{
  // Rules out the obvious candidate improvement BEFORE any effort is spent prototyping it.
  //
  // "Augment the recovery with the singular flux" is only meaningful if the singular
  // coefficient is recovered INDEPENDENTLY. If instead the singular flux is re-added at the
  // coefficient the finite-element solution already carries, it cancels exactly:
  //
  //   (E_std + E_sing) - (R(E_std) + E_sing) = E_std - R(E_std)
  //
  // which is precisely the sliced production indicator. So that form cannot possibly differ
  // from production, however the recovery is tuned. A genuine augmented recovery must
  // project the full field into a global RT + S_flux space with independently recovered
  // singular coefficients, with globally H(div)-conforming singular modes rather than
  // element-local copies, preferably orthogonalized against RT.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4);
  const int order = GENERATE(1, 2);
  CAPTURE(n, order);

  const auto report = MmsSolve(n, order, 20, true);
  REQUIRE(report.indicator_naive_augmented.size() == report.indicator_standard.size());

  // Identical to round-off, element by element. This is an algebraic identity, so the
  // tolerance is a round-off tolerance rather than a modelling one.
  double worst = 0.0;
  for (std::size_t element = 0; element < report.indicator_standard.size(); element++)
  {
    const double scale = std::max(report.indicator_standard[element], 1.0e-30);
    worst = std::max(worst, std::abs(report.indicator_naive_augmented[element] -
                                     report.indicator_standard[element]) /
                                scale);
  }
  CAPTURE(worst);
  CHECK(worst < 1.0e-10);

  // And it differs from the full-field indicator, confirming the comparison is not vacuous:
  // the variants are not all trivially equal.
  double full_difference = 0.0;
  for (std::size_t element = 0; element < report.indicator_standard.size(); element++)
  {
    full_difference =
        std::max(full_difference, std::abs(report.indicator_full[element] -
                                           report.indicator_standard[element]));
  }
  CHECK(full_difference > 0.0);
}

TEST_CASE("Manufactured L-shape production marking is concentrated and mis-ranked",
          "[singularmms][singularelements][Serial]")
{
  // Step 1 of the estimator-improvement plan: quantify HOW production's marking differs
  // from the oracle's, pass by pass, so a candidate estimator has a concrete target.
  //
  // Both indicators are evaluated on the SAME mesh sequence (the oracle's), so the
  // comparison is not confounded by the two loops diverging onto different meshes.
  //
  // CRITICAL: the two Dorfler sets have DIFFERENT CARDINALITIES, so raw capture and Jaccard
  // overlap conflate two distinct defects. Jaccard penalizes a size mismatch automatically,
  // and an element outside the minimum-cardinality oracle set is not thereby wasted. Three
  // metrics separate the defects:
  //
  //   same_budget_ranking  production's capture versus the best k_p elements by EXACT
  //   error,
  //                        where k_p = |production|. Isolates RANKING quality at
  //                        production's own budget: 1.0 means production picked the k_p
  //                        best elements.
  //   same_budget_extended capture of production's top k_o indicator elements versus the
  //                        oracle set, where k_o = |oracle|. Asks whether production has
  //                        the right broader ordering but simply STOPS TOO EARLY.
  //   concentration        k_p / k_o. The under-selection metric, reported separately.
  //
  // MEASURED, order 2, theta = 0.5, on the oracle's mesh sequence:
  //
  //   pass  ne   prod_cap  oracle_cap  best@k_p  ext@k_o   RANK   EXTEND  k_p/k_o
  //     0    24    0.249      0.503     0.371    0.296    0.670   0.587    0.600
  //     1    57    0.199      0.507     0.381    0.288    0.523   0.568    0.500
  //     2    90    0.178      0.506     0.211    0.318    0.843   0.629    0.231
  //     3   169    0.164      0.516     0.199    0.422    0.823   0.818    0.222
  //   means:                                              0.715   0.651    0.388
  //
  // VERDICT: BOTH defects are present and neither dominates. Were it only concentration,
  // RANK would be about 1.0; it is 0.715, so production does not pick the best k_p elements
  // even at its own budget. Were it only an early stop, EXTEND would be about 1.0; it is
  // 0.651, so extending production's own ordering to the oracle's budget still does not
  // recover the oracle's capture. The raw capture ratio of 0.389 therefore conflates a ~28%
  // ranking loss with a ~2.6x under-selection, which is exactly why the single-metric
  // description this test replaced was inadequate.
  //
  // The trends differ, which matters for what to fix: RANK IMPROVES with refinement
  // (0.670 -> 0.823) while concentration DEGRADES (0.600 -> 0.222). So on finer meshes the
  // dominant defect is increasingly under-selection rather than mis-ranking.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 4;
  constexpr double theta = 0.5;

  auto mesh = MmsLShapeMesh(2);
  std::vector<double> overlaps, capture_ratios, outside_oracle_shares, same_budget_ranking,
      same_budget_extended, concentrations;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto report = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);
    const auto production = MmsDorflerMark(report.indicator_standard, theta);
    const auto oracle = MmsDorflerMark(report.element_error, theta);
    REQUIRE(!production.empty());
    REQUIRE(!oracle.empty());

    const std::set<std::size_t> production_set(production.begin(), production.end());
    const std::set<std::size_t> oracle_set(oracle.begin(), oracle.end());
    std::vector<std::size_t> intersection, union_set;
    std::set_intersection(production_set.begin(), production_set.end(), oracle_set.begin(),
                          oracle_set.end(), std::back_inserter(intersection));
    std::set_union(production_set.begin(), production_set.end(), oracle_set.begin(),
                   oracle_set.end(), std::back_inserter(union_set));
    overlaps.push_back(static_cast<double>(intersection.size()) /
                       static_cast<double>(union_set.size()));

    double error_total = 0.0;
    for (double e : report.element_error)
    {
      error_total += e;
    }
    const auto capture = [&](const std::vector<std::size_t> &marked)
    {
      double captured = 0.0;
      for (std::size_t element : marked)
      {
        captured += report.element_error[element];
      }
      return captured / error_total;
    };
    // Top-k elements of a value vector, largest first.
    const auto top_k = [](const std::vector<double> &value, std::size_t count)
    {
      std::vector<std::size_t> order_index(value.size());
      std::iota(order_index.begin(), order_index.end(), std::size_t{0});
      count = std::min(count, value.size());
      std::partial_sort(order_index.begin(), order_index.begin() + count, order_index.end(),
                        [&value](std::size_t a, std::size_t b)
                        { return value[a] > value[b]; });
      return std::vector<std::size_t>(order_index.begin(), order_index.begin() + count);
    };

    const double production_capture = capture(production);
    const double oracle_capture = capture(oracle);
    capture_ratios.push_back(production_capture / oracle_capture);

    // SAME BUDGET, production's own cardinality: how much could the best possible choice of
    // k_p elements have captured? This is pure ranking quality, with size held fixed.
    const double best_at_production_budget =
        capture(top_k(report.element_error, production.size()));
    same_budget_ranking.push_back(production_capture / best_at_production_budget);

    // SAME BUDGET, the oracle's cardinality: extend production's own ordering to k_o
    // elements. If this recovers the oracle's capture, production ranks correctly but stops
    // too early, and the fix is the stopping rule rather than the indicator.
    const double extended_capture =
        capture(top_k(report.indicator_standard, oracle.size()));
    same_budget_extended.push_back(extended_capture / oracle_capture);

    concentrations.push_back(static_cast<double>(production.size()) /
                             static_cast<double>(oracle.size()));

    // Elements production marks that lie outside the oracle set. NOT "false positives":
    // the oracle set is minimum-cardinality for the Dorfler target, not a ground-truth
    // classification, so an element outside it may still carry real error and be worth
    // refining.
    outside_oracle_shares.push_back(1.0 - static_cast<double>(intersection.size()) /
                                              static_cast<double>(production_set.size()));

    CAPTURE(pass, report.elements, production.size(), oracle.size(), production_capture,
            oracle_capture, best_at_production_budget, extended_capture, overlaps.back(),
            capture_ratios.back(), same_budget_ranking.back(), same_budget_extended.back(),
            concentrations.back(), outside_oracle_shares.back());

    // The oracle set is minimum-cardinality for the Dorfler target theta, so its capture
    // must reach theta while production's need not. Note this is NOT an upper bound on
    // production's capture: a LARGER set can capture more than the oracle's, which is why
    // the same-budget metrics above exist.
    CHECK(oracle_capture >= theta - 1.0e-12);

    if (pass + 1 == passes)
    {
      break;
    }
    mfem::Array<mfem::Refinement> refinements;
    for (std::size_t element : oracle)
    {
      refinements.Append(mfem::Refinement(static_cast<int>(element)));
    }
    mesh.GeneralRefinement(refinements);
  }

  const auto mean = [](const std::vector<double> &value)
  { return std::accumulate(value.begin(), value.end(), 0.0) / value.size(); };
  const double mean_overlap = mean(overlaps), mean_capture = mean(capture_ratios);
  const double mean_ranking = mean(same_budget_ranking);
  const double mean_extended = mean(same_budget_extended);
  const double mean_concentration = mean(concentrations);
  CAPTURE(overlaps, capture_ratios, same_budget_ranking, same_budget_extended,
          concentrations, outside_oracle_shares, mean_overlap, mean_capture, mean_ranking,
          mean_extended, mean_concentration);

  // Production's marked set must differ substantially from the oracle's, or the stage-7
  // efficiency gap would have no explanation in the marking itself.
  CHECK(mean_overlap < 0.9);
  CHECK(mean_capture < 0.98);

  // Under-selection: production marks materially fewer elements than the oracle.
  CHECK(mean_concentration < 0.75);

  // Both defects are present, which is why the earlier single-metric description was
  // insufficient. Ranking is imperfect at production's own budget (so it is not ONLY
  // concentration), and extending production's ordering to the oracle's budget does not
  // recover the oracle's capture (so it is not ONLY an early stop).
  CHECK(mean_ranking < 0.95);
  CHECK(mean_extended < 0.95);
}

TEST_CASE(
    "Manufactured L-shape oracle AMR beats quasi-uniform refinement but the production "
    "indicator does not",
    "[singularmms][singularelements][Serial]")
{
  // Stage 7, reached only now that the static-mesh certification passes.
  //
  // THREE sequences, because two would not separate the cause. The production estimator
  // does not see the enrichment at all -- electrostaticsolver.cpp:1592 slices the field
  // down to GetNDSpace().GetTrueVSize() before calling AddErrorIndicator -- so showing that
  // FULL-FIELD flux recovery is inefficient with enrichment says nothing yet about what
  // production actually does. The sliced sequence is what settles that.
  //
  //   1. standard         no enrichment, ordinary AMR. Validates the loop itself.
  //   2. enriched-full    enrichment, indicator from the complete enriched field.
  //   3. enriched-sliced  enrichment, indicator from the standard block only = PRODUCTION.
  //
  // MEASURED, order 2, theta = 0.5, 5 passes, against a DENSE quasi-uniform envelope
  // n in {2,3,4,5,6,8} giving enriched (71, 1.683e-2), (139, 4.873e-3), (231, 2.305e-3),
  // (347, 1.402e-3), (487, 9.801e-4), (839, 5.955e-4). The dense envelope is essential:
  // with only {2,4,8} the reference jumped 231 -> 839 DOFs, and an adaptive point near 345
  // looked "non-dominated" merely because nothing comparable existed to compare it against.
  //
  // Final points, and cost-matched ratios against the envelope interpolated to the same DOF
  // count:
  //
  //   standard         312 dof, e2 2.593e-3    dominated at NO point
  //   enriched-oracle  576 dof, e2 4.322e-4    0.51x uniform  -- BEATS uniform
  //   enriched-sliced  345 dof, e2 2.178e-3    1.54x uniform  -- LOSES to uniform
  //   enriched-full    353 dof, e2 4.145e-3    3.01x uniform  -- LOSES badly
  //
  // Three conclusions, in the order they were established:
  //
  // 1. The production (sliced) indicator is NOT efficient here. Its 345-DOF point is
  //    "non-dominated" only on a knife edge: the nearest uniform point has 347 DOFs, i.e.
  //    two more DOFs for 1.54x LESS error. Cost-matched, production trails uniform
  //    refinement at every pass (1.56x, 2.11x, 2.58x, 1.54x). An earlier version of this
  //    comment claimed "production does not fail"; that was an artifact of the sparse
  //    envelope and is wrong.
  //
  // 2. Production is nevertheless much better than full-field recovery (1.54x
  // versus 3.01x),
  //    so the full-field redundancy story does not transfer to production. That part
  //    stands.
  //
  // 3. The ORACLE, marked by exact element errors, DOES beat quasi-uniform refinement
  //    (0.88x, 0.80x, 0.51x at 213, 378, 576 DOFs). This is the load-bearing result:
  //    adaptive refinement of the enriched space is genuinely worthwhile, so the estimator
  //    IS the limiting factor. Had the oracle also lost, no estimator work could have
  //    helped and the right conclusion would have been that AMR does not pay here at all.
  //    Cost-matched, production trails the ORACLE by 1.39x to 3.00x.
  //
  // Corner accounting uses a FIXED-RADIUS neighbourhood, not "corner-incident elements".
  // The incident set shrinks geometrically under refinement, so its share of anything must
  // fall, and reading that as "the corner is resolved" would be circular. With the fixed
  // region the standard sequence gives 0.612, 0.428, 0.297, 0.226, 0.321 -- note it rises
  // again on the last pass, which the shrinking-region accounting could never show.
  // record.incident_true_share retains the shrinking version purely for that contrast.
  //
  // The loop is deliberately conforming: MmsLShapeMesh produces a conforming triangulation
  // and mfem::Mesh::GeneralRefinement on a conforming mesh performs bisection with
  // conforming closure, so no hanging-node constraints enter and the enriched space stays
  // valid. That is the same conforming baseline the wider AMR work uses.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 5;
  constexpr double theta = 0.5;

  // Uniform reference sequences. Enriched and standard differ in DOF count, so each
  // adaptive sequence is compared against the uniform sequence built in the SAME space.
  // DENSE quasi-uniform envelope. The earlier {2,4,8} reference jumped from 231 to 839
  // enriched DOFs, so an adaptive point near 345 was "non-dominated" only because nothing
  // comparable existed to compare it against. n = 5 lands at ~347 enriched DOFs, almost
  // exactly on the 345/353 adaptive points, which is the decisive comparison.
  std::map<bool, std::vector<std::pair<long long, double>>> uniform;
  for (bool enriched : {false, true})
  {
    for (int n : {2, 3, 4, 5, 6, 8})
    {
      const auto report = MmsSolve(n, order, 20, enriched);
      uniform[enriched].emplace_back(report.dofs, report.error_direct);
    }
  }

  for (const auto variant :
       {MmsAmrVariant::STANDARD, MmsAmrVariant::ENRICHED_FULL,
        MmsAmrVariant::ENRICHED_SLICED, MmsAmrVariant::ENRICHED_ORACLE})
  {
    const bool enriched = variant != MmsAmrVariant::STANDARD;
    const std::string name = MmsAmrVariantName(variant);
    const auto history = MmsRunAmr(variant, order, passes, theta);
    REQUIRE(history.size() == static_cast<std::size_t>(passes));

    std::vector<long long> dofs;
    std::vector<double> errors, corner_true, corner_indicator, marked_capture;
    std::vector<std::size_t> marked, marked_corner;
    for (const auto &pass : history)
    {
      dofs.push_back(pass.dofs);
      errors.push_back(pass.error);
      corner_true.push_back(pass.corner_true_share);
      corner_indicator.push_back(pass.corner_indicator_share);
      marked.push_back(pass.marked);
      marked_corner.push_back(pass.marked_corner);
      marked_capture.push_back(pass.marked_true_capture);
    }
    CAPTURE(name, dofs, errors, corner_true, corner_indicator, marked, marked_corner,
            marked_capture, uniform[enriched]);

    // Convergence is required of every variant; losing to uniform refinement is not the
    // same as failing to converge.
    for (std::size_t i = 0; i + 1 < errors.size(); i++)
    {
      CHECK(errors[i + 1] < errors[i]);
    }

    // DIRECT Pareto comparison: is there a uniform point that dominates an adaptive point,
    // i.e. achieves lower error at no more DOFs? Comparing the last adaptive point against
    // a much larger uniform one would be a weaker and less honest test.
    std::size_t dominated = 0;
    for (const auto &adaptive_point : {std::pair{dofs.back(), errors.back()},
                                       std::pair{dofs[passes - 2], errors[passes - 2]}})
    {
      for (const auto &uniform_point : uniform[enriched])
      {
        if (MmsDominates(uniform_point, adaptive_point))
        {
          dominated++;
          break;
        }
      }
    }
    CAPTURE(dominated);

    // Expected domination counts over the final two adaptive points, one per variant.
    switch (variant)
    {
      case MmsAmrVariant::STANDARD:
        // Textbook expectation, and it holds: no uniform point dominates any adaptive
        // point. This validates the loop and the marking, so the enriched results below
        // cannot be dismissed as a harness artifact.
        CHECK(dominated == 0);
        break;
      case MmsAmrVariant::ENRICHED_ORACLE:
        // Exact-error marking beats the dense quasi-uniform envelope outright. This is what
        // establishes that adaptivity is worth having in the enriched space at all, and
        // hence that the estimator is the limiting factor rather than AMR itself.
        CHECK(dominated == 0);
        break;
      case MmsAmrVariant::ENRICHED_FULL:
        CHECK(dominated == 2);
        break;
      default:
        // Production's sliced indicator. Its FINAL point escapes domination only by two
        // DOFs, so the domination count alone overstates it; the cost-matched check below
        // is the one that matters.
        CHECK(dominated == 1);
        break;
    }

    // COST-MATCHED efficiency, which is the honest comparison: interpolate the
    // quasi-uniform envelope in log-log to each adaptive point's own DOF count. This is
    // immune to the knife-edge artifact that made a sparse envelope look favourable to
    // production.
    const auto envelope_at = [&uniform, enriched](long long at_dofs) -> double
    {
      const auto &points = uniform[enriched];
      for (std::size_t i = 0; i + 1 < points.size(); i++)
      {
        if (points[i].first <= at_dofs && at_dofs <= points[i + 1].first)
        {
          const double t = (std::log(static_cast<double>(at_dofs)) -
                            std::log(static_cast<double>(points[i].first))) /
                           (std::log(static_cast<double>(points[i + 1].first)) -
                            std::log(static_cast<double>(points[i].first)));
          return std::exp(std::log(points[i].second) + t * (std::log(points[i + 1].second) -
                                                            std::log(points[i].second)));
        }
      }
      return std::numeric_limits<double>::quiet_NaN();
    };
    const double final_ratio = errors.back() / envelope_at(dofs.back());
    CAPTURE(final_ratio, envelope_at(dofs.back()));
    REQUIRE(std::isfinite(final_ratio));

    switch (variant)
    {
      case MmsAmrVariant::ENRICHED_ORACLE:
        // Oracle AMR is genuinely more efficient than quasi-uniform refinement: 0.51x.
        CHECK(final_ratio < 0.7);
        break;
      case MmsAmrVariant::ENRICHED_SLICED:
        // Production trails uniform refinement cost-matched: 1.54x. Pinned as a LOSS, which
        // is the corrected finding.
        CHECK(final_ratio > 1.2);
        CHECK(final_ratio < 2.0);
        break;
      case MmsAmrVariant::ENRICHED_FULL:
        // Full-field recovery is worse still: 3.01x.
        CHECK(final_ratio > 2.5);
        break;
      default:
        break;
    }

    // Production must be clearly better than full-field but clearly worse than the oracle.
    // Recorded here as the ordering that motivates the next experiment.
    if (variant == MmsAmrVariant::ENRICHED_SLICED)
    {
      CHECK(errors.back() < 0.7 * 4.145e-3);   // better than full-field at ~350 dof
      CHECK(errors.back() > 1.5 * 1.0229e-3);  // worse than oracle at ~378 dof
    }

    // Mechanism, over the FIXED-RADIUS corner region so passes are comparable. Both
    // enriched indicators assign the corner a larger share of the indicator than it holds
    // of the exact error, i.e. they over-attribute to a region the enrichment already
    // resolves.
    if (enriched && variant != MmsAmrVariant::ENRICHED_ORACLE)
    {
      for (std::size_t i = 0; i < history.size(); i++)
      {
        CAPTURE(i, corner_true[i], corner_indicator[i], history[i].corner_true_share_tight,
                history[i].corner_indicator_share_tight);
        CHECK(corner_indicator[i] > corner_true[i]);
        // The same over-attribution must hold at a SECOND, tighter radius, which is what
        // shows the mechanism is not an artifact of one radius choice. Skipped on meshes
        // too coarse to place any element centroid inside the tight radius: at h ~ 0.5 the
        // region is empty and both shares are exactly zero, which carries no information.
        if (history[i].corner_true_share_tight > 0.0)
        {
          CHECK(history[i].corner_indicator_share_tight >
                history[i].corner_true_share_tight);
        }
      }
    }
  }
}

TEST_CASE("Manufactured L-shape exponent sensitivity sweep", "[.mmsdiag][Serial]")
{
  // Diagnostic only. How much of the enrichment gain is exponent-specific, as a function of
  // refinement? A generic local basis can capture much of the corner error on coarse
  // meshes, so this separates "six extra functions help" from "the r^nu basis is right".
  for (int order : {1, 2})
  {
    for (int n : {2, 4, 8, 16})
    {
      const auto standard = MmsSolve(n, order, 20, false);
      const auto correct = MmsSolve(n, order, 20, true);
      CAPTURE(order, n);
      WARN("BASE std=" << standard.error_direct << " correct=" << correct.error_direct
                       << " reduction=" << standard.error_direct / correct.error_direct);
      for (const double nu : {0.40, 0.50, 0.55, 0.60, 0.63, 0.65, 0.68, 0.70, 0.75, 0.85})
      {
        const auto wrong = MmsSolveWrongExponent(n, order, nu);
        WARN("SWEEP nu=" << nu << " excess=" << wrong.error_direct / correct.error_direct);
      }
    }
  }
}

namespace
{

// Feasibility measurements for augmenting an RT flux space with the singular mode.
struct MmsFluxAugmentationReport
{
  double extracted_nu = 0.0;

  // Largest relative jump in the mode's normal trace across any interior face.
  double worst_normal_jump = 0.0;

  // Energy of the mode that RT CANNOT represent, as a fraction of its total energy:
  //   schur / M_SS  with  schur = M_SS - M_RS^T M_RR^{-1} M_RS.
  // Near zero means the mode is almost inside RT and augmenting buys nothing; near one
  // means the mode is nearly orthogonal to RT and the augmentation is well posed.
  double independent_energy_fraction = 0.0;

  // Synthetic round trip: recover a field built from KNOWN coefficients.
  double worst_rt_coefficient_error = 0.0;
  double singular_coefficient_error = 0.0;

  bool mass_matrix_spd = false;
};

MmsFluxAugmentationReport MmsAssessFluxAugmentation(int n, int order, int quadrature_order)
{
  auto mesh = MmsLShapeMesh(n);
  MmsFluxAugmentationReport report;

  // Provenance: nu comes from the extractor, i.e. from geometry and materials.
  const auto features =
      fem::singular::ExtractSerialLineFeatures(mesh, {MmsReentrantAttribute}, {{1, 1.0}});
  REQUIRE(features.vertices.size() == 1);
  report.extracted_nu = features.vertices[0].nu;
  const double nu = report.extracted_nu;

  // (2) Normal-trace continuity across every interior face, checked by evaluating the mode
  // through BOTH incident elements' transformations. This catches an orientation or
  // element-local evaluation bug, which is the realistic failure mode.
  {
    mfem::Vector normal(2), point_1(3), point_2(3);
    for (int face = 0; face < mesh.GetNumFaces(); face++)
    {
      const auto info = mesh.GetFaceInformation(face);
      if (!info.IsInterior())
      {
        continue;
      }
      auto *transformation = mesh.GetFaceElementTransformations(face);
      const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, quadrature_order);
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        transformation->SetAllIntPoints(&ip);
        const auto &eip1 = transformation->GetElement1IntPoint();
        const auto &eip2 = transformation->GetElement2IntPoint();
        transformation->Elem1->Transform(eip1, point_1);
        transformation->Elem2->Transform(eip2, point_2);
        // The two elements must agree on the physical point, or the comparison is
        // meaningless.
        if (std::hypot(point_1(0) - point_2(0), point_1(1) - point_2(1)) > 1.0e-10)
        {
          continue;
        }
        // Skip points essentially at the singularity, where the mode is unbounded.
        if (std::hypot(point_1(0), point_1(1)) < 1.0e-8)
        {
          continue;
        }
        mfem::CalcOrtho(transformation->Jacobian(), normal);
        const double length = normal.Norml2();
        if (!(length > 0.0))
        {
          continue;
        }
        const auto side_1 = MmsSingularFluxMode(point_1(0), point_1(1), nu);
        const auto side_2 = MmsSingularFluxMode(point_2(0), point_2(1), nu);
        const double flux_1 = (side_1[0] * normal(0) + side_1[1] * normal(1)) / length;
        const double flux_2 = (side_2[0] * normal(0) + side_2[1] * normal(1)) / length;
        const double scale = std::max({std::abs(flux_1), std::abs(flux_2), 1.0e-12});
        report.worst_normal_jump =
            std::max(report.worst_normal_jump, std::abs(flux_1 - flux_2) / scale);
      }
    }
  }

  // Combined space RT + span{mode}. Assemble M_RR, M_RS, M_SS and the load for a synthetic
  // field with KNOWN coefficients.
  mfem::RT_FECollection flux_collection(std::max(order - 1, 0), 2);
  mfem::FiniteElementSpace flux_space(&mesh, &flux_collection);
  const int flux_size = flux_space.GetVSize();

  mfem::BilinearForm mass(&flux_space);
  mass.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  mass.Assemble();
  mass.Finalize();

  // Known RT coefficients, deterministic and not derived from the exact solution.
  mfem::Vector known_rt(flux_size);
  for (int i = 0; i < flux_size; i++)
  {
    known_rt(i) = std::sin(0.7 * i + 0.3) + 0.5 * std::cos(0.11 * i);
  }
  constexpr double known_alpha = 0.8137;

  mfem::Vector cross(flux_size), load(flux_size);
  cross = 0.0;
  load = 0.0;
  double mode_mass = 0.0, mode_load = 0.0;
  {
    mfem::GridFunction known_field(&flux_space, known_rt.GetData());
    mfem::Array<int> flux_dofs;
    mfem::DenseMatrix shape;
    mfem::Vector point(3), rt_value(2);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &fe = *flux_space.GetFE(element);
      auto &T = *mesh.GetElementTransformation(element);
      mfem::DofTransformation dof_transformation;
      flux_space.GetElementVDofs(element, flux_dofs, dof_transformation);
      mfem::Vector local_cross(fe.GetDof()), local_load(fe.GetDof());
      local_cross = 0.0;
      local_load = 0.0;
      shape.SetSize(fe.GetDof(), 2);
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::DUFFY,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            T.Transform(ip, point);
            const double weight = weight_ref * T.Weight();
            fe.CalcPhysVShape(T, shape);
            const auto mode = MmsSingularFluxMode(point(0), point(1), nu);
            known_field.GetVectorValue(T, ip, rt_value);
            // Synthetic field with known coefficients.
            const double field_x = rt_value(0) + known_alpha * mode[0];
            const double field_y = rt_value(1) + known_alpha * mode[1];
            mode_mass += weight * (mode[0] * mode[0] + mode[1] * mode[1]);
            mode_load += weight * (mode[0] * field_x + mode[1] * field_y);
            for (int i = 0; i < fe.GetDof(); i++)
            {
              local_cross(i) += weight * (shape(i, 0) * mode[0] + shape(i, 1) * mode[1]);
              local_load(i) += weight * (shape(i, 0) * field_x + shape(i, 1) * field_y);
            }
          });
      dof_transformation.TransformDual(local_cross);
      dof_transformation.TransformDual(local_load);
      cross.AddElementVector(flux_dofs, local_cross);
      load.AddElementVector(flux_dofs, local_load);
    }
  }

  // (3) Schur complement of the mode against RT: how much of the mode's energy lies outside
  // RT. This is simultaneously the "orthogonalized mode is still nonzero" check and the
  // conditioning check on the augmented block.
  const auto solve_mass = [&](const mfem::Vector &rhs)
  {
    mfem::Vector x(rhs.Size());
    x = 0.0;
    mfem::GSSmoother preconditioner(mass.SpMat());
    mfem::CGSolver cg;
    cg.SetOperator(mass.SpMat());
    cg.SetPreconditioner(preconditioner);
    cg.SetRelTol(1.0e-14);
    cg.SetAbsTol(1.0e-30);
    cg.SetMaxIter(50000);
    cg.SetPrintLevel(-1);
    cg.Mult(rhs, x);
    REQUIRE(cg.GetConverged());
    return x;
  };
  const auto mass_inverse_cross = solve_mass(cross);
  const double schur = mode_mass - mfem::InnerProduct(cross, mass_inverse_cross);
  report.independent_energy_fraction = schur / mode_mass;
  // The combined mass matrix is SPD exactly when M_RR is SPD (CG converged above) and the
  // Schur complement is positive.
  report.mass_matrix_spd = schur > 0.0;

  // (4) Synthetic round trip. Solve the 2x2 block system by Schur complement and compare
  // against the known coefficients.
  if (report.mass_matrix_spd)
  {
    const auto mass_inverse_load = solve_mass(load);
    const double alpha = (mode_load - mfem::InnerProduct(cross, mass_inverse_load)) / schur;
    mfem::Vector recovered_rt(mass_inverse_load);
    recovered_rt.Add(-alpha, mass_inverse_cross);
    report.singular_coefficient_error =
        std::abs(alpha - known_alpha) / std::abs(known_alpha);
    double worst = 0.0;
    for (int i = 0; i < flux_size; i++)
    {
      worst = std::max(worst, std::abs(recovered_rt(i) - known_rt(i)));
    }
    report.worst_rt_coefficient_error = worst / known_rt.Normlinf();
  }
  return report;
}

}  // namespace

namespace
{

// Global two-level hierarchical benchmark. Solve in V_p (+) S and V_{p+1} (+) S on the SAME
// mesh, inject the coarse solution exactly with diag(P, I), and use
//
//   eta_K^2 = integral_K |grad(delta u)|^2,   delta u = u_{p+1} - I u_p
//
// as the candidate indicator. The injection is exactly diag(P, I) because the enriched
// space is a DIRECT SUM whose singular block does not depend on p:
// BuildSerialTriangleDofTopology enumerates one node-gradient pair per singular-vertex
// incidence at order one regardless of the standard polynomial order, so the two problems
// share an identical enrichment DOF set and the singular block of the injection is the
// identity. Only the standard block needs MFEM's order-p to order-(p+1) prolongation.
//
// Unlike flux recovery this benchmark sees BOTH unresolved smooth error and an incorrect
// singular AMPLITUDE, because delta u is free to adjust the singular coefficient.

using MmsSparseColumn = fem::hierarchical::SparseColumn;
using MmsLocalElementData = fem::hierarchical::LocalOperatorContribution;

struct MmsLocalFineSystem
{
  int standard_size = 0;
  int enrichment_size = 0;
  std::vector<bool> essential;
  fem::singular::TriangleDofTopology topology;
  std::vector<MmsLocalElementData> elements;
};

MmsLocalFineSystem MmsAssembleLocalFineSystem(mfem::Mesh &mesh, int order,
                                              int quadrature_order)
{
  MmsLocalFineSystem system;
  mfem::H1_FECollection collection(order, 2);
  mfem::FiniteElementSpace space(&mesh, &collection);
  mfem::ND_FECollection nd_collection(order, 2);
  mfem::FiniteElementSpace nd_space(&mesh, &nd_collection);
  system.standard_size = space.GetVSize();
  const auto features =
      fem::singular::ExtractSerialLineFeatures(mesh, {MmsReentrantAttribute}, {{1, 1.0}});
  system.topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
  system.enrichment_size = static_cast<int>(system.topology.h1_dofs.size());
  const int combined_size = system.standard_size + system.enrichment_size;
  system.essential.assign(combined_size, false);
  mfem::Array<int> standard_essential;
  space.GetBoundaryTrueDofs(standard_essential);
  for (int dof : standard_essential)
  {
    system.essential[dof] = true;
  }
  const auto numbering =
      fem::singular::BuildParallelDofNumbering(Mpi::World(), system.topology);
  std::set<std::array<fem::singular::GlobalVertexId, 2>> boundary_segments;
  mfem::Array<int> vertices;
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    mesh.GetBdrElementVertices(boundary, vertices);
    std::array<fem::singular::GlobalVertexId, 2> segment{vertices[0], vertices[1]};
    std::sort(segment.begin(), segment.end());
    boundary_segments.insert(segment);
  }
  const auto singular_essential = fem::singular::GetEssentialTriangleH1TrueDofsOnSegments(
      Mpi::World(),
      std::vector<std::array<fem::singular::GlobalVertexId, 2>>(boundary_segments.begin(),
                                                                boundary_segments.end()),
      system.topology, numbering);
  for (int dof : singular_essential)
  {
    system.essential[system.standard_size + dof] = true;
  }

  mfem::DiffusionIntegrator diffusion;
  MmsSourceCoefficient source;
  mfem::DomainLFIntegrator load(source);
  load.SetIntRule(&mfem::IntRules.Get(mfem::Geometry::TRIANGLE, quadrature_order));
  const fem::singular::AdaptiveAssemblyOptions options{8, 1.0e-10, 1.0e-10, 8};
  system.elements.resize(mesh.GetNE());
  mfem::Array<int> standard_dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &data = system.elements[element];
    data.support_element = element;
    const auto &singular_dofs = system.topology.elements[element].h1;
    space.GetElementDofs(element, standard_dofs);
    mfem::FiniteElementSpace::AdjustVDofs(standard_dofs);
    for (int dof : standard_dofs)
    {
      data.dofs.push_back(dof);
    }
    for (const auto &dof : singular_dofs)
    {
      data.dofs.push_back(system.standard_size + static_cast<int>(dof.dof));
    }
    const int standard_count = standard_dofs.Size();
    const int local_size = static_cast<int>(data.dofs.size());
    data.matrix.SetSize(local_size);
    data.matrix = 0.0;
    data.rhs.SetSize(local_size);
    data.rhs = 0.0;
    mfem::DenseMatrix standard_matrix;
    mfem::Vector standard_rhs;
    const auto &fe = *space.GetFE(element);
    auto &T = *mesh.GetElementTransformation(element);
    diffusion.AssembleElementMatrix(fe, T, standard_matrix);
    load.AssembleRHSElementVect(fe, T, standard_rhs);
    for (int row = 0; row < standard_count; row++)
    {
      data.rhs(row) = standard_rhs(row);
      for (int column = 0; column < standard_count; column++)
      {
        data.matrix(row, column) = standard_matrix(row, column);
      }
    }
    if (!singular_dofs.empty())
    {
      const auto singular = fem::singular::AssembleTriangleElementEnrichmentMatrices(
          system.topology.elements[element], T, options);
      const auto coupling =
          fem::singular::AssembleTriangleElementStandardEnrichmentMatrices(
              system.topology.elements[element], fe, *nd_space.GetFE(element), T, options);
      REQUIRE(coupling.h1_standard_enrichment.Height() == standard_count);
      REQUIRE(coupling.h1_standard_enrichment.Width() ==
              static_cast<int>(singular_dofs.size()));
      for (int row = 0; row < standard_count; row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(row, standard_count + column) =
              coupling.h1_standard_enrichment(row, column);
          data.matrix(standard_count + column, row) =
              coupling.h1_enrichment_standard(column, row);
        }
      }
      for (int row = 0; row < static_cast<int>(singular_dofs.size()); row++)
      {
        for (int column = 0; column < static_cast<int>(singular_dofs.size()); column++)
        {
          data.matrix(standard_count + row, standard_count + column) =
              singular.h1_diffusion(row, column);
        }
      }
      mfem::Vector point(3);
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            T.Transform(ip, point);
            const double weight = weight_ref * T.Weight();
            const double value = MmsF(point(0), point(1));
            for (int local = 0; local < static_cast<int>(singular_dofs.size()); local++)
            {
              const auto &basis = singular_dofs[local].basis;
              data.rhs(standard_count + local) +=
                  weight * value *
                  fem::singular::EvaluateTriangleNodeGradientPotential(
                      lambda, basis.nodes[0], basis.nodes[1], basis.nu);
            }
          });
    }
  }
  return system;
}

struct MmsSparseLocalReport
{
  std::vector<double> indicator;
  double energy = 0.0;
  double work = 0.0;
  double maximum_patch_residual = 0.0;
  double maximum_patch_condition = 0.0;
  double maximum_actual_to_owner_support = 0.0;
  int edge_patches = 0;
  int interior_patches = 0;
  int owned_modes = 0;
  int maximum_patch_dimension = 0;
  int maximum_owner_elements = 0;
  int maximum_actual_support_elements = 0;
  int maximum_element_overlap = 0;
  double mean_element_overlap = 0.0;
  long long sparse_basis_nonzeros = 0;
  long long dense_basis_entries_avoided = 0;
  long long touched_matrix_nonzeros = 0;
  long long sum_patch_dimension_squared = 0;
  long long sum_patch_dimension_cubed = 0;
  long long fine_global_solves = 0;
  long long fine_global_matrix_assemblies = 0;
};

std::vector<double>
Mms2DCorrectionIndicator(mfem::Mesh &mesh, int order, int quadrature_order,
                         const mfem::Vector &correction, int enrichment_size,
                         MmsCornerQuadrature scheme = MmsCornerQuadrature::GRADED)
{
  mfem::H1_FECollection collection(order, 2);
  mfem::FiniteElementSpace space(&mesh, &collection);
  const int standard_size = space.GetVSize();
  REQUIRE(correction.Size() == standard_size + enrichment_size);
  const auto features =
      fem::singular::ExtractSerialLineFeatures(mesh, {MmsReentrantAttribute}, {{1, 1.0}});
  const auto topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
  REQUIRE(static_cast<int>(topology.h1_dofs.size()) == enrichment_size);
  mfem::GridFunction standard(&space);
  for (int i = 0; i < standard_size; i++)
  {
    standard(i) = correction(i);
  }
  mfem::Vector singular(std::max(enrichment_size, 1));
  singular = 0.0;
  for (int i = 0; i < enrichment_size; i++)
  {
    singular(i) = correction(standard_size + i);
  }
  std::vector<double> indicator(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, scheme,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          mfem::Vector gradient(2);
          standard.GetGradient(T, gradient);
          double gx = gradient(0), gy = gradient(1);
          if (!topology.elements[element].h1.empty())
          {
            double determinant;
            const auto grad_lambda =
                fem::singular::GetAffineTriangleBarycentricGradients(T, determinant);
            const auto value = fem::singular::EvaluateElementTriangleH1Enrichment(
                topology.elements[element], singular, lambda, grad_lambda);
            gx += value.gradient[0];
            gy += value.gradient[1];
          }
          indicator[element] += weight_ref * T.Weight() * (gx * gx + gy * gy);
        });
  }
  return indicator;
}

MmsSparseLocalReport MmsBuildSparseLocalCandidate(mfem::Mesh &mesh, int order,
                                                  int quadrature_order,
                                                  const MmsSolveReport &coarse)
{
  const long long assemblies_before = MmsGlobalFineMatrixAssemblies;
  const long long solves_before = MmsGlobalFineSolves;
  // Assemble uneliminated element contributions only. The candidate API has no fine
  // solution or global fine matrix, and both zero counters are regression invariants.
  const auto fine = MmsAssembleLocalFineSystem(mesh, order + 1, quadrature_order);
  MmsSparseLocalReport report;
  const int enrichment = static_cast<int>(coarse.enrichment_dofs);
  REQUIRE(fine.enrichment_size == coarse.enrichment_dofs);
  mfem::H1_FECollection coarse_collection(order, 2), fine_collection(order + 1, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  const int fine_combined = fine_standard + enrichment;

  // Assemble the p-injection from element transfer matrices. Duplicate element writes to a
  // conforming global row must agree; no global unit-vector transfer applications are used.
  std::vector<std::map<int, double>> injection_map(coarse_standard);
  mfem::IsoparametricTransformation identity;
  identity.SetIdentityTransformation(mfem::Geometry::TRIANGLE);
  mfem::DenseMatrix local_transfer;
  mfem::Array<int> coarse_dofs, fine_dofs;
  double maximum_transfer_mismatch = 0.0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    coarse_space.GetElementDofs(element, coarse_dofs);
    fine_space.GetElementDofs(element, fine_dofs);
    mfem::FiniteElementSpace::AdjustVDofs(coarse_dofs);
    mfem::FiniteElementSpace::AdjustVDofs(fine_dofs);
    fine_space.GetFE(element)->GetTransferMatrix(*coarse_space.GetFE(element), identity,
                                                 local_transfer);
    REQUIRE(local_transfer.Height() == fine_dofs.Size());
    REQUIRE(local_transfer.Width() == coarse_dofs.Size());
    for (int local_column = 0; local_column < coarse_dofs.Size(); local_column++)
    {
      auto &column = injection_map[coarse_dofs[local_column]];
      for (int local_row = 0; local_row < fine_dofs.Size(); local_row++)
      {
        const double value = local_transfer(local_row, local_column);
        if (std::abs(value) <= 1.0e-14)
        {
          continue;
        }
        const int row = fine_dofs[local_row];
        const auto [it, inserted] = column.emplace(row, value);
        if (!inserted)
        {
          maximum_transfer_mismatch =
              std::max(maximum_transfer_mismatch, std::abs(it->second - value));
        }
      }
    }
  }
  CHECK(maximum_transfer_mismatch < 1.0e-12);
  std::vector<MmsSparseColumn> injection(coarse_standard);
  for (int column = 0; column < coarse_standard; column++)
  {
    for (const auto &[row, value] : injection_map[column])
    {
      injection[column].dofs.push_back(row);
      injection[column].values.push_back(value);
    }
  }

  mfem::Vector injected(fine_combined);
  injected = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection[column].dofs.size(); entry++)
    {
      injected(injection[column].dofs[entry]) +=
          coarse.standard_coefficients(column) * injection[column].values[entry];
    }
  }
  for (int e = 0; e < enrichment; e++)
  {
    injected(fine_standard + e) = coarse.enrichment_block(e);
  }

  // Full combined residual assembled from uneliminated element contributions, including
  // standard-singular coupling. This is linear work and constructs no global fine matrix.
  mfem::Vector residual = fem::hierarchical::AssembleResidual(fine_combined, fine.elements,
                                                              injected, fine.essential);
  std::vector<std::set<int>> dof_elements(fine_combined);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    for (int dof : fine.elements[element].dofs)
    {
      dof_elements[dof].insert(element);
    }
  }

  const auto &topology = fine.topology;
  mfem::Array<int> element_dofs;

  const auto sparse_entity_complement =
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
      for (int local_row = 0; local_row < static_cast<int>(rows.size()); local_row++)
      {
        const auto found = injection_map[column].find(rows[local_row]);
        if (found != injection_map[column].end())
        {
          vector(local_row) = found->second;
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
    std::vector<MmsSparseColumn> result;
    for (const auto &local : complement)
    {
      MmsSparseColumn column;
      for (int row = 0; row < static_cast<int>(rows.size()); row++)
      {
        if (std::abs(local(row)) > 1.0e-14)
        {
          column.dofs.push_back(rows[row]);
          column.values.push_back(local(row));
        }
      }
      result.push_back(std::move(column));
    }
    return result;
  };

  struct Patch
  {
    int owned = 0;
    std::vector<MmsSparseColumn> basis;
    std::vector<int> coarse_guests;
    std::vector<int> singular_guests;
    mfem::Vector coefficients;
    std::set<int> support_elements;
  };
  std::vector<Patch> patches;
  std::vector<int> element_overlap(mesh.GetNE(), 0);

  const auto add_patch = [&](std::vector<MmsSparseColumn> owned,
                             const std::vector<int> &owner_elements, bool edge_patch)
  {
    if (owned.empty())
    {
      return;
    }
    Patch patch;
    patch.owned = static_cast<int>(owned.size());
    patch.basis = std::move(owned);
    report.owned_modes += patch.owned;
    report.edge_patches += edge_patch;
    report.interior_patches += !edge_patch;
    report.maximum_owner_elements =
        std::max(report.maximum_owner_elements, static_cast<int>(owner_elements.size()));
    std::set<int> coarse_guests, singular_guests;
    for (int element : owner_elements)
    {
      coarse_space.GetElementDofs(element, element_dofs);
      mfem::FiniteElementSpace::AdjustVDofs(element_dofs);
      for (int dof : element_dofs)
      {
        if (!coarse.essential_mask[dof])
        {
          coarse_guests.insert(dof);
        }
      }
      for (const auto &dof : topology.elements[element].h1)
      {
        if (!fine.essential[fine_standard + static_cast<int>(dof.dof)])
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
      MmsSparseColumn column;
      column.dofs.push_back(fine_standard + dof);
      column.values.push_back(1.0);
      patch.singular_guests.push_back(dof);
      patch.basis.push_back(std::move(column));
    }
    const int dimension = static_cast<int>(patch.basis.size());
    report.maximum_patch_dimension = std::max(report.maximum_patch_dimension, dimension);
    report.sum_patch_dimension_squared += 1LL * dimension * dimension;
    report.sum_patch_dimension_cubed += 1LL * dimension * dimension * dimension;
    bool touches_essential = false;
    for (const auto &column : patch.basis)
    {
      report.sparse_basis_nonzeros += column.dofs.size();
      report.dense_basis_entries_avoided += fine_combined - column.dofs.size();
      for (int dof : column.dofs)
      {
        touches_essential = touches_essential || fine.essential[dof];
        patch.support_elements.insert(dof_elements[dof].begin(), dof_elements[dof].end());
      }
    }
    REQUIRE_FALSE(touches_essential);
    report.maximum_actual_support_elements =
        std::max(report.maximum_actual_support_elements,
                 static_cast<int>(patch.support_elements.size()));
    report.maximum_actual_to_owner_support =
        std::max(report.maximum_actual_to_owner_support,
                 static_cast<double>(patch.support_elements.size()) /
                     std::max<std::size_t>(owner_elements.size(), 1));
    for (int element : patch.support_elements)
    {
      element_overlap[element]++;
    }

    // Element-local congruence assembly Q_K^T A_K Q_K over the complete measured support
    // union. Complete conforming guests are never truncated to the nominal owner star.
    mfem::DenseMatrix restricted;
    mfem::Vector patch_rhs;
    report.touched_matrix_nonzeros += fem::hierarchical::AssembleRestrictedOperator(
        fine.elements, patch.support_elements, patch.basis, residual, restricted,
        patch_rhs);
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
    mfem::Array<int> fine_entity, coarse_entity;
    fine_space.GetEdgeInteriorDofs(edge, fine_entity);
    coarse_space.GetEdgeInteriorDofs(edge, coarse_entity);
    add_patch(sparse_entity_complement(fine_entity, coarse_entity), edge_elements[edge],
              true);
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::Array<int> fine_entity, coarse_entity;
    fine_space.GetElementInteriorDofs(element, fine_entity);
    coarse_space.GetElementInteriorDofs(element, coarse_entity);
    add_patch(sparse_entity_complement(fine_entity, coarse_entity), {element}, false);
  }
  int coarse_free = 0, fine_free = 0;
  for (int i = 0; i < coarse_standard; i++)
  {
    coarse_free += !coarse.essential_mask[i];
  }
  for (int i = 0; i < fine_standard; i++)
  {
    fine_free += !fine.essential[i];
  }
  REQUIRE(report.owned_modes == fine_free - coarse_free);
  long long overlap_sum = 0;
  for (int overlap : element_overlap)
  {
    report.maximum_element_overlap = std::max(report.maximum_element_overlap, overlap);
    overlap_sum += overlap;
  }
  report.mean_element_overlap =
      static_cast<double>(overlap_sum) / std::max(mesh.GetNE(), 1);

  // Hierarchical-key reconstruction: own new modes once and average repeated guest
  // coefficients. A constant number of global-length vectors is used for the injected
  // coarse field, assembled residual, and final correction; patch storage remains sparse.
  mfem::Vector raw(fine_combined);
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
    if (coarse_count[dof] == 0)
    {
      continue;
    }
    const double coefficient = coarse_sum[dof] / coarse_count[dof];
    for (std::size_t entry = 0; entry < injection[dof].dofs.size(); entry++)
    {
      raw(injection[dof].dofs[entry]) += coefficient * injection[dof].values[entry];
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    if (singular_count[dof] > 0)
    {
      raw(fine_standard + dof) += singular_sum[dof] / singular_count[dof];
    }
  }
  const double raw_energy = fem::hierarchical::Energy(fine.elements, raw);
  const double raw_work = mfem::InnerProduct(raw, residual);
  const double alpha = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
  raw *= alpha;
  report.energy = alpha * alpha * raw_energy;
  report.work = alpha * raw_work;
  report.indicator =
      Mms2DCorrectionIndicator(mesh, order + 1, quadrature_order, raw, enrichment);
  report.fine_global_matrix_assemblies = MmsGlobalFineMatrixAssemblies - assemblies_before;
  report.fine_global_solves = MmsGlobalFineSolves - solves_before;
  return report;
}

struct MmsTwoLevelReport
{
  std::vector<double> indicator;   // eta_K^2 = int_K |grad(delta u)|^2
  double galerkin_residual = 0.0;  // max |a(delta u, v_p)| over the coarse basis, relative
  double energy_identity = 0.0;    // |(||delta u||_a^2) - (e_p^2 - e_{p+1}^2)| relative
  double delta_energy = 0.0;
  double coarse_error = 0.0;
  double fine_error = 0.0;

  // ABLATIONS that isolate the coarse response. In hierarchical coordinates the correction
  // solves
  //   [A_cc A_cn][dc]   [0]
  //   [A_nc A_nn][dn] = [r]
  // so dc = -A_cc^{-1} A_cn dn: the EXISTING coarse standard and singular coefficients
  // change even though the enrichment DOF set is identical across p. A localized scheme
  // that solves only for new DOFs omits exactly that response, and a patch method is really
  // approximating the Schur complement A_nn - A_nc A_cc^{-1} A_cn. These two ablations
  // measure what is lost.
  //
  // frozen_enrichment: p+1 solve with the singular coefficients held at their coarse
  // values,
  //                    i.e. the singular part of dc suppressed.
  // new_dofs_only:     correction restricted to the genuinely new polynomial DOFs, i.e. all
  // of
  //                    dc suppressed.
  std::vector<double> indicator_frozen_enrichment;
  std::vector<double> indicator_new_dofs_only;
  double frozen_enrichment_energy = 0.0;
  double new_dofs_only_energy = 0.0;

  // Practical candidate: independent edge/interior hierarchical patches. Each patch owns
  // its entity's new polynomial modes and includes the incident coarse standard and
  // singular basis functions as local response directions. The patch residual is always
  // restricted from the FULL combined residual, including standard-singular coupling.
  std::vector<double> indicator_local_patches;
  double algebraic_delta_energy = 0.0;
  double global_correction_residual = 0.0;
  double global_correction_energy_identity = 0.0;
  double maximum_patch_residual = 0.0;
  double maximum_patch_energy_identity = 0.0;
  double maximum_patch_energy_ratio = 0.0;
  double overlap_partition_error = 0.0;
  int minimum_overlap = 0;
  int maximum_overlap = 0;
  int edge_patches = 0;
  int interior_patches = 0;
  int owned_hierarchical_modes = 0;
  int maximum_patch_dimension = 0;
  int maximum_patch_elements = 0;
  int uncovered_coarse_guests = 0;
  int uncovered_singular_guests = 0;
  double maximum_patch_condition = 0.0;
  double local_patch_energy = 0.0;
  double local_patch_captured_fraction = 0.0;
  double local_indicator_energy_identity = 0.0;
};

MmsTwoLevelReport MmsTwoLevelBenchmark(mfem::Mesh &mesh, int order, int quadrature_order)
{
  MmsTwoLevelReport report;
  const auto coarse =
      MmsSolveOnMesh(mesh, order, quadrature_order, true, MmsCornerQuadrature::GRADED);
  const auto fine =
      MmsSolveOnMesh(mesh, order + 1, quadrature_order, true, MmsCornerQuadrature::GRADED);
  report.coarse_error = coarse.error_direct;
  report.fine_error = fine.error_direct;
  // Same mesh, and the enrichment block is p-independent, so the DOF sets must match.
  REQUIRE(coarse.enrichment_dofs == fine.enrichment_dofs);

  mfem::H1_FECollection coarse_collection(order, 2), fine_collection(order + 1, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);

  // Exact injection of the STANDARD block, p -> p+1 on the same mesh.
  mfem::GridFunction coarse_standard(&coarse_space);
  for (int i = 0; i < coarse_space.GetVSize(); i++)
  {
    coarse_standard(i) = coarse.standard_coefficients(i);
  }
  mfem::GridFunction injected(&fine_space);
  {
    mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
    transfer.Mult(coarse_standard, injected);
  }
  // Verify the injection really is exact: the injected function must agree pointwise with
  // the coarse one. Otherwise delta u would contain injection error rather than
  // discretization error.
  {
    double worst = 0.0;
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, 8);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      auto &T = *mesh.GetElementTransformation(element);
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        T.SetIntPoint(&ip);
        worst = std::max(
            worst, std::abs(injected.GetValue(T, ip) - coarse_standard.GetValue(T, ip)));
      }
    }
    const double scale = std::max(coarse_standard.Normlinf(), 1.0e-30);
    CAPTURE(worst / scale);
    CHECK(worst / scale < 1.0e-10);
  }

  // delta u = u_{p+1} - I u_p, in the fine standard space plus the shared singular block.
  mfem::GridFunction delta_standard(&fine_space);
  for (int i = 0; i < fine_space.GetVSize(); i++)
  {
    delta_standard(i) = fine.standard_coefficients(i) - injected(i);
  }
  mfem::Vector delta_enrichment(std::max(coarse.enrichment_dofs, 1LL));
  delta_enrichment = 0.0;
  for (int i = 0; i < coarse.enrichment_dofs; i++)
  {
    delta_enrichment(i) = fine.enrichment_block(i) - coarse.enrichment_block(i);
  }

  // Singular topology on this mesh, shared by both orders.
  const auto features =
      fem::singular::ExtractSerialLineFeatures(mesh, {MmsReentrantAttribute}, {{1, 1.0}});
  const auto topology = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);

  // grad(delta u) at a quadrature point, standard part plus singular part.
  const auto delta_gradient = [&](int element, mfem::ElementTransformation &T,
                                  const mfem::IntegrationPoint &ip,
                                  const fem::singular::TriangleBarycentricPoint &lambda)
  {
    mfem::Vector standard(2);
    delta_standard.GetGradient(T, standard);
    std::array<double, 2> value{standard(0), standard(1)};
    if (!topology.elements[element].h1.empty())
    {
      double jacobian_determinant;
      const auto grad_lambda =
          fem::singular::GetAffineTriangleBarycentricGradients(T, jacobian_determinant);
      T.SetIntPoint(&ip);
      const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
          topology.elements[element], delta_enrichment, lambda, grad_lambda);
      value[0] += singular.gradient[0];
      value[1] += singular.gradient[1];
    }
    return value;
  };

  report.indicator.assign(mesh.GetNE(), 0.0);
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &T = *mesh.GetElementTransformation(element);
    double local = 0.0;
    MmsForEachQuadraturePoint(
        mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
        [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
        {
          mfem::IntegrationPoint ip;
          ip.Set2(lambda[1], lambda[2]);
          T.SetIntPoint(&ip);
          const double weight = weight_ref * T.Weight();
          const auto gradient = delta_gradient(element, T, ip, lambda);
          local += weight * (gradient[0] * gradient[0] + gradient[1] * gradient[1]);
        });
    report.indicator[element] = local;
    report.delta_energy += local;
  }

  // ABLATION 1: p+1 solve with the enrichment coefficients FROZEN at their coarse values.
  // This suppresses the singular part of the coarse response dc, leaving the standard part.
  {
    mfem::Vector frozen(std::max(coarse.enrichment_dofs, 1LL));
    frozen = 0.0;
    for (int i = 0; i < coarse.enrichment_dofs; i++)
    {
      frozen(i) = coarse.enrichment_block(i);
    }
    const auto frozen_fine = MmsSolveOnMesh(mesh, order + 1, quadrature_order, true,
                                            MmsCornerQuadrature::GRADED, 0.0, &frozen);
    mfem::GridFunction frozen_delta(&fine_space);
    for (int i = 0; i < fine_space.GetVSize(); i++)
    {
      frozen_delta(i) = frozen_fine.standard_coefficients(i) - injected(i);
    }
    mfem::Vector zero_enrichment(std::max(coarse.enrichment_dofs, 1LL));
    zero_enrichment = 0.0;
    report.indicator_frozen_enrichment.assign(mesh.GetNE(), 0.0);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      auto &T = *mesh.GetElementTransformation(element);
      double local = 0.0;
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            const double weight = weight_ref * T.Weight();
            mfem::Vector gradient(2);
            frozen_delta.GetGradient(T, gradient);
            local += weight * (gradient(0) * gradient(0) + gradient(1) * gradient(1));
          });
      report.indicator_frozen_enrichment[element] = local;
      report.frozen_enrichment_energy += local;
    }
  }

  // ABLATION 2: correction restricted to a hierarchical complement of the injected coarse
  // space. Done properly this time.
  //
  // The previous version was WRONG in two ways, and the symptom was that its energy
  // exceeded the full correction by 4 to 5 times -- impossible, because for the true
  // residual and ANY restricted space W the SPD variational bound gives
  //
  //   eta_W^2 = r_W^T A_W^{-1} r_W <= eta_full^2.
  //
  // The two defects were: (i) the residual was formed from the standard stiffness alone, so
  // it omitted the singular coupling action A_se u_s; (ii) "new" fine DOFs were identified
  // by looking for a prolongation entry equal to one, but the p-refinement range is not a
  // coordinate subset of the fine space, so that never constructed a valid complement.
  //
  // Correct construction: build the combined injection P_c = diag(P, I) explicitly, form an
  // explicit complement Q whose columns span a subspace complementary to range(P_c), then
  //
  //   r_f = b_f - A_f P_c u_c,   A_nn = Q^T A_f Q,   r_n = Q^T r_f,   delta_new = Q
  //   A_nn^{-1} r_n.
  {
    const int coarse_standard = coarse_space.GetVSize();
    const int fine_standard = fine_space.GetVSize();
    const int enrichment = static_cast<int>(coarse.enrichment_dofs);
    const int fine_combined = fine_standard + enrichment;

    // Dense injection P_c = diag(P, I), obtained by injecting each coarse unit vector.
    // Dense is acceptable at these sizes and avoids assuming anything about the transfer's
    // sparsity.
    mfem::DenseMatrix injection(fine_combined, coarse_standard + enrichment);
    injection = 0.0;
    {
      mfem::PRefinementTransferOperator transfer(coarse_space, fine_space);
      mfem::GridFunction unit_coarse(&coarse_space), image(&fine_space);
      for (int c = 0; c < coarse_standard; c++)
      {
        unit_coarse = 0.0;
        unit_coarse(c) = 1.0;
        transfer.Mult(unit_coarse, image);
        for (int f = 0; f < fine_standard; f++)
        {
          injection(f, c) = image(f);
        }
      }
      for (int e = 0; e < enrichment; e++)
      {
        injection(fine_standard + e, coarse_standard + e) = 1.0;
      }
    }

    // Free fine DOFs. The complement is built inside the free block so that essential DOFs
    // never enter the restricted solve.
    std::vector<int> free_fine, fine_to_free(fine_combined, -1);
    for (int i = 0; i < fine_combined; i++)
    {
      if (!fine.essential_mask[i])
      {
        fine_to_free[i] = static_cast<int>(free_fine.size());
        free_fine.push_back(i);
      }
    }
    const int free_size = static_cast<int>(free_fine.size());

    // Injected coarse solution as a full fine-space vector.
    mfem::Vector coarse_combined(coarse_standard + enrichment);
    for (int i = 0; i < coarse_standard; i++)
    {
      coarse_combined(i) = coarse.standard_coefficients(i);
    }
    for (int e = 0; e < enrichment; e++)
    {
      coarse_combined(coarse_standard + e) = coarse.enrichment_block(e);
    }
    mfem::Vector injected_combined(fine_combined);
    injection.Mult(coarse_combined, injected_combined);

    // TRUE residual r_f = b_f - A_f P_c u_c, using the FULL combined operator so the
    // singular coupling A_se u_s is included.
    mfem::Vector residual(fine_combined);
    fine.combined_matrix->Mult(injected_combined, residual);
    for (int i = 0; i < fine_combined; i++)
    {
      residual(i) = fine.combined_rhs(i) - residual(i);
    }

    // Algebraic global correction. This is the reference for every restricted-space energy
    // bound below. Unlike report.delta_energy, both sides here use the SAME assembled
    // operator, so the variational bound is an exact algebraic statement rather than a
    // comparison between assembly and independent quadrature.
    mfem::Vector global_correction(fine.combined_solution_vector);
    global_correction -= injected_combined;
    {
      mfem::Vector applied(fine_combined), solve_residual(fine_combined);
      fine.combined_matrix->Mult(global_correction, applied);
      solve_residual = applied;
      solve_residual -= residual;
      for (int i = 0; i < fine_combined; i++)
      {
        if (fine.essential_mask[i])
        {
          // Essential rows are not part of the variational problem.
          solve_residual(i) = 0.0;
        }
      }
      const double residual_norm = [&]()
      {
        double sum = 0.0;
        for (int i : free_fine)
        {
          sum += residual(i) * residual(i);
        }
        return std::sqrt(sum);
      }();
      report.global_correction_residual =
          solve_residual.Norml2() / std::max(residual_norm, 1.0e-30);
      report.algebraic_delta_energy = mfem::InnerProduct(global_correction, applied);
      const double work = mfem::InnerProduct(global_correction, residual);
      report.global_correction_energy_identity =
          std::abs(report.algebraic_delta_energy - work) /
          std::max({std::abs(report.algebraic_delta_energy), std::abs(work), 1.0e-30});
    }

    // Explicit complement Q of range(P_c) within the free block, by Gram-Schmidt: project
    // the free coordinate directions against an orthonormal basis of the injected range and
    // keep whatever survives.
    mfem::DenseMatrix range_basis(free_size, coarse_standard + enrichment);
    range_basis = 0.0;
    for (int c = 0; c < coarse_standard + enrichment; c++)
    {
      for (int row = 0; row < free_size; row++)
      {
        range_basis(row, c) = injection(free_fine[row], c);
      }
    }
    std::vector<mfem::Vector> range_orthonormal;
    for (int c = 0; c < range_basis.Width(); c++)
    {
      mfem::Vector column(free_size);
      for (int row = 0; row < free_size; row++)
      {
        column(row) = range_basis(row, c);
      }
      for (const auto &basis : range_orthonormal)
      {
        column.Add(-mfem::InnerProduct(column, basis), basis);
      }
      const double norm = column.Norml2();
      if (norm > 1.0e-10)
      {
        column /= norm;
        range_orthonormal.push_back(column);
      }
    }
    std::vector<mfem::Vector> complement;
    for (int direction = 0; direction < free_size; direction++)
    {
      mfem::Vector column(free_size);
      column = 0.0;
      column(direction) = 1.0;
      for (const auto &basis : range_orthonormal)
      {
        column.Add(-mfem::InnerProduct(column, basis), basis);
      }
      for (const auto &basis : complement)
      {
        column.Add(-mfem::InnerProduct(column, basis), basis);
      }
      const double norm = column.Norml2();
      if (norm > 1.0e-8)
      {
        column /= norm;
        complement.push_back(column);
      }
    }
    const int complement_size = static_cast<int>(complement.size());
    // Dimensions must add up: free = dim(range within free) + dim(complement).
    CAPTURE(free_size, range_orthonormal.size(), complement_size);
    CHECK(static_cast<std::size_t>(complement_size) + range_orthonormal.size() ==
          static_cast<std::size_t>(free_size));

    mfem::Vector new_correction(fine_combined);
    new_correction = 0.0;
    if (complement_size > 0)
    {
      // A_nn = Q^T A_f Q and r_n = Q^T r_f, both restricted to the free block.
      mfem::DenseMatrix restricted(complement_size);
      restricted = 0.0;
      mfem::Vector restricted_rhs(complement_size);
      restricted_rhs = 0.0;
      std::vector<mfem::Vector> applied(complement_size);
      for (int j = 0; j < complement_size; j++)
      {
        mfem::Vector expanded(fine_combined), product(fine_combined);
        expanded = 0.0;
        for (int row = 0; row < free_size; row++)
        {
          expanded(free_fine[row]) = complement[j](row);
        }
        fine.combined_matrix->Mult(expanded, product);
        applied[j].SetSize(free_size);
        for (int row = 0; row < free_size; row++)
        {
          applied[j](row) = product(free_fine[row]);
        }
        double dot = 0.0;
        for (int row = 0; row < free_size; row++)
        {
          dot += complement[j](row) * residual(free_fine[row]);
        }
        restricted_rhs(j) = dot;
      }
      for (int i = 0; i < complement_size; i++)
      {
        for (int j = 0; j < complement_size; j++)
        {
          restricted(i, j) = mfem::InnerProduct(complement[i], applied[j]);
        }
      }
      mfem::DenseMatrixInverse inverse(restricted);
      mfem::Vector coefficients(complement_size);
      inverse.Mult(restricted_rhs, coefficients);
      for (int j = 0; j < complement_size; j++)
      {
        for (int row = 0; row < free_size; row++)
        {
          new_correction(free_fine[row]) += coefficients(j) * complement[j](row);
        }
      }
    }

    // Energy of the restricted correction in the a-norm, via the operator rather than by
    // quadrature, so the SPD bound below is checked on exactly the quantity it bounds.
    {
      mfem::Vector product(fine_combined);
      fine.combined_matrix->Mult(new_correction, product);
      report.new_dofs_only_energy = mfem::InnerProduct(new_correction, product);
    }

    // Elementwise indicator from the restricted correction, standard plus singular parts.
    mfem::GridFunction new_standard(&fine_space);
    for (int i = 0; i < fine_standard; i++)
    {
      new_standard(i) = new_correction(i);
    }
    mfem::Vector new_enrichment(std::max(enrichment, 1));
    new_enrichment = 0.0;
    for (int e = 0; e < enrichment; e++)
    {
      new_enrichment(e) = new_correction(fine_standard + e);
    }
    report.indicator_new_dofs_only.assign(mesh.GetNE(), 0.0);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      auto &T = *mesh.GetElementTransformation(element);
      double local = 0.0;
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            const double weight = weight_ref * T.Weight();
            mfem::Vector standard(2);
            new_standard.GetGradient(T, standard);
            double gx = standard(0), gy = standard(1);
            if (!topology.elements[element].h1.empty())
            {
              double jacobian_determinant;
              const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
                  T, jacobian_determinant);
              T.SetIntPoint(&ip);
              const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
                  topology.elements[element], new_enrichment, lambda, grad_lambda);
              gx += singular.gradient[0];
              gy += singular.gradient[1];
            }
            local += weight * (gx * gx + gy * gy);
          });
      report.indicator_new_dofs_only[element] = local;
    }

    // LOCAL EDGE/INTERIOR HIERARCHICAL PATCHES. The genuinely new p+1 modes are owned once
    // by their mesh entity. Each patch also contains all free coarse standard basis
    // functions incident on the owning entity's elements, plus incident singular basis
    // functions. Including those guest directions approximates the coarse response which
    // the new-DOFs-only ablation drops.
    struct LocalPatch
    {
      std::vector<mfem::Vector> basis;
      int owned_modes = 0;
      std::vector<int> guest_coarse_dofs;
      std::vector<int> guest_singular_dofs;
      mfem::Vector coefficients;
      mfem::Vector correction;
      double energy = 0.0;
    };
    std::vector<LocalPatch> patches;

    // Build an entity-local Euclidean complement of the injected coarse entity-interior
    // range inside the fine entity-interior coordinates. For p -> p+1 this gives exactly
    // one new mode per free edge and p-1 new modes per element interior. Modified
    // Gram-Schmidt is repeated because these tiny nodal trace blocks can be poorly scaled.
    const auto entity_complement =
        [&](mfem::Array<int> fine_entity, mfem::Array<int> coarse_entity)
    {
      mfem::FiniteElementSpace::AdjustVDofs(fine_entity);
      mfem::FiniteElementSpace::AdjustVDofs(coarse_entity);
      std::vector<int> fine_rows, coarse_columns;
      for (int dof : fine_entity)
      {
        if (!fine.essential_mask[dof])
        {
          fine_rows.push_back(dof);
        }
      }
      for (int dof : coarse_entity)
      {
        if (!coarse.essential_mask[dof])
        {
          coarse_columns.push_back(dof);
        }
      }
      std::vector<mfem::Vector> range;
      for (int column : coarse_columns)
      {
        mfem::Vector vector(static_cast<int>(fine_rows.size()));
        for (int row = 0; row < vector.Size(); row++)
        {
          vector(row) = injection(fine_rows[row], column);
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
      for (int direction = 0; direction < static_cast<int>(fine_rows.size()); direction++)
      {
        mfem::Vector vector(static_cast<int>(fine_rows.size()));
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
      std::vector<mfem::Vector> expanded;
      for (const auto &local : complement)
      {
        mfem::Vector vector(fine_combined);
        vector = 0.0;
        for (int row = 0; row < static_cast<int>(fine_rows.size()); row++)
        {
          vector(fine_rows[row]) = local(row);
        }
        expanded.push_back(vector);
      }
      return expanded;
    };

    // Element incidence of every global edge.
    std::vector<std::vector<int>> edge_elements(mesh.GetNEdges());
    {
      mfem::Array<int> edges, orientations;
      for (int element = 0; element < mesh.GetNE(); element++)
      {
        mesh.GetElementEdges(element, edges, orientations);
        for (int edge : edges)
        {
          edge_elements[edge].push_back(element);
        }
      }
    }

    const auto add_patch = [&](std::vector<mfem::Vector> owned_modes,
                               const std::vector<int> &elements, bool edge_patch)
    {
      if (owned_modes.empty())
      {
        return;
      }
      LocalPatch patch;
      patch.basis = std::move(owned_modes);
      patch.owned_modes = static_cast<int>(patch.basis.size());
      report.owned_hierarchical_modes += patch.owned_modes;
      if (edge_patch)
      {
        report.edge_patches++;
      }
      else
      {
        report.interior_patches++;
      }

      // Guest coarse standard directions: complete conforming injected basis columns, not
      // element-truncated copies. Truncation would invalidate the restricted-space bound.
      std::set<int> incident_coarse;
      std::set<int> incident_enrichment;
      mfem::Array<int> element_dofs;
      for (int element : elements)
      {
        coarse_space.GetElementDofs(element, element_dofs);
        mfem::FiniteElementSpace::AdjustVDofs(element_dofs);
        for (int dof : element_dofs)
        {
          if (!coarse.essential_mask[dof])
          {
            incident_coarse.insert(dof);
          }
        }
        for (const auto &dof : topology.elements[element].h1)
        {
          const int index = fine_standard + static_cast<int>(dof.dof);
          if (!fine.essential_mask[index])
          {
            incident_enrichment.insert(static_cast<int>(dof.dof));
          }
        }
      }
      for (int coarse_dof : incident_coarse)
      {
        mfem::Vector column(fine_combined);
        for (int row = 0; row < fine_combined; row++)
        {
          column(row) = injection(row, coarse_dof);
        }
        patch.guest_coarse_dofs.push_back(coarse_dof);
        patch.basis.push_back(column);
      }
      for (int singular_dof : incident_enrichment)
      {
        mfem::Vector column(fine_combined);
        column = 0.0;
        column(fine_standard + singular_dof) = 1.0;
        patch.guest_singular_dofs.push_back(singular_dof);
        patch.basis.push_back(column);
      }

      // Essential constraints are imposed on the patch SPACE, before solving. Every basis
      // column must vanish there; overwriting a solved correction afterwards would solve a
      // different variational problem.
      double essential_maximum = 0.0;
      for (const auto &column : patch.basis)
      {
        for (int row = 0; row < fine_combined; row++)
        {
          if (fine.essential_mask[row])
          {
            essential_maximum = std::max(essential_maximum, std::abs(column(row)));
          }
        }
      }
      CHECK(essential_maximum < 1.0e-13);

      const int patch_size = static_cast<int>(patch.basis.size());
      report.maximum_patch_dimension = std::max(report.maximum_patch_dimension, patch_size);
      report.maximum_patch_elements =
          std::max(report.maximum_patch_elements, static_cast<int>(elements.size()));
      mfem::DenseMatrix restricted(patch_size);
      restricted = 0.0;
      mfem::Vector restricted_rhs(patch_size);
      restricted_rhs = 0.0;
      std::vector<mfem::Vector> applied(patch_size);
      for (int column = 0; column < patch_size; column++)
      {
        applied[column].SetSize(fine_combined);
        fine.combined_matrix->Mult(patch.basis[column], applied[column]);
        restricted_rhs(column) = mfem::InnerProduct(patch.basis[column], residual);
      }
      for (int row = 0; row < patch_size; row++)
      {
        for (int column = 0; column < patch_size; column++)
        {
          restricted(row, column) = mfem::InnerProduct(patch.basis[row], applied[column]);
        }
      }
      // Every restricted operator must remain SPD and reasonably conditioned. Requesting
      // Cholesky factors (spd=true) catches a rank-deficient entity complement or duplicate
      // guest key. ||A||_F ||A^{-1}||_F is an inexpensive condition-number upper proxy that
      // is available in this no-LAPACK unit-test build.
      mfem::DenseMatrixInverse inverse(restricted, true);
      mfem::DenseMatrix inverse_matrix;
      inverse.GetInverseMatrix(inverse_matrix);
      report.maximum_patch_condition = std::max(
          report.maximum_patch_condition, restricted.FNorm() * inverse_matrix.FNorm());
      patch.coefficients.SetSize(patch_size);
      inverse.Mult(restricted_rhs, patch.coefficients);

      mfem::Vector patch_residual(patch_size);
      restricted.Mult(patch.coefficients, patch_residual);
      patch_residual -= restricted_rhs;
      const double relative_residual =
          patch_residual.Norml2() / std::max(restricted_rhs.Norml2(), 1.0e-30);
      report.maximum_patch_residual =
          std::max(report.maximum_patch_residual, relative_residual);

      patch.correction.SetSize(fine_combined);
      patch.correction = 0.0;
      for (int column = 0; column < patch_size; column++)
      {
        patch.correction.Add(patch.coefficients(column), patch.basis[column]);
      }
      mfem::Vector applied_correction(fine_combined);
      fine.combined_matrix->Mult(patch.correction, applied_correction);
      patch.energy = mfem::InnerProduct(patch.correction, applied_correction);
      const double work = mfem::InnerProduct(patch.correction, residual);
      // The 1e-4 denominator floor makes a 1e-10 normalized check equivalent to the
      // scale-aware identity |E-W| <= 1e-10 max(|E|,|W|) + 1e-14. Nearly inactive
      // patches otherwise turn harmless roundoff into a large relative error.
      const double identity = std::abs(patch.energy - work) /
                              (std::max(std::abs(patch.energy), std::abs(work)) + 1.0e-4);
      report.maximum_patch_energy_identity =
          std::max(report.maximum_patch_energy_identity, identity);
      report.maximum_patch_energy_ratio =
          std::max(report.maximum_patch_energy_ratio,
                   patch.energy / std::max(report.algebraic_delta_energy, 1.0e-30));
      CHECK(patch.energy >= -1.0e-14);
      CHECK(patch.energy <= report.algebraic_delta_energy * (1.0 + 1.0e-8) + 1.0e-14);

      REQUIRE(patch.owned_modes + static_cast<int>(patch.guest_coarse_dofs.size()) +
                  static_cast<int>(patch.guest_singular_dofs.size()) ==
              patch_size);
      patches.push_back(std::move(patch));
    };

    // Edge patches own trace modes; boundary-edge modes disappear when essential fine
    // coordinates are filtered. Element patches own only interior modes, so ownership is
    // disjoint by construction.
    for (int edge = 0; edge < mesh.GetNEdges(); edge++)
    {
      mfem::Array<int> fine_edge, coarse_edge;
      fine_space.GetEdgeInteriorDofs(edge, fine_edge);
      coarse_space.GetEdgeInteriorDofs(edge, coarse_edge);
      add_patch(entity_complement(fine_edge, coarse_edge), edge_elements[edge], true);
    }
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      mfem::Array<int> fine_interior, coarse_interior;
      fine_space.GetElementInteriorDofs(element, fine_interior);
      coarse_space.GetElementInteriorDofs(element, coarse_interior);
      add_patch(entity_complement(fine_interior, coarse_interior), {element}, false);
    }

    int coarse_free_standard = 0, fine_free_standard = 0;
    for (int dof = 0; dof < coarse_standard; dof++)
    {
      coarse_free_standard += !coarse.essential_mask[dof];
    }
    for (int dof = 0; dof < fine_standard; dof++)
    {
      fine_free_standard += !fine.essential_mask[dof];
    }
    CAPTURE(report.edge_patches, report.interior_patches, report.owned_hierarchical_modes,
            coarse_free_standard, fine_free_standard);
    CHECK(report.owned_hierarchical_modes == fine_free_standard - coarse_free_standard);
    REQUIRE(!patches.empty());

    // Partition of unity in HIERARCHICAL COORDINATES, not fine nodal rows. Each owned
    // edge/interior complement coefficient is inserted exactly once. Repeated response
    // coefficients are averaged by their stable coarse/singular DOF key and only then
    // reconstructed through diag(P,I). This distinction matters: an injected coarse guest
    // column can be nonzero on the fine rows used by another patch's owned mode, so
    // nodal-row averaging would incorrectly attenuate the uniquely owned hierarchical
    // contribution.
    std::vector<int> coarse_multiplicity(coarse_standard, 0);
    std::vector<int> singular_multiplicity(enrichment, 0);
    std::vector<double> coarse_sum(coarse_standard, 0.0);
    std::vector<double> singular_sum(enrichment, 0.0);
    mfem::Vector raw_correction(fine_combined);
    raw_correction = 0.0;
    for (const auto &patch : patches)
    {
      REQUIRE(patch.coefficients.Size() == static_cast<int>(patch.basis.size()));
      for (int owned = 0; owned < patch.owned_modes; owned++)
      {
        raw_correction.Add(patch.coefficients(owned), patch.basis[owned]);
      }
      int coefficient = patch.owned_modes;
      for (int coarse_dof : patch.guest_coarse_dofs)
      {
        coarse_sum[coarse_dof] += patch.coefficients(coefficient++);
        coarse_multiplicity[coarse_dof]++;
      }
      for (int singular_dof : patch.guest_singular_dofs)
      {
        singular_sum[singular_dof] += patch.coefficients(coefficient++);
        singular_multiplicity[singular_dof]++;
      }
      REQUIRE(coefficient == patch.coefficients.Size());
    }

    report.minimum_overlap = std::numeric_limits<int>::max();
    report.maximum_overlap = 0;
    for (int coarse_dof = 0; coarse_dof < coarse_standard; coarse_dof++)
    {
      const int count = coarse_multiplicity[coarse_dof];
      if (count == 0)
      {
        continue;
      }
      report.minimum_overlap = std::min(report.minimum_overlap, count);
      report.maximum_overlap = std::max(report.maximum_overlap, count);
      const double coefficient = coarse_sum[coarse_dof] / static_cast<double>(count);
      for (int row = 0; row < fine_combined; row++)
      {
        raw_correction(row) += coefficient * injection(row, coarse_dof);
      }
      double weight_sum = 0.0;
      for (const auto &patch : patches)
      {
        weight_sum +=
            static_cast<double>(std::count(patch.guest_coarse_dofs.begin(),
                                           patch.guest_coarse_dofs.end(), coarse_dof)) /
            static_cast<double>(count);
      }
      report.overlap_partition_error =
          std::max(report.overlap_partition_error, std::abs(weight_sum - 1.0));
    }
    for (int coarse_dof = 0; coarse_dof < coarse_standard; coarse_dof++)
    {
      if (!coarse.essential_mask[coarse_dof] && coarse_multiplicity[coarse_dof] == 0)
      {
        report.uncovered_coarse_guests++;
      }
    }
    for (int singular_dof = 0; singular_dof < enrichment; singular_dof++)
    {
      if (!coarse.essential_mask[coarse_standard + singular_dof] &&
          singular_multiplicity[singular_dof] == 0)
      {
        report.uncovered_singular_guests++;
      }
      const int count = singular_multiplicity[singular_dof];
      if (count == 0)
      {
        continue;
      }
      report.minimum_overlap = std::min(report.minimum_overlap, count);
      report.maximum_overlap = std::max(report.maximum_overlap, count);
      raw_correction(fine_standard + singular_dof) +=
          singular_sum[singular_dof] / static_cast<double>(count);
      double weight_sum = 0.0;
      for (const auto &patch : patches)
      {
        weight_sum +=
            static_cast<double>(std::count(patch.guest_singular_dofs.begin(),
                                           patch.guest_singular_dofs.end(), singular_dof)) /
            static_cast<double>(count);
      }
      report.overlap_partition_error =
          std::max(report.overlap_partition_error, std::abs(weight_sum - 1.0));
    }
    REQUIRE(report.minimum_overlap != std::numeric_limits<int>::max());

    // The partition-of-unity sum is an admissible direction, not itself a Galerkin solve.
    // An exact scalar line search gives the best correction along that direction and a
    // certified captured energy in [0, E_global]. The scalar does not change element
    // ranking.
    mfem::Vector raw_applied(fine_combined);
    fine.combined_matrix->Mult(raw_correction, raw_applied);
    const double raw_energy = mfem::InnerProduct(raw_correction, raw_applied);
    const double raw_work = mfem::InnerProduct(raw_correction, residual);
    const double alpha = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
    mfem::Vector local_correction(raw_correction);
    local_correction *= alpha;
    mfem::Vector local_applied(fine_combined);
    fine.combined_matrix->Mult(local_correction, local_applied);
    report.local_patch_energy = mfem::InnerProduct(local_correction, local_applied);
    const double local_work = mfem::InnerProduct(local_correction, residual);
    report.local_patch_captured_fraction =
        report.local_patch_energy / std::max(report.algebraic_delta_energy, 1.0e-30);
    CHECK(report.local_patch_energy >= -1.0e-14);
    CHECK(report.local_patch_energy <=
          report.algebraic_delta_energy * (1.0 + 1.0e-8) + 1.0e-14);
    CHECK_THAT(
        report.local_patch_energy,
        WithinAbs(local_work, 1.0e-10 * std::max({std::abs(report.local_patch_energy),
                                                  std::abs(local_work), 1.0}) +
                                  1.0e-14));

    mfem::GridFunction local_standard(&fine_space);
    for (int i = 0; i < fine_standard; i++)
    {
      local_standard(i) = local_correction(i);
    }
    mfem::Vector local_enrichment(std::max(enrichment, 1));
    local_enrichment = 0.0;
    for (int e = 0; e < enrichment; e++)
    {
      local_enrichment(e) = local_correction(fine_standard + e);
    }
    report.indicator_local_patches.assign(mesh.GetNE(), 0.0);
    double integrated_local_energy = 0.0;
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      auto &T = *mesh.GetElementTransformation(element);
      double local = 0.0;
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            const double weight = weight_ref * T.Weight();
            mfem::Vector standard(2);
            local_standard.GetGradient(T, standard);
            double gx = standard(0), gy = standard(1);
            if (!topology.elements[element].h1.empty())
            {
              double jacobian_determinant;
              const auto grad_lambda = fem::singular::GetAffineTriangleBarycentricGradients(
                  T, jacobian_determinant);
              T.SetIntPoint(&ip);
              const auto singular = fem::singular::EvaluateElementTriangleH1Enrichment(
                  topology.elements[element], local_enrichment, lambda, grad_lambda);
              gx += singular.gradient[0];
              gy += singular.gradient[1];
            }
            local += weight * (gx * gx + gy * gy);
          });
      report.indicator_local_patches[element] = local;
      integrated_local_energy += local;
    }
    report.local_indicator_energy_identity =
        std::abs(integrated_local_energy - report.local_patch_energy) /
        std::max(report.local_patch_energy, 1.0e-30);
  }

  // (a) Galerkin orthogonality a(delta u, v_p) = 0 for every coarse basis function. Because
  // the coarse space is nested in the fine one, both u_p and u_{p+1} satisfy the same
  // variational problem against v_p, so their difference is a-orthogonal to all of V_p (+)
  // S. Tested by integrating grad(delta u) . grad(v_p) for each coarse standard basis
  // function.
  {
    double worst = 0.0, scale = 0.0;
    mfem::Array<int> coarse_dofs;
    mfem::Vector coarse_shape_gradient;
    mfem::DenseMatrix coarse_dshape;
    std::vector<double> functional(coarse_space.GetVSize(), 0.0);
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &fe = *coarse_space.GetFE(element);
      auto &T = *mesh.GetElementTransformation(element);
      coarse_space.GetElementDofs(element, coarse_dofs);
      coarse_dshape.SetSize(fe.GetDof(), 2);
      MmsForEachQuadraturePoint(
          mesh, element, quadrature_order, MmsCornerQuadrature::GRADED,
          [&](const fem::singular::TriangleBarycentricPoint &lambda, double weight_ref)
          {
            mfem::IntegrationPoint ip;
            ip.Set2(lambda[1], lambda[2]);
            T.SetIntPoint(&ip);
            const double weight = weight_ref * T.Weight();
            fe.CalcPhysDShape(T, coarse_dshape);
            const auto gradient = delta_gradient(element, T, ip, lambda);
            for (int i = 0; i < fe.GetDof(); i++)
            {
              functional[coarse_dofs[i]] += weight * (coarse_dshape(i, 0) * gradient[0] +
                                                      coarse_dshape(i, 1) * gradient[1]);
            }
          });
    }
    // Exclude constrained boundary DOFs: delta u is not orthogonal to those, and they are
    // not part of the trial space.
    mfem::Array<int> essential;
    coarse_space.GetBoundaryTrueDofs(essential);
    std::set<int> essential_set(essential.begin(), essential.end());
    for (int dof = 0; dof < coarse_space.GetVSize(); dof++)
    {
      if (essential_set.count(dof))
      {
        continue;
      }
      worst = std::max(worst, std::abs(functional[dof]));
      scale = std::max(scale, std::abs(functional[dof]));
    }
    report.galerkin_residual = worst / std::max(report.delta_energy, 1.0e-30);
  }

  // (b) Energy identity ||delta u||_a^2 = ||u - u_p||_a^2 - ||u - u_{p+1}||_a^2, which
  // holds exactly for nested spaces. This is what certifies delta u as a measure of the
  // error REMOVED by the p increment rather than an arbitrary difference.
  {
    const double predicted = report.coarse_error - report.fine_error;
    report.energy_identity =
        std::abs(report.delta_energy - predicted) / std::max(std::abs(predicted), 1.0e-30);
  }
  return report;
}

}  // namespace

TEST_CASE(
    "Manufactured L-shape two-level hierarchical difference is a certified error measure",
    "[singularmms][singularelements][Serial]")
{
  // Certify the two-level benchmark BEFORE using it as an indicator. Both identities below
  // follow from nesting, so a violation means the injection or the enrichment bookkeeping
  // is wrong, not that the estimator is poor.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4);
  CAPTURE(n);
  auto mesh = MmsLShapeMesh(n);

  const auto report = MmsTwoLevelBenchmark(mesh, 1, 20);
  CAPTURE(report.delta_energy, report.algebraic_delta_energy, report.coarse_error,
          report.fine_error, report.galerkin_residual, report.energy_identity,
          report.global_correction_residual, report.global_correction_energy_identity,
          report.maximum_patch_residual, report.maximum_patch_energy_identity,
          report.maximum_patch_energy_ratio, report.overlap_partition_error,
          report.minimum_overlap, report.maximum_overlap, report.edge_patches,
          report.interior_patches, report.owned_hierarchical_modes,
          report.maximum_patch_dimension, report.maximum_patch_elements,
          report.maximum_patch_condition, report.uncovered_coarse_guests,
          report.uncovered_singular_guests, report.local_patch_energy,
          report.local_patch_captured_fraction, report.local_indicator_energy_identity);

  // The p increment must actually reduce the error, or the identity below is vacuous.
  CHECK(report.fine_error < report.coarse_error);

  // (a) a(delta u, v_p) = 0.
  CHECK(report.galerkin_residual < 1.0e-8);

  // (b) ||delta u||_a^2 = e_p^2 - e_{p+1}^2.
  CHECK(report.energy_identity < 1.0e-3);

  // The indicator must sum to the delta energy.
  const double summed =
      std::accumulate(report.indicator.begin(), report.indicator.end(), 0.0);
  CHECK_THAT(summed, WithinRel(report.delta_energy, 1.0e-12));

  // Patch correctness invariants. The full and every restricted solve use the same
  // assembled operator and full combined residual, so these are algebraic requirements,
  // not estimator-quality claims.
  CHECK(report.global_correction_residual < 1.0e-10);
  CHECK(report.global_correction_energy_identity < 1.0e-10);
  CHECK(report.maximum_patch_residual < 1.0e-10);
  CHECK(report.maximum_patch_energy_identity < 1.0e-10);
  CHECK(report.maximum_patch_energy_ratio <= 1.0 + 1.0e-8);
  CHECK(report.overlap_partition_error < 1.0e-14);
  CHECK(report.minimum_overlap >= 1);
  CHECK(report.maximum_overlap >= report.minimum_overlap);
  CHECK(report.edge_patches > 0);
  // This p=1 -> p=2 triangle increment gains no interior modes.
  CHECK(report.interior_patches == 0);
  CHECK(report.owned_hierarchical_modes > 0);
  CHECK(report.maximum_patch_dimension > 0);
  CHECK(report.maximum_patch_elements <= 2);
  CHECK(std::isfinite(report.maximum_patch_condition));
  CHECK(report.uncovered_coarse_guests == 0);
  CHECK(report.uncovered_singular_guests == 0);
  CHECK(report.local_patch_energy >= 0.0);
  CHECK(report.local_patch_captured_fraction >= 0.0);
  CHECK(report.local_patch_captured_fraction <= 1.0 + 1.0e-8);
  CHECK(report.local_indicator_energy_identity < 1.0e-7);
  const double local_summed = std::accumulate(report.indicator_local_patches.begin(),
                                              report.indicator_local_patches.end(), 0.0);
  CHECK_THAT(local_summed, WithinRel(report.local_patch_energy, 1.0e-7));
}

TEST_CASE(
    "Manufactured L-shape two-level hierarchical marking is compared against production",
    "[singularmms][singularelements][Serial]")
{
  // With the two-level difference certified, evaluate it as an INDICATOR using the same
  // diagnostics as production and the oracle: RANK at its own budget, EXTEND at the
  // oracle's budget, and cardinality. Measured on the oracle mesh sequence so no AMR loop
  // is needed to decide whether local machinery is worth building.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 3;
  constexpr double theta = 0.5;

  auto mesh = MmsLShapeMesh(2);
  std::vector<double> rank_production, rank_hierarchical, extend_production,
      extend_hierarchical, concentration_production, concentration_hierarchical,
      rank_new_only, extend_new_only, concentration_new_only, energy_ratio_new_only;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto report = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);
    const auto two_level = MmsTwoLevelBenchmark(mesh, order, 20);
    REQUIRE(two_level.indicator.size() == report.element_error.size());
    // The identities must still hold on every adapted mesh, not only the structured ones.
    CHECK(two_level.galerkin_residual < 1.0e-6);
    CHECK(two_level.energy_identity < 1.0e-2);

    double error_total = 0.0;
    for (double e : report.element_error)
    {
      error_total += e;
    }
    const auto capture = [&](const std::vector<std::size_t> &marked)
    {
      double captured = 0.0;
      for (std::size_t element : marked)
      {
        captured += report.element_error[element];
      }
      return captured / error_total;
    };
    const auto top_k = [](const std::vector<double> &value, std::size_t count)
    {
      std::vector<std::size_t> index(value.size());
      std::iota(index.begin(), index.end(), std::size_t{0});
      count = std::min(count, value.size());
      std::partial_sort(index.begin(), index.begin() + count, index.end(),
                        [&value](std::size_t a, std::size_t b)
                        { return value[a] > value[b]; });
      return std::vector<std::size_t>(index.begin(), index.begin() + count);
    };
    const auto oracle = MmsDorflerMark(report.element_error, theta);
    const double oracle_capture = capture(oracle);
    const auto measure = [&](const std::vector<double> &indicator, double &rank,
                             double &extend, double &concentration)
    {
      const auto marked = MmsDorflerMark(indicator, theta);
      REQUIRE(!marked.empty());
      rank = capture(marked) / capture(top_k(report.element_error, marked.size()));
      extend = capture(top_k(indicator, oracle.size())) / oracle_capture;
      concentration =
          static_cast<double>(marked.size()) / static_cast<double>(oracle.size());
    };
    double rank_p, extend_p, concentration_p, rank_h, extend_h, concentration_h;
    measure(report.indicator_standard, rank_p, extend_p, concentration_p);
    measure(two_level.indicator, rank_h, extend_h, concentration_h);
    double rank_f, extend_f, concentration_f, rank_nd, extend_nd, concentration_nd;
    measure(two_level.indicator_frozen_enrichment, rank_f, extend_f, concentration_f);
    measure(two_level.indicator_new_dofs_only, rank_nd, extend_nd, concentration_nd);
    CAPTURE(two_level.frozen_enrichment_energy, two_level.new_dofs_only_energy, rank_f,
            rank_nd, extend_f, extend_nd, concentration_f, concentration_nd);

    // ABLATION 1: freezing the singular coefficients costs almost nothing here. Energy
    // stays within 0.23, 0.07 and 0.008 percent of the full correction and every marking
    // metric is identical to six digits. So on THIS problem the singular part of the coarse
    // response is negligible: the p increment barely readjusts the singular amplitude, and
    // hierarchical marking does NOT succeed because it corrects that amplitude. Recorded
    // because it contradicts the natural assumption, and because it means a localized
    // scheme need not reproduce the singular response to match the benchmark HERE -- though
    // a fluxonium-class problem with a strongly hybridized singular mode could behave
    // differently.
    CHECK_THAT(two_level.frozen_enrichment_energy,
               WithinRel(two_level.delta_energy, 1.0e-2));
    CHECK_THAT(rank_f, WithinAbs(rank_h, 1.0e-5));
    CHECK_THAT(concentration_f, WithinAbs(concentration_h, 1.0e-12));

    // ABLATION 2, with the corrected complement construction. The SPD variational bound
    //   eta_W^2 = r_W^T A_W^{-1} r_W <= eta_full^2
    // holds for the true residual and ANY restricted space W, so this invariant is a hard
    // correctness check on the ablation itself, not a property of the problem. The previous
    // implementation violated it by 4 to 5 times, which is what exposed its two bugs: a
    // residual missing the singular coupling A_se u_s, and a "new DOF" test that assumed
    // the p-refinement range is a coordinate subset of the fine space.
    CHECK(two_level.new_dofs_only_energy <=
          two_level.delta_energy * (1.0 + 1.0e-8) + 1.0e-14);

    rank_production.push_back(rank_p);
    rank_hierarchical.push_back(rank_h);
    extend_production.push_back(extend_p);
    extend_hierarchical.push_back(extend_h);
    concentration_production.push_back(concentration_p);
    concentration_hierarchical.push_back(concentration_h);
    rank_new_only.push_back(rank_nd);
    extend_new_only.push_back(extend_nd);
    concentration_new_only.push_back(concentration_nd);
    energy_ratio_new_only.push_back(two_level.new_dofs_only_energy /
                                    two_level.delta_energy);
    if (pass + 1 == passes)
    {
      break;
    }
    mfem::Array<mfem::Refinement> refinements;
    for (std::size_t element : oracle)
    {
      refinements.Append(mfem::Refinement(static_cast<int>(element)));
    }
    mesh.GeneralRefinement(refinements);
  }

  const auto mean = [](const std::vector<double> &v)
  { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); };
  const double mean_rank_p = mean(rank_production), mean_rank_h = mean(rank_hierarchical);
  const double mean_extend_p = mean(extend_production);
  const double mean_extend_h = mean(extend_hierarchical);
  const double mean_conc_p = mean(concentration_production);
  const double mean_conc_h = mean(concentration_hierarchical);

  CAPTURE(rank_production, rank_hierarchical, extend_production, extend_hierarchical,
          concentration_production, concentration_hierarchical, mean_rank_p, mean_rank_h,
          mean_extend_p, mean_extend_h, mean_conc_p, mean_conc_h);

  // MEASURED, order 2, theta = 0.5, on the oracle mesh sequence. Oracle-optimal is 1.0:
  //
  //   metric      production   hierarchical
  //   RANK          0.679        0.9999
  //   EXTEND        0.595        1.0000
  //   k_p/k_o       0.444        1.199
  //
  // The hierarchical difference is essentially oracle-quality on all three, and in
  // particular it removes the under-selection: it marks slightly MORE than the minimal
  // oracle set rather than a fifth of it.
  //
  // HOW MUCH OF THIS IS BY CONSTRUCTION. The certified identity is
  // ||delta u||_a^2 = e_p^2 - e_{p+1}^2. On these order-2 runs the squared errors are
  //
  //   ne = 24   e_p^2 = 1.683e-2   e_{p+1}^2 = 1.434e-3
  //   ne = 57   e_p^2 = 5.724e-3   e_{p+1}^2 = 6.352e-4
  //   ne = 90   e_p^2 = 2.280e-3   e_{p+1}^2 = 2.425e-4
  //
  // so delta u captures 88.9 to 91.5 percent of the SQUARED error, while the remaining
  // ratio in ENERGY NORM is 29 to 33 percent. Those are different numbers because the
  // energy norm is the square root; quoting the squared ratio as though it were the
  // energy-norm one understates how much error delta u misses. delta u is therefore a good
  // but NOT near-exact proxy for u - u_p, so the high scores are only partly expected.
  //
  // This benchmark establishes an UPPER BOUND on what any hierarchical scheme can achieve
  // here; it is NOT itself a practical estimator, because the p+1 solve costs more than the
  // p solve it is estimating. Its value is that it proves the DIRECTION is right, which
  // augmented recovery did not: marking built from a hierarchical correction does not
  // under-select.
  //
  // CORRECTION to an earlier reading: I attributed that success partly to the correction
  // adjusting the SINGULAR AMPLITUDE. The frozen-enrichment ablation disproves it --
  // holding the singular coefficients fixed changes the energy by under 0.25 percent and
  // leaves every marking metric identical to six digits. On this problem the p increment
  // barely readjusts the singular amplitude, so the benefit comes from resolving the smooth
  // remainder, not from the singular response. A strongly hybridized singular mode could
  // behave differently.
  //
  // Both gates below are therefore comparisons against production, not claims of
  // practicality.
  CHECK(mean_rank_h > mean_rank_p);
  CHECK(mean_conc_h > mean_conc_p);

  // ABLATION 2, corrected. Means over the three passes:
  //
  //   metric    production   new-DOFs-only   full hierarchical
  //   RANK        0.679          0.763            0.99992
  //   EXTEND      0.595          0.704            1.00000
  //   CONC        0.444          0.749            1.199
  //   energy         --          0.475 of full        1.0
  //
  // TWO EARLIER CLAIMS WITHDRAWN, both artifacts of the buggy first implementation:
  //   * "new-only overshoots the full energy by 4 to 5.4x because it is un-relaxed against
  //   the
  //     coarse space" -- FALSE. The corrected ratio is 0.475, respecting the SPD bound.
  //   * "new-only reproduces production exactly" -- FALSE. It sits strictly BETWEEN
  //   production
  //     and the full correction on all three metrics.
  //
  // What the corrected ablation does show: the coarse response matters, but not
  // all-or-nothing. Dropping it entirely still recovers about half the correction energy
  // and roughly 60 percent of the RANK gap and 68 percent of the concentration gap. So a
  // new-DOF-only patch scheme would beat production while falling well short of the bound.
  // Retaining the coarse response remains the design requirement -- because it buys the
  // remaining half, not because omitting it is catastrophic.
  const double mean_rank_nd = mean(rank_new_only);
  const double mean_conc_nd = mean(concentration_new_only);
  const double mean_energy_ratio = mean(energy_ratio_new_only);
  CAPTURE(mean_rank_nd, mean_conc_nd, mean_energy_ratio);
  CHECK(mean_energy_ratio < 1.0);
  CHECK(mean_energy_ratio > 0.2);
  CHECK(mean_rank_nd > mean_rank_p);
  CHECK(mean_rank_nd < mean_rank_h);
  CHECK(mean_conc_nd > mean_conc_p);
  CHECK(mean_conc_nd < mean_conc_h);
  // Pin the near-oracle quality so a localized approximation can be measured against it as
  // a fraction of this bound rather than in isolation.
  CHECK(mean_rank_h > 0.95);
  CHECK(mean_extend_h > 0.95);
  CHECK(mean_conc_h > 0.8);
}

TEST_CASE("Manufactured L-shape sparse local candidate is independent and support-bounded",
          "[singularmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2;
  auto mesh = MmsLShapeMesh(2);
  const auto coarse = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);

  // Build before the global benchmark: the candidate sees a coarse solution and an
  // assembly-only p+1 system, never a fine solution.
  const auto sparse = MmsBuildSparseLocalCandidate(mesh, order, 20, coarse);
  CAPTURE(sparse.energy, sparse.maximum_patch_residual, sparse.maximum_patch_condition,
          sparse.maximum_owner_elements, sparse.maximum_actual_support_elements,
          sparse.maximum_actual_to_owner_support, sparse.maximum_element_overlap,
          sparse.mean_element_overlap, sparse.maximum_patch_dimension,
          sparse.sparse_basis_nonzeros, sparse.dense_basis_entries_avoided,
          sparse.touched_matrix_nonzeros, sparse.sum_patch_dimension_squared,
          sparse.sum_patch_dimension_cubed);
  CHECK(sparse.fine_global_solves == 0);
  CHECK(sparse.fine_global_matrix_assemblies == 0);
  CHECK(sparse.maximum_patch_residual < 1.0e-10);
  CHECK(sparse.maximum_patch_condition < 1.0e5);
  CHECK(sparse.maximum_owner_elements <= 2);
  CHECK(sparse.maximum_actual_support_elements > sparse.maximum_owner_elements);
  CHECK(sparse.maximum_actual_support_elements < 20);
  CHECK(sparse.maximum_element_overlap < 100);
  CHECK(sparse.maximum_patch_dimension < 20);
  CHECK(sparse.dense_basis_entries_avoided > 5 * sparse.sparse_basis_nonzeros);
  CHECK(sparse.touched_matrix_nonzeros > 0);
  CHECK_THAT(sparse.energy, WithinRel(sparse.work, 1.0e-10));
  CHECK_THAT(std::accumulate(sparse.indicator.begin(), sparse.indicator.end(), 0.0),
             WithinRel(sparse.energy, 1.0e-7));

  // Tiny-oracle congruence: independently scatter the element-local p+1 assembly and
  // compare every matrix/load/constraint entry with the established global assembly path.
  const auto local_system = MmsAssembleLocalFineSystem(mesh, order + 1, 20);
  const auto global_system = MmsSolveOnMesh(
      mesh, order + 1, 20, true, MmsCornerQuadrature::GRADED, 0.0, nullptr, false);
  const int combined_size = local_system.standard_size + local_system.enrichment_size;
  REQUIRE(global_system.combined_matrix->Height() == combined_size);
  mfem::DenseMatrix local_matrix(combined_size), global_matrix(combined_size);
  local_matrix = 0.0;
  global_matrix = 0.0;
  mfem::Vector local_rhs(combined_size);
  local_rhs = 0.0;
  for (const auto &data : local_system.elements)
  {
    for (int row = 0; row < static_cast<int>(data.dofs.size()); row++)
    {
      local_rhs(data.dofs[row]) += data.rhs(row);
      for (int column = 0; column < static_cast<int>(data.dofs.size()); column++)
      {
        local_matrix(data.dofs[row], data.dofs[column]) += data.matrix(row, column);
      }
    }
  }
  for (int row = 0; row < combined_size; row++)
  {
    for (int entry = global_system.combined_matrix->GetI()[row];
         entry < global_system.combined_matrix->GetI()[row + 1]; entry++)
    {
      global_matrix(row, global_system.combined_matrix->GetJ()[entry]) =
          global_system.combined_matrix->GetData()[entry];
    }
  }
  double maximum_matrix_difference = 0.0, maximum_rhs_difference = 0.0;
  for (int row = 0; row < combined_size; row++)
  {
    REQUIRE(local_system.essential[row] == global_system.essential_mask[row]);
    maximum_rhs_difference = std::max(
        maximum_rhs_difference, std::abs(local_rhs(row) - global_system.combined_rhs(row)));
    for (int column = 0; column < combined_size; column++)
    {
      maximum_matrix_difference =
          std::max(maximum_matrix_difference,
                   std::abs(local_matrix(row, column) - global_matrix(row, column)));
    }
  }
  CAPTURE(maximum_matrix_difference, maximum_rhs_difference);
  CHECK(maximum_matrix_difference < 1.0e-11);
  CHECK(maximum_rhs_difference < 1.0e-11);

  const auto global = MmsTwoLevelBenchmark(mesh, order, 20);
  const auto comparison =
      MmsCompareToOracle(sparse.indicator, coarse.element_error, global.indicator, 0.5);
  CAPTURE(comparison.rank, comparison.extend, comparison.concentration,
          comparison.correction_capture);
  CHECK(comparison.rank > 0.9);
  CHECK(comparison.extend > 0.9);
  CHECK(comparison.correction_capture > 0.4);

  // Rebuilding after the fine oracle solve must be bitwise-independent of that solution.
  const auto rebuilt = MmsBuildSparseLocalCandidate(mesh, order, 20, coarse);
  REQUIRE(rebuilt.indicator.size() == sparse.indicator.size());
  CHECK_THAT(rebuilt.energy, WithinAbs(sparse.energy, 1.0e-14));
  for (std::size_t element = 0; element < sparse.indicator.size(); element++)
  {
    CHECK_THAT(rebuilt.indicator[element], WithinAbs(sparse.indicator[element], 1.0e-14));
  }
}

TEST_CASE("Manufactured L-shape sparse local work scales at fixed order",
          "[singularmms][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2;
  std::vector<int> elements, maximum_support, maximum_overlap, maximum_dimension;
  std::vector<double> basis_per_element, matrix_touches_per_element, cubic_work_per_element;
  for (int n : {2, 3, 4, 6})
  {
    auto mesh = MmsLShapeMesh(n);
    const auto coarse = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);
    const auto local = MmsBuildSparseLocalCandidate(mesh, order, 20, coarse);
    elements.push_back(mesh.GetNE());
    maximum_support.push_back(local.maximum_actual_support_elements);
    maximum_overlap.push_back(local.maximum_element_overlap);
    maximum_dimension.push_back(local.maximum_patch_dimension);
    basis_per_element.push_back(static_cast<double>(local.sparse_basis_nonzeros) /
                                mesh.GetNE());
    matrix_touches_per_element.push_back(
        static_cast<double>(local.touched_matrix_nonzeros) / mesh.GetNE());
    cubic_work_per_element.push_back(static_cast<double>(local.sum_patch_dimension_cubed) /
                                     mesh.GetNE());
    CHECK(local.fine_global_solves == 0);
    CHECK(local.fine_global_matrix_assemblies == 0);
  }
  CAPTURE(elements, maximum_support, maximum_overlap, maximum_dimension, basis_per_element,
          matrix_touches_per_element, cubic_work_per_element);
  CHECK(*std::max_element(maximum_support.begin(), maximum_support.end()) <= 20);
  CHECK(*std::max_element(maximum_overlap.begin(), maximum_overlap.end()) <= 45);
  CHECK(*std::max_element(maximum_dimension.begin(), maximum_dimension.end()) <= 12);
  CHECK(*std::max_element(basis_per_element.begin(), basis_per_element.end()) < 200.0);
  CHECK(*std::max_element(matrix_touches_per_element.begin(),
                          matrix_touches_per_element.end()) < 1600.0);
  CHECK(*std::max_element(cubic_work_per_element.begin(), cubic_work_per_element.end()) <
        1400.0);
}

TEST_CASE("Manufactured extruded L-shape 3D edge-face patches capture the hierarchy",
          "[.singularmms3d][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  auto mesh = MmsExtrudedLShapeMesh(1, 2);
  const auto hierarchical = Mms3DHierarchicalBenchmark(mesh, 2);
  CAPTURE(hierarchical.global_energy, hierarchical.local_energy,
          hierarchical.captured_fraction, hierarchical.global_residual,
          hierarchical.maximum_patch_residual, hierarchical.maximum_patch_condition,
          hierarchical.overlap_error, hierarchical.edge_patches, hierarchical.face_patches,
          hierarchical.interior_patches, hierarchical.owned_modes,
          hierarchical.maximum_patch_dimension);
  CHECK(hierarchical.global_residual < 1.0e-10);
  CHECK(hierarchical.maximum_patch_residual < 1.0e-10);
  CHECK(hierarchical.maximum_patch_condition < 1.0e7);
  CHECK(hierarchical.overlap_error < 1.0e-14);
  CHECK(hierarchical.edge_patches > 0);
  CHECK(hierarchical.face_patches > 0);
  CHECK(hierarchical.interior_patches == 0);
  CHECK(hierarchical.owned_modes > 0);
  CHECK(hierarchical.maximum_patch_dimension < 40);
  CHECK(hierarchical.local_energy >= 0.0);
  CHECK(hierarchical.local_energy <= hierarchical.global_energy * (1.0 + 1.0e-8));
  CHECK(hierarchical.captured_fraction > 0.5);
  const double global_integrated = std::accumulate(
      hierarchical.indicator_global.begin(), hierarchical.indicator_global.end(), 0.0);
  const double local_integrated = std::accumulate(hierarchical.indicator_local.begin(),
                                                  hierarchical.indicator_local.end(), 0.0);
  CHECK_THAT(global_integrated, WithinRel(hierarchical.global_energy, 2.0e-2));
  CHECK_THAT(local_integrated, WithinRel(hierarchical.local_energy, 2.0e-2));

  const auto local =
      MmsCompareToOracle(hierarchical.indicator_local, hierarchical.element_error,
                         hierarchical.indicator_global, 0.5);
  const auto global =
      MmsCompareToOracle(hierarchical.indicator_global, hierarchical.element_error,
                         hierarchical.indicator_global, 0.5);
  CAPTURE(local.rank, local.extend, local.concentration, local.correction_capture,
          global.rank, global.extend, global.concentration, global.correction_capture);
  CHECK(local.rank > 0.7);
  CHECK(local.extend > 0.7);
  CHECK(local.correction_capture > 0.3);
  CHECK(local.rank <= global.rank + 1.0e-12);
}

TEST_CASE("Manufactured extruded L-shape 3D local hierarchical AMR beats uniform",
          "[.singularmms3d][singularelements][Serial]")
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  // Use p=1 -> p=2 in the repeated loop to keep this serial qualification bounded; the
  // separate fixed-mesh p=2 -> p=3 test above explicitly exercises 3D face ownership.
  constexpr int order = 1, passes = 3;
  constexpr double theta = 0.5;
  auto mesh = MmsExtrudedLShapeMesh(1, 2);
  std::vector<long long> dofs;
  std::vector<double> errors, capture;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto hierarchical = Mms3DHierarchicalBenchmark(mesh, order);
    dofs.push_back(hierarchical.coarse_dofs);
    errors.push_back(hierarchical.coarse_error);
    capture.push_back(hierarchical.captured_fraction);
    if (pass + 1 < passes)
    {
      mfem::Array<mfem::Refinement> refinements;
      for (std::size_t element : MmsDorflerMark(hierarchical.indicator_local, theta))
      {
        refinements.Append(mfem::Refinement(static_cast<int>(element)));
      }
      REQUIRE(refinements.Size() > 0);
      mesh.GeneralRefinement(refinements, -1, 1);
      REQUIRE(mesh.Conforming());
    }
  }

  std::vector<std::pair<long long, double>> uniform;
  for (int n : {1, 2, 3})
  {
    auto uniform_mesh = MmsExtrudedLShapeMesh(n, 2 * n);
    const auto solve = MmsSolve3D(uniform_mesh, order, true);
    uniform.emplace_back(solve.dofs, solve.error);
  }
  const auto envelope_at = [&uniform](long long at_dofs)
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
    const double envelope = envelope_at(dofs[pass]);
    REQUIRE(std::isfinite(envelope));
    ratios.push_back(errors[pass] / envelope);
    if (pass > 0)
    {
      CHECK(errors[pass] < errors[pass - 1]);
    }
  }
  CAPTURE(dofs, errors, capture, uniform, ratios);
  CHECK(ratios.back() < 1.0);
  CHECK(*std::min_element(capture.begin(), capture.end()) > 0.5);
}

TEST_CASE("Manufactured L-shape local hierarchical patches approach the global benchmark",
          "[singularmms][singularelements][Serial]")
{
  // Identical oracle meshes isolate estimator quality from path dependence. The global
  // p+1 difference is the upper benchmark; the candidate consists of overlapping
  // edge/interior corrections with incident coarse response and explicit hierarchical-key
  // PoU weighting.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 3;
  constexpr double theta = 0.5;

  auto mesh = MmsLShapeMesh(2);
  std::vector<double> production_rank, local_rank, global_rank, production_extend,
      local_extend, global_extend, production_concentration, local_concentration,
      global_concentration, production_correction_capture, local_correction_capture,
      global_correction_capture, aggregate_captured_fraction;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto solve = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);
    const auto hierarchical = MmsTwoLevelBenchmark(mesh, order, 20);
    REQUIRE(hierarchical.indicator_local_patches.size() == solve.element_error.size());

    const auto production = MmsCompareToOracle(
        solve.indicator_standard, solve.element_error, hierarchical.indicator, theta);
    const auto local =
        MmsCompareToOracle(hierarchical.indicator_local_patches, solve.element_error,
                           hierarchical.indicator, theta);
    const auto global = MmsCompareToOracle(hierarchical.indicator, solve.element_error,
                                           hierarchical.indicator, theta);
    CAPTURE(pass, solve.elements, production.rank, local.rank, global.rank,
            production.extend, local.extend, global.extend, production.concentration,
            local.concentration, global.concentration, production.correction_capture,
            local.correction_capture, global.correction_capture,
            hierarchical.local_patch_captured_fraction, hierarchical.maximum_patch_residual,
            hierarchical.maximum_patch_energy_identity,
            hierarchical.maximum_patch_energy_ratio, hierarchical.overlap_partition_error,
            hierarchical.minimum_overlap, hierarchical.maximum_overlap,
            hierarchical.maximum_patch_dimension, hierarchical.maximum_patch_elements,
            hierarchical.maximum_patch_condition, hierarchical.uncovered_coarse_guests,
            hierarchical.uncovered_singular_guests);

    // Algebraic patch invariants remain mandatory on every adapted oracle mesh.
    CHECK(hierarchical.global_correction_residual < 1.0e-10);
    CHECK(hierarchical.global_correction_energy_identity < 1.0e-10);
    CHECK(hierarchical.maximum_patch_residual < 1.0e-10);
    CHECK(hierarchical.maximum_patch_energy_identity < 1.0e-10);
    CHECK(hierarchical.maximum_patch_energy_ratio <= 1.0 + 1.0e-8);
    CHECK(hierarchical.overlap_partition_error < 1.0e-14);
    CHECK(hierarchical.edge_patches > 0);
    CHECK(hierarchical.interior_patches == mesh.GetNE());
    CHECK(hierarchical.maximum_patch_dimension <= 20);
    CHECK(hierarchical.maximum_patch_elements <= 2);
    CHECK(hierarchical.maximum_patch_condition < 1.0e5);
    CHECK(hierarchical.uncovered_coarse_guests == 0);
    CHECK(hierarchical.uncovered_singular_guests == 0);
    CHECK(hierarchical.local_patch_captured_fraction >= 0.0);
    CHECK(hierarchical.local_patch_captured_fraction <= 1.0 + 1.0e-8);
    CHECK(hierarchical.local_indicator_energy_identity < 1.0e-7);

    for (const auto *comparison : {&production, &local, &global})
    {
      CHECK(comparison->rank >= 0.0);
      CHECK(comparison->rank <= 1.0 + 1.0e-12);
      CHECK(comparison->extend >= 0.0);
      CHECK(comparison->extend <= 1.0 + 1.0e-12);
      CHECK(comparison->concentration > 0.0);
      // Dörfler controls candidate-indicator mass, not exact-error mass. True capture is a
      // quality metric and need only remain a valid fraction here.
      CHECK(comparison->true_capture >= 0.0);
      CHECK(comparison->true_capture <= 1.0 + 1.0e-12);
      CHECK(comparison->correction_capture >= 0.0);
      CHECK(comparison->correction_capture <= 1.0 + 1.0e-12);
    }

    production_rank.push_back(production.rank);
    local_rank.push_back(local.rank);
    global_rank.push_back(global.rank);
    production_extend.push_back(production.extend);
    local_extend.push_back(local.extend);
    global_extend.push_back(global.extend);
    production_concentration.push_back(production.concentration);
    local_concentration.push_back(local.concentration);
    global_concentration.push_back(global.concentration);
    production_correction_capture.push_back(production.correction_capture);
    local_correction_capture.push_back(local.correction_capture);
    global_correction_capture.push_back(global.correction_capture);
    aggregate_captured_fraction.push_back(hierarchical.local_patch_captured_fraction);

    if (pass + 1 < passes)
    {
      mfem::Array<mfem::Refinement> refinements;
      for (std::size_t element : MmsDorflerMark(solve.element_error, theta))
      {
        refinements.Append(mfem::Refinement(static_cast<int>(element)));
      }
      mesh.GeneralRefinement(refinements);
      REQUIRE(mesh.Conforming());
    }
  }

  const auto mean = [](const std::vector<double> &value)
  { return std::accumulate(value.begin(), value.end(), 0.0) / value.size(); };
  CAPTURE(production_rank, local_rank, global_rank, production_extend, local_extend,
          global_extend, production_concentration, local_concentration,
          global_concentration, production_correction_capture, local_correction_capture,
          global_correction_capture, aggregate_captured_fraction);

  // MEASURED on the three oracle meshes at theta=0.5 after hierarchical-key PoU:
  //
  //   metric                  production    local patches    global hierarchical
  //   RANK mean                  0.679           0.987              0.9999
  //   EXTEND mean                0.595           0.980              1.0000
  //   concentration mean         0.444           1.174              1.199
  //   marked correction share    0.149           0.508              0.532
  //
  // The PoU-assembled local correction itself captures 0.886 to 0.923 (mean 0.900) of the
  // global correction energy. Candidate-quality gates are comparisons, not assumptions of
  // equality with the global solve.
  CHECK(mean(local_rank) > 0.95);
  CHECK(mean(local_extend) > 0.95);
  CHECK(mean(local_concentration) > mean(production_concentration));
  CHECK(mean(local_correction_capture) > 3.0 * mean(production_correction_capture));
  CHECK(mean(local_rank) <= mean(global_rank) + 1.0e-12);
  CHECK(mean(local_extend) <= mean(global_extend) + 1.0e-12);
  CHECK(mean(local_concentration) < mean(global_concentration));
  CHECK(mean(aggregate_captured_fraction) > 0.85);
}

TEST_CASE("Manufactured L-shape local hierarchical AMR beats the dense uniform envelope",
          "[singularmms][singularelements][Serial]")
{
  // Qualification uses the existing production theta without tuning. The candidate performs
  // no global p+1 solve or global p+1 matrix assembly, stores sparse support-sized basis
  // columns, and assembles every restricted operator from complete element support unions.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 5;
  constexpr double theta = 0.5;

  std::vector<std::pair<long long, double>> uniform;
  for (int n : {2, 3, 4, 5, 6, 8})
  {
    const auto report = MmsSolve(n, order, 20, true);
    uniform.emplace_back(report.dofs, report.error_direct);
  }
  const auto envelope_at = [&uniform](long long at_dofs)
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

  std::vector<int> maximum_support, maximum_overlap, maximum_dimension;
  std::vector<long long> sparse_nonzeros, avoided_dense_entries, touched_matrix_nonzeros;
  const auto history = MmsRunAmrWithIndicator(
      true, order, passes, theta,
      [&](mfem::Mesh &mesh, const MmsSolveReport &report)
      {
        return MmsSelectSingularOnlyIndicator(
            report,
            [&]()
            {
              const auto local = MmsBuildSparseLocalCandidate(mesh, order, 20, report);
              REQUIRE(local.fine_global_solves == 0);
              REQUIRE(local.fine_global_matrix_assemblies == 0);
              maximum_support.push_back(local.maximum_actual_support_elements);
              maximum_overlap.push_back(local.maximum_element_overlap);
              maximum_dimension.push_back(local.maximum_patch_dimension);
              sparse_nonzeros.push_back(local.sparse_basis_nonzeros);
              avoided_dense_entries.push_back(local.dense_basis_entries_avoided);
              touched_matrix_nonzeros.push_back(local.touched_matrix_nonzeros);
              return local.indicator;
            });
      });
  REQUIRE(history.size() == static_cast<std::size_t>(passes));
  REQUIRE(maximum_support.size() == history.size());

  std::vector<long long> dofs;
  std::vector<double> errors, ratios;
  for (const auto &pass : history)
  {
    dofs.push_back(pass.dofs);
    errors.push_back(pass.error);
    const double envelope = envelope_at(pass.dofs);
    REQUIRE(std::isfinite(envelope));
    ratios.push_back(pass.error / envelope);
  }
  CAPTURE(dofs, errors, ratios, maximum_support, maximum_overlap, maximum_dimension,
          sparse_nonzeros, avoided_dense_entries, touched_matrix_nonzeros, uniform);
  for (std::size_t pass = 0; pass + 1 < errors.size(); pass++)
  {
    CHECK(errors[pass + 1] < errors[pass]);
  }

  // MEASURED final point: 673 solution DOFs, e^2 = 2.974e-4, 0.408x the dense uniform
  // envelope. Production previously reached 2.178e-3 at 345 DOFs (1.54x uniform), while
  // the exact-error oracle reached 4.322e-4 at 576 DOFs (0.51x uniform). The local sequence
  // therefore reaches the oracle-quality envelope without retuning theta. Its cost-matched
  // ratios after the common initial mesh are 1.24, 0.848, 0.689, and 0.408; it beats dense
  // uniform from the third solve onward rather than only at a cherry-picked final pass.
  CHECK(ratios.back() < 0.5);
  CHECK(ratios.back() < ratios.front());
  CHECK(ratios[2] < 1.0);
  CHECK(ratios[3] < 1.0);
  CHECK(*std::max_element(maximum_support.begin(), maximum_support.end()) < 25);
  CHECK(*std::max_element(maximum_overlap.begin(), maximum_overlap.end()) < 100);
  CHECK(*std::max_element(maximum_dimension.begin(), maximum_dimension.end()) < 20);
  for (std::size_t pass = 0; pass < history.size(); pass++)
  {
    CHECK(avoided_dense_entries[pass] > 5 * sparse_nonzeros[pass]);
    CHECK(touched_matrix_nonzeros[pass] > 0);
  }
}

TEST_CASE("Manufactured L-shape free-coefficient augmentation barely changes the marking",
          "[singularmms][singularelements][Serial]")
{
  // Step 1 of the augmented-recovery evaluation, using the RESIDUAL-relative fraction rho
  // rather than the mode-relative one.
  //
  //   rho = |<r, s_perp>|^2 / (||r||^2 ||s_perp||^2),   r = E - P_RT E,  s_perp = s - P_RT
  //   s
  //
  // rho is the fraction of the squared RECOVERY RESIDUAL that one free singular mode can
  // remove, which is what the estimator actually sees. It is NOT ||s_perp||^2/||s||^2: a
  // mode holding well under one percent of the total field energy can still be strongly
  // aligned with the residual and remove most of it, so the earlier mode-relative reading
  // was the wrong basis for a go/no-go decision.
  //
  // The diagnostics are recomputed on the SAME oracle mesh sequence as the production
  // comparison, so no additional AMR loop is needed to decide whether this candidate is
  // worth pursuing.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int order = 2, passes = 4;
  constexpr double theta = 0.5;

  auto mesh = MmsLShapeMesh(2);
  std::vector<double> rhos, residual_reductions, rank_production, rank_augmented,
      extend_production, extend_augmented, concentration_production,
      concentration_augmented;
  for (int pass = 0; pass < passes; pass++)
  {
    const auto report = MmsSolveOnMesh(mesh, order, 20, true, MmsCornerQuadrature::GRADED);
    REQUIRE(report.indicator_augmented.size() == report.element_error.size());
    rhos.push_back(report.augmentation_rho);
    residual_reductions.push_back(1.0 - report.residual_energy_after /
                                            report.residual_energy_before);

    double error_total = 0.0;
    for (double e : report.element_error)
    {
      error_total += e;
    }
    const auto capture = [&](const std::vector<std::size_t> &marked)
    {
      double captured = 0.0;
      for (std::size_t element : marked)
      {
        captured += report.element_error[element];
      }
      return captured / error_total;
    };
    const auto top_k = [](const std::vector<double> &value, std::size_t count)
    {
      std::vector<std::size_t> index(value.size());
      std::iota(index.begin(), index.end(), std::size_t{0});
      count = std::min(count, value.size());
      std::partial_sort(index.begin(), index.begin() + count, index.end(),
                        [&value](std::size_t a, std::size_t b)
                        { return value[a] > value[b]; });
      return std::vector<std::size_t>(index.begin(), index.begin() + count);
    };

    const auto oracle = MmsDorflerMark(report.element_error, theta);
    const double oracle_capture = capture(oracle);
    // Same three metrics as the production diagnosis, for both indicators.
    const auto measure = [&](const std::vector<double> &indicator, double &rank,
                             double &extend, double &concentration)
    {
      const auto marked = MmsDorflerMark(indicator, theta);
      REQUIRE(!marked.empty());
      rank = capture(marked) / capture(top_k(report.element_error, marked.size()));
      extend = capture(top_k(indicator, oracle.size())) / oracle_capture;
      concentration =
          static_cast<double>(marked.size()) / static_cast<double>(oracle.size());
    };
    double rank_p, extend_p, concentration_p, rank_a, extend_a, concentration_a;
    measure(report.indicator_standard, rank_p, extend_p, concentration_p);
    measure(report.indicator_augmented, rank_a, extend_a, concentration_a);
    rank_production.push_back(rank_p);
    rank_augmented.push_back(rank_a);
    extend_production.push_back(extend_p);
    extend_augmented.push_back(extend_a);
    concentration_production.push_back(concentration_p);
    concentration_augmented.push_back(concentration_a);

    CAPTURE(pass, report.elements, report.augmentation_rho, residual_reductions.back(),
            rank_p, rank_a, extend_p, extend_a, concentration_p, concentration_a);

    // The augmented recovery must be a genuine improvement in the LEAST-SQUARES sense,
    // since it minimizes over a strictly larger space. If this failed, the fit would be
    // wrong.
    CHECK(report.residual_energy_after <= report.residual_energy_before * (1.0 + 1.0e-12));

    // The measured energy reduction must equal rho exactly, since the free-coefficient
    // optimum removes |<r,s_perp>|^2 / ||s_perp||^2 by construction. This is an independent
    // check on the algebra: rho is computed from inner products, the reduction from
    // elementwise quadrature of two different recovered fields.
    CHECK_THAT(residual_reductions.back(), WithinAbs(report.augmentation_rho, 1.0e-9));

    // The augmented indicator must genuinely DIFFER elementwise, or the unchanged marking
    // below would be a wiring bug rather than a finding.
    double worst_difference = 0.0;
    for (std::size_t element = 0; element < report.element_error.size(); element++)
    {
      worst_difference =
          std::max(worst_difference, std::abs(report.indicator_augmented[element] -
                                              report.indicator_standard[element]));
    }
    CAPTURE(worst_difference);
    CHECK(worst_difference > 0.0);

    if (pass + 1 == passes)
    {
      break;
    }
    mfem::Array<mfem::Refinement> refinements;
    for (std::size_t element : oracle)
    {
      refinements.Append(mfem::Refinement(static_cast<int>(element)));
    }
    mesh.GeneralRefinement(refinements);
  }

  const auto mean = [](const std::vector<double> &v)
  { return std::accumulate(v.begin(), v.end(), 0.0) / v.size(); };
  CAPTURE(rhos, residual_reductions, rank_production, rank_augmented, extend_production,
          extend_augmented, concentration_production, concentration_augmented);
  const double mean_rho = mean(rhos);
  const double mean_rank_gain = mean(rank_augmented) - mean(rank_production);
  const double mean_extend_gain = mean(extend_augmented) - mean(extend_production);
  const double mean_concentration_gain =
      mean(concentration_augmented) - mean(concentration_production);
  CAPTURE(mean_rho, mean_rank_gain, mean_extend_gain, mean_concentration_gain);

  // rho must be a valid squared cosine.
  for (double rho : rhos)
  {
    CHECK(rho >= 0.0);
    CHECK(rho <= 1.0 + 1.0e-12);
  }

  // MEASURED, order 2, theta = 0.5, on the oracle mesh sequence:
  //
  //   pass  ne    rho     residual reduction   RANK p/a        EXTEND p/a      k_p/k_o p/a
  //     0    24   0.0817       0.0817        0.670 / 0.670   0.587 / 0.587   0.600 / 0.600
  //     1    57   0.0788       0.0788        0.523 / 0.523   0.568 / 0.568   0.500 / 0.500
  //     2    90   0.0887       0.0887        0.843 / 0.843   0.629 / 0.629   0.231 / 0.231
  //     3   169   0.0847       0.0847        0.823 / 0.823   0.818 / 0.838   0.222 / 0.222
  //   means      0.0835                      gain -1e-16     gain +0.0051    gain 0.000
  //
  // rho equals the measured residual reduction to 1e-12 in every pass, which independently
  // confirms the algebra: rho is formed from inner products while the reduction is
  // quadrature of two separately recovered fields.
  //
  // VERDICT: NO-GO for augmented recovery. One free singular mode removes about 8 percent
  // of the recovery residual -- far more than the 0.3 to 3 percent the mode-relative
  // fraction suggested, so the residual-relative measure was indeed the right one -- but it
  // does not REORDER the elements. RANK is unchanged to 1e-16, concentration is unchanged
  // exactly, and EXTEND moves by +0.005 in one of four passes. The augmented indicator does
  // differ elementwise (asserted above, so this is not a wiring artifact); the differences
  // are simply near-uniform across the elements that compete for marking, so every Dorfler
  // set comes out the same.
  //
  // Consequence: the marking defects diagnosed earlier, imperfect ranking (RANK 0.715) and
  // excessive concentration (k_p/k_o falling to 0.222), are NOT caused by the recovery
  // space missing the singular direction. Augmented recovery work stops here, and the
  // hierarchical local-correction estimator becomes the primary candidate rather than a
  // comparison.
  //
  // Gates are on the MARKING diagnostics, not on rho: a large rho that does not reorder
  // elements is useless for AMR, which is precisely what happened.
  CHECK(std::abs(mean_rank_gain) < 0.01);
  CHECK(std::abs(mean_extend_gain) < 0.05);
  CHECK(std::abs(mean_concentration_gain) < 0.01);

  // rho is substantial, so the no-go is specifically "does not reorder", not "mode is
  // irrelevant to the residual". Pinned so the distinction survives.
  CHECK(mean_rho > 0.02);
}

TEST_CASE("Manufactured L-shape singular flux modes are admissible for RT augmentation",
          "[singularmms][singularelements][Serial]")
{
  // FEASIBILITY GATE, before any recovery machinery is built on top of these modes. Each
  // requirement below is one the augmented recovery silently depends on, so each is checked
  // rather than assumed.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int n = GENERATE(2, 4);
  const int order = GENERATE(1, 2);
  CAPTURE(n, order);

  const auto report = MmsAssessFluxAugmentation(n, order, 12);
  CAPTURE(report.extracted_nu, report.worst_normal_jump, report.independent_energy_fraction,
          report.mass_matrix_spd, report.worst_rt_coefficient_error,
          report.singular_coefficient_error);

  // (1) Provenance: the exponent came from the extractor, i.e. from the wedge geometry and
  // materials, and it agrees with the closed-form 3 pi / 2 Dirichlet wedge value. The
  // manufactured solution was not consulted.
  CHECK_THAT(report.extracted_nu, WithinRel(2.0 / 3.0, 1.0e-12));

  // (2) Normal trace is continuous across every interior face.
  CHECK(report.worst_normal_jump < 1.0e-10);

  // (3) The mode is not almost inside RT, so orthogonalizing against RT leaves something,
  // and the augmented block is well conditioned rather than nearly singular.
  //
  // MEASURED, and this is a WARNING about step 2 rather than a clean pass:
  //   n=2 order 1  0.0331      n=4 order 1  0.0150
  //   n=2 order 2  0.0071      n=4 order 2  0.0028
  // RT already represents 96.7 to 99.7 percent of the mode's energy, and the independent
  // fraction DECAYS with both refinement (factor ~2.3) and order (factor ~4.7). So the
  // augmentation adds only a small and shrinking direction: the augmented recovery has
  // limited room to differ from plain RT, and less of it exactly on the finer, higher-order
  // meshes that matter. The gate below confirms well-posedness, but it is a weak prior for
  // the augmented recovery actually fixing the marking, which strengthens the case for
  // treating the hierarchical/residual estimator as a genuine comparison rather than a
  // fallback.
  CHECK(report.mass_matrix_spd);
  CHECK(report.independent_energy_fraction > 1.0e-3);
  CHECK(report.independent_energy_fraction <= 1.0);

  // (4) A field assembled from KNOWN RT and singular coefficients is recovered back to
  // quadrature tolerance. Without this, a later "the recovery works" claim would be
  // untestable. Measured: 1e-13 or better on both.
  CHECK(report.singular_coefficient_error < 1.0e-8);
  CHECK(report.worst_rt_coefficient_error < 1.0e-8);

  // (5) is structural rather than numerical and is enforced by construction:
  // MmsSingularFluxMode provides ONE mode per singular FEATURE, keyed to the extracted
  // exponent, not one per element-local enrichment gradient. An unrestricted per-element
  // set could absorb ordinary polynomial approximation error and make the recovery look
  // good for the wrong reason.
  //
  // (6) conductor-sheet faces must stay selectively broken when this reaches 3D: the normal
  // trace is intentionally DISCONTINUOUS across a zero-thickness conductor, so the
  // continuity check above must exempt exactly those faces and no others. Not exercised
  // here because the 2D reentrant corner has no internal sheet; it is a prerequisite for
  // the 3D step.
}

TEST_CASE("Manufactured L-shape enriched quadrature attribution", "[.mmsdiag][Serial]")
{
  // Diagnostic only. Reports how the stage-2(ii) gap behaves as the quadrature rule is
  // refined, which distinguishes quadrature error from an assembly inconsistency.
  for (const auto scheme : {MmsCornerQuadrature::GRADED, MmsCornerQuadrature::DUFFY})
  {
    for (int quadrature_order : {20, 24, 30, 36})
    {
      const auto report = MmsSolve(2, 2, quadrature_order, true, scheme);
      const std::string label = scheme == MmsCornerQuadrature::DUFFY ? "duffy" : "graded";
      CAPTURE(label, quadrature_order);
      WARN("a-l=" << (report.a_u_uh - report.l_uh) << " direct=" << report.error_direct
                  << " expanded=" << report.error_expanded
                  << " galerkin=" << report.error_galerkin);
    }
  }
}

TEST_CASE("Manufactured L-shape enrichment beats the standard space on identical meshes",
          "[singularmms][singularelements][Serial]")
{
  // Stage 5b. With both spaces certified by the same algebra gates, compare them on
  // IDENTICAL meshes. The enrichment adds a handful of DOFs, so any error reduction is
  // attributable to the singular basis resolving r^(2/3) rather than to refinement.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  CAPTURE(order);

  std::vector<double> standard_errors, enriched_errors;
  std::vector<long long> standard_dofs, enriched_dofs;
  for (int n : {2, 4, 8})
  {
    const auto standard = MmsSolve(n, order, 20, false);
    const auto enriched = MmsSolve(n, order, 20, true);
    CAPTURE(n, standard.dofs, standard.error_direct, enriched.dofs,
            enriched.enrichment_dofs, enriched.error_direct);

    // Identical mesh, and the enriched space is the standard space plus a few functions.
    CHECK(standard.elements == enriched.elements);
    CHECK(enriched.dofs == standard.dofs + enriched.enrichment_dofs);

    // The enriched space contains the standard one, so the Galerkin solution can only
    // improve. This is a mathematical guarantee, not a heuristic, and therefore a real
    // check on the implementation.
    CHECK(enriched.error_direct < standard.error_direct);

    standard_errors.push_back(standard.error_direct);
    enriched_errors.push_back(enriched.error_direct);
    standard_dofs.push_back(standard.dofs);
    enriched_dofs.push_back(enriched.dofs);
  }
  CAPTURE(standard_dofs, standard_errors, enriched_dofs, enriched_errors);

  // Both sequences must converge.
  for (std::size_t i = 0; i + 1 < standard_errors.size(); i++)
  {
    CHECK(standard_errors[i + 1] < standard_errors[i]);
    CHECK(enriched_errors[i + 1] < enriched_errors[i]);
  }
}

TEST_CASE("Manufactured L-shape enrichment is near-optimal in its singular exponent",
          "[singularmms][singularelements][Serial]")
{
  // Stage 5c, the control on stage 5b. Six extra functions can only help, whatever their
  // exponent, so "enriched beats standard" alone does not establish that the SINGULAR basis
  // is doing the work. Overriding nu isolates that.
  //
  // MEASURED, and worth recording because it is not the naive expectation: most of the
  // stage-5b gain is NOT exponent-specific. A deliberately wrong nu still recovers 92 to
  // 96 percent of it, because the dominant effect is simply having extra functions
  // concentrated at the corner. What the exponent controls is the remaining few percent,
  // and that part is real, grows with polynomial order, and is strongly asymmetric: at
  // order 2 / n = 8 the error relative to nu = 2/3 is 1.192 at nu = 0.40 and 1.084 at
  // nu = 0.85, but flat to within 0.3 percent across nu in [0.60, 0.68].
  //
  // So the honest claim is near-optimality within a shallow minimum, not a large penalty
  // for getting nu wrong. Asserting the latter would be asserting something false.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  constexpr int n = 8, order = 2;

  const auto standard = MmsSolve(n, order, 20, false);
  const auto correct = MmsSolve(n, order, 20, true);
  REQUIRE(correct.error_direct < standard.error_direct);

  // Exponents far from 2/3 on either side must be measurably worse. Being two-sided is what
  // makes this a real check: a monotone trend would instead just mean "larger nu is
  // better".
  for (const auto &[nu, penalty] :
       std::vector<std::pair<double, double>>{{0.40, 1.10}, {0.85, 1.05}})
  {
    const auto wrong = MmsSolveWrongExponent(n, order, nu);
    CAPTURE(nu, penalty, standard.error_direct, correct.error_direct, wrong.error_direct);
    // Same mesh and same DOF count, so this is purely the exponent's effect.
    CHECK(wrong.dofs == correct.dofs);
    // Still a superset of the standard space, hence never worse than it.
    CHECK(wrong.error_direct <= standard.error_direct * (1.0 + 1.0e-12));
    CHECK(wrong.error_direct > penalty * correct.error_direct);
  }

  // Near the true exponent the error is flat, and nothing in a neighbourhood beats 2/3 by
  // more than a fraction of a percent. The measured optimum sits near nu = 0.63 rather than
  // exactly 2/3 because u = B s is only asymptotically r^nu; the smooth factor B shifts the
  // best single-exponent fit slightly. That is a property of the manufactured solution, not
  // a defect in the enrichment.
  for (const double nu : {0.60, 0.63, 0.70})
  {
    const auto near = MmsSolveWrongExponent(n, order, nu);
    CAPTURE(nu, correct.error_direct, near.error_direct);
    CHECK(near.error_direct > 0.99 * correct.error_direct);
    CHECK(near.error_direct < 1.02 * correct.error_direct);
  }
}

TEST_CASE("Manufactured L-shape standard-space error decays under refinement",
          "[singularmms][singularelements][Serial]")
{
  // Stage 4. With the algebra checks satisfied, the error must actually decrease under
  // uniform refinement, and the three expressions must remain in agreement at every level.
  // This is what makes the sequence usable as an estimator benchmark.
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int order = GENERATE(1, 2);
  CAPTURE(order);

  std::vector<long long> dofs;
  std::vector<double> errors;
  for (int n : {2, 4, 8})
  {
    const auto report = MmsSolveStandard(n, order);
    CAPTURE(n, report.dofs, report.error_direct, report.error_expanded,
            report.error_galerkin, report.algebraic_residual);
    // Gates from stages 2-3 must hold at every refinement level, not just the coarsest.
    CHECK(report.algebraic_residual < 1.0e-10);
    CHECK_THAT(report.error_expanded,
               WithinAbs(report.error_direct, 1.0e-5 * report.error_direct));
    CHECK_THAT(report.error_galerkin,
               WithinAbs(report.error_direct, 1.0e-5 * report.error_direct));
    dofs.push_back(report.dofs);
    errors.push_back(report.error_direct);
  }

  // Monotone decrease, and a genuine rate rather than stagnation.
  REQUIRE(errors.size() == 3);
  CAPTURE(dofs, errors);
  CHECK(errors[1] < errors[0]);
  CHECK(errors[2] < errors[1]);
  // Each uniform refinement halves h; e^2 must contract by a clear factor. For P1 on a
  // reentrant corner the energy-norm rate is corner-limited, so require only a factor of
  // two per level rather than the smooth-problem factor of four.
  CHECK(errors[0] / errors[1] > 2.0);
  CHECK(errors[1] / errors[2] > 2.0);
}

}  // namespace palace
