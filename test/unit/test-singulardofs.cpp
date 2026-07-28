// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>

#include "fem/singulardofs.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

mfem::Mesh RectangleSheetMesh(bool reverse_elements = false)
{
  mfem::Mesh mesh(3, 6, 4, 2, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);

  const std::array<std::array<int, 4>, 4> elements{
      std::array<int, 4>{0, 1, 2, 4}, std::array<int, 4>{0, 2, 1, 5},
      std::array<int, 4>{1, 3, 2, 4}, std::array<int, 4>{1, 2, 3, 5}};
  if (reverse_elements)
  {
    for (auto element = elements.rbegin(); element != elements.rend(); ++element)
    {
      mesh.AddTet(element->data(), 1);
    }
  }
  else
  {
    for (const auto &element : elements)
    {
      mesh.AddTet(element.data(), 1);
    }
  }
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, 8);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh InternalLineTipMesh(bool reverse_elements = false)
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
  const std::array<std::array<int, 3>, 8> elements{
      std::array<int, 3>{0, 1, 4}, std::array<int, 3>{0, 4, 3}, std::array<int, 3>{3, 4, 7},
      std::array<int, 3>{3, 7, 6}, std::array<int, 3>{1, 2, 5}, std::array<int, 3>{1, 5, 4},
      std::array<int, 3>{4, 5, 8}, std::array<int, 3>{4, 8, 7}};
  if (reverse_elements)
  {
    for (auto element = elements.rbegin(); element != elements.rend(); ++element)
    {
      mesh.AddTriangle(element->data(), 1);
    }
  }
  else
  {
    for (const auto &element : elements)
    {
      mesh.AddTriangle(element.data(), 1);
    }
  }
  mesh.AddBdrSegment(3, 4, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

bool SameBasis(const fem::singular::HigherOrderBasis &a,
               const fem::singular::HigherOrderBasis &b)
{
  return a.family == b.family && a.nodes == b.nodes &&
         a.interpolation_indices == b.interpolation_indices && a.order == b.order &&
         a.nu == b.nu;
}

bool SameBasis(const fem::singular::TriangleBasis &a, const fem::singular::TriangleBasis &b)
{
  return a.family == b.family && a.nodes == b.nodes && a.order == b.order && a.nu == b.nu;
}

fem::singular::Vector3 Cross(const fem::singular::Vector3 &a,
                             const fem::singular::Vector3 &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double Dot(const fem::singular::Vector3 &a, const fem::singular::Vector3 &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

fem::singular::Vector3 MeshEdge(const mfem::Mesh &mesh, int first, int second)
{
  fem::singular::Vector3 result;
  const double *x = mesh.GetVertex(first);
  const double *y = mesh.GetVertex(second);
  for (int d = 0; d < 3; d++)
  {
    result[d] = y[d] - x[d];
  }
  return result;
}

fem::singular::BarycentricGradients PhysicalBarycentricGradients(const mfem::Mesh &mesh,
                                                                 int element)
{
  const int *vertices = mesh.GetElement(element)->GetVertices();
  const auto edge_1 = MeshEdge(mesh, vertices[0], vertices[1]);
  const auto edge_2 = MeshEdge(mesh, vertices[0], vertices[2]);
  const auto edge_3 = MeshEdge(mesh, vertices[0], vertices[3]);
  const double determinant = Dot(edge_1, Cross(edge_2, edge_3));
  REQUIRE(std::abs(determinant) > 0.0);

  fem::singular::BarycentricGradients gradients;
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

fem::singular::BarycentricPoint FacePoint(const mfem::Element &element,
                                          const std::array<int, 3> &face)
{
  constexpr std::array<double, 3> weights{0.217, 0.331, 0.452};
  fem::singular::BarycentricPoint lambda{};
  for (int i = 0; i < 3; i++)
  {
    const int *vertex =
        std::find(element.GetVertices(), element.GetVertices() + 4, face[i]);
    REQUIRE(vertex != element.GetVertices() + 4);
    lambda[vertex - element.GetVertices()] = weights[i];
  }
  return lambda;
}

bool IsSupportedOnFace(const fem::singular::EntityKey &support,
                       const std::array<int, 3> &face)
{
  return support.size <= 3 &&
         std::all_of(support.vertices.begin(), support.vertices.begin() + support.size,
                     [&face](auto vertex)
                     { return std::find(face.begin(), face.end(), vertex) != face.end(); });
}

double EvaluateGradientPotential(const fem::singular::BarycentricPoint &lambda,
                                 const fem::singular::HigherOrderBasis &basis)
{
  using Family = fem::singular::HigherOrderBasisFamily;
  if (basis.family == Family::NODE_GRADIENT)
  {
    return fem::singular::EvaluateHigherOrderNodeGradientPotential(
        lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
        basis.interpolation_indices, basis.order, basis.nu);
  }
  REQUIRE(basis.family == Family::EDGE_GRADIENT);
  return fem::singular::EvaluateHigherOrderEdgeGradientPotential(
      lambda, basis.nodes[0], basis.nodes[1], basis.nodes[2], basis.nodes[3],
      basis.interpolation_indices, basis.order, basis.nu);
}

int InterpolationDenominator(const fem::singular::DofKey &key)
{
  using Family = fem::singular::HigherOrderBasisFamily;
  switch (key.family)
  {
    case Family::NODE_GRADIENT:
      return key.order + 1;
    case Family::NODE_ROTATIONAL:
    case Family::EDGE_GRADIENT:
      return key.order + 2;
    case Family::EDGE_ROTATIONAL:
      return key.order + 3;
  }
  return -1;
}

std::map<std::array<int, 4>,
         std::pair<std::vector<fem::singular::DofKey>, std::vector<fem::singular::DofKey>>>
GlobalElementDofs(const mfem::Mesh &mesh, const fem::singular::DofTopology &topology)
{
  using ElementKeys =
      std::pair<std::vector<fem::singular::DofKey>, std::vector<fem::singular::DofKey>>;
  std::map<std::array<int, 4>, ElementKeys> result;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    std::array<int, 4> element_key{vertices[0], vertices[1], vertices[2], vertices[3]};
    std::sort(element_key.begin(), element_key.end());
    auto &[h1, nd] = result[element_key];
    for (const auto &dof : topology.elements[element].h1)
    {
      h1.push_back(topology.h1_dofs[dof.dof]);
    }
    for (const auto &dof : topology.elements[element].nd)
    {
      nd.push_back(topology.nd_dofs[dof.dof]);
    }
    std::sort(h1.begin(), h1.end());
    std::sort(nd.begin(), nd.end());
  }
  return result;
}

using StableElementKey = std::array<fem::singular::GlobalVertexId, 4>;
using ElementKeys =
    std::pair<std::vector<fem::singular::DofKey>, std::vector<fem::singular::DofKey>>;

std::map<StableElementKey, ElementKeys>
StableElementDofs(const mfem::Mesh &mesh,
                  const std::vector<fem::singular::GlobalVertexId> &vertex_ids,
                  const fem::singular::DofTopology &topology)
{
  REQUIRE(vertex_ids.size() == static_cast<std::size_t>(mesh.GetNV()));
  REQUIRE(topology.elements.size() == static_cast<std::size_t>(mesh.GetNE()));
  std::map<StableElementKey, ElementKeys> result;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const int *vertices = mesh.GetElement(element)->GetVertices();
    StableElementKey element_key{vertex_ids[vertices[0]], vertex_ids[vertices[1]],
                                 vertex_ids[vertices[2]], vertex_ids[vertices[3]]};
    std::sort(element_key.begin(), element_key.end());
    auto &[h1, nd] = result[element_key];
    for (const auto &dof : topology.elements[element].h1)
    {
      h1.push_back(topology.h1_dofs[dof.dof]);
    }
    for (const auto &dof : topology.elements[element].nd)
    {
      nd.push_back(topology.nd_dofs[dof.dof]);
    }
    std::sort(h1.begin(), h1.end());
    std::sort(nd.begin(), nd.end());
  }
  return result;
}

}  // namespace

TEST_CASE("Triangular singular DOFs preserve the exact enriched gradient",
          "[singulardofs][triangle][Serial]")
{
  auto mesh = InternalLineTipMesh();
  const auto features = fem::singular::ExtractSerialLineTipFeatures(mesh, {7});
  const auto dofs = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);

  REQUIRE(dofs.elements.size() == static_cast<std::size_t>(mesh.GetNE()));
  REQUIRE(dofs.h1_dofs.size() == 6);
  REQUIRE(dofs.nd_dofs.size() == 12);
  REQUIRE(dofs.h1_to_nd.size() == dofs.h1_dofs.size());
  const auto feature_membership =
      fem::singular::BuildTriangleH1DofFeatureMembership(features, dofs);
  REQUIRE(feature_membership.size() == dofs.h1_dofs.size());
  for (std::size_t h1 = 0; h1 < dofs.h1_dofs.size(); h1++)
  {
    REQUIRE(dofs.h1_to_nd[h1] < dofs.nd_dofs.size());
    CHECK(dofs.h1_dofs[h1] == dofs.nd_dofs[dofs.h1_to_nd[h1]]);
    CHECK(dofs.h1_dofs[h1].family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT);
    CHECK(dofs.h1_dofs[h1].singular_entity.size == 1);
    CHECK(dofs.h1_dofs[h1].support_entity.size == 2);
    CHECK(dofs.h1_dofs[h1].component_entity == dofs.h1_dofs[h1].support_entity);
    CHECK(dofs.h1_dofs[h1].interpolation_weights == std::array<int, 4>{1, 1, 0, 0});
    CHECK(feature_membership[h1] == std::vector<std::size_t>{0});
  }

  std::size_t enriched_elements = 0;
  std::vector<int> nd_usage(dofs.nd_dofs.size());
  for (std::size_t element = 0; element < dofs.elements.size(); element++)
  {
    const auto &element_dofs = dofs.elements[element];
    if (features.elements[element].nodes.empty())
    {
      CHECK(element_dofs.h1.empty());
      CHECK(element_dofs.nd.empty());
      continue;
    }
    enriched_elements++;
    REQUIRE(element_dofs.h1.size() == 2);
    REQUIRE(element_dofs.nd.size() == 3);
    for (const auto &h1 : element_dofs.h1)
    {
      const std::size_t nd_global = dofs.h1_to_nd[h1.dof];
      const auto nd =
          std::find_if(element_dofs.nd.begin(), element_dofs.nd.end(),
                       [nd_global](const auto &entry) { return entry.dof == nd_global; });
      REQUIRE(nd != element_dofs.nd.end());
      CHECK(SameBasis(h1.basis, nd->basis));
    }
    for (const auto &nd : element_dofs.nd)
    {
      nd_usage[nd.dof]++;
    }
  }
  CHECK(enriched_elements == 6);
  for (std::size_t nd = 0; nd < dofs.nd_dofs.size(); nd++)
  {
    const auto &key = dofs.nd_dofs[nd];
    if (key.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
    {
      CHECK(nd_usage[nd] == 2);
    }
    else
    {
      CHECK(key.family == fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL);
      CHECK(key.support_entity.size == 3);
      CHECK(key.component_entity.size == 2);
      CHECK(key.interpolation_weights == std::array<int, 4>{1, 1, 1, 0});
      CHECK(nd_usage[nd] == 1);
    }
  }
}

TEST_CASE("Triangular singular DOF keys ignore element traversal order",
          "[singulardofs][triangle][Serial]")
{
  auto forward_mesh = InternalLineTipMesh(false);
  auto reverse_mesh = InternalLineTipMesh(true);
  const auto forward_features =
      fem::singular::ExtractSerialLineTipFeatures(forward_mesh, {7});
  const auto reverse_features =
      fem::singular::ExtractSerialLineTipFeatures(reverse_mesh, {7});
  const auto forward =
      fem::singular::BuildSerialTriangleDofTopology(forward_mesh, forward_features, 1);
  const auto reverse =
      fem::singular::BuildSerialTriangleDofTopology(reverse_mesh, reverse_features, 1);
  CHECK(reverse.h1_dofs == forward.h1_dofs);
  CHECK(reverse.nd_dofs == forward.nd_dofs);
  CHECK(reverse.h1_to_nd == forward.h1_to_nd);
}

TEST_CASE("Triangular singular H1 essential DOFs follow selected PEC traces",
          "[singulardofs][triangle][Serial]")
{
  auto mesh = InternalLineTipMesh();
  auto features = fem::singular::ExtractSerialLineTipFeatures(mesh, {7});
  const auto dofs = fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1);
  const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), dofs);
  const auto essential = fem::singular::GetEssentialTriangleH1TrueDofs(
      Mpi::World(), features, dofs, numbering);
  const auto essential_nd = fem::singular::GetEssentialTriangleNDTrueDofs(
      Mpi::World(), features, dofs, numbering);
  REQUIRE(essential.Size() == 1);
  REQUIRE(essential_nd.Size() == 1);

  const auto &key = dofs.h1_dofs[essential[0]];
  const auto &nd_key = dofs.nd_dofs[essential_nd[0]];
  CHECK(nd_key == key);
  CHECK(nd_key.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT);
  CHECK(key.support_entity.size == 2);
  CHECK(key.support_entity.vertices[0] == 3);
  CHECK(key.support_entity.vertices[1] == 4);

  CHECK_THROWS_AS(fem::singular::BuildSerialTriangleDofTopology(mesh, features, 2),
                  std::invalid_argument);
  REQUIRE_FALSE(features.elements[0].nodes.empty());
  features.elements[0].nodes[0].vertex = features.vertices.size();
  CHECK_THROWS_AS(fem::singular::BuildSerialTriangleDofTopology(mesh, features, 1),
                  std::invalid_argument);
}

TEST_CASE("Triangular singular DOFs preserve canonical ownership across partitions",
          "[singulardofs][triangle][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Triangle ownership is exercised by the two-rank test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);

  auto serial_mesh = InternalLineTipMesh();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});
  const auto serial_dofs =
      fem::singular::BuildSerialTriangleDofTopology(serial_mesh, serial_features, 1);
  const std::array<int, 8> partition{0, 0, 0, 0, 1, 1, 1, 1};
  mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh, partition.data());
  const auto vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
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
      serial_features, parallel_mesh, vertex_ids, source_element_ids);
  const auto local_dofs = fem::singular::BuildLocalTriangleDofTopology(
      parallel_mesh, local_features, vertex_ids, 1);
  const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), local_dofs);

  CHECK(numbering.h1.global_size == static_cast<HYPRE_BigInt>(serial_dofs.h1_dofs.size()));
  CHECK(numbering.nd.global_size == static_cast<HYPRE_BigInt>(serial_dofs.nd_dofs.size()));
  REQUIRE(local_dofs.elements.size() == source_element_ids.size());
  for (std::size_t local_element = 0; local_element < source_element_ids.size();
       local_element++)
  {
    const std::size_t serial_element = source_element_ids[local_element];
    std::vector<fem::singular::DofKey> expected_h1, expected_nd, actual_h1, actual_nd;
    for (const auto &dof : serial_dofs.elements[serial_element].h1)
    {
      expected_h1.push_back(serial_dofs.h1_dofs[dof.dof]);
    }
    for (const auto &dof : serial_dofs.elements[serial_element].nd)
    {
      expected_nd.push_back(serial_dofs.nd_dofs[dof.dof]);
    }
    for (const auto &dof : local_dofs.elements[local_element].h1)
    {
      actual_h1.push_back(local_dofs.h1_dofs[dof.dof]);
    }
    for (const auto &dof : local_dofs.elements[local_element].nd)
    {
      actual_nd.push_back(local_dofs.nd_dofs[dof.dof]);
    }
    std::sort(expected_h1.begin(), expected_h1.end());
    std::sort(expected_nd.begin(), expected_nd.end());
    std::sort(actual_h1.begin(), actual_h1.end());
    std::sort(actual_nd.begin(), actual_nd.end());
    CHECK(actual_h1 == expected_h1);
    CHECK(actual_nd == expected_nd);
  }

  int off_rank_owner = 0;
  for (std::size_t h1 = 0; h1 < local_dofs.h1_dofs.size(); h1++)
  {
    const std::size_t nd = local_dofs.h1_to_nd[h1];
    REQUIRE(nd < local_dofs.nd_dofs.size());
    CHECK(local_dofs.h1_dofs[h1] == local_dofs.nd_dofs[nd]);
    CHECK(numbering.h1_to_nd_true[h1] == numbering.nd.local_to_true[nd]);
    if (numbering.h1.owner[h1] != Mpi::Rank(Mpi::World()))
    {
      off_rank_owner++;
    }
  }
  Mpi::GlobalSum(1, &off_rank_owner, Mpi::World());
  CHECK(off_rank_owner > 0);

  const auto essential = fem::singular::GetEssentialTriangleH1TrueDofs(
      Mpi::World(), local_features, local_dofs, numbering);
  const auto essential_nd = fem::singular::GetEssentialTriangleNDTrueDofs(
      Mpi::World(), local_features, local_dofs, numbering);
  int essential_count = essential.Size();
  int essential_nd_count = essential_nd.Size();
  Mpi::GlobalSum(1, &essential_count, Mpi::World());
  Mpi::GlobalSum(1, &essential_nd_count, Mpi::World());
  CHECK(essential_count == 1);
  CHECK(essential_nd_count == 1);
}

TEST_CASE("First-order singular DOFs preserve the enriched H1 to ND gradient",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8}, 0.5);
  const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, 1);

  REQUIRE(dofs.elements.size() == static_cast<std::size_t>(mesh.GetNE()));
  REQUIRE(dofs.h1_to_nd.size() == dofs.h1_dofs.size());
  for (std::size_t h1 = 0; h1 < dofs.h1_dofs.size(); h1++)
  {
    REQUIRE(dofs.h1_to_nd[h1] < dofs.nd_dofs.size());
    CHECK(dofs.h1_dofs[h1] == dofs.nd_dofs[dofs.h1_to_nd[h1]]);
    CHECK(fem::singular::IsGradientFamily(dofs.h1_dofs[h1].family));
  }

  for (int element = 0; element < mesh.GetNE(); element++)
  {
    const auto &feature_incidence = features.elements[element];
    const auto &element_dofs = dofs.elements[element];
    CHECK(element_dofs.h1.size() ==
          3 * feature_incidence.nodes.size() + 2 * feature_incidence.edges.size());
    CHECK(element_dofs.nd.size() ==
          6 * feature_incidence.nodes.size() + 3 * feature_incidence.edges.size());

    std::set<std::size_t> local_h1, local_nd;
    for (const auto &h1 : element_dofs.h1)
    {
      CHECK(local_h1.insert(h1.dof).second);
      const std::size_t nd_global = dofs.h1_to_nd[h1.dof];
      const auto nd =
          std::find_if(element_dofs.nd.begin(), element_dofs.nd.end(),
                       [nd_global](const auto &entry) { return entry.dof == nd_global; });
      REQUIRE(nd != element_dofs.nd.end());
      CHECK(SameBasis(h1.basis, nd->basis));
    }
    for (const auto &nd : element_dofs.nd)
    {
      CHECK(local_nd.insert(nd.dof).second);
    }
  }
}

TEST_CASE("Singular DOF keys encode retained interpolation entities",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8}, 2.0 / 3.0);

  for (int order = 1; order <= 4; order++)
  {
    CAPTURE(order);
    const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, order);
    CHECK(std::is_sorted(dofs.h1_dofs.begin(), dofs.h1_dofs.end()));
    CHECK(std::adjacent_find(dofs.h1_dofs.begin(), dofs.h1_dofs.end()) ==
          dofs.h1_dofs.end());
    CHECK(std::is_sorted(dofs.nd_dofs.begin(), dofs.nd_dofs.end()));
    CHECK(std::adjacent_find(dofs.nd_dofs.begin(), dofs.nd_dofs.end()) ==
          dofs.nd_dofs.end());

    for (const auto &key : dofs.nd_dofs)
    {
      CAPTURE(static_cast<int>(key.family), key.support_entity.size,
              key.singular_entity.size, key.component_entity.size);
      const bool node_family =
          key.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT ||
          key.family == fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL;
      CHECK(key.singular_entity.size == (node_family ? 1 : 2));
      CHECK(key.support_entity.size >= 2);
      CHECK(key.support_entity.size <= 4);
      int weight_sum = 0;
      for (std::size_t i = 0; i < key.support_entity.size; i++)
      {
        CHECK(key.interpolation_weights[i] > 0);
        weight_sum += key.interpolation_weights[i];
      }
      for (std::size_t i = key.support_entity.size; i < 4; i++)
      {
        CHECK(key.interpolation_weights[i] == 0);
      }
      CHECK(weight_sum == InterpolationDenominator(key));
    }

    std::vector<int> nd_usage(dofs.nd_dofs.size());
    for (const auto &element : dofs.elements)
    {
      for (const auto &dof : element.nd)
      {
        nd_usage[dof.dof]++;
      }
    }
    for (std::size_t i = 0; i < dofs.nd_dofs.size(); i++)
    {
      CHECK(nd_usage[i] >= 1);
      if (dofs.nd_dofs[i].support_entity.size == 4)
      {
        CHECK(nd_usage[i] == 1);
      }
    }
  }
}

TEST_CASE("Singular H1 DOFs have deterministic overlapping feature membership",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});
  REQUIRE(features.features.size() == 4);

  for (int order = 1; order <= 3; order++)
  {
    CAPTURE(order);
    const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, order);
    const auto membership = fem::singular::BuildH1DofFeatureMembership(features, dofs);
    REQUIRE(membership.size() == dofs.h1_dofs.size());

    std::vector<int> feature_coverage(features.features.size());
    for (std::size_t i = 0; i < dofs.h1_dofs.size(); i++)
    {
      const auto &key = dofs.h1_dofs[i];
      const auto &features_for_dof = membership[i];
      CAPTURE(i, static_cast<int>(key.family), key.singular_entity.vertices);
      CHECK(std::is_sorted(features_for_dof.begin(), features_for_dof.end()));
      CHECK(std::adjacent_find(features_for_dof.begin(), features_for_dof.end()) ==
            features_for_dof.end());
      if (key.family == fem::singular::HigherOrderBasisFamily::EDGE_GRADIENT)
      {
        REQUIRE(features_for_dof.size() == 1);
      }
      else
      {
        REQUIRE(key.family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT);
        REQUIRE(features_for_dof.size() == 2);
      }
      for (std::size_t feature : features_for_dof)
      {
        REQUIRE(feature < feature_coverage.size());
        feature_coverage[feature]++;
      }
    }
    CHECK(std::all_of(feature_coverage.begin(), feature_coverage.end(),
                      [](int coverage) { return coverage > 0; }));
  }

  auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, 1);
  REQUIRE_FALSE(dofs.h1_dofs.empty());
  dofs.h1_dofs.front().singular_entity.size = 3;
  CHECK_THROWS_AS(fem::singular::BuildH1DofFeatureMembership(features, dofs),
                  std::invalid_argument);
}

TEST_CASE("Singular DOF numbering ignores element traversal order",
          "[singulardofs][Serial]")
{
  auto forward_mesh = RectangleSheetMesh(false);
  auto reverse_mesh = RectangleSheetMesh(true);
  const auto forward_features =
      fem::singular::ExtractSerialSheetFeatures(forward_mesh, {7, 8});
  const auto reverse_features =
      fem::singular::ExtractSerialSheetFeatures(reverse_mesh, {7, 8});

  for (int order = 1; order <= 3; order++)
  {
    CAPTURE(order);
    const auto forward =
        fem::singular::BuildSerialDofTopology(forward_mesh, forward_features, order);
    const auto reverse =
        fem::singular::BuildSerialDofTopology(reverse_mesh, reverse_features, order);
    CHECK(reverse.h1_dofs == forward.h1_dofs);
    CHECK(reverse.nd_dofs == forward.nd_dofs);
    CHECK(reverse.h1_to_nd == forward.h1_to_nd);
    CHECK(GlobalElementDofs(reverse_mesh, reverse) ==
          GlobalElementDofs(forward_mesh, forward));
  }
}

TEST_CASE("Singular DOF keys produce matching H1 and tangential ND traces",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});

  for (int order = 1; order <= 3; order++)
  {
    CAPTURE(order);
    const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, order);
    for (int first = 0; first < mesh.GetNE(); first++)
    {
      std::set<int> first_vertices(mesh.GetElement(first)->GetVertices(),
                                   mesh.GetElement(first)->GetVertices() + 4);
      for (int second = first + 1; second < mesh.GetNE(); second++)
      {
        std::vector<int> intersection;
        const int *second_vertices = mesh.GetElement(second)->GetVertices();
        for (int i = 0; i < 4; i++)
        {
          if (first_vertices.count(second_vertices[i]) != 0)
          {
            intersection.push_back(second_vertices[i]);
          }
        }
        if (intersection.size() != 3)
        {
          continue;
        }
        std::sort(intersection.begin(), intersection.end());
        const std::array<int, 3> face{intersection[0], intersection[1], intersection[2]};
        const auto lambda_first = FacePoint(*mesh.GetElement(first), face);
        const auto lambda_second = FacePoint(*mesh.GetElement(second), face);
        const auto gradients_first = PhysicalBarycentricGradients(mesh, first);
        const auto gradients_second = PhysicalBarycentricGradients(mesh, second);
        const auto tangent_1 = MeshEdge(mesh, face[0], face[1]);
        const auto tangent_2 = MeshEdge(mesh, face[0], face[2]);

        using Trace = std::array<double, 2>;
        std::map<fem::singular::DofKey, Trace> first_nd, second_nd;
        std::map<fem::singular::DofKey, double> first_h1, second_h1;
        auto add_traces = [&](int element, const fem::singular::BarycentricPoint &lambda,
                              const fem::singular::BarycentricGradients &gradients,
                              auto &h1_traces, auto &nd_traces)
        {
          for (const auto &dof : dofs.elements[element].h1)
          {
            const auto &key = dofs.h1_dofs[dof.dof];
            if (IsSupportedOnFace(key.support_entity, face))
            {
              h1_traces.emplace(key, EvaluateGradientPotential(lambda, dof.basis));
            }
          }
          for (const auto &dof : dofs.elements[element].nd)
          {
            const auto &key = dofs.nd_dofs[dof.dof];
            if (IsSupportedOnFace(key.support_entity, face))
            {
              const auto value =
                  fem::singular::EvaluateHigherOrderBasis(lambda, gradients, dof.basis);
              nd_traces.emplace(
                  key, Trace{Dot(value.value, tangent_1), Dot(value.value, tangent_2)});
            }
          }
        };
        add_traces(first, lambda_first, gradients_first, first_h1, first_nd);
        add_traces(second, lambda_second, gradients_second, second_h1, second_nd);

        REQUIRE(first_h1.size() == second_h1.size());
        REQUIRE(first_nd.size() == second_nd.size());
        for (const auto &[key, value] : first_h1)
        {
          const auto other = second_h1.find(key);
          REQUIRE(other != second_h1.end());
          CHECK(std::abs(value - other->second) <= 2.0e-11 * (1.0 + std::abs(value)));
        }
        for (const auto &[key, value] : first_nd)
        {
          const auto other = second_nd.find(key);
          REQUIRE(other != second_nd.end());
          for (int component = 0; component < 2; component++)
          {
            CAPTURE(first, second, component, static_cast<int>(key.family));
            CHECK(std::abs(value[component] - other->second[component]) <=
                  2.0e-10 * (1.0 + std::abs(value[component])));
          }
        }
      }
    }
  }
}

TEST_CASE("Random singular ND fields are globally Hcurl conforming",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  const auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});
  constexpr std::array<std::array<double, 3>, 5> face_weights{
      std::array<double, 3>{0.217, 0.331, 0.452},
      std::array<double, 3>{0.113, 0.521, 0.366},
      std::array<double, 3>{0.587, 0.172, 0.241},
      std::array<double, 3>{0.302, 0.619, 0.079},
      std::array<double, 3>{0.431, 0.293, 0.276}};

  for (int order = 1; order <= 3; order++)
  {
    CAPTURE(order);
    const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, order);
    std::vector<double> coefficients(dofs.nd_dofs.size());
    for (std::size_t i = 0; i < coefficients.size(); i++)
    {
      coefficients[i] = std::sin(0.731 * static_cast<double>(i + 1)) +
                        0.37 * std::cos(1.193 * static_cast<double>(i + 1));
    }

    int shared_faces = 0;
    for (int first = 0; first < mesh.GetNE(); first++)
    {
      std::set<int> first_vertices(mesh.GetElement(first)->GetVertices(),
                                   mesh.GetElement(first)->GetVertices() + 4);
      for (int second = first + 1; second < mesh.GetNE(); second++)
      {
        std::vector<int> intersection;
        const int *second_vertices = mesh.GetElement(second)->GetVertices();
        for (int i = 0; i < 4; i++)
        {
          if (first_vertices.count(second_vertices[i]) != 0)
          {
            intersection.push_back(second_vertices[i]);
          }
        }
        if (intersection.size() != 3)
        {
          continue;
        }
        shared_faces++;
        std::sort(intersection.begin(), intersection.end());
        const std::array<int, 3> face{intersection[0], intersection[1], intersection[2]};
        const auto gradients_first = PhysicalBarycentricGradients(mesh, first);
        const auto gradients_second = PhysicalBarycentricGradients(mesh, second);
        const std::array<fem::singular::Vector3, 2> tangents{
            MeshEdge(mesh, face[0], face[1]), MeshEdge(mesh, face[0], face[2])};

        const auto face_point =
            [&face](const mfem::Element &element, const std::array<double, 3> &weights)
        {
          fem::singular::BarycentricPoint lambda{};
          for (int i = 0; i < 3; i++)
          {
            const int *vertex =
                std::find(element.GetVertices(), element.GetVertices() + 4, face[i]);
            REQUIRE(vertex != element.GetVertices() + 4);
            lambda[vertex - element.GetVertices()] = weights[i];
          }
          return lambda;
        };
        const auto evaluate = [&](int element,
                                  const fem::singular::BarycentricPoint &lambda,
                                  const fem::singular::BarycentricGradients &gradients)
        {
          fem::singular::Vector3 value{};
          for (const auto &dof : dofs.elements[element].nd)
          {
            const auto basis =
                fem::singular::EvaluateHigherOrderBasis(lambda, gradients, dof.basis);
            for (int d = 0; d < 3; d++)
            {
              value[d] += coefficients[dof.dof] * basis.value[d];
            }
          }
          return value;
        };

        for (const auto &weights : face_weights)
        {
          const auto lambda_first = face_point(*mesh.GetElement(first), weights);
          const auto lambda_second = face_point(*mesh.GetElement(second), weights);
          const auto value_first = evaluate(first, lambda_first, gradients_first);
          const auto value_second = evaluate(second, lambda_second, gradients_second);
          for (int component = 0; component < 2; component++)
          {
            const double trace_first = Dot(value_first, tangents[component]);
            const double trace_second = Dot(value_second, tangents[component]);
            const double scale =
                1.0 + std::max(std::abs(trace_first), std::abs(trace_second));
            CAPTURE(first, second, component, weights);
            CHECK(std::abs(trace_first - trace_second) <= 5.0e-10 * scale);
          }
        }
      }
    }
    CHECK(shared_faces > 0);
  }
}

TEST_CASE("Singular H1 essential DOFs are classified by PEC sheet trace",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});

  for (int order = 1; order <= 3; order++)
  {
    CAPTURE(order);
    const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, order);
    const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), dofs);
    const auto essential =
        fem::singular::GetEssentialH1TrueDofs(Mpi::World(), features, dofs, numbering);
    const auto essential_nd =
        fem::singular::GetEssentialNDTrueDofs(Mpi::World(), features, dofs, numbering);
    std::vector<bool> classified(dofs.h1_dofs.size());
    std::vector<bool> classified_nd(dofs.nd_dofs.size());
    for (int true_dof : essential)
    {
      REQUIRE(true_dof >= 0);
      REQUIRE(true_dof < static_cast<int>(classified.size()));
      classified[true_dof] = true;
    }
    for (int true_dof : essential_nd)
    {
      REQUIRE(true_dof >= 0);
      REQUIRE(true_dof < static_cast<int>(classified_nd.size()));
      classified_nd[true_dof] = true;
    }

    std::vector<bool> nonzero_sheet_trace(dofs.h1_dofs.size());
    std::vector<bool> nonzero_sheet_tangential_trace(dofs.nd_dofs.size());
    for (int element = 0; element < mesh.GetNE(); element++)
    {
      const auto &tetrahedron = *mesh.GetElement(element);
      const auto grad_lambda = PhysicalBarycentricGradients(mesh, element);
      std::set<int> element_vertices(tetrahedron.GetVertices(),
                                     tetrahedron.GetVertices() + 4);
      for (const auto &face : features.sheet_faces)
      {
        if (!std::all_of(face.mesh_vertices.begin(), face.mesh_vertices.end(),
                         [&element_vertices](auto vertex)
                         { return element_vertices.count(static_cast<int>(vertex)) != 0; }))
        {
          continue;
        }
        const std::array<int, 3> face_vertices{static_cast<int>(face.mesh_vertices[0]),
                                               static_cast<int>(face.mesh_vertices[1]),
                                               static_cast<int>(face.mesh_vertices[2])};
        const auto lambda = FacePoint(tetrahedron, face_vertices);
        for (const auto &dof : dofs.elements[element].h1)
        {
          const double trace = EvaluateGradientPotential(lambda, dof.basis);
          if (std::abs(trace) > 1.0e-12)
          {
            nonzero_sheet_trace[dof.dof] = true;
          }
        }
        const auto tangent_0 = MeshEdge(mesh, face_vertices[0], face_vertices[1]);
        const auto tangent_1 = MeshEdge(mesh, face_vertices[0], face_vertices[2]);
        for (const auto &dof : dofs.elements[element].nd)
        {
          const auto trace =
              fem::singular::EvaluateHigherOrderBasis(lambda, grad_lambda, dof.basis);
          if (std::abs(Dot(trace.value, tangent_0)) > 1.0e-12 ||
              std::abs(Dot(trace.value, tangent_1)) > 1.0e-12)
          {
            nonzero_sheet_tangential_trace[dof.dof] = true;
          }
        }
      }
    }

    REQUIRE(numbering.h1.local_to_true.size() == dofs.h1_dofs.size());
    for (std::size_t local = 0; local < dofs.h1_dofs.size(); local++)
    {
      REQUIRE(numbering.h1.local_to_true[local] >= 0);
      REQUIRE(numbering.h1.local_to_true[local] <
              static_cast<HYPRE_BigInt>(classified.size()));
      const std::size_t true_dof = numbering.h1.local_to_true[local];
      CAPTURE(local, true_dof, dofs.h1_dofs[local].support_entity.size);
      CHECK(classified[true_dof] == nonzero_sheet_trace[local]);
    }
    REQUIRE(numbering.nd.local_to_true.size() == dofs.nd_dofs.size());
    for (std::size_t local = 0; local < dofs.nd_dofs.size(); local++)
    {
      REQUIRE(numbering.nd.local_to_true[local] >= 0);
      REQUIRE(numbering.nd.local_to_true[local] <
              static_cast<HYPRE_BigInt>(classified_nd.size()));
      const std::size_t true_dof = numbering.nd.local_to_true[local];
      CAPTURE(local, true_dof, dofs.nd_dofs[local].support_entity.size,
              static_cast<int>(dofs.nd_dofs[local].family));
      CHECK(classified_nd[true_dof] == nonzero_sheet_tangential_trace[local]);
    }
  }

  features.sheet_faces.push_back(features.sheet_faces.front());
  const auto dofs = fem::singular::BuildSerialDofTopology(mesh, features, 1);
  const auto numbering = fem::singular::BuildParallelDofNumbering(Mpi::World(), dofs);
  CHECK_THROWS_AS(
      fem::singular::GetEssentialH1TrueDofs(Mpi::World(), features, dofs, numbering),
      std::invalid_argument);
  CHECK_THROWS_AS(
      fem::singular::GetEssentialNDTrueDofs(Mpi::World(), features, dofs, numbering),
      std::invalid_argument);
}

TEST_CASE("Singular DOF enumeration rejects inconsistent topology",
          "[singulardofs][Serial]")
{
  auto mesh = RectangleSheetMesh();
  auto features = fem::singular::ExtractSerialSheetFeatures(mesh, {7, 8});
  CHECK_THROWS_AS(fem::singular::BuildSerialDofTopology(mesh, features, 0),
                  std::invalid_argument);

  REQUIRE_FALSE(features.elements[0].nodes.empty());
  features.elements[0].nodes[0].vertex = features.vertices.size();
  CHECK_THROWS_AS(fem::singular::BuildSerialDofTopology(mesh, features, 1),
                  std::invalid_argument);
}

TEST_CASE("Original serial vertex IDs survive different mesh partitions",
          "[singulardofs][Parallel]")
{
  if (Mpi::Size(Mpi::World()) == 1)
  {
    SUCCEED("Partition-ID propagation is exercised by the [Parallel] test run.");
    return;
  }
  REQUIRE(Mpi::Size(Mpi::World()) == 2);
  auto serial_mesh = RectangleSheetMesh();
  const auto serial_features =
      fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7, 8});
  std::vector<fem::singular::GlobalVertexId> serial_vertex_ids(serial_mesh.GetNV());
  for (int vertex = 0; vertex < serial_mesh.GetNV(); vertex++)
  {
    serial_vertex_ids[vertex] = vertex;
  }
  const std::array<std::array<int, 4>, 2> partitions{std::array<int, 4>{0, 0, 1, 1},
                                                     std::array<int, 4>{0, 1, 0, 1}};

  for (const auto &partition : partitions)
  {
    CAPTURE(partition);
    mfem::ParMesh parallel_mesh(Mpi::World(), serial_mesh, partition.data());
    const auto vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
        serial_mesh, parallel_mesh, partition.data());
    REQUIRE(vertex_ids.size() == static_cast<std::size_t>(parallel_mesh.GetNV()));
    const auto local_features = fem::singular::DistributeSerialSheetFeatures(
        serial_mesh, serial_features, parallel_mesh, vertex_ids);
    std::vector<fem::singular::GlobalVertexId> source_element_ids;
    for (int element = 0; element < serial_mesh.GetNE(); element++)
    {
      if (partition[element] == Mpi::Rank(Mpi::World()))
      {
        source_element_ids.push_back(element);
      }
    }
    const auto mapped_features = fem::singular::DistributeSerialSheetFeatures(
        serial_features, parallel_mesh, vertex_ids, source_element_ids);
    CHECK(mapped_features.elements.size() == local_features.elements.size());
    CHECK(local_features.vertices.size() == serial_features.vertices.size());
    CHECK(local_features.segments.size() == serial_features.segments.size());
    CHECK(local_features.features.size() == serial_features.features.size());
    CHECK(local_features.sheet_faces.size() == serial_features.sheet_faces.size());
    CHECK(local_features.elements.size() ==
          static_cast<std::size_t>(parallel_mesh.GetNE()));

    int local_element = 0;
    for (int serial_element = 0; serial_element < serial_mesh.GetNE(); serial_element++)
    {
      if (partition[serial_element] != Mpi::Rank(Mpi::World()))
      {
        continue;
      }
      REQUIRE(local_element < parallel_mesh.GetNE());
      const int *serial_vertices = serial_mesh.GetElement(serial_element)->GetVertices();
      const int *parallel_vertices =
          parallel_mesh.GetElement(local_element++)->GetVertices();
      for (int i = 0; i < 4; i++)
      {
        CHECK(vertex_ids[parallel_vertices[i]] == serial_vertices[i]);
      }
    }
    CHECK(local_element == parallel_mesh.GetNE());

    auto duplicate_vertex_ids = vertex_ids;
    REQUIRE(duplicate_vertex_ids.size() >= 2);
    duplicate_vertex_ids[1] = duplicate_vertex_ids[0];
    CHECK_THROWS_AS(fem::singular::DistributeSerialSheetFeatures(
                        serial_mesh, serial_features, parallel_mesh, duplicate_vertex_ids),
                    std::invalid_argument);

    for (int order = 1; order <= 3; order++)
    {
      CAPTURE(order);
      const auto serial_dofs =
          fem::singular::BuildSerialDofTopology(serial_mesh, serial_features, order);
      const auto local_dofs = fem::singular::BuildLocalDofTopology(
          parallel_mesh, local_features, vertex_ids, order);
      const auto mapped_dofs = fem::singular::BuildLocalDofTopology(
          parallel_mesh, mapped_features, vertex_ids, order);
      CHECK(mapped_dofs.h1_dofs == local_dofs.h1_dofs);
      CHECK(mapped_dofs.nd_dofs == local_dofs.nd_dofs);
      CHECK(mapped_dofs.h1_to_nd == local_dofs.h1_to_nd);
      const auto numbering =
          fem::singular::BuildParallelDofNumbering(Mpi::World(), local_dofs);
      const auto essential = fem::singular::GetEssentialH1TrueDofs(
          Mpi::World(), local_features, local_dofs, numbering);
      const auto serial_element_dofs =
          StableElementDofs(serial_mesh, serial_vertex_ids, serial_dofs);
      const auto local_element_dofs =
          StableElementDofs(parallel_mesh, vertex_ids, local_dofs);
      for (const auto &[element, dofs] : local_element_dofs)
      {
        const auto serial = serial_element_dofs.find(element);
        REQUIRE(serial != serial_element_dofs.end());
        CHECK(dofs == serial->second);
      }

      int off_rank_owners = 0;
      auto check_parallel_map = [&](const std::vector<fem::singular::DofKey> &serial_keys,
                                    const std::vector<fem::singular::DofKey> &local_keys,
                                    const fem::singular::TrueDofMap &true_dofs, bool h1)
      {
        std::vector<std::set<int>> expected_ranks(serial_keys.size());
        for (int element = 0; element < serial_mesh.GetNE(); element++)
        {
          const auto &element_dofs =
              h1 ? serial_dofs.elements[element].h1 : serial_dofs.elements[element].nd;
          for (const auto &dof : element_dofs)
          {
            expected_ranks[dof.dof].insert(partition[element]);
          }
        }

        std::vector<int> expected_owner(serial_keys.size());
        std::vector<HYPRE_BigInt> expected_true(serial_keys.size());
        std::vector<HYPRE_BigInt> owned_size(Mpi::Size(Mpi::World()));
        for (std::size_t i = 0; i < serial_keys.size(); i++)
        {
          REQUIRE_FALSE(expected_ranks[i].empty());
          expected_owner[i] = *expected_ranks[i].begin();
          owned_size[expected_owner[i]]++;
        }
        std::vector<HYPRE_BigInt> owned_offset(owned_size.size());
        for (std::size_t rank = 1; rank < owned_size.size(); rank++)
        {
          owned_offset[rank] = owned_offset[rank - 1] + owned_size[rank - 1];
        }
        auto owner_position = std::vector<HYPRE_BigInt>(owned_size.size());
        for (std::size_t i = 0; i < serial_keys.size(); i++)
        {
          expected_true[i] =
              owned_offset[expected_owner[i]] + owner_position[expected_owner[i]]++;
        }

        CHECK(true_dofs.global_size == static_cast<HYPRE_BigInt>(serial_keys.size()));
        std::vector<HYPRE_BigInt> local_sizes(Mpi::Size(Mpi::World()));
        const HYPRE_BigInt local_size = static_cast<HYPRE_BigInt>(local_keys.size());
        Mpi::Allgather(1, &local_size, local_sizes.data(), Mpi::World());
        HYPRE_BigInt expected_local_offset = 0;
        HYPRE_BigInt expected_global_local_size = 0;
        for (int rank = 0; rank < Mpi::Size(Mpi::World()); rank++)
        {
          if (rank < Mpi::Rank(Mpi::World()))
          {
            expected_local_offset += local_sizes[rank];
          }
          expected_global_local_size += local_sizes[rank];
        }
        CHECK(true_dofs.local_size == local_size);
        CHECK(true_dofs.local_offset == expected_local_offset);
        CHECK(true_dofs.global_local_size == expected_global_local_size);
        CHECK(true_dofs.owned_offset == owned_offset[Mpi::Rank(Mpi::World())]);
        CHECK(true_dofs.owned_size == owned_size[Mpi::Rank(Mpi::World())]);
        REQUIRE(true_dofs.owner.size() == local_keys.size());
        REQUIRE(true_dofs.local_to_true.size() == local_keys.size());

        std::vector<fem::singular::DofKey> expected_local;
        std::vector<int> rank_count(serial_keys.size());
        for (std::size_t i = 0; i < serial_keys.size(); i++)
        {
          if (expected_ranks[i].count(Mpi::Rank(Mpi::World())) != 0)
          {
            expected_local.push_back(serial_keys[i]);
          }
          if (std::binary_search(local_keys.begin(), local_keys.end(), serial_keys[i]))
          {
            rank_count[i] = 1;
          }
        }
        CHECK(local_keys == expected_local);

        std::vector<HYPRE_BigInt> minimum_true(serial_keys.size(),
                                               std::numeric_limits<HYPRE_BigInt>::max());
        std::vector<HYPRE_BigInt> maximum_true(serial_keys.size(), -1);
        std::vector<int> owned_true(serial_keys.size());
        for (std::size_t local = 0; local < local_keys.size(); local++)
        {
          const auto serial =
              std::lower_bound(serial_keys.begin(), serial_keys.end(), local_keys[local]);
          REQUIRE(serial != serial_keys.end());
          REQUIRE(*serial == local_keys[local]);
          const std::size_t i = serial - serial_keys.begin();
          CHECK(true_dofs.owner[local] == expected_owner[i]);
          CHECK(true_dofs.local_to_true[local] == expected_true[i]);
          minimum_true[i] = true_dofs.local_to_true[local];
          maximum_true[i] = true_dofs.local_to_true[local];
          if (true_dofs.owner[local] == Mpi::Rank(Mpi::World()))
          {
            REQUIRE(true_dofs.local_to_true[local] >= 0);
            REQUIRE(true_dofs.local_to_true[local] <
                    static_cast<HYPRE_BigInt>(owned_true.size()));
            owned_true[true_dofs.local_to_true[local]] = 1;
          }
          else
          {
            off_rank_owners++;
          }
        }

        Mpi::GlobalSum(static_cast<int>(rank_count.size()), rank_count.data(),
                       Mpi::World());
        Mpi::GlobalMin(static_cast<int>(minimum_true.size()), minimum_true.data(),
                       Mpi::World());
        Mpi::GlobalMax(static_cast<int>(maximum_true.size()), maximum_true.data(),
                       Mpi::World());
        Mpi::GlobalSum(static_cast<int>(owned_true.size()), owned_true.data(),
                       Mpi::World());
        for (std::size_t i = 0; i < serial_keys.size(); i++)
        {
          CAPTURE(h1, i);
          CHECK(rank_count[i] == static_cast<int>(expected_ranks[i].size()));
          CHECK(minimum_true[i] == expected_true[i]);
          CHECK(maximum_true[i] == expected_true[i]);
          CHECK(owned_true[i] == 1);
        }
      };
      check_parallel_map(serial_dofs.h1_dofs, local_dofs.h1_dofs, numbering.h1, true);
      check_parallel_map(serial_dofs.nd_dofs, local_dofs.nd_dofs, numbering.nd, false);
      Mpi::GlobalSum(1, &off_rank_owners, Mpi::World());
      CHECK(off_rank_owners > 0);

      int local_essential_count = essential.Size();
      Mpi::GlobalSum(1, &local_essential_count, Mpi::World());
      int expected_essential_count = 0;
      for (const auto &key : serial_dofs.h1_dofs)
      {
        const bool expected = std::any_of(
            serial_features.sheet_faces.begin(), serial_features.sheet_faces.end(),
            [&key](const auto &face)
            {
              return key.support_entity.size <= face.mesh_vertices.size() &&
                     std::all_of(
                         key.support_entity.vertices.begin(),
                         key.support_entity.vertices.begin() + key.support_entity.size,
                         [&face](auto vertex)
                         {
                           return std::binary_search(face.mesh_vertices.begin(),
                                                     face.mesh_vertices.end(), vertex);
                         });
            });
        expected_essential_count += expected;
      }
      CHECK(local_essential_count == expected_essential_count);
      for (std::size_t local = 0; local < local_dofs.h1_dofs.size(); local++)
      {
        const auto &key = local_dofs.h1_dofs[local];
        const bool expected = std::any_of(
            local_features.sheet_faces.begin(), local_features.sheet_faces.end(),
            [&key](const auto &face)
            {
              return key.support_entity.size <= face.mesh_vertices.size() &&
                     std::all_of(
                         key.support_entity.vertices.begin(),
                         key.support_entity.vertices.begin() + key.support_entity.size,
                         [&face](auto vertex)
                         {
                           return std::binary_search(face.mesh_vertices.begin(),
                                                     face.mesh_vertices.end(), vertex);
                         });
            });
        const HYPRE_BigInt owned_true =
            numbering.h1.local_to_true[local] - numbering.h1.owned_offset;
        const bool present = owned_true >= 0 && owned_true < numbering.h1.owned_size &&
                             std::find(essential.begin(), essential.end(),
                                       static_cast<int>(owned_true)) != essential.end();
        CHECK(present ==
              (expected && numbering.h1.owner[local] == Mpi::Rank(Mpi::World())));
      }

      REQUIRE(local_dofs.h1_to_nd.size() == local_dofs.h1_dofs.size());
      REQUIRE(numbering.h1_to_nd_true.size() == local_dofs.h1_dofs.size());
      for (std::size_t h1 = 0; h1 < local_dofs.h1_dofs.size(); h1++)
      {
        REQUIRE(local_dofs.h1_to_nd[h1] < local_dofs.nd_dofs.size());
        CHECK(local_dofs.h1_dofs[h1] == local_dofs.nd_dofs[local_dofs.h1_to_nd[h1]]);
        CHECK(numbering.h1_to_nd_true[h1] ==
              numbering.nd.local_to_true[local_dofs.h1_to_nd[h1]]);
        CHECK(numbering.h1.owner[h1] == numbering.nd.owner[local_dofs.h1_to_nd[h1]]);
      }

      if (order == 1)
      {
        auto invalid_dofs = local_dofs;
        REQUIRE_FALSE(invalid_dofs.h1_to_nd.empty());
        invalid_dofs.h1_to_nd.pop_back();
        CHECK_THROWS_AS(
            fem::singular::BuildParallelDofNumbering(Mpi::World(), invalid_dofs),
            std::invalid_argument);
      }
    }
  }
}

}  // namespace palace
