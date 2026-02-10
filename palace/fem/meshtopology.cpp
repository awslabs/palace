// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "meshtopology.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>
#include <mfem.hpp>

namespace palace
{

int MeshTopology::GetMidpoint(int v0, int v1)
{
  auto key = EdgeKey(v0, v1);
  auto it = edge_midpoints.find(key);
  if (it != edge_midpoints.end())
  {
    return it->second;
  }
  // Create new midpoint vertex.
  int mid = static_cast<int>(vertices.size());
  std::array<double, 3> coords;
  for (int d = 0; d < 3; d++)
  {
    coords[d] = 0.5 * (vertices[v0][d] + vertices[v1][d]);
  }
  vertices.push_back(coords);
  edge_midpoints[key] = mid;
  return mid;
}

int MeshTopology::GetFaceCenter(int v0, int v1, int v2, int v3)
{
  std::array<int, 4> key = {v0, v1, v2, v3};
  std::sort(key.begin(), key.end());
  auto it = face_centers.find(key);
  if (it != face_centers.end())
  {
    return it->second;
  }
  int fc = static_cast<int>(vertices.size());
  std::array<double, 3> coords;
  for (int d = 0; d < 3; d++)
  {
    coords[d] =
        0.25 * (vertices[v0][d] + vertices[v1][d] + vertices[v2][d] + vertices[v3][d]);
  }
  vertices.push_back(coords);
  face_centers[key] = fc;
  return fc;
}

int MeshTopology::GetBodyCenter(const std::vector<int> &hex_vertices)
{
  std::array<int, 8> key;
  std::copy(hex_vertices.begin(), hex_vertices.end(), key.begin());
  std::sort(key.begin(), key.end());
  auto it = body_centers.find(key);
  if (it != body_centers.end())
  {
    return it->second;
  }
  int bc = static_cast<int>(vertices.size());
  std::array<double, 3> coords = {0.0, 0.0, 0.0};
  for (int i = 0; i < 8; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      coords[d] += vertices[hex_vertices[i]][d];
    }
  }
  for (int d = 0; d < 3; d++)
  {
    coords[d] *= 0.125;
  }
  vertices.push_back(coords);
  body_centers[key] = bc;
  return bc;
}

std::array<int, 8> MeshTopology::RefineHexElement(int elem_idx)
{
  MFEM_ASSERT(elements[elem_idx].active,
              "Cannot refine an inactive (already refined) element!");
  MFEM_ASSERT(elements[elem_idx].vertices.size() == 8,
              "RefineHexElement requires a hex element with 8 vertices!");

  // Capture parent data before any push_back (which may reallocate).
  const auto parent_verts = elements[elem_idx].vertices;
  int attr = elements[elem_idx].attribute;
  int level = elements[elem_idx].level;

  // MFEM hex vertex ordering:
  //      7--------6
  //     /|       /|
  //    4--------5 |
  //    | |      | |
  //    | 3------|-2
  //    |/       |/
  //    0--------1
  int v0 = parent_verts[0], v1 = parent_verts[1];
  int v2 = parent_verts[2], v3 = parent_verts[3];
  int v4 = parent_verts[4], v5 = parent_verts[5];
  int v6 = parent_verts[6], v7 = parent_verts[7];

  // 12 edge midpoints.
  int e01 = GetMidpoint(v0, v1);
  int e12 = GetMidpoint(v1, v2);
  int e23 = GetMidpoint(v2, v3);
  int e30 = GetMidpoint(v3, v0);
  int e45 = GetMidpoint(v4, v5);
  int e56 = GetMidpoint(v5, v6);
  int e67 = GetMidpoint(v6, v7);
  int e74 = GetMidpoint(v7, v4);
  int e04 = GetMidpoint(v0, v4);
  int e15 = GetMidpoint(v1, v5);
  int e26 = GetMidpoint(v2, v6);
  int e37 = GetMidpoint(v3, v7);

  // 6 face centers.
  int fc_bot = GetFaceCenter(v0, v1, v2, v3);  // bottom z=0
  int fc_top = GetFaceCenter(v4, v5, v6, v7);  // top z=1
  int fc_fro = GetFaceCenter(v0, v1, v5, v4);  // front y=0
  int fc_bac = GetFaceCenter(v2, v3, v7, v6);  // back y=1
  int fc_lef = GetFaceCenter(v0, v3, v7, v4);  // left x=0
  int fc_rig = GetFaceCenter(v1, v2, v6, v5);  // right x=1

  // 1 body center.
  int bc = GetBodyCenter(parent_verts);

  // Reserve space for 8 children to avoid reallocation during push_back.
  elements.reserve(elements.size() + 8);

  // Deactivate parent BEFORE push_back (push_back may reallocate and invalidate refs).
  elements[elem_idx].active = false;

  int child0_idx = static_cast<int>(elements.size());

  // 8 children, each a hex with MFEM vertex ordering.
  // Child 0: corner v0
  elements.push_back(
      {{v0, e01, fc_bot, e30, e04, fc_fro, bc, fc_lef}, attr, elem_idx, level + 1, true});
  // Child 1: corner v1
  elements.push_back(
      {{e01, v1, e12, fc_bot, fc_fro, e15, fc_rig, bc}, attr, elem_idx, level + 1, true});
  // Child 2: corner v2
  elements.push_back(
      {{fc_bot, e12, v2, e23, bc, fc_rig, e26, fc_bac}, attr, elem_idx, level + 1, true});
  // Child 3: corner v3
  elements.push_back(
      {{e30, fc_bot, e23, v3, fc_lef, bc, fc_bac, e37}, attr, elem_idx, level + 1, true});
  // Child 4: corner v4
  elements.push_back(
      {{e04, fc_fro, bc, fc_lef, v4, e45, fc_top, e74}, attr, elem_idx, level + 1, true});
  // Child 5: corner v5
  elements.push_back(
      {{fc_fro, e15, fc_rig, bc, e45, v5, e56, fc_top}, attr, elem_idx, level + 1, true});
  // Child 6: corner v6
  elements.push_back(
      {{bc, fc_rig, e26, fc_bac, fc_top, e56, v6, e67}, attr, elem_idx, level + 1, true});
  // Child 7: corner v7
  elements.push_back(
      {{fc_lef, bc, fc_bac, e37, e74, fc_top, e67, v7}, attr, elem_idx, level + 1, true});

  std::array<int, 8> children;
  for (int i = 0; i < 8; i++)
  {
    children[i] = child0_idx + i;
  }

  // Record history with sentinel values for bisection-specific fields.
  history.push_back({elem_idx,
                     {children[0], children[1], children[2], children[3], children[4],
                      children[5], children[6], children[7]},
                     -1,
                     -1,
                     -1});

  return children;
}

std::array<int, 4> MeshTopology::RefineQuadElement(int elem_idx)
{
  MFEM_ASSERT(elements[elem_idx].active,
              "Cannot refine an inactive (already refined) element!");
  MFEM_ASSERT(elements[elem_idx].vertices.size() == 4,
              "RefineQuadElement requires a quad element with 4 vertices!");

  // Capture parent data before any push_back.
  const auto parent_verts = elements[elem_idx].vertices;
  int attr = elements[elem_idx].attribute;
  int level = elements[elem_idx].level;

  // MFEM quad vertex ordering:
  //  3--------2
  //  |        |
  //  |        |
  //  0--------1
  int v0 = parent_verts[0], v1 = parent_verts[1];
  int v2 = parent_verts[2], v3 = parent_verts[3];

  // 4 edge midpoints.
  int e01 = GetMidpoint(v0, v1);
  int e12 = GetMidpoint(v1, v2);
  int e23 = GetMidpoint(v2, v3);
  int e30 = GetMidpoint(v3, v0);

  // Face center.
  int fc = GetFaceCenter(v0, v1, v2, v3);

  // Reserve space for 4 children.
  elements.reserve(elements.size() + 4);

  // Deactivate parent.
  elements[elem_idx].active = false;

  int child0_idx = static_cast<int>(elements.size());

  // 4 children, each a quad with MFEM vertex ordering.
  elements.push_back({{v0, e01, fc, e30}, attr, elem_idx, level + 1, true});
  elements.push_back({{e01, v1, e12, fc}, attr, elem_idx, level + 1, true});
  elements.push_back({{fc, e12, v2, e23}, attr, elem_idx, level + 1, true});
  elements.push_back({{e30, fc, e23, v3}, attr, elem_idx, level + 1, true});

  std::array<int, 4> children;
  for (int i = 0; i < 4; i++)
  {
    children[i] = child0_idx + i;
  }

  // Record history with sentinel values for bisection-specific fields.
  history.push_back(
      {elem_idx, {children[0], children[1], children[2], children[3]}, -1, -1, -1});

  return children;
}

std::pair<int, int> MeshTopology::FindLongestEdge(int elem_idx) const
{
  const auto &elem = elements[elem_idx];
  int nv = static_cast<int>(elem.vertices.size());
  double max_len2 = -1.0;
  int best_a = 0, best_b = 1;

  for (int i = 0; i < nv; i++)
  {
    for (int j = i + 1; j < nv; j++)
    {
      int vi = elem.vertices[i], vj = elem.vertices[j];
      double len2 = 0.0;
      for (int d = 0; d < sdim; d++)
      {
        double diff = vertices[vi][d] - vertices[vj][d];
        len2 += diff * diff;
      }
      if (len2 > max_len2)
      {
        max_len2 = len2;
        best_a = i;
        best_b = j;
      }
    }
  }
  return {best_a, best_b};
}

std::array<int, 2> MeshTopology::BisectElement(int elem_idx, int force_va, int force_vb)
{
  MFEM_ASSERT(elements[elem_idx].active,
              "Cannot bisect an inactive (already refined) element!");

  int la, lb;
  if (force_va >= 0 && force_vb >= 0)
  {
    // Find the local indices of the forced edge vertices.
    const auto &v = elements[elem_idx].vertices;
    la = static_cast<int>(std::find(v.begin(), v.end(), force_va) - v.begin());
    lb = static_cast<int>(std::find(v.begin(), v.end(), force_vb) - v.begin());
    MFEM_ASSERT(la < static_cast<int>(v.size()) && lb < static_cast<int>(v.size()),
                "Forced bisection edge vertices not found in element!");
  }
  else
  {
    std::tie(la, lb) = FindLongestEdge(elem_idx);
  }
  int va = elements[elem_idx].vertices[la];
  int vb = elements[elem_idx].vertices[lb];
  int attr = elements[elem_idx].attribute;
  int level = elements[elem_idx].level;

  // Get or create midpoint.
  int mid = GetMidpoint(va, vb);

  // Build child vertex lists by replacing one endpoint with the midpoint.
  auto verts0 = elements[elem_idx].vertices;
  auto verts1 = elements[elem_idx].vertices;
  verts0[lb] = mid;
  verts1[la] = mid;

  int child0_idx = static_cast<int>(elements.size());
  int child1_idx = child0_idx + 1;

  // Deactivate parent BEFORE push_back (push_back may reallocate and invalidate refs).
  elements[elem_idx].active = false;

  elements.push_back({std::move(verts0), attr, elem_idx, level + 1, true});
  elements.push_back({std::move(verts1), attr, elem_idx, level + 1, true});

  // Record history.
  history.push_back({elem_idx, {child0_idx, child1_idx}, va, vb, mid});

  return {child0_idx, child1_idx};
}

void MeshTopology::UpdateBdrElements(const std::vector<int> &refined_elements)
{
  // For each boundary element, check if volume element refinement split its edges.
  // Triangles and segments are handled iteratively (bisection may split one edge at a
  // time). Quads are split into 4 child quads when all 4 edge midpoints exist (they
  // always will after hex refinement since all 12 edges are split).
  bool changed = true;
  while (changed)
  {
    changed = false;
    std::vector<BdrElement> new_bdr;

    for (const auto &bdr : bdr_elements)
    {
      int nv = static_cast<int>(bdr.vertices.size());

      if (nv == 4)
      {
        // Quad boundary face: check if all 4 edge midpoints exist.
        // If so, compute face center and split into 4 child quads.
        int v0 = bdr.vertices[0], v1 = bdr.vertices[1];
        int v2 = bdr.vertices[2], v3 = bdr.vertices[3];
        auto k01 = EdgeKey(v0, v1), k12 = EdgeKey(v1, v2);
        auto k23 = EdgeKey(v2, v3), k30 = EdgeKey(v3, v0);
        auto it01 = edge_midpoints.find(k01), it12 = edge_midpoints.find(k12);
        auto it23 = edge_midpoints.find(k23), it30 = edge_midpoints.find(k30);

        if (it01 != edge_midpoints.end() && it12 != edge_midpoints.end() &&
            it23 != edge_midpoints.end() && it30 != edge_midpoints.end())
        {
          int e01 = it01->second, e12 = it12->second;
          int e23 = it23->second, e30 = it30->second;
          int fc = GetFaceCenter(v0, v1, v2, v3);

          new_bdr.push_back({{v0, e01, fc, e30}, bdr.attribute});
          new_bdr.push_back({{e01, v1, e12, fc}, bdr.attribute});
          new_bdr.push_back({{fc, e12, v2, e23}, bdr.attribute});
          new_bdr.push_back({{e30, fc, e23, v3}, bdr.attribute});
          changed = true;
        }
        else
        {
          new_bdr.push_back(bdr);
        }
      }
      else
      {
        // Triangle or segment boundary: split along the first edge with a midpoint.
        bool split = false;
        for (int i = 0; i < nv && !split; i++)
        {
          for (int j = i + 1; j < nv && !split; j++)
          {
            auto key = EdgeKey(bdr.vertices[i], bdr.vertices[j]);
            auto it = edge_midpoints.find(key);
            if (it != edge_midpoints.end())
            {
              int mid = it->second;
              if (nv == 3)
              {
                // Boundary triangle: split into 2 triangles along this edge.
                auto verts0 = bdr.vertices;
                auto verts1 = bdr.vertices;
                verts0[j] = mid;
                verts1[i] = mid;
                new_bdr.push_back({std::move(verts0), bdr.attribute});
                new_bdr.push_back({std::move(verts1), bdr.attribute});
              }
              else
              {
                // Boundary segment: split into 2 edges.
                new_bdr.push_back({{bdr.vertices[i], mid}, bdr.attribute});
                new_bdr.push_back({{mid, bdr.vertices[j]}, bdr.attribute});
              }
              split = true;
              changed = true;
            }
          }
        }

        if (!split)
        {
          new_bdr.push_back(bdr);
        }
      }
    }

    bdr_elements = std::move(new_bdr);
  }
}

MeshTopology MeshTopology::FromMFEM(const mfem::Mesh &mesh)
{
  MeshTopology topo;
  topo.dim = mesh.Dimension();
  topo.sdim = mesh.SpaceDimension();

  // Import vertices.
  int nv = mesh.GetNV();
  topo.vertices.resize(nv);
  for (int i = 0; i < nv; i++)
  {
    const double *v = mesh.GetVertex(i);
    for (int d = 0; d < topo.sdim; d++)
    {
      topo.vertices[i][d] = v[d];
    }
    for (int d = topo.sdim; d < 3; d++)
    {
      topo.vertices[i][d] = 0.0;
    }
  }

  // Import elements.
  int ne = mesh.GetNE();
  topo.elements.resize(ne);
  for (int i = 0; i < ne; i++)
  {
    mfem::Array<int> v;
    mesh.GetElementVertices(i, v);
    auto geom = mesh.GetElementGeometry(i);
    MFEM_VERIFY(geom == mfem::Geometry::TETRAHEDRON || geom == mfem::Geometry::TRIANGLE ||
                    geom == mfem::Geometry::CUBE || geom == mfem::Geometry::SQUARE,
                "MeshTopology only supports tet/tri/hex/quad elements, got geometry type "
                    << geom);
    topo.elements[i].vertices.assign(v.begin(), v.end());
    topo.elements[i].attribute = mesh.GetAttribute(i);
    topo.elements[i].parent = -1;
    topo.elements[i].level = 0;
    topo.elements[i].active = true;
  }

  // Import boundary elements.
  int nbe = mesh.GetNBE();
  topo.bdr_elements.resize(nbe);
  for (int i = 0; i < nbe; i++)
  {
    mfem::Array<int> v;
    mesh.GetBdrElementVertices(i, v);
    topo.bdr_elements[i].vertices.assign(v.begin(), v.end());
    topo.bdr_elements[i].attribute = mesh.GetBdrAttribute(i);
  }

  return topo;
}

std::unique_ptr<mfem::Mesh> MeshTopology::ToMFEM() const
{
  // Collect active elements.
  std::vector<int> active_indices;
  for (int i = 0; i < static_cast<int>(elements.size()); i++)
  {
    if (elements[i].active)
    {
      active_indices.push_back(i);
    }
  }

  int ne = static_cast<int>(active_indices.size());
  int nbe = static_cast<int>(bdr_elements.size());

  auto mesh =
      std::make_unique<mfem::Mesh>(dim, static_cast<int>(vertices.size()), ne, nbe, sdim);

  // Add vertices.
  for (int i = 0; i < static_cast<int>(vertices.size()); i++)
  {
    mesh->AddVertex(vertices[i].data());
  }

  // Determine element type from vertex count.
  auto GetElemType = [this](int nv) -> mfem::Element::Type
  {
    if (dim == 3)
    {
      return (nv == 8) ? mfem::Element::HEXAHEDRON : mfem::Element::TETRAHEDRON;
    }
    else
    {
      return (nv == 4) ? mfem::Element::QUADRILATERAL : mfem::Element::TRIANGLE;
    }
  };

  auto GetBdrElemType = [this](int nv) -> mfem::Element::Type
  {
    if (dim == 3)
    {
      return (nv == 4) ? mfem::Element::QUADRILATERAL : mfem::Element::TRIANGLE;
    }
    else
    {
      return mfem::Element::SEGMENT;
    }
  };

  // Add active elements.
  for (int idx : active_indices)
  {
    const auto &elem = elements[idx];
    mfem::Element *e = mesh->NewElement(GetElemType(elem.vertices.size()));
    e->SetVertices(elem.vertices.data());
    e->SetAttribute(elem.attribute);
    mesh->AddElement(e);
  }

  // Add boundary elements.
  for (int i = 0; i < nbe; i++)
  {
    const auto &bdr = bdr_elements[i];
    mfem::Element *e = mesh->NewElement(GetBdrElemType(bdr.vertices.size()));
    e->SetVertices(bdr.vertices.data());
    e->SetAttribute(bdr.attribute);
    mesh->AddBdrElement(e);
  }

  mesh->FinalizeTopology();
  mesh->Finalize(true, true);

  return mesh;
}

void MeshTopology::UniformRefinement()
{
  // Collect all active element indices first (they'll be deactivated during refinement).
  std::vector<int> to_refine;
  for (int i = 0; i < static_cast<int>(elements.size()); i++)
  {
    if (elements[i].active)
    {
      to_refine.push_back(i);
    }
  }

  // Clear caches for this refinement pass.
  edge_midpoints.clear();
  face_centers.clear();
  body_centers.clear();

  // Refine each element based on its type.
  for (int idx : to_refine)
  {
    int nv = static_cast<int>(elements[idx].vertices.size());
    if (nv == 8)
    {
      RefineHexElement(idx);
    }
    else if (nv == 4 && dim == 2)
    {
      RefineQuadElement(idx);
    }
    else
    {
      // Simplex element (tet or tri): bisect.
      BisectElement(idx);
    }
  }

  // Update boundary elements for any edges that were split.
  UpdateBdrElements(to_refine);
}

void MeshTopology::Refine(const std::vector<int> &marked_elements)
{
  edge_midpoints.clear();
  face_centers.clear();
  body_centers.clear();

  // Phase 1: Refine all marked elements.
  for (int idx : marked_elements)
  {
    if (!elements[idx].active)
    {
      continue;
    }
    int nv = static_cast<int>(elements[idx].vertices.size());
    if (nv == 8)
    {
      RefineHexElement(idx);
    }
    else if (nv == 4 && dim == 2)
    {
      RefineQuadElement(idx);
    }
    else
    {
      // Simplex element: bisect along longest edge.
      BisectElement(idx);
    }
  }

  // Phase 2: Closure — only needed for simplex elements.
  // Hex/quad refinement does not require closure (hanging nodes are acceptable
  // or the mesh is uniformly refined). Closure only applies to simplex meshes.
  // Scan active simplex elements for edges that have midpoints in the cache.
  bool changed = true;
  while (changed)
  {
    changed = false;
    int ne = static_cast<int>(elements.size());
    for (int i = 0; i < ne; i++)
    {
      if (!elements[i].active)
      {
        continue;
      }
      int nv = static_cast<int>(elements[i].vertices.size());
      // Only simplex elements need closure.
      if ((nv == 4 && dim == 3) || nv == 3)
      {
        for (int a = 0; a < nv; a++)
        {
          for (int b = a + 1; b < nv; b++)
          {
            int va = elements[i].vertices[a];
            int vb = elements[i].vertices[b];
            auto key = EdgeKey(va, vb);
            if (edge_midpoints.count(key))
            {
              // Hanging vertex found — bisect along this edge.
              elements.reserve(elements.size() + 2);
              BisectElement(i, va, vb);
              changed = true;
              goto next_element;
            }
          }
        }
      }
    next_element:;
    }
  }

  // Phase 3: Update boundary elements for all edges that were split.
  UpdateBdrElements(marked_elements);
}

void MeshTopology::RefineDistributed(MPI_Comm comm, const std::vector<int> &local_marks,
                                     const std::vector<int> &local_to_global_elem)
{
  // Map local marks to global element indices.
  std::vector<int> global_marks;
  global_marks.reserve(local_marks.size());
  for (int local_idx : local_marks)
  {
    MFEM_ASSERT(local_idx >= 0 && local_idx < static_cast<int>(local_to_global_elem.size()),
                "Local element index out of range in RefineDistributed!");
    global_marks.push_back(local_to_global_elem[local_idx]);
  }

  // Allgather: collect all global marks from all ranks.
  int local_count = static_cast<int>(global_marks.size());
  int world_size;
  MPI_Comm_size(comm, &world_size);

  std::vector<int> recv_counts(world_size);
  MPI_Allgather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

  std::vector<int> displacements(world_size);
  int total_count = 0;
  for (int r = 0; r < world_size; r++)
  {
    displacements[r] = total_count;
    total_count += recv_counts[r];
  }

  std::vector<int> all_marks(total_count);
  MPI_Allgatherv(global_marks.data(), local_count, MPI_INT, all_marks.data(),
                 recv_counts.data(), displacements.data(), MPI_INT, comm);

  // Deduplicate (elements could theoretically be marked by multiple ranks if they're
  // in ghost layers, though in practice each element is owned by exactly one rank).
  std::sort(all_marks.begin(), all_marks.end());
  all_marks.erase(std::unique(all_marks.begin(), all_marks.end()), all_marks.end());

  // Run the global refinement with closure. Same on all ranks → deterministic result.
  Refine(all_marks);
}

void MeshTopology::Coarsen(const std::vector<double> &elem_error, double threshold)
{
  // Walk history backwards. For each refinement record, if all children are active
  // and their combined error is below the threshold, merge them back into the parent.
  for (int h = static_cast<int>(history.size()) - 1; h >= 0; h--)
  {
    const auto &rec = history[h];

    // Check if all children are active.
    bool all_active = true;
    for (int child_idx : rec.children)
    {
      if (!elements[child_idx].active)
      {
        all_active = false;
        break;
      }
    }
    if (!all_active)
    {
      continue;  // One or more children have been further refined.
    }

    // Check error: all children must have low combined error.
    double combined_error = 0.0;
    for (int child_idx : rec.children)
    {
      if (child_idx < static_cast<int>(elem_error.size()))
      {
        combined_error += elem_error[child_idx];
      }
    }

    if (combined_error < threshold)
    {
      // Merge: deactivate children, reactivate parent.
      for (int child_idx : rec.children)
      {
        elements[child_idx].active = false;
      }
      elements[rec.parent].active = true;

      // Note: The midpoint/face center/body center vertices stay in the vertex list
      // (they're just unused). Boundary elements are not updated here — a full
      // boundary rebuild would be needed for production use.
    }
  }
}

int MeshTopology::GetNE() const
{
  return GetActiveElementCount();
}

int MeshTopology::GetActiveElementCount() const
{
  int count = 0;
  for (const auto &elem : elements)
  {
    if (elem.active)
    {
      count++;
    }
  }
  return count;
}

}  // namespace palace
