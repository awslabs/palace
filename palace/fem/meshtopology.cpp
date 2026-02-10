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
  // For each boundary element, check if volume element bisection split any of its edges.
  // A boundary triangle may need multiple splits if multiple edges have midpoints.
  // We iterate until no more splits are needed.
  bool changed = true;
  while (changed)
  {
    changed = false;
    std::vector<BdrElement> new_bdr;

    for (const auto &bdr : bdr_elements)
    {
      int nv = static_cast<int>(bdr.vertices.size());
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
            if (dim == 3)
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
              // Boundary edge: split into 2 edges.
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
    MFEM_VERIFY(geom == mfem::Geometry::TETRAHEDRON || geom == mfem::Geometry::TRIANGLE,
                "MeshTopology only supports simplex elements (tet/tri), got geometry type "
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
  int nv_per_elem = (dim == 3) ? 4 : 3;
  int nv_per_bdr = (dim == 3) ? 3 : 2;
  auto elem_type = (dim == 3) ? mfem::Element::TETRAHEDRON : mfem::Element::TRIANGLE;
  auto bdr_type = (dim == 3) ? mfem::Element::TRIANGLE : mfem::Element::SEGMENT;

  auto mesh =
      std::make_unique<mfem::Mesh>(dim, static_cast<int>(vertices.size()), ne, nbe, sdim);

  // Add vertices.
  for (int i = 0; i < static_cast<int>(vertices.size()); i++)
  {
    mesh->AddVertex(vertices[i].data());
  }

  // Add active elements.
  for (int idx : active_indices)
  {
    const auto &elem = elements[idx];
    mfem::Element *e = mesh->NewElement(elem_type);
    e->SetVertices(elem.vertices.data());
    e->SetAttribute(elem.attribute);
    mesh->AddElement(e);
  }

  // Add boundary elements.
  for (int i = 0; i < nbe; i++)
  {
    mfem::Element *e = mesh->NewElement(bdr_type);
    e->SetVertices(bdr_elements[i].vertices.data());
    e->SetAttribute(bdr_elements[i].attribute);
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

  // Clear edge midpoint cache for this refinement pass.
  edge_midpoints.clear();

  // Reserve space for children to avoid repeated reallocation.
  elements.reserve(elements.size() + 2 * to_refine.size());

  // Bisect each element.
  for (int idx : to_refine)
  {
    BisectElement(idx);
  }

  // Update boundary elements for any edges that were split.
  UpdateBdrElements(to_refine);
}

void MeshTopology::Refine(const std::vector<int> &marked_elements)
{
  edge_midpoints.clear();

  // Phase 1: Bisect all marked elements along their longest edges.
  for (int idx : marked_elements)
  {
    if (elements[idx].active)
    {
      BisectElement(idx);
    }
  }

  // Phase 2: Closure — resolve hanging vertices to maintain conformity.
  // Scan active elements for edges that have midpoints in the cache. If found,
  // that element has a hanging vertex and must be bisected along that edge.
  // Repeat until no hanging vertices remain.
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
  // Walk history backwards. For each refinement record, if both children are active
  // and their combined error is below the threshold, merge them back into the parent.
  for (int h = static_cast<int>(history.size()) - 1; h >= 0; h--)
  {
    const auto &rec = history[h];
    auto &child0 = elements[rec.children[0]];
    auto &child1 = elements[rec.children[1]];

    if (!child0.active || !child1.active)
    {
      continue;  // One or both children have been further refined.
    }

    // Check error: both children must have low error.
    double combined_error = 0.0;
    if (rec.children[0] < static_cast<int>(elem_error.size()))
    {
      combined_error += elem_error[rec.children[0]];
    }
    if (rec.children[1] < static_cast<int>(elem_error.size()))
    {
      combined_error += elem_error[rec.children[1]];
    }

    if (combined_error < threshold)
    {
      // Merge: deactivate children, reactivate parent.
      child0.active = false;
      child1.active = false;
      elements[rec.parent].active = true;

      // Note: The midpoint vertex stays in the vertex list (it's just unused).
      // Boundary elements are not updated here — a full boundary rebuild would be
      // needed for production use.
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
