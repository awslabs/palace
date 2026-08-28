// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FACE_SAMPLING_PLAN_HPP
#define PALACE_FEM_FACE_SAMPLING_PLAN_HPP

#include <atomic>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class Mesh;

// Immutable setup data shared by compatible boundary point-field evaluators. The plan
// deliberately contains only mesh/sampling information: source slots, material kernels,
// and face-neighbor source masks remain properties of each consuming evaluator.
class FaceSamplingPlan
{
public:
  struct Entry
  {
    int bdr_elem = -1;
    int attr = 0;
    int elem_a = -1, elem_b = -1;
    int face_nbr = -1;
    int ghost_attr = 0;
    bool flip = false;
    bool ghost_b = false;
    mfem::Geometry::Type bdr_geom = mfem::Geometry::INVALID;
    mfem::Geometry::Type vol_geom_a = mfem::Geometry::INVALID;
    mfem::Geometry::Type vol_geom_b = mfem::Geometry::INVALID;

    // Canonical mapped volume reference clouds. canonical_to_original routes a canonical
    // q-index back to the established boundary-refinement-rule output ordering.
    std::vector<mfem::IntegrationPoint> points_a, points_b;
    std::vector<int> canonical_to_original_a, canonical_to_original_b;
    std::vector<long long> point_key_a, point_key_b;

    // True when this side is the coarse owner of a nonconforming leaf face. These
    // flags are setup facts only; union_group_* is assigned only if several leaves of
    // that owner actually share a mapped reference coordinate.
    bool nc_coarse_a = false, nc_coarse_b = false;
    // The containing volume reference-face identity, used only to avoid combining
    // unrelated faces of the same NC element which meet at a vertex or edge.
    int reference_face_a = -1, reference_face_b = -1;
    int union_group_a = -1, union_group_b = -1;

    // Setup-only boundary Jacobians in original boundary-rule q order, with
    // [component + space_dimension * reference_dimension][q] layout. Consumers route
    // them with canonical_to_original_*.
    std::vector<double> face_jacobians;
  };

  FaceSamplingPlan(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker, int lod);

  // A selected nonconforming coarse owner is evaluated at these coordinates once;
  // routes then duplicate the result into each listed leaf's canonical trace tuple.
  // The coordinates are normalized and deduplicated by exact coordinate bits only.
  // In particular, quadrature weights are deliberately not part of this union.
  struct UnionRoute
  {
    int entry = -1;
    std::vector<int> canonical_to_union;
  };
  struct UnionGroup
  {
    bool ghost = false;
    int owner = -1;  // Local element, or face-neighbor element for a ghost owner.
    int side = -1;
    mfem::Geometry::Type geometry = mfem::Geometry::INVALID;
    std::vector<mfem::IntegrationPoint> points;
    std::vector<UnionRoute> routes;
  };

  FaceSamplingPlan(const FaceSamplingPlan &) = delete;
  FaceSamplingPlan &operator=(const FaceSamplingPlan &) = delete;

  bool Matches(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker, int lod) const;
  int LOD() const { return lod; }
  int NumPoints() const { return num_points; }
  const std::vector<int> &BufferBases() const { return buffer_bases; }
  const std::vector<Entry> &Entries() const { return entries; }
  const std::vector<UnionGroup> &UnionGroups() const { return union_groups; }
  bool HasGhostUnion() const { return has_ghost_union; }

  // Counters are process-local and intentionally inexpensive. They are useful in focused
  // tests and, with PALACE_FACE_SAMPLING_PLAN_PROFILE set, in setup profiling output.
  static long long BuildCount();
  static long long ReuseCount();
  static void NoteReuse();

private:
  const Mesh *mesh;
  long mesh_sequence;
  long nodes_sequence;
  int lod;
  int num_points = 0;
  std::vector<int> marker;
  std::vector<int> buffer_bases;
  std::vector<Entry> entries;
  std::vector<UnionGroup> union_groups;
  bool has_ghost_union = false;

  static std::atomic<long long> builds;
  static std::atomic<long long> reuses;
};

}  // namespace palace

#endif  // PALACE_FEM_FACE_SAMPLING_PLAN_HPP
