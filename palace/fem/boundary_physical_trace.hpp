// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BOUNDARY_PHYSICAL_TRACE_HPP
#define PALACE_FEM_BOUNDARY_PHYSICAL_TRACE_HPP

#include <atomic>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/ceed_group_operator.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FaceSamplingPlan;
class Mesh;
class FaceNbrFieldExchange;
class SurfaceFunctional;

// Device-resident, side-separated physical traces for one fixed boundary sampling plan.
//
// The two trace buffers for an entry are stored in canonical mapped-point order and in
// component-major layout.  They are deliberately distinct from FaceNbrFieldExchange's
// point-major wire payload: an exchange is used only while materializing a ghost-side
// trace, never as the storage consumed by a boundary point evaluator.
class BoundaryPhysicalTraceCache
{
public:
  struct Traces
  {
    const Vector *side_a = nullptr;
    const Vector *side_b = nullptr;
    int num_comp = 0;
  };

  BoundaryPhysicalTraceCache(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                             std::shared_ptr<const FaceSamplingPlan> sampling_plan);
  ~BoundaryPhysicalTraceCache();

  BoundaryPhysicalTraceCache(const BoundaryPhysicalTraceCache &) = delete;
  BoundaryPhysicalTraceCache &operator=(const BoundaryPhysicalTraceCache &) = delete;

  // Start a new ParaView save generation.  A source Vector has no mutation generation
  // of its own, so this explicit boundary is what prevents same-address solution vectors
  // from reusing a trace from a previous save.
  void Invalidate();

  // Fail fast if AMR/rebalance or nodal-coordinate updates invalidated the immutable
  // sampling/evaluator state.
  void VerifyCurrent() const;

  // Materialize (or reuse) physical traces for u.  Collective when the plan contains
  // face-neighbor routes.  The source-space identity is part of the cache key.
  Traces Get(const mfem::ParFiniteElementSpace &fespace, const Vector &u);

  const FaceSamplingPlan &SamplingPlan() const;
  int NumPoints() const;
  int NumComponents(const mfem::ParFiniteElementSpace &fespace) const;
  static bool Supports(const mfem::ParFiniteElementSpace &fespace);

  // Setup/evaluation profile counters.  They are intentionally process-local, matching
  // FaceSamplingPlan's counters, and are printed only when explicitly requested.
  static long long TraceGroupCount();
  static long long OperatorApplyCount();
  static long long CacheHitCount();
  static long long EvaluatedPointCount();
  static long long EvaluatedByteCount();

  // Exact setup counters for selected NC coarse-face leaf unions. A route denotes one
  // canonical leaf trace point scattered from a union value; public leaf tuples remain
  // duplicated even when their source coordinate is shared.
  static long long LeafPointCount();
  static long long UnionPointCount();
  static long long DuplicatePointCount();
  static long long LocalUnionGroupCount();
  static long long LocalUnionOperatorCount();
  static long long GhostUnionGroupCount();
  static long long RequestCountBefore();
  static long long RequestCountAfter();
  static long long UnionRouteCount();
  static long long RoutingOperatorCount();

private:
  struct Entry
  {
    struct RoutingOperator
    {
      Ceed ceed = nullptr;
      CeedOperator op = nullptr;
      CeedVector in_vec = nullptr, out_vec = nullptr;
      bool imported_input = false;
      int side = -1;
    };

    struct UnionValues
    {
      std::vector<fem::CeedGroupOperator> local_evaluators;
      Vector local_values;
      std::unique_ptr<FaceNbrFieldExchange> ghost_exchange;
      std::vector<RoutingOperator> routes;
    };

    struct SourceValues
    {
      Vector side_a, side_b;
      long generation = -1;
    };

    const mfem::ParFiniteElementSpace *fespace = nullptr;
    PointFieldKind kind = PointFieldKind::FIELD_E;
    int num_comp = 0;
    std::unique_ptr<SurfaceFunctional> side_a_eval, side_b_eval;
    UnionValues union_values;
    // Keep separately materialized real/imaginary (or otherwise distinct) source
    // vectors. Invalidate clears this map at the save boundary, bounding its lifetime.
    std::map<const Vector *, SourceValues> source_values;
  };

  const Mesh *mesh;
  std::vector<int> marker;
  std::shared_ptr<const FaceSamplingPlan> sampling_plan;
  long generation = 0;
  std::vector<Entry> entries;

  Entry &GetEntry(const mfem::ParFiniteElementSpace &fespace);
  void BuildUnionRoutes(Entry &entry);
  void ApplyUnionRoutes(Entry &entry, const Vector &u, Entry::SourceValues &values) const;
  static void DestroyUnionRoutes(Entry::UnionValues &union_values);
  static PointFieldKind SourceKind(const mfem::ParFiniteElementSpace &fespace);

  static std::atomic<long long> trace_groups;
  static std::atomic<long long> operator_applies;
  static std::atomic<long long> cache_hits;
  static std::atomic<long long> evaluated_points;
  static std::atomic<long long> evaluated_bytes;
  static std::atomic<long long> leaf_points;
  static std::atomic<long long> union_points;
  static std::atomic<long long> duplicate_points;
  static std::atomic<long long> local_union_groups;
  static std::atomic<long long> local_union_operators;
  static std::atomic<long long> ghost_union_groups;
  static std::atomic<long long> requests_before;
  static std::atomic<long long> requests_after;
  static std::atomic<long long> union_routes;
  static std::atomic<long long> routing_operators;
};

}  // namespace palace

#endif  // PALACE_FEM_BOUNDARY_PHYSICAL_TRACE_HPP
