// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_BOUNDARY_DERIVED_FIELD_BUNDLE_HPP
#define PALACE_FEM_BOUNDARY_DERIVED_FIELD_BUNDLE_HPP

#include <atomic>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class BoundaryPhysicalTraceCache;
class FaceSamplingPlan;
class GridFunction;
class MaterialOperator;
class Mesh;

// Fused, device-resident boundary postprocessing for the 3D E/B output family.
//
// BoundaryPhysicalTraceCache owns physical, side-separated E and B traces.  This bundle
// consumes those traces once per canonical sampling-plan route and produces a packed
// component-major result for all compatible derived fields. For registered complex E/B,
// one route application binds both physical phases, stores phase-resolved linear slices,
// and stores the combined quadratic slices. The packed result is keyed by source-vector
// identity and the explicit save generation. Direct VTU output can expose non-owning
// slice views while PointFieldEvaluator retains Copy() for non-bundle callers.
class BoundaryDerivedFieldBundle
{
public:
  enum class Slice
  {
    FIELD_E,
    FIELD_B,
    FLUX_Q,
    CURRENT_J,
    ENERGY_E,
    ENERGY_M,
    POYNTING
  };

  BoundaryDerivedFieldBundle(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                             const MaterialOperator &mat_op,
                             const mfem::ParFiniteElementSpace &nd_fespace,
                             const mfem::ParFiniteElementSpace &rt_fespace,
                             std::shared_ptr<const FaceSamplingPlan> sampling_plan,
                             std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache,
                             const GridFunction &E, const GridFunction &B,
                             double electric_scaling, double magnetic_scaling);
  ~BoundaryDerivedFieldBundle();

  BoundaryDerivedFieldBundle(const BoundaryDerivedFieldBundle &) = delete;
  BoundaryDerivedFieldBundle &operator=(const BoundaryDerivedFieldBundle &) = delete;

  bool IsValid() const { return valid; }
  bool Supports(PointFieldKind kind, const mfem::ParFiniteElementSpace *fespace_a,
                const mfem::ParFiniteElementSpace *fespace_b = nullptr) const;
  const FaceSamplingPlan *SamplingPlan() const;
  int BufferSize(PointFieldKind kind) const;
  int BufferNumComp(PointFieldKind kind) const;
  const std::vector<int> &BufferBases() const;

  // Expose the registered E/B packed result only to the narrow boundary VTU batching
  // path. GetRegisteredPacked() evaluates (or returns the cached result for) the current
  // save generation. RegisteredSlice() gives exact scalar offsets/sizes for real or
  // imaginary linear fields and the shared combined quadratic fields.
  struct PackedSlice
  {
    int offset = 0;
    int size = 0;
    int num_comp = 0;
  };
  const Vector &GetRegisteredPacked() const;
  int RegisteredPackedSize() const;
  PackedSlice RegisteredSlice(PointFieldKind kind, bool imaginary_phase = false) const;

  // Advance the save generation. The trace cache must be invalidated at the same save
  // boundary by its owner; keeping both explicit makes reuse of solver vector addresses
  // unambiguous.
  void Invalidate();

  // Copy or accumulate a requested packed slice. The single-source overload identifies
  // a real/imaginary phase from either registered E or B source and obtains its partner.
  void Copy(PointFieldKind kind, const Vector &source, Vector &buffer,
            bool zero = true) const;
  void Copy(PointFieldKind kind, const Vector &E, const Vector &B, Vector &buffer,
            bool zero = true) const;

  static long long BundleCount();
  static long long FieldCount();
  static long long QFunctionApplyCount();
  // Number of logical per-route phases evaluated. A batched complex application counts
  // both its real and imaginary phases even though it launches one operator.
  static long long PhaseCount();
  // Per-route operator launches saved by evaluating the two complex phases together.
  static long long AvoidedPhaseOperatorApplyCount();
  static long long AvoidedIndependentApplyCount();
  static long long PackedByteCount();
  static long long CacheHitCount();

  // True only for the registered complex E/B pair. This lets PointFieldEvaluator retain
  // the independent-source fallback while copying its combined quadratic slices once.
  bool BatchesComplexPhases(const Vector &E_source, const Vector &B_source) const;
  bool BatchesComplexPhases(const Vector &source) const;
  void PrintProfile() const;

private:
  struct CachedPhase
  {
    Vector packed;
    long generation = -1;
    bool batched_complex = false;
  };

  const Mesh *mesh = nullptr;
  const MaterialOperator *mat_op = nullptr;
  const mfem::ParFiniteElementSpace *nd_fespace = nullptr;
  const mfem::ParFiniteElementSpace *rt_fespace = nullptr;
  const GridFunction *E = nullptr;
  const GridFunction *B = nullptr;
  std::shared_ptr<const FaceSamplingPlan> sampling_plan;
  std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache;
  bool valid = false;
  long generation = 0;
  // Kept opaque in the header to avoid exposing libCEED implementation details at every
  // PointFieldEvaluator include site.
  struct Impl;
  std::unique_ptr<Impl> impl;
  mutable std::map<std::pair<const Vector *, const Vector *>, CachedPhase> phases;

  const CachedPhase &Get(const Vector &E_source, const Vector &B_source) const;
  std::pair<const Vector *, const Vector *> PartnerSources(const Vector &source) const;
  static Slice ToSlice(PointFieldKind kind);
  static int SliceOffset(Slice slice, int num_points, bool batched_complex,
                         bool imaginary_phase = false);
  static int SliceNumComp(Slice slice);
  bool IsImaginaryPhase(const Vector &E_source, const Vector &B_source) const;

  static std::atomic<long long> bundle_count;
  static std::atomic<long long> field_count;
  static std::atomic<long long> qfunction_applies;
  static std::atomic<long long> phase_count;
  static std::atomic<long long> avoided_phase_operator_applies;
  static std::atomic<long long> avoided_independent_applies;
  static std::atomic<long long> packed_bytes;
  static std::atomic<long long> cache_hits;
};

}  // namespace palace

#endif  // PALACE_FEM_BOUNDARY_DERIVED_FIELD_BUNDLE_HPP
