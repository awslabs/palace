// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_POINT_FIELD_EVALUATOR_HPP
#define PALACE_FEM_POINT_FIELD_EVALUATOR_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class BoundaryDerivedFieldBundle;
class BoundaryPhysicalTraceCache;
class DomainPointFieldEvaluator;
class FaceSamplingPlan;
class GridFunction;
class MaterialOperator;
class Mesh;
class SurfaceFunctional;

// Non-reducing libCEED evaluator for visualization point fields. Unlike reduction
// functionals, this produces one value per local visualization point and is intended to
// feed visualization point data. Domain buffers are point-major; boundary buffers are
// component-major. The mesh entity location selects domain elements or boundary elements;
// the operation remains pointwise, not an MPI-reduced functional.
class PointFieldEvaluator
{
public:
  using Kind = PointFieldKind;

private:
  MeshEntityType location;
  Kind kind;
  bool valid = false;
  std::unique_ptr<DomainPointFieldEvaluator> domain_eval;
  std::unique_ptr<SurfaceFunctional> boundary_eval;
  std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache;
  // Optional 3D E/B fused boundary bundle. It exposes one packed field slice per
  // callback and replaces independent trace-consumer applies for its supported kinds.
  std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle;
  const mfem::ParFiniteElementSpace *trace_fespace_a = nullptr;
  const mfem::ParFiniteElementSpace *trace_fespace_b = nullptr;

public:
  // Domain pointwise evaluator for U_e, U_m, or S at VTU visualization points. The
  // target space is retained only to preserve the existing GridFunction export path.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                      const mfem::ParFiniteElementSpace *nd_fespace,
                      const mfem::ParFiniteElementSpace *rt_fespace,
                      const mfem::ParFiniteElementSpace &target_fespace, double scaling,
                      bool build_gridfunction = true, bool build_buffer = true);

  // Boundary field evaluator for E or B at boundary VTU visualization points.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &fespace, int lod,
                      std::shared_ptr<const FaceSamplingPlan> sampling_plan = nullptr,
                      std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr,
                      std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle = nullptr);

  // Boundary material-dependent evaluator for Q_s, J_s, U_e, or U_m.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &fespace,
                      const MaterialOperator &mat_op, int lod, double scaling,
                      std::shared_ptr<const FaceSamplingPlan> sampling_plan = nullptr,
                      std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr,
                      std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle = nullptr);

  // Boundary Poynting vector evaluator.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &nd_fespace,
                      const mfem::ParFiniteElementSpace &rt_fespace,
                      const MaterialOperator &mat_op, int lod, double scaling,
                      std::shared_ptr<const FaceSamplingPlan> sampling_plan = nullptr,
                      std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache = nullptr,
                      std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle = nullptr);
  ~PointFieldEvaluator();

  PointFieldEvaluator(const PointFieldEvaluator &) = delete;
  PointFieldEvaluator &operator=(const PointFieldEvaluator &) = delete;

  bool IsValid() const;
  MeshEntityType Location() const { return location; }
  Kind GetKind() const { return kind; }

  int BufferSize() const;
  int BufferNumComp() const;
  const std::vector<int> &BufferBases() const;

  // Non-null only for boundary evaluators. This exposes identity for focused setup
  // tests while keeping evaluation resources immutable and shared.
  const FaceSamplingPlan *SamplingPlan() const;
  const BoundaryPhysicalTraceCache *TraceCache() const { return trace_cache.get(); }

  // Domain-only compatibility path for grid-function output.
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;

  // Fill the visualization buffer. Domain buffers are point-major; boundary buffers are
  // component-major. The Vector overload is for already split
  // real/imaginary fields; the GridFunction overload accumulates real and imaginary
  // contributions for quadratic single-field quantities. For POYNTING pass E and B to
  // the two-field overload.
  void EvalBuffer(const Vector &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction *E, const GridFunction *B, Vector &buffer) const;
};

}  // namespace palace

#endif  // PALACE_FEM_POINT_FIELD_EVALUATOR_HPP
