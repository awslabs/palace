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

class DomainPointFieldEvaluator;
class GridFunction;
class MaterialOperator;
class Mesh;
class SurfaceFunctional;

// Non-reducing libCEED evaluator for visualization point fields. Unlike reduction
// functionals, this produces one value per local visualization point (component-major in
// the output buffer) and is intended to feed ParaView point data. The mesh entity
// location selects domain elements or boundary elements; the operation remains
// pointwise, not an MPI-reduced functional.
class PointFieldEvaluator
{
public:
  using Kind = PointFieldKind;

private:
  enum class BoundaryEvaluatorType
  {
    NONE,
    FIELD,
    MATERIAL,
    POYNTING
  };

  MeshEntityType location;
  Kind kind;
  bool valid = false;
  std::unique_ptr<DomainPointFieldEvaluator> domain_eval;
  mutable std::unique_ptr<SurfaceFunctional> boundary_eval;

  // Boundary evaluators may be released after metadata extraction/evaluation and rebuilt
  // on demand to avoid retaining many large libCEED AtPoints operators simultaneously.
  BoundaryEvaluatorType boundary_type = BoundaryEvaluatorType::NONE;
  const Mesh *boundary_mesh = nullptr;
  mfem::Array<int> boundary_marker;
  const mfem::ParFiniteElementSpace *boundary_fespace = nullptr;
  const mfem::ParFiniteElementSpace *boundary_nd_fespace = nullptr;
  const mfem::ParFiniteElementSpace *boundary_rt_fespace = nullptr;
  const MaterialOperator *boundary_mat_op = nullptr;
  int boundary_lod = 0;
  double boundary_scaling = 1.0;
  int boundary_buffer_size = 0;
  std::vector<int> boundary_buffer_bases;
  mutable int retain_boundary_eval_count = 0;

  static int NumComponents(Kind kind);
  void EnsureBoundaryEvaluator() const;
  void CacheBoundaryMetadata();
  bool ReleaseBoundaryEvaluatorAfterUse() const;
  void MaybeReleaseBoundaryEvaluator() const;

public:
  // Domain pointwise evaluator for U_e, U_m, or S at VTU visualization points. The
  // target space is retained only to preserve the existing GridFunction export path.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                      const mfem::ParFiniteElementSpace *nd_fespace,
                      const mfem::ParFiniteElementSpace *rt_fespace,
                      const mfem::ParFiniteElementSpace &target_fespace, double scaling);

  // Boundary field evaluator for E or B at boundary VTU visualization points.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &fespace, int lod);

  // Boundary material-dependent evaluator for Q_s, J_s, U_e, or U_m.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &fespace,
                      const MaterialOperator &mat_op, int lod, double scaling);

  // Boundary Poynting vector evaluator.
  PointFieldEvaluator(Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
                      const mfem::ParFiniteElementSpace &nd_fespace,
                      const mfem::ParFiniteElementSpace &rt_fespace,
                      const MaterialOperator &mat_op, int lod, double scaling);
  ~PointFieldEvaluator();

  PointFieldEvaluator(const PointFieldEvaluator &) = delete;
  PointFieldEvaluator &operator=(const PointFieldEvaluator &) = delete;

  bool IsValid() const;
  MeshEntityType Location() const { return location; }
  Kind GetKind() const { return kind; }

  int BufferSize() const;
  int BufferNumComp() const { return NumComponents(kind); }
  const std::vector<int> &BufferBases() const;

  // Keep the next boundary evaluator instance alive for a small fixed number of
  // following EvalBuffer calls. This lets ParaView reuse the same libCEED operator for
  // adjacent mode/component fields without retaining it across unrelated field kinds.
  void RetainBoundaryEvaluatorFor(int count) const { retain_boundary_eval_count = count; }

  // Domain-only compatibility path for grid-function output.
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;

  // Fill the component-major visualization buffer. The Vector overload is for already
  // split real/imaginary boundary fields; the GridFunction overload accumulates real and
  // imaginary contributions for quadratic single-field quantities. For POYNTING pass E
  // and B to the two-field overload.
  void EvalBuffer(const Vector &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction *E, const GridFunction *B, Vector &buffer) const;
};

}  // namespace palace

#endif  // PALACE_FEM_POINT_FIELD_EVALUATOR_HPP
