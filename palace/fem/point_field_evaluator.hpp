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
// functionals, this produces one value per local visualization point and is intended to
// feed ParaView point data. Domain direct-VTU buffers are point-major; boundary buffers
// use the existing component-major writer path. The mesh entity location selects domain
// elements or boundary elements; the operation remains pointwise, not an MPI-reduced
// functional.
class PointFieldEvaluator
{
public:
  using Kind = PointFieldKind;

private:
  MeshEntityType location;
  Kind kind;
  std::unique_ptr<DomainPointFieldEvaluator> domain_eval;
  std::unique_ptr<SurfaceFunctional> boundary_eval;

  static int NumComponents(Kind kind);

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
  int BufferNumComp() const;
  const std::vector<int> &BufferBases() const;

  // Domain-only compatibility path for grid-function output.
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;

  // Fill the visualization buffer. Domain buffers are point-major for direct VTU output;
  // boundary buffers are component-major. The Vector overload is for already split
  // real/imaginary fields; the GridFunction overload accumulates real and imaginary
  // contributions for quadratic single-field quantities. For POYNTING pass E and B to
  // the two-field overload.
  void EvalBuffer(const Vector &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction &u, Vector &buffer) const;
  void EvalBuffer(const GridFunction *E, const GridFunction *B, Vector &buffer) const;
};

}  // namespace palace

#endif  // PALACE_FEM_POINT_FIELD_EVALUATOR_HPP
