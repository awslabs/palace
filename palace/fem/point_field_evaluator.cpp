// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/point_field_evaluator.hpp"

#include "fem/domain_point_field_evaluator.hpp"
#include "fem/gridfunction.hpp"
#include "fem/output_functionals.hpp"

namespace palace
{

namespace
{

DomainPointFieldEvaluator::Kind ToDomainKind(PointFieldEvaluator::Kind kind)
{
  switch (kind)
  {
    case PointFieldEvaluator::Kind::ENERGY_E:
      return DomainPointFieldEvaluator::Kind::ENERGY_E;
    case PointFieldEvaluator::Kind::ENERGY_M:
      return DomainPointFieldEvaluator::Kind::ENERGY_M;
    case PointFieldEvaluator::Kind::POYNTING:
      return DomainPointFieldEvaluator::Kind::POYNTING;
    default:
      MFEM_ABORT("Unsupported domain point field kind!");
  }
}


}  // namespace

int PointFieldEvaluator::NumComponents(Kind kind)
{
  return (kind == Kind::FLUX_Q || kind == Kind::ENERGY_E || kind == Kind::ENERGY_M) ? 1
                                                                                     : 3;
}

PointFieldEvaluator::PointFieldEvaluator(
    Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
    const mfem::ParFiniteElementSpace *nd_fespace,
    const mfem::ParFiniteElementSpace *rt_fespace,
    const mfem::ParFiniteElementSpace &target_fespace, double scaling)
  : location(MeshEntityType::Domain), kind(kind)
{
  MFEM_VERIFY(kind == Kind::ENERGY_E || kind == Kind::ENERGY_M || kind == Kind::POYNTING,
              "Unsupported domain point field kind!");
  domain_eval = std::make_unique<DomainPointFieldEvaluator>(
      ToDomainKind(kind), mesh, mat_op, nd_fespace, rt_fespace, target_fespace, scaling);
  valid = domain_eval->IsValid();
}

PointFieldEvaluator::PointFieldEvaluator(Kind kind, const Mesh &mesh,
                                         const mfem::Array<int> &bdr_attr_marker,
                                         const mfem::ParFiniteElementSpace &fespace,
                                         int lod)
  : location(MeshEntityType::Boundary), kind(kind),
    boundary_type(BoundaryEvaluatorType::FIELD), boundary_mesh(&mesh),
    boundary_marker(bdr_attr_marker), boundary_fespace(&fespace), boundary_lod(lod)
{
  MFEM_VERIFY(kind == Kind::FIELD_E || kind == Kind::FIELD_B,
              "Unsupported boundary field point kind!");
  EnsureBoundaryEvaluator();
  CacheBoundaryMetadata();
  boundary_eval.reset();
}

PointFieldEvaluator::PointFieldEvaluator(Kind kind, const Mesh &mesh,
                                         const mfem::Array<int> &bdr_attr_marker,
                                         const mfem::ParFiniteElementSpace &fespace,
                                         const MaterialOperator &mat_op, int lod,
                                         double scaling)
  : location(MeshEntityType::Boundary), kind(kind),
    boundary_type(BoundaryEvaluatorType::MATERIAL), boundary_mesh(&mesh),
    boundary_marker(bdr_attr_marker), boundary_fespace(&fespace),
    boundary_mat_op(&mat_op), boundary_lod(lod), boundary_scaling(scaling)
{
  MFEM_VERIFY(kind == Kind::FLUX_Q || kind == Kind::CURRENT_J ||
                  kind == Kind::ENERGY_E || kind == Kind::ENERGY_M,
              "Unsupported boundary material point field kind!");
  EnsureBoundaryEvaluator();
  CacheBoundaryMetadata();
  boundary_eval.reset();
}

PointFieldEvaluator::PointFieldEvaluator(Kind kind, const Mesh &mesh,
                                         const mfem::Array<int> &bdr_attr_marker,
                                         const mfem::ParFiniteElementSpace &nd_fespace,
                                         const mfem::ParFiniteElementSpace &rt_fespace,
                                         const MaterialOperator &mat_op, int lod,
                                         double scaling)
  : location(MeshEntityType::Boundary), kind(kind),
    boundary_type(BoundaryEvaluatorType::POYNTING), boundary_mesh(&mesh),
    boundary_marker(bdr_attr_marker), boundary_nd_fespace(&nd_fespace),
    boundary_rt_fespace(&rt_fespace), boundary_mat_op(&mat_op), boundary_lod(lod),
    boundary_scaling(scaling)
{
  MFEM_VERIFY(kind == Kind::POYNTING, "Unsupported two-field boundary point kind!");
  EnsureBoundaryEvaluator();
  CacheBoundaryMetadata();
  boundary_eval.reset();
}

PointFieldEvaluator::~PointFieldEvaluator() = default;

void PointFieldEvaluator::EnsureBoundaryEvaluator() const
{
  if (location != MeshEntityType::Boundary || boundary_eval)
  {
    return;
  }
  MFEM_VERIFY(boundary_mesh, "Missing mesh for boundary point field evaluator!");
  switch (boundary_type)
  {
    case BoundaryEvaluatorType::FIELD:
      boundary_eval.reset(new SurfaceFunctional(kind, *boundary_mesh, boundary_marker,
                                                *boundary_fespace, boundary_lod));
      break;
    case BoundaryEvaluatorType::MATERIAL:
      boundary_eval.reset(new SurfaceFunctional(kind, *boundary_mesh, boundary_marker,
                                                *boundary_fespace, *boundary_mat_op,
                                                boundary_lod, boundary_scaling));
      break;
    case BoundaryEvaluatorType::POYNTING:
      boundary_eval.reset(new SurfaceFunctional(kind, *boundary_mesh, boundary_marker,
                                                *boundary_nd_fespace,
                                                *boundary_rt_fespace, *boundary_mat_op,
                                                boundary_lod, boundary_scaling));
      break;
    case BoundaryEvaluatorType::NONE:
      MFEM_ABORT("Missing boundary point field evaluator type!");
  }
}

void PointFieldEvaluator::CacheBoundaryMetadata()
{
  MFEM_VERIFY(location == MeshEntityType::Boundary && boundary_eval,
              "Boundary point field metadata requires an assembled evaluator!");
  valid = boundary_eval->IsValid();
  if (valid)
  {
    boundary_buffer_size = boundary_eval->BufferSize();
    boundary_buffer_bases = boundary_eval->BufferBases();
  }
}

bool PointFieldEvaluator::ReleaseBoundaryEvaluatorAfterUse() const
{
  if (retain_boundary_eval_count > 0)
  {
    --retain_boundary_eval_count;
    return false;
  }
  return kind != Kind::FIELD_E;
}

void PointFieldEvaluator::MaybeReleaseBoundaryEvaluator() const
{
  if (ReleaseBoundaryEvaluatorAfterUse())
  {
    boundary_eval.reset();
  }
}

bool PointFieldEvaluator::IsValid() const
{
  return valid;
}

int PointFieldEvaluator::BufferSize() const
{
  return location == MeshEntityType::Domain ? domain_eval->BufferSize()
                                            : boundary_buffer_size;
}

const std::vector<int> &PointFieldEvaluator::BufferBases() const
{
  return location == MeshEntityType::Domain ? domain_eval->BufferBases()
                                            : boundary_buffer_bases;
}

void PointFieldEvaluator::Eval(const GridFunction *E, const GridFunction *B,
                               Vector &out) const
{
  MFEM_VERIFY(location == MeshEntityType::Domain && domain_eval,
              "GridFunction output is only supported for domain point fields!");
  domain_eval->Eval(E, B, out);
}

void PointFieldEvaluator::EvalBuffer(const Vector &u, Vector &buffer) const
{
  MFEM_VERIFY(IsValid() && location == MeshEntityType::Boundary && kind != Kind::POYNTING,
              "Vector EvalBuffer is only valid for single-field boundary point fields!");
  EnsureBoundaryEvaluator();
  boundary_eval->EvalBuffer(u, buffer);
  MaybeReleaseBoundaryEvaluator();
}

void PointFieldEvaluator::EvalBuffer(const GridFunction &u, Vector &buffer) const
{
  MFEM_VERIFY(IsValid() && location == MeshEntityType::Boundary && kind != Kind::POYNTING,
              "GridFunction EvalBuffer is only valid for single-field boundary point "
              "fields!");
  EnsureBoundaryEvaluator();
  boundary_eval->EvalBuffer(u, buffer);
  MaybeReleaseBoundaryEvaluator();
}

void PointFieldEvaluator::EvalBuffer(const GridFunction *E, const GridFunction *B,
                                     Vector &buffer) const
{
  MFEM_VERIFY(IsValid(), "EvalBuffer called on an invalid point field evaluator!");
  if (location == MeshEntityType::Domain)
  {
    domain_eval->EvalBuffer(E, B, buffer);
    return;
  }

  EnsureBoundaryEvaluator();
  MFEM_VERIFY(boundary_eval, "Missing boundary point field evaluator!");
  if (kind == Kind::POYNTING)
  {
    MFEM_VERIFY(E && B, "Boundary Poynting point field requires E and B fields!");
    boundary_eval->EvalBuffer(*E, *B, buffer);
  }
  else
  {
    const GridFunction *u = E ? E : B;
    MFEM_VERIFY(u, "Boundary point field requires a source field!");
    boundary_eval->EvalBuffer(*u, buffer);
  }
  MaybeReleaseBoundaryEvaluator();
}

}  // namespace palace
