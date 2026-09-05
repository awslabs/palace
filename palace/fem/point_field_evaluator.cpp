// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/point_field_evaluator.hpp"

#include "fem/boundary_derived_field_bundle.hpp"
#include "fem/boundary_physical_trace.hpp"
#include "fem/domain_point_field_evaluator.hpp"
#include "fem/face_sampling_plan.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"

namespace palace
{

namespace
{

bool SupportsBoundarySamplingPlan(const Mesh &mesh)
{
  return (mesh.Dimension() == 2 && mesh.SpaceDimension() == 2) ||
         (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3);
}

DomainPointFieldEvaluator::Kind ToDomainKind(PointFieldEvaluator::Kind kind)
{
  switch (kind)
  {
    case PointFieldEvaluator::Kind::FIELD_E:
      return DomainPointFieldEvaluator::Kind::FIELD_E;
    case PointFieldEvaluator::Kind::FIELD_B:
      return DomainPointFieldEvaluator::Kind::FIELD_B;
    case PointFieldEvaluator::Kind::FIELD_H1:
      return DomainPointFieldEvaluator::Kind::FIELD_H1;
    case PointFieldEvaluator::Kind::ENERGY_E:
      return DomainPointFieldEvaluator::Kind::ENERGY_E;
    case PointFieldEvaluator::Kind::ENERGY_M:
      return DomainPointFieldEvaluator::Kind::ENERGY_M;
    case PointFieldEvaluator::Kind::POYNTING:
      return DomainPointFieldEvaluator::Kind::POYNTING;
    case PointFieldEvaluator::Kind::MODE_SN:
      return DomainPointFieldEvaluator::Kind::MODE_SN;
    default:
      MFEM_ABORT("Unsupported domain point field kind!");
  }
}

}  // namespace

PointFieldEvaluator::PointFieldEvaluator(Kind kind, const Mesh &mesh,
                                         const MaterialOperator &mat_op,
                                         const mfem::ParFiniteElementSpace *nd_fespace,
                                         const mfem::ParFiniteElementSpace *rt_fespace,
                                         const mfem::ParFiniteElementSpace &target_fespace,
                                         double scaling, bool build_gridfunction,
                                         bool build_buffer)
  : location(MeshEntityType::Domain), kind(kind)
{
  MFEM_VERIFY(kind == Kind::FIELD_E || kind == Kind::FIELD_B || kind == Kind::FIELD_H1 ||
                  kind == Kind::ENERGY_E || kind == Kind::ENERGY_M ||
                  kind == Kind::POYNTING || kind == Kind::MODE_SN,
              "Unsupported domain point field kind!");
  domain_eval = std::make_unique<DomainPointFieldEvaluator>(
      ToDomainKind(kind), mesh, mat_op, nd_fespace, rt_fespace, target_fespace, scaling,
      build_gridfunction, build_buffer);
  valid = true;
}

PointFieldEvaluator::PointFieldEvaluator(
    Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
    const mfem::ParFiniteElementSpace &fespace, int lod,
    std::shared_ptr<const FaceSamplingPlan> sampling_plan,
    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache_,
    std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle_)
  : location(MeshEntityType::Boundary), kind(kind), trace_cache(std::move(trace_cache_)),
    derived_bundle(std::move(derived_bundle_))
{
  MFEM_VERIFY(kind == Kind::FIELD_E || kind == Kind::FIELD_B || kind == Kind::FIELD_H1,
              "Unsupported boundary field point kind!");
  if (!sampling_plan && SupportsBoundarySamplingPlan(mesh))
  {
    sampling_plan = std::make_shared<FaceSamplingPlan>(mesh, bdr_attr_marker, lod);
  }
  if (trace_cache && !BoundaryPhysicalTraceCache::Supports(fespace))
  {
    trace_cache.reset();
  }
  if (derived_bundle && !derived_bundle->Supports(kind, &fespace))
  {
    derived_bundle.reset();
  }
  if (derived_bundle)
  {
    valid = true;
    return;
  }
  if (trace_cache)
  {
    trace_fespace_a = &fespace;
  }
  boundary_eval.reset(new SurfaceFunctional(kind, mesh, bdr_attr_marker, fespace, lod,
                                            std::move(sampling_plan), trace_cache));
  valid = boundary_eval->IsValid();
}

PointFieldEvaluator::PointFieldEvaluator(
    Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
    const mfem::ParFiniteElementSpace &fespace, const MaterialOperator &mat_op, int lod,
    double scaling, std::shared_ptr<const FaceSamplingPlan> sampling_plan,
    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache_,
    std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle_)
  : location(MeshEntityType::Boundary), kind(kind), trace_cache(std::move(trace_cache_)),
    derived_bundle(std::move(derived_bundle_))
{
  MFEM_VERIFY(kind == Kind::FLUX_Q || kind == Kind::CURRENT_J || kind == Kind::ENERGY_E ||
                  kind == Kind::ENERGY_M,
              "Unsupported boundary material point field kind!");
  if (!sampling_plan && SupportsBoundarySamplingPlan(mesh))
  {
    sampling_plan = std::make_shared<FaceSamplingPlan>(mesh, bdr_attr_marker, lod);
  }
  if (trace_cache && !BoundaryPhysicalTraceCache::Supports(fespace))
  {
    trace_cache.reset();
  }
  if (derived_bundle && !derived_bundle->Supports(kind, &fespace))
  {
    derived_bundle.reset();
  }
  if (derived_bundle)
  {
    valid = true;
    return;
  }
  if (trace_cache)
  {
    trace_fespace_a = &fespace;
  }
  boundary_eval.reset(new SurfaceFunctional(kind, mesh, bdr_attr_marker, fespace, mat_op,
                                            lod, scaling, std::move(sampling_plan),
                                            trace_cache));
  valid = boundary_eval->IsValid();
}

PointFieldEvaluator::PointFieldEvaluator(
    Kind kind, const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker,
    const mfem::ParFiniteElementSpace &nd_fespace,
    const mfem::ParFiniteElementSpace &rt_fespace, const MaterialOperator &mat_op, int lod,
    double scaling, std::shared_ptr<const FaceSamplingPlan> sampling_plan,
    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache_,
    std::shared_ptr<BoundaryDerivedFieldBundle> derived_bundle_)
  : location(MeshEntityType::Boundary), kind(kind), trace_cache(std::move(trace_cache_)),
    derived_bundle(std::move(derived_bundle_))
{
  MFEM_VERIFY(kind == Kind::POYNTING, "Unsupported two-field boundary point kind!");
  if (!sampling_plan && SupportsBoundarySamplingPlan(mesh))
  {
    sampling_plan = std::make_shared<FaceSamplingPlan>(mesh, bdr_attr_marker, lod);
  }
  if (trace_cache && (!BoundaryPhysicalTraceCache::Supports(nd_fespace) ||
                      !BoundaryPhysicalTraceCache::Supports(rt_fespace)))
  {
    trace_cache.reset();
  }
  if (derived_bundle && !derived_bundle->Supports(kind, &nd_fespace, &rt_fespace))
  {
    derived_bundle.reset();
  }
  if (derived_bundle)
  {
    valid = true;
    return;
  }
  if (trace_cache)
  {
    trace_fespace_a = &nd_fespace;
    trace_fespace_b = &rt_fespace;
  }
  boundary_eval.reset(new SurfaceFunctional(kind, mesh, bdr_attr_marker, nd_fespace,
                                            rt_fespace, mat_op, lod, scaling,
                                            std::move(sampling_plan), trace_cache));
  valid = boundary_eval->IsValid();
}

PointFieldEvaluator::~PointFieldEvaluator() = default;

bool PointFieldEvaluator::IsValid() const
{
  return valid;
}

int PointFieldEvaluator::BufferSize() const
{
  if (location == MeshEntityType::Domain)
  {
    return domain_eval->BufferSize();
  }
  return derived_bundle ? derived_bundle->BufferSize(kind) : boundary_eval->BufferSize();
}

int PointFieldEvaluator::BufferNumComp() const
{
  if (location == MeshEntityType::Domain)
  {
    return domain_eval->BufferNumComp();
  }
  return derived_bundle ? derived_bundle->BufferNumComp(kind)
                        : boundary_eval->BufferNumComp();
}

const std::vector<int> &PointFieldEvaluator::BufferBases() const
{
  if (location == MeshEntityType::Domain)
  {
    return domain_eval->BufferBases();
  }
  return derived_bundle ? derived_bundle->BufferBases() : boundary_eval->BufferBases();
}

const FaceSamplingPlan *PointFieldEvaluator::SamplingPlan() const
{
  if (location != MeshEntityType::Boundary)
  {
    return nullptr;
  }
  return derived_bundle ? derived_bundle->SamplingPlan()
                        : boundary_eval->sampling_plan.get();
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
  MFEM_VERIFY(IsValid() && kind != Kind::POYNTING,
              "Vector EvalBuffer is only valid for single-field point fields!");
  if (location == MeshEntityType::Domain)
  {
    domain_eval->EvalBuffer(u, buffer);
  }
  else if (derived_bundle)
  {
    derived_bundle->Copy(kind, u, buffer);
  }
  else if (trace_cache)
  {
    MFEM_VERIFY(trace_fespace_a, "Missing source space for boundary physical trace!");
    const auto traces = trace_cache->Get(*trace_fespace_a, u);
    boundary_eval->EvalBufferTrace({traces.side_a, traces.side_b, nullptr, nullptr},
                                   buffer);
  }
  else
  {
    boundary_eval->EvalBuffer(u, buffer);
  }
}

void PointFieldEvaluator::EvalBuffer(const GridFunction &u, Vector &buffer) const
{
  MFEM_VERIFY(IsValid() && location == MeshEntityType::Boundary && kind != Kind::POYNTING,
              "GridFunction EvalBuffer is only valid for single-field boundary point "
              "fields!");
  if (derived_bundle)
  {
    derived_bundle->Copy(kind, u.Real(), buffer);
    const bool quadratic = kind == Kind::ENERGY_E || kind == Kind::ENERGY_M;
    // A registered complex bundle already stores U_e/U_m as the exact real-plus-imag
    // sum. Linear callers still select their individual source phase and retain the
    // established GridFunction accumulation behavior for external source fallbacks.
    if (u.HasImag() && !(quadratic && derived_bundle->BatchesComplexPhases(u.Real())))
    {
      derived_bundle->Copy(kind, u.Imag(), buffer, false);
    }
    return;
  }
  if (!trace_cache)
  {
    boundary_eval->EvalBuffer(u, buffer);
    return;
  }
  MFEM_VERIFY(trace_fespace_a, "Missing source space for boundary physical trace!");
  const auto real_traces = trace_cache->Get(*trace_fespace_a, u.Real());
  boundary_eval->EvalBufferTrace({real_traces.side_a, real_traces.side_b, nullptr, nullptr},
                                 buffer);
  if (u.HasImag())
  {
    const auto imag_traces = trace_cache->Get(*trace_fespace_a, u.Imag());
    boundary_eval->EvalBufferTrace(
        {imag_traces.side_a, imag_traces.side_b, nullptr, nullptr}, buffer, false);
  }
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

  if (kind == Kind::POYNTING)
  {
    MFEM_VERIFY(E && B, "Boundary Poynting point field requires E and B fields!");
    if (derived_bundle)
    {
      MFEM_VERIFY(E->HasImag() == B->HasImag(),
                  "Mismatched complex E/B fields for boundary derived bundle!");
      derived_bundle->Copy(kind, E->Real(), B->Real(), buffer);
      // The registered complex derived bundle emits the same-phase Poynting sum in one
      // packed slice. Retain two calls only for an arbitrary-source fallback pair.
      if (E->HasImag() && !derived_bundle->BatchesComplexPhases(E->Real(), B->Real()))
      {
        derived_bundle->Copy(kind, E->Imag(), B->Imag(), buffer, false);
      }
      return;
    }
    MFEM_VERIFY(boundary_eval, "Missing boundary point field evaluator!");
    if (!trace_cache)
    {
      boundary_eval->EvalBuffer(*E, *B, buffer);
      return;
    }
    MFEM_VERIFY(trace_fespace_a && trace_fespace_b,
                "Missing Poynting source space for boundary physical trace!");
    MFEM_VERIFY(E->HasImag() == B->HasImag(),
                "Mismatched complex E/B fields for boundary physical trace!");
    const auto e_real = trace_cache->Get(*trace_fespace_a, E->Real());
    const auto b_real = trace_cache->Get(*trace_fespace_b, B->Real());
    boundary_eval->EvalBufferTrace(
        {e_real.side_a, e_real.side_b, b_real.side_a, b_real.side_b}, buffer);
    if (E->HasImag())
    {
      const auto e_imag = trace_cache->Get(*trace_fespace_a, E->Imag());
      const auto b_imag = trace_cache->Get(*trace_fespace_b, B->Imag());
      boundary_eval->EvalBufferTrace(
          {e_imag.side_a, e_imag.side_b, b_imag.side_a, b_imag.side_b}, buffer, false);
    }
  }
  else
  {
    MFEM_VERIFY(boundary_eval || derived_bundle, "Missing boundary point field evaluator!");
    const GridFunction *u = E ? E : B;
    MFEM_VERIFY(u, "Boundary point field requires a source field!");
    EvalBuffer(*u, buffer);
  }
}

}  // namespace palace
