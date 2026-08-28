// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/boundary_derived_field_bundle.hpp"

#include <cstdlib>
#include <cstring>
#include <deque>
#include <map>
#include <utility>
#include "fem/boundary_physical_trace.hpp"
#include "fem/ceed_group_operator.hpp"
#include "fem/coefficient.hpp"
#include "fem/face_sampling_plan.hpp"
#include "fem/gridfunction.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/32/surf_32_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

using RouteKey = std::vector<long long>;

long long EncodeDouble(double x)
{
  static_assert(sizeof(long long) == sizeof(double));
  long long bits;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

void AppendRule(RouteKey &key, const std::vector<mfem::IntegrationPoint> &points)
{
  key.push_back(static_cast<long long>(points.size()));
  for (const auto &ip : points)
  {
    key.push_back(EncodeDouble(ip.x));
    key.push_back(EncodeDouble(ip.y));
    key.push_back(EncodeDouble(ip.z));
    key.push_back(EncodeDouble(ip.weight));
  }
}

// The temporary references are released after libCEED takes its own operator references.
struct CeedAssemblyScratch
{
  Ceed ceed;
  std::vector<CeedVector> vecs;
  std::vector<CeedElemRestriction> restrictions;

  explicit CeedAssemblyScratch(Ceed ceed_) : ceed(ceed_) {}
  ~CeedAssemblyScratch()
  {
    for (auto &vec : vecs)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&vec));
    }
    for (auto &restriction : restrictions)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&restriction));
    }
  }
};

int SideAttribute(const mfem::ParMesh &pmesh, const FaceSamplingPlan::Entry &entry,
                  int side)
{
  if (side == 0)
  {
    MFEM_VERIFY(entry.elem_a >= 0, "Boundary derived bundle is missing side-A element!");
    return pmesh.GetAttribute(entry.elem_a);
  }
  MFEM_VERIFY(entry.elem_b >= 0, "Boundary derived bundle is missing side-B element!");
  return entry.ghost_b ? entry.ghost_attr : pmesh.GetAttribute(entry.elem_b);
}

}  // namespace

struct BoundaryDerivedFieldBundle::Impl
{
  // groups retains the one-phase operator for real-valued runs and callers supplying an
  // arbitrary source Vector. complex_groups is only used for the registered complex E/B
  // pair and binds all real/imaginary physical traces in one derived-field launch.
  std::vector<fem::CeedGroupOperator> groups, complex_groups;
  // Device-capable fixed geometry and attribute inputs must outlive the operators.
  std::deque<Vector> fixed_inputs;
};

std::atomic<long long> BoundaryDerivedFieldBundle::bundle_count{0};
std::atomic<long long> BoundaryDerivedFieldBundle::field_count{0};
std::atomic<long long> BoundaryDerivedFieldBundle::qfunction_applies{0};
std::atomic<long long> BoundaryDerivedFieldBundle::phase_count{0};
std::atomic<long long> BoundaryDerivedFieldBundle::avoided_phase_operator_applies{0};
std::atomic<long long> BoundaryDerivedFieldBundle::avoided_independent_applies{0};
std::atomic<long long> BoundaryDerivedFieldBundle::packed_bytes{0};
std::atomic<long long> BoundaryDerivedFieldBundle::cache_hits{0};

BoundaryDerivedFieldBundle::BoundaryDerivedFieldBundle(
    const Mesh &mesh_, const mfem::Array<int> &bdr_attr_marker,
    const MaterialOperator &mat_op_, const mfem::ParFiniteElementSpace &nd_fespace_,
    const mfem::ParFiniteElementSpace &rt_fespace_,
    std::shared_ptr<const FaceSamplingPlan> sampling_plan_,
    std::shared_ptr<BoundaryPhysicalTraceCache> trace_cache_, const GridFunction &E_,
    const GridFunction &B_, double electric_scaling, double magnetic_scaling)
  : mesh(&mesh_), mat_op(&mat_op_), nd_fespace(&nd_fespace_), rt_fespace(&rt_fespace_),
    E(&E_), B(&B_), sampling_plan(std::move(sampling_plan_)),
    trace_cache(std::move(trace_cache_)), impl(std::make_unique<Impl>())
{
  // The 2D integral-mapped scalar B path has a deliberately different current/energy
  // contract. Keep its established trace-backed consumers rather than approximating it
  // with this 3D vector bundle.
  if (mesh->Dimension() != 3 || mesh->SpaceDimension() != 3 || !sampling_plan ||
      !trace_cache || !BoundaryPhysicalTraceCache::Supports(*nd_fespace) ||
      !BoundaryPhysicalTraceCache::Supports(*rt_fespace))
  {
    return;
  }
  MFEM_VERIFY(sampling_plan->Matches(mesh_, bdr_attr_marker, sampling_plan->LOD()),
              "Boundary derived bundle received an incompatible sampling plan!");
  MFEM_VERIFY(E->ParFESpace() == nd_fespace && B->ParFESpace() == rt_fespace,
              "Boundary derived bundle sources do not match their trace spaces!");

  const mfem::ParMesh &pmesh = mesh->Get();
  Ceed ceed = ceed::internal::GetCeedObjects()[0];

  struct Route
  {
    int side = 0;
    double normal_sign = 1.0;
    double average_scale = 1.0;
    std::vector<const FaceSamplingPlan::Entry *> entries;
  };
  std::map<RouteKey, Route> routes;
  for (const auto &entry : sampling_plan->Entries())
  {
    for (int side = 0; side < 2; side++)
    {
      if (side == 1 && entry.elem_b < 0)
      {
        continue;
      }
      const auto &points = side ? entry.points_b : entry.points_a;
      const double normal_sign =
          side ? (entry.flip ? 1.0 : -1.0) : (entry.flip ? -1.0 : 1.0);
      // Only interfaces average physical fields. Exterior faces retain their single
      // physical side, while charge/current are oriented in both cases.
      const double average_scale = entry.elem_b >= 0 ? 0.5 : 1.0;
      RouteKey key;
      // Keep a separate fused apply for every canonical plan group, physical side, and
      // ghost route. The fixed EVAL_NONE restrictions below carry only those routes.
      key.push_back(side);
      key.push_back(side && entry.ghost_b);
      key.push_back(static_cast<long long>(side ? entry.vol_geom_b : entry.vol_geom_a));
      key.push_back(normal_sign < 0.0 ? -1 : 1);
      key.push_back(average_scale == 0.5 ? 1 : 0);
      AppendRule(key, points);
      auto [it, inserted] = routes.emplace(std::move(key), Route{});
      if (inserted)
      {
        it->second.side = side;
        it->second.normal_sign = normal_sign;
        it->second.average_scale = average_scale;
      }
      it->second.entries.push_back(&entry);
    }
  }

  MaterialPropertyCoefficient epsilon(mat_op->GetAttributeToMaterial(),
                                      mat_op->GetPermittivityReal());
  MaterialPropertyCoefficient invmu(mat_op->GetAttributeToMaterial(),
                                    mat_op->GetCurlCurlInvPermeability());
  const auto epsilon_ctx = ceed::PopulateCoefficientContext(3, &epsilon);
  const auto invmu_ctx = ceed::PopulateCoefficientContext(3, &invmu);

  for (const auto &item : routes)
  {
    const auto &route = item.second;  // Remove in C++20.
    const std::size_t num_elem = route.entries.size();
    MFEM_VERIFY(num_elem > 0, "Empty boundary derived route!");
    const auto &first = *route.entries.front();
    const auto &points = route.side ? first.points_b : first.points_a;
    const auto &permutation =
        route.side ? first.canonical_to_original_b : first.canonical_to_original_a;
    const int nq = static_cast<int>(points.size());
    MFEM_VERIFY(nq > 0 && permutation.size() == static_cast<std::size_t>(nq),
                "Invalid canonical boundary route in derived bundle!");

    auto &face_jac = impl->fixed_inputs.emplace_back(static_cast<int>(num_elem) * 6 * nq);
    auto &attrs = impl->fixed_inputs.emplace_back(static_cast<int>(num_elem) * nq);
    face_jac = 0.0;
    attrs = 0.0;
    std::vector<CeedInt> geom_offsets(num_elem * nq), attr_offsets(num_elem * nq),
        trace_offsets(num_elem * nq), out_offsets(num_elem * nq);
    for (std::size_t e = 0; e < num_elem; e++)
    {
      const auto &entry = *route.entries[e];
      const auto &entry_perm =
          route.side ? entry.canonical_to_original_b : entry.canonical_to_original_a;
      MFEM_VERIFY(entry_perm.size() == static_cast<std::size_t>(nq) &&
                      entry.face_jacobians.size() == static_cast<std::size_t>(6 * nq),
                  "Invalid face Jacobian route in boundary derived bundle!");
      const int point_base = sampling_plan->BufferBases()[entry.bdr_elem];
      const int attr = SideAttribute(pmesh, entry, route.side);
      // The libCEED material contexts are indexed by the same process-local attribute
      // numbering supplied by Mesh::GetCeedAttributes() in the legacy physical-trace
      // path. This is essential on ranks where global material attributes are sparse.
      const int ceed_attr = mesh->GetCeedAttributes().at(attr);
      for (int q = 0; q < nq; q++)
      {
        const int original_q = entry_perm[static_cast<std::size_t>(q)];
        MFEM_VERIFY(original_q >= 0 && original_q < nq,
                    "Invalid canonical output route in boundary derived bundle!");
        geom_offsets[e * nq + q] = static_cast<CeedInt>(e * 6 * nq + q);
        attr_offsets[e * nq + q] = static_cast<CeedInt>(e * nq + q);
        trace_offsets[e * nq + q] = static_cast<CeedInt>(point_base + q);
        out_offsets[e * nq + q] = static_cast<CeedInt>(point_base + original_q);
        attrs[e * nq + q] = ceed_attr;
        for (int c = 0; c < 6; c++)
        {
          face_jac[e * 6 * nq + c * nq + q] = entry.face_jacobians[c * nq + original_q];
        }
      }
    }
    face_jac.UseDevice(true);
    attrs.UseDevice(true);
    auto AssembleRoute =
        [&](bool batched_complex, std::vector<fem::CeedGroupOperator> &groups)
    {
      CeedAssemblyScratch scratch(ceed);
      std::vector<ceed::CeedFunctionalFieldInput> inputs;
      std::vector<std::pair<std::string, int>> field_sources;
      auto AddFixed = [&](const std::string &name, Vector &data, int num_comp,
                          CeedInt comp_stride, CeedSize lsize,
                          const std::vector<CeedInt> &offsets)
      {
        CeedElemRestriction restriction;
        PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem),
                                                       nq, num_comp, comp_stride, lsize,
                                                       CEED_MEM_HOST, CEED_COPY_VALUES,
                                                       offsets.data(), &restriction));
        CeedVector vector;
        ceed::InitCeedVector(data, ceed, &vector);
        inputs.push_back({name, vector, restriction, nullptr, ceed::EvalMode::None});
        scratch.vecs.push_back(vector);
        scratch.restrictions.push_back(restriction);
      };
      AddFixed("grad_x_f", face_jac, 6, nq, face_jac.Size(), geom_offsets);
      AddFixed("attr", attrs, 1, nq, attrs.Size(), attr_offsets);

      auto AddTrace = [&](const std::string &name, int source)
      {
        CeedElemRestriction restriction;
        const int trace_points = trace_cache->NumPoints();
        PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                                 ceed, static_cast<CeedInt>(num_elem), nq, 3, trace_points,
                                 static_cast<CeedSize>(3) * trace_points, CEED_MEM_HOST,
                                 CEED_COPY_VALUES, trace_offsets.data(), &restriction));
        CeedVector vector;
        PalaceCeedCall(
            ceed, CeedVectorCreate(ceed, static_cast<CeedSize>(3) * trace_points, &vector));
        inputs.push_back({name, vector, restriction, nullptr, ceed::EvalMode::None});
        field_sources.emplace_back(name, source);
        scratch.vecs.push_back(vector);
        scratch.restrictions.push_back(restriction);
      };
      if (batched_complex)
      {
        AddTrace("e_real", route.side);
        AddTrace("b_real", 2 + route.side);
        AddTrace("e_imag", 4 + route.side);
        AddTrace("b_imag", 6 + route.side);
      }
      else
      {
        AddTrace("e", route.side);
        AddTrace("b", 2 + route.side);
      }

      const int output_comp = batched_complex ? 25 : 15;
      CeedElemRestriction out_restriction;
      const int trace_points = trace_cache->NumPoints();
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreate(
                    ceed, static_cast<CeedInt>(num_elem), nq, output_comp, trace_points,
                    static_cast<CeedSize>(output_comp) * trace_points, CEED_MEM_HOST,
                    CEED_COPY_VALUES, out_offsets.data(), &out_restriction));
      scratch.restrictions.push_back(out_restriction);

      std::vector<CeedIntScalar> ctx(9);
      ctx[0].second = route.normal_sign;
      ctx[1].second = electric_scaling;  // Surface charge.
      ctx[2].second = magnetic_scaling;  // Surface current.
      ctx[3].second = electric_scaling;  // Electric energy.
      ctx[4].second = magnetic_scaling;  // Magnetic energy.
      ctx[5].second = magnetic_scaling;  // Poynting vector.
      ctx[6].second = route.average_scale;
      ctx[7].first = 9;
      ctx.insert(ctx.end(), epsilon_ctx.begin(), epsilon_ctx.end());
      ctx[8].first = static_cast<CeedInt>(ctx.size());
      ctx.insert(ctx.end(), invmu_ctx.begin(), invmu_ctx.end());

      ceed::CeedQFunctionInfo info;
      info.apply_qf = batched_complex ? f_eval_bdr_derived_trace_complex_32
                                      : f_eval_bdr_derived_trace_32;
      info.apply_qf_path =
          batched_complex
              ? PalaceQFunctionRelativePath(f_eval_bdr_derived_trace_complex_32_loc)
              : PalaceQFunctionRelativePath(f_eval_bdr_derived_trace_32_loc);
      CeedOperator op;
      ceed::AssembleCeedPointEvaluator(info, ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                                       ceed, inputs, output_comp, out_restriction, &op);
      groups.push_back({ceed, op, std::move(field_sources)});
      fem::CacheGroupOperatorFieldVectors(groups.back());
      // Physical trace vectors are intentionally unbound at setup. At apply time the
      // group takes and repoints the handles to the source identities for this save.
      groups.back().field_vectors_detached = true;
    };

    // Keep the original one-phase route for real-valued solves and public callers that
    // provide an arbitrary source Vector. Registered complex E/B phases use only the
    // four-trace operator below, never two applications of this fallback route.
    AssembleRoute(false, impl->groups);
    if (E->HasImag() && B->HasImag())
    {
      AssembleRoute(true, impl->complex_groups);
    }
  }

  valid = true;
  ++bundle_count;
  field_count += 7;
}

BoundaryDerivedFieldBundle::~BoundaryDerivedFieldBundle()
{
  if (impl)
  {
    fem::DestroyGroupOperators(impl->groups);
    fem::DestroyGroupOperators(impl->complex_groups);
  }
}

bool BoundaryDerivedFieldBundle::Supports(
    PointFieldKind kind, const mfem::ParFiniteElementSpace *fespace_a,
    const mfem::ParFiniteElementSpace *fespace_b) const
{
  if (!valid)
  {
    return false;
  }
  switch (kind)
  {
    case PointFieldKind::FIELD_E:
    case PointFieldKind::FLUX_Q:
    case PointFieldKind::ENERGY_E:
      return fespace_a == nd_fespace && !fespace_b;
    case PointFieldKind::FIELD_B:
    case PointFieldKind::CURRENT_J:
    case PointFieldKind::ENERGY_M:
      return fespace_a == rt_fespace && !fespace_b;
    case PointFieldKind::POYNTING:
      return fespace_a == nd_fespace && fespace_b == rt_fespace;
    default:
      return false;
  }
}

const FaceSamplingPlan *BoundaryDerivedFieldBundle::SamplingPlan() const
{
  return sampling_plan.get();
}

int BoundaryDerivedFieldBundle::BufferSize(PointFieldKind kind) const
{
  return BufferNumComp(kind) * trace_cache->NumPoints();
}

int BoundaryDerivedFieldBundle::BufferNumComp(PointFieldKind kind) const
{
  return SliceNumComp(ToSlice(kind));
}

const std::vector<int> &BoundaryDerivedFieldBundle::BufferBases() const
{
  return sampling_plan->BufferBases();
}

const Vector &BoundaryDerivedFieldBundle::GetRegisteredPacked() const
{
  MFEM_VERIFY(valid, "Invalid boundary-derived bundle packed output request!");
  return Get(E->Real(), B->Real()).packed;
}

int BoundaryDerivedFieldBundle::RegisteredPackedSize() const
{
  MFEM_VERIFY(valid, "Invalid boundary-derived bundle packed output size request!");
  return (E->HasImag() && B->HasImag() ? 25 : 15) * trace_cache->NumPoints();
}

BoundaryDerivedFieldBundle::PackedSlice
BoundaryDerivedFieldBundle::RegisteredSlice(PointFieldKind kind, bool imaginary_phase) const
{
  MFEM_VERIFY(valid, "Invalid boundary-derived bundle packed slice request!");
  const bool batched_complex = E->HasImag() && B->HasImag();
  MFEM_VERIFY(!imaginary_phase || batched_complex,
              "Imaginary boundary-derived packed slice requires registered complex E/B!");
  const Slice slice = ToSlice(kind);
  MFEM_VERIFY(!imaginary_phase || (slice != Slice::ENERGY_E && slice != Slice::ENERGY_M &&
                                   slice != Slice::POYNTING),
              "Quadratic boundary-derived fields have one combined complex slice, not "
              "an imaginary-phase slice!");
  const int points = trace_cache->NumPoints();
  return {SliceOffset(slice, points, batched_complex, imaginary_phase),
          SliceNumComp(slice) * points, SliceNumComp(slice)};
}

void BoundaryDerivedFieldBundle::Invalidate()
{
  ++generation;
  phases.clear();
}

BoundaryDerivedFieldBundle::Slice BoundaryDerivedFieldBundle::ToSlice(PointFieldKind kind)
{
  switch (kind)
  {
    case PointFieldKind::FIELD_E:
      return Slice::FIELD_E;
    case PointFieldKind::FIELD_B:
      return Slice::FIELD_B;
    case PointFieldKind::FLUX_Q:
      return Slice::FLUX_Q;
    case PointFieldKind::CURRENT_J:
      return Slice::CURRENT_J;
    case PointFieldKind::ENERGY_E:
      return Slice::ENERGY_E;
    case PointFieldKind::ENERGY_M:
      return Slice::ENERGY_M;
    case PointFieldKind::POYNTING:
      return Slice::POYNTING;
    default:
      MFEM_ABORT("Unsupported boundary-derived bundle output kind!");
  }
}

int BoundaryDerivedFieldBundle::SliceOffset(Slice slice, int num_points,
                                            bool batched_complex, bool imaginary_phase)
{
  if (!batched_complex)
  {
    switch (slice)
    {
      case Slice::FIELD_E:
        return 0;
      case Slice::FIELD_B:
        return 3 * num_points;
      case Slice::FLUX_Q:
        return 6 * num_points;
      case Slice::CURRENT_J:
        return 7 * num_points;
      case Slice::ENERGY_E:
        return 10 * num_points;
      case Slice::ENERGY_M:
        return 11 * num_points;
      case Slice::POYNTING:
        return 12 * num_points;
    }
  }
  else
  {
    // Complex layout is component-major and phase-major for the linear slices:
    // [E_r(3), B_r(3), Q_r, J_r(3), E_i(3), B_i(3), Q_i, J_i(3), Ue, Um, S(3)].
    switch (slice)
    {
      case Slice::FIELD_E:
        return (imaginary_phase ? 10 : 0) * num_points;
      case Slice::FIELD_B:
        return (imaginary_phase ? 13 : 3) * num_points;
      case Slice::FLUX_Q:
        return (imaginary_phase ? 16 : 6) * num_points;
      case Slice::CURRENT_J:
        return (imaginary_phase ? 17 : 7) * num_points;
      case Slice::ENERGY_E:
        return 20 * num_points;
      case Slice::ENERGY_M:
        return 21 * num_points;
      case Slice::POYNTING:
        return 22 * num_points;
    }
  }
  return 0;
}

int BoundaryDerivedFieldBundle::SliceNumComp(Slice slice)
{
  switch (slice)
  {
    case Slice::FIELD_E:
    case Slice::FIELD_B:
    case Slice::CURRENT_J:
    case Slice::POYNTING:
      return 3;
    case Slice::FLUX_Q:
    case Slice::ENERGY_E:
    case Slice::ENERGY_M:
      return 1;
  }
  return 0;
}

bool BoundaryDerivedFieldBundle::BatchesComplexPhases(const Vector &E_source,
                                                      const Vector &B_source) const
{
  if (!E->HasImag() || !B->HasImag())
  {
    return false;
  }
  return (&E_source == &E->Real() && &B_source == &B->Real()) ||
         (&E_source == &E->Imag() && &B_source == &B->Imag());
}

bool BoundaryDerivedFieldBundle::BatchesComplexPhases(const Vector &source) const
{
  return E->HasImag() && B->HasImag() &&
         (&source == &E->Real() || &source == &B->Real() || &source == &E->Imag() ||
          &source == &B->Imag());
}

bool BoundaryDerivedFieldBundle::IsImaginaryPhase(const Vector &E_source,
                                                  const Vector &B_source) const
{
  return BatchesComplexPhases(E_source, B_source) && &E_source == &E->Imag();
}

const BoundaryDerivedFieldBundle::CachedPhase &
BoundaryDerivedFieldBundle::Get(const Vector &E_source, const Vector &B_source) const
{
  trace_cache->VerifyCurrent();
  const bool batched_complex = BatchesComplexPhases(E_source, B_source);
  // The registered real and imaginary source pairs are one cache identity because one
  // application consumes both. Arbitrary source vectors retain the previous pair key.
  const auto key = batched_complex ? std::make_pair(&E->Real(), &B->Real())
                                   : std::make_pair(&E_source, &B_source);
  auto [it, inserted] = phases.emplace(key, CachedPhase{});
  auto &phase = it->second;
  if (!inserted && phase.generation == generation)
  {
    ++cache_hits;
    return phase;
  }

  const auto e_traces =
      trace_cache->Get(*nd_fespace, batched_complex ? E->Real() : E_source);
  const auto b_traces =
      trace_cache->Get(*rt_fespace, batched_complex ? B->Real() : B_source);
  MFEM_VERIFY(e_traces.num_comp == 3 && b_traces.num_comp == 3,
              "3D derived boundary bundle requires vector physical traces!");
  const int points = trace_cache->NumPoints();
  phase.batched_complex = batched_complex;
  phase.packed.SetSize((batched_complex ? 25 : 15) * points);
  phase.packed.UseDevice(true);
  phase.packed = 0.0;
  if (batched_complex)
  {
    const auto e_imag_traces = trace_cache->Get(*nd_fespace, E->Imag());
    const auto b_imag_traces = trace_cache->Get(*rt_fespace, B->Imag());
    MFEM_VERIFY(e_imag_traces.num_comp == 3 && b_imag_traces.num_comp == 3,
                "3D derived boundary bundle requires vector physical traces!");
    fem::ApplyAddComplexGroupOperators(impl->complex_groups,
                                       {e_traces.side_a, e_traces.side_b, b_traces.side_a,
                                        b_traces.side_b, e_imag_traces.side_a,
                                        e_imag_traces.side_b, b_imag_traces.side_a,
                                        b_imag_traces.side_b},
                                       phase.packed);
    const long long applies = static_cast<long long>(impl->complex_groups.size());
    qfunction_applies += applies;
    phase_count += 2 * applies;
    avoided_phase_operator_applies += applies;
    // Eleven public outputs (eight phase-resolved linear and three combined quadratic)
    // share one application for this route.
    avoided_independent_applies += 10 * applies;
  }
  else
  {
    fem::ApplyAddGroupOperators(
        impl->groups, {e_traces.side_a, e_traces.side_b, b_traces.side_a, b_traces.side_b},
        phase.packed);
    const long long applies = static_cast<long long>(impl->groups.size());
    qfunction_applies += applies;
    phase_count += applies;
    avoided_independent_applies += 6 * applies;
  }
  phase.generation = generation;
  packed_bytes += static_cast<long long>(phase.packed.Size()) * sizeof(double);
  return phase;
}

std::pair<const Vector *, const Vector *>
BoundaryDerivedFieldBundle::PartnerSources(const Vector &source) const
{
  if (&source == &E->Real() || &source == &B->Real())
  {
    return {&E->Real(), &B->Real()};
  }
  if (E->HasImag() && B->HasImag() && (&source == &E->Imag() || &source == &B->Imag()))
  {
    return {&E->Imag(), &B->Imag()};
  }
  MFEM_ABORT("Boundary-derived output was invoked with an unregistered source vector!");
}

void BoundaryDerivedFieldBundle::Copy(PointFieldKind kind, const Vector &source,
                                      Vector &buffer, bool zero) const
{
  if (&source == &E->Real() || &source == &B->Real() ||
      (E->HasImag() && B->HasImag() && (&source == &E->Imag() || &source == &B->Imag())))
  {
    const auto [e_source, b_source] = PartnerSources(source);
    Copy(kind, *e_source, *b_source, buffer, zero);
    return;
  }

  // Public PointFieldEvaluator callers may evaluate another Vector on the same FE space.
  // Pair the unused family with a registered source; only the requested E- or B-family
  // slice is exposed, so this preserves the former independent-evaluator contract.
  switch (kind)
  {
    case PointFieldKind::FIELD_E:
    case PointFieldKind::FLUX_Q:
    case PointFieldKind::ENERGY_E:
      Copy(kind, source, B->Real(), buffer, zero);
      return;
    case PointFieldKind::FIELD_B:
    case PointFieldKind::CURRENT_J:
    case PointFieldKind::ENERGY_M:
      Copy(kind, E->Real(), source, buffer, zero);
      return;
    default:
      MFEM_ABORT("Unsupported single-source boundary-derived bundle output!");
  }
}

void BoundaryDerivedFieldBundle::Copy(PointFieldKind kind, const Vector &E_source,
                                      const Vector &B_source, Vector &buffer,
                                      bool zero) const
{
  const auto &phase = Get(E_source, B_source);
  const int points = trace_cache->NumPoints();
  const Slice slice_id = ToSlice(kind);
  const int size = SliceNumComp(slice_id) * points;
  MFEM_VERIFY(buffer.Size() == size,
              "Invalid output buffer size for boundary-derived field bundle!");
  Vector slice;
  // The packed source and the writer buffer are both MFEM device-aware vectors. Vector
  // assignment/AXPY preserve their device path while exposing only this public slice.
  slice.MakeRef(const_cast<Vector &>(phase.packed),
                SliceOffset(slice_id, points, phase.batched_complex,
                            IsImaginaryPhase(E_source, B_source)),
                size);
  if (zero)
  {
    buffer = slice;
  }
  else
  {
    buffer += slice;
  }
}

long long BoundaryDerivedFieldBundle::BundleCount()
{
  return bundle_count.load();
}

long long BoundaryDerivedFieldBundle::FieldCount()
{
  return field_count.load();
}

long long BoundaryDerivedFieldBundle::QFunctionApplyCount()
{
  return qfunction_applies.load();
}

long long BoundaryDerivedFieldBundle::PhaseCount()
{
  return phase_count.load();
}

long long BoundaryDerivedFieldBundle::AvoidedPhaseOperatorApplyCount()
{
  return avoided_phase_operator_applies.load();
}

long long BoundaryDerivedFieldBundle::AvoidedIndependentApplyCount()
{
  return avoided_independent_applies.load();
}

long long BoundaryDerivedFieldBundle::PackedByteCount()
{
  return packed_bytes.load();
}

long long BoundaryDerivedFieldBundle::CacheHitCount()
{
  return cache_hits.load();
}

void BoundaryDerivedFieldBundle::PrintProfile() const
{
  if (!std::getenv("PALACE_BOUNDARY_DERIVED_PROFILE"))
  {
    return;
  }
  Mpi::Print(mesh->GetComm(),
             "BoundaryDerivedFieldBundle profile bundles={} fields={} qfunction_applies={} "
             "phases={} avoided_phase_operator_applies={} avoided_independent_applies={} "
             "packed_bytes={} cache_hits={} generation={}\n",
             BundleCount(), FieldCount(), QFunctionApplyCount(), PhaseCount(),
             AvoidedPhaseOperatorApplyCount(), AvoidedIndependentApplyCount(),
             PackedByteCount(), CacheHitCount(), generation);
}

}  // namespace palace
