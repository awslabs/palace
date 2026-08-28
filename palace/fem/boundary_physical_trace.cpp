// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/boundary_physical_trace.hpp"

#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <utility>
#include "fem/face_sampling_plan.hpp"
#include "fem/facenbrexchange.hpp"
#include "fem/fespace.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/22/eval_22_qf.h"
#include "fem/qfunctions/33/eval_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

long long EncodeDouble(double x)
{
  static_assert(sizeof(long long) == sizeof(double));
  long long bits;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

std::vector<long long> MakeUnionRuleKey(mfem::Geometry::Type geometry,
                                        const std::vector<mfem::IntegrationPoint> &points)
{
  std::vector<long long> key;
  key.reserve(2 + 3 * points.size());
  key.push_back(static_cast<long long>(geometry));
  key.push_back(static_cast<long long>(points.size()));
  for (const auto &point : points)
  {
    key.push_back(EncodeDouble(point.x));
    key.push_back(EncodeDouble(point.y));
    key.push_back(EncodeDouble(point.z));
  }
  return key;
}

const mfem::IntegrationRule *
GetRegisteredUnionIr(mfem::Geometry::Type geometry,
                     const std::vector<mfem::IntegrationPoint> &points)
{
  // MFEM finite elements cache point tabulations by IntegrationRule address. Retain each
  // coordinate-only union rule for process lifetime, as the existing mapped-rule paths
  // do. Every union point has the fixed weight 1, so leaf quadrature weights cannot
  // alter this identity.
  static std::map<std::vector<long long>, std::unique_ptr<mfem::IntegrationRule>> registry;
  static std::mutex registry_mutex;
  auto key = MakeUnionRuleKey(geometry, points);
  std::lock_guard<std::mutex> lock(registry_mutex);
  auto it = registry.find(key);
  if (it == registry.end())
  {
    auto ir = std::make_unique<mfem::IntegrationRule>(static_cast<int>(points.size()));
    for (std::size_t q = 0; q < points.size(); q++)
    {
      ir->IntPoint(static_cast<int>(q)) = points[q];
    }
    it = registry.emplace(std::move(key), std::move(ir)).first;
  }
  return it->second.get();
}

ceed::CeedQFunctionInfo GetProbeQFunction(const mfem::ParFiniteElementSpace &fespace,
                                          int dim, int sdim)
{
  ceed::CeedQFunctionInfo info;
  const auto map_type = fespace.FEColl()->GetMapType(dim);
  if (map_type == mfem::FiniteElement::VALUE)
  {
    MFEM_VERIFY(fespace.GetVDim() == 1, "Scalar physical trace source must be scalar!");
    if (sdim == 2)
    {
      info.apply_qf = f_eval_probe_l2_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_l2_22_loc);
    }
    else
    {
      info.apply_qf = f_eval_probe_l2_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_l2_33_loc);
    }
  }
  else if (map_type == mfem::FiniteElement::H_CURL)
  {
    if (sdim == 2)
    {
      info.apply_qf = f_eval_probe_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hcurl_22_loc);
    }
    else
    {
      info.apply_qf = f_eval_probe_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hcurl_33_loc);
    }
  }
  else
  {
    MFEM_VERIFY(map_type == mfem::FiniteElement::H_DIV,
                "Unsupported physical trace source map type!");
    if (sdim == 2)
    {
      info.apply_qf = f_eval_probe_hdiv_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hdiv_22_loc);
    }
    else
    {
      info.apply_qf = f_eval_probe_hdiv_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hdiv_33_loc);
    }
  }
  return info;
}

}  // namespace

std::atomic<long long> BoundaryPhysicalTraceCache::trace_groups{0};
std::atomic<long long> BoundaryPhysicalTraceCache::operator_applies{0};
std::atomic<long long> BoundaryPhysicalTraceCache::cache_hits{0};
std::atomic<long long> BoundaryPhysicalTraceCache::evaluated_points{0};
std::atomic<long long> BoundaryPhysicalTraceCache::evaluated_bytes{0};
std::atomic<long long> BoundaryPhysicalTraceCache::leaf_points{0};
std::atomic<long long> BoundaryPhysicalTraceCache::union_points{0};
std::atomic<long long> BoundaryPhysicalTraceCache::duplicate_points{0};
std::atomic<long long> BoundaryPhysicalTraceCache::local_union_groups{0};
std::atomic<long long> BoundaryPhysicalTraceCache::local_union_operators{0};
std::atomic<long long> BoundaryPhysicalTraceCache::ghost_union_groups{0};
std::atomic<long long> BoundaryPhysicalTraceCache::requests_before{0};
std::atomic<long long> BoundaryPhysicalTraceCache::requests_after{0};
std::atomic<long long> BoundaryPhysicalTraceCache::union_routes{0};
std::atomic<long long> BoundaryPhysicalTraceCache::routing_operators{0};

BoundaryPhysicalTraceCache::BoundaryPhysicalTraceCache(
    const Mesh &mesh_, const mfem::Array<int> &bdr_attr_marker,
    std::shared_ptr<const FaceSamplingPlan> sampling_plan_)
  : mesh(&mesh_), marker(bdr_attr_marker.begin(), bdr_attr_marker.end()),
    sampling_plan(std::move(sampling_plan_))
{
  MFEM_VERIFY(sampling_plan, "Boundary physical traces require a sampling plan!");
  MFEM_VERIFY(sampling_plan->Matches(mesh_, bdr_attr_marker, sampling_plan->LOD()),
              "Boundary physical trace cache received an incompatible sampling plan!");
}

BoundaryPhysicalTraceCache::~BoundaryPhysicalTraceCache()
{
  for (auto &entry : entries)
  {
    DestroyUnionRoutes(entry.union_values);
  }
}

bool BoundaryPhysicalTraceCache::Supports(const mfem::ParFiniteElementSpace &fespace)
{
  const auto map_type = fespace.FEColl()->GetMapType(fespace.GetParMesh()->Dimension());
  // Integral-mapped 2D magnetic fields have a scalar H(div) representation.  They use
  // the established evaluator path for now; forcing them through a vector trace would
  // change the current-density/energy QFunction contract.
  return map_type == mfem::FiniteElement::H_CURL ||
         map_type == mfem::FiniteElement::H_DIV || map_type == mfem::FiniteElement::VALUE;
}

PointFieldKind
BoundaryPhysicalTraceCache::SourceKind(const mfem::ParFiniteElementSpace &fespace)
{
  const auto map_type = fespace.FEColl()->GetMapType(fespace.GetParMesh()->Dimension());
  if (map_type == mfem::FiniteElement::H_CURL)
  {
    return PointFieldKind::FIELD_E;
  }
  if (map_type == mfem::FiniteElement::H_DIV)
  {
    return PointFieldKind::FIELD_B;
  }
  MFEM_VERIFY(map_type == mfem::FiniteElement::VALUE,
              "Unsupported scalar physical trace source map type!");
  return PointFieldKind::FIELD_H1;
}

int BoundaryPhysicalTraceCache::NumComponents(
    const mfem::ParFiniteElementSpace &fespace) const
{
  const auto kind = SourceKind(fespace);
  return kind == PointFieldKind::FIELD_H1 ? 1 : fespace.GetParMesh()->SpaceDimension();
}

BoundaryPhysicalTraceCache::Entry &
BoundaryPhysicalTraceCache::GetEntry(const mfem::ParFiniteElementSpace &fespace)
{
  for (auto &entry : entries)
  {
    if (entry.fespace == &fespace)
    {
      return entry;
    }
  }
  MFEM_VERIFY(Supports(fespace), "Unsupported finite element space for physical trace!");

  mfem::Array<int> bdr_attr_marker(static_cast<int>(marker.size()));
  for (int i = 0; i < bdr_attr_marker.Size(); i++)
  {
    bdr_attr_marker[i] = marker[static_cast<std::size_t>(i)];
  }

  Entry entry;
  entry.fespace = &fespace;
  entry.kind = SourceKind(fespace);
  entry.num_comp = NumComponents(fespace);
  // Each trace producer evaluates exactly one physical side. Its output restriction
  // retains canonical q order; the two buffers therefore never average/jump a face.
  entry.side_a_eval.reset(new SurfaceFunctional(entry.kind, *mesh, bdr_attr_marker, fespace,
                                                sampling_plan->LOD(), sampling_plan,
                                                nullptr, 0));
  entry.side_b_eval.reset(new SurfaceFunctional(entry.kind, *mesh, bdr_attr_marker, fespace,
                                                sampling_plan->LOD(), sampling_plan,
                                                nullptr, 1));
  MFEM_VERIFY(entry.side_a_eval->IsValid() && entry.side_b_eval->IsValid(),
              "Unable to assemble boundary physical trace evaluator!");
  trace_groups += static_cast<long long>(entry.side_a_eval->TraceGroupCount() +
                                         entry.side_b_eval->TraceGroupCount());
  BuildUnionRoutes(entry);
  entries.emplace_back(std::move(entry));
  return entries.back();
}

void BoundaryPhysicalTraceCache::BuildUnionRoutes(Entry &entry)
{
  const auto &plan_groups = sampling_plan->UnionGroups();

  // Union planning is rank-local, but construction of the ghost exchange is
  // collective. Ranks with no local union groups must not return before agreeing
  // whether any rank needs that exchange.
  const mfem::ParMesh &pmesh = mesh->Get();
  const int dim = pmesh.Dimension();
  const int sdim = pmesh.SpaceDimension();
  const int num_comp = entry.num_comp;
  const int trace_points = sampling_plan->NumPoints();
  std::vector<int> local_base(plan_groups.size(), -1);
  std::vector<int> ghost_request(plan_groups.size(), -1);
  int local_points = 0;
  int local_groups = 0, ghost_groups = 0;
  long long selected_leaf_points = 0, selected_union_points = 0;
  long long selected_routes = 0, selected_requests_before = 0;
  for (std::size_t g = 0; g < plan_groups.size(); g++)
  {
    const auto &group = plan_groups[g];
    const long long group_leaf_points = [&]()
    {
      long long count = 0;
      for (const auto &route : group.routes)
      {
        count += static_cast<long long>(route.canonical_to_union.size());
      }
      return count;
    }();
    selected_leaf_points += group_leaf_points;
    selected_union_points += static_cast<long long>(group.points.size());
    selected_routes += group_leaf_points;
    if (group.ghost)
    {
      ghost_groups++;
      selected_requests_before += static_cast<long long>(group.routes.size());
    }
    else
    {
      local_base[g] = local_points;
      local_points += static_cast<int>(group.points.size());
      local_groups++;
    }
  }
  if (local_points > 0)
  {
    entry.union_values.local_values.SetSize(local_points * num_comp);
    entry.union_values.local_values.UseDevice(true);
    entry.union_values.local_values = 0.0;
  }

  // Owners with the same geometry and exact union rule share one ordinary fixed-rule
  // operator. Grouping is essential here: one operator per coarse owner would replace
  // leaf fragmentation with owner fragmentation on large NC meshes. The output remains
  // component-major, then the small EVAL_NONE routes below duplicate values into each
  // canonical leaf trace tuple.
  std::map<std::vector<long long>, std::vector<std::size_t>> local_patterns;
  for (std::size_t g = 0; g < plan_groups.size(); g++)
  {
    const auto &group = plan_groups[g];
    if (!group.ghost)
    {
      MFEM_VERIFY(group.owner >= 0 && group.owner < pmesh.GetNE() && !group.points.empty(),
                  "Invalid local NC coarse-face union owner!");
      local_patterns[MakeUnionRuleKey(group.geometry, group.points)].push_back(g);
    }
  }
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();
  const auto probe_info = GetProbeQFunction(*entry.fespace, dim, sdim);
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  for (const auto &[pattern, group_indices] : local_patterns)
  {
    (void)pattern;
    MFEM_VERIFY(!group_indices.empty(), "Empty NC union evaluator pattern!");
    const auto &representative = plan_groups[group_indices.front()];
    const int nq = static_cast<int>(representative.points.size());
    std::vector<int> owners;
    std::vector<CeedInt> out_offsets;
    owners.reserve(group_indices.size());
    out_offsets.reserve(group_indices.size() * static_cast<std::size_t>(nq));
    for (const std::size_t g : group_indices)
    {
      const auto &group = plan_groups[g];
      MFEM_VERIFY(group.geometry == representative.geometry &&
                      group.points.size() == representative.points.size(),
                  "Inconsistent NC union evaluator pattern!");
      owners.push_back(group.owner);
      for (int q = 0; q < nq; q++)
      {
        out_offsets.push_back(local_base[g] + q);
      }
    }

    const mfem::IntegrationRule *ir =
        GetRegisteredUnionIr(representative.geometry, representative.points);
    CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
        mesh_fespace, ceed, representative.geometry, owners, /*is_interp*/ true);
    const mfem::FiniteElement *mesh_fe =
        mesh_fespace.FEColl()->FiniteElementForGeometry(representative.geometry);
    MFEM_VERIFY(mesh_fe, "Unable to get mesh finite element for NC union route!");
    CeedBasis mesh_basis;
    ceed::InitCachedBasisFromRule(*mesh_fe, *ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
    CeedVector mesh_nodes_vec;
    ceed::InitCeedVector(*pmesh.GetNodes(), ceed, &mesh_nodes_vec);

    CeedElemRestriction field_restr;
    ceed::InitRestriction(*entry.fespace, owners, false, /*is_interp*/ true, false, ceed,
                          &field_restr);
    const mfem::FiniteElement *field_fe =
        entry.fespace->FEColl()->FiniteElementForGeometry(representative.geometry);
    MFEM_VERIFY(field_fe, "Unable to get field finite element for NC union route!");
    CeedBasis field_basis;
    ceed::InitCachedBasisFromRule(*field_fe, *ir, entry.fespace->GetVDim(), ceed,
                                  &field_basis);
    CeedVector field_vec;
    PalaceCeedCall(ceed,
                   CeedVectorCreate(ceed, static_cast<CeedSize>(entry.fespace->GetVSize()),
                                    &field_vec));

    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionCreate(
                       ceed, static_cast<CeedInt>(owners.size()), nq, num_comp,
                       local_points, static_cast<CeedSize>(local_points) * num_comp,
                       CEED_MEM_HOST, CEED_COPY_VALUES, out_offsets.data(), &out_restr));
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    inputs.push_back({"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
    inputs.push_back({"u_1", field_vec, field_restr, field_basis, ceed::EvalMode::Interp});
    CeedOperator op;
    ceed::AssembleCeedPointEvaluator(probe_info, nullptr, 0, ceed, inputs, num_comp,
                                     out_restr, &op);
    entry.union_values.local_evaluators.push_back({ceed, op, {{"u_1", 0}}});
    auto &evaluator = entry.union_values.local_evaluators.back();
    evaluator.mesh_nodes = pmesh.GetNodes();
    evaluator.mesh_node_fields = {"grad_x"};
    // The field vector was intentionally created without construction storage.
    evaluator.field_vectors_detached = true;
    fem::CacheGroupOperatorFieldVectors(evaluator);

    PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
    PalaceCeedCall(ceed, CeedVectorDestroy(&field_vec));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_restr));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
    PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
    PalaceCeedCall(ceed, CeedBasisDestroy(&field_basis));
  }

  // One ghost request now serves all leaf tuples owned by that ghost coarse element.
  // This exchange is deliberately constructed after the regular trace producers, and
  // applied after them, preserving the collective call order of the existing paths.
  // HasGhostUnion is communicator-wide, so ranks with no local request still construct
  // the empty peer exchange instead of returning before a collective constructor.
  if (sampling_plan->HasGhostUnion())
  {
    std::array<const mfem::ParFiniteElementSpace *, FaceNbrFieldExchange::MaxSources>
        fespaces = {entry.fespace, nullptr, nullptr, nullptr};
    std::vector<FaceNbrFieldExchange::Request> requests;
    for (std::size_t g = 0; g < plan_groups.size(); g++)
    {
      const auto &group = plan_groups[g];
      if (!group.ghost)
      {
        continue;
      }
      ghost_request[g] = static_cast<int>(requests.size());
      auto &request = requests.emplace_back();
      request.face_nbr_elem = group.owner;
      request.source_mask = 0b1u;
      // The union identity carries only geometry and normalized coordinate bits. The
      // fixed unit weights in group.points are a transport detail, not a leaf identity.
      request.point_key.reserve(2 + 3 * group.points.size());
      request.point_key.push_back(static_cast<long long>(group.geometry));
      request.point_key.push_back(static_cast<long long>(group.points.size()));
      for (const auto &point : group.points)
      {
        request.point_key.push_back(EncodeDouble(point.x));
        request.point_key.push_back(EncodeDouble(point.y));
        request.point_key.push_back(EncodeDouble(point.z));
      }
      request.pts = group.points;
    }
    entry.union_values.ghost_exchange =
        std::make_unique<FaceNbrFieldExchange>(*mesh, fespaces, requests);
  }

  // Route all owners sharing the same input layout, side, and leaf rule through one
  // EVAL_NONE operator. Otherwise a route launch per coarse owner would undo the union
  // evaluator grouping above.
  struct RouteRef
  {
    const FaceSamplingPlan::UnionGroup *group;
    const FaceSamplingPlan::UnionRoute *route;
    int input_base;
  };
  using RoutingKey = std::tuple<bool, int, int>;  // Imported input, side, leaf nq.
  std::map<RoutingKey, std::vector<RouteRef>> routing_patterns;
  for (std::size_t g = 0; g < plan_groups.size(); g++)
  {
    const auto &group = plan_groups[g];
    int input_base = local_base[g];
    if (group.ghost)
    {
      MFEM_VERIFY(entry.union_values.ghost_exchange && ghost_request[g] >= 0,
                  "Missing NC coarse-face ghost union exchange!");
      input_base = entry.union_values.ghost_exchange->ImportOffset(ghost_request[g], 0);
    }
    for (const auto &route : group.routes)
    {
      const int nq = static_cast<int>(route.canonical_to_union.size());
      MFEM_VERIFY(nq > 0, "Empty NC coarse-face union route!");
      routing_patterns[{group.ghost, group.side, nq}].push_back(
          {&group, &route, input_base});
    }
  }

  for (const auto &[routing_key, routes] : routing_patterns)
  {
    const auto [imported_input, side, nq] = routing_key;
    const int input_size = imported_input ? entry.union_values.ghost_exchange->ImportSize()
                                          : entry.union_values.local_values.Size();
    const int input_component_stride = imported_input ? 1 : local_points;
    std::vector<CeedInt> in_offsets(routes.size() * nq), out_offsets(routes.size() * nq);
    for (std::size_t e = 0; e < routes.size(); e++)
    {
      const auto &ref = routes[e];
      const auto &sampling_entry =
          sampling_plan->Entries()[static_cast<std::size_t>(ref.route->entry)];
      const int trace_base = sampling_plan->BufferBases()[sampling_entry.bdr_elem];
      MFEM_VERIFY(trace_base >= 0, "Invalid NC coarse-face trace route base!");
      for (int q = 0; q < nq; q++)
      {
        const int union_q = ref.route->canonical_to_union[static_cast<std::size_t>(q)];
        MFEM_VERIFY(union_q >= 0 && union_q < static_cast<int>(ref.group->points.size()),
                    "Invalid NC coarse-face union route index!");
        in_offsets[e * nq + q] =
            imported_input ? ref.input_base + num_comp * union_q : ref.input_base + union_q;
        out_offsets[e * nq + q] = trace_base + q;
      }
    }

    Entry::RoutingOperator route_op;
    route_op.ceed = ceed;
    route_op.imported_input = imported_input;
    route_op.side = side;
    CeedElemRestriction in_restr, out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, static_cast<CeedInt>(routes.size()), nq, num_comp,
                             input_component_stride, input_size, CEED_MEM_HOST,
                             CEED_COPY_VALUES, in_offsets.data(), &in_restr));
    PalaceCeedCall(ceed,
                   CeedElemRestrictionCreate(
                       ceed, static_cast<CeedInt>(routes.size()), nq, num_comp,
                       trace_points, static_cast<CeedSize>(trace_points) * num_comp,
                       CEED_MEM_HOST, CEED_COPY_VALUES, out_offsets.data(), &out_restr));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, input_size, &route_op.in_vec));
    PalaceCeedCall(ceed,
                   CeedVectorCreate(ceed, static_cast<CeedSize>(trace_points) * num_comp,
                                    &route_op.out_vec));
    CeedQFunction qf;
    PalaceCeedCall(ceed, CeedQFunctionCreateIdentity(ceed, num_comp, CEED_EVAL_NONE,
                                                     CEED_EVAL_NONE, &qf));
    PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf, nullptr, nullptr, &route_op.op));
    PalaceCeedCall(ceed, CeedQFunctionDestroy(&qf));
    PalaceCeedCall(ceed, CeedOperatorSetField(route_op.op, "input", in_restr,
                                              CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));
    PalaceCeedCall(ceed, CeedOperatorSetField(route_op.op, "output", out_restr,
                                              CEED_BASIS_NONE, CEED_VECTOR_ACTIVE));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(route_op.op));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&in_restr));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
    entry.union_values.routes.push_back(std::move(route_op));
  }

  leaf_points += selected_leaf_points;
  union_points += selected_union_points;
  duplicate_points += selected_leaf_points - selected_union_points;
  local_union_groups += local_groups;
  local_union_operators +=
      static_cast<long long>(entry.union_values.local_evaluators.size());
  ghost_union_groups += ghost_groups;
  requests_before += selected_requests_before;
  requests_after += ghost_groups;
  union_routes += selected_routes;
  routing_operators += static_cast<long long>(entry.union_values.routes.size());
  trace_groups += static_cast<long long>(entry.union_values.local_evaluators.size());
}

void BoundaryPhysicalTraceCache::ApplyUnionRoutes(Entry &entry, const Vector &u,
                                                  Entry::SourceValues &values) const
{
  auto &union_values = entry.union_values;
  if (!union_values.local_evaluators.empty())
  {
    union_values.local_values = 0.0;
    fem::ApplyAddGroupOperators(union_values.local_evaluators, {&u},
                                union_values.local_values);
  }
  if (union_values.ghost_exchange)
  {
    union_values.ghost_exchange->Exchange({&u, nullptr, nullptr, nullptr});
  }
  for (auto &route : union_values.routes)
  {
    const Vector &input = route.imported_input ? union_values.ghost_exchange->Imported()
                                               : union_values.local_values;
    Vector &output = route.side == 0 ? values.side_a : values.side_b;
    // Routing vectors are created unbound. Detach after every apply below so the
    // next source (real/imaginary or a new save generation) can attach directly.
    ceed::InitCeedVector(input, route.ceed, &route.in_vec, false, false);
    CeedMemType mem;
    PalaceCeedCall(route.ceed, CeedGetPreferredMemType(route.ceed, &mem));
    if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
    {
      mem = CEED_MEM_HOST;
    }
    double *output_data = output.ReadWrite(mem == CEED_MEM_DEVICE);
    PalaceCeedCall(route.ceed,
                   CeedVectorSetArray(route.out_vec, mem, CEED_USE_POINTER, output_data));
    PalaceCeedCall(route.ceed, CeedOperatorApplyAdd(route.op, route.in_vec, route.out_vec,
                                                    CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(route.ceed, CeedVectorTakeArray(route.in_vec, mem, nullptr));
    PalaceCeedCall(route.ceed, CeedVectorTakeArray(route.out_vec, mem, nullptr));
  }
}

void BoundaryPhysicalTraceCache::DestroyUnionRoutes(Entry::UnionValues &union_values)
{
  fem::DestroyGroupOperators(union_values.local_evaluators);
  for (auto &route : union_values.routes)
  {
    if (route.in_vec)
    {
      PalaceCeedCall(route.ceed, CeedVectorDestroy(&route.in_vec));
    }
    if (route.out_vec)
    {
      PalaceCeedCall(route.ceed, CeedVectorDestroy(&route.out_vec));
    }
    if (route.op)
    {
      PalaceCeedCall(route.ceed, CeedOperatorDestroy(&route.op));
    }
  }
  union_values.routes.clear();
  union_values.ghost_exchange.reset();
  union_values.local_values.Destroy();
}

void BoundaryPhysicalTraceCache::Invalidate()
{
  ++generation;
  // Do not retain arbitrary solver scratch addresses across saves. This both bounds the
  // cache and makes a same-address Vector unambiguously belong to the new generation.
  for (auto &entry : entries)
  {
    entry.source_values.clear();
  }
}

void BoundaryPhysicalTraceCache::VerifyCurrent() const
{
  mfem::Array<int> bdr_attr_marker(static_cast<int>(marker.size()));
  for (int i = 0; i < bdr_attr_marker.Size(); i++)
  {
    bdr_attr_marker[i] = marker[static_cast<std::size_t>(i)];
  }
  MFEM_VERIFY(sampling_plan->Matches(*mesh, bdr_attr_marker, sampling_plan->LOD()),
              "Boundary physical trace cache is stale after a mesh update; rebuild its "
              "sampling plan and evaluators!");
}

BoundaryPhysicalTraceCache::Traces
BoundaryPhysicalTraceCache::Get(const mfem::ParFiniteElementSpace &fespace, const Vector &u)
{
  VerifyCurrent();
  auto &entry = GetEntry(fespace);
  auto [it, inserted] = entry.source_values.emplace(&u, Entry::SourceValues{});
  auto &values = it->second;
  if (!inserted && values.generation == generation)
  {
    ++cache_hits;
    return {&values.side_a, &values.side_b, entry.num_comp};
  }

  values.side_a.SetSize(sampling_plan->NumPoints() * entry.num_comp);
  values.side_b.SetSize(sampling_plan->NumPoints() * entry.num_comp);
  values.side_a.UseDevice(true);
  values.side_b.UseDevice(true);
  entry.side_a_eval->EvalBuffer(u, values.side_a);
  entry.side_b_eval->EvalBuffer(u, values.side_b);
  ApplyUnionRoutes(entry, u, values);
  values.generation = generation;
  operator_applies += static_cast<long long>(
      entry.side_a_eval->TraceGroupCount() + entry.side_b_eval->TraceGroupCount() +
      entry.union_values.local_evaluators.size() + entry.union_values.routes.size());
  const long long points = 2LL * sampling_plan->NumPoints();
  evaluated_points += points;
  evaluated_bytes += points * entry.num_comp * static_cast<long long>(sizeof(double));
  return {&values.side_a, &values.side_b, entry.num_comp};
}

const FaceSamplingPlan &BoundaryPhysicalTraceCache::SamplingPlan() const
{
  return *sampling_plan;
}

int BoundaryPhysicalTraceCache::NumPoints() const
{
  return sampling_plan->NumPoints();
}

long long BoundaryPhysicalTraceCache::TraceGroupCount()
{
  return trace_groups.load();
}

long long BoundaryPhysicalTraceCache::OperatorApplyCount()
{
  return operator_applies.load();
}

long long BoundaryPhysicalTraceCache::CacheHitCount()
{
  return cache_hits.load();
}

long long BoundaryPhysicalTraceCache::EvaluatedPointCount()
{
  return evaluated_points.load();
}

long long BoundaryPhysicalTraceCache::EvaluatedByteCount()
{
  return evaluated_bytes.load();
}

long long BoundaryPhysicalTraceCache::LeafPointCount()
{
  return leaf_points.load();
}

long long BoundaryPhysicalTraceCache::UnionPointCount()
{
  return union_points.load();
}

long long BoundaryPhysicalTraceCache::DuplicatePointCount()
{
  return duplicate_points.load();
}

long long BoundaryPhysicalTraceCache::LocalUnionGroupCount()
{
  return local_union_groups.load();
}

long long BoundaryPhysicalTraceCache::LocalUnionOperatorCount()
{
  return local_union_operators.load();
}

long long BoundaryPhysicalTraceCache::GhostUnionGroupCount()
{
  return ghost_union_groups.load();
}

long long BoundaryPhysicalTraceCache::RequestCountBefore()
{
  return requests_before.load();
}

long long BoundaryPhysicalTraceCache::RequestCountAfter()
{
  return requests_after.load();
}

long long BoundaryPhysicalTraceCache::UnionRouteCount()
{
  return union_routes.load();
}

long long BoundaryPhysicalTraceCache::RoutingOperatorCount()
{
  return routing_operators.load();
}

void BoundaryPhysicalTraceCache::PrintProfile() const
{
  if (!std::getenv("PALACE_BOUNDARY_TRACE_PROFILE"))
  {
    return;
  }
  Mpi::Print(mesh->GetComm(),
             "BoundaryPhysicalTrace profile trace_groups={} operator_applies={} "
             "cache_hits={} evaluated_points={} evaluated_bytes={} generation={} "
             "leaf_points={} union_points={} duplicates_elided={} local_union_groups={} "
             "local_union_operators={} ghost_union_groups={} requests_before={} "
             "requests_after={} routes={} routing_operators={}\n",
             TraceGroupCount(), OperatorApplyCount(), CacheHitCount(),
             EvaluatedPointCount(), EvaluatedByteCount(), generation, LeafPointCount(),
             UnionPointCount(), DuplicatePointCount(), LocalUnionGroupCount(),
             LocalUnionOperatorCount(), GhostUnionGroupCount(), RequestCountBefore(),
             RequestCountAfter(), UnionRouteCount(), RoutingOperatorCount());
}

}  // namespace palace
