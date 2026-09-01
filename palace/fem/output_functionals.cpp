// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "output_functionals.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include "fem/coefficient.hpp"
#include "fem/facenbrexchange.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/21/surf_21_qf.h"
#include "fem/qfunctions/32/surf_32_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

// Key identifying a group of boundary elements that share a fixed mapped integration
// rule. Canonical mapped point clouds, not face orientation/local-face metadata, form
// the trace identity; exact point signatures keep the finite registry reusable.
using FaceConfigKey = std::vector<long long>;

long long EncodePointRuleDouble(double x)
{
  // Point rules are immutable fixed tabulations. Preserve every IEEE-754 bit, including
  // signed zero, so the registry can never reuse a rule whose supplied point data differ.
  static_assert(sizeof(long long) == sizeof(double));
  long long bits;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

void AppendPointRuleSignature(FaceConfigKey &key,
                              const std::vector<mfem::IntegrationPoint> &pts)
{
  // Coordinate and weight bits are ordered by quadrature point and component. This is
  // part of the process-lifetime registry identity, not merely a collision guard.
  key.push_back(static_cast<long long>(pts.size()));
  for (const auto &ip : pts)
  {
    key.push_back(EncodePointRuleDouble(ip.x));
    key.push_back(EncodePointRuleDouble(ip.y));
    key.push_back(EncodePointRuleDouble(ip.z));
    key.push_back(EncodePointRuleDouble(ip.weight));
  }
}

void VerifyMappedRuleMatches(const mfem::IntegrationRule &ir,
                             const std::vector<mfem::IntegrationPoint> &pts)
{
  MFEM_VERIFY(ir.GetNPoints() == static_cast<int>(pts.size()),
              "Canonical trace-rule key reused with a different number of mapped points!");
  for (int q = 0; q < ir.GetNPoints(); q++)
  {
    const auto &a = ir.IntPoint(q);
    const auto &b = pts[static_cast<std::size_t>(q)];
    const double err = std::max({std::abs(a.x - b.x), std::abs(a.y - b.y),
                                 std::abs(a.z - b.z), std::abs(a.weight - b.weight)});
    MFEM_VERIFY(err <= 1.0e-12,
                "Canonical trace-rule key reused for different mapped reference points; "
                "refine the finite MFEM/NCMesh face-map key before using this group!");
  }
}

// The face rules used here are symmetric under MFEM's finite set of face
// orientations. Sort the mapped volume points to remove only that q-order variation.
// The resulting canonical_q -> original_q map is retained for output routing, so this
// never changes the public boundary point ordering.
struct CanonicalMappedRule
{
  std::vector<mfem::IntegrationPoint> pts;
  std::vector<int> canonical_to_original;
};

void NormalizeReferencePoint(mfem::IntegrationPoint &ip)
{
  // A transformation can introduce -0.0 while its orientation counterpart produces
  // +0.0. They identify the same reference point, so do not let signed zero fragment a
  // finite fixed-rule registry.
  if (ip.x == 0.0)
  {
    ip.x = 0.0;
  }
  if (ip.y == 0.0)
  {
    ip.y = 0.0;
  }
  if (ip.z == 0.0)
  {
    ip.z = 0.0;
  }
  if (ip.weight == 0.0)
  {
    ip.weight = 0.0;
  }
}

bool CanonicalPointLess(const mfem::IntegrationPoint &a, const mfem::IntegrationPoint &b)
{
  if (a.x != b.x)
  {
    return a.x < b.x;
  }
  if (a.y != b.y)
  {
    return a.y < b.y;
  }
  if (a.z != b.z)
  {
    return a.z < b.z;
  }
  return a.weight < b.weight;
}

CanonicalMappedRule
CanonicalizeSymmetricMappedRule(const std::vector<mfem::IntegrationPoint> &original)
{
  CanonicalMappedRule rule;
  rule.canonical_to_original.resize(original.size());
  std::iota(rule.canonical_to_original.begin(), rule.canonical_to_original.end(), 0);
  for (const auto &ip : original)
  {
    MFEM_VERIFY(std::isfinite(ip.x) && std::isfinite(ip.y) && std::isfinite(ip.z) &&
                    std::isfinite(ip.weight),
                "Non-finite mapped reference point cannot be canonically routed!");
  }
  std::stable_sort(rule.canonical_to_original.begin(), rule.canonical_to_original.end(),
                   [&](int i, int j)
                   {
                     mfem::IntegrationPoint a = original[static_cast<std::size_t>(i)];
                     mfem::IntegrationPoint b = original[static_cast<std::size_t>(j)];
                     NormalizeReferencePoint(a);
                     NormalizeReferencePoint(b);
                     return CanonicalPointLess(a, b);
                   });
  rule.pts.resize(original.size());
  std::vector<bool> seen(original.size(), false);
  for (std::size_t q = 0; q < rule.canonical_to_original.size(); q++)
  {
    const int original_q = rule.canonical_to_original[q];
    MFEM_VERIFY(original_q >= 0 && original_q < static_cast<int>(original.size()) &&
                    !seen[static_cast<std::size_t>(original_q)],
                "Mapped face rule does not define a bijective canonical q permutation!");
    seen[static_cast<std::size_t>(original_q)] = true;
    rule.pts[q] = original[static_cast<std::size_t>(original_q)];
    NormalizeReferencePoint(rule.pts[q]);
  }
  for (const bool used : seen)
  {
    MFEM_VERIFY(used,
                "Mapped face rule does not define a bijective canonical q permutation!");
  }
  // This is the symmetric-face-rule contract used by grouping: every canonical point
  // is exactly an original mapped point and the permutation covers that point set once.
  // Equal canonical clouds may therefore share an ordinary fixed libCEED tabulation.
  for (std::size_t q = 0; q < rule.pts.size(); q++)
  {
    mfem::IntegrationPoint original_ip =
        original[static_cast<std::size_t>(rule.canonical_to_original[q])];
    NormalizeReferencePoint(original_ip);
    const auto &canonical_ip = rule.pts[q];
    MFEM_VERIFY(canonical_ip.x == original_ip.x && canonical_ip.y == original_ip.y &&
                    canonical_ip.z == original_ip.z &&
                    canonical_ip.weight == original_ip.weight,
                "Canonical face rule is not an exact permutation of mapped points!");
  }
  return rule;
}

std::vector<int> IdentityPointPermutation(int nq)
{
  std::vector<int> perm(nq);
  std::iota(perm.begin(), perm.end(), 0);
  return perm;
}

// Registry of mapped face integration rules. The IntegrationRule objects must have
// application lifetime: mfem::FiniteElement::GetDofToQuad caches tabulations keyed by
// the IntegrationRule pointer inside the (global, shared) FiniteElement objects, so
// destroying an IntegrationRule which was used for tabulation would leave a dangling
// cache entry.
const mfem::IntegrationRule *
GetRegisteredMappedIr(const FaceConfigKey &key,
                      const std::vector<mfem::IntegrationPoint> &pts)
{
  FaceConfigKey registry_key = key;
  AppendPointRuleSignature(registry_key, pts);
  static std::map<FaceConfigKey, std::unique_ptr<mfem::IntegrationRule>> registry;
  static std::mutex registry_mutex;
  std::lock_guard<std::mutex> lock(registry_mutex);
  auto it = registry.find(registry_key);
  if (it == registry.end())
  {
    auto ir = std::make_unique<mfem::IntegrationRule>(static_cast<int>(pts.size()));
    for (std::size_t q = 0; q < pts.size(); q++)
    {
      ir->IntPoint(static_cast<int>(q)) = pts[q];
    }
    it = registry.emplace(std::move(registry_key), std::move(ir)).first;
  }
  else
  {
    VerifyMappedRuleMatches(*it->second, pts);
  }
  return it->second.get();
}

// Evaluation plan for a single marked boundary element: which volume element(s) the
// field is evaluated from, and the face quadrature point positions mapped into the
// volume element reference space(s). elem_b >= 0 indicates two-sided evaluation.
struct ElemPlan
{
  int bdr_elem;
  int elem_a = -1, elem_b = -1;
  std::vector<mfem::IntegrationPoint> pts_a, pts_b;
  bool flip = false;  // Boundary element and face have the same orientation (o % 2 == 0),
                      // legacy coefficients invert the boundary element normal
  // A side whose volume element is a face neighbor (ghost) on another process: its field
  // is not in the local vector, so it is pulled via FaceNbrFieldExchange and fed to the
  // operator at the mapped points (only Elem2 can be a ghost, so at most one side).
  bool ghost_a = false, ghost_b = false;
  int face_nbr = -1;   // Face neighbor element index (Elem2No - ParMesh::GetNE())
  int ghost_attr = 0;  // Ghost element mesh attribute (material lookup on the requester)
  mfem::Geometry::Type ghost_geom = mfem::Geometry::INVALID;  // Ghost element geometry
  // Canonical identities of ghost-side point clouds: volume/face geometry plus exact
  // canonical mapped points and weights. These are passed to FaceNbrFieldExchange so
  // remote fixed-rule evaluation uses the same canonical q order.
  std::vector<long long> point_key_a, point_key_b;
};

// Data for one group of boundary elements sharing the same face configuration.
struct FaceGroup
{
  mfem::Geometry::Type bdr_geom;
  mfem::Geometry::Type vol_geom_a = mfem::Geometry::INVALID;
  mfem::Geometry::Type vol_geom_b = mfem::Geometry::INVALID;
  // Separable two-sided operations use two one-sided fixed-rule operators which
  // ApplyAdd into the same output slots. normal_scale is the effective oriented sign,
  // including the boundary-element flip when relevant.
  double normal_scale = 1.0;
  const mfem::IntegrationRule *mapped_ir_a = nullptr;  // Registry, application lifetime
  const mfem::IntegrationRule *mapped_ir_b = nullptr;
  std::vector<int> bdr_indices, vol_indices_a, vol_indices_b;
  std::vector<int> out_slots;  // Output vector slot for each boundary element
  // Per entry, maps canonical mapped-rule q to the legacy boundary-rule q. It is
  // encoded directly in the EVAL_NONE output restriction for pointwise buffers.
  std::vector<std::vector<int>> canonical_to_original;
  // Ghost (face neighbor) side, if any: at most one of ghost_a/ghost_b is set (only
  // Elem2 can be a ghost). For the ghost side, vol_indices are empty; face_nbr and
  // ghost_attr hold the per-element face neighbor index and mesh attribute, and req_idx
  // the FaceNbrFieldExchange request index (filled when the exchange is built).
  bool ghost_a = false, ghost_b = false;
  std::vector<int> face_nbr, ghost_attr, req_idx;
  std::vector<std::vector<long long>> face_nbr_point_keys;
  std::vector<std::vector<mfem::IntegrationPoint>> face_nbr_points;
};

long long EncodeDiscreteHalfScale(double scale)
{
  // Separable face operations use only +/-1 and +/-1/2. Keep that exact discrete state
  // in the group key instead of allowing rounded floating-point scale values to alias.
  const long long encoded = std::llround(2.0 * scale);
  MFEM_VERIFY(std::abs(2.0 * scale - static_cast<double>(encoded)) <= 1.0e-12,
              "Unexpected non-half-integer face side/normal scale!");
  return encoded;
}

std::vector<long long>
MakeCanonicalTracePointKey(mfem::Geometry::Type vol_geom, mfem::Geometry::Type face_geom,
                           const mfem::IntegrationRule &face_ir,
                           const std::vector<mfem::IntegrationPoint> &canonical_pts)
{
  // Reference coordinates and weights are the complete, finite identity of a mapped
  // local face/subface rule. In particular, local-face and orientation ids are excluded:
  // they are represented by the canonical cloud and (for pointwise output) its routing
  // permutation. Distinct local faces and NC subface maps remain distinct because their
  // mapped volume point clouds differ.
  FaceConfigKey key;
  key.reserve(4 + 4 * canonical_pts.size());
  key.push_back(static_cast<long long>(vol_geom));
  key.push_back(static_cast<long long>(face_geom));
  key.push_back(static_cast<long long>(face_ir.GetOrder()));
  AppendPointRuleSignature(key, canonical_pts);
  return key;
}

// Holds libCEED object references created during operator assembly for destruction once
// the assembled operator owns them.
struct CeedAssemblyScratch
{
  Ceed ceed;
  std::vector<CeedVector> vecs;
  std::vector<CeedElemRestriction> restrs;
  std::vector<CeedBasis> bases;

  CeedAssemblyScratch(Ceed ceed) : ceed(ceed) {}
  CeedAssemblyScratch(const CeedAssemblyScratch &) = delete;
  CeedAssemblyScratch &operator=(const CeedAssemblyScratch &) = delete;

  ~CeedAssemblyScratch()
  {
    for (auto &v : vecs)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
    for (auto &r : restrs)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
    }
    for (auto &b : bases)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&b));
    }
  }
};

}  // namespace

SurfaceFunctional::KernelKind SurfaceFunctional::ToKernelKind(Kind kind)
{
  switch (kind)
  {
    case Kind::AREA:
      return KernelKind::AREA;
    case Kind::HCURL_NORM2:
      return KernelKind::HCURL_NORM2;
    case Kind::INTERFACE_EPR:
      return KernelKind::INTERFACE_EPR;
    case Kind::SURFACE_FLUX:
      return KernelKind::SURFACE_FLUX;
    case Kind::FARFIELD:
      return KernelKind::FARFIELD;
  }
  MFEM_ABORT("Unknown surface functional kind!");
}

const char *SurfaceFunctional::KindName(KernelKind kind)
{
  switch (kind)
  {
    case KernelKind::AREA:
      return "AREA";
    case KernelKind::HCURL_NORM2:
      return "HCURL_NORM2";
    case KernelKind::INTERFACE_EPR:
      return "INTERFACE_EPR";
    case KernelKind::SURFACE_FLUX:
      return "SURFACE_FLUX";
    case KernelKind::FARFIELD:
      return "FARFIELD";
    case KernelKind::MODE_OVERLAP:
      return "MODE_OVERLAP";
  }
  return "UNKNOWN";
}

SurfaceFunctional::SurfaceFunctional(Kind kind, const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace *fespace)
  : kind(ToKernelKind(kind)), nd_fespace(fespace), rt_fespace(nullptr), mat_op(nullptr),
    comm(mesh.GetComm())
{
  MFEM_VERIFY(kind == Kind::AREA || kind == Kind::HCURL_NORM2,
              "Invalid SurfaceFunctional constructor for the requested functional kind!");
  MFEM_VERIFY(kind == Kind::AREA || fespace,
              "SurfaceFunctional requires a field finite element space for functionals "
              "with field inputs!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const MaterialOperator &mat_op,
                                     InterfaceDielectric type, double t_i, double epsilon_i)
  : kind(KernelKind::INTERFACE_EPR), epr_type(type), epr_t(t_i), epr_epsilon(epsilon_i),
    nd_fespace(&nd_fespace), rt_fespace(nullptr), mat_op(&mat_op), comm(mesh.GetComm())
{
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace *nd_fespace,
                                     const mfem::ParFiniteElementSpace *rt_fespace,
                                     const MaterialOperator &mat_op, SurfaceFlux type,
                                     bool two_sided, const mfem::Vector &x0)
  : kind(KernelKind::SURFACE_FLUX), flux_type(type), flux_two_sided(two_sided), flux_x0(x0),
    nd_fespace(nd_fespace), rt_fespace(rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  MFEM_VERIFY(
      (nd_fespace || (type != SurfaceFlux::ELECTRIC && type != SurfaceFlux::POWER)) &&
          (rt_fespace || (type != SurfaceFlux::MAGNETIC && type != SurfaceFlux::POWER)),
      "Missing finite element space for surface flux functional!");
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::Array<int> &bdr_attr_marker,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const mfem::ParFiniteElementSpace &rt_fespace,
                                     const MaterialOperator &mat_op,
                                     const std::vector<std::array<double, 3>> &r_naughts)
  : kind(KernelKind::FARFIELD), farfield_dirs(r_naughts), nd_fespace(&nd_fespace),
    rt_fespace(&rt_fespace), mat_op(&mat_op), comm(mesh.GetComm())
{
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::SurfaceFunctional(const Mesh &mesh,
                                     const mfem::ParFiniteElementSpace &nd_fespace,
                                     const std::vector<SurfaceModeCoefficient> &mode_coeffs)
  : kind(KernelKind::MODE_OVERLAP), nd_fespace(&nd_fespace), rt_fespace(nullptr),
    mat_op(nullptr), comm(mesh.GetComm())
{
  const mfem::ParMesh &pmesh = mesh.Get();
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker(bdr_attr_max);
  bdr_attr_marker = 0;
  for (const auto &coeff : mode_coeffs)
  {
    for (int attr : coeff.attr_list)
    {
      MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                  "Invalid boundary attribute for surface mode-overlap functional!");
      bdr_attr_marker[attr - 1] = 1;
      const bool inserted = mode_coeff_by_attr.emplace(attr, coeff).second;
      MFEM_VERIFY(inserted,
                  "Surface mode-overlap functional does not support overlapping mode "
                  "coefficient attributes!");
    }
  }
  Assemble(mesh, bdr_attr_marker);
}

SurfaceFunctional::~SurfaceFunctional()
{
  fem::DestroyGroupOperators(groups);
}

void SurfaceFunctional::Assemble(const Mesh &mesh, const mfem::Array<int> &bdr_attr_marker)
{
  AssembleLocal(mesh, bdr_attr_marker);

  // The validity decision must be globally consistent: a rank may locally encounter an
  // unsupported configuration while others do not. All ranks must agree on whether the
  // libCEED path is valid so that collective evaluation calls match; model-level callers
  // decide whether invalid means an explicit unsupported-case legacy path or an error for
  // a supported path.
  bool global_valid = valid;
  Mpi::GlobalAnd(1, &global_valid, comm);
  if (!global_valid)
  {
    // Discard any partially assembled local state before the caller selects the explicit
    // unsupported-dimension fallback or reports an invalid supported configuration.
    fem::DestroyGroupOperators(groups);
    face_nbr_exchange.reset();
    elem_attrs.clear();
    field_staging.SetSize(0);
    local_out.SetSize(0);
    local_out_attrs.clear();
    valid = false;
  }
  // Passive field vectors are re-pointed at caller data on every apply. Detach their
  // borrowed construction arrays before releasing the staging allocation.
  fem::DetachGroupOperatorFieldVectors(groups);
  field_staging.Destroy();
  MFEM_ASSERT(field_staging.Capacity() == 0,
              "Surface functional staging allocation was not released!");
}

std::vector<CeedIntScalar> SurfaceFunctional::BuildBaseContext(int dim, bool is_2d) const
{
  // The reduction-only base currently has no 2D-specific contexts.
  (void)is_2d;
  std::vector<CeedIntScalar> base_ctx;
  if (kind == KernelKind::INTERFACE_EPR)
  {
    // CeedIntScalar is a union, so the runtime integrand selector (epr_type) needs its
    // own slot: [0].first = epr_type (0 = DEFAULT, 1 = MA, 2 = MS, 3 = SA), then
    // [1].second = scale0, [2].second = scale1, then (MS only) the material context. The
    // shared kernel passes ctx + 1 to the per-type helpers so their relative layout
    // (scale0, scale1, material) is unchanged.
    base_ctx.resize(3);
    base_ctx[1].second = 0.0;
    base_ctx[2].second = 0.0;
    switch (epr_type)
    {
      case InterfaceDielectric::DEFAULT:
        base_ctx[0].first = 0;
        base_ctx[1].second = 0.5 * epr_t * epr_epsilon;
        break;
      case InterfaceDielectric::MA:
        base_ctx[0].first = 1;
        base_ctx[1].second = 0.5 * epr_t / epr_epsilon;
        break;
      case InterfaceDielectric::MS:
        {
          base_ctx[0].first = 2;
          base_ctx[1].second = 0.5 * epr_t / epr_epsilon;
          MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                                   mat_op->GetPermittivityReal());
          auto mat_ctx = ceed::PopulateCoefficientContext(dim, &epsilon_func);
          base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
        }
        break;
      case InterfaceDielectric::SA:
        base_ctx[0].first = 3;
        base_ctx[1].second = 0.5 * epr_t * epr_epsilon;
        base_ctx[2].second = 0.5 * epr_t / epr_epsilon;
        break;
    }
  }
  else if (kind == KernelKind::SURFACE_FLUX)
  {
    base_ctx.resize(5);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].first = flux_two_sided;
    for (int d = 0; d < 3; d++)
    {
      base_ctx[2 + d].second = (flux_x0.Size() > d) ? flux_x0(d) : 0.0;
    }
    if (flux_type == SurfaceFlux::ELECTRIC)
    {
      MaterialPropertyCoefficient epsilon_func(mat_op->GetAttributeToMaterial(),
                                               mat_op->GetPermittivityReal());
      auto mat_ctx = ceed::PopulateCoefficientContext(3, &epsilon_func);
      base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
    }
    else if (flux_type == SurfaceFlux::POWER)
    {
      MaterialPropertyCoefficient invmu_func(mat_op->GetAttributeToMaterial(),
                                             mat_op->GetInvPermeability());
      auto mat_ctx = ceed::PopulateCoefficientContext(3, &invmu_func);
      base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
    }
  }
  else if (kind == KernelKind::MODE_OVERLAP)
  {
    // Group-specific entries are filled after selecting the boundary attribute:
    // [0].first = type, [1].second = scale, [2..4].second = direction/origin.
    base_ctx.resize(5);
  }
  else if (kind == KernelKind::FARFIELD)
  {
    const int N = static_cast<int>(farfield_dirs.size());
    const auto b_map_type = rt_fespace->FEColl()->GetMapType(dim);
    MFEM_VERIFY(b_map_type == mfem::FiniteElement::H_CURL ||
                    b_map_type == mfem::FiniteElement::H_DIV,
                "Far-field postprocessing requires an H(curl) or H(div) magnetic field "
                "space!");
    base_ctx.resize(5 + 3 * N);
    base_ctx[0].second = 1.0;  // Normal sign, set per group
    base_ctx[1].second = farfield_omega.real();
    base_ctx[2].second = farfield_omega.imag();
    base_ctx[3].first = N;
    base_ctx[4].first = (b_map_type == mfem::FiniteElement::H_DIV) ? 1 : 0;
    for (int d = 0; d < N; d++)
    {
      for (int c = 0; c < 3; c++)
      {
        base_ctx[5 + 3 * d + c].second = farfield_dirs[d][c];
      }
    }
    MaterialPropertyCoefficient c0_func(mat_op->GetAttributeToMaterial(),
                                        mat_op->GetLightSpeed());
    auto mat_ctx = ceed::PopulateCoefficientContext(3, &c0_func);
    base_ctx.insert(base_ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }
  else
  {
    base_ctx.resize(2);
    base_ctx[0].second = 0.0;
    base_ctx[1].second = 1.0;
  }

  return base_ctx;
}

void SurfaceFunctional::ConfigureGroupQFunction(bool is_2d, bool has_b, double normal_scale,
                                                int bdr_attr,
                                                const std::vector<CeedIntScalar> &base_ctx,
                                                std::vector<CeedIntScalar> &ctx,
                                                ceed::CeedQFunctionInfo &info) const
{
  ctx = base_ctx;
  switch (kind)
  {
    case KernelKind::AREA:
      info.apply_qf = is_2d ? f_integ_surf_area_21 : f_integ_surf_area_32;
      info.apply_qf_path = is_2d ? PalaceQFunctionRelativePath(f_integ_surf_area_21_loc)
                                 : PalaceQFunctionRelativePath(f_integ_surf_area_32_loc);
      break;
    case KernelKind::HCURL_NORM2:
      info.apply_qf = is_2d ? f_integ_surf_hcurl_norm2_21 : f_integ_surf_hcurl_norm2_32;
      info.apply_qf_path =
          is_2d ? PalaceQFunctionRelativePath(f_integ_surf_hcurl_norm2_21_loc)
                : PalaceQFunctionRelativePath(f_integ_surf_hcurl_norm2_32_loc);
      break;
    case KernelKind::INTERFACE_EPR:
      // All four interface types share one kernel; epr_type is set in base_ctx[0].first.
      info.apply_qf = is_2d ? (has_b ? f_integ_surf_epr_2_21 : f_integ_surf_epr_1_21)
                            : (has_b ? f_integ_surf_epr_2_32 : f_integ_surf_epr_1_32);
      info.apply_qf_path =
          is_2d ? PalaceQFunctionRelativePath(has_b ? f_integ_surf_epr_2_21_loc
                                                    : f_integ_surf_epr_1_21_loc)
                : PalaceQFunctionRelativePath(has_b ? f_integ_surf_epr_2_32_loc
                                                    : f_integ_surf_epr_1_32_loc);
      break;
    case KernelKind::FARFIELD:
      ctx[0].second = normal_scale;
      info.apply_qf = f_integ_surf_farfield_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_farfield_32_loc);
      break;
    case KernelKind::MODE_OVERLAP:
      {
        auto mode_it = mode_coeff_by_attr.find(bdr_attr);
        MFEM_VERIFY(mode_it != mode_coeff_by_attr.end(),
                    "Missing mode coefficient for marked boundary attribute!");
        const auto &mode = mode_it->second;
        ctx[0].first = (mode.type == SurfaceModeCoefficient::Type::COAXIAL) ? 1 : 0;
        ctx[1].second = mode.scale;
        const auto &data = (mode.type == SurfaceModeCoefficient::Type::COAXIAL)
                               ? mode.origin
                               : mode.direction;
        for (int d = 0; d < 3; d++)
        {
          ctx[2 + d].second = data[d];
        }
        info.apply_qf = f_integ_surf_mode_32;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_integ_surf_mode_32_loc);
      }
      break;
    case KernelKind::SURFACE_FLUX:
      // A split flux side uses the one-sided kernel and contributes through the
      // effective context normal scale, retaining difference/average semantics.
      ctx[0].second = normal_scale;
      switch (flux_type)
      {
        case SurfaceFlux::ELECTRIC:
          info.apply_qf = has_b ? f_integ_surf_flux_e_2_32 : f_integ_surf_flux_e_1_32;
          info.apply_qf_path = PalaceQFunctionRelativePath(
              has_b ? f_integ_surf_flux_e_2_32_loc : f_integ_surf_flux_e_1_32_loc);
          break;
        case SurfaceFlux::MAGNETIC:
          info.apply_qf = has_b ? f_integ_surf_flux_m_2_32 : f_integ_surf_flux_m_1_32;
          info.apply_qf_path = PalaceQFunctionRelativePath(
              has_b ? f_integ_surf_flux_m_2_32_loc : f_integ_surf_flux_m_1_32_loc);
          break;
        case SurfaceFlux::POWER:
          info.apply_qf = has_b ? f_integ_surf_flux_p_2_32 : f_integ_surf_flux_p_1_32;
          info.apply_qf_path = PalaceQFunctionRelativePath(
              has_b ? f_integ_surf_flux_p_2_32_loc : f_integ_surf_flux_p_1_32_loc);
          break;
      }
      break;
  }
}

void SurfaceFunctional::AssembleLocal(const Mesh &mesh,
                                      const mfem::Array<int> &bdr_attr_marker)
{
  const int dim = mesh.Dimension();
  const int sdim = mesh.SpaceDimension();
  const bool is_2d = (dim == 2 && sdim == 2);
  const bool is_3d = (dim == 3 && sdim == 3);
  if (!is_2d && !is_3d)
  {
    valid = false;
    return;
  }
  if (is_2d && kind != KernelKind::AREA && kind != KernelKind::HCURL_NORM2 &&
      kind != KernelKind::INTERFACE_EPR)
  {
    // Initial 2D support covers line integrals needed by interface dielectric
    // postprocessing. Other boundary-output, surface-flux, and mode-overlap kinds still
    // use the legacy paths until their 2D scalar/vector conventions are implemented.
    valid = false;
    return;
  }
  const mfem::ParMesh &pmesh = mesh.Get();
  MFEM_VERIFY(pmesh.GetNodes(), "The mesh has no nodal FE space!");
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();
  const bool need_field = (kind != KernelKind::AREA);
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  // Plan the evaluation for each marked boundary element and group it by its fixed
  // mapped face/subface transformation. The process-lifetime registry key also carries
  // the exact mapped point-rule signature, so equal canonical maps reuse tabulations
  // across evaluator rebuilds while distinct supplied rules never alias.
  std::map<FaceConfigKey, FaceGroup> face_groups;
  int num_marked = 0;
  {
    constexpr double threshold = 1.0 - 1.0e-6;
    mfem::FaceElementTransformations FET;
    mfem::IsoparametricTransformation T1, T2;
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const int attr = pmesh.GetBdrAttribute(i);
      if (!bdr_attr_marker[attr - 1])
      {
        continue;
      }

      // Get the face and attached element transformations following the same
      // conventions as the legacy BdrGridFunctionCoefficient evaluation path
      // (FET.Elem1 always exists and is local, FET.Elem2 exists for interior
      // boundaries and may correspond to a face neighbor element on another process
      // for boundaries on parallel interfaces).
      const bool flip = BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
          i, pmesh, FET, T1, T2);
      const bool has_elem2 = (FET.Elem2 != nullptr);
      const bool elem2_local = has_elem2 && (FET.Elem2No < pmesh.GetNE());

      // Decide which side(s) the field is evaluated from following the conventions of
      // the legacy coefficients (see InterfaceDielectricCoefficient and
      // BdrSurfaceFluxCoefficient).
      ElemPlan plan;
      plan.bdr_elem = i;
      plan.flip = flip;
      if (!need_field)
      {
        // No field inputs (side selection does not apply).
      }
      else if (kind == KernelKind::HCURL_NORM2 || kind == KernelKind::MODE_OVERLAP)
      {
        plan.elem_a = FET.Elem1No;
      }
      else if (kind == KernelKind::FARFIELD)
      {
        if (has_elem2)
        {
          // Far-field computations are only supported on external boundaries. Record
          // local invalidity but keep participating in subsequent collectives; another
          // rank may own no element of this interior surface.
          valid = false;
          continue;
        }
        plan.elem_a = FET.Elem1No;
      }
      else if (kind == KernelKind::SURFACE_FLUX ||
               (kind == KernelKind::INTERFACE_EPR &&
                epr_type == InterfaceDielectric::DEFAULT))
      {
        plan.elem_a = FET.Elem1No;
        if (has_elem2)
        {
          plan.elem_b = FET.Elem2No;
        }
      }
      else  // INTERFACE_EPR with MA, MS, or SA
      {
        // Single-sided evaluation on the vacuum (MA, SA) or substrate (MS) side, with
        // averaging if both sides qualify, skipping the element if neither does.
        const bool vacuum_side = (epr_type != InterfaceDielectric::MS);
        auto OnSide = [&](int elem_attr)
        {
          const double ls = mat_op->GetLightSpeedMax(elem_attr);
          return vacuum_side ? (ls >= threshold) : (ls < threshold);
        };
        const bool use_elem1 = OnSide(FET.Elem1->Attribute);
        const bool use_elem2 = has_elem2 && OnSide(FET.Elem2->Attribute);
        if (use_elem1)
        {
          plan.elem_a = FET.Elem1No;
          if (use_elem2)
          {
            plan.elem_b = FET.Elem2No;
          }
        }
        else if (use_elem2)
        {
          plan.elem_a = FET.Elem2No;
        }
        else
        {
          // Neither side qualifies: the element contributes zero, skip it.
          continue;
        }
      }

      // A two-sided interior boundary on a parallel interface: the Elem2 side is a face
      // neighbor (ghost) element on another process. Its field is not in the local
      // vector; it is pulled via FaceNbrFieldExchange (built below from the mapped
      // points) and fed to the ghost side of the operator, so the libCEED path handles
      // it without falling back to the legacy coefficients. Only Elem2 can be a ghost.
      if (plan.elem_a >= 0 && plan.elem_a >= pmesh.GetNE())
      {
        plan.ghost_a = true;
      }
      if (plan.elem_b >= 0 && !elem2_local)
      {
        plan.ghost_b = true;
      }
      if (plan.ghost_a || plan.ghost_b)
      {
        plan.face_nbr = FET.Elem2No - pmesh.GetNE();
        plan.ghost_attr = FET.Elem2->Attribute;
        // The ghost element index is out of range for pmesh.GetElementGeometry; capture
        // its geometry here while the face neighbor transformation is in scope.
        plan.ghost_geom = FET.Elem2->GetGeometryType();
      }

      // Map the face quadrature points to the volume element reference space(s):
      // boundary element reference coordinates -> face reference coordinates
      // (TransformBdrElementToFace with the boundary element to face orientation) ->
      // element 1/2 reference coordinates (FET.Loc1/Loc2).
      const auto bdr_geom = pmesh.GetBdrElementGeometry(i);
      const mfem::IntegrationRule &face_ir =
          mfem::IntRules.Get(bdr_geom, fem::DefaultIntegrationOrder::Get(pmesh, bdr_geom));
      const int nq = face_ir.GetNPoints();
      int f, o;
      pmesh.GetBdrElementFace(i, &f, &o);
      auto MapPoints = [&](mfem::IntegrationPointTransformation &loc,
                           std::vector<mfem::IntegrationPoint> &pts)
      {
        pts.resize(nq);
        for (int q = 0; q < nq; q++)
        {
          mfem::IntegrationPoint fip = mfem::Mesh::TransformBdrElementToFace(
              FET.GetGeometryType(), o, face_ir.IntPoint(q));
          loc.Transform(fip, pts[q]);
          pts[q].weight = face_ir.IntPoint(q).weight;
        }
      };
      if (plan.elem_a >= 0)
      {
        MapPoints(plan.elem_a == FET.Elem1No ? FET.Loc1 : FET.Loc2, plan.pts_a);
      }
      if (plan.elem_b >= 0)
      {
        MapPoints(FET.Loc2, plan.pts_b);
      }

      // Canonicalize the mapped volume cloud before deriving its finite registry key.
      // The boundary element's face orientation can only reorder this symmetric rule;
      // local faces and NC subfaces with genuinely different maps retain distinct clouds.
      const auto vol_geom_a = plan.ghost_a         ? plan.ghost_geom
                              : (plan.elem_a >= 0) ? pmesh.GetElementGeometry(plan.elem_a)
                                                   : mfem::Geometry::INVALID;
      const auto vol_geom_b = plan.ghost_b         ? plan.ghost_geom
                              : (plan.elem_b >= 0) ? pmesh.GetElementGeometry(plan.elem_b)
                                                   : mfem::Geometry::INVALID;
      if (plan.elem_a >= 0)
      {
        const auto canonical = CanonicalizeSymmetricMappedRule(plan.pts_a);
        plan.point_key_a =
            MakeCanonicalTracePointKey(vol_geom_a, bdr_geom, face_ir, canonical.pts);
      }
      if (plan.elem_b >= 0)
      {
        const auto canonical = CanonicalizeSymmetricMappedRule(plan.pts_b);
        plan.point_key_b =
            MakeCanonicalTracePointKey(vol_geom_b, bdr_geom, face_ir, canonical.pts);
      }

      const int out_slot = num_marked++;
      local_out_attrs.push_back(attr);

      auto AddGroup = [&](int elem_a, bool ghost_a, mfem::Geometry::Type geom_a,
                          const std::vector<mfem::IntegrationPoint> &pts_a,
                          const std::vector<long long> &point_key_a, int elem_b,
                          bool ghost_b, mfem::Geometry::Type geom_b,
                          const std::vector<mfem::IntegrationPoint> &pts_b,
                          const std::vector<long long> &point_key_b, double normal_scale)
      {
        // INTERFACE_EPR is the only non-separable paired operation. Keep its two sides
        // in their common legacy q order. Every one-sided group (including the split
        // buffer and flux paths) instead uses the canonical mapped cloud and retains a
        // canonical_q -> original_q output permutation.
        const bool canonicalize = elem_b < 0;
        const auto Canonicalize = [&](const std::vector<mfem::IntegrationPoint> &pts)
        {
          if (canonicalize)
          {
            return CanonicalizeSymmetricMappedRule(pts);
          }
          CanonicalMappedRule rule;
          rule.pts = pts;
          rule.canonical_to_original =
              IdentityPointPermutation(static_cast<int>(pts.size()));
          return rule;
        };
        const CanonicalMappedRule routed_a =
            elem_a >= 0 ? Canonicalize(pts_a) : CanonicalMappedRule{};
        const CanonicalMappedRule routed_b =
            elem_b >= 0 ? Canonicalize(pts_b) : CanonicalMappedRule{};

        FaceConfigKey key;
        key.reserve(12 + 4 * (routed_a.pts.size() + routed_b.pts.size()));
        key.push_back(static_cast<long long>(bdr_geom));
        key.push_back(static_cast<long long>(geom_a));
        key.push_back(static_cast<long long>(geom_b));
        key.push_back(static_cast<long long>(ghost_a));
        key.push_back(static_cast<long long>(ghost_b));
        key.push_back(static_cast<long long>(nq));
        // Only normal-dependent kernels need the boundary orientation/sign in their
        // context. Normal-independent averages otherwise fragment into redundant groups.
        const bool uses_normal =
            kind == KernelKind::SURFACE_FLUX || kind == KernelKind::FARFIELD;
        const double effective_normal_scale =
            uses_normal ? (plan.flip ? -normal_scale : normal_scale) : 1.0;
        key.push_back(EncodeDiscreteHalfScale(effective_normal_scale));
        if (kind == KernelKind::MODE_OVERLAP)
        {
          // Mode coefficients are boundary-attribute dependent.
          key.push_back(static_cast<long long>(attr));
        }
        const bool ghost_only_group = elem_a >= 0 && ghost_a && elem_b < 0;
        auto AppendMappedRule = [&key](const std::vector<mfem::IntegrationPoint> &pts)
        { AppendPointRuleSignature(key, pts); };
        // A one-sided ghost group consumes point values imported with EVAL_NONE. Its
        // local operation depends only on the rule size; individual canonical keys are
        // retained for the remote FaceNbrFieldExchange requests below. All locally
        // tabulated sides are keyed by their canonical ordered volume point cloud, not
        // by boundary-face orientation or local-face metadata.
        if (elem_a >= 0 && !ghost_only_group)
        {
          AppendMappedRule(routed_a.pts);
        }
        else
        {
          key.push_back(0);
        }
        if (elem_b >= 0)
        {
          AppendMappedRule(routed_b.pts);
        }
        else
        {
          key.push_back(0);
        }

        auto it = face_groups.find(key);
        if (it == face_groups.end())
        {
          FaceGroup group;
          group.bdr_geom = bdr_geom;
          group.vol_geom_a = geom_a;
          group.vol_geom_b = geom_b;
          group.normal_scale = effective_normal_scale;
          group.ghost_a = ghost_a;
          group.ghost_b = ghost_b;
          if (elem_a >= 0)
          {
            FaceConfigKey key_a = key;
            key_a.push_back(0);
            group.mapped_ir_a = GetRegisteredMappedIr(key_a, routed_a.pts);
          }
          if (elem_b >= 0)
          {
            FaceConfigKey key_b = key;
            key_b.push_back(1);
            group.mapped_ir_b = GetRegisteredMappedIr(key_b, routed_b.pts);
          }
          it = face_groups.emplace(key, std::move(group)).first;
        }
        else
        {
          if (elem_a >= 0 && !ghost_only_group)
          {
            MFEM_VERIFY(it->second.mapped_ir_a,
                        "Missing representative side-A mapped rule for trace group!");
            VerifyMappedRuleMatches(*it->second.mapped_ir_a, routed_a.pts);
          }
          if (elem_b >= 0)
          {
            MFEM_VERIFY(it->second.mapped_ir_b,
                        "Missing representative side-B mapped rule for trace group!");
            VerifyMappedRuleMatches(*it->second.mapped_ir_b, routed_b.pts);
          }
        }
        it->second.bdr_indices.push_back(i);
        if (elem_a >= 0)
        {
          if (ghost_a)
          {
            it->second.face_nbr.push_back(plan.face_nbr);
            it->second.ghost_attr.push_back(plan.ghost_attr);
            it->second.face_nbr_point_keys.push_back(point_key_a);
            it->second.face_nbr_points.push_back(routed_a.pts);
          }
          else
          {
            it->second.vol_indices_a.push_back(elem_a);
          }
        }
        if (elem_b >= 0)
        {
          if (ghost_b)
          {
            it->second.face_nbr.push_back(plan.face_nbr);
            it->second.ghost_attr.push_back(plan.ghost_attr);
            it->second.face_nbr_point_keys.push_back(point_key_b);
            it->second.face_nbr_points.push_back(routed_b.pts);
          }
          else
          {
            it->second.vol_indices_b.push_back(elem_b);
          }
        }
        it->second.out_slots.push_back(out_slot);
        it->second.canonical_to_original.push_back(routed_a.canonical_to_original);
      };

      const bool two_sided = plan.elem_a >= 0 && plan.elem_b >= 0;
      if (kind == KernelKind::SURFACE_FLUX && two_sided)
      {
        // Surface-flux sides are separable. The normal scale preserves the established
        // two-sided difference or one-sided-value average while both operators reduce
        // into the same boundary-element output slot.
        const double scale_a = flux_two_sided ? 1.0 : 0.5;
        const double scale_b = flux_two_sided ? -1.0 : 0.5;
        AddGroup(plan.elem_a, plan.ghost_a, vol_geom_a, plan.pts_a, plan.point_key_a, -1,
                 false, mfem::Geometry::INVALID, {}, {}, scale_a);
        AddGroup(plan.elem_b, plan.ghost_b, vol_geom_b, plan.pts_b, plan.point_key_b, -1,
                 false, mfem::Geometry::INVALID, {}, {}, scale_b);
      }
      else
      {
        // Keep non-separable two-sided interface EPR paired. All other unsplit groups
        // are ordinary one-sided fixed mapped-rule operators.
        AddGroup(plan.elem_a, plan.ghost_a, vol_geom_a, plan.pts_a, plan.point_key_a,
                 plan.elem_b, plan.ghost_b, vol_geom_b, plan.pts_b, plan.point_key_b, 1.0);
      }
    }
  }

  // Initialize the local output vector and field staging vector. Far-field operators
  // produce 6 values (Re/Im of a 3-vector) per direction per element.
  const int num_out =
      (kind == KernelKind::FARFIELD) ? 6 * static_cast<int>(farfield_dirs.size()) : 1;
  local_out.SetSize(num_marked * num_out);
  local_out.UseDevice(true);
  if (need_field)
  {
    const int max_vsize = std::max(nd_fespace ? nd_fespace->GetVSize() : 0,
                                   rt_fespace ? rt_fespace->GetVSize() : 0);
    field_staging.SetSize(max_vsize);
    field_staging.UseDevice(true);
    field_staging = 0.0;
  }

  // Build the (group independent part of the) QFunction context for the integrand.
  std::vector<CeedIntScalar> base_ctx = BuildBaseContext(dim, is_2d);
  // Build the face neighbor field exchange for any ghost (face neighbor) sides of
  // two-sided interior boundaries on parallel interfaces. Collective: all processes
  // participate (those without ghost faces pose no requests), so the decision is reduced
  // globally. Each ghost face requests the neighbor's volume field at the mapped face
  // quadrature points; the returned physical-space values feed the ghost side below.
  {
    int any_ghost = 0;
    for (const auto &[k, g] : face_groups)
    {
      if (g.ghost_a || g.ghost_b)
      {
        any_ghost = 1;
        break;
      }
    }
    Mpi::GlobalSum(1, &any_ghost, comm);
    if (any_ghost > 0)
    {
      // Source slots match the field inputs: SURFACE_FLUX uses slot 0 = E (nd_fespace)
      // and slot 1 = B (rt_fespace); all other kinds carry a single field at slot 0.
      std::array<const mfem::ParFiniteElementSpace *, FaceNbrFieldExchange::MaxSources>
          ex_fes = {nullptr, nullptr, nullptr, nullptr};
      unsigned int source_mask;
      if (kind == KernelKind::SURFACE_FLUX)
      {
        ex_fes[0] = nd_fespace;
        ex_fes[1] = rt_fespace;
        source_mask = (flux_type == SurfaceFlux::ELECTRIC)   ? 0b01u
                      : (flux_type == SurfaceFlux::MAGNETIC) ? 0b10u
                                                             : 0b11u;
      }
      else
      {
        ex_fes[0] = rt_fespace ? rt_fespace : nd_fespace;
        source_mask = 0b01u;
      }
      std::vector<FaceNbrFieldExchange::Request> requests;
      for (auto &[k, g] : face_groups)
      {
        if (!(g.ghost_a || g.ghost_b))
        {
          continue;
        }
        const mfem::IntegrationRule &ir = g.ghost_a ? *g.mapped_ir_a : *g.mapped_ir_b;
        std::vector<mfem::IntegrationPoint> pts(ir.GetNPoints());
        for (int q = 0; q < ir.GetNPoints(); q++)
        {
          pts[q] = ir.IntPoint(q);
        }
        MFEM_VERIFY(g.face_nbr_point_keys.size() == g.face_nbr.size() &&
                        g.face_nbr_points.size() == g.face_nbr.size(),
                    "Invalid face-neighbor point-key layout for surface functional!");
        g.req_idx.resize(g.face_nbr.size());
        for (std::size_t e = 0; e < g.face_nbr.size(); e++)
        {
          g.req_idx[e] = static_cast<int>(requests.size());
          auto &req = requests.emplace_back();
          req.face_nbr_elem = g.face_nbr[e];
          req.source_mask = source_mask;
          req.point_key = g.face_nbr_point_keys[e];
          req.pts = g.face_nbr_points[e].empty() ? pts : g.face_nbr_points[e];
        }
      }
      face_nbr_exchange = std::make_unique<FaceNbrFieldExchange>(mesh, ex_fes, requests);
    }
  }

  // Assemble a libCEED operator for each group. For now, all operators are constructed
  // on a single Ceed context (no OpenMP parallel assembly or application; correctness
  // first, this can be extended with the thread partitioning of fem/mesh.cpp later).
  for (auto &face_group : face_groups)
  {
    auto &group = face_group.second;
    const std::size_t num_elem = group.bdr_indices.size();
    const bool has_b = !group.vol_indices_b.empty() || group.ghost_b;
    const mfem::IntegrationRule &face_ir = mfem::IntRules.Get(
        group.bdr_geom, fem::DefaultIntegrationOrder::Get(pmesh, group.bdr_geom));

    // Objects are owned by the assembled operator; scratch destroys our references.
    CeedAssemblyScratch scratch(ceed);

    // Assemble the inputs in the order expected by the QFunctions: quadrature weights
    // and boundary element mesh node gradients (surface measure and normal), per-side
    // volume element attributes (material lookup) and mesh node gradients at the mapped
    // points (Piola transformations), coordinates (SURFACE_FLUX only), then the field
    // inputs (side 1 fields, then side 2 fields). The geometry is computed on the fly
    // in the QFunctions (no stored geometry factor data).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    std::vector<std::string> mesh_node_fields;
    const bool needs_face_geom = true;
    const bool needs_elem_attr = kind != KernelKind::MODE_OVERLAP;

    // Fixed mapped volume rules use canonical q order. Face-only quantities must use
    // that same order, but a boundary-element basis is naturally tabulated in the
    // legacy boundary-rule order. Precompute the small face geometry/coordinate arrays
    // once at setup and pass them through ordinary EVAL_NONE restrictions. This is a
    // setup-only CPU loop; all recurring field evaluation and routing remains device
    // libCEED work.
    const int face_nq = face_ir.GetNPoints();
    MFEM_VERIFY(group.canonical_to_original.size() == num_elem,
                "Missing canonical permutation for face-only point data!");
    auto AddFixedFacePointInput = [&](const std::string &name, Vector &data, int num_comp)
    {
      std::vector<CeedInt> offsets(num_elem * face_nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int q = 0; q < face_nq; q++)
        {
          offsets[e * face_nq + q] = static_cast<CeedInt>(e * num_comp * face_nq + q);
        }
      }
      CeedElemRestriction restr;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem), face_nq,
                                          num_comp, face_nq, data.Size(), CEED_MEM_HOST,
                                          CEED_COPY_VALUES, offsets.data(), &restr));
      CeedVector vec;
      ceed::InitCeedVector(data, ceed, &vec);
      inputs.push_back({name, vec, restr, nullptr, ceed::EvalMode::None});
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
    };
    auto OriginalFaceQ = [&](std::size_t e, int canonical_q)
    {
      const auto &perm = group.canonical_to_original[e];
      MFEM_VERIFY(perm.size() == static_cast<std::size_t>(face_nq),
                  "Invalid canonical permutation for face-only point data!");
      const int original_q = perm[static_cast<std::size_t>(canonical_q)];
      MFEM_VERIFY(original_q >= 0 && original_q < face_nq,
                  "Invalid canonical face-only point permutation entry!");
      return original_q;
    };
    if (needs_face_geom && kind == KernelKind::AREA)
    {
      // AREA has no volume basis from which the reduction helper can obtain Q. It is
      // order-invariant, so retain its ordinary fixed face-basis pathway.
      const mfem::FiniteElement *face_mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(group.bdr_geom);
      if (!face_mesh_fe)
      {
        face_mesh_fe = mesh_fespace.FEColl()->TraceFiniteElementForGeometry(group.bdr_geom);
      }
      CeedBasis face_mesh_basis;
      ceed::InitCachedBasisFromRule(*face_mesh_fe, face_ir, mesh_fespace.GetVDim(), ceed,
                                    &face_mesh_basis);
      CeedElemRestriction face_mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, group.bdr_geom, group.bdr_indices, /*is_interp*/ true);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back({"qw", nullptr, nullptr, face_mesh_basis, ceed::EvalMode::Weight});
      inputs.push_back(
          {"x_f", mesh_nodes_vec, face_mesh_restr, face_mesh_basis, ceed::EvalMode::Grad});
      mesh_node_fields.push_back("grad_x_f");
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(face_mesh_restr);
      scratch.bases.push_back(face_mesh_basis);
    }
    else if (needs_face_geom)
    {
      auto &qw = elem_attrs.emplace_back(num_elem * face_nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int q = 0; q < face_nq; q++)
        {
          qw[e * face_nq + q] = face_ir.IntPoint(OriginalFaceQ(e, q)).weight;
        }
      }
      AddFixedFacePointInput("qw", qw, 1);
      const int face_jac_size = sdim * (dim - 1);
      auto &face_geom = elem_attrs.emplace_back(num_elem * face_jac_size * face_nq);
      face_geom = 0.0;
      for (std::size_t e = 0; e < num_elem; e++)
      {
        mfem::IsoparametricTransformation T;
        pmesh.GetBdrElementTransformation(group.bdr_indices[e], &T);
        for (int q = 0; q < face_nq; q++)
        {
          const auto &ip = face_ir.IntPoint(OriginalFaceQ(e, q));
          T.SetIntPoint(&ip);
          const mfem::DenseMatrix &J = T.Jacobian();
          for (int d = 0; d < J.Width(); d++)
          {
            for (int c = 0; c < J.Height(); c++)
            {
              face_geom[e * face_jac_size * face_nq + (c + sdim * d) * face_nq + q] =
                  J(c, d);
            }
          }
        }
      }
      AddFixedFacePointInput("grad_x_f", face_geom, face_jac_size);
    }
    auto AddVolGeomInputs = [&](const std::string &suffix, const std::vector<int> &indices,
                                mfem::Geometry::Type geom, const mfem::IntegrationRule &ir)
    {
      const int num_pts = ir.GetNPoints();
      CeedElemRestriction attr_restr = nullptr;
      CeedBasis attr_basis = nullptr;
      CeedVector attr_vec = nullptr;
      if (needs_elem_attr)
      {
        const auto &loc_attr = mesh.GetCeedAttributes();
        auto &elem_attr = elem_attrs.emplace_back(indices.size());
        for (std::size_t k = 0; k < indices.size(); k++)
        {
          elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
        }
        PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                                 ceed, indices.size(), 1, 1, indices.size(),
                                 CEED_STRIDES_BACKEND, &attr_restr));
        {
          // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
          mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
          Bt = 1.0;
          Gt = 0.0;
          qX = 0.0;
          qW = 0.0;
          PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                                 Bt.GetData(), Gt.GetData(), qX.GetData(),
                                                 qW.GetData(), &attr_basis));
        }
        ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
        inputs.push_back(
            {"attr_" + suffix, attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      }
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitCachedBasisFromRule(*mesh_fe, ir, mesh_fespace.GetVDim(), ceed,
                                    &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back(
          {"x_" + suffix, mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      if (needs_elem_attr)
      {
        scratch.vecs.push_back(attr_vec);
        scratch.restrs.push_back(attr_restr);
        scratch.bases.push_back(attr_basis);
      }
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(mesh_basis);
    };
    // Ghost (face neighbor) side variant: the volume element lives on another process,
    // so its attributes come from group.ghost_attr and the Piola Jacobian is replaced by
    // a constant identity. The exchanged field values are already physical (the owning
    // process applied the Piola map), so Piola(I) = I in the shared kernels reproduces
    // them. Same input names/positions as AddVolGeomInputs ("grad_x_<suffix>" via
    // EVAL_NONE matches the EVAL_GRAD naming), so the kernels are unchanged.
    auto AddVolGeomInputsGhost =
        [&](const std::string &suffix, const mfem::IntegrationRule &ir)
    {
      const int num_pts = ir.GetNPoints();
      CeedElemRestriction attr_restr = nullptr;
      CeedBasis attr_basis = nullptr;
      CeedVector attr_vec = nullptr;
      if (needs_elem_attr)
      {
        auto &elem_attr = elem_attrs.emplace_back(group.ghost_attr.size());
        const auto &loc_attr = mesh.GetCeedAttributes();
        for (std::size_t k = 0; k < group.ghost_attr.size(); k++)
        {
          elem_attr[k] = loc_attr.at(group.ghost_attr[k]);
        }
        PalaceCeedCall(ceed,
                       CeedElemRestrictionCreateStrided(ceed, group.ghost_attr.size(), 1, 1,
                                                        group.ghost_attr.size(),
                                                        CEED_STRIDES_BACKEND, &attr_restr));
        {
          mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
          Bt = 1.0;
          Gt = 0.0;
          qX = 0.0;
          qW = 0.0;
          PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                                 Bt.GetData(), Gt.GetData(), qX.GetData(),
                                                 qW.GetData(), &attr_basis));
        }
        ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
        inputs.push_back(
            {"attr_" + suffix, attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      }

      // Constant identity Jacobian (component-major [elem][comp][pt]), passed directly
      // to the kernel's grad_x_<suffix> input (EVAL_NONE). It is 2x2 for 2D line
      // integrals and 3x3 for 3D surface integrals.
      const int geom_comp = dim * sdim;
      auto &ident = elem_attrs.emplace_back(num_elem * geom_comp * num_pts);
      ident = 0.0;
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int d = 0; d < dim; d++)
        {
          const int c = d * (sdim + 1);
          for (int i = 0; i < num_pts; i++)
          {
            ident[e * geom_comp * num_pts + c * num_pts + i] = 1.0;
          }
        }
      }
      CeedElemRestriction ident_restr;
      const CeedInt strides[3] = {1, num_pts, geom_comp * num_pts};
      PalaceCeedCall(ceed,
                     CeedElemRestrictionCreateStrided(
                         ceed, static_cast<CeedInt>(num_elem), num_pts, geom_comp,
                         (CeedSize)num_elem * geom_comp * num_pts, strides, &ident_restr));
      CeedVector ident_vec;
      ceed::InitCeedVector(ident, ceed, &ident_vec);
      inputs.push_back(
          {"grad_x_" + suffix, ident_vec, ident_restr, nullptr, ceed::EvalMode::None});
      if (needs_elem_attr)
      {
        scratch.vecs.push_back(attr_vec);
        scratch.restrs.push_back(attr_restr);
        scratch.bases.push_back(attr_basis);
      }
      scratch.vecs.push_back(ident_vec);
      scratch.restrs.push_back(ident_restr);
    };
    if (group.ghost_a)
    {
      AddVolGeomInputsGhost("1", *group.mapped_ir_a);
    }
    else if (!group.vol_indices_a.empty())
    {
      AddVolGeomInputs("1", group.vol_indices_a, group.vol_geom_a, *group.mapped_ir_a);
    }
    if (group.ghost_b)
    {
      AddVolGeomInputsGhost("2", *group.mapped_ir_b);
    }
    else if (!group.vol_indices_b.empty())
    {
      AddVolGeomInputs("2", group.vol_indices_b, group.vol_geom_b, *group.mapped_ir_b);
    }
    if (kind == KernelKind::SURFACE_FLUX || kind == KernelKind::FARFIELD ||
        kind == KernelKind::MODE_OVERLAP)
    {
      // Like grad_x_f above, coordinates must follow the canonical volume q order.
      // Store component-major per-element point data for the ordinary EVAL_NONE
      // restriction, rather than reintroducing a runtime mapped-point basis.
      auto &x_pts = elem_attrs.emplace_back(num_elem * sdim * face_nq);
      x_pts = 0.0;
      for (std::size_t e = 0; e < num_elem; e++)
      {
        mfem::IsoparametricTransformation T;
        pmesh.GetBdrElementTransformation(group.bdr_indices[e], &T);
        for (int q = 0; q < face_nq; q++)
        {
          const auto &ip = face_ir.IntPoint(OriginalFaceQ(e, q));
          mfem::Vector x(sdim);
          T.Transform(ip, x);
          for (int d = 0; d < sdim; d++)
          {
            x_pts[e * sdim * face_nq + d * face_nq + q] = x(d);
          }
        }
      }
      AddFixedFacePointInput("x", x_pts, sdim);
    }
    auto AddFieldInput = [&](const std::string &name, int source,
                             const mfem::ParFiniteElementSpace &fespace,
                             const std::vector<int> &indices, mfem::Geometry::Type geom,
                             const mfem::IntegrationRule &ir)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "Unable to get field finite element for surface functional!");
      ceed::InitCachedBasisFromRule(*fe, ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    // Ghost (face neighbor) field input: the neighbor's physical field values at the
    // mapped points, imported via FaceNbrFieldExchange. The indexed restriction slices
    // the shared imported vector per element (base = ImportOffset) and transposes its
    // point-major layout ([pt][comp]) to the component-major layout the kernel reads
    // (offset base + sdim*pt, comp_stride 1). Re-pointed at the imported vector on each
    // apply via source index 4 (see fem::ApplyAddGroupOperators).
    auto AddFieldInputGhost =
        [&](const std::string &name, int slot, const mfem::IntegrationRule &ir)
    {
      const int nq = ir.GetNPoints();
      std::vector<CeedInt> offsets(num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        const int base = face_nbr_exchange->ImportOffset(group.req_idx[e], slot);
        for (int i = 0; i < nq; i++)
        {
          offsets[e * nq + i] = base + sdim * i;
        }
      }
      CeedElemRestriction restr;
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(num_elem), nq, sdim, 1,
                               (CeedSize)face_nbr_exchange->ImportSize(), CEED_MEM_HOST,
                               CEED_COPY_VALUES, offsets.data(), &restr));
      CeedVector vec;
      ceed::InitCeedVector(face_nbr_exchange->Imported(), ceed, &vec);
      inputs.push_back({name, vec, restr, nullptr, ceed::EvalMode::None});
      field_sources.emplace_back(name, 4);  // 4 -> imported (see ApplyAddGroupOperators)
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
    };
    if (kind == KernelKind::HCURL_NORM2 || kind == KernelKind::INTERFACE_EPR ||
        kind == KernelKind::MODE_OVERLAP)
    {
      const auto &field_fespace = rt_fespace ? *rt_fespace : *nd_fespace;
      if (group.ghost_a)
      {
        AddFieldInputGhost("u_1", 0, *group.mapped_ir_a);
      }
      else
      {
        AddFieldInput("u_1", 0, field_fespace, group.vol_indices_a, group.vol_geom_a,
                      *group.mapped_ir_a);
      }
      if (has_b)
      {
        if (group.ghost_b)
        {
          AddFieldInputGhost("u_2", 0, *group.mapped_ir_b);
        }
        else
        {
          AddFieldInput("u_2", 0, field_fespace, group.vol_indices_b, group.vol_geom_b,
                        *group.mapped_ir_b);
        }
      }
    }
    else if (kind == KernelKind::FARFIELD)
    {
      const mfem::IntegrationRule &ir = *group.mapped_ir_a;
      AddFieldInput("u_1", 0, *nd_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_2", 1, *nd_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_3", 2, *rt_fespace, group.vol_indices_a, group.vol_geom_a, ir);
      AddFieldInput("u_4", 3, *rt_fespace, group.vol_indices_a, group.vol_geom_a, ir);
    }
    else if (kind == KernelKind::SURFACE_FLUX)
    {
      int count = 0;
      auto AddSide = [&](const std::vector<int> &indices, mfem::Geometry::Type geom,
                         const mfem::IntegrationRule &ir, bool ghost)
      {
        if (flux_type == SurfaceFlux::ELECTRIC || flux_type == SurfaceFlux::POWER)
        {
          const std::string nm = "u_" + std::to_string(++count);
          if (ghost)
          {
            AddFieldInputGhost(nm, 0, ir);
          }
          else
          {
            AddFieldInput(nm, 0, *nd_fespace, indices, geom, ir);
          }
        }
        if (flux_type == SurfaceFlux::MAGNETIC || flux_type == SurfaceFlux::POWER)
        {
          const std::string nm = "u_" + std::to_string(++count);
          if (ghost)
          {
            AddFieldInputGhost(nm, 1, ir);
          }
          else
          {
            AddFieldInput(nm, 1, *rt_fespace, indices, geom, ir);
          }
        }
      };
      AddSide(group.vol_indices_a, group.vol_geom_a, *group.mapped_ir_a, group.ghost_a);
      if (has_b)
      {
        AddSide(group.vol_indices_b, group.vol_geom_b, *group.mapped_ir_b, group.ghost_b);
      }
    }

    // Output restriction: num_out slots per boundary element in the local output vector
    // (component stride num_marked).
    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, static_cast<CeedInt>(num_elem), 1, num_out, num_marked,
                             (CeedSize)num_marked * num_out, CEED_MEM_HOST,
                             CEED_COPY_VALUES, group.out_slots.data(), &out_restr));
    scratch.restrs.push_back(out_restr);

    // Select the QFunction and finalize the (group dependent) context.
    int bdr_attr = 0;
    if (kind == KernelKind::MODE_OVERLAP)
    {
      MFEM_VERIFY(!group.bdr_indices.empty(),
                  "Empty boundary group for mode-overlap functional!");
      bdr_attr = pmesh.GetBdrAttribute(group.bdr_indices.front());
    }
    std::vector<CeedIntScalar> ctx;
    ceed::CeedQFunctionInfo info;
    ConfigureGroupQFunction(is_2d, has_b, group.normal_scale, bdr_attr, base_ctx, ctx,
                            info);

    // Assemble the operator.
    CeedOperator op;
    CeedQFunctionContext op_ctx = nullptr;
    ceed::AssembleCeedSurfaceFunctional(
        info, ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, inputs, num_out,
        out_restr, &op, (kind == KernelKind::FARFIELD) ? &op_ctx : nullptr);
    groups.push_back({ceed, op, std::move(field_sources), op_ctx});
    if (!mesh_node_fields.empty())
    {
      groups.back().mesh_nodes = pmesh.GetNodes();
      groups.back().mesh_node_fields = std::move(mesh_node_fields);
    }
    fem::CacheGroupOperatorFieldVectors(groups.back());
  }
}

void SurfaceFunctional::ApplyAdd(const std::array<const Vector *, 4> &srcs) const
{
  if (face_nbr_exchange)
  {
    face_nbr_exchange->Exchange(srcs);
    fem::ApplyAddGroupOperators(groups, srcs, local_out, &face_nbr_exchange->Imported());
  }
  else
  {
    fem::ApplyAddGroupOperators(groups, srcs, local_out);
  }
}

double SurfaceFunctional::EvalLocal(const std::array<const Vector *, 4> &srcs) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) SurfaceFunctional!");
  // A process holding the face neighbor exchange must still apply (the exchange is
  // collective: it may need to export its local field for a neighbor's ghost request
  // even with no local marked elements). ApplyAdd is a no-op on the empty group set.
  if (local_out.Size() == 0 && !face_nbr_exchange)
  {
    return 0.0;
  }
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
  }
  ApplyAdd(srcs);
  return (local_out.Size() > 0) ? linalg::LocalSum(local_out) : 0.0;
}

void SurfaceFunctional::BinLocalOutByAttribute(const mfem::Array<int> &attr_to_bin,
                                               int num_bins, std::vector<double> &bins,
                                               double scale) const
{
  if (local_out.Size() == 0)
  {
    return;
  }
  MFEM_VERIFY(local_out_attrs.size() == static_cast<std::size_t>(local_out.Size()),
              "SurfaceFunctional attribute bins require one output slot per element!");
  const double *vals = local_out.HostRead();
  for (int i = 0; i < local_out.Size(); i++)
  {
    const int attr = local_out_attrs[i];
    const int bin = (attr > 0 && attr <= attr_to_bin.Size()) ? attr_to_bin[attr - 1] : -1;
    if (bin >= 0)
    {
      MFEM_VERIFY(bin < num_bins, "SurfaceFunctional attribute bin out of range!");
      bins[bin] += scale * vals[i];
    }
  }
}

double SurfaceFunctional::Eval(const Vector *u) const
{
  MFEM_VERIFY(kind == KernelKind::AREA || u,
              "SurfaceFunctional::Eval requires a field vector for functionals with "
              "field inputs!");
  double dot = EvalLocal({u, nullptr});
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

double SurfaceFunctional::Eval(const GridFunction &u) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) SurfaceFunctional!");
  MFEM_VERIFY(kind == KernelKind::HCURL_NORM2 || kind == KernelKind::INTERFACE_EPR,
              "SurfaceFunctional::Eval with a grid function is only valid for functionals "
              "quadratic in a single field!");
  double dot = 0.0;
  if (local_out.Size() > 0 || face_nbr_exchange)
  {
    if (local_out.Size() > 0)
    {
      local_out = 0.0;
    }
    ApplyAdd({&u.Real(), nullptr});
    if (u.HasImag())
    {
      ApplyAdd({&u.Imag(), nullptr});
    }
    if (local_out.Size() > 0)
    {
      dot = linalg::LocalSum(local_out);
    }
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

std::complex<double> SurfaceFunctional::EvalModeOverlap(const GridFunction &E) const
{
  MFEM_VERIFY(kind == KernelKind::MODE_OVERLAP,
              "SurfaceFunctional::EvalModeOverlap is only valid for mode-overlap "
              "functionals!");
  std::complex<double> dot(EvalLocal({&E.Real(), nullptr}), 0.0);
  if (E.HasImag())
  {
    dot.imag(EvalLocal({&E.Imag(), nullptr}));
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

std::vector<std::complex<double>> SurfaceFunctional::EvalModeOverlapByAttribute(
    const GridFunction &E, const mfem::Array<int> &attr_to_bin, int num_bins) const
{
  MFEM_VERIFY(kind == KernelKind::MODE_OVERLAP,
              "SurfaceFunctional::EvalModeOverlapByAttribute is only valid for "
              "mode-overlap functionals!");
  MFEM_VERIFY(num_bins >= 0, "Invalid number of output bins!");
  auto AccumulateBins = [&](const Vector &u, std::vector<double> &bins)
  {
    if (local_out.Size() > 0)
    {
      local_out = 0.0;
    }
    // Keep the apply collective even on ranks with no local marked elements, matching
    // the per-port EvalModeOverlap path and the binned power evaluator.
    ApplyAdd({&u, nullptr});
    BinLocalOutByAttribute(attr_to_bin, num_bins, bins, 1.0);
  };

  std::vector<double> real(num_bins, 0.0), imag(num_bins, 0.0);
  AccumulateBins(E.Real(), real);
  if (E.HasImag())
  {
    AccumulateBins(E.Imag(), imag);
  }

  std::vector<double> packed(2 * num_bins, 0.0);
  for (int i = 0; i < num_bins; i++)
  {
    packed[2 * i + 0] = real[i];
    packed[2 * i + 1] = imag[i];
  }
  Mpi::GlobalSum(static_cast<int>(packed.size()), packed.data(), comm);

  std::vector<std::complex<double>> result(num_bins);
  for (int i = 0; i < num_bins; i++)
  {
    result[i] = {packed[2 * i + 0], packed[2 * i + 1]};
  }
  return result;
}

std::complex<double> SurfaceFunctional::EvalFlux(const GridFunction *E,
                                                 const GridFunction *B) const
{
  MFEM_VERIFY(kind == KernelKind::SURFACE_FLUX,
              "SurfaceFunctional::EvalFlux is only valid for surface flux functionals!");
  MFEM_VERIFY(
      (E || (flux_type != SurfaceFlux::ELECTRIC && flux_type != SurfaceFlux::POWER)) &&
          (B || (flux_type != SurfaceFlux::MAGNETIC && flux_type != SurfaceFlux::POWER)),
      "Missing E or B field grid function for surface flux evaluation!");

  // For complex-valued fields, output the separate real and imaginary parts for the
  // time-harmonic quantity. For power flux (Poynting vector), output only the
  // stationary real part and not the part which has double the frequency (the real and
  // imaginary part contributions add).
  const bool has_imag = E ? E->HasImag() : B->HasImag();
  std::complex<double> dot(EvalLocal({E ? &E->Real() : nullptr, B ? &B->Real() : nullptr}),
                           0.0);
  if (has_imag)
  {
    const double doti = EvalLocal({E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr});
    if (flux_type == SurfaceFlux::POWER)
    {
      dot += doti;
    }
    else
    {
      dot.imag(doti);
    }
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

std::vector<std::array<std::complex<double>, 3>>
SurfaceFunctional::EvalFarField(const GridFunction &E, const GridFunction &B,
                                std::complex<double> omega)
{
  MFEM_VERIFY(kind == KernelKind::FARFIELD && E.HasImag() && B.HasImag(),
              "SurfaceFunctional::EvalFarField requires a far-field functional and "
              "complex-valued fields!");
  MFEM_VERIFY(valid, "EvalFarField called on an invalid (unassembled) SurfaceFunctional!");

  // The frequency enters only the QFunction context (the omega slots); update it in
  // place rather than reassembling the operators. Reassembly would rebuild the bases,
  // restrictions, and on-the-fly geometry inputs and re-JIT the (expensive) far-field
  // kernel on every frequency -- none of which depend on omega. FARFIELD context layout
  // (see Assemble): [0] normal sign, [1] omega_re, [2] omega_im, [3] N,
  // [4] B Piola map, [5..] directions, then the material context.
  if (omega != farfield_omega)
  {
    farfield_omega = omega;
    for (auto &group : groups)
    {
      if (!group.ctx)
      {
        continue;
      }
      CeedIntScalar *data;
      PalaceCeedCall(group.ceed,
                     CeedQFunctionContextGetData(group.ctx, CEED_MEM_HOST, &data));
      data[1].second = omega.real();
      data[2].second = omega.imag();
      PalaceCeedCall(group.ceed, CeedQFunctionContextRestoreData(group.ctx, &data));
    }
  }

  // Integrate, reduce each component over the local elements and all processes, and
  // apply the final cross products (following GetFarFieldrE).
  const int N = static_cast<int>(farfield_dirs.size());
  const int num_marked = local_out.Size() / std::max(6 * N, 1);
  std::vector<double> integrals(6 * N, 0.0);
  if (local_out.Size() > 0)
  {
    local_out = 0.0;
    fem::ApplyAddGroupOperators(groups, {&E.Real(), &E.Imag(), &B.Real(), &B.Imag()},
                                local_out);
    Vector slice;
    slice.UseDevice(true);
    for (int c = 0; c < 6 * N; c++)
    {
      slice.MakeRef(local_out, c * num_marked, num_marked);
      integrals[c] = linalg::LocalSum(slice);
    }
  }
  Mpi::GlobalSum(6 * N, integrals.data(), comm);

  std::vector<std::array<std::complex<double>, 3>> result(N);
  for (int d = 0; d < N; d++)
  {
    const auto &r = farfield_dirs[d];
    const double *Ir = integrals.data() + 6 * d, *Ii = integrals.data() + 6 * d + 3;
    const double cr[3] = {r[1] * Ir[2] - r[2] * Ir[1], r[2] * Ir[0] - r[0] * Ir[2],
                          r[0] * Ir[1] - r[1] * Ir[0]};
    const double ci[3] = {r[1] * Ii[2] - r[2] * Ii[1], r[2] * Ii[0] - r[0] * Ii[2],
                          r[0] * Ii[1] - r[1] * Ii[0]};
    for (int c = 0; c < 3; c++)
    {
      result[d][c] = {cr[c], ci[c]};
    }
  }
  return result;
}

std::complex<double> SurfaceFunctional::EvalComplexPower(const GridFunction &E,
                                                         const GridFunction &B) const
{
  MFEM_VERIFY(kind == KernelKind::SURFACE_FLUX && flux_type == SurfaceFlux::POWER &&
                  flux_two_sided,
              "SurfaceFunctional::EvalComplexPower is only valid for two-sided POWER "
              "flux functionals!");
  MFEM_VERIFY(E.HasImag() == B.HasImag(),
              "Mismatch between real- and complex-valued E and B fields in port power "
              "calculation!");

  // Following LumpedPortData::GetPower: P = ∫ (E x H*) ⋅ n dS =
  // -∫ E ⋅ (n x H*) dS, with H = μ⁻¹ B and n the normal oriented into element 1
  // (contributions from both sides of an interior boundary add). With
  // S(e, b) = ∫ (e x μ⁻¹ b) ⋅ n dS (the two-sided POWER flux functional),
  //   Re{P} = S(E_re, B_re) + S(E_im, B_im)
  //   Im{P} = S(E_im, B_re) - S(E_re, B_im) .
  const bool has_imag = E.HasImag();
  std::complex<double> dot(EvalLocal({&E.Real(), &B.Real()}), 0.0);
  if (has_imag)
  {
    dot += EvalLocal({&E.Imag(), &B.Imag()});
    dot.imag(EvalLocal({&E.Imag(), &B.Real()}) - EvalLocal({&E.Real(), &B.Imag()}));
  }
  Mpi::GlobalSum(1, &dot, comm);
  return dot;
}

std::vector<std::complex<double>>
SurfaceFunctional::EvalComplexPowerByAttribute(const GridFunction &E, const GridFunction &B,
                                               const mfem::Array<int> &attr_to_bin,
                                               int num_bins) const
{
  MFEM_VERIFY(kind == KernelKind::SURFACE_FLUX && flux_type == SurfaceFlux::POWER &&
                  flux_two_sided,
              "SurfaceFunctional::EvalComplexPowerByAttribute is only valid for "
              "two-sided POWER flux functionals!");
  MFEM_VERIFY(E.HasImag() == B.HasImag(),
              "Mismatch between real- and complex-valued E and B fields in batched "
              "port power calculation!");
  MFEM_VERIFY(num_bins >= 0, "Invalid number of output bins!");

  auto AccumulateBins = [&](const std::array<const Vector *, 4> &srcs,
                            std::vector<double> &bins, double scale)
  {
    if (local_out.Size() > 0)
    {
      local_out = 0.0;
    }
    // Keep the apply collective even on ranks with no local marked elements, since a
    // face-neighbor exchange may need this rank to export field data for another rank's
    // processor-boundary side.
    ApplyAdd(srcs);
    BinLocalOutByAttribute(attr_to_bin, num_bins, bins, scale);
  };

  std::vector<double> real(num_bins, 0.0), imag(num_bins, 0.0);
  AccumulateBins({&E.Real(), &B.Real()}, real, 1.0);
  if (E.HasImag())
  {
    AccumulateBins({&E.Imag(), &B.Imag()}, real, 1.0);
    AccumulateBins({&E.Imag(), &B.Real()}, imag, 1.0);
    AccumulateBins({&E.Real(), &B.Imag()}, imag, -1.0);
  }

  std::vector<double> packed(2 * num_bins, 0.0);
  for (int i = 0; i < num_bins; i++)
  {
    packed[2 * i + 0] = real[i];
    packed[2 * i + 1] = imag[i];
  }
  Mpi::GlobalSum(static_cast<int>(packed.size()), packed.data(), comm);

  std::vector<std::complex<double>> result(num_bins);
  for (int i = 0; i < num_bins; i++)
  {
    result[i] = {packed[2 * i + 0], packed[2 * i + 1]};
  }
  return result;
}

}  // namespace palace
