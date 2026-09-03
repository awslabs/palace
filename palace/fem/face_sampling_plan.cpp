// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/face_sampling_plan.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <numeric>
#include <tuple>
#include "fem/coefficient.hpp"
#include "fem/mesh.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

long long EncodePointRuleDouble(double x)
{
  static_assert(sizeof(long long) == sizeof(double));
  long long bits;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

using CoordinateKey = std::array<long long, 3>;

void NormalizeReferencePoint(mfem::IntegrationPoint &ip);

CoordinateKey MakeCoordinateKey(mfem::IntegrationPoint ip)
{
  // Reference-coordinate signed zero is not a distinct evaluation point. We retain
  // every other IEEE-754 coordinate bit exactly and intentionally do not inspect the
  // quadrature weight: it is irrelevant to a physical point-value evaluation.
  NormalizeReferencePoint(ip);
  return {EncodePointRuleDouble(ip.x), EncodePointRuleDouble(ip.y),
          EncodePointRuleDouble(ip.z)};
}

void AppendPointRuleSignature(std::vector<long long> &key,
                              const std::vector<mfem::IntegrationPoint> &pts)
{
  key.push_back(static_cast<long long>(pts.size()));
  for (const auto &ip : pts)
  {
    key.push_back(EncodePointRuleDouble(ip.x));
    key.push_back(EncodePointRuleDouble(ip.y));
    key.push_back(EncodePointRuleDouble(ip.z));
    key.push_back(EncodePointRuleDouble(ip.weight));
  }
}

void NormalizeReferencePoint(mfem::IntegrationPoint &ip)
{
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

int ReferenceFacet(mfem::Geometry::Type geometry,
                   const std::vector<mfem::IntegrationPoint> &points)
{
  constexpr double tol = 1.0e-12;
  auto OnFacet = [&](auto predicate)
  {
    return !points.empty() && std::all_of(points.begin(), points.end(),
                                          [&](const mfem::IntegrationPoint &point)
                                          { return predicate(point); });
  };
  switch (geometry)
  {
    case mfem::Geometry::TRIANGLE:
      if (OnFacet([&](const auto &point) { return std::abs(point.x) <= tol; }))
        return 0;
      if (OnFacet([&](const auto &point) { return std::abs(point.y) <= tol; }))
        return 1;
      if (OnFacet([&](const auto &point)
                  { return std::abs(point.x + point.y - 1.0) <= tol; }))
        return 2;
      break;
    case mfem::Geometry::SQUARE:
      for (int d = 0; d < 2; d++)
      {
        if (OnFacet([&](const auto &point)
                    { return std::abs((d == 0 ? point.x : point.y)) <= tol; }))
          return 2 * d;
        if (OnFacet([&](const auto &point)
                    { return std::abs((d == 0 ? point.x : point.y) - 1.0) <= tol; }))
          return 2 * d + 1;
      }
      break;
    case mfem::Geometry::TETRAHEDRON:
      for (int d = 0; d < 3; d++)
      {
        if (OnFacet(
                [&](const auto &point)
                {
                  return std::abs(d == 0 ? point.x : d == 1 ? point.y : point.z) <= tol;
                }))
          return d;
      }
      if (OnFacet([&](const auto &point)
                  { return std::abs(point.x + point.y + point.z - 1.0) <= tol; }))
        return 3;
      break;
    case mfem::Geometry::CUBE:
      for (int d = 0; d < 3; d++)
      {
        if (OnFacet(
                [&](const auto &point)
                {
                  return std::abs(d == 0 ? point.x : d == 1 ? point.y : point.z) <= tol;
                }))
          return 2 * d;
        if (OnFacet(
                [&](const auto &point)
                {
                  return std::abs((d == 0   ? point.x
                                   : d == 1 ? point.y
                                            : point.z) -
                                  1.0) <= tol;
                }))
          return 2 * d + 1;
      }
      break;
    case mfem::Geometry::PRISM:
      // Reference wedge = reference triangle x [0,1]. The returned numbers need only
      // distinguish facets of one owner; they are not MFEM local-face identifiers.
      if (OnFacet([&](const auto &point) { return std::abs(point.z) <= tol; }))
        return 0;
      if (OnFacet([&](const auto &point) { return std::abs(point.z - 1.0) <= tol; }))
        return 1;
      if (OnFacet([&](const auto &point) { return std::abs(point.x) <= tol; }))
        return 2;
      if (OnFacet([&](const auto &point) { return std::abs(point.y) <= tol; }))
        return 3;
      if (OnFacet([&](const auto &point)
                  { return std::abs(point.x + point.y - 1.0) <= tol; }))
        return 4;
      break;
    case mfem::Geometry::PYRAMID:
      // MFEM's reference pyramid has base z=0 and apex (0,0,1), hence its four side
      // planes are x=0, y=0, x+z=1, and y+z=1.
      if (OnFacet([&](const auto &point) { return std::abs(point.z) <= tol; }))
        return 0;
      if (OnFacet([&](const auto &point) { return std::abs(point.x) <= tol; }))
        return 1;
      if (OnFacet([&](const auto &point) { return std::abs(point.y) <= tol; }))
        return 2;
      if (OnFacet([&](const auto &point)
                  { return std::abs(point.x + point.z - 1.0) <= tol; }))
        return 3;
      if (OnFacet([&](const auto &point)
                  { return std::abs(point.y + point.z - 1.0) <= tol; }))
        return 4;
      break;
    default:
      break;
  }
  return -1;
}

void Canonicalize(const std::vector<mfem::IntegrationPoint> &original,
                  std::vector<mfem::IntegrationPoint> &canonical,
                  std::vector<int> &canonical_to_original)
{
  canonical_to_original.resize(original.size());
  std::iota(canonical_to_original.begin(), canonical_to_original.end(), 0);
  for (const auto &ip : original)
  {
    MFEM_VERIFY(std::isfinite(ip.x) && std::isfinite(ip.y) && std::isfinite(ip.z) &&
                    std::isfinite(ip.weight),
                "Non-finite mapped reference point cannot be canonically routed!");
  }
  std::stable_sort(canonical_to_original.begin(), canonical_to_original.end(),
                   [&](int i, int j)
                   {
                     mfem::IntegrationPoint a = original[static_cast<std::size_t>(i)];
                     mfem::IntegrationPoint b = original[static_cast<std::size_t>(j)];
                     NormalizeReferencePoint(a);
                     NormalizeReferencePoint(b);
                     return CanonicalPointLess(a, b);
                   });
  canonical.resize(original.size());
  for (std::size_t q = 0; q < canonical.size(); q++)
  {
    canonical[q] = original[static_cast<std::size_t>(canonical_to_original[q])];
    NormalizeReferencePoint(canonical[q]);
  }
}

std::vector<long long>
MakeCanonicalTracePointKey(mfem::Geometry::Type vol_geom, mfem::Geometry::Type face_geom,
                           const mfem::IntegrationRule &face_ir,
                           const std::vector<mfem::IntegrationPoint> &canonical_pts)
{
  std::vector<long long> key;
  key.reserve(4 + 4 * canonical_pts.size());
  key.push_back(static_cast<long long>(vol_geom));
  key.push_back(static_cast<long long>(face_geom));
  key.push_back(static_cast<long long>(face_ir.GetOrder()));
  AppendPointRuleSignature(key, canonical_pts);
  return key;
}

}  // namespace

std::atomic<long long> FaceSamplingPlan::builds{0};
std::atomic<long long> FaceSamplingPlan::reuses{0};

FaceSamplingPlan::FaceSamplingPlan(const Mesh &mesh_,
                                   const mfem::Array<int> &bdr_attr_marker, int lod_)
  : mesh(&mesh_), mesh_sequence(mesh_.Get().GetSequence()),
    nodes_sequence(mesh_.Get().GetNodesSequence()), lod(lod_)
{
  MFEM_VERIFY(lod >= 0, "Boundary face sampling level of detail must be non-negative!");
  marker.assign(bdr_attr_marker.begin(), bdr_attr_marker.end());

  const mfem::ParMesh &pmesh = mesh->Get();
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    const int attr = pmesh.GetBdrAttribute(i);
    MFEM_VERIFY(attr > 0 && attr <= static_cast<int>(marker.size()),
                "Boundary face sampling marker does not cover a boundary attribute!");
    if (!marker[static_cast<std::size_t>(attr - 1)])
    {
      continue;
    }

    Entry entry;
    entry.bdr_elem = i;
    entry.attr = attr;
    entry.flip = BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
        i, pmesh, FET, T1, T2);
    entry.elem_a = FET.Elem1No;
    entry.elem_b = FET.Elem2 ? FET.Elem2No : -1;
    entry.ghost_b = entry.elem_b >= pmesh.GetNE();
    entry.bdr_geom = pmesh.GetBdrElementGeometry(i);
    entry.vol_geom_a = pmesh.GetElementGeometry(entry.elem_a);
    entry.vol_geom_b = entry.elem_b < 0 ? mfem::Geometry::INVALID
                       : entry.ghost_b  ? FET.Elem2->GetGeometryType()
                                        : pmesh.GetElementGeometry(entry.elem_b);
    if (entry.ghost_b)
    {
      entry.face_nbr = entry.elem_b - pmesh.GetNE();
      entry.ghost_attr = FET.Elem2->Attribute;
    }

    const mfem::IntegrationRule &face_ir =
        mfem::GlobGeometryRefiner.Refine(entry.bdr_geom, lod, 1)->RefPts;
    const int nq = face_ir.GetNPoints();
    int f, orientation;
    pmesh.GetBdrElementFace(i, &f, &orientation);
    auto MapPoints = [&](mfem::IntegrationPointTransformation &loc,
                         std::vector<mfem::IntegrationPoint> &points)
    {
      points.resize(nq);
      for (int q = 0; q < nq; q++)
      {
        mfem::IntegrationPoint fip = mfem::Mesh::TransformBdrElementToFace(
            FET.GetGeometryType(), orientation, face_ir.IntPoint(q));
        loc.Transform(fip, points[q]);
        points[q].weight = face_ir.IntPoint(q).weight;
      }
    };
    std::vector<mfem::IntegrationPoint> original_a, original_b;
    MapPoints(FET.Loc1, original_a);
    Canonicalize(original_a, entry.points_a, entry.canonical_to_original_a);
    entry.reference_face_a = ReferenceFacet(entry.vol_geom_a, entry.points_a);
    entry.point_key_a = MakeCanonicalTracePointKey(entry.vol_geom_a, entry.bdr_geom,
                                                   face_ir, entry.points_a);
    if (entry.elem_b >= 0)
    {
      MapPoints(FET.Loc2, original_b);
      Canonicalize(original_b, entry.points_b, entry.canonical_to_original_b);
      entry.reference_face_b = ReferenceFacet(entry.vol_geom_b, entry.points_b);
      entry.point_key_b = MakeCanonicalTracePointKey(entry.vol_geom_b, entry.bdr_geom,
                                                     face_ir, entry.points_b);
    }

    // Face Jacobians are shared setup state. Keep their legacy q order so all public
    // output ordering is preserved by the stored canonical-to-original permutations.
    const int face_jac_size = pmesh.SpaceDimension() * (pmesh.Dimension() - 1);
    entry.face_jacobians.resize(face_jac_size * nq);
    mfem::IsoparametricTransformation T;
    pmesh.GetBdrElementTransformation(i, &T);
    for (int q = 0; q < nq; q++)
    {
      const auto &ip = face_ir.IntPoint(q);
      T.SetIntPoint(&ip);
      const mfem::DenseMatrix &J = T.Jacobian();
      for (int d = 0; d < J.Width(); d++)
      {
        for (int c = 0; c < J.Height(); c++)
        {
          entry.face_jacobians[(c + pmesh.SpaceDimension() * d) * nq + q] = J(c, d);
        }
      }
    }

    if (buffer_bases.empty())
    {
      buffer_bases.resize(pmesh.GetNBE(), -1);
    }
    buffer_bases[i] = num_points;
    num_points += nq;
    entries.push_back(std::move(entry));
  }

  // Coalesce only NC coarse-owner leaf routes. A group is activated only when the
  // coordinate-only union actually removes at least one duplicate point, so conforming
  // and no-benefit configurations retain the simpler per-leaf fixed-rule route.
  using OwnerKey = std::tuple<int, int, int, int, int>;
  std::map<OwnerKey, std::vector<int>> candidates;
  if (pmesh.Nonconforming())
  {
    for (std::size_t e = 0; e < entries.size(); e++)
    {
      const auto &entry = entries[e];
      auto AddCandidate = [&](int side, bool ghost, int owner,
                              mfem::Geometry::Type geometry, int reference_face)
      {
        // The retained boundary element describes an NC leaf while FaceInformation can
        // describe a different active mesh face. Equal owner plus reference facet is
        // the stable coarse-owner signature; it also rejects unrelated faces sharing a
        // vertex or edge.
        if (owner >= 0 && reference_face >= 0)
        {
          candidates[{static_cast<int>(ghost), owner, side, static_cast<int>(geometry),
                      reference_face}]
              .push_back(static_cast<int>(e));
        }
      };
      AddCandidate(0, false, entry.elem_a, entry.vol_geom_a, entry.reference_face_a);
      AddCandidate(1, entry.ghost_b, entry.ghost_b ? entry.face_nbr : entry.elem_b,
                   entry.vol_geom_b, entry.reference_face_b);
    }
  }
  for (const auto &[key, candidate_entries] : candidates)
  {
    if (candidate_entries.size() < 2)
    {
      continue;
    }
    const bool ghost = std::get<0>(key) != 0;
    const int owner = std::get<1>(key);
    const int side = std::get<2>(key);
    const auto geometry = static_cast<mfem::Geometry::Type>(std::get<3>(key));
    std::map<CoordinateKey, mfem::IntegrationPoint> coordinate_points;
    int leaf_points = 0;
    for (const int entry_index : candidate_entries)
    {
      const auto &entry = entries[static_cast<std::size_t>(entry_index)];
      const auto &points = side == 0 ? entry.points_a : entry.points_b;
      leaf_points += static_cast<int>(points.size());
      for (auto ip : points)
      {
        const auto coordinate_key = MakeCoordinateKey(ip);
        NormalizeReferencePoint(ip);
        // Union point values are independent of quadrature weights. Use a fixed
        // representative so a later face-neighbor request also has a deterministic
        // rule signature without letting leaf weights fragment the union.
        ip.weight = 1.0;
        coordinate_points.emplace(coordinate_key, ip);
      }
    }
    if (static_cast<int>(coordinate_points.size()) >= leaf_points)
    {
      continue;
    }

    UnionGroup group;
    group.ghost = ghost;
    group.owner = owner;
    group.side = side;
    group.geometry = geometry;
    std::map<CoordinateKey, int> coordinate_to_union;
    for (const auto &[coordinate_key, point] : coordinate_points)
    {
      coordinate_to_union.emplace(coordinate_key, static_cast<int>(group.points.size()));
      group.points.push_back(point);
    }
    for (const int entry_index : candidate_entries)
    {
      const auto &entry = entries[static_cast<std::size_t>(entry_index)];
      const auto &points = side == 0 ? entry.points_a : entry.points_b;
      UnionRoute route;
      route.entry = entry_index;
      route.canonical_to_union.reserve(points.size());
      for (const auto &point : points)
      {
        const auto it = coordinate_to_union.find(MakeCoordinateKey(point));
        MFEM_VERIFY(it != coordinate_to_union.end(),
                    "Missing NC coarse-face union coordinate route!");
        route.canonical_to_union.push_back(it->second);
      }
      group.routes.push_back(std::move(route));
    }
    const int group_index = static_cast<int>(union_groups.size());
    for (const int entry_index : candidate_entries)
    {
      auto &entry = entries[static_cast<std::size_t>(entry_index)];
      if (side == 0)
      {
        entry.nc_coarse_a = true;
        entry.union_group_a = group_index;
      }
      else
      {
        entry.nc_coarse_b = true;
        entry.union_group_b = group_index;
      }
    }
    union_groups.push_back(std::move(group));
  }

  // Ghost unions are rank-local plan entries, but their exchange is collective. Record
  // the communicator-wide fact once while all ranks build the shared sampling plan so
  // each source-space trace entry can construct an empty peer exchange where needed
  // without another setup collective or a rank-local early return.
  has_ghost_union = std::any_of(union_groups.begin(), union_groups.end(),
                                [](const auto &group) { return group.ghost; });
  Mpi::GlobalOr(1, &has_ghost_union, mesh->GetComm());

  for (const auto &entry : entries)
  {
    MFEM_VERIFY(entry.bdr_elem >= 0 && entry.bdr_elem < pmesh.GetNBE() &&
                    entry.elem_a >= 0 && entry.elem_a < pmesh.GetNE(),
                "Invalid boundary entry after NC coarse-face union planning!");
  }

  ++builds;
}

bool FaceSamplingPlan::Matches(const Mesh &other, const mfem::Array<int> &bdr_attr_marker,
                               int other_lod) const
{
  if (mesh != &other || mesh_sequence != other.Get().GetSequence() ||
      nodes_sequence != other.Get().GetNodesSequence() || lod != other_lod ||
      marker.size() != static_cast<std::size_t>(bdr_attr_marker.Size()))
  {
    return false;
  }
  return std::equal(marker.begin(), marker.end(), bdr_attr_marker.begin());
}

long long FaceSamplingPlan::BuildCount()
{
  return builds.load();
}

long long FaceSamplingPlan::ReuseCount()
{
  return reuses.load();
}

void FaceSamplingPlan::NoteReuse()
{
  ++reuses;
}

}  // namespace palace
