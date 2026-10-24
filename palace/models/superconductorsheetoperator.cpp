// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superconductorsheetoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>
#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

// The two faces of a cracked two-sided film are paired by their vertex coordinates,
// quantized to 1e-4 of the shortest film edge, far coarser than the crack displacement
// between them. A face key is the vertex count followed by the vertex keys in canonical
// order.
using VertexKey = std::array<std::int64_t, 3>;
using FaceKey = std::array<std::int64_t, 13>;

double GetKeyTolerance(const mfem::ParMesh &mesh, const std::unordered_set<int> &attrs)
{
  double h = mfem::infinity();
  mfem::Array<int> edges, orient, ev;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (!attrs.contains(mesh.GetBdrAttribute(be)))
    {
      continue;
    }
    mesh.GetBdrElementEdges(be, edges, orient);
    for (int e : edges)
    {
      mesh.GetEdgeVertices(e, ev);
      const double *x0 = mesh.GetVertex(ev[0]), *x1 = mesh.GetVertex(ev[1]);
      h = std::min(h, std::sqrt((x1[0] - x0[0]) * (x1[0] - x0[0]) +
                                (x1[1] - x0[1]) * (x1[1] - x0[1]) +
                                (x1[2] - x0[2]) * (x1[2] - x0[2])));
    }
  }
  Mpi::GlobalMin(1, &h, mesh.GetComm());
  return 1.0e-4 * ((h < mfem::infinity()) ? h : 1.0);
}

VertexKey GetVertexKey(const double *x, double tol)
{
  return {std::llround(x[0] / tol), std::llround(x[1] / tol), std::llround(x[2] / tol)};
}

// The rank that pairs records with the given key (FNV-1a, identical on every rank).
template <std::size_t N>
int GetKeyRank(const std::array<std::int64_t, N> &key, int nranks)
{
  std::uint64_t h = 14695981039346656037ULL;
  for (auto k : key)
  {
    h = (h ^ static_cast<std::uint64_t>(k)) * 1099511628211ULL;
  }
  return static_cast<int>(h % static_cast<std::uint64_t>(nranks));
}

// Key of a film face, with its vertices in a canonical order shared by the two faces of a
// pair: the smallest vertex key first, then (quadrilaterals) around the face toward the
// smaller neighbor. canon[k] is the local vertex in canonical position k.
FaceKey GetFaceKey(const mfem::ParMesh &mesh, int be, double tol, std::array<int, 4> &canon)
{
  mfem::Array<int> vids;
  mesh.GetBdrElementVertices(be, vids);
  const int nv = vids.Size();
  MFEM_VERIFY(nv == 3 || nv == 4,
              "Two-sided superconductor sheet faces must be triangles or quadrilaterals!");
  std::array<VertexKey, 4> vk{};
  for (int k = 0; k < nv; k++)
  {
    vk[k] = GetVertexKey(mesh.GetVertex(vids[k]), tol);
  }
  if (nv == 3)
  {
    canon = {0, 1, 2, 3};
    std::sort(canon.begin(), canon.begin() + 3,
              [&vk](int p, int q) { return vk[p] < vk[q]; });
  }
  else
  {
    const int s = static_cast<int>(std::min_element(vk.begin(), vk.end()) - vk.begin());
    const int step = (vk[(s + 1) % 4] < vk[(s + 3) % 4]) ? 1 : 3;
    for (int k = 0; k < 4; k++)
    {
      canon[k] = (s + step * k) % 4;
    }
  }
  FaceKey key{};
  key[0] = nv;
  for (int k = 0; k < nv; k++)
  {
    std::copy(vk[canon[k]].begin(), vk[canon[k]].end(), key.begin() + 1 + 3 * k);
  }
  return key;
}

// Element reference coordinates of a point given in the canonical frame of a face. The
// canonical order is a symmetry of the reference face, so the map is affine with unit
// Jacobian and the two faces of a pair map it to the same physical point.
mfem::IntegrationPoint CanonicalToElement(mfem::Geometry::Type geom,
                                          const std::array<int, 4> &canon,
                                          const mfem::IntegrationPoint &ip)
{
  const mfem::IntegrationRule &verts = *mfem::Geometries.GetVertices(geom);
  const double x = ip.x, y = ip.y;
  const std::array<double, 4> N =
      (geom == mfem::Geometry::TRIANGLE)
          ? std::array<double, 4>{1.0 - x - y, x, y, 0.0}
          : std::array<double, 4>{(1.0 - x) * (1.0 - y), x * (1.0 - y), x * y,
                                  (1.0 - x) * y};
  mfem::IntegrationPoint ipe;
  ipe.Init(0);
  for (int k = 0; k < verts.GetNPoints(); k++)
  {
    ipe.x += N[k] * verts.IntPoint(canon[k]).x;
    ipe.y += N[k] * verts.IntPoint(canon[k]).y;
  }
  ipe.weight = ip.weight;
  return ipe;
}

// Send send[r] to rank r for every rank r, returning everything received.
std::vector<double> Exchange(const std::vector<std::vector<double>> &send, MPI_Comm comm)
{
  const int nranks = Mpi::Size(comm);
  std::vector<int> scounts(nranks), sdispls(nranks), rcounts(nranks), rdispls(nranks);
  std::vector<double> sbuf;
  for (int r = 0; r < nranks; r++)
  {
    scounts[r] = static_cast<int>(send[r].size());
    sdispls[r] = static_cast<int>(sbuf.size());
    sbuf.insert(sbuf.end(), send[r].begin(), send[r].end());
  }
  Mpi::Alltoall(1, scounts.data(), rcounts.data(), comm);
  int rtotal = 0;
  for (int r = 0; r < nranks; r++)
  {
    rdispls[r] = rtotal;
    rtotal += rcounts[r];
  }
  std::vector<double> rbuf(rtotal);
  Mpi::Alltoallv(sbuf.data(), scounts.data(), sdispls.data(), rbuf.data(), rcounts.data(),
                 rdispls.data(), comm);
  return rbuf;
}

// Barycentric test that p lies in the face with vertex coordinates x (orthogonal
// projection onto each triangle of the face).
bool FaceContains(const std::vector<std::array<double, 3>> &x, const double *p)
{
  auto in_triangle = [p](const auto &a, const auto &b, const auto &c)
  {
    double v0[3], v1[3], v2[3];
    for (int d = 0; d < 3; d++)
    {
      v0[d] = b[d] - a[d];
      v1[d] = c[d] - a[d];
      v2[d] = p[d] - a[d];
    }
    auto dot = [](const double *u, const double *v)
    { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
    const double d00 = dot(v0, v0), d01 = dot(v0, v1), d11 = dot(v1, v1), d20 = dot(v2, v0),
                 d21 = dot(v2, v1), det = d00 * d11 - d01 * d01;
    const double b1 = (d11 * d20 - d01 * d21) / det, b2 = (d00 * d21 - d01 * d20) / det;
    constexpr double eps = 1.0e-8;
    return b1 >= -eps && b2 >= -eps && b1 + b2 <= 1.0 + eps;
  };
  return in_triangle(x[0], x[1], x[2]) || (x.size() == 4 && in_triangle(x[0], x[2], x[3]));
}

}  // namespace

SuperconductorSheetOperator::SuperconductorSheetOperator(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(superconductor, cracked_attributes, mesh);
  PrintBoundaryInfo(units, mesh);
}

SuperconductorSheetOperator::SuperconductorSheetOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         const mfem::ParMesh &mesh)
  : SuperconductorSheetOperator(iodata.boundaries.superconductor,
                                iodata.boundaries.cracked_attributes, iodata.units, mat_op,
                                mesh)
{
}

void SuperconductorSheetOperator::SetUpBoundaryProperties(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const mfem::ParMesh &mesh)
{
  // Check that superconductor sheet boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!superconductor.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, superconductor, "superconductor sheet");
  }

  // Kinetic sheet inductance L_ksq [H/sq]: supplied directly, or from (lambda, d) via the
  // finite-thickness London sheet L_ksq = lambda*coth(d/lambda) (nondimensional, mu0
  // absorbed). Reduces to the thin-film Pearl limit lambda^2/d for d << lambda and
  // saturates at lambda for d >> lambda, capturing the through-thickness current profile
  // without resolving the film thickness.
  boundaries.reserve(superconductor.size());
  for (const auto &data : superconductor)
  {
    // Two-sided (two-port) sheet: the diagonal self-term inductance is λ·tanh(d/λ), so
    // 1/Ls = coth(d/λ)/λ; the cross-face term is added separately via
    // BuildTwoPortCoupling. Otherwise the single-sheet L_ksq applies.
    double Ls;
    if (data.two_sided)
    {
      Ls = data.lambda_L * std::tanh(data.thickness / data.lambda_L);
      has_two_port_ = true;
    }
    else
    {
      Ls =
          (data.Ls > 0.0) ? data.Ls : KineticSheetInductance(data.lambda_L, data.thickness);
    }
    MFEM_VERIFY(Ls > 0.0,
                "Superconductor sheet has non-positive kinetic sheet inductance!");
    auto &bdr = boundaries.emplace_back();
    bdr.Ls = Ls;
    bdr.lambda_L = data.lambda_L;
    bdr.thickness = data.thickness;
    bdr.two_sided = data.two_sided;
    bdr.attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      if (!mesh_attrs.Contains(attr))
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
      // Per-attribute scaling for the doubled area of a cracked interface. A two-sided
      // sheet keeps each face at its full self-term (the two faces are physically
      // distinct), so it is not halved even though it cracks.
      bdr.attr_scaling[attr] =
          (!data.two_sided && cracked_attributes.contains(attr)) ? 2.0 : 1.0;
    }
  }
  MFEM_VERIFY(!has_two_port_ || mesh.Dimension() == 3,
              "Superconductor boundary \"TwoSided\" requires a 3D mesh!");
}

void SuperconductorSheetOperator::PrintBoundaryInfo(const Units &units,
                                                    const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  fmt::format_to(out, "\nConfiguring superconductor sheet BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      if (bdr.two_sided)
      {
        // Each face's self-term, and the single sheet ½λ·coth(d/2λ) the pair reduces to
        // when both faces carry the same field.
        fmt::format_to(
            out,
            " {:d}: L_ksq = {:.3e} H/sq per face, {:.3e} H/sq common mode, n = ({:+.1f})\n",
            attr, units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
            units.Dimensionalize<VT::INDUCTANCE>(
                0.5 * KineticSheetInductance(bdr.lambda_L, 0.5 * bdr.thickness)),
            fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
        continue;
      }
      fmt::format_to(out, " {:d}: L_ksq = {:.3e} H/sq, n = ({:+.1f})\n", attr,
                     units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buffer));
}

double SuperconductorSheetOperator::KineticSheetInductance(double lambda_L,
                                                           double thickness)
{
  // Finite-thickness London sheet inductance L_ksq = lambda * coth(d/lambda), with
  // coth(x) = 1/tanh(x). Reduces to lambda^2/d for d << lambda; saturates at lambda for
  // d >> lambda.
  return lambda_L / std::tanh(thickness / lambda_L);
}

std::unique_ptr<mfem::HypreParMatrix> SuperconductorSheetOperator::BuildTwoPortCoupling(
    mfem::ParFiniteElementSpace &nd_fespace) const
{
  if (!has_two_port_)
  {
    return nullptr;
  }
  mfem::ParMesh &pmesh = *nd_fespace.GetParMesh();
  MPI_Comm comm = nd_fespace.GetComm();
  const int sd = pmesh.SpaceDimension(), rank = Mpi::Rank(comm), nranks = Mpi::Size(comm);

  // Cross-face coupling coefficient -1/(λ·sinh(d/λ)) per two-sided attribute.
  std::unordered_map<int, double> attr_off;
  std::unordered_set<int> attrs;
  for (const auto &bdr : boundaries)
  {
    if (!bdr.two_sided)
    {
      continue;
    }
    const double off = -1.0 / (bdr.lambda_L * std::sinh(bdr.thickness / bdr.lambda_L));
    for (auto attr : bdr.attr_list)
    {
      attr_off[attr] = off;
      attrs.insert(attr);
    }
  }
  const double tol = GetKeyTolerance(pmesh, attrs);

  // Fixed-size face records, padded to the largest face element and quadrature rule:
  // key | off | nd | nqp | rank | local face | nd global L-dofs | nd signs | nqp·nd·sd
  // basis | nqp weights. The quadrature is set in the canonical frame of each face, so the
  // two faces of a pair evaluate their bases at the same physical points.
  const int order = nd_fespace.GetMaxElementOrder();
  int nd_max = 0, nqp_max = 0;
  for (auto geom : {mfem::Geometry::TRIANGLE, mfem::Geometry::SQUARE})
  {
    nd_max =
        std::max(nd_max, nd_fespace.FEColl()->FiniteElementForGeometry(geom)->GetDof());
    nqp_max = std::max(nqp_max, mfem::IntRules.Get(geom, 2 * order + 2).GetNPoints());
  }
  const int nk = std::tuple_size_v<FaceKey>, o_off = nk, o_nd = nk + 1, o_nqp = nk + 2,
            o_rank = nk + 3, o_face = nk + 4, o_dof = nk + 5, o_sgn = o_dof + nd_max,
            o_vsh = o_sgn + nd_max, o_w = o_vsh + nqp_max * nd_max * sd,
            rec = o_w + nqp_max;

  // Pack the local film faces and route each to the rank that pairs its key.
  const HYPRE_BigInt ldof_offset = nd_fespace.GetMyDofOffset();
  std::vector<double> local;
  std::vector<std::vector<double>> send(nranks);
  mfem::DenseMatrix vsh;
  mfem::Array<int> dofs;
  mfem::DofTransformation doftrans;
  std::array<int, 4> canon;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    auto it = attr_off.find(pmesh.GetBdrAttribute(be));
    if (it == attr_off.end())
    {
      continue;
    }
    const FaceKey key = GetFaceKey(pmesh, be, tol, canon);
    const mfem::FiniteElement &fe = *nd_fespace.GetBE(be);
    const mfem::IntegrationRule &ir = mfem::IntRules.Get(fe.GetGeomType(), 2 * order + 2);
    const int nd = fe.GetDof(), nqp = ir.GetNPoints();
    mfem::ElementTransformation &T = *pmesh.GetBdrElementTransformation(be);
    nd_fespace.GetBdrElementDofs(be, dofs, doftrans);

    const std::size_t r = local.size();
    local.resize(r + rec, 0.0);
    std::copy(key.begin(), key.end(), local.begin() + r);
    local[r + o_off] = it->second;
    local[r + o_nd] = nd;
    local[r + o_nqp] = nqp;
    local[r + o_rank] = rank;
    local[r + o_face] = static_cast<double>(r / rec);
    for (int a = 0; a < nd; a++)
    {
      const int ld = (dofs[a] >= 0) ? dofs[a] : -1 - dofs[a];
      local[r + o_dof + a] = static_cast<double>(ldof_offset + ld);
      local[r + o_sgn + a] = (dofs[a] >= 0) ? 1.0 : -1.0;
    }
    vsh.SetSize(nd, sd);
    for (int q = 0; q < nqp; q++)
    {
      const mfem::IntegrationPoint ipe =
          CanonicalToElement(fe.GetGeomType(), canon, ir.IntPoint(q));
      T.SetIntPoint(&ipe);
      fe.CalcVShape(T, vsh);
      // Map the local basis to the global DOF orientation, as in assembly.
      for (int d = 0; d < sd; d++)
      {
        doftrans.TransformDual(vsh.GetColumn(d));
      }
      local[r + o_w + q] = ipe.weight * T.Weight();
      for (int a = 0; a < nd; a++)
      {
        for (int d = 0; d < sd; d++)
        {
          local[r + o_vsh + (q * nd + a) * sd + d] = vsh(a, d);
        }
      }
    }
    auto &buf = send[GetKeyRank(key, nranks)];
    buf.insert(buf.end(), local.begin() + r, local.end());
  }

  // Pair records by key and send each face its twin's record, prefixed with the face's
  // local index on its own rank.
  const std::vector<double> paired = Exchange(send, comm);
  std::map<FaceKey, std::vector<std::size_t>> by_key;
  for (std::size_t r = 0; r < paired.size(); r += rec)
  {
    FaceKey key;
    std::copy(paired.begin() + r, paired.begin() + r + nk, key.begin());
    by_key[key].push_back(r);
  }
  int n_unpaired = 0;
  for (auto &v : send)
  {
    v.clear();
  }
  for (const auto &[key, recs] : by_key)
  {
    if (recs.size() != 2)
    {
      n_unpaired += static_cast<int>(recs.size());
      continue;
    }
    for (int i = 0; i < 2; i++)
    {
      const std::size_t R = recs[i], S = recs[1 - i];
      auto &buf = send[static_cast<int>(paired[R + o_rank])];
      buf.push_back(paired[R + o_face]);
      buf.insert(buf.end(), paired.begin() + S, paired.begin() + S + rec);
    }
  }
  Mpi::GlobalSum(1, &n_unpaired, comm);
  MFEM_VERIFY(n_unpaired == 0, "Two-sided superconductor sheet has "
                                   << n_unpaired
                                   << " face(s) without a coincident twin; the film must "
                                      "be fully interior!");
  const std::vector<double> twins = Exchange(send, comm);

  // Each rank assembles C(row, col) = off · s_a · s_b · ∫ φ_a·φ_b over its own faces, rows
  // from the face and columns from its twin, as a parallel matrix on the L-dofs.
  std::vector<std::tuple<int, HYPRE_BigInt, double>> entries;
  for (std::size_t t = 0; t < twins.size(); t += rec + 1)
  {
    const double *R = local.data() + static_cast<std::size_t>(twins[t]) * rec;
    const double *S = twins.data() + t + 1;
    const int nd = static_cast<int>(R[o_nd]), nqp = static_cast<int>(R[o_nqp]);
    for (int a = 0; a < nd; a++)
    {
      const int row =
          static_cast<int>(static_cast<HYPRE_BigInt>(R[o_dof + a]) - ldof_offset);
      for (int b = 0; b < nd; b++)
      {
        double m_ab = 0.0;
        for (int q = 0; q < nqp; q++)
        {
          double dot = 0.0;
          for (int d = 0; d < sd; d++)
          {
            dot += R[o_vsh + (q * nd + a) * sd + d] * S[o_vsh + (q * nd + b) * sd + d];
          }
          m_ab += R[o_w + q] * dot;
        }
        entries.emplace_back(row, static_cast<HYPRE_BigInt>(S[o_dof + b]),
                             R[o_off] * R[o_sgn + a] * S[o_sgn + b] * m_ab);
      }
    }
  }
  std::sort(entries.begin(), entries.end());
  const int nrows = nd_fespace.GetVSize();
  std::vector<int> I(nrows + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> data;
  for (std::size_t k = 0; k < entries.size(); k++)
  {
    const auto &[row, col, val] = entries[k];
    if (k > 0 && std::get<0>(entries[k - 1]) == row && std::get<1>(entries[k - 1]) == col)
    {
      data.back() += val;
      continue;
    }
    I[row + 1]++;
    J.push_back(col);
    data.push_back(val);
  }
  for (int i = 0; i < nrows; i++)
  {
    I[i + 1] += I[i];
  }
  const HYPRE_BigInt nglob = nd_fespace.GlobalVSize();
  mfem::HypreParMatrix C_l(comm, nrows, nglob, nglob, I.data(), J.data(), data.data(),
                           nd_fespace.GetDofOffsets(), nd_fespace.GetDofOffsets());

  // C = Pᵀ C_l P maps to the true DOFs, including shared and nonconforming (slave) DOFs.
  return std::unique_ptr<mfem::HypreParMatrix>(
      mfem::RAP(&C_l, nd_fespace.Dof_TrueDof_Matrix()));
}

int SuperconductorSheetOperator::MatchTwoSidedRefinement(
    const std::vector<config::SuperconductorData> &superconductor, mfem::ParMesh &mesh,
    int nc_limit)
{
  std::unordered_set<int> attrs;
  for (const auto &data : superconductor)
  {
    if (data.two_sided)
    {
      attrs.insert(data.attributes.begin(), data.attributes.end());
    }
  }
  if (attrs.empty())
  {
    return 0;
  }
  MPI_Comm comm = mesh.GetComm();
  const int rank = Mpi::Rank(comm), nranks = Mpi::Size(comm);
  const auto initial_elem_count = mesh.GetGlobalNE();

  // Refining the elements of one face leaves its twin coarser. A film face with no
  // coincident twin is coarse if a smaller face of the other side lies inside it: such a
  // face touches one of its vertices, so faces are routed by vertex key and compared
  // there. Record layout: vertex key | face key | nv | 4 × 3 vertex coordinates | rank |
  // element.
  constexpr int nk = std::tuple_size_v<FaceKey>, o_face = 3, o_nv = o_face + nk,
                o_x = o_nv + 1, o_rank = o_x + 12, o_elem = o_rank + 1, rec = o_elem + 1;
  std::array<int, 4> canon;
  mfem::Array<int> vids;
  while (true)
  {
    const double tol = GetKeyTolerance(mesh, attrs);
    std::vector<std::vector<double>> send(nranks);
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      if (!attrs.contains(mesh.GetBdrAttribute(be)))
      {
        continue;
      }
      const FaceKey key = GetFaceKey(mesh, be, tol, canon);
      mesh.GetBdrElementVertices(be, vids);
      int elem, info;
      mesh.GetBdrElementAdjacentElement(be, elem, info);
      std::vector<double> r(rec, 0.0);
      std::copy(key.begin(), key.end(), r.begin() + o_face);
      r[o_nv] = vids.Size();
      for (int k = 0; k < vids.Size(); k++)
      {
        std::copy(mesh.GetVertex(vids[k]), mesh.GetVertex(vids[k]) + 3,
                  r.begin() + o_x + 3 * k);
      }
      r[o_rank] = rank;
      r[o_elem] = elem;
      for (int k = 0; k < vids.Size(); k++)
      {
        const VertexKey vk = GetVertexKey(mesh.GetVertex(vids[k]), tol);
        std::copy(vk.begin(), vk.end(), r.begin());
        auto &buf = send[GetKeyRank(vk, nranks)];
        buf.insert(buf.end(), r.begin(), r.end());
      }
    }
    const std::vector<double> recv = Exchange(send, comm);

    // At each vertex, mark every unmatched face that contains a smaller face of the other
    // side. Unmatched faces are counted once, at their first canonical vertex.
    std::map<VertexKey, std::vector<std::size_t>> by_vertex;
    for (std::size_t r = 0; r < recv.size(); r += rec)
    {
      by_vertex[{static_cast<std::int64_t>(recv[r]), static_cast<std::int64_t>(recv[r + 1]),
                 static_cast<std::int64_t>(recv[r + 2])}]
          .push_back(r);
    }
    auto face_key = [&recv](std::size_t r)
    {
      FaceKey key;
      std::copy(recv.begin() + r + o_face, recv.begin() + r + o_face + nk, key.begin());
      return key;
    };
    // Centroid c and normal n of a face, with |n| twice its area: the cross product of two
    // edges of a triangle, or of the diagonals of a quadrilateral.
    auto face_geometry = [&recv](std::size_t r, double *c, double *n)
    {
      const int nv = static_cast<int>(recv[r + o_nv]);
      const double *x = recv.data() + r + o_x;
      double u[3], v[3];
      for (int d = 0; d < 3; d++)
      {
        u[d] = (nv == 4) ? x[6 + d] - x[d] : x[3 + d] - x[d];
        v[d] = (nv == 4) ? x[9 + d] - x[3 + d] : x[6 + d] - x[d];
        c[d] = 0.0;
        for (int k = 0; k < nv; k++)
        {
          c[d] += x[3 * k + d] / nv;
        }
      }
      n[0] = u[1] * v[2] - u[2] * v[1];
      n[1] = u[2] * v[0] - u[0] * v[2];
      n[2] = u[0] * v[1] - u[1] * v[0];
    };
    auto norm = [](const double *n)
    { return std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]); };
    std::vector<std::vector<double>> marks(nranks);
    int n_unmatched = 0;
    for (const auto &[vk, recs] : by_vertex)
    {
      for (auto G : recs)
      {
        const FaceKey key_G = face_key(G);
        bool matched = false;
        for (auto F : recs)
        {
          matched |= (F != G && face_key(F) == key_G);
        }
        if (matched)
        {
          continue;
        }
        n_unmatched += (std::equal(vk.begin(), vk.end(), key_G.begin() + 1));
        double c_G[3], n_G[3];
        face_geometry(G, c_G, n_G);
        const int nv = static_cast<int>(recv[G + o_nv]);
        std::vector<std::array<double, 3>> x_G(nv);
        for (int k = 0; k < nv; k++)
        {
          std::copy(recv.begin() + G + o_x + 3 * k, recv.begin() + G + o_x + 3 * k + 3,
                    x_G[k].begin());
        }
        for (auto F : recs)
        {
          double c_F[3], n_F[3];
          face_geometry(F, c_F, n_F);
          const double a_F = norm(n_F), a_G = norm(n_G),
                       cos_FG = (n_F[0] * n_G[0] + n_F[1] * n_G[1] + n_F[2] * n_G[2]) /
                                (a_F * a_G);
          if (F != G && a_F < 0.75 * a_G && std::abs(cos_FG) > 0.5 &&
              FaceContains(x_G, c_F))
          {
            marks[static_cast<int>(recv[G + o_rank])].push_back(recv[G + o_elem]);
            break;
          }
        }
      }
    }
    Mpi::GlobalSum(1, &n_unmatched, comm);
    if (n_unmatched == 0)
    {
      break;
    }
    const std::vector<double> marked = Exchange(marks, comm);
    std::set<int> unique(marked.begin(), marked.end());
    mfem::Array<int> refine;
    for (int e : unique)
    {
      refine.Append(e);
    }
    int n_marked = refine.Size();
    Mpi::GlobalSum(1, &n_marked, comm);
    MFEM_VERIFY(n_marked > 0, "Two-sided superconductor sheet has "
                                  << n_unmatched
                                  << " face(s) whose twin is refined differently; the two "
                                     "faces of the film must refine alike!");
    mesh.GeneralRefinement(refine, -1, nc_limit);
  }
  return static_cast<int>(mesh.GetGlobalNE() - initial_elem_count);
}

mfem::Array<int> SuperconductorSheetOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

void SuperconductorSheetOperator::AddStiffnessBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fb) const
{
  // Kinetic sheet inductance boundaries: add (coeff / L_ksq) as a tangential surface mass.
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      const double s = bdr.attr_scaling.at(attr);
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Ls * s));
    }
  }
}

}  // namespace palace
