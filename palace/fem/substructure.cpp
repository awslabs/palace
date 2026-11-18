// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructure.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <set>
#include <nlohmann/json.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/outputdir.hpp"

namespace palace
{

namespace
{

// MarkInterfaceTrueDofs for the element DOFs given by element_dofs(e, dofs) into an
// L-vector of size vsize, with the prolongation P.
template <typename ElementDofs>
int MarkInterfaceTrueDofs(const mfem::ParMesh &mesh, int vsize,
                          const mfem::HypreParMatrix &P, const ElementDofs &element_dofs,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker)
{
  auto in = [](const mfem::Array<int> &s, int a)
  {
    for (int x : s)
    {
      if (x == a)
      {
        return true;
      }
    }
    return false;
  };

  // L-vector attribute markers over element DOFs.
  mfem::Vector reg_l(vsize), env_l(vsize);
  reg_l = 0.0;
  env_l = 0.0;
  mfem::Array<int> dofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int attr = mesh.GetAttribute(e);
    const bool is_reg = in(region_attrs, attr), is_env = in(environment_attrs, attr);
    if (!is_reg && !is_env)
    {
      continue;
    }
    element_dofs(e, dofs);
    for (int i = 0; i < dofs.Size(); i++)
    {
      const int d = dofs[i] >= 0 ? dofs[i] : -1 - dofs[i];
      if (is_reg)
      {
        reg_l(d) = 1.0;
      }
      else
      {
        env_l(d) = 1.0;
      }
    }
  }

  // Reduce to true DOFs with cross-rank accumulation via |P|^T.
  const int nt = P.Width();
  mfem::Vector reg_t(nt), env_t(nt);
  P.AbsMultTranspose(1.0, reg_l, 0.0, reg_t);
  P.AbsMultTranspose(1.0, env_l, 0.0, env_t);

  region_marker.SetSize(nt);
  env_marker.SetSize(nt);
  const double *rt = reg_t.HostRead(), *et = env_t.HostRead();
  int local_iface = 0;
  for (int i = 0; i < nt; i++)
  {
    region_marker[i] = (rt[i] > 0.5) ? 1 : 0;
    env_marker[i] = (et[i] > 0.5) ? 1 : 0;
    local_iface += region_marker[i] && env_marker[i];
  }
  int global_iface = 0;
  MPI_Allreduce(&local_iface, &global_iface, 1, MPI_INT, MPI_SUM, mesh.GetComm());
  return global_iface;
}

}  // namespace

void CheckSubstructuringAttributes(const mfem::ParMesh &mesh,
                                   const std::vector<int> &region_attributes,
                                   const std::vector<int> &environment_attributes)
{
  // The attributes of a ParMesh are those of all processes.
  auto in = [](const std::vector<int> &s, int a)
  { return std::ranges::find(s, a) != s.end(); };
  for (int a : mesh.attributes)
  {
    MFEM_VERIFY(in(region_attributes, a) || in(environment_attributes, a),
                "Domain attribute " << a
                                    << " of the mesh is in neither the substructuring "
                                       "\"Region\" nor the \"Environment\"!");
  }
  for (int a : region_attributes)
  {
    MFEM_VERIFY(mesh.attributes.Find(a) >= 0, "Substructuring \"Region\" attribute "
                                                  << a
                                                  << " is not a domain attribute of "
                                                     "the mesh!");
  }
}

int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker)
{
  auto element_dofs = [&](int e, mfem::Array<int> &dofs)
  { parent_fespace.GetElementDofs(e, dofs); };
  return MarkInterfaceTrueDofs(*parent_fespace.GetParMesh(), parent_fespace.GetVSize(),
                               *parent_fespace.Dof_TrueDof_Matrix(), element_dofs,
                               region_attrs, environment_attrs, region_marker, env_marker);
}

int MarkInterfaceTrueDofs(const FiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker)
{
  mfem::DofTransformation dof_trans;
  auto element_dofs = [&](int e, mfem::Array<int> &dofs)
  { parent_fespace.GetElementDofs(e, dofs, dof_trans); };
  return MarkInterfaceTrueDofs(parent_fespace.GetParMesh(), parent_fespace.GetVSize(),
                               *parent_fespace.Dof_TrueDof_Matrix(), element_dofs,
                               region_attrs, environment_attrs, region_marker, env_marker);
}

namespace
{

bool IsNedelec(const mfem::ParFiniteElementSpace &fespace)
{
  return dynamic_cast<const mfem::ND_FECollection *>(fespace.FEColl()) != nullptr;
}

}  // namespace

int SignatureWidth(const mfem::ParFiniteElementSpace &fespace)
{
  return IsNedelec(fespace) ? 12 : 3;
}

namespace
{

// The signatures of the local true DOFs i with index[i] >= 0, at rows index[i] of an
// n x width array (zero elsewhere).
std::vector<double> LocalTrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                           const std::vector<int> &index, int n)
{
  const int w = SignatureWidth(fespace), nt = fespace.GetTrueVSize();
  std::vector<double> loc(static_cast<std::size_t>(n) * w, 0.0);
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(fespace);
  mfem::ParGridFunction gf(&fes);
  mfem::Vector td(nt);
  auto stamp = [&](int slot)
  {
    gf.GetTrueDofs(td);
    const double *t = td.HostRead();
    for (int i = 0; i < nt; i++)
    {
      if (index[i] >= 0)
      {
        loc[static_cast<std::size_t>(index[i]) * w + slot] = t[i];
      }
    }
  };
  if (!IsNedelec(fespace))
  {
    for (int d = 0; d < fes.GetParMesh()->SpaceDimension(); d++)
    {
      mfem::FunctionCoefficient xc([d](const mfem::Vector &x) { return x(d); });
      gf.ProjectCoefficient(xc);
      stamp(d);
    }
  }
  else
  {
    for (int b = 0; b < 3; b++)
    {
      mfem::Vector e(3);
      e = 0.0;
      e(b) = 1.0;
      mfem::VectorConstantCoefficient ec(e);
      gf.ProjectCoefficient(ec);
      stamp(b);
    }
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        mfem::VectorFunctionCoefficient xc(3,
                                           [a, b](const mfem::Vector &x, mfem::Vector &v)
                                           {
                                             v = 0.0;
                                             v(b) = x(a);
                                           });
        gf.ProjectCoefficient(xc);
        stamp(3 + a * 3 + b);
      }
    }
  }
  return loc;
}

}  // namespace

std::vector<double> TrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<int> &index, int n)
{
  const std::vector<double> loc = LocalTrueDofSignatures(fespace, index, n);
  std::vector<double> glob(loc.size(), 0.0);
  MPI_Allreduce(loc.data(), glob.data(), static_cast<int>(loc.size()), MPI_DOUBLE, MPI_SUM,
                fespace.GetComm());
  return glob;
}

namespace
{

// The point p of each DOF from its signature: the coordinates of an H1 DOF, and for a
// point-tangent DOF u(p).t with the signature (t_b, p_a t_b) (and a side moment), p = Q t /
// |t|^2 with Q_ab = p_a t_b.
std::vector<std::array<double, 3>> SignaturePoints(const std::vector<double> &sig,
                                                   int width)
{
  const int n = static_cast<int>(sig.size() / width);
  std::vector<std::array<double, 3>> p(n, {0.0, 0.0, 0.0});
  for (int g = 0; g < n; g++)
  {
    const double *v = &sig[static_cast<std::size_t>(g) * width];
    if (width == 3)
    {
      p[g] = {v[0], v[1], v[2]};
      continue;
    }
    const double t2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    for (int a = 0; a < 3 && t2 > 0.0; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        p[g][a] += v[3 + 3 * a + b] * v[b] / t2;
      }
    }
  }
  return p;
}

// The matching tolerance for the signatures of a saved set.
double SignatureTolerance(const std::vector<double> &saved)
{
  double scale = 1.0;
  for (double v : saved)
  {
    scale = std::max(scale, std::abs(v));
  }
  return 1.0e-8 * scale;
}

}  // namespace

std::vector<char> MatchSignaturePoints(const std::vector<double> &candidates, int width,
                                       const std::vector<double> &saved, int saved_width)
{
  auto known = [](int w) { return w == 3 || w == 12 || w == 13; };
  MFEM_VERIFY(known(width) && known(saved_width) && (width == 3) == (saved_width == 3),
              "Unknown DOF signature widths " << width << " and " << saved_width << "!");
  const double tol = SignatureTolerance(saved);
  const auto pc = SignaturePoints(candidates, width),
             ps = SignaturePoints(saved, saved_width);
  // The saved points by their first coordinate, searched within the tolerance.
  std::vector<int> order(ps.size());
  for (std::size_t s = 0; s < order.size(); s++)
  {
    order[s] = static_cast<int>(s);
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) { return ps[a][0] < ps[b][0]; });
  std::vector<char> match(pc.size(), 0);
  for (std::size_t g = 0; g < pc.size(); g++)
  {
    auto it = std::lower_bound(order.begin(), order.end(), pc[g][0] - tol,
                               [&](int s, double x) { return ps[s][0] < x; });
    for (; it != order.end() && ps[*it][0] <= pc[g][0] + tol && !match[g]; ++it)
    {
      const auto &q = ps[*it];
      match[g] = std::hypot(pc[g][0] - q[0], pc[g][1] - q[1], pc[g][2] - q[2]) <= tol;
    }
  }
  return match;
}

std::vector<int> SignatureMap::SavedRows(int g0, int m) const
{
  std::vector<int> saved;
  for (int g = g0; g < g0 + m; g++)
  {
    for (const auto &[a, c] : rows[g])
    {
      saved.push_back(a);
    }
  }
  std::sort(saved.begin(), saved.end());
  saved.erase(std::unique(saved.begin(), saved.end()), saved.end());
  return saved;
}

SignatureMap MatchSignatureBasis(const std::vector<double> &current,
                                 const std::vector<double> &saved, int width)
{
  MFEM_VERIFY(width == 3 || width == 12 || width == 13,
              "Unknown DOF signature width " << width << "!");
  MFEM_VERIFY(current.size() == saved.size(),
              "Saved and current DOF sets of different sizes cannot be matched!");
  const int n = static_cast<int>(current.size() / width);
  const bool nedelec = (width != 3);
  const double tol = SignatureTolerance(saved);
  auto row = [width](const std::vector<double> &sig, int g)
  { return &sig[static_cast<std::size_t>(g) * width]; };

  // One-to-one matches, signed for H(curl), among the saved DOFs whose key |w . s|
  // (unchanged by a flip) is within |w| tol of the current one, as for any signature within
  // tol of it: a fixed generic w makes the candidates few.
  std::vector<double> w(width);
  double wn = 0.0;
  for (int q = 0; q < width; q++)
  {
    w[q] = std::sqrt(q + 2.0);
    wn += w[q] * w[q];
  }
  const double window = std::sqrt(wn) * tol;
  auto key = [&](const double *v)
  {
    double k = 0.0;
    for (int q = 0; q < width; q++)
    {
      k += w[q] * v[q];
    }
    return std::abs(k);
  };
  std::vector<double> ks(n);
  std::vector<int> order(n);
  for (int s = 0; s < n; s++)
  {
    ks[s] = key(row(saved, s));
    order[s] = s;
  }
  std::sort(order.begin(), order.end(), [&](int x, int y) { return ks[x] < ks[y]; });
  SignatureMap map;
  map.rows.assign(n, {});
  std::vector<char> used(n, 0), matched(n, 0);
  for (int g = 0; g < n; g++)
  {
    const double *c = row(current, g);
    const double kc = key(c);
    int best = -1;
    double bd = tol * tol, bs = 1.0;
    auto it = std::lower_bound(order.begin(), order.end(), kc - window,
                               [&](int s, double x) { return ks[s] < x; });
    for (; it != order.end() && ks[*it] <= kc + window; ++it)
    {
      const double *v = row(saved, *it);
      double dp = 0.0, dm = 0.0;
      for (int q = 0; q < width; q++)
      {
        dp += (c[q] - v[q]) * (c[q] - v[q]);
        dm += (c[q] + v[q]) * (c[q] + v[q]);
      }
      if (!nedelec)
      {
        dm = dp;  // H1 DOFs are not flipped
      }
      if (!used[*it] && std::min(dp, dm) <= bd)
      {
        best = *it;
        bd = std::min(dp, dm);
        bs = (dp <= dm) ? 1.0 : -1.0;
      }
    }
    if (best >= 0)
    {
      used[best] = matched[g] = 1;
      map.rows[g] = {{best, bs}};
    }
  }

  // The others with all DOFs at their point (for second-order Nédélec face DOFs on
  // tetrahedra, the two tangents at the face center, one of which may have matched): the
  // block M_g with current = M_g saved, M_g = C S^T (S S^T)^-1 for the signature blocks.
  bool ok = true;
  double worst = 0.0;  // largest block residual
  if (std::find(matched.begin(), matched.end(), 0) != matched.end())
  {
    const auto pc = SignaturePoints(current, width), ps = SignaturePoints(saved, width);
    // The DOFs at a point within the tolerance, searched by their first coordinate.
    auto by_x = [n](const std::vector<std::array<double, 3>> &pts)
    {
      std::vector<int> o(n);
      for (int g = 0; g < n; g++)
      {
        o[g] = g;
      }
      std::sort(o.begin(), o.end(), [&](int x, int y) { return pts[x][0] < pts[y][0]; });
      return o;
    };
    const std::vector<int> oc = by_x(pc), os = by_x(ps);
    auto at = [tol](const std::vector<std::array<double, 3>> &pts,
                    const std::vector<int> &o, const std::array<double, 3> &q)
    {
      std::vector<int> found;
      auto it = std::lower_bound(o.begin(), o.end(), q[0] - tol,
                                 [&](int g, double x) { return pts[g][0] < x; });
      for (; it != o.end() && pts[*it][0] <= q[0] + tol; ++it)
      {
        const auto &p = pts[*it];
        if (std::hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) <= tol)
        {
          found.push_back(*it);
        }
      }
      return found;
    };
    std::vector<char> done(n, 0);
    for (int g0 = 0; g0 < n && ok; g0++)
    {
      if (matched[g0] || done[g0])
      {
        continue;
      }
      if (width == 13 && row(current, g0)[12] != 0.0)
      {
        ok = false;  // a version of a split DOF matches one to one
        break;
      }
      const std::vector<int> C = at(pc, oc, pc[g0]), Sv = at(ps, os, pc[g0]);
      const int k = static_cast<int>(C.size());
      if (k != static_cast<int>(Sv.size()))
      {
        ok = false;
        break;
      }
      mfem::DenseMatrix Cm(k, width), Sm(k, width), SSt(k, k), CSt(k, k), Mg(k, k),
          R(k, width);
      for (int i = 0; i < k; i++)
      {
        for (int q = 0; q < width; q++)
        {
          Cm(i, q) = row(current, C[i])[q];
          Sm(i, q) = row(saved, Sv[i])[q];
        }
      }
      mfem::MultABt(Sm, Sm, SSt);
      mfem::MultABt(Cm, Sm, CSt);
      SSt.Invert();
      mfem::Mult(CSt, SSt, Mg);
      mfem::Mult(Mg, Sm, R);
      R -= Cm;
      worst = std::max(worst, R.MaxMaxNorm());
      ok = (R.MaxMaxNorm() <= tol);  // also false for a NaN residual
      if (!ok)
      {
        break;  // M_g can be singular: do not invert it
      }
      // Row i of M_g^-T: column i of M_g^-1.
      mfem::DenseMatrix Mi(Mg);
      Mi.Invert();
      for (int i = 0; i < k; i++)
      {
        map.rows[C[i]].clear();
        for (int j = 0; j < k; j++)
        {
          map.rows[C[i]].push_back({Sv[j], Mi(j, i)});
        }
        done[C[i]] = 1;
      }
    }
  }
  // A bijection: every current DOF mapped, every saved DOF used.
  std::fill(used.begin(), used.end(), 0);
  for (const auto &r : map.rows)
  {
    ok = ok && !r.empty();
    for (const auto &[s, c] : r)
    {
      used[s] = 1;
    }
  }
  ok = ok && std::find(used.begin(), used.end(), 0) == used.end();
  const auto unmatched = std::count(matched.begin(), matched.end(), 0);
  MFEM_VERIFY(ok, "Online interface DOFs do not match the saved model ("
                      << unmatched << " of " << n
                      << " not matched one to one, largest block residual " << worst
                      << "): the interface (Gamma) must be identical!");
  return map;
}

std::vector<char> BoundaryTrueDofs(mfem::ParFiniteElementSpace &fespace,
                                   const std::vector<int> &bdr_attributes)
{
  // Mark the local DOFs of the boundary elements, then the true DOFs they reach (|P|^T).
  const auto &mesh = *fespace.GetParMesh();
  mfem::Vector lmark(fespace.GetVSize()), tmark(fespace.GetTrueVSize());
  lmark = 0.0;
  mfem::Array<int> dofs;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (std::find(bdr_attributes.begin(), bdr_attributes.end(), mesh.GetBdrAttribute(be)) !=
        bdr_attributes.end())
    {
      fespace.GetBdrElementDofs(be, dofs);
      for (int d : dofs)
      {
        lmark(d >= 0 ? d : -1 - d) = 1.0;
      }
    }
  }
  fespace.Dof_TrueDof_Matrix()->AbsMultTranspose(1.0, lmark, 0.0, tmark);
  const double *t = tmark.HostRead();
  std::vector<char> mark(fespace.GetTrueVSize(), 0);
  for (int i = 0; i < fespace.GetTrueVSize(); i++)
  {
    mark[i] = t[i] > 0.0;
  }
  return mark;
}

std::unique_ptr<mfem::HypreParMatrix> EliminateToIdentity(const mfem::HypreParMatrix &A,
                                                          const mfem::Array<int> &dofs)
{
  mfem::HypreParMatrix B(A);
  B.EliminateBC(dofs, mfem::Operator::DiagonalPolicy::DIAG_ZERO);
  mfem::SparseMatrix diag(A.Height());
  for (int d : dofs)
  {
    diag.Set(d, d, 1.0);
  }
  diag.Finalize();
  mfem::HypreParMatrix I(A.GetComm(), A.GetGlobalNumRows(), A.GetRowStarts(), &diag);
  return std::unique_ptr<mfem::HypreParMatrix>(mfem::ParAdd(&B, &I));
}

std::string RegionMeshPath(const std::string &model_path)
{
  std::filesystem::path p(model_path);
  return p.replace_filename(p.stem().string() + "_region.mesh").string();
}

namespace
{

// A serial mesh of the region, written without the boundary elements duplicated by the
// cracking of interior boundaries (the online run cracks them again).
class SerialRegionMesh : public mfem::Mesh
{
public:
  explicit SerialRegionMesh(mfem::Mesh &&mesh) : mfem::Mesh(std::move(mesh)) {}

  // Shift the attributes above shift back and drop the duplicated boundary elements (for
  // writing only: the boundary tables are not updated). With a plain attribute (nonzero),
  // its boundary elements inside the mesh are dropped, and those the submesh generated
  // (with the next attribute) take it.
  void PrepareForPrint(int shift, int plain)
  {
    for (int e = 0; e < GetNE(); e++)
    {
      if (GetAttribute(e) > shift)
      {
        SetAttribute(e, GetAttribute(e) - shift);
      }
    }
    std::set<std::vector<int>> seen;  // attribute, then the sorted vertices
    int nbe = 0;
    for (int be = 0; be < NumOfBdrElements; be++)
    {
      if (plain > 0 && boundary[be]->GetAttribute() == plain &&
          GetFaceInformation(GetBdrElementFaceIndex(be)).IsInterior())
      {
        FreeElement(boundary[be]);
        continue;
      }
      if (plain > 0 && boundary[be]->GetAttribute() == plain + 1)
      {
        boundary[be]->SetAttribute(plain);
      }
      mfem::Array<int> v;
      boundary[be]->GetVertices(v);
      std::vector<int> key(v.begin(), v.end());
      std::sort(key.begin(), key.end());
      key.insert(key.begin(), boundary[be]->GetAttribute());
      if (!seen.insert(key).second)
      {
        FreeElement(boundary[be]);
        continue;
      }
      boundary[nbe++] = boundary[be];
    }
    NumOfBdrElements = nbe;
    boundary.SetSize(nbe);
    SetAttributes();
  }
};

}  // namespace

std::vector<char> WriteRegionMesh(const std::string &path,
                                  mfem::ParFiniteElementSpace &fespace,
                                  const std::vector<int> &region_attributes,
                                  const std::vector<int> &environment_attributes,
                                  const std::vector<char> &keep, double length_scale,
                                  const std::vector<char> &plain_bdr)
{
  // On a copy of the mesh, the environment elements touching the kept DOFs take shifted
  // attributes for the domain submesh, and get theirs back in the written mesh.
  const auto &parent = *fespace.GetParMesh();
  MFEM_VERIFY(!parent.Nonconforming(),
              "The mesh of the region is written from a conforming mesh only!");
  mfem::Vector tmark(fespace.GetTrueVSize()), lmark(fespace.GetVSize());
  for (int i = 0; i < fespace.GetTrueVSize(); i++)
  {
    tmark(i) = keep[i] ? 1.0 : 0.0;
  }
  fespace.Dof_TrueDof_Matrix()->Mult(tmark, lmark);
  const double *lm = lmark.HostRead();
  mfem::ParMesh mesh(parent);
  const int shift = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
  auto in = [](const std::vector<int> &s, int a)
  { return std::find(s.begin(), s.end(), a) != s.end(); };
  std::vector<char> written(mesh.GetNE(), 0);
  mfem::Array<int> dofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (!in(environment_attributes, parent.GetAttribute(e)))
    {
      continue;
    }
    fespace.GetElementDofs(e, dofs);
    if (std::any_of(dofs.begin(), dofs.end(),
                    [&](int d) { return lm[d >= 0 ? d : -1 - d] != 0.0; }))
    {
      mesh.SetAttribute(e, parent.GetAttribute(e) + shift);
      written[e] = 1;
    }
  }
  // The plain boundary elements take the attribute of the boundary faces the submesh
  // generates (the one after the largest).
  int plain = std::any_of(plain_bdr.begin(), plain_bdr.end(), [](char c) { return c; });
  Mpi::GlobalMax(1, &plain, parent.GetComm());
  if (plain)
  {
    plain = (parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0) + 1;
    for (int be = 0; be < static_cast<int>(plain_bdr.size()); be++)
    {
      if (plain_bdr[be])
      {
        mesh.SetBdrAttribute(be, plain);
      }
    }
  }
  mesh.SetAttributes();
  mfem::Array<int> attrs;
  for (int a : region_attributes)
  {
    attrs.Append(a);
  }
  for (int a : environment_attributes)
  {
    attrs.Append(a + shift);
  }

  // The submesh is extracted from the whole mesh on rank 0 (gathering a parallel submesh
  // can duplicate entities shared between ranks that hold few of its elements).
  MPI_Comm comm = parent.GetComm();
  auto write = [&](mfem::Mesh &smesh)
  {
    WriteRootOutputFile(path, comm,
                        [&](std::ostream &out)
                        {
                          out << std::scientific;
                          out.precision(std::numeric_limits<double>::max_digits10);
                          mesh::DimensionalizeMesh(smesh, length_scale);
                          smesh.Mesh::Print(out);
                        });
  };
  mfem::Mesh smesh = mesh.GetSerialMesh(0);
  SerialRegionMesh sub(Mpi::Root(comm)
                           ? mfem::Mesh(mfem::SubMesh::CreateFromDomain(smesh, attrs))
                           : mfem::Mesh());
  if (Mpi::Root(comm))
  {
    sub.PrepareForPrint(shift, plain);
  }
  write(sub);
  return written;
}

namespace
{

using json = nlohmann::json;

// Add an element of volume (or area) v and center c to the moments of AttributeMoments.
void AddMoments(std::array<double, 6> &m, double v, const mfem::Vector &c)
{
  m[0] += 1.0;
  m[1] += v;
  m[2] += v * c.Norml2();
  for (int d = 0; d < c.Size() && d < 3; d++)
  {
    m[3 + d] += v * c(d);
  }
}

// The moments of each of n attributes, summed over the processes.
void SumMoments(std::vector<std::array<double, 6>> &m, MPI_Comm comm)
{
  if (!m.empty())
  {
    Mpi::GlobalSum(static_cast<int>(6 * m.size()), m.front().data(), comm);
  }
}

// Call f(b) for each boundary element b next to an element with one of the given domain
// attributes (on either side, for an interior boundary), or with within, whose neighbors
// all have one of them.
template <typename F>
void ForEachBoundaryElementNextTo(mfem::ParMesh &mesh, const std::vector<int> &attributes,
                                  bool within, F &&f)
{
  mesh.ExchangeFaceNbrData();
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  auto in = [&](const mfem::ElementTransformation *T)
  { return T && std::ranges::find(attributes, T->Attribute) != attributes.end(); };
  for (int b = 0; b < mesh.GetNBE(); b++)
  {
    BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(b, mesh, FET, T1, T2);
    if (within ? (in(FET.Elem1) && (!FET.Elem2 || in(FET.Elem2)))
               : (in(FET.Elem1) || in(FET.Elem2)))
    {
      f(b);
    }
  }
}

// An entry of the configuration cut to the given sorted attributes (see
// EnvironmentConfigHash), or for an array its entries listing one of them, sorted, and for
// a group of entries those of its values; null if none. Without the excitations, the
// indices are left out.
json CutToAttributes(const json &node, const std::vector<int> &attrs, bool excitations)
{
  if (node.is_array())
  {
    std::vector<json> entries;
    for (const auto &x : node)
    {
      if (json e = CutToAttributes(x, attrs, excitations); !e.is_null())
      {
        entries.push_back(std::move(e));
      }
    }
    std::ranges::sort(entries, {}, [](const json &e) { return e.dump(); });
    return entries.empty() ? json() : json(entries);
  }
  if (!node.is_object())
  {
    return json();
  }
  const bool is_entry = std::ranges::any_of(node.items(),
                                            [](const auto &item)
                                            {
                                              return item.value().is_array() &&
                                                     (item.key().ends_with("Attributes") ||
                                                      item.key() == "Elements");
                                            });
  if (!is_entry)
  {
    json group = json::object();
    for (const auto &[key, value] : node.items())
    {
      if (json cut = CutToAttributes(value, attrs, excitations); !cut.is_null())
      {
        group[key] = std::move(cut);
      }
    }
    return group.empty() ? json() : group;
  }
  json entry = node;
  if (!excitations)
  {
    entry.erase("Index");
  }
  bool listed = false;
  for (auto &[key, value] : entry.items())
  {
    if (key.ends_with("Attributes") && value.is_array())
    {
      std::vector<int> cut;
      for (const auto &a : value)
      {
        if (a.is_number_integer() && std::ranges::binary_search(attrs, a.get<int>()))
        {
          cut.push_back(a.get<int>());
        }
      }
      std::ranges::sort(cut);
      cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
      listed = listed || !cut.empty();
      value = cut;
    }
    else if (key == "Elements" && value.is_array())
    {
      json elements = CutToAttributes(value, attrs, excitations);
      listed = listed || !elements.is_null();
      value = elements.is_null() ? json::array() : elements;
    }
  }
  return listed ? entry : json();
}

// The entries of a section of the configuration, but postprocessing (and without the
// excitations, the sources), cut to the given attributes.
json SectionOn(const json &config, const char *section, std::vector<int> attrs,
               bool excitations)
{
  std::ranges::sort(attrs);
  json entries = json::object();
  const auto it = config.find(section);
  if (it == config.end() || !it->is_object())
  {
    return entries;
  }
  for (const auto &[key, value] : it->items())
  {
    if (key == "Postprocessing" ||
        (!excitations && (key == "SurfaceCurrent" || key == "FluxLoop")))
    {
      continue;
    }
    if (json cut = CutToAttributes(value, attrs, excitations); !cut.is_null())
    {
      entries[key] = std::move(cut);
    }
  }
  return entries;
}

}  // namespace

std::vector<std::array<double, 6>> AttributeMoments(mfem::ParMesh &mesh,
                                                    const std::vector<int> &attributes,
                                                    const std::vector<char> &mark)
{
  std::vector<std::array<double, 6>> m(attributes.size(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const auto it = std::ranges::find(attributes, mesh.GetAttribute(e));
    if ((!mark.empty() && !mark[e]) || it == attributes.end())
    {
      continue;
    }
    mesh.GetElementCenter(e, c);
    AddMoments(m[it - attributes.begin()], mesh.GetElementVolume(e), c);
  }
  SumMoments(m, mesh.GetComm());
  return m;
}

std::vector<int> BoundaryAttributesNextTo(mfem::ParMesh &mesh,
                                          const std::vector<int> &attributes)
{
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::vector<int> next(bdr_attr_max, 0);
  ForEachBoundaryElementNextTo(mesh, attributes, false,
                               [&](int b) { next[mesh.GetBdrAttribute(b) - 1] = 1; });
  Mpi::GlobalMax(bdr_attr_max, next.data(), mesh.GetComm());
  std::vector<int> bdr_attributes;
  for (int a = 0; a < bdr_attr_max; a++)
  {
    if (next[a])
    {
      bdr_attributes.push_back(a + 1);
    }
  }
  return bdr_attributes;
}

std::vector<std::array<double, 6>>
BoundaryAttributeMoments(mfem::ParMesh &mesh, const std::vector<int> &bdr_attributes,
                         const std::vector<int> &attributes)
{
  std::vector<std::array<double, 6>> m(bdr_attributes.size(),
                                       {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  mfem::Vector c;
  ForEachBoundaryElementNextTo(
      mesh, attributes, true,
      [&](int b)
      {
        const auto it = std::ranges::find(bdr_attributes, mesh.GetBdrAttribute(b));
        if (it == bdr_attributes.end())
        {
          return;
        }
        const auto geom = mesh.GetBdrElementGeometry(b);
        mfem::ElementTransformation &T = *mesh.GetBdrElementTransformation(b);
        double area = 0.0;
        for (const auto &ip : mfem::IntRules.Get(geom, 2 * T.Order() + 2))
        {
          T.SetIntPoint(&ip);
          area += ip.weight * T.Weight();
        }
        T.Transform(mfem::Geometries.GetCenter(geom), c);
        AddMoments(m[it - bdr_attributes.begin()], area, c);
      });
  SumMoments(m, mesh.GetComm());
  return m;
}

std::uint64_t EnvironmentConfigHash(const nlohmann::json &config,
                                    const std::vector<int> &attributes,
                                    const std::vector<int> &bdr_attributes,
                                    bool excitations)
{
  auto setting = [&config](const char *section, const char *key)
  {
    const auto it = config.find(section);
    return (it != config.end() && it->is_object()) ? it->value(key, json()) : json();
  };
  const json described = {
      {"Problem", {{"Type", setting("Problem", "Type")}}},
      {"Model",
       {{"CrackInternalBoundaryElements",
         setting("Model", "CrackInternalBoundaryElements")}}},
      {"Solver",
       {{"Order", setting("Solver", "Order")},
        {"QuadratureOrderJacobian", setting("Solver", "QuadratureOrderJacobian")},
        {"QuadratureOrderExtra", setting("Solver", "QuadratureOrderExtra")}}},
      {"Domains", SectionOn(config, "Domains", attributes, excitations)},
      {"Boundaries", SectionOn(config, "Boundaries", bdr_attributes, excitations)}};
  // FNV-1a of the description.
  std::uint64_t hash = 14695981039346656037ULL;
  for (const unsigned char c : described.dump())
  {
    hash = (hash ^ c) * 1099511628211ULL;
  }
  return hash;
}

EnvironmentRecord EnvironmentRecord::Make(mfem::ParMesh &mesh,
                                          const std::vector<int> &attributes,
                                          const nlohmann::json &config,
                                          const std::vector<char> &written)
{
  EnvironmentRecord r;
  r.attributes = attributes;
  std::ranges::sort(r.attributes);
  r.attributes.erase(std::unique(r.attributes.begin(), r.attributes.end()),
                     r.attributes.end());
  r.bdr_attributes = BoundaryAttributesNextTo(mesh, r.attributes);
  r.moments = AttributeMoments(mesh, r.attributes);
  r.region_moments = AttributeMoments(mesh, r.attributes, written);
  r.bdr_moments = BoundaryAttributeMoments(mesh, r.bdr_attributes, r.attributes);
  r.config_hash = EnvironmentConfigHash(config, r.attributes, r.bdr_attributes, false);
  r.excitations_hash = EnvironmentConfigHash(config, r.attributes, r.bdr_attributes, true);
  return r;
}

std::vector<char> EnvironmentRecord::Serialize() const
{
  std::vector<char> bytes;
  auto put = [&bytes]<typename T>(const T *x, std::size_t n)
  {
    const auto *p = reinterpret_cast<const char *>(x);
    bytes.insert(bytes.end(), p, p + sizeof(T) * n);
  };
  const int n[2] = {static_cast<int>(attributes.size()),
                    static_cast<int>(bdr_attributes.size())};
  put(n, 2);
  put(attributes.data(), attributes.size());
  put(bdr_attributes.data(), bdr_attributes.size());
  put(moments.data(), moments.size());
  put(region_moments.data(), region_moments.size());
  put(bdr_moments.data(), bdr_moments.size());
  put(&config_hash, 1);
  put(&excitations_hash, 1);
  return bytes;
}

std::optional<EnvironmentRecord>
EnvironmentRecord::Deserialize(const std::vector<char> &bytes)
{
  std::size_t at = 0;
  auto get = [&]<typename T>(T *x, std::size_t n)
  {
    if (bytes.size() - at < sizeof(T) * n)
    {
      return false;
    }
    std::memcpy(x, bytes.data() + at, sizeof(T) * n);
    at += sizeof(T) * n;
    return true;
  };
  int n[2];
  if (bytes.size() < sizeof(n) || !get(n, 2) || n[0] < 0 || n[1] < 0)
  {
    return std::nullopt;
  }
  EnvironmentRecord r;
  r.attributes.resize(n[0]);
  r.bdr_attributes.resize(n[1]);
  r.moments.resize(n[0]);
  r.region_moments.resize(n[0]);
  r.bdr_moments.resize(n[1]);
  if (!get(r.attributes.data(), r.attributes.size()) ||
      !get(r.bdr_attributes.data(), r.bdr_attributes.size()) ||
      !get(r.moments.data(), r.moments.size()) ||
      !get(r.region_moments.data(), r.region_moments.size()) ||
      !get(r.bdr_moments.data(), r.bdr_moments.size()) || !get(&r.config_hash, 1) ||
      !get(&r.excitations_hash, 1) || at != bytes.size())
  {
    return std::nullopt;
  }
  return r;
}

bool EnvironmentRecord::Check(mfem::ParMesh &mesh, const nlohmann::json &config,
                              const std::vector<int> &environment_attributes, double tol,
                              bool region_mesh) const
{
  auto match = [tol](const std::vector<std::array<double, 6>> &a,
                     const std::vector<std::array<double, 6>> &b)
  {
    return std::ranges::equal(a, b, [tol](const auto &x, const auto &y)
                              { return MomentsMatch(x, y, tol); });
  };
  const auto current = AttributeMoments(mesh, attributes);
  const bool full = match(current, moments), region = match(current, region_moments);
  if (region_mesh)
  {
    MFEM_VERIFY(region,
                (full ? "\"RegionMesh\" is set, but the mesh holds the whole "
                        "environment of the saved substructuring model: unset it on "
                        "the full mesh!"
                      : "The mesh is not the mesh of the region written with the "
                        "saved substructuring model, nor a re-meshing of its region!"));
  }
  else
  {
    MFEM_VERIFY(full, (region ? "The mesh is the mesh of the region written with the saved "
                                "substructuring model: set \"RegionMesh\": true!"
                              : kEnvironmentMeshDiffers));
    MFEM_VERIFY(
        BoundaryAttributesNextTo(mesh, attributes) == bdr_attributes &&
            match(BoundaryAttributeMoments(mesh, bdr_attributes, attributes), bdr_moments),
        kEnvironmentMeshDiffers);
  }
  std::vector<int> env = environment_attributes;
  std::ranges::sort(env);
  env.erase(std::unique(env.begin(), env.end()), env.end());
  MFEM_VERIFY((region_mesh ? std::ranges::includes(attributes, env) : env == attributes) &&
                  EnvironmentConfigHash(config, attributes, bdr_attributes, false) ==
                      config_hash,
              kEnvironmentConfigDiffers);
  const bool excitations =
      EnvironmentConfigHash(config, attributes, bdr_attributes, true) == excitations_hash;
  MFEM_VERIFY(!region_mesh || excitations, kEnvironmentExcitationsDiffer);
  return excitations;
}

bool MomentsMatch(const std::array<double, 6> &a, const std::array<double, 6> &b,
                  double tol)
{
  bool match = (a[0] == b[0]);
  for (int q = 1; q < 6; q++)
  {
    const int s = std::min(q, 2);  // the scale of the moment
    match =
        match && std::abs(a[q] - b[q]) <= tol * std::max(std::abs(a[s]), std::abs(b[s]));
  }
  return match;
}

std::vector<char> RegionMeshInterface(mfem::ParFiniteElementSpace &parent_fespace,
                                      const std::vector<char> &candidate,
                                      const std::vector<int> &parent,
                                      const std::vector<double> &saved_signatures,
                                      int saved_width)
{
  const auto &mesh = *parent_fespace.GetParMesh();
  MPI_Comm comm = parent_fespace.GetComm();
  const int nt = static_cast<int>(candidate.size());
  auto parent_of = [&](int t) { return parent.empty() ? t : parent[t]; };
  const std::vector<char> bdr =
      BoundaryTrueDofs(parent_fespace, std::vector<int>(mesh.bdr_attributes.begin(),
                                                        mesh.bdr_attributes.end()));
  // The candidates' parent DOFs on the boundary, numbered over the processes.
  std::vector<int> index(parent_fespace.GetTrueVSize(), -1);
  int nc = 0;
  for (int t = 0; t < nt; t++)
  {
    const int t0 = parent_of(t);
    if (bdr[t0] && candidate[t] && index[t0] < 0)
    {
      index[t0] = nc++;
    }
  }
  int off = 0;
  MPI_Exscan(&nc, &off, 1, MPI_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  Mpi::GlobalSum(1, &nc, comm);
  for (auto &g : index)
  {
    g += (g >= 0) ? off : 0;
  }
  const auto match =
      MatchSignaturePoints(TrueDofSignatures(parent_fespace, index, nc),
                           SignatureWidth(parent_fespace), saved_signatures, saved_width);
  std::vector<char> cut(nt, 0);
  for (int t = 0; t < nt; t++)
  {
    const int t0 = parent_of(t);
    cut[t] = candidate[t] && index[t0] >= 0 && match[index[t0]];
  }
  return cut;
}

namespace
{

// The scalar polynomial k of FingerprintField at the point (X, Y, Z).
double FingerprintPolynomial(int k, double X, double Y, double Z)
{
  return (k == 0) ? X : (k == 1) ? Y + Z * Z : X * Y * Z;
}

}  // namespace

mfem::Vector FingerprintField(const mfem::ParFiniteElementSpace &fespace, int k,
                              const std::array<double, 3> &center)
{
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(fespace);
  mfem::ParGridFunction gf(&fes);
  if (fes.FEColl()->GetRangeType(fes.GetMesh()->Dimension()) == mfem::FiniteElement::VECTOR)
  {
    mfem::VectorFunctionCoefficient c(
        3,
        [k, &center](const mfem::Vector &x, mfem::Vector &v)
        {
          const double X = x(0) - center[0], Y = x(1) - center[1], Z = x(2) - center[2];
          const double f[3][3] = {{Z, X, Y}, {Y * Y, Z * Z, X * X}, {X * Y, Y * Z, Z * X}};
          v.SetSize(3);
          for (int d = 0; d < 3; d++)
          {
            v(d) = f[k][d];
          }
        });
    gf.ProjectCoefficient(c);
  }
  else
  {
    mfem::FunctionCoefficient c(
        [k, &center](const mfem::Vector &x)
        {
          return FingerprintPolynomial(k, x(0) - center[0], x(1) - center[1],
                                       (x.Size() > 2 ? x(2) : 0.0) - center[2]);
        });
    gf.ProjectCoefficient(c);
  }
  mfem::Vector r(fes.GetTrueVSize());
  gf.GetTrueDofs(r);
  return r;
}

}  // namespace palace
