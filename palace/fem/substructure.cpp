// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructure.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/outputdir.hpp"

namespace palace
{

int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker,
                          mfem::Array<int> &interface_marker)
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
  mfem::Vector reg_l(parent_fespace.GetVSize()), env_l(parent_fespace.GetVSize());
  reg_l = 0.0;
  env_l = 0.0;
  mfem::Array<int> dofs;
  auto *mesh = parent_fespace.GetParMesh();
  for (int e = 0; e < mesh->GetNE(); e++)
  {
    const int attr = mesh->GetAttribute(e);
    const bool is_reg = in(region_attrs, attr), is_env = in(environment_attrs, attr);
    if (!is_reg && !is_env)
    {
      continue;
    }
    parent_fespace.GetElementDofs(e, dofs);
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
  const mfem::HypreParMatrix *P = parent_fespace.Dof_TrueDof_Matrix();
  const int nt = parent_fespace.GetTrueVSize();
  mfem::Vector reg_t(nt), env_t(nt);
  P->AbsMultTranspose(1.0, reg_l, 0.0, reg_t);
  P->AbsMultTranspose(1.0, env_l, 0.0, env_t);

  region_marker.SetSize(nt);
  env_marker.SetSize(nt);
  interface_marker.SetSize(nt);
  int local_iface = 0;
  for (int i = 0; i < nt; i++)
  {
    region_marker[i] = (reg_t(i) > 0.5) ? 1 : 0;
    env_marker[i] = (env_t(i) > 0.5) ? 1 : 0;
    interface_marker[i] = (region_marker[i] && env_marker[i]) ? 1 : 0;
    local_iface += interface_marker[i];
  }
  int global_iface = 0;
  MPI_Allreduce(&local_iface, &global_iface, 1, MPI_INT, MPI_SUM, parent_fespace.GetComm());
  return global_iface;
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

std::vector<double> TrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<int> &index, int n)
{
  const int w = SignatureWidth(fespace), nt = fespace.GetTrueVSize();
  std::vector<double> loc(static_cast<std::size_t>(n) * w, 0.0),
      glob(static_cast<std::size_t>(n) * w, 0.0);
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(fespace);
  mfem::ParGridFunction gf(&fes);
  mfem::Vector td(nt);
  auto stamp = [&](int slot)
  {
    gf.GetTrueDofs(td);
    for (int i = 0; i < nt; i++)
    {
      if (index[i] >= 0)
      {
        loc[static_cast<std::size_t>(index[i]) * w + slot] = td(i);
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
  MPI_Allreduce(loc.data(), glob.data(), n * w, MPI_DOUBLE, MPI_SUM, fespace.GetComm());
  return glob;
}

namespace
{

// The point p of each DOF from its signature: the coordinates of an H1 DOF, and for a
// point-tangent DOF u(p).t with the signature (t_b, p_a t_b), p = Q t / |t|^2 with
// Q_ab = p_a t_b.
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

std::vector<char> MatchSignaturePoints(const std::vector<double> &candidates,
                                       const std::vector<double> &saved, int width)
{
  MFEM_VERIFY(width == 3 || width == 12, "Unknown DOF signature width " << width << "!");
  const double tol = SignatureTolerance(saved);
  const auto pc = SignaturePoints(candidates, width), ps = SignaturePoints(saved, width);
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

SignatureMap MatchSignatureBasis(const std::vector<double> &current,
                                 const std::vector<double> &saved, int width)
{
  MFEM_VERIFY(width == 12, "MatchSignatureBasis expects Nédélec signatures!");
  MFEM_VERIFY(current.size() == saved.size(),
              "Saved and current DOF sets of different sizes cannot be matched!");
  const int n = static_cast<int>(current.size() / width);
  const double tol = SignatureTolerance(saved);
  auto row = [width](const std::vector<double> &sig, int g) { return &sig[g * width]; };

  // Signed one-to-one matches.
  SignatureMap map;
  map.rows.assign(n, {});
  std::vector<char> used(n, 0), matched(n, 0);
  for (int g = 0; g < n; g++)
  {
    const double *c = row(current, g);
    int best = -1;
    double bd = tol * tol, bs = 1.0;
    for (int s = 0; s < n; s++)
    {
      const double *v = row(saved, s);
      double dp = 0.0, dm = 0.0;
      for (int q = 0; q < width; q++)
      {
        dp += (c[q] - v[q]) * (c[q] - v[q]);
        dm += (c[q] + v[q]) * (c[q] + v[q]);
      }
      if (!used[s] && std::min(dp, dm) <= bd)
      {
        best = s;
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
    auto near = [tol](const std::array<double, 3> &p, const std::array<double, 3> &q)
    { return std::hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) <= tol; };
    std::vector<char> done(n, 0);
    for (int g0 = 0; g0 < n && ok; g0++)
    {
      if (matched[g0] || done[g0])
      {
        continue;
      }
      std::vector<int> C, Sv;
      for (int g = 0; g < n; g++)
      {
        if (near(pc[g], pc[g0]))
        {
          C.push_back(g);
        }
        if (near(ps[g], pc[g0]))
        {
          Sv.push_back(g);
        }
      }
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
      ok = (R.MaxMaxNorm() <= tol);  // also false for a singular block (NaN)
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
  MFEM_VERIFY(ok, "Online interface DOFs do not match the saved model (max mismatch "
                      << worst << "); the interface Gamma must be identical.");
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
  std::vector<char> mark(fespace.GetTrueVSize(), 0);
  for (int i = 0; i < fespace.GetTrueVSize(); i++)
  {
    mark[i] = tmark(i) > 0.0;
  }
  return mark;
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
  // writing only: the boundary tables are not updated).
  void Finalize(int shift)
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

void WriteRegionMesh(const std::string &path, mfem::ParFiniteElementSpace &fespace,
                     const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes,
                     const std::vector<char> &keep, double length_scale)
{
  // On a copy of the mesh, the environment elements touching the kept DOFs take shifted
  // attributes for the domain submesh, and get theirs back in the written mesh.
  const auto &parent = *fespace.GetParMesh();
  mfem::Vector tmark(fespace.GetTrueVSize()), lmark(fespace.GetVSize());
  for (int i = 0; i < fespace.GetTrueVSize(); i++)
  {
    tmark(i) = keep[i] ? 1.0 : 0.0;
  }
  fespace.Dof_TrueDof_Matrix()->Mult(tmark, lmark);
  mfem::ParMesh mesh(parent);
  const int shift = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
  auto in = [](const std::vector<int> &s, int a)
  { return std::find(s.begin(), s.end(), a) != s.end(); };
  mfem::Array<int> dofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (!in(environment_attributes, parent.GetAttribute(e)))
    {
      continue;
    }
    fespace.GetElementDofs(e, dofs);
    if (std::any_of(dofs.begin(), dofs.end(),
                    [&](int d) { return lmark(d >= 0 ? d : -1 - d) != 0.0; }))
    {
      mesh.SetAttribute(e, parent.GetAttribute(e) + shift);
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
  if (mesh.Nonconforming())
  {
    mfem::Array<int> rank0(mesh.GetNE());
    rank0 = 0;
    mesh.Rebalance(rank0);
    auto sub = mfem::ParSubMesh::CreateFromDomain(mesh, attrs);
    for (int e = 0; e < sub.GetNE(); e++)
    {
      if (sub.GetAttribute(e) > shift)
      {
        sub.SetAttribute(e, sub.GetAttribute(e) - shift);
      }
    }
    sub.SetAttributes();
    write(sub);
  }
  else
  {
    mfem::Mesh smesh = mesh.GetSerialMesh(0);
    SerialRegionMesh sub(Mpi::Root(comm)
                             ? mfem::Mesh(mfem::SubMesh::CreateFromDomain(smesh, attrs))
                             : mfem::Mesh());
    if (Mpi::Root(comm))
    {
      sub.Finalize(shift);
    }
    write(sub);
  }
}

}  // namespace palace
