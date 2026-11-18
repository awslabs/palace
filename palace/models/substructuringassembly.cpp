// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <cmath>
#include <optional>
#include <set>
#include <utility>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

namespace
{

// The sum of the given matrices (the null ones left out).
std::unique_ptr<mfem::HypreParMatrix>
Sum(std::initializer_list<const mfem::HypreParMatrix *> terms)
{
  std::unique_ptr<mfem::HypreParMatrix> A;
  for (const auto *B : terms)
  {
    if (B)
    {
      A.reset(A ? mfem::Add(1.0, *A, 1.0, *B) : new mfem::HypreParMatrix(*B));
    }
  }
  return A;
}

}  // namespace

std::unique_ptr<ParOperator>
SubstructuringSolver::Impl::GetTerm(const std::vector<int> *attrs, Term term,
                                    bool assemble) const
{
  if (!magnetostatic)
  {
    MFEM_ASSERT(term == Term::VOLUME, "Electrostatic operators have the volume term only!");
    std::optional<LaplaceOperator::AssemblyRestriction> restriction;
    if (attrs)
    {
      restriction.emplace(*laplace_op, *attrs);
    }
    return laplace_op->GetUnconstrainedStiffnessMatrix(assemble);
  }
  if (term == Term::MASS)
  {
    // The regularization mag_eps (1/mu) M on the domains (the exactly zero entries of the
    // other elements are not stored), weighted as the curl-curl operator. The coefficient
    // is built for every attribute and then restricted, so that it is a matrix one also on
    // a process without an element of the domains.
    MFEM_ASSERT(attrs, "The regularization is assembled on the region or the environment!");
    const auto &mat_op = curlcurl_op->GetMaterialOp();
    MaterialPropertyCoefficient eps(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvPermeability(), mag_eps);
    eps.RestrictCoefficient(mat_op.GetCeedAttributes(*attrs));
    BilinearForm m(SolutionSpace());
    m.AddDomainIntegrator<VectorFEMassIntegrator>(eps);
    return std::make_unique<ParOperator>(
        assemble ? std::unique_ptr<Operator>(m.FullAssemble(true)) : m.PartialAssemble(),
        SolutionSpace());
  }
  if (term == Term::SHEETS && !has_sheets)
  {
    return nullptr;
  }
  std::optional<CurlCurlOperator::AssemblyRestriction> restriction;
  if (attrs)
  {
    restriction.emplace(*curlcurl_op, *attrs);
  }
  return curlcurl_op->GetUnconstrainedStiffnessMatrix(term == Term::VOLUME,
                                                      term == Term::SHEETS, assemble);
}

std::unique_ptr<mfem::HypreParMatrix>
SubstructuringSolver::Impl::AssembleTerm(const std::vector<int> *attrs, Term term) const
{
  auto op = GetTerm(attrs, term, true);
  if (!op)
  {
    return nullptr;
  }
  auto A = op->StealParallelAssemble();
  if (broken)
  {
    // The transpose product Pᵀ (A P) can store zeros on one side only of the pattern of a
    // broken space. The sparse direct solvers are given a symmetric pattern (STRUMPACK
    // orders the graph of the pattern as it is): symmetrize.
    std::unique_ptr<mfem::HypreParMatrix> At(A->Transpose());
    A.reset(mfem::Add(0.5, *A, 0.5, *At));
  }
  return A;
}

void SubstructuringSolver::Impl::AssembleRegion()
{
  const auto *region = &iodata.solver.substructuring->region_attributes;
  if (!magnetostatic)
  {
    A_region = AssembleTerm(region, Term::VOLUME);
    K_region_e = A_region.get();
    return;
  }
  // The energy operator is the volume term (see SourceEnergyMatrix for the London sheets).
  A_region_energy = AssembleTerm(region, Term::VOLUME);
  K_region_e = A_region_energy.get();
  M_sheet_region = AssembleTerm(region, Term::SHEETS);
  A_region = Sum({A_region_energy.get(), M_sheet_region.get(),
                  AssembleTerm(region, Term::MASS).get()});
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Aenv()
{
  if (!A_env)
  {
    A_env = magnetostatic ? Sum({&Kenv(), has_sheets ? &MsheetEnv() : nullptr, &Denv()})
                          : AssembleTerm(&env_attributes, Term::VOLUME);
  }
  return *A_env;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Kenv()
{
  if (!magnetostatic)
  {
    return Aenv();
  }
  if (!K_env)
  {
    K_env = AssembleTerm(&env_attributes, Term::VOLUME);
  }
  return *K_env;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Denv()
{
  if (!D_env)
  {
    D_env = AssembleTerm(&env_attributes, Term::MASS);
  }
  return *D_env;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Msheet()
{
  if (!M_sheet)
  {
    M_sheet = AssembleTerm(nullptr, Term::SHEETS);
  }
  return *M_sheet;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::MsheetEnv()
{
  if (!M_sheet_env)
  {
    M_sheet_env = AssembleTerm(&env_attributes, Term::SHEETS);
  }
  return *M_sheet_env;
}

bool SubstructuringSolver::Impl::EnvSheets()
{
  if (!has_sheets)
  {
    return false;
  }
  // hypre computes the norm on the host, the entries of a GPU build being on the device.
  const mfem::HypreParMatrix &M = MsheetEnv();
  M.HostRead();
  const double norm = M.FNorm();
  M.HypreRead();
  return norm > 0.0;
}

std::vector<mfem::Vector>
SubstructuringSolver::Impl::EnvApply(const std::vector<const mfem::Vector *> &X,
                                     bool volume) const
{
  // The assembled operator if there is one (offline), else the operator applied
  // matrix-free, set up on first use.
  const mfem::HypreParMatrix *M = volume ? A_env.get() : M_sheet_env.get();
  auto &op = volume ? env_op : env_sheet_op;
  if (!M && !op)
  {
    const auto *env = &env_attributes;
    auto sum = std::make_unique<SumOperator>(nt, nt);
    for (auto term : volume ? std::vector<Term>{Term::VOLUME, Term::SHEETS, Term::MASS}
                            : std::vector<Term>{Term::SHEETS})
    {
      if (!magnetostatic && term != Term::VOLUME)
      {
        continue;
      }
      if (auto t = GetTerm(env, term, false))
      {
        sum->AddOperator(std::move(t));
      }
    }
    op = std::move(sum);
  }
  std::vector<mfem::Vector> Y(X.size());
  Vector x(nt);
  x.UseDevice(true);
  for (std::size_t k = 0; k < X.size(); k++)
  {
    x = *X[k];
    Y[k].SetSize(nt);
    Y[k].UseDevice(true);
    if (M)
    {
      M->Mult(x, Y[k]);
    }
    else
    {
      op->Mult(x, Y[k]);
    }
  }
  return Y;
}

const FiniteElementSpace &SubstructuringSolver::Impl::SolutionSpace() const
{
  return broken ? *broken : parent_space;
}

const mfem::Vector &SubstructuringSolver::Impl::FirstVersion(const mfem::Vector &x,
                                                             mfem::Vector &y) const
{
  if (!broken)
  {
    return x;
  }
  y.SetSize(first_version->Height());
  first_version->Mult(x, y);
  return y;
}

void SubstructuringSolver::Impl::SetUpBrokenSpace()
{
  if (!curlcurl_op->IsBroken())
  {
    return;
  }
  broken = &curlcurl_op->GetSolutionNDSpace();
  first_version = broken->BuildVersionMap(FiniteElementSpace::VersionMap::FIRST);
  all_versions = broken->BuildVersionMap(FiniteElementSpace::VersionMap::ALL);
  version_parent.resize(broken->GetTrueVSize());
  for (int t = 0; t < parent_fes.GetTrueVSize(); t++)
  {
    for (int v = 0; v < broken->GetNumVersions(t); v++)
    {
      version_parent[broken->GetVersionTrueDof(t, v)] = t;
    }
  }
}

void SubstructuringSolver::Impl::CheckSheets()
{
  // The faces of the one-sided sheets, of the film faces between the region and the
  // environment, and of the film faces in the region. The native operators assign a sheet
  // face to the side of its neighbor with the lower attribute (see
  // Mesh::GetCeedBdrAttributesByNeighbor), which for a face on Gamma may be either side.
  std::set<int> sheets, films;
  for (const auto &sc : iodata.boundaries.superconductor)
  {
    (sc.two_sided ? films : sheets).insert(sc.attributes.begin(), sc.attributes.end());
  }
  std::array<int, 3> count = {0, 0, 0};
  parent.ExchangeFaceNbrData();
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  gamma_sheet_faces.clear();
  gamma_sheet_sides = {};
  constexpr double k[3][3] = {{1.3, 2.9, 0.7}, {2.3, 0.5, 3.1}, {0.9, 1.7, 2.6}};
  mfem::Array<int> vertices;
  for (int be = 0; be < parent.GetNBE(); be++)
  {
    const int attr = parent.GetBdrAttribute(be);
    const bool sheet = sheets.contains(attr), film = films.contains(attr);
    if (!sheet && !film)
    {
      continue;
    }
    BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(be, parent, FET, T1,
                                                                     T2);
    const bool region1 = ra_arr.Find(FET.Elem1->Attribute) >= 0,
               region2 = FET.Elem2 && ra_arr.Find(FET.Elem2->Attribute) >= 0;
    if (film)
    {
      count[1] += (region1 != region2);
      count[2] += region1;
      continue;
    }
    if (FET.Elem2 && region1 != region2)
    {
      gamma_sheet_faces.push_back(be);
      if ((FET.Elem1->Attribute < FET.Elem2->Attribute) ? region1 : region2)
      {
        std::array<double, 3> c = {0.0, 0.0, 0.0};
        parent.GetBdrElementVertices(be, vertices);
        for (int v : vertices)
        {
          for (int d = 0; d < 3; d++)
          {
            c[d] += parent.GetVertex(v)[d] / vertices.Size();
          }
        }
        gamma_sheet_sides[0] += 1.0;
        for (int j = 0; j < 3; j++)
        {
          gamma_sheet_sides[j + 1] +=
              std::cos(k[j][0] * c[0] + k[j][1] * c[1] + k[j][2] * c[2]);
        }
      }
    }
    count[0]++;
  }
  Mpi::GlobalSum(3, count.data(), parent.GetComm());
  Mpi::GlobalSum(4, gamma_sheet_sides.data(), parent.GetComm());
  has_sheets = (count[0] > 0) || broken;
  MFEM_VERIFY(count[1] == 0, "Substructuring does not support a \"TwoSided\" "
                             "superconductor film between the region and the environment!");
  MFEM_VERIFY(count[2] == 0 || iodata.model.refinement.max_it == 0,
              "Adaptive refinement with substructuring does not support \"TwoSided\" "
              "superconductor films in the region!");
}

mfem::Vector SubstructuringSolver::Impl::EnvSource(const Source &src,
                                                   std::array<double, 3> &fp)
{
  mfem::Vector b(nt);
  b = 0.0;
  if (src.HasA())
  {
    b = std::move(EnvApply({src.a}, false)[0]);
  }
  double *bb = b.HostReadWrite();
  if (src.HasJ())
  {
    const double *J = src.J->HostRead();
    for (int i = 0; i < nt; i++)
    {
      if (is_env_int[i])
      {
        bb[i] += J[i];
      }
    }
  }
  for (int d = 0; d < dbc_tdofs.Size(); d++)
  {
    bb[dbc_tdofs[d]] = 0.0;
  }
  double loc[3] = {b.Normlinf(), CachedFingerprintField(0) * b,
                   CachedFingerprintField(1) * b},
         glob[3];
  MPI_Allreduce(loc, glob, 1, MPI_DOUBLE, MPI_MAX, parent_fes.GetComm());
  MPI_Allreduce(loc + 1, glob + 1, 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
  fp = {glob[0] > 0.0 ? 1.0 : 0.0, glob[1], glob[2]};
  return b;
}

mfem::Vector SubstructuringSolver::Impl::RegionSource(const Source &src) const
{
  mfem::Vector b(nt);
  b = 0.0;
  if (src.HasA() && M_sheet_region)
  {
    M_sheet_region->Mult(*src.a, b);
  }
  if (src.HasJ())
  {
    const double *J = src.J->HostRead();
    double *bb = b.HostReadWrite();
    for (int i = 0; i < nt; i++)
    {
      if (is_region_free[i])
      {
        bb[i] += J[i];
      }
    }
  }
  return b;
}

std::vector<HYPRE_BigInt> SubstructuringSolver::Impl::InterfaceTrueDofs() const
{
  const HYPRE_BigInt tstart = SolutionSpace().GetTrueDofOffsets()[0];
  std::vector<HYPRE_BigInt> mine(gamma_nloc), g2t(nG_global);
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      mine[gamma_global[i] - gamma_off] = tstart + i;
    }
  }
  std::vector<int> cnt, disp;
  RowLayout(cnt, disp);
  MPI_Allgatherv(mine.data(), gamma_nloc, HYPRE_MPI_BIG_INT, g2t.data(), cnt.data(),
                 disp.data(), HYPRE_MPI_BIG_INT, parent_fes.GetComm());
  return g2t;
}

bool SubstructuringSolver::Impl::InterfaceSplit() const
{
  int split = 0;
  for (int t = 0; broken && t < nt && !split; t++)
  {
    split = is_gamma[t] && broken->GetNumVersions(ParentTrueDof(t)) > 1;
  }
  Mpi::GlobalMax(1, &split, parent.GetComm());
  return split > 0;
}

std::vector<double> SubstructuringSolver::Impl::InterfaceSignatures(int width) const
{
  if (!broken)
  {
    return TrueDofSignatures(parent_fes, gamma_global, nG_global);
  }
  // The row of an interface DOF is that of its parent DOF, stored at the interface index
  // of the first of its versions on Gamma.
  MPI_Comm comm = parent_fes.GetComm();
  const int w0 = SignatureWidth(parent_fes);
  std::vector<int> index(parent_fes.GetTrueVSize(), -1), mine(gamma_nloc), row(nG_global);
  for (int t = 0; t < nt; t++)
  {
    if (is_gamma[t])
    {
      int &i = index[ParentTrueDof(t)];
      i = (i < 0) ? gamma_global[t] : i;
      mine[gamma_global[t] - gamma_off] = i;
    }
  }
  std::vector<int> cnt, disp;
  RowLayout(cnt, disp);
  MPI_Allgatherv(mine.data(), gamma_nloc, MPI_INT, row.data(), cnt.data(), disp.data(),
                 MPI_INT, comm);
  const std::vector<double> sig0 = TrueDofSignatures(parent_fes, index, nG_global);
  std::vector<double> sig(static_cast<std::size_t>(nG_global) * width, 0.0);
  for (int g = 0; g < nG_global; g++)
  {
    std::copy_n(&sig0[static_cast<std::size_t>(row[g]) * w0], w0,
                &sig[static_cast<std::size_t>(g) * width]);
  }
  if (width == w0)
  {
    return sig;
  }

  // The sums Σ_f (c_f - m) × n_f of the versions: on the L-vector, from the region
  // element of each boundary face of the region (with no region element across it, also
  // on another process), then on the true DOFs.
  auto bounds_region = [this](int e, int f)
  {
    const auto info = parent.GetFaceInformation(f);
    if (info.IsBoundary())
    {
      return true;
    }
    const int other =
        info.IsShared()
            ? parent.face_nbr_elements[info.element[1].index]->GetAttribute()
            : parent.GetAttribute((info.element[0].index == e) ? info.element[1].index
                                                               : info.element[0].index);
    return ra_arr.Find(other) < 0;
  };
  std::array<mfem::Vector, 3> WL, W;
  for (auto &v : WL)
  {
    v.SetSize(broken->GetVSize());
    v = 0.0;
  }
  mfem::Array<int> faces, orient, fedges, fo, eedges, eo, vdofs, ents;
  mfem::DofTransformation doftrans;
  mfem::Vector ce(3), c(3), n(3), m(3);
  for (int e = 0; e < parent.GetNE(); e++)
  {
    if (ra_arr.Find(parent.GetAttribute(e)) < 0)
    {
      continue;
    }
    const auto geom = parent.GetElementGeometry(e);
    const int nv = mfem::Geometry::NumVerts[geom];
    parent.GetElementFaces(e, faces, orient);
    bool first = true;
    for (int f : faces)
    {
      if (!bounds_region(e, f))
      {
        continue;
      }
      if (first)
      {
        parent.GetElementTransformation(e)->Transform(mfem::Geometries.GetCenter(geom), ce);
        parent.GetElementEdges(e, eedges, eo);
        broken->GetElementDofs(e, vdofs, doftrans);
        fem::GetElementDofEntities(parent_fes, e, ents);
        first = false;
      }
      mfem::ElementTransformation &T = *parent.GetFaceTransformation(f);
      const mfem::IntegrationPoint &ip = mfem::Geometries.GetCenter(T.GetGeometryType());
      T.SetIntPoint(&ip);
      T.Transform(ip, c);
      mfem::CalcOrtho(T.Jacobian(), n);
      n *= ((c(0) - ce(0)) * n(0) + (c(1) - ce(1)) * n(1) + (c(2) - ce(2)) * n(2) < 0.0)
               ? -1.0 / n.Norml2()
               : 1.0 / n.Norml2();
      parent.GetFaceEdges(f, fedges, fo);
      for (int j = 0; j < vdofs.Size(); j++)
      {
        const int b = ents[j] - nv;
        if (b < 0 || b >= eedges.Size() || fedges.Find(eedges[b]) < 0)
        {
          continue;
        }
        mfem::ElementTransformation &Te = *parent.GetEdgeTransformation(eedges[b]);
        Te.Transform(mfem::Geometries.GetCenter(mfem::Geometry::SEGMENT), m);
        const double d[3] = {c(0) - m(0), c(1) - m(1), c(2) - m(2)};
        const int l = (vdofs[j] >= 0) ? vdofs[j] : -1 - vdofs[j];
        WL[0](l) += d[1] * n(2) - d[2] * n(1);
        WL[1](l) += d[2] * n(0) - d[0] * n(2);
        WL[2](l) += d[0] * n(1) - d[1] * n(0);
      }
    }
  }
  const mfem::HypreParMatrix &P = *broken->Dof_TrueDof_Matrix();
  std::array<const double *, 3> w;
  for (int k = 0; k < 3; k++)
  {
    W[k].SetSize(nt);
    P.AbsMultTranspose(1.0, WL[k], 0.0, W[k]);
    w[k] = W[k].HostRead();
  }
  std::vector<double> q(nG_global, 0.0);
  int missing = 0;
  for (int t = 0; t < nt; t++)
  {
    if (is_gamma[t] && broken->GetNumVersions(ParentTrueDof(t)) > 1)
    {
      const std::size_t g = gamma_global[t];
      for (int b = 0; b < 3; b++)
      {
        q[g] += sig[g * width + b] * w[b][t];
      }
      missing += (w[0][t] == 0.0 && w[1][t] == 0.0 && w[2][t] == 0.0);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, q.data(), nG_global, MPI_DOUBLE, MPI_SUM, comm);
  Mpi::GlobalSum(1, &missing, comm);
  MFEM_VERIFY(missing == 0, "Substructuring found "
                                << missing
                                << " interface DOF(s) on \"TwoSided\" superconductor "
                                   "films without a boundary face of the region on their "
                                   "side!");
  for (int g = 0; g < nG_global; g++)
  {
    sig[static_cast<std::size_t>(g) * width + w0] = q[g];
  }
  return sig;
}

mfem::Vector SubstructuringSolver::Impl::TerminalMode(int idx) const
{
  mfem::Vector x(nt);
  x = 0.0;
  for (int d : terminal_tdofs.at(idx))
  {
    x(d) = 1.0;
  }
  return x;
}

std::unique_ptr<mfem::HypreParMatrix>
SubstructuringSolver::Impl::AssembleSelection(const std::vector<char> &mark) const
{
  MPI_Comm comm = parent_fes.GetComm();
  int nloc = 0;
  for (int i = 0; i < nt; i++)
  {
    nloc += mark[i];
  }
  int off = 0, nglob = 0;
  MPI_Exscan(&nloc, &off, 1, MPI_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, comm);
  std::vector<int> I(nt + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> V;
  for (int i = 0; i < nt; i++)
  {
    if (mark[i])
    {
      J.push_back(off + static_cast<int>(J.size()));
      V.push_back(1.0);
    }
    I[i + 1] = static_cast<int>(J.size());
  }
  // Column partitioning in the layout HYPRE expects for this build.
  std::vector<HYPRE_BigInt> cols;
  if (HYPRE_AssumedPartitionCheck())
  {
    cols = {off, off + nloc, nglob};
  }
  else
  {
    const int nranks = Mpi::Size(comm);
    std::vector<int> cnt(nranks);
    MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    cols.assign(nranks + 1, 0);
    for (int r = 0; r < nranks; r++)
    {
      cols[r + 1] = cols[r] + cnt[r];
    }
  }
  return std::make_unique<mfem::HypreParMatrix>(
      comm, nt, SolutionSpace().GlobalTrueVSize(), static_cast<HYPRE_BigInt>(nglob),
      I.data(), J.data(), V.data(), SolutionSpace().GetTrueDofOffsets(), cols.data());
}

std::unique_ptr<mfem::HypreParMatrix>
SubstructuringSolver::Impl::AssembleDenseInterface() const
{
  const int nG = nG_global, nloc = gamma_nloc;
  const std::vector<HYPRE_BigInt> g2t = InterfaceTrueDofs();
  std::vector<int> I(nt + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> V;
  J.reserve(static_cast<std::size_t>(nloc) * nG);
  V.reserve(static_cast<std::size_t>(nloc) * nG);
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      const double *row =
          &S_rows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * nG];
      for (int j = 0; j < nG; j++)
      {
        J.push_back(g2t[j]);
        V.push_back(row[j]);
      }
    }
    I[i + 1] = static_cast<int>(J.size());
  }
  const HYPRE_BigInt glob = SolutionSpace().GlobalTrueVSize();
  return std::make_unique<mfem::HypreParMatrix>(
      parent_fes.GetComm(), nt, glob, glob, I.data(), J.data(), V.data(),
      SolutionSpace().GetTrueDofOffsets(), SolutionSpace().GetTrueDofOffsets());
}

void SubstructuringSolver::Impl::EnumerateInterface()
{
  MPI_Comm comm = parent_fes.GetComm();
  int nloc = 0;
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      nloc++;
    }
  }
  int off = 0;
  MPI_Exscan(&nloc, &off, 1, MPI_INT, MPI_SUM, comm);
  off = Mpi::Root(comm) ? 0 : off;
  nG_global = 0;
  MPI_Allreduce(&nloc, &nG_global, 1, MPI_INT, MPI_SUM, comm);
  gamma_global.assign(nt, -1);
  {
    int c = off;
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        gamma_global[i] = c++;
      }
    }
  }

  gamma_off = off;
  gamma_nloc = nloc;
  // The condensation fills S_rows: free that of a previous interface first.
  std::vector<double>().swap(S_rows);
}

}  // namespace palace
