// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <utility>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

namespace
{

// Coefficient given per (domain or boundary) element.
class ElementCoefficient : public mfem::Coefficient
{
  const std::vector<double> &value;

public:
  explicit ElementCoefficient(const std::vector<double> &value) : value(value) {}
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &) override
  {
    return value[T.ElementNo];
  }
};

// MFEM integrator with the quadrature order of Palace's own operators, so the condensed
// discretization matches the native one also on curved meshes.
template <typename Integrator>
class NativeQuadrature : public Integrator
{
public:
  using Integrator::Integrator;
  void AssembleElementMatrix(const mfem::FiniteElement &el, mfem::ElementTransformation &T,
                             mfem::DenseMatrix &M) override
  {
    this->IntRule =
        &mfem::IntRules.Get(el.GetGeomType(), fem::DefaultIntegrationOrder::Get(T));
    Integrator::AssembleElementMatrix(el, T, M);
  }
};

// Piecewise-constant (per element attribute) matrix coefficient for anisotropic materials.
// Attributes absent from the map contribute a zero tensor (e.g. environment attributes when
// assembling the region operator).
class PWMatrixCoefficient : public mfem::MatrixCoefficient
{
  const std::map<int, mfem::DenseMatrix> &mats;

public:
  PWMatrixCoefficient(int dim, const std::map<int, mfem::DenseMatrix> &m)
    : mfem::MatrixCoefficient(dim), mats(m)
  {
  }
  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &) override
  {
    auto it = mats.find(T.Attribute);
    if (it != mats.end())
    {
      K = it->second;
    }
    else
    {
      K.SetSize(width);
      K = 0.0;
    }
  }
};

}  // namespace

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Aenv()
{
  if (!A_env)
  {
    A_env = AssembleParent(env_eps, magnetostatic, true, &sheet_env);
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
    K_env = AssembleParent(env_eps, false);
  }
  return *K_env;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Denv()
{
  if (!D_env)
  {
    D_env = AssembleParent(env_eps, true, false);
  }
  return *D_env;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::Msheet()
{
  if (!M_sheet)
  {
    M_sheet = AssembleParent({}, false, false, &sheet_all);
  }
  return *M_sheet;
}

const mfem::HypreParMatrix &SubstructuringSolver::Impl::MsheetEnv()
{
  if (!M_sheet_env)
  {
    M_sheet_env = AssembleParent({}, false, false, &sheet_env);
  }
  return *M_sheet_env;
}

void SubstructuringSolver::Impl::ForEachElement(
    const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass, bool with_curl,
    const std::vector<double> *sheet,
    const std::function<bool(const mfem::Array<int> &)> &use,
    const std::function<void(const mfem::Array<int> &, const mfem::DenseMatrix &)> &f) const
{
  PWMatrixCoefficient coef(parent.Dimension(), coef_by_attr);
  mfem::Vector mass(parent.attributes.Size() ? parent.attributes.Max() : 1);
  mass = mag_eps;
  mfem::PWConstCoefficient mcoef(mass);
  NativeQuadrature<mfem::CurlCurlIntegrator> curlcurl(coef);
  NativeQuadrature<mfem::VectorFEMassIntegrator> regularization(mcoef);
  NativeQuadrature<mfem::DiffusionIntegrator> diffusion(coef);
  mfem::Array<int> vdofs;
  mfem::DofTransformation doftrans;
  mfem::DenseMatrix Ae, Me;
  const bool volume = with_curl || (magnetostatic && with_mass);
  for (int e = 0; volume && e < parent.GetNE(); e++)
  {
    if (!coef_by_attr.contains(parent.GetAttribute(e)))
    {
      continue;
    }
    SolutionSpace().GetElementDofs(e, vdofs, doftrans);
    if (!use(vdofs))
    {
      continue;
    }
    const mfem::FiniteElement &fe = *parent_fes.GetFE(e);
    mfem::ElementTransformation &T = *parent_fes.GetElementTransformation(e);
    if (!magnetostatic)
    {
      diffusion.AssembleElementMatrix(fe, T, Ae);
    }
    else if (with_curl)
    {
      curlcurl.AssembleElementMatrix(fe, T, Ae);
      if (with_mass)
      {
        regularization.AssembleElementMatrix(fe, T, Me);
        Ae += Me;
      }
    }
    else
    {
      regularization.AssembleElementMatrix(fe, T, Ae);
    }
    doftrans.TransformDual(Ae);
    f(vdofs, Ae);
  }
  if (sheet && magnetostatic && has_sheets)
  {
    ElementCoefficient sheet_coef(*sheet);
    NativeQuadrature<mfem::VectorFEMassIntegrator> integ(sheet_coef);
    for (int be = 0; be < parent.GetNBE(); be++)
    {
      if ((*sheet)[be] == 0.0)
      {
        continue;
      }
      SolutionSpace().GetBdrElementDofs(be, vdofs, doftrans);
      if (!use(vdofs))
      {
        continue;
      }
      integ.AssembleElementMatrix(*parent_fes.GetBE(be),
                                  *parent_fes.GetBdrElementTransformation(be), Ae);
      doftrans.TransformDual(Ae);
      f(vdofs, Ae);
    }
  }
}

std::vector<mfem::Vector>
SubstructuringSolver::Impl::EnvApply(const std::vector<const mfem::Vector *> &X,
                                     bool volume, bool sheets) const
{
  const int n = static_cast<int>(X.size());
  const mfem::Operator &P = *SolutionSpace().GetProlongationMatrix();
  std::vector<mfem::Vector> XL(n, mfem::Vector(SolutionSpace().GetVSize())),
      YL(n, mfem::Vector(SolutionSpace().GetVSize()));
  for (int k = 0; k < n; k++)
  {
    P.Mult(*X[k], XL[k]);
    YL[k] = 0.0;
  }
  std::vector<char> touched(n, 0);
  std::vector<const double *> xl(n);
  for (int k = 0; k < n; k++)
  {
    xl[k] = XL[k].HostRead();
  }
  auto use = [&](const mfem::Array<int> &vdofs)
  {
    bool any = false;
    for (int k = 0; k < n; k++)
    {
      touched[k] = 0;
      for (int d : vdofs)
      {
        if (xl[k][d >= 0 ? d : -1 - d] != 0.0)
        {
          touched[k] = 1;
          break;
        }
      }
      any = any || touched[k];
    }
    return any;
  };
  mfem::Vector xe, ye;
  ForEachElement(env_eps, volume, volume, sheets ? &sheet_env.coef : nullptr, use,
                 [&](const mfem::Array<int> &vdofs, const mfem::DenseMatrix &Ae)
                 {
                   for (int k = 0; k < n; k++)
                   {
                     if (touched[k])
                     {
                       XL[k].GetSubVector(vdofs, xe);
                       ye.SetSize(xe.Size());
                       Ae.Mult(xe, ye);
                       YL[k].AddElementVector(vdofs, ye);
                     }
                   }
                 });
  std::vector<mfem::Vector> Y(n, mfem::Vector(nt));
  for (int k = 0; k < n; k++)
  {
    P.MultTranspose(YL[k], Y[k]);
    if (sheets && sheet_env.two_sided)
    {
      sheet_env.two_sided->Mult(1.0, *X[k], 1.0, Y[k]);
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

std::array<std::unique_ptr<mfem::HypreParMatrix>, 2>
SubstructuringSolver::Impl::SetUpBrokenSpace()
{
  std::vector<config::SuperconductorData> films;
  for (const auto &sc : iodata.boundaries.superconductor)
  {
    if (sc.two_sided)
    {
      auto &film = films.emplace_back(sc);
      std::erase_if(film.attributes,
                    [this](int a) { return parent.bdr_attributes.Find(a) < 0; });
    }
  }
  if (films.empty())
  {
    return {};
  }
  MaterialOperator mat_op(iodata, mesh);
  SuperconductorSheetOperator sheet_op(films, iodata.boundaries.cracked_attributes,
                                       iodata.units, mat_op, parent, false);
  const auto attrs = sheet_op.GetTwoSidedAttributes();
  if (attrs.empty())
  {
    return {};
  }
  const auto &sides = mesh.GetCrackSides(attrs);
  broken = std::make_unique<FiniteElementSpace>(parent_space, sides);
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
  const auto &sub = *iodata.solver.substructuring;
  return {sheet_op.AssembleBrokenSheetMatrix(parent_space, *broken, sides,
                                             &sub.region_attributes),
          sheet_op.AssembleBrokenSheetMatrix(parent_space, *broken, sides,
                                             &sub.environment_attributes)};
}

void SubstructuringSolver::Impl::SetUpSheets(
    const mfem::Array<int> &ra,
    std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> two_sided)
{
  const int nbe = parent.GetNBE();
  sheet_region.coef.assign(nbe, 0.0);
  sheet_env.coef.assign(nbe, 0.0);
  sheet_all.coef.assign(nbe, 0.0);
  std::map<int, double> coef;
  std::set<int> films;
  for (const auto &sc : iodata.boundaries.superconductor)
  {
    if (sc.two_sided)
    {
      films.insert(sc.attributes.begin(), sc.attributes.end());
      continue;
    }
    const double Ls = (sc.Ls > 0.0) ? sc.Ls
                                    : SuperconductorSheetOperator::KineticSheetInductance(
                                          sc.lambda_L, sc.thickness);
    for (int a : sc.attributes)
    {
      coef[a] = 1.0 / (Ls * (iodata.boundaries.cracked_attributes.contains(a) ? 2.0 : 1.0));
    }
  }
  // Counts of the sheet faces, of the film faces between the region and the environment,
  // and of the film faces in the region.
  std::array<int, 3> count = {0, 0, 0};
  parent.ExchangeFaceNbrData();
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  for (int be = 0; be < nbe; be++)
  {
    const int attr = parent.GetBdrAttribute(be);
    auto it = coef.find(attr);
    const bool film = films.contains(attr);
    if (it == coef.end() && !film)
    {
      continue;
    }
    BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(be, parent, FET, T1,
                                                                     T2);
    const bool region1 = ra.Find(FET.Elem1->Attribute) >= 0,
               region2 = FET.Elem2 && ra.Find(FET.Elem2->Attribute) >= 0;
    if (film)
    {
      count[1] += (region1 != region2);
      count[2] += region1;
      continue;
    }
    (region1 || region2 ? sheet_region : sheet_env).coef[be] = it->second;
    sheet_all.coef[be] = it->second;  // all sheet faces, for Msheet()
    count[0]++;
  }
  Mpi::GlobalSum(3, count.data(), parent.GetComm());
  has_sheets = (count[0] > 0) || two_sided[0];
  MFEM_VERIFY(count[1] == 0, "Substructuring does not support a \"TwoSided\" "
                             "superconductor film between the region and the environment!");
  MFEM_VERIFY(count[2] == 0 || iodata.model.refinement.max_it == 0,
              "Adaptive refinement with substructuring does not support \"TwoSided\" "
              "superconductor films in the region!");
  if (!two_sided[0])
  {
    return;
  }
  sheet_all.two_sided.reset(mfem::Add(1.0, *two_sided[0], 1.0, *two_sided[1]));
  sheet_region.two_sided = std::move(two_sided[0]);
  sheet_env.two_sided = std::move(two_sided[1]);
}

mfem::Vector SubstructuringSolver::Impl::EnvSource(const Source &src,
                                                   std::array<double, 3> &fp)
{
  mfem::Vector b(nt);
  b = 0.0;
  if (src.HasA())
  {
    b = std::move(EnvApply({src.a}, false, true)[0]);
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

std::unique_ptr<mfem::HypreParMatrix> SubstructuringSolver::Impl::AssembleParent(
    const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass, bool with_curl,
    const Sheets *sheet)
{
  if (broken)
  {
    // On the broken space: the terms element by element, in the DOFs of the side of each
    // element, and the two-port terms of the films.
    auto A_l = std::make_unique<mfem::SparseMatrix>(broken->GetVSize());
    ForEachElement(
        coef_by_attr, with_mass, with_curl, sheet ? &sheet->coef : nullptr,
        [](const mfem::Array<int> &) { return true; },
        [&](const mfem::Array<int> &vdofs, const mfem::DenseMatrix &Ae)
        { A_l->AddSubMatrix(vdofs, vdofs, Ae, 0); });
    A_l->Finalize();
    auto A =
        linalg::RAP(*linalg::BlockDiagonalParMatrix(
                        parent_fes.GetComm(), broken->GlobalVSize(), broken->GlobalVSize(),
                        broken->GetDofOffsets(), broken->GetDofOffsets(), std::move(A_l)),
                    *broken->Dof_TrueDof_Matrix());
    if (sheet && sheet->two_sided)
    {
      A.reset(mfem::Add(1.0, *A, 1.0, *sheet->two_sided));
    }
    // The transpose product Pᵀ (A P) can store zeros on one side only of the pattern of a
    // broken space. The sparse direct solvers are given a symmetric pattern (STRUMPACK
    // orders the graph of the pattern as it is): symmetrize.
    std::unique_ptr<mfem::HypreParMatrix> At(A->Transpose());
    return std::unique_ptr<mfem::HypreParMatrix>(mfem::Add(0.5, *A, 0.5, *At));
  }
  PWMatrixCoefficient coef(parent.Dimension(), coef_by_attr);
  mfem::ParBilinearForm a(&parent_fes);
  std::unique_ptr<ElementCoefficient> sheet_coef;
  if (magnetostatic && sheet && has_sheets)
  {
    // London kinetic sheet term (1/L_ksq) A_t . v_t on the given boundary elements.
    sheet_coef = std::make_unique<ElementCoefficient>(sheet->coef);
    a.AddBoundaryIntegrator(
        new NativeQuadrature<mfem::VectorFEMassIntegrator>(*sheet_coef));
  }
  // Domain integrators only on the attributes of coef_by_attr (zero elsewhere).
  const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
  mfem::Array<int> marker(max_attr);
  marker = 0;
  for (const auto &[attr, m] : coef_by_attr)
  {
    marker[attr - 1] = 1;
  }
  mfem::Vector mass(max_attr);
  mass = mag_eps;
  mfem::PWConstCoefficient mcoef(mass);
  if (magnetostatic)
  {
    if (with_curl)
    {
      a.AddDomainIntegrator(new NativeQuadrature<mfem::CurlCurlIntegrator>(coef), marker);
    }
    if (with_mass)
    {
      a.AddDomainIntegrator(new NativeQuadrature<mfem::VectorFEMassIntegrator>(mcoef),
                            marker);
    }
  }
  else
  {
    a.AddDomainIntegrator(new NativeQuadrature<mfem::DiffusionIntegrator>(coef), marker);
  }
  a.Assemble();
  a.Finalize();
  return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
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
  RowLayout(1, cnt, disp);
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
  RowLayout(1, cnt, disp);
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

std::unique_ptr<mfem::HypreParMatrix> SubstructuringSolver::Impl::AssembleEnvSubmesh()
{
  MFEM_VERIFY(!has_sheets, "Superconductor sheets in substructuring need a direct "
                           "environment factorization: a sparse direct solver and an "
                           "environment within \"Solver.Substructuring.DirectMaxSize\"!");
  PWMatrixCoefficient coef(parent.Dimension(), env_eps);
  mfem::Vector mass(env_submesh->attributes.Size() ? env_submesh->attributes.Max() : 1);
  mass = 0.0;
  for (const auto &[attr, t] : env_eps)
  {
    mass(attr - 1) = mag_eps;
  }
  mfem::PWConstCoefficient mcoef(mass);
  mfem::ParBilinearForm a(env_sfes.get());
  if (magnetostatic)
  {
    a.AddDomainIntegrator(new NativeQuadrature<mfem::CurlCurlIntegrator>(coef));
    a.AddDomainIntegrator(new NativeQuadrature<mfem::VectorFEMassIntegrator>(mcoef));
  }
  else
  {
    a.AddDomainIntegrator(new NativeQuadrature<mfem::DiffusionIntegrator>(coef));
  }
  a.Assemble();
  a.Finalize();
  return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
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
  S_rows.assign(static_cast<std::size_t>(nloc) * nG_global, 0.0);
}

}  // namespace palace
