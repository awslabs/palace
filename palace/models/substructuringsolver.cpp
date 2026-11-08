// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/amg.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Preconditioner forwarding to a callable (the region factor or region-submesh multigrid,
// both built up front, so SetOperator is a no-op).
class CallableSolver : public Solver<Operator>
{
public:
  CallableSolver(int h, std::function<void(const mfem::Vector &, mfem::Vector &)> f)
    : Solver<Operator>(false), apply(std::move(f))
  {
    height = width = h;
  }
  void SetOperator(const Operator &) override {}
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { apply(x, y); }

private:
  std::function<void(const mfem::Vector &, mfem::Vector &)> apply;
};

}  // namespace

SubstructuringSolver::Impl::Impl(const IoData &iodata, Mesh &mesh)
  : iodata(iodata), mesh(mesh), parent(mesh.Get()),
    magnetostatic(iodata.problem.type == ProblemType::MAGNETOSTATIC),
    fec(magnetostatic
            ? std::unique_ptr<mfem::FiniteElementCollection>(
                  new mfem::ND_FECollection(iodata.solver.order, parent.Dimension()))
            : std::unique_ptr<mfem::FiniteElementCollection>(
                  new mfem::H1_FECollection(iodata.solver.order, parent.Dimension()))),
    parent_space(mesh, fec.get()), parent_fes(parent_space.Get())
{
  MFEM_VERIFY(iodata.solver.substructuring,
              "SubstructuringSolver requires a Solver.Substructuring configuration!");
  const auto &sub = *iodata.solver.substructuring;
  CheckSubstructuringAttributes(parent, sub.region_attributes, sub.environment_attributes);
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> two_sided;
  if (magnetostatic)
  {
    two_sided = SetUpBrokenSpace();
  }
  nt = SolutionSpace().GetTrueVSize();
  for (int a : sub.region_attributes)
  {
    ra_arr.Append(a);
  }
  for (int a : sub.environment_attributes)
  {
    ea_arr.Append(a);
  }

  // Interface / region / environment true-DOF markers.
  mfem::Array<int> rm, em, im;
  MarkInterfaceTrueDofs(SolutionSpace(), ra_arr, ea_arr, rm, em, im);

  // Dirichlet true DOFs: the terminals (per terminal) and the grounded boundaries. The
  // set is fixed across excitations, so the interface/interior partition is too.
  const auto &terminals = iodata.boundaries.terminal;
  mfem::Array<int> dir_mark(nt);
  dir_mark = 0;
  {
    const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
    for (const auto &[idx, term] : terminals)
    {
      mfem::Array<int> ess_bdr(maxb), ess;
      ess_bdr = 0;
      for (int a : term.attributes)
      {
        if (a >= 1 && a <= maxb)
        {
          ess_bdr[a - 1] = 1;
        }
      }
      SolutionSpace().GetEssentialTrueDofs(ess_bdr, ess);
      auto &list = terminal_tdofs[idx];
      for (int i = 0; i < ess.Size(); i++)
      {
        dir_mark[ess[i]] = 1;
        list.push_back(ess[i]);
      }
    }
  }
  // PEC boundaries are fixed at zero, as in the native operators (flux-loop films are
  // London sheets, not Dirichlet boundaries).
  {
    const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
    mfem::Array<int> ess_bdr(maxb), ess;
    ess_bdr = 0;
    auto mark = [&](int a)
    {
      if (a >= 1 && a <= maxb)
      {
        ess_bdr[a - 1] = 1;
      }
    };
    for (int a : iodata.boundaries.pec.attributes)
    {
      mark(a);
    }
    SolutionSpace().GetEssentialTrueDofs(ess_bdr, ess);
    for (int i = 0; i < ess.Size(); i++)
    {
      dir_mark[ess[i]] = 1;
    }
  }
  for (int i = 0; i < nt; i++)
  {
    if (dir_mark[i])
    {
      dbc_tdofs.Append(i);
    }
  }
  region_dbc.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    region_dbc[i] = rm[i] && dir_mark[i];
  }

  // Online, the head of the model and whether the mesh is that of the region.
  std::vector<char> cut(nt, 0);
  double model_eps = 0.0;
  if (sub.mode == SubstructuringMode::ONLINE && !sub.save_model.empty())
  {
    {
      BlockTimer bt(Timer::IO);
      saved = ReadSavedHead(sub.save_model, parent_fes.GetComm());
    }
    const std::string &model = sub.save_model;
    MFEM_VERIFY(saved.version >= 0, "Cannot read the saved substructuring model \""
                                        << model << "\" (run in \"Offline\" mode first)!");
    MFEM_VERIFY(saved.version > 0, "\"" << model
                                        << "\" is not a saved substructuring model: run "
                                           "in \"Offline\" mode to write one!");
    MFEM_VERIFY(saved.version == kModelVersion,
                "The saved substructuring model \""
                    << model
                    << "\" was written by another version of Palace: rerun in "
                       "\"Offline\" mode!");
    MFEM_VERIFY(saved.width > 0, "The saved substructuring model \""
                                     << model
                                     << "\" is incomplete: rerun in \"Offline\" mode!");
    MFEM_VERIFY((saved.sig_type == 1) != magnetostatic,
                "The saved substructuring model \""
                    << model << "\" is "
                    << (magnetostatic ? "electrostatic" : "magnetostatic")
                    << ", unlike this run!");
    MFEM_VERIFY(saved.order == iodata.solver.order,
                "The saved substructuring model \""
                    << model << "\" has \"Order\": " << saved.order << ", unlike this run ("
                    << iodata.solver.order << ")!");
    // The mesh holds the environment of the model, or its part in the mesh of the region.
    const auto moments = ElementMoments(parent, sub.environment_attributes);
    region_only = !MomentsMatch(moments, saved.env_moments, kFingerprintTol);
    MFEM_VERIFY(!region_only ||
                    MomentsMatch(moments, saved.region_moments, kFingerprintTol),
                kEnvironmentDiffers);
    if (region_only)
    {
      model_eps = saved.regularization;
      MFEM_VERIFY(iodata.model.refinement.max_it == 0,
                  "Adaptive refinement with substructuring needs the full mesh, not the "
                  "mesh of the region!");
      std::vector<char> bdr =
          BoundaryTrueDofs(parent_fes, std::vector<int>(parent.bdr_attributes.begin(),
                                                        parent.bdr_attributes.end()));
      // The parent DOFs on the boundary with a version on the region only, at the points
      // of the saved interface DOFs: those versions are interface DOFs (both sides of a
      // two-sided film crossing Gamma).
      auto region_side = [&](int t) { return rm[t] && !em[t] && !dir_mark[t]; };
      std::vector<int> index(parent_fes.GetTrueVSize(), -1);
      int nc = 0;
      for (int t = 0; t < nt; t++)
      {
        const int t0 = ParentTrueDof(t);
        if (bdr[t0] && region_side(t) && index[t0] < 0)
        {
          index[t0] = nc++;
        }
      }
      int off = 0;
      MPI_Exscan(&nc, &off, 1, MPI_INT, MPI_SUM, parent_fes.GetComm());
      off = Mpi::Root(parent_fes.GetComm()) ? 0 : off;
      Mpi::GlobalSum(1, &nc, parent_fes.GetComm());
      for (auto &g : index)
      {
        g += (g >= 0) ? off : 0;
      }
      const auto match =
          MatchSignaturePoints(TrueDofSignatures(parent_fes, index, nc),
                               SignatureWidth(parent_fes), saved.sig, saved.width);
      for (int t = 0; t < nt; t++)
      {
        const int t0 = ParentTrueDof(t);
        cut[t] = region_side(t) && (index[t0] >= 0) && match[index[t0]];
      }
    }
  }

  is_gamma.assign(nt, 0);
  is_env_int.assign(nt, 0);
  is_region_free.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    if (dir_mark[i])
    {
      env_dbc_local += em[i] ? 1 : 0;
      continue;
    }
    const bool r = rm[i], e = em[i];
    if ((r && e) || cut[i])
    {
      is_gamma[i] = 1;
    }
    else if (e)
    {
      is_env_int[i] = 1;
    }
    if (r)
    {
      is_region_free[i] = 1;  // region-free includes the interface
    }
  }
  for (int i = 0; i < nt; i++)
  {
    if (is_region_free[i])
    {
      rf_idx.Append(i);
    }
    if (rm[i])
    {
      rc_idx.Append(i);
    }
    if (is_gamma[i])
    {
      gamma_idx.Append(i);
    }
  }

  // Domain-restricted material coefficients (zero outside the subdomain).
  const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
  const int dim = parent.Dimension();
  // Material tensor from its eigendecomposition, M_ij = sum_k s[k] v[k]_i v[k]_j
  // (inverted for the curl-curl inverse permeability).
  auto tensor = [dim](const config::SymmetricMatrixData<3> &prop, bool invert)
  {
    mfem::DenseMatrix e(dim);
    e = 0.0;
    for (int k = 0; k < 3; k++)
    {
      for (int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          e(i, j) += prop.s[k] * prop.v[k][i] * prop.v[k][j];
        }
      }
    }
    if (invert)
    {
      e.Invert();
    }
    return e;
  };
  for (const auto &mat : iodata.domains.materials)
  {
    for (int a : mat.attributes)
    {
      if (a < 1 || a > max_attr)
      {
        continue;
      }
      if (ra_arr.Find(a) >= 0)
      {
        region_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
      }
      else if (ea_arr.Find(a) >= 0)
      {
        env_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
      }
    }
  }
  if (magnetostatic)
  {
    SetUpSheets(ra_arr, std::move(two_sided));
    // The region solve factors A_region_free + S_E with the dense S_E (see
    // EnsureRegionFactor), which the compression discards.
    MFEM_VERIFY(!has_sheets || sub.interface_offdiag_tol == 0.0,
                "Substructuring with London superconductor sheets needs the dense "
                "interface operator (\"InterfaceOffdiagTol\": 0)!");
    mag_eps = (model_eps > 0.0) ? model_eps
              : has_sheets      ? kMagRegularizationSheets
                                : kMagRegularization;
  }
  A_region = AssembleParent(region_eps, magnetostatic, true, &sheet_region);
  if (magnetostatic)
  {
    // Energy operator: the pure curl-curl (see SourceEnergyMatrix for the London sheets).
    A_region_energy = AssembleParent(region_eps, false);
    K_region_e = A_region_energy.get();
    if (has_sheets)
    {
      M_sheet_region = AssembleParent({}, false, false, &sheet_region);
    }
  }
  else
  {
    K_region_e = A_region.get();
  }
}

void SubstructuringSolver::Impl::EnsureRegionFactor(int n)
{
  if (!region_direct || reg_lu)
  {
    return;
  }
  BlockTimer bt(Timer::KSP_SETUP);
  // The exact factor solves each batch directly (one multi-RHS solve), the other one
  // preconditions CG per excitation, which converges to its tolerance. With London
  // sheets, the exact factor is used whenever it fits.
  region_exact =
      region_dense_eligible &&
      (has_sheets ||
       static_cast<long long>(std::max(n, 1)) * kDenseInterfacePerExcitation >= nG_global);
  long long nrf = rf_idx.Size();
  MPI_Allreduce(MPI_IN_PLACE, &nrf, 1, MPI_LONG_LONG, MPI_SUM, parent_fes.GetComm());
  const bool serial = (nrf <= kSerialRegionFactorDofs);
  reg_lu = std::make_unique<substructuring::BlockedDirectSolver>(
      substructuring::MakeDirectSolver(iodata, parent_fes.GetComm(), serial),
      region_exact
          ? std::min(substructuring::BlockedDirectSolver::kMaxBlock, std::max(n, 1))
          : 1);
  std::unique_ptr<mfem::HypreParMatrix> E = AssembleSelection(is_region_free);
  std::unique_ptr<mfem::HypreParMatrix> P;
  if (region_exact)
  {
    std::unique_ptr<mfem::HypreParMatrix> S_mat = AssembleDenseInterface();
    std::unique_ptr<mfem::HypreParMatrix> A(mfem::ParAdd(A_region_free.get(), S_mat.get()));
    S_mat.reset();
    P.reset(mfem::RAP(A.get(), E.get()));
  }
  else
  {
    P.reset(mfem::RAP(A_region_free.get(), E.get()));
  }
  reg_lu->SetOperator(*P);  // the factor keeps its own copy of P
  reg_r.SetSize(rf_idx.Size());
  reg_z.SetSize(rf_idx.Size());
}

void SubstructuringSolver::Impl::ApplyRegionFactor(const mfem::Vector &r,
                                                   mfem::Vector &z) const
{
  z = r;
  r.GetSubVector(rf_idx, reg_r);
  reg_lu->Mult({&reg_r}, {&reg_z});
  z.SetSubVector(rf_idx, reg_z);
}

void SubstructuringSolver::Impl::SolveRegion(std::vector<mfem::Vector> &x) const
{
  BlockTimer bt(Timer::KSP);
  const int n = static_cast<int>(x.size());
  if (region_exact)
  {
    std::vector<const mfem::Vector *> X(n);
    std::vector<mfem::Vector *> Y(n);
    for (int k = 0; k < n; k++)
    {
      X[k] = &x[k];
      Y[k] = &x[k];
    }
    reg_lu->Mult(X, Y);
    return;
  }
  mfem::Vector b(nt), u(nt);
  b.UseDevice(true);
  u.UseDevice(true);
  for (int k = 0; k < n; k++)
  {
    b = 0.0;
    b.SetSubVector(rf_idx, x[k]);
    u = 0.0;
    region_ksp->Mult(b, u);
    MFEM_VERIFY(region_ksp->GetConverged() || !has_sheets,
                "Substructuring region solve with London sheets did not converge!");
    u.GetSubVector(rf_idx, x[k]);
  }
}

mfem::Vector SubstructuringSolver::Impl::DirichletField(const mfem::Vector &xd) const
{
  mfem::Vector u(nt);
  u = 0.0;
  u.SetSubVector(dbc_tdofs, xd);
  return u;
}

mfem::Vector SubstructuringSolver::Impl::RegionField(const mfem::Vector &x_rf,
                                                     const mfem::Vector *xd) const
{
  mfem::Vector u(nt);
  u = 0.0;
  u.SetSubVector(rf_idx, x_rf);
  if (xd)
  {
    u.SetSubVector(dbc_tdofs, *xd);
  }
  return u;
}

bool SubstructuringSolver::Impl::BuildRegionGmg()
{
  const int order = iodata.solver.order;
  const int dim = parent.Dimension();
  const int mg_levels = iodata.solver.linear.mg_max_levels;
  if (magnetostatic || order <= 1 || mg_levels <= 1)
  {
    return false;
  }

  reg_mesh_vec.clear();
  reg_mesh_vec.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(parent, ra_arr))));
  reg_submesh = dynamic_cast<mfem::ParSubMesh *>(&reg_mesh_vec[0]->Get());
  reg_mesh_vec[0]->RebuildCeedAttributes();
  reg_pgf.SetSpace(&parent_fes);

  // Essential submesh boundary attributes: the region Dirichlet ones (Gamma free).
  mfem::Vector marker(nt);
  marker = 0.0;
  for (int i = 0; i < dbc_tdofs.Size(); i++)
  {
    marker(dbc_tdofs[i]) = 1.0;
  }
  mfem::Array<int> ess_attr = EssentialSubmeshAttributes(*reg_submesh, marker);
  if (ess_attr.Size() == 0)
  {
    return false;  // pure-Neumann region block: keep the single-level preconditioner
  }

  reg_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      order, dim, mg_levels, iodata.solver.linear.mg_coarsening, false);
  reg_dbc_lists.clear();
  reg_hierarchy = std::make_unique<FiniteElementSpaceHierarchy>(
      fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
          mg_levels, reg_mesh_vec, reg_fecs, &ess_attr, &reg_dbc_lists));
  if (reg_hierarchy->GetNumLevels() < 2)
  {
    reg_hierarchy.reset();
    reg_fecs.clear();
    return false;
  }
  reg_solve_fes = &reg_hierarchy->GetFinestFESpace().Get();

  MaterialOperator reg_mat_op(iodata, *reg_mesh_vec[0]);
  MaterialPropertyCoefficient coef(reg_mat_op.GetAttributeToMaterial(),
                                   reg_mat_op.GetPermittivityReal());
  BilinearForm a(reg_hierarchy->GetFinestFESpace());
  a.AddDomainIntegrator<DiffusionIntegrator>(coef);
  auto a_vec = a.Assemble(*reg_hierarchy, false);

  const std::size_t nl = reg_hierarchy->GetNumLevels();
  reg_mg_op = std::make_unique<MultigridOperator>(nl);
  for (std::size_t l = 0; l < nl; l++)
  {
    auto &fes_l = reg_hierarchy->GetFESpaceAtLevel(l);
    auto A_l = std::make_unique<ParOperator>(std::move(a_vec[l]), fes_l);
    A_l->SetEssentialTrueDofs(reg_dbc_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    reg_mg_op->AddOperator(std::move(A_l));
  }
  reg_ess = reg_dbc_lists.back();

  auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
      std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
  amg->SetDropSmallEntries(false);
  reg_gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
      iodata, reg_submesh->GetComm(), std::move(amg),
      reg_hierarchy->GetProlongationOperators());
  reg_gmg->SetOperator(*reg_mg_op);

  reg_sgf.SetSpace(reg_solve_fes);
  reg_srhs.SetSize(reg_solve_fes->GetTrueVSize());
  reg_ssol.SetSize(reg_solve_fes->GetTrueVSize());
  reg_srhs.UseDevice(true);
  reg_ssol.UseDevice(true);
  return true;
}

void SubstructuringSolver::Impl::ApplyRegionGmgPc(const mfem::Vector &r,
                                                  mfem::Vector &z) const
{
  z.SetSize(nt);
  {
    const double *rr = r.HostRead();
    double *zz = z.HostWrite();
    for (int i = 0; i < nt; i++)
    {
      zz[i] = is_region_free[i] ? 0.0 : rr[i];
    }
  }
  reg_pgf.SetFromTrueDofs(r);
  reg_sgf = 0.0;
  reg_submesh->Transfer(reg_pgf, reg_sgf);
  reg_sgf.GetTrueDofs(reg_srhs);
  linalg::SetSubVector(reg_srhs, reg_ess, 0.0);
  reg_ssol = 0.0;
  reg_gmg->Mult(reg_srhs, reg_ssol);
  reg_sgf.SetFromTrueDofs(reg_ssol);
  reg_pgf = 0.0;
  reg_submesh->Transfer(reg_sgf, reg_pgf);
  mfem::Vector zr(nt);
  reg_pgf.GetTrueDofs(zr);
  const double *zrr = zr.HostRead();
  double *zz = z.HostReadWrite();
  for (int i = 0; i < nt; i++)
  {
    if (is_region_free[i])
    {
      zz[i] = zrr[i];
    }
  }
}

std::string SubstructuringSolver::Impl::SetUpRegionSolver()
{
  const int nG = nG_global;
  MPI_Comm comm = parent_fes.GetComm();
  const auto &subcfg = *iodata.solver.substructuring;
  std::string region_solver;
  std::unique_ptr<Solver<Operator>> pc;
  // Direct factorization when the region fits.
  long long reg_glob = 0;
  if (substructuring::kHasDirectSolver)
  {
    long long reg_loc = 0;
    for (int i = 0; i < nt; i++)
    {
      reg_loc += is_region_free[i];
    }
    MPI_Allreduce(&reg_loc, &reg_glob, 1, MPI_LONG_LONG, MPI_SUM, comm);
  }
  MFEM_VERIFY(!broken ||
                  (substructuring::kHasDirectSolver && reg_glob <= subcfg.direct_max_size),
              "\"TwoSided\" superconductor films in substructuring need a direct region "
              "factorization: a sparse direct solver and a region within "
              "\"Solver.Substructuring.DirectMaxSize\"!");
  if (substructuring::kHasDirectSolver && reg_glob <= subcfg.direct_max_size)
  {
    // The exact condensed operator A_region_free + S_E (dense Gamma block; a sparsified
    // S_E is a worse preconditioner than none) or A_region_free, chosen at the first
    // solve batch (EnsureRegionFactor).
    region_direct = true;
    region_dense_eligible = !hodlr && nG > 0 && nG <= kDirectCondensedMaxInterface;
    pc = std::make_unique<CallableSolver>(A_region_free->Height(),
                                          [this](const mfem::Vector &r, mfem::Vector &z)
                                          { ApplyRegionFactor(r, z); });
    region_solver = fmt::format("{} factorization", substructuring::kDirectSolverName);
  }
  else if (BuildRegionGmg())
  {
    pc = std::make_unique<CallableSolver>(A_region_free->Height(),
                                          [this](const mfem::Vector &r, mfem::Vector &z)
                                          { ApplyRegionGmgPc(r, z); });
    region_solver = "conjugate gradients with geometric multigrid";
  }
  else if (magnetostatic)
  {
    region_solver = "conjugate gradients with AMS";
    auto ams = std::make_unique<mfem::HypreAMS>(&parent_fes);
    ams->SetPrintLevel(0);
    pc = std::make_unique<MfemWrapperSolver<Operator>>(
        std::move(ams), /*save_assembled=*/true, /*complex_matrix=*/false,
        /*drop_small_entries=*/false);
  }
  else
  {
    region_solver = "conjugate gradients with BoomerAMG";
    pc = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, true, 0), /*save_assembled=*/true,
        /*complex_matrix=*/false, /*drop_small_entries=*/false);
  }
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(kRegionCgTol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(kRegionCgMaxIts);
  region_op = std::make_unique<substructuring::RegionCondensedOperator>(
      *A_region_free, *mat_dtn, is_region_free);
  region_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
  region_ksp->SetOperators(*region_op, *A_region_free);
  return region_solver;
}

SubstructuringSolver::SubstructuringSolver(const IoData &iodata,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh)
  : impl(std::make_unique<Impl>(iodata, *mesh.back()))
{
}

SubstructuringSolver::~SubstructuringSolver() = default;

void SubstructuringSolver::CondenseEnvironment()
{
  // A_region restricted to region-free true DOFs (non-region-free pinned to identity).
  mfem::Array<int> non_rfree;
  for (int i = 0; i < impl->nt; i++)
  {
    if (!impl->is_region_free[i])
    {
      non_rfree.Append(i);
    }
  }
  impl->A_region_free = EliminateToIdentity(*impl->A_region, non_rfree);

  impl->EnumerateInterface();
  MPI_Comm comm = impl->parent_fes.GetComm();
  const int nG = impl->nG_global;

  // Interface true-DOF geometric signature (replicated), to re-order a saved S_E onto the
  // current interface after re-meshing or re-partitioning the region (with the side moments
  // of the versions of the DOFs split by two-sided films, see InterfaceSignatures).
  MFEM_VERIFY(SignatureWidth(impl->parent_fes) == (impl->magnetostatic ? 12 : 3),
              "Unexpected finite element space for the interface signatures!");
  const int sig_w = impl->magnetostatic ? (impl->InterfaceSplit() ? 13 : 12) : 3;

  // Online with a saved model: load S_E (and the model's sections); otherwise materialize
  // S_E, and save it when a path is set.
  const auto &subcfg = *impl->iodata.solver.substructuring;
  const bool online = (subcfg.mode == SubstructuringMode::ONLINE);
  const std::string &model_path = subcfg.save_model;
  bool loaded = false;
  if (online && !model_path.empty())
  {
    impl->LoadModel(model_path, sig_w);
    loaded = true;
  }
  // MUMPS Schur materialization when available and the environment fits a direct
  // factorization; not when S^K is also materialized (its correction needs the lifted
  // columns, which the batched back-solve path produces).
  bool mumps_done = false;
#if defined(MFEM_USE_MUMPS)
  if (!loaded && nG > 0 && !(impl->EnergyDiffers() && !model_path.empty()) &&
      impl->EnvDirectSize() <= subcfg.direct_max_size)
  {
    impl->MaterializeMumps();
    mumps_done = true;
  }
#endif
  if (!loaded && !mumps_done && impl->iodata.solver.substructuring->factorization_tol > 0.0)
  {
    Mpi::Warning(comm, "Solver.Substructuring.FactorizationTol is only used with the MUMPS "
                       "environment factorization; the environment is factored exactly.\n");
  }
  if (!loaded)
  {
    impl->EnsureEnv();  // no-op after the MUMPS materialization (it is the env solver)
    if (!mumps_done)
    {
      impl->MaterializeByBackSolves(impl->EnergyDiffers() && !model_path.empty());
    }
    if (!impl->magnetostatic)
    {
      // Terminal modes (the flux-loop source modes are computed by SourceEnergyMatrix).
      std::vector<int> ids;
      std::vector<Vector> lifts;
      for (const auto &[idx, dofs] : impl->terminal_tdofs)
      {
        ids.push_back(idx);
        lifts.push_back(impl->TerminalMode(idx));
      }
      impl->ComputeModes(ids, lifts);
    }
    if (!model_path.empty())
    {
      impl->SaveModel(model_path, sig_w);
    }
  }
  // Optional HODLR compression of the interface operators.
  if (subcfg.interface_offdiag_tol > 0.0 && nG > 0)
  {
    impl->CompressInterface(subcfg.interface_offdiag_tol);
  }

  impl->mat_dtn = std::make_unique<substructuring::MaterializedDtN>(
      impl->S_rows, impl->gamma_off, impl->gamma_global, impl->nG_global, comm,
      impl->hodlr.get());
  if (impl->have_sk)
  {
    impl->mat_dtn_K = std::make_unique<substructuring::MaterializedDtN>(
        impl->SK_rows, impl->gamma_off, impl->gamma_global, impl->nG_global, comm,
        impl->hodlr_K.get());
  }

  const std::string region_solver = impl->SetUpRegionSolver();

  // Summary: the sizes, and how the environment was condensed and the region is solved.
  std::string env_solver;
  if (loaded)
  {
    env_solver = fmt::format("loaded from {}", model_path);
  }
  else
  {
    if (mumps_done)
    {
      env_solver = "MUMPS partial factorization";
      if (subcfg.factorization_tol > 0.0)
      {
        env_solver +=
            fmt::format(" (block low-rank, tolerance {:.1e})", subcfg.factorization_tol);
      }
    }
    else if (impl->env_parent_direct)
    {
      env_solver = fmt::format("{} factorization", substructuring::kDirectSolverName);
    }
    else
    {
      env_solver = fmt::format("conjugate gradients with {}",
                               impl->env_hierarchy   ? "geometric multigrid"
                               : impl->magnetostatic ? "AMS"
                                                     : "BoomerAMG");
    }
    if (!model_path.empty())
    {
      env_solver += fmt::format(", saved to {}", model_path);
    }
  }
  long long counts[2] = {0, 0};  // region and environment interiors
  for (int i = 0; i < impl->nt; i++)
  {
    counts[0] += impl->is_region_free[i] && !impl->is_gamma[i];
    counts[1] += impl->is_env_int[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_LONG_LONG, MPI_SUM, comm);
  Mpi::Print(comm,
             "\nSubstructuring, number of global unknowns:\n Region: {:d}, environment: "
             "{:d}, interface: {:d}\n Environment: {}\n Region: {}\n",
             counts[0], counts[1], nG, env_solver, region_solver);
}

namespace
{

// Terminal excitation Dirichlet field: the driven terminal at 1 V, all others grounded.
Vector TerminalDbc(const std::map<int, std::vector<int>> &terminal_tdofs, int drive_idx,
                   int nt)
{
  Vector dbc(nt);
  dbc = 0.0;
  for (const auto &[idx, dofs] : terminal_tdofs)
  {
    const double value = (idx == drive_idx) ? 1.0 : 0.0;
    for (int d : dofs)
    {
      dbc(d) = value;
    }
  }
  return dbc;
}

}  // namespace

Vector SubstructuringSolver::SolveExcitation(int drive_idx)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  return SolveDirichletBatch({TerminalDbc(impl->terminal_tdofs, drive_idx, impl->nt)})[0];
}

std::vector<Vector>
SubstructuringSolver::SolveDirichlets(const std::vector<Vector> &dbc_values)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  // Prescribe an arbitrary field on the Dirichlet DOFs.
  std::vector<Vector> dbcs(dbc_values.size(), Vector(impl->nt));
  Vector xd;
  for (std::size_t k = 0; k < dbc_values.size(); k++)
  {
    dbcs[k] = 0.0;
    dbc_values[k].GetSubVector(impl->dbc_tdofs, xd);
    dbcs[k].SetSubVector(impl->dbc_tdofs, xd);
  }
  return SolveDirichletBatch(dbcs);
}

std::vector<Vector>
SubstructuringSolver::SolveDirichletBatch(const std::vector<Vector> &dbcs)
{
  impl->EnsureEnv();  // general Dirichlet data needs environment solves
  const int nt = impl->nt;
  const int n = static_cast<int>(dbcs.size());
  impl->EnsureRegionFactor(n);
  MPI_Comm comm = impl->parent_fes.GetComm();

  // g_E per excitation: the interface response to the Dirichlet data,
  //   g_E = (A_env x)|_Gamma - (A_env A_EE^-1 (A_env x)|_E)|_Gamma.
  // The environment solve is skipped (exactly) when (A_env x)|_E vanishes, i.e. the
  // Dirichlet data does not touch the environment (e.g. a region terminal); the others run
  // as one batched multi-RHS solve.
  std::vector<Vector> t(n, Vector(nt)), gE(n, Vector(nt));
  std::vector<int> touches(n, 0);
  for (int k = 0; k < n; k++)
  {
    impl->Aenv().Mult(dbcs[k], t[k]);
    const double *tt = t[k].HostRead();
    for (int i = 0; i < nt && !touches[k]; i++)
    {
      if (impl->is_env_int[i] && tt[i] != 0.0)
      {
        touches[k] = 1;
      }
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, touches.data(), n, MPI_INT, MPI_MAX, comm);
  }
  {
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    std::vector<int> ks;
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        X.push_back(&t[k]);
        Y.push_back(&gE[k]);
        ks.push_back(k);
      }
    }
    if (!X.empty())
    {
      impl->ApplyAeeInvMulti(X, Y);
    }
    Vector t2(nt);
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        impl->Aenv().Mult(gE[k], t2);
        gE[k] = t[k];
        gE[k] -= t2;
      }
      else
      {
        gE[k] = t[k];
      }
    }
  }

  // Region-condensed solve per excitation using the materialized S_E (no environment solves
  // in the loop). RHS: region Dirichlet elimination + environment load g_E on the
  // interface. Both are complete (assembled) true-DOF vectors, so each rank reads its own
  // entries.
  std::vector<Vector> x(n, Vector(impl->rf_idx.Size())), u(n);
  {
    Vector tr(nt), g;
    for (int k = 0; k < n; k++)
    {
      impl->A_region->Mult(dbcs[k], tr);
      gE[k].GetSubVector(impl->gamma_idx, g);
      tr.AddElementVector(impl->gamma_idx, g);
      tr.Neg();
      tr.GetSubVector(impl->rf_idx, x[k]);
    }
  }
  impl->SolveRegion(x);
  for (int k = 0; k < n; k++)
  {
    Vector xd(impl->dbc_tdofs.Size());
    dbcs[k].GetSubVector(impl->dbc_tdofs, xd);
    u[k] = impl->RegionField(x[k], &xd);
  }

  impl->RecoverEnvInterior(u);  // batched: u_E = -A_EE^-1 (A_env u)|_E
  return u;
}

Vector SubstructuringSolver::SolveSource(const Vector &f)
{
  return SolveSourceBatch({f})[0];
}

std::vector<Vector> SubstructuringSolver::SolveSourceBatch(const std::vector<Vector> &f)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  impl->EnsureEnv();
  const int nt = impl->nt, n = static_cast<int>(f.size());
  impl->EnsureRegionFactor(n);

  // Environment interior source responses w = A_EE^-1 f_E, in one batch.
  std::vector<Vector> fE(n, Vector(nt)), w(n, Vector(nt));
  std::vector<const mfem::Vector *> X(n);
  std::vector<mfem::Vector *> Y(n);
  for (int k = 0; k < n; k++)
  {
    const double *ff = f[k].HostRead();
    double *fe = fE[k].HostWrite();
    for (int i = 0; i < nt; i++)
    {
      fe[i] = impl->is_env_int[i] ? ff[i] : 0.0;
    }
    w[k] = 0.0;
    X[k] = &fE[k];
    Y[k] = &w[k];
  }
  if (n > 0)
  {
    impl->ApplyAeeInvMulti(X, Y);
  }

  // Region-condensed loads: the region and interface sources minus the environment source
  // correction (A_env w)|_Gamma (assembled: each rank reads its own interface entries).
  std::vector<Vector> x(n, Vector(impl->rf_idx.Size()));
  Vector Aw(nt), b(nt);
  for (int k = 0; k < n; k++)
  {
    impl->Aenv().Mult(w[k], Aw);
    const double *ff = f[k].HostRead(), *aw = Aw.HostRead();
    double *bb = b.HostWrite();
    for (int i = 0; i < nt; i++)
    {
      bb[i] = (impl->is_region_free[i] ? ff[i] : 0.0) - (impl->is_gamma[i] ? aw[i] : 0.0);
    }
    b.GetSubVector(impl->rf_idx, x[k]);
  }
  impl->SolveRegion(x);
  std::vector<Vector> u(n);
  for (int k = 0; k < n; k++)
  {
    u[k] = impl->RegionField(x[k]);
  }
  impl->RecoverEnvInterior(u, &fE);  // u_E = A_EE^-1 (f_E - (A_env u)|_E)
  return u;
}

mfem::DenseMatrix SubstructuringSolver::MagnetostaticEnergyMatrix(
    const std::vector<int> &current_ids, const std::vector<Vector> &J,
    const std::vector<int> &flux_ids, const std::vector<Vector> &a,
    const std::vector<Vector> &apertures, mfem::DenseMatrix *linked_flux,
    std::vector<Vector> *fields, int n_fields, std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(impl->magnetostatic,
              "MagnetostaticEnergyMatrix needs a magnetostatic model!");
  MFEM_VERIFY(current_ids.size() == J.size() && flux_ids.size() == a.size(),
              "MagnetostaticEnergyMatrix: one id per excitation!");
  MFEM_VERIFY(a.empty() || impl->has_sheets,
              "Flux loops need London superconductor sheets!");
  MFEM_VERIFY(apertures.empty() || (apertures.size() == J.size() && linked_flux),
              "MagnetostaticEnergyMatrix: one aperture per surface current!");
  const int nc = static_cast<int>(J.size()), nf = static_cast<int>(a.size());
  std::vector<int> ids(current_ids);
  ids.insert(ids.end(), flux_ids.begin(), flux_ids.end());
  std::vector<Source> src;
  for (const auto &Jk : J)
  {
    src.push_back({nullptr, &Jk});
  }
  for (const auto &ak : a)
  {
    src.push_back({&ak, nullptr});
  }
  // Aperture modes have the negative ids of their ports (indices are positive and unique
  // across surface currents and flux loops).
  Functionals fun;
  for (int k = 0; k < static_cast<int>(apertures.size()); k++)
  {
    fun.ids.push_back(-current_ids[k]);
    fun.l.push_back(&apertures[k]);
  }
  std::vector<double> work;
  mfem::DenseMatrix G = SourceEnergyMatrix(ids, src, fields, n_fields, &work,
                                           fun.l.empty() ? nullptr : &fun, region_fields);

  // Surface currents: G_ij = u_i^T L u_j - J_i^T u_j - J_j^T u_i = -E_ij at the solution,
  // stationary in the regularization. A current that does not close has a null-space part,
  // amplified by 1/eps: its regularization energy J^T u - u^T L u then dominates.
  for (int k = 0; k < nc; k++)
  {
    const double uLu = G(k, k) + 2.0 * work[k], uDu = work[k] - uLu;
    MFEM_VERIFY(uDu <= Impl::kMaxRegularizationFraction * std::abs(uLu),
                "Surface current excitation "
                    << current_ids[k] << " does not close: its current must start and end "
                    << "on conductors or superconducting films (regularization energy "
                    << uDu / std::abs(uLu) << " of the magnetic energy)!");
  }
  // Flux loops: G is the London energy. The two kinds of states are energy-orthogonal, so
  // their cross block is zero.
  mfem::DenseMatrix E(nc + nf);
  E = 0.0;
  for (int i = 0; i < nc + nf; i++)
  {
    for (int j = 0; j < nc + nf; j++)
    {
      if ((i < nc) == (j < nc))
      {
        E(i, j) = (i < nc) ? -G(i, j) : G(i, j);
      }
    }
  }
  if (linked_flux)
  {
    linked_flux->SetSize(nc, nf);
    *linked_flux = 0.0;
    for (int c = 0; c < static_cast<int>(fun.l.size()); c++)
    {
      for (int f = 0; f < nf; f++)
      {
        (*linked_flux)(c, f) = fun.values(c, nc + f);
      }
    }
  }
  return E;
}

mfem::DenseMatrix SubstructuringSolver::SourceEnergyMatrix(
    const std::vector<int> &ids, const std::vector<Source> &src,
    std::vector<Vector> *fields, int n_fields, std::vector<double> *work, Functionals *fun,
    std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  const int nt = impl->nt;
  const int n = static_cast<int>(src.size());
  const int nl = fun ? static_cast<int>(fun->l.size()) : 0;
  const int nf = fields ? std::max(0, std::min(n_fields, n)) : 0;
  MPI_Comm comm = impl->parent_fes.GetComm();
  if (fun)
  {
    fun->values.SetSize(nl, n);
  }

  // Energy of the sources b_k = M_sheet a_k + J_k, with L = K + M_sheet:
  //   G_ij = u_i^T K u_j + (u_i - a_i)^T M_sheet (u_j - a_j) - J_i^T u_j - J_j^T u_i,
  // the native energy for London flux states (J = 0). The kinetic term is formed from the
  // differences (u - a is small for a stiff sheet, and the expanded form would cancel).
  auto current = [&](const Vector &ui, const Vector &uj, const Source &si, const Source &sj)
  { return (si.HasJ() ? *si.J * uj : 0.0) + (sj.HasJ() ? *sj.J * ui : 0.0); };

  // Without the energy interface operator S^K (no model saved or loaded): full solves with
  // the environment, in one batch, and each operator applied once per source.
  if (!impl->mat_dtn_K)
  {
    MFEM_VERIFY(!region_fields, "Region solutions need a saved substructuring model!");
    std::vector<Vector> b(n, Vector(nt));
    for (int k = 0; k < n; k++)
    {
      b[k] = 0.0;
      if (src[k].HasA())
      {
        impl->Msheet().Mult(*src[k].a, b[k]);
      }
      if (src[k].HasJ())
      {
        b[k] += *src[k].J;
      }
      b[k].SetSubVector(impl->dbc_tdofs, 0.0);
    }
    const std::vector<Vector> u = SolveSourceBatch(b);
    std::vector<Vector> d(u), Md(n, Vector(nt));
    for (int j = 0; j < n; j++)
    {
      if (src[j].HasA())
      {
        d[j] -= *src[j].a;
      }
      Md[j] = 0.0;
      if (impl->has_sheets)
      {
        impl->Msheet().Mult(d[j], Md[j]);
      }
    }
    // The kinetic and current terms, the source works b_k^T u_k and the functionals,
    // reduced together.
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    std::vector<double> loc(nn + n + static_cast<std::size_t>(nl) * n, 0.0);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        loc[i * n + j] = (d[i] * Md[j]) - current(u[i], u[j], src[i], src[j]);
      }
      loc[nn + i] = b[i] * u[i];
      for (int l = 0; l < nl; l++)
      {
        loc[nn + n + l * n + i] = *fun->l[l] * u[i];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, loc.data(), static_cast<int>(loc.size()), MPI_DOUBLE,
                  MPI_SUM, comm);
    mfem::DenseMatrix E = impl->EnergyGram(u);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        E(i, j) += loc[i * n + j];
      }
      for (int l = 0; l < nl; l++)
      {
        fun->values(l, i) = loc[nn + n + l * n + i];
      }
    }
    if (work)
    {
      work->assign(loc.begin() + nn, loc.begin() + nn + n);
    }
    if (fields)
    {
      fields->assign(u.begin(), u.begin() + nf);
    }
    return E;
  }

  // Source modes of the sources and functionals with an environment part (one without,
  // such as a port inside the region, has zero modes): reuse the saved / current ones if
  // they match (id + environment source fingerprint), else compute them (needs the
  // environment once; appended to the model in an offline run that saves one).
  std::vector<int> mode_ids(ids);
  std::vector<Source> mode_src(src);
  for (int l = 0; l < nl; l++)
  {
    mode_ids.push_back(fun->ids[l]);
    mode_src.push_back({nullptr, fun->l[l]});
  }
  std::vector<int> col;
  std::vector<char> need;
  if (!impl->SourceModesMatch(mode_ids, mode_src, col, need))
  {
    std::vector<int> req_ids;
    std::vector<Source> req_src;
    for (std::size_t j = 0; j < mode_ids.size(); j++)
    {
      col[j] = need[j] ? static_cast<int>(req_ids.size()) : -1;
      if (need[j])
      {
        req_ids.push_back(mode_ids[j]);
        req_src.push_back(mode_src[j]);
      }
    }
    impl->ComputeSourceModes(req_ids, req_src);
    const auto &subcfg = *impl->iodata.solver.substructuring;
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.save_model.empty())
    {
      impl->AppendSourceModes(subcfg.save_model);
    }
  }
  const int K = static_cast<int>(impl->src_ids.size());
  auto row = [&](int i)
  { return static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K; };
  // Modes of request j (zero without an environment part).
  auto g_of = [&](int i, int j)
  { return col[j] < 0 ? 0.0 : impl->SG_rows[row(i) + col[j]]; };
  auto h_of = [&](int i, int j)
  { return col[j] < 0 ? 0.0 : impl->SH_rows[row(i) + col[j]]; };
  auto pair_of = [&](const std::vector<double> &C, int i, int j)
  {
    return (col[i] < 0 || col[j] < 0) ? 0.0
                                      : C[static_cast<std::size_t>(col[i]) * K + col[j]];
  };

  // Region-condensed solve per source (environment interior left at zero): the region
  // source plus the environment's interface load g_k.
  impl->EnsureRegionFactor(n);
  const auto &rc = impl->rc_idx, &gm = impl->gamma_idx;
  std::vector<Vector> x(n, Vector(impl->rf_idx.Size()));
  for (int j = 0; j < n; j++)
  {
    Vector bR = impl->RegionSource(src[j]);
    double *br = bR.HostReadWrite();
    for (int i : gm)
    {
      br[i] += g_of(i, j);
    }
    bR.GetSubVector(impl->rf_idx, x[j]);
  }
  // Source work b_j^T u_j = (region-condensed load)^T u_rf + B_jj.
  std::vector<Vector> loads;
  if (work)
  {
    loads = x;
  }
  impl->SolveRegion(x);
  if (work)
  {
    work->assign(n, 0.0);
    for (int j = 0; j < n; j++)
    {
      (*work)[j] = loads[j] * x[j];
    }
    loads.clear();
    MPI_Allreduce(MPI_IN_PLACE, work->data(), n, MPI_DOUBLE, MPI_SUM, comm);
    for (int j = 0; j < n; j++)
    {
      (*work)[j] += pair_of(impl->SB, j, j);
    }
  }
  // Functionals: l^T u_j = l_rf^T u_rf + g_l^T u_Gamma + B_lj.
  if (fun)
  {
    Vector l_rf(impl->rf_idx.Size());
    for (int l = 0; l < nl; l++)
    {
      fun->l[l]->GetSubVector(impl->rf_idx, l_rf);
      for (int j = 0; j < n; j++)
      {
        fun->values(l, j) = l_rf * x[j];
      }
    }
  }

  // Per source, on the region elements' DOFs: u_j, K_R u_j, d_j = u_j - a_j, M_sheet^R d_j
  // and the region part of J_j; on Gamma: u_j and S^K u_j. Full fields only for the saved
  // ones (and all of them for region_fields).
  std::vector<Vector> uc(n), Kuc(n), dc(n), Mdc(n), Jc(n), ug(n), Sug(n);
  if (fields)
  {
    fields->clear();
  }
  if (region_fields)
  {
    region_fields->clear();
  }
  {
    Vector Ku(nt), Su(nt), d(nt), Md(nt), JR(nt);
    for (int j = 0; j < n; j++)
    {
      Vector u = impl->RegionField(x[j]);
      x[j].Destroy();
      impl->K_region_e->Mult(u, Ku);
      impl->mat_dtn_K->Mult(u, Su);
      d = u;
      if (src[j].HasA())
      {
        d -= *src[j].a;
      }
      Md = 0.0;
      if (impl->M_sheet_region)
      {
        impl->M_sheet_region->Mult(d, Md);
      }
      JR = 0.0;
      if (src[j].HasJ())
      {
        const double *J = src[j].J->HostRead();
        double *jr = JR.HostWrite();
        for (int i = 0; i < nt; i++)
        {
          jr[i] = impl->is_region_free[i] ? J[i] : 0.0;
        }
      }
      u.GetSubVector(rc, uc[j]);
      Ku.GetSubVector(rc, Kuc[j]);
      d.GetSubVector(rc, dc[j]);
      Md.GetSubVector(rc, Mdc[j]);
      JR.GetSubVector(rc, Jc[j]);
      u.GetSubVector(gm, ug[j]);
      Su.GetSubVector(gm, Sug[j]);
      if (region_fields)
      {
        region_fields->push_back(u);
      }
      if (j < nf)
      {
        fields->push_back(std::move(u));
      }
    }
  }

  // G_ij = [u_i^T K_R u_j + (u_i - a_i)^T M_sheet^R (u_j - a_j) - J_i^R^T u_j
  //         - J_j^R^T u_i]                                          (region)
  //      + u_i,G^T S^K u_j,G + u_i,G^T h_j + h_i^T u_j,G + c_ij     (environment).
  mfem::DenseMatrix E(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double e = (uc[i] * Kuc[j]) + (dc[i] * Mdc[j]) + (ug[i] * Sug[j]) - (Jc[i] * uc[j]) -
                 (Jc[j] * uc[i]);
      for (int q = 0; q < gm.Size(); q++)
      {
        e += ug[i](q) * h_of(gm[q], j) + h_of(gm[q], i) * ug[j](q);
      }
      E(i, j) = e;
    }
  }
  for (int l = 0; l < nl; l++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int q = 0; q < gm.Size(); q++)
      {
        fun->values(l, j) += g_of(gm[q], n + l) * ug[j](q);
      }
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, E.GetData(), n * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  if (nl > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, fun->values.GetData(), nl * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      E(i, j) += pair_of(impl->SC, i, j);
    }
  }
  for (int l = 0; l < nl; l++)
  {
    for (int j = 0; j < n; j++)
    {
      fun->values(l, j) += pair_of(impl->SB, n + l, j);
    }
  }

  // Optional full fields: the environment interior recovered from its sources.
  if (nf > 0)
  {
    std::vector<mfem::Vector> bE(nf);
    std::array<double, 3> fp;
    for (int k = 0; k < nf; k++)
    {
      bE[k] = impl->EnvSource(src[k], fp);
    }
    impl->RecoverEnvInterior(*fields, &bE);
  }
  return E;
}

mfem::DenseMatrix SubstructuringSolver::EnergyMatrix(const std::vector<int> &ids,
                                                     const std::vector<Vector> &lifts,
                                                     std::vector<Vector> *fields,
                                                     int n_fields)
{
  MFEM_VERIFY(ids.size() == lifts.size(), "EnergyMatrix: one id per lift!");
  // Dirichlet data restricted to the Dirichlet DOF set (as SolveDirichlets does).
  std::vector<Vector> xd(lifts.size(), Vector(impl->dbc_tdofs.Size()));
  for (std::size_t k = 0; k < lifts.size(); k++)
  {
    lifts[k].GetSubVector(impl->dbc_tdofs, xd[k]);
  }
  return EnergyMatrixDbc(ids, xd, fields, n_fields, nullptr);
}

mfem::DenseMatrix SubstructuringSolver::EnergyMatrixDbc(const std::vector<int> &ids,
                                                        const std::vector<Vector> &xd,
                                                        std::vector<Vector> *fields,
                                                        int n_fields,
                                                        std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  const int nt = impl->nt;
  const int n = static_cast<int>(xd.size());
  const int nf = fields ? std::max(0, std::min(n_fields, n)) : 0;
  MPI_Comm comm = impl->parent_fes.GetComm();
  auto full = [&]()
  {
    std::vector<Vector> x;
    for (const auto &v : xd)
    {
      x.push_back(impl->DirichletField(v));
    }
    return x;
  };

  // Magnetostatics without the energy interface operator S^K: environment path (full
  // fields, then the energy u_i^T K u_j).
  if (impl->EnergyDiffers() && !impl->mat_dtn_K)
  {
    std::vector<Vector> u = SolveDirichletBatch(full());
    mfem::DenseMatrix E = impl->EnergyGram(u);
    if (fields)
    {
      fields->assign(u.begin(), u.begin() + nf);
    }
    return E;
  }

  // Modes of these lifts: reuse the saved / current ones if they match (id + environment
  // fingerprint), else compute them (needs the environment once; appended to the model in
  // an offline run that saves one).
  std::vector<int> col;
  if (!impl->ModesMatch(ids, xd, col))
  {
    impl->ComputeModes(ids, full());
    col.resize(n);
    std::iota(col.begin(), col.end(), 0);
    const auto &subcfg = *impl->iodata.solver.substructuring;
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.save_model.empty())
    {
      impl->AppendModes(subcfg.save_model);
    }
  }
  const int K = static_cast<int>(impl->mode_ids.size());
  const std::vector<double> &GE = impl->have_gk ? impl->GK_rows : impl->G_rows;
  const substructuring::MaterializedDtN &SK =
      impl->mat_dtn_K ? *impl->mat_dtn_K : *impl->mat_dtn;
  auto row = [&](int i)
  { return static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K; };

  // Region-condensed solve per lift (environment interior left at zero): the interface load
  // is the precomputed g_k, so no environment solve is needed.
  impl->EnsureRegionFactor(n);
  const auto &rc = impl->rc_idx, &gm = impl->gamma_idx;
  std::vector<Vector> xr(n, Vector(impl->rf_idx.Size()));
  {
    Vector tr(nt);
    for (int j = 0; j < n; j++)
    {
      impl->A_region->Mult(impl->DirichletField(xd[j]), tr);
      double *trr = tr.HostReadWrite();
      for (int i : gm)
      {
        trr[i] += impl->G_rows[row(i) + col[j]];
      }
      tr.Neg();
      tr.GetSubVector(impl->rf_idx, xr[j]);
    }
  }
  impl->SolveRegion(xr);

  // Per lift, on the region elements' DOFs: u_j and K_R u_j; on Gamma: u_j and S^K u_j.
  // Full fields only for the saved ones (and all of them for region_fields).
  std::vector<Vector> uc(n), Kuc(n), ug(n), Sug(n);
  if (fields)
  {
    fields->clear();
  }
  if (region_fields)
  {
    region_fields->clear();
  }
  {
    Vector Ku(nt), Su(nt);
    for (int j = 0; j < n; j++)
    {
      Vector u = impl->RegionField(xr[j], &xd[j]);
      xr[j].Destroy();
      impl->K_region_e->Mult(u, Ku);
      SK.Mult(u, Su);
      u.GetSubVector(rc, uc[j]);
      Ku.GetSubVector(rc, Kuc[j]);
      u.GetSubVector(gm, ug[j]);
      Su.GetSubVector(gm, Sug[j]);
      if (j < nf)
      {
        fields->push_back(u);
      }
      if (region_fields)
      {
        region_fields->push_back(std::move(u));
      }
    }
  }

  // E_ij = u_i^T K_R u_j + [u_G,i; e_i]^T [S^K, G^K; G^K^T, Cmode] [u_G,j; e_j]: the region
  // energy plus the environment's exact condensed energy.
  mfem::DenseMatrix E(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double a = (uc[i] * Kuc[j]) + (ug[i] * Sug[j]);
      for (int q = 0; q < gm.Size(); q++)
      {
        a += ug[i](q) * GE[row(gm[q]) + col[j]] + GE[row(gm[q]) + col[i]] * ug[j](q);
      }
      E(i, j) = a;
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, E.GetData(), n * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      E(i, j) += impl->Cmode[static_cast<std::size_t>(col[i]) * K + col[j]];
    }
  }

  // Optional full fields (environment interior recovered on demand).
  if (nf > 0)
  {
    impl->RecoverEnvInterior(*fields);
  }
  return E;
}

mfem::DenseMatrix
SubstructuringSolver::CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                        std::vector<Vector> *fields, int n_fields,
                                        std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(!impl->magnetostatic,
              "CapacitanceMatrix is for electrostatic substructuring problems!");
  std::vector<Vector> xd;
  for (int idx : terminal_indices)
  {
    xd.emplace_back(impl->dbc_tdofs.Size());
    TerminalLift(idx).GetSubVector(impl->dbc_tdofs, xd.back());
  }
  return EnergyMatrixDbc(terminal_indices, xd, fields, n_fields, region_fields);
}

ErrorIndicator
SubstructuringSolver::RegionErrorIndicator(const std::vector<Vector> &region_fields,
                                           const mfem::DenseMatrix &E) const
{
  BlockTimer bt(Timer::ESTIMATION);
  MFEM_VERIFY(E.Height() == static_cast<int>(region_fields.size()),
              "RegionErrorIndicator: one energy per field!");
  const auto &iodata = impl->iodata;
  const int order = iodata.solver.order, dim = impl->parent.Dimension();
  const double tol = iodata.solver.linear.estimator_tol;
  const int max_it = iodata.solver.linear.estimator_max_it;
  const auto coarsening = iodata.solver.linear.mg_coarsening;

  // Palace's flux estimators on the region submesh: their flux recovery sees only the
  // region, so the unknown environment field does not pollute the indicators near the
  // interface.
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(impl->parent, impl->ra_arr))));
  mesh[0]->RebuildCeedAttributes();
  auto &submesh = static_cast<mfem::ParSubMesh &>(mesh[0]->Get());
  MaterialOperator mat_op(iodata, *mesh[0]);

  // Each field on the submesh, its flux s D u, and its indicator.
  mfem::ParGridFunction pgf(&impl->parent_fes);
  mfem::Vector u0;
  ErrorIndicator region;
  auto estimate = [&](const auto &estimator, const auto &D, double s,
                      FiniteElementSpace &fespace, FiniteElementSpace &flux_fespace)
  {
    mfem::ParGridFunction sgf(&fespace.Get());
    Vector u(fespace.GetTrueVSize()), f(flux_fespace.GetTrueVSize());
    u.UseDevice(true);
    f.UseDevice(true);
    for (std::size_t k = 0; k < region_fields.size(); k++)
    {
      pgf.SetFromTrueDofs(impl->FirstVersion(region_fields[k], u0));
      sgf = 0.0;
      submesh.Transfer(pgf, sgf);
      sgf.GetTrueDofs(u);
      f = 0.0;
      D.AddMult(u, f, s);
      const double energy = 0.5 * E(static_cast<int>(k), static_cast<int>(k));
      estimator.AddErrorIndicator(f, energy, region);
    }
  };
  if (impl->magnetostatic)
  {
    // B = curl A.
    auto nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(order, dim, 1,
                                                                      coarsening, false);
    auto nd_fespaces =
        fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(1, mesh, nd_fecs);
    auto &nd_fespace = nd_fespaces.GetFinestFESpace();
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace rt_fespace(*mesh[0], &rt_fec);
    CurlFluxErrorEstimator<Vector> estimator(mat_op, rt_fespace, nd_fespaces, tol, max_it,
                                             0, false);
    estimate(estimator, rt_fespace.GetDiscreteInterpolator(nd_fespace), 1.0, nd_fespace,
             rt_fespace);
  }
  else
  {
    // E = -grad V.
    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(*mesh[0], &h1_fec), nd_fespace(*mesh[0], &nd_fec);
    auto rt_fecs = fem::ConstructFECollections<mfem::RT_FECollection>(order - 1, dim, 1,
                                                                      coarsening, false);
    auto rt_fespaces =
        fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(1, mesh, rt_fecs);
    GradFluxErrorEstimator<Vector> estimator(mat_op, nd_fespace, rt_fespaces, tol, max_it,
                                             0, false);
    estimate(estimator, nd_fespace.GetDiscreteInterpolator(h1_fespace), -1.0, h1_fespace,
             nd_fespace);
  }

  // Region element indicators onto the parent mesh (zero on the environment).
  Vector local(impl->parent.GetNE());
  local = 0.0;
  if (region.Local().Size() > 0)
  {
    const auto &parent_id = submesh.GetParentElementIDMap();
    const double *r = region.Local().HostRead();
    for (int i = 0; i < parent_id.Size(); i++)
    {
      local(parent_id[i]) = r[i];
    }
  }
  return ErrorIndicator(std::move(local));
}

Vector SubstructuringSolver::TerminalLift(int terminal_index) const
{
  MFEM_VERIFY(impl->terminal_tdofs.contains(terminal_index),
              "Unknown terminal index " << terminal_index << "!");
  return impl->TerminalMode(terminal_index);
}

bool SubstructuringSolver::EnvironmentFactored() const
{
  return impl->env_built;
}

std::vector<int> SubstructuringSolver::TerminalIndices() const
{
  std::vector<int> idx;
  for (const auto &[i, dofs] : impl->terminal_tdofs)
  {
    idx.push_back(i);
  }
  return idx;
}

double SubstructuringSolver::MutualEnergy(const Vector &ui, const Vector &uj) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->K_region_e->Mult(uj, t);
  impl->Kenv().Mult(uj, t2);
  t += t2;
  return mfem::InnerProduct(impl->parent_fes.GetComm(), ui, t);
}

long long int SubstructuringSolver::GlobalTrueVSize() const
{
  return impl->SolutionSpace().GlobalTrueVSize();
}

void SubstructuringSolver::WriteParaView(const std::string &dir,
                                         const std::vector<int> &ids,
                                         const std::vector<Vector> &fields) const
{
  BlockTimer bt(Timer::POSTPRO_PARAVIEW);
  MFEM_VERIFY(ids.size() == fields.size(), "WriteParaView requires one field per id!");
  const int order = impl->iodata.solver.order;
  mfem::ParGridFunction phi(&impl->parent_fes);
  mfem::ParaViewDataCollection pv("paraview", &impl->parent);
  pv.SetPrefixPath(dir);
  pv.SetLevelsOfDetail(order);
  pv.SetHighOrderOutput(order > 1);
  pv.SetDataFormat(mfem::VTKFormat::BINARY);
  pv.RegisterField(impl->magnetostatic ? "A" : "V", &phi);
  mfem::Vector u0;
  for (std::size_t j = 0; j < fields.size(); j++)
  {
    phi.SetFromTrueDofs(impl->FirstVersion(fields[j], u0));
    phi.HostRead();  // the writer reads the host copy
    pv.SetCycle(static_cast<int>(j));
    pv.SetTime(static_cast<double>(ids[j]));
    pv.Save();
  }
}

}  // namespace palace
