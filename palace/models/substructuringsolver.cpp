// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include "fem/gridfunction.hpp"
#include "fem/multigrid.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/iterative.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/geodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Preconditioner forwarding to a callable (the region factor, built on the first solve
// batch, so SetOperator is a no-op).
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

// The solver configuration of the native operator of the substructuring: the finest mesh
// and space only (the factorizations need no multigrid hierarchy).
config::SolverData SingleLevel(const config::SolverData &solver)
{
  config::SolverData single(solver);
  single.linear.mg_max_levels = 1;
  return single;
}

std::unique_ptr<LaplaceOperator>
MakeLaplaceOperator(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  // The terminals of the environment are not on the mesh of the region (their modes are in
  // the model), which the operator would reject.
  config::BoundaryData boundaries(iodata.boundaries);
  const auto &bdr_attributes = mesh.back()->Get().bdr_attributes;
  for (auto it = boundaries.terminal.begin(); it != boundaries.terminal.end();)
  {
    std::erase_if(it->second.attributes, [&](int a) { return bdr_attributes.Find(a) < 0; });
    it = it->second.attributes.empty() ? boundaries.terminal.erase(it) : std::next(it);
  }
  return std::make_unique<LaplaceOperator>(boundaries, SingleLevel(iodata.solver),
                                           iodata.domains.materials, iodata.problem.type,
                                           mesh);
}

// The entrywise absolute value of a matrix, zero on the rows not marked (all rows are
// marked when the marker is empty).
std::unique_ptr<mfem::HypreParMatrix> AbsRows(const mfem::HypreParMatrix &A,
                                              const std::vector<bool> &marker)
{
  auto B = std::make_unique<mfem::HypreParMatrix>(A);
  B->HostReadWrite();
  hypre_ParCSRMatrix *h = *B;
  for (hypre_CSRMatrix *M : {hypre_ParCSRMatrixDiag(h), hypre_ParCSRMatrixOffd(h)})
  {
    const HYPRE_Int *I = hypre_CSRMatrixI(M);
    HYPRE_Real *data = hypre_CSRMatrixData(M);
    for (int i = 0; i < hypre_CSRMatrixNumRows(M); i++)
    {
      const bool keep = marker.empty() || marker[i];
      for (int k = I[i]; k < I[i + 1]; k++)
      {
        data[k] = keep ? std::abs(data[k]) : 0.0;
      }
    }
  }
  B->HypreReadWrite();
  return B;
}

// The connected components of the graph of the nonzero entries of a symmetric matrix,
// among its rows with a nonzero diagonal entry: each such row is labeled by the least
// global row of its component, the others by -1. Labels propagate along the local rows
// until none changes, then across processes through the off-diagonal columns, until no
// label changes on any process. Collective.
std::vector<HYPRE_BigInt> ConnectedComponents(const mfem::HypreParMatrix &A)
{
  A.HostRead();
  hypre_ParCSRMatrix *h = A;
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(h),
                        *offd = hypre_ParCSRMatrixOffd(h);
  const HYPRE_Int *di = hypre_CSRMatrixI(diag), *dj = hypre_CSRMatrixJ(diag),
                  *oi = hypre_CSRMatrixI(offd), *oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Real *dv = hypre_CSRMatrixData(diag), *ov = hypre_CSRMatrixData(offd);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(h);
  const int n = hypre_CSRMatrixNumRows(diag), n_offd = hypre_CSRMatrixNumCols(offd);
  const HYPRE_BigInt first = hypre_ParCSRMatrixFirstRowIndex(h);
  MPI_Comm comm = A.GetComm();
  const int np = Mpi::Size(comm);

  // The owners of the off-diagonal columns (in increasing order), and the rows each process
  // reads from this one.
  std::vector<HYPRE_BigInt> starts(np + 1);
  MPI_Allgather(&first, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  starts[np] = A.GetGlobalNumRows();
  std::vector<int> recv_count(np, 0), send_count(np), recv_disp(np + 1, 0),
      send_disp(np + 1, 0);
  for (int k = 0; k < n_offd; k++)
  {
    recv_count[std::upper_bound(starts.begin(), starts.end(), cmap[k]) - starts.begin() -
               1]++;
  }
  MPI_Alltoall(recv_count.data(), 1, MPI_INT, send_count.data(), 1, MPI_INT, comm);
  for (int p = 0; p < np; p++)
  {
    recv_disp[p + 1] = recv_disp[p] + recv_count[p];
    send_disp[p + 1] = send_disp[p] + send_count[p];
  }
  std::vector<HYPRE_BigInt> send_rows(send_disp[np]);
  MPI_Alltoallv(cmap, recv_count.data(), recv_disp.data(), HYPRE_MPI_BIG_INT,
                send_rows.data(), send_count.data(), send_disp.data(), HYPRE_MPI_BIG_INT,
                comm);

  constexpr HYPRE_BigInt none = std::numeric_limits<HYPRE_BigInt>::max();
  std::vector<HYPRE_BigInt> label(n, none), halo(n_offd), send(send_rows.size());
  for (int i = 0; i < n; i++)
  {
    for (int k = di[i]; k < di[i + 1]; k++)
    {
      if (dj[k] == i && dv[k] != 0.0)
      {
        label[i] = first + i;
      }
    }
  }
  bool changed = true;
  while (changed)
  {
    for (std::size_t k = 0; k < send_rows.size(); k++)
    {
      send[k] = label[send_rows[k] - first];
    }
    MPI_Alltoallv(send.data(), send_count.data(), send_disp.data(), HYPRE_MPI_BIG_INT,
                  halo.data(), recv_count.data(), recv_disp.data(), HYPRE_MPI_BIG_INT,
                  comm);
    changed = false;
    for (bool sweep = true; sweep;)
    {
      sweep = false;
      for (int i = 0; i < n; i++)
      {
        if (label[i] == none)
        {
          continue;
        }
        HYPRE_BigInt m = label[i];
        for (int k = di[i]; k < di[i + 1]; k++)
        {
          m = (dv[k] != 0.0) ? std::min(m, label[dj[k]]) : m;
        }
        for (int k = oi[i]; k < oi[i + 1]; k++)
        {
          m = (ov[k] != 0.0) ? std::min(m, halo[oj[k]]) : m;
        }
        if (m < label[i])
        {
          label[i] = m;
          sweep = changed = true;
        }
      }
    }
    Mpi::GlobalOr(1, &changed, comm);
  }
  for (auto &l : label)
  {
    l = (l == none) ? -1 : l;
  }
  A.HypreRead();
  return label;
}

}  // namespace

SubstructuringSolver::Impl::Impl(const IoData &iodata,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh_)
  : iodata(iodata), mesh(*mesh_.back()), parent(mesh.Get()),
    magnetostatic(iodata.problem.type == ProblemType::MAGNETOSTATIC),
    laplace_op(magnetostatic ? nullptr : MakeLaplaceOperator(iodata, mesh_)),
    curlcurl_op(magnetostatic ? std::make_unique<CurlCurlOperator>(
                                    iodata.boundaries, SingleLevel(iodata.solver),
                                    iodata.domains.materials, iodata.problem.type,
                                    iodata.units, mesh_)
                              : nullptr),
    parent_space(magnetostatic ? curlcurl_op->GetNDSpace() : laplace_op->GetH1Space()),
    parent_fes(parent_space.Get())
{
  MFEM_VERIFY(iodata.solver.substructuring,
              "SubstructuringSolver requires a Solver.Substructuring configuration!");
  const auto &sub = *iodata.solver.substructuring;
  env_attributes = sub.environment_attributes;
  if (env_attributes.empty())
  {
    for (int a : parent.attributes)
    {
      if (std::ranges::find(sub.region_attributes, a) == sub.region_attributes.end())
      {
        env_attributes.push_back(a);
      }
    }
  }
  MFEM_VERIFY(!env_attributes.empty() || sub.region_mesh,
              "The substructuring \"Region\" holds every domain attribute of the mesh: the "
              "environment is empty!");
  CheckSubstructuringAttributes(parent, sub.region_attributes, env_attributes);
  if (magnetostatic)
  {
    SetUpBrokenSpace();
  }
  nt = SolutionSpace().GetTrueVSize();
  for (int a : sub.region_attributes)
  {
    ra_arr.Append(a);
  }
  for (int a : env_attributes)
  {
    ea_arr.Append(a);
  }

  // Interface / region / environment true-DOF markers.
  mfem::Array<int> rm, em;
  MarkInterfaceTrueDofs(SolutionSpace(), ra_arr, ea_arr, rm, em);

  // Dirichlet true DOFs: the terminals (per terminal) and the grounded boundaries. The
  // set is fixed across excitations, so the interface/interior partition is too.
  // Terminals are Dirichlet boundaries in electrostatics only: the regular curl-curl
  // operator ignores them (CurlCurlOperator::SetUpBoundaryProperties).
  const std::map<int, config::TerminalData> no_terminals;
  const auto &terminals = magnetostatic ? no_terminals : iodata.boundaries.terminal;
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
  if (sub.mode == SubstructuringMode::ONLINE && !sub.model_file.empty())
  {
    {
      BlockTimer bt(Timer::IO);
      saved = ReadSavedHead(sub.model_file, parent_fes.GetComm());
    }
    const std::string &model = sub.model_file;
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
    // The mesh holds the environment of the model, or with "RegionMesh" its part in the
    // mesh of the region, which cannot be refined.
    region_only = sub.region_mesh;
    if (region_only)
    {
      MFEM_VERIFY(iodata.model.refinement.max_it == 0,
                  "Adaptive refinement with substructuring needs the full mesh, not the "
                  "mesh of the region!");
      const auto &refinement = iodata.model.refinement;
      MFEM_VERIFY(refinement.uniform_ref_levels == 0 &&
                      refinement.ser_uniform_ref_levels == 0 &&
                      refinement.GetBoxes().empty() && refinement.GetSpheres().empty(),
                  kRegionMeshRefined);
    }
    saved.environment.Check(parent, *iodata.resolved_config, env_attributes,
                            kFingerprintTol, region_only);
    // S_E was condensed with the regularization of the model: the region uses the same.
    model_eps = saved.regularization;
    if (region_only)
    {
      // The versions on the region side only of the DOFs on the boundary, at the points of
      // the saved interface DOFs, are interface DOFs (both sides of a two-sided film
      // crossing Gamma).
      std::vector<char> region_side(nt);
      for (int t = 0; t < nt; t++)
      {
        region_side[t] = rm[t] && !em[t] && !dir_mark[t];
      }
      cut = RegionMeshInterface(parent_fes, region_side,
                                broken ? version_parent : std::vector<int>{}, saved.sig,
                                saved.width);
    }
  }
  // As in the regular solve, the terminal attributes are boundary attributes of the mesh
  // (but those of the environment terminals are not on the mesh of the region).
  if (!region_only)
  {
    for (const auto &[idx, term] : terminals)
    {
      for (int a : term.attributes)
      {
        MFEM_VERIFY(parent.bdr_attributes.Find(a) >= 0,
                    "Unknown terminal boundary attribute " << a << "!");
      }
    }
  }

  is_gamma.assign(nt, 0);
  is_env_int.assign(nt, 0);
  is_region_free.assign(nt, 0);
  env_dbc.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    if (dir_mark[i])
    {
      env_dbc[i] = em[i];
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

  if (magnetostatic)
  {
    CheckSheets();
    if (!saved.sheet_sides.empty())
    {
      bool same = std::abs(gamma_sheet_sides[0] - saved.sheet_sides[0]) < 0.5;
      for (int j = 1; j < 4; j++)
      {
        same = same && std::abs(gamma_sheet_sides[j] - saved.sheet_sides[j]) <=
                           kFingerprintTol * std::max(1.0, saved.sheet_sides[0]);
      }
      MFEM_VERIFY(same, "The sheet faces on the interface are not on the same sides as in "
                        "the saved substructuring model: a face belongs to the side of its "
                        "neighbor with the lower attribute, and the attributes next to the "
                        "interface changed. Keep them, or rerun in \"Offline\" mode!");
    }
    mag_eps = model_eps;
    if (!(mag_eps > 0.0))
    {
      mfem::Vector bbmin, bbmax;
      mesh::GetAxisAlignedBoundingBox(parent, bbmin, bbmax);
      bbmax -= bbmin;
      const double L = bbmax.Max();
      mag_eps = (has_sheets ? kMagRegularizationSheets : kMagRegularization) / (L * L);
    }
  }
  AssembleRegion();
}

void SubstructuringSolver::Impl::EnsureRegionFactor(int n)
{
  if (reg_lu)
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
  const bool serial = !region_exact && nrf <= kSerialRegionFactorDofs;
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

void SubstructuringSolver::Impl::SolveRegion(std::vector<mfem::Vector> &x)
{
  const int n = static_cast<int>(x.size());
  EnsureRegionFactor(n);
  BlockTimer bt(Timer::KSP);
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

const Operator &SubstructuringSolver::Impl::ClosureGradient() const
{
  return SolutionSpace().GetDiscreteInterpolator(
      curlcurl_op->GetSolutionH1Spaces().GetFinestFESpace());
}

void SubstructuringSolver::Impl::SetUpClosure()
{
  // |G| (G is assembled for a broken space).
  const Operator &G_op = ClosureGradient();
  const auto *G = dynamic_cast<const mfem::HypreParMatrix *>(&G_op);
  if (!G)
  {
    G = &dynamic_cast<const ParOperator &>(G_op).ParallelAssemble();
  }
  closure_G_abs = AbsRows(*G, {});

  // The conductor DOFs of the ND space: those of the PEC boundaries and of the sheets (the
  // rows of the sheet mass on both sides with a nonzero diagonal), and the connected
  // conductors: the components of the graph of |G_C|^T |G_C|, G_C the rows of G on them.
  const int n = SolutionSpace().GetTrueVSize();
  std::vector<bool> on(n, false);
  for (int t : dbc_tdofs)
  {
    on[t] = true;
  }
  if (has_sheets)
  {
    Vector d;
    Msheet().GetDiag(d);
    const double *h_d = d.HostRead();
    for (int t = 0; t < n; t++)
    {
      on[t] = on[t] || (h_d[t] != 0.0);
    }
  }
  {
    auto G_C = AbsRows(*G, on);
    std::unique_ptr<mfem::HypreParMatrix> G_Ct(G_C->Transpose());
    std::unique_ptr<mfem::HypreParMatrix> H(mfem::ParMult(G_Ct.get(), G_C.get(), true));
    conductor = ConnectedComponents(*H);
  }
  if (!region_only)
  {
    return;
  }

  // On the mesh of the region: the H1 DOFs of the environment elements it holds (those of
  // their ND DOFs off the region elements and the Dirichlet DOFs), whose currents close in
  // the environment of the model, and the conductors touching them.
  Vector g(n), t(closure_G_abs->Width());
  g.UseDevice(true);
  t.UseDevice(true);
  g = 1.0;
  g.SetSubVector(rc_idx, 0.0);
  g.SetSubVector(dbc_tdofs, 0.0);
  closure_G_abs->MultTranspose(g, t);
  const double *h_t = t.HostRead();
  std::vector<HYPRE_BigInt> touched;
  for (std::size_t i = 0; i < conductor.size(); i++)
  {
    if (h_t[i] > 0.0)
    {
      if (conductor[i] >= 0)
      {
        touched.push_back(conductor[i]);
      }
      conductor[i] = kUnchecked;
    }
  }
  MPI_Comm comm = parent_fes.GetComm();
  int count = static_cast<int>(touched.size());
  std::vector<int> counts(Mpi::Size(comm)), displs(Mpi::Size(comm) + 1, 0);
  MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  for (int p = 0; p < Mpi::Size(comm); p++)
  {
    displs[p + 1] = displs[p] + counts[p];
  }
  std::vector<HYPRE_BigInt> all(displs.back());
  MPI_Allgatherv(touched.data(), count, HYPRE_MPI_BIG_INT, all.data(), counts.data(),
                 displs.data(), HYPRE_MPI_BIG_INT, comm);
  const std::set<HYPRE_BigInt> unchecked(all.begin(), all.end());
  for (auto &c : conductor)
  {
    c = (c >= 0 && unchecked.contains(c)) ? kUnchecked : c;
  }
}

bool SubstructuringSolver::Impl::CurrentCloses(const Vector &J)
{
  if (!closure_G_abs)
  {
    SetUpClosure();
  }
  const int n = static_cast<int>(conductor.size());
  Vector r(n), c(n), J_abs(J);
  r.UseDevice(true);
  c.UseDevice(true);
  J_abs.Abs();
  ClosureGradient().MultTranspose(J, r);
  closure_G_abs->MultTranspose(J_abs, c);
  const double *h_r = r.HostRead(), *h_c = c.HostRead();
  MPI_Comm comm = parent_fes.GetComm();
  double scale = (n > 0) ? *std::max_element(h_c, h_c + n) : 0.0;
  Mpi::GlobalMax(1, &scale, comm);

  // The DOFs off the conductors, and the sums over each conductor on this process.
  bool closes = true;
  std::map<HYPRE_BigInt, std::array<double, 2>> sums;
  for (int i = 0; i < n; i++)
  {
    if (conductor[i] == kOffConductor)
    {
      closes = closes && std::abs(h_r[i]) <= kClosureTol * scale;
    }
    else if (conductor[i] != kUnchecked)
    {
      auto &s = sums[conductor[i]];
      s[0] += h_r[i];
      s[1] += h_c[i];
    }
  }
  bool open = !closes;
  Mpi::GlobalOr(1, &open, comm);

  // The sums over each conductor, from all processes.
  std::vector<double> loc;
  for (const auto &[label, s] : sums)
  {
    loc.insert(loc.end(), {static_cast<double>(label), s[0], s[1]});
  }
  int count = static_cast<int>(loc.size());
  std::vector<int> counts(Mpi::Size(comm)), displs(Mpi::Size(comm) + 1, 0);
  MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  for (int p = 0; p < Mpi::Size(comm); p++)
  {
    displs[p + 1] = displs[p] + counts[p];
  }
  std::vector<double> all(displs.back());
  MPI_Allgatherv(loc.data(), count, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
                 MPI_DOUBLE, comm);
  std::map<double, std::array<double, 2>> total;
  for (std::size_t k = 0; k < all.size(); k += 3)
  {
    total[all[k]][0] += all[k + 1];
    total[all[k]][1] += all[k + 2];
  }
  for (const auto &[label, s] : total)
  {
    open = open || std::abs(s[0]) > kClosureTol * std::max(scale, s[1]);
  }
  return !open;
}

std::string SubstructuringSolver::Impl::SetUpRegionSolver()
{
  MPI_Comm comm = parent_fes.GetComm();
  // The exact condensed operator A_region_free + S_E (dense Gamma block; a sparsified S_E
  // is a worse preconditioner than none) or A_region_free, chosen at the first solve batch
  // (EnsureRegionFactor).
  region_dense_eligible = nG_global > 0 && nG_global <= kDirectCondensedMaxInterface;
  auto pc = std::make_unique<CallableSolver>(A_region_free->Height(),
                                             [this](const mfem::Vector &r, mfem::Vector &z)
                                             { ApplyRegionFactor(r, z); });
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(kRegionCgTol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(kRegionCgMaxIts);
  region_op = std::make_unique<substructuring::RegionCondensedOperator>(
      *A_region_free, *mat_dtn, is_region_free);
  region_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
  region_ksp->SetOperators(*region_op, *A_region_free);
  return fmt::format("{} factorization", substructuring::DirectSolverName(iodata));
}

SubstructuringSolver::SubstructuringSolver(const IoData &iodata,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh)
  : impl(std::make_unique<Impl>(iodata, mesh))
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
  const int sig_w = impl->magnetostatic ? (impl->InterfaceSplit() ? 13 : 12) : 3;

  // Online with a saved model: load S_E (and the model's sections); otherwise materialize
  // S_E, and save it when a path is set.
  const auto &subcfg = *impl->iodata.solver.substructuring;
  const bool online = (subcfg.mode == SubstructuringMode::ONLINE);
  const std::string &model_path = subcfg.model_file;
  if (!online && !model_path.empty())
  {
    // Fail before the condensation if the model cannot be written (an existing one is
    // kept until it is replaced).
    substructuring::ModelFile::CheckWritable(model_path, impl->parent_fes.GetComm());
  }
  // The versions of the interface DOFs on a two-sided film crossing Gamma are told apart
  // by the faces of their sides (InterfaceSignatures), of which a coarse face of a
  // nonconforming mesh has several on one side.
  MFEM_VERIFY(
      model_path.empty() || !impl->InterfaceSplit() || !impl->parent.Nonconforming(),
      "A saved substructuring model with \"TwoSided\" superconductor films crossing "
      "the interface needs a conforming mesh!");
  bool loaded = false;
  if (online && !model_path.empty())
  {
    impl->LoadModel(model_path, sig_w);
    loaded = true;
  }
  // S_E from a partial factorization when one is built, else from solves with A_EE (with
  // S^K materialized, for a magnetostatic model to save).
  const bool with_sk = impl->EnergyDiffers() && !model_path.empty();
  if (!loaded && nG > 0)
  {
    impl->MaterializePartial();
  }
  if (!loaded)
  {
    impl->EnsureEnv();  // no-op after a partial factorization (it is the env solver)
    if (with_sk)
    {
      impl->MaterializeWithEnergy();
    }
    else
    {
      impl->S_rows = impl->env_schur->ReleaseSchurRows();
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
  impl->mat_dtn = std::make_unique<substructuring::MaterializedDtN>(
      impl->S_rows, impl->gamma_off, impl->gamma_global, impl->nG_global, comm);
  if (impl->have_sk)
  {
    impl->mat_dtn_K = std::make_unique<substructuring::MaterializedDtN>(
        impl->SK_rows, impl->gamma_off, impl->gamma_global, impl->nG_global, comm);
  }
  if (impl->have_skin)
  {
    impl->mat_dtn_Kin = std::make_unique<substructuring::MaterializedDtN>(
        impl->SKin_rows, impl->gamma_off, impl->gamma_global, impl->nG_global, comm);
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
    if (!impl->env_partial.empty())
    {
      env_solver = fmt::format("{} partial factorization", impl->env_partial);
    }
    else
    {
      env_solver =
          fmt::format("{} factorization", substructuring::DirectSolverName(impl->iodata));
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
SubstructuringSolver::SolveDirichletBatch(const std::vector<Vector> &dbcs)
{
  impl->EnsureEnv();  // general Dirichlet data needs environment solves
  const int nt = impl->nt;
  const int n = static_cast<int>(dbcs.size());
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
    std::vector<Vector> *fields, int n_fields, std::vector<Vector> *region_fields,
    std::vector<double> *volume_energy)
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
  // A current that does not close has no solution.
  for (int k = 0; k < nc; k++)
  {
    MFEM_VERIFY(impl->CurrentCloses(J[k]),
                "Surface current excitation "
                    << current_ids[k] << " does not close: its current must start and end "
                    << "on conductors or superconducting films!");
  }
  std::vector<double> work, kinetic;
  mfem::DenseMatrix G =
      SourceEnergyMatrix(ids, src, fields, n_fields, &work, fun.l.empty() ? nullptr : &fun,
                         region_fields, volume_energy ? &kinetic : nullptr);

  // Surface currents: G_ij = u_i^T L u_j - J_i^T u_j - J_j^T u_i = -E_ij at the solution,
  // stationary in the regularization. The part of a current that its quadrature leaves
  // incompatible with the null space is amplified by 1/eps, and twice its regularization
  // energy J^T u - u^T L u bounds the error of the energy.
  for (int k = 0; k < nc; k++)
  {
    const double uLu = G(k, k) + 2.0 * work[k], uDu = work[k] - uLu;
    MFEM_VERIFY(uDu <= Impl::kMaxRegularizationFraction * std::abs(uLu),
                "Surface current excitation "
                    << current_ids[k] << " closes too inexactly for the regularization of "
                    << "the solve (regularization energy " << uDu / std::abs(uLu)
                    << " of the magnetic energy): refine its mesh!");
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
  if (volume_energy)
  {
    volume_energy->resize(nc + nf);
    for (int k = 0; k < nc + nf; k++)
    {
      (*volume_energy)[k] = 0.5 * (E(k, k) - kinetic[k]);
    }
  }
  return E;
}

mfem::DenseMatrix SubstructuringSolver::SourceEnergyMatrix(
    const std::vector<int> &ids, const std::vector<Source> &src,
    std::vector<Vector> *fields, int n_fields, std::vector<double> *work, Functionals *fun,
    std::vector<Vector> *region_fields, std::vector<double> *kinetic)
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
    const std::size_t nn = static_cast<std::size_t>(n) * n,
                      nk = nn + n + static_cast<std::size_t>(nl) * n;
    std::vector<double> loc(nk + n, 0.0);
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
      loc[nk + i] = d[i] * Md[i];
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
    if (kinetic)
    {
      kinetic->assign(loc.begin() + nk, loc.end());
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
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.model_file.empty())
    {
      impl->AppendSourceModes(subcfg.model_file);
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
  auto hkin_of = [&](int i, int j)
  { return col[j] < 0 ? 0.0 : impl->SHkin_rows[row(i) + col[j]]; };
  auto pair_of = [&](const std::vector<double> &C, int i, int j)
  {
    return (col[i] < 0 || col[j] < 0) ? 0.0
                                      : C[static_cast<std::size_t>(col[i]) * K + col[j]];
  };

  // Region-condensed solve per source (environment interior left at zero): the region
  // source plus the environment's interface load g_k.
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
  if (kinetic)
  {
    // (u_k - a_k)^T M_sheet^R (u_k - a_k) in the region, and in the environment u_k,G^T
    // S^Kin u_k,G + 2 u_k,G^T h^kin_k + c^kin_kk.
    kinetic->assign(n, 0.0);
    Vector u(nt), Qu(nt), Qug;
    u = 0.0;
    for (int k = 0; k < n; k++)
    {
      double e = dc[k] * Mdc[k];
      if (impl->mat_dtn_Kin)
      {
        u.SetSubVector(gm, ug[k]);  // S^Kin reads the interface entries only
        impl->mat_dtn_Kin->Mult(u, Qu);
        Qu.GetSubVector(gm, Qug);
        e += ug[k] * Qug;
        for (int q = 0; q < gm.Size(); q++)
        {
          e += 2.0 * ug[k](q) * hkin_of(gm[q], k);
        }
      }
      (*kinetic)[k] = e;
    }
    MPI_Allreduce(MPI_IN_PLACE, kinetic->data(), n, MPI_DOUBLE, MPI_SUM, comm);
    if (impl->mat_dtn_Kin)
    {
      for (int k = 0; k < n; k++)
      {
        (*kinetic)[k] += pair_of(impl->SCkin, k, k);
      }
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

mfem::DenseMatrix SubstructuringSolver::EnergyMatrixDbc(const std::vector<int> &ids,
                                                        const std::vector<Vector> &xd,
                                                        std::vector<Vector> *fields,
                                                        int n_fields,
                                                        std::vector<Vector> *region_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  MFEM_VERIFY(!impl->magnetostatic,
              "Dirichlet-lift energies are for electrostatic substructuring problems!");
  const int nt = impl->nt;
  const int n = static_cast<int>(xd.size());
  const int nf = fields ? std::max(0, std::min(n_fields, n)) : 0;
  MPI_Comm comm = impl->parent_fes.GetComm();

  // Modes of the lifts that reach the environment: reuse the saved / current ones if they
  // match (id + environment fingerprint), else compute them (needs the environment once;
  // appended to the model in an offline run that saves one). A lift that does not reach the
  // environment, such as a terminal inside the region, has zero modes.
  std::vector<int> col;
  std::vector<char> need;
  if (!impl->ModesMatch(ids, xd, col, need))
  {
    std::vector<int> req_ids;
    std::vector<Vector> req_lifts;
    for (int j = 0; j < n; j++)
    {
      col[j] = need[j] ? static_cast<int>(req_ids.size()) : -1;
      if (need[j])
      {
        req_ids.push_back(ids[j]);
        req_lifts.push_back(impl->DirichletField(xd[j]));
      }
    }
    impl->ComputeModes(req_ids, req_lifts);
    const auto &subcfg = *impl->iodata.solver.substructuring;
    if (subcfg.mode != SubstructuringMode::ONLINE && !subcfg.model_file.empty())
    {
      impl->AppendModes(subcfg.model_file);
    }
  }
  const int K = static_cast<int>(impl->mode_ids.size());
  auto row = [&](int i)
  { return static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K; };
  // Modes of lift j (zero without an environment part).
  auto g_of = [&](int i, int j)
  { return col[j] < 0 ? 0.0 : impl->G_rows[row(i) + col[j]]; };

  // Region-condensed solve per lift (environment interior left at zero): the interface load
  // is the precomputed g_k, so no environment solve is needed.
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
        trr[i] += g_of(i, j);
      }
      tr.Neg();
      tr.GetSubVector(impl->rf_idx, xr[j]);
    }
  }
  impl->SolveRegion(xr);

  // Per lift, on the region elements' DOFs: u_j and K_R u_j; on Gamma: u_j and S_E u_j.
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
      impl->mat_dtn->Mult(u, Su);
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

  // E_ij = u_i^T K_R u_j + [u_G,i; e_i]^T [S_E, G; G^T, Cmode] [u_G,j; e_j]: the region
  // energy plus the environment's exact condensed energy.
  mfem::DenseMatrix E(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double a = (uc[i] * Kuc[j]) + (ug[i] * Sug[j]);
      for (int q = 0; q < gm.Size(); q++)
      {
        a += ug[i](q) * g_of(gm[q], j) + g_of(gm[q], i) * ug[j](q);
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
      if (col[i] >= 0 && col[j] >= 0)
      {
        E(i, j) += impl->Cmode[static_cast<std::size_t>(col[i]) * K + col[j]];
      }
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
                                           const std::vector<double> &energies) const
{
  BlockTimer bt(Timer::ESTIMATION);
  MFEM_VERIFY(energies.size() == region_fields.size(),
              "RegionErrorIndicator: one energy per field!");
  const auto &iodata = impl->iodata;
  const int order = iodata.solver.order, dim = impl->parent.Dimension();
  const double tol = iodata.solver.linear.estimator_tol;
  const int max_it = iodata.solver.linear.estimator_max_it;
  const auto coarsening = iodata.solver.linear.mg_coarsening;

  // Palace's flux estimators on the region submesh: their flux recovery sees only the
  // region, so the unknown environment field does not pollute the indicators near the
  // interface, and is discontinuous across the sheets, as in the regular solves.
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
      estimator.AddErrorIndicator(f, energies[k], region);
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
                                             0, false,
                                             iodata.boundaries.GetSheetAttributes());
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
                                             0, false,
                                             iodata.boundaries.GetSheetAttributes());
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

double SubstructuringSolver::MutualSystemEnergy(const Vector &ui, const Vector &uj) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->A_region->Mult(uj, t);
  impl->Aenv().Mult(uj, t2);
  t += t2;
  return mfem::InnerProduct(impl->parent_fes.GetComm(), ui, t);
}

long long int SubstructuringSolver::GlobalTrueVSize() const
{
  return impl->SolutionSpace().GlobalTrueVSize();
}

const FiniteElementSpaceHierarchy &SubstructuringSolver::GetFESpaces() const
{
  return impl->magnetostatic ? impl->curlcurl_op->GetNDSpaces()
                             : impl->laplace_op->GetH1Spaces();
}

void SubstructuringSolver::WriteParaView(const std::string &dir,
                                         const std::vector<int> &ids,
                                         const std::vector<Vector> &fields) const
{
  BlockTimer bt(Timer::POSTPRO_PARAVIEW);
  MFEM_VERIFY(ids.size() == fields.size(), "WriteParaView requires one field per id!");
  const int order = impl->iodata.solver.order;
  const auto &units = impl->iodata.units;
  GridFunction phi(impl->parent_fes);
  mfem::ParaViewDataCollection pv("paraview", &impl->parent);
  pv.SetPrefixPath(dir);
  pv.SetLevelsOfDetail(order);
  pv.SetHighOrderOutput(order > 1);
  pv.SetDataFormat(mfem::VTKFormat::BINARY);
  pv.RegisterField(impl->magnetostatic ? "A" : "V", &phi.Real());
  // The mesh and the fields in the units of the input, as PostOperator writes them.
  const double mesh_Lc0 = units.GetMeshLengthRelativeScale();
  mesh::DimensionalizeMesh(impl->parent, mesh_Lc0);
  mfem::Vector u0;
  for (std::size_t j = 0; j < fields.size(); j++)
  {
    phi.Real().SetFromTrueDofs(impl->FirstVersion(fields[j], u0));
    ScaleGridFunction(mesh_Lc0, phi);
    if (impl->magnetostatic)
    {
      units.DimensionalizeInPlace<Units::ValueType::CURRENT>(phi.Real());
    }
    else
    {
      units.DimensionalizeInPlace<Units::ValueType::VOLTAGE>(phi.Real());
    }
    phi.Real().HostRead();  // the writer reads the host copy
    pv.SetCycle(static_cast<int>(j));
    pv.SetTime(static_cast<double>(ids[j]));
    pv.Save();
  }
  mesh::NondimensionalizeMesh(impl->parent, mesh_Lc0);
}

}  // namespace palace
