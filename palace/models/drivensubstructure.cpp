// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <utility>
#include "fem/substructure.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"

extern "C"
{
  void zsytrf_(const char *, const int *, std::complex<double> *, const int *, int *,
               std::complex<double> *, const int *, int *);
  void zsytrs_(const char *, const int *, const int *, const std::complex<double> *,
               const int *, const int *, std::complex<double> *, const int *, int *);
}

namespace palace
{

#if !defined(MFEM_USE_MUMPS)
class MumpsSchurSolver
{
};
#endif

namespace
{

// BLR tolerance of the factorizations, at round-off: MUMPS's BLR factorization of these
// symmetric indefinite operators is several times faster than its full-rank one even
// without compression, at no loss of accuracy.
constexpr double kBlrTol = 1.0e-14;

// The assembled real (imag = false) or imaginary part of an operator, if any.
const mfem::HypreParMatrix *Part(const ComplexOperator *A, bool imag)
{
  const auto *part =
      A ? dynamic_cast<const ParOperator *>(imag ? A->Imag() : A->Real()) : nullptr;
  return part ? &part->ParallelAssemble() : nullptr;
}

// sum_k a_k X_k over the given (coefficient, matrix) terms (null matrices skipped).
std::unique_ptr<mfem::HypreParMatrix>
Sum(const std::vector<std::pair<double, const mfem::HypreParMatrix *>> &terms)
{
  std::unique_ptr<mfem::HypreParMatrix> sum;
  for (const auto &[a, X] : terms)
  {
    if (!X)
    {
      continue;
    }
    if (!sum)
    {
      sum = std::make_unique<mfem::HypreParMatrix>(*X);
      *sum *= a;
    }
    else
    {
      sum.reset(mfem::Add(1.0, *sum, a, *X));
    }
  }
  return sum;
}

}  // namespace

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes)
  : space_op(space_op), region_attrs(region_attributes), env_attrs(environment_attributes)
{
  const auto &fes = space_op.GetNDSpace().Get();
  const auto &mesh = *fes.GetParMesh();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int a = mesh.GetAttribute(e);
    MFEM_VERIFY(std::ranges::find(region_attrs, a) != region_attrs.end() ||
                    std::ranges::find(env_attrs, a) != env_attrs.end(),
                "Domain attribute " << a
                                    << " is neither in the region nor in the "
                                       "environment!");
  }

  // Interface DOFs: on region and environment elements; environment interior: on
  // environment elements only; region-free: on region elements (all without the Dirichlet
  // DOFs).
  mfem::Array<int> rm, em, im;
  {
    mfem::Array<int> ra(region_attrs.data(), static_cast<int>(region_attrs.size())),
        ea(env_attrs.data(), static_cast<int>(env_attrs.size()));
    MarkInterfaceTrueDofs(const_cast<mfem::ParFiniteElementSpace &>(fes), ra, ea, rm, em,
                          im);
  }
  const int nt = fes.GetTrueVSize();
  std::vector<char> dbc(nt, 0);
  for (int d : space_op.GetNDDbcTDofLists().back())
  {
    dbc[d] = 1;
  }
  is_env_int.assign(nt, 0);
  is_gamma.assign(nt, 0);
  is_region_int.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    is_gamma[i] = !dbc[i] && rm[i] && em[i];
    is_env_int[i] = !dbc[i] && em[i] && !rm[i];
    is_region_int[i] = !dbc[i] && rm[i] && !em[i];
    if (!is_gamma[i] && !is_env_int[i])
    {
      other_env.Append(i);
    }
    if (!is_gamma[i] && !is_region_int[i])
    {
      other_region.Append(i);
    }
  }

  // Interface order: by owner rank, then local index.
  std::vector<HYPRE_BigInt> mine;
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      mine.push_back(tstart + i);
    }
  }
  MPI_Comm comm = fes.GetComm();
  const int nranks = Mpi::Size(comm), nloc = static_cast<int>(mine.size());
  gamma_cnt.assign(nranks, 0);
  gamma_disp.assign(nranks, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, gamma_cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    gamma_disp[r] = gamma_disp[r - 1] + gamma_cnt[r - 1];
  }
  gamma_tdofs.resize(gamma_disp[nranks - 1] + gamma_cnt[nranks - 1]);
  MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, gamma_tdofs.data(), gamma_cnt.data(),
                 gamma_disp.data(), HYPRE_MPI_BIG_INT, comm);

  // The frequency-independent operators of each side.
  {
    SpaceOperator::AssemblyRestriction environment(space_op, env_attrs);
    K_env = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    C_env = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    M_env = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }
  {
    SpaceOperator::AssemblyRestriction region(space_op, region_attrs);
    K_region = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    C_region = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    M_region = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }
}

DrivenSubstructure::~DrivenSubstructure() = default;

std::unique_ptr<mfem::HypreParMatrix>
DrivenSubstructure::BlockOperator(double omega, const std::vector<int> &attrs,
                                  const ComplexOperator *K, const ComplexOperator *C,
                                  const ComplexOperator *M, const mfem::Array<int> &pinned)
{
  // K + iω C - ω² M + A2(ω): real part Kr - ω Ci - ω² Mr + A2r, imaginary part
  // Ki + ω Cr - ω² Mi + A2i.
  std::unique_ptr<ComplexOperator> A2;
  {
    SpaceOperator::AssemblyRestriction side(space_op, attrs);
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  }
  const double w2 = omega * omega;
  auto Ar = Sum({{1.0, Part(K, false)},
                 {-omega, Part(C, true)},
                 {-w2, Part(M, false)},
                 {1.0, Part(A2.get(), false)}});
  auto Ai = Sum({{1.0, Part(K, true)},
                 {omega, Part(C, false)},
                 {-w2, Part(M, true)},
                 {1.0, Part(A2.get(), true)}});
  MFEM_VERIFY(Ar, "Missing real part of a substructure operator!");
  if (!Ai)
  {
    // A lossless side: a zero imaginary part with the pattern of the real one.
    Ai = std::make_unique<mfem::HypreParMatrix>(*Ar);
    *Ai = 0.0;
  }
  // Pin the given DOFs. A side-restricted operator stores no entries on the DOFs no element
  // of its side touches, so the unit diagonal is added rather than set.
  Ar->EliminateBC(pinned, Operator::DIAG_ZERO);
  Ai->EliminateBC(pinned, Operator::DIAG_ZERO);
  {
    const auto &fes = space_op.GetNDSpace().Get();
    mfem::SparseMatrix diag(fes.GetTrueVSize());
    for (int i : pinned)
    {
      diag.Set(i, i, 1.0);
    }
    diag.Finalize();
    mfem::HypreParMatrix I(
        fes.GetComm(), fes.GlobalTrueVSize(),
        const_cast<mfem::ParFiniteElementSpace &>(fes).GetTrueDofOffsets(), &diag);
    Ar.reset(mfem::Add(1.0, *Ar, 1.0, I));
  }
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  mfem::Array2D<double> coeffs(2, 2);
  blocks(0, 0) = Ar.get();
  blocks(0, 1) = Ai.get();
  blocks(1, 0) = Ai.get();
  blocks(1, 1) = Ar.get();
  coeffs(0, 0) = coeffs(0, 1) = coeffs(1, 0) = 1.0;
  coeffs(1, 1) = -1.0;
  return std::unique_ptr<mfem::HypreParMatrix>(
      mfem::HypreParMatrixFromBlocks(blocks, &coeffs));
}

std::vector<std::complex<double>>
DrivenSubstructure::ComplexSchur(const MumpsSchurSolver &schur) const
{
  // The Schur complement of [[Ar, Ai], [Ai, -Ar]] is [[Sr, Si], [Si, -Sr]].
  std::vector<std::complex<double>> X;
#if defined(MFEM_USE_MUMPS)
  if (Mpi::Root(space_op.GetComm()))
  {
    const int nG = InterfaceSize();
    const auto &SB = schur.Schur();
    const std::size_t ld = 2 * static_cast<std::size_t>(nG);
    X.resize(static_cast<std::size_t>(nG) * nG);
    for (int b = 0; b < nG; b++)
    {
      for (int a = 0; a < nG; a++)
      {
        X[static_cast<std::size_t>(b) * nG + a] = {SB[b * ld + a], SB[(nG + b) * ld + a]};
      }
    }
  }
#endif
  return X;
}

void DrivenSubstructure::Condense(double omega)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MPI_Comm comm = space_op.GetComm();
  const int nG = InterfaceSize();
  env_op =
      BlockOperator(omega, env_attrs, K_env.get(), C_env.get(), M_env.get(), other_env);
  region_op = BlockOperator(omega, region_attrs, K_region.get(), C_region.get(),
                            M_region.get(), other_region);
  if (!env_schur)
  {
    // Schur variables of the real forms: the real, then the imaginary interface parts.
    // Their local rows are the local real rows, then the local imaginary rows.
    const auto &fes = space_op.GetNDSpace().Get();
    const HYPRE_BigInt tstart = fes.GetMyTDofOffset(), bstart = env_op->GetRowStarts()[0];
    const int nt = fes.GetTrueVSize();
    schur_vars.assign(2 * nG, -1);
    for (int g = 0; g < nG; g++)
    {
      const HYPRE_BigInt i = gamma_tdofs[g] - tstart;  // local only on its owner
      if (i >= 0 && i < nt)
      {
        schur_vars[g] = bstart + i;
        schur_vars[nG + g] = bstart + nt + i;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, schur_vars.data(), 2 * nG, HYPRE_MPI_BIG_INT, MPI_MAX,
                  comm);
    env_schur = std::make_unique<MumpsSchurSolver>(*env_op, schur_vars, kBlrTol, false,
                                                   true, false);
    region_schur = std::make_unique<MumpsSchurSolver>(*region_op, schur_vars, kBlrTol,
                                                      false, true, false);
  }
  else
  {
    env_schur->Refactor(*env_op);
    region_schur->Refactor(*region_op);
  }

  // The interface system S_R + S_E, factored on rank 0 (complex symmetric).
  S = ComplexSchur(*env_schur);
  T = ComplexSchur(*region_schur);
  if (Mpi::Root(comm))
  {
    for (std::size_t k = 0; k < T.size(); k++)
    {
      T[k] += S[k];
    }
    T_piv.resize(nG);
    int info = 0, lwork = -1;
    std::complex<double> wq;
    zsytrf_("L", &nG, T.data(), &nG, T_piv.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq.real()));
    std::vector<std::complex<double>> work(lwork);
    zsytrf_("L", &nG, T.data(), &nG, T_piv.data(), work.data(), &lwork, &info);
    MFEM_VERIFY(info == 0, "Factorization of the interface system failed: info = " << info);
  }
#endif
}

void DrivenSubstructure::Solve(const std::vector<const ComplexVector *> &rhs,
                               std::vector<ComplexVector> &u)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MFEM_VERIFY(env_schur && region_schur, "Condense must be called before Solve!");
  MPI_Comm comm = space_op.GetComm();
  const int n = static_cast<int>(rhs.size()), nt = space_op.GetNDSpace().GetTrueVSize();
  const int nG = InterfaceSize(), nloc = gamma_cnt[Mpi::Rank(comm)];
  // Vectors of the real forms [Re(x); -Im(x)] act as x; right-hand sides are
  // [Re(b); Im(b)]. Local layout: the real rows, then the imaginary rows.
  auto rhs_block = [&](const ComplexVector &b, const std::vector<char> &mask)
  {
    mfem::Vector y(2 * nt);
    const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
    for (int i = 0; i < nt; i++)
    {
      y(i) = mask[i] ? br[i] : 0.0;
      y(nt + i) = mask[i] ? bi[i] : 0.0;
    }
    return y;
  };
  auto solve = [](MumpsSchurSolver &lu, std::vector<mfem::Vector> &x)
  {
    std::vector<const mfem::Vector *> X(x.size());
    std::vector<mfem::Vector *> Y(x.size());
    for (std::size_t k = 0; k < x.size(); k++)
    {
      X[k] = Y[k] = &x[k];
    }
    lu.SolveInternal(X, Y);
  };
  mfem::Vector t(2 * nt);
  auto subtract_rows = [&](const mfem::HypreParMatrix &A, const mfem::Vector &x,
                           const std::vector<char> &mask, mfem::Vector &y)
  {
    A.Mult(x, t);
    for (int i = 0; i < nt; i++)
    {
      if (mask[i])
      {
        y(i) -= t(i);
        y(nt + i) -= t(nt + i);
      }
    }
  };

  // Interior solves with the interface held at zero: environment and region sources.
  std::vector<mfem::Vector> w(n), y(n), r(n);
  for (int k = 0; k < n; k++)
  {
    w[k] = rhs_block(*rhs[k], is_env_int);
    y[k] = rhs_block(*rhs[k], is_region_int);
  }
  solve(*env_schur, w);
  solve(*region_schur, y);

  // Interface right-hand sides b_Γ - A_ΓE w - A_ΓR y, gathered on rank 0 (complex).
  std::vector<std::complex<double>> rG(Mpi::Root(comm) ? static_cast<std::size_t>(nG) * n
                                                       : 0),
      mine(nloc);
  for (int k = 0; k < n; k++)
  {
    r[k] = rhs_block(*rhs[k], is_gamma);
    subtract_rows(*env_op, w[k], is_gamma, r[k]);
    subtract_rows(*region_op, y[k], is_gamma, r[k]);
    for (int i = 0, g = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        mine[g++] = {r[k](i), r[k](nt + i)};
      }
    }
    MPI_Gatherv(mine.data(), nloc, MPI_C_DOUBLE_COMPLEX,
                Mpi::Root(comm) ? rG.data() + static_cast<std::size_t>(k) * nG : nullptr,
                gamma_cnt.data(), gamma_disp.data(), MPI_C_DOUBLE_COMPLEX, 0, comm);
  }
  if (Mpi::Root(comm) && nG > 0)
  {
    int info = 0;
    zsytrs_("L", &nG, &n, T.data(), &nG, T_piv.data(), rG.data(), &nG, &info);
    MFEM_VERIFY(info == 0, "Solve of the interface system failed: info = " << info);
  }

  // The interface solution in real form, then both interiors: A_II^-1 (b_I - A_IΓ u_Γ).
  for (int k = 0; k < n; k++)
  {
    MPI_Scatterv(Mpi::Root(comm) ? rG.data() + static_cast<std::size_t>(k) * nG : nullptr,
                 gamma_cnt.data(), gamma_disp.data(), MPI_C_DOUBLE_COMPLEX, mine.data(),
                 nloc, MPI_C_DOUBLE_COMPLEX, 0, comm);
    r[k] = 0.0;
    for (int i = 0, g = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        r[k](i) = mine[g].real();
        r[k](nt + i) = -mine[g].imag();
        g++;
      }
    }
    w[k] = rhs_block(*rhs[k], is_env_int);
    y[k] = rhs_block(*rhs[k], is_region_int);
    subtract_rows(*env_op, r[k], is_env_int, w[k]);
    subtract_rows(*region_op, r[k], is_region_int, y[k]);
  }
  solve(*env_schur, w);
  solve(*region_schur, y);
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    double *ur = u[k].Real().HostWrite(), *ui = u[k].Imag().HostWrite();
    for (int i = 0; i < nt; i++)
    {
      const mfem::Vector &x = is_env_int[i] ? w[k] : (is_region_int[i] ? y[k] : r[k]);
      ur[i] = x(i);
      ui[i] = -x(nt + i);
    }
  }
#endif
}

}  // namespace palace
