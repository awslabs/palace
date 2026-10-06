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

namespace palace
{

#if !defined(MFEM_USE_MUMPS)
class MumpsSchurSolver
{
};
#endif

namespace
{

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
  is_region_free.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    is_gamma[i] = !dbc[i] && rm[i] && em[i];
    is_env_int[i] = !dbc[i] && em[i] && !rm[i];
    is_region_free[i] = !dbc[i] && rm[i];
    if (!is_gamma[i] && !is_env_int[i])
    {
      other_env.Append(i);
    }
    if (!is_region_free[i])
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
  const int nranks = Mpi::Size(comm);
  gamma_nloc = static_cast<int>(mine.size());
  std::vector<int> cnt(nranks), disp(nranks, 0);
  MPI_Allgather(&gamma_nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    disp[r] = disp[r - 1] + cnt[r - 1];
  }
  gamma_off = disp[Mpi::Rank(comm)];
  gamma_tdofs.resize(disp[nranks - 1] + cnt[nranks - 1]);
  MPI_Allgatherv(mine.data(), gamma_nloc, HYPRE_MPI_BIG_INT, gamma_tdofs.data(), cnt.data(),
                 disp.data(), HYPRE_MPI_BIG_INT, comm);

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

std::unique_ptr<mfem::HypreParMatrix> DrivenSubstructure::BlockOperator(
    double omega, const std::vector<int> &attrs, const ComplexOperator *K,
    const ComplexOperator *C, const ComplexOperator *M, const mfem::Array<int> &pinned,
    const mfem::HypreParMatrix *Xr, const mfem::HypreParMatrix *Xi)
{
  // K + iω C - ω² M + A2(ω) (+ X): real part Kr - ω Ci - ω² Mr + A2r (+ Xr), imaginary part
  // Ki + ω Cr - ω² Mi + A2i (+ Xi).
  std::unique_ptr<ComplexOperator> A2;
  {
    SpaceOperator::AssemblyRestriction side(space_op, attrs);
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  }
  const double w2 = omega * omega;
  auto Ar = Sum({{1.0, Part(K, false)},
                 {-omega, Part(C, true)},
                 {-w2, Part(M, false)},
                 {1.0, Part(A2.get(), false)},
                 {1.0, Xr}});
  auto Ai = Sum({{1.0, Part(K, true)},
                 {omega, Part(C, false)},
                 {-w2, Part(M, true)},
                 {1.0, Part(A2.get(), true)},
                 {1.0, Xi}});
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

std::pair<std::unique_ptr<mfem::HypreParMatrix>, std::unique_ptr<mfem::HypreParMatrix>>
DrivenSubstructure::InterfaceMatrices() const
{
  const auto &fes = space_op.GetNDSpace().Get();
  const int nt = fes.GetTrueVSize(), nG = InterfaceSize();
  std::vector<int> I(nt + 1, 0);
  std::vector<HYPRE_BigInt> J;
  std::vector<double> Vr, Vi;
  J.reserve(static_cast<std::size_t>(gamma_nloc) * nG);
  Vr.reserve(J.capacity());
  Vi.reserve(J.capacity());
  for (int i = 0, g = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      for (int j = 0; j < nG; j++)
      {
        const auto v = S_rows[static_cast<std::size_t>(g) * nG + j];
        J.push_back(gamma_tdofs[j]);
        Vr.push_back(v.real());
        Vi.push_back(v.imag());
      }
      g++;
    }
    I[i + 1] = static_cast<int>(J.size());
  }
  const HYPRE_BigInt glob = fes.GlobalTrueVSize();
  auto *offsets = const_cast<mfem::ParFiniteElementSpace &>(fes).GetTrueDofOffsets();
  auto make = [&](std::vector<double> &V)
  {
    return std::make_unique<mfem::HypreParMatrix>(fes.GetComm(), nt, glob, glob, I.data(),
                                                  J.data(), V.data(), offsets, offsets);
  };
  return {make(Vr), make(Vi)};
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
  if (!env_schur)
  {
    // Schur variables of the real form: the real, then the imaginary interface parts. Its
    // local rows are the local real rows, then the local imaginary rows.
    const auto &fes = space_op.GetNDSpace().Get();
    const HYPRE_BigInt tstart = fes.GetMyTDofOffset(), bstart = env_op->GetRowStarts()[0];
    const int nt = fes.GetTrueVSize();
    std::vector<HYPRE_BigInt> vars(2 * nG);
    for (int g = 0; g < nG; g++)
    {
      const HYPRE_BigInt i = gamma_tdofs[g] - tstart;  // local only on its owner
      vars[g] = vars[nG + g] = -1;
      if (i >= 0 && i < nt)
      {
        vars[g] = bstart + i;
        vars[nG + g] = bstart + nt + i;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, vars.data(), 2 * nG, HYPRE_MPI_BIG_INT, MPI_MAX, comm);
    env_schur = std::make_unique<MumpsSchurSolver>(*env_op, vars, 0.0, false, true);
  }
  else
  {
    env_schur->Refactor(*env_op);
  }

  // The Schur complement of [[Ar, Ai], [Ai, -Ar]] is [[Sr, Si], [Si, -Sr]].
  S.clear();
  if (Mpi::Root(comm))
  {
    const auto &SB = env_schur->Schur();
    const std::size_t ld = 2 * static_cast<std::size_t>(nG);
    S.resize(static_cast<std::size_t>(nG) * nG);
    for (int b = 0; b < nG; b++)
    {
      for (int a = 0; a < nG; a++)
      {
        S[static_cast<std::size_t>(b) * nG + a] = {SB[b * ld + a], SB[(nG + b) * ld + a]};
      }
    }
  }

  // Each rank's interface rows of S_E (symmetric: the columns of rank 0's column-major
  // array), then the region operator with S_E on Γ.
  {
    const int nranks = Mpi::Size(comm);
    std::vector<int> cnt(nranks), disp(nranks);
    const int mine = gamma_nloc * nG;
    MPI_Allgather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 0, d = 0; r < nranks; d += cnt[r], r++)
    {
      disp[r] = d;
    }
    S_rows.resize(mine);
    MPI_Scatterv(S.data(), cnt.data(), disp.data(), MPI_C_DOUBLE_COMPLEX, S_rows.data(),
                 mine, MPI_C_DOUBLE_COMPLEX, 0, comm);
  }
  auto [Sr, Si] = InterfaceMatrices();
  auto region_op = BlockOperator(omega, region_attrs, K_region.get(), C_region.get(),
                                 M_region.get(), other_region, Sr.get(), Si.get());
  if (!region_lu)
  {
    region_lu = std::make_unique<MumpsSchurSolver>(*region_op, std::vector<HYPRE_BigInt>{},
                                                   0.0, false, true);
  }
  else
  {
    region_lu->Refactor(*region_op);
  }
#endif
}

void DrivenSubstructure::Solve(const std::vector<const ComplexVector *> &rhs,
                               std::vector<ComplexVector> &u)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MFEM_VERIFY(env_op && region_lu, "Condense must be called before Solve!");
  const int n = static_cast<int>(rhs.size()), nt = space_op.GetNDSpace().GetTrueVSize();
  // Vectors of the real form [Re(x); -Im(x)] act as x; right-hand sides are
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

  // Environment sources b_E: interface loads -(A_E A_EE^-1 b_E)|_Γ.
  std::vector<mfem::Vector> w(n), x(n);
  for (int k = 0; k < n; k++)
  {
    w[k] = rhs_block(*rhs[k], is_env_int);
  }
  solve(*env_schur, w);
  mfem::Vector t(2 * nt);
  for (int k = 0; k < n; k++)
  {
    x[k] = rhs_block(*rhs[k], is_region_free);
    env_op->Mult(w[k], t);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        x[k](i) -= t(i);
        x[k](nt + i) -= t(nt + i);
      }
    }
  }
  solve(*region_lu, x);

  // Environment interior: A_EE^-1 (b_E - A_EΓ u_Γ).
  mfem::Vector xg(2 * nt);
  for (int k = 0; k < n; k++)
  {
    xg = 0.0;
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        xg(i) = x[k](i);
        xg(nt + i) = x[k](nt + i);
      }
    }
    env_op->Mult(xg, t);
    w[k] = rhs_block(*rhs[k], is_env_int);
    for (int i = 0; i < nt; i++)
    {
      if (is_env_int[i])
      {
        w[k](i) -= t(i);
        w[k](nt + i) -= t(nt + i);
      }
    }
  }
  solve(*env_schur, w);
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    double *ur = u[k].Real().HostWrite(), *ui = u[k].Imag().HostWrite();
    for (int i = 0; i < nt; i++)
    {
      const bool env = is_env_int[i];
      ur[i] = env ? w[k](i) : x[k](i);
      ui[i] = -(env ? w[k](nt + i) : x[k](nt + i));
    }
  }
#endif
}

}  // namespace palace
