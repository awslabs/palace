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
template <typename T>
class MumpsSchurSolverT
{
};
#endif

namespace
{

// BLR tolerance of the factorizations, at round-off: MUMPS's BLR factorization of these
// complex symmetric operators is several times faster than its full-rank one even
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

DrivenSubstructure::Parts
DrivenSubstructure::SideOperator(double omega, const std::vector<int> &attrs,
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
  Parts A;
  A.real = Sum({{1.0, Part(K, false)},
                {-omega, Part(C, true)},
                {-w2, Part(M, false)},
                {1.0, Part(A2.get(), false)}});
  A.imag = Sum({{1.0, Part(K, true)},
                {omega, Part(C, false)},
                {-w2, Part(M, true)},
                {1.0, Part(A2.get(), true)}});
  MFEM_VERIFY(A.real, "Missing real part of a substructure operator!");
  // Pin the given DOFs. A side-restricted operator stores no entries on the DOFs no element
  // of its side touches, so the unit diagonal is added rather than set.
  A.real->EliminateBC(pinned, Operator::DIAG_ZERO);
  if (A.imag)
  {
    A.imag->EliminateBC(pinned, Operator::DIAG_ZERO);
  }
  const auto &fes = space_op.GetNDSpace().Get();
  mfem::SparseMatrix diag(fes.GetTrueVSize());
  for (int i : pinned)
  {
    diag.Set(i, i, 1.0);
  }
  diag.Finalize();
  mfem::HypreParMatrix I(fes.GetComm(), fes.GlobalTrueVSize(),
                         const_cast<mfem::ParFiniteElementSpace &>(fes).GetTrueDofOffsets(),
                         &diag);
  A.real.reset(mfem::Add(1.0, *A.real, 1.0, I));
  return A;
}

void DrivenSubstructure::Condense(double omega)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  MPI_Comm comm = space_op.GetComm();
  const int nG = InterfaceSize();
  env_op = SideOperator(omega, env_attrs, K_env.get(), C_env.get(), M_env.get(), other_env);
  region_op = SideOperator(omega, region_attrs, K_region.get(), C_region.get(),
                           M_region.get(), other_region);
  if (!env_schur)
  {
    using Solver = MumpsSchurSolverT<std::complex<double>>;
    env_schur = std::make_unique<Solver>(*env_op.real, gamma_tdofs, kBlrTol, false, true,
                                         false, env_op.imag.get());
    region_schur = std::make_unique<Solver>(*region_op.real, gamma_tdofs, kBlrTol, false,
                                            true, false, region_op.imag.get());
  }
  else
  {
    env_schur->Refactor(*env_op.real, env_op.imag.get());
    region_schur->Refactor(*region_op.real, region_op.imag.get());
  }

  // MUMPS keeps its own copy of the entries.
  env_op = {};
  region_op = {};
  // The interface system S_R + S_E, factored on rank 0 (complex symmetric).
  if (Mpi::Root(comm))
  {
    S = env_schur->Schur();
    T = region_schur->Schur();
    S.resize(static_cast<std::size_t>(nG) * nG);
    T.resize(static_cast<std::size_t>(nG) * nG);
  }
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
  const int nG = InterfaceSize();
  auto masked = [&](const ComplexVector &b, const std::vector<char> &mask, bool gamma)
  {
    ComplexVector y(nt);
    const double *br = b.Real().HostRead(), *bi = b.Imag().HostRead();
    double *yr = y.Real().HostWrite(), *yi = y.Imag().HostWrite();
    for (int i = 0; i < nt; i++)
    {
      const bool keep = mask[i] || (gamma && is_gamma[i]);
      yr[i] = keep ? br[i] : 0.0;
      yi[i] = keep ? bi[i] : 0.0;
    }
    return y;
  };

  // Condensation of the sources onto Γ: the environment's with the interface loads b_Γ,
  // b_Γ - A_ΓE A_EE^-1 b_E, and the region's, -A_ΓR A_RR^-1 b_R.
  std::vector<ComplexVector> w(n), y(n);
  std::vector<const ComplexVector *> W(n), Y(n);
  for (int k = 0; k < n; k++)
  {
    w[k] = masked(*rhs[k], is_env_int, true);
    y[k] = masked(*rhs[k], is_region_int, false);
    W[k] = &w[k];
    Y[k] = &y[k];
  }
  std::vector<std::complex<double>> rG, rR;
  env_schur->Reduce(W, rG);
  region_schur->Reduce(Y, rR);
  if (Mpi::Root(comm) && nG > 0)
  {
    for (std::size_t q = 0; q < rG.size(); q++)
    {
      rG[q] += rR[q];
    }
    int info = 0;
    zsytrs_("L", &nG, &n, T.data(), &nG, T_piv.data(), rG.data(), &nG, &info);
    MFEM_VERIFY(info == 0, "Solve of the interface system failed: info = " << info);
  }

  // Both interiors from the interface solution, A_II^-1 (b_I - A_IΓ u_Γ).
  std::vector<ComplexVector *> Wo(n), Yo(n);
  for (int k = 0; k < n; k++)
  {
    Wo[k] = &w[k];
    Yo[k] = &y[k];
  }
  env_schur->Expand(rG, Wo);
  region_schur->Expand(rG, Yo);
  u.resize(n);
  for (int k = 0; k < n; k++)
  {
    u[k].SetSize(nt);
    u[k].UseDevice(true);
    double *ur = u[k].Real().HostWrite(), *ui = u[k].Imag().HostWrite();
    const double *wr = w[k].Real().HostRead(), *wi = w[k].Imag().HostRead();
    const double *yr = y[k].Real().HostRead(), *yi = y[k].Imag().HostRead();
    for (int i = 0; i < nt; i++)
    {
      // The interface values come out of either expansion; the Dirichlet DOFs are 0.
      const bool region = is_region_int[i];
      ur[i] = (is_env_int[i] || is_gamma[i]) ? wr[i] : (region ? yr[i] : 0.0);
      ui[i] = (is_env_int[i] || is_gamma[i]) ? wi[i] : (region ? yi[i] : 0.0);
    }
  }
#endif
}

}  // namespace palace
