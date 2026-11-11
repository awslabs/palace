// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "rap.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <typeinfo>
#include <vector>

#include "fem/bilinearform.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/hypre.hpp"
#include "utils/communication.hpp"
#include "utils/perfprobes.hpp"
#include "utils/perfstats.hpp"

namespace palace
{

namespace
{

// Inclusive tallies for the parallel prolongation P and its transpose in forward operator
// applications (local work plus the halo exchange inside them).
int HaloTallyP()
{
  static const int id = perf::RegisterPhaseTally("Halo/P");
  return id;
}

int HaloTallyPT()
{
  static const int id = perf::RegisterPhaseTally("Halo/PT");
  return id;
}

const hypre::MaskedProlongation *GetMasked(const hypre::MaskedProlongation *&masked_P,
                                           const FiniteElementSpace &trial_fespace,
                                           const FiniteElementSpace &test_fespace,
                                           bool use_R,
                                           const mfem::Array<int> &dbc_tdof_list)
{
  if (&trial_fespace != &test_fespace || use_R)
  {
    return nullptr;
  }
  if (!masked_P ||
      static_cast<const Operator *>(&masked_P->GetProlongation()) !=
          trial_fespace.GetProlongationMatrix() ||
      !masked_P->Matches(dbc_tdof_list))
  {
    masked_P = trial_fespace.GetMaskedProlongation(dbc_tdof_list);
  }
  return masked_P;
}

// Apply the prolongation (or its transpose) of a space to the real and imaginary parts of a
// complex vector, at once when possible.
void ProlongationMult(const FiniteElementSpace &fespace, const Vector &xr, const Vector &xi,
                      Vector &yr, Vector &yi)
{
  const Operator *P = fespace.GetProlongationMatrix();
  if (const auto *fP = dynamic_cast<const hypre::HypreProlongation *>(P))
  {
    fP->Mult(xr, xi, yr, yi);
    return;
  }
  if (const auto *cP = dynamic_cast<const ConformingProlongation *>(P))
  {
    cP->Mult(xr, xi, yr, yi);
    return;
  }
  P->Mult(xr, yr);
  P->Mult(xi, yi);
}

void ProlongationMultTranspose(const FiniteElementSpace &fespace, const Vector &xr,
                               const Vector &xi, Vector &yr, Vector &yi)
{
  const Operator *P = fespace.GetProlongationMatrix();
  if (const auto *fP = dynamic_cast<const hypre::HypreProlongation *>(P))
  {
    fP->MultTranspose(xr, xi, yr, yi);
    return;
  }
  if (const auto *cP = dynamic_cast<const ConformingProlongation *>(P))
  {
    cP->MultTranspose(xr, xi, yr, yi);
    return;
  }
  P->MultTranspose(xr, yr);
  P->MultTranspose(xi, yi);
}

// A += B for an assembled square operator, keeping the diagonal entry first in each row.
void AddToAssembled(std::unique_ptr<mfem::HypreParMatrix> &A, const mfem::HypreParMatrix &B)
{
  A.reset(mfem::ParAdd(A.get(), &B));
  hypre_ParCSRMatrixSetNumNonzeros(*A);
  hypre_CSRMatrixReorder(hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)*A));
}

// The diagonal of a local operator in L-vector coordinates, and its entries between coupled
// L-dofs (see ceed::AssembleExactDiagonal), for libCEED operators and assembled local
// matrices (whose diagonal is always exact). Returns false for other operators.
bool AssembleLocalExactDiagonal(const Operator &A, const std::vector<char> *coupled,
                                bool exact_local, Vector &d,
                                const std::function<void(int, int, double)> &add_entry)
{
  if (const auto *cA = dynamic_cast<const ceed::Operator *>(&A))
  {
    return ceed::AssembleExactDiagonal(*cA, coupled, exact_local, d, add_entry);
  }
  if (const auto *sA = dynamic_cast<const hypre::HypreCSRMatrix *>(&A))
  {
    hypre_CSRMatrix *mat = *sA, *host = nullptr;
    if (hypre_CSRMatrixMemoryLocation(mat) != HYPRE_MEMORY_HOST)
    {
      host = hypre_CSRMatrixClone_v2(mat, 1, HYPRE_MEMORY_HOST);
      mat = host;
    }
    const HYPRE_Int n = hypre_CSRMatrixNumRows(mat);
    const HYPRE_Int *I = hypre_CSRMatrixI(mat), *J = hypre_CSRMatrixJ(mat);
    const double *a = hypre_CSRMatrixData(mat);
    d.SetSize(n);
    double *dd = d.HostWrite();
    for (HYPRE_Int i = 0; i < n; i++)
    {
      dd[i] = 0.0;
      for (HYPRE_Int k = I[i]; k < I[i + 1]; k++)
      {
        const HYPRE_Int j = J[k];
        if (j == i)
        {
          dd[i] += a[k];
        }
        else if (coupled && (*coupled)[i] && (*coupled)[j])
        {
          add_entry(i, j, a[k]);
        }
      }
    }
    if (host)
    {
      hypre_CSRMatrixDestroy(host);
    }
    return true;
  }
  return false;
}

// The exact diagonal of Pᵀ A P for a prolongation P assembled as a HypreParMatrix (as on
// nonconforming meshes), whose rows (L-dofs) can have several entries and whose columns
// (true dofs) can be shared by several rows of a process. With
//   (Pᵀ A P)ⱼⱼ = Σₖ Qₖⱼ,  Qₖⱼ = Pₖⱼ (Pₖⱼ Aₖₖ + Σₗ≠ₖ Aₖₗ Pₗⱼ),
// a matrix with the sparsity of P, the diagonal is Qᵀ 1 (with the communication of Pᵀ). The
// sums over l ≠ k are nonzero only for "coupled" L-dofs, whose row has several entries or
// shares a column with another row of the process, and only through the entries Aₖₗ between
// coupled L-dofs (the elements with at least two coupled L-dofs). Returns false if the
// local operator is not supported.
bool AssembleProlongatedDiagonal(const Operator &A, const mfem::HypreParMatrix &P,
                                 bool exact_local, Vector &lx, Vector &diag)
{
  P.HostRead();
  hypre_ParCSRMatrix *hP = P;
  const hypre_CSRMatrix *Pd = hypre_ParCSRMatrixDiag(hP), *Po = hypre_ParCSRMatrixOffd(hP);
  const HYPRE_Int n = hypre_CSRMatrixNumRows(Pd), nd = hypre_CSRMatrixNumCols(Pd),
                  no = hypre_CSRMatrixNumCols(Po);
  MFEM_VERIFY(n == lx.Size(), "Prolongation and local operator size mismatch!");
  const HYPRE_Int *di = hypre_CSRMatrixI(Pd), *dj = hypre_CSRMatrixJ(Pd),
                  *oi = hypre_CSRMatrixI(Po), *oj = hypre_CSRMatrixJ(Po);
  const double *da = hypre_CSRMatrixData(Pd), *oa = hypre_CSRMatrixData(Po);
  auto OffdBegin = [&](HYPRE_Int k) { return oi ? oi[k] : 0; };
  auto OffdEnd = [&](HYPRE_Int k) { return oi ? oi[k + 1] : 0; };

  // References of the columns (local true dofs, then off-diagonal ones) by rows with one
  // entry and by rows with several entries, and the coupled rows: 2 for rows with several
  // entries or sharing their column with another row with one entry, 1 for rows with one
  // entry sharing their column only with rows with several entries. Two L-dofs whose
  // supports intersect are coupled, and at least one of them is marked 2.
  auto RowSize = [&](HYPRE_Int k)
  { return (di[k + 1] - di[k]) + (OffdEnd(k) - OffdBegin(k)); };
  auto SingleColumn = [&](HYPRE_Int k)
  { return (di[k + 1] > di[k]) ? dj[di[k]] : nd + oj[OffdBegin(k)]; };
  std::vector<int> refs_single(static_cast<std::size_t>(nd) + no, 0),
      refs_multi(refs_single.size(), 0);
  for (HYPRE_Int k = 0; k < n; k++)
  {
    auto &refs = (RowSize(k) > 1) ? refs_multi : refs_single;
    for (HYPRE_Int p = di[k]; p < di[k + 1]; p++)
    {
      refs[dj[p]]++;
    }
    for (HYPRE_Int p = OffdBegin(k); p < OffdEnd(k); p++)
    {
      refs[nd + oj[p]]++;
    }
  }
  std::vector<char> coupled(n, 0);
  std::vector<int> start(n + 1, 0);
  bool any = false;
  for (HYPRE_Int k = 0; k < n; k++)
  {
    const HYPRE_Int nnz = RowSize(k);
    char c = (nnz > 1) ? 2 : 0;
    if (nnz == 1)
    {
      const HYPRE_Int j = SingleColumn(k);
      c = (refs_single[j] > 1) ? 2 : (refs_multi[j] > 0) ? 1 : 0;
    }
    coupled[k] = c;
    any = any || c;
    start[k + 1] = start[k] + (c ? nnz : 0);
  }

  // Supports of the coupled rows sorted by column, with the position of each entry in the
  // values of P (diagonal block, then off-diagonal block).
  struct Entry
  {
    int col;
    double val;
    int slot;
  };
  const HYPRE_Int nnz_d = di[n];
  std::vector<Entry> support(start[n]);
  for (HYPRE_Int k = 0; k < n; k++)
  {
    if (!coupled[k])
    {
      continue;
    }
    Entry *e = support.data() + start[k];
    for (HYPRE_Int p = di[k]; p < di[k + 1]; p++)
    {
      *e++ = {static_cast<int>(dj[p]), da[p], static_cast<int>(p)};
    }
    for (HYPRE_Int p = OffdBegin(k); p < OffdEnd(k); p++)
    {
      *e++ = {static_cast<int>(nd + oj[p]), oa[p], static_cast<int>(nnz_d + p)};
    }
    std::sort(support.data() + start[k], e,
              [](const Entry &a, const Entry &b) { return a.col < b.col; });
  }
  std::vector<double> cross(static_cast<std::size_t>(nnz_d) + (oi ? oi[n] : 0), 0.0);
  auto AddEntry = [&](int r, int c, double v)
  {
    // Σₗ≠ₖ Aₖₗ Pₗⱼ for the common columns j of rows r = k and c = l.
    const Entry *a = support.data() + start[r], *a_end = support.data() + start[r + 1];
    const Entry *b = support.data() + start[c], *b_end = support.data() + start[c + 1];
    while (a != a_end && b != b_end)
    {
      if (a->col < b->col)
      {
        a++;
      }
      else if (b->col < a->col)
      {
        b++;
      }
      else
      {
        cross[a->slot] += v * b->val;
        a++;
        b++;
      }
    }
  };
  // The product with Q is collective: fall back on every process if the local operator of
  // any is not supported.
  MPI_Comm comm = hypre_ParCSRMatrixComm(hP);
  bool supported =
      AssembleLocalExactDiagonal(A, any ? &coupled : nullptr, exact_local, lx, AddEntry);
  Mpi::GlobalAnd(1, &supported, comm);
  if (!supported)
  {
    return false;
  }
  // Q on the sparsity of P, and the diagonal Qᵀ 1.
  const double *d = lx.HostRead();
  mfem::HypreParMatrix Q(P);
  Q.HostReadWrite();
  hypre_ParCSRMatrix *hQ = Q;
  double *qd = hypre_CSRMatrixData(hypre_ParCSRMatrixDiag(hQ));
  double *qo = hypre_CSRMatrixData(hypre_ParCSRMatrixOffd(hQ));
  for (HYPRE_Int k = 0; k < n; k++)
  {
    for (HYPRE_Int p = di[k]; p < di[k + 1]; p++)
    {
      qd[p] = da[p] * (da[p] * d[k] + cross[p]);
    }
    for (HYPRE_Int p = OffdBegin(k); p < OffdEnd(k); p++)
    {
      qo[p] = oa[p] * (oa[p] * d[k] + cross[nnz_d + p]);
    }
  }
  Vector ones(n);
  ones.UseDevice(true);
  ones = 1.0;
  Q.MultTranspose(1.0, ones, 0.0, diag);
  return true;
}

}  // namespace

ParOperator::ParOperator(std::unique_ptr<Operator> &&dA, const Operator *pA,
                         const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace, bool test_restrict)
  : Operator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    data_A(std::move(dA)), A((data_A != nullptr) ? data_A.get() : pA),
    trial_fespace(trial_fespace), test_fespace(test_fespace), use_R(test_restrict),
    diag_policy(DiagonalPolicy::DIAG_ONE), RAP(nullptr)
{
  MFEM_VERIFY(A, "Cannot construct ParOperator from an empty matrix!");
  // Owned square libCEED operators may apply small simplex leaves through their element
  // matrices (CPU only, decided by the operator).
  if (auto *ceed_A = dynamic_cast<ceed::Operator *>(data_A.get()))
  {
    if (&trial_fespace == &test_fespace &&
        (typeid(*ceed_A) == typeid(ceed::Operator) ||
         typeid(*ceed_A) == typeid(ceed::SymmetricOperator)))
    {
      ceed_A->UseElementAssembly();
    }
  }
}

ParOperator::ParOperator(std::unique_ptr<Operator> &&A,
                         const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace, bool test_restrict)
  : ParOperator(std::move(A), nullptr, trial_fespace, test_fespace, test_restrict)
{
}

ParOperator::ParOperator(const Operator &A, const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace, bool test_restrict)
  : ParOperator(nullptr, &A, trial_fespace, test_fespace, test_restrict)
{
}

const hypre::MaskedProlongation *ParOperator::GetMaskedProlongation() const
{
  return GetMasked(masked_P, trial_fespace, test_fespace, use_R, dbc_tdof_list);
}

bool ParOperator::MappedMult(const Operator &local, const Vector &x, Vector &y, bool add,
                             double a) const
{
  // The same operations as the products with P, the local operator on the L-vectors and
  // Pᵀ, without forming the L-vectors (see ceed::AddMultMapped).
  const auto *mP = GetMaskedProlongation();
  if (!mP || B || x.HostRead() == y.HostRead() || !ceed::SupportsMapped(local))
  {
    return false;
  }
  const int nx = mP->NumExtra();
  mapped_x.SetSize(nx);
  mapped_y.SetSize(nx);
  auto &ty = test_fespace.GetTVector<Vector>();
  Vector &yt = add ? ty : y;
  {
    perf::ScopedTally t(HaloTallyP());
    mP->ExtraMult(x, mapped_x);
  }
  yt = 0.0;
  mapped_y = 0.0;
  ceed::AddMultMapped(local, {mP->LocalDofMap().data(), mP->Key()}, x.HostRead(),
                      mapped_x.HostRead(), yt.HostReadWrite(), mapped_y.HostReadWrite());
  {
    perf::ScopedTally t(HaloTallyPT());
    mP->TransposeComplete(yt, mapped_y);
  }
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(yt, dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(yt, dbc_tdof_list, 0.0);
    }
  }
  if (add)
  {
    y.Add(a, yt);
  }
  return true;
}

class ParOperator::SinglePrecisionOperator : public Operator
{
private:
  const ParOperator &op;

public:
  SinglePrecisionOperator(const ParOperator &op) : Operator(op.Height(), op.Width()), op(op)
  {
  }

  void Mult(const Vector &x, Vector &y) const override { op.MultWith(*op.A_single, x, y); }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    op.AddMultWith(*op.A_single, x, y, a);
  }

  void MultTranspose(const Vector &x, Vector &y) const override { op.MultTranspose(x, y); }

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    op.AddMultTranspose(x, y, a);
  }
};

std::size_t ParOperator::UseSinglePrecision()
{
  MFEM_VERIFY(&trial_fespace == &test_fespace && !use_R,
              "Single-precision approximation requires a square ParOperator!");
  std::size_t bytes = 0;
  if (!A_single && !RAP)
  {
    A_single = ceed::CreateSinglePrecisionOperator(*A, &bytes);
    if (A_single)
    {
      op_single = std::make_unique<SinglePrecisionOperator>(*this);
    }
  }
  return bytes;
}

void ParOperator::SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                                       DiagonalPolicy policy)
{
  MFEM_VERIFY(policy == DiagonalPolicy::DIAG_ONE || policy == DiagonalPolicy::DIAG_ZERO,
              "Essential boundary condition true dof elimination for ParOperator supports "
              "only DiagonalPolicy::DIAG_ONE or DiagonalPolicy::DIAG_ZERO!");
  MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                               "for rectangular ParOperator!");
  MFEM_VERIFY(!B,
              "Set essential true dofs before adding a parallel matrix to ParOperator!");
  tdof_list.Read();
  dbc_tdof_list.MakeRef(tdof_list);
  diag_policy = policy;
}

void ParOperator::AddParallelMatrix(const mfem::HypreParMatrix &B_)
{
  MFEM_VERIFY(!B && height == width && B_.Height() == height && B_.Width() == width,
              "Invalid parallel matrix for square ParOperator!");
  // The diagonal policy at the essential true dofs is kept by RAP, so B contributes zero.
  B = std::make_unique<mfem::HypreParMatrix>(B_);
  B->EliminateBC(dbc_tdof_list, DiagonalPolicy::DIAG_ZERO);
  if (RAP)
  {
    AddToAssembled(RAP, *B);
  }
}

Operator::DiagonalPolicy ParOperator::GetDiagonalPolicy() const
{
  MFEM_VERIFY(dbc_tdof_list.Size() > 0,
              "There is no DiagonalPolicy if no essential dofs have been set!");
  return diag_policy;
}

void ParOperator::EliminateRHS(const Vector &x, Vector &b) const
{
  MFEM_VERIFY(A, "No local matrix available for ParOperator::EliminateRHS!");
  MFEM_VERIFY(!B, "ParOperator::EliminateRHS does not support an added parallel matrix!");
  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  {
    auto &tx = trial_fespace.GetTVector<Vector>();
    tx = 0.0;
    linalg::SetSubVector(tx, dbc_tdof_list, x);
    trial_fespace.GetProlongationMatrix()->Mult(tx, lx);
  }

  // Apply the unconstrained operator.
  A->Mult(lx, ly);

  auto &ty = test_fespace.GetTVector<Vector>();
  RestrictionMatrixMult(ly, ty);
  b.Add(-1.0, ty);
  if (diag_policy == DiagonalPolicy::DIAG_ONE)
  {
    linalg::SetSubVector(b, dbc_tdof_list, x);
  }
  else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
  {
    linalg::SetSubVector(b, dbc_tdof_list, 0.0);
  }
}

mfem::HypreParMatrix &ParOperator::ParallelAssemble(bool skip_zeros) const
{
  if (RAP)
  {
    return *RAP;
  }
  perf::ScopedTally t(perf::AssemblyTally(perf::AssemblyKind::PARALLEL));

  // Build the square or rectangular assembled HypreParMatrix.
  const auto *sA = dynamic_cast<const hypre::HypreCSRMatrix *>(A);
  std::unique_ptr<hypre::HypreCSRMatrix> data_sA;
  if (!sA)
  {
    const auto *cA = dynamic_cast<const ceed::Operator *>(A);
    MFEM_VERIFY(cA,
                "ParOperator::ParallelAssemble requires A as an hypre::HypreCSRMatrix or "
                "ceed::Operator!");
    data_sA = BilinearForm::FullAssemble(*cA, skip_zeros, use_R);
    sA = data_sA.get();
  }

  hypre_ParCSRMatrix *hA = hypre_ParCSRMatrixCreate(
      trial_fespace.GetComm(), test_fespace.GlobalVSize(), trial_fespace.GlobalVSize(),
      test_fespace.GetDofOffsets(), trial_fespace.GetDofOffsets(), 0, sA->NNZ(), 0);
  hypre_CSRMatrix *hA_diag = hypre_ParCSRMatrixDiag(hA);
  hypre_ParCSRMatrixDiag(hA) = *const_cast<hypre::HypreCSRMatrix *>(sA);
  hypre_ParCSRMatrixInitialize(hA);

  const mfem::HypreParMatrix *P = trial_fespace.Dof_TrueDof_Matrix();
  if (!use_R)
  {
    const mfem::HypreParMatrix *Rt = test_fespace.Dof_TrueDof_Matrix();
    if (trial_fespace.IsBroken() || test_fespace.IsBroken())
    {
      // hypre_ParCSRMatrixRAPKT reads out of bounds when a process has more columns of Rt
      // (true DOFs) than rows of A (L-DOFs) and no off-diagonal columns of A P, which can
      // happen for broken spaces, where a process owns the copies of split true DOFs which
      // are read only by other processes. Compute Rtᵀ (A P) instead.
      hypre_ParCSRMatrix *AP = hypre_ParCSRMatMat(hA, *P);
      RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRTMatMatKT(*Rt, AP, 1), true);
      hypre_ParCSRMatrixDestroy(AP);
    }
    else
    {
      RAP = std::make_unique<mfem::HypreParMatrix>(
          hypre_ParCSRMatrixRAPKT(*Rt, hA, *P, 1, 1), true);
    }
  }
  else
  {
    mfem::HypreParMatrix hR(
        test_fespace.GetComm(), test_fespace.GlobalTrueVSize(), test_fespace.GlobalVSize(),
        test_fespace.GetTrueDofOffsets(), test_fespace.GetDofOffsets(),
        const_cast<mfem::SparseMatrix *>(test_fespace.GetRestrictionMatrix()));
    hypre_ParCSRMatrix *AP = hypre_ParCSRMatMat(hA, *P);
    RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatMat(hR, AP), true);
    hypre_ParCSRMatrixDestroy(AP);
  }

  hypre_ParCSRMatrixDiag(hA) = hA_diag;
  hypre_ParCSRMatrixDestroy(hA);
  hypre_ParCSRMatrixSetNumNonzeros(*RAP);
  if (&trial_fespace == &test_fespace)
  {
    // Make sure that the first entry in each row is the diagonal one, for a square matrix.
    hypre_CSRMatrixReorder(hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)*RAP));
  }

  // Eliminate boundary conditions on the assembled (square) matrix.
  if (&trial_fespace == &test_fespace)
  {
    RAP->EliminateBC(dbc_tdof_list, diag_policy);
  }
  else
  {
    MFEM_VERIFY(dbc_tdof_list.Size() == 0,
                "Essential BC elimination is only available for square ParOperator!");
  }
  if (B)
  {
    AddToAssembled(RAP, *B);
  }

  return *RAP;
}

void ParOperator::AssembleDiagonal(Vector &diag) const
{
  diag.UseDevice(true);
  if (RAP)
  {
    RAP->AssembleDiagonal(diag);
    return;
  }
  perf::ScopedTally t(perf::AssemblyTally(perf::AssemblyKind::DIAGONAL));

  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "Diagonal assembly is only available for square ParOperator!");
  auto &lx = trial_fespace.GetLVector<Vector>();
  const Operator *P = test_fespace.GetProlongationMatrix();
  const auto *hP = dynamic_cast<const mfem::HypreParMatrix *>(P);
  if (const auto *fP = dynamic_cast<const hypre::HypreProlongation *>(P))
  {
    hP = &fP->GetMatrix();
  }

  // The exact diagonal of Pᵀ A P: for a prolongation assembled as a HypreParMatrix (as on
  // nonconforming meshes) through its entries, otherwise (each L-dof a copy of one true
  // dof) as Pᵀ dₗ with the exact local diagonal dₗ.
  bool exact = false;
  if (hP)
  {
    exact = AssembleProlongatedDiagonal(*A, *hP, true, lx, diag);
  }
  else if (AssembleLocalExactDiagonal(*A, nullptr, true, lx, {}))
  {
    P->MultTranspose(lx, diag);
    exact = true;
  }
  if (!exact)
  {
    // For an AMR mesh, a convergent diagonal is assembled with |P|ᵀ dₗ, where |P| has
    // entry-wise absolute values of the conforming prolongation operator.
    A->AssembleDiagonal(lx);
    if (hP)
    {
      hP->AbsMultTranspose(1.0, lx, 0.0, diag);
    }
    else
    {
      P->MultTranspose(lx, diag);
    }
  }

  // Eliminate essential true dofs.
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(diag, dbc_tdof_list, 1.0);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(diag, dbc_tdof_list, 0.0);
    }
  }
  if (B)
  {
    // Search each row for its diagonal entry: B need not store one, let alone first.
    B->HypreRead();
    hypre_CSRMatrixExtractDiagonal(hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)*B),
                                   diag.ReadWrite(), 10);
  }
}

void ParOperator::Mult(const Vector &x, Vector &y) const
{
  MultWith(*A, x, y);
}

void ParOperator::MultWith(const Operator &local, const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::Mult!");
  if (RAP)
  {
    RAP->Mult(x, y);
    return;
  }

  if (MappedMult(local, x, y, false, 0.0))
  {
    return;
  }
  const auto *mP = GetMaskedProlongation();
  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  {
    perf::ScopedTally t(HaloTallyP());
    if (mP)
    {
      mP->Mult(x, lx);
    }
    else if (dbc_tdof_list.Size())
    {
      auto &tx = trial_fespace.GetTVector<Vector>();
      tx = x;
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
      trial_fespace.GetProlongationMatrix()->Mult(tx, lx);
    }
    else
    {
      trial_fespace.GetProlongationMatrix()->Mult(x, lx);
    }
  }

  // Apply the operator on the L-vector.
  local.Mult(lx, ly);

  {
    perf::ScopedTally t(HaloTallyPT());
    RestrictionMatrixMult(ly, y);
  }
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, dbc_tdof_list, 0.0);
    }
  }
  if (B)
  {
    B->Mult(1.0, x, 1.0, y);
  }
}

void ParOperator::MultTranspose(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ParOperator::MultTranspose!");
  if (RAP)
  {
    RAP->MultTranspose(x, y);
    return;
  }

  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<Vector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  trial_fespace.GetProlongationMatrix()->MultTranspose(lx, y);
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, dbc_tdof_list, 0.0);
    }
  }
  if (B)
  {
    B->MultTranspose(1.0, x, 1.0, y);
  }
}

void ParOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  AddMultWith(*A, x, y, a);
}

void ParOperator::AddMultWith(const Operator &local, const Vector &x, Vector &y,
                              const double a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::AddMult!");
  if (RAP)
  {
    RAP->AddMult(x, y, a);
    return;
  }

  if (MappedMult(local, x, y, true, a))
  {
    return;
  }
  const auto *mP = GetMaskedProlongation();
  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  {
    perf::ScopedTally t(HaloTallyP());
    if (mP)
    {
      mP->Mult(x, lx);
    }
    else if (dbc_tdof_list.Size())
    {
      auto &tx = trial_fespace.GetTVector<Vector>();
      tx = x;
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
      trial_fespace.GetProlongationMatrix()->Mult(tx, lx);
    }
    else
    {
      trial_fespace.GetProlongationMatrix()->Mult(x, lx);
    }
  }

  // Apply the operator on the L-vector.
  local.Mult(lx, ly);

  if (mP)
  {
    perf::ScopedTally t(HaloTallyPT());
    mP->AddMultTranspose(ly, y, a,
                         (diag_policy == DiagonalPolicy::DIAG_ONE) ? &x : nullptr);
  }
  else
  {
    auto &ty = test_fespace.GetTVector<Vector>();
    {
      perf::ScopedTally t(HaloTallyPT());
      RestrictionMatrixMult(ly, ty);
    }
    if (dbc_tdof_list.Size())
    {
      if (diag_policy == DiagonalPolicy::DIAG_ONE)
      {
        linalg::SetSubVector(ty, dbc_tdof_list, x);
      }
      else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
      {
        linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
      }
    }
    y.Add(a, ty);
  }
  if (B)
  {
    B->Mult(a, x, 1.0, y);
  }
}

void ParOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ParOperator::AddMultTranspose!");
  if (RAP)
  {
    RAP->AddMultTranspose(x, y, a);
    return;
  }

  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<Vector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  auto &tx = trial_fespace.GetTVector<Vector>();
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx, tx);
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
    }
  }
  y.Add(a, tx);
  if (B)
  {
    B->MultTranspose(a, x, 1.0, y);
  }
}

void ParOperator::RestrictionMatrixMult(const Vector &ly, Vector &ty) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->MultTranspose(ly, ty);
  }
  else
  {
    if (const auto *R = test_fespace.GetRestrictionRuns())
    {
      R->Mult(ly, ty);
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->Mult(ly, ty);
    }
  }
}

void ParOperator::RestrictionMatrixMultTranspose(const Vector &ty, Vector &ly) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->Mult(ty, ly);
  }
  else
  {
    if (const auto *R = test_fespace.GetRestrictionRuns())
    {
      R->MultTranspose(ty, ly);
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(ty, ly);
    }
  }
}

Vector &ParOperator::GetTestLVector() const
{
  return (&trial_fespace == &test_fespace) ? trial_fespace.GetLVector2<Vector>()
                                           : test_fespace.GetLVector<Vector>();
}

ComplexParOperator::ComplexParOperator(std::unique_ptr<Operator> &&dAr,
                                       std::unique_ptr<Operator> &&dAi, const Operator *pAr,
                                       const Operator *pAi,
                                       const FiniteElementSpace &trial_fespace,
                                       const FiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexOperator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    data_A((dAr != nullptr || dAi != nullptr)
               ? ceed::CreateComplexOperator(std::move(dAr), std::move(dAi), trial_fespace)
               : std::make_unique<ComplexWrapperOperator>(pAr, pAi)),
    A(data_A.get()), trial_fespace(trial_fespace), test_fespace(test_fespace),
    use_R(test_restrict), diag_policy(Operator::DiagonalPolicy::DIAG_ONE),
    RAPr(A->Real()
             ? std::make_unique<ParOperator>(*A->Real(), trial_fespace, test_fespace, use_R)
             : nullptr),
    RAPi(A->Imag()
             ? std::make_unique<ParOperator>(*A->Imag(), trial_fespace, test_fespace, use_R)
             : nullptr)
{
  // We use the non-owning constructors for real and imaginary part ParOperators, since we
  // construct A as a ComplexWrapperOperator which has separate access to the real and
  // imaginary components.
}

ComplexParOperator::ComplexParOperator(std::unique_ptr<Operator> &&Ar,
                                       std::unique_ptr<Operator> &&Ai,
                                       const FiniteElementSpace &trial_fespace,
                                       const FiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexParOperator(std::move(Ar), std::move(Ai), nullptr, nullptr, trial_fespace,
                       test_fespace, test_restrict)
{
}

ComplexParOperator::ComplexParOperator(const Operator *Ar, const Operator *Ai,
                                       const FiniteElementSpace &trial_fespace,
                                       const FiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexParOperator(nullptr, nullptr, Ar, Ai, trial_fespace, test_fespace, test_restrict)
{
}

const hypre::MaskedProlongation *ComplexParOperator::GetMaskedProlongation() const
{
  return GetMasked(masked_P, trial_fespace, test_fespace, use_R, dbc_tdof_list);
}

bool ComplexParOperator::MappedMult(const ComplexOperator &local, const ComplexVector &x,
                                    ComplexVector &y, bool add, double a) const
{
  // See ParOperator::MappedMult.
  const auto *mP = GetMaskedProlongation();
  if (!mP || x.Real().HostRead() == y.Real().HostRead() || !ceed::SupportsMapped(local))
  {
    return false;
  }
  const int nx = mP->NumExtra();
  mapped_x.SetSize(nx);
  mapped_y.SetSize(nx);
  auto &ty = test_fespace.GetTVector<ComplexVector>();
  ComplexVector &yt = add ? ty : y;
  {
    perf::ScopedTally t(HaloTallyP());
    mP->ExtraMult(x.Real(), x.Imag(), mapped_x.Real(), mapped_x.Imag());
  }
  yt = 0.0;
  mapped_y = 0.0;
  ceed::AddMultMapped(local, {mP->LocalDofMap().data(), mP->Key()}, x.Real().HostRead(),
                      x.Imag().HostRead(), mapped_x.Real().HostRead(),
                      mapped_x.Imag().HostRead(), yt.Real().HostReadWrite(),
                      yt.Imag().HostReadWrite(), mapped_y.Real().HostReadWrite(),
                      mapped_y.Imag().HostReadWrite());
  {
    perf::ScopedTally t(HaloTallyPT());
    mP->TransposeComplete(yt.Real(), yt.Imag(), mapped_y.Real(), mapped_y.Imag());
  }
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(yt, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(yt, dbc_tdof_list, 0.0);
    }
  }
  if (add)
  {
    y.AXPY(a, yt);
  }
  return true;
}

class ComplexParOperator::SinglePrecisionOperator : public ComplexOperator
{
private:
  const ComplexParOperator &op;

public:
  SinglePrecisionOperator(const ComplexParOperator &op)
    : ComplexOperator(op.Height(), op.Width()), op(op)
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    op.MultWith(*op.A_single, x, y);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    op.AddMultWith(*op.A_single, x, y, a);
  }

  void MultTranspose(const ComplexVector &x, ComplexVector &y) const override
  {
    op.MultTranspose(x, y);
  }

  void MultHermitianTranspose(const ComplexVector &x, ComplexVector &y) const override
  {
    op.MultHermitianTranspose(x, y);
  }

  void AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                        const std::complex<double> a = 1.0) const override
  {
    op.AddMultTranspose(x, y, a);
  }

  void AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a = 1.0) const override
  {
    op.AddMultHermitianTranspose(x, y, a);
  }
};

std::size_t
ComplexParOperator::UseSinglePrecision(int *dropped,
                                       const ceed::ElementMatrixCombination *comb)
{
  MFEM_VERIFY(&trial_fespace == &test_fespace && !use_R,
              "Single-precision approximation requires a square ComplexParOperator!");
  std::size_t bytes = 0;
  if (dropped)
  {
    *dropped = 0;
  }
  if (!A_single)
  {
    A_single = ceed::CreateSinglePrecisionOperator(A->Real(), A->Imag(), trial_fespace,
                                                   &bytes, dropped, comb);
    if (A_single)
    {
      op_single = std::make_unique<SinglePrecisionOperator>(*this);
    }
  }
  return bytes;
}

void ComplexParOperator::SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                                              Operator::DiagonalPolicy policy)
{
  MFEM_VERIFY(policy == Operator::DiagonalPolicy::DIAG_ONE ||
                  policy == Operator::DiagonalPolicy::DIAG_ZERO,
              "Essential boundary condition true dof elimination for ComplexParOperator "
              "supports only DiagonalPolicy::DIAG_ONE or DiagonalPolicy::DIAG_ZERO!");
  MFEM_VERIFY(
      policy != Operator::DiagonalPolicy::DIAG_ONE || RAPr,
      "DiagonalPolicy::DIAG_ONE specified for ComplexParOperator with no real part!");
  MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                               "for rectangular ComplexParOperator!");
  tdof_list.Read();
  dbc_tdof_list.MakeRef(tdof_list);
  diag_policy = policy;
  if (RAPr)
  {
    RAPr->SetEssentialTrueDofs(tdof_list, policy);
  }
  if (RAPi)
  {
    RAPi->SetEssentialTrueDofs(tdof_list, Operator::DiagonalPolicy::DIAG_ZERO);
  }
}

Operator::DiagonalPolicy ComplexParOperator::GetDiagonalPolicy() const
{
  MFEM_VERIFY(dbc_tdof_list.Size() > 0,
              "There is no DiagonalPolicy if no essential dofs have been set!");
  return diag_policy;
}

void ComplexParOperator::AssembleDiagonal(ComplexVector &diag) const
{
  diag.UseDevice(true);
  if (assembled_diag)
  {
    diag = *assembled_diag;
    return;
  }
  diag = 0.0;
  if (RAPr)
  {
    RAPr->AssembleDiagonal(diag.Real());
  }
  if (RAPi)
  {
    RAPi->AssembleDiagonal(diag.Imag());
  }
}

void ComplexParOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  MultWith(*A, x, y);
}

void ComplexParOperator::MultWith(const ComplexOperator &local, const ComplexVector &x,
                                  ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::Mult!");

  if (MappedMult(local, x, y, false, 0.0))
  {
    return;
  }
  const auto *mP = GetMaskedProlongation();
  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  {
    perf::ScopedTally t(HaloTallyP());
    if (mP)
    {
      mP->Mult(x.Real(), x.Imag(), lx.Real(), lx.Imag());
    }
    else if (dbc_tdof_list.Size())
    {
      auto &tx = trial_fespace.GetTVector<ComplexVector>();
      tx = x;
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
      ProlongationMult(trial_fespace, tx.Real(), tx.Imag(), lx.Real(), lx.Imag());
    }
    else
    {
      ProlongationMult(trial_fespace, x.Real(), x.Imag(), lx.Real(), lx.Imag());
    }
  }

  // Apply the operator on the L-vector.
  local.Mult(lx, ly);

  {
    perf::ScopedTally t(HaloTallyPT());
    RestrictionMatrixMult(ly, y);
  }
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::MultTranspose(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::MultTranspose!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<ComplexVector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  ProlongationMultTranspose(trial_fespace, lx.Real(), lx.Imag(), y.Real(), y.Imag());
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::MultHermitianTranspose(const ComplexVector &x,
                                                ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::MultHermitianTranspose!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<ComplexVector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultHermitianTranspose(ly, lx);

  ProlongationMultTranspose(trial_fespace, lx.Real(), lx.Imag(), y.Real(), y.Imag());
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a) const
{
  AddMultWith(*A, x, y, a);
}

void ComplexParOperator::AddMultWith(const ComplexOperator &local, const ComplexVector &x,
                                     ComplexVector &y, const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::AddMult!");

  if (a.imag() == 0.0 && MappedMult(local, x, y, true, a.real()))
  {
    return;
  }
  const auto *mP = GetMaskedProlongation();
  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  {
    perf::ScopedTally t(HaloTallyP());
    if (mP)
    {
      mP->Mult(x.Real(), x.Imag(), lx.Real(), lx.Imag());
    }
    else if (dbc_tdof_list.Size())
    {
      auto &tx = trial_fespace.GetTVector<ComplexVector>();
      tx = x;
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
      ProlongationMult(trial_fespace, tx.Real(), tx.Imag(), lx.Real(), lx.Imag());
    }
    else
    {
      ProlongationMult(trial_fespace, x.Real(), x.Imag(), lx.Real(), lx.Imag());
    }
  }

  // Apply the operator on the L-vector.
  local.Mult(lx, ly);

  if (mP && a.imag() == 0.0)
  {
    perf::ScopedTally t(HaloTallyPT());
    const bool one = (diag_policy == Operator::DiagonalPolicy::DIAG_ONE);
    mP->AddMultTranspose(ly.Real(), ly.Imag(), y.Real(), y.Imag(), a.real(),
                         one ? &x.Real() : nullptr, one ? &x.Imag() : nullptr);
    return;
  }
  auto &ty = test_fespace.GetTVector<ComplexVector>();
  {
    perf::ScopedTally t(HaloTallyPT());
    RestrictionMatrixMult(ly, ty);
  }
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    }
  }
  y.AXPY(a, ty);
}

void ComplexParOperator::AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                          const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultTranspose!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<ComplexVector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  auto &tx = trial_fespace.GetTVector<ComplexVector>();
  ProlongationMultTranspose(trial_fespace, lx.Real(), lx.Imag(), tx.Real(), tx.Imag());
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
    }
  }
  y.AXPY(a, tx);
}

void ComplexParOperator::AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                                   const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultHermitianTranspose!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &ty = test_fespace.GetTVector<ComplexVector>();
    ty = x;
    linalg::SetSubVector(ty, dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultHermitianTranspose(ly, lx);

  auto &tx = trial_fespace.GetTVector<ComplexVector>();
  ProlongationMultTranspose(trial_fespace, lx.Real(), lx.Imag(), tx.Real(), tx.Imag());
  if (dbc_tdof_list.Size())
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
    }
  }
  y.AXPY(a, tx);
}

void ComplexParOperator::RestrictionMatrixMult(const ComplexVector &ly,
                                               ComplexVector &ty) const
{
  if (!use_R)
  {
    ProlongationMultTranspose(test_fespace, ly.Real(), ly.Imag(), ty.Real(), ty.Imag());
  }
  else
  {
    if (const auto *R = test_fespace.GetRestrictionRuns())
    {
      R->Mult(ly.Real(), ty.Real());
      R->Mult(ly.Imag(), ty.Imag());
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->Mult(ly.Real(), ty.Real());
      test_fespace.GetRestrictionMatrix()->Mult(ly.Imag(), ty.Imag());
    }
  }
}

void ComplexParOperator::RestrictionMatrixMultTranspose(const ComplexVector &ty,
                                                        ComplexVector &ly) const
{
  if (!use_R)
  {
    ProlongationMult(test_fespace, ty.Real(), ty.Imag(), ly.Real(), ly.Imag());
  }
  else
  {
    if (const auto *R = test_fespace.GetRestrictionRuns())
    {
      R->MultTranspose(ty.Real(), ly.Real());
      R->MultTranspose(ty.Imag(), ly.Imag());
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(ty.Real(), ly.Real());
      test_fespace.GetRestrictionMatrix()->MultTranspose(ty.Imag(), ly.Imag());
    }
  }
}

ComplexVector &ComplexParOperator::GetTestLVector() const
{
  return (&trial_fespace == &test_fespace) ? trial_fespace.GetLVector2<ComplexVector>()
                                           : test_fespace.GetLVector<ComplexVector>();
}

// Helper that checks if two containers (Vector or Array<T>) are actually references to the
// same underlying data.
template <typename C>
bool ReferencesSameMemory(const C &c1, const C &c2)
{
  const auto &m1 = c1.GetMemory();
  const auto &m2 = c2.GetMemory();
  return (m1.HostIsValid() && m2.HostIsValid() && c1.HostRead() == c2.HostRead()) ||
         (m1.DeviceIsValid() && m2.DeviceIsValid() && c1.Read() == c2.Read());
}

// Apply the essential dofs shared by the operators to their sum, with the highest diagonal
// policy among them in enum order.
template <typename Ops, typename SumOperType>
void SetSumEssentialTrueDofs(const Ops &ops, SumOperType &O)
{
  // Extract essential dof pointer from first operator with one.
  auto it_ess = std::ranges::find_if(ops, [](auto p)
                                     { return p != nullptr && p->GetEssentialTrueDofs(); });
  if (it_ess == ops.end())
  {
    return;
  }
  const auto *ess_dofs = (*it_ess)->GetEssentialTrueDofs();

  // Check other existent essential dof arrays are references.
  MFEM_VERIFY(std::ranges::all_of(ops,
                                  [&](auto p)
                                  {
                                    if (p == nullptr)
                                    {
                                      return true;
                                    }
                                    auto p_ess_dofs = p->GetEssentialTrueDofs();
                                    return p_ess_dofs == nullptr ||
                                           ReferencesSameMemory(*ess_dofs, *p_ess_dofs);
                                  }),
              "If essential dofs are set, all suboperators must agree on them!");

  // Use implied ordering of enumeration.
  Operator::DiagonalPolicy policy = Operator::DiagonalPolicy::DIAG_ZERO;
  for (auto p : ops)
  {
    policy = (p && p->GetEssentialTrueDofs()) ? std::max(policy, p->GetDiagonalPolicy())
                                              : policy;
  }
  O.SetEssentialTrueDofs(*ess_dofs, policy);
}

// Combine a collection of ParOperator into a weighted summation. If set_essential is true,
// extract the essential dofs from the operator array, and apply to the summed operator.
template <std::size_t N>
std::unique_ptr<ParOperator>
BuildParSumOperator(const std::array<double, N> &coeff,
                    const std::array<const ParOperator *, N> &ops, bool set_essential)
{
  auto it = std::ranges::find_if(ops, [](auto p) { return p != nullptr; });
  MFEM_VERIFY(it != ops.end(),
              "BuildParSumOperator requires at least one valid ParOperator!");
  const auto first_op = *it;
  const auto &fespace = first_op->TrialFiniteElementSpace();
  MFEM_VERIFY(std::ranges::all_of(
                  ops, [&fespace](auto p)
                  { return p == nullptr || &p->TrialFiniteElementSpace() == &fespace; }),
              "All ComplexParOperators must have the same FiniteElementSpace!");

  auto sum = std::make_unique<SumOperator>(first_op->LocalOperator().Height(),
                                           first_op->LocalOperator().Width());
  for (std::size_t i = 0; i < coeff.size(); i++)
  {
    if (ops[i] && coeff[i] != 0)
    {
      sum->AddOperator(ops[i]->LocalOperator(), coeff[i]);
    }
  }

  auto O = std::make_unique<ParOperator>(std::move(sum), fespace);
  if (set_essential)
  {
    SetSumEssentialTrueDofs(ops, *O);
  }
  return O;
}

// Combine a collection of ComplexParOperator into a weighted summation. If set_essential is
// true, extract the essential dofs from the operator array, and apply to the summed
// operator.
template <std::size_t N>
std::unique_ptr<ComplexParOperator>
BuildParSumOperator(const std::array<std::complex<double>, N> &coeff,
                    const std::array<const ComplexParOperator *, N> &ops,
                    bool set_essential)
{
  return BuildParSumOperator(
      std::vector<std::complex<double>>(coeff.begin(), coeff.end()),
      std::vector<const ComplexParOperator *>(ops.begin(), ops.end()), set_essential);
}

// Runtime-sized overload for a caller-determined number of terms.
std::unique_ptr<ComplexParOperator>
BuildParSumOperator(const std::vector<std::complex<double>> &coeff,
                    const std::vector<const ComplexParOperator *> &ops, bool set_essential)
{
  MFEM_VERIFY(coeff.size() == ops.size(),
              "BuildParSumOperator requires matching coefficient and operator counts!");
  auto it = std::ranges::find_if(ops, [](auto p) { return p != nullptr; });
  MFEM_VERIFY(it != ops.end(),
              "BuildParSumOperator requires at least one valid ComplexParOperator!");
  const auto first_op = *it;
  const auto &fespace = first_op->TrialFiniteElementSpace();
  MFEM_VERIFY(std::ranges::all_of(
                  ops, [&fespace](auto p)
                  { return p == nullptr || &p->TrialFiniteElementSpace() == &fespace; }),
              "All ComplexParOperators must have the same FiniteElementSpace!");

  auto sumr = std::make_unique<SumOperator>(first_op->LocalOperator().Height(),
                                            first_op->LocalOperator().Width());
  auto sumi = std::make_unique<SumOperator>(first_op->LocalOperator().Height(),
                                            first_op->LocalOperator().Width());
  for (std::size_t i = 0; i < coeff.size(); i++)
  {
    if (ops[i] && coeff[i].real() != 0)
    {
      if (ops[i]->LocalOperator().Real())
      {
        sumr->AddOperator(*ops[i]->LocalOperator().Real(), coeff[i].real());
      }
      if (ops[i]->LocalOperator().Imag())
      {
        sumi->AddOperator(*ops[i]->LocalOperator().Imag(), coeff[i].real());
      }
    }
    if (ops[i] && coeff[i].imag() != 0)
    {
      if (ops[i]->LocalOperator().Imag())
      {
        sumr->AddOperator(*ops[i]->LocalOperator().Imag(), -coeff[i].imag());
      }
      if (ops[i]->LocalOperator().Real())
      {
        sumi->AddOperator(*ops[i]->LocalOperator().Real(), coeff[i].imag());
      }
    }
  }
  auto O = std::make_unique<ComplexParOperator>(std::move(sumr), std::move(sumi), fespace);
  if (set_essential)
  {
    SetSumEssentialTrueDofs(ops, *O);
  }
  return O;
}

std::unique_ptr<ComplexParOperator>
BuildParSumOperator(const std::vector<std::complex<double>> &coeff,
                    const std::vector<const ComplexOperator *> &ops, bool set_essential)
{
  std::vector<const ComplexParOperator *> par_ops(ops.size());
  std::ranges::transform(ops, par_ops.begin(), [](const ComplexOperator *op)
                         { return dynamic_cast<const ComplexParOperator *>(op); });
  return BuildParSumOperator(coeff, par_ops, set_essential);
}

std::unique_ptr<ComplexOperator>
BuildOperatorWithA2(const std::vector<std::complex<double>> &coeff,
                    const std::vector<const ComplexOperator *> &ops,
                    const ComplexOperator *A2)
{
  if (A2 && !dynamic_cast<const ComplexParOperator *>(A2))
  {
    auto sum = BuildParSumOperator(coeff, ops, true);
    return std::make_unique<SumComplexOperator>(std::move(sum), *A2);
  }
  std::vector<std::complex<double>> c(coeff);
  std::vector<const ComplexOperator *> o(ops);
  c.emplace_back(1.0, 0.0);
  o.push_back(A2);
  return BuildParSumOperator(c, o, true);
}

template <std::size_t N>
std::unique_ptr<ComplexParOperator>
BuildParSumOperator(std::complex<double> (&&coeff_in)[N],
                    const ComplexParOperator *(&&ops_in)[N], bool set_essential)
{
  return BuildParSumOperator(std::to_array<std::complex<double>>(std::move(coeff_in)),
                             std::to_array<const ComplexParOperator *>(std::move(ops_in)),
                             set_essential);
}

template <std::size_t N, typename ScalarType, typename OperType>
std::unique_ptr<std::conditional_t<std::is_base_of_v<ComplexOperator, OperType>,
                                   ComplexParOperator, ParOperator>>
BuildParSumOperator(ScalarType (&&coeff_in)[N], const OperType *(&&ops_in)[N],
                    bool set_essential)
{
  using ParOperType =
      typename std::conditional_t<std::is_base_of_v<ComplexOperator, OperType>,
                                  ComplexParOperator, ParOperator>;

  std::array<const ParOperType *, N> par_ops;
  std::transform(ops_in, ops_in + N, par_ops.begin(),
                 [](const OperType *op) { return dynamic_cast<const ParOperType *>(op); });

  return BuildParSumOperator(std::to_array<ScalarType>(std::move(coeff_in)),
                             std::move(par_ops), set_essential);
}

// Explicit instantiation.
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[1],
                                                          const Operator *(&&)[1], bool);
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[2],
                                                          const Operator *(&&)[2], bool);
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[3],
                                                          const Operator *(&&)[3], bool);
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[4],
                                                          const Operator *(&&)[4], bool);

template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[1], const ComplexOperator *(&&)[1], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[2], const ComplexOperator *(&&)[2], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[3], const ComplexOperator *(&&)[3], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[4], const ComplexOperator *(&&)[4], bool);

namespace linalg
{

std::unique_ptr<mfem::HypreParMatrix>
BlockDiagonalParMatrix(MPI_Comm comm, HYPRE_BigInt global_rows, HYPRE_BigInt global_cols,
                       HYPRE_BigInt *row_starts, HYPRE_BigInt *col_starts,
                       std::unique_ptr<mfem::SparseMatrix> &&diag)
{
  diag->Finalize();
  auto *offd = new mfem::SparseMatrix(diag->Height(), 0, 0);
  offd->Finalize();
  return std::make_unique<mfem::HypreParMatrix>(comm, global_rows, global_cols, row_starts,
                                                col_starts, diag.release(), offd, nullptr,
                                                true);
}

std::unique_ptr<mfem::HypreParMatrix> RAP(const mfem::HypreParMatrix &A,
                                          const mfem::HypreParMatrix &P)
{
  hypre_ParCSRMatrix *AP = hypre_ParCSRMatMat(const_cast<mfem::HypreParMatrix &>(A),
                                              const_cast<mfem::HypreParMatrix &>(P));
  auto PtAP = std::make_unique<mfem::HypreParMatrix>(
      hypre_ParCSRTMatMatKT(const_cast<mfem::HypreParMatrix &>(P), AP, 1), true);
  hypre_ParCSRMatrixDestroy(AP);
  hypre_ParCSRMatrixSetNumNonzeros(*PtAP);
  return PtAP;
}

}  // namespace linalg

}  // namespace palace
