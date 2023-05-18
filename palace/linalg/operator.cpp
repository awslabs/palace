// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <general/forall.hpp>
#include "linalg/complex.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

ParOperator::ParOperator(std::unique_ptr<Operator> &&A,
                         const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         bool test_restrict)
  : Operator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()), A_(std::move(A)),
    trial_fespace_(trial_fespace), test_fespace_(test_fespace), use_R_(test_restrict),
    trial_dbc_tdof_list_(nullptr), test_dbc_tdof_list_(nullptr),
    diag_policy_(DiagonalPolicy::DIAG_ONE), RAP_(nullptr), save_A_(false)
{
  MFEM_VERIFY(A_, "Cannot construct ParOperator from an empty matrix!");
  lx_.SetSize(A_->Width());
  ly_.SetSize(A_->Height());
  tx_.SetSize(width);
  if (height != width)
  {
    ty_.SetSize(height);
  }
  else
  {
    ty_.MakeRef(tx_, 0, height);
  }
}

void ParOperator::EliminateRHS(const Vector &x, Vector &b) const
{
  if (!trial_dbc_tdof_list_ || !test_dbc_tdof_list_)
  {
    return;
  }

  MFEM_VERIFY(A_, "No local matrix available for ParOperator::EliminateRHS!");
  tx_ = 0.0;
  {
    const int N = trial_dbc_tdof_list_->Size();
    const auto *idx = trial_dbc_tdof_list_->Read();
    const auto *X = x.Read();
    auto *TX = tx_.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const int id = idx[i];
                   TX[id] = X[id];
                 });
  }

  // Apply the unconstrained operator.
  trial_fespace_.GetProlongationMatrix()->Mult(tx_, lx_);
  A_->Mult(lx_, ly_);
  if (!use_R_)
  {
    test_fespace_.GetProlongationMatrix()->AddMultTranspose(ly_, b, -1.0);
  }
  else
  {
    test_fespace_.GetRestrictionMatrix()->AddMult(ly_, b, -1.0);
  }

  {
    if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = test_dbc_tdof_list_->Size();
      const auto *idx = test_dbc_tdof_list_->Read();
      const auto *X = x.Read();
      auto *B = b.ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     B[id] = X[id];
                   });
    }
    else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
    {
      b.SetSubVector(*test_dbc_tdof_list_, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
  }
}

void ParOperator::AssembleDiagonal(Vector &diag) const
{
  if (RAP_)
  {
    RAP_->GetDiag(diag);
    return;
  }

  // For an AMR mesh, a convergent diagonal is assembled with |P|ᵀ dₗ, where |P| has
  // entry-wise absolute values of the conforming prolongation operator.
  MFEM_VERIFY(&trial_fespace_ == &test_fespace_,
              "Diagonal assembly is only available for square ParOperator!");
  if (auto *bfA = dynamic_cast<mfem::BilinearForm *>(A_.get()))
  {
    if (bfA->HasSpMat())
    {
      bfA->SpMat().GetDiag(ly_);
    }
    else if (bfA->HasExt())
    {
      bfA->Ext().AssembleDiagonal(ly_);
    }
    else
    {
      MFEM_ABORT("Unable to assemble the local operator diagonal of BilinearForm!");
    }
  }
  else if (auto *sA = dynamic_cast<mfem::SparseMatrix *>(A_.get()))
  {
    sA->GetDiag(ly_);
  }
  else
  {
    MFEM_ABORT("ParOperator::AssembleDiagonal requires A as a BilinearForm or "
               "SparseMatrix!");
  }

  const Operator *P = test_fespace_.GetProlongationMatrix();
  if (const auto *hP = dynamic_cast<const mfem::HypreParMatrix *>(P))
  {
    hP->AbsMultTranspose(1.0, ly_, 0.0, diag);
  }
  else
  {
    P->MultTranspose(ly_, diag);
  }

  if (test_dbc_tdof_list_)
  {
    if (diag_policy_ == DiagonalPolicy::DIAG_ONE)
    {
      diag.SetSubVector(*test_dbc_tdof_list_, 1.0);
    }
    else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO)
    {
      diag.SetSubVector(*test_dbc_tdof_list_, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
  }
}

mfem::HypreParMatrix &ParOperator::ParallelAssemble()
{
  if (RAP_)
  {
    return *RAP_;
  }

  // XX TODO: For mfem::AssemblyLevel::PARTIAL, we cannot use CeedOperatorFullAssemble for
  //          a ND space with p > 1. We should throw an error here that the user needs to
  //          use AssemblyLevel::LEGACY in this case.

  // Build the square or rectangular RᵀAP HypreParMatrix.
  if (&trial_fespace_ == &test_fespace_)
  {
    mfem::SparseMatrix *lA;
    bool own_lA = false;
    if (auto *bfA = dynamic_cast<mfem::BilinearForm *>(A_.get()))
    {
#ifdef MFEM_USE_CEED
      if (bfA->HasSpMat())
      {
        lA = &bfA->SpMat();
      }
      else if (bfA->HasExt())
      {
        lA = mfem::ceed::CeedOperatorFullAssemble(*bfA);
        own_lA = true;
      }
      else
      {
        MFEM_ABORT("Unable to assemble the local operator for parallel assembly of "
                   "BilinearForm!");
      }
#else
      MFEM_VERIFY(bfA->HasSpMat(),
                  "Missing assembled SparseMatrix for parallel assembly of BilinearForm!");
      lA = &bfA->SpMat();
#endif
    }
    else if (auto *sA = dynamic_cast<mfem::SparseMatrix *>(A_.get()))
    {
      lA = sA;
    }
    else
    {
      MFEM_ABORT("ParOperator::ParallelAssemble requires A as a BilinearForm or "
                 "SparseMatrix!");
      lA = nullptr;
    }
    mfem::HypreParMatrix *hA =
        new mfem::HypreParMatrix(trial_fespace_.GetComm(), trial_fespace_.GlobalVSize(),
                                 trial_fespace_.GetDofOffsets(), lA);
    const mfem::HypreParMatrix *P = trial_fespace_.Dof_TrueDof_Matrix();
    RAP_ = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*P, *hA, *P), true);
    delete hA;
    if (own_lA)
    {
      delete lA;
    }
  }
  else
  {
    mfem::SparseMatrix *lA;
    bool own_lA = false;
    if (auto *mbfA = dynamic_cast<mfem::MixedBilinearForm *>(A_.get()))
    {
#ifdef MFEM_USE_CEED
      if (mbfA->HasSpMat())
      {
        lA = &mbfA->SpMat();
      }
      else if (bfA->HasExt())
      {
        lA = mfem::ceed::CeedOperatorFullAssemble(*bfA);
        own_lA = true;
      }
      else
      {
        MFEM_ABORT("Unable to assemble the local operator for parallel assembly of "
                   "MixedBilinearForm!");
      }
#else
      MFEM_VERIFY(
          mbfA->HasSpMat(),
          "Missing assembled SparseMatrix for parallel assembly of MixedBilinearForm!");
      lA = &mbfA->SpMat();
#endif
    }
    else if (auto *sA = dynamic_cast<mfem::SparseMatrix *>(A_.get()))
    {
      lA = sA;
    }
    else
    {
      MFEM_ABORT("ParOperator::ParallelAssemble requires A as a MixedBilinearForm or "
                 "SparseMatrix!");
      lA = nullptr;
    }
    mfem::HypreParMatrix *hA = new mfem::HypreParMatrix(
        trial_fespace_.GetComm(), test_fespace_.GlobalVSize(), trial_fespace_.GlobalVSize(),
        test_fespace_.GetDofOffsets(), trial_fespace_.GetDofOffsets(), lA);
    const mfem::HypreParMatrix *P = trial_fespace_.Dof_TrueDof_Matrix();
    if (!use_R_)
    {
      const mfem::HypreParMatrix *Rt = test_fespace_.Dof_TrueDof_Matrix();
      RAP_ =
          std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*Rt, *hA, *P), true);
    }
    else
    {
      mfem::SparseMatrix *sRt = mfem::Transpose(*test_fespace_.GetRestrictionMatrix());
      mfem::HypreParMatrix *hRt = new mfem::HypreParMatrix(
          test_fespace_.GetComm(), test_fespace_.GlobalVSize(),
          test_fespace_.GlobalTrueVSize(), test_fespace_.GetDofOffsets(),
          test_fespace_.GetTrueDofOffsets(), sRt);
      RAP_ = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*hRt, *hA, *P),
                                                    true);
      delete sRt;
      delete hRt;
    }
    delete hA;
    if (own_lA)
    {
      delete lA;
    }
  }
  hypre_ParCSRMatrixSetNumNonzeros(*RAP_);

  // Delete the original local operator.
  if (!save_A_)
  {
    A_.reset();
  }

  // Eliminate boundary conditions on the assembled matrix.
  if (test_dbc_tdof_list_ || trial_dbc_tdof_list_)
  {
    if (test_dbc_tdof_list_ == trial_dbc_tdof_list_)
    {
      // Elimination for a square operator.
      MFEM_VERIFY(
          &trial_fespace_ == &test_fespace_,
          "Only square ParOperator should have same trial and test eliminated tdofs!");
      RAP_->EliminateBC(*trial_dbc_tdof_list_, diag_policy_);
    }
    else
    {
      // Rectangular elimination sets all eliminated rows/columns to zero.
      if (test_dbc_tdof_list_)
      {
        RAP_->EliminateRows(*test_dbc_tdof_list_);
      }
      if (trial_dbc_tdof_list_)
      {
        mfem::HypreParMatrix *RAPe = RAP_->EliminateCols(*trial_dbc_tdof_list_);
        delete RAPe;
      }
    }
  }
  return *RAP_;
}

void ParOperator::Mult(const Vector &x, Vector &y) const
{
  if (RAP_)
  {
    RAP_->Mult(x, y);
    return;
  }
  y = 0.0;
  AddMult(x, y);
}

void ParOperator::MultTranspose(const Vector &x, Vector &y) const
{
  if (RAP_)
  {
    RAP_->MultTranspose(x, y);
    return;
  }
  y = 0.0;
  AddMultTranspose(x, y);
}

void ParOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  if (RAP_)
  {
    RAP_->AddMult(x, y, a);
    return;
  }
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::AddMult!");
  if (trial_dbc_tdof_list_)
  {
    tx_ = x;
    tx_.SetSubVector(*trial_dbc_tdof_list_, 0.0);
  }
  trial_fespace_.GetProlongationMatrix()->Mult(trial_dbc_tdof_list_ ? tx_ : x, lx_);

  // Apply the operator on the L-vector.
  A_->Mult(lx_, ly_);

  if (test_dbc_tdof_list_)
  {
    if (!use_R_)
    {
      test_fespace_.GetProlongationMatrix()->MultTranspose(ly_, ty_);
    }
    else
    {
      test_fespace_.GetRestrictionMatrix()->Mult(ly_, ty_);
    }
    if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = test_dbc_tdof_list_->Size();
      const auto *idx = test_dbc_tdof_list_->Read();
      const auto *X = x.Read();
      auto *TY = ty_.ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TY[id] = X[id];
                   });
    }
    else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty_.SetSubVector(*test_dbc_tdof_list_, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    y.Add(a, ty_);
  }
  else
  {
    if (!use_R_)
    {
      test_fespace_.GetProlongationMatrix()->AddMultTranspose(ly_, y, a);
    }
    else
    {
      test_fespace_.GetRestrictionMatrix()->AddMult(ly_, y, a);
    }
  }
}

void ParOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  if (RAP_)
  {
    RAP_->AddMultTranspose(x, y, a);
    return;
  }
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ParOperator::AddMultTranspose!");
  if (test_dbc_tdof_list_)
  {
    ty_ = x;
    ty_.SetSubVector(*test_dbc_tdof_list_, 0.0);
  }
  if (!use_R_)
  {
    test_fespace_.GetProlongationMatrix()->Mult(test_dbc_tdof_list_ ? ty_ : x, ly_);
  }
  else
  {
    test_fespace_.GetRestrictionMatrix()->MultTranspose(test_dbc_tdof_list_ ? ty_ : x, ly_);
  }

  // Apply the operator on the L-vector.
  A_->MultTranspose(ly_, lx_);

  if (trial_dbc_tdof_list_)
  {
    trial_fespace_.GetProlongationMatrix()->MultTranspose(lx_, tx_);
    if (diag_policy_ == DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = trial_dbc_tdof_list_->Size();
      const auto *idx = trial_dbc_tdof_list_->Read();
      const auto *X = x.Read();
      auto *TX = tx_.ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TX[id] = X[id];
                   });
    }
    else if (diag_policy_ == DiagonalPolicy::DIAG_ZERO || height != width)
    {
      tx_.SetSubVector(*test_dbc_tdof_list_, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    y.Add(a, tx_);
  }
  else
  {
    trial_fespace_.GetProlongationMatrix()->AddMultTranspose(lx_, y, a);
  }
}

void DiagonalOperator::Mult(const Vector &x, Vector &y) const
{
  const int N = height;
  const auto *D = d_.Read();
  const auto *X = x.Read();
  auto *Y = y.Write();
  mfem::forall(N, [=] MFEM_HOST_DEVICE(int i) { Y[i] = D[i] * X[i]; });
}

namespace linalg
{

double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym, double tol, int max_it)
{
  // The SumOperator does not take ownership of A and allows the ComplexWrapperOperator
  // to own its input.
  ComplexWrapperOperator Ar(std::make_unique<SumOperator>(A, 1.0), nullptr);
  return SpectralNorm(comm, Ar, sym, tol, max_it);
}

double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm, double tol,
                    int max_it)
{
  // XX TODO: Use ARPACK or SLEPc for this when configured.
#if defined(PALACE_WITH_SLEPC)
  return slepc::GetMaxSingularValue(comm, A, herm, tol, max_it);
#else
  // Power iteration loop: ||A||₂² = λₙ(Aᴴ A).
  int it = 0;
  double res = 0.0;
  double l, l0 = 0.0;
  ComplexVector u(A.Height()), v(A.Height());
  SetRandom(comm, u);
  Normalize(comm, u);
  while (it < max_it)
  {
    A.Mult(u, v);
    if (herm)
    {
      u = v;
    }
    else
    {
      A.MultHermitianTranspose(v, u);
    }
    l = Normalize(comm, u);
    if (it > 0)
    {
      res = std::abs(l - l0) / l0;
      if (res < tol)
      {
        break;
      }
    }
    l0 = l;
    it++;
  }
  if (it >= max_it)
  {
    Mpi::Warning(comm,
                 "Power iteration did not converge in {:d} iterations, res = {:.3e}, "
                 "lambda = {:.3e}!\n",
                 it, res, l);
  }
  return herm ? l : std::sqrt(l);
#endif
}

}  // namespace linalg

}  // namespace palace
