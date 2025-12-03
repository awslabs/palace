// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "rap.hpp"

#include "fem/bilinearform.hpp"
#include "linalg/hypre.hpp"

namespace palace
{

ParOperator::ParOperator(std::unique_ptr<Operator> &&dA, const Operator *pA,
                         const FiniteElementSpace &trial_fespace,
                         const FiniteElementSpace &test_fespace, bool test_restrict)
  : Operator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    data_A(std::move(dA)), A((data_A != nullptr) ? data_A.get() : pA),
    trial_fespace(trial_fespace), test_fespace(test_fespace), use_R(test_restrict),
    diag_policy(DiagonalPolicy::DIAG_ONE), RAP(nullptr)
{
  MFEM_VERIFY(A, "Cannot construct ParOperator from an empty matrix!");
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

void ParOperator::SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                                       DiagonalPolicy policy)
{
  MFEM_VERIFY(policy == DiagonalPolicy::DIAG_ONE || policy == DiagonalPolicy::DIAG_ZERO,
              "Essential boundary condition true dof elimination for ParOperator supports "
              "only DiagonalPolicy::DIAG_ONE or DiagonalPolicy::DIAG_ZERO!");
  MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                               "for rectangular ParOperator!");
  tdof_list.Read();
  dbc_tdof_list.MakeRef(tdof_list);
  diag_policy = policy;
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
      test_fespace.Get().GetDofOffsets(), trial_fespace.Get().GetDofOffsets(), 0, sA->NNZ(),
      0);
  hypre_CSRMatrix *hA_diag = hypre_ParCSRMatrixDiag(hA);
  hypre_ParCSRMatrixDiag(hA) = *const_cast<hypre::HypreCSRMatrix *>(sA);
  hypre_ParCSRMatrixInitialize(hA);

  const mfem::HypreParMatrix *P = trial_fespace.Get().Dof_TrueDof_Matrix();
  if (!use_R)
  {
    const mfem::HypreParMatrix *Rt = test_fespace.Get().Dof_TrueDof_Matrix();
    RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAPKT(*Rt, hA, *P, 1),
                                                 true);
  }
  else
  {
    mfem::HypreParMatrix *hR = new mfem::HypreParMatrix(
        test_fespace.GetComm(), test_fespace.GlobalTrueVSize(), test_fespace.GlobalVSize(),
        test_fespace.Get().GetTrueDofOffsets(), test_fespace.Get().GetDofOffsets(),
        const_cast<mfem::SparseMatrix *>(test_fespace.GetRestrictionMatrix()));
    hypre_ParCSRMatrix *AP = hypre_ParCSRMatMat(hA, *P);
    RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatMat(*hR, AP), true);
    hypre_ParCSRMatrixDestroy(AP);
    delete hR;
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

  // For an AMR mesh, a convergent diagonal is assembled with |P|ᵀ dₗ, where |P| has
  // entry-wise absolute values of the conforming prolongation operator.
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "Diagonal assembly is only available for square ParOperator!");
  auto &lx = trial_fespace.GetLVector<Vector>();
  A->AssembleDiagonal(lx);

  // Parallel assemble and eliminate essential true dofs.
  const Operator *P = test_fespace.GetProlongationMatrix();
  if (const auto *hP = dynamic_cast<const mfem::HypreParMatrix *>(P))
  {
    hP->AbsMultTranspose(1.0, lx, 0.0, diag);
  }
  else
  {
    P->MultTranspose(lx, diag);
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
}

void ParOperator::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::Mult!");
  if (RAP)
  {
    RAP->Mult(x, y);
    return;
  }

  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
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

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  RestrictionMatrixMult(ly, y);
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
}

void ParOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::AddMult!");
  if (RAP)
  {
    RAP->AddMult(x, y, a);
    return;
  }

  auto &lx = trial_fespace.GetLVector<Vector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
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

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  auto &ty = test_fespace.GetTVector<Vector>();
  RestrictionMatrixMult(ly, ty);
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
}

void ParOperator::RestrictionMatrixMult(const Vector &ly, Vector &ty) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->MultTranspose(ly, ty);
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->Mult(ly, ty);
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
    test_fespace.GetRestrictionMatrix()->MultTranspose(ty, ly);
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
               ? std::make_unique<ComplexWrapperOperator>(std::move(dAr), std::move(dAi))
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
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::Mult!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &tx = trial_fespace.GetTVector<ComplexVector>();
    tx = x;
    linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(tx.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(tx.Imag(), lx.Imag());
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(x.Imag(), lx.Imag());
  }

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  RestrictionMatrixMult(ly, y);
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

  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), y.Real());
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), y.Imag());
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

  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), y.Real());
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), y.Imag());
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
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::AddMult!");

  auto &lx = trial_fespace.GetLVector<ComplexVector>();
  auto &ly = GetTestLVector();
  if (dbc_tdof_list.Size())
  {
    auto &tx = trial_fespace.GetTVector<ComplexVector>();
    tx = x;
    linalg::SetSubVector(tx, dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(tx.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(tx.Imag(), lx.Imag());
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(x.Imag(), lx.Imag());
  }

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  auto &ty = test_fespace.GetTVector<ComplexVector>();
  RestrictionMatrixMult(ly, ty);
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
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), tx.Real());
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), tx.Imag());
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
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), tx.Real());
  trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), tx.Imag());
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
    test_fespace.GetProlongationMatrix()->MultTranspose(ly.Real(), ty.Real());
    test_fespace.GetProlongationMatrix()->MultTranspose(ly.Imag(), ty.Imag());
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->Mult(ly.Real(), ty.Real());
    test_fespace.GetRestrictionMatrix()->Mult(ly.Imag(), ty.Imag());
  }
}

void ComplexParOperator::RestrictionMatrixMultTranspose(const ComplexVector &ty,
                                                        ComplexVector &ly) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->Mult(ty.Real(), ly.Real());
    test_fespace.GetProlongationMatrix()->Mult(ty.Imag(), ly.Imag());
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->MultTranspose(ty.Real(), ly.Real());
    test_fespace.GetRestrictionMatrix()->MultTranspose(ty.Imag(), ly.Imag());
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

// Combine a collection of ParOperator into a weighted summation. If set_essential is true,
// extract the essential dofs from the operator array, and apply to the summed operator.
template <std::size_t N>
std::unique_ptr<ParOperator>
BuildParSumOperator(const std::array<double, N> &coeff,
                    const std::array<const ParOperator *, N> &ops, bool set_essential)
{
  auto it = std::find_if(ops.begin(), ops.end(), [](auto p) { return p != nullptr; });
  MFEM_VERIFY(it != ops.end(),
              "BuildParSumOperator requires at least one valid ParOperator!");
  const auto first_op = *it;
  const auto &fespace = first_op->TrialFiniteElementSpace();
  MFEM_VERIFY(
      std::all_of(ops.begin(), ops.end(), [&fespace](auto p)
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
    // Extract essential dof pointer from first operator with one.
    auto it_ess = std::find_if(ops.begin(), ops.end(), [](auto p)
                               { return p != nullptr && p->GetEssentialTrueDofs(); });
    if (it_ess == ops.end())
    {
      return O;
    }
    const auto *ess_dofs = (*it_ess)->GetEssentialTrueDofs();

    // Check other existent essential dof arrays are references.
    MFEM_VERIFY(std::all_of(ops.begin(), ops.end(),
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
    O->SetEssentialTrueDofs(*ess_dofs, policy);
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
  auto it = std::find_if(ops.begin(), ops.end(), [](auto p) { return p != nullptr; });
  MFEM_VERIFY(it != ops.end(),
              "BuildParSumOperator requires at least one valid ComplexParOperator!");
  const auto first_op = *it;
  const auto &fespace = first_op->TrialFiniteElementSpace();
  MFEM_VERIFY(
      std::all_of(ops.begin(), ops.end(), [&fespace](auto p)
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
    // Extract essential dof pointer from first operator with one.
    auto it_ess = std::find_if(ops.begin(), ops.end(), [](auto p)
                               { return p != nullptr && p->GetEssentialTrueDofs(); });
    if (it_ess == ops.end())
    {
      return O;
    }
    const auto *ess_dofs = (*it_ess)->GetEssentialTrueDofs();

    // Check other existent essential dof arrays are references.
    MFEM_VERIFY(std::all_of(ops.begin(), ops.end(),
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
    O->SetEssentialTrueDofs(*ess_dofs, policy);
  }
  return O;
}

// TODO: replace with std::to_array in c++20.
namespace detail
{
// Helper for conversion to std::array.
template <class T, std::size_t N, std::size_t... I>
constexpr std::array<std::remove_cv_t<T>, N> to_array_impl(T (&&a)[N],
                                                           std::index_sequence<I...>)
{
  return {{std::move(a[I])...}};
}
}  // namespace detail

template <class T, std::size_t N>
constexpr std::array<std::remove_cv_t<T>, N> to_array(T (&&a)[N])
{
  return detail::to_array_impl(std::move(a), std::make_index_sequence<N>{});
}

template <std::size_t N>
std::unique_ptr<ComplexParOperator>
BuildParSumOperator(std::complex<double> (&&coeff_in)[N],
                    const ComplexParOperator *(&&ops_in)[N], bool set_essential)
{
  return BuildParSumOperator(to_array<std::complex<double>>(std::move(coeff_in)),
                             to_array<const ComplexParOperator *>(std::move(ops_in)),
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

  return BuildParSumOperator(to_array<ScalarType>(std::move(coeff_in)), std::move(par_ops),
                             set_essential);
}

// Explicit instantiation.
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[2],
                                                          const Operator *(&&)[2], bool);
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[3],
                                                          const Operator *(&&)[3], bool);
template std::unique_ptr<ParOperator> BuildParSumOperator(double (&&)[4],
                                                          const Operator *(&&)[4], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[2], const ComplexOperator *(&&)[2], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[3], const ComplexOperator *(&&)[3], bool);
template std::unique_ptr<ComplexParOperator>
    BuildParSumOperator(std::complex<double> (&&)[4], const ComplexOperator *(&&)[4], bool);

}  // namespace palace
