// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "rap.hpp"

namespace palace
{

ParOperator::ParOperator(std::unique_ptr<Operator> &&dA, Operator *pA,
                         const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         bool test_restrict)
  : Operator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    data_A(std::move(dA)), A((data_A != nullptr) ? data_A.get() : pA),
    trial_fespace(trial_fespace), test_fespace(test_fespace), use_R(test_restrict),
    dbc_tdof_list(nullptr), diag_policy(DiagonalPolicy::DIAG_ONE), RAP(nullptr)
{
  MFEM_VERIFY(A, "Cannot construct ParOperator from an empty matrix!");
  lx.SetSize(A->Width());
  ly.SetSize(A->Height());
  ty.SetSize(width);
}

ParOperator::ParOperator(std::unique_ptr<Operator> &&A,
                         const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         bool test_restrict)
  : ParOperator(std::move(A), nullptr, trial_fespace, test_fespace, test_restrict)
{
}

ParOperator::ParOperator(Operator &A, const mfem::ParFiniteElementSpace &trial_fespace,
                         const mfem::ParFiniteElementSpace &test_fespace,
                         bool test_restrict)
  : ParOperator(nullptr, &A, trial_fespace, test_fespace, test_restrict)
{
}

const Operator &ParOperator::LocalOperator() const
{
  MFEM_ASSERT(A, "No local matrix available for ParOperator::LocalOperator!");
  return *A;
}

Operator &ParOperator::LocalOperator()
{
  MFEM_ASSERT(A, "No local matrix available for ParOperator::LocalOperator!");
  return *A;
}

void ParOperator::SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                                       DiagonalPolicy policy)
{
  MFEM_VERIFY(policy == DiagonalPolicy::DIAG_ONE || policy == DiagonalPolicy::DIAG_ZERO,
              "Essential boundary condition true dof elimination for ParOperator supports "
              "only DiagonalPolicy::DIAG_ONE or DiagonalPolicy::DIAG_ZERO!");
  MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                               "for rectangular ParOperator!");
  dbc_tdof_list = &tdof_list;
  diag_policy = policy;
}

void ParOperator::AssembleDiagonal(Vector &diag) const
{
  // For an AMR mesh, a convergent diagonal is assembled with |P|ᵀ dₗ, where |P| has
  // entry-wise absolute values of the conforming prolongation operator.
  MFEM_VERIFY(&trial_fespace == &test_fespace,
              "Diagonal assembly is only available for square ParOperator!");
  if (const auto *bfA = dynamic_cast<const mfem::BilinearForm *>(A))
  {
    // if (bfA->HasSpMat())   //XX TODO WIP FOR LIBCEED
    // {
    //   bfA->SpMat().GetDiag(ly);
    // }
    // else if (bfA->HasExt())
    // {
    //   bfA->Ext().AssembleDiagonal(ly);
    // }
    // else
    {
      MFEM_ABORT("Unable to assemble the local operator diagonal of BilinearForm!");
    }
  }
  else if (const auto *sA = dynamic_cast<const mfem::SparseMatrix *>(A))
  {
    sA->GetDiag(ly);
  }
  else
  {
    MFEM_ABORT("ParOperator::AssembleDiagonal requires A as a BilinearForm or "
               "SparseMatrix!");
  }

  // Parallel assemble and eliminate essential true dofs.
  const Operator *P = test_fespace.GetProlongationMatrix();
  if (const auto *hP = dynamic_cast<const mfem::HypreParMatrix *>(P))
  {
    hP->AbsMultTranspose(1.0, ly, 0.0, diag);
  }
  else
  {
    P->MultTranspose(ly, diag);
  }
  if (dbc_tdof_list)
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(diag, *dbc_tdof_list, 1.0);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(diag, *dbc_tdof_list, 0.0);
    }
  }
}

mfem::HypreParMatrix &ParOperator::ParallelAssemble() const
{
  if (RAP)
  {
    return *RAP;
  }

  // Build the square or rectangular assembled HypreParMatrix.
  mfem::SparseMatrix *sA = dynamic_cast<mfem::SparseMatrix *>(A);
  MFEM_VERIFY(sA, "ParOperator::ParallelAssemble requires A as a SparseMatrix!");
  if (&trial_fespace == &test_fespace)
  {
    mfem::HypreParMatrix *hA =
        new mfem::HypreParMatrix(trial_fespace.GetComm(), trial_fespace.GlobalVSize(),
                                 trial_fespace.GetDofOffsets(), sA);
    const mfem::HypreParMatrix *P = trial_fespace.Dof_TrueDof_Matrix();
    RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*P, *hA, *P), true);
    delete hA;
  }
  else
  {
    mfem::HypreParMatrix *hA = new mfem::HypreParMatrix(
        trial_fespace.GetComm(), test_fespace.GlobalVSize(), trial_fespace.GlobalVSize(),
        test_fespace.GetDofOffsets(), trial_fespace.GetDofOffsets(), sA);
    const mfem::HypreParMatrix *P = trial_fespace.Dof_TrueDof_Matrix();
    if (!use_R)
    {
      const mfem::HypreParMatrix *Rt = test_fespace.Dof_TrueDof_Matrix();
      RAP =
          std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*Rt, *hA, *P), true);
    }
    else
    {

      // XX TODO WIP MIGHT NEED TO EDIT WITHOUT MFEM PATCH...

      mfem::SparseMatrix *sRt = mfem::Transpose(*test_fespace.GetRestrictionMatrix());
      mfem::HypreParMatrix *hRt = new mfem::HypreParMatrix(
          test_fespace.GetComm(), test_fespace.GlobalVSize(),
          test_fespace.GlobalTrueVSize(), test_fespace.GetDofOffsets(),
          test_fespace.GetTrueDofOffsets(), sRt);
      RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*hRt, *hA, *P),
                                                   true);
      delete sRt;
      delete hRt;
    }
    delete hA;
  }
  hypre_ParCSRMatrixSetNumNonzeros(*RAP);

  // Eliminate boundary conditions on the assembled (square) matrix.
  if (dbc_tdof_list)
  {
    RAP->EliminateBC(*dbc_tdof_list, diag_policy);
  }
  return *RAP;
}

void ParOperator::EliminateRHS(const Vector &x, Vector &b) const
{
  if (!dbc_tdof_list)
  {
    return;
  }

  MFEM_VERIFY(A, "No local matrix available for ParOperator::EliminateRHS!");
  ty = 0.0;
  linalg::SetSubVector(ty, *dbc_tdof_list, x);
  trial_fespace.GetProlongationMatrix()->Mult(ty, lx);

  // Apply the unconstrained operator.
  A->Mult(lx, ly);

  RestrictionMatrixAddMult(ly, b, -1.0);
  if (diag_policy == DiagonalPolicy::DIAG_ONE)
  {
    linalg::SetSubVector(b, *dbc_tdof_list, x);
  }
  else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
  {
    linalg::SetSubVector(b, *dbc_tdof_list, 0.0);
  }
}

void ParOperator::Mult(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::Mult!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(ty, lx);
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x, lx);
  }

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  RestrictionMatrixMult(ly, y);
  if (dbc_tdof_list)
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, 0.0);
    }
  }
}

void ParOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::AddMult!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(ty, lx);
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x, lx);
  }

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  if (dbc_tdof_list)
  {
    RestrictionMatrixMult(ly, ty);
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    }
    y.Add(a, ty);
  }
  else
  {
    RestrictionMatrixAddMult(ly, y, a);
  }
}

void ParOperator::MultTranspose(const Vector &x, Vector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ParOperator::MultTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  trial_fespace.GetProlongationMatrix()->MultTranspose(lx, y);
  if (dbc_tdof_list)
  {
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, 0.0);
    }
  }
}

void ParOperator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ParOperator::AddMultTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx, ty);
    if (diag_policy == DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, x);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    }
    y.Add(a, ty);
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx, y, a);
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
    test_fespace.GetRestrictionMatrix()->Mult(ly, ty);
  }
}

void ParOperator::RestrictionMatrixAddMult(const Vector &ly, Vector &ty,
                                           const double a) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->AddMultTranspose(ly, ty, a);
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->AddMult(ly, ty, a);
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

ComplexParOperator::ComplexParOperator(std::unique_ptr<Operator> &&dAr,
                                       std::unique_ptr<Operator> &&dAi, Operator *pAr,
                                       Operator *pAi,
                                       const mfem::ParFiniteElementSpace &trial_fespace,
                                       const mfem::ParFiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexOperator(test_fespace.GetTrueVSize(), trial_fespace.GetTrueVSize()),
    data_A((dAr != nullptr || dAi != nullptr)
               ? std::make_unique<ComplexWrapperOperator>(std::move(dAr), std::move(dAi))
               : std::make_unique<ComplexWrapperOperator>(pAr, pAi)),
    A(data_A.get()), trial_fespace(trial_fespace), test_fespace(test_fespace),
    use_R(test_restrict), dbc_tdof_list(nullptr),
    diag_policy(Operator::DiagonalPolicy::DIAG_ONE),
    RAPr(A->HasReal()
             ? std::make_unique<ParOperator>(*A->Real(), trial_fespace, test_fespace, use_R)
             : nullptr),
    RAPi(A->HasImag()
             ? std::make_unique<ParOperator>(*A->Imag(), trial_fespace, test_fespace, use_R)
             : nullptr)
{
  // We use the non-owning constructors for real and imaginary part ParOperators. We know A
  // is a ComplexWrapperOperator which has separate access to the real and imaginary
  // components.
  lx.SetSize(A->Width());
  ly.SetSize(A->Height());
  ty.SetSize(width);
}

ComplexParOperator::ComplexParOperator(std::unique_ptr<Operator> &&Ar,
                                       std::unique_ptr<Operator> &&Ai,
                                       const mfem::ParFiniteElementSpace &trial_fespace,
                                       const mfem::ParFiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexParOperator(std::move(Ar), std::move(Ai), nullptr, nullptr, trial_fespace,
                       test_fespace, test_restrict)
{
}

ComplexParOperator::ComplexParOperator(Operator *Ar, Operator *Ai,
                                       const mfem::ParFiniteElementSpace &trial_fespace,
                                       const mfem::ParFiniteElementSpace &test_fespace,
                                       bool test_restrict)
  : ComplexParOperator(nullptr, nullptr, Ar, Ai, trial_fespace, test_fespace, test_restrict)
{
}

const ComplexOperator &ComplexParOperator::LocalOperator() const
{
  MFEM_ASSERT(A, "No local matrix available for ComplexParOperator::LocalOperator!");
  return *A;
}

ComplexOperator &ComplexParOperator::LocalOperator()
{
  MFEM_ASSERT(A, "No local matrix available for ComplexParOperator::LocalOperator!");
  return *A;
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
  dbc_tdof_list = &tdof_list;
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

void ComplexParOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::Mult!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(ty.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(ty.Imag(), lx.Imag());
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(x.Imag(), lx.Imag());
  }

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  RestrictionMatrixMult(ly, y);
  if (dbc_tdof_list)
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ComplexParOperator::AddMult!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    trial_fespace.GetProlongationMatrix()->Mult(ty.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(ty.Imag(), lx.Imag());
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->Mult(x.Real(), lx.Real());
    trial_fespace.GetProlongationMatrix()->Mult(x.Imag(), lx.Imag());
  }

  // Apply the operator on the L-vector.
  ly = 0.0;
  A->AddMult(lx, ly, a);

  if (dbc_tdof_list)
  {
    RestrictionMatrixMult(ly, ty);
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    }
    y += ty;
  }
  else
  {
    RestrictionMatrixAddMult(ly, y, 1.0);
  }
}

void ComplexParOperator::MultTranspose(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::MultTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
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
  if (dbc_tdof_list)
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                          const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  lx = 0.0;
  A->AddMultTranspose(ly, lx, a);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), ty.Real());
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), ty.Imag());
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    }
    y += ty;
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Real(), y.Real());
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Imag(), y.Imag());
  }
}

void ComplexParOperator::MultHermitianTranspose(const ComplexVector &x,
                                                ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::MultHermitianTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
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
  if (dbc_tdof_list)
  {
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, 0.0);
    }
  }
}

void ComplexParOperator::AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                                   const std::complex<double> a) const
{
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultHermitianTranspose!");
  if (dbc_tdof_list)
  {
    ty = x;
    linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    RestrictionMatrixMultTranspose(ty, ly);
  }
  else
  {
    RestrictionMatrixMultTranspose(x, ly);
  }

  // Apply the operator on the L-vector.
  lx = 0.0;
  A->AddMultHermitianTranspose(ly, lx, a);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), ty.Real());
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), ty.Imag());
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, x);
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO)
    {
      linalg::SetSubVector(ty, *dbc_tdof_list, 0.0);
    }
    y += ty;
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Real(), y.Real());
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Imag(), y.Imag());
  }
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

void ComplexParOperator::RestrictionMatrixAddMult(const ComplexVector &ly,
                                                  ComplexVector &ty, const double a) const
{
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->AddMultTranspose(ly.Real(), ty.Real(), a);
    test_fespace.GetProlongationMatrix()->AddMultTranspose(ly.Imag(), ty.Imag(), a);
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->AddMult(ly.Real(), ty.Real(), a);
    test_fespace.GetRestrictionMatrix()->AddMult(ly.Imag(), ty.Imag(), a);
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

}  // namespace palace
