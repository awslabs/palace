// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "rap.hpp"

#include <general/forall.hpp>

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
  if (auto *bfA = dynamic_cast<const mfem::BilinearForm *>(A))
  {
    if (bfA->HasSpMat())
    {
      bfA->SpMat().GetDiag(ly);
    }
    else if (bfA->HasExt())
    {
      bfA->Ext().AssembleDiagonal(ly);
    }
    else
    {
      MFEM_ABORT("Unable to assemble the local operator diagonal of BilinearForm!");
    }
  }
  else if (auto *sA = dynamic_cast<const mfem::SparseMatrix *>(A))
  {
    sA->GetDiag(ly);
  }
  else
  {
    MFEM_ABORT("ParOperator::AssembleDiagonal requires A as a BilinearForm or "
               "SparseMatrix!");
  }

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
      diag.SetSubVector(*dbc_tdof_list, 1.0);
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO)
    {
      diag.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
  }
}

mfem::HypreParMatrix &ParOperator::ParallelAssemble()
{
  if (RAP)
  {
    return *RAP;
  }

  // XX TODO: For mfem::AssemblyLevel::PARTIAL, we cannot use CeedOperatorFullAssemble for
  //          a ND space with p > 1. We should throw an error here that the user needs to
  //          use AssemblyLevel::LEGACY in this case.

  // Build the square or rectangular RAP HypreParMatrix.
  if (&trial_fespace == &test_fespace)
  {
    mfem::SparseMatrix *lA;
    bool own_lA = false;
    if (auto *bfA = dynamic_cast<mfem::BilinearForm *>(A))
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
    else if (auto *sA = dynamic_cast<mfem::SparseMatrix *>(A))
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
        new mfem::HypreParMatrix(trial_fespace.GetComm(), trial_fespace.GlobalVSize(),
                                 trial_fespace.GetDofOffsets(), lA);
    const mfem::HypreParMatrix *P = trial_fespace.Dof_TrueDof_Matrix();
    RAP = std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*P, *hA, *P), true);
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
    if (auto *mbfA = dynamic_cast<mfem::MixedBilinearForm *>(A))
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
    else if (auto *sA = dynamic_cast<mfem::SparseMatrix *>(A))
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
        trial_fespace.GetComm(), test_fespace.GlobalVSize(), trial_fespace.GlobalVSize(),
        test_fespace.GetDofOffsets(), trial_fespace.GetDofOffsets(), lA);
    const mfem::HypreParMatrix *P = trial_fespace.Dof_TrueDof_Matrix();
    if (!use_R)
    {
      const mfem::HypreParMatrix *Rt = test_fespace.Dof_TrueDof_Matrix();
      RAP =
          std::make_unique<mfem::HypreParMatrix>(hypre_ParCSRMatrixRAP(*Rt, *hA, *P), true);
    }
    else
    {
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
    if (own_lA)
    {
      delete lA;
    }
  }
  hypre_ParCSRMatrixSetNumNonzeros(*RAP);

  // Eliminate boundary conditions on the assembled (square) matrix.
  if (dbc_tdof_list)
  {
    MFEM_VERIFY(
        &trial_fespace == &test_fespace,
        "Only square ParOperator should have same trial and test eliminated tdofs!");
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
  {
    const int N = dbc_tdof_list->Size();
    const auto *idx = dbc_tdof_list->Read();
    const auto *X = x.Read();
    auto *TY = ty.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const int id = idx[i];
                   TY[id] = X[id];
                 });
  }

  // Apply the unconstrained operator.
  trial_fespace.GetProlongationMatrix()->Mult(ty, lx);
  A->Mult(lx, ly);

  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->AddMultTranspose(ly, b, -1.0);
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->AddMult(ly, b, -1.0);
  }
  if (diag_policy == DiagonalPolicy::DIAG_ONE && height == width)
  {
    const int N = dbc_tdof_list->Size();
    const auto *idx = dbc_tdof_list->Read();
    const auto *X = x.Read();
    auto *B = b.ReadWrite();
    mfem::forall(N,
                 [=] MFEM_HOST_DEVICE(int i)
                 {
                   const int id = idx[i];
                   B[id] = X[id];
                 });
  }
  else if (diag_policy == DiagonalPolicy::DIAG_ZERO || height != width)
  {
    b.SetSubVector(*dbc_tdof_list, 0.0);
  }
  else
  {
    MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
  }
}

void ParOperator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Incompatible dimensions for ParOperator::AddMult!");
  if (dbc_tdof_list)
  {
    ty = x;
    ty.SetSubVector(*dbc_tdof_list, 0.0);
  }
  trial_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty : x, lx);

  // Apply the operator on the L-vector.
  A->Mult(lx, ly);

  if (dbc_tdof_list)
  {
    if (!use_R)
    {
      test_fespace.GetProlongationMatrix()->MultTranspose(ly, ty);
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->Mult(ly, ty);
    }
    if (diag_policy == DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = dbc_tdof_list->Size();
      const auto *idx = dbc_tdof_list->Read();
      const auto *X = x.Read();
      auto *TY = ty.ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TY[id] = X[id];
                   });
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    y.Add(a, ty);
  }
  else
  {
    if (!use_R)
    {
      test_fespace.GetProlongationMatrix()->AddMultTranspose(ly, y, a);
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->AddMult(ly, y, a);
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
    ty.SetSubVector(*dbc_tdof_list, 0.0);
  }
  if (!use_R)
  {
    test_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty : x, ly);
  }
  else
  {
    test_fespace.GetRestrictionMatrix()->MultTranspose(dbc_tdof_list ? ty : x, ly);
  }

  // Apply the operator on the L-vector.
  A->MultTranspose(ly, lx);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx, ty);
    if (diag_policy == DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = dbc_tdof_list->Size();
      const auto *idx = dbc_tdof_list->Read();
      const auto *X = x.Read();
      auto *TY = ty.ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TY[id] = X[id];
                   });
    }
    else if (diag_policy == DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    y.Add(a, ty);
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx, y, a);
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

void ComplexParOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  MFEM_ASSERT(xr.Size() == width && xi.Size() == width && yr.Size() == height &&
                  yi.Size() == height,
              "Incompatible dimensions for ComplexParOperator::AddMult!");
  if (dbc_tdof_list)
  {
    ty.Real() = xr;
    ty.Imag() = xi;
    ty.SetSubVector(*dbc_tdof_list, 0.0);
  }
  if (!zero_real)
  {
    trial_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Real() : xr, lx.Real());
  }
  if (!zero_imag)
  {
    trial_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Imag() : xi, lx.Imag());
  }

  // Apply the operator on the L-vector.
  ly = 0.0;
  A->AddMult(lx, ly, a);

  if (dbc_tdof_list)
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
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = dbc_tdof_list->Size();
      const auto *idx = dbc_tdof_list->Read();
      const auto *XR = xr.Read();
      const auto *XI = xi.Read();
      auto *TYR = ty.Real().ReadWrite();
      auto *TYI = ty.Imag().ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TYR[id] = XR[id];
                     TYI[id] = XI[id];
                   });
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    yr += ty.Real();
    yi += ty.Imag();
  }
  else
  {
    if (!use_R)
    {
      test_fespace.GetProlongationMatrix()->MultTranspose(ly.Real(), yr);
      test_fespace.GetProlongationMatrix()->MultTranspose(ly.Imag(), yi);
    }
    else
    {
      test_fespace.GetRestrictionMatrix()->Mult(ly.Real(), yr);
      test_fespace.GetRestrictionMatrix()->Mult(ly.Imag(), yi);
    }
  }
}

void ComplexParOperator::AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                                          const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  MFEM_ASSERT(xr.Size() == height && xi.Size() == height && yr.Size() == width &&
                  yi.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultTranspose!");
  if (dbc_tdof_list)
  {
    ty.Real() = xr;
    ty.Imag() = xi;
    ty.SetSubVector(*dbc_tdof_list, 0.0);
  }
  if (!use_R)
  {
    if (!zero_real)
    {
      test_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Real() : xr, ly.Real());
    }
    if (!zero_imag)
    {
      test_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Imag() : xi, ly.Imag());
    }
  }
  else
  {
    if (!zero_real)
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(dbc_tdof_list ? ty.Real() : xr,
                                                         ly.Real());
    }
    if (!zero_imag)
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(dbc_tdof_list ? ty.Imag() : xi,
                                                         ly.Imag());
    }
  }

  // Apply the operator on the L-vector.
  lx = 0.0;
  A->AddMultTranspose(ly, lx, a);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), ty.Real());
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), ty.Imag());
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = dbc_tdof_list->Size();
      const auto *idx = dbc_tdof_list->Read();
      const auto *XR = xr.Read();
      const auto *XI = xi.Read();
      auto *TYR = ty.Real().ReadWrite();
      auto *TYI = ty.Imag().ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TYR[id] = XR[id];
                     TYI[id] = XI[id];
                   });
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    yr += ty.Real();
    yi += ty.Imag();
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Real(), yr);
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Imag(), yi);
  }
}

void ComplexParOperator::AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                                   const std::complex<double> a) const
{
  constexpr bool zero_real = false;
  constexpr bool zero_imag = false;
  const Vector &xr = x.Real();
  const Vector &xi = x.Imag();
  Vector &yr = y.Real();
  Vector &yi = y.Imag();
  MFEM_ASSERT(xr.Size() == height && xi.Size() == height && yr.Size() == width &&
                  yi.Size() == width,
              "Incompatible dimensions for ComplexParOperator::AddMultHermitianTranspose!");
  if (dbc_tdof_list)
  {
    ty.Real() = xr;
    ty.Imag() = xi;
    ty.SetSubVector(*dbc_tdof_list, 0.0);
  }
  if (!use_R)
  {
    if (!zero_real)
    {
      test_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Real() : xr, ly.Real());
    }
    if (!zero_imag)
    {
      test_fespace.GetProlongationMatrix()->Mult(dbc_tdof_list ? ty.Imag() : xi, ly.Imag());
    }
  }
  else
  {
    if (!zero_real)
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(dbc_tdof_list ? ty.Real() : xr,
                                                         ly.Real());
    }
    if (!zero_imag)
    {
      test_fespace.GetRestrictionMatrix()->MultTranspose(dbc_tdof_list ? ty.Imag() : xi,
                                                         ly.Imag());
    }
  }

  // Apply the operator on the L-vector.
  lx = 0.0;
  A->AddMultHermitianTranspose(ly, lx, a);

  if (dbc_tdof_list)
  {
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Real(), ty.Real());
    trial_fespace.GetProlongationMatrix()->MultTranspose(lx.Imag(), ty.Imag());
    if (diag_policy == Operator::DiagonalPolicy::DIAG_ONE && height == width)
    {
      const int N = dbc_tdof_list->Size();
      const auto *idx = dbc_tdof_list->Read();
      const auto *XR = xr.Read();
      const auto *XI = xi.Read();
      auto *TYR = ty.Real().ReadWrite();
      auto *TYI = ty.Imag().ReadWrite();
      mfem::forall(N,
                   [=] MFEM_HOST_DEVICE(int i)
                   {
                     const int id = idx[i];
                     TYR[id] = XR[id];
                     TYI[id] = XI[id];
                   });
    }
    else if (diag_policy == Operator::DiagonalPolicy::DIAG_ZERO || height != width)
    {
      ty.SetSubVector(*dbc_tdof_list, 0.0);
    }
    else
    {
      MFEM_ABORT("Unsupported Operator::DiagonalPolicy for ParOperator!");
    }
    yr += ty.Real();
    yi += ty.Imag();
  }
  else
  {
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Real(), yr);
    trial_fespace.GetProlongationMatrix()->AddMultTranspose(lx.Imag(), yi);
  }
}

}  // namespace palace
