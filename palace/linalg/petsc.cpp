// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "petsc.hpp"

#include <petsc.h>
#include <petscblaslapack.h>
#include <general/forall.hpp>
// #include "linalg/hypre.hpp"
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"

static PetscErrorCode __mat_shell_init(Mat);
static PetscErrorCode __mat_shell_destroy(Mat);
static PetscErrorCode __mat_shell_apply(Mat, Vec, Vec);
static PetscErrorCode __mat_shell_apply_transpose(Mat, Vec, Vec);
static PetscErrorCode __mat_shell_apply_hermitian_transpose(Mat, Vec, Vec);
static PetscErrorCode __mat_shell_apply_add(Mat, Vec, Vec);
static PetscErrorCode __mat_shell_apply_transpose_add(Mat, Vec, Vec);
static PetscErrorCode __mat_shell_apply_hermitian_transpose_add(Mat, Vec, Vec);
#if defined(PETSC_USE_COMPLEX)
static PetscErrorCode __mat_shell_apply(Mat, const mfem::Vector &, Vec);
static PetscErrorCode __mat_shell_apply_transpose(Mat, const mfem::Vector &, Vec);
static PetscErrorCode __mat_shell_apply_hermitian_transpose(Mat, const mfem::Vector &, Vec);
#endif
static PetscErrorCode __mat_shell_get_diagonal(Mat, Vec);
// static PetscErrorCode __mat_shell_shift(Mat, PetscScalar);
// static PetscErrorCode __mat_shell_scale(Mat, PetscScalar);
// static PetscErrorCode __mat_shell_conj(Mat);
// static PetscErrorCode __mat_shell_axpy(Mat, PetscScalar, Mat, MatStructure);
// static PetscErrorCode __mat_shell_norm(Mat, NormType, PetscReal *);
static PetscErrorCode __mat_shell_real_part(Mat);
static PetscErrorCode __mat_shell_imag_part(Mat);
static PetscErrorCode __mat_convert_hypreParCSR_AIJ(hypre_ParCSRMatrix *, Mat *);
static PetscErrorCode __array_container_destroy(void *);

namespace palace::petsc
{

using mfem::ForallWrap;

void Initialize(int &argc, char **&argv, const char rc_file[], const char help[])
{
  PalacePetscCall(PetscInitialize(&argc, &argv, rc_file, help));
}

void Finalize()
{
  PalacePetscCall(PetscFinalize());
}

// PetscScatter methods
PetscScatter::PetscScatter(PetscScatter::Type type, const PetscParVector &x,
                           std::unique_ptr<PetscParVector> &y)
{
  Vec yy;
  if (type == Type::TO_ZERO)
  {
    PalacePetscCall(VecScatterCreateToZero(x, &ctx, &yy));
  }
  else  // type == Type::TO_ALL
  {
    PalacePetscCall(VecScatterCreateToAll(x, &ctx, &yy));
  }
  y = std::make_unique<PetscParVector>(yy, false);
}

PetscScatter::~PetscScatter()
{
  PalacePetscCall(VecScatterDestroy(&ctx));
}

void PetscScatter::Forward(const PetscParVector &x, PetscParVector &y)
{
  PalacePetscCall(VecScatterBegin(ctx, x, y, INSERT_VALUES, SCATTER_FORWARD));
  PalacePetscCall(VecScatterEnd(ctx, x, y, INSERT_VALUES, SCATTER_FORWARD));
}

void PetscScatter::Reverse(const PetscParVector &x, PetscParVector &y)
{
  PalacePetscCall(VecScatterBegin(ctx, x, y, INSERT_VALUES, SCATTER_REVERSE));
  PalacePetscCall(VecScatterEnd(ctx, x, y, INSERT_VALUES, SCATTER_REVERSE));
}

// PetscParVector methods

PetscParVector::PetscParVector(const PetscParMatrix &A, bool transpose)
{
  if (!transpose)
  {
    PalacePetscCall(MatCreateVecs(A, &x, nullptr));
  }
  else
  {
    PalacePetscCall(MatCreateVecs(A, nullptr, &x));
  }
}

PetscParVector::PetscParVector(MPI_Comm comm, const mfem::Vector &y)
{
  PalacePetscCall(VecCreate(comm, &x));
  PalacePetscCall(VecSetSizes(x, y.Size(), PETSC_DECIDE));
  PalacePetscCall(VecSetType(x, VECSTANDARD));
  SetFromVector(y);
}

PetscParVector::PetscParVector(const mfem::Vector &y)
{
  PalacePetscCall(VecCreateSeq(PETSC_COMM_SELF, y.Size(), &x));
  SetFromVector(y);
}

#if defined(PETSC_USE_COMPLEX)
PetscParVector::PetscParVector(MPI_Comm comm, const mfem::Vector &yr,
                               const mfem::Vector &yi)
{
  MFEM_VERIFY(yr.Size() == yi.Size(),
              "Mismatch in size of real and imaginary vector parts!");
  PalacePetscCall(VecCreate(comm, &x));
  PalacePetscCall(VecSetSizes(x, yr.Size(), PETSC_DECIDE));
  PalacePetscCall(VecSetType(x, VECSTANDARD));
  SetFromVectors(yr, yi);
}

PetscParVector::PetscParVector(const mfem::Vector &yr, const mfem::Vector &yi)
{
  MFEM_VERIFY(yr.Size() == yi.Size(),
              "Mismatch in size of real and imaginary vector parts!");
  PalacePetscCall(VecCreateSeq(PETSC_COMM_SELF, yr.Size(), &x));
  SetFromVectors(yr, yi);
}
#endif

PetscParVector::PetscParVector(MPI_Comm comm, PetscInt n, PetscInt N)
{
  PalacePetscCall(VecCreateMPI(comm, n, N, &x));
}

// PetscParVector::PetscParVector(PetscInt n)
// {
//   PalacePetscCall(VecCreateSeq(PETSC_COMM_SELF, n, &x));
// }

PetscParVector::PetscParVector(MPI_Comm comm, PetscInt n, PetscInt N, PetscScalar *data)
{
  PalacePetscCall(VecCreateMPIWithArray(comm, 1, n, N, data, &x));
}

PetscParVector::PetscParVector(PetscInt n, PetscScalar *data)
{
  PalacePetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, data, &x));
}

PetscParVector::PetscParVector(const PetscParVector &y)
{
  PalacePetscCall(VecDuplicate(y, &x));
  Copy(y);
}

PetscParVector::PetscParVector(Vec y, bool ref)
{
  x = y;
  if (ref)
  {
    PalacePetscCall(PetscObjectReference(reinterpret_cast<PetscObject>(y)));
  }
}

PetscParVector::~PetscParVector()
{
  PalacePetscCall(VecDestroy(&x));
}

void PetscParVector::Copy(const PetscParVector &y)
{
  MFEM_VERIFY(GetSize() == y.GetSize(), "Invalid size!");
  PalacePetscCall(VecCopy(y, x));
}

void PetscParVector::GetToVector(mfem::Vector &v, PetscInt start, PetscInt end) const
{
  const PetscScalar *xv;
  if (start < 0)
  {
    start = 0;
  }
  if (end < 0)
  {
    end = GetSize();
  }
  MFEM_VERIFY(0 <= start && start <= end && end <= GetSize() && v.Size() == end - start,
              "Invalid start/end indices for vector extraction!");
  PalacePetscCall(VecGetArrayRead(x, &xv));
  auto vv = v.Write();
#if defined(PETSC_USE_COMPLEX)
  MFEM_FORALL(i, end - start, { vv[i] = PetscRealPart(xv[i + start]); });
#else
  MFEM_FORALL(i, end - start, { vv[i] = xv[i + start]; });
#endif
  PalacePetscCall(VecRestoreArrayRead(x, &xv));
}

void PetscParVector::SetFromVector(const mfem::Vector &v)
{
  PetscScalar *xv;
  MFEM_VERIFY(GetSize() == v.Size(), "Invalid size!");
  PalacePetscCall(VecGetArray(x, &xv));
  const auto vv = v.Read();
  MFEM_FORALL(i, GetSize(), { xv[i] = vv[i]; });
  PalacePetscCall(VecRestoreArray(x, &xv));
}

void PetscParVector::AddFromVector(const mfem::Vector &v)
{
  PetscScalar *xv;
  MFEM_VERIFY(GetSize() == v.Size(), "Invalid size!");
  PalacePetscCall(VecGetArray(x, &xv));
  const auto vv = v.Read();
  MFEM_FORALL(i, GetSize(), { xv[i] += vv[i]; });
  PalacePetscCall(VecRestoreArray(x, &xv));
}

#if defined(PETSC_USE_COMPLEX)
void PetscParVector::GetToVectors(mfem::Vector &vr, mfem::Vector &vi, PetscInt start,
                                  PetscInt end) const
{
  const PetscScalar *xv;
  if (start < 0)
  {
    start = 0;
  }
  if (end < 0)
  {
    end = GetSize();
  }
  MFEM_VERIFY(0 <= start && start <= end && end <= GetSize() && vr.Size() == end - start &&
                  vi.Size() == end - start,
              "Invalid start/end indices for vector extraction!");
  PalacePetscCall(VecGetArrayRead(x, &xv));
  auto vvr = vr.Write();
  auto vvi = vi.Write();
  MFEM_FORALL(i, end - start, {
    vvr[i] = PetscRealPart(xv[i + start]);
    vvi[i] = PetscImaginaryPart(xv[i + start]);
  });
  PalacePetscCall(VecRestoreArrayRead(x, &xv));
}

ComplexVector PetscParVector::GetToVectors(PetscInt start, PetscInt end) const
{
  ComplexVector cv;
  GetToVectors(cv.real, cv.imag, start, end);
  return cv;
}

void PetscParVector::SetFromVectors(const mfem::Vector &vr, const mfem::Vector &vi)
{
  PetscScalar *xv;
  MFEM_VERIFY(GetSize() == vr.Size() || GetSize() == vi.Size(), "Invalid size!");
  PalacePetscCall(VecGetArray(x, &xv));
  const auto vvr = vr.Read();
  const auto vvi = vi.Read();
  MFEM_FORALL(i, GetSize(), {
    // xv[i] = vvr[i] + PETSC_i * vvi[i];
    reinterpret_cast<PetscReal *>(&xv[i])[0] = vvr[i];
    reinterpret_cast<PetscReal *>(&xv[i])[1] = vvi[i];
  });
  PalacePetscCall(VecRestoreArray(x, &xv));
}

void PetscParVector::AddFromVectors(const mfem::Vector &vr, const mfem::Vector &vi)
{
  PetscScalar *xv;
  MFEM_VERIFY(GetSize() == vr.Size() || GetSize() == vi.Size(), "Invalid size!");
  PalacePetscCall(VecGetArray(x, &xv));
  const auto vvr = vr.Read();
  const auto vvi = vi.Read();
  MFEM_FORALL(i, GetSize(), {
    // xv[i] += vvr[i] + PETSC_i * vvi[i];
    reinterpret_cast<PetscReal *>(&xv[i])[0] += vvr[i];
    reinterpret_cast<PetscReal *>(&xv[i])[1] += vvi[i];
  });
  PalacePetscCall(VecRestoreArray(x, &xv));
}
#endif

PetscScalar *PetscParVector::GetArray()
{
  PetscScalar *data;
  PalacePetscCall(VecGetArray(x, &data));
  return data;
}

const PetscScalar *PetscParVector::GetArrayRead() const
{
  const PetscScalar *data;
  PalacePetscCall(VecGetArrayRead(x, &data));
  return data;
}

void PetscParVector::RestoreArray(PetscScalar *data)
{
  PalacePetscCall(VecRestoreArray(x, &data));
}

void PetscParVector::RestoreArrayRead(const PetscScalar *data) const
{
  PalacePetscCall(VecRestoreArrayRead(x, &data));
}

void PetscParVector::PlaceArray(const PetscScalar *data)
{
  PalacePetscCall(VecPlaceArray(x, data));
}

void PetscParVector::ResetArray()
{
  PalacePetscCall(VecResetArray(x));
}

PetscInt PetscParVector::GetSize() const
{
  PetscInt n;
  PalacePetscCall(VecGetLocalSize(x, &n));
  return n;
}

PetscInt PetscParVector::GetGlobalSize() const
{
  PetscInt N;
  PalacePetscCall(VecGetSize(x, &N));
  return N;
}

void PetscParVector::Resize(PetscInt n, bool copy)
{
  Vec y;
  const PetscScalar *xv;
  PetscScalar *yv;
  PetscInt n0 = GetSize();
  VecType type;
  if (n0 == n)
  {
    return;
  }
  PalacePetscCall(VecGetType(x, &type));
  PalacePetscCall(VecCreate(GetComm(), &y));
  PalacePetscCall(VecSetSizes(y, n, PETSC_DECIDE));
  PalacePetscCall(VecSetType(y, type));
  if (copy)
  {
    PalacePetscCall(VecGetArrayRead(x, &xv));
    PalacePetscCall(VecGetArray(y, &yv));
    MFEM_FORALL(i, std::min(n, n0), { yv[i] = xv[i]; });
    PalacePetscCall(VecRestoreArrayRead(x, &xv));
    PalacePetscCall(VecRestoreArray(y, &yv));
  }
  PalacePetscCall(VecDestroy(&x));
  x = y;
}

void PetscParVector::SetZero()
{
  PalacePetscCall(VecZeroEntries(x));
}

void PetscParVector::SetRandom()
{
  PetscRandom rand;
  MPI_Comm comm = GetComm();
  PalacePetscCall(PetscRandomCreate(comm, &rand));
#if defined(PETSC_USE_COMPLEX)
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0 - PETSC_i, 1.0 + PETSC_i));
#else
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0, 1.0));
#endif
  PalacePetscCall(VecSetRandom(x, rand));
  PalacePetscCall(PetscRandomDestroy(&rand));
}

#if defined(PETSC_USE_COMPLEX)
void PetscParVector::SetRandomReal()
{
  PetscRandom rand;
  MPI_Comm comm = GetComm();
  PalacePetscCall(PetscRandomCreate(comm, &rand));
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0, 1.0));
  PalacePetscCall(VecSetRandom(x, rand));
  PalacePetscCall(PetscRandomDestroy(&rand));
}
#endif

void PetscParVector::SetRandomSign(bool init)
{
  PetscScalar *xv;
  if (!init)
  {
    SetRandomReal();
  }
  PalacePetscCall(VecGetArray(x, &xv));
  MFEM_FORALL(i, GetSize(), {
    // Leave zeros alone.
    xv[i] =
        (PetscRealPart(xv[i]) > 0.0) ? 1.0 : ((PetscRealPart(xv[i]) < 0.0) ? -1.0 : 0.0);
  });
  PalacePetscCall(VecRestoreArray(x, &xv));
}

PetscParVector &PetscParVector::operator=(PetscScalar s)
{
  PalacePetscCall(VecSet(x, s));
  return *this;
}

void PetscParVector::Scale(PetscScalar s)
{
  PalacePetscCall(VecScale(x, s));
}

void PetscParVector::Shift(PetscScalar s)
{
  PalacePetscCall(VecShift(x, s));
}

void PetscParVector::Abs()
{
  PalacePetscCall(VecAbs(x));
}

void PetscParVector::SqrtAbs()
{
  PalacePetscCall(VecSqrtAbs(x));
}

void PetscParVector::Inv()
{
  PalacePetscCall(VecReciprocal(x));
}

void PetscParVector::InvSqrt()
{
  PalacePetscCall(VecPow(x, -0.5));
}

#if defined(PETSC_USE_COMPLEX)
void PetscParVector::Conj()
{
  PalacePetscCall(VecConjugate(x));
}

void PetscParVector::GetRealPart()
{
  PalacePetscCall(VecRealPart(x));
}

void PetscParVector::GetImagPart()
{
  PalacePetscCall(VecImaginaryPart(x));
}
#endif

PetscReal PetscParVector::Normalize()
{
  PetscReal norm;
  PalacePetscCall(VecNormalize(x, &norm));
  return norm;
}

PetscReal PetscParVector::Normalize(const PetscParMatrix &B, PetscParVector &Bx)
{
  B.Mult(*this, Bx);
  PetscReal norm =
      PetscSqrtReal(PetscAbsScalar(Bx.Dot(*this)));  // For SPD B, xᴴ B x is real
  Scale(1.0 / norm);
  return norm;
}

PetscReal PetscParVector::Norml2() const
{
  PetscReal norm;
  PalacePetscCall(VecNorm(x, NORM_2, &norm));
  return norm;
}

PetscReal PetscParVector::Normlinf() const
{
  PetscReal norm;
  PalacePetscCall(VecNorm(x, NORM_INFINITY, &norm));
  return norm;
}

void PetscParVector::ZeroRows(const mfem::Array<int> &rows)
{
  PetscScalar *xv;
  PalacePetscCall(VecGetArray(x, &xv));
  MFEM_FORALL(i, rows.Size(), { xv[rows[i]] = 0.0; });
  PalacePetscCall(VecRestoreArray(x, &xv));
}

void PetscParVector::PointwiseMult(const PetscParVector &y, bool replace_zeros)
{
  MFEM_VERIFY(GetSize() == y.GetSize(), "Invalid size!");
  if (replace_zeros)
  {
    PetscScalar *yv;
    PalacePetscCall(VecGetArray(y, &yv));
    MFEM_FORALL(i, GetSize(), {
      if (yv[i] == 0.0)
      {
        yv[i] = 1.0;
      }
    });
    PalacePetscCall(VecRestoreArray(y, &yv));
  }
  PalacePetscCall(VecPointwiseMult(x, x, y));
}

void PetscParVector::AXPY(PetscScalar alpha, const PetscParVector &y)
{
  MFEM_VERIFY(GetSize() == y.GetSize(), "Invalid size!");
  PalacePetscCall(VecAXPY(x, alpha, y));
}

void PetscParVector::AXPBY(PetscScalar alpha, const PetscParVector &y, PetscScalar beta)
{
  MFEM_VERIFY(GetSize() == y.GetSize(), "Invalid size!");
  PalacePetscCall(VecAXPBY(x, alpha, beta, y));
}

void PetscParVector::AXPBYPCZ(PetscScalar alpha, const PetscParVector &y, PetscScalar beta,
                              const PetscParVector &z, PetscScalar gamma)
{
  MFEM_VERIFY(GetSize() == y.GetSize() && GetSize() == z.GetSize(), "Invalid size!");
  PalacePetscCall(VecAXPBYPCZ(x, alpha, beta, gamma, y, z));
}

PetscScalar PetscParVector::Dot(const PetscParVector &y) const
{
  PetscScalar val;
  PalacePetscCall(VecDot(x, y, &val));
  return val;
}

PetscScalar PetscParVector::TransposeDot(const PetscParVector &y) const
{
  PetscScalar val;
  PalacePetscCall(VecTDot(x, y, &val));
  return val;
}

void PetscParVector::Print(const char *fname, bool binary) const
{
  if (fname)
  {
    PetscViewer view;
    if (binary)
    {
      PalacePetscCall(
          PetscViewerBinaryOpen(PetscObjectComm(reinterpret_cast<PetscObject>(x)), fname,
                                FILE_MODE_WRITE, &view));
    }
    else
    {
      PalacePetscCall(PetscViewerASCIIOpen(
          PetscObjectComm(reinterpret_cast<PetscObject>(x)), fname, &view));
    }
    PalacePetscCall(VecView(x, view));
    PalacePetscCall(PetscViewerDestroy(&view));
  }
  else
  {
    PalacePetscCall(VecView(x, nullptr));
  }
}

MPI_Comm PetscParVector::GetComm() const
{
  return x ? PetscObjectComm(reinterpret_cast<PetscObject>(x)) : MPI_COMM_NULL;
}

// PetscParMatrix methods

PetscParMatrix::PetscParMatrix(const PetscParMatrix &B)
{
  PalacePetscCall(MatDuplicate(B, MAT_COPY_VALUES, &A));
}

PetscParMatrix::PetscParMatrix(Mat B, bool ref)
{
  A = B;
  if (ref)
  {
    PalacePetscCall(PetscObjectReference(reinterpret_cast<PetscObject>(B)));
  }
}

PetscParMatrix::~PetscParMatrix()
{
  MPI_Comm comm;
  PalacePetscCall(PetscObjectGetComm(reinterpret_cast<PetscObject>(A), &comm));
  PalacePetscCall(MatDestroy(&A));
}

void PetscParMatrix::SetSymmetric(bool sym)
{
  PalacePetscCall(MatSetOption(A, MAT_SYMMETRIC, sym ? PETSC_TRUE : PETSC_FALSE));
  PalacePetscCall(MatSetOption(A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
}

void PetscParMatrix::SetHermitian(bool herm)
{
  PalacePetscCall(MatSetOption(A, MAT_HERMITIAN, herm ? PETSC_TRUE : PETSC_FALSE));
  PalacePetscCall(MatSetOption(A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
}

bool PetscParMatrix::GetSymmetric() const
{
  PetscBool flg, sym;
  PalacePetscCall(MatIsSymmetricKnown(A, &flg, &sym));
  return (flg == PETSC_TRUE && sym == PETSC_TRUE);
}

bool PetscParMatrix::GetHermitian() const
{
  PetscBool flg, herm;
  PalacePetscCall(MatIsHermitianKnown(A, &flg, &herm));
  return (flg == PETSC_TRUE && herm == PETSC_TRUE);
}

#if defined(PETSC_USE_COMPLEX)
void PetscParMatrix::SetRealSymmetric()
{
  PalacePetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE));
  PalacePetscCall(MatSetOption(A, MAT_HERMITIAN, PETSC_TRUE));
  PalacePetscCall(MatSetOption(A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
}
#endif

void PetscParMatrix::CopySymmetry(const PetscParMatrix &B)
{
  PalacePetscCall(MatPropagateSymmetryOptions(B, A));
}

PetscInt PetscParMatrix::GetNumRows() const
{
  PetscInt m;
  PalacePetscCall(MatGetLocalSize(A, &m, nullptr));
  return m;
}

PetscInt PetscParMatrix::GetNumCols() const
{
  PetscInt n;
  PalacePetscCall(MatGetLocalSize(A, nullptr, &n));
  return n;
}

PetscInt PetscParMatrix::GetGlobalNumRows() const
{
  PetscInt M;
  PalacePetscCall(MatGetSize(A, &M, nullptr));
  return M;
}

PetscInt PetscParMatrix::GetGlobalNumCols() const
{
  PetscInt N;
  PalacePetscCall(MatGetSize(A, nullptr, &N));
  return N;
}

PetscInt PetscParMatrix::NNZ() const
{
  MatInfo info;
  PalacePetscCall(MatGetInfo(A, MAT_GLOBAL_SUM, &info));
  return (PetscInt)info.nz_used;
}

PetscReal PetscParMatrix::NormF() const
{
  PetscReal norm;
  PalacePetscCall(MatNorm(A, NORM_FROBENIUS, &norm));
  return norm;
}

PetscReal PetscParMatrix::NormInf() const
{
  PetscReal norm;
  PalacePetscCall(MatNorm(A, NORM_INFINITY, &norm));
  return norm;
}

PetscReal PetscParMatrix::Norm2(PetscReal tol, PetscInt maxits) const
{
  // XX TODO: Add separate if condition using ARPACK estimate before reverting to power
  //          iteration.
  if (tol == PETSC_DEFAULT)
  {
    tol = 1.0e-4;
  }
  if (maxits == PETSC_DEFAULT)
  {
    maxits = 100;
  }
#if defined(PALACE_WITH_SLEPC)
  return slepc::GetMaxSingularValue(*this, tol, maxits);
#else
  // Power iteration loop: ||A||₂² = λₙ(Aᴴ A) .
  PetscInt it = 0;
  PetscReal res = 0.0;
  PetscReal l, l0 = 0.0;
  PetscParVector u(*this), v(*this);
  u.SetRandom();
  u.Normalize();
  while (it < maxits)
  {
    Mult(u, v);
    if (GetHermitian())
    {
      u.Copy(v);
    }
    else
    {
      MultHermitianTranspose(v, u);
    }
    l = u.Normalize();
    if (it > 0)
    {
      res = PetscAbsReal(l - l0) / PetscAbsReal(l0);
      if (res < tol)
      {
        break;
      }
    }
    l0 = l;
    it++;
  }
  if (it >= maxits)
  {
    Mpi::Warning(GetComm(),
                 "Power iteration did not converge in {:d} "
                 "iterations, res = {:.3e}, lambda = {:.3e}!\n",
                 it, res, l);
  }
  return GetHermitian() ? l : PetscSqrtReal(l);
#endif
}

void PetscParMatrix::Scale(PetscScalar s)
{
  PalacePetscCall(MatScale(A, s));
}

#if defined(PETSC_USE_COMPLEX)
void PetscParMatrix::Conj()
{
  PalacePetscCall(MatConjugate(A));
}

void PetscParMatrix::GetRealPart()
{
  PalacePetscCall(MatRealPart(A));
}

void PetscParMatrix::GetImagPart()
{
  PalacePetscCall(MatImaginaryPart(A));
}
#endif

void PetscParMatrix::AXPY(PetscScalar alpha, const PetscParMatrix &B,
                          PetscParMatrix::NNZStructure struc)
{
  switch (struc)
  {
    case NNZStructure::DIFFERENT:
      PalacePetscCall(MatAXPY(A, alpha, B, DIFFERENT_NONZERO_PATTERN));
      break;
    case NNZStructure::SAME:
      PalacePetscCall(MatAXPY(A, alpha, B, SAME_NONZERO_PATTERN));
      break;
    case NNZStructure::SUBSET:
      PalacePetscCall(MatAXPY(A, alpha, B, SUBSET_NONZERO_PATTERN));
      break;
    default:
      MFEM_ABORT("Nonzero type for MatAXPY not implemented!");
      break;
  }
}

void PetscParMatrix::Mult(const PetscParVector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(::MatMult(A, x, y));
}

void PetscParMatrix::MultAdd(const PetscParVector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(MatMultAdd(A, x, y, y));
}

void PetscParMatrix::MultTranspose(const PetscParVector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumRows() && y.GetSize() == GetNumCols(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(::MatMultTranspose(A, (Vec)x, (Vec)y));
}

void PetscParMatrix::MultTransposeAdd(const PetscParVector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumRows() && y.GetSize() == GetNumCols(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(MatMultTransposeAdd(A, x, y, y));
}

void PetscParMatrix::MultHermitianTranspose(const PetscParVector &x,
                                            PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumRows() && y.GetSize() == GetNumCols(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(MatMultHermitianTranspose(A, x, y));
}

void PetscParMatrix::MultHermitianTransposeAdd(const PetscParVector &x,
                                               PetscParVector &y) const
{
  MFEM_VERIFY(x.GetSize() == GetNumRows() && y.GetSize() == GetNumCols(),
              "Incorrect vector sizes for matrix-vector product!");
  PalacePetscCall(MatMultHermitianTransposeAdd(A, x, y, y));
}

#if defined(PETSC_USE_COMPLEX)
void PetscParMatrix::Mult(const mfem::Vector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  PetscParVector xx(GetComm(), x);
  Mult(xx, y);
}

void PetscParMatrix::MultTranspose(const mfem::Vector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  PetscParVector xx(GetComm(), x);
  MultTranspose(xx, y);
}

void PetscParMatrix::MultHermitianTranspose(const mfem::Vector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  PetscParVector xx(GetComm(), x);
  MultHermitianTranspose(xx, y);
}
#endif

void PetscParMatrix::Print(const char *fname, bool binary) const
{
  if (fname)
  {
    PetscViewer view;
    if (binary)
    {
      PalacePetscCall(
          PetscViewerBinaryOpen(PetscObjectComm(reinterpret_cast<PetscObject>(A)), fname,
                                FILE_MODE_WRITE, &view));
    }
    else
    {
      PalacePetscCall(PetscViewerASCIIOpen(
          PetscObjectComm(reinterpret_cast<PetscObject>(A)), fname, &view));
    }
    PalacePetscCall(MatView(A, view));
    PalacePetscCall(PetscViewerDestroy(&view));
  }
  else
  {
    PalacePetscCall(MatView(A, nullptr));
  }
}

std::unique_ptr<mfem::HypreParMatrix>
#if defined(PETSC_USE_COMPLEX)
PetscParMatrix::GetHypreParMatrix(PetscParMatrix::ExtractStructure struc) const
#else
PetscParMatrix::GetHypreParMatrix() const
#endif
{
  HYPRE_BigInt M = GetGlobalNumRows();
  HYPRE_BigInt N = GetGlobalNumCols();
  std::unique_ptr<HYPRE_BigInt[]> rows, cols;
  if (HYPRE_AssumedPartitionCheck())
  {
    PetscInt start, end;
    rows = std::make_unique<HYPRE_BigInt[]>(2);
    PalacePetscCall(MatGetOwnershipRange(A, &start, &end));
    rows[0] = start;
    rows[1] = end;
    if (M != N)
    {
      cols = std::make_unique<HYPRE_BigInt[]>(2);
      PalacePetscCall(MatGetOwnershipRangeColumn(A, &start, &end));
      cols[0] = start;
      cols[1] = end;
    }
  }
  else
  {
    PetscMPIInt comm_size;
    const PetscInt *ranges;
    MPI_Comm_size(GetComm(), &comm_size);
    rows = std::make_unique<HYPRE_BigInt[]>(comm_size + 1);
    PalacePetscCall(MatGetOwnershipRanges(A, &ranges));
    for (PetscMPIInt i = 0; i < comm_size + 1; i++)
    {
      rows[i] = ranges[i];
    }
    if (M != N)
    {
      cols = std::make_unique<HYPRE_BigInt[]>(comm_size + 1);
      PalacePetscCall(MatGetOwnershipRangesColumn(A, &ranges));
      for (PetscMPIInt i = 0; i < comm_size + 1; i++)
      {
        cols[i] = ranges[i];
      }
    }
  }

  // Count nonzeros.
  MatInfo info;
  PalacePetscCall(MatGetInfo(A, MAT_LOCAL, &info));
  PetscInt nnz = (PetscInt)info.nz_used;

  // Copy local CSR block of rows (columns in global numbering).
  PetscInt rstart, rend, n;
  const PetscInt *jj;
  const PetscScalar *vals;
  PalacePetscCall(MatGetOwnershipRange(A, &rstart, &rend));

  int m = rend - rstart;
  std::unique_ptr<int[]> II = std::make_unique<int[]>(m + 1);
  std::unique_ptr<HYPRE_BigInt[]> JJ = std::make_unique<HYPRE_BigInt[]>(nnz);
  std::unique_ptr<double[]> data = std::make_unique<double[]>(nnz);
  nnz = 0;

  for (PetscInt i = rstart; i < rend; i++)
  {
    PalacePetscCall(MatGetRow(A, i, &n, &jj, &vals));
    II[i - rstart] = nnz;
    for (PetscInt j = 0; j < n; j++)
    {
#if defined(PETSC_USE_COMPLEX)
      if (struc == ExtractStructure::REAL)
      {
        data[nnz] = PetscRealPart(vals[j]);
      }
      else if (struc == ExtractStructure::IMAGINARY)
      {
        data[nnz] = PetscImaginaryPart(vals[j]);
      }
      else  // struc == ExtractStructure::SUM
      {
        data[nnz] = PetscRealPart(vals[j]) + PetscImaginaryPart(vals[j]);
      }
#else
      data[nnz] = vals[j];
#endif
      JJ[nnz++] = jj[j];
    }
    PalacePetscCall(MatRestoreRow(A, i, &n, &jj, &vals));
  }
  II[m] = nnz;

  // Create the HypreParMatrix (copies all inputs so memory of local variables is released
  // after return).
  if (M == N)
  {
    return std::make_unique<mfem::HypreParMatrix>(GetComm(), m, M, N, II.get(), JJ.get(),
                                                  data.get(), rows.get(), rows.get());
  }
  else
  {
    return std::make_unique<mfem::HypreParMatrix>(GetComm(), m, M, N, II.get(), JJ.get(),
                                                  data.get(), rows.get(), cols.get());
  }
}

PetscErrorCode Convert_Array_IS(MPI_Comm comm, bool islist, const mfem::Array<int> &list,
                                PetscInt start, IS *is)
{
  // Converts from a list (or a marked Array if islist is false) to an IS. The offset where
  // to start numbering is given as start.
  PetscInt n = list.Size(), *idxs;
  const auto *data = list.HostRead();
  PetscFunctionBeginUser;

  if (islist)
  {
    PetscCall(PetscMalloc1(n, &idxs));
    for (PetscInt i = 0; i < n; i++)
    {
      idxs[i] = data[i] + start;
    }
  }
  else
  {
    PetscInt cum = 0;
    for (PetscInt i = 0; i < n; i++)
    {
      if (data[i])
      {
        cum++;
      }
    }
    PetscCall(PetscMalloc1(cum, &idxs));
    cum = 0;
    for (PetscInt i = 0; i < n; i++)
    {
      if (data[i])
      {
        idxs[cum++] = i + start;
      }
    }
    n = cum;
  }
  PetscCall(ISCreateGeneral(comm, n, idxs, PETSC_OWN_POINTER, is));
  PetscFunctionReturn(0);
}

std::unique_ptr<PetscParMatrix> PetscParMatrix::GetSubMatrix(const mfem::Array<int> &rows,
                                                             const mfem::Array<int> &cols)
{
  PetscInt rst, cst;
  IS row_is, col_is;
  Mat B;
  PalacePetscCall(MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE));
  // Rows need to be in global numbering.
  PalacePetscCall(MatGetOwnershipRange(A, &rst, nullptr));
  PalacePetscCall(MatGetOwnershipRange(A, &cst, nullptr));
  PalacePetscCall(Convert_Array_IS(GetComm(), true, rows, rst, &row_is));
  PalacePetscCall(Convert_Array_IS(GetComm(), true, cols, cst, &col_is));
  PalacePetscCall(MatCreateSubMatrix(A, row_is, col_is, MAT_INITIAL_MATRIX, &B));
  PalacePetscCall(ISDestroy(&row_is));
  PalacePetscCall(ISDestroy(&col_is));
  return std::make_unique<PetscParMatrix>(B, false);
}

std::unique_ptr<PetscParMatrix> PetscParMatrix::GetSequentialMatrix(bool create)
{
  IS row_is, col_is;
  PetscInt nmat = create ? 1 : 0;
  Mat *pB = nullptr, B = nullptr;
  if (create)
  {
    PetscInt M = GetGlobalNumRows(), N = GetGlobalNumCols();
    PalacePetscCall(ISCreateStride(PETSC_COMM_SELF, M, 0, 1, &row_is));
    PalacePetscCall(ISCreateStride(PETSC_COMM_SELF, N, 0, 1, &col_is));
  }
  PalacePetscCall(MatCreateSubMatrices(A, nmat, &row_is, &col_is, MAT_INITIAL_MATRIX, &pB));
  if (create)
  {
    PalacePetscCall(ISDestroy(&row_is));
    PalacePetscCall(ISDestroy(&col_is));
    B = pB[0];
  }
  PalacePetscCall(PetscFree(pB));
  return (B) ? std::make_unique<PetscParMatrix>(B, false) : nullptr;
}

MPI_Comm PetscParMatrix::GetComm() const
{
  return A ? PetscObjectComm(reinterpret_cast<PetscObject>(A)) : MPI_COMM_NULL;
}

PetscShellMatrix::PetscShellMatrix(MPI_Comm comm, std::unique_ptr<mfem::Operator> &&B)
{
  // Wrap the MFEM Operator as a PETSc shell, which inherets the underlying matrix storage
  // (when the PETSc matrix is destroyed, so is the Hypre one).
  MFEM_VERIFY(B, "Cannot construct PETSc shell from an empty matrix!");
  PetscInt m = (PetscInt)B->Height();
  PetscInt n = (PetscInt)B->Width();

  PetscMatShellCtx *ctx = new PetscMatShellCtx;
  ctx->Ar = std::move(B);
#if defined(PETSC_USE_COMPLEX)
  ctx->Ai = nullptr;
  ctx->x.SetSize(2 * n);
  ctx->y.SetSize(2 * m);
#else
  ctx->x.SetSize(n);
  ctx->y.SetSize(m);
#endif

  PalacePetscCall(MatCreateShell(comm, m, n, PETSC_DECIDE, PETSC_DECIDE, (void *)ctx, &A));
  __mat_shell_init(A);
}

#if defined(PETSC_USE_COMPLEX)
PetscShellMatrix::PetscShellMatrix(MPI_Comm comm, std::unique_ptr<mfem::Operator> &&Br,
                                   std::unique_ptr<mfem::Operator> &&Bi)
{
  MFEM_VERIFY(Br || Bi, "Cannot construct PETSc shell from an empty matrix!");
  MFEM_VERIFY((!Br || !Bi) || (Br->Height() == Bi->Height() && Br->Width() == Bi->Width()),
              "Mismatch in dimension of real and imaginary matrix parts!");
  PetscInt m, n;
  if (Br)
  {
    m = (PetscInt)Br->Height();
    n = (PetscInt)Br->Width();
  }
  else
  {
    m = (PetscInt)Bi->Height();
    n = (PetscInt)Bi->Width();
  }

  PetscMatShellCtx *ctx = new PetscMatShellCtx;
  ctx->Ar = std::move(Br);
  ctx->Ai = std::move(Bi);
  ctx->x.SetSize(2 * n);
  ctx->y.SetSize(2 * m);

  PalacePetscCall(MatCreateShell(comm, m, n, PETSC_DECIDE, PETSC_DECIDE, (void *)ctx, &A));
  __mat_shell_init(A);
}
#endif

PetscMatShellCtx *PetscShellMatrix::GetContext() const
{
  PetscMatShellCtx *ctx;
  PalacePetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  return ctx;
}

PetscInt PetscShellMatrix::NNZ() const
{
  HYPRE_BigInt nnz;
  PetscMatShellCtx *ctx = GetContext();
#if defined(PETSC_USE_COMPLEX)
  MFEM_VERIFY(!(ctx->Ar && ctx->Ai), "Use NNZReal/NNZImag methods for complex matrices!");
  nnz = (ctx->Ar) ? dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).NNZ()
                  : ((ctx->Ai) ? dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai).NNZ() : 0);
#else
  nnz = (ctx->Ar) ? dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).NNZ() : 0;
#endif
  return (PetscInt)nnz;
}

#if defined(PETSC_USE_COMPLEX)
PetscInt PetscShellMatrix::NNZReal() const
{
  HYPRE_BigInt nnz;
  PetscMatShellCtx *ctx = GetContext();
  nnz = (ctx->Ar) ? dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).NNZ() : 0;
  return (PetscInt)nnz;
}

PetscInt PetscShellMatrix::NNZImag() const
{
  HYPRE_BigInt nnz;
  PetscMatShellCtx *ctx = GetContext();
  nnz = (ctx->Ai) ? dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai).NNZ() : 0;
  return (PetscInt)nnz;
}
#endif

#if defined(PETSC_USE_COMPLEX)
PetscReal PetscShellMatrix::NormFReal() const
{
  HYPRE_Real norm;
  PetscMatShellCtx *ctx = GetContext();
  norm = (ctx->Ar) ? hypre_ParCSRMatrixFnorm(dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar))
                   : 0.0;
  return norm;
}

PetscReal PetscShellMatrix::NormFImag() const
{
  HYPRE_Real norm;
  PetscMatShellCtx *ctx = GetContext();
  norm = (ctx->Ai) ? hypre_ParCSRMatrixFnorm(dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai))
                   : 0.0;
  return norm;
}

PetscReal PetscShellMatrix::NormInfReal() const
{
  HYPRE_Real norm;
  PetscMatShellCtx *ctx = GetContext();
  if (ctx->Ar)
  {
    hypre_ParCSRMatrixInfNorm(dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar), &norm);
  }
  else
  {
    norm = 0.0;
  }
  return norm;
}

PetscReal PetscShellMatrix::NormInfImag() const
{
  HYPRE_Real norm;
  PetscMatShellCtx *ctx = GetContext();
  if (ctx->Ai)
  {
    hypre_ParCSRMatrixInfNorm(dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai), &norm);
  }
  else
  {
    norm = 0.0;
  }
  return norm;
}
#endif

#if defined(PETSC_USE_COMPLEX)
void PetscShellMatrix::Mult(const mfem::Vector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  __mat_shell_apply(A, x, y);
}

void PetscShellMatrix::MultTranspose(const mfem::Vector &x, PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  __mat_shell_apply_transpose(A, x, y);
}

void PetscShellMatrix::MultHermitianTranspose(const mfem::Vector &x,
                                              PetscParVector &y) const
{
  MFEM_VERIFY(x.Size() == GetNumCols() && y.GetSize() == GetNumRows(),
              "Incorrect vector sizes for matrix-vector product!");
  __mat_shell_apply_hermitian_transpose(A, x, y);
}
#endif

void PetscShellMatrix::Print(const char *fname, bool binary) const
{
  MFEM_VERIFY(
      fname && !binary,
      "PetscShellMatrix::Print only works with a specified filename and binary = false!")
  PetscMatShellCtx *ctx = GetContext();
#if defined(PETSC_USE_COMPLEX)
  MFEM_VERIFY(!(ctx->Ar && ctx->Ai),
              "Use PrintReal/PrintImag methods for complex matrices!");
  if (ctx->Ar)
  {
    dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).Print(fname);
  }
  else if (ctx->Ai)
  {
    dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai).Print(fname);
  }
#else
  if (ctx->Ar)
  {
    dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).Print(fname);
  }
#endif
}

#if defined(PETSC_USE_COMPLEX)
void PetscShellMatrix::PrintReal(const char *fname) const
{
  PetscMatShellCtx *ctx = GetContext();
  if (ctx->Ar)
  {
    dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ar).Print(fname);
  }
}

void PetscShellMatrix::PrintImag(const char *fname) const
{
  PetscMatShellCtx *ctx = GetContext();
  if (ctx->Ai)
  {
    dynamic_cast<mfem::HypreParMatrix &>(*ctx->Ai).Print(fname);
  }
}
#endif

#if defined(PETSC_USE_COMPLEX)
bool PetscShellMatrix::HasReal() const
{
  PetscMatShellCtx *ctx = GetContext();
  return (ctx->Ar != nullptr);
}

bool PetscShellMatrix::HasImag() const
{
  PetscMatShellCtx *ctx = GetContext();
  return (ctx->Ai != nullptr);
}
#endif

const mfem::Operator *
#if defined(PETSC_USE_COMPLEX)
PetscShellMatrix::GetOperator(PetscParMatrix::ExtractStructure struc) const
#else
PetscShellMatrix::GetOperator() const
#endif
{
  PetscMatShellCtx *ctx = GetContext();
#if defined(PETSC_USE_COMPLEX)
  if (struc == ExtractStructure::REAL)
  {
    MFEM_VERIFY(ctx->Ar, "Invalid use of GetOperator, no real matrix component defined!");
    return ctx->Ar.get();
  }
  else if (struc == ExtractStructure::IMAGINARY)
  {
    MFEM_VERIFY(ctx->Ai,
                "Invalid use of GetOperator, no imaginary matrix component defined!");
    return ctx->Ai.get();
  }
  MFEM_ABORT("ExtractStructure::SUM is not implemented for PetscShellMatrix!");
  return nullptr;
#else
  MFEM_VERIFY(ctx->Ar, "Invalid use of GetOperator, no matrix defined!");
  return ctx->Ar.get();
#endif
}

PetscAijMatrix::PetscAijMatrix(const mfem::Operator &B)
{
  auto hB = dynamic_cast<const mfem::HypreParMatrix *>(&B);
  MFEM_VERIFY(hB, "PetscAijMatrix constructor requires Operator of type HypreParMatrix!");
  PalacePetscCall(__mat_convert_hypreParCSR_AIJ(*hB, &A));
}

#if defined(PETSC_USE_COMPLEX)
PetscAijMatrix::PetscAijMatrix(const mfem::Operator &Br, const mfem::Operator &Bi)
{
  Mat Ai;
  auto hBr = dynamic_cast<const mfem::HypreParMatrix *>(&Br);
  auto hBi = dynamic_cast<const mfem::HypreParMatrix *>(&Bi);
  MFEM_VERIFY(hBr && hBi,
              "PetscAijMatrix constructor requires Operator of type HypreParMatrix!");
  PalacePetscCall(__mat_convert_hypreParCSR_AIJ(*hBr, &A));
  PalacePetscCall(__mat_convert_hypreParCSR_AIJ(*hBi, &Ai));
  PalacePetscCall(MatAXPY(A, PETSC_i, Ai, UNKNOWN_NONZERO_PATTERN));
  PalacePetscCall(MatDestroy(&Ai));
}
#endif

PetscDenseMatrix::PetscDenseMatrix(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt M,
                                   PetscInt N, PetscScalar *data)
{
  PalacePetscCall(MatCreateDense(comm, m, n, M, N, data, &A));
}

PetscDenseMatrix::PetscDenseMatrix(PetscInt m, PetscInt n, PetscScalar *data)
{
  PalacePetscCall(MatCreateSeqDense(PETSC_COMM_SELF, m, n, data, &A));
}

void PetscDenseMatrix::Resize(PetscInt m, PetscInt n, bool copy)
{
  Mat B;
  PetscScalar *Aj, *Bj;
  PetscInt m0 = GetNumRows(), n0 = GetNumCols();
  if (m0 == m && n0 == n)
  {
    return;
  }
  PalacePetscCall(MatCreateDense(GetComm(), m, n, PETSC_DECIDE, PETSC_DECIDE, nullptr, &B));
  if (copy)
  {
    for (PetscInt j = 0; j < std::min(n, n0); j++)
    {
      PalacePetscCall(MatDenseGetColumn(A, j, &Aj));
      PalacePetscCall(MatDenseGetColumn(B, j, &Bj));
      for (PetscInt i = 0; i < std::min(m, m0); i++)
      {
        Bj[i] = Aj[i];
      }
      PalacePetscCall(MatDenseRestoreColumn(A, &Aj));
      PalacePetscCall(MatDenseRestoreColumn(B, &Bj));
    }
  }
  PalacePetscCall(MatPropagateSymmetryOptions(A, B));
  PalacePetscCall(MatDestroy(&A));
  A = B;
}

PetscParVector PetscDenseMatrix::GetColumn(PetscInt j)
{
  MFEM_VERIFY(j >= 0 && j < GetGlobalNumCols(), "Column index out of range!");
  Vec v;
  PalacePetscCall(MatDenseGetColumnVec(A, j, &v));
  return PetscParVector(v, true);
}

const PetscParVector PetscDenseMatrix::GetColumnRead(PetscInt j) const
{
  MFEM_VERIFY(j >= 0 && j < GetGlobalNumCols(), "Column index out of range!");
  Vec v;
  PalacePetscCall(MatDenseGetColumnVecRead(A, j, &v));
  return PetscParVector(v, true);
}

void PetscDenseMatrix::RestoreColumn(PetscInt j, PetscParVector &v)
{
  MFEM_VERIFY(j >= 0 && j < GetGlobalNumCols(), "Column index out of range!");
  Vec u = v;
  PalacePetscCall(MatDenseRestoreColumnVec(A, j, &u));
}

void PetscDenseMatrix::RestoreColumnRead(PetscInt j, const PetscParVector &v) const
{
  MFEM_VERIFY(j >= 0 && j < GetGlobalNumCols(), "Column index out of range!");
  Vec u = v;
  PalacePetscCall(MatDenseRestoreColumnVecRead(A, j, &u));
}

PetscScalar *PetscDenseMatrix::GetArray()
{
  PetscScalar *data;
  PalacePetscCall(MatDenseGetArray(A, &data));
  return data;
}

const PetscScalar *PetscDenseMatrix::GetArrayRead() const
{
  const PetscScalar *data;
  PalacePetscCall(MatDenseGetArrayRead(A, &data));
  return data;
}

void PetscDenseMatrix::RestoreArray(PetscScalar *data)
{
  PalacePetscCall(MatDenseRestoreArray(A, &data));
}

void PetscDenseMatrix::RestoreArrayRead(const PetscScalar *data) const
{
  PalacePetscCall(MatDenseRestoreArrayRead(A, &data));
}

void PetscDenseMatrix::SetRandom(PetscInt start, PetscInt end)
{
  PetscRandom rand;
  MPI_Comm comm = GetComm();
  PalacePetscCall(PetscRandomCreate(comm, &rand));
#if defined(PETSC_USE_COMPLEX)
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0 - PETSC_i, 1.0 + PETSC_i));
#else
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0, 1.0));
#endif
  if (start < 0)
  {
    start = 0;
  }
  if (end < 0)
  {
    end = GetGlobalNumCols();
  }
  MFEM_VERIFY(0 <= start && start <= end && end <= GetGlobalNumCols(),
              "Invalid start/end columns for SetRandom!");
  for (PetscInt j = start; j < end; j++)
  {
    PetscParVector v = GetColumn(j);
    PalacePetscCall(VecSetRandom(v, rand));
    RestoreColumn(j, v);
  }
  PalacePetscCall(PetscRandomDestroy(&rand));
}

#if defined(PETSC_USE_COMPLEX)
void PetscDenseMatrix::SetRandomReal(PetscInt start, PetscInt end)
{
  PetscRandom rand;
  MPI_Comm comm = GetComm();
  PalacePetscCall(PetscRandomCreate(comm, &rand));
  PalacePetscCall(PetscRandomSetInterval(rand, -1.0, 1.0));
  if (start < 0)
  {
    start = 0;
  }
  if (end < 0)
  {
    end = GetGlobalNumCols();
  }
  MFEM_VERIFY(0 <= start && start <= end && end <= GetGlobalNumCols(),
              "Invalid start/end columns for SetRandom!");
  for (PetscInt j = start; j < end; j++)
  {
    PetscParVector v = GetColumn(j);
    PalacePetscCall(VecSetRandom(v, rand));
    RestoreColumn(j, v);
  }
  PalacePetscCall(PetscRandomDestroy(&rand));
}
#endif

void PetscDenseMatrix::SetRandomSign(PetscInt start, PetscInt end, bool init)
{
  if (start < 0)
  {
    start = 0;
  }
  if (end < 0)
  {
    end = GetGlobalNumCols();
  }
  MFEM_VERIFY(0 <= start && start <= end && end <= GetGlobalNumCols(),
              "Invalid start/end columns for SetRandom!");
  if (!init)
  {
    SetRandomReal(start, end);
  }
  for (PetscInt j = start; j < end; j++)
  {
    PetscParVector v = GetColumn(j);
    v.SetRandomSign(true);
    RestoreColumn(j, v);
  }
}

PetscReal PetscDenseMatrix::OrthonormalizeColumn(PetscInt j, bool mgs, bool cgs2)
{
  auto Dot = [](const PetscParVector &v, const PetscParVector &w) -> PetscScalar
  { return v.Dot(w); };
  auto VecDot = [](const PetscParVector &v, const PetscParMatrix &A,
                   PetscParVector &dot) -> void { A.MultHermitianTranspose(v, dot); };
  auto Normalize = [](PetscParVector &v) -> PetscReal { return v.Normalize(); };
  return OrthonormalizeColumnInternal(j, mgs, cgs2, Dot, VecDot, Normalize);
}

PetscReal PetscDenseMatrix::OrthonormalizeColumn(PetscInt j, bool mgs, bool cgs2,
                                                 const PetscParMatrix &B,
                                                 PetscParVector &Bv)
{
  MFEM_VERIFY(Bv.GetSize() == B.GetNumRows(),
              "Workspace error for B-matrix orthonormalization!");
  auto Dot = [&B, &Bv](const PetscParVector &v, const PetscParVector &w) -> PetscScalar
  {
    B.Mult(v, Bv);
    return Bv.Dot(w);
  };
  auto VecDot = [&B, &Bv](const PetscParVector &v, const PetscParMatrix &A,
                          PetscParVector &dot) -> void
  {
    B.Mult(v, Bv);
    A.MultHermitianTranspose(Bv, dot);
  };
  auto Normalize = [&B, &Bv](PetscParVector &v) -> PetscReal { return v.Normalize(B, Bv); };
  return OrthonormalizeColumnInternal(j, mgs, cgs2, Dot, VecDot, Normalize);
}

PetscReal PetscDenseMatrix::OrthonormalizeColumnInternal(
    PetscInt j, bool mgs, bool cgs2,
    const std::function<PetscScalar(PetscParVector &, PetscParVector &)> &Dot,
    const std::function<void(PetscParVector &, PetscDenseMatrix &, PetscParVector &)>
        &VecDot,
    const std::function<PetscReal(PetscParVector &)> &Normalize)
{
  MFEM_VERIFY(j >= 0 && j < GetGlobalNumCols(), "Column index out of range!");
  PetscParVector v = GetColumn(j);
  if (j > 0)
  {
    if (mgs)
    {
      // We can't call GetColumn twice.
      PetscScalar *pA = GetArray();
      for (int i = 0; i < j; i++)
      {
        PetscParVector w(GetComm(), GetNumRows(), PETSC_DECIDE, pA + i * GetNumRows());
        PetscScalar dot = Dot(v, w);
        v.AXPY(-dot, w);
      }
      RestoreArray(pA);
    }
    else
    {
      int refine = (cgs2) ? 2 : 1;
      PetscScalar *pA = GetArray();
      for (int l = 0; l < refine; l++)
      {
        PetscDenseMatrix Aj(GetComm(), GetNumRows(), PETSC_DECIDE, PETSC_DECIDE, j, pA);
        PetscParVector dot(Aj);
        VecDot(v, Aj, dot);
        dot.Scale(-1.0);
        Aj.MultAdd(dot, v);
      }
      RestoreArray(pA);
    }
  }
  PetscReal norm = Normalize(v);
  MFEM_VERIFY(norm > 0.0,
              "Linearly dependent column encountered during vector orthonormalization!");
  RestoreColumn(j, v);
  // {
  //   // Debug
  //   Mpi::Print(GetComm(), "Orthogonality error (j = {:d}):\n", j);
  //   for (int ii = 0; ii <= j; ii++)
  //   {
  //     PetscParVector vv = GetColumn(ii);
  //     PetscScalar err = Dot(vv, vv);
  //     Mpi::Print(GetComm(), "  ({:d}, {:d}): {:e}{:+e}i\n", ii, ii, PetscRealPart(err),
  //                 PetscImaginaryPart(err));
  //     PetscScalar *pA = GetArray();
  //     for (int jj = ii + 1; jj <= j; jj++)
  //     {
  //       // We can't call GetColumn twice.
  //       PetscParVector ww(GetComm(), GetNumRows(), PETSC_DECIDE, pA + jj * GetNumRows());
  //       err = Dot(vv, ww);
  //       Mpi::Print(GetComm(), "  ({:d}, {:d}): {:e}{:+e}i\n", ii, jj, PetscRealPart(err),
  //                   PetscImaginaryPart(err));
  //     }
  //     RestoreArray(pA);
  //     RestoreColumn(ii, vv);
  //   }
  // }
  return norm;
}

void PetscDenseMatrix::MatMult(const PetscDenseMatrix &X, PetscDenseMatrix &Y) const
{
  MFEM_VERIFY(X.GetNumRows() == GetNumCols() && Y.GetNumRows() == GetNumRows(),
              "Incorrect matrix sizes for matrix-matrix product!");
  MFEM_VERIFY(Mpi::Size(GetComm()) == 1,
              "PetscDenseMatrix::MatMult is only implemented for sequential "
              "matrices!");
  const PetscScalar *pA, *pX;
  PetscScalar *pY;
  PetscInt lda;
  PetscBLASInt m, k, n, ldaA, ldaX, ldaY;
  PetscScalar One = 1.0, Zero = 0.0;
  PetscBLASIntCast(Y.GetNumRows(), &m);
  PetscBLASIntCast(Y.GetNumCols(), &n);
  PetscBLASIntCast(GetNumCols(), &k);

  PalacePetscCall(MatDenseGetLDA(A, &lda));
  PetscBLASIntCast(lda, &ldaA);
  PalacePetscCall(MatDenseGetLDA(X, &lda));
  PetscBLASIntCast(lda, &ldaX);
  PalacePetscCall(MatDenseGetLDA(Y, &lda));
  PetscBLASIntCast(lda, &ldaY);

  PalacePetscCall(MatDenseGetArrayRead(A, &pA));
  PalacePetscCall(MatDenseGetArrayRead(X, &pX));
  PalacePetscCall(MatDenseGetArrayWrite(Y, &pY));
  BLASgemm_("N", "N", &m, &n, &k, &One, pA, &ldaA, pX, &ldaX, &Zero, pY, &ldaY);
  PalacePetscCall(MatDenseRestoreArrayRead(A, &pA));
  PalacePetscCall(MatDenseRestoreArrayRead(X, &pX));
  PalacePetscCall(MatDenseRestoreArrayWrite(Y, &pY));
}

void PetscDenseMatrix::MatMultTranspose(const PetscDenseMatrix &X,
                                        PetscDenseMatrix &Y) const
{
  MFEM_VERIFY(X.GetNumCols() == GetNumCols() && Y.GetNumRows() == GetNumRows(),
              "Incorrect matrix sizes for matrix-matrix product!");
  MFEM_VERIFY(Mpi::Size(GetComm()) == 1,
              "PetscDenseMatrix::MatMultTranspose is only implemented for "
              "sequential matrices!");
  const PetscScalar *pA, *pX;
  PetscScalar *pY;
  PetscInt lda;
  PetscBLASInt m, k, n, ldaA, ldaX, ldaY;
  PetscScalar One = 1.0, Zero = 0.0;
  PetscBLASIntCast(Y.GetNumRows(), &m);
  PetscBLASIntCast(Y.GetNumCols(), &n);
  PetscBLASIntCast(GetNumCols(), &k);

  PalacePetscCall(MatDenseGetLDA(A, &lda));
  PetscBLASIntCast(lda, &ldaA);
  PalacePetscCall(MatDenseGetLDA(X, &lda));
  PetscBLASIntCast(lda, &ldaX);
  PalacePetscCall(MatDenseGetLDA(Y, &lda));
  PetscBLASIntCast(lda, &ldaY);

  PalacePetscCall(MatDenseGetArrayRead(A, &pA));
  PalacePetscCall(MatDenseGetArrayRead(X, &pX));
  PalacePetscCall(MatDenseGetArrayWrite(Y, &pY));
  BLASgemm_("N", "T", &m, &n, &k, &One, pA, &ldaA, pX, &ldaX, &Zero, pY, &ldaY);
  PalacePetscCall(MatDenseRestoreArrayRead(A, &pA));
  PalacePetscCall(MatDenseRestoreArrayRead(X, &pX));
  PalacePetscCall(MatDenseRestoreArrayWrite(Y, &pY));
}

void PetscDenseMatrix::MatTransposeMult(const PetscDenseMatrix &X,
                                        PetscDenseMatrix &Y) const
{
  MFEM_VERIFY(X.GetNumRows() == GetNumRows() && Y.GetNumRows() == GetNumCols(),
              "Incorrect matrix sizes for matrix-matrix product!");
  MFEM_VERIFY(Mpi::Size(GetComm()) == 1,
              "PetscDenseMatrix::MatTransposeMult is only implemented for "
              "sequential matrices!");
  const PetscScalar *pA, *pX;
  PetscScalar *pY;
  PetscInt lda;
  PetscBLASInt m, k, n, ldaA, ldaX, ldaY;
  PetscScalar One = 1.0, Zero = 0.0;
  PetscBLASIntCast(Y.GetNumRows(), &m);
  PetscBLASIntCast(Y.GetNumCols(), &n);
  PetscBLASIntCast(GetNumRows(), &k);

  PalacePetscCall(MatDenseGetLDA(A, &lda));
  PetscBLASIntCast(lda, &ldaA);
  PalacePetscCall(MatDenseGetLDA(X, &lda));
  PetscBLASIntCast(lda, &ldaX);
  PalacePetscCall(MatDenseGetLDA(Y, &lda));
  PetscBLASIntCast(lda, &ldaY);

  PalacePetscCall(MatDenseGetArrayRead(A, &pA));
  PalacePetscCall(MatDenseGetArrayRead(X, &pX));
  PalacePetscCall(MatDenseGetArrayWrite(Y, &pY));
  BLASgemm_("T", "N", &m, &n, &k, &One, pA, &ldaA, pX, &ldaX, &Zero, pY, &ldaY);
  PalacePetscCall(MatDenseRestoreArrayRead(A, &pA));
  PalacePetscCall(MatDenseRestoreArrayRead(X, &pX));
  PalacePetscCall(MatDenseRestoreArrayWrite(Y, &pY));
}

}  // namespace palace::petsc

PetscErrorCode __mat_shell_init(Mat A)
{
  PetscFunctionBeginUser;

  PalacePetscCall(MatShellSetManageScalingShifts(A));
  PalacePetscCall(MatShellSetOperation(A, MATOP_DESTROY, (void (*)())__mat_shell_destroy));
  PetscCall(MatShellSetOperation(
      A, MATOP_MULT,
      (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_shell_apply)));
  PetscCall(
      MatShellSetOperation(A, MATOP_MULT_TRANSPOSE,
                           (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                               &__mat_shell_apply_transpose)));
  PetscCall(
      MatShellSetOperation(A, MATOP_MULT_HERMITIAN_TRANSPOSE,
                           (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                               &__mat_shell_apply_hermitian_transpose)));
  PetscCall(MatShellSetOperation(
      A, MATOP_MULT_ADD,
      (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(&__mat_shell_apply_add)));
  PetscCall(
      MatShellSetOperation(A, MATOP_MULT_TRANSPOSE_ADD,
                           (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                               &__mat_shell_apply_transpose_add)));
  PetscCall(
      MatShellSetOperation(A, MATOP_MULT_HERMITIAN_TRANS_ADD,
                           (void (*)()) static_cast<PetscErrorCode (*)(Mat, Vec, Vec)>(
                               &__mat_shell_apply_hermitian_transpose_add)));
  PetscCall(
      MatShellSetOperation(A, MATOP_GET_DIAGONAL, (void (*)())__mat_shell_get_diagonal));
  // PetscCall(MatShellSetOperation(A, MATOP_SHIFT, (void (*)())__mat_shell_shift));
  // PetscCall(MatShellSetOperation(A, MATOP_SCALE, (void (*)())__mat_shell_scale));
  // PetscCall(MatShellSetOperation(A, MATOP_CONJUGATE, (void (*)())__mat_shell_conj));
  // PetscCall(MatShellSetOperation(A, MATOP_AXPY, (void (*)())__mat_shell_axpy));
  // PetscCall(MatShellSetOperation(A, MATOP_NORM, (void (*)())__mat_shell_norm));
  PetscCall(MatShellSetOperation(A, MATOP_REAL_PART, (void (*)())__mat_shell_real_part));
  PetscCall(
      MatShellSetOperation(A, MATOP_IMAGINARY_PART, (void (*)())__mat_shell_imag_part));
  PetscCall(MatSetUp(A));
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_destroy(Mat A)
{
  palace::petsc::PetscMatShellCtx *ctx;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  delete ctx;
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply_add(Mat A, Vec x, Vec y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
#if defined(PETSC_USE_COMPLEX)
  {
    mfem::Vector xr, xi, yr, yi;
    xr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    xi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    yr.MakeRef(ctx->y, 0, ctx->y.Size() / 2);
    yi.MakeRef(ctx->y, ctx->y.Size() / 2, ctx->y.Size() / 2);
    xx.GetToVectors(xr, xi);
    if (ctx->Ar)
    {
      ctx->Ar->Mult(xr, yr);
      ctx->Ar->Mult(xi, yi);
    }
    else
    {
      yr = 0.0;
      yi = 0.0;
    }
    if (ctx->Ai)
    {
      ctx->Ai->AddMult(xi, yr, -1.0);
      ctx->Ai->AddMult(xr, yi, 1.0);
    }
    yy.AddFromVectors(yr, yi);
  }
#else
  {
    xx.GetToVector(ctx->x);
    if (ctx->Ar)
    {
      ctx->Ar->Mult(ctx->x, ctx->y);
    }
    else
    {
      ctx->y = 0.0;
    }
    yy.AddFromVector(ctx->y);
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply(Mat A, Vec x, Vec y)
{
  PetscCall(VecZeroEntries(y));
  return __mat_shell_apply_add(A, x, y);
}

PetscErrorCode __mat_shell_apply_transpose_add(Mat A, Vec x, Vec y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscBool flg, sym;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  PetscCall(MatIsSymmetricKnown(A, &flg, &sym));
  if (flg == PETSC_TRUE && sym == PETSC_TRUE)
  {
    PetscCall(__mat_shell_apply_add(A, x, y));
    PetscFunctionReturn(0);
  }
#if defined(PETSC_USE_COMPLEX)
  {
    mfem::Vector xr, xi, yr, yi;
    xr.MakeRef(ctx->y, 0, ctx->y.Size() / 2);
    xi.MakeRef(ctx->y, ctx->y.Size() / 2, ctx->y.Size() / 2);
    yr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    yi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    xx.GetToVectors(xr, xi);
    if (ctx->Ar)
    {
      ctx->Ar->MultTranspose(xr, yr);
      ctx->Ar->MultTranspose(xi, yi);
    }
    else
    {
      yr = 0.0;
      yi = 0.0;
    }
    if (ctx->Ai)
    {
      ctx->Ai->AddMultTranspose(xi, yr, -1.0);
      ctx->Ai->AddMultTranspose(xr, yi, 1.0);
    }
    yy.AddFromVectors(yr, yi);
  }
#else
  {
    xx.GetToVector(ctx->y);
    if (ctx->Ar)
    {
      ctx->Ar->MultTranspose(ctx->y, ctx->x);
    }
    else
    {
      ctx->x = 0.0;
    }
    yy.AddFromVector(ctx->x);
  }
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply_transpose(Mat A, Vec x, Vec y)
{
  PetscCall(VecZeroEntries(y));
  return __mat_shell_apply_transpose_add(A, x, y);
}

PetscErrorCode __mat_shell_apply_hermitian_transpose_add(Mat A, Vec x, Vec y)
{
#if defined(PETSC_USE_COMPLEX)
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector xx(x, true), yy(y, true);
  PetscBool flg, sym;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  PetscCall(MatIsHermitianKnown(A, &flg, &sym));
  if (flg == PETSC_TRUE && sym == PETSC_TRUE)
  {
    PetscCall(__mat_shell_apply_add(A, x, y));
    PetscFunctionReturn(0);
  }
  if (!ctx->Ai)
  {
    PetscCall(__mat_shell_apply_transpose_add(A, x, y));
    PetscFunctionReturn(0);
  }
  PetscCall(MatIsSymmetricKnown(A, &flg, &sym));
  {
    mfem::Vector xr, xi, yr, yi;
    xr.MakeRef(ctx->y, 0, ctx->y.Size() / 2);
    xi.MakeRef(ctx->y, ctx->y.Size() / 2, ctx->y.Size() / 2);
    yr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    yi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    xx.GetToVectors(xr, xi);
    if (ctx->Ar)
    {
      if (flg == PETSC_TRUE && sym == PETSC_TRUE)
      {
        ctx->Ar->Mult(xr, yr);
        ctx->Ar->Mult(xi, yi);
      }
      else
      {
        ctx->Ar->MultTranspose(xr, yr);
        ctx->Ar->MultTranspose(xi, yi);
      }
    }
    else
    {
      yr = 0.0;
      yi = 0.0;
    }
    if (ctx->Ai)
    {
      if (flg == PETSC_TRUE && sym == PETSC_TRUE)
      {
        ctx->Ai->AddMult(xi, yr, 1.0);
        ctx->Ai->AddMult(xr, yi, -1.0);
      }
      else
      {
        ctx->Ai->AddMultTranspose(xi, yr, 1.0);
        ctx->Ai->AddMultTranspose(xr, yi, -1.0);
      }
    }
    yy.AddFromVectors(yr, yi);
  }
#else
  PetscCall(__mat_shell_apply_transpose_add(A, x, y));
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply_hermitian_transpose(Mat A, Vec x, Vec y)
{
  PetscCall(VecZeroEntries(y));
  return __mat_shell_apply_hermitian_transpose_add(A, x, y);
}

#if defined(PETSC_USE_COMPLEX)
PetscErrorCode __mat_shell_apply(Mat A, const mfem::Vector &x, Vec y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  {
    mfem::Vector yr, yi;
    yr.MakeRef(ctx->y, 0, ctx->y.Size() / 2);
    yi.MakeRef(ctx->y, ctx->y.Size() / 2, ctx->y.Size() / 2);
    if (ctx->Ar)
    {
      ctx->Ar->Mult(x, yr);
    }
    else
    {
      yr = 0.0;
    }
    if (ctx->Ai)
    {
      ctx->Ai->Mult(x, yi);
    }
    else
    {
      yi = 0.0;
    }
    yy.SetFromVectors(yr, yi);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply_transpose(Mat A, const mfem::Vector &x, Vec y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector yy(y, true);
  PetscBool flg, sym;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  PetscCall(MatIsSymmetricKnown(A, &flg, &sym));
  if (flg == PETSC_TRUE && sym == PETSC_TRUE)
  {
    PetscCall(__mat_shell_apply(A, x, y));
    PetscFunctionReturn(0);
  }
  {
    mfem::Vector yr, yi;
    yr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    yi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    if (ctx->Ar)
    {
      ctx->Ar->MultTranspose(x, yr);
    }
    else
    {
      yr = 0.0;
    }
    if (ctx->Ai)
    {
      ctx->Ai->MultTranspose(x, yi);
    }
    else
    {
      yi = 0.0;
    }
    yy.SetFromVectors(yr, yi);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_apply_hermitian_transpose(Mat A, const mfem::Vector &x, Vec y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector yy(y, true);
  PetscBool flg, sym;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
  PetscCall(MatIsHermitianKnown(A, &flg, &sym));
  if (flg == PETSC_TRUE && sym == PETSC_TRUE)
  {
    PetscCall(__mat_shell_apply(A, x, y));
    PetscFunctionReturn(0);
  }
  if (!ctx->Ai)
  {
    PetscCall(__mat_shell_apply_transpose(A, x, y));
    PetscFunctionReturn(0);
  }
  {
    mfem::Vector yr, yi;
    yr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    yi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    PetscCall(MatIsSymmetricKnown(A, &flg, &sym));
    if (ctx->Ar)
    {
      if (flg == PETSC_TRUE && sym == PETSC_TRUE)
      {
        ctx->Ar->Mult(x, yr);
      }
      else
      {
        ctx->Ar->MultTranspose(x, yr);
      }
    }
    else
    {
      yr = 0.0;
    }
    if (ctx->Ai)
    {
      if (flg == PETSC_TRUE && sym == PETSC_TRUE)
      {
        ctx->Ai->Mult(x, yi);
      }
      else
      {
        ctx->Ai->MultTranspose(x, yi);
      }
      yi.Neg();
    }
    else
    {
      yi = 0.0;
    }
    yy.SetFromVectors(yr, yi);
  }
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode __mat_shell_get_diagonal(Mat A, Vec diag)
{
  palace::petsc::PetscMatShellCtx *ctx;
  palace::petsc::PetscParVector ddiag(diag, true);
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(A, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
#if defined(PETSC_USE_COMPLEX)
  {
    mfem::Vector xr, xi;
    xr.MakeRef(ctx->x, 0, ctx->x.Size() / 2);
    xi.MakeRef(ctx->x, ctx->x.Size() / 2, ctx->x.Size() / 2);
    if (ctx->Ar)
    {
      ctx->Ar->AssembleDiagonal(xr);
    }
    else
    {
      xr = 0.0;
    }
    if (ctx->Ai)
    {
      ctx->Ai->AssembleDiagonal(xi);
    }
    else
    {
      xi = 0.0;
    }
    ddiag.SetFromVectors(xr, xi);
  }
#else
  {
    if (ctx->Ar)
    {
      ctx->Ar->AssembleDiagonal(ctx->x);
    }
    else
    {
      ctx->x = 0.0;
    }
    ddiag.SetFromVector(ctx->x);
  }
#endif
  PetscFunctionReturn(0);
}

// PetscErrorCode __mat_shell_shift(Mat Y, PetscScalar a)
// {
//   palace::petsc::PetscMatShellCtx *ctx;
//   HYPRE_Real as;
//   PetscFunctionBeginUser;

//   PetscCall(MatShellGetContext(Y, (void **)&ctx));
//   MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
//   as = PetscRealPart(a);  // Does nothing if not PETSC_USE_COMPLEX
//   if (std::abs(as) > 0.0)
//   {
//     auto hAr = dynamic_cast<mfem::HypreParMatrix *>(ctx->Ar.get());
//     MFEM_VERIFY(hAr, "Invalid real shift with no real matrix part!");
//     int n = hAr->Height();
//     const hypre_ParCSRMatrix *A = *hAr;
//     const HYPRE_Int *A_diag_i = A->diag->i;
//     HYPRE_Real *A_diag_d = A->diag->data;
//     for (int j = 0; j < n; j++)
//     {
//       A_diag_d[A_diag_i[j]] += as;
//     }
//   }
// #if defined(PETSC_USE_COMPLEX)
//   as = PetscImaginaryPart(a);
//   if (std::abs(as) > 0.0)
//   {
//     auto hAi = dynamic_cast<mfem::HypreParMatrix *>(ctx->Ai.get());
//     MFEM_VERIFY(hAi, "Invalid imaginary shift with no imaginary matrix part!");
//     int n = hAi->Height();
//     const hypre_ParCSRMatrix *A = *hAi;
//     const HYPRE_Int *A_diag_i = A->diag->i;
//     HYPRE_Real *A_diag_d = A->diag->data;
//     for (int j = 0; j < n; j++)
//     {
//       A_diag_d[A_diag_i[j]] += as;
//     }
//   }
// #endif
//   PetscFunctionReturn(0);
// }

// PetscErrorCode __mat_shell_scale(Mat Y, PetscScalar a)
// {
//   palace::petsc::PetscMatShellCtx *ctx;
//   PetscFunctionBeginUser;

//   PetscCall(MatShellGetContext(Y, (void **)&ctx));
//   MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
//   if (a == 0.0)
//   {
//     ctx->Ar.reset();
// #if defined(PETSC_USE_COMPLEX)
//     ctx->Ai.reset();
// #endif
//   }
//   else
//   {
// #if defined(PETSC_USE_COMPLEX)
//     HYPRE_Real ar, ai;
//     ar = PetscRealPart(a);
//     ai = PetscImaginaryPart(a);
//     if (std::abs(ar) > 0.0 && std::abs(ai) == 0.0)
//     {
//       if (ctx->Ar)
//       {
//         *ctx->Ar *= ar;
//       }
//       if (ctx->Ai)
//       {
//         *ctx->Ai *= ar;
//       }
//     }
//     else if (std::abs(ai) > 0.0 && std::abs(ar) == 0.0)
//     {
//       ctx->Ar.swap(ctx->Ai);
//       if (ctx->Ar)
//       {
//         *ctx->Ar *= -ai;
//       }
//       if (ctx->Ai)
//       {
//         *ctx->Ai *= ai;
//       }
//     }
//     else
//     {
//       // General complex coefficient case.
//       mfem::HypreParMatrix *aYr, *aYi;
//       if (ctx->Ar && ctx->Ai)
//       {
//         aYr = mfem::Add(ar, *ctx->Ar, -ai, *ctx->Ai);
//         aYi = mfem::Add(ai, *ctx->Ar, ar, *ctx->Ai);
//         ctx->Ar.reset(aYr);
//         ctx->Ai.reset(aYi);
//       }
//       else if (!ctx->Ar)
//       {
//         ctx->Ar = std::make_unique<mfem::HypreParMatrix>(*ctx->Ai);
//         *ctx->Ar *= -ai;
//         *ctx->Ai *= ar;
//       }
//       else  // !ctx->Ai
//       {
//         ctx->Ai = std::make_unique<mfem::HypreParMatrix>(*ctx->Ar);
//         *ctx->Ar *= ar;
//         *ctx->Ai *= ai;
//       }
//     }
// #else
//     if (ctx->Ar)
//     {
//       *ctx->Ar *= a;
//     }
// #endif
//   }
//   PetscFunctionReturn(0);
// }

// PetscErrorCode __mat_shell_conj(Mat Y)
// {
//   palace::petsc::PetscMatShellCtx *ctx;
//   PetscFunctionBeginUser;

//   PetscCall(MatShellGetContext(Y, (void **)&ctx));
//   MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
// #if defined(PETSC_USE_COMPLEX)
//   if (ctx->Ai)
//   {
//     *ctx->Ai *= -1.0;
//   }
// #endif
//   PetscFunctionReturn(0);
// }

// PetscErrorCode __mat_shell_axpy(Mat Y, PetscScalar a, Mat X, MatStructure str)
// {
//   palace::petsc::PetscMatShellCtx *ctxY, *ctxX;
// #if defined(PETSC_USE_COMPLEX)
//   HYPRE_Real ar, ai;
// #endif
//   PetscFunctionBeginUser;

//   auto Add = [&str](std::unique_ptr<mfem::HypreParMatrix> &Y, HYPRE_Real a,
//                     const std::unique_ptr<mfem::HypreParMatrix> &X)
//   {
//     if (Y)
//     {
//       if (str == SAME_NONZERO_PATTERN)
//       {
//         Y->Add(a, *X);
//       }
//       else
//       {
//         Y.reset(mfem::Add(1.0, *Y, a, *X));
//       }
//     }
//     else
//     {
//       Y = std::unique_ptr<mfem::HypreParMatrix>(*X);
//       *Y *= a;
//     }
//   };
//   PetscCall(MatShellGetContext(Y, (void **)&ctxY));
//   PetscCall(MatShellGetContext(X, (void **)&ctxX));
//   MFEM_VERIFY(ctxY && ctxX, "Invalid PETSc shell matrix contexts!");
// #if defined(PETSC_USE_COMPLEX)
//   ar = PetscRealPart(a);
//   ai = PetscImaginaryPart(a);
//   if (std::abs(ar) > 0.0)
//   {
//     if (ctxX->Ar)
//     {
//       Add(ctxY->Ar, ar, ctxX->Ar);
//     }
//     if (ctxX->Ai)
//     {
//       Add(ctxY->Ai, ar, ctxX->Ai);
//     }
//   }
//   else if (std::abs(ai) > 0.0)
//   {
//     if (ctxX->Ai)
//     {
//       Add(ctxY->Ar, -ai, ctxX->Ai);
//     }
//     if (ctxX->Ar)
//     {
//       Add(ctxY->Ai, ai, ctxX->Ar);
//     }
//   }
// #else
//   if (std::abs(a) > 0.0 && ctxX->Ar)
//   {
//     Add(ctxY->Ar, a, ctxX->Ar);
//   }
// #endif
//   PetscFunctionReturn(0);
// }

// PetscErrorCode __mat_shell_norm(Mat A, NormType type, PetscReal *norm)
// {
//   palace::petsc::PetscMatShellCtx *ctx;
//   PetscFunctionBeginUser;

//   PetscCall(MatShellGetContext(A, (void **)&ctx));
//   MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
//   switch (type)
//   {
//     case NORM_FROBENIUS:
// #if defined(PETSC_USE_COMPLEX)
//       *norm = std::hypot((ctx->Ar) ? hypre_ParCSRMatrixFnorm(*ctx->Ar) : 0.0,
//                          (ctx->Ai) ? hypre_ParCSRMatrixFnorm(*ctx->Ai) : 0.0);
// #else
//       *norm = (ctx->Ar) ? hypre_ParCSRMatrixFnorm(*ctx->Ar) : 0.0;
// #endif
//       break;
//     case NORM_INFINITY:  // Max absolute row sum
// #if defined(PETSC_USE_COMPLEX)
//       if (!ctx->Ar && !ctx->Ai)
//       {
//         *norm = 0.0;
//       }
//       else if (ctx->Ar && !ctx->Ai)
//       {
//         hypre_ParCSRMatrixInfNorm(*ctx->Ar, norm);
//       }
//       else if (ctx->Ai && !ctx->Ar)
//       {
//         hypre_ParCSRMatrixInfNorm(*ctx->Ai, norm);
//       }
//       else
//       {
//         // Need to consider both real and imaginary parts of the matrix.
//         hypre::hypreParCSRInfNorm(*ctx->Ar, *ctx->Ai, norm);
//       }
// #else
//       if (ctx->Ar)
//       {
//         hypre_ParCSRMatrixInfNorm(*ctx->Ar, norm);
//       }
//       else
//       {
//         *norm = 0.0;
//       }
// #endif
//       break;
//     case NORM_1:  // Max absolute column sum (not supported yet)
//     default:
//       MFEM_ABORT("Unsupported matrix norm type!");
//   }
//   PetscFunctionReturn(0);
// }

PetscErrorCode __mat_shell_real_part(Mat Y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(Y, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
#if defined(PETSC_USE_COMPLEX)
  ctx->Ai.reset();
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_shell_imag_part(Mat Y)
{
  palace::petsc::PetscMatShellCtx *ctx;
  PetscFunctionBeginUser;

  PetscCall(MatShellGetContext(Y, (void **)&ctx));
  MFEM_VERIFY(ctx, "Invalid PETSc shell matrix context!");
#if defined(PETSC_USE_COMPLEX)
  ctx->Ar = std::move(ctx->Ai);
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode __mat_convert_hypreParCSR_AIJ(hypre_ParCSRMatrix *hA, Mat *pA)
{
  // Same as PETSc's MatConvert_HYPRE_AIJ function with mtype == MATAIJ, reuse ==
  // MAT_INITIAL_MATRIX, and sameint = true. Slightly modified to allow for using real
  // Hypre matrices (double type) to construct a PETSc matrix with general complex entires
  // (if PETSC_USE_COMPLEX is enabled). See also MFEM's MatConvert_hypreParCSR_AIJ which is
  // a copy of the PETSc version.
  hypre_CSRMatrix *hdiag, *hoffd;
  MPI_Comm comm;
  HYPRE_Int dnnz, onnz, m, n;
  PetscScalar *da, *oa, *aptr;
  PetscInt *dii, *djj, *oii, *ojj, *iptr;
  PetscInt i;
  PetscMPIInt size;
  PetscBool sameint = (PetscBool)(sizeof(PetscInt) == sizeof(HYPRE_Int));
  PetscFunctionBeginUser;

  comm = hypre_ParCSRMatrixComm(hA);
  MPI_Comm_size(comm, &size);
  hdiag = hypre_ParCSRMatrixDiag(hA);
  hoffd = hypre_ParCSRMatrixOffd(hA);
  m = hypre_CSRMatrixNumRows(hdiag);
  n = hypre_CSRMatrixNumCols(hdiag);
  dnnz = hypre_CSRMatrixNumNonzeros(hdiag);
  onnz = hypre_CSRMatrixNumNonzeros(hoffd);
  PetscCall(PetscMalloc1(m + 1, &dii));
  PetscCall(PetscMalloc1(dnnz, &djj));
  PetscCall(PetscMalloc1(dnnz, &da));
  // MFEM_VERIFY(sizeof(HYPRE_Int) == sizeof(PetscInt),
  //             "Index size mismatch inf Hypre-PETSc MatConvert!");
  if (sameint)
  {
    PetscCall(PetscArraycpy(dii, hypre_CSRMatrixI(hdiag), m + 1));
    PetscCall(PetscArraycpy(djj, hypre_CSRMatrixJ(hdiag), dnnz));
  }
  else
  {
    for (i = 0; i < m + 1; i++)
    {
      dii[i] = (PetscInt)(hypre_CSRMatrixI(hdiag)[i]);
    }
    for (i = 0; i < dnnz; i++)
    {
      djj[i] = (PetscInt)(hypre_CSRMatrixJ(hdiag)[i]);
    }
  }
  // This loop replaces the call to PetscArraycpy to convert HYPRE_Complex to PetscScalar
  // values.
  for (i = 0; i < dnnz; i++)
  {
    da[i] = (PetscScalar)(hypre_CSRMatrixData(hdiag)[i]);
  }
  iptr = djj;
  aptr = da;
  for (i = 0; i < m; i++)
  {
    PetscInt nc = dii[i + 1] - dii[i];
    PetscCall(PetscSortIntWithScalarArray(nc, iptr, aptr));
    iptr += nc;
    aptr += nc;
  }
  if (size > 1)
  {
    HYPRE_BigInt *coffd;
    PetscCall(PetscMalloc1(m + 1, &oii));
    PetscCall(PetscMalloc1(onnz, &ojj));
    PetscCall(PetscMalloc1(onnz, &oa));
    if (sameint)
    {
      PetscCall(PetscArraycpy(oii, hypre_CSRMatrixI(hoffd), m + 1));
    }
    else
    {
      for (i = 0; i < m + 1; i++)
      {
        oii[i] = (PetscInt)(hypre_CSRMatrixI(hoffd)[i]);
      }
    }
    coffd = hypre_ParCSRMatrixColMapOffd(hA);
    for (i = 0; i < onnz; i++)
    {
      ojj[i] = (PetscInt)coffd[hypre_CSRMatrixJ(hoffd)[i]];
    }
    for (i = 0; i < onnz; i++)
    {
      oa[i] = (PetscScalar)(hypre_CSRMatrixData(hoffd)[i]);
    }
    iptr = ojj;
    aptr = oa;
    for (i = 0; i < m; i++)
    {
      PetscInt nc = oii[i + 1] - oii[i];
      PetscCall(PetscSortIntWithScalarArray(nc, iptr, aptr));
      iptr += nc;
      aptr += nc;
    }
    PetscCall(MatCreateMPIAIJWithSplitArrays(comm, m, n, PETSC_DECIDE, PETSC_DECIDE, dii,
                                             djj, da, oii, ojj, oa, pA));
  }
  else
  {
    oii = ojj = nullptr;
    oa = nullptr;
    PetscCall(MatCreateSeqAIJWithArrays(comm, m, n, dii, djj, da, pA));
  }
  // We are responsible to free the CSR arrays. However, since we can take references of a
  // PetscParMatrix but we cannot take reference of PETSc arrays, we need to create a
  // PetscContainer object to take reference of these arrays in reference objects.
  void *ptrs[6] = {dii, djj, da, oii, ojj, oa};
  const char *names[6] = {"_csr_dii", "_csr_djj", "_csr_da",
                          "_csr_oii", "_csr_ojj", "_csr_oa"};
  for (i = 0; i < 6; i++)
  {
    PetscContainer c;
    PetscCall(PetscContainerCreate(comm, &c));
    PetscCall(PetscContainerSetPointer(c, ptrs[i]));
    PetscCall(PetscContainerSetUserDestroy(c, __array_container_destroy));
    PetscCall(PetscObjectCompose(reinterpret_cast<PetscObject>(*pA), names[i],
                                 reinterpret_cast<PetscObject>(c)));
    PetscCall(PetscContainerDestroy(&c));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode __array_container_destroy(void *ptr)
{
  PetscFunctionBeginUser;

  PetscCall(PetscFree(ptr));
  PetscFunctionReturn(0);
}
