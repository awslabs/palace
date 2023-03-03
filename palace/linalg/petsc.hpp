// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_PETSC_HPP
#define PALACE_PETSC_HPP

#include <petscsys.h>
#include <mfem.hpp>

#include "complexvector.hpp"

#if defined(PETSC_HAVE_HYPRE)
#error \
    "PETSc should be built without Hypre to avoid conflicts with MFEM's Hypre dependency!"
#endif
#if !defined(PETSC_USE_REAL_DOUBLE)
#error "PETSc should be compiled with double precision!"
#endif
#if defined(PETSC_USE_64BIT_INDICES) && !(defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT))
#warning "Mismatch between big HYPRE (32bit) and PETSc (64bit) integer types!"
#endif
#if !defined(PETSC_USE_64BIT_INDICES) && (defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT))
#warning "Mismatch between big HYPRE (64bit) and PETSc (32bit) integer types!"
#endif

#include <functional>
#include <memory>

// Forward declarations of PETSc objects.
typedef struct _p_PetscSF *VecScatter;
typedef struct _p_Vec *Vec;
typedef struct _p_Mat *Mat;
typedef struct _p_KSP *KSP;
typedef struct _p_PC *PC;

// Error handling similar to Petsc's PetscCallAbort but always aborts on the global
// PETSC_COMM_WORLD communicator.
#define PalacePetscCall(...) PetscCallAbort(PETSC_COMM_WORLD, __VA_ARGS__)

namespace palace::petsc
{

//
// A minimal implementation of MFEM's PETSc wrappers to support PETSc built with complex
// numbers.
//

class PetscParMatrix;
class PetscParVector;

// Wrappers for PetscInitialize/PetscFinalize.
void Initialize(int &argc, char **&argv, const char rc_file[], const char help[]);
void Finalize();

// Wrapper for PETSc's vector scatter class.
class PetscScatter
{
public:
  enum class Type
  {
    TO_ZERO,
    TO_ALL
  };

private:
  // The actual PETSc object.
  VecScatter ctx;

public:
  // Creates a scatter context that copies all entries from the parallel vector to either
  // all processes or to the root process. Allocates the
  PetscScatter(Type type, const PetscParVector &x, std::unique_ptr<PetscParVector> &y);

  // Calls PETSc's destroy function.
  ~PetscScatter();

  // Routines for forward/reverse scattering.
  void Forward(const PetscParVector &x, PetscParVector &y);
  void Reverse(const PetscParVector &x, PetscParVector &y);
};

// Wrapper for PETSc's vector class.
class PetscParVector
{
private:
  // The actual PETSc object.
  Vec x;

public:
  // Creates vector compatible with (i.e. in the domain of) A or Aᵀ.
  PetscParVector(const PetscParMatrix &A, bool transpose = false);

  // Parallel and serial copy constructors from MFEM's Vector object.
  PetscParVector(MPI_Comm comm, const mfem::Vector &y);
  PetscParVector(const mfem::Vector &y);
#if defined(PETSC_USE_COMPLEX)
  PetscParVector(MPI_Comm comm, const mfem::Vector &yr, const mfem::Vector &yi);
  PetscParVector(const mfem::Vector &yr, const mfem::Vector &yi);
  PetscParVector(MPI_Comm comm, const ComplexVector &cv)
  : PetscParVector(comm, cv.real, cv.imag) {};
  PetscParVector(const ComplexVector &cv)
  : PetscParVector(cv.real, cv.imag) {};
#endif

  // Create a parallel or sequential PETSc vector with the provided dimension.
  PetscParVector(MPI_Comm comm, PetscInt n, PetscInt N);
  // PetscParVector(PetscInt n);

  // Create a parallel or sequential PETSc vector with a data array.
  PetscParVector(MPI_Comm comm, PetscInt n, PetscInt N, PetscScalar *data);
  PetscParVector(PetscInt n, PetscScalar *data);

  // Copy constructor, calls VecDuplicate.
  PetscParVector(const PetscParVector &y);

  // Constructor which wraps an existing PETSc Vec object and takes over ownership unless
  // ref is true.
  PetscParVector(Vec y, bool ref);

  // Calls PETSc's destroy function.
  virtual ~PetscParVector();

  // Copy to/from MFEM's Vector type.
  void GetToVector(mfem::Vector &v, PetscInt start = -1, PetscInt end = -1) const;
  void SetFromVector(const mfem::Vector &v);
  void AddFromVector(const mfem::Vector &v);
#if defined(PETSC_USE_COMPLEX)
  void GetToVectors(mfem::Vector &vr, mfem::Vector &vi, PetscInt start = -1,
                    PetscInt end = -1) const;
  ComplexVector GetToVectors(PetscInt start = -1, PetscInt end = -1) const;
  void SetFromVectors(const mfem::Vector &vr, const mfem::Vector &vi);
  void AddFromVectors(const mfem::Vector &vr, const mfem::Vector &vi);
  void SetFromVectors(const ComplexVector &v) { SetFromVectors(v.real, v.imag); }
  void AddFromVectors(const ComplexVector &v) { AddFromVectors(v.real, v.imag); }
#endif

  // Access the data array of the vector.
  PetscScalar *GetArray();
  const PetscScalar *GetArrayRead() const;
  void RestoreArray(PetscScalar *data);
  void RestoreArrayRead(const PetscScalar *data) const;

  // Temporarily replace the data array of the vector.
  void PlaceArray(const PetscScalar *data);
  void ResetArray();

  // Copy entries of y to x.
  void Copy(const PetscParVector &y);

  // Returns the local vector size.
  PetscInt GetSize() const;

  // Returns the global vector size.
  PetscInt GetGlobalSize() const;

  // Set the (local) vector dimension to n, copying previous contents to the upper block.
  void Resize(PetscInt n, bool copy = false);

  // Zero all entries of the vector.
  void SetZero();

  // Sets all entries of the vector to random numbers sampled from the range[-1-i, 1+i], or
  // [-1, 1].
  void SetRandom();
#if defined(PETSC_USE_COMPLEX)
  void SetRandomReal();
#else
  void SetRandomReal() { SetRandom(); }
#endif
  void SetRandomSign(bool init = false);

  // Set all entries to s.
  PetscParVector &operator=(PetscScalar s);

  // Scale all entries by s.
  void Scale(PetscScalar s);

  // Shift all entries by +s.
  void Shift(PetscScalar s);

  // Compute pointwise |x|.
  void Abs();

  // Compute pointwise sqrt(|x|).
  void SqrtAbs();

  // Compute pointwise 1/x.
  void Inv();

  // Compute pointwise 1/sqrt(x).
  void InvSqrt();

#if defined(PETSC_USE_COMPLEX)
  // Replace entries with complex conjugate.
  void Conj();

  // Zero the imaginary part of the vector.
  void GetRealPart();

  // Move the imaginary part to the real part of the vector.
  void GetImagPart();
#endif

  // Normalize the vector.
  PetscReal Normalize();
  PetscReal Normalize(const PetscParMatrix &B, PetscParVector &Bv);

  // Calculate the vector 2-norm.
  PetscReal Norml2() const;

  // Calculate the vector infinity-norm.
  PetscReal Normlinf() const;

  // Zero specified (local) rows of the vector.
  void ZeroRows(const mfem::Array<int> &rows);

  // Pointwise multiplication x *= y.
  void PointwiseMult(const PetscParVector &y, bool replace_zeros);

  // In-place addition x += alpha * y.
  void AXPY(PetscScalar alpha, const PetscParVector &y);

  // In-place addition x = alpha * y + beta * x.
  void AXPBY(PetscScalar alpha, const PetscParVector &y, PetscScalar beta);

  // In-place addition x = alpha * y + beta * z + gamma * x.
  void AXPBYPCZ(PetscScalar alpha, const PetscParVector &y, PetscScalar beta,
                const PetscParVector &z, PetscScalar gamma);

  // Vector dot product (yᴴ x) or indefinite dot product (yᵀ x) for complex vectors.
  PetscScalar Dot(const PetscParVector &y) const;
  PetscScalar TransposeDot(const PetscParVector &y) const;

  // Prints the vector (to stdout if fname is nullptr).
  void Print(const char *fname = nullptr, bool binary = false) const;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const;

  // Typecasting to PETSc's Vec type.
  operator Vec() const { return x; }

  // Typecasting to PETSc object.
  operator PetscObject() const { return reinterpret_cast<PetscObject>(x); }
};

// Base wrapper for PETSc's matrix class.
class PetscParMatrix
{
public:
  enum class NNZStructure
  {
    DIFFERENT,
    SAME,
    SUBSET
  };

#if defined(PETSC_USE_COMPLEX)
  enum class ExtractStructure
  {
    REAL,
    IMAGINARY,
    SUM
  };
#endif

protected:
  // The actual PETSc object.
  Mat A;

  // Default constructor for derived classes.
  PetscParMatrix() : A(nullptr) {}

public:
  // Copy constructor, calls MatDuplicate.
  PetscParMatrix(const PetscParMatrix &B);

  // Constructor which wraps an existing PETSc Mat object and takes over ownership unless
  // ref is true.
  PetscParMatrix(Mat B, bool ref);

  // Calls PETSc's destroy function.
  virtual ~PetscParMatrix();

  // Get/set symmetric or Hermitian flags for the matrix. When setting the flags, it is
  // assumed the structure does not change for the lifetime of the matrix(unless explicitly
  // set again).
  void SetSymmetric(bool sym = true);
  void SetHermitian(bool herm = true);
  bool GetSymmetric() const;
  bool GetHermitian() const;
#if defined(PETSC_USE_COMPLEX)
  void SetRealSymmetric();
#endif
  void CopySymmetry(const PetscParMatrix &B);

  // Returns the local number of rows.
  PetscInt GetNumRows() const;
  PetscInt Height() const { return GetNumRows(); }

  // Returns the local number of columns.
  PetscInt GetNumCols() const;
  PetscInt Width() const { return GetNumCols(); }

  // Returns the global number of rows.
  PetscInt GetGlobalNumRows() const;

  // Returns the global number of columns.
  PetscInt GetGlobalNumCols() const;

  // Returns the number of nonzeros.
  virtual PetscInt NNZ() const;
#if defined(PETSC_USE_COMPLEX)
  virtual PetscInt NNZReal() const
  {
    MFEM_ABORT("NNZReal is not supported for base class PetscParMatrix!");
    return 0;
  }
  virtual PetscInt NNZImag() const
  {
    MFEM_ABORT("NNZImag is not supported for base class PetscParMatrix!");
    return 0;
  }
#endif

  // Calculate matrix Frobenius and infinity norms.
  PetscReal NormF() const;
  PetscReal NormInf() const;
#if defined(PETSC_USE_COMPLEX)
  virtual PetscReal NormFReal() const
  {
    MFEM_ABORT("NormFReal is not supported for base class PetscParMatrix!");
    return 0.0;
  }
  virtual PetscReal NormFImag() const
  {
    MFEM_ABORT("NormFImag is not supported for base class PetscParMatrix!");
    return 0.0;
  }
  virtual PetscReal NormInfReal() const
  {
    MFEM_ABORT("NormInfReal is not supported for base class PetscParMatrix!");
    return 0.0;
  }
  virtual PetscReal NormInfImag() const
  {
    MFEM_ABORT("NormInfImag is not supported for base class PetscParMatrix!");
    return 0.0;
  }
#endif

  // Estimate matrix 2-norm (spectral norm) using power iteration.
  PetscReal Norm2(PetscReal tol = PETSC_DEFAULT, PetscInt maxits = PETSC_DEFAULT) const;

  // Scale all entries by s.
  void Scale(PetscScalar s);

#if defined(PETSC_USE_COMPLEX)
  // Replace entries with complex conjugate.
  void Conj();

  // Zero the imaginary part of the matrix.
  void GetRealPart();

  // Move the imaginary part to the real part of the matrix.
  void GetImagPart();
#endif

  // In-place addition A += alpha * B.
  void AXPY(PetscScalar alpha, const PetscParMatrix &B, NNZStructure struc);

  // Matrix-vector multiplication.
  void Mult(const PetscParVector &x, PetscParVector &y) const;
  void MultAdd(const PetscParVector &x, PetscParVector &y) const;
  void MultTranspose(const PetscParVector &x, PetscParVector &y) const;
  void MultTransposeAdd(const PetscParVector &x, PetscParVector &y) const;
  void MultHermitianTranspose(const PetscParVector &x, PetscParVector &y) const;
  void MultHermitianTransposeAdd(const PetscParVector &x, PetscParVector &y) const;

#if defined(PETSC_USE_COMPLEX)
  // Multiplication with a real-valued vector.
  virtual void Mult(const mfem::Vector &x, PetscParVector &y) const;
  virtual void MultTranspose(const mfem::Vector &x, PetscParVector &y) const;
  virtual void MultHermitianTranspose(const mfem::Vector &x, PetscParVector &y) const;
#endif

  // Prints the matrix (to stdout if fname is nullptr).
  virtual void Print(const char *fname = nullptr, bool binary = false) const;
#if defined(PETSC_USE_COMPLEX)
  virtual void PrintReal(const char *fname) const
  {
    MFEM_ABORT("PrintReal is not supported for base class PetscParMatrix!");
  }
  virtual void PrintImag(const char *fname) const
  {
    MFEM_ABORT("PrintImag is not supported for base class PetscParMatrix!");
  }
#endif

  // Returns a (real) MFEM Operator from the underlying shell matrix data. When complex
  // scalars are used, the parameter controls which part of the matrix to extract.
#if defined(PETSC_USE_COMPLEX)
  virtual const mfem::Operator *GetOperator(ExtractStructure struc) const
#else
  virtual const mfem::Operator *GetOperator() const
#endif
  {
    MFEM_ABORT("GetOperator is not supported for base class PetscParMatrix!");
    return nullptr;
  }

  // Test whether or not a shell matrix has a real or imaginary parts.
#if defined(PETSC_USE_COMPLEX)
  virtual bool HasReal() const
  {
    MFEM_ABORT("HasReal is not supported for base class PetscParMatrix!");
    return false;
  }
  virtual bool HasImag() const
  {
    MFEM_ABORT("HasImag is not supported for base class PetscParMatrix!");
    return false;
  }
#endif

  // Constructs a (real) HypreParMatrix from the PETSc matrix data. When complex scalars
  // are used, the parameter controls which part of the matrix to extract.
#if defined(PETSC_USE_COMPLEX)
  virtual std::unique_ptr<mfem::HypreParMatrix>
  GetHypreParMatrix(ExtractStructure struc) const;
#else
  virtual std::unique_ptr<mfem::HypreParMatrix> GetHypreParMatrix() const;
#endif

  // Create a submatrix on the same number of processors as the original matrix,
  // corresponding to the provided rows and columns which are the selected(local) indices.
  virtual std::unique_ptr<PetscParMatrix> GetSubMatrix(const mfem::Array<int> &rows,
                                                       const mfem::Array<int> &cols);

  // Create a sequential gathered matrix corresponding to the parallel matrix. All processes
  // on the original communicator must call this function, but if the argument is false, no
  // matrix is created (returned pointer is nullptr).
  virtual std::unique_ptr<PetscParMatrix> GetSequentialMatrix(bool create);

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const;

  // Typecasting to PETSc's Mat type.
  operator Mat() const { return A; }

  // Typecasting to PETSc object.
  operator PetscObject() const { return reinterpret_cast<PetscObject>(A); }
};

// Context data for PETSc shell matrices. These store complex matrices as
// Ar + i Ai and perform matrix-vector products.
struct PetscMatShellCtx
{
  std::unique_ptr<mfem::Operator> Ar;
  mfem::Vector x, y;
#if defined(PETSC_USE_COMPLEX)
  std::unique_ptr<mfem::Operator> Ai;
#endif
};

// Wrapper for PETSc's MATSHELL matrix class.
class PetscShellMatrix : public PetscParMatrix
{
private:
  // Returns the shell matrix context.
  PetscMatShellCtx *GetContext() const;

public:
  // Create a PETSc shell matrix wrapping an MFEM Operator. Ownership of the operator is
  // transfered to the PETSc shell. When PETSc is compiled with complex numbers support,
  // the shell matrix wraps the real and imaginary parts to act on complex PETSc Vec
  // objects.
  PetscShellMatrix(MPI_Comm comm, std::unique_ptr<mfem::Operator> &&B);
#if defined(PETSC_USE_COMPLEX)
  PetscShellMatrix(MPI_Comm comm, std::unique_ptr<mfem::Operator> &&Br,
                   std::unique_ptr<mfem::Operator> &&Bi);
#endif

  // Returns the number of nonzeros.
  PetscInt NNZ() const override;
#if defined(PETSC_USE_COMPLEX)
  PetscInt NNZReal() const override;
  PetscInt NNZImag() const override;
#endif

  // Calculate matrix Frobenius and infinity norms.
#if defined(PETSC_USE_COMPLEX)
  PetscReal NormFReal() const override;
  PetscReal NormFImag() const override;
  PetscReal NormInfReal() const override;
  PetscReal NormInfImag() const override;
#endif

#if defined(PETSC_USE_COMPLEX)
  // Multiplication with a real-valued vector.
  void Mult(const mfem::Vector &x, PetscParVector &y) const override;
  void MultTranspose(const mfem::Vector &x, PetscParVector &y) const override;
  void MultHermitianTranspose(const mfem::Vector &x, PetscParVector &y) const override;
#endif

  // Prints the locally owned matrix rows in parallel.
  void Print(const char *fname = nullptr, bool binary = false) const override;
#if defined(PETSC_USE_COMPLEX)
  void PrintReal(const char *fname) const override;
  void PrintImag(const char *fname) const override;
#endif

  // Test whether or not a shell matrix has a real or imaginary parts.
#if defined(PETSC_USE_COMPLEX)
  bool HasReal() const override;
  bool HasImag() const override;
#endif

  // Returns a (real) MFEM Operator from the underlying shell matrix data. When complex
  // scalars are used, the parameter controls which part of the matrix to extract.
#if defined(PETSC_USE_COMPLEX)
  const mfem::Operator *GetOperator(ExtractStructure struc) const override;
#else
  const mfem::Operator *GetOperator() const override;
#endif

  // These methods are not supported for MATSHELL.
#if defined(PETSC_USE_COMPLEX)
  std::unique_ptr<mfem::HypreParMatrix>
  GetHypreParMatrix(ExtractStructure struc) const override
#else
  std::unique_ptr<mfem::HypreParMatrix> GetHypreParMatrix() const override
#endif
  {
    MFEM_ABORT("GetHypreParMatrix is not supported for PetscShellMatrix!");
    return {};
  }
  std::unique_ptr<PetscParMatrix> GetSubMatrix(const mfem::Array<int> &,
                                               const mfem::Array<int> &) override
  {
    MFEM_ABORT("GetSubMatrix is not supported for PetscShellMatrix!");
    return {};
  }
  std::unique_ptr<PetscParMatrix> GetSequentialMatrix(bool) override
  {
    MFEM_ABORT("GetSequentialMatrix is not supported for PetscShellMatrix!");
    return {};
  }
};

// Wrapper for PETSc's MATIJ matrix class.
class PetscAijMatrix : public PetscParMatrix
{
public:
  // Create a PETSc matrix explicitly converted from an MFEM Operator.
  PetscAijMatrix(const mfem::Operator &B);
#if defined(PETSC_USE_COMPLEX)
  PetscAijMatrix(const mfem::Operator &Br, const mfem::Operator &Bi);
#endif
};

// Wrapper for PETSc's MATDENSE matrix class.
class PetscDenseMatrix : public PetscParMatrix
{
private:
  // Helper method for column orthonormalization.
  PetscReal OrthonormalizeColumnInternal(
      PetscInt j, bool mgs, bool cgs2,
      const std::function<PetscScalar(PetscParVector &, PetscParVector &)> &Dot,
      const std::function<void(PetscParVector &, PetscDenseMatrix &, PetscParVector &)>
          &VecDot,
      const std::function<PetscReal(PetscParVector &)> &Normalize);

public:
  // Create a parallel or sequential PETSc dense matrix. Option to specify an existing data
  // array.
  PetscDenseMatrix(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt M, PetscInt N,
                   PetscScalar *data);
  PetscDenseMatrix(PetscInt m, PetscInt n, PetscScalar *data);

  // Set the (local) matrix dimensions to m x n, copying previous contents to the upper-left
  // block.
  void Resize(PetscInt m, PetscInt n, bool copy = false);

  // Access methods for columns of the dense matrix.
  PetscParVector GetColumn(PetscInt j);
  const PetscParVector GetColumnRead(PetscInt j) const;
  void RestoreColumn(PetscInt j, PetscParVector &v);
  void RestoreColumnRead(PetscInt j, const PetscParVector &v) const;

  // Access the data array of the dense matrix.
  PetscScalar *GetArray();
  const PetscScalar *GetArrayRead() const;
  void RestoreArray(PetscScalar *data);
  void RestoreArrayRead(const PetscScalar *data) const;

  // Sets all entries of the vector to random numbers sampled from the range[-1-i, 1+i], or
  // [-1, 1].
  void SetRandom(PetscInt start = -1, PetscInt end = -1);
#if defined(PETSC_USE_COMPLEX)
  void SetRandomReal(PetscInt start = -1, PetscInt end = -1);
#else
  void SetRandomReal(PetscInt start = -1, PetscInt end = -1) { SetRandom(start, end); }
#endif
  void SetRandomSign(PetscInt start = -1, PetscInt end = -1, bool init = false);

  // Orthonormalize column j of the matrix against the preceeding columns, using classical
  // or modified Gram-Schmidt.
  PetscReal OrthonormalizeColumn(PetscInt j, bool mgs, bool cgs2);
  PetscReal OrthonormalizeColumn(PetscInt j, bool mgs, bool cgs2, const PetscParMatrix &B,
                                 PetscParVector &Bv);

  // Dense matrix-matrix multiplication.
  void MatMult(const PetscDenseMatrix &X, PetscDenseMatrix &Y) const;
  void MatMultTranspose(const PetscDenseMatrix &X, PetscDenseMatrix &Y) const;
  void MatTransposeMult(const PetscDenseMatrix &X, PetscDenseMatrix &Y) const;
};

}  // namespace palace::petsc

#endif  // PALACE_PETSC_HPP
