// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_RAP_HPP
#define PALACE_LINALG_RAP_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// A parallel operators represented by RAP constructed through the actions of R, A, and P,
// usually with R = Pᵀ, and with possible eliminated essential BC.
//

// Real-valued RAP operator.
class ParOperator : public Operator
{
private:
  // Storage and access for the local operator.
  std::unique_ptr<Operator> data_A;
  Operator *A;

  // Finite element spaces for parallel prolongation and restriction.
  const mfem::ParFiniteElementSpace &trial_fespace, &test_fespace;
  const bool use_R;

  // Lists of constrained essential boundary true dofs for elimination.
  const mfem::Array<int> *dbc_tdof_list;

  // Diagonal policy for constrained true dofs.
  DiagonalPolicy diag_policy;

  // Assembled operator as a parallel Hypre matrix. If assembled, the local operator is not
  // deleted.
  std::unique_ptr<mfem::HypreParMatrix> RAP;

  // Temporary storage for operator application.
  mutable Vector lx, ly, ty;

  ParOperator(std::unique_ptr<Operator> &&data_A, Operator *A,
              const mfem::ParFiniteElementSpace &trial_fespace,
              const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);

public:
  // Construct the parallel operator, inheriting ownership of the local operator.
  ParOperator(std::unique_ptr<Operator> &&A,
              const mfem::ParFiniteElementSpace &trial_fespace,
              const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);
  ParOperator(std::unique_ptr<Operator> &&A, const mfem::ParFiniteElementSpace &fespace)
    : ParOperator(std::move(A), fespace, fespace, false)
  {
  }

  // Non-owning constructors.
  ParOperator(Operator *A, const mfem::ParFiniteElementSpace &trial_fespace,
              const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);
  ParOperator(Operator *A, const mfem::ParFiniteElementSpace &fespace)
    : ParOperator(A, fespace, fespace, false)
  {
  }

  // Get access to the underlying local (L-vector) operator.
  const Operator &LocalOperator() const;

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list, DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace.GetComm(); }

  // Assemble the diagonal for the parallel operator.
  void AssembleDiagonal(Vector &diag) const override;

  // Assemble the operator as a parallel sparse matrix. The memory associated with the
  // local operator is not freed.
  mfem::HypreParMatrix &ParallelAssemble();

  // Steal the assembled parallel sparse matrix.
  std::unique_ptr<mfem::HypreParMatrix> StealParallelAssemble()
  {
    ParallelAssemble();
    return std::move(RAP);
  }

  // Eliminate essential true dofs from the RHS vector b, using the essential boundary
  // condition values in x.
  void EliminateRHS(const Vector &x, Vector &b) const;

  void Mult(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    AddMult(x, y);
  }

  void MultTranspose(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    AddMultTranspose(x, y);
  }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// Complex-valued RAP operator.
class ComplexParOperator : public ComplexOperator
{
private:
  // Storage and access for the local operator.
  std::unique_ptr<ComplexOperator> data_A;
  ComplexOperator *A;

  // Finite element spaces for parallel prolongation and restriction.
  const mfem::ParFiniteElementSpace &trial_fespace, &test_fespace;
  const bool use_R;

  // Lists of constrained essential boundary true dofs for elimination.
  mutable const mfem::Array<int> *dbc_tdof_list;

  // Diagonal policy for constrained true dofs.
  DiagonalPolicy diag_policy;

  // Real and imaginary parts of the operator as non-owning ParOperator objects.
  std::unique_ptr<ParOperator> RAPr, RAPi;

  // Temporary storage for operator application.
  mutable ComplexVector lx, ly, ty;

public:
  // Construct the complex-valued parallel operator, inheriting ownership of the local
  // operator.
  ComplexParOperator(std::unique_ptr<ComplexOperator> &&A,
                     const mfem::ParFiniteElementSpace &trial_fespace,
                     const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);
  ComplexParOperator(std::unique_ptr<ComplexOperator> &&A,
                     const mfem::ParFiniteElementSpace &fespace)
    : ComplexParOperator(std::move(Ar), std::move(A), fespace, fespace, false)
  {
  }
  ComplexParOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai,
                     const mfem::ParFiniteElementSpace &trial_fespace,
                     const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);
  ComplexParOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai,
                     const mfem::ParFiniteElementSpace &fespace)
    : ComplexParOperator(std::move(Ar), std::move(Ai), fespace, fespace, false)
  {
  }

  // Get access to the underlying local (L-vector) operator.
  const ComplexOperator &LocalOperator() const;

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list, DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace.GetComm(); }

  bool IsReal() const override { return A->IsReal(); }
  bool IsImag() const override { return A->IsImag(); }
  bool HasReal() const override { return RAPr != nullptr; }
  bool HasImag() const override { return RAPi != nullptr; }
  const Operator *Real() const override { return RAPr.get(); }
  Operator *Real() override { return RAPr.get(); }
  const Operator *Imag() const override { return RAPi.get(); }
  Operator *Imag() override { return RAPi.get(); }

  using ComplexOperator::AddMult;
  using ComplexOperator::AddMultHermitianTranspose;
  using ComplexOperator::AddMultTranspose;
  using ComplexOperator::Mult;
  using ComplexOperator::MultHermitianTranspose;
  using ComplexOperator::MultTranspose;

  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
            bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMult(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                     bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void MultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                              bool zero_real = false, bool zero_imag = false) const override
  {
    yr = 0.0;
    yi = 0.0;
    AddMultHermitianTranspose(xr, xi, yr, yi, 1.0, zero_real, zero_imag);
  }

  void AddMult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
               const std::complex<double> a = 1.0, bool zero_real = false,
               bool zero_imag = false) const override;

  void AddMultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                        const std::complex<double> a = 1.0, bool zero_real = false,
                        bool zero_imag = false) const override;

  void AddMultHermitianTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi,
                                 const std::complex<double> a = 1.0, bool zero_real = false,
                                 bool zero_imag = false) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_RAP_HPP
