// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_RAP_HPP
#define PALACE_LINALG_RAP_HPP

#include <memory>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// A parallel operator represented by RAP constructed through the actions of R, A, and P,
// usually with R = Pᵀ, and with possible eliminated essential BC. Here R and P are the
// parallel restriction and prolongation matrices.
//

// Real-valued RAP operator.
class ParOperator : public Operator
{
private:
  // Storage and access for the local operator.
  std::unique_ptr<Operator> data_A;
  const Operator *A;

  // Finite element spaces for parallel prolongation and restriction.
  const FiniteElementSpace &trial_fespace, &test_fespace;
  const bool use_R;

  // Lists of constrained essential boundary true dofs for elimination.
  mfem::Array<int> dbc_tdof_list;

  // Diagonal policy for constrained true dofs.
  DiagonalPolicy diag_policy;

  // Assembled operator as a parallel Hypre matrix. If assembled, the local operator is not
  // deleted.
  mutable std::unique_ptr<mfem::HypreParMatrix> RAP;

  // Helper methods for operator application.
  void RestrictionMatrixMult(const Vector &ly, Vector &ty) const;
  void RestrictionMatrixMultTranspose(const Vector &ty, Vector &ly) const;
  Vector &GetTestLVector() const;

  ParOperator(std::unique_ptr<Operator> &&dA, const Operator *pA,
              const FiniteElementSpace &trial_fespace,
              const FiniteElementSpace &test_fespace, bool test_restrict);

public:
  // Construct the parallel operator, inheriting ownership of the local operator.
  ParOperator(std::unique_ptr<Operator> &&A, const FiniteElementSpace &trial_fespace,
              const FiniteElementSpace &test_fespace, bool test_restrict);
  ParOperator(std::unique_ptr<Operator> &&A, const FiniteElementSpace &fespace)
    : ParOperator(std::move(A), fespace, fespace, false)
  {
  }

  // Non-owning constructors.
  ParOperator(const Operator &A, const FiniteElementSpace &trial_fespace,
              const FiniteElementSpace &test_fespace, bool test_restrict);
  ParOperator(const Operator &A, const FiniteElementSpace &fespace)
    : ParOperator(A, fespace, fespace, false)
  {
  }

  // Get access to the underlying local (L-vector) operator.
  const Operator &LocalOperator() const { return *A; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace.GetComm(); }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list, DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    return dbc_tdof_list.Size() ? &dbc_tdof_list : nullptr;
  }

  // Eliminate essential true dofs from the RHS vector b, using the essential boundary
  // condition values in x.
  void EliminateRHS(const Vector &x, Vector &b) const;

  // Assemble the operator as a parallel sparse matrix. The memory associated with the
  // local operator is free'd.
  mfem::HypreParMatrix &ParallelAssemble(bool skip_zeros = false) const;

  // Steal the assembled parallel sparse matrix.
  std::unique_ptr<mfem::HypreParMatrix> StealParallelAssemble(bool skip_zeros = false) const
  {
    ParallelAssemble(skip_zeros);
    return std::move(RAP);
  }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// Complex-valued RAP operator.
class ComplexParOperator : public ComplexOperator
{
private:
  // Storage and access for the local operator.
  std::unique_ptr<ComplexWrapperOperator> data_A;
  const ComplexWrapperOperator *A;

  // Finite element spaces for parallel prolongation and restriction.
  const FiniteElementSpace &trial_fespace, &test_fespace;
  const bool use_R;

  // Lists of constrained essential boundary true dofs for elimination.
  mfem::Array<int> dbc_tdof_list;

  // Diagonal policy for constrained true dofs.
  Operator::DiagonalPolicy diag_policy;

  // Real and imaginary parts of the operator as non-owning ParOperator objects.
  std::unique_ptr<ParOperator> RAPr, RAPi;

  // Helper methods for operator application.
  void RestrictionMatrixMult(const ComplexVector &ly, ComplexVector &ty) const;
  void RestrictionMatrixMultTranspose(const ComplexVector &ty, ComplexVector &ly) const;
  ComplexVector &GetTestLVector() const;

  ComplexParOperator(std::unique_ptr<Operator> &&dAr, std::unique_ptr<Operator> &&dAi,
                     const Operator *pAr, const Operator *pAi,
                     const FiniteElementSpace &trial_fespace,
                     const FiniteElementSpace &test_fespace, bool test_restrict);

public:
  // Construct the complex-valued parallel operator from the separate real and imaginary
  // parts, inheriting ownership of the local operator.
  ComplexParOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai,
                     const FiniteElementSpace &trial_fespace,
                     const FiniteElementSpace &test_fespace, bool test_restrict);
  ComplexParOperator(std::unique_ptr<Operator> &&Ar, std::unique_ptr<Operator> &&Ai,
                     const FiniteElementSpace &fespace)
    : ComplexParOperator(std::move(Ar), std::move(Ai), fespace, fespace, false)
  {
  }

  // Non-owning constructors.
  ComplexParOperator(const Operator *Ar, const Operator *Ai,
                     const FiniteElementSpace &trial_fespace,
                     const FiniteElementSpace &test_fespace, bool test_restrict);
  ComplexParOperator(const Operator *Ar, const Operator *Ai,
                     const FiniteElementSpace &fespace)
    : ComplexParOperator(Ar, Ai, fespace, fespace, false)
  {
  }

  const Operator *Real() const override { return RAPr.get(); }
  const Operator *Imag() const override { return RAPi.get(); }

  // Get access to the underlying local (L-vector) operator.
  const ComplexOperator &LocalOperator() const { return *A; }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace.GetComm(); }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                            Operator::DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    return dbc_tdof_list.Size() ? &dbc_tdof_list : nullptr;
  }

  void AssembleDiagonal(ComplexVector &diag) const override;

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

  void MultTranspose(const ComplexVector &x, ComplexVector &y) const override;

  void MultHermitianTranspose(const ComplexVector &x, ComplexVector &y) const override;

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;

  void AddMultTranspose(const ComplexVector &x, ComplexVector &y,
                        const std::complex<double> a = 1.0) const override;

  void AddMultHermitianTranspose(const ComplexVector &x, ComplexVector &y,
                                 const std::complex<double> a = 1.0) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_RAP_HPP
