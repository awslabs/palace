// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_OPERATOR_HPP
#define PALACE_LINALG_OPERATOR_HPP

#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class ComplexOperator;

using Operator = mfem::Operator;
using Vector = mfem::Vector;

//
// Derived operator classes extending mfem::Operator from MFEM.
//

// A parallel operator represented by RᵀAP constructed through the actions of Rᵀ, A, and P
// with possible eliminated essential BC.
class ParOperator : public Operator
{
private:
  std::unique_ptr<Operator> A_;
  const mfem::ParFiniteElementSpace &trial_fespace_, &test_fespace_;
  const bool use_R_;

  // Lists of constrained essential boundary true dofs for elimination.
  const mfem::Array<int> *trial_dbc_tdof_list_, *test_dbc_tdof_list_;

  // Diagonal policy for constrained true dofs.
  DiagonalPolicy diag_policy_;

  // Assembled operator as a parallel Hypre matrix. If the save flag is true, calls to
  // ParallelAssemble will not delete the local operator. This is useful for later on calls
  // to EliminateRHS, for example.
  std::unique_ptr<mfem::HypreParMatrix> RAP_;
  bool save_A_;

  // Temporary storage for operator application.
  mutable Vector lx_, ly_, tx_, ty_;

public:
  // Construct the parallel operator, inheriting ownership of the local operator.
  ParOperator(std::unique_ptr<Operator> &&A,
              const mfem::ParFiniteElementSpace &trial_fespace,
              const mfem::ParFiniteElementSpace &test_fespace, bool test_restrict);
  ParOperator(std::unique_ptr<Operator> &&A, const mfem::ParFiniteElementSpace &fespace)
    : ParOperator(std::move(A), fespace, fespace, false)
  {
  }

  // Get access to the underlying local (L-vector) operator.
  const Operator &LocalOperator() const
  {
    MFEM_ASSERT(A_, "No local matrix available for ParOperator::LocalOperator!");
    return *A_;
  }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &dbc_tdof_list,
                            DiagonalPolicy diag_policy)
  {
    MFEM_VERIFY(height == width, "Set essential true dofs for both test and trial spaces "
                                 "for rectangular ParOperator!");
    trial_dbc_tdof_list_ = &dbc_tdof_list;
    test_dbc_tdof_list_ = &dbc_tdof_list;
    diag_policy_ = diag_policy;
  }

  // Set essential boundary condition true dofs for rectangular operators.
  void SetEssentialTrueDofs(const mfem::Array<int> *trial_dbc_tdof_list,
                            const mfem::Array<int> *test_dbc_tdof_list,
                            DiagonalPolicy diag_policy)
  {
    MFEM_VERIFY(diag_policy == DiagonalPolicy::DIAG_ZERO,
                "Essential boundary condition true dof elimination for rectangular "
                "ParOperator only supports DiagonalPolicy::DIAG_ZERO!");
    trial_dbc_tdof_list_ = trial_dbc_tdof_list;
    test_dbc_tdof_list_ = test_dbc_tdof_list;
    diag_policy_ = diag_policy;
  }

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    MFEM_VERIFY(trial_dbc_tdof_list_ == test_dbc_tdof_list_ && height == width,
                "GetEssentialTrueDofs should only be used for square ParOperator!");
    return trial_dbc_tdof_list_;
  }

  // Get access to the finite element spaces associated with the operator.
  const mfem::ParFiniteElementSpace &GetFESpace() const
  {
    MFEM_VERIFY(&trial_fespace_ == &test_fespace_ && height == width,
                "GetFESpace should only be used for square ParOperator!");
    return trial_fespace_;
  }

  // A call to ParallelAssemble will typically free the memory associated with the local
  // operator as it is no longer required. When the save flag is set, the local operator
  // will not be deleted during parallel assembly.
  void SaveLocalOperator() { save_A_ = true; }

  // Eliminate essential true dofs from the RHS vector b, using the essential boundary
  // condition values in x.
  void EliminateRHS(const Vector &x, Vector &b) const;

  // Assemble the diagonal for the parallel operator.
  void AssembleDiagonal(Vector &diag) const override;

  // Assemble the operator as a parallel sparse matrix. This frees the memory associated
  // with the local operator.
  mfem::HypreParMatrix &ParallelAssemble();

  // Steal the assembled parallel sparse matrix. The local operator is saved so that this
  // object still can perform operations after this is called.
  std::unique_ptr<mfem::HypreParMatrix> StealParallelAssemble()
  {
    SaveLocalOperator();
    ParallelAssemble();
    return std::move(RAP_);
  }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return trial_fespace_.GetComm(); }

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// Wrap a sequence of operators of the same dimensions and optional coefficients.
class SumOperator : public Operator
{
private:
  std::vector<std::pair<const Operator *, double>> ops_;

public:
  SumOperator(int s) : Operator(s) {}
  SumOperator(int h, int w) : Operator(h, w) {}
  SumOperator(const Operator &op, double c = 1.0) : Operator(op.Height(), op.Width())
  {
    AddOperator(op, c);
  }

  void AddOperator(const Operator &op, double c = 1.0)
  {
    MFEM_VERIFY(op.Height() == height && op.Width() == width,
                "Invalid Operator dimensions for SumOperator!");
    ops_.emplace_back(&op, c);
  }

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

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    for (const auto &[op, c] : ops_)
    {
      op->AddMult(x, y, a * c);
    }
  }

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    for (const auto &[op, c] : ops_)
    {
      op->AddMultTranspose(x, y, a * c);
    }
  }
};

// Wraps two symmetric operators such that: (AB)ᵀ = BᵀAᵀ = BA.
class SymmetricProductOperator : public Operator
{
private:
  const Operator &A_, &B_;
  mutable Vector z_;

public:
  SymmetricProductOperator(const Operator &A, const Operator &B)
    : Operator(A.Height(), B.Width()), A_(A), B_(B), z_(B_.Height())
  {
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    B_.Mult(x, z_);
    A_.Mult(z_, y);
  }

  void MultTranspose(const Vector &x, Vector &y) const override
  {
    A_.Mult(x, z_);
    B_.Mult(z_, y);
  }
};

// Applies the simple (symmetric) operator: diag(d).
class DiagonalOperator : public Operator
{
private:
  const Vector &d_;

public:
  DiagonalOperator(const Vector &d) : Operator(d.Size()), d_(d) {}

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

namespace linalg
{

// Estimate operator 2-norm (spectral norm) using power iteration. Assumes the operator is
// not symmetric or Hermitian unless specified.
double SpectralNorm(MPI_Comm comm, const Operator &A, bool sym = false, double tol = 1.0e-4,
                    int max_it = 100);
double SpectralNorm(MPI_Comm comm, const ComplexOperator &A, bool herm = false,
                    double tol = 1.0e-4, int max_it = 100);

}  // namespace linalg

}  // namespace palace

namespace mfem
{

// A symmetric bilinear form operator which replaces *MultTranspose with *Mult.
class SymmetricBilinearForm : public BilinearForm
{
public:
  using BilinearForm::BilinearForm;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
  void AddMultTranspose(const Vector &x, Vector &y, double c = 1.0) const override
  {
    AddMult(x, y, c);
  }
};

}  // namespace mfem

#endif  // PALACE_LINALG_OPERATOR_HPP
