// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_RAP_HPP
#define PALACE_LINALG_RAP_HPP

#include <array>
#include <memory>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace detail
{
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
  DiagonalPolicy diag_policy = DiagonalPolicy::DIAG_ZERO;

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

  // Accessor for trial finite element space.
  const FiniteElementSpace &TrialFiniteElementSpace() const { return trial_fespace; }

  // Accessor for test finite element space.
  const FiniteElementSpace &TestFiniteElementSpace() const { return test_fespace; }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list, DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    return dbc_tdof_list.Size() ? &dbc_tdof_list : nullptr;
  }

  // Get the diagonal policy that was most recently used. If there are no essential dofs,
  // and thus no valid policy, will error.
  DiagonalPolicy GetDiagonalPolicy() const;

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
  Operator::DiagonalPolicy diag_policy = Operator::DiagonalPolicy::DIAG_ZERO;
  ;

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

  // Accessor for trial finite element space.
  const FiniteElementSpace &TrialFiniteElementSpace() const { return trial_fespace; }

  // Accessor for test finite element space.
  const FiniteElementSpace &TestFiniteElementSpace() const { return test_fespace; }

  // Set essential boundary condition true dofs for square operators.
  void SetEssentialTrueDofs(const mfem::Array<int> &tdof_list,
                            Operator::DiagonalPolicy policy);

  // Get the essential boundary condition true dofs associated with the operator. May be
  // nullptr.
  const mfem::Array<int> *GetEssentialTrueDofs() const
  {
    return dbc_tdof_list.Size() ? &dbc_tdof_list : nullptr;
  }

  // Get the diagonal policy that was most recently used. If there are no essential dofs,
  // and thus no valid policy, will error.
  Operator::DiagonalPolicy GetDiagonalPolicy() const;

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
auto BuildParSumOperator(const std::array<double, N> &coeff,
                         const std::array<const ParOperator *, N> &ops,
                         bool set_essential = true)
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

    // Check other existant essential dof arrays are references.
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
auto BuildParSumOperator(const std::array<std::complex<double>, N> &coeff,
                         const std::array<const ComplexParOperator *, N> &ops,
                         bool set_essential = true)
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

    // Check other existant essential dof arrays are references.
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

// Dispatch for ParOperators which have been type erased.
template <typename OperType, std::size_t N>
auto BuildParSumOperator(
    const std::array<
        typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                  std::complex<double>, double>::type,
        N> &coeff,
    const std::array<const OperType *, N> &ops)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  std::array<const ParOperType *, N> par_ops;
  std::transform(ops.begin(), ops.end(), par_ops.begin(),
                 [](const OperType *op) { return dynamic_cast<const ParOperType *>(op); });
  return BuildParSumOperator(coeff, std::move(par_ops));
}

// Dispatcher to convert initializer list or C arrays into std::array whilst deducing sizes
// and types.
template <std::size_t N, typename ScalarType, typename OperatorType>
auto BuildParSumOperator(ScalarType (&&coeff_in)[N], const OperatorType *(&&ops_in)[N],
                         bool set_essential = true)
{
  return BuildParSumOperator(to_array<ScalarType>(std::move(coeff_in)),
                             to_array<const OperatorType *>(std::move(ops_in)),
                             set_essential);
}

}  // namespace palace

#endif  // PALACE_LINALG_RAP_HPP
