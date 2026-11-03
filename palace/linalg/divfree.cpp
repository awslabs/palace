// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "divfree.hpp"

#include <limits>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "models/materialoperator.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

template <OperatorType OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a,
                           const FiniteElementSpace &fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&a,
                                     const FiniteElementSpace &fespace)
{
  return std::make_unique<ParOperator>(std::move(a), fespace);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&a,
                                            const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(a), nullptr, fespace);
}

// Apply an operator, or its transpose, to a vector or to the parts of a complex vector.
void MultParts(const Operator &A, const Vector &x, Vector &y, double a = 0.0)
{
  if (a == 0.0)
  {
    A.Mult(x, y);
  }
  else
  {
    A.AddMult(x, y, a);
  }
}

void MultParts(const Operator &A, const ComplexVector &x, ComplexVector &y, double a = 0.0)
{
  MultParts(A, x.Real(), y.Real(), a);
  MultParts(A, x.Imag(), y.Imag(), a);
}

void MultTransposeParts(const Operator &A, const Vector &x, Vector &y)
{
  A.MultTranspose(x, y);
}

void MultTransposeParts(const Operator &A, const ComplexVector &x, ComplexVector &y)
{
  A.MultTranspose(x.Real(), y.Real());
  A.MultTranspose(x.Imag(), y.Imag());
}

// The operator Gᵀ B G on the potentials, with identity rows and columns for the essential
// true DOFs.
template <OperatorType OperType>
class GradientOperator : public OperType
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;

private:
  const OperType &B;
  const Operator &Grad;
  const mfem::Array<int> *dbc_tdof_list;
  mutable VecType x_free, u, Bu;

public:
  GradientOperator(const OperType &B, const Operator &Grad,
                   const mfem::Array<int> *dbc_tdof_list)
    : OperType(Grad.Width()), B(B), Grad(Grad), dbc_tdof_list(dbc_tdof_list),
      x_free(Grad.Width()), u(Grad.Height()), Bu(Grad.Height())
  {
    x_free.UseDevice(true);
    u.UseDevice(true);
    Bu.UseDevice(true);
  }

  void Mult(const VecType &x, VecType &y) const override
  {
    x_free = x;
    if (dbc_tdof_list)
    {
      linalg::SetSubVector(x_free, *dbc_tdof_list, 0.0);
    }
    MultParts(Grad, x_free, u);
    B.Mult(u, Bu);
    MultTransposeParts(Grad, Bu, y);
    if (dbc_tdof_list)
    {
      linalg::SetSubVector(y, *dbc_tdof_list, x);
    }
  }
};

// A linear solver as a preconditioner.
template <OperatorType OperType>
class KspPreconditioner : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  const BaseKspSolver<OperType> &ksp;

public:
  KspPreconditioner(const BaseKspSolver<OperType> &ksp) : ksp(ksp) {}

  void SetOperator(const OperType &) override {}

  void Mult(const VecType &x, VecType &y) const override { ksp.Mult(x, y); }
};

}  // namespace

template <typename VecType>
DivFreeSolver<VecType>::DivFreeSolver(
    const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
    FiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists, double tol, int max_it,
    int print)
{
  BlockTimer bt(Timer::DIV_FREE);
  comm = nd_fespace.GetComm();

  // If no boundaries on the mesh have been marked, add a single degree of freedom
  // constraint so the system for the projection is not singular. This amounts to enforcing
  // a scalar potential of 0 at a point in space if it is otherwise completely
  // unconstrained.
  const auto *ptr_h1_bdr_tdof_lists = &h1_bdr_tdof_lists;
  {
    MFEM_VERIFY(
        !h1_bdr_tdof_lists.empty(),
        "Unexpected empty list of boundary true dofs for finite element space hierarchy!");
    HYPRE_BigInt coarse_bdr_tdofs = h1_bdr_tdof_lists[0].Size();
    MPI_Comm comm = h1_fespaces.GetFESpaceAtLevel(0).GetComm();
    Mpi::GlobalSum(1, &coarse_bdr_tdofs, comm);
    if (coarse_bdr_tdofs == 0)
    {
      int root = (h1_fespaces.GetFESpaceAtLevel(0).GetTrueVSize() == 0) ? Mpi::Size(comm)
                                                                        : Mpi::Rank(comm);
      Mpi::GlobalMin(1, &root, comm);
      MFEM_VERIFY(root < Mpi::Size(comm),
                  "No root process found for single true dof constraint!");
      if (root == Mpi::Rank(comm))
      {
        aux_tdof_lists.reserve(h1_fespaces.GetNumLevels());
        for (std::size_t l = 0; l < h1_fespaces.GetNumLevels(); l++)
        {
          auto &tdof_list = aux_tdof_lists.emplace_back(1);
          tdof_list[0] = 0;
        }
        ptr_h1_bdr_tdof_lists = &aux_tdof_lists;
      }
    }
  }

  // Create the mass and weak divergence operators for divergence-free projection.
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  {
    constexpr bool skip_zeros = false;
    BilinearForm m(h1_fespaces.GetFinestFESpace());
    m.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
    // m.AssembleQuadratureData();
    auto m_vec = m.Assemble(h1_fespaces, skip_zeros);
    auto M_mg =
        std::make_unique<BaseMultigridOperator<OperType>>(h1_fespaces.GetNumLevels());
    for (std::size_t l = 0; l < h1_fespaces.GetNumLevels(); l++)
    {
      const auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
      auto M_l = BuildLevelParOperator<OperType>(std::move(m_vec[l]), h1_fespace_l);
      M_l->SetEssentialTrueDofs((*ptr_h1_bdr_tdof_lists)[l],
                                Operator::DiagonalPolicy::DIAG_ONE);
      if (l == h1_fespaces.GetNumLevels() - 1)
      {
        bdr_tdof_list_M = M_l->GetEssentialTrueDofs();
      }
      M_mg->AddOperator(std::move(M_l));
    }
    M = std::move(M_mg);
  }
  {
    // Weak divergence operator is always partially assembled.
    BilinearForm weakdiv(nd_fespace, h1_fespaces.GetFinestFESpace());
    weakdiv.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(epsilon_func);
    WeakDiv = std::make_unique<ParOperator>(weakdiv.PartialAssemble(), nd_fespace,
                                            h1_fespaces.GetFinestFESpace(), false);
  }
  Grad = &nd_fespace.GetDiscreteInterpolator(h1_fespaces.GetFinestFESpace());

  // The system matrix for the projection is real and SPD.
  auto amg = std::make_unique<MfemWrapperSolver<OperType>>(
      std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
  amg->SetDropSmallEntries(false);
  std::unique_ptr<Solver<OperType>> pc;
  if (h1_fespaces.GetNumLevels() > 1)
  {
    const int mg_smooth_order =
        std::max(h1_fespaces.GetFinestFESpace().GetMaxElementOrder(), 2);
    pc = std::make_unique<GeometricMultigridSolver<OperType>>(
        h1_fespaces.GetFinestFESpace().GetComm(), std::move(amg),
        h1_fespaces.GetProlongationOperators(), nullptr, 1, 1, mg_smooth_order, 1.0, 0.0,
        true);
  }
  else
  {
    pc = std::move(amg);
  }

  auto pcg =
      std::make_unique<CgSolver<OperType>>(h1_fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(pc));
  ksp->SetOperators(*M, *M);

  psi.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  rhs.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  psi.UseDevice(true);
  rhs.UseDevice(true);
}

template <typename VecType>
void DivFreeSolver<VecType>::Mult(VecType &y) const
{
  BlockTimer bt(Timer::DIV_FREE);

  // Compute the divergence of y.
  if constexpr (std::is_same_v<VecType, ComplexVector>)
  {
    WeakDiv->Mult(y.Real(), rhs.Real());
    WeakDiv->Mult(y.Imag(), rhs.Imag());
  }
  else
  {
    WeakDiv->Mult(y, rhs);
  }

  // Apply essential BC and solve the linear system.
  if (bdr_tdof_list_M)
  {
    linalg::SetSubVector(rhs, *bdr_tdof_list_M, 0.0);
  }
  ksp->Mult(rhs, psi);

  // Compute the irrotational portion of y and subtract.
  if constexpr (std::is_same_v<VecType, ComplexVector>)
  {
    Grad->AddMult(psi.Real(), y.Real(), 1.0);
    Grad->AddMult(psi.Imag(), y.Imag(), 1.0);
  }
  else
  {
    Grad->AddMult(psi, y, 1.0);
  }
}

template <typename VecType>
void DivFreeSolver<VecType>::Mult(VecType &y, double tol) const
{
  const double tol0 = ksp->GetRelTol();
  ksp->SetRelTol(tol);
  Mult(y);
  ksp->SetRelTol(tol0);
}

template <typename VecType>
void DivFreeSolver<VecType>::Mult(const OperType &B, VecType &y) const
{
  BlockTimer bt(Timer::DIV_FREE);

  // Solve Gᵀ B G ψ = Gᵀ B y with GMRES, preconditioned with the solver for Gᵀ M G, and
  // subtract the gradient of ψ.
  VecType By(y.Size());
  By.UseDevice(true);
  B.Mult(y, By);
  MultTransposeParts(*Grad, By, rhs);
  if (bdr_tdof_list_M)
  {
    linalg::SetSubVector(rhs, *bdr_tdof_list_M, 0.0);
  }
  GradientOperator<OperType> A(B, *Grad, bdr_tdof_list_M);
  KspPreconditioner<OperType> pc(*ksp);
  GmresSolver<OperType> gmres(comm, 0);
  gmres.SetOperator(A);
  gmres.SetPreconditioner(pc);
  gmres.SetRelTol(1.0e-12);
  gmres.SetMaxIter(100);
  gmres.SetRestartDim(100);
  psi = 0.0;
  gmres.Mult(rhs, psi);
  if (!gmres.GetConverged())
  {
    Mpi::Warning(comm,
                 "Divergence-free correction did not converge in {:d} iterations "
                 "(relative residual {:.3e})!\n",
                 gmres.GetNumIterations(), gmres.GetFinalRes() / gmres.GetInitialRes());
  }
  MultParts(*Grad, psi, y, -1.0);
}

template class DivFreeSolver<Vector>;
template class DivFreeSolver<ComplexVector>;

}  // namespace palace
