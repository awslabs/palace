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
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

template <typename OperType>
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

}  // namespace

template <typename VecType>
DivFreeSolver<VecType>::DivFreeSolver(
    const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
    FiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists, double tol, int max_it,
    int print)
{
  BlockTimer bt(Timer::DIV_FREE);

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
  if constexpr (std::is_same<VecType, ComplexVector>::value)
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
  if constexpr (std::is_same<VecType, ComplexVector>::value)
  {
    Grad->AddMult(psi.Real(), y.Real(), 1.0);
    Grad->AddMult(psi.Imag(), y.Imag(), 1.0);
  }
  else
  {
    Grad->AddMult(psi, y, 1.0);
  }
}

template class DivFreeSolver<Vector>;
template class DivFreeSolver<ComplexVector>;

}  // namespace palace
