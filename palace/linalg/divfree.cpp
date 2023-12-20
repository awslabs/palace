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

namespace palace
{

DivFreeSolver::DivFreeSolver(const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                             AuxiliaryFiniteElementSpaceHierarchy &h1_fespaces,
                             const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists,
                             double tol, int max_it, int print)
{
  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  {
    auto M_mg = std::make_unique<MultigridOperator>(h1_fespaces.GetNumLevels());
    for (std::size_t l = 0; l < h1_fespaces.GetNumLevels(); l++)
    {
      // Force coarse level operator to be fully assembled always.
      const auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
      BilinearForm m(h1_fespace_l);
      m.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
      auto M_l = std::make_unique<ParOperator>(m.Assemble(skip_zeros), h1_fespace_l);
      M_l->SetEssentialTrueDofs(h1_bdr_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      M_mg->AddOperator(std::move(M_l));
    }
    M = std::move(M_mg);
  }
  {
    BilinearForm weakdiv(nd_fespace, h1_fespaces.GetFinestFESpace());
    weakdiv.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(epsilon_func);
    WeakDiv = std::make_unique<ParOperator>(weakdiv.Assemble(skip_zeros), nd_fespace,
                                            h1_fespaces.GetFinestFESpace(), false);
  }
  Grad = &h1_fespaces.GetFinestFESpace().GetDiscreteInterpolator();
  bdr_tdof_list_M = &h1_bdr_tdof_lists.back();

  // The system matrix for the projection is real and SPD.
  auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
      std::make_unique<BoomerAmgSolver>(1, 1, 0));
  std::unique_ptr<Solver<Operator>> pc;
  if (h1_fespaces.GetNumLevels() > 1)
  {
    const int mg_smooth_order =
        std::max(h1_fespaces.GetFinestFESpace().GetMaxElementOrder(), 2);
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        h1_fespaces.GetFinestFESpace().GetComm(),

        std::move(amg), h1_fespaces.GetProlongationOperators(), nullptr, 1, 1,
        mg_smooth_order, 1.0, 0.0, true);
  }
  else
  {
    pc = std::move(amg);
  }

  auto pcg =
      std::make_unique<CgSolver<Operator>>(h1_fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
  ksp->SetOperators(*M, *M);

  psi.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  rhs.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
}

}  // namespace palace
