// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hcurl.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

WeightedHCurlNormSolver::WeightedHCurlNormSolver(
    const MaterialOperator &mat_op, FiniteElementSpaceHierarchy &nd_fespaces,
    AuxiliaryFiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
    const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists, double tol, int max_it,
    int print)
{
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == nd_fespaces.GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");
  const auto n_levels = nd_fespaces.GetNumLevels();
  {
    constexpr bool skip_zeros = false;
    MaterialPropertyCoefficient muinv_func(mat_op, mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    MaterialPropertyCoefficient epsilon_func(mat_op, mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    auto A_mg = std::make_unique<MultigridOperator>(n_levels);
    for (bool aux : {false, true})
    {
      for (std::size_t l = 0; l < n_levels; l++)
      {
        // Force coarse level operator to be fully assembled always.
        const auto &fespace_l =
            aux ? h1_fespaces.GetFESpaceAtLevel(l) : nd_fespaces.GetFESpaceAtLevel(l);
        const auto &dbc_tdof_lists_l = aux ? h1_dbc_tdof_lists[l] : nd_dbc_tdof_lists[l];
        BilinearForm a(fespace_l);
        if (aux)
        {
          a.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
        }
        else
        {
          a.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_func, epsilon_func);
        }
        auto A_l = std::make_unique<ParOperator>(a.Assemble(skip_zeros), fespace_l);
        A_l->SetEssentialTrueDofs(dbc_tdof_lists_l, Operator::DiagonalPolicy::DIAG_ONE);
        if (aux)
        {
          A_mg->AddAuxiliaryOperator(std::move(A_l));
        }
        else
        {
          A_mg->AddOperator(std::move(A_l));
        }
      }
    }
    A = std::move(A_mg);
  }

  // The system matrix K + M is real and SPD. We use Hypre's AMS solver as the coarse-level
  // multigrid solve.
  auto ams = std::make_unique<WrapperSolver<Operator>>(std::make_unique<HypreAmsSolver>(
      nd_fespaces.GetFESpaceAtLevel(0), h1_fespaces.GetFESpaceAtLevel(0), 1, 1, 1, false,
      false, 0));
  std::unique_ptr<Solver<Operator>> pc;
  if (n_levels > 1)
  {
    const auto G = h1_fespaces.GetDiscreteInterpolators();
    const int mg_smooth_order =
        std::max(nd_fespaces.GetFinestFESpace().GetMaxElementOrder(), 2);
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        std::move(ams), nd_fespaces.GetProlongationOperators(), &G, 1, 1, mg_smooth_order,
        1.0, 0.0, true);
  }
  else
  {
    pc = std::move(ams);
  }

  auto pcg =
      std::make_unique<CgSolver<Operator>>(nd_fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
  ksp->SetOperators(*A, *A);
}

}  // namespace palace
