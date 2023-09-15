// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hcurl.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

WeightedHCurlNormSolver::WeightedHCurlNormSolver(
    const MaterialOperator &mat_op, const mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
    const mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
    const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists, double tol, int max_it,
    int print, int pa_order_threshold, bool pa_discrete_interp)
{
  constexpr bool skip_zeros = false;
  constexpr auto MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  constexpr auto MatTypeEps = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
  MaterialPropertyCoefficient<MatTypeEps> epsilon_func(mat_op);
  {
    auto A_mg = std::make_unique<MultigridOperator>(nd_fespaces.GetNumLevels());
    for (int s = 0; s < 2; s++)
    {
      const auto &fespaces = (s == 0) ? nd_fespaces : h1_fespaces;
      const auto &dbc_tdof_lists = (s == 0) ? nd_dbc_tdof_lists : h1_dbc_tdof_lists;
      for (int l = 0; l < fespaces.GetNumLevels(); l++)
      {
        // Force coarse level operator to be fully assembled always.
        const auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
        BilinearForm a(fespace_l);
        if (s == 0)
        {
          a.AddDomainIntegrator(
              std::make_unique<CurlCurlMassIntegrator>(muinv_func, epsilon_func));
        }
        else
        {
          a.AddDomainIntegrator(std::make_unique<DiffusionIntegrator>(epsilon_func));
        }
        auto A_l = std::make_unique<ParOperator>(
            a.Assemble((l > 0) ? pa_order_threshold : 99, skip_zeros), fespace_l);
        A_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
        if (s == 0)
        {
          A_mg->AddOperator(std::move(A_l));
        }
        else
        {
          A_mg->AddAuxiliaryOperator(std::move(A_l));
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
  if (nd_fespaces.GetNumLevels() > 1)
  {
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        std::move(ams), nd_fespaces, &h1_fespaces, 1, 1, 2, 1.0, 0.0, true,
        pa_order_threshold, pa_discrete_interp);
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
