// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "divfree.hpp"

#include <limits>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/multigrid.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

DivFreeSolver::DivFreeSolver(const MaterialOperator &mat_op,
                             mfem::ParFiniteElementSpace &nd_fespace,
                             mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                             const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists,
                             double tol, int max_it, int print, int pa_order_threshold)
{
  constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
  {
    auto M_mg = std::make_unique<MultigridOperator>(h1_fespaces.GetNumLevels());
    for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
    {
      auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
      auto m = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace_l);
      m->AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon_func));
      m->SetAssemblyLevel(
          utils::GetAssemblyLevel(h1_fespace_l.GetMaxElementOrder(), pa_order_threshold));
      m->Assemble(0);
      m->Finalize(0);
      auto M_l = std::make_unique<ParOperator>(
          utils::AssembleOperator(std::move(m), pa_order_threshold), h1_fespace_l);
      M_l->SetEssentialTrueDofs(h1_bdr_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      M_mg->AddOperator(std::move(M_l));
    }
    M = std::move(M_mg);
  }
  {
    // Partial assembly for this operator is only available with libCEED backend.
    auto weakdiv = std::make_unique<mfem::MixedBilinearForm>(
        &nd_fespace, &h1_fespaces.GetFinestFESpace());
    weakdiv->AddDomainIntegrator(
        new mfem::MixedVectorWeakDivergenceIntegrator(epsilon_func));
    weakdiv->SetAssemblyLevel(utils::GetAssemblyLevel());
    weakdiv->Assemble();
    weakdiv->Finalize();
    WeakDiv = std::make_unique<ParOperator>(
        utils::AssembleOperator(std::move(weakdiv), pa_order_threshold), nd_fespace,
        h1_fespaces.GetFinestFESpace(), false);
  }
  {
    // XX TODO: Separate interpolator partial assembly option?
    auto grad = std::make_unique<mfem::DiscreteLinearOperator>(
        &h1_fespaces.GetFinestFESpace(), &nd_fespace);
    grad->AddDomainInterpolator(new mfem::GradientInterpolator);
    grad->SetAssemblyLevel(
        utils::GetAssemblyLevel(nd_fespace.GetMaxElementOrder(), pa_order_threshold));
    grad->Assemble();
    grad->Finalize();
    Grad = std::make_unique<ParOperator>(
        utils::AssembleOperator(std::move(grad), pa_order_threshold),
        h1_fespaces.GetFinestFESpace(), nd_fespace, true);
  }
  bdr_tdof_list_M = &h1_bdr_tdof_lists.back();

  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg =
      std::make_unique<WrapperSolver<Operator>>(std::make_unique<BoomerAmgSolver>(1, 1, 0));
  auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
      std::move(amg), h1_fespaces, nullptr, 1, 1, 2, pa_order_threshold);

  auto pcg =
      std::make_unique<CgSolver<Operator>>(h1_fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
  ksp->SetOperators(*M, *M);

  psi.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  rhs.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
}

}  // namespace palace
