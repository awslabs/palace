// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "divfree.hpp"

#include <limits>
#include "fem/coefficient.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

DivFreeSolver::DivFreeSolver(const MaterialOperator &mat_op,
                             mfem::ParFiniteElementSpace &nd_fespace,
                             mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                             const std::vector<mfem::Array<int>> &h1_bdr_tdof_lists,
                             double tol, int max_it, int print)
  : mfem::Solver(nd_fespace.GetTrueVSize())
{
  constexpr MaterialPropertyType MatType = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
  {
    M.clear();
    M.reserve(h1_fespaces.GetNumLevels());
    for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
    {
      auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
      auto m = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace_l);
      m->AddDomainIntegrator(new mfem::MixedGradGradIntegrator(epsilon_func));
      // XX TODO: Partial assembly option?
      m->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
      m->Assemble(0);
      m->Finalize(0);
      M.push_back(std::make_unique<ParOperator>(std::move(m), h1_fespace_l, h1_fespace_l));
      M.back()->SetEssentialTrueDofs(h1_bdr_tdof_lists[l],
                                     Operator::DiagonalPolicy::DIAG_ONE);
    }
  }
  {
    // XX TODO: Partial assembly option?
    auto weakDiv = std::make_unique<mfem::MixedBilinearForm>(
        &nd_fespace, &h1_fespaces.GetFinestFESpace());
    weakDiv->AddDomainIntegrator(
        new mfem::MixedVectorWeakDivergenceIntegrator(epsilon_func));
    weakDiv->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    weakDiv->Assemble();
    weakDiv->Finalize();
    WeakDiv = std::make_unique<ParOperator>(std::move(weakDiv), nd_fespace,
                                            h1_fespaces.GetFinestFESpace());
  }
  {
    // XX TODO: Partial assembly option?
    auto grad = std::make_unique<mfem::DiscreteLinearOperator>(
        &h1_fespaces.GetFinestFESpace(), &nd_fespace);
    grad->AddDomainInterpolator(new mfem::GradientInterpolator);
    grad->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    grad->Assemble();
    grad->Finalize();
    Grad = std::make_unique<ParOperator>(std::move(grad), h1_fespaces.GetFinestFESpace(),
                                         nd_fespace, true);
  }

  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg = std::make_unique<BoomerAmgSolver>(1, 1, 0);
  auto gmg = std::make_unique<GeometricMultigridSolver>(std::move(amg), h1_fespaces,
                                                        nullptr, 1, 1, 2);

  auto pcg = std::make_unique<mfem::CGSolver>(h1_fespaces.GetFinestFESpace().GetComm());
  pcg->iterative_mode = false;
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  pcg->SetPrintLevel(print);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(amg));
  ksp->SetOperator(*M.back(), M);

  psi.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  rhs.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
}

}  // namespace palace
