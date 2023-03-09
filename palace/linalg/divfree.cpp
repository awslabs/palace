// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "divfree.hpp"

#include <limits>
#include "fem/materialoperator.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "utils/mfemcoefficients.hpp"

namespace palace
{

DivFreeSolver::DivFreeSolver(const MaterialOperator &mat_op,
                             const mfem::Array<int> &bdr_marker,
                             mfem::ParFiniteElementSpace &nd_fespace,
                             mfem::ParFiniteElementSpaceHierarchy &h1_fespaces, double tol,
                             int max_it, int print)
  : mfem::Solver(nd_fespace.GetTrueVSize())
{
  MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL> epsilon_func(mat_op);
  MFEM_VERIFY(bdr_marker.Size() ==
                  h1_fespaces.GetFinestFESpace().GetParMesh()->bdr_attributes.Max(),
              "Invalid boundary marker for divergence-free solver!");
  M.reserve(h1_fespaces.GetNumLevels());
  for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    auto &h1_fespace_l = h1_fespaces.GetFESpaceAtLevel(l);
    mfem::Array<int> dbc_tdof_list_l;
    h1_fespace_l.GetEssentialTrueDofs(bdr_marker, dbc_tdof_list_l);

    mfem::ParBilinearForm m(&h1_fespace_l);
    m.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(epsilon_func));
    // m.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    m.Assemble();
    m.Finalize();
    mfem::HypreParMatrix *hM = m.ParallelAssemble();
    hM->EliminateBC(dbc_tdof_list_l, mfem::Operator::DiagonalPolicy::DIAG_ONE);
    M.emplace_back(hM);
  }
  {
    mfem::ParMixedBilinearForm weakDiv(&nd_fespace, &h1_fespaces.GetFinestFESpace());
    weakDiv.AddDomainIntegrator(
        new mfem::MixedVectorWeakDivergenceIntegrator(epsilon_func));
    // weakDiv.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    weakDiv.Assemble();
    weakDiv.Finalize();
    WeakDiv.reset(weakDiv.ParallelAssemble());
  }
  {
    mfem::ParDiscreteLinearOperator grad(&h1_fespaces.GetFinestFESpace(), &nd_fespace);
    grad.AddDomainInterpolator(new mfem::GradientInterpolator);
    // grad.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    grad.Assemble();
    grad.Finalize();
    Grad.reset(grad.ParallelAssemble());
  }
  h1_fespaces.GetFinestFESpace().GetEssentialTrueDofs(bdr_marker, h1_bdr_tdof_list);

  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg = std::make_unique<BoomerAmgSolver>();
  amg->SetCoarseRelaxType(8);
  auto gmg = std::make_unique<GeometricMultigridSolver>(std::move(amg), bdr_marker,
                                                        h1_fespaces, nullptr, 1, 1, 2);
  gmg->SetOperator(M);
  pc = std::move(gmg);

  ksp = std::make_unique<mfem::CGSolver>(h1_fespaces.GetFinestFESpace().GetComm());
  ksp->SetRelTol(tol);
  ksp->SetAbsTol(std::numeric_limits<double>::epsilon());
  ksp->SetMaxIter(max_it);
  ksp->SetPrintLevel(print);
  ksp->SetOperator(*M.back());
  ksp->SetPreconditioner(*pc);

  psi.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  rhs.SetSize(h1_fespaces.GetFinestFESpace().GetTrueVSize());
  xr.SetSize(height);
  xi.SetSize(height);
}

}  // namespace palace
