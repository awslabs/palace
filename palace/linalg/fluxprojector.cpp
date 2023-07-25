// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fluxprojector.hpp"

#include <limits>
#include "fem/coefficient.hpp"
#include "fem/multigrid.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

// Given a finite element space hierarchy, construct a vector of mass matrix
// operators corresponding to each level.
std::unique_ptr<Operator> BuildMassMatrixOperator(mfem::ParFiniteElementSpaceHierarchy &h, int pa_order_threshold)
{
  constexpr int skip_zeros = 0;
  const bool is_scalar_FE_space =
      h.GetFESpaceAtLevel(0).GetFE(0)->GetRangeType() == mfem::FiniteElement::SCALAR;

  // Assemble the bilinear form operator
  auto M = std::make_unique<MultigridOperator>(h.GetNumLevels());
  for (int l = 0; l < h.GetNumLevels(); ++l)
  {
    auto &h_l = h.GetFESpaceAtLevel(l);
    auto m = std::make_unique<mfem::SymmetricBilinearForm>(&h_l);

    if (is_scalar_FE_space)
    {
      MFEM_ASSERT(h_l.GetVDim() == 1,
                  "Scalar mass matrix hierarchy assumes a component-wise solve.");
      m->AddDomainIntegrator(new mfem::MassIntegrator);
    }
    else
    {
      m->AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    }

    auto M_l = std::make_unique<ParOperator>(
      utils::AssembleOperator(std::move(m), true, (l > 0) ? pa_order_threshold : 100, skip_zeros), h_l);

    // Set the essential dofs (none).
    M->AddOperator(std::move(M_l));
  }
  return M;
}

FluxProjector::FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fes,
                             double tol, int max_it, int print_level, int pa_order_threshold)
  : M(BuildMassMatrixOperator(smooth_flux_fes, pa_order_threshold))
{
  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg = std::make_unique<WrapperSolver<Operator>>(std::make_unique<BoomerAmgSolver>(1, 1, 0));
  auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
      std::move(amg), smooth_flux_fes, nullptr, 1, 1, 2, 1.0, 0.0, true,
      pa_order_threshold);

  auto pcg = std::make_unique<CgSolver<Operator>>(
      smooth_flux_fes.GetFinestFESpace().GetComm(), print_level);

  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
  ksp->SetOperators(*M, *M);

  tmp.SetSize(smooth_flux_fes.GetFinestFESpace().GetTrueVSize());
}

}  // namespace palace
