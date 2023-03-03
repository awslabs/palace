// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fluxprojector.hpp"

#include <limits>
#include "fem/materialoperator.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "utils/mfemcoefficients.hpp"

namespace palace
{

// Given a finite element space hierarchy, construct a vector of mass matrix
// operators corresponding to each level.
auto BuildMassMatrixHierarchy(mfem::ParFiniteElementSpaceHierarchy &h)
{
  std::vector<std::unique_ptr<mfem::Operator>> M;

  const bool is_scalar_FE_space =
      h.GetFESpaceAtLevel(0).GetFE(0)->GetRangeType() == mfem::FiniteElement::SCALAR;

  // Assemble the bilinear form operator
  M.reserve(h.GetNumLevels());
  for (int l = 0; l < h.GetNumLevels(); l++)
  {
    auto &h_l = h.GetFESpaceAtLevel(l);

    mfem::ParBilinearForm m(&h_l);

    if (is_scalar_FE_space)
    {
      auto *vmass = new mfem::VectorMassIntegrator;
      vmass->SetVDim(h_l.GetVDim());
      m.AddDomainIntegrator(new mfem::VectorMassIntegrator);
    }
    else
    {
      m.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    }

    m.Assemble();
    m.Finalize();
    M.emplace_back(m.ParallelAssemble());
  }

  return M;
}

FluxProjector::FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fes,
                             double tol, int max_it, int print)
  : mfem::Solver(smooth_flux_fes.GetFinestFESpace().GetTrueVSize())
{
  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg = std::make_unique<BoomerAmgSolver>();
  amg->SetCoarseRelaxType(8);
  auto gmg = std::make_unique<GeometricMultigridSolver>(std::move(amg), mfem::Array<int>(),
                                                        smooth_flux_fes, nullptr, 1, 1, 2);
  gmg->SetOperator(M);
  pc = std::move(gmg);

  ksp = std::make_unique<mfem::CGSolver>(smooth_flux_fes.GetFinestFESpace().GetComm());
  ksp->SetRelTol(tol);
  ksp->SetAbsTol(std::numeric_limits<double>::epsilon());
  ksp->SetMaxIter(max_it);
  ksp->SetPrintLevel(print);
  ksp->SetOperator(*M.back());
  ksp->SetPreconditioner(*pc);
}

}  // namespace palace
