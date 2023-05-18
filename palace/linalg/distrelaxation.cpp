// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "distrelaxation.hpp"

#include <general/forall.hpp>
#include "linalg/chebyshev.hpp"

namespace palace
{

DistRelaxationSmoother::DistRelaxationSmoother(mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace,
                                               int smooth_it, int cheby_smooth_it,
                                               int cheby_order)
  : mfem::Solver(), pc_it(smooth_it), A(nullptr), A_G(nullptr)
{
  // Construct discrete gradient matrix for the auxiliary space.
  {
    // XX TODO: Partial assembly option?
    auto grad = std::make_unique<mfem::DiscreteLinearOperator>(&h1_fespace, &nd_fespace);
    grad->AddDomainInterpolator(new mfem::GradientInterpolator);
    grad->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    grad->Assemble();
    grad->Finalize();
    G = std::make_unique<ParOperator>(std::move(grad), h1_fespace, nd_fespace, true);
  }

  // Initialize smoothers.
  B = std::make_unique<ChebyshevSmoother>(cheby_smooth_it, cheby_order);
  B_G = std::make_unique<ChebyshevSmoother>(cheby_smooth_it, cheby_order);
  B_G->iterative_mode = false;
}

void DistRelaxationSmoother::SetOperator(const Operator &op, const Operator &op_G)
{
  MFEM_VERIFY(op.Height() == G->Height() && op.Width() == G->Height() &&
                  op_G.Height() == G->Width() && op_G.Width() == G->Width(),
              "Invalid operator sizes for DistRelaxationSmoother!");
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  const auto *PtAP_G = dynamic_cast<const ParOperator *>(&op_G);
  MFEM_VERIFY(PtAP && PtAP_G, "ChebyshevSmoother requires a ParOperator operator!");
  A = PtAP;
  A_G = PtAP_G;

  height = A->Height();
  width = A->Width();
  r.SetSize(height);
  x_G.SetSize(A_G->Height());
  y_G.SetSize(A_G->Height());

  // Set up smoothers for A and A_G.
  B->SetOperator(*A);
  B_G->SetOperator(*A_G);
}

void DistRelaxationSmoother::Mult(const Vector &x, Vector &y) const
{
  // Apply smoother.
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + B (x - A y)
    B->iterative_mode = (iterative_mode || it > 0);
    B->Mult(x, y);

    // y = y + G B_G Gᵀ (x - A y)
    A->Mult(y, r);
    subtract(x, r, r);
    G->MultTranspose(r, x_G);
    if (A_G->GetEssentialTrueDofs())
    {
      x_G.SetSubVector(*A_G->GetEssentialTrueDofs(), 0.0);
    }
    B_G->Mult(x_G, y_G);
    G->AddMult(y_G, y, 1.0);
  }
}

void DistRelaxationSmoother::MultTranspose(const Vector &x, Vector &y) const
{
  // Apply transpose.
  B->iterative_mode = true;
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + G B_Gᵀ Gᵀ (x - A y)
    if (iterative_mode || it > 0)
    {
      A->Mult(y, r);
      subtract(x, r, r);
      G->MultTranspose(r, x_G);
    }
    else
    {
      y = 0.0;
    }
    if (A_G->GetEssentialTrueDofs())
    {
      x_G.SetSubVector(*A_G->GetEssentialTrueDofs(), 0.0);
    }
    B_G->MultTranspose(x_G, y_G);
    G->AddMult(y_G, y, 1.0);

    // y = y + Bᵀ (x - A y)
    B->MultTranspose(x, y);
  }
}

}  // namespace palace
