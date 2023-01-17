// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "distrelaxation.hpp"

#include "linalg/chebyshev.hpp"

namespace palace
{

DistRelaxationSmoother::DistRelaxationSmoother(mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace,
                                               const mfem::Array<int> &dbc_marker,
                                               int smooth_it, int cheby_smooth_it,
                                               int cheby_order)
  : mfem::Solver(), A(nullptr), A_G(nullptr), pc_it(smooth_it)
{
  // Construct discrete gradient matrix for the auxiliary space.
  {
    mfem::ParDiscreteLinearOperator grad(&h1_fespace, &nd_fespace);
    grad.AddDomainInterpolator(new mfem::GradientInterpolator);
    // grad.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    grad.Assemble();
    grad.Finalize();
    G.reset(grad.ParallelAssemble());
  }

  // Initialize smoothers.
  mfem::Array<int> nd_dbc_tdof_list;
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  B = std::make_unique<ChebyshevSmoother>(nd_fespace.GetComm(), nd_dbc_tdof_list,
                                          cheby_smooth_it, cheby_order);
  B_G = std::make_unique<ChebyshevSmoother>(h1_fespace.GetComm(), h1_dbc_tdof_list,
                                            cheby_smooth_it, cheby_order);
  B_G->iterative_mode = false;
}

void DistRelaxationSmoother::InitVectors(int nrhs) const
{
  if (nrhs * A->Height() == r.Size())
  {
    return;
  }
  DestroyVectors();
  r.SetSize(nrhs * A->Height());
  x_G.SetSize(nrhs * A_G->Height());
  y_G.SetSize(nrhs * A_G->Height());
  R.SetSize(nrhs);
  X_G.SetSize(nrhs);
  Y_G.SetSize(nrhs);
  for (int j = 0; j < nrhs; j++)
  {
    R[j] = new mfem::Vector(r, j * A->Height(), A->Height());
    X_G[j] = new mfem::Vector(x_G, j * A_G->Height(), A_G->Height());
    Y_G[j] = new mfem::Vector(y_G, j * A_G->Height(), A_G->Height());
  }
}

void DistRelaxationSmoother::DestroyVectors() const
{
  for (int j = 0; j < R.Size(); j++)
  {
    delete R[j];
    delete X_G[j];
    delete Y_G[j];
  }
}

void DistRelaxationSmoother::SetOperator(const mfem::Operator &op,
                                         const mfem::Operator &op_G)
{
  A = &op;
  A_G = &op_G;
  MFEM_VERIFY(A->Height() == G->Height() && A->Width() == G->Height() &&
                  A_G->Height() == G->Width() && A_G->Width() == G->Width(),
              "Invalid operator sizes for DistRelaxationSmoother!");
  height = A->Height();
  width = A->Width();

  // Set up smoothers for A and A_G.
  B->SetOperator(*A);
  B_G->SetOperator(*A_G);
}

}  // namespace palace
