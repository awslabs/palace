// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "distrelaxation.hpp"

#include <mfem.hpp>
#include <general/forall.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/rap.hpp"

namespace palace
{

template <typename OperType>
DistRelaxationSmoother<OperType>::DistRelaxationSmoother(
    mfem::ParFiniteElementSpace &nd_fespace, mfem::ParFiniteElementSpace &h1_fespace,
    int smooth_it, int cheby_smooth_it, int cheby_order)
  : Solver<OperType>(), pc_it(smooth_it), A(nullptr), A_G(nullptr), dbc_tdof_list_G(nullptr)
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
  B = std::make_unique<ChebyshevSmoother<OperType>>(cheby_smooth_it, cheby_order);
  B_G = std::make_unique<ChebyshevSmoother<OperType>>(cheby_smooth_it, cheby_order);
  B_G->SetInitialGuess(false);
}

template <typename OperType>
void DistRelaxationSmoother<OperType>::SetOperators(const OperType &op,
                                                    const OperType &op_G)
{
  typedef typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                    ComplexParOperator, ParOperator>::type ParOperType;

  MFEM_VERIFY(op.Height() == G->Height() && op.Width() == G->Height() &&
                  op_G.Height() == G->Width() && op_G.Width() == G->Width(),
              "Invalid operator sizes for DistRelaxationSmoother!");
  A = &op;
  A_G = &op_G;

  const auto *PtAP_G = dynamic_cast<const ParOperType *>(&op_G);
  MFEM_VERIFY(PtAP_G,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  dbc_tdof_list_G = PtAP_G->GetEssentialTrueDofs();

  r.SetSize(A->Height());
  x_G.SetSize(A_G->Height());
  y_G.SetSize(A_G->Height());

  // Set up smoothers for A and A_G.
  B->SetOperator(*A);
  B_G->SetOperator(*A_G);
}

namespace
{

inline void RealAddMult(Operator &op, const Vector &x, Vector &y)
{
  op.AddMult(x, y, 1.0);
}

inline void RealAddMult(Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.AddMult(x.Real(), y.Real(), 1.0);
  op.AddMult(x.Imag(), y.Imag(), 1.0);
}

inline void RealMultTranspose(Operator &op, const Vector &x, Vector &y)
{
  op.MultTranspose(x, y);
}

inline void RealMultTranspose(Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.MultTranspose(x.Real(), y.Real());
  op.MultTranspose(x.Imag(), y.Imag());
}

}  // namespace

template <typename OperType>
void DistRelaxationSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  // Apply smoother.
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + B (x - A y)
    B->SetInitialGuess(this->initial_guess || it > 0);
    B->Mult(x, y);

    // y = y + G B_G Gᵀ (x - A y)
    A->Mult(y, r);
    linalg::AXPBY(1.0, x, -1.0, r);
    RealMultTranspose(*G, r, x_G);
    if (dbc_tdof_list_G)
    {
      x_G.SetSubVector(*dbc_tdof_list_G, 0.0);
    }
    B_G->Mult(x_G, y_G);
    RealAddMult(*G, y_G, y);
  }
}

template <typename OperType>
void DistRelaxationSmoother<OperType>::MultTranspose(const VecType &x, VecType &y) const
{
  // Apply transpose.
  B->SetInitialGuess(true);
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + G B_Gᵀ Gᵀ (x - A y)
    if (this->initial_guess || it > 0)
    {
      A->Mult(y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
      RealMultTranspose(*G, r, x_G);
    }
    else
    {
      y = 0.0;
      RealMultTranspose(*G, x, x_G);
    }
    if (dbc_tdof_list_G)
    {
      x_G.SetSubVector(*dbc_tdof_list_G, 0.0);
    }
    B_G->MultTranspose(x_G, y_G);
    RealAddMult(*G, y_G, y);

    // y = y + Bᵀ (x - A y)
    B->MultTranspose(x, y);
  }
}

template class DistRelaxationSmoother<Operator>;
template class DistRelaxationSmoother<ComplexOperator>;

}  // namespace palace
