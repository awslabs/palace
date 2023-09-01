// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "distrelaxation.hpp"

#include <mfem.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/rap.hpp"

namespace palace
{

template <typename OperType>
DistRelaxationSmoother<OperType>::DistRelaxationSmoother(
    MPI_Comm comm, const Operator &G, int smooth_it, int cheby_smooth_it, int cheby_order,
    double cheby_sf_max, double cheby_sf_min, bool cheby_4th_kind)
  : Solver<OperType>(), pc_it(smooth_it), G(&G), A(nullptr), A_G(nullptr),
    dbc_tdof_list_G(nullptr)
{
  // Initialize smoothers.
  if (cheby_4th_kind)
  {
    B = std::make_unique<ChebyshevSmoother<OperType>>(comm, cheby_smooth_it, cheby_order,
                                                      cheby_sf_max);
    B_G = std::make_unique<ChebyshevSmoother<OperType>>(comm, cheby_smooth_it, cheby_order,
                                                        cheby_sf_max);
  }
  else
  {
    B = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
        comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_sf_min);
    B_G = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
        comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_sf_min);
  }
  B_G->SetInitialGuess(false);
}

template <typename OperType>
void DistRelaxationSmoother<OperType>::SetOperators(const OperType &op,
                                                    const OperType &op_G)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  MFEM_VERIFY(op.Height() == G->Height() && op.Width() == G->Height() &&
                  op_G.Height() == G->Width() && op_G.Width() == G->Width(),
              "Invalid operator sizes for DistRelaxationSmoother!");
  A = &op;
  A_G = &op_G;
  x_G.SetSize(op_G.Height());
  y_G.SetSize(op_G.Height());
  r_G.SetSize(op_G.Height());
  x_G.UseDevice(true);
  y_G.UseDevice(true);
  r_G.UseDevice(true);

  const auto *PtAP_G = dynamic_cast<const ParOperType *>(&op_G);
  MFEM_VERIFY(PtAP_G,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  dbc_tdof_list_G = PtAP_G->GetEssentialTrueDofs();

  // Set up smoothers for A and A_G.
  B->SetOperator(op);
  B_G->SetOperator(op_G);

  this->height = op.Height();
  this->width = op.Width();
}

namespace
{

inline void RealAddMult(const Operator &op, const Vector &x, Vector &y)
{
  op.AddMult(x, y, 1.0);
}

inline void RealAddMult(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.AddMult(x.Real(), y.Real(), 1.0);
  op.AddMult(x.Imag(), y.Imag(), 1.0);
}

inline void RealMultTranspose(const Operator &op, const Vector &x, Vector &y)
{
  op.MultTranspose(x, y);
}

inline void RealMultTranspose(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.MultTranspose(x.Real(), y.Real());
  op.MultTranspose(x.Imag(), y.Imag());
}

}  // namespace

template <typename OperType>
void DistRelaxationSmoother<OperType>::Mult2(const VecType &x, VecType &y, VecType &r) const
{
  // Apply smoother.
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + B (x - A y)
    B->SetInitialGuess(this->initial_guess || it > 0);
    B->Mult2(x, y, r);

    // y = y + G B_G Gᵀ (x - A y)
    A->Mult(y, r);
    linalg::AXPBY(1.0, x, -1.0, r);
    RealMultTranspose(*G, r, x_G);
    if (dbc_tdof_list_G)
    {
      linalg::SetSubVector(x_G, *dbc_tdof_list_G, 0.0);
    }
    B_G->Mult2(x_G, y_G, r_G);
    RealAddMult(*G, y_G, y);
  }
}

template <typename OperType>
void DistRelaxationSmoother<OperType>::MultTranspose2(const VecType &x, VecType &y,
                                                      VecType &r) const
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
      linalg::SetSubVector(x_G, *dbc_tdof_list_G, 0.0);
    }
    B_G->MultTranspose2(x_G, y_G, r_G);
    RealAddMult(*G, y_G, y);

    // y = y + Bᵀ (x - A y)
    B->MultTranspose2(x, y, r);
  }
}

template class DistRelaxationSmoother<Operator>;
template class DistRelaxationSmoother<ComplexOperator>;

}  // namespace palace
