// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "distrelaxation.hpp"

#include <limits>
#include <tuple>
#include <utility>
#include <mfem.hpp>
#include "linalg/chebyshev.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/perfstats.hpp"

namespace palace
{

template <OperatorType OperType>
DistRelaxationSmoother<OperType>::DistRelaxationSmoother(
    MPI_Comm comm, const Operator &G, int smooth_it, int cheby_smooth_it, int cheby_order,
    double cheby_sf_max, double cheby_sf_min, bool cheby_4th_kind, int cheby_max_order,
    int cheby_order_aux, MultigridAuxiliaryResidual aux_residual)
  : Solver<OperType>(), pc_it(smooth_it), G(&G), A(nullptr), A_G(nullptr),
    dbc_tdof_list_G(nullptr), order_aux(cheby_order_aux),
    aux_residual(MultigridAuxiliaryResidual::EXACT),
    aux_residual_req(MultigridAuxiliaryResidual::EXACT), A_dir(nullptr), comm(comm),
    grad_error(-1.0), grad_error_aux(-1.0), exact_aux_input(false), check_aux(true)
{
  // Initialize smoothers. Unless given its own order, the auxiliary space smoother uses the
  // order chosen for the primary space.
  const int order_G = (order_aux > 0) ? order_aux : cheby_order;
  if (cheby_4th_kind)
  {
    B = std::make_unique<ChebyshevSmoother<OperType>>(comm, cheby_smooth_it, cheby_order,
                                                      cheby_sf_max, cheby_max_order);
    B_G = std::make_unique<ChebyshevSmoother<OperType>>(comm, cheby_smooth_it, order_G,
                                                        cheby_sf_max);
  }
  else
  {
    B = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
        comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_sf_min);
    B_G = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
        comm, cheby_smooth_it, order_G, cheby_sf_max, cheby_sf_min);
  }
  B_G->SetInitialGuess(false);
  SetAuxiliaryResidual(aux_residual);
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetOperators(const OperType &op,
                                                    const OperType &op_G)
{
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;

  MFEM_VERIFY(op.Height() == G->Height() && op.Width() == G->Height() &&
                  op_G.Height() == G->Width() && op_G.Width() == G->Width(),
              "Invalid operator sizes for DistRelaxationSmoother!");
  A = &op;
  A_G = &op_G;
  A_dir = &GetDirectionOperator(op);
  x_G.SetSize(op_G.Height());
  y_G.SetSize(op_G.Height());
  r_G.SetSize(op_G.Height());
  x_G.UseDevice(true);
  y_G.UseDevice(true);
  r_G.UseDevice(true);
  if (aux_residual_req != MultigridAuxiliaryResidual::EXACT && z.Size() != op.Height())
  {
    z.SetSize(op.Height());
    z.UseDevice(true);
  }

  const auto *PtAP_G = dynamic_cast<const ParOperType *>(&op_G);
  MFEM_VERIFY(PtAP_G,
              "ChebyshevSmoother requires a ParOperator or ComplexParOperator operator!");
  dbc_tdof_list_G = PtAP_G->GetEssentialTrueDofs();

  // Set up smoothers for A and A_G.
  B->SetOperator(op);
  const auto *cheby = dynamic_cast<const ChebyshevSmoother<OperType> *>(B.get());
  if (cheby && order_aux <= 0)
  {
    static_cast<ChebyshevSmoother<OperType> &>(*B_G).SetOrder(cheby->GetOrder());
  }
  B_G->SetOperator(op_G);

  // The accuracy of the residual updates depends on the operator (for example, on the
  // frequency of a driven solve): measure it again.
  grad_error = grad_error_aux = -1.0;
  UpdateAuxiliaryResidual();

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetSmootherOrder(int order)
{
  auto *cheby = dynamic_cast<ChebyshevSmoother<OperType> *>(B.get());
  MFEM_VERIFY(cheby, "SetSmootherOrder requires a 4th-kind Chebyshev smoother!");
  cheby->SetOrder(order);
  if (order_aux <= 0)
  {
    static_cast<ChebyshevSmoother<OperType> &>(*B_G).SetOrder(order);
  }
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetAuxiliarySmootherOrder(int order)
{
  auto *cheby_G = dynamic_cast<ChebyshevSmoother<OperType> *>(B_G.get());
  MFEM_VERIFY(cheby_G, "SetAuxiliarySmootherOrder requires a 4th-kind Chebyshev smoother!");
  order_aux = order;
  if (order > 0)
  {
    cheby_G->SetOrder(order);
  }
  else if (const auto *cheby = dynamic_cast<const ChebyshevSmoother<OperType> *>(B.get()))
  {
    cheby_G->SetOrder(cheby->GetOrder());
  }
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetAuxiliaryResidual(
    MultigridAuxiliaryResidual residual)
{
  // An unresolved default keeps the exact residuals.
  aux_residual_req = (residual == MultigridAuxiliaryResidual::DEFAULT)
                         ? MultigridAuxiliaryResidual::EXACT
                         : residual;
  MFEM_VERIFY(aux_residual_req == MultigridAuxiliaryResidual::EXACT || B->UpdatesResidual(),
              "Updated or lagged residuals require a smoother which updates its residual!");
  if (aux_residual_req != MultigridAuxiliaryResidual::EXACT && A && z.Size() != A->Height())
  {
    z.SetSize(A->Height());
    z.UseDevice(true);
  }
  UpdateAuxiliaryResidual();
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetAuxiliaryGradientCheck(bool check)
{
  if (check != check_aux)
  {
    check_aux = check;
    grad_error = grad_error_aux = -1.0;
    UpdateAuxiliaryResidual();
  }
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::UpdateAuxiliaryResidual()
{
  aux_residual = aux_residual_req;
  exact_aux_input = false;
  if (aux_residual == MultigridAuxiliaryResidual::EXACT || !A || !A_G)
  {
    return;
  }
  if (grad_error < 0.0)
  {
    std::tie(grad_error, grad_error_aux) = GradientResidualError();
  }
  if (!(grad_error <= max_grad_error))
  {
    aux_residual = MultigridAuxiliaryResidual::EXACT;
  }
  else if (!(grad_error_aux <= max_grad_error_aux))
  {
    exact_aux_input = true;
  }
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

inline void RealMult(const Operator &op, const Vector &x, Vector &y)
{
  op.Mult(x, y);
}

inline void RealMult(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.Mult(x.Real(), y.Real());
  op.Mult(x.Imag(), y.Imag());
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

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::SetPerfTallies(const std::string &prefix,
                                                      Timer::Index parent)
{
  perf_main = perf::RegisterTally(prefix + "/Main", parent);
  perf_aux = perf::RegisterTally(prefix + "/Aux", parent);
  perf_res = perf::RegisterTally(prefix + "/Residual", parent);
  perf_gt = perf::RegisterTally(prefix + "/GT", parent);
  perf_g = perf::RegisterTally(prefix + "/G", parent);
  B->SetPerfTallies(prefix + "/Main", parent);
  B_G->SetPerfTallies(prefix + "/Aux", parent);
}

template <OperatorType OperType>
std::pair<double, double> DistRelaxationSmoother<OperType>::GradientResidualError() const
{
  // The single-precision approximation is local (a process without elements it applies to
  // has none) and the measurement is collective.
  int single = (A_dir != A);
  Mpi::GlobalMax(1, &single, comm);
  if (!single)
  {
    return {0.0, check_aux ? 0.0 : -1.0};
  }
  VecType v(A_G->Height()), a(A->Height()), b(A->Height());
  v.UseDevice(true);
  a.UseDevice(true);
  b.UseDevice(true);
  linalg::SetRandomReal(comm, v);
  if (dbc_tdof_list_G)
  {
    linalg::SetSubVector(v, *dbc_tdof_list_G, 0.0);
  }
  RealMult(*G, v, z);
  A->Mult(z, a);
  A_dir->Mult(z, b);
  const double norm_a = linalg::Norml2(comm, a);
  linalg::AXPY(-1.0, a, b);
  const double norm_e = linalg::Norml2(comm, b);
  constexpr double inf = std::numeric_limits<double>::infinity();
  const double rho = (norm_a > 0.0) ? norm_e / norm_a : inf;

  // The norms above are dominated by the largest terms of A on gradients, for example the
  // inductive boundary terms of impedance boundaries, while on the remaining gradients the
  // error is the rounding of the curl-curl part against the much smaller mass part. Measure
  // the error also on the auxiliary space, where each entry is scaled by the inverse
  // diagonal of GᵀAG (the scaling of the auxiliary space smoother which amplifies it): the
  // root mean square of the local errors.
  if (!check_aux)
  {
    return {rho, -1.0};
  }
  const VecType *dinv_G = nullptr;
  if (const auto *c = dynamic_cast<const ChebyshevSmoother<OperType> *>(B_G.get()))
  {
    dinv_G = &c->GetInverseDiagonal();
  }
  else if (const auto *c =
               dynamic_cast<const ChebyshevSmoother1stKind<OperType> *>(B_G.get()))
  {
    dinv_G = &c->GetInverseDiagonal();
  }
  if (!dinv_G || dinv_G->Size() != A_G->Height())
  {
    return {rho, -1.0};
  }
  const BaseDiagonalOperator<OperType> Dinv_G(*dinv_G);
  auto ScaledNorm = [&](const VecType &x)
  {
    RealMultTranspose(*G, x, x_G);
    if (dbc_tdof_list_G)
    {
      linalg::SetSubVector(x_G, *dbc_tdof_list_G, 0.0);
    }
    Dinv_G.Mult(x_G, y_G);
    return linalg::Norml2(comm, y_G);
  };
  const double norm_a_G = ScaledNorm(a), norm_e_G = ScaledNorm(b);
  const double rho_G = (norm_a_G > 0.0) ? norm_e_G / norm_a_G : inf;
  return {rho, rho_G};
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::Mult2(const VecType &x, VecType &y, VecType &r) const
{
  // Apply smoother.
  const bool update =
                 (aux_residual == MultigridAuxiliaryResidual::UPDATED) && !exact_aux_input,
             lagged =
                 (aux_residual == MultigridAuxiliaryResidual::LAGGED) && !exact_aux_input;
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + B (x - A y). With updated residuals, B returns the residual of its result,
    // starting from the one updated by the previous iteration. Otherwise, r holds the
    // residual before the last step of B.
    B->SetInitialGuess(this->initial_guess || it > 0);
    B->SetResidual(update && it > 0, update);
    {
      perf::ScopedTally t(perf_main);
      B->Mult2(x, y, r);
    }

    // y = y + G B_G Gᵀ (x - A y), with the lagged residual or after computing it.
    if (!update && !lagged)
    {
      perf::ScopedTally t(perf_res);
      A->Mult(y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    {
      perf::ScopedTally t(perf_gt);
      RealMultTranspose(*G, r, x_G);
      if (dbc_tdof_list_G)
      {
        linalg::SetSubVector(x_G, *dbc_tdof_list_G, 0.0);
      }
    }
    {
      perf::ScopedTally t(perf_aux);
      B_G->Mult2(x_G, y_G, r_G);
    }
    if (update && it < pc_it - 1)
    {
      {
        perf::ScopedTally t(perf_g);
        RealMult(*G, y_G, z);
        y += z;
      }
      perf::ScopedTally t(perf_res);
      A_dir->AddMult(z, r, -1.0);
    }
    else
    {
      perf::ScopedTally t(perf_g);
      RealAddMult(*G, y_G, y);
    }
  }
}

template <OperatorType OperType>
void DistRelaxationSmoother<OperType>::MultTranspose2(const VecType &x, VecType &y,
                                                      VecType &r) const
{
  // Apply transpose.
  const bool update = (aux_residual == MultigridAuxiliaryResidual::UPDATED),
             lagged = (aux_residual == MultigridAuxiliaryResidual::LAGGED);
  for (int it = 0; it < pc_it; it++)
  {
    // y = y + G B_Gᵀ Gᵀ (x - A y). With updated residuals, the residual is the one returned
    // by B in the previous iteration.
    if (this->initial_guess || it > 0)
    {
      if (!update || it == 0 || exact_aux_input)
      {
        perf::ScopedTally t(perf_res);
        A->Mult(y, r);
        linalg::AXPBY(1.0, x, -1.0, r);
      }
      perf::ScopedTally t(perf_gt);
      RealMultTranspose(*G, r, x_G);
    }
    else
    {
      y = 0.0;
      if (update || lagged)
      {
        r = x;
      }
      perf::ScopedTally t(perf_gt);
      RealMultTranspose(*G, x, x_G);
    }
    if (dbc_tdof_list_G)
    {
      linalg::SetSubVector(x_G, *dbc_tdof_list_G, 0.0);
    }
    {
      perf::ScopedTally t(perf_aux);
      B_G->MultTranspose2(x_G, y_G, r_G);
    }

    // y = y + Bᵀ (x - A y). The correction G y_G is added to the iterate and its residual
    // before B, or by B together with its first search direction for lagged residuals.
    B->SetInitialGuess(true);
    if (update || lagged)
    {
      {
        perf::ScopedTally t(perf_g);
        RealMult(*G, y_G, z);
      }
      if (update)
      {
        y += z;
        perf::ScopedTally t(perf_res);
        A_dir->AddMult(z, r, -1.0);
        B->SetResidual(true, it < pc_it - 1);
      }
      else
      {
        B->SetResidual(true, false, &z);
      }
    }
    else
    {
      perf::ScopedTally t(perf_g);
      RealAddMult(*G, y_G, y);
      B->SetResidual(false, false);
    }
    perf::ScopedTally t(perf_main);
    B->MultTranspose2(x, y, r);
  }
}

template class DistRelaxationSmoother<Operator>;
template class DistRelaxationSmoother<ComplexOperator>;

}  // namespace palace
