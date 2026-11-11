// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
#define PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP

#include <memory>
#include <utility>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace mfem
{

template <typename T>
class Array;

}  // namespace mfem

namespace palace
{

//
// Hiptmair distributive relaxation smoother applying smoothers to both the operator in the
// primary space as well as its projection into an auxiliary space.
// Reference: Hiptmair, Multigrid method for Maxwell's equations, SIAM J. Numer. Anal.
//            (1998).
//
// Each smoothing step applies the smoother of the primary space and the auxiliary space
// correction y += G B_G Gᵀ (x - A y), in this order for Mult and in the reverse order for
// MultTranspose, which need the residual of the iterate in between. With EXACT residuals
// it is computed with the operator. UPDATED residuals are carried through the applications
// of the operator (its single-precision approximation when it has one) by the Chebyshev
// smoother of the primary space and to the auxiliary space correction: the smoother is the
// same, up to the rounding of these applications. With LAGGED residuals, the auxiliary
// space correction of Mult uses the residual before the last step of the Chebyshev
// smoother, and the first Chebyshev step of MultTranspose uses the residual before the
// auxiliary space correction (applied together with that step): this saves both operator
// applications per smoothing step, but MultTranspose is not the transpose of Mult.
//
// The residual updates of the UPDATED and LAGGED residuals apply the single-precision
// approximation of A to the auxiliary space correction G y_G (and the residual of LAGGED
// Mult comes from such applications to the Chebyshev search directions), but rounding the
// curl-curl part of A does not preserve its null space: the error is about ε₃₂ ‖K‖ ‖G y_G‖,
// while the update is ‖A G y_G‖, the much smaller mass part at low frequencies (small k h),
// where y_G is large. Each level with requested residual updates measures the relative
// error of the single-precision operator on gradients when its operators are set, and
// keeps EXACT residuals when it is not small. The error is measured both globally and on
// the auxiliary space scaled by the inverse diagonal of GᵀAG: large terms on gradients
// (such as the inductance of impedance boundaries) dominate the global norms, but not the
// error on the remaining gradients. The auxiliary space error is the root mean square of
// the local errors, which the smallest elements dominate, so its bound is larger. Where
// only it exceeds its bound, Mult computes the residual for the auxiliary space correction
// exactly and MultTranspose keeps the requested residuals: the error breaks down the
// multigrid preconditioner through the residual of the Chebyshev search directions which
// Mult projects on the auxiliary space, while the update after the correction in
// MultTranspose is harmless. The bound 100 comes from driven simulations at low
// frequencies, lagged against exact residuals: up to auxiliary space errors of about 20,
// the lagged residuals needed the same iterations, and one or two more up to 100; larger
// errors increased the iterations, and much larger ones broke down the preconditioner. The
// auxiliary space check can be turned off (SetAuxiliaryGradientCheck): the shift-and-invert
// solves of eigenmode simulations converged with lagged residuals in fewer iterations than
// with exact residuals up to auxiliary space errors of about 1e4.
//
template <OperatorType OperType>
class DistRelaxationSmoother : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of smoother iterations.
  const int pc_it;

  // Discrete gradient matrix (not owned).
  const Operator *G;

  // System matrix and its projection GᵀAG (not owned).
  const OperType *A, *A_G;
  const mfem::Array<int> *dbc_tdof_list_G;

  // Point smoother objects for each matrix.
  mutable std::unique_ptr<Solver<OperType>> B;
  std::unique_ptr<Solver<OperType>> B_G;

  // Chebyshev order of the auxiliary space smoother (<= 0: the order of the primary one).
  int order_aux;

  // Residuals between the primary and auxiliary space smoothers, and the operator for the
  // residual updates of the auxiliary space corrections (A or its single-precision
  // approximation).
  MultigridAuxiliaryResidual aux_residual, aux_residual_req;
  const OperType *A_dir;
  MPI_Comm comm;

  // Relative errors ‖(A_dir - A) G v‖ / ‖A G v‖ and ‖D_G⁻¹ Gᵀ (A_dir - A) G v‖ /
  // ‖D_G⁻¹ Gᵀ A G v‖ (D_G the diagonal of the auxiliary space operator) of the operator for
  // the residual updates on the range of the discrete gradient (negative when not
  // measured), and their bounds for residual updates.
  double grad_error, grad_error_aux;
  static constexpr double max_grad_error = 0.1, max_grad_error_aux = 100.0;

  // Whether Mult computes the residual for the auxiliary space correction exactly while
  // the requested residuals are kept elsewhere (when only grad_error_aux exceeds its
  // bound), and whether grad_error_aux is measured.
  bool exact_aux_input, check_aux;

  // Measure grad_error and grad_error_aux (negative when the auxiliary space smoother has
  // no diagonal) for a random auxiliary space vector v.
  std::pair<double, double> GradientResidualError() const;

  // Resolve the residuals used from the requested ones and grad_error.
  void UpdateAuxiliaryResidual();

  // Temporary vectors for smoother application.
  mutable VecType x_G, y_G, r_G, r, z;

  // Performance tally ids: main and auxiliary smoothers, residual, Gᵀ and G transfers.
  int perf_main = -1, perf_aux = -1, perf_res = -1, perf_gt = -1, perf_g = -1;

public:
  DistRelaxationSmoother(
      MPI_Comm comm, const Operator &G, int smooth_it, int cheby_smooth_it, int cheby_order,
      double cheby_sf_max, double cheby_sf_min, bool cheby_4th_kind,
      int cheby_max_order = 0, int cheby_order_aux = -1,
      MultigridAuxiliaryResidual aux_residual = MultigridAuxiliaryResidual::EXACT);

  void SetPerfTallies(const std::string &prefix, Timer::Index parent) override;

  void SetOperator(const OperType &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "DistRelaxationSmoother, use the two argument signature instead!");
  }

  void SetOperators(const OperType &op, const OperType &op_G);

  const Solver<OperType> *GetSmoother() const { return B.get(); }
  const Solver<OperType> *GetAuxiliarySmoother() const { return B_G.get(); }

  // Set the order of the Chebyshev smoother of the primary space, and of the auxiliary
  // space unless it has its own order.
  void SetSmootherOrder(int order);

  // Order of the Chebyshev smoother of the auxiliary space: its own order (positive), or
  // that of the primary space (zero or negative).
  int GetAuxiliarySmootherOrder() const { return order_aux; }
  void SetAuxiliarySmootherOrder(int order);

  // Residuals between the primary and auxiliary space smoothers (EXACT, UPDATED or
  // LAGGED). The residuals used are EXACT on levels where the single-precision operator is
  // not accurate on gradients (see GetGradientError).
  void SetAuxiliaryResidual(MultigridAuxiliaryResidual residual);
  MultigridAuxiliaryResidual GetAuxiliaryResidual() const { return aux_residual; }

  // Measure the error on the auxiliary space too (the default) or only globally.
  void SetAuxiliaryGradientCheck(bool check);
  bool GetExactAuxiliaryInput() const { return exact_aux_input; }

  // Relative error of the operator for the residual updates on the range of the discrete
  // gradient (zero without a single-precision approximation, negative when not measured),
  // globally and on the auxiliary space.
  double GetGradientError() const { return grad_error; }
  double GetAuxiliaryGradientError() const { return grad_error_aux; }

  void Mult(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    Mult2(x, y, r);
  }

  void MultTranspose(const VecType &x, VecType &y) const override
  {
    if (r.Size() != y.Size())
    {
      r.SetSize(y.Size());
      r.UseDevice(true);
    }
    MultTranspose2(x, y, r);
  }

  void Mult2(const VecType &x, VecType &y, VecType &r) const override;

  void MultTranspose2(const VecType &x, VecType &y, VecType &r) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
