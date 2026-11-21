// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_AMG_HPP
#define PALACE_LINALG_AMG_HPP

#include <memory>
#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for Hypre's BoomerAMG solver.
//
class BoomerAmgSolver : public mfem::HypreBoomerAMG
{
private:
  struct ParVectorDeleter
  {
    void operator()(hypre_ParVector *v) const { hypre_ParVectorDestroy(v); }
  };

  // Number of V-cycles per application.
  const int cycle_it;

  // Multicomponent right-hand side and solution for ArrayMult.
  mutable std::unique_ptr<hypre_ParVector, ParVectorDeleter> B_mv, X_mv;

public:
  // Use agg_levels levels of aggressive coarsening (0 to disable).
  BoomerAmgSolver(int cycle_it = 1, int smooth_it = 1, int agg_levels = 1, int print = 0);
  BoomerAmgSolver(const IoData &iodata, bool coarse_solver, int print)
    : BoomerAmgSolver(coarse_solver ? 1 : iodata.solver.linear.mg_cycle_it,
                      iodata.solver.linear.mg_smooth_it,
                      iodata.solver.linear.amg_agg_coarsen, print)
  {
  }

  // Unless the initial guess is used, zero it and mark it as zero, which lets the first
  // relaxation skip the exchange of its off-processor values (same result).
  void Mult(const mfem::HypreParVector &b, mfem::HypreParVector &x) const override;
  using mfem::HypreBoomerAMG::Mult;

  void SetOperator(const mfem::Operator &op) override;

  // Apply the V-cycle to several right-hand sides at once (for example the real and
  // imaginary parts of a complex vector), with one message per neighbor for all of them,
  // when MultiComponentSupported(). Each component sees the same operations as in Mult.
  // Otherwise, apply it to each right-hand side separately.
  void ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                 mfem::Array<mfem::Vector *> &Y) const override;

  // Whether hypre applies BoomerAMG to multicomponent vectors with the relaxation and
  // coarse solve used here: with the hypre patch for them, on the host and with one
  // thread.
  static bool MultiComponentSupported();
};

}  // namespace palace

#endif  // PALACE_LINALG_AMG_HPP
