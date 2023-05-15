// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
#define PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Hiptmair distributive relaxation smoother applying smoothers to both the operator in the
// primary space as well as its projection into an auxiliary space.
// Reference: Hiptmair, Multigrid method for Maxwell's equations, SIAM J. Numer. Anal.
//            (1998).
//
class DistRelaxationSmoother : public mfem::Solver
{
private:
  // Number of smoother iterations.
  const int pc_it;

  // System matrix and its projection Gᵀ A G (not owned).
  const ParOperator *A, *A_G;

  // Discrete gradient matrix.
  std::unique_ptr<ParOperator> G;

  // Point smoother objects for each matrix.
  mutable std::unique_ptr<mfem::Solver> B;
  std::unique_ptr<mfem::Solver> B_G;

  // Temporary vectors for smoother application.
  mutable Vector r, x_G, y_G;

public:
  DistRelaxationSmoother(mfem::ParFiniteElementSpace &nd_fespace,
                         mfem::ParFiniteElementSpace &h1_fespace, int smooth_it,
                         int cheby_smooth_it, int cheby_order);

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for "
               "DistRelaxationSmoother, use the two argument signature instead!");
  }
  void SetOperator(const ParOperator &op, const ParOperator &op_G);

  void Mult(const Vector &x, Vector &y) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_DIST_RELAXATION_SMOOTHER_HPP
