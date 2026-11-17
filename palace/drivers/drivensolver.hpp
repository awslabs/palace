// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SOLVER_HPP
#define PALACE_DRIVERS_DRIVEN_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/operator.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class SpaceOperator;
class SurfacePostGeometry;
class SurfaceResponseGeometry;

namespace internal
{

// The system operator of one frequency of the uniform driven sweep. A is applied through
// the finest level of the preconditioner matrix P when SpaceOperator::GetSystemOperator
// provides it (the matrix-free path), else it is assembled from K, C, M and the extra
// system operator A2 (K, C, M are built on first use into the caller's cache). With a
// response mass dM (and stiffness dK), corrected = A - omega^2 dM (+ dK) wraps whichever A
// was built, so the self-consistent corrected solve sees the same system on both paths.
// Members are declared in dependency order: corrected refers to A, A to A2 and P.
struct DrivenSystemOperator
{
  std::unique_ptr<ComplexOperator> A2;
  std::unique_ptr<ComplexOperator> A;
  std::unique_ptr<SumComplexOperator> corrected;
  bool assembled = false;
};

DrivenSystemOperator BuildDrivenSystemOperator(SpaceOperator &space_op, double omega,
                                               const ComplexOperator &P,
                                               std::unique_ptr<ComplexOperator> &K,
                                               std::unique_ptr<ComplexOperator> &C,
                                               std::unique_ptr<ComplexOperator> &M,
                                               const ComplexOperator *response_mass,
                                               const ComplexOperator *response_stiffness);

}  // namespace internal

//
// Driver class for driven terminal simulations.
//
class DrivenSolver : public BaseSolver
{
private:
  mutable std::shared_ptr<const SurfacePostGeometry> surface_post_geometry;
  mutable std::shared_ptr<const SurfaceResponseGeometry> response_geometry;

  ErrorIndicator SweepUniform(SpaceOperator &space_op) const;

  ErrorIndicator SweepAdaptive(SpaceOperator &space_op) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SOLVER_HPP
