// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_EIGEN_SOLVER_HPP
#define PALACE_DRIVERS_EIGEN_SOLVER_HPP

#include <complex>
#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/eps.hpp"
#include "linalg/operator.hpp"
#include "utils/labels.hpp"

namespace palace
{

class ErrorIndicator;
class Mesh;
class SurfacePostGeometry;
class SurfaceResponseGeometry;

namespace internal
{

struct ResponseCorrectedMass
{
  std::unique_ptr<Operator> real;
  std::unique_ptr<ComplexOperator> op;
};

struct EigenvalueTarget
{
  std::complex<double> shift;
  EigenvalueSolver::WhichType which;
};

// Construct M + R while preserving separate access to its real and imaginary parts.
ResponseCorrectedMass BuildResponseCorrectedMass(const ComplexOperator &mass,
                                                 const Operator &response);

// Return the shift and spectrum selector for an eigenproblem represented in either
// λ = iω or μ = -λ² = ω².
EigenvalueTarget GetEigenvalueTarget(double target, bool lambda_eigenproblem,
                                     EigenSolverBackend backend,
                                     NonlinearEigenSolver nonlinear_type);

// Convert the solver eigenvalue to angular frequency.
std::complex<double> EigenvalueToAngularFrequency(std::complex<double> value,
                                                  bool lambda_eigenproblem);

// Return the maximum-total-weight one-to-one column assignment for a rectangular matrix
// with no more rows than columns.
std::vector<int>
MaximumWeightModeAssignment(const std::vector<std::vector<double>> &weights);

}  // namespace internal

//
// Driver class for eigenmode simulations.
//
class EigenSolver : public BaseSolver
{
private:
  mutable std::shared_ptr<const SurfacePostGeometry> surface_post_geometry;
  mutable std::shared_ptr<const SurfaceResponseGeometry> response_geometry;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
