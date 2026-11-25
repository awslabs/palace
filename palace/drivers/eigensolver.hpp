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
class VectorFieldTransfer;

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

// Construct K + dK (the magnetic inductive defect, mx-sc-magnetic DESIGN-v2 E3): the twin
// of BuildResponseCorrectedMass on the stiffness side (the essential diagonal of K kept).
ResponseCorrectedMass BuildResponseCorrectedStiffness(const ComplexOperator &stiffness,
                                                      const Operator &response_stiffness);

// The response-corrected twins of the stiffness operators: K + dK and, when the HYBRID
// polynomial seed exists (Kp = K + A2_0), Kp + dK. Kp is null without that seed (no A2
// term); then no corrected Kp is built (OI-UPSYNC-NULL-KP: it was built from the null Kp).
struct ResponseCorrectedStiffnesses
{
  ResponseCorrectedMass K, Kp;
};
ResponseCorrectedStiffnesses
BuildResponseCorrectedStiffnesses(const ComplexOperator &stiffness,
                                  const ComplexOperator *seed_stiffness,
                                  const Operator &response_stiffness);

// Refuse the self-consistent response-corrected eigenmode solve with the magnetic defect dK
// under the SLP nonlinear eigensolver (the exact A2 of wave ports or absorbing boundaries
// beside K + dK): no case exercises that operator (OI-EIGEN-SLP-MAGNETIC-SELF-CONSISTENT).
void VerifySelfConsistentMagneticNonlinearSolver(bool magnetic_defect,
                                                 NonlinearEigenSolver nonlinear_type);

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
  // The impedance-metal rule candidate's own geometry cache (RECORD only; unused under
  // Off).
  mutable std::shared_ptr<const SurfaceResponseGeometry> impedance_candidate_geometry;

  // Combination of the eigenvectors of the previous adaptation iteration (and of λ x for
  // quadratic problems), kept through the refinement and rebalancing of the finest mesh to
  // start the next eigenvalue solve, and the number of eigenvectors in it.
  mutable std::unique_ptr<VectorFieldTransfer> prev_eigenvectors;
  mutable int prev_num_eigenvectors = 0;

  // Number of Krylov cycles of the last solve without a warm start, and number of solves so
  // far (the initial one, then one per adaptation iteration).
  mutable int last_cold_its = -1;
  mutable int num_solves = 0;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

  void UpdateAfterMeshChange(Mesh &mesh) const override;

public:
  EigenSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
              const char *git_tag = nullptr);
  ~EigenSolver() override;

  // CUT-PLAN step 4b (lane v1step4b): PreflightMaxwellResolutionCensus on the first mesh.
  void PreflightResolutionCensus(std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_EIGEN_SOLVER_HPP
