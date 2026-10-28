// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SPACE_EIGEN_SOLVER_HPP
#define PALACE_MODELS_SPACE_EIGEN_SOLVER_HPP

#include <complex>
#include <functional>
#include <memory>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class Interpolation;
class IoData;
class SpaceOperator;

template <typename VecType>
class DivFreeSolver;

//
// Eigenvalue solve for the frequency-domain system assembled by SpaceOperator,
//
//   T(λ) u = (K + λ C + λ² M + A2(λ)) u = 0,   λ = iω,
//
// for the eigenpairs nearest config["Solver"]["Eigenmode"]["Target"]. Without
// frequency-dependent boundaries this is the linear problem K u = ω² M u or, with damping,
// the quadratic polynomial problem. Frequency-dependent boundaries (wave ports,
// second-order farfield, surface conductivity, rational surface impedance) contribute
// A2(λ) and make the problem nonlinear. It is then solved with the configured
// "NonlinearType": Hybrid (a polynomial seed pencil from a Newton interpolation of A2 over
// [Target, TargetUpper], refined by quasi-Newton iterations on the exact T(λ)) or SLP
// (SLEPc successive linear problems). Shared by the eigenmode simulation driver and the
// eigenmode enrichment of the circuit-synthesis reduced-order model, so that both solve the
// same eigenvalue problem.
//
class SpaceEigenSolver
{
private:
  using A2Function = std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>;
  using PreconditionerFunction = std::function<std::unique_ptr<ComplexOperator>(
      std::complex<double>, std::complex<double>, std::complex<double>,
      std::complex<double>)>;

  const IoData &iodata;
  SpaceOperator &space_op;

  // System matrices and linear solver for the spectral transformation, owned by the caller.
  const ComplexOperator &K, *C, &M;
  ComplexKspSolver &ksp;

  // A2(λ) as the sparse boundary matrix (the polynomial seed interpolates it), and with the
  // matrix-free wave-port modal correction W (the exact operator for SLP and quasi-Newton).
  A2Function funcA2, funcA2_full;
  PreconditionerFunction funcP;

  NonlinearEigenSolver nonlinear_type;
  bool has_A2, quadratic, nonlinear_slp;

  // Operators referenced by the eigenvalue solver. Declared before it so that it is
  // destroyed first; A2_0/A2_1/A2_2 reference the interpolation operator DAG and Kp/Cp/Mp
  // reference them.
  std::unique_ptr<ComplexOperator> A2;
  std::unique_ptr<Interpolation> interp_op;
  std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2, Kp, Cp, Mp;
  std::unique_ptr<Operator> KM;
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  std::unique_ptr<ComplexOperator> A2_shift, A, P;
  std::unique_ptr<EigenvalueSolver> eigen;

public:
  // K, C (null without damping), M, and the linear solver are owned by the caller and must
  // outlive this object. The constructor freezes the wave-port modal reference at the
  // target and configures the eigenvalue solver, its spectral transformation, and the
  // linear solver operators; Solve runs it. The linear solver's operators are reset by the
  // solve.
  SpaceEigenSolver(const IoData &iodata, SpaceOperator &space_op, const ComplexOperator &K,
                   const ComplexOperator *C, const ComplexOperator &M,
                   ComplexKspSolver &ksp);
  ~SpaceEigenSolver();

  // Solve the eigenvalue problem, including the quasi-Newton refinement of the Hybrid
  // nonlinear strategy, and return the number of converged eigenpairs.
  int Solve();

  // Whether the eigenvalues follow the λ = iω convention (quadratic or nonlinear problem)
  // rather than μ = ω² (linear problem); see EigenvalueToOmega.
  bool IsQuadratic() const { return quadratic; }

  // Access the converged eigenpairs, as for EigenvalueSolver.
  std::complex<double> GetEigenvalue(int i) const { return eigen->GetEigenvalue(i); }
  std::complex<double> GetOmega(int i) const
  {
    return EigenvalueToOmega(GetEigenvalue(i), quadratic);
  }
  void GetEigenvector(int i, ComplexVector &x) const { eigen->GetEigenvector(i, x); }
  double GetError(int i, EigenvalueSolver::ErrorType type) const
  {
    return eigen->GetError(i, type);
  }

  // Normalize the first num_eig eigenvectors to unit M-norm (unit electric field energy),
  // unless the eigenvalue solver already used the M-inner product ("MassOrthogonal").
  void NormalizeEigenvectors(int num_eig);
};

}  // namespace palace

#endif  // PALACE_MODELS_SPACE_EIGEN_SOLVER_HPP
