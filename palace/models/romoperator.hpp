// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

#ifdef MFEM_USE_SINGLE
using EigenVectorX = Eigen::VectorXf;
using EigenVectorXc = Eigen::VectorXcf;
using EigenMatrixX = Eigen::MatrixXf;
using EigenMatrixXc = Eigen::MatrixXcf;
#else
using EigenVectorX = Eigen::VectorXd;
using EigenVectorXc = Eigen::VectorXcd;
using EigenMatrixX = Eigen::MatrixXd;
using EigenMatrixXc = Eigen::MatrixXcd;
#endif

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
private:
  // Reference to HDM discretization (not owned).
  SpaceOperator &space_op;

  // HDM system matrices and excitation RHS.
  std::unique_ptr<ComplexOperator> K, M, C, A2;
  ComplexVector RHS1, RHS2, r;
  bool has_A2, has_RHS1, has_RHS2;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // PROM matrices and vectors.
  EigenMatrixXc Kr, Mr, Cr, Ar;
  EigenVectorXc RHS1r, RHSr;

  // PROM reduced-order basis (real-valued) and active dimension.
  std::vector<Vector> V;
  std::size_t dim_V;
  GmresSolverBase::OrthogType orthog_type;

  // (Complex-valued) upper-trianglar matrix R from orthogonalization of the HDM samples.
  // Minimal rational interpolant (MRI) defined by the vector q of interpolation weights and
  // support points z is used as an error indicator.
  std::vector<ComplexVector> Q;
  std::size_t dim_Q;
  EigenMatrixXc R;
  EigenVectorXc q;
  std::vector<mfem::real_t> z;

public:
  RomOperator(const IoData &iodata, SpaceOperator &space_op, int max_size);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension.
  int GetReducedDimension() const { return dim_V; }

  // Return set of sampled parameter points for basis construction.
  const auto &GetSamplePoints() const { return z; }

  // Assemble and solve the HDM at the specified frequency.
  void SolveHDM(mfem::real_t omega, ComplexVector &u);

  // Add the solution vector to the reduced-order basis and update the PROM.
  void UpdatePROM(mfem::real_t omega, const ComplexVector &u);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional space.
  void SolvePROM(mfem::real_t omega, ComplexVector &u);

  // Compute the location(s) of the maximum error in the range of the previously sampled
  // parameter points.
  std::vector<mfem::real_t> FindMaxError(int N = 1) const;

  // Compute eigenvalue estimates for the current PROM system.
  std::vector<std::complex<mfem::real_t>> ComputeEigenvalueEstimates() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
