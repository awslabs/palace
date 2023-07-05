// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <memory>
#include <random>
#include <set>
#include <vector>
#include <Eigen/Dense>
#include "linalg/curlcurl.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
private:
  // Reference to HDM discretization (not owned).
  SpaceOperator &spaceop;

  // HDM system matrices and excitation RHS.
  std::unique_ptr<ComplexOperator> K, M, C, A2;
  ComplexVector RHS1, RHS2;
  bool has_A2, has_RHS2;

  // Working storage for HDM vectors.
  ComplexVector r, w, z;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // Linear solver for inner product solves for error metric.
  std::unique_ptr<CurlCurlMassSolver> kspKM;

  // PROM matrices and vectors.
  Eigen::MatrixXcd Kr, Mr, Cr, Ar;
  Eigen::VectorXcd RHS1r, RHSr;

  // PROM reduced-order basis (real-valued) and active dimension.
  std::vector<Vector> V;
  int dim_V;
  bool orthog_mgs;

  // Data structures for parameter domain sampling.
  std::set<double> PS, P_m_PS;
  std::default_random_engine engine;

public:
  RomOperator(const IoData &iodata, SpaceOperator &sp);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension.
  int GetReducedDimension() const { return dim_V; }

  // Return set of sampled parameter points for basis construction.
  const std::set<double> &GetSampleFrequencies() const { return PS; }

  // Initialize the parameter domain P = {ω_L, ω_L + δ, ..., ω_R}. Also sets the maximum
  // number of sample points for the PROM construction.
  void Initialize(double start, double delta, int num_steps, int max_dim);

  // Assemble and solve the HDM at the specified frequency.
  void SolveHDM(double omega, ComplexVector &e);

  // Add the solution vector to the reduced-order basis and update the PROM.
  void AddHDMSample(double omega, ComplexVector &e);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional solution space.
  void AssemblePROM(double omega);
  void SolvePROM(ComplexVector &e);

  // Compute the error metric for the PROM at the specified frequency.
  double ComputeError(double omega);

  // Compute the maximum error over a randomly sampled set of candidate points. Returns the
  // maximum error and its correcponding frequency, as well as the number of candidate
  // points used (if fewer than those availble in the unsampled parameter domain).
  double ComputeMaxError(int num_cand, double &omega_star);
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
