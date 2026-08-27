// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_REDUCED_MODEL_HPP
#define PALACE_MODELS_WAVE_PORT_REDUCED_MODEL_HPP

#include <complex>
#include <memory>
#include <optional>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/modeoperatorassembly.hpp"

namespace palace
{

// Per-wave-port projection model for the boundary-mode generalized eigenproblem. Exact
// solve orchestration and parent-communicator coordination remain in ModeEigenSolver.
class WavePortReducedModel
{
public:
  struct Stats
  {
    std::size_t exact_solves = 0;
    std::size_t reduced_solves = 0;
    std::size_t fallbacks = 0;
    std::size_t periodic_checks = 0;
    std::size_t offline_basis_rank = 0;
    std::size_t online_basis_cap = 0;
    double worst_residual = 0.0;
  };

  WavePortReducedModel(int num_modes, double eig_tol, int nd_size, int h1_size,
                       MPI_Comm solver_comm, bool wave_port,
                       const mode_assembly::ModeOperatorModel &mode_op_model,
                       const ComplexOperator &opB);

  void ConfigureTraining(std::size_t max_basis_size);
  void Enable(bool enable, double adaptive_tol, const ComplexOperator *truth_operator,
              std::complex<double> truth_omega, double truth_sigma);
  void ResetSolution() { has_solution = false; }

  bool IsEnabled() const { return enabled; }
  bool IsReady() const
  {
    return enabled && affine_model_ready && solver_comm != MPI_COMM_NULL &&
           basis.size() >= static_cast<std::size_t>(num_modes);
  }
  bool ExactCheckDue() const { return solves_since_exact >= EXACT_CHECK_INTERVAL; }
  bool TrySolve(double omega, double sigma);

  void ObserveExactEigenvectors(int num_converged, const EigenvalueSolver &eigen,
                                const std::vector<int> &mode_perm);
  void RecordExactSolve();
  void RecordReducedSolve();
  void RecordFallback() { stats.fallbacks++; }
  void RecordPeriodicCheck() { stats.periodic_checks++; }

  const Stats &GetStats() const { return stats; }
  std::size_t GetBasisSize() const { return basis.size(); }
  double GetTolerance() const { return tolerance; }
  bool HasBasis() const { return !basis.empty(); }
  bool HasSolution() const { return has_solution; }
  std::complex<double> GetEigenvalue(int i) const { return eigenvalues.at(i); }
  void GetEigenvector(int i, ComplexVector &x) const { x = eigenvectors.at(i); }
  double GetBackwardError(int i) const { return errors.at(i); }

private:
  struct ProjectedComponent
  {
    const mode_assembly::OperatorComponent *source;
    std::vector<ComplexVector> actions;
    Eigen::MatrixXcd projection;
  };

  int num_modes;
  double eig_tol;
  int nd_size, h1_size;
  MPI_Comm solver_comm;
  bool wave_port;
  const mode_assembly::ModeOperatorModel &mode_op_model;
  const ComplexOperator &opB;

  bool enabled = false;
  bool has_solution = false;
  bool affine_model_ready = false;
  bool basis_cap_warned = false;
  bool gram_direct_self_checked = false;
  bool gram_residual_trusted = true;
  std::size_t basis_cap = 0;
  std::size_t solves_since_exact = 0;
  static constexpr std::size_t EXACT_CHECK_INTERVAL = 20;
  double tolerance = 0.0;
  Stats stats;

  std::vector<ComplexVector> basis;
  std::vector<std::complex<double>> eigenvalues;
  std::vector<ComplexVector> eigenvectors;
  std::vector<double> errors;
  std::vector<ProjectedComponent> components;
  std::vector<ComplexVector> B_actions;
  Eigen::MatrixXcd B_projection;
  Eigen::MatrixXcd action_gram;

  void BuildAffineModel(const ComplexOperator *truth_operator,
                        std::complex<double> truth_omega, double truth_sigma);
  void UpdateProjection(std::size_t old_basis_size);
  void BuildActionGram();
  double ValidateAffineModel(const ComplexOperator &truth_operator, double omega,
                             double sigma) const;
  std::optional<double>
  EvaluateGramResidual(const Eigen::VectorXcd &y, std::complex<double> lambda,
                       const std::vector<std::complex<double>> &coefficients);
  double
  EvaluateDirectResidual(const Eigen::VectorXcd &y, std::complex<double> lambda,
                         const std::vector<std::complex<double>> &coefficients) const;
  bool SolveFromGram(double sigma, const Eigen::MatrixXcd &Ar,
                     const std::vector<std::complex<double>> &coefficients);
  bool AddBasisVector(const ComplexVector &x);
};

}  // namespace palace

#endif  // PALACE_MODELS_WAVE_PORT_REDUCED_MODEL_HPP
