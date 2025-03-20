// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <complex>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/strongtype.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

// Class for handling minimal-rational interpolation of solutions in frequency space. Used
// as an error indicator and to efficiently selecting next frequency sample points in PROM
// construciton. Each excitation gets a separate MRI, so sample frequencies are not shared.
class MinimalRationInterpolation
{
private:
  // (Complex-valued) upper-trianglar matrix R from orthogonalization of the HDM samples.
  // Minimal rational interpolant (MRI) defined by the vector q of interpolation weights and
  // support points z is used as an error indicator.
  std::vector<ComplexVector> Q;
  std::size_t dim_Q = 0;
  Eigen::MatrixXcd R;
  Eigen::VectorXcd q;
  std::vector<double> z;

public:
  MinimalRationInterpolation(int max_size);
  void AddSolutionSample(double omega, const ComplexVector &u,
                         const SpaceOperator &space_op,
                         GmresSolverBase::OrthogType orthog_type);
  std::vector<double> FindMaxError(int N) const;

  const auto &GetSamplePoints() const { return z; }
};

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
private:
  // Reference to HDM discretization (not owned).
  SpaceOperator &space_op;

  // Used for constructing & resuse of RHS1
  ExcitationIdx excitation_idx_cache = ExcitationIdx(0);

  // HDM system matrices and excitation RHS.
  std::unique_ptr<ComplexOperator> K, M, C, A2;
  ComplexVector RHS1, RHS2, r;
  // Defaults: will be toggeled by SetExcitationIndex & SolveHDM
  bool has_A2 = true;
  bool has_RHS1 = true;
  bool has_RHS2 = true;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // PROM matrices and vectors.
  Eigen::MatrixXcd Kr, Mr, Cr, Ar;
  Eigen::VectorXcd RHS1r;
  Eigen::VectorXcd RHSr;
  Eigen::VectorXd voltage_norm_H;
  std::vector<std::string> v_node_label;

  // PROM reduced-order basis (real-valued) and active dimension.
  std::vector<Vector> V;
  std::size_t dim_V = 0;
  GmresSolverBase::OrthogType orthog_type;

  // MRIs: one for each excitation index
  std::map<ExcitationIdx, MinimalRationInterpolation> mri;

public:
  RomOperator(const IoData &iodata, SpaceOperator &space_op, int max_size_per_excitation);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension.
  auto GetReducedDimension() const { return dim_V; }

  // Return set of sampled parameter points for basis construction.
  const auto &GetSamplePoints(ExcitationIdx excitation_idx) const
  {
    return mri.at(excitation_idx).GetSamplePoints();
  }

  std::tuple<const Eigen::MatrixXcd &, const Eigen::MatrixXcd &, const Eigen::MatrixXcd &,
             const Eigen::VectorXd &, const std::vector<std::string> &>
  GetReducedMatrices() const
  {
    return std::tie(Kr, Mr, Cr, voltage_norm_H, v_node_label);
  }

  // Set excitation index to build corresponding RHS vector (linear in frequency part)
  void SetExcitationIndex(ExcitationIdx excitation_idx);

  // Assemble and solve the HDM at the specified frequency.
  void SolveHDM(ExcitationIdx excitation_idx, double omega, ComplexVector &u);

  // Add field configuration to the reduced-order basis and update the PROM.
  void UpdatePROM(const ComplexVector &u, std::string_view node_label);

  // Add solution u to the minimal-rational interpolation for error estimation. MRI are
  // separated by excitation index.
  void UpdateMRI(ExcitationIdx excitation_idx, double omega, const ComplexVector &u);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional space.
  void SolvePROM(ExcitationIdx excitation_idx, double omega, ComplexVector &u);

  // Compute the location(s) of the maximum error in the range of the previously sampled
  // parameter points.
  std::vector<double> FindMaxError(ExcitationIdx excitation_idx, int N = 1) const
  {
    return mri.at(excitation_idx).FindMaxError(N);
  }

  // Compute eigenvalue estimates for the current PROM system.
  std::vector<std::complex<double>> ComputeEigenvalueEstimates() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
