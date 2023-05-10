// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#if 0  // XX TODO DISABLE ROM FOR NOW

#include <memory>
#include <random>
#include <vector>
#include <mfem.hpp>
#include "linalg/curlcurl.hpp"
#include "linalg/ksp.hpp"
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

//
// A class handling PROM construction and use for adaptive fast frequency sweeps.
//
class RomOperator
{
private:
  // Reference to HDM discretization (not owned).
  SpaceOperator &spaceop;

  // HDM system matrices and excitation RHS.
  std::unique_ptr<petsc::PetscParMatrix> K, M, C;
  std::unique_ptr<petsc::PetscParVector> RHS1;

  // HDM storage for terms with non-polynomial frequency dependence.
  std::vector<std::unique_ptr<petsc::PetscParMatrix>> A2;
  std::vector<std::unique_ptr<petsc::PetscParVector>> RHS2;
  bool init2, hasA2, hasRHS2;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<KspSolver> ksp0;
  std::unique_ptr<KspPreconditioner> pc0;

  // Working storage for HDM vectors.
  std::unique_ptr<petsc::PetscParVector> E0, R0, T0;

  // PROM matrices, vectors, and linear solver.
  std::unique_ptr<petsc::PetscDenseMatrix> Kr, Mr, Cr, Ar;
  std::unique_ptr<petsc::PetscParVector> RHS1r, RHSr, Er;
  std::unique_ptr<KspSolver> ksp;

  // Linear solver for inner product solves for error metric.
  std::unique_ptr<CurlCurlMassSolver> kspKM;
  std::unique_ptr<petsc::PetscParMatrix> opKM;

  // PROM reduced-order basis and parameter domain samplings.
  int dim;
  std::unique_ptr<petsc::PetscDenseMatrix> V;
  std::vector<double> Ps, PmPs;
  double omega_min, delta_omega;
  std::default_random_engine engine;

  // Compute the error metric for the PROM solution (computed internally) at the specified
  // frequency.
  double ComputeError(double omega);

  // Helper functions for reduced-order matrix or vector construction/update.
  void BVMatProjectInternal(petsc::PetscDenseMatrix &V, petsc::PetscParMatrix &A,
                            petsc::PetscDenseMatrix &Ar, petsc::PetscParVector &r, int n0,
                            int n);
  void BVDotVecInternal(petsc::PetscDenseMatrix &V, petsc::PetscParVector &b,
                        petsc::PetscParVector &br, int n0, int n);

public:
  RomOperator(const IoData &iodata, SpaceOperator &sp, int nmax);

  // Return set of sampled parameter points for basis construction.
  const std::vector<double> &GetSampleFrequencies() const { return Ps; }

  // Return PROM dimension.
  int GetReducedDimension() const { return dim; }

  // Return number of HDM linear solves and linear solver iterations performed during
  // offline training.
  int GetTotalKspMult() const { return ksp0->GetTotalNumMult(); }
  int GetTotalKspIter() const { return ksp0->GetTotalNumIter(); }

  // Initialize the solution basis with HDM samples at the minimum and maximum frequencies.
  void Initialize(int steps, double start, double delta);

  // Assemble and solve the HDM at the specified frequency, adding the solution vector to
  // the reduced-order basis.
  void SolveHDM(double omega, petsc::PetscParVector &E, bool print = false);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional solution space.
  void AssemblePROM(double omega);
  void SolvePROM(petsc::PetscParVector &E);

  // Compute the maximum error over a randomly sampled set of candidate points. Returns the
  // maximum error and its correcponding frequency, as well as the number of candidate
  // points used (if fewer than those availble in the unsampled parameter domain).
  double ComputeMaxError(int Nc, double &omega_star);
};

}  // namespace palace

#endif

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
