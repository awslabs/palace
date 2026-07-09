// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eps.hpp"

#include "linalg/arpack.hpp"
#include "linalg/slepc.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

namespace palace
{

using namespace std::complex_literals;

std::unique_ptr<EigenvalueSolver>
BuildEigenvalueSolver(MPI_Comm comm, const config::EigenSolverData &eigenmode,
                      Orthogonalization gs_orthog, int verbose, bool quadratic,
                      bool nonlinear_slp)
{
  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenvalue solver requires building with ARPACK or SLEPc!"
#endif
  MFEM_VERIFY(!nonlinear_slp || type != EigenSolverBackend::ARPACK,
              "SLP nonlinear eigenvalue solver requires the SLEPc backend!");
  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (quadratic)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(comm, verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(comm, verbose);
    }
#endif
  }
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (nonlinear_slp)
    {
      slepc = std::make_unique<slepc::SlepcNEPSolver>(comm, verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
    }
    else
    {
      if (quadratic)
      {
        if (!eigenmode.pep_linear)
        {
          slepc = std::make_unique<slepc::SlepcPEPSolver>(comm, verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
        }
        else
        {
          slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(comm, verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
        }
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcEPSSolver>(comm, verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    }
    slepc->SetOrthogonalization(gs_orthog == Orthogonalization::MGS,
                                gs_orthog == Orthogonalization::CGS2);
    eigen = std::move(slepc);
#endif
  }
  MFEM_VERIFY(eigen,
              "Failed to configure an eigenvalue solver for the requested backend!");
  eigen->SetNumModes(eigenmode.n, eigenmode.max_size);
  eigen->SetTol(eigenmode.tol);
  eigen->SetMaxIter(eigenmode.max_it);
  return eigen;
}

void SetEigenSolverShiftInvert(EigenvalueSolver &eigen, EigenSolverBackend backend,
                               double target, bool quadratic, bool nonlinear_slp)
{
  if (quadratic || nonlinear_slp)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen.SetShiftInvert(1i * target);
    if (nonlinear_slp)
    {
      eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    }
    else
    {
      eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    eigen.SetShiftInvert(target * target);
    if (backend == EigenSolverBackend::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen.SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }
}

std::complex<double> EigenvalueToOmega(std::complex<double> lambda, bool quadratic)
{
  if (quadratic)
  {
    // Quadratic (and nonlinear) EVP solves for eigenvalue λ = iω.
    return lambda / 1i;
  }
  // Linear EVP has eigenvalue μ = -λ² = ω².
  return std::sqrt(lambda);
}

}  // namespace palace
