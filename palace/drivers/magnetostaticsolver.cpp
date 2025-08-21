// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "magnetostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"
#include "drivers/surfacecurlsolver.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
MagnetostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs.
  BlockTimer bt0(Timer::CONSTRUCT);
  CurlCurlOperator curlcurl_op(iodata, mesh);
  auto K = curlcurl_op.GetStiffnessMatrix();
  const auto &Curl = curlcurl_op.GetCurlMatrix();
  SaveMetadata(curlcurl_op.GetNDSpaces());

  // Set up the linear solver.
  KspSolver ksp(iodata, curlcurl_op.GetNDSpaces(), &curlcurl_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Terminal indices are the set of boundaries over which to compute the inductance matrix.
  PostOperator<ProblemType::MAGNETOSTATIC> post_op(iodata, curlcurl_op);
  int n_current_steps = static_cast<int>(curlcurl_op.GetSurfaceCurrentOp().Size());
  int n_flux_steps = static_cast<int>(iodata.boundaries.fluxloop.size());
  int n_step = n_current_steps + n_flux_steps;
  
  MFEM_VERIFY(n_step > 0,
              "No surface current boundaries or flux loops specified for magnetostatic simulation!");

  // Source term and solution vector storage.
  Vector RHS(Curl.Width()), B(Curl.Height());
  std::vector<Vector> A(n_step);
  std::vector<double> I_inc(n_step);
  std::vector<double> Phi_inc(n_step);

  // Initialize structures for storing and reducing the results of error estimation.
  CurlFluxErrorEstimator estimator(
      curlcurl_op.GetMaterialOp(), curlcurl_op.GetRTSpace(), curlcurl_op.GetNDSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Unified loop over all excitation sources (current + flux loops)
  Mpi::Print("\nComputing magnetostatic fields for {:d} source {}\n", n_step,
             (n_step > 1) ? "boundaries" : "boundary");
  auto t0 = Timer::Now();
  
  for (int step = 0; step < n_step; step++)
  {
    A[step].SetSize(RHS.Size());
    A[step].UseDevice(true);
    A[step] = 0.0;
    
    if (step < n_current_steps)
    {
      // Current source excitation
      const auto &[idx, data] = *std::next(curlcurl_op.GetSurfaceCurrentOp().begin(), step);
      Mpi::Print("\nIt {:d}/{:d}: Current Index = {:d} (elapsed time = {:.2e} s)\n", 
                step + 1, n_step, idx, Timer::Duration(Timer::Now() - t0).count());
      
      curlcurl_op.GetExcitationVector(idx, RHS);
      I_inc[step] = data.GetExcitationCurrent();
      Phi_inc[step] = 0.0; // Zero flux for current sources
    }
    else
    {
      // Flux loop excitation
      int flux_idx = step - n_current_steps;
      const auto &[idx, flux_data] = *std::next(iodata.boundaries.fluxloop.begin(), flux_idx);
      Mpi::Print("\nIt {:d}/{:d}: FluxLoop Index = {:d} (elapsed time = {:.2e} s)\n", 
                step + 1, n_step, idx, Timer::Duration(Timer::Now() - t0).count());
      
      // Solve 2D surface curl problem for this specific flux loop
      auto flux_solution = curlcurl_op.SolveSurfaceCurlProblem(idx);
      curlcurl_op.GetFluxLoopExcitationVector(*flux_solution, RHS);
      I_inc[step] = 0.0; // Zero current for flux loops
      Phi_inc[step] = flux_data.flux_amounts[0]; // Store prescribed flux
    }
    
    // Solve 3D magnetostatic problem
    ksp.Mult(RHS, A[step]);
    
    // Common post-processing
    Curl.Mult(A[step], B);
    
    // Flux verification for flux loops
    if (step >= n_current_steps)
    {
      int flux_idx = step - n_current_steps;
      const auto &[idx, flux_data] = *std::next(iodata.boundaries.fluxloop.begin(), flux_idx);
      auto &B_gf = post_op.GetBGridFunction().Real();
      B_gf.SetFromTrueDofs(B);
      VerifyFluxThroughHoles(B_gf, flux_data.hole_attributes, flux_data.flux_amounts, 
                            curlcurl_op.GetMesh(), curlcurl_op.GetComm());
    }
    
    // Energy calculation and error estimation
    int terminal_idx = (step < n_current_steps) ? 
      std::next(curlcurl_op.GetSurfaceCurrentOp().begin(), step)->first :
      std::next(iodata.boundaries.fluxloop.begin(), step - n_current_steps)->first;
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, A[step], B, terminal_idx);
    estimator.AddErrorIndicator(B, total_domain_energy, indicator);
  }

  // Postprocess the inductance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(post_op, curlcurl_op.GetSurfaceCurrentOp(), A, I_inc, Phi_inc);
  post_op.MeasureFinalize(indicator);
  return {indicator, curlcurl_op.GlobalTrueVSize()};
}

void MagnetostaticSolver::PostprocessTerminals(
    PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
    const SurfaceCurrentOperator &surf_j_op, 
    const std::vector<Vector> &A,
    const std::vector<double> &I_inc,
    const std::vector<double> &Phi_inc) const
{
  int n_current = static_cast<int>(surf_j_op.Size());
  int n_flux = A.size() - n_current;
  int n = A.size();
  std::vector<bool> is_flux_loop(n, false);
  
  for (int i = n_current; i < n; i++)
    is_flux_loop[i] = true;

  mfem::DenseMatrix M(n), Mm(n);
  mfem::DenseMatrix R(n); // Temporary for computation only
  
  // Compute cross-energy matrix and diagonals
  mfem::DenseMatrix cross_energy(n);
  for (int i = 0; i < n; i++)
  {
    auto &A_gf = post_op.GetAGridFunction().Real();
    auto &H_gf = post_op.GetDomainPostOp().H;
    A_gf.SetFromTrueDofs(A[i]);
    post_op.GetDomainPostOp().M_mag->Mult(A_gf, H_gf);
    cross_energy(i, i) = linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);
    
    // Diagonal terms from energy
    if (is_flux_loop[i])
      R(i, i) = cross_energy(i, i) / (Phi_inc[i] * Phi_inc[i]);
    else
      M(i, i) = cross_energy(i, i) / (I_inc[i] * I_inc[i]);
    
    // Off-diagonal cross-energies
    for (int j = i + 1; j < n; j++)
    {
      A_gf.SetFromTrueDofs(A[j]);
      cross_energy(i, j) = cross_energy(j, i) = 
        linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);
    }
  }
  
  // System of equations for mixed terms
  if (n_current > 0 && n_flux > 0)
  {
    int n_off = n * (n - 1) / 2; // Number of off-diagonal elements
    mfem::DenseMatrix A_sys(n * n, 2 * n_off);
    mfem::Vector b_sys(n * n), x_sol(2 * n_off);
    A_sys = 0.0; b_sys = 0.0;
    
    // Set up M×R = I constraint equations
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        int eq = i * n + j;
        b_sys[eq] = (i == j) ? 1.0 : 0.0;
        
        for (int k = 0; k < n; k++)
        {
          if (i != k && k != j) // Off-diagonal terms
          {
            int M_idx = (i < k) ? i * n + k - (i + 1) * (i + 2) / 2 : 
                                  k * n + i - (k + 1) * (k + 2) / 2;
            int R_idx = (k < j) ? k * n + j - (k + 1) * (k + 2) / 2 : 
                                  j * n + k - (j + 1) * (j + 2) / 2;
            A_sys(eq, M_idx) += (k == j) ? 1.0 : 0.0;
            A_sys(eq, n_off + R_idx) += (i == k) ? 1.0 : 0.0;
          }
        }
        // Diagonal contributions
        if (i == j)
          b_sys[eq] -= M(i, i) * R(j, j);
      }
    }
    
    // Solve system
    mfem::DenseMatrixInverse A_inv(A_sys);
    A_inv.Mult(b_sys, x_sol);
    
    // Extract off-diagonal elements
    int idx = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        M(i, j) = M(j, i) = x_sol[idx];
        R(i, j) = R(j, i) = x_sol[n_off + idx];
        idx++;
      }
    }
  }
  else
  {
    // Pure current or pure flux - use simplified approach
    for (int i = 0; i < n; i++)
    {
      for (int j = i + 1; j < n; j++)
      {
        if (is_flux_loop[i] && is_flux_loop[j])
        {
          R(i, j) = R(j, i) = cross_energy(i, j) / (Phi_inc[i] * Phi_inc[j]);
        }
        else if (!is_flux_loop[i] && !is_flux_loop[j])
        {
          M(i, j) = M(j, i) = cross_energy(i, j) / (I_inc[i] * I_inc[j]);
        }
      }
    }
    
    // Convert flux reluctance to inductance
    if (n_flux > 0)
    {
      mfem::DenseMatrix R_flux(n_flux), M_flux(n_flux);
      for (int i = 0; i < n_flux; i++)
        for (int j = 0; j < n_flux; j++)
          R_flux(i, j) = R(n_current + i, n_current + j);
      
      M_flux = R_flux;
      M_flux.Invert();
      
      for (int i = 0; i < n_flux; i++)
        for (int j = 0; j < n_flux; j++)
          M(n_current + i, n_current + j) = M_flux(i, j);
    }
  }
  
  // Compute Mm matrix only
  for (int i = 0; i < n; i++)
  {
    Mm(i, i) = M(i, i);
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        Mm(i, j) = -M(i, j);
        Mm(i, i) += M(i, j);
      }
    }
  }

  mfem::DenseMatrix Minv(M);
  Minv.Invert();
  
  // Only root writes to disk
  if (!root)
  {
    return;
  }
  using fmt::format;

  // Write matrix data using existing pattern
  auto PrintMatrix = [&surf_j_op, this, n_current, n_flux](const std::string &file, const std::string &name,
                                  const std::string &unit,
                                  const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    
    auto AddTerminal = [&](int idx, int col_pos) {
      output.table.insert(format("i2{}", idx), format("{}[i][{}] {}", name, idx, unit));
      output.table["i"] << idx;
      auto &col = output.table[format("i2{}", idx)];
      for (int i = 0; i < mat.Height(); i++)
      {
        col << mat(i, col_pos) * scale;
      }
    };
    
    // Add current sources
    int j = 0;
    for (const auto &[idx, data] : surf_j_op)
      AddTerminal(idx, j++);
    
    // Add flux loops  
    for (const auto &[idx, flux_data] : iodata.boundaries.fluxloop)
      AddTerminal(idx, j++);
    
    output.WriteFullTableTrunc();
  };

  const double H = iodata.units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  PrintMatrix("terminal-M.csv", "M", "(H)", M, H);
  PrintMatrix("terminal-Minv.csv", "M⁻¹", "(1/H)", Minv, 1.0 / H);
  PrintMatrix("terminal-Mm.csv", "M_m", "(H)", Mm, H);

  // Write out a file with source current excitations.
  {
    TableWithCSVFile terminal_I(post_dir / "terminal-I.csv");
    terminal_I.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_I.table.insert("Iinc", "I_inc[i] (A)");
    int i = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      terminal_I.table["i"] << double(idx);
      terminal_I.table["Iinc"] << iodata.units.Dimensionalize<Units::ValueType::CURRENT>(
          I_inc[i]);
      i++;
    }
    terminal_I.WriteFullTableTrunc();
  }

  // Write out a file with flux loop excitations.
  if (n_flux > 0)
  {
    TableWithCSVFile terminal_Phi(post_dir / "terminal-Phi.csv");
    terminal_Phi.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_Phi.table.insert("Phiinc", "Phi_inc[i] (unit flux)");
    int i = n_current;
    for (const auto &[idx, data] : iodata.boundaries.fluxloop)
    {
      terminal_Phi.table["i"] << double(idx);
      terminal_Phi.table["Phiinc"] << Phi_inc[i];  // Keep as dimensionless
      i++;
    }
    terminal_Phi.WriteFullTableTrunc();
  }
}

}  // namespace palace
