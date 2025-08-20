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
  // Determine which excitations are flux loops vs current sources
  int n_current = static_cast<int>(surf_j_op.Size());
  int n_flux = A.size() - n_current;
  std::vector<bool> is_flux_loop(A.size(), false);
  
  // Mark flux loop excitations (they come after current excitations)
  for (int i = n_current; i < static_cast<int>(A.size()); i++)
  {
    is_flux_loop[i] = true;
  }

  mfem::DenseMatrix M(A.size()), Mm(A.size()), R(A.size()), Rm(A.size());
  
  for (int i = 0; i < M.Height(); i++)
  {
    auto &A_gf = post_op.GetAGridFunction().Real();
    auto &H_gf = post_op.GetDomainPostOp().H;
    A_gf.SetFromTrueDofs(A[i]);
    post_op.GetDomainPostOp().M_mag->Mult(A_gf, H_gf);
    double energy_ii = linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);
    
    // Diagonal terms
    if (is_flux_loop[i])
    {
      // Flux-based: R_ii = 2 U_m(A_i) / Φ_i²
      R(i, i) = Rm(i, i) = energy_ii / (Phi_inc[i] * Phi_inc[i]);
    }
    else
    {
      // Current-based: M_ii = 2 U_m(A_i) / I_i²
      M(i, i) = Mm(i, i) = energy_ii / (I_inc[i] * I_inc[i]);
    }
    
    // Off-diagonal terms
    for (int j = i + 1; j < M.Width(); j++)
    {
      A_gf.SetFromTrueDofs(A[j]);
      double cross_energy = linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf);
      
      if (is_flux_loop[i] && is_flux_loop[j])
      {
        // Flux-flux: R_ij = (A_j^T K A_i) / (Φ_i Φ_j)
        R(i, j) = cross_energy / (Phi_inc[i] * Phi_inc[j]);
        Rm(i, j) = -R(i, j);
        Rm(i, i) -= Rm(i, j);
      }
      else if (!is_flux_loop[i] && !is_flux_loop[j])
      {
        // Current-current: M_ij = (A_j^T K A_i) / (I_i I_j)
        M(i, j) = cross_energy / (I_inc[i] * I_inc[j]);
        Mm(i, j) = -M(i, j);
        Mm(i, i) -= Mm(i, j);
      }
      // Mixed terms: skip for now
    }
    
    // Copy to lower triangle
    for (int j = 0; j < i; j++)
    {
      if (is_flux_loop[i] && is_flux_loop[j])
      {
        R(i, j) = R(j, i);
        Rm(i, j) = Rm(j, i);
        Rm(i, i) -= Rm(i, j);
      }
      else if (!is_flux_loop[i] && !is_flux_loop[j])
      {
        M(i, j) = M(j, i);
        Mm(i, j) = Mm(j, i);
        Mm(i, i) -= Mm(i, j);
      }
    }
  }
  
  // For flux loops: compute inductance from reluctance
  if (n_flux > 0)
  {
    mfem::DenseMatrix R_flux(n_flux), M_flux(n_flux);
    for (int i = 0; i < n_flux; i++)
    {
      for (int j = 0; j < n_flux; j++)
      {
        R_flux(i, j) = R(n_current + i, n_current + j);
      }
    }
    M_flux = R_flux;
    M_flux.Invert(); // M_flux = R_flux^(-1)
    
    // Copy back to main matrix
    for (int i = 0; i < n_flux; i++)
    {
      for (int j = 0; j < n_flux; j++)
      {
        M(n_current + i, n_current + j) = M_flux(i, j);
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
  
  if (n_flux > 0)
  {
    PrintMatrix("terminal-R.csv", "R", "(1/H)", R, 1.0 / H);
    PrintMatrix("terminal-Rm.csv", "R_m", "(1/H)", Rm, 1.0 / H);
  }

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
