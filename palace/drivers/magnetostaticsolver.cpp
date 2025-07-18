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
  int n_step = static_cast<int>(curlcurl_op.GetSurfaceCurrentOp().Size());
  MFEM_VERIFY(n_step > 0,
              "No surface current boundaries specified for magnetostatic simulation!");

  // Source term and solution vector storage.
  Vector RHS(Curl.Width()), B(Curl.Height());
  std::vector<Vector> A(n_step);
  std::vector<double> I_inc(n_step);

  // Initialize structures for storing and reducing the results of error estimation.
  CurlFluxErrorEstimator estimator(
      curlcurl_op.GetMaterialOp(), curlcurl_op.GetRTSpace(), curlcurl_op.GetNDSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over current source boundaries.
  Mpi::Print("\nComputing magnetostatic fields for {:d} source {}\n", n_step,
             (n_step > 1) ? "boundaries" : "boundary");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : curlcurl_op.GetSurfaceCurrentOp())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed current on the specified source.
    Mpi::Print("\n");
    A[step].SetSize(RHS.Size());
    A[step].UseDevice(true);
    A[step] = 0.0;
    curlcurl_op.GetExcitationVector(idx, RHS);
    ksp.Mult(RHS, A[step]);

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||A|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(curlcurl_op.GetComm(), A[step]),
               linalg::Norml2(curlcurl_op.GetComm(), RHS));

    // Compute B = ∇ x A on the true dofs.
    Curl.Mult(A[step], B);

    // Save excitation current for inductance matrix calculation.
    I_inc[step] = data.GetExcitationCurrent();

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, A[step], B, idx);

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(B, total_domain_energy, indicator);

    // Next source.
    step++;
  }

  // Postprocess the inductance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(post_op, curlcurl_op.GetSurfaceCurrentOp(), A, I_inc);
  post_op.MeasureFinalize(indicator);
  return {indicator, curlcurl_op.GlobalTrueVSize()};
}

void MagnetostaticSolver::PostprocessTerminals(
    PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
    const SurfaceCurrentOperator &surf_j_op, const std::vector<Vector> &A,
    const std::vector<double> &I_inc) const
{
  // Postprocess the Maxwell inductance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the magnetic field energy based on a current
  // excitation for each port. Alternatively, we could compute the resulting loop fluxes to
  // get M directly as:
  //                         Φ_i = ∫ B ⋅ n_j dS
  // and M_ij = Φ_i/I_j. The energy formulation avoids having to locally integrate B =
  // ∇ x A.
  mfem::DenseMatrix M(A.size()), Mm(A.size());
  for (int i = 0; i < M.Height(); i++)
  {
    // Diagonal: Mᵢᵢ = 2 Uₘ(Aᵢ) / Iᵢ² = (Aᵢᵀ K Aᵢ) / Iᵢ²
    auto &A_gf = post_op.GetAGridFunction().Real();
    auto &H_gf = post_op.GetDomainPostOp().H;
    A_gf.SetFromTrueDofs(A[i]);
    post_op.GetDomainPostOp().M_mag->Mult(A_gf, H_gf);
    M(i, i) = Mm(i, i) =
        linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf) / (I_inc[i] * I_inc[i]);

    // Off-diagonals: Mᵢⱼ = Uₘ(Aᵢ + Aⱼ) / (Iᵢ Iⱼ) - 1/2 (Iᵢ/Iⱼ Mᵢᵢ + Iⱼ/Iᵢ Mⱼⱼ)
    //                    = (Aⱼᵀ K Aᵢ) / (Iᵢ Iⱼ)
    for (int j = i + 1; j < M.Width(); j++)
    {
      A_gf.SetFromTrueDofs(A[j]);
      M(i, j) = linalg::Dot<Vector>(post_op.GetComm(), A_gf, H_gf) / (I_inc[i] * I_inc[j]);
      Mm(i, j) = -M(i, j);
      Mm(i, i) -= Mm(i, j);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      M(i, j) = M(j, i);
      Mm(i, j) = Mm(j, i);
      Mm(i, i) -= Mm(i, j);
    }
  }
  mfem::DenseMatrix Minv(M);
  Minv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using fmt::format;

  // Write inductance matrix data.
  auto PrintMatrix = [&surf_j_op, this](const std::string &file, const std::string &name,
                                        const std::string &unit,
                                        const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    int j = 0;
    for (const auto &[idx2, data2] : surf_j_op)
    {
      output.table.insert(format("i2{}", idx2), format("{}[i][{}] {}", name, idx2, unit));
      // Use the fact that iterator over i and j is the same span.
      output.table["i"] << idx2;

      auto &col = output.table[format("i2{}", idx2)];
      for (std::size_t i = 0; i < surf_j_op.Size(); i++)
      {
        col << mat(i, j) * scale;
      }
      j++;
    }
    output.WriteFullTableTrunc();
  };
  const double H = iodata.units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
  PrintMatrix("terminal-M.csv", "M", "(H)", M, H);
  PrintMatrix("terminal-Minv.csv", "M⁻¹", "(1/H)", Minv, 1.0 / H);
  PrintMatrix("terminal-Mm.csv", "M_m", "(H)", Mm, H);

  // Also write out a file with source current excitations.
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
}

}  // namespace palace
