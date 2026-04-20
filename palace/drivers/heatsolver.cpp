// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "heatsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
HeatSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  LaplaceOperator laplace_op(iodata, mesh);
  auto K = laplace_op.GetStiffnessMatrix();
  const auto &Grad = laplace_op.GetGradMatrix();
  SaveMetadata(laplace_op.GetH1Spaces());

  // Set up the linear solver.
  KspSolver ksp(iodata, laplace_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Terminal indices are the set of boundaries over which to compute the temperature
  // solutions. Terminal boundaries are aliases for prescribed temperature surfaces.
  PostOperator<ProblemType::HEAT> post_op(iodata, laplace_op);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0,
              "No terminal boundaries specified for heat equation simulation!");

  // Right-hand side term and solution vector storage.
  Vector RHS(Grad.Width()), Q(Grad.Height());
  std::vector<Vector> T(n_step);

  // Initialize structures for storing and reducing the results of error estimation.
  GradFluxErrorEstimator estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing heat equation fields for {:d} terminal {}\n", n_step,
             (n_step > 1) ? "boundaries" : "boundary");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero temperature on the
    // specified terminal.
    Mpi::Print("\n");
    laplace_op.GetExcitationVector(idx, *K, T[step], RHS);
    ksp.Mult(RHS, T[step]);

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||T|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), T[step]),
               linalg::Norml2(laplace_op.GetComm(), RHS));

    // Compute Q = -grad(T) on the true dofs.
    Q = 0.0;
    Grad.AddMult(T[step], Q, -1.0);

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, T[step], Q, idx);

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(Q, total_domain_energy, indicator);

    // Next terminal.
    step++;
  }

  // Postprocess the thermal conductance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(post_op, laplace_op.GetSources(), T);
  post_op.MeasureFinalize(indicator);
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void HeatSolver::PostprocessTerminals(
    PostOperator<ProblemType::HEAT> &post_op,
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &T) const
{
  // Postprocess the thermal conductance matrix. Analogous to the capacitance matrix in
  // electrostatics: G_ij = (T_j^T K T_i) / (T_i * T_j), where T_i is the temperature
  // field for unit temperature on terminal i.
  mfem::DenseMatrix G(T.size()), Gm(T.size());
  for (int i = 0; i < G.Height(); i++)
  {
    auto &T_gf = post_op.GetVGridFunction().Real();
    auto &D_gf = post_op.GetDomainPostOp().D;
    T_gf.SetFromTrueDofs(T[i]);
    post_op.GetDomainPostOp().M_elec->Mult(T_gf, D_gf);
    G(i, i) = Gm(i, i) = linalg::Dot<Vector>(post_op.GetComm(), T_gf, D_gf);

    for (int j = i + 1; j < G.Width(); j++)
    {
      T_gf.SetFromTrueDofs(T[j]);
      G(i, j) = linalg::Dot<Vector>(post_op.GetComm(), T_gf, D_gf);
      Gm(i, j) = -G(i, j);
      Gm(i, i) -= Gm(i, j);
    }

    for (int j = 0; j < i; j++)
    {
      G(i, j) = G(j, i);
      Gm(i, j) = Gm(j, i);
      Gm(i, i) -= Gm(i, j);
    }
  }
  mfem::DenseMatrix Ginv(G);
  Ginv.Invert();

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using VT = Units::ValueType;

  // Write conductance matrix data.
  auto PrintMatrix = [&terminal_sources, this](const std::string &file,
                                               const std::string &name,
                                               const std::string &unit,
                                               const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    int j = 0;
    for (const auto &[idx2, data2] : terminal_sources)
    {
      output.table.insert(fmt::format("i2{}", idx2),
                          fmt::format("{}[i][{}] {}", name, idx2, unit));
      output.table["i"] << idx2;

      auto &col = output.table[fmt::format("i2{}", idx2)];
      for (std::size_t i = 0; i < terminal_sources.size(); i++)
      {
        col << mat(i, j) * scale;
      }
      j++;
    }
    output.WriteFullTableTrunc();
  };
  const double F = iodata.units.Dimensionalize<VT::CAPACITANCE>(1.0);
  PrintMatrix("terminal-G.csv", "G", "(F)", G, F);
  PrintMatrix("terminal-Ginv.csv", "G⁻¹", "(1/F)", Ginv, 1.0 / F);
  PrintMatrix("terminal-Gm.csv", "G_m", "(F)", Gm, F);

  // Also write out a file with terminal temperature excitations.
  {
    TableWithCSVFile terminal_T(post_dir / "terminal-T.csv");
    terminal_T.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_T.table.insert("Tinc", "T_inc[i] (V)");
    for (const auto &[idx, data] : terminal_sources)
    {
      terminal_T.table["i"] << double(idx);
      terminal_T.table["Tinc"] << iodata.units.Dimensionalize<VT::VOLTAGE>(1.0);
    }
    terminal_T.WriteFullTableTrunc();
  }
}

}  // namespace palace
