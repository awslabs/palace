// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

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
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
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

  // Terminal indices are the set of boundaries over which to compute the capacitance
  // matrix. Terminal boundaries are aliases for ports.
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0, "No terminal boundaries specified for electrostatic simulation!");

  // Right-hand side term and solution vector storage.
  Vector RHS(Grad.Width()), E(Grad.Height());
  std::vector<Vector> V(n_step);

  // Initialize structures for storing and reducing the results of error estimation.
  GradFluxErrorEstimator estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing electrostatic fields for {:d} terminal {}\n", n_step,
             (n_step > 1) ? "boundaries" : "boundary");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero voltage on the specified
    // terminal.
    Mpi::Print("\n");
    laplace_op.GetExcitationVector(idx, *K, V[step], RHS);
    ksp.Mult(RHS, V[step]);

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), V[step]),
               linalg::Norml2(laplace_op.GetComm(), RHS));

    // Compute E = -∇V on the true dofs.
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, V[step], E, idx);

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, total_domain_energy, indicator);

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(laplace_op, *K, laplace_op.GetSources(), V);
  post_op.MeasureFinalize(indicator);
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void ElectrostaticSolver::PostprocessTerminals(
    const LaplaceOperator &laplace_op, const Operator &K,
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &V) const
{
  // Postprocess the Maxwell capacitance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the electric field energy based on a unit voltage
  // excitation for each terminal. Alternatively, we could compute the resulting terminal
  // charges from the prescribed voltage to get C directly as:
  //         Q_i = ∫ ρ dV = ∫ ∇ ⋅ (ε E) dV = ∫ (ε E) ⋅ n dS
  // and C_ij = Q_i/V_j. The energy formulation avoids having to locally integrate E = -∇V.
  auto energy_op = LaplaceOperator::GetUnconstrainedStiffnessOperator(K);
  mfem::DenseMatrix C =
      LaplaceOperator::ComputeCapacitanceMatrix(laplace_op.GetComm(), *energy_op, V);
  mfem::DenseMatrix Cm(C);
  for (int i = 0; i < C.Height(); i++)
  {
    for (int j = 0; j < C.Width(); j++)
    {
      if (i != j)
      {
        Cm(i, j) = -C(i, j);
        Cm(i, i) += C(i, j);
      }
    }
  }
  mfem::DenseMatrix Cinv(C);
  Cinv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  std::vector<int> terminal_indices;
  terminal_indices.reserve(terminal_sources.size());
  for (const auto &[idx, data] : terminal_sources)
  {
    terminal_indices.push_back(idx);
  }
  using VT = Units::ValueType;

  // Write capacitance matrix data.
  const double F = iodata.units.Dimensionalize<VT::CAPACITANCE>(1.0);
  LaplaceOperator::WriteTerminalMatrix(laplace_op.GetComm(), post_dir, "terminal-C.csv",
                                       "C", "(F)", terminal_indices, C, F);
  LaplaceOperator::WriteTerminalMatrix(laplace_op.GetComm(), post_dir, "terminal-Cinv.csv",
                                       "C⁻¹", "(1/F)", terminal_indices, Cinv, 1.0 / F);
  LaplaceOperator::WriteTerminalMatrix(laplace_op.GetComm(), post_dir, "terminal-Cm.csv",
                                       "C_m", "(F)", terminal_indices, Cm, F);

  // Also write out a file with terminal voltage excitations.
  if (!root)
  {
    return;
  }
  {
    TableWithCSVFile terminal_V(post_dir / "terminal-V.csv");
    terminal_V.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_V.table.insert("Vinc", "V_inc[i] (V)");
    for (const auto &[idx, data] : terminal_sources)
    {
      terminal_V.table["i"] << double(idx);
      terminal_V.table["Vinc"] << iodata.units.Dimensionalize<VT::VOLTAGE>(1.0);
    }
    terminal_V.WriteFullTableTrunc();
  }
}

}  // namespace palace
