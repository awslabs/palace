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
  PostOperator post_op(iodata, laplace_op, "electrostatic");
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0, "No terminal boundaries specified for electrostatic simulation!");
  PostprocessPrintResults post_results(root, post_dir, post_op,
                                       iodata.solver.electrostatic.n_post);

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
  Mpi::Print("\nComputing electrostatic fields for {:d} terminal boundar{}\n", n_step,
             (n_step > 1) ? "ies" : "y");
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

    // Compute E = -∇V on the true dofs, and set the internal GridFunctions in PostOperator
    // for all postprocessing operations.
    BlockTimer bt2(Timer::POSTPRO);
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);
    post_op.SetVGridFunction(V[step]);
    post_op.SetEGridFunction(E);
    const double E_elec = post_op.GetEFieldEnergy();
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), V[step]),
               linalg::Norml2(laplace_op.GetComm(), RHS));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy E = {:.3e} J\n", E_elec * J);
    }

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(E, E_elec, indicator);

    // Postprocess field solutions and optionally write solution to disk.
    post_results.PostprocessStep(iodata, post_op, step, idx, E_elec);

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(post_op, laplace_op.GetSources(), V);
  post_results.PostprocessFinal(post_op, indicator);
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void ElectrostaticSolver::PostprocessTerminals(
    PostOperator &post_op, const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &V) const
{
  // Postprocess the Maxwell capacitance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the electric field energy based on a unit voltage
  // excitation for each terminal. Alternatively, we could compute the resulting terminal
  // charges from the prescribed voltage to get C directly as:
  //         Q_i = ∫ ρ dV = ∫ ∇ ⋅ (ε E) dV = ∫ (ε E) ⋅ n dS
  // and C_ij = Q_i/V_j. The energy formulation avoids having to locally integrate E = -∇V.
  mfem::DenseMatrix C(V.size()), Cm(V.size());
  for (int i = 0; i < C.Height(); i++)
  {
    // Diagonal: Cᵢᵢ = 2 Uₑ(Vᵢ) / Vᵢ² = (Vᵢᵀ K Vᵢ) / Vᵢ² (with ∀i, Vᵢ = 1)
    auto &V_gf = post_op.GetVGridFunction().Real();
    auto &D_gf = post_op.GetDomainPostOp().D;
    V_gf.SetFromTrueDofs(V[i]);
    post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    C(i, i) = Cm(i, i) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);

    // Off-diagonals: Cᵢⱼ = Uₑ(Vᵢ + Vⱼ) / (Vᵢ Vⱼ) - 1/2 (Vᵢ/Vⱼ Cᵢᵢ + Vⱼ/Vᵢ Cⱼⱼ)
    //                    = (Vⱼᵀ K Vᵢ) / (Vᵢ Vⱼ)
    for (int j = i + 1; j < C.Width(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      C(i, j) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);
      Cm(i, j) = -C(i, j);
      Cm(i, i) -= Cm(i, j);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      C(i, j) = C(j, i);
      Cm(i, j) = Cm(j, i);
      Cm(i, i) -= Cm(i, j);
    }
  }
  mfem::DenseMatrix Cinv(C);
  Cinv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Write capactance matrix data.
  auto PrintMatrix = [&terminal_sources, this](const std::string &file,
                                               const std::string &name,
                                               const std::string &unit,
                                               const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert_column(Column("i", "i", 0, {}, {}, ""));
    int j = 0;
    for (const auto &[idx2, data2] : terminal_sources)
    {
      output.table.insert_column(format("i2{}", idx2),
                                 format("{}[i][{}] {}", name, idx2, unit));
      // Use the fact that iterator over i and j is the same span
      output.table["i"] << idx2;

      auto &col = output.table[format("i2{}", idx2)];
      for (int i = 0; i < terminal_sources.size(); i++)
      {
        col << mat(i, j) * scale;
      }
      j++;
    }
    output.WriteFullTableTrunc();
  };
  const double F = iodata.DimensionalizeValue(VT::CAPACITANCE, 1.0);
  PrintMatrix("terminal-C.csv", "C", "(F)", C, F);
  PrintMatrix("terminal-Cinv.csv", "C⁻¹", "(1/F)", Cinv, 1.0 / F);
  PrintMatrix("terminal-Cm.csv", "C_m", "(F)", Cm, F);

  // Also write out a file with terminal voltage excitations.
  {
    TableWithCSVFile terminal_V(post_dir / "terminal-V.csv");
    terminal_V.table.insert_column(Column("i", "i", 0, {}, {}, ""));
    terminal_V.table.insert_column("Vinc", "V_inc[i] (V)");
    int i = 0;
    for (const auto &[idx, data] : terminal_sources)
    {
      terminal_V.table["i"] << double(idx);
      terminal_V.table["Vinc"] << iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
      i++;
    }
    terminal_V.WriteFullTableTrunc();
  }
}

ElectrostaticSolver::PostprocessPrintResults::PostprocessPrintResults(
    bool root, const fs::path &post_dir, const PostOperator &post_op, int n_post_)
  : n_post(n_post_), write_paraview_fields(n_post_ > 0),
    domains{true, root, post_dir, post_op.GetDomainPostOp(), "i", n_post},
    surfaces{true, root, post_dir, post_op, "i", n_post},
    probes{true, root, post_dir, post_op, "i", n_post},
    error_indicator{true, root, post_dir}
{
}

void ElectrostaticSolver::PostprocessPrintResults::PostprocessStep(
    const IoData &iodata, const PostOperator &post_op, int step, int idx, double E_elec)
{
  domains.AddMeasurement(idx, post_op, E_elec, 0.0, 0.0, 0.0, iodata);
  surfaces.AddMeasurement(idx, post_op, E_elec, 0.0, iodata);
  probes.AddMeasurement(idx, post_op, iodata);
  // The internal GridFunctions in PostOperator have already been set from V:
  if (write_paraview_fields && step < n_post)
  {
    Mpi::Print("\n");
    post_op.WriteFields(step, idx);
    Mpi::Print(" Wrote fields to disk for source {:d}\n", idx);
  }
}

void ElectrostaticSolver::PostprocessPrintResults::PostprocessFinal(
    const PostOperator &post_op, const ErrorIndicator &indicator)
{
  BlockTimer bt0(Timer::POSTPRO);
  error_indicator.PrintIndicatorStatistics(post_op, indicator);
  if (write_paraview_fields)
  {
    post_op.WriteFieldsFinal(&indicator);
  }
}

}  // namespace palace
