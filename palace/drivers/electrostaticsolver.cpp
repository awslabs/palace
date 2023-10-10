// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
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

ErrorIndicator
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  LaplaceOperator laplaceop(iodata, mesh);
  auto K = laplaceop.GetStiffnessMatrix();
  SaveMetadata(laplaceop.GetH1Spaces());

  // Set up the linear solver.
  KspSolver ksp(iodata, laplaceop.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Terminal indices are the set of boundaries over which to compute the capacitance
  // matrix. Terminal boundaries are aliases for ports.
  PostOperator postop(iodata, laplaceop, "electrostatic");
  int nstep = static_cast<int>(laplaceop.GetSources().size());
  MFEM_VERIFY(nstep > 0, "No terminal boundaries specified for electrostatic simulation!");

  // Right-hand side term and solution vector storage.
  Vector RHS(K->Height());
  std::vector<Vector> V(nstep);

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing electrostatic fields for {:d} terminal boundar{}\n", nstep,
             (nstep > 1) ? "ies" : "y");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplaceop.GetSources())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, nstep,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero voltage on the specified
    // terminal.
    Mpi::Print("\n");
    laplaceop.GetExcitationVector(idx, *K, V[step], RHS);

    BlockTimer bt1(Timer::SOLVE);
    ksp.Mult(RHS, V[step]);

    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplaceop.GetComm(), V[step]),
               linalg::Norml2(laplaceop.GetComm(), RHS));

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  return Postprocess(laplaceop, postop, V);
}

ErrorIndicator ElectrostaticSolver::Postprocess(LaplaceOperator &laplaceop,
                                                PostOperator &postop,
                                                const std::vector<Vector> &V) const
{
  // Postprocess the Maxwell capacitance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the electric field energy based on a unit voltage
  // excitation for each terminal. Alternatively, we could compute the resulting terminal
  // charges from the prescribed voltage to get C directly as:
  //         Q_i = ∫ ρ dV = ∫ ∇ ⋅ (ε E) dV = ∫ (ε E) ⋅ n dS
  // and C_ij = Q_i/V_j. The energy formulation avoids having to locally integrate E = -∇V.
  // Additionally compute error estimates for each terminal.
  auto Grad = laplaceop.GetGradMatrix();
  const std::map<int, mfem::Array<int>> &terminal_sources = laplaceop.GetSources();
  int nstep = static_cast<int>(terminal_sources.size());
  mfem::DenseMatrix C(nstep), Cm(nstep);
  Vector E(Grad->Height()), Vij(Grad->Width());
  if (iodata.solver.electrostatic.n_post > 0)
  {
    Mpi::Print("\n");
  }

  // Calculate and record the error indicators.
  GradFluxErrorEstimator estimator(
      laplaceop.GetMaterialOp(), laplaceop.GetH1Spaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it,
      iodata.problem.verbose, iodata.solver.pa_order_threshold);
  ErrorIndicator indicator;
  for (int i = 0; i < nstep; i++)
  {
    estimator.AddErrorIndicator(V[i], indicator);
  }

  int i = 0;
  for (const auto &[idx, data] : terminal_sources)
  {
    // Compute E = -∇V on the true dofs, and set the internal GridFunctions in PostOperator
    // for all postprocessing operations.
    E = 0.0;
    Grad->AddMult(V[i], E, -1.0);
    postop.SetEGridFunction(E);
    postop.SetVGridFunction(V[i]);
    double Ue = postop.GetEFieldEnergy();
    PostprocessDomains(postop, "i", i, idx, Ue, 0.0, 0.0, 0.0);
    PostprocessSurfaces(postop, "i", i, idx, Ue, 0.0, 1.0, 0.0);
    PostprocessProbes(postop, "i", i, idx);
    if (i < iodata.solver.electrostatic.n_post)
    {
      PostprocessFields(postop, i, idx, (i == 0) ? &indicator : nullptr);
      Mpi::Print(" Wrote fields to disk for terminal {:d}\n", idx);
    }
    if (i == 0)
    {
      PostprocessErrorIndicator(postop, indicator);
    }

    // Diagonal: C_ii = 2 U_e(V_i) / V_i².
    C(i, i) = Cm(i, i) = 2.0 * Ue;
    i++;
  }

  // Off-diagonals: C_ij = U_e(V_i + V_j) / (V_i V_j) - 1/2 (V_i/V_j C_ii + V_j/V_i C_jj).
  for (i = 0; i < C.Height(); i++)
  {
    for (int j = 0; j < C.Width(); j++)
    {
      if (j < i)
      {
        // Copy lower triangle from already computed upper triangle.
        C(i, j) = C(j, i);
        Cm(i, j) = Cm(j, i);
        Cm(i, i) -= Cm(i, j);
      }
      else if (j > i)
      {
        linalg::AXPBYPCZ(1.0, V[i], 1.0, V[j], 0.0, Vij);
        E = 0.0;
        Grad->AddMult(Vij, E, -1.0);
        postop.SetEGridFunction(E);
        double Ue = postop.GetEFieldEnergy();
        C(i, j) = Ue - 0.5 * (C(i, i) + C(j, j));
        Cm(i, j) = -C(i, j);
        Cm(i, i) -= Cm(i, j);
      }
    }
  }
  mfem::DenseMatrix Cinv(C);
  Cinv.Invert();  // In-place, uses LAPACK (when available) and should be cheap
  PostprocessTerminals(terminal_sources, C, Cinv, Cm);
  return indicator;
}

void ElectrostaticSolver::PostprocessTerminals(
    const std::map<int, mfem::Array<int>> &terminal_sources, const mfem::DenseMatrix &C,
    const mfem::DenseMatrix &Cinv, const mfem::DenseMatrix &Cm) const
{
  // Only root writes to disk (every process has full matrices).
  if (!root || post_dir.length() == 0)
  {
    return;
  }

  // Write capactance matrix data.
  auto PrintMatrix = [&terminal_sources, this](const std::string &file,
                                               const std::string &name,
                                               const std::string &unit,
                                               const mfem::DenseMatrix &mat, double scale)
  {
    std::string path = post_dir + file;
    auto output = OutputFile(path, false);
    output.print("{:>{}s},", "i", table.w1);
    for (const auto &[idx2, data2] : terminal_sources)
    {
      // clang-format off
      output.print("{:>{}s}{}",
                   name + "[i][" + std::to_string(idx2) + "] " + unit, table.w,
                   (idx2 == terminal_sources.rbegin()->first) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
    int i = 0;
    for (const auto &[idx, data] : terminal_sources)
    {
      int j = 0;
      output.print("{:{}.{}e},", static_cast<double>(idx), table.w1, table.p1);
      for (const auto &[idx2, data2] : terminal_sources)
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     mat(i, j) * scale, table.w, table.p,
                     (idx2 == terminal_sources.rbegin()->first) ? "" : ",");
        // clang-format on
        j++;
      }
      output.print("\n");
      i++;
    }
  };
  const double F = iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, 1.0);
  PrintMatrix("terminal-C.csv", "C", "(F)", C, F);
  PrintMatrix("terminal-Cinv.csv", "C⁻¹", "(1/F)", Cinv, 1.0 / F);
  PrintMatrix("terminal-Cm.csv", "C_m", "(F)", Cm, F);
}

}  // namespace palace
