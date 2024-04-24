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
  LaplaceOperator laplaceop(iodata, mesh);
  auto K = laplaceop.GetStiffnessMatrix();
  const auto &Grad = laplaceop.GetGradMatrix();
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
  Vector RHS(Grad.Width()), E(Grad.Height());
  std::vector<Vector> V(nstep);

  // Initialize structures for storing and reducing the results of error estimation.
  GradFluxErrorEstimator estimator(
      laplaceop.GetMaterialOp(), laplaceop.GetH1Space(), laplaceop.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

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
    ksp.Mult(RHS, V[step]);

    // Compute E = -∇V on the true dofs, and set the internal GridFunctions in PostOperator
    // for all postprocessing operations.
    BlockTimer bt2(Timer::POSTPRO);
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);
    postop.SetVGridFunction(V[step]);
    postop.SetEGridFunction(E);
    double E_elec = postop.GetEFieldEnergy();
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplaceop.GetComm(), V[step]),
               linalg::Norml2(laplaceop.GetComm(), RHS));
    {
      const double J = iodata.DimensionalizeValue(IoData::ValueType::ENERGY, 1.0);
      Mpi::Print(" Field energy E = {:.3e} J\n", E_elec * J);
    }

    // Calculate and record the error indicators.
    Mpi::Print(" Updating solution error estimates\n");
    estimator.AddErrorIndicator(V[step], indicator);

    // Postprocess field solutions and optionally write solution to disk.
    Postprocess(postop, step, idx, E_elec, (step == nstep - 1) ? &indicator : nullptr);

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  PostprocessTerminals(postop, laplaceop.GetSources(), V);
  return {indicator, laplaceop.GlobalTrueVSize()};
}

void ElectrostaticSolver::Postprocess(const PostOperator &postop, int step, int idx,
                                      double E_elec, const ErrorIndicator *indicator) const
{
  // The internal GridFunctions for PostOperator have already been set from the V solution
  // in the main loop.
  PostprocessDomains(postop, "i", step, idx, E_elec, 0.0, 0.0, 0.0);
  PostprocessSurfaces(postop, "i", step, idx, E_elec, 0.0, 1.0, 0.0);
  PostprocessProbes(postop, "i", step, idx);
  if (step < iodata.solver.electrostatic.n_post)
  {
    PostprocessFields(postop, step, idx);
    Mpi::Print(" Wrote fields to disk for terminal {:d}\n", idx);
  }
  if (indicator)
  {
    PostprocessErrorIndicator(postop, *indicator, iodata.solver.electrostatic.n_post > 0);
  }
}

void ElectrostaticSolver::PostprocessTerminals(
    PostOperator &postop, const std::map<int, mfem::Array<int>> &terminal_sources,
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
    // Diagonal: Cᵢᵢ = 2 Uₑ(Vᵢ) / Vᵢ² = (Vᵢᵀ K Vᵢ) / Vᵢ²
    auto &V_gf = postop.GetVGridFunction().Real();
    auto &D_gf = postop.GetDomainPostOp().D;
    V_gf.SetFromTrueDofs(V[i]);
    postop.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    C(i, i) = Cm(i, i) = linalg::Dot<Vector>(postop.GetComm(), V_gf, D_gf);

    // Off-diagonals: Cᵢⱼ = Uₑ(Vᵢ + Vⱼ) / (Vᵢ Vⱼ) - 1/2 (Vᵢ/Vⱼ Cᵢᵢ + Vⱼ/Vᵢ Cⱼⱼ)
    //                    = (Vⱼᵀ K Vᵢ) / (Vᵢ Vⱼ)
    for (int j = i + 1; j < C.Width(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      C(i, j) = linalg::Dot<Vector>(postop.GetComm(), V_gf, D_gf);
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
