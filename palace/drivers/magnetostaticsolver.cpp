// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "magnetostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
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

ErrorIndicator
MagnetostaticSolver::Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs.
  BlockTimer bt0(Timer::CONSTRUCT);
  CurlCurlOperator curlcurlop(iodata, mesh);
  auto K = curlcurlop.GetStiffnessMatrix();
  SaveMetadata(curlcurlop.GetNDSpaces());

  // Set up the linear solver.
  KspSolver ksp(iodata, curlcurlop.GetNDSpaces(), &curlcurlop.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Terminal indices are the set of boundaries over which to compute the inductance matrix.
  PostOperator postop(iodata, curlcurlop, "magnetostatic");
  int nstep = static_cast<int>(curlcurlop.GetSurfaceCurrentOp().Size());
  MFEM_VERIFY(nstep > 0,
              "No surface current boundaries specified for magnetostatic simulation!");

  // Source term and solution vector storage.
  Vector RHS(K->Height());
  std::vector<Vector> A(nstep);

  // Main loop over current source boundaries.
  Mpi::Print("\nComputing magnetostatic fields for {:d} source boundar{}\n", nstep,
             (nstep > 1) ? "ies" : "y");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : curlcurlop.GetSurfaceCurrentOp())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, nstep,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed current on the specified source.
    Mpi::Print("\n");
    A[step].SetSize(RHS.Size());
    A[step] = 0.0;
    curlcurlop.GetExcitationVector(idx, RHS);

    BlockTimer bt1(Timer::SOLVE);
    ksp.Mult(RHS, A[step]);

    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||A|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(curlcurlop.GetComm(), A[step]),
               linalg::Norml2(curlcurlop.GetComm(), RHS));

    // Next source.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  return Postprocess(curlcurlop, postop, A);
}

ErrorIndicator MagnetostaticSolver::Postprocess(CurlCurlOperator &curlcurlop,
                                                PostOperator &postop,
                                                const std::vector<Vector> &A) const
{
  // Postprocess the Maxwell inductance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the magnetic field energy based on a current
  // excitation for each port. Alternatively, we could compute the resulting loop fluxes to
  // get M directly as:
  //                         Φ_i = ∫ B ⋅ n_j dS
  // and M_ij = Φ_i/I_j. The energy formulation avoids having to locally integrate B =
  // ∇ x A.
  const auto &Curl = curlcurlop.GetCurlMatrix();
  const SurfaceCurrentOperator &surf_j_op = curlcurlop.GetSurfaceCurrentOp();
  int nstep = static_cast<int>(surf_j_op.Size());
  mfem::DenseMatrix M(nstep), Mm(nstep);
  Vector B(Curl.Height()), Aij(Curl.Width());
  Vector Iinc(nstep);
  if (iodata.solver.magnetostatic.n_post > 0)
  {
    Mpi::Print("\n");
  }

  // Calculate and record the error indicators.
  CurlFluxErrorEstimator<Vector> estimator(
      curlcurlop.GetMaterialOp(), curlcurlop.GetNDSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.pa_order_threshold);
  ErrorIndicator indicator;
  for (int i = 0; i < nstep; i++)
  {
    estimator.AddErrorIndicator(A[i], indicator);
  }

  int i = 0;
  for (const auto &[idx, data] : surf_j_op)
  {
    // Get the magnitude of the current excitations (unit J_s,inc, but circuit current I is
    // the integral of J_s,inc over port).
    Iinc(i) = data.GetExcitationCurrent();
    MFEM_VERIFY(Iinc(i) > 0.0, "Zero current excitation for magnetostatic solver!");

    // Compute B = ∇ x A on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    Curl.Mult(A[i], B);
    postop.SetBGridFunction(B);
    postop.SetAGridFunction(A[i]);
    double Um = postop.GetHFieldEnergy();
    PostprocessDomains(postop, "i", i, idx, 0.0, Um, 0.0, 0.0);
    PostprocessSurfaces(postop, "i", i, idx, 0.0, Um, 0.0, Iinc(i));
    PostprocessProbes(postop, "i", i, idx);
    if (i < iodata.solver.magnetostatic.n_post)
    {
      PostprocessFields(postop, i, idx, (i == 0) ? &indicator : nullptr);
      Mpi::Print(" Wrote fields to disk for terminal {:d}\n", idx);
    }
    if (i == 0)
    {
      PostprocessErrorIndicator(postop, indicator);
    }

    // Diagonal: M_ii = 2 U_m(A_i) / I_i².
    M(i, i) = Mm(i, i) = 2.0 * Um / (Iinc(i) * Iinc(i));
    i++;
  }

  // Off-diagonals: M_ij = U_m(A_i + A_j) / (I_i I_j) - 1/2 (I_i/I_j M_ii + I_j/I_i M_jj).
  for (i = 0; i < M.Height(); i++)
  {
    for (int j = 0; j < M.Width(); j++)
    {
      if (j < i)
      {
        // Copy lower triangle from already computed upper triangle.
        M(i, j) = M(j, i);
        Mm(i, j) = Mm(j, i);
        Mm(i, i) -= Mm(i, j);
      }
      else if (j > i)
      {
        linalg::AXPBYPCZ(1.0, A[i], 1.0, A[j], 0.0, Aij);
        Curl.Mult(Aij, B);
        postop.SetBGridFunction(B);
        double Um = postop.GetHFieldEnergy();
        M(i, j) = Um / (Iinc(i) * Iinc(j)) -
                  0.5 * (M(i, i) * Iinc(i) / Iinc(j) + M(j, j) * Iinc(j) / Iinc(i));
        Mm(i, j) = -M(i, j);
        Mm(i, i) -= Mm(i, j);
      }
    }
  }
  mfem::DenseMatrix Minv(M);
  Minv.Invert();  // In-place, uses LAPACK (when available) and should be cheap
  PostprocessTerminals(surf_j_op, M, Minv, Mm);
  return indicator;
}

void MagnetostaticSolver::PostprocessTerminals(const SurfaceCurrentOperator &surf_j_op,
                                               const mfem::DenseMatrix &M,
                                               const mfem::DenseMatrix &Minv,
                                               const mfem::DenseMatrix &Mm) const
{
  // Only root writes to disk (every process has full matrices).
  if (!root || post_dir.length() == 0)
  {
    return;
  }

  // Write inductance matrix data.
  auto PrintMatrix = [&surf_j_op, this](const std::string &file, const std::string &name,
                                        const std::string &unit,
                                        const mfem::DenseMatrix &mat, double scale)
  {
    std::string path = post_dir + file;
    auto output = OutputFile(path, false);
    output.print("{:>{}s},", "i", table.w1);
    for (const auto &[idx2, data2] : surf_j_op)
    {
      // clang-format off
      output.print("{:>{}s}{}",
                   name + "[i][" + std::to_string(idx2) + "] " + unit, table.w,
                   (idx2 == surf_j_op.rbegin()->first) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
    int i = 0;
    for (const auto &[idx, data] : surf_j_op)
    {
      int j = 0;
      output.print("{:{}.{}e},", static_cast<double>(idx), table.w1, table.p1);
      for (const auto &[idx2, data2] : surf_j_op)
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     mat(i, j) * scale, table.w, table.p,
                     (idx2 == surf_j_op.rbegin()->first) ? "" : ",");
        // clang-format on
        j++;
      }
      output.print("\n");
      i++;
    }
  };
  const double H = iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, 1.0);
  PrintMatrix("terminal-M.csv", "M", "(H)", M, H);
  PrintMatrix("terminal-Minv.csv", "M⁻¹", "(1/H)", Minv, 1.0 / H);
  PrintMatrix("terminal-Mm.csv", "M_m", "(H)", Mm, H);
}

}  // namespace palace
