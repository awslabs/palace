// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "magnetostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/curlcurloperator.hpp"
#include "fem/postoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"
#include "linalg/gmg.hpp"
#include "linalg/pc.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

BaseSolver::SolveOutput
MagnetostaticSolver::Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                           Timer &timer) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs.
  timer.Lap();
  std::vector<std::unique_ptr<mfem::Operator>> K;
  CurlCurlOperator curlcurlop(iodata, mesh);
  curlcurlop.GetStiffnessMatrix(K);
  SaveMetadata(curlcurlop.GetNDSpace());

  // Set up the linear solver.
  std::unique_ptr<mfem::Solver> pc =
      ConfigurePreconditioner(iodata, curlcurlop.GetDbcMarker(), curlcurlop.GetNDSpaces());
  auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pc.get());
  if (gmg)
  {
    gmg->SetOperator(K);
  }
  else
  {
    pc->SetOperator(*K.back());
  }

  mfem::IterativeSolver::PrintLevel print =
      mfem::IterativeSolver::PrintLevel().Warnings().Errors();
  if (iodata.problem.verbose > 0)
  {
    print.Summary();
    if (iodata.problem.verbose > 1)
    {
      print.Iterations();
      if (iodata.problem.verbose > 2)
      {
        print.All();
      }
    }
  }
  mfem::CGSolver pcg(mesh.back()->GetComm());
  pcg.SetRelTol(iodata.solver.linear.tol);
  pcg.SetMaxIter(iodata.solver.linear.max_it);
  pcg.SetPrintLevel(print);
  pcg.SetOperator(*K.back());  // Call before SetPreconditioner, PC operator set separately
  pcg.SetPreconditioner(*pc);
  if (iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::DEFAULT &&
      iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::CG)
  {
    Mpi::Warning("Magnetostatic problem type always uses CG as the Krylov solver!\n");
  }

  // Terminal indices are the set of boundaries over which to compute the inductance matrix.
  PostOperator postop(iodata, curlcurlop, "magnetostatic");
  int nstep = static_cast<int>(curlcurlop.GetSurfaceCurrentOp().Size());
  MFEM_VERIFY(nstep > 0,
              "No surface current boundaries specified for magnetostatic simulation!");

  // Source term and solution vector storage.
  mfem::Vector RHS(K.back()->Height());
  std::vector<mfem::Vector> A(nstep);
  timer.construct_time += timer.Lap();

  // Main loop over current source boundaries.
  Mpi::Print("\nComputing magnetostatic fields for {:d} source boundar{}\n", nstep,
             (nstep > 1) ? "ies" : "y");
  int step = 0, ksp_it = 0;
  auto t0 = timer.Now();
  for (const auto &[idx, data] : curlcurlop.GetSurfaceCurrentOp())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, nstep,
               idx, Timer::Duration(timer.Now() - t0).count());

    // Form and solve the linear system for a prescribed current on the specified source.
    Mpi::Print("\n");
    A[step].SetSize(RHS.Size());
    A[step] = 0.0;
    curlcurlop.GetExcitationVector(idx, RHS);
    timer.construct_time += timer.Lap();

    pcg.Mult(RHS, A[step]);
    if (!pcg.GetConverged())
    {
      Mpi::Warning("Linear solver did not converge in {:d} iterations!\n",
                   pcg.GetNumIterations());
    }
    ksp_it += pcg.GetNumIterations();
    timer.solve_time += timer.Lap();

    // A[step]->Print();
    Mpi::Print(" Sol. ||A|| = {:.6e} (||RHS|| = {:.6e})\n",
               std::sqrt(mfem::InnerProduct(mesh.back()->GetComm(), A[step], A[step])),
               std::sqrt(mfem::InnerProduct(mesh.back()->GetComm(), RHS, RHS)));
    timer.postpro_time += timer.Lap();

    // Next source.
    step++;
  }

  // Postprocess the capacitance matrix from the computed field solutions.
  const auto io_time_prev = timer.io_time;
  SaveMetadata(nstep, ksp_it);
  Postprocess(curlcurlop, postop, A, timer);
  timer.postpro_time += timer.Lap() - (timer.io_time - io_time_prev);

  return BaseSolver::SolveOutput();
}

void MagnetostaticSolver::Postprocess(CurlCurlOperator &curlcurlop, PostOperator &postop,
                                      const std::vector<mfem::Vector> &A,
                                      Timer &timer) const
{
  // Postprocess the Maxwell inductance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the magnetic field energy based on a current
  // excitation for each port. Alternatively, we could compute the resulting loop fluxes to
  // get M directly as:
  //                         Φ_i = ∫ B ⋅ n_j dS
  // and M_ij = Φ_i/I_j. The energy formulation avoids having to locally integrate B =
  // ∇ x A.
  std::unique_ptr<mfem::Operator> Curl = curlcurlop.GetCurlMatrix();
  const SurfaceCurrentOperator &surf_j_op = curlcurlop.GetSurfaceCurrentOp();
  int nstep = static_cast<int>(surf_j_op.Size());
  mfem::DenseMatrix M(nstep), Mm(nstep);
  mfem::Vector B(Curl->Height()), Aij(Curl->Width());
  mfem::Vector Iinc(nstep);
  if (iodata.solver.magnetostatic.n_post > 0)
  {
    Mpi::Print("\n");
  }
  int i = 0;
  for (const auto &[idx, data] : surf_j_op)
  {
    // Get the magnitude of the current excitations (unit J_s,inc, but circuit current I is
    // the integral of J_s,inc over port).
    Iinc(i) = data.GetExcitationCurrent();
    MFEM_VERIFY(Iinc(i) > 0.0, "Zero current excitation for magnetostatic solver!");

    // Set the internal GridFunctions in PostOperator for all postprocessing operations.
    PostOperator::GetBField(*Curl, A[i], B);
    postop.SetBGridFunction(B);
    postop.SetAGridFunction(A[i]);
    double Um = postop.GetHFieldEnergy();
    PostprocessDomains(postop, "i", i, idx, 0.0, Um, 0.0, 0.0);
    PostprocessSurfaces(postop, "i", i, idx, 0.0, Um, 0.0, Iinc(i));
    PostprocessProbes(postop, "i", i, idx);
    if (i < iodata.solver.magnetostatic.n_post)
    {
      auto t0 = timer.Now();
      PostprocessFields(postop, i, idx);
      Mpi::Print(" Wrote fields to disk for terminal {:d}\n", idx);
      timer.io_time += timer.Now() - t0;
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
        add(A[i], A[j], Aij);
        PostOperator::GetBField(*Curl, Aij, B);
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
