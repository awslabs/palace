// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <mfem.hpp>
#include "linalg/arpack.hpp"
#include "linalg/complex.hpp"
#include "linalg/divfree.hpp"
#include "linalg/feast.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/slepc.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

void EigenSolver::Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                        Timer &timer) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  timer.Lap();
  SpaceOperator spaceop(iodata, mesh);
  std::unique_ptr<ComplexParOperator> K = spaceop.GetComplexSystemMatrix(
      SpaceOperator::OperatorType::STIFFNESS, Operator::DIAG_ONE);
  std::unique_ptr<ComplexParOperator> M = spaceop.GetComplexSystemMatrix(
      SpaceOperator::OperatorType::MASS, Operator::DIAG_ZERO);
  std::unique_ptr<ComplexParOperator> C = spaceop.GetComplexSystemMatrix(
      SpaceOperator::OperatorType::DAMPING, Operator::DIAG_ZERO);
  std::unique_ptr<ComplexParOperator> Curl = spaceop.GetComplexCurlMatrix();
  SaveMetadata(spaceop.GetNDSpace());

  // Configure objects for postprocessing.
  PostOperator postop(iodata, spaceop, "eigenmode");
  ComplexVector E(Curl->Width()), B(Curl->Height());

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. A shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ. In general, the system matrices are complex and
  // symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  config::EigenSolverData::Type type = iodata.solver.eigenmode.type;
#if defined(PALACE_WITH_ARPACK) && defined(PALACE_WITH_SLEPC)
  if (type == config::EigenSolverData::Type::DEFAULT)
  {
    type = config::EigenSolverData::Type::SLEPC;
  }
#elif defined(PALACE_WITH_ARPACK)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::SLEPC)
  {
    Mpi::Warning("SLEPc eigensolver not available, using ARPACK!\n");
  }
  else if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::FEAST)
  {
    Mpi::Warning("FEAST eigensolver requires SLEPc, using ARPACK!\n");
  }
  type = config::EigenSolverData::Type::ARPACK;
#elif defined(PALACE_WITH_SLEPC)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::ARPACK)
  {
    Mpi::Warning("ARPACK eigensolver not available, using SLEPc!\n");
  }
  type = config::EigenSolverData::Type::SLEPC;
#else
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == config::EigenSolverData::Type::FEAST)
  {
    MFEM_ABORT("FEAST eigenvalue solver is currently not supported!");
#if defined(PALACE_WITH_SLEPC)
    // Mpi::Print("\nConfiguring FEAST eigenvalue solver\n");
    // if (C)
    // {
    //   eigen = std::make_unique<feast::FeastPEPSolver>(
    //       K->GetComm(), iodata, spaceop, iodata.solver.eigenmode.feast_contour_np,
    //       iodata.problem.verbose);
    // }
    // else
    // {
    //   eigen = std::make_unique<feast::FeastEPSSolver>(
    //       K->GetComm(), iodata, spaceop, iodata.solver.eigenmode.feast_contour_np,
    //       iodata.problem.verbose);
    // }
#endif
  }
  else if (type == config::EigenSolverData::Type::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver\n");
    if (C)
    {
      eigen =
          std::make_unique<arpack::ArpackPEPSolver>(K->GetComm(), iodata.problem.verbose);
    }
    else
    {
      eigen =
          std::make_unique<arpack::ArpackEPSSolver>(K->GetComm(), iodata.problem.verbose);
    }
#endif
  }
  else  // config::EigenSolverData::Type::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver\n");
    std::unique_ptr<slepc::SlepcEigenSolver> slepc;
    if (C)
    {
      if (!iodata.solver.eigenmode.pep_linear)
      {
        slepc =
            std::make_unique<slepc::SlepcPEPSolver>(K->GetComm(), iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenSolver::Type::TOAR);
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(K->GetComm(),
                                                              iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenSolver::Type::KRYLOVSCHUR);
      }
    }
    else
    {
      slepc = std::make_unique<slepc::SlepcEPSSolver>(K->GetComm(), iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenSolver::Type::KRYLOVSCHUR);
    }
    slepc->SetProblemType(slepc::SlepcEigenSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetOrthogonalization(iodata.solver.linear.orthog_mgs,
                                iodata.solver.linear.orthog_cgs2);
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (C)
  {
    eigen->SetOperators(*K, *C, *M, scale);
  }
  else
  {
    eigen->SetOperators(*K, *M, scale);
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  eigen->SetTol(iodata.solver.eigenmode.tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  const double target = iodata.solver.eigenmode.target;
  const double f_target = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, target);
  std::unique_ptr<ComplexParOperator> A;
  std::vector<std::unique_ptr<ParOperator>> P, AuxP;
  std::unique_ptr<ComplexKspSolver> ksp;
  // #if defined(PALACE_WITH_SLEPC)
  //   auto *feast = dynamic_cast<feast::FeastEigenSolver *>(eigen.get());
  //   if (feast)
  //   {
  //     // Configure the FEAST integration contour. The linear solvers are set up inside
  //     the
  //     // solver.
  //     if (iodata.solver.eigenmode.feast_contour_np > 1)
  //     {
  //       double contour_ub = iodata.solver.eigenmode.feast_contour_ub;
  //       double f_contour_ub =
  //           iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, contour_ub);
  //       double contour_ar = iodata.solver.eigenmode.feast_contour_ar;
  //       MFEM_VERIFY(contour_ub > target,
  //                   "FEAST eigensolver requires a specified upper frequency target!");
  //       MFEM_VERIFY(
  //           contour_ar >= 0.0 && contour_ar <= 1.0,
  //           "Contour aspect ratio for FEAST eigenvalue solver must be in range
  //           [0.0, 1.0]!");
  //       Mpi::Print(" FEAST search contour: σ_lower = {:.3e} GHz ({:.3e})\n"
  //                  "                       σ_upper = {:.3e} GHz ({:.3e})\n"
  //                  "                       AR = {:.1e}\n",
  //                  f_target, target, f_contour_ub, contour_ub, contour_ar);
  //       if (C)
  //       {
  //         // Search for eigenvalues in the range λ = iσₗₒ to iσₕᵢ.
  //         double h = (contour_ub - target) * contour_ar;
  //         feast->SetContour(-0.5 * h, target, 0.5 * h, contour_ub, false, true);
  //       }
  //       else
  //       {
  //         // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues from μ =
  //         σₗₒ² to
  //         // σₕᵢ².
  //         double h = (contour_ub * contour_ub - target * target) * contour_ar;
  //         feast->SetContour(target * target, -0.5 * h, contour_ub * contour_ub, 0.5 * h);
  //       }
  //     }
  //     else
  //     {
  //       Mpi::Print(" FEAST search target: σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  //       if (C)
  //       {
  //         feast->SetContour(0.0, target, 0.0, target, false, true);
  //       }
  //       else
  //       {
  //         feast->SetContour(target * target, 0.0, target * target, 0.0);
  //       }
  //     }
  //   }
  //   else
  // #endif
  {
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
    if (C)
    {
      // Search for eigenvalues closest to λ = iσ.
      eigen->SetShiftInvert(0.0, target);
      if (type == config::EigenSolverData::Type::ARPACK)
      {
        // ARPACK searches based on eigenvalues of the transformed problem. The eigenvalue
        // 1 / (λ - σ) will be a large-magnitude negative imaginary number for an eigenvalue
        // λ with frequency close to but not below the target σ.
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_IMAGINARY);
      }
      else
      {
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
      }
    }
    else
    {
      // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
      eigen->SetShiftInvert(target * target, 0.0);
      if (type == config::EigenSolverData::Type::ARPACK)
      {
        // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
        // will be a large-magnitude positive real number for an eigenvalue μ with frequency
        // close to but below the target σ².
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
      }
      else
      {
        eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
      }
    }

    // Set up the linear solver required for solving systems involving the shifted operator
    // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
    // preconditioner for complex linear systems is constructed from a real approximation
    // to the complex system matrix.
    A = spaceop.GetComplexSystemMatrix(1.0, 1i * target, -target * target, K.get(), C.get(),
                                       M.get(), nullptr);

    spaceop.GetPreconditionerMatrix(1.0, target, -target * target, target, P, AuxP);

    ksp = std::make_unique<ComplexKspSolver>(iodata, spaceop.GetNDSpaces(),
                                             &spaceop.GetH1Spaces());
    ksp->SetOperator(*A, P, &AuxP);
    eigen->SetLinearSolver(*ksp);

    // XX TODO REVISIT FOR LOG PRINT FORMATTING...
  }

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed).
  std::unique_ptr<Operator> Mr;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    // Mpi::Print(" Basis uses M-inner product\n");
    // Mr = std::make_unique<SumOperator>(M->Real(), 1.0);
    // eigen->SetBMat(*Mr);

    Mpi::Print(" Basis uses (K + M)-inner product\n");
    auto KM = std::make_unique<SumOperator>(M->Real(), 1.0);
    KM->AddOperator(K->Real(), 1.0);
    Mr = std::move(KM);
    eigen->SetBMat(*Mr);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver> divfree;
  if (iodata.solver.linear.divfree_max_it > 0)
  {
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver>(
        spaceop.GetMaterialOp(), spaceop.GetNDSpace(), spaceop.GetH1Spaces(),
        spaceop.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      spaceop.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      spaceop.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // auto Grad = spaceop.GetComplexGradMatrix();
    // ComplexVector r0(Grad->Width());
    // Grad->MultTranspose(v0, r0);
    // r0.Print();
  }
  timer.construct_time += timer.Lap();

  // Eigenvalue problem solve.
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
#if defined(PALACE_WITH_SLEPC)
  // auto *feast = dynamic_cast<feast::FeastEigenSolver *>(eigen.get());
  // if (feast)
  // {
  //   SaveMetadata(feast->GetLinearSolver());
  // }
  // else
#endif
  {
    SaveMetadata(*ksp);
  }
  timer.solve_time += timer.Lap();

  // Postprocess the results.
  const auto io_time_prev = timer.io_time;
  for (int i = 0; i < num_conv; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    double error1 = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error2 = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (!C)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }
    if (i == 0)
    {
      Mpi::Print(" Found {:d} converged eigenvalue{} (first = {:.3e}{:+.3e}i)\n", num_conv,
                 (num_conv > 1) ? "s" : "", omega.real(), omega.imag());
      Mpi::Print("\n");
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    eigen->GetEigenvector(i, E);
    Curl->Mult(E, B);
    B *= -1.0 / (1i * omega);
    postop.SetEGridFunction(E);
    postop.SetBGridFunction(B);
    postop.UpdatePorts(spaceop.GetLumpedPortOp(), omega.real());

    // Postprocess the mode.
    Postprocess(postop, spaceop.GetLumpedPortOp(), i, omega, error1, error2, num_conv,
                timer);
  }
  timer.postpro_time += timer.Lap() - (timer.io_time - io_time_prev);
}

void EigenSolver::Postprocess(const PostOperator &postop,
                              const LumpedPortOperator &lumped_port_op, int i,
                              std::complex<double> omega, double error1, double error2,
                              int num_conv, Timer &timer) const
{
  // The internal GridFunctions for PostOperator have already been set from the E and B
  // solutions in the main loop over converged eigenvalues. Note: The energies output are
  // nondimensional (they can be dimensionalized using the scaling μ₀ * H₀² * L₀³, which
  // are the free space permeability, characteristic magnetic field strength, and
  // characteristic length scale, respectively).
  double E_elec = postop.GetEFieldEnergy();
  double E_mag = postop.GetHFieldEnergy();
  double E_cap = postop.GetLumpedCapacitorEnergy(lumped_port_op);
  double E_ind = postop.GetLumpedInductorEnergy(lumped_port_op);
  PostprocessEigen(i, omega, error1, error2, num_conv);
  PostprocessEPR(postop, lumped_port_op, i, omega, E_elec + E_cap);
  PostprocessDomains(postop, "m", i, i + 1, E_elec, E_mag, E_cap, E_ind);
  PostprocessSurfaces(postop, "m", i, i + 1, E_elec + E_cap, E_mag + E_ind, 1.0, 1.0);
  PostprocessProbes(postop, "m", i, i + 1);
  if (i < iodata.solver.eigenmode.n_post)
  {
    auto t0 = timer.Now();
    PostprocessFields(postop, i, i + 1);
    Mpi::Print(" Wrote mode {:d} to disk\n", i + 1);
    timer.io_time += timer.Now() - t0;
  }
}

namespace
{

struct EprLData
{
  const int idx;    // Lumped inductor index
  const double pj;  // Inductor energy-participation ratio
};

struct EprIOData
{
  const int idx;    // Lumped resistor index
  const double Ql;  // Quality factor
  const double Kl;  // κ for loss rate
};

}  // namespace

void EigenSolver::PostprocessEigen(int i, std::complex<double> omega, double error1,
                                   double error2, int num_conv) const
{
  // Dimensionalize the result and print in a nice table of frequencies and Q-factors. Save
  // to file if user has specified.
  const std::complex<double> f = {
      iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.real()),
      iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega.imag())};
  const double Q =
      (f.imag() == 0.0) ? mfem::infinity() : 0.5 * std::abs(f) / std::abs(f.imag());

  // Print table to stdout.
  {
    const int int_width = 1 + static_cast<int>(std::log10(num_conv));
    constexpr int p = 6;
    constexpr int w = 6 + p + 7;  // Column spaces + precision + extra for table
    if (i == 0)
    {
      // clang-format off
      Mpi::Print("{:>{}s}{:>{}s}{:>{}s}{:>{}s}{:>{}s}\n{}\n",
                 "m", int_width,
                 "Re{ω}/2π (GHz)", w,
                 "Im{ω}/2π (GHz)", w,
                 "Bkwd. Error", w,
                 "Abs. Error", w,
                 std::string(int_width + 4 * w, '='));
      // clang-format on
    }
    // clang-format off
    Mpi::Print("{:{}d}{:+{}.{}e}{:+{}.{}e}{:+{}.{}e}{:+{}.{}e}\n",
               i + 1, int_width,
               f.real(), w, p,
               f.imag(), w, p,
               error1, w, p,
               error2, w, p);
    // clang-format on
  }

  // Print table to file.
  if (root && post_dir.length() > 0)
  {
    std::string path = post_dir + "eig.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      // clang-format off
      output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s}\n",
                   "m", table.w1,
                   "Re{f} (GHz)", table.w,
                   "Im{f} (GHz)", table.w,
                   "Q", table.w);
      // clang-format on
    }
    // clang-format off
    output.print("{:{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}\n",
                 static_cast<double>(i + 1), table.w1, table.p1,
                 f.real(), table.w, table.p,
                 f.imag(), table.w, table.p,
                 Q, table.w, table.p);
    // clang-format on
  }
}

void EigenSolver::PostprocessEPR(const PostOperator &postop,
                                 const LumpedPortOperator &lumped_port_op, int i,
                                 std::complex<double> omega, double Em) const
{
  // If ports have been specified in the model, compute the corresponding energy-
  // participation ratios (EPR) and write out to disk.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the mode EPR for lumped inductor elements.
  std::vector<EprLData> epr_L_data;
  epr_L_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.GetL()) > 0.0)
    {
      const double pj = postop.GetInductorParticipation(lumped_port_op, idx, Em);
      epr_L_data.push_back({idx, pj});
    }
  }
  if (root && !epr_L_data.empty())
  {
    std::string path = post_dir + "port-EPR.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      output.print("{:>{}s},", "m", table.w1);
      for (const auto &data : epr_L_data)
      {
        // clang-format off
        output.print("{:>{}s}{}",
                     "p[" + std::to_string(data.idx) + "]", table.w,
                     (data.idx == epr_L_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", static_cast<double>(i + 1), table.w1, table.p1);
    for (const auto &data : epr_L_data)
    {
      // clang-format off
      output.print("{:+{}.{}e}{}",
                   data.pj, table.w, table.p,
                   (data.idx == epr_L_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }

  // Write the mode EPR for lumped resistor elements.
  std::vector<EprIOData> epr_IO_data;
  epr_IO_data.reserve(lumped_port_op.Size());
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.GetR()) > 0.0)
    {
      const double Kl = postop.GetExternalKappa(lumped_port_op, idx, Em);
      const double Ql = (Kl == 0.0) ? mfem::infinity() : omega.real() / std::abs(Kl);
      epr_IO_data.push_back(
          {idx, Ql, iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, Kl)});
    }
  }
  if (root && !epr_IO_data.empty())
  {
    std::string path = post_dir + "port-Q.csv";
    auto output = OutputFile(path, (i > 0));
    if (i == 0)
    {
      output.print("{:>{}s},", "m", table.w1);
      for (const auto &data : epr_IO_data)
      {
        // clang-format off
        output.print("{:>{}s},{:>{}s}{}",
                     "Q_ext[" + std::to_string(data.idx) + "]", table.w,
                     "κ_ext[" + std::to_string(data.idx) + "] (GHz)", table.w,
                     (data.idx == epr_IO_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", static_cast<double>(i + 1), table.w1, table.p1);
    for (const auto &data : epr_IO_data)
    {
      // clang-format off
      output.print("{:+{}.{}e},{:+{}.{}e}{}",
                   data.Ql, table.w, table.p,
                   data.Kl, table.w, table.p,
                   (data.idx == epr_IO_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

}  // namespace palace
