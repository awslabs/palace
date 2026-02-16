// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeanalysissolver.hpp"

#include <complex>

#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/mumps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/modeanalysisoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

std::pair<ErrorIndicator, long long int>
ModeAnalysisSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  MFEM_VERIFY(mesh.back()->Dimension() == 2,
              "ModeAnalysis solver requires a 2D mesh (waveguide cross-section)!");

  const auto &ma_data = iodata.solver.mode_analysis;
  const double freq_GHz = ma_data.freq;
  const int num_modes = ma_data.n;
  const double tol = ma_data.tol;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz (omega = "
             "{:.6e})\n",
             freq_GHz, omega);

  // Construct FE spaces: ND for tangential E, H1 for normal (out-of-plane) E.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                         mesh.back()->Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                         mesh.back()->Dimension());
  FiniteElementSpace nd_fespace(*mesh.back(), nd_fec.get());
  FiniteElementSpace h1_fespace(*mesh.back(), h1_fec.get());

  // PEC essential BCs.
  const auto &pmesh = mesh.back()->Get();
  int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> dbc_marker =
      palace::mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  const int nd_size = nd_fespace.GetTrueVSize();
  const int h1_size = h1_fespace.GetTrueVSize();
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             nd_fespace.GlobalTrueVSize(), h1_fespace.GlobalTrueVSize(),
             nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize());

  // Material operator, mode analysis operator wrapper, and PostOperator.
  MaterialOperator mat_op(iodata, *mesh.back());
  ModeAnalysisOperator mode_op(mat_op, nd_fespace, h1_fespace, *mesh.back(),
                               iodata.solver.order);
  PostOperator<ProblemType::MODEANALYSIS> post_op(iodata, mode_op);
  double c_min = mat_op.GetLightSpeedMax().Min();
  Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light in ModeAnalysisSolver!");
  const double mu_eps_max = 1.0 / (c_min * c_min) * 1.1;
  const double sigma = -omega * omega * mu_eps_max;

  // Assemble block system matrices.
  MaterialPropertyCoefficient muinv_cc(mat_op.GetAttributeToMaterial(),
                                       mat_op.GetCurlCurlInvPermeability());
  MaterialPropertyCoefficient eps_shifted(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetPermittivityReal(), -omega * omega);
  eps_shifted.AddCoefficient(mat_op.GetAttributeToMaterial(),
                             mat_op.GetInvPermeability(), -sigma);
  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc, eps_shifted);
  auto Att = ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  MaterialPropertyCoefficient muinv_neg(mat_op.GetAttributeToMaterial(),
                                        mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_neg);
  auto Atn = ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
                 .StealParallelAssemble();

  MaterialPropertyCoefficient eps_pos(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetPermittivityReal(), 1.0);
  BilinearForm ant(nd_fespace, h1_fespace);
  ant.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(eps_pos);
  auto Ant = ParOperator(ant.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
                 .StealParallelAssemble();

  MaterialPropertyCoefficient eps_neg(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetPermittivityScalar(), -1.0);
  BilinearForm ann(h1_fespace);
  ann.AddDomainIntegrator<MassIntegrator>(eps_neg);
  auto Ann = ParOperator(ann.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  MaterialPropertyCoefficient muinv_pos(mat_op.GetAttributeToMaterial(),
                                        mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_pos);
  auto Btt = ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  // Store Btt in ModeAnalysisOperator for impedance postprocessing in PostOperator.
  // Make a copy since the original is also needed in the block B matrix below.
  {
    auto Btt_copy = std::make_unique<mfem::HypreParMatrix>(*Btt);
    mode_op.SetBttMatrix(std::move(Btt_copy));
  }

  // Block A matrix with essential BC elimination.
  std::unique_ptr<mfem::HypreParMatrix> opAr;
  {
    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    blocks(0, 0) = Att.get();
    blocks(0, 1) = Atn.get();
    blocks(1, 0) = Ant.get();
    blocks(1, 1) = Ann.get();
    opAr.reset(mfem::HypreParMatrixFromBlocks(blocks));
    mfem::Array<int> dbc_tdof_list;
    dbc_tdof_list.Append(nd_dbc_tdof_list);
    for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
    {
      dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
    }
    opAr->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  }

  // Block B matrix.
  std::unique_ptr<mfem::HypreParMatrix> opBr;
  {
    mfem::SparseMatrix zero_sp(h1_size, h1_size);
    zero_sp.Finalize();
    mfem::HypreParMatrix zero_h1(h1_fespace.GetComm(), h1_fespace.GlobalTrueVSize(),
                                  h1_fespace.GlobalTrueVSize(),
                                  h1_fespace.Get().GetTrueDofOffsets(),
                                  h1_fespace.Get().GetTrueDofOffsets(), &zero_sp);
    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    blocks(0, 0) = Btt.get();
    blocks(0, 1) = nullptr;
    blocks(1, 0) = nullptr;
    blocks(1, 1) = &zero_h1;
    opBr.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  auto opA = std::make_unique<ComplexWrapperOperator>(std::move(opAr), nullptr);
  auto opB = std::make_unique<ComplexWrapperOperator>(std::move(opBr), nullptr);

  // Eigenvalue solver.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving for {:d} propagation mode(s)...\n", num_modes);

  std::unique_ptr<EigenvalueSolver> eigen;
#if defined(PALACE_WITH_SLEPC)
  {
    auto slepc = std::make_unique<slepc::SlepcEPSSolver>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetNumModes(num_modes, std::max(2 * num_modes + 1, 20));
    slepc->SetTol(tol);
    slepc->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    eigen = std::move(slepc);
  }
#elif defined(PALACE_WITH_ARPACK)
  {
    auto arpack = std::make_unique<arpack::ArpackEPSSolver>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    arpack->SetNumModes(num_modes, std::max(2 * num_modes + 1, 20));
    arpack->SetTol(tol);
    arpack->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    eigen = std::move(arpack);
  }
#else
  MFEM_ABORT("ModeAnalysis solver requires SLEPc or ARPACK!");
#endif

  // Linear solver for shift-and-invert. Uses GMRES preconditioned with a sparse direct
  // solver (same pattern as the wave port solver). The direct solver type is read from the
  // linear solver config, defaulting to whatever sparse direct solver is available.
  // Note: multigrid (AMS, BoomerAMG) is not applicable for the block eigenvalue problem
  // since the combined ND+H1 system doesn't have a standard multigrid hierarchy.
  std::unique_ptr<ComplexKspSolver> ksp;
  {
    const auto &linear = iodata.solver.linear;
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    gmres->SetInitialGuess(false);
    gmres->SetRelTol(linear.tol);
    gmres->SetMaxIter(linear.max_it);
    gmres->SetRestartDim(linear.max_it);

    // Select sparse direct solver from config (or default to whatever is available).
    LinearSolver pc_type = linear.type;
    if (pc_type == LinearSolver::DEFAULT || pc_type == LinearSolver::AMS ||
        pc_type == LinearSolver::BOOMER_AMG)
    {
      // Multigrid/AMS not applicable for the block system; fall back to sparse direct.
#if defined(MFEM_USE_SUPERLU)
      pc_type = LinearSolver::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      pc_type = LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      pc_type = LinearSolver::MUMPS;
#else
      MFEM_ABORT(
          "ModeAnalysis solver requires SuperLU_DIST, STRUMPACK, or MUMPS for the "
          "shift-and-invert linear solver!");
#endif
    }
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        [&]() -> std::unique_ptr<mfem::Solver>
        {
          if (pc_type == LinearSolver::SUPERLU)
          {
#if defined(MFEM_USE_SUPERLU)
            return std::make_unique<SuperLUSolver>(
                mesh.back()->GetComm(), linear.sym_factorization,
                linear.superlu_3d, true, iodata.problem.verbose - 1);
#endif
          }
          else if (pc_type == LinearSolver::STRUMPACK ||
                   pc_type == LinearSolver::STRUMPACK_MP)
          {
#if defined(MFEM_USE_STRUMPACK)
            return std::make_unique<StrumpackSolver>(
                mesh.back()->GetComm(), linear.sym_factorization,
                linear.strumpack_compression_type, linear.strumpack_lr_tol,
                linear.strumpack_lossy_precision, linear.strumpack_butterfly_l, true,
                iodata.problem.verbose - 1);
#endif
          }
          else if (pc_type == LinearSolver::MUMPS)
          {
#if defined(MFEM_USE_MUMPS)
            return std::make_unique<MumpsSolver>(
                mesh.back()->GetComm(), mfem::MUMPSSolver::UNSYMMETRIC,
                linear.sym_factorization, linear.strumpack_lr_tol, true,
                iodata.problem.verbose - 1);
#endif
          }
          MFEM_ABORT("Unsupported linear solver type for ModeAnalysis!");
          return {};
        }());
    pc->SetSaveAssembled(false);
    pc->SetDropSmallEntries(false);
    ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
  }
  ComplexWrapperOperator opP(opA->Real(), nullptr);
  ksp->SetOperators(*opA, opP);
  eigen->SetLinearSolver(*ksp);
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);

  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  // Postprocessing: extract propagation constants, effective indices, and impedance.
  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  Mpi::Print("\nComputing mode analysis results and performing postprocessing\n\n");
  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    std::complex<double> lambda = eigen->GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);

    // Extract the tangential E eigenvector for PostOperator (only the ND portion).
    ComplexVector et(nd_size);
    {
      ComplexVector e0(nd_size + h1_size);
      eigen->GetEigenvector(i, e0);
      std::copy_n(e0.Real().begin(), nd_size, et.Real().begin());
      std::copy_n(e0.Imag().begin(), nd_size, et.Imag().begin());
    }

    // PostOperator handles all measurements, printing, and CSV output.
    post_op.MeasureAndPrintAll(i, et, kn, omega, n_print);
  }
  Mpi::Print("\n");

  ErrorIndicator indicator;
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
