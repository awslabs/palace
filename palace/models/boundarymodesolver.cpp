// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/arpack.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"

namespace palace
{

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

BoundaryModeSolver::BoundaryModeSolver(const BoundaryModeSolverConfig &config,
                                       const FiniteElementSpace &nd_fespace,
                                       const FiniteElementSpace &h1_fespace,
                                       const mfem::Array<int> &dbc_tdof_list,
                                       MPI_Comm solver_comm)
  : config(config), nd_fespace(nd_fespace), h1_fespace(h1_fespace),
    dbc_tdof_list(dbc_tdof_list)
{
  nd_size = nd_fespace.GetTrueVSize();
  h1_size = h1_fespace.GetTrueVSize();

  // Assemble frequency-independent matrices. These use the FE space communicator
  // (all processes that share the mesh), NOT solver_comm.
  std::tie(Atnr, Atni) = AssembleAtn();
  std::tie(Antr, Anti) = AssembleAnt();
  std::tie(Annr, Anni) = AssembleAnn();

  // Assemble Btt and build the block B matrix.
  {
    auto [bttr, btti] = AssembleBtt();
    Bttr = std::make_unique<mfem::HypreParMatrix>(*bttr);  // Keep a copy for external use

    // Construct the zero diagonal block for the H1 portion of B.
    Vector d(h1_size);
    d.UseDevice(false);  // SparseMatrix constructor uses Vector on host
    d = 0.0;
    mfem::SparseMatrix diag(d);
    auto Dnn = std::make_unique<mfem::HypreParMatrix>(
        h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
        h1_fespace.Get().GetTrueDofOffsets(), &diag);

    auto [Br, Bi] = BuildSystemMatrixB(bttr.get(), btti.get(), Dnn.get());
    opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));
  }

  // Configure linear and eigenvalue solvers.
  // - ModeAnalysis: solver_comm == MPI_COMM_NULL → use FE space communicator (all procs).
  // - WavePort on port procs: solver_comm is a valid subset communicator → use it.
  // - WavePort on non-port procs: solver_comm == MPI_COMM_NULL but FE space has no local
  //   DOFs → skip solver setup (ksp and eigen remain null).
  if (solver_comm != MPI_COMM_NULL)
  {
    SetUpLinearSolver(solver_comm);
    SetUpEigenSolver(solver_comm);
  }
  else if (nd_size > 0)
  {
    // Standalone mode (ModeAnalysis): all procs have DOFs, use FE space comm.
    SetUpLinearSolver(nd_fespace.GetComm());
    SetUpEigenSolver(nd_fespace.GetComm());
  }
  // else: non-port process with empty FE space — no solvers needed.
}

BoundaryModeSolver::~BoundaryModeSolver()
{
  // Free the solvers before any communicator they depend on is freed externally.
  ksp.reset();
  eigen.reset();
}

void BoundaryModeSolver::AssembleFrequencyDependent(double omega, double sigma)
{
  // Assemble frequency-dependent Att and build the block A matrix. Impedance/absorbing/
  // conductivity boundary terms enter Att as a tangential ND mass from the et component
  // of the impedance BC: (iω/Zs)(et·t̂, ft·t̂)_γ.
  //
  // The impedance BC also has an en component: (iω/Zs)(en, fn)_γ, which would enter Ann.
  // However, with the variable substitution en_code = ikn × en_phys, this becomes
  // (iω/Zs)/(ikn) × (en_code, fn)_γ, introducing kn-dependence (nonlinear eigenvalue).
  // For now, Ann boundary terms are omitted and H1 Dirichlet is enforced on all conductor
  // boundaries (en=0), which approximates the PEC condition for en.
  auto [Attr, Atti] = AssembleAtt(omega, sigma);
  auto [Ar, Ai] = BuildSystemMatrixA(Attr.get(), Atti.get(), Atnr.get(), Atni.get(),
                                     Antr.get(), Anti.get(), Annr.get(), Anni.get());
  opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

BoundaryModeSolver::SolveResult
BoundaryModeSolver::Solve(double omega, double sigma, const ComplexVector *initial_space)
{
  // Assemble frequency-dependent matrices (MPI collective on FE space communicator).
  AssembleFrequencyDependent(omega, sigma);

  // Solve the eigenvalue problem. All processes must have solvers for this path.
  MFEM_VERIFY(ksp && eigen,
              "BoundaryModeSolver::Solve called on process without solvers! "
              "Use SolveSplit() for wave port mode where only a subset has solvers.");
  ComplexWrapperOperator opP(opA->Real(), nullptr);  // Non-owning constructor
  ksp->SetOperators(*opA, opP);
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);

  if (initial_space)
  {
    eigen->SetInitialSpace(*initial_space);
  }

  int num_conv = eigen->Solve();

  return {num_conv, sigma};
}

BoundaryModeSolver::SolveResult
BoundaryModeSolver::SolveSplit(double omega, double sigma, bool has_solver,
                               const ComplexVector *initial_space)
{
  // Assemble frequency-dependent matrices (MPI collective on all processes).
  AssembleFrequencyDependent(omega, sigma);

  // Only processes with solvers (port_comm) do the eigenvalue solve.
  if (!has_solver)
  {
    return {0, sigma};
  }

  ComplexWrapperOperator opP(opA->Real(), nullptr);
  ksp->SetOperators(*opA, opP);
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);

  if (initial_space)
  {
    eigen->SetInitialSpace(*initial_space);
  }

  int num_conv = eigen->Solve();

  return {num_conv, sigma};
}

std::complex<double> BoundaryModeSolver::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(i);
}

void BoundaryModeSolver::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(i, x);
}

// Stiffness matrix (shifted): Att = (mu_cc^{-1} curl_t x u, curl_t x v)
//                                   - omega^2 (eps u, v)
//                                   - sigma (mu^{-1} u, v).
BoundaryModeSolver::ComplexHypreParMatrix
BoundaryModeSolver::AssembleAtt(double omega, double sigma) const
{
  MaterialPropertyCoefficient muinv_cc_func(*config.attr_to_material,
                                            *config.curlcurl_inv_permeability);
  if (config.normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*config.normal);
  }

  MaterialPropertyCoefficient eps_shifted_func(*config.attr_to_material,
                                               *config.permittivity_real, -omega * omega);
  eps_shifted_func.AddCoefficient(*config.attr_to_material, *config.inv_permeability,
                                  -sigma);

  // Assemble Att real part: domain + boundary contributions in a single BilinearForm,
  // following the same pattern as SpaceOperator. Boundary impedance/absorbing terms use
  // Palace's ceed-based BilinearForm with VectorFEMassIntegrator.
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);

  // Add boundary contributions from impedance, absorbing, and conductivity operators.
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);        // Ls -> real
    config.surf_z_op->AddDampingBdrCoefficients(omega, fbi);        // Rs -> imag (ω/Rs)
    config.surf_z_op->AddMassBdrCoefficients(-omega * omega, fbr);  // Cs -> real (-ω²Cs)
  }
  if (config.farfield_op)
  {
    config.farfield_op->AddDampingBdrCoefficients(omega, fbi);  // 1st-order ABC
  }
  if (config.surf_sigma_op)
  {
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  }

  Mpi::Print(" [BMS] Att boundary: surf_z_op={}, farfield_op={}, surf_sigma_op={}, "
             "max_bdr_attr={}, fbr.empty={}, fbi.empty={}\n",
             config.surf_z_op != nullptr, config.farfield_op != nullptr,
             config.surf_sigma_op != nullptr, max_bdr_attr, fbr.empty(), fbi.empty());

  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_func);
  if (!fbr.empty())
  {
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }

  // Assemble domain + boundary contributions via Palace's BilinearForm (libCEED).
  auto Attr_assembled =
      ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
  {
    double norm_r =
        hypre_ParCSRMatrixFnorm(static_cast<hypre_ParCSRMatrix *>(*Attr_assembled));
    Mpi::Print(" [BMS] Attr Frobenius norm = {:.6e}\n", norm_r);
  }

  // Assemble imaginary part: domain (loss tangent, conductivity) + boundary (Rs, ABC).
  std::unique_ptr<mfem::HypreParMatrix> Atti_assembled;
  {
    bool has_imag = config.has_loss_tangent ||
                    (config.has_conductivity && config.conductivity) || !fbi.empty();
    if (has_imag)
    {
      BilinearForm atti(nd_fespace);
      if (config.has_loss_tangent)
      {
        MaterialPropertyCoefficient negepstandelta_func(
            *config.attr_to_material, *config.permittivity_imag, -omega * omega);
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
      }
      if (config.has_conductivity && config.conductivity)
      {
        MaterialPropertyCoefficient fi_domain(*config.attr_to_material,
                                              *config.conductivity, omega);
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(fi_domain);
      }
      if (!fbi.empty())
      {
        atti.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbi);
      }
      Atti_assembled =
          ParOperator(atti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
    }
  }

  return {std::move(Attr_assembled), std::move(Atti_assembled)};
}

// Coupling matrix: Atn = -(mu^{-1} grad_t u, v).
BoundaryModeSolver::ComplexHypreParMatrix BoundaryModeSolver::AssembleAtn() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material, *config.inv_permeability,
                                         -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

// Coupling matrix: Ant = -(eps u, grad_t v).
BoundaryModeSolver::ComplexHypreParMatrix BoundaryModeSolver::AssembleAnt() const
{
  MaterialPropertyCoefficient epsilon_func(*config.attr_to_material,
                                           *config.permittivity_real, 1.0);
  BilinearForm antr(nd_fespace, h1_fespace);
  antr.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(epsilon_func);

  // Contribution for loss tangent: eps -> eps * (1 - i tan(delta)).
  if (!config.has_loss_tangent)
  {
    return {ParOperator(antr.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
                .StealParallelAssemble(),
            nullptr};
  }

  MaterialPropertyCoefficient negepstandelta_func(*config.attr_to_material,
                                                  *config.permittivity_imag, 1.0);
  BilinearForm anti(nd_fespace, h1_fespace);
  anti.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(negepstandelta_func);
  return {ParOperator(antr.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
              .StealParallelAssemble(),
          ParOperator(anti.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
              .StealParallelAssemble()};
}

// Mass matrix: Ann = -(eps u, v), with optional normal projection for 3D boundaries.
BoundaryModeSolver::ComplexHypreParMatrix BoundaryModeSolver::AssembleAnn() const
{
  MaterialPropertyCoefficient epsilon_func(*config.attr_to_material,
                                           *config.permittivity_scalar, -1.0);
  if (config.normal)
  {
    epsilon_func.NormalProjectedCoefficient(*config.normal);
  }
  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<MassIntegrator>(epsilon_func);

  auto Annr_assembled =
      ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  // Contribution for loss tangent: eps -> eps * (1 - i tan(delta)).
  if (!config.has_loss_tangent)
  {
    return {std::move(Annr_assembled), nullptr};
  }

  MaterialPropertyCoefficient negepstandelta_func(*config.attr_to_material,
                                                  *config.permittivity_imag, -1.0);
  if (config.normal)
  {
    negepstandelta_func.NormalProjectedCoefficient(*config.normal);
  }
  BilinearForm anni(h1_fespace);
  anni.AddDomainIntegrator<MassIntegrator>(negepstandelta_func);
  return {std::move(Annr_assembled),
          ParOperator(anni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble()};
}

// Mass matrix: Btt = (mu^{-1} u, v).
BoundaryModeSolver::ComplexHypreParMatrix BoundaryModeSolver::AssembleBtt() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material,
                                         *config.inv_permeability);
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

BoundaryModeSolver::ComplexHypreParMatrix BoundaryModeSolver::BuildSystemMatrixA(
    const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
    const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
    const mfem::HypreParMatrix *Antr, const mfem::HypreParMatrix *Anti,
    const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Attr;
  blocks(0, 1) = Atnr;
  blocks(1, 0) = Antr;
  blocks(1, 1) = Annr;
  std::unique_ptr<mfem::HypreParMatrix> Ar(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti)
  {
    blocks(0, 0) = Atti;
    blocks(0, 1) = Atni;
    blocks(1, 0) = Anti;
    blocks(1, 1) = Anni;
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs constrained by Dirichlet BCs.
  Ar->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  if (Ai)
  {
    Ai->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Ar), std::move(Ai)};
}

BoundaryModeSolver::ComplexHypreParMatrix
BoundaryModeSolver::BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                       const mfem::HypreParMatrix *Btti,
                                       const mfem::HypreParMatrix *Dnn) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Bttr;
  blocks(0, 1) = nullptr;
  blocks(1, 0) = nullptr;
  blocks(1, 1) = Dnn;
  std::unique_ptr<mfem::HypreParMatrix> Br(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Bi;
  if (Btti)
  {
    blocks(0, 0) = Btti;
    Bi.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs constrained by Dirichlet BCs.
  Br->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (Bi)
  {
    Bi->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Br), std::move(Bi)};
}

void BoundaryModeSolver::SetUpLinearSolver(MPI_Comm comm)
{
  // GMRES iterative solver preconditioned with a sparse direct solver.
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, config.verbose);
  gmres->SetInitialGuess(false);
  gmres->SetRelTol(config.linear->tol);
  gmres->SetMaxIter(config.linear->max_it);
  gmres->SetRestartDim(config.linear->max_it);

  // Select sparse direct solver type. Multigrid (AMS, BoomerAMG) is not applicable for
  // the combined ND+H1 block system, so fall back to a sparse direct solver.
  LinearSolver pc_type = config.linear->type;
  if (pc_type == LinearSolver::DEFAULT || pc_type == LinearSolver::AMS ||
      pc_type == LinearSolver::BOOMER_AMG)
  {
#if defined(MFEM_USE_SUPERLU)
    pc_type = LinearSolver::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
    pc_type = LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
    pc_type = LinearSolver::MUMPS;
#else
    MFEM_ABORT("Boundary mode solver requires building with SuperLU_DIST, STRUMPACK, "
               "or MUMPS!");
#endif
  }
  else if (pc_type == LinearSolver::SUPERLU)
  {
#if !defined(MFEM_USE_SUPERLU)
    MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
               "different solver!");
#endif
  }
  else if (pc_type == LinearSolver::STRUMPACK || pc_type == LinearSolver::STRUMPACK_MP)
  {
#if !defined(MFEM_USE_STRUMPACK)
    MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
               "different solver!");
#endif
  }
  else if (pc_type == LinearSolver::MUMPS)
  {
#if !defined(MFEM_USE_MUMPS)
    MFEM_ABORT("Solver was not built with MUMPS support, please choose a "
               "different solver!");
#endif
  }

  const auto &linear = *config.linear;
  auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
      [&]() -> std::unique_ptr<mfem::Solver>
      {
        if (pc_type == LinearSolver::SUPERLU)
        {
#if defined(MFEM_USE_SUPERLU)
          return std::make_unique<SuperLUSolver>(
              comm, linear.sym_factorization, linear.superlu_3d, true, config.verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::STRUMPACK ||
                 pc_type == LinearSolver::STRUMPACK_MP)
        {
#if defined(MFEM_USE_STRUMPACK)
          return std::make_unique<StrumpackSolver>(
              comm, linear.sym_factorization, linear.strumpack_compression_type,
              linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
              linear.strumpack_lossy_precision, true, config.verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::MUMPS)
        {
#if defined(MFEM_USE_MUMPS)
          return std::make_unique<MumpsSolver>(
              comm, mfem::MUMPSSolver::UNSYMMETRIC, linear.sym_factorization,
              linear.strumpack_lr_tol, true, config.verbose - 1);
#endif
        }
        MFEM_ABORT("Unsupported linear solver type for boundary mode solver!");
        return {};
      }());
  pc->SetSaveAssembled(false);
  pc->SetDropSmallEntries(false);
  ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
}

void BoundaryModeSolver::SetUpEigenSolver(MPI_Comm comm)
{
  constexpr int print = 0;
  EigenSolverBackend type = config.eigen_backend;

  if (type == EigenSolverBackend::SLEPC)
  {
#if !defined(PALACE_WITH_SLEPC)
    MFEM_ABORT("Solver was not built with SLEPc support, please choose a "
               "different solver!");
#endif
  }
  else if (type == EigenSolverBackend::ARPACK)
  {
#if !defined(PALACE_WITH_ARPACK)
    MFEM_ABORT("Solver was not built with ARPACK support, please choose a "
               "different solver!");
#endif
  }
  else  // Default choice
  {
#if defined(PALACE_WITH_SLEPC)
    type = EigenSolverBackend::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
    type = EigenSolverBackend::ARPACK;
#else
#error "Boundary mode solver requires building with ARPACK or SLEPc!"
#endif
  }

  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    auto arpack = std::make_unique<arpack::ArpackEPSSolver>(comm, print);
    arpack->SetNumModes(config.num_modes, config.num_vec);
    arpack->SetTol(config.eig_tol);
    arpack->SetWhichEigenpairs(config.which_eig);
    arpack->SetLinearSolver(*ksp);
    eigen = std::move(arpack);
#endif
  }
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    auto slepc = std::make_unique<slepc::SlepcEPSSolver>(comm, print);
    slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetNumModes(config.num_modes, config.num_vec);
    slepc->SetTol(config.eig_tol);
    slepc->SetWhichEigenpairs(config.which_eig);
    slepc->SetLinearSolver(*ksp);
    eigen = std::move(slepc);
#endif
  }

  // Note: gradient null-space modes (kn = 0, et = grad phi) are separated from physical
  // propagating modes by the shift-and-invert spectral transformation. Physical modes have
  // large |lambda|, gradient modes have small |lambda| ~ 1/|sigma|, so LARGEST_MAGNITUDE
  // targeting naturally skips them. A DivFreeSolver is NOT appropriate here: it would
  // project et to tangential divergence-free and destroy physical modes whose divergence
  // couples to en through the block equation (div_t(eps et) = eps en != 0).
}

}  // namespace palace
