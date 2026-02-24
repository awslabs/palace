// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "newboundarymodesolver.hpp"

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

NewBoundaryModeSolver::NewBoundaryModeSolver(const BoundaryModeSolverConfig &config,
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
  //
  // Atn: gradient coupling -(mu^{-1} grad_t u, v), same as BoundaryModeSolver.
  std::tie(Atnr, Atni) = AssembleAtn();

  // Btn: NEGATED transpose of Atn. Atn = -(mu^{-1} grad ·, ·), so Atn^T = -(mu^{-1} ·, grad·).
  // But the physical Btn from Eq 2 is +(mu^{-1} ·, grad·) (positive), so Btn = -Atn^T.
  Btnr.reset(Atnr->Transpose());
  *Btnr *= -1.0;
  if (Atni)
  {
    Btni.reset(Atni->Transpose());
    *Btni *= -1.0;
  }

  // Assemble Btt (negative mu^{-1} ND mass) and build the block B matrix.
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

    auto [Br, Bi] = BuildSystemMatrixB(bttr.get(), btti.get(), Btnr.get(), Btni.get(),
                                       Dnn.get());
    opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));
  }

  // Configure linear and eigenvalue solvers.
  // - ModeAnalysis: solver_comm == MPI_COMM_NULL -> use FE space communicator (all procs).
  // - WavePort on port procs: solver_comm is a valid subset communicator -> use it.
  // - WavePort on non-port procs: solver_comm == MPI_COMM_NULL but FE space has no local
  //   DOFs -> skip solver setup (ksp and eigen remain null).
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
  // else: non-port process with empty FE space -- no solvers needed.
}

NewBoundaryModeSolver::~NewBoundaryModeSolver()
{
  // Free the solvers before any communicator they depend on is freed externally.
  ksp.reset();
  eigen.reset();
}

void NewBoundaryModeSolver::AssembleFrequencyDependent(double omega, double sigma)
{
  // Assemble frequency-dependent Att (same as BoundaryModeSolver) and Ann (NEW: H1
  // stiffness + omega^2 eps mass + BC-n impedance), then build the block A matrix.
  // The shift-and-invert transformation requires A_shifted = A - sigma * B. Since B has
  // an off-diagonal Btn block, the shifted A has a (2,1) block: -sigma * Btn.
  // The Att shift is handled inside AssembleAtt (via the sigma parameter).
  auto [Attr, Atti] = AssembleAtt(omega, sigma);
  auto [Annr_local, Anni_local] = AssembleAnn(omega);

  // Compute the shifted (2,1) block: -sigma * Btn_r. For sigma = -kn_target^2, this
  // gives +kn_target^2 * Btn_r. Btn is real-only (no imaginary part).
  std::unique_ptr<mfem::HypreParMatrix> shifted_Btnr;
  if (Btnr && std::abs(sigma) > 0.0)
  {
    shifted_Btnr.reset(mfem::Add(-sigma, *Btnr, 0.0, *Btnr));
  }

  auto [Ar, Ai] = BuildSystemMatrixA(Attr.get(), Atti.get(), Atnr.get(), Atni.get(),
                                     Annr_local.get(), Anni_local.get(),
                                     shifted_Btnr.get());
  opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

NewBoundaryModeSolver::SolveResult
NewBoundaryModeSolver::Solve(double omega, double sigma, const ComplexVector *initial_space)
{
  // Assemble frequency-dependent matrices (MPI collective on FE space communicator).
  AssembleFrequencyDependent(omega, sigma);

  // Solve the eigenvalue problem. All processes must have solvers for this path.
  MFEM_VERIFY(ksp && eigen,
              "NewBoundaryModeSolver::Solve called on process without solvers! "
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

NewBoundaryModeSolver::SolveResult
NewBoundaryModeSolver::SolveSplit(double omega, double sigma, bool has_solver,
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

std::complex<double> NewBoundaryModeSolver::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(i);
}

void NewBoundaryModeSolver::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(i, x);
}

// Stiffness matrix (shifted): Att = (mu_cc^{-1} curl_t x u, curl_t x v)
//                                   - omega^2 (eps u, v)
//                                   - sigma (mu^{-1} u, v)
//                                   + BC-t impedance/absorbing/conductivity.
// SAME as BoundaryModeSolver::AssembleAtt.
NewBoundaryModeSolver::ComplexHypreParMatrix
NewBoundaryModeSolver::AssembleAtt(double omega, double sigma) const
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
  // London superconductor contribution: +(1/lambda_L^2)(et, ft).
  if (config.has_london_depth)
  {
    eps_shifted_func.AddCoefficient(*config.attr_to_material, *config.inv_london_depth,
                                    1.0);
  }

  // Assemble Att real part: domain + boundary contributions in a single BilinearForm,
  // following the same pattern as SpaceOperator. Boundary impedance/absorbing terms use
  // Palace's ceed-based BilinearForm with VectorFEMassIntegrator.
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);

  // Add boundary contributions from impedance, absorbing, and conductivity operators.
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);        // Ls -> real
    config.surf_z_op->AddDampingBdrCoefficients(omega, fbi);        // Rs -> imag
    config.surf_z_op->AddMassBdrCoefficients(-omega * omega, fbr);  // Cs -> real
  }
  if (config.farfield_op)
  {
    config.farfield_op->AddDampingBdrCoefficients(omega, fbi);  // 1st-order ABC
  }
  if (config.surf_sigma_op)
  {
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  }

  Mpi::Print(" [NBMS] Att boundary: surf_z_op={}, farfield_op={}, surf_sigma_op={}, "
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
    Mpi::Print(" [NBMS] Attr Frobenius norm = {:.6e}\n", norm_r);
  }

  // Assemble imaginary part: domain (loss tangent, conductivity) + boundary (Rs, ABC).
  std::unique_ptr<mfem::HypreParMatrix> Atti_assembled;
  {
    bool has_imag = config.has_loss_tangent ||
                    (config.has_conductivity && config.conductivity) || !fbi.empty();
    if (has_imag)
    {
      // Coefficients must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = config.attr_to_material->Size();
      MaterialPropertyCoefficient negepstandelta_func(n_attr);
      MaterialPropertyCoefficient fi_domain(n_attr);
      if (config.has_loss_tangent)
      {
        negepstandelta_func.AddCoefficient(*config.attr_to_material,
                                           *config.permittivity_imag, -omega * omega);
      }
      if (config.has_conductivity && config.conductivity)
      {
        fi_domain.AddCoefficient(*config.attr_to_material, *config.conductivity, omega);
      }
      BilinearForm atti(nd_fespace);
      if (!negepstandelta_func.empty())
      {
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
      }
      if (!fi_domain.empty())
      {
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
// SAME as BoundaryModeSolver::AssembleAtn.
NewBoundaryModeSolver::ComplexHypreParMatrix NewBoundaryModeSolver::AssembleAtn() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material, *config.inv_permeability,
                                         -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

// H1 block: Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n impedance.
//
// This replaces BoundaryModeSolver's Ann (which is just -eps mass) with the full normal
// curl-curl equation from Equation 2. The stiffness (diffusion) is NEGATIVE (IBP of
// Laplacian), the mass is POSITIVE (+omega^2 eps), and impedance boundary terms use
// NEGATED signs relative to the tangential Att terms (from the double cross product
// in the impedance BC: [n x (n x E)]_z = -En).
NewBoundaryModeSolver::ComplexHypreParMatrix
NewBoundaryModeSolver::AssembleAnn(double omega) const
{
  // H1 stiffness: -(mu^{-1} grad u, grad v) -- negative sign from IBP of Laplacian.
  MaterialPropertyCoefficient neg_muinv_func(*config.attr_to_material,
                                             *config.inv_permeability, -1.0);
  if (config.normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*config.normal);
  }

  // H1 mass: +omega^2(eps u, v) -- positive sign in the normal equation.
  MaterialPropertyCoefficient poseps_h1_func(*config.attr_to_material,
                                             *config.permittivity_scalar, omega * omega);
  if (config.normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*config.normal);
  }

  // London superconductor contribution: +(1/lambda_L^2)(en, fn).
  // In the normal curl-curl equation, the London term is +(1/lambda_L^2) En, which
  // enters as a positive H1 mass (same sign as the omega^2 eps mass).
  // The inv_london_depth tensor is NxN for ND; for H1 we need scalar (0,0) component.
  if (config.has_london_depth)
  {
    const auto &ild = *config.inv_london_depth;
    mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
    for (int k = 0; k < ild.SizeK(); k++)
    {
      ild_scalar(0, 0, k) = ild(0, 0, k);
    }
    poseps_h1_func.AddCoefficient(*config.attr_to_material, ild_scalar);
  }

  // Boundary impedance for en: -(iw/Zs)(en, fn)_gamma -- NEGATIVE sign from the
  // double cross product in the impedance BC: [n x (n x E)]_z = -En, giving
  // -(iw/Zs)(En, fn) on the boundary. Note: opposite sign from et boundary term.
  // This follows the same pattern as QuadBoundaryModeSolver's M0(2,2) assembly.
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr), nn_fbi(max_bdr_attr);
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    config.surf_z_op->AddDampingBdrCoefficients(-omega, nn_fbi);
    config.surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (config.farfield_op && config.farfield_op->GetAttrList().Size() > 0)
  {
    // The farfield operator's AddDampingBdrCoefficients adds a tensor-valued inverse
    // impedance (for VectorFEMassIntegrator on ND). For the scalar H1 MassIntegrator,
    // we need the scalar inverse impedance sqrt(eps/mu). Extract the (0,0) component
    // of the tensor for each boundary material and add as scalar coefficients.
    const auto &farfield_attrs = config.farfield_op->GetAttrList();
    const auto &inv_z = config.mat_op->GetInvImpedance();
    const auto &bdr_attr_to_mat = config.mat_op->GetBdrAttributeToMaterial();
    for (auto attr : farfield_attrs)
    {
      int mat_idx =
          (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
      double inv_z0_scalar = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
      auto ceed_attrs = config.mat_op->GetCeedBdrAttributes(attr);
      if (ceed_attrs.Size() > 0)
      {
        nn_fbi.AddMaterialProperty(ceed_attrs, inv_z0_scalar, -omega);
      }
    }
  }
  if (config.surf_sigma_op)
  {
    // For conductivity: negate the coefficients (opposite sign from et boundary term).
    MaterialPropertyCoefficient cond_r(max_bdr_attr), cond_i(max_bdr_attr);
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
    if (!cond_i.empty())
    {
      cond_i *= -1.0;
      nn_fbi.AddCoefficient(cond_i.GetAttributeToMaterial(),
                            cond_i.GetMaterialProperties());
    }
  }

  Mpi::Print(" [NBMS] Ann boundary: nn_fbr.empty={}, nn_fbi.empty={}\n", nn_fbr.empty(),
             nn_fbi.empty());

  // Assemble Ann real: H1 stiffness (negative Laplacian) + mass (positive) + boundary.
  // Use DiffusionMassIntegrator to combine stiffness and mass in a single form.
  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
  if (!nn_fbr.empty())
  {
    annr.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
  }
  auto Annr_assembled =
      ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  // Assemble Ann imaginary: loss tangent + boundary impedance.
  std::unique_ptr<mfem::HypreParMatrix> Anni_assembled;
  {
    bool has_imag = config.has_loss_tangent || !nn_fbi.empty();
    if (has_imag)
    {
      // Coefficient must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = config.attr_to_material->Size();
      MaterialPropertyCoefficient posepsi_h1_func(n_attr);
      if (config.has_loss_tangent)
      {
        posepsi_h1_func.AddCoefficient(*config.attr_to_material,
                                       *config.permittivity_imag, omega * omega);
        if (config.normal)
        {
          posepsi_h1_func.NormalProjectedCoefficient(*config.normal);
        }
      }
      BilinearForm anni(h1_fespace);
      if (!posepsi_h1_func.empty())
      {
        anni.AddDomainIntegrator<MassIntegrator>(posepsi_h1_func);
      }
      if (!nn_fbi.empty())
      {
        anni.AddBoundaryIntegrator<MassIntegrator>(nn_fbi);
      }
      Anni_assembled =
          ParOperator(anni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
    }
  }

  return {std::move(Annr_assembled), std::move(Anni_assembled)};
}

// Mass matrix: Btt = -(mu^{-1} u, v).
// SAME as BoundaryModeSolver::AssembleBtt.
NewBoundaryModeSolver::ComplexHypreParMatrix NewBoundaryModeSolver::AssembleBtt() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material,
                                         *config.inv_permeability);
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

NewBoundaryModeSolver::ComplexHypreParMatrix NewBoundaryModeSolver::BuildSystemMatrixA(
    const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
    const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
    const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
    const mfem::HypreParMatrix *shifted_Btnr) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  // The (1,0) block is -sigma * Btn from the shift-and-invert transformation.
  // Without shift (sigma=0), this block is zero (upper block-triangular).
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Attr;
  blocks(0, 1) = Atnr;
  blocks(1, 0) = shifted_Btnr;  // -sigma * Btn (nullptr when sigma=0)
  blocks(1, 1) = Annr;
  std::unique_ptr<mfem::HypreParMatrix> Ar(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti || Atni || Anni)
  {
    blocks(0, 0) = Atti;
    blocks(0, 1) = Atni;
    blocks(1, 0) = nullptr;  // Shifted Btn is real-only
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

NewBoundaryModeSolver::ComplexHypreParMatrix
NewBoundaryModeSolver::BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                          const mfem::HypreParMatrix *Btti,
                                          const mfem::HypreParMatrix *Btnr,
                                          const mfem::HypreParMatrix *Btni,
                                          const mfem::HypreParMatrix *Dnn) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  // B = [Btt, 0; Btn, Dnn] where Btn = Atn^T and Dnn is the zero diagonal.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Bttr;
  blocks(0, 1) = nullptr;
  blocks(1, 0) = Btnr;
  blocks(1, 1) = Dnn;
  std::unique_ptr<mfem::HypreParMatrix> Br(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Bi;
  if (Btti || Btni)
  {
    blocks(0, 0) = Btti;
    blocks(0, 1) = nullptr;
    blocks(1, 0) = Btni;
    blocks(1, 1) = nullptr;
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

void NewBoundaryModeSolver::SetUpLinearSolver(MPI_Comm comm)
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
    MFEM_ABORT("NewBoundaryModeSolver requires building with SuperLU_DIST, STRUMPACK, "
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

void NewBoundaryModeSolver::SetUpEigenSolver(MPI_Comm comm)
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
#error "NewBoundaryModeSolver requires building with ARPACK or SLEPc!"
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
}

}  // namespace palace
