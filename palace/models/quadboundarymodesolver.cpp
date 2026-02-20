// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "quadboundarymodesolver.hpp"

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

QuadBoundaryModeSolver::QuadBoundaryModeSolver(const BoundaryModeSolverConfig &config,
                                               const FiniteElementSpace &nd_fespace,
                                               const FiniteElementSpace &h1_fespace,
                                               const mfem::Array<int> &dbc_tdof_list,
                                               double omega, MPI_Comm solver_comm)
  : config(config), nd_fespace(nd_fespace), h1_fespace(h1_fespace),
    dbc_tdof_list(dbc_tdof_list)
{
  nd_size = nd_fespace.GetTrueVSize();
  h1_size = h1_fespace.GetTrueVSize();

  // Assemble the three coefficient matrices for the QEP:
  //   (M0 + lambda * M1 + lambda^2 * M2) [et; en] = 0,  lambda = ikn
  auto [M0r, M0i] = AssembleM0(omega);
  opM0 = std::make_unique<ComplexWrapperOperator>(std::move(M0r), std::move(M0i));

  auto [M1r, M1i] = AssembleM1();
  opM1 = std::make_unique<ComplexWrapperOperator>(std::move(M1r), std::move(M1i));

  auto [M2r, M2i] = AssembleM2();
  opM2 = std::make_unique<ComplexWrapperOperator>(std::move(M2r), std::move(M2i));

  // Configure solvers.
  if (solver_comm != MPI_COMM_NULL)
  {
    SetUpLinearSolver(solver_comm);
    SetUpEigenSolver(solver_comm);
  }
  else if (nd_size > 0)
  {
    SetUpLinearSolver(nd_fespace.GetComm());
    SetUpEigenSolver(nd_fespace.GetComm());
  }
}

QuadBoundaryModeSolver::~QuadBoundaryModeSolver()
{
  ksp.reset();
  eigen.reset();
}

QuadBoundaryModeSolver::SolveResult QuadBoundaryModeSolver::Solve(double kn_target)
{
  MFEM_VERIFY(ksp && eigen,
              "QuadBoundaryModeSolver::Solve called on process without solvers!");

  // The PEP eigenvalue is lambda = ikn. Target: lambda_target = i * kn_target.
  std::complex<double> lambda_target(0.0, kn_target);
  double kn2 = kn_target * kn_target;

  // Build the shifted polynomial P(sigma) = M0 + sigma*M1 + sigma^2*M2 for the KSP.
  // With sigma = i*kn_target: P = M0 + ikn*M1 - kn^2*M2.
  //   P_real = M0_real - kn^2 * M2_real
  //   P_imag = M0_imag + kn * M1_real
  {
    const auto *M0r = dynamic_cast<const mfem::HypreParMatrix *>(opM0->Real());
    const auto *M0i = dynamic_cast<const mfem::HypreParMatrix *>(opM0->Imag());
    const auto *M1r = dynamic_cast<const mfem::HypreParMatrix *>(opM1->Real());
    const auto *M2r = dynamic_cast<const mfem::HypreParMatrix *>(opM2->Real());

    // P_real = M0_real - kn^2 * M2_real.
    std::unique_ptr<mfem::HypreParMatrix> Pr(mfem::Add(1.0, *M0r, -kn2, *M2r));

    // P_imag = M0_imag + kn * M1_real (M0_imag may be null).
    std::unique_ptr<mfem::HypreParMatrix> Pi;
    if (M0i && M1r)
    {
      Pi.reset(mfem::Add(1.0, *M0i, kn_target, *M1r));
    }
    else if (M1r)
    {
      Pi.reset(mfem::Add(kn_target, *M1r, 0.0, *M1r));  // kn * M1r + 0
    }
    else if (M0i)
    {
      Pi = std::make_unique<mfem::HypreParMatrix>(*M0i);
    }

    auto opP = std::make_unique<ComplexWrapperOperator>(std::move(Pr), std::move(Pi));
    ksp->SetOperators(*opP, *opP);  // Full complex P(sigma) as preconditioner

    // The KSP now owns the reference to opP, but opP goes out of scope. We need to keep
    // it alive. Store it as a member temporarily.
    opP_shifted = std::move(opP);
  }

  // Set PEP operators: M0 + lambda * M1 + lambda^2 * M2.
  eigen->SetOperators(*opM0, *opM1, *opM2, EigenvalueSolver::ScaleType::NONE);
  eigen->SetShiftInvert(lambda_target);

  int num_conv = eigen->Solve();
  Mpi::Print(" Found {:d} converged eigenvalue{} (first lambda = {:.3e}{:+.3e}i)\n",
             num_conv, (num_conv != 1) ? "s" : "",
             (num_conv > 0) ? GetEigenvalue(0).real() : 0.0,
             (num_conv > 0) ? GetEigenvalue(0).imag() : 0.0);

  return {num_conv, lambda_target};
}

std::complex<double> QuadBoundaryModeSolver::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(i);
}

void QuadBoundaryModeSolver::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(i, x);
}

std::complex<double> QuadBoundaryModeSolver::GetKn(int i) const
{
  // lambda = ikn  =>  kn = lambda / i = -i * lambda
  std::complex<double> lambda = GetEigenvalue(i);
  return std::complex<double>(lambda.imag(), -lambda.real());
}

// ---------------------------------------------------------------------------
// M0: frequency-dependent block matrix
//
// M0(1,1) = mu_cc^{-1}(curl_t u, curl_t v) - omega^2 (eps u, v)
//           + (iw/Zs)(u_t, v_t)_gamma                          [ND]
// M0(2,2) = mu^{-1}(grad_t u, grad_t v) - omega^2 (eps u, v)
//           + (iw/Zs)(u, v)_gamma                                [H1]
// M0(1,2) = M0(2,1) = 0
// ---------------------------------------------------------------------------
QuadBoundaryModeSolver::ComplexHypreParMatrix
QuadBoundaryModeSolver::AssembleM0(double omega) const
{
  // --- M0(1,1): ND curl-curl + mass + impedance boundary ---
  // Domain: mu_cc^{-1}(curl u, curl v) - omega^2 (eps u, v)
  MaterialPropertyCoefficient muinv_cc_func(*config.attr_to_material,
                                            *config.curlcurl_inv_permeability);
  if (config.normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*config.normal);
  }
  MaterialPropertyCoefficient negeps_func(*config.attr_to_material,
                                          *config.permittivity_real, -omega * omega);
  // London superconductor contribution: +(1/lambda_L^2)(et, ft).
  if (config.has_london_depth)
  {
    negeps_func.AddCoefficient(*config.attr_to_material, *config.inv_london_depth, 1.0);
  }

  // Boundary impedance for et: (iw/Zs)(et.t, ft.t)_gamma
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient tt_fbr(max_bdr_attr), tt_fbi(max_bdr_attr);
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(1.0, tt_fbr);
    config.surf_z_op->AddDampingBdrCoefficients(omega, tt_fbi);
    config.surf_z_op->AddMassBdrCoefficients(-omega * omega, tt_fbr);
  }
  if (config.farfield_op)
  {
    config.farfield_op->AddDampingBdrCoefficients(omega, tt_fbi);
  }
  if (config.surf_sigma_op)
  {
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, tt_fbr, tt_fbi);
  }

  // Assemble M0(1,1) real: domain + boundary.
  BilinearForm m0tt(nd_fespace);
  m0tt.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, negeps_func);
  if (!tt_fbr.empty())
  {
    m0tt.AddBoundaryIntegrator<VectorFEMassIntegrator>(tt_fbr);
  }
  auto M0ttr =
      ParOperator(m0tt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  // Assemble M0(1,1) imaginary: loss tangent + domain conductivity + boundary.
  HyprePtr M0tti;
  {
    bool has_imag = config.has_loss_tangent ||
                    (config.has_conductivity && config.conductivity) || !tt_fbi.empty();
    if (has_imag)
    {
      BilinearForm m0tti(nd_fespace);
      if (config.has_loss_tangent)
      {
        MaterialPropertyCoefficient negepsi_func(*config.attr_to_material,
                                                 *config.permittivity_imag, -omega * omega);
        m0tti.AddDomainIntegrator<VectorFEMassIntegrator>(negepsi_func);
      }
      if (config.has_conductivity && config.conductivity)
      {
        MaterialPropertyCoefficient cond_func(*config.attr_to_material,
                                              *config.conductivity, omega);
        m0tti.AddDomainIntegrator<VectorFEMassIntegrator>(cond_func);
      }
      if (!tt_fbi.empty())
      {
        m0tti.AddBoundaryIntegrator<VectorFEMassIntegrator>(tt_fbi);
      }
      M0tti =
          ParOperator(m0tti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
    }
  }

  // --- M0(2,2): H1 stiffness + mass + impedance boundary ---
  // From the z-component of the curl-curl equation tested with fn (IBP of Laplacian):
  //   -(mu^{-1} grad en, grad fn) + omega^2 (eps en, fn) + impedance_bdr = 0
  // Note: stiffness is NEGATIVE (IBP of Laplacian), mass is POSITIVE.
  MaterialPropertyCoefficient neg_muinv_func(*config.attr_to_material,
                                             *config.inv_permeability, -1.0);
  if (config.normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*config.normal);
  }
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
    // Extract scalar (0,0) component from the tensor for each material.
    mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
    for (int k = 0; k < ild.SizeK(); k++)
    {
      ild_scalar(0, 0, k) = ild(0, 0, k);
    }
    poseps_h1_func.AddCoefficient(*config.attr_to_material, ild_scalar);
  }

  // Boundary impedance for en: -(iw/Zs)(en, fn)_gamma — NEGATIVE sign from the
  // double cross product in the impedance BC: [n x (n x E)]_z = -En, giving
  // -(iw/Zs)(En, fn) on the boundary. Note: opposite sign from et boundary term.
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
    // For conductivity: negate the coefficients.
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

  Mpi::Print(
      " [QEP] M0: tt_fbr.empty={}, tt_fbi.empty={}, nn_fbr.empty={}, nn_fbi.empty={}\n",
      tt_fbr.empty(), tt_fbi.empty(), nn_fbr.empty(), nn_fbi.empty());

  // Assemble M0(2,2) real: H1 stiffness (negative Laplacian) + mass (positive) + boundary.
  BilinearForm m0nn(h1_fespace);
  m0nn.AddDomainIntegrator<DiffusionIntegrator>(neg_muinv_func);
  m0nn.AddDomainIntegrator<MassIntegrator>(poseps_h1_func);
  if (!nn_fbr.empty())
  {
    m0nn.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
  }
  auto M0nnr =
      ParOperator(m0nn.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  // Assemble M0(2,2) imaginary: loss tangent + boundary.
  HyprePtr M0nni;
  {
    bool has_imag = config.has_loss_tangent || !nn_fbi.empty();
    if (has_imag)
    {
      BilinearForm m0nni(h1_fespace);
      if (config.has_loss_tangent)
      {
        // Loss tangent contribution: +omega^2 * eps_imag (same positive sign as real mass).
        MaterialPropertyCoefficient posepsi_h1_func(
            *config.attr_to_material, *config.permittivity_imag, omega * omega);
        if (config.normal)
        {
          posepsi_h1_func.NormalProjectedCoefficient(*config.normal);
        }
        m0nni.AddDomainIntegrator<MassIntegrator>(posepsi_h1_func);
      }
      if (!nn_fbi.empty())
      {
        m0nni.AddBoundaryIntegrator<MassIntegrator>(nn_fbi);
      }
      M0nni =
          ParOperator(m0nni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
    }
  }

  // Build 2x2 block M0 = [M0tt, 0; 0, M0nn] with Dirichlet BC elimination.
  return BuildBlockMatrix(M0ttr.get(), nullptr, nullptr, M0nnr.get(), M0tti.get(), nullptr,
                          nullptr, M0nni.get(), Operator::DIAG_ONE);
}

// ---------------------------------------------------------------------------
// M1: frequency-independent coupling
//
// M1(1,2) = -mu^{-1}(grad_t u, v)  [MixedVectorGrad, H1→ND]
// M1(2,1) = -mu^{-1}(u, grad_t v)  [transpose of M1(1,2)]
// M1(1,1) = M1(2,2) = 0
// ---------------------------------------------------------------------------
QuadBoundaryModeSolver::ComplexHypreParMatrix QuadBoundaryModeSolver::AssembleM1() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material, *config.inv_permeability,
                                         -1.0);
  BilinearForm grad(h1_fespace, nd_fespace);
  grad.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  auto Gr = ParOperator(grad.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
                .StealParallelAssemble();

  // M1(2,1) = transpose of M1(1,2).
  auto GrT = std::unique_ptr<mfem::HypreParMatrix>(Gr->Transpose());

  // Build 2x2 block M1 = [0, G; G^T, 0] with Dirichlet BC elimination.
  return BuildBlockMatrix(nullptr, Gr.get(), GrT.get(), nullptr, nullptr, nullptr, nullptr,
                          nullptr, Operator::DIAG_ZERO);
}

// ---------------------------------------------------------------------------
// M2: frequency-independent ND mass
//
// M2(1,1) = -mu^{-1}(u, v)  [ND mass, negative sign from lambda^2 = -kn^2]
// M2(1,2) = M2(2,1) = M2(2,2) = 0
// ---------------------------------------------------------------------------
QuadBoundaryModeSolver::ComplexHypreParMatrix QuadBoundaryModeSolver::AssembleM2()
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material, *config.inv_permeability,
                                         -1.0);
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  auto Bttr_local =
      ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  // Save a POSITIVE copy of Btt for postprocessing (impedance, power normalization).
  // The M2 assembly uses -mu^{-1}, but postprocessing needs +mu^{-1}.
  {
    MaterialPropertyCoefficient muinv_pos(*config.attr_to_material,
                                          *config.inv_permeability);
    BilinearForm btt_pos(nd_fespace);
    btt_pos.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_pos);
    Bttr =
        ParOperator(btt_pos.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
  }

  // M2(2,2) = 0: construct zero H1 diagonal.
  Vector d(h1_size);
  d.UseDevice(false);
  d = 0.0;
  mfem::SparseMatrix diag(d);
  auto Dnn = std::make_unique<mfem::HypreParMatrix>(
      h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
      h1_fespace.Get().GetTrueDofOffsets(), &diag);

  // Build 2x2 block M2 = [Btt, 0; 0, 0] with Dirichlet BC elimination.
  return BuildBlockMatrix(Bttr_local.get(), nullptr, nullptr, Dnn.get(), nullptr, nullptr,
                          nullptr, nullptr, Operator::DIAG_ZERO);
}

// ---------------------------------------------------------------------------
// Build a 2x2 block HypreParMatrix and eliminate Dirichlet BCs.
// ---------------------------------------------------------------------------
QuadBoundaryModeSolver::ComplexHypreParMatrix QuadBoundaryModeSolver::BuildBlockMatrix(
    const mfem::HypreParMatrix *block00r, const mfem::HypreParMatrix *block01r,
    const mfem::HypreParMatrix *block10r, const mfem::HypreParMatrix *block11r,
    const mfem::HypreParMatrix *block00i, const mfem::HypreParMatrix *block01i,
    const mfem::HypreParMatrix *block10i, const mfem::HypreParMatrix *block11i,
    Operator::DiagonalPolicy diag_policy) const
{
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = block00r;
  blocks(0, 1) = block01r;
  blocks(1, 0) = block10r;
  blocks(1, 1) = block11r;
  auto Ar = HyprePtr(mfem::HypreParMatrixFromBlocks(blocks));

  HyprePtr Ai;
  if (block00i || block01i || block10i || block11i)
  {
    blocks(0, 0) = block00i;
    blocks(0, 1) = block01i;
    blocks(1, 0) = block10i;
    blocks(1, 1) = block11i;
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  Ar->EliminateBC(dbc_tdof_list, diag_policy);
  if (Ai)
  {
    Ai->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Ar), std::move(Ai)};
}

// ---------------------------------------------------------------------------
// Linear solver: GMRES + sparse direct preconditioner (same as BoundaryModeSolver).
// ---------------------------------------------------------------------------
void QuadBoundaryModeSolver::SetUpLinearSolver(MPI_Comm comm)
{
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, config.verbose);
  gmres->SetInitialGuess(false);
  gmres->SetRelTol(config.linear->tol);
  gmres->SetMaxIter(config.linear->max_it);
  gmres->SetRestartDim(config.linear->max_it);

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
    MFEM_ABORT("QuadBoundaryModeSolver requires SuperLU_DIST, STRUMPACK, or MUMPS!");
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
        MFEM_ABORT("Unsupported linear solver type!");
        return {};
      }());
  pc->SetSaveAssembled(false);
  pc->SetDropSmallEntries(false);
  pc->SetComplexMatrix(linear.complex_coarse_solve);
  ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
}

// ---------------------------------------------------------------------------
// Eigenvalue solver: SLEPc PEP (TOAR or linearized Krylov-Schur).
// ---------------------------------------------------------------------------
void QuadBoundaryModeSolver::SetUpEigenSolver(MPI_Comm comm)
{
  constexpr int print = 0;
  EigenSolverBackend type = config.eigen_backend;

  if (type == EigenSolverBackend::SLEPC)
  {
#if !defined(PALACE_WITH_SLEPC)
    MFEM_ABORT("Solver was not built with SLEPc support!");
#endif
  }
  else if (type == EigenSolverBackend::ARPACK)
  {
#if !defined(PALACE_WITH_ARPACK)
    MFEM_ABORT("Solver was not built with ARPACK support!");
#endif
  }
  else
  {
#if defined(PALACE_WITH_SLEPC)
    type = EigenSolverBackend::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
    type = EigenSolverBackend::ARPACK;
#else
#error "QuadBoundaryModeSolver requires ARPACK or SLEPc!"
#endif
  }

  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    auto arpack = std::make_unique<arpack::ArpackPEPSolver>(comm, print);
    arpack->SetNumModes(config.num_modes, config.num_vec);
    arpack->SetTol(config.eig_tol);
    arpack->SetWhichEigenpairs(config.which_eig);
    arpack->SetLinearSolver(*ksp);
    eigen = std::move(arpack);
#endif
  }
  else
  {
#if defined(PALACE_WITH_SLEPC)
    // Use linearized PEP solver (companion linearization + Krylov-Schur).
    auto slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(comm, print);
    slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetNumModes(config.num_modes, config.num_vec);
    slepc->SetTol(config.eig_tol);
    slepc->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    slepc->SetLinearSolver(*ksp);
    eigen = std::move(slepc);
#endif
  }
}

}  // namespace palace
