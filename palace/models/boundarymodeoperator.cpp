// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include <algorithm>
#include <numeric>

#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
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

BoundaryModeOperator::BoundaryModeOperator(const BoundaryModeOperatorConfig &config,
                                           const FiniteElementSpace &nd_fespace,
                                           const FiniteElementSpace &h1_fespace,
                                           const mfem::Array<int> &dbc_tdof_list,
                                           MPI_Comm solver_comm,
                                           const BoundaryModeMultigridConfig *mg_config)
  : config(config), nd_fespace(nd_fespace), h1_fespace(h1_fespace),
    dbc_tdof_list(dbc_tdof_list), mg_config(mg_config)
{
  nd_size = nd_fespace.GetTrueVSize();
  h1_size = h1_fespace.GetTrueVSize();

  // Assemble frequency-independent matrices. These use the FE space communicator
  // (all processes that share the mesh), NOT solver_comm.
  //
  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  std::tie(Atnr, Atni) = AssembleAtn();

  // Btn: NEGATED transpose of Atn. Atn = -(mu^{-1} grad ·, ·), so Atn^T = -(mu^{-1} ·,
  // grad·). The physical Btn from Eq 2 is +(mu^{-1} ·, grad·) (positive), so Btn = -Atn^T.
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

    auto [Br, Bi] =
        BuildSystemMatrixB(bttr.get(), btti.get(), Btnr.get(), Btni.get(), Dnn.get());
    opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));
  }

  // Configure linear and eigenvalue solvers.
  // - BoundaryMode: solver_comm == MPI_COMM_NULL -> use FE space communicator (all procs).
  // - WavePort on port procs: solver_comm is a valid subset communicator -> use it.
  // - WavePort on non-port procs: solver_comm == MPI_COMM_NULL but FE space has no local
  //   DOFs -> skip solver setup (ksp and eigen remain null).
  const bool use_mg =
      mg_config && mg_config->nd_fespaces && mg_config->nd_fespaces->GetNumLevels() > 1;
  if (solver_comm != MPI_COMM_NULL)
  {
    if (use_mg)
    {
      SetUpMultigridLinearSolver(solver_comm);
    }
    else
    {
      SetUpLinearSolver(solver_comm);
    }
    SetUpEigenSolver(solver_comm);
  }
  else if (nd_size > 0)
  {
    // Standalone mode (BoundaryMode): all procs have DOFs, use FE space comm.
    if (use_mg)
    {
      SetUpMultigridLinearSolver(nd_fespace.GetComm());
    }
    else
    {
      SetUpLinearSolver(nd_fespace.GetComm());
    }
    SetUpEigenSolver(nd_fespace.GetComm());
  }
  // else: non-port process with empty FE space -- no solvers needed.
}

BoundaryModeOperator::~BoundaryModeOperator()
{
  // Free the solvers before any communicator they depend on is freed externally.
  ksp.reset();
  eigen.reset();
}

void BoundaryModeOperator::AssembleFrequencyDependent(double omega, double sigma)
{
  auto [Attr, Atti] = AssembleAtt(omega, sigma);
  auto [Annr_local, Anni_local] = AssembleAnn(omega);

  // Compute the shifted (2,1) block: -sigma * Btn_r. For sigma = -kn_target^2, this
  // gives +kn_target^2 * Btn_r. Btn is real-only (no imaginary part).
  std::unique_ptr<mfem::HypreParMatrix> shifted_Btnr;
  if (Btnr && std::abs(sigma) > 0.0)
  {
    shifted_Btnr = std::make_unique<mfem::HypreParMatrix>(*Btnr);
    *shifted_Btnr *= -sigma;
  }

  auto [Ar, Ai] =
      BuildSystemMatrixA(Attr.get(), Atti.get(), Atnr.get(), Atni.get(), Annr_local.get(),
                         Anni_local.get(), shifted_Btnr.get());
  opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

BoundaryModeOperator::SolveResult
BoundaryModeOperator::Solve(double omega, double sigma, bool has_solver,
                            const ComplexVector *initial_space)
{
  // Assemble frequency-dependent matrices (MPI collective on FE space communicator).
  AssembleFrequencyDependent(omega, sigma);

  // Only processes with solvers participate in the eigenvalue solve. For BoundaryMode
  // (has_solver=true on all processes) all processes solve. For WavePort, only port
  // processes solve while non-port processes skip.
  if (!has_solver)
  {
    return {0, sigma};
  }

  MFEM_VERIFY(ksp && eigen,
              "BoundaryModeOperator::Solve called on process without solvers!");

  if (block_pc_ptr)
  {
    // Multigrid path: assemble preconditioner operators at all levels and set on the
    // block-diagonal preconditioner. The outer Krylov solver gets the monolithic opA.
    att_mg_op = AssembleAttPreconditioner(omega, sigma);
    ann_mg_op = AssembleAnnPreconditioner(omega);
    block_pc_ptr->SetBlockOperators(*att_mg_op, *ann_mg_op);

    // Set the off-diagonal operator -sigma*Btn for block lower-triangular preconditioning.
    // This captures the shift-and-invert coupling that dominates the off-diagonal.
    if (Btnr && std::abs(sigma) > 0.0)
    {
      auto sBtnr = std::make_unique<mfem::HypreParMatrix>(*Btnr);
      *sBtnr *= -sigma;
      shifted_Btn_op = std::make_unique<ComplexWrapperOperator>(std::move(sBtnr), nullptr);
      block_pc_ptr->SetOffDiagonalOperator(shifted_Btn_op.get());
    }
    else
    {
      block_pc_ptr->SetOffDiagonalOperator(nullptr);
    }

    ksp->SetOperators(*opA, *opA);  // opA passed twice; pc uses block_pc_ptr
  }
  else
  {
    // Sparse direct path: precondition with real part of the full block system.
    ComplexWrapperOperator opP(opA->Real(), nullptr);
    ksp->SetOperators(*opA, opP);
  }
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);

  if (initial_space)
  {
    eigen->SetInitialSpace(*initial_space);
  }

  int num_conv = eigen->Solve();

  // Build a permutation sorted by descending Re{kn} so that mode ordering is consistent
  // across eigensolver backends (ARPACK vs SLEPc sort eigenvalues differently). Descending
  // order puts the mode closest to the shift target first.
  mode_perm.resize(num_conv);
  std::iota(mode_perm.begin(), mode_perm.end(), 0);
  std::sort(mode_perm.begin(), mode_perm.end(),
            [this, sigma](int a, int b)
            {
              auto kn_a = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(a));
              auto kn_b = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(b));
              return kn_a.real() > kn_b.real();
            });

  return {num_conv, sigma};
}

std::complex<double> BoundaryModeOperator::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(mode_perm[i]);
}

double BoundaryModeOperator::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  return eigen->GetError(mode_perm[i], type);
}

void BoundaryModeOperator::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(mode_perm[i], x);
}

// Stiffness matrix (shifted): Att = (mu_cc^{-1} curl_t x u, curl_t x v)
//                                   - omega^2 (eps u, v)
//                                   - sigma (mu^{-1} u, v)
//                                   + BC-t impedance/absorbing/conductivity.
BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAtt(double omega, double sigma) const
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

  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_func);
  if (!fbr.empty())
  {
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }

  // Assemble domain + boundary contributions via Palace's BilinearForm (libCEED).
  auto Attr_assembled =
      ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
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
BoundaryModeOperator::ComplexHypreParMatrix BoundaryModeOperator::AssembleAtn() const
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
// This replaces BoundaryModeOperator's Ann (which is just -eps mass) with the full normal
// curl-curl equation from Equation 2. The stiffness (diffusion) is NEGATIVE (IBP of
// Laplacian), the mass is POSITIVE (+omega^2 eps), and impedance boundary terms use
// NEGATED signs relative to the tangential Att terms (from the double cross product
// in the impedance BC: [n x (n x E)]_z = -En).
BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAnn(double omega) const
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
  if (config.has_london_depth)
  {
    if (config.inv_london_depth_scalar)
    {
      // 2D mode analysis: use pre-computed scalar out-of-plane London depth.
      poseps_h1_func.AddCoefficient(*config.attr_to_material,
                                    *config.inv_london_depth_scalar);
    }
    else if (config.inv_london_depth && config.normal)
    {
      // 3D wave port: use normal projection of the full tensor (handled by caller via
      // NormalProjectedCoefficient on the assembled form).
      const auto &ild = *config.inv_london_depth;
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(*config.attr_to_material, ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*config.normal);
    }
  }

  // Boundary impedance for en: -(iw/Zs)(en, fn)_gamma -- NEGATIVE sign from the
  // double cross product in the impedance BC: [n x (n x E)]_z = -En, giving
  // -(iw/Zs)(En, fn) on the boundary. Note: opposite sign from et boundary term.
  // This follows the same pattern as the quadratic formulation's M0(2,2) assembly.
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
        if (config.normal)
        {
          // 3D wave port: project the full tensor using the surface normal.
          posepsi_h1_func.AddCoefficient(*config.attr_to_material,
                                         *config.permittivity_imag, omega * omega);
          posepsi_h1_func.NormalProjectedCoefficient(*config.normal);
        }
        else if (config.permittivity_imag_scalar)
        {
          // 2D mode analysis: use pre-computed scalar imaginary permittivity.
          posepsi_h1_func.AddCoefficient(*config.attr_to_material,
                                         *config.permittivity_imag_scalar, omega * omega);
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

// Mass matrix: Btt = (mu^{-1} u, v).
BoundaryModeOperator::ComplexHypreParMatrix BoundaryModeOperator::AssembleBtt() const
{
  MaterialPropertyCoefficient muinv_func(*config.attr_to_material,
                                         *config.inv_permeability);
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

BoundaryModeOperator::ComplexHypreParMatrix BoundaryModeOperator::BuildSystemMatrixA(
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

BoundaryModeOperator::ComplexHypreParMatrix BoundaryModeOperator::BuildSystemMatrixB(
    const mfem::HypreParMatrix *Bttr, const mfem::HypreParMatrix *Btti,
    const mfem::HypreParMatrix *Btnr, const mfem::HypreParMatrix *Btni,
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

void BoundaryModeOperator::SetUpLinearSolver(MPI_Comm comm)
{
  // GMRES iterative solver preconditioned with a sparse direct solver.
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, config.verbose);
  gmres->SetInitialGuess(false);
  gmres->SetRelTol(config.linear->tol);
  gmres->SetMaxIter(config.linear->max_it);
  gmres->SetRestartDim(config.linear->max_it);
  gmres->EnableTimer();

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
    MFEM_ABORT("BoundaryModeOperator requires building with SuperLU_DIST, STRUMPACK, "
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

void BoundaryModeOperator::SetUpMultigridLinearSolver(MPI_Comm comm)
{
  MFEM_VERIFY(mg_config && mg_config->nd_fespaces && mg_config->h1_fespaces &&
                  mg_config->h1_aux_fespaces,
              "Multigrid linear solver requires ND, H1, and H1 auxiliary space "
              "hierarchies!");
  const auto &linear = *config.linear;
  const int print = config.verbose - 1;

  // Determine coarse solver types for the ND and H1 blocks from the config. The default
  // type is resolved in IoData (sparse direct for frequency-domain problems, AMS
  // otherwise). For the H1 block, AMS is not applicable — use BoomerAMG instead.
  LinearSolver nd_pc_type = linear.type;
  LinearSolver h1_pc_type = linear.type;

  if (nd_pc_type == LinearSolver::BOOMER_AMG)
  {
    Mpi::Warning(comm,
                 "BoomerAMG is not well-suited for the Nedelec system matrix, consider "
                 "using another solver.\n");
  }
  if (h1_pc_type == LinearSolver::AMS)
  {
    h1_pc_type = LinearSolver::BOOMER_AMG;
    Mpi::Print(" Multigrid coarse solve: AMS for ND block, BoomerAMG for H1 block\n");
  }

  // Helper to create a coarse solver from the type.
  auto MakeCoarseSolver =
      [&](LinearSolver type) -> std::unique_ptr<MfemWrapperSolver<ComplexOperator>>
  {
    switch (type)
    {
      case LinearSolver::AMS:
        MFEM_VERIFY(mg_config->h1_aux_fespaces,
                    "AMS coarse solver requires auxiliary H1 space!");
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<HypreAmsSolver>(
                mg_config->nd_fespaces->GetFESpaceAtLevel(0),
                mg_config->h1_aux_fespaces->GetFESpaceAtLevel(0), linear.ams_max_it,
                linear.mg_smooth_it, linear.ams_vector_interp, linear.ams_singular_op,
                linear.amg_agg_coarsen, print),
            true, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
      case LinearSolver::BOOMER_AMG:
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<BoomerAmgSolver>(1, linear.mg_smooth_it,
                                              linear.amg_agg_coarsen, print),
            true, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
      case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<SuperLUSolver>(comm, linear.sym_factorization,
                                            linear.superlu_3d, true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with SuperLU_DIST support!");
        return {};
#endif
      case LinearSolver::STRUMPACK:
      case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<StrumpackSolver>(
                comm, linear.sym_factorization, linear.strumpack_compression_type,
                linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
                linear.strumpack_lossy_precision, true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with STRUMPACK support!");
        return {};
#endif
      case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<MumpsSolver>(comm, mfem::MUMPSSolver::UNSYMMETRIC,
                                          linear.sym_factorization, linear.strumpack_lr_tol,
                                          true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with MUMPS support!");
        return {};
#endif
      default:
        MFEM_ABORT("Unsupported coarse solver type for multigrid boundary mode solver!");
        return {};
    }
  };

  // ND block: p-multigrid with Hiptmair distributive relaxation smoothing.
  const auto nd_P = mg_config->nd_fespaces->GetProlongationOperators();
  const auto nd_G =
      mg_config->nd_fespaces->GetDiscreteInterpolators(*mg_config->h1_aux_fespaces);
  auto nd_gmg = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
      comm, MakeCoarseSolver(nd_pc_type), nd_P, &nd_G, linear.mg_cycle_it,
      linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
      linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  nd_gmg->EnableTimer();

  // H1 block: p-multigrid with Chebyshev smoothing.
  const auto h1_P = mg_config->h1_fespaces->GetProlongationOperators();
  auto h1_gmg = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
      comm, MakeCoarseSolver(h1_pc_type), h1_P, nullptr, linear.mg_cycle_it,
      linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
      linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  h1_gmg->EnableTimer();

  // Combine into block-diagonal preconditioner.
  auto block_pc = std::make_unique<BlockDiagonalPreconditioner<ComplexOperator>>(
      nd_size, std::move(nd_gmg), std::move(h1_gmg));
  block_pc_ptr = block_pc.get();

  // Outer Krylov solver — use the user-configured type from Solver.Linear.KSPType.
  std::unique_ptr<IterativeSolver<ComplexOperator>> krylov;
  switch (linear.krylov_solver)
  {
    case KrylovSolver::CG:
      krylov = std::make_unique<CgSolver<ComplexOperator>>(comm, config.verbose);
      break;
    case KrylovSolver::FGMRES:
      {
        auto fgmres = std::make_unique<FgmresSolver<ComplexOperator>>(comm, config.verbose);
        fgmres->SetRestartDim(linear.max_size);
        krylov = std::move(fgmres);
      }
      break;
    case KrylovSolver::GMRES:
    case KrylovSolver::DEFAULT:
    default:
      {
        auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, config.verbose);
        gmres->SetRestartDim(linear.max_size);
        krylov = std::move(gmres);
      }
      break;
  }
  krylov->SetInitialGuess(linear.initial_guess);
  krylov->SetRelTol(linear.tol);
  krylov->SetMaxIter(linear.max_it);
  krylov->EnableTimer();

  ksp = std::make_unique<ComplexKspSolver>(std::move(krylov), std::move(block_pc));
}

std::unique_ptr<ComplexMultigridOperator>
BoundaryModeOperator::AssembleAttPreconditioner(double omega, double sigma) const
{
  MFEM_VERIFY(mg_config && mg_config->nd_fespaces && mg_config->h1_aux_fespaces,
              "AssembleAttPreconditioner requires ND and H1 auxiliary hierarchies!");
  const auto n_levels = mg_config->nd_fespaces->GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients (same at all levels — indexed by element attribute, not by p).
  MaterialPropertyCoefficient muinv_cc_func(*config.attr_to_material,
                                            *config.curlcurl_inv_permeability);
  if (config.normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*config.normal);
  }

  // Preconditioner mass coefficient: -omega^2 eps - sigma/mu (+ London). When
  // pc_mat_shifted is enabled (config), use absolute values to ensure well-conditioning
  // near the shift target (where -omega^2 eps + kn^2/mu ≈ 0). This follows
  // SpaceOperator's pc_mat_shifted pattern.
  const bool shifted = config.linear->pc_mat_shifted;
  const double mass_coeff = shifted ? std::abs(omega * omega) : (-omega * omega);
  const double shift_coeff = shifted ? std::abs(sigma) : (-sigma);

  MaterialPropertyCoefficient eps_shifted_pc(*config.attr_to_material,
                                             *config.permittivity_real, mass_coeff);
  eps_shifted_pc.AddCoefficient(*config.attr_to_material, *config.inv_permeability,
                                shift_coeff);
  if (config.has_london_depth)
  {
    eps_shifted_pc.AddCoefficient(*config.attr_to_material, *config.inv_london_depth, 1.0);
  }

  // Boundary coefficients for preconditioner (real part only).
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient fbr(max_bdr_attr);
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);
    config.surf_z_op->AddMassBdrCoefficients(
        shifted ? std::abs(omega * omega) : (-omega * omega), fbr);
  }
  if (config.surf_sigma_op)
  {
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbr);
  }

  // Assemble ND operators at all levels using the hierarchy-aware assembly pattern
  // (matching SpaceOperator::AssembleOperators). Creates a BilinearForm on the finest
  // space and uses Assemble(fespaces) to produce operators at all levels.
  constexpr bool assemble_q_data = false;
  {
    BilinearForm att(mg_config->nd_fespaces->GetFinestFESpace());
    att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_pc);
    if (!fbr.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
    }
    auto att_ops = att.Assemble(*mg_config->nd_fespaces, skip_zeros);

    // Auxiliary H1 operators for Hiptmair smoothing (matching SpaceOperator::
    // AssembleAuxOperators — uses DiffusionIntegrator with the mass coefficient).
    BilinearForm att_aux(mg_config->h1_aux_fespaces->GetFinestFESpace());
    att_aux.AddDomainIntegrator<DiffusionIntegrator>(eps_shifted_pc);
    if (!fbr.empty())
    {
      att_aux.AddBoundaryIntegrator<DiffusionIntegrator>(fbr);
    }
    auto att_aux_ops = att_aux.Assemble(*mg_config->h1_aux_fespaces, skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &nd_fespace_l = mg_config->nd_fespaces->GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(att_ops[l]), nullptr,
                                                      nd_fespace_l);
      if (mg_config->nd_dbc_tdof_lists && l < mg_config->nd_dbc_tdof_lists->size())
      {
        B_l->SetEssentialTrueDofs((*mg_config->nd_dbc_tdof_lists)[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));

      const auto &h1_aux_l = mg_config->h1_aux_fespaces->GetFESpaceAtLevel(l);
      auto B_aux_l = std::make_unique<ComplexParOperator>(std::move(att_aux_ops[l]),
                                                          nullptr, h1_aux_l);
      if (mg_config->h1_aux_dbc_tdof_lists && l < mg_config->h1_aux_dbc_tdof_lists->size())
      {
        B_aux_l->SetEssentialTrueDofs((*mg_config->h1_aux_dbc_tdof_lists)[l],
                                      Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddAuxiliaryOperator(std::move(B_aux_l));
    }
  }

  return B;
}

std::unique_ptr<ComplexMultigridOperator>
BoundaryModeOperator::AssembleAnnPreconditioner(double omega) const
{
  MFEM_VERIFY(mg_config && mg_config->h1_fespaces,
              "AssembleAnnPreconditioner requires H1 hierarchy!");
  const auto n_levels = mg_config->h1_fespaces->GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients matching AssembleAnn (real part only). The negative diffusion
  // sign from the IBP is preserved — this matches the actual operator.
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
  if (config.has_london_depth)
  {
    if (config.inv_london_depth_scalar)
    {
      poseps_h1_func.AddCoefficient(*config.attr_to_material,
                                    *config.inv_london_depth_scalar);
    }
    else if (config.inv_london_depth && config.normal)
    {
      const auto &ild = *config.inv_london_depth;
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(*config.attr_to_material, ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*config.normal);
    }
  }

  // Boundary coefficients (real part only).
  int max_bdr_attr = config.mat_op ? config.mat_op->MaxCeedBdrAttribute() : 0;
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr);
  if (config.surf_z_op)
  {
    config.surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    config.surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (config.surf_sigma_op)
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr);
    config.surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_r);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
  }

  // Assemble H1 operators at all levels using hierarchy-aware assembly.
  {
    BilinearForm ann(mg_config->h1_fespaces->GetFinestFESpace());
    ann.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
    if (!nn_fbr.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
    }
    auto ann_ops = ann.Assemble(*mg_config->h1_fespaces, skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &h1_fespace_l = mg_config->h1_fespaces->GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(ann_ops[l]), nullptr,
                                                      h1_fespace_l);
      if (mg_config->h1_dbc_tdof_lists && l < mg_config->h1_dbc_tdof_lists->size())
      {
        B_l->SetEssentialTrueDofs((*mg_config->h1_dbc_tdof_lists)[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));
    }
  }

  return B;
}

void BoundaryModeOperator::SetUpEigenSolver(MPI_Comm comm)
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
#error "BoundaryModeOperator requires building with ARPACK or SLEPc!"
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
