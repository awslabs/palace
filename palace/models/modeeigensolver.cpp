// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeeigensolver.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

ModeEigenSolver::ModeEigenSolver(
    const MaterialOperator &mat_op, const mfem::Vector *normal,
    SurfaceImpedanceOperator *surf_z_op, FarfieldBoundaryOperator *farfield_op,
    SurfaceConductivityOperator *surf_sigma_op, const FiniteElementSpace &nd_fespace,
    const FiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_tdof_list,
    int num_modes, int num_vec, double eig_tol, EigenvalueSolver::WhichType which_eig,
    const config::LinearSolverData &linear, EigenSolverBackend eigen_backend, int verbose,
    MPI_Comm solver_comm, BoundaryModeOperator *bmo)
  : num_modes(num_modes), num_vec(num_vec), eig_tol(eig_tol), which_eig(which_eig),
    linear(linear), eigen_backend(eigen_backend), verbose(verbose), mat_op(mat_op),
    normal(normal), surf_z_op(surf_z_op), farfield_op(farfield_op),
    surf_sigma_op(surf_sigma_op), nd_fespace(nd_fespace), h1_fespace(h1_fespace),
    dbc_tdof_list(dbc_tdof_list), bmo(bmo)
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

  // Pick a communicator for the solvers: solver_comm if supplied (WavePort port ranks),
  // else the FE space comm (BoundaryMode, all ranks). Ranks with no DOFs and no
  // solver_comm (WavePort non-port ranks) skip solver setup entirely and return
  // num_converged = 0 from Solve.
  const bool use_mg = bmo && bmo->GetNDSpaceHierarchy().GetNumLevels() > 1;
  MPI_Comm configure_comm = (solver_comm != MPI_COMM_NULL) ? solver_comm
                                                           : (nd_size > 0)
                                                                 ? nd_fespace.GetComm()
                                                                 : MPI_COMM_NULL;
  if (configure_comm != MPI_COMM_NULL)
  {
    if (use_mg)
    {
      SetUpMultigridLinearSolver(configure_comm);
    }
    else
    {
      SetUpLinearSolver(configure_comm);
    }
    SetUpEigenSolver(configure_comm);
  }
}

void ModeEigenSolver::AssembleFrequencyDependent(double omega, double sigma)
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

ModeEigenSolver::SolveResult
ModeEigenSolver::Solve(double omega, double sigma, const ComplexVector *initial_space)
{
  sigma_cached = sigma;

  // Frequency-dependent matrices assemble on the FE space communicator.
  AssembleFrequencyDependent(omega, sigma);

  // Ranks configured without a solver (wave port non-port ranks) return after assembly.
  if (!ksp || !eigen)
  {
    return {0, sigma};
  }

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

  // Build a permutation sorted by proximity to the shift target so that mode ordering is
  // consistent across eigensolver backends (ARPACK vs SLEPc sort eigenvalues differently).
  // The shift sigma = -kn_target^2, so kn_target = sqrt(-sigma). Sorting by ascending
  // |Re{kn} - kn_target| puts the mode closest to the target first.
  const double kn_target = std::sqrt(-sigma);
  mode_perm.resize(num_conv);
  std::iota(mode_perm.begin(), mode_perm.end(), 0);
  std::sort(mode_perm.begin(), mode_perm.end(),
            [this, sigma, kn_target](int a, int b)
            {
              auto kn_a = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(a));
              auto kn_b = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(b));
              return std::abs(kn_a.real() - kn_target) < std::abs(kn_b.real() - kn_target);
            });

  return {num_conv, sigma};
}

std::complex<double> ModeEigenSolver::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(mode_perm[i]);
}

double ModeEigenSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  return eigen->GetError(mode_perm[i], type);
}

void ModeEigenSolver::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(mode_perm[i], x);
}

std::complex<double> ModeEigenSolver::GetPropagationConstant(int i) const
{
  return std::sqrt(-sigma_cached - 1.0 / eigen->GetEigenvalue(mode_perm[i]));
}

void ModeEigenSolver::ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn,
                                           ComplexVector &et, ComplexVector &en) const
{
  et.Real().MakeRef(e0.Real(), 0, nd_size);
  et.Imag().MakeRef(e0.Imag(), 0, nd_size);
  en.Real().MakeRef(e0.Real(), nd_size, h1_size);
  en.Imag().MakeRef(e0.Imag(), nd_size, h1_size);
  const auto ikn_inv = 1.0 / (std::complex<double>(0.0, 1.0) * kn);
  ComplexVector::AXPBY(ikn_inv, en.Real(), en.Imag(), 0.0, en.Real(), en.Imag());
}

std::complex<double> ModeEigenSolver::ComputePoyntingPower(double omega,
                                                           std::complex<double> kn,
                                                           const ComplexVector &et,
                                                           const ComplexVector &en) const
{
  if (!Bttr)
  {
    return 0.0;
  }
  auto comm = nd_fespace.GetComm();
  std::complex<double> P = 0.5 * std::conj(kn) / omega * linalg::Dot(comm, et, *Bttr, et);
  if (Atnr && en.Size() == h1_size)
  {
    ComplexWrapperOperator Atn(const_cast<mfem::HypreParMatrix *>(Atnr.get()),
                               const_cast<mfem::HypreParMatrix *>(Atni.get()));
    P += std::complex<double>(0.0, 1.0) / (2.0 * omega) * linalg::Dot(comm, en, Atn, et);
  }
  return P;
}

bool ModeEigenSolver::IsPropagating(std::complex<double> kn)
{
  return std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
}

// Stiffness matrix (shifted): Att = (mu_cc^{-1} curl_t x u, curl_t x v)
//                                   - omega^2 (eps u, v)
//                                   - sigma (mu^{-1} u, v)
//                                   + BC-t impedance/absorbing/conductivity.
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAtt(double omega,
                                                                    double sigma) const
{
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient eps_shifted_func(
      mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), -omega * omega);
  eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvPermeability(), -sigma);
  // London superconductor contribution: +(1/lambda_L^2)(et, ft).
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepth(), 1.0);
  }

  // Assemble Att real part: domain + boundary contributions in a single BilinearForm,
  // following the same pattern as SpaceOperator. Boundary impedance/absorbing terms use
  // Palace's ceed-based BilinearForm with VectorFEMassIntegrator.
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);

  // Add boundary contributions from impedance, absorbing, and conductivity operators.
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);        // Ls -> real
    surf_z_op->AddDampingBdrCoefficients(omega, fbi);        // Rs -> imag
    surf_z_op->AddMassBdrCoefficients(-omega * omega, fbr);  // Cs -> real
  }
  if (farfield_op)
  {
    farfield_op->AddDampingBdrCoefficients(omega, fbi);  // 1st-order ABC
  }
  if (surf_sigma_op)
  {
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbi);
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
    bool has_imag = mat_op.HasLossTangent() || mat_op.HasConductivity() || !fbi.empty();
    if (has_imag)
    {
      // Coefficients must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient negepstandelta_func(n_attr);
      MaterialPropertyCoefficient fi_domain(n_attr);
      if (mat_op.HasLossTangent())
      {
        negepstandelta_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityImag(), -omega * omega);
      }
      if (mat_op.HasConductivity())
      {
        fi_domain.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(),
                                 omega);
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
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAtn() const
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

// H1 block: Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n impedance.
//
// This replaces ModeEigenSolver's Ann (which is just -eps mass) with the full normal
// curl-curl equation from Equation 2. The stiffness (diffusion) is NEGATIVE (IBP of
// Laplacian), the mass is POSITIVE (+omega^2 eps), and impedance boundary terms use
// NEGATED signs relative to the tangential Att terms (from the double cross product
// in the impedance BC: [n x (n x E)]_z = -En).
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAnn(double omega) const
{
  // H1 stiffness: -(mu^{-1} grad u, grad v) -- negative sign from IBP of Laplacian.
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  // H1 mass: +omega^2(eps u, v) -- positive sign in the normal equation.
  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             normal ? mat_op.GetPermittivityReal()
                                                    : mat_op.GetPermittivityScalar(),
                                             omega * omega);
  if (normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*normal);
  }

  // London superconductor contribution: +(1/lambda_L^2)(en, fn).
  // In the normal curl-curl equation, the London term is +(1/lambda_L^2) En, which
  // enters as a positive H1 mass (same sign as the omega^2 eps mass).
  if (mat_op.HasLondonDepth())
  {
    if (!normal)
    {
      // 2D mode analysis: use pre-computed scalar out-of-plane London depth.
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepthScalar());
    }
    else if (normal)
    {
      // 3D wave port: use normal projection of the full tensor (handled by caller via
      // NormalProjectedCoefficient on the assembled form).
      const auto &ild = mat_op.GetInvLondonDepth();
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(), ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*normal);
    }
  }

  // Boundary impedance for en: -(iw/Zs)(en, fn)_gamma -- NEGATIVE sign from the
  // double cross product in the impedance BC: [n x (n x E)]_z = -En, giving
  // -(iw/Zs)(En, fn) on the boundary. Note: opposite sign from et boundary term.
  // This follows the same pattern as the quadratic formulation's M0(2,2) assembly.
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr), nn_fbi(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    surf_z_op->AddDampingBdrCoefficients(-omega, nn_fbi);
    surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (farfield_op && farfield_op->GetAttrList().Size() > 0)
  {
    // The farfield operator's AddDampingBdrCoefficients adds a tensor-valued inverse
    // impedance (for VectorFEMassIntegrator on ND). For the scalar H1 MassIntegrator,
    // we need the scalar inverse impedance sqrt(eps/mu). Extract the (0,0) component
    // of the tensor for each boundary material and add as scalar coefficients.
    const auto &farfield_attrs = farfield_op->GetAttrList();
    const auto &inv_z = mat_op.GetInvImpedance();
    const auto &bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
    for (auto attr : farfield_attrs)
    {
      int mat_idx =
          (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
      double inv_z0_scalar = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
      auto ceed_attrs = mat_op.GetCeedBdrAttributes(attr);
      if (ceed_attrs.Size() > 0)
      {
        nn_fbi.AddMaterialProperty(ceed_attrs, inv_z0_scalar, -omega);
      }
    }
  }
  if (surf_sigma_op)
  {
    // For conductivity: negate the coefficients (opposite sign from et boundary term).
    MaterialPropertyCoefficient cond_r(max_bdr_attr), cond_i(max_bdr_attr);
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
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
    bool has_imag = mat_op.HasLossTangent() || !nn_fbi.empty();
    if (has_imag)
    {
      // Coefficient must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient posepsi_h1_func(n_attr);
      if (mat_op.HasLossTangent())
      {
        if (normal)
        {
          // 3D wave port: project the full tensor using the surface normal.
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImag(), omega * omega);
          posepsi_h1_func.NormalProjectedCoefficient(*normal);
        }
        else if (!normal)
        {
          // 2D mode analysis: use pre-computed scalar imaginary permittivity.
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImagScalar(), omega * omega);
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
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleBtt() const
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::BuildSystemMatrixA(
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
    // HypreParMatrixFromBlocks requires at least one non-null block per row and column
    // to determine sizes. Since (1,0) is always null (shifted Btn is real-only), add
    // zero diagonal placeholders when an entire block row or column would be null.
    std::unique_ptr<mfem::HypreParMatrix> Dtt_zero, Dnn_zero;
    if (!Atti && !Atni)
    {
      Vector d(nd_size);
      d.UseDevice(false);
      d = 0.0;
      mfem::SparseMatrix diag(d);
      Dtt_zero = std::make_unique<mfem::HypreParMatrix>(
          nd_fespace.Get().GetComm(), nd_fespace.Get().GlobalTrueVSize(),
          nd_fespace.Get().GetTrueDofOffsets(), &diag);
    }
    if (!Anni)
    {
      Vector d(h1_size);
      d.UseDevice(false);
      d = 0.0;
      mfem::SparseMatrix diag(d);
      Dnn_zero = std::make_unique<mfem::HypreParMatrix>(
          h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
          h1_fespace.Get().GetTrueDofOffsets(), &diag);
    }
    blocks(0, 0) = Atti ? Atti : Dtt_zero.get();
    blocks(0, 1) = Atni;
    blocks(1, 0) = nullptr;  // Shifted Btn is real-only
    blocks(1, 1) = Anni ? Anni : Dnn_zero.get();
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

ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::BuildSystemMatrixB(
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
    // NOTE: Currently unreachable (Btt, Btn are real for real permeability). If complex
    // permeability is added, zero placeholder blocks would be needed here too (same as
    // the imaginary A block above) to prevent HypreParMatrixFromBlocks sizing errors.
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

void ModeEigenSolver::SetUpLinearSolver(MPI_Comm comm)
{
  // GMRES iterative solver preconditioned with a sparse direct solver.
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, verbose);
  gmres->SetInitialGuess(false);
  gmres->SetRelTol(linear.tol);
  gmres->SetMaxIter(linear.max_it);
  gmres->SetRestartDim(linear.max_size);
  gmres->EnableTimer();

  // Select sparse direct solver type. Multigrid (AMS, BoomerAMG) is not applicable for
  // the combined ND+H1 block system, so fall back to a sparse direct solver.
  LinearSolver pc_type = linear.type;
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
    MFEM_ABORT("ModeEigenSolver requires building with SuperLU_DIST, STRUMPACK, "
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

  auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
      [&]() -> std::unique_ptr<mfem::Solver>
      {
        if (pc_type == LinearSolver::SUPERLU)
        {
#if defined(MFEM_USE_SUPERLU)
          return std::make_unique<SuperLUSolver>(comm, linear.sym_factorization,
                                                 linear.superlu_3d, true, verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::STRUMPACK ||
                 pc_type == LinearSolver::STRUMPACK_MP)
        {
#if defined(MFEM_USE_STRUMPACK)
          return std::make_unique<StrumpackSolver>(
              comm, linear.sym_factorization, linear.strumpack_compression_type,
              linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
              linear.strumpack_lossy_precision, true, verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::MUMPS)
        {
#if defined(MFEM_USE_MUMPS)
          return std::make_unique<MumpsSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
                                               linear.sym_factorization,
                                               linear.strumpack_lr_tol, true, verbose - 1);
#endif
        }
        MFEM_ABORT("Unsupported linear solver type for boundary mode solver!");
        return {};
      }());
  pc->SetSaveAssembled(false);
  pc->SetDropSmallEntries(false);
  ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
}

void ModeEigenSolver::SetUpMultigridLinearSolver(MPI_Comm comm)
{
  MFEM_VERIFY(bmo,
              "Multigrid linear solver requires a BoundaryModeOperator reference!");
  const int print = verbose - 1;

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
        MFEM_VERIFY(bmo, "AMS coarse solver requires auxiliary H1 space!");
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<HypreAmsSolver>(
                bmo->GetNDSpaceHierarchy().GetFESpaceAtLevel(0),
                bmo->GetH1AuxSpaceHierarchy().GetFESpaceAtLevel(0), linear.ams_max_it,
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
            std::make_unique<MumpsSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
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
  const auto nd_P = bmo->GetNDSpaceHierarchy().GetProlongationOperators();
  const auto nd_G =
      bmo->GetNDSpaceHierarchy().GetDiscreteInterpolators(bmo->GetH1AuxSpaceHierarchy());
  auto nd_gmg = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
      comm, MakeCoarseSolver(nd_pc_type), nd_P, &nd_G, linear.mg_cycle_it,
      linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
      linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  nd_gmg->EnableTimer();

  // H1 block: p-multigrid with Chebyshev smoothing.
  const auto h1_P = bmo->GetH1SpaceHierarchy().GetProlongationOperators();
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
      krylov = std::make_unique<CgSolver<ComplexOperator>>(comm, verbose);
      break;
    case KrylovSolver::FGMRES:
      {
        auto fgmres = std::make_unique<FgmresSolver<ComplexOperator>>(comm, verbose);
        fgmres->SetRestartDim(linear.max_size);
        krylov = std::move(fgmres);
      }
      break;
    case KrylovSolver::GMRES:
    case KrylovSolver::DEFAULT:
    default:
      {
        auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, verbose);
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
ModeEigenSolver::AssembleAttPreconditioner(double omega, double sigma) const
{
  MFEM_VERIFY(bmo, "AssembleAttPreconditioner requires a BoundaryModeOperator!");
  const auto n_levels = bmo->GetNDSpaceHierarchy().GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients (same at all levels — indexed by element attribute, not by p).
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  // Preconditioner mass coefficient: -omega^2 eps - sigma/mu (+ London). When
  // pc_mat_shifted is enabled (config), use absolute values to ensure well-conditioning
  // near the shift target (where -omega^2 eps + kn^2/mu ≈ 0). This follows
  // SpaceOperator's pc_mat_shifted pattern.
  const bool shifted = linear.pc_mat_shifted;
  const double mass_coeff = shifted ? std::abs(omega * omega) : (-omega * omega);
  const double shift_coeff = shifted ? std::abs(sigma) : (-sigma);

  MaterialPropertyCoefficient eps_shifted_pc(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal(), mass_coeff);
  eps_shifted_pc.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                mat_op.GetInvPermeability(), shift_coeff);
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_pc.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvLondonDepth(), 1.0);
  }

  // Boundary coefficients for preconditioner (real part only, matching SpaceOperator).
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);
    surf_z_op->AddMassBdrCoefficients(shifted ? std::abs(omega * omega) : (-omega * omega),
                                      fbr);
  }
  if (farfield_op)
  {
    farfield_op->AddDampingBdrCoefficients(omega, fbr);
  }
  if (surf_sigma_op)
  {
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbr);
  }

  // Assemble ND operators at all levels using the hierarchy-aware assembly pattern
  // (matching SpaceOperator::AssembleOperators). Creates a BilinearForm on the finest
  // space and uses Assemble(fespaces) to produce operators at all levels.
  constexpr bool assemble_q_data = false;
  {
    BilinearForm att(bmo->GetNDSpaceHierarchy().GetFinestFESpace());
    att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_pc);
    if (!fbr.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
    }
    auto att_ops = att.Assemble(bmo->GetNDSpaceHierarchy(), skip_zeros);

    // Auxiliary H1 operators for Hiptmair smoothing (matching SpaceOperator::
    // AssembleAuxOperators — uses DiffusionIntegrator with the mass coefficient).
    BilinearForm att_aux(bmo->GetH1AuxSpaceHierarchy().GetFinestFESpace());
    att_aux.AddDomainIntegrator<DiffusionIntegrator>(eps_shifted_pc);
    if (!fbr.empty())
    {
      att_aux.AddBoundaryIntegrator<DiffusionIntegrator>(fbr);
    }
    auto att_aux_ops = att_aux.Assemble(bmo->GetH1AuxSpaceHierarchy(), skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &nd_fespace_l = bmo->GetNDSpaceHierarchy().GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(att_ops[l]), nullptr,
                                                      nd_fespace_l);
      if (l < bmo->GetNDDbcTDofLists().size())
      {
        B_l->SetEssentialTrueDofs(bmo->GetNDDbcTDofLists()[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));

      const auto &h1_aux_l = bmo->GetH1AuxSpaceHierarchy().GetFESpaceAtLevel(l);
      auto B_aux_l = std::make_unique<ComplexParOperator>(std::move(att_aux_ops[l]),
                                                          nullptr, h1_aux_l);
      if (l < bmo->GetH1AuxDbcTDofLists().size())
      {
        B_aux_l->SetEssentialTrueDofs(bmo->GetH1AuxDbcTDofLists()[l],
                                      Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddAuxiliaryOperator(std::move(B_aux_l));
    }
  }

  return B;
}

std::unique_ptr<ComplexMultigridOperator>
ModeEigenSolver::AssembleAnnPreconditioner(double omega) const
{
  MFEM_VERIFY(bmo, "AssembleAnnPreconditioner requires a BoundaryModeOperator!");
  const auto n_levels = bmo->GetH1SpaceHierarchy().GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients matching AssembleAnn (real part only). The negative diffusion
  // sign from the IBP is preserved — this matches the actual operator.
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             normal ? mat_op.GetPermittivityReal()
                                                    : mat_op.GetPermittivityScalar(),
                                             omega * omega);
  if (normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*normal);
  }
  if (mat_op.HasLondonDepth())
  {
    if (!normal)
    {
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepthScalar());
    }
    else if (normal)
    {
      const auto &ild = mat_op.GetInvLondonDepth();
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(), ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*normal);
    }
  }

  // Boundary coefficients (real part only).
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (surf_sigma_op)
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr);
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_r);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
  }

  // Assemble H1 operators at all levels using hierarchy-aware assembly.
  {
    BilinearForm ann(bmo->GetH1SpaceHierarchy().GetFinestFESpace());
    ann.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
    if (!nn_fbr.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
    }
    auto ann_ops = ann.Assemble(bmo->GetH1SpaceHierarchy(), skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &h1_fespace_l = bmo->GetH1SpaceHierarchy().GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(ann_ops[l]), nullptr,
                                                      h1_fespace_l);
      if (l < bmo->GetH1DbcTDofLists().size())
      {
        B_l->SetEssentialTrueDofs(bmo->GetH1DbcTDofLists()[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));
    }
  }

  return B;
}

void ModeEigenSolver::SetUpEigenSolver(MPI_Comm comm)
{
  constexpr int print = 0;
  EigenSolverBackend type = eigen_backend;

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
#error "ModeEigenSolver requires building with ARPACK or SLEPc!"
#endif
  }

  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    auto arpack = std::make_unique<arpack::ArpackEPSSolver>(comm, print);
    arpack->SetNumModes(num_modes, num_vec);
    arpack->SetTol(eig_tol);
    arpack->SetWhichEigenpairs(which_eig);
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
    slepc->SetNumModes(num_modes, num_vec);
    slepc->SetTol(eig_tol);
    slepc->SetWhichEigenpairs(which_eig);
    slepc->SetLinearSolver(*ksp);
    eigen = std::move(slepc);
#endif
  }
}

}  // namespace palace
