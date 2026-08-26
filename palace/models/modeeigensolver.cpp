// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeeigensolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/cudss.hpp"
#include "linalg/gmg.hpp"
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
#include "models/surfacerationalimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

extern "C"
{
  void zggev_(char *, char *, int *, std::complex<double> *, int *, std::complex<double> *,
              int *, std::complex<double> *, std::complex<double> *, std::complex<double> *,
              int *, std::complex<double> *, int *, std::complex<double> *, int *, double *,
              int *);
}

namespace palace
{

namespace mode_assembly
{

namespace
{
constexpr bool skip_zeros = false;
}  // namespace

ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op)
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op)
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix AssembleAtt(
    const FiniteElementSpace &nd_fespace, const MaterialOperator &mat_op,
    const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
    FarfieldBoundaryOperator &farfield_op, SurfaceConductivityOperator &surf_sigma_op,
    SurfaceRationalImpedanceOperator &surf_rz_op, std::complex<double> omega, double sigma)
{
  // Complex-ω decomposition. The transverse block carries the second-order operator
  // K + iωC - ω²M (domain + boundary). For real ω the imaginary parts below are all
  // zero and this reduces bit-for-bit to the original real-ω assembly. For complex ω,
  // ω² = (w2r + i·w2i) so the -ω²·ε_real domain mass term acquires an imaginary part
  // -w2i·ε_real (Atti), and the surface RLC / farfield multipliers (iω, -ω²) pick up
  // real-part contributions that move to the real block. Cross terms between Im(ω) and
  // material loss (ε_imag, σ) are included for completeness but are second order.
  const double wr = omega.real(), wi = omega.imag();
  const std::complex<double> w2 = omega * omega;
  const double w2r = w2.real(), w2i = w2.imag();
  const bool complex_omega = (wi != 0.0);

  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient eps_shifted_func(mat_op.GetAttributeToMaterial(),
                                               mat_op.GetPermittivityReal(), -w2r);
  eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvPermeability(), -sigma);
  if (complex_omega && mat_op.HasLossTangent())
  {
    // Cross term: Re(-ω²·(i·ε_imag)) = +w2i·ε_imag.
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetPermittivityImag(), w2i);
  }
  if (complex_omega && mat_op.HasConductivity())
  {
    // Cross term: Re(+iω·σ) = -wi·σ.
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetConductivity(), -wi);
  }
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepth(), 1.0);
  }

  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);
  surf_z_op.AddStiffnessBdrCoefficients(1.0, fbr);  // K: 1/Ls (real, ω-independent)
  surf_z_op.AddDampingBdrCoefficients(wr, fbi);     // C: Im(iω)/Rs = wr
  surf_z_op.AddMassBdrCoefficients(-w2r, fbr);      // M: Re(-ω²)·Cs
  farfield_op.AddDampingBdrCoefficients(wr, fbi);
  if (complex_omega)
  {
    surf_z_op.AddDampingBdrCoefficients(-wi, fbr);  // C: Re(iω)/Rs = -wi
    surf_z_op.AddMassBdrCoefficients(-w2i, fbi);    // M: Im(-ω²)·Cs
    farfield_op.AddDampingBdrCoefficients(-wi, fbr);
  }
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_func);
  if (!fbr.empty())
  {
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }
  auto Attr_assembled =
      ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  std::unique_ptr<mfem::HypreParMatrix> Atti_assembled;
  {
    const bool has_imag = mat_op.HasLossTangent() || mat_op.HasConductivity() ||
                          !fbi.empty() || complex_omega;
    if (has_imag)
    {
      // Coefficients must outlive the BilinearForm (integrators hold raw pointers).
      const int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient negepstandelta_func(n_attr);
      MaterialPropertyCoefficient fi_domain(n_attr);
      if (complex_omega)
      {
        // New term: Im(-ω²·ε_real) = -w2i·ε_real.
        negepstandelta_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(), -w2i);
      }
      if (mat_op.HasLossTangent())
      {
        // Im(-ω²·(i·ε_imag)) = -w2r·ε_imag.
        negepstandelta_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityImag(), -w2r);
      }
      if (mat_op.HasConductivity())
      {
        // Im(+iω·σ) = +wr·σ.
        fi_domain.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(),
                                 wr);
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

ComplexHypreParMatrix
AssembleAnn(const FiniteElementSpace &h1_fespace, const MaterialOperator &mat_op,
            const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
            FarfieldBoundaryOperator &farfield_op,
            SurfaceConductivityOperator &surf_sigma_op,
            SurfaceRationalImpedanceOperator &surf_rz_op, std::complex<double> omega)
{
  // Complex-ω decomposition for the normal (H1) block. Sign convention is OPPOSITE the
  // transverse block: the ε mass term enters as +ω²·ε (vs -ω²·ε in Att), surf-Z damping
  // as -iω, surf-Z mass as +ω², farfield as -iω. For real ω the imaginary parts are
  // zero and this reduces bit-for-bit to the original assembly.
  const double wr = omega.real(), wi = omega.imag();
  const std::complex<double> w2 = omega * omega;
  const double w2r = w2.real(), w2i = w2.imag();
  const bool complex_omega = (wi != 0.0);

  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient poseps_h1_func(
      mat_op.GetAttributeToMaterial(),
      normal ? mat_op.GetPermittivityReal() : mat_op.GetPermittivityScalar(), w2r);
  if (complex_omega && mat_op.HasLossTangent())
  {
    // Cross term: Re(+ω²·(i·ε_imag)) = -w2i·ε_imag — mirror of the +w2i·ε_imag term in
    // AssembleAtt. Accumulate before the (destructive) normal projection below;
    // projection is linear so this equals projecting the sum.
    poseps_h1_func.AddCoefficient(
        mat_op.GetAttributeToMaterial(),
        normal ? mat_op.GetPermittivityImag() : mat_op.GetPermittivityImagScalar(), -w2i);
  }
  if (complex_omega && mat_op.HasConductivity())
  {
    // Cross term: Re(-iω·σ) = +wi·σ — mirror of the -wi·σ term in AssembleAtt (the
    // conduction current enters the normal block with the opposite sign, -iωσ vs +iωσ,
    // following the +ω²ε_c vs -ω²ε_c convention).
    poseps_h1_func.AddCoefficient(
        mat_op.GetAttributeToMaterial(),
        normal ? mat_op.GetConductivity() : mat_op.GetConductivityScalar(), wi);
  }
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
    else
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

  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr), nn_fbi(max_bdr_attr);
  surf_z_op.AddStiffnessBdrCoefficients(-1.0, nn_fbr);
  surf_z_op.AddDampingBdrCoefficients(-wr, nn_fbi);  // Im(-iω)/Rs = -wr
  surf_z_op.AddMassBdrCoefficients(w2r, nn_fbr);
  if (complex_omega)
  {
    surf_z_op.AddDampingBdrCoefficients(wi, nn_fbr);  // Re(-iω)/Rs = wi
    surf_z_op.AddMassBdrCoefficients(w2i, nn_fbi);
  }
  if (farfield_op.GetAttrList().Size() > 0)
  {
    // Farfield boundary: scalar inverse impedance for the H1 mass integrator. The
    // multiplier is -iω·(1/Z0): Im → nn_fbi (·-wr), Re → nn_fbr (·wi).
    const auto &farfield_attrs = farfield_op.GetAttrList();
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
        nn_fbi.AddMaterialProperty(ceed_attrs, inv_z0_scalar, -wr);
        if (complex_omega)
        {
          nn_fbr.AddMaterialProperty(ceed_attrs, inv_z0_scalar, wi);
        }
      }
    }
  }
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr), cond_i(max_bdr_attr);
    surf_sigma_op.AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
    surf_rz_op.AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
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

  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
  if (!nn_fbr.empty())
  {
    annr.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
  }
  auto Annr_assembled =
      ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  std::unique_ptr<mfem::HypreParMatrix> Anni_assembled;
  {
    const bool has_imag = mat_op.HasLossTangent() || mat_op.HasConductivity() ||
                          !nn_fbi.empty() || complex_omega;
    if (has_imag)
    {
      const int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient posepsi_h1_func(n_attr);
      // Accumulate both contributions into the (possibly tensor) coefficient first, then
      // normal-project ONCE at the end: NormalProjectedCoefficient is destructive (it
      // collapses the tensor to a scalar), and projection is linear so the sum of
      // projections equals the projection of the sum.
      if (complex_omega)
      {
        // New term: Im(+ω²·ε_real) = +w2i·ε_real.
        posepsi_h1_func.AddCoefficient(
            mat_op.GetAttributeToMaterial(),
            normal ? mat_op.GetPermittivityReal() : mat_op.GetPermittivityScalar(), w2i);
      }
      if (mat_op.HasLossTangent())
      {
        // Im(+ω²·(i·ε_imag)) = +w2r·ε_imag.
        posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                       normal ? mat_op.GetPermittivityImag()
                                              : mat_op.GetPermittivityImagScalar(),
                                       w2r);
      }
      if (mat_op.HasConductivity())
      {
        // Im(-iω·σ) = -wr·σ: the conduction-current part of the complex permittivity
        // ε_c = ε - iσ/ω in the normal block's +ω²·ε_c term (mirror of the +wr·σ term in
        // AssembleAtt, opposite sign convention).
        posepsi_h1_func.AddCoefficient(
            mat_op.GetAttributeToMaterial(),
            normal ? mat_op.GetConductivity() : mat_op.GetConductivityScalar(), -wr);
      }
      if (normal && !posepsi_h1_func.empty())
      {
        posepsi_h1_func.NormalProjectedCoefficient(*normal);
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

void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en)
{
  et.Real().MakeRef(e0.Real(), 0, nd_size);
  et.Imag().MakeRef(e0.Imag(), 0, nd_size);
  en.Real().MakeRef(e0.Real(), nd_size, h1_size);
  en.Imag().MakeRef(e0.Imag(), nd_size, h1_size);
  const auto ikn_inv = 1.0 / (std::complex<double>(0.0, 1.0) * kn);
  ComplexVector::AXPBY(ikn_inv, en.Real(), en.Imag(), 0.0, en.Real(), en.Imag());
}

}  // namespace mode_assembly

using namespace std::complex_literals;

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

ModeEigenSolver::ModeEigenSolver(
    const MaterialOperator &mat_op, const mfem::Vector *normal,
    SurfaceImpedanceOperator &surf_z_op, FarfieldBoundaryOperator &farfield_op,
    SurfaceConductivityOperator &surf_sigma_op,
    SurfaceRationalImpedanceOperator &surf_rz_op, const FiniteElementSpace &nd_fespace,
    const FiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_tdof_list,
    int num_modes, int num_vec, double eig_tol, EigenvalueSolver::WhichType which_eig,
    const config::LinearSolverData &linear, EigenSolverBackend eigen_backend, int verbose,
    MPI_Comm solver_comm)
  : num_modes(num_modes), num_vec(num_vec), eig_tol(eig_tol), which_eig(which_eig),
    linear(linear), eigen_backend(eigen_backend), verbose(verbose), mat_op(mat_op),
    normal(normal), surf_z_op(surf_z_op), farfield_op(farfield_op),
    surf_sigma_op(surf_sigma_op), surf_rz_op(surf_rz_op), nd_fespace(nd_fespace),
    h1_fespace(h1_fespace), dbc_tdof_list(dbc_tdof_list)
{
  // Assemble Atn, Btn = -Atn^T, Btt locally (no BMO available on this path).
  std::tie(owned_Atnr, owned_Atni) =
      mode_assembly::AssembleAtn(nd_fespace, h1_fespace, mat_op);
  owned_Btnr.reset(owned_Atnr->Transpose());
  *owned_Btnr *= -1.0;
  if (owned_Atni)
  {
    owned_Btni.reset(owned_Atni->Transpose());
    *owned_Btni *= -1.0;
  }
  auto [bttr_tmp, btti_tmp] = mode_assembly::AssembleBtt(nd_fespace, mat_op);
  owned_Bttr = std::move(bttr_tmp);

  Atnr = owned_Atnr.get();
  Atni = owned_Atni.get();
  Btnr = owned_Btnr.get();
  Btni = owned_Btni.get();
  Bttr = owned_Bttr.get();

  Init(solver_comm);
}

ModeEigenSolver::ModeEigenSolver(BoundaryModeOperator &bmo,
                                 const mfem::Array<int> &dbc_tdof_list, int num_modes,
                                 int num_vec, double eig_tol,
                                 EigenvalueSolver::WhichType which_eig,
                                 const config::LinearSolverData &linear,
                                 EigenSolverBackend eigen_backend, int verbose)
  : num_modes(num_modes), num_vec(num_vec), eig_tol(eig_tol), which_eig(which_eig),
    linear(linear), eigen_backend(eigen_backend), verbose(verbose),
    mat_op(bmo.GetMaterialOp()), normal(nullptr), surf_z_op(bmo.GetSurfZOp()),
    farfield_op(bmo.GetFarfieldOp()), surf_sigma_op(bmo.GetSurfSigmaOp()),
    surf_rz_op(bmo.GetSurfRZOp()), nd_fespace(bmo.GetNDSpace()),
    h1_fespace(bmo.GetH1Space()), bmo(&bmo), dbc_tdof_list(dbc_tdof_list)
{
  // Alias BMO-owned frequency-independent matrices; no local assembly.
  Atnr = bmo.GetAtnr();
  Atni = bmo.GetAtni();
  Btnr = bmo.GetBtnr();
  Btni = bmo.GetBtni();
  Bttr = bmo.GetBtt();

  Init(nd_fespace.GetComm());
}

void ModeEigenSolver::Init(MPI_Comm solver_comm)
{
  nd_size = nd_fespace.GetTrueVSize();
  h1_size = h1_fespace.GetTrueVSize();

  // Build the frequency-independent block B matrix. Btn is real-only; Btt and the zero
  // H1 diagonal placeholder make up the rest.
  Vector d(h1_size);
  d.UseDevice(false);
  d = 0.0;
  mfem::SparseMatrix diag(d);
  auto Dnn = std::make_unique<mfem::HypreParMatrix>(
      h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
      h1_fespace.Get().GetTrueDofOffsets(), &diag);
  auto [Br, Bi] = BuildSystemMatrixB(Bttr, nullptr, Btnr, Btni, Dnn.get());
  opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));

  // Pick a communicator for the solvers: `solver_comm` if non-null (WavePort sub-
  // communicator), else the FE space comm. Ranks with no DOFs and no solver_comm
  // (WavePort non-port ranks) skip solver setup entirely.
  const bool use_mg = bmo && bmo->GetNDSpaceHierarchy().GetNumLevels() > 1;
  MPI_Comm configure_comm = (solver_comm != MPI_COMM_NULL) ? solver_comm
                            : (nd_size > 0)                ? nd_fespace.GetComm()
                                                           : MPI_COMM_NULL;
  this->solver_comm = configure_comm;
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

void ModeEigenSolver::AssembleFrequencyDependent(std::complex<double> omega, double sigma)
{
  reduced_stats.full_operator_assemblies++;
  last_assembled_omega = omega;
  last_assembled_sigma = sigma;

  // Frequency-dependent Att/Ann: delegate to BMO on the 2D domain path; otherwise
  // assemble locally via the shared free functions.
  std::unique_ptr<mfem::HypreParMatrix> Attr, Atti, Annr_local, Anni_local;
  if (bmo)
  {
    std::tie(Attr, Atti) = bmo->AssembleAtt(omega, sigma);
    std::tie(Annr_local, Anni_local) = bmo->AssembleAnn(omega);
  }
  else
  {
    std::tie(Attr, Atti) =
        mode_assembly::AssembleAtt(nd_fespace, mat_op, normal, surf_z_op, farfield_op,
                                   surf_sigma_op, surf_rz_op, omega, sigma);
    std::tie(Annr_local, Anni_local) =
        mode_assembly::AssembleAnn(h1_fespace, mat_op, normal, surf_z_op, farfield_op,
                                   surf_sigma_op, surf_rz_op, omega);
  }

  // Shifted (1,0) block: -sigma * Btn_r (real-only).
  std::unique_ptr<mfem::HypreParMatrix> shifted_Btnr;
  if (Btnr && std::abs(sigma) > 0.0)
  {
    shifted_Btnr = std::make_unique<mfem::HypreParMatrix>(*Btnr);
    *shifted_Btnr *= -sigma;
  }

  auto [Ar, Ai] = BuildSystemMatrixA(Attr.get(), Atti.get(), Atnr, Atni, Annr_local.get(),
                                     Anni_local.get(), shifted_Btnr.get());
  opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

ModeEigenSolver::SolveResult ModeEigenSolver::Solve(std::complex<double> omega,
                                                    double sigma,
                                                    const ComplexVector *initial_space)
{
  sigma_cached = sigma;
  reduced_solution = false;

  const bool real_frequency = (omega.imag() == 0.0);
  const bool can_reduce = real_frequency && reduced_enabled && affine_model_ready &&
                          reduced_basis.size() >= static_cast<std::size_t>(num_modes);
  if (can_reduce && reduced_solves_since_exact < REDUCED_EXACT_CHECK_INTERVAL)
  {
    // The accepted online path evaluates only preprojected exact affine components. It does
    // not assemble the full finite-element operator.
    if (TryAffineReducedSolve(omega.real(), sigma))
    {
      reduced_solution = true;
      reduced_stats.reduced_solves++;
      reduced_stats.affine_reduced_solves++;
      reduced_solves_since_exact++;
      return {static_cast<int>(reduced_eigenvalues.size()), sigma};
    }
    reduced_stats.reduced_fallbacks++;
  }
  else if (can_reduce)
  {
    // Periodically refresh the truth subspace even when every reduced residual passes.
    // A residual can certify a represented Ritz pair but cannot prove that a new mode has
    // not entered the requested rank from outside the current reduced basis.
    reduced_stats.periodic_exact_checks++;
  }

  // Truth/fallback and all complex-frequency queries retain the original exact assembly.
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
    // The preconditioner is a real approximation built at Re(omega); the imaginary parts
    // of the complex-ω operator are carried only by the monolithic opA above, so the
    // matching system matrix still sees the full complex frequency.
    att_mg_op = AssembleAttPreconditioner(omega.real(), sigma);
    ann_mg_op = AssembleAnnPreconditioner(omega.real());
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

  if (real_frequency && warm_start.Size() > 0)
  {
    eigen->SetInitialSpace(warm_start);
  }
  else if (initial_space)
  {
    eigen->SetInitialSpace(*initial_space);
  }

  int num_conv = eigen->Solve();
  if (real_frequency)
  {
    reduced_stats.exact_solves++;
    reduced_solves_since_exact = 0;
  }
  else
  {
    reduced_stats.complex_exact_solves++;
  }

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

  if (real_frequency)
  {
    EnrichReducedBasis(num_conv);
    if (affine_model_ready)
    {
      const double discrepancy = ValidateAffineModel(omega.real(), sigma);
      reduced_stats.worst_affine_discrepancy =
          std::max(reduced_stats.worst_affine_discrepancy, discrepancy);
      if (!(discrepancy <= 1.0e-9))
      {
        affine_model_ready = false;
        Mpi::Warning(solver_comm,
                     "Disabling affine wave-port reduced operator after exact refresh: "
                     "component/action discrepancy {:.3e} exceeds 1e-9.\n",
                     discrepancy);
      }
    }
  }
  return {num_conv, sigma};
}

void ModeEigenSolver::SetReducedModelTraining(bool enable, std::size_t max_basis_size)
{
  reduced_training = enable;
  reduced_basis_cap = std::max(max_basis_size, static_cast<std::size_t>(num_modes));
}

void ModeEigenSolver::EnableReducedModel(double adaptive_tol)
{
  reduced_training = true;
  reduced_enabled = !reduced_basis.empty();
  // Do not ask the reduced model to be more accurate than the truth eigensolves used to
  // construct it. For loose adaptive sweeps, cap the port backward error at 1e-4.
  reduced_tol = std::max(10.0 * eig_tol, std::min(1.0e-4, 0.01 * adaptive_tol));
  if (reduced_enabled && !bmo)
  {
    BuildAffineModel();
  }
}

void ModeEigenSolver::BuildAffineModel()
{
  affine_model_ready = false;
  affine_components.clear();
  affine_BV.clear();
  affine_Br.resize(0, 0);

  // The explicit component construction is currently for the bare wave-port path. The
  // BoundaryMode path owns a separate hierarchy/operator assembly and never enables the
  // adaptive wave-port reduced model.
  if (bmo || reduced_basis.empty())
  {
    return;
  }

  auto assemble_nd = [&](BilinearForm &form)
  {
    return ParOperator(form.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
  };
  auto assemble_h1 = [&](BilinearForm &form)
  {
    return ParOperator(form.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
  };
  auto add_component =
      [&](AffineCoefficientType type, int index, std::unique_ptr<mfem::HypreParMatrix> Attr,
          std::unique_ptr<mfem::HypreParMatrix> Atti, const mfem::HypreParMatrix *Atn_real,
          const mfem::HypreParMatrix *Atn_imag, std::unique_ptr<mfem::HypreParMatrix> Annr,
          std::unique_ptr<mfem::HypreParMatrix> Anni,
          std::unique_ptr<mfem::HypreParMatrix> lower_left,
          Operator::DiagonalPolicy diag_policy)
  {
    auto [Ar, Ai] =
        BuildSystemMatrixA(Attr.get(), Atti.get(), Atn_real, Atn_imag, Annr.get(),
                           Anni.get(), lower_left.get(), diag_policy);
    if (Ar || Ai)
    {
      AffineComponent component;
      component.type = type;
      component.index = index;
      component.op = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
      affine_components.push_back(std::move(component));
    }
  };

  const auto &attr_to_mat = mat_op.GetAttributeToMaterial();
  const int n_attr = attr_to_mat.Size();
  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();

  // Constant block: curl-curl / diffusion, London mass, surface inductance, and Atn.
  {
    MaterialPropertyCoefficient muinv_cc(attr_to_mat,
                                         normal ? mat_op.GetInvPermeability()
                                                : mat_op.GetCurlCurlInvPermeability());
    if (normal)
    {
      muinv_cc.NormalProjectedCoefficient(*normal);
    }
    BilinearForm att(nd_fespace);
    att.AddDomainIntegrator<CurlCurlIntegrator>(muinv_cc);
    MaterialPropertyCoefficient london_t(n_attr);
    if (mat_op.HasLondonDepth())
    {
      london_t.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepth());
      att.AddDomainIntegrator<VectorFEMassIntegrator>(london_t);
    }
    MaterialPropertyCoefficient surf_k_t(max_bdr_attr);
    surf_z_op.AddStiffnessBdrCoefficients(1.0, surf_k_t);
    if (!surf_k_t.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(surf_k_t);
    }
    auto Attr = assemble_nd(att);

    MaterialPropertyCoefficient neg_muinv(attr_to_mat, mat_op.GetInvPermeability(), -1.0);
    if (normal)
    {
      neg_muinv.NormalProjectedCoefficient(*normal);
    }
    BilinearForm ann(h1_fespace);
    ann.AddDomainIntegrator<DiffusionIntegrator>(neg_muinv);
    MaterialPropertyCoefficient london_n(n_attr);
    if (mat_op.HasLondonDepth())
    {
      if (!normal)
      {
        london_n.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepthScalar());
      }
      else
      {
        const auto &ild = mat_op.GetInvLondonDepth();
        mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
        for (int k = 0; k < ild.SizeK(); k++)
        {
          ild_scalar(0, 0, k) = ild(0, 0, k);
        }
        london_n.AddCoefficient(attr_to_mat, ild_scalar);
        london_n.NormalProjectedCoefficient(*normal);
      }
      ann.AddDomainIntegrator<MassIntegrator>(london_n);
    }
    MaterialPropertyCoefficient surf_k_n(max_bdr_attr);
    surf_z_op.AddStiffnessBdrCoefficients(-1.0, surf_k_n);
    if (!surf_k_n.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(surf_k_n);
    }
    auto Annr = assemble_h1(ann);

    add_component(AffineCoefficientType::CONSTANT, -1, std::move(Attr), nullptr, Atnr, Atni,
                  std::move(Annr), nullptr, nullptr, Operator::DIAG_ONE);
  }

  // Linear-in-omega block: bulk conductivity plus resistive/farfield damping. The factor i
  // is stored in the component's imaginary matrix, so the online scalar is real omega.
  if (mat_op.HasConductivity() || surf_z_op.GetRsAttrList().Size() > 0 ||
      farfield_op.GetAttrList().Size() > 0)
  {
    BilinearForm atti(nd_fespace);
    MaterialPropertyCoefficient cond_t(n_attr);
    if (mat_op.HasConductivity())
    {
      cond_t.AddCoefficient(attr_to_mat, mat_op.GetConductivity());
      atti.AddDomainIntegrator<VectorFEMassIntegrator>(cond_t);
    }
    MaterialPropertyCoefficient damp_t(max_bdr_attr);
    surf_z_op.AddDampingBdrCoefficients(1.0, damp_t);
    farfield_op.AddDampingBdrCoefficients(1.0, damp_t);
    if (!damp_t.empty())
    {
      atti.AddBoundaryIntegrator<VectorFEMassIntegrator>(damp_t);
    }
    auto Atti = assemble_nd(atti);

    BilinearForm anni(h1_fespace);
    MaterialPropertyCoefficient cond_n(n_attr);
    if (mat_op.HasConductivity())
    {
      cond_n.AddCoefficient(
          attr_to_mat, normal ? mat_op.GetConductivity() : mat_op.GetConductivityScalar(),
          -1.0);
      if (normal)
      {
        cond_n.NormalProjectedCoefficient(*normal);
      }
      anni.AddDomainIntegrator<MassIntegrator>(cond_n);
    }
    MaterialPropertyCoefficient damp_n(max_bdr_attr);
    surf_z_op.AddDampingBdrCoefficients(-1.0, damp_n);
    if (farfield_op.GetAttrList().Size() > 0)
    {
      const auto &inv_z = mat_op.GetInvImpedance();
      const auto &bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
      for (auto attr : farfield_op.GetAttrList())
      {
        int mat_idx =
            (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
        const double inv_z0 = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
        auto ceed_attrs = mat_op.GetCeedBdrAttributes(attr);
        if (ceed_attrs.Size() > 0)
        {
          damp_n.AddMaterialProperty(ceed_attrs, inv_z0, -1.0);
        }
      }
    }
    if (!damp_n.empty())
    {
      anni.AddBoundaryIntegrator<MassIntegrator>(damp_n);
    }
    auto Anni = assemble_h1(anni);

    add_component(AffineCoefficientType::OMEGA, -1, nullptr, std::move(Atti), nullptr,
                  nullptr, nullptr, std::move(Anni), nullptr, Operator::DIAG_ZERO);
  }

  // Quadratic block: real permittivity/surface capacitance and dielectric loss. Signs in
  // the transverse and normal blocks follow the original VD formulation.
  {
    BilinearForm attr(nd_fespace), atti(nd_fespace);
    MaterialPropertyCoefficient neg_eps(attr_to_mat, mat_op.GetPermittivityReal(), -1.0);
    attr.AddDomainIntegrator<VectorFEMassIntegrator>(neg_eps);
    MaterialPropertyCoefficient surf_m_t(max_bdr_attr);
    surf_z_op.AddMassBdrCoefficients(-1.0, surf_m_t);
    if (!surf_m_t.empty())
    {
      attr.AddBoundaryIntegrator<VectorFEMassIntegrator>(surf_m_t);
    }
    std::unique_ptr<mfem::HypreParMatrix> Atti;
    if (mat_op.HasLossTangent())
    {
      MaterialPropertyCoefficient neg_eps_i(attr_to_mat, mat_op.GetPermittivityImag(),
                                            -1.0);
      atti.AddDomainIntegrator<VectorFEMassIntegrator>(neg_eps_i);
      Atti = assemble_nd(atti);
    }
    auto Attr = assemble_nd(attr);

    BilinearForm annr(h1_fespace), anni(h1_fespace);
    MaterialPropertyCoefficient eps_n(attr_to_mat, normal ? mat_op.GetPermittivityReal()
                                                          : mat_op.GetPermittivityScalar());
    if (normal)
    {
      eps_n.NormalProjectedCoefficient(*normal);
    }
    annr.AddDomainIntegrator<MassIntegrator>(eps_n);
    MaterialPropertyCoefficient surf_m_n(max_bdr_attr);
    surf_z_op.AddMassBdrCoefficients(1.0, surf_m_n);
    if (!surf_m_n.empty())
    {
      annr.AddBoundaryIntegrator<MassIntegrator>(surf_m_n);
    }
    std::unique_ptr<mfem::HypreParMatrix> Anni;
    if (mat_op.HasLossTangent())
    {
      MaterialPropertyCoefficient eps_i_n(attr_to_mat,
                                          normal ? mat_op.GetPermittivityImag()
                                                 : mat_op.GetPermittivityImagScalar());
      if (normal)
      {
        eps_i_n.NormalProjectedCoefficient(*normal);
      }
      anni.AddDomainIntegrator<MassIntegrator>(eps_i_n);
      Anni = assemble_h1(anni);
    }
    auto Annr = assemble_h1(annr);

    add_component(AffineCoefficientType::OMEGA_SQUARED, -1, std::move(Attr),
                  std::move(Atti), nullptr, nullptr, std::move(Annr), std::move(Anni),
                  nullptr, Operator::DIAG_ZERO);
  }

  // Algebraic shift block: -sigma*Btt in Att and -sigma*Btn in the lower-left block.
  {
    auto Attr = std::make_unique<mfem::HypreParMatrix>(*Bttr);
    *Attr *= -1.0;
    auto lower_left = std::make_unique<mfem::HypreParMatrix>(*Btnr);
    *lower_left *= -1.0;
    add_component(AffineCoefficientType::SHIFT, -1, std::move(Attr), nullptr, nullptr,
                  nullptr, nullptr, nullptr, std::move(lower_left), Operator::DIAG_ZERO);
  }

  auto add_boundary_component =
      [&](AffineCoefficientType type, int index, MaterialPropertyCoefficient &unit)
  {
    BilinearForm att(nd_fespace), ann(h1_fespace);
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(unit);
    auto Attr = assemble_nd(att);
    unit *= -1.0;
    ann.AddBoundaryIntegrator<MassIntegrator>(unit);
    auto Annr = assemble_h1(ann);
    add_component(type, index, std::move(Attr), nullptr, nullptr, nullptr, std::move(Annr),
                  nullptr, nullptr, Operator::DIAG_ZERO);
  };

  // Non-polynomial scalar laws remain exact: one fixed boundary mass per configured group.
  for (std::size_t g = 0; g < surf_sigma_op.Size(); g++)
  {
    if (surf_sigma_op.IsActive(g))
    {
      MaterialPropertyCoefficient unit(max_bdr_attr);
      surf_sigma_op.AddBoundaryMassBdrCoefficients(g, unit);
      add_boundary_component(AffineCoefficientType::SURFACE_CONDUCTIVITY,
                             static_cast<int>(g), unit);
    }
  }
  for (int b = 0; b < surf_rz_op.GetNumBoundaries(); b++)
  {
    MaterialPropertyCoefficient unit(max_bdr_attr);
    surf_rz_op.AddUnitBdrCoefficients(b, unit);
    add_boundary_component(AffineCoefficientType::RATIONAL_IMPEDANCE, b, unit);
  }

  reduced_stats.affine_model_builds++;
  UpdateAffineProjection(0);
  affine_model_ready = solver_comm == MPI_COMM_NULL || !affine_components.empty();

  // Validate the decomposition against the last truth operator assembled during offline
  // sampling. Any mismatch disables the affine path rather than silently omitting physics.
  if (solver_comm != MPI_COMM_NULL && opA && last_assembled_omega.imag() == 0.0)
  {
    const double discrepancy =
        ValidateAffineModel(last_assembled_omega.real(), last_assembled_sigma);
    reduced_stats.worst_affine_discrepancy =
        std::max(reduced_stats.worst_affine_discrepancy, discrepancy);
    if (!(discrepancy <= 1.0e-9))
    {
      affine_model_ready = false;
      Mpi::Warning(solver_comm,
                   "Disabling affine wave-port reduced operator: component/action "
                   "discrepancy {:.3e} exceeds 1e-9.\n",
                   discrepancy);
    }
  }
}

std::complex<double>
ModeEigenSolver::EvaluateAffineCoefficient(const AffineComponent &component, double omega,
                                           double sigma) const
{
  switch (component.type)
  {
    case AffineCoefficientType::CONSTANT:
      return 1.0;
    case AffineCoefficientType::OMEGA:
      return omega;
    case AffineCoefficientType::OMEGA_SQUARED:
      return omega * omega;
    case AffineCoefficientType::SHIFT:
      return sigma;
    case AffineCoefficientType::SURFACE_CONDUCTIVITY:
      return surf_sigma_op.EvaluateScalar(static_cast<std::size_t>(component.index),
                                          std::complex<double>(omega, 0.0));
    case AffineCoefficientType::RATIONAL_IMPEDANCE:
      return surf_rz_op.EvalRobinCoefficient(component.index,
                                             std::complex<double>(0.0, omega));
  }
  MFEM_ABORT("Unknown affine mode-operator coefficient type!");
  return 0.0;
}

void ModeEigenSolver::UpdateAffineProjection(std::size_t old_basis_size)
{
  if (solver_comm == MPI_COMM_NULL)
  {
    return;
  }
  const std::size_t n = reduced_basis.size();
  MFEM_VERIFY(old_basis_size <= n, "Invalid old affine basis size!");

  affine_BV.resize(n);
  for (std::size_t j = old_basis_size; j < n; j++)
  {
    affine_BV[j].SetSize(opB->Height());
    affine_BV[j].UseDevice(true);
    opB->Mult(reduced_basis[j], affine_BV[j]);
  }
  affine_Br.resize(n, n);
  for (std::size_t j = 0; j < n; j++)
  {
    for (std::size_t i = 0; i < n; i++)
    {
      affine_Br(i, j) = linalg::Dot(solver_comm, affine_BV[j], reduced_basis[i]);
    }
  }

  for (auto &component : affine_components)
  {
    component.AV.resize(n);
    for (std::size_t j = old_basis_size; j < n; j++)
    {
      component.AV[j].SetSize(component.op->Height());
      component.AV[j].UseDevice(true);
      component.op->Mult(reduced_basis[j], component.AV[j]);
    }
    component.Ar.resize(n, n);
    for (std::size_t j = 0; j < n; j++)
    {
      for (std::size_t i = 0; i < n; i++)
      {
        component.Ar(i, j) = linalg::Dot(solver_comm, component.AV[j], reduced_basis[i]);
      }
    }
  }
  if (old_basis_size > 0 && n > old_basis_size)
  {
    reduced_stats.affine_projection_extensions++;
  }
}

double ModeEigenSolver::ValidateAffineModel(double omega, double sigma) const
{
  if (solver_comm == MPI_COMM_NULL || !opA || affine_components.empty())
  {
    return 0.0;
  }
  ComplexVector x(opA->Width()), exact(opA->Height()), affine(opA->Height()),
      component_action(opA->Height());
  x.UseDevice(true);
  exact.UseDevice(true);
  affine.UseDevice(true);
  component_action.UseDevice(true);
  double discrepancy = 0.0;
  // Deterministic random probes measure operator-level agreement without the misleading
  // relative amplification obtained from an eigenvector for which A*x is nearly null.
  for (int probe = 0; probe < 2; probe++)
  {
    linalg::SetRandom(solver_comm, x, 9137 + probe);
    linalg::SetSubVector(x, dbc_tdof_list, 0.0);
    linalg::Normalize(solver_comm, x);
    opA->Mult(x, exact);
    affine = 0.0;
    for (const auto &component : affine_components)
    {
      component.op->Mult(x, component_action);
      affine.Add(EvaluateAffineCoefficient(component, omega, sigma), component_action);
    }
    affine.Add(-1.0, exact);
    const double relative = linalg::Norml2(solver_comm, affine) /
                            std::max(linalg::Norml2(solver_comm, exact), 1.0e-300);
    discrepancy = std::max(discrepancy, relative);
  }
  return discrepancy;
}

bool ModeEigenSolver::AddReducedBasisVector(const ComplexVector &x)
{
  if (reduced_basis.size() >= reduced_basis_cap || solver_comm == MPI_COMM_NULL)
  {
    return false;
  }

  ComplexVector q(x);
  const double norm0 = linalg::Norml2(solver_comm, q);
  if (!(norm0 > 0.0))
  {
    return false;
  }
  // Two-pass modified Gram-Schmidt is sufficient for the deliberately small basis and
  // avoids making the port ROM depend on the real-valued 3D PROM orthogonalization code.
  for (int pass = 0; pass < 2; pass++)
  {
    for (const auto &v : reduced_basis)
    {
      const std::complex<double> alpha = linalg::Dot(solver_comm, q, v);
      q.Add(-alpha, v);
    }
  }
  const double norm = linalg::Norml2(solver_comm, q);
  if (!(norm > 1.0e-10 * norm0))
  {
    return false;
  }
  q *= 1.0 / norm;
  reduced_basis.push_back(std::move(q));
  return true;
}

void ModeEigenSolver::EnrichReducedBasis(int num_converged)
{
  if (num_converged <= 0)
  {
    return;
  }

  // The selected highest requested mode is the best portable single-vector warm start for
  // the next exact ARPACK/SLEPc solve. Capture it before WavePortData applies the VD
  // back-transform to its separate output copy.
  const int selected = std::min(num_modes - 1, num_converged - 1);
  warm_start.SetSize(nd_size + h1_size);
  warm_start.UseDevice(true);
  eigen->GetEigenvector(mode_perm[selected], warm_start);
  linalg::Normalize(solver_comm, warm_start);

  if (!reduced_training)
  {
    return;
  }
  const std::size_t old_basis_size = reduced_basis.size();
  const int n = std::min(num_modes, num_converged);
  ComplexVector x(nd_size + h1_size);
  x.UseDevice(true);
  for (int i = 0; i < n && reduced_basis.size() < reduced_basis_cap; i++)
  {
    eigen->GetEigenvector(mode_perm[i], x);
    AddReducedBasisVector(x);
  }
  if (!affine_components.empty() && reduced_basis.size() > old_basis_size)
  {
    UpdateAffineProjection(old_basis_size);
  }
}

bool ModeEigenSolver::TryAssembledReducedSolve(double sigma)
{
  int n = static_cast<int>(reduced_basis.size());
  if (n < num_modes || solver_comm == MPI_COMM_NULL)
  {
    return false;
  }

  Eigen::MatrixXcd Ar(n, n), Br(n, n);
  std::vector<ComplexVector> AV, BV;
  AV.reserve(n);
  BV.reserve(n);
  for (int j = 0; j < n; j++)
  {
    auto &av = AV.emplace_back(opA->Height());
    auto &bv = BV.emplace_back(opB->Height());
    av.UseDevice(true);
    bv.UseDevice(true);
    opA->Mult(reduced_basis[j], av);
    opB->Mult(reduced_basis[j], bv);
    for (int i = 0; i < n; i++)
    {
      // linalg::Dot(comm, x, y) computes y^H x.
      Ar(i, j) = linalg::Dot(solver_comm, av, reduced_basis[i]);
      Br(i, j) = linalg::Dot(solver_comm, bv, reduced_basis[i]);
    }
  }

  return TryReducedSolveFromActions(sigma, Ar, Br, AV, BV);
}

bool ModeEigenSolver::TryAffineReducedSolve(double omega, double sigma)
{
  const int n = static_cast<int>(reduced_basis.size());
  if (!affine_model_ready || n < num_modes || solver_comm == MPI_COMM_NULL)
  {
    return false;
  }

  Eigen::MatrixXcd Ar = Eigen::MatrixXcd::Zero(n, n);
  std::vector<ComplexVector> AV(n);
  for (int j = 0; j < n; j++)
  {
    AV[j].SetSize(nd_size + h1_size);
    AV[j].UseDevice(true);
    AV[j] = 0.0;
  }
  for (const auto &component : affine_components)
  {
    const std::complex<double> coeff = EvaluateAffineCoefficient(component, omega, sigma);
    Ar.noalias() += coeff * component.Ar;
    for (int j = 0; j < n; j++)
    {
      AV[j].Add(coeff, component.AV[j]);
    }
  }
  return TryReducedSolveFromActions(sigma, Ar, affine_Br, AV, affine_BV);
}

bool ModeEigenSolver::TryReducedSolveFromActions(double sigma, const Eigen::MatrixXcd &Ar,
                                                 const Eigen::MatrixXcd &Br,
                                                 const std::vector<ComplexVector> &AV,
                                                 const std::vector<ComplexVector> &BV)
{
  int n = static_cast<int>(reduced_basis.size());
  MFEM_VERIFY(Ar.rows() == n && Ar.cols() == n && Br.rows() == n && Br.cols() == n &&
                  AV.size() == static_cast<std::size_t>(n) &&
                  BV.size() == static_cast<std::size_t>(n),
              "Invalid reduced mode-operator dimensions!");

  // The full eigensolver solves B x = lambda A_sigma x. Preserve that orientation in the
  // dense QZ solve and avoid forming A_sigma^{-1} B because the projected pair can be
  // singular near cutoff.
  Eigen::MatrixXcd Bq = Br, Aq = Ar;
  Eigen::VectorXcd alpha(n), beta(n);
  Eigen::MatrixXcd vl_dummy(1, 1), vr(n, n);
  int ldvl = 1, ldvr = n, lwork = std::max(4 * n, 1), info = 0;
  Eigen::VectorXcd work(lwork);
  Eigen::VectorXd rwork(std::max(8 * n, 1));
  char job_n = 'N', job_v = 'V';
  zggev_(&job_n, &job_v, &n, Bq.data(), &n, Aq.data(), &n, alpha.data(), beta.data(),
         vl_dummy.data(), &ldvl, vr.data(), &ldvr, work.data(), &lwork, rwork.data(),
         &info);
  if (info != 0)
  {
    return false;
  }

  struct Candidate
  {
    std::complex<double> lambda;
    std::complex<double> kn;
    ComplexVector vector;
    double residual;
    double distance;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(n);
  const double kn_target = std::sqrt(-sigma);
  double best_rejected_residual = std::numeric_limits<double>::infinity();
  for (int k = 0; k < n; k++)
  {
    const double scale = std::max(std::abs(alpha(k)), std::abs(beta(k)));
    if (!(scale > 0.0) || std::abs(beta(k)) <= 1.0e-13 * scale ||
        std::abs(alpha(k)) <= 1.0e-13 * scale)
    {
      continue;
    }
    const std::complex<double> lambda = alpha(k) / beta(k);
    const std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    if (!std::isfinite(lambda.real()) || !std::isfinite(lambda.imag()) ||
        !std::isfinite(kn.real()) || !std::isfinite(kn.imag()) ||
        std::abs(kn) <= 1.0e-10 * std::max(kn_target, 1.0))
    {
      continue;
    }

    ComplexVector x(nd_size + h1_size);
    x.UseDevice(true);
    x = 0.0;
    for (int j = 0; j < n; j++)
    {
      x.Add(vr(j, k), reduced_basis[j]);
    }
    const double xnorm = linalg::Norml2(solver_comm, x);
    if (!(xnorm > 0.0))
    {
      continue;
    }
    x *= 1.0 / xnorm;

    ComplexVector ax(nd_size + h1_size), bx(nd_size + h1_size), residual(nd_size + h1_size);
    ax.UseDevice(true);
    bx.UseDevice(true);
    residual.UseDevice(true);
    ax = 0.0;
    bx = 0.0;
    for (int j = 0; j < n; j++)
    {
      ax.Add(vr(j, k), AV[j]);
      bx.Add(vr(j, k), BV[j]);
    }
    residual = bx;
    residual.Add(-lambda, ax);
    const double denom = linalg::Norml2(solver_comm, bx) +
                         std::abs(lambda) * linalg::Norml2(solver_comm, ax);
    const double eta = linalg::Norml2(solver_comm, residual) /
                       std::max(denom, std::numeric_limits<double>::min());
    if (!std::isfinite(eta))
    {
      continue;
    }
    if (eta > reduced_tol)
    {
      // Projection can introduce spurious Ritz roots, especially for the singular,
      // non-normal block pencil. Filter them by the reconstructed full-space residual
      // before applying the physical mode-ordering rule; otherwise a spurious root closer
      // to the shift target can mask an accurate represented mode.
      best_rejected_residual = std::min(best_rejected_residual, eta);
      continue;
    }
    candidates.push_back({lambda, kn, std::move(x), eta, std::abs(kn.real() - kn_target)});
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &a, const Candidate &b) { return a.distance < b.distance; });
  if (candidates.size() < static_cast<std::size_t>(num_modes))
  {
    reduced_stats.last_residual = std::isfinite(best_rejected_residual)
                                      ? best_rejected_residual
                                      : std::numeric_limits<double>::infinity();
    return false;
  }
  // A nearly tied mode just outside the requested set makes rank-based identification
  // ambiguous. Let the exact eigensolver resolve and enrich that local mode cluster.
  if (candidates.size() > static_cast<std::size_t>(num_modes))
  {
    const double gap = candidates[num_modes].distance - candidates[num_modes - 1].distance;
    if (gap <= 1.0e-10 * std::max(kn_target, 1.0))
    {
      return false;
    }
  }

  reduced_eigenvalues.clear();
  reduced_eigenvectors.clear();
  reduced_errors.clear();
  for (int i = 0; i < num_modes; i++)
  {
    reduced_stats.last_residual = candidates[i].residual;
    reduced_stats.worst_accepted_residual =
        std::max(reduced_stats.worst_accepted_residual, candidates[i].residual);
    reduced_eigenvalues.push_back(candidates[i].lambda);
    reduced_eigenvectors.push_back(std::move(candidates[i].vector));
    reduced_errors.push_back(candidates[i].residual);
  }
  return true;
}

std::complex<double> ModeEigenSolver::GetEigenvalue(int i) const
{
  return reduced_solution ? reduced_eigenvalues.at(i) : eigen->GetEigenvalue(mode_perm[i]);
}

double ModeEigenSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  if (reduced_solution)
  {
    return reduced_errors.at(i);
  }
  return eigen->GetError(mode_perm[i], type);
}

void ModeEigenSolver::GetEigenvector(int i, ComplexVector &x) const
{
  if (reduced_solution)
  {
    x = reduced_eigenvectors.at(i);
  }
  else
  {
    eigen->GetEigenvector(mode_perm[i], x);
  }
}

std::complex<double> ModeEigenSolver::GetPropagationConstant(int i) const
{
  return std::sqrt(-sigma_cached - 1.0 / GetEigenvalue(i));
}
bool ModeEigenSolver::IsPropagating(std::complex<double> kn)
{
  return std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
}
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::BuildSystemMatrixA(
    const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
    const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
    const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
    const mfem::HypreParMatrix *shifted_Btnr, Operator::DiagonalPolicy diag_policy) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  // The (1,0) block is -sigma * Btn from the shift-and-invert transformation.
  // Without shift (sigma=0), this block is zero (upper block-triangular).
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);

  // Zero diagonal placeholders let affine components contain only one block while keeping
  // HypreParMatrixFromBlocks informed of both block-row/column dimensions. Their backing
  // storage must outlive every FromBlocks call below.
  std::unique_ptr<mfem::HypreParMatrix> Dtt_zero, Dnn_zero;
  Vector dtt(nd_size), dnn(h1_size);
  dtt.UseDevice(false);
  dnn.UseDevice(false);
  dtt = 0.0;
  dnn = 0.0;
  auto diag_tt = std::make_unique<mfem::SparseMatrix>(dtt);
  auto diag_nn = std::make_unique<mfem::SparseMatrix>(dnn);
  Dtt_zero = std::make_unique<mfem::HypreParMatrix>(
      nd_fespace.Get().GetComm(), nd_fespace.Get().GlobalTrueVSize(),
      nd_fespace.Get().GetTrueDofOffsets(), diag_tt.get());
  Dnn_zero = std::make_unique<mfem::HypreParMatrix>(
      h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
      h1_fespace.Get().GetTrueDofOffsets(), diag_nn.get());

  std::unique_ptr<mfem::HypreParMatrix> Ar;
  if (Attr || Atnr || Annr || shifted_Btnr)
  {
    blocks(0, 0) = Attr ? Attr : Dtt_zero.get();
    blocks(0, 1) = Atnr;
    blocks(1, 0) = shifted_Btnr;  // -sigma * Btn (nullptr when sigma=0)
    blocks(1, 1) = Annr ? Annr : Dnn_zero.get();
    Ar.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti || Atni || Anni)
  {
    // HypreParMatrixFromBlocks requires at least one non-null block per row and column
    // to determine sizes. Since (1,0) is always null (shifted Btn is real-only), add
    // zero diagonal placeholders when an entire block row or column would be null.
    //
    // The 4-arg HypreParMatrix(comm, glob, row_starts, &diag) constructor does NOT
    // deep-copy `diag`: it aliases the SparseMatrix's CSR arrays (CopyCSR with
    // mem_owner=false). The backing Vector + SparseMatrix must therefore outlive both
    // the placeholder HypreParMatrix AND the HypreParMatrixFromBlocks call below (which
    // reads each block's diag via hypre_MergeDiagAndOffd). Keep them in this outer scope
    // — declaring them inside the `if` blocks would free them before FromBlocks runs,
    // leaving the placeholder pointing at freed memory (use-after-free → segfault).
    blocks(0, 0) = Atti ? Atti : Dtt_zero.get();
    blocks(0, 1) = Atni;
    blocks(1, 0) = nullptr;  // Shifted Btn is real-only
    blocks(1, 1) = Anni ? Anni : Dnn_zero.get();
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs constrained by Dirichlet BCs.
  if (Ar)
  {
    Ar->EliminateBC(dbc_tdof_list, diag_policy);
  }
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
#elif defined(MFEM_USE_CUDSS)
    pc_type = LinearSolver::CUDSS;
#else
    MFEM_ABORT("ModeEigenSolver requires building with SuperLU_DIST, STRUMPACK, "
               "MUMPS, or cuDSS!");
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
  else if (pc_type == LinearSolver::CUDSS)
  {
#if !defined(MFEM_USE_CUDSS)
    MFEM_ABORT("Solver was not built with cuDSS support, please choose a "
               "different solver");
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
        else if (pc_type == LinearSolver::CUDSS)
        {
#if defined(MFEM_USE_CUDSS)
          return std::make_unique<CuDSSSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
                                               linear.sym_factorization, true, verbose - 1);
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
  MFEM_VERIFY(bmo, "Multigrid linear solver requires BMO ctor (2D domain path)!");
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
        MFEM_VERIFY(bmo, "AMS coarse solver requires BMO ctor (2D domain path)!");
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
      case LinearSolver::CUDSS:
#if defined(MFEM_USE_CUDSS)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<CuDSSSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
                                          linear.sym_factorization, true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with cuDSS support!");
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
  MFEM_VERIFY(bmo, "AssembleAttPreconditioner requires BMO ctor (2D domain path)!");
  const auto n_levels = bmo->GetNDSpaceHierarchy().GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients (same at all levels — indexed by element attribute, not by p).
  // The 2D MaterialOperator exposes the scalar out-of-plane μ⁻¹ directly via
  // GetCurlCurlInvPermeability(); material rotation (if any) was baked into iodata before
  // the operator was built.
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            mat_op.GetCurlCurlInvPermeability());

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
  surf_z_op.AddStiffnessBdrCoefficients(1.0, fbr);
  surf_z_op.AddMassBdrCoefficients(shifted ? std::abs(omega * omega) : (-omega * omega),
                                   fbr);
  farfield_op.AddDampingBdrCoefficients(omega, fbr);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbr);
  surf_rz_op.AddExtraSystemBdrCoefficients(omega, fbr, fbr);

  // Assemble ND operators at all levels from one BilinearForm. Unlike SpaceOperator's
  // configured-coarse/active-fine split, these domain coefficients are always nonzero.
  constexpr bool assemble_q_data = false;
  {
    BilinearForm att(bmo->GetNDSpaceHierarchy().GetFinestFESpace());
    att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_pc);
    if (!fbr.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
    }
    auto att_ops = att.Assemble(bmo->GetNDSpaceHierarchy(), skip_zeros);

    // Assemble the auxiliary H1 operators for Hiptmair smoothing from the same always-
    // nonzero mass coefficient.
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
  MFEM_VERIFY(bmo, "AssembleAnnPreconditioner requires BMO ctor (2D domain path)!");
  const auto n_levels = bmo->GetH1SpaceHierarchy().GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients matching AssembleAnn (real part only). The negative diffusion
  // sign from the IBP is preserved — this matches the actual operator. In-plane tensors
  // (μ⁻¹, ε) and scalar out-of-plane components are already in the local frame.
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);

  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityScalar(), omega * omega);
  if (mat_op.HasLondonDepth())
  {
    poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvLondonDepthScalar());
  }

  // Boundary coefficients (real part only).
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr);
  surf_z_op.AddStiffnessBdrCoefficients(-1.0, nn_fbr);
  surf_z_op.AddMassBdrCoefficients(omega * omega, nn_fbr);
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr);
    surf_sigma_op.AddExtraSystemBdrCoefficients(omega, cond_r, cond_r);
    surf_rz_op.AddExtraSystemBdrCoefficients(omega, cond_r, cond_r);
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
