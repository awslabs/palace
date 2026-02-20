// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "quadmodeanalysissolver.hpp"

#include <complex>

#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodesolver.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeanalysisoperator.hpp"
#include "models/postoperator.hpp"
#include "models/quadboundarymodesolver.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
QuadModeAnalysisSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  MFEM_VERIFY(mesh.back()->Dimension() == 2,
              "QuadModeAnalysis solver requires a 2D mesh (waveguide cross-section)!");

  const auto &ma_data = iodata.solver.mode_analysis;
  const double freq_GHz = ma_data.freq;
  const int num_modes = ma_data.n;
  const double tol = ma_data.tol;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide QEP mode analysis at f = {:.3e} GHz (omega = "
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

  // Essential (Dirichlet) BCs: PEC boundaries get Dirichlet on both ND (et×n=0) and
  // H1 (en=0). Impedance/absorbing/conductivity boundaries: NO Dirichlet — the impedance
  // Robin BC enters the M0 matrix as boundary mass on both et and en.
  const auto &pmesh = mesh.back()->Get();
  int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  // Both ND and H1 use the same Dirichlet attributes: PEC only.
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

  // Material operator and PostOperator.
  MaterialOperator mat_op(iodata, *mesh.back());
  ModeAnalysisOperator mode_op(mat_op, nd_fespace, h1_fespace, *mesh.back(),
                               iodata.solver.order);
  PostOperator<ProblemType::MODEANALYSIS> post_op(iodata, mode_op);

  // Error estimator setup.
  auto rt_fec = std::make_unique<mfem::RT_FECollection>(iodata.solver.order - 1,
                                                        mesh.back()->Dimension());
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(*mesh.back(), nd_fec.get()));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(*mesh.back(), rt_fec.get()));
  FiniteElementSpaceHierarchy h1_fespaces_est(
      std::make_unique<FiniteElementSpace>(*mesh.back(), h1_fec.get()));
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      mat_op, nd_fespaces, rt_fespaces, iodata.solver.linear.estimator_tol,
      iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg,
      &mode_op.GetCurlSpace(), &h1_fespaces_est);
  ErrorIndicator indicator;

  // Determine kn_target from user-specified n_eff target or material properties.
  double kn_target;
  if (ma_data.target > 0.0)
  {
    kn_target = ma_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", ma_data.target, kn_target);
  }
  else
  {
    double c_min = mat_op.GetLightSpeedMax().Min();
    Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * 1.05;  // Slightly above light line
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
  }

  // Construct boundary operators for impedance, absorbing, and conductivity BCs.
  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, *mesh.back());
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, *mesh.back());
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, *mesh.back());

  // Configure solver.
  BoundaryModeSolverConfig config;
  config.attr_to_material = &mat_op.GetAttributeToMaterial();
  config.inv_permeability = &mat_op.GetInvPermeability();
  config.curlcurl_inv_permeability = &mat_op.GetCurlCurlInvPermeability();
  config.permittivity_real = &mat_op.GetPermittivityReal();
  config.permittivity_scalar = &mat_op.GetPermittivityScalar();
  config.permittivity_imag =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImag() : nullptr;
  config.has_loss_tangent = mat_op.HasLossTangent();
  config.conductivity = mat_op.HasConductivity() ? &mat_op.GetConductivity() : nullptr;
  config.has_conductivity = mat_op.HasConductivity();
  config.mat_op = &mat_op;
  config.surf_z_op = &surf_z_op;
  config.farfield_op = &farfield_op;
  config.surf_sigma_op = &surf_sigma_op;
  config.normal = nullptr;  // 2D domain mesh
  config.num_modes = num_modes;
  config.num_vec = std::max(2 * num_modes + 1, 20);
  config.eig_tol = tol;
  config.which_eig = EigenvalueSolver::WhichType::TARGET_MAGNITUDE;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = ma_data.type;
  config.verbose = iodata.problem.verbose;

  // Build combined dbc_tdof_list for the block system (PEC only on both ND and H1).
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  // Construct the QEP boundary mode solver.
  QuadBoundaryModeSolver mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list, omega,
                                     mesh.back()->GetComm());

  // Store Btt for impedance postprocessing.
  mode_op.SetBttMatrix(std::make_unique<mfem::HypreParMatrix>(*mode_solver.GetBtt()));

  // Solve the QEP.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving QEP for {:d} propagation mode(s)...\n", num_modes);

  auto result = mode_solver.Solve(kn_target);
  int num_conv = result.num_converged;

  // Print ALL converged eigenvalues for debugging.
  for (int i = 0; i < num_conv; i++)
  {
    auto kn_i = mode_solver.GetKn(i);
    auto lam_i = mode_solver.GetEigenvalue(i);
    double neff_r = kn_i.real() / omega;
    double neff_i = kn_i.imag() / omega;
    Mpi::Print(" [QEP] eig {:d}: lambda = {:.6e}{:+.6e}i, kn = {:.6e}{:+.6e}i, "
               "n_eff = {:.6e}{:+.6e}i\n",
               i, lam_i.real(), lam_i.imag(), kn_i.real(), kn_i.imag(), neff_r, neff_i);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = mode_solver.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }

  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  auto &l2_curl_fespace = mode_op.GetCurlSpace();
  const auto &CurlOp = l2_curl_fespace.GetDiscreteInterpolator(nd_fespace);
  const int l2_size = l2_curl_fespace.GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    std::complex<double> kn = mode_solver.GetKn(i);

    // Extract et (ND) and en (H1) from eigenvector.
    ComplexVector et(nd_size), en(h1_size);
    {
      ComplexVector e0(nd_size + h1_size);
      mode_solver.GetEigenvector(i, e0);
      std::copy_n(e0.Real().begin(), nd_size, et.Real().begin());
      std::copy_n(e0.Imag().begin(), nd_size, et.Imag().begin());
      std::copy_n(e0.Real().begin() + nd_size, h1_size, en.Real().begin());
      std::copy_n(e0.Imag().begin() + nd_size, h1_size, en.Imag().begin());
    }

    // Power-normalize eigenvector.
    {
      const auto *Btt = mode_op.GetBtt();
      if (Btt)
      {
        Vector Btt_etr(nd_size), Btt_eti(nd_size);
        Btt->Mult(et.Real(), Btt_etr);
        Btt->Mult(et.Imag(), Btt_eti);
        double p_rr = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_etr);
        double p_ii = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_eti);
        double p_ri = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_eti);
        double p_ir = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_etr);
        std::complex<double> etH_Btt_et(p_rr + p_ii, p_ri - p_ir);
        std::complex<double> P = 0.5 * std::conj(kn) / omega * etH_Btt_et;
        double P_abs = std::abs(P);
        if (P_abs > 0.0)
        {
          double scale = 1.0 / std::sqrt(P_abs);
          et.Real() *= scale;
          et.Imag() *= scale;
          en.Real() *= scale;
          en.Imag() *= scale;
        }
        else
        {
          double norm2 = p_rr + p_ii;
          if (norm2 > 0.0)
          {
            double scale = 1.0 / std::sqrt(norm2);
            et.Real() *= scale;
            et.Imag() *= scale;
            en.Real() *= scale;
            en.Imag() *= scale;
          }
        }
      }
    }

    auto total_domain_energy = post_op.MeasureAndPrintAll(i, et, en, kn, omega, n_print);

    const bool is_propagating =
        std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
    if (i < num_modes && is_propagating)
    {
      ComplexVector bz(l2_size);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        bz.Real() = curl_eti;
        bz.Real() *= 1.0 / omega;
        bz.Imag() = curl_etr;
        bz.Imag() *= -1.0 / omega;
      }
      estimator.AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
