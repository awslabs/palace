// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeanalysissolver.hpp"

#include <complex>

#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodesolver.hpp"
#include "models/materialoperator.hpp"
#include "models/modeanalysisoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

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

  // Construct RT FE space for error estimation (smooth D-field recovery) and wrap in a
  // single-level hierarchy. For the 2D case, the curl flux estimator also needs the L2
  // curl space (from ModeAnalysisOperator) and an H1 hierarchy.
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
  double c_min = mat_op.GetLightSpeedMax().Min();
  Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light in ModeAnalysisSolver!");
  double sigma;
  if (ma_data.target > 0.0)
  {
    // User-specified target effective index: search for modes near n_eff = target.
    // kn_target = n_eff * omega (nondimensional), sigma = -kn_target^2.
    const double kn_target = ma_data.target * omega;
    sigma = -kn_target * kn_target;
    Mpi::Print(" Target n_eff = {:.6e}, sigma = {:.6e}\n", ma_data.target, sigma);
  }
  else
  {
    // Auto: target the light line from material properties.
    const double mu_eps_max = 1.0 / (c_min * c_min) * 1.1;
    sigma = -omega * omega * mu_eps_max;
  }

  // Configure and construct the boundary mode solver.
  BoundaryModeSolverConfig config;
  config.attr_to_material = &mat_op.GetAttributeToMaterial();
  config.inv_permeability = &mat_op.GetInvPermeability();
  config.curlcurl_inv_permeability = &mat_op.GetCurlCurlInvPermeability();
  config.permittivity_real = &mat_op.GetPermittivityReal();
  config.permittivity_scalar = &mat_op.GetPermittivityScalar();
  config.permittivity_imag =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImag() : nullptr;
  config.has_loss_tangent = mat_op.HasLossTangent();
  config.normal = nullptr;  // 2D domain mesh, no normal projection
  config.num_modes = num_modes;
  config.num_vec = std::max(2 * num_modes + 1, 20);
  config.eig_tol = tol;
  config.which_eig = EigenvalueSolver::WhichType::LARGEST_MAGNITUDE;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = ma_data.type;
  config.verbose = iodata.problem.verbose;

  // Build combined dbc_tdof_list for the block system.
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  BoundaryModeSolver mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list,
                                 mesh.back()->GetComm());

  // Store Btt in ModeAnalysisOperator for impedance postprocessing in PostOperator.
  mode_op.SetBttMatrix(std::make_unique<mfem::HypreParMatrix>(*mode_solver.GetBtt()));

  // Solve the eigenvalue problem.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving for {:d} propagation mode(s)...\n", num_modes);

  auto result = mode_solver.Solve(omega, sigma);
  int num_conv = result.num_converged;
  {
    std::complex<double> lambda = (num_conv > 0) ? mode_solver.GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  // Postprocessing: extract propagation constants, effective indices, and impedance.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = mode_solver.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }

  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  // Get the discrete curl operator for computing Bz = curl_t(Et) / (iω) on the L2 space,
  // needed for the curl flux error estimator.
  auto &l2_curl_fespace = mode_op.GetCurlSpace();
  const auto &CurlOp = l2_curl_fespace.GetDiscreteInterpolator(nd_fespace);
  const int l2_size = l2_curl_fespace.GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    std::complex<double> lambda = mode_solver.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);

    // Extract the tangential (ND) and normal (H1) E eigenvector components.
    ComplexVector et(nd_size), en(h1_size);
    {
      ComplexVector e0(nd_size + h1_size);
      mode_solver.GetEigenvector(i, e0);
      std::copy_n(e0.Real().begin(), nd_size, et.Real().begin());
      std::copy_n(e0.Imag().begin(), nd_size, et.Imag().begin());
      std::copy_n(e0.Real().begin() + nd_size, h1_size, en.Real().begin());
      std::copy_n(e0.Imag().begin() + nd_size, h1_size, en.Imag().begin());
    }

    // Power-normalize eigenvector to unit transmitted Poynting power:
    //   P = (1/2) Re{kn*/omega} * (et^H Btt et)
    // SLEPc/ARPACK normalize to some operator norm, giving arbitrary field magnitudes.
    // Normalizing to P = 1 (nondimensional) makes field magnitudes physically meaningful
    // and comparable to COMSOL (which normalizes to 1 W/m transmitted power).
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
          // Evanescent mode (kn purely imaginary): P = 0 since Re{kn*/omega} = 0.
          // Normalize to unit et^H Btt et instead so field magnitudes are well-defined.
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

    // PostOperator handles all measurements, printing, and CSV output.
    auto total_domain_energy = post_op.MeasureAndPrintAll(i, et, en, kn, omega, n_print);

    // Calculate and record the error indicators. Compute Bz = curl_t(Et) / (iω) on the
    // L2 curl space for the curl flux part of the estimator.
    if (i < num_modes)
    {
      ComplexVector bz(l2_size);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        // Bz = curl(Et) / (iω) = curl(Et) * (-i/ω)
        // Bz_r = curl(Et_i) / ω,  Bz_i = -curl(Et_r) / ω
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
