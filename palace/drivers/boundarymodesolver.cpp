// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include "linalg/errorestimator.hpp"
#include "linalg/modeeigensolver.hpp"
#include "linalg/operator.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
BoundaryModeSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &bm_data = iodata.solver.boundary_mode;
  const int num_modes = bm_data.n;
  const double freq_GHz = bm_data.freq;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  // Construct the boundary mode operator (mesh, spaces, materials, boundary ops).
  BlockTimer bt0(Timer::CONSTRUCT);
  BoundaryModeOperator mode_op(iodata, mesh);

  // Construct the eigenvalue solver on top of the operator's FE context.
  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (int i = 0; i < mode_op.GetH1DbcTDofLists().back().Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + mode_op.GetH1DbcTDofLists().back()[i]);
  }

  const bool have_mg = (mode_op.GetNDSpaceHierarchy().GetNumLevels() > 1);
  if (have_mg)
  {
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(mode_op.GetNDSpaceHierarchy().GetNumLevels()));
  }
  const auto which_eig = (bm_data.target > 0.0)
                             ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                             : EigenvalueSolver::WhichType::LARGEST_REAL;
  ModeEigenSolver eig(mode_op.GetMaterialOp(), mode_op.GetSubmeshNormal(),
                      mode_op.GetSurfZOp(), mode_op.GetFarfieldOp(),
                      mode_op.GetSurfSigmaOp(), mode_op.GetNDSpace(), mode_op.GetH1Space(),
                      dbc_tdof_list, num_modes, bm_data.max_size, bm_data.tol, which_eig,
                      iodata.solver.linear, bm_data.type, iodata.problem.verbose,
                      mode_op.GetComm(), have_mg ? &mode_op : nullptr);

  // Construct the error estimator. BoundaryModeOperator owns the RT hierarchy for
  // flux recovery (analogous to SpaceOperator::rt_fespaces).
  BoundaryModeFluxErrorEstimator<ComplexVector> estimator(
      mode_op.GetMaterialOp(), mode_op.GetNDSpaceHierarchy(),
      mode_op.GetRTSpaceHierarchy(), mode_op.GetCurlSpace(), mode_op.GetH1SpaceHierarchy(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);

  // Construct PostOperator.
  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);
  if (mode_op.IsFromSubmesh())
  {
    post_op.ProjectImpedancePaths(mode_op.GetSubmeshCentroid(), mode_op.GetSubmeshE1(),
                                  mode_op.GetSubmeshE2());
    post_op.ProjectVoltagePaths(mode_op.GetSubmeshCentroid(), mode_op.GetSubmeshE1(),
                                mode_op.GetSubmeshE2());
  }

  ErrorIndicator indicator;

  // Determine kn_target.
  double kn_target;
  if (bm_data.target > 0.0)
  {
    kn_target = bm_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", bm_data.target, kn_target);
  }
  else
  {
    double c_min = mode_op.GetMaterialOp().GetLightSpeedMax().Min();
    Mpi::GlobalMin(1, &c_min, mode_op.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * std::sqrt(1.1);
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
  }

  // Solve the eigenvalue problem.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);
  const double sigma = -kn_target * kn_target;
  auto result = eig.Solve(omega, sigma);
  int num_conv = result.num_converged;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", result.sigma);

  for (int i = 0; i < num_conv; i++)
  {
    auto kn = eig.GetPropagationConstant(i);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i, kn.real(),
               kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = eig.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  const int h1_size = mode_op.GetH1TrueVSize();
  const int l2_size = mode_op.GetCurlSpace().GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    ComplexVector et, en;
    std::complex<double> kn = eig.GetPhysicalMode(i, omega, e0, et, en);
    double error_bkwd = eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eig.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    // Poynting power P for impedance postprocessing (|P| ≈ 1 after normalization).
    std::complex<double> P = eig.ComputePoyntingPower(omega, kn, et, en);

    auto total_domain_energy = post_op.MeasureAndPrintAll(i, et, en, kn, P, omega,
                                                          error_abs, error_bkwd, n_print);

    if (i < num_modes && ModeEigenSolver::IsPropagating(kn))
    {
      // Bz = curl(Et) / (i·omega) = (Im(curl Et) - i·Re(curl Et)) / omega.
      ComplexVector bz(l2_size);
      bz.UseDevice(true);
      const auto &CurlOp =
          mode_op.GetCurlSpace().GetDiscreteInterpolator(mode_op.GetNDSpace());
      Vector curl_etr(l2_size), curl_eti(l2_size);
      curl_etr.UseDevice(true);
      curl_eti.UseDevice(true);
      CurlOp.Mult(et.Real(), curl_etr);
      CurlOp.Mult(et.Imag(), curl_eti);
      bz.Real() = curl_eti;
      bz.Real() *= 1.0 / omega;
      bz.Imag() = curl_etr;
      bz.Imag() *= -1.0 / omega;
      estimator.AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);

  return {indicator,
          mode_op.GetNDSpace().GlobalTrueVSize() + mode_op.GetH1Space().GlobalTrueVSize()};
}

}  // namespace palace
