// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

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

  // Construct the boundary mode operator (mesh, spaces, materials, eigen solver).
  BlockTimer bt0(Timer::CONSTRUCT);
  BoundaryModeOperator mode_op(iodata, mesh);

  // Construct PostOperator.
  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);
  if (mode_op.IsFromSubmesh())
  {
    post_op.ProjectImpedancePaths(mode_op.GetSubmeshCentroid(), mode_op.GetSubmeshE1(),
                                  mode_op.GetSubmeshE2());
    post_op.ProjectVoltagePaths(mode_op.GetSubmeshCentroid(), mode_op.GetSubmeshE1(),
                                mode_op.GetSubmeshE2());
  }

  // Error indicator (estimator is owned by mode_op).
  auto &nd_fespace = mode_op.GetNDSpace();
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
    Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * std::sqrt(1.1);
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
  }

  // Solve the eigenvalue problem.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);
  auto result = mode_op.Solve(omega, kn_target);
  int num_conv = result.num_converged;
  double sigma = result.sigma;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", sigma);

  for (int i = 0; i < num_conv; i++)
  {
    std::complex<double> lambda = mode_op.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i, kn.real(),
               kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = mode_op.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  const int nd_size = mode_op.GetNDTrueVSize();
  const int h1_size = mode_op.GetH1TrueVSize();
  auto &l2_curl_fespace = mode_op.GetCurlSpace();
  const auto &CurlOp = l2_curl_fespace.GetDiscreteInterpolator(nd_fespace);
  const int l2_size = l2_curl_fespace.GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    std::complex<double> lambda = mode_op.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    double error_bkwd = mode_op.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = mode_op.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    mode_op.GetEigenvector(i, e0);
    ComplexVector et, en;
    et.Real().MakeRef(e0.Real(), 0, nd_size);
    et.Imag().MakeRef(e0.Imag(), 0, nd_size);
    en.Real().MakeRef(e0.Real(), nd_size, h1_size);
    en.Imag().MakeRef(e0.Imag(), nd_size, h1_size);

    // Power-normalize eigenvector using the full Poynting integral:
    //   P = (1/2) Re{conj(kn)/ω × et^H Btt et} + Re{1/(2ωkn) × et^H Atn ẽn}
    // The second term accounts for the Et·∇t(En) cross-coupling. At this point en = ẽn
    // (VD variable, not yet back-transformed to physical En).
    {
      const auto *Btt = mode_op.GetBtt();
      if (Btt)
      {
        auto etH_Btt_et = linalg::Dot(nd_fespace.GetComm(), et, *Btt, et);
        std::complex<double> P = 0.5 * std::conj(kn) / omega * etH_Btt_et;

        // Cross-term: Re{1/(2ωkn) × et^H Atn ẽn}.
        const auto *Atnr = mode_op.GetAtnr();
        if (Atnr)
        {
          ComplexWrapperOperator Atn(Atnr, mode_op.GetAtni());
          P += 1.0 / (2.0 * omega * kn) * linalg::Dot(nd_fespace.GetComm(), en, Atn, et);
        }

        double P_abs = std::abs(P);
        if (P_abs > 0.0)
        {
          e0 *= 1.0 / std::sqrt(P_abs);
        }
        else
        {
          // Fallback to Btt mass norm (etH_Btt_et is real for symmetric Btt).
          double norm2 = etH_Btt_et.real();
          if (norm2 > 0.0)
          {
            e0 *= 1.0 / std::sqrt(norm2);
          }
        }
      }
    }

    // Back-transform en from VD variable ẽn = ikn·En to physical En = ẽn/(ikn).
    {
      auto ikn_inv = 1.0 / (std::complex<double>(0.0, 1.0) * kn);
      ComplexVector::AXPBY(ikn_inv, en.Real(), en.Imag(), 0.0, en.Real(), en.Imag());
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, et, en, kn, omega, error_abs, error_bkwd, n_print);

    const bool is_propagating =
        std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
    if (i < num_modes && is_propagating)
    {
      ComplexVector bz(l2_size);
      bz.UseDevice(true);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        curl_etr.UseDevice(true);
        curl_eti.UseDevice(true);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        bz.Real() = curl_eti;
        bz.Real() *= 1.0 / omega;
        bz.Imag() = curl_etr;
        bz.Imag() *= -1.0 / omega;
      }
      mode_op.AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);

  return {indicator,
          mode_op.GetNDSpace().GlobalTrueVSize() + mode_op.GetH1Space().GlobalTrueVSize()};
}

}  // namespace palace
