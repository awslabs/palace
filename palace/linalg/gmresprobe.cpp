// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmresprobe.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mfem.hpp>

namespace palace
{

template <typename VecType>
ProbeGmres<VecType>::ProbeGmres(MPI_Comm comm, const VecType &v0, int max_its)
  : comm(comm), V(max_its + 1), H(Eigen::MatrixXcd::Zero(max_its + 1, max_its))
{
  V[0].SetSize(v0.Size());
  V[0].UseDevice(true);
  V[0] = v0;
}

template <typename VecType>
void ProbeGmres<VecType>::Step(const std::function<void(const VecType &, VecType &)> &apply)
{
  MFEM_VERIFY(its + 1 < static_cast<int>(V.size()) && !Stopped(),
              "No more iterations for the GMRES probe run!");
  const int k = its;
  auto &w = V[k + 1];
  w.SetSize(V[k].Size());
  w.UseDevice(true);
  apply(V[k], w);
  for (int pass = 0; pass < 2; pass++)
  {
    for (int i = 0; i <= k; i++)
    {
      const auto h = linalg::Dot(comm, w, V[i]);
      H(i, k) += h;
      linalg::AXPY(-h, V[i], w);
    }
  }
  const double norm = linalg::Norml2(comm, w);
  H(k + 1, k) = norm;
  its = k + 1;

  // Breakdown when M A v vanishes (the residual is unchanged), or when the new direction
  // lies in the Krylov space to round-off relative to ‖M A v‖.
  constexpr double breakdown_tol = 1.0e-12;
  const Eigen::MatrixXcd Hk = H.topLeftCorner(k + 2, k + 1);
  const double norm_mav = Hk.col(k).norm();
  if (!Hk.allFinite())
  {
    failed = true;
    return;
  }
  if (!(norm_mav > breakdown_tol))
  {
    stopped = true;
    return;
  }
  Eigen::VectorXcd e1 = Eigen::VectorXcd::Zero(k + 2);
  e1(0) = 1.0;
  res = (e1 - Hk * Hk.colPivHouseholderQr().solve(e1)).norm();
  if (!std::isfinite(res))
  {
    failed = true;
  }
  else if (!(norm > breakdown_tol * norm_mav) || res == 0.0)
  {
    stopped = true;
  }
  else
  {
    w *= 1.0 / norm;
  }
}

bool UseReferencePreconditioner(double res_candidate, double res_reference,
                                bool failed_candidate, bool failed_reference,
                                double max_ratio)
{
  return failed_candidate ||
         (!failed_reference && res_reference < std::pow(res_candidate, max_ratio));
}

template <typename VecType>
GmresProbeResult
CompareWithGmresProbe(MPI_Comm comm, const VecType &v0,
                      const std::function<void(bool)> &select,
                      const std::function<void(const VecType &, VecType &)> &apply,
                      const GmresProbeOptions &opts)
{
  ProbeGmres<VecType> run_c(comm, v0, 2 * opts.max_its), run_r(comm, v0, 2 * opts.max_its);
  const double q = opts.max_ratio, d = opts.log_tol;
  while (true)
  {
    if (!run_c.Stopped())
    {
      select(false);
      run_c.Step(apply);
    }
    if (!run_r.Stopped())
    {
      select(true);
      run_r.Step(apply);
    }
    const double r_c = run_c.Residual(), r_r = run_r.Residual();
    const double l_c = -std::log(std::max(r_c, std::numeric_limits<double>::min())),
                 l_r = -std::log(std::max(r_r, std::numeric_limits<double>::min()));
    const int its = std::max(run_c.Iterations(), run_r.Iterations());
    const bool keep = its >= opts.min_keep_its && l_c > d && l_r + d < q * (l_c - d);
    const bool reference = l_r - d > q * (l_c + d);
    const double lo = std::min(1.0, q * q), hi = std::max(1.0, q * q);
    const bool ambiguous =
        l_r >= lo * l_c && l_r <= hi * l_c && l_r >= opts.min_log_decrease;
    bool flags[5] = {run_c.Stopped() && !run_c.Failed(), run_c.Failed(),
                     run_r.Stopped() && !run_r.Failed(), run_r.Failed(),
                     keep || reference || r_c <= opts.tol || r_r <= opts.tol ||
                         (its >= opts.max_its && !ambiguous) || its >= 2 * opts.max_its};
    Mpi::GlobalOr(5, flags, comm);
    run_c.SetFlags(flags[0], flags[1]);
    run_r.SetFlags(flags[2], flags[3]);
    if (flags[4] || flags[1] || flags[3] || (run_c.Stopped() && run_r.Stopped()))
    {
      break;
    }
  }
  GmresProbeResult result;
  result.its = run_c.Iterations();
  result.its_reference = run_r.Iterations();
  result.res_candidate = run_c.Residual();
  result.res_reference = run_r.Residual();
  result.failed_candidate = run_c.Failed();
  result.failed_reference = run_r.Failed();
  result.use_reference = UseReferencePreconditioner(
      result.res_candidate, result.res_reference, result.failed_candidate,
      result.failed_reference, opts.max_ratio);
  Mpi::GlobalOr(1, &result.use_reference, comm);
  result.ratio = (result.res_candidate > 0.0 && result.res_candidate < 1.0 &&
                  result.res_reference > 0.0)
                     ? std::log(result.res_reference) / std::log(result.res_candidate)
                     : std::nan("");
  select(result.use_reference);
  return result;
}

template class ProbeGmres<Vector>;
template class ProbeGmres<ComplexVector>;

template GmresProbeResult
CompareWithGmresProbe(MPI_Comm, const Vector &, const std::function<void(bool)> &,
                      const std::function<void(const Vector &, Vector &)> &,
                      const GmresProbeOptions &);
template GmresProbeResult
CompareWithGmresProbe(MPI_Comm, const ComplexVector &, const std::function<void(bool)> &,
                      const std::function<void(const ComplexVector &, ComplexVector &)> &,
                      const GmresProbeOptions &);

}  // namespace palace
