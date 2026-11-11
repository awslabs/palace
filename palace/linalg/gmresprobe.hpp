// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GMRES_PROBE_HPP
#define PALACE_LINALG_GMRES_PROBE_HPP

#include <functional>
#include <vector>
#include <Eigen/Dense>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// One GMRES run of a probe of a preconditioned operator: the Arnoldi process from a unit
// start vector, orthogonalized by classical Gram-Schmidt applied twice, and the relative
// residual of the least squares problem min ‖e₁ - H y‖. The run stops (breakdown) when
// M A v vanishes for the unit vector v (a start vector in the null space; the residual is
// unchanged) or when the new direction is in the Krylov space to round-off, and fails on
// non-finite values.
//
template <typename VecType>
class ProbeGmres
{
private:
  MPI_Comm comm;
  std::vector<VecType> V;
  Eigen::MatrixXcd H;
  int its = 0;
  double res = 1.0;
  bool stopped = false, failed = false;

public:
  ProbeGmres(MPI_Comm comm, const VecType &v0, int max_its);

  int Iterations() const { return its; }
  double Residual() const { return res; }

  // The run cannot continue: it broke down (its residual is final) or failed.
  bool Stopped() const { return stopped || failed; }
  bool Failed() const { return failed; }

  // Set the flags, for example to their reduction over the processes.
  void SetFlags(bool stopped, bool failed)
  {
    this->stopped = stopped;
    this->failed = failed;
  }

  // One iteration, where apply(x, y) computes y = M A x.
  void Step(const std::function<void(const VecType &, VecType &)> &apply);
};

//
// Comparison of a candidate preconditioner with a reference one by interleaved GMRES runs
// of the preconditioned operator (see CompareWithGmresProbe). With the relative residuals
// r_c and r_r of the two runs after the same iterations, the candidate is predicted to need
// ratio = ln(r_r) / ln(r_c) times the iterations of the reference, and the reference is
// used when ratio > max_ratio, that is r_r < r_c^max_ratio (which also holds for zero or
// unit residuals). The probe stops at the first of:
//  - the prediction stays on the same side of max_ratio when both logarithms of the
//    residuals change by log_tol (about their spread between random start vectors); for
//    the candidate only after min_keep_its iterations, since the prediction tends to grow
//    with the iterations when the runs converge slowly;
//  - a residual decreased by tol;
//  - max_its iterations, unless the prediction is ambiguous, between 1 and max_ratio², and
//    the reference residual decreased by at least min_log_decrease in logarithm (not total
//    stagnation): the probe then continues while it stays so, up to 2 max_its iterations.
//
struct GmresProbeOptions
{
  double max_ratio = 1.25;
  int max_its = 20;
  double tol = 1.0e-3;
  double log_tol = 0.3;
  int min_keep_its = 3;
  double min_log_decrease = 1.0e-3;
};

struct GmresProbeResult
{
  // Iterations of the candidate run, final residuals, failures, the decision, and the
  // predicted iteration ratio (NaN when not defined).
  int its = 0;
  double res_candidate = 1.0, res_reference = 1.0;
  bool failed_candidate = false, failed_reference = false;
  bool use_reference = false;
  double ratio = 0.0;
};

// Decision of the probe from its final state: the reference if the candidate failed, the
// candidate if only the reference failed, else from the predicted iteration ratio.
bool UseReferencePreconditioner(double res_candidate, double res_reference,
                                bool failed_candidate, bool failed_reference,
                                double max_ratio);

// Run the probe from the unit start vector v0, where select(reference) makes apply(x, y)
// compute y = M A x with the reference (true) or candidate (false) preconditioner. The
// stopping tests and the decision are reduced over the processes, so that all of them run
// the same iterations and take the same decision. On return, the preconditioner decided
// is selected.
template <typename VecType>
GmresProbeResult
CompareWithGmresProbe(MPI_Comm comm, const VecType &v0,
                      const std::function<void(bool)> &select,
                      const std::function<void(const VecType &, VecType &)> &apply,
                      const GmresProbeOptions &opts = {});

}  // namespace palace

#endif  // PALACE_LINALG_GMRES_PROBE_HPP
