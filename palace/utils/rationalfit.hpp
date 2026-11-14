// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_RATIONAL_FIT_HPP
#define PALACE_UTILS_RATIONAL_FIT_HPP

#include <array>
#include <complex>
#include <functional>
#include <vector>

namespace palace::utils
{

//
// Rational fits of the scalar functions f(λ) of the frequency-dependent terms f(λ) B of an
// eigenvalue problem which are not rational, for the linearized eigensolvers (see
// RationalTerms). The fit is made and checked in a window of the search, in the
// nondimensional angular frequency ω = λ / i: the quarter disc
//   ω = ω_t + ρ e^{iφ},  0 ≤ ρ ≤ R,  -m ≤ φ ≤ π/2 + m ,
// of the damped modes above the target ω_t (within the margin angle m, also the modes
// slightly below it and slightly growing).
//
struct FitWindow
{
  double omega_t, radius, margin = 0.1;
  bool Contains(std::complex<double> omega) const;
};

// A fit f(λ) ≈ poly[0] + poly[1] λ + poly[2] λ² + Σⱼ residues[j] / (λ - poles[j]), valid in
// the window except within the distance `strip` (in ω) of the singular points of f.
struct ScalarFit
{
  std::array<std::complex<double>, 3> poly = {0.0, 0.0, 0.0};
  std::vector<std::complex<double>> poles, residues;
  bool success = false;
  // The evaluations of f in the fit, and the points of its check.
  int num_samples = 0, num_checks = 0;
  // The largest error of the fit on the check points, relative to the largest |f| there.
  double error = 0.0;
  // The half-width of the strip around the singular points where the fit is not valid
  // (zero if there is none in the window), and the singular points.
  double strip = 0.0;
  std::vector<std::complex<double>> singular;

  std::complex<double> Eval(std::complex<double> l) const;

  // Whether the fit is valid at ω: in the window, and not in the strip.
  bool IsValid(const FitWindow &window, std::complex<double> omega) const;
};

// Fit f(λ) = α λᵐ g(λ), m = -1, 0, 1, or 2, in the window. g is analytic away from the
// singular points (values of ω, the branch points first, then dense along each branch
// cut). The fit error is analytic in the window outside a strip around the singular points
// when no pole of the fit is there, so its largest value there is on the boundary of that
// region. g is fitted with AAA on samples of this boundary, at most a quarter of their
// distance to the singular points apart, and checked on 8 times as many points of it,
// offset from the samples, and on a grid inside. The strip is the narrowest of 0 (only
// without singular points in the window) and R / 1024 to R / 8 for which the fit has at
// most max_poles poles and an error of at most tol.
ScalarFit FitScalar(const std::function<std::complex<double>(std::complex<double>)> &g,
                    std::complex<double> alpha, int m,
                    const std::vector<std::complex<double>> &singular,
                    const FitWindow &window, double tol, int max_poles = 24);

// Points (as ω = λ / i) at most `spacing` apart along the branch cut of √e(λ), for
// e(λ) = 1 + Σₖ rₖ / (λ - pₖ): the set where e is real and not positive, the roots of
// D(λ) + u P(λ) for u ∈ [0, 1], with D = Πₖ (λ - pₖ) and P = (e - 1) D. It runs from the
// zeros of e (u = 1) to its poles (u = 0), the branch points, which are the first points.
std::vector<std::complex<double>> SquareRootCut(const std::vector<std::complex<double>> &p,
                                                const std::vector<std::complex<double>> &r,
                                                double spacing);

}  // namespace palace::utils

#endif  // PALACE_UTILS_RATIONAL_FIT_HPP
