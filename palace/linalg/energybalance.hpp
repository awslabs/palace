// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ENERGY_BALANCE_HPP
#define PALACE_LINALG_ENERGY_BALANCE_HPP

#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// Energy balance of a computed eigenpair (λ, x) of the eigenvalue problem
//
//                    T(λ) x = (K + λ C + λ² M + A2(λ)) x = 0,   λ = iω,
//
// the root λ_EB of the scalar equation xᴴ T(λ_EB) x = 0 closest to λ. For an exact
// eigenpair it is the eigenvalue. Its loss rate, Re{λ_EB} = -Im{ω_EB}, is the power
// dissipated by the eigenvector x over its stored energy, which depends on x only. Errors
// of the eigenvalue solve, for example from inexact linear solves in the spectral
// transformation, perturb the eigenvalue λ with an absolute error shared by its real and
// imaginary parts, so that the relative error of the quality factor Q = |ω| / (2 |Im{ω}|)
// computed from λ is amplified by 2Q, while the relative error of the quality factor of
// λ_EB is that of the eigenvector (to first order) and is not amplified by Q.
//

using ExtraSystemMatrixFunction =
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)>;

// Root of m λ² + c λ + k = 0 closest to λ₀ (the root of c λ + k = 0 if m = 0).
std::complex<double> QuadraticRootNearest(std::complex<double> m, std::complex<double> c,
                                          std::complex<double> k,
                                          std::complex<double> lambda0);

// Root near λ₀ of k + λ c + λ² m + a(λ) = 0, where a(λ) is analytic and varies slowly with
// λ compared to the quadratic terms. The secant method is applied to the fixed-point map λ
// ↦ root of m z² + c z + (k + a(λ)) nearest to λ. Returns the root and the number of
// evaluations of a(λ), negative if the iteration did not converge in max_it evaluations.
std::pair<std::complex<double>, int>
EnergyBalanceRoot(std::complex<double> k, std::complex<double> c, std::complex<double> m,
                  const std::function<std::complex<double>(std::complex<double>)> &a,
                  std::complex<double> lambda0, int max_it = 20);

// Energy balance eigenvalue λ_EB of the computed eigenpair (λ, x) of the quadratic or
// nonlinear eigenvalue problem (K + λ C + λ² M + A2(λ)) x = 0, λ = iω. The damping operator
// C and the frequency-dependent operator function A2 may be null. The vector y is used as
// workspace. If given, num_it returns the number of evaluations of A2 (negative if the root
// did not converge).
std::complex<double> GetEnergyBalanceEigenvalue(
    MPI_Comm comm, const ComplexVector &x, std::complex<double> lambda,
    const ComplexOperator &K, const ComplexOperator *C, const ComplexOperator &M,
    const ExtraSystemMatrixFunction *funcA2, ComplexVector &y, int *num_it = nullptr);

// Energy balance eigenvalue μ_EB = xᴴ K x / xᴴ M x of the computed eigenpair (μ, x) of the
// linear eigenvalue problem (K - μ M) x = 0, μ = -λ² = ω². The vector y is used as
// workspace.
std::complex<double> GetEnergyBalanceEigenvalue(MPI_Comm comm, const ComplexVector &x,
                                                const ComplexOperator &K,
                                                const ComplexOperator &M, ComplexVector &y);

// Quality factor Q = |ω| / (2 |Im{ω}|) of the angular frequency ω (infinite for ω = 0 or
// real ω).
double GetQualityFactor(std::complex<double> omega);

// Relative difference |Q - Q_EB| / Q_EB of the quality factors of the eigenvalue ω and of
// the energy balance ω_EB of its eigenvector, evaluated from the loss rates 1 / Q. It is
// zero when both loss rates are below round-off (lossless modes).
double GetQualityFactorError(std::complex<double> omega, std::complex<double> omega_eb);

}  // namespace palace

#endif  // PALACE_LINALG_ENERGY_BALANCE_HPP
