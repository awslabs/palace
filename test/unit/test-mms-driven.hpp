// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_DRIVEN_HPP
#define PALACE_TEST_MMS_DRIVEN_HPP

// Shared declarations for the driven (frequency-domain) MMS test backends.

#include <array>
#include <functional>

#include <mfem.hpp>

namespace palace::testing::mms::driven
{

using VectorFunction = std::function<void(const mfem::Vector &, mfem::Vector &)>;

// Distinct anisotropic diagonal μ_r and ε_r, so that swapping them, dropping either, or
// applying one component to the wrong direction changes the right-hand side.
inline constexpr std::array<double, 3> kAnisotropicMuR = {2.0, 3.0, 5.0};
inline constexpr std::array<double, 3> kAnisotropicEpsilonR = {4.0, 2.0, 3.0};

// Nondimensional angular frequency. The mass term ω²εE is then of the same size as the
// curl-curl term ∇×(μ⁻¹∇×E) for the trigonometric field (both O(1) to O(10)), so an error
// in either term is visible.
//
// ω² lies below the first PEC cavity eigenvalue of the unit cube. For diagonal tensors
// the modes have E_x ∝ cos mπx sin nπy sin lπz (and cyclically) with wave vector
// k = π(m, n, l), and ω² an eigenvalue of ε^(−1/2) [k×]ᵀ μ⁻¹ [k×] ε^(−1/2). If one index
// is zero, only the component along that axis survives. The lowest mode is
// E = x̂ sin πy sin πz with
//   ω₁² = π² (1/μ_y + 1/μ_z) / ε_x = π² (1/3 + 1/5) / 4 ≈ 1.32
// (next: the lower (1, 1, 1) mode at ≈ 2.37, then ẑ sin πx sin πy at ≈ 2.74), so
// ω²/ω₁² ≈ 0.76 and the error constant is amplified by roughly 1/(1 − 0.76) ≈ 4.
// The system K − ω²M is nonetheless indefinite: K vanishes on the discrete gradients in
// ND_p, where the operator is −ω²M. That gradient block is what the shifted AMS
// preconditioner handles, so it is exercised at any ω > 0. Placing ω² above ω₁² would
// only add a few negative curl modes, at the cost of a worse-conditioned solve.
inline constexpr double kOmega = 1.0;

// A manufactured electric field with n×E = 0 on the boundary, its curl ∇×E, the
// curl-curl term ∇×(μ⁻¹∇×E) of the strong-form source (for this case's μ_r), and the
// constant diagonal material properties and frequency of the driven solve.
struct MmsCase
{
  VectorFunction electric_field;
  VectorFunction curl;
  VectorFunction curl_curl;
  std::array<double, 3> mu_r;
  std::array<double, 3> epsilon_r;
  double omega;
};

namespace analytic
{

const MmsCase &PolynomialCase();
const MmsCase &TrigCase();

}  // namespace analytic

// Test 1 (mass projection): only electric_field, curl, and epsilon_r are used.
void CheckMassProjectionExactness(const MmsCase &polynomial);
void CheckMassProjectionConvergence(const MmsCase &smooth);

// Test 3 (curl-curl + mass).
void CheckPolynomialExactness(const MmsCase &polynomial);
void CheckOptimalConvergence(const MmsCase &smooth);
void CheckBackendAgreement(const MmsCase &candidate, const MmsCase &reference);

}  // namespace palace::testing::mms::driven

#endif  // PALACE_TEST_MMS_DRIVEN_HPP
