// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_MAGNETOSTATIC_HPP
#define PALACE_TEST_MMS_MAGNETOSTATIC_HPP

// Shared declarations for the analytic and Enzyme magnetostatic MMS test backends.

#include <array>
#include <functional>

#include <mfem.hpp>

namespace palace::testing::mms::magnetostatic
{

using VectorFunction = std::function<void(const mfem::Vector &, mfem::Vector &)>;

// A non-unit, anisotropic μ makes a mismatch between the operator's μ⁻¹ and the one
// applied to the right-hand side visible.
inline constexpr std::array<double, 3> kAnisotropicMuR = {2.0, 3.0, 5.0};

// A manufactured vector potential with n×A = 0 on the boundary, its flux density
// B = ∇×A, and the constant diagonal relative permeability used by the magnetostatic
// solve.
struct MmsCase
{
  VectorFunction potential;
  VectorFunction magnetic_flux_density;
  std::array<double, 3> mu_r;
};

namespace analytic
{

const MmsCase &PolynomialCase();
const MmsCase &TrigCase();

}  // namespace analytic

void CheckPolynomialExactness(const MmsCase &polynomial);
void CheckOptimalConvergence(const MmsCase &smooth);
void CheckBackendAgreement(const MmsCase &candidate, const MmsCase &reference);

}  // namespace palace::testing::mms::magnetostatic

#endif  // PALACE_TEST_MMS_MAGNETOSTATIC_HPP
