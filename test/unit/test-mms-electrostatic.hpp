// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_ELECTROSTATIC_HPP
#define PALACE_TEST_MMS_ELECTROSTATIC_HPP

// Shared declarations for the analytic and Enzyme electrostatic MMS test backends.

#include <array>
#include <functional>
#include <vector>

#include <mfem.hpp>

namespace palace::testing::mms::electrostatic
{

using ScalarFunction = std::function<double(const mfem::Vector &)>;
using VectorFunction = std::function<void(const mfem::Vector &, mfem::Vector &)>;

inline constexpr std::array<double, 3> kIsotropicEpsilonR = {1.0, 1.0, 1.0};
inline constexpr std::array<double, 3> kAnisotropicEpsilonR = {2.0, 3.0, 5.0};

// A manufactured potential, its derived quantities, and the constant diagonal material
// used by the electrostatic solve. The Neumann flux is optional.
struct MmsCase
{
  ScalarFunction potential;
  VectorFunction electric_field;
  ScalarFunction charge_density;
  std::array<double, 3> epsilon_r;
  bool nonzero_dirichlet;
  ScalarFunction neumann_flux;
};

namespace analytic
{

const MmsCase &HomogeneousCase();
const MmsCase &NonHomogeneousCase();
const MmsCase &PolynomialCase();
const MmsCase &CurvedCylinderCase();

}  // namespace analytic

void CheckPointwiseAgreement(const MmsCase &enzyme, const MmsCase &analytic,
                             const std::vector<std::array<double, 3>> &points,
                             bool check_neumann = false);
void CheckSmallErrors(const MmsCase &homogeneous, const MmsCase &nonhomogeneous);
void CheckPolynomialExactness(const MmsCase &polynomial);
void CheckNeumannBoundary(const MmsCase &polynomial);
void CheckOptimalConvergence(const MmsCase &nonhomogeneous);
void CheckCurvedConvergence(const MmsCase &curved);

}  // namespace palace::testing::mms::electrostatic

#endif  // PALACE_TEST_MMS_ELECTROSTATIC_HPP
