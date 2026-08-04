// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/hierarchicalerrorestimator.hpp"

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("Hierarchical patch algebra includes domain and boundary contributions",
          "[hierarchicalerrorestimator][Serial]")
{
  using fem::hierarchical::LocalOperatorContribution;
  using fem::hierarchical::SparseColumn;
  LocalOperatorContribution domain, boundary, exterior;
  domain.support_element = boundary.support_element = 0;
  exterior.support_element = 1;
  domain.dofs = boundary.dofs = exterior.dofs = {0, 1};
  domain.matrix.SetSize(2);
  boundary.matrix.SetSize(2);
  exterior.matrix.SetSize(2);
  domain.matrix(0, 0) = 3.0;
  domain.matrix(0, 1) = domain.matrix(1, 0) = -1.0;
  domain.matrix(1, 1) = 2.0;
  boundary.matrix = 0.0;
  boundary.matrix(0, 0) = 0.5;
  boundary.matrix(1, 1) = 1.5;
  exterior.matrix = 0.0;
  exterior.matrix(0, 0) = 100.0;
  exterior.matrix(1, 1) = 100.0;
  domain.rhs.SetSize(2);
  boundary.rhs.SetSize(2);
  exterior.rhs.SetSize(2);
  domain.rhs(0) = 2.0;
  domain.rhs(1) = 1.0;
  boundary.rhs(0) = 0.25;
  boundary.rhs(1) = -0.5;
  exterior.rhs = 0.0;
  const std::vector<LocalOperatorContribution> contributions{domain, boundary, exterior};

  mfem::Vector injected(2);
  injected(0) = 0.2;
  injected(1) = -0.1;
  const auto residual =
      fem::hierarchical::AssembleResidual(2, contributions, injected, {false, false});
  // AssembleResidual includes every contribution, including the exterior element.
  CHECK_THAT(residual(0), WithinAbs(2.25 - 3.5 * 0.2 + 1.0 * -0.1 - 20.0, 1.0e-14));
  CHECK_THAT(residual(1), WithinAbs(0.5 + 1.0 * 0.2 - 3.5 * -0.1 + 10.0, 1.0e-14));

  SparseColumn first{{0}, {1.0}}, second{{0, 1}, {0.5, 1.0}};
  mfem::DenseMatrix restricted;
  mfem::Vector restricted_rhs;
  const long long touched = fem::hierarchical::AssembleRestrictedOperator(
      contributions, {0}, {first, second}, residual, restricted, restricted_rhs);
  // The support filter excludes exterior while retaining both the domain and boundary
  // facet.
  mfem::DenseMatrix expected(2);
  expected(0, 0) = 3.5;
  expected(0, 1) = expected(1, 0) = 0.75;
  expected(1, 1) = 3.375;
  CHECK(touched == 8);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK_THAT(restricted(i, j), WithinAbs(expected(i, j), 1.0e-14));
    }
  }

  mfem::Vector vector(2);
  vector(0) = -0.3;
  vector(1) = 0.7;
  const double energy = fem::hierarchical::Energy({domain, boundary}, vector);
  CHECK_THAT(energy, WithinAbs(2.45, 1.0e-14));
}

}  // namespace palace
