// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>

#include "fem/singularassembly.hpp"
#include "fem/singularelements.hpp"

namespace palace
{

namespace
{

using Basis = fem::singular::HigherOrderBasis;
using BasisFamily = fem::singular::HigherOrderBasisFamily;
using BasisKey =
    std::tuple<int, std::array<int, 4>, fem::singular::InterpolationIndices, int>;
using PairKey = std::pair<BasisKey, BasisKey>;

enum class SingularFeature
{
  NODE,
  EDGE
};

BasisKey MakeKey(const Basis &basis)
{
  return {static_cast<int>(basis.family), basis.nodes, basis.interpolation_indices,
          basis.order};
}

bool IsGradient(const Basis &basis)
{
  return basis.family == BasisFamily::NODE_GRADIENT ||
         basis.family == BasisFamily::EDGE_GRADIENT;
}

std::vector<Basis> EnumerateBases(SingularFeature feature, int order, double nu)
{
  constexpr std::array<int, 4> Nodes{0, 1, 2, 3};
  std::vector<Basis> bases;
  const auto append = [&](const std::vector<Basis> &entries)
  { bases.insert(bases.end(), entries.begin(), entries.end()); };
  if (feature == SingularFeature::NODE)
  {
    append(fem::singular::EnumerateHigherOrderNodeGradientBases(Nodes, order, nu));
    append(fem::singular::EnumerateHigherOrderNodeRotationalBases(Nodes, order, nu));
  }
  else
  {
    append(fem::singular::EnumerateHigherOrderEdgeGradientBases(Nodes, order, nu));
    append(fem::singular::EnumerateHigherOrderEdgeRotationalBases(Nodes, order, nu));
  }
  return bases;
}

std::map<PairKey, std::pair<Basis, Basis>> CanonicalPairs(const std::vector<Basis> &bases,
                                                          std::size_t &retained_pair_count)
{
  std::map<PairKey, std::pair<Basis, Basis>> pairs;
  retained_pair_count = 0;
  for (std::size_t row = 0; row < bases.size(); row++)
  {
    for (std::size_t column = row; column < bases.size(); column++)
    {
      retained_pair_count++;
      const auto canonical = fem::singular::CanonicalizeReferenceBasisPair(
          fem::singular::ReferenceBasis{bases[row]},
          fem::singular::ReferenceBasis{bases[column]});
      const auto &canonical_row = std::get<Basis>(canonical.row_basis);
      const auto &canonical_column = std::get<Basis>(canonical.column_basis);
      const PairKey key{MakeKey(canonical_row), MakeKey(canonical_column)};
      const auto [entry, inserted] =
          pairs.emplace(key, std::pair{canonical_row, canonical_column});
      if (!inserted)
      {
        REQUIRE(MakeKey(entry->second.first) == MakeKey(canonical_row));
        REQUIRE(MakeKey(entry->second.second) == MakeKey(canonical_column));
      }
    }
  }
  return pairs;
}

double CertificationRadialPower(double nu)
{
  if (nu == 0.5)
  {
    return 2.0;
  }
  REQUIRE(nu == 2.0 / 3.0);
  return 3.0;
}

int ExactDuffyOrder(SingularFeature feature, int singular_order, double nu,
                    double radial_power)
{
  // The minimal grading q = 2 for nu = 1/2 and q = 3 for nu = 2/3 makes every
  // transformed integrand polynomial without placing high-order node points
  // needlessly close to the singular feature. The largest radial degree comes
  // from the regular part of a rotational-basis mass product. Including the
  // Duffy Jacobian gives q*(2*s + 2*nu + 3) - 1 at a node and
  // q*(2*s + 2*nu + 7) - 1 at an edge. Angular degrees are lower. MFEM's
  // segment rule parameter is its polynomial exactness order.
  const double degree = radial_power * (2 * singular_order + 2 * nu +
                                        (feature == SingularFeature::NODE ? 3.0 : 7.0)) -
                        1.0;
  REQUIRE(std::abs(degree - std::round(degree)) <=
          16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(degree)));
  return static_cast<int>(std::lround(degree));
}

fem::singular::ReferenceIntegral Integrate(SingularFeature feature,
                                           const fem::singular::ReferenceBasis &row,
                                           const fem::singular::ReferenceBasis &column,
                                           int order, double radial_power)
{
  if (feature == SingularFeature::NODE)
  {
    return fem::singular::ComputeNodeDuffyReferenceIntegral(row, column, order,
                                                            radial_power);
  }
  return fem::singular::ComputeEdgeDuffyReferenceIntegral(row, column, order, radial_power);
}

void CertifyHigherOrderTensors(SingularFeature feature, int singular_order)
{
  constexpr int OrderSeparation = 8;
  constexpr int AuditOrderIncrement = 4;
  constexpr double ErrorSafetyFactor = 8.0;
  constexpr double RelativeCertificationTolerance = 5.0e-12;

  std::size_t retained_pairs_total = 0;
  std::size_t canonical_pairs_total = 0;
  std::size_t audited_coefficients = 0;
  int minimum_comparison_order = std::numeric_limits<int>::max();
  int maximum_audit_order = 0;

  for (const double nu : {0.5, 2.0 / 3.0})
  {
    const double radial_power = CertificationRadialPower(nu);
    const int reference_order = ExactDuffyOrder(feature, singular_order, nu, radial_power);
    const int comparison_order = reference_order - OrderSeparation;
    const int audit_order = reference_order + AuditOrderIncrement;
    minimum_comparison_order = std::min(minimum_comparison_order, comparison_order);
    maximum_audit_order = std::max(maximum_audit_order, audit_order);
    const auto bases = EnumerateBases(feature, singular_order, nu);
    REQUIRE_FALSE(bases.empty());
    std::size_t retained_pair_count;
    const auto pairs = CanonicalPairs(bases, retained_pair_count);
    REQUIRE_FALSE(pairs.empty());
    REQUIRE(pairs.size() <= retained_pair_count);
    retained_pairs_total += retained_pair_count;
    canonical_pairs_total += pairs.size();

    for (const auto &[key, pair] : pairs)
    {
      const auto key_copy = key;
      const auto &[row, column] = pair;
      const fem::singular::ReferenceBasis row_basis{row};
      const fem::singular::ReferenceBasis column_basis{column};
      const auto comparison =
          Integrate(feature, row_basis, column_basis, comparison_order, radial_power);
      const auto reference =
          Integrate(feature, row_basis, column_basis, reference_order, radial_power);
      const auto audit =
          Integrate(feature, row_basis, column_basis, audit_order, radial_power);
      REQUIRE(reference.quadrature_order == reference_order);
      REQUIRE(reference.radial_power == radial_power);

      const bool curl_free = IsGradient(row) || IsGradient(column);
      for (int u = 0; u < 3; u++)
      {
        for (int v = 0; v < 3; v++)
        {
          const auto check_entry = [&](double comparison_value, double reference_value,
                                       double audit_value, bool must_be_zero,
                                       const char *quantity)
          {
            const double scale =
                std::max({1.0, std::abs(comparison_value), std::abs(reference_value),
                          std::abs(audit_value)});
            const double roundoff = 256.0 * std::numeric_limits<double>::epsilon() * scale;
            const double estimate =
                ErrorSafetyFactor * std::abs(reference_value - comparison_value) + roundoff;
            CAPTURE(static_cast<int>(feature), singular_order, nu, key_copy, u, v, quantity,
                    comparison_order, reference_order, audit_order, comparison_value,
                    reference_value, audit_value, estimate);
            REQUIRE(std::isfinite(comparison_value));
            REQUIRE(std::isfinite(reference_value));
            REQUIRE(std::isfinite(audit_value));
            REQUIRE(std::abs(reference_value - audit_value) <= estimate);
            REQUIRE(std::abs(reference_value - audit_value) <=
                    RelativeCertificationTolerance * scale);
            if (must_be_zero)
            {
              REQUIRE(comparison_value == 0.0);
              REQUIRE(reference_value == 0.0);
              REQUIRE(audit_value == 0.0);
            }
            audited_coefficients++;
          };
          check_entry(comparison.mass[u][v], reference.mass[u][v], audit.mass[u][v], false,
                      "mass");
          check_entry(comparison.curl_curl[u][v], reference.curl_curl[u][v],
                      audit.curl_curl[u][v], curl_free, "curl-curl");
        }
      }
    }
  }

  REQUIRE(audited_coefficients == 18 * canonical_pairs_total);
  std::cout << "Certified singular full-wave tensors: feature = "
            << (feature == SingularFeature::NODE ? "node" : "edge")
            << ", singular order = " << singular_order
            << ", retained pairs = " << retained_pairs_total
            << ", canonical pairs = " << canonical_pairs_total
            << ", tensor coefficients = " << audited_coefficients
            << ", Duffy order range = " << minimum_comparison_order << ".."
            << maximum_audit_order << '\n';
}

}  // namespace

// These exhaustive tests are hidden from the routine unit-test suite because
// their purpose is release certification, not per-change smoke testing. Run
// them explicitly with the [singulartensorcertification] tag.
TEST_CASE("Certify order-2 node singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::NODE, 2);
}

TEST_CASE("Certify order-2 edge singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::EDGE, 2);
}

TEST_CASE("Certify order-3 node singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::NODE, 3);
}

TEST_CASE("Certify order-3 edge singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::EDGE, 3);
}

TEST_CASE("Certify order-4 node singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::NODE, 4);
}

TEST_CASE("Certify order-4 edge singular full-wave tensors",
          "[.singulartensorcertification][singularelements][Serial]")
{
  CertifyHigherOrderTensors(SingularFeature::EDGE, 4);
}

}  // namespace palace
