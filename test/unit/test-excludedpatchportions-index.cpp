// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <tuple>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "models/surfaceresponseoperator.hpp"

using namespace palace;

namespace
{

using UncoveredPortionData =
    config::ElectrostaticSolverData::ResponseCorrectionData::UncoveredPortionData;

// = Norm(Subtract(a, b)) of surfaceresponseoperator.cpp in the same expression shape (the
// same rounding / contraction), so the reference's lengths and tolerance compares are the
// implementation's to the bit.
double Distance(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
  const std::array<double, 3> d = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

// The former form of AppendExcludedPatchPortion (main f7912d28cd,
// surfaceresponseoperator.cpp:27080-27113): the linear scan over EVERY portion for the
// first co-located interval of the same feature. The reference the indexed form must equal.
void AppendExcludedPatchPortionLinear(ExcludedPatchPortions &result, std::size_t patch_idx,
                                      const UncoveredPortionData &portion,
                                      bool translational, double tolerance)
{
  const double length = Distance(portion.p1, portion.p0);
  if (length <= tolerance)
  {
    return;
  }
  auto SameInterval = [&](const UncoveredPortionData &a)
  {
    return (Distance(a.p0, portion.p0) <= tolerance &&
            Distance(a.p1, portion.p1) <= tolerance) ||
           (Distance(a.p0, portion.p1) <= tolerance &&
            Distance(a.p1, portion.p0) <= tolerance);
  };
  for (std::size_t i = 0; i < result.portions.size(); i++)
  {
    if (result.portions[i].feature == portion.feature && SameInterval(result.portions[i]))
    {
      result.patches[i].push_back(patch_idx);
      result.duplicate_patches++;
      return;
    }
  }
  result.portions.push_back(portion);
  result.patches.push_back({patch_idx});
  result.translational.push_back(translational);
  result.length += length;
  if (translational)
  {
    result.geometric_cells++;
  }
}

void CheckEqual(const ExcludedPatchPortions &indexed, const ExcludedPatchPortions &linear)
{
  REQUIRE(indexed.portions.size() == linear.portions.size());
  REQUIRE(indexed.patches.size() == linear.patches.size());
  REQUIRE(indexed.translational.size() == linear.translational.size());
  for (std::size_t i = 0; i < indexed.portions.size(); i++)
  {
    CHECK(indexed.portions[i].feature == linear.portions[i].feature);
    CHECK(indexed.portions[i].topology == linear.portions[i].topology);
    CHECK(indexed.portions[i].segment == linear.portions[i].segment);
    CHECK(indexed.portions[i].p0 == linear.portions[i].p0);
    CHECK(indexed.portions[i].p1 == linear.portions[i].p1);
    CHECK(indexed.patches[i] == linear.patches[i]);
    CHECK(indexed.translational[i] == linear.translational[i]);
  }
  CHECK(indexed.geometric_cells == linear.geometric_cells);
  CHECK(indexed.duplicate_patches == linear.duplicate_patches);
  // The length is a sum in insertion order on both sides: bit-equal.
  CHECK(indexed.length == linear.length);
}

}  // namespace

// Lane fix-response-memory (decisions 1439 / 1446; DESIGN s2.6 / s3 (M3), commit C2): the
// feature-indexed co-located-interval scan of AppendExcludedPatchPortion gives the same
// result as the former linear scan over every portion, to the bit: the same first match,
// the same portions / patches / translational / counters / length.
TEST_CASE("AppendExcludedPatchPortion feature index equals the linear scan",
          "[surfaceresponse][memfix][Serial]")
{
  const double tolerance = 1.0e-3;
  // Co-located intervals across features (the same interval appears under features 3 and
  // 7: distinct entries), reversed intervals (p0 <-> p1: the same entry), near-duplicates
  // within and beyond the tolerance, zero-length portions (dropped), translational cells,
  // and a long random tail with many repeats.
  std::vector<std::tuple<std::size_t, UncoveredPortionData, bool>> appends;
  auto Portion = [](int feature, const char *topology, int segment, double x0, double y0,
                    double x1, double y1)
  {
    UncoveredPortionData portion;
    portion.feature = feature;
    portion.topology = topology;
    portion.segment = segment;
    portion.p0 = {x0, y0, 0.0};
    portion.p1 = {x1, y1, 0.0};
    return portion;
  };
  appends.emplace_back(0, Portion(3, "A", 0, 0.0, 0.0, 1.0, 0.0), true);
  appends.emplace_back(1, Portion(7, "B", 0, 0.0, 0.0, 1.0, 0.0), false);  // other feature
  appends.emplace_back(2, Portion(3, "A", 0, 1.0, 0.0, 0.0, 0.0), true);   // reversed: dup
  appends.emplace_back(3, Portion(3, "A", 1, 0.0, 0.0004, 1.0, 0.0), false);  // within tol
  appends.emplace_back(4, Portion(3, "A", 1, 0.0, 0.002, 1.0, 0.0), false);   // beyond tol
  appends.emplace_back(5, Portion(3, "A", 2, 2.0, 0.0, 2.0, 0.0), true);      // zero length
  appends.emplace_back(6, Portion(7, "B", 0, 1.0, 0.0, 0.0, 0.0), true);      // dup of 1
  appends.emplace_back(7, Portion(-1, "C", 5, 5.0, 5.0, 6.0, 5.0), false);
  appends.emplace_back(8, Portion(-1, "C", 5, 6.0, 5.0, 5.0, 5.0), false);  // dup of 7
  std::mt19937 generator(20261109);
  std::uniform_int_distribution<int> feature(0, 11), cell(0, 24), flip(0, 1);
  for (std::size_t patch = 9; patch < 2000; patch++)
  {
    const int f = feature(generator);
    const double x = 0.25 * cell(generator), y = 0.1 * f;
    auto portion = flip(generator) ? Portion(f, "R", f, x, y, x + 0.25, y)
                                   : Portion(f, "R", f, x + 0.25, y, x, y);
    if (flip(generator) && flip(generator))
    {
      // A tolerance-scale perturbation (within tol / 2): still the same interval.
      portion.p0[1] += 0.0004;
    }
    appends.emplace_back(patch, portion, flip(generator) == 1);
  }

  ExcludedPatchPortions indexed, linear;
  for (const auto &[patch, portion, translational] : appends)
  {
    AppendExcludedPatchPortion(indexed, patch, portion, translational, tolerance);
    AppendExcludedPatchPortionLinear(linear, patch, portion, translational, tolerance);
  }
  CheckEqual(indexed, linear);
  CHECK(indexed.duplicate_patches > 100);
  CHECK(indexed.portions.size() > 50);
  CHECK(indexed.indexed_portions == indexed.portions.size());

  // Portions pushed directly (the DomainBoundary ledger's geometry-fallback cell portions
  // form) are found by later appends: the index catches up.
  {
    const auto direct = Portion(42, "Direct", 9, 9.0, 9.0, 10.0, 9.0);
    for (auto *portions : {&indexed, &linear})
    {
      portions->portions.push_back(direct);
      portions->patches.push_back({});
      portions->translational.push_back(true);
      portions->geometric_cells++;
      portions->length += 1.0;
    }
    AppendExcludedPatchPortion(indexed, 5000, Portion(42, "Direct", 9, 10.0, 9.0, 9.0, 9.0),
                               false, tolerance);
    AppendExcludedPatchPortionLinear(
        linear, 5000, Portion(42, "Direct", 9, 10.0, 9.0, 9.0, 9.0), false, tolerance);
    CheckEqual(indexed, linear);
    CHECK(indexed.patches.back() == std::vector<std::size_t>{5000});
  }

  // A copy keeps a consistent index (the struct is copied into the operator's members).
  {
    ExcludedPatchPortions copy = indexed;
    AppendExcludedPatchPortion(copy, 6000, Portion(3, "A", 0, 0.0, 0.0, 1.0, 0.0), false,
                               tolerance);
    AppendExcludedPatchPortionLinear(linear, 6000, Portion(3, "A", 0, 0.0, 0.0, 1.0, 0.0),
                                     false, tolerance);
    CheckEqual(copy, linear);
  }
}
