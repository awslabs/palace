// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Regression cases: full Palace solves diffed against test data references.
// Each case carries [Serial][Parallel][GPU][Regression] (the category
// tag is orthogonal to execution style — every case is valid at any
// rank count and on either CPU or GPU). [Long] modifies a regression
// case to be skipped from the default sweep and run only under the
// long-tests CI workflow.

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "regression_helpers.hpp"
#include "utils/communication.hpp"
#include "utils/tablecsv.hpp"

namespace
{

// ---------------------------------------------------------------------------
// Shared custom-check factories. Each returns a callable matching
// palace::test::CustomCheck so TEST_CASEs can register them per-CSV.
// ---------------------------------------------------------------------------

// Magnitude-only farfield check. Phase of the complex E components is not
// reproducible across MPI partitions, so we compare
//   |E|(row) = sqrt(|Ex|^2 + |Ey|^2 + |Ez|^2)
// at the case's rtol. The six trailing columns are r*Re{E_x}, r*Im{E_x}, ...,
// r*Im{E_z} for both driven and eigenmode output (the leading columns differ).
palace::test::CustomCheck TestFarfield(double rtol)
{
  return [rtol](palace::Table &a, palace::Table &r)
  {
    CHECK(a.n_cols() >= 6);
    CHECK(a.n_cols() == r.n_cols());
    if (a.n_cols() < 6 || a.n_cols() != r.n_cols())
    {
      return;
    }
    const std::size_t e0 = a.n_cols() - 6;
    const std::size_t n_rows = std::min(a.n_rows(), r.n_rows());
    for (std::size_t i = 0; i < n_rows; ++i)
    {
      auto mag_sq = [&](palace::Table &t)
      {
        auto sq = [](double x, double y) { return x * x + y * y; };
        return sq(t[e0].data[i], t[e0 + 1].data[i]) +
               sq(t[e0 + 2].data[i], t[e0 + 3].data[i]) +
               sq(t[e0 + 4].data[i], t[e0 + 5].data[i]);
      };
      const double m_a = std::sqrt(mag_sq(a));
      const double m_r = std::sqrt(mag_sq(r));
      INFO("row " << i + 1 << ": |E| actual=" << m_a << " ref=" << m_r);
      const double thresh = rtol * std::max(m_a, m_r);
      CHECK(std::abs(m_a - m_r) <= thresh);
    }
  };
}

// Generalised magnitude comparison for any CSV with paired
// "Re{X} (unit)" / "Im{X} (unit)" complex columns. For each Re column
// it pairs the matching Im (looked up by header_text), then compares
// |X| = sqrt(Re^2 + Im^2) at the case's tolerances. Real-only columns
// (frequency, V_inc, etc.) fall through to the same WithinRel ||
// WithinAbs check the default comparator uses. Used for adaptive
// frequency sweeps where the ROM-interpolated phase wobbles at zero
// crossings but |X| tracks the converged response.
palace::test::CustomCheck CompareComplexMagnitudes(double rtol, double atol)
{
  return [rtol, atol](palace::Table &actual, palace::Table &reference)
  {
    auto find_by_header = [](palace::Table &t, const std::string &h) -> palace::Column *
    {
      for (auto &c : t)
      {
        if (c.header_text == h)
        {
          return &c;
        }
      }
      return nullptr;
    };

    const std::size_t n_cols = std::min(actual.n_cols(), reference.n_cols());
    const std::size_t n_rows = std::min(actual.n_rows(), reference.n_rows());
    for (std::size_t c = 0; c < n_cols; ++c)
    {
      const std::string &hdr = actual[c].header_text;
      INFO("column " << c << ": '" << hdr << "'");
      CHECK(hdr == reference[c].header_text);

      // Im{...} columns are checked alongside their Re partner; skip
      // them here so we don't double-count or compare the imaginary
      // half on its own.
      if (hdr.rfind("Im{", 0) == 0)
      {
        continue;
      }

      if (hdr.rfind("Re{", 0) == 0)
      {
        std::string im_hdr = hdr;
        im_hdr.replace(0, 3, "Im{");  // "Re{X} (unit)" -> "Im{X} (unit)"
        palace::Column *im_a = find_by_header(actual, im_hdr);
        palace::Column *im_r = find_by_header(reference, im_hdr);
        CHECK(im_a != nullptr);
        CHECK(im_r != nullptr);
        if (!im_a || !im_r)
        {
          continue;
        }
        for (std::size_t r = 0; r < n_rows; ++r)
        {
          const double re_av = actual[c].data[r];
          const double im_av = im_a->data[r];
          const double re_rv = reference[c].data[r];
          const double im_rv = im_r->data[r];
          const double mag_a = std::sqrt(re_av * re_av + im_av * im_av);
          const double mag_r = std::sqrt(re_rv * re_rv + im_rv * im_rv);
          INFO("row " << r + 1 << " |" << hdr.substr(3, hdr.find('}', 3) - 3) << "|");
          CHECK_THAT(mag_a, Catch::Matchers::WithinRel(mag_r, rtol) ||
                                Catch::Matchers::WithinAbs(mag_r, atol));
        }
        continue;
      }

      // Real-valued column (e.g. f (GHz), V_inc, |S| (dB)): compare
      // values directly with the same matchers.
      for (std::size_t r = 0; r < n_rows; ++r)
      {
        CHECK_THAT(actual[c].data[r],
                   Catch::Matchers::WithinRel(reference[c].data[r], rtol) ||
                       Catch::Matchers::WithinAbs(reference[c].data[r], atol));
      }
    }
  };
}

// Compare a canonical-keyed diagnostic against a reference subset. The final
// `numeric_columns` columns are compared numerically; preceding columns form
// the lookup key except for explicitly ignored partition-dependent columns.
palace::test::CustomCheck CompareCanonicalRows(std::size_t numeric_columns,
                                               std::set<std::string> ignored_key_headers,
                                               double rtol, double atol,
                                               bool compare_complex_magnitude = false)
{
  return [numeric_columns, ignored_key_headers = std::move(ignored_key_headers), rtol, atol,
          compare_complex_magnitude](palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    REQUIRE(actual.n_cols() >= numeric_columns);
    const std::size_t numeric_begin = actual.n_cols() - numeric_columns;
    for (std::size_t c = 0; c < actual.n_cols(); c++)
    {
      CHECK(actual[c].header_text == reference[c].header_text);
    }
    for (std::size_t rr = 0; rr < reference.n_rows(); rr++)
    {
      std::vector<std::size_t> matches;
      for (std::size_t ar = 0; ar < actual.n_rows(); ar++)
      {
        bool match = true;
        for (std::size_t c = 0; c < numeric_begin; c++)
        {
          if (ignored_key_headers.count(actual[c].header_text) == 0 &&
              actual[c].data[ar] != reference[c].data[rr])
          {
            match = false;
            break;
          }
        }
        if (match)
        {
          matches.push_back(ar);
        }
      }
      INFO("canonical reference row " << rr + 1);
      REQUIRE(matches.size() == 1);
      const std::size_t ar = matches.front();
      if (compare_complex_magnitude)
      {
        REQUIRE(numeric_columns == 2);
        const double actual_magnitude =
            std::hypot(actual[numeric_begin].data[ar], actual[numeric_begin + 1].data[ar]);
        const double reference_magnitude = std::hypot(
            reference[numeric_begin].data[rr], reference[numeric_begin + 1].data[rr]);
        CHECK_THAT(actual_magnitude,
                   Catch::Matchers::WithinRel(reference_magnitude, rtol) ||
                       Catch::Matchers::WithinAbs(reference_magnitude, atol));
        continue;
      }
      for (std::size_t c = numeric_begin; c < actual.n_cols(); c++)
      {
        INFO("column '" << actual[c].header_text << "'");
        CHECK_THAT(actual[c].data[ar],
                   Catch::Matchers::WithinRel(reference[c].data[rr], rtol) ||
                       Catch::Matchers::WithinAbs(reference[c].data[rr], atol));
      }
    }
  };
}

palace::test::CustomCheck CompareCanonicalKeySets(std::size_t numeric_columns,
                                                  std::set<std::string> ignored_key_headers)
{
  return [numeric_columns, ignored_key_headers = std::move(ignored_key_headers)](
             palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    REQUIRE(actual.n_cols() >= numeric_columns);
    const std::size_t numeric_begin = actual.n_cols() - numeric_columns;
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      CHECK(actual[column].header_text == reference[column].header_text);
    }
    const auto Keys = [&](palace::Table &table)
    {
      std::set<std::vector<double>> keys;
      for (std::size_t row = 0; row < table.n_rows(); row++)
      {
        std::vector<double> key;
        for (std::size_t column = 0; column < numeric_begin; column++)
        {
          if (ignored_key_headers.count(table[column].header_text) == 0)
          {
            key.push_back(table[column].data[row]);
          }
        }
        REQUIRE(keys.insert(std::move(key)).second);
        for (std::size_t column = numeric_begin; column < table.n_cols(); column++)
        {
          CHECK(std::isfinite(table[column].data[row]));
        }
      }
      return keys;
    };
    CHECK(Keys(actual) == Keys(reference));
  };
}

palace::test::CustomCheck CompareBoundaryModeSpectrum(double rtol)
{
  return [rtol](palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    REQUIRE(actual.n_rows() == reference.n_rows());
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      CHECK(actual[column].header_text == reference[column].header_text);
    }
    auto FindColumn = [](palace::Table &table, std::string_view header)
    {
      const auto column =
          std::find_if(table.begin(), table.end(),
                       [header](const auto &entry) { return entry.header_text == header; });
      REQUIRE(column != table.end());
      return static_cast<std::size_t>(std::distance(table.begin(), column));
    };
    const std::size_t kn_real = FindColumn(actual, "Re{kn} (1/m)");
    const std::size_t kn_imag = FindColumn(actual, "Im{kn} (1/m)");
    const std::size_t neff_real = FindColumn(actual, "Re{n_eff}");
    const std::size_t neff_imag = FindColumn(actual, "Im{n_eff}");
    const std::size_t backward_error = FindColumn(actual, "Error (Bkwd.)");
    const std::size_t absolute_error = FindColumn(actual, "Error (Abs.)");
    std::vector<std::size_t> actual_order(actual.n_rows());
    std::vector<std::size_t> reference_order(reference.n_rows());
    std::iota(actual_order.begin(), actual_order.end(), 0);
    std::iota(reference_order.begin(), reference_order.end(), 0);
    const auto SortByPropagationConstant = [kn_real](palace::Table &table, auto &order)
    {
      std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right)
                { return table[kn_real].data[left] < table[kn_real].data[right]; });
    };
    SortByPropagationConstant(actual, actual_order);
    SortByPropagationConstant(reference, reference_order);
    for (std::size_t mode = 0; mode < actual_order.size(); mode++)
    {
      const std::size_t a = actual_order[mode];
      const std::size_t r = reference_order[mode];
      CAPTURE(mode);
      CHECK_THAT(actual[kn_real].data[a],
                 Catch::Matchers::WithinRel(reference[kn_real].data[r], rtol));
      CHECK_THAT(actual[neff_real].data[a],
                 Catch::Matchers::WithinRel(reference[neff_real].data[r], rtol));
      CHECK(std::abs(actual[kn_imag].data[a]) < 1.0e-3);
      CHECK(std::abs(actual[neff_imag].data[a]) < 1.0e-5);
      CHECK(std::abs(actual[backward_error].data[a]) < 1.0e-7);
      CHECK(std::abs(actual[absolute_error].data[a]) < 1.0e-7);
    }
  };
}

palace::test::CustomCheck CompareBoundaryModeSubspaceDiagnostics(double rtol)
{
  return [rtol](palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    REQUIRE(actual.n_rows() == reference.n_rows());
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      CHECK(actual[column].header_text == reference[column].header_text);
      for (double value : actual[column].data)
      {
        CHECK(std::isfinite(value));
      }
    }
    auto FindColumn = [](palace::Table &table, std::string_view header)
    {
      const auto column =
          std::find_if(table.begin(), table.end(),
                       [header](const auto &entry) { return entry.header_text == header; });
      REQUIRE(column != table.end());
      return static_cast<std::size_t>(std::distance(table.begin(), column));
    };
    const std::size_t power_real = FindColumn(actual, "Re{P} (normalized)");
    const std::size_t power_imag = FindColumn(actual, "Im{P} (normalized)");
    const std::size_t electric = FindColumn(actual, "Total electric field energy (J)");
    const std::size_t magnetic = FindColumn(actual, "Total magnetic field energy (J)");
    for (std::size_t mode = 0; mode < actual.n_rows(); mode++)
    {
      CAPTURE(mode);
      CHECK(std::abs(actual[power_real].data[mode] - 1.0) < 1.0e-8);
      CHECK(std::abs(actual[power_imag].data[mode]) < 1.0e-5);
      CHECK(actual[electric].data[mode] > 0.0);
      CHECK(actual[magnetic].data[mode] > 0.0);
      CHECK(std::abs(actual[electric].data[mode] - actual[magnetic].data[mode]) /
                std::max(actual[electric].data[mode], actual[magnetic].data[mode]) <
            2.0e-3);
    }
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      if (actual[column].header_text.find("field energy") == std::string::npos)
      {
        continue;
      }
      const double actual_trace =
          std::accumulate(actual[column].data.begin(), actual[column].data.end(), 0.0);
      const double reference_trace = std::accumulate(reference[column].data.begin(),
                                                     reference[column].data.end(), 0.0);
      INFO("subspace trace of '" << actual[column].header_text << "'");
      CHECK_THAT(actual_trace, Catch::Matchers::WithinRel(reference_trace, rtol) ||
                                   Catch::Matchers::WithinAbs(reference_trace, 1.0e-12));
    }
  };
}

palace::test::CustomCheck CompareBoundaryModeSurfaceSubspace(double rtol)
{
  return [rtol](palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    REQUIRE(actual.n_rows() == reference.n_rows());
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      CHECK(actual[column].header_text == reference[column].header_text);
    }
    REQUIRE(actual.n_cols() == 3);
    double actual_trace = 0.0;
    double reference_trace = 0.0;
    for (std::size_t mode = 0; mode < actual.n_rows(); mode++)
    {
      const double participation = actual[1].data[mode];
      const double quality_factor = actual[2].data[mode];
      CAPTURE(mode, participation, quality_factor);
      CHECK(participation > 0.0);
      CHECK(quality_factor > 0.0);
      CHECK(std::abs(participation * quality_factor - 1.0) < 1.0e-8);
      actual_trace += participation;
      reference_trace += reference[1].data[mode];
    }
    CHECK_THAT(actual_trace, Catch::Matchers::WithinRel(reference_trace, rtol));
  };
}

palace::test::CustomCheck ValidateBoundaryModeTipSlopes()
{
  return [](palace::Table &actual, palace::Table &reference)
  {
    REQUIRE(actual.n_cols() == reference.n_cols());
    for (std::size_t column = 0; column < actual.n_cols(); column++)
    {
      CHECK(actual[column].header_text == reference[column].header_text);
    }
    auto FindColumn = [](palace::Table &table, std::string_view header)
    {
      const auto column =
          std::find_if(table.begin(), table.end(),
                       [header](const auto &entry) { return entry.header_text == header; });
      REQUIRE(column != table.end());
      return static_cast<std::size_t>(std::distance(table.begin(), column));
    };
    const std::size_t valid = FindColumn(actual, "fit_valid");
    const std::size_t expected = FindColumn(actual, "expected_slope");
    const std::size_t fitted = FindColumn(actual, "fitted_slope");
    const std::size_t r_squared = FindColumn(actual, "R_squared");
    for (std::size_t row = 0; row < actual.n_rows(); row++)
    {
      CAPTURE(row);
      CHECK(actual[valid].data[row] == 1.0);
      CHECK(std::abs(actual[fitted].data[row] - actual[expected].data[row]) < 3.0e-2);
      CHECK(actual[r_squared].data[row] > 0.999);
    }
  };
}

// Standard "drop per-element extrema + eigenmode error columns" list.
const std::vector<std::string> kEigenExcluded = {"Maximum", "Minimum", "Mean",
                                                 "Error (Bkwd.)", "Error (Abs.)"};
constexpr auto kForceDefaultSolver = palace::test::SolverOverridePolicy::ForceDefault;

palace::Table LoadRegressionOutputTable(std::string_view case_dir,
                                        std::string_view postpro_subdir,
                                        std::string_view filename)
{
  const auto path =
      palace::test::GetRegressionOutputPath(case_dir, postpro_subdir) / filename;
  palace::TableWithCSVFile wrapped(path.string(), /*load_existing_file=*/true);
  return std::move(wrapped.table);
}

double GetTableValue(palace::Table &table, std::string_view header, std::size_t row = 0)
{
  for (auto &column : table)
  {
    if (column.header_text == header)
    {
      REQUIRE(row < column.data.size());
      return column.data[row];
    }
  }
  FAIL("Missing CSV column '" << header << "'");
  return std::numeric_limits<double>::quiet_NaN();
}

double RelativeDifference(double first, double second)
{
  return std::abs(first - second) / std::max(std::abs(first), std::abs(second));
}

double RelativeComplexDifference(double first_real, double first_imag, double second_real,
                                 double second_imag)
{
  return std::hypot(first_real - second_real, first_imag - second_imag) /
         std::max(std::hypot(first_real, first_imag), std::hypot(second_real, second_imag));
}

// Floquet-port S-parameters: compare only the |S[...]| (dB) magnitude columns
// (phase isn't reproducible). NaN entries (evanescent modes) and signals below
// -200 dB (negligible) count as matches. Ports the Julia test_floquet_sparams.
palace::test::CustomCheck TestFloquetSParams(double rtol, double atol)
{
  return [rtol, atol](palace::Table &actual, palace::Table &reference)
  {
    const std::size_t n_cols = std::min(actual.n_cols(), reference.n_cols());
    const std::size_t n_rows = std::min(actual.n_rows(), reference.n_rows());
    for (std::size_t c = 0; c < n_cols; ++c)
    {
      const std::string &hdr = reference[c].header_text;
      if (hdr.find("|S[") == std::string::npos || hdr.find("(dB)") == std::string::npos)
      {
        continue;
      }
      for (std::size_t r = 0; r < n_rows; ++r)
      {
        const double v_new = actual[c].data[r];
        const double v_ref = reference[c].data[r];
        if ((std::isnan(v_new) && std::isnan(v_ref)) || v_ref < -200.0)
        {
          continue;
        }
        INFO("row " << r + 1 << " column '" << hdr << "'");
        CHECK_THAT(v_new, Catch::Matchers::WithinRel(v_ref, rtol) ||
                              Catch::Matchers::WithinAbs(v_ref, atol));
      }
    }
  };
}

}  // namespace

// ===========================================================================
// 3D cases. Order follows the retired Julia regression harness.
// ===========================================================================

TEST_CASE("singular_sheet_electrostatic", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  // Combined-field coefficients and near-edge samples retain small
  // partition-dependent Krylov/reduction differences.
  opts.rtol = 5.0e-9;
  opts.atol = 1.0e-20;
  opts.excluded_columns = {"true_dof"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_sheet", "singular_sheet.json", "", opts);
}

TEST_CASE("singular_sheet_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_sheet_eigenmode",
                                  "singular_sheet_eigenmode.json", "eigenmode", opts);
}

TEST_CASE("singular_wedge_electrostatic", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-20;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-coefficients.csv"] =
      CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
  opts.custom_checks["singular-edge-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_electrostatic.json",
                                  "electrostatic", opts);
}

TEST_CASE("singular_wedge_electrostatic_multigrid", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-20;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-coefficients.csv"] =
      CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
  opts.custom_checks["singular-edge-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_electrostatic_p2.json",
                                  "electrostatic_p2", opts);
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_electrostatic_mg2.json",
                                  "electrostatic_p2", opts);
}

TEST_CASE("singular_finite_metal_fixed_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 5.0e-12;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_finite_metal_fixed",
                                  "singular_finite_metal_fixed_regression.json",
                                  "eigenmode", opts);
}

TEST_CASE("singular_wedge_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_eigenmode.json",
                                  "eigenmode", opts);
}

TEST_CASE("singular_wedge_eigenmode_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-12;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  if (palace::Mpi::Size(palace::Mpi::World()) > 1)
  {
    // Conforming closure is partition-dependent even for identical stable seed IDs. Keep a
    // decomposition-specific reference instead of weakening one shared tolerance.
    opts.reference_subdir = "amr_eigenmode_mpi";
  }
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_eigenmode_amr.json",
                                  "amr_eigenmode", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    for (const std::string prefix : {"", "iteration1/"})
    {
      auto diagnostics = LoadRegressionOutputTable("singular_wedge", "amr_eigenmode",
                                                   prefix + "singular-eigenmode.csv");
      const double electric = GetTableValue(diagnostics, "Electric field energy (J)");
      const double magnetic = GetTableValue(diagnostics, "Magnetic field energy (J)");
      const double mismatch = GetTableValue(diagnostics, "Relative energy mismatch");
      const double divergence = GetTableValue(diagnostics, "Relative weak divergence");
      CAPTURE(prefix, electric, magnetic, mismatch, divergence);
      CHECK(electric > 0.0);
      CHECK(magnetic > 0.0);
      CHECK(mismatch < 1.0e-10);
      CHECK(divergence < 1.0e-10);
    }
  }
}

TEST_CASE("singular_wedge_loss_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_loss_eigenmode.json",
                                  "loss_eigenmode", opts);
}

TEST_CASE("singular_wedge_loss_eigenmode_ams", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-12;
  opts.excluded_columns = {"Error (", "Relative energy mismatch",
                           "Relative weak divergence"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_loss_eigenmode_ams.json",
                                  "loss_eigenmode_ams", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto eigenvalues =
        LoadRegressionOutputTable("singular_wedge", "loss_eigenmode_ams", "eig.csv");
    auto diagnostics = LoadRegressionOutputTable("singular_wedge", "loss_eigenmode_ams",
                                                 "singular-eigenmode.csv");
    const double backward_error = GetTableValue(eigenvalues, "Error (Bkwd.)");
    const double absolute_error = GetTableValue(eigenvalues, "Error (Abs.)");
    const double mismatch = GetTableValue(diagnostics, "Relative energy mismatch");
    const double divergence = GetTableValue(diagnostics, "Relative weak divergence");
    CAPTURE(backward_error, absolute_error, mismatch, divergence);
    CHECK(backward_error < 1.0e-9);
    CHECK(absolute_error < 1.0e-8);
    CHECK(mismatch < 5.0e-4);
    CHECK(divergence < 1.0e-10);
  }
}

TEST_CASE("singular_wedge_loss_eigenmode_agreement", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error ("};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_loss_compare_eigenmode.json",
                                  "loss_compare_eigenmode", opts);
  auto standard_opts = opts;
  standard_opts.excluded_columns.push_back("Minimum");
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_loss_standard_eigenmode.json",
                                  "loss_standard_eigenmode", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_wedge", "loss_compare_eigenmode", "eig.csv");
    auto standard =
        LoadRegressionOutputTable("singular_wedge", "loss_standard_eigenmode", "eig.csv");
    auto singular_energy = LoadRegressionOutputTable(
        "singular_wedge", "loss_compare_eigenmode", "singular-eigenmode.csv");
    auto standard_energy = LoadRegressionOutputTable(
        "singular_wedge", "loss_standard_eigenmode", "domain-E.csv");

    const double singular_frequency_real = GetTableValue(singular, "Re{f} (GHz)");
    const double standard_frequency_real = GetTableValue(standard, "Re{f} (GHz)");
    const double singular_frequency_imag = GetTableValue(singular, "Im{f} (GHz)");
    const double standard_frequency_imag = GetTableValue(standard, "Im{f} (GHz)");
    const double singular_q = GetTableValue(singular, "Q");
    const double standard_q = GetTableValue(standard, "Q");
    CAPTURE(singular_frequency_real, standard_frequency_real, singular_frequency_imag,
            standard_frequency_imag, singular_q, standard_q);
    CHECK(RelativeDifference(singular_frequency_real, standard_frequency_real) < 5.0e-4);
    CHECK(RelativeDifference(singular_frequency_imag, standard_frequency_imag) < 5.0e-4);
    CHECK(RelativeDifference(singular_q, standard_q) < 1.0e-8);

    const double singular_electric =
        GetTableValue(singular_energy, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_energy, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard_energy, "E_mag (J)");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 1.0e-8);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 1.0e-8);
  }
}

TEST_CASE("singular_wedge_driven", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_driven.json", "driven",
                                  opts);
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_standard_driven.json",
                                  "standard_driven", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_wedge", "driven", "singular-driven.csv");
    auto standard =
        LoadRegressionOutputTable("singular_wedge", "standard_driven", "domain-E.csv");
    auto singular_surface =
        LoadRegressionOutputTable("singular_wedge", "driven", "surface-Q.csv");
    auto standard_surface =
        LoadRegressionOutputTable("singular_wedge", "standard_driven", "surface-Q.csv");

    const double singular_electric = GetTableValue(singular, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard, "E_elec (J)");
    const double singular_magnetic = GetTableValue(singular, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard, "E_mag (J)");
    const double singular_residual = GetTableValue(singular, "Relative residual");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic,
            singular_residual);
    CHECK(singular_residual < 1.0e-7);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 5.0e-3);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 1.5e-2);
    CHECK(RelativeDifference(singular_electric + singular_magnetic,
                             standard_electric + standard_magnetic) < 5.0e-3);

    const double singular_participation = GetTableValue(singular_surface, "p_surf[1]");
    const double standard_participation = GetTableValue(standard_surface, "p_surf[1]");
    CAPTURE(singular_participation, standard_participation);
    CHECK(singular_participation / standard_participation > 1.5);
  }
}

TEST_CASE("singular_wedge_driven_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-14;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  if (palace::Mpi::Size(palace::Mpi::World()) > 1)
  {
    opts.reference_subdir = "amr_driven_mpi";
  }
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_driven_amr.json",
                                  "amr_driven", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    for (const std::string prefix : {"", "iteration1/"})
    {
      auto diagnostics = LoadRegressionOutputTable("singular_wedge", "amr_driven",
                                                   prefix + "singular-driven.csv");
      const double electric = GetTableValue(diagnostics, "Electric field energy (J)");
      const double magnetic = GetTableValue(diagnostics, "Magnetic field energy (J)");
      const double residual = GetTableValue(diagnostics, "Relative residual");
      CAPTURE(prefix, electric, magnetic, residual);
      CHECK(electric > 0.0);
      CHECK(magnetic > 0.0);
      CHECK(residual < 1.0e-8);
    }
  }
}

TEST_CASE("singular_wedge_loss_driven", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_loss_driven.json",
                                  "loss_driven", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_wedge", "loss_driven", "singular-driven.csv");
    const double singular_residual = GetTableValue(singular, "Relative residual");
    CAPTURE(singular_residual);
    CHECK(singular_residual < 1.0e-7);
  }
}

TEST_CASE("singular_wedge_impedance_driven", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_impedance_driven.json",
                                  "impedance_driven", opts);
  auto standard_opts = opts;
  standard_opts.excluded_columns.push_back("Minimum");
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_impedance_standard_driven.json",
                                  "impedance_standard_driven", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular = LoadRegressionOutputTable("singular_wedge", "impedance_driven",
                                              "singular-driven.csv");
    auto standard = LoadRegressionOutputTable("singular_wedge", "impedance_standard_driven",
                                              "domain-E.csv");
    const double singular_electric = GetTableValue(singular, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard, "E_elec (J)");
    const double singular_magnetic = GetTableValue(singular, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard, "E_mag (J)");
    const double residual = GetTableValue(singular, "Relative residual");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic,
            residual);
    CHECK(residual < 1.0e-7);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 2.0e-2);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 1.0e-2);
    CHECK(RelativeDifference(singular_electric + singular_magnetic,
                             standard_electric + standard_magnetic) < 1.0e-2);
  }
}

TEST_CASE("singular_wedge_impedance_driven_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-14;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_impedance_driven_amr.json",
                                  "impedance_amr_driven", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    for (const std::string prefix : {"", "iteration1/"})
    {
      auto diagnostics = LoadRegressionOutputTable("singular_wedge", "impedance_amr_driven",
                                                   prefix + "singular-driven.csv");
      const double electric = GetTableValue(diagnostics, "Electric field energy (J)");
      const double magnetic = GetTableValue(diagnostics, "Magnetic field energy (J)");
      const double residual = GetTableValue(diagnostics, "Relative residual");
      CAPTURE(prefix, electric, magnetic, residual);
      CHECK(electric > 0.0);
      CHECK(magnetic > 0.0);
      CHECK(residual < 1.0e-7);
    }
  }
}

TEST_CASE("singular_wedge_impedance_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error (", "Relative energy mismatch",
                           "Relative weak divergence"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_impedance_eigenmode.json",
                                  "impedance_eigenmode", opts);
  auto standard_opts = opts;
  standard_opts.excluded_columns.insert(standard_opts.excluded_columns.end(),
                                        {"Minimum", "Im{f}", "Q"});
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_impedance_standard_eigenmode.json",
                                  "impedance_standard_eigenmode", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular_mode =
        LoadRegressionOutputTable("singular_wedge", "impedance_eigenmode", "eig.csv");
    auto standard_mode = LoadRegressionOutputTable(
        "singular_wedge", "impedance_standard_eigenmode", "eig.csv");
    auto singular_diagnostics = LoadRegressionOutputTable(
        "singular_wedge", "impedance_eigenmode", "singular-eigenmode.csv");
    auto standard_energy = LoadRegressionOutputTable(
        "singular_wedge", "impedance_standard_eigenmode", "domain-E.csv");

    const double singular_frequency = GetTableValue(singular_mode, "Re{f} (GHz)");
    const double standard_frequency = GetTableValue(standard_mode, "Re{f} (GHz)");
    const double backward_error = GetTableValue(singular_mode, "Error (Bkwd.)");
    const double absolute_error = GetTableValue(singular_mode, "Error (Abs.)");
    const double energy_mismatch =
        GetTableValue(singular_diagnostics, "Relative energy mismatch");
    CAPTURE(singular_frequency, standard_frequency, backward_error, absolute_error,
            energy_mismatch);
    CHECK(RelativeDifference(singular_frequency, standard_frequency) < 1.0e-3);
    CHECK(backward_error < 1.0e-10);
    CHECK(absolute_error < 1.0e-8);
    CHECK(energy_mismatch < 1.0e-9);

    const double singular_electric =
        GetTableValue(singular_diagnostics, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_diagnostics, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard_energy, "E_mag (J)");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 1.0e-3);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 1.0e-3);
  }
}

TEST_CASE("singular_wedge_lumped_driven", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_lumped_driven.json",
                                  "lumped_driven", opts);
  auto standard_opts = opts;
  standard_opts.excluded_columns.push_back("Minimum");
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_standard_driven.json",
                                  "lumped_standard_driven", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_wedge", "lumped_driven", "singular-driven.csv");
    auto standard = LoadRegressionOutputTable("singular_wedge", "lumped_standard_driven",
                                              "domain-E.csv");
    auto standard_voltage =
        LoadRegressionOutputTable("singular_wedge", "lumped_standard_driven", "port-V.csv");
    auto standard_current =
        LoadRegressionOutputTable("singular_wedge", "lumped_standard_driven", "port-I.csv");
    auto standard_scattering =
        LoadRegressionOutputTable("singular_wedge", "lumped_standard_driven", "port-S.csv");
    auto singular_surface =
        LoadRegressionOutputTable("singular_wedge", "lumped_driven", "surface-Q.csv");
    auto standard_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_driven", "surface-Q.csv");

    const double singular_electric = GetTableValue(singular, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard, "E_elec (J)");
    const double singular_magnetic = GetTableValue(singular, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard, "E_mag (J)");
    const double singular_residual = GetTableValue(singular, "Relative residual");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic,
            singular_residual);
    CHECK(singular_residual < 1.0e-7);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 3.0e-3);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 5.0e-4);
    CHECK(RelativeDifference(singular_electric + singular_magnetic,
                             standard_electric + standard_magnetic) < 5.0e-4);

    const double singular_voltage_real = GetTableValue(singular, "Re{V[1]} (V)");
    const double singular_voltage_imag = GetTableValue(singular, "Im{V[1]} (V)");
    const double standard_voltage_real = GetTableValue(standard_voltage, "Re{V[1]} (V)");
    const double standard_voltage_imag = GetTableValue(standard_voltage, "Im{V[1]} (V)");
    CAPTURE(singular_voltage_real, singular_voltage_imag, standard_voltage_real,
            standard_voltage_imag);
    CHECK(RelativeComplexDifference(singular_voltage_real, singular_voltage_imag,
                                    standard_voltage_real, standard_voltage_imag) < 2.0e-4);

    const double singular_current_real = GetTableValue(singular, "Re{I[1]} (A)");
    const double singular_current_imag = GetTableValue(singular, "Im{I[1]} (A)");
    const double standard_current_real = GetTableValue(standard_current, "Re{I[1]} (A)");
    const double standard_current_imag = GetTableValue(standard_current, "Im{I[1]} (A)");
    CAPTURE(singular_current_real, singular_current_imag, standard_current_real,
            standard_current_imag);
    CHECK(RelativeComplexDifference(singular_current_real, singular_current_imag,
                                    standard_current_real, standard_current_imag) < 2.0e-4);

    const double standard_scattering_db =
        GetTableValue(standard_scattering, "|S[1][1]| (dB)");
    const double standard_scattering_phase =
        GetTableValue(standard_scattering, "arg(S[1][1]) (deg.)") * M_PI / 180.0;
    const double standard_scattering_magnitude =
        std::pow(10.0, standard_scattering_db / 20.0);
    const double standard_scattering_real =
        standard_scattering_magnitude * std::cos(standard_scattering_phase);
    const double standard_scattering_imag =
        standard_scattering_magnitude * std::sin(standard_scattering_phase);
    const double singular_scattering_real = GetTableValue(singular, "Re{S[1]}");
    const double singular_scattering_imag = GetTableValue(singular, "Im{S[1]}");
    CAPTURE(singular_scattering_real, singular_scattering_imag, standard_scattering_real,
            standard_scattering_imag);
    CHECK(RelativeComplexDifference(singular_scattering_real, singular_scattering_imag,
                                    standard_scattering_real,
                                    standard_scattering_imag) < 2.0e-4);

    const double singular_participation = GetTableValue(singular_surface, "p_surf[1]");
    const double standard_participation = GetTableValue(standard_surface, "p_surf[1]");
    CAPTURE(singular_participation, standard_participation);
    CHECK(singular_participation / standard_participation > 2.0);
  }
}

TEST_CASE("singular_wedge_lumped_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-7;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Error (",
                           "Relative energy mismatch",
                           "Relative weak divergence",
                           "Im{V[1]}",
                           "Re{I[1]}",
                           "p_ind[1]"};
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("singular_wedge", "singular_wedge_lumped_eigenmode.json",
                                  "lumped_eigenmode", opts);

  auto standard_opts = opts;
  standard_opts.excluded_columns.push_back("Minimum");
  standard_opts.excluded_columns.push_back("Q");
  standard_opts.excluded_columns.push_back("p[1]");
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_standard_eigenmode.json",
                                  "lumped_standard_eigenmode", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular_mode =
        LoadRegressionOutputTable("singular_wedge", "lumped_eigenmode", "eig.csv");
    auto standard_mode =
        LoadRegressionOutputTable("singular_wedge", "lumped_standard_eigenmode", "eig.csv");
    auto singular_diagnostics = LoadRegressionOutputTable(
        "singular_wedge", "lumped_eigenmode", "singular-eigenmode.csv");
    auto standard_energy = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_eigenmode", "domain-E.csv");
    auto standard_voltage = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_eigenmode", "port-V.csv");
    auto standard_current = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_eigenmode", "port-I.csv");
    auto standard_epr = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_eigenmode", "port-EPR.csv");
    auto singular_surface =
        LoadRegressionOutputTable("singular_wedge", "lumped_eigenmode", "surface-Q.csv");
    auto standard_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_standard_eigenmode", "surface-Q.csv");

    const double singular_frequency = GetTableValue(singular_mode, "Re{f} (GHz)");
    const double standard_frequency = GetTableValue(standard_mode, "Re{f} (GHz)");
    const double singular_frequency_imaginary = GetTableValue(singular_mode, "Im{f} (GHz)");
    const double standard_frequency_imaginary = GetTableValue(standard_mode, "Im{f} (GHz)");
    const double singular_backward_error = GetTableValue(singular_mode, "Error (Bkwd.)");
    const double singular_absolute_error = GetTableValue(singular_mode, "Error (Abs.)");
    const double singular_energy_mismatch =
        GetTableValue(singular_diagnostics, "Relative energy mismatch");
    CAPTURE(singular_frequency, standard_frequency, singular_frequency_imaginary,
            standard_frequency_imaginary, singular_backward_error, singular_absolute_error,
            singular_energy_mismatch);
    CHECK(RelativeDifference(singular_frequency, standard_frequency) < 5.0e-4);
    CHECK(std::abs(singular_frequency_imaginary) < 1.0e-12);
    CHECK(std::abs(standard_frequency_imaginary) < 1.0e-12);
    CHECK(singular_backward_error < 1.0e-10);
    CHECK(singular_absolute_error < 1.0e-9);
    CHECK(singular_energy_mismatch < 1.0e-10);

    const double singular_electric =
        GetTableValue(singular_diagnostics, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_diagnostics, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard_energy, "E_mag (J)");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 2.0e-3);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 2.0e-3);

    const double singular_voltage =
        std::hypot(GetTableValue(singular_diagnostics, "Re{V[1]} (V)"),
                   GetTableValue(singular_diagnostics, "Im{V[1]} (V)"));
    const double standard_voltage_magnitude =
        std::hypot(GetTableValue(standard_voltage, "Re{V[1]} (V)"),
                   GetTableValue(standard_voltage, "Im{V[1]} (V)"));
    const double singular_current =
        std::hypot(GetTableValue(singular_diagnostics, "Re{I[1]} (A)"),
                   GetTableValue(singular_diagnostics, "Im{I[1]} (A)"));
    const double standard_current_magnitude =
        std::hypot(GetTableValue(standard_current, "Re{I[1]} (A)"),
                   GetTableValue(standard_current, "Im{I[1]} (A)"));
    CAPTURE(singular_voltage, standard_voltage_magnitude, singular_current,
            standard_current_magnitude);
    CHECK(RelativeDifference(singular_voltage, standard_voltage_magnitude) < 1.5e-2);
    CHECK(RelativeDifference(singular_current, standard_current_magnitude) < 1.5e-2);

    const double singular_epr = std::abs(GetTableValue(singular_diagnostics, "p_ind[1]"));
    const double standard_epr_magnitude = std::abs(GetTableValue(standard_epr, "p[1]"));
    CAPTURE(singular_epr, standard_epr_magnitude);
    CHECK(RelativeDifference(singular_epr, standard_epr_magnitude) < 3.0e-2);

    const double singular_participation = GetTableValue(singular_surface, "p_surf[1]");
    const double standard_participation = GetTableValue(standard_surface, "p_surf[1]");
    CAPTURE(singular_participation, standard_participation);
    CHECK(RelativeDifference(singular_participation, standard_participation) > 5.0e-2);
  }
}

TEST_CASE("singular_wedge_lumped_resistive_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-6;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Error (",
                           "Relative energy mismatch",
                           "Relative weak divergence",
                           "Re{V[1]}",
                           "Im{V[1]}",
                           "Re{I[1]}",
                           "Im{I[1]}",
                           "p_ind[1]"};
  opts.paraview_fields = false;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_resistive_eigenmode.json",
                                  "lumped_resistive_eigenmode", opts);

  auto standard_opts = opts;
  standard_opts.excluded_columns.push_back("Minimum");
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_resistive_standard_eigenmode.json",
                                  "lumped_resistive_standard_eigenmode", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular_mode = LoadRegressionOutputTable("singular_wedge",
                                                   "lumped_resistive_eigenmode", "eig.csv");
    auto standard_mode = LoadRegressionOutputTable(
        "singular_wedge", "lumped_resistive_standard_eigenmode", "eig.csv");
    auto singular_diagnostics = LoadRegressionOutputTable(
        "singular_wedge", "lumped_resistive_eigenmode", "singular-eigenmode.csv");
    auto standard_energy = LoadRegressionOutputTable(
        "singular_wedge", "lumped_resistive_standard_eigenmode", "domain-E.csv");
    auto singular_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_resistive_eigenmode", "surface-Q.csv");
    auto standard_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_resistive_standard_eigenmode", "surface-Q.csv");

    const double singular_frequency_real = GetTableValue(singular_mode, "Re{f} (GHz)");
    const double singular_frequency_imaginary = GetTableValue(singular_mode, "Im{f} (GHz)");
    const double singular_quality = GetTableValue(singular_mode, "Q");
    const double standard_frequency_real = GetTableValue(standard_mode, "Re{f} (GHz)");
    const double standard_frequency_imaginary = GetTableValue(standard_mode, "Im{f} (GHz)");
    const double standard_quality = GetTableValue(standard_mode, "Q");
    CAPTURE(singular_frequency_real, singular_frequency_imaginary, singular_quality,
            standard_frequency_real, standard_frequency_imaginary, standard_quality);
    CHECK(singular_frequency_real > 0.0);
    CHECK(singular_frequency_imaginary > 0.0);
    CHECK(standard_frequency_real > 0.0);
    CHECK(standard_frequency_imaginary > 0.0);
    CHECK(RelativeDifference(
              singular_quality,
              0.5 * std::hypot(singular_frequency_real, singular_frequency_imaginary) /
                  singular_frequency_imaginary) < 1.0e-10);
    CHECK(RelativeDifference(
              standard_quality,
              0.5 * std::hypot(standard_frequency_real, standard_frequency_imaginary) /
                  standard_frequency_imaginary) < 1.0e-10);
    CHECK(RelativeDifference(singular_frequency_real, standard_frequency_real) < 1.0e-2);
    CHECK(RelativeDifference(singular_frequency_imaginary, standard_frequency_imaginary) <
          2.0e-2);
    CHECK(RelativeDifference(singular_quality, standard_quality) < 2.0e-2);

    const double singular_backward_error = GetTableValue(singular_mode, "Error (Bkwd.)");
    const double singular_absolute_error = GetTableValue(singular_mode, "Error (Abs.)");
    const double standard_backward_error = GetTableValue(standard_mode, "Error (Bkwd.)");
    const double standard_absolute_error = GetTableValue(standard_mode, "Error (Abs.)");
    const double energy_mismatch =
        GetTableValue(singular_diagnostics, "Relative energy mismatch");
    const double weak_divergence =
        GetTableValue(singular_diagnostics, "Relative weak divergence");
    CAPTURE(singular_backward_error, singular_absolute_error, standard_backward_error,
            standard_absolute_error, energy_mismatch, weak_divergence);
    CHECK(singular_backward_error < 1.0e-9);
    CHECK(singular_absolute_error < 1.0e-7);
    CHECK(standard_backward_error < 1.0e-9);
    CHECK(standard_absolute_error < 1.0e-7);
    CHECK(energy_mismatch < 1.0e-9);
    CHECK(std::isfinite(weak_divergence));
    CHECK(weak_divergence < 5.0e-1);

    const double singular_electric =
        GetTableValue(singular_diagnostics, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_diagnostics, "Magnetic field energy (J)");
    const double standard_magnetic = GetTableValue(standard_energy, "E_mag (J)");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, standard_magnetic);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 1.0e-6);
    CHECK(RelativeDifference(singular_magnetic, standard_magnetic) < 1.0e-6);

    const double singular_participation = GetTableValue(singular_surface, "p_surf[1]");
    const double standard_participation = GetTableValue(standard_surface, "p_surf[1]");
    const double participation_ratio = singular_participation / standard_participation;
    CAPTURE(singular_participation, standard_participation, participation_ratio);
    CHECK(participation_ratio > 1.5);
    CHECK(participation_ratio < 2.5);
  }
}

TEST_CASE("singular_wedge_lumped_rlc_divfree_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-6;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Error (",
                           "Relative energy mismatch",
                           "Relative weak divergence",
                           "Re{V[1]}",
                           "Im{V[1]}",
                           "Re{I[1]}",
                           "Im{I[1]}",
                           "p_ind[1]"};
  opts.paraview_fields = false;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_rlc_divfree_eigenmode.json",
                                  "lumped_rlc_divfree_eigenmode", opts);
  palace::test::RunRegressionCase("singular_wedge",
                                  "singular_wedge_lumped_rlc_no_divfree_eigenmode.json",
                                  "lumped_rlc_no_divfree_eigenmode", opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto projected_mode = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_divfree_eigenmode", "eig.csv");
    auto unprojected_mode = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_no_divfree_eigenmode", "eig.csv");
    auto projected_diagnostics = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_divfree_eigenmode", "singular-eigenmode.csv");
    auto unprojected_diagnostics = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_no_divfree_eigenmode", "singular-eigenmode.csv");
    auto projected_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_divfree_eigenmode", "surface-Q.csv");
    auto unprojected_surface = LoadRegressionOutputTable(
        "singular_wedge", "lumped_rlc_no_divfree_eigenmode", "surface-Q.csv");

    for (const auto header : {"Re{f} (GHz)", "Im{f} (GHz)", "Q"})
    {
      const double projected = GetTableValue(projected_mode, header);
      const double unprojected = GetTableValue(unprojected_mode, header);
      CAPTURE(header, projected, unprojected);
      CHECK(RelativeDifference(projected, unprojected) < 2.0e-6);
    }
    for (auto *mode : {&projected_mode, &unprojected_mode})
    {
      CHECK(GetTableValue(*mode, "Error (Bkwd.)") < 1.0e-8);
      CHECK(GetTableValue(*mode, "Error (Abs.)") < 1.0e-6);
    }
    for (auto *diagnostics : {&projected_diagnostics, &unprojected_diagnostics})
    {
      CHECK(GetTableValue(*diagnostics, "Relative energy mismatch") < 1.0e-8);
      CHECK(std::isfinite(GetTableValue(*diagnostics, "Relative weak divergence")));
    }
    for (const auto header : {"Electric field energy (J)", "Magnetic field energy (J)",
                              "Relative weak divergence"})
    {
      const double projected = GetTableValue(projected_diagnostics, header);
      const double unprojected = GetTableValue(unprojected_diagnostics, header);
      CAPTURE(header, projected, unprojected);
      CHECK(RelativeDifference(projected, unprojected) < 2.0e-6);
    }

    const double projected_voltage =
        std::hypot(GetTableValue(projected_diagnostics, "Re{V[1]} (V)"),
                   GetTableValue(projected_diagnostics, "Im{V[1]} (V)"));
    const double unprojected_voltage =
        std::hypot(GetTableValue(unprojected_diagnostics, "Re{V[1]} (V)"),
                   GetTableValue(unprojected_diagnostics, "Im{V[1]} (V)"));
    const double projected_current =
        std::hypot(GetTableValue(projected_diagnostics, "Re{I[1]} (A)"),
                   GetTableValue(projected_diagnostics, "Im{I[1]} (A)"));
    const double unprojected_current =
        std::hypot(GetTableValue(unprojected_diagnostics, "Re{I[1]} (A)"),
                   GetTableValue(unprojected_diagnostics, "Im{I[1]} (A)"));
    const double projected_participation = GetTableValue(projected_surface, "p_surf[1]");
    const double unprojected_participation =
        GetTableValue(unprojected_surface, "p_surf[1]");
    CAPTURE(projected_voltage, unprojected_voltage, projected_current, unprojected_current,
            projected_participation, unprojected_participation);
    CHECK(RelativeDifference(projected_voltage, unprojected_voltage) < 2.0e-6);
    CHECK(RelativeDifference(projected_current, unprojected_current) < 2.0e-6);
    CHECK(RelativeDifference(projected_participation, unprojected_participation) < 2.0e-6);
  }
}

TEST_CASE("spheres", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.gridfunction_fields = true;
  palace::test::RunRegressionCase("spheres", "spheres.json", "", opts);
}

TEST_CASE("rings", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("rings", "rings.json", "", opts);
}

TEST_CASE("circular_hole_flux_loop", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("circular_hole", "circular_hole.json", "", opts);
}

// --- cylinder: reltol=1e-4, abstol=1e-16 for the four eigen-style cases ---

TEST_CASE("cylinder_cavity_pec", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cylinder", "cavity_pec.json", "cavity_pec", opts);
}

TEST_CASE("cylinder_cavity_impedance", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cylinder", "cavity_impedance.json", "cavity_impedance",
                                  opts);
}

TEST_CASE("cylinder_waveguide", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cylinder", "waveguide.json", "waveguide", opts);
}

TEST_CASE("cylinder_floquet", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cylinder", "floquet.json", "floquet", opts);
}

// Looser rtol than the cylinder eigen cases above (partition-sensitive).
TEST_CASE("cylinder_driven_wave", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  palace::test::RunRegressionCase("cylinder", "driven_wave.json", "driven_wave", opts);
}

// Floquet-port dielectric grating: structure + Floquet S-parameter magnitudes
// (phase, evanescent-NaN and negligible entries skipped). Tolerances from the
// driven_wave block.
TEST_CASE("dielectric_grating_uniform", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.custom_checks["port-floquet-S.csv"] = TestFloquetSParams(opts.rtol, opts.atol);
  palace::test::RunRegressionCase("dielectric_grating", "dielectric_grating_uniform.json",
                                  "uniform", opts);
}

// Mixed Floquet + lumped/wave port S-parameter tests: validate the sqrt(2) power
// normalization bridge for cross-type observations. Use custom checks to verify power
// balance (sum |S|^2 = 1 for the lossless structure). Floquet NaN and negligible
// entries handled as in dielectric_grating.
TEST_CASE("floquet_lumped", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.paraview_fields = false;
  opts.custom_checks["port-floquet-S.csv"] = TestFloquetSParams(opts.rtol, opts.atol);
  palace::test::RunRegressionCase("floquet_lumped", "floquet_lumped.json", "floquet_lumped",
                                  opts);
}

TEST_CASE("floquet_wave", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.paraview_fields = false;
  opts.custom_checks["port-floquet-S.csv"] = TestFloquetSParams(opts.rtol, opts.atol);
  palace::test::RunRegressionCase("floquet_wave", "floquet_wave.json", "floquet_wave",
                                  opts);
}

// --- antenna: reltol=2e-2, atol=50*1e-10 = 5e-9 for all three ---

TEST_CASE("antenna_halfwave_dipole", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 5.0e-9;
  palace::test::RunRegressionCase("antenna", "antenna_halfwave_dipole.json",
                                  "antenna_halfwave_dipole", opts);
}

TEST_CASE("antenna_short_dipole", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 5.0e-9;
  opts.custom_checks["farfield-rE.csv"] = TestFarfield(opts.rtol);
  palace::test::RunRegressionCase("antenna", "antenna_short_dipole.json",
                                  "antenna_short_dipole", opts);
}

TEST_CASE("antenna_halfwave_dipole_surfacecurrent", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 5.0e-9;
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("antenna", "antenna_halfwave_dipole_surfacecurrent.json",
                                  "antenna_halfwave_dipole_surfacecurrent", opts);
}

// --- coaxial + cpw non-adaptive + cpw eigen: reltol still 2e-2 (from  ---
// --- antenna block), abstol drops to 1e-11 before the coaxial block.  ---

TEST_CASE("coaxial_open", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("coaxial", "coaxial_open.json", "open", opts);
}

TEST_CASE("coaxial_matched", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("coaxial", "coaxial_matched.json", "matched", opts);
}

// Mixed lumped + wave port S-parameters (added in PR #743). Same tolerances as the
// other coaxial cases; port-Z.csv / port-S.csv are diffed like any other CSV.
TEST_CASE("coaxial_lumped_wave", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("coaxial", "coaxial_lumped_wave.json", "lumped_wave",
                                  opts);
}

TEST_CASE("rational_impedance", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("rational_impedance", "parallel_rlc_rational.json",
                                  "rational", opts);
}

// Eigenmode with a series-RLC rational impedance.
TEST_CASE("rational_impedance_eigen", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("rational_impedance", "series_rlc18_rational_eigen.json",
                                  "series_rlc18_eigen", opts);
}

TEST_CASE("cpw_lumped_uniform", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.custom_checks["farfield-rE.csv"] = TestFarfield(opts.rtol);
  palace::test::RunRegressionCase("cpw", "cpw_lumped_uniform.json", "lumped_uniform", opts);
}

TEST_CASE("cpw_wave_uniform", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("cpw", "cpw_wave_uniform.json", "wave_uniform", opts);
}

// Adaptive frequency sweeps run as smoke tests: structure (output tree,
// CSV set, headers, row/column counts) is enforced, numeric comparison
// is disabled via infinite tolerances.
TEST_CASE("cpw_lumped_adaptive", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = std::numeric_limits<double>::infinity();
  opts.atol = std::numeric_limits<double>::infinity();
  palace::test::RunRegressionCase("cpw", "cpw_lumped_adaptive.json", "lumped_adaptive",
                                  opts);
}

// cpw_wave_adaptive omits [GPU] (awslabs/palace#375).
TEST_CASE("cpw_wave_adaptive", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = std::numeric_limits<double>::infinity();
  opts.atol = std::numeric_limits<double>::infinity();
  palace::test::RunRegressionCase("cpw", "cpw_wave_adaptive.json", "wave_adaptive", opts);
}

TEST_CASE("cpw_lumped_eigen", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  // Drop phase-sensitive port columns and Q_ext / κ_ext normalisations;
  // keep the farfield magnitude check.
  opts.excluded_columns = {"Maximum",      "Minimum", "Mean",      "Error (Bkwd.)",
                           "Error (Abs.)", "Re{V[",   "Im{V[",     "Re{I[",
                           "Im{I[",        "Q_ext",   "\u03ba_ext"};
  opts.skip_rowcount = true;
  opts.custom_checks["farfield-rE.csv"] = TestFarfield(opts.rtol);
  palace::test::RunRegressionCase("cpw", "cpw_lumped_eigen.json", "lumped_eigen", opts);
}

TEST_CASE("cpw_wave_eigen", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cpw", "cpw_wave_eigen.json", "wave_eigen", opts);
}

// --- adapter: same reltol / atol as the cpw eigen block above ---

TEST_CASE("adapter_hybrid", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("adapter", "hybrid.json", "hybrid", opts);
}

// Adaptive driven sweep + circuit synthesis with wave ports whose kn(ω) dispersion
// forces the Augmented (rational + auxiliary state) synthesis regime. The raw
// synthesized circuit matrices (rom-Linv/Rinv/C and per-port loads) are excluded from
// comparison: their dimension follows the adaptive sample count and their entries
// depend on the orthogonalized basis, so runs that represent the same physical system
// can differ freely. The physically meaningful regression signal is
// rom-eigenvalues.csv — the eigenvalues of the synthesized L⁻¹/R⁻¹/C system inside the
// trained band (the adapter resonance at 10.526 GHz, Q 30.5, matches the eigenmode
// reference of adapter_hybrid) — together with the S-parameters and the port
// reference table. Omits [GPU] like cpw_wave_adaptive (wave-port adaptive sweeps are
// skipped on GPU, awslabs/palace#375).
TEST_CASE("adapter_driven_synth", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  // The synthesized-eigenvalue rows track the adaptive greedy sampling: modes are
  // only reproducible where the PROM converged, so allow row-count drift and compare
  // the leading (sorted, in-band) modes. Require at least 5 rows (reference has 7
  // modes) so a partial loss of trailing modes cannot pass as a leading-subset match.
  // The HDM eigenpair error columns in rom-eigenvalues.csv are residual diagnostics,
  // not regression targets (same exclusion as the eigenmode cases).
  opts.skip_rowcount = true;
  opts.min_rows = 5;
  opts.excluded_columns = {"Error (Bkwd.)", "Error (Abs.)"};
  // The raw synthesis matrices and the eigenvectors (which live in the same
  // basis-dependent node coordinates) vary with the greedy sampling and MPI partition;
  // only their presence is checked. The eigenvalues are partition-independent and are
  // compared numerically.
  opts.excluded_files = {"rom-Linv",
                         "rom-Rinv",
                         "rom-C-",
                         "rom-portload-",
                         "rom-orthogonalization-matrix-R",
                         "rom-eigenvectors"};
  // No field output is requested in the config.
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("adapter", "driven_synth.json", "driven_synth", opts);
}

// --- transmon: heavy eigen cases. Tagged `[Regression][Long]` so the
// default `ctest -L "^regression$"` sweep skips them (the regression
// block uses TEST_SPEC `[Regression]~[Long]`); CI exercises them via
// the `trigger-long-tests` PR label, which runs `ctest -L "^long$"`.
// Each case takes ~10 minutes at np=2 — hence the dedicated `[Long]`
// modifier and longer per-test TIMEOUT in the CMakeLists registration.

TEST_CASE("transmon_coarse", "[Serial][Parallel][GPU][Regression][Long]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-2;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean",  "Error (Bkwd.)", "Error (Abs.)",
                           "Re{V[",   "Im{V[",   "Re{I[", "Im{I["};
  opts.abs_columns = {"\u03ba_ext"};
  opts.skip_rowcount = true;
  opts.gridfunction_fields = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("transmon", "transmon_coarse.json", "transmon_coarse",
                                  opts);
}

TEST_CASE("transmon_amr", "[Serial][Parallel][GPU][Regression][Long]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-2;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean",  "Error (Bkwd.)", "Error (Abs.)",
                           "Re{V[",   "Im{V[",   "Re{I[", "Im{I["};
  opts.abs_columns = {"\u03ba_ext"};
  opts.skip_rowcount = true;
  opts.gridfunction_fields = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("transmon", "transmon_amr.json", "transmon_amr", opts);
}

// adapter_slp is a dormant entry in runtests.jl: the case is listed in
// the conditional block but neither examples/adapter/slp.json nor
// No adapter/slp reference data exists on this branch. Omit until the
// fixtures are added.

// ===========================================================================
// 2D cases.
// ===========================================================================

TEST_CASE("singular_line_electrostatic", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-9;
  opts.atol = 1.0e-20;
  opts.excluded_columns = {"true_dof"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-coefficients.csv"] =
      CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
  opts.custom_checks["singular-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_line_electrostatic",
                                  "singular_line_electrostatic.json", "", opts);
}

TEST_CASE("singular_line_electrostatic_multigrid", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-9;
  opts.atol = 1.0e-20;
  opts.excluded_columns = {"true_dof"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-coefficients.csv"] =
      CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
  opts.custom_checks["singular-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_line_electrostatic",
                                  "singular_line_electrostatic_p2.json", "p2", opts);
  palace::test::RunRegressionCase("singular_line_electrostatic",
                                  "singular_line_electrostatic_mg2.json", "p2", opts);
}

TEST_CASE("singular_corner_electrostatic", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-8;
  opts.atol = 1.0e-20;
  opts.excluded_columns = {"true_dof"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-coefficients.csv"] =
      CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
  opts.custom_checks["singular-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_corner_electrostatic",
                                  "singular_corner_electrostatic.json", "baseline", opts);
}

TEST_CASE("singular_corner_electrostatic_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-18;
  opts.excluded_columns = {"true_dof"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  for (const std::string prefix : {"", "iteration1/"})
  {
    opts.custom_checks[prefix + "singular-coefficients.csv"] =
        CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
    opts.custom_checks[prefix + "singular-tip-slopes.csv"] =
        CompareCanonicalRows(8, {"source_element"}, opts.rtol, opts.atol);
  }
  palace::test::RunRegressionCase("singular_corner_electrostatic",
                                  "singular_corner_electrostatic_amr.json", "amr", opts);
}

TEST_CASE("singular_corner_electrostatic_amr_nonconforming",
          "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-18;
  opts.excluded_columns = {"true_dof"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  for (const std::string prefix : {"", "iteration1/", "iteration2/"})
  {
    opts.custom_checks[prefix + "singular-coefficients.csv"] =
        CompareCanonicalRows(2, {"true_dof"}, opts.rtol, opts.atol);
    opts.custom_checks[prefix + "singular-tip-slopes.csv"] =
        CompareCanonicalRows(8, {"source_element"}, opts.rtol, opts.atol);
  }
  palace::test::RunRegressionCase("singular_corner_electrostatic",
                                  "singular_corner_electrostatic_amr_nc.json", "amr_nc",
                                  opts);
}

TEST_CASE("standard_corner_electrostatic_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-18;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_corner_electrostatic",
                                  "singular_corner_electrostatic_standard_amr.json",
                                  "standard_amr", opts);
}

TEST_CASE("singular_line_boundarymode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-4;
  opts.atol = 2.0e-9;
  opts.excluded_columns = {"enrichment_true_dof", "Im{kn}", "Im{n_eff}", "Error (", "Im{P}",
                           "coefficient norm"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-mode-nd-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 3.0e-4, 2.0e-3, true);
  opts.custom_checks["singular-mode-h1-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 3.0e-4, 2.0e-3, true);
  opts.custom_checks["singular-mode-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_boundarymode.json", "lossless", opts);
  auto standard_opts = opts;
  standard_opts.rtol = 5.0e-7;
  standard_opts.atol = 1.0e-13;
  standard_opts.excluded_columns = {"Minimum", "Error (", "Im{kn}", "Im{n_eff}"};
  standard_opts.custom_checks.clear();
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_standard_boundarymode.json",
                                  "standard_lossless", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular_mode =
        LoadRegressionOutputTable("singular_line_boundarymode", "lossless", "mode-kn.csv");
    auto standard_mode = LoadRegressionOutputTable("singular_line_boundarymode",
                                                   "standard_lossless", "mode-kn.csv");
    auto singular_diagnostics = LoadRegressionOutputTable(
        "singular_line_boundarymode", "lossless", "singular-mode-diagnostics.csv");
    auto standard_energy = LoadRegressionOutputTable("singular_line_boundarymode",
                                                     "standard_lossless", "domain-E.csv");
    auto singular_surface = LoadRegressionOutputTable("singular_line_boundarymode",
                                                      "lossless", "surface-Q.csv");
    auto standard_surface = LoadRegressionOutputTable("singular_line_boundarymode",
                                                      "standard_lossless", "surface-Q.csv");

    const double singular_kn = GetTableValue(singular_mode, "Re{kn} (1/m)");
    const double standard_kn = GetTableValue(standard_mode, "Re{kn} (1/m)");
    const double singular_neff = GetTableValue(singular_mode, "Re{n_eff}");
    const double standard_neff = GetTableValue(standard_mode, "Re{n_eff}");
    CAPTURE(singular_kn, standard_kn, singular_neff, standard_neff);
    CHECK(RelativeDifference(singular_kn, standard_kn) < 5.0e-4);
    CHECK(RelativeDifference(singular_neff, standard_neff) < 5.0e-4);

    const double singular_electric =
        GetTableValue(singular_diagnostics, "Total electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_diagnostics, "Total magnetic field energy (J)");
    const double singular_power = GetTableValue(singular_diagnostics, "Re{P} (normalized)");
    CAPTURE(singular_electric, standard_electric, singular_magnetic, singular_power);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 5.0e-4);
    CHECK(RelativeDifference(singular_electric, singular_magnetic) < 1.0e-10);
    CHECK(std::abs(singular_power - 1.0) < 1.0e-10);

    const double singular_participation = GetTableValue(singular_surface, "p_surf[1]");
    const double standard_participation = GetTableValue(standard_surface, "p_surf[1]");
    CAPTURE(singular_participation, standard_participation);
    CHECK(std::max(singular_participation, standard_participation) /
              std::min(singular_participation, standard_participation) >
          2.0);
  }
}

TEST_CASE("singular_line_boundarymode_loss", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 3.0e-4;
  opts.atol = 2.0e-9;
  opts.excluded_columns = {"enrichment_true_dof", "Error (", "coefficient norm"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-mode-nd-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 3.0e-4, 2.0e-3, true);
  opts.custom_checks["singular-mode-h1-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 3.0e-4, 2.0e-3, true);
  opts.custom_checks["singular-mode-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, opts.rtol, opts.atol);
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_boundarymode_loss.json", "loss", opts);
  auto standard_opts = opts;
  standard_opts.rtol = 5.0e-7;
  standard_opts.atol = 1.0e-13;
  standard_opts.excluded_columns = {"Minimum", "Error ("};
  standard_opts.custom_checks.clear();
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_standard_boundarymode_loss.json",
                                  "standard_loss", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_line_boundarymode", "loss", "mode-kn.csv");
    auto standard = LoadRegressionOutputTable("singular_line_boundarymode", "standard_loss",
                                              "mode-kn.csv");
    auto singular_diagnostics = LoadRegressionOutputTable(
        "singular_line_boundarymode", "loss", "singular-mode-diagnostics.csv");
    auto standard_energy = LoadRegressionOutputTable("singular_line_boundarymode",
                                                     "standard_loss", "domain-E.csv");

    const double singular_kn_real = GetTableValue(singular, "Re{kn} (1/m)");
    const double standard_kn_real = GetTableValue(standard, "Re{kn} (1/m)");
    const double singular_kn_imag = GetTableValue(singular, "Im{kn} (1/m)");
    const double standard_kn_imag = GetTableValue(standard, "Im{kn} (1/m)");
    const double singular_neff_real = GetTableValue(singular, "Re{n_eff}");
    const double standard_neff_real = GetTableValue(standard, "Re{n_eff}");
    const double singular_neff_imag = GetTableValue(singular, "Im{n_eff}");
    const double standard_neff_imag = GetTableValue(standard, "Im{n_eff}");
    CAPTURE(singular_kn_real, standard_kn_real, singular_kn_imag, standard_kn_imag,
            singular_neff_real, standard_neff_real, singular_neff_imag, standard_neff_imag);
    CHECK(RelativeDifference(singular_kn_real, standard_kn_real) < 5.0e-4);
    CHECK(RelativeDifference(singular_kn_imag, standard_kn_imag) < 5.0e-4);
    CHECK(RelativeDifference(singular_neff_real, standard_neff_real) < 5.0e-4);
    CHECK(RelativeDifference(singular_neff_imag, standard_neff_imag) < 5.0e-4);

    const double singular_electric =
        GetTableValue(singular_diagnostics, "Total electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    CAPTURE(singular_electric, standard_electric);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 5.0e-4);
  }
}

TEST_CASE("singular_line_boundarymode_multigrid", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 2.0e-9;
  opts.excluded_columns = {"enrichment_true_dof", "Im{kn}", "Im{n_eff}", "Error (",
                           "coefficient norm"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.custom_checks["singular-mode-nd-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 2.0e-3, 3.0e-1, true);
  opts.custom_checks["singular-mode-h1-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 2.0e-3, 2.0e-2, true);
  opts.custom_checks["singular-mode-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, 2.0e-3, opts.atol);
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_boundarymode_mg2.json", "mg2", opts);
}

TEST_CASE("singular_line_boundarymode_amr", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-5;
  opts.atol = 1.0e-12;
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  for (const std::string prefix : {"", "iteration1/"})
  {
    opts.custom_checks[prefix + "mode-kn.csv"] = CompareBoundaryModeSpectrum(2.0e-5);
    opts.custom_checks[prefix + "singular-mode-diagnostics.csv"] =
        CompareBoundaryModeSubspaceDiagnostics(1.0e-3);
    opts.custom_checks[prefix + "singular-mode-nd-coefficients.csv"] =
        CompareCanonicalKeySets(2, {"enrichment_true_dof"});
    opts.custom_checks[prefix + "singular-mode-h1-coefficients.csv"] =
        CompareCanonicalKeySets(2, {"enrichment_true_dof"});
    opts.custom_checks[prefix + "singular-mode-tip-slopes.csv"] =
        ValidateBoundaryModeTipSlopes();
    opts.custom_checks[prefix + "surface-Q.csv"] =
        CompareBoundaryModeSurfaceSubspace(5.0e-2);
  }
  palace::test::RunRegressionCase("singular_line_boundarymode",
                                  "singular_line_boundarymode_amr.json", "amr", opts);
}

TEST_CASE("singular_line_eigenmode", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-8;
  opts.atol = 1.0e-13;
  opts.excluded_columns = {"Error ("};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("singular_line_fullwave", "singular_line_eigenmode.json",
                                  "eigenmode", opts);
  auto standard_opts = opts;
  standard_opts.excluded_columns.insert(standard_opts.excluded_columns.end(),
                                        {"Q", "Norm", "Minimum", "Maximum", "Mean"});
  palace::test::RunRegressionCase("singular_line_fullwave",
                                  "singular_line_standard_eigenmode.json",
                                  "standard_eigenmode", standard_opts);

  if (palace::Mpi::Root(palace::Mpi::World()))
  {
    auto singular =
        LoadRegressionOutputTable("singular_line_fullwave", "eigenmode", "eig.csv");
    auto standard = LoadRegressionOutputTable("singular_line_fullwave",
                                              "standard_eigenmode", "eig.csv");
    auto singular_energy = LoadRegressionOutputTable("singular_line_fullwave", "eigenmode",
                                                     "singular-eigenmode.csv");
    auto standard_energy = LoadRegressionOutputTable("singular_line_fullwave",
                                                     "standard_eigenmode", "domain-E.csv");

    const double singular_frequency = GetTableValue(singular, "Re{f} (GHz)");
    const double standard_frequency = GetTableValue(standard, "Re{f} (GHz)");
    const double singular_damping = GetTableValue(singular, "Im{f} (GHz)");
    const double standard_damping = GetTableValue(standard, "Im{f} (GHz)");
    CAPTURE(singular_frequency, standard_frequency, singular_damping, standard_damping);
    CHECK(RelativeDifference(singular_frequency, standard_frequency) < 1.0e-2);
    CHECK(std::abs(singular_damping) / singular_frequency < 1.0e-10);
    CHECK(std::abs(standard_damping) / standard_frequency < 1.0e-10);

    const double singular_electric =
        GetTableValue(singular_energy, "Electric field energy (J)");
    const double standard_electric = GetTableValue(standard_energy, "E_elec (J)");
    const double singular_magnetic =
        GetTableValue(singular_energy, "Magnetic field energy (J)");
    const double singular_energy_mismatch =
        GetTableValue(singular_energy, "Relative energy mismatch");
    CAPTURE(singular_electric, standard_electric, singular_magnetic,
            singular_energy_mismatch);
    CHECK(RelativeDifference(singular_electric, standard_electric) < 1.0e-10);
    CHECK(RelativeDifference(singular_electric, singular_magnetic) < 1.0e-10);
    CHECK(singular_energy_mismatch < 1.0e-10);
  }
}

// cavity2d eigenmode uses the tight 1e-4 / 1e-16 tolerance.
TEST_CASE("cavity2d_eigenmode", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d.json", "eigenmode", opts);
}

// Coarse 2D driven tolerances (partition- and platform-sensitive).
TEST_CASE("cavity2d_driven", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-8;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d_driven.json", "driven", opts);
}

TEST_CASE("cavity2d_electrostatic", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = {"Maximum", "Minimum"};
  palace::test::RunRegressionCase("cavity2d", "cavity2d_electrostatic.json",
                                  "electrostatic", opts);
}

TEST_CASE("cavity2d_magnetostatic", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d_magnetostatic.json",
                                  "magnetostatic", opts);
}

TEST_CASE("cavity2d_transient", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d_transient.json", "transient", opts);
}

// cpw2d boundary-mode cases: reltol=1e-2 because impedance depends on
// voltage path integration (cross-platform reproducibility).
TEST_CASE("cpw2d_thin", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-2;
  opts.atol = 1.0e-10;  // abstol from the preceding cavity2d block
  opts.excluded_columns = {"Maximum",      "Minimum",      "Mean",     "Error (Bkwd.)",
                           "Error (Abs.)", "Im{kn} (1/m)", "Im{n_eff}"};
  opts.skip_rowcount = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["mode-V.csv"] = CompareComplexMagnitudes(opts.rtol, opts.atol);
  palace::test::RunRegressionCase("cpw2d", "cpw2d_thin.json", "thin", opts);
}

TEST_CASE("cpw2d_thin_driven_singular", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 5.0e-3;
  opts.atol = 1.0e-15;
  opts.excluded_columns = {"Relative residual"};
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cpw2d", "cpw2d_thin_driven_singular.json",
                                  "thin_driven_singular", opts);
}

TEST_CASE("cpw2d_thick_impedance", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-2;
  opts.atol = 1.0e-10;
  opts.excluded_columns = {"Maximum",      "Minimum",      "Mean",     "Error (Bkwd.)",
                           "Error (Abs.)", "Im{kn} (1/m)", "Im{n_eff}"};
  opts.skip_rowcount = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cpw2d", "cpw2d_thick_impedance.json", "thick_impedance",
                                  opts);
}

TEST_CASE("cpw2d_thick_impedance_singular", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 3.0e-4;
  opts.atol = 2.0e-9;
  opts.excluded_columns = {"enrichment_true_dof", "Im{kn}",
                           "Im{n_eff}",           "Error (",
                           "coefficient norm",    "Im{P} (normalized)"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.custom_checks["singular-mode-nd-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 1.0e-2, 2.0e-3, true);
  opts.custom_checks["singular-mode-h1-coefficients.csv"] =
      CompareCanonicalRows(2, {"enrichment_true_dof"}, 1.0e-2, 2.0e-3, true);
  opts.custom_checks["singular-mode-tip-slopes.csv"] =
      CompareCanonicalRows(8, {}, 1.0e-2, opts.atol);
  palace::test::RunRegressionCase("cpw2d", "cpw2d_thick_impedance_singular.json",
                                  "thick_impedance_singular", opts);
}

// 2D mode analysis from a 3D mesh: back to reltol=1e-4, atol still 1e-10.
TEST_CASE("cpw_wave_2dmode", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = kEigenExcluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("cpw", "cpw_wave_2dmode.json", "wave_2dmode", opts);
}
