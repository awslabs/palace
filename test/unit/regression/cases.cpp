// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Catch2 ports of the Julia regression cases in test/examples/runtests.jl.
// Each case is tagged [Serial][Parallel][GPU][Regression] (regression
// is a category, orthogonal to the [Serial]/[Parallel]/[GPU] execution
// styles; the same case is valid at any rank count and on either CPU
// or GPU). The two transmon eigen cases additionally carry [Long] —
// the default regression sweep is registered with TEST_SPEC
// `[Regression]~[Long]` so they run only under the long-tests CI
// workflow.
//
// Per-case rtol / atol / excluded_columns / skip_rowcount track the
// Julia values verbatim, with variable names from the original file
// preserved in comments where they change between cases.

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "regression_helpers.hpp"
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
// at the case's rtol. Column layout in farfield-rE.csv:
//   3=r*Re{Ex}, 4=r*Im{Ex}, 5=r*Re{Ey}, 6=r*Im{Ey}, 7=r*Re{Ez}, 8=r*Im{Ez}
palace::test::CustomCheck TestFarfield(double rtol)
{
  return [rtol](palace::Table &a, palace::Table &r)
  {
    CHECK(a.n_cols() >= 9);
    CHECK(r.n_cols() >= 9);
    if (a.n_cols() < 9 || r.n_cols() < 9)
    {
      return;
    }
    const std::size_t n_rows = std::min(a.n_rows(), r.n_rows());
    for (std::size_t i = 0; i < n_rows; ++i)
    {
      auto mag_sq = [&](palace::Table &t)
      {
        auto sq = [](double x, double y) { return x * x + y * y; };
        return sq(t[3].data[i], t[4].data[i]) + sq(t[5].data[i], t[6].data[i]) +
               sq(t[7].data[i], t[8].data[i]);
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

// Standard "drop per-element extrema + eigenmode error columns" list.
const std::vector<std::string> kEigenExcluded = {"Maximum", "Minimum", "Mean",
                                                 "Error (Bkwd.)", "Error (Abs.)"};
constexpr auto kForceDefaultSolver = palace::test::SolverOverridePolicy::ForceDefault;

}  // namespace

// ===========================================================================
// 3D cases. Order follows test/examples/runtests.jl.
// ===========================================================================

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

// cylinder/driven_wave bumps reltol to 1e-3 in Julia.
TEST_CASE("cylinder_driven_wave", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  palace::test::RunRegressionCase("cylinder", "driven_wave.json", "driven_wave", opts);
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

TEST_CASE("cpw_lumped_uniform", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = {"Maximum", "Minimum"};
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

// Julia kept adaptive frequency sweeps as smoke tests only: enforce the
// output tree, CSV set, headers, and row/column structure, but skip all
// numeric comparisons with infinite tolerances.
TEST_CASE("cpw_lumped_adaptive", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = std::numeric_limits<double>::infinity();
  opts.atol = std::numeric_limits<double>::infinity();
  palace::test::RunRegressionCase("cpw", "cpw_lumped_adaptive.json", "lumped_adaptive",
                                  opts);
}

TEST_CASE("cpw_wave_adaptive", "[Serial][Parallel][GPU][Regression]")
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
  // Match Julia parity: drop eigenmode phase-sensitive port columns and
  // Q_ext / κ_ext normalisation values; keep the farfield magnitude check.
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
  opts.skip_rowcount = true;
  opts.gridfunction_fields = true;
  opts.linear_solver_policy = kForceDefaultSolver;
  opts.eigen_solver_policy = kForceDefaultSolver;
  palace::test::RunRegressionCase("transmon", "transmon_amr.json", "transmon_amr", opts);
}

// adapter_slp is a dormant entry in runtests.jl: the case is listed in
// the conditional block but neither examples/adapter/slp.json nor
// test/examples/ref/adapter/slp/ exists on this branch. Omit until the
// fixtures are added.

// ===========================================================================
// 2D cases. Live on hughcars/boundarymode-2d; not yet on origin/main.
// ===========================================================================

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

// Julia switches back to reltol=1e-4, abstol=1e-10 for the remaining
// cavity2d cases.
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
