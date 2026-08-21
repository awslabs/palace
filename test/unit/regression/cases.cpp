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
#include <map>
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

// Standard "drop per-element extrema + eigenmode error columns" list.
const std::vector<std::string> kEigenExcluded = {"Maximum", "Minimum", "Mean",
                                                 "Error (Bkwd.)", "Error (Abs.)"};
constexpr auto kForceDefaultSolver = palace::test::SolverOverridePolicy::ForceDefault;

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

// Power-conservation check for a driven wave-port run on a lossless structure: for
// each excitation j, require Sum_i |S_ij|^2 == 1 within rtol. Reference-free (asserts
// the invariant directly), so it is insensitive to MPI partitioning. Guards the
// wave-port n x H fix (GetModalCorrectionOperator): without it Sum|S|^2 > 1 for modes
// with E_n != 0 (cylinder TM01, |S11| ~ 1.07). Column headers are "|S[i][j]| (dB)".
palace::test::CustomCheck TestWavePortLossless(double rtol)
{
  return [rtol](palace::Table &actual, palace::Table &)
  {
    const std::size_t n_rows = actual.n_rows();
    for (std::size_t r = 0; r < n_rows; ++r)
    {
      std::map<std::string, double> power_by_excitation;
      for (auto &c : actual)
      {
        const std::string &hdr = c.header_text;
        if (hdr.rfind("|S[", 0) != 0 || hdr.find("(dB)") == std::string::npos)
        {
          continue;
        }
        // Extract the excitation index j (second [..]) from "|S[i][j]| (dB)".
        const auto b1 = hdr.find(']');
        const auto b2 = hdr.find('[', b1);
        const auto b2e = hdr.find(']', b2);
        if (b1 == std::string::npos || b2 == std::string::npos || b2e == std::string::npos)
        {
          continue;
        }
        const std::string j = hdr.substr(b2 + 1, b2e - b2 - 1);
        const double s_lin = std::pow(10.0, c.data[r] / 20.0);  // dB -> linear |S|
        power_by_excitation[j] += s_lin * s_lin;
      }
      for (const auto &[j, sum_sq] : power_by_excitation)
      {
        INFO("row " << r + 1 << ", excitation " << j << ": Sum_i |S[i][" << j
                    << "]|^2 = " << sum_sq);
        CHECK_THAT(sum_sq, Catch::Matchers::WithinAbs(1.0, rtol));
      }
    }
  };
}

}  // namespace

// ===========================================================================
// 3D cases. Order follows the retired Julia regression harness.
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

// Drives the TM01 port mode (Mode 3) of a lossless PEC-shorted circular waveguide.
// TM01 carries a longitudinal field (E_n != 0 on the port face), so a correct wave
// port must still reflect all incident power: |S11| = 1. TestWavePortLossless asserts
// Sum_i |S_i1|^2 = 1 directly. The dominant TE11 mode (Mode 1, E_n = 0) is a separate
// case (cylinder_driven_wave).
TEST_CASE("cylinder_driven_wave_tm", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.custom_checks["port-S.csv"] = TestWavePortLossless(1.0e-3);
  palace::test::RunRegressionCase("cylinder", "driven_wave_tm.json", "driven_wave_tm",
                                  opts);
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
