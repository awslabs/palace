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
#include <complex>
#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>
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
  return [rtol](palace::Table &a, palace::Table &r, const std::filesystem::path &)
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
  return [rtol, atol](palace::Table &actual, palace::Table &reference,
                      const std::filesystem::path &)
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
const std::vector<std::string> eigen_excluded = {"Maximum", "Minimum", "Mean",
                                                 "Error (Bkwd.)", "Error (Abs.)"};
constexpr auto force_default_solver = palace::test::SolverOverridePolicy::ForceDefault;

// Floquet-port S-parameters: compare only the |S[...]| (dB) magnitude columns
// (phase isn't reproducible). NaN entries (evanescent modes) and signals below
// -200 dB (negligible) count as matches. Ports the Julia test_floquet_sparams.
palace::test::CustomCheck TestFloquetSParams(double rtol, double atol)
{
  return [rtol, atol](palace::Table &actual, palace::Table &reference,
                      const std::filesystem::path &)
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

// S-parameter check for a driven wave-port run on a lossless structure. Two guards:
//   (1) Reference-free power conservation: for each excitation j, Sum_i |S_ij|^2 == 1
//       within rtol. Insensitive to MPI partitioning; catches the wave-port n x H fix
//       (GetModalCorrectionOperator) regressing (without it Sum|S|^2 > 1 for E_n != 0
//       modes, e.g. cylinder TM01 |S11| ~ 1.07).
//   (2) Reference diff of the S-parameter PHASE columns "arg(S[i][j]) (deg.)" within
//       rtol/atol. For a lossless structure power conservation already pins the
//       magnitudes, so the phase is the discriminating quantity (it encodes the port
//       reactance/impedance): a wrong-phase result that still conserves power is caught
//       here. The magnitude column is skipped — |S|=1 is ~0 dB, where a relative diff is
//       meaningless.
palace::test::CustomCheck TestWavePortLossless(double rtol, double atol = 1.0e-16)
{
  return [rtol, atol](palace::Table &actual, palace::Table &reference,
                      const std::filesystem::path &)
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
    // (2) Reference diff of the phase columns only.
    const std::size_t n_cols = std::min(actual.n_cols(), reference.n_cols());
    const std::size_t n_cmp = std::min(actual.n_rows(), reference.n_rows());
    for (std::size_t c = 0; c < n_cols; ++c)
    {
      const palace::Column &ac = actual[c];
      const palace::Column &rc = reference[c];
      if (rc.header_text.rfind("arg(S", 0) != 0)
      {
        continue;
      }
      for (std::size_t r = 0; r < n_cmp; ++r)
      {
        INFO("row " << r + 1 << " column '" << rc.header_text << "'");
        CHECK_THAT(ac.data[r], Catch::Matchers::WithinRel(rc.data[r], rtol) ||
                                   Catch::Matchers::WithinAbs(rc.data[r], atol));
      }
    }
  };
}

// Assert the synthesized model reproduces the eigenmode resonant frequency. Among the
// physical resonances (Q > 1, which excludes the aux-realization's critically-damped
// spurious roots), the one nearest f_re_eigen must match it in Re{f} within rtol. Robust to
// the basis/partition-dependent root count and ordering.
//
// q_eigen > 0 additionally asserts that same root's Q within q_rtol. Pass it when the
// pole is extractable from a real-frequency fit (a moderate-to-high-Q resonance sampled by
// a band that brackets it) and the common-pole AAA fit resolves it as a converged HDM
// eigenpair (backward error ~1e-9), so its Q is toolchain/partition-robust. For a broad
// off-axis pole that stays under-determined by real-axis samples, leave q_eigen < 0 to
// skip.
palace::test::CustomCheck TestRomEigenvalueMatchesEigenmode(double f_re_eigen, double rtol,
                                                            double q_eigen = -1.0,
                                                            double q_rtol = 0.0)
{
  return [f_re_eigen, rtol, q_eigen, q_rtol](palace::Table &actual, palace::Table &,
                                             const std::filesystem::path &)
  {
    int col_re = -1, col_q = -1;
    for (std::size_t c = 0; c < actual.n_cols(); ++c)
    {
      const std::string &hdr = actual[c].header_text;
      if (hdr.find("Re{f}") != std::string::npos)
      {
        col_re = static_cast<int>(c);
      }
      else if (hdr.find("Q") != std::string::npos)
      {
        col_q = static_cast<int>(c);
      }
    }
    REQUIRE(col_re >= 0);
    REQUIRE(col_q >= 0);
    const auto &re = actual[col_re].data;
    const auto &q = actual[col_q].data;
    double best = -1.0, best_re = 0.0, best_q = 0.0;
    for (std::size_t r = 0; r < re.size(); ++r)
    {
      if (q[r] <= 1.0)  // skip critically-damped spurious roots
      {
        continue;
      }
      const double d = std::abs(re[r] - f_re_eigen);
      if (best < 0.0 || d < best)
      {
        best = d;
        best_re = re[r];
        best_q = q[r];
      }
    }
    REQUIRE(best >= 0.0);  // at least one physical resonance present
    INFO("nearest physical (Q>1) synth root Re{f} = " << best_re << " GHz (Q = " << best_q
                                                      << ") vs eigenmode " << f_re_eigen
                                                      << " GHz (Q = " << q_eigen << ")");
    CHECK_THAT(best_re, Catch::Matchers::WithinRel(f_re_eigen, rtol));
    if (q_eigen > 0.0)
    {
      CHECK_THAT(best_q, Catch::Matchers::WithinRel(q_eigen, q_rtol));
    }
  };
}

// End-to-end de-embed/reconnect round-trip for a synthesized wave-port network. The swept
// port-S.csv is field-derived (modal overlap of the ROM-reconstructed field); the
// rom-{Linv,Rinv,C} pencil is the independently synthesized circuit. This check closes the
// loop: it reduces the total pencil Y(iω) = L⁻¹/(iω) + R⁻¹ + iω·C to the physical port
// terminals by a Schur complement onto the internal (basis + aux) states, de-embeds each
// port's own reference load y_ref (baked into the port's diagonal load block), renormalizes
// to the real modal reference y_ref via a Kurokawa power-wave transform, and checks the
// resulting |S| against Palace's field-derived port-S. Agreement validates that the
// exported circuit reproduces the solver's own S-parameters. The W modal correction lives
// in the port loads, so a regressed W shifts |S| by dB-scale amounts. Frequencies and
// per-port y_ref are read from the sibling rom-port-reference.csv (same sweep grid).
palace::test::CustomCheck TestWavePortSRoundTrip(double atol_lin)
{
  return [atol_lin](palace::Table &, palace::Table &reference,
                    const std::filesystem::path &actual_path)
  {
    namespace fs = std::filesystem;
    const fs::path dir = actual_path.parent_path();

    // Read a synthesized matrix CSV. Table layout (see RomOperator::print_table): column c
    // has header = node label, data[r] = M(r, c). A zero real/imag part is not written, so
    // a missing file yields an empty matrix (treated as zero by the caller).
    auto read_real = [](const fs::path &p, std::vector<std::string> *lab) -> Eigen::MatrixXd
    {
      if (!fs::is_regular_file(p))
      {
        return Eigen::MatrixXd(0, 0);
      }
      palace::TableWithCSVFile w(p.string(), /*load_existing_file=*/true);
      palace::Table &t = w.table;
      const long n = static_cast<long>(t.n_cols());
      const long m = static_cast<long>(t.n_rows());
      Eigen::MatrixXd M(m, n);
      for (long c = 0; c < n; ++c)
      {
        if (lab)
        {
          lab->push_back(t[c].header_text);
        }
        for (long r = 0; r < m; ++r)
        {
          M(r, c) = t[c].data[r];
        }
      }
      return M;
    };

    // L⁻¹ real part is always populated (α₀ + HDM K), so it defines the label set and size.
    std::vector<std::string> labels;
    const Eigen::MatrixXd linv_re = read_real(dir / "rom-Linv-re.csv", &labels);
    REQUIRE(linv_re.rows() > 0);
    const long n_total = linv_re.rows();
    auto term = [&](const std::string &stem) -> Eigen::MatrixXcd
    {
      const Eigen::MatrixXd re = read_real(dir / (stem + "-re.csv"), nullptr);
      const Eigen::MatrixXd im = read_real(dir / (stem + "-im.csv"), nullptr);
      Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(n_total, n_total);
      if (re.rows() == n_total && re.cols() == n_total)
      {
        M.real() = re;
      }
      if (im.rows() == n_total && im.cols() == n_total)
      {
        M.imag() = im;
      }
      return M;
    };
    const Eigen::MatrixXcd linv = term("rom-Linv");
    const Eigen::MatrixXcd rinv = term("rom-Rinv");
    const Eigen::MatrixXcd cap = term("rom-C");

    // Port terminals and their real modal reference y_ref(ω), read from the sibling
    // rom-port-reference.csv. Its "Re{Y_ref[<label>]} (S)" columns give the port order and
    // the de-embed/renormalization reference; its "f (GHz)" column gives the sweep grid.
    palace::TableWithCSVFile ref_w((dir / "rom-port-reference.csv").string(),
                                   /*load_existing_file=*/true);
    palace::Table &ref_t = ref_w.table;
    std::vector<std::string> port_labels;
    std::vector<int> re_yref_col;
    int f_col = -1;
    for (std::size_t c = 0; c < ref_t.n_cols(); ++c)
    {
      const std::string &h = ref_t[c].header_text;
      if (h.rfind("f (GHz)", 0) == 0)
      {
        f_col = static_cast<int>(c);
      }
      else if (h.rfind("Re{Y_ref[", 0) == 0)
      {
        const auto lb = h.find('[');
        const auto rb = h.find(']', lb);
        port_labels.push_back(h.substr(lb + 1, rb - lb - 1));
        re_yref_col.push_back(static_cast<int>(c));
      }
    }
    REQUIRE(f_col >= 0);
    const long n_ports = static_cast<long>(port_labels.size());
    REQUIRE(n_ports >= 1);

    // Map each port terminal to its pencil row, and everything else to the internal block.
    std::vector<long> pidx;
    for (const auto &pl : port_labels)
    {
      const auto it = std::find(labels.begin(), labels.end(), pl);
      REQUIRE(it != labels.end());
      pidx.push_back(static_cast<long>(it - labels.begin()));
    }
    std::vector<long> iidx;
    for (long i = 0; i < n_total; ++i)
    {
      if (std::find(pidx.begin(), pidx.end(), i) == pidx.end())
      {
        iidx.push_back(i);
      }
    }
    const long n_int = static_cast<long>(iidx.size());

    // Sweep grid + per-port y_ref from the reference table.
    const std::size_t n_freq = ref_t.n_rows();
    const std::complex<double> j1(0.0, 1.0);
    for (std::size_t rr = 0; rr < n_freq; ++rr)
    {
      const double f_ghz = ref_t[f_col].data[rr];
      const std::complex<double> s = j1 * (2.0 * M_PI * f_ghz * 1.0e9);

      // Full nodal admittance, then Schur-reduce onto the n_ports terminals.
      const Eigen::MatrixXcd y = linv / s + rinv + s * cap;
      Eigen::MatrixXcd ypp(n_ports, n_ports), ypi(n_ports, n_int), yip(n_int, n_ports),
          yii(n_int, n_int);
      for (long a = 0; a < n_ports; ++a)
      {
        for (long b = 0; b < n_ports; ++b)
        {
          ypp(a, b) = y(pidx[a], pidx[b]);
        }
        for (long b = 0; b < n_int; ++b)
        {
          ypi(a, b) = y(pidx[a], iidx[b]);
          yip(b, a) = y(iidx[b], pidx[a]);
        }
      }
      for (long a = 0; a < n_int; ++a)
      {
        for (long b = 0; b < n_int; ++b)
        {
          yii(a, b) = y(iidx[a], iidx[b]);
        }
      }
      Eigen::MatrixXcd yports = ypp;
      if (n_int > 0)
      {
        yports -= ypi * yii.fullPivLu().solve(yip);
      }

      // De-embed each port's own reference load (baked into the diagonal), then renormalize
      // to y_ref via the Kurokawa power-wave transform S = √Y0 (I − Z0·Yd)(I + Z0·Yd)⁻¹ √Z0
      // with Z0 = diag(1/y_ref) (real reference).
      Eigen::VectorXd yref(n_ports);
      for (long a = 0; a < n_ports; ++a)
      {
        yref(a) = ref_t[re_yref_col[a]].data[rr];
        yports(a, a) -= yref(a);
      }
      Eigen::MatrixXcd z0 = Eigen::MatrixXcd::Zero(n_ports, n_ports);
      Eigen::MatrixXcd rootY = Eigen::MatrixXcd::Zero(n_ports, n_ports);
      Eigen::MatrixXcd rootZ = Eigen::MatrixXcd::Zero(n_ports, n_ports);
      for (long a = 0; a < n_ports; ++a)
      {
        z0(a, a) = 1.0 / yref(a);
        rootY(a, a) = std::sqrt(yref(a));
        rootZ(a, a) = 1.0 / std::sqrt(yref(a));
      }
      const Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(n_ports, n_ports);
      const Eigen::MatrixXcd smat =
          rootY * (id - z0 * yports) * (id + z0 * yports).fullPivLu().inverse() * rootZ;

      // Compare |S[i][j]| against the field-derived reference (dB → linear).
      for (long i = 0; i < n_ports; ++i)
      {
        for (long jj = 0; jj < n_ports; ++jj)
        {
          const std::string hdr =
              "|S[" + std::to_string(i + 1) + "][" + std::to_string(jj + 1) + "]| (dB)";
          int sc = -1;
          for (std::size_t c = 0; c < reference.n_cols(); ++c)
          {
            if (reference[c].header_text == hdr)
            {
              sc = static_cast<int>(c);
              break;
            }
          }
          if (sc < 0 || rr >= reference.n_rows())
          {
            continue;
          }
          const double meas_lin = std::pow(10.0, reference[sc].data[rr] / 20.0);
          const double pen_lin = std::abs(smat(i, jj));
          INFO("f = " << f_ghz << " GHz, |S[" << i + 1 << "][" << jj + 1 << "]| pencil "
                      << pen_lin << " vs field " << meas_lin);
          CHECK_THAT(pen_lin, Catch::Matchers::WithinAbs(meas_lin, atol_lin));
        }
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
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("rings", "rings.json", "", opts);
}

TEST_CASE("rings_multiring_inactive_ports", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("multiring", "multiring_inactive_ports.json",
                                  "multiring_inactive_ports", opts);
}

TEST_CASE("circular_hole_flux_loop", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("circular_hole", "circular_hole.json", "", opts);
}

// Mixed current-flux excitation. The aperture integral recovering M[1][2] is
// reduced over surfaces the partitioner may split, so this case catches a
// double-counted contribution.
TEST_CASE("ring_disk_mixed_current_flux", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("ring_disk", "ring_disk.json", "", opts);
}

// --- cylinder: reltol=1e-4, abstol=1e-16 for the four eigen-style cases ---

TEST_CASE("cylinder_cavity_pec", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("cylinder", "cavity_pec.json", "cavity_pec", opts);
}

TEST_CASE("cylinder_cavity_impedance", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cylinder", "cavity_impedance.json", "cavity_impedance",
                                  opts);
}

TEST_CASE("cylinder_waveguide", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cylinder", "waveguide.json", "waveguide", opts);
}

TEST_CASE("cylinder_floquet", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-16;
  opts.excluded_columns = eigen_excluded;
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

// TM01 (E_n != 0) driven through a uniaxial-anisotropic mu (mu_t != mu_z); mu_x = mu_y
// keeps TM01 clean so the anisotropic mu^-1 path must still give Sum_i |S_i1|^2 = 1.
// Rotated mu is covered by the WavePortData reconstruction unit test.
TEST_CASE("cylinder_driven_wave_tm_aniso", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.custom_checks["port-S.csv"] = TestWavePortLossless(1.0e-3);
  palace::test::RunRegressionCase("cylinder", "driven_wave_tm_aniso.json",
                                  "driven_wave_tm_aniso", opts);
}

// TM01 (E_n != 0) swept just above cutoff (~2.9 GHz): small k_n stresses the
// modal-correction conditioning while the port propagates, so Sum_i |S_i1|^2 = 1 must hold
// across the sweep.
TEST_CASE("cylinder_driven_wave_tm_cutoff", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.custom_checks["port-S.csv"] = TestWavePortLossless(1.0e-3);
  palace::test::RunRegressionCase("cylinder", "driven_wave_tm_cutoff.json",
                                  "driven_wave_tm_cutoff", opts);
}

// Adaptive (PROM) counterpart of cylinder_driven_wave_tm: the same non-TEM TM01 mode swept
// 3.1-3.6 GHz, exercising the wave-port modal correction W on the reduced operator.
// TestWavePortLossless asserts Sum_i |S_i1|^2 = 1 (reference-free); partition-dependent
// error-indicators.csv / domain-E.csv are kept for the file-set check but not diffed.
// Omits [GPU] like cpw_wave_adaptive (awslabs/palace#375).
TEST_CASE("cylinder_driven_wave_tm_adaptive", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = {"Maximum", "Minimum", "Mean"};
  opts.excluded_files = {"domain-E.csv", "error-indicators.csv"};
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  opts.custom_checks["port-S.csv"] = TestWavePortLossless(1.0e-3);
  palace::test::RunRegressionCase("cylinder", "driven_wave_tm_adaptive.json",
                                  "driven_wave_tm_adaptive", opts);
}

// Eigenmode counterpart of iris_filter_driven_wave_synth. Both wave ports of this two-port
// iris-coupled slab-loaded guide carry the hybrid/LSM mode (Mode 3), so the modal
// correction W is rank>=2 and enters the nonlinear eigensolver on the LHS. Dropping or
// breaking W shifts the cavity resonance Q (Re{f} ~ 6.690 GHz, Q ~ 22), so the diffed
// eig.csv is the W-sensitivity guard.
TEST_CASE("iris_filter_wave_eigen", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-16;
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("iris_filter", "wave_eigen.json", "wave_eigen", opts);
}

// Circuit synthesis of the two-port iris-coupled slab-loaded guide (hybrid/LSM port, W
// rank>=2). The inductive irises make the cavity a strongly-coupled |S21| resonance, so the
// greedy real-frequency samples capture the pole and the synthesized root reproduces the
// eigenmode Re{f} (6.690 GHz) AND Q (~22.4) across partitions, asserted by
// TestRomEigenvalueMatchesEigenmode. The common-pole AAA + signed-residue fit resolves the
// pole as a converged HDM eigenpair (backward error ~1e-9), so Q is partition-stable
// (np=1,2,4,8 all give Q=22.40); swept S and the pencil matrices are partition-dependent
// and presence-checked only.
TEST_CASE("iris_filter_driven_wave_synth", "[Serial][Parallel][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-3;
  opts.atol = 1.0e-11;
  opts.skip_rowcount = true;
  opts.min_rows = 1;
  // The physical (W-inclusive) Y_ref carries a tiny partition-dependent residue in its
  // imaginary part (Im{Y}~1e-9, Im{Z}~1e-3 Ohm) that the old scalar-kn reference had as
  // exactly zero; the real parts carry the W signal and are partition-stable, so diff those
  // and drop the imaginary columns.
  opts.excluded_columns = {"Error (Bkwd.)", "Error (Abs.)", "Maximum", "Minimum",
                           "Mean",          "Im{Y_ref",     "Im{Z_ref"};
  opts.excluded_files = {"rom-Linv", "rom-Rinv", "rom-C-", "rom-portload-",
                         "rom-orthogonalization-matrix-R", "rom-eigenvectors",
                         // Sharp resonance: swept S and the per-element error-estimator
                         // extrema are partition/arithmetic-sensitive near the pole, so the
                         // W-dependent signal is the synthesized eigenvalue (custom check)
                         // rather than a pointwise S diff.
                         "port-S", "error-indicators.csv"};
  opts.custom_checks["rom-eigenvalues.csv"] =
      TestRomEigenvalueMatchesEigenmode(6.690, 1.0e-3, 22.409, 3.0e-3);
  opts.paraview_fields = false;
  palace::test::RunRegressionCase("iris_filter", "driven_wave_synth.json",
                                  "driven_wave_synth", opts);
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
  opts.excluded_columns = eigen_excluded;
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
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  palace::test::RunRegressionCase("cpw", "cpw_wave_eigen.json", "wave_eigen", opts);
}

// --- adapter: same reltol / atol as the cpw eigen block above ---

TEST_CASE("adapter_hybrid", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-11;
  opts.excluded_columns = eigen_excluded;
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
  // End-to-end de-embed/reconnect: reduce the total synthesized pencil to the port
  // terminals, de-embed the reference loads, and confirm the resulting S reproduces the
  // solver's field-derived port-S (a stronger check than the pointwise field-S diff, and
  // one that also catches the dB-scale shift a regressed W would produce). The physical S
  // is partition-independent: |S| matches to ~6e-8 at np=1 and np=2, so a 1e-6 linear
  // tolerance is robust.
  opts.custom_checks["port-S.csv"] = TestWavePortSRoundTrip(1.0e-6);
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
  opts.linear_solver_policy = force_default_solver;
  opts.eigen_solver_policy = force_default_solver;
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
  opts.linear_solver_policy = force_default_solver;
  opts.eigen_solver_policy = force_default_solver;
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
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d.json", "eigenmode", opts);
}

// Coarse 2D driven tolerances (partition- and platform-sensitive).
TEST_CASE("cavity2d_driven", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 2.0e-2;
  opts.atol = 1.0e-8;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = force_default_solver;
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
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("cavity2d", "cavity2d_magnetostatic.json",
                                  "magnetostatic", opts);
}

TEST_CASE("cavity2d_transient", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = {"Maximum", "Minimum"};
  opts.linear_solver_policy = force_default_solver;
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
  opts.linear_solver_policy = force_default_solver;
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
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("cpw2d", "cpw2d_thick_impedance.json", "thick_impedance",
                                  opts);
}

// 2D mode analysis from a 3D mesh: back to reltol=1e-4, atol still 1e-10.
TEST_CASE("cpw_wave_2dmode", "[Serial][Parallel][GPU][Regression]")
{
  palace::test::RegressionOptions opts;
  opts.rtol = 1.0e-4;
  opts.atol = 1.0e-10;
  opts.excluded_columns = eigen_excluded;
  opts.skip_rowcount = true;
  opts.linear_solver_policy = force_default_solver;
  palace::test::RunRegressionCase("cpw", "cpw_wave_2dmode.json", "wave_2dmode", opts);
}
