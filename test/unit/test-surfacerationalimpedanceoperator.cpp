// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <complex>
#include <optional>
#include <unordered_set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Returns the scalar coefficient value assigned to bdr_attr, or std::nullopt if the
// attribute is not present on this rank.
std::optional<double> GetBdrCoeffValue(const MaterialPropertyCoefficient &fb,
                                       const Mesh &mesh, int bdr_attr)
{
  auto ceed_attrs = mesh.GetCeedBdrAttributes(bdr_attr);
  if (ceed_attrs.Size() == 0)
  {
    return std::nullopt;
  }
  int ceed_attr = ceed_attrs[0];
  const auto &attr_mat = fb.GetAttributeToMaterial();
  int mat_idx = attr_mat[ceed_attr - 1];
  if (mat_idx < 0)
  {
    return std::nullopt;
  }
  return fb.GetMaterialProperties()(0, 0, mat_idx);
}

}  // namespace

TEST_CASE("SurfaceRationalImpedanceOperator",
          "[surfacerationalimpedanceoperator][Serial][Parallel]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  if (Mpi::Size(Mpi::World()) > serial_mesh->GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  int min_ne = palace_mesh.GetNE();
  Mpi::GlobalMin(1, &min_ne, Mpi::World());
  REQUIRE(min_ne > 0);

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  Units units(1.0, 1.0);

  // Nondimensional parallel-RLC surface impedance: Zs(s) = N(s)/D(s) with
  // N = [Ls*Rs, 0] and D = [Ls*Rs*Cs, Ls, Rs]. The Robin coefficient is
  // iω·Ys = iω/Rs + 1/Ls - ω²Cs, i.e. real part 1/Ls - ω²Cs, imaginary part ω/Rs.
  // Resonance at ω0 = 1/sqrt(Ls*Cs) = 4 with these values.
  const double Rs = 2.0, Ls = 0.5, Cs = 0.125;
  const std::vector<double> num = {Ls * Rs, 0.0};
  const std::vector<double> den = {Ls * Rs * Cs, Ls, Rs};

  auto require_global_coverage = [](bool local_v1, bool local_v2)
  {
    bool flags[2] = {local_v1, local_v2};
    Mpi::GlobalOr(2, flags, Mpi::World());
    REQUIRE(flags[0]);
    REQUIRE(flags[1]);
  };

  auto make_data = [&](const std::vector<double> &n, const std::vector<double> &d)
  {
    config::RationalImpedanceData data;
    data.num = n;
    data.den = d;
    data.attributes = {1, 2};
    return data;
  };

  SECTION("Matches the parallel-RLC admittance at finite frequency")
  {
    std::unordered_set<int> cracked = {};
    SurfaceRationalImpedanceOperator op({make_data(num, den)}, cracked, ProblemType::DRIVEN,
                                        units, mat_op, palace_mesh);

    for (double omega : {1.0, 4.0, 8.0})
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

      const double re = 1.0 / Ls - omega * omega * Cs;
      const double im = omega / Rs;
      auto vr = GetBdrCoeffValue(fbr, palace_mesh, 1);
      auto vi = GetBdrCoeffValue(fbi, palace_mesh, 1);
      if (vr)
      {
        CHECK_THAT(*vr, WithinRel(re, 1e-12) || WithinAbs(re, 1e-12));
      }
      if (vi)
      {
        CHECK_THAT(*vi, WithinRel(im, 1e-12));
      }
      require_global_coverage(vr.has_value(), vi.has_value());
    }
  }

  SECTION("DC limit: transmission zero at s = 0 gives iω·Ys → 1/Ls")
  {
    // N(0) = 0 with a simple zero, so lim_{s→0} s·D/N = D(0)/N'(0) = Rs/(Ls·Rs) = 1/Ls.
    std::unordered_set<int> cracked = {};
    SurfaceRationalImpedanceOperator op({make_data(num, den)}, cracked, ProblemType::DRIVEN,
                                        units, mat_op, palace_mesh);

    CHECK_THAT(op.EvalRobinCoefficient(0, 0.0).real(), WithinRel(1.0 / Ls, 1e-12));
    CHECK_THAT(op.EvalRobinCoefficient(0, 0.0).imag(), WithinAbs(0.0, 1e-15));

    MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
    MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
    op.AddExtraSystemBdrCoefficients(0.0, fbr, fbi);

    auto vr = GetBdrCoeffValue(fbr, palace_mesh, 1);
    auto vi = GetBdrCoeffValue(fbi, palace_mesh, 1);
    if (vr)
    {
      CHECK_THAT(*vr, WithinRel(1.0 / Ls, 1e-12));
    }
    if (vi)
    {
      CHECK_THAT(*vi, WithinAbs(0.0, 1e-15));
    }
    require_global_coverage(vr.has_value(), vi.has_value());
  }

  SECTION("DC limit: finite impedance at DC contributes nothing")
  {
    // Purely resistive Zs = Rs: N(0) != 0, so the ω = 0 contribution is exactly zero and
    // no boundary coefficient is assembled at all.
    std::unordered_set<int> cracked = {};
    SurfaceRationalImpedanceOperator op({make_data({Rs}, {1.0})}, cracked,
                                        ProblemType::DRIVEN, units, mat_op, palace_mesh);

    CHECK_THAT(op.EvalRobinCoefficient(0, 0.0).real(), WithinAbs(0.0, 1e-15));
    CHECK_THAT(op.EvalRobinCoefficient(0, 0.0).imag(), WithinAbs(0.0, 1e-15));

    MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
    MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
    op.AddExtraSystemBdrCoefficients(0.0, fbr, fbi);

    bool flags[2] = {GetBdrCoeffValue(fbr, palace_mesh, 1).has_value(),
                     GetBdrCoeffValue(fbr, palace_mesh, 2).has_value()};
    Mpi::GlobalOr(2, flags, Mpi::World());
    CHECK(!flags[0]);
    CHECK(!flags[1]);
  }

  SECTION("Impedance pole: contribution vanishes without inf/nan")
  {
    // Lossless parallel LC: Zs = sLs/(s²LsCs + 1) has a pole at ω0 = 1/sqrt(Ls*Cs) = 4,
    // where D(iω0) = 0. The direct admittance form gives coef = s·D/N = 0 there.
    std::unordered_set<int> cracked = {};
    SurfaceRationalImpedanceOperator op({make_data({Ls, 0.0}, {Ls * Cs, 0.0, 1.0})},
                                        cracked, ProblemType::DRIVEN, units, mat_op,
                                        palace_mesh);

    MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
    MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
    op.AddExtraSystemBdrCoefficients(4.0, fbr, fbi);

    auto vr = GetBdrCoeffValue(fbr, palace_mesh, 1);
    auto vi = GetBdrCoeffValue(fbi, palace_mesh, 1);
    if (vr)
    {
      CHECK(std::isfinite(*vr));
      CHECK_THAT(*vr, WithinAbs(0.0, 1e-12));
    }
    if (vi)
    {
      CHECK(std::isfinite(*vi));
      CHECK_THAT(*vi, WithinAbs(0.0, 1e-12));
    }
    require_global_coverage(vr.has_value(), vi.has_value());
  }

  SECTION("Per-attribute scaling, mixed cracked/uncracked")
  {
    std::unordered_set<int> cracked = {2};
    SurfaceRationalImpedanceOperator op({make_data(num, den)}, cracked, ProblemType::DRIVEN,
                                        units, mat_op, palace_mesh);

    const double omega = 1.0;
    MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
    MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
    op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

    const double re = 1.0 / Ls - omega * omega * Cs;
    auto v1 = GetBdrCoeffValue(fbr, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fbr, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(re / 1.0, 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(re / 2.0, 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }

  SECTION("Pole-residue form of the Robin coefficient")
  {
    std::unordered_set<int> cracked = {};
    std::array<std::complex<double>, 3> poly;
    std::vector<std::complex<double>> poles, residues;

    // Parallel RLC: g(s) = 1/Ls + s/Rs + Cs s², a polynomial (N has a zero at s = 0 which
    // s cancels, and no pole).
    SurfaceRationalImpedanceOperator op_par(
        {make_data(num, den)}, cracked, ProblemType::EIGENMODE, units, mat_op, palace_mesh);
    REQUIRE(op_par.GetRobinPoleResidue(0, poly, poles, residues));
    CHECK(poles.empty());
    CHECK(std::abs(poly[0] - 1.0 / Ls) <= 1.0e-12 / Ls);
    CHECK(std::abs(poly[1] - 1.0 / Rs) <= 1.0e-12 / Rs);
    CHECK(std::abs(poly[2] - Cs) <= 1.0e-12 * Cs);

    // Series RLC: Zs = Rs + s Ls + 1/(s Cs), g(s) = Cs s² / N(s), N(s) = Ls Cs s² + Rs Cs s
    // + 1: g = 1/Ls + Σ Cs p² / N'(p) / (s - p) at the two roots p of N.
    SurfaceRationalImpedanceOperator op_ser({make_data({Ls * Cs, Rs * Cs, 1.0}, {Cs, 0.0})},
                                            cracked, ProblemType::EIGENMODE, units, mat_op,
                                            palace_mesh);
    REQUIRE(op_ser.GetRobinPoleResidue(0, poly, poles, residues));
    CHECK(std::abs(poly[0] - 1.0 / Ls) <= 1.0e-12 / Ls);
    CHECK(poly[1] == 0.0);
    CHECK(poly[2] == 0.0);
    REQUIRE(poles.size() == 2);
    const std::complex<double> disc =
        std::sqrt(std::complex<double>(Rs * Rs * Cs * Cs - 4.0 * Ls * Cs));
    for (const auto p :
         {(-Rs * Cs + disc) / (2.0 * Ls * Cs), (-Rs * Cs - disc) / (2.0 * Ls * Cs)})
    {
      const auto j = (std::abs(poles[0] - p) < std::abs(poles[1] - p)) ? 0 : 1;
      CAPTURE(p, poles[j], residues[j]);
      CHECK(std::abs(poles[j] - p) <= 1.0e-12 * std::abs(p));
      const auto c = Cs * p * p / (2.0 * Ls * Cs * p + Rs * Cs);
      CHECK(std::abs(residues[j] - c) <= 1.0e-12 * std::abs(c));
    }

    // A repeated zero of Zs, N(s) = (s + 1)², and a Robin coefficient s⁴ of degree above
    // two have no such form.
    SurfaceRationalImpedanceOperator op_double({make_data({1.0, 2.0, 1.0}, {1.0})}, cracked,
                                               ProblemType::EIGENMODE, units, mat_op,
                                               palace_mesh);
    CHECK_FALSE(op_double.GetRobinPoleResidue(0, poly, poles, residues));
    SurfaceRationalImpedanceOperator op_quartic({make_data({1.0}, {1.0, 0.0, 0.0, 0.0})},
                                                cracked, ProblemType::EIGENMODE, units,
                                                mat_op, palace_mesh);
    CHECK_FALSE(op_quartic.GetRobinPoleResidue(0, poly, poles, residues));
  }

  SECTION("Invalid inputs are rejected")
  {
    std::unordered_set<int> cracked = {};

    // Transmission zero (Zs = 0) at the evaluation frequency: series LC with
    // N = [Ls*Cs, 0, 1], zero at ω0 = 4.
    SurfaceRationalImpedanceOperator op_zero({make_data({Ls * Cs, 0.0, 1.0}, {Cs, 0.0})},
                                             cracked, ProblemType::DRIVEN, units, mat_op,
                                             palace_mesh);
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      CHECK_THROWS(op_zero.AddExtraSystemBdrCoefficients(4.0, fbr, fbi));
    }

    // Double transmission zero at s = 0: the DC limit diverges.
    SurfaceRationalImpedanceOperator op_dc({make_data({1.0, 0.0, 0.0}, {1.0})}, cracked,
                                           ProblemType::DRIVEN, units, mat_op, palace_mesh);
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      CHECK_THROWS(op_dc.AddExtraSystemBdrCoefficients(0.0, fbr, fbi));
    }

    // Rational impedance boundaries are frequency-domain only: transient (and static)
    // simulation types are rejected, while driven, eigenmode, and boundary mode are
    // allowed.
    CHECK_THROWS(SurfaceRationalImpedanceOperator({make_data(num, den)}, cracked,
                                                  ProblemType::TRANSIENT, units, mat_op,
                                                  palace_mesh));
    CHECK_NOTHROW(SurfaceRationalImpedanceOperator({make_data(num, den)}, cracked,
                                                   ProblemType::EIGENMODE, units, mat_op,
                                                   palace_mesh));
  }
}

}  // namespace palace