// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <optional>
#include <unordered_set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "models/surfaceterms.hpp"
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

// The frequency-dependent boundary mass coefficients (real and imaginary parts) of the
// surface terms at ω, for the tangential trace (vector) or the out-of-plane (scalar) field.
using SummedMethod = void (SurfaceTerms::*)(double, MaterialPropertyCoefficient &,
                                            MaterialPropertyCoefficient &);
const SummedMethod summed_methods[] = {
    &SurfaceTerms::AddExtraSystemBdrMassCoefficients,
    &SurfaceTerms::AddScalarExtraSystemBdrMassCoefficients};

}  // namespace

TEST_CASE("SurfaceRationalImpedanceOperator",
          "[surfacerationalimpedanceoperator][Serial][Parallel]")
{
  // The frequency-dependent boundary mass coefficients of a rational impedance boundary,
  // as SurfaceTerms adds them to the coefficient of each field (from the scale and unit
  // coefficient of its term).
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

  // The surface terms of the rational impedance boundary, which add its frequency-dependent
  // boundary mass coefficients.
  auto make_terms =
      [&](const config::RationalImpedanceData &data, const std::unordered_set<int> &cracked)
  {
    config::BoundaryData boundaries;
    boundaries.rational_impedance = {data};
    boundaries.split_impedance_attributes = cracked;
    return std::make_unique<SurfaceTerms>(boundaries, ProblemType::DRIVEN, units, mat_op,
                                          palace_mesh.Get());
  };

  SECTION("Matches the parallel-RLC admittance at finite frequency")
  {
    auto terms = make_terms(make_data(num, den), {});

    for (double omega : {1.0, 4.0, 8.0})
    {
      for (auto summed : summed_methods)
      {
        MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
        MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
        ((*terms).*summed)(omega, fbr, fbi);

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
  }

  SECTION("Matches the parallel-RLC admittance at complex frequency")
  {
    // The analytic continuation iω·Ys = iω/Rs + 1/Ls - ω²Cs, at the complex frequencies of
    // the eigenmode solve.
    auto terms = make_terms(make_data(num, den), {});

    for (std::complex<double> omega :
         {std::complex<double>(1.0, 0.2), std::complex<double>(4.0, -0.5)})
    {
      MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
          dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
          fbi(mat_op.MaxCeedBdrAttribute()), frbr(mat_op.MaxCeedBdrAttribute()),
          frbi(mat_op.MaxCeedBdrAttribute());
      terms->AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, nullptr,
                                           nullptr);

      const std::complex<double> coef =
          std::complex<double>(0.0, 1.0) * omega / Rs + 1.0 / Ls - omega * omega * Cs;
      auto vr = GetBdrCoeffValue(fbr, palace_mesh, 1);
      auto vi = GetBdrCoeffValue(fbi, palace_mesh, 1);
      if (vr)
      {
        CHECK_THAT(*vr, WithinRel(coef.real(), 1e-12));
      }
      if (vi)
      {
        CHECK_THAT(*vi, WithinRel(coef.imag(), 1e-12));
      }
      require_global_coverage(vr.has_value(), vi.has_value());
    }
  }

  SECTION("Impedance pole: contribution vanishes without inf/nan")
  {
    // Lossless parallel LC: Zs = sLs/(s²LsCs + 1) has a pole at ω0 = 1/sqrt(Ls*Cs) = 4,
    // where D(iω0) = 0. The direct admittance form gives coef = s·D/N = 0 there.
    auto terms = make_terms(make_data({Ls, 0.0}, {Ls * Cs, 0.0, 1.0}), {});

    for (auto summed : summed_methods)
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      ((*terms).*summed)(4.0, fbr, fbi);

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
  }

  SECTION("Per-attribute scaling, mixed cracked/uncracked")
  {
    auto terms = make_terms(make_data(num, den), {2});

    const double omega = 1.0;
    for (auto summed : summed_methods)
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      ((*terms).*summed)(omega, fbr, fbi);

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
  }

  SECTION("Invalid inputs are rejected")
  {
    std::unordered_set<int> cracked = {};

    // Transmission zero (Zs = 0) at the evaluation frequency: series LC with
    // N = [Ls*Cs, 0, 1], zero at ω0 = 4.
    auto terms_zero = make_terms(make_data({Ls * Cs, 0.0, 1.0}, {Cs, 0.0}), cracked);
    for (auto summed : summed_methods)
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      MaterialPropertyCoefficient fbi(mat_op.MaxCeedBdrAttribute());
      CHECK_THROWS_WITH(((*terms_zero).*summed)(4.0, fbr, fbi),
                        ContainsSubstring("Rational impedance boundary has a transmission "
                                          "zero (Zs = 0) at the evaluation frequency"));
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