// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Unit tests for Floquet port computational functions: the static methods in
// FloquetPortData, and the S-parameters of a port which is not driven.

#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <memory>
#include <numbers>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/floquetportoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

using namespace Catch::Matchers;
using namespace std::complex_literals;

// =============================================================================
// Reciprocal lattice: FloquetPortData::ComputeReciprocalLattice
// =============================================================================
TEST_CASE("Floquet Reciprocal Lattice", "[floquetport][Serial]")
{
  SECTION("Square cell: a₁ = (1,0,0), a₂ = (0,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 1.0;
    a2 = 0.0;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
    CHECK_THAT(b1.Norml2(), WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(b2.Norml2(), WithinRel(2.0 * std::numbers::pi, 1e-10));
  }

  SECTION("Rectangular cell: a₁ = (4,0,0), a₂ = (0,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 4.0;
    a2 = 0.0;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
    // |b₁| = 2π/4 = π/2, |b₂| = 2π
    CHECK_THAT(b1.Norml2(), WithinRel(std::numbers::pi / 2.0, 1e-10));
    CHECK_THAT(b2.Norml2(), WithinRel(2.0 * std::numbers::pi, 1e-10));
    // b₁ along x, b₂ along y
    CHECK_THAT(b1(0), WithinRel(std::numbers::pi / 2.0, 1e-10));
    CHECK_THAT(std::abs(b1(1)), WithinAbs(0.0, 1e-10));
    CHECK_THAT(b2(1), WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(std::abs(b2(0)), WithinAbs(0.0, 1e-10));
  }

  SECTION("Non-axis-aligned cell: a₁ = (1,0,0), a₂ = (0.5,1,0)")
  {
    mfem::Vector a1(3), a2(3), b1(3), b2(3);
    a1 = 0.0;
    a1(0) = 1.0;
    a2 = 0.0;
    a2(0) = 0.5;
    a2(1) = 1.0;

    FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);

    CHECK_THAT(a1 * b1, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(a2 * b2, WithinRel(2.0 * std::numbers::pi, 1e-10));
    CHECK_THAT(std::abs(a1 * b2), WithinAbs(0.0, 1e-10));
    CHECK_THAT(std::abs(a2 * b1), WithinAbs(0.0, 1e-10));
  }
}

// =============================================================================
// BZ wrapping offset: FloquetPortData::ComputeBZOffset
// =============================================================================
TEST_CASE("Floquet BZ Wrapping Offset", "[floquetport][Serial]")
{
  // Period a = 4 in x → b₁ = 2π/4 = π/2 in x. BZ: [-π/4, π/4].
  mfem::Vector a1(3), a2(3), b1(3), b2(3);
  a1 = 0.0;
  a1(0) = 4.0;
  a2 = 0.0;
  a2(1) = 1.0;
  FloquetPortData::ComputeReciprocalLattice(a1, a2, b1, b2);
  double b1_sq = b1 * b1;

  SECTION("k_F inside BZ: no wrapping")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = 0.5;
    kF_w = kF;  // No wrapping needed

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == 0);
  }

  SECTION("k_F outside BZ: wrapping by +1 b₁")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = 1.05;  // Outside π/4 ≈ 0.785
    kF_w = 0.0;
    kF_w(0) = std::remainder(1.05, 2.0 * std::numbers::pi / 4.0);  // Wrapped

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == 1);
    // Verify roundtrip: kF_w + bz_m * b₁ ≈ kF
    CHECK_THAT(kF_w(0) + bz_m * b1(0), WithinRel(kF(0), 1e-12));
  }

  SECTION("Negative k_F outside BZ: wrapping by -1 b₁")
  {
    mfem::Vector kF(3), kF_w(3);
    kF = 0.0;
    kF(0) = -1.05;
    kF_w = 0.0;
    kF_w(0) = std::remainder(-1.05, 2.0 * std::numbers::pi / 4.0);

    int bz_m = FloquetPortData::ComputeBZOffset(kF, kF_w, b1, b1_sq);
    CHECK(bz_m == -1);
    CHECK_THAT(kF_w(0) + bz_m * b1(0), WithinRel(kF(0), 1e-12));
  }
}

// =============================================================================
// A port which is not driven, with its specular order below cutoff
// =============================================================================
TEST_CASE_METHOD(test::SharedTempDir, "Floquet port terminating an evanescent field",
                 "[floquetport][Serial][Parallel]")
{
  // Total internal reflection of a plane wave incident at 45° from a dielectric, ε = 4, at
  // the excited port 1 (z = 0) on vacuum (z > L / 2), which port 2 (z = L) terminates. The
  // tangential wave vector k_F = 2 k₀ sin 45° = √2 k₀ exceeds the vacuum wavenumber k₀, so
  // the specular order is below cutoff at port 2, which only terminates the evanescent
  // field, γ₂ = −i k₀. Port 2 is not excited and needs no incident wave: it has no order
  // which carries power, and so no S-parameters. The reflection is total,
  // S11 = r exp(−2i γ₁ L / 2) with γ₁ = √2 k₀ and the Fresnel coefficient
  // r = (λ₁ − λ₂) / (λ₁ + λ₂) of the DtN eigenvalues λ = γ (TE) and ω²ε/γ (TM), |r| = 1.
  using json = nlohmann::json;
  const bool is_te = GENERATE(true, false);
  CAPTURE(is_te);
  MPI_Comm comm = Mpi::World();
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  iodata.model.L0 = 1.0e-3;
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  iodata.solver.order = 3;
  iodata.solver.driven.sample_f = {100.0};  // GHz
  // Precondition with a sparse direct solve of the complex-valued system, without the
  // low-rank port terms: GMRES converges in a few iterations. (With multigrid and its
  // real-valued coarse solve, the TM solve at oblique incidence converges slowly.)
  iodata.solver.linear.mg_max_levels = 1;
  iodata.solver.linear.complex_coarse_solve = true;
  iodata.solver.linear.tol = 1.0e-10;
  // Wavenumbers in 1/mm, the mesh length unit: k₀ of vacuum and the tangential k_F.
  constexpr double eps1 = 4.0, eps2 = 1.0;
  const double k0 = 2.0 * std::numbers::pi * 1.0e11 * 1.0e-3 / electromagnetics::c0_;
  const double kF = std::sqrt(eps1) * k0 * std::sin(std::numbers::pi / 4.0);
  iodata.domains = config::DomainData(
      json{{"Materials",
            json::array({{{"Attributes", json::array({1})}, {"Permittivity", eps1}},
                         {{"Attributes", json::array({2})}, {"Permittivity", eps2}}})}});
  // MakeCartesian3D boundary attributes: 1 (z = 0), 6 (z = L), 5 and 3 (x faces), 2 and
  // 4 (y faces).
  iodata.boundaries = config::BoundaryData(json{
      {"Periodic",
       {{"FloquetWaveVector", json::array({kF, 0.0, 0.0})},
        {"BoundaryPairs", json::array({{{"DonorAttributes", json::array({5})},
                                        {"ReceiverAttributes", json::array({3})}},
                                       {{"DonorAttributes", json::array({2})},
                                        {"ReceiverAttributes", json::array({4})}}})}}},
      {"FloquetPort",
       json::array({{{"Index", 1},
                     {"Attributes", json::array({1})},
                     {"Excitation", 1},
                     {"IncidentPolarization", is_te ? "TE" : "TM"},
                     {"MaxOrder", 0}},
                    {{"Index", 2}, {"Attributes", json::array({6})}, {"MaxOrder", 0}}})}});
  iodata.CheckConfiguration();
  constexpr double a = 0.5, L = 0.25;
  iodata.model.mesh = (temp_dir / "floquet_tir.mesh").string();
  if (Mpi::Root(comm))
  {
    // The dielectric (attribute 1) below z = L / 2, vacuum (attribute 2) above.
    auto box = mfem::Mesh::MakeCartesian3D(3, 3, 8, mfem::Element::HEXAHEDRON, a, a, L);
    mfem::Vector center(3);
    for (int e = 0; e < box.GetNE(); e++)
    {
      box.GetElementCenter(e, center);
      box.SetAttribute(e, (center(2) < 0.5 * L) ? 1 : 2);
    }
    box.SetAttributes();
    std::ofstream ofs(iodata.model.mesh);
    box.Print(ofs);
  }
  Mpi::Barrier(comm);
  auto smesh = mesh::Load(iodata, comm);
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
  SpaceOperator space_op(iodata, mesh);

  const double omega = iodata.solver.driven.sample_f.front();
  REQUIRE_THAT(omega, WithinRel(k0, 1.0e-12));
  auto &floquet_op = space_op.GetFloquetPortOp();
  floquet_op.Initialize(omega);
  // The (0,0) order of port 2 carries no power, and its Robin reference is zero.
  REQUIRE(floquet_op.GetPort(2).GetGamma0() == 0.0);

  ComplexVector E;
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
    auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                      K.get(), C.get(), M.get(), A2.get());
    auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
        1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
    ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
    ksp.SetOperators(*A, *P);
    ComplexVector RHS;
    RHS.UseDevice(true);
    E.UseDevice(true);
    REQUIRE(space_op.GetExcitationVector(1, omega, RHS));
    E.SetSize(RHS.Size());
    E = 0.0;
    ksp.Mult(RHS, E);
  }
  GridFunction E_gf(space_op.GetNDSpace(), true);
  E_gf.Real().SetFromTrueDofs(E.Real());
  E_gf.Imag().SetFromTrueDofs(E.Imag());

  std::map<std::tuple<int, int, bool>, std::complex<double>> S1, S2;
  REQUIRE_NOTHROW(S1 = floquet_op.GetPort(1).GetAllSParameters(E_gf, true));
  REQUIRE_NOTHROW(S2 = floquet_op.GetPort(2).GetAllSParameters(E_gf));
  CHECK(S2.empty());

  const double gamma1 = std::sqrt(eps1 * k0 * k0 - kF * kF);
  const std::complex<double> gamma2(0.0, -std::sqrt(kF * kF - eps2 * k0 * k0));
  const std::complex<double> lambda1 = is_te ? gamma1 : eps1 * k0 * k0 / gamma1;
  const std::complex<double> lambda2 = is_te ? gamma2 : eps2 * k0 * k0 / gamma2;
  const auto r = (lambda1 - lambda2) / (lambda1 + lambda2);
  REQUIRE_THAT(std::abs(r), WithinRel(1.0, 1.0e-14));
  const auto S11_expected = r * std::exp(-1i * gamma1 * L);
  REQUIRE(S1.size() == 2);
  const auto S11 = S1.at({0, 0, is_te}), S11_cross = S1.at({0, 0, !is_te});
  CAPTURE(S11, S11_expected, S11_cross);
  CHECK(std::abs(S11 - S11_expected) < 1.0e-8);
  CHECK(std::abs(S11_cross) < 1.0e-8);
}

}  // namespace palace
