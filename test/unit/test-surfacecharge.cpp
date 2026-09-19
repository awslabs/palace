// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

// Electric surface flux Φ = ∫ D ⋅ n dS and surface charge for a lossy dielectric. Palace's
// constitutive relation for a material with loss tangent is the complex permittivity
// ε = ε' (1 - i tan δ) (Im{ε} = -ε' tan δ, e^{iωt} convention), so for a complex
// peak-phasor field D = ε E has D_r = ε' E_r - ε'' E_i and D_i = ε' E_i + ε'' E_r with ε''
// = -ε' tan δ. The expected values below are written directly from this definition
// (internal units, ε₀ = 1): for a uniform field E = (0, 0, E_z) through the face z = s_z of
// the cube with outward normal +z, Φ = ε_r (1 - i tan δ) E_z s_x s_y.

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("Electric surface flux of a lossy dielectric uses the complex permittivity",
          "[surfacepostoperator][Serial][Parallel]")
{
  constexpr int order = 1;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  const double tandelta = GENERATE(0.0, 0.3);
  constexpr double eps_r = 4.0;
  const std::complex<double> Ez(1.0, 0.5);  // Uniform complex peak phasor, internal units.
  CAPTURE(tandelta);

  // Cube [0, sx] x [0, sy] x [0, sz] in internal units; MakeCartesian3D boundary attribute
  // 6 is the top face z = sz.
  constexpr double sx = 1.1, sy = 2.5, sz = 3.8;
  constexpr int top_attr = 6;
  MPI_Comm comm = Mpi::World();
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, sx, sy, sz));
  const int dim = serial_mesh->Dimension();
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, *serial_mesh));
  serial_mesh.reset();

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {eps_r, eps_r, eps_r};
  material.tandelta.s = {tandelta, tandelta, tandelta};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);

  GridFunction E(nd_fespace, true);
  {
    mfem::Vector vr(3), vi(3);
    vr = 0.0;
    vi = 0.0;
    vr(2) = Ez.real();
    vi(2) = Ez.imag();
    mfem::VectorConstantCoefficient cr(vr), ci(vi);
    E.Real().ProjectCoefficient(cr);
    E.Imag().ProjectCoefficient(ci);
  }

  config::BoundaryPostData postpro;
  config::SurfaceFluxData flux;
  flux.type = SurfaceFlux::ELECTRIC;
  flux.two_sided = false;
  flux.center = {0.5 * sx, 0.5 * sy, 0.5 * sz};
  flux.no_center = false;
  flux.attributes = {top_attr};
  postpro.flux.emplace(1, flux);
  SurfacePostOperator surf_post_op(postpro, ProblemType::DRIVEN, mat_op, h1_fespace,
                                   nd_fespace);

  // Φ = ∫ D ⋅ n dS with D = ε_r (1 - i tan δ) E and n = +z on the top face.
  const std::complex<double> eps(eps_r, -eps_r * tandelta);
  const std::complex<double> expected = eps * Ez * sx * sy;
  const auto flux_value = surf_post_op.GetSurfaceFlux(1, &E, nullptr);
  CHECK_THAT(flux_value.real(), WithinRel(expected.real(), 1.0e-10));
  CHECK_THAT(flux_value.imag(), WithinRel(expected.imag(), 1.0e-10));
}

TEST_CASE("Lossy parallel-plate capacitor conserves charge at the driven port",
          "[surfacepostoperator][Serial][Parallel]")
{
  // A 1 mm cube of lossy dielectric between PEC plates at z = 0 and z = 1 mm (a section of
  // an infinite parallel-plate capacitor; the remaining faces are PMC), driven through a
  // 50 Ω lumped port bridging the plates on the face x = 0 at 1 GHz (electrically small).
  // Charge conservation on the top electrode relates the surface charge Q = ∫ D ⋅ n dS to
  // the net current delivered by the port, I = (2 V_inc - V) / R (Thevenin source 2 V_inc
  // behind R; V the total port voltage): i ω Q = I. This holds independently of how the
  // constitutive relation is written, so it discriminates between D = ε' E and the correct
  // D = ε' (1 - i tan δ) E for a lossy dielectric. The port admittance is additionally
  // compared with the lossy capacitor i ω C (1 - i tan δ), C = ε_r ε₀ A / d.
  using namespace std::complex_literals;
  using json = nlohmann::json;
  const double tandelta = GENERATE(0.0, 0.3);
  CAPTURE(tandelta);
  MPI_Comm comm = Mpi::World();

  // Mesh units are millimetres; L0 = Lc = 1 mm so the mesh is already in internal units.
  IoData iodata{Units(1.0e-3, 1.0e-3)};
  iodata.model.L0 = 1.0e-3;
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  iodata.solver.order = 3;
  iodata.solver.driven.sample_f = {1.0};  // GHz, nondimensionalized below
  iodata.domains = config::DomainData(
      json{{"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                                    {"Permittivity", 4.0},
                                                    {"LossTan", tandelta}})})}});
  // MakeCartesian3D boundary attributes: 1 (z = 0), 6 (z = sz), 5 (x = 0).
  iodata.boundaries = config::BoundaryData(json{
      {"PEC", json::object({{"Attributes", json::array({1, 6})}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", 1},
                                                {"Attributes", json::array({5})},
                                                {"Direction", "+Z"}})})},
      {"Postprocessing",
       json::object(
           {{"SurfaceFlux",
             json::array({json::object({{"Index", 1},
                                        {"Type", "Electric"},
                                        {"Attributes", json::array({6})},
                                        {"Center", json::array({0.5, 0.5, 0.5})}})})}})}});
  iodata.CheckConfiguration();

  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
  SpaceOperator space_op(iodata, mesh);

  // Solve (K + i ω C - ω² M) E = RHS at the single frequency.
  const double omega = iodata.solver.driven.sample_f.front();
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ksp.SetOperators(*A, *P);
  ComplexVector RHS, E;
  RHS.UseDevice(true);
  E.UseDevice(true);
  REQUIRE(space_op.GetExcitationVector(1, omega, RHS));
  E.SetSize(RHS.Size());
  E = 0.0;
  ksp.Mult(RHS, E);

  GridFunction E_gf(space_op.GetNDSpace(), true);
  E_gf.Real().SetFromTrueDofs(E.Real());
  E_gf.Imag().SetFromTrueDofs(E.Imag());
  E_gf.Real().ExchangeFaceNbrData();
  E_gf.Imag().ExchangeFaceNbrData();

  // Net current into the top plate from the port and the plate charge.
  const auto &port = space_op.GetLumpedPortOp().GetPort(1);
  const std::complex<double> V = port.GetVoltage(E_gf);
  const double V_inc = port.GetExcitationVoltage();
  const std::complex<double> I_net = (2.0 * V_inc - V) / port.R;
  SurfacePostOperator surf_post_op(iodata.boundaries.postpro, ProblemType::DRIVEN,
                                   space_op.GetMaterialOp(), space_op.GetH1Space(),
                                   space_op.GetNDSpace());
  const std::complex<double> Q = surf_post_op.GetSurfaceFlux(1, &E_gf, nullptr);
  const std::complex<double> ratio = 1i * omega * Q / I_net;
  CAPTURE(V, V_inc, I_net, Q, ratio);
  CHECK_THAT(ratio.real(), WithinRel(1.0, 1.0e-3));
  CHECK_THAT(ratio.imag(), WithinAbs(0.0, 1.0e-3));

  // Setup control: the port sees the lossy capacitor Y = i ω C (1 - i tan δ) with
  // C = ε_r A / d in internal units (ε₀ = 1, A = d² for the unit cube).
  const std::complex<double> Y = I_net / V;
  const std::complex<double> Y_expected = 1i * omega * 4.0 * (1.0 - 1i * tandelta);
  CAPTURE(Y, Y_expected);
  CHECK_THAT(std::abs(Y / Y_expected - 1.0), WithinAbs(0.0, 1.0e-2));
}

}  // namespace palace
