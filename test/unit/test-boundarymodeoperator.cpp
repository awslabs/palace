// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <functional>
#include <numeric>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/singularfeatures.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Solve for modes of a 2D rectangular waveguide cross-section using
// ModeEigenSolver. Returns eigenvalues as complex kn values.
// MakeCartesian2D boundary attributes: bottom=1, right=2, top=3, left=4.
struct ModeResult
{
  std::vector<std::complex<double>> kn;
  int num_converged;
};

mfem::Mesh InternalPecLineMesh()
{
  mfem::Mesh mesh(2, 9, 8, 9, 2);
  mesh.AddVertex(-1.0, -1.0);
  mesh.AddVertex(0.0, -1.0);
  mesh.AddVertex(1.0, -1.0);
  mesh.AddVertex(-1.0, 0.0);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(-1.0, 1.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);
  mesh.AddTriangle(0, 1, 4, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(3, 4, 7, 1);
  mesh.AddTriangle(3, 7, 6, 1);
  mesh.AddTriangle(1, 2, 5, 1);
  mesh.AddTriangle(1, 5, 4, 1);
  mesh.AddTriangle(4, 5, 8, 1);
  mesh.AddTriangle(4, 8, 7, 1);
  mesh.AddBdrSegment(0, 1, 1);
  mesh.AddBdrSegment(1, 2, 1);
  mesh.AddBdrSegment(2, 5, 1);
  mesh.AddBdrSegment(5, 8, 1);
  mesh.AddBdrSegment(8, 7, 1);
  mesh.AddBdrSegment(7, 6, 1);
  mesh.AddBdrSegment(6, 3, 1);
  mesh.AddBdrSegment(3, 0, 1);
  mesh.AddBdrSegment(3, 4, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

void SetQuadraticCurvedGeometry(mfem::Mesh &mesh)
{
  mesh.SetCurvature(2, false, 2, mfem::Ordering::byVDIM);
  mfem::VectorFunctionCoefficient geometry(2,
                                           [](const mfem::Vector &x, mfem::Vector &value)
                                           {
                                             value.SetSize(2);
                                             value = x;
                                             value[1] += 0.1 * x[0] * x[1];
                                           });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

double RelativeNorm(const mfem::Vector &residual, const mfem::Vector &reference)
{
  return residual.Norml2() / std::max(1.0, reference.Norml2());
}

ModeResult SolveRectangularModes(double width, double height, double freq_ghz,
                                 double epsilon_r, int order, int num_modes,
                                 const std::function<void(IoData &)> &configure_bcs)
{
  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;

  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {epsilon_r, epsilon_r, epsilon_r};

  // Default: PEC on all boundaries.
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};

  // Let the caller configure specific BCs.
  configure_bcs(iodata);

  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = freq_ghz;
  iodata.solver.boundary_mode.n = num_modes;
  iodata.solver.boundary_mode.tol = 1.0e-8;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 200;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(10, 5, mfem::Element::TRIANGLE, false, width, height));
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.CheckConfiguration();
  Mesh palace_mesh(std::move(par_mesh));

  auto nd_fec = std::make_unique<mfem::ND_FECollection>(order, palace_mesh.Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(order, palace_mesh.Dimension());
  FiniteElementSpace nd_fespace(palace_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(palace_mesh, h1_fec.get());
  MaterialOperator mat_op(iodata, palace_mesh);

  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, palace_mesh.Get());
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, palace_mesh.Get());
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, palace_mesh.Get());
  SurfaceRationalImpedanceOperator surf_rz_op(iodata, mat_op, palace_mesh.Get());

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = palace_mesh.Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }

  int nd_size = nd_fespace.GetTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
  double kn_target = omega * std::sqrt(1.1 * mat_op.GetMaxMuEpsilon());

  // ModeEigenSolver requires a positive Krylov subspace size (num_vec). Mirror the
  // formula used by IoData::CheckConfiguration for eigenmode.max_size.
  const int num_vec = std::max(2 * num_modes, num_modes + 15);
  ModeEigenSolver mode_solver(mat_op, nullptr, surf_z_op, farfield_op, surf_sigma_op,
                              surf_rz_op, nd_fespace, h1_fespace, dbc_tdof_list, num_modes,
                              num_vec, 1.0e-8, EigenvalueSolver::WhichType::LARGEST_REAL,
                              iodata.solver.linear, iodata.solver.boundary_mode.type, 0,
                              nd_fespace.GetComm());

  double sigma = -kn_target * kn_target;
  auto result = mode_solver.Solve(omega, sigma);

  ModeResult out;
  out.num_converged = result.num_converged;
  for (int i = 0; i < result.num_converged; i++)
  {
    auto lambda = mode_solver.GetEigenvalue(i);
    out.kn.push_back(std::sqrt(-sigma - 1.0 / lambda));
  }
  return out;
}

}  // namespace

TEST_CASE("Singular BoundaryMode rejects unsupported material and boundary physics",
          "[boundarymodeoperator][singularelements][Serial]")
{
  Units units(1.0, 1.0);
  const auto make_input = [&]()
  {
    IoData iodata(units);
    iodata.problem.type = ProblemType::BOUNDARYMODE;
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    iodata.boundaries.pec.attributes = {1, 7};
    iodata.solver.linear.mg_max_levels = 1;
    iodata.solver.singular_elements.attributes = {7};
    return iodata;
  };

  SECTION("Anisotropic material")
  {
    auto iodata = make_input();
    iodata.domains.materials[0].epsilon_r.s = {1.0, 1.0, 2.0};
    CHECK_THROWS(iodata.CheckConfiguration());
  }
  SECTION("Lossy or London material")
  {
    auto iodata = make_input();
    iodata.domains.materials[0].tandelta.s = {1.0e-4, 1.0e-4, 1.0e-4};
    CHECK_THROWS(iodata.CheckConfiguration());

    iodata = make_input();
    iodata.domains.materials[0].lambda_L = 1.0e-7;
    CHECK_THROWS(iodata.CheckConfiguration());
  }
  SECTION("Non-PEC boundary operator")
  {
    auto iodata = make_input();
    auto &impedance = iodata.boundaries.impedance.emplace_back();
    impedance.attributes = {2};
    impedance.Ls = 1.0e-9;
    CHECK_THROWS(iodata.CheckConfiguration());
  }
}

TEST_CASE("Singular BoundaryMode blocks preserve the complete exact sequence",
          "[boundarymodeoperator][singularelements][curved][Serial]")
{
  MPI_Comm comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) == 1);

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {2.3, 2.3, 2.3};
  material.mu_r.s = {1.7, 1.7, 1.7};
  iodata.boundaries.pec.attributes = {1, 7};
  iodata.solver.order = 1;
  iodata.solver.boundary_mode.freq = 10.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.linear.mg_max_levels = 1;
  iodata.solver.singular_elements.attributes = {7};
  iodata.solver.singular_elements.order = 1;
  iodata.solver.singular_elements.abs_tol = 2.0e-10;
  iodata.solver.singular_elements.rel_tol = 2.0e-10;

  auto serial_mesh = std::make_unique<mfem::Mesh>(InternalPecLineMesh());
  SetQuadraticCurvedGeometry(*serial_mesh);
  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(*serial_mesh, {7});
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(serial_mesh->GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(serial_mesh->GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);

  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *par_mesh, source_vertex_ids, source_element_ids);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  MaterialOperator mat_op(iodata, *meshes.back());
  BoundaryModeOperator mode_op(iodata, meshes, mat_op, &local_features, &source_vertex_ids);

  REQUIRE(mode_op.HasSingularEnrichment());
  REQUIRE(mode_op.GetCombinedGradient());
  REQUIRE(mode_op.GetBtt());
  REQUIRE(mode_op.GetAtnr());
  REQUIRE(mode_op.GetBtnr());
  CHECK(mode_op.GetNDTrueVSize() > mode_op.GetNDSpace().GetTrueVSize());
  CHECK(mode_op.GetH1TrueVSize() > mode_op.GetH1Space().GetTrueVSize());
  CHECK(mode_op.GetCombinedDbcTDofList().Size() >
        mode_op.GetNDDbcTDofLists().back().Size() +
            mode_op.GetH1DbcTDofLists().back().Size());

  const auto &gradient = *mode_op.GetCombinedGradient();
  const auto &btt = *mode_op.GetBtt();
  const auto &atn = *mode_op.GetAtnr();
  const auto &btn = *mode_op.GetBtnr();
  mfem::Vector h1(mode_op.GetH1TrueVSize());
  mfem::Vector nd(mode_op.GetNDTrueVSize());
  for (int i = 0; i < h1.Size(); i++)
  {
    h1[i] = 0.17 + 0.031 * i;
  }
  for (int i = 0; i < nd.Size(); i++)
  {
    nd[i] = -0.23 + 0.019 * i;
  }

  mfem::Vector gradient_h1(nd.Size()), btt_gradient_h1(nd.Size()), atn_h1(nd.Size());
  gradient.Mult(h1, gradient_h1);
  btt.Mult(gradient_h1, btt_gradient_h1);
  atn.Mult(h1, atn_h1);
  atn_h1 += btt_gradient_h1;
  CHECK(RelativeNorm(atn_h1, btt_gradient_h1) < 2.0e-11);

  mfem::Vector btn_nd(h1.Size()), atn_transpose_nd(h1.Size());
  btn.Mult(nd, btn_nd);
  atn.MultTranspose(nd, atn_transpose_nd);
  btn_nd += atn_transpose_nd;
  CHECK(RelativeNorm(btn_nd, atn_transpose_nd) < 2.0e-13);

  auto [att, att_imag] = mode_op.AssembleAtt(0.0, 0.0);
  auto [ann, ann_imag] = mode_op.AssembleAnn(0.0);
  REQUIRE(att);
  REQUIRE(ann);
  CHECK_FALSE(att_imag);
  CHECK_FALSE(ann_imag);
  REQUIRE(att->Height() == nd.Size());
  REQUIRE(ann->Height() == h1.Size());

  mfem::Vector curl_gradient(nd.Size());
  att->Mult(gradient_h1, curl_gradient);
  CHECK(RelativeNorm(curl_gradient, gradient_h1) < 2.0e-10);

  mfem::Vector ann_h1(h1.Size()), gradient_transpose_btt_gradient(h1.Size());
  ann->Mult(h1, ann_h1);
  gradient.MultTranspose(btt_gradient_h1, gradient_transpose_btt_gradient);
  ann_h1 += gradient_transpose_btt_gradient;
  CHECK(RelativeNorm(ann_h1, gradient_transpose_btt_gradient) < 2.0e-10);

  constexpr int num_modes = 4;
  constexpr int num_vectors = 19;
  const double omega = 2.0 * M_PI *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(
                           iodata.solver.boundary_mode.freq);
  const double kn_target =
      omega * std::sqrt(1.1 * mode_op.GetMaterialOp().GetMaxMuEpsilon());
  ModeEigenSolver eig(mode_op, mode_op.GetCombinedDbcTDofList(), num_modes, num_vectors,
                      1.0e-8, EigenvalueSolver::WhichType::LARGEST_REAL,
                      iodata.solver.linear, iodata.solver.boundary_mode.type, 0);
  const auto result = eig.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= 1);

  double maximum_nd_gradient = 0.0;
  double maximum_nd_rotational = 0.0;
  double maximum_h1_gradient = 0.0;
  double maximum_electric_energy = 0.0;
  double maximum_magnetic_normal_energy = 0.0;
  for (int i = 0; i < std::min(result.num_converged, num_modes); i++)
  {
    const auto kn = eig.GetPropagationConstant(i);
    CHECK(std::isfinite(kn.real()));
    CHECK(std::isfinite(kn.imag()));
    CHECK(eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD) < 1.0e-6);

    ComplexVector eigenvector(mode_op.GetNDTrueVSize() + mode_op.GetH1TrueVSize());
    eig.GetEigenvector(i, eigenvector);
    ComplexVector et, en;
    mode_op.ApplyVDBackTransform(eigenvector, kn, et, en);
    if (i == 0)
    {
      const ComplexVector et_before(et);
      const ComplexVector en_before(en);
      const auto power_before = mode_op.ComputePoyntingPower(omega, kn, et, en);
      constexpr double scale = 1.75;
      eigenvector *= scale;
      for (int j = 0; j < et.Size(); j++)
      {
        CHECK(std::abs(et.Real()[j] - scale * et_before.Real()[j]) <=
              2.0e-14 * std::max(1.0, std::abs(scale * et_before.Real()[j])));
        CHECK(std::abs(et.Imag()[j] - scale * et_before.Imag()[j]) <=
              2.0e-14 * std::max(1.0, std::abs(scale * et_before.Imag()[j])));
      }
      for (int j = 0; j < en.Size(); j++)
      {
        CHECK(std::abs(en.Real()[j] - scale * en_before.Real()[j]) <=
              2.0e-14 * std::max(1.0, std::abs(scale * en_before.Real()[j])));
        CHECK(std::abs(en.Imag()[j] - scale * en_before.Imag()[j]) <=
              2.0e-14 * std::max(1.0, std::abs(scale * en_before.Imag()[j])));
      }
      const auto power_after = mode_op.ComputePoyntingPower(omega, kn, et, en);
      CHECK(std::abs(power_after - scale * scale * power_before) <=
            2.0e-12 * std::max({1.0, std::abs(power_after),
                                std::abs(scale * scale * power_before)}));
      eigenvector *= 1.0 / scale;
    }
    const auto norms = mode_op.ComputeSingularCoefficientNorms(et, en);
    maximum_nd_gradient = std::max(maximum_nd_gradient, norms.nd_gradient);
    maximum_nd_rotational = std::max(maximum_nd_rotational, norms.nd_rotational);
    maximum_h1_gradient = std::max(maximum_h1_gradient, norms.h1_gradient);
    const auto energies = mode_op.ComputeSingularFieldEnergies(omega, kn, et, en);
    CHECK(std::isfinite(energies.electric_transverse));
    CHECK(std::isfinite(energies.electric_normal));
    CHECK(std::isfinite(energies.magnetic_transverse));
    CHECK(std::isfinite(energies.magnetic_normal));
    CHECK(energies.electric_transverse >= 0.0);
    CHECK(energies.electric_normal >= 0.0);
    CHECK(energies.magnetic_transverse >= 0.0);
    CHECK(energies.magnetic_normal >= 0.0);
    maximum_electric_energy = std::max(
        maximum_electric_energy, energies.electric_transverse + energies.electric_normal);
    maximum_magnetic_normal_energy =
        std::max(maximum_magnetic_normal_energy,
                 energies.magnetic_transverse + energies.magnetic_normal);
    const auto power = mode_op.ComputePoyntingPower(omega, kn, et, en);
    CHECK(std::isfinite(power.real()));
    CHECK(std::isfinite(power.imag()));
  }
  CHECK(maximum_nd_gradient > 1.0e-10);
  CHECK(maximum_nd_rotational > 1.0e-10);
  CHECK(maximum_h1_gradient > 1.0e-10);
  CHECK(maximum_electric_energy > 1.0e-10);
  CHECK(maximum_magnetic_normal_energy > 1.0e-10);
}

TEST_CASE("ModeEigenSolver PEC", "[boundarymodeoperator][Serial]")
{
  // Rectangular waveguide: 1000×500 μm (L0=1e-6), ε=4, f=500 GHz.
  // Analytical kn for TE10 mode:
  //   kc = π / (a * L0) = π / 1e-3 ≈ 3141.6 1/m
  //   ω = 2π * 500e9 ≈ 3.1416e12 rad/s
  //   kn = sqrt(ω²ε/c² - kc²) = sqrt(4*(π*1e12/c)² - (π/1e-3)²)
  //   In nondimensional units (Lc = a*L0 = 1e-3 m):
  //     kn_nd = kn * Lc
  auto result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  REQUIRE(result.num_converged >= 1);

  double kn_real = result.kn[0].real();
  double kn_imag = result.kn[0].imag();
  CAPTURE(kn_real, kn_imag, result.num_converged);

  // First mode should be propagating (real kn, negligible imaginary part).
  CHECK(kn_real > 0.0);
  CHECK(std::abs(kn_imag) < 1.0e-6 * std::abs(kn_real));

  // Analytical kn for TE10 mode of rectangular waveguide with PEC walls:
  //   a = 1000 μm = 1e-3 m, ε_r = 4, f = 500 GHz
  //   kc = π / a = π / 1e-3 m
  //   kn = sqrt(ω²ε_r/c² - kc²)
  //   kn ≈ 20708 1/m → nondimensional (×Lc where Lc = 1e-6 m) ≈ 0.02071
  // Allow 5% tolerance for the coarse 10×5 mesh at order 2.
  CHECK_THAT(kn_real, WithinRel(0.02071, 0.05));
}

TEST_CASE("ModeEigenSolver Impedance shifts kn", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  // Use a large enough inductance so the impedance shift is well above numerical
  // noise across different BLAS/LAPACK implementations.
  auto imp_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                                          [](IoData &iodata)
                                          {
                                            iodata.boundaries.pec.attributes = {1, 3, 4};
                                            auto &imp =
                                                iodata.boundaries.impedance.emplace_back();
                                            imp.attributes = {2};
                                            imp.Ls = 1.0e-8;
                                          });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(imp_result.num_converged >= 1);

  CHECK(imp_result.kn[0].real() > pec_result.kn[0].real());
}

TEST_CASE("ModeEigenSolver Conductivity adds loss", "[boundarymodeoperator][Serial]")
{
  auto pec_result = SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3, [](IoData &) {});

  auto cond_result =
      SolveRectangularModes(1000.0, 500.0, 500.0, 4.0, 2, 3,
                            [](IoData &iodata)
                            {
                              iodata.boundaries.pec.attributes = {1, 3, 4};
                              auto &cond = iodata.boundaries.conductivity.emplace_back();
                              cond.attributes = {2};
                              cond.sigma = 5.0e7;
                              cond.h = 0.001;
                            });

  REQUIRE(pec_result.num_converged >= 1);
  REQUIRE(cond_result.num_converged >= 1);

  CHECK(std::abs(cond_result.kn[0].imag()) > std::abs(pec_result.kn[0].imag()));
}

}  // namespace palace
