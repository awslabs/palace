// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "drivers/singularsolver.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/singulardofs.hpp"
#include "fem/singularfeatures.hpp"
#include "linalg/vector.hpp"
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

void FillVector(mfem::Vector &vector, double phase)
{
  for (int i = 0; i < vector.Size(); i++)
  {
    vector[i] = std::sin(phase + 0.37 * (i + 1));
  }
}

double RelativeGlobalNorm(const mfem::Vector &residual, const mfem::Vector &reference)
{
  return linalg::Norml2(Mpi::World(), residual) /
         std::max(1.0, linalg::Norml2(Mpi::World(), reference));
}

std::complex<double> SolveFundamentalMode(const IoData &iodata,
                                          BoundaryModeOperator &mode_op)
{
  constexpr int num_modes = 4;
  constexpr int num_vectors = 19;
  const double omega = 2.0 * M_PI *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(
                           iodata.solver.boundary_mode.freq);
  const double kn_target =
      iodata.solver.boundary_mode.target > 0.0
          ? iodata.solver.boundary_mode.target * omega
          : omega * std::sqrt(1.1 * mode_op.GetMaterialOp().GetMaxMuEpsilon());
  const auto which = iodata.solver.boundary_mode.target > 0.0
                         ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                         : EigenvalueSolver::WhichType::LARGEST_REAL;
  ModeEigenSolver eig(mode_op, mode_op.GetCombinedDbcTDofList(), num_modes, num_vectors,
                      1.0e-9, which, iodata.solver.linear, iodata.solver.boundary_mode.type,
                      0);
  const auto result = eig.Solve(omega, -kn_target * kn_target);
  REQUIRE(result.num_converged >= 1);
  int selected = -1;
  double selected_distance = std::numeric_limits<double>::infinity();
  for (int mode = 0; mode < result.num_converged; mode++)
  {
    const auto kn = eig.GetPropagationConstant(mode);
    if (ModeEigenSolver::IsPropagating(kn))
    {
      const double distance = std::abs(kn - kn_target);
      if (distance < selected_distance)
      {
        selected = mode;
        selected_distance = distance;
      }
    }
  }
  REQUIRE(selected >= 0);
  CHECK(eig.GetError(selected, EigenvalueSolver::ErrorType::BACKWARD) < 1.0e-7);
  return eig.GetPropagationConstant(selected);
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

TEST_CASE("Singular BoundaryMode preserves resolved propagation constants",
          "[boundarymodeoperator][singularelements][Serial]")
{
  MPI_Comm comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) == 1);

  Units units(1.0, 1.0);
  const auto make_input = [&]()
  {
    IoData iodata(units);
    iodata.model.L0 = 1.0;
    iodata.model.Lc = 2.0;
    iodata.problem.type = ProblemType::BOUNDARYMODE;
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
    material.epsilon_r.s = {2.3, 2.3, 2.3};
    material.mu_r.s = {1.7, 1.7, 1.7};
    iodata.boundaries.pec.attributes = {1, 7};
    iodata.solver.order = 1;
    iodata.solver.boundary_mode.freq = 10.0;
    iodata.solver.boundary_mode.n = 1;
    iodata.solver.boundary_mode.target = 1.0;
    iodata.solver.boundary_mode.tol = 1.0e-9;
    iodata.solver.linear.mg_max_levels = 1;
    iodata.solver.linear.tol = 1.0e-11;
    iodata.solver.linear.max_it = 300;
    return iodata;
  };

  auto singular_iodata = make_input();
  singular_iodata.solver.singular_elements.attributes = {7};
  singular_iodata.solver.singular_elements.order = 1;
  singular_iodata.solver.singular_elements.abs_tol = 2.0e-8;
  singular_iodata.solver.singular_elements.rel_tol = 2.0e-8;

  auto serial_mesh = std::make_unique<mfem::Mesh>(InternalPecLineMesh());
  auto standard_serial_mesh = std::make_unique<mfem::Mesh>(*serial_mesh);
  singular_iodata.NondimensionalizeInputs(serial_mesh);
  singular_iodata.CheckConfiguration();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(*serial_mesh, {7});
  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(serial_mesh->GetNV());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::vector<fem::singular::GlobalVertexId> source_element_ids(serial_mesh->GetNE());
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);

  auto singular_par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *singular_par_mesh, source_vertex_ids, source_element_ids);
  std::vector<std::unique_ptr<Mesh>> singular_meshes;
  singular_meshes.push_back(std::make_unique<Mesh>(std::move(singular_par_mesh)));
  MaterialOperator singular_material(singular_iodata, *singular_meshes.back());
  BoundaryModeOperator singular_op(singular_iodata, singular_meshes, singular_material,
                                   &local_features, &source_vertex_ids);

  auto standard_iodata = make_input();
  standard_iodata.NondimensionalizeInputs(standard_serial_mesh);
  standard_iodata.CheckConfiguration();
  auto standard_par_mesh = std::make_unique<mfem::ParMesh>(comm, *standard_serial_mesh);
  std::vector<std::unique_ptr<Mesh>> standard_meshes;
  standard_meshes.push_back(std::make_unique<Mesh>(std::move(standard_par_mesh)));
  MaterialOperator standard_material(standard_iodata, *standard_meshes.back());
  BoundaryModeOperator standard_op(standard_iodata, standard_meshes, standard_material);

  const auto singular_kn = SolveFundamentalMode(singular_iodata, singular_op);
  const auto standard_kn = SolveFundamentalMode(standard_iodata, standard_op);
  const double omega = 2.0 * M_PI *
                       singular_iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(
                           singular_iodata.solver.boundary_mode.freq);
  const auto singular_neff = singular_kn / omega;
  const auto standard_neff = standard_kn / omega;
  const double kn_relative_difference =
      std::abs(singular_kn - standard_kn) / std::abs(singular_kn);
  const double neff_relative_difference =
      std::abs(singular_neff - standard_neff) / std::abs(singular_neff);
  CAPTURE(singular_kn, standard_kn, singular_neff, standard_neff, kn_relative_difference,
          neff_relative_difference);
  CHECK(kn_relative_difference < 5.0e-4);
  CHECK(neff_relative_difference < 5.0e-4);
}

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
  SECTION("Isotropic loss tangent is supported")
  {
    auto iodata = make_input();
    iodata.domains.materials[0].tandelta.s = {1.0e-4, 1.0e-4, 1.0e-4};
    CHECK_NOTHROW(iodata.CheckConfiguration());
  }
  SECTION("Anisotropic loss tangent or London material")
  {
    auto iodata = make_input();
    iodata.domains.materials[0].tandelta.s = {1.0e-4, 2.0e-4, 1.0e-4};
    CHECK_THROWS(iodata.CheckConfiguration());

    iodata = make_input();
    iodata.domains.materials[0].lambda_L = 1.0e-7;
    CHECK_THROWS(iodata.CheckConfiguration());
  }
  SECTION("Surface impedance is admitted")
  {
    auto iodata = make_input();
    auto &impedance = iodata.boundaries.impedance.emplace_back();
    impedance.attributes = {2};
    impedance.Ls = 1.0e-9;
    CHECK_NOTHROW(iodata.CheckConfiguration());
  }
  SECTION("Surface conductivity remains unsupported")
  {
    auto iodata = make_input();
    auto &conductivity = iodata.boundaries.conductivity.emplace_back();
    conductivity.attributes = {2};
    conductivity.sigma = 5.0e7;
    conductivity.h = 1.0e-7;
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
  constexpr double loss_tangent = 1.0e-4;
  material.tandelta.s = {loss_tangent, loss_tangent, loss_tangent};
  iodata.boundaries.pec.attributes = {1, 7};
  iodata.solver.order = 1;
  iodata.solver.boundary_mode.freq = 10.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.linear.mg_max_levels = 1;
  iodata.solver.linear.tol = 1.0e-12;
  iodata.solver.linear.max_it = 200;
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
  // The complete Btt * G product has the correct action but creates structural fill
  // in the standard-standard block through element-level cancellation. Keep the
  // blockwise assembly sparse without relying on a numerical drop tolerance.
  std::unique_ptr<mfem::HypreParMatrix> full_product(mfem::ParMult(&btt, &gradient, true));
  REQUIRE(full_product);
  CAPTURE(atn.NNZ(), full_product->NNZ());
  CHECK(atn.NNZ() < full_product->NNZ());

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
  REQUIRE(att_imag);
  REQUIRE(ann_imag);
  REQUIRE(att->Height() == nd.Size());
  REQUIRE(ann->Height() == h1.Size());

  mfem::Vector curl_gradient(nd.Size());
  att->Mult(gradient_h1, curl_gradient);
  CHECK(RelativeNorm(curl_gradient, gradient_h1) < 2.0e-10);
  mfem::Vector zero_imaginary_action_nd(nd.Size()), zero_imaginary_action_h1(h1.Size());
  att_imag->Mult(nd, zero_imaginary_action_nd);
  ann_imag->Mult(h1, zero_imaginary_action_h1);
  CHECK(zero_imaginary_action_nd.Norml2() == 0.0);
  CHECK(zero_imaginary_action_h1.Norml2() == 0.0);

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
  auto [att_loss_real, att_loss_imag] = mode_op.AssembleAtt(omega, 0.0);
  auto [ann_loss_real, ann_loss_imag] = mode_op.AssembleAnn(omega);
  REQUIRE(att_loss_real);
  REQUIRE(att_loss_imag);
  REQUIRE(ann_loss_real);
  REQUIRE(ann_loss_imag);
  mfem::Vector att_zero_action(nd.Size()), att_loss_action(nd.Size()),
      att_imag_action(nd.Size()), ann_zero_action(h1.Size()), ann_loss_action(h1.Size()),
      ann_imag_action(h1.Size());
  att->Mult(nd, att_zero_action);
  att_loss_real->Mult(nd, att_loss_action);
  att_loss_imag->Mult(nd, att_imag_action);
  att_loss_action -= att_zero_action;
  att_imag_action.Add(loss_tangent, att_loss_action);
  CHECK(RelativeNorm(att_imag_action, att_loss_action) < 2.0e-10);
  ann->Mult(h1, ann_zero_action);
  ann_loss_real->Mult(h1, ann_loss_action);
  ann_loss_imag->Mult(h1, ann_imag_action);
  ann_loss_action -= ann_zero_action;
  ann_imag_action.Add(loss_tangent, ann_loss_action);
  CHECK(RelativeNorm(ann_imag_action, ann_loss_action) < 2.0e-10);

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

TEST_CASE("Singular BoundaryMode constrains thin-sheet impedance enrichment traces",
          "[boundarymodeoperator][singularelements][impedance][curved][Serial]")
{
  MPI_Comm comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) == 1);

  IoData iodata(Units(1.0, 1.0));
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {2.3, 2.3, 2.3};
  material.mu_r.s = {1.0, 1.0, 1.0};
  iodata.boundaries.pec.attributes = {1};
  auto &impedance = iodata.boundaries.impedance.emplace_back();
  impedance.attributes = {7};
  impedance.Ls = 1.0;
  iodata.solver.order = 1;
  iodata.solver.boundary_mode.freq = 1.0;
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
  REQUIRE(GetConstrainedSingularImpedanceAttributes(iodata, serial_features) ==
          std::set<int>{7});

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
  CHECK(mode_op.GetCombinedDbcTDofList().Size() >
        mode_op.GetNDDbcTDofLists().back().Size() +
            mode_op.GetH1DbcTDofLists().back().Size());
  const double omega = 2.0 * M_PI *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(
                           iodata.solver.boundary_mode.freq);
  auto [att, att_imaginary] = mode_op.AssembleAtt(omega, 0.0);
  auto [ann, ann_imaginary] = mode_op.AssembleAnn(omega);
  REQUIRE(att);
  REQUIRE(ann);
  CHECK_FALSE(att_imaginary);
  CHECK_FALSE(ann_imaginary);

  mfem::Vector input_nd(att->Width()), output_nd(att->Height());
  mfem::Vector input_h1(ann->Width()), output_h1(ann->Height());
  input_nd = 1.0;
  input_h1 = 1.0;
  att->Mult(input_nd, output_nd);
  ann->Mult(input_h1, output_h1);
  CHECK(std::isfinite(output_nd.Norml2()));
  CHECK(std::isfinite(output_h1.Norml2()));
}

TEST_CASE("Singular BoundaryMode multigrid hierarchy commutes",
          "[boundarymodeoperator][singularelements][Serial][Parallel]")
{
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = {2.3, 2.3, 2.3};
  material.mu_r.s = {1.7, 1.7, 1.7};
  iodata.boundaries.pec.attributes = {1, 7};
  iodata.solver.order = 2;
  iodata.solver.boundary_mode.freq = 10.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.linear.type = LinearSolver::AMS;
  iodata.solver.linear.mg_max_levels = 2;
  iodata.solver.linear.pc_mat_real = true;
  iodata.solver.linear.pc_mat_shifted = true;
  iodata.solver.linear.complex_coarse_solve = false;
  iodata.solver.singular_elements.attributes = {7};
  iodata.solver.singular_elements.order = 1;
  iodata.solver.singular_elements.abs_tol = 2.0e-8;
  iodata.solver.singular_elements.rel_tol = 2.0e-8;

  auto serial_mesh = std::make_unique<mfem::Mesh>(InternalPecLineMesh());
  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  const auto serial_features =
      fem::singular::ExtractSerialLineTipFeatures(*serial_mesh, {7});
  std::vector<int> partition(serial_mesh->GetNE());
  for (int element = 0; element < serial_mesh->GetNE(); element++)
  {
    partition[element] = std::min(Mpi::Size(Mpi::World()) - 1,
                                  element * Mpi::Size(Mpi::World()) / serial_mesh->GetNE());
  }
  auto par_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh, partition.data());
  auto source_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
      *serial_mesh, *par_mesh, partition.data());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(par_mesh->GetNE());
  for (int element = 0; element < par_mesh->GetNE(); element++)
  {
    source_element_ids[element] = par_mesh->GetGlobalElementNum(element);
  }
  const auto local_features = fem::singular::DistributeSerialLineTipFeatures(
      serial_features, *par_mesh, source_vertex_ids, source_element_ids);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  MaterialOperator mat_op(iodata, *meshes.back());
  BoundaryModeOperator mode_op(iodata, meshes, mat_op, &local_features, &source_vertex_ids);

  const auto nd_prolongations = mode_op.GetCombinedNDProlongationOperators();
  const auto h1_prolongations = mode_op.GetCombinedH1ProlongationOperators();
  const auto gradients = mode_op.GetCombinedGradientOperators();
  const auto &nd_spaces = mode_op.GetNDSpaceHierarchy();
  const auto &h1_spaces = mode_op.GetH1SpaceHierarchy();
  const auto &nd_essential = mode_op.GetCombinedNDDbcTDofLists();
  const auto &h1_essential = mode_op.GetCombinedH1DbcTDofLists();
  REQUIRE(nd_spaces.GetNumLevels() == 2);
  REQUIRE(h1_spaces.GetNumLevels() == nd_spaces.GetNumLevels());
  REQUIRE(nd_prolongations.size() + 1 == nd_spaces.GetNumLevels());
  REQUIRE(h1_prolongations.size() == nd_prolongations.size());
  REQUIRE(gradients.size() == nd_spaces.GetNumLevels());
  REQUIRE(nd_essential.size() == nd_spaces.GetNumLevels());
  REQUIRE(h1_essential.size() == h1_spaces.GetNumLevels());

  for (std::size_t level = 0; level < gradients.size(); level++)
  {
    const int standard_nd_size = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    const int standard_h1_size = h1_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    CHECK(gradients[level]->Height() - standard_nd_size ==
          static_cast<int>(mode_op.GetSingularDofNumbering().nd.owned_size));
    CHECK(gradients[level]->Width() - standard_h1_size ==
          static_cast<int>(mode_op.GetSingularDofNumbering().h1.owned_size));
    CHECK(std::is_sorted(nd_essential[level].begin(), nd_essential[level].end()));
    CHECK(std::is_sorted(h1_essential[level].begin(), h1_essential[level].end()));
  }

  const auto check_identity_enrichment = [](const mfem::HypreParMatrix &prolongation,
                                            int coarse_standard_size,
                                            int fine_standard_size)
  {
    mfem::Vector coarse(prolongation.Width()), fine(prolongation.Height());
    coarse = 0.0;
    for (int i = coarse_standard_size; i < coarse.Size(); i++)
    {
      coarse[i] = 0.13 * (i - coarse_standard_size + 1);
    }
    prolongation.Mult(coarse, fine);
    for (int i = 0; i < fine_standard_size; i++)
    {
      CHECK(fine[i] == 0.0);
    }
    for (int i = fine_standard_size; i < fine.Size(); i++)
    {
      CHECK(fine[i] == coarse[coarse_standard_size + i - fine_standard_size]);
    }
  };

  for (std::size_t level = 0; level < nd_prolongations.size(); level++)
  {
    const auto *nd_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations[level]);
    const auto *h1_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(h1_prolongations[level]);
    const auto *coarse_gradient =
        dynamic_cast<const mfem::HypreParMatrix *>(gradients[level]);
    const auto *fine_gradient =
        dynamic_cast<const mfem::HypreParMatrix *>(gradients[level + 1]);
    REQUIRE(nd_prolongation);
    REQUIRE(h1_prolongation);
    REQUIRE(coarse_gradient);
    REQUIRE(fine_gradient);
    CHECK(nd_prolongation->Width() == coarse_gradient->Height());
    CHECK(nd_prolongation->Height() == fine_gradient->Height());
    CHECK(h1_prolongation->Width() == coarse_gradient->Width());
    CHECK(h1_prolongation->Height() == fine_gradient->Width());

    check_identity_enrichment(*nd_prolongation,
                              nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize(),
                              nd_spaces.GetFESpaceAtLevel(level + 1).GetTrueVSize());
    check_identity_enrichment(*h1_prolongation,
                              h1_spaces.GetFESpaceAtLevel(level).GetTrueVSize(),
                              h1_spaces.GetFESpaceAtLevel(level + 1).GetTrueVSize());

    mfem::Vector coarse_h1(coarse_gradient->Width()), fine_h1(fine_gradient->Width()),
        fine_gradient_h1(fine_gradient->Height()),
        coarse_gradient_h1(coarse_gradient->Height()),
        prolonged_gradient(nd_prolongation->Height());
    FillVector(coarse_h1, 0.83);
    h1_prolongation->Mult(coarse_h1, fine_h1);
    fine_gradient->Mult(fine_h1, fine_gradient_h1);
    coarse_gradient->Mult(coarse_h1, coarse_gradient_h1);
    nd_prolongation->Mult(coarse_gradient_h1, prolonged_gradient);
    fine_gradient_h1 -= prolonged_gradient;
    CHECK(RelativeGlobalNorm(fine_gradient_h1, prolonged_gradient) < 2.0e-12);
  }

  const double omega = 2.0 * M_PI *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(
                           iodata.solver.boundary_mode.freq);
  auto finest_att = mode_op.AssembleAttPreconditioner(omega, -omega * omega);
  auto finest_ann = mode_op.AssembleAnnPreconditioner(omega);
  REQUIRE(finest_att);
  REQUIRE(finest_ann);
  finest_att->EliminateBC(nd_essential.back(), Operator::DIAG_ONE);
  finest_ann->EliminateBC(h1_essential.back(), Operator::DIAG_ONE);

  const auto check_galerkin = [](const mfem::HypreParMatrix &fine,
                                 const mfem::HypreParMatrix &prolongation,
                                 const mfem::Array<int> &coarse_essential)
  {
    std::unique_ptr<mfem::HypreParMatrix> coarse(
        mfem::RAP(&prolongation, &fine, &prolongation));
    REQUIRE(coarse);
    coarse->EliminateBC(coarse_essential, Operator::DIAG_ONE);
    mfem::Vector coarse_x(coarse->Width()), coarse_action(coarse->Height()),
        fine_x(fine.Width()), fine_action(fine.Height()),
        restricted_action(coarse->Height());
    FillVector(coarse_x, 0.41);
    linalg::SetSubVector(coarse_x, coarse_essential, 0.0);
    coarse->Mult(coarse_x, coarse_action);
    prolongation.Mult(coarse_x, fine_x);
    fine.Mult(fine_x, fine_action);
    prolongation.MultTranspose(fine_action, restricted_action);
    linalg::SetSubVector(restricted_action, coarse_essential, 0.0);
    restricted_action -= coarse_action;
    CHECK(RelativeGlobalNorm(restricted_action, coarse_action) < 2.0e-11);
  };

  const auto *nd_prolongation =
      dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations.front());
  const auto *h1_prolongation =
      dynamic_cast<const mfem::HypreParMatrix *>(h1_prolongations.front());
  REQUIRE(nd_prolongation);
  REQUIRE(h1_prolongation);
  check_galerkin(*finest_att, *nd_prolongation, nd_essential.front());
  check_galerkin(*finest_ann, *h1_prolongation, h1_essential.front());
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
