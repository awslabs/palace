// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// A PrescribedPotential boundary (a coupon potential trace and its one-volt conductor) is
// an essential (Dirichlet) conductor sheet, like a Terminal: when it lies in the interior
// of the domain, the fields are discontinuous across it but it is not cracked by default.
// It is therefore a sheet attribute (BoundaryData::GetSheetAttributes: broken flux recovery
// of the error estimator, two-sided interface dielectric energies) and is not in the
// default crack set (BoundaryData::GetMeshCrackAttributes), only in the legacy one.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

// The unit cube mesh of mfem::Mesh::MakeCartesian3D has the exterior boundary attributes
// 1 to 6; the interior sheet on the plane x = 0.5 gets the next one.
constexpr int sheet_attr = 7;
constexpr double x_sheet = 0.5;

// Unit cube mesh with interior boundary elements of attribute sheet_attr on the plane
// x = x_sheet (the mesh is not cracked: the faces are shared by the elements on both
// sides).
mfem::Mesh MakeInteriorSheetCubeMesh(int n, mfem::Element::Type type)
{
  auto base = mfem::Mesh::MakeCartesian3D(n, n, n, type, 1.0, 1.0, 1.0);
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    if (e2 < 0)
    {
      continue;
    }
    base.GetFaceVertices(f, fv);
    bool on_plane = true;
    for (auto v : fv)
    {
      on_plane = on_plane && (std::abs(base.GetVertex(v)[0] - x_sheet) < 1.0e-12);
    }
    if (on_plane)
    {
      sheet_faces.push_back(f);
    }
  }
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&mesh));
  }
  for (auto f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(sheet_attr);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.SetAttribute(e, 1);
  }
  mesh.SetAttributes();
  return mesh;
}

// Partition such that process boundaries both coincide with the sheet and cross it.
Mesh MakeParMesh(MPI_Comm comm, mfem::Mesh &smesh)
{
  const int num_procs = Mpi::Size(comm);
  std::vector<int> part(smesh.GetNE());
  mfem::Vector c(3);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    smesh.GetElementCenter(e, c);
    part[e] = ((c(0) > x_sheet) + 2 * (c(1) > 0.5) + 4 * (c(2) > 0.5)) % num_procs;
  }
  smesh.EnsureNodes();
  return Mesh(std::make_unique<mfem::ParMesh>(comm, smesh, part.data()));
}

// The boundaries of a coupon with a PrescribedPotential sheet on sheet_attr: either the
// trace itself (Attributes) with the exterior boundary grounded, or its one-volt conductor
// (TerminalAttributes) with the trace on the exterior boundary. The data file is only named
// (it is read by the Laplace operator, not by the configuration).
json PrescribedSheetBoundaries(bool sheet_is_terminal)
{
  const std::vector<int> exterior = {1, 2, 3, 4, 5, 6};
  if (sheet_is_terminal)
  {
    return {{"PrescribedPotential",
             {{{"Index", 1},
               {"Attributes", exterior},
               {"DataFile", "trace.csv"},
               {"TerminalAttributes", {sheet_attr}}}}}};
  }
  return {{"Ground", {{"Attributes", exterior}}},
          {"PrescribedPotential",
           {{{"Index", 1}, {"Attributes", {sheet_attr}}, {"DataFile", "trace.csv"}}}}};
}

mfem::Array<int> AttrToMarker(int max_attr, const std::vector<int> &attr_list)
{
  mfem::Array<int> marker(max_attr);
  marker = 0;
  for (auto attr : attr_list)
  {
    if (attr <= max_attr)
    {
      marker[attr - 1] = 1;
    }
  }
  return marker;
}

// The electric field of a surface charge on the sheet: a normal component which flips sign
// across it, E = -grad |x - x_sheet| (piecewise constant: recovered exactly by a broken
// flux space, not by a continuous one), and the same field without the jump.
void JumpE(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = (x(0) > x_sheet) ? -1.0 : 1.0;
  E(1) = 0.0;
  E(2) = 0.0;
}

void ConstantE(const mfem::Vector &x, mfem::Vector &E)
{
  E(0) = 1.0;
  E(1) = 0.0;
  E(2) = 0.0;
}

// An electrostatic coupon read from its configuration, as the electrostatic driver sets it
// up: the mesh loaded by mesh::ReadMesh with the crack form of the configuration, the
// Laplace operator, and the error estimator of ElectrostaticSolver (its recovery space
// broken across the sheet attributes of the configuration where the mesh has interior faces
// there).
struct PrescribedSheetCoupon
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> meshes;
  std::unique_ptr<LaplaceOperator> laplace_op;
  std::unique_ptr<GradFluxErrorEstimator<Vector>> estimator;

  explicit PrescribedSheetCoupon(const json &config) : iodata(config, false)
  {
    meshes.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
    laplace_op = std::make_unique<LaplaceOperator>(iodata, meshes);
    const auto &linear = iodata.solver.linear;
    estimator = std::make_unique<GradFluxErrorEstimator<Vector>>(
        laplace_op->GetMaterialOp(), laplace_op->GetNDSpace(), laplace_op->GetRTSpaces(),
        linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg,
        iodata.boundaries.GetSheetAttributes());
  }

  // The largest error indicator of the field that jumps across the sheet.
  double JumpIndicator(const GradFluxErrorEstimator<Vector> &jump_estimator) const
  {
    mfem::VectorFunctionCoefficient coeff(3, JumpE);
    mfem::ParGridFunction gf(&laplace_op->GetNDSpace().Get());
    gf.ProjectCoefficient(coeff);
    Vector E(laplace_op->GetNDSpace().GetTrueVSize());
    gf.ParallelProject(E);
    ErrorIndicator indicator;
    jump_estimator.AddErrorIndicator(E, 0.5, indicator);
    return indicator.Max(Mpi::World());
  }

  // The interface dielectric energies of the jump field (SurfacePostOperator from the
  // configuration, as PostOperator builds it), by entry index.
  std::vector<double> JumpInterfaceEnergies(const std::vector<int> &indices) const
  {
    SurfacePostOperator surf_post_op(iodata, laplace_op->GetMaterialOp(),
                                     laplace_op->GetH1Space(), laplace_op->GetNDSpace());
    GridFunction E(laplace_op->GetNDSpace());
    mfem::VectorFunctionCoefficient f(3, JumpE);
    E.Real().ProjectCoefficient(f);
    Vector e(laplace_op->GetNDSpace().GetTrueVSize());
    E.Real().GetTrueDofs(e);
    E.Real().SetFromTrueDofs(e);
    E.Real().ExchangeFaceNbrData();
    std::vector<double> energies;
    for (const int idx : indices)
    {
      energies.push_back(surf_post_op.GetInterfaceElectricFieldEnergy(idx, E, 0.0));
    }
    return energies;
  }
};

}  // namespace

TEST_CASE("PrescribedPotential sheet classification",
          "[prescribedpotential][config][Serial]")
{
  // The trace attributes and the terminal attributes of a PrescribedPotential entry are
  // sheet attributes, as Ground (PEC) and Terminal; the default crack set has only the
  // boundary conditions which apply to each side separately (ZeroCharge here); the
  // deprecated legacy option cracks every boundary condition attribute (but lumped ports).
  json boundaries = {{"Ground", {{"Attributes", {1}}}},
                     {"ZeroCharge", {{"Attributes", {2}}}},
                     {"PrescribedPotential",
                      {{{"Index", 1},
                        {"Attributes", {3}},
                        {"DataFile", "trace.csv"},
                        {"TerminalAttributes", {4}}},
                       {{"Index", 2}, {"Attributes", {5}}, {"DataFile", "trace2.csv"}}}}};
  config::BoundaryData data(boundaries);
  REQUIRE(data.prescribed_potential.size() == 2);
  CHECK(data.GetSheetAttributes() == std::vector<int>{1, 3, 4, 5});
  CHECK(data.GetMeshCrackAttributes(std::nullopt) == std::vector<int>{2});
  CHECK(data.GetMeshCrackAttributes(true) == std::vector<int>{1, 2, 3, 4, 5});
  CHECK(data.GetMeshCrackAttributes(false).empty());

  // A Terminal sheet is classified the same way (it is not a legacy crack attribute,
  // unlike a PrescribedPotential one, which the legacy option cracked before the default
  // changed).
  json terminal_boundaries = {{"Ground", {{"Attributes", {1}}}},
                              {"ZeroCharge", {{"Attributes", {2}}}},
                              {"Terminal", {{{"Index", 1}, {"Attributes", {3}}}}}};
  config::BoundaryData terminal_data(terminal_boundaries);
  CHECK(terminal_data.GetSheetAttributes() == std::vector<int>{1, 3});
  CHECK(terminal_data.GetMeshCrackAttributes(std::nullopt) == std::vector<int>{2});
  CHECK(terminal_data.GetMeshCrackAttributes(true) == std::vector<int>{1, 2});
}

TEST_CASE("PrescribedPotential sheet on an uncracked mesh",
          "[prescribedpotential][errorestimator][surfacefunctional][Serial][Parallel]")
{
  // On an uncracked mesh, an interior PrescribedPotential sheet (trace or terminal) gets
  // the broken flux recovery of the error estimator and the two-sided interface dielectric
  // energy, through the sheet attributes of the configuration.
  const auto comm = MPI_COMM_WORLD;
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = GENERATE(1, 2);
  const bool sheet_is_terminal = GENERATE(false, true);
  CAPTURE(static_cast<int>(type), order, sheet_is_terminal);
  fem::DefaultIntegrationOrder::p_trial = order;

  const config::BoundaryData boundaries(PrescribedSheetBoundaries(sheet_is_terminal));
  const auto sheet_attr_list = boundaries.GetSheetAttributes();
  REQUIRE(std::find(sheet_attr_list.begin(), sheet_attr_list.end(), sheet_attr) !=
          sheet_attr_list.end());
  const auto crack_attr_list = boundaries.GetMeshCrackAttributes(std::nullopt);
  REQUIRE(std::find(crack_attr_list.begin(), crack_attr_list.end(), sheet_attr) ==
          crack_attr_list.end());
  const std::vector<int> no_attr_list = {};

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 2.0, 2.0};
  config::PeriodicBoundaryData periodic;
  auto smesh = MakeInteriorSheetCubeMesh(4, type);
  auto mesh = MakeParMesh(comm, smesh);
  MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));

  SECTION("Broken flux recovery")
  {
    // The flux of the jump field is recovered exactly on the space broken across the
    // sheet attributes (indicator at round-off), not on the continuous space.
    auto Estimate = [&](const std::vector<int> &attr_list)
    {
      constexpr double tol = 1.0e-14;
      constexpr int max_it = 10000, print = 0;
      GradFluxErrorEstimator<Vector> estimator(mat_op, nd_fespace, rt_fespaces, tol, max_it,
                                               print, false, attr_list);
      mfem::VectorFunctionCoefficient coeff(3, JumpE);
      mfem::ParGridFunction gf(&nd_fespace.Get());
      gf.ProjectCoefficient(coeff);
      Vector E(nd_fespace.GetTrueVSize());
      gf.ParallelProject(E);
      ErrorIndicator indicator;
      estimator.AddErrorIndicator(E, 0.5, indicator);
      return indicator.Max(comm);
    };
    const double broken = Estimate(sheet_attr_list), continuous = Estimate(no_attr_list);
    CAPTURE(broken, continuous);
    CHECK(continuous > 1.0e-2);
    CHECK(broken < 1.0e-6 * continuous);
  }

  SECTION("Two-sided interface dielectric energy")
  {
    // On the sheet, the interface energy of the jump field sums both sides: it equals that
    // of the field without the jump (the same |E| on each side), which is twice the
    // one-sided value; the average of the two sides (no sum) cancels the normal field.
    auto Project = [&](void (*F)(const mfem::Vector &, mfem::Vector &))
    {
      auto E = std::make_unique<GridFunction>(nd_fespace);
      mfem::VectorFunctionCoefficient f(3, F);
      E->Real().ProjectCoefficient(f);
      return E;
    };
    auto E_jump = Project(JumpE);
    auto E_const = Project(ConstantE);
    const int bdr_attr_max = mesh.Get().bdr_attributes.Max();
    const auto marker = AttrToMarker(bdr_attr_max, {sheet_attr});
    const auto sum_sides_marker = AttrToMarker(bdr_attr_max, sheet_attr_list);
    const double t_i = 2.0e-3, epsilon_i = 10.0;
    SurfaceFunctional sum_sides(mesh, marker, nd_fespace.Get(), mat_op,
                                InterfaceDielectric::DEFAULT, t_i, epsilon_i,
                                &sum_sides_marker);
    SurfaceFunctional average(mesh, marker, nd_fespace.Get(), mat_op,
                              InterfaceDielectric::DEFAULT, t_i, epsilon_i, nullptr);
    REQUIRE(sum_sides.IsValid());
    REQUIRE(average.IsValid());
    const double sum_jump = sum_sides.Eval(*E_jump, 0.0),
                 sum_const = sum_sides.Eval(*E_const, 0.0),
                 avg_jump = average.Eval(*E_jump, 0.0),
                 avg_const = average.Eval(*E_const, 0.0);
    CAPTURE(sum_jump, sum_const, avg_jump, avg_const);
    CHECK(avg_const > 0.0);
    CHECK_THAT(sum_const, WithinRel(2.0 * avg_const, 1.0e-10));
    CHECK_THAT(sum_jump, WithinRel(sum_const, 1.0e-10));
    CHECK(std::abs(avg_jump) < 1.0e-10 * avg_const);
  }
}

TEST_CASE("PrescribedPotential sheet end to end in both crack forms",
          "[prescribedpotential][errorestimator][surfacepostoperator][Serial][Parallel]")
{
  // P3-31 (SL2 review m-4 (b); SD-201; PLAN-ERRATUM-1560 E-1): an electrostatic coupon with
  // an interior PrescribedPotential sheet (its trace, or its one-volt conductor), read from
  // its configuration under the legacy key (the sheet is cracked) and the new default (not
  // cracked), through the electrostatic driver's wiring: the driver's error estimator
  // recovers the flux of a field that jumps across the sheet exactly in both forms (a
  // broken recovery space on the uncracked mesh, the crack on the legacy one), and the
  // interface dielectric energies of SurfacePostOperator built from the configuration sum
  // both sides of the uncracked sheet, equal to the per-copy sum on the cracked one.
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const bool sheet_is_terminal = GENERATE(false, true);
  CAPTURE(static_cast<int>(type), sheet_is_terminal);
  test::SharedTempDir temp;
  const auto mesh_path = temp.temp_dir / "sheet-cube.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    auto smesh = MakeInteriorSheetCubeMesh(4, type);
    std::ofstream output(mesh_path);
    smesh.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  json boundaries = PrescribedSheetBoundaries(sheet_is_terminal);
  boundaries["Postprocessing"] = {{"Dielectric", json::array({{{"Index", 1},
                                                               {"Attributes", {sheet_attr}},
                                                               {"Type", "Default"},
                                                               {"Thickness", 2.0e-3},
                                                               {"Permittivity", 10.0}},
                                                              {{"Index", 2},
                                                               {"Attributes", {sheet_attr}},
                                                               {"Type", "MS"},
                                                               {"Thickness", 2.0e-3},
                                                               {"Permittivity", 10.0}}})}};
  auto Config = [&](bool legacy_key)
  {
    json config = {
        {"Problem",
         {{"Type", "Electrostatic"}, {"Verbose", 0}, {"Output", temp.temp_dir.string()}}},
        {"Model", {{"Mesh", mesh_path.string()}}},
        {"Domains",
         {{"Materials", json::array({{{"Attributes", {1}}, {"Permittivity", 2.0}}})}}},
        {"Boundaries", boundaries},
        {"Solver",
         {{"Order", 1},
          {"Linear", {{"EstimatorTol", 1.0e-14}, {"EstimatorMaxIts", 10000}}}}}};
    if (legacy_key)
    {
      config["Model"]["CrackInternalBoundaryElements"] = true;
    }
    return config;
  };
  PrescribedSheetCoupon legacy(Config(true)), uncracked(Config(false));
  REQUIRE(legacy.iodata.boundaries.cracked_attributes.contains(sheet_attr));
  REQUIRE(!uncracked.iodata.boundaries.cracked_attributes.contains(sheet_attr));

  // The recovery space is broken across the uncracked sheet only; the flux of the jump
  // field is recovered in both forms (indicator at round-off), not on a continuous space.
  CHECK(uncracked.estimator->GetFluxSpace().GlobalTrueVSize() >
        uncracked.laplace_op->GetRTSpace().GlobalTrueVSize());
  CHECK(legacy.estimator->GetFluxSpace().GlobalTrueVSize() ==
        legacy.laplace_op->GetRTSpace().GlobalTrueVSize());
  const auto &linear = uncracked.iodata.solver.linear;
  GradFluxErrorEstimator<Vector> continuous(
      uncracked.laplace_op->GetMaterialOp(), uncracked.laplace_op->GetNDSpace(),
      uncracked.laplace_op->GetRTSpaces(), linear.estimator_tol, linear.estimator_max_it, 0,
      linear.estimator_mg);
  const double reference = uncracked.JumpIndicator(continuous);
  const double broken = uncracked.JumpIndicator(*uncracked.estimator);
  const double cracked = legacy.JumpIndicator(*legacy.estimator);
  CAPTURE(reference, broken, cracked);
  CHECK(reference > 1.0e-2);
  CHECK(broken < 1.0e-6 * reference);
  CHECK(cracked < 1.0e-6 * reference);

  // The two-sided interface energies (Default and MS: the permittivity 2 on both sides
  // makes both substrate sides) equal the per-copy sums of the cracked sheet.
  const auto legacy_energies = legacy.JumpInterfaceEnergies({1, 2});
  const auto uncracked_energies = uncracked.JumpInterfaceEnergies({1, 2});
  for (std::size_t i = 0; i < legacy_energies.size(); i++)
  {
    INFO("interface " << i + 1);
    CHECK(legacy_energies[i] > 0.0);
    CHECK_THAT(uncracked_energies[i], WithinRel(legacy_energies[i], 1.0e-10));
  }
}

}  // namespace palace
