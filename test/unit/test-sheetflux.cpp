// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The interface dielectric postprocessing on interior metal sheets under both crack forms
// of the integration (SD-000): the legacy key (Model.CrackInternalBoundaryElements true:
// interior metal sheets cracked, every sheet face one-sided) and the new default (the key
// unset: PEC / Ground sheets uncracked, their faces shared by the elements on both sides).

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/coefficient.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using json = nlohmann::json;
using namespace Catch::Matchers;
namespace fs = std::filesystem;

namespace
{

// The unit cube (4 x 4 x 4 tetrahedra or hexahedra) split by the plane x = 0.5: element
// attribute 1 for x < 0.5 and 2 above; the interior faces of the plane with y and z in
// [0.25, 0.75] form the metal island (boundary attribute island_attr), the other faces of
// the plane the open interface (open_attr), and the faces of the plane x = 0.75 facing the
// island the probe plane (probe_attr, vacuum on both sides); mfem::Mesh::MakeCartesian3D
// numbers the exterior boundary 1 to 6.
constexpr int island_attr = 7, open_attr = 8, probe_attr = 9;

mfem::Mesh MakeSerialIslandMesh(mfem::Element::Type type = mfem::Element::TETRAHEDRON)
{
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, type, 1.0, 1.0, 1.0);
  for (int face = 0; face < mesh.GetNumFaces(); face++)
  {
    int e1, e2;
    mesh.GetFaceElements(face, &e1, &e2);
    if (e1 < 0 || e2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    mesh.GetFaceVertices(face, vertices);
    bool on_plane = true, on_probe = true, inside = true;
    for (const int v : vertices)
    {
      const double *x = mesh.GetVertex(v);
      on_plane = on_plane && std::abs(x[0] - 0.5) < 1.0e-12;
      on_probe = on_probe && std::abs(x[0] - 0.75) < 1.0e-12;
      inside = inside && x[1] > 0.25 - 1.0e-12 && x[1] < 0.75 + 1.0e-12 &&
               x[2] > 0.25 - 1.0e-12 && x[2] < 0.75 + 1.0e-12;
    }
    if (on_plane || (on_probe && inside))
    {
      mesh.AddBdrElement(mesh.GetFace(face)->Duplicate(&mesh));
      mesh.SetBdrAttribute(mesh.GetNBE() - 1,
                           on_probe ? probe_attr : (inside ? island_attr : open_attr));
    }
  }
  mfem::Vector center(3);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, center);
    mesh.SetAttribute(e, (center(0) < 0.5) ? 1 : 2);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();
  return mesh;
}

// The island mesh written once (root) into the shared temporary directory.
fs::path WriteIslandMesh(const fs::path &dir,
                         mfem::Element::Type type = mfem::Element::TETRAHEDRON)
{
  const auto path = dir / ("island-" + std::to_string(int(type)) + ".mesh");
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream output(path);
    MakeSerialIslandMesh(type).Print(output);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

// An electrostatic configuration of the island mesh: the island is a Ground (PEC) sheet,
// the exterior face z = 0 a terminal, the permittivity of each side as given (1 = vacuum);
// one interface dielectric entry of the type given on the island, with automatic metal
// edges and the localized edge energy of the response matrix. Legacy states the legacy
// crack key, Default leaves it unset (the island is then not cracked).
json IslandConfig(const fs::path &mesh_path, const fs::path &output, bool legacy,
                  double eps_left, double eps_right, const std::string &type)
{
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Verbose", 0}, {"Output", output.string()}}},
      {"Model", {{"Mesh", mesh_path.string()}}},
      {"Domains",
       {{"Materials", json::array({{{"Attributes", {1}}, {"Permittivity", eps_left}},
                                   {{"Attributes", {2}}, {"Permittivity", eps_right}}})}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {island_attr}}}},
        {"Terminal", json::array({{{"Index", 1}, {"Attributes", {1}}}})},
        {"Postprocessing",
         {{"Dielectric", json::array({{{"Index", 1},
                                       {"Attributes", {island_attr}},
                                       {"Type", type},
                                       {"Thickness", 2.0e-3},
                                       {"Permittivity", 10.0},
                                       {"AutomaticEdges", true},
                                       {"EdgeDistances", {0.1, 0.2}},
                                       {"LocalizeEdgeEnergy", true}}})}}}}},
      {"Solver", {{"Order", 1}}}};
  if (legacy)
  {
    config["Model"]["CrackInternalBoundaryElements"] = true;
  }
  return config;
}

// A partition of the island mesh in the file format of Model.Partitioning, written once
// (root): the rank of an element = (x > 0.5) + 2 (y > 0.5) + (z > 0.5) of its center,
// modulo the process count, so that the process boundaries lie ON the sheet plane and
// ACROSS the island (the planes y = 0.5 and z = 0.5 cut it).
fs::path WriteIslandPartition(const fs::path &dir, mfem::Element::Type type)
{
  const int num_procs = Mpi::Size(Mpi::World());
  const auto path = dir / ("island-" + std::to_string(int(type)) + "-np" +
                           std::to_string(num_procs) + ".partition");
  if (Mpi::Root(Mpi::World()))
  {
    auto mesh = MakeSerialIslandMesh(type);
    std::ofstream output(path);
    output << "number_of_elements " << mesh.GetNE() << "\n"
           << "number_of_processors " << num_procs << "\n";
    mfem::Vector center(3);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementCenter(e, center);
      output << ((center(0) > 0.5) + 2 * (center(1) > 0.5) + (center(2) > 0.5)) % num_procs
             << "\n";
    }
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

// An electrostatic configuration of the island mesh with FluxRecovery entries NEXT to the
// island (never on it: our method rule keeps FluxRecovery off interior sheets, ruling
// 1600): 1 = SA on the open interface (the rest of the plane x = 0.5, which meets the
// island along its perimeter), 2 = MA on the exterior face x = 1 (vacuum side), 3 = MS on
// the exterior face x = 0 (substrate side), 4 = entry 1 without FluxRecovery, 5 = MA on the
// probe plane (it reads the normal flux of the elements that face the island). Substrate
// (permittivity 4) for x < 0.5, vacuum above; second order; a tight flux recovery; the
// crack without the bisection of RefineCrackElements in both forms (the same elements), and
// the partition.
json FluxRecoveryIslandConfig(const fs::path &mesh_path, const fs::path &partition,
                              const fs::path &output, bool legacy)
{
  auto Entry = [](int index, int attr, const char *type, bool flux_recovery)
  {
    return json{{"Index", index},       {"Attributes", {attr}},
                {"Type", type},         {"Thickness", 2.0e-3},
                {"Permittivity", 10.0}, {"FluxRecovery", flux_recovery}};
  };
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Verbose", 0}, {"Output", output.string()}}},
      {"Model",
       {{"Mesh", mesh_path.string()},
        {"Partitioning", partition.string()},
        {"RefineCrackElements", false}}},
      {"Domains",
       {{"Materials", json::array({{{"Attributes", {1}}, {"Permittivity", 4.0}},
                                   {{"Attributes", {2}}, {"Permittivity", 1.0}}})}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {island_attr}}}},
        {"Terminal", json::array({{{"Index", 1}, {"Attributes", {1}}}})},
        {"Postprocessing",
         {{"Dielectric",
           json::array({Entry(1, open_attr, "SA", true), Entry(2, 3, "MA", true),
                        Entry(3, 5, "MS", true), Entry(4, open_attr, "SA", false),
                        Entry(5, probe_attr, "MA", true)})}}}}},
      {"Solver",
       {{"Order", 2},
        {"Linear", {{"EstimatorTol", 1.0e-14}, {"EstimatorMaxIts", 10000}}}}}};
  if (legacy)
  {
    config["Model"]["CrackInternalBoundaryElements"] = true;
  }
  return config;
}

// The FluxRecovery energies of the entries 1-5 for a smooth field (no jump, so that the
// cracked and the uncracked mesh carry the same field element by element) with D recovered
// by the electrostatic error estimator as the driver builds it; on a broken flux space also
// the probe entry 5 with the CONTINUOUS recovery (the RT space, not broken across the
// island), i.e. what a D on the RT space gives there.
// The per-side D read AT the island faces, one-sided (MS reads the substrate side, MA the
// vacuum side; no sum-sides marker), directly through the coefficient. The
// SurfacePostOperator refuses FluxRecovery on an interior sheet (decision 1600, R-b), so
// its FluxRecovery energies never read an element's DOFs on the sheet faces: these reads
// are the ones that fail when the per-side holder hands both sides the same version of a
// split DOF (P3-B review MAJOR-1).
class RecoveredFluxPostOperator : public PostOperator<ProblemType::ELECTROSTATIC>
{
public:
  using PostOperator<ProblemType::ELECTROSTATIC>::PostOperator;
  const GridFunction *RecoveredD() const { return this->D_recovered.get(); }
};

template <InterfaceDielectric Type>
double IslandSideEnergy(const GridFunction &E, const GridFunction *D,
                        const MaterialOperator &mat_op)
{
  auto &mesh = *E.ParFESpace()->GetParMesh();
  mesh.ExchangeFaceNbrData();
  InterfaceDielectricCoefficient<Type> coeff(E, mat_op, 2.0e-3, 10.0, nullptr, D);
  double sum = 0.0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (mesh.GetBdrAttribute(be) != island_attr)
    {
      continue;
    }
    auto *T = mesh.GetBdrElementTransformation(be);
    const auto &ir = mfem::IntRules.Get(T->GetGeometryType(), 6);
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      sum += coeff.Eval(*T, ip) * ip.weight * T->Weight();
    }
  }
  Mpi::GlobalSum(1, &sum, mesh.GetComm());
  return sum;
}

struct IslandFluxEnergies
{
  double island_ms = 0.0, island_ma = 0.0;
  std::vector<double> energies;
  double continuous_probe = 0.0;
  bool broken = false;
  long long flux_true_dofs = 0;
};

IslandFluxEnergies RecoverIslandFluxEnergies(const json &config)
{
  IoData iodata(config, false);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
  LaplaceOperator laplace_op(iodata, mesh);
  RecoveredFluxPostOperator post_op(iodata, laplace_op);
  REQUIRE(post_op.NeedsRecoveredElectricFlux());
  mfem::FunctionCoefficient V_coeff(
      [](const mfem::Vector &x) { return 1.0 - x(2) + 0.3 * x(0) * x(1) - 0.2 * x(1); });
  mfem::VectorFunctionCoefficient E_coeff(3,
                                          [](const mfem::Vector &x, mfem::Vector &v)
                                          {
                                            v(0) = -0.3 * x(1);
                                            v(1) = -0.3 * x(0) + 0.2;
                                            v(2) = 1.0;
                                          });
  GridFunction V(laplace_op.GetH1Space()), E(laplace_op.GetNDSpace());
  V.Real().ProjectCoefficient(V_coeff);
  E.Real().ProjectCoefficient(E_coeff);
  Vector v(laplace_op.GetH1Space().GetTrueVSize()),
      e(laplace_op.GetNDSpace().GetTrueVSize());
  V.Real().GetTrueDofs(v);
  E.Real().GetTrueDofs(e);
  E.Real().SetFromTrueDofs(e);
  E.Real().ExchangeFaceNbrData();
  const auto &linear = iodata.solver.linear;
  auto Estimator = [&](const std::vector<int> &crack_attr_list)
  {
    return GradFluxErrorEstimator<Vector>(laplace_op.GetMaterialOp(),
                                          laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
                                          linear.estimator_tol, linear.estimator_max_it, 0,
                                          linear.estimator_mg, crack_attr_list);
  };
  const auto estimator = Estimator(iodata.boundaries.GetSheetAttributes());
  IslandFluxEnergies result;
  result.broken = estimator.GetFluxSpace().IsBroken();
  result.flux_true_dofs = estimator.GetFluxSpace().GlobalTrueVSize();
  REQUIRE(post_op.GetRecoveredElectricFluxTrueVSize() ==
          estimator.GetFluxSpace().GetTrueVSize());
  Vector d(estimator.GetFluxSpace().GetTrueVSize());
  estimator.RecoverFlux(e, d);
  const auto energies = post_op.GetElectrostaticEnergies(v, e, &d);
  for (int idx : {1, 2, 3, 4, 5})
  {
    result.energies.push_back(energies.interfaces.at(idx).energy);
  }
  result.island_ms = IslandSideEnergy<InterfaceDielectric::MS>(E, post_op.RecoveredD(),
                                                               laplace_op.GetMaterialOp());
  result.island_ma = IslandSideEnergy<InterfaceDielectric::MA>(E, post_op.RecoveredD(),
                                                               laplace_op.GetMaterialOp());
  if (result.broken)
  {
    const auto continuous = Estimator({});
    Vector dc(continuous.GetFluxSpace().GetTrueVSize());
    continuous.RecoverFlux(e, dc);
    GridFunction Dc(laplace_op.GetRTSpace());
    Dc.Real().SetFromTrueDofs(dc);
    Dc.Real().ExchangeFaceNbrData();
    result.continuous_probe =
        post_op.GetSurfacePostOp().GetInterfaceElectricFieldEnergy(5, E, &Dc);
  }
  return result;
}

}  // namespace

// P3-16: the aggregated response matrix (the library build, also the response archive's
// writer and reducer) and the archive estimate use batched interface amplitudes, which
// cannot represent the summed energies of both sides of an interior sheet face (SD-401
// (d)). Such a run is refused at the construction of its postprocessing, before any solve,
// where the island (uncracked under the default) has two qualifying sides: vacuum on both
// sides for MA, substrate on both sides for MS. One qualifying side, the cracked island of
// the legacy key (one-sided copies) and the non-aggregated response matrix (the per-pair
// energies, which sum the sides) are not refused.
TEST_CASE_METHOD(test::SharedTempDir,
                 "Batched interface response refused before the solve on a summed sheet",
                 "[sheetflux][Serial][Parallel]")
{
  const auto mesh_path = WriteIslandMesh(temp_dir);
  struct Island
  {
    double eps_left, eps_right;
    const char *type;
    bool two_qualifying_sides;
  };
  const std::vector<Island> islands = {{1.0, 1.0, "MA", true},  {1.0, 1.0, "MS", false},
                                       {4.0, 4.0, "MS", true},  {4.0, 4.0, "MA", false},
                                       {4.0, 1.0, "MA", false}, {4.0, 1.0, "MS", false}};
  const std::string refusal =
      "cannot represent the summed energies of both sides of an interior sheet";
  for (const bool legacy : {false, true})
  {
    for (const auto &island : islands)
    {
      INFO((legacy ? "legacy" : "default")
           << " crack form, " << island.type << ", permittivity " << island.eps_left
           << " | " << island.eps_right);
      IoData iodata(IslandConfig(mesh_path, temp_dir / "postpro", legacy, island.eps_left,
                                 island.eps_right, island.type),
                    false);
      std::vector<std::unique_ptr<Mesh>> mesh;
      mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
      REQUIRE(iodata.boundaries.cracked_attributes.contains(island_attr) == legacy);
      LaplaceOperator laplace_op(iodata, mesh);
      const bool refused = island.two_qualifying_sides && !legacy;

      // The archive estimate's call, whatever the response-matrix keys.
      {
        PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
        if (refused)
        {
          CHECK_THROWS_WITH(post_op.VerifyBatchedInterfaceResponse(),
                            ContainsSubstring(refusal));
        }
        else
        {
          CHECK_NOTHROW(post_op.VerifyBatchedInterfaceResponse());
        }
      }

      // The aggregated response matrix: refused at the construction.
      iodata.solver.electrostatic.response_matrix = true;
      iodata.solver.electrostatic.aggregate_response_matrix = true;
      if (refused)
      {
        CHECK_THROWS_WITH(PostOperator<ProblemType::ELECTROSTATIC>(iodata, laplace_op),
                          ContainsSubstring(refusal));
      }
      else
      {
        CHECK_NOTHROW(PostOperator<ProblemType::ELECTROSTATIC>(iodata, laplace_op));
      }

      // The per-pair response matrix evaluates the summed energies: not refused.
      iodata.solver.electrostatic.aggregate_response_matrix = false;
      CHECK_NOTHROW(PostOperator<ProblemType::ELECTROSTATIC>(iodata, laplace_op));
    }
  }
}

// P3-14 (SD-3202, T1): with FluxRecovery next to an interior sheet that is NOT cracked (the
// new default), the drivers' recovered D lives on the flux space broken across the sheet,
// and the postprocessing reads it per side: each element sees the version of its own side.
// That is the function space of the same elements CRACKED along the sheet (the legacy key,
// without the bisection), so the FluxRecovery energies of the two forms agree to the
// recovery's tolerance, with process boundaries on and across the sheet (the other side
// through the face-neighbour data); the continuous recovery (a D on the RT space, as
// before the broken space) differs next to the sheet.
TEST_CASE_METHOD(test::SharedTempDir,
                 "Per-side recovered flux on an uncracked sheet equals the cracked mesh's",
                 "[sheetflux][Serial][Parallel]")
{
  constexpr auto type = mfem::Element::HEXAHEDRON;
  const auto mesh_path = WriteIslandMesh(temp_dir, type);
  const auto partition = WriteIslandPartition(temp_dir, type);
  const auto uncracked = RecoverIslandFluxEnergies(
      FluxRecoveryIslandConfig(mesh_path, partition, temp_dir / "postpro-default", false));
  const auto cracked = RecoverIslandFluxEnergies(
      FluxRecoveryIslandConfig(mesh_path, partition, temp_dir / "postpro-legacy", true));
  CHECK(uncracked.broken);
  CHECK(!cracked.broken);
  // The broken view and the RT space of the cracked mesh have the same true DOFs.
  CHECK(uncracked.flux_true_dofs == cracked.flux_true_dofs);
  for (std::size_t i = 0; i < uncracked.energies.size(); i++)
  {
    INFO("entry " << i + 1 << ": uncracked (per side) " << std::setprecision(17)
                  << uncracked.energies[i] << ", cracked " << cracked.energies[i]
                  << ", relative difference "
                  << std::abs(uncracked.energies[i] - cracked.energies[i]) /
                         std::abs(cracked.energies[i]));
    CHECK(cracked.energies[i] > 0.0);
    // MEASURED (np 1 / 2 / 3): <= 1e-14 relative; frozen at 1e-12.
    CHECK_THAT(uncracked.energies[i], WithinRel(cracked.energies[i], 1.0e-12));
  }
  INFO("island MS per side " << std::setprecision(17) << uncracked.island_ms << " cracked "
                             << cracked.island_ms << "; island MA per side "
                             << uncracked.island_ma << " cracked " << cracked.island_ma);
  CHECK(cracked.island_ms > 0.0);
  CHECK(cracked.island_ma > 0.0);
  CHECK_THAT(uncracked.island_ms, WithinRel(cracked.island_ms, 1.0e-10));
  CHECK_THAT(uncracked.island_ma, WithinRel(cracked.island_ma, 1.0e-10));
  // The probe plane reads the normal flux of the elements facing the island: there the
  // per-side recovery differs from the continuous one. (The open interface does not show
  // it on this Cartesian hexahedral mesh: the normal flux of its faces is decoupled from
  // the island's in the projection, MEASURED bitwise equal.)
  INFO("entry 5 with the continuous recovery " << std::setprecision(17)
                                               << uncracked.continuous_probe);
  CHECK(std::abs(uncracked.continuous_probe - uncracked.energies[4]) >
        1.0e-6 * uncracked.energies[4]);
}

// P3-14 (SD-3202, T6): the response-matrix paths with FluxRecovery read the per-side D too:
// the batched Gram (the in-run response matrix) and the streamed rows (the archive
// reduction, whose FLUX vectors have the flux space's size, P3-21) agree with a D on the
// broken flux space, and the uncracked (per-side) matrices equal the cracked mesh's.
TEST_CASE_METHOD(
    test::SharedTempDir,
    "Batched and streamed interface responses read the per-side recovered flux",
    "[sheetflux][Serial][Parallel]")
{
  constexpr auto type = mfem::Element::HEXAHEDRON;
  const auto mesh_path = WriteIslandMesh(temp_dir, type);
  const auto partition = WriteIslandPartition(temp_dir, type);
  using ResponseMatrices =
      std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>>;
  auto Responses = [&](bool legacy)
  {
    auto config = FluxRecoveryIslandConfig(
        mesh_path, partition, temp_dir / (legacy ? "rm-legacy" : "rm-default"), legacy);
    auto &entry = config["Boundaries"]["Postprocessing"]["Dielectric"][0];
    entry["AutomaticEdges"] = true;
    entry["EdgeDistances"] = {0.1, 0.2};
    entry["LocalizeEdgeEnergy"] = true;
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
    LaplaceOperator laplace_op(iodata, mesh);
    PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
    const auto &linear = iodata.solver.linear;
    GradFluxErrorEstimator<Vector> estimator(
        laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
        linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg,
        iodata.boundaries.GetSheetAttributes());
    CHECK(estimator.GetFluxSpace().IsBroken() == !legacy);
    // Two basis fields and their recovered fluxes.
    std::vector<Vector> e(2), d(2);
    for (int i = 0; i < 2; i++)
    {
      mfem::VectorFunctionCoefficient f(3,
                                        [i](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = (i == 0) ? -0.3 * x(1) : 0.2 + 0.1 * x(2);
                                          v(1) = (i == 0) ? -0.3 * x(0) + 0.2 : 0.5 * x(0);
                                          v(2) = (i == 0) ? 1.0 : -0.4 + 0.2 * x(1);
                                        });
      GridFunction E(laplace_op.GetNDSpace());
      E.Real().ProjectCoefficient(f);
      e[i].SetSize(laplace_op.GetNDSpace().GetTrueVSize());
      E.Real().GetTrueDofs(e[i]);
      d[i].SetSize(estimator.GetFluxSpace().GetTrueVSize());
      estimator.RecoverFlux(e[i], d[i]);
    }
    auto batched = post_op.GetInterfaceElectricFieldEnergyMatrices(e, &d);
    ResponseMatrices streamed;
    for (const auto &samples : post_op.CacheInterfaceResponseSamples())
    {
      std::vector<double> rows(2 * samples.RowSize());
      for (int i = 0; i < 2; i++)
      {
        post_op.SetInterfaceResponseField(e[i], &d[i]);
        post_op.EvaluateInterfaceResponseRow(samples, rows.data() + i * samples.RowSize());
      }
      streamed.emplace(samples.interface_index,
                       post_op.AssembleInterfaceResponseMatrices(samples, rows.data(), 2));
    }
    return std::make_pair(std::move(batched), std::move(streamed));
  };
  auto Compare =
      [](const ResponseMatrices &a, const ResponseMatrices &b, double tol, const char *what)
  {
    REQUIRE(a.size() == 1);
    REQUIRE(b.size() == 1);
    const auto &ma = a.at(1), &mb = b.at(1);
    REQUIRE(ma.size() == mb.size());
    REQUIRE(!ma.empty());
    for (std::size_t k = 0; k < ma.size(); k++)
    {
      for (int i = 0; i < 2; i++)
      {
        for (int j = 0; j < 2; j++)
        {
          INFO(what << ": radius " << k << ", (" << i << ", " << j << "): total "
                    << std::setprecision(17) << ma[k].energy_total(i, j) << " vs "
                    << mb[k].energy_total(i, j) << " (relative difference "
                    << std::abs(ma[k].energy_total(i, j) - mb[k].energy_total(i, j)) /
                           std::abs(mb[k].energy_total(i, j))
                    << "); inside " << ma[k].energy_inside(i, j) << " vs "
                    << mb[k].energy_inside(i, j));
          CHECK_THAT(ma[k].energy_total(i, j),
                     WithinRel(mb[k].energy_total(i, j), tol) ||
                         WithinAbs(mb[k].energy_total(i, j), 1.0e-30));
          CHECK_THAT(ma[k].energy_inside(i, j),
                     WithinRel(mb[k].energy_inside(i, j), tol) ||
                         WithinAbs(mb[k].energy_inside(i, j), 1.0e-30));
        }
      }
    }
  };
  const auto [batched_default, streamed_default] = Responses(false);
  const auto [batched_legacy, streamed_legacy] = Responses(true);
  // MEASURED (np 1 / 2 / 3): <= 1.4e-15 relative (all three comparisons); frozen at 1e-12.
  Compare(batched_default, streamed_default, 1.0e-12, "default: batched vs streamed");
  Compare(batched_legacy, streamed_legacy, 1.0e-12, "legacy: batched vs streamed");
  Compare(batched_default, batched_legacy, 1.0e-12,
          "per side (default) vs cracked (legacy)");
}

// P3-14 (SD-3202, T1; P3-B review MAJOR-1): tetrahedra, Order 2, the WHOLE plane x = 0.5 a
// Ground sheet (uncracked, the new default), a uniform E = (1, 0, 0): the recovered D on
// the flux space broken across the sheet is exactly 4 E below and E above; the per-side D
// read at the island faces one-sided (MS substrate side, MA vacuum side) must equal the
// read of the exact constant fields (4, 0, 0) / (1, 0, 0) on the RT space.
TEST_CASE_METHOD(
    test::SharedTempDir,
    "Per-side recovered flux read one-sided at an uncracked sheet on tetrahedra",
    "[sheetflux][Serial][Parallel]")
{
  const auto mesh_path = WriteIslandMesh(temp_dir, mfem::Element::TETRAHEDRON);
  json config = {
      {"Problem",
       {{"Type", "Electrostatic"},
        {"Verbose", 0},
        {"Output", (temp_dir / "postpro-tet2").string()}}},
      {"Model", {{"Mesh", mesh_path.string()}}},
      {"Domains",
       {{"Materials", json::array({{{"Attributes", {1}}, {"Permittivity", 4.0}},
                                   {{"Attributes", {2}}, {"Permittivity", 1.0}}})}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {island_attr, open_attr}}}},
        {"Terminal", json::array({{{"Index", 1}, {"Attributes", {1}}}})},
        {"Postprocessing",
         {{"Dielectric", json::array({{{"Index", 1},
                                       {"Attributes", {3}},
                                       {"Type", "MA"},
                                       {"Thickness", 2.0e-3},
                                       {"Permittivity", 10.0},
                                       {"FluxRecovery", true}}})}}}}},
      {"Solver",
       {{"Order", 2},
        {"Linear", {{"EstimatorTol", 1.0e-14}, {"EstimatorMaxIts", 10000}}}}}};
  IoData iodata(config, false);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
  LaplaceOperator laplace_op(iodata, mesh);
  RecoveredFluxPostOperator post_op(iodata, laplace_op);
  REQUIRE(post_op.NeedsRecoveredElectricFlux());
  auto X = [](double dx)
  {
    mfem::Vector v(3);
    v = 0.0;
    v(0) = dx;
    return v;
  };
  mfem::VectorConstantCoefficient E_coeff(X(1.0));
  GridFunction E(laplace_op.GetNDSpace());
  E.Real().ProjectCoefficient(E_coeff);
  Vector e(laplace_op.GetNDSpace().GetTrueVSize());
  E.Real().GetTrueDofs(e);
  E.Real().SetFromTrueDofs(e);
  E.Real().ExchangeFaceNbrData();
  const auto &linear = iodata.solver.linear;
  GradFluxErrorEstimator<Vector> estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg,
      iodata.boundaries.GetSheetAttributes());
  REQUIRE(estimator.GetFluxSpace().IsBroken());
  Vector d(estimator.GetFluxSpace().GetTrueVSize());
  estimator.RecoverFlux(e, d);
  post_op.SetRecoveredElectricFlux(d);
  auto Exact = [&](double dx)
  {
    auto D = std::make_unique<GridFunction>(laplace_op.GetRTSpace());
    mfem::VectorConstantCoefficient c(X(dx));
    D->Real().ProjectCoefficient(c);
    Vector t(laplace_op.GetRTSpace().GetTrueVSize());
    D->Real().GetTrueDofs(t);
    D->Real().SetFromTrueDofs(t);
    D->Real().ExchangeFaceNbrData();
    return D;
  };
  const auto D_sub = Exact(4.0), D_vac = Exact(1.0);
  const auto &mat_op = laplace_op.GetMaterialOp();
  const double ms =
      IslandSideEnergy<InterfaceDielectric::MS>(E, post_op.RecoveredD(), mat_op);
  const double ms_exact = IslandSideEnergy<InterfaceDielectric::MS>(E, D_sub.get(), mat_op);
  const double ma =
      IslandSideEnergy<InterfaceDielectric::MA>(E, post_op.RecoveredD(), mat_op);
  const double ma_exact = IslandSideEnergy<InterfaceDielectric::MA>(E, D_vac.get(), mat_op);
  INFO("tetrahedra MS per side " << std::setprecision(17) << ms << " exact " << ms_exact
                                 << "; MA per side " << ma << " exact " << ma_exact);
  CHECK(ms_exact > 0.0);
  CHECK(ma_exact > 0.0);
  CHECK_THAT(ms, WithinRel(ms_exact, 1.0e-8));
  CHECK_THAT(ma, WithinRel(ma_exact, 1.0e-8));
}

}  // namespace palace
