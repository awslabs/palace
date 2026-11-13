// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The interface dielectric postprocessing on interior metal sheets under both crack forms
// of the integration (SD-000): the legacy key (Model.CrackInternalBoundaryElements true:
// interior metal sheets cracked, every sheet face one-sided) and the new default (the key
// unset: PEC / Ground sheets uncracked, their faces shared by the elements on both sides).

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/laplaceoperator.hpp"
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

// The unit cube (4 x 4 x 4 tetrahedra) split by the plane x = 0.5: element attribute 1
// for x < 0.5 and 2 above; the interior faces of the plane with y and z in [0.25, 0.75]
// form the metal island (boundary attribute island_attr), the other faces of the plane the
// open interface (open_attr); mfem::Mesh::MakeCartesian3D numbers the exterior boundary 1
// to 6.
constexpr int island_attr = 7, open_attr = 8;

mfem::Mesh MakeSerialIslandMesh()
{
  mfem::Mesh mesh =
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
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
    bool on_plane = true, inside = true;
    for (const int v : vertices)
    {
      const double *x = mesh.GetVertex(v);
      on_plane = on_plane && std::abs(x[0] - 0.5) < 1.0e-12;
      inside = inside && x[1] > 0.25 - 1.0e-12 && x[1] < 0.75 + 1.0e-12 &&
               x[2] > 0.25 - 1.0e-12 && x[2] < 0.75 + 1.0e-12;
    }
    if (on_plane)
    {
      mesh.AddBdrElement(mesh.GetFace(face)->Duplicate(&mesh));
      mesh.SetBdrAttribute(mesh.GetNBE() - 1, inside ? island_attr : open_attr);
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
fs::path WriteIslandMesh(const fs::path &dir)
{
  const auto path = dir / "island.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream output(path);
    MakeSerialIslandMesh().Print(output);
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

}  // namespace palace
