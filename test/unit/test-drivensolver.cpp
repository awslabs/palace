// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/mesh.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"

using json = nlohmann::json;
using namespace palace;

TEST_CASE("PortOrthogonalityAdaptivePROM", "[driven_solver][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  // Recycle previous test input which has non-zero port overlap at edges.
  auto filename =
      fmt::format("{}/{}", PALACE_TEST_DIR, "postoperatorcsv_restart/restart.json");
  IoData iodata{filename.c_str(), false};

  iodata.model.mesh = fs::path(PALACE_TEST_DIR) / "mesh/fichera-tet.mesh";
  iodata.solver.driven.adaptive_tol = 1e-02;
  iodata.solver.driven.adaptive_circuit_synthesis = true;
  iodata.problem.output = "output/driven_adaptive/prom_ortho";

  MakeOutputFolder(iodata, world_comm);

  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::ReadMesh(iodata, world_comm));
    iodata.NondimensionalizeInputs(*mfem_mesh[0]);
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }

  DrivenSolver solver{iodata, Mpi::Root(world_comm), Mpi::Size(world_comm),
                      utils::ConfigureOmp(), ""};
  CHECK_THROWS(solver.SolveEstimateMarkRefine(mesh_),
               Catch::Matchers::ContainsSubstring(
                   "Lumped port modes should have exactly zero overlap"));
}
