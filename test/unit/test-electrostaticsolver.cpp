// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <memory>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "drivers/electrostaticsolver.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using json = nlohmann::json;

namespace
{

// Load, nondimensionalize, partition, and wrap the mesh referenced by the config. Mirrors
// LoadScaleParMesh in test-waveportoperator.cpp.
auto LoadScaleParMesh(IoData &iodata, MPI_Comm comm)
{
  std::vector<std::unique_ptr<Mesh>> mesh;
  auto smesh = mesh::Load(iodata, comm);
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
  }
  iodata.NondimensionalizeInputs(smesh);
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
  mesh::RefineMesh(iodata, mfem_mesh);
  for (auto &m : mfem_mesh)
  {
    mesh.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  return mesh;
}

// Read the cavity2d electrostatic regression config, absolutizing the mesh path so it
// resolves regardless of cwd. The solver is constructed as a non-root instance below to
// suppress metadata output.
json LoadCavity2dElectrostaticConfig()
{
  auto dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "cavity2d";
  std::ifstream f(dir / "cavity2d_electrostatic.json");
  REQUIRE(f.good());
  json setup = json::parse(f, /*cb=*/nullptr, /*allow_exceptions=*/true,
                           /*ignore_comments=*/true);
  setup["Model"]["Mesh"] = (dir / setup["Model"]["Mesh"].get<std::string>()).string();
  setup["Problem"]["Output"] = "";
  return setup;
}

}  // namespace

// Exercise the solution-exposing inner Solve overload (added for MMS): build the operator
// here, call Solve(V, laplace_op) via the test-only subclass, and check the returned field.
TEST_CASE("ElectrostaticSolver exposes per-terminal potentials",
          "[electrostaticsolver][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(LoadCavity2dElectrostaticConfig(), /*print=*/false);
  auto mesh = LoadScaleParMesh(iodata, comm);

  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  LaplaceOperator laplace_op(iodata, mesh);

  std::vector<Vector> V;
  solver.Solve(V, laplace_op);

  // The contract: one solved potential per terminal.
  const std::size_t n_terminal = laplace_op.GetSources().size();
  REQUIRE(n_terminal > 0);
  REQUIRE(V.size() == n_terminal);

  // Each potential should be a genuine (nonzero) solution over the H1 space.
  REQUIRE(laplace_op.GlobalTrueVSize() > 0);
  for (const auto &Vi : V)
  {
    CHECK(Vi.Size() > 0);
    CHECK(linalg::Norml2(comm, Vi) > 0.0);
  }
}
