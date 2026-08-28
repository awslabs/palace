// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <memory>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;

namespace
{

// Helper to create a file with given content.
void CreateFile(const fs::path &path, const std::string &content = "test")
{
  std::ofstream f(path);
  f << content;
}

// Helper to read file content.
std::string ReadFile(const fs::path &path)
{
  std::ifstream f(path);
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// Minimal concrete BaseSolver so BaseSolver::SaveMetadata can be exercised directly; the
// Solve override is never invoked by these tests.
class TestSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;
  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &) const override
  {
    return {};
  }
};

// Verify the emitted true entity counts reproduce the EXACT conforming DOF count of the
// H1 and ND spaces (used by Palace's solvers) across orders, via the standard per-entity
// decomposition on tetrahedra. This is what lets the counts size a re-run at any order --
// a stronger check than true < leaf (it also catches an edge/face swap).
void CheckDofDecomposition(mfem::Mesh &mesh, const mesh::MeshEntityCounts &counts)
{
  const long long V = counts.true_vertices, E = counts.true_edges, F = counts.true_faces,
                  T = counts.elements;
  const int dim = mesh.Dimension();
  // H1 order p on tets: V + E(p-1) + F(p-1)(p-2)/2 + T(p-1)(p-2)(p-3)/6. Order 4 exercises
  // the tetrahedron-interior term.
  for (int p = 1; p <= 4; p++)
  {
    const long long expected = V + E * (p - 1) + F * ((p - 1) * (p - 2) / 2) +
                               T * ((p - 1) * (p - 2) * (p - 3) / 6);
    mfem::H1_FECollection fec(p, dim);
    CHECK(mfem::FiniteElementSpace(&mesh, &fec).GetTrueVSize() == expected);
  }
  // ND (Nedelec, first kind) order p on tets: E*p + F*p(p-1) + T*p(p-1)(p-2)/2. Order 3
  // exercises the tetrahedron-interior term.
  for (int p = 1; p <= 3; p++)
  {
    const long long expected = E * p + F * (p * (p - 1)) + T * (p * (p - 1) * (p - 2) / 2);
    mfem::ND_FECollection fec(p, dim);
    CHECK(mfem::FiniteElementSpace(&mesh, &fec).GetTrueVSize() == expected);
  }
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration moves files and leaves symlinks", "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Create test files in the output directory.
  CreateFile(temp_dir / "palace.json", R"({"key": "value"})");
  CreateFile(temp_dir / "domain-E.csv", "f,E\n1.0,2.0\n");
  CreateFile(temp_dir / "port-S.csv", "f,S11\n1.0,0.5\n");

  SaveIteration(comm, temp_dir, 1, 2);

  auto iter_dir = temp_dir / "iteration01";

  SECTION("Iteration subfolder is created")
  {
    CHECK(fs::is_directory(iter_dir));
  }

  SECTION("Regular files are moved to iteration subfolder")
  {
    CHECK(fs::is_regular_file(iter_dir / "domain-E.csv"));
    CHECK(ReadFile(iter_dir / "domain-E.csv") == "f,E\n1.0,2.0\n");
    CHECK(fs::is_regular_file(iter_dir / "port-S.csv"));
    CHECK(ReadFile(iter_dir / "port-S.csv") == "f,S11\n1.0,0.5\n");
  }

  SECTION("Symlinks are left behind for moved files")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(fs::is_symlink(temp_dir / "port-S.csv"));
    // Symlinks should resolve to the moved files.
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "f,E\n1.0,2.0\n");
    CHECK(ReadFile(temp_dir / "port-S.csv") == "f,S11\n1.0,0.5\n");
  }

  SECTION("Symlinks are relative")
  {
    auto target = fs::read_symlink(temp_dir / "domain-E.csv");
    CHECK(target.is_relative());
    CHECK(target == fs::path("iteration01") / "domain-E.csv");
  }

  SECTION("palace.json is copied, not moved or symlinked")
  {
    CHECK(fs::is_regular_file(temp_dir / "palace.json"));
    CHECK(!fs::is_symlink(temp_dir / "palace.json"));
    CHECK(fs::is_regular_file(iter_dir / "palace.json"));
    CHECK(ReadFile(temp_dir / "palace.json") == R"({"key": "value"})");
    CHECK(ReadFile(iter_dir / "palace.json") == R"({"key": "value"})");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SaveIteration handles directories",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  CreateFile(temp_dir / "palace.json", "{}");
  fs::create_directories(temp_dir / "paraview" / "driven");
  CreateFile(temp_dir / "paraview" / "driven" / "fields.vtu", "<VTK/>");

  SaveIteration(comm, temp_dir, 1, 1);

  auto iter_dir = temp_dir / "iteration1";

  SECTION("Directory is moved to iteration subfolder")
  {
    CHECK(fs::is_directory(iter_dir / "paraview" / "driven"));
    CHECK(ReadFile(iter_dir / "paraview" / "driven" / "fields.vtu") == "<VTK/>");
  }

  SECTION("Symlink is left behind for directory")
  {
    CHECK(fs::is_symlink(temp_dir / "paraview"));
    CHECK(fs::is_directory(temp_dir / "paraview" / "driven"));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration handles two consecutive iterations", "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Simulate first iteration output.
  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "iter1");

  SaveIteration(comm, temp_dir, 1, 1);

  // After first save: domain-E.csv is a symlink to iteration1/domain-E.csv.
  CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
  CHECK(ReadFile(temp_dir / "domain-E.csv") == "iter1");

  // Simulate second iteration: solver overwrites symlink with a real file.
  fs::remove(temp_dir / "domain-E.csv");
  CreateFile(temp_dir / "domain-E.csv", "iter2");
  CHECK(!fs::is_symlink(temp_dir / "domain-E.csv"));

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("First iteration data is preserved")
  {
    CHECK(ReadFile(temp_dir / "iteration1" / "domain-E.csv") == "iter1");
  }

  SECTION("Second iteration has new data")
  {
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "iter2");
  }

  SECTION("Symlink now points to second iteration")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(fs::read_symlink(temp_dir / "domain-E.csv") ==
          fs::path("iteration2") / "domain-E.csv");
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "iter2");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration preserves old symlinks for files not reproduced by solver",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // First iteration produces two files.
  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "iter1");
  CreateFile(temp_dir / "extra.csv", "extra_data");

  SaveIteration(comm, temp_dir, 1, 1);

  // Simulate second iteration that only produces domain-E.csv (not extra.csv).
  fs::remove(temp_dir / "domain-E.csv");
  CreateFile(temp_dir / "domain-E.csv", "iter2");

  // extra.csv is still a symlink from iteration 1.
  CHECK(fs::is_symlink(temp_dir / "extra.csv"));

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("Old symlink is preserved, still accessible")
  {
    CHECK(fs::is_symlink(temp_dir / "extra.csv"));
    CHECK(ReadFile(temp_dir / "extra.csv") == "extra_data");
  }

  SECTION("New file is moved to iteration2")
  {
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "iter2");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveAdaptMesh preserves an archived mesh between iterations",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  auto serial_mesh = SingleTetMesh();
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  mesh::RebalanceMesh(iodata, parallel_mesh);

  // Save the first adapted mesh contents for comparison.
  auto mesh_file = temp_dir / "model.mesh";
  auto first_mesh = ReadFile(mesh_file);
  REQUIRE_FALSE(first_mesh.empty());

  // Archive the first mesh, leaving a top-level symlink to it.
  SaveIteration(comm, temp_dir, 1, 1);
  REQUIRE(fs::is_symlink(mesh_file));
  REQUIRE(ReadFile(temp_dir / "iteration1" / "model.mesh") == first_mesh);

  // Simulate the next AMR iteration with a changed mesh.
  parallel_mesh->GetVertex(0)[0] = -1.0;
  mesh::RebalanceMesh(iodata, parallel_mesh);

  // The new mesh replaces the symlink without changing the archived mesh.
  CHECK_FALSE(fs::is_symlink(mesh_file));
  CHECK(ReadFile(mesh_file) != first_mesh);
  CHECK(ReadFile(temp_dir / "iteration1" / "model.mesh") == first_mesh);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh reports the true topological entity counts",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  auto serial_mesh = SingleTetMesh();
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);

  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);

  // The counts are computed on the root rank, where the gathered serial mesh is complete.
  const long long global_ne = parallel_mesh->GetGlobalNE();
  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    // A single tetrahedron has 4 vertices, 6 edges, 4 faces, and 1 element. These are the
    // order-independent conforming entity counts (H1(1)/ND(1)/RT(0)/L2(0) true sizes).
    CHECK(counts.true_vertices == 4);
    CHECK(counts.true_edges == 6);
    CHECK(counts.true_faces == 4);
    CHECK(counts.elements == 1);
    // The element count must agree with the mesh's own global element count.
    CHECK(counts.elements == global_ne);
    // SingleTetMesh tags its one domain 1 and its four boundary faces 1..4.
    CHECK(counts.domain_attributes == std::vector<int>{1});
    CHECK(counts.boundary_attributes == std::vector<int>{1, 2, 3, 4});
    // The counts compose into higher-order true sizes: an order-2 (P2) H1 space on a
    // tetrahedron has one dof per vertex and one per edge (no face/interior dofs at p=2),
    // so its true size must equal true_vertices + true_edges.
    mfem::H1_FECollection h1_p2(2, 3);
    mfem::FiniteElementSpace fes_p2(&serial_mesh, &h1_p2);
    CHECK(fes_p2.GetTrueVSize() == counts.true_vertices + counts.true_edges);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh discounts hanging entities on a nonconforming mesh",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  // Build a multi-element tetrahedral mesh and refine a subset nonconformally so that
  // hanging nodes appear at the refined/unrefined interfaces.
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  serial_mesh.EnsureNCMesh(true);
  mfem::Array<int> refs;
  for (int i = 0; i < serial_mesh.GetNE(); i += 3)
  {
    refs.Append(i);
  }
  serial_mesh.GeneralRefinement(refs, 1);  // 1 => nonconforming
  REQUIRE(serial_mesh.Nonconforming());

  // Raw leaf entity counts of the (global) serial mesh, before distribution.
  const long long leaf_vertices = serial_mesh.GetNV();
  const long long leaf_edges = serial_mesh.GetNEdges();
  const long long leaf_faces = serial_mesh.GetNFaces();
  const long long leaf_elements = serial_mesh.GetNE();

  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    // The true (conforming) counts discount the constrained hanging vertices, edges, and
    // faces, so they are strictly smaller than the raw leaf counts. This is the whole
    // reason leaf counts cannot be used directly for sizing a nonconforming mesh.
    CHECK(counts.true_vertices > 0);
    CHECK(counts.true_edges > 0);
    CHECK(counts.true_faces > 0);
    CHECK(counts.true_vertices < leaf_vertices);
    CHECK(counts.true_edges < leaf_edges);
    CHECK(counts.true_faces < leaf_faces);
    // Elements are never constrained, so the element count is not discounted.
    CHECK(counts.elements == leaf_elements);
    // The true entity counts must reproduce the exact conforming DOF count of the H1/ND
    // spaces Palace solves with, across orders, on this nonconforming mesh -- the property
    // that lets them size a re-run at any order (stronger than the leaf-count
    // inequalities).
    CheckDofDecomposition(serial_mesh, counts);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh leaves entity counts unset without SaveAdaptMesh",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = false;

  auto serial_mesh = SingleTetMesh();
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);

  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);

  // With SaveAdaptMesh disabled no mesh is written, so no counts are produced. Also assert
  // the save side-effect did not happen, tying valid == false to its actual cause rather
  // than merely to the struct's default.
  CHECK_FALSE(counts.valid);
  CHECK_FALSE(fs::exists(temp_dir / "model.mesh"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveAdaptMesh writes the Mesh block to palace.json",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  auto serial_mesh = SingleTetMesh();
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);

  // The BaseSolver constructor writes the initial palace.json (root only); SaveMetadata
  // then reads it, adds the "Mesh" block, and writes it back -- exactly the production
  // path.
  TestSolver solver(iodata, Mpi::Root(comm));
  if (counts.valid)
  {
    solver.SaveMetadata(counts);
  }

  if (Mpi::Root(comm))
  {
    const auto meta = nlohmann::json::parse(ReadFile(temp_dir / "palace.json"));
    REQUIRE(meta.contains("Mesh"));
    const auto &m = meta.at("Mesh");
    CHECK(m.at("Dimension").get<int>() == 3);
    CHECK(m.at("TrueVertices").get<long long>() == 4);
    CHECK(m.at("TrueEdges").get<long long>() == 6);
    CHECK(m.at("TrueFaces").get<long long>() == 4);
    CHECK(m.at("Elements").get<long long>() == 1);
    CHECK(m.at("DomainAttributes").get<std::vector<int>>() == std::vector<int>{1});
    CHECK(m.at("BoundaryAttributes").get<std::vector<int>>() ==
          std::vector<int>{1, 2, 3, 4});
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveIteration handles dirty output directory from previous run",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Simulate a previous run that left files in the iteration folder.
  CreateFile(temp_dir / "palace.json", "{}");
  fs::create_directories(temp_dir / "iteration1" / "paraview" / "driven");
  CreateFile(temp_dir / "iteration1" / "paraview" / "driven" / "old.vtu", "old");
  CreateFile(temp_dir / "iteration1" / "domain-E.csv", "old_data");
  fs::create_directories(temp_dir / "iteration1" / ".palace-archive-claim");

  // New run produces fresh output.
  CreateFile(temp_dir / "domain-E.csv", "new_data");
  fs::create_directories(temp_dir / "paraview" / "driven");
  CreateFile(temp_dir / "paraview" / "driven" / "new.vtu", "new");

  // Should not throw despite destination files existing.
  SaveIteration(comm, temp_dir, 1, 1);

  SECTION("New data replaces old in iteration folder")
  {
    CHECK(ReadFile(temp_dir / "iteration1" / "domain-E.csv") == "new_data");
    CHECK(fs::exists(temp_dir / "iteration1" / "paraview" / "driven" / "new.vtu"));
    CHECK(!fs::exists(temp_dir / "iteration1" / "paraview" / "driven" / "old.vtu"));
  }

  SECTION("Symlinks point to new data")
  {
    CHECK(fs::is_symlink(temp_dir / "domain-E.csv"));
    CHECK(ReadFile(temp_dir / "domain-E.csv") == "new_data");
  }

  SECTION("Stale archive claims are removed")
  {
    CHECK_FALSE(fs::exists(temp_dir / "iteration1" / ".palace-archive-claim"));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SaveIteration skips iteration subfolders",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "domain-E.csv", "data");

  // Create a pre-existing iteration folder.
  fs::create_directories(temp_dir / "iteration1");
  CreateFile(temp_dir / "iteration1" / "old.csv", "old");

  SaveIteration(comm, temp_dir, 2, 1);

  SECTION("Pre-existing iteration folder is untouched")
  {
    CHECK(fs::is_directory(temp_dir / "iteration1"));
    CHECK(ReadFile(temp_dir / "iteration1" / "old.csv") == "old");
  }

  SECTION("New iteration folder is created")
  {
    CHECK(fs::is_directory(temp_dir / "iteration2"));
    CHECK(ReadFile(temp_dir / "iteration2" / "domain-E.csv") == "data");
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SaveIteration is MPI collective",
                 "[basesolver][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if (Mpi::Root(comm))
  {
    CreateFile(temp_dir / "palace.json", "{}");
    fs::create_directories(temp_dir / "gridfunction" / "electrostatic");
    CreateFile(temp_dir / "gridfunction" / "electrostatic" / "V.gf.000000", "field");
  }
  Mpi::Barrier(comm);

  bool setup_complete =
      fs::is_regular_file(temp_dir / "gridfunction" / "electrostatic" / "V.gf.000000");
  Mpi::GlobalAnd(1, &setup_complete, comm);
  REQUIRE(setup_complete);

  SaveIteration(comm, temp_dir, 1, 1);

  CHECK(fs::is_symlink(temp_dir / "gridfunction"));
  CHECK(fs::is_regular_file(temp_dir / "iteration1" / "gridfunction" / "electrostatic" /
                            "V.gf.000000"));
}
