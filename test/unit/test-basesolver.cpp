// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/ncmeshio.hpp"

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

// Reconstruct the exact conforming DOF count (`GetTrueVSize`) from the per-geometry
// true entity counts. Agreement with `GetTrueVSize` verifies that constrained entities
// are excluded and face geometries are classified correctly on nonconforming meshes.
long long Reconstruct(const mesh::MeshEntityCounts &c, mfem::FiniteElementCollection &fec)
{
  long long n = c.true_vertices * fec.DofForGeometry(mfem::Geometry::POINT) +
                c.true_edges * fec.DofForGeometry(mfem::Geometry::SEGMENT);
  for (const auto &[g, cnt] : c.true_faces)
  {
    n += cnt * fec.DofForGeometry(g);
  }
  for (const auto &[g, cnt] : c.cells)
  {
    n += cnt * fec.DofForGeometry(g);
  }
  return n;
}

// Assert the per-geometry true counts reproduce the EXACT conforming true DOF size of the
// H1, ND, and RT spaces Palace sizes with, across orders, on `serial_mesh`. RT (p=0..2) is
// the ONLY space with FACE dofs, so its match directly validates the tri/quad TRUE-face
// split -- the quantity a prior RT-solve corrupted. `counts` come from a gathered serial
// copy whose topology is identical to `serial_mesh`, so the true sizes must match exactly.
// Runs on the caller's rank; call inside a root guard (counts are valid only on root).
void CheckReconstructionOracle(const mesh::MeshEntityCounts &counts,
                               mfem::Mesh &serial_mesh)
{
  const int dim = serial_mesh.Dimension();
  for (int p = 1; p <= 3; p++)
  {
    mfem::H1_FECollection fec(p, dim);
    CHECK(Reconstruct(counts, fec) ==
          mfem::FiniteElementSpace(&serial_mesh, &fec).GetTrueVSize());
  }
  for (int p = 1; p <= 3; p++)
  {
    mfem::ND_FECollection fec(p, dim);
    CHECK(Reconstruct(counts, fec) ==
          mfem::FiniteElementSpace(&serial_mesh, &fec).GetTrueVSize());
  }
  if (dim == 3)
  {
    for (int p = 0; p <= 2; p++)
    {
      mfem::RT_FECollection fec(p, dim);
      CHECK(Reconstruct(counts, fec) ==
            mfem::FiniteElementSpace(&serial_mesh, &fec).GetTrueVSize());
    }
  }
}

// Distribute `serial_mesh`, gather geometry-resolved counts while writing the saved mesh,
// then collectively complete nonconforming true vertex and edge counts.
mesh::MeshEntityCounts ComputeCounts(MPI_Comm comm, const fs::path &out,
                                     mfem::Mesh &serial_mesh)
{
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = out.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);
  mesh::CompleteMeshEntityCounts(*parallel_mesh, counts);
  return counts;
}

// Build a nonconforming prism (WEDGE) mesh: extrude a 2D triangular mesh into wedges, mark
// it nonconforming, and refine a subset so hanging nodes appear. A prism carries BOTH
// triangular (top/bottom) and quadrilateral (side) faces, so it exercises the mixed
// tri/quad TRUE-face split with a genuinely nonzero answer for both geometries.
mfem::Mesh BuildNCPrismMesh()
{
  auto base2d = mfem::Mesh::MakeCartesian2D(3, 3, mfem::Element::TRIANGLE);
  std::unique_ptr<mfem::Mesh> extruded(mfem::Extrude2D(&base2d, 3, 1.0));
  mfem::Mesh serial_mesh(*extruded);
  serial_mesh.EnsureNCMesh(true);
  mfem::Array<int> refs;
  for (int i = 0; i < serial_mesh.GetNE(); i += 3)
  {
    refs.Append(i);
  }
  serial_mesh.GeneralRefinement(refs, 1);  // 1 => nonconforming
  return serial_mesh;
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

  auto iter_dir = SaveIteration(comm, temp_dir, 1, 2);

  SECTION("Iteration subfolder is created and returned")
  {
    CHECK(iter_dir == temp_dir / "iteration01");
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

  // RebalanceMesh gzip-compresses the saved mesh (".meshgz") when MFEM is built with zlib,
  // and writes it uncompressed (".mesh") otherwise; the archived bytes match the write
  // either way, so this test is agnostic to which path ran.
#if defined(MFEM_USE_ZLIB)
  const std::string saved_mesh = "model.meshgz";
#else
  const std::string saved_mesh = "model.mesh";
#endif

  // Save the first adapted mesh contents for comparison.
  auto mesh_file = temp_dir / saved_mesh;
  auto first_mesh = ReadFile(mesh_file);
  REQUIRE_FALSE(first_mesh.empty());

  // Archive the first mesh, leaving a top-level symlink to it.
  SaveIteration(comm, temp_dir, 1, 1);
  REQUIRE(fs::is_symlink(mesh_file));
  REQUIRE(ReadFile(temp_dir / "iteration1" / saved_mesh) == first_mesh);

  // Simulate the next AMR iteration with a changed mesh.
  parallel_mesh->GetVertex(0)[0] = -1.0;
  mesh::RebalanceMesh(iodata, parallel_mesh);

  // The new mesh replaces the symlink without changing the archived mesh.
  CHECK_FALSE(fs::is_symlink(mesh_file));
  CHECK(ReadFile(mesh_file) != first_mesh);
  CHECK(ReadFile(temp_dir / "iteration1" / saved_mesh) == first_mesh);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Mesh loader reads a gzip-compressed mesh (.meshgz round-trip)",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto serial_mesh = SingleTetMesh();

  // Write the same mesh two ways: gzip-compressed (the form SaveAdaptMesh emits under
  // MFEM_USE_ZLIB) and uncompressed. Loading the compressed mesh must yield the same
  // topology as the uncompressed one; before LoadMesh read through named_ifgzstream, a
  // ".meshgz" aborted with "Unknown input mesh format" (the reader saw raw gzip bytes).
  const auto gz_path = temp_dir / "mesh.meshgz";
  const auto plain_path = temp_dir / "mesh.mesh";
  {
    mfem::ofgzstream os(gz_path.string(), /*compression=*/true);
    serial_mesh.Print(os);
  }
  {
    std::ofstream os(plain_path);
    serial_mesh.Print(os);
  }

  Units units(1.0, 1.0);
  auto load = [&](const fs::path &p)
  {
    IoData iodata(units);
    iodata.model.mesh = p.string();
    iodata.model.clean_unused_elements = false;
    return mesh::Load(iodata, comm);
  };
  auto gz = load(gz_path);
  auto plain = load(plain_path);
  REQUIRE(gz);
  REQUIRE(plain);
  CHECK(gz->Dimension() == plain->Dimension());
  CHECK(gz->GetNE() == plain->GetNE());
  CHECK(gz->GetNV() == plain->GetNV());

#if defined(MFEM_USE_ZLIB)
  // Confirm the ".meshgz" really is gzip (magic bytes 0x1f 0x8b), so the load above
  // exercised decompression rather than a silently-uncompressed file.
  std::ifstream raw(gz_path, std::ios::binary);
  unsigned char magic[2] = {0, 0};
  raw.read(reinterpret_cast<char *>(magic), 2);
  CHECK(magic[0] == 0x1f);
  CHECK(magic[1] == 0x8b);
#endif
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
    // A single tetrahedron has 4 vertices, 6 edges, 4 triangular faces, and is one
    // TETRAHEDRON cell. These are the order-independent conforming entity counts
    // (H1(1)/ND(1)/RT(0)/L2(0) true sizes), reported per geometry.
    CHECK(counts.true_vertices == 4);
    CHECK(counts.true_edges == 6);
    CHECK(counts.true_faces.at(mfem::Geometry::TRIANGLE) == 4);
    // A tet has no quadrilateral faces, so the SQUARE key must be absent (not
    // present-zero).
    CHECK_FALSE(counts.true_faces.contains(mfem::Geometry::SQUARE));
    CHECK(counts.cells.at(mfem::Geometry::TETRAHEDRON) == 1);
    CHECK(counts.cells.size() == 1);
    // The (single) cell count must agree with the mesh's own global element count.
    CHECK(counts.cells.at(mfem::Geometry::TETRAHEDRON) == global_ne);
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
    CHECK_FALSE(counts.valid);
    CHECK(counts.true_vertices == 0);
    CHECK(counts.true_edges == 0);
  }
  mesh::CompleteMeshEntityCounts(*parallel_mesh, counts);

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    // The true (conforming) counts discount the constrained hanging vertices, edges, and
    // faces, so they are strictly smaller than the raw leaf counts. This is the whole
    // reason leaf counts cannot be used directly for sizing a nonconforming mesh. This
    // all-tet mesh has only triangular faces, so the SQUARE key must be absent.
    CHECK(counts.true_vertices > 0);
    CHECK(counts.true_edges > 0);
    CHECK(counts.true_faces.at(mfem::Geometry::TRIANGLE) > 0);
    CHECK_FALSE(counts.true_faces.contains(mfem::Geometry::SQUARE));
    CHECK(counts.true_vertices < leaf_vertices);
    CHECK(counts.true_edges < leaf_edges);
    // Cheap sanity check that the face discount happened.
    CHECK(counts.true_faces.at(mfem::Geometry::TRIANGLE) < leaf_faces);
    // Cells are never constrained, so the raw per-geometry cell count is not discounted.
    CHECK(counts.cells.at(mfem::Geometry::TETRAHEDRON) == leaf_elements);
    // No entity count may ever be negative. Guards the class of bug where the prior
    // nonconforming RT true-size solve underflowed into a negative face count (a real AMR
    // run produced TrueFaces.quadrilateral = -3447 on a pure-tet mesh). The direct NC
    // face-list count that replaced the solve cannot go negative by construction.
    for (const auto &[g, c] : counts.true_faces)
    {
      CHECK(c >= 0);
    }
    for (const auto &[g, c] : counts.cells)
    {
      CHECK(c >= 0);
    }
    CHECK(counts.true_vertices >= 0);
    CHECK(counts.true_edges >= 0);
    // Reconstruction-vs-GetTrueVSize oracle across H1, ND, and RT (RT validates the
    // tri/quad TRUE-face split that the old solve corrupted). See
    // CheckReconstructionOracle.
    CheckReconstructionOracle(counts, serial_mesh);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh reports exact true counts on a nonconforming prism mesh",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // A refined nonconforming prism (wedge) mesh carries BOTH triangular and quadrilateral
  // faces, so it exercises the mixed tri/quad TRUE-face split with a nonzero answer for
  // both geometries -- the case the removed RT(0)/RT(1) solve got wrong.
  auto serial_mesh = BuildNCPrismMesh();
  REQUIRE(serial_mesh.Nonconforming());

  const auto counts = ComputeCounts(comm, temp_dir, serial_mesh);

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    // Prisms bear both triangular and quadrilateral faces; both must be present and > 0.
    CHECK(counts.true_faces.at(mfem::Geometry::TRIANGLE) > 0);
    CHECK(counts.true_faces.at(mfem::Geometry::SQUARE) > 0);
    // Cells are all prisms.
    CHECK(counts.cells.count(mfem::Geometry::PRISM) == 1);
    for (const auto &[g, c] : counts.true_faces)
    {
      CHECK(c >= 0);
    }
    // The full H1/ND/RT sweep. RT directly validates the tri/quad split on a mixed mesh.
    CheckReconstructionOracle(counts, serial_mesh);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh reports exact true counts on a nonconforming hex mesh",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // A pure-hexahedral nonconforming mesh must report quadrilateral true faces only.
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON);
  serial_mesh.EnsureNCMesh();
  mfem::Array<int> refs;
  for (int i = 0; i < serial_mesh.GetNE(); i += 3)
  {
    refs.Append(i);
  }
  serial_mesh.GeneralRefinement(refs, 1);  // 1 => nonconforming
  REQUIRE(serial_mesh.Nonconforming());

  const auto counts = ComputeCounts(comm, temp_dir, serial_mesh);

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    // A hex mesh has only quadrilateral faces: SQUARE present-positive, TRIANGLE absent.
    CHECK(counts.true_faces.at(mfem::Geometry::SQUARE) > 0);
    CHECK_FALSE(counts.true_faces.contains(mfem::Geometry::TRIANGLE));
    CHECK(counts.cells.count(mfem::Geometry::CUBE) == 1);
    CheckReconstructionOracle(counts, serial_mesh);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh counts match between a conforming mesh and its unrefined "
                 "nonconforming twin",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // For each element geometry: compute counts on the conforming mesh (raw-count path), then
  // flag the same mesh nonconforming WITHOUT refining (topologically identical, routed
  // through the NC FE-space + NC-face-list path) and compute again. The two must agree
  // exactly -- the robust guard that the NC counting reduces to the conforming answer when
  // there are no hanging entities. Each geometry uses its own subdirectory so the written
  // saved-mesh files do not collide.
  auto check_twin = [&](const std::string &name, const mfem::Mesh &base)
  {
    const auto conf_dir = temp_dir / (name + "-conf");
    const auto nc_dir = temp_dir / (name + "-nc");
    fs::create_directories(conf_dir);
    fs::create_directories(nc_dir);

    mfem::Mesh conf_mesh(base);
    const auto conf = ComputeCounts(comm, conf_dir, conf_mesh);

    mfem::Mesh nc_mesh(base);
    nc_mesh.EnsureNCMesh(true);
    REQUIRE(nc_mesh.Nonconforming());
    const auto nc = ComputeCounts(comm, nc_dir, nc_mesh);

    if (Mpi::Root(comm))
    {
      INFO("geometry: " << name);
      REQUIRE(conf.valid);
      REQUIRE(nc.valid);
      CHECK(conf.dim == nc.dim);
      CHECK(conf.true_vertices == nc.true_vertices);
      CHECK(conf.true_edges == nc.true_edges);
      CHECK(conf.true_faces == nc.true_faces);
      CHECK(conf.cells == nc.cells);
    }
  };

  check_twin("tet", mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
  check_twin("hex", mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON));
  {
    auto base2d = mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::TRIANGLE);
    std::unique_ptr<mfem::Mesh> prism(mfem::Extrude2D(&base2d, 2, 1.0));
    check_twin("prism", *prism);
  }
  check_twin("tri2d", mfem::Mesh::MakeCartesian2D(3, 3, mfem::Element::TRIANGLE));
  check_twin("quad2d", mfem::Mesh::MakeCartesian2D(3, 3, mfem::Element::QUADRILATERAL));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "RebalanceMesh entity counts reflect only the final mesh across iterations",
    "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  // Regression for the accumulation bug: one MeshEntityCounts reused across AMR iterations
  // must be reset each call, so the accumulating members (cells, true_faces) reflect ONLY
  // the last mesh, never the sum over iterations.
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);

  mesh::MeshEntityCounts counts;

  // Iteration 1.
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);
  const long long ne1 = parallel_mesh->GetGlobalNE();

  // Iteration 2: refine (still conforming) so the mesh genuinely changes, reusing `counts`.
  parallel_mesh->UniformRefinement();
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);
  const long long ne2 = parallel_mesh->GetGlobalNE();

  // Build the expected final serial mesh (same single uniform refinement) for the oracle.
  serial_mesh.UniformRefinement();

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    REQUIRE(ne2 > ne1);
    // Cells must equal ONLY the final element count, not ne1 + ne2. Without the per-call
    // reset this would be inflated to the running sum.
    CHECK(counts.cells.at(mfem::Geometry::TETRAHEDRON) == ne2);
    long long total_cells = 0;
    for (const auto &[g, c] : counts.cells)
    {
      total_cells += c;
    }
    CHECK(total_cells == ne2);
    // The reconstruction oracle on the final mesh also fails if true_faces (or any true
    // count) had accumulated across the two iterations.
    CheckReconstructionOracle(counts, serial_mesh);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveMetadata writes no SavedAdaptedMesh block when counts are invalid",
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

  // With SaveAdaptMesh disabled (a zero-AMR-iteration run also lands here), counts are
  // invalid and the driver skips the block; palace.json must contain no SavedAdaptedMesh.
  TestSolver solver(iodata, Mpi::Root(comm));
  if (counts.valid)
  {
    solver.SaveMetadata(counts);
  }
  if (Mpi::Root(comm))
  {
    CHECK_FALSE(counts.valid);
    const auto meta = nlohmann::json::parse(ReadFile(temp_dir / "palace.json"));
    CHECK_FALSE(meta.contains("SavedAdaptedMesh"));
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
  // than merely to the struct's default. Check both possible names (compressed and not) so
  // the guard holds regardless of whether MFEM was built with zlib.
  CHECK_FALSE(counts.valid);
  CHECK_FALSE(fs::exists(temp_dir / "model.mesh"));
  CHECK_FALSE(fs::exists(temp_dir / "model.meshgz"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SaveAdaptMesh writes a SavedAdaptedMesh block to palace.json",
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
  // then reads it, adds the "SavedAdaptedMesh" block, and writes it back -- exactly the
  // path.
  TestSolver solver(iodata, Mpi::Root(comm));
  if (counts.valid)
  {
    solver.SaveMetadata(counts);
  }

  if (Mpi::Root(comm))
  {
    const auto meta = nlohmann::json::parse(ReadFile(temp_dir / "palace.json"));
    REQUIRE(meta.contains("SavedAdaptedMesh"));
    const auto &m = meta.at("SavedAdaptedMesh");
    CHECK(m.at("Dimension").get<int>() == 3);
    CHECK(m.at("TrueVertices").get<long long>() == 4);
    CHECK(m.at("TrueEdges").get<long long>() == 6);
    // Per-geometry blocks: a single tet has four triangular faces and is one tetrahedron.
    CHECK(m.at("TrueFaces").at("triangle").get<long long>() == 4);
    // A tet has no quadrilateral faces, so that key must be absent (not present-zero).
    CHECK_FALSE(m.at("TrueFaces").contains("quadrilateral"));
    CHECK(m.at("Cells").at("tetrahedron").get<long long>() == 1);
    CHECK(m.at("DomainAttributes").get<std::vector<int>>() == std::vector<int>{1});
    CHECK(m.at("BoundaryAttributes").get<std::vector<int>>() ==
          std::vector<int>{1, 2, 3, 4});
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh reports the true topological entity counts in parallel",
                 "[basesolver][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.problem.output = temp_dir.string();
  iodata.model.mesh = "model.msh";
  iodata.model.refinement.save_adapt_mesh = true;

  // Build a nonconforming tet mesh identically on every rank, distribute it across the
  // ranks, then RebalanceMesh gathers geometry-resolved counts on root and the collective
  // completion step obtains true vertex and edge totals from distributed spaces. The
  // reported counts must equal the single-rank answer regardless of partitioning, while
  // non-root count records remain invalid.
  auto serial_mesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  serial_mesh.EnsureNCMesh(true);
  mfem::Array<int> refs;
  for (int i = 0; i < serial_mesh.GetNE(); i += 3)
  {
    refs.Append(i);
  }
  serial_mesh.GeneralRefinement(refs, 1);  // 1 => nonconforming
  REQUIRE(serial_mesh.Nonconforming());

  auto parallel_mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  mesh::MeshEntityCounts counts;
  mesh::RebalanceMesh(iodata, parallel_mesh, &counts);
  mesh::CompleteMeshEntityCounts(*parallel_mesh, counts);

  if (Mpi::Root(comm))
  {
    REQUIRE(counts.valid);
    CHECK(counts.dim == 3);
    CHECK(counts.true_faces.at(mfem::Geometry::TRIANGLE) > 0);
    CHECK_FALSE(counts.true_faces.contains(mfem::Geometry::SQUARE));
    // Partition-invariant: the gathered root counts reproduce the serial true DOF sizes.
    CheckReconstructionOracle(counts, serial_mesh);
  }
  else
  {
    // Only root owns the gathered geometry-resolved count record.
    CHECK_FALSE(counts.valid);
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

namespace
{

// Compare, on the root, the meshes read from two files, up to the order of their elements
// and boundary elements: entity counts, roots, hanging entities, the elements matched by
// attribute and centroid (vertices and nodes), the boundary, one refinement and one
// derefinement; and the entity counts of the parallel save.
void CompareSavedMeshes(const fs::path &file_a, const fs::path &file_b,
                        const mesh::MeshEntityCounts &counts)
{
  mfem::named_ifgzstream fa(file_a.string()), fb(file_b.string());
  mfem::Mesh a(fa), b(fb);
  CHECK((a.Nonconforming() && b.Nonconforming()));
  CHECK(a.GetNV() == b.GetNV());
  CHECK(a.GetNEdges() == b.GetNEdges());
  CHECK(a.GetNFaces() == b.GetNFaces());
  CHECK(a.ncmesh->GetNumRootElements() == b.ncmesh->GetNumRootElements());
  for (auto list : {&mfem::NCMesh::GetEdgeList, &mfem::NCMesh::GetFaceList})
  {
    const auto &la = (a.ncmesh->*list)(), &lb = (b.ncmesh->*list)();
    CHECK(la.conforming.Size() == lb.conforming.Size());
    CHECK(la.masters.Size() == lb.masters.Size());
    CHECK(la.slaves.Size() == lb.slaves.Size());
  }
  std::map<mfem::Geometry::Type, long long> cells, true_faces;
  for (int e = 0; e < a.GetNE(); e++)
  {
    cells[a.GetElementBaseGeometry(e)]++;
  }
  for (int f = 0; a.Dimension() == 3 && f < a.GetNFaces(); f++)
  {
    true_faces[a.GetFaceGeometry(f)] +=
        a.ncmesh->GetFaceList().GetMeshIdType(f) != mfem::NCMesh::NCList::MeshIdType::SLAVE;
  }
  std::erase_if(true_faces, [](const auto &c) { return c.second == 0; });
  CHECK(counts.dim == a.Dimension());
  CHECK(counts.cells == cells);
  CHECK(counts.true_faces == true_faces);
  CHECK(counts.domain_attributes ==
        std::vector<int>(a.attributes.begin(), a.attributes.end()));
  CHECK(counts.boundary_attributes ==
        std::vector<int>(a.bdr_attributes.begin(), a.bdr_attributes.end()));
  auto Order = [](mfem::Mesh &m, bool bdr)
  {
    std::vector<std::pair<std::vector<long long>, int>> keys;
    for (int e = 0; e < (bdr ? m.GetNBE() : m.GetNE()); e++)
    {
      mfem::Array<int> v;
      bdr ? m.GetBdrElementVertices(e, v) : m.GetElementVertices(e, v);
      std::vector<long long> key = {bdr ? m.GetBdrAttribute(e) : m.GetAttribute(e)};
      for (int d = 0; d < m.SpaceDimension(); d++)
      {
        double c = 0.0;
        for (int i : v)
        {
          c += m.GetVertex(i)[d] / v.Size();
        }
        key.push_back(std::llround(c * 1e9));
      }
      keys.emplace_back(key, e);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
  };
  const auto ka = Order(a, false), kb = Order(b, false);
  const auto ba = Order(a, true), bb = Order(b, true);
  auto SameKeys = [](const auto &x, const auto &y)
  {
    return std::equal(x.begin(), x.end(), y.begin(), y.end(),
                      [](const auto &p, const auto &q) { return p.first == q.first; });
  };
  const bool same = SameKeys(ka, kb);
  CHECK(same);
  CHECK(SameKeys(ba, bb));
  if (!same)
  {
    return;  // the elements do not correspond (the other processes go on)
  }
  double max_diff = 0.0;
  mfem::Array<int> refs_a, refs_b;
  for (std::size_t i = 0; i < ka.size(); i++)
  {
    const int ea = ka[i].second, eb = kb[i].second;
    mfem::Array<int> va, vb;
    a.GetElementVertices(ea, va);
    b.GetElementVertices(eb, vb);
    for (int j = 0; j < va.Size(); j++)
    {
      for (int d = 0; d < a.SpaceDimension(); d++)
      {
        max_diff =
            std::max(max_diff, std::abs(a.GetVertex(va[j])[d] - b.GetVertex(vb[j])[d]));
      }
    }
    if (a.GetNodes())
    {
      mfem::Array<int> da, db;
      mfem::Vector xa, xb;
      a.GetNodes()->FESpace()->GetElementVDofs(ea, da);
      b.GetNodes()->FESpace()->GetElementVDofs(eb, db);
      a.GetNodes()->GetSubVector(da, xa);
      b.GetNodes()->GetSubVector(db, xb);
      xa -= xb;
      max_diff = std::max(max_diff, xa.Normlinf());
    }
    if (i % 5 == 0)
    {
      refs_a.Append(ea);
      refs_b.Append(eb);
    }
  }
  CHECK(max_diff < 1e-11);
  mfem::Mesh ra(a), rb(b);
  ra.GeneralRefinement(refs_a, 1);
  rb.GeneralRefinement(refs_b, 1);
  CHECK(ra.GetNE() == rb.GetNE());
  CHECK(ra.GetNV() == rb.GetNV());
  mfem::Vector zero(a.GetNE());
  zero = 0.0;
  a.DerefineByError(zero, 1.0);
  b.DerefineByError(zero, 1.0);
  CHECK(a.GetNE() == b.GetNE());
  CHECK(a.GetNV() == b.GetNV());
}

// Refine and scatter nonconforming meshes, and compare their parallel saves (each saved
// mesh is also read back and saved again) with a gathered save.
void TestParallelNCSave(MPI_Comm comm, const fs::path &temp_dir, mfem::Element::Type type,
                        int order)
{
  const int P = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const bool is3d =
      (type != mfem::Element::TRIANGLE && type != mfem::Element::QUADRILATERAL);
  auto serial_mesh = is3d ? mfem::Mesh::MakeCartesian3D(3, 3, 2, type)
                          : mfem::Mesh::MakeCartesian2D(5, 4, type);
  serial_mesh.EnsureNCMesh(true);
  if (order > 0)
  {
    serial_mesh.SetCurvature(order, false, -1,
                             order == 2 ? mfem::Ordering::byNODES : mfem::Ordering::byVDIM);
    serial_mesh.Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(0) += 0.05 * std::sin(3.0 * x(1) + x(x.Size() - 1));
          y(1) += 0.04 * std::cos(2.0 * x(0));
        });
  }
  auto mesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  const double L = 2.5;  // the length scale of the saved mesh
  mesh::MeshEntityCounts counts;
  auto Check = [&](const std::string &name, bool parallel)
  {
    const auto gathered = temp_dir / (name + "-gathered.mesh");
    auto file = temp_dir / (name + ".meshgz");
    {
      mfem::ParMesh smesh(*mesh);
      mfem::Array<int> part(smesh.GetNE());
      part = 0;
      smesh.Rebalance(part);
      if (Mpi::Root(comm))
      {
        mesh::DimensionalizeMesh(smesh, L);
        std::ofstream fo(gathered);
        fo << std::scientific;
        fo.precision(std::numeric_limits<double>::max_digits10);
        smesh.Mesh::Print(fo);
      }
    }
    std::string reason;
    if (parallel)
    {
      REQUIRE(mesh::PrintParallelNCMesh(*mesh, L, file.string(), order != 1, &counts,
                                        &reason));  // also uncompressed
    }
    else
    {
      // Through RebalanceMesh, which rebalances a nonconforming mesh, then saves it in
      // parallel (not gathered: as PrintParallelNCMesh saves the rebalanced mesh).
      Units units(1.0, L);
      IoData iodata(units);
      iodata.problem.output = temp_dir.string();
      iodata.model.mesh = name + ".msh";
      iodata.model.refinement.save_adapt_mesh = true;
      mesh::RebalanceMesh(iodata, mesh, &counts);
      const bool compressed = fs::exists(file);
      file.replace_extension(compressed ? ".meshgz" : ".mesh");
      const auto direct = temp_dir / (name + "-direct.mesh");
      REQUIRE(mesh::PrintParallelNCMesh(*mesh, L, direct.string(), compressed));
      if (Mpi::Root(comm))
      {
        CHECK(ReadFile(file) == ReadFile(direct));
      }
    }
    if (Mpi::Root(comm))
    {
      CompareSavedMeshes(gathered, file, counts);
    }
    Mpi::Barrier(comm);
    return file;
  };
  fs::path file;
  for (int pass = 0; pass < 3; pass++)
  {
    // Refine some elements, then, after the first pass, scatter the leaves over the
    // processes, so that roots are split between them.
    mfem::Array<int> refs;
    for (int e = 0; e < mesh->GetNE(); e++)
    {
      if ((e * 7 + rank * 3 + pass) % 4 == 0)
      {
        refs.Append(e);
      }
    }
    mesh->GeneralRefinement(refs, 1, pass);
    if (pass > 0)
    {
      mfem::Array<int> part(mesh->GetNE());
      for (int e = 0; e < mesh->GetNE(); e++)
      {
        part[e] = (pass == 1) ? (e / 3 + rank) % P : (e * 13 + rank + pass) % P;
      }
      mesh->Rebalance(part);
    }
    file = Check("pass" + std::to_string(pass), true);
  }

  // Read the last mesh back in parallel from an uncompressed file and from the last file
  // (compressed but for order 1), as a restart does, and check that they are the same. Then
  // refine it and the mesh alike (with a limit on the hanging nodes, which keeps their
  // differences), and compare them.
  {
    const auto plain = temp_dir / "plain.mesh", gathered = temp_dir / "pass2-gathered.mesh";
    REQUIRE(mesh::PrintParallelNCMesh(*mesh, L, plain.string(), false));
    std::string text;
    CHECK((order == 0 ||
           !mesh::ReadNCMeshCoarse(comm, gathered.string(), true,
                                   text)));  // not in the order of PrintParallelNCMesh
    auto Read = [&](const fs::path &path)
    {
      auto kept = mesh::ReadNCMeshCoarse(comm, path.string(), true, text);
      REQUIRE(kept);
      std::istringstream is(text);
      auto coarse = std::make_unique<mfem::Mesh>(is);
      mesh::NondimensionalizeMesh(*coarse, L);
      auto read = mesh::PartitionNCMeshFile(comm, *kept, std::move(coarse), L, 1.1);
      REQUIRE(read);
      REQUIRE(mesh::PrintParallelNCMesh(*read, L, path.string() + ".read", false));
      return read;
    };
    Read(file);
    auto saved = std::move(mesh);
    mesh = Read(plain);
    if (Mpi::Root(comm))
    {
      CHECK(ReadFile(file.string() + ".read") == ReadFile(plain.string() + ".read"));
    }

    // A file which does not match its coarse mesh is not read, rather than aborting.
    {
      auto kept = mesh::ReadNCMeshCoarse(comm, plain.string(), true, text);
      REQUIRE(kept);
      std::istringstream is(text);
      auto coarse = std::make_unique<mfem::Mesh>(is);
      coarse->UniformRefinement();
      CHECK_FALSE(mesh::PartitionNCMeshFile(comm, *kept, std::move(coarse), L, 1.1));
    }

    // Through mesh::Load and mesh::Partition, as a restart reads its mesh: in parallel from
    // the files of PrintParallelNCMesh, and from the gathered file in parallel or whole.
    auto Restart = [&](const fs::path &path)
    {
      IoData iodata(Units(1.0, L));
      iodata.model.mesh = path.string();
      iodata.model.clean_unused_elements = false;
      auto smesh = mesh::Load(iodata, comm);
      const bool parallel = static_cast<bool>(iodata.nc_mesh_file);
      if (smesh)
      {
        mesh::NondimensionalizeMesh(*smesh, L);
      }
      auto read = mesh::Partition(iodata, std::move(smesh), comm);
      REQUIRE(mesh::PrintParallelNCMesh(*read, L, path.string() + ".restart", false));
      if (Mpi::Root(comm))
      {
        CompareSavedMeshes(gathered, path.string() + ".restart", counts);
      }
      return parallel;
    };
    CHECK(Restart(plain));
    CHECK(Restart(file));
    Restart(gathered);
    for (auto *m : {saved.get(), mesh.get()})
    {
      mfem::Array<int> refs;
      for (int e = 0; e < m->GetNE(); e++)
      {
        mfem::Vector c;
        m->GetElementCenter(e, c);
        if (c(0) < 0.4)
        {
          refs.Append(e);
        }
      }
      m->GeneralRefinement(refs, 1, 1);
    }
    Check("read-refined", true);
    mesh = std::move(saved);
    Check("refined", true);
    if (Mpi::Root(comm))
    {
      CompareSavedMeshes(temp_dir / "refined-gathered.mesh",
                         temp_dir / "read-refined-gathered.mesh", counts);
    }
  }

  // Read the last file back (as a restart does), scatter it, and save it again through
  // RebalanceMesh.
  mfem::named_ifgzstream fi(file.string());
  mfem::Mesh reread(fi);
  mesh::NondimensionalizeMesh(reread, L);
  std::vector<int> part(reread.GetNE());
  for (int e = 0; e < reread.GetNE(); e++)
  {
    part[e] = (e % 3 == 0) ? 0 : (e * 7919) % P;
  }
  mesh = std::make_unique<mfem::ParMesh>(comm, reread, part.data());
  Check("reread", false);
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh parallel NC save matches the gathered save",
                 "[basesolver][Serial][Parallel]")
{
  const auto type =
      GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON, mfem::Element::WEDGE,
               mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  const int order = GENERATE(0, 1, 2, 3);
  TestParallelNCSave(MPI_COMM_WORLD, temp_dir, type, order);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh parallel NC save matches the gathered save, more processes",
                 "[basesolver][Parallel][Regression]")
{
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  TestParallelNCSave(MPI_COMM_WORLD, temp_dir, type, 2);
}
