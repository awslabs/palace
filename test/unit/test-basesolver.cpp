// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/dorfler.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using json = nlohmann::json;

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

// P3-24 (SL5-RESOLUTION s10; the P3 stub of
// resolutions/SL5/palace__drivers__basesolver.cpp.md): WriteMetadata = upstream's
// palace.json.tmp + rename (#906) with ours' checked streamed write BEFORE the rename. A
// short write (here: palace.json.tmp is a symbolic link to /dev/full, so the flush fails
// with ENOSPC) aborts and leaves the previous palace.json byte-identical; the partial
// temporary file is left behind, never renamed into place. The fault is injected THROUGH
// the symbolic link, so the case relies on WriteMetadata opening its temporary file through
// an existing link: if that writer is hardened to remove links before writing (as the other
// outputs, SD-302), the injection must change (P3-C review m-5). REQUIRE (not CHECK) stops
// the case on the defect it guards: a renamed link to /dev/full is never read back (P3-C
// review m-7: the read would not end).
TEST_CASE_METHOD(palace::test::SharedTempDir, "WriteMetadata short write keeps palace.json",
                 "[basesolver][Serial]")
{
#if !defined(MFEM_USE_EXCEPTIONS)
  SKIP("The short-write abort is observed as an exception (MFEM_USE_EXCEPTIONS)");
#else
  if (!fs::exists("/dev/full"))
  {
    SKIP("No /dev/full on this host");
  }
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
  REQUIRE(counts.valid);

  // The BaseSolver constructor writes the initial palace.json; SaveMetadata rewrites it.
  TestSolver solver(iodata, Mpi::Root(comm));
  const fs::path path = temp_dir / "palace.json";
  const fs::path tmp_path = temp_dir / "palace.json.tmp";
  REQUIRE(fs::is_regular_file(path));
  const std::string before = ReadFile(path);
  REQUIRE_FALSE(before.empty());
  fs::create_symlink("/dev/full", tmp_path);
  REQUIRE_THROWS_WITH(solver.SaveMetadata(counts),
                      Catch::Matchers::ContainsSubstring("a short write"));
  REQUIRE_FALSE(fs::is_symlink(path));
  CHECK(ReadFile(path) == before);
  CHECK(fs::is_symlink(tmp_path));
  fs::remove(tmp_path);

  // The control: the same call without the failing stream writes the block.
  solver.SaveMetadata(counts);
  CHECK(nlohmann::json::parse(ReadFile(path)).contains("SavedAdaptedMesh"));
  CHECK_FALSE(fs::exists(tmp_path));
#endif
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

TEST_CASE("PredictRefinedSize scales the last solved size by the element growth",
          "[basesolver][Serial]")
{
  // Stage-2 S5 thin P4 (order 4) AMR cycles: the c9 -> c10 refinement grew the mesh
  // 9,034,822 -> 25,852,252 elements after a c9 solve of 67,992,598 unknowns; c10 assembled
  // 196,931,656 unknowns (the prediction reads 1.2% low on this graded nonconforming mesh).
  CHECK(PredictRefinedSize(67992598, 9034822, 25852252) == 194554113);
  // c8 -> c9 of the same run: 25,343,372 unknowns, 3,435,515 -> 9,034,822 elements
  // (actual c9: 67,992,598).
  CHECK(PredictRefinedSize(25343372, 3435515, 9034822) == 66648772);
  // The ratio is rounded up and an empty refinement step leaves the size unchanged.
  CHECK(PredictRefinedSize(10, 3, 4) == 14);
  CHECK(PredictRefinedSize(10, 3, 3) == 10);
  CHECK(PredictRefinedSize(0, 3, 4) == 0);
}

TEST_CASE("ExceedsPredictedMaxSize stops before an oversized solve", "[basesolver][Serial]")
{
  config::RefinementData refinement;
  refinement.max_size = 100000000;

  SECTION("Predicted rule (default): the S5 P4 c9 -> c10 step stops, c8 -> c9 does not")
  {
    CHECK(refinement.max_size_predicted);
    // The last solved size (67,992,598) is below MaxSize, so the legacy loop-head check
    // would have let the 196,931,656-unknown c10 solve run (x2.9 MaxSize).
    CHECK(ExceedsPredictedMaxSize(refinement, 67992598, 9034822, 25852252));
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 25343372, 3435515, 9034822));
    // A prediction exactly at MaxSize is allowed; one above is not.
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 50000000, 1000, 2000));
    CHECK(ExceedsPredictedMaxSize(refinement, 50000000, 1000, 2001));
  }

  SECTION("Legacy rule (MaxSizePredicted = false) never stops before the solve")
  {
    refinement.max_size_predicted = false;
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 67992598, 9034822, 25852252));
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 99999999, 1, 1000));
  }

  SECTION("No MaxSize (<= 0) means no size constraint")
  {
    refinement.max_size = 0;
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 67992598, 9034822, 25852252));
    refinement.max_size = -1;
    CHECK_FALSE(ExceedsPredictedMaxSize(refinement, 67992598, 9034822, 25852252));
  }
}

TEST_CASE("MarkEligibleElements without a floor marks the Dörfler set",
          "[basesolver][Serial][Parallel]")
{
  // Every element eligible: the marked set, the threshold and the marked error are those
  // of the unrestricted Dörfler marking on the same indicators (the default code path).
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm);
  Vector e(6);
  for (int i = 0; i < e.Size(); i++)
  {
    e[i] = 0.5 + i + 0.1 * rank;
  }
  const std::vector<bool> eligible(e.Size(), true);
  const double fraction = 0.7;
  const auto marking = MarkEligibleElements(comm, e, eligible, fraction);
  const auto [threshold, marked_error] = utils::ComputeDorflerThreshold(comm, e, fraction);
  CHECK(marking.markable);
  CHECK(marking.threshold == threshold);
  CHECK(marking.marked_error == marked_error);
  CHECK(marking.ineligible_elements == 0);
  CHECK(marking.ineligible_error == 0.0);
  mfem::Array<int> expected;
  for (int i = 0; i < e.Size(); i++)
  {
    if (e[i] >= threshold)
    {
      expected.Append(i);
    }
  }
  CHECK(marking.marked_elements == expected);
  CHECK(marking.marked_elements.Size() > 0);
}

TEST_CASE("MarkEligibleElements excludes the ineligible elements from the threshold",
          "[basesolver][Serial][Parallel]")
{
  // The two largest indicators sit on ineligible elements: they are never marked, do not
  // enter the threshold, and are reported with their share of the squared indicator.
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  Vector e(6);
  e[0] = 10.0;
  e[1] = 9.0;
  e[2] = 1.0;
  e[3] = 2.0;
  e[4] = 3.0;
  e[5] = 4.0;
  const std::vector<bool> eligible = {false, false, true, true, true, true};
  const double fraction = 0.7;
  const auto marking = MarkEligibleElements(comm, e, eligible, fraction);
  CHECK(marking.markable);

  Vector eligible_e(4);
  eligible_e[0] = 1.0;
  eligible_e[1] = 2.0;
  eligible_e[2] = 3.0;
  eligible_e[3] = 4.0;
  const auto [threshold, marked_error] =
      utils::ComputeDorflerThreshold(comm, eligible_e, fraction);
  CHECK(marking.threshold == threshold);
  CHECK(marking.marked_error == marked_error);
  CHECK(marking.ineligible_elements == 2 * size);
  CHECK_THAT(marking.ineligible_error,
             Catch::Matchers::WithinRel(
                 (100.0 + 81.0) / (100.0 + 81.0 + 1.0 + 4.0 + 9.0 + 16.0), 1.0e-14));
  // 70% of the eligible squared error 30 needs {4, 3} (25) : {4} alone is 16 < 21.
  mfem::Array<int> expected({4, 5});
  CHECK(marking.marked_elements == expected);
  for (int i : marking.marked_elements)
  {
    CHECK(eligible[i]);
  }
}

TEST_CASE("MarkEligibleElements with an empty eligible set marks nothing",
          "[basesolver][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);
  Vector e(3);
  e[0] = 1.0;
  e[1] = 2.0;
  e[2] = 3.0;

  SECTION("No eligible element on any process")
  {
    const std::vector<bool> eligible(e.Size(), false);
    const auto marking = MarkEligibleElements(comm, e, eligible, 0.7);
    CHECK(!marking.markable);
    CHECK(marking.marked_elements.Size() == 0);
    CHECK(marking.ineligible_elements == 3 * size);
    CHECK(marking.ineligible_error == 1.0);
    CHECK(marking.threshold == 0.0);
    CHECK(marking.marked_error == 0.0);
  }

  SECTION("Eligible elements with a zero indicator only")
  {
    // All the error sits on the ineligible elements: there is no Dörfler threshold over
    // the eligible set (0 / 0), so the marking reads like the empty eligible set.
    e[0] = 0.0;
    const std::vector<bool> eligible = {true, false, false};
    const auto marking = MarkEligibleElements(comm, e, eligible, 0.7);
    CHECK(!marking.markable);
    CHECK(marking.marked_elements.Size() == 0);
    CHECK(marking.ineligible_elements == 2 * size);
    CHECK(marking.ineligible_error == 1.0);
    CHECK(marking.threshold == 0.0);
    CHECK(marking.marked_error == 0.0);
  }

  SECTION("Eligible elements on the root process only")
  {
    // The marking is collective: the other processes take part in the threshold and mark
    // nothing; the root marks eligible elements only, among them the largest indicator,
    // and the marked share of the eligible error (4 + 9) reaches the fraction.
    const bool root = Mpi::Root(comm);
    const std::vector<bool> eligible = {false, root, root};
    const auto marking = MarkEligibleElements(comm, e, eligible, 0.7);
    CHECK(marking.markable);
    CHECK(marking.ineligible_elements == 3 * size - 2);
    CHECK(marking.marked_error >= 0.7);
    CHECK(marking.marked_elements.Size() == (root ? marking.marked_elements.Size() : 0));
    if (root)
    {
      REQUIRE(marking.marked_elements.Size() >= 1);
      CHECK(marking.marked_elements[marking.marked_elements.Size() - 1] == 2);
      for (int i : marking.marked_elements)
      {
        CHECK(eligible[i]);
      }
    }
  }
}

TEST_CASE("ElementsEligibleForMarking reads the element size of the mesh statistics",
          "[basesolver][Serial][Parallel]")
{
  // Two 1.5 x 1 x 1 hexahedra; one is refined once, so its 8 children have half its
  // volume-based size h = (1.5)^(1/3).
  auto serial_mesh =
      mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::HEXAHEDRON, 3.0, 1.0, 1.0);
  serial_mesh.EnsureNCMesh(true);
  mfem::Array<int> refine_first(1);
  refine_first[0] = 0;
  serial_mesh.GeneralRefinement(refine_first, 1);
  REQUIRE(serial_mesh.GetNE() == 9);
  mfem::ParMesh mesh(Mpi::World(), serial_mesh);
  const double h_parent = std::cbrt(1.5), h_child = 0.5 * h_parent;

  // The same measure as the h of the mesh statistics (Mesh::GetCharacteristics).
  double h_min, h_max, kappa_min, kappa_max;
  mesh.GetCharacteristics(h_min, h_max, kappa_min, kappa_max);
  double h_min_elements = mfem::infinity();
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    h_min_elements = std::min(h_min_elements, mesh.GetElementSize(i, 0));
  }
  Mpi::GlobalMin(1, &h_min_elements, mesh.GetComm());
  CHECK(h_min == h_min_elements);
  CHECK_THAT(h_min, Catch::Matchers::WithinRel(h_child, 1.0e-12));
  CHECK_THAT(h_max, Catch::Matchers::WithinRel(h_parent, 1.0e-12));

  // Floor at h_child / 2: every element is eligible (h >= 2 x floor = h_child).
  for (bool eligible : ElementsEligibleForMarking(mesh, 0.5 * h_child * (1.0 - 1.0e-12)))
  {
    CHECK(eligible);
  }
  // Floor just above h_child / 2: the children are ineligible, the parent eligible.
  {
    const auto eligible = ElementsEligibleForMarking(mesh, 0.5 * h_child * (1.0 + 1.0e-12));
    int num_eligible = 0;
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      CHECK(eligible[i] == (mesh.GetElementSize(i, 0) > h_child * (1.0 + 1.0e-13)));
      num_eligible += eligible[i];
    }
    Mpi::GlobalSum(1, &num_eligible, mesh.GetComm());
    CHECK(num_eligible == 1);
  }

  // Elements below a floor (the closure count after a refinement step).
  CHECK(CountElementsBelowSize(mesh, h_child * (1.0 - 1.0e-12)) == 0);
  CHECK(CountElementsBelowSize(mesh, h_child * (1.0 + 1.0e-12)) == 8);
  CHECK(CountElementsBelowSize(mesh, h_parent * (1.0 + 1.0e-12)) == 9);
}

// OI-ADAPT-SAVEMESH-OVERWRITE (idfix DESIGN section 7; decision 1072 Q1): the mesh writer
// never writes THROUGH the symlink SaveIteration left.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RemoveSymlinkBeforeWrite leaves the archived target untouched",
                 "[basesolver][savemesh][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  CreateFile(temp_dir / "palace.json", "{}");
  CreateFile(temp_dir / "box.mesh", "iteration1 mesh");
  SaveIteration(comm, temp_dir, 1, 2);
  REQUIRE(fs::is_symlink(temp_dir / "box.mesh"));
  if (Mpi::Root(comm))
  {
    CHECK(mesh::RemoveSymlinkBeforeWrite(temp_dir / "box.mesh"));
    CHECK(!fs::exists(fs::symlink_status(temp_dir / "box.mesh")));
    CreateFile(temp_dir / "box.mesh", "iteration2 mesh");
    // The top-level file is real, the archived iteration keeps its own bytes.
    CHECK(!fs::is_symlink(temp_dir / "box.mesh"));
    CHECK(ReadFile(temp_dir / "box.mesh") == "iteration2 mesh");
    CHECK(ReadFile(temp_dir / "iteration01" / "box.mesh") == "iteration1 mesh");
    // A real file or a missing path is left alone.
    CHECK(!mesh::RemoveSymlinkBeforeWrite(temp_dir / "box.mesh"));
    CHECK(ReadFile(temp_dir / "box.mesh") == "iteration2 mesh");
    CHECK(!mesh::RemoveSymlinkBeforeWrite(temp_dir / "absent.mesh"));
  }
  Mpi::Barrier(comm);
}

// OI-MAGNETIC-CSV-WRITE-THROUGH-LINK (fix round P4 re-freeze; decision 1368): the magnetic
// table writer (PostOperator::PrintSurfaceResponseMagnetic, "rewritten whole at every
// print") opens surface-response-magnetic.csv for truncation; after SaveIteration that path
// is a symlink into iteration<N>/, so the writer removes the link first (the same form as
// the mesh writer). The writer's open sequence, twice across a SaveIteration: the archived
// first print keeps its bytes and the top-level file is real.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "surface-response-magnetic.csv is rewritten as a real file after "
                 "SaveIteration (OI-MAGNETIC-CSV-WRITE-THROUGH-LINK)",
                 "[basesolver][savemesh][magnetic][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const fs::path table = temp_dir / "surface-response-magnetic.csv";
  auto WriteWhole = [&](const std::string &text)
  {
    // The writer's sequence: remove a SaveIteration symlink, then truncate a real file.
    mesh::RemoveSymlinkBeforeWrite(table);
    std::ofstream output(table, std::ios::trunc);
    REQUIRE(output);
    output << text;
  };
  CreateFile(temp_dir / "palace.json", "{}");
  if (Mpi::Root(comm))
  {
    WriteWhole("kind,model\nmodel,0\n");
  }
  Mpi::Barrier(comm);
  SaveIteration(comm, temp_dir, 1, 2);
  REQUIRE(fs::is_symlink(fs::symlink_status(table)));
  if (Mpi::Root(comm))
  {
    WriteWhole("kind,model\nmodel,0\ngate,0\n");
    CHECK(!fs::is_symlink(fs::symlink_status(table)));
    CHECK(fs::is_regular_file(fs::symlink_status(table)));
    CHECK(ReadFile(table) == "kind,model\nmodel,0\ngate,0\n");
    CHECK(ReadFile(temp_dir / "iteration01" / "surface-response-magnetic.csv") ==
          "kind,model\nmodel,0\n");
  }
  Mpi::Barrier(comm);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RebalanceMesh with SaveAdaptMesh keeps the archived iteration mesh",
                 "[basesolver][savemesh][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const fs::path output = temp_dir / "out";
  if (Mpi::Root(comm))
  {
    fs::create_directories(output);
  }
  Mpi::Barrier(comm);
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", output.string()}}},
                 {"Model",
                  {{"Mesh", (temp_dir / "box.mesh").string()},
                   {"Refinement", {{"SaveAdaptMesh", true}}}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", {{"Ground", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}}},
                 {"Solver", {{"Order", 1}, {"Electrostatic", json::object()}}}};
  IoData iodata(config, false);
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON));
  serial_mesh->EnsureNCMesh(true);
  auto mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  serial_mesh.reset();

  // The saved mesh takes upstream's name rule (SD-302): gzip-compressed ".meshgz" when MFEM
  // is built with zlib, ".mesh" otherwise.
#if defined(MFEM_USE_ZLIB)
  const std::string saved_mesh = "box.meshgz";
#else
  const std::string saved_mesh = "box.mesh";
#endif

  // Iteration 1: the writer saves out/<saved_mesh>, SaveIteration archives it and leaves
  // the symlink out/<saved_mesh> -> iteration01/<saved_mesh>.
  mesh::RebalanceMesh(iodata, mesh);
  REQUIRE(fs::is_regular_file(fs::symlink_status(output / saved_mesh)));
  if (Mpi::Root(comm))
  {
    CreateFile(output / "palace.json", "{}");
  }
  SaveIteration(comm, output, 1, 2);
  REQUIRE(fs::is_symlink(output / saved_mesh));
  const std::string archived = ReadFile(output / "iteration01" / saved_mesh);
  REQUIRE(!archived.empty());

  // Iteration 2: the mesh is refined and saved again. The archived iteration-1 mesh is
  // byte-identical before and after; the top-level file is a real file with the new mesh.
  mesh->UniformRefinement();
  mesh::RebalanceMesh(iodata, mesh);
  Mpi::Barrier(comm);
  CHECK(!fs::is_symlink(output / saved_mesh));
  CHECK(fs::is_regular_file(output / saved_mesh));
  CHECK(ReadFile(output / "iteration01" / saved_mesh) == archived);
  CHECK(ReadFile(output / saved_mesh) != archived);
  Mpi::Barrier(comm);
}
