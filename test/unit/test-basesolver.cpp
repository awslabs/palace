// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "drivers/basesolver.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/dorfler.hpp"
#include "utils/filesystem.hpp"

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
                 "SaveIteration handles dirty output directory from previous run",
                 "[basesolver][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Simulate a previous run that left files in the iteration folder.
  CreateFile(temp_dir / "palace.json", "{}");
  fs::create_directories(temp_dir / "iteration1" / "paraview" / "driven");
  CreateFile(temp_dir / "iteration1" / "paraview" / "driven" / "old.vtu", "old");
  CreateFile(temp_dir / "iteration1" / "domain-E.csv", "old_data");

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
