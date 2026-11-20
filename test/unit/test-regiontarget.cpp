// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "driver.hpp"
#include "surfaceresponse-fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

// Ruling 2251 (OI-REGION-TARGET-CORRECTED-NOT-FILTERED): the corrected columns of a target
// entry with a Region count s_p of each patch's fabricated energy f_p, s_p = the share of
// the patch's own within-R raw energy inside the region (the same 0 / 1 test at the same
// face quadrature points as the raw terms). On the convex island of the unit box (x, z in
// [0.25, 0.75] on y = 0.5, h = 0.125, R = 0.2; the default patch construction, whose
// patches carry raw claims) the half-space box x < 0.375 holds the patches of the edge x =
// 0.25 whole, cuts the corner patches at x = 0.25 (their claims run R = 0.2 along the edges
// z = 0.25 / 0.75, across the wall) and leaves out those of the edge x = 0.75; its
// complement x >= 0.375 tiles the island with it (half-open; the wall is a mesh plane, no
// quadrature point on it).

namespace palace
{
namespace
{

namespace fs = std::filesystem;
using json = nlohmann::json;

constexpr double kBoxWall = 0.375;

// The Region of the box x < kBoxWall (lower = true) or of its complement.
json HalfSpaceRegion(bool lower)
{
  return {{"BoxMin", {lower ? -10.0 : kBoxWall, -10.0, -10.0}},
          {"BoxMax", {lower ? kBoxWall : 10.0, 10.0, 10.0}}};
}

// The values of `header` on the rows where every (header, value) of `keys` holds (integer
// key columns; a key value below 0 matches every row).
std::vector<double> Select(const Table &table, const std::string &header,
                           const std::vector<std::pair<std::string, int>> &keys)
{
  const Column &column = test::ColumnByHeader(table, header);
  std::vector<const Column *> key_columns;
  for (const auto &key : keys)
  {
    key_columns.push_back(&test::ColumnByHeader(table, key.first));
  }
  std::vector<double> values;
  for (std::size_t i = 0; i < column.data.size(); i++)
  {
    bool match = true;
    for (std::size_t k = 0; k < keys.size(); k++)
    {
      match = match && (keys[k].second < 0 ||
                        std::lround(key_columns[k]->data[i]) == keys[k].second);
    }
    if (match)
    {
      values.push_back(column.data[i]);
    }
  }
  return values;
}

double Single(const Table &table, const std::string &header,
              const std::vector<std::pair<std::string, int>> &keys)
{
  const auto values = Select(table, header, keys);
  REQUIRE(values.size() == 1);
  return values[0];
}

// sum_p s_p f_p of the patch table rows (source, evaluation) of target 4 (source below 0:
// the table's first source), with s_p from the region-share rows (excitation, share
// evaluation); without a share file s_p = 1 (the window-wide sum). A patch the share file
// does not list has s_p = 0.
double RecomposedFabricated(const fs::path &dir, int source, int evaluation, bool region,
                            int excitation, int share_evaluation)
{
  std::map<long, double> shares;
  if (region)
  {
    const Table table = test::LoadCsv(dir / "surface-response-region-share.csv");
    const std::vector<std::pair<std::string, int>> keys = {
        {"excitation", excitation}, {"evaluation", share_evaluation}, {"target", 4}};
    const auto patches = Select(table, "patch", keys);
    const auto values = Select(table, "share", keys);
    REQUIRE(!patches.empty());
    for (std::size_t i = 0; i < patches.size(); i++)
    {
      CHECK(values[i] >= 0.0);
      CHECK(values[i] <= 1.0);
      shares[std::lround(patches[i])] = values[i];
    }
  }
  const Table table = test::LoadCsv(dir / "surface-response-patch-energy.csv");
  if (source < 0)
  {
    REQUIRE(table.n_rows() > 0);
    source = static_cast<int>(std::lround(test::ColumnByHeader(table, "source").data[0]));
  }
  const std::vector<std::pair<std::string, int>> keys = {{"source", source},
                                                         {"evaluation", evaluation}};
  const auto patches = Select(table, "patch", keys);
  const auto energies = Select(table, "fabricated surface energy[4] (J)", keys);
  REQUIRE(!patches.empty());
  double sum = 0.0;
  for (std::size_t i = 0; i < patches.size(); i++)
  {
    const auto share = shares.find(std::lround(patches[i]));
    sum += (!region ? 1.0 : share != shares.end() ? share->second : 0.0) * energies[i];
  }
  return sum;
}

// The convex library in the default patch construction (as the corner-arm trim case
// writes it: every model on SA coupon 1, the isolated edge on a 4-point basis square);
// written by the root rank.
fs::path WriteDefaultConstructionLibrary(const fs::path &convex_library,
                                         const fs::path &dir)
{
  const fs::path basis_path = dir / "region-target-basis-points.csv";
  const fs::path library_path = dir / "fabrication-process-region-target-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream basis(basis_path);
    basis << "x,y,z\n-0.2,-0.2,0.0\n0.2,-0.2,0.0\n0.2,0.2,0.0\n-0.2,0.2,0.0\n";
    std::ifstream input(convex_library);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-region-target-3d";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      if (model["Topology"] == "IsolatedEdge")
      {
        model["BasisPoints"] = basis_path.string();
      }
    }
    std::ofstream output(library_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return library_path;
}

// The default-construction correction of a fixture config (the corner-arm trim case's; the
// surface mortar on the electrostatic path only: Maxwell refuses it).
void UseDefaultConstruction(json &correction, const fs::path &library, bool mortar)
{
  correction.erase("PatchConstruction");
  correction["Library"] = library.string();
  if (mortar)
  {
    correction["TraceCoupling"] = "SurfaceMortar";
    correction["MortarOversampling"] = 2;
  }
}

json RegionRecord(const fs::path &dir)
{
  std::ifstream input(dir / "palace.json");
  const json meta = json::parse(input);
  const auto &response = meta.at("SurfaceResponse");
  if (!response.contains("Diagnostics") ||
      !response["Diagnostics"].contains("RegionTargets"))
  {
    return json();
  }
  return response["Diagnostics"]["RegionTargets"];
}

}  // namespace

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Region target corrected energies count each patch's share inside the box "
                 "(electrostatic)",
                 "[regiontarget][electrostaticsolver][surfaceresponseoperator][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "region-target-island.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  json config = ConvexIslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  auto &correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  UseDefaultConstruction(
      correction, WriteDefaultConstructionLibrary(convex_library_3d_path, temp.temp_dir),
      true);
  correction["CorrectionMode"] = "Both";
  correction["PatchEnergy"] = true;
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 400}};
  auto Run = [&](const std::string &stem, const json &region)
  {
    json run = config;
    if (!region.is_null())
    {
      run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Region"] = region;
    }
    const fs::path output = temp.temp_dir / stem;
    run["Problem"]["Output"] = output.string();
    test::RunElectrostatic(run);
    Mpi::Barrier(Mpi::World());
    return output;
  };
  const fs::path window = Run("region-window", json());
  const fs::path lower = Run("region-lower", HalfSpaceRegion(true));
  const fs::path upper = Run("region-upper", HalfSpaceRegion(false));
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  auto Corrected = [](const fs::path &dir, const std::string &column)
  {
    return Single(test::LoadCsv(dir / "surface-Q-corrected.csv"),
                  fmt::format("E_surf {}[4] (J)", column), {{"i", 1}});
  };
  auto Outside = [](const fs::path &dir)
  {
    return Single(test::LoadCsv(dir / "surface-Q-edge.csv"), "E_out (J)",
                  {{"i", 1}, {"interface", 4}});
  };
  // An entry without a Region: no record, no share file, the window rule (no kept term on
  // this island: E_ft = E_out + sum_p f_p).
  CHECK(RegionRecord(window).is_null());
  CHECK_FALSE(fs::exists(window / "surface-response-region-share.csv"));
  for (const fs::path &dir : {window, lower, upper})
  {
    for (const char *kept : {"surface-response-uncovered-energy.csv",
                             "surface-response-domain-boundary-energy.csv",
                             "surface-response-conductor-consistency-energy.csv"})
    {
      CHECK_FALSE(fs::exists(dir / kept));
    }
  }
  const double window_ft = RecomposedFabricated(window, 1, 0, false, 0, 0);
  CHECK(Corrected(window, "postprocessed fixed-trace") ==
        Catch::Approx(Outside(window) + window_ft).epsilon(1.0e-10));
  // The patch table (f_p and the excluded set) is the window's on every run.
  CHECK(test::ReadFile(lower / "surface-response-patch-energy.csv") ==
        test::ReadFile(window / "surface-response-patch-energy.csv"));
  CHECK(test::ReadFile(upper / "surface-response-patch-energy.csv") ==
        test::ReadFile(window / "surface-response-patch-energy.csv"));
  for (const fs::path &dir : {lower, upper})
  {
    INFO("run " << dir.filename().string());
    // The recomposition from the open tables: E_out of the box + sum_p s_p f_p.
    const double box_ft = RecomposedFabricated(dir, 1, 0, true, 1, 0);
    const double box_ff = RecomposedFabricated(dir, 1, 1, true, 1, 0);
    CHECK(Corrected(dir, "postprocessed fixed-trace") ==
          Catch::Approx(Outside(dir) + box_ft).epsilon(1.0e-10));
    CHECK(Corrected(dir, "postprocessed fixed-flux") ==
          Catch::Approx(Outside(dir) + box_ff).epsilon(1.0e-10));
    // The record: the raw and the self-consistent evaluation of terminal 1; the box holds,
    // cuts and leaves out patches; its box / window ratio is the tables' one.
    const json record = RegionRecord(dir);
    REQUIRE(!record.is_null());
    REQUIRE(record.at("Evaluations").size() == 2);
    CHECK(record["Evaluations"][0]["Evaluation"] == "Raw");
    CHECK(record["Evaluations"][1]["Evaluation"] == "SelfConsistent");
    const auto &target = record["Evaluations"][0]["Targets"][0];
    CHECK(record["Evaluations"][0]["Excitation"].get<int>() == 1);
    CHECK(target["Index"].get<int>() == 4);
    CHECK(target["PatchesInside"].get<int>() > 0);
    CHECK(target["PatchesCut"].get<int>() > 0);
    CHECK(target["PatchesOutside"].get<int>() > 0);
    CHECK(target["CutShare"].get<double>() > 0.0);
    CHECK(target["CutShare"].get<double>() < 1.0);
    CHECK(target["FabricatedEnergyRegion"].get<double>() /
              target["FabricatedEnergyWindow"].get<double>() ==
          Catch::Approx(box_ft / window_ft).epsilon(1.0e-10));
  }
  // Boxes that tile the island add up to the entry without a Region (every corrected
  // column; the raw participation too).
  for (const char *column :
       {"postprocessed fixed-trace", "postprocessed fixed-flux", "corrected"})
  {
    INFO("column " << column);
    CHECK(Corrected(lower, column) + Corrected(upper, column) ==
          Catch::Approx(Corrected(window, column)).epsilon(1.0e-10));
  }
  auto Raw = [](const fs::path &dir)
  { return Single(test::LoadCsv(dir / "surface-Q.csv"), "p_surf[4]", {{"i", 1}}); };
  CHECK(Raw(lower) + Raw(upper) == Catch::Approx(Raw(window)).epsilon(1.0e-10));
#endif
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Region target corrected energies count each patch's share inside the box "
                 "(eigenmode)",
                 "[regiontarget][eigensolver][surfaceresponseoperator][3d][Serial]"
                 "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  const fs::path mesh_path = temp.temp_dir / "region-target-maxwell-island.mesh";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    std::ofstream output(mesh_path);
    serial.Print(output);
  }
  Mpi::Barrier(comm);
  json config = ConvexMaxwellIslandConfig();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  UseDefaultConstruction(
      config["Solver"]["SurfaceResponseCorrection"],
      WriteDefaultConstructionLibrary(convex_library_3d_path, temp.temp_dir), false);
  config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
  config["Solver"]["SurfaceResponseCorrection"]["PatchEnergy"] = true;
  // The lowest mode of the 1 m PEC box.
  config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 1}};
  auto Run = [&](const std::string &stem, const json &region)
  {
    json run = config;
    if (!region.is_null())
    {
      run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Region"] = region;
    }
    const fs::path output = temp.temp_dir / stem;
    run["Problem"]["Output"] = output.string();
    IoData iodata(run, false);
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    palace::Run(iodata, comm, utils::ConfigureOmp(), /*git_tag=*/nullptr);
    Mpi::Barrier(comm);
    return output;
  };
  const fs::path window = Run("region-maxwell-window", json());
  const fs::path lower = Run("region-maxwell-lower", HalfSpaceRegion(true));
  const fs::path upper = Run("region-maxwell-upper", HalfSpaceRegion(false));
  if (!Mpi::Root(comm))
  {
    return;
  }
  auto Corrected = [](const fs::path &dir, const std::string &column)
  {
    return Single(test::LoadCsv(dir / "surface-Q-corrected.csv"),
                  fmt::format("E_surf {} (J)", column), {{"m", 1}, {"interface", 4}});
  };
  auto Outside = [](const fs::path &dir)
  {
    return Single(test::LoadCsv(dir / "surface-Q-edge.csv"), "E_out (J)",
                  {{"m", 1}, {"interface", 4}});
  };
  CHECK(RegionRecord(window).is_null());
  CHECK_FALSE(fs::exists(window / "surface-response-region-share.csv"));
  const double window_ft = RecomposedFabricated(window, -1, 0, false, 0, 0);
  // surface-Q-corrected.csv of the Maxwell path prints its energies in J with 16 decimals:
  // its resolution is 1e-16 J, not a relative one (the island's energies are ~1e-10 J).
  constexpr double kPrinted = 1.0e-15;
  CHECK(Corrected(window, "corrected") ==
        Catch::Approx(Outside(window) + window_ft).epsilon(1.0e-10).margin(kPrinted));
  CHECK(test::ReadFile(lower / "surface-response-patch-energy.csv") ==
        test::ReadFile(window / "surface-response-patch-energy.csv"));
  CHECK(test::ReadFile(upper / "surface-response-patch-energy.csv") ==
        test::ReadFile(window / "surface-response-patch-energy.csv"));
  for (const fs::path &dir : {lower, upper})
  {
    INFO("run " << dir.filename().string());
    // Mode 1 is the construction's first measurement: excitation 0 of the share rows.
    const double box_ft = RecomposedFabricated(dir, -1, 0, true, 0, 0);
    const double box_ff = RecomposedFabricated(dir, -1, 1, true, 0, 0);
    CHECK(Corrected(dir, "corrected") ==
          Catch::Approx(Outside(dir) + box_ft).epsilon(1.0e-10).margin(kPrinted));
    CHECK(Corrected(dir, "corrected fixed-flux") ==
          Catch::Approx(Outside(dir) + box_ff).epsilon(1.0e-10).margin(kPrinted));
    // Every measured mode has a raw and a self-consistent record (the solver may measure
    // more converged modes than N); mode 1 = excitation 0.
    const json record = RegionRecord(dir);
    REQUIRE(!record.is_null());
    REQUIRE(record.at("Evaluations").size() >= 2);
    REQUIRE(record["Evaluations"].size() % 2 == 0);
    CHECK(record["Evaluations"][0]["Evaluation"] == "Raw");
    CHECK(record["Evaluations"][1]["Evaluation"] == "SelfConsistent");
    CHECK(record["Evaluations"][1]["Excitation"].get<int>() == 0);
    const auto &target = record["Evaluations"][0]["Targets"][0];
    CHECK(record["Evaluations"][0]["Excitation"].get<int>() == 0);
    CHECK(target["Index"].get<int>() == 4);
    CHECK(target["PatchesInside"].get<int>() > 0);
    CHECK(target["PatchesCut"].get<int>() > 0);
    CHECK(target["PatchesOutside"].get<int>() > 0);
    CHECK(target["FabricatedEnergyRegion"].get<double>() /
              target["FabricatedEnergyWindow"].get<double>() ==
          Catch::Approx(box_ft / window_ft).epsilon(1.0e-10));
  }
  for (const char *column : {"corrected", "corrected fixed-flux", "self-consistent"})
  {
    INFO("column " << column);
    CHECK(Corrected(lower, column) + Corrected(upper, column) ==
          Catch::Approx(Corrected(window, column)).epsilon(1.0e-10).margin(kPrinted));
  }
  CHECK(Corrected(lower, "raw") + Corrected(upper, "raw") ==
        Catch::Approx(Corrected(window, "raw")).epsilon(1.0e-10).margin(kPrinted));
#endif
}

}  // namespace palace
