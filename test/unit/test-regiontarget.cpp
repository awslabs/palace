// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
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
// patches carry raw claims: 4 corners whose claims run R along both edges, 4 cells of 0.025
// per edge in [0.45, 0.55]) the half-space box x < 0.375 cuts the cells of the edge x =
// 0.25 that have points (their within-R points reach x = 0.45) and the corner patches at x
// = 0.25 (their
// claims run R along the edges z = 0.25 / 0.75, across the wall) and leaves out those of
// the edge x = 0.75; its complement x >= 0.375 tiles the island with it (half-open; the
// wall is a mesh plane, no quadrature point on it). The wall x = 0.5 lies beyond R of both
// edges x = 0.25 / 0.75: each of them counts whole on one side of it.

namespace palace
{
namespace
{

namespace fs = std::filesystem;
using json = nlohmann::json;

constexpr double kBoxWall = 0.375;
constexpr double kMidWall = 0.5;

// The Region of the box x < wall (lower = true) or of its complement.
json HalfSpaceRegion(bool lower, double wall = kBoxWall)
{
  return {{"BoxMin", {lower ? -10.0 : wall, -10.0, -10.0}},
          {"BoxMax", {lower ? wall : 10.0, 10.0, 10.0}}};
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

// The share rows of target 4 at (excitation, evaluation), each checked against ruling
// 2251's tiers: the tier is the first of (1) the points no kept class holds, (2) every
// point of the claim by energy, (3) by area that has a denominator, (4) none (the reference
// point; a tier-4 patch is listed only when that point lies in the box: share 1), and the
// share is that tier's ratio. Returns the share per listed patch (1-based).
std::map<long, double> CheckedShareRows(const fs::path &dir, int excitation, int evaluation)
{
  const Table table = test::LoadCsv(dir / "surface-response-region-share.csv");
  const std::vector<std::pair<std::string, int>> keys = {
      {"excitation", excitation}, {"evaluation", evaluation}, {"target", 4}};
  const auto patches = Select(table, "patch", keys);
  const auto tiers = Select(table, "tier", keys);
  const auto shares = Select(table, "share", keys);
  const auto open = Select(table, "open_energy (J)", keys);
  const auto open_in = Select(table, "open_energy_in (J)", keys);
  const auto energy = Select(table, "energy (J)", keys);
  const auto energy_in = Select(table, "energy_in (J)", keys);
  const auto area = Select(table, "area (m^2)", keys);
  const auto area_in = Select(table, "area_in (m^2)", keys);
  REQUIRE(!patches.empty());
  std::map<long, double> result;
  for (std::size_t i = 0; i < patches.size(); i++)
  {
    INFO("patch " << std::lround(patches[i]) << ", tier " << std::lround(tiers[i]));
    const auto Ratio = [](double part, double whole)
    { return Catch::Approx(part / whole).epsilon(1.0e-12).margin(1.0e-14); };
    if (open[i] != 0.0)
    {
      CHECK(std::lround(tiers[i]) == 1);
      CHECK(shares[i] == Ratio(open_in[i], open[i]));
    }
    else if (energy[i] != 0.0)
    {
      CHECK(std::lround(tiers[i]) == 2);
      CHECK(shares[i] == Ratio(energy_in[i], energy[i]));
    }
    else if (area[i] > 0.0)
    {
      CHECK(std::lround(tiers[i]) == 3);
      CHECK(shares[i] == Ratio(area_in[i], area[i]));
    }
    else
    {
      CHECK(std::lround(tiers[i]) == 4);
      CHECK(shares[i] == 1.0);
    }
    result[std::lround(patches[i])] = shares[i];
  }
  return result;
}

// The patches (1-based) of the patch table whose origin lies on the island edge x = `edge`
// between its corners (the translational cells of that edge; mesh units = m).
std::set<long> EdgeCells(const fs::path &dir, double edge)
{
  const Table table = test::LoadCsv(dir / "surface-response-patch-energy.csv");
  const auto &patch = test::ColumnByHeader(table, "patch");
  const auto &x = test::ColumnByHeader(table, "origin x (m)");
  const auto &z = test::ColumnByHeader(table, "origin z (m)");
  std::set<long> cells;
  for (std::size_t i = 0; i < patch.data.size(); i++)
  {
    if (std::abs(x.data[i] - edge) < 1.0e-9 && z.data[i] > 0.25 + 1.0e-9 &&
        z.data[i] < 0.75 - 1.0e-9)
    {
      cells.insert(std::lround(patch.data[i]));
    }
  }
  return cells;
}

// The shares themselves on the two boxes of the wall x = kMidWall (raw and self-consistent
// rows of `excitation`): every row is its tier's ratio; the cells of the edge x = 0.25
// count whole in the lower box and not at all in the upper one, those of x = 0.75 the
// reverse (a patch the share file does not list has s_p = 0).
void CheckEdgeShares(const fs::path &window, const fs::path &lower, const fs::path &upper,
                     int excitation)
{
  const auto near = EdgeCells(window, 0.25), far = EdgeCells(window, 0.75);
  REQUIRE(near.size() == 4);
  REQUIRE(far.size() == 4);
  for (const int evaluation : {0, 2})
  {
    INFO("share evaluation " << evaluation);
    const auto low = CheckedShareRows(lower, excitation, evaluation);
    const auto up = CheckedShareRows(upper, excitation, evaluation);
    auto Share = [](const std::map<long, double> &shares, long patch)
    {
      const auto found = shares.find(patch);
      return found != shares.end() ? found->second : 0.0;
    };
    for (const long patch : near)
    {
      INFO("cell " << patch << " of the edge x = 0.25");
      CHECK(Share(low, patch) == 1.0);
      CHECK(Share(up, patch) == 0.0);
    }
    for (const long patch : far)
    {
      INFO("cell " << patch << " of the edge x = 0.75");
      CHECK(Share(low, patch) == 0.0);
      CHECK(Share(up, patch) == 1.0);
    }
  }
}

// The 2D automatic library with its response matrices scaled by 1e-6 (the energy columns
// from `first_energy` on; the form of test-correctedcapacitance.cpp, whose device needs it:
// the unscaled correction exceeds the raw energy of the 2D device), written by the root.
void WriteScaledMatrix(const fs::path &source, const fs::path &target,
                       std::size_t first_energy, double scale)
{
  std::ifstream input(source);
  REQUIRE(input);
  std::ofstream output(target);
  std::string line;
  std::getline(input, line);
  output << line << "\n";
  while (std::getline(input, line))
  {
    if (line.empty())
    {
      continue;
    }
    std::vector<std::string> entries;
    std::size_t begin = 0;
    while (begin <= line.size())
    {
      const std::size_t end = line.find(',', begin);
      entries.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
      if (end == std::string::npos)
      {
        break;
      }
      begin = end + 1;
    }
    for (std::size_t i = 0; i < entries.size(); i++)
    {
      output << (i >= first_energy ? fmt::format("{:.17e}", std::stod(entries[i]) * scale)
                                   : entries[i])
             << (i + 1 < entries.size() ? "," : "\n");
    }
  }
}

fs::path WriteScaled2DLibrary(const fs::path &library_path, const fs::path &dir)
{
  const fs::path scaled_path = dir / "fabrication-process-region-target-2d.json";
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream input(library_path);
    REQUIRE(input);
    json library = json::parse(input);
    REQUIRE(library["Models"].size() == 1);
    auto &model = library["Models"][0];
    for (const auto &[key, first_energy] :
         {std::pair<const char *, std::size_t>{"FabricatedMatrix", 2},
          std::pair<const char *, std::size_t>{"ThinMatrix", 2},
          std::pair<const char *, std::size_t>{"FabricatedSurfaceMatrix", 5},
          std::pair<const char *, std::size_t>{"ThinSurfaceMatrix", 5}})
    {
      const fs::path source = model[key].get<std::string>();
      const fs::path target = dir / ("region-target-2d-" + source.filename().string());
      WriteScaledMatrix(source, target, first_energy, 1.0e-6);
      model[key] = target.string();
    }
    std::ofstream output(scaled_path);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());
  return scaled_path;
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

// The text a callable writes to the process's standard output (the root's log; the form of
// test-conductorconsistency.cpp).
template <typename F>
std::string CaptureStdout(const fs::path &file, F &&run)
{
  std::fflush(stdout);
  const int saved = dup(fileno(stdout));
  std::FILE *target = std::fopen(file.c_str(), "w");
  REQUIRE(target);
  dup2(fileno(target), fileno(stdout));
  run();
  std::fflush(stdout);
  dup2(saved, fileno(stdout));
  close(saved);
  std::fclose(target);
  std::ifstream input(file);
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

// The number of occurrences of `text` in `log`.
int Occurrences(const std::string &log, const std::string &text)
{
  int count = 0;
  for (auto at = log.find(text); at != std::string::npos; at = log.find(text, at + 1))
  {
    count++;
  }
  return count;
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
  // Both crack forms of the island sheet (Legacy: cracked; Default: the key unset).
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  json config = ConvexIslandConfig(crack_policy);
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
    const fs::path output =
        temp.temp_dir / (stem + "-" + test::CrackPolicyName(crack_policy));
    run["Problem"]["Output"] = output.string();
    test::RunElectrostatic(run);
    Mpi::Barrier(Mpi::World());
    return output;
  };
  const fs::path window = Run("region-window", json());
  const fs::path lower = Run("region-lower", HalfSpaceRegion(true));
  const fs::path upper = Run("region-upper", HalfSpaceRegion(false));
  const fs::path mid_lower = Run("region-mid-lower", HalfSpaceRegion(true, kMidWall));
  const fs::path mid_upper = Run("region-mid-upper", HalfSpaceRegion(false, kMidWall));
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
  // The shares themselves: each row its tier's ratio on the boxes that cut patches; the
  // edges whole on their side of the wall x = kMidWall.
  for (const fs::path &dir : {lower, upper})
  {
    INFO("run " << dir.filename().string());
    CheckedShareRows(dir, 1, 0);
    CheckedShareRows(dir, 1, 2);
  }
  CheckEdgeShares(window, mid_lower, mid_upper, 1);
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
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  json config = ConvexMaxwellIslandConfig(crack_policy);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  UseDefaultConstruction(
      config["Solver"]["SurfaceResponseCorrection"],
      WriteDefaultConstructionLibrary(convex_library_3d_path, temp.temp_dir), false);
  config["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both";
  config["Solver"]["SurfaceResponseCorrection"]["PatchEnergy"] = true;
  // The two lowest modes of the 1 m PEC box (as the staged / unstaged case: with N 1 an
  // 8-rank run measured a converged mode the self-consistent step had no field for); the
  // checks read mode 1.
  config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}};
  auto Run = [&](const std::string &stem, const json &region)
  {
    json run = config;
    if (!region.is_null())
    {
      run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Region"] = region;
    }
    const fs::path output =
        temp.temp_dir / (stem + "-" + test::CrackPolicyName(crack_policy));
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
  const fs::path mid_lower =
      Run("region-maxwell-mid-lower", HalfSpaceRegion(true, kMidWall));
  const fs::path mid_upper =
      Run("region-maxwell-mid-upper", HalfSpaceRegion(false, kMidWall));
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
    CHECK(target["CutShare"].get<double>() > 0.0);
    CHECK(target["FabricatedEnergyRegion"].get<double>() /
              target["FabricatedEnergyWindow"].get<double>() ==
          Catch::Approx(box_ft / window_ft).epsilon(1.0e-10));
  }
  // The shares themselves (mode 1 = excitation 0), as in the electrostatic case.
  for (const fs::path &dir : {lower, upper})
  {
    INFO("run " << dir.filename().string());
    CheckedShareRows(dir, 0, 0);
    CheckedShareRows(dir, 0, 2);
  }
  CheckEdgeShares(window, mid_lower, mid_upper, 0);
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

// Ruling 2251 Q4 where a kept class holds part of an applied patch's claim (the situation
// in which the first two tiers differ; the default construction lays kept and applied
// claims end to end, so no driver case forms it). The operator-level island of the cases
// above with the Region x < kBoxWall on entry 4: the potential x[2] - 0.5 (inconsistent
// with the island's conductor state) makes the conductor-consistency gate exclude the 8
// cells of the edges x = 0.25 / 0.75 (their whole claims kept: tier 2 where they have
// points), and one more kept portion covers x in [0.25, 0.3125] of corner 1's claim along z
// = 0.25 (inside the box; the wall cuts that arm). Corner 1's share is then the in-box part
// of the points of its claim no kept class holds, not of all its points.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Region target share of a patch whose claim is partly kept raw counts only "
    "its points no kept class holds (electrostatic)",
    "[regiontarget][surfaceresponseoperator][conductorconsistency][3d][Serial]"
    "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  json config = ConvexIslandConfig(crack_policy);
  UseDefaultConstruction(
      config["Solver"]["Electrostatic"]["ResponseCorrection"],
      WriteDefaultConstructionLibrary(convex_library_3d_path, temp.temp_dir), true);
  config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Region"] = HalfSpaceRegion(true);
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
  LaplaceOperator laplace(iodata, meshes);
  SurfaceResponseOperator response(iodata, laplace);
  REQUIRE(response.HasRegionClaims());
  REQUIRE(response.GetPatchCount() == 20);
  mfem::ParGridFunction potential(&laplace.GetH1Space().Get());
  mfem::FunctionCoefficient coefficient([](const mfem::Vector &x) { return x[2] - 0.5; });
  potential.ProjectCoefficient(coefficient);
  Vector V;
  potential.GetTrueDofs(V);
  std::set<int> excluded;
  for (const auto &record : response.ApplyConductorConsistencyGate(V, 1))
  {
    if (record.excluded)
    {
      excluded.insert(record.patch);
    }
  }
  REQUIRE(excluded.size() == 8);
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> gate_kept;
  for (const auto &portion : response.GetConductorConsistencyPortions())
  {
    gate_kept.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
  }
  REQUIRE(!gate_kept.empty());
  const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> corner_kept = {
      {{0.25, 0.5, 0.25}, {0.3125, 0.5, 0.25}, "kept", -1}};
  const std::vector<SurfacePostOperator::CoveragePartitionClass> kept = {
      {"ConductorConsistency", &gate_kept}, {"Kept", &corner_kept}};
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace);
  mfem::ParDiscreteLinearOperator gradient(&laplace.GetH1Space().Get(),
                                           &laplace.GetNDSpace().Get());
  gradient.AddDomainInterpolator(new mfem::GradientInterpolator());
  gradient.Assemble();
  GridFunction field(laplace.GetNDSpace());
  gradient.Mult(potential, field.Real());
  field.Real() *= -1.0;
  Vector e(laplace.GetNDSpace().GetTrueVSize());
  field.Real().GetTrueDofs(e);
  field.Real().SetFromTrueDofs(e);
  field.Real().ExchangeFaceNbrData();
  const auto footprints = post_op.GetSurfacePostOp().GetInterfaceRegionPatchFootprints(
      {4}, field, nullptr, kept, response.GetRegionClaimPortions(),
      response.GetRegionClaimPatches(), response.GetRegionReferencePoints());
  const auto shares = SurfaceResponseOperator::ComputeRegionPatchShares(footprints);
  const auto &list = footprints.at(4);
  const auto &target = shares.at(4);
  REQUIRE(list.size() == 20);
  // Every patch: the first tier with a denominator, and that tier's ratio.
  for (std::size_t p = 0; p < list.size(); p++)
  {
    INFO("patch " << p + 1);
    const auto &f = list[p];
    const auto Ratio = [](double part, double whole)
    { return Catch::Approx(part / whole).epsilon(1.0e-14).margin(1.0e-15); };
    const int tier = f.open_energy != 0.0 ? 1 : f.energy != 0.0 ? 2 : f.area > 0.0 ? 3 : 4;
    CHECK(target.tier[p] == tier);
    CHECK(target.share[p] == (tier == 1   ? Ratio(f.open_energy_in, f.open_energy)
                              : tier == 2 ? Ratio(f.energy_in, f.energy)
                              : tier == 3 ? Ratio(f.area_in, f.area)
                                          : Ratio(f.reference_in ? 1.0 : 0.0, 1.0)));
    if (excluded.count(static_cast<int>(p)))
    {
      CHECK(f.open_energy == 0.0);  // the gate's class holds its whole claim
    }
  }
  CHECK(std::count(target.tier.begin(), target.tier.end(), 2) >= 4);
  // Corner 1 (origin (0.25, 0.5, 0.25)): cut by the wall, part of its claim kept, and a
  // first-tier share that differs from the energy ratio of all its points (found by its
  // origin among the reference points: every placed patch's, on every rank).
  const auto &origins = response.GetRegionReferencePoints();
  const auto corner = std::find_if(origins.begin(), origins.end(),
                                   [](const auto &origin)
                                   {
                                     return std::abs(origin[0] - 0.25) < 1.0e-12 &&
                                            std::abs(origin[2] - 0.25) < 1.0e-12;
                                   });
  REQUIRE(corner != origins.end());
  const auto c = static_cast<std::size_t>(std::distance(origins.begin(), corner));
  const auto &f = list[c];
  CHECK(target.tier[c] == 1);
  CHECK(target.placement[c] == SurfaceResponseOperator::RegionPlacement::CUT);
  CHECK(f.open_energy > 0.0);
  CHECK(f.open_energy < f.energy);
  CHECK(target.share[c] > 0.0);
  CHECK(target.share[c] < 1.0);
  CHECK(target.share[c] ==
        Catch::Approx(f.open_energy_in / f.open_energy).epsilon(1.0e-14));
  CHECK(std::abs(target.share[c] - f.energy_in / f.energy) > 1.0e-3);
#endif
}

// Decision 2257 MAJOR-2 (the text): on a two-dimensional mesh the warning says that every
// patch counts whole where its reference point lies and that the run is not refused; on a
// three-dimensional mesh it is due only when the wall cuts patches.
TEST_CASE("Region target warning text", "[regiontarget][surfaceresponseoperator][Serial]")
{
  using Catch::Matchers::ContainsSubstring;
  SurfaceResponseOperator::RegionTargetResponse response;
  response.patches_inside = 1;
  response.patches_outside = 1;
  response.patches_by_tier = {0, 0, 0, 2};
  const std::string text =
      SurfaceResponseOperator::RegionTargetWarning(4, response, 1, "Raw", true);
  CHECK_THAT(text, ContainsSubstring("REGION TARGET 4 ON A TWO-DIMENSIONAL MESH"));
  CHECK_THAT(text, ContainsSubstring("every patch counts whole where its reference point "
                                     "(its origin) lies"));
  CHECK_THAT(text, ContainsSubstring("a wall within the matching distance R of an edge "
                                     "splits the raw terms but not the patches"));
  CHECK_THAT(text, ContainsSubstring("The run is not refused."));
  CHECK(SurfaceResponseOperator::RegionTargetWarning(4, response, 1, "Raw", false).empty());
  response.patches_cut = 2;
  response.cut_share = 0.25;
  CHECK_THAT(SurfaceResponseOperator::RegionTargetWarning(4, response, 1, "Raw", false),
             ContainsSubstring("its region wall cuts 2 patch(es) carrying 2.500e-01"));
}

// Decision 2257 MAJOR-2 (the run): the 2D automatic problem (MakeAutomatic2DMesh: the metal
// line y = 0.5, x in [0.25, 0.75]; its SA entry 4 at R = 0.1) with entry 4 in the box x <
// 0.3 (a wall within R of the line end x = 0.25) or x < 0.5 (the corrected-capacitance
// device: the line as terminal 2, the library scaled by 1e-6). No claims exist: every patch
// is of tier 4 and the one at x = 0.25 counts whole in both boxes, while the raw energy of
// entry 4 moves with the wall; the warning prints once per construction; the runs end.
TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Region target on a two-dimensional mesh counts every patch whole where its "
    "reference point lies, with a warning",
    "[regiontarget][electrostaticsolver][surfaceresponseoperator][2d][Serial]"
    "[Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  using Catch::Matchers::ContainsSubstring;
  const auto crack_policy = GENERATE(test::CrackPolicy::Legacy, test::CrackPolicy::Default);
  INFO("crack policy: " << test::CrackPolicyName(crack_policy));
  const fs::path mesh_path = temp.temp_dir / "region-target-2d.mesh";
  {
    mfem::Mesh serial = MakeAutomatic2DMesh();
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream output(mesh_path);
      serial.Print(output);
    }
  }
  Mpi::Barrier(Mpi::World());
  json config = AutomaticConfig2D(crack_policy);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      WriteScaled2DLibrary(library_path, temp.temp_dir).string();
  config["Boundaries"]["Ground"] = {{"Attributes", {1, 3, 4}}};
  config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {2}}},
                                      {{"Index", 2}, {"Attributes", {9, 10}}}};
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 100}};
  auto Run = [&](const std::string &stem, double wall)
  {
    json run = config;
    run["Boundaries"]["Postprocessing"]["Dielectric"][0]["Region"] =
        HalfSpaceRegion(true, wall);
    const std::string name = stem + "-" + test::CrackPolicyName(crack_policy);
    const fs::path output = temp.temp_dir / name;
    run["Problem"]["Output"] = output.string();
    const std::string log = CaptureStdout(
        temp.temp_dir / fmt::format("{}-{}.log", name, Mpi::Rank(Mpi::World())),
        [&]() { test::RunElectrostatic(run); });
    Mpi::Barrier(Mpi::World());
    return std::make_pair(output, log);
  };
  const auto near = Run("region-2d-near", 0.3);
  const auto far = Run("region-2d-far", 0.5);
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  auto Fabricated = [](const fs::path &dir)
  {
    const json record = RegionRecord(dir);
    REQUIRE(!record.is_null());
    REQUIRE(record.at("Evaluations").size() >= 1);
    for (const auto &evaluation : record["Evaluations"])
    {
      const auto &target = evaluation["Targets"][0];
      CHECK(target["Index"].get<int>() == 4);
      CHECK(target["PatchesCut"].get<int>() == 0);
      const auto tiers = target["PatchesByTier"].get<std::array<int, 4>>();
      CHECK(tiers[0] + tiers[1] + tiers[2] == 0);
      CHECK(tiers[3] > 0);
    }
    return record["Evaluations"][0]["Targets"][0]["FabricatedEnergyRegion"].get<double>();
  };
  for (const auto &[dir, log] : {near, far})
  {
    INFO("run " << dir.filename().string());
    CHECK(Occurrences(log, "ON A TWO-DIMENSIONAL MESH") == 1);
    CHECK_THAT(log,
               ContainsSubstring("every patch counts whole where its reference point"));
  }
  const double near_fabricated = Fabricated(near.first);
  CHECK(near_fabricated != 0.0);
  CHECK(near_fabricated == Catch::Approx(Fabricated(far.first)).epsilon(1.0e-12));
  auto Raw = [](const fs::path &dir)
  { return Single(test::LoadCsv(dir / "surface-Q.csv"), "p_surf[4]", {{"i", 1}}); };
  CHECK(Raw(near.first) != Catch::Approx(Raw(far.first)).epsilon(1.0e-6));
#endif
}

}  // namespace palace
