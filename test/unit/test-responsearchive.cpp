// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "drivers/electrostaticsolver.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;

namespace
{

using ResponseArchiveMode = config::ElectrostaticSolverData::ResponseArchiveData::Mode;

json ElectrostaticConfig(const json &electrostatic)
{
  return {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
          {"Model", {{"Mesh", "test.msh"}}},
          {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
          {"Boundaries", json::object()},
          {"Solver", {{"Electrostatic", electrostatic}}}};
}

const char *const response_archive_variables[] = {"PALACE_RESPONSE_ARCHIVE_ONLY",
                                                  "PALACE_RESPONSE_REDUCE_ONLY",
                                                  "PALACE_RESPONSE_ESTIMATE_ONLY",
                                                  "PALACE_RESPONSE_ARCHIVE_DIR",
                                                  "PALACE_RESPONSE_BLOCK_SIZE",
                                                  "PALACE_RESPONSE_SOURCE_TIMING",
                                                  "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS",
                                                  "PALACE_RESPONSE_ESTIMATE_REQUEST"};

// Unsets the alias variables on this rank at construction and destruction.
struct ResponseArchiveEnvironmentGuard
{
  ResponseArchiveEnvironmentGuard() { Clear(); }
  ~ResponseArchiveEnvironmentGuard() { Clear(); }
  static void Clear()
  {
    for (const char *name : response_archive_variables)
    {
      unsetenv(name);
    }
  }
};

config::ElectrostaticSolverData Keys(const json &archive)
{
  return config::ElectrostaticSolverData(json{{"ResponseArchive", archive}});
}

}  // namespace

TEST_CASE("Response archive keys: schema, parse and resolved configuration",
          "[responsearchive][config][Serial]")
{
  const json full = {{"Mode", "Archive"},
                     {"Directory", "archive"},
                     {"BlockSize", 48},
                     {"SourceTiming", true},
                     {"RecycleInitialGuess", false},
                     {"EstimateRequest", "request.json"}};
  CHECK(ValidateConfig(ElectrostaticConfig({{"ResponseArchive", full}})).empty());
  CHECK(ValidateConfig(ElectrostaticConfig({{"ResponseArchive", json::object()}})).empty());
  const std::vector<std::pair<std::string, json>> invalid = {
      {"Mode", "Stream"},      {"BlockSize", 0},      {"Directory", ""},
      {"EstimateRequest", ""}, {"SourceTiming", "1"}, {"Unknown", true}};
  for (const auto &[key, value] : invalid)
  {
    CHECK(!ValidateConfig(ElectrostaticConfig({{"ResponseArchive", {{key, value}}}}))
               .empty());
  }

  // Absent block: no keys (the alias variables decide alone).
  CHECK(!config::ElectrostaticSolverData(json::object()).response_archive);
  const auto empty = Keys(json::object()).response_archive;
  REQUIRE(empty);
  CHECK(!empty->mode);
  CHECK(!empty->directory);
  CHECK(!empty->block_size);
  CHECK(!empty->source_timing);
  CHECK(!empty->recycle_initial_guess);
  CHECK(!empty->estimate_request);
  const auto parsed = Keys(full).response_archive;
  REQUIRE(parsed);
  CHECK(parsed->mode == ResponseArchiveMode::ARCHIVE);
  CHECK(*parsed->directory == "archive");
  CHECK(*parsed->block_size == 48);
  CHECK(*parsed->source_timing);
  CHECK(!*parsed->recycle_initial_guess);
  CHECK(*parsed->estimate_request == "request.json");
  CHECK(Keys({{"Mode", "None"}}).response_archive->mode == ResponseArchiveMode::NONE);
  CHECK(Keys({{"Mode", "Reduce"}}).response_archive->mode == ResponseArchiveMode::REDUCE);
  CHECK(Keys({{"Mode", "Estimate"}}).response_archive->mode ==
        ResponseArchiveMode::ESTIMATE);
  CHECK_THROWS(Keys({{"Mode", "Stream"}}));
  CHECK_THROWS(Keys({{"BlockSize", 0}}));
  CHECK_THROWS(Keys({{"Directory", ""}}));
  CHECK_THROWS(Keys({{"EstimateRequest", ""}}));
  CHECK_THROWS(Keys(json::array()));

  // The resolved configuration: unchanged without the block; its defaults with it.
  {
    const auto config = ElectrostaticConfig({{"Save", 0}});
    IoData iodata(config, false);
    const auto resolved = IoData::ConcretizeDefaults(iodata, config);
    CHECK(!resolved["Solver"]["Electrostatic"].contains("ResponseArchive"));
  }
  {
    const auto config = ElectrostaticConfig(
        {{"ResponseArchive", {{"Mode", "Reduce"}, {"Directory", "a"}}}});
    IoData iodata(config, false);
    const auto resolved = IoData::ConcretizeDefaults(iodata, config);
    CHECK(resolved["Solver"]["Electrostatic"]["ResponseArchive"] ==
          json({{"Mode", "Reduce"},
                {"Directory", "a"},
                {"BlockSize", 3},
                {"SourceTiming", false},
                {"RecycleInitialGuess", false}}));
  }
}

// The alias rule of Solver.Electrostatic.ResponseArchive: key only -> the key, variable
// only
// -> the variable (the job scripts of record keep their meaning), both equal -> that value,
// both different -> an error; the variables are the root's on every rank.
TEST_CASE("Response archive keys and alias variables resolve collectively",
          "[responsearchive][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const bool root = Mpi::Root(comm);
  ResponseArchiveEnvironmentGuard guard;
  const auto Path = [](const std::string &value)
  { return std::filesystem::absolute(value).lexically_normal(); };
  const auto Set = [](const char *name, const char *value) { setenv(name, value, 1); };

  // Nothing given: the ordinary solve.
  {
    const auto settings = ResolveResponseArchive(config::ElectrostaticSolverData(), comm);
    CHECK(!settings.archive_only);
    CHECK(!settings.reduce_only);
    CHECK(!settings.estimate_only);
    CHECK(!settings.directory);
    CHECK(settings.block_size == 3);
    CHECK(!settings.source_timing);
    CHECK(!settings.recycle_initial_guess);
    CHECK(!settings.estimate_request);
    CHECK(settings.estimate_conflicts.empty());
  }

  // The worker of record by its variables, by its keys, and by both (equal).
  const auto worker_keys =
      Keys({{"Mode", "Archive"}, {"Directory", "archive/./w"}, {"SourceTiming", true}});
  const auto CheckWorker = [&](const ResponseArchiveSettings &settings)
  {
    CHECK(settings.archive_only);
    CHECK(!settings.reduce_only);
    CHECK(!settings.estimate_only);
    REQUIRE(settings.directory);
    CHECK(*settings.directory == Path("archive/w"));
    CHECK(settings.source_timing);
    CHECK(!settings.recycle_initial_guess);
  };
  CheckWorker(ResolveResponseArchive(worker_keys, comm));
  Set("PALACE_RESPONSE_ARCHIVE_ONLY", "1");
  Set("PALACE_RESPONSE_SOURCE_TIMING", "1");
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "archive/w");
  CheckWorker(ResolveResponseArchive(config::ElectrostaticSolverData(), comm));
  CheckWorker(ResolveResponseArchive(worker_keys, comm));
  // Both set and different: refused for every key.
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Mode", "Reduce"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Mode", "None"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Directory", "archive/v"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"SourceTiming", false}}), comm));
  Set("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "0");
  CHECK_THROWS(ResolveResponseArchive(Keys({{"RecycleInitialGuess", true}}), comm));
  CHECK(!ResolveResponseArchive(Keys({{"RecycleInitialGuess", false}}), comm)
             .recycle_initial_guess);
  ResponseArchiveEnvironmentGuard::Clear();

  // The reducer of record; BLOCK_SIZE is parsed by a reduction only.
  Set("PALACE_RESPONSE_BLOCK_SIZE", "not a number");
  CHECK(ResolveResponseArchive(config::ElectrostaticSolverData(), comm).block_size == 3);
  Set("PALACE_RESPONSE_REDUCE_ONLY", "1");
  CHECK_THROWS(ResolveResponseArchive(config::ElectrostaticSolverData(), comm));
  Set("PALACE_RESPONSE_BLOCK_SIZE", "0");
  CHECK_THROWS(ResolveResponseArchive(config::ElectrostaticSolverData(), comm));
  Set("PALACE_RESPONSE_BLOCK_SIZE", "48");
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "archive/w");
  const auto reducer_keys =
      Keys({{"Mode", "Reduce"}, {"Directory", "archive/w"}, {"BlockSize", 48}});
  for (const auto &keys : {config::ElectrostaticSolverData(), reducer_keys})
  {
    const auto settings = ResolveResponseArchive(keys, comm);
    CHECK(settings.reduce_only);
    CHECK(!settings.archive_only);
    CHECK(settings.block_size == 48);
    CHECK(*settings.directory == Path("archive/w"));
  }
  CHECK_THROWS(ResolveResponseArchive(Keys({{"BlockSize", 8}}), comm));
  ResponseArchiveEnvironmentGuard::Clear();
  {
    const auto settings = ResolveResponseArchive(reducer_keys, comm);
    CHECK(settings.reduce_only);
    CHECK(settings.block_size == 48);
    CHECK(*settings.directory == Path("archive/w"));
  }

  // An empty ARCHIVE_DIR is unset; a key then decides alone.
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "");
  CHECK(!ResolveResponseArchive(config::ElectrostaticSolverData(), comm).directory);
  CHECK(*ResolveResponseArchive(Keys({{"Directory", "d"}}), comm).directory == Path("d"));
  ResponseArchiveEnvironmentGuard::Clear();

  // The estimate: the refused settings, set at all, in the order of record.
  Set("PALACE_RESPONSE_ESTIMATE_ONLY", "1");
  Set("PALACE_RESPONSE_ESTIMATE_REQUEST", "request.json");
  Set("PALACE_RESPONSE_BLOCK_SIZE", "4");
  Set("PALACE_RESPONSE_ARCHIVE_ONLY", "0");
  Set("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "0");
  {
    const auto settings = ResolveResponseArchive(config::ElectrostaticSolverData(), comm);
    CHECK(settings.estimate_only);
    CHECK(!settings.archive_only);
    CHECK(*settings.estimate_request == "request.json");
    CHECK(settings.estimate_conflicts ==
          std::vector<std::string>{"PALACE_RESPONSE_ARCHIVE_ONLY",
                                   "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS",
                                   "PALACE_RESPONSE_BLOCK_SIZE"});
  }
  CHECK_THROWS(ResolveResponseArchive(Keys({{"EstimateRequest", "other.json"}}), comm));
  ResponseArchiveEnvironmentGuard::Clear();
  {
    const auto settings = ResolveResponseArchive(
        Keys({{"Mode", "Estimate"}, {"EstimateRequest", "request.json"}, {"BlockSize", 4}}),
        comm);
    CHECK(settings.estimate_only);
    CHECK(*settings.estimate_request == "request.json");
    CHECK(settings.estimate_conflicts ==
          std::vector<std::string>{"Solver.Electrostatic.ResponseArchive.BlockSize"});
  }

  // The root's variables on every rank: set on the root only, then on the others only.
  if (root)
  {
    Set("PALACE_RESPONSE_ARCHIVE_ONLY", "1");
    Set("PALACE_RESPONSE_ARCHIVE_DIR", "archive/w");
  }
  {
    const auto settings = ResolveResponseArchive(config::ElectrostaticSolverData(), comm);
    CHECK(settings.archive_only);
    REQUIRE(settings.directory);
    CHECK(*settings.directory == Path("archive/w"));
    CHECK_NOTHROW(ResolveResponseArchive(worker_keys, comm));
  }
  ResponseArchiveEnvironmentGuard::Clear();
  if (!root)
  {
    Set("PALACE_RESPONSE_REDUCE_ONLY", "1");
  }
  CHECK(!ResolveResponseArchive(config::ElectrostaticSolverData(), comm).reduce_only);
  CHECK_NOTHROW(ResolveResponseArchive(worker_keys, comm));
}
