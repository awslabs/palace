// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "drivers/electrostaticsolver.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/jsonschema.hpp"

using json = nlohmann::json;
using namespace palace;
using Catch::Matchers::ContainsSubstring;

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

  // The resolved configuration: unchanged without the block; the block as given with it
  // (an absent key stays absent: its alias variable decides).
  {
    const auto config = ElectrostaticConfig({{"Save", 0}});
    IoData iodata(config, false);
    const auto resolved = IoData::ConcretizeDefaults(iodata, config);
    CHECK(!resolved["Solver"]["Electrostatic"].contains("ResponseArchive"));
  }
  for (const json &given : {json::object(), json({{"Directory", "a"}}),
                            json({{"Mode", "Reduce"}, {"Directory", "a"}})})
  {
    const auto config = ElectrostaticConfig({{"ResponseArchive", given}});
    IoData iodata(config, false);
    const auto resolved = IoData::ConcretizeDefaults(iodata, config);
    CHECK(resolved["Solver"]["Electrostatic"]["ResponseArchive"] == given);
  }
}

// The alias rule of Solver.Electrostatic.ResponseArchive: key only -> the key, variable
// only -> the variable, both equal -> that value, both different -> an error; the
// variables are the root's on every rank.
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
    CHECK(!settings.Active());
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

  // The archive worker by its variables, by its keys, and by both (equal).
  const auto worker_keys =
      Keys({{"Mode", "Archive"}, {"Directory", "archive/./w"}, {"SourceTiming", true}});
  const auto CheckWorker = [&](const ResponseArchiveSettings &settings)
  {
    CHECK(settings.Active());
    CHECK(settings.archive_only);
    CHECK(!settings.reduce_only);
    CHECK(!settings.estimate_only);
    REQUIRE(settings.directory);
    CHECK(*settings.directory == Path("archive/w"));
    CHECK(settings.source_timing);
    CHECK(!settings.recycle_initial_guess);
  };
  {
    const auto settings = ResolveResponseArchive(worker_keys, comm);
    CheckWorker(settings);
    CHECK(settings.mode_source == ResponseArchiveSource::KEY);
    CHECK(settings.directory_source == ResponseArchiveSource::KEY);
    CHECK(settings.block_size_source == ResponseArchiveSource::DEFAULT);
  }
  Set("PALACE_RESPONSE_ARCHIVE_ONLY", "1");
  Set("PALACE_RESPONSE_SOURCE_TIMING", "1");
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "archive/w");
  {
    const auto settings = ResolveResponseArchive(config::ElectrostaticSolverData(), comm);
    CheckWorker(settings);
    CHECK(settings.mode_source == ResponseArchiveSource::VARIABLE);
    CHECK(settings.source_timing_source == ResponseArchiveSource::VARIABLE);
  }
  {
    const auto settings = ResolveResponseArchive(worker_keys, comm);
    CheckWorker(settings);
    CHECK(settings.mode_source == ResponseArchiveSource::KEY_AND_VARIABLE);
    CHECK(settings.directory_source == ResponseArchiveSource::KEY_AND_VARIABLE);
    CHECK(settings.source_timing_source == ResponseArchiveSource::KEY_AND_VARIABLE);
  }
  // A partial block: the absent keys are decided by their variables, and the log line
  // names each setting's source.
  {
    const auto settings = ResolveResponseArchive(Keys({{"Directory", "archive/w"}}), comm);
    CheckWorker(settings);
    CHECK(ResponseArchiveSummary(settings) ==
          "Response-archive settings: Mode=Archive (variable), Directory=" +
              Path("archive/w").string() +
              " (key and variable), BlockSize=3 (default), SourceTiming=true (variable), "
              "RecycleInitialGuess=false (default), EstimateRequest=- (default)");
  }
  // A relative variable and an absolute key naming the same directory agree (compared on
  // the root).
  CheckWorker(ResolveResponseArchive(
      Keys({{"Mode", "Archive"}, {"Directory", Path("archive/w").string()}}), comm));
  // Both set and different: refused for every key.
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Mode", "Reduce"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Mode", "None"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"Directory", "archive/v"}}), comm));
  CHECK_THROWS(ResolveResponseArchive(Keys({{"SourceTiming", false}}), comm));
  Set("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "0");
  CHECK_THROWS(ResolveResponseArchive(Keys({{"RecycleInitialGuess", true}}), comm));
  CHECK(!ResolveResponseArchive(Keys({{"RecycleInitialGuess", false}}), comm)
             .recycle_initial_guess);
  // Recycling is turned on by its variable and by its key.
  Set("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "1");
  CHECK(ResolveResponseArchive(config::ElectrostaticSolverData(), comm)
            .recycle_initial_guess);
  ResponseArchiveEnvironmentGuard::Clear();
  CHECK(ResolveResponseArchive(Keys({{"RecycleInitialGuess", true}}), comm)
            .recycle_initial_guess);
  // The Reduce and Estimate variables against a different Mode key.
  Set("PALACE_RESPONSE_REDUCE_ONLY", "1");
  CHECK_THROWS_WITH(ResolveResponseArchive(Keys({{"Mode", "Archive"}}), comm),
                    ContainsSubstring("PALACE_RESPONSE_REDUCE_ONLY"));
  ResponseArchiveEnvironmentGuard::Clear();
  Set("PALACE_RESPONSE_ESTIMATE_ONLY", "1");
  CHECK_THROWS_WITH(ResolveResponseArchive(Keys({{"Mode", "None"}}), comm),
                    ContainsSubstring("PALACE_RESPONSE_ESTIMATE_ONLY"));
  ResponseArchiveEnvironmentGuard::Clear();

  // The archive reducer; BLOCK_SIZE is parsed by a reduction only.
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

  // The estimate: the refused settings, set at all, in a fixed order.
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

// The startup checks of every run (ValidateArchiveEstimateOptions): on the key path, a key
// with no effect in the resolved mode and a key-given Archive / Reduce mode without a
// Directory are refused naming the key; the variables keep their silent semantics; an
// archive estimate refuses its conflicting settings.
TEST_CASE("Response archive startup checks", "[responsearchive][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  ResponseArchiveEnvironmentGuard guard;
  const auto Set = [](const char *name, const char *value) { setenv(name, value, 1); };
  const auto Validate = [&](const json &electrostatic)
  {
    const IoData iodata(ElectrostaticConfig(electrostatic), false);
    ValidateArchiveEstimateOptions(iodata, comm, false);
  };
  const auto Block = [](const json &archive) { return json{{"ResponseArchive", archive}}; };

  // The key path.
  CHECK_NOTHROW(Validate(json::object()));
  CHECK_NOTHROW(Validate(Block(json::object())));
  CHECK_NOTHROW(
      Validate(Block({{"Mode", "Archive"}, {"Directory", "a"}, {"SourceTiming", true}})));
  CHECK_NOTHROW(
      Validate(Block({{"Mode", "Reduce"}, {"Directory", "a"}, {"BlockSize", 8}})));
  CHECK_NOTHROW(Validate(Block({{"Mode", "None"}, {"Directory", "a"}})));
  CHECK_THROWS_WITH(Validate(Block({{"SourceTiming", true}})),
                    ContainsSubstring("\"SourceTiming\" has no effect"));
  CHECK_THROWS_WITH(
      Validate(
          Block({{"Mode", "None"}, {"Directory", "a"}, {"RecycleInitialGuess", false}})),
      ContainsSubstring("\"RecycleInitialGuess\" has no effect"));
  CHECK_THROWS_WITH(
      Validate(Block({{"Mode", "Archive"}, {"Directory", "a"}, {"BlockSize", 8}})),
      ContainsSubstring("\"BlockSize\" has no effect"));
  CHECK_THROWS_WITH(
      Validate(
          Block({{"Mode", "Reduce"}, {"Directory", "a"}, {"EstimateRequest", "r.json"}})),
      ContainsSubstring("\"EstimateRequest\" has no effect"));
  CHECK_THROWS_WITH(Validate(Block({{"Mode", "Archive"}})),
                    ContainsSubstring("requires a \"Directory\""));
  CHECK_THROWS_WITH(Validate(Block({{"Mode", "Reduce"}, {"BlockSize", 8}})),
                    ContainsSubstring("requires a \"Directory\""));

  // A key made effective by a variable's mode, a key-given mode with the variable's
  // Directory.
  Set("PALACE_RESPONSE_ARCHIVE_ONLY", "1");
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "a");
  CHECK_NOTHROW(Validate(Block({{"SourceTiming", true}})));
  ResponseArchiveEnvironmentGuard::Clear();
  Set("PALACE_RESPONSE_ARCHIVE_DIR", "a");
  CHECK_NOTHROW(Validate(Block({{"Mode", "Archive"}})));
  ResponseArchiveEnvironmentGuard::Clear();

  // The variables path: no startup refusal of a variable without effect, nor of a Reduce
  // mode without a Directory (the solve refuses it, as before the keys existed).
  Set("PALACE_RESPONSE_SOURCE_TIMING", "1");
  Set("PALACE_RESPONSE_BLOCK_SIZE", "4");
  Set("PALACE_RESPONSE_ESTIMATE_REQUEST", "r.json");
  CHECK_NOTHROW(Validate(json::object()));
  ResponseArchiveEnvironmentGuard::Clear();
  Set("PALACE_RESPONSE_REDUCE_ONLY", "1");
  CHECK_NOTHROW(Validate(json::object()));
  ResponseArchiveEnvironmentGuard::Clear();

  // The archive estimate refuses a conflicting variable set at all, and its BlockSize key.
  Set("PALACE_RESPONSE_ESTIMATE_ONLY", "1");
  Set("PALACE_RESPONSE_ARCHIVE_ONLY", "0");
  CHECK_THROWS_WITH(Validate(json::object()),
                    ContainsSubstring("conflicts with PALACE_RESPONSE_ARCHIVE_ONLY"));
  ResponseArchiveEnvironmentGuard::Clear();
  CHECK_THROWS_WITH(
      Validate(
          Block({{"Mode", "Estimate"}, {"EstimateRequest", "r.json"}, {"BlockSize", 4}})),
      ContainsSubstring("conflicts with Solver.Electrostatic.ResponseArchive.BlockSize"));
}
