// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <streambuf>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "utils/metadatafile.hpp"

using namespace palace;
using json = nlohmann::json;

namespace
{

// The fields every driver writes (BaseSolver: the constructor, the spaces, the solver, the
// timer and the memory records).
json CommonMetadata()
{
  json meta;
  meta["GitTag"] = "v0.13.0-memfix";
  meta["Problem"]["MPISize"] = 48;
  meta["Problem"]["OpenMPThreads"] = 1;
  meta["Problem"]["MeshElements"] = 1234567;
  meta["Problem"]["DegreesOfFreedom"] = 9876543210LL;
  meta["Problem"]["MultigridDegreesOfFreedom"] = {1024, 8192, 9876543210LL};
  meta["LinearSolver"]["TotalSolves"] = 12;
  meta["LinearSolver"]["TotalIts"] = 345;
  for (const char *key : {"Initialization", "Solve", "Postprocessing", "IO", "Total"})
  {
    meta["ElapsedTime"]["Durations"][key] = 1.0 / 3.0;
    meta["ElapsedTime"]["MinimumDurations"][key] = 1.0e-300;
    meta["ElapsedTime"]["MaximumDurations"][key] = 6.02214076e23;
    meta["ElapsedTime"]["AverageDurations"][key] = -0.0;
    meta["ElapsedTime"]["Counts"][key] = 7;
    meta["PeakMemoryGrowthMegabytes"]["Min"][key] = 0.1;
  }
  meta["PeakMemoryMegabytes"] = {
      {"Min", 1.5}, {"Max", 2048.25}, {"Average", 100.125}, {"Total", 4.0e6}};
  meta["PeakNodeMemoryMegabytes"] = meta["PeakMemoryMegabytes"];
  return meta;
}

// A gate record entry in the form of FinishConductorConsistencyGate (probe snaps, the
// self-consistent pre-gate flag, a Maxwell frequency).
json GateRecord(int patch, int source, bool maxwell, bool pre_gate, int snaps)
{
  json entry = {{"Patch", patch},
                {"Model", "isolated-edge \"quoted\" / unicode \xc2\xb5m"},
                {"ModelIndex", patch % 5},
                {"Source", source},
                {"Amplitude", 1.0e-17 * patch},
                {"State", 0.5},
                {"FloorApplied", patch % 2 == 0},
                {"MaxRatio", 3.0e-9},
                {"ProbeKind", "TranslationalEdge"},
                {"Maxwell", maxwell},
                {"Frequency", maxwell ? 5.123456789e9 : 0.0},
                {"SelfConsistentPreGate", pre_gate},
                {"RecordOnly", false},
                {"WorstPoint", {-1.25e-3, 2.5e-3, 0.0}},
                {"Excluded", false}};
  json knots = json::array();
  for (int i = 0; i < snaps; i++)
  {
    knots.push_back({{"Knot", i}, {"Segment", 2 * i}, {"Arc", -1}, {"Displacement", 1e-9}});
  }
  entry["ProbeSnaps"] = {
      {"Count", snaps}, {"MaxDisplacement", 1e-9 * snaps}, {"Knots", std::move(knots)}};
  return entry;
}

json SurfaceResponseMetadata(bool maxwell, int sources)
{
  json response = {{"Version", 1}};
  response["ModelCatalog"] = json::array({{{"Index", 0}, {"Name", "isolated-edge"}}});
  auto &gate = response["Diagnostics"]["ConductorConsistency"];
  gate["Count"] = 0;
  gate["ExcludedPatches"] = json::array();
  gate["Records"] = json::array();
  for (int source = 1; source <= sources; source++)
  {
    for (int patch = 0; patch < 3; patch++)
    {
      gate["Records"].push_back(GateRecord(patch, source, maxwell, source == 1, patch));
    }
  }
  response["Diagnostics"]["EdgeExclusion"] = json::object();
  return response;
}

// An electrostatic, an eigenmode and a driven-form metadata object.
json ElectrostaticMetadata()
{
  json meta = CommonMetadata();
  meta["SurfaceResponse"] = SurfaceResponseMetadata(false, 4);
  meta["SurfaceResponse"]["CorrectedSolver"]["Electrostatic"] = {{"LinearSolves", 4},
                                                                 {"LinearIterations", 99}};
  return meta;
}

json EigenmodeMetadata()
{
  json meta = CommonMetadata();
  meta["SurfaceResponse"] = SurfaceResponseMetadata(true, 3);
  meta["AdaptiveMeshRefinement"]["MaxSizeStop"] = nullptr;
  return meta;
}

json DrivenMetadata()
{
  json meta = CommonMetadata();
  meta["Excitations"] = {{"1", {{"LumpedPort", {1, 2}}, {"WavePort", json::array()}}},
                         {"2", {{"LumpedPort", {3}}}}};
  meta["SurfaceResponse"] = SurfaceResponseMetadata(true, 2);
  meta["Empty"] = {{"Object", json::object()}, {"Array", json::array()}};
  return meta;
}

// A stream buffer that accepts `capacity` bytes and then refuses (a full disk).
class BoundedBuffer : public std::streambuf
{
public:
  explicit BoundedBuffer(std::size_t capacity) : capacity(capacity) {}
  std::string data;

protected:
  int_type overflow(int_type c) override
  {
    if (traits_type::eq_int_type(c, traits_type::eof()))
    {
      return traits_type::not_eof(c);
    }
    if (data.size() >= capacity)
    {
      return traits_type::eof();
    }
    data.push_back(traits_type::to_char_type(c));
    return c;
  }
  std::streamsize xsputn(const char *s, std::streamsize n) override
  {
    const auto room = static_cast<std::streamsize>(capacity - data.size());
    const auto taken = n < room ? n : room;
    data.append(s, static_cast<std::size_t>(taken));
    return taken;
  }

private:
  std::size_t capacity;
};

}  // namespace

// Lane fix-response-memory round 2a R1 (decision 1475 (1)): palace.json is streamed, not
// built as a dump() string first; its bytes must be dump(2) + '\n' exactly (the file is a
// record compared byte for byte), and a short write must be detected (fail loud).
TEST_CASE("Streamed metadata file is byte-identical to dump(2)",
          "[metadatafile][memfix][Serial]")
{
  const auto dir = std::filesystem::temp_directory_path() / "palace-metadatafile-test";
  std::filesystem::create_directories(dir);
  for (const auto &[form, meta] :
       {std::pair<const char *, json>{"electrostatic", ElectrostaticMetadata()},
        {"eigenmode", EigenmodeMetadata()},
        {"driven", DrivenMetadata()}})
  {
    INFO(form);
    const std::string expected = meta.dump(2) + "\n";
    std::ostringstream os;
    CHECK(StreamJsonMetadata(os, meta));
    CHECK(os.str() == expected);

    const auto path = dir / (std::string(form) + ".json");
    {
      std::ofstream fo(path);
      REQUIRE(fo.is_open());
      CHECK(StreamJsonMetadata(fo, meta));
      fo.close();
      CHECK(!fo.fail());
    }
    std::ifstream fi(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(fi)),
                            std::istreambuf_iterator<char>());
    CHECK(bytes == expected);
    CHECK(json::parse(bytes) == meta);

    // A stream that takes only part of the text (a full disk) reports the short write.
    for (std::size_t capacity :
         {std::size_t(0), std::size_t(1), expected.size() / 2, expected.size() - 1})
    {
      BoundedBuffer buffer(capacity);
      std::ostream bounded(&buffer);
      CHECK_FALSE(StreamJsonMetadata(bounded, meta));
      CHECK(buffer.data == expected.substr(0, capacity));
    }
    BoundedBuffer exact(expected.size());
    std::ostream fits(&exact);
    CHECK(StreamJsonMetadata(fits, meta));
    CHECK(exact.data == expected);
  }
  std::filesystem::remove_all(dir);
}

namespace
{

// The form of a spilled palace.json: the Records entries moved from `meta` to a spill (one
// Flush per source / mode, as the gate passes do); `meta` keeps the empty placeholder.
json SpillRecords(json meta, JsonArraySpill &spill)
{
  auto &records = meta["SurfaceResponse"]["Diagnostics"]["ConductorConsistency"]["Records"];
  int source = records.empty() ? 0 : records.front()["Source"].get<int>();
  for (const auto &entry : records)
  {
    if (entry["Source"].get<int>() != source)
    {
      REQUIRE(spill.Flush());
      source = entry["Source"].get<int>();
    }
    REQUIRE(spill.Append(entry));
  }
  REQUIRE(spill.Flush());
  records = json::array();
  return meta;
}

std::string Streamed(const json &meta, const JsonArraySpill *spill)
{
  std::ostringstream os;
  REQUIRE(StreamJsonMetadata(os, meta, spill));
  return os.str();
}

}  // namespace

// Lane fix-response-memory round 2a R2 (decision 1475 (1)): the root's
// conductor-consistency Records are spilled to a scratch file in their palace.json text and
// spliced at the metadata write; the file must be byte-identical to the in-memory form:
// multi-source ES, multi-mode eigen, driven, probe snaps, the self-consistent pre-gate
// records, an EMPTY record set, no record place at all, and successive AMR iterations in
// one directory (each written, re-read and re-written as the later SaveMetadata calls do).
TEST_CASE("Spilled metadata Records are spliced byte-identically",
          "[metadatafile][memfix][Serial]")
{
  const auto dir = std::filesystem::temp_directory_path() / "palace-metadatafile-spill";
  std::filesystem::create_directories(dir);
  const std::string path = (dir / ".records.spill").string();
  const std::vector<std::string> keys = {"SurfaceResponse", "Diagnostics",
                                         "ConductorConsistency", "Records"};
  json empty = ElectrostaticMetadata();
  empty["SurfaceResponse"]["Diagnostics"]["ConductorConsistency"]["Records"] =
      json::array();
  json no_place = CommonMetadata();
  for (const auto &[form, meta] :
       {std::pair<const char *, json>{"electrostatic", ElectrostaticMetadata()},
        {"eigenmode", EigenmodeMetadata()},
        {"driven", DrivenMetadata()},
        {"empty", empty},
        {"no-place", no_place}})
  {
    INFO(form);
    const std::string expected = meta.dump(2) + "\n";
    JsonArraySpill spill(path, keys);
    CHECK(spill.EntryIndent() == 10);
    const json placeholder =
        meta.contains("SurfaceResponse") ? SpillRecords(meta, spill) : meta;
    CHECK(std::filesystem::exists(path) == (spill.Size() > 0));
    CHECK(Streamed(placeholder, &spill) == expected);
    // Without the spill the placeholder writes the empty array (the records are the
    // file's).
    if (spill.Size() > 0)
    {
      CHECK(Streamed(placeholder, nullptr) != expected);
    }
    CHECK(spill.Remove());
    CHECK(!std::filesystem::exists(path));
  }

  // AMR: two iterations, each operator's spill at the same path; each palace.json is then
  // re-read and re-written (the timer / memory SaveMetadata calls) and copied aside.
  for (int iteration = 0; iteration < 2; iteration++)
  {
    INFO(iteration);
    json meta = iteration == 0 ? ElectrostaticMetadata() : EigenmodeMetadata();
    meta["AdaptiveMeshRefinement"]["Iteration"] = iteration;
    JsonArraySpill spill(path, keys);
    const json placeholder = SpillRecords(meta, spill);
    const auto file = dir / "palace.json";
    {
      std::ofstream fo(file);
      REQUIRE(StreamJsonMetadata(fo, placeholder, &spill));
    }
    CHECK(spill.Remove());
    std::ifstream fi(file);
    json reread = json::parse(fi);
    CHECK(reread == meta);
    reread["PeakMemoryMegabytes"]["Max"] = 1.0 + iteration;
    meta["PeakMemoryMegabytes"]["Max"] = 1.0 + iteration;
    CHECK(Streamed(reread, nullptr) == meta.dump(2) + "\n");
  }

  // A scratch file that cannot be created reports the failed write.
  JsonArraySpill unwritable((dir / "missing" / "x.spill").string(), keys);
  CHECK_FALSE(unwritable.Append(json::object()));
  std::filesystem::remove_all(dir);
}
