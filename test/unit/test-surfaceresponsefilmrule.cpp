// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceresponsefilmrule.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace
{

using ResponseCorrectionData = config::ElectrostaticSolverData::ResponseCorrectionData;

// The boundary block of the film-rule cases: one attribute per law.
json FilmRuleBoundaries()
{
  return {{"PEC", {{"Attributes", {1}}}},
          {"Impedance",
           {{{"Attributes", {2}}, {"Ls", 4.5e-14}},
            {{"Attributes", {3}}, {"Rs", 1.0}},
            {{"Attributes", {4}}, {"Cs", 4.0e-17}},
            {{"Attributes", {7}}, {"Ls", 2.0e-12}}}},
          {"Conductivity", {{{"Attributes", {5}}, {"Conductivity", 5.8e7}}}},
          {"RationalImpedance",
           {{{"Attributes", {6}}, {"Numerator", {1.0, 0.0}}, {"Denominator", {1.0}}}}}};
}

bool HasLine(const FilmRuleResolution &film, const std::string &text)
{
  return std::any_of(film.lines.begin(), film.lines.end(), [&](const std::string &line)
                     { return line.find(text) != std::string::npos; });
}

}  // namespace

// t1 (STAGEA-NOTE s6): the film rule of decision 2006 and its overrides (ruling 2012).
TEST_CASE("Impedance-metal rule film resolution (decisions 2006 / 2012)",
          "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  const config::BoundaryData boundaries(FilmRuleBoundaries());
  const Units units(1.0, 1.0);
  ResponseCorrectionData request;
  request.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::RECORD;

  SECTION("the default: PEC and a purely inductive sheet are the film; every other law "
          "is not served, named")
  {
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(1).service == FilmService::FILM_PEC);
    CHECK(film.attributes.at(2).service == FilmService::FILM_INDUCTIVE);
    CHECK(film.attributes.at(7).service == FilmService::FILM_INDUCTIVE);  // any Ls
    CHECK(film.attributes.at(3).service == FilmService::NOT_DECLARED);    // Rs only
    CHECK(film.attributes.at(4).service == FilmService::NOT_DECLARED);    // Cs only
    CHECK(film.attributes.at(5).service == FilmService::NOT_DECLARED);    // Conductivity
    CHECK(film.attributes.at(6).service == FilmService::NOT_DECLARED);    // Rational
    CHECK(film.candidate_differs);
    REQUIRE(!film.lines.empty());
    CHECK_THAT(film.lines.front(),
               ContainsSubstring("film = the library's film: an ASSUMPTION"));
    CHECK_THAT(film.lines.front(), ContainsSubstring("PEC [1]"));
    CHECK(HasLine(film, "not served (Conductivity): attributes [5]"));
    CHECK(HasLine(film, "not served (RationalImpedance): attributes [6]"));
    CHECK(HasLine(film, "not served (Impedance Rs 1 Ls 0 Cs 0): attributes [3]"));
    CHECK(HasLine(film, "not served (Impedance Rs 0 Ls 0 Cs 4e-17): attributes [4]"));
    // A segment is served only when every supporting attribute is.
    CHECK(film.SegmentService({1}) == FilmService::FILM_PEC);
    CHECK(film.SegmentService({2}) == FilmService::FILM_INDUCTIVE);
    CHECK(film.SegmentService({1, 2}) == FilmService::FILM_INDUCTIVE);
    CHECK(film.SegmentService({2, 3}) == FilmService::NOT_DECLARED);
    CHECK(film.SegmentService({42}) == FilmService::NOT_DECLARED);
  }

  SECTION("overrides in both directions, any law, each printed")
  {
    request.impedance_metal_film_overrides = {{{3, 5}, true}, {{2}, false}, {{1}, true}};
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(3).service == FilmService::SERVED_BY_OVERRIDE);
    CHECK(film.attributes.at(5).service == FilmService::SERVED_BY_OVERRIDE);
    CHECK(film.attributes.at(2).service == FilmService::NOT_SERVED_BY_OVERRIDE);
    CHECK(film.attributes.at(1).service == FilmService::FILM_PEC);
    CHECK(film.attributes.at(3).overridden);
    CHECK(HasLine(film,
                  "attribute 3 (Impedance Rs 1 Ls 0 Cs 0) served by override: not "
                  "a purely inductive sheet; only the measured check stands behind it"));
    CHECK(HasLine(film, "attribute 5 (Conductivity) served by override: not a purely "
                        "inductive sheet"));
    CHECK(HasLine(film, "attribute 2 (Impedance Rs 0 Ls 4.5e-14 Cs 0) NOT served"));
    CHECK(HasLine(film, "attribute 1 (PEC) served: it is the library's film by default"));
    CHECK(film.SegmentService({2, 1}) == FilmService::NOT_SERVED_BY_OVERRIDE);
    CHECK(film.SegmentService({3}) == FilmService::SERVED_BY_OVERRIDE);
    CHECK(film.SegmentService({1, 3}) == FilmService::SERVED_BY_OVERRIDE);
  }

  SECTION("an override naming an attribute that is not a metal boundary is refused")
  {
    request.impedance_metal_film_overrides = {{{8}, true}};
    CHECK_THROWS_WITH(ResolveImpedanceMetalFilm(boundaries, request, units, false),
                      ContainsSubstring("attribute 8 is not a metal boundary"));
  }

  SECTION("perfect conductors and laws outside the rule only: nothing substituted")
  {
    const config::BoundaryData pec_only(
        json{{"PEC", {{"Attributes", {1}}}},
             {"Impedance", {{{"Attributes", {3}}, {"Rs", 1.0}}}},
             {"Conductivity", {{{"Attributes", {5}}, {"Conductivity", 5.8e7}}}}});
    const auto film = ResolveImpedanceMetalFilm(pec_only, request, units, false);
    CHECK_FALSE(film.candidate_differs);
  }
}

// t1 (STAGEA-NOTE s1): the configuration keys.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule configuration keys (ImpedanceMetalRule, "
                 "ImpedanceMetalFilmOverrides)",
                 "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  json config = ConvexMaxwellIslandConfig();
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  {
    const IoData iodata(config, false);
    CHECK(iodata.solver.surface_response_correction->impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::OFF);
    CHECK(
        iodata.solver.surface_response_correction->impedance_metal_film_overrides.empty());
  }
  correction["ImpedanceMetalRule"] = "Record";
  correction["ImpedanceMetalFilmOverrides"] = {{{"Attributes", {9}}, {"Serve", false}}};
  {
    const IoData iodata(config, false);
    const auto &request = *iodata.solver.surface_response_correction;
    CHECK(request.impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::RECORD);
    REQUIRE(request.impedance_metal_film_overrides.size() == 1);
    CHECK(request.impedance_metal_film_overrides[0].attributes == std::vector<int>{9});
    CHECK_FALSE(request.impedance_metal_film_overrides[0].serve);
  }
  json apply = config;
  apply["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "Apply";
  CHECK_THROWS_WITH(IoData(apply, false), ContainsSubstring("\"Apply\" is reserved"));
  json other = config;
  other["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "On";
  CHECK_THROWS_WITH(IoData(other, false),
                    ContainsSubstring("must be \"Off\" or \"Record\""));
  json off = config;
  off["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "Off";
  CHECK_THROWS_WITH(IoData(off, false),
                    ContainsSubstring("requires \"ImpedanceMetalRule\""));
  json twice = config;
  twice["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalFilmOverrides"] = {
      {{"Attributes", {9}}, {"Serve", true}}, {{"Attributes", {9}}, {"Serve", false}}};
  CHECK_THROWS_WITH(IoData(twice, false), ContainsSubstring("listed more than once"));
  json unknown = config;
  unknown["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalFilmOverrides"] = {
      {{"Attributes", {9}}, {"Serve", true}, {"Basis", "x"}}};
  CHECK_THROWS_WITH(IoData(unknown, false), ContainsSubstring("unknown key \"Basis\""));
}

// t2 / t3 / t5 and the negative laws (STAGEA-NOTE s6): the candidate operator on the convex
// island (attribute 9; Maxwell, Features) against its perfect-conductor twin.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule candidate operator, statuses and law reading",
                 "[surfaceresponsefilmrule][impedancemetalrule][maxwell][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // The island's metal (9; 10 too on the diagonal / neighbouring meshes) with `laws`:
  // {attribute, Impedance json or "PEC"}.
  auto Config = [&](const std::map<int, json> &laws, const fs::path &library)
  {
    json config = ConvexMaxwellIslandConfig();
    std::vector<int> ground = {1, 2, 3, 4, 5, 6};
    json impedance = json::array();
    std::vector<int> metal;
    for (const auto &[attribute, law] : laws)
    {
      metal.push_back(attribute);
      if (law.is_string())
      {
        ground.push_back(attribute);
      }
      else
      {
        json entry = law;
        entry["Attributes"] = {attribute};
        impedance.push_back(entry);
      }
    }
    config["Boundaries"]["Ground"]["Attributes"] = ground;
    if (!impedance.empty())
    {
      config["Boundaries"]["Impedance"] = impedance;
    }
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = metal;
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    correction.erase("PatchConstruction");
    correction["Library"] = library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction["TraceCoupling"] = "Collocated";
    correction["CorrectionMode"] = "PostprocessOnly";
    correction["ImpedanceMetalRule"] = "Record";
    return config;
  };
  struct Run
  {
    std::unique_ptr<IoData> iodata;
    std::vector<std::unique_ptr<Mesh>> meshes;
    std::unique_ptr<SpaceOperator> space;
  };
  // Held by pointer: the space operator keeps references to the IoData and the mesh list.
  auto Build = [&](const json &config, std::unique_ptr<mfem::ParMesh> mesh)
  {
    auto run = std::make_unique<Run>();
    run->iodata = std::make_unique<IoData>(config, false);
    run->iodata->boundaries.cracked_attributes.insert(9);
    run->iodata->boundaries.cracked_attributes.insert(10);
    run->meshes.push_back(std::make_unique<Mesh>(std::move(mesh)));
    run->space = std::make_unique<SpaceOperator>(*run->iodata, run->meshes);
    return run;
  };
  auto Film = [](const Run &run)
  {
    return ResolveImpedanceMetalFilm(
        run.iodata->boundaries, *run.iodata->solver.surface_response_correction,
        run.iodata->units, run.iodata->InputsNondimensionalized());
  };
  // A field with a tangential component on the metal plane y = 0.5 (x and z) and a
  // normal one: the response and the on-metal legs read it.
  auto Field = [](SpaceOperator &space, double ex, double ey, double ez)
  {
    GridFunction field(space.GetNDSpace(), true);
    mfem::Vector value(3);
    value(0) = ex;
    value(1) = ey;
    value(2) = ez;
    mfem::VectorConstantCoefficient coefficient(value);
    field.Real().ProjectCoefficient(coefficient);
    field.Imag() = 0.0;
    return field;
  };
  auto TrueField = [](const GridFunction &field)
  {
    ComplexVector E(field.Real().ParFESpace()->GetTrueVSize());
    field.Real().GetTrueDofs(E.Real());
    field.Imag().GetTrueDofs(E.Imag());
    return E;
  };
  auto Statuses = [](const SurfaceResponseOperator &op)
  {
    std::map<std::string, std::map<std::string, int>> by_type;
    for (const auto &feature : op.GetImpedanceLawFeatures())
    {
      by_type[feature.type][feature.status]++;
    }
    return by_type;
  };
  auto SameResponse = [](const SurfaceResponseOperator::MaxwellResponse &a,
                         const SurfaceResponseOperator::MaxwellResponse &b)
  {
    REQUIRE(a.fabricated_surface_energy.size() == b.fabricated_surface_energy.size());
    for (const auto &[interface, energy] : a.fabricated_surface_energy)
    {
      CHECK(b.fabricated_surface_energy.at(interface) ==
            Catch::Approx(energy).epsilon(1.0e-12).margin(1.0e-300));
    }
    CHECK(b.domain_correction ==
          Catch::Approx(a.domain_correction).epsilon(1.0e-12).margin(1.0e-300));
  };
  const json film_law = {{"Ls", 1.0e-13}};
  // The convex library with every model on the SA interface only (the island's one target
  // interface; the r0a translational-probe cases' form).
  const fs::path sa_library_path = temp.temp_dir / "impedance-rule-convex-sa-3d.json";
  if (Mpi::Root(Mpi::World()))
  {
    json library = json::parse(std::ifstream(convex_library_3d_path));
    library["Name"] = "unit-test-process-impedance-rule-convex-sa-3d";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream(sa_library_path) << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  // The perfect-conductor twin of record (the island PEC).
  auto twin = Build(Config({{9, "PEC"}}, sa_library_path), MakeIslandMesh());
  SurfaceResponseOperator twin_op(*twin->iodata, *twin->space);
  REQUIRE(twin_op.GetPatchCount() > 0);
  const auto twin_field = Field(*twin->space, 0.3, 1.0, 0.7);
  const auto twin_response = twin_op.GetMaxwellResponse(twin_field, 0.0);

  SECTION("RECORD on the perfect conductor: nothing substituted; the law reading on "
          "separate lines leaves every value of the operator unchanged")
  {
    CHECK_FALSE(Film(*twin).candidate_differs);
    SurfaceResponseOperator reading_op(*twin->iodata, *twin->space, nullptr, nullptr, true);
    CHECK(reading_op.LawReadingEnabled());
    CHECK_FALSE(reading_op.IsImpedanceRuleCandidate());
    CHECK(reading_op.GetPatchCount() == twin_op.GetPatchCount());
    const auto response = reading_op.GetMaxwellResponse(twin_field, 0.0);
    for (const auto &[interface, energy] : twin_response.fabricated_surface_energy)
    {
      CHECK(response.fabricated_surface_energy.at(interface) == energy);  // bitwise
    }
    CHECK(response.domain_correction == twin_response.domain_correction);
    CHECK(response.loop_residual == twin_response.loop_residual);
    if (twin_op.HasConductorConsistencyProbes())
    {
      // The gate's records are the same with the reading on (it refuses nothing).
      SurfaceResponseOperator gate_op(*twin->iodata, *twin->space);
      const auto E = TrueField(twin_field);
      const auto records = gate_op.ApplyConductorConsistencyGate(E, 0, 1.0, false);
      const auto reading_records =
          reading_op.ApplyConductorConsistencyGate(E, 0, 1.0, false);
      REQUIRE(records.size() == reading_records.size());
      for (std::size_t i = 0; i < records.size(); i++)
      {
        CHECK(records[i].max_ratio == reading_records[i].max_ratio);
        CHECK(records[i].excluded == reading_records[i].excluded);
      }
      CHECK(gate_op.GetLawReadings().empty());
      const auto &readings = reading_op.GetLawReadings();
      CHECK(static_cast<int>(readings.size()) <= reading_op.GetPatchCount());
      CHECK(!readings.empty());
      // t5: a field normal to the sheet reads zero on every on-metal leg; a tangential
      // field reads its drop.
      const auto normal_field = Field(*twin->space, 0.0, 1.0, 0.0);
      reading_op.ApplyConductorConsistencyGate(TrueField(normal_field), 1, 1.0, false);
      int legs = 0;
      for (const auto &reading : reading_op.GetLawReadings())
      {
        if (reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG)
        {
          legs++;
          CHECK(reading.legs > 0);
          CHECK(reading.max_leg <= 1.0e-12);
        }
      }
      const auto tangential_field = Field(*twin->space, 1.0, 0.0, 1.0);
      reading_op.ApplyConductorConsistencyGate(TrueField(tangential_field), 2, 1.0, false);
      double largest = 0.0;
      for (const auto &reading : reading_op.GetLawReadings())
      {
        if (reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG)
        {
          largest = std::max(largest, reading.max_leg);
        }
      }
      if (legs > 0)
      {
        CHECK(largest > 1.0e-3);
      }
    }
  }

  SECTION("a purely inductive island: the record is raw, the candidate is the PEC twin "
          "(SUBSTITUTED)")
  {
    auto run = Build(Config({{9, film_law}}, sa_library_path), MakeIslandMesh());
    SurfaceResponseOperator record_op(*run->iodata, *run->space);
    CHECK(record_op.IsRawFallback());  // no model of the Impedance law: raw, as today
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.IsImpedanceRuleCandidate());
    REQUIRE_FALSE(candidate.IsRawFallback());
    CHECK(candidate.GetPatchCount() == twin_op.GetPatchCount());
    const auto statuses = Statuses(candidate);
    for (const auto &[type, counts] : statuses)
    {
      CHECK(counts.size() == 1);
      CHECK(counts.count("SUBSTITUTED") == 1);
    }
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      CHECK(feature.device_hash != feature.served_hash);
      CHECK(feature.override_status == "none");
    }
    SameResponse(twin_response,
                 candidate.GetMaxwellResponse(Field(*run->space, 0.3, 1.0, 0.7), 0.0));
  }

  SECTION("the negative laws: Rs only, Cs only, a Conductivity law: not served (nothing "
          "substituted); an override serves them, named")
  {
    for (const json &law :
         {json{{"Rs", 1.0}}, json{{"Cs", 4.0e-17}}, json{{"Rs", 2.0}, {"Ls", 1.0e-13}}})
    {
      auto run = Build(Config({{9, law}}, sa_library_path), MakeIslandMesh());
      const auto film = Film(*run);
      CHECK_FALSE(film.candidate_differs);
      CHECK(film.attributes.at(9).service == FilmService::NOT_DECLARED);
    }
    {
      json config = Config({{9, "PEC"}}, sa_library_path);
      config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
      config["Boundaries"]["Conductivity"] = {
          {{"Attributes", {9}}, {"Conductivity", 5.8e7}, {"Thickness", 1.0e-7}}};
      auto run = Build(config, MakeIslandMesh());
      const auto film = Film(*run);
      CHECK_FALSE(film.candidate_differs);
      CHECK(film.attributes.at(9).service == FilmService::NOT_DECLARED);
      CHECK(HasLine(film, "not served (Conductivity): attributes [9]"));
    }
    json config = Config({{9, json{{"Rs", 1.0}}}}, sa_library_path);
    config["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalFilmOverrides"] = {
        {{"Attributes", {9}}, {"Serve", true}}};
    auto run = Build(config, MakeIslandMesh());
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.GetPatchCount() == twin_op.GetPatchCount());
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      CHECK(feature.status == "SERVED-BY-OVERRIDE");
      CHECK_THAT(feature.override_status,
                 ContainsSubstring("not a purely inductive sheet"));
    }
  }

  SECTION("an override that un-serves a perfect conductor: Missing in the candidate")
  {
    json config = Config({{9, "PEC"}}, sa_library_path);
    config["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalFilmOverrides"] = {
        {{"Attributes", {9}}, {"Serve", false}}};
    auto run = Build(config, MakeIslandMesh());
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.IsRawFallback());
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      CHECK(feature.status == "NOT-SERVED-BY-OVERRIDE");
      CHECK(feature.model.empty());
    }
  }

  SECTION("a law change along one conductor and corners at it: LAW-JOINT, Missing whole "
          "(ERRATUM-1974 E-2)")
  {
    // The diagonal island: attribute 9 below x + z = 1 (PEC), 10 above (the film sheet):
    // one conductor whose law changes at the corners (0.25, 0.75) and (0.75, 0.25).
    auto run = Build(Config({{9, "PEC"}, {10, film_law}}, sa_library_path),
                     MakeDiagonalIslandMesh());
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    const auto statuses = Statuses(candidate);
    REQUIRE(statuses.count("ConvexCorner"));
    const auto &corners = statuses.at("ConvexCorner");
    INFO(
        "LAW-JOINT corners " << (corners.count("LAW-JOINT") ? corners.at("LAW-JOINT") : 0));
    REQUIRE(corners.count("LAW-JOINT"));
    CHECK(corners.at("LAW-JOINT") >= 2);
    CHECK(corners.count("PEC"));
    CHECK(corners.count("SUBSTITUTED"));
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      if (feature.status == "LAW-JOINT")
      {
        CHECK(feature.model.empty());
        CHECK(feature.device_hash == "MIXED");
        CHECK_THAT(feature.device_law, ContainsSubstring("|"));
      }
    }
  }

  SECTION("a mixed library (t3; a TEST-ONLY copy of the isolated model relabelled with "
          "the film law): exactly that key's features EXACT")
  {
    json library = json::parse(std::ifstream(sa_library_path));
    json copy;
    for (const auto &model : library["Models"])
    {
      if (model["Topology"] == "IsolatedEdge")
      {
        copy = model;
      }
    }
    REQUIRE(copy.is_object());
    copy["Name"] = "isolated-film-copy";
    copy["BoundaryCondition"] = {{"Type", "Impedance"}, {"Ls", 1.0e-13}};
    copy["BoundaryLawQualification"] = {{"Version", 1}, {"Status", "Qualified"}};
    library["Name"] = "unit-test-process-impedance-rule-convex-sa-3d-mixed";
    library["Models"].push_back(copy);
    const fs::path mixed_path = temp.temp_dir / "impedance-rule-mixed-library.json";
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream(mixed_path) << library.dump(2) << "\n";
    }
    Mpi::Barrier(Mpi::World());
    auto run = Build(Config({{9, film_law}}, mixed_path), MakeIslandMesh());
    const auto film = Film(*run);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.GetPatchCount() == twin_op.GetPatchCount());
    int exact = 0;
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      if (feature.type == "IsolatedEdge")
      {
        CHECK(feature.status == "EXACT");
        CHECK(feature.model == "isolated-film-copy");
        CHECK(feature.device_hash != feature.served_hash);
        exact++;
      }
      else
      {
        CHECK(feature.status == "SUBSTITUTED");
      }
    }
    CHECK(exact > 0);
  }
#endif
}
