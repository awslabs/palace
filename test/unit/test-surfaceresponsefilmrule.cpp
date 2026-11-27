// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unistd.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "models/surfaceresponsefilmrule.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/metadatafile.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"
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

// The text a callable writes to the process's standard output (the root's log).
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

// Decision 2206 (a) (DESIGN-APPLY s2 / s4 / s5): "Apply" serves PEC and purely inductive
// sheets within the Ls range [0, ImpedanceMetalMaxInductance], names each sheet's form,
// and refuses an override that would serve anything else.
TEST_CASE("Impedance-metal rule Apply film resolution, Ls range check and refusals",
          "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  config::BoundaryData boundaries(FilmRuleBoundaries());
  const Units units(1.0, 1.0);
  ResponseCorrectionData request;
  request.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::APPLY;

  SECTION("the default bound 9e-14 H/sq: the designers' 4.5e-14 served, 2e-12 not")
  {
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.apply);
    CHECK(film.max_inductance == 9.0e-14);
    CHECK(film.max_inductance_is_default);
    CHECK(film.attributes.at(1).service == FilmService::FILM_PEC);
    CHECK(film.attributes.at(2).service == FilmService::FILM_INDUCTIVE);
    CHECK(film.attributes.at(7).service == FilmService::INDUCTIVE_ABOVE_RANGE);
    CHECK(film.attributes.at(3).service == FilmService::NOT_DECLARED);  // Rs only
    CHECK(film.attributes.at(4).service == FilmService::NOT_DECLARED);  // Cs only
    CHECK(film.attributes.at(5).service == FilmService::NOT_DECLARED);  // Conductivity
    CHECK(film.attributes.at(6).service == FilmService::NOT_DECLARED);  // Rational
    CHECK(std::string(FilmServiceName(FilmService::INDUCTIVE_ABOVE_RANGE)) ==
          "ABOVE-RANGE");
    CHECK(film.served_inductances == std::vector<double>{4.5e-14});
    REQUIRE(!film.lines.empty());
    const auto &line = film.lines.front();
    CHECK_THAT(line, ContainsSubstring("Impedance-metal rule (APPLY): treated as the "
                                       "library's perfect-conductor film (an ASSUMPTION)"));
    CHECK_THAT(line, ContainsSubstring("PEC [1]"));
    CHECK_THAT(line, ContainsSubstring("Impedance Rs 0 Ls 4.5e-14 Cs 0, uncracked (one "
                                       "sheet of total Ls) [2]"));
    CHECK_THAT(line, ContainsSubstring("Ls range served: 0 < configured Ls <= 9e-14 H/sq "
                                       "in every sheet form (ImpedanceMetalMaxInductance, "
                                       "the default)"));
    // CUT-PLAN step 3: the regime label replaced the line "experimental and not validated":
    // ONE regime line with the bound and the validated range; no line says "not validated".
    CHECK(std::none_of(film.lines.begin(), film.lines.end(), [](const std::string &text)
                       { return text.find("not validated") != std::string::npos; }));
    CHECK(std::count_if(film.lines.begin(), film.lines.end(), [](const std::string &text)
                        { return text.find("(APPLY): regime: ") != std::string::npos; }) ==
          1);
    CHECK(HasLine(film, "Impedance-metal rule (APPLY): regime: purely inductive sheets at "
                        "configured Ls <= 9e-14 H/sq are served by the perfect-conductor "
                        "models; validated at the designers' 4.5e-14 H/sq in driven and "
                        "eigenmode runs, cracked and uncracked"));
    CHECK(HasLine(film, kImpedanceApplyValidatedRange));
    CHECK_THAT(line, ContainsSubstring("served Ls met: 4.5e-14 H/sq"));
    // Not served (as ruled for the default: raw and named, no abort): the line names the
    // bound and the validated range in H/sq.
    CHECK(HasLine(film,
                  "not served (Impedance Rs 0 Ls 2e-12 Cs 0, uncracked (one sheet "
                  "of total Ls): ABOVE the Ls range <= 9e-14 H/sq): attributes [7]: "
                  "not served by the perfect-conductor models: raw unless the library "
                  "holds a model of this law (law status EXACT, corrected by that "
                  "model); the perfect-conductor models serve Ls sheets at or below "
                  "9e-14 H/sq (validated at the designers' 4.5e-14 H/sq"));
    CHECK(HasLine(film, "not served (Impedance Rs 1 Ls 0 Cs 0, uncracked (one sheet of "
                        "total Ls)): attributes [3]: not served by the perfect-conductor "
                        "models: raw unless the library holds a model of this law"));
    CHECK(HasLine(film, "not served (Conductivity): attributes [5]"));
    // An above-range sheet serves no segment; a served one does.
    CHECK(film.SegmentService({7}) == FilmService::NOT_DECLARED);
    CHECK(film.SegmentService({2, 7}) == FilmService::NOT_DECLARED);
    CHECK(film.SegmentService({1, 2}) == FilmService::FILM_INDUCTIVE);
  }

  SECTION("at, below and above a configured bound")
  {
    request.impedance_metal_max_inductance = 4.5e-14;
    request.impedance_metal_max_inductance_given = true;
    auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(2).service == FilmService::FILM_INDUCTIVE);  // at the bound
    CHECK_FALSE(film.max_inductance_is_default);
    CHECK_THAT(film.lines.front(), ContainsSubstring("(ImpedanceMetalMaxInductance, "
                                                     "configured)"));
    request.impedance_metal_max_inductance = 4.4e-14;
    film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(2).service == FilmService::INDUCTIVE_ABOVE_RANGE);
    CHECK_THAT(film.lines.front(), ContainsSubstring("served Ls met: none"));
    request.impedance_metal_max_inductance = 3.0e-12;
    film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(7).service == FilmService::FILM_INDUCTIVE);
    CHECK(film.served_inductances == std::vector<double>{4.5e-14, 2.0e-12});
    CHECK_THAT(film.lines.front(), ContainsSubstring("served Ls met: 4.5e-14, 2e-12 H/sq"));
  }

  // Lane v1lsdefault (bench-driven-ls NOTES 21): the bound is inclusive in the H/sq the
  // designer configures. Main compared the scaled Ls times its scale, fl(fl(Ls / s) s): at
  // this characteristic length (tools/roundtrip_scan_v1.py of the lane) that lands one unit
  // in the last place above 9e-14, so a sheet at exactly the bound was refused.
  SECTION("the bound is inclusive in the configured H/sq: at it, one ulp below, one above")
  {
    const Units scaled(1.0e-6, 1.0011520131719906e-6);
    const double s = scaled.GetScaleFactor<Units::ValueType::INDUCTANCE>();
    const double bound = ResponseCorrectionData::kDefaultImpedanceMetalMaxInductance;
    REQUIRE(bound == 9.0e-14);
    REQUIRE((bound / s) * s > bound);  // the premise: main's round trip refused the bound
    const double below = std::nextafter(bound, 0.0);
    const double above = std::nextafter(bound, 1.0);
    const json sheets = {{"Impedance",
                          {{{"Attributes", {2}}, {"Ls", bound}},
                           {{"Attributes", {3}}, {"Ls", below}},
                           {{"Attributes", {4}}, {"Ls", above}}}}};
    config::BoundaryData si(sheets), nondimensional(sheets);
    for (auto &data : nondimensional.impedance)
    {
      config::Nondimensionalize(scaled, data);
    }
    CHECK(nondimensional.impedance[0].configured_Ls == bound);
    CHECK(nondimensional.impedance[0].Ls == bound / s);
    for (const bool scaled_data : {true, false})
    {
      INFO("nondimensionalized " << scaled_data);
      const auto film = ResolveImpedanceMetalFilm(scaled_data ? nondimensional : si,
                                                  request, scaled, scaled_data);
      CHECK(film.attributes.at(2).service == FilmService::FILM_INDUCTIVE);
      CHECK(film.attributes.at(3).service == FilmService::FILM_INDUCTIVE);
      CHECK(film.attributes.at(4).service == FilmService::INDUCTIVE_ABOVE_RANGE);
      // The recorded Ls values are the configured numbers.
      CHECK(film.attributes.at(2).inductance == bound);
      CHECK(film.served_inductances == std::vector<double>{below, bound});
      // The configured Ls is printed exactly (review v1lsdefault m3): the sheet one ulp
      // above the bound does not read as the bound in its not-served line, and the served
      // Ls list tells one ulp below from the bound.
      REQUIRE(fmt::format("{}", above) == "9.000000000000001e-14");
      REQUIRE(fmt::format("{}", below) == "8.999999999999998e-14");
      CHECK(HasLine(film, "not served (Impedance Rs 0 Ls 9.000000000000001e-14 Cs 0, "
                          "uncracked (one sheet of total Ls): ABOVE the Ls range <= 9e-14 "
                          "H/sq): attributes [4]"));
      CHECK_FALSE(HasLine(film, "Ls 9e-14 Cs 0, uncracked (one sheet of total Ls): ABOVE"));
      CHECK_THAT(film.lines.front(),
                 ContainsSubstring("served Ls met: 8.999999999999998e-14, 9e-14 H/sq"));
    }
  }

  SECTION("the sheet form of each attribute: uncracked, cracked, the legacy split")
  {
    boundaries.cracked_attributes = {2, 7};
    boundaries.split_impedance_attributes = {7};
    request.impedance_metal_max_inductance = 3.0e-12;
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(2).form == "cracked (two faces of Ls each)");
    CHECK(film.attributes.at(7).form == "cracked, legacy split (two faces of 2 Ls each)");
    CHECK_THAT(film.lines.front(),
               ContainsSubstring("Ls 4.5e-14 Cs 0, cracked (two faces of Ls each) [2]"));
    CHECK_THAT(film.lines.front(),
               ContainsSubstring("Ls 2e-12 Cs 0, cracked, legacy split (two faces of 2 Ls "
                                 "each) [7]"));
    // The range check compares the CONFIGURED Ls in every form (decision 2208 (3)): the
    // same 2e-12 is served cracked and split alike under a 3e-12 bound.
    CHECK(film.attributes.at(7).inductance == 2.0e-12);
    // Record keeps the branch's lines (no form, no range).
    request.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::RECORD;
    const auto record = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(record.attributes.at(7).service == FilmService::FILM_INDUCTIVE);
    CHECK_THAT(record.lines.front(),
               ContainsSubstring("Impedance-metal rule (RECORD): film = the library's "
                                 "film: an ASSUMPTION: Impedance Rs 0 Ls 2e-12 Cs 0 [7]"));
  }

  SECTION("overrides: un-serve any metal; serve only what the rule serves by default")
  {
    request.impedance_metal_film_overrides = {{{2}, false}, {{1}, true}};
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    CHECK(film.attributes.at(2).service == FilmService::NOT_SERVED_BY_OVERRIDE);
    CHECK(film.attributes.at(1).service == FilmService::FILM_PEC);
    CHECK(film.served_inductances.empty());
    CHECK(HasLine(film, "attribute 1 (PEC) served: it is the library's film by default"));
    request.impedance_metal_film_overrides = {{{3}, true}};
    CHECK_THROWS_WITH(
        ResolveImpedanceMetalFilm(boundaries, request, units, false),
        ContainsSubstring("attribute 3 (Impedance Rs 1 Ls 0 Cs 0) cannot be served by an "
                          "override under \"ImpedanceMetalRule\": \"Apply\": Apply serves "
                          "perfect conductors and purely inductive sheets (Rs = 0, Cs = 0) "
                          "only"));
    request.impedance_metal_film_overrides = {{{5}, true}};
    CHECK_THROWS_WITH(ResolveImpedanceMetalFilm(boundaries, request, units, false),
                      ContainsSubstring("attribute 5 (Conductivity) cannot be served"));
    request.impedance_metal_film_overrides = {{{7}, true}};
    CHECK_THROWS_WITH(ResolveImpedanceMetalFilm(boundaries, request, units, false),
                      ContainsSubstring("a larger ImpedanceMetalMaxInductance serves it"));
  }

  SECTION("the lower end: Ls = 0 with Rs = Cs = 0 is no sheet: not served, named; an "
          "override cannot serve it; Record keeps it")
  {
    config::BoundaryData zero(json{
        {"PEC", {{"Attributes", {1}}}},
        {"Impedance",
         {{{"Attributes", {8}}, {"Ls", 0.0}}, {{"Attributes", {2}}, {"Ls", 1.0e-15}}}}});
    const auto film = ResolveImpedanceMetalFilm(zero, request, units, false);
    CHECK(film.attributes.at(8).service == FilmService::NOT_DECLARED);
    CHECK(film.attributes.at(8).no_sheet);
    CHECK(film.attributes.at(2).service == FilmService::FILM_INDUCTIVE);
    CHECK_FALSE(film.attributes.at(2).no_sheet);
    CHECK(film.served_inductances == std::vector<double>{1.0e-15});
    CHECK_THAT(film.lines.front(), ContainsSubstring("0 < configured Ls <= 9e-14 H/sq"));
    CHECK(HasLine(film,
                  "not served (Impedance Rs 0 Ls 0 Cs 0, uncracked (one sheet of total "
                  "Ls): no sheet (Rs, Ls and Cs are all 0)): attributes [8]"));
    CHECK(film.SegmentService({8}) == FilmService::NOT_DECLARED);
    request.impedance_metal_film_overrides = {{{8}, true}};
    CHECK_THROWS_WITH(
        ResolveImpedanceMetalFilm(zero, request, units, false),
        ContainsSubstring("attribute 8 (Impedance Rs 0 Ls 0 Cs 0) cannot be "
                          "served by an override under \"ImpedanceMetalRule\": "
                          "\"Apply\": it is no sheet (Rs, Ls and Cs are all 0)"));
    request.impedance_metal_film_overrides.clear();
    request.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::RECORD;
    const auto record = ResolveImpedanceMetalFilm(zero, request, units, false);
    CHECK(record.attributes.at(8).service == FilmService::FILM_INDUCTIVE);
    CHECK_FALSE(record.attributes.at(8).no_sheet);
  }

  SECTION("the palace.json / manifest record of an Apply run")
  {
    const auto film = ResolveImpedanceMetalFilm(boundaries, request, units, false);
    // n 4 (the rows on an impedance sheet), s 3 (their law status serves the sheet), k 2
    // (corrected by a model), m 1 (SUBSTITUTED without a model: Missing).
    std::vector<ResponseCorrectionData::ImpedanceLawFeatureData> features(5);
    features[0].status = "SUBSTITUTED";
    features[0].impedance_sheet = true;
    features[0].model = "isolated";
    features[0].length = 2.0;
    features[1].status = "PEC";
    features[2].status = "MIXED-LAW";
    features[2].impedance_sheet = true;
    features[2].model = "convex-corner-90";
    features[2].length = 0.5;
    features[3].status = "SUBSTITUTED";
    features[3].impedance_sheet = true;
    features[3].length = 4.0;
    features[4].status = "NOT-DECLARED";
    features[4].impedance_sheet = true;
    features[4].length = 8.0;
    const auto counts = CountImpedanceApplyFeatures(features);
    CHECK(counts.sheet == 4);
    CHECK(counts.substituted == 3);
    CHECK(counts.corrected == 2);
    CHECK(counts.missing == 1);
    CHECK(counts.corrected_length == 2.5);
    // CUT-PLAN step 3: the features on a metal the rule does not serve (the not-served
    // flags word of surface-participation.csv): NOT-DECLARED here; LAW-JOINT and
    // NOT-SERVED-BY-OVERRIDE below.
    CHECK(counts.not_served == 1);
    auto joint = features;
    joint[1].status = "LAW-JOINT";
    joint[3].status = "NOT-SERVED-BY-OVERRIDE";
    CHECK(CountImpedanceApplyFeatures(joint).not_served == 3);
    CHECK(CountImpedanceApplyFeatures(
              std::vector<ResponseCorrectionData::ImpedanceLawFeatureData>(1, features[0]))
              .not_served == 0);
    CHECK(ImpedanceApplyCountText(counts) ==
          "2 of 4 feature(s) on impedance sheets are corrected by perfect-conductor models "
          "(1 Missing, as for a perfect conductor)");
    // ServedImpedanceLength in metres: the rows' nondimensional length times Lc.
    const Units scaled_units(1.0e-6, 2.0e-3);
    const auto record = ImpedanceApplyRecord(film, features, scaled_units, false);
    CHECK(record.at("Mode") == "Apply");
    CHECK(record.at("ImpedanceSheetFeatures") == 4);
    CHECK(record.at("SubstitutedImpedanceFeatures") == 3);
    CHECK(record.at("ServedImpedanceFeatures") == 2);
    CHECK(record.at("MissingImpedanceFeatures") == 1);
    CHECK(record.at("ServedImpedanceLength").get<double>() ==
          Catch::Approx(2.5 * 2.0e-3).epsilon(1.0e-14));
    CHECK(record.at("LawStatusCounts").at("PEC") == 1);
    CHECK(record.at("LawStatusCounts").at("SUBSTITUTED") == 2);
    CHECK(record.at("MaxInductance").get<double>() == 9.0e-14);
    CHECK(record.at("MaxInductanceIsDefault") == true);
    CHECK(record.at("Lines").size() == film.lines.size());
    CHECK_FALSE(record.contains("Features"));
    CHECK(ImpedanceApplyRecord(film, features, units, true).at("Features").size() == 5);
    // An all-perfect-conductor run: no feature on an impedance sheet (n 0).
    const std::vector<ResponseCorrectionData::ImpedanceLawFeatureData> pec(2, features[1]);
    CHECK(ImpedanceApplyCountText(CountImpedanceApplyFeatures(pec)) ==
          "0 of 0 feature(s) on impedance sheets are corrected by perfect-conductor models "
          "(0 Missing, as for a perfect conductor)");
  }

  // CUT-PLAN step 3 ("Apply" is the Maxwell default): the rule acts only where there is a
  // sheet to serve or an override; a perfect-conductor run is "Off" by construction.
  SECTION("Apply acts only on a run with an Impedance boundary or an override")
  {
    CHECK(ImpedanceApplyActs(boundaries, request));
    const config::BoundaryData no_sheet(
        json{{"PEC", {{"Attributes", {1}}}},
             {"Conductivity", {{{"Attributes", {5}}, {"Conductivity", 5.8e7}}}}});
    CHECK_FALSE(ImpedanceApplyActs(no_sheet, request));
    auto overridden = request;
    overridden.impedance_metal_film_overrides = {{{1}, false}};
    CHECK(ImpedanceApplyActs(no_sheet, overridden));
    for (const auto rule : {ResponseCorrectionData::ImpedanceMetalRule::OFF,
                            ResponseCorrectionData::ImpedanceMetalRule::RECORD})
    {
      auto other = request;
      other.impedance_metal_rule = rule;
      CHECK_FALSE(ImpedanceApplyActs(boundaries, other));
    }
  }
}

// Decision 2215: an override refused under "Apply" stops every rank (the refusal follows
// from the configuration alone); the root carries the reason, the other ranks one short
// line.
TEST_CASE("Impedance-metal rule Apply override refusal names its reason on the root only",
          "[surfaceresponsefilmrule][impedancemetalrule][Serial][Parallel]")
{
  config::BoundaryData boundaries(FilmRuleBoundaries());
  const Units units(1.0, 1.0);
  ResponseCorrectionData request;
  request.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::APPLY;
  request.impedance_metal_film_overrides = {{{3}, true}};
  const std::string reason = "attribute 3 (Impedance Rs 1 Ls 0 Cs 0) cannot be served";
  const std::string pointer = "(the reason is printed by the root process)";
  if (Mpi::Root(Mpi::World()))
  {
    CHECK_THROWS_WITH(ResolveImpedanceMetalFilm(boundaries, request, units, false),
                      ContainsSubstring(reason) && !ContainsSubstring(pointer));
  }
  else
  {
    CHECK_THROWS_WITH(ResolveImpedanceMetalFilm(boundaries, request, units, false),
                      ContainsSubstring(pointer) && !ContainsSubstring(reason));
  }
}

// Decision 2208 / DESIGN-APPLY s8: an electrostatic run solves no impedance boundary, so
// "Apply" serves nothing differently there: the operator is the one of "Off" (no seam, the
// same patches, the same corrected response).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule Apply in an electrostatic run serves nothing "
                 "differently",
                 "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  auto Response = [&](const std::string &rule)
  {
    json config = ConvexIslandConfig();
    config["Solver"]["Electrostatic"]["ResponseCorrection"]["ImpedanceMetalRule"] = rule;
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    LaplaceOperator laplace(iodata, meshes);
    auto K = laplace.GetStiffnessMatrix();
    KspSolver ksp(iodata, laplace.GetH1Spaces());
    ksp.SetOperators(*K, *K);
    Vector V, RHS;
    laplace.GetExcitationVector(1, *K, V, RHS);
    ksp.Mult(RHS, V);
    SurfaceResponseOperator response(iodata, laplace);
    CHECK(response.GetImpedanceApplyRule() == nullptr);
    CHECK_FALSE(response.IsImpedanceRuleCandidate());
    return std::pair{response.GetPatchCount(), response.GetElectrostaticResponse(V)};
  };
  const auto [off_patches, off] = Response("Off");
  const auto [apply_patches, apply] = Response("Apply");
  REQUIRE(off_patches > 0);
  CHECK(apply_patches == off_patches);
  CHECK(apply.domain_correction == off.domain_correction);
  REQUIRE(apply.fabricated_surface_energy.size() == off.fabricated_surface_energy.size());
  for (const auto &[interface, energy] : off.fabricated_surface_energy)
  {
    CHECK(apply.fabricated_surface_energy.at(interface) == energy);
  }
#endif
}

// Review r2206 MINOR-4 (decision 2215 (B)): the Legacy patch construction has no film seam,
// so an "Apply" preflight there writes the manifest of "Off" (its requirements carry the
// device law, the law of the run).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule Apply preflight on the Legacy patch construction "
                 "writes the manifest of Off",
                 "[surfaceresponsefilmrule][impedancemetalrule][maxwell][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  Mesh island_mesh(MakeIslandMesh());
  auto StripWallTime = [](auto &self, json &node) -> void
  {
    if (node.is_object())
    {
      node.erase("WallTime");
    }
    if (node.is_object() || node.is_array())
    {
      for (auto &child : node)
      {
        self(self, child);
      }
    }
  };
  auto Preflight = [&](const std::string &rule)
  {
    json config = ConvexMaxwellIslandConfig();  // PatchConstruction "Legacy"
    config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
    config["Boundaries"]["Impedance"] =
        json::array({{{"Attributes", {9}}, {"Ls", 4.5e-14}}});
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    REQUIRE(correction.at("PatchConstruction") == "Legacy");
    correction["UnmatchedPolicy"] = "Warn";
    correction["ImpedanceMetalRule"] = rule;
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto dir = temp.temp_dir / ("impedance-legacy-preflight-" + rule);
    fs::create_directories(dir);
    const auto manifest_path = dir / "surface-response-requirements.json";
    const std::string log = CaptureStdout(
        dir / "preflight.log", [&]
        { WriteSurfaceResponseRequirements(iodata, island_mesh, manifest_path.string()); });
    std::ifstream input(manifest_path);
    REQUIRE(input);
    json manifest = json::parse(input);
    StripWallTime(StripWallTime, manifest);
    return std::pair<json, std::string>(std::move(manifest), log);
  };
  const auto [off, off_log] = Preflight("Off");
  const auto [apply, apply_log] = Preflight("Apply");
  CHECK(apply == off);
  CHECK_FALSE(apply.contains("ImpedanceMetalRule"));
  CHECK_THAT(apply_log,
             ContainsSubstring("Impedance-metal rule (APPLY): no seam on this path "
                               "(PatchConstruction Legacy): as \"Off\""));
  CHECK_THAT(off_log, !ContainsSubstring("Impedance-metal rule"));
  // The requirements carry the device law (the Impedance sheet), as under "Off".
  bool impedance_requirement = false;
  for (const auto &requirement : apply.at("Requirements"))
  {
    impedance_requirement = impedance_requirement ||
                            requirement.at("BoundaryCondition").at("Type") == "Impedance";
  }
  CHECK(impedance_requirement);
#endif
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
    // CUT-PLAN step 3: "Apply" by default in the Maxwell block, "Off" in the
    // electrostatic block (ConvexIslandConfig).
    const IoData iodata(config, false);
    CHECK(iodata.solver.surface_response_correction->impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::APPLY);
    CHECK(
        iodata.solver.surface_response_correction->impedance_metal_film_overrides.empty());
    const IoData electrostatic(ConvexIslandConfig(), false);
    CHECK(electrostatic.solver.electrostatic.response_correction->impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::OFF);
    CHECK(ResponseCorrectionData().impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::OFF);  // a request built in code
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
  // "Apply" (decision 2206 (a)): accepted, with its bound (default 9e-14 H/sq) and the
  // overrides; the bound is refused without "Apply" and when not a positive number.
  json apply = config;
  apply["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "Apply";
  {
    const IoData iodata(apply, false);
    const auto &request = *iodata.solver.surface_response_correction;
    CHECK(request.impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::APPLY);
    CHECK(request.impedance_metal_max_inductance == 9.0e-14);
    CHECK_FALSE(request.impedance_metal_max_inductance_given);
    REQUIRE(request.impedance_metal_film_overrides.size() == 1);
  }
  apply["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalMaxInductance"] = 2.5e-13;
  {
    const IoData iodata(apply, false);
    const auto &request = *iodata.solver.surface_response_correction;
    CHECK(request.impedance_metal_max_inductance == 2.5e-13);
    CHECK(request.impedance_metal_max_inductance_given);
  }
  for (const json &bad : {json(0.0), json(-1.0e-14), json("9e-14")})
  {
    json invalid = apply;
    invalid["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalMaxInductance"] = bad;
    CHECK_THROWS_WITH(IoData(invalid, false),
                      ContainsSubstring("must be a finite positive surface inductance"));
  }
  json bound_record = config;
  bound_record["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalMaxInductance"] =
      1.0e-13;
  CHECK_THROWS_WITH(IoData(bound_record, false),
                    ContainsSubstring("\"ImpedanceMetalMaxInductance\" requires "
                                      "\"ImpedanceMetalRule\": \"Apply\""));
  // "Apply" writes no candidate value: it runs under either coverage policy.
  json apply_marked = apply;
  apply_marked["Solver"]["SurfaceResponseCorrection"]["Coverage"]["Policy"] =
      "MarkNotApplied";
  CHECK_NOTHROW(IoData(apply_marked, false));
  json other = config;
  other["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "On";
  CHECK_THROWS_WITH(IoData(other, false),
                    ContainsSubstring("must be \"Off\", \"Record\" or \"Apply\""));
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
  // Ruling 2065 (1): "Record" runs only under Coverage.Policy "Warn" (the measuring mode
  // prints candidate values beside the values of record; MarkNotApplied withholds them).
  const std::string policy_refusal =
      "\"Record\" runs only with \"Coverage\": {\"Policy\": \"Warn\"}";
  json marked = config;
  marked["Solver"]["SurfaceResponseCorrection"]["Coverage"]["Policy"] = "MarkNotApplied";
  CHECK_THROWS_WITH(IoData(marked, false), ContainsSubstring(policy_refusal));
  json warned = config;
  warned["Solver"]["SurfaceResponseCorrection"]["Coverage"]["Policy"] = "Warn";
  CHECK_NOTHROW(IoData(warned, false));
  json marked_off = marked;
  marked_off["Solver"]["SurfaceResponseCorrection"].erase("ImpedanceMetalFilmOverrides");
  marked_off["Solver"]["SurfaceResponseCorrection"]["ImpedanceMetalRule"] = "Off";
  CHECK_NOTHROW(IoData(marked_off, false));
  {
    // The refusal reads the RESOLVED policy: with the key absent the parse resolves
    // CoverageData::kDefaultCoveragePolicy (a constant: not patchable here), so the check
    // is called with each value the key-absent parse can resolve to, the flipped default
    // MARK_NOT_APPLIED included.
    const IoData absent(config, false);
    const auto &request = *absent.solver.surface_response_correction;
    CHECK(request.coverage.policy ==
          ResponseCorrectionData::CoverageData::kDefaultCoveragePolicy);
    ResponseCorrectionData flipped = request;
    flipped.coverage.policy =
        ResponseCorrectionData::CoverageData::Policy::MARK_NOT_APPLIED;
    CHECK_THROWS_WITH(config::VerifyImpedanceMetalRuleCoveragePolicy(flipped),
                      ContainsSubstring(policy_refusal));
    flipped.coverage.policy = ResponseCorrectionData::CoverageData::Policy::WARN;
    CHECK_NOTHROW(config::VerifyImpedanceMetalRuleCoveragePolicy(flipped));
    flipped.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::OFF;
    flipped.coverage.policy =
        ResponseCorrectionData::CoverageData::Policy::MARK_NOT_APPLIED;
    CHECK_NOTHROW(config::VerifyImpedanceMetalRuleCoveragePolicy(flipped));
    flipped.impedance_metal_rule = ResponseCorrectionData::ImpedanceMetalRule::APPLY;
    CHECK_NOTHROW(config::VerifyImpedanceMetalRuleCoveragePolicy(flipped));
  }
}

// CUT-PLAN step 3 (lane v1lsdefault; the supervisor's rulings): "Record",
// "ImpedanceMetalFilmOverrides" and a bound above the validated 9e-14 H/sq are developer
// options; at or below it any positive bound is accepted (a stricter one).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule developer options and the validated bound",
                 "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  json config = ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction.erase("Developer");
  {
    const IoData iodata(config, false);
    const auto &request = *iodata.solver.surface_response_correction;
    CHECK(request.impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::APPLY);
    CHECK(request.developer_options.empty());
  }
  const double bound = ResponseCorrectionData::kDefaultImpedanceMetalMaxInductance;
  auto WithKeys = [&](const json &keys, bool developer)
  {
    json run = config;
    for (const auto &[key, value] : keys.items())
    {
      run["Solver"]["SurfaceResponseCorrection"][key] = value;
    }
    if (developer)
    {
      run["Solver"]["SurfaceResponseCorrection"]["Developer"] = true;
    }
    return run;
  };
  const std::string switch_text = "is a parked developer option";
  const json record = {{"ImpedanceMetalRule", "Record"},
                       {"Coverage", {{"Policy", "Warn"}}}};
  const json overrides = {
      {"ImpedanceMetalFilmOverrides", {{{"Attributes", {9}}, {"Serve", false}}}}};
  const json larger = {{"ImpedanceMetalMaxInductance", std::nextafter(bound, 1.0)}};
  CHECK_THROWS_WITH(IoData(WithKeys(record, false), false),
                    ContainsSubstring("\"ImpedanceMetalRule\": \"Record\" " + switch_text));
  CHECK_THROWS_WITH(IoData(WithKeys(overrides, false), false),
                    ContainsSubstring("\"ImpedanceMetalFilmOverrides\" " + switch_text));
  CHECK_THROWS_WITH(IoData(WithKeys(larger, false), false),
                    ContainsSubstring("\"ImpedanceMetalMaxInductance\": "
                                      "9.000000000000001e-14 H/sq (above the validated "
                                      "bound 9e-14 H/sq) " +
                                      switch_text));
  for (const auto &[keys, option] :
       {std::pair{record, "\"ImpedanceMetalRule\": \"Record\""},
        std::pair{overrides, "\"ImpedanceMetalFilmOverrides\""},
        std::pair{larger, "\"ImpedanceMetalMaxInductance\": 9.000000000000001e-14 H/sq"}})
  {
    INFO(keys.dump());
    const IoData iodata(WithKeys(keys, true), false);
    const auto &options = iodata.solver.surface_response_correction->developer_options;
    REQUIRE(options.size() == 1);
    CHECK_THAT(options[0], ContainsSubstring(option));
  }
  // At or below the validated bound: accepted without the switch, as a stricter bound.
  for (const double value : {bound, std::nextafter(bound, 0.0), 4.5e-14})
  {
    const IoData iodata(WithKeys({{"ImpedanceMetalMaxInductance", value}}, false), false);
    const auto &request = *iodata.solver.surface_response_correction;
    CHECK(request.impedance_metal_max_inductance == value);
    CHECK(request.impedance_metal_max_inductance_given);
    CHECK(request.developer_options.empty());
  }
  // The bound still needs "Apply": refused under an explicit "Off" and in the
  // electrostatic block, whose default stays "Off".
  CHECK_THROWS_WITH(IoData(WithKeys({{"ImpedanceMetalRule", "Off"},
                                     {"ImpedanceMetalMaxInductance", bound}},
                                    false),
                           false),
                    ContainsSubstring("requires \"ImpedanceMetalRule\": \"Apply\""));
  json electrostatic = ConvexIslandConfig(test::CrackPolicy::Default);
  electrostatic["Solver"]["Electrostatic"]["ResponseCorrection"]
               ["ImpedanceMetalMaxInductance"] = bound;
  CHECK_THROWS_WITH(IoData(electrostatic, false),
                    ContainsSubstring("requires \"ImpedanceMetalRule\": \"Apply\""));
  // "Off" stays allowed in the Maxwell block.
  const IoData off(WithKeys({{"ImpedanceMetalRule", "Off"}}, false), false);
  CHECK(off.solver.surface_response_correction->impedance_metal_rule ==
        ResponseCorrectionData::ImpedanceMetalRule::OFF);
  // CUT-PLAN step 4b (lane v1step4b) with step 3: PreflightOnly refuses and accepts these
  // keys exactly as the solve form, the same message without "Developer" and the same
  // developer options with it.
  for (const auto &keys : {record, overrides, larger})
  {
    INFO(keys.dump());
    json preflight_keys = keys;
    preflight_keys["PreflightOnly"] = true;
    auto Refusal = [&](const json &run)
    {
      try
      {
        IoData iodata(run, false);
      }
      catch (const std::exception &error)
      {
        return std::string(error.what());
      }
      return std::string();
    };
    const std::string refusal = Refusal(WithKeys(keys, false));
    CHECK_THAT(refusal, ContainsSubstring(switch_text));
    CHECK(Refusal(WithKeys(preflight_keys, false)) == refusal);
    const IoData solve(WithKeys(keys, true), false);
    const IoData preflight(WithKeys(preflight_keys, true), false);
    CHECK(preflight.solver.surface_response_correction->preflight_only);
    CHECK(preflight.solver.surface_response_correction->developer_options ==
          solve.solver.surface_response_correction->developer_options);
  }
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

  SECTION("RECORD on the perfect conductor: the candidate is a twin of the operator of "
          "record (ruling 2024); the law reading runs on the twin only")
  {
    const auto film = Film(*twin);
    CHECK_FALSE(film.candidate_differs);
    SurfaceResponseOperator candidate(*twin->iodata, *twin->space, nullptr, &film, true);
    CHECK(candidate.LawReadingEnabled());
    CHECK(candidate.IsImpedanceRuleCandidate());
    CHECK_FALSE(candidate.ImpedanceRuleActs());
    CHECK_FALSE(twin_op.LawReadingEnabled());
    CHECK(candidate.GetPatchCount() == twin_op.GetPatchCount());
    const auto response = candidate.GetMaxwellResponse(twin_field, 0.0);
    for (const auto &[interface, energy] : twin_response.fabricated_surface_energy)
    {
      CHECK(response.fabricated_surface_energy.at(interface) == energy);  // bitwise
    }
    CHECK(response.domain_correction == twin_response.domain_correction);
    CHECK(response.loop_residual == twin_response.loop_residual);
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      CHECK(feature.status == "PEC");
      CHECK(feature.device_hash == feature.served_hash);
    }
    if (twin_op.HasConductorConsistencyProbes())
    {
      // The gate's records are the same on the twin with the reading on (it refuses
      // nothing); the operator of record takes no reading.
      SurfaceResponseOperator gate_op(*twin->iodata, *twin->space);
      const auto E = TrueField(twin_field);
      const auto records = gate_op.ApplyConductorConsistencyGate(E, 0, 1.0, false);
      const auto twin_records = candidate.ApplyConductorConsistencyGate(E, 0, 1.0, false);
      REQUIRE(records.size() == twin_records.size());
      for (std::size_t i = 0; i < records.size(); i++)
      {
        CHECK(records[i].max_ratio == twin_records[i].max_ratio);
        CHECK(records[i].excluded == twin_records[i].excluded);
      }
      CHECK(gate_op.GetLawReadings().empty());
      const auto &readings = candidate.GetLawReadings();
      CHECK(static_cast<int>(readings.size()) <= candidate.GetPatchCount());
      CHECK(!readings.empty());
      // MINOR-8 (decision 2063): LAW-UNPROBED on real geometry: the island's four corner
      // coupons carry no probe (no conductor vertex on their trace mesh): their rows are
      // written with no leg read and no L; every isolated-edge cell reads an on-metal LEG.
      int unprobed = 0, isolated_legs = 0;
      for (const auto &reading : readings)
      {
        const std::string &topology = candidate.ModelTopology(reading.model);
        INFO(topology);
        if (reading.leg_class == SurfaceResponseOperator::LawLegClass::UNPROBED)
        {
          unprobed++;
          CHECK(topology == "convex corner");
          CHECK(reading.legs == 0);
          CHECK(std::isnan(reading.ratio));
        }
        else if (topology == "isolated edge")
        {
          CHECK(reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG);
          isolated_legs++;
        }
      }
      CHECK(unprobed == 4);
      CHECK(isolated_legs > 0);
      // t5: a field normal to the sheet reads zero on every on-metal leg; a tangential
      // field reads its drop. A FRESH twin per field: the gate's exclusion is sticky, so a
      // patch the first field excluded would not be read by a later one.
      auto Readings = [&](double ex, double ey, double ez)
      {
        SurfaceResponseOperator fresh(*twin->iodata, *twin->space, nullptr, &film, true);
        const auto field = Field(*twin->space, ex, ey, ez);
        fresh.ApplyConductorConsistencyGate(TrueField(field), 0, 1.0, false);
        return fresh.GetLawReadings();
      };
      int legs = 0;
      for (const auto &reading : Readings(0.0, 1.0, 0.0))
      {
        if (reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG)
        {
          legs++;
          CHECK(reading.legs > 0);
          CHECK(reading.max_leg <= 1.0e-12);
        }
      }
      CHECK(legs > 0);
      double largest = 0.0;
      int tangential_legs = 0;
      for (const auto &reading : Readings(1.0, 0.0, 1.0))
      {
        if (reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG)
        {
          tangential_legs++;
          largest = std::max(largest, reading.max_leg);
        }
      }
      CHECK(tangential_legs > 0);
      CHECK(largest > 1.0e-3);
    }
  }

  SECTION("the candidate's conductor-consistency Records spill to its own scratch file, "
          "never to memory (MINOR-3)")
  {
    const auto film = Film(*twin);
    SurfaceResponseOperator candidate(*twin->iodata, *twin->space, nullptr, &film, true);
    candidate.SetPatchRecordResidency(
        SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY);
    const fs::path directory = temp.temp_dir / "candidate-records-spill";
    if (Mpi::Root(Mpi::World()))
    {
      fs::create_directories(directory);
    }
    Mpi::Barrier(Mpi::World());
    const std::string name =
        ".palace-impedance-candidate-conductor-consistency-records.spill";
    candidate.SpillConductorConsistencyRecords(directory.string(), name);
    REQUIRE(candidate.HasConductorConsistencyProbes());
    const auto records =
        candidate.ApplyConductorConsistencyGate(TrueField(twin_field), 0, 1.0, false);
    if (Mpi::Root(Mpi::World()))
    {
      const auto *spill = candidate.GetConductorConsistencyRecordsSpill();
      REQUIRE(spill);
      CHECK(fs::path(spill->Path()) == directory / name);
      CHECK(fs::exists(spill->Path()));
      CHECK(spill->Size() == records.size());
      CHECK(!records.empty());
      const auto statistics = candidate.GetStatistics();
      CHECK(statistics["Diagnostics"]["ConductorConsistency"]["Records"].empty());
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

  SECTION("RECORD through the post-operator of a three-dimensional eigenmode run: the "
          "candidate is set up and the two-dimensional line is not printed")
  {
    auto run = Build(Config({{9, film_law}}, sa_library_path), MakeIslandMesh());
    MPI_Comm comm = Mpi::World();
    MakeOutputFolder(*run->iodata, comm);
    const std::string log = CaptureStdout(
        temp.temp_dir / "impedance-record-post.log",
        [&] { PostOperator<ProblemType::EIGENMODE> post_op(*run->iodata, *run->space); });
    CHECK_THAT(log, ContainsSubstring("Impedance-metal rule (RECORD): film = the library's "
                                      "film: an ASSUMPTION"));
    CHECK_THAT(log, !ContainsSubstring("two-dimensional boundary-mode runs"));
    CHECK_THAT(log, !ContainsSubstring("Impedance-metal rule (APPLY)"));
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
    if (candidate.HasConductorConsistencyProbes())
    {
      // MINOR-8 (decision 2063): LAW-NO-LEG on real geometry. An EXACT feature's placement
      // reads its device law (not a perfect conductor): its Maxwell anchor is the edge's
      // own foot, so the probe chain holds no on-metal leg and the row reads LAW-NO-LEG
      // (legs 0, no L); the SUBSTITUTED cells (served as PEC) keep their LEG.
      SurfaceResponseOperator fresh(*run->iodata, *run->space, nullptr, &film, true);
      fresh.ApplyConductorConsistencyGate(TrueField(Field(*run->space, 0.0, 1.0, 0.0)), 0,
                                          1.0, false);
      int no_leg = 0;
      for (const auto &reading : fresh.GetLawReadings())
      {
        const int status = fresh.ImpedanceLawStatusOfFeature(reading.feature);
        if (status == static_cast<int>(ImpedanceLawStatus::EXACT))
        {
          CHECK(reading.leg_class == SurfaceResponseOperator::LawLegClass::NO_LEG);
          CHECK(reading.legs == 0);
          CHECK(std::isnan(reading.ratio));
          no_leg++;
        }
        else if (fresh.ModelTopology(reading.model) == "isolated edge")
        {
          CHECK(reading.leg_class == SurfaceResponseOperator::LawLegClass::LEG);
        }
      }
      CHECK(no_leg > 0);
    }
  }

  // The island split MID-EDGE at x = 0.5 (attribute 9 for x < 0.5, 10 beyond): the two
  // sides along x (z = 0.25, 0.75) change law in their middle.
  auto MidEdgeSplitIslandMesh = [&]()
  {
    mfem::Mesh serial = MakeSerialIslandMesh();
    for (int be = 0; be < serial.GetNBE(); be++)
    {
      if (serial.GetBdrAttribute(be) != 9)
      {
        continue;
      }
      mfem::Array<int> vertices;
      serial.GetBdrElementVertices(be, vertices);
      double x = 0.0;
      for (const int vertex : vertices)
      {
        x += serial.GetVertex(vertex)[0];
      }
      if (x > 0.5 * static_cast<double>(vertices.Size()) * (1.0 + 1.0e-9))
      {
        serial.SetBdrAttribute(be, 10);
      }
    }
    serial.SetAttributes();
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };

  SECTION("a MID-EDGE law change (MINOR-8; MINOR-2 named): the spanning sides LAW-JOINT, "
          "the corners too (the check reaches one neighbouring segment: fail-closed)")
  {
    auto run = Build(Config({{9, "PEC"}, {10, film_law}}, sa_library_path),
                     MidEdgeSplitIslandMesh());
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.ImpedanceRuleActs());
    const auto statuses = Statuses(candidate);
    REQUIRE(statuses.count("IsolatedEdge"));
    const auto &sides = statuses.at("IsolatedEdge");
    INFO("IsolatedEdge LAW-JOINT "
         << (sides.count("LAW-JOINT") ? sides.at("LAW-JOINT") : 0) << ", PEC "
         << (sides.count("PEC") ? sides.at("PEC") : 0) << ", SUBSTITUTED "
         << (sides.count("SUBSTITUTED") ? sides.at("SUBSTITUTED") : 0));
    REQUIRE(sides.count("LAW-JOINT"));
    CHECK(sides.at("LAW-JOINT") >= 2);  // the two sides the law change splits
    // MINOR-2 (not changed; DESIGN-NOTE): every corner is next to a split side's segment
    // of the other law and is refused although it lies wholly in one law.
    REQUIRE(statuses.count("ConvexCorner"));
    const auto &corners = statuses.at("ConvexCorner");
    CHECK(corners.size() == 1);
    CHECK(corners.count("LAW-JOINT") == 1);
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      if (feature.status == "LAW-JOINT")
      {
        CHECK(feature.model.empty());
        CHECK(feature.device_hash == "MIXED");
      }
      else
      {
        // A served or plain feature carries one device law.
        CHECK(feature.device_law.find('|') == std::string::npos);
      }
    }
  }

  SECTION("two different inductances side by side (MINOR-8): MIXED-LAW, served")
  {
    const json second_law = {{"Ls", 2.0e-13}};
    auto run = Build(Config({{9, film_law}, {10, second_law}}, sa_library_path),
                     MakeIslandMesh(false, false, false, true));
    const auto film = Film(*run);
    REQUIRE(film.candidate_differs);
    CHECK(film.attributes.at(9).service == FilmService::FILM_INDUCTIVE);
    CHECK(film.attributes.at(10).service == FilmService::FILM_INDUCTIVE);
    SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
    CHECK(candidate.ImpedanceRuleActs());
    // The two islands within reach form ONE feature (the identification's cluster over both
    // conductors): its two device laws lie on different conductors, so it is MIXED-LAW (a
    // served status), never LAW-JOINT.
    int mixed = 0;
    for (const auto &feature : candidate.GetImpedanceLawFeatures())
    {
      INFO(feature.type << " " << feature.status << " " << feature.device_law);
      CHECK(feature.status != "LAW-JOINT");
      CHECK(feature.status != "NOT-DECLARED");
      if (feature.status == "MIXED-LAW")
      {
        mixed++;
        CHECK(feature.device_hash == "MIXED");
        // Its device laws, both purely inductive sheets, named.
        CHECK_THAT(feature.device_law, ContainsSubstring("|"));
      }
    }
    CHECK(mixed > 0);
  }

  SECTION("the identification-only candidate (MINOR-1): decided on the BUILT statuses; a "
          "perfect conductor and a non-film law on one conductor write the candidate "
          "manifest with its LAW-JOINT refusals; a perfect conductor alone writes none")
  {
    auto Preflight = [&](const Run &run, const std::string &name)
    {
      const fs::path directory = temp.temp_dir / name;
      if (Mpi::Root(Mpi::World()))
      {
        fs::create_directories(directory);
      }
      Mpi::Barrier(Mpi::World());
      const auto path = directory / "surface-response-requirements.json";
      const auto outcome =
          WriteSurfaceResponseRequirements(*run.iodata, *run.meshes[0], path.string());
      Mpi::Barrier(Mpi::World());
      return std::make_pair(
          outcome, directory / "surface-response-impedance-candidate-requirements.json");
    };
    {
      // C7c of review r2038: PEC | Rs 1 split mid-edge along one conductor.
      auto run = Build(Config({{9, "PEC"}, {10, json{{"Rs", 1.0}}}}, sa_library_path),
                       MidEdgeSplitIslandMesh());
      const auto film = Film(*run);
      CHECK_FALSE(film.candidate_differs);  // the film rule alone substitutes nothing
      SurfaceResponseOperator candidate(*run->iodata, *run->space, nullptr, &film, true);
      CHECK(candidate.ImpedanceRuleActs());  // ... but the built candidate refuses
      const auto [outcome, candidate_path] = Preflight(*run, "c7c-preflight");
      REQUIRE(fs::exists(candidate_path));
      CHECK(std::find(outcome.written.begin(), outcome.written.end(),
                      candidate_path.string()) != outcome.written.end());
      const json manifest = json::parse(std::ifstream(candidate_path));
      CHECK(manifest["Candidate"] == "CANDIDATE: NOT OF RECORD");
      const auto &block = manifest["ImpedanceRuleCandidate"];
      CHECK(block["RuleActs"]["Substituted"].get<int>() == 0);
      CHECK(block["RuleActs"]["Refused"].get<int>() > 0);
      CHECK(block["LawStatusCounts"].contains("LAW-JOINT"));
    }
    {
      auto run = Build(Config({{9, "PEC"}}, sa_library_path), MakeIslandMesh());
      const auto [outcome, candidate_path] = Preflight(*run, "pec-preflight");
      CHECK_FALSE(fs::exists(candidate_path));
      CHECK(outcome.written.size() == 2);  // the manifest and the patch dry run
    }
  }
#endif
}

// Ruling 2024 (4): every runtime model topology name the code emits decides the leg class;
// an unknown name fails closed as CONDITIONAL-UNTESTED (never a tested class).
TEST_CASE("Impedance-metal rule law reading leg class of every runtime topology name",
          "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  using LawLegClass = SurfaceResponseOperator::LawLegClass;
  const auto names = SurfaceResponseOperator::LawReadingTopologyNames();
  CHECK(names.size() == 16);
  const std::set<std::string> translational = {"isolated edge",
                                               "same-conductor gap",
                                               "different-conductor gap",
                                               "same-conductor strip",
                                               "parallel-edge cluster",
                                               "curved edge",
                                               "curved same-conductor gap",
                                               "curved different-conductor gap",
                                               "curved same-conductor strip"};
  const std::set<std::string> derived = {"convex corner", "ConvexCorner"};
  const std::set<std::string> conditional = {
      "concave corner", "ConcaveCorner", "spatial edge cluster", "endpoint", "junction"};
  for (const auto &name : names)
  {
    INFO(name);
    bool known = false;
    const auto with_leg =
        SurfaceResponseOperator::LawLegClassOfTopology(name, true, true, &known);
    CHECK(known);
    const auto without_leg =
        SurfaceResponseOperator::LawLegClassOfTopology(name, true, false, nullptr);
    CHECK(SurfaceResponseOperator::LawLegClassOfTopology(name, false, true, nullptr) ==
          LawLegClass::UNPROBED);
    if (translational.count(name))
    {
      CHECK(with_leg == LawLegClass::LEG);
      CHECK(without_leg == LawLegClass::NO_LEG);
    }
    else if (derived.count(name))
    {
      CHECK(with_leg == LawLegClass::LEG_DERIVED);
      CHECK(without_leg == LawLegClass::LEG_DERIVED);
    }
    else
    {
      CHECK(conditional.count(name) == 1);
      CHECK(with_leg == LawLegClass::CONDITIONAL_UNTESTED);
      CHECK(without_leg == LawLegClass::CONDITIONAL_UNTESTED);
    }
  }
  CHECK(translational.size() + derived.size() + conditional.size() == names.size());
  for (const std::string name : {"Explicit", "IsolatedEdge", "bogus", ""})
  {
    INFO(name);
    bool known = true;
    CHECK(SurfaceResponseOperator::LawLegClassOfTopology(name, true, true, &known) ==
          LawLegClass::CONDITIONAL_UNTESTED);
    CHECK_FALSE(known);
  }
}

// MINOR-4 (decision 2063): the candidate's phases book into its own timer categories; the
// record's categories (the nested phases' own) are not counted while a redirect lives.
TEST_CASE("Impedance-metal rule candidate timers book into their own categories",
          "[surfaceresponsefilmrule][impedancemetalrule][Serial]")
{
  CHECK(Timer::descriptions.size() == static_cast<std::size_t>(Timer::NUM_TIMINGS));
  CHECK(Timer::descriptions[Timer::CONSTRUCT_IMPEDANCE_RULE] ==
        "  Impedance Candidate Setup");
  CHECK(Timer::descriptions[Timer::POSTPRO_IMPEDANCE_RULE] == "  Impedance Candidate");
  const auto &timer = BlockTimer::GlobalTimer();
  const int response = timer.Counts(Timer::POSTPRO_RESPONSE);
  const int energy = timer.Counts(Timer::POSTPRO_RESPONSE_ENERGY);
  const int setup = timer.Counts(Timer::CONSTRUCT_RESPONSE);
  const int candidate = timer.Counts(Timer::POSTPRO_IMPEDANCE_RULE);
  const int candidate_setup = timer.Counts(Timer::CONSTRUCT_IMPEDANCE_RULE);
  {
    RedirectingBlockTimer redirect(Timer::POSTPRO_IMPEDANCE_RULE);
    BlockTimer outer(Timer::POSTPRO_RESPONSE);
    {
      BlockTimer inner(Timer::POSTPRO_RESPONSE_ENERGY);
    }
  }
  {
    RedirectingBlockTimer redirect(Timer::CONSTRUCT_IMPEDANCE_RULE);
    BlockTimer inner(Timer::CONSTRUCT_RESPONSE);
  }
  CHECK(timer.Counts(Timer::POSTPRO_RESPONSE) == response);
  CHECK(timer.Counts(Timer::POSTPRO_RESPONSE_ENERGY) == energy);
  CHECK(timer.Counts(Timer::CONSTRUCT_RESPONSE) == setup);
  CHECK(timer.Counts(Timer::POSTPRO_IMPEDANCE_RULE) == candidate + 3);
  CHECK(timer.Counts(Timer::CONSTRUCT_IMPEDANCE_RULE) == candidate_setup + 2);
  {
    // The redirect ends with its scope: a phase books into its own category again.
    BlockTimer after(Timer::POSTPRO_RESPONSE);
  }
  CHECK(timer.Counts(Timer::POSTPRO_RESPONSE) == response + 1);
}

namespace
{

// CUT-PLAN step 3 (lane v1lsdefault): the configuration of a designer's eigenmode run of
// the island whose metal (attribute 9) is an uncracked Impedance sheet of `ls` H/sq, with
// no ImpedanceMetalRule key (the Maxwell default "Apply"), no developer switch and the
// library of the convex island on its SA interface (mesh and library written by the root).
json DefaultSheetEigenmodeConfig(test::SurfaceResponseFiles &files, const std::string &stem,
                                 double ls, const std::string &unmatched_policy)
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path = files.temp.temp_dir / "impedance-default-island.mesh";
  const auto library_path = files.temp.temp_dir / "impedance-default-convex-sa-3d.json";
  if (Mpi::Root(comm))
  {
    mfem::Mesh serial = test::SurfaceResponseFiles::MakeSerialIslandMesh();
    std::ofstream mesh_output(mesh_path);
    serial.Print(mesh_output);
    json library = json::parse(std::ifstream(files.convex_library_3d_path));
    library["Name"] = "unit-test-process-impedance-default-convex-sa-3d";
    for (auto &model : library["Models"])
    {
      model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
    }
    std::ofstream(library_path) << library.dump(2) << "\n";
  }
  Mpi::Barrier(comm);
  json config = files.ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
  config["Model"]["Mesh"] = mesh_path.string();
  config["Model"]["L0"] = 1.0;
  config["Model"]["Lc"] = 1.0;
  config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
  config["Boundaries"]["Impedance"] = json::array({{{"Attributes", {9}}, {"Ls", ls}}});
  config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9};
  auto &correction = config["Solver"]["SurfaceResponseCorrection"];
  correction.erase("PatchConstruction");
  correction.erase("Developer");
  correction["Library"] = library_path.string();
  correction["UnmatchedPolicy"] = unmatched_policy;
  correction["GeometryFailurePolicy"] = "RecordRaw";
  config["Solver"]["Eigenmode"] = {{"Target", 0.05}, {"N", 2}};
  config["Problem"]["Output"] =
      (files.temp.temp_dir / ("impedance-default-" + stem)).string();
  return config;
}

// The run of DefaultSheetEigenmodeConfig through the driver; the output directory.
std::filesystem::path RunDefaultSheetEigenmode(test::SurfaceResponseFiles &files,
                                               const std::string &stem, double ls,
                                               const std::string &unmatched_policy)
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(DefaultSheetEigenmodeConfig(files, stem, ls, unmatched_policy), false);
  REQUIRE(iodata.solver.surface_response_correction->impedance_metal_rule ==
          ResponseCorrectionData::ImpedanceMetalRule::APPLY);
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  palace::Run(iodata, comm, utils::ConfigureOmp(), /*git_tag=*/nullptr);
  Mpi::Barrier(comm);
  return iodata.problem.output;
}

// The words under `header` of the rows of surface-participation.csv that carry a value
// (every row but the NOT-TARGET ones, by their flags).
std::vector<std::string> ValuedRowWords(const std::filesystem::path &dir,
                                        const std::string &header = "flags")
{
  const Table rows = test::LoadCsv(dir / kSurfaceParticipationFile);
  const auto &flags = test::ColumnByHeader(rows, "flags").text;
  // A column whose cells are all empty loads as a number column (no text): empty words.
  const auto &words = test::ColumnByHeader(rows, header).text;
  REQUIRE((words.empty() || words.size() == flags.size()));
  std::vector<std::string> valued;
  for (std::size_t i = 0; i < flags.size(); i++)
  {
    if (flags[i] != "NOT-TARGET")
    {
      valued.push_back(words.empty() ? std::string() : words[i]);
    }
  }
  return valued;
}

}  // namespace

// CUT-PLAN step 3 (lane v1lsdefault; the supervisor's rulings (b) and (c)): with no rule
// key an Ls sheet at the designers' 4.5e-14 H/sq is served: every valued row of the
// designer table carries the served word and palace.json the regime line (the log prints
// the same lines). Above the bound (100 x) with UnmatchedPolicy "Warn" the run goes on raw:
// nothing matches (RawFallback) and every valued row carries the not-served word.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule default in an eigenmode run: the designer table's "
                 "words for a served and a not-served sheet",
                 "[surfaceresponsefilmrule][impedancemetalrule][surfaceparticipation]"
                 "[eigensolver][3d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  fs::path served;
  const std::string served_log = CaptureStdout(
      temp.temp_dir / fmt::format("served-{}.log", Mpi::Rank(Mpi::World())),
      [&] { served = RunDefaultSheetEigenmode(*this, "served", 4.5e-14, "Warn"); });
  const auto above = RunDefaultSheetEigenmode(*this, "above", 4.5e-12, "Warn");
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto ReadMetadata = [](const fs::path &dir)
  { return json::parse(std::ifstream(dir / "palace.json")); };
  {
    // Lane v1wordscode: the served word is an information word; flags hold bad words only.
    // The Maxwell gates of the value reach the rows and the log: this coarse test device's
    // kR is far above 0.1, so every served row names it and reads CHECK; the trace-closure
    // gates are one information line of the run.
    const auto flags = ValuedRowWords(served);
    const auto information = ValuedRowWords(served, "information");
    const auto summary = ValuedRowWords(served, "summary");
    REQUIRE(!flags.empty());
    for (std::size_t i = 0; i < flags.size(); i++)
    {
      CHECK_THAT(information[i], ContainsSubstring(kParticipationSheetServedFlag));
      CHECK_THAT(flags[i], !ContainsSubstring(kParticipationSheetServedFlag));
      CHECK_THAT(flags[i], !ContainsSubstring(kParticipationSheetNotServedFlag));
      CHECK_THAT(flags[i], ContainsSubstring("kR ") && ContainsSubstring(" above 0.1"));
      CHECK(summary[i] == kParticipationSummaryCheck);
    }
    CHECK_THAT(served_log, ContainsSubstring("Maxwell surface response: a gate of the "
                                             "postprocessing-only value failed: kR "));
    const std::string line = "the gap between two closures, information, not a verdict";
    CHECK(served_log.find(line) != std::string::npos);
    CHECK(served_log.find(line) == served_log.rfind(line));
    const auto record = ReadMetadata(served).at("SurfaceResponse").at("ImpedanceMetalRule");
    CHECK(record.at("ServedImpedanceFeatures").get<int>() > 0);
    const auto lines = record.at("Lines").get<std::vector<std::string>>();
    CHECK(std::count_if(lines.begin(), lines.end(),
                        [](const std::string &line)
                        {
                          return line.find("Impedance-metal rule (APPLY): regime: purely "
                                           "inductive sheets at configured Ls <= 9e-14 "
                                           "H/sq are served by the perfect-conductor "
                                           "models; validated") != std::string::npos;
                        }) == 1);
  }
  {
    const auto flags = ValuedRowWords(above);
    const auto information = ValuedRowWords(above, "information");
    const auto summary = ValuedRowWords(above, "summary");
    REQUIRE(!flags.empty());
    for (std::size_t i = 0; i < flags.size(); i++)
    {
      CHECK_THAT(flags[i], ContainsSubstring(kParticipationSheetNotServedFlag));
      CHECK_THAT(information[i], !ContainsSubstring(kParticipationSheetServedFlag));
      CHECK(summary[i] == kParticipationSummaryDoNotUse);  // raw where a correction is due
    }
    const Table rows = test::LoadCsv(above / kSurfaceParticipationFile);
    for (const auto &reason : test::ColumnByHeader(rows, "reason").text)
    {
      CHECK_THAT(reason, ContainsSubstring("RawFallback"));
    }
    const auto record = ReadMetadata(above).at("SurfaceResponse").at("ImpedanceMetalRule");
    CHECK(record.at("ServedImpedanceFeatures").get<int>() == 0);
    const auto lines = record.at("Lines").get<std::vector<std::string>>();
    CHECK(std::count_if(lines.begin(), lines.end(),
                        [](const std::string &line)
                        {
                          return line.find("ABOVE the Ls range <= 9e-14 H/sq): attributes "
                                           "[9]: not served by the perfect-conductor "
                                           "models: raw unless the library holds a model "
                                           "of this law") != std::string::npos;
                        }) == 1);
  }
#endif
}

// The same sheet above the bound with UnmatchedPolicy "Error": the run stops naming its
// features (as on main for every Missing feature; the warning before it names the law).
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule default: a sheet above the bound stops a run under "
                 "UnmatchedPolicy Error",
                 "[surfaceresponsefilmrule][impedancemetalrule][eigensolver][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  std::string message;
  const std::string log = CaptureStdout(
      temp.temp_dir / "impedance-default-error.log",
      [&]
      {
        IoData iodata(DefaultSheetEigenmodeConfig(*this, "error", 4.5e-12, "Error"), false);
        std::vector<std::unique_ptr<Mesh>> meshes;
        meshes.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
        SpaceOperator space_op(iodata, meshes);
        try
        {
          SurfaceResponseOperator response(iodata, space_op);
        }
        catch (const std::exception &error)
        {
          message = error.what();
        }
      });
  CHECK_THAT(message, ContainsSubstring("identified feature(s) have no library model"));
  CHECK_THAT(message, ContainsSubstring("(UnmatchedPolicy = Error)"));
  CHECK_THAT(log, ContainsSubstring("ABOVE the Ls range <= 9e-14 H/sq): attributes [9]"));
  CHECK_THAT(log, ContainsSubstring("under another metal boundary law: device metal"));
#endif
}

// The supervisor's ruling on lane v1lsdefault's finding: a sheet above the bound is not
// served by the perfect-conductor models; its values are raw UNLESS the library holds a
// model of its own law (law status EXACT: corrected by that model). The not-served line
// names both outcomes (it precedes the matching); the feature statuses and the log's
// law-status line say which one the run met, and only the raw one gives the designer
// table's not-served word. Ls 1e-13 H/sq (> 9e-14), no rule key, the Features construction:
// with the perfect-conductor library (no model of the law) every feature is NOT-DECLARED
// and nothing is placed; with the library of that law (the finite-impedance convex
// library, the six value tests' library) every feature is EXACT and corrected by a model.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Impedance-metal rule default: a sheet above the bound is raw unless the "
                 "library holds a model of its law",
                 "[surfaceresponsefilmrule][impedancemetalrule][maxwell][3d][Serial]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  MPI_Comm comm = Mpi::World();
  // Both libraries keyed on the island's one target interface (SA), as the Features
  // construction reads them.
  auto SaLibrary = [&](const fs::path &source, const std::string &name)
  {
    const auto path = temp.temp_dir / (name + ".json");
    if (Mpi::Root(comm))
    {
      json library = json::parse(std::ifstream(source));
      library["Name"] = name;
      for (auto &model : library["Models"])
      {
        model["Interfaces"] = {{{"Type", "SA"}, {"Coupon", 1}}};
      }
      std::ofstream(path) << library.dump(2) << "\n";
    }
    Mpi::Barrier(comm);
    return path;
  };
  const auto pec_library =
      SaLibrary(convex_library_3d_path, "unit-test-process-above-bound-pec-sa-3d");
  const auto law_library = SaLibrary(finite_impedance_convex_library_3d_path,
                                     "unit-test-process-above-bound-impedance-sa-3d");
  const std::string line_text =
      "not served (Impedance Rs 0 Ls 1e-13 Cs 0, cracked (two faces of Ls each): ABOVE the "
      "Ls range <= 9e-14 H/sq): attributes [9]: not served by the perfect-conductor "
      "models: "
      "raw unless the library holds a model of this law (law status EXACT, corrected by "
      "that model); the perfect-conductor models serve Ls sheets at or below 9e-14 H/sq";
  auto Check = [&](const fs::path &library, const std::string &stem, bool exact)
  {
    INFO(stem);
    json config = ConvexMaxwellIslandConfig(test::CrackPolicy::Default);
    config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
    config["Boundaries"]["Impedance"] =
        json::array({{{"Attributes", {9}}, {"Ls", 1.0e-13}}});
    config["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = {9};
    auto &correction = config["Solver"]["SurfaceResponseCorrection"];
    correction.erase("PatchConstruction");
    correction.erase("Developer");
    correction["Library"] = library.string();
    correction["UnmatchedPolicy"] = "Warn";
    correction["GeometryFailurePolicy"] = "RecordRaw";
    IoData iodata(config, false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(MakeIslandMesh()));
    SpaceOperator space(iodata, meshes);
    std::unique_ptr<SurfaceResponseOperator> response;
    const std::string log = CaptureStdout(
        temp.temp_dir / (stem + ".log"),
        [&] { response = std::make_unique<SurfaceResponseOperator>(iodata, space); });
    const auto *film = response->GetImpedanceApplyRule();
    REQUIRE(film != nullptr);
    CHECK(film->attributes.at(9).service == FilmService::INDUCTIVE_ABOVE_RANGE);
    CHECK(HasLine(*film, line_text));
    CHECK_THAT(log, ContainsSubstring(line_text));
    CHECK_THAT(log, ContainsSubstring(exact ? "feature law statuses: EXACT"
                                            : "feature law statuses: NOT-DECLARED"));
    const auto &features = response->GetImpedanceLawFeatures();
    REQUIRE(!features.empty());
    for (const auto &feature : features)
    {
      CHECK(feature.status == (exact ? "EXACT" : "NOT-DECLARED"));
      CHECK(feature.model.empty() == !exact);
    }
    const auto counts = CountImpedanceApplyFeatures(features);
    CHECK(counts.sheet == static_cast<int>(features.size()));
    CHECK(counts.not_served == (exact ? 0 : counts.sheet));
    CHECK((response->GetPatchCount() > 0) == exact);
  };
  Check(pec_library, "above-bound-pec-library", false);
  Check(law_library, "above-bound-law-library", true);
#endif
}
