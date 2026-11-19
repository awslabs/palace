// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP

#include <map>
#include <string>
#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "utils/configfile.hpp"

namespace palace
{

class Units;

//
// The impedance-metal rule's film resolution (decisions 2006 / 2012 / 2206; EXPERIMENTAL:
// read under Solver.SurfaceResponseCorrection.ImpedanceMetalRule "Record" or "Apply"). By
// default, with no tag from the designer, a metal boundary is treated as the library's film
// iff it is a perfect conductor or a PURELY INDUCTIVE Impedance sheet (Rs = 0 and Cs = 0;
// any Ls under "Record"; 0 <= Ls <= ImpedanceMetalMaxInductance under "Apply"): an
// ASSUMPTION, printed on every run. A boundary with Rs != 0 or Cs != 0, a Conductivity or a
// RationalImpedance law is not served. ImpedanceMetalFilmOverrides serves or un-serves an
// attribute: under "Record" in either direction, whatever its law, under "Apply" serving
// only what the rule serves by default; every override is printed.
//

// What the candidate does with the metal of one boundary attribute.
enum class FilmService : char
{
  FILM_PEC,            // a perfect conductor: the library's film, served as today
  FILM_INDUCTIVE,      // a purely inductive sheet: served as the film (substituted)
  NOT_DECLARED,        // another law: not served (as today)
  SERVED_BY_OVERRIDE,  // another law, served by an override (not a purely inductive sheet)
  NOT_SERVED_BY_OVERRIDE,  // a film law the designer marked "not the library's film"
  INDUCTIVE_ABOVE_RANGE    // Apply: a purely inductive sheet above the Ls range: not served
};

// The name of a service in the candidate tables.
const char *FilmServiceName(FilmService service);

struct FilmRuleAttribute
{
  int attribute = 0;
  // The device law in plain words ("PEC", "Impedance Rs 0 Ls 4.5e-14 Cs 0", ...; SI units).
  std::string law;
  FilmService service = FilmService::NOT_DECLARED;
  bool overridden = false;
  // A purely inductive sheet's configured Ls (H/sq, SI; 0 for every other law: the value
  // Apply's range check compares) and the sheet form of the attribute in the mesh
  // ("uncracked (one sheet of total Ls)", "cracked (two faces of Ls each)" or "cracked,
  // legacy split (two faces of 2 Ls each)").
  double inductance = 0.0;
  std::string form;
};

struct FilmRuleResolution
{
  std::map<int, FilmRuleAttribute> attributes;
  // "Apply" (decision 2206 (a)): the film rule is the seam of the operator of record, with
  // the Ls range [0, max_inductance] (H/sq); false = "Record" (the candidate's seam).
  bool apply = false;
  double max_inductance = 0.0;
  bool max_inductance_is_default = true;
  // The Ls values (H/sq, sorted, distinct) of the purely inductive sheets served.
  std::vector<double> served_inductances;
  // The candidate's seam differs from the record's on some attribute (a served law that is
  // not PEC, or an override that changes the default): a candidate operator is built.
  bool candidate_differs = false;
  // The printed lines: the film line ("film = the library's film: an ASSUMPTION: ..."),
  // then one line per override and one per law not served.
  std::vector<std::string> lines;

  // The service of a perimeter segment supported by these metal attributes: served (film
  // or override) only when every attribute is served; an attribute un-served by an
  // override wins; otherwise NOT_DECLARED. An attribute outside the run's metal boundaries
  // is NOT_DECLARED.
  FilmService SegmentService(const std::vector<int> &metal_attributes) const;
};

// The law status of one CANDIDATE feature (DESIGN-IMPK-v3 s5.1 as ruled in 2012 Q2 / Q3):
// formed by the film rule at the seam, classified from the DEVICE laws of the feature's
// segments and of the segments sharing a vertex with them.
enum class ImpedanceLawStatus : signed char
{
  PEC,                     // every device law PEC (served as today)
  EXACT,                   // a verified library model of the device law matches by key
  SUBSTITUTED,             // one film law, served by the PEC model of the served key
  MIXED_LAW,               // different film laws on different conductors: served
  LAW_JOINT,               // one conductor carries >= 2 device laws: Missing, named
  NOT_DECLARED,            // a law outside the film rule: as today
  NOT_SERVED_BY_OVERRIDE,  // the designer marked a metal "not the library's film"
  SERVED_BY_OVERRIDE       // served, and a metal is a non-film law served by override
};

const char *ImpedanceLawStatusName(ImpedanceLawStatus status);

// Whether a service serves the metal (the candidate canonicalises it to the PEC law).
inline bool FilmServiceServes(FilmService service)
{
  return service == FilmService::FILM_PEC || service == FilmService::FILM_INDUCTIVE ||
         service == FilmService::SERVED_BY_OVERRIDE;
}

// The resolution of the rule on the run's metal boundaries (PEC, AuxPEC, Conductivity,
// Impedance, RationalImpedance), in the mode of `request` ("Record" or "Apply").
// `nondimensionalized`: the boundary data hold scaled values (printed and range-checked in
// SI through `units`). An override naming an attribute that is not a metal boundary of the
// run is refused; under "Apply" so is an override that serves what the rule does not.
FilmRuleResolution ResolveImpedanceMetalFilm(
    const config::BoundaryData &boundaries,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
    const Units &units, bool nondimensionalized);

// Whether a feature's law status (ImpedanceLawStatusName) is served on a metal that is not
// a perfect conductor (SUBSTITUTED, MIXED-LAW, SERVED-BY-OVERRIDE).
bool ImpedanceLawStatusServesSheet(const std::string &status);

// The count of an "Apply" run's features on impedance sheets (from the root's feature
// rows): `sheet` (n) the features whose device law is an Impedance sheet, served or not;
// `substituted` (s) the features whose law status serves the sheet; of those, `corrected`
// (k) the features corrected by a perfect-conductor model (matched and placed) and
// `missing` (m = s - k) the features without one (Missing, as for a perfect conductor);
// `corrected_length` the k features' length (the rows' nondimensional length).
struct ImpedanceApplyCounts
{
  int sheet = 0;
  int substituted = 0;
  int corrected = 0;
  int missing = 0;
  double corrected_length = 0.0;
};

ImpedanceApplyCounts CountImpedanceApplyFeatures(
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceLawFeatureData>
        &features);

// The count line of an "Apply" run: "k of n feature(s) on impedance sheets are corrected
// by perfect-conductor models (m Missing, as for a perfect conductor)".
std::string ImpedanceApplyCountText(const ImpedanceApplyCounts &counts);

// The "ImpedanceMetalRule" record of an "Apply" run (palace.json SurfaceResponse and the
// preflight manifest of record; decision 2206 (a)): the mode, the printed lines, the Ls
// range, the Ls values met, the law-status counts and, from the feature rows (`features`:
// the root's), ImpedanceSheetFeatures n, SubstitutedImpedanceFeatures s,
// ServedImpedanceFeatures k (corrected), MissingImpedanceFeatures m and
// ServedImpedanceLength (the k features' length in metres, through `units`); the rows
// themselves when `with_features` (the lines include the one experimental-status line,
// decision 2208 (1)).
nlohmann::json ImpedanceApplyRecord(
    const FilmRuleResolution &film,
    const std::vector<
        config::ElectrostaticSolverData::ResponseCorrectionData::ImpedanceLawFeatureData>
        &features,
    const Units &units, bool with_features);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP
