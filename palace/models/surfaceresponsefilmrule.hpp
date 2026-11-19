// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP

#include <map>
#include <string>
#include <vector>
#include "utils/configfile.hpp"

namespace palace
{

class Units;

//
// The impedance-metal rule's film resolution (decisions 2006 / 2012; EXPERIMENTAL, stage A:
// read only under Solver.SurfaceResponseCorrection.ImpedanceMetalRule = "Record"). By
// default, with no tag from the designer, a metal boundary is treated as the library's film
// iff it is a perfect conductor or a PURELY INDUCTIVE Impedance sheet (Rs = 0 and Cs = 0,
// any Ls): an ASSUMPTION, printed on every run. A boundary with Rs != 0 or Cs != 0, a
// Conductivity or a RationalImpedance law is not served. ImpedanceMetalFilmOverrides serves
// or un-serves an attribute in either direction, whatever its law; every override is
// printed.
//

// What the candidate does with the metal of one boundary attribute.
enum class FilmService : char
{
  FILM_PEC,            // a perfect conductor: the library's film, served as today
  FILM_INDUCTIVE,      // a purely inductive sheet: served as the film (substituted)
  NOT_DECLARED,        // another law: not served (as today)
  SERVED_BY_OVERRIDE,  // another law, served by an override (not a purely inductive sheet)
  NOT_SERVED_BY_OVERRIDE  // a film law the designer marked "not the library's film"
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
};

struct FilmRuleResolution
{
  std::map<int, FilmRuleAttribute> attributes;
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
// Impedance, RationalImpedance). `nondimensionalized`: the boundary data hold scaled
// values (printed in SI through `units`). An override naming an attribute that is not a
// metal boundary of the run is refused.
FilmRuleResolution ResolveImpedanceMetalFilm(
    const config::BoundaryData &boundaries,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
    const Units &units, bool nondimensionalized);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_FILM_RULE_HPP
