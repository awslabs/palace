// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_LIBRARY_CHECK_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_LIBRARY_CHECK_HPP

#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

class IoData;

//
// THE CHECKS OF A CONFIGURATION AGAINST ITS FABRICATION-PROCESS RESPONSE LIBRARY. They run
// once, on the root, at load, before the mesh (CheckResponseLibraryAtLoad,
// surfaceresponseoperator.hpp), except the substrate's, which needs the mesh (the domains
// next to the corrected metal) and runs where the library is matched; both print the rows
// of this file in one message.
//
// THE RULE: a mismatch is tolerated while its ESTIMATED EFFECT on a corrected value stays
// within kLibraryMismatchEffectLimit; inside it one warning lists every mismatch with its
// size and effect, beyond it the run stops with the same list (no override: the developer
// switch does not relax it). The effect of one mismatch is its factor
// (LibraryMismatchFactor) times |config / library - 1|; the effects that act on one target
// interface add (its own layer's permittivity and thickness, and the largest of the
// substrate rows), and the limit applies to that sum. |config / library - 1| <=
// kLibraryMismatchEqualTolerance is equality. The other rows stop at any difference (R, a
// missing InterfaceLayers entry, the library reader's own refusals, the length unit, a
// parked metal law without the developer switch, the target set) or warn that nothing can
// be verified (no InterfaceLayers). An effect is beyond the limit only when it exceeds it
// by more than kLibraryMismatchLimitTolerance (relative): a slip of exactly the limit,
// written in decimal, is tolerated. THE VACUUM SIDE: Palace takes the vacuum side of MA and
// SA interfaces only where the material is vacuum by its own test (IsVacuumMaterial); a
// material that is not, and has a permittivity component nearer to 1 than to the library's
// substrate, is neither the substrate nor vacuum (never judged as a substrate slip): named
// at load, it stops the run after the mesh if it touches the corrected metal; with an MA or
// SA target and no vacuum material at all, the run stops at load.
//
inline constexpr double kLibraryMismatchEffectLimit = 0.01;
inline constexpr double kLibraryMismatchEqualTolerance = 1.0e-10;
inline constexpr double kLibraryMismatchLimitTolerance = 1.0e-9;

// Palace's own vacuum test of a configured material (MaterialOperator::IsVacuum, which
// InterfaceDielectricCoefficient uses for the vacuum side of MA and SA interfaces), on the
// full 3 x 3 tensors (the mesh's dimension is not known at load).
bool IsVacuumMaterial(const config::MaterialData &material);

enum class LibraryMismatchQuantity : char
{
  LAYER_PERMITTIVITY,
  LAYER_THICKNESS,
  SUBSTRATE_PERMITTIVITY
};

// The effect on a corrected value of Type `affected` per unit relative mismatch of
// `quantity` (of the layer of Type `layer`; any `layer` for the substrate). A layer's
// mismatch acts on its own Type only (0 otherwise). Where the numbers come from: the .cpp.
double LibraryMismatchFactor(LibraryMismatchQuantity quantity, InterfaceDielectric layer,
                             InterfaceDielectric affected);

// What the check compares, as the run's own library reader reads it.
struct LibraryCheckLibrary
{
  std::string name, path;
  // The reader's refusal (the file, its JSON, its Version, an entry): nothing else is read.
  std::string refusal;
  double matching_radius = 0.0;  // mesh length units
  struct Layer
  {
    double thickness = 0.0;  // mesh length units
    double permittivity = 0.0;
  };
  std::map<InterfaceDielectric, Layer> interface_layers;
  std::optional<double> substrate_permittivity;
  std::optional<std::string> length_unit;  // Fabrication.LengthUnit
  // The fabricated profile of record (decision 2388): key and value, in the library's
  // length unit and degrees (absent: std::nullopt).
  std::vector<std::pair<std::string, std::optional<double>>> profile;
  // The models of a parked metal law (conductivity, rational impedance) the device carries.
  std::vector<std::string> parked_models;
};

// A target interface as the matching selects it (BuildAutomaticResponseData3D): a
// Boundaries.Postprocessing.Dielectric entry of an explicit Type with EdgeDistances,
// listed in TargetInterfaces when that list is given.
struct LibraryCheckTarget
{
  int index = 0;
  InterfaceDielectric type = InterfaceDielectric::DEFAULT;
  std::vector<int> attributes;
  double thickness = 0.0, permittivity = 0.0, radius = 0.0;  // mesh length units
};
std::vector<LibraryCheckTarget>
LibraryCheckTargets(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request);

struct LibraryCheckRow
{
  // EFFECT: judged by the rule; STOP: stops at any difference; WARNING: cannot be verified.
  enum class Kind : char
  {
    EFFECT,
    STOP,
    WARNING
  };
  Kind kind = Kind::STOP;
  std::string quantity, entry, config_value, library_value, tolerance;
  // EFFECT rows: config / library - 1, the two numbers, and the estimated effect on every
  // target interface the mismatch acts on (index -> factor x |mismatch|). `substrate`: a
  // row of the substrate (per target the largest of them counts: one patch lies on one
  // substrate, so a target's error is a weighted mean of the materials' effects).
  double mismatch = std::numeric_limits<double>::quiet_NaN();
  double config_number = std::numeric_limits<double>::quiet_NaN();
  double library_number = std::numeric_limits<double>::quiet_NaN();
  std::map<int, double> effect;
  bool substrate = false;
};

struct LibraryCheckReport
{
  std::vector<LibraryCheckRow> rows;
  // Per target interface: its summed effect and its Type.
  std::map<int, double> target_effect;
  std::map<int, InterfaceDielectric> target_type;
  bool beyond_rule = false;  // a summed effect beyond kLibraryMismatchEffectLimit
  bool stop = false;         // a STOP row, or beyond the rule (no override)
  bool Empty() const { return rows.empty(); }
};

// The rows at load: everything the configuration and the library alone decide.
std::vector<LibraryCheckRow>
FindLibraryLoadRows(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
                    const std::vector<LibraryCheckTarget> &targets,
                    const LibraryCheckLibrary &library);

// The substrate at load, judged WITH the rows found at load (`load_rows`): a configured
// material of the substrate class can pass only if, for every target, the target's layer
// effect plus the substrate factor of its Type times the material's |mismatch| stays within
// the limit. The line names the materials that can pass, those that cannot and those that
// are neither the substrate nor vacuum; the last two stop the run after the mesh if they
// touch the corrected metal (`conditional`: the load verdict is then conditional). When
// none can pass (and the closest is not equal to the library's substrate): the closest
// material's row, which the judge then stops.
struct LibrarySubstrateLine
{
  std::string text;
  std::vector<std::string> conditional;  // e.g. "Domains.Materials[1]"
};
std::optional<LibraryCheckRow>
FindSubstrateLoadRow(const IoData &iodata, const LibraryCheckLibrary &library,
                     const std::vector<LibraryCheckTarget> &targets,
                     const std::vector<LibraryCheckRow> &load_rows,
                     LibrarySubstrateLine *line);

// The substrate after the mesh: one row per material next to the corrected metal that is
// neither the vacuum side nor the library's substrate
// (FindSubstratePermittivityMismatches): an EFFECT row for the substrate class, a STOP row
// for a material that is neither.
std::vector<LibraryCheckRow>
SubstrateMeshRows(const IoData &iodata, const std::map<int, std::vector<int>> &mismatches,
                  double substrate_permittivity,
                  const std::vector<LibraryCheckTarget> &targets);

// The verdict of the rule on the rows: per target the layer rows' effects plus the largest
// substrate row's; beyond the limit the run stops (no override).
LibraryCheckReport JudgeLibraryRows(std::vector<LibraryCheckRow> rows,
                                    const std::vector<LibraryCheckTarget> &targets);

// The one message: the library, the rule, the profile, the substrate line, every row with
// its tolerance and action, the summed effect per Type, and the verdict (conditional on the
// line's `conditional` materials when the run goes on).
std::string LibraryCheckMessage(const LibraryCheckReport &report,
                                const LibraryCheckLibrary &library, const std::string &when,
                                const LibrarySubstrateLine &substrate_line);

// The tolerated mismatches of `report` into the request's record (empty: none).
void RecordLibraryMismatch(
    const LibraryCheckReport &report,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request);

// The record as palace.json writes it (SurfaceResponse.LibraryMismatch: the rule, every
// tolerated mismatch with its size and effect per Type, the summed effect per Type) and
// back (the root's record to every rank).
nlohmann::json LibraryMismatchJson(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record);
config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
LibraryMismatchFromJson(const nlohmann::json &json);

// The record's target interfaces and Types with an effect (LibraryMismatchActsOn's sets,
// surfaceparticipationtable.hpp).
std::set<int> LibraryMismatchInterfaces(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record);
std::set<std::string> LibraryMismatchTypes(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_LIBRARY_CHECK_HPP
