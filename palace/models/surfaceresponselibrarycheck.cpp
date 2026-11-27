// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceresponselibrarycheck.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mfem.hpp>
#include "linalg/densematrix.hpp"
#include "models/surfaceparticipationtable.hpp"
#include "utils/enum_string.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

std::string Percent(double fraction)
{
  return fmt::format("{:+.3g} %", 100.0 * fraction);
}

std::string Number(double value)
{
  return fmt::format("{:.10g}", value);
}

bool Equal(double config, double library)
{
  return std::abs(config - library) <=
         kLibraryMismatchEqualTolerance * std::max(std::abs(config), std::abs(library));
}

// The one test of the rule's limit (at load and after the mesh alike): an effect of
// exactly the limit (a slip of exactly 1 % written in decimal) is within it.
bool BeyondLimit(double effect)
{
  return effect > kLibraryMismatchEffectLimit * (1.0 + kLibraryMismatchLimitTolerance);
}

mfem::DenseMatrix Dense(const config::SymmetricMatrixData<3> &data)
{
  mfem::DenseMatrix M(3, 3);
  M = 0.0;
  mfem::Vector V(3);
  for (std::size_t i = 0; i < 3; i++)
  {
    for (std::size_t j = 0; j < 3; j++)
    {
      V(j) = data.v[i][j];
    }
    AddMult_a_VVt(data.s[i], V, M);
  }
  return M;
}

}  // namespace

bool IsVacuumMaterial(const config::MaterialData &material)
{
  // MaterialOperator::IsVacuum: the largest light speed, a singular value of (mu
  // eps)^(-1/2), at least 1 - 1e-6, and no frequency-dependent permittivity.
  constexpr double threshold = 1.0 - 1.0e-6;
  mfem::DenseMatrix T(3, 3);
  Mult(Dense(material.mu_r), Dense(material.epsilon_r), T);
  return linalg::SingularValueMax(linalg::MatrixPow(T, -0.5)) >= threshold &&
         !material.HasFrequencyDependentPermittivity();
}

// THE FACTORS. A layer's permittivity or thickness: 1.0 on its own Type (the error of a
// corrected value is the library-served share of it times the mismatch, at most the
// mismatch; independent of a design). The substrate's permittivity changes the field
// itself, so its factor is measured: on two windows of two designs (lane v1loadcheck,
// cluster runs on the container of main fe67940f30): the merge tier's electrostatic
// member on the CTX-002 flat window FC1v2cg3 (Order 3; reads/eps-cl-moves.tsv) and its
// electrostatic member on the CTX-003 window C2 (Order 4, one fixed mesh;
// reads/fix-sub-c2-moves.tsv), the substrate moved by +0.175 % and +1 %, the library's
// models unchanged; effect = |the move of the corrected value - the move of the raw value|,
// its largest over the rows of the Type, over the mismatch (linear in it); per Type the
// larger of the two windows, rounded up to one decimal.
double LibraryMismatchFactor(LibraryMismatchQuantity quantity, InterfaceDielectric layer,
                             InterfaceDielectric affected)
{
  const auto Of = [affected](double ma, double ms, double sa)
  {
    return affected == InterfaceDielectric::MA   ? ma
           : affected == InterfaceDielectric::MS ? ms
           : affected == InterfaceDielectric::SA ? sa
                                                 : 0.0;
  };
  switch (quantity)
  {
    case LibraryMismatchQuantity::LAYER_PERMITTIVITY:
    case LibraryMismatchQuantity::LAYER_THICKNESS:
      return layer == affected ? Of(1.0, 1.0, 1.0) : 0.0;
    case LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY:
      // CTX-002 / CTX-003: MA 0.157 / 0.230, MS 1.858 / 1.795, SA 0.102 / 0.019.
      return Of(0.3, 1.9, 0.2);
  }
  return 0.0;
}

namespace
{

std::string TargetEntry(const LibraryCheckTarget &target)
{
  return fmt::format("Boundaries.Postprocessing.Dielectric Index {} ({}; attributes {})",
                     target.index, ToString(target.type),
                     fmt::join(target.attributes, ", "));
}

// The tolerance of one quantity alone as the rule reads it (the limit over the largest
// factor of that quantity), printed beside each row.
std::string AloneTolerance(LibraryMismatchQuantity quantity, InterfaceDielectric layer)
{
  double factor = 0.0;
  for (const auto affected :
       {InterfaceDielectric::MA, InterfaceDielectric::MS, InterfaceDielectric::SA})
  {
    factor = std::max(factor, LibraryMismatchFactor(quantity, layer, affected));
  }
  return fmt::format("effect = {:.2f} x |mismatch| on {}; alone within {:.3g} %", factor,
                     quantity == LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY
                         ? std::string("the most affected Type")
                         : std::string(ToString(layer)),
                     100.0 * kLibraryMismatchEffectLimit / factor);
}

}  // namespace

std::vector<LibraryCheckTarget>
LibraryCheckTargets(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request)
{
  const double scale = iodata.units.GetMeshLengthRelativeScale();
  const std::set<int> filter(request.target_interfaces.begin(),
                             request.target_interfaces.end());
  std::vector<LibraryCheckTarget> targets;
  for (const auto &[index, dielectric] : iodata.boundaries.postpro.dielectric)
  {
    if ((!filter.empty() && filter.count(index) == 0) ||
        dielectric.type == InterfaceDielectric::DEFAULT ||
        dielectric.edge_distances.empty())
    {
      continue;
    }
    targets.push_back({index, dielectric.type, dielectric.attributes, dielectric.t * scale,
                       dielectric.epsilon_r, dielectric.edge_distances.back() * scale});
  }
  return targets;
}

namespace
{

// The metres of a Fabrication.LengthUnit word (0: not a unit this check knows). The micro
// forms: "um", the micro sign U+00B5, the Greek mu U+03BC, "micron".
double LengthUnitMetres(const std::string &unit)
{
  const std::map<std::string, double> units = {
      {"m", 1.0},          {"cm", 1.0e-2},      {"mm", 1.0e-3},     {"um", 1.0e-6},
      {"\u00b5m", 1.0e-6}, {"\u03bcm", 1.0e-6}, {"micron", 1.0e-6}, {"nm", 1.0e-9}};
  const auto it = units.find(unit);
  return it != units.end() ? it->second : 0.0;
}

LibraryCheckRow StopRow(std::string quantity, std::string entry, std::string config_value,
                        std::string library_value, std::string tolerance = "exact")
{
  LibraryCheckRow row;
  row.kind = LibraryCheckRow::Kind::STOP;
  row.quantity = std::move(quantity);
  row.entry = std::move(entry);
  row.config_value = std::move(config_value);
  row.library_value = std::move(library_value);
  row.tolerance = std::move(tolerance);
  return row;
}

// A row judged by the rule: the relative mismatch and its effect on the listed targets.
LibraryCheckRow EffectRow(LibraryMismatchQuantity quantity, InterfaceDielectric layer,
                          std::string name, std::string entry, double config,
                          double library, const std::string &library_key,
                          const std::vector<LibraryCheckTarget> &affected)
{
  LibraryCheckRow row;
  row.kind = LibraryCheckRow::Kind::EFFECT;
  row.quantity = std::move(name);
  row.entry = std::move(entry);
  row.config_number = config;
  row.library_number = library;
  row.config_value = Number(config);
  row.library_value = Number(library) + " (" + library_key + ")";
  row.mismatch = config / library - 1.0;
  row.tolerance = AloneTolerance(quantity, layer);
  row.substrate = quantity == LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY;
  for (const auto &target : affected)
  {
    row.effect[target.index] =
        LibraryMismatchFactor(quantity, layer, target.type) * std::abs(row.mismatch);
  }
  return row;
}

// The vacuum side at load (below): an MA or SA target with no material vacuum at all.
constexpr const char *kVacuumSideQuantity = "vacuum side (MA and SA targets)";
std::optional<LibraryCheckRow>
VacuumSideRow(const IoData &iodata, const LibraryCheckLibrary &library,
              const std::vector<LibraryCheckTarget> &targets);

}  // namespace

std::vector<LibraryCheckRow>
FindLibraryLoadRows(const IoData &iodata,
                    const config::ElectrostaticSolverData::ResponseCorrectionData &request,
                    const std::vector<LibraryCheckTarget> &targets,
                    const LibraryCheckLibrary &library)
{
  std::vector<LibraryCheckRow> rows;
  const std::string library_entry = "Library \"" + request.library + "\"";
  if (!library.refusal.empty())
  {
    rows.push_back(StopRow("library file", library_entry, "-", library.refusal, "-"));
    return rows;
  }
  // The target set (the matching's own refusals, unchanged in action).
  if (targets.empty())
  {
    rows.push_back(StopRow("target set", "Boundaries.Postprocessing.Dielectric",
                           "no entry of an explicit Type with EdgeDistances",
                           "at least one target interface", "-"));
  }
  std::set<int> missing(request.target_interfaces.begin(), request.target_interfaces.end());
  // (AutomaticEdges / EdgeAttributes depend on the mesh's dimension: checked with it.)
  for (const auto &target : targets)
  {
    missing.erase(target.index);
  }
  for (const int index : missing)
  {
    rows.push_back(StopRow("target set", fmt::format("TargetInterfaces entry {}", index),
                           "missing, untyped, or without EdgeDistances",
                           "a target interface", "-"));
  }
  // R: the library's models replace the energy inside it, so any other R is another model.
  for (const auto &target : targets)
  {
    if (!Equal(target.radius, library.matching_radius))
    {
      rows.push_back(StopRow("R (the largest EdgeDistances)", TargetEntry(target),
                             Number(target.radius) + " mesh length units",
                             Number(library.matching_radius) + " (MatchingRadius)"));
    }
  }
  if (library.length_unit)
  {
    const double metres = LengthUnitMetres(*library.length_unit);
    if (!Equal(iodata.model.L0, metres))
    {
      rows.push_back(StopRow(
          "length unit", "Model.L0", Number(iodata.model.L0) + " m",
          "\"" + *library.length_unit + "\" (Fabrication.LengthUnit)" +
              (metres > 0.0 ? std::string()
                            : std::string(": unknown LengthUnit word (known: m, cm, mm, "
                                          "um, \u00b5m, \u03bcm, micron, nm)"))));
    }
  }
  if (!request.developer)
  {
    for (const auto &model : library.parked_models)
    {
      rows.push_back(StopRow(
          "metal law", "model \"" + model + "\"", "the device carries its law",
          "a conductivity or rational-impedance coupon (parked)", "\"Developer\": true"));
    }
  }
  // The layers (decision 1009's InterfaceLayers): a library without them cannot be checked.
  if (library.interface_layers.empty())
  {
    LibraryCheckRow row = StopRow("InterfaceLayers", library_entry, "-",
                                  "no InterfaceLayers metadata: the target thicknesses "
                                  "and permittivities cannot be verified",
                                  "-");
    row.kind = LibraryCheckRow::Kind::WARNING;
    rows.push_back(std::move(row));
    return rows;
  }
  for (const auto &target : targets)
  {
    const std::string type = ToString(target.type);
    const auto layer = library.interface_layers.find(target.type);
    if (layer == library.interface_layers.end())
    {
      rows.push_back(StopRow(type + " InterfaceLayers entry", TargetEntry(target), "-",
                             "none for " + type, "-"));
      continue;
    }
    if (!Equal(target.thickness, layer->second.thickness))
    {
      rows.push_back(EffectRow(
          LibraryMismatchQuantity::LAYER_THICKNESS, target.type, type + " layer thickness",
          TargetEntry(target), target.thickness, layer->second.thickness,
          "Fabrication.InterfaceLayers." + type + ".Thickness", {target}));
    }
    if (!Equal(target.permittivity, layer->second.permittivity))
    {
      rows.push_back(EffectRow(LibraryMismatchQuantity::LAYER_PERMITTIVITY, target.type,
                               type + " layer permittivity", TargetEntry(target),
                               target.permittivity, layer->second.permittivity,
                               "Fabrication.InterfaceLayers." + type + ".Permittivity",
                               {target}));
    }
  }
  if (auto row = VacuumSideRow(iodata, library, targets))
  {
    rows.push_back(std::move(*row));
  }
  return rows;
}

namespace
{

// A material against the library's substrate (THE VACUUM SIDE, the .hpp): the vacuum side
// (IsVacuumMaterial); NEITHER the substrate nor vacuum (a permittivity component nearer to
// 1 than to the substrate); else the substrate class, its mismatch the component farthest
// from the substrate, as config / library - 1.
struct MaterialClass
{
  enum class Kind : char
  {
    VACUUM,
    SUBSTRATE,
    NEITHER
  };
  Kind kind = Kind::VACUUM;
  double mismatch = 0.0;
};

// The words of a material that is neither (the load line, the rows).
constexpr const char *kNotVacuumWords =
    "not the library's substrate and not vacuum (Palace takes the vacuum side of MA and SA "
    "interfaces only where the material is vacuum: permittivity and permeability exactly "
    "1, no frequency dependence; if it is the air, set its permittivity to exactly 1.0; if "
    "it is a dielectric, use a library made for it)";

MaterialClass ClassifyMaterial(const config::MaterialData &material, double substrate)
{
  const auto &s = material.epsilon_r.s;
  if (IsVacuumMaterial(material))
  {
    return {MaterialClass::Kind::VACUUM};
  }
  if (std::any_of(s.begin(), s.end(),
                  [&](double e) { return std::abs(e - 1.0) < std::abs(e - substrate); }))
  {
    return {MaterialClass::Kind::NEITHER};
  }
  double mismatch = 0.0;
  for (const double e : s)
  {
    if (!Equal(e, substrate) && std::abs(e / substrate - 1.0) > std::abs(mismatch))
    {
      mismatch = e / substrate - 1.0;
    }
  }
  return {MaterialClass::Kind::SUBSTRATE, mismatch};
}

std::string MaterialValue(const std::array<double, 3> &s)
{
  return s[0] == s[1] && s[1] == s[2]
             ? Number(s[0])
             : fmt::format("{} / {} / {}", Number(s[0]), Number(s[1]), Number(s[2]));
}

// With an MA or SA target and no configured material vacuum by Palace's test, no interface
// has a vacuum side (its raw MA and SA read 0): a STOP row naming the materials nearer to 1
// than to the library's substrate (every material when the library states none).
std::optional<LibraryCheckRow> VacuumSideRow(const IoData &iodata,
                                             const LibraryCheckLibrary &library,
                                             const std::vector<LibraryCheckTarget> &targets)
{
  const auto &materials = iodata.domains.materials;
  if (std::none_of(targets.begin(), targets.end(),
                   [](const LibraryCheckTarget &target)
                   {
                     return target.type == InterfaceDielectric::MA ||
                            target.type == InterfaceDielectric::SA;
                   }) ||
      std::any_of(materials.begin(), materials.end(), IsVacuumMaterial))
  {
    return std::nullopt;
  }
  std::vector<std::string> named;
  for (std::size_t m = 0; m < materials.size(); m++)
  {
    if (!library.substrate_permittivity ||
        ClassifyMaterial(materials[m], *library.substrate_permittivity).kind ==
            MaterialClass::Kind::NEITHER)
    {
      named.push_back(fmt::format("Domains.Materials[{}] (attributes {}) {}", m,
                                  fmt::join(materials[m].attributes, ", "),
                                  MaterialValue(materials[m].epsilon_r.s)));
    }
  }
  return StopRow(
      kVacuumSideQuantity, "Domains.Materials: no material is vacuum",
      named.empty() ? std::string("none nearer to 1 than to the substrate")
                    : fmt::format("{}", fmt::join(named, "; ")),
      "a vacuum material",
      "Palace takes the vacuum side of MA and SA interfaces only where the "
      "material is vacuum (permittivity and permeability exactly 1, no frequency "
      "dependence): set the air's permittivity to exactly 1.0");
}

// Per target: the summed effect of the layer rows (every EFFECT row but the substrate's).
std::map<int, double> LayerEffects(const std::vector<LibraryCheckRow> &rows)
{
  std::map<int, double> effects;
  for (const auto &row : rows)
  {
    if (row.substrate)
    {
      continue;
    }
    for (const auto &[index, effect] : row.effect)
    {
      effects[index] += effect;
    }
  }
  return effects;
}

// A substrate material of relative mismatch `mismatch` passes the rule beside the layer
// rows: for every target, its layer effect + the substrate factor of its Type x |mismatch|.
bool SubstrateCanPass(double mismatch, const std::vector<LibraryCheckTarget> &targets,
                      const std::map<int, double> &layer_effects)
{
  return std::all_of(
      targets.begin(), targets.end(),
      [&](const LibraryCheckTarget &target)
      {
        const auto layer = layer_effects.find(target.index);
        return !BeyondLimit(
            (layer != layer_effects.end() ? layer->second : 0.0) +
            LibraryMismatchFactor(LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY,
                                  InterfaceDielectric::DEFAULT, target.type) *
                std::abs(mismatch));
      });
}

}  // namespace

std::optional<LibraryCheckRow>
FindSubstrateLoadRow(const IoData &iodata, const LibraryCheckLibrary &library,
                     const std::vector<LibraryCheckTarget> &targets,
                     const std::vector<LibraryCheckRow> &load_rows,
                     LibrarySubstrateLine *line)
{
  *line = {};
  if (!library.substrate_permittivity || targets.empty())
  {
    return std::nullopt;
  }
  const double substrate = *library.substrate_permittivity;
  // A material of the substrate class could be the substrate if, beside the layer rows, it
  // passes the rule; the others (failing, or neither) are named with their condition.
  const auto layer_effects = LayerEffects(load_rows);
  const auto &materials = iodata.domains.materials;
  std::vector<std::string> names, failing, neither;
  std::optional<std::size_t> closest;
  double closest_mismatch = 0.0;
  for (std::size_t m = 0; m < materials.size(); m++)
  {
    const auto &s = materials[m].epsilon_r.s;
    const auto material = ClassifyMaterial(materials[m], substrate);
    const std::string entry = fmt::format("Domains.Materials[{}]", m);
    const std::string named =
        fmt::format("{} (attributes {}) {}", entry,
                    fmt::join(materials[m].attributes, ", "), MaterialValue(s));
    if (material.kind == MaterialClass::Kind::VACUUM)
    {
      continue;
    }
    if (material.kind == MaterialClass::Kind::NEITHER)
    {
      neither.push_back(named);
      line->conditional.push_back(entry);
      continue;
    }
    if (SubstrateCanPass(material.mismatch, targets, layer_effects))
    {
      names.push_back(fmt::format("{} ({})", named, Percent(material.mismatch)));
    }
    else
    {
      failing.push_back(fmt::format("{} ({})", named, Percent(material.mismatch)));
      line->conditional.push_back(entry);
    }
    if (!closest || std::abs(material.mismatch) < std::abs(closest_mismatch))
    {
      closest = m;
      closest_mismatch = material.mismatch;
    }
  }
  const std::string stops = ": it stops the run after the mesh if it touches the corrected "
                            "metal";
  line->text = fmt::format(
      "substrate (checked after the mesh, on the domains next to the corrected metal, with "
      "the layer rows above): library {} (Fabrication.SubstratePermittivity); could be the "
      "substrate: {}",
      Number(substrate),
      names.empty() ? std::string("none") : fmt::format("{}", fmt::join(names, "; ")));
  for (const auto &named : failing)
  {
    line->text += "; " + named + ": beyond the rule with the layer rows" + stops;
  }
  for (const auto &named : neither)
  {
    line->text += fmt::format("; {}: {}{}", named, kNotVacuumWords, stops);
  }
  if (!closest)
  {
    // No material of the substrate class (row (b) of the table at load): with an MS or SA
    // target none can be the library's substrate. Not a substrate slip: no mismatch.
    if (std::none_of(targets.begin(), targets.end(),
                     [](const LibraryCheckTarget &target)
                     {
                       return target.type == InterfaceDielectric::MS ||
                              target.type == InterfaceDielectric::SA;
                     }))
    {
      return std::nullopt;
    }
    return StopRow("substrate (MS and SA targets)",
                   "Domains.Materials: no configured material can be the library's "
                   "substrate",
                   neither.empty() ? std::string("none (vacuum only)")
                                   : fmt::format("{}", fmt::join(neither, "; ")),
                   Number(substrate) + " (Fabrication.SubstratePermittivity)",
                   "a material nearer to the library's substrate than to 1");
  }
  // A material equal to the library's substrate is no row (the layer rows alone are then
  // beyond the rule and stop the run).
  if (!names.empty() || closest_mismatch == 0.0)
  {
    return std::nullopt;
  }
  // No material of the substrate class can pass the check after the mesh: judged now, on
  // the closest (beyond the rule with the layer rows, by construction: the sum grows with
  // |mismatch|).
  const auto &s = materials[*closest].epsilon_r.s;
  auto row = EffectRow(
      LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY, InterfaceDielectric::DEFAULT,
      "substrate permittivity (no material of the configuration within "
      "the rule beside the layer rows: the check after the mesh cannot pass)",
      fmt::format("Domains.Materials[{}] (attributes {}; the closest)", *closest,
                  fmt::join(materials[*closest].attributes, ", ")),
      (1.0 + closest_mismatch) * substrate, substrate, "Fabrication.SubstratePermittivity",
      targets);
  row.config_value = MaterialValue(s);
  return row;
}

std::vector<LibraryCheckRow>
SubstrateMeshRows(const IoData &iodata, const std::map<int, std::vector<int>> &mismatches,
                  double substrate_permittivity,
                  const std::vector<LibraryCheckTarget> &targets)
{
  std::vector<LibraryCheckRow> rows;
  for (const auto &[material, attributes] : mismatches)
  {
    const auto &s = iodata.domains.materials[material].epsilon_r.s;
    const auto side =
        ClassifyMaterial(iodata.domains.materials[material], substrate_permittivity);
    const std::string entry =
        fmt::format("Domains.Materials[{}] (attributes next to the corrected metal: {})",
                    material, fmt::join(attributes, ", "));
    if (side.kind == MaterialClass::Kind::VACUUM)
    {
      continue;
    }
    if (side.kind == MaterialClass::Kind::NEITHER)
    {
      rows.push_back(StopRow(
          "material next to the corrected metal", entry, MaterialValue(s),
          fmt::format("the substrate {} (Fabrication.SubstratePermittivity) or vacuum",
                      Number(substrate_permittivity)),
          kNotVacuumWords));
      continue;
    }
    auto row = EffectRow(
        LibraryMismatchQuantity::SUBSTRATE_PERMITTIVITY, InterfaceDielectric::DEFAULT,
        "substrate permittivity", entry, (1.0 + side.mismatch) * substrate_permittivity,
        substrate_permittivity, "Fabrication.SubstratePermittivity", targets);
    row.config_value = MaterialValue(s);
    rows.push_back(std::move(row));
  }
  return rows;
}

LibraryCheckReport JudgeLibraryRows(std::vector<LibraryCheckRow> rows,
                                    const std::vector<LibraryCheckTarget> &targets)
{
  LibraryCheckReport report;
  for (const auto &target : targets)
  {
    report.target_type[target.index] = target.type;
  }
  // Per target: the layer rows add; of the substrate rows the largest counts.
  report.target_effect = LayerEffects(rows);
  std::map<int, double> substrate_effect;
  for (const auto &row : rows)
  {
    report.stop = report.stop || row.kind == LibraryCheckRow::Kind::STOP;
    for (const auto &[index, effect] : row.effect)
    {
      if (row.substrate)
      {
        substrate_effect[index] = std::max(substrate_effect[index], effect);
      }
    }
  }
  for (const auto &[index, effect] : substrate_effect)
  {
    report.target_effect[index] += effect;
  }
  for (const auto &[index, effect] : report.target_effect)
  {
    report.beyond_rule = report.beyond_rule || BeyondLimit(effect);
  }
  report.stop = report.stop || report.beyond_rule;
  report.rows = std::move(rows);
  return report;
}

namespace
{

std::string RowAction(const LibraryCheckRow &row, const LibraryCheckReport &report)
{
  if (row.kind != LibraryCheckRow::Kind::EFFECT)
  {
    return row.kind == LibraryCheckRow::Kind::STOP ? "STOP" : "WARNING";
  }
  const bool beyond =
      std::any_of(row.effect.begin(), row.effect.end(), [&](const auto &entry)
                  { return BeyondLimit(report.target_effect.at(entry.first)); });
  return beyond ? "STOP (beyond the rule)" : "TOLERATED";
}

std::string RowLine(const LibraryCheckRow &row, const LibraryCheckReport &report)
{
  std::string line =
      fmt::format("    {}: {}; {}; config {}; library {}", RowAction(row, report),
                  row.quantity, row.entry, row.config_value, row.library_value);
  if (row.kind == LibraryCheckRow::Kind::EFFECT)
  {
    double largest = 0.0;
    for (const auto &[index, effect] : row.effect)
    {
      largest = std::max(largest, effect);
    }
    line += fmt::format("; mismatch {}; estimated effect {:.3g} %", Percent(row.mismatch),
                        100.0 * largest);
  }
  return line + "; tolerance: " + row.tolerance + "\n";
}

}  // namespace

std::string LibraryCheckMessage(const LibraryCheckReport &report,
                                const LibraryCheckLibrary &library, const std::string &when,
                                const LibrarySubstrateLine &substrate_line)
{
  std::string message = fmt::format(
      "Response-correction library check {}: library \"{}\" ({})\n  the rule: a mismatch "
      "is tolerated while its estimated effect on a corrected value stays within {:.3g} % "
      "and stops the run beyond it, with no override (effect = the factor x |config / "
      "library - 1|, summed over the layer mismatches acting on one target interface plus "
      "the largest substrate one); R, a missing InterfaceLayers entry, the library "
      "file, the length unit, a parked metal law and the target set stop at any "
      "difference\n",
      when, library.name, library.path, 100.0 * kLibraryMismatchEffectLimit);
  if (!library.profile.empty())
  {
    std::vector<std::string> profile;
    for (const auto &[key, value] : library.profile)
    {
      profile.push_back(key + " " + (value ? Number(*value) : std::string("absent")));
    }
    message +=
        fmt::format("  the library's profile (Fabrication; length unit {}, "
                    "degrees): {}\n",
                    library.length_unit.value_or("not stated"), fmt::join(profile, ", "));
  }
  if (!substrate_line.text.empty())
  {
    message += "  " + substrate_line.text + "\n";
  }
  // A run that goes on is tolerated only conditionally while a configured material would
  // stop it after the mesh if it touched the corrected metal (named in the line above).
  const std::string unless =
      substrate_line.conditional.empty()
          ? std::string()
          : fmt::format(" unless {} touches the corrected metal (then the run stops after "
                        "the mesh)",
                        fmt::join(substrate_line.conditional, " or "));
  if (report.Empty())
  {
    return message + (unless.empty() ? std::string("  no mismatch\n")
                                     : "  no row at load\n  -> TOLERATED" + unless + "\n");
  }
  message += fmt::format("  {} row(s):\n", report.rows.size());
  for (const auto &row : report.rows)
  {
    message += RowLine(row, report);
  }
  // Per Type: the largest summed effect of one of its target interfaces.
  std::vector<std::string> per_type;
  for (const auto type :
       {InterfaceDielectric::MA, InterfaceDielectric::MS, InterfaceDielectric::SA})
  {
    std::optional<std::pair<int, double>> largest;
    for (const auto &[index, effect] : report.target_effect)
    {
      if (report.target_type.at(index) == type && (!largest || effect > largest->second))
      {
        largest = std::make_pair(index, effect);
      }
    }
    if (largest && largest->second > 0.0)
    {
      per_type.push_back(fmt::format("{} {:.3g} % (Index {})", ToString(type),
                                     100.0 * largest->second, largest->first));
    }
  }
  if (!per_type.empty())
  {
    message += fmt::format("  estimated effect per Type (the largest sum over its target "
                           "interfaces; the limit {:.3g} %): {}\n",
                           100.0 * kLibraryMismatchEffectLimit, fmt::join(per_type, ", "));
  }
  const bool effect_rows =
      std::any_of(report.rows.begin(), report.rows.end(), [](const LibraryCheckRow &row)
                  { return row.kind == LibraryCheckRow::Kind::EFFECT; });
  const bool parked_law =
      std::any_of(report.rows.begin(), report.rows.end(),
                  [](const LibraryCheckRow &row) { return row.quantity == "metal law"; });
  // No library can serve air that is not vacuum (its raw MA and SA are 0): that remedy
  // first, the library's for the other rows after it.
  const bool vacuum_side =
      std::any_of(report.rows.begin(), report.rows.end(), [](const LibraryCheckRow &row)
                  { return row.quantity == kVacuumSideQuantity; });
  if (report.stop)
  {
    message += vacuum_side ? "  -> STOP: set the air's permittivity to exactly 1.0" : "";
    if (!vacuum_side || report.rows.size() > 1)
    {
      message +=
          vacuum_side ? "; set each other entry above" : "  -> STOP: set each entry above";
      message += " to the library's value, or use a library made for the configuration's "
                 "values";
    }
    if (parked_law)
    {
      message += ", or set \"Developer\": true for the parked coupons named above";
    }
    message += report.beyond_rule
                   ? " (no override: an effect beyond the rule stops the run with or "
                     "without \"Developer\": true)\n"
                   : "\n";
  }
  else if (effect_rows)
  {
    message +=
        "  -> TOLERATED" + unless +
        ": the run goes on; the outputs name it (the flags word \"" +
        std::string(kParticipationLibraryMismatchFlag) +
        "\" of surface-participation.csv, palace.json SurfaceResponse.LibraryMismatch)\n";
  }
  else
  {
    message += "  -> WARNING: the run goes on unverified" + unless + "\n";
  }
  return message;
}

void RecordLibraryMismatch(
    const LibraryCheckReport &report,
    const config::ElectrostaticSolverData::ResponseCorrectionData &request)
{
  auto &record = request.library_mismatch;
  record = {};
  if (report.stop)
  {
    return;
  }
  for (const auto &row : report.rows)
  {
    if (row.kind != LibraryCheckRow::Kind::EFFECT)
    {
      continue;
    }
    config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord::Entry
        entry;
    entry.quantity = row.quantity;
    entry.entry = row.entry;
    entry.config_value = row.config_number;
    entry.library_value = row.library_number;
    entry.mismatch = row.mismatch;
    for (const auto &[index, effect] : row.effect)
    {
      auto &value = entry.effect[ToString(report.target_type.at(index))];
      value = std::max(value, effect);
    }
    record.entries.push_back(std::move(entry));
  }
  for (const auto &[index, effect] : report.target_effect)
  {
    if (effect > 0.0)
    {
      auto &value = record.effect_per_type[ToString(report.target_type.at(index))];
      value = std::max(value, effect);
      record.interfaces.push_back(index);
    }
  }
}

nlohmann::json LibraryMismatchJson(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  if (record.Empty())
  {
    return nlohmann::json::object();
  }
  nlohmann::json json;
  json["Rule"] = fmt::format(
      "tolerated while the estimated effect on a corrected value stays within {:.3g} % "
      "(effect = the factor x |config / library - 1|, summed per target interface, the "
      "largest substrate mismatch counting); beyond it the run stops (no override)",
      100.0 * kLibraryMismatchEffectLimit);
  json["Mismatches"] = nlohmann::json::array();
  for (const auto &entry : record.entries)
  {
    json["Mismatches"].push_back({{"Quantity", entry.quantity},
                                  {"Entry", entry.entry},
                                  {"Config", entry.config_value},
                                  {"Library", entry.library_value},
                                  {"Mismatch", entry.mismatch},
                                  {"EstimatedEffect", entry.effect}});
  }
  json["EstimatedEffectPerType"] = record.effect_per_type;
  json["Interfaces"] = record.interfaces;
  return json;
}

config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
LibraryMismatchFromJson(const nlohmann::json &json)
{
  config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord record;
  if (!json.contains("Mismatches"))
  {
    return record;
  }
  for (const auto &item : json.at("Mismatches"))
  {
    record.entries.push_back(
        {item.at("Quantity").get<std::string>(), item.at("Entry").get<std::string>(),
         item.at("Config").get<double>(), item.at("Library").get<double>(),
         item.at("Mismatch").get<double>(),
         item.at("EstimatedEffect").get<std::map<std::string, double>>()});
  }
  record.effect_per_type =
      json.at("EstimatedEffectPerType").get<std::map<std::string, double>>();
  record.interfaces = json.at("Interfaces").get<std::vector<int>>();
  return record;
}

std::set<int> LibraryMismatchInterfaces(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  return {record.interfaces.begin(), record.interfaces.end()};
}

std::set<std::string> LibraryMismatchTypes(
    const config::ElectrostaticSolverData::ResponseCorrectionData::LibraryMismatchRecord
        &record)
{
  std::set<std::string> types;
  for (const auto &[type, effect] : record.effect_per_type)
  {
    types.insert(type);
  }
  return types;
}

}  // namespace palace
