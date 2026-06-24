// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_DOMAINS_HPP
#define PALACE_SCHEMA_TYPES_DOMAINS_HPP

// Mirrors palace::config::{Domain, Material, DomainEnergy,
// Probe, CurrentDipole, DomainPost} (palace/utils/configfile.hpp).
//
// Deviations from PR 716's domains.json for Phase 1:
//
//   - Permeability / Permittivity / LossTan / Conductivity accept either a
//     scalar or a 3-array. Palace's runtime expands scalars to isotropic
//     material data, so the schema mirrors that wire format with a variant.
//
//   - `CurrentDipole.Direction` uses the `palace::schema::DipoleDirection`
//     alias (Cartesian axis-keyword string OR 3-array) from `common.hpp` —
//     Dirac dipole sources have no notion of a cylindrical frame, so the
//     keyword set excludes `R`/`r` variants.
//
//   - Palace's runtime uses `std::map<int, T>` keyed by Index for Energy /
//     Probe entries. The JSON wire format is always an array with Index as an
//     element member, which is what we model here.

#include <array>
#include <string>
#include <vector>

#include <rfl.hpp>

#include "common.hpp"
#include "schema/utils/annotations.hpp"

namespace palace::schema
{

struct Material
{
  PALACE_SCHEMA_DESC_REQUIRED(Attributes,
                              "Integer array of mesh domain this object applies to.",
                              AttributeList) = {};

  PALACE_SCHEMA_DESC(Permeability,
                     "Relative permeability for this material. Scalar or vector of 3 "
                     "coefficients corresponding to each of `\"MaterialAxes\"`.",
                     ScalarOrVector3) = 1.0;

  PALACE_SCHEMA_DESC(Permittivity,
                     "Relative permittivity for this material. Scalar or vector of 3 "
                     "coefficients corresponding to each of `\"MaterialAxes\"`.",
                     ScalarOrVector3) = 1.0;

  PALACE_SCHEMA_DESC(LossTan,
                     "Loss tangent for this material. Scalar or vector of 3 coefficients "
                     "corresponding to each of `\"MaterialAxes\"`.",
                     ScalarOrVector3) = 0.0;

  PALACE_SCHEMA_DESC(Conductivity,
                     "Electrical conductivity for this material, S/m. Activates the Ohmic "
                     "loss model in this domain. Scalar or vector of 3 coefficients "
                     "corresponding to each of `\"MaterialAxes\"`.",
                     ScalarOrVector3) = 0.0;

  PALACE_SCHEMA_DESC(LondonDepth,
                     "London penetration depth for this material, specified in mesh length "
                     "units. Activates the London equations-based model relating "
                     "superconducting current and electromagnetic fields in this domain.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(
      MaterialAxes,
      "Axes directions for specification of anisotropic material properties. "
      "Required to be unit length and orthogonal.",
      std::array<Vector3, 3>) = std::array<Vector3, 3> {
    {{{1.0, 0.0, 0.0}}, {{0.0, 1.0, 0.0}}, {{0.0, 0.0, 1.0}}}
  };
};

struct CurrentDipole
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index,
      "Index of this current dipole source, used in postprocessing output "
      "files.",
      palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(Moment, "Current dipole moment magnitude, A·m.",
                              double) = 0.0;

  PALACE_SCHEMA_DESC_REQUIRED(
      Center,
      "Coordinates of the dipole center position `[x, y, z]`, in mesh "
      "length units.",
      Vector3) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC_REQUIRED(
      Direction,
      "Direction of the Dirac current source specifying the dipole. "
      "Axis-aligned directions can be specified using keywords: `\"+X\"`, "
      "`\"-X\"`, `\"+Y\"`, `\"-Y\"`, `\"+Z\"`, `\"-Z\"`. The direction can "
      "alternatively be specified as a normalized array of three values, "
      "for example `[0.0, 1.0, 0.0]`.",
      DipoleDirection) = {};
};

struct Energy
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index, "Index of this energy postprocessing domain, used in output files.",
      palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(Attributes,
                              "Integer array of mesh domain this object applies to.",
                              AttributeList) = {};
};

struct Probe
{
  PALACE_SCHEMA_DESC_REQUIRED(Index,
                              "Index of this probe, used in postprocessing output files.",
                              palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(
      Center, "Coordinates of this probe `[x, y, z]`, in mesh length units.", Vector3) = {
    {0.0, 0.0, 0.0}
  };
};

struct DomainPostprocessing
{
  PALACE_SCHEMA_DESC(Energy,
                     "Array of objects for postprocessing domain energies. Postprocesses "
                     "the electric and magnetic field energy inside a given domain. "
                     "Results are written to `domain-E.csv` in the output directory.",
                     std::vector<::palace::schema::Energy>) = {};

  PALACE_SCHEMA_DESC(Probe,
                     "Array of probe points for evaluating field values at specified "
                     "locations in space. The electric field **E** and magnetic flux "
                     "density **B** are probed and written to `probe-E.csv` and "
                     "`probe-B.csv` in the output directory.",
                     std::vector<::palace::schema::Probe>) = {};
};

struct Domains
{
  PALACE_SCHEMA_DESC_REQUIRED(Materials, "Array of material property objects.",
                              std::vector<Material>) = {};

  PALACE_SCHEMA_DESC(CurrentDipole, "Array of current dipole source excitations.",
                     std::vector<::palace::schema::CurrentDipole>) = {};

  PALACE_SCHEMA_DESC(Postprocessing, "Configuration for domain postprocessing.",
                     DomainPostprocessing) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_DOMAINS_HPP
