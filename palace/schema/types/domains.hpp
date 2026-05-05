// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_DOMAINS_HPP
#define PALACE_SCHEMA_TYPES_DOMAINS_HPP

// Mirrors palace::config::{DomainData, MaterialData, DomainEnergyData,
// ProbeData, CurrentDipoleData, DomainPostData} (palace/utils/configfile.hpp).
//
// Deviations from PR 716's domains.json for Phase 1:
//
//   - Permeability / Permittivity / LossTan / Conductivity: PR 716 accepts
//     scalar-or-3-array via `oneOf`. reflect-cpp cannot emit a `oneOf` of
//     unlike primitive vs array shapes from a single C++ type; we standardize
//     on the 3-array form here. Users migrating from scalar notation must
//     expand `2.0` to `[2.0, 2.0, 2.0]`.
//
//   - MaterialAxes and the CurrentDipole.Direction `anyOf` (string-or-array)
//     are not modelled in Phase 1 for the same reason; Direction is modelled
//     as a string and the array form is unavailable until Phase 1.5 wires up
//     a `std::variant`-aware schema pass.
//
//   - Palace's runtime uses `std::map<int, T>` keyed by Index for Energy /
//     Probe entries. The JSON wire format is always an array with Index as an
//     element member, which is what we model here.

#include <array>
#include <optional>
#include <string>
#include <vector>

#include <rfl.hpp>

#include "common.hpp"
#include "schema/utils/annotations.hpp"

namespace palace::schema
{

struct MaterialData
{
  PALACE_SCHEMA_DESC(Attributes, "Integer array of mesh domain this object applies to.",
                     AttributeList) = {};

  PALACE_SCHEMA_DESC(Permeability,
                     "Relative permeability for this material. Scalar or vector of 3 "
                     "coefficients corresponding to each of `\"MaterialAxes\"`.",
                     std::array<double, 3>) = {
    {1.0, 1.0, 1.0}
  };

  PALACE_SCHEMA_DESC(Permittivity,
                     "Relative permittivity for this material. Scalar or vector of 3 "
                     "coefficients corresponding to each of `\"MaterialAxes\"`.",
                     std::array<double, 3>) = {
    {1.0, 1.0, 1.0}
  };

  PALACE_SCHEMA_DESC(LossTan,
                     "Loss tangent for this material. Scalar or vector of 3 coefficients "
                     "corresponding to each of `\"MaterialAxes\"`.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC(Conductivity,
                     "Electrical conductivity for this material, S/m. Activates the Ohmic "
                     "loss model in this domain. Scalar or vector of 3 coefficients "
                     "corresponding to each of `\"MaterialAxes\"`.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC(LondonDepth,
                     "London penetration depth for this material, specified in mesh length "
                     "units. Activates the London equations-based model relating "
                     "superconducting current and electromagnetic fields in this domain.",
                     double) = 0.0;
};

struct CurrentDipoleData
{
  PALACE_SCHEMA_DESC(Index,
                     "Index of this current dipole source, used in postprocessing output "
                     "files.",
                     palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(Moment, "Current dipole moment magnitude, A·m.", double) = 0.0;

  PALACE_SCHEMA_DESC(Center,
                     "Coordinates of the dipole center position `[x, y, z]`, in mesh "
                     "length units.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC(Direction,
                     "Direction of the Dirac current source specifying the dipole. "
                     "Axis-aligned directions can be specified using keywords: `\"+X\"`, "
                     "`\"-X\"`, `\"+Y\"`, `\"-Y\"`, `\"+Z\"`, `\"-Z\"`. The direction can "
                     "alternatively be specified as a normalized array of three values, "
                     "for example `[0.0, 1.0, 0.0]`.",
                     std::string) = "";
};

struct DomainEnergyData
{
  PALACE_SCHEMA_DESC(Index,
                     "Index of this energy postprocessing domain, used in output files.",
                     palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(Attributes, "Integer array of mesh domain this object applies to.",
                     AttributeList) = {};
};

struct ProbeData
{
  PALACE_SCHEMA_DESC(Index, "Index of this probe, used in postprocessing output files.",
                     palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(Center, "Coordinates of this probe `[x, y, z]`, in mesh length units.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };
};

struct DomainPostData
{
  PALACE_SCHEMA_DESC(Energy,
                     "Array of objects for postprocessing domain energies. Postprocesses "
                     "the electric and magnetic field energy inside a given domain. "
                     "Results are written to `domain-E.csv` in the output directory.",
                     std::vector<DomainEnergyData>) = {};

  PALACE_SCHEMA_DESC(Probe,
                     "Array of probe points for evaluating field values at specified "
                     "locations in space. The electric field **E** and magnetic flux "
                     "density **B** are probed and written to `probe-E.csv` and "
                     "`probe-B.csv` in the output directory.",
                     std::vector<ProbeData>) = {};
};

struct DomainData
{
  PALACE_SCHEMA_DESC(Materials, "Array of material property objects.",
                     std::vector<MaterialData>) = {};

  PALACE_SCHEMA_DESC(CurrentDipole, "Array of current dipole source excitations.",
                     std::vector<CurrentDipoleData>) = {};

  PALACE_SCHEMA_DESC(Postprocessing, "Configuration for domain postprocessing.",
                     DomainPostData) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_DOMAINS_HPP
