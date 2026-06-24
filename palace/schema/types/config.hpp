// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_CONFIG_HPP
#define PALACE_SCHEMA_TYPES_CONFIG_HPP

// Top-level Palace configuration: aggregates all five sections Palace parses
// from a single JSON config file (Problem, Model, Domains, Boundaries,
// Solver). PR 716 marks all five sections as required at the top level; the
// annotated types here carry defaults so `PalaceConfig{}` is default-constructible
// (needed by the schema::utils compile-time type-graph walk). `SchemaVersion`
// remains optional for backwards compatibility with existing Palace inputs.

#include <string>

#include <rfl.hpp>

#include "boundaries.hpp"
#include "domains.hpp"
#include "model.hpp"
#include "problem.hpp"
#include "schema/utils/annotations.hpp"
#include "solver.hpp"

namespace palace::schema
{

struct PalaceConfiguration
{
  PALACE_SCHEMA_DESC(SchemaVersion,
                     "Schema version in SchemaVer format: `MODEL-REVISION-ADDITION` "
                     "(three non-negative integers, hyphen-separated, e.g. "
                     "`1-0-0`). See https://docs.snowplow.io/docs/api-reference/"
                     "iglu/common-architecture/schemaver/.",
                     PALACE_SCHEMA_PATTERN("[0-9]+-[0-9]+-[0-9]+", "schema-ver")) = "1-0-0";

  PALACE_SCHEMA_DESC_REQUIRED(Problem,
                              "Top-level configuration for the simulation type and output.",
                              ::palace::schema::Problem) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Model, "Mesh and model configuration.",
                              ::palace::schema::Model) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Domains, "Material and domain configuration.",
                              ::palace::schema::Domains) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Boundaries, "Boundary condition configuration.",
                              ::palace::schema::Boundaries) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Solver, "Solver configuration for all simulation types.",
                              ::palace::schema::Solver) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_CONFIG_HPP
