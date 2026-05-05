// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_CONFIG_HPP
#define PALACE_SCHEMA_TYPES_CONFIG_HPP

// Top-level Palace configuration: aggregates all five sections Palace parses
// from a single JSON config file (Problem, Model, Domains, Boundaries,
// Solver). PR 716 marks all five as required at the top level; the annotated
// types here carry defaults so `PalaceConfig{}` is default-constructible
// (needed by the schema::utils compile-time type-graph walk). The emitted schema's
// `required` list still names all five because of the pruning rule — Palace
// enforces presence via MFEM_VERIFY at SetUp time.

#include <rfl.hpp>

#include "schema/utils/annotations.hpp"
#include "boundaries.hpp"
#include "domains.hpp"
#include "model.hpp"
#include "problem.hpp"
#include "solver.hpp"

namespace palace::schema {

struct PalaceConfig {
    PALACE_SCHEMA_DESC(Problem,
             "Top-level configuration for the simulation type and output.",
             ProblemData) = {};

    PALACE_SCHEMA_DESC(Model, "Mesh and model configuration.", ModelData) = {};

    PALACE_SCHEMA_DESC(Domains, "Material and domain configuration.", DomainData) = {};

    PALACE_SCHEMA_DESC(Boundaries, "Boundary condition configuration.",
             BoundaryData) = {};

    PALACE_SCHEMA_DESC(Solver, "Solver configuration for all simulation types.",
             SolverData) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_CONFIG_HPP
