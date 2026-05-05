// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_PROBLEM_HPP
#define PALACE_SCHEMA_TYPES_PROBLEM_HPP

// Mirrors palace::config::ProblemData (palace/utils/configfile.hpp). Field
// names match PR 716's scripts/schema/config/problem.json exactly; C++ members
// are PascalCase so reflect-cpp emits the right JSON keys with no rfl::Rename.

#include <string>

#include <rfl.hpp>

#include "schema/utils/annotations.hpp"
#include "common.hpp"

namespace palace::schema {

struct OutputFormatsData {
    PALACE_SCHEMA_DESC(Paraview,
             "Set to `true` to output fields in [ParaView](https://www.paraview.org/) format.",
             bool) = true;

    PALACE_SCHEMA_DESC(GridFunction,
             "Set to `true` to output fields in MFEM grid function format for "
             "visualization with [GLVis](https://glvis.org/).",
             bool) = false;
};

struct ProblemData {
    PALACE_SCHEMA_DESC(Type, "Controls the simulation type.", ProblemType) = ProblemType::Driven;

    PALACE_SCHEMA_DESC(Verbose, "Controls the level of log file printing.",
             palace::schema::utils::Min<int, 0>) = 1;

    PALACE_SCHEMA_DESC(Output, "Directory path for saving postprocessing outputs.",
             std::string) = "";

    PALACE_SCHEMA_DESC(OutputFormats, "Configures the field output formats.",
             OutputFormatsData) = {};
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_PROBLEM_HPP
