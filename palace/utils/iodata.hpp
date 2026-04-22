// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_IODATA_HPP
#define PALACE_UTILS_IODATA_HPP

#include <sstream>
#include <nlohmann/json_fwd.hpp>
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace mfem
{

class ParMesh;

}  // namespace mfem

namespace palace
{

std::stringstream PreprocessFile(const char *filename);

//
// A parser class for processing the configuration file which controls runtime options.
//
class IoData
{
public:
  // Configuration file objects.
  config::ProblemData problem;
  config::ModelData model;
  config::DomainData domains;
  config::BoundaryData boundaries;
  config::SolverData solver;

  // Class that holds mesh scale and converts between SI quantities and normalized values.
  Units units;

private:
  bool init;

public:
  // Check configuration file options and compatibility with requested problem type.
  // Exposed for testing; not intended for general use.
  void CheckConfiguration();

  explicit IoData(const Units &units);

  // Construct from a pre-parsed and validated JSON config object.
  IoData(const nlohmann::json &config, bool print);

  // Parse configuration file, validate, and construct IoData.
  IoData(const char *filename, bool print);

  // Parse and validate a configuration file, returning the JSON object.
  static nlohmann::json ParseAndValidate(const char *filename);

  // Write resolved IoData values back into the JSON config. Call after
  // CheckConfiguration() to produce a self-describing config with all sentinels resolved.
  static void ConcretizeDefaults(const IoData &iodata, nlohmann::json &config);

  // Nondimensionalize input values for use in the solver, including the mesh coordinates.
  void NondimensionalizeInputs(mfem::ParMesh &mesh);
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
