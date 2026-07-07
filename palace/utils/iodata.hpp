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

// Reflect-cpp schema types, parsed from the configuration and converted to the runtime
// config:: structs. Forward-declared (used only as parameter/return types here) so
// reflect-cpp stays out of this widely-included header; definitions live in schema/types.
namespace schema
{

struct Problem;
struct PalaceConfiguration;

}  // namespace schema

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

  // Shared body of the JSON constructors: reads Model/Domains/Boundaries/Solver from
  // `config`, rejects unknown sections, and runs CheckConfiguration(). `problem` must
  // already be populated by the caller.
  void SetUpFromJson(const nlohmann::json &config, bool print);

public:
  // Check configuration file options and compatibility with requested problem type.
  // Exposed for testing; not intended for general use.
  void CheckConfiguration();

  explicit IoData(const Units &units);

  // Construct from a pre-parsed and validated JSON config object together with the schema
  // configuration object parsed from the same source. The Problem section is initialized
  // from `pconfig`; the remaining sections are still read from `config` (to be migrated to
  // `pconfig` in follow-ups). This is the seam for retiring the nlohmann::json path.
  IoData(const nlohmann::json &config, const schema::PalaceConfiguration &pconfig,
         bool print);

  // Parse configuration file, validate, and construct IoData.
  IoData(const char *filename, bool print);

  // Parse and validate a configuration file, returning the JSON object.
  static nlohmann::json ParseAndValidate(const char *filename);

  // Parse a configuration directly into the reflect-cpp schema configuration object. The
  // file overload preprocesses (strips comments, expands integer ranges) independently of
  // nlohmann::json so the JSON path can eventually be removed; the json overload parses an
  // already-loaded config object (used where the JSON has already been read/validated).
  static schema::PalaceConfiguration ParsePalaceConfiguration(const char *filename);
  static schema::PalaceConfiguration ParsePalaceConfiguration(const nlohmann::json &config);

  // Derive the runtime ProblemData from a parsed "Problem" schema object. Exposed for
  // testing the schema bridge; not intended for general use.
  static config::ProblemData FromSchema(const schema::Problem &problem);

  // Return the user's config with any entries that were absent filled in from the
  // resolved IoData. User-provided values pass through untouched; only missing keys
  // are added. Call after CheckConfiguration() so the filled values are concrete.
  static nlohmann::json ConcretizeDefaults(const IoData &iodata, nlohmann::json config);

  // Concretize defaults from `raw_config` and write the result to
  // `<problem.output>/<input_stem>_resolved.json` so users have a self-contained record
  // of every Palace decision. The `_resolved` suffix (derived from `input_config_path`)
  // ensures the sidecar never overwrites the user's input config, even when `"Output"`
  // resolves to the input's own directory. The caller is responsible for ensuring this is
  // invoked on a single process (e.g. MPI rank 0). Aborts on filesystem failure.
  void WriteResolvedConfig(const nlohmann::json &raw_config,
                           const std::string &input_config_path) const;

  // Nondimensionalize input values and mesh coordinates. Requires model.Lc > 0 (caller
  // populates it from the config or via mesh::ComputeReferenceLength). `mesh` may be
  // null on ranks that don't hold the serial mesh; iodata is scaled on every rank. No
  // MPI.
  void NondimensionalizeInputs(std::unique_ptr<mfem::Mesh> &mesh);
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
