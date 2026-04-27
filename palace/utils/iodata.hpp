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

  // Take parsed json and override options defaults.
  explicit IoData(nlohmann::json &&config, bool print = false);

  // Parse command line arguments and override options defaults.
  explicit IoData(const char *filename, bool print);

  // Nondimensionalize input values for use in the solver, including the mesh coordinates.
  // Requires model.Lc > 0; the caller is responsible for populating Lc (from the config
  // or via mesh::ComputeReferenceLength) so that every rank agrees on the scale before
  // this is called. The mesh pointer may be null on ranks that don't hold a copy of the
  // serial mesh in the byte-string distribution path; iodata is still scaled everywhere.
  // Pure on the MPI side: no collectives.
  void NondimensionalizeInputs(std::unique_ptr<mfem::Mesh> &mesh);
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
