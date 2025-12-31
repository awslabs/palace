// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_IODATA_HPP
#define PALACE_UTILS_IODATA_HPP

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
  bool explicit_units;  // Track if units were explicitly provided via constructor

public:
  // Check configuration file options and compatibility with requested problem type. Should
  // not be called by user, but temporarily made public for testing.
  void CheckConfiguration();

  IoData(const Units &units) : units(units), init(false), explicit_units(true) {}

  // Parse command line arguments and override options defaults.
  IoData(const char *filename, bool print);

  // Nondimensionalize input values for use in the solver, including the mesh coordinates.
  void NondimensionalizeInputs(mfem::ParMesh &mesh);
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
