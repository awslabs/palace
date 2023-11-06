// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_IODATA_HPP
#define PALACE_UTILS_IODATA_HPP

#include "utils/configfile.hpp"

namespace mfem
{

class ParMesh;

}  // namespace mfem

namespace palace
{

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

private:
  // Characteristic reference length [m] and time [ns] for nondimensionalization.
  double Lc, tc;
  bool init;

  // Check configuration file options and compatibility with requested problem type.
  void CheckConfiguration();

public:
  // Parse command line arguments and override options defaults.
  IoData(const char *filename, bool print);

  // Nondimensionalize input values for use in the solver, including the mesh coordinates.
  void NondimensionalizeInputs(mfem::ParMesh &mesh);

  // Return the mesh scaling factor in units model.L0 x [m] for mesh IO.
  double GetLengthScale() const { return Lc / model.L0; }

  // Redimensionalize values for output.
  enum class ValueType
  {
    TIME,         // [ns]
    FREQUENCY,    // [GHz]
    LENGTH,       // [m]
    IMPEDANCE,    // [Ω]
    INDUCTANCE,   // [H]
    CAPACITANCE,  // [F]
    CONDUCTIVITY  // [S/m]
  };
  double DimensionalizeValue(ValueType type, double v) const;
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
