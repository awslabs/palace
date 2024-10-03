// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_IODATA_HPP
#define PALACE_UTILS_IODATA_HPP

#include <complex>
#include "utils/configfile.hpp"

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
  double GetMeshLengthScale() const { return Lc / model.L0; }

  // Redimensionalize values for output. Outputs which depend on the fields assume a
  // characteristic reference magnetic field strength Hc such that Pc = 1 W, where Pc is the
  // characteristic reference power.
  enum class ValueType
  {
    TIME,          // [ns]
    FREQUENCY,     // [GHz]
    LENGTH,        // [m]
    IMPEDANCE,     // [Ω]
    INDUCTANCE,    // [H] = [Ωs]
    CAPACITANCE,   // [F] = [s/Ω]
    CONDUCTIVITY,  // [S/m]
    VOLTAGE,       // [V]
    CURRENT,       // [A]
    POWER,         // [W]
    ENERGY,        // [J]
    FIELD_E,       // [V/m]
    FIELD_D,       // [C/m²] = [A⋅s/m²]
    FIELD_H,       // [A/m]
    FIELD_B        // [Wb/m²] = [V⋅s/m²]
  };
  template <typename T>
  T DimensionalizeValue(ValueType type, T v) const;
  template <typename T>
  std::complex<T> DimensionalizeValue(ValueType type, std::complex<T> v) const
  {
    return {DimensionalizeValue(type, v.real()), DimensionalizeValue(type, v.imag())};
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_IODATA_HPP
