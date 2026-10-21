// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GRIDFUNCTION_IO_HPP
#define PALACE_UTILS_GRIDFUNCTION_IO_HPP

#include <array>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class Mesh;

namespace fieldio
{

//
// Reading of named scalar field data from standard mesh files into grid functions on a
// solve mesh. MFEM's own mesh readers discard field data, so the named arrays are extracted
// here (Gmsh v2.2 $NodeData node data, or VTK/VTU point data) and mapped to the solve mesh
// by coordinate matching.
//

// Named scalar field as (coordinates, values) pairs.
struct FieldData
{
  std::vector<std::array<double, 3>> points;
  std::vector<double> values;
};

// A scalar grid function together with the finite element collection and space that own it.
struct FieldGridFunction
{
  std::unique_ptr<mfem::H1_FECollection> fec;
  std::unique_ptr<mfem::ParFiniteElementSpace> fes;
  std::unique_ptr<mfem::ParGridFunction> gf;
};

// Read a named scalar array from a Gmsh v2.2 (.msh) file or a VTK/VTU point-data file.
FieldData ReadFieldFile(const std::string &filename, const std::string &name);

// Build a P1 (nodal) grid function on `mesh` from `data` by matching vertex coordinates.
// Aborts if any local vertex has no matching point.
FieldGridFunction BuildGridFunction(const Mesh &mesh, const FieldData &data);

}  // namespace fieldio
}  // namespace palace

#endif  // PALACE_UTILS_GRIDFUNCTION_IO_HPP
