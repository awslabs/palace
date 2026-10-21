// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/gridfunctionio.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

double FieldFunction(double x, double y) { return 2.0 + 0.5 * x - 0.25 * y; }

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir, "GridFunction field file reader",
                 "[gridfunctionio][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));

  const int nv = mesh.Get().GetNV();

  // Gmsh v2.2 file with a named $NodeData array over the mesh vertices.
  {
    std::ofstream out((temp_dir / "field.msh").string());
    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << nv << "\n";
    for (int v = 0; v < nv; v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      out << v + 1 << " " << x[0] << " " << x[1] << " 0.0\n";
    }
    out << "$EndNodes\n$NodeData\n1\n\"eps\"\n1\n0\n2\n1 " << nv << "\n";
    for (int v = 0; v < nv; v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      out << v + 1 << " " << FieldFunction(x[0], x[1]) << "\n";
    }
    out << "$EndNodeData\n";
  }

  // VTU file with the same named point data.
  {
    std::ofstream out((temp_dir / "field.vtu").string());
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\">\n"
           "<UnstructuredGrid>\n<Piece NumberOfPoints=\""
        << nv
        << "\">\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
           "format=\"ascii\">\n";
    for (int v = 0; v < nv; v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      out << x[0] << " " << x[1] << " 0.0 ";
    }
    out << "\n</DataArray>\n</Points>\n<PointData>\n<DataArray type=\"Float64\" "
           "Name=\"eps\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int v = 0; v < nv; v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      out << FieldFunction(x[0], x[1]) << " ";
    }
    out << "\n</DataArray>\n</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  }

  for (const std::string file : {"field.msh", "field.vtu"})
  {
    CAPTURE(file);
    const auto data = fieldio::ReadFieldFile((temp_dir / file).string(), "eps");
    REQUIRE(data.points.size() == static_cast<std::size_t>(nv));
    auto grid = fieldio::BuildGridFunction(mesh, data);
    mfem::Array<int> dofs;
    for (int v = 0; v < nv; v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      grid.fes->GetVertexDofs(v, dofs);
      REQUIRE(dofs.Size() == 1);
      CHECK((*grid.gf)(dofs[0]) == Catch::Approx(FieldFunction(x[0], x[1])));
    }
  }
}

TEST_CASE("MaterialFields config parsing", "[gridfunctionio][Serial]")
{
  auto setup = nlohmann::json::parse(R"({
    "Problem": { "Type": "BoundaryMode" },
    "Model": { "Mesh": "x.msh" },
    "Domains": { "Materials": [{ "Attributes": [1] }] },
    "Boundaries": { "PEC": { "Attributes": [1, 2, 3, 4] } },
    "Solver": {
      "BoundaryMode": {
        "Freq": 1.0,
        "MaterialFields": {
          "Permittivity": { "File": "eps.msh", "Name": "eps" },
          "Permeability": { "File": "mu.vtu" }
        }
      }
    }
  })");
  IoData iodata(setup, /*print=*/false);
  CHECK(iodata.solver.boundary_mode.permittivity_field.file == "eps.msh");
  CHECK(iodata.solver.boundary_mode.permittivity_field.name == "eps");
  CHECK(iodata.solver.boundary_mode.permeability_field.file == "mu.vtu");
  CHECK(iodata.solver.boundary_mode.permeability_field.name.empty());
  CHECK_FALSE(iodata.solver.boundary_mode.permittivity_field.empty());
  CHECK_FALSE(iodata.solver.boundary_mode.permeability_field.empty());
}

}  // namespace palace
