// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

using namespace Catch::Matchers;

namespace palace
{

namespace
{

double MeasureConstantTopSurface(double edge_cutoff)
{
  fem::DefaultIntegrationOrder::p_trial = 2;
  // The selected unit-square top face is split into a 4 x 4 conforming triangle mesh.
  // With cutoff 1/4, the retained domain is exactly the central 1/2 x 1/2 square.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 4, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0));
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(parallel_mesh));

  mfem::H1_FECollection h1_collection(2, 3);
  mfem::ND_FECollection nd_collection(2, 3);
  FiniteElementSpace h1_fespace(mesh, &h1_collection);
  FiniteElementSpace nd_fespace(mesh, &nd_collection);

  GridFunction electric_field(nd_fespace);
  mfem::Vector constant(3);
  constant = 0.0;
  constant[0] = 1.0;
  mfem::VectorConstantCoefficient electric_coefficient(constant);
  electric_field.Real().ProjectCoefficient(electric_coefficient);

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator material_operator({material}, periodic, ProblemType::ELECTROSTATIC,
                                     mesh);

  config::InterfaceDielectricData interface;
  interface.attributes = {6};  // z = 1 top face of MakeCartesian3D.
  interface.type = InterfaceDielectric::DEFAULT;
  interface.t = 2.0;
  interface.epsilon_r = 1.0;
  interface.tandelta = 0.0;
  interface.edge_cutoff = edge_cutoff;
  config::BoundaryPostData postprocessing;
  postprocessing.dielectric.emplace(1, interface);

  SurfacePostOperator surface_operator(postprocessing, ProblemType::ELECTROSTATIC,
                                       material_operator, h1_fespace, nd_fespace);
  return surface_operator.GetInterfaceElectricFieldEnergy(1, electric_field);
}

}  // namespace

TEST_CASE("Standard surface EdgeCutoff excludes a physical boundary band",
          "[surfacepostoperator][Serial][Parallel]")
{
  // 0.5 * thickness * epsilon * |E|^2 = 1, so energy equals retained area.
  CHECK_THAT(MeasureConstantTopSurface(0.0), WithinAbs(1.0, 1.0e-12));
  CHECK_THAT(MeasureConstantTopSurface(0.25), WithinAbs(0.25, 1.0e-11));
  // The cutoff crosses element interiors here, exercising the fixed high-order mask rather
  // than a cutoff aligned to element boundaries. The analytic retained area is 0.6^2.
  CHECK_THAT(MeasureConstantTopSurface(0.2), WithinAbs(0.36, 5.0e-3));
}

}  // namespace palace
