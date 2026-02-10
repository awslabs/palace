// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/coefficient.hpp"
#include "fem/mesh.hpp"

using namespace palace;

TEST_CASE("BdrFieldCoefficient interior average", "[coefficient][Serial]")
{
  // 2x1 quad mesh: two elements sharing an interior face at x=1.
  mfem::Mesh smesh =
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);

  // Add interior boundary at x=1 (the shared face between elements).
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)  // Interior face.
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      smesh.AddBdrElement(face_elem);
      smesh.SetBdrAttribute(smesh.GetNBE() - 1, 99);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  palace::Mesh mesh(std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, smesh));

  // H1 space, order 1.
  mfem::H1_FECollection fec(1, 2);
  mfem::ParFiniteElementSpace fes(&mesh.Get(), &fec);
  mfem::ParGridFunction u(&fes);

  // Set u = x. At interior face (x=1): elem0 centroid ~0.5, elem1 centroid ~1.5.
  // Proper average = 1.0. Bug (comma op) would give 0.5 * (value from elem1).
  mfem::FunctionCoefficient x_coord([](const mfem::Vector &x) { return x(0); });
  u.ProjectCoefficient(x_coord);

  BdrFieldCoefficient coeff(mesh, u);

  // Evaluate on interior boundary (attr 99).
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (mesh.GetBdrAttribute(be) == 99)
    {
      auto *T = mesh.GetBdrElementTransformation(be);
      mfem::IntegrationPoint ip;
      ip.Set1w(0.5, 1.0);  // Midpoint of the face.
      double val = coeff.Eval(*T, ip);
      // At x=1 face, both elements should give u=1.0 (since u=x and face is at x=1).
      // Average is 1.0. With comma bug: 0.5 * 1.0 = 0.5.
      CHECK_THAT(val, Catch::Matchers::WithinAbs(1.0, 1e-12));
      return;
    }
  }
  FAIL("Interior boundary not found");
}
