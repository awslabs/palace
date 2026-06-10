// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/surfacefunctional.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{

namespace
{

std::unique_ptr<Mesh> LoadTestMesh(MPI_Comm comm, const std::string &file)
{
  mfem::Mesh smesh(std::string(PALACE_TEST_DATA_DIR) + "/mesh/" + file, 1, 1);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Reference (host, mfem-based) computation of the surface area with the same
// quadrature order as SurfaceFunctional (mesh::GetSurfaceArea uses a different
// quadrature order, which only matters for curved boundary elements).
double RefSurfaceArea(mfem::ParMesh &pmesh, const mfem::Array<int> &marker)
{
  double area = 0.0;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = pmesh.GetBdrElementTransformation(i);
    const auto &ir = mfem::IntRules.Get(pmesh.GetBdrElementGeometry(i),
                                        fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      area += ip.weight * T->Weight();
    }
  }
  Mpi::GlobalSum(1, &area, pmesh.GetComm());
  return area;
}

// Reference (host, mfem-based) computation of ∫ |E|² dS over the marked boundary
// elements, evaluating E from the attached volume element at each quadrature point in
// the same way as the legacy BdrGridFunctionCoefficient-based postprocessing.
double RefSurfaceHCurlNorm2(mfem::ParMesh &pmesh, const mfem::ParGridFunction &E,
                            const mfem::Array<int> &marker)
{
  double sum = 0.0;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  mfem::Vector Ev(pmesh.SpaceDimension());
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = pmesh.GetBdrElementTransformation(i);
    const auto &ir = mfem::IntRules.Get(pmesh.GetBdrElementGeometry(i),
                                        fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(i, pmesh, FET, T1,
                                                                       T2, &ip);
      E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), Ev);
      sum += ip.weight * T->Weight() * (Ev * Ev);
    }
  }
  Mpi::GlobalSum(1, &sum, pmesh.GetComm());
  return sum;
}

}  // namespace

TEST_CASE("SurfaceFunctional Area", "[surfacefunctional][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  auto file = GENERATE(std::string("fichera-tet.mesh"), std::string("fichera-hex.mesh"),
                       std::string("fichera-mixed-p2.mesh"));
  CAPTURE(file);
  auto mesh = LoadTestMesh(comm, file);
  auto &pmesh = mesh->Get();

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  SECTION("All boundary attributes")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    SurfaceFunctional area(SurfaceFunctional::Kind::AREA, *mesh, marker);
    CHECK(area.Eval() == Catch::Approx(RefSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    if (file != "fichera-mixed-p2.mesh")
    {
      // For meshes with straight-sided boundary elements, also cross-check against the
      // independent implementation in mesh::GetSurfaceArea (exact for any quadrature
      // order, so quadrature rule differences don't matter).
      CHECK(area.Eval() ==
            Catch::Approx(mesh::GetSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    }
  }

  SECTION("Single boundary attribute")
  {
    for (auto attr : pmesh.bdr_attributes)
    {
      CAPTURE(attr);
      mfem::Array<int> marker(bdr_attr_max);
      marker = 0;
      marker[attr - 1] = 1;
      SurfaceFunctional area(SurfaceFunctional::Kind::AREA, *mesh, marker);
      CHECK(area.Eval() == Catch::Approx(RefSurfaceArea(pmesh, marker)).epsilon(1.0e-12));
    }
  }
}

TEST_CASE("SurfaceFunctional HCurl Norm", "[surfacefunctional][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto file = GENERATE(std::string("fichera-tet.mesh"), std::string("fichera-hex.mesh"),
                       std::string("fichera-mixed-p2.mesh"));
  auto order = GENERATE(1, 2);
  CAPTURE(file, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  auto mesh = LoadTestMesh(comm, file);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  // Project a non-trivial smooth vector field onto the ND space. The reference value is
  // computed from the same projected grid function, so the comparison is exact up to
  // quadrature evaluation differences (not projection error).
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient f(3,
                                    [](const mfem::Vector &x, mfem::Vector &v)
                                    {
                                      v(0) = std::sin(x(1)) + x(2) * x(2);
                                      v(1) = std::cos(x(2)) + x(0);
                                      v(2) = x(0) * x(1) + 1.0;
                                    });
  E.ProjectCoefficient(f);

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  SECTION("All boundary attributes")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    SurfaceFunctional norm2(SurfaceFunctional::Kind::HCURL_NORM2, *mesh, marker,
                            &nd_fespace);
    const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
    CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
  }

  SECTION("Single boundary attribute")
  {
    for (auto attr : pmesh.bdr_attributes)
    {
      CAPTURE(attr);
      mfem::Array<int> marker(bdr_attr_max);
      marker = 0;
      marker[attr - 1] = 1;
      SurfaceFunctional norm2(SurfaceFunctional::Kind::HCURL_NORM2, *mesh, marker,
                              &nd_fespace);
      const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
      CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
    }
  }
}

}  // namespace palace
