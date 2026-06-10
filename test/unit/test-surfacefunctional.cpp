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
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/surfacefunctional.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/labels.hpp"

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

// Reference (host, mfem-based) integration of a legacy boundary coefficient over the
// marked boundary elements, mirroring SurfacePostOperator::GetLocalSurfaceIntegral.
double RefSurfaceCoefficientIntegral(mfem::ParMesh &pmesh, mfem::Coefficient &f,
                                     const mfem::Array<int> &marker)
{
  double sum = 0.0;
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
      sum += ip.weight * T->Weight() * f.Eval(*T, ip);
    }
  }
  Mpi::GlobalSum(1, &sum, pmesh.GetComm());
  return sum;
}

// Build a 3D mesh with two material regions (attribute 1: z < 0.5, vacuum; attribute
// 2: z >= 0.5, dielectric) and interior boundary elements (attribute 7) added on the
// material interface at z = 0.5.
std::unique_ptr<Mesh> MakeInterfaceMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  // Add interior boundary elements on the material interface.
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
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

TEST_CASE("SurfaceFunctional Interface Dielectric", "[surfacefunctional][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  // Materials: vacuum for attribute 1, dielectric for attribute 2.
  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  // Project non-trivial smooth vector fields onto the (complex-valued) ND space.
  GridFunction E(nd_fespace, complex);
  mfem::VectorFunctionCoefficient fr(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1)) + x(2) * x(2);
                                       v(1) = std::cos(x(2)) + x(0);
                                       v(2) = x(0) * x(1) + 1.0;
                                     });
  E.Real().ProjectCoefficient(fr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fi(3,
                                       [](const mfem::Vector &x, mfem::Vector &v)
                                       {
                                         v(0) = x(1) * x(2) - 0.5;
                                         v(1) = std::sin(x(0)) - x(2);
                                         v(2) = std::cos(x(1)) + x(0) * x(0);
                                       });
    E.Imag().ProjectCoefficient(fi);
  }

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](InterfaceDielectric type, const mfem::Array<int> &marker)
  {
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
    auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
    {
      switch (type)
      {
        case InterfaceDielectric::DEFAULT:
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(E, mat_op, t_i,
                                                                            epsilon_i);
        case InterfaceDielectric::MA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::MS:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(
              E, mat_op, t_i, epsilon_i);
        case InterfaceDielectric::SA:
          return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>(
              E, mat_op, t_i, epsilon_i);
      }
      return {};
    };
    auto legacy = MakeLegacy();
    const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
    const double val = epr.Eval(E);
    CAPTURE(ref, val);
    if (std::abs(ref) > 0.0)
    {
      CHECK(val == Catch::Approx(ref).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(val) < 1.0e-12);
    }
  };

  SECTION("Exterior boundary")
  {
    // Boundary attribute 1 is an exterior boundary (single-sided evaluation). The
    // bottom boundary (z = 0) neighbors vacuum, the top (z = 1) neighbors dielectric.
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 6})  // z = 0 (vacuum side), z = 1 (dielectric side)
      {
        CAPTURE(attr);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
    // Boundary attribute 7 is the interior vacuum-dielectric interface (two-sided,
    // side selection by material light speed for MA/MS/SA, averaging for DEFAULT).
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[7 - 1] = 1;
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      TestType(type, marker);
    }
  }
}

}  // namespace palace
