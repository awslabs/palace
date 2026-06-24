// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/coefficient.hpp"
#include "fem/domain_field_evaluator.hpp"
#include "fem/facenbrexchange.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/interpolator.hpp"
#include "fem/mesh.hpp"
#include "fem/output_functionals.hpp"
#include "fixtures.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/strattonchu.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/labels.hpp"
#include "utils/units.hpp"

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

// Build a 2D mesh with two material regions (attribute 1: y < 0.5, vacuum; attribute
// 2: y >= 0.5, dielectric) and interior boundary elements (attribute 7) added on the
// material interface at y = 0.5.
std::unique_ptr<Mesh> MakeInterfaceMesh2D(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian2D(2, 2, elem_type, false, 1.0, 1.0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(2);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(1) < 0.5) ? 1 : 2);
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

class TestBoundaryModePostOperator : public PostOperator<ProblemType::BOUNDARYMODE>
{
public:
  using PostOperator<ProblemType::BOUNDARYMODE>::PostOperator;
  const Measurement &Cache() const { return measurement_cache; }
};

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

// Build a 3D mesh with two material regions and interior boundary elements (attribute
// 7) on the material interface at z = 0.5, then apply nonconformal refinement to a
// subset of elements so that marked surfaces (exterior and interior) contain
// nonconformal (master/slave) faces.
std::unique_ptr<Mesh> MakeNCInterfaceMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
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
  smesh.EnsureNCMesh(true);

  // Nonconformally refine the elements in the x < 0.5 half (crosses both the interior
  // interface and the exterior boundaries), creating master/slave faces on the marked
  // surfaces.
  mfem::Array<int> refs;
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    if (center(0) < 0.5)
    {
      refs.Append(e);
    }
  }
  smesh.GeneralRefinement(refs, 1);

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

}  // namespace

TEST_CASE("SurfaceFunctional Area", "[surfacefunctional][Serial][Parallel][GPU]")
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

TEST_CASE("SurfaceFunctional HCurl Norm", "[surfacefunctional][Serial][Parallel][GPU]")
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
                            &nd_fespace.Get());
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
                              &nd_fespace.Get());
      const double ref = RefSurfaceHCurlNorm2(pmesh, E, marker);
      CHECK(norm2.Eval(&E) == Catch::Approx(ref).epsilon(1.0e-10));
    }
  }
}

TEST_CASE("SurfaceFunctional Interface Dielectric", "[surfacefunctional][Serial][GPU]")
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

TEST_CASE("SurfaceFunctional Interface Dielectric 2D", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh2D(comm, elem_type);
  auto &pmesh = mesh->Get();

  // Materials: vacuum for attribute 1, dielectric for attribute 2. Use diagonal
  // anisotropy to exercise the 2x2 coefficient unpacking path in the libCEED kernels.
  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 3.2;
  dielectric.epsilon_r.s[2] = 5.1;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);

  GridFunction E(nd_fespace, complex);
  mfem::VectorFunctionCoefficient fr(2,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(1.7 * x(1)) + 0.25 * x(0);
                                       v(1) = std::cos(0.3 + x(0)) + x(1) * x(1);
                                     });
  E.Real().ProjectCoefficient(fr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fi(2,
                                       [](const mfem::Vector &x, mfem::Vector &v)
                                       {
                                         v(0) = x(0) * x(1) - 0.2;
                                         v(1) = std::sin(x(0) + 0.5 * x(1));
                                       });
    E.Imag().ProjectCoefficient(fi);
  }

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](InterfaceDielectric type, const mfem::Array<int> &marker)
  {
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
    REQUIRE(epr.IsValid());
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
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 3})  // bottom (vacuum side), top (dielectric side)
      {
        CAPTURE(attr);
        REQUIRE(pmesh.bdr_attributes.Find(attr) >= 0);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
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

TEST_CASE_METHOD(test::SharedTempDir,
                 "PostOperator boundary mode 2D interface dielectric matches legacy",
                 "[postoperator][surfacefunctional][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  REQUIRE(Mpi::Size(comm) == 1);
  constexpr int order = 2;

  IoData iodata(Units(1.0, 1.0));
  iodata.problem.type = ProblemType::BOUNDARYMODE;
  iodata.problem.verbose = 0;
  iodata.problem.output = temp_dir.string();
  iodata.problem.output_formats.paraview = false;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;
  iodata.solver.order = order;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.n_post = 1;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 50;

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 4.0;
  dielectric.epsilon_r.s[1] = 6.0;
  dielectric.epsilon_r.s[2] = 8.0;
  iodata.domains.materials = {vacuum, dielectric};
  iodata.boundaries.pec.attributes = {1, 2, 3, 4};

  config::InterfaceDielectricData interface;
  interface.attributes = {7};
  interface.type = InterfaceDielectric::SA;
  interface.t = 2.0e-3;
  interface.epsilon_r = 10.0;
  interface.tandelta = 3.0e-4;
  iodata.boundaries.postpro.dielectric.emplace(1, interface);

  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(MakeInterfaceMesh2D(comm, mfem::Element::TRIANGLE));
  auto &pmesh = mesh.front()->Get();
  MaterialOperator mat_op(iodata, *mesh.front());
  BoundaryModeOperator fem_op(iodata, mesh, mat_op);
  TestBoundaryModePostOperator post_op(iodata, fem_op);

  GridFunction E(fem_op.GetNDSpace(), true), En(fem_op.GetH1Space(), true);
  mfem::VectorFunctionCoefficient etr(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + 0.3 * x(0);
                                        v(1) = std::cos(x(0)) + x(0) * x(1);
                                      });
  mfem::VectorFunctionCoefficient eti(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(0) - 0.2 * x(1);
                                        v(1) = std::sin(x(0) + x(1));
                                      });
  mfem::FunctionCoefficient enr(
      [](const mfem::Vector &x) { return 0.5 + x(0) * x(1); });
  mfem::FunctionCoefficient eni(
      [](const mfem::Vector &x) { return std::cos(x(0)) - 0.1 * x(1); });
  E.Real().ProjectCoefficient(etr);
  E.Imag().ProjectCoefficient(eti);
  En.Real().ProjectCoefficient(enr);
  En.Imag().ProjectCoefficient(eni);

  ComplexVector et(fem_op.GetNDTrueVSize()), en(fem_op.GetH1TrueVSize());
  E.Real().GetTrueDofs(et.Real());
  E.Imag().GetTrueDofs(et.Imag());
  En.Real().GetTrueDofs(en.Real());
  En.Imag().GetTrueDofs(en.Imag());

  // Public PostOperator entry point under test: BoundaryMode computes its derived B
  // fields internally, then measures configured interface dielectric participation.
  const double omega = 2.0;
  const std::complex<double> kn(0.7, 0.05);
  post_op.MeasureAndPrintAll(0, et, en, kn, omega, 0.0, 0.0, 1);

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 0;
  marker[7 - 1] = 1;
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> legacy(
      E, mat_op, interface.t, interface.epsilon_r);
  const double interface_energy = RefSurfaceCoefficientIntegral(pmesh, legacy, marker);
  DomainPostOperator domain_post(iodata.domains.postpro, mat_op, fem_op.GetNDSpace(),
                                 fem_op.GetCurlSpace());
  const double domain_energy = domain_post.GetElectricFieldEnergy(E);
  const double participation = interface_energy / domain_energy;

  const auto &cache = post_op.Cache();
  REQUIRE(cache.interface_eps_i.size() == 1);
  const auto &measured = cache.interface_eps_i.front();
  CHECK(measured.idx == 1);
  CHECK(measured.energy == Catch::Approx(interface_energy).epsilon(1.0e-10));
  CHECK(measured.energy_participation == Catch::Approx(participation).epsilon(1.0e-10));
  CHECK(measured.quality_factor ==
        Catch::Approx(1.0 / (participation * interface.tandelta)).epsilon(1.0e-10));
}

TEST_CASE("SurfaceFunctional Surface Flux", "[surfacefunctional][Serial][GPU]")
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

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  // Project non-trivial smooth vector fields onto the (complex-valued) ND and RT
  // spaces.
  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  mfem::Vector x0(3);
  x0(0) = 0.4;
  x0(1) = 0.6;
  x0(2) = -0.2;  // Off-center so the orientation sign flip logic is exercised
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  auto TestType = [&](SurfaceFlux type, bool two_sided, const mfem::Array<int> &marker)
  {
    const GridFunction *E_use =
        (type == SurfaceFlux::ELECTRIC || type == SurfaceFlux::POWER) ? &E : nullptr;
    const GridFunction *B_use =
        (type == SurfaceFlux::MAGNETIC || type == SurfaceFlux::POWER) ? &B : nullptr;
    SurfaceFunctional flux(*mesh, marker, E_use ? &nd_fespace.Get() : nullptr,
                           B_use ? &rt_fespace.Get() : nullptr, mat_op, type, two_sided,
                           x0);

    // Legacy reference following SurfacePostOperator::GetSurfaceFlux.
    auto MakeLegacy =
        [&](const mfem::ParGridFunction *Er,
            const mfem::ParGridFunction *Br) -> std::unique_ptr<mfem::Coefficient>
    {
      switch (type)
      {
        case SurfaceFlux::ELECTRIC:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
              Er, Br, mat_op, two_sided, x0);
        case SurfaceFlux::MAGNETIC:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>(
              Er, Br, mat_op, two_sided, x0);
        case SurfaceFlux::POWER:
          return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
              Er, Br, mat_op, two_sided, x0);
      }
      return {};
    };
    auto legacy_re = MakeLegacy(E_use ? &E.Real() : nullptr, B_use ? &B.Real() : nullptr);
    std::complex<double> ref(RefSurfaceCoefficientIntegral(pmesh, *legacy_re, marker), 0.0);
    if (complex)
    {
      auto legacy_im = MakeLegacy(E_use ? &E.Imag() : nullptr, B_use ? &B.Imag() : nullptr);
      const double ref_im = RefSurfaceCoefficientIntegral(pmesh, *legacy_im, marker);
      if (type == SurfaceFlux::POWER)
      {
        ref += ref_im;
      }
      else
      {
        ref.imag(ref_im);
      }
    }

    const std::complex<double> val = flux.EvalFlux(E_use, B_use);
    CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
    auto CheckPart = [](double v, double r)
    {
      if (std::abs(r) > 1.0e-12)
      {
        CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
      }
      else
      {
        CHECK(std::abs(v) < 1.0e-10);
      }
    };
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), ref.imag());
  };

  SECTION("Exterior boundary")
  {
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      CAPTURE(static_cast<int>(type));
      for (int attr : {1, 6})  // z = 0 (vacuum side), z = 1 (dielectric side)
      {
        CAPTURE(attr);
        mfem::Array<int> marker(bdr_attr_max);
        marker = 0;
        marker[attr - 1] = 1;
        TestType(type, false, marker);
      }
    }
  }

  SECTION("Interior material interface")
  {
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[7 - 1] = 1;
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      CAPTURE(static_cast<int>(type));
      for (bool two_sided : {false, true})
      {
        CAPTURE(two_sided);
        TestType(type, two_sided, marker);
      }
    }
  }
}

TEST_CASE("SurfaceFunctional AtPoints surface flux matches mapped path",
          "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  constexpr auto elem_type = mfem::Element::TETRAHEDRON;
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 3.1;
  dielectric.epsilon_r.s[2] = 2.4;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.8;
  dielectric.mu_r.s[2] = 2.2;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  mfem::Array<int> marker(pmesh.bdr_attributes.Max());
  marker = 0;
  marker[7 - 1] = 1;
  mfem::Vector x0(3);
  x0(0) = 0.4;
  x0(1) = 0.6;
  x0(2) = -0.2;

  const char *old_disable = std::getenv("PALACE_SURFACE_DISABLE_ATPOINTS");
  auto RestoreEnv = [&]()
  {
    if (old_disable)
    {
      setenv("PALACE_SURFACE_DISABLE_ATPOINTS", old_disable, 1);
    }
    else
    {
      unsetenv("PALACE_SURFACE_DISABLE_ATPOINTS");
    }
  };
  auto Eval = [&](SurfaceFlux type, bool two_sided, bool disable_at_points)
  {
    if (disable_at_points)
    {
      setenv("PALACE_SURFACE_DISABLE_ATPOINTS", "1", 1);
    }
    else
    {
      unsetenv("PALACE_SURFACE_DISABLE_ATPOINTS");
    }
    const GridFunction *E_use =
        (type == SurfaceFlux::ELECTRIC || type == SurfaceFlux::POWER) ? &E : nullptr;
    const GridFunction *B_use =
        (type == SurfaceFlux::MAGNETIC || type == SurfaceFlux::POWER) ? &B : nullptr;
    SurfaceFunctional flux(*mesh, marker, E_use ? &nd_fespace.Get() : nullptr,
                           B_use ? &rt_fespace.Get() : nullptr, mat_op, type, two_sided,
                           x0);
    REQUIRE(flux.IsValid());
    return flux.EvalFlux(E_use, B_use);
  };

  for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
  {
    for (bool two_sided : {false, true})
    {
      CAPTURE(static_cast<int>(type), two_sided);
      const auto mapped = Eval(type, two_sided, true);
      const auto maybe_at_points = Eval(type, two_sided, false);
      CAPTURE(mapped.real(), mapped.imag(), maybe_at_points.real(), maybe_at_points.imag());
      CHECK(maybe_at_points.real() ==
            Catch::Approx(mapped.real()).epsilon(1.0e-10).margin(1.0e-12));
      CHECK(maybe_at_points.imag() ==
            Catch::Approx(mapped.imag()).epsilon(1.0e-10).margin(1.0e-12));
    }
  }
  RestoreEnv();
}

TEST_CASE("SurfaceFunctional Complex Power", "[surfacefunctional][Serial][GPU]")
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

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  // Legacy reference: replicate LumpedPortData::GetPower exactly (linear form from
  // BdrSurfaceCurrentVectorCoefficient dotted with the E field expansion).
  auto LegacyPower = [&](int attr, const mfem::Array<int> &attr_marker)
  {
    mfem::Array<int> attr_list(1);
    attr_list[0] = attr;
    std::complex<double> dot;
    {
      RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbr_c(
          attr_list, B.Real(), mat_op);
      mfem::LinearForm pr(&nd_fespace.Get());
      pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr_c),
                               const_cast<mfem::Array<int> &>(attr_marker));
      pr.UseFastAssembly(false);
      pr.UseDevice(false);
      pr.Assemble();
      pr.UseDevice(true);
      dot =
          -(pr * E.Real()) + (complex ? std::complex<double>(0.0, -(pr * E.Imag())) : 0.0);
    }
    if (complex)
    {
      RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbi_c(
          attr_list, B.Imag(), mat_op);
      mfem::LinearForm pi(&nd_fespace.Get());
      pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi_c),
                               const_cast<mfem::Array<int> &>(attr_marker));
      pi.UseFastAssembly(false);
      pi.UseDevice(false);
      pi.Assemble();
      pi.UseDevice(true);
      dot += std::complex<double>(-(pi * E.Imag()), pi * E.Real());
    }
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  };

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Vector x0(3);
  x0 = 0.0;

  auto CheckPart = [](double v, double r)
  {
    if (std::abs(r) > 1.0e-12)
    {
      CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(v) < 1.0e-10);
    }
  };

  for (int attr : {1, 6, 7})  // Exterior (vacuum side), exterior (dielectric), interior
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    const std::complex<double> ref = LegacyPower(attr, marker);
    const std::complex<double> val = power.EvalComplexPower(E, B);
    CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), ref.imag());
  }

  SECTION("Batched by boundary attribute")
  {
    mfem::Array<int> marker(bdr_attr_max), attr_to_bin(bdr_attr_max);
    marker = 0;
    attr_to_bin = -1;
    std::vector<int> attrs = {1, 6, 7};
    for (std::size_t i = 0; i < attrs.size(); i++)
    {
      marker[attrs[i] - 1] = 1;
      attr_to_bin[attrs[i] - 1] = static_cast<int>(i);
    }
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    const auto values = power.EvalComplexPowerByAttribute(
        E, B, attr_to_bin, static_cast<int>(attrs.size()));
    REQUIRE(values.size() == attrs.size());
    for (std::size_t i = 0; i < attrs.size(); i++)
    {
      mfem::Array<int> single_marker(bdr_attr_max);
      single_marker = 0;
      single_marker[attrs[i] - 1] = 1;
      const std::complex<double> ref = LegacyPower(attrs[i], single_marker);
      CAPTURE(attrs[i], ref.real(), ref.imag(), values[i].real(), values[i].imag());
      CheckPart(values[i].real(), ref.real());
      CheckPart(values[i].imag(), ref.imag());
    }
  }
}

// Benchmark comparing the legacy mfem::Coefficient-based measurement paths against the
// libCEED surface functionals. Hidden by default; run explicitly with:
//   ./palace-unit-tests "[surfacefunctional-bench]" [--device cuda]
// Mesh size can be controlled with the PALACE_BENCH_N environment variable (default
// 20, i.e. a 20x20x20 hex mesh).
TEST_CASE("SurfaceFunctional Benchmark", "[.][surfacefunctional-bench][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  const int order = 2;
  const int N =
      std::getenv("PALACE_BENCH_N") ? std::atoi(std::getenv("PALACE_BENCH_N")) : 20;
  const int n_reps = 10;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Two-material N x N x N hex mesh with an interior interface (attribute 7).
  auto mesh = [&]()
  {
    mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(N, N, N, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(3);
      smesh.GetElementCenter(e, center);
      smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
    }
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
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
    return std::make_unique<Mesh>(std::move(pmesh));
  }();
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec),
      h1_fespace(*mesh, &h1_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  const int bdr_attr_max = pmesh.bdr_attributes.Max();
  mfem::Array<int> marker(bdr_attr_max);
  marker = 0;
  marker[7 - 1] = 1;
  const double t_i = 2.0e-3, epsilon_i = 10.0;

  auto TimeIt = [&](auto &&f, int reps)
  {
    // Warm up (also triggers any JIT compilation on GPU backends).
    f();
    double t_min = 1.0e30, t_sum = 0.0;
    for (int r = 0; r < reps; r++)
    {
      auto t0 = std::chrono::steady_clock::now();
      f();
      auto t1 = std::chrono::steady_clock::now();
      const double t = std::chrono::duration<double, std::milli>(t1 - t0).count();
      t_min = std::min(t_min, t);
      t_sum += t;
    }
    return std::make_pair(t_min, t_sum / reps);
  };

  std::cout << "\n=== SurfaceFunctional benchmark: N = " << N << " (" << pmesh.GetNE()
            << " elements, " << N * N << " interface boundary elements, ND order " << order
            << ", " << nd_fespace.GlobalTrueVSize() << " ND dofs) ===\n";

  // 1. Interface dielectric EPR (MA type) over the interior interface.
  {
    volatile double sink;
    auto legacy = [&]()
    {
      // Replicates SurfacePostOperator::GetInterfaceElectricFieldEnergy +
      // GetLocalSurfaceIntegral.
      InterfaceDielectricCoefficient<InterfaceDielectric::MA> f(E, mat_op, t_i, epsilon_i);
      mfem::LinearForm s(&h1_fespace.Get());
      s.AddBoundaryIntegrator(new BoundaryLFIntegrator(f),
                              const_cast<mfem::Array<int> &>(marker));
      s.UseFastAssembly(false);
      s.UseDevice(false);
      s.Assemble();
      s.UseDevice(true);
      double dot = linalg::LocalSum(s);
      Mpi::GlobalSum(1, &dot, comm);
      sink = dot;
    };
    auto t0 = std::chrono::steady_clock::now();
    SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, InterfaceDielectric::MA, t_i,
                          epsilon_i);
    auto t1 = std::chrono::steady_clock::now();
    const double t_constr = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto ceed_path = [&]() { sink = epr.Eval(E); };

    // Verify agreement before timing.
    InterfaceDielectricCoefficient<InterfaceDielectric::MA> f(E, mat_op, t_i, epsilon_i);
    const double ref = RefSurfaceCoefficientIntegral(pmesh, f, marker);
    const double val = epr.Eval(E);
    REQUIRE(val == Catch::Approx(ref).epsilon(1.0e-10));

    auto [legacy_min, legacy_avg] = TimeIt(legacy, n_reps);
    auto [ceed_min, ceed_avg] = TimeIt(ceed_path, n_reps);
    std::cout << "Interface EPR (MA), complex E, per measurement:\n"
              << "  legacy coefficient path: min " << legacy_min << " ms, avg "
              << legacy_avg << " ms\n"
              << "  libCEED functional:      min " << ceed_min << " ms, avg " << ceed_avg
              << " ms  (construction, once: " << t_constr << " ms)\n"
              << "  speedup (avg): " << legacy_avg / ceed_avg << "x\n";
  }

  // 2. Port power over the interior interface (driven solver inner loop quantity).
  {
    volatile double sink;
    mfem::Array<int> attr_list(1);
    attr_list[0] = 7;
    auto legacy = [&]()
    {
      // Replicates LumpedPortData::GetPower.
      std::complex<double> dot;
      {
        RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbr_c(
            attr_list, B.Real(), mat_op);
        mfem::LinearForm pr(&nd_fespace.Get());
        pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr_c),
                                 const_cast<mfem::Array<int> &>(marker));
        pr.UseFastAssembly(false);
        pr.UseDevice(false);
        pr.Assemble();
        pr.UseDevice(true);
        dot = std::complex<double>(-(pr * E.Real()), -(pr * E.Imag()));
      }
      {
        RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient> fbi_c(
            attr_list, B.Imag(), mat_op);
        mfem::LinearForm pi(&nd_fespace.Get());
        pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi_c),
                                 const_cast<mfem::Array<int> &>(marker));
        pi.UseFastAssembly(false);
        pi.UseDevice(false);
        pi.Assemble();
        pi.UseDevice(true);
        dot += std::complex<double>(-(pi * E.Imag()), pi * E.Real());
      }
      Mpi::GlobalSum(1, &dot, comm);
      sink = dot.real();
    };
    mfem::Vector x0(3);
    x0 = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    auto t1 = std::chrono::steady_clock::now();
    const double t_constr = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto ceed_path = [&]() { sink = power.EvalComplexPower(E, B).real(); };

    auto [legacy_min, legacy_avg] = TimeIt(legacy, n_reps);
    auto [ceed_min, ceed_avg] = TimeIt(ceed_path, n_reps);
    std::cout << "Port power, complex E and B, per measurement:\n"
              << "  legacy linear form path: min " << legacy_min << " ms, avg "
              << legacy_avg << " ms\n"
              << "  libCEED functional:      min " << ceed_min << " ms, avg " << ceed_avg
              << " ms  (construction, once: " << t_constr << " ms)\n"
              << "  speedup (avg): " << legacy_avg / ceed_avg << "x\n";
  }
}

TEST_CASE("SurfaceFunctional Nonconformal Mesh", "[surfacefunctional][Serial][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();
  REQUIRE(pmesh.Nonconforming());

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Vector x0(3);
  x0(0) = 0.4;
  x0(1) = 0.6;
  x0(2) = -0.2;

  auto CheckPart = [](double v, double r)
  {
    if (std::abs(r) > 1.0e-12)
    {
      CHECK(v == Catch::Approx(r).epsilon(1.0e-10));
    }
    else
    {
      CHECK(std::abs(v) < 1.0e-10);
    }
  };

  // The exterior boundary attribute 1 (z = 0 in the refined half and beyond) and the
  // interior interface attribute 7 both contain a mix of conformal and nonconformal
  // (slave) faces after the partial refinement.
  for (int attr : {1, 6, 7})
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;

    // Interface dielectric (all variants).
    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA,
                      InterfaceDielectric::MS, InterfaceDielectric::SA})
    {
      CAPTURE(static_cast<int>(type));
      SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);
      REQUIRE(epr.IsValid());
      auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
      {
        switch (type)
        {
          case InterfaceDielectric::DEFAULT:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(
                E, mat_op, t_i, epsilon_i);
          case InterfaceDielectric::MA:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(E, mat_op, t_i,
                                                                         epsilon_i);
          case InterfaceDielectric::MS:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::MS>>(E, mat_op, t_i,
                                                                         epsilon_i);
          case InterfaceDielectric::SA:
            return std::make_unique<
                InterfaceDielectricCoefficient<InterfaceDielectric::SA>>(E, mat_op, t_i,
                                                                         epsilon_i);
        }
        return {};
      };
      auto legacy = MakeLegacy();
      const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
      const double val = epr.Eval(E);
      CAPTURE(ref, val);
      CheckPart(val, ref);
    }

    // Surface flux (all types, both two_sided settings).
    for (auto type : {SurfaceFlux::ELECTRIC, SurfaceFlux::MAGNETIC, SurfaceFlux::POWER})
    {
      for (bool two_sided : {false, true})
      {
        CAPTURE(static_cast<int>(type), two_sided);
        const GridFunction *E_use =
            (type == SurfaceFlux::ELECTRIC || type == SurfaceFlux::POWER) ? &E : nullptr;
        const GridFunction *B_use =
            (type == SurfaceFlux::MAGNETIC || type == SurfaceFlux::POWER) ? &B : nullptr;
        SurfaceFunctional flux(*mesh, marker, E_use ? &nd_fespace.Get() : nullptr,
                               B_use ? &rt_fespace.Get() : nullptr, mat_op, type, two_sided,
                               x0);
        REQUIRE(flux.IsValid());
        auto MakeLegacy =
            [&](const mfem::ParGridFunction *Er,
                const mfem::ParGridFunction *Br) -> std::unique_ptr<mfem::Coefficient>
        {
          switch (type)
          {
            case SurfaceFlux::ELECTRIC:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
                  Er, Br, mat_op, two_sided, x0);
            case SurfaceFlux::MAGNETIC:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>(
                  Er, Br, mat_op, two_sided, x0);
            case SurfaceFlux::POWER:
              return std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
                  Er, Br, mat_op, two_sided, x0);
          }
          return {};
        };
        auto legacy_re =
            MakeLegacy(E_use ? &E.Real() : nullptr, B_use ? &B.Real() : nullptr);
        std::complex<double> ref(RefSurfaceCoefficientIntegral(pmesh, *legacy_re, marker),
                                 0.0);
        auto legacy_im =
            MakeLegacy(E_use ? &E.Imag() : nullptr, B_use ? &B.Imag() : nullptr);
        const double ref_im = RefSurfaceCoefficientIntegral(pmesh, *legacy_im, marker);
        if (type == SurfaceFlux::POWER)
        {
          ref += ref_im;
        }
        else
        {
          ref.imag(ref_im);
        }
        const std::complex<double> val = flux.EvalFlux(E_use, B_use);
        CAPTURE(ref.real(), ref.imag(), val.real(), val.imag());
        CheckPart(val.real(), ref.real());
        CheckPart(val.imag(), ref.imag());
      }
    }
  }
}

TEST_CASE("SurfaceFunctional Nonconformal Parallel", "[surfacefunctional][Parallel][GPU]")
{
  // On parallel (possibly nonconformal) meshes, two-sided evaluation across process
  // boundaries uses FaceNbrFieldExchange to pull remote-side physical field values.
  // This test makes that processor-boundary behavior a committed regression: exterior
  // and interior interface surfaces must all assemble through libCEED and match the
  // legacy coefficient oracles.
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  const int order = 2;
  CAPTURE(elem_type);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeNCInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  // The legacy reference evaluates fields on remote sides of shared faces via face
  // neighbor data.
  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();
  B.Imag().ExchangeFaceNbrData();

  const double t_i = 2.0e-3, epsilon_i = 10.0;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

  for (int attr : {1, 6, 7})
  {
    CAPTURE(attr);
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    marker[attr - 1] = 1;
    auto CheckPart = [](double val, double ref)
    {
      if (std::abs(ref) > 1.0e-12)
      {
        CHECK(val == Catch::Approx(ref).epsilon(1.0e-10));
      }
      else
      {
        CHECK(std::abs(val) < 1.0e-10);
      }
    };

    for (auto type : {InterfaceDielectric::DEFAULT, InterfaceDielectric::MA})
    {
      CAPTURE(static_cast<int>(type));
      SurfaceFunctional epr(*mesh, marker, nd_fespace, mat_op, type, t_i, epsilon_i);

      bool valid = epr.IsValid();
      bool valid_and = valid, valid_or = valid;
      Mpi::GlobalAnd(1, &valid_and, comm);
      Mpi::GlobalOr(1, &valid_or, comm);
      REQUIRE(valid_and == valid_or);
      REQUIRE(valid);

      auto MakeLegacy = [&]() -> std::unique_ptr<mfem::Coefficient>
      {
        if (type == InterfaceDielectric::DEFAULT)
        {
          return std::make_unique<
              InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>(
              E, mat_op, t_i, epsilon_i);
        }
        return std::make_unique<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>(
            E, mat_op, t_i, epsilon_i);
      };
      auto legacy = MakeLegacy();
      const double ref = RefSurfaceCoefficientIntegral(pmesh, *legacy, marker);
      const double val = epr.Eval(E);
      CAPTURE(ref, val);
      CheckPart(val, ref);
    }

    mfem::Vector x0(3);
    x0 = 0.0;
    SurfaceFunctional power(*mesh, marker, &nd_fespace.Get(), &rt_fespace.Get(), mat_op,
                            SurfaceFlux::POWER, /*two_sided*/ true, x0);
    bool power_valid = power.IsValid();
    bool power_valid_and = power_valid, power_valid_or = power_valid;
    Mpi::GlobalAnd(1, &power_valid_and, comm);
    Mpi::GlobalOr(1, &power_valid_or, comm);
    REQUIRE(power_valid_and == power_valid_or);
    REQUIRE(power_valid);

    auto MakePowerLegacy = [&](const mfem::ParGridFunction &Er,
                               const mfem::ParGridFunction &Br)
    {
      auto coeff = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>(
          &Er, &Br, mat_op, /*two_sided*/ true, x0);
      return RefSurfaceCoefficientIntegral(pmesh, *coeff, marker);
    };
    const std::complex<double> ref(MakePowerLegacy(E.Real(), B.Real()) +
                                       MakePowerLegacy(E.Imag(), B.Imag()),
                                   0.0);
    const std::complex<double> val = power.EvalFlux(&E, &B);
    CAPTURE(ref.real(), val.real(), val.imag());
    CheckPart(val.real(), ref.real());
    CheckPart(val.imag(), 0.0);
  }
}

TEST_CASE("DomainFieldEvaluator", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, complex), B(rt_fespace, complex);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(1) * x(2) - 0.5;
                                          v(1) = std::sin(x(0)) - x(2);
                                          v(2) = std::cos(x(1)) + x(0) * x(0);
                                        });
    mfem::VectorFunctionCoefficient fbi(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::cos(x(2)) - 0.2;
                                          v(1) = x(0) * x(2) + 0.1;
                                          v(2) = std::sin(x(1)) - x(0);
                                        });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  // Interpolatory L2 output spaces (legacy ProjectCoefficient on these spaces evaluates
  // the coefficient at the nodal points, which is exactly the libCEED evaluator
  // semantics, at any order).
  mfem::L2_FECollection viz_fec(order, 3);
  mfem::ParFiniteElementSpace viz_scalar(&pmesh, &viz_fec), viz_vector(&pmesh, &viz_fec, 3);

  const double scaling = 2.5;
  auto CheckField = [](const mfem::ParGridFunction &val, const mfem::ParGridFunction &ref)
  {
    // HostRead to sync from device (the evaluator fills on device).
    const double *v = val.HostRead();
    const double *r = ref.HostRead();
    double max_diff = 0.0, max_ref = 0.0;
    for (int i = 0; i < ref.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(v[i] - r[i]));
      max_ref = std::max(max_ref, std::abs(r[i]));
    }
    CAPTURE(max_diff, max_ref);
    CHECK(max_diff <= 1.0e-11 * std::max(max_ref, 1.0));
  };

  SECTION("Electric energy density")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::ENERGY_E, *mesh, mat_op,
                              E.ParFESpace(), nullptr, viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(&E, nullptr, val);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Magnetic energy density")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::ENERGY_M, *mesh, mat_op, nullptr,
                              B.ParFESpace(), viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(nullptr, &B, val);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Poynting vector")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::POYNTING, *mesh, mat_op,
                              E.ParFESpace(), B.ParFESpace(), viz_vector, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_vector), ref(&viz_vector);
    eval.Eval(&E, &B, val);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }
}

TEST_CASE("DomainFieldEvaluator 2D", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto order = GENERATE(1, 2);
  auto complex = GENERATE(false, true);
  CAPTURE(elem_type, order, complex);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(3, 2, elem_type, false, 1.2, 0.7));
  smesh->EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh->GetNE());
  auto pmesh_ptr = std::make_unique<mfem::ParMesh>(comm, *smesh);
  Mesh mesh(std::move(pmesh_ptr));
  auto &pmesh = mesh.Get();

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s[0] = 2.0;
  material.epsilon_r.s[1] = 3.0;
  material.epsilon_r.s[2] = 4.0;
  material.mu_r.s[0] = 1.0;
  material.mu_r.s[1] = 1.5;
  material.mu_r.s[2] = 2.0;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, mesh);

  mfem::ND_FECollection nd_fec(order, 2);
  mfem::L2_FECollection l2_fec(order, 2);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), l2_fespace(mesh, &l2_fec);

  GridFunction E(nd_fespace, complex), B(l2_fespace, complex);
  mfem::VectorFunctionCoefficient fer(2,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + 0.2 * x(0);
                                        v(1) = std::cos(x(0)) + x(1) * x(1);
                                      });
  mfem::FunctionCoefficient fbr(
      [](const mfem::Vector &x) { return 0.3 + x(0) - 0.7 * x(1); });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  if (complex)
  {
    mfem::VectorFunctionCoefficient fei(2,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = x(0) * x(1) - 0.1;
                                          v(1) = std::sin(x(0) + x(1));
                                        });
    mfem::FunctionCoefficient fbi(
        [](const mfem::Vector &x) { return std::cos(x(0)) + 0.4 * x(1); });
    E.Imag().ProjectCoefficient(fei);
    B.Imag().ProjectCoefficient(fbi);
  }

  mfem::L2_FECollection viz_fec(order, 2);
  mfem::ParFiniteElementSpace viz_scalar(&pmesh, &viz_fec), viz_vector(&pmesh, &viz_fec, 2);

  const double scaling = 1.7;
  auto CheckField = [](const mfem::ParGridFunction &val, const mfem::ParGridFunction &ref)
  {
    const double *v = val.HostRead();
    const double *r = ref.HostRead();
    double max_diff = 0.0, max_ref = 0.0;
    for (int i = 0; i < ref.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(v[i] - r[i]));
      max_ref = std::max(max_ref, std::abs(r[i]));
    }
    CAPTURE(max_diff, max_ref);
    CHECK(max_diff <= 1.0e-11 * std::max(max_ref, 1.0));
  };

  SECTION("Electric energy density")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::ENERGY_E, mesh, mat_op,
                              E.ParFESpace(), nullptr, viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(&E, nullptr, val);
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Magnetic energy density")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::ENERGY_M, mesh, mat_op, nullptr,
                              B.ParFESpace(), viz_scalar, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_scalar), ref(&viz_scalar);
    eval.Eval(nullptr, &B, val);
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }

  SECTION("Poynting vector")
  {
    DomainFieldEvaluator eval(DomainFieldEvaluator::Kind::POYNTING, mesh, mat_op,
                              E.ParFESpace(), B.ParFESpace(), viz_vector, scaling);
    REQUIRE(eval.IsValid());
    mfem::ParGridFunction val(&viz_vector), ref(&viz_vector);
    eval.Eval(&E, &B, val);
    PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
    ref.ProjectCoefficient(legacy);
    CheckField(val, ref);
  }
}

TEST_CASE("SurfaceFunctional FarField", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};  // Isotropic (far-field requirement)
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  mfem::VectorFunctionCoefficient fbi(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::cos(x(2)) - 0.2;
                                        v(1) = x(0) * x(2) + 0.1;
                                        v(2) = std::sin(x(1)) - x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);
  B.Real().ProjectCoefficient(fbr);
  B.Imag().ProjectCoefficient(fbi);

  // Observation directions and (complex) frequency.
  std::vector<std::array<double, 3>> r_naughts;
  for (auto [theta, phi] : {std::pair{0.3, 0.7}, {1.2, 2.1}, {2.4, 4.5}})
  {
    r_naughts.push_back({std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
                         std::cos(theta)});
  }
  const double omega_re = 2.7, omega_im = 0.15;

  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 0;
  marker[1 - 1] = 1;  // Exterior boundary (z = 0)

  // Legacy reference following GetFarFieldrE.
  std::vector<std::array<double, 3>> integrals_r(r_naughts.size()),
      integrals_i(r_naughts.size());
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(pmesh).GetBdrElementTransformation(i);
    const auto *fe = nd_fespace.Get().GetBE(i);
    const auto *ir =
        &mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));
    AddStrattonChuIntegrandAtElement(E, B, mat_op, omega_re, omega_im, r_naughts, *T, *ir,
                                     integrals_r, integrals_i);
  }
  Mpi::GlobalSum(3 * r_naughts.size(), integrals_r.data()->data(), comm);
  Mpi::GlobalSum(3 * r_naughts.size(), integrals_i.data()->data(), comm);

  SurfaceFunctional farfield(*mesh, marker, nd_fespace, rt_fespace, mat_op, r_naughts);
  REQUIRE(farfield.IsValid());
  auto result = farfield.EvalFarField(E, B, {omega_re, omega_im});
  REQUIRE(result.size() == r_naughts.size());
  for (std::size_t d = 0; d < r_naughts.size(); d++)
  {
    const auto &r = r_naughts[d];
    const auto &Ir = integrals_r[d];
    const auto &Ii = integrals_i[d];
    const std::array<double, 3> cr = {r[1] * Ir[2] - r[2] * Ir[1],
                                      r[2] * Ir[0] - r[0] * Ir[2],
                                      r[0] * Ir[1] - r[1] * Ir[0]};
    const std::array<double, 3> ci = {r[1] * Ii[2] - r[2] * Ii[1],
                                      r[2] * Ii[0] - r[0] * Ii[2],
                                      r[0] * Ii[1] - r[1] * Ii[0]};
    for (int c = 0; c < 3; c++)
    {
      CAPTURE(d, c, cr[c], ci[c], result[d][c].real(), result[d][c].imag());
      CHECK(result[d][c].real() == Catch::Approx(cr[c]).epsilon(1.0e-10).margin(1.0e-14));
      CHECK(result[d][c].imag() == Catch::Approx(ci[c]).epsilon(1.0e-10).margin(1.0e-14));
    }
  }

  // The AtPoints-specialized far-field path must agree with the mapped integration-rule
  // path independently for every observation direction.
  const char *old_disable = std::getenv("PALACE_SURFACE_DISABLE_ATPOINTS");
  setenv("PALACE_SURFACE_DISABLE_ATPOINTS", "1", 1);
  SurfaceFunctional mapped_farfield(*mesh, marker, nd_fespace, rt_fespace, mat_op,
                                    r_naughts);
  REQUIRE(mapped_farfield.IsValid());
  auto mapped_result = mapped_farfield.EvalFarField(E, B, {omega_re, omega_im});
  if (old_disable)
  {
    setenv("PALACE_SURFACE_DISABLE_ATPOINTS", old_disable, 1);
  }
  else
  {
    unsetenv("PALACE_SURFACE_DISABLE_ATPOINTS");
  }
  for (std::size_t d = 0; d < r_naughts.size(); d++)
  {
    for (int c = 0; c < 3; c++)
    {
      CAPTURE(d, c, result[d][c].real(), result[d][c].imag(),
              mapped_result[d][c].real(), mapped_result[d][c].imag());
      CHECK(result[d][c].real() ==
            Catch::Approx(mapped_result[d][c].real()).epsilon(1.0e-10).margin(1.0e-14));
      CHECK(result[d][c].imag() ==
            Catch::Approx(mapped_result[d][c].imag()).epsilon(1.0e-10).margin(1.0e-14));
    }
  }

  // Changing the frequency must reassemble and still agree (different result).
  auto result2 = farfield.EvalFarField(E, B, {2.0 * omega_re, 0.0});
  CHECK(std::abs(result2[0][0] - result[0][0]) > 0.0);
}

#if defined(MFEM_USE_GSLIB)
TEST_CASE("InterpolationOperator Ceed Probes", "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  GridFunction E(nd_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);

  // Probe points (element interiors, in both material regions; points on element
  // borders are avoided since interpolated values of H(curl) fields are multi-valued
  // there and the GSLIB reference resolves the donor element rank-dependently).
  std::map<int, config::ProbeData> probes;
  const std::array<std::array<double, 3>, 3> pts = {
      {{0.21, 0.37, 0.23}, {0.74, 0.52, 0.81}, {0.53, 0.48, 0.52}}};
  for (std::size_t i = 0; i < pts.size(); i++)
  {
    config::ProbeData data;
    data.center = pts[i];
    probes.emplace(static_cast<int>(i + 1), data);
  }
  Units units(1.0, 1.0);
  InterpolationOperator interp(probes, units, nd_fespace);

  // Reference: GSLIB interpolation at the same points (byVDIM).
  const int npts = static_cast<int>(pts.size());
  mfem::Vector xyz(npts * 3), vals_r(npts * 3), vals_i(npts * 3);
  for (int i = 0; i < npts; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      xyz(d * npts + i) = pts[i][d];
    }
  }
  fem::InterpolateFunction(xyz, E.Real(), vals_r, mfem::Ordering::byNODES);
  fem::InterpolateFunction(xyz, E.Imag(), vals_i, mfem::Ordering::byNODES);

  auto vals = interp.ProbeField(E);
  REQUIRE(static_cast<int>(vals.size()) == npts * 3);
  for (int i = 0; i < npts; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      // ProbeField returns byVDIM; the InterpolateFunction reference returns with the
      // requested (byNODES) ordering.
      const auto val = vals[3 * i + d];
      const double ref_r = vals_r(d * npts + i), ref_i = vals_i(d * npts + i);
      CAPTURE(i, d, ref_r, ref_i);
      CHECK(val.real() == Catch::Approx(ref_r).epsilon(1.0e-9).margin(1.0e-12));
      CHECK(val.imag() == Catch::Approx(ref_i).epsilon(1.0e-9).margin(1.0e-12));
    }
  }
}
#endif

TEST_CASE("SurfaceFunctional Boundary Viz Fields",
          "[surfacefunctional][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s[0] = 11.7;
  dielectric.epsilon_r.s[1] = 11.7;
  dielectric.epsilon_r.s[2] = 11.7;
  dielectric.mu_r.s[0] = 1.4;
  dielectric.mu_r.s[1] = 1.4;
  dielectric.mu_r.s[2] = 1.4;
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, false), B(rt_fespace, false);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fbr(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) - 0.3 * x(2);
                                        v(1) = std::sin(x(2)) + 0.5;
                                        v(2) = std::cos(x(0)) - x(1) * x(2);
                                      });
  E.Real().ProjectCoefficient(fer);
  B.Real().ProjectCoefficient(fbr);
  E.Real().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();

  const int lod = order;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 1;

  auto TestKind = [&](SurfaceFunctional::Kind kind, const mfem::ParFiniteElementSpace &fes,
                      const mfem::ParGridFunction &U)
  {
    SurfaceFunctional viz(kind, *mesh, marker, fes, lod);

    // The validity decision must be identical on all ranks; interior surfaces split
    // across processes fall back (consistently).
    bool valid = viz.IsValid();
    bool valid_and = valid, valid_or = valid;
    Mpi::GlobalAnd(1, &valid_and, comm);
    Mpi::GlobalOr(1, &valid_or, comm);
    REQUIRE(valid_and == valid_or);
    if (!valid)
    {
      return;
    }

    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    viz.EvalBuffer(U, buffer);
    const double *buf = buffer.HostRead();
    const auto &bases = viz.BufferBases();
    const int component_stride = viz.BufferSize() / 3;

    BdrFieldVectorCoefficient legacy(U);
    mfem::Vector V(3);
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const auto &RefG =
          *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
      auto *T = pmesh.GetBdrElementTransformation(i);
      for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
      {
        const auto &ip = RefG.RefPts.IntPoint(j);
        T->SetIntPoint(&ip);
        legacy.Eval(V, *T, ip);
        for (int c = 0; c < 3; c++)
        {
          const double val = buf[bases[i] + j + c * component_stride];
          CAPTURE(i, j, c, V(c), val);
          CHECK(val == Catch::Approx(V(c)).epsilon(1.0e-10).margin(1.0e-13));
        }
      }
    }
  };

  TestKind(SurfaceFunctional::Kind::BDR_FIELD_E, nd_fespace, E.Real());
  TestKind(SurfaceFunctional::Kind::BDR_FIELD_B, rt_fespace, B.Real());

  // Material-dependent boundary visualization kinds (surface charge, surface current,
  // boundary energy densities) against the corresponding legacy coefficients.
  const double scaling = 1.7;
  auto TestKindCoeff = [&](SurfaceFunctional::Kind kind,
                           const mfem::ParFiniteElementSpace &fes,
                           mfem::Coefficient *legacy_s, mfem::VectorCoefficient *legacy_v)
  {
    const int nc = SurfaceFunctional::BufferNumComp(kind);
    SurfaceFunctional viz(kind, *mesh, marker, fes, mat_op, lod, scaling);
    bool valid = viz.IsValid();
    bool valid_and = valid, valid_or = valid;
    Mpi::GlobalAnd(1, &valid_and, comm);
    Mpi::GlobalOr(1, &valid_or, comm);
    REQUIRE(valid_and == valid_or);
    if (!valid)
    {
      return;
    }
    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    viz.EvalBuffer(kind == SurfaceFunctional::Kind::BDR_CURRENT_J ||
                           kind == SurfaceFunctional::Kind::BDR_ENERGY_M
                       ? B.Real()
                       : E.Real(),
                   buffer);
    const double *buf = buffer.HostRead();
    const auto &bases = viz.BufferBases();
    const int component_stride = viz.BufferSize() / nc;
    mfem::Vector V(3);
    for (int i = 0; i < pmesh.GetNBE(); i++)
    {
      const auto &RefG =
          *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
      auto *T = pmesh.GetBdrElementTransformation(i);
      for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
      {
        const auto &ip = RefG.RefPts.IntPoint(j);
        T->SetIntPoint(&ip);
        if (nc == 1)
        {
          const double ref = legacy_s->Eval(*T, ip);
          const double val = buf[bases[i] + j];
          CAPTURE(i, j, ref, val);
          CHECK(val == Catch::Approx(ref).epsilon(1.0e-10).margin(1.0e-13));
        }
        else
        {
          legacy_v->Eval(V, *T, ip);
          for (int c = 0; c < 3; c++)
          {
            const double val = buf[bases[i] + j + c * component_stride];
            CAPTURE(i, j, c, V(c), val);
            CHECK(val == Catch::Approx(V(c)).epsilon(1.0e-10).margin(1.0e-13));
          }
        }
      }
    }
  };
  {
    BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> q_legacy(
        &E.Real(), nullptr, mat_op, true, mfem::Vector(), scaling);
    TestKindCoeff(SurfaceFunctional::Kind::BDR_FLUX_Q, nd_fespace, &q_legacy, nullptr);
  }
  {
    BdrSurfaceCurrentVectorCoefficient j_legacy(B.Real(), mat_op, scaling);
    TestKindCoeff(SurfaceFunctional::Kind::BDR_CURRENT_J, rt_fespace, nullptr, &j_legacy);
  }
  {
    EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> ue_legacy(E, mat_op, scaling);
    TestKindCoeff(SurfaceFunctional::Kind::BDR_ENERGY_E, nd_fespace, &ue_legacy, nullptr);
  }
  {
    EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> um_legacy(B, mat_op, scaling);
    TestKindCoeff(SurfaceFunctional::Kind::BDR_ENERGY_M, rt_fespace, &um_legacy, nullptr);
  }
}

TEST_CASE("FaceNbrFieldExchange", "[surfacefunctional][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal);

  auto mesh = nonconformal ? MakeNCInterfaceMesh(comm, elem_type)
                           : MakeInterfaceMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  // Project non-trivial smooth fields. The reference values are computed from the same
  // projected grid functions through the legacy mfem face neighbor dof exchange, so the
  // comparison is exact up to roundoff (not projection error).
  mfem::ParGridFunction E(&nd_fespace.Get()), B(&rt_fespace.Get());
  mfem::VectorFunctionCoefficient fe(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1)) + x(2) * x(2);
                                       v(1) = std::cos(x(2)) + x(0);
                                       v(2) = x(0) * x(1) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient fb(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) * x(2) - 0.5;
                                       v(1) = std::sin(x(0)) - x(2);
                                       v(2) = std::cos(x(1)) + x(0) * x(0);
                                     });
  E.ProjectCoefficient(fe);
  B.ProjectCoefficient(fb);

  // Request E (slot 0) at a few reference points of every ghost element, and both E
  // and B (slot 1) for every other ghost element (exercising the value layouts). The
  // points are valid reference coordinates for both tetrahedra and hexahedra.
  const int num_ghost = pmesh.GetNFaceNeighborElements();
  std::vector<FaceNbrFieldExchange::Request> requests;
  for (int fn = 0; fn < num_ghost; fn++)
  {
    auto &req = requests.emplace_back();
    req.face_nbr_elem = fn;
    req.source_mask = (fn % 2 == 0) ? 0b01u : 0b11u;
    req.point_key = {static_cast<long long>(elem_type), static_cast<long long>(order),
                     4};
    req.pts.resize(4);
    req.pts[0].Set3(0.1, 0.2, 0.3);
    req.pts[1].Set3(0.25, 0.25, 0.25);
    req.pts[2].Set3(0.05, 0.1, 0.7);
    req.pts[3].Set3(0.3, 0.05, 0.05);
  }
  FaceNbrFieldExchange exchange(
      *mesh, {&nd_fespace.Get(), &rt_fespace.Get(), nullptr, nullptr}, requests);
  exchange.Exchange({&E, &B, nullptr, nullptr});

  // Reference: evaluate the ghost elements through the legacy dof exchange.
  E.ExchangeFaceNbrData();
  B.ExchangeFaceNbrData();
  std::vector<double> vals(exchange.Imported().HostRead(),
                           exchange.Imported().HostRead() + exchange.ImportSize());
  mfem::Vector ref(3);
  int num_checked = 0;
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    for (int s = 0; s < 2; s++)
    {
      const int offset = exchange.ImportOffset(static_cast<int>(r), s);
      if (!(req.source_mask & (1u << s)))
      {
        CHECK(offset < 0);
        continue;
      }
      const auto &U = (s == 0) ? E : B;
      for (std::size_t j = 0; j < req.pts.size(); j++)
      {
        U.GetVectorValue(pmesh.GetNE() + req.face_nbr_elem, req.pts[j], ref);
        for (int c = 0; c < 3; c++)
        {
          CAPTURE(r, s, j, c);
          CHECK(vals[offset + 3 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
          num_checked++;
        }
      }
    }
  }
  // With more than one process, the partition must produce at least one ghost element
  // somewhere (the exchange is the point of the test).
  int num_global = num_checked;
  Mpi::GlobalSum(1, &num_global, comm);
  CHECK((Mpi::Size(comm) == 1 || num_global > 0));

  // The field inputs are re-pointed at the sources on each call: scaling the field
  // scales the exchanged values.
  E *= 2.0;
  exchange.Exchange({&E, &B, nullptr, nullptr});
  const double *vals2 = exchange.Imported().HostRead();
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const int offset = exchange.ImportOffset(static_cast<int>(r), 0);
    for (std::size_t j = 0; j < 3 * requests[r].pts.size(); j++)
    {
      CHECK(vals2[offset + j] == Catch::Approx(2.0 * vals[offset + j]).margin(1.0e-12));
    }
  }
}

}  // namespace palace
