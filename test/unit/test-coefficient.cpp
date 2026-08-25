// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

using namespace palace;

TEST_CASE("BdrFieldCoefficient interior average", "[coefficient][Serial]")
{
  // 2x1 quad mesh: two elements sharing an interior face at x=1.
  mfem::Mesh mesh =
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);

  // Add interior boundary at x=1 (the shared face between elements).
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)  // Interior face.
    {
      auto *face_elem = mesh.GetFace(f)->Duplicate(&mesh);
      mesh.AddBdrElement(face_elem);
      mesh.SetBdrAttribute(mesh.GetNBE() - 1, 99);
    }
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  auto pmesh = mfem::ParMesh(MPI_COMM_WORLD, mesh);

  // H1 space, order 1.
  mfem::H1_FECollection fec(1, 2);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec);
  mfem::ParGridFunction u(&fes);

  // Set u = x. At interior face (x=1): elem0 centroid ~0.5, elem1 centroid ~1.5.
  // Proper average = 1.0. Bug (comma op) would give 0.5 * (value from elem1).
  mfem::FunctionCoefficient x_coord([](const mfem::Vector &x) { return x(0); });
  u.ProjectCoefficient(x_coord);

  BdrFieldCoefficient coeff(u);

  // Evaluate on interior boundary (attr 99).
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (pmesh.GetBdrAttribute(be) == 99)
    {
      auto *T = pmesh.GetBdrElementTransformation(be);
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

TEST_CASE("InterfaceDielectricCoefficient recovered normal flux", "[coefficient][Serial]")
{
  // Two materials separated by an internal boundary at x=1. Element 1 is substrate
  // (epsilon=4), and element 2 is vacuum.
  mfem::Mesh serial_mesh =
      mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::QUADRILATERAL, false, 2.0, 1.0);
  serial_mesh.SetAttribute(0, 1);
  serial_mesh.SetAttribute(1, 2);
  for (int f = 0; f < serial_mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    serial_mesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)
    {
      serial_mesh.AddBdrElement(serial_mesh.GetFace(f)->Duplicate(&serial_mesh));
      serial_mesh.SetBdrAttribute(serial_mesh.GetNBE() - 1, 99);
    }
  }
  serial_mesh.FinalizeTopology();
  serial_mesh.Finalize();

  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  Mesh mesh(std::move(par_mesh));
  config::MaterialData substrate, vacuum;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  vacuum.attributes = {2};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({substrate, vacuum}, periodic, ProblemType::ELECTROSTATIC, mesh);

  mfem::ND_FECollection nd_fec(1, mesh.Dimension());
  mfem::RT_FECollection rt_fec(0, mesh.Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpace rt_fespace(mesh, &rt_fec);
  GridFunction E(nd_fespace, true), D(rt_fespace, true);

  mfem::VectorFunctionCoefficient Er(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 3.0;
                                       v[1] = 4.0;
                                     });
  mfem::VectorFunctionCoefficient Ei(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 5.0;
                                       v[1] = 6.0;
                                     });
  mfem::VectorFunctionCoefficient Dr(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 7.0;
                                       v[1] = 8.0;
                                     });
  mfem::VectorFunctionCoefficient Di(2,
                                     [](const mfem::Vector &, mfem::Vector &v)
                                     {
                                       v[0] = 9.0;
                                       v[1] = 10.0;
                                     });
  E.Real().ProjectCoefficient(Er);
  E.Imag().ProjectCoefficient(Ei);
  D.Real().ProjectCoefficient(Dr);
  D.Imag().ProjectCoefficient(Di);
  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  D.Real().ExchangeFaceNbrData();
  D.Imag().ExchangeFaceNbrData();

  mfem::ElementTransformation *T = nullptr;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (mesh.Get().GetBdrAttribute(be) == 99)
    {
      T = mesh.Get().GetBdrElementTransformation(be);
      break;
    }
  }
  REQUIRE(T);
  mfem::IntegrationPoint ip;
  ip.Set1w(0.5, 1.0);
  T->SetIntPoint(&ip);

  constexpr double thickness = 2.0;
  constexpr double interface_epsilon = 5.0;
  constexpr double raw_normal_squared = 3.0 * 3.0 + 5.0 * 5.0;
  constexpr double recovered_normal_squared = 7.0 * 7.0 + 9.0 * 9.0;
  constexpr double tangent_squared = 4.0 * 4.0 + 6.0 * 6.0;

  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_recovered(
      E, mat_op, thickness, interface_epsilon, &D);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_normal(
      E, mat_op, thickness, interface_epsilon, &D, InterfaceDielectricComponent::NORMAL);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_tangential(
      E, mat_op, thickness, interface_epsilon, &D,
      InterfaceDielectricComponent::TANGENTIAL);
  CHECK_THAT(ma_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * raw_normal_squared, 1.0e-12));
  CHECK_THAT(ma_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * recovered_normal_squared, 1.0e-12));
  CHECK_THAT(ma_normal.Eval(*T, ip),
             Catch::Matchers::WithinAbs(ma_recovered.Eval(*T, ip), 1.0e-12));
  CHECK_THAT(ma_tangential.Eval(*T, ip), Catch::Matchers::WithinAbs(0.0, 1.0e-12));

  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_recovered(
      E, mat_op, thickness, interface_epsilon, &D);
  CHECK_THAT(ms_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness / interface_epsilon * 4.0 * 4.0 *
                                            raw_normal_squared,
                                        1.0e-12));
  CHECK_THAT(ms_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(
                 0.5 * thickness / interface_epsilon * recovered_normal_squared, 1.0e-12));

  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_raw(E, mat_op, thickness,
                                                                 interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_recovered(
      E, mat_op, thickness, interface_epsilon, &D);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_normal(
      E, mat_op, thickness, interface_epsilon, &D, InterfaceDielectricComponent::NORMAL);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_tangential(
      E, mat_op, thickness, interface_epsilon, &D,
      InterfaceDielectricComponent::TANGENTIAL);
  CHECK_THAT(sa_raw.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness *
                                            (interface_epsilon * tangent_squared +
                                             raw_normal_squared / interface_epsilon),
                                        1.0e-12));
  CHECK_THAT(sa_recovered.Eval(*T, ip),
             Catch::Matchers::WithinAbs(0.5 * thickness *
                                            (interface_epsilon * tangent_squared +
                                             recovered_normal_squared / interface_epsilon),
                                        1.0e-12));
  CHECK_THAT(sa_normal.Eval(*T, ip) + sa_tangential.Eval(*T, ip),
             Catch::Matchers::WithinAbs(sa_recovered.Eval(*T, ip), 1.0e-12));

  // The real-valued amplitudes used by batched response-matrix assembly reproduce the
  // ordinary quadratic coefficients exactly.
  GridFunction E_batch(nd_fespace), D_batch(rt_fespace);
  E_batch.Real().ProjectCoefficient(Er);
  D_batch.Real().ProjectCoefficient(Dr);
  E_batch.Real().ExchangeFaceNbrData();
  D_batch.Real().ExchangeFaceNbrData();
  mfem::Vector field;
  auto CheckBatched = [&](auto &coefficient)
  {
    coefficient.EvalEnergyField(*T, ip, field);
    CHECK_THAT(0.5 * (field * field),
               Catch::Matchers::WithinAbs(coefficient.Eval(*T, ip), 1.0e-12));
  };
  InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT> default_batch(
      E_batch, mat_op, thickness, interface_epsilon);
  InterfaceDielectricCoefficient<InterfaceDielectric::MA> ma_batch(
      E_batch, mat_op, thickness, interface_epsilon, &D_batch);
  InterfaceDielectricCoefficient<InterfaceDielectric::MS> ms_batch(
      E_batch, mat_op, thickness, interface_epsilon, &D_batch);
  InterfaceDielectricCoefficient<InterfaceDielectric::SA> sa_batch(
      E_batch, mat_op, thickness, interface_epsilon, &D_batch);
  CheckBatched(default_batch);
  CheckBatched(ma_batch);
  CheckBatched(ms_batch);
  CheckBatched(sa_batch);
}
