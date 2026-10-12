// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <functional>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/materialfield.hpp"
#include "fem/mesh.hpp"
#include "linalg/hypre.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Build a per-attribute scalar coefficient using the mesh's local libCEED attribute
// mapping. `value(global_attr)` returns the coefficient for a global domain attribute.
MaterialPropertyCoefficient MakeAttributeCoefficient(
    const Mesh &mesh, const std::function<double(int)> &value)
{
  const auto &loc_attr = mesh.GetCeedAttributes();
  mfem::Array<int> attr_mat(loc_attr.size());
  mfem::DenseTensor mat_coeff;
  mat_coeff.SetSize(1, 1, loc_attr.size());
  for (const auto &[global_attr, ceed_attr] : loc_attr)
  {
    attr_mat[ceed_attr - 1] = ceed_attr - 1;
    mat_coeff(0, 0, ceed_attr - 1) = value(global_attr);
  }
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

// Relative error of the operators A and B applied to a deterministic test vector.
double RelativeOperatorError(const Operator &A, const Operator &B, int size, MPI_Comm comm)
{
  Vector x(size), y1(size), y2(size);
  for (int i = 0; i < size; i++)
  {
    x[i] = 0.3 + std::sin(1.7 * (i + 1)) + 0.1 * (i % 7);
  }
  A.Mult(x, y1);
  B.Mult(x, y2);
  y1 -= y2;
  double err = y1.Norml2(), ref = y2.Norml2();
  Mpi::GlobalSum(1, &err, comm);
  Mpi::GlobalSum(1, &ref, comm);
  return (ref > 0.0) ? err / ref : err;
}

// Initialize the libCEED integration order, normally configured by
// IoData::CheckConfiguration.
void SetIntegrationOrder(int p)
{
  fem::DefaultIntegrationOrder::p_trial = p;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
}

}  // namespace

// A grid function that is exactly constant on every element must reproduce the
// attribute-based constant material assembly exactly. Every element is assigned its own
// attribute and an arbitrary value, exercising the full attribute/material mapping.
TEST_CASE("MaterialField elementwise-constant H1 mass equivalence",
          "[materialfield][Serial][Parallel]")
{
  constexpr int nx = 6, ny = 4;
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, false, 2.0, 1.0));
  for (int e = 0; e < serial_mesh->GetNE(); e++)
  {
    serial_mesh->SetAttribute(e, 1 + (e % 11));
  }
  serial_mesh->SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));

  constexpr int order = 1;
  SetIntegrationOrder(order);
  auto fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
  FiniteElementSpace fespace(mesh, fec.get());

  auto value = [](int attr) { return 0.75 + 0.37 * attr - 0.01 * attr * attr; };
  MaterialPropertyCoefficient Q = MakeAttributeCoefficient(mesh, value);

  BilinearForm a_ref(fespace);
  a_ref.AddDomainIntegrator<MassIntegrator>(Q);
  auto A_ref = a_ref.FullAssemble(false);

  // Elementwise-constant (L2, order 0) material field matching the per-element values.
  mfem::L2_FECollection l2_fec(0, mesh.Dimension());
  mfem::ParFiniteElementSpace l2_fes(&mesh.Get(), &l2_fec);
  mfem::ParGridFunction gf(&l2_fes);
  mfem::Array<int> dofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    l2_fes.GetElementDofs(e, dofs);
    gf(dofs[0]) = value(mesh.Get().GetAttribute(e));
  }
  MaterialField field(gf);

  BilinearForm a_field(fespace);
  a_field.AddDomainIntegrator<MassIntegrator>().SetMaterialField(&field);
  auto A_field = a_field.FullAssemble(false);

  CHECK(RelativeOperatorError(*A_ref, *A_field, fespace.GetVSize(), Mpi::World()) < 1.0e-13);
}

// Verify that the field data cached for a smooth higher-order grid function matches direct
// evaluation at the libCEED quadrature points for every element geometry.
TEST_CASE("MaterialField smooth quadrature sampling", "[materialfield][Serial][Parallel]")
{
  const auto type = GENERATE(mfem::Element::QUADRILATERAL, mfem::Element::TRIANGLE);
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(5, 3, type, false, 1.5, 1.0));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh mesh(std::move(par_mesh));

  constexpr int order = 2;
  SetIntegrationOrder(order);
  auto fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
  FiniteElementSpace fespace(mesh, fec.get());

  mfem::ParGridFunction gf(&fespace.Get());
  mfem::FunctionCoefficient coeff(
      [](const mfem::Vector &x) { return 1.25 + 0.5 * x(0) - 0.75 * x(1) * x(1); });
  gf.ProjectCoefficient(coeff);
  MaterialField field(gf);

  double max_err = 0.0;
  for (Ceed ceed : ceed::internal::GetCeedObjects())
  {
    for (const auto &[geom, geom_data] : mesh.GetCeedGeomFactorData(ceed))
    {
      if (mfem::Geometry::Dimension[geom] != mesh.Dimension())
      {
        continue;
      }
      const CeedMaterialFieldData &data = field.GetData(ceed, geom, geom_data.indices);
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);
      const int num_qpts = ir.GetNPoints();
      const CeedScalar *values;
      PalaceCeedCall(ceed, CeedVectorGetArrayRead(data.data, CEED_MEM_HOST, &values));
      for (std::size_t e = 0; e < geom_data.indices.size(); e++)
      {
        for (int q = 0; q < num_qpts; q++)
        {
          const double expected = gf.GetValue(geom_data.indices[e], ir.IntPoint(q));
          max_err = std::max(max_err, std::abs(values[e * num_qpts + q] - expected));
        }
      }
      PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(data.data, &values));
    }
  }
  CHECK(max_err == 0.0);
}

}  // namespace palace
