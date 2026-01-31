// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

//
// Unit tests and benchmarks for IntegrateFunction and IntegrateFunctionLocal.
//
// These functions replace MFEM's LinearForm + BoundaryLFIntegrator pattern for computing
// surface integrals. The key advantage is that they take a callback for quadrature order,
// allowing the caller to specify the order without relying on global state.
//

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"

namespace palace
{
namespace
{

// Simple constant coefficient for testing.
// This is thread-safe as it has no mutable state.
class ConstantCoefficient : public mfem::Coefficient
{
private:
  double value;

public:
  ConstantCoefficient(double v) : value(v) {}
  double Eval(mfem::ElementTransformation &, const mfem::IntegrationPoint &) override
  {
    return value;
  }
};

// Coordinate-dependent coefficient for testing: f(x,y,z) = x + 2y + 3z.
// This is thread-safe as it has no mutable state.
class LinearCoefficient : public mfem::Coefficient
{
public:
  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector x;
    T.Transform(ip, x);
    return x[0] + 2.0 * x[1] + (x.Size() > 2 ? 3.0 * x[2] : 0.0);
  }
};

// Helper to compute integral using MFEM's LinearForm with BoundaryLFIntegrator.
double IntegrateWithLinearForm(mfem::ParFiniteElementSpace &fespace,
                               const mfem::Array<int> &marker, mfem::Coefficient &Q,
                               int q_order)
{
  // Set global quadrature order for MFEM integrator.
  fem::DefaultIntegrationOrder::p_trial = q_order;
  fem::DefaultIntegrationOrder::q_order_jac = false;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  mfem::LinearForm lf(&fespace);
  lf.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(Q),
                           const_cast<mfem::Array<int> &>(marker));
  lf.UseFastAssembly(false);
  lf.UseDevice(false);
  lf.Assemble();
  lf.UseDevice(true);

  // Sum all entries (this is what the original GetLocalSurfaceIntegral did).
  double sum = 0.0;
  for (int i = 0; i < lf.Size(); i++)
  {
    sum += lf[i];
  }

  // Global sum across MPI ranks.
  double global_sum = sum;
  MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, fespace.GetComm());
  return global_sum;
}

// Create a simple 2D mesh for testing.
std::unique_ptr<mfem::ParMesh> CreateTestMesh2D(MPI_Comm comm)
{
  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL, true, 1.0, 1.0));
  return std::make_unique<mfem::ParMesh>(comm, *mesh);
}

// Create a simple 3D mesh for testing.
std::unique_ptr<mfem::ParMesh> CreateTestMesh3D(MPI_Comm comm)
{
  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  return std::make_unique<mfem::ParMesh>(comm, *mesh);
}

TEST_CASE("IntegrateFunction", "[integrate][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  SECTION("2D boundary integral with constant coefficient")
  {
    auto mesh = CreateTestMesh2D(comm);
    mfem::H1_FECollection fec(1, mesh->Dimension());
    mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);

    // Mark all boundary elements.
    int bdr_attr_max = mesh->bdr_attributes.Max();
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;

    ConstantCoefficient coeff(2.5);
    int q_order = 2;

    // Compute using IntegrateFunction.
    auto GetOrder = [q_order](const mfem::ElementTransformation &) { return q_order; };
    double result_new = fem::IntegrateFunction(*mesh, marker, true, coeff, GetOrder);

    // Compute using MFEM LinearForm.
    double result_mfem = IntegrateWithLinearForm(fespace, marker, coeff, q_order / 2);

    // For a unit square, boundary length = 4, so integral of 2.5 = 10.0.
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(10.0, 1e-10));
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(result_mfem, 1e-10));
  }

  SECTION("3D boundary integral with constant coefficient")
  {
    auto mesh = CreateTestMesh3D(comm);
    mfem::H1_FECollection fec(1, mesh->Dimension());
    mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);

    // Mark all boundary elements.
    int bdr_attr_max = mesh->bdr_attributes.Max();
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;

    ConstantCoefficient coeff(1.0);
    int q_order = 2;

    // Compute using IntegrateFunction.
    auto GetOrder = [q_order](const mfem::ElementTransformation &) { return q_order; };
    double result_new = fem::IntegrateFunction(*mesh, marker, true, coeff, GetOrder);

    // Compute using MFEM LinearForm.
    double result_mfem = IntegrateWithLinearForm(fespace, marker, coeff, q_order / 2);

    // For a unit cube, surface area = 6.
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(6.0, 1e-10));
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(result_mfem, 1e-10));
  }

  SECTION("3D boundary integral with linear coefficient")
  {
    auto mesh = CreateTestMesh3D(comm);
    mfem::H1_FECollection fec(2, mesh->Dimension());  // Higher order for accuracy
    mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);

    // Mark all boundary elements.
    int bdr_attr_max = mesh->bdr_attributes.Max();
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;

    LinearCoefficient coeff;
    int q_order = 4;  // Higher order for linear coefficient

    // Compute using IntegrateFunction.
    auto GetOrder = [q_order](const mfem::ElementTransformation &) { return q_order; };
    double result_new = fem::IntegrateFunction(*mesh, marker, true, coeff, GetOrder);

    // Compute using MFEM LinearForm.
    double result_mfem = IntegrateWithLinearForm(fespace, marker, coeff, q_order / 2);

    // Results should match between the two methods.
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(result_mfem, 1e-10));
  }

  SECTION("Partial boundary marking")
  {
    auto mesh = CreateTestMesh3D(comm);
    mfem::H1_FECollection fec(1, mesh->Dimension());
    mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);

    // Mark only some boundary attributes (e.g., first half).
    int bdr_attr_max = mesh->bdr_attributes.Max();
    mfem::Array<int> marker(bdr_attr_max);
    marker = 0;
    for (int i = 0; i < bdr_attr_max / 2; i++)
    {
      marker[i] = 1;
    }

    ConstantCoefficient coeff(1.0);
    int q_order = 2;

    // Compute using IntegrateFunction.
    auto GetOrder = [q_order](const mfem::ElementTransformation &) { return q_order; };
    double result_new = fem::IntegrateFunction(*mesh, marker, true, coeff, GetOrder);

    // Compute using MFEM LinearForm.
    double result_mfem = IntegrateWithLinearForm(fespace, marker, coeff, q_order / 2);

    // Results should match between the two methods.
    REQUIRE_THAT(result_new, Catch::Matchers::WithinRel(result_mfem, 1e-10));
  }
}

TEST_CASE("IntegrateFunction Benchmark", "[Benchmark][integrate][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  // Create a larger mesh for meaningful benchmarks.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(10, 10, 10, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  auto mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  serial_mesh.reset();

  mfem::H1_FECollection fec(2, mesh->Dimension());
  mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);

  int bdr_attr_max = mesh->bdr_attributes.Max();
  mfem::Array<int> marker(bdr_attr_max);
  marker = 1;

  LinearCoefficient coeff;
  int q_order = 4;

  BENCHMARK("IntegrateFunction")
  {
    auto GetOrder = [q_order](const mfem::ElementTransformation &) { return q_order; };
    return fem::IntegrateFunction(*mesh, marker, true, coeff, GetOrder);
  };

  BENCHMARK("MFEM LinearForm")
  {
    return IntegrateWithLinearForm(fespace, marker, coeff, q_order / 2);
  };
}

}  // namespace
}  // namespace palace
