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
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"

namespace palace
{
namespace
{

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

TEST_CASE("InnerProductCoefficient", "[coefficient][Serial]")
{
  SECTION("Computes dot product correctly")
  {
    // Create two constant vector coefficients: a = (1, 2, 3), b = (4, 5, 6).
    // Expected: a · b = 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32.
    mfem::Vector a_vec(3), b_vec(3);
    a_vec[0] = 1.0;
    a_vec[1] = 2.0;
    a_vec[2] = 3.0;
    b_vec[0] = 4.0;
    b_vec[1] = 5.0;
    b_vec[2] = 6.0;

    mfem::VectorConstantCoefficient a_coeff(a_vec), b_coeff(b_vec);
    InnerProductCoefficient ip_coeff(a_coeff, b_coeff);

    // Evaluate at a dummy transformation (the constant coefficients ignore it).
    mfem::IsoparametricTransformation T;
    mfem::IntegrationPoint ip;
    ip.Set3(0.5, 0.5, 0.5);

    double result = ip_coeff.Eval(T, ip);
    REQUIRE_THAT(result, Catch::Matchers::WithinRel(32.0, 1e-14));
  }
}

TEST_CASE("BdrInnerProductCoefficient", "[coefficient][Serial]")
{
  MPI_Comm comm = MPI_COMM_WORLD;

  SECTION("Computes boundary inner product correctly")
  {
    // Create a simple 3D mesh.
    auto serial_mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    auto mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
    serial_mesh.reset();

    // Create Nedelec space and grid function.
    mfem::ND_FECollection fec(1, mesh->Dimension());
    mfem::ParFiniteElementSpace fespace(mesh.get(), &fec);
    mfem::ParGridFunction gf(&fespace);

    // Set grid function to a constant vector field (1, 0, 0).
    mfem::Vector const_vec(3);
    const_vec[0] = 1.0;
    const_vec[1] = 0.0;
    const_vec[2] = 0.0;
    mfem::VectorConstantCoefficient const_coeff(const_vec);
    gf.ProjectCoefficient(const_coeff);

    // Create a coefficient that returns (2, 3, 4).
    mfem::Vector coeff_vec(3);
    coeff_vec[0] = 2.0;
    coeff_vec[1] = 3.0;
    coeff_vec[2] = 4.0;
    mfem::VectorConstantCoefficient vec_coeff(coeff_vec);

    // BdrInnerProductCoefficient should compute gf · vec_coeff = (1,0,0) · (2,3,4) = 2.
    BdrInnerProductCoefficient bdr_ip(gf, vec_coeff);

    // Mark all boundary elements.
    int bdr_attr_max = mesh->bdr_attributes.Max();
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;

    // Integrate over boundary - should be 2 * surface_area = 2 * 6 = 12.
    auto GetOrder = [](const mfem::ElementTransformation &) { return 2; };
    double result = fem::IntegrateFunction(*mesh, marker, true, bdr_ip, GetOrder);

    REQUIRE_THAT(result, Catch::Matchers::WithinRel(12.0, 1e-10));
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
