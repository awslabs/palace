// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "test-helpers.hpp"

#include "fem/bilinearform.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"

namespace palace
{
using json = nlohmann::json;
using namespace std::complex_literals;
using namespace Catch;
namespace fs = std::filesystem;

TEST_CASE("BuildParSumOperator", "[rap][Serial][Parallel]")
{
  // Verify that BuildParSumOperator can assemble collections of ParOperators.
  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;

  auto comm = Mpi::World();

  mfem::Mesh mfem_mesh(SingleTetMesh());
  while (mfem_mesh.GetNE() < Mpi::Size(comm))
  {
    mfem_mesh.UniformRefinement();  // avoid empty partition
  }
  Mesh mesh(comm, mfem_mesh);

  constexpr int order = 2, dim = 3;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, dim, mfem::BasisType::GaussLobatto,
                               mfem::BasisType::GaussLegendre);
  FiniteElementSpace nd_fes(mesh, &nd_fec);
  MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, mesh);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute());

  df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvPermeability(), 1.0);
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), 1.0);
  SECTION("ParOperator")
  {
    int empty = {df.empty() && f.empty()};
    Mpi::GlobalMin(1, &empty, comm);
    REQUIRE(empty == 0);  // There must be a non-empty.

    BilinearForm da(nd_fes), a(nd_fes);
    da.AddDomainIntegrator<CurlCurlIntegrator>(df);
    a.AddDomainIntegrator<VectorFEMassIntegrator>(df);

    constexpr bool skip_zeros = false;
    std::unique_ptr<Operator> DA =
        std::make_unique<ParOperator>(da.Assemble(skip_zeros), nd_fes);
    std::unique_ptr<Operator> A =
        std::make_unique<ParOperator>(a.Assemble(skip_zeros), nd_fes);

    REQUIRE(DA);
    REQUIRE(A);

    constexpr double c1 = 1.1, c2 = 2.3;
    auto sum = BuildParSumOperator({c1, 2.3, c2},
                                   {DA.get(), std::unique_ptr<Operator>().get(), A.get()});

    // Check that the action of the sum operator is equivalent to the weighted sum of the
    // actions of the components.
    Vector v0(DA->Height()), x1(DA->Height()), x2(DA->Height());
    v0 = 1.5;
    x1 = 0.0;
    x2 = 0.0;

    sum->Mult(v0, x1);
    DA->AddMult(v0, x2, c1);
    A->AddMult(v0, x2, c2);

    x1 -= x2;
    x1.Abs();
    double constexpr tol = 1e-12;
    CHECK(x1.Max() < tol);
    CHECK(x1.Max() < tol);

    // Exercise BaseSumOperator's unit-weight first-term initialization path directly.
    SumOperator direct_sum(DA->Height());
    direct_sum.AddOperator(*DA);
    direct_sum.AddOperator(*A, c2);
    x1 = 100.0;  // The first term must fully overwrite the output.
    direct_sum.Mult(v0, x1);
    DA->Mult(v0, x2);
    A->AddMult(v0, x2, c2);
    x1 -= x2;
    CHECK(x1.Normlinf() < tol);
    x1 = 100.0;
    direct_sum.MultTranspose(v0, x1);
    DA->MultTranspose(v0, x2);
    A->AddMultTranspose(v0, x2, c2);
    x1 -= x2;
    CHECK(x1.Normlinf() < tol);

    // An empty sum zeroes its output.
    SumOperator empty_sum(DA->Height());
    empty_sum.Mult(v0, x1);
    CHECK(x1.Normlinf() == 0.0);
  }

  SECTION("ComplexParOperator")
  {
    int empty = {df.empty() && f.empty()};
    Mpi::GlobalMin(1, &empty, comm);
    REQUIRE(empty == 0);  // There must be a non-empty.

    BilinearForm dar(nd_fes), dai(nd_fes), ar(nd_fes), ai(nd_fes);
    dar.AddDomainIntegrator<CurlCurlIntegrator>(df);
    dai.AddDomainIntegrator<VectorFEMassIntegrator>(f);
    ar.AddDomainIntegrator<VectorFEMassIntegrator>(df);
    ai.AddDomainIntegrator<CurlCurlIntegrator>(f);

    constexpr bool skip_zeros = false;
    std::unique_ptr<ComplexOperator> DA = std::make_unique<ComplexParOperator>(
        dar.Assemble(skip_zeros), dai.Assemble(skip_zeros), nd_fes);
    std::unique_ptr<ComplexOperator> A = std::make_unique<ComplexParOperator>(
        ar.Assemble(skip_zeros), ai.Assemble(skip_zeros), nd_fes);

    REQUIRE(DA);
    REQUIRE(A);

    const std::complex<double> c1 = 1.1 + 3.4i, c2 = 2.3 + 0.3i;
    auto sum = BuildParSumOperator(
        {c1, 3.1 + 0i, c2}, {DA.get(), std::unique_ptr<ComplexOperator>().get(), A.get()});

    // Check that the action of the sum operator is equivalent to the weighted sum of the
    // actions of the components.
    ComplexVector v0(DA->Height()), x1(DA->Height()), x2(DA->Height());
    v0.Real() = 1.5;
    v0.Imag() = 0.6;
    x1 = 0.0;
    x2 = 0.0;

    sum->Mult(v0, x1);
    DA->AddMult(v0, x2, c1);
    A->AddMult(v0, x2, c2);

    x1 -= x2;
    x1.Abs();
    double constexpr tol = 1e-12;
    CHECK(x1.Real().Max() < tol);
    CHECK(x1.Imag().Max() < tol);

    SumComplexOperator direct_sum(DA->Height());
    direct_sum.AddOperator(*DA);
    direct_sum.AddOperator(*A, c2);
    x1.Real() = 100.0;  // The first term must fully overwrite the output.
    x1.Imag() = 100.0;
    direct_sum.Mult(v0, x1);
    DA->Mult(v0, x2);
    A->AddMult(v0, x2, c2);
    x1 -= x2;
    CHECK(std::max(x1.Real().Normlinf(), x1.Imag().Normlinf()) < tol);
    x1.Real() = 100.0;
    x1.Imag() = 100.0;
    direct_sum.MultTranspose(v0, x1);
    DA->MultTranspose(v0, x2);
    A->AddMultTranspose(v0, x2, c2);
    x1 -= x2;
    CHECK(std::max(x1.Real().Normlinf(), x1.Imag().Normlinf()) < tol);

    // An empty sum zeroes its output.
    SumComplexOperator empty_sum(DA->Height());
    empty_sum.Mult(v0, x1);
    CHECK(std::max(x1.Real().Normlinf(), x1.Imag().Normlinf()) == 0.0);
  }
}

}  // namespace palace
