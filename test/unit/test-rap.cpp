// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "test-helpers.hpp"

#include "drivers/eigensolver.hpp"
#include "fem/bilinearform.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
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
  }

  SECTION("Essential diagonal wrapper")
  {
    BilinearForm ar(nd_fes), ai(nd_fes);
    ar.AddDomainIntegrator<VectorFEMassIntegrator>(df);
    ai.AddDomainIntegrator<CurlCurlIntegrator>(f);

    constexpr bool skip_zeros = false;
    ComplexParOperator A(ar.Assemble(skip_zeros), ai.Assemble(skip_zeros), nd_fes);
    mfem::Array<int> essential_tdofs(1);
    essential_tdofs[0] = 0;
    A.SetEssentialTrueDofs(essential_tdofs, Operator::DiagonalPolicy::DIAG_ONE);

    auto real_sum = std::make_unique<SumOperator>(*A.Real(), *A.Real());
    auto corrected_real = ApplyEssentialDiagonal(std::move(real_sum), *A.Real());
    ComplexWrapperOperator corrected(corrected_real.get(), A.Imag());

    ComplexVector x(A.Width()), y(A.Height());
    x.Real() = 1.5;
    x.Imag() = 0.6;
    corrected.Mult(x, y);
    CHECK(y.Real()[0] == Catch::Approx(x.Real()[0]));
    CHECK(y.Imag()[0] == Catch::Approx(x.Imag()[0]));

    Vector yr(A.Height()), yi(A.Height());
    corrected.Real()->Mult(x.Real(), yr);
    corrected.Imag()->Mult(x.Real(), yi);
    CHECK(yr[0] == Catch::Approx(x.Real()[0]));
    CHECK(yi[0] == Catch::Approx(0.0));
  }
}

namespace
{

// A rank-one complex operator y = u (v . x): a stand-in for the matrix-free wave-port modal
// correction W in the extra system operator A2. u and v vanish on the essential true dofs.
class RankOneComplexOperator : public ComplexOperator
{
private:
  ComplexVector u, v;
  MPI_Comm comm;

public:
  RankOneComplexOperator(const ComplexVector &u, const ComplexVector &v, MPI_Comm comm)
    : ComplexOperator(u.Size()), u(u), v(v), comm(comm)
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y = u;
    y *= linalg::Dot(comm, v, x);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    y.Add(a * linalg::Dot(comm, v, x), u);
  }
};

}  // namespace

// SD-402 / SD-408 (P3-08): BuildOperatorWithA2 keeps every term when a stiffness, damping
// or mass operator is not a ComplexParOperator (the response-corrected wrappers of the
// corrected eigenproblem) or A2 is matrix-free, and the sum keeps the essential-row policy
// of the stiffness: the identity on the essential rows (K DIAG_ONE; C, M, A2 zero there).
TEST_CASE("BuildOperatorWithA2 abstract terms", "[rap][Serial][Parallel]")
{
  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  auto comm = Mpi::World();
  mfem::Mesh mfem_mesh(SingleTetMesh());
  while (mfem_mesh.GetNE() < 2 * Mpi::Size(comm))
  {
    mfem_mesh.UniformRefinement();  // avoid empty partitions
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

  // Every exterior face is essential.
  mfem::Array<int> ess_bdr(mesh.Get().bdr_attributes.Max()), ess_tdofs;
  ess_bdr = 1;
  nd_fes.Get().GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  int num_ess = ess_tdofs.Size(), num_tdofs = nd_fes.GetTrueVSize();
  Mpi::GlobalSum(1, &num_ess, comm);
  Mpi::GlobalSum(1, &num_tdofs, comm);
  REQUIRE(num_ess > 0);
  REQUIRE(num_ess < num_tdofs);

  constexpr bool skip_zeros = false;
  BilinearForm kf(nd_fes), cf(nd_fes), mf(nd_fes), a2rf(nd_fes), a2if(nd_fes);
  kf.AddDomainIntegrator<CurlCurlIntegrator>(df);
  cf.AddDomainIntegrator<VectorFEMassIntegrator>(f);
  mf.AddDomainIntegrator<VectorFEMassIntegrator>(df);
  a2rf.AddDomainIntegrator<VectorFEMassIntegrator>(f);
  a2if.AddDomainIntegrator<CurlCurlIntegrator>(df);
  ComplexParOperator K(kf.Assemble(skip_zeros), nullptr, nd_fes);
  ComplexParOperator C(cf.Assemble(skip_zeros), nullptr, nd_fes);
  ComplexParOperator M(mf.Assemble(skip_zeros), nullptr, nd_fes);
  ComplexParOperator A2_sparse(a2rf.Assemble(skip_zeros), a2if.Assemble(skip_zeros),
                               nd_fes);
  K.SetEssentialTrueDofs(ess_tdofs, Operator::DiagonalPolicy::DIAG_ONE);
  C.SetEssentialTrueDofs(ess_tdofs, Operator::DiagonalPolicy::DIAG_ZERO);
  M.SetEssentialTrueDofs(ess_tdofs, Operator::DiagonalPolicy::DIAG_ZERO);
  A2_sparse.SetEssentialTrueDofs(ess_tdofs, Operator::DiagonalPolicy::DIAG_ZERO);

  const int n = K.Height();
  ComplexVector u(n), v(n), x(n), y(n), expected(n);
  for (auto *w : {&u, &v, &x, &y, &expected})
  {
    w->UseDevice(true);
  }
  linalg::SetRandom(comm, u, 11);
  linalg::SetRandom(comm, v, 12);
  linalg::SetRandom(comm, x, 13);
  linalg::SetSubVector(u, ess_tdofs, 0.0);
  linalg::SetSubVector(v, ess_tdofs, 0.0);
  RankOneComplexOperator W(u, v, comm);
  SumComplexOperator A2(A2_sparse);  // An abstract A2: the sparse part + W.
  A2.AddOperator(W);

  const std::complex<double> l = 0.3 + 1.7i;
  const std::vector<std::complex<double>> coeff = {1.0, l, l * l};

  // The action of op against the term-by-term sum of the terms, every row, and the
  // identity on the essential rows.
  using Terms = std::vector<std::pair<std::complex<double>, const ComplexOperator *>>;
  auto CheckSum = [&](const ComplexOperator &op, const Terms &terms)
  {
    ComplexVector term_x(n);
    term_x.UseDevice(true);
    expected = 0.0;
    for (const auto &[c, term] : terms)
    {
      term->Mult(x, term_x);
      expected.Add(c, term_x);
    }
    op.Mult(x, y);
    ComplexVector difference(y);
    linalg::AXPY(-1.0, expected, difference);
    CHECK(linalg::Norml2(comm, difference) / linalg::Norml2(comm, expected) < 1.0e-12);
    difference = y;
    linalg::AXPY(-1.0, x, difference);
    const auto *dr = difference.Real().HostRead();
    const auto *di = difference.Imag().HostRead();
    double essential_error = 0.0;
    for (int i = 0; i < ess_tdofs.Size(); i++)
    {
      essential_error +=
          dr[ess_tdofs[i]] * dr[ess_tdofs[i]] + di[ess_tdofs[i]] * di[ess_tdofs[i]];
    }
    Mpi::GlobalSum(1, &essential_error, comm);
    CHECK(std::sqrt(essential_error) <= 1.0e-12 * linalg::Norml2(comm, x));
  };

  SECTION("(i) parallel K, C, M and an abstract A2")
  {
    auto op = BuildOperatorWithA2(coeff, {&K, &C, &M}, &A2);
    REQUIRE(op);
    CheckSum(*op, {{1.0, &K}, {l, &C}, {l * l, &M}, {1.0, &A2}});
    auto reference =
        BuildComplexSumOperator({{1.0, &K}, {l, &C}, {l * l, &M}, {1.0, &A2}}, K);
    CheckSum(*reference, {{1.0, &K}, {l, &C}, {l * l, &M}, {1.0, &A2}});
  }

  SECTION("(ii) an abstract M is included")
  {
    auto M_real =
        ApplyEssentialDiagonal(std::make_unique<SumOperator>(*M.Real()), *M.Real());
    ComplexWrapperOperator M_abstract(M_real.get(), nullptr);
    auto op = BuildOperatorWithA2(coeff, {&K, &C, &M_abstract}, &A2);
    REQUIRE(op);
    CheckSum(*op, {{1.0, &K}, {l, &C}, {l * l, &M}, {1.0, &A2}});
    // The M term is not negligible here: dropping it (the unextended dispatcher's
    // dynamic_cast) would fail the check above.
    M.Mult(x, expected);
    CHECK(std::abs(l * l) * linalg::Norml2(comm, expected) >
          1.0e-3 * linalg::Norml2(comm, y));
  }

  SECTION("(iii) response-corrected K and M and an abstract A2")
  {
    auto K_corrected = internal::BuildResponseCorrectedStiffness(K, *M.Real());
    auto M_corrected = internal::BuildResponseCorrectedMass(M, *M.Real());
    REQUIRE(K_corrected.op);
    REQUIRE(M_corrected.op);
    auto op =
        BuildOperatorWithA2(coeff, {K_corrected.op.get(), &C, M_corrected.op.get()}, &A2);
    REQUIRE(op);
    const Terms terms = {
        {1.0, K_corrected.op.get()}, {l, &C}, {l * l, M_corrected.op.get()}, {1.0, &A2}};
    CheckSum(*op, terms);
    auto reference = BuildComplexSumOperator(
        {{1.0, K_corrected.op.get()}, {l, &C}, {l * l, M_corrected.op.get()}, {1.0, &A2}},
        *K_corrected.op);
    CheckSum(*reference, terms);
  }
}

}  // namespace palace
