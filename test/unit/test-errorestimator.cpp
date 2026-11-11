// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <complex>
#include <memory>
#include <random>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/entityjacobi.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

auto LoadMesh(MPI_Comm comm, const std::string &input)
{
  mfem::Mesh smesh(input, 1, 1);
  smesh.EnsureNodes();
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    smesh.SetAttribute(i, 1);
  }
  smesh.SetAttributes();
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  return Mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
}

void FillRandom(Vector &v, unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto *h = v.HostWrite();
  for (int i = 0; i < v.Size(); i++)
  {
    h[i] = dist(gen);
  }
}

}  // namespace

// The error estimators integrate the flux error with one libCEED sub-operator per element
// geometry type, all reading the same field vectors. libCEED re-gathers a passive input
// only when its vector's state changes, and for complex fields the vectors are re-pointed
// from the real to the imaginary part between the two integration passes, so every
// sub-operator must see those updates: on a mixed mesh, a second estimate and the imaginary
// pass must not reuse the data of the first estimate for the geometry types other than the
// first.
TEST_CASE("Flux error estimators on mixed-geometry meshes",
          "[errorestimator][Serial][Parallel][GPU]")
{
  const auto comm = MPI_COMM_WORLD;
  auto mesh =
      LoadMesh(comm, std::string(PALACE_TEST_DATA_DIR "/mesh/fichera-mixed-p2.mesh"));
  const int dim = mesh.Dimension();
  REQUIRE(dim == 3);
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;

  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {2.0, 3.0, 4.0};
  material.mu_r.s = {1.0, 1.5, 2.0};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::EIGENMODE, mesh);

  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  auto &rt_fespace = rt_fespaces.GetFinestFESpace();

  // Tight flux projection so that the complex (paired) and real solves agree well beyond
  // the tolerance used below.
  constexpr double tol = 1.0e-14;
  constexpr int max_it = 1000, print = 0;
  constexpr bool use_mg = false;

  auto Compare = [&](const Vector &actual, const Vector &expected, const char *what)
  {
    REQUIRE(actual.Size() == mesh.GetNE());
    REQUIRE(expected.Size() == actual.Size());
    const auto *ha = actual.HostRead();
    const auto *he = expected.HostRead();
    int num_bad = 0;
    double max_rel = 0.0;
    for (int i = 0; i < actual.Size(); i++)
    {
      const double rel = std::abs(ha[i] - he[i]) / std::max(std::abs(he[i]), 1.0e-300);
      max_rel = std::max(max_rel, rel);
      num_bad += (rel > 1.0e-8);
    }
    INFO(what << ": max relative difference " << max_rel);
    CHECK(num_bad == 0);
  };

  // MakeEstimator() constructs a fresh estimator of the requested type; a fresh estimator's
  // first estimate is always correct and serves as the reference.
  auto Run = [&](auto MakeEstimator, const FiniteElementSpace &fespace)
  {
    ComplexVector X(fespace.GetTrueVSize());
    X.UseDevice(true);
    FillRandom(X.Real(), 7);
    FillRandom(X.Imag(), 11);

    ErrorIndicator ref_re, ref_im;
    MakeEstimator(Vector()).AddErrorIndicator(X.Real(), 1.0, ref_re);
    MakeEstimator(Vector()).AddErrorIndicator(X.Imag(), 1.0, ref_im);
    {
      // The two parts give different estimates, so that reading the wrong part is detected.
      const auto *hre = ref_re.Local().HostRead();
      const auto *him = ref_im.Local().HostRead();
      int distinct = 0;
      for (int i = 0; i < ref_re.Local().Size(); i++)
      {
        distinct += (std::abs(hre[i] - him[i]) > 1.0e-6 * std::max(hre[i], him[i]));
      }
      Mpi::GlobalSum(1, &distinct, comm);
      REQUIRE(distinct > 0);
    }

    SECTION("Repeated real estimates")
    {
      auto estimator = MakeEstimator(Vector());
      ErrorIndicator first, second;
      estimator.AddErrorIndicator(X.Real(), 1.0, first);
      estimator.AddErrorIndicator(X.Imag(), 1.0, second);
      Compare(first.Local(), ref_re.Local(), "first real estimate");
      Compare(second.Local(), ref_im.Local(), "second real estimate");
    }

    SECTION("Complex estimate")
    {
      auto estimator = MakeEstimator(ComplexVector());
      ErrorIndicator ind_c;
      estimator.AddErrorIndicator(X, 1.0, ind_c);
      Vector expected(ref_re.Local().Size());
      {
        const auto *hre = ref_re.Local().HostRead();
        const auto *him = ref_im.Local().HostRead();
        auto *he = expected.HostWrite();
        for (int i = 0; i < expected.Size(); i++)
        {
          he[i] = std::sqrt(hre[i] * hre[i] + him[i] * him[i]);
        }
      }
      Compare(ind_c.Local(), expected, "complex estimate");
    }
  };

  SECTION("Curl flux estimator")
  {
    Run(
        [&](auto vec_type)
        {
          return CurlFluxErrorEstimator<decltype(vec_type)>(mat_op, rt_fespace, nd_fespaces,
                                                            tol, max_it, print, use_mg);
        },
        rt_fespace);
  }

  SECTION("Gradient flux estimator")
  {
    Run(
        [&](auto vec_type)
        {
          return GradFluxErrorEstimator<decltype(vec_type)>(mat_op, nd_fespace, rt_fespaces,
                                                            tol, max_it, print, use_mg);
        },
        nd_fespace);
  }
}

// With the entity blocks, the flux projection stops on the residual norm of Jacobi
// preconditioning, so that on flat elements its error at a given tolerance is about that of
// Jacobi preconditioning, at most twice here (stopping on the residual norm of the blocks,
// it is two orders of magnitude larger).
TEST_CASE("Flux projection accuracy on flat elements", "[errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  smesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(2) *= 1.0e-3;
      });
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 3, max_it = 100000;
  constexpr double tol = 1.0e-6;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  MaterialPropertyCoefficient coeff(mesh.MaxCeedAttribute());
  coeff.AddMaterialProperty(mesh.GetCeedAttributes(1), 2.0);

  auto Check = [&](const mfem::FiniteElementCollection &smooth_fec,
                   const mfem::FiniteElementCollection &rhs_fec)
  {
    FiniteElementSpaceHierarchy smooth_fespaces(
        std::make_unique<FiniteElementSpace>(mesh, &smooth_fec));
    const auto &smooth_fespace = smooth_fespaces.GetFinestFESpace();
    FiniteElementSpace rhs_fespace(mesh, &rhs_fec);
    if (!EntityBlockJacobiSolver<Operator>::Supported(smooth_fespace))
    {
      SKIP("The flux projection does not use entity blocks with this backend");
    }
    BilinearForm m(smooth_fespace), flux(rhs_fespace, smooth_fespace);
    m.AddDomainIntegrator<VectorFEMassIntegrator>();
    flux.AddDomainIntegrator<VectorFEMassIntegrator>(coeff);
    ParOperator M(m.PartialAssemble(), smooth_fespace);
    ParOperator Flux(flux.PartialAssemble(), rhs_fespace, smooth_fespace, false);
    const int n = smooth_fespace.GetTrueVSize();
    Vector x(rhs_fespace.GetTrueVSize()), b(n), y(n), y_jacobi(n), y_ref(n), Mv(n);
    FillRandom(x, Mpi::Rank(comm) + 3);
    Flux.Mult(x, b);
    auto MNorm = [&](const Vector &v)
    {
      M.Mult(v, Mv);
      return std::sqrt(linalg::Dot(comm, Mv, v));
    };

    FluxProjector<Vector> proj(coeff, smooth_fespaces, rhs_fespace, tol, max_it, 0, false);
    const int its = proj.NumTotalMultIterations();
    proj.Mult(x, y);
    const int its_ebj = proj.NumTotalMultIterations() - its;
    FluxProjector<Vector> ref(coeff, smooth_fespaces, rhs_fespace, 1.0e-12, max_it, 0,
                              false);
    ref.Mult(x, y_ref);
    JacobiSmoother<Operator> jacobi(comm);
    jacobi.SetOperator(M);
    CgSolver<Operator> pcg(comm, 0);
    pcg.SetOperator(M);
    pcg.SetPreconditioner(jacobi);
    pcg.SetRelTol(tol);
    pcg.SetMaxIter(max_it);
    y_jacobi = 0.0;
    pcg.Mult(b, y_jacobi);
    REQUIRE(pcg.GetConverged());

    const double norm_ref = MNorm(y_ref);
    y -= y_ref;
    y_jacobi -= y_ref;
    const double err = MNorm(y) / norm_ref, err_jacobi = MNorm(y_jacobi) / norm_ref;
    INFO("relative M-norm errors " << err << " (Jacobi " << err_jacobi << "), iterations "
                                   << its_ebj << " (Jacobi " << pcg.GetNumIterations()
                                   << ")");
    CHECK(err <= 2.0 * err_jacobi);
    CHECK(2 * its_ebj <= pcg.GetNumIterations());
  };
  Check(rt_fec, nd_fec);
  Check(nd_fec, rt_fec);
}

// In 2D, B = ∇ × E is a scalar in an L2 space with INTEGRAL map type (B = 1 / det(J) ̂B).
// Its smooth H1 reconstruction is exact for a uniform B, so the estimate vanishes,
// including on a graded mesh where det(J) varies between elements.
TEST_CASE("2D curl flux error estimator", "[errorestimator][Serial][Parallel][GPU]")
{
  const auto comm = MPI_COMM_WORLD;
  const auto elem_type = GENERATE(mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL);
  auto smesh = mfem::Mesh::MakeCartesian2D(4, 4, elem_type);
  for (int i = 0; i < smesh.GetNV(); i++)
  {
    double *v = smesh.GetVertex(i);
    v[0] = v[0] * v[0];
    v[1] = 0.5 * v[1] * (1.0 + v[1]);
  }
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;

  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {2.0, 2.0, 2.0};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::EIGENMODE, mesh);

  mfem::ND_FECollection nd_fec(order, 2);
  mfem::H1_FECollection h1_fec(order, 2);
  mfem::L2_FECollection l2_fec(order - 1, 2, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), l2_fespace(mesh, &l2_fec);
  FiniteElementSpaceHierarchy h1_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &h1_fec));

  // B = ∇ × E = 1 for E = (-y, x) / 2.
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient E_coeff(2,
                                          [](const mfem::Vector &x, mfem::Vector &v)
                                          {
                                            v(0) = -0.5 * x(1);
                                            v(1) = 0.5 * x(0);
                                          });
  E.ProjectCoefficient(E_coeff);
  Vector e(nd_fespace.GetTrueVSize()), b(l2_fespace.GetTrueVSize());
  e.UseDevice(true);
  b.UseDevice(true);
  E.GetTrueDofs(e);
  l2_fespace.GetDiscreteInterpolator(nd_fespace).Mult(e, b);

  CurlFluxErrorEstimator<Vector> estimator(mat_op, l2_fespace, h1_fespaces, 1.0e-14, 1000,
                                           0, false);
  ErrorIndicator indicator;
  estimator.AddErrorIndicator(b, 1.0, indicator);
  CHECK(indicator.Max(comm) < 1.0e-10);

  // The error QFunction requires the INTEGRAL map type.
  mfem::L2_FECollection l2_value_fec(order - 1, 2);
  FiniteElementSpace l2_value_fespace(mesh, &l2_value_fec);
  CHECK_THROWS(CurlFluxErrorEstimator<Vector>(mat_op, l2_value_fespace, h1_fespaces,
                                              1.0e-14, 1000, 0, false));
}

// Consecutive flux projections into the same output vector start from the previous
// solution, scaled to the closest multiple of it: a flux proportional to the previous one
// (with a complex factor) needs no iterations, an unrelated one about as many as from zero
// (with Jacobi preconditioning, up to 1% more), and the result has the same accuracy.
TEST_CASE("Flux projection initial guess", "[errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  smesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(2) *= 0.1 + x(0);
      });
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  const auto &rt_fespace = rt_fespaces.GetFinestFESpace();
  MaterialPropertyCoefficient coeff(mesh.MaxCeedAttribute());
  coeff.AddMaterialProperty(mesh.GetCeedAttributes(1), 2.0);
  BilinearForm m(rt_fespace);
  m.AddDomainIntegrator<VectorFEMassIntegrator>();
  ParOperator M(m.PartialAssemble(), rt_fespace);
  constexpr int max_it = 10000, print = 0;

  auto Check = [&](auto vec_type, std::complex<double> s)
  {
    using VecType = decltype(vec_type);
    const int n = nd_fespace.GetTrueVSize(), nr = rt_fespace.GetTrueVSize();
    auto Random = [n](unsigned seed)
    {
      VecType x(n);
      if constexpr (std::is_same_v<VecType, ComplexVector>)
      {
        FillRandom(x.Real(), seed);
        FillRandom(x.Imag(), seed + 1);
      }
      else
      {
        FillRandom(x, seed);
      }
      return x;
    };
    auto MNorm = [&](const VecType &v)
    {
      Vector Mv(nr);
      if constexpr (std::is_same_v<VecType, ComplexVector>)
      {
        M.Mult(v.Real(), Mv);
        double vMv = linalg::Dot(comm, Mv, v.Real());
        M.Mult(v.Imag(), Mv);
        return std::sqrt(vMv + linalg::Dot(comm, Mv, v.Imag()));
      }
      else
      {
        M.Mult(v, Mv);
        return std::sqrt(linalg::Dot(comm, Mv, v));
      }
    };
    auto Solve = [](const FluxProjector<VecType> &proj, const VecType &x, VecType &y)
    {
      const int its = proj.NumTotalMultIterations();
      proj.Mult(x, y);
      return proj.NumTotalMultIterations() - its;
    };
    const VecType x1 = Random(Mpi::Rank(comm) + 1), x3 = Random(Mpi::Rank(comm) + 7);
    VecType x2(x1);
    if constexpr (std::is_same_v<VecType, ComplexVector>)
    {
      x2 *= s;
    }
    else
    {
      x2 *= s.real();
    }
    FluxProjector<VecType> proj(coeff, rt_fespaces, nd_fespace, 1.0e-10, max_it, print,
                                false);
    VecType y(nr), y_zero(nr), y_ref(nr);
    const int its1 = Solve(proj, x1, y);
    const int its2 = Solve(proj, x2, y);
    const int its3 = Solve(proj, x3, y);
    FluxProjector<VecType> zero(coeff, rt_fespaces, nd_fespace, 1.0e-10, max_it, print,
                                false);
    const int its3_zero = Solve(zero, x3, y_zero);
    FluxProjector<VecType> ref(coeff, rt_fespaces, nd_fespace, 1.0e-14, max_it, print,
                               false);
    Solve(ref, x3, y_ref);
    INFO("iterations " << its1 << ", proportional " << its2 << ", unrelated " << its3
                       << " (from zero " << its3_zero << ")");
    REQUIRE(its1 >= 20);
    CHECK(4 * its2 <= its1);
    CHECK(its3 <= its3_zero + its3_zero / 100 + 2);
    y -= y_ref;
    y_zero -= y_ref;
    const double norm_ref = MNorm(y_ref), err = MNorm(y) / norm_ref,
                 err_zero = MNorm(y_zero) / norm_ref;
    INFO("relative M-norm errors " << err << " (from zero " << err_zero << ")");
    CHECK(err <= 1.0e-8);
    CHECK(err_zero <= 1.0e-8);
  };
  Check(Vector(), 1.7);
  Check(ComplexVector(), {0.3, -1.2});
}

TEST_CASE("Flux projection initial guess from several solutions",
          "[errorestimator][Serial][Parallel]")
{
  // With the last num_prev solutions kept, a flux in the span of the fluxes of the last
  // num_prev projections is recovered from its initial guess in a few iterations, also
  // after the oldest solutions are dropped, and to the same accuracy.
  const auto comm = MPI_COMM_WORLD;
  auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  smesh.Transform(
      [](const mfem::Vector &x, mfem::Vector &y)
      {
        y = x;
        y(2) *= 0.1 + x(0);
      });
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  const auto &rt_fespace = rt_fespaces.GetFinestFESpace();
  MaterialPropertyCoefficient coeff(mesh.MaxCeedAttribute());
  coeff.AddMaterialProperty(mesh.GetCeedAttributes(1), 2.0);
  BilinearForm m(rt_fespace);
  m.AddDomainIntegrator<VectorFEMassIntegrator>();
  ParOperator M(m.PartialAssemble(), rt_fespace);
  constexpr int max_it = 10000, print = 0;
  const int n = nd_fespace.GetTrueVSize(), nr = rt_fespace.GetTrueVSize();
  auto Random = [n](unsigned seed)
  {
    ComplexVector x(n);
    FillRandom(x.Real(), seed);
    FillRandom(x.Imag(), seed + 1);
    return x;
  };
  auto MNorm = [&](const ComplexVector &v)
  {
    Vector Mv(nr);
    M.Mult(v.Real(), Mv);
    double vMv = linalg::Dot(comm, Mv, v.Real());
    M.Mult(v.Imag(), Mv);
    return std::sqrt(vMv + linalg::Dot(comm, Mv, v.Imag()));
  };
  auto Solve =
      [](const FluxProjector<ComplexVector> &proj, const ComplexVector &x, ComplexVector &y)
  {
    const int its = proj.NumTotalMultIterations();
    proj.Mult(x, y);
    return proj.NumTotalMultIterations() - its;
  };
  const int rank = Mpi::Rank(comm);
  const ComplexVector x1 = Random(rank + 1), x2 = Random(rank + 11), x3 = Random(rank + 21);
  ComplexVector x23(x2), x3s(x3);
  x23 *= std::complex<double>(0.4, 0.9);
  x3s *= std::complex<double>(-1.3, 0.2);
  x23 += x3s;  // In the span of the fluxes of the last two projections, but not of x3 alone
  FluxProjector<ComplexVector> proj(coeff, rt_fespaces, nd_fespace, 1.0e-10, max_it, print,
                                    false, 2);
  FluxProjector<ComplexVector> proj1(coeff, rt_fespaces, nd_fespace, 1.0e-10, max_it, print,
                                     false, 1);
  ComplexVector y(nr), y1(nr), y_ref(nr);
  const int its1 = Solve(proj, x1, y);
  Solve(proj, x2, y);
  Solve(proj, x3, y);  // x1 is dropped
  const int its23 = Solve(proj, x23, y);
  for (const auto *x : {&x1, &x2, &x3})
  {
    Solve(proj1, *x, y1);
  }
  const int its23_1 = Solve(proj1, x23, y1);
  FluxProjector<ComplexVector> ref(coeff, rt_fespaces, nd_fespace, 1.0e-14, max_it, print,
                                   false);
  Solve(ref, x23, y_ref);
  INFO("iterations " << its1 << ", in the span of the last two " << its23 << " (with one "
                     << its23_1 << ")");
  REQUIRE(its1 >= 20);
  CHECK(4 * its23 <= its1);
  CHECK(its23 < its23_1);
  y -= y_ref;
  y1 -= y_ref;
  const double norm_ref = MNorm(y_ref);
  INFO("relative M-norm errors " << MNorm(y) / norm_ref << " (with one "
                                 << MNorm(y1) / norm_ref << ")");
  CHECK(MNorm(y) <= 1.0e-8 * norm_ref);
  CHECK(MNorm(y1) <= 1.0e-8 * norm_ref);
}

TEST_CASE("ErrorIndicator weighted samples", "[errorestimator][Serial]")
{
  // An indicator of weight w counts as w samples of the running root mean square, and the
  // weight 1 gives the unweighted mean.
  Vector a(3), b(3);
  a(0) = 1.0;
  a(1) = 2.0;
  a(2) = 0.5;
  b(0) = 3.0;
  b(1) = 0.0;
  b(2) = 1.5;
  ErrorIndicator weighted, repeated, unweighted;
  weighted.AddIndicator(a);
  weighted.AddIndicator(b, 3.0);
  repeated.AddIndicator(a);
  for (int k = 0; k < 3; k++)
  {
    repeated.AddIndicator(b);
  }
  unweighted.AddIndicator(a);
  unweighted.AddIndicator(b);
  const auto *w = weighted.Local().HostRead(), *r = repeated.Local().HostRead(),
             *u = unweighted.Local().HostRead();
  for (int i = 0; i < 3; i++)
  {
    CHECK_THAT(w[i], WithinRel(r[i], 1.0e-14));
    CHECK(u[i] == std::sqrt(0.5 * (a(i) * a(i) + b(i) * b(i))));
  }
}

}  // namespace palace
