// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <random>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
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

// Diagonal tensor coefficient, a^p on attribute 1 and b^p on attribute 2.
struct PiecewiseDiagonal
{
  mfem::MatrixConstantCoefficient c1, c2;
  mfem::PWMatrixCoefficient pw;

  static mfem::DenseMatrix Diag(const double (&a)[3], double p)
  {
    mfem::DenseMatrix m(3);
    m = 0.0;
    for (int i = 0; i < 3; i++)
    {
      m(i, i) = std::pow(a[i], p);
    }
    return m;
  }

  PiecewiseDiagonal(const double (&a)[3], double b, double p)
    : c1(Diag(a, p)), c2(Diag({b, b, b}, p)),
      pw(3, mfem::Array<int>({1, 2}), mfem::Array<mfem::MatrixCoefficient *>({&c1, &c2}))
  {
  }
};

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
          return CurlFluxErrorEstimator<decltype(vec_type)>(
              mat_op, rt_fespace, nd_fespaces, tol, 0.0, max_it, print, use_mg);
        },
        rt_fespace);
  }

  SECTION("Gradient flux estimator")
  {
    Run(
        [&](auto vec_type)
        {
          return GradFluxErrorEstimator<decltype(vec_type)>(
              mat_op, nd_fespace, rt_fespaces, tol, 0.0, max_it, print, use_mg);
        },
        nd_fespace);
  }
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

  CurlFluxErrorEstimator<Vector> estimator(mat_op, l2_fespace, h1_fespaces, 1.0e-14, 0.0,
                                           1000, 0, false);
  ErrorIndicator indicator;
  estimator.AddErrorIndicator(b, 1.0, indicator);
  CHECK(indicator.Max(comm) < 1.0e-10);

  // The error QFunction requires the INTEGRAL map type.
  mfem::L2_FECollection l2_value_fec(order - 1, 2);
  FiniteElementSpace l2_value_fespace(mesh, &l2_value_fec);
  CHECK_THROWS(CurlFluxErrorEstimator<Vector>(mat_op, l2_value_fespace, h1_fespaces,
                                              1.0e-14, 0.0, 1000, 0, false));
}

// The flux error is measured in the energy norm, ‖ε^{-1/2} (D − ε E)‖ for the gradient flux
// and ‖μ^{1/2} (H − μ⁻¹ B)‖ for the curl flux, so the recovered flux must be the projection
// which minimizes it: (ψ, ε⁻¹ D) = (ψ, E) and (ψ, μ H) = (ψ, B) for every ψ in the smooth
// space. Across a material interface this differs from the unweighted projection
// (ψ, D) = (ψ, ε E), which is checked as well so that the test can tell the two apart.
TEST_CASE("Flux recovery minimizes the energy-norm error across material interfaces",
          "[errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  mfem::Mesh smesh =
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
  smesh.EnsureNodes();
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    mfem::Vector c;
    smesh.GetElementCenter(i, c);
    smesh.SetAttribute(i, (c(0) < 0.5) ? 1 : 2);
  }
  smesh.SetAttributes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;

  // Material 1 is anisotropic so that the inverse of a full tensor coefficient is
  // exercised.
  const double eps1[3] = {2.0, 3.0, 4.0}, mu1[3] = {1.0, 1.5, 2.0};
  constexpr double eps2 = 10.0, mu2 = 0.5;
  std::vector<config::MaterialData> materials(2);
  materials[0].attributes = {1};
  materials[0].epsilon_r.s = {eps1[0], eps1[1], eps1[2]};
  materials[0].mu_r.s = {mu1[0], mu1[1], mu1[2]};
  materials[1].attributes = {2};
  materials[1].epsilon_r.s.fill(eps2);
  materials[1].mu_r.s.fill(mu2);
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::EIGENMODE, mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  auto &rt_fespace = rt_fespaces.GetFinestFESpace();
  constexpr double tol = 1.0e-14;
  constexpr int max_it = 1000, print = 0;
  constexpr bool use_mg = false;

  // Relative residuals of the weighted and unweighted projection equations, assembled
  // directly with MFEM, for the flux G recovered in smooth_fespace from F in rhs_fespace.
  auto Residuals = [&](const FiniteElementSpace &smooth_fespace,
                       const FiniteElementSpace &rhs_fespace, const Vector &F,
                       const Vector &G, mfem::MatrixCoefficient &coeff,
                       mfem::MatrixCoefficient &inv_coeff)
  {
    auto Mass = [&](mfem::MatrixCoefficient *Q)
    {
      mfem::ParBilinearForm m(
          const_cast<mfem::ParFiniteElementSpace *>(&smooth_fespace.Get()));
      m.AddDomainIntegrator(Q ? new mfem::VectorFEMassIntegrator(*Q)
                              : new mfem::VectorFEMassIntegrator());
      m.Assemble();
      m.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(m.ParallelAssemble());
    };
    auto Mixed = [&](mfem::MatrixCoefficient *Q)
    {
      mfem::ParMixedBilinearForm b(
          const_cast<mfem::ParFiniteElementSpace *>(&rhs_fespace.Get()),
          const_cast<mfem::ParFiniteElementSpace *>(&smooth_fespace.Get()));
      b.AddDomainIntegrator(Q ? new mfem::VectorFEMassIntegrator(*Q)
                              : new mfem::VectorFEMassIntegrator());
      b.Assemble();
      b.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(b.ParallelAssemble());
    };
    auto Residual = [&](mfem::HypreParMatrix &M, mfem::HypreParMatrix &B)
    {
      mfem::Vector r(G.Size()), b(G.Size());
      B.Mult(F, b);
      M.Mult(G, r);
      r -= b;
      return std::sqrt(mfem::InnerProduct(comm, r, r) / mfem::InnerProduct(comm, b, b));
    };
    return std::pair{Residual(*Mass(&inv_coeff), *Mixed(nullptr)),
                     Residual(*Mass(nullptr), *Mixed(&coeff))};
  };

  SECTION("Gradient flux")
  {
    PiecewiseDiagonal eps(eps1, eps2, 1.0), inv_eps(eps1, eps2, -1.0);
    FluxProjector<Vector> projector(
        MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetPermittivityReal()),
        rt_fespaces, nd_fespace, tol, 0.0, max_it, print, use_mg);
    Vector E(nd_fespace.GetTrueVSize()), D(rt_fespace.GetTrueVSize());
    FillRandom(E, 7);
    projector.Mult(E, D, 0.0);  // No test relative to the estimate
    const auto [weighted, unweighted] =
        Residuals(rt_fespace, nd_fespace, E, D, eps.pw, inv_eps.pw);
    CHECK(weighted < 1.0e-10);
    CHECK(unweighted > 1.0e-3);
  }

  SECTION("Curl flux")
  {
    PiecewiseDiagonal mu(mu1, mu2, 1.0), inv_mu(mu1, mu2, -1.0);
    FluxProjector<Vector> projector(
        MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetCurlCurlInvPermeability()),
        nd_fespaces, rt_fespace, tol, 0.0, max_it, print, use_mg);
    Vector B(rt_fespace.GetTrueVSize()), H(nd_fespace.GetTrueVSize());
    FillRandom(B, 11);
    projector.Mult(B, H, 0.0);
    // The projector's coefficient is μ⁻¹, so its inverse weighting is μ.
    const auto [weighted, unweighted] =
        Residuals(nd_fespace, rt_fespace, B, H, inv_mu.pw, mu.pw);
    CHECK(weighted < 1.0e-10);
    CHECK(unweighted > 1.0e-3);
  }
}

// The estimate of any G in the recovery space is ‖G - D‖² = ‖G* - D‖² + ‖G - G*‖²_M for the
// projection G* (in the norm of the estimate, with the projector's mass matrix M), so the
// estimate of an iterate of the projection exceeds that of the projection by its squared
// error. The projection stopped at a tolerance tol has an error of at most tol times the
// square root of its estimate. The estimates are assembled directly with MFEM.
TEST_CASE("Flux projection error relative to the estimate",
          "[errorestimator][Serial][Parallel]")
{
  const auto comm = MPI_COMM_WORLD;
  mfem::Mesh smesh =
      mfem::Mesh::MakeCartesian3D(4, 4, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 0.01);
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    mfem::Vector c;
    smesh.GetElementCenter(i, c);
    smesh.SetAttribute(i, (c(0) < 0.5) ? 1 : 2);
  }
  smesh.SetAttributes();
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  const double eps1[3] = {2.0, 3.0, 4.0}, mu1[3] = {1.0, 1.5, 2.0};
  constexpr double eps2 = 10.0, mu2 = 0.5;
  std::vector<config::MaterialData> materials(2);
  materials[0].attributes = {1};
  materials[0].epsilon_r.s = {eps1[0], eps1[1], eps1[2]};
  materials[0].mu_r.s = {mu1[0], mu1[1], mu1[2]};
  materials[1].attributes = {2};
  materials[1].epsilon_r.s.fill(eps2);
  materials[1].mu_r.s.fill(mu2);
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::EIGENMODE, mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  auto &rt_fespace = rt_fespaces.GetFinestFESpace();

  // Projection of the flux of F in rhs_fespace onto the finest space of smooth_fespaces,
  // with the projector's coefficient coeff, c as an MFEM coefficient and inv_c its inverse.
  auto Check = [&](const MaterialPropertyCoefficient &coeff,
                   const FiniteElementSpaceHierarchy &smooth_fespaces,
                   const FiniteElementSpace &rhs_fespace, mfem::MatrixCoefficient &c,
                   mfem::MatrixCoefficient &inv_c)
  {
    const auto &smooth_fespace = smooth_fespaces.GetFinestFESpace();
    auto *smooth = const_cast<mfem::ParFiniteElementSpace *>(&smooth_fespace.Get());
    auto *rhs = const_cast<mfem::ParFiniteElementSpace *>(&rhs_fespace.Get());
    auto Assemble = [](auto &form)
    {
      form.Assemble();
      form.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(form.ParallelAssemble());
    };
    mfem::ParBilinearForm m(smooth), k(rhs);
    m.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(inv_c));
    k.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(c));
    mfem::ParMixedBilinearForm mixed(rhs, smooth);
    mixed.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
    auto M = Assemble(m), K = Assemble(k), B = Assemble(mixed);

    Vector F(rhs_fespace.GetTrueVSize()), b(smooth_fespace.GetTrueVSize()),
        t(rhs_fespace.GetTrueVSize());
    FillRandom(F, 17);
    K->Mult(F, t);
    B->Mult(F, b);
    const double flux_norm2 = mfem::InnerProduct(comm, F, t);
    // The estimate ‖G - D‖² = ‖D‖² - 2 (F, G) + (coeff⁻¹ G, G), and the squared error.
    auto Estimate = [&](const Vector &G)
    {
      Vector MG(G.Size());
      M->Mult(G, MG);
      return flux_norm2 - 2.0 * mfem::InnerProduct(comm, b, G) +
             mfem::InnerProduct(comm, G, MG);
    };
    auto Error2 = [&](const Vector &G, const Vector &G0)
    {
      Vector e(G), Me(G.Size());
      e -= G0;
      M->Mult(e, Me);
      return mfem::InnerProduct(comm, e, Me);
    };
    // The projection with the residual tolerance res_tol and the tolerance error_tol
    // relative to the estimate.
    auto Project = [&](double res_tol, double error_tol)
    {
      Vector G(smooth_fespace.GetTrueVSize());
      FluxProjector<Vector>(coeff, smooth_fespaces, rhs_fespace, res_tol, error_tol, 10000,
                            0, false)
          .Mult(F, G, flux_norm2);
      return G;
    };

    // The projection to a relative residual of 10⁻¹³, and the Pythagoras identity for an
    // iterate and for an arbitrary G.
    const Vector G0 = Project(1.0e-13, 0.0);
    const double eta2 = Estimate(G0);
    REQUIRE(eta2 > 1.0e-3 * flux_norm2);
    Vector Gr(G0.Size());
    FillRandom(Gr, 19);
    for (const auto &G : {Project(0.0, 1.0e-1), Gr})
    {
      const double excess = Estimate(G) - eta2, err2 = Error2(G, G0);
      CHECK_THAT(excess, WithinRel(err2, 1.0e-6) || WithinAbs(err2, 1.0e-12 * flux_norm2));
    }

    // The error of the projection stopped at the tolerance relative to the estimate
    // (without the residual test), relative to its estimate, and the tolerance stops the
    // iterations before the rounding error.
    for (const double tol : {1.0e-1, 3.0e-2, 1.0e-2, 1.0e-3})
    {
      const Vector G = Project(0.0, tol);
      const double eta2_tol = Estimate(G), err2 = Error2(G, G0);
      INFO("tol = " << tol
                    << ": error² / (tol² estimate) = " << err2 / (tol * tol * eta2_tol));
      CHECK(err2 <= tol * tol * eta2_tol);
      CHECK(eta2_tol >= eta2 * (1.0 - 1.0e-12));
      if (tol == 1.0e-1)
      {
        CHECK(err2 >= 1.0e-8 * eta2_tol);
      }
    }
  };

  SECTION("Gradient flux")
  {
    PiecewiseDiagonal eps(eps1, eps2, 1.0), inv_eps(eps1, eps2, -1.0);
    Check(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetPermittivityReal()),
          rt_fespaces, nd_fespace, eps.pw, inv_eps.pw);
  }

  SECTION("Curl flux")
  {
    // The projector's coefficient is μ⁻¹, so its inverse weighting is μ.
    PiecewiseDiagonal mu(mu1, mu2, 1.0), inv_mu(mu1, mu2, -1.0);
    Check(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetCurlCurlInvPermeability()),
          nd_fespaces, rt_fespace, inv_mu.pw, mu.pw);
  }
}

// The estimator stops the projection relative to the estimate (here without the residual
// test): its estimate is never below that of the converged projection and exceeds it by at
// most tol² (relative), and the element indicators differ by at most tol times the estimate
// in the ℓ² norm.
TEST_CASE("Flux error estimator with the projection tolerance relative to the estimate",
          "[errorestimator][Serial][Parallel][GPU]")
{
  const auto comm = MPI_COMM_WORLD;
  auto smesh =
      mfem::Mesh::MakeCartesian3D(4, 4, 2, mfem::Element::TETRAHEDRON, 1.0, 1.0, 0.01);
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    mfem::Vector c;
    smesh.GetElementCenter(i, c);
    smesh.SetAttribute(i, (c(0) < 0.5) ? 1 : 2);
  }
  smesh.SetAttributes();
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  std::vector<config::MaterialData> materials(2);
  materials[0].attributes = {1};
  materials[0].epsilon_r.s.fill(1.0);
  materials[1].attributes = {2};
  materials[1].epsilon_r.s.fill(10.0);
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::ELECTROSTATIC, mesh);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(mesh, &rt_fec));
  Vector E(nd_fespace.GetTrueVSize());
  E.UseDevice(true);
  FillRandom(E, 13);

  ErrorIndicator ind_ref;
  GradFluxErrorEstimator<Vector>(mat_op, nd_fespace, rt_fespaces, 1.0e-13, 0.0, 10000, 0,
                                 false)
      .AddErrorIndicator(E, 1.0, ind_ref);
  const double eta_ref = ind_ref.Norml2(comm);
  for (const double tol : {1.0e-1, 1.0e-2})
  {
    ErrorIndicator ind;
    GradFluxErrorEstimator<Vector>(mat_op, nd_fespace, rt_fespaces, 0.0, tol, 10000, 0,
                                   false)
        .AddErrorIndicator(E, 1.0, ind);
    const double eta = ind.Norml2(comm);
    Vector diff(ind.Local());
    diff -= ind_ref.Local();
    INFO("tol = " << tol << ": (estimate² - reference²) / (tol² estimate²) = "
                  << (eta * eta - eta_ref * eta_ref) / (tol * tol * eta * eta)
                  << ", ‖indicator - reference‖ / (tol estimate) = "
                  << linalg::Norml2(comm, diff) / (tol * eta));
    CHECK(eta >= eta_ref * (1.0 - 1.0e-12));
    CHECK(eta * eta - eta_ref * eta_ref <= tol * tol * eta * eta);
    CHECK(linalg::Norml2(comm, diff) <= tol * eta);
    if (tol == 1.0e-1)
    {
      // The tolerance stops the iterations before the rounding error.
      CHECK(eta * eta - eta_ref * eta_ref >= 1.0e-8 * eta * eta);
    }
  }
}

}  // namespace palace
