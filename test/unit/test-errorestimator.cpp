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
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
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
        rt_fespaces, nd_fespace, tol, max_it, print, use_mg);
    Vector E(nd_fespace.GetTrueVSize()), D(rt_fespace.GetTrueVSize());
    FillRandom(E, 7);
    projector.Mult(E, D);
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
        nd_fespaces, rt_fespace, tol, max_it, print, use_mg);
    Vector B(rt_fespace.GetTrueVSize()), H(nd_fespace.GetTrueVSize());
    FillRandom(B, 11);
    projector.Mult(B, H);
    // The projector's coefficient is μ⁻¹, so its inverse weighting is μ.
    const auto [weighted, unweighted] =
        Residuals(nd_fespace, rt_fespace, B, H, inv_mu.pw, mu.pw);
    CHECK(weighted < 1.0e-10);
    CHECK(unweighted > 1.0e-3);
  }
}

}  // namespace palace
