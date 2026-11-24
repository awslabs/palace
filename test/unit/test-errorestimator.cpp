// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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
      // Complex electric fields require their frequency; this constant material does not
      // depend on it, so any fixture frequency gives the same estimate.
      constexpr double omega = 1.0;
      if constexpr (requires { estimator.AddErrorIndicator(X, 1.0, omega, ind_c); })
      {
        estimator.AddErrorIndicator(X, 1.0, omega, ind_c);
      }
      else
      {
        estimator.AddErrorIndicator(X, 1.0, ind_c);  // Magnetic flux: no frequency
      }
      Vector expected(ref_re.Local().Size());
      {
        const auto *hre = ref_re.Local().HostRead();
        const auto *him = ref_im.Local().HostRead();
        auto *he = expected.HostWrite();
        for (int i = 0; i < expected.Size(); i++)
        {
          // At the same physical energy denominator, a peak phasor's squared
          // error has half the weight of each instantaneous real-field result.
          he[i] = std::sqrt(0.5 * (hre[i] * hre[i] + him[i] * him[i]));
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

namespace
{

// The unit square or cube (4 simplices per direction) with an interior boundary sheet of
// attribute sheet_attr on the plane x = 0.5 which is NOT cracked (its faces are shared by
// the elements on both sides, as for the new default crack policy), the exterior boundary
// attributes of MakeCartesian2D / MakeCartesian3D (1-4 / 1-6), the element attribute 1 for
// y < 0.5 and 2 above, and a partition with process boundaries on and across the sheet.
constexpr int sheet_attr = 7;

Mesh MakeSheetMesh(MPI_Comm comm, int dim)
{
  auto base =
      (dim == 2)
          ? mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0)
          : mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
  std::vector<int> sheet_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    base.GetFaceVertices(f, fv);
    if (e2 >= 0 && std::all_of(fv.begin(), fv.end(), [&](int v)
                               { return std::abs(base.GetVertex(v)[0] - 0.5) < 1.0e-12; }))
    {
      sheet_faces.push_back(f);
    }
  }
  mfem::Mesh smesh(dim, base.GetNV(), base.GetNE(),
                   base.GetNBE() + static_cast<int>(sheet_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    smesh.AddVertex(base.GetVertex(v));
  }
  mfem::Vector c(dim);
  for (int e = 0; e < base.GetNE(); e++)
  {
    base.GetElementCenter(e, c);
    auto *el = base.GetElement(e)->Duplicate(&smesh);
    el->SetAttribute((c(1) < 0.5) ? 1 : 2);
    smesh.AddElement(el);
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    smesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&smesh));
  }
  for (const int f : sheet_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&smesh);
    el->SetAttribute(sheet_attr);
    smesh.AddBdrElement(el);
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  const int num_procs = Mpi::Size(comm);
  if (num_procs > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  std::vector<int> part(smesh.GetNE());
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    smesh.GetElementCenter(e, c);
    part[e] = ((c(0) > 0.5) + 2 * (c(1) > 0.5) + (dim > 2 && c(2) > 0.5)) % num_procs;
  }
  smesh.EnsureNodes();
  return Mesh(std::make_unique<mfem::ParMesh>(comm, smesh, part.data()));
}

std::vector<config::MaterialData> SheetMaterials()
{
  config::MaterialData substrate, vacuum;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  vacuum.attributes = {2};
  return {substrate, vacuum};
}

// The estimators' spaces as the drivers build them (ND of order p, RT of order p - 1; in 2D
// also H1 and the INTEGRAL L2 space of B = curl E, as the boundary mode driver).
struct SheetEstimatorSetup
{
  Mesh mesh;
  MaterialOperator mat_op;
  mfem::ND_FECollection nd_fec;
  mfem::RT_FECollection rt_fec;
  mfem::H1_FECollection h1_fec;
  mfem::L2_FECollection l2_fec;
  FiniteElementSpaceHierarchy nd_fespaces, rt_fespaces, h1_fespaces;
  FiniteElementSpace l2_fespace;

  SheetEstimatorSetup(MPI_Comm comm, int dim, int order)
    : mesh(MakeSheetMesh(comm, dim)),
      mat_op(SheetMaterials(), config::PeriodicBoundaryData(), ProblemType::EIGENMODE,
             mesh),
      nd_fec(order, dim), rt_fec(order - 1, dim), h1_fec(order, dim),
      l2_fec(order - 1, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL),
      nd_fespaces(std::make_unique<FiniteElementSpace>(mesh, &nd_fec)),
      rt_fespaces(std::make_unique<FiniteElementSpace>(mesh, &rt_fec)),
      h1_fespaces(std::make_unique<FiniteElementSpace>(mesh, &h1_fec)),
      l2_fespace(mesh, &l2_fec)
  {
  }
};

// The element-wise indicators of two forms agree to round-off (relative to the largest).
void CheckSameIndicator(MPI_Comm comm, const Vector &actual, const Vector &expected,
                        const char *what)
{
  REQUIRE(actual.Size() == expected.Size());
  double max_expected = expected.Size() ? expected.Max() : 0.0, max_diff = 0.0;
  const auto *ha = actual.HostRead();
  const auto *he = expected.HostRead();
  for (int i = 0; i < actual.Size(); i++)
  {
    max_diff = std::max(max_diff, std::abs(ha[i] - he[i]));
  }
  Mpi::GlobalMax(1, &max_expected, comm);
  Mpi::GlobalMax(1, &max_diff, comm);
  INFO(what << ": max. difference " << max_diff << " (max. indicator " << max_expected
            << ")");
  CHECK(max_expected > 0.0);
  CHECK(max_diff <= 1.0e-12 * max_expected);
}

double GlobalMaxEntry(MPI_Comm comm, const Vector &v)
{
  double max = v.Size() ? v.Max() : 0.0;
  Mpi::GlobalMax(1, &max, comm);
  return max;
}

constexpr const char *flux_size_refusal =
    "Recovered flux size does not match the flux recovery space";

}  // namespace

// The drivers with an interface "FluxRecovery" recover the electric flux D once and reuse
// it in the error estimate (SD-507 / SD-408): AddErrorIndicator(E, [B,] D) with the D of
// RecoverFlux / RecoverElectricFlux gives the indicator of AddErrorIndicator(E, [B]), on a
// recovery space which is unbroken (no interior boundary in the list) or broken across an
// interior sheet (SD-403: D lives on GetFluxSpace() / GetElectricFluxSpace()), for the
// gradient estimator (electrostatics) and the drivers' estimators (the time-dependent one
// in 3D, the boundary mode one in 2D). A D of another size is refused.
TEST_CASE("Reused recovered D gives the recomputed indicator",
          "[errorestimator][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  const int dim = GENERATE(2, 3);
  const bool broken = GENERATE(false, true);
  INFO((dim == 2 ? "2D" : "3D")
       << (broken ? ", sheet in the list" : ", no sheet in the list"));
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  SheetEstimatorSetup setup(comm, dim, order);
  auto &nd_fespace = setup.nd_fespaces.GetFinestFESpace();
  auto &rt_fespace = setup.rt_fespaces.GetFinestFESpace();
  const std::vector<int> crack_attr_list =
      broken ? std::vector<int>{1, sheet_attr} : std::vector<int>{1};
  constexpr double tol = 1.0e-14;
  constexpr int max_it = 10000, print = 0;
  constexpr bool use_mg = false;

  GradFluxErrorEstimator<Vector> grad(setup.mat_op, nd_fespace, setup.rt_fespaces, tol,
                                      max_it, print, use_mg, crack_attr_list);
  const auto &flux_space = grad.GetFluxSpace();
  CHECK(flux_space.IsBroken() == broken);
  if (broken)
  {
    CHECK(flux_space.GlobalTrueVSize() > rt_fespace.GlobalTrueVSize());
  }
  else
  {
    CHECK(&flux_space == &rt_fespace);
  }

  Vector E(nd_fespace.GetTrueVSize());
  E.UseDevice(true);
  FillRandom(E, 13);
  ErrorIndicator recomputed, reused;
  grad.AddErrorIndicator(E, 1.0, recomputed);
  Vector D(flux_space.GetTrueVSize());
  D.UseDevice(true);
  grad.RecoverFlux(E, D);
  grad.AddErrorIndicator(E, D, 1.0, reused);
  CheckSameIndicator(comm, reused.Local(), recomputed.Local(), "gradient estimator");

  Vector D_wrong(flux_space.GetTrueVSize() + 1);
  D_wrong.UseDevice(true);
  D_wrong = 0.0;
  ErrorIndicator unused;
  CHECK_THROWS_WITH(grad.RecoverFlux(E, D_wrong), ContainsSubstring(flux_size_refusal));
  CHECK_THROWS_WITH(grad.AddErrorIndicator(E, D_wrong, 1.0, unused),
                    ContainsSubstring(flux_size_refusal));

  // The drivers' form: complex E and B = curl E (on the RT space in 3D, on the L2 space in
  // 2D), D recovered by RecoverElectricFlux on GetElectricFluxSpace().
  auto CheckDriverForm = [&](const auto &estimator, const FiniteElementSpace &b_fespace)
  {
    const auto &d_space = estimator.GetElectricFluxSpace();
    CHECK(d_space.IsBroken() == broken);
    CHECK(d_space.GlobalTrueVSize() == flux_space.GlobalTrueVSize());
    ComplexVector Ec(nd_fespace.GetTrueVSize()), Bc(b_fespace.GetTrueVSize()),
        Dc(d_space.GetTrueVSize()), Dc_wrong(d_space.GetTrueVSize() + 1);
    Ec.UseDevice(true);
    Bc.UseDevice(true);
    Dc.UseDevice(true);
    Dc_wrong.UseDevice(true);
    FillRandom(Ec.Real(), 17);
    FillRandom(Ec.Imag(), 19);
    FillRandom(Bc.Real(), 23);
    FillRandom(Bc.Imag(), 29);
    Dc_wrong = 0.0;
    ErrorIndicator recomputed_c, reused_c, unused_c;
    estimator.AddErrorIndicator(Ec, Bc, 1.0, 1.0, recomputed_c);
    estimator.RecoverElectricFlux(Ec, Dc);
    estimator.AddErrorIndicator(Ec, Bc, Dc, 1.0, reused_c);
    CheckSameIndicator(comm, reused_c.Local(), recomputed_c.Local(), "drivers' estimator");
    CHECK_THROWS_WITH(estimator.RecoverElectricFlux(Ec, Dc_wrong),
                      ContainsSubstring(flux_size_refusal));
    CHECK_THROWS_WITH(estimator.AddErrorIndicator(Ec, Bc, Dc_wrong, 1.0, unused_c),
                      ContainsSubstring(flux_size_refusal));
  };
  if (dim == 3)
  {
    TimeDependentFluxErrorEstimator<ComplexVector> estimator(
        setup.mat_op, setup.nd_fespaces, setup.rt_fespaces, tol, max_it, print, use_mg,
        crack_attr_list);
    CheckDriverForm(estimator, rt_fespace);
  }
  else
  {
    BoundaryModeFluxErrorEstimator<ComplexVector> estimator(
        setup.mat_op, setup.nd_fespaces, setup.rt_fespaces, setup.l2_fespace,
        setup.h1_fespaces, tol, max_it, print, use_mg, crack_attr_list);
    CheckDriverForm(estimator, setup.l2_fespace);
  }
}

// On a mesh with an uncracked interior sheet the recovered electric flux lives on the
// recovery space broken across the sheet (SD-403) and carries the normal flux of EACH side:
// for a field whose normal component jumps across the sheet (with D = eps E in the
// recovery space on each side), the recovered D reproduces eps E on both sides, so that the
// indicator with this D vanishes, whereas a recovery space continuous across the sheet
// cannot follow the jump. The drivers' estimators recover the same D on
// GetElectricFluxSpace().
TEST_CASE("Recovered electric flux on the broken recovery space",
          "[errorestimator][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  const int dim = GENERATE(2, 3);
  INFO("dimension " << dim);
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  SheetEstimatorSetup setup(comm, dim, order);
  auto &nd_fespace = setup.nd_fespaces.GetFinestFESpace();
  constexpr double tol = 1.0e-14;
  constexpr int max_it = 10000, print = 0;
  constexpr bool use_mg = false;

  // E_x jumps across the sheet x = 0.5, and eps E_y = 1 is continuous across the dielectric
  // interface y = 0.5 (eps_r 4 below, 1 above): D = eps E is piecewise constant with a jump
  // of its normal component across the sheet only.
  mfem::VectorFunctionCoefficient E_coeff(dim,
                                          [](const mfem::Vector &x, mfem::Vector &v)
                                          {
                                            v(0) = (x(0) > 0.5) ? 2.0 : 1.0;
                                            v(1) = (x(1) > 0.5) ? 1.0 : 0.25;
                                            if (v.Size() > 2)
                                            {
                                              v(2) = -0.5;
                                            }
                                          });
  mfem::ParGridFunction E_gf(&nd_fespace.Get());
  E_gf.ProjectCoefficient(E_coeff);
  Vector E(nd_fespace.GetTrueVSize());
  E.UseDevice(true);
  E_gf.ParallelProject(E);

  auto Recover = [&](const std::vector<int> &crack_attr_list, bool expect_broken, Vector &D)
  {
    GradFluxErrorEstimator<Vector> grad(setup.mat_op, nd_fespace, setup.rt_fespaces, tol,
                                        max_it, print, use_mg, crack_attr_list);
    const auto &flux_space = grad.GetFluxSpace();
    CHECK(flux_space.IsBroken() == expect_broken);
    D.SetSize(flux_space.GetTrueVSize());
    D.UseDevice(true);
    grad.RecoverFlux(E, D);
    ErrorIndicator indicator;
    grad.AddErrorIndicator(E, D, 1.0, indicator);
    return std::make_pair(GlobalMaxEntry(comm, indicator.Local()),
                          flux_space.GlobalTrueVSize());
  };
  Vector D_continuous, D_broken;
  const auto [continuous_max, continuous_size] = Recover({}, false, D_continuous);
  const auto [broken_max, broken_size] = Recover({sheet_attr}, true, D_broken);
  INFO("max. indicator: continuous recovery " << continuous_max << ", broken across the "
                                              << "sheet " << broken_max);
  CHECK(broken_size > continuous_size);
  CHECK(continuous_max > 1.0e-3);
  CHECK(broken_max < 1.0e-8 * continuous_max);

  // The drivers' estimators (complex fields): the same D on GetElectricFluxSpace().
  auto CheckDriverForm = [&](const auto &estimator)
  {
    const auto &d_space = estimator.GetElectricFluxSpace();
    CHECK(d_space.IsBroken());
    REQUIRE(d_space.GetTrueVSize() == D_broken.Size());
    ComplexVector Ec(E.Size()), Dc(d_space.GetTrueVSize());
    Ec.UseDevice(true);
    Dc.UseDevice(true);
    Ec.Real() = E;
    Ec.Imag() = 0.0;
    estimator.RecoverElectricFlux(Ec, Dc);
    double max_d = D_broken.Size() ? D_broken.Normlinf() : 0.0, max_diff = 0.0,
           max_imag = Dc.Imag().Size() ? Dc.Imag().Normlinf() : 0.0;
    const auto *hd = D_broken.HostRead();
    const auto *hc = Dc.Real().HostRead();
    for (int i = 0; i < D_broken.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(hc[i] - hd[i]));
    }
    Mpi::GlobalMax(1, &max_d, comm);
    Mpi::GlobalMax(1, &max_diff, comm);
    Mpi::GlobalMax(1, &max_imag, comm);
    INFO("drivers' D: max. difference " << max_diff << ", max. imaginary part " << max_imag
                                        << " (max. |D| " << max_d << ")");
    CHECK(max_d > 0.0);
    CHECK(max_diff <= 1.0e-10 * max_d);
    CHECK(max_imag <= 1.0e-10 * max_d);
  };
  const std::vector<int> crack_attr_list = {sheet_attr};
  if (dim == 3)
  {
    CheckDriverForm(TimeDependentFluxErrorEstimator<ComplexVector>(
        setup.mat_op, setup.nd_fespaces, setup.rt_fespaces, tol, max_it, print, use_mg,
        crack_attr_list));
  }
  else
  {
    CheckDriverForm(BoundaryModeFluxErrorEstimator<ComplexVector>(
        setup.mat_op, setup.nd_fespaces, setup.rt_fespaces, setup.l2_fespace,
        setup.h1_fespaces, tol, max_it, print, use_mg, crack_attr_list));
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
