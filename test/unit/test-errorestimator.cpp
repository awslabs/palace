// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

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

namespace
{

// Expose the work vectors, to check that ReleaseWorkspace actually frees them and that the
// next estimate reallocates them.
template <typename VecType>
class ProbeFluxProjector : public FluxProjector<VecType>
{
public:
  using FluxProjector<VecType>::FluxProjector;

  int RhsSize() const { return this->rhs.Size(); }
};

template <typename VecType>
class ProbeGradFluxErrorEstimator : public GradFluxErrorEstimator<VecType>
{
public:
  using GradFluxErrorEstimator<VecType>::GradFluxErrorEstimator;

  std::array<int, 3> WorkVectorSizes() const
  {
    return {this->E_gf.Size(), this->D.Size(), this->D_gf.Size()};
  }
};

template <typename VecType>
class ProbeCurlFluxErrorEstimator : public CurlFluxErrorEstimator<VecType>
{
public:
  using CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator;

  std::array<int, 3> WorkVectorSizes() const
  {
    return {this->B_gf.Size(), this->H.Size(), this->H_gf.Size()};
  }
};

// Small driven-problem setup on a Cartesian tetrahedral mesh, providing the material
// operator and the ND/RT space hierarchies the estimators are built on.
class DrivenFixture
{
public:
  Units units{0.496, 1.453};
  IoData iodata{units};
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<SpaceOperator> space_op;

  DrivenFixture()
  {
    iodata.domains.materials.emplace_back().attributes = {1};
    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.driven.sample_f = {1.0};
    iodata.solver.linear.estimator_tol = 1.0e-10;
    iodata.solver.linear.estimator_max_it = 100;
    const nlohmann::json boundaries = {{"LumpedPort",
                                        {{{"Attributes", {2}},
                                          {"Index", 1},
                                          {"R", 50.0},
                                          {"Direction", "+X"},
                                          {"Excitation", true}}}}};
    iodata.boundaries.lumpedport = config::BoundaryData(boundaries).lumpedport;
    iodata.CheckConfiguration();  // Initializes the quadrature orders

    constexpr int resolution = 2;
    auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
        resolution, resolution, resolution, mfem::Element::TETRAHEDRON));
    const auto comm = Mpi::World();
    iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, comm);
    iodata.NondimensionalizeInputs(serial_mesh);
    mesh.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));
    space_op = std::make_unique<SpaceOperator>(iodata, mesh);
  }

  // Random electric field and the matching magnetic flux density on the true dofs.
  void GetFields(ComplexVector &e, ComplexVector &b) const
  {
    const auto &Curl = space_op->GetCurlMatrix();
    e.SetSize(Curl.Width());
    b.SetSize(Curl.Height());
    e.UseDevice(true);
    b.UseDevice(true);
    e.Real().Randomize(1 + Mpi::Rank(Mpi::World()));
    e.Imag().Randomize(101 + Mpi::Rank(Mpi::World()));
    Curl.Mult(e.Real(), b.Real());
    Curl.Mult(e.Imag(), b.Imag());
  }
};

// Releasing the workspace may only change the result at roundoff level: GPU reductions are
// not required to be bitwise reproducible after reallocating their work arrays.
void CheckClose(const Vector &x, const Vector &x_ref)
{
  REQUIRE(x.Size() == x_ref.Size());
  Vector diff(x);
  diff -= x_ref;
  const double norm_ref = linalg::Norml2(Mpi::World(), x_ref);
  const double norm_diff = linalg::Norml2(Mpi::World(), diff);
  CAPTURE(norm_ref, norm_diff);
  CHECK(norm_diff <= 1.0e-13 + 1.0e-12 * norm_ref);
}

std::unique_ptr<Mesh> LoadTestMesh(MPI_Comm comm, const std::string &file)
{
  mfem::Mesh smesh(std::string(PALACE_TEST_DATA_DIR) + "/mesh/" + file, 1, 1);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  return std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, smesh));
}

}  // namespace

TEST_CASE("Flux projector workspace release reproduces the projection",
          "[errorestimator][Serial][Parallel]")
{
  DrivenFixture fix;
  const auto &mat_op = fix.space_op->GetMaterialOp();
  ProbeFluxProjector<ComplexVector> projector(
      MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetPermittivityReal()),
      fix.space_op->GetRTSpaces(), fix.space_op->GetNDSpace(), 1.0e-12, 100, 0, false);

  ComplexVector e, b;
  fix.GetFields(e, b);
  ComplexVector d_ref(b.Size()), d(b.Size());
  d_ref.UseDevice(true);
  d.UseDevice(true);
  projector.Mult(e, d_ref);
  const int rhs_size = projector.RhsSize();
  REQUIRE(rhs_size > 0);

  projector.ReleaseWorkspace();
  CHECK(projector.RhsSize() == 0);

  projector.Mult(e, d);
  CHECK(projector.RhsSize() == rhs_size);
  CheckClose(d.Real(), d_ref.Real());
  CheckClose(d.Imag(), d_ref.Imag());
}

TEST_CASE("Gradient flux estimator workspace release reproduces the indicator",
          "[errorestimator][Serial][Parallel]")
{
  DrivenFixture fix;
  ProbeGradFluxErrorEstimator<ComplexVector> estimator(
      fix.space_op->GetMaterialOp(), fix.space_op->GetNDSpace(),
      fix.space_op->GetRTSpaces(), fix.iodata.solver.linear.estimator_tol,
      fix.iodata.solver.linear.estimator_max_it, 0, false);

  ComplexVector e, b;
  fix.GetFields(e, b);
  ErrorIndicator indicator_ref, indicator;
  estimator.AddErrorIndicator(e, 1.0, indicator_ref);
  const auto sizes = estimator.WorkVectorSizes();
  REQUIRE(sizes[0] > 0);
  REQUIRE(sizes[1] > 0);
  REQUIRE(sizes[2] > 0);

  estimator.ReleaseWorkspace();
  CHECK(estimator.WorkVectorSizes() == std::array<int, 3>{0, 0, 0});

  estimator.AddErrorIndicator(e, 1.0, indicator);
  CHECK(estimator.WorkVectorSizes() == sizes);
  CheckClose(indicator.Local(), indicator_ref.Local());
}

TEST_CASE("Curl flux estimator workspace release reproduces the indicator",
          "[errorestimator][Serial][Parallel]")
{
  DrivenFixture fix;
  ProbeCurlFluxErrorEstimator<ComplexVector> estimator(
      fix.space_op->GetMaterialOp(), fix.space_op->GetCurlSpace(),
      fix.space_op->GetNDSpaces(), fix.iodata.solver.linear.estimator_tol,
      fix.iodata.solver.linear.estimator_max_it, 0, false);

  ComplexVector e, b;
  fix.GetFields(e, b);
  ErrorIndicator indicator_ref, indicator;
  estimator.AddErrorIndicator(b, 1.0, indicator_ref);
  const auto sizes = estimator.WorkVectorSizes();
  REQUIRE(sizes[0] > 0);
  REQUIRE(sizes[1] > 0);
  REQUIRE(sizes[2] > 0);

  estimator.ReleaseWorkspace();
  CHECK(estimator.WorkVectorSizes() == std::array<int, 3>{0, 0, 0});

  estimator.AddErrorIndicator(b, 1.0, indicator);
  CHECK(estimator.WorkVectorSizes() == sizes);
  CheckClose(indicator.Local(), indicator_ref.Local());
}

// The driven solver releases the workspace of the combined estimator after every frequency
// or PROM sample, so repeated estimates of the same fields must agree to roundoff.
TEST_CASE("Combined flux estimator workspace release reproduces the indicator",
          "[errorestimator][Serial][Parallel]")
{
  DrivenFixture fix;
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      fix.space_op->GetMaterialOp(), fix.space_op->GetNDSpaces(),
      fix.space_op->GetRTSpaces(), fix.iodata.solver.linear.estimator_tol,
      fix.iodata.solver.linear.estimator_max_it, 0, false);

  ComplexVector e, b;
  fix.GetFields(e, b);
  ErrorIndicator indicator_ref, indicator;
  estimator.AddErrorIndicator(e, b, 1.0, indicator_ref);
  estimator.ReleaseWorkspace();
  estimator.AddErrorIndicator(e, b, 1.0, indicator);
  CheckClose(indicator.Local(), indicator_ref.Local());

  // A second release without an intervening estimate is a no-op.
  estimator.ReleaseWorkspace();
  estimator.ReleaseWorkspace();
  ErrorIndicator indicator_again;
  estimator.AddErrorIndicator(e, b, 1.0, indicator_again);
  CheckClose(indicator_again.Local(), indicator_ref.Local());
}

// Every element geometry sub-operator shares one passive input-vector pair per libCEED
// context, so a mixed mesh can release and reallocate the grid function work vectors too.
TEST_CASE("Mixed mesh flux estimator workspace release reproduces the indicator",
          "[errorestimator][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  auto mesh = LoadTestMesh(comm, "fichera-mixed-p2.mesh");
  REQUIRE(mesh->Get().GetNumGeometries(mesh->Dimension()) > 1);

  config::MaterialData material;
  material.attributes = {1};
  const std::vector<config::MaterialData> materials = {material};
  const config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(1, mesh->Dimension());
  mfem::RT_FECollection rt_fec(0, mesh->Dimension());
  FiniteElementSpaceHierarchy nd_fespaces(
      std::make_unique<FiniteElementSpace>(*mesh, &nd_fec));
  FiniteElementSpaceHierarchy rt_fespaces(
      std::make_unique<FiniteElementSpace>(*mesh, &rt_fec));
  ProbeGradFluxErrorEstimator<ComplexVector> estimator(
      mat_op, nd_fespaces.GetFinestFESpace(), rt_fespaces, 1.0e-10, 100, 0, false);

  ComplexVector e(nd_fespaces.GetFinestFESpace().GetTrueVSize());
  e.UseDevice(true);
  e.Real().Randomize(1 + Mpi::Rank(comm));
  e.Imag().Randomize(101 + Mpi::Rank(comm));
  ErrorIndicator indicator_ref, indicator;
  estimator.AddErrorIndicator(e, 1.0, indicator_ref);
  const auto sizes = estimator.WorkVectorSizes();

  estimator.ReleaseWorkspace();
  CHECK(estimator.WorkVectorSizes() == std::array<int, 3>{0, 0, 0});

  estimator.AddErrorIndicator(e, 1.0, indicator);
  CHECK(estimator.WorkVectorSizes() == sizes);
  CheckClose(indicator.Local(), indicator_ref.Local());
}
}  // namespace palace
