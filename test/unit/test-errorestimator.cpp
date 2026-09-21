// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
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

// Releasing the workspace may not change the computed indicator at all, since none of the
// freed storage carries state between estimates.
void CheckExactlyEqual(const Vector &x, const Vector &x_ref)
{
  REQUIRE(x.Size() == x_ref.Size());
  const auto *px = x.HostRead();
  const auto *px_ref = x_ref.HostRead();
  bool equal = true;
  for (int i = 0; i < x.Size(); i++)
  {
    equal = equal && (px[i] == px_ref[i]);
  }
  CHECK(equal);
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
  CheckExactlyEqual(d.Real(), d_ref.Real());
  CheckExactlyEqual(d.Imag(), d_ref.Imag());
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
  CheckExactlyEqual(indicator.Local(), indicator_ref.Local());
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
  CheckExactlyEqual(indicator.Local(), indicator_ref.Local());
}

// The driven solver releases the workspace of the combined estimator after every frequency
// or PROM sample, so repeated estimates of the same fields have to agree exactly.
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
  CheckExactlyEqual(indicator.Local(), indicator_ref.Local());

  // A second release without an intervening estimate is a no-op.
  estimator.ReleaseWorkspace();
  estimator.ReleaseWorkspace();
  ErrorIndicator indicator_again;
  estimator.AddErrorIndicator(e, b, 1.0, indicator_again);
  CheckExactlyEqual(indicator_again.Local(), indicator_ref.Local());
}

// A mesh with more than one element geometry gives the error integration operator several
// sub-operators, each with its own passive input vectors referencing the grid function work
// vectors. Only the first sub-operator's are re-pointed at each estimate, so those vectors
// stay allocated while the rest of the workspace is released.
TEST_CASE("Mixed mesh flux estimator workspace release keeps the passive inputs",
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
  CHECK(estimator.WorkVectorSizes() == std::array<int, 3>{sizes[0], 0, sizes[2]});

  estimator.AddErrorIndicator(e, 1.0, indicator);
  CHECK(estimator.WorkVectorSizes() == sizes);
  CheckExactlyEqual(indicator.Local(), indicator_ref.Local());
}

}  // namespace palace
