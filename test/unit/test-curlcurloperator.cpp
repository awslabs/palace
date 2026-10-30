// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

using namespace Catch::Matchers;

// A uniform field A = a has no curl, so AᵀK·A is the London energy |a|² V_L/(μλ²) of the
// London domain alone, as is AᵀM_L·A.
TEST_CASE("CurlCurlOperator London domains", "[curlcurloperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(Units(1.0e-6, 1.0e-6));
  iodata.problem.type = ProblemType::MAGNETOSTATIC;
  const double mu2 = 2.0, lambda = 0.3;  // London domain, in mesh length units
  {
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {1};
  }
  {
    auto &material = iodata.domains.materials.emplace_back();
    material.attributes = {2};
    material.mu_r = mu2;
    material.lambda_L = lambda;
  }
  iodata.solver.order = 2;
  iodata.CheckConfiguration();

  // Unit cube split at x = 1/2 into a normal (1) and a London (2) domain.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 2, 2, mfem::Element::HEXAHEDRON));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(i, center);
    serial_mesh->SetAttribute(i, center(0) < 0.5 ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, comm);
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *serial_mesh)));

  // Nondimensional domain volumes and London depth.
  double volume[2] = {0.0, 0.0};
  auto &pmesh = mesh.back()->Get();
  for (int i = 0; i < pmesh.GetNE(); i++)
  {
    volume[pmesh.GetAttribute(i) - 1] += pmesh.GetElementVolume(i);
  }
  Mpi::GlobalSum(2, volume, comm);
  const double lambda_nd = iodata.domains.materials[1].lambda_L;

  CurlCurlOperator curlcurl_op(iodata, mesh);
  REQUIRE(curlcurl_op.HasLondonDomains());

  // Interpolate a uniform field, exactly representable in the Nédélec space.
  mfem::Vector a(3);
  a(0) = 1.0;
  a(1) = -2.0;
  a(2) = 3.0;
  mfem::VectorConstantCoefficient a_coeff(a);
  GridFunction A_gf(curlcurl_op.GetNDSpace());
  A_gf.Real().ProjectCoefficient(a_coeff);
  Vector A(curlcurl_op.GetNDSpace().GetTrueVSize()), y(A.Size());
  A.UseDevice(true);
  y.UseDevice(true);
  A_gf.Real().GetTrueDofs(A);
  const double london = (a * a) * volume[1] / (mu2 * lambda_nd * lambda_nd);

  auto K = curlcurl_op.GetStiffnessMatrix();
  K->Mult(A, y);
  CHECK_THAT(linalg::Dot(comm, A, y), WithinRel(london, 1.0e-10));

  curlcurl_op.ApplyLondonMass(A, y);
  CHECK_THAT(linalg::Dot(comm, A, y), WithinRel(london, 1.0e-10));

  // The preconditioner adds the gauge shift · (1/μ) ∫|A|² over the whole domain.
  auto P = curlcurl_op.GetPreconditionerMatrix();
  REQUIRE(P);
  P->Mult(A, y);
  const double shift =
      iodata.solver.linear.london_pc_shift * (a * a) * (volume[0] + volume[1] / mu2);
  CHECK_THAT(linalg::Dot(comm, A, y), WithinRel(london + shift, 1.0e-10));
}

}  // namespace palace
