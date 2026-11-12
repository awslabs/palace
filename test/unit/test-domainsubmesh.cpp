// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
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

// Hex cube [0, 1]³ with n³ cells, attribute 1 for x < 0.5 and 2 for x > 0.5, optionally
// with all attribute 2 elements on process 0 and with nonconforming refinement on both
// sides of x = 0.5.
std::unique_ptr<mfem::ParMesh> MakeSplitCube(MPI_Comm comm, int n, bool attr2_on_root,
                                             bool refine)
{
  auto serial = mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::HEXAHEDRON);
  mfem::Vector c;
  for (int e = 0; e < serial.GetNE(); e++)
  {
    serial.GetElementCenter(e, c);
    serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
  }
  serial.SetAttributes();
  if (refine)
  {
    serial.EnsureNCMesh();
  }
  std::unique_ptr<int[]> partitioning;
  if (attr2_on_root)
  {
    const int size = Mpi::Size(comm);
    partitioning = std::make_unique<int[]>(serial.GetNE());
    for (int e = 0, k = 0; e < serial.GetNE(); e++)
    {
      partitioning[e] = (serial.GetAttribute(e) == 2) ? 0 : 1 + (k++ % (size - 1));
    }
  }
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial, partitioning.get());
  if (refine)
  {
    mfem::Array<int> marked;
    for (int e = 0; e < pmesh->GetNE(); e++)
    {
      pmesh->GetElementCenter(e, c);
      if ((c(0) > 0.3 && c(0) < 0.5 && c(1) < 0.5) ||
          (c(0) > 0.5 && c(0) < 0.7 && c(1) > 0.5))
      {
        marked.Append(e);
      }
    }
    pmesh->GeneralRefinement(marked, 1);
  }
  return pmesh;
}

// Domain submesh with libCEED data built from the submesh itself.
std::unique_ptr<Mesh> MakeDomainSubMesh(const mfem::ParMesh &parent,
                                        const mfem::Array<int> &attributes)
{
  auto mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromDomain(parent, attributes)));
  mesh->RebuildCeedAttributes();
  return mesh;
}

}  // namespace

TEST_CASE("Gradient flux error estimator on a domain submesh with an empty process",
          "[errorestimator][Parallel]")
{
  const auto comm = Mpi::World();
  if (Mpi::Size(comm) < 2)
  {
    SKIP("Needs at least two processes");
  }
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;

  // The global estimate must not change when a process has no submesh elements.
  auto Estimate = [&](bool attr2_on_root)
  {
    auto parent = MakeSplitCube(comm, 4, attr2_on_root, false);
    auto mesh = MakeDomainSubMesh(*parent, mfem::Array<int>({2}));
    if (attr2_on_root)
    {
      REQUIRE((mesh->GetNE() > 0) == (Mpi::Rank(comm) == 0));
    }

    config::MaterialData material;
    material.attributes = {2};
    material.epsilon_r.s = {2.0, 3.0, 4.0};
    config::PeriodicBoundaryData periodic;
    MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, *mesh);
    mfem::ND_FECollection nd_fec(order, 3);
    mfem::RT_FECollection rt_fec(order - 1, 3);
    FiniteElementSpace nd_fespace(*mesh, &nd_fec);
    FiniteElementSpaceHierarchy rt_fespaces(
        std::make_unique<FiniteElementSpace>(*mesh, &rt_fec));
    GradFluxErrorEstimator<Vector> estimator(mat_op, nd_fespace, rt_fespaces, 1.0e-14, 1000,
                                             0, false);

    mfem::VectorFunctionCoefficient f(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(0) + x(1));
                                        v(1) = std::cos(x(1) * x(2));
                                        v(2) = x(0) * x(2);
                                      });
    mfem::ParGridFunction gf(&nd_fespace.Get());
    gf.ProjectCoefficient(f);
    Vector E(nd_fespace.GetTrueVSize());
    E.UseDevice(true);
    gf.GetTrueDofs(E);
    ErrorIndicator indicator;
    estimator.AddErrorIndicator(E, 1.0, indicator);
    REQUIRE(indicator.Local().Size() == mesh->GetNE());
    return indicator.Norml2(comm);
  };

  const double ref = Estimate(false);
  REQUIRE(ref > 0.0);
  CHECK_THAT(Estimate(true), WithinRel(ref, 1.0e-8));
}

TEST_CASE("libCEED operators on a nonconforming domain submesh",
          "[libCEED][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  const int order = GENERATE(1, 2);
  fem::DefaultIntegrationOrder::p_trial = order;

  // Diffusion with a coefficient jump across a nonconforming interface, also between
  // processes, against MFEM's assembly.
  auto parent = MakeSplitCube(comm, 4, false, true);
  auto mesh = MakeDomainSubMesh(*parent, mfem::Array<int>({1, 2}));
  REQUIRE(mesh->Get().Nonconforming());
  if (Mpi::Size(comm) > 1)
  {
    REQUIRE(mesh->Get().GetNSharedFaces() > 0);
  }

  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace fespace(*mesh, &h1_fec);
  mfem::Vector coeff({1.0, 10.0});
  mfem::Array<int> attr_mat(mesh->MaxCeedAttribute());
  attr_mat = -1;
  mfem::DenseTensor mat_coeff(1, 1, 2);
  for (int k = 0; k < 2; k++)
  {
    mat_coeff(k)(0, 0) = coeff(k);
    for (auto attr : mesh->GetCeedAttributes(k + 1))
    {
      attr_mat[attr - 1] = k;
    }
  }
  MaterialPropertyCoefficient Q(attr_mat, mat_coeff);
  BilinearForm a_test(fespace);
  a_test.AddDomainIntegrator<DiffusionIntegrator>(Q);
  auto op_test = a_test.PartialAssemble();

  mfem::PWConstCoefficient Q_ref(coeff);
  mfem::ParBilinearForm a_ref(&fespace.Get());
  a_ref.AddDomainIntegrator(new mfem::DiffusionIntegrator(Q_ref));
  a_ref.Assemble();
  a_ref.Finalize();

  Vector x(op_test->Width()), y_test(op_test->Height()), y_ref(op_test->Height());
  x.UseDevice(true);
  y_test.UseDevice(true);
  y_ref.UseDevice(true);
  linalg::SetRandom(comm, x);
  op_test->Mult(x, y_test);
  a_ref.SpMat().Mult(x, y_ref);
  REQUIRE(linalg::Norml2(comm, y_ref) > 0.0);
  y_test -= y_ref;
  CHECK(linalg::Norml2(comm, y_test) < 1.0e-12 * linalg::Norml2(comm, y_ref));
}

}  // namespace palace
