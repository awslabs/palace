// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/shapederivativeform.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "test-helpers.hpp"

namespace palace
{

namespace
{

// Load a test mesh with alternating domain and boundary attributes. Its nodes are moved by
// a smooth map, so that the elements are not affine (and are curved for high-order nodes).
std::unique_ptr<Mesh> LoadMesh(MPI_Comm comm, const std::string &name)
{
  mfem::Mesh smesh(TestMeshPath(name), 1, 1);
  smesh.EnsureNodes();
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    smesh.SetAttribute(i, 1 + i % 2);
  }
  for (int i = 0; i < smesh.GetNBE(); i++)
  {
    smesh.SetBdrAttribute(i, 1 + i % 3);
  }
  smesh.SetAttributes();
  mfem::GridFunction &nodes = *smesh.GetNodes();
  const auto &nodal_fespace = *nodes.FESpace();
  const int dim = smesh.SpaceDimension();
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    double x[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < dim; d++)
    {
      x[d] = nodes(nodal_fespace.DofToVDof(i, d));
    }
    for (int d = 0; d < dim; d++)
    {
      nodes(nodal_fespace.DofToVDof(i, d)) +=
          0.02 * std::sin(1.3 * x[(d + 1) % dim] + 0.7 * x[(d + 2) % dim] + d);
    }
  }
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  return std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, smesh));
}

// A piecewise constant coefficient on the domain or boundary attributes of the mesh, with
// a different value for each attribute: scalar, or an anisotropic matrix (symmetric or
// not).
MaterialPropertyCoefficient BuildCoefficient(const Mesh &mesh, bool bdr, int dim,
                                             bool symmetric)
{
  const auto &attributes = bdr ? mesh.Get().bdr_attributes : mesh.Get().attributes;
  const int num_mat = attributes.Max();
  mfem::DenseTensor mat_coeff(dim, dim, num_mat);
  for (int k = 0; k < num_mat; k++)
  {
    for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
      {
        mat_coeff(i, j, k) =
            (i == j) ? 1.0 + 0.5 * k + 0.3 * i : 0.1 * (1 + k) * (symmetric ? 1.0 : i - j);
      }
    }
  }
  mfem::Array<int> attr_mat(bdr ? mesh.MaxCeedBdrAttribute() : mesh.MaxCeedAttribute());
  attr_mat = 0;
  for (int k = 0; k < num_mat; k++)
  {
    for (auto attr : bdr ? mesh.GetCeedBdrAttributes(k + 1) : mesh.GetCeedAttributes(k + 1))
    {
      attr_mat[attr - 1] = k;
    }
  }
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

Vector RandomVector(int size, std::mt19937 &gen)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Vector v(size);
  for (int i = 0; i < size; i++)
  {
    v(i) = dist(gen);
  }
  return v;
}

// The nodal displacement V(X) = A X + b on the true DOFs of the mesh nodal space.
Vector AffineDisplacement(const mfem::ParFiniteElementSpace &nodal_fespace,
                          const mfem::GridFunction &nodes, const mfem::DenseMatrix &A,
                          const mfem::Vector &b)
{
  const int dim = nodal_fespace.GetVDim();
  mfem::Vector V_l(nodal_fespace.GetVSize());
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    for (int d = 0; d < dim; d++)
    {
      double v = b(d);
      for (int e = 0; e < dim; e++)
      {
        v += A(d, e) * nodes(nodal_fespace.DofToVDof(i, e));
      }
      V_l(nodal_fespace.DofToVDof(i, d)) = v;
    }
  }
  Vector V(nodal_fespace.GetTrueVSize());
  nodal_fespace.GetRestrictionMatrix()->Mult(V_l, V);
  return V;
}

// The form a(u, w) = wᵀ A u, for true DOF vectors.
double EvaluateForm(const BilinearForm &a, const FiniteElementSpace &fespace,
                    const Vector &u, const Vector &w)
{
  Vector u_l(fespace.GetVSize()), w_l(fespace.GetVSize()), y_l(fespace.GetVSize());
  fespace.GetProlongationMatrix()->Mult(u, u_l);
  fespace.GetProlongationMatrix()->Mult(w, w_l);
  a.PartialAssemble()->Mult(u_l, y_l);
  return linalg::Dot<Vector>(fespace.GetComm(), w_l, y_l);
}

// Move the mesh nodes by h V, for V on the true DOFs of the mesh nodal space.
void MoveNodes(Mesh &mesh, const mfem::ParFiniteElementSpace &nodal_fespace,
               const Vector &V, double h)
{
  mfem::Vector V_l(nodal_fespace.GetVSize());
  nodal_fespace.GetProlongationMatrix()->Mult(V, V_l);
  mesh.Get().GetNodes()->Add(h, V_l);
  mesh.Get().NodesUpdated();
  mesh.ResetCeedObjects();
}

enum class FormType
{
  Diffusion,
  DiffusionIsotropic,
  DiffusionNonsymmetric,
  Mass,
  BoundaryMass
};

void TestShapeDerivative(const std::string &mesh_name, int order, FormType type)
{
  auto mesh = LoadMesh(Mpi::World(), mesh_name);
  const int dim = mesh->Dimension();
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  mfem::H1_FECollection fec(order, dim);
  FiniteElementSpace fespace(*mesh, &fec);
  const bool bdr = (type == FormType::BoundaryMass);
  const bool diffusion =
      (type == FormType::Diffusion || type == FormType::DiffusionIsotropic ||
       type == FormType::DiffusionNonsymmetric);
  auto Q = BuildCoefficient(*mesh, bdr, diffusion ? dim : 1,
                            type != FormType::DiffusionNonsymmetric);
  BilinearForm a(fespace);
  ShapeDerivativeForm da(fespace);
  switch (type)
  {
    case FormType::Diffusion:
    case FormType::DiffusionNonsymmetric:
      a.AddDomainIntegrator<DiffusionIntegrator>(Q);
      da.AddDomainIntegrator<DiffusionShapeIntegrator>(Q);
      break;
    case FormType::DiffusionIsotropic:
      a.AddDomainIntegrator<DiffusionIntegrator>();
      da.AddDomainIntegrator<DiffusionShapeIntegrator>();
      break;
    case FormType::Mass:
      a.AddDomainIntegrator<MassIntegrator>(Q);
      da.AddDomainIntegrator<MassShapeIntegrator>(Q);
      break;
    case FormType::BoundaryMass:
      a.AddBoundaryIntegrator<MassIntegrator>(Q);
      da.AddBoundaryIntegrator<MassShapeIntegrator>(Q);
      break;
  }

  std::mt19937 gen(1234 + Mpi::Rank(Mpi::World()));
  const Vector u = RandomVector(fespace.GetTrueVSize(), gen);
  const Vector w = RandomVector(fespace.GetTrueVSize(), gen);
  const auto da_op = da.Assemble();
  const auto &nodal_fespace = da_op->GetNodalSpace();
  Vector g;
  da_op->Mult(u, w, g);
  REQUIRE(g.Size() == nodal_fespace.GetTrueVSize());
  const double a0 = EvaluateForm(a, fespace, u, w);
  const auto &nodes = *mesh->Get().GetNodes();
  auto Contract = [&](const Vector &V) { return linalg::Dot<Vector>(Mpi::World(), g, V); };
  double g_scale = 0.0;
  for (int i = 0; i < g.Size(); i++)
  {
    g_scale += std::abs(g(i));
  }
  Mpi::GlobalSum(1, &g_scale, Mpi::World());

  // Translations leave the form unchanged, and a dilation by (1 + t) scales it by
  // (1 + t)^n with n = dim - 2 for diffusion, dim for domain mass, dim - 1 for boundary
  // mass. Both are exact for the discrete form.
  mfem::DenseMatrix A(dim), I(dim);
  mfem::Vector b(dim), zero(dim);
  A = 0.0;
  I = 0.0;
  zero = 0.0;
  for (int d = 0; d < dim; d++)
  {
    I(d, d) = 1.0;
    b = 0.0;
    b(d) = 1.0;
    CHECK(std::abs(Contract(AffineDisplacement(nodal_fespace, nodes, A, b))) <=
          1.0e-12 * g_scale);
  }
  const int n = diffusion ? dim - 2 : (bdr ? dim - 1 : dim);
  CHECK(std::abs(Contract(AffineDisplacement(nodal_fespace, nodes, I, zero)) - n * a0) <=
        1.0e-11 * (std::abs(a0) + g_scale));

  // Rotations leave forms with isotropic coefficients unchanged.
  if (type != FormType::Diffusion && type != FormType::DiffusionNonsymmetric)
  {
    mfem::DenseMatrix W(dim);
    W = 0.0;
    W(0, 1) = 1.0;
    W(1, 0) = -1.0;
    if (dim == 3)
    {
      W(0, 2) = 0.5;
      W(2, 0) = -0.5;
      W(1, 2) = -0.3;
      W(2, 1) = 0.3;
    }
    CHECK(std::abs(Contract(AffineDisplacement(nodal_fespace, nodes, W, zero))) <=
          1.0e-11 * g_scale);
  }

  // A random displacement, against central finite differences of the assembled form on
  // moved meshes, with Richardson extrapolation.
  const Vector V = RandomVector(nodal_fespace.GetTrueVSize(), gen);
  auto CentralDifference = [&](double h)
  {
    MoveNodes(*mesh, nodal_fespace, V, h);
    const double ap = EvaluateForm(a, fespace, u, w);
    MoveNodes(*mesh, nodal_fespace, V, -2.0 * h);
    const double am = EvaluateForm(a, fespace, u, w);
    MoveNodes(*mesh, nodal_fespace, V, h);
    return (ap - am) / (2.0 * h);
  };
  const double h = 1.0e-3;
  const double fd = (4.0 * CentralDifference(0.5 * h) - CentralDifference(h)) / 3.0;
  const double dv = Contract(V);
  INFO("a = " << a0 << ", g ⋅ V = " << dv << ", finite difference = " << fd);
  CHECK(std::abs(dv - fd) <= 1.0e-7 * std::abs(fd));
}

}  // namespace

TEST_CASE("Shape derivatives of H1 forms", "[shapederivative][Serial][Parallel]")
{
  auto mesh =
      GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh", "fichera-hex.mesh",
               "fichera-tet.mesh", "fichera-mixed-p2.mesh", "tinyzoo-3d.mesh");
  auto order = GENERATE(1, 2, 3);
  auto type =
      GENERATE(FormType::Diffusion, FormType::DiffusionIsotropic,
               FormType::DiffusionNonsymmetric, FormType::Mass, FormType::BoundaryMass);
  INFO("Mesh: " << mesh << ", order: " << order << ", form: " << static_cast<int>(type));
  TestShapeDerivative(mesh, order, type);
}

}  // namespace palace
