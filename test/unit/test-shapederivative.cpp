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
#include "fem/coefficient.hpp"
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
// With amr, it is refined nonconformally, with hanging nodes.
std::unique_ptr<Mesh> LoadMesh(MPI_Comm comm, const std::string &name, bool amr)
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
  if (amr)
  {
    smesh.EnsureNCMesh(true);
  }
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  if (amr)
  {
    pmesh->RandomRefinement(0.5);
  }
  return std::make_unique<Mesh>(std::move(pmesh));
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

// The linear form b(w) = ∫ f ⋅ w ds on the boundary, with f constant on each attribute.
double EvaluateLinearForm(FiniteElementSpace &fespace, mfem::VectorCoefficient &f,
                          const Vector &w)
{
  mfem::LinearForm b(&fespace.Get());
  b.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(f));
  b.UseFastAssembly(false);
  b.Assemble();
  Vector w_l(fespace.GetVSize());
  fespace.GetProlongationMatrix()->Mult(w, w_l);
  return linalg::Dot<Vector>(fespace.GetComm(), w_l, b);
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
  BoundaryMass,
  VectorMass,
  VectorMassNonsymmetric,
  BoundaryVectorMass,
  CurlCurl,
  BoundaryCurlCurl,
  BoundaryLinearForm
};

bool IsVectorForm(FormType type)
{
  return type >= FormType::VectorMass;
}

bool IsBoundaryForm(FormType type)
{
  return type == FormType::BoundaryMass || type == FormType::BoundaryVectorMass ||
         type == FormType::BoundaryCurlCurl || type == FormType::BoundaryLinearForm;
}

// The exponent n of the scaling (1 + t)^n of the form under a dilation by (1 + t): dim - 2
// for diffusion, dim for H1 mass, dim - 2 for H(curl) mass (covariant fields), dim - 4 for
// curl-curl (contravariant curls), one less on boundaries, and dim - 2 for the boundary
// linear form with a fixed f.
int DilationExponent(FormType type, int dim)
{
  switch (type)
  {
    case FormType::Diffusion:
    case FormType::DiffusionIsotropic:
    case FormType::DiffusionNonsymmetric:
    case FormType::VectorMass:
    case FormType::VectorMassNonsymmetric:
    case FormType::BoundaryLinearForm:
      return dim - 2;
    case FormType::Mass:
      return dim;
    case FormType::BoundaryMass:
      return dim - 1;
    case FormType::BoundaryVectorMass:
      return dim - 3;
    case FormType::CurlCurl:
      return dim - 4;
    case FormType::BoundaryCurlCurl:
      return dim - 5;
  }
  return 0;
}

void TestShapeDerivative(const std::string &mesh_name, bool amr, int order, FormType type)
{
  auto mesh = LoadMesh(Mpi::World(), mesh_name, amr);
  const int dim = mesh->Dimension();
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  std::unique_ptr<mfem::FiniteElementCollection> fec;
  if (IsVectorForm(type))
  {
    fec = std::make_unique<mfem::ND_FECollection>(order, dim);
  }
  else
  {
    fec = std::make_unique<mfem::H1_FECollection>(order, dim);
  }
  FiniteElementSpace fespace(*mesh, fec.get());
  const bool bdr = IsBoundaryForm(type);
  const bool symmetric =
      (type != FormType::DiffusionNonsymmetric && type != FormType::VectorMassNonsymmetric);
  int coeff_dim = 1;
  if (type == FormType::Diffusion || type == FormType::DiffusionNonsymmetric ||
      type == FormType::VectorMass || type == FormType::VectorMassNonsymmetric ||
      type == FormType::BoundaryVectorMass || type == FormType::BoundaryLinearForm ||
      (type == FormType::CurlCurl && dim == 3))
  {
    coeff_dim = dim;
  }
  auto Q = BuildCoefficient(*mesh, bdr, coeff_dim, symmetric);
  BilinearForm a(fespace);
  ShapeDerivativeForm da(fespace);
  SumVectorCoefficient f(dim);
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
    case FormType::VectorMass:
    case FormType::VectorMassNonsymmetric:
      a.AddDomainIntegrator<VectorFEMassIntegrator>(Q);
      da.AddDomainIntegrator<VectorFEMassShapeIntegrator>(Q);
      break;
    case FormType::BoundaryVectorMass:
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(Q);
      da.AddBoundaryIntegrator<VectorFEMassShapeIntegrator>(Q);
      break;
    case FormType::CurlCurl:
      a.AddDomainIntegrator<CurlCurlIntegrator>(Q);
      da.AddDomainIntegrator<CurlCurlShapeIntegrator>(Q);
      break;
    case FormType::BoundaryCurlCurl:
      a.AddBoundaryIntegrator<CurlCurlIntegrator>(Q);
      da.AddBoundaryIntegrator<CurlCurlShapeIntegrator>(Q);
      break;
    case FormType::BoundaryLinearForm:
      // The vector of each attribute is the first column of its coefficient.
      for (int k = 0; k < mesh->Get().bdr_attributes.Max(); k++)
      {
        mfem::Array<int> attr_list(1);
        attr_list[0] = k + 1;
        mfem::Vector f_k(dim);
        for (int d = 0; d < dim; d++)
        {
          f_k(d) = Q.GetMaterialProperties()(d, 0, k);
        }
        f.AddCoefficient(
            std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
                attr_list, f_k));
      }
      da.AddBoundaryIntegrator<VectorFEBoundaryLFShapeIntegrator>(Q);
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
  auto Evaluate = [&]()
  {
    return (type == FormType::BoundaryLinearForm) ? EvaluateLinearForm(fespace, f, w)
                                                  : EvaluateForm(a, fespace, u, w);
  };
  const double a0 = Evaluate();
  const auto &nodes = *mesh->Get().GetNodes();
  auto Contract = [&](const Vector &V) { return linalg::Dot<Vector>(Mpi::World(), g, V); };
  double g_scale = 0.0;
  for (int i = 0; i < g.Size(); i++)
  {
    g_scale += std::abs(g(i));
  }
  Mpi::GlobalSum(1, &g_scale, Mpi::World());

  // Translations leave the form unchanged, and a dilation by (1 + t) scales it by
  // (1 + t)^n. Both are exact for the discrete form.
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
  const int n = DilationExponent(type, dim);
  CHECK(std::abs(Contract(AffineDisplacement(nodal_fespace, nodes, I, zero)) - n * a0) <=
        1.0e-11 * (std::abs(a0) + g_scale));

  // Rotations leave forms with isotropic coefficients unchanged.
  if (coeff_dim == 1 || type == FormType::DiffusionIsotropic)
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
    const double ap = Evaluate();
    MoveNodes(*mesh, nodal_fespace, V, -2.0 * h);
    const double am = Evaluate();
    MoveNodes(*mesh, nodal_fespace, V, h);
    return (ap - am) / (2.0 * h);
  };
  const double h = 2.5e-4;
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
  auto amr = GENERATE(false, true);
  if (amr && std::string(mesh) == "tinyzoo-3d.mesh")
  {
    return;  // The forward forms are NaN on nonconforming refinements of its pyramid
  }
  auto order = GENERATE(1, 2, 3);
  auto type =
      GENERATE(FormType::Diffusion, FormType::DiffusionIsotropic,
               FormType::DiffusionNonsymmetric, FormType::Mass, FormType::BoundaryMass);
  INFO("Mesh: " << mesh << ", AMR: " << amr << ", order: " << order
                << ", form: " << static_cast<int>(type));
  TestShapeDerivative(mesh, amr, order, type);
}

TEST_CASE("Shape derivatives of H(curl) forms", "[shapederivative][Serial][Parallel]")
{
  auto mesh =
      GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh", "fichera-hex.mesh",
               "fichera-tet.mesh", "fichera-mixed-p2.mesh", "tinyzoo-3d.mesh");
  auto amr = GENERATE(false, true);
  if (amr && std::string(mesh) == "tinyzoo-3d.mesh")
  {
    return;  // The forward forms are NaN on nonconforming refinements of its pyramid
  }
  auto order = GENERATE(1, 2, 3);
  auto type = GENERATE(FormType::VectorMass, FormType::VectorMassNonsymmetric,
                       FormType::BoundaryVectorMass, FormType::CurlCurl,
                       FormType::BoundaryCurlCurl, FormType::BoundaryLinearForm);
  if (std::string(mesh).starts_with("star") && type == FormType::BoundaryCurlCurl)
  {
    return;  // The curl of a trace is only defined on the boundaries of 3D meshes
  }
  INFO("Mesh: " << mesh << ", AMR: " << amr << ", order: " << order
                << ", form: " << static_cast<int>(type));
  TestShapeDerivative(mesh, amr, order, type);
}

}  // namespace palace
