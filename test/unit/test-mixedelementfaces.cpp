// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Palace's libCEED operators on meshes of a few elements of different types sharing faces,
// at every orientation-preserving vertex ordering of two of the elements (every third one
// at order 4), in serial and on two ranks: the discrete gradient (H1 -> H(curl)), the
// discrete curl (H(curl) -> H(div)) and the order prolongation of H(curl) reproduce fields
// of the spaces' degrees exactly; the element DOFs of the H(curl) restrictions equal MFEM's
// and the field evaluated from them is tangentially continuous; the H(curl) mass and
// curl-curl and the H(div) mass operators equal MFEM's element assembly. From order 2, the
// H(curl) DOFs on a triangular face need the element's DofTransformation, which the
// restrictions must take from the space (MFEM's pyramid does not return it to the
// collection).

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

struct MixedElement
{
  mfem::Geometry::Type geom;
  std::vector<int> vertices;
};

// Elements in mesh order, their ranks when run on two ranks, and the two elements whose
// vertex orderings are varied (the others keep theirs).
struct MixedCase
{
  std::string name;
  std::vector<std::array<double, 3>> coords;
  std::vector<MixedElement> elements;
  std::vector<int> ranks;
  int vary0, vary1;
};

// A pyramid and a prism sharing triangles with tetrahedra; a pyramid and a prism sharing a
// quadrilateral with a hexahedron (and a triangle with a tetrahedron); a pyramid sharing a
// triangle with a prism (and one with a tetrahedron); a tetrahedron | tetrahedron control.
// On two ranks one shared face is on the rank boundary and the pyramid or prism shares a
// second face with an element of its own rank which precedes it, so its local face
// orientation is not the identity there either. The coordinates are generic.
std::vector<MixedCase> MixedCases()
{
  using G = mfem::Geometry;
  return {
      {"tet|pyramid|tet",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {1.1, 1.0, 0.05},
        {-0.05, 0.9, 0.0},
        {0.45, 0.4, 0.95},
        {1.6, 0.7, 0.8},
        {0.5, -0.8, 0.4}},
       {{G::TETRAHEDRON, {1, 2, 4, 5}},
        {G::TETRAHEDRON, {0, 1, 4, 6}},
        {G::PYRAMID, {0, 1, 2, 3, 4}}},
       {0, 1, 1},
       2,
       0},
      {"tet|prism|tet",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {0.1, 0.9, 0.05},
        {0.05, 0.0, 1.0},
        {1.05, 0.1, 1.1},
        {0.0, 1.0, 0.95},
        {0.35, 0.4, 1.9},
        {0.3, 0.35, -0.9}},
       {{G::TETRAHEDRON, {3, 4, 5, 6}},
        {G::TETRAHEDRON, {0, 2, 1, 7}},
        {G::PRISM, {0, 1, 2, 3, 4, 5}}},
       {0, 1, 1},
       2,
       0},
      {"tet|pyramid|hex",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {1.1, 1.0, 0.05},
        {-0.05, 0.9, 0.0},
        {0.05, 0.0, 1.0},
        {1.0, 0.05, 1.1},
        {1.0, 1.05, 0.95},
        {0.0, 1.0, 1.0},
        {0.45, 0.55, 1.9},
        {0.5, -0.6, 1.6}},
       {{G::TETRAHEDRON, {4, 5, 8, 9}},
        {G::PYRAMID, {4, 5, 6, 7, 8}},
        {G::CUBE, {0, 1, 2, 3, 4, 5, 6, 7}}},
       {0, 0, 1},
       1,
       2},
      {"tet|prism|hex",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {1.1, 1.0, 0.05},
        {-0.05, 0.9, 0.0},
        {0.05, 0.0, 1.0},
        {1.0, 0.05, 1.1},
        {1.0, 1.05, 0.95},
        {0.0, 1.0, 1.0},
        {1.9, 0.1, 0.5},
        {1.95, 1.0, 0.45},
        {1.3, -0.7, 0.5}},
       {{G::TETRAHEDRON, {1, 8, 5, 10}},
        {G::PRISM, {1, 8, 5, 2, 9, 6}},
        {G::CUBE, {0, 1, 2, 3, 4, 5, 6, 7}}},
       {0, 0, 1},
       1,
       2},
      {"prism|pyramid|tet",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {1.1, 1.0, 0.05},
        {-0.05, 0.9, 0.0},
        {0.45, 0.4, 0.95},
        {1.9, 0.15, 0.1},
        {2.0, 1.05, 0.15},
        {1.35, 0.45, 1.05},
        {0.5, -0.8, 0.4}},
       {{G::PRISM, {1, 2, 4, 5, 6, 7}},
        {G::TETRAHEDRON, {0, 1, 4, 8}},
        {G::PYRAMID, {0, 1, 2, 3, 4}}},
       {0, 1, 1},
       2,
       0},
      {"tet|tet",
       {{0.0, 0.0, 0.0},
        {1.0, 0.1, 0.0},
        {0.1, 0.9, 0.05},
        {0.3, 0.35, 0.9},
        {0.2, 0.25, -0.95}},
       {{G::TETRAHEDRON, {0, 1, 2, 3}}, {G::TETRAHEDRON, {0, 2, 1, 4}}},
       {0, 1},
       0,
       1},
  };
}

// The vertex permutations of a reference element that map its edge set onto itself (its
// symmetries, both handednesses; the caller keeps the orientation-preserving ones).
std::vector<std::vector<int>> ElementSymmetries(mfem::Geometry::Type geom)
{
  std::unique_ptr<mfem::Element> el;
  switch (geom)
  {
    case mfem::Geometry::TETRAHEDRON:
      el = std::make_unique<mfem::Tetrahedron>();
      break;
    case mfem::Geometry::PYRAMID:
      el = std::make_unique<mfem::Pyramid>();
      break;
    case mfem::Geometry::PRISM:
      el = std::make_unique<mfem::Wedge>();
      break;
    default:
      el = std::make_unique<mfem::Hexahedron>();
      break;
  }
  const int nv = el->GetNVertices();
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < el->GetNEdges(); i++)
  {
    const int *ev = el->GetEdgeVertices(i);
    edges.emplace_back(std::min(ev[0], ev[1]), std::max(ev[0], ev[1]));
  }
  std::sort(edges.begin(), edges.end());
  std::vector<std::vector<int>> symmetries;
  std::vector<int> perm(nv);
  std::iota(perm.begin(), perm.end(), 0);
  do
  {
    std::vector<std::pair<int, int>> mapped;
    for (const auto &[a, b] : edges)
    {
      mapped.emplace_back(std::min(perm[a], perm[b]), std::max(perm[a], perm[b]));
    }
    std::sort(mapped.begin(), mapped.end());
    if (mapped == edges)
    {
      symmetries.push_back(perm);
    }
  } while (std::next_permutation(perm.begin(), perm.end()));
  return symmetries;
}

// A polynomial of total degree p with fixed generic coefficients (distinct per seed), and
// its gradient.
double Poly(const mfem::Vector &x, int p, int seed, mfem::Vector *grad)
{
  double v = 0.0;
  if (grad)
  {
    grad->SetSize(3);
    *grad = 0.0;
  }
  int n = 0;
  for (int i = 0; i <= p; i++)
  {
    for (int j = 0; i + j <= p; j++)
    {
      for (int k = 0; i + j + k <= p; k++)
      {
        const double c = 0.3 + 0.7 * std::sin(1.3 * (++n) + 0.9 * seed);
        const double xi = std::pow(x(0), i), yj = std::pow(x(1), j), zk = std::pow(x(2), k);
        v += c * xi * yj * zk;
        if (grad)
        {
          (*grad)(0) += (i > 0) ? c * i * std::pow(x(0), i - 1) * yj * zk : 0.0;
          (*grad)(1) += (j > 0) ? c * j * xi * std::pow(x(1), j - 1) * zk : 0.0;
          (*grad)(2) += (k > 0) ? c * k * xi * yj * std::pow(x(2), k - 1) : 0.0;
        }
      }
    }
  }
  return v;
}

// A vector field with components of total degree p and its curl.
void PolyField(const mfem::Vector &x, int p, mfem::Vector &F, mfem::Vector *curl)
{
  mfem::Vector g[3];
  F.SetSize(3);
  for (int c = 0; c < 3; c++)
  {
    F(c) = Poly(x, p, c + 1, &g[c]);
  }
  if (curl)
  {
    curl->SetSize(3);
    (*curl)(0) = g[2](1) - g[1](2);
    (*curl)(1) = g[0](2) - g[2](0);
    (*curl)(2) = g[1](0) - g[0](1);
  }
}

// The mesh of a case with the given vertex orderings of its two varied elements, or nullptr
// if an ordering inverts its element.
std::unique_ptr<mfem::Mesh> MakeCaseMesh(const MixedCase &mixed,
                                         const std::vector<int> &perm0,
                                         const std::vector<int> &perm1)
{
  const int ne = static_cast<int>(mixed.elements.size());
  auto smesh =
      std::make_unique<mfem::Mesh>(3, static_cast<int>(mixed.coords.size()), ne, 0, 3);
  for (const auto &x : mixed.coords)
  {
    smesh->AddVertex(x[0], x[1], x[2]);
  }
  for (int i = 0; i < ne; i++)
  {
    const auto &el = mixed.elements[i];
    std::vector<int> v(el.vertices);
    const std::vector<int> *perm =
        (i == mixed.vary0) ? &perm0 : ((i == mixed.vary1) ? &perm1 : nullptr);
    for (std::size_t j = 0; perm && j < v.size(); j++)
    {
      v[j] = el.vertices[(*perm)[j]];
    }
    mfem::Element *e = smesh->NewElement(el.geom);
    e->SetVertices(v.data());
    e->SetAttribute(i + 1);
    smesh->AddElement(e);
  }
  smesh->FinalizeTopology();
  smesh->Finalize(/*refine*/ true, /*fix_orientation*/ false);
  for (int e = 0; e < ne; e++)
  {
    auto *T = smesh->GetElementTransformation(e);
    T->SetIntPoint(&mfem::Geometries.GetCenter(smesh->GetElementBaseGeometry(e)));
    if (T->Jacobian().Det() <= 0.0)
    {
      return nullptr;
    }
  }
  return smesh;
}

// A smooth field that no polynomial space represents exactly.
void SmoothField(const mfem::Vector &x, mfem::Vector &F)
{
  F.SetSize(3);
  F(0) = std::cos(0.8 * x(1) + 0.3 * x(2));
  F(1) = std::exp(0.4 * x(0)) * std::sin(1.2 * x(2));
  F(2) = std::sin(0.6 * x(0) - 0.9 * x(1)) + x(0) * x(1);
}

// The view of an H(curl) field that every libCEED operator has: the element DOFs of
// Palace's libCEED element restrictions. Returns the largest difference from MFEM's element
// DOFs (GetElementDofs, then the inverse of its DofTransformation) and the largest
// tangential jump across the interior faces of the field evaluated element by element from
// those DOFs, both relative to the largest value.
std::array<double, 2> RestrictionErrors(const FiniteElementSpace &fespace,
                                        const mfem::ParGridFunction &u, int q_order)
{
  auto &pmesh = const_cast<mfem::ParMesh &>(fespace.GetParMesh());
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  std::vector<mfem::Vector> ceed_dofs(pmesh.GetNE());
  for (const auto &[geom, data] : fespace.GetMesh().GetCeedGeomFactorData(ceed))
  {
    if (mfem::Geometry::Dimension[geom] != pmesh.Dimension())
    {
      continue;
    }
    CeedElemRestriction restr = fespace.GetCeedElemRestriction(ceed, geom, data);
    CeedVector l_vec, e_vec;
    PalaceCeedCall(ceed, CeedElemRestrictionCreateVector(restr, &l_vec, &e_vec));
    PalaceCeedCall(ceed, CeedVectorSetArray(l_vec, CEED_MEM_HOST, CEED_COPY_VALUES,
                                            const_cast<double *>(u.HostRead())));
    PalaceCeedCall(ceed, CeedElemRestrictionApply(restr, CEED_NOTRANSPOSE, l_vec, e_vec,
                                                  CEED_REQUEST_IMMEDIATE));
    const CeedScalar *e_data;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(e_vec, CEED_MEM_HOST, &e_data));
    const int P = fespace.Get().GetFE(data.indices[0])->GetDof();
    for (std::size_t i = 0; i < data.indices.size(); i++)
    {
      auto &v = ceed_dofs[data.indices[i]];
      v.SetSize(P);
      std::copy(e_data + i * P, e_data + (i + 1) * P, v.GetData());
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(e_vec, &e_data));
    PalaceCeedCall(ceed, CeedVectorDestroy(&l_vec));
    PalaceCeedCall(ceed, CeedVectorDestroy(&e_vec));
  }
  double err[3] = {0.0, 0.0, 0.0};  // DOF difference, tangential jump, reference
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    fespace.Get().GetElementDofs(e, dofs, dof_trans);
    mfem::Vector mfem_dofs;
    u.GetSubVector(dofs, mfem_dofs);
    dof_trans.InvTransformPrimal(mfem_dofs);
    err[2] = std::max(err[2], mfem_dofs.Normlinf());
    mfem_dofs -= ceed_dofs[e];
    err[0] = std::max(err[0], mfem_dofs.Normlinf());
  }
  double ref = 0.0;
  for (int f = 0; f < pmesh.GetNumFaces(); f++)
  {
    auto *FET = pmesh.GetInteriorFaceTransformations(f);
    if (!FET)
    {
      continue;
    }
    const auto &ir = mfem::IntRules.Get(FET->GetGeometryType(), q_order);
    mfem::Vector v[2], n(3);
    mfem::DenseMatrix shape;
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      FET->SetAllIntPoints(&ir.IntPoint(q));
      mfem::CalcOrtho(FET->Jacobian(), n);
      n /= n.Norml2();
      for (int k = 0; k < 2; k++)
      {
        auto &T = (k == 0) ? *FET->Elem1 : *FET->Elem2;
        const auto &fe = *fespace.Get().GetFE(T.ElementNo);
        shape.SetSize(fe.GetDof(), 3);
        fe.CalcVShape(T, shape);
        v[k].SetSize(3);
        shape.MultTranspose(ceed_dofs[T.ElementNo], v[k]);
        ref = std::max(ref, v[k].Normlinf());
      }
      v[0] -= v[1];
      v[0].Add(-(v[0] * n), n);
      err[1] = std::max(err[1], v[0].Normlinf());
    }
  }
  double out[4] = {err[0], err[2], err[1], ref};
  Mpi::GlobalMax(4, out, fespace.GetComm());
  return {out[0] / out[1], out[2] / std::max(out[3], 1.0e-300)};
}

struct MixedErrors
{
  double gradient = 0.0, curl = 0.0, prolongation = 0.0, element_dofs = 0.0,
         element_jump = 0.0, nd_mass = 0.0, nd_curlcurl = 0.0, rt_mass = 0.0;

  void Max(const MixedErrors &o)
  {
    gradient = std::max(gradient, o.gradient);
    curl = std::max(curl, o.curl);
    prolongation = std::max(prolongation, o.prolongation);
    element_dofs = std::max(element_dofs, o.element_dofs);
    element_jump = std::max(element_jump, o.element_jump);
    nd_mass = std::max(nd_mass, o.nd_mass);
    nd_curlcurl = std::max(nd_curlcurl, o.nd_curlcurl);
    rt_mass = std::max(rt_mass, o.rt_mass);
  }
};

// Largest relative difference of Palace's operator from MFEM's element assembly on the
// local DOFs, for a random local vector.
template <typename T>
double OperatorDifference(const FiniteElementSpace &fespace,
                          mfem::BilinearFormIntegrator *ref)
{
  BilinearForm a_test(fespace);
  a_test.AddDomainIntegrator<T>();
  mfem::BilinearForm a_ref(const_cast<mfem::ParFiniteElementSpace *>(&fespace.Get()));
  a_ref.AddDomainIntegrator(ref);
  a_ref.Assemble();
  a_ref.Finalize();
  Vector x(fespace.GetVSize()), y_test(fespace.GetVSize()), y_ref(fespace.GetVSize());
  x.Randomize(1);
  a_test.PartialAssemble()->Mult(x, y_test);
  a_ref.SpMat().Mult(x, y_ref);
  y_test -= y_ref;
  double diff[2] = {y_test.Normlinf(), y_ref.Normlinf()};
  Mpi::GlobalMax(2, diff, fespace.GetComm());
  return diff[0] / diff[1];
}

// The errors of one mesh at order p: pointwise errors relative to the largest exact value
// of the gradient of a degree-p potential, the curl of a degree-(p - 1) field and the
// prolongation from order p - 1 of a degree-(p - 2) field; the libCEED view of a projected
// smooth field (RestrictionErrors); operator differences.
// The parallel mesh of a case (on one rank, or on the case's two ranks), with the default
// quadrature order of the libCEED operators (set by IoData in a simulation) set first: of
// MFEM's default integrators for the mass (jac = true, extra = 0) or for the curl-curl
// (jac = false, extra = -2) integrator, so that the operators can be compared with MFEM's.
std::unique_ptr<Mesh> MakeCaseParMesh(MPI_Comm comm, const MixedCase &mixed,
                                      mfem::Mesh &smesh, int p, bool q_order_jac,
                                      int q_order_extra_pk)
{
  fem::DefaultIntegrationOrder::p_trial = p;
  fem::DefaultIntegrationOrder::q_order_jac = q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = q_order_extra_pk;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  std::vector<int> partitioning(mixed.ranks);
  if (Mpi::Size(comm) == 1)
  {
    std::fill(partitioning.begin(), partitioning.end(), 0);
  }
  return std::make_unique<Mesh>(
      std::make_unique<mfem::ParMesh>(comm, smesh, partitioning.data()));
}

MixedErrors CaseErrors(MPI_Comm comm, const MixedCase &mixed, mfem::Mesh &smesh, int p)
{
  auto mesh_ptr = MakeCaseParMesh(comm, mixed, smesh, p, true, 0);
  auto &mesh = *mesh_ptr;
  mfem::H1_FECollection h1_fec(p, 3);
  mfem::ND_FECollection nd_fec(p, 3), nd_coarse_fec(std::max(p - 1, 1), 3);
  mfem::RT_FECollection rt_fec(p - 1, 3);
  FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec),
      nd_coarse_fespace(mesh, &nd_coarse_fec), rt_fespace(mesh, &rt_fec);
  auto apply = [](const Operator &op, const FiniteElementSpace &trial,
                  const FiniteElementSpace &test, mfem::Coefficient *s_coeff,
                  mfem::VectorCoefficient *v_coeff, mfem::ParGridFunction &out)
  {
    mfem::ParGridFunction in(const_cast<mfem::ParFiniteElementSpace *>(&trial.Get()));
    s_coeff ? in.ProjectCoefficient(*s_coeff) : in.ProjectCoefficient(*v_coeff);
    Vector in_t(trial.GetTrueVSize()), out_t(test.GetTrueVSize());
    in.GetTrueDofs(in_t);
    op.Mult(in_t, out_t);
    out.SetFromTrueDofs(out_t);
  };
  const auto &grad = nd_fespace.GetDiscreteInterpolator(h1_fespace);
  const auto &curl = rt_fespace.GetDiscreteInterpolator(nd_fespace);
  DiscreteLinearOperator prolongation(nd_coarse_fespace, nd_fespace);
  prolongation.AddDomainInterpolator<IdentityInterpolator>();
  ParOperator P(prolongation.PartialAssemble(), nd_coarse_fespace, nd_fespace, true);
  mfem::FunctionCoefficient potential([p](const mfem::Vector &x)
                                      { return Poly(x, p, 0, nullptr); });
  mfem::VectorFunctionCoefficient field(3, [p](const mfem::Vector &x, mfem::Vector &F)
                                        { PolyField(x, p - 1, F, nullptr); }),
      coarse_field(3, [p](const mfem::Vector &x, mfem::Vector &F)
                   { PolyField(x, std::max(p - 2, 0), F, nullptr); }),
      smooth_field(3, SmoothField);
  mfem::ParGridFunction E(&nd_fespace.Get()), B(&rt_fespace.Get()), U(&nd_fespace.Get()),
      U_smooth(&nd_fespace.Get());
  apply(grad, h1_fespace, nd_fespace, &potential, nullptr, E);
  apply(curl, nd_fespace, rt_fespace, nullptr, &field, B);
  apply(P, nd_coarse_fespace, nd_fespace, nullptr, &coarse_field, U);
  U_smooth.ProjectCoefficient(smooth_field);

  // Pointwise errors [0, 3) and reference magnitudes [3, 6): gradient, curl, prolongation.
  double err[6] = {0.0};
  auto &pmesh = mesh.Get();
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const auto &ir = mfem::IntRules.Get(pmesh.GetElementGeometry(e), 2 * p + 2);
    auto *T = pmesh.GetElementTransformation(e);
    mfem::Vector x, g, F, c, val;
    auto add = [&](int i, const mfem::ParGridFunction &u, const mfem::Vector &exact)
    {
      u.GetVectorValue(*T, T->GetIntPoint(), val);
      val -= exact;
      err[i] = std::max(err[i], val.Normlinf());
      err[i + 3] = std::max(err[i + 3], exact.Normlinf());
    };
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      T->SetIntPoint(&ir.IntPoint(q));
      T->Transform(ir.IntPoint(q), x);
      Poly(x, p, 0, &g);
      add(0, E, g);
      PolyField(x, p - 1, F, &c);
      add(1, B, c);
      err[4] = std::max(err[4], F.Normlinf());
      PolyField(x, std::max(p - 2, 0), F, nullptr);
      add(2, U, F);
    }
  }
  Mpi::GlobalMax(6, err, comm);

  MixedErrors errors;
  errors.gradient = err[0] / err[3];
  errors.curl = err[1] / err[4];
  errors.prolongation = err[2] / err[5];
  const auto restriction = RestrictionErrors(nd_fespace, U_smooth, 2 * p + 2);
  errors.element_dofs = restriction[0];
  errors.element_jump = restriction[1];
  errors.nd_mass = OperatorDifference<VectorFEMassIntegrator>(
      nd_fespace, new mfem::VectorFEMassIntegrator);
  errors.rt_mass = OperatorDifference<VectorFEMassIntegrator>(
      rt_fespace, new mfem::VectorFEMassIntegrator);
  auto cc_mesh = MakeCaseParMesh(comm, mixed, smesh, p, false, -2);
  FiniteElementSpace cc_fespace(*cc_mesh, &nd_fec);
  errors.nd_curlcurl =
      OperatorDifference<CurlCurlIntegrator>(cc_fespace, new mfem::CurlCurlIntegrator);
  return errors;
}

}  // namespace

TEST_CASE("Mixed-element faces: H(curl) and H(div) operators",
          "[MixedElementFaces][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  if (Mpi::Size(comm) > 2)
  {
    SKIP("The mixed-element cases run on one or two ranks");
  }
  const auto cases = MixedCases();
  const auto &mixed = cases[GENERATE(0, 1, 2, 3, 4, 5)];
  const int p = GENERATE(1, 2, 3, 4);
  MixedErrors worst;
  int orderings = 0, k = 0;
  const int stride = (p < 4) ? 1 : 3;  // Every third ordering at order 4, for time.
  for (const auto &perm0 : ElementSymmetries(mixed.elements[mixed.vary0].geom))
  {
    for (const auto &perm1 : ElementSymmetries(mixed.elements[mixed.vary1].geom))
    {
      auto smesh = (k++ % stride == 0) ? MakeCaseMesh(mixed, perm0, perm1) : nullptr;
      if (smesh)
      {
        worst.Max(CaseErrors(comm, mixed, *smesh, p));
        orderings++;
      }
    }
  }
  INFO(mixed.name << " order " << p << ", " << orderings << " orderings: gradient "
                  << worst.gradient << ", curl " << worst.curl << ", prolongation "
                  << worst.prolongation << ", element DOFs " << worst.element_dofs
                  << ", element tangential jump " << worst.element_jump << ", ND mass "
                  << worst.nd_mass << ", ND curl-curl " << worst.nd_curlcurl << ", RT mass "
                  << worst.rt_mass);
  CHECK(orderings > 0);
  CHECK(worst.gradient < 1.0e-10);
  CHECK(worst.curl < 1.0e-10);
  CHECK(worst.prolongation < 1.0e-10);
  CHECK(worst.element_dofs < 1.0e-12);
  CHECK(worst.element_jump < 1.0e-10);
  CHECK(worst.nd_mass < 1.0e-12);
  CHECK(worst.nd_curlcurl < 1.0e-12);
  CHECK(worst.rt_mass < 1.0e-12);
}

namespace
{

// A space that withholds the DofTransformation of element 0 from GetElementDofs, as MFEM's
// ND_FECollection withholds that of pyramids: a restriction whose first element is element
// 0 is then built by the sign-only branch, which the face transformations of the other
// elements must not reach.
class FirstElementDofTransformationWithheld : public mfem::FiniteElementSpace
{
public:
  using mfem::FiniteElementSpace::FiniteElementSpace;
  using mfem::FiniteElementSpace::GetElementDofs;
  void GetElementDofs(int elem, mfem::Array<int> &dofs,
                      mfem::DofTransformation &doftrans) const override
  {
    mfem::FiniteElementSpace::GetElementDofs(elem, dofs, doftrans);
    if (elem == 0)
    {
      doftrans.SetDofTransformation(
          static_cast<const mfem::StatelessDofTransformation *>(nullptr));
    }
  }
};

void BuildElementRestriction(const mfem::FiniteElementSpace &fespace,
                             const std::vector<int> &indices)
{
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  CeedElemRestriction restr = nullptr;
  ceed::InitRestriction(fespace, indices, false, false, false, ceed, &restr);
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&restr));
}

}  // namespace

TEST_CASE("Mixed-element faces: no face transformation reaches a sign-only restriction",
          "[MixedElementFaces][Serial]")
{
#if defined(MFEM_USE_OPENMP)
  SKIP("An exception thrown inside an OpenMP parallel region ends the program");
#endif
  // Tetrahedra at order 2: every element has a non-identity face transformation.
  auto mesh = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON);
  REQUIRE(mesh.GetNE() > 1);
  mfem::ND_FECollection fec(2, mesh.Dimension());
  mfem::FiniteElementSpace fespace(&mesh, &fec);
  REQUIRE(ceed::HasElementDofTransformation(fespace, 0));
  REQUIRE(ceed::HasElementDofTransformation(fespace, 1));
  CHECK_NOTHROW(BuildElementRestriction(fespace, {0, 1}));

  FirstElementDofTransformationWithheld withheld(&mesh, &fec);
  REQUIRE_FALSE(ceed::HasElementDofTransformation(withheld, 0));
  REQUIRE(ceed::HasElementDofTransformation(withheld, 1));
  CHECK_NOTHROW(BuildElementRestriction(withheld, {0}));
  CHECK_THROWS_WITH(BuildElementRestriction(withheld, {0, 1}),
                    Catch::Matchers::ContainsSubstring(
                        "Unexpected DofTransformation for oriented element restriction"));
}

}  // namespace palace
