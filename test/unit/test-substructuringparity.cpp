// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "fixtures.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/substructuringsolver.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

// Parity of the substructured operators with the operators of the regular electrostatic and
// magnetostatic solves: the region and environment operators of SubstructuringSolver must
// add up to the stiffness matrix of LaplaceOperator or CurlCurlOperator (with the sheet
// terms of the superconductor sheets and the two-sided films) and, in magnetostatics, the
// regularization eps M, assembled by Palace on the same mesh and space. The bilinear forms
// are compared on random vectors which vanish on the Dirichlet DOFs, so that the
// elimination of the regular operators does not enter.

using json = nlohmann::json;

namespace palace
{

namespace
{

// The mass regularization of the magnetostatic substructured operator (see
// SubstructuringSolver::Impl::kMagRegularization and kMagRegularizationSheets).
constexpr double kMagRegularization = 1.0e-3, kMagRegularizationSheets = 1.0e-7;

// SubstructuringSolver's scale of the magnetostatic regularization eps (1/mu) M: 1/L^2, L
// the largest side of the bounding box of the mesh.
double RegularizationScale(const mfem::ParMesh &mesh)
{
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  bbmax -= bbmin;
  return 1.0 / (bbmax.Max() * bbmax.Max());
}

// Domain attributes of the cube: the region is x < 1/2 (attribute 1, and 3 in the block y,
// z < 1/2), the environment x > 1/2 (attribute 2, and 4 in the block).
const std::vector<int> kRegion = {1, 3}, kEnvironment = {2, 4};

// Boundary attributes: 1 at x = 0, 2 at x = 1, 3 on the other faces of the cube, and the
// interior patches 5 on y = 1/4 and 6 on y = 3/4 (both 1/4 < x, z < 3/4, so they cross the
// interface x = 1/2), and 7 on the interface itself (3/8 < y < 5/8, 1/4 < z < 3/4).
constexpr int kSheet = 5, kFilm = 6, kInterfaceSheet = 7;

struct ParityMesh
{
  bool tet = false;   // split the hexahedra into tetrahedra
  int curvature = 0;  // order of a curved mesh (0: straight)
  bool port = false;  // the faces y = 0 and y = 1 with attributes 8 and 9, not 3
  bool nc = false;    // refine next to Gamma, nonconformingly, on both sides
};

// The unit cube with n^3 hexahedra (n a multiple of 8), see kRegion and kSheet.
mfem::Mesh MakeParitySerial(int n, const ParityMesh &opt)
{
  MFEM_VERIFY(n % 8 == 0, "The parity mesh needs n a multiple of 8!");
  auto vid = [n](int i, int j, int k) { return (i * (n + 1) + j) * (n + 1) + k; };
  auto in_patch = [n](int a) { return 4 * a >= n && 4 * a < 3 * n; };
  mfem::Mesh mesh(3, (n + 1) * (n + 1) * (n + 1), n * n * n, 6 * n * n + 5 * n * n / 8, 3);
  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      for (int k = 0; k <= n; k++)
      {
        mesh.AddVertex(static_cast<double>(i) / n, static_cast<double>(j) / n,
                       static_cast<double>(k) / n);
      }
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int k = 0; k < n; k++)
      {
        const int v[8] = {vid(i, j, k),
                          vid(i + 1, j, k),
                          vid(i + 1, j + 1, k),
                          vid(i, j + 1, k),
                          vid(i, j, k + 1),
                          vid(i + 1, j, k + 1),
                          vid(i + 1, j + 1, k + 1),
                          vid(i, j + 1, k + 1)};
        const bool block = (2 * j < n) && (2 * k < n);
        mesh.AddHex(v, (2 * i < n) ? (block ? 3 : 1) : (block ? 4 : 2));
      }
    }
  }
  // The faces of the cube, and of the interior patches (quadrilaterals in the planes i, j
  // or k = const, with the first index of their cell a, b in the two other directions).
  auto quad = [&](int dir, int c, int a, int b, int attr)
  {
    std::array<int, 4> q;
    if (dir == 0)
    {
      q = {vid(c, a, b), vid(c, a + 1, b), vid(c, a + 1, b + 1), vid(c, a, b + 1)};
    }
    else if (dir == 1)
    {
      q = {vid(a, c, b), vid(a + 1, c, b), vid(a + 1, c, b + 1), vid(a, c, b + 1)};
    }
    else
    {
      q = {vid(a, b, c), vid(a + 1, b, c), vid(a + 1, b + 1, c), vid(a, b + 1, c)};
    }
    mesh.AddBdrQuad(q.data(), attr);
  };
  for (int a = 0; a < n; a++)
  {
    for (int b = 0; b < n; b++)
    {
      quad(0, 0, a, b, 1);
      quad(0, n, a, b, 2);
      for (int dir = 1; dir < 3; dir++)
      {
        quad(dir, 0, a, b, (opt.port && dir == 1) ? 8 : 3);
        quad(dir, n, a, b, (opt.port && dir == 1) ? 9 : 3);
      }
      if (in_patch(a) && in_patch(b))
      {
        quad(1, n / 4, a, b, kSheet);
        quad(1, 3 * n / 4, a, b, kFilm);
        if (8 * a >= 3 * n && 8 * a < 5 * n)
        {
          quad(0, n / 2, a, b, kInterfaceSheet);
        }
      }
    }
  }
  mesh.FinalizeHexMesh(1, 0, true);
  if (opt.nc)
  {
    mesh.EnsureNCMesh(true);
    mfem::Array<int> marked;
    mfem::Vector c;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementCenter(e, c);
      if ((c(0) > 0.37 && c(0) < 0.5 && c(1) < 0.5) ||
          (c(0) > 0.5 && c(0) < 0.63 && c(1) > 0.5))
      {
        marked.Append(e);
      }
    }
    mesh.GeneralRefinement(marked, 1);
  }
  if (opt.tet)
  {
    mesh = mfem::Mesh::MakeSimplicial(mesh);
  }
  mesh.SetAttributes();
  if (opt.curvature > 0)
  {
    // A smooth deformation fixing the faces of the cube: the interface and the patches are
    // curved.
    mesh.SetCurvature(opt.curvature, false, 3, mfem::Ordering::byVDIM);
    mesh.Transform(
        [](const mfem::Vector &x, mfem::Vector &y)
        {
          const double b =
              std::sin(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::sin(M_PI * x(2));
          y = x;
          y(0) += 0.06 * b * std::cos(M_PI * x(1));
          y(1) += 0.05 * b * std::cos(M_PI * x(2));
          y(2) += 0.04 * b * std::cos(M_PI * x(0));
        });
  }
  return mesh;
}

// An orthonormal frame for anisotropic materials.
json MaterialAxes(double a, double b)
{
  return {{std::cos(a), std::sin(a), 0.0},
          {-std::sin(a) * std::cos(b), std::cos(a) * std::cos(b), std::sin(b)},
          {std::sin(a) * std::sin(b), -std::cos(a) * std::sin(b), std::cos(b)}};
}

// Anisotropic materials, with a rotated frame in each side.
json Materials(const std::string &property)
{
  return {{{"Attributes", {1}},
           {property, {2.0, 3.0, 5.0}},
           {"MaterialAxes", MaterialAxes(0.3, 0.7)}},
          {{"Attributes", {3}}, {property, 4.0}},
          {{"Attributes", {2}},
           {property, {7.0, 1.5, 2.0}},
           {"MaterialAxes", MaterialAxes(-0.5, 1.1)}},
          {{"Attributes", {4}}, {property, {1.0, 6.0, 2.5}}}};
}

json Substructuring()
{
  return {{"Region", {{"Attributes", kRegion}}},
          {"Environment", {{"Attributes", kEnvironment}}},
          {"SaveModel", false}};
}

json ElectrostaticConfig(int order)
{
  return {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
          {"Model", {{"Mesh", "test.msh"}}},
          {"Domains", {{"Materials", Materials("Permittivity")}}},
          {"Boundaries",
           {{"Terminal",
             {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
          {"Solver", {{"Order", order}, {"Substructuring", Substructuring()}}}};
}

// Superconductor sheets: one-sided sheets crossing the interface and on it, and a
// two-sided film crossing the interface; London domains: a domain on each side.
json MagnetostaticConfig(int order, bool sheets, bool film, bool london = false)
{
  json sc = json::array();
  if (sheets)
  {
    sc.push_back(
        {{"Attributes", {kSheet}}, {"PenetrationDepth", 0.2}, {"Thickness", 0.05}});
    sc.push_back({{"Attributes", {kInterfaceSheet}}, {"KineticInductance", 0.5}});
  }
  if (film)
  {
    sc.push_back({{"Attributes", {kFilm}},
                  {"PenetrationDepth", 0.1},
                  {"Thickness", 0.2},
                  {"TwoSided", true}});
  }
  json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", Materials("Permeability")}}},
                 {"Boundaries", {{"PEC", {{"Attributes", {1, 2}}}}}},
                 {"Solver", {{"Order", order}, {"Substructuring", Substructuring()}}}};
  if (!sc.empty())
  {
    config["Boundaries"]["Superconductor"] = sc;
  }
  if (london)
  {
    for (auto &material : config["Domains"]["Materials"])
    {
      const int a = material["Attributes"][0].get<int>();
      if (a == kRegion.back() || a == kEnvironment.back())
      {
        material["LondonDepth"] = 0.5;
      }
    }
  }
  return config;
}

// A random vector on the true DOFs, zero on the given ones.
Vector RandomVector(int n, const mfem::Array<int> &zero, unsigned seed)
{
  std::mt19937 gen(seed + 7919u * static_cast<unsigned>(Mpi::Rank(Mpi::World())));
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Vector x(n);
  double *xx = x.HostWrite();
  for (int i = 0; i < n; i++)
  {
    xx[i] = dist(gen);
  }
  const int *z = zero.HostRead();
  for (int k = 0; k < zero.Size(); k++)
  {
    xx[z[k]] = 0.0;
  }
  x.UseDevice(true);
  return x;
}

double Form(const Operator &A, const Vector &x, const Vector &y)
{
  Vector t(y.Size());
  t.UseDevice(true);
  A.Mult(y, t);
  return linalg::Dot(Mpi::World(), x, t);
}

// The largest difference of the forms a and b over the pairs of the given vectors,
// relative to the scale sqrt(b(x, x) b(y, y)) of each pair (b is positive definite on the
// vectors).
template <typename A, typename B>
double FormDifference(const std::vector<Vector> &v, const A &a, const B &b)
{
  std::vector<double> diag(v.size());
  for (std::size_t i = 0; i < v.size(); i++)
  {
    diag[i] = b(v[i], v[i]);
    REQUIRE(diag[i] > 0.0);
  }
  double diff = 0.0;
  for (std::size_t i = 0; i < v.size(); i++)
  {
    for (std::size_t j = i; j < v.size(); j++)
    {
      const double bij = (i == j) ? diag[i] : b(v[i], v[j]);
      diff = std::max(diff, std::abs(a(v[i], v[j]) - bij) / std::sqrt(diag[i] * diag[j]));
    }
  }
  return diff;
}

// The true DOFs of a space on the given boundary attributes.
mfem::Array<int> EssentialTrueDofs(const FiniteElementSpace &fespace,
                                   const std::vector<int> &attrs)
{
  const auto &mesh = fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker = mesh::AttrToMarker(bdr_attr_max, attrs), dofs;
  fespace.GetEssentialTrueDofs(marker, dofs);
  return dofs;
}

// The rounding of a bilinear form of the assembled operators on random vectors.
constexpr double kFormTol = 1.0e-12;

// The magnetostatic form on gradients which vanish on the sheets, where only the
// regularization remains: the rounding of the curl-curl term of the assembled operators,
// relative to the regularization (eps h^2 of it, eps down to 1e-7), bounds the agreement.
constexpr double kRegularizationTol = 1.0e-5;

void CheckElectrostatic(IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  SubstructuringSolver ss(iodata, mesh);
  LaplaceOperator laplace_op(iodata, mesh);
  auto K = laplace_op.GetStiffnessMatrix();
  const auto &fespace = laplace_op.GetH1Space();
  REQUIRE(ss.GlobalTrueVSize() == fespace.GlobalTrueVSize());
  const auto dbc = EssentialTrueDofs(fespace, {1, 2});
  std::vector<Vector> v;
  for (unsigned s = 0; s < 3; s++)
  {
    v.push_back(RandomVector(fespace.GetTrueVSize(), dbc, s));
  }
  auto native = [&](const Vector &x, const Vector &y) { return Form(*K, x, y); };
  const double d_energy = FormDifference(
      v, [&](const Vector &x, const Vector &y) { return ss.MutualEnergy(x, y); }, native);
  const double d_system = FormDifference(
      v, [&](const Vector &x, const Vector &y) { return ss.MutualSystemEnergy(x, y); },
      native);
  CAPTURE(d_energy, d_system);
  CHECK(d_energy < kFormTol);
  CHECK(d_system < kFormTol);
}

void CheckMagnetostatic(IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh,
                        bool with_sheets)
{
  SubstructuringSolver ss(iodata, mesh);
  CurlCurlOperator curlcurl_op(iodata, mesh);
  const auto &nd_fespace = curlcurl_op.GetSolutionNDSpace();
  const auto &h1_fespace = curlcurl_op.GetSolutionH1Spaces().GetFinestFESpace();
  REQUIRE(ss.GlobalTrueVSize() == nd_fespace.GlobalTrueVSize());
  const auto &mat_op = curlcurl_op.GetMaterialOp();

  // The stiffness matrix of the regular solve (curl-curl, London and sheet terms), the
  // volume terms alone (the energy operator), and the regularization eps M off the London
  // domains.
  auto K = curlcurl_op.GetStiffnessMatrix();
  const bool london = mat_op.HasLondonDepth();
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetCurlCurlInvPermeability());
  MaterialPropertyCoefficient london_func(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetInvLondonDepth());
  BilinearForm kcc(nd_fespace);
  kcc.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  if (london)
  {
    kcc.AddDomainIntegrator<VectorFEMassIntegrator>(london_func);
  }
  ParOperator Kcc(kcc.Assemble(false), nd_fespace);
  std::vector<int> domain_attrs(kRegion);
  domain_attrs.insert(domain_attrs.end(), kEnvironment.begin(), kEnvironment.end());
  mfem::DenseTensor muinv(mat_op.GetInvPermeability());
  for (int i = 0; i < muinv.SizeK(); i++)
  {
    if (mat_op.GetInvLondonDepth()(i).MaxMaxNorm() > 0.0)
    {
      muinv(i) = 0.0;
    }
  }
  MaterialPropertyCoefficient eps_func(
      mat_op.GetAttributeToMaterial(), muinv,
      ((with_sheets || london) ? kMagRegularizationSheets : kMagRegularization) *
          RegularizationScale(static_cast<const mfem::ParMesh &>(nd_fespace.GetMesh())));
  eps_func.RestrictCoefficient(mat_op.GetCeedAttributes(domain_attrs));
  BilinearForm m(nd_fespace);
  m.AddDomainIntegrator<VectorFEMassIntegrator>(eps_func);
  ParOperator M(m.Assemble(false), nd_fespace);
  const auto dbc = EssentialTrueDofs(nd_fespace, {1, 2});
  M.SetEssentialTrueDofs(dbc, Operator::DiagonalPolicy::DIAG_ZERO);
  SumOperator A(nd_fespace.GetTrueVSize(), nd_fespace.GetTrueVSize());
  A.AddOperator(*K);
  A.AddOperator(M);

  // Random vectors; gradients of random scalars zero on the PEC boundary, on which the
  // curl-curl term vanishes; and gradients of random scalars also zero on the sheets and
  // films, on which only the regularization remains (lifted to all the versions of the
  // DOFs on a broken space).
  const auto h1_dbc = EssentialTrueDofs(h1_fespace, {1, 2});
  const auto &G = nd_fespace.GetDiscreteInterpolator(h1_fespace);
  const auto h1_sheet_dbc =
      EssentialTrueDofs(curlcurl_op.GetH1Space(), {1, 2, kSheet, kFilm, kInterfaceSheet});
  const auto &G0 =
      curlcurl_op.GetNDSpace().GetDiscreteInterpolator(curlcurl_op.GetH1Space());
  const auto all_versions =
      curlcurl_op.IsBroken()
          ? nd_fespace.BuildVersionMap(FiniteElementSpace::VersionMap::ALL)
          : nullptr;
  std::vector<Vector> v, g, r;
  for (unsigned s = 0; s < 3; s++)
  {
    v.push_back(RandomVector(nd_fespace.GetTrueVSize(), dbc, s));
    const Vector phi = RandomVector(h1_fespace.GetTrueVSize(), h1_dbc, 10 + s);
    Vector &gs = g.emplace_back(nd_fespace.GetTrueVSize());
    gs.UseDevice(true);
    G.Mult(phi, gs);
    const Vector psi =
        RandomVector(curlcurl_op.GetH1Space().GetTrueVSize(), h1_sheet_dbc, 20 + s);
    Vector rs(curlcurl_op.GetNDSpace().GetTrueVSize());
    rs.UseDevice(true);
    G0.Mult(psi, rs);
    if (all_versions)
    {
      Vector &rb = r.emplace_back(nd_fespace.GetTrueVSize());
      rb.UseDevice(true);
      all_versions->MultTranspose(rs, rb);
    }
    else
    {
      r.push_back(std::move(rs));
    }
  }
  auto system = [&](const Vector &x, const Vector &y)
  { return ss.MutualSystemEnergy(x, y); };
  auto native = [&](const Vector &x, const Vector &y) { return Form(A, x, y); };
  const double d_energy = FormDifference(
      v, [&](const Vector &x, const Vector &y) { return ss.MutualEnergy(x, y); },
      [&](const Vector &x, const Vector &y) { return Form(Kcc, x, y); });
  const double d_system = FormDifference(v, system, native);
  const double d_gradient = FormDifference(g, system, native);
  const double d_regularization = FormDifference(r, system, native);
  CAPTURE(d_energy, d_system, d_gradient, d_regularization);
  CHECK(d_energy < kFormTol);
  CHECK(d_system < kFormTol);
  if (with_sheets)
  {
    CHECK(d_gradient < kFormTol);
  }
  CHECK(d_regularization < kRegularizationTol);
}

// With by_side, the processes own slabs of x (process floor(np x^2)), so that the last ones
// own no element of x < 1/2, nor a neighbor of one. With first_only, the first process owns
// every element, and the others none.
std::vector<std::unique_ptr<Mesh>> MakeMesh(mfem::Mesh &&serial, bool by_side = false,
                                            bool first_only = false)
{
  const int np = Mpi::Size(Mpi::World());
  std::vector<int> part;
  if ((by_side || first_only) && np > 1)
  {
    part.assign(serial.GetNE(), 0);
    mfem::Vector c;
    for (int e = 0; by_side && e < serial.GetNE(); e++)
    {
      serial.GetElementCenter(e, c);
      part[e] = std::min(np - 1, static_cast<int>(np * c(0) * c(0)));
    }
  }
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(
      Mpi::World(), serial, part.empty() ? nullptr : part.data())));
  return mesh;
}

// The mixed mesh fichera-mixed-p2 (tetrahedra, hexahedra and prisms, so that the elements
// of a geometry are not numbered contiguously), refined twice, with the domain attributes
// of the cube: 1 and 3 for x < 0, 2 and 4 otherwise, 3 and 4 for y < 0.
mfem::Mesh MakeMixedSerial()
{
  mfem::Mesh mesh(PALACE_TEST_DATA_DIR "/mesh/fichera-mixed-p2.mesh", 1, 1);
  mesh.UniformRefinement();
  mesh.UniformRefinement();
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    mesh.SetAttribute(e, (c(0) < 0.0) ? ((c(1) < 0.0) ? 3 : 1) : ((c(1) < 0.0) ? 4 : 2));
  }
  mesh.SetAttributes();
  return mesh;
}

// A coefficient with the value v(attr) on the libCEED attributes of each domain attribute.
MaterialPropertyCoefficient DomainCoefficient(const Mesh &mesh,
                                              const std::function<double(int)> &v)
{
  const int n = static_cast<int>(mesh.MaxCeedAttribute());
  mfem::Array<int> attr_mat(n);
  mfem::DenseTensor mat_coeff(1, 1, n);
  for (const auto &[attr, ceed_attr] : mesh.GetCeedAttributes())
  {
    attr_mat[ceed_attr - 1] = ceed_attr - 1;
    mat_coeff(0, 0, ceed_attr - 1) = v(attr);
  }
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

// A coefficient with the value v(attr, nbr_attr) on the libCEED boundary attribute of each
// boundary attribute and neighboring domain attribute.
MaterialPropertyCoefficient BoundaryCoefficient(const Mesh &mesh,
                                                const std::function<double(int, int)> &v)
{
  const int n = static_cast<int>(mesh.MaxCeedBdrAttribute());
  mfem::Array<int> attr_mat(n);
  mfem::DenseTensor mat_coeff(1, 1, n);
  for (const auto &[attr, bdr_attr_map] : mesh.GetCeedBdrAttributes())
  {
    for (const auto &[nbr_attr, ceed_attr] : bdr_attr_map)
    {
      attr_mat[ceed_attr - 1] = ceed_attr - 1;
      mat_coeff(0, 0, ceed_attr - 1) = v(attr, nbr_attr);
    }
  }
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

// The number of elements of the libCEED sub-operators of an operator, summed over the
// processes.
long long NumAssembledElements(const ceed::Operator &op)
{
  long long n = 0;
  for (std::size_t i = 0; i < op.Size(); i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    CeedInt num_sub_ops;
    CeedOperator *sub_ops;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[i], &num_sub_ops));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[i], &sub_ops));
    for (CeedInt j = 0; j < num_sub_ops; j++)
    {
      CeedInt num_elem;
      PalaceCeedCall(ceed, CeedOperatorGetNumElements(sub_ops[j], &num_elem));
      n += num_elem;
    }
  }
  Mpi::GlobalSum(1, &n, Mpi::World());
  return n;
}

// The number of elements with one of the given domain attributes, summed over the
// processes.
long long NumElements(const mfem::ParMesh &mesh, const std::vector<int> &attrs)
{
  long long n = 0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    n += (std::find(attrs.begin(), attrs.end(), mesh.GetAttribute(e)) != attrs.end());
  }
  Mpi::GlobalSum(1, &n, Mpi::World());
  return n;
}

// The number of entries of the local matrix A which differ from those of B (in the sparsity
// pattern, or in value by more than tol relative to the largest entry of B, over the
// processes), checking that the matrices are not zero.
long long NumDifferentEntries(hypre::HypreCSRMatrix &A, hypre::HypreCSRMatrix &B,
                              double tol)
{
  hypre_CSRMatrixMigrate(A, HYPRE_MEMORY_HOST);
  hypre_CSRMatrixMigrate(B, HYPRE_MEMORY_HOST);
  REQUIRE(A.Height() == B.Height());
  REQUIRE(A.NNZ() == B.NNZ());
  double scale = 0.0;
  for (int k = 0; k < B.NNZ(); k++)
  {
    scale = std::max(scale, std::abs(B.GetData()[k]));
  }
  Mpi::GlobalMax(1, &scale, Mpi::World());
  REQUIRE(scale > 0.0);
  long long diff = 0;
  for (int i = 0; i <= A.Height(); i++)
  {
    diff += (A.GetI()[i] != B.GetI()[i]);
  }
  for (int k = 0; k < A.NNZ(); k++)
  {
    // Written so that a NaN counts as a difference.
    diff += (A.GetJ()[k] != B.GetJ()[k]) ||
            !(std::abs(A.GetData()[k] - B.GetData()[k]) <= tol * scale);
  }
  Mpi::GlobalSum(1, &diff, Mpi::World());
  return diff;
}

}  // namespace

TEST_CASE("Substructured electrostatic operators match LaplaceOperator",
          "[substructure][Serial][Parallel][GPU]")
{
  const auto [order, opt] = GENERATE(
      table<int, ParityMesh>({{1, {}}, {2, {}}, {3, {}}, {2, {false, 3}}, {2, {true, 2}}}));
  CAPTURE(order, opt.tet, opt.curvature);
  IoData iodata(ElectrostaticConfig(order), false);
  CheckElectrostatic(iodata, MakeMesh(MakeParitySerial(8, opt)));
}

TEST_CASE("Substructured magnetostatic operators match CurlCurlOperator",
          "[substructure][Serial][Parallel][GPU]")
{
  // Without sheets, with one-sided sheets (crossing the interface and on it), with a
  // two-sided film crossing the interface, and with both.
  const auto [order, sheets, film, opt] = GENERATE(
      table<int, bool, bool, ParityMesh>({{1, false, false, {}},
                                          {2, false, false, {}},
                                          {2, false, false, {false, 0, false, true}},
                                          {1, true, false, {}},
                                          {2, true, false, {}},
                                          {2, true, false, {false, 2}},
                                          {1, false, true, {}},
                                          {2, false, true, {}},
                                          {1, true, true, {}},
                                          {1, true, true, {true, 2}}}));
  CAPTURE(order, sheets, film, opt.tet, opt.curvature, opt.nc);
  IoData iodata(MagnetostaticConfig(order, sheets, film), false);
  CheckMagnetostatic(iodata, MakeMesh(MakeParitySerial(8, opt)), sheets || film);
}

TEST_CASE(
    "Substructured magnetostatic operators match CurlCurlOperator with London domains",
    "[substructure][Serial][Parallel][GPU]")
{
  // A London domain on each side, whose London mass is in the energy operator and lifts
  // the regularization, with and without sheets.
  const auto [order, sheets] =
      GENERATE(table<int, bool>({{1, false}, {2, false}, {1, true}}));
  CAPTURE(order, sheets);
  IoData iodata(MagnetostaticConfig(order, sheets, false, true), false);
  CheckMagnetostatic(iodata, MakeMesh(MakeParitySerial(8, {})), sheets);
}

TEST_CASE("Substructured third-order magnetostatic operators match CurlCurlOperator",
          "[substructure][Serial][Parallel][GPU][Long]")
{
  // As above, without sheets and with both kinds.
  const bool sheets = GENERATE(false, true);
  CAPTURE(sheets);
  IoData iodata(MagnetostaticConfig(3, sheets, sheets), false);
  CheckMagnetostatic(iodata, MakeMesh(MakeParitySerial(8, {})), sheets);
}

TEST_CASE("Substructured operators match on processes owning one side only",
          "[substructure][Parallel]")
{
  // A process without an element of a side, or a neighbor of one, assembles nothing of
  // its operators, but takes part in their parallel assembly.
  IoData electrostatic(ElectrostaticConfig(2), false);
  CheckElectrostatic(electrostatic, MakeMesh(MakeParitySerial(8, {}), true));
  IoData magnetostatic(MagnetostaticConfig(1, true, true), false);
  CheckMagnetostatic(magnetostatic, MakeMesh(MakeParitySerial(8, {}), true), true);
}

TEST_CASE("BilinearForm assembles the elements of its assembly domains only",
          "[substructure][Serial][Parallel][GPU]")
{
  // A form restricted to the region or to the environment assembles the elements of these
  // domains where its coefficients are nonzero, and the boundary elements whose neighboring
  // domain (the one with the lower attribute, for the interior boundaries) is one of them,
  // which is the operator of the form with its coefficients restricted to them. The
  // coefficient vanishes on attribute 3, a domain of the region. The meshes: the cube in
  // hexahedra and in tetrahedra, and a mixed mesh.
  const int mesh_type = GENERATE(0, 1, 2);
  const bool by_side = GENERATE(false, true);
  const bool nd = GENERATE(false, true);
  CAPTURE(mesh_type, by_side, nd);
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  auto mesh = MakeMesh(
      (mesh_type < 2) ? MakeParitySerial(8, {mesh_type == 1}) : MakeMixedSerial(), by_side);
  const auto v = [](int attr) { return (attr == 3) ? 0.0 : 1.0 + attr; };
  const auto vb = [](int attr, int nbr_attr) { return 1.0 + 0.1 * attr + 0.01 * nbr_attr; };
  mfem::H1_FECollection h1_fec(order, 3);
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace fespace(*mesh[0],
                             nd ? static_cast<mfem::FiniteElementCollection *>(&nd_fec)
                                : static_cast<mfem::FiniteElementCollection *>(&h1_fec));

  // Restricting the elements keeps the order of the sums of the full assembly on the CPU.
  const double tol = (ceed::Print().find("/gpu/") == std::string::npos)
                         ? 0.0
                         : 8.0 * std::numeric_limits<double>::epsilon();
  for (const auto &attrs : {kRegion, kEnvironment})
  {
    CAPTURE(attrs);
    auto f = DomainCoefficient(*mesh[0], v), fb = BoundaryCoefficient(*mesh[0], vb);
    auto f_ref = f, fb_ref = fb;
    f_ref.RestrictCoefficient(mesh[0]->GetCeedAttributes(attrs));
    fb_ref.RestrictCoefficient(mesh[0]->GetCeedBdrAttributesByNeighbor(attrs));
    BilinearForm a(fespace), a_ref(fespace), a_dom(fespace);
    if (nd)
    {
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(f, f);
      a_ref.AddDomainIntegrator<CurlCurlMassIntegrator>(f_ref, f_ref);
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb);
      a_ref.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb_ref);
      a_dom.AddDomainIntegrator<CurlCurlMassIntegrator>(f, f);
    }
    else
    {
      a.AddDomainIntegrator<DiffusionIntegrator>(f);
      a_ref.AddDomainIntegrator<DiffusionIntegrator>(f_ref);
      a.AddBoundaryIntegrator<MassIntegrator>(fb);
      a_ref.AddBoundaryIntegrator<MassIntegrator>(fb_ref);
      a_dom.AddDomainIntegrator<DiffusionIntegrator>(f);
    }
    a.SetAssemblyDomains(attrs);
    a_dom.SetAssemblyDomains(attrs);
    constexpr bool skip_zeros = true;
    auto A = a.FullAssemble(skip_zeros), A_ref = a_ref.FullAssemble(skip_zeros);
    CHECK(NumDifferentEntries(*A, *A_ref, tol) == 0);

    // The same operator applied matrix-free (the sums of shared DOFs may be reordered).
    auto op = a.PartialAssemble(), op_ref = a_ref.PartialAssemble();
    Vector x(op_ref->Width()), y(op_ref->Height()), y_ref(op_ref->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    y_ref.UseDevice(true);
    x.Randomize(1);
    op->Mult(x, y);
    op_ref->Mult(x, y_ref);
    const double norm_ref = linalg::Norml2(Mpi::World(), y_ref);
    REQUIRE(norm_ref > 0.0);
    y -= y_ref;
    CHECK(linalg::Norml2(Mpi::World(), y) <= 1.0e-14 * norm_ref);

    // Only the elements of the domains with a nonzero coefficient are assembled.
    std::vector<int> nonzero_attrs;
    std::copy_if(attrs.begin(), attrs.end(), std::back_inserter(nonzero_attrs),
                 [&](int attr) { return v(attr) != 0.0; });
    CHECK(NumAssembledElements(*a_dom.PartialAssemble()) ==
          NumElements(mesh[0]->Get(), nonzero_attrs));
  }
}

TEST_CASE("AssemblyDomains restricts the assembly while a Restriction is in scope",
          "[substructure][Serial]")
{
  AssemblyDomains domains;
  REQUIRE(domains.Get() == nullptr);
  {
    AssemblyDomains::Restriction region(domains, kRegion);
    REQUIRE(domains.Get() != nullptr);
    CHECK(*domains.Get() == kRegion);
    {
      AssemblyDomains::Restriction environment(domains, kEnvironment);
      REQUIRE(domains.Get() != nullptr);
      CHECK(*domains.Get() == kEnvironment);
    }
    REQUIRE(domains.Get() != nullptr);
    CHECK(*domains.Get() == kRegion);
  }
  CHECK(domains.Get() == nullptr);

  // A form restricted after the scopes ended assembles every element.
  fem::DefaultIntegrationOrder::p_trial = 1;
  auto mesh = MakeMesh(MakeParitySerial(8, {}));
  mfem::H1_FECollection fec(1, 3);
  FiniteElementSpace fespace(*mesh[0], &fec);
  auto f = DomainCoefficient(*mesh[0], [](int attr) { return 1.0 + attr; });
  BilinearForm a(fespace);
  a.AddDomainIntegrator<DiffusionIntegrator>(f);
  domains.Restrict(a);
  CHECK(NumAssembledElements(*a.PartialAssemble()) ==
        NumElements(mesh[0]->Get(), {1, 2, 3, 4}));
}

TEST_CASE(
    "Substructured magnetostatic operators match CurlCurlOperator with cracked sheets",
    "[substructure][Serial][Parallel][GPU]")
{
  // The sheets of a mesh cracked along them, of which each face carries half of the sheet
  // term. The mesh preprocessing does not crack superconductor sheets, so the attributes
  // are marked cracked on the uncracked mesh.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  IoData iodata(MagnetostaticConfig(order, true, false), false);
  iodata.boundaries.cracked_attributes = {kSheet, kInterfaceSheet};
  CheckMagnetostatic(iodata, MakeMesh(MakeParitySerial(8, {})), true);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Substructured sheet faces on the interface belong to one side",
                 "[substructure][Serial][Parallel][GPU]")
{
  // A face of a sheet on Gamma belongs to the side of its neighbor with the lower
  // attribute: the region with the region attributes numbered below those of the
  // environment, the environment with them numbered above. Either way, the inductance of a
  // surface current (on y = 0, from the PEC face x = 0 to x = 1, which the PEC face y = 1
  // connects) of the offline run must be that of an online run on the full mesh and of one
  // on the mesh of the region, which keeps the environment neighbors of these faces. Also
  // partitioned by side, and with the first process owning every element: the processes
  // without true DOFs take part in the collectives of the others.
  // Not structured bindings, which the lambdas below cannot capture with OpenMP (Clang 18).
  bool region_below, by_side, first_only;
  int order;
  std::tie(region_below, order, by_side, first_only) =
      GENERATE(table<bool, int, bool, bool>({{true, 1, false, false},
                                             {false, 1, false, false},
                                             {true, 2, false, false},
                                             {false, 2, false, false},
                                             {true, 1, true, false},
                                             {false, 1, true, false},
                                             {true, 1, false, true}}));
  CAPTURE(region_below, order, by_side, first_only);
  const std::string model_path = (temp_dir / "sheet_side" / "environment.model").string();
  auto make_config = [&](const std::string &mode)
  {
    json config = MagnetostaticConfig(order, true, false);
    config["Boundaries"]["PEC"]["Attributes"] = {1, 2, 9};
    config["Boundaries"]["SurfaceCurrent"] = {
        {{"Index", 1}, {"Attributes", {8}}, {"Direction", "+X"}}};
    auto &sub = config["Solver"]["Substructuring"];
    sub["Region"]["Attributes"] = region_below ? kRegion : kEnvironment;
    sub["Environment"]["Attributes"] = region_below ? kEnvironment : kRegion;
    test::SetModelRun(config, mode, model_path);
    return IoData(config, false);
  };
  // The energy of the surface current, and whether the environment was solved.
  auto energy = [&](const IoData &iodata, mfem::Mesh &&serial)
  {
    auto mesh = MakeMesh(std::move(serial), by_side, first_only);
    Vector J;
    CurlCurlOperator(iodata, mesh).GetCurrentExcitationVector(1, J);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    const double E = ss.MagnetostaticEnergyMatrix({1}, {J}, {}, {}, {}, nullptr)(0, 0);
    CHECK(ss.EnvironmentFactored() == !iodata.solver.substructuring->Online());
    return E;
  };
  const ParityMesh opt{false, 0, true};
  const double E_off = energy(make_config("Offline"), MakeParitySerial(8, opt));
  const double E_on = energy(make_config("Online"), MakeParitySerial(8, opt));
  mfem::Mesh region_mesh(RegionMeshPath(model_path).c_str());
  const int region_ne = region_mesh.GetNE();
  const double E_region = energy(make_config("Online"), std::move(region_mesh));
  CAPTURE(E_off, E_on, E_region, region_ne);
  CHECK(region_ne < 8 * 8 * 8);  // the mesh of the region, not the full mesh
  CHECK(E_off > 0.0);
  CHECK(std::abs(E_on - E_off) <= 1.0e-9 * E_off);
  CHECK(std::abs(E_region - E_off) <= 1.0e-9 * E_off);
  // Refined in the region next to the sheet on Gamma, nonconformingly as by the adaptation
  // of an online run, its faces split and keep their side.
  mfem::Mesh refined = MakeParitySerial(8, opt);
  refined.EnsureNCMesh(true);
  mfem::Array<int> marked;
  mfem::Vector c;
  for (int e = 0; e < refined.GetNE(); e++)
  {
    refined.GetElementCenter(e, c);
    if ((region_below ? (c(0) > 0.37 && c(0) < 0.5) : (c(0) > 0.5 && c(0) < 0.63)) &&
        c(1) > 0.37 && c(1) < 0.63)
    {
      marked.Append(e);
    }
  }
  refined.GeneralRefinement(marked, 1, 0);
  const double E_refined = energy(make_config("Online"), std::move(refined));
  CAPTURE(E_refined);
  CHECK(std::abs(E_refined - E_off) <= 1.0e-2 * E_off);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Substructuring tells apart the copies of a crack across the interface",
                 "[substructure][Serial][Parallel]")
{
  // A cracked boundary across the interface (the patch y = 1/4, of zero charge) doubles its
  // interface DOFs, which a saved model matches by their positions: displaced as by
  // default, the copies are told apart and an online run gives the capacitance of the
  // offline one, and not displaced, the run is rejected.
  const std::string mesh_path = (temp_dir / "cracked.mesh").string();
  if (Mpi::Root(Mpi::World()))
  {
    MakeParitySerial(8, {}).Save(mesh_path);
  }
  Mpi::Barrier(Mpi::World());
  const std::string model_path = (temp_dir / "cracked" / "environment.model").string();
  auto capacitance = [&](const std::string &mode, double displacement)
  {
    json config = ElectrostaticConfig(1);
    config["Model"] = {{"Mesh", mesh_path}, {"CrackDisplacementFactor", displacement}};
    config["Boundaries"]["ZeroCharge"] = {{"Attributes", {kSheet}}};
    test::SetModelRun(config, mode, model_path);
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  CHECK_THROWS_WITH(
      capacitance("Offline", 0.0),
      Catch::Matchers::ContainsSubstring("crosses the substructuring interface"));
  const mfem::DenseMatrix C_off = capacitance("Offline", 1.0e-12);
  mfem::DenseMatrix C_on = capacitance("Online", 1.0e-12);
  C_on -= C_off;
  CHECK(C_on.MaxMaxNorm() <= 1.0e-10 * C_off.MaxMaxNorm());
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Substructured sheet faces on the interface keep their sides online",
                 "[substructure][Serial][Parallel]")
{
  // Renumbering the region above the environment on the mesh of the region moves the sheet
  // faces on Gamma, which belong to the side of their neighbor with the lower attribute, to
  // the environment, which the online run does not assemble: the run is rejected.
  const std::string model_path = (temp_dir / "sheet_moved" / "environment.model").string();
  constexpr int shift = 10;
  auto make_config = [&](const std::string &mode, bool renumbered)
  {
    json config = MagnetostaticConfig(1, true, false);
    test::SetModelRun(config, mode, model_path);
    auto &sub = config["Solver"]["Substructuring"];
    if (renumbered)
    {
      std::vector<int> region;
      for (int a : kRegion)
      {
        region.push_back(a + shift);
      }
      sub["Region"]["Attributes"] = region;
      for (auto &m : config["Domains"]["Materials"])
      {
        const int a = m["Attributes"][0].get<int>();
        if (std::find(kRegion.begin(), kRegion.end(), a) != kRegion.end())
        {
          m["Attributes"][0] = a + shift;
        }
      }
    }
    return IoData(config, false);
  };
  {
    IoData iodata = make_config("Offline", false);
    auto mesh = MakeMesh(MakeParitySerial(8, {}));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  }
  auto online = [&](bool renumbered)
  {
    mfem::Mesh region_mesh(RegionMeshPath(model_path).c_str());
    if (renumbered)
    {
      for (int e = 0; e < region_mesh.GetNE(); e++)
      {
        const int a = region_mesh.GetAttribute(e);
        if (std::find(kRegion.begin(), kRegion.end(), a) != kRegion.end())
        {
          region_mesh.SetAttribute(e, a + shift);
        }
      }
      region_mesh.SetAttributes();
    }
    IoData iodata = make_config("Online", renumbered);
    auto mesh = MakeMesh(std::move(region_mesh));
    SubstructuringSolver ss(iodata, mesh);
  };
  CHECK_NOTHROW(online(false));
#if defined(MFEM_USE_EXCEPTIONS)
  CHECK_THROWS_WITH(online(true),
                    Catch::Matchers::ContainsSubstring("not on the same sides"));
#endif
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Substructuring model with a two-sided film needs a conforming mesh",
                 "[substructure][Serial][Parallel]")
{
  // A model with a two-sided film crossing the interface on a nonconforming mesh, whose
  // coarse faces have several neighbors, is rejected before the condensation.
  json config = MagnetostaticConfig(1, false, true);
  test::SetModelRun(config, "Offline",
                    (temp_dir / "film_nc" / "environment.model").string());
  IoData iodata(config, false);
  auto mesh = MakeMesh(MakeParitySerial(8, {false, 0, false, true}));
  SubstructuringSolver ss(iodata, mesh);
  CHECK_THROWS_WITH(ss.CondenseEnvironment(),
                    Catch::Matchers::ContainsSubstring("needs a conforming mesh"));
}

}  // namespace palace
