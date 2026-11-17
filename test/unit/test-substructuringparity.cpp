// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <filesystem>
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
          {"Environment", {{"Attributes", kEnvironment}}}};
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
// two-sided film crossing the interface.
json MagnetostaticConfig(int order, bool sheets, bool film)
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
  for (int i : zero)
  {
    xx[i] = 0.0;
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

  // The stiffness matrix of the regular solve (curl-curl and sheet terms), the curl-curl
  // term alone (the energy operator), and the regularization eps M.
  auto K = curlcurl_op.GetStiffnessMatrix();
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetCurlCurlInvPermeability());
  BilinearForm kcc(nd_fespace);
  kcc.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  ParOperator Kcc(kcc.Assemble(false), nd_fespace);
  std::vector<int> domain_attrs(kRegion);
  domain_attrs.insert(domain_attrs.end(), kEnvironment.begin(), kEnvironment.end());
  MaterialPropertyCoefficient eps_func(mat_op.MaxCeedAttribute());
  eps_func.AddMaterialProperty(mat_op.GetCeedAttributes(domain_attrs),
                               with_sheets ? kMagRegularizationSheets : kMagRegularization);
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
// own no element of x < 1/2, nor a neighbor of one.
std::vector<std::unique_ptr<Mesh>> MakeMesh(mfem::Mesh &&serial, bool by_side = false)
{
  const int np = Mpi::Size(Mpi::World());
  std::vector<int> part;
  if (by_side && np > 1)
  {
    part.resize(serial.GetNE());
    mfem::Vector c;
    for (int e = 0; e < serial.GetNE(); e++)
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
  // on the mesh of the region, which keeps the environment neighbors of these faces.
  // Not structured bindings, which the lambdas below cannot capture with OpenMP (Clang 18).
  bool region_below, by_side;
  int order;
  std::tie(region_below, order, by_side) =
      GENERATE(table<bool, int, bool>({{true, 1, false},
                                       {false, 1, false},
                                       {true, 2, false},
                                       {false, 2, false},
                                       {true, 1, true},
                                       {false, 1, true}}));
  CAPTURE(region_below, order, by_side);
  const std::string model_path = (temp_dir / "sheet_side.model").string();
  auto make_config = [&](const std::string &mode)
  {
    json config = MagnetostaticConfig(order, true, false);
    config["Boundaries"]["PEC"]["Attributes"] = {1, 2, 9};
    config["Boundaries"]["SurfaceCurrent"] = {
        {{"Index", 1}, {"Attributes", {8}}, {"Direction", "+X"}}};
    auto &sub = config["Solver"]["Substructuring"];
    sub["Region"]["Attributes"] = region_below ? kRegion : kEnvironment;
    sub["Environment"]["Attributes"] = region_below ? kEnvironment : kRegion;
    sub["Mode"] = mode;
    sub["SaveModel"] = model_path;
    return IoData(config, false);
  };
  // The energy of the surface current, and whether the environment was solved.
  auto energy = [&](const IoData &iodata, mfem::Mesh &&serial)
  {
    auto mesh = MakeMesh(std::move(serial), by_side);
    Vector J;
    CurlCurlOperator(iodata, mesh).GetCurrentExcitationVector(1, J);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    const double E = ss.MagnetostaticEnergyMatrix({1}, {J}, {}, {}, {}, nullptr)(0, 0);
    CHECK(ss.EnvironmentFactored() ==
          (iodata.solver.substructuring->mode == SubstructuringMode::OFFLINE));
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
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Substructuring model with a two-sided film needs a conforming mesh",
                 "[substructure][Serial][Parallel]")
{
  // A model with a two-sided film crossing the interface on a nonconforming mesh, whose
  // coarse faces have several neighbors, is rejected before the condensation.
  json config = MagnetostaticConfig(1, false, true);
  config["Solver"]["Substructuring"]["SaveModel"] = (temp_dir / "film_nc.model").string();
  IoData iodata(config, false);
  auto mesh = MakeMesh(MakeParitySerial(8, {false, 0, false, true}));
  SubstructuringSolver ss(iodata, mesh);
  CHECK_THROWS_WITH(ss.CondenseEnvironment(),
                    Catch::Matchers::ContainsSubstring("needs a conforming mesh"));
}

}  // namespace palace
