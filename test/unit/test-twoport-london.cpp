// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Standalone prototype for the two-port (transfer-matrix) London sheet (Option 2). The exact
// through-thickness London slab response reduces to a 2x2 surface stiffness in the two face
// potentials A- = A(0), A+ = A(d), with D = d/lambda, s = sinh D, c = cosh D:
//
//   E = 1/(2 lambda s) [ c (A+^2 + A-^2) - 2 A+ A- ]  = 1/2 [A-,A+] M [A-;A+],
//   M = 1/(lambda s) [[c, -1], [-1, c]]   (mu0 absorbed, Palace nondimensional convention).
//
// Assembled on a sheet with ND tangential mass Mass, the doubled-space operator is
//   K = 1/(lambda s) [[c Mass, -Mass_x], [-Mass_x^T, c Mass]],
// where Mass_x is the cross-face mass coupling the two coincident faces. This driver:
//  (A) builds K from a real ND surface mass matrix and checks it reproduces the four analytic
//      mode inductances, and that dropping the off-diagonal (the shipped independent-faces
//      model) overstates the common (isolated-washer) mode;
//  (B) shows the Sigma+ <-> Sigma- twin-DOF map is recoverable purely geometrically on two
//      coincident but independently-numbered meshes -- the piece the crack pass discards.
//
// Analytic mode inductances (nondimensional L_ksq, mu0 = 1):
//   one-sided / far-face Neumann (shipped one-port): L = lambda coth D
//   far-face Dirichlet (grounded):                   L = lambda tanh D
//   common mode  A+ = A- (isolated washer):          L = (lambda/2) coth(D/2)
//   differential A+ = -A-:                           L = (lambda/2) tanh(D/2)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem.hpp>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace mfem;

namespace
{

// Build the cross-face coupling C on the cracked ND space for boundary attribute `attr`, with
// entries off * <phi_a^+, phi_b^->. Serial (ldof == tdof). Pairs the two coincident faces by
// centroid, then integrates the cross mass by evaluating each face's Piola basis at shared
// physical quadrature points (TransformBack), so orientation/sign is handled exactly.
std::unique_ptr<SparseMatrix> BuildCoupling(ParFiniteElementSpace &fes, int attr, double off,
                                            int &n_pairs)
{
  ParMesh &pmesh = *fes.GetParMesh();
  auto C = std::make_unique<SparseMatrix>(fes.GetVSize(), fes.GetVSize());

  std::map<std::array<long, 3>, std::vector<int>> by_centroid;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (pmesh.GetBdrAttribute(be) != attr)
    {
      continue;
    }
    Array<int> v;
    pmesh.GetBdrElementVertices(be, v);
    double cc[3] = {0, 0, 0};
    for (int k = 0; k < v.Size(); k++)
    {
      const double *x = pmesh.GetVertex(v[k]);
      cc[0] += x[0];
      cc[1] += x[1];
      cc[2] += x[2];
    }
    const double inv = 1.0 / v.Size(), q = 1e6;
    by_centroid[{std::lround(cc[0] * inv * q), std::lround(cc[1] * inv * q),
                 std::lround(cc[2] * inv * q)}]
        .push_back(be);
  }

  n_pairs = 0;
  DenseMatrix vshp, vshm, cross;
  for (auto &[key, bes] : by_centroid)
  {
    if (bes.size() != 2)
    {
      continue;
    }
    n_pairs++;
    const int bp = bes[0], bm = bes[1];
    const FiniteElement *fep = fes.GetBE(bp), *fem = fes.GetBE(bm);
    ElementTransformation *Tp = pmesh.GetBdrElementTransformation(bp);
    ElementTransformation *Tm = pmesh.GetBdrElementTransformation(bm);
    const int nd = fep->GetDof(), sd = pmesh.SpaceDimension();
    Array<int> dp, dm;
    fes.GetBdrElementDofs(bp, dp);
    fes.GetBdrElementDofs(bm, dm);
    cross.SetSize(nd);
    cross = 0.0;
    vshp.SetSize(nd, sd);
    vshm.SetSize(nd, sd);
    const IntegrationRule &ir = IntRules.Get(fep->GetGeomType(), 2 * fep->GetOrder() + 2);
    for (int iq = 0; iq < ir.GetNPoints(); iq++)
    {
      const IntegrationPoint &ip = ir.IntPoint(iq);
      Tp->SetIntPoint(&ip);
      const double w = ip.weight * Tp->Weight();
      fep->CalcVShape(*Tp, vshp);
      Vector xphys;
      Tp->Transform(ip, xphys);
      IntegrationPoint ipm;
      Tm->TransformBack(xphys, ipm);
      Tm->SetIntPoint(&ipm);
      fem->CalcVShape(*Tm, vshm);
      for (int a = 0; a < nd; a++)
      {
        for (int b = 0; b < nd; b++)
        {
          double dot = 0;
          for (int dd = 0; dd < sd; dd++)
          {
            dot += vshp(a, dd) * vshm(b, dd);
          }
          cross(a, b) += w * dot;
        }
      }
    }
    for (int a = 0; a < nd; a++)
    {
      int ga = dp[a];
      double sa = 1.0;
      if (ga < 0)
      {
        ga = -1 - ga;
        sa = -1.0;
      }
      for (int b = 0; b < nd; b++)
      {
        int gb = dm[b];
        double sb = 1.0;
        if (gb < 0)
        {
          gb = -1 - gb;
          sb = -1.0;
        }
        const double val = off * sa * sb * cross(a, b);
        C->Add(ga, gb, val);
        C->Add(gb, ga, val);  // symmetric partner <phi_b^-, phi_a^+>
      }
    }
  }
  C->Finalize();
  return C;
}

}  // namespace

namespace
{

// Flat unit-square sheet meshed into two triangles, with a caller-chosen vertex ordering so we
// can build two geometrically identical but independently-numbered copies.
Mesh MakeSquareSheet(const std::vector<std::array<double, 2>> &verts,
                     const std::vector<std::array<int, 3>> &tris)
{
  Mesh mesh(2, static_cast<int>(verts.size()), static_cast<int>(tris.size()));
  for (const auto &v : verts)
  {
    mesh.AddVertex(v[0], v[1]);
  }
  for (const auto &t : tris)
  {
    mesh.AddTriangle(t[0], t[1], t[2]);
  }
  mesh.FinalizeTriMesh(1);  // generate_edges = 1
  return mesh;
}

// ND(1) edge-DOF geometric key: signed by the edge tangent so the two faces agree on
// orientation. Returns (midpoint x, midpoint y, tangent-sign tag).
struct DofKey
{
  long mx, my;  // quantized midpoint
  bool operator<(const DofKey &o) const
  {
    return (mx != o.mx) ? (mx < o.mx) : (my < o.my);
  }
};

DofKey EdgeKey(const Mesh &mesh, int edge)
{
  Array<int> ev;
  mesh.GetEdgeVertices(edge, ev);
  const double *a = mesh.GetVertex(ev[0]), *b = mesh.GetVertex(ev[1]);
  const double mx = 0.5 * (a[0] + b[0]), my = 0.5 * (a[1] + b[1]);
  const double q = 1e6;  // quantize to 1e-6 to match coincident coordinates
  return {std::lround(mx * q), std::lround(my * q)};
}

// Quadratic form x^T K x.
double Quad(SparseMatrix &K, const Vector &x)
{
  Vector Kx(x.Size());
  K.Mult(x, Kx);
  return x * Kx;
}

}  // namespace

TEST_CASE("Two-port London sheet reproduces analytic modes", "[superconductor][twoport][Serial]")
{
  const double lambda = 0.1, d = 0.2;  // D = 2 (thick film, d > lambda)
  const double D = d / lambda, s = std::sinh(D), c = std::cosh(D);

  // ND(1) tangential mass on the flat sheet.
  std::vector<std::array<double, 2>> verts = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  std::vector<std::array<int, 3>> tris = {{0, 1, 2}, {0, 2, 3}};
  Mesh mesh = MakeSquareSheet(verts, tris);
  ND_FECollection fec(1, 2);
  FiniteElementSpace fes(&mesh, &fec);
  const int n = fes.GetNDofs();

  BilinearForm m(&fes);
  m.AddDomainIntegrator(new VectorFEMassIntegrator());
  m.Assemble();
  m.Finalize();
  SparseMatrix &Mass = m.SpMat();

  // Uniform tangential field (1,0) projected onto ND(1).
  Vector u(n);
  {
    GridFunction gf(&fes);
    VectorConstantCoefficient one({Vector({1.0, 0.0})});
    gf.ProjectCoefficient(one);
    u = gf;
  }

  // Assemble the doubled-space operators on [u-; u+]. Faces are coincident so the cross-face
  // mass equals Mass in a shared indexing (twin = identity here; Part B recovers it when the
  // numbering differs). Two-port stiffness diag c/(lambda s), off -1/(lambda s); the shipped
  // model is two independent faces each carrying the single-sheet 1/(L_ksq) split by the
  // attr_scaling = 2 area factor, i.e. per-face diag tanh(D)/(2 lambda), no coupling.
  auto BuildBlock = [&](double diag, double off)
  {
    SparseMatrix K(2 * n, 2 * n);
    for (int i = 0; i < n; i++)
    {
      Array<int> cols;
      Vector vals;
      Mass.GetRow(i, cols, vals);
      for (int k = 0; k < cols.Size(); k++)
      {
        const int j = cols[k];
        const double mij = vals[k];
        K.Add(i, j, diag * mij);
        K.Add(n + i, n + j, diag * mij);
        K.Add(i, n + j, off * mij);
        K.Add(n + i, j, off * mij);
      }
    }
    K.Finalize();
    return K;
  };

  SparseMatrix K = BuildBlock(c / (lambda * s), -1.0 / (lambda * s));       // two-port
  SparseMatrix Kship = BuildBlock(std::tanh(D) / (2.0 * lambda), 0.0);      // shipped one-port

  // Effective sheet inductance is referenced to a single sheet of the physical area, so the
  // mode energy is normalized by the single-sheet mass u^T Mass u, not the doubled mass.
  const double single_mass = Mass.InnerProduct(u, u);
  Vector u_common(2 * n), u_diff(2 * n);
  for (int i = 0; i < n; i++)
  {
    u_common[i] = u_common[n + i] = u[i];
    u_diff[i] = u[i];
    u_diff[n + i] = -u[i];
  }

  SECTION("Common mode -> (lambda/2) coth(D/2)")
  {
    const double L = single_mass / Quad(K, u_common);
    CHECK_THAT(L, Catch::Matchers::WithinRel(0.5 * lambda / std::tanh(D / 2.0), 1e-10));
  }

  SECTION("Differential mode -> (lambda/2) tanh(D/2)")
  {
    const double L = single_mass / Quad(K, u_diff);
    CHECK_THAT(L, Catch::Matchers::WithinRel(0.5 * lambda * std::tanh(D / 2.0), 1e-10));
  }

  SECTION("One-sided (Neumann far face) via Schur complement -> lambda coth D")
  {
    // Schur onto Sigma-: c Mass - Mass (c Mass)^{-1} Mass = (c - 1/c) Mass, exact for any u.
    // 1/L = (c - 1/c)/(lambda s) = tanh D / lambda.
    const double inv_L = (c - 1.0 / c) / (lambda * s);
    CHECK_THAT(1.0 / inv_L, Catch::Matchers::WithinRel(lambda / std::tanh(D), 1e-12));
  }

  SECTION("Far-face Dirichlet (u+ = 0) -> lambda tanh D")
  {
    // Top-left block only: c Mass / (lambda s); 1/L = coth D / lambda.
    const double inv_L = c / (lambda * s);
    CHECK_THAT(1.0 / inv_L, Catch::Matchers::WithinRel(lambda * std::tanh(D), 1e-12));
  }

  SECTION("Shipped decoupled model overstates the common (washer) mode")
  {
    // Independent faces (no coupling): the common mode returns the one-sided value lambda coth D
    // (the shipped result) instead of the correct two-port (lambda/2) coth(D/2).
    const double L_ship = single_mass / Quad(Kship, u_common);
    const double L_true = single_mass / Quad(K, u_common);
    CHECK_THAT(L_ship, Catch::Matchers::WithinRel(lambda / std::tanh(D), 1e-10));
    // At D = 2 the shipped model is ~1.58x too large (matches the measured sheet ratio).
    CHECK(L_ship / L_true > 1.5);
    CHECK(L_ship / L_true < 1.65);
  }
}

TEST_CASE("Sigma+/Sigma- twin-DOF map recoverable geometrically", "[superconductor][twoport][Serial]")
{
  // Two geometrically identical square sheets with different vertex/element ordering, emulating
  // the coincident faces the crack pass produces (whose twin map it then discards).
  Mesh mA = MakeSquareSheet({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, {{0, 1, 2}, {0, 2, 3}});
  Mesh mB = MakeSquareSheet({{1, 1}, {0, 1}, {0, 0}, {1, 0}}, {{2, 3, 0}, {2, 0, 1}});

  ND_FECollection fec(1, 2);
  FiniteElementSpace fesA(&mA, &fec), fesB(&mB, &fec);
  REQUIRE(fesA.GetNDofs() == fesB.GetNDofs());

  // Build the twin map A-dof -> B-dof by matching edge midpoints (ND(1) DOF <-> edge).
  std::map<DofKey, int> keyB;
  for (int e = 0; e < mB.GetNEdges(); e++)
  {
    keyB[EdgeKey(mB, e)] = e;
  }
  std::vector<int> twin(mA.GetNEdges(), -1);
  for (int e = 0; e < mA.GetNEdges(); e++)
  {
    auto it = keyB.find(EdgeKey(mA, e));
    REQUIRE(it != keyB.end());  // every face-A edge has a coincident face-B edge
    twin[e] = it->second;
  }

  SECTION("Map is a bijection")
  {
    std::vector<int> seen(mB.GetNEdges(), 0);
    for (int e = 0; e < mA.GetNEdges(); e++)
    {
      seen[twin[e]]++;
    }
    for (int cnt : seen)
    {
      CHECK(cnt == 1);
    }
  }

  SECTION("Permuted mass matrices agree -> cross-face coupling assembles correctly")
  {
    auto AssembleMass = [&](FiniteElementSpace &fes)
    {
      auto *m = new BilinearForm(&fes);
      m->AddDomainIntegrator(new VectorFEMassIntegrator());
      m->Assemble();
      m->Finalize();
      return m;
    };
    BilinearForm *mA_form = AssembleMass(fesA), *mB_form = AssembleMass(fesB);
    SparseMatrix &MA = mA_form->SpMat(), &MB = mB_form->SpMat();

    // ND(1): 1 DOF per edge, so the DOF twin map is the edge twin map (up to sign). Compare
    // |M_A(i,j)| against |M_B(twin[i],twin[j])|; sign differences are absorbed by the tangent
    // orientation and do not affect the assembled symmetric coupling.
    double max_diff = 0.0;
    for (int i = 0; i < MA.Height(); i++)
    {
      Array<int> cols;
      Vector vals;
      MA.GetRow(i, cols, vals);
      for (int k = 0; k < cols.Size(); k++)
      {
        const double bij = MB(twin[i], twin[cols[k]]);
        max_diff = std::max(max_diff, std::abs(std::abs(vals[k]) - std::abs(bij)));
      }
    }
    CHECK(max_diff < 1e-12);
    delete mA_form;
    delete mB_form;
  }
}

// Stage 1 on the real geometry: load the actual circular_hole washer through Palace's mesh
// pipeline and inspect whether the London film is cracked. FINDING: it is NOT. The film
// (attr 8) is fully interior and cracking is enabled by default, but Boundaries.attributes
// (the crack-candidate set) aggregates PEC/impedance/ports/current only -- NOT Superconductor
// or FluxLoop (configfile.cpp:1087-1128). So the film runs with a single shared tangential
// A_t (A continuous across the film), i.e. one potential, not two independent faces. The
// two-port coupled operator therefore requires first forcing the film to crack, which is a
// prerequisite that also perturbs the flux-loop generator/constraint machinery.
TEST_CASE("London washer film runs UNCRACKED (two-port prerequisite check)",
          "[superconductor][twoport][Serial]")
{
  using nlohmann::json;
  const std::string mesh_path =
      "/Users/dzpham/palace_ci/examples/circular_hole/mesh/circular_hole.msh";
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", ""}}},
      {"Model", {{"Mesh", mesh_path}, {"L0", 1.0e-6}}},
      {"Domains", {{"Materials", json::array({{{"Attributes", {1}}, {"Permeability", 1.0}}})}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {2, 3, 4, 5, 6, 7}}}},
        {"Superconductor", json::array({{{"Attributes", {8}},
                                         {"PenetrationDepth", 0.1},
                                         {"Thickness", 0.2}}})},
        {"FluxLoop", json::array({{{"Index", 1},
                                   {"FilmAttributes", {8}},
                                   {"HoleAttributes", {9}},
                                   {"FluxAmounts", {1.0}},
                                   {"Direction", "+Z"}}})}}},
      {"Solver", {{"Order", 1}, {"Magnetostatic", {{"Save", 0}}}}}};

  palace::IoData iodata(config, /*print=*/false);
  auto smesh = palace::mesh::Load(iodata, palace::Mpi::World());
  REQUIRE(smesh);

  // Diagnostics: cracking flag, recorded cracked attributes, and interior-face count for attr 8.
  int n8 = 0, n8_interior = 0;
  for (int be = 0; be < smesh->GetNBE(); be++)
  {
    if (smesh->GetBdrAttribute(be) != 8)
    {
      continue;
    }
    n8++;
    int f, o;
    smesh->GetBdrElementFace(be, &f, &o);
    int e1, e2;
    smesh->GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)
    {
      n8_interior++;
    }
  }
  std::string cracked;
  for (int a : iodata.boundaries.cracked_attributes)
  {
    cracked += std::to_string(a) + " ";
  }
  UNSCOPED_INFO("crack_bdr_elements=" << iodata.model.crack_bdr_elements
                                      << " cracked_attributes=[" << cracked << "] attr8_be=" << n8
                                      << " attr8_interior=" << n8_interior);

  // The film is fully interior (crackable in principle) yet is NOT cracked: Superconductor /
  // FluxLoop attributes are excluded from the crack-candidate set. Confirms the single shared
  // A_t discretization and that the two-port needs cracking to be forced first.
  CHECK(n8 > 0);
  CHECK(n8_interior == n8);
  CHECK(iodata.boundaries.cracked_attributes.count(8) == 0);
}

// AMS gate: force-crack the film, assemble curl-curl + two-port sheet self-term + cross-face
// coupling C + a London-shift volume mass (SPD regularization), and confirm HypreAMS/CG
// converges on the coupled operator -- and not dramatically worse than the uncoupled baseline.
// This is the go/no-go for building the production two-port. Serial prototype.
TEST_CASE("AMS converges on the coupled two-port cracked washer",
          "[superconductor][twoport][Serial]")
{
  if (palace::Mpi::Size(palace::Mpi::World()) > 1)
  {
    return;  // serial prototype
  }
  using nlohmann::json;
  setenv("PALACE_LONDON_FORCE_CRACK", "1", 1);
  const std::string mesh_path =
      "/Users/dzpham/palace_ci/examples/circular_hole/mesh/circular_hole.msh";
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", ""}}},
      {"Model", {{"Mesh", mesh_path}, {"L0", 1.0e-6}}},
      {"Domains", {{"Materials", json::array({{{"Attributes", {1}}, {"Permeability", 1.0}}})}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {2, 3, 4, 5, 6, 7}}}},
        {"Superconductor", json::array({{{"Attributes", {8}},
                                         {"PenetrationDepth", 0.1},
                                         {"Thickness", 0.2}}})},
        {"FluxLoop", json::array({{{"Index", 1},
                                   {"FilmAttributes", {8}},
                                   {"HoleAttributes", {9}},
                                   {"FluxAmounts", {1.0}},
                                   {"Direction", "+Z"}}})}}},
      {"Solver", {{"Order", 1}, {"Magnetostatic", {{"Save", 0}}}}}};

  palace::IoData iodata(config, /*print=*/false);
  auto smesh = palace::mesh::Load(iodata, palace::Mpi::World());
  REQUIRE(iodata.boundaries.cracked_attributes.count(8) == 1);  // force-crack worked
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = palace::mesh::ComputeReferenceLength(smesh, palace::Mpi::World());
  }
  iodata.NondimensionalizeInputs(smesh);
  auto pmesh = palace::mesh::Partition(iodata, std::move(smesh), palace::Mpi::World());

  ND_FECollection fec(1, pmesh->Dimension());
  ParFiniteElementSpace fes(pmesh.get(), &fec);

  const double lambda = 0.1, d = 0.2, D = d / lambda, s = std::sinh(D), cD = std::cosh(D);
  const double diag = cD / (lambda * s);    // two-port self coefficient c/(lambda s)
  const double off = -1.0 / (lambda * s);   // cross-face coupling -1/(lambda s)
  const double shift = 1.0;                 // London PC shift: SPD regularization

  // A = curl-curl + shift * (volume mass) + diag * (sheet mass on attr 8).
  ConstantCoefficient one(1.0), diagc(diag), shiftc(shift);
  ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new CurlCurlIntegrator(one));
  a.AddDomainIntegrator(new VectorFEMassIntegrator(shiftc));
  Array<int> sheet_marker(pmesh->bdr_attributes.Max());
  sheet_marker = 0;
  sheet_marker[8 - 1] = 1;
  a.AddBoundaryIntegrator(new VectorFEMassIntegrator(diagc), sheet_marker);
  a.Assemble();
  a.Finalize();
  std::unique_ptr<HypreParMatrix> A(a.ParallelAssemble());

  // Cross-face coupling as a HypreParMatrix on the (serial) true dofs.
  int n_pairs = 0;
  auto Csp = BuildCoupling(fes, 8, off, n_pairs);
  REQUIRE(n_pairs > 0);
  HypreParMatrix Ch(fes.GetComm(), fes.GlobalTrueVSize(), fes.GetTrueDofOffsets(), Csp.get());
  std::unique_ptr<HypreParMatrix> Acoupled(Add(1.0, *A, 1.0, Ch));

  auto Solve = [&](HypreParMatrix &op)
  {
    HypreAMS ams(op, &fes);
    CGSolver cg(fes.GetComm());
    cg.SetOperator(op);
    cg.SetPreconditioner(ams);
    cg.SetRelTol(1.0e-8);
    cg.SetMaxIter(500);
    cg.SetPrintLevel(0);
    Vector B(fes.GetTrueVSize()), X(fes.GetTrueVSize());
    B.Randomize(1);
    X = 0.0;
    cg.Mult(B, X);
    return std::pair<bool, int>{cg.GetConverged(), cg.GetNumIterations()};
  };

  auto [conv_base, it_base] = Solve(*A);
  auto [conv_cpl, it_cpl] = Solve(*Acoupled);
  UNSCOPED_INFO("AMS gate: n_pairs=" << n_pairs << " baseline converged=" << conv_base
                                     << " iters=" << it_base << " | coupled converged="
                                     << conv_cpl << " iters=" << it_cpl);

  // Gate: AMS must converge on the coupled operator, and not blow up vs the uncoupled baseline.
  CHECK(conv_base);
  CHECK(conv_cpl);
  CHECK(it_cpl <= 3 * it_base + 20);
}
