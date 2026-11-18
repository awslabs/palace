// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/bilinearform.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/rap.hpp"
#include "models/curlcurloperator.hpp"
#include "models/materialoperator.hpp"
#include "models/substructuringsolver.hpp"
#include "models/substructuringsolverimpl.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using json = nlohmann::json;

namespace palace
{

namespace
{

// The magnetostatic regularization eps (1/mu) M / L^2 of SubstructuringSolver (see
// SubstructuringSolver::Impl::kMagRegularization), and its scale 1/L^2, L the largest side
// of the bounding box of the mesh.
constexpr double kMagRegularization = 1.0e-3;

double RegularizationScale(const mfem::ParMesh &mesh)
{
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  bbmax -= bbmin;
  return 1.0 / (bbmax.Max() * bbmax.Max());
}

// Unit-cube ParMesh split at x=0.5 (domain attrs 1/2), boundary attrs 1 (x=0), 2 (x=1),
// 3 (other faces).
std::unique_ptr<mfem::ParMesh> MakeSplitCube(int nx)
{
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
  for (int e = 0; e < serial.GetNE(); e++)
  {
    mfem::Vector c;
    serial.GetElementCenter(e, c);
    serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
  }
  for (int b = 0; b < serial.GetNBE(); b++)
  {
    mfem::Array<int> vtx;
    serial.GetBdrElementVertices(b, vtx);
    double xc = 0.0;
    for (int j = 0; j < vtx.Size(); j++)
    {
      xc += serial.GetVertex(vtx[j])[0];
    }
    xc /= vtx.Size();
    serial.SetBdrAttribute(b, xc < 1e-9 ? 1 : (xc > 1.0 - 1e-9 ? 2 : 3));
  }
  serial.SetAttributes();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// An unstructured tetrahedral cube split by a wavy (non-planar) interface into region (attr
// 1) and environment (attr 2). Tets + a curved interface exercise the DOF maps and
// interface identification on a more realistic mesh than the structured hex half-space
// split, while keeping the x=0 / x=1 terminal-face convention (bdr attr 1 / 2, sides 3).
// With reorder, the same tetrahedra come in another element and vertex order (so other face
// orientations, as another partition gives).
mfem::Mesh MakeWavyTetSerial(int nx, bool reorder = false)
{
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::TETRAHEDRON);
  if (reorder)
  {
    const int nv = serial.GetNV(), step = 7;
    MFEM_VERIFY(std::gcd(nv, step) == 1, "Reordering needs nv coprime with the step!");
    std::vector<int> inv(nv);
    mfem::Mesh perm(3, nv, serial.GetNE(), serial.GetNBE());
    for (int v = 0; v < nv; v++)
    {
      perm.AddVertex(serial.GetVertex((v * step) % nv));
      inv[(v * step) % nv] = v;
    }
    // Reversed element order, and even vertex permutations (which keep the orientation).
    for (int e = serial.GetNE() - 1; e >= 0; e--)
    {
      const int *ev = serial.GetElement(e)->GetVertices();
      perm.AddTet(inv[ev[1]], inv[ev[2]], inv[ev[0]], inv[ev[3]]);
    }
    for (int b = 0; b < serial.GetNBE(); b++)
    {
      const int *bv = serial.GetBdrElement(b)->GetVertices();
      perm.AddBdrTriangle(inv[bv[1]], inv[bv[2]], inv[bv[0]]);
    }
    perm.FinalizeTetMesh(1, 0, true);
    serial = std::move(perm);
  }
  auto iface = [](double y, double z)
  { return 0.5 + 0.15 * std::sin(M_PI * y) * std::sin(M_PI * z); };
  for (int e = 0; e < serial.GetNE(); e++)
  {
    mfem::Vector c;
    serial.GetElementCenter(e, c);
    serial.SetAttribute(e, (c(0) < iface(c(1), c(2))) ? 1 : 2);
  }
  for (int b = 0; b < serial.GetNBE(); b++)
  {
    mfem::Array<int> vtx;
    serial.GetBdrElementVertices(b, vtx);
    double xc = 0.0;
    for (int j = 0; j < vtx.Size(); j++)
    {
      xc += serial.GetVertex(vtx[j])[0];
    }
    xc /= vtx.Size();
    serial.SetBdrAttribute(b, xc < 1e-9 ? 1 : (xc > 1.0 - 1e-9 ? 2 : 3));
  }
  serial.SetAttributes();
  return serial;
}

// With nonconforming, ready for nonconforming refinement, as Palace prepares a mesh for it.
std::unique_ptr<mfem::ParMesh> MakeWavyTetSplit(int nx, bool reorder = false,
                                                bool nonconforming = false)
{
  mfem::Mesh serial = MakeWavyTetSerial(nx, reorder);
  if (nonconforming)
  {
    serial.EnsureNCMesh(true);
  }
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// Optional design features inside the region of MakeGradedSplit: a conductor patch (the
// region terminal, boundary attribute 1, covers only the x = 0 faces with center y <
// patch_y; the rest of that face gets the natural attribute 3) and a dielectric block
// (region hexes with center (x, y, z) < block get element attribute 3).
struct RegionDesign
{
  double patch_y = 1.0;
  std::array<double, 3> block = {0.0, 0.0, 0.0};
  bool sheet = false;  // interior boundary (attribute 4) on the plane y = 0.5, all x
};

// A hex cube whose region (x<0.5, attr 1) and environment (x>0.5, attr 2) have independent
// x-resolutions (a and b cells) but share the same y-z grid (n cells), so the interface at
// x=0.5 has identical nodes regardless of a. Re-meshing the region (varying a, fixing b and
// n) keeps Gamma and the environment fixed -- exactly the offline/online region-redesign
// workflow. `design` adds features inside the region (see RegionDesign).
mfem::Mesh MakeGradedSplitSerial(int a, int b, int n, const RegionDesign &design = {})
{
  std::vector<double> xs;
  for (int i = 0; i <= a; i++)
  {
    xs.push_back(0.5 * i / a);
  }
  for (int i = 1; i <= b; i++)
  {
    xs.push_back(0.5 + 0.5 * i / b);
  }
  const int nx = static_cast<int>(xs.size());
  auto vid = [&](int i, int j, int k) { return (i * (n + 1) + j) * (n + 1) + k; };
  mfem::Mesh serial(3, nx * (n + 1) * (n + 1), (nx - 1) * n * n,
                    2 * n * n + 4 * (nx - 1) * n + (design.sheet ? (nx - 1) * n : 0), 3);
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      for (int k = 0; k <= n; k++)
      {
        serial.AddVertex(xs[i], static_cast<double>(j) / n, static_cast<double>(k) / n);
      }
    }
  }
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int k = 0; k < n; k++)
      {
        int v[8] = {vid(i, j, k),
                    vid(i + 1, j, k),
                    vid(i + 1, j + 1, k),
                    vid(i, j + 1, k),
                    vid(i, j, k + 1),
                    vid(i + 1, j, k + 1),
                    vid(i + 1, j + 1, k + 1),
                    vid(i, j + 1, k + 1)};
        const double cx = 0.5 * (xs[i] + xs[i + 1]), cy = (j + 0.5) / n, cz = (k + 0.5) / n;
        const bool in_block =
            cx < design.block[0] && cy < design.block[1] && cz < design.block[2];
        serial.AddHex(v, (cx < 0.5) ? (in_block ? 3 : 1) : 2);
      }
    }
  }
  for (int j = 0; j < n; j++)
  {
    for (int k = 0; k < n; k++)
    {
      int q0[4] = {vid(0, j, k), vid(0, j + 1, k), vid(0, j + 1, k + 1), vid(0, j, k + 1)};
      serial.AddBdrQuad(q0, ((j + 0.5) / n < design.patch_y) ? 1 : 3);
      int q1[4] = {vid(nx - 1, j, k), vid(nx - 1, j, k + 1), vid(nx - 1, j + 1, k + 1),
                   vid(nx - 1, j + 1, k)};
      serial.AddBdrQuad(q1, 2);
    }
  }
  for (int i = 0; i < nx - 1; i++)
  {
    for (int k = 0; k < n; k++)
    {
      int y0[4] = {vid(i, 0, k), vid(i, 0, k + 1), vid(i + 1, 0, k + 1), vid(i + 1, 0, k)};
      serial.AddBdrQuad(y0, 3);
      int y1[4] = {vid(i, n, k), vid(i + 1, n, k), vid(i + 1, n, k + 1), vid(i, n, k + 1)};
      serial.AddBdrQuad(y1, 3);
    }
  }
  for (int i = 0; i < nx - 1; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int z0[4] = {vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0)};
      serial.AddBdrQuad(z0, 3);
      int z1[4] = {vid(i, j, n), vid(i, j + 1, n), vid(i + 1, j + 1, n), vid(i + 1, j, n)};
      serial.AddBdrQuad(z1, 3);
    }
  }
  if (design.sheet)
  {
    MFEM_VERIFY(n % 2 == 0, "The interior sheet needs an even n!");
    for (int i = 0; i < nx - 1; i++)
    {
      for (int k = 0; k < n; k++)
      {
        const int j = n / 2;
        int q[4] = {vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j, k + 1), vid(i, j, k + 1)};
        serial.AddBdrQuad(q, 4);
      }
    }
  }
  serial.FinalizeHexMesh(1, 0, true);
  serial.SetAttributes();
  return serial;
}

// With nonconforming, ready for nonconforming refinement.
std::unique_ptr<mfem::ParMesh> MakeGradedSplit(int a, int b, int n,
                                               const RegionDesign &design = {},
                                               bool nonconforming = false)
{
  mfem::Mesh serial = MakeGradedSplitSerial(a, b, n, design);
  if (nonconforming)
  {
    serial.EnsureNCMesh();  // a ParMesh cannot be made nonconforming afterwards
  }
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

// The energy 1/2 u^T K u of a field with the energy operator of the whole domain.
double Energy(const SubstructuringSolver &ss, const Vector &u)
{
  return 0.5 * ss.MutualEnergy(u, u);
}

// The energy matrices of London flux states and of surface currents alone.
mfem::DenseMatrix SheetEnergies(SubstructuringSolver &ss, const std::vector<int> &ids,
                                const std::vector<Vector> &a,
                                std::vector<Vector> *fields = nullptr, int n_fields = 0,
                                std::vector<double> *volume_energy = nullptr)
{
  return ss.MagnetostaticEnergyMatrix({}, {}, ids, a, {}, nullptr, fields, n_fields,
                                      nullptr, volume_energy);
}
mfem::DenseMatrix CurrentEnergies(SubstructuringSolver &ss, const std::vector<int> &ids,
                                  const std::vector<Vector> &J,
                                  std::vector<Vector> *fields = nullptr, int n_fields = 0)
{
  return ss.MagnetostaticEnergyMatrix(ids, J, {}, {}, {}, nullptr, fields, n_fields);
}

// The closed current J = curl(0, 0, psi), psi = (x(x1-x) y(1-y) z(1-z))^2 for x < x1 and
// zero beyond (zero normal component on the faces of the unit cube and on x = x1), as a
// dual vector of the Nédélec space of the given order.
Vector ClosedCurrent(mfem::ParMesh &mesh, int order, double x1 = 1.0)
{
  mfem::VectorFunctionCoefficient jc(3,
                                     [x1](const mfem::Vector &x, mfem::Vector &j)
                                     {
                                       const double a = (x(0) < x1) ? x(0) * (x1 - x(0))
                                                                    : 0.0,
                                                    b = x(1) * (1.0 - x(1)),
                                                    c = x(2) * (1.0 - x(2));
                                       j(0) = 2.0 * a * a * b * (1.0 - 2.0 * x(1)) * c * c;
                                       j(1) = -2.0 * a * (x1 - 2.0 * x(0)) * b * b * c * c;
                                       j(2) = 0.0;
                                     });
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&mesh, &fec);
  mfem::ParLinearForm lf(&pfes);
  auto *lfi = new mfem::VectorFEDomainLFIntegrator(jc);
  lfi->SetIntRule(&mfem::IntRules.Get(mesh.GetTypicalElementGeometry(), 14));
  lf.AddDomainIntegrator(lfi);
  lf.Assemble();
  Vector J(pfes.GetTrueVSize());
  lf.ParallelAssemble(J);
  return J;
}
Vector ClosedCurrent(Mesh &mesh, int order, double x1 = 1.0)
{
  return ClosedCurrent(mesh.Get(), order, x1);
}

}  // namespace

TEST_CASE("SubstructuringSolver reproduces full-domain electrostatics",
          "[substructure][Serial][Parallel]")
{
  auto run = [](double eps_r, double eps_e, int order,
                const std::function<std::unique_ptr<mfem::ParMesh>()> &make_mesh)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", eps_r}},
            {{"Attributes", {2}}, {"Permittivity", eps_e}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);

    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(make_mesh()));

    // Region-condensed solve (parent true-DOF solution).
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    Vector u = ss.SolveExcitation(ss.TerminalIndices()[0]);

    // Reuse: a second solve reuses the materialized environment DtN (no re-condensation)
    // and must give an identical result.
    Vector u2 = ss.SolveExcitation(ss.TerminalIndices()[0]);
    {
      Vector d(u2);
      d -= u;
      CHECK(d.Norml2() <= 1.0e-12 * (u.Norml2() + 1.0e-30));
    }

    // Full-domain reference on the same parent space, parallel CG.
    auto &pmesh = mesh.back()->Get();
    mfem::H1_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    const int max_attr = pmesh.attributes.Max();
    mfem::Vector eps_by_attr(max_attr);
    eps_by_attr = 1.0;
    for (const auto &m : iodata.domains.materials)
    {
      for (int a : m.attributes)
      {
        eps_by_attr(a - 1) = m.epsilon_r.s[0];
      }
    }
    mfem::PWConstCoefficient eps(eps_by_attr);
    mfem::ParBilinearForm a(&pfes);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
    a.Assemble();
    mfem::ParGridFunction xgf(&pfes);
    xgf = 0.0;
    const int maxb = pmesh.bdr_attributes.Max();
    mfem::Array<int> m1(maxb), m2(maxb);
    m1 = 0;
    m2 = 0;
    m1[0] = 1;  // attr 1 -> 1 V
    m2[1] = 1;  // attr 2 -> 0 V
    mfem::ConstantCoefficient one(1.0), zero(0.0);
    xgf.ProjectBdrCoefficient(one, m1);
    xgf.ProjectBdrCoefficient(zero, m2);
    mfem::Array<int> ess_bdr(maxb), ess_tdofs;
    ess_bdr = 0;
    ess_bdr[0] = 1;
    ess_bdr[1] = 1;
    pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    mfem::ParLinearForm bform(&pfes);
    bform = 0.0;
    bform.Assemble();
    mfem::OperatorPtr A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdofs, xgf, bform, A, X, B);
    mfem::HypreParMatrix *Ah = A.As<mfem::HypreParMatrix>();
    mfem::HypreBoomerAMG amg(*Ah);
    amg.SetPrintLevel(0);
    mfem::HyprePCG pcg(*Ah);
    pcg.SetTol(1e-12);
    pcg.SetMaxIter(500);
    pcg.SetPrintLevel(0);
    pcg.SetPreconditioner(amg);
    pcg.Mult(B, X);
    a.RecoverFEMSolution(X, bform, xgf);
    mfem::Vector u_full;
    xgf.GetTrueDofs(u_full);

    // Compare the full reconstructed field (region + recovered environment) to the monolith
    // on all true DOFs.
    double num = 0.0, den = 0.0;
    for (int i = 0; i < pfes.GetTrueVSize(); i++)
    {
      double e = u(i) - u_full(i);
      num += e * e;
      den += u_full(i) * u_full(i);
    }
    double gnum = 0.0, gden = 0.0;
    MPI_Allreduce(&num, &gnum, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
    MPI_Allreduce(&den, &gden, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
    CHECK(std::sqrt(gnum / gden) < 1.0e-8);

    // Electrostatic energy (QoI) must match the monolith. Monolith energy = 1/2 X^T (B -
    // r), computed here directly as 1/2 u_full^T A_full u_full via the substructuring
    // accessor on the reference field.
    const double e_sub = Energy(ss, u);
    const double e_ref = Energy(ss, u_full);
    CHECK(std::abs(e_sub - e_ref) <= 1.0e-8 * std::abs(e_ref));

    // Capacitance sweep: solve each terminal excitation and form C_ij = phi_i^T K phi_j.
    // The Maxwell capacitance matrix must be symmetric and, with no grounded conductor, its
    // rows must sum to zero (K annihilates the constant vector).
    const std::vector<int> terms = ss.TerminalIndices();
    REQUIRE(terms.size() == 2);
    std::vector<Vector> phi(terms.size());
    for (std::size_t j = 0; j < terms.size(); j++)
    {
      phi[j] = ss.SolveExcitation(terms[j]);
    }
    const double c00 = ss.MutualEnergy(phi[0], phi[0]);
    const double c01 = ss.MutualEnergy(phi[0], phi[1]);
    const double c10 = ss.MutualEnergy(phi[1], phi[0]);
    const double c11 = ss.MutualEnergy(phi[1], phi[1]);
    CHECK(std::abs(c01 - c10) <= 1.0e-9 * std::abs(c00));
    CHECK(std::abs(c00 + c01) <= 1.0e-8 * std::abs(c00));
    CHECK(std::abs(c11 + c10) <= 1.0e-8 * std::abs(c11));
    // u drives the lowest terminal, matching phi[0].
    CHECK(std::abs(2.0 * e_sub - c00) <= 1.0e-9 * std::abs(c00));
  };

  SECTION("uniform permittivity, order 1")
  {
    run(1.0, 1.0, 1, [] { return MakeSplitCube(6); });
  }
  SECTION("contrast across interface, order 1")
  {
    run(1.0, 10.0, 1, [] { return MakeSplitCube(6); });
  }
  SECTION("contrast across interface, order 2")
  {
    run(1.0, 10.0, 2, [] { return MakeSplitCube(6); });
  }
  SECTION("contrast across interface, order 3")
  {
    run(1.0, 10.0, 3, [] { return MakeSplitCube(6); });
  }
  // Unstructured tets + a curved (wavy) interface, vs the monolith.
  SECTION("wavy tet interface, contrast, order 1")
  {
    run(1.0, 10.0, 1, [] { return MakeWavyTetSplit(8); });
  }
  SECTION("wavy tet interface, contrast, order 2")
  {
    run(1.0, 10.0, 2, [] { return MakeWavyTetSplit(8); });
  }
}

TEST_CASE("SubstructuringSolver rejects attribute sets that do not split the mesh",
          "[substructure][Serial][Parallel]")
{
  auto make = [](const std::vector<int> &region, const std::vector<int> &environment)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1, 2}}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Substructuring",
           {{"Region", {{"Attributes", region}}},
            {"Environment", {{"Attributes", environment}}}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(2)));
    SubstructuringSolver ss(iodata, mesh);
  };
  // An attribute in neither set would be left out of both operators.
  CHECK_THROWS_WITH(make({1}, {3}),
                    Catch::Matchers::ContainsSubstring("Domain attribute 2 of the mesh is "
                                                       "in neither"));
  CHECK_THROWS_WITH(make({1, 4}, {2}), Catch::Matchers::ContainsSubstring(
                                           "\"Region\" attribute 4 is not a domain"));
  CHECK_NOTHROW(make({1}, {2}));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver mesh of the region takes no a priori refinement",
                 "[substructure][Serial][Parallel][GPU]")
{
  // The offline run writes the mesh of the region refined: an online run on it with a
  // priori refinement would refine the region again, its interface included, so it is
  // rejected (on the full mesh, both runs refine alike).
  const std::string model_path = (temp_dir / "substruct_region_refined.bin").string();
  auto make_config =
      [&](const std::string &mode, const json &refinement, bool ground = false)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}, {"Refinement", refinement}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    if (ground)
    {
      config["Boundaries"]["Ground"] = {{"Attributes", {3}}};
    }
    return IoData(config, false);
  };
  auto condense = [](const IoData &iodata, std::unique_ptr<mfem::ParMesh> &&pmesh)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  };
  auto region_mesh = [&]()
  {
    mfem::Mesh serial(RegionMeshPath(model_path).c_str());
    return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
  };
  const json none = json::object(), uniform = {{"UniformLevels", 1}},
             boxes = {{"Boxes",
                       {{{"Levels", 1},
                         {"BoundingBoxMin", {0.0, 0.0, 0.0}},
                         {"BoundingBoxMax", {0.25, 1.0, 1.0}}}}}};
  condense(make_config("Offline", none), MakeSplitCube(4));
  CHECK_NOTHROW(condense(make_config("Online", none), region_mesh()));
  CHECK_NOTHROW(condense(make_config("Online", uniform), MakeSplitCube(4)));
  for (const json &refinement : {uniform, boxes})
  {
    CHECK_THROWS_WITH(condense(make_config("Online", refinement), region_mesh()),
                      Catch::Matchers::ContainsSubstring("is written refined"));
  }

  // With the ground reaching the interface, environment elements are written with the
  // region: once refined, they are told from a changed environment by the moments that
  // refinement keeps.
  auto refined = [](std::unique_ptr<mfem::ParMesh> &&pmesh)
  {
    pmesh->UniformRefinement();
    return std::move(pmesh);
  };
  condense(make_config("Offline", none, true), MakeSplitCube(4));
  CHECK_THROWS_WITH(condense(make_config("Online", uniform, true), refined(region_mesh())),
                    Catch::Matchers::ContainsSubstring("is written refined"));
  CHECK_THROWS_WITH(
      condense(make_config("Online", uniform, true), refined(MakeSplitCube(4))),
      Catch::Matchers::ContainsSubstring("The environment differs"));
}

#if defined(MFEM_USE_STRUMPACK)
TEST_CASE("SubstructuringSolver factors without STRUMPACK compression",
          "[substructure][Serial][Parallel]")
{
  // The factorizations of substructuring are exact, as no Krylov iteration corrects the
  // error of a compressed factor: the block low-rank compression of Solver.Linear leaves
  // the capacitance unchanged (on a mesh with fronts large enough to compress).
  auto capacitance = [](const std::string &compression)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 2},
          {"Linear",
           {{"Type", "STRUMPACK"},
            {"STRUMPACKCompressionType", compression},
            {"STRUMPACKCompressionTol", 1.0e-2}}},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(16)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  const mfem::DenseMatrix C0 = capacitance("None"), C1 = capacitance("BLR");
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CAPTURE(i, j, C0(i, j), C1(i, j));
      CHECK(std::abs(C1(i, j) - C0(i, j)) <= 1.0e-12 * std::abs(C0(i, i)));
    }
  }
}
#endif

TEST_CASE("SubstructuringSolver rejects unknown terminal attributes",
          "[substructure][Serial][Parallel]")
{
  // As in the regular solve, every terminal attribute is a boundary attribute of the mesh.
  auto make = [](const std::vector<int> &attributes)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1, 2}}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", attributes}},
            {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(2)));
    SubstructuringSolver ss(iodata, mesh);
  };
  CHECK_THROWS_WITH(make({1, 12}), Catch::Matchers::ContainsSubstring(
                                       "Unknown terminal boundary attribute 12"));
  CHECK_NOTHROW(make({1}));
}

TEST_CASE("SubstructuringSolver anisotropic permittivity",
          "[substructure][Serial][Parallel]")
{
  // The cube is layered along x with terminals at x = 0 and x = 1, so the field is purely
  // x-directed and only the xx permittivity component matters. A diagonal anisotropic
  // tensor diag(eps_x, *, *) must therefore give the same capacitance as the scalar-eps_x
  // series capacitor: C11 = 1 / (0.5/eps_x_region + 0.5/eps_x_env).
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", {1.0, 5.0, 7.0}}},
          {{"Attributes", {2}}, {"Permittivity", {10.0, 3.0, 2.0}}}}}}},
      {"Boundaries",
       {{"Terminal",
         {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  Vector phi0 = ss.SolveExcitation(1);
  Vector phi1 = ss.SolveExcitation(2);
  const double c00 = ss.MutualEnergy(phi0, phi0);
  const double c01 = ss.MutualEnergy(phi0, phi1);
  const double c_series = 1.0 / (0.5 / 1.0 + 0.5 / 10.0);
  CHECK(std::abs(c00 - c_series) <= 1.0e-6 * c_series);
  CHECK(std::abs(c00 + c01) <= 1.0e-8 * c_series);
}

TEST_CASE("SubstructuringSolver reproduces full-domain magnetostatic energy",
          "[substructure][Serial][Parallel]")
{
  // Gauge-free H(curl) magnetostatics: the fields are gauge-ambiguous, but the magnetic
  // energy is gauge-invariant. A closed current (divergence-free, with no normal component
  // on the boundary) is a valid magnetostatic excitation; the region-condensed energy must
  // match a monolith of the same regularized operator.
  const int order = 1;
  const double mu_r = 1.0, mu_e = 4.0;
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
          {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
      {"Boundaries", {}},
      {"Solver",
       {{"Order", order},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();

  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const Vector f = ClosedCurrent(pmesh, order);
  Vector u = ss.SolveSource(f);
  const double e_sub = Energy(ss, u);  // 1/2 u^T K_curlcurl u (magnetic energy)

  // Monolith reference with the identical approach: solve curl-curl + small mass (definite,
  // AMS-convergent), measure energy with the pure curl-curl operator.
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector nu_by_attr(max_attr), mass_by_attr(max_attr);
  nu_by_attr = 0.0;
  mass_by_attr = 0.0;
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  const double eps = kMagRegularization * RegularizationScale(pmesh);
  mass_by_attr(0) = eps / mu_r;
  mass_by_attr(1) = eps / mu_e;
  mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
  mfem::ParBilinearForm asolve(&pfes);
  asolve.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
  asolve.Assemble();
  asolve.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Ksolve(asolve.ParallelAssemble());
  mfem::ParBilinearForm aenergy(&pfes);
  aenergy.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  aenergy.Assemble();
  aenergy.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Kpure(aenergy.ParallelAssemble());
  mfem::HypreAMS ams(*Ksolve, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ksolve);
  pcg.SetTol(1e-13);
  pcg.SetMaxIter(2000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  Vector u_full(pfes.GetTrueVSize());
  u_full = 0.0;
  pcg.Mult(f, u_full);
  Vector t(pfes.GetTrueVSize());
  Kpure->Mult(u_full, t);
  double local = 0.0;
  for (int i = 0; i < pfes.GetTrueVSize(); i++)
  {
    local += u_full(i) * t(i);
  }
  double e_mono = 0.0;
  MPI_Allreduce(&local, &e_mono, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
  e_mono *= 0.5;

  CHECK(e_sub > 0.0);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-10 * std::abs(e_mono));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver offline/online model reuse",
                 "[substructure][Serial][Parallel][GPU]")
{
  const std::string model_path = (temp_dir / "substruct_model_roundtrip.bin").string();
  auto make_config = [](const std::string &mode, const std::string &path)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", path}}}}}};
    return IoData(config, false);
  };

  // Offline: condense the environment and write the model to disk.
  IoData iodata_off = make_config("Offline", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  Vector u_off = off.SolveExcitation(1);

  // Online: load the saved model (no environment materialization) and solve again.
  IoData iodata_on = make_config("Online", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  Vector u_on = on.SolveExcitation(1);

  Vector d(u_on);
  d -= u_off;
  CHECK(d.Norml2() <= 1.0e-12 * (u_off.Norml2() + 1.0e-30));
}

#if defined(MFEM_USE_EXCEPTIONS)
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver saved model rejects a changed environment",
                 "[substructure][Serial][Parallel][GPU]")
{
  // The model records a fingerprint of the environment operator: loading it with a changed
  // environment (material or mesh) must fail instead of silently reusing a stale S_E, while
  // changes confined to the region are accepted.
  const std::string model_path = (temp_dir / "substruct_env_fingerprint.bin").string();
  auto make_config =
      [&model_path](const std::string &mode, double eps_region, double eps_env)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", eps_region}},
            {{"Attributes", {2}}, {"Permittivity", eps_env}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto condense = [](const IoData &iodata, int a, int b)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(a, b, 4)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  };
  condense(make_config("Offline", 1.0, 10.0), 3, 4);
  // Region-only changes (material, re-meshed region) are accepted.
  CHECK_NOTHROW(condense(make_config("Online", 2.5, 10.0), 5, 4));
  // A changed environment material or environment mesh is rejected.
  CHECK_THROWS_WITH(condense(make_config("Online", 1.0, 12.0), 3, 4),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
  CHECK_THROWS_WITH(condense(make_config("Online", 1.0, 10.0), 3, 5),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
  // Also with fewer environment DOFs than the model, like the mesh of the region.
  CHECK_THROWS_WITH(condense(make_config("Online", 1.0, 10.0), 3, 3),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
  // So is a changed boundary condition on the environment (its terminal removed).
  IoData no_env_terminal = make_config("Online", 1.0, 10.0);
  no_env_terminal.boundaries.terminal.erase(2);
  CHECK_THROWS_WITH(condense(no_env_terminal, 3, 4),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver saved model rejects mirrored environment changes",
                 "[substructure][Serial][Parallel][GPU]")
{
  // An environment of two halves mirrored about y = 0: swapping their materials, or moving
  // a grounded boundary to its mirror image (with as many Dirichlet DOFs), changes the
  // environment, which its fingerprint must see although its counts are unchanged.
  // Not structured bindings, which make_config cannot capture with OpenMP (Clang 18).
  std::string type;
  int order;
  std::tie(type, order) = GENERATE(table<std::string, int>({{"Electrostatic", 1},
                                                            {"Electrostatic", 2},
                                                            {"Magnetostatic", 1},
                                                            {"Magnetostatic", 2}}));
  CAPTURE(type, order);
  const bool electrostatic = (type == "Electrostatic");
  const std::string model_path = (temp_dir / "substruct_mirror.bin").string();
  auto make_config = [&](const std::string &mode, double m2, double m3, int grounded)
  {
    const char *property = electrostatic ? "Permittivity" : "Permeability";
    json boundaries = {{electrostatic ? "Ground" : "PEC", {{"Attributes", {grounded}}}}};
    if (electrostatic)
    {
      boundaries["Terminal"] = {{{"Index", 1}, {"Attributes", {1}}},
                                {{"Index", 2}, {"Attributes", {2}}}};
    }
    json config = {{"Problem", {{"Type", type}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {property, 1.0}},
                       {{"Attributes", {2}}, {property, m2}},
                       {{"Attributes", {3}}, {property, m3}}}}}},
                   {"Boundaries", boundaries},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2, 3}}}},
                       {"Mode", mode},
                       {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  // The unit cube moved to -0.5 < y < 0.5: the region x < 0.5 (attribute 1), the
  // environment halves y < 0 and y > 0 (2 and 3), the boundaries x = 0 and x = 1 (1 and 2),
  // the faces y = -0.5 and y = 0.5 with x > 0.75 (4 and 5, away from the interface), and
  // the others (3).
  auto condense = [](const IoData &iodata)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
    for (int v = 0; v < serial.GetNV(); v++)
    {
      serial.GetVertex(v)[1] -= 0.5;
    }
    mfem::Vector c;
    for (int e = 0; e < serial.GetNE(); e++)
    {
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : (c(1) < 0.0) ? 2 : 3);
    }
    for (int b = 0; b < serial.GetNBE(); b++)
    {
      mfem::Array<int> vtx;
      serial.GetBdrElementVertices(b, vtx);
      double xc[2] = {0.0, 0.0};
      for (int j = 0; j < vtx.Size(); j++)
      {
        xc[0] += serial.GetVertex(vtx[j])[0] / vtx.Size();
        xc[1] += serial.GetVertex(vtx[j])[1] / vtx.Size();
      }
      const int attr = (xc[0] < 1e-9)                          ? 1
                       : (xc[0] > 1.0 - 1e-9)                  ? 2
                       : (xc[0] > 0.75 && xc[1] < -0.5 + 1e-9) ? 4
                       : (xc[0] > 0.75 && xc[1] > 0.5 - 1e-9)  ? 5
                                                               : 3;
      serial.SetBdrAttribute(b, attr);
    }
    serial.SetAttributes();
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  };
  condense(make_config("Offline", 2.0, 5.0, 4));
  CHECK_NOTHROW(condense(make_config("Online", 2.0, 5.0, 4)));
  CHECK_THROWS_WITH(condense(make_config("Online", 5.0, 2.0, 4)),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
  CHECK_THROWS_WITH(condense(make_config("Online", 2.0, 5.0, 5)),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver saved model of another problem or order",
                 "[substructure][Serial][Parallel]")
{
  // A model is loaded by a run of its problem type and order only, and a file that is not
  // one is reported as such.
  const std::string model_path = (temp_dir / "substruct_problem.bin").string();
  auto make_config =
      [&model_path](const std::string &type, const std::string &mode, int order)
  {
    const json terminals = {{{"Index", 1}, {"Attributes", {1}}},
                            {{"Index", 2}, {"Attributes", {2}}}};
    json config = {{"Problem", {{"Type", type}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains", {{"Materials", {{{"Attributes", {1, 2}}}}}}},
                   {"Boundaries", (type == "Electrostatic")
                                      ? json{{"Terminal", terminals}}
                                      : json{{"PEC", {{"Attributes", {1, 2}}}}}},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto condense = [](const IoData &iodata)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(3, 4, 4)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  };
  condense(make_config("Electrostatic", "Offline", 1));
  CHECK_NOTHROW(condense(make_config("Electrostatic", "Online", 1)));
  CHECK_THROWS_WITH(
      condense(make_config("Magnetostatic", "Online", 1)),
      Catch::Matchers::ContainsSubstring("is electrostatic, unlike this run"));
  CHECK_THROWS_WITH(
      condense(make_config("Electrostatic", "Online", 2)),
      Catch::Matchers::ContainsSubstring("has \"Order\": 1, unlike this run (2)"));
  // The model is read by the root process: a model of another format version, and a
  // truncated one.
  if (Mpi::Root(Mpi::World()))
  {
    std::fstream f(model_path, std::ios::binary | std::ios::in | std::ios::out);
    const int version = 1000;
    f.seekp(sizeof(int));
    f.write(reinterpret_cast<const char *>(&version), sizeof(int));
  }
  CHECK_THROWS_WITH(condense(make_config("Electrostatic", "Online", 1)),
                    Catch::Matchers::ContainsSubstring("written by another version"));
  // Cut inside its last section (the terminal modes), and inside S_E.
  condense(make_config("Electrostatic", "Offline", 1));
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::resize_file(model_path, std::filesystem::file_size(model_path) - 8);
  }
  CHECK_THROWS_WITH(condense(make_config("Electrostatic", "Online", 1)),
                    Catch::Matchers::ContainsSubstring("is incomplete"));
  condense(make_config("Electrostatic", "Offline", 1));
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::resize_file(model_path, std::filesystem::file_size(model_path) / 2);
  }
  CHECK_THROWS_WITH(condense(make_config("Electrostatic", "Online", 1)),
                    Catch::Matchers::ContainsSubstring("is incomplete"));
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream(model_path, std::ios::binary) << "not a substructuring model";
  }
  CHECK_THROWS_WITH(
      condense(make_config("Electrostatic", "Online", 1)),
      Catch::Matchers::ContainsSubstring("is not a saved substructuring model"));
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::remove(model_path);
  }
  CHECK_THROWS_WITH(condense(make_config("Electrostatic", "Online", 1)),
                    Catch::Matchers::ContainsSubstring("Cannot read the saved"));
  // A model that cannot be written, found before the condensation.
  IoData unwritable = make_config("Electrostatic", "Offline", 1);
  unwritable.solver.substructuring->save_model =
      (temp_dir / "missing" / "model.bin").string();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(3, 4, 4)));
  SubstructuringSolver ss(unwritable, mesh);
  CHECK_THROWS_WITH(ss.CondenseEnvironment(),
                    Catch::Matchers::ContainsSubstring("Writing the substructuring model"));
  CHECK_FALSE(ss.EnvironmentFactored());
}
#endif

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver environment-free capacitance matrix",
                 "[substructure][Serial][Parallel]")
{
  // The electrostatic capacitance matrix from the condensed environment (S_E +
  // terminal-mode couplings g_k, Cmode), with no environment solve, must equal the
  // full-field energy C_ij = u_i^T K u_j. Terminal 1 lies in the region and terminal 2 in
  // the environment, so both mode kinds are exercised. The online run (model loaded,
  // including a re-meshed region) must reproduce it WITHOUT factoring the environment;
  // saved fields are recovered on demand.
  const int order = GENERATE(1, 2);
  const bool remesh = GENERATE(false, true);
  CAPTURE(order, remesh);
  const std::string model_path = (temp_dir / "substruct_capacitance_model.bin").string();
  auto make_config = [order](const std::string &mode, const std::string &path)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", path}}}}}};
    return IoData(config, false);
  };
  const std::vector<int> terms = {1, 2};
  auto full_field_C = [&terms](SubstructuringSolver &ss)
  {
    std::vector<Vector> u;
    for (int t : terms)
    {
      u.push_back(ss.SolveExcitation(t));
    }
    mfem::DenseMatrix C(2);
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        C(i, j) = ss.MutualEnergy(u[i], u[j]);
      }
    }
    return C;
  };
  auto rel_diff = [](const mfem::DenseMatrix &A, const mfem::DenseMatrix &B)
  {
    double d = 0.0, m = 0.0;
    for (int i = 0; i < A.Height(); i++)
    {
      for (int j = 0; j < A.Width(); j++)
      {
        d = std::max(d, std::abs(A(i, j) - B(i, j)));
        m = std::max(m, std::abs(B(i, j)));
      }
    }
    return d / m;
  };

  // Offline: condensed capacitance vs. the full-field energy in the same run.
  IoData iodata_off = make_config("Offline", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const mfem::DenseMatrix C_off = off.CapacitanceMatrix(terms);
  const mfem::DenseMatrix C_ref_off = full_field_C(off);
  CHECK(rel_diff(C_off, C_ref_off) <= 1.0e-9);
  CHECK(std::abs(C_off(0, 1) - C_off(1, 0)) <= 1.0e-9 * std::abs(C_off(0, 0)));
  CHECK(C_off(0, 1) < 0.0);  // Maxwell capacitance: negative mutual term

  // Online (same or re-meshed region): no environment factorization for the matrix.
  IoData iodata_on = make_config("Online", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(remesh ? 7 : 4, 5, 6)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix C_on = on.CapacitanceMatrix(terms);
  CHECK_FALSE(on.EnvironmentFactored());
  if (!remesh)
  {
    CHECK(rel_diff(C_on, C_off) <= 1.0e-9);
  }

  // Recover one saved field: builds the environment once and matches the full solve.
  std::vector<Vector> fields;
  const mfem::DenseMatrix C_on2 = on.CapacitanceMatrix(terms, &fields, 1);
  CHECK(on.EnvironmentFactored());
  REQUIRE(fields.size() == 1);
  CHECK(rel_diff(C_on2, C_on) <= 1.0e-12);
  const Vector u_full = on.SolveExcitation(1);
  Vector d(fields[0]);
  d -= u_full;
  const double un = std::sqrt(mfem::InnerProduct(Mpi::World(), u_full, u_full));
  CHECK(std::sqrt(mfem::InnerProduct(Mpi::World(), d, d)) <= 1.0e-8 * un);
  // The online condensed matrix equals the online full-field energy (re-meshed or not).
  CHECK(rel_diff(C_on, full_field_C(on)) <= 1.0e-9);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver surface-current energies across face orientations",
                 "[substructure][Serial][Parallel]")
{
  // As above, for a surface-current excitation: J_s = cos(pi y) e_x on the side walls of
  // the wavy tetrahedral split, across the interface, closing through the x = 0 and x = 1
  // PEC walls. Its environment part lies on boundary face DOFs, whose basis changes with
  // the face orientation; the online energy needs no environment solve and equals the
  // offline one.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_current_orient_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
        {"Boundaries",
         {{"PEC", {{"Attributes", {1, 2}}}},
          {"SurfaceCurrent", {{{"Index", 1}, {"Attributes", {3}}, {"Direction", "+X"}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto current = [order](Mesh &mesh)
  {
    auto &pmesh = mesh.Get();
    mfem::ND_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    mfem::Array<int> pec_marker(pmesh.bdr_attributes.Max()), pec_tdofs,
        wall_marker(pmesh.bdr_attributes.Max());
    pec_marker = 0;
    pec_marker[0] = pec_marker[1] = 1;
    pfes.GetEssentialTrueDofs(pec_marker, pec_tdofs);
    wall_marker = 0;
    wall_marker[2] = 1;
    mfem::VectorFunctionCoefficient c(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v = 0.0;
                                        v(0) = std::cos(M_PI * x(1));
                                      });
    mfem::ParLinearForm f(&pfes);
    f.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(c), wall_marker);
    f.Assemble();
    Vector J(pfes.GetTrueVSize());
    f.ParallelAssemble(J);
    for (int i : pec_tdofs)
    {
      J(i) = 0.0;
    }
    return std::vector<Vector>{J};
  };

  IoData iodata_off = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const mfem::DenseMatrix E_off = CurrentEnergies(off, {1}, current(*mesh_off.back()));

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4, true)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix E_on = CurrentEnergies(on, {1}, current(*mesh_on.back()));
  CHECK_FALSE(on.EnvironmentFactored());
  CAPTURE(E_on(0, 0), E_off(0, 0));
  CHECK(std::abs(E_on(0, 0) - E_off(0, 0)) <= 1.0e-9 * std::abs(E_off(0, 0)));
}

#if defined(MFEM_USE_MUMPS)
TEST_CASE("SubstructuringSolver block low-rank environment factorization",
          "[substructure][Serial][Parallel]")
{
  // Solver.Substructuring.FactorizationTol > 0 selects a block low-rank (BLR) MUMPS
  // environment factorization (the matrix is scaled so the tolerance is relative; the Schur
  // and the solves are rescaled). On fronts large enough to compress, it changes the
  // capacitance matrix by about the tolerance, and not more.
  auto capacitance = [](double tol)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 2},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"FactorizationTol", tol}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(16)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  const mfem::DenseMatrix C0 = capacitance(0.0), C1 = capacitance(1.0e-10);
  double d = 0.0, m = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      d = std::max(d, std::abs(C1(i, j) - C0(i, j)));
      m = std::max(m, std::abs(C0(i, j)));
    }
  }
  CAPTURE(d / m);
  CHECK(d > 1.0e-12 * m);  // compressed (rounding changes it by about 1e-14)
  CHECK(d <= 1.0e-6 * m);
}
#endif

TEST_CASE_METHOD(palace::test::SharedTempDir, "SubstructuringSolver column orderings",
                 "[substructure][Serial][Parallel]")
{
  // Solver.Linear.ColumnOrdering orders the MUMPS factorizations: the environment factor of
  // a saved magnetostatic model has no Schur variables, so ParMETIS analyzes it in
  // parallel. The energy of a closed current does not depend on the ordering.
  const std::string model_path = (temp_dir / "substruct_ordering_model.bin").string();
  auto energy = [&](const std::string &ordering)
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", 1.0}},
                       {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
                   {"Boundaries", {}},
                   {"Solver",
                    {{"Order", 1},
                     {"Linear", {{"ColumnOrdering", ordering}}},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", "Offline"},
                       {"SaveModel", model_path}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4)));
    const Vector J = ClosedCurrent(*mesh.back(), 1);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return CurrentEnergies(ss, {1}, {J})(0, 0);
  };
  const double E0 = energy("Default");
  for (const std::string ordering : {"ParMETIS", "AMD"})
  {
    const double E = energy(ordering);
    CAPTURE(ordering, E, E0);
    CHECK(std::abs(E - E0) <= 1.0e-10 * std::abs(E0));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver online London sheets added to or removed from the "
                 "region",
                 "[substructure][Serial][Parallel][GPU]")
{
  // The environment of a model was condensed with the regularization of the model, which
  // depends on whether there are London sheets: an online run that adds the first sheet of
  // the region, or removes its last, keeps it. The full mesh then accepts the model and
  // gives the energy of a current closed in the region that the mesh of the region gives.
  const bool offline_sheet = GENERATE(false, true);
  CAPTURE(offline_sheet);
  const std::string model_path = (temp_dir / "substruct_sheet_eps.bin").string();
  auto make_config = [&](const std::string &mode, bool sheet)
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", 1.0}},
                       {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
                   {"Boundaries", sheet ? json{{"Superconductor",
                                                {{{"Attributes", {1}},
                                                  {"PenetrationDepth", 0.2},
                                                  {"Thickness", 0.05}}}}}
                                        : json::object()},
                   {"Solver",
                    {{"Order", 1},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto energy = [](const IoData &iodata, std::unique_ptr<mfem::ParMesh> &&pmesh)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    const Vector J = ClosedCurrent(*mesh.back(), 1, 0.5);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return CurrentEnergies(ss, {1}, {J})(0, 0);
  };
  (void)energy(make_config("Offline", offline_sheet), MakeSplitCube(4));
  {
    // With sheets in the region only, the model holds no kinetic operator S^Kin.
    std::ifstream in(model_path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
    CHECK(bytes.find("SKIN") == std::string::npos);
  }
  mfem::Mesh serial(RegionMeshPath(model_path).c_str());
  const double E_full = energy(make_config("Online", !offline_sheet), MakeSplitCube(4));
  const double E_region = energy(make_config("Online", !offline_sheet),
                                 std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
  CAPTURE(E_full, E_region);
  CHECK(E_full > 0.0);
  CHECK(std::abs(E_region - E_full) <= 1.0e-9 * E_full);
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver on the mesh of the region needs the model's terminals",
    "[substructure][Serial][Parallel]")
{
  // On the mesh of the region, a terminal on the environment's boundary (attribute 2, at
  // x = 1) has no DOFs: the model must hold it, so a renamed one needs the full mesh (and
  // would otherwise give a zero row of C).
  const std::string model_path = (temp_dir / "substruct_terminals.bin").string();
  auto make_config = [&](const std::string &mode, int env_index)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1, 2}}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Ground", {{"Attributes", {3}}}},
          {"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}},
            {{"Index", env_index}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto capacitance = [](const IoData &iodata, std::unique_ptr<mfem::ParMesh> &&pmesh, int k)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, k});
  };
  (void)capacitance(make_config("Offline", 2), MakeSplitCube(4), 2);
  mfem::Mesh serial(RegionMeshPath(model_path).c_str());
  auto region = [&]() { return std::make_unique<mfem::ParMesh>(Mpi::World(), serial); };
  CHECK(capacitance(make_config("Online", 2), region(), 2)(1, 1) > 0.0);
  CHECK_THROWS_WITH(capacitance(make_config("Online", 5), region(), 5),
                    Catch::Matchers::ContainsSubstring("Excitation 5 has no DOFs"));
}

TEST_CASE("SubstructuringSolver magnetostatics ignores terminals",
          "[substructure][Serial][Parallel]")
{
  // As the regular curl-curl operator, magnetostatic substructuring does not ground
  // "Terminal" boundaries: one leaves the energy of a closed current unchanged.
  auto energy = [](bool terminal)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
        {"Boundaries", terminal ? json{{"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}}}
                                : json::object()},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
    const Vector J = ClosedCurrent(*mesh.back(), 1);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return CurrentEnergies(ss, {1}, {J})(0, 0);
  };
  const double E0 = energy(false), E1 = energy(true);
  CHECK(E0 > 0.0);
  CHECK(std::abs(E1 - E0) <= 1.0e-12 * E0);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver online terminal added in the region",
                 "[substructure][Serial][Parallel][GPU]")
{
  // A terminal added inside the region does not reach the environment, so its modes are
  // zero: online runs on the full mesh and on the mesh of the region solve no environment
  // and give the capacitance of a fresh condensation with it.
  const std::string model_path = (temp_dir / "substruct_new_terminal.bin").string();
  auto make_config =
      [](const std::string &mode, bool region_terminal, const std::string &path)
  {
    json terminals = {{{"Index", 2}, {"Attributes", {2}}}};
    if (region_terminal)
    {
      terminals.push_back({{"Index", 1}, {"Attributes", {1}}});
    }
    json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permittivity", 1.0}},
                       {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
                   {"Boundaries", {{"Terminal", terminals}}},
                   {"Solver",
                    {{"Order", 1},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", path}}}}}};
    return IoData(config, false);
  };
  auto capacitance =
      [](const IoData &iodata, std::unique_ptr<mfem::ParMesh> &&pmesh, bool *factored)
  {
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    mfem::DenseMatrix C;
    if (factored)
    {
      C = ss.CapacitanceMatrix({1, 2});
      *factored = ss.EnvironmentFactored();
    }
    return C;
  };
  (void)capacitance(make_config("Offline", false, model_path), MakeSplitCube(4), nullptr);
  bool factored = true;
  const mfem::DenseMatrix C_ref =
      capacitance(make_config("Offline", true, (temp_dir / "substruct_ref.bin").string()),
                  MakeSplitCube(4), &factored);
  const mfem::DenseMatrix C_full =
      capacitance(make_config("Online", true, model_path), MakeSplitCube(4), &factored);
  CHECK_FALSE(factored);
  mfem::Mesh serial(RegionMeshPath(model_path).c_str());
  const mfem::DenseMatrix C_region =
      capacitance(make_config("Online", true, model_path),
                  std::make_unique<mfem::ParMesh>(Mpi::World(), serial), &factored);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CAPTURE(i, j, C_ref(i, j), C_full(i, j), C_region(i, j));
      CHECK(std::abs(C_full(i, j) - C_ref(i, j)) <= 1.0e-10 * C_ref(i, i));
      CHECK(std::abs(C_region(i, j) - C_ref(i, j)) <= 1.0e-10 * C_ref(i, i));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver incomplete model on every process",
                 "[substructure][Serial][Parallel][GPU]")
{
  // A model cut after its reference length is incomplete: every process reads the same of
  // it (the length is not used), and the online run says that it is incomplete.
  const std::string model_path = (temp_dir / "substruct_incomplete_model.bin").string();
  auto make_config = [&model_path](const std::string &mode, double Lc)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    if (Lc > 0.0)
    {
      config["Model"]["Lc"] = Lc;
    }
    return IoData(config, false);
  };
  {
    IoData offline = make_config("Offline", 2.0);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
    SubstructuringSolver off(offline, mesh);
    off.CondenseEnvironment();
  }
  // Keep the records up to the reference length.
  const auto file = substructuring::ModelFile::Read(model_path, Mpi::World());
  const auto *length = file.Find(substructuring::ModelFile::Tag("LENC"));
  REQUIRE(length);
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::resize_file(model_path, length->pos + length->size);
  }
  Mpi::Barrier(Mpi::World());
  IoData online = make_config("Online", 0.0);
  SubstructuringSolver::UseModelLength(online, Mpi::World());
  double Lc[2] = {online.model.Lc, -online.model.Lc};
  Mpi::GlobalMin(2, Lc, Mpi::World());
  CHECK(Lc[0] == -Lc[1]);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
  CHECK_THROWS_WITH(SubstructuringSolver(online, mesh),
                    Catch::Matchers::ContainsSubstring("is incomplete"));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver writes no mesh of the region of a nonconforming mesh",
    "[substructure][Serial][Parallel][GPU]")
{
  // On a nonconforming mesh, whose refined environment elements next to the interface are
  // among those written with the mesh of the region (they touch the grounded faces that
  // reach the interface), no mesh of the region is written; the model still serves an
  // online run on the full mesh.
  const std::string model_path = (temp_dir / "substruct_nc_model.bin").string();
  auto make_config = [&model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}},
          {"Ground", {{"Attributes", {3}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto capacitance = [](const IoData &iodata)
  {
    auto pmesh = MakeGradedSplit(4, 4, 4, {}, true);
    mfem::Array<int> refs;
    mfem::Vector c;
    for (int e = 0; e < pmesh->GetNE(); e++)
    {
      pmesh->GetElementCenter(e, c);
      if (c(0) > 0.5 && c(0) < 0.75)
      {
        refs.Append(e);
      }
    }
    pmesh->GeneralRefinement(refs, 1);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    REQUIRE(mesh.back()->Get().Nonconforming());
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  const mfem::DenseMatrix C_off = capacitance(make_config("Offline"));
  CHECK_FALSE(std::filesystem::exists(RegionMeshPath(model_path)));
  const mfem::DenseMatrix C_on = capacitance(make_config("Online"));
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(std::abs(C_on(i, j) - C_off(i, j)) <= 1.0e-10 * C_off(i, i));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver capacitance from a model without terminal modes",
                 "[substructure][Serial][Parallel][GPU]")
{
  // Modes missing from a saved model are recomputed on demand (this needs the environment
  // once) and give the same matrix.
  const std::string model_path = (temp_dir / "substruct_nomodes_model.bin").string();
  auto make_config = [&model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  const std::vector<int> terms = {1, 2};
  IoData iodata_off = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const mfem::DenseMatrix C_off = off.CapacitanceMatrix(terms);

  // Strip the record of the terminal modes, the last one.
  const auto file = substructuring::ModelFile::Read(model_path, Mpi::World());
  REQUIRE(!file.records.empty());
  REQUIRE(file.records.back().tag == substructuring::ModelFile::Tag("MODE"));
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::resize_file(model_path, file.records.back().pos - 12);
  }
  Mpi::Barrier(Mpi::World());

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  CHECK_FALSE(on.EnvironmentFactored());
  const mfem::DenseMatrix C_on = on.CapacitanceMatrix(terms);
  CHECK(on.EnvironmentFactored());  // modes recomputed on demand
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(std::abs(C_on(i, j) - C_off(i, j)) <= 1.0e-9 * std::abs(C_off(0, 0)));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver cross-run region re-meshing",
                 "[substructure][Serial][Parallel][GPU]")
{
  // Offline: condense the environment on one region mesh and save S_E. Online: RE-MESH the
  // region (different x-resolution) while keeping the interface Gamma and the environment
  // fixed, load and geometrically re-order S_E onto the new interface DOFs, and solve. The
  // re-meshed region-condensed energy must match a monolith on the online mesh (S_E is
  // exact for the fixed environment).
  const std::string model_path = (temp_dir / "substruct_remesh_model.bin").string();
  // Several (offline, online) region resolutions: re-mesh finer, coarser, and identical
  // (the identity-reorder sanity case). All must match the monolith on the online mesh.
  const auto res = GENERATE(std::make_pair(4, 8), std::make_pair(8, 4),
                            std::make_pair(5, 9), std::make_pair(6, 6));
  const int a_off = res.first, a_on = res.second;
  const int order = GENERATE(1, 2);
  CAPTURE(a_off, a_on, order);
  auto make_config = [order](const std::string &mode, const std::string &path)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", path}}}}}};
    return IoData(config, false);
  };

  // Offline: coarse region (a=4), environment b=5, interface n=6.
  IoData iodata_off = make_config("Offline", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeGradedSplit(a_off, 5, 6)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  (void)off.SolveExcitation(1);

  // Online: re-meshed region (a=8), same environment (b=5) and interface (n=6).
  IoData iodata_on = make_config("Online", model_path);
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(a_on, 5, 6)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  Vector u_on = on.SolveExcitation(1);
  const double e_sub = Energy(on, u_on);

  // Monolith on the online mesh, same excitation (attr 1 -> 1 V, attr 2 -> 0).
  auto &pmesh = mesh_on.back()->Get();
  mfem::H1_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector eps_by_attr(max_attr);
  eps_by_attr = 1.0;
  eps_by_attr(1) = 10.0;  // environment (attr 2)
  mfem::PWConstCoefficient eps(eps_by_attr);
  mfem::ParBilinearForm k(&pfes);
  k.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
  k.Assemble();
  k.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> K(k.ParallelAssemble());
  mfem::ParBilinearForm a(&pfes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
  a.Assemble();
  mfem::ParGridFunction xgf(&pfes);
  xgf = 0.0;
  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> m1(maxb), m2(maxb);
  m1 = 0;
  m2 = 0;
  m1[0] = 1;
  m2[1] = 1;
  mfem::ConstantCoefficient one(1.0), zero(0.0);
  xgf.ProjectBdrCoefficient(one, m1);
  xgf.ProjectBdrCoefficient(zero, m2);
  mfem::Array<int> ess_bdr(maxb), ess_tdofs;
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  mfem::ParLinearForm bform(&pfes);
  bform = 0.0;
  bform.Assemble();
  mfem::OperatorPtr A;
  mfem::Vector B, X;
  a.FormLinearSystem(ess_tdofs, xgf, bform, A, X, B);
  mfem::HypreParMatrix *Ah = A.As<mfem::HypreParMatrix>();
  mfem::HypreBoomerAMG amg(*Ah);
  amg.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ah);
  pcg.SetTol(1e-12);
  pcg.SetMaxIter(500);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(amg);
  pcg.Mult(B, X);
  a.RecoverFEMSolution(X, bform, xgf);
  mfem::Vector u_full, t(pfes.GetTrueVSize());
  xgf.GetTrueDofs(u_full);
  K->Mult(u_full, t);
  double local = 0.0;
  for (int i = 0; i < pfes.GetTrueVSize(); i++)
  {
    local += u_full(i) * t(i);
  }
  double e_mono = 0.0;
  MPI_Allreduce(&local, &e_mono, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
  e_mono *= 0.5;

  CHECK(e_sub > 1.0e-8);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-10 * std::abs(e_mono));
}

namespace
{

// Monolithic reference capacitance C_ij = u_i^T K u_j (terminal i at 1 V on boundary
// attribute i + 1, the others at 0 V), with MFEM directly (independent of substructuring).
// eps_by_attr: permittivity by element attribute.
mfem::DenseMatrix MonolithCapacitance(mfem::ParMesh &pmesh, int order,
                                      const mfem::Vector &eps_by_attr, int n_terminals,
                                      const std::vector<int> &grounded = {})
{
  mfem::H1_FECollection fec(order, pmesh.Dimension());
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::PWConstCoefficient eps(eps_by_attr);
  mfem::ParBilinearForm a(&pfes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> K(a.ParallelAssemble());
  const int maxb = pmesh.bdr_attributes.Max();
  std::vector<mfem::Vector> u(n_terminals);
  for (int i = 0; i < n_terminals; i++)
  {
    mfem::ParGridFunction x(&pfes);
    x = 0.0;
    mfem::Array<int> drive(maxb), ess_bdr(maxb), ess_tdofs;
    drive = 0;
    drive[i] = 1;
    ess_bdr = 0;
    for (int t = 0; t < n_terminals; t++)
    {
      ess_bdr[t] = 1;
    }
    for (int a : grounded)
    {
      ess_bdr[a - 1] = 1;
    }
    mfem::ConstantCoefficient one(1.0);
    x.ProjectBdrCoefficient(one, drive);
    pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    mfem::ParLinearForm rhs(&pfes);
    rhs = 0.0;
    mfem::OperatorPtr A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdofs, x, rhs, A, X, B);
    mfem::HypreBoomerAMG amg(*A.As<mfem::HypreParMatrix>());
    amg.SetPrintLevel(0);
    mfem::HyprePCG pcg(*A.As<mfem::HypreParMatrix>());
    pcg.SetTol(1.0e-14);
    pcg.SetMaxIter(2000);
    pcg.SetPrintLevel(0);
    pcg.SetPreconditioner(amg);
    pcg.Mult(B, X);
    a.RecoverFEMSolution(X, rhs, x);
    x.GetTrueDofs(u[i]);
  }
  mfem::DenseMatrix C(n_terminals);
  mfem::Vector Ku(pfes.GetTrueVSize());
  for (int j = 0; j < n_terminals; j++)
  {
    K->Mult(u[j], Ku);
    for (int i = 0; i < n_terminals; i++)
    {
      C(i, j) = mfem::InnerProduct(pmesh.GetComm(), u[i], Ku);
    }
  }
  return C;
}

// Nonconforming refinement of the graded split (as produced by AMR): refined elements on
// both sides of the interface x = 0.5 (so Gamma has hanging nodes from either side) and
// away from it.
std::unique_ptr<mfem::ParMesh> MakeRefinedGradedSplit(int a, int b, int n,
                                                      const RegionDesign &design = {})
{
  auto pmesh = MakeGradedSplit(a, b, n, design, true);
  mfem::Array<int> marked;
  mfem::Vector c;
  for (int e = 0; e < pmesh->GetNE(); e++)
  {
    pmesh->GetElementCenter(e, c);
    if ((c(0) > 0.37 && c(0) < 0.5 && c(1) < 0.5) ||  // region side of Gamma
        (c(0) > 0.5 && c(0) < 0.62 && c(1) > 0.5) ||  // environment side of Gamma
        (c(0) < 0.2 && c(2) < 0.3) || (c(0) > 0.85 && c(2) > 0.7))
    {
      marked.Append(e);
    }
  }
  pmesh->GeneralRefinement(marked, 1);
  return pmesh;
}

}  // namespace

TEST_CASE("SubstructuringSolver region solve by conjugate gradients",
          "[substructure][Serial][Parallel]")
{
  // With one excitation and an interface of over 1000 DOFs, the region is solved by CG,
  // preconditioned by the factor of the region operator without S_E, also for a region
  // without a Dirichlet boundary, whose region operator is singular (its factor still
  // preconditions the nonsingular condensed operator). The capacitance is that of a
  // monolith, of the terminal at x = 1 against the grounded sides of the environment, with
  // the region grounded at x = 0 or not.
  const bool grounded = GENERATE(true, false);
  CAPTURE(grounded);
  const std::vector<int> ground = grounded ? std::vector<int>{2, 4} : std::vector<int>{4};
  json boundaries = {{"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
                     {"Ground", {{"Attributes", ground}}}};
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 1.0}},
          {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
      {"Boundaries", boundaries},
      {"Solver",
       {{"Order", 2},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);
  // The split cube with the boundary attributes 1 at x = 1, 2 at x = 0, 4 on the sides of
  // the environment away from the interface (x > 9 / 16) and 3 elsewhere.
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(16, 16, 16, mfem::Element::HEXAHEDRON);
  mfem::Vector c;
  for (int e = 0; e < serial.GetNE(); e++)
  {
    serial.GetElementCenter(e, c);
    serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
  }
  for (int b = 0; b < serial.GetNBE(); b++)
  {
    mfem::Array<int> vtx;
    serial.GetBdrElementVertices(b, vtx);
    double xc = 0.0;
    for (int j = 0; j < vtx.Size(); j++)
    {
      xc += serial.GetVertex(vtx[j])[0] / vtx.Size();
    }
    serial.SetBdrAttribute(b, (xc > 1.0 - 1e-9) ? 1
                              : (xc < 1e-9)     ? 2
                              : (xc > 0.5625)   ? 4
                                                : 3);
  }
  serial.SetAttributes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  const double C_sub = ss.CapacitanceMatrix({1})(0, 0);
  mfem::Vector eps(2);
  eps(0) = 1.0;
  eps(1) = 10.0;
  const double C_ref = MonolithCapacitance(mesh.back()->Get(), 2, eps, 1, ground)(0, 0);
  CAPTURE(C_sub, C_ref);
  CHECK(std::abs(C_sub - C_ref) <= 1.0e-8 * C_ref);
}

TEST_CASE("SubstructuringSolver grounded PEC boundaries",
          "[substructure][Serial][Parallel]")
{
  // PEC boundaries are grounded (0 V) Dirichlet boundaries, as in the native electrostatic
  // operator: here on the environment's far face and on the side faces, which cross the
  // interface. The capacitance of the one terminal must match a monolith with those
  // boundaries grounded.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 1.0}},
          {{"Attributes", {3}}, {"Permittivity", 4.0}},
          {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
      {"Boundaries",
       {{"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
        {"PEC", {{"Attributes", {2, 3}}}}}},
      {"Solver",
       {{"Order", order},
        {"Substructuring",
         {{"Region", {{"Attributes", {1, 3}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);
  const RegionDesign design{0.5, {0.25, 0.5, 0.5}};
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  const mfem::DenseMatrix C = ss.CapacitanceMatrix({1});
  mfem::Vector eps_by_attr(3);
  eps_by_attr(0) = 1.0;
  eps_by_attr(1) = 10.0;
  eps_by_attr(2) = 4.0;
  const mfem::DenseMatrix C_mono =
      MonolithCapacitance(mesh.back()->Get(), order, eps_by_attr, 1, {2, 3});
  CAPTURE(C(0, 0), C_mono(0, 0));
  CHECK(std::abs(C(0, 0) - C_mono(0, 0)) <= 1.0e-9 * std::abs(C_mono(0, 0)));
}

TEST_CASE("SubstructuringSolver sparse direct solver choice",
          "[substructure][Serial][Parallel]")
{
  // Solver.Linear.Type names the sparse direct solver of the factorizations, else Palace
  // prefers STRUMPACK, MUMPS, then SuperLU_DIST. The environment is condensed by the
  // partial factorization of MUMPS, or else of STRUMPACK with its Schur complement
  // interface, and otherwise by solves. The capacitance matrix does not depend on it.
  auto config = [](const std::string &type)
  {
    return json{{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                {"Model", {{"Mesh", "test.msh"}}},
                {"Domains",
                 {{"Materials",
                   {{{"Attributes", {1}}, {"Permittivity", 1.0}},
                    {{"Attributes", {3}}, {"Permittivity", 4.0}},
                    {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
                {"Boundaries",
                 {{"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
                  {"PEC", {{"Attributes", {2, 3}}}}}},
                {"Solver",
                 {{"Order", 2},
                  {"Linear", {{"Type", type}}},
                  {"Substructuring",
                   {{"Region", {{"Attributes", {1, 3}}}},
                    {"Environment", {{"Attributes", {2}}}}}}}}};
  };
  std::vector<std::string> types;
#if defined(MFEM_USE_STRUMPACK)
  types.push_back("STRUMPACK");
#endif
#if defined(MFEM_USE_MUMPS)
  types.push_back("MUMPS");
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back("SuperLU");
#endif

  // The solver of each choice, and its partial factorization of a small SPD matrix.
  mfem::Mesh serial = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL);
  mfem::ParMesh pmesh(Mpi::World(), serial);
  mfem::H1_FECollection fec(1, 2);
  mfem::ParFiniteElementSpace fes(&pmesh, &fec);
  mfem::ConstantCoefficient one(1.0);
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.AddDomainIntegrator(new mfem::MassIntegrator(one));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  const std::vector<HYPRE_BigInt> schur_vars = {0, 3};
  std::vector<std::string> choices = {"Default"};
  choices.insert(choices.end(), types.begin(), types.end());
  for (const auto &type : choices)
  {
    CAPTURE(type);
    IoData iodata(config(type), false);
    const std::string direct = substructuring::DirectSolverName(iodata);
    std::string partial;
    const bool built = (substructuring::MakePartialFactorization(iodata, *A, schur_vars,
                                                                 partial) != nullptr);
    std::string direct_ref = (type == "SuperLU") ? "SuperLU_DIST" : type, partial_ref;
#if defined(PALACE_WITH_STRUMPACK_SCHUR)
    partial_ref = (type == "Default" || type == "STRUMPACK") ? "STRUMPACK" : "";
#endif
#if defined(MFEM_USE_MUMPS)
    partial_ref = (type == "Default" || type == "MUMPS") ? "MUMPS" : partial_ref;
#endif
    if (type == "Default")
    {
      direct_ref = (types[0] == "SuperLU") ? "SuperLU_DIST" : types[0];
    }
    CHECK(direct == direct_ref);
    const auto solver = substructuring::MakeDirectSolver(iodata, Mpi::World());
    std::string made;
#if defined(MFEM_USE_STRUMPACK)
    made = dynamic_cast<StrumpackSolver *>(solver.get()) ? "STRUMPACK" : made;
#endif
#if defined(MFEM_USE_MUMPS)
    made = dynamic_cast<substructuring::MumpsDirectSolver *>(solver.get()) ? "MUMPS" : made;
#endif
#if defined(MFEM_USE_SUPERLU)
    made = dynamic_cast<SuperLUSolver *>(solver.get()) ? "SuperLU_DIST" : made;
#endif
    CHECK(made == direct_ref);
    CHECK(built == !partial_ref.empty());
    CHECK(partial == partial_ref);
  }

  // The capacitance of each choice, against a monolith: the grounded boundaries of the
  // environment cross the interface, so its Dirichlet rows couple to it.
  const RegionDesign design{0.5, {0.25, 0.5, 0.5}};
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  mfem::Vector eps_by_attr(3);
  eps_by_attr(0) = 1.0;
  eps_by_attr(1) = 10.0;
  eps_by_attr(2) = 4.0;
  const mfem::DenseMatrix C_mono =
      MonolithCapacitance(mesh.back()->Get(), 2, eps_by_attr, 1, {2, 3});
  for (const auto &type : choices)
  {
    IoData iodata(config(type), false);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    const double C = ss.CapacitanceMatrix({1})(0, 0);
    CAPTURE(type, C, C_mono(0, 0));
    CHECK(std::abs(C - C_mono(0, 0)) <= 1.0e-9 * std::abs(C_mono(0, 0)));
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "SubstructuringSolver London sheet energies",
                 "[substructure][Serial][Parallel]")
{
  // Magnetostatics with a London superconductor sheet (finite penetration depth) on an
  // interior plane crossing the interface: region and environment sheet faces, meeting
  // Gamma along their edges. Two sheet generators a_k drive the sources M_sheet a_k with
  // the sheet DOFs free; the energy matrix E_ij = u_i^T K_cc u_j + (u_i - a_i)^T M_sheet
  // (u_j - a_j) must match a monolith (curl-curl + sheet, with a refined regularization),
  // offline and after reloading the model, and so must the energies in the volume, without
  // the kinetic energy (which normalize the error indicators), also without a model.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const double mu_r = 1.0, mu_e = 4.0, lambda = 0.2, thickness = 0.05;
  const std::string model_path = (temp_dir / "substruct_london_model.bin").string();
  auto make_config =
      [&](const std::string &mode, bool save = true, const std::string &path = "")
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", mu_r}},
                       {{"Attributes", {2}}, {"Permeability", mu_e}}}}}},
                   {"Boundaries",
                    {{"PEC", {{"Attributes", {1, 2}}}},
                     {"Superconductor",
                      {{{"Attributes", {4}},
                        {"PenetrationDepth", lambda},
                        {"Thickness", thickness}}}}}},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", path.empty() ? model_path : path}}}}}};
    if (!save)
    {
      config["Solver"]["Substructuring"].erase("SaveModel");
    }
    return IoData(config, false);
  };
  RegionDesign design;
  design.sheet = true;
  IoData iodata = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  auto &pmesh = mesh.back()->Get();

  // Sheet generators with curl on the sheet (an in-plane rotation and a shear, like
  // flux-loop generators; a gradient would be a zero-energy state), interpolated on the
  // sheet DOFs.
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> sheet_marker(maxb), sheet_tdofs;
  sheet_marker = 0;
  sheet_marker[3] = 1;
  pfes.GetEssentialTrueDofs(sheet_marker, sheet_tdofs);
  std::vector<Vector> a(2, Vector(pfes.GetTrueVSize()));
  for (int k = 0; k < 2; k++)
  {
    mfem::VectorFunctionCoefficient c(3,
                                      [k](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v = 0.0;
                                        if (k == 0)
                                        {
                                          v(0) = -(x(2) - 0.5);
                                          v(2) = x(0) - 0.5;
                                        }
                                        else
                                        {
                                          v(2) = x(0);
                                        }
                                      });
    mfem::ParGridFunction g(&pfes);
    g.ProjectCoefficient(c);
    Vector full;
    g.GetTrueDofs(full);
    a[k] = 0.0;
    for (int i : sheet_tdofs)
    {
      a[k](i) = full(i);
    }
  }

  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  std::vector<double> V;
  const mfem::DenseMatrix E = SheetEnergies(ss, {1, 2}, a, nullptr, 0, &V);
  // The full-field path (fields requested) gives the same energies.
  std::vector<Vector> fields;
  const mfem::DenseMatrix E_full = SheetEnergies(ss, {1, 2}, a, &fields, 2);

  // Monolith: A = K_cc + M_sheet + eps M, u = A^-1 M_sheet a (PEC pinned) with one
  // refinement step for the regularization, and the energy formula above.
  const double L_ksq = SuperconductorSheetOperator::KineticSheetInductance(
      iodata.boundaries.superconductor[0].lambda_L,
      iodata.boundaries.superconductor[0].thickness);
  mfem::Vector nu_by_attr(2), eps_by_attr(2), sheet_by_attr(maxb);
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  eps_by_attr = 1.0e-3;
  sheet_by_attr = 0.0;
  sheet_by_attr(3) = 1.0 / L_ksq;
  mfem::PWConstCoefficient nu(nu_by_attr), epsc(eps_by_attr), sheetc(sheet_by_attr);
  auto assemble = [&](bool curl, bool mass, bool sheet)
  {
    mfem::ParBilinearForm f(&pfes);
    if (curl)
    {
      f.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
    }
    if (mass)
    {
      f.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsc));
    }
    if (sheet)
    {
      f.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(sheetc));
    }
    f.Assemble();
    f.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(f.ParallelAssemble());
  };
  auto A = assemble(true, true, true), Kcc = assemble(true, false, false),
       Ms = assemble(false, false, true), Dm = assemble(false, true, false);
  mfem::Array<int> pec_marker(maxb), pec_tdofs;
  pec_marker = 0;
  pec_marker[0] = pec_marker[1] = 1;
  pfes.GetEssentialTrueDofs(pec_marker, pec_tdofs);
  std::unique_ptr<mfem::HypreParMatrix> Ae(A->EliminateRowsCols(pec_tdofs));
  mfem::HypreAMS ams(*A, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*A);
  pcg.SetTol(1.0e-14);
  pcg.SetMaxIter(5000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  auto solve = [&](Vector rhs)
  {
    for (int i : pec_tdofs)
    {
      rhs(i) = 0.0;
    }
    Vector x(rhs.Size());
    x = 0.0;
    pcg.Mult(rhs, x);
    return x;
  };
  std::vector<Vector> u(2), d(2), Md(2);
  for (int k = 0; k < 2; k++)
  {
    Vector b(a[k].Size());
    Ms->Mult(a[k], b);
    u[k] = solve(b);
    Vector r(b.Size());
    Dm->Mult(u[k], r);
    u[k] += solve(r);
    d[k] = u[k];
    d[k] -= a[k];
    Md[k].SetSize(b.Size());
    Ms->Mult(d[k], Md[k]);
  }
  mfem::DenseMatrix E_mono(2);
  Vector t(pfes.GetTrueVSize());
  for (int j = 0; j < 2; j++)
  {
    Kcc->Mult(u[j], t);
    for (int i = 0; i < 2; i++)
    {
      E_mono(i, j) = mfem::InnerProduct(pmesh.GetComm(), u[i], t) +
                     mfem::InnerProduct(pmesh.GetComm(), d[i], Md[j]);
    }
  }
  double diff = 0.0, m = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      diff = std::max(diff, std::abs(E(i, j) - E_mono(i, j)));
      m = std::max(m, std::abs(E_mono(i, j)));
    }
  }
  CAPTURE(E(0, 0), E_mono(0, 0), E(0, 1), E_mono(0, 1), E(1, 1), E_mono(1, 1));
  CHECK(std::min(E_mono(0, 0), E_mono(1, 1)) >= 1.0e-3 * m);  // no zero-energy state
  double dfull = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      dfull = std::max(dfull, std::abs(E_full(i, j) - E_mono(i, j)));
    }
  }
  CAPTURE(E_full(0, 0), E_full(1, 1));
  CHECK(dfull <= 1.0e-9 * m);
  CHECK(diff <= 1.0e-9 * m);

  // The energies in the volume, (u_k^T K_cc u_k) / 2, with a kinetic share of the E_kk
  // large enough to tell them apart. Without a model (full solves with the environment),
  // those of the monolith up to the regularization of their fields (to first order, unlike
  // the stationary E_kk); with the model, those without.
  std::vector<double> V_mono(2), V_none;
  for (int k = 0; k < 2; k++)
  {
    Kcc->Mult(u[k], t);
    V_mono[k] = 0.5 * mfem::InnerProduct(pmesh.GetComm(), u[k], t);
    CHECK(E_mono(k, k) - 2.0 * V_mono[k] >= 0.1 * E_mono(k, k));
  }
  {
    IoData iodata_none = make_config("Offline", false);
    std::vector<std::unique_ptr<Mesh>> mesh_none;
    mesh_none.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
    SubstructuringSolver none(iodata_none, mesh_none);
    none.CondenseEnvironment();
    SheetEnergies(none, {1, 2}, a, nullptr, 0, &V_none);
  }
  REQUIRE(V_none.size() == 2);
  for (int k = 0; k < 2; k++)
  {
    CAPTURE(k, V_none[k], V_mono[k]);
    CHECK(std::abs(V_none[k] - V_mono[k]) <= 1.0e-7 * V_mono[k]);
  }
  auto check_volume = [&](const std::vector<double> &Vk)
  {
    REQUIRE(Vk.size() == 2);
    for (int k = 0; k < 2; k++)
    {
      CAPTURE(k, Vk[k], V_none[k]);
      CHECK(std::abs(Vk[k] - V_none[k]) <= 1.0e-9 * m);
    }
  };
  check_volume(V);

  // Online reuse of the saved model.
  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  std::vector<double> V_on;
  const mfem::DenseMatrix E_on = SheetEnergies(on, {1, 2}, a, nullptr, 0, &V_on);
  CHECK_FALSE(on.EnvironmentFactored());
  check_volume(V_on);
  double don = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      don = std::max(don, std::abs(E_on(i, j) - E(i, j)));
    }
  }
  CHECK(don <= 1.0e-9 * m);

  // Online fields: the region solution with the environment interior recovered from its
  // sources, as in the full-field path. They differ only in the zero-energy gauge (the
  // regularization fixes it to rounding), so they are compared in the energy norm.
  std::vector<Vector> fields_on;
  SheetEnergies(on, {1, 2}, a, &fields_on, 2);
  REQUIRE(fields_on.size() == 2);
  auto energy_norm = [&](const Vector &v)
  {
    Vector t1(v.Size()), t2(v.Size());
    Kcc->Mult(v, t1);
    Ms->Mult(v, t2);
    t1 += t2;
    return std::sqrt(std::max(0.0, mfem::InnerProduct(pmesh.GetComm(), v, t1)));
  };
  for (int k = 0; k < 2; k++)
  {
    Vector d(fields_on[k]);
    d -= fields[k];
    const double nd = energy_norm(d), nf = energy_norm(fields[k]);
    CAPTURE(k, nd, nf);
    CHECK(nd <= 1.0e-10 * nf);
  }

  // A model cut in the kinetic parts of its source modes (its last section) is incomplete.
  const std::string cut_path = (temp_dir / "substruct_london_cut.bin").string();
  if (Mpi::Root(Mpi::World()))
  {
    std::ifstream in(model_path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
    CHECK(bytes.find("SKIN") != std::string::npos);  // the kinetic operator S^Kin
    std::ofstream(cut_path, std::ios::binary)
        .write(bytes.data(), static_cast<std::streamsize>(bytes.size() - sizeof(double)));
  }
  Mpi::Barrier(Mpi::World());
  IoData iodata_cut = make_config("Online", true, cut_path);
  std::vector<std::unique_ptr<Mesh>> mesh_cut;
  mesh_cut.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  CHECK_THROWS_WITH(SubstructuringSolver(iodata_cut, mesh_cut),
                    Catch::Matchers::ContainsSubstring("is incomplete"));
}

namespace
{

// Magnetostatics in a box with PEC walls, driven by surface currents on an interior plane
// crossing the interface (from the x = 0 wall to the x = 1 wall, returning through the
// walls): region and environment parts of the excitation. The energies must match a
// monolith, offline and after reloading the model (also with an added port in the
// region), with and without a London sheet on the plane; a current that does not close is
// rejected.
void CheckSurfaceCurrentEnergies(const std::filesystem::path &temp_dir, int order,
                                 bool london)
{
  CAPTURE(order, london);
  const double mu_r = 1.0, mu_e = 4.0, lambda = 0.2, thickness = 0.05;
  const std::string model_path = (temp_dir / "substruct_current_model.bin").string();
  auto make_config = [&](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", mu_r}},
            {{"Attributes", {2}}, {"Permeability", mu_e}}}}}},
        {"Boundaries",
         {{"PEC", {{"Attributes", {1, 2, 3}}}},
          {"SurfaceCurrent", {{{"Index", 1}, {"Attributes", {4}}, {"Direction", "+X"}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    if (london)
    {
      config["Boundaries"].erase("SurfaceCurrent");  // a port may not be a Superconductor
      config["Boundaries"]["Superconductor"] = {
          {{"Attributes", {4}}, {"PenetrationDepth", lambda}, {"Thickness", thickness}}};
    }
    return IoData(config, false);
  };
  RegionDesign design;
  design.sheet = true;
  IoData iodata = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  auto &pmesh = mesh.back()->Get();

  // Excitations (the configured port is only validated), zero on the walls: sheet currents
  // J_s = x and J_s = x z^2 on the plane (both sides), J_s = z on its strip x < 3/8 (region
  // only, from the z = 0 to the z = 1 wall), and J_s = x on that strip, which ends on the
  // free sheet and closes only through a London sheet.
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> pec_marker(maxb), pec_tdofs, sheet_marker(maxb);
  pec_marker = 0;
  pec_marker[0] = pec_marker[1] = pec_marker[2] = 1;
  pfes.GetEssentialTrueDofs(pec_marker, pec_tdofs);
  sheet_marker = 0;
  sheet_marker[3] = 1;
  std::vector<Vector> J(4, Vector(pfes.GetTrueVSize()));
  for (int k = 0; k < 4; k++)
  {
    mfem::VectorFunctionCoefficient c(3,
                                      [k](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v = 0.0;
                                        const bool strip = x(0) < 0.375;
                                        if (k < 2)
                                        {
                                          v(0) = (k == 0) ? 1.0 : x(2) * x(2);
                                        }
                                        else if (strip)
                                        {
                                          v((k == 2) ? 2 : 0) = 1.0;
                                        }
                                      });
    mfem::ParLinearForm f(&pfes);
    f.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(c), sheet_marker);
    f.Assemble();
    f.ParallelAssemble(J[k]);
    for (int i : pec_tdofs)
    {
      J[k](i) = 0.0;
    }
  }
  const std::vector<Vector> J12 = {J[0], J[1]}, J123 = {J[0], J[1], J[2]};

  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  const mfem::DenseMatrix E = CurrentEnergies(ss, {1, 2}, J12);
  std::vector<Vector> fields;
  const mfem::DenseMatrix E_full = CurrentEnergies(ss, {1, 2}, J12, &fields, 2);

  // Monolith: L = K + M_sheet, u = L^+ J from A = L + eps M with two refinement steps
  // (PEC pinned), and E_ij = u_i^T L u_j.
  mfem::Vector nu_by_attr(2), eps_by_attr(2), sheet_by_attr(maxb);
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  eps_by_attr = 1.0e-3;
  sheet_by_attr = 0.0;
  if (london)
  {
    sheet_by_attr(3) = 1.0 / SuperconductorSheetOperator::KineticSheetInductance(
                                 iodata.boundaries.superconductor[0].lambda_L,
                                 iodata.boundaries.superconductor[0].thickness);
  }
  mfem::PWConstCoefficient nu(nu_by_attr), epsc(eps_by_attr), sheetc(sheet_by_attr);
  auto assemble = [&](bool mass)
  {
    mfem::ParBilinearForm f(&pfes);
    f.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
    if (mass)
    {
      f.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsc));
    }
    f.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(sheetc));
    f.Assemble();
    f.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(f.ParallelAssemble());
  };
  auto A = assemble(true), L = assemble(false);
  std::unique_ptr<mfem::HypreParMatrix> Ae(A->EliminateRowsCols(pec_tdofs));
  mfem::HypreAMS ams(*A, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*A);
  pcg.SetTol(1.0e-14);
  pcg.SetMaxIter(5000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  std::unique_ptr<mfem::HypreParMatrix> Le(L->EliminateRowsCols(pec_tdofs));
  std::vector<Vector> u(3);
  for (int k = 0; k < 3; k++)
  {
    u[k].SetSize(J[k].Size());
    u[k] = 0.0;
    Vector r(J[k]), du(J[k].Size()), Lu(J[k].Size());
    for (int it = 0; it < 3; it++)
    {
      du = 0.0;
      pcg.Mult(r, du);
      u[k] += du;
      L->Mult(u[k], Lu);
      r = J[k];
      r -= Lu;
    }
  }
  mfem::DenseMatrix E_mono(3);
  Vector t(pfes.GetTrueVSize());
  for (int j = 0; j < 3; j++)
  {
    L->Mult(u[j], t);
    for (int i = 0; i < 3; i++)
    {
      E_mono(i, j) = mfem::InnerProduct(pmesh.GetComm(), u[i], t);
    }
  }
  auto max_diff = [](const mfem::DenseMatrix &X, const mfem::DenseMatrix &Y)
  {
    double d = 0.0;
    for (int i = 0; i < X.Height(); i++)
    {
      for (int j = 0; j < X.Width(); j++)
      {
        d = std::max(d, std::abs(X(i, j) - Y(i, j)));
      }
    }
    return d;
  };
  mfem::DenseMatrix E_mono12(2);
  E_mono12.CopyMN(E_mono, 2, 2, 0, 0);
  // The regularization eps M of the substructured solve gives an energy error of order
  // (eps / lambda_1)^2, with lambda_1 the smallest nonzero eigenvalue of L relative to M:
  // about 1e-8 for eps = 1e-3 in this unit box, and below 1e-9 for eps = 1e-7 with sheets.
  const double m = E_mono.MaxMaxNorm(), tol = london ? 1.0e-9 : 1.0e-7;
  CAPTURE(E(0, 0), E_mono(0, 0), E(0, 1), E_mono(0, 1), E(1, 1), E_mono(1, 1));
  CHECK(std::min({E_mono(0, 0), E_mono(1, 1), E_mono(2, 2)}) >= 1.0e-3 * m);
  CHECK(std::abs(E_mono(0, 1)) >= 1.0e-3 * m);
  CHECK(max_diff(E, E_mono12) <= tol * m);
  CHECK(max_diff(E_full, E_mono12) <= tol * m);

  // Online reuse of the saved model, without the environment, also with an added port in
  // the region (no environment part, so no modes).
  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix E_on = CurrentEnergies(on, {1, 2}, J12);
  const mfem::DenseMatrix E_on3 = CurrentEnergies(on, {1, 2, 3}, J123);
  CHECK_FALSE(on.EnvironmentFactored());
  CHECK(max_diff(E_on, E) <= 1.0e-9 * m);
  CAPTURE(E_on3(2, 2), E_mono(2, 2), E_on3(0, 2), E_mono(0, 2));
  CHECK(max_diff(E_on3, E_mono) <= tol * m);

  // Online fields agree with the full-field path in the energy norm (they differ only in
  // the zero-energy gauge).
  std::vector<Vector> fields_on;
  CurrentEnergies(on, {1, 2}, J12, &fields_on, 2);
  REQUIRE(fields_on.size() == 2);
  for (int k = 0; k < 2; k++)
  {
    Vector d(fields_on[k]), Ld(d.Size()), Lf(d.Size());
    d -= fields[k];
    L->Mult(d, Ld);
    L->Mult(fields[k], Lf);
    const double nd = std::sqrt(std::max(0.0, mfem::InnerProduct(pmesh.GetComm(), d, Ld))),
                 nf = std::sqrt(mfem::InnerProduct(pmesh.GetComm(), fields[k], Lf));
    CAPTURE(k, nd, nf);
    CHECK(nd <= 1.0e-10 * nf);
  }

  // A reversed excitation does not match the saved source modes.
  std::vector<Vector> J_rev = {J[0]};
  J_rev[0].Neg();
  const mfem::DenseMatrix E_rev = CurrentEnergies(on, {1}, J_rev);
  CHECK(std::abs(E_rev(0, 0) - E(0, 0)) <= 1.0e-9 * m);

  // A current ending on the free sheet does not close (with PEC only).
  if (london)
  {
    CHECK_NOTHROW(CurrentEnergies(on, {4}, {J[3]}));
  }
#if defined(MFEM_USE_EXCEPTIONS)
  else
  {
    CHECK_THROWS_WITH(CurrentEnergies(on, {4}, {J[3]}),
                      Catch::Matchers::ContainsSubstring("does not close"));
  }
#endif

  // The surface flux functional matches ComputeFluxThroughSurface (the sheet crosses the
  // interface, so it has shared faces in parallel).
  if (!london)
  {
    CurlCurlOperator curlcurl_op(iodata, mesh);
    const auto &rt = curlcurl_op.GetRTSpace();
    mfem::Vector dir(3);
    dir = 0.0;
    dir(1) = 1.0;
    const Vector f = FluxThroughSurfaceFunctional(rt, {4}, dir);
    Vector B(rt.GetTrueVSize());
    B.Randomize(1 + Mpi::Rank(pmesh.GetComm()));
    mfem::ParGridFunction B_gf(&const_cast<mfem::ParFiniteElementSpace &>(rt.Get()));
    B_gf.SetFromTrueDofs(B);
    const double flux =
        ComputeFluxThroughSurface(B_gf, {4}, curlcurl_op.GetMesh(),
                                  curlcurl_op.GetMaterialOp(), dir, pmesh.GetComm());
    CHECK(std::abs(mfem::InnerProduct(pmesh.GetComm(), f, B) - flux) <=
          1.0e-12 * std::abs(flux));
  }

  // Mixed currents and London flux states, with functionals l_c (here the excitations J2,
  // with an environment part, and J3, without): the current block as above, the London
  // energy and l_c^T u_f against the monolith, offline and online without the environment.
  if (london)
  {
    mfem::Array<int> sheet_tdofs;
    pfes.GetEssentialTrueDofs(sheet_marker, sheet_tdofs);
    Vector a(pfes.GetTrueVSize()), full;
    {
      mfem::VectorFunctionCoefficient c(3,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v.SetSize(3);
                                          v = 0.0;
                                          v(0) = -(x(2) - 0.5);
                                          v(2) = x(0) - 0.5;
                                        });
      mfem::ParGridFunction g(&pfes);
      g.ProjectCoefficient(c);
      g.GetTrueDofs(full);
      a = 0.0;
      for (int i : sheet_tdofs)
      {
        a(i) = full(i);
      }
    }
    auto Ms = [&]()
    {
      mfem::ParBilinearForm m(&pfes);
      m.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(sheetc));
      m.Assemble();
      m.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(m.ParallelAssemble());
    }();
    Vector b(a.Size()), ua(a.Size()), r(a.Size()), du(a.Size()), Lu(a.Size());
    Ms->Mult(a, b);
    for (int i : pec_tdofs)
    {
      b(i) = 0.0;
    }
    ua = 0.0;
    r = b;
    for (int it = 0; it < 3; it++)
    {
      du = 0.0;
      pcg.Mult(r, du);
      ua += du;
      L->Mult(ua, Lu);
      r = b;
      r -= Lu;
    }
    // u^T K u + (u - a)^T M_sheet (u - a) = u^T L u - 2 a^T M_sheet u + a^T M_sheet a.
    Vector Ma(a.Size());
    Ms->Mult(a, Ma);
    L->Mult(ua, Lu);
    const double E_ff = mfem::InnerProduct(pmesh.GetComm(), ua, Lu) -
                        2.0 * mfem::InnerProduct(pmesh.GetComm(), Ma, ua) +
                        mfem::InnerProduct(pmesh.GetComm(), Ma, a);
    const std::vector<Vector> Jm = {J[1], J[2]};
    for (const bool online : {false, true})
    {
      IoData io = make_config(online ? "Online" : "Offline");
      std::vector<std::unique_ptr<Mesh>> m_mix;
      m_mix.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
      SubstructuringSolver mix(io, m_mix);
      mix.CondenseEnvironment();
      mfem::DenseMatrix linked;
      const mfem::DenseMatrix E_mix =
          mix.MagnetostaticEnergyMatrix({2, 3}, Jm, {5}, {a}, Jm, &linked);
      CAPTURE(online, E_mix(2, 2), E_ff, linked(0, 0), linked(1, 0));
      CHECK(mix.EnvironmentFactored() == !online);
      CHECK(std::abs(E_mix(0, 0) - E_mono(1, 1)) <= tol * m);
      CHECK(std::abs(E_mix(1, 1) - E_mono(2, 2)) <= tol * m);
      CHECK(std::abs(E_mix(0, 1) - E_mono(1, 2)) <= tol * m);
      CHECK(E_mix(0, 2) == 0.0);
      CHECK(std::abs(E_mix(2, 2) - E_ff) <= 1.0e-9 * std::abs(E_ff));
      for (int c = 0; c < 2; c++)
      {
        const double ref = mfem::InnerProduct(pmesh.GetComm(), Jm[c], ua);
        CHECK(std::abs(linked(c, 0) - ref) <= 1.0e-9 * std::sqrt(m * std::abs(E_ff)));
      }
    }
  }
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver surface-current energies",
                 "[substructure][Serial][Parallel]")
{
  CheckSurfaceCurrentEnergies(temp_dir, 1, GENERATE(false, true));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver second-order surface-current energies",
                 "[substructure][Serial][Parallel][Long]")
{
  CheckSurfaceCurrentEnergies(temp_dir, 2, GENERATE(false, true));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver on a nonconforming mesh",
                 "[substructure][Serial][Parallel]")
{
  // Offline condensation on an adapted (nonconforming) mesh must reproduce the monolith on
  // the same mesh, including hanging nodes on the interface from either side. A partial
  // conductor and a dielectric block in the region make the field vary across the hanging
  // faces (a field linear in x would satisfy the constraints trivially).
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_nc_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {3}}, {"Permittivity", 4.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1, 3}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  const RegionDesign design{0.5, {0.25, 0.5, 0.5}};
  IoData iodata = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeRefinedGradedSplit(4, 5, 6, design)));
  REQUIRE(mesh.back()->Get().Nonconforming());
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  const mfem::DenseMatrix C = ss.CapacitanceMatrix({1, 2});
  mfem::Vector eps_by_attr(3);
  eps_by_attr(0) = 1.0;
  eps_by_attr(1) = 10.0;
  eps_by_attr(2) = 4.0;
  const mfem::DenseMatrix C_mono =
      MonolithCapacitance(mesh.back()->Get(), order, eps_by_attr, 2);
  double d = 0.0, m = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      d = std::max(d, std::abs(C(i, j) - C_mono(i, j)));
      m = std::max(m, std::abs(C_mono(i, j)));
    }
  }
  CAPTURE(C(0, 0), C_mono(0, 0), C(0, 1), C_mono(0, 1));
  CHECK(d <= 1.0e-9 * m);

  // The refinement changes the capacitance, so matching the refined monolith is meaningful.
  const auto coarse = MakeGradedSplit(4, 5, 6, design);
  const mfem::DenseMatrix C_coarse = MonolithCapacitance(*coarse, order, eps_by_attr, 2);
  CHECK(std::abs(C_coarse(0, 0) - C_mono(0, 0)) >= 1.0e-4 * m);

  // Online reuse of the model on the same nonconforming mesh (interface matched through the
  // master DOFs of the hanging faces) is exact and needs no environment factorization.
  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeRefinedGradedSplit(4, 5, 6, design)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix C_on = on.CapacitanceMatrix({1, 2});
  CHECK_FALSE(on.EnvironmentFactored());
  double don = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      don = std::max(don, std::abs(C_on(i, j) - C(i, j)));
    }
  }
  CHECK(don <= 1.0e-12 * m);
}

namespace
{

// H(curl) on an adapted (nonconforming) mesh with hanging edges on the interface: the
// condensed magnetic energy of a closed current must match a monolith of the operators of
// the regular solve on the same mesh.
void CheckNonconformingMagnetostaticEnergy(const std::filesystem::path &temp_dir, int order)
{
  CAPTURE(order);
  const double mu_r = 1.0, mu_e = 4.0;
  const std::string model_path = (temp_dir / "substruct_nc_mag_model.bin").string();
  auto make_config = [order, mu_r, mu_e, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
        {"Boundaries", {}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  IoData iodata = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeRefinedGradedSplit(4, 5, 6)));
  REQUIRE(mesh.back()->Get().Nonconforming());
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();

  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const Vector f = ClosedCurrent(pmesh, order);
  const double e_sub = Energy(ss, ss.SolveSource(f));

  const double mag_eps = kMagRegularization;
  CurlCurlOperator curlcurl_op(iodata, mesh);
  const auto &nd_fespace = curlcurl_op.GetSolutionNDSpace();
  const auto &mat_op = curlcurl_op.GetMaterialOp();
  auto K = curlcurl_op.GetUnconstrainedStiffnessMatrix(true, false);
  MaterialPropertyCoefficient eps_func(mat_op.GetAttributeToMaterial(),
                                       mat_op.GetInvPermeability(),
                                       mag_eps * RegularizationScale(pmesh));
  BilinearForm m(nd_fespace);
  m.AddDomainIntegrator<VectorFEMassIntegrator>(eps_func);
  ParOperator M(m.Assemble(false), nd_fespace);
  const mfem::HypreParMatrix &Kpure = K->ParallelAssemble();
  std::unique_ptr<mfem::HypreParMatrix> Ksolve(
      mfem::Add(1.0, Kpure, 1.0, M.ParallelAssemble()));
  mfem::HypreAMS ams(*Ksolve, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ksolve);
  pcg.SetTol(1.0e-13);
  pcg.SetMaxIter(4000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  Vector u(pfes.GetTrueVSize()), t(pfes.GetTrueVSize());
  u = 0.0;
  pcg.Mult(f, u);
  Kpure.Mult(u, t);
  const double e_mono = 0.5 * mfem::InnerProduct(pmesh.GetComm(), u, t);
  CAPTURE(e_sub, e_mono);
  CHECK(e_sub > 0.0);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-10 * std::abs(e_mono));

  // The refinement changes the energy by far more than the mismatch above, so the match is
  // meaningful.
  {
    mfem::Vector nu_by_attr(2), mass_by_attr(2);
    auto coarse = MakeGradedSplit(4, 5, 6);
    nu_by_attr(0) = 1.0 / mu_r;
    nu_by_attr(1) = 1.0 / mu_e;
    mass_by_attr = nu_by_attr;
    mass_by_attr *= mag_eps * RegularizationScale(*coarse);
    mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
    mfem::ParFiniteElementSpace cfes(coarse.get(), &fec);
    const Vector cf = ClosedCurrent(*coarse, order);
    mfem::ParBilinearForm cs(&cfes), ce(&cfes);
    cs.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
    cs.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
    ce.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
    cs.Assemble();
    cs.Finalize();
    ce.Assemble();
    ce.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> CS(cs.ParallelAssemble()),
        CE(ce.ParallelAssemble());
    mfem::HypreAMS cams(*CS, &cfes);
    cams.SetPrintLevel(0);
    mfem::HyprePCG cpcg(*CS);
    cpcg.SetTol(1.0e-13);
    cpcg.SetMaxIter(4000);
    cpcg.SetPrintLevel(0);
    cpcg.SetPreconditioner(cams);
    Vector cu(cfes.GetTrueVSize()), ct(cfes.GetTrueVSize());
    cu = 0.0;
    cpcg.Mult(cf, cu);
    CE->Mult(cu, ct);
    const double e_coarse = 0.5 * mfem::InnerProduct(coarse->GetComm(), cu, ct);
    CHECK(std::abs(e_coarse - e_mono) >= 100.0 * std::abs(e_sub - e_mono));
  }

  // Online reuse on the same nonconforming mesh: the saved S_E is matched onto the
  // interface edge DOFs, which include the master edges of hanging faces.
  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeRefinedGradedSplit(4, 5, 6)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const double e_on = Energy(on, on.SolveSource(f));
  CAPTURE(e_on);
  CHECK(std::abs(e_on - e_sub) <= 1.0e-10 * std::abs(e_sub));
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver magnetostatic energy on a nonconforming mesh",
                 "[substructure][Serial][Parallel]")
{
  CheckNonconformingMagnetostaticEnergy(temp_dir, 1);
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver second-order magnetostatic energy on a nonconforming mesh",
    "[substructure][Serial][Parallel][Long]")
{
  CheckNonconformingMagnetostaticEnergy(temp_dir, 2);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver online region adaptation",
                 "[substructure][Serial][Parallel]")
{
  // One step of online adaptive refinement of the region only: the region error indicators
  // are zero on the environment; nonconforming refinement of the marked region elements
  // leaves the environment and the interface DOFs as condensed, so the saved model still
  // applies without factoring the environment, matches a monolith on the refined mesh, and
  // the estimated region error decreases.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_amr_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {3}}, {"Permittivity", 4.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1, 3}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  const RegionDesign design{0.5, {0.25, 0.5, 0.5}};
  const std::vector<int> terms = {1, 2};
  {
    IoData iodata = make_config("Offline");
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design, true)));
    SubstructuringSolver off(iodata, mesh);
    off.CondenseEnvironment();
    (void)off.CapacitanceMatrix(terms);
  }
  IoData iodata_on = make_config("Online");
  auto count_env = [](const mfem::ParMesh &m)
  {
    long c = 0;
    for (int e = 0; e < m.GetNE(); e++)
    {
      c += (m.GetAttribute(e) == 2);
    }
    MPI_Allreduce(MPI_IN_PLACE, &c, 1, MPI_LONG, MPI_SUM, m.GetComm());
    return c;
  };

  // Online solve and region indicators on the initial mesh.
  std::vector<std::unique_ptr<Mesh>> mesh0;
  mesh0.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design, true)));
  auto &pm0 = mesh0.back()->Get();
  SubstructuringSolver on0(iodata_on, mesh0);
  on0.CondenseEnvironment();
  std::vector<Vector> region_fields;
  const mfem::DenseMatrix C0 = on0.CapacitanceMatrix(terms, nullptr, 0, &region_fields);
  auto half_diagonal = [](const mfem::DenseMatrix &C)
  {
    std::vector<double> e(C.Height());
    for (int i = 0; i < C.Height(); i++)
    {
      e[i] = 0.5 * C(i, i);
    }
    return e;
  };
  const ErrorIndicator ind0 = on0.RegionErrorIndicator(region_fields, half_diagonal(C0));
  CHECK_FALSE(on0.EnvironmentFactored());
  REQUIRE(ind0.Local().Size() == pm0.GetNE());
  double env_max = 0.0, reg_max = 0.0;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    (pm0.GetAttribute(e) == 2 ? env_max : reg_max) =
        std::max(pm0.GetAttribute(e) == 2 ? env_max : reg_max, ind0.Local()(e));
  }
  MPI_Allreduce(MPI_IN_PLACE, &env_max, 1, MPI_DOUBLE, MPI_MAX, pm0.GetComm());
  MPI_Allreduce(MPI_IN_PLACE, &reg_max, 1, MPI_DOUBLE, MPI_MAX, pm0.GetComm());
  CHECK(env_max == 0.0);
  CHECK(reg_max > 0.0);

  // Refine the region elements with the largest indicators (nonconforming, no level
  // constraint), on an identical copy of the mesh.
  auto pm1 = MakeGradedSplit(4, 5, 6, design, true);
  mfem::Array<int> marked;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    if (ind0.Local()(e) >= 0.3 * reg_max)
    {
      marked.Append(e);
    }
  }
  const long n_env = count_env(*pm1);
  pm1->GeneralRefinement(marked, 1, 0);
  CHECK(count_env(*pm1) == n_env);  // the environment is not refined

  std::vector<std::unique_ptr<Mesh>> mesh1;
  mesh1.push_back(std::make_unique<Mesh>(std::move(pm1)));
  SubstructuringSolver on1(iodata_on, mesh1);
  on1.CondenseEnvironment();
  const mfem::DenseMatrix C1 = on1.CapacitanceMatrix(terms, nullptr, 0, &region_fields);
  const ErrorIndicator ind1 = on1.RegionErrorIndicator(region_fields, half_diagonal(C1));
  CHECK_FALSE(on1.EnvironmentFactored());
  mfem::Vector eps_by_attr(3);
  eps_by_attr(0) = 1.0;
  eps_by_attr(1) = 10.0;
  eps_by_attr(2) = 4.0;
  const mfem::DenseMatrix C_mono =
      MonolithCapacitance(mesh1.back()->Get(), order, eps_by_attr, 2);
  double d = 0.0, m = 0.0, change = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      d = std::max(d, std::abs(C1(i, j) - C_mono(i, j)));
      m = std::max(m, std::abs(C_mono(i, j)));
      change = std::max(change, std::abs(C1(i, j) - C0(i, j)));
    }
  }
  CAPTURE(C0(0, 0), C1(0, 0), C_mono(0, 0));
  CHECK(d <= 1.0e-9 * m);
  CHECK(change >= 1.0e3 * d);  // the refinement matters at the tolerance of the check
  const MPI_Comm comm = mesh1.back()->GetComm();
  CAPTURE(ind0.Norml2(pm0.GetComm()), ind1.Norml2(comm));
  CHECK(ind1.Norml2(comm) < ind0.Norml2(pm0.GetComm()));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver online magnetostatic region adaptation",
                 "[substructure][Serial][Parallel]")
{
  // One step of online adaptive refinement of the region only, for a closed current on
  // tetrahedra: the region error indicators are zero on the environment; nonconforming
  // refinement of the marked region elements leaves the environment and the interface DOFs
  // as condensed, so the saved model still applies without factoring the environment,
  // matches a monolith on the refined mesh, and the estimated region error decreases.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_mag_amr_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", 1.0}},
                       {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
                   {"Boundaries", {}},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  // J = curl (0, 0, phi) with phi = (x (1 - x) y (1 - y) z (1 - z))^2: divergence-free with
  // no normal component, and integrated exactly, so that it is orthogonal to the discrete
  // gradients (it closes).
  mfem::VectorFunctionCoefficient jc(3,
                                     [](const mfem::Vector &x, mfem::Vector &j)
                                     {
                                       const double a = x(0) * (1.0 - x(0)),
                                                    b = x(1) * (1.0 - x(1)),
                                                    c = x(2) * (1.0 - x(2));
                                       j(0) = 2.0 * a * a * b * (1.0 - 2.0 * x(1)) * c * c;
                                       j(1) = -2.0 * a * (1.0 - 2.0 * x(0)) * b * b * c * c;
                                       j(2) = 0.0;
                                     });
  mfem::ND_FECollection fec(order, 3);
  auto current = [&](mfem::ParMesh &pmesh)
  {
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    mfem::ParLinearForm lf(&pfes);
    auto *lfi = new mfem::VectorFEDomainLFIntegrator(jc);
    lfi->SetIntRule(&mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, 2 * order + 12));
    lf.AddDomainIntegrator(lfi);
    lf.Assemble();
    Vector J(pfes.GetTrueVSize());
    lf.ParallelAssemble(J);
    return J;
  };
  {
    IoData iodata = make_config("Offline");
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4)));
    SubstructuringSolver off(iodata, mesh);
    off.CondenseEnvironment();
    (void)CurrentEnergies(off, {1}, {current(mesh.back()->Get())});
  }
  IoData iodata_on = make_config("Online");
  auto count_env = [](const mfem::ParMesh &m)
  {
    long c = 0;
    for (int e = 0; e < m.GetNE(); e++)
    {
      c += (m.GetAttribute(e) == 2);
    }
    MPI_Allreduce(MPI_IN_PLACE, &c, 1, MPI_LONG, MPI_SUM, m.GetComm());
    return c;
  };
  auto solve = [&](std::vector<std::unique_ptr<Mesh>> &mesh, double &energy)
  {
    SubstructuringSolver on(iodata_on, mesh);
    on.CondenseEnvironment();
    std::vector<Vector> region_fields;
    const mfem::DenseMatrix E =
        on.MagnetostaticEnergyMatrix({1}, {current(mesh.back()->Get())}, {}, {}, {},
                                     nullptr, nullptr, 0, &region_fields);
    CHECK_FALSE(on.EnvironmentFactored());
    energy = E(0, 0);
    return on.RegionErrorIndicator(region_fields, {0.5 * E(0, 0)});
  };

  // Online solve and region indicators on the initial mesh.
  std::vector<std::unique_ptr<Mesh>> mesh0;
  mesh0.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4, false, true)));
  auto &pm0 = mesh0.back()->Get();
  double e0 = 0.0;
  const ErrorIndicator ind0 = solve(mesh0, e0);
  REQUIRE(ind0.Local().Size() == pm0.GetNE());
  double env_max = 0.0, reg_max = 0.0;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    double &m = (pm0.GetAttribute(e) == 2) ? env_max : reg_max;
    m = std::max(m, ind0.Local()(e));
  }
  Mpi::GlobalMax(1, &env_max, pm0.GetComm());
  Mpi::GlobalMax(1, &reg_max, pm0.GetComm());
  CHECK(env_max == 0.0);
  CHECK(reg_max > 0.0);

  // Refine the region elements with the largest indicators (nonconforming, no level
  // constraint), on an identical copy of the mesh.
  auto pm1 = MakeWavyTetSplit(4, false, true);
  mfem::Array<int> marked;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    if (ind0.Local()(e) >= 0.3 * reg_max)
    {
      marked.Append(e);
    }
  }
  const long n_env = count_env(*pm1);
  pm1->GeneralRefinement(marked, 1, 0);
  CHECK(count_env(*pm1) == n_env);  // the environment is not refined
  std::vector<std::unique_ptr<Mesh>> mesh1;
  mesh1.push_back(std::make_unique<Mesh>(std::move(pm1)));
  double e1 = 0.0;
  const ErrorIndicator ind1 = solve(mesh1, e1);

  // Monolith on the refined mesh: u = (K + eps M)^-1 J and E = 2 J^T u - u^T K u, with
  // SubstructuringSolver's regularization.
  auto &pm = mesh1.back()->Get();
  mfem::ParFiniteElementSpace pfes(&pm, &fec);
  const Vector J = current(pm);
  mfem::Vector nu_by_attr(2), mass_by_attr(2);
  nu_by_attr(0) = 1.0;
  nu_by_attr(1) = 0.25;
  mass_by_attr = nu_by_attr;
  mass_by_attr *= kMagRegularization * RegularizationScale(pm);
  mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
  mfem::ParBilinearForm asolve(&pfes), aenergy(&pfes);
  asolve.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
  aenergy.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.Assemble();
  asolve.Finalize();
  aenergy.Assemble();
  aenergy.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Ksolve(asolve.ParallelAssemble()),
      Kpure(aenergy.ParallelAssemble());
  mfem::HypreAMS ams(*Ksolve, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ksolve);
  pcg.SetTol(1.0e-14);
  pcg.SetMaxIter(4000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  Vector u(pfes.GetTrueVSize()), t(pfes.GetTrueVSize());
  u = 0.0;
  pcg.Mult(J, u);
  Kpure->Mult(u, t);
  const double e_mono =
      2.0 * mfem::InnerProduct(pm.GetComm(), J, u) - mfem::InnerProduct(pm.GetComm(), u, t);
  CAPTURE(e0, e1, e_mono);
  CHECK(std::abs(e1 - e_mono) <= 1.0e-8 * std::abs(e_mono));
  CHECK(std::abs(e1 - e0) >= 1.0e3 * std::abs(e1 - e_mono));  // the refinement matters
  CAPTURE(ind0.Norml2(pm0.GetComm()), ind1.Norml2(pm.GetComm()));
  CHECK(ind1.Norml2(pm.GetComm()) < ind0.Norml2(pm0.GetComm()));
}

TEST_CASE("SubstructuringSolver region error indicators across sheets",
          "[substructure][Serial][Parallel]")
{
  // As in the regular solve, the flux recovery of the region indicators is discontinuous
  // across the superconductor sheets: those of a London sheet on y = 0.5 equal the
  // indicators of Palace's estimator, cracked at the sheet, on a serial mesh of the region.
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permeability", 3.0}},
          {{"Attributes", {2}}, {"Permeability", 7.0}}}}}},
      {"Boundaries",
       {{"Superconductor",
         {{{"Attributes", {4}}, {"PenetrationDepth", 0.2}, {"Thickness", 0.05}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Linear", {{"EstimatorTol", 1.0e-13}, {"EstimatorMaxIts", 1000}}},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);
  RegionDesign design;
  design.sheet = true;
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 4, 4, design)));
  SubstructuringSolver sub(iodata, mesh);
  mfem::VectorFunctionCoefficient ac(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1) + 2.0 * x(2));
                                       v(1) = std::cos(3.0 * x(0)) * x(2);
                                       v(2) = std::exp(x(0)) * x(1);
                                     });
  auto project = [&](mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction gf(&fes);
    gf.ProjectCoefficient(ac);
    Vector u(fes.GetTrueVSize());
    gf.GetTrueDofs(u);
    return u;
  };
  mfem::ND_FECollection fec(1, 3);
  mfem::ParFiniteElementSpace pfes(&mesh.back()->Get(), &fec);
  const ErrorIndicator ind = sub.RegionErrorIndicator({project(pfes)}, {1.0});

  // The reference, on every rank.
  mfem::Mesh serial = MakeGradedSplitSerial(4, 4, 4, design);
  mfem::Array<int> attrs({1});
  mfem::Mesh region_serial = mfem::SubMesh::CreateFromDomain(serial, attrs);
  std::vector<std::unique_ptr<Mesh>> rmesh;
  rmesh.push_back(std::make_unique<Mesh>(MPI_COMM_SELF, region_serial));
  MaterialOperator mat_op(iodata, *rmesh[0]);
  FiniteElementSpace rfes(*rmesh[0], &fec);
  auto nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      1, 3, 1, MultigridCoarsening::LINEAR, false);
  auto nd_fespaces =
      fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(1, rmesh, nd_fecs);
  mfem::RT_FECollection rt_fec(0, 3);
  FiniteElementSpace rt_fespace(*rmesh[0], &rt_fec);
  CurlFluxErrorEstimator<Vector> estimator(mat_op, rt_fespace, nd_fespaces, 1.0e-13, 1000,
                                           0, false, {4});
  Vector b(rt_fespace.GetTrueVSize());
  rt_fespace.GetDiscreteInterpolator(rfes).Mult(project(rfes.Get()), b);
  ErrorIndicator ref;
  estimator.AddErrorIndicator(b, 1.0, ref);

  // Element by element, matched by their centers.
  auto &rm = rmesh[0]->Get();
  auto &pm = mesh.back()->Get();
  double err = 0.0, scale = 0.0;
  int matched = 0;
  mfem::Vector c, rc;
  for (int e = 0; e < pm.GetNE(); e++)
  {
    if (pm.GetAttribute(e) != 1)
    {
      continue;
    }
    pm.GetElementCenter(e, c);
    for (int r = 0; r < rm.GetNE(); r++)
    {
      rm.GetElementCenter(r, rc);
      rc -= c;
      if (rc.Normlinf() < 1.0e-12)
      {
        err = std::max(err, std::abs(ind.Local()(e) - ref.Local()(r)));
        scale = std::max(scale, std::abs(ref.Local()(r)));
        matched++;
        break;
      }
    }
  }
  Mpi::GlobalMax(1, &err, Mpi::World());
  Mpi::GlobalMax(1, &scale, Mpi::World());
  Mpi::GlobalSum(1, &matched, Mpi::World());
  CHECK(matched == rm.GetNE());
  CHECK(scale > 0.0);
  CHECK(err <= 1.0e-10 * scale);
}

TEST_CASE("SubstructuringSolver region error indicators",
          "[substructure][Serial][Parallel]")
{
  // The region's indicators (on the region submesh) must equal those of Palace's estimator
  // on a serial mesh of the region alone, on any partition, also with the face DOFs of
  // second-order Nédélec and first-order Raviart-Thomas spaces (in parallel, MFEM's
  // submesh spaces need mfem/mfem#5518).
  const int order = GENERATE(1, 2);
  const bool magnetostatic = GENERATE(false, true);
  CAPTURE(order, magnetostatic);
  json config = {
      {"Problem",
       {{"Type", magnetostatic ? "Magnetostatic" : "Electrostatic"},
        {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"Permeability", 3.0}},
          {{"Attributes", {2}}, {"Permittivity", 5.0}, {"Permeability", 7.0}}}}}},
      {"Boundaries", {}},
      {"Solver",
       {{"Order", order},
        {"Linear", {{"EstimatorTol", 1.0e-13}, {"EstimatorMaxIts", 1000}}},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);
  const int nx = 4;
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(nx)));
  SubstructuringSolver sub(iodata, mesh);

  // A smooth potential, V or A, on the parent and on a serial mesh of the region.
  mfem::FunctionCoefficient vc(
      [](const mfem::Vector &x)
      { return std::sin(x(0) + 2.0 * x(1)) * std::cos(3.0 * x(2)); });
  mfem::VectorFunctionCoefficient ac(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1) + 2.0 * x(2));
                                       v(1) = std::cos(3.0 * x(0)) * x(2);
                                       v(2) = std::exp(x(0)) * x(1);
                                     });
  auto make_fec = [&]() -> std::unique_ptr<mfem::FiniteElementCollection>
  {
    if (magnetostatic)
    {
      return std::make_unique<mfem::ND_FECollection>(order, 3);
    }
    return std::make_unique<mfem::H1_FECollection>(order, 3);
  };
  auto project = [&](mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction gf(&fes);
    if (magnetostatic)
    {
      gf.ProjectCoefficient(ac);
    }
    else
    {
      gf.ProjectCoefficient(vc);
    }
    Vector u(fes.GetTrueVSize());
    gf.GetTrueDofs(u);
    return u;
  };
  const auto fec = make_fec();
  mfem::ParFiniteElementSpace pfes(&mesh.back()->Get(), fec.get());
  const ErrorIndicator ind = sub.RegionErrorIndicator({project(pfes)}, {1.0});

  // The reference, on every rank.
  mfem::Mesh serial = MakeWavyTetSerial(nx);
  mfem::Array<int> attrs({1});
  mfem::Mesh region_serial = mfem::SubMesh::CreateFromDomain(serial, attrs);
  std::vector<std::unique_ptr<Mesh>> rmesh;
  rmesh.push_back(std::make_unique<Mesh>(MPI_COMM_SELF, region_serial));
  MaterialOperator mat_op(iodata, *rmesh[0]);
  ErrorIndicator ref;
  const auto rfec = make_fec();
  FiniteElementSpace rfes(*rmesh[0], rfec.get());
  const Vector ur = project(rfes.Get());
  if (magnetostatic)
  {
    auto nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
        order, 3, 1, MultigridCoarsening::LINEAR, false);
    auto nd_fespaces =
        fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(1, rmesh, nd_fecs);
    mfem::RT_FECollection rt_fec(order - 1, 3);
    FiniteElementSpace rt_fespace(*rmesh[0], &rt_fec);
    CurlFluxErrorEstimator<Vector> estimator(mat_op, rt_fespace, nd_fespaces, 1.0e-13, 1000,
                                             0, false);
    Vector b(rt_fespace.GetTrueVSize());
    rt_fespace.GetDiscreteInterpolator(rfes).Mult(ur, b);
    estimator.AddErrorIndicator(b, 1.0, ref);
  }
  else
  {
    mfem::ND_FECollection nd_fec(order, 3);
    FiniteElementSpace nd_fespace(*rmesh[0], &nd_fec);
    auto rt_fecs = fem::ConstructFECollections<mfem::RT_FECollection>(
        order - 1, 3, 1, MultigridCoarsening::LINEAR, false);
    auto rt_fespaces =
        fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(1, rmesh, rt_fecs);
    GradFluxErrorEstimator<Vector> estimator(mat_op, nd_fespace, rt_fespaces, 1.0e-13, 1000,
                                             0, false);
    Vector e(nd_fespace.GetTrueVSize());
    e = 0.0;
    nd_fespace.GetDiscreteInterpolator(rfes).AddMult(ur, e, -1.0);
    estimator.AddErrorIndicator(e, 1.0, ref);
  }

  // Element by element, matched by their centers.
  auto center = [](const mfem::Mesh &m, int e)
  {
    mfem::Vector c;
    const_cast<mfem::Mesh &>(m).GetElementCenter(e, c);
    return std::array<double, 3>{c(0), c(1), c(2)};
  };
  const auto &rm = rmesh[0]->Get();
  const auto &pm = mesh.back()->Get();
  double err = 0.0, scale = 0.0, env = 0.0, unmatched = 0.0;
  for (int e = 0; e < pm.GetNE(); e++)
  {
    if (pm.GetAttribute(e) != 1)
    {
      env = std::max(env, std::abs(ind.Local()(e)));
      continue;
    }
    const auto c = center(pm, e);
    int match = -1;
    for (int r = 0; r < rm.GetNE() && match < 0; r++)
    {
      const auto q = center(rm, r);
      match = (std::hypot(c[0] - q[0], c[1] - q[1], c[2] - q[2]) < 1.0e-10) ? r : -1;
    }
    if (match < 0)
    {
      unmatched += 1.0;
      continue;
    }
    err = std::max(err, std::abs(ind.Local()(e) - ref.Local()(match)));
    scale = std::max(scale, ref.Local()(match));
  }
  Mpi::GlobalMax(1, &err, Mpi::World());
  Mpi::GlobalMax(1, &scale, Mpi::World());
  Mpi::GlobalMax(1, &env, Mpi::World());
  Mpi::GlobalSum(1, &unmatched, Mpi::World());
  CAPTURE(err, scale);
  CHECK(unmatched == 0.0);
  CHECK(scale > 0.0);
  CHECK(env == 0.0);
  CHECK(err <= 1.0e-8 * scale);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver geometric redesign of the region",
                 "[substructure][Serial][Parallel][GPU]")
{
  // A real design change inside the region between the offline and online runs: the region
  // conductor (terminal 1) changes shape, a dielectric block inside the region changes
  // shape and material, and the region is re-meshed at a different resolution, while the
  // environment and the interface Gamma are unchanged. The online run reuses the saved
  // model (without factoring the environment) and must reproduce a monolithic solve on the
  // redesigned mesh.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_redesign_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode, double eps_block)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {3}}, {"Permittivity", eps_block}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1, 3}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  const std::vector<int> terms = {1, 2};
  const int b = 5, n = 6;  // environment and interface: fixed across designs
  // Design A (offline) and design B (online): different conductor, dielectric and mesh.
  const RegionDesign design_a{0.5, {0.25, 0.5, 0.5}}, design_b{0.84, {0.34, 0.34, 0.67}};
  const double eps_a = 4.0, eps_b = 6.0;

  IoData iodata_off = make_config("Offline", eps_a);
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, b, n, design_a)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const mfem::DenseMatrix C_a = off.CapacitanceMatrix(terms);

  IoData iodata_on = make_config("Online", eps_b);
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(6, b, n, design_b)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix C_b = on.CapacitanceMatrix(terms);
  CHECK_FALSE(on.EnvironmentFactored());

  // Monolithic reference on the redesigned mesh: C_ij = u_i^T K u_j, terminal i at 1 V.
  auto &pmesh = mesh_on.back()->Get();
  mfem::H1_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::Vector eps_by_attr(pmesh.attributes.Max());
  eps_by_attr = 1.0;
  eps_by_attr(1) = 10.0;   // environment (attr 2)
  eps_by_attr(2) = eps_b;  // dielectric block (attr 3)
  mfem::PWConstCoefficient eps(eps_by_attr);
  mfem::ParBilinearForm a(&pfes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> K(a.ParallelAssemble());
  const int maxb = pmesh.bdr_attributes.Max();
  std::vector<mfem::Vector> u(2);
  for (int i = 0; i < 2; i++)
  {
    mfem::ParGridFunction x(&pfes);
    x = 0.0;
    mfem::Array<int> drive(maxb), ess_bdr(maxb), ess_tdofs;
    drive = 0;
    drive[i] = 1;  // terminal i+1 is boundary attribute i+1
    ess_bdr = 0;
    ess_bdr[0] = ess_bdr[1] = 1;
    mfem::ConstantCoefficient one(1.0);
    x.ProjectBdrCoefficient(one, drive);
    pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
    mfem::ParLinearForm rhs(&pfes);
    rhs = 0.0;
    mfem::OperatorPtr A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdofs, x, rhs, A, X, B);
    mfem::HypreBoomerAMG amg(*A.As<mfem::HypreParMatrix>());
    amg.SetPrintLevel(0);
    mfem::HyprePCG pcg(*A.As<mfem::HypreParMatrix>());
    pcg.SetTol(1.0e-14);
    pcg.SetMaxIter(1000);
    pcg.SetPrintLevel(0);
    pcg.SetPreconditioner(amg);
    pcg.Mult(B, X);
    a.RecoverFEMSolution(X, rhs, x);
    x.GetTrueDofs(u[i]);
  }
  mfem::DenseMatrix C_mono(2);
  mfem::Vector Ku(pfes.GetTrueVSize());
  for (int j = 0; j < 2; j++)
  {
    K->Mult(u[j], Ku);
    for (int i = 0; i < 2; i++)
    {
      C_mono(i, j) = mfem::InnerProduct(Mpi::World(), u[i], Ku);
    }
  }
  double d = 0.0, m = 0.0, change = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      d = std::max(d, std::abs(C_b(i, j) - C_mono(i, j)));
      m = std::max(m, std::abs(C_mono(i, j)));
      change = std::max(change, std::abs(C_b(i, j) - C_a(i, j)));
    }
  }
  CHECK(d <= 1.0e-9 * m);
  CHECK(change >= 1.0e-2 * m);  // the redesign really changed the capacitance
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver surface currents close through connected conductors",
                 "[substructure][Serial][Parallel]")
{
  // The sheet current J_s = x on the plane y = 1/2 runs from the PEC wall x = 0 to the PEC
  // wall x = 1. It closes when the other faces are PEC too, connecting the two walls, or
  // when the plane is a London sheet, and not when only the two walls are PEC: then it
  // leaves one conductor for another. On its strip x < 3/8 of the region, it ends in free
  // space, also on the mesh of the region.
  const int order = 1;
  const std::string model_path = (temp_dir / "substruct_closure_model.bin").string();
  RegionDesign design;
  design.sheet = true;
  auto energy = [&](const std::vector<int> &pec, bool london, bool strip = false,
                    const std::string &mode = "", mfem::Mesh *serial = nullptr)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
        {"Boundaries", {{"PEC", {{"Attributes", pec}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    if (london)
    {
      config["Boundaries"]["Superconductor"] = {
          {{"Attributes", {4}}, {"PenetrationDepth", 0.2}, {"Thickness", 0.05}}};
    }
    else
    {
      config["Boundaries"]["SurfaceCurrent"] = {
          {{"Index", 1}, {"Attributes", {4}}, {"Direction", "+X"}}};
    }
    if (!mode.empty())
    {
      config["Solver"]["Substructuring"]["Mode"] = mode;
      config["Solver"]["Substructuring"]["SaveModel"] = model_path;
    }
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(serial ? std::make_unique<Mesh>(
                                std::make_unique<mfem::ParMesh>(Mpi::World(), *serial))
                          : std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
    auto &pmesh = mesh.back()->Get();
    mfem::ND_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    const int maxb = pmesh.bdr_attributes.Max();
    mfem::Array<int> pec_marker(maxb), pec_tdofs, sheet_marker(maxb);
    pec_marker = 0;
    for (int a : pec)
    {
      pec_marker[a - 1] = 1;
    }
    pfes.GetEssentialTrueDofs(pec_marker, pec_tdofs);
    sheet_marker = 0;
    sheet_marker[3] = 1;
    mfem::VectorFunctionCoefficient c(3,
                                      [strip](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v = 0.0;
                                        v(0) = (!strip || x(0) < 0.375) ? 1.0 : 0.0;
                                      });
    mfem::ParLinearForm f(&pfes);
    f.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(c), sheet_marker);
    f.Assemble();
    Vector J(pfes.GetTrueVSize());
    f.ParallelAssemble(J);
    J.SetSubVector(pec_tdofs, 0.0);
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return CurrentEnergies(ss, {1}, {J})(0, 0);
  };
  CHECK(energy({1, 2, 3}, false) > 0.0);
  CHECK(energy({1, 2}, true) > 0.0);
  CHECK(energy({1, 2, 3}, false, false, "Offline") > 0.0);
  mfem::Mesh region_mesh(RegionMeshPath(model_path).c_str());
  CHECK(energy({1, 2, 3}, false, false, "Online", &region_mesh) > 0.0);
#if defined(MFEM_USE_EXCEPTIONS)
  CHECK_THROWS_WITH(energy({1, 2}, false),
                    Catch::Matchers::ContainsSubstring("does not close"));
  CHECK_THROWS_WITH(energy({1, 2, 3}, false, true),
                    Catch::Matchers::ContainsSubstring("does not close"));
  CHECK_THROWS_WITH(energy({1, 2, 3}, false, true, "Online", &region_mesh),
                    Catch::Matchers::ContainsSubstring("does not close"));
#endif
}

TEST_CASE("SubstructuringSolver magnetostatic energies scale with the permeability and "
          "the mesh",
          "[substructure][Serial][Parallel]")
{
  // For a fixed dual vector J, scaling the permeability by m scales the curl-curl operator
  // by 1/m and the energy by m, and scaling the mesh by s scales the curl-curl operator by
  // 1/s and the mass by s, so the energy by s. The regularization scales as the curl-curl
  // operator, so this holds whatever the scales.
  const int order = 1;
  auto energy = [order](double m, double s, Vector &J)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", m}},
            {{"Attributes", {2}}, {"Permeability", 4.0 * m}}}}}},
        {"Boundaries", json::object()},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);
    auto pmesh = MakeSplitCube(4);
    mfem::VectorFunctionCoefficient scale(3, [s](const mfem::Vector &x, mfem::Vector &y)
                                          { add(s, x, 0.0, x, y); });
    pmesh->Transform(scale);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    if (J.Size() == 0)
    {
      J = ClosedCurrent(*mesh.back(), order);
    }
    return CurrentEnergies(ss, {1}, {J})(0, 0);
  };
  Vector J;
  const double E = energy(1.0, 1.0, J), E_m = energy(1.0e5, 1.0, J),
               E_s = energy(1.0, 1.0e3, J);
  CAPTURE(E, E_m, E_s);
  CHECK(E > 0.0);
  CHECK(std::abs(E_m - 1.0e5 * E) <= 1.0e-10 * 1.0e5 * E);
  CHECK(std::abs(E_s - 1.0e3 * E) <= 1.0e-10 * 1.0e3 * E);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver magnetostatic cross-run re-meshing",
                 "[substructure][Serial][Parallel][GPU]")
{
  // H(curl) region re-meshing: the environment (attr 2, b=5) is identical between a coarse-
  // region offline run and a re-meshed (finer) online run, so the loaded S_E, S^K and
  // source modes, reordered by signed edge-signature matching, give the energy of a closed
  // current without an environment solve, and its field, as a fresh condensation on the
  // online mesh.
  const auto res = GENERATE(std::make_pair(4, 8), std::make_pair(6, 6));
  const int a_off = res.first, a_on = res.second;
  const int order = GENERATE(1, 2);
  CAPTURE(a_off, a_on, order);
  auto make_config = [&](const std::string &mode, const std::string &path)
  {
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", 1.0}},
                       {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"Mode", mode},
                       {"SaveModel", path}}}}}};
    return IoData(config, false);
  };
  auto solve = [&](const std::string &mode, const std::string &path, int a, Vector *u)
  {
    IoData io = make_config(mode, path);
    std::vector<std::unique_ptr<Mesh>> m;
    m.push_back(std::make_unique<Mesh>(MakeGradedSplit(a, 5, 6)));
    SubstructuringSolver s(io, m);
    s.CondenseEnvironment();
    std::vector<Vector> fields;
    const double E = CurrentEnergies(s, {1}, {ClosedCurrent(*m.back(), order)}, &fields,
                                     u ? 1 : 0)(0, 0);
    if (u)
    {
      *u = fields[0];
    }
    return E;
  };

  // Offline on the coarse region (the model with the source modes), a fresh condensation on
  // the re-meshed region, and the online run on it with the coarse-region model.
  const std::string model_path = (temp_dir / "substruct_mag_remesh.bin").string();
  (void)solve("Offline", model_path, a_off, nullptr);
  Vector u_ref, u_on;
  const double E_ref =
      solve("Offline", (temp_dir / "substruct_mag_ref.bin").string(), a_on, &u_ref);
  const double E_on = solve("Online", model_path, a_on, &u_on);
  CHECK(E_ref > 0.0);
  CHECK(std::abs(E_on - E_ref) <= 1.0e-10 * E_ref);
  Vector d(u_on);
  d -= u_ref;
  CHECK(d.Norml2() <= 1.0e-9 * u_ref.Norml2());
}

TEST_CASE("Rows of an interface operator pass through rank 0 in chunks",
          "[substructure][Serial][Parallel]")
{
  // Rank q holds q + 2 rows of width 3; a chunk of 1 or 6 doubles holds 1 or 2 rows, the
  // largest one all of them.
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), nranks = Mpi::Size(comm), w = 3;
  const std::size_t chunk = GENERATE(1, 6, 1000);
  auto value = [](int q, int r, int c) { return 100.0 * q + 10.0 * r + c; };
  std::vector<int> cnt(nranks);
  for (int q = 0; q < nranks; q++)
  {
    cnt[q] = q + 2;
  }
  std::vector<double> rows(static_cast<std::size_t>(cnt[rank]) * w), expected;
  for (int r = 0; r < cnt[rank]; r++)
  {
    for (int c = 0; c < w; c++)
    {
      rows[r * w + c] = value(rank, r, c);
    }
  }

  // To rank 0: every row once, in rank order, in chunks of at most the chunk.
  std::vector<double> got;
  int most = 0;
  substructuring::RowsToRoot(comm, rows, w, cnt, chunk,
                             [&](const double *x, int n)
                             {
                               got.insert(got.end(), x, x + n * w);
                               most = std::max(most, n);
                             });
  for (int q = 0; Mpi::Root(comm) && q < nranks; q++)
  {
    for (int r = 0; r < cnt[q]; r++)
    {
      for (int c = 0; c < w; c++)
      {
        expected.push_back(value(q, r, c));
      }
    }
  }
  CHECK(got == expected);
  CHECK(most <= std::max<int>(1, static_cast<int>(chunk) / w));

  // From rank 0: rank q takes the rows 3 i + q (i < 2 q + 3) of a source, in that order.
  auto need = [](int q)
  {
    std::vector<int> a;
    for (int i = 0; i < 2 * q + 3; i++)
    {
      a.push_back(3 * i + q);
    }
    return a;
  };
  auto read = [](int a, double *x)
  {
    for (int c = 0; c < 3; c++)
    {
      x[c] = a + 0.25 * c;
    }
  };
  substructuring::RowsFromRoot from(comm, w, chunk);
  if (Mpi::Root(comm))
  {
    from.Serve(need, read);
  }
  else
  {
    // Each chunk holds the next rows of the chunk, or the rest (probed before Next takes
    // it).
    const std::vector<int> a = need(rank);
    const std::size_t per = std::max<std::size_t>(1, chunk / w);
    for (std::size_t k = 0; k < a.size(); k++)
    {
      if (k % per == 0)
      {
        MPI_Status status;
        int n = 0;
        MPI_Probe(0, 0, comm, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &n);
        CHECK(n == static_cast<int>(std::min(per, a.size() - k)) * w);
      }
      std::vector<double> x(w), ref(w);
      from.Next(x.data());
      read(a[k], ref.data());
      CHECK(x == ref);
    }
  }
}

}  // namespace palace
