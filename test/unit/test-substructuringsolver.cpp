// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "drivers/drivensolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "fixtures.hpp"
#include "linalg/directschur.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/lowertriangle.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/curlcurloperator.hpp"
#include "models/drivensubstructure.hpp"
#include "models/environmentoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/regionsolver.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/substructuringsolver.hpp"
#include "models/substructuringsolverimpl.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

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

TEST_CASE("SubstructuringSolver environment defaults to the attributes not in the region",
          "[substructure][Serial][Parallel]")
{
  // Without "Environment", the environment is the domain attributes of the mesh not in the
  // region: the capacitance matrix is that with the environment given, and a region with
  // every attribute leaves no environment.
  auto capacitance = [](const std::vector<int> &region, bool environment)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 4.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", region}}}, {"SaveModel", false}}}}}};
    if (environment)
    {
      config["Solver"]["Substructuring"]["Environment"] = {{"Attributes", {2}}};
    }
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  const mfem::DenseMatrix C = capacitance({1}, true), C_default = capacitance({1}, false);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(std::abs(C_default(i, j) - C(i, j)) <= 1.0e-12 * std::abs(C(i, i)));
    }
  }
#if defined(MFEM_USE_EXCEPTIONS)
  CHECK_THROWS_WITH(capacitance({1, 2}, false),
                    Catch::Matchers::ContainsSubstring("the environment is empty"));
#endif
}

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
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"SaveModel", false}}}}}};
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
    const double *h_u = u.HostRead(), *h_full = u_full.HostRead();
    for (int i = 0; i < pfes.GetTrueVSize(); i++)
    {
      double e = h_u[i] - h_full[i];
      num += e * e;
      den += h_full[i] * h_full[i];
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
            {"Environment", {{"Attributes", environment}}},
            {"SaveModel", false}}}}}};
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
  const std::string model_path =
      (temp_dir / "substruct_region_refined" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  auto no_ground = region_mesh();
  condense(make_config("Offline", none, true), MakeSplitCube(4));
  CHECK_THROWS_WITH(condense(make_config("Online", uniform, true), refined(region_mesh())),
                    Catch::Matchers::ContainsSubstring("is written refined"));
  // A mesh which is neither, as that of the region without the ground, is rejected.
  CHECK_THROWS_WITH(condense(make_config("Online", none, true), std::move(no_ground)),
                    Catch::Matchers::ContainsSubstring("its mesh changed"));
  CHECK_THROWS_WITH(
      condense(make_config("Online", uniform, true), refined(MakeSplitCube(4))),
      Catch::Matchers::ContainsSubstring("The environment differs"));
}

#if defined(MFEM_USE_STRUMPACK)
TEST_CASE("SubstructuringSolver factors without STRUMPACK compression",
          "[substructure][Serial][Parallel][Long]")
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
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"SaveModel", false}}}}}};
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
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"SaveModel", false}}}}}};
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
         {{"Region", {{"Attributes", {1}}}},
          {"Environment", {{"Attributes", {2}}}},
          {"SaveModel", false}}}}}};
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
         {{"Region", {{"Attributes", {1}}}},
          {"Environment", {{"Attributes", {2}}}},
          {"SaveModel", false}}}}}};
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
  // On the device the solution and its product are read there (not as host entries).
  const double e_mono = 0.5 * mfem::InnerProduct(Mpi::World(), u_full, t);

  CHECK(e_sub > 0.0);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-10 * std::abs(e_mono));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver offline/online model reuse",
                 "[substructure][Serial][Parallel][GPU]")
{
  const std::string model_path =
      (temp_dir / "substruct_model_roundtrip" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path);
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
  const std::string model_path =
      (temp_dir / "substruct_env_fingerprint" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_mirror" / "environment.model").string();
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
                       {"Environment", {{"Attributes", {2, 3}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  // An environment half moved to the region, on the same mesh and materials.
  IoData half = make_config("Online", 2.0, 5.0, 4);
  half.solver.substructuring->region_attributes = {1, 3};
  half.solver.substructuring->environment_attributes = {2};
  CHECK_THROWS_WITH(condense(half),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver mesh of the region rejects a changed environment configuration",
    "[substructure][Serial][Parallel][GPU]")
{
  // The mesh of the region holds little of the environment, so its configuration tells a
  // changed environment there: another material or boundary condition on it is rejected,
  // while the region may change.
  const std::string model_path =
      (temp_dir / "substruct_region_config" / "environment.model").string();
  auto make_config =
      [&](const std::string &mode, double eps_region, double eps_env, bool ground)
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  condense(make_config("Offline", 1.0, 10.0, false), MakeSplitCube(4));
  CHECK_NOTHROW(condense(make_config("Online", 1.0, 10.0, false), region_mesh()));
  CHECK_NOTHROW(condense(make_config("Online", 2.5, 10.0, false), region_mesh()));
  const auto changed = Catch::Matchers::ContainsSubstring(
      "its attributes, materials, boundary conditions or discretization changed");
  CHECK_THROWS_WITH(condense(make_config("Online", 1.0, 12.0, false), region_mesh()),
                    changed);
  CHECK_THROWS_WITH(condense(make_config("Online", 1.0, 10.0, true), region_mesh()),
                    changed);
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver saved model rejects re-tagged environment boundaries",
    "[substructure][Serial][Parallel][GPU]")
{
  // The grounded faces of the environment extended, with the configuration unchanged: the
  // moments of the boundary attributes next to the environment tell the mesh apart.
  const std::string model_path =
      (temp_dir / "substruct_retag" / "environment.model").string();
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
         {{"Ground", {{"Attributes", {4}}}},
          {"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
    return IoData(config, false);
  };
  // The unit cube split at x = 1/2 (attributes 1 and 2), with the faces x = 0 (1), x = 1
  // (2), y = 0 with x > x0 (4), and the others (3).
  auto condense = [](const IoData &iodata, double x0)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
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
      double xc[2] = {0.0, 0.0};
      for (int j = 0; j < vtx.Size(); j++)
      {
        xc[0] += serial.GetVertex(vtx[j])[0] / vtx.Size();
        xc[1] += serial.GetVertex(vtx[j])[1] / vtx.Size();
      }
      serial.SetBdrAttribute(b, (xc[0] < 1e-9)                 ? 1
                                : (xc[0] > 1.0 - 1e-9)         ? 2
                                : (xc[1] < 1e-9 && xc[0] > x0) ? 4
                                                               : 3);
    }
    serial.SetAttributes();
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
  };
  condense(make_config("Offline"), 0.75);
  CHECK_NOTHROW(condense(make_config("Online"), 0.75));
  CHECK_THROWS_WITH(condense(make_config("Online"), 0.5),
                    Catch::Matchers::ContainsSubstring("its mesh changed"));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "SubstructuringSolver saved model takes a boundary attribute of the region",
    "[substructure][Serial][Parallel]")
{
  // Palace adds boundary elements on the interface, numbered from the largest boundary
  // attribute of the mesh: an attribute added to the region renumbers them, and the model
  // compares the boundary attributes the configuration names only. Online, the grounded
  // faces y = 0 of the region take the new attribute 4.
  const std::string model_path =
      (temp_dir / "substruct_added" / "environment.model").string();
  auto write_mesh = [&](bool added)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
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
      double xc[2] = {0.0, 0.0};
      for (int v : vtx)
      {
        xc[0] += serial.GetVertex(v)[0] / vtx.Size();
        xc[1] += serial.GetVertex(v)[1] / vtx.Size();
      }
      serial.SetBdrAttribute(b, (xc[0] < 1e-9)                           ? 1
                                : (xc[0] > 1.0 - 1e-9)                   ? 2
                                : (added && xc[1] < 1e-9 && xc[0] < 0.5) ? 4
                                                                         : 3);
    }
    serial.SetAttributes();
    const auto path = temp_dir / (added ? "added.mesh" : "base.mesh");
    if (Mpi::Root(Mpi::World()))
    {
      serial.Save(path.string());
    }
    Mpi::Barrier(Mpi::World());
    return path.string();
  };
  auto capacitance = [&](const std::string &mode, const std::string &mesh_path,
                         const std::vector<int> &ground)
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", mesh_path}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries",
         {{"Ground", {{"Attributes", ground}}},
          {"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.CapacitanceMatrix({1, 2});
  };
  const mfem::DenseMatrix C_off = capacitance("Offline", write_mesh(false), {3});
  mfem::DenseMatrix C_on;
  REQUIRE_NOTHROW(C_on = capacitance("Online", write_mesh(true), {3, 4}));
  C_on -= C_off;
  CHECK(C_on.MaxMaxNorm() <= 1.0e-10 * C_off.MaxMaxNorm());
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver saved model of another problem or order",
                 "[substructure][Serial][Parallel]")
{
  // A model is loaded by a run of its problem type and order only, and a file that is not
  // one is reported as such.
  const std::string model_path =
      (temp_dir / "substruct_problem" / "environment.model").string();
  auto make_config =
      [&model_path](const std::string &type, const std::string &mode, int order)
  {
    const json terminals = {{{"Index", 1}, {"Attributes", {1}}},
                            {{"Index", 2}, {"Attributes", {2}}}};
    json config = {
        {"Problem", {{"Type", type}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1, 2}}}}}}},
        {"Boundaries", (type == "Electrostatic") ? json{{"Terminal", terminals}}
                                                 : json{{"PEC", {{"Attributes", {1, 2}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  unwritable.problem.output = (temp_dir / "missing").string();
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
  const std::string model_path =
      (temp_dir / "substruct_capacitance_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path);
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
  const std::string model_path =
      (temp_dir / "substruct_current_orient_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
    J.SetSubVector(pec_tdofs, 0.0);
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

TEST_CASE_METHOD(palace::test::SharedTempDir, "SubstructuringSolver column orderings",
                 "[substructure][Serial][Parallel]")
{
  // Solver.Linear.ColumnOrdering orders the factorizations: with MUMPS, the environment
  // factor of a saved magnetostatic model has no Schur variables, so ParMETIS analyzes it
  // in parallel, and the partial factorization of STRUMPACK takes METIS for it. The energy
  // of a closed current does not depend on the ordering.
  const std::string model_path =
      (temp_dir / "substruct_ordering_model" / "environment.model").string();
  auto energy = [&](const std::string &ordering)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, "Offline", model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_sheet_eps" / "environment.model").string();
  auto make_config = [&](const std::string &mode, bool sheet)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  // x = 1) has no DOFs: the model must hold it (the terminal would otherwise give a zero
  // row of C). A renamed one is an excitation of the environment which differs from the
  // model's, and so is one of a model without its terminal modes: both need the full mesh.
  // An attribute in neither the mesh nor the environment of the model is unknown.
  const std::string model_path =
      (temp_dir / "substruct_terminals" / "environment.model").string();
  auto make_config =
      [&](const std::string &mode, int env_index, const std::vector<int> &attrs = {1})
  {
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains", {{"Materials", {{{"Attributes", {1, 2}}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Ground", {{"Attributes", {3}}}},
          {"Terminal",
           {{{"Index", 1}, {"Attributes", attrs}},
            {{"Index", env_index}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  CHECK_THROWS_WITH(
      capacitance(make_config("Online", 2, {1, 99}), region(), 2),
      Catch::Matchers::ContainsSubstring("Unknown terminal boundary attribute 99"));
  CHECK_THROWS_WITH(
      capacitance(make_config("Online", 5), region(), 5),
      Catch::Matchers::ContainsSubstring("The excitations of the environment differ"));
  // Strip the record of the terminal modes, the last one.
  const auto file = substructuring::ModelFile::Read(model_path, Mpi::World());
  REQUIRE(!file.records.empty());
  REQUIRE(file.records.back().tag == substructuring::ModelFile::Tag("MODE"));
  if (Mpi::Root(Mpi::World()))
  {
    std::filesystem::resize_file(model_path, file.records.back().pos - 12);
  }
  Mpi::Barrier(Mpi::World());
  CHECK_THROWS_WITH(capacitance(make_config("Online", 2), region(), 2),
                    Catch::Matchers::ContainsSubstring("Excitation 2 has no DOFs"));
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
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"SaveModel", false}}}}}};
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
  const std::string model_path =
      (temp_dir / "substruct_new_terminal" / "environment.model").string();
  auto make_config =
      [](const std::string &mode, bool region_terminal, const std::string &path)
  {
    json terminals = {{{"Index", 2}, {"Attributes", {2}}}};
    if (region_terminal)
    {
      terminals.push_back({{"Index", 1}, {"Attributes", {1}}});
    }
    json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
        {"Boundaries", {{"Terminal", terminals}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path);
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
      capacitance(make_config("Offline", true,
                              (temp_dir / "substruct_ref" / "environment.model").string()),
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
  const std::string model_path =
      (temp_dir / "substruct_incomplete_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const double L = substructuring::InterfaceRecords::ReadLength(model_path, Mpi::World());
  double Lc[2] = {L, -L};
  Mpi::GlobalMin(2, Lc, Mpi::World());
  CHECK(L == 2.0);
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
  const std::string model_path =
      (temp_dir / "substruct_nc_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_nomodes_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_remesh_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path);
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
  // On the device the solution and its product are read there (not as host entries).
  const double e_mono = 0.5 * mfem::InnerProduct(Mpi::World(), u_full, t);

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

TEST_CASE_METHOD(palace::test::SharedTempDir, "SubstructuringSolver boundaries of the mesh",
                 "[substructure][Serial]")
{
  // The operators on the mesh of a run leave out the terminal attributes it does not hold
  // and, on the mesh of the region, the flux loops of the environment that the model lists
  // (whatever attributes the mesh holds) and the superconductor attributes it does not
  // hold. The mesh holds the boundary attributes 1 to 6.
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON);
  mfem::ParMesh mesh(Mpi::World(), serial);
  const std::string model = (temp_dir / "substruct_boundaries.model").string();
  using substructuring::ModelFile;
  ModelFile::Create(model, Mpi::World());
  const int env_loop = 2;
  ModelFile::Append(model, Mpi::World(), ModelFile::Tag("ELPS"), sizeof(int),
                    [&](std::ostream &f)
                    { f.write(reinterpret_cast<const char *>(&env_loop), sizeof(int)); });
  auto boundaries = [&](const std::string &type, const json &bdr, bool region_mesh)
  {
    const json sub = {{"Region", {{"Attributes", {1}}}},
                      {"Environment", {{"Attributes", {2}}}},
                      {"ModelFile", model}};
    const json config = {{"Problem", {{"Type", type}, {"Output", "test_output"}}},
                         {"Model", {{"Mesh", "test.msh"}}},
                         {"Domains", {{"Materials", {{{"Attributes", {1, 2}}}}}}},
                         {"Boundaries", bdr},
                         {"Solver", {{"Substructuring", sub}}}};
    return SubstructuringSolver::MeshBoundaries(IoData(config, false), mesh, region_mesh);
  };
  const json terminals = {
      {"Terminal",
       {{{"Index", 1}, {"Attributes", {1, 9}}}, {{"Index", 2}, {"Attributes", {9}}}}},
      {"Ground", {{"Attributes", {3}}}}};
  for (const bool region_mesh : {false, true})
  {
    CAPTURE(region_mesh);
    const auto b = boundaries("Electrostatic", terminals, region_mesh);
    CHECK(b.terminal.size() == 1);
    CHECK(b.terminal.contains(1));
    CHECK((b.terminal.contains(1) && b.terminal.at(1).attributes == std::vector<int>{1}));
  }
  auto loop = [](int index, int film, int hole)
  {
    return json{{"Index", index},
                {"FilmAttributes", {film}},
                {"HoleAttributes", {hole}},
                {"FluxAmounts", {1.0}},
                {"Direction", "+Z"}};
  };
  const json films = {
      {"PEC", {{"Attributes", {6}}}},
      {"FluxLoop", {loop(1, 2, 3), loop(2, 4, 5)}},
      {"Superconductor",
       {{{"Attributes", {2, 9}}, {"PenetrationDepth", 0.1}, {"Thickness", 0.05}},
        {{"Attributes", {8}}, {"PenetrationDepth", 0.1}, {"Thickness", 0.05}}}}};
  auto attributes = [](const config::BoundaryData &b)
  {
    std::vector<int> attrs;
    for (const auto &data : b.superconductor)
    {
      attrs.insert(attrs.end(), data.attributes.begin(), data.attributes.end());
    }
    std::ranges::sort(attrs);
    return attrs;
  };
  // With the sheet of the film of loop 2, which has no superconductor entry.
  const auto full = boundaries("Magnetostatic", films, false);
  CHECK(full.fluxloop.size() == 2);
  CHECK(full.superconductor.size() == 3);
  CHECK(attributes(full) == std::vector<int>{2, 4, 8, 9});
  const auto region = boundaries("Magnetostatic", films, true);
  CHECK(region.fluxloop.size() == 1);
  CHECK(region.fluxloop.contains(1));
  CHECK(region.superconductor.size() == 2);
  CHECK(attributes(region) == std::vector<int>{2, 4});
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver postprocessing of the region",
                 "[substructure][Serial][Parallel]")
{
  // The postprocessing of the region keeps the energies of its domains, the probes in it
  // (also one on the interface), and the surface fluxes and interface dielectrics on its
  // surfaces, and leaves out the others: of the environment, across the interface, and a
  // probe outside the mesh. On the mesh of the region, the sides of the cube (a natural
  // boundary, along which no environment element is written) seem to border the region
  // only: they are left out there too.
  auto flux = [](int index, int attr)
  { return json{{"Index", index}, {"Attributes", {attr}}, {"Type", "Electric"}}; };
  auto dielectric = [](int index, int attr)
  {
    return json{{"Index", index},
                {"Attributes", {attr}},
                {"Thickness", 1.0e-3},
                {"Permittivity", 4.0}};
  };
  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 1.0}},
          {{"Attributes", {2}}, {"Permittivity", 4.0}}}},
        {"Postprocessing",
         {{"Energy",
           {{{"Index", 1}, {"Attributes", {1}}},
            {{"Index", 2}, {"Attributes", {2}}},
            {{"Index", 3}, {"Attributes", {1, 2}}}}}}}}},
      {"Boundaries",
       {{"Terminal",
         {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}},
        {"Postprocessing",
         {{"SurfaceFlux", {flux(1, 1), flux(2, 2), flux(3, 3)}},
          {"Dielectric", {dielectric(1, 1), dielectric(2, 3)}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}},
          {"Environment", {{"Attributes", {2}}}},
          {"SaveModel", false}}}}}};
#if defined(MFEM_USE_GSLIB)
  // Probes need GSLIB.
  config["Domains"]["Postprocessing"]["Probe"] = {
      {{"Index", 1}, {"Center", {0.25, 0.5, 0.5}}},
      {{"Index", 2}, {"Center", {0.75, 0.5, 0.5}}},
      {{"Index", 3}, {"Center", {2.0, 0.5, 0.5}}},
      {{"Index", 4}, {"Center", {0.5, 0.5, 0.5}}}};
#endif
  IoData iodata(config, false);
  CHECK(SubstructuringSolver::Postprocesses(iodata));
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
  SubstructuringSolver ss(iodata, mesh);
  const IoData post = ss.RegionPostprocessing(iodata);
  auto indices = [](const auto &entries)
  {
    std::vector<int> idx;
    for (const auto &[i, data] : entries)
    {
      idx.push_back(i);
    }
    return idx;
  };
  CHECK(indices(post.domains.postpro.energy) == std::vector<int>{1});
#if defined(MFEM_USE_GSLIB)
  CHECK(indices(post.domains.postpro.probe) == std::vector<int>{1, 4});
#endif
  CHECK(indices(post.boundaries.postpro.flux) == std::vector<int>{1});
  CHECK(indices(post.boundaries.postpro.dielectric) == std::vector<int>{1});
  const std::string model_path =
      (temp_dir / "substruct_postprocessing_model" / "environment.model").string();
  test::SetModelRun(config, "Offline", model_path);
  {
    std::vector<std::unique_ptr<Mesh>> full;
    full.push_back(std::make_unique<Mesh>(MakeSplitCube(4)));
    const IoData saving(config, false);
    SubstructuringSolver offline(saving, full);
    offline.CondenseEnvironment();
  }
  test::SetModelRun(config, "Online", model_path);
  const IoData online(config, false);
  mfem::Mesh serial(RegionMeshPath(model_path).c_str());
  std::vector<std::unique_ptr<Mesh>> region;
  region.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  const IoData region_post =
      SubstructuringSolver(online, region).RegionPostprocessing(online);
  CHECK(indices(region_post.boundaries.postpro.flux) == std::vector<int>{1});
  CHECK(indices(region_post.boundaries.postpro.dielectric) == std::vector<int>{1});
  config["Domains"].erase("Postprocessing");
  config["Boundaries"].erase("Postprocessing");
  CHECK_FALSE(SubstructuringSolver::Postprocesses(IoData(config, false)));
}

TEST_CASE("SubstructuringSolver region factored with the environment operator",
          "[substructure][Serial][Parallel][Long]")
{
  // The region is factored with S_E as one sparse system, here with one excitation and an
  // interface of over 1000 DOFs, the dense block of the factorization. Without a Dirichlet
  // boundary, the region operator alone is singular, but not with S_E. The capacitance is
  // that of a monolith, of the terminal at x = 1 against the grounded sides of the
  // environment, with the region grounded at x = 0 or not.
  const bool grounded = GENERATE(true, false);
  CAPTURE(grounded);
  const std::vector<int> ground = grounded ? std::vector<int>{2, 4} : std::vector<int>{4};
  json boundaries = {{"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}},
                     {"Ground", {{"Attributes", ground}}}};
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains",
                  {{"Materials",
                    {{{"Attributes", {1}}, {"Permittivity", 1.0}},
                     {{"Attributes", {2}}, {"Permittivity", 10.0}}}}}},
                 {"Boundaries", boundaries},
                 {"Solver",
                  {{"Order", 2},
                   {"Substructuring",
                    {{"Region", {{"Attributes", {1}}}},
                     {"Environment", {{"Attributes", {2}}}},
                     {"SaveModel", false}}}}}};
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
  json config = {{"Problem", {{"Type", "Electrostatic"}, {"Output", "test_output"}}},
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
                    {{"Region", {{"Attributes", {1, 3}}}},
                     {"Environment", {{"Attributes", {2}}}},
                     {"SaveModel", false}}}}}};
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
                    {"Environment", {{"Attributes", {2}}}},
                    {"SaveModel", false}}}}}};
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
                 "[substructure][Serial][Parallel][Long]")
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
  const std::string model_path =
      (temp_dir / "substruct_london_model" / "environment.model").string();
  auto make_config =
      [&](const std::string &mode, bool save = true, const std::string &path = "")
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path.empty() ? model_path : path);
    if (!save)
    {
      config["Solver"]["Substructuring"]["SaveModel"] = false;
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
    Vector on_sheets;
    full.GetSubVector(sheet_tdofs, on_sheets);
    a[k].SetSubVector(sheet_tdofs, on_sheets);
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
    rhs.SetSubVector(pec_tdofs, 0.0);
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

  // On the mesh of the region, a source without DOFs in it is of the environment, whose
  // modes the model must hold (a model without them reads as complete).
  {
    mfem::Mesh serial(RegionMeshPath(model_path).c_str());
    std::vector<std::unique_ptr<Mesh>> region;
    region.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
    SubstructuringSolver on_region(iodata_on, region);
    on_region.CondenseEnvironment();
    Vector zero(on_region.GetCurlCurlOperator().GetNDSpace().GetTrueVSize());
    zero.UseDevice(true);
    zero = 0.0;
    CHECK_THROWS_WITH(SheetEnergies(on_region, {7}, {zero}),
                      Catch::Matchers::ContainsSubstring(
                          "Excitation 7 has no DOFs in the mesh of the region"));
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
  const std::string model_path =
      (temp_dir / "substruct_current_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
    J[k].SetSubVector(pec_tdofs, 0.0);
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
      Vector on_sheets;
      full.GetSubVector(sheet_tdofs, on_sheets);
      a.SetSubVector(sheet_tdofs, on_sheets);
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

TEST_CASE("SpaceOperator assembly restricted to region and environment",
          "[substructure][Serial][Parallel]")
{
  // The driven operators assembled on the region and on the environment add up to the
  // full operators, and neither side couples to the other side's interior DOFs: lumped
  // ports on both sides, a lossy dielectric region, a conducting environment, a
  // second-order absorbing boundary on both sides and an impedance sheet crossing Γ.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  json config = {
      {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"LossTan", 0.01}},
          {{"Attributes", {2}}, {"Permittivity", 4.0}, {"Conductivity", 1.0}}}}}},
      {"Boundaries",
       {{"LumpedPort",
         {{{"Index", 1},
           {"R", 1.0},
           {"Attributes", {1}},
           {"Direction", "+Y"},
           {"Excitation", true}},
          {{"Index", 2}, {"R", 1.0}, {"Attributes", {2}}, {"Direction", "+Y"}}}},
        {"Absorbing", {{"Attributes", {3}}, {"Order", 2}}},
        {"Impedance", {{{"Attributes", {4}}, {"Rs", 1.0}, {"Ls", 2.0}, {"Cs", 0.5}}}}}},
      {"Solver",
       {{"Order", order},
        {"Device", "CPU"},
        {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 1.0}}}}}};
  IoData iodata(config, false);
  RegionDesign design;
  design.sheet = true;
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(2, 3, 4, design)));
  SpaceOperator space_op(iodata, mesh);
  auto &pmesh = mesh.back()->Get();
  const auto &fes = space_op.GetNDSpace().Get();

  // True DOFs of region (attribute 1) and environment (attribute 2) elements.
  mfem::Array<int> marks[2];
  for (int side = 0; side < 2; side++)
  {
    mfem::Vector m(fes.GetVSize()), mt(fes.GetTrueVSize());
    m = 0.0;
    mfem::Array<int> vdofs;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      if (pmesh.GetAttribute(e) == side + 1)
      {
        fes.GetElementVDofs(e, vdofs);
        for (int d : vdofs)
        {
          m(d >= 0 ? d : -1 - d) = 1.0;
        }
      }
    }
    fes.GetProlongationMatrix()->MultTranspose(m, mt);
    marks[side].SetSize(mt.Size());
    const double *t = mt.HostRead();
    for (int i = 0; i < mt.Size(); i++)
    {
      marks[side][i] = (t[i] > 0.0);
    }
  }

  const std::vector<int> region = {1}, environment = {2};
  const double omega = 1.3;  // nondimensional (the configuration is not nondimensionalized)
  auto assemble = [&](const std::vector<int> *domains)
  {
    std::optional<AssemblyDomains::Restriction> restriction;
    if (domains)
    {
      restriction.emplace(space_op.GetAssemblyDomains(), *domains);
    }
    std::vector<std::unique_ptr<ComplexOperator>> ops;
    ops.push_back(space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO));
    ops.push_back(space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO));
    ops.push_back(space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO));
    ops.push_back(
        space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO));
    return ops;
  };
  auto full = assemble(nullptr), reg = assemble(&region), env = assemble(&environment);

  // y = A x, for one part (real or imaginary) of an operator (zero if absent).
  auto apply = [](const ComplexOperator *A, bool imag, const Vector &x, Vector &y)
  {
    y = 0.0;
    const Operator *part = A ? (imag ? A->Imag() : A->Real()) : nullptr;
    if (part)
    {
      part->Mult(x, y);
    }
  };
  const int nt = fes.GetTrueVSize();
  Vector x(nt), y_full(nt), y_reg(nt), y_env(nt), x_eint(nt), x_rint(nt);
  x.Randomize(1 + Mpi::Rank(pmesh.GetComm()));
  for (int i = 0; i < nt; i++)
  {
    x_eint(i) = (marks[1][i] && !marks[0][i]) ? x(i) : 0.0;  // environment interior
    x_rint(i) = (marks[0][i] && !marks[1][i]) ? x(i) : 0.0;  // region interior
  }
  const char *names[4] = {"K", "C", "M", "A2"};
  for (int k = 0; k < 4; k++)
  {
    for (const bool imag : {false, true})
    {
      CAPTURE(names[k], imag);
      apply(full[k].get(), imag, x, y_full);
      apply(reg[k].get(), imag, x, y_reg);
      apply(env[k].get(), imag, x, y_env);
      const double n_full = linalg::Norml2(pmesh.GetComm(), y_full);
      y_reg += y_env;
      y_reg -= y_full;
      CHECK(linalg::Norml2(pmesh.GetComm(), y_reg) <= 1.0e-12 * (n_full + 1.0e-300));
      apply(reg[k].get(), imag, x_eint, y_reg);
      apply(env[k].get(), imag, x_rint, y_env);
      CHECK(linalg::Norml2(pmesh.GetComm(), y_reg) <= 1.0e-12 * (n_full + 1.0e-300));
      CHECK(linalg::Norml2(pmesh.GetComm(), y_env) <= 1.0e-12 * (n_full + 1.0e-300));
    }
  }
  // Each side carries its own part: nonzero region and environment stiffness and damping.
  for (int k : {0, 1})
  {
    apply(reg[k].get(), false, x, y_reg);
    apply(env[k].get(), false, x, y_env);
    CHECK(linalg::Norml2(pmesh.GetComm(), y_reg) > 0.0);
    CHECK(linalg::Norml2(pmesh.GetComm(), y_env) > 0.0);
  }

  // So do the preconditioners, at every level and on the auxiliary spaces.
  auto precond = [&](const std::vector<int> *domains)
  {
    std::optional<AssemblyDomains::Restriction> restriction;
    if (domains)
    {
      restriction.emplace(space_op.GetAssemblyDomains(), *domains);
    }
    return space_op.GetPreconditionerMatrix<ComplexOperator>(
        std::complex<double>(1.0), std::complex<double>(0.0, omega),
        std::complex<double>(-omega * omega), omega);
  };
  const auto P_full = precond(nullptr), P_reg = precond(&region),
             P_env = precond(&environment);
  auto mg = [](const ComplexOperator &P) -> const ComplexMultigridOperator &
  { return dynamic_cast<const ComplexMultigridOperator &>(P); };
  for (const bool aux : {false, true})
  {
    for (std::size_t l = 0; l < mg(*P_full).GetNumLevels(); l++)
    {
      auto level = [&](const ComplexOperator &P)
      {
        return aux ? &mg(P).GetAuxiliaryOperatorAtLevel(l) : &mg(P).GetOperatorAtLevel(l);
      };
      const int m = level(*P_full)->Width();
      Vector xl(m), yl_full(m), yl_reg(m), yl_env(m);
      xl.Randomize(3 + Mpi::Rank(pmesh.GetComm()));
      for (const bool imag : {false, true})
      {
        CAPTURE(aux, l, imag);
        apply(level(*P_full), imag, xl, yl_full);
        apply(level(*P_reg), imag, xl, yl_reg);
        apply(level(*P_env), imag, xl, yl_env);
        const double n_full = linalg::Norml2(pmesh.GetComm(), yl_full);
        yl_reg += yl_env;
        yl_reg -= yl_full;
        CHECK(linalg::Norml2(pmesh.GetComm(), yl_reg) <= 1.0e-12 * (n_full + 1.0e-300));
      }
    }
  }
  auto finest = [&](const ComplexOperator &P) { return &mg(P).GetFinestOperator(); };
  for (const bool imag : {false, true})
  {
    CAPTURE(imag);
    apply(finest(*P_reg), imag, x_eint, y_reg);
    apply(finest(*P_env), imag, x_rint, y_env);
    const double n = linalg::Norml2(pmesh.GetComm(), x);
    CHECK(linalg::Norml2(pmesh.GetComm(), y_reg) <= 1.0e-12 * n);
    CHECK(linalg::Norml2(pmesh.GetComm(), y_env) <= 1.0e-12 * n);
  }
  // The assembled matrix of the region couples no row of the environment's interior (which
  // keeps its zero diagonal entry): its elements are left out, not assembled with zero
  // coefficients, which a sparse direct solver would factor.
  const auto *part = dynamic_cast<const ParOperator *>(finest(*P_reg)->Real());
  REQUIRE(part);
  auto *A = (hypre_ParCSRMatrix *)part->ParallelAssemble();
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(A),
                        *offd = hypre_ParCSRMatrixOffd(A);
  int coupled = 0;
  for (int i = 0; i < nt; i++)
  {
    if (marks[1][i] && !marks[0][i])
    {
      for (HYPRE_Int k = hypre_CSRMatrixI(diag)[i]; k < hypre_CSRMatrixI(diag)[i + 1]; k++)
      {
        coupled += (hypre_CSRMatrixJ(diag)[k] != i);
      }
      coupled += hypre_CSRMatrixI(offd)[i + 1] - hypre_CSRMatrixI(offd)[i];
    }
  }
  Mpi::GlobalSum(1, &coupled, pmesh.GetComm());
  CHECK(coupled == 0);
}

TEST_CASE("Barycentric interpolation of a rational vector function",
          "[substructure][Serial]")
{
  // x(ω) = Σ_k r_k / (ω - p_k) + c_0 + c_1 ω on [1, 2], with two near-real poles in the
  // band and one outside it: greedy sampling at the least denominator converges with a few
  // samples, interpolates them, and the interpolant matches x between them. (Once the fit
  // is exact to round-off, further samples leave the weights undetermined.)
  constexpr int n = 40;
  const std::complex<double> p[3] = {{1.3, 1.0e-3}, {1.7, 2.0e-3}, {2.5, 0.1}};
  std::vector<std::complex<double>> r[3], c[2];
  std::uint32_t seed = 7;
  auto rnd = [&]()
  {
    seed = 1664525u * seed + 1013904223u;
    return static_cast<double>(seed) / 4294967296.0 - 0.5;
  };
  for (auto *v : {&r[0], &r[1], &r[2], &c[0], &c[1]})
  {
    for (int q = 0; q < n; q++)
    {
      v->emplace_back(rnd(), rnd());
    }
  }
  auto x = [&](double omega)
  {
    std::vector<std::complex<double>> y(n);
    for (int q = 0; q < n; q++)
    {
      y[q] = c[0][q] + omega * c[1][q];
      for (int k = 0; k < 3; k++)
      {
        y[q] += r[k][q] / (omega - p[k]);
      }
    }
    return y;
  };
  auto rel = [](const std::vector<std::complex<double>> &a,
                const std::vector<std::complex<double>> &b)
  {
    double d = 0.0, m = 0.0;
    for (std::size_t q = 0; q < a.size(); q++)
    {
      d += std::norm(a[q] - b[q]);
      m += std::norm(b[q]);
    }
    return std::sqrt(d / m);
  };
  BarycentricInterpolant fit;
  fit.AddSample(1.0, x(1.0));
  fit.AddSample(2.0, x(2.0));
  int memory = 0;
  while (fit.Samples().size() < 20 && memory < 2)
  {
    const double omega = fit.FindMaxError();
    const auto exact = x(omega);
    memory = (rel(fit.Evaluate(omega), exact) < 1.0e-8) ? memory + 1 : 0;
    fit.AddSample(omega, exact);
  }
  CAPTURE(fit.Samples().size());
  CHECK(memory == 2);
  CHECK(fit.Samples().size() <= 12);
  for (double omega : fit.Samples())
  {
    CHECK(rel(fit.Evaluate(omega), x(omega)) <= 1.0e-12);
  }
  double worst = 0.0;
  for (int i = 0; i <= 1000; i++)
  {
    const double omega = 1.0 + i / 1000.0;
    worst = std::max(worst, rel(fit.Evaluate(omega), x(omega)));
  }
  CAPTURE(worst);
  CHECK(worst <= 1.0e-7);
  const auto a =
      BarycentricInterpolant::Coefficients(fit.Samples(), fit.Weights(), fit.Samples()[2]);
  CHECK(a[2] == 1.0);
}

TEST_CASE("Passivity of a condensed environment", "[substructure][Serial]")
{
  // S = Re + i Im with Im = Q diag(d) Q^T (Q a rotation): the least eigenvalue of Im,
  // relative to ‖S‖_F, for any scaling of S, and the Cholesky test of passivity to a
  // tolerance, about that eigenvalue.
  const double c = std::cos(0.3), s = std::sin(0.3);
  auto lower = [&](double d0, double d1, double scale)
  {
    // Lower triangle by columns of the 2 x 2 S: (0,0), (1,0), (1,1).
    const double i00 = c * c * d0 + s * s * d1, i10 = c * s * (d0 - d1),
                 i11 = s * s * d0 + c * c * d1;
    return std::vector<std::complex<double>>{scale * std::complex<double>(2.0, i00),
                                             scale * std::complex<double>(0.5, i10),
                                             scale * std::complex<double>(1.0, i11)};
  };
  auto norm = [](const std::vector<std::complex<double>> &S)
  { return std::sqrt(std::norm(S[0]) + 2.0 * std::norm(S[1]) + std::norm(S[2])); };
  for (double scale : {1.0, 1.0e-6})
  {
    const auto passive = lower(0.3, 0.1, scale), active = lower(0.3, -0.2, scale);
    CHECK(DrivenSubstructureModel::Passivity(passive.data(), 2) ==
          Catch::Approx(0.1 * scale / norm(passive)).epsilon(1.0e-12));
    CHECK(DrivenSubstructureModel::Passivity(active.data(), 2) ==
          Catch::Approx(-0.2 * scale / norm(active)).epsilon(1.0e-12));
    const double violation = 0.2 * scale / norm(active);
    CHECK(DrivenSubstructureModel::Passive(passive.data(), 2, 1.0e-12));
    CHECK_FALSE(DrivenSubstructureModel::Passive(active.data(), 2, 0.5 * violation));
    CHECK(DrivenSubstructureModel::Passive(active.data(), 2, 2.0 * violation));
  }

  // A least eigenvalue repeated four times, whose whole cluster the bisection of LAPACK
  // returns: S = I + i diag(-0.5, -0.5, -0.5, -0.5, 1, 1).
  const int n = 6;
  std::vector<std::complex<double>> S;
  for (int j = 0; j < n; j++)
  {
    for (int i = j; i < n; i++)
    {
      S.emplace_back((i == j) ? 1.0 : 0.0, (i == j) ? ((j < 4) ? -0.5 : 1.0) : 0.0);
    }
  }
  const double fro = std::sqrt(4.0 * 1.25 + 2.0 * 2.0);
  CHECK(DrivenSubstructureModel::Passivity(S.data(), n) ==
        Catch::Approx(-0.5 / fro).epsilon(1.0e-12));
}

TEST_CASE("Saved driven excitations match across partitions", "[substructure][Serial]")
{
  // Source fingerprints of a lumped port of the CPW resonator grid example (in the plane
  // z = 0, along x), saved on 192 processes and recomputed on 8: the pairing with the first
  // fingerprint field (z, x, y) vanishes, up to rounding that differs between partitions.
  using Model = DrivenSubstructureModel;
  const std::vector<double> online = {1.0,
                                      0.0,
                                      -5.2511761922216104e-16,
                                      0.0,
                                      9.2812996599596332e-04,
                                      0.0,
                                      -1.9606511512972598e-03},
                            saved = {1.0,
                                     0.0,
                                     -4.8624843274134622e-16,
                                     0.0,
                                     9.2812996599596874e-04,
                                     0.0,
                                     -1.9606511512972719e-03},
                            other = {1.0,
                                     0.0,
                                     -1.2238455583072405e-15,
                                     0.0,
                                     9.2812996599578388e-04,
                                     0.0,
                                     1.9606511512968894e-03};  // the port at the other end
  REQUIRE(online.size() == static_cast<std::size_t>(Model::kSourceFp));
  CHECK(Model::SameSource(online.data(), saved.data()));
  CHECK_FALSE(Model::SameSource(online.data(), other.data()));
  auto no_env = online;
  no_env[0] = 0.0;
  CHECK_FALSE(Model::SameSource(no_env.data(), saved.data()));

  // A wave port's source follows its mode, reproducible to an eigensolver tolerance only:
  // a change at that level matches to a tolerance above it, a different source does not.
  auto perturbed = saved;
  for (int q = 1; q < Model::kSourceFp; q++)
  {
    perturbed[q] *= 1.0 + 1.0e-7 * q;
  }
  CHECK_FALSE(Model::SameSource(online.data(), perturbed.data()));
  CHECK(Model::SameSource(online.data(), perturbed.data(), 1.0e-4));
  CHECK_FALSE(Model::SameSource(online.data(), other.data(), 1.0e-4));

  // Environment fingerprints: counts, then quadratic forms with a vanishing imaginary part.
  const std::vector<double> env = {100.0, 20.0, 3.0, 0.0, -2.0, 1.0e-3, 5.0, 0.0};
  auto env2 = env;
  env2[2] *= 1.0 + 1.0e-13;
  CHECK(Model::SameEnvironment(env, env2));
  env2[6] *= 1.0 + 1.0e-6;
  CHECK_FALSE(Model::SameEnvironment(env, env2));
  env2 = env;
  env2[1] += 1.0;
  CHECK_FALSE(Model::SameEnvironment(env, env2));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "A cut or corrupted driven model header is rejected",
                 "[substructure][Serial]")
{
  // The header is read part by part, each only if its record holds it: a header cut short,
  // one whose sizes exceed its record (also a negative one), one without frequencies, and a
  // model without the record of its environment, are incomplete; another format version is
  // rejected, and a file of another kind is not a model.
  DrivenSubstructureModel m;
  m.nG = 2;
  m.sig_w = 1;
  m.signatures = {0.25, 0.75};
  m.env_fp = {100.0, 20.0, 3.0, 0.0};
  m.omega = {1.0};
  m.fp_centers = {{{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}}};
  m.converged = false;
  m.error = 0.5;
  m.tol = 1.0e-8;
  m.passivity = -1.0e-3;
  const auto path = (temp_dir / "driven-model-cut.bin").string();
  // The position of the payload of the last record with a tag in the file.
  using substructuring::ModelFile;
  auto payload = [&path](std::uint32_t tag)
  { return ModelFile::Read(path, Mpi::World()).Find(tag)->pos; };
  const auto header = ModelFile::Tag("DHED"), signatures = ModelFile::Tag("SIGS");
  m.WriteHeader(path);
  DrivenSubstructureModel saved;
  saved.ReadHeader(path, Mpi::World());
  CHECK(saved.signatures == m.signatures);
  CHECK(saved.fp_centers == m.fp_centers);
  CHECK_FALSE(saved.converged);
  CHECK(saved.error == m.error);
  CHECK(saved.tol == m.tol);
  CHECK(saved.passivity == m.passivity);
  using Catch::Matchers::ContainsSubstring;
  std::filesystem::resize_file(path, std::filesystem::file_size(path) / 2);
  CHECK_THROWS_WITH(saved.ReadHeader(path, Mpi::World()),
                    ContainsSubstring("is incomplete"));
  // The interface size (the first entry of the interface record) beyond the record or
  // negative, and more frequencies (head entry 4) than the record holds, or none (an
  // adaptive offline run stopped before its first sample).
  for (const auto [tag, entry, value] : {std::tuple{signatures, 0, 1000},
                                         {signatures, 0, -5},
                                         {header, 4, 7},
                                         {header, 4, 0}})
  {
    m.WriteHeader(path);
    const auto pos = payload(tag);
    {
      std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
      f.seekp(pos + entry * static_cast<long long>(sizeof(int)));
      f.write(reinterpret_cast<const char *>(&value), sizeof(int));
    }
    CAPTURE(entry, value);
    CHECK_THROWS_WITH(saved.ReadHeader(path, Mpi::World()),
                      ContainsSubstring("is incomplete"));
  }
  m.WriteHeader(path);
  {
    using substructuring::ModelFile;
    const auto file = ModelFile::Read(path, Mpi::World());
    const auto *e = file.Find(ModelFile::Tag("ENVR"));
    REQUIRE(e);
    std::filesystem::resize_file(path, e->pos - 12);
  }
  CHECK_THROWS_WITH(saved.ReadHeader(path, Mpi::World()),
                    ContainsSubstring("is incomplete"));
  m.WriteHeader(path);
  {
    const auto pos = payload(header);
    std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
    const int version = DrivenSubstructureModel::kVersion + 1;
    f.seekp(pos);
    f.write(reinterpret_cast<const char *>(&version), sizeof(int));
  }
  CHECK_THROWS_WITH(saved.ReadHeader(path, Mpi::World()),
                    ContainsSubstring("written by another version"));
  std::ofstream(path, std::ios::binary) << "not a model";
  CHECK_THROWS_WITH(saved.ReadHeader(path, Mpi::World()),
                    ContainsSubstring("Cannot read the driven substructuring model"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "A saved driven model continues from its complete records",
                 "[substructure][Serial]")
{
  // A model of two interface unknowns, an excitation and a port at three frequencies,
  // interrupted while writing its third record, then cut to its first record by a restart
  // at the second frequency.
  using Model = DrivenSubstructureModel;
  Model m;
  m.nG = 2;
  m.sig_w = 1;
  m.signatures = {0.25, 0.75};
  m.env_fp = {100.0, 20.0, 3.0, 0.0, -2.0, 1.0e-3};
  m.excitations = {1};
  m.exc_fp = {1.0, 0.0, 0.5, 0.0, 0.25, 0.0, -0.125};
  m.ports = {1};
  m.region_bdr_attributes = {5, 7};
  m.environment.attributes = {2};
  m.environment.bdr_attributes = {3, 4};
  m.environment.moments = {{4.0, 1.0, 0.5, 0.25, 0.75, 0.3}};
  m.environment.region_moments = {{2.0, 0.5, 0.25, 0.1, 0.4, 0.1}};
  m.environment.bdr_moments = {{2.0, 1.0, 0.5, 0.5, 1.0, 0.25},
                               {1.0, 0.5, 0.0, 0.25, 0.5, 0.1}};
  m.environment.config_hash = 11;
  m.environment.excitations_hash = 13;
  m.omega = {1.0, 2.0, 3.0};
  m.length_scale = 1.0e-6;
  auto record = [](double s)
  {
    Model::Record r;
    r.S = {s, 2.0 * s, 2.0 * s, 3.0 * s};
    r.g = {s, -s};
    r.h = {std::complex<double>(0.0, s), s};
    r.c = {s};
    return r;
  };
  const auto path = (temp_dir / "driven-model-restart.bin").string();
  m.WriteHeader(path);
  m.AppendRecord(path, record(1.0));
  m.AppendRecord(path, record(2.0));
  {
    std::ofstream f(path, std::ios::binary | std::ios::app);
    const double part[3] = {1.0, 2.0, 3.0};
    f.write(reinterpret_cast<const char *>(part), sizeof(part));
  }
  CHECK(m.NumRecords(path) == 2);
  Model saved;
  saved.ReadHeader(path, Mpi::World());
  CHECK(saved.region_bdr_attributes == m.region_bdr_attributes);
  CHECK(saved.environment.bdr_moments == m.environment.bdr_moments);
  CHECK(saved.environment.config_hash == m.environment.config_hash);
  CHECK(m.SameHeader(saved, 1.0e-10));
  m.TruncateRecords(path, 1);
  CHECK(m.NumRecords(path) == 1);
  CHECK(m.ReadFlatRecord(path, 0) == m.Flatten(record(1.0)));
  m.AppendRecord(path, record(4.0));
  CHECK(m.NumRecords(path) == 2);
  CHECK(m.ReadFlatRecord(path, 1) == m.Flatten(record(4.0)));

  // The model of another sweep does not continue: other frequencies, interface, boundaries,
  // ports, environment (beyond the rounding of its moments) or sources (beyond their
  // tolerance).
  auto other = m;
  other.region_bdr_attributes = {5};
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.ports = {2};
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.environment.config_hash++;
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.environment.region_moments[0][1] *= 1.0 + 1.0e-6;
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other.environment.region_moments[0][1] =
      m.environment.region_moments[0][1] * (1.0 + 1.0e-13);
  CHECK(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.omega[2] = 3.5;
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.signatures[1] = 0.5;
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  other = m;
  other.exc_fp[6] *= 1.0 + 1.0e-6;
  CHECK_FALSE(other.SameHeader(saved, 1.0e-10));
  CHECK(other.SameHeader(saved, 1.0e-4));
}

#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)
TEST_CASE("DirectSchurSolver condenses exactly", "[substructure][Serial][Parallel]")
{
  // The complex symmetric shifted Laplacian of the MumpsSchurSolver tests, bordered by an
  // unknown on the last rank (as by a wave port's modal term), with the Schur variables on
  // a face of the cube and every other entry passed by another rank than that of its row
  // (some split in two): the Schur complement and a condensed solve against dense ones,
  // also after a refactorization with new values, by STRUMPACK and by SuperLU_DIST.
  std::vector<LinearSolver> types;
#if defined(MFEM_USE_STRUMPACK)
  types.push_back(LinearSolver::STRUMPACK);
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back(LinearSolver::SUPERLU);
#endif
  const LinearSolver type = GENERATE_COPY(from_range(types));
  CAPTURE(static_cast<int>(type));
  auto pmesh = MakeSplitCube(4);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fes(pmesh.get(), &fec);
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.AddDomainIntegrator(new mfem::MassIntegrator);
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize(), nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const bool last = (rank == nranks - 1);
  const HYPRE_BigInt n_dof = A->GetGlobalNumRows(), off = fes.GetMyTDofOffset();
  const int n_glob = static_cast<int>(n_dof) + 1;
  mfem::Array<int> bdr(pmesh->bdr_attributes.Max()), face;
  bdr = 0;
  bdr[0] = 1;
  fes.GetEssentialTrueDofs(bdr, face);
  std::vector<HYPRE_BigInt> mine(face.Size()), schur_vars;
  for (int i = 0; i < face.Size(); i++)
  {
    mine[i] = off + face[i];
  }
  std::vector<int> cnt(nranks), disp(nranks, 0);
  {
    int n = face.Size();
    MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    schur_vars.resize(disp.back() + cnt.back());
    MPI_Allgatherv(mine.data(), n, HYPRE_MPI_BIG_INT, schur_vars.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
  }
  REQUIRE(!schur_vars.empty());

  // The entries: c A by its lower triangle, and the border unknown (row n_dof) with entries
  // b against the first local DOFs of each rank and -2 on its diagonal.
  const std::complex<double> c(1.0, 0.1), b(0.5, 0.2);
  DirectSchurSolver::Coo coo;
  {
    std::vector<std::vector<double>> vals;
    std::vector<int> row_ptr;
    LowerTrianglePattern({A.get()}, coo.irn, coo.jcn, vals, row_ptr);
    for (double v : vals[0])
    {
      coo.val.push_back(c * v);
    }
  }
  for (int i = 0; i < std::min(3, nt); i++)
  {
    coo.irn.push_back(n_glob);
    coo.jcn.push_back(static_cast<int>(off) + i + 1);
    coo.val.push_back(b);
  }
  if (last)
  {
    coo.irn.push_back(n_glob);
    coo.jcn.push_back(n_glob);
    coo.val.push_back(-2.0);
  }
  {
    DirectSchurSolver::Coo keep, give;
    auto push = [](DirectSchurSolver::Coo &x, int i, int j, std::complex<double> v)
    {
      x.irn.push_back(i);
      x.jcn.push_back(j);
      x.val.push_back(v);
    };
    for (std::size_t e = 0; e < coo.irn.size(); e++)
    {
      const bool split = (e % 4 == 3);
      push((e % 2) ? give : keep, coo.irn[e], coo.jcn[e], (split ? 0.5 : 1.0) * coo.val[e]);
      if (split)
      {
        push(keep, coo.irn[e], coo.jcn[e], 0.5 * coo.val[e]);
      }
    }
    const int dst = (rank + 1) % nranks, src = (rank + nranks - 1) % nranks;
    int ns = static_cast<int>(give.irn.size()), nr = 0;
    MPI_Sendrecv(&ns, 1, MPI_INT, dst, 0, &nr, 1, MPI_INT, src, 0, comm, MPI_STATUS_IGNORE);
    std::vector<int> ri(nr), rj(nr);
    std::vector<std::complex<double>> rv(nr);
    MPI_Sendrecv(give.irn.data(), ns, MPI_INT, dst, 1, ri.data(), nr, MPI_INT, src, 1, comm,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(give.jcn.data(), ns, MPI_INT, dst, 2, rj.data(), nr, MPI_INT, src, 2, comm,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(give.val.data(), ns, MPI_C_DOUBLE_COMPLEX, dst, 3, rv.data(), nr,
                 MPI_C_DOUBLE_COMPLEX, src, 3, comm, MPI_STATUS_IGNORE);
    for (int e = 0; e < nr; e++)
    {
      push(keep, ri[e], rj[e], rv[e]);
    }
    coo = std::move(keep);
  }

  // The dense matrix on rank 0, from all the entries.
  Eigen::MatrixXcd M;
  {
    int ne = static_cast<int>(coo.irn.size());
    std::vector<int> ecnt(nranks), edisp(nranks, 0);
    MPI_Gather(&ne, 1, MPI_INT, ecnt.data(), 1, MPI_INT, 0, comm);
    for (int r = 1; r < nranks; r++)
    {
      edisp[r] = edisp[r - 1] + ecnt[r - 1];
    }
    const int ntot = edisp.back() + ecnt.back();
    std::vector<int> gi(rank == 0 ? ntot : 0), gj(gi.size());
    std::vector<std::complex<double>> gv(gi.size());
    MPI_Gatherv(coo.irn.data(), ne, MPI_INT, gi.data(), ecnt.data(), edisp.data(), MPI_INT,
                0, comm);
    MPI_Gatherv(coo.jcn.data(), ne, MPI_INT, gj.data(), ecnt.data(), edisp.data(), MPI_INT,
                0, comm);
    MPI_Gatherv(coo.val.data(), ne, MPI_C_DOUBLE_COMPLEX, gv.data(), ecnt.data(),
                edisp.data(), MPI_C_DOUBLE_COMPLEX, 0, comm);
    if (rank == 0)
    {
      M = Eigen::MatrixXcd::Zero(n_glob, n_glob);
      for (int e = 0; e < ntot; e++)
      {
        M(gi[e] - 1, gj[e] - 1) += gv[e];
        if (gi[e] != gj[e])
        {
          M(gj[e] - 1, gi[e] - 1) += gv[e];
        }
      }
    }
  }
  ComplexVector x(nt);
  x.Real().Randomize(1 + rank);
  x.Imag().Randomize(2 + rank);
  // A vector of the true DOFs on rank 0.
  auto gather = [&](const ComplexVector &v)
  {
    std::vector<int> vcnt(nranks), vdisp(nranks, 0);
    int n = nt;
    MPI_Gather(&n, 1, MPI_INT, vcnt.data(), 1, MPI_INT, 0, comm);
    for (int r = 1; r < nranks; r++)
    {
      vdisp[r] = vdisp[r - 1] + vcnt[r - 1];
    }
    std::vector<std::complex<double>> loc(nt), glob(rank == 0 ? n_dof : 0);
    for (int i = 0; i < nt; i++)
    {
      loc[i] = {v.Real()[i], v.Imag()[i]};
    }
    MPI_Gatherv(loc.data(), nt, MPI_C_DOUBLE_COMPLEX, glob.data(), vcnt.data(),
                vdisp.data(), MPI_C_DOUBLE_COMPLEX, 0, comm);
    return glob;
  };
  std::vector<int> sch, in;
  Eigen::MatrixXcd S_ref;
  Eigen::VectorXcd u_ref, rhs_in;
  if (rank == 0)
  {
    std::vector<char> is_sch(n_glob, 0);
    for (auto g : schur_vars)
    {
      is_sch[g] = 1;
    }
    for (auto g : schur_vars)
    {
      sch.push_back(static_cast<int>(g));
    }
    for (int g = 0; g < n_glob; g++)
    {
      if (!is_sch[g])
      {
        in.push_back(g);
      }
    }
    const Eigen::MatrixXcd Aii = M(in, in), Ais = M(in, sch), Ass = M(sch, sch);
    S_ref = Ass - Ais.transpose() * Aii.partialPivLu().solve(Ais);
  }
  const auto xg = gather(x);
  if (rank == 0)
  {
    Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n_glob);
    for (HYPRE_BigInt g = 0; g < n_dof; g++)
    {
      rhs(g) = xg[g];
    }
    u_ref = M.partialPivLu().solve(rhs);
    rhs_in = rhs(in);
  }

  std::vector<int> rows(nt + (last ? 1 : 0), -1);
  for (int i = 0; i < nt; i++)
  {
    rows[i] = i;
  }
  DirectSchurSolver solver(comm, n_glob, static_cast<int>(rows.size()), std::move(coo),
                           schur_vars, rows, type, SymbolicFactorization::AMD);
  auto dist = [](const auto &a, const auto &b)
  { return (a - b).cwiseAbs().maxCoeff() / b.cwiseAbs().maxCoeff(); };
  for (int k = 0; k < 2; k++)
  {
    const double s = (k == 0) ? 1.0 : 2.0;
    if (k == 1)
    {
      for (auto &v : solver.Values())
      {
        v *= s;
      }
      solver.Refactor();
    }
    const auto S = solver.ReleaseSchur();
    std::vector<std::complex<double>> red;
    solver.Reduce({&x}, red);
    std::vector<std::complex<double>> uS;
    if (rank == 0)
    {
      const int nS = static_cast<int>(sch.size());
      const Eigen::Map<const Eigen::MatrixXcd> Sm(S.data(), nS, nS);
      CHECK(dist(Sm, s * S_ref) <= 1.0e-10);
      const Eigen::VectorXcd u =
          Sm.partialPivLu().solve(Eigen::Map<const Eigen::VectorXcd>(red.data(), nS));
      uS.assign(u.data(), u.data() + nS);
    }
    ComplexVector y(nt);
    solver.Expand(uS, {&y});
    const auto yg = gather(y);
    if (rank == 0)
    {
      const Eigen::Map<const Eigen::VectorXcd> ym(yg.data(), static_cast<int>(n_dof));
      CHECK(dist(ym, u_ref.head(n_dof) / s) <= 1.0e-10);
    }
    solver.SolveInternal({&x}, {&y});
    const auto zg = gather(y);
    if (rank == 0)
    {
      Eigen::VectorXcd z = Eigen::VectorXcd::Zero(n_glob);
      const Eigen::MatrixXcd Aii = M(in, in);
      const Eigen::VectorXcd zi = Aii.partialPivLu().solve(rhs_in);
      z(in) = zi;
      const Eigen::Map<const Eigen::VectorXcd> zm(zg.data(), static_cast<int>(n_dof));
      CHECK(dist(zm, z.head(n_dof) / s) <= 1.0e-10);
    }
  }

  // The rows of S (of the values scaled by 2) on the ranks holding its variables, from
  // another factorization.
  solver.Refactor();
  const int nS = static_cast<int>(schur_vars.size()), m_loc = static_cast<int>(rows.size());
  if (rank != 0)
  {
    S_ref.resize(nS, nS);
  }
  MPI_Bcast(S_ref.data(), nS * nS, MPI_C_DOUBLE_COMPLEX, 0, comm);
  HYPRE_BigInt n = m_loc, first = 0;
  MPI_Exscan(&n, &first, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
  first = (rank == 0) ? 0 : first;
  std::vector<int> held;
  for (int p = 0; p < nS; p++)
  {
    if (schur_vars[p] >= first && schur_vars[p] < first + m_loc)
    {
      held.push_back(p);
    }
  }
  const auto R = solver.ReleaseSchurRows();
  CHECK(R.size() == held.size() * nS);
  double err = 0.0;
  for (std::size_t r = 0; r < held.size() && R.size() == held.size() * nS; r++)
  {
    for (int col = 0; col < nS; col++)
    {
      err = std::max(err, std::abs(R[r * nS + col] - 2.0 * S_ref(held[r], col)));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, comm);
  CHECK(err <= 1.0e-10 * 2.0 * S_ref.cwiseAbs().maxCoeff());
}
#endif

// The sweeps run in process condense the environment offline, by MUMPS or STRUMPACK.
#if defined(MFEM_USE_MUMPS) || defined(MFEM_USE_STRUMPACK)

namespace
{

// The model of the offline driven runs below, saved in their output directory.
std::filesystem::path DrivenModel(const std::filesystem::path &dir)
{
  return dir / "postpro_Offline" / "environment.model";
}

// A driven sweep of the substructuring_driven mesh (Γ the plane x = 0.5, the sheet 8 at
// z = 0.5 crossing it, lumped port 1 in the region and 2 in the environment, each its own
// excitation), run in process. Collective.
void RunDrivenSheetSquare(const std::filesystem::path &dir, const std::string &mode,
                          const std::function<void(json &)> &edit = {})
{
  const auto mesh = std::filesystem::path(PALACE_TEST_DATA_DIR) / "regression" / "input" /
                    "substructuring_magnetostatic" / "mesh" / "sheet_square.msh";
  json config = {
      {"Problem", {{"Type", "Driven"}, {"Output", (dir / ("postpro_" + mode)).string()}}},
      {"Model", {{"Mesh", mesh.string()}, {"L0", 1.0e-6}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 4.0}},
          {{"Attributes", {2}}, {"Permittivity", 2.0}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {2, 3, 4, 5, 6, 7}}}},
        {"LumpedPort",
         {{{"Index", 1},
           {"R", 50.0},
           {"Attributes", {9}},
           {"Direction", "+X"},
           {"Excitation", 1}},
          {{"Index", 2},
           {"R", 50.0},
           {"Attributes", {10}},
           {"Direction", "+X"},
           {"Excitation", 2}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Device", "CPU"},
        {"Driven", {{"MinFreq", 1000.0}, {"MaxFreq", 2000.0}, {"FreqStep", 1000.0}}},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  if (mode == "Online")
  {
    config["Solver"]["Substructuring"]["ModelFile"] = DrivenModel(dir).string();
  }
  if (edit)
  {
    edit(config);
  }
  IoData iodata(config, false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  palace::Run(iodata, comm, 1);
}

// The header (trimmed) and the rows of a CSV file written by Palace.
std::pair<std::vector<std::string>, std::vector<std::vector<double>>>
ReadCsv(const std::filesystem::path &path)
{
  std::ifstream f(path);
  REQUIRE(f.good());
  auto split = [](const std::string &line)
  {
    std::vector<std::string> cells;
    std::stringstream ss(line);
    for (std::string cell; std::getline(ss, cell, ',');)
    {
      cells.push_back(cell);
    }
    return cells;
  };
  std::string line;
  std::getline(f, line);
  std::vector<std::string> header;
  for (auto cell : split(line))
  {
    cell.erase(0, cell.find_first_not_of(' '));
    cell.erase(cell.find_last_not_of(' ') + 1);
    header.push_back(cell);
  }
  std::vector<std::vector<double>> rows;
  while (std::getline(f, line))
  {
    std::vector<double> row;
    for (const auto &cell : split(line))
    {
      row.push_back(std::stod(cell));
    }
    rows.push_back(row);
  }
  return {header, rows};
}

// The configuration of a driven sweep of the box of MakeInteriorBoundaryBox (Γ the interior
// boundary x = 1 between the region x < 1 and the environment, lumped port 1 on the face
// x = 0 and the face x = 2 absorbing), whose mesh it writes. Collective.
json DrivenBoxConfig(const std::filesystem::path &dir, const std::string &mode,
                     const std::function<void(json &)> &edit = {})
{
  const auto mesh = dir / "box.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    MakeInteriorBoundaryBox().Save(mesh.string());
  }
  Mpi::Barrier(Mpi::World());
  json config = {
      {"Problem", {{"Type", "Driven"}, {"Output", (dir / ("postpro_" + mode)).string()}}},
      {"Model", {{"Mesh", mesh.string()}, {"L0", 1.0e-3}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 4.0}},
          {{"Attributes", {2}}, {"Permittivity", 2.0}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {3}}}},
        {"Absorbing", {{"Attributes", {2}}, {"Order", 1}}},
        {"LumpedPort",
         {{{"Index", 1},
           {"R", 50.0},
           {"Attributes", {1}},
           {"Direction", "+Y"},
           {"Excitation", 1}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Device", "CPU"},
        {"Driven", {{"MinFreq", 20.0}, {"MaxFreq", 20.0}, {"FreqStep", 1.0}}},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  if (mode == "Online")
  {
    config["Solver"]["Substructuring"]["ModelFile"] = DrivenModel(dir).string();
  }
  if (edit)
  {
    edit(config);
  }
  return config;
}

// The sweep of DrivenBoxConfig, run in process. Collective.
void RunDrivenBox(const std::filesystem::path &dir, const std::string &mode,
                  const std::function<void(json &)> &edit = {})
{
  IoData iodata(DrivenBoxConfig(dir, mode, edit), false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  palace::Run(iodata, comm, 1);
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring region solver inverts the region's system",
                 "[substructure][Serial][Parallel]")
{
  // Online, the RegionSolver of the region of the box (a lumped port, so that the
  // preconditioner matrix is the system's but for S_E) solves the system of the sweep, with
  // S_E, to rounding, also after a refactorization: with the dense matrices on Γ on the
  // grid of their size, and on a grid of all the processes, and with each sparse direct
  // solver built for the Schur complement (MUMPS's distributed, STRUMPACK's or
  // SuperLU_DIST's on rank 0), and the default one. Also with second-order absorbing walls
  // crossing Γ instead of the PEC ones: their term keeps its sign in the region's
  // preconditioner matrix (it is flipped on the coarse level of a multigrid preconditioner,
  // applied to a real-valued approximation).
  const bool absorbing = GENERATE(false, true);
  auto walls = [=](json &config)
  {
    if (absorbing)
    {
      config["Boundaries"].erase("PEC");
      config["Boundaries"]["Absorbing"] = {{"Attributes", {2, 3}}, {"Order", 2}};
    }
  };
  RunDrivenBox(temp_dir, "Offline", walls);
  MPI_Comm comm = Mpi::World();
  const int grid_procs = GENERATE_COPY(0, Mpi::Size(comm));
  std::vector<std::string> types = {"Default"};
#if defined(MFEM_USE_MUMPS)
  types.emplace_back("MUMPS");
#endif
#if defined(MFEM_USE_STRUMPACK)
  types.emplace_back("STRUMPACK");
#endif
#if defined(MFEM_USE_SUPERLU)
  types.emplace_back("SuperLU");
#endif
  const std::string type = GENERATE_COPY(Catch::Generators::from_range(types));
  IoData iodata(DrivenBoxConfig(temp_dir, "Online",
                                [&](json &config)
                                {
                                  walls(config);
                                  config["Solver"]["Linear"]["Type"] = type;
                                }),
                false);
  MakeOutputFolder(iodata, comm);
  DrivenSolver driven(iodata, Mpi::Root(comm), Mpi::Size(comm), 1, nullptr);
  auto smesh = mesh::Load(iodata, comm);
  driven.Preprocess(iodata, smesh, comm);
  std::vector<std::unique_ptr<mfem::ParMesh>> pmesh;
  pmesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
  mesh::RefineMesh(iodata, pmesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(pmesh[0])));
  SpaceOperator space_op(iodata, mesh);
  EnvironmentOperator env(iodata, space_op, grid_procs);
  space_op.SetEnvironment(env);

  const double omega = iodata.solver.driven.sample_f[0];
  const std::complex<double> a0(1.0, 0.0), a1(0.0, omega), a2(-omega * omega, 0.0);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  auto A = space_op.GetSystemMatrix(a0, a1, a2, K.get(), C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, omega);
  RegionSolver solver(iodata, space_op, env);
  const int n = A->Height();
  ComplexVector x(n), y(n), r(n);
  x.UseDevice(true);
  y.UseDevice(true);
  r.UseDevice(true);
  x.Real().Randomize(1);
  x.Imag().Randomize(2);
  for (int factorization = 0; factorization < 2; factorization++)
  {
    solver.SetOperator(dynamic_cast<const BaseMultigridOperator<ComplexOperator> &>(*P)
                           .GetFinestOperator());
    solver.Mult(x, y);
    A->Mult(y, r);
    linalg::AXPY(-1.0, x, r);
    CHECK(linalg::Norml2(comm, r) <= 1.0e-10 * linalg::Norml2(comm, x));
  }
}

namespace
{

// A RomOperator which exposes the truncation of its basis, and the numbers of rows of the
// pairings of its basis with the sources of the environment, one per excitation.
class RegionRomOperator : public RomOperator
{
public:
  using RomOperator::RomOperator;
  using RomOperator::TruncatePROM;

  std::vector<Eigen::Index> EnvironmentPairingRows() const
  {
    std::vector<Eigen::Index> rows;
    for (const auto &pairings : env_RHSr)
    {
      rows.push_back(pairings.rows());
    }
    return rows;
  }
};

// A driven PostOperator which exposes its last measurement.
class RegionPostOperator : public PostOperator<ProblemType::DRIVEN>
{
public:
  using PostOperator<ProblemType::DRIVEN>::PostOperator;
  const Measurement &GetMeasurement() const { return measurement_cache; }
};

// The box of DrivenBoxConfig with a lumped port in the environment (excitation 2), swept
// adaptively over 10 to 30 GHz.
void EnvironmentPort(json &config)
{
  config["Boundaries"]["LumpedPort"].push_back({{"Index", 2},
                                                {"R", 50.0},
                                                {"Attributes", {5}},
                                                {"Direction", "+Y"},
                                                {"Excitation", 2}});
  config["Solver"]["Driven"] = {
      {"MinFreq", 10.0}, {"MaxFreq", 30.0}, {"FreqStep", 10.0}, {"AdaptiveTol", 1.0e-6}};
}

// The region of an online configuration with its environment, from the saved model.
// Collective.
struct OnlineRegion
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> meshes;
  std::unique_ptr<SpaceOperator> space_op;
  std::unique_ptr<EnvironmentOperator> env;

  explicit OnlineRegion(const json &config) : iodata(config, false)
  {
    MPI_Comm comm = Mpi::World();
    MakeOutputFolder(iodata, comm);
    DrivenSolver driven(iodata, Mpi::Root(comm), Mpi::Size(comm), 1, nullptr);
    auto smesh = mesh::Load(iodata, comm);
    driven.Preprocess(iodata, smesh, comm);
    std::vector<std::unique_ptr<mfem::ParMesh>> pmesh;
    pmesh.push_back(mesh::Partition(iodata, std::move(smesh), comm));
    mesh::RefineMesh(iodata, pmesh);
    meshes.push_back(std::make_unique<Mesh>(std::move(pmesh[0])));
    space_op = std::make_unique<SpaceOperator>(iodata, meshes);
    env = std::make_unique<EnvironmentOperator>(iodata, *space_op);
    space_op->SetEnvironment(*env);
  }
};

// A random field, zero on the Dirichlet DOFs.
ComplexVector RandomField(SpaceOperator &space_op, int seed)
{
  ComplexVector u(space_op.GetNDSpace().GetTrueVSize());
  u.UseDevice(true);
  linalg::SetRandom(space_op.GetComm(), u.Real(), seed);
  linalg::SetRandom(space_op.GetComm(), u.Imag(), seed + 1);
  linalg::SetSubVector(u, space_op.GetNDDbcTDofLists().back(), 0.0);
  return u;
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring PROM truncation also truncates the environment's "
                 "sources",
                 "[substructure][Serial][Parallel]")
{
  // TruncatePROM(n), with which AddSamplePROM removes a sample to add it again refined,
  // gives back the PROM of the first n basis vectors also online: with a lumped port in the
  // environment of a model saved by an adaptive sweep, the pairings of the basis with its
  // sources lose the rows of the removed vectors, and the reduced solution is that before
  // the sample was added.
  RunDrivenBox(temp_dir, "Offline", EnvironmentPort);
  OnlineRegion region(DrivenBoxConfig(temp_dir, "Online", EnvironmentPort));
  RegionRomOperator prom_op(region.iodata, *region.space_op, 4);
  auto pairing_rows = [&](Eigen::Index n)
  {
    const auto rows = prom_op.EnvironmentPairingRows();
    CAPTURE(rows, n);
    REQUIRE(std::count(rows.begin(), rows.end(), n) > 0);
    for (const auto r : rows)
    {
      REQUIRE((r == 0 || r == n));
    }
  };
  const double omega = region.iodata.solver.driven.sample_f[0];
  prom_op.UpdatePROM(RandomField(*region.space_op, 71), "s1");
  const auto n = static_cast<Eigen::Index>(prom_op.GetReducedDimension());
  pairing_rows(n);
  prom_op.SolvePROM(omega);
  const Eigen::VectorXcd x_n = prom_op.GetReducedSolution(2);
  REQUIRE(x_n.norm() > 0.0);

  prom_op.UpdatePROM(RandomField(*region.space_op, 73), "s2");
  REQUIRE(static_cast<Eigen::Index>(prom_op.GetReducedDimension()) == n + 2);
  pairing_rows(n + 2);
  prom_op.TruncatePROM(n);
  REQUIRE(static_cast<Eigen::Index>(prom_op.GetReducedDimension()) == n);
  pairing_rows(n);
  prom_op.SolvePROM(omega);
  CHECK((prom_op.GetReducedSolution(2) - x_n).norm() <= 1.0e-12 * x_n.norm());
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring reduced postprocessing matches the full evaluation",
                 "[substructure][Serial][Parallel]")
{
  // Online, the reduced postprocessing of the adaptive sweep takes the values of the
  // environment's ports (a lumped port in it, and with a wave port at its end, above the
  // cutoff of the box, an inactive wave port) from the values on Γ of the reduced
  // solution: they are those of the expanded solution for both excitations, and with them
  // the reduced measurements, the first of which the full evaluation checks, are the full
  // ones. An online adaptive sweep then postprocesses every frequency in reduced form.
  const bool wave = GENERATE(false, true);
  CAPTURE(wave);
  auto edit = [wave](json &config)
  {
    EnvironmentPort(config);
    if (wave)
    {
      config["Boundaries"].erase("Absorbing");
      config["Boundaries"]["WavePort"] = {
          {{"Index", 3}, {"Attributes", {2}}, {"Mode", 1}, {"Offset", 0.0}}};
      config["Solver"]["Driven"]["MinFreq"] = 150.0;
      config["Solver"]["Driven"]["MaxFreq"] = 250.0;
      config["Solver"]["Driven"]["FreqStep"] = 50.0;
    }
  };
  RunDrivenBox(temp_dir, "Offline", edit);
  {
    OnlineRegion region(DrivenBoxConfig(temp_dir, "Online", edit));
    auto &space_op = *region.space_op;
    auto &env = *region.env;
    RomOperator prom_op(region.iodata, space_op, 4);
    RegionPostOperator post_op(region.iodata, space_op);
    prom_op.UpdatePROM(RandomField(space_op, 81), "s1");
    post_op.ConfigureReducedPostprocessing(prom_op,
                                           region.iodata.solver.driven.adaptive_tol);
    REQUIRE(post_op.HasReducedPostprocessing());

    // A reduced solution y, its expansion E = V y, and its values on Γ from those of the
    // basis vectors.
    const auto &basis = prom_op.GetBasis();
    Eigen::VectorXcd y(static_cast<Eigen::Index>(basis.size()));
    for (Eigen::Index i = 0; i < y.size(); i++)
    {
      y(i) = std::complex<double>(1.0 / (i + 1.0), -0.25 / (i + 1.0));
    }
    ComplexVector E(space_op.GetNDSpace().GetTrueVSize());
    E.UseDevice(true);
    E = 0.0;
    std::vector<std::complex<double>> e_gamma;
    for (std::size_t i = 0; i < basis.size(); i++)
    {
      const auto yi = y(static_cast<Eigen::Index>(i));
      linalg::AXPY(yi.real(), basis[i], E.Real());
      linalg::AXPY(yi.imag(), basis[i], E.Imag());
      const auto v_gamma = env.LocalInterfaceValues(basis[i]);
      e_gamma.resize(v_gamma.size(), 0.0);
      for (std::size_t a = 0; a < v_gamma.size(); a++)
      {
        e_gamma[a] += yi * v_gamma[a];
      }
    }
    const double omega = region.iodata.solver.driven.sample_f[0];
    ComplexVector B(space_op.GetCurlMatrix().Height());
    B.UseDevice(true);
    space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
    space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
    B *= -1.0 / (std::complex<double>(0.0, 1.0) * omega);

    auto close = [](std::complex<double> x, std::complex<double> ref, double tol)
    { return std::abs(x - ref) <= tol * std::max(std::abs(ref), 1.0e-300); };
    for (const int excitation_idx : {1, 2})
    {
      CAPTURE(excitation_idx);
      std::map<int, std::complex<double>> V, S, V_gamma, S_gamma;
      env.PortValues(excitation_idx, omega, E, V, S);
      env.PortValues(excitation_idx, omega, e_gamma, V_gamma, S_gamma);
      REQUIRE(V.size() == 1);
      REQUIRE(V_gamma.size() == 1);
      REQUIRE(std::abs(V.at(2)) > 0.0);
      CHECK(close(V_gamma.at(2), V.at(2), 1.0e-12));
      REQUIRE(S.size() == (wave ? 1 : 0));
      REQUIRE(S_gamma.size() == S.size());
      if (wave)
      {
        REQUIRE(std::abs(S.at(3)) > 0.0);
        CHECK(close(S_gamma.at(3), S.at(3), 1.0e-12));
      }

      post_op.SetLumpedPortVoltages(V);
      post_op.SetWavePortOverlaps(S);
      post_op.MeasureAndPrintAll(excitation_idx, 0, E, B, omega);
      const Measurement full = post_op.GetMeasurement();
      post_op.SetLumpedPortVoltages(V_gamma);
      post_op.SetWavePortOverlaps(S_gamma);
      post_op.MeasureAndPrintReduced(excitation_idx, 1, &E, omega, y, prom_op);
      REQUIRE(post_op.HasReducedPostprocessing());
      post_op.MeasureAndPrintReduced(excitation_idx, 2, nullptr, omega, y, prom_op);
      const Measurement reduced = post_op.GetMeasurement();
      REQUIRE(reduced.lumped_port_vi.size() == full.lumped_port_vi.size());
      for (const auto &[idx, vi] : full.lumped_port_vi)
      {
        CAPTURE(idx);
        CHECK(close(reduced.lumped_port_vi.at(idx).V, vi.V, 1.0e-10));
        CHECK(close(reduced.lumped_port_vi.at(idx).S, vi.S, 1.0e-10));
      }
      REQUIRE(reduced.wave_port_vi.size() == full.wave_port_vi.size());
      for (const auto &[idx, vi] : full.wave_port_vi)
      {
        CAPTURE(idx);
        CHECK(close(reduced.wave_port_vi.at(idx).S, vi.S, 1.0e-10));
      }
      REQUIRE(reduced.domain_E_field_energy_i.size() ==
              full.domain_E_field_energy_i.size());
      for (std::size_t i = 0; i < full.domain_E_field_energy_i.size(); i++)
      {
        CHECK(close(*reduced.domain_E_field_energy_i[i].energy,
                    *full.domain_E_field_energy_i[i].energy, 1.0e-10));
      }
    }
  }
  // The full online postprocessing prints the norm of the expanded solution at each
  // frequency, which the reduced one does not compute.
  const auto out = CaptureStdout([&]() { RunDrivenBox(temp_dir, "Online", edit); });
  CHECK(out.find("Sol. ||E|| = ") == std::string::npos);
  CHECK(out.find("Reduced postprocessing disagrees") == std::string::npos);
  if (Mpi::Root(Mpi::World()))
  {
    CHECK(out.find("Beginning fast frequency sweep online phase") != std::string::npos);
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring online with a process holding no region element",
                 "[substructure][Parallel]")
{
  // An online sweep partitioned with every region element on the first process, so that
  // the others hold environment elements only, whose DOFs the region eliminates, gives the
  // S-parameters of the offline sweep.
  const int nranks = Mpi::Size(Mpi::World());
  if (nranks < 2)
  {
    return;
  }
  RunDrivenBox(temp_dir, "Offline");
  const auto part = temp_dir / "box_region_on_root.part";
  if (Mpi::Root(Mpi::World()))
  {
    const mfem::Mesh serial = MakeInteriorBoundaryBox();
    std::ofstream f(part);
    f << "number_of_elements " << serial.GetNE() << "\nnumber_of_processors " << nranks
      << "\n";
    for (int e = 0; e < serial.GetNE(); e++)
    {
      f << ((serial.GetAttribute(e) == 1) ? 0 : 1 + e % (nranks - 1)) << "\n";
    }
  }
  Mpi::Barrier(Mpi::World());
  RunDrivenBox(temp_dir, "Online",
               [&](json &config) { config["Model"]["Partitioning"] = part.string(); });
  const auto offline = ReadCsv(temp_dir / "postpro_Offline" / "port-S.csv");
  const auto online = ReadCsv(temp_dir / "postpro_Online" / "port-S.csv");
  REQUIRE(online.first == offline.first);
  REQUIRE(online.second.size() == offline.second.size());
  for (std::size_t r = 0; r < online.second.size(); r++)
  {
    for (std::size_t c = 0; c < online.second[r].size(); c++)
    {
      CAPTURE(online.first[c]);
      CHECK(std::abs(online.second[r][c] - offline.second[r][c]) <=
            1.0e-8 * (1.0 + std::abs(offline.second[r][c])));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring needs the sources of a model away from Γ",
                 "[substructure][Serial][Parallel]")
{
  // A surface current on Γ, the interior boundary between the domains of the box, loads Γ,
  // which a model cannot hold: an offline sweep saving a model and an online one from it
  // reject it.
  auto current = [](json &config)
  {
    config["Boundaries"]["SurfaceCurrent"] = {
        {{"Index", 2}, {"Attributes", {4}}, {"Direction", "+Y"}}};
  };
  using Catch::Matchers::ContainsSubstring;
  CHECK_THROWS_WITH(RunDrivenBox(temp_dir, "Offline", current),
                    ContainsSubstring("has sources on the substructuring interface"));
  RunDrivenBox(temp_dir, "Offline");
  CHECK_THROWS_WITH(RunDrivenBox(temp_dir, "Online", current),
                    ContainsSubstring("has sources on the substructuring interface"));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring tells the ports of the environment by their attributes",
    "[substructure][Serial][Parallel]")
{
  // The ports 2 and 3 of the environment with their indices swapped online, which leaves
  // the configuration of the environment and the excitations as they are: the saved port
  // functionals are of the other port.
  auto ports = [](bool swapped)
  {
    return [swapped](json &config)
    {
      config["Boundaries"]["PEC"]["Attributes"] = {2, 3, 4, 6, 7, 8};
      config["Boundaries"]["LumpedPort"] = {{{"Index", 1},
                                             {"R", 50.0},
                                             {"Attributes", {9}},
                                             {"Direction", "+X"},
                                             {"Excitation", 1}},
                                            {{"Index", swapped ? 3 : 2},
                                             {"R", 50.0},
                                             {"Attributes", {10}},
                                             {"Direction", "+X"}},
                                            {{"Index", swapped ? 2 : 3},
                                             {"R", 75.0},
                                             {"Attributes", {5}},
                                             {"Direction", "+Y"}}};
    };
  };
  RunDrivenSheetSquare(temp_dir, "Offline", ports(false));
  CHECK_NOTHROW(RunDrivenSheetSquare(temp_dir, "Online", ports(false)));
  CHECK_THROWS_WITH(RunDrivenSheetSquare(temp_dir, "Online", ports(true)),
                    Catch::Matchers::ContainsSubstring("in their indices or attributes"));
}

#if defined(MFEM_USE_GSLIB)
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring writes the probes it interpolates in the region",
                 "[substructure][Serial][Parallel][Long]")
{
  // Two probes on Γ and one in the region. Online on the full mesh, all three are written,
  // those on Γ with the field of the region, on its side of Γ, where the normal E jumps
  // (the permittivities of the sides differ): the values of a regular sweep with the probes
  // on Γ moved into the region by 1e-9.
  auto probes = [](double shift)
  {
    return [shift](json &config)
    {
      config["Domains"]["Postprocessing"]["Probe"] = {
          {{"Index", 1}, {"Center", {0.5 - shift, -3.1, 4.7}}},
          {{"Index", 2}, {"Center", {0.5 - shift, 0.3, 0.2}}},
          {{"Index", 3}, {"Center", {-6.0, 0.0, -3.0}}}};
    };
  };
  RunDrivenSheetSquare(temp_dir, "Offline", probes(0.0));
  RunDrivenSheetSquare(temp_dir, "Online", probes(0.0));
  RunDrivenSheetSquare(temp_dir, "Regular",
                       [&](json &config)
                       {
                         probes(1.0e-9)(config);
                         config["Solver"].erase("Substructuring");
                         config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 1000}};
                       });
  for (const std::string csv : {"probe-E.csv", "probe-B.csv"})
  {
    CAPTURE(csv);
    const auto online = ReadCsv(temp_dir / "postpro_Online" / csv);
    const auto regular = ReadCsv(temp_dir / "postpro_Regular" / csv);
    REQUIRE(online.second.size() == regular.second.size());
    double scale = 0.0;
    for (const auto &row : regular.second)
    {
      for (std::size_t c = 1; c < row.size(); c++)
      {
        scale = std::max(scale, std::abs(row[c]));
      }
    }
    std::set<std::string> written;
    for (std::size_t c = 1; c < online.first.size(); c++)
    {
      const auto &name = online.first[c];
      CAPTURE(name);
      written.insert(name.substr(name.find('[') + 1, 1));
      const auto it = std::ranges::find(regular.first, name);
      REQUIRE(it != regular.first.end());
      const auto k = static_cast<std::size_t>(it - regular.first.begin());
      for (std::size_t r = 0; r < online.second.size(); r++)
      {
        CHECK(std::abs(online.second[r][c] - regular.second[r][k]) <= 1.0e-8 * scale);
      }
    }
    CHECK(written == std::set<std::string>{"1", "2", "3"});
  }
}
#endif

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring writes each excitation's fields to its own collection",
    "[substructure][Serial][Parallel]")
{
  // Two excitations at two frequencies: as in a regular sweep, the ParaView collections of
  // each excitation hold its fields at both frequencies (the last one also the final mesh
  // partitioning, at a later time), and not those of an earlier sweep in the same folder.
  auto sweep = [&](double min_freq, double max_freq)
  {
    RunDrivenBox(temp_dir, "Offline",
                 [=](json &config)
                 {
                   config["Boundaries"]["LumpedPort"].push_back({{"Index", 2},
                                                                 {"R", 50.0},
                                                                 {"Attributes", {5}},
                                                                 {"Direction", "+Y"},
                                                                 {"Excitation", 2}});
                   config["Solver"]["Driven"] = {{"MinFreq", min_freq},
                                                 {"MaxFreq", max_freq},
                                                 {"FreqStep", 10.0},
                                                 {"SaveStep", 1}};
                 });
  };
  sweep(10.0, 10.0);
  sweep(20.0, 30.0);
  const auto root = temp_dir / "postpro_Offline" / "paraview";
  for (const std::string collection : {"driven", "driven_boundary"})
  {
    for (int excitation : {1, 2})
    {
      const std::string folder = fmt::format("excitation_{}", excitation);
      CAPTURE(collection, excitation);
      std::ifstream pvd(root / collection / folder / (folder + ".pvd"));
      REQUIRE(pvd.good());
      const std::string contents(std::istreambuf_iterator<char>(pvd), {});
      const std::regex timestep("timestep=\"([^\"]+)\"");
      std::vector<double> freq;
      for (auto it = std::sregex_iterator(contents.begin(), contents.end(), timestep);
           it != std::sregex_iterator(); ++it)
      {
        const double t = std::stod((*it)[1].str());
        if (t < 90.0)
        {
          freq.push_back(t);
        }
      }
      REQUIRE(freq.size() == 2);
      CHECK_THAT(freq[0], Catch::Matchers::WithinRel(20.0, 1.0e-6));
      CHECK_THAT(freq[1], Catch::Matchers::WithinRel(30.0, 1.0e-6));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring tests its rational model against its samples",
                 "[substructure][Serial][Parallel][Long]")
{
  // A second-order absorbing boundary makes Im S_E indefinite at the samples themselves:
  // the rational model is passive between them if no less so than at them, so that the
  // sampling stops at its tolerance, without passivity samples or a warning.
  const auto dir = std::filesystem::path(PALACE_TEST_DATA_DIR) / "regression" / "input" /
                   "substructuring_driven";
  std::ifstream fi(dir / "driven_offline_absorbing.json");
  json config = json::parse(fi);
  config["Model"]["Mesh"] = (dir / config["Model"]["Mesh"].get<std::string>()).string();
  config["Problem"]["Output"] = (temp_dir / "absorbing").string();
  config["Boundaries"]["Absorbing"]["Order"] = 2;
  config["Solver"]["Driven"]["FreqStep"] = 250.0;
  config["Solver"]["Driven"]["AdaptiveTol"] = 1.0e-8;
  config["Solver"]["Substructuring"]["SaveModel"] = false;
  IoData iodata(config, false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  const auto out = CaptureStdout([&]() { palace::Run(iodata, comm, 1); });
  if (Mpi::Root(comm))
  {
    CHECK(out.find("Passive between the samples") != std::string::npos);
    CHECK(out.find("Passivity sample") == std::string::npos);
    CHECK(out.find("not passive") == std::string::npos);
  }
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring gives no voltage of the environment's wave port online",
    "[substructure][Serial][Parallel][Long]")
{
  // The voltage and the impedance of a wave port need its field: offline they are those of
  // the environment's port, online, where the port's field is not known, they are NaN, not
  // 0 (its S-parameter comes from the model).
  const auto dir = std::filesystem::path(PALACE_TEST_DATA_DIR) / "regression" / "input" /
                   "substructuring_driven_wave";
  auto run = [&](const std::string &file)
  {
    std::ifstream fi(dir / file);
    json config = json::parse(fi);
    config["Model"]["Mesh"] = (dir / config["Model"]["Mesh"].get<std::string>()).string();
    config["Problem"]["Output"] = (temp_dir / file).string();
    config["Solver"]["Linear"]["Type"] = "Default";
    if (config["Solver"]["Substructuring"].contains("ModelFile"))
    {
      config["Solver"]["Substructuring"]["ModelFile"] =
          (temp_dir / "driven_wave_inactive_offline.json" / "environment.model").string();
    }
    config["Boundaries"]["WavePort"][1]["VoltagePath"] = {{0.0, 0.0, 0.0},
                                                          {2.74, 0.0, 0.0}};
    IoData iodata(config, false);
    MPI_Comm comm = Mpi::World();
    MakeOutputFolder(iodata, comm);
    BlockTimer::Reset();
    palace::Run(iodata, comm, 1);
  };
  // The values of the columns of a CSV file whose header holds name.
  auto values =
      [&](const std::string &file, const std::string &csv, const std::string &name)
  {
    const auto csv_data = ReadCsv(temp_dir / file / csv);
    const auto &header = csv_data.first;
    std::vector<double> v;
    for (const auto &row : csv_data.second)
    {
      for (std::size_t c = 0; c < header.size() && c < row.size(); c++)
      {
        if (header[c].find(name) != std::string::npos)
        {
          v.push_back(row[c]);
        }
      }
    }
    return v;
  };
  run("driven_wave_inactive_offline.json");
  run("driven_wave_inactive_online.json");
  for (const std::string name : {"V_wp[2]", "Z[2]"})
  {
    CAPTURE(name);
    const auto off = values("driven_wave_inactive_offline.json",
                            (name == "Z[2]") ? "port-Z.csv" : "port-V.csv", name);
    const auto on = values("driven_wave_inactive_online.json",
                           (name == "Z[2]") ? "port-Z.csv" : "port-V.csv", name);
    REQUIRE(!off.empty());
    REQUIRE(on.size() == off.size());
    CHECK(std::ranges::all_of(off, [](double x) { return std::isfinite(x); }));
    CHECK(std::ranges::any_of(off, [](double x) { return x != 0.0; }));
    CHECK(std::ranges::all_of(on, [](double x) { return std::isnan(x); }));
  }
}

TEST_CASE("DrivenSubstructure finds the boundary attributes on Γ",
          "[substructure][Serial][Parallel]")
{
  // In the box, the interior boundary between the domains is Γ, with its DOFs on the PEC
  // faces around it; the PEC faces themselves, all Dirichlet, and the interior boundary
  // inside the environment are not.
  json config = {{"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1, 2}}}}}}},
                 {"Boundaries",
                  {{"PEC", {{"Attributes", {3}}}},
                   {"LumpedPort",
                    {{{"Index", 1},
                      {"R", 50.0},
                      {"Attributes", {1}},
                      {"Direction", "+Y"},
                      {"Excitation", true}}}}}},
                 {"Solver",
                  {{"Device", "CPU"},
                   {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 1.0}}}}}};
  IoData iodata(config, false);
  mfem::Mesh serial = MakeInteriorBoundaryBox();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  SpaceOperator space_op(iodata, mesh);
  DrivenSubstructure ds(space_op, {1}, {2});
  CHECK(ds.InterfaceBoundaryAttributes({3, 4, 5}) == std::vector<int>{4});
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring on the mesh of the region rejects boundary terms on Γ",
    "[substructure][Serial][Parallel]")
{
  // An impedance (also a rational one) or a conducting sheet on Γ, the interior boundary
  // between the domains of the box, would be the region's alone on the mesh of the region,
  // while the full mesh runs from the model. An impedance on the interior boundary inside
  // the environment is the model's.
  const std::string region_mesh = RegionMeshPath(DrivenModel(temp_dir).string());
  using Edit = std::function<void(json &)>;
  const std::vector<std::tuple<std::string, Edit, bool>> terms = {
      {"impedance", [](json &config)
       { config["Boundaries"]["Impedance"] = {{{"Attributes", {4}}, {"Rs", 377.0}}}; },
       true},
      {"rational impedance",
       [](json &config)
       {
         config["Boundaries"]["RationalImpedance"] = {
             {{"Attributes", {4}}, {"Numerator", {377.0}}, {"Denominator", {1.0}}}};
       },
       true},
      {"conductivity",
       [](json &config)
       {
         config["Boundaries"]["Conductivity"] = {
             {{"Attributes", {4}}, {"Conductivity", 1.0e6}}};
       },
       true},
      {"impedance in the environment", [](json &config)
       { config["Boundaries"]["Impedance"] = {{{"Attributes", {5}}, {"Rs", 377.0}}}; },
       false}};
  using Catch::Matchers::ContainsSubstring;
  for (const auto &t : terms)
  {
    const std::string &name = std::get<0>(t);
    const Edit &term = std::get<1>(t);
    const bool on_gamma = std::get<2>(t);
    CAPTURE(name);
    RunDrivenBox(temp_dir, "Offline", term);
    RunDrivenBox(temp_dir, "Online", term);
    auto on_region = [&](json &config)
    {
      term(config);
      config["Model"]["Mesh"] = region_mesh;
    };
    if (on_gamma)
    {
      CHECK_THROWS_WITH(RunDrivenBox(temp_dir, "Online", on_region),
                        ContainsSubstring("boundary conditions on the substructuring "
                                          "interface (attributes 4)"));
    }
    else
    {
      RunDrivenBox(temp_dir, "Online", on_region);
    }
  }
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring needs the full mesh for a source of the environment",
    "[substructure][Serial][Parallel]")
{
  // A current dipole in the environment, away from Γ, lies outside the mesh of the region,
  // which holds no source of it: a sweep on it is rejected, and one on the full mesh runs.
  auto dipole = [](json &config)
  {
    config["Domains"]["CurrentDipole"] = {{{"Index", 1},
                                           {"Moment", 1.0e-3},
                                           {"Center", {1.75, 0.25, 0.75}},
                                           {"Direction", {0.0, 0.0, 1.0}}}};
    config["Boundaries"]["LumpedPort"][0].erase("Excitation");
  };
  RunDrivenBox(temp_dir, "Offline", dipole);
  CHECK_NOTHROW(RunDrivenBox(temp_dir, "Online", dipole));
  const std::string region_mesh = RegionMeshPath(DrivenModel(temp_dir).string());
  CHECK_THROWS_WITH(
      RunDrivenBox(temp_dir, "Online",
                   [&](json &config)
                   {
                     dipole(config);
                     config["Model"]["Mesh"] = region_mesh;
                   }),
      Catch::Matchers::ContainsSubstring("the mesh of the region does not hold"));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring checks the environment against the model",
                 "[substructure][Serial][Parallel][Long]")
{
  // The record of the environment tells the mesh of the region from the full mesh, and
  // rejects an environment refined, with another material, or with another resistance of
  // its port on the mesh of the region (which holds no environment to fingerprint). The
  // excitation of the environment's port is not its configuration: a sweep without it runs,
  // as does one with it on the mesh of the region, which fingerprints its sources there.
  RunDrivenSheetSquare(temp_dir, "Offline");
  const std::string region_mesh = RegionMeshPath(DrivenModel(temp_dir).string());
  auto on_region = [&](json &config) { config["Model"]["Mesh"] = region_mesh; };
  auto online = [&](const std::function<void(json &)> &edit)
  { RunDrivenSheetSquare(temp_dir, "Online", edit); };
  using Catch::Matchers::ContainsSubstring;
  CHECK_THROWS_WITH(
      online([](json &config) { config["Model"]["Refinement"]["UniformLevels"] = 1; }),
      ContainsSubstring("its mesh changed"));
  CHECK_THROWS_WITH(
      online([](json &config) { config["Domains"]["Materials"][1]["Permittivity"] = 3.0; }),
      ContainsSubstring("discretization changed"));
  CHECK_THROWS_WITH(online(
                        [&](json &config)
                        {
                          on_region(config);
                          config["Boundaries"]["LumpedPort"][1]["R"] = 60.0;
                        }),
                    ContainsSubstring("discretization changed"));
  online(
      [&](json &config)
      {
        on_region(config);
        config["Boundaries"]["LumpedPort"][1].erase("Excitation");
      });
  online(on_region);
}

TEST_CASE("The record of a driven environment covers its ports", "[substructure][Serial]")
{
  // The configuration hash of the environment covers every key of its ports but their
  // excitations, so also the direction of a lumped port and the mode target of a wave port,
  // whose functionals give the voltages of an online sweep.
  const json base = json::parse(R"({
    "Problem": {"Type": "Driven", "Output": "postpro"},
    "Model": {"Mesh": "mesh.msh"},
    "Domains": {"Materials": [{"Attributes": [1, 2], "Permittivity": 2.0}]},
    "Boundaries": {
      "LumpedPort": [{"Index": 1, "R": 50.0, "Attributes": [3], "Direction": "+X",
                      "Excitation": true}],
      "WavePort": [{"Index": 2, "Attributes": [4], "Mode": 1, "Target": 0.65}]},
    "Solver": {"Order": 2}})");
  const std::vector<int> attrs = {2}, bdr = {3, 4};
  const auto hash = EnvironmentConfigHash(base, attrs, bdr, false);
  auto differs = [&](const std::function<void(json &)> &change)
  {
    json config = base;
    change(config);
    return EnvironmentConfigHash(config, attrs, bdr, false) != hash;
  };
  CHECK(differs([](json &c) { c["Boundaries"]["LumpedPort"][0]["Direction"] = "-X"; }));
  CHECK(differs([](json &c) { c["Boundaries"]["WavePort"][0]["Target"] = 1.1; }));
  CHECK_FALSE(
      differs([](json &c) { c["Boundaries"]["LumpedPort"][0].erase("Excitation"); }));
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring warns about a rational model which did not converge",
    "[substructure][Serial][Parallel][Long]")
{
  // An online sweep warns about a rational model whose sampling stopped before it met
  // "AdaptiveTol": at "AdaptiveMaxSamples", or in an offline run cut short, which leaves
  // the header of its last sample.
  auto offline = [&](int max_samples)
  {
    RunDrivenSheetSquare(temp_dir, "Offline",
                         [&](json &config)
                         {
                           config["Solver"]["Driven"] = {
                               {"MinFreq", 1000.0},
                               {"MaxFreq", 3000.0},
                               {"FreqStep", 250.0},
                               {"AdaptiveTol", 1.0e-8},
                               {"AdaptiveMaxSamples", max_samples}};
                         });
  };
  auto online = [&]()
  {
    return CaptureStdout(
        [&]()
        {
          RunDrivenSheetSquare(temp_dir, "Online",
                               [](json &config)
                               {
                                 config["Solver"]["Driven"] = {{"MinFreq", 1500.0},
                                                               {"MaxFreq", 2500.0},
                                                               {"FreqStep", 1000.0}};
                               });
        });
  };
  offline(3);
  const auto capped = online();
  offline(20);
  const auto converged = online();
  const bool root = Mpi::Root(Mpi::World());
  if (root)
  {
    // The model without its final header.
    using substructuring::ModelFile;
    const auto path = DrivenModel(temp_dir).string();
    long long pos = 0;
    for (const auto &r : ModelFile::Read(path, MPI_COMM_SELF).records)
    {
      pos = (r.tag == ModelFile::Tag("DHED")) ? r.pos : pos;
    }
    CHECK(pos > 0);
    if (pos > 0)
    {
      std::filesystem::resize_file(
          path, pos - static_cast<long long>(sizeof(std::uint32_t) + sizeof(long long)));
    }
  }
  Mpi::Barrier(Mpi::World());
  const auto cut = online();
  if (root)
  {
    using Catch::Matchers::ContainsSubstring;
    CHECK_THAT(capped, ContainsSubstring("did not converge"));
    CHECK_THAT(converged, !ContainsSubstring("did not converge"));
    CHECK_THAT(cut, ContainsSubstring("did not converge"));
  }
}

TEST_CASE_METHOD(
    palace::test::SharedTempDir,
    "Driven substructuring rejects an environment reflected through the origin",
    "[substructure][Serial][Parallel]")
{
  // A box [-2, 2] x [-1, 1]^2 of 8 x 4 x 4 hexahedra about the origin: the region x > 1
  // (attribute 1) and the environment, which holds two elements reflected into each other
  // through the origin (3 and 4), an object of four elements centered at the origin, away
  // from the boundary (5), and the rest (2). Swapping the materials of 3 and 4 changes the
  // configuration of the environment. Reflecting the object through the origin moves the
  // boundary elements added on its interfaces, but changes no moment of an attribute and
  // no energy of a field centered at the origin: without the interface elements, only
  // the fingerprint of the operator, with fields centered off the origin, tells it.
  auto write = [&](bool reflect)
  {
    const auto path = temp_dir / (reflect ? "reflected.mesh" : "mirror.mesh");
    if (Mpi::Root(Mpi::World()))
    {
      auto mesh =
          mfem::Mesh::MakeCartesian3D(8, 4, 4, mfem::Element::HEXAHEDRON, 4.0, 2.0, 2.0);
      for (int v = 0; v < mesh.GetNV(); v++)
      {
        double *x = mesh.GetVertex(v);
        x[0] -= 2.0;
        x[1] -= 1.0;
        x[2] -= 1.0;
      }
      const double s = reflect ? -1.0 : 1.0;
      const std::vector<std::array<double, 4>> parts = {
          // element center, attribute
          {-0.75, 0.25, -0.25, 3},
          {0.75, -0.25, 0.25, 4},
          {0.75 * s, 0.25 * s, 0.25 * s, 5},
          {-0.25 * s, -0.25 * s, -0.25 * s, 5},
          {-0.25 * s, 0.25 * s, -0.25 * s, 5},
          {-0.25 * s, -0.25 * s, 0.25 * s, 5}};
      mfem::Vector c;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
        mesh.GetElementCenter(e, c);
        int attr = (c(0) > 1.0) ? 1 : 2;
        for (const auto &p : parts)
        {
          if (std::abs(c(0) - p[0]) + std::abs(c(1) - p[1]) + std::abs(c(2) - p[2]) <
              1.0e-9)
          {
            attr = static_cast<int>(p[3]);
          }
        }
        mesh.SetAttribute(e, attr);
      }
      // The faces x = 2 (lumped port 1), x = -2 (absorbing) and the others (PEC).
      mfem::Array<int> v;
      for (int b = 0; b < mesh.GetNBE(); b++)
      {
        mesh.GetBdrElementVertices(b, v);
        double x = 0.0;
        for (int i : v)
        {
          x += mesh.GetVertex(i)[0] / v.Size();
        }
        mesh.SetBdrAttribute(b, (x > 2.0 - 1.0e-9) ? 1 : (x < -2.0 + 1.0e-9) ? 2 : 3);
      }
      mesh.SetAttributes();
      mesh.Save(path.string());
    }
    Mpi::Barrier(Mpi::World());
    return path.string();
  };
  const auto mesh = write(false), reflected = write(true);
  auto run = [&](const std::string &mode, const std::string &m, double eps3, double eps4,
                 bool interfaces = true)
  {
    RunDrivenBox(
        temp_dir, mode,
        [&](json &config)
        {
          config["Model"]["Mesh"] = m;
          config["Model"]["AddInterfaceBoundaryElements"] = interfaces;
          config["Domains"]["Materials"] = {{{"Attributes", {1}}, {"Permittivity", 4.0}},
                                            {{"Attributes", {2}}, {"Permittivity", 2.0}},
                                            {{"Attributes", {3}}, {"Permittivity", eps3}},
                                            {{"Attributes", {4}}, {"Permittivity", eps4}},
                                            {{"Attributes", {5}}, {"Permittivity", 8.0}}};
          config["Solver"]["Substructuring"]["Environment"]["Attributes"] = {2, 3, 4, 5};
        });
  };
  run("Offline", mesh, 3.0, 6.0);
  using Catch::Matchers::ContainsSubstring;
  CHECK_THROWS_WITH(run("Online", mesh, 6.0, 3.0),
                    ContainsSubstring("materials, boundary conditions or discretization"));
  CHECK_THROWS_WITH(run("Online", reflected, 3.0, 6.0),
                    ContainsSubstring("its mesh changed"));
  run("Online", mesh, 3.0, 6.0);
  run("Offline", mesh, 3.0, 6.0, false);
  CHECK_THROWS_WITH(run("Online", reflected, 3.0, 6.0, false),
                    ContainsSubstring("(environment mesh, materials"));
  run("Online", mesh, 3.0, 6.0, false);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Driven substructuring mesh of the region takes no a priori refinement",
                 "[substructure][Serial][Parallel][Long]")
{
  // As in the static substructuring, the mesh of the region is written refined: an online
  // run on it takes no a priori refinement, also of the region only (which leaves the
  // written environment elements as they were).
  RunDrivenSheetSquare(temp_dir, "Offline");
  const std::string region_mesh = RegionMeshPath(DrivenModel(temp_dir).string());
  const json uniform = {{"UniformLevels", 1}},
             boxes = {{"Nonconformal", true},
                      {"MaxNCLevels", 0},
                      {"Boxes",
                       {{{"Levels", 1},
                         {"BoundingBoxMin", {-8.4, -7.8, -7.1}},
                         {"BoundingBoxMax", {-7.5, -6.9, -6.2}}}}}};
  for (const json &refinement : {uniform, boxes})
  {
    CHECK_THROWS_WITH(RunDrivenSheetSquare(temp_dir, "Online",
                                           [&](json &config)
                                           {
                                             config["Model"]["Mesh"] = region_mesh;
                                             config["Model"]["Refinement"] = refinement;
                                           }),
                      Catch::Matchers::ContainsSubstring("is written refined"));
  }
}

#endif

#if defined(MFEM_USE_MUMPS) || defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK)

#if defined(MFEM_USE_MUMPS)
TEST_CASE("MumpsSchurSolver refactors with new values", "[substructure][Serial][Parallel]")
{
  // A complex symmetric shifted Laplacian with Schur variables on a face of the cube,
  // refactored with its values doubled: the Schur complement doubles, and the internal
  // solves halve.
  auto pmesh = MakeSplitCube(4);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fes(pmesh.get(), &fec);
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.AddDomainIntegrator(new mfem::MassIntegrator);
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize(), nranks = Mpi::Size(comm);
  mfem::Array<int> bdr(pmesh->bdr_attributes.Max()), face;
  bdr = 0;
  bdr[0] = 1;
  fes.GetEssentialTrueDofs(bdr, face);
  std::vector<HYPRE_BigInt> mine(face.Size()), schur_vars;
  for (int i = 0; i < face.Size(); i++)
  {
    mine[i] = fes.GetMyTDofOffset() + face[i];
  }
  {
    int n = face.Size();
    std::vector<int> cnt(nranks), disp(nranks, 0);
    MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    schur_vars.resize(disp.back() + cnt.back());
    MPI_Allgatherv(mine.data(), n, HYPRE_MPI_BIG_INT, schur_vars.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
  }
  REQUIRE(!schur_vars.empty());
  ComplexVector x(nt);
  x.Real().Randomize(1 + Mpi::Rank(comm));
  x.Imag().Randomize(2 + Mpi::Rank(comm));

  ComplexMumpsSchurSolver::Coo coo;
  std::vector<std::vector<double>> vals;
  std::vector<int> row_ptr;
  LowerTrianglePattern({A.get()}, coo.irn, coo.jcn, vals, row_ptr);
  for (double v : vals[0])
  {
    coo.val.push_back(std::complex<double>(1.0, 0.1) * v);
  }
  ComplexMumpsSchurSolver solver(comm, A->GetGlobalNumRows(), nt, std::move(coo),
                                 schur_vars, 0.0, true, false);
  std::vector<std::complex<double>> S[2];
  ComplexVector y[2];
  for (int k = 0; k < 2; k++)
  {
    if (k == 1)
    {
      for (auto &v : solver.Values())
      {
        v *= 2.0;
      }
      solver.Refactor();
    }
    S[k] = solver.ReleaseSchur();
    solver.SolveInternal({&x}, {&y[k]});
  }
  if (Mpi::Root(comm))
  {
    double d = 0.0, m = 0.0;
    for (std::size_t q = 0; q < S[0].size(); q++)
    {
      d = std::max(d, std::abs(S[1][q] - 2.0 * S[0][q]));
      m = std::max(m, std::abs(2.0 * S[0][q]));
    }
    CHECK(d <= 1.0e-12 * m);
  }
  ComplexVector d(y[1]);
  d *= 2.0;
  d -= y[0];
  CHECK(linalg::Norml2(comm, d) <= 1.0e-12 * linalg::Norml2(comm, y[0]));
}

TEST_CASE("MumpsSchurSolver keeps its tolerance relative when refactored",
          "[substructure][Serial][Parallel]")
{
  // The block low-rank (BLR) tolerance is relative to the largest diagonal entry of the
  // values factored: refactored with values 1e6 times smaller, the Schur complement keeps
  // the accuracy of the first factorization (with the scale of the first values, it loses
  // three digits here).
  auto pmesh = MakeSplitCube(12);
  mfem::H1_FECollection fec(2, 3);
  mfem::ParFiniteElementSpace fes(pmesh.get(), &fec);
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
  a.AddDomainIntegrator(new mfem::MassIntegrator);
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  MPI_Comm comm = fes.GetComm();
  const int nt = fes.GetTrueVSize(), nranks = Mpi::Size(comm);
  mfem::Array<int> bdr(pmesh->bdr_attributes.Max()), face;
  bdr = 0;
  bdr[0] = 1;
  fes.GetEssentialTrueDofs(bdr, face);
  std::vector<HYPRE_BigInt> mine(face.Size()), schur_vars;
  for (int i = 0; i < face.Size(); i++)
  {
    mine[i] = fes.GetMyTDofOffset() + face[i];
  }
  {
    int n = face.Size();
    std::vector<int> cnt(nranks), disp(nranks, 0);
    MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    schur_vars.resize(disp.back() + cnt.back());
    MPI_Allgatherv(mine.data(), n, HYPRE_MPI_BIG_INT, schur_vars.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
  }
  auto schur = [&](double blr_tol)
  {
    ComplexMumpsSchurSolver::Coo coo;
    std::vector<std::vector<double>> vals;
    std::vector<int> row_ptr;
    LowerTrianglePattern({A.get()}, coo.irn, coo.jcn, vals, row_ptr);
    for (double v : vals[0])
    {
      coo.val.push_back(std::complex<double>(1.0, 0.1) * v);
    }
    const auto values = coo.val;
    ComplexMumpsSchurSolver solver(comm, A->GetGlobalNumRows(), nt, std::move(coo),
                                   schur_vars, blr_tol, true, false);
    std::array<std::vector<std::complex<double>>, 2> S;
    S[0] = solver.ReleaseSchur();
    for (std::size_t q = 0; q < values.size(); q++)
    {
      solver.Values()[q] = 1.0e-6 * values[q];
    }
    solver.Refactor();
    S[1] = solver.ReleaseSchur();
    return S;
  };
  const auto exact = schur(0.0), blr = schur(1.0e-8);
  if (Mpi::Root(comm))
  {
    for (int k = 0; k < 2; k++)
    {
      double d = 0.0, m = 0.0;
      for (std::size_t q = 0; q < exact[k].size(); q++)
      {
        d = std::max(d, std::abs(blr[k][q] - exact[k][q]));
        m = std::max(m, std::abs(exact[k][q]));
      }
      CAPTURE(k, d, m);
      CHECK(d <= 1.0e-8 * m);
    }
  }
}
#endif

TEST_CASE("DrivenSubstructure rejects attribute sets that do not split the mesh",
          "[substructure][Serial][Parallel]")
{
  json config = {{"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
                 {"Model", {{"Mesh", "test.msh"}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1, 2}}}}}}},
                 {"Boundaries",
                  {{"PEC", {{"Attributes", {1}}}},
                   {"LumpedPort",
                    {{{"Index", 1},
                      {"R", 1.0},
                      {"Attributes", {2}},
                      {"Direction", "+Y"},
                      {"Excitation", true}}}}}},
                 {"Solver",
                  {{"Device", "CPU"},
                   {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 1.0}}}}}};
  IoData iodata(config, false);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeGradedSplit(2, 3, 2)));
  SpaceOperator space_op(iodata, mesh);
  CHECK_THROWS_WITH(DrivenSubstructure(space_op, {1}, {3}),
                    Catch::Matchers::ContainsSubstring("Domain attribute 2 of the mesh is "
                                                       "in neither"));
  CHECK_THROWS_WITH(
      DrivenSubstructure(space_op, {1, 4}, {2}),
      Catch::Matchers::ContainsSubstring("\"Region\" attribute 4 is not a domain"));
}

TEST_CASE("DrivenSubstructure condenses the environment exactly",
          "[substructure][Serial][Parallel][Long]")
{
  // S_E(ω) from the partial factorization of the environment operator (alone, and with the
  // region) against a dense condensation of the same operator, and the substructured solve
  // against the full system: PEC in the region, a lumped port in the environment, a
  // second-order absorbing boundary on both sides, an impedance sheet crossing Γ, a lossy
  // dielectric region and a conducting environment. A second frequency reuses the analyses.
  // Hex meshes with an impedance sheet, conforming and nonconforming (hanging nodes on Γ
  // and on the sheet from either side), and tet meshes with a non-planar interface (where
  // second-order Nédélec face DOFs shared between ranks combine with signs).
  enum class MeshKind
  {
    HEX,
    NC_HEX,
    TET
  };
  const auto kind = GENERATE(MeshKind::HEX, MeshKind::NC_HEX, MeshKind::TET);
  const int order = GENERATE(1, 2);
  // The default sparse direct solvers, and each one Palace is built with.
  std::vector<LinearSolver> types = {LinearSolver::DEFAULT};
#if defined(MFEM_USE_STRUMPACK)
  types.push_back(LinearSolver::STRUMPACK);
#endif
#if defined(MFEM_USE_MUMPS)
  types.push_back(LinearSolver::MUMPS);
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back(LinearSolver::SUPERLU);
#endif
  const LinearSolver type = GENERATE_COPY(from_range(types));
  const bool tet = (kind == MeshKind::TET);
  CAPTURE(order, static_cast<int>(kind), static_cast<int>(type));
  // SuperLU_DIST's static pivoting, in the real equivalent form, is accurate to a few
  // 1e-10 here: by default, it factors the sides without MUMPS and STRUMPACK, and the
  // online region without STRUMPACK.
#if defined(MFEM_USE_SUPERLU) && !defined(MFEM_USE_STRUMPACK)
  const bool superlu = (type == LinearSolver::SUPERLU || type == LinearSolver::DEFAULT);
#else
  const bool superlu = (type == LinearSolver::SUPERLU);
#endif
  const double exact_tol = superlu ? 1.0e-9 : 1.0e-10;
  json config = {
      {"Problem", {{"Type", "Driven"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"LossTan", 0.01}},
          {{"Attributes", {2}}, {"Permittivity", 4.0}, {"Conductivity", 0.5}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", {1}}}},
        {"LumpedPort",
         {{{"Index", 1},
           {"R", 1.0},
           {"Attributes", {2}},
           {"Direction", "+Y"},
           {"Excitation", true}}}},
        {"Absorbing", {{"Attributes", {3}}, {"Order", 2}}},
        {"Impedance", {{{"Attributes", {4}}, {"Rs", 1.0}, {"Ls", 2.0}, {"Cs", 0.5}}}}}},
      {"Solver",
       {{"Order", order},
        {"Device", "CPU"},
        {"Driven", {{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 1.0}}}}}};
  if (tet)
  {
    config["Boundaries"].erase("Impedance");
  }
  IoData iodata(config, false);
  RegionDesign design;
  design.sheet = true;
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(tet ? MakeWavyTetSplit(3)
                                        : (kind == MeshKind::NC_HEX)
                                            ? MakeRefinedGradedSplit(2, 3, 2, design)
                                            : MakeGradedSplit(2, 3, 2, design)));
  REQUIRE(mesh.back()->Get().Nonconforming() == (kind == MeshKind::NC_HEX));
  SpaceOperator space_op(iodata, mesh);
  const auto &fes = space_op.GetNDSpace().Get();
  MPI_Comm comm = fes.GetComm();
  const bool root = Mpi::Root(comm);
  DrivenSubstructure ds(space_op, {1}, {2}, false, nullptr, false, type);
  const int nG = ds.InterfaceSize();
  REQUIRE(nG > 0);

  // Global environment-interior and interface true DOFs (replicated), and the environment
  // operator parts, for the dense reference.
  const int nt = fes.GetTrueVSize(), nranks = Mpi::Size(comm);
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  std::vector<HYPRE_BigInt> E_loc, E;
  for (int i = 0; i < nt; i++)
  {
    if (ds.EnvironmentInterior()[i])
    {
      E_loc.push_back(tstart + i);
    }
  }
  {
    int n = static_cast<int>(E_loc.size());
    std::vector<int> cnt(nranks), disp(nranks, 0);
    MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
    E.resize(disp[nranks - 1] + cnt[nranks - 1]);
    MPI_Allgatherv(E_loc.data(), n, HYPRE_MPI_BIG_INT, E.data(), cnt.data(), disp.data(),
                   HYPRE_MPI_BIG_INT, comm);
  }
  const std::vector<HYPRE_BigInt> &G = ds.InterfaceTrueDofs();
  std::vector<HYPRE_BigInt> idx(E);
  idx.insert(idx.end(), G.begin(), G.end());
  const int nE = static_cast<int>(E.size()), n = nE + nG;
  const std::vector<int> environment = {2};
  std::unique_ptr<ComplexOperator> K, C, M;
  {
    AssemblyDomains::Restriction restriction(space_op.GetAssemblyDomains(), environment);
    K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  }

  // Dense A_E(ω) on (E, Γ) at rank 0 from the operators applied to unit vectors.
  std::vector<int> cnt(nranks), disp(nranks, 0);
  MPI_Allgather(&nt, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    disp[r] = disp[r - 1] + cnt[r - 1];
  }
  auto dense_operator = [&](double omega)
  {
    std::unique_ptr<ComplexOperator> A2;
    {
      AssemblyDomains::Restriction restriction(space_op.GetAssemblyDomains(), environment);
      A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    }
    std::vector<std::complex<double>> A(root ? static_cast<std::size_t>(n) * n : 0);
    Vector x(nt), y(nt), yr(nt), yi(nt),
        full(root ? disp[nranks - 1] + cnt[nranks - 1] : 0), full_i(full.Size());
    auto add = [&](const ComplexOperator *X, std::complex<double> a)
    {
      for (const bool imag : {false, true})
      {
        const Operator *part = X ? (imag ? X->Imag() : X->Real()) : nullptr;
        if (part)
        {
          part->Mult(x, y);
          const std::complex<double> c = imag ? a * std::complex<double>(0.0, 1.0) : a;
          yr.Add(c.real(), y);
          yi.Add(c.imag(), y);
        }
      }
    };
    for (int j = 0; j < n; j++)
    {
      x = 0.0;
      if (idx[j] >= tstart && idx[j] < tstart + nt)
      {
        x.HostReadWrite()[idx[j] - tstart] = 1.0;
      }
      yr = 0.0;
      yi = 0.0;
      add(K.get(), 1.0);
      add(C.get(), {0.0, omega});
      add(M.get(), -omega * omega);
      add(A2.get(), 1.0);
      MPI_Gatherv(yr.HostRead(), nt, MPI_DOUBLE, full.HostWrite(), cnt.data(), disp.data(),
                  MPI_DOUBLE, 0, comm);
      MPI_Gatherv(yi.HostRead(), nt, MPI_DOUBLE, full_i.HostWrite(), cnt.data(),
                  disp.data(), MPI_DOUBLE, 0, comm);
      for (int i = 0; root && i < n; i++)
      {
        A[static_cast<std::size_t>(j) * n + i] = {full(idx[i]), full_i(idx[i])};
      }
    }
    return A;
  };

  // S_ref = A_GG - A_GE A_EE^-1 A_EG with the standard real representation [[Re, -Im],
  // [Im, Re]] of the complex blocks.
  auto dense_schur = [&](const std::vector<std::complex<double>> &A)
  {
    auto real_form = [&](int r0, int nr, int c0, int nc)
    {
      mfem::DenseMatrix T(2 * nr, 2 * nc);
      for (int i = 0; i < nr; i++)
      {
        for (int j = 0; j < nc; j++)
        {
          const auto a = A[static_cast<std::size_t>(c0 + j) * n + (r0 + i)];
          T(i, j) = T(nr + i, nc + j) = a.real();
          T(i, nc + j) = -a.imag();
          T(nr + i, j) = a.imag();
        }
      }
      return T;
    };
    mfem::DenseMatrix Aee = real_form(0, nE, 0, nE), Aeg = real_form(0, nE, nE, nG),
                      Age = real_form(nE, nG, 0, nE), Agg = real_form(nE, nG, nE, nG);
    mfem::DenseMatrixInverse inv(Aee);
    mfem::DenseMatrix X(2 * nE, 2 * nG), Y(2 * nG, 2 * nG);
    inv.Mult(Aeg, X);
    mfem::Mult(Age, X, Y);
    Agg -= Y;
    std::vector<std::complex<double>> S(static_cast<std::size_t>(nG) * nG);
    for (int b = 0; b < nG; b++)
    {
      for (int a = 0; a < nG; a++)
      {
        S[static_cast<std::size_t>(b) * nG + a] = {Agg(a, b), Agg(nG + a, b)};
      }
    }
    return S;
  };

  // Online, with the region only: given S_E and the environment's source condensation
  // (offline), the region and interface solution must be the offline one.
  DrivenSubstructure ds_online(space_op, {1}, {2}, true, nullptr, false, type);
  REQUIRE(ds_online.InterfaceSize() == nG);

  for (const double omega : {1.3, 0.7})
  {
    CAPTURE(omega);
    const auto A = dense_operator(omega);
    const auto S_ref = root ? dense_schur(A) : std::vector<std::complex<double>>();
    for (const bool region : {false, true})  // the environment alone, then both sides
    {
      ds.Condense(omega, region);
      if (!root)
      {
        continue;
      }
      const auto &S = ds.Schur();
      REQUIRE(S.size() == S_ref.size());
      double d = 0.0, m = 0.0, asym = 0.0;
      for (int b = 0; b < nG; b++)
      {
        for (int a = 0; a < nG; a++)
        {
          const std::size_t ab = static_cast<std::size_t>(b) * nG + a;
          d = std::max(d, std::abs(S[ab] - S_ref[ab]));
          m = std::max(m, std::abs(S_ref[ab]));
          asym = std::max(asym, std::abs(S[ab] - S[static_cast<std::size_t>(a) * nG + b]));
        }
      }
      CAPTURE(nG, nE, d, m, asym);
      CHECK(d <= exact_tol * m);
      CHECK(asym <= exact_tol * m);
    }

    // Substructured solves for the port excitation (an environment source) and random
    // right-hand sides (sources on both sides; the third has the second's environment
    // part): residuals of the full system.
    std::vector<ComplexVector> b(3, ComplexVector(nt)), u;
    space_op.GetExcitationVector(1, omega, b[0]);
    b[1].Real().Randomize(3 + Mpi::Rank(comm));
    b[1].Imag().Randomize(5 + Mpi::Rank(comm));
    b[2].Real().Randomize(7 + Mpi::Rank(comm));
    b[2].Imag().Randomize(11 + Mpi::Rank(comm));
    for (int i = 0; i < nt; i++)
    {
      if (!ds.RegionInterior()[i])
      {
        b[2].Real()(i) = b[1].Real()(i);
        b[2].Imag()(i) = b[1].Imag()(i);
      }
    }
    for (int d : space_op.GetNDDbcTDofLists().back())
    {
      b[1].Real()(d) = b[1].Imag()(d) = b[2].Real()(d) = b[2].Imag()(d) = 0.0;
    }
    ds.Solve({&b[0], &b[1], &b[2]}, u);
    const auto g_env = ds.EnvironmentSourceCondensation();
    const auto u_gamma = ds.InterfaceSolution();

    // A condensed environment functional (as for the environment's port voltages):
    // l^T u_k = c + h^T u_Γ,k with c fixed by the environment sources, so the value for the
    // third right-hand side follows from the second's.
    {
      ComplexVector l(nt);
      l.Real().Randomize(13 + Mpi::Rank(comm));
      l.Imag() = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (!ds.EnvironmentInterior()[i])
        {
          l.Real()(i) = 0.0;
        }
      }
      const auto h = ds.CondenseEnvironment({&l});
      CHECK(ds.CondenseEnvironment({}).empty());  // an environment without ports
      std::vector<ComplexVector> x;  // the environment interior solution for b[1]
      ds.CondenseEnvironment({&b[1]}, &x);
      double lx[2] = {mfem::InnerProduct(l.Real(), x[0].Real()),
                      mfem::InnerProduct(l.Real(), x[0].Imag())};
      Mpi::GlobalSum(2, lx, comm);
      std::complex<double> V[2];
      for (int k = 0; k < 2; k++)
      {
        double v[2] = {mfem::InnerProduct(l.Real(), u[k + 1].Real()),
                       mfem::InnerProduct(l.Real(), u[k + 1].Imag())};
        Mpi::GlobalSum(2, v, comm);
        V[k] = {v[0], v[1]};
      }
      if (root)
      {
        std::complex<double> hu[2] = {0.0, 0.0};
        for (int k = 0; k < 2; k++)
        {
          for (int a = 0; a < nG; a++)
          {
            hu[k] += h[a] * u_gamma[static_cast<std::size_t>(k + 1) * nG + a];
          }
        }
        const double dv = std::abs(V[1] - (V[0] - hu[0] + hu[1]));
        CAPTURE(dv, std::abs(V[1]));
        CHECK(dv <= exact_tol * std::abs(V[1]));
        // and l^T u = l^T x + h^T u_Γ, with x the environment interior solution.
        const double dc = std::abs(V[0] - (std::complex<double>(lx[0], lx[1]) + hu[0]));
        CAPTURE(dc);
        CHECK(dc <= exact_tol * std::abs(V[0]));
      }
    }

    // Online: the region and interface values of the offline solution, 0 in the
    // environment interior.
    {
      std::vector<std::complex<double>> S_env;
      if (root)
      {
        S_env = ds.Schur();
      }
      ds_online.Condense(omega, std::move(S_env));
      std::vector<ComplexVector> u_online;
      ds_online.Solve({&b[0], &b[1], &b[2]}, u_online, &g_env);
      for (int k = 0; k < 3; k++)
      {
        double dmax[2] = {0.0, 0.0};
        const double *xr = u_online[k].Real().HostRead(),
                     *xi = u_online[k].Imag().HostRead();
        const double *yr = u[k].Real().HostRead(), *yi = u[k].Imag().HostRead();
        for (int i = 0; i < nt; i++)
        {
          const std::complex<double> x(xr[i], xi[i]), y(yr[i], yi[i]);
          if (ds.EnvironmentInterior()[i])
          {
            dmax[1] = std::max(dmax[1], std::abs(x));
          }
          else
          {
            dmax[0] = std::max(dmax[0], std::abs(x - y));
          }
        }
        Mpi::GlobalMax(2, dmax, comm);
        const double unorm = linalg::Norml2(comm, u[k]);
        CAPTURE(k, dmax[0], dmax[1], unorm);
        CHECK(dmax[0] <= exact_tol * unorm);
        CHECK(dmax[1] == 0.0);
      }
    }
    auto K1 = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C1 = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M1 = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A21 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    auto A_full = space_op.GetSystemMatrix(
        std::complex<double>(1.0, 0.0), std::complex<double>(0.0, omega),
        std::complex<double>(-omega * omega, 0.0), K1.get(), C1.get(), M1.get(), A21.get());
    for (int k = 0; k < 3; k++)
    {
      ComplexVector r(nt);
      A_full->Mult(u[k], r);
      r -= b[k];
      const double res = linalg::Norml2(comm, r) / linalg::Norml2(comm, b[k]);
      CAPTURE(k, res);
      CHECK(res <= exact_tol);
    }
  }
}
#endif

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
  const std::string model_path =
      (temp_dir / "substruct_nc_model" / "environment.model").string();
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
            {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_nc_mag_model" / "environment.model").string();
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
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_amr_model" / "environment.model").string();
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
            {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const double *h_ind0 = ind0.Local().HostRead();
  double env_max = 0.0, reg_max = 0.0;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    (pm0.GetAttribute(e) == 2 ? env_max : reg_max) =
        std::max(pm0.GetAttribute(e) == 2 ? env_max : reg_max, h_ind0[e]);
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
    if (h_ind0[e] >= 0.3 * reg_max)
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
                 "[substructure][Serial][Parallel][Long]")
{
  // One step of online adaptive refinement of the region only, for a closed current on
  // tetrahedra: the region error indicators are zero on the environment; nonconforming
  // refinement of the marked region elements leaves the environment and the interface DOFs
  // as condensed, so the saved model still applies without factoring the environment,
  // matches a monolith on the refined mesh, and the estimated region error decreases.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path =
      (temp_dir / "substruct_mag_amr_model" / "environment.model").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
        {"Boundaries", {}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const double *h_ind0 = ind0.Local().HostRead();
  double env_max = 0.0, reg_max = 0.0;
  for (int e = 0; e < pm0.GetNE(); e++)
  {
    double &m = (pm0.GetAttribute(e) == 2) ? env_max : reg_max;
    m = std::max(m, h_ind0[e]);
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
    if (h_ind0[e] >= 0.3 * reg_max)
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
         {{"Region", {{"Attributes", {1}}}},
          {"Environment", {{"Attributes", {2}}}},
          {"SaveModel", false}}}}}};
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
  const double *h_ind = ind.Local().HostRead(), *h_ref = ref.Local().HostRead();
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
        err = std::max(err, std::abs(h_ind[e] - h_ref[r]));
        scale = std::max(scale, std::abs(h_ref[r]));
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
         {{"Region", {{"Attributes", {1}}}},
          {"Environment", {{"Attributes", {2}}}},
          {"SaveModel", false}}}}}};
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
  const double *h_ind = ind.Local().HostRead(), *h_ref = ref.Local().HostRead();
  for (int e = 0; e < pm.GetNE(); e++)
  {
    if (pm.GetAttribute(e) != 1)
    {
      env = std::max(env, std::abs(h_ind[e]));
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
    err = std::max(err, std::abs(h_ind[e] - h_ref[match]));
    scale = std::max(scale, h_ref[match]);
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
  const std::string model_path =
      (temp_dir / "substruct_redesign_model" / "environment.model").string();
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
            {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, model_path);
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
  const std::string model_path =
      (temp_dir / "substruct_closure_model" / "environment.model").string();
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
    if (mode.empty())
    {
      config["Solver"]["Substructuring"]["SaveModel"] = false;
    }
    else
    {
      test::SetModelRun(config, mode, model_path);
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
    json config = {{"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
                   {"Model", {{"Mesh", "test.msh"}}},
                   {"Domains",
                    {{"Materials",
                      {{{"Attributes", {1}}, {"Permeability", m}},
                       {{"Attributes", {2}}, {"Permeability", 4.0 * m}}}}}},
                   {"Boundaries", json::object()},
                   {"Solver",
                    {{"Order", order},
                     {"Substructuring",
                      {{"Region", {{"Attributes", {1}}}},
                       {"Environment", {{"Attributes", {2}}}},
                       {"SaveModel", false}}}}}};
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
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}}}}}},
        {"Boundaries", json::object()},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    test::SetModelRun(config, mode, path);
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
  const std::string model_path =
      (temp_dir / "substruct_mag_remesh" / "environment.model").string();
  (void)solve("Offline", model_path, a_off, nullptr);
  Vector u_ref, u_on;
  const double E_ref =
      solve("Offline", (temp_dir / "substruct_mag_ref" / "environment.model").string(),
            a_on, &u_ref);
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
