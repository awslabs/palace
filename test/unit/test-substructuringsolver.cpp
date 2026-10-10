// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/errorindicator.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "fixtures.hpp"
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/curlcurloperator.hpp"
#include "models/drivensubstructure.hpp"
#include "models/spaceoperator.hpp"
#include "models/substructuringsolver.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

using json = nlohmann::json;

namespace palace
{

namespace
{

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
std::unique_ptr<mfem::ParMesh> MakeWavyTetSplit(int nx, bool reorder = false)
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
std::unique_ptr<mfem::ParMesh> MakeGradedSplit(int a, int b, int n,
                                               const RegionDesign &design = {},
                                               bool nonconforming = false)
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
  if (nonconforming)
  {
    serial.EnsureNCMesh();  // a ParMesh cannot be made nonconforming afterwards
  }
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
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
    const double e_sub = ss.ElectrostaticEnergy(u);
    const double e_ref = ss.ElectrostaticEnergy(u_full);
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
  // Gauge-free H(curl) magnetostatics: pure curl-curl (singular). Fields are
  // gauge-ambiguous, but the magnetic energy is gauge-invariant. A divergence-free
  // (rotational) current source is a valid magnetostatic excitation; the region-condensed
  // energy must match a monolith AMS-singular (min-norm) solve of the identical pure
  // curl-curl operator.
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

  // Divergence-free rotational current source on the parent ND space.
  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  auto jfun = [](const mfem::Vector &x, mfem::Vector &j)
  {
    j.SetSize(3);
    j(0) = -(x(1) - 0.5);
    j(1) = (x(0) - 0.5);
    j(2) = 0.0;
  };
  mfem::VectorFunctionCoefficient jc(3, jfun);
  mfem::ParLinearForm lf(&pfes);
  lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(jc));
  lf.Assemble();
  Vector f(pfes.GetTrueVSize());
  lf.ParallelAssemble(f);

  Vector u = ss.SolveSource(f);
  const double e_sub = ss.ElectrostaticEnergy(u);  // 1/2 u^T K_curlcurl u (magnetic energy)

  // Monolith reference with the identical approach: solve curl-curl + small mass (definite,
  // AMS-convergent), measure energy with the pure curl-curl operator.
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector nu_by_attr(max_attr), mass_by_attr(max_attr);
  nu_by_attr = 0.0;
  mass_by_attr = 0.0;
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  mass_by_attr(0) = 1.0e-3;  // match SubstructuringSolver::Impl::kMagRegularization
  mass_by_attr(1) = 1.0e-3;
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

  CHECK(e_sub > 1.0e-6);  // nontrivial magnetic energy
  CHECK(std::abs(e_sub - e_mono) <=
        1.0e-6 * std::abs(e_mono));  // substructuring == monolith (same operator)
}

TEST_CASE("SubstructuringSolver magnetostatic Dirichlet excitation",
          "[substructure][Serial][Parallel]")
{
  // A magnetostatic Dirichlet-lift excitation (prescribed tangential A on a boundary) gives
  // the same gauge-invariant energy region-condensed and monolithic.
  const int order = 1;
  const double mu_r = 1.0, mu_e = 4.0;
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
          {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
      {"Boundaries",
       {{"Terminal",
         {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
      {"Solver",
       {{"Order", order},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  Vector u = ss.SolveExcitation(1);  // tangential A = 1 on attr-1 boundary, 0 on attr-2
  const double e_sub = ss.ElectrostaticEnergy(u);

  // Monolith: same Dirichlet data, curl-curl + small mass (definite), pure-curl-curl
  // energy.
  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector nu_by_attr(max_attr), mass_by_attr(max_attr);
  nu_by_attr = 0.0;
  mass_by_attr = 0.0;
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  mass_by_attr(0) = 1.0e-3;  // match kMagRegularization
  mass_by_attr(1) = 1.0e-3;
  mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
  mfem::ParBilinearForm asolve(&pfes);
  asolve.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
  asolve.Assemble();
  mfem::ParBilinearForm aenergy(&pfes);
  aenergy.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  aenergy.Assemble();
  aenergy.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Kpure(aenergy.ParallelAssemble());

  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> ess_bdr(maxb), ess_tdofs, t1_bdr(maxb), t1_tdofs;
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  t1_bdr = 0;
  t1_bdr[0] = 1;
  pfes.GetEssentialTrueDofs(t1_bdr, t1_tdofs);
  mfem::ParGridFunction xgf(&pfes);
  {
    mfem::Vector td(pfes.GetTrueVSize());
    td = 0.0;
    for (int i = 0; i < t1_tdofs.Size(); i++)
    {
      td(t1_tdofs[i]) = 1.0;
    }
    xgf.SetFromTrueDofs(td);
  }
  mfem::ParLinearForm bform(&pfes);
  bform = 0.0;
  bform.Assemble();
  mfem::OperatorPtr A;
  mfem::Vector B, X;
  asolve.FormLinearSystem(ess_tdofs, xgf, bform, A, X, B);
  mfem::HypreParMatrix *Ah = A.As<mfem::HypreParMatrix>();
  mfem::HypreAMS ams(*Ah, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ah);
  pcg.SetTol(1e-13);
  pcg.SetMaxIter(2000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  pcg.Mult(B, X);
  asolve.RecoverFEMSolution(X, bform, xgf);
  mfem::Vector u_full;
  xgf.GetTrueDofs(u_full);
  mfem::Vector t(pfes.GetTrueVSize());
  Kpure->Mult(u_full, t);
  double local = 0.0;
  for (int i = 0; i < pfes.GetTrueVSize(); i++)
  {
    local += u_full(i) * t(i);
  }
  double e_mono = 0.0;
  MPI_Allreduce(&local, &e_mono, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
  e_mono *= 0.5;

  CHECK(e_sub > 1.0e-8);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-6 * std::abs(e_mono));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver offline/online model reuse",
                 "[substructure][Serial][Parallel]")
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
                 "[substructure][Serial][Parallel]")
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
  // So is a changed boundary condition on the environment (its terminal removed).
  IoData no_env_terminal = make_config("Online", 1.0, 10.0);
  no_env_terminal.boundaries.terminal.erase(2);
  CHECK_THROWS_WITH(condense(no_env_terminal, 3, 4),
                    Catch::Matchers::ContainsSubstring("The environment differs"));
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
    std::vector<Vector> u = ss.SolveExcitations(terms);
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
                 "SubstructuringSolver environment-free magnetostatic energies",
                 "[substructure][Serial][Parallel]")
{
  // Magnetostatics: the energy operator (pure curl-curl K) differs from the solve operator
  // (curl-curl + eps mass, A). The condensed energy uses S^K = S_E - P^T D P (D = A - K)
  // and the D-corrected lift couplings. An offline run that saves a model must reproduce
  // the full-field energies u_i^T K u_j, and the online run must reproduce them without
  // factoring the environment (lifts matched by id + fingerprint). Tolerance: the synthetic
  // unit lifts are mostly curl-free, so the order-1 energies are tiny (~1e-6) and the full-
  // field reference itself is only good to ~5e-9 (its own asymmetry); omitting the
  // eps-correction errs by orders of magnitude more.
  const int order = GENERATE(1, 2);
  const bool tet = GENERATE(false, true);
  CAPTURE(order, tet);
  const std::string model_path = (temp_dir / "substruct_mag_energy_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  auto make_mesh = [tet]() { return tet ? MakeWavyTetSplit(4) : MakeSplitCube(4); };
  const std::vector<int> ids = {1, 2};
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
  auto full_field_E = [&ids](SubstructuringSolver &ss)
  {
    std::vector<Vector> lifts = {ss.TerminalLift(1), ss.TerminalLift(2)};
    std::vector<Vector> u = ss.SolveDirichlets(lifts);
    mfem::DenseMatrix E(2);
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        E(i, j) = ss.MutualEnergy(u[i], u[j]);
      }
    }
    return E;
  };

  IoData iodata_off = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(make_mesh()));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const std::vector<Vector> lifts_off = {off.TerminalLift(1), off.TerminalLift(2)};
  const mfem::DenseMatrix E_off = off.EnergyMatrix(ids, lifts_off);
  const mfem::DenseMatrix E_ref = full_field_E(off);
  CHECK(rel_diff(E_off, E_ref) <= 1.0e-7);

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(make_mesh()));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const std::vector<Vector> lifts_on = {on.TerminalLift(1), on.TerminalLift(2)};
  std::vector<Vector> fields;
  const mfem::DenseMatrix E_on = on.EnergyMatrix(ids, lifts_on);
  CHECK_FALSE(on.EnvironmentFactored());
  // Same arithmetic up to the run-to-run rounding of the parallel direct solvers, which the
  // tiny order-1 energies amplify to ~1e-8 (as for the full-field reference above).
  CHECK(rel_diff(E_on, E_off) <= 5.0e-8);
  // Recover a field on demand; it matches the full solve.
  (void)on.EnergyMatrix(ids, lifts_on, &fields, 1);
  CHECK(on.EnvironmentFactored());
  const Vector u_full = on.SolveDirichlet(lifts_on[0]);
  Vector d(fields[0]);
  d -= u_full;
  const double un = std::sqrt(mfem::InnerProduct(Mpi::World(), u_full, u_full));
  CHECK(std::sqrt(mfem::InnerProduct(Mpi::World(), d, d)) <= 1.0e-8 * un);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver magnetostatic energies across face orientations",
                 "[substructure][Serial][Parallel]")
{
  // A model saved on tetrahedra, reused online on the same tetrahedra in another element
  // and vertex order: other face orientations (as another partition gives), so second-order
  // Nédélec face DOFs on Gamma have another basis, not a signed permutation of the saved
  // one. The lifts are the same fields on both meshes, projected on all DOFs so that they
  // reach Gamma (and their fingerprints see it); the online energies need no environment
  // solve and equal the offline ones.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const std::string model_path = (temp_dir / "substruct_mag_orient_model.bin").string();
  auto make_config = [order, &model_path](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", model_path}}}}}};
    return IoData(config, false);
  };
  const std::vector<int> ids = {1, 2};
  auto lifts = [order, &ids](Mesh &mesh)
  {
    mfem::ND_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace fes(&mesh.Get(), &fec);
    std::vector<Vector> x;
    for (int id : ids)
    {
      mfem::VectorFunctionCoefficient f(3,
                                        [id](const mfem::Vector &p, mfem::Vector &v)
                                        {
                                          v.SetSize(3);
                                          const double k = 3.0 * id;
                                          v(0) = std::sin(k * p(1)) * p(2);
                                          v(1) = std::cos(k * p(2)) * p(0);
                                          v(2) = std::sin(k * p(0) + p(1));
                                        });
      mfem::ParGridFunction gf(&fes);
      gf.ProjectCoefficient(f);
      Vector t(fes.GetTrueVSize());
      gf.GetTrueDofs(t);
      x.push_back(std::move(t));
    }
    return x;
  };

  IoData iodata_off = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  const mfem::DenseMatrix E_off = off.EnergyMatrix(ids, lifts(*mesh_off.back()));

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4, true)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix E_on = on.EnergyMatrix(ids, lifts(*mesh_on.back()));
  CHECK_FALSE(on.EnvironmentFactored());
  double d = 0.0, m = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      d = std::max(d, std::abs(E_on(i, j) - E_off(i, j)));
      m = std::max(m, std::abs(E_off(i, j)));
    }
  }
  CAPTURE(d, m);
  CHECK(d <= 1.0e-9 * m);
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
  const mfem::DenseMatrix E_off = off.CurrentEnergyMatrix({1}, current(*mesh_off.back()));

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeWavyTetSplit(4, true)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix E_on = on.CurrentEnergyMatrix({1}, current(*mesh_on.back()));
  CHECK_FALSE(on.EnvironmentFactored());
  CAPTURE(E_on(0, 0), E_off(0, 0));
  CHECK(std::abs(E_on(0, 0) - E_off(0, 0)) <= 1.0e-9 * std::abs(E_off(0, 0)));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver magnetostatic energies on a re-meshed region",
                 "[substructure][Serial][Parallel]")
{
  // Offline on one region mesh, online on a re-meshed region (fixed Gamma + environment,
  // signed H(curl) re-ordering of S_E, S^K and the lift couplings): the online condensed
  // energies, computed without an environment solve, must equal the online full-field
  // energies.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  auto make_config = [this, order](const std::string &mode)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", 1.0}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", 4.0}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", order},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"Mode", mode},
            {"SaveModel", (temp_dir / "substruct_mag_remesh_energy.bin").string()}}}}}};
    return IoData(config, false);
  };
  const std::vector<int> ids = {1, 2};
  IoData iodata_off = make_config("Offline");
  std::vector<std::unique_ptr<Mesh>> mesh_off;
  mesh_off.push_back(std::make_unique<Mesh>(MakeGradedSplit(3, 4, 4)));
  SubstructuringSolver off(iodata_off, mesh_off);
  off.CondenseEnvironment();
  (void)off.EnergyMatrix(ids, {off.TerminalLift(1), off.TerminalLift(2)});

  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(5, 4, 4)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const std::vector<Vector> lifts = {on.TerminalLift(1), on.TerminalLift(2)};
  const mfem::DenseMatrix E_on = on.EnergyMatrix(ids, lifts);
  CHECK_FALSE(on.EnvironmentFactored());
  std::vector<Vector> u = on.SolveDirichlets(lifts);
  double d = 0.0, m = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      const double e = on.MutualEnergy(u[i], u[j]);
      d = std::max(d, std::abs(E_on(i, j) - e));
      m = std::max(m, std::abs(e));
    }
  }
  CHECK(d <= 1.0e-7 * m);  // see the tolerance note above
}

TEST_CASE("SubstructuringSolver block low-rank environment factorization",
          "[substructure][Serial][Parallel]")
{
  // Solver.Substructuring.FactorizationTol > 0 selects a block low-rank (BLR) MUMPS
  // environment factorization (the matrix is scaled so the tolerance is relative; the Schur
  // and the solves are rescaled). The capacitance matrix must stay within ~the tolerance of
  // the exact factorization. Without MUMPS the tolerance is ignored (exact).
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  auto capacitance = [order](double tol)
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
            {"FactorizationTol", tol}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(8)));
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
  CHECK(d <= 1.0e-6 * m);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver capacitance from a model without terminal modes",
                 "[substructure][Serial][Parallel]")
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

  // Strip the terminal-mode section: keep the header (2 ints), signature (nG x 3), S_E (nG
  // x nG) and environment fingerprint (1 int + 5 doubles).
  if (Mpi::Root(Mpi::World()))
  {
    int nG = 0;
    {
      std::ifstream f(model_path, std::ios::binary);
      f.read(reinterpret_cast<char *>(&nG), sizeof(int));
    }
    const auto kept = static_cast<std::uintmax_t>(3 * sizeof(int)) +
                      static_cast<std::uintmax_t>(sizeof(double)) * (nG * 3 + nG * nG + 5);
    REQUIRE(std::filesystem::file_size(model_path) > kept);
    std::filesystem::resize_file(model_path, kept);
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

TEST_CASE("SubstructuringSolver magnetostatic inductance matrix",
          "[substructure][Serial][Parallel]")
{
  // Inductance matrix M = R^-1 of two Dirichlet-lift excitations, with the reluctance
  // R_ij = A_i^T K A_j / (Phi_i Phi_j) (pure curl-curl energy): the region-condensed matrix
  // must match a monolith and be symmetric.
  const int order = GENERATE(1, 2);
  const bool tet = GENERATE(false, true);
  CAPTURE(order);
  CAPTURE(tet);
  const double mu_r = 1.0, mu_e = 4.0;
  json config = {
      {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
      {"Model", {{"Mesh", "test.msh"}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
          {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
      {"Boundaries",
       {{"Terminal",
         {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
      {"Solver",
       {{"Order", order},
        {"Substructuring",
         {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
  IoData iodata(config, false);

  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(tet ? MakeWavyTetSplit(6) : MakeSplitCube(6)));
  SubstructuringSolver ss(iodata, mesh);
  ss.CondenseEnvironment();
  std::vector<Vector> As = {ss.SolveExcitation(1), ss.SolveExcitation(2)};
  auto invert2 = [](mfem::DenseMatrix &R)
  {
    mfem::DenseMatrix M(R);
    M.Invert();
    return M;
  };
  mfem::DenseMatrix R_sub(2);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      R_sub(i, j) = ss.MutualEnergy(As[i], As[j]);  // Phi = 1
    }
  }
  mfem::DenseMatrix M_sub = invert2(R_sub);
  // Symmetric up to iterative-solver residual (cross-energies are analytically symmetric).
  CHECK(std::abs(M_sub(0, 1) - M_sub(1, 0)) <=
        1.0e-5 * std::max(std::abs(M_sub(0, 0)), std::abs(M_sub(1, 1))));

  // Monolith: same two Dirichlet excitations, curl-curl + small mass, pure-curl-curl
  // cross-energies, then reluctance inversion.
  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector nu_by_attr(max_attr), mass_by_attr(max_attr);
  nu_by_attr = 0.0;
  mass_by_attr = 0.0;
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  mass_by_attr(0) = 1.0e-3;
  mass_by_attr(1) = 1.0e-3;
  mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
  mfem::ParBilinearForm asolve(&pfes);
  asolve.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
  asolve.Assemble();
  mfem::ParBilinearForm aenergy(&pfes);
  aenergy.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  aenergy.Assemble();
  aenergy.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Kpure(aenergy.ParallelAssemble());
  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> ess_bdr(maxb), ess_tdofs;
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  auto solve_dir = [&](int drive_attr)
  {
    mfem::Array<int> d_bdr(maxb), d_tdofs;
    d_bdr = 0;
    d_bdr[drive_attr - 1] = 1;
    pfes.GetEssentialTrueDofs(d_bdr, d_tdofs);
    mfem::ParGridFunction xgf(&pfes);
    mfem::Vector td(pfes.GetTrueVSize());
    td = 0.0;
    for (int i = 0; i < d_tdofs.Size(); i++)
    {
      td(d_tdofs[i]) = 1.0;
    }
    xgf.SetFromTrueDofs(td);
    mfem::ParLinearForm bform(&pfes);
    bform = 0.0;
    bform.Assemble();
    mfem::OperatorPtr A;
    mfem::Vector B, X;
    asolve.FormLinearSystem(ess_tdofs, xgf, bform, A, X, B);
    mfem::HypreParMatrix *Ah = A.As<mfem::HypreParMatrix>();
    mfem::HypreAMS amsp(*Ah, &pfes);
    amsp.SetPrintLevel(0);
    mfem::HyprePCG pcg(*Ah);
    pcg.SetTol(1e-13);
    pcg.SetMaxIter(2000);
    pcg.SetPrintLevel(0);
    pcg.SetPreconditioner(amsp);
    pcg.Mult(B, X);
    asolve.RecoverFEMSolution(X, bform, xgf);
    mfem::Vector uu;
    xgf.GetTrueDofs(uu);
    return uu;
  };
  std::vector<mfem::Vector> Am = {solve_dir(1), solve_dir(2)};
  auto cross = [&](const mfem::Vector &a, const mfem::Vector &b)
  {
    mfem::Vector t(pfes.GetTrueVSize());
    Kpure->Mult(b, t);
    double loc = 0.0;
    for (int i = 0; i < pfes.GetTrueVSize(); i++)
    {
      loc += a(i) * t(i);
    }
    double g = 0.0;
    MPI_Allreduce(&loc, &g, 1, MPI_DOUBLE, MPI_SUM, Mpi::World());
    return g;
  };
  mfem::DenseMatrix R_mono(2);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      R_mono(i, j) = cross(Am[i], Am[j]);
    }
  }
  mfem::DenseMatrix M_mono = invert2(R_mono);
  // Scale the comparison by the largest inductance entry (M(0,0) is tiny at higher order,
  // so scaling by it alone would make the tolerance meaningless for the other entries).
  double scale = 0.0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      scale = std::max(scale, std::abs(M_mono(i, j)));
    }
  }
  // Known limitation: order >= 2 H(curl) on tetrahedra is exact serially but wrong in
  // parallel (the higher-order tetrahedral edge/face DOF orientation across a partition cut
  // is mishandled in the interface identification -- order-1 tets and order-2 hexes are
  // fine). Skip that combination in parallel until the signed higher-order simplex map is
  // fixed.
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(std::abs(M_sub(i, j) - M_mono(i, j)) <= 1.0e-5 * scale);
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver cross-run region re-meshing",
                 "[substructure][Serial][Parallel]")
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
  const double e_sub = on.ElectrostaticEnergy(u_on);

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
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-6 * std::abs(e_mono));
}

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

TEST_CASE_METHOD(palace::test::SharedTempDir, "SubstructuringSolver London sheet energies",
                 "[substructure][Serial][Parallel]")
{
  // Magnetostatics with a London superconductor sheet (finite penetration depth) on an
  // interior plane crossing the interface: region and environment sheet faces, meeting
  // Gamma along their edges. Two sheet generators a_k drive the sources M_sheet a_k with
  // the sheet DOFs free; the energy matrix E_ij = u_i^T K_cc u_j + (u_i - a_i)^T M_sheet
  // (u_j - a_j) must match a monolith (curl-curl + sheet, with a refined regularization),
  // offline and after reloading the model.
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  const double mu_r = 1.0, mu_e = 4.0, lambda = 0.2, thickness = 0.05;
  const std::string model_path = (temp_dir / "substruct_london_model.bin").string();
  auto make_config = [&](const std::string &mode)
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
                       {"SaveModel", model_path}}}}}};
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
  REQUIRE(ss.HasSheets());
  const mfem::DenseMatrix E = ss.SheetEnergyMatrix({1, 2}, a);
  // The full-field path (fields requested) gives the same energies.
  std::vector<Vector> fields;
  const mfem::DenseMatrix E_full = ss.SheetEnergyMatrix({1, 2}, a, &fields, 2);

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

  // Online reuse of the saved model.
  IoData iodata_on = make_config("Online");
  std::vector<std::unique_ptr<Mesh>> mesh_on;
  mesh_on.push_back(std::make_unique<Mesh>(MakeGradedSplit(4, 5, 6, design)));
  SubstructuringSolver on(iodata_on, mesh_on);
  on.CondenseEnvironment();
  const mfem::DenseMatrix E_on = on.SheetEnergyMatrix({1, 2}, a);
  CHECK_FALSE(on.EnvironmentFactored());
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
  on.SheetEnergyMatrix({1, 2}, a, &fields_on, 2);
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
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver surface-current energies",
                 "[substructure][Serial][Parallel]")
{
  // Magnetostatics in a box with PEC walls, driven by surface currents on an interior plane
  // crossing the interface (from the x = 0 wall to the x = 1 wall, returning through the
  // walls): region and environment parts of the excitation. The energies must match a
  // monolith, offline and after reloading the model (also with an added port in the
  // region), with and without a London sheet on the plane; a current that does not close is
  // rejected.
  const int order = GENERATE(1, 2);
  const bool london = GENERATE(false, true);
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
  CHECK(ss.HasSheets() == london);
  const mfem::DenseMatrix E = ss.CurrentEnergyMatrix({1, 2}, J12);
  std::vector<Vector> fields;
  const mfem::DenseMatrix E_full = ss.CurrentEnergyMatrix({1, 2}, J12, &fields, 2);

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
  const mfem::DenseMatrix E_on = on.CurrentEnergyMatrix({1, 2}, J12);
  const mfem::DenseMatrix E_on3 = on.CurrentEnergyMatrix({1, 2, 3}, J123);
  CHECK_FALSE(on.EnvironmentFactored());
  CHECK(max_diff(E_on, E) <= 1.0e-9 * m);
  CAPTURE(E_on3(2, 2), E_mono(2, 2), E_on3(0, 2), E_mono(0, 2));
  CHECK(max_diff(E_on3, E_mono) <= tol * m);

  // Online fields agree with the full-field path in the energy norm (they differ only in
  // the zero-energy gauge).
  std::vector<Vector> fields_on;
  on.CurrentEnergyMatrix({1, 2}, J12, &fields_on, 2);
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
  const mfem::DenseMatrix E_rev = on.CurrentEnergyMatrix({1}, J_rev);
  CHECK(std::abs(E_rev(0, 0) - E(0, 0)) <= 1.0e-9 * m);

  // A current ending on the free sheet does not close (with PEC only).
  if (london)
  {
    CHECK_NOTHROW(on.CurrentEnergyMatrix({4}, {J[3]}));
  }
#if defined(MFEM_USE_EXCEPTIONS)
  else
  {
    CHECK_THROWS_WITH(on.CurrentEnergyMatrix({4}, {J[3]}),
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

TEST_CASE("SpaceOperator assembly restricted to region and environment",
          "[substructure][Serial][Parallel]")
{
  // The driven operators assembled on the region and on the environment add up to the
  // full operators, and neither side couples to the other side's interior DOFs: lumped
  // ports on both sides, a lossy dielectric region, a conducting environment, a
  // second-order absorbing boundary on both sides and an impedance sheet crossing Gamma.
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
    for (int i = 0; i < mt.Size(); i++)
    {
      marks[side][i] = (mt(i) > 0.0);
    }
  }

  const std::vector<int> region = {1}, environment = {2};
  const double omega = 1.3;  // nondimensional (the configuration is not nondimensionalized)
  auto assemble = [&](const std::vector<int> *domains)
  {
    std::optional<SpaceOperator::AssemblyRestriction> restriction;
    if (domains)
    {
      restriction.emplace(space_op, *domains);
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
  // relative to ‖S‖_F, for any scaling of S.
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
  }
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

#if defined(MFEM_USE_MUMPS)
TEST_CASE("MumpsSchurSolver factors on a subset of the ranks",
          "[substructure][Serial][Parallel]")
{
  // A complex symmetric shifted Laplacian with Schur variables on a face of the cube,
  // factored on all ranks and on fewer: the same Schur complement and internal solves,
  // also after a refactorization with new values.
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

  struct Result
  {
    std::vector<std::complex<double>> S[2];
    ComplexVector y[2];
  };
  auto run = [&](int procs)
  {
    ComplexMumpsSchurSolver::Coo coo;
    std::vector<std::vector<double>> vals;
    std::vector<int> row_ptr;
    LowerTrianglePattern({A.get()}, coo.irn, coo.jcn, vals, row_ptr);
    for (double v : vals[0])
    {
      coo.val.push_back(std::complex<double>(1.0, 0.1) * v);
    }
    ComplexMumpsSchurSolver solver(comm, A->GetGlobalNumRows(), nt, std::move(coo),
                                   schur_vars, 0.0, procs, true, false);
    Result res;
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
      res.S[k] = solver.Schur();
      solver.SolveInternal({&x}, {&res.y[k]});
    }
    return res;
  };
  auto dist = [](const std::vector<std::complex<double>> &a,
                 const std::vector<std::complex<double>> &b, std::complex<double> c)
  {
    double d = 0.0, m = 0.0;
    for (std::size_t q = 0; q < b.size(); q++)
    {
      d = std::max(d, std::abs(a[q] - c * b[q]));
      m = std::max(m, std::abs(c * b[q]));
    }
    return d / m;
  };
  const Result ref = run(0);
  if (Mpi::Root(comm))
  {
    CHECK(dist(ref.S[1], ref.S[0], 2.0) <= 1.0e-12);
  }
  for (int procs : {1, 2})
  {
    CAPTURE(procs, nranks);
    const Result res = run(procs);
    for (int k = 0; k < 2; k++)
    {
      if (Mpi::Root(comm))
      {
        CHECK(dist(res.S[k], ref.S[k], 1.0) <= 1.0e-12);
      }
      ComplexVector d(res.y[k]);
      d -= ref.y[k];
      CHECK(linalg::Norml2(comm, d) <= 1.0e-12 * linalg::Norml2(comm, ref.y[k]));
    }
  }
}

TEST_CASE("DrivenSubstructure condenses the environment exactly",
          "[substructure][Serial][Parallel]")
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
  const bool tet = (kind == MeshKind::TET);
  CAPTURE(order, static_cast<int>(kind));
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
  DrivenSubstructure ds(space_op, {1}, {2});
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
    SpaceOperator::AssemblyRestriction restriction(space_op, environment);
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
      SpaceOperator::AssemblyRestriction restriction(space_op, environment);
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
        x(static_cast<int>(idx[j] - tstart)) = 1.0;
      }
      yr = 0.0;
      yi = 0.0;
      add(K.get(), 1.0);
      add(C.get(), {0.0, omega});
      add(M.get(), -omega * omega);
      add(A2.get(), 1.0);
      MPI_Gatherv(yr.GetData(), nt, MPI_DOUBLE, full.GetData(), cnt.data(), disp.data(),
                  MPI_DOUBLE, 0, comm);
      MPI_Gatherv(yi.GetData(), nt, MPI_DOUBLE, full_i.GetData(), cnt.data(), disp.data(),
                  MPI_DOUBLE, 0, comm);
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
  DrivenSubstructure ds_online(space_op, {1}, {2}, true);
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
      CHECK(d <= 1.0e-10 * m);
      CHECK(asym <= 1.0e-10 * m);
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
        CHECK(dv <= 1.0e-10 * std::abs(V[1]));
        // and l^T u = l^T x + h^T u_Γ, with x the environment interior solution.
        const double dc = std::abs(V[0] - (std::complex<double>(lx[0], lx[1]) + hu[0]));
        CAPTURE(dc);
        CHECK(dc <= 1.0e-10 * std::abs(V[0]));
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
        for (int i = 0; i < nt; i++)
        {
          const std::complex<double> x(u_online[k].Real()(i), u_online[k].Imag()(i)),
              y(u[k].Real()(i), u[k].Imag()(i));
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
        CHECK(dmax[0] <= 1.0e-10 * unorm);
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
      CHECK(res <= 1.0e-10);
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

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SubstructuringSolver magnetostatic energy on a nonconforming mesh",
                 "[substructure][Serial][Parallel]")
{
  // H(curl) on an adapted (nonconforming) mesh with hanging edges on the interface: the
  // condensed magnetic energy of a divergence-free current source must match a monolith of
  // the same operator on the same mesh.
  const int order = GENERATE(1, 2);
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
  auto jfun = [](const mfem::Vector &x, mfem::Vector &j)
  {
    j.SetSize(3);
    j(0) = -(x(1) - 0.5);
    j(1) = (x(0) - 0.5);
    j(2) = 0.0;
  };
  mfem::VectorFunctionCoefficient jc(3, jfun);
  mfem::ParLinearForm lf(&pfes);
  lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(jc));
  lf.Assemble();
  Vector f(pfes.GetTrueVSize());
  lf.ParallelAssemble(f);
  const double e_sub = ss.ElectrostaticEnergy(ss.SolveSource(f));

  mfem::Vector nu_by_attr(2), mass_by_attr(2);
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  mass_by_attr = 1.0e-3;  // SubstructuringSolver's regularization
  mfem::PWConstCoefficient nu(nu_by_attr), massc(mass_by_attr);
  mfem::ParBilinearForm asolve(&pfes), aenergy(&pfes);
  asolve.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(massc));
  aenergy.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  asolve.Assemble();
  asolve.Finalize();
  aenergy.Assemble();
  aenergy.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> Ksolve(asolve.ParallelAssemble());
  std::unique_ptr<mfem::HypreParMatrix> Kpure(aenergy.ParallelAssemble());
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
  Kpure->Mult(u, t);
  const double e_mono = 0.5 * mfem::InnerProduct(pmesh.GetComm(), u, t);
  CAPTURE(e_sub, e_mono);
  CHECK(e_sub > 1.0e-6);
  CHECK(std::abs(e_sub - e_mono) <= 1.0e-6 * std::abs(e_mono));

  // The refinement changes the energy by far more than the mismatch above, so the match is
  // meaningful.
  {
    auto coarse = MakeGradedSplit(4, 5, 6);
    mfem::ParFiniteElementSpace cfes(coarse.get(), &fec);
    mfem::ParLinearForm clf(&cfes);
    clf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(jc));
    clf.Assemble();
    Vector cf(cfes.GetTrueVSize());
    clf.ParallelAssemble(cf);
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
  const double e_on = on.ElectrostaticEnergy(on.SolveSource(f));
  CAPTURE(e_on);
  CHECK(std::abs(e_on - e_sub) <= 1.0e-7 * std::abs(e_sub));  // iterative region solves
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
  const ErrorIndicator ind0 = on0.RegionErrorIndicator(region_fields, C0);
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
  const ErrorIndicator ind1 = on1.RegionErrorIndicator(region_fields, C1);
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
                 "SubstructuringSolver geometric redesign of the region",
                 "[substructure][Serial][Parallel]")
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
                 "SubstructuringSolver magnetostatic cross-run re-meshing",
                 "[substructure][Serial][Parallel]")
{
  // H(curl) region re-meshing: the environment (attr 2, b=5) is identical between a coarse-
  // region offline run and a re-meshed (finer) online run, so the loaded+reordered S_E (via
  // signed edge-signature matching) must reproduce a freshly materialized S_E on the online
  // mesh. Compare the recovered field of an Online (loaded) solve to an Offline (fresh)
  // solve on the same re-meshed mesh.
  const auto res = GENERATE(std::make_pair(4, 8), std::make_pair(6, 6));
  const int a_off = res.first, a_on = res.second;
  const int order = GENERATE(1, 2);
  CAPTURE(a_off, a_on, order);
  const double mu_r = 1.0, mu_e = 4.0;
  auto make_config = [&](const std::string &mode, const std::string &path)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
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

  const std::string model_path = (temp_dir / "substruct_mag_remesh.bin").string();
  // Offline on the coarse region: materialize + save S_E (signed edge signature).
  {
    IoData io = make_config("Offline", model_path);
    std::vector<std::unique_ptr<Mesh>> m;
    m.push_back(std::make_unique<Mesh>(MakeGradedSplit(a_off, 5, 6)));
    SubstructuringSolver off(io, m);
    off.CondenseEnvironment();
    (void)off.SolveExcitation(1);
  }
  // Fresh reference on the re-meshed online mesh (materialize S_E on the same environment).
  Vector ref;
  {
    IoData io = make_config("Offline", (temp_dir / "substruct_mag_ref.bin").string());
    std::vector<std::unique_ptr<Mesh>> m;
    m.push_back(std::make_unique<Mesh>(MakeGradedSplit(a_on, 5, 6)));
    SubstructuringSolver s(io, m);
    s.CondenseEnvironment();
    ref = s.SolveExcitation(1);
  }
  // Online on the re-meshed mesh: load the coarse-region model and reorder S_E onto the new
  // interface edge DOFs (with orientation signs).
  Vector got;
  {
    IoData io = make_config("Online", model_path);
    std::vector<std::unique_ptr<Mesh>> m;
    m.push_back(std::make_unique<Mesh>(MakeGradedSplit(a_on, 5, 6)));
    SubstructuringSolver s(io, m);
    s.CondenseEnvironment();
    got = s.SolveExcitation(1);
  }
  Vector d(got);
  d -= ref;
  CHECK(ref.Norml2() > 1.0e-30);
  // The offline and online runs materialize S_E independently (different region meshes/
  // partitions) with iterative solves, so the fields agree to ~solver tolerance, not to
  // machine precision.
  CHECK(d.Norml2() <= 1.0e-5 * (ref.Norml2() + 1.0e-30));
}

TEST_CASE("SubstructuringSolver HODLR off-diagonal compression accuracy vs tolerance",
          "[substructure][Serial]")
{
  // HODLR compression of S_E: the region energy converges to the dense one as the
  // compression tolerance decreases.
  auto energy_at = [](double tol)
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
            {"InterfaceOffdiagTol", tol}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(10)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.ElectrostaticEnergy(ss.SolveExcitation(1));
  };
  const double e_exact = energy_at(0.0);
  CHECK(e_exact > 1.0e-12);
  double relerr = 1.0;
  for (double tol : {1e-2, 1e-4, 1e-6, 1e-8})
  {
    relerr = std::abs(energy_at(tol) - e_exact) / std::abs(e_exact);
  }
  // Tightest tolerance should recover the exact region energy closely.
  CHECK(relerr <= 1.0e-6);
}

TEST_CASE("SubstructuringSolver HODLR compressed apply matches dense (parallel)",
          "[substructure][Parallel]")
{
  // Exercises the parallel build (gather -> build on rank 0 -> broadcast) and the
  // replicated compressed apply: the region energy with a tight-tolerance compressed S_E
  // must match the dense S_E at any process count.
  auto energy_at = [](double tol)
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
            {"InterfaceOffdiagTol", tol}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(8)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    return ss.ElectrostaticEnergy(ss.SolveExcitation(1));
  };
  const double e_dense = energy_at(0.0);
  const double e_comp = energy_at(1.0e-8);
  CHECK(e_dense > 1.0e-12);
  CHECK(std::abs(e_comp - e_dense) / std::abs(e_dense) <= 1.0e-6);
}

TEST_CASE("SubstructuringSolver HODLR magnetostatic compressed apply matches dense",
          "[substructure][Parallel]")
{
  // HODLR compression on a Nedelec (H(curl)) interface: the clustering coordinate is the
  // edge midpoint. The region reluctance (self-energy) with a tight-tolerance compressed
  // S_E must match the dense operator at any process count.
  const double mu_r = 1.0, mu_e = 4.0;
  auto energy_at = [&](double tol)
  {
    json config = {
        {"Problem", {{"Type", "Magnetostatic"}, {"Output", "test_output"}}},
        {"Model", {{"Mesh", "test.msh"}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permeability", mu_r}, {"Permittivity", 1.0}},
            {{"Attributes", {2}}, {"Permeability", mu_e}, {"Permittivity", 1.0}}}}}},
        {"Boundaries",
         {{"Terminal",
           {{{"Index", 1}, {"Attributes", {1}}}, {{"Index", 2}, {"Attributes", {2}}}}}}},
        {"Solver",
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}},
            {"Environment", {{"Attributes", {2}}}},
            {"InterfaceOffdiagTol", tol}}}}}};
    IoData iodata(config, false);
    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    const Vector a = ss.SolveExcitation(1);
    return ss.MutualEnergy(a, a);
  };
  const double e_dense = energy_at(0.0);
  const double e_comp = energy_at(1.0e-8);
  const double relerr = std::abs(e_comp - e_dense) / std::abs(e_dense);
  CHECK(std::abs(e_dense) > 1.0e-12);
  CHECK(relerr <= 1.0e-5);
}

}  // namespace palace
