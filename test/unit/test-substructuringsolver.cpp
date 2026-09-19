// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "fem/substructure.hpp"
#include "models/substructuringsolver.hpp"
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

}  // namespace

TEST_CASE("SubstructuringSolver reproduces full-domain electrostatics",
          "[substructure][Serial][Parallel]")
{
  auto run = [](double eps_r, double eps_e, int order)
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
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));

    // Region-condensed solve (parent true-DOF solution).
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    Vector u = ss.SolveRegion();

    // Reuse: a second solve reuses the materialized environment DtN (no re-condensation)
    // and must give an identical result.
    Vector u2 = ss.SolveRegion();
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
    // The default SolveRegion excitation drives the lowest terminal, matching phi[0].
    CHECK(std::abs(2.0 * e_sub - c00) <= 1.0e-9 * std::abs(c00));
  };

  SECTION("uniform permittivity, order 1")
  {
    run(1.0, 1.0, 1);
  }
  SECTION("contrast across interface, order 1")
  {
    run(1.0, 10.0, 1);
  }
  SECTION("contrast across interface, order 2")
  {
    run(1.0, 10.0, 2);
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

TEST_CASE("SubstructuringSolver reproduces full-domain magnetostatics",
          "[substructure][Serial][Parallel]")
{
  // H(curl) region-condensed solve of the regularized curl-curl operator vs a monolith
  // using the identical operator (1/mu curl-curl + subdomain unit mass regularization).
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
  Vector u = ss.SolveRegion();

  // Monolith reference on the same ND parent space with the identical operator.
  auto &pmesh = mesh.back()->Get();
  mfem::ND_FECollection fec(order, 3);
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  const int max_attr = pmesh.attributes.Max();
  mfem::Vector nu_by_attr(max_attr), mass_by_attr(max_attr);
  nu_by_attr = 0.0;
  mass_by_attr = 1.0;
  nu_by_attr(0) = 1.0 / mu_r;
  nu_by_attr(1) = 1.0 / mu_e;
  mfem::PWConstCoefficient nu(nu_by_attr), mass(mass_by_attr);
  mfem::ParBilinearForm a(&pfes);
  a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(nu));
  a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mass));
  a.Assemble();
  mfem::ParGridFunction xgf(&pfes);
  xgf = 0.0;
  const int maxb = pmesh.bdr_attributes.Max();
  mfem::Array<int> ess_bdr(maxb), ess_tdofs, t1_bdr(maxb), t1_tdofs;
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  pfes.GetEssentialTrueDofs(ess_bdr, ess_tdofs);
  // Match the substructuring Dirichlet data exactly: raw edge DOF = 1 on the driven
  // terminal (attr 1), 0 elsewhere.
  t1_bdr = 0;
  t1_bdr[0] = 1;
  pfes.GetEssentialTrueDofs(t1_bdr, t1_tdofs);
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
  a.FormLinearSystem(ess_tdofs, xgf, bform, A, X, B);
  mfem::HypreParMatrix *Ah = A.As<mfem::HypreParMatrix>();
  mfem::HypreAMS ams(*Ah, &pfes);
  ams.SetPrintLevel(0);
  mfem::HyprePCG pcg(*Ah);
  pcg.SetTol(1e-12);
  pcg.SetMaxIter(1000);
  pcg.SetPrintLevel(0);
  pcg.SetPreconditioner(ams);
  pcg.Mult(B, X);
  a.RecoverFEMSolution(X, bform, xgf);
  mfem::Vector u_full;
  xgf.GetTrueDofs(u_full);

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
  CHECK(std::sqrt(gnum / gden) < 1.0e-7);
}

TEST_CASE("SubstructuringSolver offline/online model reuse",
          "[substructure][Serial][Parallel]")
{
  const std::string model_path = "substruct_model_roundtrip.bin";
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

}  // namespace palace
