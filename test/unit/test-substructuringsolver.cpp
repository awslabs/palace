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

}  // namespace palace
