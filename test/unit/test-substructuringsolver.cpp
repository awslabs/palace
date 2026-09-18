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

// Build a unit-cube ParMesh split at x=0.5 (domain attrs 1/2) with boundary attrs 1 (x=0),
// 2 (x=1), 3 (all other faces).
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
    int a = 3;
    if (xc < 1e-9)
    {
      a = 1;
    }
    else if (xc > 1.0 - 1e-9)
    {
      a = 2;
    }
    serial.SetBdrAttribute(b, a);
  }
  serial.SetAttributes();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

}  // namespace

TEST_CASE("SubstructuringSolver reproduces full-domain electrostatics",
          "[substructure][Serial]")
{
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("SubstructuringSolver test is serial-only");
  }
  auto run = [](double eps_r, double eps_e)
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
         {{"Order", 1},
          {"Substructuring",
           {{"Region", {{"Attributes", {1}}}}, {"Environment", {{"Attributes", {2}}}}}}}}};
    IoData iodata(config, false);

    std::vector<std::unique_ptr<Mesh>> mesh;
    mesh.push_back(std::make_unique<Mesh>(MakeSplitCube(6)));

    // Region-condensed solve.
    SubstructuringSolver ss(iodata, mesh);
    ss.CondenseEnvironment();
    Vector u_region = ss.SolveRegion();

    // Full-domain reference: assemble grad(eps grad) with terminal Dirichlet and solve.
    auto &pmesh = mesh.back()->Get();
    mfem::H1_FECollection fec(1, 3);
    mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
    const int N = pfes.GetVSize();
    mfem::Vector eps_by_attr(pmesh.attributes.Max());
    eps_by_attr = 1.0;
    for (const auto &m : iodata.domains.materials)
    {
      for (int a : m.attributes)
      {
        eps_by_attr(a - 1) = m.epsilon_r.s[0];
      }
    }
    mfem::PWConstCoefficient eps(eps_by_attr);
    mfem::BilinearForm a(&pfes);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(eps));
    a.Assemble();
    a.Finalize();
    mfem::SparseMatrix &A = a.SpMat();
    // Dirichlet: attr 1 -> 1 V, attr 2 -> 0 V (by vertex x-coordinate, order-1 H1).
    std::vector<char> dir(N, 0);
    std::vector<double> udir(N, 0.0);
    for (int v = 0; v < pmesh.GetNV(); v++)
    {
      const double *x = pmesh.GetVertex(v);
      if (x[0] < 1e-9)
      {
        dir[v] = 1;
        udir[v] = 1.0;
      }
      else if (x[0] > 1.0 - 1e-9)
      {
        dir[v] = 1;
        udir[v] = 0.0;
      }
    }
    std::vector<int> fl(N, -1), freed;
    for (int p = 0; p < N; p++)
    {
      if (!dir[p])
      {
        fl[p] = (int)freed.size();
        freed.push_back(p);
      }
    }
    const int nf = freed.size();
    mfem::DenseMatrix Kf(nf);
    Kf = 0.0;
    mfem::Vector bf(nf);
    bf = 0.0;
    for (int af = 0; af < nf; af++)
    {
      int p = freed[af];
      const int *cols = A.GetRowColumns(p);
      const double *vals = A.GetRowEntries(p);
      for (int k = 0; k < A.RowSize(p); k++)
      {
        int q = cols[k];
        double v = vals[k];
        if (fl[q] >= 0)
        {
          Kf(af, fl[q]) += v;
        }
        else if (dir[q])
        {
          bf(af) -= v * udir[q];
        }
      }
    }
    mfem::DenseMatrix Kfi(Kf);
    Kfi.Invert();
    mfem::Vector uf(nf);
    Kfi.Mult(bf, uf);
    std::vector<double> u_full(N, 0.0);
    for (int p = 0; p < N; p++)
    {
      if (dir[p])
      {
        u_full[p] = udir[p];
      }
    }
    for (int af = 0; af < nf; af++)
    {
      u_full[freed[af]] = uf(af);
    }

    // Map region submesh DOFs to parent for comparison.
    mfem::Array<int> ra(1);
    ra[0] = 1;
    Substructure region(pfes, ra, fec);
    const auto &par = region.GetParentDof();
    double num = 0.0, den = 0.0;
    for (int i = 0; i < u_region.Size(); i++)
    {
      double e = u_region(i) - u_full[par[i]];
      num += e * e;
      den += u_full[par[i]] * u_full[par[i]];
    }
    CHECK(std::sqrt(num / den) < 1.0e-10);
  };

  SECTION("uniform permittivity")
  {
    run(1.0, 1.0);
  }
  SECTION("contrast across interface")
  {
    run(1.0, 10.0);
  }
}

}  // namespace palace
