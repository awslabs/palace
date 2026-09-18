// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/substructure.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Scatter a submesh element-space sparse matrix into the parent DOF space with the signed
// map: A_par(par[i], par[j]) += sign[i]*sign[j]*A(i,j).
void AddToParent(const mfem::SparseMatrix &A, const std::vector<int> &par,
                 const std::vector<double> &sgn, mfem::DenseMatrix &A_par)
{
  for (int i = 0; i < A.Height(); i++)
  {
    const int *cols = A.GetRowColumns(i);
    const double *vals = A.GetRowEntries(i);
    for (int k = 0; k < A.RowSize(i); k++)
    {
      A_par(par[i], par[cols[k]]) += sgn[i] * sgn[cols[k]] * vals[k];
    }
  }
}

void DensifyInto(const mfem::SparseMatrix &A, mfem::DenseMatrix &D)
{
  for (int i = 0; i < A.Height(); i++)
  {
    const int *cols = A.GetRowColumns(i);
    const double *vals = A.GetRowEntries(i);
    for (int k = 0; k < A.RowSize(i); k++)
    {
      D(i, cols[k]) += vals[k];
    }
  }
}

// Assemble region+environment independently on submeshes and check the signed composition
// into parent-DOF space reproduces a direct parent assembly to machine precision.
double CompositionError(bool nedelec, int nx)
{
  mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
  for (int e = 0; e < serial.GetNE(); e++)
  {
    mfem::Vector c;
    serial.GetElementCenter(e, c);
    serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
  }
  serial.SetAttributes();
  mfem::ParMesh mesh(Mpi::World(), serial);

  std::unique_ptr<mfem::FiniteElementCollection> fec;
  if (nedelec)
  {
    fec = std::make_unique<mfem::ND_FECollection>(1, 3);
  }
  else
  {
    fec = std::make_unique<mfem::H1_FECollection>(1, 3);
  }
  mfem::ParFiniteElementSpace pfes(&mesh, fec.get());
  const int N = pfes.GetVSize();

  auto assemble = [&](mfem::ParFiniteElementSpace &fes)
  {
    auto *a = new mfem::BilinearForm(&fes);
    if (nedelec)
    {
      a->AddDomainIntegrator(new mfem::CurlCurlIntegrator);
      a->AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    }
    else
    {
      a->AddDomainIntegrator(new mfem::DiffusionIntegrator);
      a->AddDomainIntegrator(new mfem::MassIntegrator);
    }
    a->Assemble();
    a->Finalize();
    return a;
  };

  // Direct parent assembly (reference).
  std::unique_ptr<mfem::BilinearForm> ap(assemble(pfes));
  mfem::DenseMatrix A_ref(N);
  A_ref = 0.0;
  DensifyInto(ap->SpMat(), A_ref);

  // Independent subdomain assembly + signed composition into parent space.
  mfem::DenseMatrix A_comp(N);
  A_comp = 0.0;
  bool map_ok = true;
  for (int attr : {1, 2})
  {
    mfem::Array<int> a(1);
    a[0] = attr;
    mfem::ParSubMesh sm = mfem::ParSubMesh::CreateFromDomain(mesh, a);
    mfem::ParFiniteElementSpace fes(&sm, fec.get());
    mfem::Array<int> emap = sm.GetParentElementIDMap();
    std::vector<int> par;
    std::vector<double> sgn;
    map_ok = BuildSubMeshDofMap(fes, pfes, emap, par, sgn) && map_ok;
    std::unique_ptr<mfem::BilinearForm> ak(assemble(fes));
    AddToParent(ak->SpMat(), par, sgn, A_comp);
  }
  CHECK(map_ok);

  A_comp -= A_ref;
  double num = A_comp.FNorm(), den = A_ref.FNorm();
  return num / den;
}

}  // namespace

TEST_CASE("Substructure signed composition reproduces the monolith",
          "[substructure][Serial]")
{
  // The composition check compares dense parent-DOF-space matrices, which is a serial
  // construction; parallel interface handling is a later milestone.
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("substructure composition test is serial-only");
  }
  SECTION("H1 (scalar)")
  {
    CHECK(CompositionError(false, 6) < 1.0e-12);
    CHECK(CompositionError(false, 8) < 1.0e-12);
  }
  SECTION("H(curl) (signed edge dofs)")
  {
    CHECK(CompositionError(true, 6) < 1.0e-12);
    CHECK(CompositionError(true, 8) < 1.0e-12);
  }
}

TEST_CASE("Substructure interface partition", "[substructure][Serial]")
{
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("substructure partition test is serial-only");
  }
  // Cube split at x = 0.5 into region (attr 1) and environment (attr 2). The interface is
  // the x = 0.5 plane; its DOF count is exact and geometry-determined.
  auto interface_count = [](bool nedelec, int nx)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    if (nedelec)
    {
      fec = std::make_unique<mfem::ND_FECollection>(1, 3);
    }
    else
    {
      fec = std::make_unique<mfem::H1_FECollection>(1, 3);
    }
    mfem::ParFiniteElementSpace pfes(&mesh, fec.get());
    mfem::Array<int> ar(1), ae(1);
    ar[0] = 1;
    ae[0] = 2;
    Substructure region(pfes, ar, *fec), env(pfes, ae, *fec);
    CHECK(region.ConformingMap());
    CHECK(env.ConformingMap());
    std::vector<int> owner;
    MarkParentDofOwnership({&region, &env}, pfes.GetVSize(), owner);
    int n_iface = 0, n_covered = 0;
    for (int o : owner)
    {
      if (o == 3)
      {
        n_iface++;
      }
      if (o != 0)
      {
        n_covered++;
      }
    }
    // Region and environment together cover every parent DOF.
    CHECK(n_covered == pfes.GetVSize());
    return n_iface;
  };

  // H1: 7x7 vertices on the x=0.5 plane. H(curl): in-plane edges = 2*6*7.
  CHECK(interface_count(false, 6) == 49);
  CHECK(interface_count(true, 6) == 84);
}

TEST_CASE("DtN condensation reproduces the monolith and is reusable",
          "[substructure][Serial]")
{
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("DtN test is serial-only");
  }
  // Region (attr 1) redesigned while the environment (attr 2) is condensed once to S_E,
  // g_E. Each region design must reproduce a fresh monolithic solve, reusing the same DtN
  // operator.
  auto run = [](bool nedelec)
  {
    const int nx = 6;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    if (nedelec)
    {
      fec = std::make_unique<mfem::ND_FECollection>(1, 3);
    }
    else
    {
      fec = std::make_unique<mfem::H1_FECollection>(1, 3);
    }
    mfem::ParFiniteElementSpace pfes(&mesh, fec.get());
    const int N = pfes.GetVSize();
    mfem::Array<int> ar(1), ae(1);
    ar[0] = 1;
    ae[0] = 2;
    Substructure region(pfes, ar, *fec), env(pfes, ae, *fec);
    std::vector<int> owner, gamma_index;
    MarkParentDofOwnership({&region, &env}, N, owner);
    const int nG = BuildInterfaceIndex(owner, gamma_index);

    // Assemble an operator (+ unit load) on a given substructure. cmass scales the mass
    // term (the region "design"); the environment always uses cmass = 1.
    auto assemble = [&](mfem::ParFiniteElementSpace &fes, double cmass,
                        std::unique_ptr<mfem::SparseMatrix> &A, mfem::Vector &f)
    {
      mfem::ConstantCoefficient cm(cmass);
      mfem::BilinearForm a(&fes);
      mfem::LinearForm l(&fes);
      if (nedelec)
      {
        a.AddDomainIntegrator(new mfem::CurlCurlIntegrator);
        a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(cm));
        mfem::Vector one3(3);
        one3 = 1.0;
        auto *vc = new mfem::VectorConstantCoefficient(one3);
        l.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*vc));
        l.Assemble();
        delete vc;
      }
      else
      {
        a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
        a.AddDomainIntegrator(new mfem::MassIntegrator(cm));
        mfem::ConstantCoefficient one(1.0);
        l.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        l.Assemble();
      }
      a.Assemble();
      a.Finalize();
      A = std::make_unique<mfem::SparseMatrix>(a.SpMat());
      f.SetSize(fes.GetVSize());
      f = l;
    };

    // Environment condensed once.
    std::unique_ptr<mfem::SparseMatrix> Ae;
    mfem::Vector fe;
    assemble(env.GetFESpace(), 1.0, Ae, fe);
    DtNBoundaryOperator dtn(env, gamma_index, *Ae, fe);

    // Scatter a substructure operator/load into a compact (region-interior + interface)
    // system. Region interior DOFs get compact indices after the nG interface DOFs.
    auto region_solve = [&](double cmass, mfem::Vector &u_parent)
    {
      std::unique_ptr<mfem::SparseMatrix> Ar;
      mfem::Vector fr;
      assemble(region.GetFESpace(), cmass, Ar, fr);
      const auto &par = region.GetParentDof();
      const auto &sgn = region.GetSign();
      const int nr = Ar->Height();
      // Region interior compact index (parent dof -> local), interface uses gamma_index.
      std::vector<int> ri(N, -1);
      int nI = 0;
      for (int i = 0; i < nr; i++)
      {
        if (gamma_index[par[i]] < 0 && ri[par[i]] < 0)
        {
          ri[par[i]] = nG + nI++;
        }
      }
      const int n = nG + nI;
      mfem::DenseMatrix K(n);
      K = 0.0;
      mfem::Vector b(n);
      b = 0.0;
      auto idx = [&](int i)
      {
        int p = par[i];
        return gamma_index[p] >= 0 ? gamma_index[p] : ri[p];
      };
      for (int i = 0; i < nr; i++)
      {
        const int *cols = Ar->GetRowColumns(i);
        const double *vals = Ar->GetRowEntries(i);
        for (int k = 0; k < Ar->RowSize(i); k++)
        {
          const int j = cols[k];
          K(idx(i), idx(j)) += sgn[i] * sgn[j] * vals[k];
        }
        b(idx(i)) += sgn[i] * fr(i);
      }
      // Add environment DtN on the interface block/load.
      for (int a = 0; a < nG; a++)
      {
        b(a) += dtn.Load()(a);
        for (int bb = 0; bb < nG; bb++)
        {
          K(a, bb) += dtn.Schur()(a, bb);
        }
      }
      mfem::DenseMatrix Ki(K);
      Ki.Invert();
      mfem::Vector u(n);
      Ki.Mult(b, u);
      // Scatter back to parent DOFs. The compact system is already in parent orientation,
      // so no sign is applied here.
      u_parent.SetSize(N);
      u_parent = 0.0;
      for (int i = 0; i < nr; i++)
      {
        u_parent(par[i]) = u(idx(i));
      }
    };

    auto monolith = [&](double cmass, mfem::Vector &u_parent)
    {
      std::unique_ptr<mfem::SparseMatrix> Ar, Ae2;
      mfem::Vector fr, fe2;
      assemble(region.GetFESpace(), cmass, Ar, fr);
      assemble(env.GetFESpace(), 1.0, Ae2, fe2);
      mfem::DenseMatrix A(N);
      A = 0.0;
      mfem::Vector b(N);
      b = 0.0;
      auto add = [&](Substructure &s, mfem::SparseMatrix &As, mfem::Vector &fs)
      {
        const auto &par = s.GetParentDof();
        const auto &sgn = s.GetSign();
        for (int i = 0; i < As.Height(); i++)
        {
          const int *cols = As.GetRowColumns(i);
          const double *vals = As.GetRowEntries(i);
          for (int k = 0; k < As.RowSize(i); k++)
          {
            A(par[i], par[cols[k]]) += sgn[i] * sgn[cols[k]] * vals[k];
          }
          b(par[i]) += sgn[i] * fs(i);
        }
      };
      add(region, *Ar, fr);
      add(env, *Ae2, fe2);
      mfem::DenseMatrix Ai(A);
      Ai.Invert();
      u_parent.SetSize(N);
      Ai.Mult(b, u_parent);
    };

    for (double cmass : {1.0, 5.0})
    {
      mfem::Vector us, um;
      region_solve(cmass, us);
      monolith(cmass, um);
      // Compare on region-owned parent DOFs (region interior + interface).
      double num = 0.0, den = 0.0;
      for (int p = 0; p < N; p++)
      {
        if (owner[p] & 1)
        {
          double e = us(p) - um(p);
          num += e * e;
          den += um(p) * um(p);
        }
      }
      CHECK(std::sqrt(num / den) < 1.0e-10);
    }
  };

  SECTION("H1 (scalar)")
  {
    run(false);
  }
  SECTION("H(curl) (definite, signed)")
  {
    run(true);
  }
}

TEST_CASE("Region DtN operator solves iteratively to the monolith",
          "[substructure][Serial]")
{
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("region DtN operator test is serial-only");
  }
  // Solve (K_region + S_E on Gamma) u = b_region + g_E with CG via RegionDtNOperator, and
  // check it reproduces the monolithic solve on region DOFs.
  auto run = [](bool nedelec)
  {
    const int nx = 6;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    std::unique_ptr<mfem::FiniteElementCollection> fec;
    if (nedelec)
    {
      fec = std::make_unique<mfem::ND_FECollection>(1, 3);
    }
    else
    {
      fec = std::make_unique<mfem::H1_FECollection>(1, 3);
    }
    mfem::ParFiniteElementSpace pfes(&mesh, fec.get());
    const int N = pfes.GetVSize();
    mfem::Array<int> ar(1), ae(1);
    ar[0] = 1;
    ae[0] = 2;
    Substructure region(pfes, ar, *fec), env(pfes, ae, *fec);
    std::vector<int> owner, gamma_index;
    MarkParentDofOwnership({&region, &env}, N, owner);
    BuildInterfaceIndex(owner, gamma_index);

    auto assemble = [&](mfem::ParFiniteElementSpace &fes,
                        std::unique_ptr<mfem::SparseMatrix> &A, mfem::Vector &f)
    {
      mfem::ConstantCoefficient one(1.0);
      mfem::BilinearForm a(&fes);
      mfem::LinearForm l(&fes);
      if (nedelec)
      {
        a.AddDomainIntegrator(new mfem::CurlCurlIntegrator);
        a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
        mfem::Vector one3(3);
        one3 = 1.0;
        auto *vc = new mfem::VectorConstantCoefficient(one3);
        l.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*vc));
        l.Assemble();
        delete vc;
      }
      else
      {
        a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
        a.AddDomainIntegrator(new mfem::MassIntegrator);
        l.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        l.Assemble();
      }
      a.Assemble();
      a.Finalize();
      A = std::make_unique<mfem::SparseMatrix>(a.SpMat());
      f.SetSize(fes.GetVSize());
      f = l;
    };

    std::unique_ptr<mfem::SparseMatrix> Ae;
    mfem::Vector fe;
    assemble(env.GetFESpace(), Ae, fe);
    DtNBoundaryOperator dtn(env, gamma_index, *Ae, fe);

    // Region operator + RHS on the region submesh (local orientation).
    std::unique_ptr<mfem::SparseMatrix> Ar;
    mfem::Vector fr;
    assemble(region.GetFESpace(), Ar, fr);
    const auto &par = region.GetParentDof();
    const auto &sgn = region.GetSign();
    const int nr = Ar->Height();
    mfem::Vector b(fr);
    for (int i = 0; i < nr; i++)
    {
      const int g = gamma_index[par[i]];
      if (g >= 0)
      {
        b(i) += sgn[i] * dtn.Load()(g);  // parent -> local orientation
      }
    }

    // Combined operator K_region + S_E-on-Gamma, solved with CG.
    RegionDtNOperator dtn_op(region, gamma_index, dtn);
    mfem::SumOperator combined(Ar.get(), 1.0, &dtn_op, 1.0, false, false);
    mfem::CGSolver cg;
    cg.SetOperator(combined);
    cg.SetRelTol(1.0e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(0);
    mfem::Vector u(nr);
    u = 0.0;
    cg.Mult(b, u);
    CHECK(cg.GetConverged());

    // Monolithic reference.
    mfem::DenseMatrix A(N);
    A = 0.0;
    mfem::Vector rhs(N);
    rhs = 0.0;
    auto add = [&](Substructure &s, mfem::SparseMatrix &As, mfem::Vector &fs)
    {
      const auto &p = s.GetParentDof();
      const auto &sg = s.GetSign();
      for (int i = 0; i < As.Height(); i++)
      {
        const int *cols = As.GetRowColumns(i);
        const double *vals = As.GetRowEntries(i);
        for (int k = 0; k < As.RowSize(i); k++)
        {
          A(p[i], p[cols[k]]) += sg[i] * sg[cols[k]] * vals[k];
        }
        rhs(p[i]) += sg[i] * fs(i);
      }
    };
    add(region, *Ar, fr);
    add(env, *Ae, fe);
    mfem::DenseMatrix Ai(A);
    Ai.Invert();
    mfem::Vector um(N);
    Ai.Mult(rhs, um);

    // Compare region solution (local orientation) to the monolith (parent orientation).
    double num = 0.0, den = 0.0;
    for (int i = 0; i < nr; i++)
    {
      const double us = sgn[i] * u(i);  // local -> parent orientation
      const double e = us - um(par[i]);
      num += e * e;
      den += um(par[i]) * um(par[i]);
    }
    CHECK(std::sqrt(num / den) < 1.0e-8);
  };

  SECTION("H1 (scalar)")
  {
    run(false);
  }
  SECTION("H(curl) (definite, signed)")
  {
    run(true);
  }
}

TEST_CASE("DtN with Dirichlet terminals reproduces electrostatics",
          "[substructure][Serial]")
{
  if (Mpi::Size(Mpi::World()) > 1)
  {
    SKIP("electrostatic DtN test is serial-only");
  }
  // Pure-Laplace two-terminal capacitor: phi=1 on x=0 (region), phi=0 on x=1 (environment).
  // The environment is condensed with its ground terminal as a Dirichlet DOF (keeping A_EE
  // nonsingular); the region is solved with its drive terminal + the DtN, and must
  // reproduce the full-domain solve. Includes a permittivity contrast across the interface.
  auto run = [](double eps_r, double eps_e)
  {
    const int nx = 6;
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(nx, nx, nx, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::H1_FECollection fec(1, 3);
    mfem::ParFiniteElementSpace pfes(&mesh, &fec);
    const int N = pfes.GetVSize();

    // Dirichlet terminal DOFs by vertex x-coordinate (order-1 H1: DOF == vertex).
    std::vector<char> dir(N, 0);
    std::vector<double> udir(N, 0.0);
    for (int v = 0; v < mesh.GetNV(); v++)
    {
      const double *x = mesh.GetVertex(v);
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

    mfem::Array<int> ar(1), ae(1);
    ar[0] = 1;
    ae[0] = 2;
    Substructure region(pfes, ar, fec), env(pfes, ae, fec);
    std::vector<int> owner, gamma_index;
    MarkParentDofOwnership({&region, &env}, N, owner);
    BuildInterfaceIndex(owner, gamma_index);

    auto assemble = [&](mfem::ParFiniteElementSpace &fes, double eps,
                        std::unique_ptr<mfem::SparseMatrix> &A)
    {
      mfem::ConstantCoefficient ec(eps);
      mfem::BilinearForm a(&fes);
      a.AddDomainIntegrator(new mfem::DiffusionIntegrator(ec));
      a.Assemble();
      a.Finalize();
      A = std::make_unique<mfem::SparseMatrix>(a.SpMat());
    };

    // Environment submesh Dirichlet marker/values (map parent Dirichlet -> env submesh
    // DOF).
    auto submesh_dbc = [&](Substructure &s, std::vector<char> &m, mfem::Vector &vals)
    {
      const auto &par = s.GetParentDof();
      const int ns = par.size();
      m.assign(ns, 0);
      vals.SetSize(ns);
      vals = 0.0;
      for (int i = 0; i < ns; i++)
      {
        if (dir[par[i]])
        {
          m[i] = 1;
          vals(i) = udir[par[i]];
        }
      }
    };

    std::unique_ptr<mfem::SparseMatrix> Ae;
    assemble(env.GetFESpace(), eps_e, Ae);
    std::vector<char> em;
    mfem::Vector ev;
    submesh_dbc(env, em, ev);
    mfem::Vector fe(Ae->Height());
    fe = 0.0;
    DtNBoundaryOperator dtn(env, gamma_index, *Ae, fe, em, ev);

    // Region system on (region-interior-free + interface), region Dirichlet eliminated.
    std::unique_ptr<mfem::SparseMatrix> Ar;
    assemble(region.GetFESpace(), eps_r, Ar);
    const auto &par = region.GetParentDof();
    const int nr = Ar->Height();
    std::vector<char> rm;
    mfem::Vector rv;
    submesh_dbc(region, rm, rv);
    // Compact region free DOFs: interface (gamma) + region interior (not Dirichlet, not
    // iface).
    std::vector<int> rint;
    std::vector<int> ri(N, -1);
    const int nG = dtn.Size();
    for (int i = 0; i < nr; i++)
    {
      if (!rm[i] && gamma_index[par[i]] < 0 && ri[par[i]] < 0)
      {
        ri[par[i]] = (int)rint.size();
        rint.push_back(par[i]);
      }
    }
    const int nRi = rint.size(), ndof = nG + nRi;
    auto idx = [&](int i)
    { return gamma_index[par[i]] >= 0 ? gamma_index[par[i]] : nG + ri[par[i]]; };
    mfem::DenseMatrix K(ndof);
    K = 0.0;
    mfem::Vector b(ndof);
    b = 0.0;
    for (int i = 0; i < nr; i++)
    {
      if (rm[i])
      {
        continue;
      }
      const int *cols = Ar->GetRowColumns(i);
      const double *vals = Ar->GetRowEntries(i);
      for (int k = 0; k < Ar->RowSize(i); k++)
      {
        const int j = cols[k];
        const double v = vals[k];
        if (rm[j])
        {
          b(idx(i)) -= v * rv(j);
          continue;
        }  // eliminate region Dirichlet
        K(idx(i), idx(j)) += v;
      }
    }
    for (int a = 0; a < nG; a++)
    {
      b(a) += dtn.Load()(a);
      for (int bb = 0; bb < nG; bb++)
      {
        K(a, bb) += dtn.Schur()(a, bb);
      }
    }
    mfem::DenseMatrix Ki(K);
    Ki.Invert();
    mfem::Vector u(ndof);
    Ki.Mult(b, u);

    // Full-domain reference.
    std::unique_ptr<mfem::SparseMatrix> ArF, AeF;
    assemble(region.GetFESpace(), eps_r, ArF);
    assemble(env.GetFESpace(), eps_e, AeF);
    mfem::DenseMatrix A(N);
    A = 0.0;
    auto add = [&](Substructure &s, mfem::SparseMatrix &As)
    {
      const auto &p = s.GetParentDof();
      for (int i = 0; i < As.Height(); i++)
      {
        const int *cols = As.GetRowColumns(i);
        const double *vals = As.GetRowEntries(i);
        for (int k = 0; k < As.RowSize(i); k++)
        {
          A(p[i], p[cols[k]]) += vals[k];
        }
      }
    };
    add(region, *ArF);
    add(env, *AeF);
    std::vector<int> freed;
    std::vector<int> fl(N, -1);
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
    for (int a = 0; a < nf; a++)
    {
      int p = freed[a];
      for (int q = 0; q < N; q++)
      {
        double v = A(p, q);
        if (v == 0.0)
        {
          continue;
        }
        if (fl[q] >= 0)
        {
          Kf(a, fl[q]) += v;
        }
        else if (dir[q])
        {
          bf(a) -= v * udir[q];
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
    for (int a = 0; a < nf; a++)
    {
      u_full[freed[a]] = uf(a);
    }

    // Compare region-condensed to full-domain on region free DOFs (interior + interface).
    double num = 0.0, den = 0.0;
    for (int i = 0; i < nr; i++)
    {
      if (rm[i])
      {
        continue;
      }
      double us = u(idx(i));
      double um = u_full[par[i]];
      double e = us - um;
      num += e * e;
      den += um * um;
    }
    CHECK(std::sqrt(num / den) < 1.0e-10);

    // Recover the environment field from the interface solution and check the full
    // reconstructed field (region + environment) matches the monolith everywhere.
    mfem::Vector u_gamma(nG);
    for (int a = 0; a < nG; a++)
    {
      u_gamma(a) = u(a);  // region compact solve is parent-oriented (H1)
    }
    mfem::Vector u_env = dtn.RecoverEnvironment(u_gamma);
    const auto &epar = env.GetParentDof();
    double enum_ = 0.0, eden = 0.0;
    for (int i = 0; i < u_env.Size(); i++)
    {
      double e = u_env(i) - u_full[epar[i]];
      enum_ += e * e;
      eden += u_full[epar[i]] * u_full[epar[i]];
    }
    CHECK(std::sqrt(enum_ / eden) < 1.0e-10);
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
