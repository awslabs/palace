// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

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

}  // namespace palace
