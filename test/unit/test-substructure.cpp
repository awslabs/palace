// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "fem/substructure.hpp"
#include "utils/communication.hpp"

namespace palace
{

TEST_CASE("Interface true DOFs are identified consistently in parallel",
          "[substructure][Serial][Parallel]")
{
  // The global interface true-DOF count must be independent of the MPI partition. Cube
  // split at x=0.5, H1: 49 interface DOFs at order 1 (7x7 plane), 169 at order 2.
  auto count = [](int order)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::H1_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&mesh, &fec);
    mfem::Array<int> ra(1), ea(1);
    ra[0] = 1;
    ea[0] = 2;
    mfem::Array<int> rm, em, im;
    return MarkInterfaceTrueDofs(pfes, ra, ea, rm, em, im);
  };
  CHECK(count(1) == 49);
  CHECK(count(2) == 169);
}

TEST_CASE("Saved Nédélec DOFs are mapped across meshes and partitions",
          "[substructure][Serial][Parallel]")
{
  // The same tetrahedra in a different element and vertex order (so different face
  // orientations, as a different partition gives): the true-DOF bases of second-order
  // Nédélec elements differ beyond signs on the faces. The signature map must transform a
  // linear form and the mass matrix of one space into those of the other (dual quantities,
  // x_cur = M^-T x_saved).
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  mfem::Mesh base = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::Mesh perm(3, base.GetNV(), base.GetNE());
  const int nv = base.GetNV();
  for (int v = 0; v < nv; v++)
  {
    perm.AddVertex(base.GetVertex((v * 7) % nv));  // 7 and nv = 27 are coprime
  }
  std::vector<int> inv(nv);
  for (int v = 0; v < nv; v++)
  {
    inv[(v * 7) % nv] = v;
  }
  for (int e = base.GetNE() - 1; e >= 0; e--)
  {
    const int *ev = base.GetElement(e)->GetVertices();
    // An even permutation of the vertices keeps the orientation.
    perm.AddTet(inv[ev[1]], inv[ev[2]], inv[ev[0]], inv[ev[3]]);
  }
  perm.FinalizeTetMesh(1, 0, true);
  auto measure = [order](mfem::Mesh &serial, std::vector<double> &sig,
                         std::vector<double> &b, std::vector<double> &M)
  {
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::ND_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace fes(&mesh, &fec);
    const int nt = fes.GetTrueVSize();
    const int n = static_cast<int>(fes.GlobalTrueVSize());
    std::vector<int> index(nt);
    for (int i = 0; i < nt; i++)
    {
      index[i] = static_cast<int>(fes.GetMyTDofOffset()) + i;
    }
    sig = TrueDofSignatures(fes, index, n);
    // A linear form and the mass matrix, replicated (n is small).
    mfem::VectorFunctionCoefficient f(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v(0) = std::sin(3.0 * x(1)) + x(2);
                                        v(1) = x(0) * x(2) * x(2);
                                        v(2) = std::cos(2.0 * x(0) + x(1));
                                      });
    mfem::ParLinearForm lf(&fes);
    lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f));
    lf.Assemble();
    std::unique_ptr<mfem::HypreParVector> bt(lf.ParallelAssemble());
    mfem::ParBilinearForm a(&fes);
    a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    b.assign(n, 0.0);
    M.assign(static_cast<std::size_t>(n) * n, 0.0);
    mfem::Vector e(nt), Ae(nt);
    for (int j = 0; j < n; j++)
    {
      e = 0.0;
      const int jl = j - static_cast<int>(fes.GetMyTDofOffset());
      if (jl >= 0 && jl < nt)
      {
        e(jl) = 1.0;
      }
      A->Mult(e, Ae);
      for (int i = 0; i < nt; i++)
      {
        M[static_cast<std::size_t>(j) * n + index[i]] = Ae(i);
      }
    }
    for (int i = 0; i < nt; i++)
    {
      b[index[i]] = (*bt)(i);
    }
    MPI_Allreduce(MPI_IN_PLACE, b.data(), n, MPI_DOUBLE, MPI_SUM, Mpi::World());
    MPI_Allreduce(MPI_IN_PLACE, M.data(), n * n, MPI_DOUBLE, MPI_SUM, Mpi::World());
  };
  std::vector<double> sig_s, b_s, M_s, sig_c, b_c, M_c;
  measure(base, sig_s, b_s, M_s);
  measure(perm, sig_c, b_c, M_c);
  const SignatureMap map = MatchSignatureBasis(sig_c, sig_s, 12);
  int blocks = 0;
  for (const auto &r : map.rows)
  {
    blocks += (r.size() > 1);
  }
  if (order > 1)
  {
    CHECK(blocks > 0);  // face DOFs whose basis is not a signed permutation of the other's
  }
  const auto b = map.Dual(b_s.data());
  const auto M = map.DualMatrix(M_s);
  double db = 0.0, mb = 0.0, dM = 0.0, mM = 0.0;
  for (std::size_t i = 0; i < b.size(); i++)
  {
    db = std::max(db, std::abs(b[i] - b_c[i]));
    mb = std::max(mb, std::abs(b_c[i]));
  }
  for (std::size_t q = 0; q < M.size(); q++)
  {
    dM = std::max(dM, std::abs(M[q] - M_c[q]));
    mM = std::max(mM, std::abs(M_c[q]));
  }
  CAPTURE(blocks, db, mb, dM, mM);
  CHECK(db <= 1.0e-12 * mb);
  CHECK(dM <= 1.0e-12 * mM);
}

}  // namespace palace
