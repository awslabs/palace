// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
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

TEST_CASE("Element moments tell sets of elements apart", "[substructure][Serial][Parallel]")
{
  // A box [0, 2] x [-0.3, 0.3]^2, and the same rotated by 90 degrees about the z axis: the
  // same count, volume and distances from the origin, told apart by the centroid moments
  // (of which some vanish, up to rounding). The moments are those of any partition, and a
  // subset differs.
  auto box = [](bool rotate)
  {
    mfem::Mesh serial =
        mfem::Mesh::MakeCartesian3D(4, 2, 2, mfem::Element::HEXAHEDRON, 2.0, 0.6, 0.6);
    serial.Transform(
        [rotate](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(1) -= 0.3;
          y(2) -= 0.3;
          if (rotate)
          {
            y(0) = -(x(1) - 0.3);
            y(1) = x(0);
          }
        });
    return serial;
  };
  mfem::Mesh serial = box(false), rotated = box(true);
  mfem::ParMesh mesh(Mpi::World(), serial), mesh_rotated(Mpi::World(), rotated);
  std::vector<int> part(serial.GetNE());
  for (int e = 0; e < serial.GetNE(); e++)
  {
    part[e] = (serial.GetNE() - 1 - e) % Mpi::Size(Mpi::World());  // another partition
  }
  mfem::ParMesh mesh_repartitioned(Mpi::World(), serial, part.data());
  const auto m = ElementMoments(mesh, {1});
  CHECK(m[0] == 16.0);
  CHECK(std::abs(m[1] - 0.72) <= 1.0e-12);
  CHECK(MomentsMatch(m, ElementMoments(mesh_repartitioned, {1}), 1.0e-12));
  const auto r = ElementMoments(mesh_rotated, {1});
  CHECK(r[0] == m[0]);
  CHECK(std::abs(r[1] - m[1]) <= 1.0e-12 * m[1]);
  CHECK(std::abs(r[2] - m[2]) <= 1.0e-12 * m[2]);
  CHECK_FALSE(MomentsMatch(m, r, 1.0e-9));
  std::vector<char> half(mesh.GetNE(), 0);
  for (int e = 0; e < mesh.GetNE(); e += 2)
  {
    half[e] = 1;
  }
  CHECK_FALSE(MomentsMatch(m, ElementMoments(mesh, {1}, half), 1.0e-9));
}

TEST_CASE("Nédélec DOFs at the points of saved ones", "[substructure][Serial]")
{
  // Point-tangent DOFs u(p).t, with signatures (t_b, p_a t_b): a DOF at a saved point
  // matches whatever its tangent (another one, as a face DOF in another basis, or the
  // opposite one), a DOF off the saved points does not.
  auto add = [](std::vector<double> &sig, std::array<double, 3> p, std::array<double, 3> t)
  {
    sig.insert(sig.end(), t.begin(), t.end());
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        sig.push_back(p[a] * t[b]);
      }
    }
  };
  std::vector<double> saved, candidates;
  add(saved, {0.1, 0.2, 0.3}, {1.0, 0.0, 0.0});
  add(saved, {0.5, -0.2, 0.0}, {0.0, 0.6, 0.8});
  add(saved, {0.5, 0.4, 0.0}, {0.0, 0.0, 1.0});
  add(candidates, {0.1, 0.2, 0.3}, {0.0, 1.0, 0.0});
  add(candidates, {0.5, -0.2, 0.0}, {0.0, -0.6, -0.8});
  add(candidates, {0.5, -0.2, 1.0e-6}, {0.0, 0.6, 0.8});
  add(candidates, {0.1 + 1.0e-12, 0.2, 0.3}, {1.0, 0.0, 0.0});
  add(candidates, {0.5, 0.4, 0.0}, {0.6, 0.8, 0.0});
  const auto match = MatchSignaturePoints(candidates, 12, saved, 12);
  REQUIRE(match.size() == 5);
  CHECK(match[0]);
  CHECK(match[1]);
  CHECK_FALSE(match[2]);
  CHECK(match[3]);
  CHECK(match[4]);
}

TEST_CASE("Versions of split Nédélec DOFs are matched by their sides",
          "[substructure][Serial]")
{
  // The versions of a DOF split by a two-sided film share its 12 moments and have side
  // moments of opposite signs; a DOF away from the film has none. An orientation flip
  // negates all 13 columns. Each current version must map to the saved version of its side
  // (of the same 12 moments, which alone would pair the first one with the wrong side).
  auto add = [](std::vector<double> &sig, std::array<double, 3> p, std::array<double, 3> t,
                double side)
  {
    sig.insert(sig.end(), t.begin(), t.end());
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        sig.push_back(p[a] * t[b]);
      }
    }
    sig.push_back(side);
  };
  const std::array<double, 3> p = {0.0, 0.3, 0.0}, t = {0.0, 0.1, 0.0}, q = {0.5, 0.2, 0.1},
                              u = {0.1, 0.0, 0.05};
  const double side = 2.0e-3;
  std::vector<double> saved, current;
  add(saved, p, t, -side);
  add(saved, p, t, side);
  add(saved, q, u, 0.0);
  add(current, p, t, -side);
  add(current, q, {-u[0], -u[1], -u[2]}, 0.0);
  add(current, p, {-t[0], -t[1], -t[2]}, -side);
  const SignatureMap map = MatchSignatureBasis(current, saved, 13);
  REQUIRE(map.rows.size() == 3);
  using Entry = std::vector<std::pair<int, double>>;
  CHECK(map.rows[0] == Entry{{0, 1.0}});
  CHECK(map.rows[1] == Entry{{2, -1.0}});
  CHECK(map.rows[2] == Entry{{1, -1.0}});

  // A version whose side moment matches neither saved one is rejected, even though the
  // combinations of the saved versions reproduce it.
  current[12] *= 2.0;
  CHECK_THROWS(MatchSignatureBasis(current, saved, 13));
}

TEST_CASE("H1 DOFs at the points of saved ones", "[substructure][Serial]")
{
  // H1 signatures are the DOF coordinates: a DOF matches a saved one at its point, within
  // the matching tolerance, and not off it.
  const std::vector<double> saved = {0.1, 0.2, 0.3, 0.5, -0.2, 0.0},
                            candidates = {0.5, -0.2, 0.0,  0.1,   0.2 + 1.0e-12,
                                          0.3, 0.5,  -0.2, 1.0e-6};
  const auto match = MatchSignaturePoints(candidates, 3, saved, 3);
  REQUIRE(match.size() == 3);
  CHECK(match[0]);
  CHECK(match[1]);
  CHECK_FALSE(match[2]);
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
  const auto b = map.DualRows(b_s, 1);
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

TEST_CASE("Eliminated rows without a stored diagonal become identity rows",
          "[substructure][Serial][Parallel]")
{
  // Per rank: rows 0 and 3 store no entry (as a device-assembled matrix has for a DOF
  // without element contributions), rows 1 and 2 are coupled, rows 4 and 5 diagonal.
  // Eliminating rows 0, 1 and 3 leaves the identity on them and decouples row 2.
  MPI_Comm comm = Mpi::World();
  const int n = 6, rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  const HYPRE_BigInt off = static_cast<HYPRE_BigInt>(rank) * n,
                     glob = static_cast<HYPRE_BigInt>(size) * n;
  std::vector<HYPRE_BigInt> part;
  if (HYPRE_AssumedPartitionCheck())
  {
    part = {off, off + n};
  }
  else
  {
    for (int r = 0; r <= size; r++)
    {
      part.push_back(static_cast<HYPRE_BigInt>(r) * n);
    }
  }
  std::vector<int> I = {0, 0, 2, 4, 4, 5, 6};
  std::vector<HYPRE_BigInt> J = {off + 1, off + 2, off + 1, off + 2, off + 4, off + 5};
  std::vector<double> V = {2.0, -1.0, -1.0, 2.0, 2.0, 2.0};
  mfem::HypreParMatrix A(comm, n, glob, glob, I.data(), J.data(), V.data(), part.data(),
                         part.data());
  mfem::Array<int> dofs({0, 1, 3});
  std::unique_ptr<mfem::HypreParMatrix> B = EliminateToIdentity(A, dofs);
  mfem::Vector x(n), y(n);
  x = 1.0;
  B->Mult(x, y);
  const double *yy = y.HostRead();
  const std::array<double, 6> expected = {1.0, 1.0, 2.0, 1.0, 2.0, 2.0};
  for (int i = 0; i < n; i++)
  {
    CAPTURE(i);
    CHECK(yy[i] == expected[i]);
  }
}

}  // namespace palace
