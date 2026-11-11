// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// A unit cube of n³ tetrahedra or hexahedra whose vertex orders are randomly permuted
// (tetrahedra) or transformed by a random symmetry of the cube (hexahedra), so that the
// shared faces of neighboring elements take all relative orientations, including the
// rotations, which only meshes with mixed element orientations have.
std::unique_ptr<mfem::Mesh> MixedOrientationMesh(int n, bool hex, unsigned int seed)
{
  auto base = mfem::Mesh::MakeCartesian3D(
      n, n, n, hex ? mfem::Element::HEXAHEDRON : mfem::Element::TETRAHEDRON);
  auto mesh = std::make_unique<mfem::Mesh>(3, base.GetNV(), base.GetNE(), base.GetNBE());
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh->AddVertex(base.GetVertex(v));
  }
  std::mt19937 gen(seed);
  // Reference coordinates of the hexahedron vertices.
  constexpr int ref[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                             {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  for (int e = 0; e < base.GetNE(); e++)
  {
    const int *v = base.GetElement(e)->GetVertices();
    if (!hex)
    {
      int w[4];
      std::copy_n(v, 4, w);
      std::shuffle(w, w + 4, gen);
      mesh->AddTet(w, 1);
      continue;
    }
    int perm[3] = {0, 1, 2}, w[8];
    std::shuffle(perm, perm + 3, gen);
    const int flip = static_cast<int>(gen() % 8);
    for (int k = 0; k < 8; k++)
    {
      int c[3];
      for (int d = 0; d < 3; d++)
      {
        c[d] = (flip & (1 << d)) ? 1 - ref[k][perm[d]] : ref[k][perm[d]];
      }
      for (int j = 0; j < 8; j++)
      {
        if (ref[j][0] == c[0] && ref[j][1] == c[1] && ref[j][2] == c[2])
        {
          w[k] = v[j];
        }
      }
    }
    mesh->AddHex(w, 1);
  }
  for (int b = 0; b < base.GetNBE(); b++)
  {
    if (hex)
    {
      mesh->AddBdrQuad(base.GetBdrElement(b)->GetVertices(), 1);
    }
    else
    {
      mesh->AddBdrTriangle(base.GetBdrElement(b)->GetVertices(), 1);
    }
  }
  // Keep the vertex orders (no reorientation).
  if (hex)
  {
    mesh->FinalizeHexMesh(1, 0, false);
  }
  else
  {
    mesh->FinalizeTetMesh(1, 0, false);
  }
  return mesh;
}

// For u a polynomial in the finite element space (of degree deg = p for H1, p - 1 for ND on
// an affine mesh), returns max |P R u - u| / max |u|, which is zero up to rounding when the
// conforming interpolation P of the space is correct.
double ProlongationError(mfem::ParFiniteElementSpace &fespace, int deg)
{
  auto f = [deg](const mfem::Vector &x, int c)
  {
    double v = 0.0;
    int k = 0;
    for (int a = 0; a <= deg; a++)
    {
      for (int b = 0; a + b <= deg; b++)
      {
        for (int e = 0; a + b + e <= deg; e++, k++)
        {
          v += (std::sin(1.3 + 0.77 * k + 2.1 * c) + 0.1) * std::pow(x(0), a) *
               std::pow(x(1), b) * std::pow(x(2), e);
        }
      }
    }
    return v;
  };
  mfem::ParGridFunction u(&fespace);
  if (fespace.FEColl()->GetRangeType(3) == mfem::FiniteElement::SCALAR)
  {
    mfem::FunctionCoefficient coeff([&f](const mfem::Vector &x) { return f(x, 0); });
    u.ProjectCoefficient(coeff);
  }
  else
  {
    mfem::VectorFunctionCoefficient coeff(3,
                                          [&f](const mfem::Vector &x, mfem::Vector &v)
                                          {
                                            for (int c = 0; c < 3; c++)
                                            {
                                              v(c) = f(x, c + 1);
                                            }
                                          });
    u.ProjectCoefficient(coeff);
  }
  mfem::Vector t(fespace.GetTrueVSize()), y(fespace.GetVSize());
  fespace.GetRestrictionMatrix()->Mult(u, t);
  fespace.GetProlongationMatrix()->Mult(t, y);
  y -= u;
  double err[2] = {y.Normlinf(), u.Normlinf()};
  Mpi::GlobalMax(2, err, fespace.GetComm());
  return err[0] / err[1];
}

}  // namespace

TEST_CASE("Nonconforming interpolation on faces with rotated orientations",
          "[ncmesh][Serial][Parallel]")
{
  // The rows of the conforming interpolation which a process receives for the DOFs of a
  // shared face are reordered with the relative orientation of the face, whose inverse is
  // only different for rotations (MFEM had the inverse for nonconforming meshes). H1 has
  // reordered face DOFs from order 4 on triangles and 3 on quadrilaterals, Nédélec from
  // order 2.
  MPI_Comm comm = Mpi::World();
  for (bool hex : {false, true})
  {
    auto smesh = MixedOrientationMesh(3, hex, 1);
    smesh->EnsureNCMesh(true);
    {
      // Also have master and slave faces.
      mfem::Array<int> refs;
      refs.Append(0);
      smesh->GeneralRefinement(refs, 1);
    }

    // A random partition: most faces are shared between processes.
    std::vector<int> partitioning(smesh->GetNE());
    std::mt19937 gen(2);
    for (auto &p : partitioning)
    {
      p = static_cast<int>(gen() % Mpi::Size(comm));
    }
    for (int r = 0; r < Mpi::Size(comm); r++)
    {
      partitioning[r] = r;
    }
    mfem::ParMesh pmesh(comm, *smesh, partitioning.data());

    for (int p : {hex ? 3 : 4, 5})
    {
      mfem::H1_FECollection fec(p, 3);
      mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
      INFO((hex ? "hexahedra" : "tetrahedra") << ", H1 order " << p);
      CHECK(ProlongationError(fespace, p) < 1.0e-12);
    }
    for (int p : {2, 3})
    {
      mfem::ND_FECollection fec(p, 3);
      mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
      INFO((hex ? "hexahedra" : "tetrahedra") << ", ND order " << p);
      CHECK(ProlongationError(fespace, p - 1) < 1.0e-12);
    }
  }
}

}  // namespace palace
