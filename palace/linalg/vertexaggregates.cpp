// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vertexaggregates.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/blockutils.hpp"
#include "linalg/vector.hpp"

namespace palace
{

namespace
{

using namespace blocks;

// An element is flat when two consecutive nonzero eigenvalues of its scaled matrix have a
// ratio below the gap, and two of its vertices are tied when their entry of the projector
// onto the eigenvectors above the gap is at least the tie threshold in magnitude.
//
// The gap. The smoother is a Chebyshev polynomial of the 4th kind of order k in the
// diagonally scaled operator: an error component with the eigenvalue λ = x λ_max is
// multiplied by sin((2k+1)θ) / ((2k+1) sin θ), sin²θ = x, which is at most
// 1 / ((2k+1) √x). This bound promises a reduction only for x > 1 / (2k+1)². In a region
// of flat elements a function which is constant across the thickness and oscillates along
// it has only weak energy, is not represented on a coarser level, and its x is the ratio
// of a weak to a strong eigenvalue of the scaled element matrix (both are scaled by the
// same assembled diagonal). The point smoother thus leaves it undamped when that ratio is
// below 1 / (2k+1)²: 1/81 for k = 4, the lowest default order. The value is used for every
// order: the argument is for the nodal Laplacian, and for an edge space the blocks are
// still needed at larger ratios than the argument gives for higher orders (at order 3, an
// aspect ratio of 10 needs blocks for the curl-curl operator and not for the Laplacian).
//
// The tie threshold has no such derivation: it separates the entries of the projector
// between two vertices on opposite sides of a flat element from those between two
// vertices on the same side, which is a property of the shape of the strong eigenvector.
// For tetrahedra it lies between two values: a flat tetrahedron with one vertex above
// the center of the opposite face has the entries 1/12 between two vertices of that face
// (not to be tied, or all vertices of a layer are) and 1/4 between the vertex and each of
// them (to be tied), as has a flat tetrahedron with two vertices on each side for its
// pairs across. It is a constant which is not derived: 0.09, 0.1 and 1/6 give the same
// iterations in the cases tried, and 0.24, at the upper end, loses ties which a mesh with
// hanging nodes needs.
constexpr double STRONG_GAP = 1.0 / 81.0;
constexpr double STRONG_TIE = 1.0e-1;

// Marks the identifier of a hanging vertex, which is a hash and not a true dof number.
constexpr Key SLAVE_BIT = Key(1) << 62;

std::shared_ptr<const VertexAggregates>
ComputeVertexAggregates(const FiniteElementSpace &h1_fespace, const ceed::Operator &op)
{
  const mfem::ParFiniteElementSpace &fes = h1_fespace.Get();
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  MPI_Comm comm = fes.GetComm();
  const int nv = mesh.GetNV();
  MFEM_VERIFY(fes.GetVSize() == nv && fes.GetMaxElementOrder() == 1,
              "Vertex aggregates require the lowest order nodal space of the mesh!");
  const int dim = mesh.Dimension();

  // Element matrices of the volume sub-operators (L-dof = vertex for an order-1 H1 space),
  // by element-local probing.
  std::vector<int> eptr(1, 0), edofs;
  std::vector<double> emat;
  ForEachSubOperator(op,
                     [&](CeedOperator sub)
                     {
                       Ceed ceed;
                       PalaceCeedCallBackend(CeedOperatorGetCeed(sub, &ceed));
                       CeedBasis basis;
                       CeedInt bdim;
                       PalaceCeedCall(ceed, CeedOperatorGetActiveBasis(sub, &basis));
                       PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &bdim));
                       PalaceCeedCall(ceed, CeedBasisDestroy(&basis));
                       if (bdim != dim)
                       {
                         return;
                       }
                       std::size_t m0 = 0;
                       ProbeSubOperator(
                           sub,
                           [&](CeedInt ne, CeedInt n, const CeedInt *offsets, CeedInt a,
                               const CeedScalar *y)
                           {
                             if (a == 0)
                             {
                               m0 = emat.size();
                               emat.resize(m0 + static_cast<std::size_t>(ne) * n * n);
                               for (CeedInt e = 0; e < ne; e++)
                               {
                                 for (CeedInt b = 0; b < n; b++)
                                 {
                                   edofs.push_back(
                                       offsets[static_cast<std::size_t>(e) * n + b]);
                                 }
                                 eptr.push_back(static_cast<int>(edofs.size()));
                               }
                             }
                             for (CeedInt e = 0; e < ne; e++)
                             {
                               for (CeedInt b = 0; b < n; b++)
                               {
                                 emat[m0 + (static_cast<std::size_t>(e) * n + b) * n + a] =
                                     y[static_cast<std::size_t>(e) * n + b];
                               }
                             }
                           });
                     });
  const int ne = static_cast<int>(eptr.size()) - 1;

  // Assembled diagonal at the L-level: exact (summed over the processes) for the true
  // vertices, the local part for the hanging ones.
  PRows P(h1_fespace.Dof_TrueDof_Matrix());
  std::vector<double> dl(nv, 0.0);
  {
    std::size_t m = 0;
    for (int e = 0; e < ne; e++)
    {
      const int n = eptr[e + 1] - eptr[e];
      for (int a = 0; a < n; a++)
      {
        dl[edofs[eptr[e] + a]] += emat[m + a * n + a];
      }
      m += static_cast<std::size_t>(n) * n;
    }
  }
  std::vector<char> slave(nv, 0);
  std::vector<double> d(nv);
  {
    Vector dls(nv), dt(fes.GetTrueVSize()), dv(nv);
    dls.UseDevice(true);
    dt.UseDevice(true);
    dv.UseDevice(true);
    auto *DLS = dls.HostWrite();
    for (int i = 0; i < nv; i++)
    {
      slave[i] = !P.IsUnit(i);
      DLS[i] = slave[i] ? 0.0 : dl[i];
    }
    h1_fespace.Dof_TrueDof_Matrix()->MultTranspose(dls, dt);
    h1_fespace.Dof_TrueDof_Matrix()->Mult(dt, dv);
    const auto *DV = dv.HostRead();
    for (int i = 0; i < nv; i++)
    {
      d[i] = slave[i] ? dl[i] : DV[i];
    }
  }

  // Global identity of the vertices: the true dof, or for a hanging vertex a key of its
  // masters and weights.
  std::vector<Key> gid(nv);
  {
    std::vector<std::pair<HYPRE_BigInt, double>> row;
    for (int i = 0; i < nv; i++)
    {
      row.clear();
      P.ForRow(i, [&](int c, bool offd, double w)
               { row.emplace_back(P.GlobalCol(c, offd), w); });
      MFEM_VERIFY(!row.empty(), "Vertex without a prolongation row!");
      if (!slave[i])
      {
        gid[i] = static_cast<Key>(row[0].first);
      }
      else
      {
        std::sort(row.begin(), row.end());
        Key h = 0x13198a2e03707344ULL;
        for (auto &[c, w] : row)
        {
          h = Combine(h, static_cast<Key>(c));
          h = Combine(h, static_cast<Key>(std::llround(w * 1048576.0)));
        }
        gid[i] = h | SLAVE_BIT;
      }
    }
  }

  // Strong ties from the element matrices.
  std::vector<Key> pairs;
  std::vector<int> tie_a, tie_b;
  std::vector<double> tie_w;
  {
    using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 8, 8>;
    Eigen::SelfAdjointEigenSolver<Mat> es;
    std::size_t m = 0;
    for (int e = 0; e < ne; e++)
    {
      const int n = eptr[e + 1] - eptr[e];
      const int *dofs = &edofs[eptr[e]];
      const double *K = &emat[m];
      m += static_cast<std::size_t>(n) * n;
      if (n < 3 || n > 8)
      {
        continue;
      }
      Mat Ks(n, n);
      bool ok = true;
      for (int a = 0; a < n && ok; a++)
      {
        ok = d[dofs[a]] > 0.0;
      }
      if (!ok)
      {
        continue;
      }
      for (int a = 0; a < n; a++)
      {
        for (int b = 0; b < n; b++)
        {
          Ks(a, b) =
              0.5 * (K[a * n + b] + K[b * n + a]) / std::sqrt(d[dofs[a]] * d[dofs[b]]);
        }
      }
      es.compute(Ks);
      const auto &w = es.eigenvalues();
      const auto &U = es.eigenvectors();
      // The eigenvalues above the lowest spectral gap of ratio >= 1 / gap among the nonzero
      // ones (w[0] is the kernel of constants) are the strong ones.
      int first = -1;
      for (int j = 1; j + 1 < n; j++)
      {
        if (w(j + 1) > 0.0 && w(j) <= STRONG_GAP * w(j + 1))
        {
          first = j + 1;
          break;
        }
      }
      if (first < 0)
      {
        continue;
      }
      for (int a = 0; a < n; a++)
      {
        for (int b = a + 1; b < n; b++)
        {
          double p = 0.0;
          for (int j = first; j < n; j++)
          {
            p += U(a, j) * U(b, j);
          }
          if (std::abs(p) >= STRONG_TIE && gid[dofs[a]] != gid[dofs[b]])
          {
            pairs.push_back(gid[dofs[a]]);
            pairs.push_back(gid[dofs[b]]);
            tie_a.push_back(dofs[a]);
            tie_b.push_back(dofs[b]);
            tie_w.push_back(p);
          }
        }
      }
    }
  }

  auto aggregates = std::make_shared<VertexAggregates>();
  aggregates->tie_a = std::move(tie_a);
  aggregates->tie_b = std::move(tie_b);
  aggregates->tie_w = std::move(tie_w);
  const auto map = GlobalClasses(comm, pairs);
  HYPRE_BigInt num_pairs = pairs.size();
  Mpi::GlobalSum(1, &num_pairs, comm);
  aggregates->any = (num_pairs > 0);
  aggregates->aggregate.resize(nv);
  aggregates->nontrivial.resize(nv);
  for (int i = 0; i < nv; i++)
  {
    const auto it = map.find(gid[i]);
    aggregates->aggregate[i] = (it != map.end()) ? it->second : gid[i];
    aggregates->nontrivial[i] = (it != map.end());
  }
  return aggregates;
}

}  // namespace

const VertexAggregates &GetVertexAggregates(const Mesh &mesh)
{
  if (const auto *aggregates = mesh.GetVertexAggregates())
  {
    return *aggregates;
  }
  // The Laplacian on the lowest order nodal space of the mesh, partially assembled. The
  // space only reads the mesh.
  mfem::H1_FECollection fec(1, mesh.Dimension());
  FiniteElementSpace h1_fespace(const_cast<Mesh &>(mesh), &fec);
  BilinearForm k(h1_fespace);
  k.AddDomainIntegrator<DiffusionIntegrator>();
  mesh.SetVertexAggregates(ComputeVertexAggregates(h1_fespace, *k.PartialAssemble()));
  return *mesh.GetVertexAggregates();
}

}  // namespace palace
