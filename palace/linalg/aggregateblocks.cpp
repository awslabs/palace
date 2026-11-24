// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "aggregateblocks.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <metis.h>
#include "fem/fespace.hpp"
#include "linalg/amg.hpp"
#include "linalg/blockutils.hpp"
#include "linalg/rap.hpp"
#include "linalg/vertexaggregates.hpp"

namespace palace
{

using namespace blocks;

namespace
{

using SpMat = Eigen::SparseMatrix<double>;

// The tie between two L-DOFs, from the ties of the vertex aggregates: every vertex of the
// mesh entity of each of the two must be a vertex of the entity of the other or be tied to
// one of its vertices across a flat element, and the weight is the weakest of these ties
// (1 for a common vertex, 0 without a tie). DOFs of entities above each other across a
// stack of flat elements are tied, DOFs of neighboring entities along it are not.
struct DofTies
{
  std::vector<int> lent, eptr, evert;
  std::unordered_map<std::uint64_t, double> tie;

  static std::uint64_t Pair(int a, int b)
  {
    return (static_cast<std::uint64_t>(std::min(a, b)) << 32) |
           static_cast<std::uint64_t>(std::max(a, b));
  }

  DofTies(const mfem::ParFiniteElementSpace &fes, const VertexAggregates &va)
    : lent(fes.GetVSize(), -1), eptr(1, 0)
  {
    const mfem::ParMesh &mesh = *fes.GetParMesh();
    mfem::Array<int> dofs, vert;
    auto SetEntity = [&]()
    {
      if (!dofs.Size())
      {
        return;
      }
      const int e = static_cast<int>(eptr.size()) - 1;
      for (int v : vert)
      {
        evert.push_back(v);
      }
      eptr.push_back(static_cast<int>(evert.size()));
      for (int dof : dofs)
      {
        lent[(dof < 0) ? -1 - dof : dof] = e;
      }
    };
    for (int v = 0; v < mesh.GetNV(); v++)
    {
      fes.GetVertexDofs(v, dofs);
      vert.SetSize(1);
      vert[0] = v;
      SetEntity();
    }
    for (int e = 0; e < mesh.GetNEdges(); e++)
    {
      fes.GetEdgeInteriorDofs(e, dofs);
      mesh.GetEdgeVertices(e, vert);
      SetEntity();
    }
    for (int f = 0; f < mesh.GetNFaces(); f++)
    {
      fes.GetFaceInteriorDofs(f, dofs);
      mesh.GetFaceVertices(f, vert);
      SetEntity();
    }
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      fes.GetElementInteriorDofs(e, dofs);
      mesh.GetElementVertices(e, vert);
      SetEntity();
    }
    // Only the ties across a flat element count (negative projector entries), with their
    // magnitude.
    for (std::size_t k = 0; k < va.tie_w.size(); k++)
    {
      if (va.tie_w[k] < 0.0)
      {
        double &w = tie[Pair(va.tie_a[k], va.tie_b[k])];
        w = std::max(w, -va.tie_w[k]);
      }
    }
  }

  double Tie(int a, int b) const
  {
    if (a == b)
    {
      return 1.0;
    }
    const auto it = tie.find(Pair(a, b));
    return (it != tie.end()) ? it->second : 0.0;
  }

  double Directed(int ei, int ej) const
  {
    double d = 1.0;
    for (int p = eptr[ei]; p < eptr[ei + 1]; p++)
    {
      double m = 0.0;
      for (int q = eptr[ej]; q < eptr[ej + 1]; q++)
      {
        m = std::max(m, Tie(evert[p], evert[q]));
      }
      d = std::min(d, m);
    }
    return d;
  }

  // The tie of L-DOFs i and j.
  double operator()(int i, int j) const
  {
    const int ei = lent[i], ej = lent[j];
    if (ei < 0 || ej < 0)
    {
      return 0.0;
    }
    return (ei == ej) ? 1.0 : std::min(Directed(ei, ej), Directed(ej, ei));
  }
};

// The entries of a matrix between rows of a common group, symmetrized.
SpMat GroupMatrix(hypre_CSRMatrix *A, const std::vector<int> &group)
{
  const int n = static_cast<int>(group.size());
  const HYPRE_Int *I = hypre_CSRMatrixI(A), *J = hypre_CSRMatrixJ(A);
  const HYPRE_Real *V = hypre_CSRMatrixData(A);
  std::vector<Eigen::Triplet<double>> entries;
  for (int i = 0; i < n; i++)
  {
    for (HYPRE_Int e = I[i]; e < I[i + 1]; e++)
    {
      const int j = J[e];
      if (j >= 0 && j < n && group[i] == group[j])
      {
        entries.emplace_back(i, j, 0.5 * V[e]);
        entries.emplace_back(j, i, 0.5 * V[e]);
      }
    }
  }
  SpMat M(n, n);
  M.setFromTriplets(entries.begin(), entries.end());
  M.makeCompressed();
  return M;
}

// The dense blocks which a process owns, in the size classes of 1, 2, 4, ..., 64 rows. A
// block of m rows is in the class of c = the smallest power of two >= m rows and is padded
// with the identity. Per class: the inverses of its blocks, one after the other, and the
// index of every row of every block (-1 for padding) into a source vector, which is
// gathered into one block vector, multiplied class by class, and scattered back.
//
// The inverse of a block is symmetric, and its lower triangle is stored, packed by rows:
// c (c + 1) / 2 entries, 4 c (c + 1) bytes per block. Next to the full block with one dot
// product per row, the product with the packed triangle (a dot product and an update per
// row) takes about 0.6 to 0.75 of the time for c = 16, 32 and 64, where the time of the
// blocks is spent, and within 15 % of it for the smaller classes (measured with -O3 on a
// 64-bit Arm core, blocks out of cache), on half the memory.
struct ClassBlocks
{
  static constexpr int NUM_CLASSES = AggregateBlocks::NUM_CLASSES;
  std::array<std::vector<double>, NUM_CLASSES> inv;
  std::array<std::vector<int>, NUM_CLASSES> row;
  std::array<long long, NUM_CLASSES> num{};
  long long num_split = 0, num_dropped = 0, num_negative = 0;
  mutable std::vector<double> x, y;

  // The pieces of the rows [begin, end) of a group: connected pieces of at most max_size
  // rows each, in the graph of the group (rows; an edge for every pair which is coupled in
  // the matrix M). A group of at most max_size rows is split into its connected components
  // (the DOFs of collinear edges, for example, are in one group and not coupled). A larger
  // one is split by recursive bisection of its graph. A subgraph
  // which needs k = ceil(n / max_size) pieces is split in the ratio floor(k / 2) :
  // ceil(k / 2), so that the pieces come out close to n / k rows each, by the bisection of
  // METIS, which minimizes the weight of the pairs it cuts: 1 for a pair without a tie
  // across the flat direction, 1 + TIE_WEIGHT tie for a tied one, so that the pairs
  // without a tie are cut first and the weakest ties next. An aggregate is connected by
  // its ties, so a cut of it always cuts some. A part which is not connected is split
  // into its connected components. piece[i] is set to the first row of the piece of row
  // i.
  template <typename TieFunc>
  static void Split(const SpMat &M, int begin, int end, int max_size, TieFunc &&Tie,
                    std::vector<int> &piece)
  {
    constexpr double TIE_WEIGHT = 50.0;
    const idx_t nv = end - begin;
    const bool cut = (nv > max_size);
    std::vector<idx_t> xadj(nv + 1, 0), adjncy, adjwgt, part(nv, 0);
    for (int j = begin; j < end; j++)
    {
      for (SpMat::InnerIterator it(M, j); it; ++it)
      {
        const int i = static_cast<int>(it.row());
        if (i != j && i >= begin && i < end)
        {
          adjncy.push_back(static_cast<idx_t>(i - begin));
          adjwgt.push_back(
              cut ? 1 + static_cast<idx_t>(std::llround(
                            TIE_WEIGHT * std::min(1.0, std::max(0.0, Tie(i, j)))))
                  : 1);
        }
      }
      xadj[j - begin + 1] = static_cast<idx_t>(adjncy.size());
    }
    idx_t num_parts = 0;
    std::vector<std::vector<idx_t>> stack;
    stack.emplace_back(nv);
    std::iota(stack.back().begin(), stack.back().end(), 0);
    std::vector<idx_t> local(nv, -1), sx, sa, sw, sp, queue;
    while (!stack.empty())
    {
      const std::vector<idx_t> verts = std::move(stack.back());
      stack.pop_back();
      const idx_t ns = static_cast<idx_t>(verts.size());
      if (ns <= max_size)
      {
        for (idx_t v : verts)
        {
          part[v] = num_parts;
        }
        num_parts++;
        continue;
      }
      // The subgraph, and its connected components.
      for (idx_t q = 0; q < ns; q++)
      {
        local[verts[q]] = q;
      }
      sx.assign(1, 0);
      sa.clear();
      sw.clear();
      for (idx_t q = 0; q < ns; q++)
      {
        for (idx_t r = xadj[verts[q]]; r < xadj[verts[q] + 1]; r++)
        {
          if (local[adjncy[r]] >= 0)
          {
            sa.push_back(local[adjncy[r]]);
            sw.push_back(adjwgt[r]);
          }
        }
        sx.push_back(static_cast<idx_t>(sa.size()));
      }
      sp.assign(ns, -1);
      idx_t num_split = 0;
      for (idx_t q = 0; q < ns; q++)
      {
        if (sp[q] >= 0)
        {
          continue;
        }
        queue.assign(1, q);
        sp[q] = num_split;
        for (std::size_t h = 0; h < queue.size(); h++)
        {
          for (idx_t r = sx[queue[h]]; r < sx[queue[h] + 1]; r++)
          {
            if (sp[sa[r]] < 0)
            {
              sp[sa[r]] = num_split;
              queue.push_back(sa[r]);
            }
          }
        }
        num_split++;
      }
      if (num_split == 1)
      {
        const idx_t k = (ns + max_size - 1) / max_size, k0 = k / 2;
        real_t tpwgts[2] = {static_cast<real_t>(k0) / k, 1.0f};
        tpwgts[1] -= tpwgts[0];
        idx_t ncon = 1, objval = 0, two = 2, nvs = ns;
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_NUMBERING] = 0;
        const int ierr = METIS_PartGraphRecursive(&nvs, &ncon, sx.data(), sa.data(),
                                                  nullptr, nullptr, sw.data(), &two, tpwgts,
                                                  nullptr, options, &objval, sp.data());
        MFEM_VERIFY(ierr == METIS_OK, "Block Jacobi: METIS failed to split a group!");
        num_split = 2;
        const idx_t n0 = static_cast<idx_t>(std::count(sp.begin(), sp.end(), 0));
        if (n0 == 0 || n0 == ns)
        {
          // No bisection: the two parts in the order of the rows.
          for (idx_t q = 0; q < ns; q++)
          {
            sp[q] = (q < (ns * k0) / k) ? 0 : 1;
          }
        }
      }
      std::vector<std::vector<idx_t>> subs(num_split);
      for (idx_t q = 0; q < ns; q++)
      {
        subs[sp[q]].push_back(verts[q]);
        local[verts[q]] = -1;
      }
      for (auto &sub : subs)
      {
        stack.push_back(std::move(sub));
      }
    }
    // The pieces: the connected components of the parts (the parts of METIS need not be
    // connected).
    std::vector<int> parent(nv);
    std::iota(parent.begin(), parent.end(), 0);
    auto Find = [&parent](int i)
    {
      while (parent[i] != i)
      {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    };
    for (idx_t v = 0; v < nv; v++)
    {
      for (idx_t r = xadj[v]; r < xadj[v + 1]; r++)
      {
        if (part[adjncy[r]] == part[v])
        {
          const int rv = Find(static_cast<int>(v)), rw = Find(static_cast<int>(adjncy[r]));
          if (rv != rw)
          {
            parent[std::max(rv, rw)] = std::min(rv, rw);
          }
        }
      }
    }
    for (idx_t v = 0; v < nv; v++)
    {
      piece[begin + v] = begin + Find(static_cast<int>(v));
    }
  }

  // The inverse of a block from its symmetric eigenvalue decomposition.
  //
  // Without A: M = V Λ Vᵀ is the block of the operator, and the inverse is V Λ⁺ Vᵀ, the
  // pseudo-inverse. An eigenvalue is taken to be zero if |λ| <= m ε max |λ|: the computed
  // eigenvalues of a symmetric m x m matrix are the exact ones of a matrix within a
  // distance of the order m ε ||M|| (the backward error of the reduction to tridiagonal
  // form and of the QR iteration), so that a smaller eigenvalue cannot be told from zero.
  // The blocks of a curl-curl operator without a mass term have such eigenvalues, for the
  // gradients which are supported in the block.
  //
  // With A: M is the block of the definite counterpart and A the block of the operator.
  // With W = V Λ^(-1/2) on the eigenvalues of M above the threshold, Wᵀ A W = Q μ Qᵀ and
  // U = W Q, so that A U = M U μ and Uᵀ M U = I, the inverse is U sign(μ) Uᵀ. A μ at the
  // rounding level of Wᵀ A W (by the same bound) gets the sign +: the relaxed operator
  // has an eigenvalue |μ| of that size in this direction with either sign.
  static void Invert(const Eigen::MatrixXd &M, const Eigen::MatrixXd *A,
                     Eigen::MatrixXd &Minv, long long &dropped, long long &negative)
  {
    const int m = static_cast<int>(M.rows());
    constexpr double eps = std::numeric_limits<double>::epsilon();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    MFEM_VERIFY(es.info() == Eigen::Success,
                "Block Jacobi: eigenvalue decomposition of a block failed!");
    const auto &lambda = es.eigenvalues();
    const double tol = m * eps * lambda.cwiseAbs().maxCoeff();
    std::vector<int> kept;
    for (int j = 0; j < m; j++)
    {
      if (std::abs(lambda(j)) > tol)
      {
        kept.push_back(j);
      }
    }
    const int r = static_cast<int>(kept.size());
    dropped += m - r;
    Eigen::MatrixXd W(m, r);
    Eigen::VectorXd sign(r);
    for (int j = 0; j < r; j++)
    {
      W.col(j) = es.eigenvectors().col(kept[j]) / std::sqrt(std::abs(lambda(kept[j])));
      sign(j) = (lambda(kept[j]) < 0.0) ? -1.0 : 1.0;
    }
    if (A && r > 0)
    {
      MFEM_VERIFY(sign.minCoeff() > 0.0,
                  "Block Jacobi: the definite counterpart of a block is not definite!");
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ea(W.transpose() * (*A) * W);
      MFEM_VERIFY(ea.info() == Eigen::Success,
                  "Block Jacobi: eigenvalue decomposition of a block failed!");
      const auto &mu = ea.eigenvalues();
      const double tol_mu = r * eps * mu.cwiseAbs().maxCoeff();
      W = W * ea.eigenvectors();
      for (int j = 0; j < r; j++)
      {
        sign(j) = (mu(j) < -tol_mu) ? -1.0 : 1.0;
      }
    }
    negative += (sign.array() < 0.0).count();
    Minv = W * sign.asDiagonal() * W.transpose();
  }

  // The classes for the rows 0, ..., n - 1 in pieces (piece[i]: the first row of the piece
  // of row i), with the matrices M and, if not null, A (see Invert), which couple only
  // rows of a common piece.
  void Setup(const SpMat &M, const SpMat *A, const std::vector<int> &piece)
  {
    const int n = static_cast<int>(piece.size());
    // The rows piece by piece, and the position of every row in its piece.
    std::vector<int> order(n), pos(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&piece](int a, int b) { return piece[a] < piece[b]; });
    Eigen::MatrixXd Mb, Ab, Minv;
    for (int b = 0; b < n;)
    {
      int e = b;
      for (; e < n && piece[order[e]] == piece[order[b]]; e++)
      {
        pos[order[e]] = e - b;
      }
      const int m = e - b;
      MFEM_VERIFY(m <= AggregateBlocks::MAX_BLOCK_SIZE,
                  "Block Jacobi: unexpected size of a block!");
      int q = 0;
      while ((1 << q) < m)
      {
        q++;
      }
      const int c = 1 << q;
      Mb.setZero(m, m);
      Ab.setZero(A ? m : 0, A ? m : 0);
      for (int k = 0; k < m; k++)
      {
        for (SpMat::InnerIterator it(M, order[b + k]); it; ++it)
        {
          Mb(pos[it.row()], k) = it.value();
        }
        if (A)
        {
          for (SpMat::InnerIterator it(*A, order[b + k]); it; ++it)
          {
            Ab(pos[it.row()], k) = it.value();
          }
        }
      }
      Invert(Mb, A ? &Ab : nullptr, Minv, num_dropped, num_negative);
      for (int i = 0; i < c; i++)
      {
        for (int j = 0; j <= i; j++)
        {
          inv[q].push_back((i < m) ? 0.5 * (Minv(i, j) + Minv(j, i)) : (i == j));
        }
        row[q].push_back((i < m) ? order[b + i] : -1);
      }
      num[q]++;
      b = e;
    }
    std::size_t total = 0;
    for (int q = 0; q < NUM_CLASSES; q++)
    {
      total += row[q].size();
    }
    x.assign(total, 0.0);
    y.assign(total, 0.0);
  }

  // The rows of the blocks (see Setup) renumbered: row i becomes map[i].
  void Renumber(const std::vector<int> &map)
  {
    for (auto &rows : row)
    {
      for (int &i : rows)
      {
        i = (i >= 0) ? map[i] : -1;
      }
    }
  }

  // y = M⁻¹ x for the nb blocks of the class of C rows.
  template <int C>
  static void Kernel(long long nb, const double *inv, const double *x, double *y)
  {
    for (long long b = 0; b < nb; b++)
    {
      const double *L = inv + b * (C * (C + 1) / 2), *xb = x + b * C;
      double *yb = y + b * C;
      for (int i = 0; i < C; i++)
      {
        // Row i of the lower triangle: its part of y_i, and its part of y_j, j < i, from
        // the entries of the upper triangle.
        const double xi = xb[i];
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
          sum += L[j] * xb[j];
          yb[j] += L[j] * xi;
        }
        yb[i] = sum + L[i] * xi;
        L += i + 1;
      }
    }
  }

  // Z[row] = M⁻¹ X[row] for the rows of the blocks: one gather, one kernel per class, one
  // scatter.
  void Apply(const double *X, double *Z) const
  {
    std::size_t o = 0;
    for (int q = 0; q < NUM_CLASSES; q++)
    {
      for (const int i : row[q])
      {
        x[o++] = (i >= 0) ? X[i] : 0.0;
      }
    }
    const double *xq = x.data();
    double *yq = y.data();
    Kernel<1>(num[0], inv[0].data(), xq, yq);
    xq += row[0].size(), yq += row[0].size();
    Kernel<2>(num[1], inv[1].data(), xq, yq);
    xq += row[1].size(), yq += row[1].size();
    Kernel<4>(num[2], inv[2].data(), xq, yq);
    xq += row[2].size(), yq += row[2].size();
    Kernel<8>(num[3], inv[3].data(), xq, yq);
    xq += row[3].size(), yq += row[3].size();
    Kernel<16>(num[4], inv[4].data(), xq, yq);
    xq += row[4].size(), yq += row[4].size();
    Kernel<32>(num[5], inv[5].data(), xq, yq);
    xq += row[5].size(), yq += row[5].size();
    Kernel<64>(num[6], inv[6].data(), xq, yq);
    o = 0;
    for (int q = 0; q < NUM_CLASSES; q++)
    {
      for (const int i : row[q])
      {
        if (i >= 0)
        {
          Z[i] = y[o];
        }
        o++;
      }
    }
  }
};

// Approximate solves of the groups of a nodal space which are too large for a block
// ("sheets"). They stay on the processes of their DOFs: their matrix A_b (the block
// diagonal matrix of all of them, a parallel sparse matrix on their DOFs, symmetric
// positive definite) is kept, and A_b⁻¹ r is replaced by a fixed number of steps of the
// Chebyshev iteration for B A_b on [a, b], with B one BoomerAMG V-cycle of A_b. The result
// S is a fixed symmetric positive definite operator. Its memory is that of A_b and of the
// AMG hierarchy.
//
// The interval is an estimate. A Lanczos (preconditioned CG) run for a random right-hand
// side gives Ritz values θ₁ ≤ ... ≤ θ_m of B A_b, each with a residual bound δ: an
// eigenvalue of B A_b lies within δ of θ. That locates an eigenvalue, not necessarily the
// smallest or the largest one (the run starts from A_b v, in which the weak modes are
// small). The interval is [θ₁ - δ₁, θ_m + δ_m]. If the spectrum extends below its lower
// end, the error polynomial of the Chebyshev iteration is still between 0 and 1 there:
// these components are reduced by less than the target, not amplified, and S stays
// symmetric positive definite. (The upper end is less exposed: the largest Ritz value
// converges first, and for a symmetric V-cycle with convergent smoothers the eigenvalues
// of B A_b are at most 1.)
//
// The degree targets a reduction on that interval. The sheet solve is one step of a
// smoothing sweep next to the Chebyshev polynomial of order k of the relaxed operator,
// whose guarantee for the components it is responsible for is a reduction by
// 1 / (2k + 1) (the bound of the 4th-kind polynomial, which is used for the polynomial of
// the 1st kind as well). The components on a sheet are the responsibility of the sheet
// solve: with the same reduction they are as small after the sweep as the others, and a
// more accurate solve does not make the sweep a better smoother. This is an argument for
// the size of the reduction, not a bound on the convergence of the multigrid cycle. The
// Chebyshev iteration of degree d on [a, b] reduces the error in the A_b-norm by
// 2 ρᵈ / (1 + ρ²ᵈ), ρ = (√κ - 1) / (√κ + 1), κ = b / a, and the degree is the smallest d
// for which this is at most the required reduction r (AggregateBlocks::SheetDegree).
//
// The work is bounded. One AMG cycle is taken to be a preconditioner of the sheets as
// long as it leaves a condition number κ of at most 1 / r², which is (2k + 1)²: the range
// λ_max / λ of the eigenvalues for which the bound of the polynomial smoother of order k
// itself promises a reduction. Two limits follow from it (AggregateBlocks::GetSheetLimits).
// The degree is at most the one for κ = 1 / r² (13 for k = 4). The Lanczos run has at
// most the number of steps after which the error bound of CG for that κ, 2 ρᵐ, is at the
// rounding level (165 for k = 4): a run which has not converged by then has a larger κ.
// Beyond the bound (a larger estimated κ, or a run which has not converged), the
// polynomial is the one of the largest degree, for the interval [b r², b]: a fixed
// symmetric positive definite operator as before, which reaches the required reduction
// on that part of the spectrum and less below it. The smoother reports when this limit
// is used.
struct SheetSolve
{
  MPI_Comm comm = MPI_COMM_NULL;
  std::vector<int> dofs;  // Local true DOFs of the sheets
  HYPRE_BigInt num_global = 0;
  std::unique_ptr<mfem::HypreParMatrix> A;
  std::unique_ptr<mfem::HypreBoomerAMG> amg;
  double a = 0.0, b = 0.0;
  int degree = 0;
  bool limited = false;      // The degree is the largest one (see above)
  HYPRE_BigInt entries = 0;  // Of A_b and of the matrices of the AMG hierarchy (global)
  mutable Vector r, x, d, w, z;

  int Size() const { return static_cast<int>(dofs.size()); }

  double Dot(const Vector &u, const Vector &v) const
  {
    double dot = linalg::LocalDot(u, v);
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  }

  // Returns whether A_b is positive definite, as far as the Lanczos run can tell: the
  // iteration is not defined otherwise.
  bool Setup(double reduction)
  {
    const auto limits = AggregateBlocks::GetSheetLimits(reduction);
    const int n = Size();
    for (Vector *v : {&r, &x, &d, &w, &z})
    {
      v->SetSize(n);
      v->UseDevice(false);
      *v = 0.0;
    }
    Vector p(n), Ap(n);
    p.UseDevice(false);
    Ap.UseDevice(false);
    {
      std::minstd_rand gen(1234 + Mpi::Rank(comm));
      std::uniform_real_distribution<double> dist(-1.0, 1.0);
      auto *P = p.HostWrite();
      for (int i = 0; i < n; i++)
      {
        P[i] = dist(gen);
      }
    }
    A->Mult(p, r);
    amg = std::make_unique<BoomerAmgSolver>(1, 1, 0, 0);
    amg->SetOperator(*A);
    HYPRE_BoomerAMGSetCumNnzAP(*amg, 1.0);  // Count the entries of the hierarchy

    // Lanczos for B A_b: CG on A_b u = A_b v with a random v. For a definite A_b (and B),
    // every step has positive curvature. The interval of the Ritz values, widened by
    // their residual bounds, is examined every CHECK steps, until the bounds are within
    // INTERVAL_TOL of the two ends (the condition number, and with it the degree, is then
    // overestimated by at most this much), or the residual is at the rounding level, or
    // the number of steps is at its limit.
    amg->Mult(r, z);
    p = z;
    double rz = Dot(r, z);
    const double rz0 = rz;
    std::vector<double> alphas, betas;
    constexpr int CHECK = 10;
    constexpr double INTERVAL_TOL = 0.1;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    auto Interval = [&]()
    {
      // The tridiagonal matrix of the Lanczos method from the coefficients of CG.
      const int k = static_cast<int>(alphas.size());
      Eigen::VectorXd diag(k), subdiag(std::max(k - 1, 0));
      for (int i = 0; i < k; i++)
      {
        diag(i) = 1.0 / alphas[i] + ((i > 0) ? betas[i - 1] / alphas[i - 1] : 0.0);
        if (i + 1 < k)
        {
          subdiag(i) = std::sqrt(betas[i]) / alphas[i];
        }
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      es.computeFromTridiagonal(diag, subdiag);
      // The residual of a Ritz pair (θ, s) of the Lanczos relation is the next
      // off-diagonal entry times the last entry of s.
      const double next = (static_cast<int>(betas.size()) == k)
                              ? std::sqrt(betas[k - 1]) / alphas[k - 1]
                              : 0.0;
      const double da = next * std::abs(es.eigenvectors()(k - 1, 0)),
                   db = next * std::abs(es.eigenvectors()(k - 1, k - 1));
      a = es.eigenvalues()(0) - da;
      b = es.eigenvalues()(k - 1) + db;
      return da <= INTERVAL_TOL * es.eigenvalues()(0) &&
             db <= INTERVAL_TOL * es.eigenvalues()(k - 1);
    };
    bool converged = false;
    for (int it = 0; it < limits.lanczos_steps && !converged; it++)
    {
      A->Mult(p, Ap);
      const double pAp = Dot(p, Ap);
      if (!(pAp > 0.0))
      {
        return false;
      }
      const double alpha = rz / pAp;
      r.Add(-alpha, Ap);
      amg->Mult(r, z);
      const double rz_new = Dot(r, z);
      alphas.push_back(alpha);
      if (!(rz_new >= 0.0))
      {
        return false;
      }
      if (rz_new <= eps * eps * rz0)
      {
        converged = Interval();
        break;
      }
      betas.push_back(rz_new / rz);
      linalg::AXPBY(1.0, z, betas.back(), p);
      rz = rz_new;
      converged = ((it + 1) % CHECK == 0) && Interval();
    }
    if (!converged)
    {
      converged = Interval();
    }
    if (!(b > 0.0) || !std::isfinite(b))
    {
      return false;
    }
    limited = !converged || !(limits.condition * a >= b);
    if (limited)
    {
      a = b / limits.condition;
    }
    degree = AggregateBlocks::SheetDegree(a, b, reduction);
    HYPRE_Real cum_nnz = 0.0;
    HYPRE_BoomerAMGGetCumNnzAP(*amg, &cum_nnz);
    entries = static_cast<HYPRE_BigInt>(cum_nnz);
    return true;
  }

  // Y[dofs] += q(B A_b) B X[dofs].
  void AddMult(const double *X, double *Y) const
  {
    const int n = Size();
    const double theta = 0.5 * (b + a), delta = 0.5 * (b - a);
    {
      auto *R = r.HostWrite();
      for (int i = 0; i < n; i++)
      {
        R[i] = X[dofs[i]];
      }
    }
    x = 0.0;
    amg->Mult(r, d);
    d *= 1.0 / theta;
    double rho = delta / theta;
    for (int k = 0; k < degree; k++)
    {
      x += d;
      if (k == degree - 1)
      {
        break;
      }
      A->Mult(d, w);
      r -= w;
      const double rho_new = 1.0 / (2.0 * theta / delta - rho);
      amg->Mult(r, z);
      linalg::AXPBY(2.0 * rho_new / delta, z, rho_new * rho, d);
      rho = rho_new;
    }
    const auto *XX = x.HostRead();
    for (int i = 0; i < n; i++)
    {
      Y[dofs[i]] += XX[i];
    }
  }
};

// Step 1: the column of every L-DOF, the key of the set of aggregates of the vertices of
// the mesh entity which carries it, and whether one of these vertices is in an aggregate
// of several vertices (in flat).
std::vector<Key> EntityColumns(const mfem::ParFiniteElementSpace &fes,
                               const VertexAggregates &va, std::vector<char> &flat)
{
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  std::vector<Key> lkey(fes.GetVSize(), 0);
  flat.assign(fes.GetVSize(), 0);
  mfem::Array<int> dofs, vert;
  std::vector<Key> ids;
  auto SetEntity = [&]()
  {
    if (!dofs.Size())
    {
      return;
    }
    ids.clear();
    bool nontrivial = false;
    for (int v : vert)
    {
      ids.push_back(va.aggregate[v]);
      nontrivial = nontrivial || va.nontrivial[v];
    }
    const Key k = SetKey(ids);
    for (int dof : dofs)
    {
      const int l = (dof < 0) ? -1 - dof : dof;
      lkey[l] = k;
      flat[l] = nontrivial;
    }
  };
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    fes.GetVertexDofs(v, dofs);
    vert.SetSize(1);
    vert[0] = v;
    SetEntity();
  }
  for (int e = 0; e < mesh.GetNEdges(); e++)
  {
    fes.GetEdgeInteriorDofs(e, dofs);
    mesh.GetEdgeVertices(e, vert);
    SetEntity();
  }
  for (int f = 0; f < mesh.GetNFaces(); f++)
  {
    fes.GetFaceInteriorDofs(f, dofs);
    mesh.GetFaceVertices(f, vert);
    SetEntity();
  }
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    fes.GetElementInteriorDofs(e, dofs);
    mesh.GetElementVertices(e, vert);
    SetEntity();
  }
  return lkey;
}

// Step 3, for nonconforming meshes: within a column the strong part of the operator
// penalizes the differences of the L-DOF values, which are the rows of P. Summed over a
// column of true DOFs, the weights of the rows of a column of L-DOFs must agree, or else
// the flat vectors of the columns involved are tied by a strong term: the columns of true
// DOFs whose weight is not the same in every row of a column of L-DOFs over flat elements
// are merged. Returns the map from the merged columns to their representative. slave marks
// the L-DOFs which are not the L-DOF of a true DOF of this process.
std::unordered_map<Key, Key> MergeColumns(MPI_Comm comm, const PRows &P,
                                          const std::vector<Key> &lkey,
                                          const std::vector<char> &flat,
                                          const std::vector<Key> &tkey,
                                          const std::vector<char> &slave)
{
  const int rank = Mpi::Rank(comm);
  const std::vector<Key> okey = P.Halo(tkey);
  // The rows (column of the L-DOF, L-DOF, column of true DOFs, weight) of the constrained
  // L-DOFs, and the columns with an unconstrained one (whose row is the unit row).
  std::vector<Key> row_col, row_id, row_blk, true_cols;
  std::vector<double> row_w;
  std::vector<std::pair<Key, double>> row;
  for (int i = 0; i < P.nrows; i++)
  {
    if (!flat[i])
    {
      continue;
    }
    if (!slave[i])
    {
      true_cols.push_back(lkey[i]);
      continue;
    }
    row.clear();
    P.ForRow(i,
             [&](int c, bool offd, double w)
             {
               const Key k = offd ? okey[c] : tkey[c];
               if (k)
               {
                 row.emplace_back(k, w);
               }
             });
    if (row.size() == 1 && row[0].first == lkey[i] && std::abs(row[0].second - 1.0) < 1e-12)
    {
      // A shared L-DOF of a true DOF of the same column on another process.
      true_cols.push_back(lkey[i]);
      continue;
    }
    std::sort(row.begin(), row.end());
    for (std::size_t a = 0; a < row.size();)
    {
      double w = 0.0;
      std::size_t b = a;
      for (; b < row.size() && row[b].first == row[a].first; b++)
      {
        w += row[b].second;
      }
      row_col.push_back(lkey[i]);
      row_id.push_back((static_cast<Key>(rank) << 32) | static_cast<Key>(i));
      row_blk.push_back(row[a].first);
      row_w.push_back(w);
      a = b;
    }
  }
  std::sort(true_cols.begin(), true_cols.end());
  true_cols.erase(std::unique(true_cols.begin(), true_cols.end()), true_cols.end());
  std::vector<int> counts;
  const std::vector<Key> all_col = GatherRoot(comm, row_col, counts);
  const std::vector<Key> all_id = GatherRoot(comm, row_id, counts);
  const std::vector<Key> all_blk = GatherRoot(comm, row_blk, counts);
  const std::vector<double> all_w = GatherRoot(comm, row_w, counts);
  const std::vector<Key> all_true = GatherRoot(comm, true_cols, counts);
  std::vector<Key> classes;
  if (rank == 0)
  {
    constexpr double tol = 1.0e-9;
    const std::unordered_set<Key> has_true(all_true.begin(), all_true.end());
    const int num_rec = static_cast<int>(all_col.size());
    std::vector<int> order(num_rec);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b)
              {
                return (all_col[a] != all_col[b]) ? all_col[a] < all_col[b]
                       : (all_id[a] != all_id[b]) ? all_id[a] < all_id[b]
                                                  : all_blk[a] < all_blk[b];
              });
    KeyUnion merge;
    std::vector<std::pair<Key, double>> ref, cur;
    std::vector<Key> differ;
    for (int a = 0; a < num_rec;)
    {
      // The rows of one column of L-DOFs, compared with its first one (the unit row if the
      // column has an unconstrained L-DOF).
      const Key col = all_col[order[a]];
      int end = a;
      while (end < num_rec && all_col[order[end]] == col)
      {
        end++;
      }
      ref.clear();
      differ.clear();
      bool first = true;
      if (has_true.count(col))
      {
        ref.emplace_back(col, 1.0);
        first = false;
      }
      for (int b = a; b < end;)
      {
        const Key id = all_id[order[b]];
        cur.clear();
        for (; b < end && all_id[order[b]] == id; b++)
        {
          cur.emplace_back(all_blk[order[b]], all_w[order[b]]);
        }
        if (first)
        {
          ref = cur;
          first = false;
          continue;
        }
        // Both are sorted by the key of the column of true DOFs.
        for (std::size_t i = 0, j = 0; i < cur.size() || j < ref.size();)
        {
          if (j == ref.size() || (i < cur.size() && cur[i].first < ref[j].first))
          {
            if (std::abs(cur[i].second) > tol)
            {
              differ.push_back(cur[i].first);
            }
            i++;
          }
          else if (i == cur.size() || ref[j].first < cur[i].first)
          {
            if (std::abs(ref[j].second) > tol)
            {
              differ.push_back(ref[j].first);
            }
            j++;
          }
          else
          {
            if (std::abs(cur[i].second - ref[j].second) > tol)
            {
              differ.push_back(cur[i].first);
            }
            i++;
            j++;
          }
        }
      }
      for (std::size_t k = 1; k < differ.size(); k++)
      {
        merge.Union(differ[0], differ[k]);
      }
      a = end;
    }
    classes = merge.Classes();
  }
  BroadcastRoot(comm, classes);
  std::unordered_map<Key, Key> merged;
  for (std::size_t k = 0; k + 1 < classes.size(); k += 2)
  {
    merged.emplace(classes[k], classes[k + 1]);
  }
  return merged;
}

// Step 4: what becomes of the DOFs of a block on this process.
struct BlockInfo
{
  enum Kind : int
  {
    NONE,   // Not a group: a single DOF
    DENSE,  // On its owner (gathered there if need be), in dense blocks
    SHEET   // Solved approximately as a whole where its DOFs are
  };
  Kind kind = NONE;
  HYPRE_BigInt slot = 0;        // First global slot of the DOFs of this process
  HYPRE_BigInt block_slot = 0;  // First global slot of the block
  int size = 0;                 // Global number of DOFs of the group
  int owner = 0;                // Process which owns the group: the one with most DOFs
};

// The census of the blocks, on the root: for the blocks (keys) with DOFs on this process
// and their numbers of DOFs here, their BlockInfo, the first global slot of every process
// in slot_starts, and the global counts.
std::unordered_map<Key, BlockInfo> CountBlocks(MPI_Comm comm, const std::vector<Key> &keys,
                                               const std::vector<int> &sizes,
                                               int max_dense_size, bool sheets,
                                               std::vector<HYPRE_BigInt> &slot_starts,
                                               AggregateBlocks::Census &census)
{
  const int rank = Mpi::Rank(comm), np = Mpi::Size(comm);
  constexpr int num_fields = 5;
  std::vector<int> counts;
  const std::vector<Key> all_keys = GatherRoot(comm, keys, counts);
  const std::vector<int> all_sizes = GatherRoot(comm, sizes, counts);
  std::vector<long long> reply(num_fields * all_keys.size(), 0), starts(np + 1, 0);
  // Blocks and DOFs of each kind, and the largest factored block.
  long long totals[7] = {0, 0, 0, 0, 0, 0, 0};  // totals[6]: groups of several processes
  if (rank == 0)
  {
    // The contributions to every block, in rank order.
    const int num_rec = static_cast<int>(all_keys.size());
    std::vector<int> all_rank(num_rec), order(num_rec);
    for (int r = 0, q = 0; r < np; r++)
    {
      for (int k = 0; k < counts[r]; k++)
      {
        all_rank[q++] = r;
      }
    }
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b)
              {
                return (all_keys[a] != all_keys[b]) ? all_keys[a] < all_keys[b]
                                                    : all_rank[a] < all_rank[b];
              });
    // First pass: size, owner and kind of every block, and the slots of every owner.
    struct Block
    {
      int begin, end, owner, size;
      BlockInfo::Kind kind;
    };
    std::vector<Block> blocks;
    for (int a = 0; a < num_rec;)
    {
      int b = a, size = 0, owner = all_rank[order[a]], most = -1;
      for (; b < num_rec && all_keys[order[b]] == all_keys[order[a]]; b++)
      {
        size += all_sizes[order[b]];
        if (all_sizes[order[b]] > most)
        {
          most = all_sizes[order[b]];
          owner = all_rank[order[b]];
        }
      }
      const auto kind = (size <= 1)                ? BlockInfo::NONE
                        : (size <= max_dense_size) ? BlockInfo::DENSE
                        : sheets                   ? BlockInfo::SHEET
                                                   : BlockInfo::NONE;
      blocks.push_back({a, b, owner, size, kind});
      if (kind == BlockInfo::DENSE)
      {
        starts[owner + 1] += size;
        totals[6] += (b - a > 1);
      }
      if (size > 1)
      {
        totals[2 * kind]++;
        totals[2 * kind + 1] += size;
      }
      a = b;
    }
    for (int r = 0; r < np; r++)
    {
      starts[r + 1] += starts[r];
    }
    std::vector<long long> next(starts.begin(), starts.end() - 1);
    for (const auto &block : blocks)
    {
      long long slot = 0, block_slot = 0;
      if (block.kind == BlockInfo::DENSE)
      {
        slot = block_slot = next[block.owner];
        next[block.owner] += block.size;
      }
      for (int q = block.begin; q < block.end; q++)
      {
        long long *r = &reply[num_fields * static_cast<std::size_t>(order[q])];
        r[0] = block.kind;
        r[1] = slot;
        r[2] = block_slot;
        r[3] = block.size;
        r[4] = block.owner;
        slot += all_sizes[order[q]];
      }
    }
  }
  const std::vector<long long> my_reply =
      ScatterRoot(comm, reply, counts, num_fields, static_cast<int>(keys.size()));
  MPI_Bcast(starts.data(), np + 1, MPI_LONG_LONG, 0, comm);
  MPI_Bcast(totals, 7, MPI_LONG_LONG, 0, comm);
  slot_starts.assign(starts.begin(), starts.end());
  census.skipped = totals[2 * BlockInfo::NONE];
  census.skipped_dofs = totals[2 * BlockInfo::NONE + 1];
  census.block_dofs = totals[2 * BlockInfo::DENSE + 1];
  census.sheets = totals[2 * BlockInfo::SHEET];
  census.sheet_dofs = totals[2 * BlockInfo::SHEET + 1];
  census.gathered = totals[6];
  std::unordered_map<Key, BlockInfo> info;
  for (std::size_t k = 0; k < keys.size(); k++)
  {
    const long long *r = &my_reply[num_fields * k];
    info.emplace(keys[k],
                 BlockInfo{static_cast<BlockInfo::Kind>(r[0]),
                           static_cast<HYPRE_BigInt>(r[1]), static_cast<HYPRE_BigInt>(r[2]),
                           static_cast<int>(r[3]), static_cast<int>(r[4])});
  }
  return info;
}

// Steps 5 and 6: the entries of the element matrices of a local operator between L-DOFs
// of a common block, summed into a sparse matrix A_f, and Pᵀ A_f P: the operator on the
// pairs of true DOFs of a common block. lkeys[lkeys_ptr[i]], ..., lkeys[lkeys_ptr[i + 1]
// - 1] are the blocks (sorted) of the true DOFs which L-DOF i depends on.
std::unique_ptr<mfem::HypreParMatrix>
AssembleBlockEntries(const Operator &local_op, const FiniteElementSpace &fespace,
                     const std::vector<int> &lkeys_ptr, const std::vector<Key> &lkeys)
{
  const int nl = fespace.GetVSize();
  auto Match = [&](int i, int j)
  {
    for (int p = lkeys_ptr[i], q = lkeys_ptr[j];
         p < lkeys_ptr[i + 1] && q < lkeys_ptr[j + 1];)
    {
      if (lkeys[p] == lkeys[q])
      {
        return true;
      }
      (lkeys[p] < lkeys[q]) ? p++ : q++;
    }
    return false;
  };
  std::vector<int> ti, tj;
  std::vector<double> tv;
  ForEachLocalEntry(local_op,
                    [&](int row, int col, double v)
                    {
                      if (Match(row, col))
                      {
                        ti.push_back(row);
                        tj.push_back(col);
                        tv.push_back(v);
                      }
                    });
  // Sort the entries by row and column, and sum the duplicates.
  const std::size_t num_entries = ti.size();
  std::vector<int> I(nl + 1, 0);
  for (std::size_t k = 0; k < num_entries; k++)
  {
    I[ti[k] + 1]++;
  }
  for (int i = 0; i < nl; i++)
  {
    I[i + 1] += I[i];
  }
  std::vector<int> order(num_entries), fill(I.begin(), I.end() - 1);
  for (std::size_t k = 0; k < num_entries; k++)
  {
    order[fill[ti[k]]++] = static_cast<int>(k);
  }
  std::vector<HYPRE_Int> cI(nl + 1, 0), cJ;
  std::vector<HYPRE_Real> cA;
  for (int i = 0; i < nl; i++)
  {
    std::sort(order.begin() + I[i], order.begin() + I[i + 1],
              [&](int a, int b) { return tj[a] < tj[b]; });
    for (int q = I[i]; q < I[i + 1]; q++)
    {
      if (q > I[i] && tj[order[q]] == tj[order[q - 1]])
      {
        cA.back() += tv[order[q]];
      }
      else
      {
        cJ.push_back(tj[order[q]]);
        cA.push_back(tv[order[q]]);
      }
    }
    cI[i + 1] = static_cast<HYPRE_Int>(cJ.size());
  }
  // The local matrix, in the memory of hypre, as the diagonal block of a parallel matrix
  // (as in ParOperator::ParallelAssemble).
  hypre::HypreCSRMatrix sA(nl, nl, static_cast<int>(cJ.size()));
  {
    const auto location = hypre_CSRMatrixMemoryLocation(static_cast<hypre_CSRMatrix *>(sA));
    hypre_TMemcpy(sA.GetI(), cI.data(), HYPRE_Int, cI.size(), location, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(sA.GetJ(), cJ.data(), HYPRE_Int, cJ.size(), location, HYPRE_MEMORY_HOST);
    hypre_TMemcpy(sA.GetData(), cA.data(), HYPRE_Real, cA.size(), location,
                  HYPRE_MEMORY_HOST);
  }
  hypre_ParCSRMatrix *hA = hypre_ParCSRMatrixCreate(
      fespace.GetComm(), fespace.GlobalVSize(), fespace.GlobalVSize(),
      fespace.GetDofOffsets(), fespace.GetDofOffsets(), 0, sA.NNZ(), 0);
  hypre_CSRMatrix *hA_diag = hypre_ParCSRMatrixDiag(hA);
  hypre_ParCSRMatrixDiag(hA) = sA;
  hypre_ParCSRMatrixInitialize(hA);
  const mfem::HypreParMatrix *P = fespace.Dof_TrueDof_Matrix();
  auto PtAP = std::make_unique<mfem::HypreParMatrix>(
      hypre_ParCSRMatrixRAPKT(*P, hA, *P, 1, 1), true);
  hypre_ParCSRMatrixDiag(hA) = hA_diag;
  hypre_ParCSRMatrixDestroy(hA);
  return PtAP;
}

// The selection matrix (true DOFs x columns) with a unit entry at (dofs[k], cols[k]), for
// columns [col_start, col_end) of this process out of num_cols.
mfem::HypreParMatrix SelectionMatrix(const FiniteElementSpace &fespace,
                                     const std::vector<int> &dofs,
                                     std::vector<HYPRE_BigInt> cols, HYPRE_BigInt col_start,
                                     HYPRE_BigInt col_end, HYPRE_BigInt num_cols)
{
  const int nt = fespace.GetTrueVSize();
  std::vector<int> I(nt + 1, 0);
  for (int dof : dofs)
  {
    I[dof + 1] = 1;
  }
  for (int t = 0; t < nt; t++)
  {
    I[t + 1] += I[t];
  }
  std::vector<double> data(std::max<std::size_t>(cols.size(), 1), 1.0);
  if (cols.empty())
  {
    cols.push_back(0);
    data[0] = 0.0;
  }
  HYPRE_BigInt row_range[2] = {fespace.GetTrueDofOffsets()[0],
                               fespace.GetTrueDofOffsets()[1]};
  HYPRE_BigInt col_range[2] = {col_start, col_end};
  return mfem::HypreParMatrix(fespace.GetComm(), nt, fespace.GlobalTrueVSize(), num_cols,
                              I.data(), cols.data(), data.data(), row_range, col_range);
}

hypre_CSRMatrix *DiagOf(mfem::HypreParMatrix *M)
{
  if (!M)
  {
    return nullptr;
  }
  M->HostRead();
  return hypre_ParCSRMatrixDiag(static_cast<hypre_ParCSRMatrix *>(*M));
}

}  // namespace

struct AggregateBlocks::Impl
{
  int nt = 0;
  Census census;
  ClassBlocks blocks;                       // Blocks owned by this process
  std::unique_ptr<mfem::HypreParMatrix> R;  // Slots x true DOFs (blocks of this process
                                            // with DOFs of others: gathered here)
  std::vector<int> dofs;                    // Local true DOFs in blocks
  SheetSolve sheets;
  mutable Vector g, y, zt;
};

AggregateBlocks::~AggregateBlocks() = default;

bool AggregateBlocks::IsSupported(const ParOperator &A)
{
  const auto &fespace = A.TrialFiniteElementSpace();
  return fespace.Dimension() == 3 && !fespace.IsBroken();
}

AggregateBlocks::SheetLimits AggregateBlocks::GetSheetLimits(double reduction)
{
  MFEM_VERIFY(reduction > 0.0 && reduction < 1.0,
              "Invalid reduction for the limits of the sheet polynomial!");
  SheetLimits limits;
  limits.condition = 1.0 / (reduction * reduction);
  limits.degree = SheetDegree(1.0, limits.condition, reduction);
  // 2 ρᵐ ≤ ε for ρ = (√κ - 1) / (√κ + 1) = (1 - r) / (1 + r).
  const double rho = (1.0 - reduction) / (1.0 + reduction);
  limits.lanczos_steps = static_cast<int>(std::ceil(
      std::log(2.0 / std::numeric_limits<double>::epsilon()) / std::log(1.0 / rho)));
  return limits;
}

int AggregateBlocks::SheetDegree(double a, double b, double reduction)
{
  MFEM_VERIFY(a > 0.0 && b >= a && reduction > 0.0,
              "Invalid interval or reduction for the degree of the sheet polynomial!");
  const double sqrt_kappa = std::sqrt(b / a);
  const double rho = (sqrt_kappa - 1.0) / (sqrt_kappa + 1.0);
  int degree = 1;
  for (double rho_d = rho; 2.0 * rho_d / (1.0 + rho_d * rho_d) > reduction; rho_d *= rho)
  {
    degree++;
  }
  return degree;
}

AggregateBlocks::AggregateBlocks(const ParOperator &A, const ParOperator *A_def,
                                 int block_size, double sheet_reduction)
  : impl(std::make_unique<Impl>())
{
  const FiniteElementSpace &fespace = A.TrialFiniteElementSpace();
  const mfem::ParFiniteElementSpace &fes = fespace.Get();
  MPI_Comm comm = fes.GetComm();
  MFEM_VERIFY(IsSupported(A), "The blocks of the block Jacobi smoother require a 3D mesh!");
  MFEM_VERIFY(block_size >= 2 && block_size <= MAX_BLOCK_SIZE &&
                  (block_size & (block_size - 1)) == 0,
              "The block size of the block Jacobi smoother is a power of two up to "
                  << MAX_BLOCK_SIZE << "!");
  const int nl = fes.GetVSize(), nt = fes.GetTrueVSize();
  const int rank = Mpi::Rank(comm);
  const auto *fec = fes.FEColl();
  const bool is_h1 = (fec->GetContType() == mfem::FiniteElementCollection::CONTINUOUS);
  impl->nt = nt;
  auto &census = impl->census;

  // 1. Columns of the L-DOFs.
  std::vector<char> lflat;
  const VertexAggregates &aggregates = GetVertexAggregates(fespace.GetMesh());
  if (!aggregates.any)
  {
    return;
  }
  const std::vector<Key> lkey = EntityColumns(fes, aggregates, lflat);

  // 2. Columns of the true DOFs. The L-DOF of a true DOF is the one selected by the
  //    restriction matrix. Every other L-DOF is a constrained ("slave") row of P for step
  //    3, including the unit rows: a hanging vertex at the midpoint of a master edge
  //    coincides with the order-2 DOF of that edge, but belongs to another column (and a
  //    shared L-DOF is a row of the owner's column).
  PRows P(fespace.Dof_TrueDof_Matrix());
  MFEM_VERIFY(P.nrows == nl && P.ncols == nt, "Block Jacobi: prolongation size mismatch!");
  std::vector<Key> tkey(nt, 0);
  std::vector<char> tflat(nt, 0), slave(nl, 1);
  std::vector<int> tldof(nt);
  {
    const mfem::SparseMatrix *Rt = fespace.GetRestrictionMatrix();
    MFEM_VERIFY(Rt && Rt->Height() == nt && Rt->Width() == nl,
                "Block Jacobi: unexpected restriction matrix!");
    const int *RI = Rt->HostReadI(), *RJ = Rt->HostReadJ();
    for (int t = 0; t < nt; t++)
    {
      MFEM_VERIFY(RI[t + 1] == RI[t] + 1,
                  "Block Jacobi: unexpected restriction matrix row!");
      tldof[t] = RJ[RI[t]];
      slave[tldof[t]] = 0;
      tkey[t] = lkey[tldof[t]];
      tflat[t] = lflat[tldof[t]];
    }
  }

  // 3. Columns merged at hanging entities, and the group of every true DOF (0: none): its
  //    column if that is over flat elements, or merged with one which is.
  if (fes.GetParMesh()->Nonconforming())
  {
    const auto merged = MergeColumns(comm, P, lkey, lflat, tkey, slave);
    for (int t = 0; t < nt; t++)
    {
      if (const auto it = merged.find(tkey[t]); it != merged.end())
      {
        tkey[t] = it->second;
        tflat[t] = 1;
      }
    }
  }
  for (int t = 0; t < nt; t++)
  {
    if (!tflat[t])
    {
      tkey[t] = 0;
    }
  }
  if (const auto *ess_list = A.GetEssentialTrueDofs())
  {
    const auto *list = ess_list->HostRead();
    for (int k = 0; k < ess_list->Size(); k++)
    {
      tkey[list[k]] = 0;
    }
  }

  // 4. Census of the groups. A group which fits into a block is gathered on the process
  //    with most of its DOFs. A larger group of a nodal space is a sheet. Of a larger
  //    group of another space, every process keeps the DOFs it owns as a group of its
  //    own, to be split there: a process boundary is one more cut of a group which is cut
  //    anyway, and nothing of it is gathered.
  std::unordered_map<Key, BlockInfo> info;
  std::vector<HYPRE_BigInt> slot_starts;
  {
    std::vector<Key> keys;
    std::vector<int> sizes;
    auto Count = [&]()
    {
      std::unordered_map<Key, int> count;
      for (Key k : tkey)
      {
        if (k)
        {
          count[k]++;
        }
      }
      keys.clear();
      sizes.clear();
      for (const auto &[k, n] : count)
      {
        keys.push_back(k);
      }
      std::sort(keys.begin(), keys.end());
      for (Key k : keys)
      {
        sizes.push_back(count[k]);
      }
    };
    Count();
    if (!is_h1)
    {
      Census whole;
      std::vector<HYPRE_BigInt> unused;
      const auto sizes_of = CountBlocks(comm, keys, sizes, std::numeric_limits<int>::max(),
                                        false, unused, whole);
      const Key salt = Combine(0x243f6a8885a308d3ULL, static_cast<Key>(rank));
      for (const auto &[k, group] : sizes_of)
      {
        census.split += (group.size > block_size && group.owner == rank);
      }
      for (Key &k : tkey)
      {
        if (k && sizes_of.at(k).size > block_size)
        {
          k = Combine(k, salt);
          k += (k == 0);
        }
      }
      Mpi::GlobalSum(1, &census.split, comm);
      Count();
    }
    const HYPRE_BigInt split = census.split;
    info =
        CountBlocks(comm, keys, sizes, is_h1 ? block_size : std::numeric_limits<int>::max(),
                    is_h1, slot_starts, census);
    census.split = split;
  }
  // The slot of every true DOF of a group in blocks, and the DOFs of the sheets.
  std::vector<HYPRE_BigInt> slots;
  for (int t = 0; t < nt; t++)
  {
    if (tkey[t])
    {
      auto &group = info[tkey[t]];
      if (group.kind == BlockInfo::DENSE)
      {
        impl->dofs.push_back(t);
        slots.push_back(group.slot++);
        continue;
      }
      if (group.kind == BlockInfo::SHEET)
      {
        impl->sheets.dofs.push_back(t);
        continue;
      }
      tkey[t] = 0;
    }
  }
  if (census.block_dofs == 0 && census.sheets == 0)
  {
    return;
  }

  // 5. Groups of the L-DOFs (through the true DOFs they depend on), and the entries of the
  //    operator and of its definite counterpart within the groups.
  std::unique_ptr<mfem::HypreParMatrix> Af, Afd;
  {
    const std::vector<Key> okey = P.Halo(tkey);
    std::vector<int> lkeys_ptr(nl + 1, 0);
    std::vector<Key> lkeys, ks;
    for (int i = 0; i < nl; i++)
    {
      ks.clear();
      P.ForRow(i,
               [&](int c, bool offd, double)
               {
                 if (const Key k = offd ? okey[c] : tkey[c])
                 {
                   ks.push_back(k);
                 }
               });
      std::sort(ks.begin(), ks.end());
      ks.erase(std::unique(ks.begin(), ks.end()), ks.end());
      lkeys.insert(lkeys.end(), ks.begin(), ks.end());
      lkeys_ptr[i + 1] = static_cast<int>(lkeys.size());
    }
    Af = AssembleBlockEntries(A.LocalOperator(), fespace, lkeys_ptr, lkeys);
    if (A_def)
    {
      Afd = AssembleBlockEntries(A_def->LocalOperator(), fespace, lkeys_ptr, lkeys);
    }
  }

  // 6. Blocks: the matrices R A_f Rᵀ of the groups on their owners, where R (slots x true
  //    DOFs) gathers the DOFs of every group on its owner; the groups which are too large
  //    are split into pieces; the inverses of the blocks.
  if (census.block_dofs > 0)
  {
    const mfem::HypreParMatrix Rt =
        SelectionMatrix(fespace, impl->dofs, slots, slot_starts[rank],
                        slot_starts[rank + 1], slot_starts.back());
    std::unique_ptr<mfem::HypreParMatrix> R(Rt.Transpose());
    std::unique_ptr<mfem::HypreParMatrix> B(mfem::RAP(Af.get(), &Rt));
    std::unique_ptr<mfem::HypreParMatrix> Bd(Afd ? mfem::RAP(Afd.get(), &Rt) : nullptr);
    // The groups owned here, in slot order, and the local true DOF of every slot (-1 for
    // a DOF of another process).
    std::vector<std::pair<HYPRE_BigInt, int>> owned;
    for (const auto &[k, group] : info)
    {
      if (group.kind == BlockInfo::DENSE && group.owner == rank)
      {
        owned.emplace_back(group.block_slot - slot_starts[rank], group.size);
      }
    }
    std::sort(owned.begin(), owned.end());
    const int num_slots = static_cast<int>(slot_starts[rank + 1] - slot_starts[rank]);
    MFEM_VERIFY(hypre_CSRMatrixNumRows(DiagOf(B.get())) == num_slots,
                "Block Jacobi: slot count mismatch!");
    std::vector<int> piece(num_slots), slot_dof(num_slots, -1);
    {
      hypre_CSRMatrix *Rd = DiagOf(R.get());
      const HYPRE_Int *RI = hypre_CSRMatrixI(Rd), *RJ = hypre_CSRMatrixJ(Rd);
      for (int i = 0; i < num_slots; i++)
      {
        if (RI[i + 1] > RI[i])
        {
          slot_dof[i] = RJ[RI[i]];
        }
      }
    }
    int next = 0;
    bool any_split = false;
    for (const auto &[slot, size] : owned)
    {
      MFEM_VERIFY(slot == next, "Block Jacobi: unexpected slot of a group!");
      std::fill(piece.begin() + next, piece.begin() + next + size, next);
      any_split = any_split || (size > block_size);
      next += size;
    }
    MFEM_VERIFY(next == num_slots, "Block Jacobi: unexpected number of slots!");
    // The matrix which gives the blocks their magnitudes: the definite counterpart if
    // there is one, else the operator.
    hypre_CSRMatrix *Bm = DiagOf(Bd ? Bd.get() : B.get());
    {
      // The ties are those of the groups which are too large, which have only DOFs of
      // this process (see step 4).
      const SpMat M = GroupMatrix(Bm, piece);
      std::unique_ptr<DofTies> ties;
      if (any_split)
      {
        ties = std::make_unique<DofTies>(fes, aggregates);
      }
      auto Tie = [&](int i, int j)
      { return (*ties)(tldof[slot_dof[i]], tldof[slot_dof[j]]); };
      for (const auto &[slot, size] : owned)
      {
        ClassBlocks::Split(M, static_cast<int>(slot), static_cast<int>(slot) + size,
                           block_size, Tie, piece);
      }
    }
    {
      const SpMat M = GroupMatrix(Bm, piece);
      if (Bd)
      {
        const SpMat As = GroupMatrix(DiagOf(B.get()), piece);
        impl->blocks.Setup(M, &As, piece);
      }
      else
      {
        impl->blocks.Setup(M, nullptr, piece);
      }
    }
    // Without gathered blocks the rows of the blocks are the local true DOFs themselves,
    // and nothing is communicated.
    if (census.gathered > 0)
    {
      impl->R = std::move(R);
      impl->g.SetSize(num_slots);
      impl->y.SetSize(num_slots);
      impl->g.UseDevice(false);
      impl->y.UseDevice(false);
      impl->zt.SetSize(nt);
      impl->zt.UseDevice(false);
    }
    else
    {
      impl->blocks.Renumber(slot_dof);
    }
    HYPRE_BigInt counts[NUM_CLASSES + 3];
    for (int q = 0; q < NUM_CLASSES; q++)
    {
      counts[q] = impl->blocks.num[q];
    }
    counts[NUM_CLASSES] = impl->blocks.num_dropped;
    counts[NUM_CLASSES + 1] = impl->blocks.num_negative;
    counts[NUM_CLASSES + 2] = 0;
    for (int q = 0; q < NUM_CLASSES; q++)
    {
      const HYPRE_BigInt c = HYPRE_BigInt(1) << q;
      counts[NUM_CLASSES + 2] += impl->blocks.num[q] * (c * (c + 1) / 2);
    }
    census.max_entries = counts[NUM_CLASSES + 2];
    Mpi::GlobalMax(1, &census.max_entries, comm);
    Mpi::GlobalSum(NUM_CLASSES + 3, counts, comm);
    std::copy(counts, counts + NUM_CLASSES, census.blocks);
    census.dropped = counts[NUM_CLASSES];
    census.negative = counts[NUM_CLASSES + 1];
    census.entries = counts[NUM_CLASSES + 2];
  }

  // 7. Sheets: their matrix Sᵀ A_f S on their DOFs, which S selects, from the definite
  //    counterpart if there is one. Without a definite matrix, their DOFs are left to the
  //    point Jacobi method.
  auto &sheets = impl->sheets;
  sheets.comm = comm;
  sheets.num_global = census.sheet_dofs;
  if (sheets.num_global > 0)
  {
    HYPRE_BigInt end = sheets.Size();
    MPI_Scan(MPI_IN_PLACE, &end, 1, HYPRE_MPI_BIG_INT, MPI_SUM, comm);
    std::vector<HYPRE_BigInt> cols(sheets.Size());
    std::iota(cols.begin(), cols.end(), end - sheets.Size());
    const mfem::HypreParMatrix S = SelectionMatrix(
        fespace, sheets.dofs, cols, end - sheets.Size(), end, sheets.num_global);
    sheets.A.reset(mfem::RAP(Afd ? Afd.get() : Af.get(), &S));
    if (sheets.Setup(sheet_reduction))
    {
      census.sheet_entries = sheets.entries;
      census.sheet_degree = sheets.degree;
      census.sheet_limited = sheets.limited;
    }
    else
    {
      census.skipped += census.sheets;
      census.skipped_dofs += census.sheet_dofs;
      census.sheets = census.sheet_dofs = 0;
      sheets = SheetSolve();
    }
  }
}

const AggregateBlocks::Census &AggregateBlocks::GetCensus() const
{
  return impl->census;
}

const std::vector<int> &AggregateBlocks::GetBlockDofs() const
{
  return impl->dofs;
}

void AggregateBlocks::MultInverse(const Vector &x, Vector &z) const
{
  if (impl->dofs.empty() && !impl->R)
  {
    return;
  }
  if (!impl->R)
  {
    impl->blocks.Apply(x.HostRead(), z.HostReadWrite());
    return;
  }
  // Blocks with DOFs of other processes: gather on the owners, apply, and return.
  impl->R->Mult(x, impl->g);
  impl->blocks.Apply(impl->g.HostRead(), impl->y.HostWrite());
  impl->R->MultTranspose(impl->y, impl->zt);
  const auto *ZT = impl->zt.HostRead();
  auto *Z = z.HostReadWrite();
  for (int i : impl->dofs)
  {
    Z[i] = ZT[i];
  }
}

bool AggregateBlocks::HasSheets() const
{
  return impl->sheets.num_global > 0;
}

const std::vector<int> &AggregateBlocks::GetSheetDofs() const
{
  return impl->sheets.dofs;
}

void AggregateBlocks::AddSheetSolve(const Vector &r, Vector &y) const
{
  impl->sheets.AddMult(r.HostRead(), y.HostReadWrite());
}

}  // namespace palace
