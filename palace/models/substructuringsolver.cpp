// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolver.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/substructure.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "linalg/superlu.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Piecewise-constant (per element attribute) matrix coefficient for anisotropic materials.
// Attributes absent from the map contribute a zero tensor (e.g. environment attributes when
// assembling the region operator).
class PWMatrixCoefficient : public mfem::MatrixCoefficient
{
  const std::map<int, mfem::DenseMatrix> &mats;

public:
  PWMatrixCoefficient(int dim, const std::map<int, mfem::DenseMatrix> &m)
    : mfem::MatrixCoefficient(dim), mats(m)
  {
  }
  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    auto it = mats.find(T.Attribute);
    if (it != mats.end())
    {
      K = it->second;
    }
    else
    {
      K.SetSize(width);
      K = 0.0;
    }
  }
};

// Implicit environment Dirichlet-to-Neumann action on parent true DOFs:
//   y|_Gamma = A_GG x - A_GE A_EE^-1 A_EG x
// via A_env matvecs and one A_EE solve (A_EE = A_env with non-interior true DOFs
// eliminated).
class ImplicitDtN : public mfem::Operator
{
public:
  using ApplyFn = std::function<void(const mfem::Vector &, mfem::Vector &)>;
  ImplicitDtN(mfem::HypreParMatrix &A_env, ApplyFn Aee_inv,
              const std::vector<char> &is_gamma, const std::vector<char> &is_env_int)
    : mfem::Operator(A_env.Height()), A_env(A_env), Aee_inv(std::move(Aee_inv)),
      is_gamma(is_gamma), is_env_int(is_env_int), t(A_env.Height()), rhs(A_env.Height()),
      ye(A_env.Height()), t2(A_env.Height())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A_env.Mult(x, t);
    rhs = 0.0;
    for (int i = 0; i < height; i++)
    {
      if (is_env_int[i])
      {
        rhs(i) = t(i);
      }
    }
    ye = 0.0;
    Aee_inv(rhs, ye);
    for (int i = 0; i < height; i++)
    {
      if (!is_env_int[i])
      {
        ye(i) = 0.0;
      }
    }
    A_env.Mult(ye, t2);
    y = 0.0;
    for (int i = 0; i < height; i++)
    {
      if (is_gamma[i])
      {
        y(i) = t(i) - t2(i);
      }
    }
  }

private:
  mfem::HypreParMatrix &A_env;
  ApplyFn Aee_inv;
  const std::vector<char> &is_gamma, &is_env_int;
  mutable mfem::Vector t, rhs, ye, t2;
};

// Region-condensed system operator: (region operator on region-free true DOFs) + implicit
// DtN on the interface, with non-region-free true DOFs pinned to identity.
class RegionCondensedOperator : public mfem::Operator
{
public:
  RegionCondensedOperator(mfem::HypreParMatrix &A_region_free, const mfem::Operator &dtn,
                          const std::vector<char> &is_region_free)
    : mfem::Operator(A_region_free.Height()), A_region_free(A_region_free), dtn(dtn),
      is_region_free(is_region_free), td(A_region_free.Height())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A_region_free.Mult(x, y);
    dtn.Mult(x, td);
    for (int i = 0; i < height; i++)
    {
      if (is_region_free[i])
      {
        y(i) += td(i);
      }
    }
  }

private:
  mfem::HypreParMatrix &A_region_free;
  const mfem::Operator &dtn;
  const std::vector<char> &is_region_free;
  mutable mfem::Vector td;
};

// Preconditioner adapter: forwards Solver::Mult to a callable (used to route the region
// preconditioner through the region-submesh geometric multigrid). SetOperator is a no-op
// because the underlying multigrid is built once, up front.
class CallableSolver : public Solver<Operator>
{
public:
  CallableSolver(int h, std::function<void(const mfem::Vector &, mfem::Vector &)> f)
    : Solver<Operator>(false), apply(std::move(f))
  {
    height = width = h;
  }
  void SetOperator(const Operator &) override {}
  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { apply(x, y); }

private:
  std::function<void(const mfem::Vector &, mfem::Vector &)> apply;
};

// (MaterializedDtN is defined after the Hodlr helpers below, since it can apply either the
// dense per-rank row blocks or a compressed hierarchical operator.)

// Symmetric eigendecomposition by cyclic Jacobi (this mfem build has no LAPACK): A (n x n,
// row-major, destroyed), eigenvectors -> V (row-major, columns are eigenvectors),
// eigenvalues
// -> lam.
inline void SymEig(std::vector<double> &A, int n, std::vector<double> &V,
                   std::vector<double> &lam)
{
  V.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; i++)
  {
    V[static_cast<std::size_t>(i) * n + i] = 1.0;
  }
  auto a = [&](int i, int j) -> double & { return A[static_cast<std::size_t>(i) * n + j]; };
  auto v = [&](int i, int j) -> double & { return V[static_cast<std::size_t>(i) * n + j]; };
  for (int sweep = 0; sweep < 100; sweep++)
  {
    double off = 0.0;
    for (int p = 0; p < n; p++)
    {
      for (int q = p + 1; q < n; q++)
      {
        off += a(p, q) * a(p, q);
      }
    }
    if (off < 1e-30)
    {
      break;
    }
    for (int p = 0; p < n; p++)
    {
      for (int q = p + 1; q < n; q++)
      {
        const double apq = a(p, q);
        if (std::abs(apq) < 1e-300)
        {
          continue;
        }
        const double phi = 0.5 * std::atan2(2.0 * apq, a(q, q) - a(p, p));
        const double c = std::cos(phi), s = std::sin(phi);
        for (int k = 0; k < n; k++)
        {
          const double kp = a(k, p), kq = a(k, q);
          a(k, p) = c * kp - s * kq;
          a(k, q) = s * kp + c * kq;
        }
        for (int k = 0; k < n; k++)
        {
          const double pk = a(p, k), qk = a(q, k);
          a(p, k) = c * pk - s * qk;
          a(q, k) = s * pk + c * qk;
        }
        for (int k = 0; k < n; k++)
        {
          const double kp = v(k, p), kq = v(k, q);
          v(k, p) = c * kp - s * kq;
          v(k, q) = s * kp + c * kq;
        }
      }
    }
  }
  lam.assign(n, 0.0);
  for (int i = 0; i < n; i++)
  {
    lam[i] = a(i, i);
  }
}

// Hierarchical off-diagonal low-rank (HODLR) representation of the symmetric interface
// operator S_E. The interface DOFs are permuted (recursive coordinate-median clustering);
// diagonal leaf blocks are stored dense, and each internal node's off-diagonal coupling
// S[I1,I2] is stored as a low-rank factor pair U V^T (symmetry gives the transpose block V
// U^T). Storage and apply are O(nG * (leaf + rank * log nG)) instead of O(nG^2). Replicated
// across ranks (small once compressed); serializable for MPI broadcast and file I/O.
struct Hodlr
{
  struct Leaf
  {
    int s, m;               // start (permuted), size
    std::vector<double> D;  // m x m, row-major
  };
  struct Block
  {
    int s1, m1, s2, m2, r;  // row/col starts+sizes (permuted), rank
    std::vector<double> U;  // m1 x r, row-major
    std::vector<double> V;  // m2 x r, row-major (B[I1,I2] ~= U V^T)
  };
  int n = 0;
  std::vector<int> perm;  // permuted position -> interface global index
  std::vector<Leaf> leaves;
  std::vector<Block> blocks;

  long long Storage() const
  {
    long long s = 0;
    for (const auto &lf : leaves)
    {
      s += static_cast<long long>(lf.m) * lf.m;
    }
    for (const auto &b : blocks)
    {
      s += static_cast<long long>(b.r) * (b.m1 + b.m2);
    }
    return s;
  }

  // y = S x, both in the permuted ordering (length n).
  void Mult(const double *xp, double *yp) const
  {
    std::fill(yp, yp + n, 0.0);
    for (const auto &lf : leaves)
    {
      for (int i = 0; i < lf.m; i++)
      {
        const double *row = &lf.D[static_cast<std::size_t>(i) * lf.m];
        double acc = 0.0;
        for (int j = 0; j < lf.m; j++)
        {
          acc += row[j] * xp[lf.s + j];
        }
        yp[lf.s + i] += acc;
      }
    }
    std::vector<double> t;
    for (const auto &b : blocks)
    {
      t.assign(b.r, 0.0);  // t = V^T x2
      for (int j = 0; j < b.m2; j++)
      {
        const double *vj = &b.V[static_cast<std::size_t>(j) * b.r];
        const double xj = xp[b.s2 + j];
        for (int k = 0; k < b.r; k++)
        {
          t[k] += vj[k] * xj;
        }
      }
      for (int i = 0; i < b.m1; i++)  // y1 += U t
      {
        const double *ui = &b.U[static_cast<std::size_t>(i) * b.r];
        double acc = 0.0;
        for (int k = 0; k < b.r; k++)
        {
          acc += ui[k] * t[k];
        }
        yp[b.s1 + i] += acc;
      }
      t.assign(b.r, 0.0);  // t = U^T x1
      for (int i = 0; i < b.m1; i++)
      {
        const double *ui = &b.U[static_cast<std::size_t>(i) * b.r];
        const double xi = xp[b.s1 + i];
        for (int k = 0; k < b.r; k++)
        {
          t[k] += ui[k] * xi;
        }
      }
      for (int j = 0; j < b.m2; j++)  // y2 += V t
      {
        const double *vj = &b.V[static_cast<std::size_t>(j) * b.r];
        double acc = 0.0;
        for (int k = 0; k < b.r; k++)
        {
          acc += vj[k] * t[k];
        }
        yp[b.s2 + j] += acc;
      }
    }
  }

  // Flatten to a double buffer (indices packed as doubles) for MPI_Bcast and file I/O.
  std::vector<double> Serialize() const
  {
    std::vector<double> b;
    b.push_back(n);
    b.push_back(static_cast<double>(leaves.size()));
    b.push_back(static_cast<double>(blocks.size()));
    for (int p : perm)
    {
      b.push_back(p);
    }
    for (const auto &lf : leaves)
    {
      b.push_back(lf.s);
      b.push_back(lf.m);
      b.insert(b.end(), lf.D.begin(), lf.D.end());
    }
    for (const auto &bl : blocks)
    {
      b.push_back(bl.s1);
      b.push_back(bl.m1);
      b.push_back(bl.s2);
      b.push_back(bl.m2);
      b.push_back(bl.r);
      b.insert(b.end(), bl.U.begin(), bl.U.end());
      b.insert(b.end(), bl.V.begin(), bl.V.end());
    }
    return b;
  }

  static Hodlr Deserialize(const std::vector<double> &b)
  {
    Hodlr h;
    std::size_t p = 0;
    h.n = static_cast<int>(b[p++]);
    const int nleaf = static_cast<int>(b[p++]);
    const int nblk = static_cast<int>(b[p++]);
    h.perm.resize(h.n);
    for (int i = 0; i < h.n; i++)
    {
      h.perm[i] = static_cast<int>(b[p++]);
    }
    for (int i = 0; i < nleaf; i++)
    {
      Leaf lf;
      lf.s = static_cast<int>(b[p++]);
      lf.m = static_cast<int>(b[p++]);
      lf.D.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(lf.m) * lf.m);
      p += static_cast<std::size_t>(lf.m) * lf.m;
      h.leaves.push_back(std::move(lf));
    }
    for (int i = 0; i < nblk; i++)
    {
      Block bl;
      bl.s1 = static_cast<int>(b[p++]);
      bl.m1 = static_cast<int>(b[p++]);
      bl.s2 = static_cast<int>(b[p++]);
      bl.m2 = static_cast<int>(b[p++]);
      bl.r = static_cast<int>(b[p++]);
      bl.U.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(bl.m1) * bl.r);
      p += static_cast<std::size_t>(bl.m1) * bl.r;
      bl.V.assign(b.begin() + p, b.begin() + p + static_cast<std::size_t>(bl.m2) * bl.r);
      p += static_cast<std::size_t>(bl.m2) * bl.r;
      h.blocks.push_back(std::move(bl));
    }
    return h;
  }
};

// Recursively build a Hodlr from a dense symmetric S (n x n, row-major): split the index
// set `idx` at the median of its widest coordinate axis, compress the off-diagonal block
// via a truncated SVD to relative tolerance `tol`, recurse on the diagonal blocks. Leaf
// blocks (size <= `leaf`) are stored dense. `base` is the permuted offset of this block.
// `coords` is n x 3 (global-index order).
inline void BuildHodlr(const std::vector<double> &S, int n, const std::vector<int> &idx,
                       int base, const std::vector<double> &coords, double tol, int leaf,
                       Hodlr &h)
{
  const int m = static_cast<int>(idx.size());
  if (m <= leaf)
  {
    Hodlr::Leaf lf;
    lf.s = base;
    lf.m = m;
    lf.D.assign(static_cast<std::size_t>(m) * m, 0.0);
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
      {
        lf.D[static_cast<std::size_t>(i) * m + j] =
            S[static_cast<std::size_t>(idx[i]) * n + idx[j]];
      }
      h.perm[base + i] = idx[i];
    }
    h.leaves.push_back(std::move(lf));
    return;
  }
  int axis = 0;
  double best = -1.0;
  for (int d = 0; d < 3; d++)
  {
    double lo = 1e300, hi = -1e300;
    for (int i : idx)
    {
      const double x = coords[static_cast<std::size_t>(i) * 3 + d];
      lo = std::min(lo, x);
      hi = std::max(hi, x);
    }
    if (hi - lo > best)
    {
      best = hi - lo;
      axis = d;
    }
  }
  std::vector<int> s(idx);
  std::sort(s.begin(), s.end(),
            [&](int i, int j)
            {
              return coords[static_cast<std::size_t>(i) * 3 + axis] <
                     coords[static_cast<std::size_t>(j) * 3 + axis];
            });
  std::vector<int> I1(s.begin(), s.begin() + m / 2), I2(s.begin() + m / 2, s.end());
  const int m1 = static_cast<int>(I1.size()), m2 = static_cast<int>(I2.size());
  // Recurse first so the final permuted order within each child range is fixed; the
  // off-diagonal factors below are then indexed consistently with `perm` (the recursion
  // reorders I1/I2 internally).
  BuildHodlr(S, n, I1, base, coords, tol, leaf, h);
  BuildHodlr(S, n, I2, base + m1, coords, tol, leaf, h);
  // Off-diagonal block B = S[rows, cols] with rows/cols the final permuted indices.
  const int *rows = &h.perm[base], *cols = &h.perm[base + m1];
  // C = B^T B (m2 x m2); singular values of B are sqrt(eig(C)).
  std::vector<double> C(static_cast<std::size_t>(m2) * m2, 0.0);
  for (int a = 0; a < m2; a++)
  {
    for (int b = a; b < m2; b++)
    {
      double acc = 0.0;
      for (int k = 0; k < m1; k++)
      {
        acc += S[static_cast<std::size_t>(rows[k]) * n + cols[a]] *
               S[static_cast<std::size_t>(rows[k]) * n + cols[b]];
      }
      C[static_cast<std::size_t>(a) * m2 + b] = acc;
      C[static_cast<std::size_t>(b) * m2 + a] = acc;
    }
  }
  std::vector<double> Vc, lam;
  SymEig(C, m2, Vc, lam);
  double lmax = 0.0;
  for (double l : lam)
  {
    lmax = std::max(lmax, l);
  }
  std::vector<int> keep;
  for (int k = 0; k < m2; k++)
  {
    if (lam[k] > tol * tol * lmax)
    {
      keep.push_back(k);
    }
  }
  const int rB = static_cast<int>(keep.size());
  // B ~= (B V_keep) V_keep^T; store U := B V_keep (m1 x rB), V := V_keep (m2 x rB).
  Hodlr::Block bl;
  bl.s1 = base;
  bl.m1 = m1;
  bl.s2 = base + m1;
  bl.m2 = m2;
  bl.r = rB;
  bl.U.assign(static_cast<std::size_t>(m1) * rB, 0.0);
  bl.V.assign(static_cast<std::size_t>(m2) * rB, 0.0);
  for (int i = 0; i < m1; i++)
  {
    for (int kk = 0; kk < rB; kk++)
    {
      double acc = 0.0;
      for (int l = 0; l < m2; l++)
      {
        acc += S[static_cast<std::size_t>(rows[i]) * n + cols[l]] *
               Vc[static_cast<std::size_t>(l) * m2 + keep[kk]];
      }
      bl.U[static_cast<std::size_t>(i) * rB + kk] = acc;
    }
  }
  for (int j = 0; j < m2; j++)
  {
    for (int kk = 0; kk < rB; kk++)
    {
      bl.V[static_cast<std::size_t>(j) * rB + kk] =
          Vc[static_cast<std::size_t>(j) * m2 + keep[kk]];
    }
  }
  h.blocks.push_back(std::move(bl));
}

// Materialized DtN: applies the precomputed interface operator S_E to a distributed
// interface vector. The global interface vector is gathered by an Allreduce, then either
// each rank multiplies its own dense row block (S_rows) or -- when a compressed Hodlr is
// supplied -- the replicated hierarchical operator is applied and each rank reads back its
// interface rows.
class MaterializedDtN : public mfem::Operator
{
public:
  MaterializedDtN(const std::vector<double> &S_rows, int row_off,
                  const std::vector<int> &gamma_global, int nG_global, MPI_Comm comm,
                  const Hodlr *hodlr = nullptr)
    : mfem::Operator(static_cast<int>(gamma_global.size())), S_rows(S_rows),
      row_off(row_off), gamma_global(gamma_global), nG_global(nG_global), comm(comm),
      hodlr(hodlr)
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    std::vector<double> xl(nG_global, 0.0), xg(nG_global, 0.0);
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        xl[gamma_global[i]] = x(i);
      }
    }
    MPI_Allreduce(xl.data(), xg.data(), nG_global, MPI_DOUBLE, MPI_SUM, comm);
    y = 0.0;
    if (hodlr)
    {
      // Compressed apply: permute to tree order, apply the replicated Hodlr, un-permute,
      // and read back this rank's interface rows.
      std::vector<double> xp(nG_global), yp(nG_global), yg(nG_global);
      for (int k = 0; k < nG_global; k++)
      {
        xp[k] = xg[hodlr->perm[k]];
      }
      hodlr->Mult(xp.data(), yp.data());
      for (int k = 0; k < nG_global; k++)
      {
        yg[hodlr->perm[k]] = yp[k];
      }
      for (int i = 0; i < height; i++)
      {
        if (gamma_global[i] >= 0)
        {
          y(i) = yg[gamma_global[i]];
        }
      }
      return;
    }
    for (int i = 0; i < height; i++)
    {
      if (gamma_global[i] >= 0)
      {
        const int r = gamma_global[i] - row_off;  // owned row block is contiguous
        const double *row = &S_rows[static_cast<std::size_t>(r) * nG_global];
        double s = 0.0;
        for (int j = 0; j < nG_global; j++)
        {
          s += row[j] * xg[j];
        }
        y(i) = s;
      }
    }
  }

private:
  const std::vector<double> &S_rows;
  int row_off;
  const std::vector<int> &gamma_global;
  int nG_global;
  MPI_Comm comm;
  const Hodlr *hodlr;
};

}  // namespace

// Parallel region-condensed static solve (electrostatic H1 or magnetostatic H(curl)):
// region/environment operators assembled on the parent space with domain-restricted
// material coefficients, interface identified in true-DOF space, environment condensed via
// a distributed DtN. Magnetostatics uses a small curl-curl mass regularization (see
// kMagRegularization).
struct SubstructuringSolver::Impl
{
  const IoData &iodata;
  mfem::ParMesh &parent;
  bool magnetostatic;
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  mfem::ParFiniteElementSpace parent_fes;
  int nt;

  std::unique_ptr<mfem::HypreParMatrix> A_region, A_env, A_region_free;
  // Environment interior solve on the environment submesh (Gamma as an essential boundary):
  // avoids the all-identity ranks of the parent-space eliminated operator and enables
  // geometric multigrid on a standard Dirichlet problem.
  mfem::Array<int> ra_arr, ea_arr;
  std::unique_ptr<mfem::ParSubMesh> env_submesh_owned;    // single-level ownership
  std::vector<std::unique_ptr<Mesh>> env_mesh_vec;        // GMG: [0] owns the ParSubMesh
  mfem::ParSubMesh *env_submesh = nullptr;                // raw ptr to the owned submesh
  std::unique_ptr<mfem::ParFiniteElementSpace> env_sfes;  // single-level solve space
  mfem::ParFiniteElementSpace *env_solve_fes = nullptr;   // fespace the solver acts on
  std::vector<std::unique_ptr<mfem::H1_FECollection>> env_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> env_hierarchy;
  std::unique_ptr<MultigridOperator> env_mg_op;
  std::vector<mfem::Array<int>>
      env_dbc_lists;  // per-level essential (ParOperator MakeRefs)
  std::unique_ptr<mfem::HypreParMatrix> env_A_ee;  // single-level submesh operator
  std::unique_ptr<KspSolver> env_ksp;
#if defined(MFEM_USE_SUPERLU)
  std::unique_ptr<SuperLUSolver> reg_lu;  // direct A_region_free factorization (region pc)
  std::unique_ptr<mfem::HypreParMatrix>
      env_A_ee_parent;  // parent-space env A_EE (eliminated)
  std::unique_ptr<SuperLUSolver>
      env_lu_parent;  // parent-space direct env solve (no transfer)
#endif
  // Block size for the batched multi-RHS S_E materialization (SuperLU triangular solves are
  // per-call-overhead bound for a single RHS; ~4x cheaper per column at 32).
  static constexpr int kMaterializeBlock = 32;
  // Largest interface for which the region preconditioner factors the dense Gamma block of
  // the exact condensed operator (dense block + its fill ~ 3 |Gamma|^2 doubles in total).
  static constexpr int kDirectCondensedMaxInterface = 10000;
  // RHS count of the direct env factor, fixed by its first solve (MFEM's SuperLU wrapper
  // does not allow changing it): kMaterializeBlock when S_E is materialized, else the size
  // of the first excitation batch. Smaller batches are zero-padded to it.
  mutable int env_block = 0;
  mutable mfem::Vector env_pad_in,
      env_pad_out;               // shared zero input / scratch output (padding)
  mfem::Array<int> non_env_int;  // parent true DOFs that are not environment-interior
  bool env_parent_direct = false;
  mfem::Array<int> env_ess;  // solve-space essential true DOFs (Gamma + env Dirichlet)
  mutable mfem::ParGridFunction env_pgf, env_sgf;  // parent / submesh transfer buffers
  mutable mfem::Vector env_srhs, env_ssol;
  // Region-solve geometric multigrid (Phase C): the region preconditioner runs on the
  // region ParSubMesh (region Dirichlet terminals essential, interface Gamma free),
  // avoiding the all-identity ranks of the parent-space region-free operator.
  std::vector<std::unique_ptr<Mesh>> reg_mesh_vec;
  mfem::ParSubMesh *reg_submesh = nullptr;
  mfem::ParFiniteElementSpace *reg_solve_fes = nullptr;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> reg_fecs;
  std::unique_ptr<FiniteElementSpaceHierarchy> reg_hierarchy;
  std::unique_ptr<MultigridOperator> reg_mg_op;
  std::vector<mfem::Array<int>> reg_dbc_lists;
  std::unique_ptr<KspSolver> reg_gmg_ksp;  // inner GMG solve on the region submesh
  mfem::Array<int> reg_ess;
  bool region_gmg = false;
  mutable mfem::ParGridFunction reg_pgf, reg_sgf;
  mutable mfem::Vector reg_srhs, reg_ssol;
  // Energy (QoI) operators. Electrostatic: alias the solve operators. Magnetostatic: pure
  // curl-curl (no mass), so the magnetic energy / inductance is physical.
  std::unique_ptr<mfem::HypreParMatrix> A_region_energy, A_env_energy;
  mfem::HypreParMatrix *K_region_e = nullptr, *K_env_e = nullptr;
  std::vector<char> is_gamma, is_env_int, is_region_free;
  mfem::Array<int> dbc_tdofs;
  std::map<int, std::vector<int>> terminal_tdofs;  // terminal index -> its true DOFs

  std::unique_ptr<RegionCondensedOperator> region_op;
  std::unique_ptr<KspSolver> region_ksp;  // Palace CG + wrapped AMG/AMS on the region block
  std::unique_ptr<ImplicitDtN> dtn;

  // Materialized (reusable) interface operator: replicated dense S_E + load g_E over a
  // global interface enumeration, computed once so region solves need no environment
  // solves.
  std::vector<int> gamma_global;  // owned parent true DOF -> global interface index, or -1
  int nG_global = 0;
  // Distributed dense interface operator S_E: each rank stores the contiguous block of rows
  // [gamma_off, gamma_off + gamma_nloc) it owns (row-major, gamma_nloc x nG_global),
  // instead of the full nG x nG matrix replicated on every rank. Cuts per-rank storage
  // O(nG^2) -> O(nG^2 / P) and parallelizes the apply.
  std::vector<double> S_rows;
  int gamma_off = 0, gamma_nloc = 0;
  // Optional compressed (HODLR) form of S_E, used in place of S_rows when interface
  // compression is requested (replicated across ranks).
  std::unique_ptr<Hodlr> hodlr;
  std::unique_ptr<MaterializedDtN> mat_dtn;

  // Electrostatic terminal modes (environment condensed onto Gamma + terminal potentials):
  // for each terminal k with unit mode x_k (1 on its Dirichlet DOFs), g_k = (A_env x_k -
  // A_env A_EE^-1 (A_env x_k)|_E)|_Gamma (the interface load, stored like S_E by owned
  // interface rows: G_rows[r * K + k]) and the replicated K x K environment energy block
  // Cmode[k][l] = x_k^T A_env x_l - (A_env x_k)|_E^T A_EE^-1 (A_env x_l)|_E. Together with
  // S_E they give the environment's exact energy contribution without an environment solve.
  std::vector<int> mode_ids;  // terminal index per mode column (terminal_tdofs key order)
  std::vector<double> G_rows, Cmode;
  bool modes_ready = false;
  bool env_built = false;  // environment interior operator factored / set up

  Impl(const IoData &iodata, mfem::ParMesh &parent)
    : iodata(iodata), parent(parent),
      magnetostatic(iodata.problem.type == ProblemType::MAGNETOSTATIC),
      fec(magnetostatic
              ? std::unique_ptr<mfem::FiniteElementCollection>(
                    new mfem::ND_FECollection(iodata.solver.order, parent.Dimension()))
              : std::unique_ptr<mfem::FiniteElementCollection>(
                    new mfem::H1_FECollection(iodata.solver.order, parent.Dimension()))),
      parent_fes(&parent, fec.get()), nt(parent_fes.GetTrueVSize())
  {
    const auto &sub = *iodata.solver.substructuring;
    mfem::Array<int> ra(static_cast<int>(sub.region_attributes.size())),
        ea(static_cast<int>(sub.environment_attributes.size()));
    std::copy(sub.region_attributes.begin(), sub.region_attributes.end(), ra.begin());
    std::copy(sub.environment_attributes.begin(), sub.environment_attributes.end(),
              ea.begin());
    ra_arr = ra;
    ea_arr = ea;

    // Interface / region / environment true-DOF markers.
    mfem::Array<int> rm, em, im;
    MarkInterfaceTrueDofs(parent_fes, ra, ea, rm, em, im);

    // Terminal Dirichlet true DOFs (per terminal). The Dirichlet DOF *set* is the union of
    // all terminals and is fixed across excitations; only the prescribed values change, so
    // the interface/interior partition below is excitation-independent.
    const auto &terminals = iodata.boundaries.terminal;
    mfem::Array<int> dir_mark(nt);
    dir_mark = 0;
    {
      const int maxb = parent.bdr_attributes.Size() ? parent.bdr_attributes.Max() : 0;
      for (const auto &[idx, term] : terminals)
      {
        mfem::Array<int> ess_bdr(maxb), ess;
        ess_bdr = 0;
        for (int a : term.attributes)
        {
          if (a >= 1 && a <= maxb)
          {
            ess_bdr[a - 1] = 1;
          }
        }
        parent_fes.GetEssentialTrueDofs(ess_bdr, ess);
        auto &list = terminal_tdofs[idx];
        for (int i = 0; i < ess.Size(); i++)
        {
          dir_mark[ess[i]] = 1;
          list.push_back(ess[i]);
        }
      }
    }
    for (int i = 0; i < nt; i++)
    {
      if (dir_mark[i])
      {
        dbc_tdofs.Append(i);
      }
    }

    is_gamma.assign(nt, 0);
    is_env_int.assign(nt, 0);
    is_region_free.assign(nt, 0);
    for (int i = 0; i < nt; i++)
    {
      if (dir_mark[i])
      {
        continue;
      }
      const bool r = rm[i], e = em[i];
      if (r && e)
      {
        is_gamma[i] = 1;
      }
      else if (e)
      {
        is_env_int[i] = 1;
      }
      if (r)
      {
        is_region_free[i] = 1;  // region-free includes the interface
      }
    }

    // Domain-restricted scalar permittivity coefficients (zero outside the subdomain).
    const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
    const int dim = parent.Dimension();
    auto in = [](const mfem::Array<int> &s, int a)
    {
      for (int x : s)
      {
        if (x == a)
        {
          return true;
        }
      }
      return false;
    };
    // Reconstruct the (possibly anisotropic) material tensor per attribute from its eigen-
    // decomposition: M_ij = sum_k s[k] v[k]_i v[k]_j. For magnetostatics the curl-curl
    // coefficient is the inverse permeability, so the reconstructed mu tensor is inverted.
    auto tensor = [dim](const config::SymmetricMatrixData<3> &prop, bool invert)
    {
      mfem::DenseMatrix e(dim);
      e = 0.0;
      for (int k = 0; k < 3; k++)
      {
        for (int i = 0; i < dim; i++)
        {
          for (int j = 0; j < dim; j++)
          {
            e(i, j) += prop.s[k] * prop.v[k][i] * prop.v[k][j];
          }
        }
      }
      if (invert)
      {
        e.Invert();
      }
      return e;
    };
    for (const auto &mat : iodata.domains.materials)
    {
      for (int a : mat.attributes)
      {
        if (a < 1 || a > max_attr)
        {
          continue;
        }
        if (in(ra, a))
        {
          region_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
        }
        else if (in(ea, a))
        {
          env_eps[a] = tensor(magnetostatic ? mat.mu_r : mat.epsilon_r, magnetostatic);
        }
      }
    }
    A_region = AssembleParent(region_eps, magnetostatic);
    A_env = AssembleParent(env_eps, magnetostatic);
    if (magnetostatic)
    {
      A_region_energy = AssembleParent(region_eps, false);  // pure curl-curl
      A_env_energy = AssembleParent(env_eps, false);
      K_region_e = A_region_energy.get();
      K_env_e = A_env_energy.get();
    }
    else
    {
      K_region_e = A_region.get();
      K_env_e = A_env.get();
    }
  }

  // Small mass regularization making the singular magnetostatic curl-curl definite (AMS-
  // convergent). Exact for divergence-free-compatible excitations (flux loops); a
  // gauge-free (pseudo-inverse) treatment for general surface currents is a follow-up.
  static constexpr double kMagRegularization = 1.0e-3;

  // Environment/region size (global true DOFs) below which a sparse-direct factorization is
  // used by default for the interior solve; above it, fall back to iterative / GMG. Chosen
  // so typical region-in-chip environments factor comfortably while very large domains do
  // not exhaust memory on the SuperLU factorization.
  static constexpr long long kDirectMaxDofs = 2000000;

  std::unique_ptr<mfem::HypreParMatrix>
  AssembleParent(const std::map<int, mfem::DenseMatrix> &coef_by_attr, bool with_mass)
  {
    PWMatrixCoefficient coef(parent.Dimension(), coef_by_attr);
    mfem::ParBilinearForm a(&parent_fes);
    if (magnetostatic)
    {
      a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(coef));
      const int max_attr = parent.attributes.Size() ? parent.attributes.Max() : 1;
      mfem::Vector mass(max_attr);
      mass = 0.0;
      if (with_mass)
      {
        for (const auto &[a_attr, m] : coef_by_attr)
        {
          mass(a_attr - 1) = kMagRegularization;
        }
      }
      mfem::PWConstCoefficient mcoef(mass);  // outlives Assemble() below
      if (with_mass)
      {
        a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mcoef));
      }
      a.Assemble();
      a.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    }
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  }

  std::map<int, mfem::DenseMatrix> region_eps, env_eps;

  // Environment interior solver: A_EE^-1 (env-interior with Gamma + env Dirichlet held
  // fixed), built on the environment submesh so every rank owns real DOFs (no all-identity
  // ranks).
  // Build (factor) the environment interior solver on first use: an online run that loads
  // S_E and the terminal modes never needs it unless environment fields are requested.
  void EnsureEnv()
  {
    if (!env_built)
    {
      BuildEnvSubmeshSolver();
      env_built = true;
    }
  }

  // Unit Dirichlet mode of terminal idx: 1 on its true DOFs.
  mfem::Vector TerminalMode(int idx) const
  {
    mfem::Vector x(nt);
    x = 0.0;
    for (int d : terminal_tdofs.at(idx))
    {
      x(d) = 1.0;
    }
    return x;
  }

  // Terminal-mode couplings g_k and Cmode (see the members) from K batched environment
  // solves; terminals whose mode does not reach the environment interior skip the solve.
  void ComputeTerminalModes()
  {
    EnsureEnv();
    MPI_Comm comm = parent_fes.GetComm();
    mode_ids.clear();
    for (const auto &[idx, dofs] : terminal_tdofs)
    {
      mode_ids.push_back(idx);
    }
    const int K = static_cast<int>(mode_ids.size());
    std::vector<mfem::Vector> x(K), t(K, mfem::Vector(nt)), w(K, mfem::Vector(nt));
    std::vector<int> touches(K, 0);
    for (int k = 0; k < K; k++)
    {
      x[k] = TerminalMode(mode_ids[k]);
      A_env->Mult(x[k], t[k]);
      for (int i = 0; i < nt && !touches[k]; i++)
      {
        touches[k] = (is_env_int[i] && t[k](i) != 0.0);
      }
    }
    if (K > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, touches.data(), K, MPI_INT, MPI_MAX, comm);
    }
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    for (int k = 0; k < K; k++)
    {
      w[k] = 0.0;
      if (touches[k])
      {
        X.push_back(&t[k]);
        Y.push_back(&w[k]);
      }
    }
    if (!X.empty())
    {
      ApplyAeeInvMulti(X, Y);
    }
    G_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
    mfem::Vector t2(nt);
    for (int k = 0; k < K; k++)
    {
      A_env->Mult(w[k], t2);
      for (int i = 0; i < nt; i++)
      {
        if (is_gamma[i])
        {
          G_rows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * K + k] =
              t[k](i) - t2(i);
        }
      }
    }
    Cmode.assign(static_cast<std::size_t>(K) * K, 0.0);
    for (int k = 0; k < K; k++)
    {
      for (int l = 0; l < K; l++)
      {
        Cmode[static_cast<std::size_t>(k) * K + l] = (x[k] * t[l]) - (t[k] * w[l]);
      }
    }
    if (K > 0)
    {
      MPI_Allreduce(MPI_IN_PLACE, Cmode.data(), K * K, MPI_DOUBLE, MPI_SUM, comm);
    }
    modes_ready = true;
  }

  // Recover the environment interior of full fields (region + Dirichlet parts set),
  // batched: u_E = -A_EE^-1 (A_env u)|_E.
  void RecoverEnvInterior(std::vector<mfem::Vector> &u)
  {
    EnsureEnv();
    const int n = static_cast<int>(u.size());
    std::vector<mfem::Vector> rhs(n, mfem::Vector(nt));
    std::vector<const mfem::Vector *> X(n);
    std::vector<mfem::Vector *> Y(n);
    for (int k = 0; k < n; k++)
    {
      A_env->Mult(u[k], rhs[k]);
      rhs[k].Neg();
      X[k] = &rhs[k];
      Y[k] = &rhs[k];  // solve in place
    }
    if (n > 0)
    {
      ApplyAeeInvMulti(X, Y);
    }
    for (int k = 0; k < n; k++)
    {
      for (int i = 0; i < nt; i++)
      {
        if (is_env_int[i])
        {
          u[k](i) = rhs[k](i);
        }
      }
    }
  }

  // Per-rank counts/displacements of the interface-row-distributed G (gamma_nloc * K each).
  void ModeRowLayout(int K, std::vector<int> &cnt, std::vector<int> &disp) const
  {
    MPI_Comm comm = parent_fes.GetComm();
    const int nranks = Mpi::Size(comm);
    cnt.assign(nranks, 0);
    disp.assign(nranks, 0);
    const int mine = gamma_nloc * K;
    MPI_Allgather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nranks; r++)
    {
      disp[r] = disp[r - 1] + cnt[r - 1];
    }
  }

  // Append the terminal modes to a saved model (collective; rank 0 writes).
  static constexpr int kModesMagic = 0x53444f4d;  // "MODS"
  void AppendModes(const std::string &path) const
  {
    MPI_Comm comm = parent_fes.GetComm();
    const int K = static_cast<int>(mode_ids.size());
    std::vector<int> cnt, disp;
    ModeRowLayout(K, cnt, disp);
    const bool root = (Mpi::Rank(comm) == 0);
    std::vector<double> G_full;
    if (root)
    {
      G_full.assign(static_cast<std::size_t>(nG_global) * K, 0.0);
    }
    MPI_Gatherv(G_rows.data(), gamma_nloc * K, MPI_DOUBLE, root ? G_full.data() : nullptr,
                cnt.data(), disp.data(), MPI_DOUBLE, 0, comm);
    if (root)
    {
      std::ofstream f(path, std::ios::binary | std::ios::app);
      f.write(reinterpret_cast<const char *>(&kModesMagic), sizeof(int));
      f.write(reinterpret_cast<const char *>(&K), sizeof(int));
      f.write(reinterpret_cast<const char *>(mode_ids.data()), sizeof(int) * K);
      f.write(reinterpret_cast<const char *>(G_full.data()),
              sizeof(double) * G_full.size());
      f.write(reinterpret_cast<const char *>(Cmode.data()), sizeof(double) * Cmode.size());
    }
  }

  // Load the terminal modes following S_E at byte offset `pos` (collective). Interface rows
  // are re-ordered onto the online numbering like S_E: row g <- saved row perm[g], times
  // sgn[g]. Leaves modes_ready false (so they are recomputed on demand) if the file
  // predates modes or the terminal set differs.
  void LoadModes(const std::string &path, std::streamoff pos, const std::vector<int> &perm,
                 const std::vector<double> &sgn)
  {
    MPI_Comm comm = parent_fes.GetComm();
    const bool root = (Mpi::Rank(comm) == 0);
    int hdr[2] = {0, 0};
    std::ifstream f;
    if (root)
    {
      f.open(path, std::ios::binary);
      f.seekg(pos);
      if (!f.read(reinterpret_cast<char *>(hdr), 2 * sizeof(int)))
      {
        hdr[0] = hdr[1] = 0;
      }
    }
    MPI_Bcast(hdr, 2, MPI_INT, 0, comm);
    const int K = hdr[1];
    if (hdr[0] != kModesMagic || K < 0)
    {
      return;
    }
    std::vector<int> ids(K);
    std::vector<double> G_full, C(static_cast<std::size_t>(K) * K);
    if (root)
    {
      std::vector<double> G_off(static_cast<std::size_t>(nG_global) * K);
      f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * K);
      f.read(reinterpret_cast<char *>(G_off.data()), sizeof(double) * G_off.size());
      f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
      G_full.assign(G_off.size(), 0.0);
      for (int g = 0; g < nG_global; g++)
      {
        for (int k = 0; k < K; k++)
        {
          G_full[static_cast<std::size_t>(g) * K + k] =
              sgn[g] * G_off[static_cast<std::size_t>(perm[g]) * K + k];
        }
      }
    }
    MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
    MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
    std::vector<int> cur;
    for (const auto &[idx, dofs] : terminal_tdofs)
    {
      cur.push_back(idx);
    }
    if (ids != cur)
    {
      return;  // different terminal set: recompute on demand
    }
    std::vector<int> cnt, disp;
    ModeRowLayout(K, cnt, disp);
    G_rows.assign(static_cast<std::size_t>(gamma_nloc) * K, 0.0);
    MPI_Scatterv(root ? G_full.data() : nullptr, cnt.data(), disp.data(), MPI_DOUBLE,
                 G_rows.data(), gamma_nloc * K, MPI_DOUBLE, 0, comm);
    mode_ids = ids;
    Cmode = C;
    modes_ready = true;
  }

  void BuildEnvSubmeshSolver()
  {
    const int mg_levels = iodata.solver.linear.mg_max_levels;
#if defined(MFEM_USE_SUPERLU)
    // Direct A_EE (default when the environment fits): factor A_env with the
    // non-env-interior true DOFs eliminated, in PARENT space. SuperLU tolerates the
    // all-identity ranks that rule out a parent-space *iterative* solve, and no ParSubMesh
    // transfer is needed per solve (measured ~40% of the S_E materialization). Order >= 2
    // H(curl) always takes this path: the transfer mishandles higher-order tetrahedral
    // edge/face orientation across a cut.
    bool direct_fits = false;
    {
      long long env_loc = 0;
      for (int i = 0; i < nt; i++)
      {
        if (is_env_int[i] || is_gamma[i])
        {
          env_loc++;
        }
      }
      long long env_glob = 0;
      MPI_Allreduce(&env_loc, &env_glob, 1, MPI_LONG_LONG, MPI_SUM, parent_fes.GetComm());
      direct_fits = (env_glob <= kDirectMaxDofs);
    }
    if (direct_fits || (magnetostatic && iodata.solver.order > 1))
    {
      non_env_int.SetSize(0);
      for (int i = 0; i < nt; i++)
      {
        if (!is_env_int[i])
        {
          non_env_int.Append(i);
        }
      }
      env_A_ee_parent = std::make_unique<mfem::HypreParMatrix>(*A_env);
      {
        std::unique_ptr<mfem::HypreParMatrix> e(
            env_A_ee_parent->EliminateRowsCols(non_env_int));
      }
      env_lu_parent = std::make_unique<SuperLUSolver>(iodata, parent_fes.GetComm(), 0);
      env_lu_parent->SetOperator(*env_A_ee_parent);
      env_parent_direct = true;
      return;
    }
#endif
    // Iterative fallback (very large environment, or no SuperLU): environment submesh
    // solve.
    const bool use_gmg = !magnetostatic && iodata.solver.order > 1 && mg_levels > 1;

    // Create the environment submesh. For the GMG path it must be owned by a Palace Mesh
    // (for the CEED attribute maps + the FE-space hierarchy); otherwise a standalone
    // ParSubMesh suffices for the MFEM assembly + transfer.
    if (use_gmg)
    {
      env_mesh_vec.clear();
      env_mesh_vec.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
          mfem::ParSubMesh::CreateFromDomain(parent, ea_arr))));
      env_submesh = dynamic_cast<mfem::ParSubMesh *>(&env_mesh_vec[0]->Get());
      env_mesh_vec[0]->RebuildCeedAttributes();
    }
    else
    {
      env_submesh_owned = std::make_unique<mfem::ParSubMesh>(
          mfem::ParSubMesh::CreateFromDomain(parent, ea_arr));
      env_submesh = env_submesh_owned.get();
    }
    env_pgf.SetSpace(&parent_fes);

    if (!(use_gmg && BuildEnvGmg()))
    {
      // Single-level submesh Dirichlet solve (order 1, magnetostatic H(curl), or GMG
      // fallback): standalone FE space + wrapped AMG / AMS.
      env_sfes = std::make_unique<mfem::ParFiniteElementSpace>(env_submesh, fec.get());
      env_solve_fes = env_sfes.get();
      env_sgf.SetSpace(env_solve_fes);
      ComputeEnvEss();
      env_A_ee = AssembleEnvSubmesh();
      {
        std::unique_ptr<mfem::HypreParMatrix> e(env_A_ee->EliminateRowsCols(env_ess));
      }
      MPI_Comm comm = env_solve_fes->GetComm();
      {
        std::unique_ptr<Solver<Operator>> pc;
        if (magnetostatic)
        {
          auto ams = std::make_unique<mfem::HypreAMS>(env_solve_fes);
          ams->SetPrintLevel(0);
          pc = std::make_unique<MfemWrapperSolver<Operator>>(std::move(ams), true, false,
                                                             false);
        }
        else
        {
          pc = std::make_unique<MfemWrapperSolver<Operator>>(
              std::make_unique<BoomerAmgSolver>(1, 1, true, 0), true, false, false);
        }
        auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
        pcg->SetInitialGuess(false);
        pcg->SetRelTol(1.0e-12);
        pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
        pcg->SetMaxIter(1000);
        env_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
        env_ksp->SetOperators(*env_A_ee, *env_A_ee);
      }
    }

    env_srhs.SetSize(env_solve_fes->GetTrueVSize());
    env_ssol.SetSize(env_solve_fes->GetTrueVSize());
  }

  // Essential submesh DOFs = the parent non-(environment-interior) DOFs (interface Gamma +
  // environment Dirichlet), mapped onto the solve space by transferring the parent marker.
  // Robust to Palace inserting material-interface boundary elements at Gamma.
  void ComputeEnvEss()
  {
    mfem::Vector t(nt);
    for (int i = 0; i < nt; i++)
    {
      t(i) = is_env_int[i] ? 0.0 : 1.0;
    }
    mfem::ParGridFunction pg(&parent_fes), sg(env_solve_fes);
    pg.SetFromTrueDofs(t);
    sg = 0.0;
    env_submesh->Transfer(pg, sg);
    mfem::Vector st(env_solve_fes->GetTrueVSize());
    sg.GetTrueDofs(st);
    env_ess.SetSize(0);
    for (int i = 0; i < st.Size(); i++)
    {
      if (std::abs(st(i)) > 0.5)  // fabs: ND transfer may flip the marker's sign
      {
        env_ess.Append(i);
      }
    }
  }

  // Higher-order H1 environment Dirichlet solve via geometric p-multigrid on the submesh.
  // Returns false (falling back to the single-level solve) if a hierarchy cannot be built.
  bool BuildEnvGmg()
  {
    const int order = iodata.solver.order;
    const int dim = parent.Dimension();
    const int mg_levels = iodata.solver.linear.mg_max_levels;

    // Identify essential boundary attributes on the submesh (Gamma + env Dirichlet) from a
    // marker transfer onto a scratch order-p space: a boundary attribute is essential iff
    // all of its true DOFs lie in the essential set.
    auto scratch = std::make_unique<mfem::ParFiniteElementSpace>(env_submesh, fec.get());
    std::set<int> ess_set;
    {
      mfem::Vector t(nt);
      for (int i = 0; i < nt; i++)
      {
        t(i) = is_env_int[i] ? 0.0 : 1.0;
      }
      mfem::ParGridFunction pg(&parent_fes), sg(scratch.get());
      pg.SetFromTrueDofs(t);
      sg = 0.0;
      env_submesh->Transfer(pg, sg);
      mfem::Vector st(scratch->GetTrueVSize());
      sg.GetTrueDofs(st);
      for (int i = 0; i < st.Size(); i++)
      {
        if (std::abs(st(i)) > 0.5)
        {
          ess_set.insert(i);
        }
      }
    }
    mfem::Array<int> ess_attr;
    // The essential-attribute decision must be identical on every rank (otherwise ranks
    // diverge between the GMG and single-level paths and deadlock in the collectives
    // below). An attribute is essential iff, globally, it has DOFs and none of them are
    // non-essential.
    MPI_Comm comm = env_submesh->GetComm();
    const int lbmax =
        env_submesh->bdr_attributes.Size() ? env_submesh->bdr_attributes.Max() : 0;
    int bmax = 0;
    MPI_Allreduce(&lbmax, &bmax, 1, MPI_INT, MPI_MAX, comm);
    std::vector<int> all_in(bmax, 1), has_dofs(bmax, 0);
    for (int a = 1; a <= bmax; a++)
    {
      mfem::Array<int> m(bmax);
      m = 0;
      m[a - 1] = 1;
      mfem::Array<int> adofs;
      scratch->GetEssentialTrueDofs(m, adofs);
      if (adofs.Size() > 0)
      {
        has_dofs[a - 1] = 1;
      }
      for (int d : adofs)
      {
        if (!ess_set.count(d))
        {
          all_in[a - 1] = 0;
          break;
        }
      }
    }
    std::vector<int> g_all_in(bmax), g_has(bmax);
    MPI_Allreduce(all_in.data(), g_all_in.data(), bmax, MPI_INT, MPI_LAND, comm);
    MPI_Allreduce(has_dofs.data(), g_has.data(), bmax, MPI_INT, MPI_LOR, comm);
    for (int a = 1; a <= bmax; a++)
    {
      if (g_all_in[a - 1] && g_has[a - 1])
      {
        ess_attr.Append(a);
      }
    }
    if (ess_attr.Size() == 0)
    {
      return false;
    }

    // p-multigrid hierarchy on the submesh, with ess_attr as the Dirichlet boundary.
    env_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
        order, dim, mg_levels, iodata.solver.linear.mg_coarsening, false);
    std::vector<mfem::Array<int>> &dbc_lists = env_dbc_lists;
    dbc_lists.clear();
    env_hierarchy = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
            mg_levels, env_mesh_vec, env_fecs, &ess_attr, &dbc_lists));
    if (env_hierarchy->GetNumLevels() < 2)
    {
      env_hierarchy.reset();
      env_fecs.clear();
      return false;
    }
    env_solve_fes = &env_hierarchy->GetFinestFESpace().Get();

    // Ceed-consistent material coefficient (cf. divfree.cpp): built via a MaterialOperator
    // on the submesh so the attribute-to-material map matches Palace's local CEED numbering
    // (a hand-built coefficient keyed by global attribute assembles the interior to zero).
    MaterialOperator env_mat_op(iodata, *env_mesh_vec[0]);
    MaterialPropertyCoefficient coef(env_mat_op.GetAttributeToMaterial(),
                                     env_mat_op.GetPermittivityReal());
    BilinearForm a(env_hierarchy->GetFinestFESpace());
    a.AddDomainIntegrator<DiffusionIntegrator>(coef);
    auto a_vec = a.Assemble(*env_hierarchy, false);

    const std::size_t nl = env_hierarchy->GetNumLevels();
    env_mg_op = std::make_unique<MultigridOperator>(nl);
    for (std::size_t l = 0; l < nl; l++)
    {
      auto &fes_l = env_hierarchy->GetFESpaceAtLevel(l);
      auto A_l = std::make_unique<ParOperator>(std::move(a_vec[l]), fes_l);
      A_l->SetEssentialTrueDofs(dbc_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      env_mg_op->AddOperator(std::move(A_l));
    }
    env_ess = dbc_lists.back();  // finest essential, consistent with the operators

    auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
    amg->SetDropSmallEntries(false);
    auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
        iodata, comm, std::move(amg), env_hierarchy->GetProlongationOperators());
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-12);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(1000);
    env_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
    env_ksp->SetOperators(*env_mg_op, *env_mg_op);

    env_sgf.SetSpace(env_solve_fes);
    return true;
  }

  std::unique_ptr<mfem::HypreParMatrix> AssembleEnvSubmesh()
  {
    PWMatrixCoefficient coef(parent.Dimension(), env_eps);
    mfem::ParBilinearForm a(env_sfes.get());
    if (magnetostatic)
    {
      a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(coef));
      const int am = env_submesh->attributes.Size() ? env_submesh->attributes.Max() : 1;
      mfem::Vector mass(am);
      mass = 0.0;
      for (const auto &[attr, t] : env_eps)
      {
        mass(attr - 1) = kMagRegularization;
      }
      mfem::PWConstCoefficient mcoef(mass);
      a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mcoef));
      a.Assemble();
      a.Finalize();
      return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
    }
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(coef));
    a.Assemble();
    a.Finalize();
    return std::unique_ptr<mfem::HypreParMatrix>(a.ParallelAssemble());
  }

  // Apply A_EE^-1 to a parent true-DOF vector (nonzero on the environment interior):
  // transfer to the submesh, solve the Dirichlet problem, transfer back, and restrict to
  // the environment interior.
  // Interface selection E_Gamma (parent true DOFs x global interface index): E(i, g_i) = 1
  // for each interface true DOF. Each rank's interface indices are its own contiguous
  // block, so the matrix is local (diag block only).
  std::unique_ptr<mfem::HypreParMatrix> AssembleInterfaceSelection() const
  {
    MPI_Comm comm = parent_fes.GetComm();
    std::vector<int> I(nt + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<double> V;
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        J.push_back(gamma_global[i]);
        V.push_back(1.0);
      }
      I[i + 1] = static_cast<int>(J.size());
    }
    // Column (interface) partitioning in the layout HYPRE expects for this build.
    std::vector<HYPRE_BigInt> cols;
    if (HYPRE_AssumedPartitionCheck())
    {
      cols = {gamma_off, gamma_off + gamma_nloc, nG_global};
    }
    else
    {
      const int nranks = Mpi::Size(comm);
      std::vector<int> cnt(nranks);
      MPI_Allgather(&gamma_nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
      cols.assign(nranks + 1, 0);
      for (int r = 0; r < nranks; r++)
      {
        cols[r + 1] = cols[r] + cnt[r];
      }
    }
    return std::make_unique<mfem::HypreParMatrix>(
        comm, nt, parent_fes.GlobalTrueVSize(), static_cast<HYPRE_BigInt>(nG_global),
        I.data(), J.data(), V.data(), parent_fes.GetTrueDofOffsets(), cols.data());
  }

  // S_E as a parent-space HypreParMatrix (dense Gamma x Gamma block, zero elsewhere), for
  // factoring the exact condensed region operator. Each rank's S_E rows are its own
  // interface true DOFs; columns map interface index -> global true DOF via a replicated
  // table.
  std::unique_ptr<mfem::HypreParMatrix> AssembleDenseInterface() const
  {
    MPI_Comm comm = parent_fes.GetComm();
    const int nG = nG_global, nloc = gamma_nloc;
    const HYPRE_BigInt tstart = parent_fes.GetMyTDofOffset();
    std::vector<HYPRE_BigInt> mine(nloc), g2t(nG);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        mine[gamma_global[i] - gamma_off] = tstart + i;
      }
    }
    const int nranks = Mpi::Size(comm);
    std::vector<int> cnt(nranks), disp(nranks);
    MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
    for (int r = 0, acc = 0; r < nranks; r++)
    {
      disp[r] = acc;
      acc += cnt[r];
    }
    MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, g2t.data(), cnt.data(),
                   disp.data(), HYPRE_MPI_BIG_INT, comm);
    std::vector<int> I(nt + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<double> V;
    J.reserve(static_cast<std::size_t>(nloc) * nG);
    V.reserve(static_cast<std::size_t>(nloc) * nG);
    for (int i = 0; i < nt; i++)
    {
      if (is_gamma[i])
      {
        const double *row =
            &S_rows[static_cast<std::size_t>(gamma_global[i] - gamma_off) * nG];
        for (int j = 0; j < nG; j++)
        {
          J.push_back(g2t[j]);
          V.push_back(row[j]);
        }
      }
      I[i + 1] = static_cast<int>(J.size());
    }
    const HYPRE_BigInt glob = parent_fes.GlobalTrueVSize();
    return std::make_unique<mfem::HypreParMatrix>(comm, nt, glob, glob, I.data(), J.data(),
                                                  V.data(), parent_fes.GetTrueDofOffsets(),
                                                  parent_fes.GetTrueDofOffsets());
  }

  // Batched environment interior solves y_k = A_EE^-1 x_k (inputs/outputs are parent
  // true-DOF vectors; only the environment interior of y_k is kept, the rest is zeroed).
  // The direct path runs multi-RHS blocks through the single factor; the iterative path
  // loops. y_k may alias x_k.
  void ApplyAeeInvMulti(const std::vector<const mfem::Vector *> &X,
                        const std::vector<mfem::Vector *> &Y) const
  {
    MFEM_ASSERT(X.size() == Y.size(), "ApplyAeeInvMulti size mismatch!");
    MFEM_VERIFY(env_built, "Environment solver used before EnsureEnv()!");
#if defined(MFEM_USE_SUPERLU)
    if (env_parent_direct)
    {
      if (env_block == 0)
      {
        env_block = std::max(1, std::min(kMaterializeBlock, static_cast<int>(X.size())));
        env_pad_in.SetSize(nt);
        env_pad_in = 0.0;
        env_pad_out.SetSize(nt);
      }
      // The factored operator is block diagonal [A_EE, 0; 0, I] (non-interior rows/cols
      // eliminated), so the interior solution depends only on the interior rhs: the inputs
      // are passed as-is and only the output is masked. The tail block is padded with a
      // shared zero rhs (SuperLU fixes the RHS count at the first solve).
      mfem::Array<const mfem::Vector *> Xp(env_block);
      mfem::Array<mfem::Vector *> Yp(env_block);
      const int n = static_cast<int>(X.size());
      for (int c0 = 0; c0 < n; c0 += env_block)
      {
        const int nb = std::min(env_block, n - c0);
        for (int k = 0; k < env_block; k++)
        {
          Xp[k] = (k < nb) ? X[c0 + k] : &env_pad_in;
          Yp[k] = (k < nb) ? Y[c0 + k] : &env_pad_out;
          Yp[k]->SetSize(nt);
        }
        if (env_block == 1)
        {
          env_lu_parent->Mult(*Xp[0], *Yp[0]);
        }
        else
        {
          env_lu_parent->ArrayMult(Xp, Yp);
        }
        for (int k = 0; k < nb; k++)
        {
          mfem::Vector &y = *Y[c0 + k];
          for (int i = 0; i < nt; i++)
          {
            if (!is_env_int[i])
            {
              y(i) = 0.0;
            }
          }
        }
      }
      return;
    }
#endif
    for (std::size_t k = 0; k < X.size(); k++)
    {
      ApplyAeeInv(*X[k], *Y[k]);
    }
  }

  void ApplyAeeInv(const mfem::Vector &x_parent, mfem::Vector &y_parent) const
  {
#if defined(MFEM_USE_SUPERLU)
    if (env_parent_direct)
    {
      // Parent-space direct env interior solve (no submesh transfer), via the batched path
      // so the factor's fixed RHS count is respected.
      ApplyAeeInvMulti({&x_parent}, {&y_parent});
      return;
    }
#endif
    env_pgf.SetFromTrueDofs(x_parent);
    env_sgf = 0.0;
    env_submesh->Transfer(env_pgf, env_sgf);
    env_sgf.GetTrueDofs(env_srhs);
    for (int i = 0; i < env_ess.Size(); i++)
    {
      env_srhs(env_ess[i]) = 0.0;
    }
    env_ssol = 0.0;
    env_ksp->Mult(env_srhs, env_ssol);
    env_sgf.SetFromTrueDofs(env_ssol);
    env_pgf = 0.0;
    env_submesh->Transfer(env_sgf, env_pgf);
    y_parent.SetSize(nt);
    env_pgf.GetTrueDofs(y_parent);
    for (int i = 0; i < nt; i++)
    {
      if (!is_env_int[i])
      {
        y_parent(i) = 0.0;
      }
    }
  }

  // Region-solve geometric multigrid preconditioner on the region submesh (order>=2 H1).
  // Returns false to fall back to the single-level wrapped BoomerAMG preconditioner.
  bool BuildRegionGmg()
  {
    const int order = iodata.solver.order;
    const int dim = parent.Dimension();
    const int mg_levels = iodata.solver.linear.mg_max_levels;
    if (magnetostatic || order <= 1 || mg_levels <= 1)
    {
      return false;
    }

    reg_mesh_vec.clear();
    reg_mesh_vec.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromDomain(parent, ra_arr))));
    reg_submesh = dynamic_cast<mfem::ParSubMesh *>(&reg_mesh_vec[0]->Get());
    reg_mesh_vec[0]->RebuildCeedAttributes();
    reg_pgf.SetSpace(&parent_fes);

    // Essential submesh DOFs = the region Dirichlet terminals (Gamma stays free), found by
    // transferring the parent Dirichlet marker onto a scratch order-p space.
    auto scratch = std::make_unique<mfem::ParFiniteElementSpace>(reg_submesh, fec.get());
    std::set<int> ess_set;
    {
      mfem::Vector t(nt);
      t = 0.0;
      for (int i = 0; i < dbc_tdofs.Size(); i++)
      {
        t(dbc_tdofs[i]) = 1.0;
      }
      mfem::ParGridFunction pg(&parent_fes), sg(scratch.get());
      pg.SetFromTrueDofs(t);
      sg = 0.0;
      reg_submesh->Transfer(pg, sg);
      mfem::Vector st(scratch->GetTrueVSize());
      sg.GetTrueDofs(st);
      for (int i = 0; i < st.Size(); i++)
      {
        if (std::abs(st(i)) > 0.5)
        {
          ess_set.insert(i);
        }
      }
    }
    // Global essential-attribute decision (identical on every rank, see BuildEnvGmg).
    MPI_Comm comm = reg_submesh->GetComm();
    const int lbmax =
        reg_submesh->bdr_attributes.Size() ? reg_submesh->bdr_attributes.Max() : 0;
    int bmax = 0;
    MPI_Allreduce(&lbmax, &bmax, 1, MPI_INT, MPI_MAX, comm);
    std::vector<int> all_in(bmax, 1), has_dofs(bmax, 0);
    for (int a = 1; a <= bmax; a++)
    {
      mfem::Array<int> m(bmax);
      m = 0;
      m[a - 1] = 1;
      mfem::Array<int> adofs;
      scratch->GetEssentialTrueDofs(m, adofs);
      if (adofs.Size() > 0)
      {
        has_dofs[a - 1] = 1;
      }
      for (int d : adofs)
      {
        if (!ess_set.count(d))
        {
          all_in[a - 1] = 0;
          break;
        }
      }
    }
    std::vector<int> g_all_in(bmax), g_has(bmax);
    MPI_Allreduce(all_in.data(), g_all_in.data(), bmax, MPI_INT, MPI_LAND, comm);
    MPI_Allreduce(has_dofs.data(), g_has.data(), bmax, MPI_INT, MPI_LOR, comm);
    mfem::Array<int> ess_attr;
    for (int a = 1; a <= bmax; a++)
    {
      if (g_all_in[a - 1] && g_has[a - 1])
      {
        ess_attr.Append(a);
      }
    }
    if (ess_attr.Size() == 0)
    {
      return false;  // pure-Neumann region block: keep the single-level preconditioner
    }

    reg_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
        order, dim, mg_levels, iodata.solver.linear.mg_coarsening, false);
    reg_dbc_lists.clear();
    reg_hierarchy = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
            mg_levels, reg_mesh_vec, reg_fecs, &ess_attr, &reg_dbc_lists));
    if (reg_hierarchy->GetNumLevels() < 2)
    {
      reg_hierarchy.reset();
      reg_fecs.clear();
      return false;
    }
    reg_solve_fes = &reg_hierarchy->GetFinestFESpace().Get();

    MaterialOperator reg_mat_op(iodata, *reg_mesh_vec[0]);
    MaterialPropertyCoefficient coef(reg_mat_op.GetAttributeToMaterial(),
                                     reg_mat_op.GetPermittivityReal());
    BilinearForm a(reg_hierarchy->GetFinestFESpace());
    a.AddDomainIntegrator<DiffusionIntegrator>(coef);
    auto a_vec = a.Assemble(*reg_hierarchy, false);

    const std::size_t nl = reg_hierarchy->GetNumLevels();
    reg_mg_op = std::make_unique<MultigridOperator>(nl);
    for (std::size_t l = 0; l < nl; l++)
    {
      auto &fes_l = reg_hierarchy->GetFESpaceAtLevel(l);
      auto A_l = std::make_unique<ParOperator>(std::move(a_vec[l]), fes_l);
      A_l->SetEssentialTrueDofs(reg_dbc_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      reg_mg_op->AddOperator(std::move(A_l));
    }
    reg_ess = reg_dbc_lists.back();

    auto amg = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, true, 0));
    amg->SetDropSmallEntries(false);
    auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
        iodata, comm, std::move(amg), reg_hierarchy->GetProlongationOperators());
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-10);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(1000);
    reg_gmg_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
    reg_gmg_ksp->SetOperators(*reg_mg_op, *reg_mg_op);

    reg_sgf.SetSpace(reg_solve_fes);
    reg_srhs.SetSize(reg_solve_fes->GetTrueVSize());
    reg_ssol.SetSize(reg_solve_fes->GetTrueVSize());
    region_gmg = true;
    return true;
  }

  // Region preconditioner apply (approximates A_region_free^-1): identity on the
  // non-region-free DOFs (A_region_free is identity there), and the region-submesh GMG
  // solve on the region-free DOFs (transfer parent -> submesh, solve, transfer back,
  // restrict).
  void ApplyRegionGmgPc(const mfem::Vector &r, mfem::Vector &z) const
  {
    z.SetSize(nt);
    for (int i = 0; i < nt; i++)
    {
      z(i) = is_region_free[i] ? 0.0 : r(i);
    }
    reg_pgf.SetFromTrueDofs(r);
    reg_sgf = 0.0;
    reg_submesh->Transfer(reg_pgf, reg_sgf);
    reg_sgf.GetTrueDofs(reg_srhs);
    for (int i = 0; i < reg_ess.Size(); i++)
    {
      reg_srhs(reg_ess[i]) = 0.0;
    }
    reg_ssol = 0.0;
    reg_gmg_ksp->Mult(reg_srhs, reg_ssol);
    reg_sgf.SetFromTrueDofs(reg_ssol);
    reg_pgf = 0.0;
    reg_submesh->Transfer(reg_sgf, reg_pgf);
    mfem::Vector zr(nt);
    reg_pgf.GetTrueDofs(zr);
    for (int i = 0; i < nt; i++)
    {
      if (is_region_free[i])
      {
        z(i) = zr(i);
      }
    }
  }
};

SubstructuringSolver::SubstructuringSolver(const IoData &iodata,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh)
  : impl(std::make_unique<Impl>(iodata, mesh.back()->Get()))
{
  MFEM_VERIFY(iodata.solver.substructuring,
              "SubstructuringSolver requires a Solver.Substructuring configuration!");
}

SubstructuringSolver::~SubstructuringSolver() = default;

void SubstructuringSolver::CondenseEnvironment()
{
  // The environment interior solver is built lazily (EnsureEnv): an online run that loads
  // S_E and the terminal modes computes capacitance without ever factoring the environment.
  Impl *pi = impl.get();
  impl->dtn = std::make_unique<ImplicitDtN>(
      *impl->A_env, [pi](const mfem::Vector &x, mfem::Vector &y) { pi->ApplyAeeInv(x, y); },
      impl->is_gamma, impl->is_env_int);

  // A_region restricted to region-free true DOFs (non-region-free pinned to identity).
  impl->A_region_free = std::make_unique<mfem::HypreParMatrix>(*impl->A_region);
  mfem::Array<int> non_rfree;
  for (int i = 0; i < impl->nt; i++)
  {
    if (!impl->is_region_free[i])
    {
      non_rfree.Append(i);
    }
  }
  {
    std::unique_ptr<mfem::HypreParMatrix> tmp(
        impl->A_region_free->EliminateRowsCols(non_rfree));
  }

  // Materialize the interface operator S_E and load g_E once, so region solves reuse them
  // without any environment solves. Global interface enumeration via MPI_Exscan.
  MPI_Comm comm = impl->parent_fes.GetComm();
  int nloc = 0;
  for (int i = 0; i < impl->nt; i++)
  {
    if (impl->is_gamma[i])
    {
      nloc++;
    }
  }
  int off = 0;
  MPI_Exscan(&nloc, &off, 1, MPI_INT, MPI_SUM, comm);
  impl->nG_global = 0;
  MPI_Allreduce(&nloc, &impl->nG_global, 1, MPI_INT, MPI_SUM, comm);
  impl->gamma_global.assign(impl->nt, -1);
  {
    int c = off;
    for (int i = 0; i < impl->nt; i++)
    {
      if (impl->is_gamma[i])
      {
        impl->gamma_global[i] = c++;
      }
    }
  }

  const int nG = impl->nG_global;
  impl->gamma_off = off;
  impl->gamma_nloc = nloc;
  impl->S_rows.assign(static_cast<std::size_t>(nloc) * nG, 0.0);
  // Per-rank row counts/displacements for gather/scatter of the distributed S_E (row
  // blocks).
  const int nranks = Mpi::Size(comm);
  std::vector<int> row_cnt(nranks), row_disp(nranks);
  {
    std::vector<int> nloc_all(nranks);
    MPI_Allgather(&nloc, 1, MPI_INT, nloc_all.data(), 1, MPI_INT, comm);
    int acc = 0;
    for (int r = 0; r < nranks; r++)
    {
      row_cnt[r] = nloc_all[r] * nG;
      row_disp[r] = acc;
      acc += row_cnt[r];
    }
  }

  // Interface true-DOF geometric signature (replicated), for re-ordering a saved S_E onto
  // the current interface after region re-meshing / re-partitioning (fixed Gamma). H1: DOF
  // coordinates (width 3). H(curl): per-edge moments dof(e_b) (edge vector) and dof(x_a
  // e_b) (edge_b * mid_a), width 12 -- position + signed orientation, matching edge DOFs up
  // to an orientation flip.
  const int sdim = impl->parent.Dimension();
  const int sig_w = impl->magnetostatic ? 12 : 3;
  auto gamma_sig = [&]()
  {
    std::vector<double> loc(static_cast<std::size_t>(nG) * sig_w, 0.0),
        glob(static_cast<std::size_t>(nG) * sig_w, 0.0);
    mfem::ParGridFunction gf(&impl->parent_fes);
    Vector td(impl->nt);
    auto stamp = [&](int slot, mfem::Coefficient *sc, mfem::VectorCoefficient *vc)
    {
      if (sc)
      {
        gf.ProjectCoefficient(*sc);
      }
      else
      {
        gf.ProjectCoefficient(*vc);
      }
      gf.GetTrueDofs(td);
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_gamma[i])
        {
          loc[static_cast<std::size_t>(impl->gamma_global[i]) * sig_w + slot] = td(i);
        }
      }
    };
    if (!impl->magnetostatic)
    {
      for (int d = 0; d < sdim; d++)
      {
        mfem::FunctionCoefficient xc([d](const mfem::Vector &x) { return x(d); });
        stamp(d, &xc, nullptr);
      }
    }
    else
    {
      // Edge (tangential-integral) moments: o_b = dof(e_b) = edge vector; and
      // M[a][b] = dof(x_a e_b) = (edge vector)_b * midpoint_a. Together (12 numbers) these
      // uniquely identify an edge DOF up to an orientation flip (both are linear in the
      // tangent, so a flip negates all of them).
      for (int b = 0; b < 3; b++)
      {
        mfem::Vector e(3);
        e = 0.0;
        e(b) = 1.0;
        mfem::VectorConstantCoefficient ec(e);
        stamp(b, nullptr, &ec);
      }
      for (int a = 0; a < 3; a++)
      {
        for (int b = 0; b < 3; b++)
        {
          mfem::VectorFunctionCoefficient xc(3,
                                             [a, b](const mfem::Vector &x, mfem::Vector &v)
                                             {
                                               v = 0.0;
                                               v(b) = x(a);
                                             });
          stamp(3 + a * 3 + b, nullptr, &xc);
        }
      }
    }
    MPI_Allreduce(loc.data(), glob.data(), nG * sig_w, MPI_DOUBLE, MPI_SUM, comm);
    return glob;
  };

  // Offline/online: the interface enumeration above is cheap and deterministic to rebuild,
  // but materializing S_E costs |Gamma| environment solves. In Online mode with a saved
  // model, load S_E instead; in Offline mode with a path set, materialize and save it.
  const auto &subcfg = *impl->iodata.solver.substructuring;
  const bool online = (subcfg.mode == SubstructuringMode::ONLINE);
  const std::string &model_path = subcfg.save_model;
  const int rank = Mpi::Rank(comm);
  bool loaded = false;
  if (online && !model_path.empty())
  {
    int nG_file = -1;
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      if (f.good())
      {
        f.read(reinterpret_cast<char *>(&nG_file), sizeof(int));
      }
    }
    MPI_Bcast(&nG_file, 1, MPI_INT, 0, comm);
    MFEM_VERIFY(nG_file == nG,
                "Saved substructuring model interface size ("
                    << nG_file << ") does not match this run (" << nG
                    << "); the interface (Gamma) must be identical between the offline and "
                       "online runs.");
    int sig_type = 0;  // 0: none (identical mesh), 1: H1 coords, 2: H(curl) edge signature
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(sizeof(int));
      f.read(reinterpret_cast<char *>(&sig_type), sizeof(int));
    }
    MPI_Bcast(&sig_type, 1, MPI_INT, 0, comm);
    const int file_w = (sig_type == 2) ? 12 : (sig_type == 1 ? 3 : 0);
    std::vector<double> saved_sig;
    if (file_w > 0)
    {
      saved_sig.assign(static_cast<std::size_t>(nG) * file_w, 0.0);
      if (rank == 0)
      {
        std::ifstream f(model_path, std::ios::binary);
        f.seekg(2 * sizeof(int));
        f.read(reinterpret_cast<char *>(saved_sig.data()), sizeof(double) * nG * file_w);
      }
      MPI_Bcast(saved_sig.data(), nG * file_w, MPI_DOUBLE, 0, comm);
    }
    // Rank 0 holds the full saved S_E transiently, re-orders it onto the current interface
    // numbering, then scatters contiguous row blocks -- no replicated full matrix persists.
    std::vector<double> S_full;
    if (rank == 0)
    {
      S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(static_cast<std::streamoff>(2 * sizeof(int) + sizeof(double) * nG * file_w));
      f.read(reinterpret_cast<char *>(S_full.data()), sizeof(double) * nG * nG);
    }
    // Signature match -> (perm, sgn): online interface index g corresponds to saved index
    // perm[g] with orientation sgn[g] (sgn=1 for H1; +/-1 for H(curl)). Computed on all
    // ranks from the replicated signatures. gamma_global stays fresh (contiguous), and S_E
    // is re-indexed to the online order: S_on[i][j] = sgn_i sgn_j S_off[perm_i][perm_j].
    std::vector<int> perm(nG);
    std::iota(perm.begin(), perm.end(), 0);  // identity unless re-ordered by signature
    std::vector<double> sgn(nG, 1.0);
    if (file_w > 0)
    {
      const std::vector<double> cur = gamma_sig();
      const bool signed_match = (sig_type == 2);
      double worst = 0.0;
      for (int g = 0; g < nG; g++)
      {
        int best = 0;
        double bs = 1.0, bd = 1e300;
        for (int s = 0; s < nG; s++)
        {
          double dp = 0.0, dm = 0.0;
          for (int d = 0; d < file_w; d++)
          {
            const double a = cur[static_cast<std::size_t>(g) * file_w + d];
            const double b = saved_sig[static_cast<std::size_t>(s) * file_w + d];
            dp += (a - b) * (a - b);
            if (signed_match)
            {
              dm += (a + b) * (a + b);
            }
          }
          if (dp < bd)
          {
            bd = dp;
            best = s;
            bs = 1.0;
          }
          if (signed_match && dm < bd)
          {
            bd = dm;
            best = s;
            bs = -1.0;
          }
        }
        perm[g] = best;
        sgn[g] = bs;
        worst = std::max(worst, bd);
      }
      double gworst = 0.0;
      MPI_Allreduce(&worst, &gworst, 1, MPI_DOUBLE, MPI_MAX, comm);
      MFEM_VERIFY(std::sqrt(gworst) < 1e-8,
                  "Online interface DOFs do not match the saved model (max mismatch "
                      << std::sqrt(gworst) << "); the interface Gamma must be identical.");
    }
    if (rank == 0 && file_w > 0)
    {
      const std::vector<double> S_off = S_full;
      for (int i = 0; i < nG; i++)
      {
        for (int j = 0; j < nG; j++)
        {
          S_full[static_cast<std::size_t>(i) * nG + j] =
              sgn[i] * sgn[j] * S_off[static_cast<std::size_t>(perm[i]) * nG + perm[j]];
        }
      }
    }
    MPI_Scatterv(rank == 0 ? S_full.data() : nullptr, row_cnt.data(), row_disp.data(),
                 MPI_DOUBLE, impl->S_rows.data(), nloc * nG, MPI_DOUBLE, 0, comm);
    loaded = true;
    if (!impl->magnetostatic)
    {
      const std::streamoff modes_pos =
          static_cast<std::streamoff>(2 * sizeof(int)) +
          static_cast<std::streamoff>(sizeof(double)) * nG * file_w +
          static_cast<std::streamoff>(sizeof(double)) * nG * nG;
      impl->LoadModes(model_path, modes_pos, perm, sgn);
    }
  }
  if (!loaded)
  {
    impl->EnsureEnv();
    // Map each owned interface column (global index in [off, off+nloc)) to its local true
    // DOF, so each unit-vector RHS is set in O(1) instead of scanning all true DOFs per
    // column.
    std::vector<int> col_to_dof(nloc, -1);
    for (int i = 0; i < impl->nt; i++)
    {
      if (impl->is_gamma[i])
      {
        col_to_dof[impl->gamma_global[i] - off] = i;
      }
    }
    // S_E e_c = (A_env e_c)|_Gamma - (A_env A_EE^-1 (A_env e_c)|_E)|_Gamma. The results are
    // complete (assembled) true-DOF vectors, so each rank stores its own S_E rows directly
    // -- no interface Allreduce (O(nG^2) communication) is needed.
    auto store_col = [&](int c, const Vector &y)
    {
      for (int r = 0; r < nloc; r++)
      {
        impl->S_rows[static_cast<std::size_t>(r) * nG + c] = y(col_to_dof[r]);
      }
    };
    auto set_unit = [&](Vector &e, int c)
    {
      e = 0.0;
      if (c >= off && c < off + nloc)
      {
        e(col_to_dof[c - off]) = 1.0;
      }
    };
#if defined(MFEM_USE_SUPERLU)
    if (impl->env_parent_direct)
    {
      // Batched direct materialization: kMaterializeBlock columns per multi-RHS back-solve
      // through the (single) environment factor. This is its first solve, which fixes the
      // factor's block size at kMaterializeBlock for the online excitation batches too.
      // Thin couplings built once: G = A_env E_Gamma (column c = A_env e_c, supported next
      // to Gamma) and R = E_Gamma^T A_env (the interface rows), so each column costs O(nnz
      // near Gamma) instead of two full A_env matvecs.
      std::unique_ptr<mfem::HypreParMatrix> E = impl->AssembleInterfaceSelection();
      std::unique_ptr<mfem::HypreParMatrix> G(mfem::ParMult(impl->A_env.get(), E.get()));
      std::unique_ptr<mfem::HypreParMatrix> Et(E->Transpose());
      std::unique_ptr<mfem::HypreParMatrix> R(mfem::ParMult(Et.get(), impl->A_env.get()));
      E.reset();
      Et.reset();
      const int B = Impl::kMaterializeBlock;
      std::vector<Vector> t(B, Vector(impl->nt));
      std::vector<const mfem::Vector *> X(B);
      std::vector<mfem::Vector *> Y(B);
      Vector ec(nloc), z(nloc);
      std::vector<double> agg(static_cast<std::size_t>(B) * nloc);  // A_GG e_c, own rows
      for (int c0 = 0; c0 < nG; c0 += B)
      {
        const int nb = std::min(B, nG - c0);
        X.resize(nb);
        Y.resize(nb);
        for (int k = 0; k < nb; k++)
        {
          const int c = c0 + k;
          ec = 0.0;
          if (c >= off && c < off + nloc)
          {
            ec(c - off) = 1.0;
          }
          // A_env e_c: its env-interior part is the solve RHS, its Gamma part is A_GG e_c.
          G->Mult(ec, t[k]);
          for (int r = 0; r < nloc; r++)
          {
            agg[static_cast<std::size_t>(k) * nloc + r] = t[k](col_to_dof[r]);
          }
          X[k] = &t[k];
          Y[k] = &t[k];  // solve in place
        }
        impl->ApplyAeeInvMulti(X, Y);
        for (int k = 0; k < nb; k++)
        {
          R->Mult(t[k], z);  // (A_env y)|_Gamma = A_GE A_EE^-1 A_EG e_c, this rank's rows
          const int c = c0 + k;
          for (int r = 0; r < nloc; r++)
          {
            impl->S_rows[static_cast<std::size_t>(r) * nG + c] =
                agg[static_cast<std::size_t>(k) * nloc + r] - z(r);
          }
        }
      }
    }
    else
#endif
    {
      Vector e(impl->nt), y(impl->nt);
      for (int c = 0; c < nG; c++)
      {
        set_unit(e, c);
        impl->dtn->Mult(e, y);
        store_col(c, y);
      }
    }
    if (!impl->magnetostatic)
    {
      impl->ComputeTerminalModes();  // K batched env solves (cheap vs. the |Gamma| columns)
    }
    if (!model_path.empty())
    {
      const int sig_type = impl->magnetostatic ? 2 : 1;
      const std::vector<double> sig = gamma_sig();  // collective: all ranks participate
      // Gather the distributed row blocks to rank 0 to write the full matrix to disk.
      std::vector<double> S_full;
      if (rank == 0)
      {
        S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
      }
      MPI_Gatherv(impl->S_rows.data(), nloc * nG, MPI_DOUBLE,
                  rank == 0 ? S_full.data() : nullptr, row_cnt.data(), row_disp.data(),
                  MPI_DOUBLE, 0, comm);
      if (rank == 0)
      {
        std::ofstream f(model_path, std::ios::binary);
        f.write(reinterpret_cast<const char *>(&nG), sizeof(int));
        f.write(reinterpret_cast<const char *>(&sig_type), sizeof(int));
        f.write(reinterpret_cast<const char *>(sig.data()), sizeof(double) * nG * sig_w);
        f.write(reinterpret_cast<const char *>(S_full.data()), sizeof(double) * nG * nG);
      }
      if (impl->modes_ready)
      {
        impl->AppendModes(model_path);  // collective
      }
    }
  }
  // Optional hierarchical (HODLR) off-diagonal compression of the assembled S_E: the DtN's
  // well-separated interface-block couplings are low-rank, so this compresses S_E storage +
  // apply to a controlled relative tolerance while the near-field / diagonal stays exact
  // (unlike a global low-rank, which fails -- S_E is full rank). Built on rank 0 from the
  // gathered dense S_E, then broadcast (compressed) and applied replicated on every rank.
  {
    const double hodlr_tol = impl->iodata.solver.substructuring->interface_offdiag_tol;
    if (hodlr_tol > 0.0 && nG > 0)
    {
      // Interface DOF coordinates (replicated) for the coordinate-median clustering. Each
      // true DOF is placed at its interpolation point: H1 -> node coordinate, Nedelec ->
      // edge midpoint. Computed per element (GetNodes mapped through the element
      // transformation), reduced onto the owning rank's true DOF, then summed to a
      // replicated array.
      std::vector<double> coords_loc(static_cast<std::size_t>(nG) * 3, 0.0),
          coords(static_cast<std::size_t>(nG) * 3, 0.0);
      std::vector<double> tdof_xyz(static_cast<std::size_t>(impl->nt) * 3, 0.0);
      std::vector<char> have(impl->nt, 0);
      mfem::Array<int> edofs;
      mfem::Vector phys;
      for (int e = 0; e < impl->parent_fes.GetNE(); e++)
      {
        const mfem::FiniteElement *fe = impl->parent_fes.GetFE(e);
        mfem::ElementTransformation *T = impl->parent.GetElementTransformation(e);
        const mfem::IntegrationRule &nodes = fe->GetNodes();
        impl->parent_fes.GetElementDofs(e, edofs);
        for (int j = 0; j < edofs.Size(); j++)
        {
          const int ldof = edofs[j] >= 0 ? edofs[j] : -1 - edofs[j];
          const int t = impl->parent_fes.GetLocalTDofNumber(ldof);
          if (t < 0 || have[t])
          {
            continue;
          }
          T->Transform(nodes.IntPoint(j), phys);
          for (int d = 0; d < phys.Size() && d < 3; d++)
          {
            tdof_xyz[static_cast<std::size_t>(t) * 3 + d] = phys(d);
          }
          have[t] = 1;
        }
      }
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_gamma[i])
        {
          for (int d = 0; d < 3; d++)
          {
            coords_loc[static_cast<std::size_t>(impl->gamma_global[i]) * 3 + d] =
                tdof_xyz[static_cast<std::size_t>(i) * 3 + d];
          }
        }
      }
      MPI_Allreduce(coords_loc.data(), coords.data(), nG * 3, MPI_DOUBLE, MPI_SUM, comm);
      // Gather the dense S_E to rank 0 and build the Hodlr there.
      std::vector<double> S_full;
      if (rank == 0)
      {
        S_full.assign(static_cast<std::size_t>(nG) * nG, 0.0);
      }
      MPI_Gatherv(impl->S_rows.data(), nloc * nG, MPI_DOUBLE,
                  rank == 0 ? S_full.data() : nullptr, row_cnt.data(), row_disp.data(),
                  MPI_DOUBLE, 0, comm);
      std::vector<double> buf;
      int buflen = 0;
      if (rank == 0)
      {
        Hodlr h;
        h.n = nG;
        h.perm.assign(nG, 0);
        std::vector<int> idx(nG);
        std::iota(idx.begin(), idx.end(), 0);
        BuildHodlr(S_full, nG, idx, 0, coords, hodlr_tol, 32, h);
        buf = h.Serialize();
        buflen = static_cast<int>(buf.size());
      }
      MPI_Bcast(&buflen, 1, MPI_INT, 0, comm);
      buf.resize(buflen);
      MPI_Bcast(buf.data(), buflen, MPI_DOUBLE, 0, comm);
      impl->hodlr = std::make_unique<Hodlr>(Hodlr::Deserialize(buf));
      // The compressed operator supplants the dense row blocks; free them to realize the
      // memory saving.
      const long long ret = impl->hodlr->Storage();
      impl->S_rows.clear();
      impl->S_rows.shrink_to_fit();
      Mpi::Print("[HODLR] tol={:.1e}: S_E storage {:d}/{:d} = {:.3f} of dense\n", hodlr_tol,
                 ret, static_cast<long long>(nG) * nG,
                 static_cast<double>(ret) / (static_cast<double>(nG) * nG));
    }
  }

  // g_E is excitation-dependent; it is computed per excitation in the region solve.
  impl->mat_dtn =
      std::make_unique<MaterializedDtN>(impl->S_rows, impl->gamma_off, impl->gamma_global,
                                        impl->nG_global, comm, impl->hodlr.get());

  // Region-condensed solver: Palace CG preconditioned by a wrapped AMS (H(curl)) or
  // BoomerAMG (H1) on the region-free block. Built once and reused across excitations (the
  // condensed interface operator S_E is fixed). The region operator carries the same mass
  // regularization as the environment, so it is definite and AMS is used without the
  // singular-problem option.
  {
    std::unique_ptr<Solver<Operator>> pc;
    bool use_direct = false;
#if defined(MFEM_USE_SUPERLU)
    // Default to a direct factorization of A_region_free when the region is small enough to
    // factor (a near-exact preconditioner -> the outer CG converges in a few iterations);
    // fall back to GMG / AMG for very large regions.
    {
      long long reg_loc = 0;
      for (int i = 0; i < impl->nt; i++)
      {
        if (impl->is_region_free[i])
        {
          reg_loc++;
        }
      }
      long long reg_glob = 0;
      MPI_Allreduce(&reg_loc, &reg_glob, 1, MPI_LONG_LONG, MPI_SUM, comm);
      use_direct = (reg_glob <= Impl::kDirectMaxDofs);
    }
#endif
#if defined(MFEM_USE_SUPERLU)
    if (use_direct)
    {
      // Direct factorization of the region preconditioner. When the interface is small
      // enough for a dense Gamma block (and S_E is dense, not HODLR), factor the EXACT
      // condensed operator A_region_free + S_E, so the outer CG converges in one iteration:
      // the non-local DtN cannot be sparsified (truncating it to the interface stencil
      // makes the preconditioner far worse than omitting it), and A_region_free alone needs
      // O(50) iterations. Otherwise fall back to A_region_free (S_E omitted).
      impl->reg_lu = std::make_unique<SuperLUSolver>(impl->iodata, comm, 0);
      if (!impl->hodlr && nG > 0 && nG <= Impl::kDirectCondensedMaxInterface)
      {
        std::unique_ptr<mfem::HypreParMatrix> S_mat = impl->AssembleDenseInterface();
        std::unique_ptr<mfem::HypreParMatrix> P(
            mfem::ParAdd(impl->A_region_free.get(), S_mat.get()));
        S_mat.reset();
        impl->reg_lu->SetOperator(*P);  // copies P into SuperLU's distributed format
      }
      else
      {
        impl->reg_lu->SetOperator(*impl->A_region_free);
      }
      pc = std::make_unique<CallableSolver>(impl->A_region_free->Height(),
                                            [pi](const mfem::Vector &r, mfem::Vector &z)
                                            { pi->reg_lu->Mult(r, z); });
    }
    else
#endif
        if (impl->BuildRegionGmg())
    {
      // Order>=2 H1: geometric multigrid on the region submesh, routed through the region
      // preconditioner apply.
      pc = std::make_unique<CallableSolver>(impl->A_region_free->Height(),
                                            [pi](const mfem::Vector &r, mfem::Vector &z)
                                            { pi->ApplyRegionGmgPc(r, z); });
    }
    else if (impl->magnetostatic)
    {
      auto ams = std::make_unique<mfem::HypreAMS>(&impl->parent_fes);
      ams->SetPrintLevel(0);
      pc = std::make_unique<MfemWrapperSolver<Operator>>(
          std::move(ams), /*save_assembled=*/true, /*complex_matrix=*/false,
          /*drop_small_entries=*/false);
    }
    else
    {
      pc = std::make_unique<MfemWrapperSolver<Operator>>(
          std::make_unique<BoomerAmgSolver>(1, 1, true, 0), /*save_assembled=*/true,
          /*complex_matrix=*/false, /*drop_small_entries=*/false);
    }
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-10);
    pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
    pcg->SetMaxIter(2000);
    impl->region_op = std::make_unique<RegionCondensedOperator>(
        *impl->A_region_free, *impl->mat_dtn, impl->is_region_free);
    impl->region_ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
    impl->region_ksp->SetOperators(*impl->region_op, *impl->A_region_free);
  }
}

namespace
{

// Terminal excitation Dirichlet field: the driven terminal at 1 V, all others grounded.
Vector TerminalDbc(const std::map<int, std::vector<int>> &terminal_tdofs, int drive_idx,
                   int nt)
{
  Vector dbc(nt);
  dbc = 0.0;
  for (const auto &[idx, dofs] : terminal_tdofs)
  {
    const double value = (idx == drive_idx) ? 1.0 : 0.0;
    for (int d : dofs)
    {
      dbc(d) = value;
    }
  }
  return dbc;
}

}  // namespace

Vector SubstructuringSolver::SolveExcitation(int drive_idx)
{
  return SolveExcitations({drive_idx})[0];
}

std::vector<Vector>
SubstructuringSolver::SolveExcitations(const std::vector<int> &drive_terminal_indices)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  std::vector<Vector> dbcs;
  dbcs.reserve(drive_terminal_indices.size());
  for (int idx : drive_terminal_indices)
  {
    dbcs.push_back(TerminalDbc(impl->terminal_tdofs, idx, impl->nt));
  }
  return SolveDirichletBatch(dbcs);
}

Vector SubstructuringSolver::SolveDirichlet(const Vector &dbc_values)
{
  return SolveDirichlets({dbc_values})[0];
}

std::vector<Vector>
SubstructuringSolver::SolveDirichlets(const std::vector<Vector> &dbc_values)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  // Prescribe an arbitrary boundary field on the Dirichlet DOF set (e.g. a flux-loop lift).
  std::vector<Vector> dbcs(dbc_values.size(), Vector(impl->nt));
  for (std::size_t k = 0; k < dbc_values.size(); k++)
  {
    dbcs[k] = 0.0;
    for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
    {
      const int d = impl->dbc_tdofs[i];
      dbcs[k](d) = dbc_values[k](d);
    }
  }
  return SolveDirichletBatch(dbcs);
}

std::vector<Vector>
SubstructuringSolver::SolveDirichletBatch(const std::vector<Vector> &dbcs)
{
  impl->EnsureEnv();  // general Dirichlet data needs environment solves
  const int nt = impl->nt;
  const int n = static_cast<int>(dbcs.size());
  MPI_Comm comm = impl->parent_fes.GetComm();

  // g_E per excitation: the interface response to the Dirichlet data,
  //   g_E = (A_env x)|_Gamma - (A_env A_EE^-1 (A_env x)|_E)|_Gamma.
  // The environment solve is skipped (exactly) when (A_env x)|_E vanishes, i.e. the
  // Dirichlet data does not touch the environment (e.g. a region terminal); the others run
  // as one batched multi-RHS solve.
  std::vector<Vector> t(n, Vector(nt)), gE(n, Vector(nt));
  std::vector<int> touches(n, 0);
  for (int k = 0; k < n; k++)
  {
    impl->A_env->Mult(dbcs[k], t[k]);
    for (int i = 0; i < nt && !touches[k]; i++)
    {
      if (impl->is_env_int[i] && t[k](i) != 0.0)
      {
        touches[k] = 1;
      }
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, touches.data(), n, MPI_INT, MPI_MAX, comm);
  }
  {
    std::vector<const mfem::Vector *> X;
    std::vector<mfem::Vector *> Y;
    std::vector<int> ks;
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        X.push_back(&t[k]);
        Y.push_back(&gE[k]);
        ks.push_back(k);
      }
    }
    if (!X.empty())
    {
      impl->ApplyAeeInvMulti(X, Y);
    }
    Vector t2(nt);
    for (int k = 0; k < n; k++)
    {
      if (touches[k])
      {
        impl->A_env->Mult(gE[k], t2);
        for (int i = 0; i < nt; i++)
        {
          gE[k](i) = t[k](i) - t2(i);
        }
      }
      else
      {
        gE[k] = t[k];
      }
    }
  }

  // Region-condensed solve per excitation using the materialized S_E (no environment solves
  // in the loop). RHS: region Dirichlet elimination + environment load g_E on the
  // interface. Both are complete (assembled) true-DOF vectors, so each rank reads its own
  // entries.
  std::vector<Vector> u(n, Vector(nt));
  {
    Vector b(nt), tr(nt);
    for (int k = 0; k < n; k++)
    {
      impl->A_region->Mult(dbcs[k], tr);
      b = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (impl->is_region_free[i])
        {
          b(i) -= tr(i);
        }
        if (impl->is_gamma[i])
        {
          b(i) -= gE[k](i);
        }
      }
      u[k] = 0.0;
      impl->region_ksp->Mult(b, u[k]);
      for (int i = 0; i < impl->dbc_tdofs.Size(); i++)
      {
        u[k](impl->dbc_tdofs[i]) = dbcs[k](impl->dbc_tdofs[i]);
      }
    }
  }

  impl->RecoverEnvInterior(u);  // batched: u_E = -A_EE^-1 (A_env u)|_E
  return u;
}

Vector SubstructuringSolver::SolveRegion()
{
  // Default single excitation: drive the lowest-index terminal.
  MFEM_VERIFY(!impl->terminal_tdofs.empty(), "No terminals configured!");
  return SolveExcitation(impl->terminal_tdofs.begin()->first);
}

Vector SubstructuringSolver::SolveSource(const Vector &f)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  impl->EnsureEnv();
  const int nt = impl->nt;

  // Environment interior source response: solve A_EE w = f_E, correction (A_env w)|_Gamma.
  Vector fE(nt), w(nt), Aw(nt);
  fE = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      fE(i) = f(i);
    }
  }
  w = 0.0;
  impl->ApplyAeeInv(fE, w);
  impl->A_env->Mult(w, Aw);  // assembled: each rank reads its own interface entries

  // RHS: region/interface source minus the environment source correction on the interface.
  Vector b(nt);
  b = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_region_free[i])
    {
      b(i) = f(i);
    }
  }
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_gamma[i])
    {
      b(i) -= Aw(i);
    }
  }

  Vector u(nt);
  u = 0.0;
  impl->region_ksp->Mult(b, u);

  // Recover environment interior: u_E = A_EE^-1 (f_E - (A_env u)|_E).
  Vector Au(nt), rhs(nt), uE(nt);
  impl->A_env->Mult(u, Au);
  rhs = 0.0;
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      rhs(i) = fE(i) - Au(i);
    }
  }
  uE = 0.0;
  impl->ApplyAeeInv(rhs, uE);
  for (int i = 0; i < nt; i++)
  {
    if (impl->is_env_int[i])
    {
      u(i) = uE(i);
    }
  }
  return u;
}

mfem::DenseMatrix
SubstructuringSolver::CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                        std::vector<Vector> *fields, int n_fields)
{
  MFEM_VERIFY(impl->mat_dtn, "CondenseEnvironment must be called before solving!");
  MFEM_VERIFY(!impl->magnetostatic,
              "CapacitanceMatrix is for electrostatic substructuring problems!");
  if (!impl->modes_ready)
  {
    impl->ComputeTerminalModes();  // model saved without modes: needs the environment once
  }
  const int nt = impl->nt;
  const int n = static_cast<int>(terminal_indices.size());
  const int K = static_cast<int>(impl->mode_ids.size());
  MPI_Comm comm = impl->parent_fes.GetComm();
  std::vector<int> col(n, -1);  // mode column of each requested terminal
  for (int j = 0; j < n; j++)
  {
    for (int k = 0; k < K; k++)
    {
      if (impl->mode_ids[k] == terminal_indices[j])
      {
        col[j] = k;
      }
    }
    MFEM_VERIFY(col[j] >= 0, "Unknown terminal index " << terminal_indices[j] << "!");
  }
  auto g = [&](int i, int k)  // g_k at owned interface true DOF i
  {
    return impl
        ->G_rows[static_cast<std::size_t>(impl->gamma_global[i] - impl->gamma_off) * K + k];
  };

  // Region-condensed solve per terminal (environment interior left at zero): the interface
  // load is the precomputed g_k, so no environment solve is needed.
  std::vector<Vector> u(n, Vector(nt)), Ku(n, Vector(nt)), Su(n, Vector(nt));
  {
    Vector b(nt), tr(nt);
    for (int j = 0; j < n; j++)
    {
      const Vector x = impl->TerminalMode(terminal_indices[j]);
      impl->A_region->Mult(x, tr);
      b = 0.0;
      for (int i = 0; i < nt; i++)
      {
        if (impl->is_region_free[i])
        {
          b(i) -= tr(i);
        }
        if (impl->is_gamma[i])
        {
          b(i) -= g(i, col[j]);
        }
      }
      u[j] = 0.0;
      impl->region_ksp->Mult(b, u[j]);
      for (int d = 0; d < impl->dbc_tdofs.Size(); d++)
      {
        const int i = impl->dbc_tdofs[d];
        u[j](i) = x(i);  // driven terminal at 1, the others grounded
      }
      impl->K_region_e->Mult(u[j], Ku[j]);
      impl->mat_dtn->Mult(u[j], Su[j]);
    }
  }

  // C_ij = u_i^T K_R u_j + [u_G,i; e_i]^T [S_E, G; G^T, Cmode] [u_G,j; e_j]: the region
  // energy plus the environment's exact condensed energy (S_E couples the interface values,
  // G the interface to the terminal potentials, Cmode the terminals through the
  // environment interior).
  mfem::DenseMatrix C(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      double a = u[i] * Ku[j];
      for (int r = 0; r < nt; r++)
      {
        if (impl->is_gamma[r])
        {
          a += u[i](r) * Su[j](r) + u[i](r) * g(r, col[j]) + g(r, col[i]) * u[j](r);
        }
      }
      C(i, j) = a;
    }
  }
  if (n > 0)
  {
    MPI_Allreduce(MPI_IN_PLACE, C.GetData(), n * n, MPI_DOUBLE, MPI_SUM, comm);
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      C(i, j) += impl->Cmode[static_cast<std::size_t>(col[i]) * K + col[j]];
    }
  }

  // Optional full fields (environment interior recovered on demand).
  if (fields)
  {
    const int nf = std::max(0, std::min(n_fields, n));
    fields->assign(u.begin(), u.begin() + nf);
    if (nf > 0)
    {
      impl->RecoverEnvInterior(*fields);
    }
  }
  return C;
}

bool SubstructuringSolver::EnvironmentFactored() const
{
  return impl->env_built;
}

std::vector<int> SubstructuringSolver::TerminalIndices() const
{
  std::vector<int> idx;
  for (const auto &[i, dofs] : impl->terminal_tdofs)
  {
    idx.push_back(i);
  }
  return idx;
}

double SubstructuringSolver::MutualEnergy(const Vector &ui, const Vector &uj) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->K_region_e->Mult(uj, t);
  impl->K_env_e->Mult(uj, t2);
  t += t2;
  double local = 0.0;
  for (int i = 0; i < impl->nt; i++)
  {
    local += ui(i) * t(i);
  }
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, impl->parent_fes.GetComm());
  return global;
}

long long int SubstructuringSolver::RegionGlobalTrueVSize() const
{
  return impl->parent_fes.GlobalTrueVSize();
}

void SubstructuringSolver::WriteParaView(const std::string &dir,
                                         const std::vector<int> &terminals,
                                         const std::vector<Vector> &fields) const
{
  MFEM_VERIFY(terminals.size() == fields.size(),
              "WriteParaView requires one field per terminal!");
  const int order = impl->iodata.solver.order;
  mfem::ParGridFunction phi(&impl->parent_fes);
  mfem::ParaViewDataCollection pv("paraview", &impl->parent);
  pv.SetPrefixPath(dir);
  pv.SetLevelsOfDetail(order);
  pv.SetHighOrderOutput(order > 1);
  pv.SetDataFormat(mfem::VTKFormat::BINARY);
  pv.RegisterField("V", &phi);
  for (std::size_t j = 0; j < fields.size(); j++)
  {
    phi.SetFromTrueDofs(fields[j]);
    pv.SetCycle(static_cast<int>(j));
    pv.SetTime(static_cast<double>(terminals[j]));
    pv.Save();
  }
}

double SubstructuringSolver::ElectrostaticEnergy(const Vector &u) const
{
  Vector t(impl->nt), t2(impl->nt);
  impl->K_region_e->Mult(u, t);
  impl->K_env_e->Mult(u, t2);
  t += t2;
  double local = 0.0;
  for (int i = 0; i < impl->nt; i++)
  {
    local += u(i) * t(i);
  }
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, impl->parent_fes.GetComm());
  return 0.5 * global;
}

}  // namespace palace
