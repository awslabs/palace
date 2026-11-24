// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_BLOCK_UTILS_HPP
#define PALACE_LINALG_BLOCK_UTILS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>
#include <_hypre_parcsr_mv.h>
#include <ceed/backend.h>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/hypre.hpp"
#include "linalg/operator.hpp"
#include "utils/communication.hpp"

//
// Helpers shared by the vertex aggregates and the blocks of the block Jacobi smoother (not
// part of the interface of the smoother): keys for sets, root-based parallel equivalence
// classes, the rows of a parallel prolongation matrix, and the element matrices of a
// libCEED operator.
//

namespace palace::blocks
{

using Key = std::uint64_t;

// Keys (64-bit hashes) for sets of aggregates

inline Key SplitMix(Key x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline Key Combine(Key h, Key v)
{
  return SplitMix(h ^ (SplitMix(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2)));
}

// Key of a set of ids (sorted and made unique in place). Never zero.
inline Key SetKey(std::vector<Key> &ids)
{
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  Key h = 0x243f6a8885a308d3ULL;
  for (Key id : ids)
  {
    h = Combine(h, id);
  }
  return h ? h : 1;
}

// ---------------------------------------------------------------------------------------
// Small parallel helpers (root-based: the data are the strong ties between mesh vertices
// and the census of the blocks)

template <typename T>
MPI_Datatype MpiType();
template <>
inline MPI_Datatype MpiType<Key>()
{
  return MPI_UINT64_T;
}
template <>
inline MPI_Datatype MpiType<double>()
{
  return MPI_DOUBLE;
}
template <>
inline MPI_Datatype MpiType<int>()
{
  return MPI_INT;
}
template <>
inline MPI_Datatype MpiType<long long>()
{
  return MPI_LONG_LONG;
}

// The counts and displacements of these helpers are of type int: verifies, on all
// processes, that the global number of entries of a gather, scatter or broadcast is in its
// range. Collective.
inline void VerifyCount(MPI_Comm comm, long long local, bool sum)
{
  long long total = local;
  if (sum)
  {
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, comm);
  }
  else
  {
    MPI_Bcast(&total, 1, MPI_LONG_LONG, 0, comm);
  }
  MFEM_VERIFY(total <= std::numeric_limits<int>::max(),
              "Block Jacobi: the setup of the blocks gathers "
                  << total
                  << " entries on one process, which exceeds the range of the MPI counts!");
}

// Gather the local vectors on the root (returns the concatenation there, in rank order, and
// the counts per rank).
template <typename T>
std::vector<T> GatherRoot(MPI_Comm comm, const std::vector<T> &loc,
                          std::vector<int> &counts)
{
  const int np = Mpi::Size(comm), rank = Mpi::Rank(comm);
  VerifyCount(comm, static_cast<long long>(loc.size()), true);
  int n = static_cast<int>(loc.size());
  counts.assign(np, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  std::vector<int> displs(np + 1, 0);
  std::vector<T> all;
  if (rank == 0)
  {
    for (int i = 0; i < np; i++)
    {
      displs[i + 1] = displs[i] + counts[i];
    }
    all.resize(displs[np]);
  }
  MPI_Gatherv(loc.data(), n, MpiType<T>(), all.data(), counts.data(), displs.data(),
              MpiType<T>(), 0, comm);
  return all;
}

// Scatter from the root, to every process stride entries for each of its counts (as
// returned by GatherRoot).
template <typename T>
std::vector<T> ScatterRoot(MPI_Comm comm, const std::vector<T> &all,
                           const std::vector<int> &counts, int stride, int nloc)
{
  const int np = Mpi::Size(comm), rank = Mpi::Rank(comm);
  VerifyCount(comm, static_cast<long long>(all.size()), false);
  std::vector<int> scounts(np, 0), displs(np + 1, 0);
  if (rank == 0)
  {
    for (int i = 0; i < np; i++)
    {
      scounts[i] = stride * counts[i];
      displs[i + 1] = displs[i] + scounts[i];
    }
    MFEM_VERIFY(static_cast<std::size_t>(displs[np]) == all.size(),
                "Block Jacobi: unexpected size of the data to scatter!");
  }
  std::vector<T> loc(stride * static_cast<std::size_t>(nloc));
  MPI_Scatterv(all.data(), scounts.data(), displs.data(), MpiType<T>(), loc.data(),
               static_cast<int>(loc.size()), MpiType<T>(), 0, comm);
  return loc;
}

template <typename T>
void BroadcastRoot(MPI_Comm comm, std::vector<T> &v)
{
  long long n = static_cast<long long>(v.size());
  MPI_Bcast(&n, 1, MPI_LONG_LONG, 0, comm);
  MFEM_VERIFY(n <= std::numeric_limits<int>::max(),
              "Block Jacobi: the setup of the blocks broadcasts "
                  << n << " entries, which exceeds the range of the MPI counts!");
  v.resize(n);
  if (n > 0)
  {
    MPI_Bcast(v.data(), static_cast<int>(n), MpiType<T>(), 0, comm);
  }
}

// Union-find on keys.
struct KeyUnion
{
  std::unordered_map<Key, Key> parent;

  Key Find(Key a)
  {
    auto it = parent.find(a);
    if (it == parent.end())
    {
      parent.emplace(a, a);
      return a;
    }
    Key r = a;
    while (true)
    {
      const Key p = parent[r];
      if (p == r)
      {
        break;
      }
      r = p;
    }
    while (a != r)
    {
      const Key p = parent[a];
      parent[a] = r;
      a = p;
    }
    return r;
  }

  void Union(Key a, Key b)
  {
    a = Find(a);
    b = Find(b);
    if (a != b)
    {
      // The smaller key is the representative.
      if (a < b)
      {
        parent[b] = a;
      }
      else
      {
        parent[a] = b;
      }
    }
  }

  // (member, representative) for the members of classes with more than one member.
  std::vector<Key> Classes()
  {
    std::unordered_map<Key, int> size;
    for (auto &[k, p] : parent)
    {
      size[Find(k)]++;
    }
    std::vector<Key> out;
    for (auto &[k, p] : parent)
    {
      const Key r = Find(k);
      if (size[r] > 1)
      {
        out.push_back(k);
        out.push_back(r);
      }
    }
    return out;
  }
};

// Classes of the union of the pairs of all processes: on every process, the map from the
// members of the classes with more than one member to their representative.
inline std::unordered_map<Key, Key> GlobalClasses(MPI_Comm comm,
                                                  const std::vector<Key> &pairs)
{
  std::vector<int> counts;
  std::vector<Key> all = GatherRoot(comm, pairs, counts);
  std::vector<Key> flat;
  if (Mpi::Root(comm))
  {
    KeyUnion uf;
    for (std::size_t k = 0; k + 1 < all.size(); k += 2)
    {
      uf.Union(all[k], all[k + 1]);
    }
    flat = uf.Classes();
  }
  BroadcastRoot(comm, flat);
  std::unordered_map<Key, Key> map;
  map.reserve(flat.size());
  for (std::size_t k = 0; k + 1 < flat.size(); k += 2)
  {
    map.emplace(flat[k], flat[k + 1]);
  }
  return map;
}

// ---------------------------------------------------------------------------------------

// Rows of the parallel prolongation (L-dofs x true dofs)

struct PRows
{
  hypre_ParCSRMatrix *h;
  const HYPRE_Int *di, *dj, *oi, *oj;
  const HYPRE_Real *da, *oa;
  const HYPRE_BigInt *cmap;
  HYPRE_BigInt col0;
  int nrows, ncols, noffd;

  PRows(const mfem::HypreParMatrix *P)
  {
    MFEM_VERIFY(P, "Missing prolongation matrix!");
    P->HostRead();
    h = *P;
    hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(h), *offd = hypre_ParCSRMatrixOffd(h);
    di = hypre_CSRMatrixI(diag);
    dj = hypre_CSRMatrixJ(diag);
    da = hypre_CSRMatrixData(diag);
    oi = hypre_CSRMatrixI(offd);
    oj = hypre_CSRMatrixJ(offd);
    oa = hypre_CSRMatrixData(offd);
    cmap = hypre_ParCSRMatrixColMapOffd(h);
    col0 = hypre_ParCSRMatrixFirstColDiag(h);
    nrows = hypre_CSRMatrixNumRows(diag);
    ncols = hypre_CSRMatrixNumCols(diag);
    noffd = hypre_CSRMatrixNumCols(offd);
  }

  int NumOffdRow(int i) const { return (noffd > 0 && oi) ? oi[i + 1] - oi[i] : 0; }

  // A row with a single unit entry: the L-dof is a true dof (of this or another process).
  bool IsUnit(int i) const
  {
    const int nd = di[i + 1] - di[i], no = NumOffdRow(i);
    if (nd + no != 1)
    {
      return false;
    }
    const double a = nd ? da[di[i]] : oa[oi[i]];
    return std::abs(a - 1.0) < 1.0e-12;
  }

  // f(local column or offd column, is_offd, weight)
  template <typename F>
  void ForRow(int i, F &&f) const
  {
    for (int k = di[i]; k < di[i + 1]; k++)
    {
      f(dj[k], false, da[k]);
    }
    for (int k = 0, n = NumOffdRow(i); k < n; k++)
    {
      f(oj[oi[i] + k], true, oa[oi[i] + k]);
    }
  }

  HYPRE_BigInt GlobalCol(int c, bool offd) const { return offd ? cmap[c] : col0 + c; }

  // Values of a true-dof quantity at the offd columns (halo exchange through the matvec
  // communication package), as two 32-bit halves.
  std::vector<Key> Halo(const std::vector<Key> &tvals) const
  {
    std::vector<Key> out(noffd, 0);
    if (!hypre_ParCSRMatrixCommPkg(h))
    {
      hypre_MatvecCommPkgCreate(h);
    }
    hypre_ParCSRCommPkg *pkg = hypre_ParCSRMatrixCommPkg(h);
    const int ns = hypre_ParCSRCommPkgNumSends(pkg);
    const int nsend = hypre_ParCSRCommPkgSendMapStart(pkg, ns);
    std::vector<HYPRE_Int> sbuf(std::max(nsend, 1)), rbuf(std::max(noffd, 1));
    for (int half = 0; half < 2; half++)
    {
      for (int k = 0; k < nsend; k++)
      {
        const Key v = tvals[hypre_ParCSRCommPkgSendMapElmt(pkg, k)];
        sbuf[k] = static_cast<HYPRE_Int>(
            static_cast<std::uint32_t>(half ? (v >> 32) : (v & 0xffffffffULL)));
      }
      hypre_ParCSRCommHandle *hd =
          hypre_ParCSRCommHandleCreate(11, pkg, sbuf.data(), rbuf.data());
      hypre_ParCSRCommHandleDestroy(hd);
      for (int k = 0; k < noffd; k++)
      {
        const Key v = static_cast<Key>(static_cast<std::uint32_t>(rbuf[k]));
        out[k] |= half ? (v << 32) : v;
      }
    }
    return out;
  }
};

// ---------------------------------------------------------------------------------------

// Element matrices of a local operator, matrix-free

// Element matrices of a libCEED sub-operator by element-local probing. An operator with the
// same QFunction, bases and passive (quadrature data) fields is built whose active fields
// live on a broken L-vector: every element has its own copy of its dofs, with the
// orientation transformations of the original element restriction. Applied to the indicator
// of local dof a of every element at once, it returns column a of every element matrix (in
// terms of the L-dofs, as an assembly would produce them). The cost is one operator
// application per element dof and the storage two vectors of (elements x element dofs): no
// element matrix is stored. f(ne, n, offsets, a, y) is called for a = 0, ..., n - 1, with
// y[e n + b] = K_e(b, a) and offsets[e n + b] the L-dof of local dof b of element e.
template <typename F>
void ProbeSubOperator(CeedOperator sub, F &&f)
{
  Ceed ceed;
  PalaceCeedCallBackend(CeedOperatorGetCeed(sub, &ceed));
  CeedQFunction qf;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(sub, &qf));
  CeedInt ni, no;
  CeedOperatorField *inputs, *outputs;
  PalaceCeedCall(ceed, CeedOperatorGetFields(sub, &ni, &inputs, &no, &outputs));

  // The active element restriction.
  CeedElemRestriction rstr = nullptr;
  for (CeedInt k = 0; k < ni && !rstr; k++)
  {
    CeedElemRestriction r;
    CeedBasis b;
    CeedVector v;
    PalaceCeedCall(ceed, CeedOperatorFieldGetData(inputs[k], nullptr, &r, &b, &v));
    if (v == CEED_VECTOR_ACTIVE)
    {
      rstr = r;
    }
    else
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
    }
    PalaceCeedCall(ceed, CeedBasisDestroy(&b));
    PalaceCeedCall(ceed, CeedVectorDestroy(&v));
  }
  MFEM_VERIFY(rstr, "Block Jacobi: sub-operator without an active field!");
  CeedInt ne, n, ncomp;
  CeedRestrictionType type;
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(rstr, &ne));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(rstr, &n));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(rstr, &ncomp));
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(rstr, &type));
  MFEM_VERIFY(ncomp == 1 &&
                  (type == CEED_RESTRICTION_STANDARD || type == CEED_RESTRICTION_ORIENTED ||
                   type == CEED_RESTRICTION_CURL_ORIENTED),
              "Block Jacobi: unsupported element restriction!");
  const CeedSize esize = static_cast<CeedSize>(ne) * n;
  if (ne == 0)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&rstr));
    PalaceCeedCall(ceed, CeedQFunctionDestroy(&qf));
    return;
  }

  // The broken restriction (identity offsets, same orientations).
  CeedElemRestriction rb;
  {
    std::vector<CeedInt> ident(esize);
    std::iota(ident.begin(), ident.end(), 0);
    if (type == CEED_RESTRICTION_CURL_ORIENTED)
    {
      const CeedInt8 *orients;
      PalaceCeedCall(ceed,
                     CeedElemRestrictionGetCurlOrientations(rstr, CEED_MEM_HOST, &orients));
      PalaceCeedCall(ceed, CeedElemRestrictionCreateCurlOriented(
                               ceed, ne, n, 1, 1, esize, CEED_MEM_HOST, CEED_COPY_VALUES,
                               ident.data(), orients, &rb));
      PalaceCeedCall(ceed, CeedElemRestrictionRestoreCurlOrientations(rstr, &orients));
    }
    else if (type == CEED_RESTRICTION_ORIENTED)
    {
      const bool *orients;
      PalaceCeedCall(ceed,
                     CeedElemRestrictionGetOrientations(rstr, CEED_MEM_HOST, &orients));
      PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                               ceed, ne, n, 1, 1, esize, CEED_MEM_HOST, CEED_COPY_VALUES,
                               ident.data(), orients, &rb));
      PalaceCeedCall(ceed, CeedElemRestrictionRestoreOrientations(rstr, &orients));
    }
    else
    {
      PalaceCeedCall(ceed,
                     CeedElemRestrictionCreate(ceed, ne, n, 1, 1, esize, CEED_MEM_HOST,
                                               CEED_COPY_VALUES, ident.data(), &rb));
    }
  }

  // The element-level operator.
  CeedOperator ope;
  PalaceCeedCall(
      ceed, CeedOperatorCreate(ceed, qf, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &ope));
  for (int io = 0; io < 2; io++)
  {
    CeedOperatorField *fields = io ? outputs : inputs;
    for (CeedInt k = 0; k < (io ? no : ni); k++)
    {
      const char *name;
      CeedElemRestriction r;
      CeedBasis b;
      CeedVector v;
      PalaceCeedCall(ceed, CeedOperatorFieldGetData(fields[k], &name, &r, &b, &v));
      if (v == CEED_VECTOR_ACTIVE)
      {
        MFEM_VERIFY(r == rstr, "Block Jacobi: active fields with different restrictions!");
        PalaceCeedCall(ceed, CeedOperatorSetField(ope, name, rb, b, CEED_VECTOR_ACTIVE));
      }
      else
      {
        PalaceCeedCall(ceed, CeedOperatorSetField(ope, name, r, b, v));
      }
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
      PalaceCeedCall(ceed, CeedBasisDestroy(&b));
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(ope));

  // Probe.
  CeedVector x, y;
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, esize, &x));
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, esize, &y));
  const CeedInt *offsets;
  PalaceCeedCall(ceed, CeedElemRestrictionGetOffsets(rstr, CEED_MEM_HOST, &offsets));
  for (CeedInt a = 0; a < n; a++)
  {
    CeedScalar *xa;
    PalaceCeedCall(ceed, CeedVectorGetArrayWrite(x, CEED_MEM_HOST, &xa));
    std::fill(xa, xa + esize, 0.0);
    for (CeedInt e = 0; e < ne; e++)
    {
      xa[static_cast<std::size_t>(e) * n + a] = 1.0;
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArray(x, &xa));
    PalaceCeedCall(ceed, CeedOperatorApply(ope, x, y, CEED_REQUEST_IMMEDIATE));
    const CeedScalar *ya;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(y, CEED_MEM_HOST, &ya));
    f(ne, n, offsets, a, ya);
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(y, &ya));
  }
  PalaceCeedCall(ceed, CeedElemRestrictionRestoreOffsets(rstr, &offsets));
  PalaceCeedCall(ceed, CeedVectorDestroy(&x));
  PalaceCeedCall(ceed, CeedVectorDestroy(&y));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&ope));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&rb));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&rstr));
  PalaceCeedCall(ceed, CeedQFunctionDestroy(&qf));
}

// The sub-operators of a libCEED operator (over threads, geometries, integrators).
template <typename F>
void ForEachSubOperator(const ceed::Operator &op, F &&f)
{
  for (std::size_t id = 0; id < op.Size(); id++)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    CeedInt nsub;
    CeedOperator *subs;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &nsub));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[id], &subs));
    for (CeedInt k = 0; k < nsub; k++)
    {
      f(subs[k]);
    }
  }
}

// f(row, col, value) for every entry of every element matrix of a local operator (libCEED:
// by element-local probing), or for every entry of an assembled local matrix.
template <typename F>
void ForEachLocalEntry(const Operator &op, F &&f)
{
  if (const auto *cop = dynamic_cast<const ceed::Operator *>(&op))
  {
    ForEachSubOperator(*cop,
                       [&](CeedOperator sub)
                       {
                         ProbeSubOperator(sub,
                                          [&](CeedInt ne, CeedInt n, const CeedInt *offsets,
                                              CeedInt a, const CeedScalar *y)
                                          {
                                            for (CeedInt e = 0; e < ne; e++)
                                            {
                                              const std::size_t base =
                                                  static_cast<std::size_t>(e) * n;
                                              const CeedInt col = offsets[base + a];
                                              for (CeedInt b = 0; b < n; b++)
                                              {
                                                f(offsets[base + b], col, y[base + b]);
                                              }
                                            }
                                          });
                       });
  }
  else if (const auto *sop = dynamic_cast<const hypre::HypreCSRMatrix *>(&op))
  {
    // A copy of the matrix in host memory, wherever the matrix is.
    hypre_CSRMatrix *host = hypre_CSRMatrixClone_v2(
        *const_cast<hypre::HypreCSRMatrix *>(sop), 1, HYPRE_MEMORY_HOST);
    const int n = hypre_CSRMatrixNumRows(host);
    const auto *I = hypre_CSRMatrixI(host);
    const auto *J = hypre_CSRMatrixJ(host);
    const auto *A = hypre_CSRMatrixData(host);
    for (int i = 0; i < n; i++)
    {
      for (int k = I[i]; k < I[i + 1]; k++)
      {
        f(i, J[k], A[k]);
      }
    }
    hypre_CSRMatrixDestroy(host);
  }
  else
  {
    MFEM_ABORT("Block Jacobi: the local operator is neither a libCEED operator nor a CSR "
               "matrix!");
  }
}

}  // namespace palace::blocks

#endif  // PALACE_LINALG_BLOCK_UTILS_HPP
