// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "entityjacobi.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <type_traits>
#include <ceed/backend.h>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/operator.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// In-place inverse of a symmetric positive definite matrix (row-major) from its Cholesky
// factorization A = L Lᵀ, A⁻¹ = L⁻ᵀ L⁻¹. If A is not numerically positive definite, A is
// replaced by the inverse of its diagonal and false is returned.
bool InvertSPD(double *A, int n, std::vector<double> &L)
{
  L.assign(std::size_t(n) * n, 0.0);
  for (int j = 0; j < n; j++)
  {
    double s = A[j * n + j];
    for (int k = 0; k < j; k++)
    {
      s -= L[j * n + k] * L[j * n + k];
    }
    if (!(s > 0.0))
    {
      for (int i = 0; i < n; i++)
      {
        for (int k = 0; k < n; k++)
        {
          A[i * n + k] = (i == k) ? 1.0 / A[i * n + i] : 0.0;
        }
      }
      return false;
    }
    L[j * n + j] = std::sqrt(s);
    for (int i = j + 1; i < n; i++)
    {
      double t = A[i * n + j];
      for (int k = 0; k < j; k++)
      {
        t -= L[i * n + k] * L[j * n + k];
      }
      L[i * n + j] = t / L[j * n + j];
    }
  }
  // Overwrite L with L⁻¹ column by column (lower triangular).
  for (int j = 0; j < n; j++)
  {
    L[j * n + j] = 1.0 / L[j * n + j];
    for (int i = j + 1; i < n; i++)
    {
      double t = 0.0;
      for (int k = j; k < i; k++)
      {
        t -= L[i * n + k] * L[k * n + j];
      }
      L[i * n + j] = t / L[i * n + i];
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      double t = 0.0;
      for (int k = i; k < n; k++)
      {
        t += L[k * n + i] * L[k * n + j];
      }
      A[i * n + j] = A[j * n + i] = t;
    }
  }
  return true;
}

// Mesh entity of each L-DOF of the space whose prolongation row (in Pd and Po) is a single
// ±1 entry (unconstrained), numbering the edges, then the faces (in 3D), then the element
// interiors, or -1. Returns the number of entities.
long long GetLDofEntities(const FiniteElementSpace &fespace, const mfem::SparseMatrix &Pd,
                          const mfem::SparseMatrix &Po, std::vector<long long> &entity)
{
  const auto &ufes = fespace.Get();  // Underlying MFEM space, also for a broken space
  const auto &mesh = fespace.GetMesh().Get();
  const long long nedge = mesh.GetNEdges(),
                  nface = (mesh.Dimension() == 3) ? mesh.GetNFaces() : 0;
  entity.assign(fespace.GetVSize(), -1);
  mfem::Array<int> dofs;
  auto Set = [&](long long id)
  {
    for (int d : dofs)
    {
      entity[(d >= 0) ? d : -1 - d] = id;
    }
  };
  for (int i = 0; i < nedge; i++)
  {
    ufes.GetEdgeInteriorDofs(i, dofs);
    Set(i);
  }
  for (int i = 0; i < nface; i++)
  {
    ufes.GetFaceInteriorDofs(i, dofs);
    Set(nedge + i);
  }
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    ufes.GetElementInteriorDofs(i, dofs);
    Set(nedge + nface + i);
  }
  for (int l = 0; l < fespace.GetVSize(); l++)
  {
    const int nd = Pd.RowSize(l), no = Po.RowSize(l);
    const double v = nd ? Pd.GetRowEntries(l)[0] : (no ? Po.GetRowEntries(l)[0] : 0.0);
    if (nd + no != 1 || std::abs(std::abs(v) - 1.0) >= 1.0e-12)
    {
      entity[l] = -1;
    }
  }
  return nedge + nface + mesh.GetNE();
}

// Calls f(p, dofs, A_e) for the element matrix A_e (p x p, row-major) of each element of
// each libCEED sub-operator (one at a time), with the L-DOFs dofs of its rows and columns.
template <typename Function>
void ForEachElementMatrix(const ceed::Operator &op, Function &&f)
{
  for (std::size_t id = 0; id < op.Size(); id++)
  {
    Ceed ceed;
    CeedInt nsub;
    CeedOperator *subs;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &nsub));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[id], &subs));
    for (CeedInt s = 0; s < nsub; s++)
    {
      CeedElemRestriction restr;
      CeedInt ne, p, ncomp;
      PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestriction(subs[s], &restr));
      PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(restr, &ne));
      PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(restr, &p));
      PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(restr, &ncomp));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&restr));
      MFEM_VERIFY(ncomp == 1, "EntityBlockJacobiSolver requires a single-component space!");
      if (ne == 0)
      {
        continue;
      }
      CeedSize nnz;
      CeedInt *rows, *cols;
      CeedVector vals;
      PalaceCeedCall(ceed, CeedOperatorLinearAssembleSymbolic(subs[s], &nnz, &rows, &cols));
      MFEM_VERIFY(nnz == CeedSize(ne) * p * p,
                  "Unexpected element matrix layout for EntityBlockJacobiSolver!");
      std::vector<int> edofs(std::size_t(ne) * p);  // L-DOF of each element matrix row
      for (std::size_t k = 0; k < edofs.size(); k++)
      {
        edofs[k] = rows[k * p];
      }
      std::free(rows);
      std::free(cols);
      PalaceCeedCall(ceed, CeedVectorCreate(ceed, nnz, &vals));
      PalaceCeedCall(ceed, CeedOperatorLinearAssemble(subs[s], vals));
      PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(subs[s]));
      const CeedScalar *data;
      PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals, CEED_MEM_HOST, &data));
      for (CeedInt e = 0; e < ne; e++)
      {
        f(p, edofs.data() + std::size_t(e) * p, data + std::size_t(e) * p * p);
      }
      PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals, &data));
      PalaceCeedCall(ceed, CeedVectorDestroy(&vals));
    }
  }
}

}  // namespace

template <OperatorType OperType>
EntityBlockJacobiSolver<OperType>::EntityBlockJacobiSolver(
    const ceed::Operator &op, const FiniteElementSpace &fespace)
  : Solver<OperType>()
{
  const mfem::HypreParMatrix *P = fespace.Dof_TrueDof_Matrix();
  P->HostRead();
  const int nl = fespace.GetVSize(), nt = fespace.GetTrueVSize();
  mfem::SparseMatrix Pd, Po;
  HYPRE_BigInt *cmap;
  P->GetDiag(Pd);
  P->GetOffd(Po, cmap);

  // Block of each owned true DOF: the entity of its unconstrained local L-DOF, with the
  // copies of a broken space in their own blocks, or a block of its own.
  std::vector<long long> entity, key(nt);
  const long long nent = GetLDofEntities(fespace, Pd, Po, entity);
  for (int t = 0; t < nt; t++)
  {
    key[t] = -1 - t;
  }
  std::vector<int> ltdof(nl, -1);  // Owned true DOF of each unconstrained L-DOF
  for (int l = 0; l < nl; l++)
  {
    if (entity[l] >= 0 && Pd.RowSize(l) == 1)
    {
      ltdof[l] = Pd.GetRowColumns(l)[0];
      key[ltdof[l]] = entity[l];
    }
  }
  if (const auto *bP = fespace.GetBrokenProlongation())
  {
    for (int t = 0; t < fespace.Get().GetTrueVSize(); t++)
    {
      for (int v = 1; key[t] >= 0 && v < bP->num_versions[t]; v++)
      {
        key[bP->GetVersionTrueDof(t, v)] = key[t] + v * nent;
      }
    }
  }
  block_dofs.resize(nt);
  std::iota(block_dofs.begin(), block_dofs.end(), 0);
  std::stable_sort(block_dofs.begin(), block_dofs.end(),
                   [&key](int s, int t) { return key[s] < key[t]; });
  std::vector<int> block(nt), pos(nt);
  block_ptr.assign(1, 0);
  for (int k = 0; k < nt; k++)
  {
    const int t = block_dofs[k];
    if (k > 0 && key[t] != key[block_dofs[k - 1]])
    {
      block_ptr.push_back(k);
    }
    block[t] = static_cast<int>(block_ptr.size()) - 1;
    pos[t] = k - block_ptr.back();
  }
  block_ptr.push_back(nt);
  const int nb = static_cast<int>(block_ptr.size()) - 1;
  inv_ptr.assign(nb + 1, 0);
  for (int b = 0; b < nb; b++)
  {
    const std::size_t n = block_ptr[b + 1] - block_ptr[b];
    inv_ptr[b + 1] = inv_ptr[b] + n * n;
  }
  inv.assign(inv_ptr[nb], 0.0);

  // The blocks of Pᵀ A_loc P from the element matrices: the entries between unconstrained
  // L-DOFs of owned true DOFs of the same entity are added to them directly, and the other
  // entries which can contribute to a block (with an L-DOF of no entity, or between L-DOFs
  // of the same entity) are assembled as by ParallelAssemble, which adds the contributions
  // of the elements of other processes.
  {
    mfem::SparseMatrix A(nl);
    mfem::DenseMatrix Ae_f;
    auto AddElementMatrix = [&](int p, int *de, const CeedScalar *Ae)
    {
      Ae_f.SetSize(p);
      Ae_f = 0.0;
      bool rest = false;
      for (int i = 0; i < p; i++)
      {
        for (int j = 0; j < p; j++)
        {
          const int ti = ltdof[de[i]], tj = ltdof[de[j]];
          const long long ei = entity[de[i]], ej = entity[de[j]];
          const double a = 0.5 * (Ae[i * p + j] + Ae[j * p + i]);
          if (ti >= 0 && tj >= 0)
          {
            if (ei == ej)
            {
              const int b = block[ti], n = block_ptr[b + 1] - block_ptr[b];
              inv[inv_ptr[b] + std::size_t(pos[ti]) * n + pos[tj]] +=
                  Pd.GetRowEntries(de[i])[0] * a * Pd.GetRowEntries(de[j])[0];
            }
          }
          else if (ei < 0 || ej < 0 || ei == ej)
          {
            Ae_f(i, j) = a;
            rest = true;
          }
        }
      }
      if (rest)
      {
        const mfem::Array<int> dofs(de, p);
        A.AddSubMatrix(dofs, dofs, Ae_f);
      }
    };
    ForEachElementMatrix(op, AddElementMatrix);
    A.Finalize();
    hypre::HypreCSRMatrix hA(A);
    auto PtAP = ParOperator(hA, fespace).StealParallelAssemble();
    mfem::SparseMatrix S;
    PtAP->GetDiag(S);
    for (int t = 0; t < nt; t++)
    {
      const int b = block[t];
      const std::size_t n = block_ptr[b + 1] - block_ptr[b];
      for (int k = S.GetI()[t]; k < S.GetI()[t + 1]; k++)
      {
        if (block[S.GetJ()[k]] == b)
        {
          inv[inv_ptr[b] + pos[t] * n + pos[S.GetJ()[k]]] += S.GetData()[k];
        }
      }
    }
  }

  // Invert the blocks, using the diagonal for those which are not numerically positive
  // definite.
  std::vector<double> L;
  long long num_diag = 0;
  for (int b = 0; b < nb; b++)
  {
    num_diag += !InvertSPD(inv.data() + inv_ptr[b], block_ptr[b + 1] - block_ptr[b], L);
  }
  Mpi::GlobalSum(1, &num_diag, fespace.GetComm());
  if (num_diag > 0)
  {
    Mpi::Print(fespace.GetComm(),
               " Entity block Jacobi preconditioner: {:d} blocks are not numerically "
               "positive definite, using their diagonal\n",
               num_diag);
  }

  // The applications read the upper triangles of the inverse blocks, which halves the
  // memory traffic of the preconditioner.
  {
    int max_n = 0;
    for (int b = 0; b < nb; b++)
    {
      max_n = std::max(max_n, block_ptr[b + 1] - block_ptr[b]);
    }
    upper = max_n <= 64;
  }
  if (upper)
  {
    // The upper triangles of the (symmetrized) inverse blocks.
    sym_ptr.assign(nb + 1, 0);
    for (int b = 0; b < nb; b++)
    {
      const std::size_t n = block_ptr[b + 1] - block_ptr[b];
      sym_ptr[b + 1] = sym_ptr[b] + n * (n + 1) / 2;
    }
    std::vector<double> packed(sym_ptr[nb]);
    for (int b = 0; b < nb; b++)
    {
      const std::size_t n = block_ptr[b + 1] - block_ptr[b];
      const double *S = inv.data() + inv_ptr[b];
      double *U = packed.data() + sym_ptr[b];
      for (std::size_t i = 0; i < n; i++)
      {
        for (std::size_t j = i; j < n; j++)
        {
          *U++ = 0.5 * (S[i * n + j] + S[j * n + i]);
        }
      }
    }
    inv = std::move(packed);
  }
}

template <OperatorType OperType>
void EntityBlockJacobiSolver<OperType>::SetOperator(const OperType &op)
{
  this->height = op.Height();
  this->width = op.Width();
  VecType d(op.Height());
  d.UseDevice(true);
  op.AssembleDiagonal(d);
  d.Reciprocal();
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    dinv = d.Real();
  }
  else
  {
    dinv = d;
  }
}

template <OperatorType OperType>
bool EntityBlockJacobiSolver<OperType>::Supported(const FiniteElementSpace &fespace)
{
  // The blocks are assembled with Dof_TrueDof_Matrix(), which is what the prolongation of
  // the space applies on the host: an IdentityOperator or ConformingProlongationOperator
  // of MFEM, Dof_TrueDof_Matrix() itself, or the prolongation of a broken space.
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK | mfem::Backend::OMP_MASK))
  {
    return false;
  }
  // At first order, every entity has a single DOF. On hexahedra and quadrilaterals, Jacobi
  // is already close to optimal for the tensor-product bases, and the blocks increase the
  // number of iterations. The blocks become large on prisms from fourth order (78 DOFs) and
  // on simplices from fifth order (60 DOFs). The largest blocks used, the element interiors
  // of third-order prisms and fourth-order tetrahedra, have 30 DOFs.
  const auto map_type = fespace.GetFEColl().GetMapType(fespace.Dimension());
  const int order = fespace.GetMaxElementOrder();
  const int simplex = 1, prism = 4;  // Element types of mfem::Mesh::MeshGenerator()
  const int types = fespace.GetMesh().Get().MeshGenerator();
  return (map_type == mfem::FiniteElement::H_DIV ||
          map_type == mfem::FiniteElement::H_CURL) &&
         order >= 2 && !(types & ~(simplex | prism)) &&
         (order <= 3 || (order == 4 && types == simplex));
}

namespace
{

// y = S x by blocks for NV real vectors at once (the real and imaginary parts of a complex
// vector), so that the inverse blocks are read once per application.
template <int NV>
void ApplyEntityBlocks(const std::vector<int> &block_ptr,
                       const std::vector<int> &block_dofs,
                       const std::vector<std::size_t> &inv_ptr,
                       const std::vector<double> &inv,
                       const std::array<const double *, NV> &X,
                       const std::array<double *, NV> &Y)
{
  for (std::size_t b = 0; b + 1 < block_ptr.size(); b++)
  {
    const int n = block_ptr[b + 1] - block_ptr[b], *db = block_dofs.data() + block_ptr[b];
    const double *S = inv.data() + inv_ptr[b];
    for (int i = 0; i < n; i++)
    {
      double t[NV] = {};
      for (int j = 0; j < n; j++)
      {
        const double s = S[i * n + j];
        for (int v = 0; v < NV; v++)
        {
          t[v] += s * X[v][db[j]];
        }
      }
      for (int v = 0; v < NV; v++)
      {
        Y[v][db[i]] = t[v];
      }
    }
  }
}

// y = S x as ApplyEntityBlocks, with the upper triangle of each symmetric block S (row i
// holding the columns j >= i): the entries above the diagonal contribute to both of their
// rows.
template <int NV>
void ApplyUpperEntityBlocks(const std::vector<int> &block_ptr,
                            const std::vector<int> &block_dofs,
                            const std::vector<std::size_t> &sym_ptr, const double *inv,
                            const std::array<const double *, NV> &X,
                            const std::array<double *, NV> &Y)
{
  constexpr int max_n = 64;
  double x[max_n][NV], t[max_n][NV];
  for (std::size_t b = 0; b + 1 < block_ptr.size(); b++)
  {
    const int n = block_ptr[b + 1] - block_ptr[b], *db = block_dofs.data() + block_ptr[b];
    const double *U = inv + sym_ptr[b];
    for (int j = 0; j < n; j++)
    {
      for (int v = 0; v < NV; v++)
      {
        x[j][v] = X[v][db[j]];
        t[j][v] = 0.0;
      }
    }
    for (int i = 0; i < n; i++)
    {
      const double d = *U++;
      double s[NV];
      for (int v = 0; v < NV; v++)
      {
        s[v] = t[i][v] + d * x[i][v];
      }
      for (int j = i + 1; j < n; j++)
      {
        const double a = *U++;
        for (int v = 0; v < NV; v++)
        {
          s[v] += a * x[j][v];
          t[j][v] += a * x[i][v];
        }
      }
      for (int v = 0; v < NV; v++)
      {
        Y[v][db[i]] = s[v];
      }
    }
  }
}

}  // namespace

template <OperatorType OperType>
template <int NV>
void EntityBlockJacobiSolver<OperType>::Apply(const std::array<const double *, NV> &X,
                                              const std::array<double *, NV> &Y) const
{
  if (!upper)
  {
    ApplyEntityBlocks<NV>(block_ptr, block_dofs, inv_ptr, inv, X, Y);
  }
  else
  {
    ApplyUpperEntityBlocks<NV>(block_ptr, block_dofs, sym_ptr, inv.data(), X, Y);
  }
}

template <OperatorType OperType>
void EntityBlockJacobiSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    Apply<2>({x.Real().HostRead(), x.Imag().HostRead()},
             {y.Real().HostWrite(), y.Imag().HostWrite()});
  }
  else
  {
    Apply<1>({x.HostRead()}, {y.HostWrite()});
  }
}

template class EntityBlockJacobiSolver<Operator>;
template class EntityBlockJacobiSolver<ComplexOperator>;

}  // namespace palace
