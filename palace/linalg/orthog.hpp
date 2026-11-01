// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>
#include "linalg/fusedvector.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/labels.hpp"

namespace palace::linalg
{

//
// Orthogonalization functions for orthogonalizing a vector against a number of basis
// vectors using modified or classical Gram-Schmidt.
//
// Assumes that the input vectors are normalized, but does not normalize the output vectors!
// If done in a loop, normalization has to be managed by hand! (TODO: Reconsider).
//

// Concept: InnerProductHelper has function InnerProduct(VecType, VecType) -> ScalarType,
// acting on local degrees of freedom. Also add MPI reduction.

// Simplest case is canonical inner product on R & C.
struct IdentityInnerProduct
{
  template <typename VecType>
  auto operator()(const VecType &x, const VecType &y) const
  {
    return LocalDot(x, y);
  }
};

namespace internal
{

// Modified Gram-Schmidt with each column update fused with the inner product for the next
// column (or the squared norm of the result, if norm2 is given) on the host: same
// arithmetic and reductions as the unfused loop, one fewer pass over w per column.
template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnMGSFused(MPI_Comm comm, const std::vector<VecType> &V,
                                        VecType &w, ScalarType *H, std::size_t m,
                                        ScalarType *norm2)
{
  H[0] = LocalDot(w, V[0]);
  Mpi::GlobalSum(1, &H[0], comm);
  for (std::size_t j = 0; j + 1 < m; j++)
  {
    H[j + 1] = fused::AddDot(-H[j], V[j], w, &V[j + 1]);
    Mpi::GlobalSum(1, &H[j + 1], comm);
  }
  if (norm2)
  {
    *norm2 = fused::AddDot(-H[m - 1], V[m - 1], w, nullptr);
    Mpi::GlobalSum(1, norm2, comm);
  }
  else
  {
    w.Add(-H[m - 1], V[m - 1]);
  }
}

}  // namespace internal

template <typename VecType, typename ScalarType,
          typename InnerProductW = IdentityInnerProduct>
inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m,
                                   const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(m <= V.size(), "Out of bounds number of columns for MGS orthogonalization!");
  if constexpr (std::is_same_v<InnerProductW, IdentityInnerProduct>)
  {
    if (m > 0 && fused::UseHost())
    {
      internal::OrthogonalizeColumnMGSFused(comm, V, w, H, m, (ScalarType *)nullptr);
      return;
    }
  }
  for (std::size_t j = 0; j < m; j++)
  {
    // Global inner product: Note order is important for complex vectors.
    H[j] = dot_op(w, V[j]);
    Mpi::GlobalSum(1, &H[j], comm);
    w.Add(-H[j], V[j]);
  }
}

namespace internal
{

template <typename VecType, typename ScalarType, typename Project>
inline void OrthogonalizeColumnCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m, bool refine,
                                   const Project &project)
{
  MFEM_ASSERT(m <= V.size(), "Out of bounds number of columns for CGS orthogonalization!");
  if (m == 0)
  {
    return;
  }
  project(w, H);
  Mpi::GlobalSum(m, H, comm);
  for (std::size_t j = 0; j < m; j++)
  {
    w.Add(-H[j], V[j]);
  }
  if (refine)
  {
    std::vector<ScalarType> dH(m);
    project(w, dH.data());
    Mpi::GlobalSum(m, dH.data(), comm);
    for (std::size_t j = 0; j < m; j++)
    {
      H[j] += dH[j];
      w.Add(-dH[j], V[j]);
    }
  }
}

}  // namespace internal

template <typename VecType, typename ScalarType,
          typename InnerProductW = IdentityInnerProduct>
inline void OrthogonalizeColumnCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m, bool refine = false,
                                   const InnerProductW &dot_op = {})
{
  internal::OrthogonalizeColumnCGS(comm, V, w, H, m, refine,
                                   [&V, m, &dot_op](const VecType &x, ScalarType *h)
                                   {
                                     for (std::size_t j = 0; j < m; j++)
                                     {
                                       h[j] = dot_op(x, V[j]);  // Local inner product
                                     }
                                   });
}

// Operator-weighted CGS needs only one application of W per pass because w is unchanged
// while computing its projections. Reapply W after the first update when refining.
// The caller owns the workspace, which must not alias w or any basis vector.
template <typename VecType, typename ScalarType, typename WeightOperator>
inline void OrthogonalizeColumnWeightedCGS(MPI_Comm comm, const std::vector<VecType> &V,
                                           VecType &w, ScalarType *H, std::size_t m,
                                           const WeightOperator &W, VecType &work,
                                           bool refine = false)
{
  internal::OrthogonalizeColumnCGS(comm, V, w, H, m, refine,
                                   [&V, m, &W, &work](const VecType &x, ScalarType *h)
                                   {
                                     work.SetSize(x.Size());
                                     work.UseDevice(x.UseDevice());
                                     W.Mult(x, work);
                                     for (std::size_t j = 0; j < m; j++)
                                     {
                                       h[j] = LocalDot(work, V[j]);
                                     }
                                   });
}

template <typename VecType, typename ScalarType,
          typename InnerProductW = IdentityInnerProduct>
inline void OrthogonalizeColumn(Orthogonalization type, MPI_Comm comm,
                                const std::vector<VecType> &V, VecType &w, ScalarType *H,
                                std::size_t m, const InnerProductW &dot_op = {})
{
  switch (type)
  {
    case Orthogonalization::MGS:
      OrthogonalizeColumnMGS(comm, V, w, H, m, dot_op);
      break;
    case Orthogonalization::CGS:
      OrthogonalizeColumnCGS(comm, V, w, H, m, false, dot_op);
      break;
    case Orthogonalization::CGS2:
      OrthogonalizeColumnCGS(comm, V, w, H, m, true, dot_op);
      break;
  }
}

// Orthogonalize w against the columns of V as OrthogonalizeColumn and return the 2-norm of
// the result, as Norml2.
template <typename VecType, typename ScalarType>
inline double OrthogonalizeColumnAndNorm(Orthogonalization type, MPI_Comm comm,
                                         const std::vector<VecType> &V, VecType &w,
                                         ScalarType *H, std::size_t m)
{
  if (type == Orthogonalization::MGS && m > 0 && fused::UseHost())
  {
    ScalarType norm2;
    internal::OrthogonalizeColumnMGSFused(comm, V, w, H, m, &norm2);
    return std::sqrt(std::abs(norm2));
  }
  OrthogonalizeColumn(type, comm, V, w, H, m);
  return Norml2(comm, w);
}

}  // namespace palace::linalg

#endif  // PALACE_LINALG_ORTHOG_HPP
