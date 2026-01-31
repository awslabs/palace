// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

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

class InnerProductStandard
{
public:
  template <typename VecType>
  static auto InnerProduct(const VecType &x, const VecType &y)
  {
    return LocalDot(x, y);
  }

  template <typename VecType>
  static auto InnerProduct(MPI_Comm comm, const VecType &x, const VecType &y)
  {
    return Dot(comm, x, y);
  }
};

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProductStandard>
inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m,
                                   const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(m <= V.size(), "Out of bounds number of columns for MGS orthogonalization!");
  for (std::size_t j = 0; j < m; j++)
  {
    // Global inner product: Note order is important for complex vectors.
    H[j] = dot_op.InnerProduct(comm, w, V[j]);
    w.Add(-H[j], V[j]);
  }
}

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProductStandard>
inline void OrthogonalizeColumnCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m, bool refine = false,
                                   const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(m <= V.size(), "Out of bounds number of columns for CGS orthogonalization!");
  if (m == 0)
  {
    return;
  }
  for (std::size_t j = 0; j < m; j++)
  {
    H[j] = dot_op.InnerProduct(w, V[j]);  // Local inner product
  }
  Mpi::GlobalSum(m, H, comm);
  for (std::size_t j = 0; j < m; j++)
  {
    w.Add(-H[j], V[j]);
  }
  if (refine)
  {
    std::vector<ScalarType> dH(m);
    for (int j = 0; j < m; j++)
    {
      dH[j] = dot_op.InnerProduct(w, V[j]);  // Local inner product
    }
    Mpi::GlobalSum(m, dH.data(), comm);
    for (std::size_t j = 0; j < m; j++)
    {
      H[j] += dH[j];
      w.Add(-dH[j], V[j]);
    }
  }
}

}  // namespace palace::linalg

#endif  // PALACE_LINALG_ORTHOG_HPP
