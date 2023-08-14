// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

#include <vector>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace::linalg
{

//
// Orthogonalization functions for orthogonalizing a vector against a number of basis
// vectors using modified or classical Gram-Schmidt.
//

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, int m)
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for MGS orthogonalization!");
  for (int j = 0; j < m; j++)
  {
    H[j] = linalg::Dot(comm, w, V[j]);  // Global inner product
    w.Add(-H[j], V[j]);
  }
}

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, int m, bool refine = false)
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for CGS orthogonalization!");
  if (m == 0)
  {
    return;
  }
  for (int j = 0; j < m; j++)
  {
    H[j] = w * V[j];  // Local inner product
  }
  Mpi::GlobalSum(m, H, comm);
  for (int j = 0; j < m; j++)
  {
    w.Add(-H[j], V[j]);
  }
  if (refine)
  {
    std::vector<ScalarType> dH(m);
    for (int j = 0; j < m; j++)
    {
      dH[j] = w * V[j];  // Local inner product
    }
    Mpi::GlobalSum(m, dH.data(), comm);
    for (int j = 0; j < m; j++)
    {
      H[j] += dH[j];
      w.Add(-dH[j], V[j]);
    }
  }
}

}  // namespace palace::linalg

#endif  // PALACE_LINALG_ORTHOG_HPP
