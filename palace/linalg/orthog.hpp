// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ORTHOG_HPP
#define PALACE_LINALG_ORTHOG_HPP

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

namespace
{

inline void MultHelper(const Operator &M, const Vector &v_in, Vector &v_out)
{
  M.Mult(v_in, v_out);
}

inline void MultHelper(const Operator &M, const ComplexVector &v_in, ComplexVector &v_out)
{
  M.Mult(v_in.Real(), v_out.Real());
  M.Mult(v_in.Imag(), v_out.Imag());
}

}  // namespace

// Orthogonalization functions with weigh matrix W assumed to be SPD. Make these separate
// functions due overhead of allocation of intermediate vector and different inner product
// call structure.
//
// We need a temporary vector `WVj` to write in the matrix-vector multiplication as the
// inner-product function on W is not available. Assume this temporary vector to be passed
// in and correct in shape, device properties, etc just like in `linalg::Norml2`.
template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnMGSWeighted(MPI_Comm comm, const Operator &weight_matrix_W,
                                           const std::vector<VecType> &V, VecType &w,
                                           ScalarType *H, int m, VecType &WVj)
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for CGS orthogonalization!");
  for (int j = 0; j < m; j++)
  {
    MultHelper(weight_matrix_W, V[j], WVj);
    H[j] = linalg::Dot(comm, w, WVj);  // Global inner product
    w.Add(-H[j], V[j]);
  }
}

template <typename VecType, typename ScalarType>
inline void OrthogonalizeColumnCGSWeighted(MPI_Comm comm, const Operator &weight_matrix_W,
                                           const std::vector<VecType> &V, VecType &w,
                                           ScalarType *H, int m, VecType &WVj,
                                           bool refine = false)
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for CGS orthogonalization!");
  for (int j = 0; j < m; j++)
  {
    MultHelper(weight_matrix_W, V[j], WVj);
    H[j] = w * WVj;  // Local inner product
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
      MultHelper(weight_matrix_W, V[j], WVj);
      dH[j] = w * WVj;  // Local inner product
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
