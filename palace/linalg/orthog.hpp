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
struct InnerProduct
{
  template <typename VecType>
  auto operator()(const VecType &x, const VecType &y) const
  {
    return LocalDot(x, y);
  }
};

class RealWeightedInnerProduct
{
  // Choose generic operator, although can improve by refining for specialized type.
  std::shared_ptr<Operator> weight_op;

  // Don't have inner product wrapper / implemented in Operator, so need to allocate a
  // vector as a workspace. (TODO: Optimize this away).
  mutable Vector v_workspace = {};

  void SetWorkspace(const Vector &blueprint) const
  {
    v_workspace.SetSize(blueprint.Size());
    v_workspace.UseDevice(blueprint.UseDevice());
  }

public:
  template <typename OpType>
  explicit RealWeightedInnerProduct(const std::shared_ptr<OpType> &weight_op_)
    : weight_op(weight_op_)
  {
    MFEM_VERIFY(weight_op->Height() == weight_op->Width(),
                "Real weight operator must be square! ("
                    << weight_op->Height() << " != " << weight_op->Width());
  }
  // Follow same conventions as Dot:  yᴴ x or yᵀ x (note y comes second in the arguments).
  double operator()(const Vector &x, const Vector &y) const
  {
    SetWorkspace(x);
    weight_op->Mult(x, v_workspace);
    return LocalDot(v_workspace, y);
  }

  std::complex<double> operator()(const ComplexVector &x, const ComplexVector &y) const
  {
    using namespace std::complex_literals;
    SetWorkspace(x.Real());

    // weight_op is real.
    weight_op->Mult(x.Real(), v_workspace);
    std::complex<double> dot = {+LocalDot(v_workspace, y.Real()),
                                -LocalDot(v_workspace, y.Imag())};

    weight_op->Mult(x.Imag(), v_workspace);
    dot += std::complex<double>{LocalDot(v_workspace, y.Imag()),
                                LocalDot(v_workspace, y.Real())};

    return dot;
  }
};

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProduct>
inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, std::size_t m,
                                   const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(m <= V.size(), "Out of bounds number of columns for MGS orthogonalization!");
  for (std::size_t j = 0; j < m; j++)
  {
    // Global inner product: Note order is important for complex vectors.
    H[j] = dot_op(w, V[j]);
    Mpi::GlobalSum(1, &H[j], comm);
    w.Add(-H[j], V[j]);
  }
}

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProduct>
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
    H[j] = dot_op(w, V[j]);  // Local inner product
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
      dH[j] = dot_op(w, V[j]);  // Local inner product
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
