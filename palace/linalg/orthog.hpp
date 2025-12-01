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

// Concept: InnerProductHelper has function InnerProduct(VecType, VecType) -> ScalarType,
// acting on local degrees of freedom. Also add MPI reduction.

// Simplest case is canonical inner product on R & C.

class InnerProductStandard
{
public:
  double InnerProduct(const Vector &x, const Vector &y) const { return LocalDot(x, y); }

  double InnerProduct(MPI_Comm comm, const Vector &x, const Vector &y) const
  {
    return Dot(comm, x, y);
  }

  std::complex<double> InnerProduct(const ComplexVector &x, const ComplexVector &y) const
  {
    return LocalDot(x, y);
  }

  std::complex<double> InnerProduct(MPI_Comm comm, const ComplexVector &x,
                                    const ComplexVector &y) const
  {
    return Dot(comm, x, y);
  }
};

class InnerProductRealWeight
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
  explicit InnerProductRealWeight(const std::shared_ptr<OpType> &weight_op_)
    : weight_op(weight_op_)
  {
  }
  // Follow same conventions as Dot:  yᴴ x or yᵀ x (not y comes second in the arguments).
  double InnerProduct(const Vector &x, const Vector &y) const
  {
    SetWorkspace(x);
    weight_op->Mult(x, v_workspace);
    return LocalDot(v_workspace, y);
  }

  double InnerProduct(MPI_Comm comm, const Vector &x, const Vector &y) const
  {
    SetWorkspace(x);
    weight_op->Mult(x, v_workspace);
    return Dot(comm, v_workspace, y);
  }

  std::complex<double> InnerProduct(const ComplexVector &x, const ComplexVector &y) const
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

  std::complex<double> InnerProduct(MPI_Comm comm, const ComplexVector &x,
                                    const ComplexVector &y) const
  {
    auto dot = InnerProduct(x, y);
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  }
};

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProductStandard>
inline void OrthogonalizeColumnMGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, int m, const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for MGS orthogonalization!");
  for (int j = 0; j < m; j++)
  {
    // Global inner product: Note order is important for complex vectors.
    H[j] = dot_op.InnerProduct(comm, w, V[j]);
    w.Add(-H[j], V[j]);
  }
}

template <typename VecType, typename ScalarType,
          typename InnerProductW = InnerProductStandard>
inline void OrthogonalizeColumnCGS(MPI_Comm comm, const std::vector<VecType> &V, VecType &w,
                                   ScalarType *H, int m, bool refine = false,
                                   const InnerProductW &dot_op = {})
{
  MFEM_ASSERT(static_cast<std::size_t>(m) <= V.size(),
              "Out of bounds number of columns for CGS orthogonalization!");
  if (m == 0)
  {
    return;
  }
  for (int j = 0; j < m; j++)
  {
    H[j] = dot_op.InnerProduct(w, V[j]);  // Local inner product
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
      dH[j] = dot_op.InnerProduct(w, V[j]);  // Local inner product
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
