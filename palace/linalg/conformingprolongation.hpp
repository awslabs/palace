// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CONFORMING_PROLONGATION_HPP
#define PALACE_LINALG_CONFORMING_PROLONGATION_HPP

#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

// Host products with the prolongation of a conforming parallel finite element space
// (mfem::ConformingProlongationOperator), with an additional fused product for the real and
// imaginary parts of a complex vector: a single exchange of interleaved values instead of
// one per part. Bitwise identical to two separate products with the MFEM operator.
class ConformingProlongation : public Operator
{
private:
  const mfem::ConformingProlongationOperator &P;

  // Group communicator for interleaved (real, imaginary) L-vector data, with the same
  // groups and the same (consistent across processes) order of shared L-DOFs as the
  // communicator of the space.
  mfem::GroupCommunicator gc2;

  // Sorted L-DOFs owned by other processes.
  mfem::Array<int> external_ldofs;

  mutable std::vector<double> buf;

public:
  ConformingProlongation(const mfem::ConformingProlongationOperator &P);

  const mfem::ConformingProlongationOperator &GetOperator() const { return P; }

  void Mult(const Vector &x, Vector &y) const override { P.Mult(x, y); }

  void MultTranspose(const Vector &x, Vector &y) const override { P.MultTranspose(x, y); }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    P.AddMult(x, y, a);
  }

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    P.AddMultTranspose(x, y, a);
  }

  // y = P x and y = Pᵀ x for both parts (r, i) of a complex vector.
  void Mult(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi) const;
  void MultTranspose(const Vector &xr, const Vector &xi, Vector &yr, Vector &yi) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_CONFORMING_PROLONGATION_HPP
