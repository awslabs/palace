// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CURL_CURL_HPP
#define PALACE_LINALG_CURL_CURL_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/complex.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MaterialOperator;
class KspSolver;

//
// This solver implements a solver for the operator K + M in a Nedelec space.
//
class CurlCurlMassSolver : public mfem::Solver
{
private:
  // H(curl) norm operator A = K + M and its projection Gᵀ A G.
  std::vector<std::unique_ptr<ParOperator>> A, AuxA;

  // Linear solver for the linear system A y = x;
  std::unique_ptr<KspSolver> ksp;

public:
  CurlCurlMassSolver(const MaterialOperator &mat_op,
                     mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
                     mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                     const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
                     const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists, double tol,
                     int max_it, int print);

  void SetOperator(const Operator &op) override {}

  void Mult(const Vector &x, Vector &y) const override { ksp->Mult(x, y); }
  void Mult(const ComplexVector &x, ComplexVector &y)
  {
    Mult(x.Real(), y.Real());
    Mult(x.Imag(), y.Imag());
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_CURL_CURL_HPP
