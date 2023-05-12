// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CURL_CURL_HPP
#define PALACE_LINALG_CURL_CURL_HPP

#include <memory>
#include <mfem.hpp>
#include "linalg/petsc.hpp"

namespace palace
{

class MaterialOperator;

//
// This solver implements a solver for the operator K + M in the Nedelec space.
//
class CurlCurlMassSolver : public mfem::Solver
{
private:
  // H(curl) norm operator A = K + M.
  std::vector<std::unique_ptr<mfem::Operator>> A, AuxA;

  // Linear solver and preconditioner for the linear system A y = x;
  std::unique_ptr<mfem::IterativeSolver> ksp;
  std::unique_ptr<mfem::Solver> pc;

  // Workspace objects for solver application.
  mutable mfem::Vector xr, xi, yr, yi;

public:
  CurlCurlMassSolver(const MaterialOperator &mat_op, const mfem::Array<int> &dbc_marker,
                     mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
                     mfem::ParFiniteElementSpaceHierarchy &h1_fespaces, double tol,
                     int max_it, int print);

  void SetOperator(const mfem::Operator &op) override {}

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override { ksp->Mult(x, y); }
  void Mult(const petsc::PetscParVector &x, petsc::PetscParVector &y) const
  {
    x.GetToVectors(xr, xi);
    Mult(xr, yr);
    Mult(xi, yi);
    y.SetFromVectors(yr, yi);
  }
  using mfem::Operator::Mult;
};

}  // namespace palace

#endif  // PALACE_LINALG_CURL_CURL_HPP
