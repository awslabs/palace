// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
#define PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP

#include <memory>
#include <mfem.hpp>
#include "drivers/basesolver.hpp"
#include "fem/fespace.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

//
// Lightweight operator wrapper for the boundary mode solver that provides the interface
// expected by PostOperator. Holds references to the FE spaces and material operator created
// by the solver. Also owns the Btt mass matrix and L2 curl space for B-field.
//
class BoundaryModeFemOp
{
private:
  MaterialOperator &mat_op;
  FiniteElementSpace &nd_fespace;
  FiniteElementSpace &h1_fespace;

  // L2 curl space for the B-field (Hz = curl_t(Et)/(iωμ), scalar in 2D).
  std::unique_ptr<mfem::L2_FECollection> l2_curl_fec;
  std::unique_ptr<FiniteElementSpace> l2_curl_fespace;

  // Btt = (mu^{-1}) mass matrix on the ND space, for impedance power computation.
  std::unique_ptr<mfem::HypreParMatrix> Btt;

public:
  BoundaryModeFemOp(MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
                    FiniteElementSpace &h1_fespace, Mesh &mesh, int order)
    : mat_op(mat_op), nd_fespace(nd_fespace), h1_fespace(h1_fespace)
  {
    // Create L2 curl space for the scalar B-field (Hz) in 2D.
    l2_curl_fec = std::make_unique<mfem::L2_FECollection>(order - 1, mesh.Dimension(),
                                                          mfem::BasisType::GaussLegendre,
                                                          mfem::FiniteElement::INTEGRAL);
    l2_curl_fespace = std::make_unique<FiniteElementSpace>(mesh, l2_curl_fec.get());
  }

  const MaterialOperator &GetMaterialOp() const { return mat_op; }
  auto &GetNDSpace() { return nd_fespace; }
  const auto &GetNDSpace() const { return nd_fespace; }
  auto &GetH1Space() { return h1_fespace; }
  const auto &GetH1Space() const { return h1_fespace; }

  // L2 curl space for B-field (scalar curl of Et in 2D).
  auto &GetCurlSpace() { return *l2_curl_fespace; }
  const auto &GetCurlSpace() const { return *l2_curl_fespace; }

  const auto &GetMesh() const { return nd_fespace.GetMesh(); }
  MPI_Comm GetComm() const { return nd_fespace.GetComm(); }

  void SetBttMatrix(std::unique_ptr<mfem::HypreParMatrix> Btt_) { Btt = std::move(Btt_); }
  const mfem::HypreParMatrix *GetBtt() const { return Btt.get(); }
};

//
// Driver class for 2D waveguide mode analysis using a linear eigenvalue formulation
// (Eq 1 + Eq 2 with VD substitution). This formulation supports full impedance BC
// handling (both BC-t and BC-n) while maintaining a standard generalized eigenvalue
// problem in kn^2 (no quadratic linearization).
//
class BoundaryModeSolver : public BaseSolver
{
public:
  using BaseSolver::BaseSolver;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BOUNDARY_MODE_SOLVER_HPP
