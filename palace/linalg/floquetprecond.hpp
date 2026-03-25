// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_FLOQUET_PRECOND_HPP
#define PALACE_LINALG_FLOQUET_PRECOND_HPP

#include <complex>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Block preconditioner for Floquet port DtN boundary conditions.
//
// Wraps an existing volume preconditioner P_inner and adds a dense boundary solve:
//   1. y = P_inner^{-1}(x)         — handles interior DOFs via AMS/direct/multigrid
//   2. y_B = S^{-1}(x_B)           — overwrites boundary DOFs with exact dense solve
//
// S is the boundary Schur complement: S = A_BB + F_BB where A_BB is the boundary-boundary
// block of the assembled system and F_BB is the dense DtN matrix. S is pre-factored via
// Eigen LU (N_bdr × N_bdr complex). The dense solve uses port_comm (sub-communicator of
// ranks owning port boundary elements) to avoid redundant work on non-port ranks.
//
class FloquetBoundaryPreconditioner : public Solver<ComplexOperator>
{
private:
  // Inner (volume) preconditioner — takes ownership.
  std::unique_ptr<Solver<ComplexOperator>> pc_inner;

  // Port sub-communicator: only ranks owning Floquet boundary DOFs participate.
  // MPI_COMM_NULL on non-port ranks.
  MPI_Comm port_comm;

  // Local boundary true DOF indices into the full vector.
  mfem::Array<int> bdr_tdof_list;
  int n_bdr_local;
  int n_bdr_global;

  // Allgatherv infrastructure on port_comm.
  std::vector<int> recv_counts, recv_displs;

  // Pre-factored boundary Schur complement.
  Eigen::PartialPivLU<Eigen::MatrixXcd> S_lu;
  bool factored = false;

  // Workspace vectors.
  mutable Eigen::VectorXcd x_bdr_local, x_bdr_global, y_bdr_global;

public:
  // Construct with the inner preconditioner (takes ownership) and boundary DOF info.
  // port_comm, bdr_tdof_list, n_bdr_global, recv_counts, recv_displs come from the
  // FloquetPortData (or are aggregated across ports by the caller).
  FloquetBoundaryPreconditioner(MPI_Comm port_comm,
                                std::unique_ptr<Solver<ComplexOperator>> &&pc_inner,
                                const mfem::Array<int> &bdr_tdof_list, int n_bdr_global,
                                const std::vector<int> &recv_counts,
                                const std::vector<int> &recv_displs);

  ~FloquetBoundaryPreconditioner() override;

  // Delegates to inner preconditioner.
  void SetOperator(const ComplexOperator &op) override;

  // Update the boundary Schur complement and factor it.
  // S should be the N_bdr_global × N_bdr_global complex boundary system matrix.
  void UpdateBoundarySchur(const Eigen::MatrixXcd &S);

  // Apply: y = P_inner^{-1}(x), then overwrite y_B = S^{-1}(x_B).
  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

}  // namespace palace

#endif  // PALACE_LINALG_FLOQUET_PRECOND_HPP
