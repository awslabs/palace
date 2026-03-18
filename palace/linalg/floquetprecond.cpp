// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "floquetprecond.hpp"

#include "utils/communication.hpp"

namespace palace
{

FloquetBoundaryPreconditioner::FloquetBoundaryPreconditioner(
    MPI_Comm port_comm, std::unique_ptr<Solver<ComplexOperator>> &&pc_inner,
    const mfem::Array<int> &bdr_tdof_list_in, int n_bdr_global,
    const std::vector<int> &recv_counts, const std::vector<int> &recv_displs)
  : Solver<ComplexOperator>(), pc_inner(std::move(pc_inner)), port_comm(port_comm),
    n_bdr_global(n_bdr_global), recv_counts(recv_counts), recv_displs(recv_displs),
    x_bdr_global(n_bdr_global), y_bdr_global(n_bdr_global)
{
  bdr_tdof_list.MakeRef(bdr_tdof_list_in);
  n_bdr_local = bdr_tdof_list.Size();
  x_bdr_local.resize(n_bdr_local);
}

FloquetBoundaryPreconditioner::~FloquetBoundaryPreconditioner() = default;

void FloquetBoundaryPreconditioner::SetOperator(const ComplexOperator &op)
{
  // Delegate to inner preconditioner. The boundary Schur complement is updated
  // separately via UpdateBoundarySchur().
  pc_inner->SetOperator(op);
  this->height = op.Height();
  this->width = op.Width();
}

void FloquetBoundaryPreconditioner::UpdateBoundarySchur(const Eigen::MatrixXcd &S)
{
  S_lu.compute(S);
  factored = true;
}

void FloquetBoundaryPreconditioner::Mult(const ComplexVector &x, ComplexVector &y) const
{
  // Step 1: Apply inner preconditioner (handles interior DOFs).
  pc_inner->Mult(x, y);

  if (!factored || port_comm == MPI_COMM_NULL)
  {
    return;  // No boundary correction (non-port rank or not yet factored).
  }

  // Step 2: Extract x at local boundary DOFs.
  x.Real().HostRead();
  x.Imag().HostRead();
  for (int i = 0; i < n_bdr_local; i++)
  {
    int dof = bdr_tdof_list[i];
    x_bdr_local(i) = std::complex<double>(x.Real()[dof], x.Imag()[dof]);
  }

  // Step 3: Gather boundary DOFs across port_comm.
  MPI_Allgatherv(x_bdr_local.data(), n_bdr_local, MPI_DOUBLE_COMPLEX, x_bdr_global.data(),
                 recv_counts.data(), recv_displs.data(), MPI_DOUBLE_COMPLEX, port_comm);

  // Step 4: Dense boundary solve.
  y_bdr_global.noalias() = S_lu.solve(x_bdr_global);

  // Step 5: Scatter back — overwrite y at local boundary DOFs.
  y.Real().HostReadWrite();
  y.Imag().HostReadWrite();
  int offset = recv_displs[Mpi::Rank(port_comm)];
  for (int i = 0; i < n_bdr_local; i++)
  {
    int dof = bdr_tdof_list[i];
    y.Real()[dof] = y_bdr_global(offset + i).real();
    y.Imag()[dof] = y_bdr_global(offset + i).imag();
  }
}

}  // namespace palace
