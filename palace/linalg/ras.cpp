// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ras.hpp"

namespace palace
{

RasSolver::RasSolver(int fill_level, int print) : mfem::HypreILU()
{
  MFEM_VERIFY(fill_level >= 0, "RAS fill level must be nonnegative!");
  MFEM_VERIFY(!mfem::Device::Allows(mfem::Backend::DEVICE_MASK),
              "The hypre RAS-ILU preconditioner does not support GPU execution in "
              "Palace because hypre requires unified memory for this solver!");

  // hypre ILU type 30 is restricted additive Schwarz with ILU(k) subdomain solves. The
  // type is selected in SetOperator because hypre RAS requires an off-process coupling.
  SetLevelOfFill(fill_level);
  SetMaxIter(1);
  SetTol(0.0);
  SetPrintLevel((print > 1) ? print - 1 : 0);
}

void RasSolver::SetOperator(const mfem::Operator &op)
{
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(hA, "RAS requires an assembled HypreParMatrix operator!");

  // RAS setup assumes that the matrix has an off-process communication package. For a
  // serial or globally block-diagonal matrix there is no interface, and one block ILU
  // solve is exactly the corresponding degenerate domain decomposition method.
  mfem::SparseMatrix offd;
  HYPRE_BigInt *cmap = nullptr;
  hA->GetOffd(offd, cmap);
  int has_offd = (offd.Width() > 0);
  MPI_Allreduce(MPI_IN_PLACE, &has_offd, 1, MPI_INT, MPI_LOR, hA->GetComm());
  SetType(has_offd ? 30 : 0);

  mfem::HypreILU::SetOperator(op);
}

}  // namespace palace
