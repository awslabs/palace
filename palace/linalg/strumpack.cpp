// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#if defined(MFEM_USE_STRUMPACK)

namespace palace
{

namespace
{

strumpack::CompressionType
GetCompressionType(config::LinearSolverData::CompressionType type)
{
  switch (type)
  {
    case config::LinearSolverData::CompressionType::HSS:
      return strumpack::CompressionType::HSS;
    case config::LinearSolverData::CompressionType::BLR:
      return strumpack::CompressionType::BLR;
    case config::LinearSolverData::CompressionType::HODLR:
      return strumpack::CompressionType::HODLR;
    case config::LinearSolverData::CompressionType::ZFP:
      return strumpack::CompressionType::LOSSY;
    case config::LinearSolverData::CompressionType::BLR_HODLR:
      return strumpack::CompressionType::BLR_HODLR;
      break;
    case config::LinearSolverData::CompressionType::ZFP_BLR_HODLR:
      return strumpack::CompressionType::ZFP_BLR_HODLR;
      break;
    case config::LinearSolverData::CompressionType::NONE:
      return strumpack::CompressionType::NONE;
  }
  return strumpack::CompressionType::NONE;  // For compiler warning
}

}  // namespace

template <typename StrumpackSolverType>
StrumpackSolverBase<StrumpackSolverType>::StrumpackSolverBase(
    MPI_Comm comm, config::LinearSolverData::SymFactType reorder,
    config::LinearSolverData::CompressionType compression, double lr_tol, int butterfly_l,
    int lossy_prec, int print)
  : StrumpackSolverType(comm), comm(comm)
{
  // Configure the solver.
  this->SetPrintFactorStatistics(print > 1);
  this->SetPrintSolveStatistics(print > 1);
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  switch (reorder)
  {
    case config::LinearSolverData::SymFactType::METIS:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
      // this->SetReorderingStrategy(strumpack::ReorderingStrategy::AND);
      break;
    case config::LinearSolverData::SymFactType::PARMETIS:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::PARMETIS);
      break;
    case config::LinearSolverData::SymFactType::SCOTCH:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::SCOTCH);
      break;
    case config::LinearSolverData::SymFactType::PTSCOTCH:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::PTSCOTCH);
      break;
    case config::LinearSolverData::SymFactType::AMD:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::AMD);
      // this->SetReorderingStrategy(strumpack::ReorderingStrategy::MMD);
      break;
    case config::LinearSolverData::SymFactType::RCM:
      this->SetReorderingStrategy(strumpack::ReorderingStrategy::RCM);
    case config::LinearSolverData::SymFactType::PORD:
    case config::LinearSolverData::SymFactType::DEFAULT:
      // Should have good default
      break;
  }
  this->SetReorderingReuse(true);  // Repeated calls use same sparsity pattern

  // Configure compression.
  this->SetCompression(GetCompressionType(compression));
  switch (compression)
  {
    case config::LinearSolverData::CompressionType::ZFP:
      if (lossy_prec <= 0)
      {
        this->SetCompression(strumpack::CompressionType::LOSSLESS);
      }
      else
      {
        this->SetCompressionLossyPrecision(lossy_prec);
      }
      break;
    case config::LinearSolverData::CompressionType::ZFP_BLR_HODLR:
      this->SetCompressionLossyPrecision(lossy_prec);
    case config::LinearSolverData::CompressionType::HODLR:
    case config::LinearSolverData::CompressionType::BLR_HODLR:
      this->SetCompressionButterflyLevels(butterfly_l);
    case config::LinearSolverData::CompressionType::HSS:
    case config::LinearSolverData::CompressionType::BLR:
      this->SetCompressionRelTol(lr_tol);
      break;
    case config::LinearSolverData::CompressionType::NONE:
      break;
  }
  // if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  // {
  //   this->EnableGPU();  // XX TODO: GPU support disabled for now
  // }
  // else
  {
    this->DisableGPU();
  }
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetOperator(const Operator &op)
{
  // Convert the input operator to a distributed STRUMPACK matrix (always assume a symmetric
  // sparsity pattern). This is very similar to the MFEM's STRUMPACKRowLocMatrix from a
  // HypreParMatrix but avoids using the communicator from the Hypre matrix in the case that
  // the solver is constructed on a different communicator.
  const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(hA && hA->GetGlobalNumRows() == hA->GetGlobalNumCols(),
              "StrumpackSolver requires a square HypreParMatrix operator!");
  auto *parcsr = (hypre_ParCSRMatrix *)const_cast<mfem::HypreParMatrix &>(*hA);
  hypre_CSRMatrix *csr = hypre_MergeDiagAndOffd(parcsr);
  hypre_CSRMatrixMigrate(csr, HYPRE_MEMORY_HOST);

  // Create the STRUMPACKRowLocMatrix by taking the internal data from a hypre_CSRMatrix.
  HYPRE_BigInt glob_n = hypre_ParCSRMatrixGlobalNumRows(parcsr);
  HYPRE_BigInt first_row = hypre_ParCSRMatrixFirstRowIndex(parcsr);
  HYPRE_Int n_loc = hypre_CSRMatrixNumRows(csr);
  HYPRE_Int *I = hypre_CSRMatrixI(csr);
  HYPRE_BigInt *J = hypre_CSRMatrixBigJ(csr);
  double *data = hypre_CSRMatrixData(csr);

  // Safe to delete the matrix since STRUMPACK copies it on input. Also clean up the Hypre
  // data structure once we are done with it.
#if !defined(HYPRE_BIGINT)
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc, first_row, glob_n, glob_n, I, J, data, true);
#else
  int n_loc_int = static_cast<int>(n_loc);
  MFEM_ASSERT(n_loc == (HYPRE_Int)n_loc_int,
              "Overflow error for local sparse matrix size!");
  mfem::Array<int> II(n_loc_int + 1);
  for (int i = 0; i <= n_loc_int; i++)
  {
    II[i] = static_cast<int>(I[i]);
    MFEM_ASSERT(I[i] == (HYPRE_Int)II[i], "Overflow error for local sparse matrix index!");
  }
  mfem::STRUMPACKRowLocMatrix A(comm, n_loc_int, first_row, glob_n, glob_n, II.HostRead(),
                                J, data, true);
#endif
  StrumpackSolverType::SetOperator(A);
  hypre_CSRMatrixDestroy(csr);
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif
