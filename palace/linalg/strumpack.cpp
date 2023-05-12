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
    case config::LinearSolverData::CompressionType::INVALID:
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
  this->SetPrintFactorStatistics((print > 1));
  this->SetPrintSolveStatistics((print > 1));
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  if (reorder == config::LinearSolverData::SymFactType::METIS)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PARMETIS)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::PARMETIS);
  }
  else if (reorder == config::LinearSolverData::SymFactType::SCOTCH)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::SCOTCH);
  }
  else if (reorder == config::LinearSolverData::SymFactType::PTSCOTCH)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::PTSCOTCH);
  }
  else
  {
    // Use default
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
    default:
      break;
  }
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetOperator(const Operator &op)
{
  // Convert the input operator to a distributed STRUMPACK matrix (always assume a symmetric
  // sparsity pattern). Safe to delete the matrix since STRUMPACK copies it on input.
  auto *PtAP = const_cast<ParOperator *>(dynamic_cast<const ParOperator *>(&op));
  MFEM_VERIFY(PtAP, "StrumpackSolver requires a ParOperator operator!");
  mfem::STRUMPACKRowLocMatrix A(PtAP->ParallelAssemble(), true);

  // Set up base class.
  StrumpackSolverType::SetOperator(A);
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
#if STRUMPACK_VERSION_MAJOR >= 6 && STRUMPACK_VERSION_MINOR >= 3 && \
    STRUMPACK_VERSION_PATCH > 1
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;
#endif

}  // namespace palace

#endif
