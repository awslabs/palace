// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "strumpack.hpp"

#if defined(MFEM_USE_STRUMPACK)

namespace palace
{

template <typename StrumpackSolverType>
StrumpackSolverBase<StrumpackSolverType>::StrumpackSolverBase(
    MPI_Comm comm, int sym_fact_type, strumpack::CompressionType comp_type, double lr_tol,
    int butterfly_l, int lossy_prec, int print)
  : StrumpackSolverType(comm)
{
  // Configure the solver.
  this->SetPrintFactorStatistics(print > 1);
  this->SetPrintSolveStatistics(print > 1);
  this->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);  // Always as a preconditioner or
                                                           // direct solver
  this->SetMatching(strumpack::MatchingJob::NONE);
  if (sym_fact_type == 2)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::PARMETIS);
  }
  else if (sym_fact_type == 1)
  {
    this->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
  }
  else
  {
    // Use default
  }
  this->SetReorderingReuse(true);  // Repeated calls use same sparsity pattern

  // Configure compression.
  this->SetCompression(comp_type);
  switch (comp_type)
  {
    case strumpack::CompressionType::LOSSLESS:
    case strumpack::CompressionType::LOSSY:
      if (lossy_prec <= 0)
      {
        this->SetCompression(strumpack::CompressionType::LOSSLESS);
      }
      else
      {
        this->SetCompression(strumpack::CompressionType::LOSSY);
        this->SetCompressionLossyPrecision(lossy_prec);
      }
      break;
    case strumpack::CompressionType::ZFP_BLR_HODLR:
      this->SetCompressionLossyPrecision(lossy_prec);
    case strumpack::CompressionType::HODLR:
    case strumpack::CompressionType::BLR_HODLR:
      this->SetCompressionButterflyLevels(butterfly_l);
    case strumpack::CompressionType::HSS:
    case strumpack::CompressionType::BLR:
      this->SetCompressionRelTol(lr_tol);
      break;
    case strumpack::CompressionType::NONE:
    default:
      break;
  }
}

template <typename StrumpackSolverType>
void StrumpackSolverBase<StrumpackSolverType>::SetOperator(const ParOperator &op)
{
  // Convert the input operator to a distributed STRUMPACK matrix (always assume a symmetric
  // sparsity pattern). Safe to delete the matrix since STRUMPACK copies it on input.
  mfem::STRUMPACKRowLocMatrix A(const_cast<ParOperator *>(&op)->ParallelAssemble(), true);

  // Set up base class.
  StrumpackSolverType::SetOperator(A);
}

template class StrumpackSolverBase<mfem::STRUMPACKSolver>;
template class StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif
