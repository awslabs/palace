// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_HPP
#define PALACE_LINALG_STRUMPACK_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include "linalg/operator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for the STRUMPACK direct solver package.
//
template <typename StrumpackSolverType>
class StrumpackSolverBase : public StrumpackSolverType
{
private:
  strumpack::CompressionType CompressionType(config::LinearSolverData::CompressionType type)
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
      default:
        return strumpack::CompressionType::NONE;
    }
  }

public:
  StrumpackSolverBase(MPI_Comm comm, int sym_fact_type,
                      strumpack::CompressionType comp_type, double lr_tol, int butterfly_l,
                      int lossy_prec, int print_lvl);

  StrumpackSolverBase(MPI_Comm comm, const IoData &iodata, int print_lvl)
    : StrumpackSolverBase(comm,
                          (iodata.solver.linear.sym_fact_type ==
                           config::LinearSolverData::SymFactType::PARMETIS)
                              ? 2
                              : ((iodata.solver.linear.sym_fact_type ==
                                  config::LinearSolverData::SymFactType::METIS)
                                     ? 1
                                     : 0),
                          CompressionType(iodata.solver.linear.strumpack_compression_type),
                          iodata.solver.linear.strumpack_lr_tol,
                          iodata.solver.linear.strumpack_butterfly_l,
                          iodata.solver.linear.strumpack_lossy_precision, print_lvl)
  {
  }

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("StrumpackSolver requires a ParOperator operator!");
  }
  void SetOperator(const ParOperator &op);
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

#if STRUMPACK_VERSION_MAJOR >= 6 && STRUMPACK_VERSION_MINOR >= 3 && \
    STRUMPACK_VERSION_PATCH > 1
using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;
#endif

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
