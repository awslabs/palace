// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MUMPS_HPP
#define PALACE_LINALG_MUMPS_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_MUMPS)

#include "utils/labels.hpp"

namespace palace
{

class IoData;

//
// A wrapper for the MUMPS direct solver package.
//
class MumpsSolver : public mfem::MUMPSSolver
{
public:
  MumpsSolver(MPI_Comm comm, MatrixSymmetry sym, SymbolicFactorization reorder,
              double blr_tol, bool reorder_reuse, int print);
  MumpsSolver(const IoData &iodata, MPI_Comm comm, int print);

  void SetReorderReuse(bool reorder_reuse);
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_MUMPS_HPP
