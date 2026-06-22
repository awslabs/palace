// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#include "cudss.hpp"

#if defined(MFEM_USE_CUDSS)

#include "utils/iodata.hpp"

namespace palace
{

namespace
{

mfem::CuDSSSolver::MatType GetCuDSSMatType(MatrixSymmetry sym)
{
  switch (sym)
  {
    case MatrixSymmetry::SPD:
      return mfem::CuDSSSolver::SYMMETRIC_POSITIVE_DEFINITE;
    case MatrixSymmetry::SYMMETRIC:
      return mfem::CuDSSSolver::SYMMETRIC_INDEFINITE;
    case MatrixSymmetry::UNSYMMETRIC:
      return mfem::CuDSSSolver::NONSYMMETRIC;
  }
  return mfem::CuDSSSolver::NONSYMMETRIC;
}

}  // namespace

CuDSSSolver::CuDSSSolver(MPI_Comm comm, MatrixSymmetry sym, SymbolicFactorization reorder,
                         bool reorder_reuse, int print)
  : mfem::CuDSSSolver(comm)
{
  SetMatrixSymType(GetCuDSSMatType(sym));
  SetReorderingReuse(reorder_reuse);
}

void CuDSSSolver::SetReorderReuse(bool reorder_reuse)
{
  SetReorderingReuse(reorder_reuse);
}

}  // namespace palace

#endif
