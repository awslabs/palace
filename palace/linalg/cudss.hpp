// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_CUDSS_HPP
#define PALACE_LINALG_CUDSS_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_CUDSS)

#include "utils/labels.hpp"

namespace palace
{

//
// A wrapper for cuDSS.
//
class CuDSSSolver : public mfem::CuDSSSolver
{
public:
  CuDSSSolver(MPI_Comm comm, MatrixSymmetry sym, SymbolicFactorization reorder,
              bool reorder_reuse, int print);

  void SetReorderReuse(bool reorder_reuse);
};

}  // namespace palace

#endif

#endif  // PALACE_LINALG_CUDSS_HPP
