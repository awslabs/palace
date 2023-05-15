// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_HYPRE_HPP
#define PALACE_LINALG_HYPRE_HPP

#include <mfem.hpp>

// XX TODO REVISIT AFTER WAVE PORT OPERATOR IF NEEDED....

namespace mfem
{

// Convenience wrapper for casting away the const on the pointers and dispatching onto the
// original function that has the argument type: mfem::Array2D<mfem::HypreParMatrix*> &.
mfem::HypreParMatrix *
HypreParMatrixFromBlocks(mfem::Array2D<const mfem::HypreParMatrix *> &blocks,
                         mfem::Array2D<double> *coeff = nullptr);

}  // namespace mfem

#endif  // PALACE_LINALG_HYPRE_HPP
