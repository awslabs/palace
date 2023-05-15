// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hypre.hpp"

namespace mfem
{

mfem::HypreParMatrix *
HypreParMatrixFromBlocks(mfem::Array2D<const mfem::HypreParMatrix *> &blocks,
                         mfem::Array2D<double> *coeff)
{
  mfem::Array2D<mfem::HypreParMatrix *> blocks_without_const(blocks.NumRows(),
                                                             blocks.NumCols());
  for (int i = 0; i < blocks.NumRows(); i++)
  {
    for (int j = 0; j < blocks.NumCols(); j++)
    {
      blocks_without_const(i, j) = const_cast<mfem::HypreParMatrix *>(blocks(i, j));
    }
  }
  return HypreParMatrixFromBlocks(blocks_without_const, coeff);
}

}  // namespace mfem
