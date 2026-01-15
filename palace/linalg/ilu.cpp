// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ilu.hpp"

namespace palace
{

ILUSolver::ILUSolver(int type, int fill_level, int print)
  : mfem::HypreILU()
{
  // Hypre ILU factorization types:
  // 0: Block-Jacobi ILUK (BJ-ILUK).
  // 1: Block-Jacobi ILUT (BJ-ILUT).
  // 10: GMRES with ILUK (GMRES-ILUK).
  // 11: GMRES with ILUT (GMRES-ILUT).
  // 20: NSH with ILUK (NSH-ILUK).
  // 21: NSH with ILUT (NSH-ILUT).
  // 30: RAS with ILUK (RAS-ILUK).
  // 31: RAS with ILUT (RAS-ILUT).
  // 40: ddPQ-GMRES with ILUK (ddPQ-GMRES-ILUK).
  // 41: ddPQ-GMRES with ILUT (ddPQ-GMRES-ILUT).
  // 50: GMRES with RAP-ILU0 with modified ILU0 (GMRES-RAP-ILU0).
  Mpi::Print("Setting ILU type: {} and fill level: {}\n", type, fill_level);
  HYPRE_ILUSetType(*this, type);
  HYPRE_ILUSetLevelOfFill(*this, fill_level); // for ILU(k)
  HYPRE_ILUSetPrintLevel(*this, (print > 1) ? print - 1 : 0);
}

}  // namespace palace
