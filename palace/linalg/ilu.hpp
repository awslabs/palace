// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_ILU_HPP
#define PALACE_LINALG_ILU_HPP

#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

//
// A wrapper for Hypre's ILU solver.
//
class ILUSolver : public mfem::HypreILU
{
public:
  ILUSolver(int type, int fill_level = 1, int print = 0);
  ILUSolver(const IoData &iodata, int print)
    : ILUSolver(iodata.solver.linear.ilu_type, iodata.solver.linear.ilu_fill_level, print)
  {
  }
};

}  // namespace palace

#endif  // PALACE_LINALG_ILU_HPP
