// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_COMMON_HPP
#define PALACE_TEST_MMS_COMMON_HPP

// Helpers shared by the Method of Manufactured Solutions (MMS) tests for the different
// operators (electrostatic Laplace, magnetostatic curl-curl, ...).

#include <cmath>
#include <cstdlib>
#include <string>

#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::testing::mms
{

// Observed order of convergence between two refinement levels, assuming e ~ C h^rate.
inline double ConvergenceRate(double h0, double e0, double h1, double e1)
{
  return std::log(e0 / e1) / std::log(h0 / h1);
}

// Print one CSV row of convergence data when PALACE_MMS_REPORT is set, for plotting. Each
// test reports its two error measures (for example potential and field for electrostatics,
// or field and curl for curl-curl) in the last two columns.
inline void ReportMmsData(const std::string &dataset, const std::string &solution,
                          const std::string &element, int order,
                          const std::string &resolution_kind, int resolution, double h,
                          HYPRE_BigInt primary_dofs, double error_0, double error_1)
{
  if (std::getenv("PALACE_MMS_REPORT"))
  {
    Mpi::Print("MMS_DATA,{},{},{},{},{},{},{:.17g},{},{:.17g},{:.17g}\n", dataset, solution,
               element, order, resolution_kind, resolution, h, primary_dofs, error_0,
               error_1);
  }
}

}  // namespace palace::testing::mms

#endif  // PALACE_TEST_MMS_COMMON_HPP
