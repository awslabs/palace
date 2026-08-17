// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_MMS_SOLVERS_HPP
#define PALACE_TEST_MMS_SOLVERS_HPP

// Test-only solver subclasses for manufactured-solution (MMS) verification. Each re-exposes
// a driver's protected Solve overloads as public so unit tests can drive the solve directly
// and inspect the solution vector, without enlarging the shipping API. Add subclasses for
// further MMS-verified solvers (magnetostatic, etc.) here as they are covered.

#include "drivers/electrostaticsolver.hpp"

namespace palace
{

// Re-exposes ElectrostaticSolver's protected Solve overloads as public for unit tests.
class ExposedElectrostaticSolver : public ElectrostaticSolver
{
public:
  using ElectrostaticSolver::ElectrostaticSolver;
  using ElectrostaticSolver::Solve;
};

}  // namespace palace

#endif  // PALACE_TEST_MMS_SOLVERS_HPP
