// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP

#include <memory>
#include <vector>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

// Phase 1 substructuring for electrostatics: condense the environment (the rest of the
// domain, selected by domain attributes) to an implicit Dirichlet-to-Neumann boundary
// operator on the shared interface, then solve the region of interest against it. Works in
// parallel (true-DOF interface identification and a distributed implicit DtN).
//
// Scope (this phase): single-excitation electrostatic. The capacitance sweep, reuse-
// optimized (materialized) DtN, and magnetostatics are later phases.
class SubstructuringSolver
{
public:
  SubstructuringSolver(const IoData &iodata,
                       const std::vector<std::unique_ptr<Mesh>> &mesh);
  ~SubstructuringSolver();

  // Assemble and condense the environment to S_E, g_E (offline). Called once; the result is
  // reused by every SolveRegion call.
  void CondenseEnvironment();

  // Solve the region-condensed electrostatic problem: the region's own Dirichlet terminals
  // plus the environment DtN term on the interface. Returns the potential on the parent H1
  // true DOFs (the region-touched entries are the region solution; environment-interior
  // entries are left zero). Requires CondenseEnvironment to have been called.
  Vector SolveRegion();

  // Global parent H1 true-DOF size, for reporting.
  long long int RegionGlobalTrueVSize() const;

private:
  // Hides the heavy internals (LaplaceOperator on each submesh, Substructure,
  // DtNBoundaryOperator, KspSolver) so the header stays decoupled from the plumbing.
  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
