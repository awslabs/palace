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
// domain, selected by domain attributes) to a Dirichlet-to-Neumann boundary operator
// S_E, g_E on the shared interface, then solve the region of interest against it with
// Palace's linear solver. The environment condensation is performed once and reused across
// region solves, enabling region redesign without re-solving the environment.
//
// Scope (this phase): single-excitation electrostatic, offline mode. The capacitance sweep
// (multiple terminals, environment-terminal excitation), online/serialized reuse, and
// magnetostatics are later phases.
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
  // plus the environment DtN term on the interface. Returns the region potential on the
  // region H1 (true) DOFs. Requires CondenseEnvironment to have been called.
  Vector SolveRegion();

  // Global (region) H1 true-DOF size, for reporting.
  long long int RegionGlobalTrueVSize() const;

private:
  // Hides the heavy internals (LaplaceOperator on each submesh, Substructure,
  // DtNBoundaryOperator, KspSolver) so the header stays decoupled from the plumbing.
  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
