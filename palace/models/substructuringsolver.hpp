// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP

#include <memory>
#include <string>
#include <vector>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

// Phase 1 substructuring for electrostatics: condense the environment (the rest of the
// domain, selected by domain attributes) to an implicit Dirichlet-to-Neumann boundary
// operator on the shared interface, then solve the region of interest against it. The
// environment DtN is materialized once and reused across terminal excitations, enabling a
// cheap capacitance sweep. Works in parallel (true-DOF interface identification and a
// distributed implicit DtN).
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
  // plus the environment DtN term on the interface, then recover the environment interior.
  // Returns the full parent H1 true-DOF potential. Requires CondenseEnvironment first.
  Vector SolveRegion();

  // Solve for a single excitation: the given terminal driven to 1 V, all other terminals
  // grounded. Reuses the materialized environment DtN. Returns the full parent field.
  Vector SolveExcitation(int drive_terminal_index);

  // Solve a Dirichlet-lift excitation with an arbitrary prescribed boundary field on the
  // Dirichlet DOF set (e.g. a magnetostatic flux-loop lift). Returns the full parent field.
  Vector SolveDirichlet(const Vector &dbc_values);

  // Solve K u = f for a full parent-space source f (magnetostatic current excitation),
  // reusing the condensed environment. Returns the full parent field.
  Vector SolveSource(const Vector &f);

  // Terminal indices (sorted), for a capacitance sweep.
  std::vector<int> TerminalIndices() const;

  // Mutual energy ui^T K uj over the full (region + environment) stiffness; the (i, j)
  // Maxwell capacitance entry for unit terminal excitations.
  double MutualEnergy(const Vector &ui, const Vector &uj) const;

  // Total electrostatic energy 1/2 phi^T K phi of a full parent-space potential, using the
  // full (region + environment) stiffness. For a single-terminal 1 V excitation this is
  // half the driven terminal's self-capacitance.
  double ElectrostaticEnergy(const Vector &u) const;

  // Global parent H1 true-DOF size, for reporting.
  long long int RegionGlobalTrueVSize() const;

  // Write the recovered full parent-space potentials to a ParaView collection under dir,
  // one time step per excitation (time = terminal index).
  void WriteParaView(const std::string &dir, const std::vector<int> &terminals,
                     const std::vector<Vector> &fields) const;

private:
  // Region-condensed Dirichlet-lift solve using the currently set impl->dbc_values.
  Vector SolveWithCurrentDbc();

  // Hides the heavy internals (LaplaceOperator on each submesh, Substructure,
  // DtNBoundaryOperator, KspSolver) so the header stays decoupled from the plumbing.
  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
