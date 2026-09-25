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

  // Batched sweeps: solve several excitations at once. Equivalent to calling the single
  // versions in turn, but the environment solves (interface load g_E and interior recovery
  // u_E) run as multi-RHS blocks through one direct factorization.
  std::vector<Vector> SolveExcitations(const std::vector<int> &drive_terminal_indices);
  std::vector<Vector> SolveDirichlets(const std::vector<Vector> &dbc_values);

  // Energy matrix E_ij = u_i^T K u_j (K the pure energy operator) of Dirichlet-lift
  // excitations: lifts[k] is a full parent true-DOF vector prescribing the Dirichlet data
  // (e.g. a terminal unit potential or a magnetostatic flux-loop lift) and ids[k] labels it
  // for model reuse. Computed from region solves against the condensed environment (S_E,
  // the energy operator S^K, and per-lift interface couplings): once the lifts' modes are
  // known (saved with the model and matched by id + an environment-side fingerprint), NO
  // environment solve is needed. Modes are otherwise computed with the environment factor
  // (and appended to the model in an offline run that saves one). Magnetostatic runs
  // without a materialized S^K (a model is not being saved/loaded) use the environment
  // path. If fields is non-null, the full fields of the first n_fields lifts are also
  // recovered (this needs the environment interior).
  mfem::DenseMatrix EnergyMatrix(const std::vector<int> &ids,
                                 const std::vector<Vector> &lifts,
                                 std::vector<Vector> *fields = nullptr, int n_fields = 0);

  // Electrostatic Maxwell capacitance matrix over the given terminals: EnergyMatrix of the
  // terminal unit potentials.
  mfem::DenseMatrix CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                      std::vector<Vector> *fields = nullptr,
                                      int n_fields = 0);

  // Unit Dirichlet lift of a terminal (1 on its true DOFs, 0 elsewhere).
  Vector TerminalLift(int terminal_index) const;

  // Whether the environment interior operator has been factored (diagnostics / tests).
  bool EnvironmentFactored() const;

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
  // Region-condensed Dirichlet-lift solves for a batch of prescribed Dirichlet fields (each
  // a full parent true-DOF vector, nonzero only on the Dirichlet DOFs).
  std::vector<Vector> SolveDirichletBatch(const std::vector<Vector> &dbcs);

  // Hides the heavy internals (LaplaceOperator on each submesh, Substructure,
  // DtNBoundaryOperator, KspSolver) so the header stays decoupled from the plumbing.
  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
