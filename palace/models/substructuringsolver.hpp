// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP

#include <memory>
#include <string>
#include <vector>
#include "fem/errorindicator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;

// Substructuring for electrostatics and magnetostatics: condense the environment (the rest
// of the domain, selected by domain attributes) to a Dirichlet-to-Neumann operator S_E on
// the interface with the region of interest, and solve the region against it. S_E is
// materialized once (and can be saved for online runs on a redesigned region), so the
// capacitance or inductance sweep needs region solves only.
class SubstructuringSolver
{
public:
  SubstructuringSolver(const IoData &iodata,
                       const std::vector<std::unique_ptr<Mesh>> &mesh);
  ~SubstructuringSolver();

  // Condense the environment to S_E (or load it from a saved model). Called once, before
  // any solve.
  void CondenseEnvironment();

  // Solve for the given terminal driven to 1 V, all other terminals grounded. Returns the
  // full parent field.
  Vector SolveExcitation(int drive_terminal_index);

  // Solve for an arbitrary prescribed field on the Dirichlet DOFs. Returns the full parent
  // field.
  Vector SolveDirichlet(const Vector &dbc_values);

  // Batched versions of the above (the environment solves run as multi-RHS blocks).
  std::vector<Vector> SolveExcitations(const std::vector<int> &drive_terminal_indices);
  std::vector<Vector> SolveDirichlets(const std::vector<Vector> &dbc_values);

  // Energy matrix E_ij = u_i^T K u_j (K the energy operator) of Dirichlet-lift excitations
  // (lifts[k] a full parent true-DOF vector with the Dirichlet data, labeled ids[k]). Once
  // the lift modes are known (from a saved model, matched by id and an environment-side
  // fingerprint), no environment solve is needed; otherwise they are computed and appended
  // to the model being saved. If fields is non-null, it receives the full fields of the
  // first n_fields lifts; if region_fields is non-null, the solutions of all lifts on the
  // region and interface (zero on the environment interior).
  mfem::DenseMatrix EnergyMatrix(const std::vector<int> &ids,
                                 const std::vector<Vector> &lifts,
                                 std::vector<Vector> *fields = nullptr, int n_fields = 0,
                                 std::vector<Vector> *region_fields = nullptr);

  // Magnetostatic energy matrix of London flux states: flux loop k with the fluxoid
  // generator a_k drives the source M_sheet a_k with the superconducting films as free
  // London sheets, and E_ij = u_i^T K u_j + (u_i - a_i)^T M_sheet (u_j - a_j) as in the
  // native solver. If fields is non-null, it receives the full fields of the first n_fields
  // states.
  mfem::DenseMatrix SheetEnergyMatrix(const std::vector<int> &ids,
                                      const std::vector<Vector> &a,
                                      std::vector<Vector> *fields = nullptr,
                                      int n_fields = 0);

  // Magnetostatic energy matrix of surface-current excitations: port k drives the assembled
  // excitation J_k (zero on the Dirichlet DOFs), and E_ij = u_i^T K u_j (+ the sheet
  // kinetic energy u_i^T M_sheet u_j), in a form stationary in the regularization of the
  // solve. A current that does not close is rejected.
  mfem::DenseMatrix CurrentEnergyMatrix(const std::vector<int> &ids,
                                        const std::vector<Vector> &J,
                                        std::vector<Vector> *fields = nullptr,
                                        int n_fields = 0);

  // Magnetostatic energy matrix of surface currents and London flux states (currents
  // first), with the blocks of CurrentEnergyMatrix and SheetEnergyMatrix and a zero cross
  // block (the two kinds of states are energy-orthogonal). With the aperture flux
  // functionals l_c of the ports, linked_flux(c, f) = l_c^T u_f for the flux states f.
  mfem::DenseMatrix MagnetostaticEnergyMatrix(
      const std::vector<int> &current_ids, const std::vector<Vector> &J,
      const std::vector<int> &flux_ids, const std::vector<Vector> &a,
      const std::vector<Vector> &apertures, mfem::DenseMatrix *linked_flux,
      std::vector<Vector> *fields = nullptr, int n_fields = 0);

  // Whether the model has London superconductor sheets (magnetostatics).
  bool HasSheets() const;

  // Maxwell capacitance matrix over the given terminals (EnergyMatrix of the terminal unit
  // potentials).
  mfem::DenseMatrix CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                      std::vector<Vector> *fields = nullptr,
                                      int n_fields = 0,
                                      std::vector<Vector> *region_fields = nullptr);

  // Error indicators for adaptive refinement of the region only (electrostatics): Palace's
  // gradient-flux estimator on the region solutions of CapacitanceMatrix, normalized by the
  // energies E_kk / 2, and zero on the environment.
  ErrorIndicator RegionErrorIndicator(const std::vector<Vector> &region_fields,
                                      const mfem::DenseMatrix &E) const;

  // Unit Dirichlet lift of a terminal (1 on its true DOFs, 0 elsewhere).
  Vector TerminalLift(int terminal_index) const;

  // Whether the environment interior solver has been set up.
  bool EnvironmentFactored() const;

  // Solve A u = f for a full parent-space source f. Returns the full parent field.
  Vector SolveSource(const Vector &f);

  // Terminal indices (sorted).
  std::vector<int> TerminalIndices() const;

  // Mutual energy u_i^T K u_j with the full (region + environment) energy operator.
  double MutualEnergy(const Vector &ui, const Vector &uj) const;

  // Energy 1/2 u^T K u with the full energy operator.
  double ElectrostaticEnergy(const Vector &u) const;

  // Global parent true-DOF size.
  long long int GlobalTrueVSize() const;

  // Write full parent-space fields (V or A) to a ParaView collection under dir, one time
  // step per excitation (time = excitation index).
  void WriteParaView(const std::string &dir, const std::vector<int> &ids,
                     const std::vector<Vector> &fields) const;

private:
  // Region-condensed solves for a batch of prescribed Dirichlet fields.
  std::vector<Vector> SolveDirichletBatch(const std::vector<Vector> &dbcs);

  // EnergyMatrix for Dirichlet data given on the Dirichlet true DOFs only.
  mfem::DenseMatrix EnergyMatrixDbc(const std::vector<int> &ids,
                                    const std::vector<Vector> &xd,
                                    std::vector<Vector> *fields, int n_fields,
                                    std::vector<Vector> *region_fields);

  // A magnetostatic source: a London flux state drives M_sheet a (a the fluxoid generator),
  // a surface current the assembled excitation J. Either may be null.
  struct Source
  {
    const Vector *a, *J;
    bool HasA() const { return a && a->Size() > 0; }
    bool HasJ() const { return J && J->Size() > 0; }
  };

  // Linear functionals l_k (with mode ids) of the solutions: values(k, j) = l_k^T u_j.
  struct Functionals
  {
    std::vector<int> ids;
    std::vector<const Vector *> l;
    mfem::DenseMatrix values;
  };

  // Energy G of magnetostatic sources (see MagnetostaticEnergyMatrix). If work is non-null,
  // it receives the source works b_k^T u_k; fun receives its values.
  mfem::DenseMatrix SourceEnergyMatrix(const std::vector<int> &ids,
                                       const std::vector<Source> &src,
                                       std::vector<Vector> *fields, int n_fields,
                                       std::vector<double> *work = nullptr,
                                       Functionals *fun = nullptr);

  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
