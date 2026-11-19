// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
#define PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "fem/errorindicator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class CurlCurlOperator;
class FiniteElementSpaceHierarchy;
class IoData;
class LaplaceOperator;
class Mesh;

namespace config
{
struct BoundaryData;
}  // namespace config

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

  // Online, with a saved model: Model.Lc from the model (that of its offline run), which a
  // configured one must equal. Collective.
  static void UseModelLength(IoData &iodata, MPI_Comm comm);

  // The boundary conditions of the configuration for the operators on this mesh: without
  // the terminal attributes it does not hold and, on the mesh of the region, without the
  // flux loops it leaves out and the superconductor attributes it does not hold (the model
  // holds their parts). Collective.
  static config::BoundaryData MeshBoundaries(const IoData &iodata,
                                             const mfem::ParMesh &mesh);

  // Condense the environment to S_E (or load it from a saved model). Called once, before
  // any solve.
  void CondenseEnvironment();

  // Maxwell capacitance matrix over the given terminals (EnergyMatrixDbc of the terminal
  // unit potentials). If fields is non-null, it receives the full fields of the first
  // n_fields terminals; if region_fields is non-null, it receives the solutions of all
  // terminals on the region and interface (zero on the environment interior).
  mfem::DenseMatrix CapacitanceMatrix(const std::vector<int> &terminal_indices,
                                      std::vector<Vector> *fields = nullptr,
                                      int n_fields = 0,
                                      std::vector<Vector> *region_fields = nullptr);

  // Magnetostatic energy matrix of surface currents and London flux states (currents
  // first). Surface current k drives the assembled excitation J_k (zero on the Dirichlet
  // DOFs), with E_ij = u_i^T K u_j (+ the sheet kinetic energy u_i^T M_sheet u_j) in a form
  // stationary in the regularization; a current that does not close is rejected. Flux loop
  // k with the fluxoid generator a_k drives M_sheet a_k with the films as free London
  // sheets, with E_ij = u_i^T K u_j + (u_i - a_i)^T M_sheet (u_j - a_j) as in the native
  // solver. The cross block is zero (the two kinds of states are energy-orthogonal). With
  // the aperture flux functionals l_c of the ports, linked_flux(c, f) = l_c^T u_f. If
  // non-null, fields receives the full fields of the first n_fields states, region_fields
  // the solutions of all states on the region and interface (zero on the environment
  // interior; the full fields without a saved model), and volume_energy their energies in
  // the volume, (E_kk - (u_k - a_k)^T M_sheet (u_k - a_k)) / 2, which normalize the error
  // indicators as in the regular solve.
  mfem::DenseMatrix MagnetostaticEnergyMatrix(
      const std::vector<int> &current_ids, const std::vector<Vector> &J,
      const std::vector<int> &flux_ids, const std::vector<Vector> &a,
      const std::vector<Vector> &apertures, mfem::DenseMatrix *linked_flux,
      std::vector<Vector> *fields = nullptr, int n_fields = 0,
      std::vector<Vector> *region_fields = nullptr,
      std::vector<double> *volume_energy = nullptr);

  // The fluxoids c_j^T a_k of the flux-loop excitations k (ids) on the loops j (loop_ids),
  // measured on the mesh: saved with the model offline, and on the mesh of the region taken
  // from the model on the loops it leaves out. Collective.
  void EnvironmentFluxoids(const std::vector<int> &ids, const std::vector<int> &loop_ids,
                           mfem::DenseMatrix &fluxoid) const;

  // Online, the excitations (index -> loops) of a single flux loop of the environment whose
  // modes the model holds, if the excitations of the environment are those of the model
  // (always on the mesh of the region): having no part in the region, they need no fluxoid
  // generator. None if an excitation drives several loops, whose fluxoids are reported.
  // After CondenseEnvironment.
  std::vector<int> ModelExcitations(const std::map<int, std::vector<int>> &excitations);

  // Error indicators for adaptive refinement of the region only: Palace's gradient-flux
  // (electrostatics) or curl-flux (magnetostatics) estimator on the region solutions of
  // CapacitanceMatrix or MagnetostaticEnergyMatrix, normalized by the given energies (as in
  // the regular solve, the energies of the fields in the volume), and zero on the
  // environment.
  ErrorIndicator RegionErrorIndicator(const std::vector<Vector> &region_fields,
                                      const std::vector<double> &energies) const;

  // Terminal indices (sorted).
  std::vector<int> TerminalIndices() const;

  // Global parent true-DOF size.
  long long int GlobalTrueVSize() const;

  // The native operator of the problem (finest mesh and space only), for its postprocessing
  // and, in magnetostatics, its excitations.
  LaplaceOperator &GetLaplaceOperator();
  CurlCurlOperator &GetCurlCurlOperator();

  // Whether the configuration asks for the postprocessing of the fields (domain energies,
  // probes, surface fluxes or interface dielectrics) or saves fields.
  static bool Postprocesses(const IoData &iodata);

  // The configuration with the postprocessing that the field of the region gives: the
  // energies of its domains, the probes in it, and the surface fluxes and interface
  // dielectrics on surfaces in it (on the mesh of the region, none next to the environment
  // on the full mesh, which is cut off there). The rest is left out with a warning.
  // Collective.
  IoData RegionPostprocessing(const IoData &iodata) const;

  // The finite element spaces of the regular operator (for the problem size).
  const FiniteElementSpaceHierarchy &GetFESpaces() const;

  // Verification: the general solves and energies against which the tests check the
  // condensation (the simulations use the above).

  // Solve for the given terminal driven to 1 V, all other terminals grounded. Returns the
  // full parent field.
  Vector SolveExcitation(int drive_terminal_index);

  // Solve A u = f for a full parent-space source f. Returns the full parent field.
  Vector SolveSource(const Vector &f);

  // Mutual energy u_i^T K u_j with the full (region + environment) energy operator.
  double MutualEnergy(const Vector &ui, const Vector &uj) const;

  // Mutual energy u_i^T A u_j with the full (region + environment) operator of the solves:
  // the energy operator and, in magnetostatics, the sheet terms and the regularization.
  double MutualSystemEnergy(const Vector &ui, const Vector &uj) const;

  // Unit Dirichlet lift of a terminal (1 on its true DOFs, 0 elsewhere).
  Vector TerminalLift(int terminal_index) const;

  // Whether the environment interior solver has been set up.
  bool EnvironmentFactored() const;

private:
  // Region-condensed solves for a batch of prescribed Dirichlet fields.
  std::vector<Vector> SolveDirichletBatch(const std::vector<Vector> &dbcs);

  // Region-condensed solves of A u = f for a batch of full parent-space sources f.
  std::vector<Vector> SolveSourceBatch(const std::vector<Vector> &f);

  // Energy matrix E_ij = u_i^T K u_j (K the energy operator) of Dirichlet-lift excitations,
  // with the Dirichlet data xd[k] on the Dirichlet true DOFs, labeled ids[k]. Once the lift
  // modes are known (from a saved model, matched by id and an environment-side
  // fingerprint), no environment solve is needed; otherwise they are computed and appended
  // to the model being saved. If fields is non-null, it receives the full fields of the
  // first n_fields lifts, and region_fields, if non-null, their solutions on the region and
  // interface.
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
  // it receives the source works b_k^T u_k; fun->values receives l_k^T u_j; kinetic
  // receives the kinetic parts (u_k - a_k)^T M_sheet (u_k - a_k) of the G_kk.
  mfem::DenseMatrix SourceEnergyMatrix(const std::vector<int> &ids,
                                       const std::vector<Source> &src,
                                       std::vector<Vector> *fields, int n_fields,
                                       std::vector<double> *work = nullptr,
                                       Functionals *fun = nullptr,
                                       std::vector<Vector> *region_fields = nullptr,
                                       std::vector<double> *kinetic = nullptr);

  struct Impl;
  std::unique_ptr<Impl> impl;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUBSTRUCTURING_SOLVER_HPP
