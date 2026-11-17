// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP
#define PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/substructure.hpp"
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class DrivenSubstructure;
struct DrivenSubstructureModel;
class IoData;
class SpaceOperator;

//
// The environment of an online driven substructuring run, condensed onto the interface Γ
// by an offline run (DrivenSubstructureModel), for the regular driven solver to solve the
// region against it: at a frequency of the model, the environment operator S_E(ω), the
// condensation g_k(ω) of the environment's sources of each excitation and the terms of the
// environment's port voltages V_jk = c_jk + h_j^T u_Γ,k, on the current interface
// (SignatureMap). With SpaceOperator::SetCondensedEnvironment, the system is the region's
// operator (assembled under a restriction to the region) with S_E(ω) added on Γ and the
// identity on the true DOFs outside the region and Γ (the environment interior and the
// Dirichlet DOFs of environment elements only), where the solution is zero, and the
// excitation is the region's with g_k(ω) on Γ. The preconditioner of that system
// (GetPreconditioner) is the coupled factorization of the online sweep: the region interior
// and Γ with S_E(ω), factored as one complex symmetric system by DrivenSubstructure, its
// pattern refilled at each frequency.
//
// This is a prototype for the review of the online path of DrivenSubstructuredSweep: exact
// models (the frequencies of the model only) and lumped ports.
//
class CondensedEnvironment
{
public:
  // Read the model and check it against this run, as the online sweep does. Construct
  // before SpaceOperator::SetCondensedEnvironment, which changes the operators and
  // excitations the checks use. Collective.
  CondensedEnvironment(const IoData &iodata, SpaceOperator &space_op);
  ~CondensedEnvironment();

  // Whether an online driven substructuring run takes this path: the prototype is chosen
  // by PALACE_SUBSTRUCTURING_DRIVENSOLVER=1 in the environment (else, or for an adaptive
  // sweep, DrivenSubstructuredSweep).
  static bool Enabled(const config::SolverData &solver);

  // The preconditioner of the region's system: solves with the coupled factorization of the
  // region interior and Γ for the entries of a vector there, zero elsewhere (where the
  // solution is zero). Its SetOperator factors the region against the S_E(ω) read last,
  // that of the system operator (GetExtraSystemOperator(ω) reads it), and ignores the
  // operator otherwise.
  std::unique_ptr<Solver<ComplexOperator>> GetPreconditioner();

  const std::vector<int> &RegionAttributes() const;

  // The essential true DOFs of the region's operators (on the finest space): the Dirichlet
  // DOFs of region elements. The other true DOFs outside the region and Γ have the identity
  // added with S_E (GetSchurOperator), since the operators of the region store no diagonal
  // entry for them: the system operator is nonsingular, with a zero solution there.
  const mfem::Array<int> &EssentialTrueDofs() const { return ess_tdofs; }

  // At ω (a frequency of the model; its record is read when ω changes), S_E(ω) as an
  // operator on the true DOFs, with the identity on the pinned true DOFs which are not
  // essential. Collective.
  std::unique_ptr<ComplexOperator> GetSchurOperator(double omega);

  // The excitation of the region at ω from the full excitation vector: RHS on the region
  // interior, g_k(ω) on Γ and zero elsewhere. Collective.
  void ApplyExcitation(int excitation_idx, double omega, ComplexVector &RHS);

  // The voltages of the environment's lumped ports at ω, for the solution E of an
  // excitation (given to PostOperator, which cannot measure them on the field, zero in the
  // environment). Collective.
  std::map<int, std::complex<double>> PortVoltages(int excitation_idx, double omega,
                                                   const ComplexVector &E);

  // The postprocessing of the domains and boundaries restricted to the region (see
  // RestrictPostprocessingToRegion).
  const config::DomainData &GetDomains() const { return domains; }
  const config::BoundaryData &GetBoundaries() const { return boundaries; }

private:
  SpaceOperator &space_op;
  const config::SubstructuringData &sub;
  MPI_Comm comm;
  std::unique_ptr<DrivenSubstructureModel> model;
  std::unique_ptr<DrivenSubstructure> ds;
  SignatureMap gamma_map;
  int nG = 0;
  // The excitations of the run, and their columns in the model (-1: no environment
  // sources); the environment's lumped ports.
  std::vector<int> ex_idx, ex_col, env_ports;
  config::DomainData domains;
  config::BoundaryData boundaries;
  // Local true DOFs: the interface ones in interface order (a contiguous range of it from
  // gamma_off on each rank, with the counts and offsets of the ranks), and the region
  // interior; the pinned true DOFs which are not essential (with the identity in
  // GetSchurOperator).
  std::vector<int> gamma_local, gamma_cnt, gamma_disp;
  int gamma_off = 0;
  std::vector<char> region_int;
  mfem::Array<int> ess_tdofs, identity_tdofs;
  // The record read last (its index in the model, -1 if none): the rows of S_E of the local
  // interface DOFs (m x nG, row-major), and g (nG x excitations of the run), h (nG x
  // ports) on the current interface and c (ports x excitations of the model), replicated.
  int record = -1;
  std::shared_ptr<const std::vector<std::complex<double>>> S_rows;
  std::vector<std::complex<double>> g, h, c;
  // For the coupled factorization: the frequency of the record and its full S_E on rank 0
  // until factored, and the record factored last.
  double omega_record = 0.0;
  std::vector<std::complex<double>> S_root;
  int factored = -1;

  void Load(double omega);

  // Factor the region against the S_E of the record read last (if not yet).
  void Factor();

  class Preconditioner;
  class SchurOperator;
};

//
// Helpers of the model of driven substructuring, for its offline and online sweeps.
//

// The lumped ports on the environment side (and its wave ports, as -index), and an abort
// for lumped ports touching Γ (their excitation and voltage would depend on both sides).
std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::map<int, std::vector<char>> &port_dofs);

// The parameters of the environment's ports that enter S_E (see
// DrivenSubstructureModel::port_data).
std::vector<double> EnvironmentPortData(const IoData &iodata,
                                        const std::vector<int> &env_ports);

// The full S_E (n x n, column-major) of a record's lower triangle (by columns), mapped to
// the current interface as M^-T S_E M^-1 with map (rows[g]: row g of M^-T), on one
// triangle.
std::vector<std::complex<double>> FullSchur(const std::complex<double> *lower, int n,
                                            const SignatureMap *map);

// The index of x in y (relative tolerance tol), -1 if none.
int FindValue(const std::vector<double> &y, double x, double tol = 1.0e-10);

// From a model, the field is known in the region (and on Γ) only: the domain and boundary
// postprocessing restricted to energies of region domains, without totals or participation
// ratios, and to the surface fluxes and probes the region holds, with a warning of what is
// skipped (and that the environment's ports, if any, have no power). mode names the run in
// the warning. Collective.
void RestrictPostprocessingToRegion(const mfem::ParMesh &mesh,
                                    const config::SubstructuringData &sub,
                                    const DrivenSubstructure &ds,
                                    const DrivenSubstructureModel &model, bool env_ports,
                                    const char *mode, config::DomainData &domains,
                                    config::BoundaryData &boundaries);

}  // namespace palace

#endif  // PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP
