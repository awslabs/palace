// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP
#define PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP

#include <array>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/substructure.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/drivensubstructure.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

//
// The environment of an online driven substructuring run, condensed onto the interface Γ by
// its saved model (see DrivenSubstructureModel): its Schur complement S_E(ω), the
// condensation g_E(ω) of the sources of the excitations, and the values of its ports, at
// the frequencies of an exact model or anywhere between the samples of a rational one. The
// region is then the driven problem on its elements with S_E(ω) added on Γ (see
// SpaceOperator::SetEnvironment).
//
// The terms at a frequency are those of the model's record there, or of the combination
// Σ_k c_k(ω) of its records with the barycentric coefficients of a rational model, on the
// current interface, each process holding its rows of Γ. They are formed when a frequency
// is first asked for, and held until another one is: the calls are collective, for the
// frequencies in turn. The records of a rational model are held whole, for the terms of a
// reduced model of the region (see RomOperator).
//
class EnvironmentOperator
{
public:
  // The records of a model, in order (flat, in the saved basis, on root).
  using RecordFunction = std::function<std::vector<std::complex<double>>(std::size_t)>;

  // Read the model and check it against the run: the environment (its record and, on the
  // full mesh, the fingerprint of its operator), the interface, the environment's ports,
  // the sources of the excitations and the frequencies of the sweep. Call before the
  // region is set up, on the operators of the full problem. Collective.
  EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op);

  // The rational model of this run's environment, sampled offline on the full mesh, with
  // its records. Collective.
  EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op,
                      const DrivenSubstructureModel &sampled, const RecordFunction &record);
  ~EnvironmentOperator();

  // Whether the mesh is the mesh of the region written with the model.
  bool RegionMesh() const { return region_mesh; }

  // The domain attributes of the region.
  const std::vector<int> &RegionAttributes() const { return region_attrs; }

  // The true DOFs the region pins, those of the environment off Γ: its interior and, on
  // the mesh of the region, the other DOFs of the environment elements it keeps.
  const mfem::Array<int> &PinnedTrueDofs() const { return pinned; }

  // The environment's lumped ports (index) and wave ports (-index).
  const std::vector<int> &Ports() const { return env_ports; }

  // The boundary attributes in the region the model records (see
  // DrivenSubstructureModel::region_bdr_attributes).
  const std::vector<int> &RegionBoundaryAttributes() const
  {
    return model.region_bdr_attributes;
  }

  // y += R_Γᵀ S_E(ω) R_Γ x.
  void AddMult(double omega, const ComplexVector &x, ComplexVector &y);

  // The real and imaginary parts of R_Γᵀ S_E(ω) R_Γ, assembled on the true DOFs.
  std::array<std::unique_ptr<mfem::HypreParMatrix>, 2> Assemble(double omega);

  // RHS += R_Γᵀ g_E(ω) of an excitation; false if it has no sources in the environment.
  bool AddExcitation(int excitation_idx, double omega, ComplexVector &RHS);

  // The environment's lumped port voltages (V) and wave port overlaps (S) for a solution E
  // of an excitation at ω.
  void PortValues(int excitation_idx, double omega, const ComplexVector &E,
                  std::map<int, std::complex<double>> &V,
                  std::map<int, std::complex<double>> &S);

  // A rational model: its number of records, their coefficients c_k(ω), y += R_Γᵀ S_k R_Γ x
  // for record k, and v = R_Γᵀ g_k of an excitation (false if it has no sources in the
  // environment).
  std::size_t NumRecords() const { return records.size(); }
  std::vector<std::complex<double>> Coefficients(double omega) const;
  void RecordAddMult(std::size_t k, const ComplexVector &x, ComplexVector &y) const;

  // For a real vector: y += the real or imaginary part of R_Γᵀ S_k R_Γ x, and the pairings
  // v_Γᵀ g_k of v with the condensed sources of an excitation of every record k (empty if
  // it has no sources in the environment).
  void RecordAddMultPart(std::size_t k, bool imag, const Vector &x, Vector &y) const;
  std::vector<std::complex<double>> RecordExcitationPairings(int excitation_idx,
                                                             const Vector &v) const;

  MPI_Comm GetComm() const { return comm; }

private:
  EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op,
                      const DrivenSubstructureModel *sampled, const RecordFunction &record);

  SpaceOperator &space_op;
  MPI_Comm comm;
  const bool root;
  const std::string model_path;
  std::vector<int> region_attrs;
  DrivenSubstructureModel model;
  std::unique_ptr<DrivenSubstructure> ds;
  bool region_mesh = false;
  int nG = 0, np = 0, ne = 0;
  std::vector<int> env_ports;
  SignatureMap gamma_map;  // the saved interface on the current one
  // The model column of each excitation with environment sources (absent: none).
  std::map<int, int> ex_col;
  mfem::Array<int> pinned;

  // The rows of Γ of this process: its local true DOFs (in interface order), and the
  // counts and offsets of all processes.
  std::vector<int> gamma_local;
  std::vector<int> gamma_cnt, gamma_disp;

  // A record on the current interface, as the rows of Γ of this process: of S_E (x nG,
  // row-major), of the condensed sources g (per model column) and of the port functionals h
  // (per port), with c (ports x model columns, column-major) on every process.
  struct Rows
  {
    std::vector<std::complex<double>> S;
    std::vector<std::vector<std::complex<double>>> g, h;
    std::vector<std::complex<double>> c;
  };

  // The records of a rational model, and the terms at omega_cur.
  std::vector<Rows> records;
  double omega_cur = std::numeric_limits<double>::quiet_NaN();
  Rows cur;

  // The rows of a flat record of the model (on root, in the saved basis). Collective.
  Rows Distribute(const std::vector<std::complex<double>> &flat) const;

  // Form the terms at ω, unless they are held.
  void Prepare(double omega);

  // The values of Γ of a vector on this process, in interface order, and the product of
  // the rows of a dense block with the interface values of x (all of Γ) added to y.
  std::vector<std::complex<double>> LocalInterfaceValues(const ComplexVector &x) const;
  void AddRows(const std::vector<std::complex<double>> &rows, const ComplexVector &x,
               ComplexVector &y) const;
  std::vector<std::complex<double>> AllInterfaceValues(const Vector &x) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP
