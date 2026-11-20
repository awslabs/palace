// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP
#define PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP

#include <array>
#include <complex>
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
// The terms at a frequency are those of the model's record there, or of the combination of
// its records, on the current interface, with the rows of Γ of each process held by it.
// They are formed when a frequency is first asked for, and held until another one is: the
// calls are collective, for the frequencies in turn.
//
class EnvironmentOperator
{
public:
  // Read the model and check it against the run: the environment (its record and, on the
  // full mesh, the fingerprint of its operator), the interface, the environment's ports,
  // the sources of the excitations and the frequencies of the sweep. Call before the
  // region is set up, on the operators of the full problem. Collective.
  EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op);
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

  MPI_Comm GetComm() const { return comm; }

private:
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
  // The records of a rational model (root).
  std::vector<std::vector<std::complex<double>>> records;
  mfem::Array<int> pinned;

  // The rows of Γ of this process: its local true DOFs (in interface order), and the
  // counts and offsets of all processes.
  std::vector<int> gamma_local;
  std::vector<int> gamma_cnt, gamma_disp;

  // The terms at omega_cur: the rows of Γ of this process of S_E (row-major, x nG) and of
  // the condensed sources g (per model column), and on root, the ports' h (per port) and
  // c (ports x model columns, column-major).
  double omega_cur = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::complex<double>> S_loc;
  std::vector<std::vector<std::complex<double>>> g_loc, h;
  std::vector<std::complex<double>> c;

  // Form the terms at ω, unless they are held.
  void Prepare(double omega);

  // The values of Γ of a vector on this process, in interface order.
  std::vector<std::complex<double>> LocalInterfaceValues(const ComplexVector &x) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP
