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
#include <optional>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/substructure.hpp"
#include "linalg/blockcyclic.hpp"
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
// region is then the driven problem on its elements with S_E(ω) added on Γ (see the
// SpaceOperator of an environment).
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
  // Record j of a model (flat, in the saved basis, on root): called on every process, in
  // turn. Collective.
  using RecordFunction = std::function<std::vector<std::complex<double>>(std::size_t)>;

  // Read the model and check it against the run: the environment (its record and, on the
  // full mesh, the fingerprint of its operator), the interface, the environment's ports,
  // the sources of the excitations and the frequencies of the sweep. Call before the
  // region is set up, on the operators of the full problem. grid_procs: the number of
  // processes of the grid of the dense matrices on Γ (by default, see BlockCyclicGrid).
  // Collective.
  EnvironmentOperator(const IoData &iodata, SpaceOperator &space_op, int grid_procs = 0);

  // A model of this run's environment, condensed offline on the full mesh: rational, with
  // its samples as records, taken here, or exact, at the frequencies of the sweep, whose
  // records are taken when the sweep reaches them. Collective.
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

  // On the mesh of the region, the boundary attributes of the environment elements it holds
  // (see EnvironmentRecord::bdr_attributes).
  const std::vector<int> &EnvironmentBoundaryAttributes() const
  {
    return model.environment.bdr_attributes;
  }

  // y += R_Γᵀ S_E(ω) R_Γ x.
  void AddMult(double omega, const ComplexVector &x, ComplexVector &y);

  // The dense matrices on Γ (S_E and its records) are 2D block cyclic on this grid, that of
  // the Schur complement of the region (see RegionSolver).
  const BlockCyclicGrid &Grid() const { return *grid; }

  // C += S_E(ω) for the local part C of a matrix on Grid().
  void AddSchurBlocks(double omega, std::vector<std::complex<double>> &C);

  // The interface Γ of the run, in its order (see SubstructureInterface).
  const SubstructureInterface &Interface() const { return iface; }

  // RHS += R_Γᵀ g_E(ω) of an excitation; false if it has no sources in the environment.
  bool AddExcitation(int excitation_idx, double omega, ComplexVector &RHS);

  // The environment's lumped port voltages (V) and wave port overlaps (S) for a solution E
  // of an excitation at ω, or for its values on Γ on this process (in interface order, as
  // LocalInterfaceValues gives them).
  void PortValues(int excitation_idx, double omega, const ComplexVector &E,
                  std::map<int, std::complex<double>> &V,
                  std::map<int, std::complex<double>> &S);
  void PortValues(int excitation_idx, double omega,
                  const std::vector<std::complex<double>> &e_gamma,
                  std::map<int, std::complex<double>> &V,
                  std::map<int, std::complex<double>> &S);

  // The values of a real vector on Γ on this process, in interface order.
  std::vector<double> LocalInterfaceValues(const Vector &x) const;

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
                      const DrivenSubstructureModel *sampled, const RecordFunction &record,
                      int grid_procs);

  MPI_Comm comm;
  const bool root;
  const std::string model_path;
  std::vector<int> region_attrs;
  DrivenSubstructureModel model;
  // The records of an exact model of this run (see the constructor).
  RecordFunction exact_records;
  bool region_mesh = false;
  int nG = 0, np = 0, ne = 0;
  std::vector<int> env_ports;
  SignatureMap gamma_map;  // the saved interface on the current one
  // The current interface DOFs each saved one adds to, with M^-T's coefficients, and the
  // positions of the records in the model file.
  std::vector<std::vector<std::pair<int, double>>> saved_to_current;
  std::vector<long long> record_pos;
  // The model column of each excitation with environment sources (absent: none).
  std::map<int, int> ex_col;
  mfem::Array<int> pinned;

  // The interface, and the local true DOFs of Γ on this process, in interface order.
  SubstructureInterface iface;
  std::vector<int> gamma_local;

  // A record on the current interface: the local part of S_E on the grid, the rows of Γ of
  // this process of the condensed sources g (per model column) and of the port functionals
  // h (per port), and c (ports x model columns, column-major) on every process.
  struct Rows
  {
    std::vector<std::complex<double>> S;
    std::vector<std::vector<std::complex<double>>> g, h;
    std::vector<std::complex<double>> c;
  };

  // The records of a rational model, and the terms at ω: S_E and g at omega_cur, c and h at
  // omega_fun.
  std::vector<Rows> records;
  double omega_cur = std::numeric_limits<double>::quiet_NaN();
  double omega_fun = std::numeric_limits<double>::quiet_NaN();
  Rows cur;

  // The grid of the dense matrices on Γ, and the global rows and columns of the local part.
  std::optional<BlockCyclicGrid> grid;
  std::vector<int> grid_rows, grid_cols;

  // A record from its flat entries on root (in the saved basis, see
  // DrivenSubstructureModel::Flatten), or record j of the model file, whose S_E each
  // process reads a share of, and the record from the parts after S_E on root. Collective.
  Rows Distribute(const std::vector<std::complex<double>> &flat) const;
  Rows ReadRecord(int j) const;
  Rows DistributeVectors(const std::complex<double> *tail) const;

  // S_E of the record whose payload starts at byte pos of the model file, on the grid: each
  // process reads an equal share of its lower triangle, and sends each entry, mapped onto
  // the current interface, to the processes holding the entries it adds to. Collective.
  std::vector<std::complex<double>> ReadSchur(long long pos) const;

  // Form the terms of the operator (S_E and g) at ω, or the port functionals (c and h),
  // unless they are held (an exact model reads them whole).
  void Prepare(double omega);
  void PrepareFunctionals(double omega);

  // The values of Γ of a vector on this process, in interface order, and those of all of
  // Γ (on every process).
  std::vector<std::complex<double>> LocalInterfaceValues(const ComplexVector &x) const;
  std::vector<std::complex<double>> AllInterfaceValues(const ComplexVector &x) const;
  std::vector<double> AllInterfaceValues(const Vector &x) const;

  // y = S x on all of Γ (on every process), for a matrix S on the grid (its local part), or
  // its real or imaginary part (for a real x). Collective.
  std::vector<std::complex<double>>
  Product(const std::vector<std::complex<double>> &S,
          const std::vector<std::complex<double>> &x) const;
  std::vector<double> ProductPart(const std::vector<std::complex<double>> &S, bool imag,
                                  const std::vector<double> &x) const;
  // y += R_Γᵀ S R_Γ x for a matrix S on the grid (its local part).
  void RecordProduct(const std::vector<std::complex<double>> &S, const ComplexVector &x,
                     ComplexVector &y) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ENVIRONMENT_OPERATOR_HPP
