// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
#define PALACE_MODELS_WAVE_PORT_OPERATOR_HPP

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <unordered_map>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/modeeigensolver.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class SumVectorCoefficient;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;
class Units;

namespace config
{

struct BoundaryData;
struct DomainData;
struct LinearSolverData;
struct SolverData;
struct WavePortData;

}  // namespace config

enum class ProblemType : char;

class FiniteElementSpace;

//
// Wave port data for a single port boundary.
//
class WavePortData
{
public:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Wave port properties.
  int mode_idx;
  double d_offset;
  int excitation;
  bool active;
  std::complex<double> kn0;
  double omega0;
  mfem::Vector port_normal;

private:
  // List of all boundary attributes making up this port boundary.
  mfem::Array<int> attr_list;

  // SubMesh data structures to define finite element spaces and grid functions on the
  // SubMesh corresponding to this port boundary.
  std::unique_ptr<Mesh> port_mesh;
  std::unique_ptr<mfem::FiniteElementCollection> port_nd_fec, port_h1_fec;
  std::unique_ptr<FiniteElementSpace> port_nd_fespace, port_h1_fespace;
  std::unique_ptr<mfem::ParTransferMap> port_nd_transfer, port_h1_transfer;
  std::unordered_map<int, int> submesh_parent_elems;
  mfem::Array<int> port_dbc_tdof_list;
  double mu_eps_max;

  // Submesh-specific material operator and boundary condition operators. Constructed on
  // the remapped submesh with rebuilt CEED data for correct attribute mapping.
  std::unique_ptr<MaterialOperator> port_mat_op;
  std::unique_ptr<SurfaceImpedanceOperator> port_surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> port_farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> port_surf_sigma_op;

  // Boundary mode eigenvalue problem solver.
  std::unique_ptr<ModeEigenSolver> mode_solver;
  ComplexVector v0, e0;

  // Communicator for processes which have elements for this port.
  MPI_Comm port_comm = MPI_COMM_NULL;
  int port_root = 0;

  // Grid functions storing the last computed electric field mode on the port, and stored
  // objects for computing functions of the port modes for use as an excitation or in
  // postprocessing.
  std::unique_ptr<GridFunction> port_E0t, port_E0n, port_S0t, port_E;
  std::unique_ptr<mfem::LinearForm> port_sr, port_si;

  // Voltage path for line integral (optional, for impedance postprocessing).
  std::vector<mfem::Vector> voltage_path;
  int voltage_n_samples;
  bool has_voltage_coords = false;

  // Reverse transfer map (port submesh → parent mesh) and parent-mesh GridFunction
  // used to evaluate line integrals of the port mode field via GSLIB on the 3D parent
  // mesh. Only allocated if the user configured a voltage path.
  std::unique_ptr<mfem::ParTransferMap> port_nd_transfer_reverse;
  std::unique_ptr<GridFunction> parent_E0t;

  // Cached GSLIB point locator for the voltage-path line integrals, Setup once on the
  // (fixed) parent mesh. The geometric Setup is the dominant cost, so reusing it across
  // every Initialize() avoids an O(num_elements) hash rebuild per frequency — critical for
  // the circuit-synthesis dispersion fit, which evaluates the mode at many frequencies.
  // Only allocated if the user configured a voltage path.
#if defined(MFEM_USE_GSLIB)
  std::unique_ptr<mfem::FindPointsGSLIB> voltage_gslib_op;
#endif

  // Optional polarity attributes (parent-mesh boundary attrs [high, low], signal
  // first, ground second). When non-zero (i.e. set by the user) the mode is flipped
  // so that E points from high to low, matching the lumped-port `+R Direction`
  // convention. Computed without GSLIB by evaluating the mode at port-submesh
  // boundary element centroids adjacent to each attribute.
  std::array<int, 2> polarity_attributes = {0, 0};

public:
  // 3D submesh constructor: extracts submesh from parent mesh.
  WavePortData(const config::WavePortData &data, const config::BoundaryData &boundaries,
               const config::DomainData &domains, ProblemType problem_type,
               const config::LinearSolverData &linear, const Units &units,
               const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
               mfem::ParFiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_attr);

  ~WavePortData();

  [[nodiscard]] constexpr bool HasExcitation() const { return excitation != 0; }
  [[nodiscard]] bool HasVoltageCoords() const { return has_voltage_coords; }

  const auto &GetAttrList() const { return attr_list; }

  void Initialize(double omega);

  // Compute the sign of the modal E-field projected on the (high → low) direction
  // implied by the given pair of parent-mesh boundary attributes (signal terminal
  // first, ground terminal second). Returns +1, -1, or 0 if attributes were not
  // resolvable on this mesh partition. Used by Initialize() to flip the mode polarity
  // without requiring GSLIB.
  [[nodiscard]] int GetModePolaritySign(int high_attr, int low_attr) const;

  HYPRE_BigInt GlobalTrueNDSize() const { return port_nd_fespace->GlobalTrueVSize(); }
  HYPRE_BigInt GlobalTrueH1Size() const { return port_h1_fespace->GlobalTrueVSize(); }

  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientReal() const;
  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientImag() const;

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeFieldCoefficientReal(double scaling = 1.0) const;
  std::unique_ptr<mfem::VectorCoefficient>
  GetModeFieldCoefficientImag(double scaling = 1.0) const;

  // Characteristic impedance Z = |V|^2 / (2P) from the port mode voltage and unit power.
  // Requires voltage coordinates to be configured.
  std::complex<double> GetCharacteristicImpedance() const;

  double GetExcitationPower() const;

  // Excitation voltage from the normalized port mode field.
  std::complex<double> GetExcitationVoltage() const;

  std::complex<double> GetPower(GridFunction &E, GridFunction &B) const;
  std::complex<double> GetSParameter(GridFunction &E) const;

  // Voltage line integral V = integral of E . dl on the port, using the 3D E field.
  // Requires voltage coordinates to be configured.
  std::complex<double> GetVoltage(GridFunction &E) const;
};

//
// A class handling wave port boundaries and their postprocessing.
//
class WavePortOperator
{
private:
  // Mapping from port index to data structure containing port information.
  std::map<int, WavePortData> ports;

  // Flag which forces no printing during WavePortData::Print().
  bool suppress_output;
  double fc, kc;

  void SetUpBoundaryProperties(const config::BoundaryData &boundaries,
                               const config::DomainData &domains,
                               const config::SolverData &solver, ProblemType problem_type,
                               const Units &units, const MaterialOperator &mat_op,
                               mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

protected:
  // Compute boundary modes for all wave port boundaries at the specified frequency.
  // Protected so test fixtures can drive this directly via a thin subclass.
  void Initialize(double omega);

public:
  WavePortOperator(const config::BoundaryData &boundaries,
                   const config::DomainData &domains, const config::SolverData &solver,
                   ProblemType problem_type, const Units &units,
                   const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace);
  WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                   mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace);

  // Access data structures for the wave port with the given index.
  const WavePortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Enable or suppress all outputs (log printing and fields to disk).
  void SetSuppressOutput(bool suppress) { suppress_output = suppress; }

  // Returns array of wave port attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix from wave ports.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Add the ω-independent boundary mass contribution from a single wave port to a material
  // property coefficient. The full system contribution from this port is
  // i·kₙ(ω)·(this) — see AddExtraSystemBdrCoefficients. Used to factor out the
  // ω-independent operator for reuse by the reduced-order model. The port_idx must refer
  // to an active port; if it does not, no contribution is added.
  void AddBoundaryMassBdrCoefficients(int port_idx, MaterialPropertyCoefficient &fb) const;

  // As above, with an explicit scale factor on the contribution. Used internally by
  // AddExtraSystemBdrCoefficients to share the per-port coefficient setup without
  // duplicating the mfem::MaterialPropertyCoefficient construction.
  void AddBoundaryMassBdrCoefficients(int port_idx, MaterialPropertyCoefficient &fb,
                                      double scale) const;

  // Compute and return the (real part of the) modal propagation constant kₙ for the
  // specified port at the given operating frequency. Triggers the per-port cross-section
  // EVP via WavePortData::Initialize; result is cached per port (per ω). Used by the
  // reduced-order model to assemble the wave-port contribution online without touching
  // any HDM-size object.
  double GetWavePortKn(int port_idx, double omega);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries.
  void AddExcitationBdrCoefficients(int excitation_idx, double omega,
                                    SumVectorCoefficient &fbr, SumVectorCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
