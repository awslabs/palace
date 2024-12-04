// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_POST_OPERATOR_HPP
#define PALACE_MODELS_POST_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/domainpostoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/filesystem.hpp"

namespace palace
{

class CurlCurlOperator;
class ErrorIndicator;
class IoData;
class LaplaceOperator;
class MaterialOperator;
class SpaceOperator;
class SurfaceCurrentOperator;
class WavePortOperator;

fs::path ParaviewPath(const IoData &iodata, size_t excitation_idx = 0,
                      size_t max_excitation = 1);

//
// A class to handle solution postprocessing.
//
class PostOperator
{
private:
  // Reference to material property operator (not owned).
  std::reference_wrapper<const MaterialOperator> mat_op;
  // Reference to mesh for paraview (re)-init
  std::reference_wrapper<mfem::ParMesh> mesh_ND;

  // Objects for grid function postprocessing from the FE solution.
  std::unique_ptr<GridFunction> E, B, V, A;
  std::unique_ptr<mfem::VectorCoefficient> S, E_sr, E_si, B_sr, B_si, A_s, J_sr, J_si;
  std::unique_ptr<mfem::Coefficient> U_e, U_m, V_s, Q_sr, Q_si;

  // Data collection for writing fields to disk for visualization
  // Paraview fields are mutable as the writing is triggered by const solver printer
  std::string name_;
  mutable mfem::ParaViewDataCollection paraview, paraview_bdr;
  double mesh_Lc0;

  // Set paraview options and register fields
  void InitializeParaviewDataCollection();

  // ----- Measurements from Fields -----

  DomainPostOperator dom_post_op;           // Energy in bulk
  SurfacePostOperator surf_post_op;         //  Dielectric Interface Energy and Flux
  mutable InterpolationOperator interp_op;  // E & B fields: mutates during measure

  // Port Contributions: not owned, view onto space_op only, it must not go out of scope
  LumpedPortOperator *lumped_port_op = nullptr;
  WavePortOperator *wave_port_op = nullptr;

  // Wave port boundary mode field postprocessing.
  struct WavePortFieldData
  {
    std::unique_ptr<mfem::VectorCoefficient> E0r, E0i;
  };
  std::map<int, WavePortFieldData> port_E0;

public:
  // Mini storage clases for cache: Definitions are public
  struct FluxData
  {
    int idx;                   // Surface index
    std::complex<double> Phi;  // Integrated flux
    SurfaceFluxType type;      // Flux type
  };

  struct InterfaceData
  {
    int idx;          // Interface index
    double energy;    // Surface ELectric Field Energy
    double tandelta;  // Dissipation tangent tan(δ)
  };

  // For both lumped and wave port
  struct PortPostData
  {
    std::complex<double> P, V, S;
    std::array<std::complex<double>, 3> I;  // Separate R, L, and C branches
  };

private:
  struct MeasurementCache
  {
    std::optional<std::complex<double>> omega = std::nullopt;

    std::optional<double> domain_E_field_energy_all = std::nullopt;
    std::optional<double> domain_H_field_energy_all = std::nullopt;

    std::optional<std::map<int, double>> domain_E_field_energy_i = std::nullopt;
    std::optional<std::map<int, double>> domain_H_field_energy_i = std::nullopt;

    std::optional<std::vector<FluxData>> surface_flux_i = std::nullopt;
    std::optional<std::vector<InterfaceData>> interface_eps_i = std::nullopt;

    std::optional<std::map<int, PortPostData>> lumped_port_vi = std::nullopt;
    std::optional<std::map<int, PortPostData>> wave_port_vi = std::nullopt;

    std::optional<double> lumped_port_inductor_energy = std::nullopt;
    std::optional<double> lumped_port_capacitor_energy = std::nullopt;

    std::optional<std::vector<std::complex<double>>> probe_E_field = std::nullopt;
    std::optional<std::vector<std::complex<double>>> probe_B_field = std::nullopt;
  };
  mutable MeasurementCache measurment_cache = {};

  void ValidateDoPortMeasurement() const;

public:
  PostOperator(const IoData &iodata, SpaceOperator &space_op, std::string name,
               const std::optional<fs::path> &paraview_path = {});
  PostOperator(const IoData &iodata, LaplaceOperator &laplace_op, std::string name,
               const std::optional<fs::path> &paraview_path = {});
  PostOperator(const IoData &iodata, CurlCurlOperator &curlcurl_op, std::string name,
               const std::optional<fs::path> &paraview_path = {});

  // Make new paraview output target for multiple excitations / prom
  // Previous output must be correctly closed / deregistered
  void SetNewParaviewOutput(const fs::path &paraview_path);

  // Access to surface and domain postprocessing objects.
  const auto &GetSurfacePostOp() const { return surf_post_op; }
  const auto &GetDomainPostOp() const { return dom_post_op; }

  // Return options for postprocessing configuration.
  bool HasE() const { return E != nullptr; }
  bool HasB() const { return B != nullptr; }
  bool HasImag() const { return HasE() ? E->HasImag() : B->HasImag(); }

  // Populate the grid function solutions for the E- and B-field using the solution vectors
  // on the true dofs. For the real-valued overload, the electric scalar potential can be
  // specified too for electrostatic simulations. The output mesh and fields are
  // nondimensionalized consistently (B ~ E (L₀ ω₀ E₀⁻¹)).
  void SetEGridFunction(const ComplexVector &e, bool exchange_face_nbr_data = true);
  void SetBGridFunction(const ComplexVector &b, bool exchange_face_nbr_data = true);
  void SetEGridFunction(const Vector &e, bool exchange_face_nbr_data = true);
  void SetBGridFunction(const Vector &b, bool exchange_face_nbr_data = true);
  void SetVGridFunction(const Vector &v, bool exchange_face_nbr_data = true);
  void SetAGridFunction(const Vector &a, bool exchange_face_nbr_data = true);

  // Access grid functions for field solutions.
  auto &GetEGridFunction()
  {
    MFEM_ASSERT(E, "Missing GridFunction object when accessing from PostOperator!");
    return *E;
  }
  auto &GetBGridFunction()
  {
    MFEM_ASSERT(B, "Missing GridFunction object when accessing from PostOperator!");
    return *B;
  }
  auto &GetVGridFunction()
  {
    MFEM_ASSERT(V, "Missing GridFunction object when accessing from PostOperator!");
    return *V;
  }
  auto &GetAGridFunction()
  {
    MFEM_ASSERT(A, "Missing GridFunction object when accessing from PostOperator!");
    return *A;
  }

  // Function that triggers all available post-processing measurements and populate cache.
  void MeasureAll();

  // Clear internal measurement caches
  void ClearAllMeasurementCache();

  // Treat the frequency, for driven and eigemode solvers, as a "measurement", that other
  // measurements can depend on. This has to be supplied during the solver loop separate
  // from the fields.
  void SetFrequency(double omega);
  void SetFrequency(std::complex<double> omega);

  // Return stored frequency that was given in SetFrequency. Always promotes to complex
  // frequency.
  std::complex<double> GetFrequency() const;

  // Postprocess the total electric and magnetic field energies in the electric and magnetic
  // fields.
  double GetEFieldEnergy() const;
  double GetHFieldEnergy() const;

  // Postprocess the electric and magnetic field energies in the domain with the given
  // index.
  double GetEFieldEnergy(int idx) const;
  double GetHFieldEnergy(int idx) const;

  // Postprocess the electric or magnetic field flux for a surface index using the computed
  // electcric field and/or magnetic flux density field solutions.
  std::vector<FluxData> GetSurfaceFluxAll() const;
  FluxData GetSurfaceFlux(int idx) const;

  // Postprocess the partitipation ratio for interface lossy dielectric losses in the
  // electric field mode.
  double GetInterfaceParticipation(int idx, double E_m) const;
  std::vector<InterfaceData> GetInterfaceEFieldEnergyAll() const;
  InterfaceData GetInterfaceEFieldEnergy(int idx) const;

  // Measure and cache port voltages and currents for lumped and wave port operators.
  void MeasureLumpedPorts() const;
  void MeasureWavePorts() const;

  // Postprocess the energy in lumped capacitor or inductor port boundaries with index in
  // the provided set.
  double GetLumpedInductorEnergy() const;
  double GetLumpedCapacitorEnergy() const;

  // Postprocess the S-parameter for recieving lumped or wave port index using the electric
  // field solution.
  std::complex<double> GetSParameter(bool is_lumped_port, int idx, int source_idx) const;

  // Postprocess the circuit voltage and current across lumped port index using the electric
  // field solution. When the internal grid functions are real-valued, the returned voltage
  // has only a nonzero real part.
  std::complex<double> GetPortPower(int idx) const;
  std::complex<double> GetPortVoltage(int idx) const;
  std::complex<double>
  GetPortCurrent(int idx,
                 LumpedPortData::Branch branch = LumpedPortData::Branch::TOTAL) const;

  // Postprocess the EPR for the electric field solution and lumped port index.
  double GetInductorParticipation(int idx, double E_m) const;

  // Postprocess the coupling rate for radiative loss to the given I-O port index.
  double GetExternalKappa(int idx, double E_m) const;

  // Write to disk the E- and B-fields extracted from the solution vectors. Note that fields
  // are not redimensionalized, to do so one needs to compute: B <= B * (μ₀ H₀), E <= E *
  // (Z₀ H₀), V <= V * (Z₀ H₀ L₀), etc.
  void WriteFields(int step, double time) const;
  void WriteFieldsFinal(const ErrorIndicator *indicator = nullptr) const;

  // Probe the E- and B-fields for their vector-values at speceified locations in space.
  // Locations of probes are set up in constructor from configuration file data. If
  // the internal grid functions are real-valued, the returned fields have only nonzero real
  // parts. Output vectors are ordered by vector dimension, that is [v1x, v1y, v1z, v2x,
  // v2y, v2z, ...].
  int GetInterpolationOpVDim() const { return interp_op.GetVDim(); }
  const auto &GetProbes() const { return interp_op.GetProbes(); }
  std::vector<std::complex<double>> ProbeEField() const;
  std::vector<std::complex<double>> ProbeBField() const;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const
  {
    return (E) ? E->ParFESpace()->GetComm() : B->ParFESpace()->GetComm();
  }
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_HPP
