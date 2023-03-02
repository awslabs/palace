// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_POST_OPERATOR_HPP
#define PALACE_POST_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/domainpostoperator.hpp"
#include "fem/interpolationoperator.hpp"
#include "fem/surfacepostoperator.hpp"

namespace palace
{

class CurlCurlOperator;
class IoData;
class LaplaceOperator;
class LumpedPortOperator;
class MaterialOperator;
class SpaceOperator;
class SurfaceCurrentOperator;
class WavePortOperator;

namespace petsc
{

class PetscParMatrix;
class PetscParVector;

}  // namespace petsc

//
// A class to handle solution postprocessing.
//
class PostOperator
{
private:
  // Shared face mapping for boundary coefficients.
  std::map<int, int> local_to_shared;

  // Reference to material property operator (not owned).
  const MaterialOperator &mat_op;

  // Surface boundary and domain postprocessors.
  const SurfacePostOperator surf_post_op;
  const DomainPostOperator dom_post_op;

  // Objects for grid function postprocessing from the FE solution.
  const bool has_imaginary;
  std::optional<mfem::ParComplexGridFunction> E, B;
  std::optional<mfem::ParGridFunction> V, A;
  std::unique_ptr<mfem::VectorCoefficient> Esr, Esi, Bsr, Bsi, As, Jsr, Jsi;
  std::unique_ptr<mfem::Coefficient> Vs, Ue, Um, Qsr, Qsi;

  // Lumped and wave port voltage and current (R, L, and C branches) caches updated when
  // the grid functions are set.
  struct PortPostData
  {
    std::complex<double> S, P, V, Z;
  };
  std::map<int, PortPostData> lumped_port_vi, wave_port_vi;
  bool lumped_port_init, wave_port_init;

  // Data collection for writing fields to disk for visualization and sampling points.
  mutable mfem::ParaViewDataCollection paraview, paraview_bdr;
  mutable InterpolationOperator interp_op;
  void InitializeDataCollection(const IoData &iodata);

public:
  PostOperator(const IoData &iodata, SpaceOperator &spaceop, const std::string &name);
  PostOperator(const IoData &iodata, LaplaceOperator &laplaceop, const std::string &name);
  PostOperator(const IoData &iodata, CurlCurlOperator &curlcurlop, const std::string &name);

  // Access to surface and domain postprocessing objects.
  const auto &GetSurfacePostOp() const { return surf_post_op; }
  const auto &GetDomainPostOp() const { return dom_post_op; }

  // Return options for postprocessing configuration.
  bool HasImaginary() const { return has_imaginary; }
  bool HasE() const { return E.has_value(); }
  bool HasB() const { return B.has_value(); }

  // Accessors for grid function representations. Performs check for initialization.
  const auto& GetE() const;
  const auto& GetB() const;
  const auto& GetV() const;
  const auto& GetA() const;

  // Compute the magnetic flux density B in RT space from electric field solution E solution
  // in ND space for the time-harmonic case: B =  -1/(iω) ∇ x E.
  static void GetBField(std::complex<double> omega, const petsc::PetscParMatrix &NegCurl,
                        const petsc::PetscParVector &e, petsc::PetscParVector &b);

  // Compute the magnetic flux density B in RT space from the magnetic vector potential
  // solution A in ND space: B = ∇ x A.
  static void GetBField(const mfem::Operator &Curl, const mfem::Vector &a, mfem::Vector &b);

  // Compute the electric field E in ND space from the scalar potential solution V in H1
  // space: E = -∇V.
  static void GetEField(const mfem::Operator &NegGrad, const mfem::Vector &v,
                        mfem::Vector &e);

  // Populate the grid function solutions for the E- and B-field using the solution vectors
  // on the true dofs. For the real-valued overload, the electric scalar potential can be
  // specified too for electrostatic simulations. The output mesh and fields are
  // nondimensionalized consistently (B ~ E (L₀ ω₀ E₀⁻¹)).
  void SetEGridFunction(const petsc::PetscParVector &e);
  void SetBGridFunction(const petsc::PetscParVector &b);
  void SetEGridFunction(const mfem::Vector &e);
  void SetBGridFunction(const mfem::Vector &b);
  void SetVGridFunction(const mfem::Vector &v);
  void SetAGridFunction(const mfem::Vector &a);

  // Update cached port voltages and currents for lumped and wave port operators.
  void UpdatePorts(const LumpedPortOperator &lumped_port_op,
                   const WavePortOperator &wave_port_op, double omega = 0.0)
  {
    UpdatePorts(lumped_port_op, omega);
    UpdatePorts(wave_port_op, omega);
  }
  void UpdatePorts(const LumpedPortOperator &lumped_port_op, double omega = 0.0);
  void UpdatePorts(const WavePortOperator &wave_port_op, double omega = 0.0);

  // Postprocess the total electric and magnetic field energies in the electric and magnetic
  // fields.
  double GetEFieldEnergy() const;
  double GetHFieldEnergy() const;

  // Postprocess the energy in lumped capacitor or inductor port boundaries with index in
  // the provided set.
  double GetLumpedInductorEnergy(const LumpedPortOperator &lumped_port_op) const;
  double GetLumpedCapacitorEnergy(const LumpedPortOperator &lumped_port_op) const;

  // Postprocess the S-parameter for recieving lumped or wave port index using the electric
  // field solution.
  std::complex<double> GetSParameter(const LumpedPortOperator &lumped_port_op, int idx,
                                     int source_idx) const;
  std::complex<double> GetSParameter(const WavePortOperator &wave_port_op, int idx,
                                     int source_idx) const;

  // Postprocess the circuit voltage and current across lumped port index using the electric
  // field solution. When has_imaginary is false, the returned voltage has only a nonzero
  // real part.
  std::complex<double> GetPortPower(const LumpedPortOperator &lumped_port_op,
                                    int idx) const;
  std::complex<double> GetPortPower(const WavePortOperator &wave_port_op, int idx) const;
  std::complex<double> GetPortVoltage(const LumpedPortOperator &lumped_port_op,
                                      int idx) const;
  std::complex<double> GetPortVoltage(const WavePortOperator &wave_port_op, int idx) const
  {
    MFEM_ABORT("GetPortVoltage is not yet implemented for wave port boundaries!");
    return 0.0;
  }
  std::complex<double> GetPortCurrent(const LumpedPortOperator &lumped_port_op,
                                      int idx) const;
  std::complex<double> GetPortCurrent(const WavePortOperator &wave_port_op, int idx) const
  {
    MFEM_ABORT("GetPortCurrent is not yet implemented for wave port boundaries!");
    return 0.0;
  }

  // Postprocess the EPR for the electric field solution and lumped port index.
  double GetInductorParticipation(const LumpedPortOperator &lumped_port_op, int idx,
                                  double Em) const;

  // Postprocess the coupling rate for radiative loss to the given I-O port index.
  double GetExternalKappa(const LumpedPortOperator &lumped_port_op, int idx,
                          double Em) const;

  // Postprocess the participation ratio or quality factor for bulk lossy dielectric losses
  // in the electric field mode.
  double GetBulkParticipation(int idx, double Em) const;
  double GetBulkQualityFactor(int idx, double Em) const;

  // Postprocess the partitipation ratio for interface lossy dielectric losses in the
  // electric field mode.
  double GetInterfaceParticipation(int idx, double Em) const;

  // Postprocess the charge or flux for a surface index using the electric field solution
  // or the magnetic flux density field solution.
  double GetSurfaceCharge(int idx) const;
  double GetSurfaceFlux(int idx) const;

  // Write to disk the E- and B-fields extracted from the solution vectors. Note that fields
  // are not redimensionalized, to do so one needs to compute: B <= B * (μ₀ H₀), E <= E *
  // (Z₀ H₀), V <= V * (Z₀ H₀ L₀), etc.
  void WriteFields(int step, double time) const;

  // Probe the E- and B-fields for their vector-values at speceified locations in space.
  // Locations of probes are set up in constructor from configuration file data. If
  // has_imaginary is false, the returned fields have only nonzero real parts. Output
  // vectors are ordered by vector dimension, that is [v1x, v1y, v1z, v2x, v2y, v2z, ...].
  const auto &GetProbes() const { return interp_op.GetProbes(); }
  std::vector<std::complex<double>> ProbeEField() const;
  std::vector<std::complex<double>> ProbeBField() const;
};

}  // namespace palace

#endif  // PALACE_POST_OPERATOR_HPP
