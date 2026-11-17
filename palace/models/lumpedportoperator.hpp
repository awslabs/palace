// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP
#define PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/lumpedelement.hpp"
#include "fem/output_functionals.hpp"
#include "linalg/operator.hpp"

namespace palace
{

class FiniteElementSpace;
class GridFunction;
class IoData;
class MaterialOperator;
class Units;
class MaterialPropertyCoefficient;
class SumVectorCoefficient;

namespace config
{

struct LumpedPortData;

}  // namespace config

//
// Helper class for lumped ports in a model.
//
class LumpedPortData
{
  friend class LumpedPortOperator;

public:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // To accommodate multielement lumped ports, a port may be made up of elements with
  // different attributes and directions which add in parallel.
  std::vector<std::unique_ptr<LumpedElementData>> elems;

  // Lumped port properties.
  double R, L, C;
  int excitation;
  bool active;
  bool include_in_synthesis;

private:
  // Complex weight a = A exp(i φ) on the unit-power incident wave of this port, from the
  // configured "ExcitationAmplitude" and "ExcitationPhase". It is dimensionless and
  // frequency independent. It is 1 for a port without excitation, so it is private: use
  // GetExcitationWeight, which is 0 for such a port.
  std::complex<double> excitation_weight = 1.0;

protected:
  // Linear forms for postprocessing integrated quantities on the port.
  mutable std::unique_ptr<mfem::LinearForm> s, v;

  // libCEED surface functionals for port power, voltage, and S-parameter computation,
  // replacing host boundary LinearForm paths when supported.
  mutable std::unique_ptr<SurfaceFunctional> power_func, v_func;

  // Boundary mass operator of the resistive sheet, for the dissipated port power.
  mutable std::unique_ptr<Operator> resistor_op;

  void InitializeLinearForms(mfem::ParFiniteElementSpace &nd_fespace) const;
  SurfaceModeCoefficient GetModeCoefficient(const LumpedElementData &elem,
                                            double coeff) const;

public:
  LumpedPortData(const config::LumpedPortData &data, const MaterialOperator &mat_op,
                 const mfem::ParMesh &mesh);

  double GetToSquare(const LumpedElementData &elem) const
  {
    return elem.GetGeometryWidth() / elem.GetGeometryLength() * elems.size();
  }

  // Normalization of tangential electric field of port corresponding ∫|E|²ds = |Z₀| *
  // ∑(Wₑ/Lₑ). In this function we set the impedance magnitude|Z₀| = 1 in internal units (=
  // Z_freespace Ohm). The actual lumped impedance is frequency dependant for ports with L,C
  // so it is more convenient to add this factor later as needed.
  double GetExcitationFieldEtNormSqWithUnityZR() const
  {
    double norm_et = 0.0;
    for (const auto &el : elems)
    {
      norm_et += el->GetGeometryWidth() / el->GetGeometryLength();
    }
    return norm_et;
  }

  // Normalization of tangential electric field of port corresponding ∫|H|²ds = 1 / (|Z₀|
  // nₑₗₑₘₛ²) * ∑(Lₑ/Wₑ). Same details as in GetExcitationFieldEtNormSqWithUnityZR above.
  double GetExcitationFieldHtNormSqWithUnityZR() const
  {
    double norm_ht = 0.0;
    for (const auto &el : elems)
    {
      norm_ht += el->GetGeometryLength() / el->GetGeometryWidth();
    }
    norm_ht /= elems.size() * elems.size();
    return norm_ht;
  }

  constexpr bool HasExcitation() const { return excitation != 0; }

  // Real reference resistance used to normalize the incident-field amplitude of an excited
  // port. For a resistive port this is R, so the incident/reflected-wave normalization and
  // the port's own termination impedance coincide and behaviour is unchanged. For a purely
  // reactive excited port (R == 0) there is no real port resistance, so we fall back to the
  // unit reference impedance |Z_R| = 1 in internal units (= Z_freespace) purely to define a
  // finite drive amplitude; the reactance acts physically through the system-matrix
  // termination, and that port's own S column is referenced to this unit impedance.
  double GetExcitationRefResistance() const { return (std::abs(R) > 0.0) ? R : 1.0; }

  enum class Branch
  {
    TOTAL,
    R,
    L,
    C
  };
  std::complex<double> GetCharacteristicImpedance(double omega = 0.0,
                                                  Branch branch = Branch::TOTAL) const;

  // UNIT-POWER incident quantities of the excited port, which do not include the excitation
  // weight: the incident power is 1 and the voltage is that of this unit-power wave. Use
  // these for normalization (for example the S-parameters). For the incident voltage and
  // current actually applied and reported, which include the weight, use
  // GetIncidentVoltage and GetIncidentCurrent below.
  double GetExcitationPower() const;
  double GetExcitationVoltage() const;

  // Complex weight a = A exp(i φ) applied to the incident wave of an excited port when the
  // excitation is assembled and when incident quantities are reported; 0 for a port
  // without excitation. The S-parameter normalization (GetSParameter) does not include it.
  std::complex<double> GetExcitationWeight() const
  {
    return HasExcitation() ? excitation_weight : std::complex<double>(0.0, 0.0);
  }

  // Whether the excitation weight has a nonzero imaginary part, such that the excitation
  // vector of the port is complex.
  bool HasComplexExcitationWeight() const
  {
    return HasExcitation() && excitation_weight.imag() != 0.0;
  }

  // WEIGHTED incident voltage a V_inc,0 and current a I_inc,0 = a P_inc,0 / V_inc,0 of an
  // excited port, in internal units, with V_inc,0 = GetExcitationVoltage() the unit-power
  // value and a = GetExcitationWeight(). These are the incident values applied to the port
  // and written to the output files.
  std::complex<double> GetIncidentVoltage() const
  {
    return GetExcitationWeight() * GetExcitationVoltage();
  }
  std::complex<double> GetIncidentCurrent() const;

  std::complex<double> GetPower(GridFunction &E, GridFunction &B) const;
  std::complex<double> GetPowerLegacy(GridFunction &E, GridFunction &B) const;
  std::complex<double> GetSParameter(GridFunction &E) const;
  std::complex<double> GetVoltage(GridFunction &E) const;

  // Add the sheet admittance 1 / R_s of each port element, scaled by coeff, to the boundary
  // coefficient. This term is zero for a port without resistance. It does not depend on
  // whether the port is active.
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb) const;

  // Time-averaged power dissipated in the resistive sheet of the port,
  //                           P = 1/2 ∫ |E_t|² / R_s dS,
  // for the field E on nd_fespace. This is the port contribution to the boundary quadratic
  // of the damping matrix. For R > 0, the modal voltage gives |V|² / (2R) ≤ P, with
  // equality only when E_t on each element matches the voltage mode of the port and all
  // elements have the same voltage: field components orthogonal to that mode do not change
  // V, but they dissipate power. The boundary operator is assembled on the first call.
  double GetResistorPower(const FiniteElementSpace &nd_fespace,
                          const GridFunction &E) const;
};

//
// A class handling lumped port boundaries and their postprocessing.
//
class LumpedPortOperator
{
private:
  // Mapping from port index to data structure containing port information and methods to
  // calculate circuit properties like voltage and current on lumped or multielement lumped
  // ports.
  std::map<int, LumpedPortData> ports;

  struct PortAttributeBins
  {
    mfem::Array<int> attr_to_port;
    std::vector<int> port_indices;
    bool unavailable = false;
  };

  // Batched libCEED surface functionals for lumped-port power and voltage. Individual
  // LumpedPortData evaluation remains available for port sets that cannot be batched.
  mutable std::unique_ptr<SurfaceFunctional> batched_power_func;
  mutable PortAttributeBins batched_power_bins;
  mutable std::unique_ptr<SurfaceFunctional> batched_voltage_func;
  mutable PortAttributeBins batched_voltage_bins;

  void SetUpBoundaryProperties(const std::map<int, config::LumpedPortData> &lumpedport,
                               const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

public:
  LumpedPortOperator(const std::map<int, config::LumpedPortData> &lumpedport,
                     const Units &units, const MaterialOperator &mat_op,
                     const mfem::ParMesh &mesh);
  LumpedPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const mfem::ParMesh &mesh);

  // Access data structures for the lumped port with the given index.
  const LumpedPortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Compute port powers, batching disjoint libCEED surface-functional evaluations when
  // possible while preserving per-port scalar outputs.
  std::map<int, std::complex<double>> GetPowers(GridFunction &E, GridFunction &B) const;

  // Compute port voltages, batching disjoint libCEED mode-overlap evaluations when
  // possible while preserving per-port scalar outputs.
  std::map<int, std::complex<double>> GetVoltages(GridFunction &E) const;

  // Returns array of lumped port attributes.
  mfem::Array<int> GetAttrList() const;
  mfem::Array<int> GetRsAttrList() const;
  mfem::Array<int> GetLsAttrList() const;
  mfem::Array<int> GetCsAttrList() const;

  // Add contributions to system matrices from lumped elements with nonzero inductance,
  // resistance, and/or capacitance.
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries, -U_inc/(iω) for the real version (versus the full -U_inc for
  // the complex one). The incident field of each port includes its excitation weight, whose
  // real part is added to fbr, usable in the time domain, and whose imaginary part is added
  // to fbi, for the frequency domain only. Both are frequency independent.
  void AddExcitationBdrCoefficients(int excitation_idx, SumVectorCoefficient &fbr,
                                    SumVectorCoefficient &fbi);

  // Same drive term for a single port selected by port index, with unit weight, independent
  // of any "Excitation" grouping and excitation weight in the configuration. Used for the
  // quasistatic anchor solves in circuit synthesis.
  void AddPortExcitationBdrCoefficient(int port_idx, SumVectorCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP
