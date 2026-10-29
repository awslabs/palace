// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class GridFunction;
class FiniteElementSpace;
class IoData;
class MaterialOperator;

namespace config
{

struct DomainPostData;

}  // namespace config

//
// Class to handle domain energy postprocessing. We use a leading factor of 1/2 instead of
// 1/4 even though the eigenmodes are peak phasors and not RMS normalized because the same
// peak phasors are used to compute the voltages/currents which are 2x the time-averaged
// values. This correctly yields an EPR of 1 in cases where expected.
//
class DomainPostOperator
{
public:
  // Temporary vectors for inner product calculations.
  mutable Vector D, H;

  // Bilinear forms for computing field energy integrals over domains.
  std::unique_ptr<Operator> M_elec, M_mag;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_i;

  // Unit mass forms of each frequency-dependent material (by material index), over the
  // entire domain and over each postprocessing domain (only the materials in it). With
  // M_elec weighted by ε∞, the electric field energy adds that of the polarization of each
  // material, its unit mass energy times the stored-energy susceptibility wₘ(ω) of its
  // model.
  using FrequencyDependentForms =
      std::vector<std::pair<std::size_t, std::unique_ptr<Operator>>>;
  FrequencyDependentForms M_elec_fd;
  std::map<int, FrequencyDependentForms> M_elec_fd_i;
  const MaterialOperator *mat_op = nullptr;

  DomainPostOperator() = default;
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace);
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);

  // Mode analysis: ND space for electric field energy only (no B-field).
  // Uses VectorFEMassIntegrator with ε (not CurlCurlIntegrator).
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace, bool electric_energy_only);

  // The electric energy of a frequency-dependent material at the real angular frequency ω
  // (unused without one): STORED is that of its model, the field and the polarization;
  // MATERIAL is the field energy at its permittivity Re ε(ω), which for a lossless mode
  // equals the magnetic energy. They coincide for a constant material.
  enum class ElectricEnergy
  {
    STORED,
    MATERIAL
  };

  // Whether the STORED electric energy is available in the entire domain, or in the
  // postprocessing domain idx: whether every frequency-dependent material in it has a
  // stored electric energy (MaterialOperator::HasStoredElectricEnergy). This is fixed for
  // the run and the same on every process. The MATERIAL energy, the magnetic energy, and
  // the energies of static problems are always available.
  bool HasStoredElectricFieldEnergy() const;
  bool HasStoredDomainElectricFieldEnergy(int idx) const;

  // Get volume integrals computing the electric or magnetic field energy in the entire
  // domain. The STORED electric energy requires HasStoredElectricFieldEnergy() (or
  // HasStoredDomainElectricFieldEnergy(idx) for a portion of the domain).
  double GetElectricFieldEnergy(const GridFunction &E, double omega,
                                ElectricEnergy type = ElectricEnergy::STORED) const;
  double GetMagneticFieldEnergy(const GridFunction &B) const;

  // The STORED (if HasStoredElectricFieldEnergy()) and MATERIAL electric energies in the
  // entire domain, from one set of field products.
  struct ElectricFieldEnergies
  {
    std::optional<double> stored;
    double material;
  };
  ElectricFieldEnergies GetElectricFieldEnergies(const GridFunction &E, double omega) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain.
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E, double omega,
                                      ElectricEnergy type = ElectricEnergy::STORED) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B) const;

  // The susceptibility weighting the unit mass form of material m: wₘ(ω) for STORED (which
  // requires MaterialOperator::HasStoredElectricEnergy(m)) and Re χₘ(iω) for MATERIAL.
  // Errors at a singularity of the material (an exact lossless resonance), as do the
  // energies with it: they are not finite there.
  double GetFrequencyDependentEnergyWeight(std::size_t material_idx, double omega,
                                           ElectricEnergy type) const;

private:
  void SetUpFrequencyDependentForms(const config::DomainPostData &postpro,
                                    const MaterialOperator &mat_op,
                                    const FiniteElementSpace &nd_fespace);
  bool HasStoredEnergy(const FrequencyDependentForms &M_fd) const;
  double ElectricFieldEnergy(const Operator &M, const FrequencyDependentForms &M_fd,
                             const GridFunction &E, double omega,
                             ElectricEnergy type) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
