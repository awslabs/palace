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
#include "utils/constants.hpp"

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
// Class to handle domain energy postprocessing. Energies are physical stored energies: for
// real-valued (instantaneous) fields E_elec = 1/2 ∫ D ⋅ E dV, while for complex peak-phasor
// fields the reported value is the time average over one period, E_elec = 1/4 Re{∫ Dᴴ E dV}
// (see electromagnetics::TimeAverageWeight). Ratios of energies (participation ratios) are
// independent of this choice; quality factors formed from time-averaged losses divide by
// the total (electric + magnetic) time-averaged energy, see PostOperator.
//
class DomainPostOperator
{
public:
  // Weight of a field-energy quadratic form; shared with reduced-coordinate evaluation.
  static constexpr double EnergyWeight(bool complex_phasor)
  {
    return 0.5 * electromagnetics::TimeAverageWeight(complex_phasor);
  }

  // Temporary vectors for inner product calculations.
  mutable Vector D, H, D_n;

  // Bilinear forms for computing field energy integrals over domains. The total energies
  // M_elec and M_mag exclude the PML regions, whose energies M_elec_pml and M_mag_pml (and
  // those of the domains M_i in the PML regions) are computed with the material properties
  // of the background medium of the layer: in the PML, the fields are those of the
  // stretched coordinates, and the integrals measure the field in the layer rather than a
  // physical energy.
  std::unique_ptr<Operator> M_elec, M_mag, M_elec_pml, M_mag_pml;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_i;

  // Mode analysis: forms of the out-of-plane E (H1) and in-plane B (ND), which complete
  // those of the in-plane E and out-of-plane B, over the entire domain and each
  // postprocessing domain.
  std::unique_ptr<Operator> M_elec_n, M_mag_t;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_n_i;

  // London domain mass ∫ (μλ²)⁻¹ A·v for magnetostatics, with ½AᵀM_london·A the kinetic
  // energy of the supercurrent (null without London domains).
  std::unique_ptr<Operator> M_london;

  // Unit mass forms of each frequency-dependent material (by material index), over the
  // entire domain and over each postprocessing domain (only the materials in it). With
  // M_elec weighted by ε∞, the electric field energy adds that of the polarization of each
  // material, its unit mass energy times the stored-energy susceptibility wₘ(ω) of its
  // model. A form on every attribute of a material is its unit mass, borrowed when given at
  // construction; the others are owned by M_elec_fd_owned.
  using FrequencyDependentForms = std::vector<std::pair<std::size_t, const Operator *>>;
  FrequencyDependentForms M_elec_fd, M_elec_n_fd;
  std::map<int, FrequencyDependentForms> M_elec_fd_i, M_elec_n_fd_i;
  std::vector<std::unique_ptr<Operator>> M_elec_fd_owned;
  const MaterialOperator *mat_op = nullptr;

  DomainPostOperator() = default;
  // The unit masses of the frequency-dependent materials, by material index (null where not
  // given), are borrowed instead of assembled: the unit masses B_m of SpaceOperator
  // (GetFrequencyDependentPermittivityLocalMassMatrix), which must outlive this object.
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace,
                     const std::vector<const Operator *> &unit_mass = {});
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);
  // Mode analysis: in-plane fields on the ND space, and out-of-plane E and B on the H1 and
  // L2 curl spaces.
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &h1_fespace,
                     const FiniteElementSpace &l2_fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);

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
  // domain (time-averaged for complex-valued fields). The STORED electric energy requires
  // HasStoredElectricFieldEnergy() (or HasStoredDomainElectricFieldEnergy(idx) for a
  // portion of the domain). For mode analysis, E and B are the in-plane E and out-of-plane
  // B, completed by the out-of-plane E (En) and the in-plane B, given as zxBt = z_hat x Bt,
  // which lies in the ND space.
  double GetElectricFieldEnergy(const GridFunction &E, double omega,
                                ElectricEnergy type = ElectricEnergy::STORED,
                                const GridFunction *En = nullptr) const;
  double GetMagneticFieldEnergy(const GridFunction &B,
                                const GridFunction *zxBt = nullptr) const;

  // The STORED (if HasStoredElectricFieldEnergy()) and MATERIAL electric energies in the
  // entire domain, from one set of field products.
  struct ElectricFieldEnergies
  {
    std::optional<double> stored;
    double material;
  };
  ElectricFieldEnergies GetElectricFieldEnergies(const GridFunction &E, double omega,
                                                 const GridFunction *En = nullptr) const;

  // The electric field energy in the entire domain with the high-frequency permittivity ε∞
  // of each frequency-dependent material, and the permittivity of the constant materials,
  // without the energy of the polarization. Available for every material, and positive for
  // a nonzero field also at a zero of a permittivity.
  double GetHighFrequencyElectricFieldEnergy(const GridFunction &E,
                                             const GridFunction *En = nullptr) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain (time-averaged for complex-valued fields).
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E, double omega,
                                      ElectricEnergy type = ElectricEnergy::STORED,
                                      const GridFunction *En = nullptr) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B,
                                      const GridFunction *zxBt = nullptr) const;

  // The magnetic field energies from a real L-vector of the space of construction, for
  // spaces without MFEM grid functions (broken spaces). Real-valued fields, so these are
  // instantaneous energies.
  double GetMagneticFieldEnergy(const Vector &a_l, MPI_Comm comm) const;
  double GetDomainMagneticFieldEnergy(int idx, const Vector &a_l, MPI_Comm comm) const;

  // Get volume integrals for the electric or magnetic field energy in the PML regions (zero
  // without PML regions).
  bool HasPML() const { return M_elec_pml != nullptr; }
  double GetPMLElectricFieldEnergy(const GridFunction &E) const;
  double GetPMLMagneticFieldEnergy(const GridFunction &B) const;

  // The kinetic energy of the London domains from a real L-vector of the space of
  // construction. A real-valued field, so an instantaneous energy.
  double GetLondonKineticEnergy(const Vector &a_l, MPI_Comm comm) const;

  // The susceptibility weighting the unit mass form of material m: wₘ(ω) for STORED (which
  // requires MaterialOperator::HasStoredElectricEnergy(m)) and Re χₘ(iω) for MATERIAL.
  // Errors at a singularity of the material (an exact lossless resonance), as do the
  // energies with it: they are not finite there.
  double GetFrequencyDependentEnergyWeight(std::size_t material_idx, double omega,
                                           ElectricEnergy type) const;

private:
  // The forms of the frequency-dependent materials on the ND or H1 space, over the entire
  // domain (M_fd) and each postprocessing domain (M_fd_i).
  void SetUpFrequencyDependentForms(const config::DomainPostData &postpro,
                                    const MaterialOperator &mat_op,
                                    const FiniteElementSpace &fespace,
                                    const std::vector<const Operator *> &unit_mass,
                                    FrequencyDependentForms &M_fd,
                                    std::map<int, FrequencyDependentForms> &M_fd_i);
  bool HasStoredEnergy(const FrequencyDependentForms &M_fd) const;
  // The electric energy of E with the forms (M, M_fd) and, for mode analysis, of En with
  // (M_n, M_n_fd).
  double ElectricFieldEnergy(const Operator &M, const FrequencyDependentForms &M_fd,
                             const GridFunction &E, const Operator *M_n,
                             const FrequencyDependentForms *M_n_fd, const GridFunction *En,
                             double omega, ElectricEnergy type) const;

  // Evaluate the energy quadratic form 1/2 uᴴ M u for a real or complex field u, with the
  // time-average weight applied for complex peak phasors.
  static double GetFieldEnergy(const Operator &M, const GridFunction &u, Vector &Mu);
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
