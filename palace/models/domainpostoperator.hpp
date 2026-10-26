// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <utility>
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

  // Bilinear forms for computing field energy integrals over domains. The total energies
  // M_elec and M_mag exclude the PML regions, whose energies M_elec_pml and M_mag_pml (and
  // those of the domains M_i in the PML regions) are computed with the material properties
  // of the background medium of the layer: in the PML, the fields are those of the
  // stretched coordinates, and the integrals measure the field in the layer rather than a
  // physical energy.
  std::unique_ptr<Operator> M_elec, M_mag, M_elec_pml, M_mag_pml;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_i;

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

  // Get volume integrals computing the electric or magnetic field energy in the entire
  // domain.
  double GetElectricFieldEnergy(const GridFunction &E) const;
  double GetMagneticFieldEnergy(const GridFunction &B) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain.
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B) const;

  // The magnetic field energies from a real L-vector of the space of construction, for
  // spaces without MFEM grid functions (broken spaces).
  double GetMagneticFieldEnergy(const Vector &a_l, MPI_Comm comm) const;
  double GetDomainMagneticFieldEnergy(int idx, const Vector &a_l, MPI_Comm comm) const;

  // Get volume integrals for the electric or magnetic field energy in the PML regions (zero
  // without PML regions).
  bool HasPML() const { return M_elec_pml != nullptr; }
  double GetPMLElectricFieldEnergy(const GridFunction &E) const;
  double GetPMLMagneticFieldEnergy(const GridFunction &B) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
