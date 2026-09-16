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
  // Temporary vectors for inner product calculations.
  mutable Vector D, H;

  // Bilinear forms for computing field energy integrals over domains.
  std::unique_ptr<Operator> M_elec, M_mag;
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
  // domain (time-averaged for complex-valued fields).
  double GetElectricFieldEnergy(const GridFunction &E) const;
  double GetMagneticFieldEnergy(const GridFunction &B) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain (time-averaged for complex-valued fields).
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B) const;

private:
  // Evaluate the energy quadratic form 1/2 uᴴ M u for a real or complex field u, with the
  // time-average weight applied for complex peak phasors.
  static double GetFieldEnergy(const Operator &M, const GridFunction &u, Vector &Mu);
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
