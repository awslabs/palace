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

  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);

  // Get volume integrals computing the electric or magnetic field energy in the entire
  // domain.
  double GetElectricFieldEnergy(const GridFunction &E) const;
  double GetMagneticFieldEnergy(const GridFunction &B) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain.
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
