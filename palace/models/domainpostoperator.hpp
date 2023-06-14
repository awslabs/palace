// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <map>
#include <optional>
#include <utility>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;

//
// A class handling domain postprocessing.
//
class DomainPostOperator
{
private:
  // Bilinear forms for computing field energy integrals over domains.
  std::optional<mfem::BilinearForm> M_ND, M_RT;
  std::map<int, std::pair<mfem::BilinearForm, mfem::BilinearForm>> M_NDi;

  // Temporary vectors for inner product calculations.
  mutable Vector D, H;

public:
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     mfem::ParFiniteElementSpace *nd_fespace,
                     mfem::ParFiniteElementSpace *rt_fespace);

  // Access underlying bulk loss postprocessing data structures (for keys).
  const auto &GetEps() const { return M_NDi; }
  auto SizeEps() const { return M_NDi.size(); }

  // Get volume integrals computing bulk electric or magnetic field energy.
  double GetElectricFieldEnergy(const mfem::ParComplexGridFunction &E) const;
  double GetElectricFieldEnergy(const mfem::ParGridFunction &E) const;
  double GetMagneticFieldEnergy(const mfem::ParComplexGridFunction &B) const;
  double GetMagneticFieldEnergy(const mfem::ParGridFunction &B) const;
  double GetDomainElectricFieldEnergy(int idx, const mfem::ParComplexGridFunction &E) const;
  double GetDomainElectricFieldEnergy(int idx, const mfem::ParGridFunction &E) const;
  double GetDomainElectricFieldEnergyLoss(int idx,
                                          const mfem::ParComplexGridFunction &E) const;
  double GetDomainElectricFieldEnergyLoss(int idx, const mfem::ParGridFunction &E) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
