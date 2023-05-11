// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <map>
#include <optional>
#include <utility>
#include <mfem.hpp>

namespace palace
{

class IoData;
class MaterialOperator;

//
// A class handling domain postprocessing (matrix-free).
//
class DomainPostOperatorMF
{
private:
  // Reference to material property operator (not owned).
  const MaterialOperator &mat_op;

  // Unit function used for computing volume integrals.
  mfem::ParGridFunction ones;

  // Mapping from domain index to marker and loss tangent for postprocessing bulk dielectic
  // loss.
  mutable std::map<int, mfem::Array<int>> attr_markers;

  double GetVolumeIntegral(mfem::Coefficient &f) const;
  double GetVolumeIntegral(mfem::Coefficient &f, mfem::Array<int> &attr_marker) const;

public:
  DomainPostOperatorMF(const IoData &iodata, const MaterialOperator &mat,
                       mfem::ParFiniteElementSpace &h1_fespace);

  // Access underlying bulk loss postprocessing data structures (for keys).
  const auto &GetEps() const { return attr_markers; }
  auto SizeEps() const { return attr_markers.size(); }

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

//
// A class handling domain postprocessing.
//
class DomainPostOperator
{
private:
  // Bilinear forms for computing field energy integrals over domains.
  std::optional<mfem::ParBilinearForm> m0ND, m0RT;
  std::map<int, std::pair<mfem::ParBilinearForm, mfem::ParBilinearForm>> m0NDi;

public:
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     mfem::ParFiniteElementSpace *nd_fespace,
                     mfem::ParFiniteElementSpace *rt_fespace);

  // Access underlying bulk loss postprocessing data structures (for keys).
  const auto &GetEps() const { return m0NDi; }
  auto SizeEps() const { return m0NDi.size(); }

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
