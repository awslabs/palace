// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <map>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;

//
// A class handling material attributes.
//
class MaterialOperator
{
private:
  // Material properties for domain attributes: relative permeability, relative
  // permittivity, and others (like electrical conductivity and London penetration depth
  // for superconductors. The i-1-th entry of each Vector is the property for mesh domain
  // attribute i. Marker arrays contain a 1 for each domain attribute labeled, and 0 else.
  std::vector<mfem::DenseMatrix> mat_muinv, mat_epsilon, mat_epsilon_imag, mat_epsilon_abs,
      mat_invz0, mat_c0, mat_sigma, mat_invLondon;
  std::vector<double> mat_c0_min, mat_c0_max;
  mfem::Array<int> losstan_marker, conductivity_marker, london_marker;
  void SetUpMaterialProperties(const IoData &iodata, mfem::ParMesh &mesh);

  // Shared face mapping for boundary coefficients.
  std::map<int, int> local_to_shared;

public:
  MaterialOperator(const IoData &iodata, mfem::ParMesh &mesh);

  int SpaceDimension() const { return mat_muinv.front().Height(); }

  const auto &GetLocalToSharedFaceMap() const { return local_to_shared; }

  const auto &GetInvPermeability(int attr) const { return mat_muinv[attr - 1]; }
  const auto &GetPermittivityReal(int attr) const { return mat_epsilon[attr - 1]; }
  const auto &GetPermittivityImag(int attr) const { return mat_epsilon_imag[attr - 1]; }
  const auto &GetPermittivityAbs(int attr) const { return mat_epsilon_abs[attr - 1]; }
  const auto &GetInvImpedance(int attr) const { return mat_invz0[attr - 1]; }
  const auto &GetLightSpeed(int attr) const { return mat_c0[attr - 1]; }
  const auto &GetLightSpeedMin(int attr) const { return mat_c0_min[attr - 1]; }
  const auto &GetLightSpeedMax(int attr) const { return mat_c0_max[attr - 1]; }
  const auto &GetConductivity(int attr) const { return mat_sigma[attr - 1]; }
  const auto &GetInvLondonDepth(int attr) const { return mat_invLondon[attr - 1]; }

  bool HasLossTangent() const { return (losstan_marker.Max() > 0); }
  bool HasConductivity() const { return (conductivity_marker.Max() > 0); }
  bool HasLondonDepth() const { return (london_marker.Max() > 0); }

  const auto &GetLossTangentMarker() const { return losstan_marker; }
  const auto &GetConductivityMarker() const { return conductivity_marker; }
  const auto &GetLondonDepthMarker() const { return london_marker; }
};

}  // namespace palace

#endif  // PALACE_MODELS_MATERIAL_OPERATOR_HPP
