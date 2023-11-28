// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <map>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

// XX TODO WIP MATERIAL PROPERTY COEFFICIENTS, ELEMENT ATTRIBUTE VECTORS

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
  // for superconductors. Marker arrays contain a 1 for each domain attribute labeled, and
  // 0 else.
  std::vector<int> mat_idx;
  std::vector<mfem::DenseMatrix> mat_muinv, mat_epsilon, mat_epsilon_imag, mat_epsilon_abs,
      mat_invz0, mat_c0, mat_sigma, mat_invLondon;
  std::vector<double> mat_c0_min, mat_c0_max;
  mfem::Array<int> losstan_marker, conductivity_marker, london_marker;
  void SetUpMaterialProperties(const IoData &iodata, const mfem::ParMesh &mesh);

  // Shared face mapping for boundary coefficients.
  std::map<int, int> local_to_shared;

  // Data structure with information used to assemble libCEED operators.
  ceed::CeedObjectMap<ceed::CeedGeomFactorData> geom_data;

public:
  MaterialOperator(const IoData &iodata, mfem::ParMesh &mesh);

  int SpaceDimension() const { return mat_muinv.front().Height(); }

  const auto &GetInvPermeability(int attr) const { return mat_muinv[mat_idx[attr - 1]]; }
  const auto &GetPermittivityReal(int attr) const { return mat_epsilon[mat_idx[attr - 1]]; }
  const auto &GetPermittivityImag(int attr) const
  {
    return mat_epsilon_imag[mat_idx[attr - 1]];
  }
  const auto &GetPermittivityAbs(int attr) const
  {
    return mat_epsilon_abs[mat_idx[attr - 1]];
  }
  const auto &GetInvImpedance(int attr) const { return mat_invz0[mat_idx[attr - 1]]; }
  const auto &GetLightSpeed(int attr) const { return mat_c0[mat_idx[attr - 1]]; }
  const auto &GetLightSpeedMin(int attr) const { return mat_c0_min[mat_idx[attr - 1]]; }
  const auto &GetLightSpeedMax(int attr) const { return mat_c0_max[mat_idx[attr - 1]]; }
  const auto &GetConductivity(int attr) const { return mat_sigma[mat_idx[attr - 1]]; }
  const auto &GetInvLondonDepth(int attr) const { return mat_invLondon[mat_idx[attr - 1]]; }

  bool HasLossTangent() const { return (losstan_marker.Max() > 0); }
  bool HasConductivity() const { return (conductivity_marker.Max() > 0); }
  bool HasLondonDepth() const { return (london_marker.Max() > 0); }

  const auto &GetLossTangentMarker() const { return losstan_marker; }
  const auto &GetConductivityMarker() const { return conductivity_marker; }
  const auto &GetLondonDepthMarker() const { return london_marker; }

  const auto &GetLocalToSharedFaceMap() const { return local_to_shared; }

  const auto &GetCeedGeomFactorData() const { return geom_data; }
  const auto &GetCeedGeomFactorData(Ceed ceed, mfem::Geometry::Type geom) const
  {
    const auto it = geom_data.find(std::make_pair(ceed, geom));
    MFEM_ASSERT(it != geom_data.end(), "Unable to geometry factor data for geometry "
                                           << mfem::Geometry::Name[geom] << "!");
    return it->second;
  }
};

}  // namespace palace

#endif  // PALACE_MODELS_MATERIAL_OPERATOR_HPP
