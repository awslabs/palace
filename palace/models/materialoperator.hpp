// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <unordered_map>
#include <mfem.hpp>

namespace palace
{

class IoData;
class Mesh;

//
// A class handling material attributes.
//
class MaterialOperator
{
private:
  // Useful references to objects from the underlying mesh (not owned).
  const std::unordered_map<int, int> &loc_attr, &loc_bdr_attr;
  const std::unordered_map<int, int> &local_to_shared;

  // Mapping from the local attribute to material index. For boundary elements, this is the
  // material of the neighboring element (for internal boundaries, use the element which
  // corresponds to the vacuum domain, or at least the one with the higher speed of light.
  mfem::Array<int> attr_mat, bdr_attr_mat;

  // Material properties: relative permeability, relative permittivity, and others (like
  // electrical conductivity and London penetration depth for superconductors.
  mfem::DenseTensor mat_muinv, mat_epsilon, mat_epsilon_imag, mat_epsilon_abs, mat_invz0,
      mat_c0, mat_sigma, mat_invLondon;
  mfem::Array<double> mat_c0_min, mat_c0_max;

  // Domain attributes with nonzero loss tangent, electrical conductivity, London
  // penetration depth.
  mfem::Array<int> losstan_attr, conductivity_attr, london_attr;

  void SetUpMaterialProperties(const IoData &iodata, const Mesh &mesh);

  const auto AttrToMat(int attr, bool bdr)
  {
    if (bdr)
    {
      MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end(),
                  "Missing local boundary attribute for attribute " << attr << "!");
      return bdr_attr_mat[loc_bdr_attr[attr] - 1];
    }
    else
    {
      MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                  "Missing local domain attribute for attribute " << attr << "!");
      return attr_mat[loc_attr[attr] - 1];
    }
  }

  const auto Wrap(const mfem::DenseTensor &data, int attr, bool bdr) const
  {
    const int k = AttrToMat(attr, bdr);
    return mfem::DenseMatrix(const_cast<double *>(data.GetData(k)), data.SizeI(),
                             data.SizeJ());
  }

public:
  MaterialOperator(const IoData &iodata, const Mesh &mesh);

  int SpaceDimension() const { return mat_muinv.SizeI(); }

  const auto GetInvPermeability(int attr, bool bdr = false) const
  {
    return Wrap(mat_muinv, attr, bdr);
  }
  const auto GetPermittivityReal(int attr, bool bdr = false) const
  {
    return Wrap(mat_epsilon, attr, bdr);
  }
  const auto GetPermittivityImag(int attr, bool bdr = false) const
  {
    return Wrap(mat_epsilon_imag, attr, bdr);
  }
  const auto GetPermittivityAbs(int attr, bool bdr = false) const
  {
    return Wrap(mat_epsilon_abs, attr, bdr);
  }
  const auto GetInvImpedance(int attr, bool bdr = false) const
  {
    return Wrap(mat_invz0, attr, bdr);
  }
  const auto GetLightSpeed(int attr, bool bdr = false) const
  {
    return Wrap(mat_c0, attr, bdr);
  }
  const auto GetConductivity(int attr, bool bdr = false) const
  {
    return Wrap(mat_sigma, attr, bdr);
  }
  const auto GetInvLondonDepth(int attr, bool bdr = false) const
  {
    return Wrap(mat_invLondon, attr, bdr);
  }

  auto GetLightSpeedMin(int attr, bool bdr = false) const
  {
    return mat_c0_min[AttrToMat(attr, bdr)];
  }
  auto GetLightSpeedMax(int attr, bool bdr = false) const
  {
    return mat_c0_max[AttrToMat(attr, bdr)];
  }

  const auto &GetInvPermeability() const { return mat_muinv; }
  const auto &GetPermittivityReal() const { return mat_epsilon; }
  const auto &GetPermittivityImag() const { return mat_epsilon_imag; }
  const auto &GetPermittivityAbs() const { return mat_epsilon_abs; }
  const auto &GetInvImpedance() const { return mat_invz0; }
  const auto &GetLightSpeed() const { return mat_c0; }
  const auto &GetConductivity() const { return mat_sigma; }
  const auto &GetInvLondonDepth() const { return mat_invLondon; }

  bool HasLossTangent() const { return (losstan_attr.Size() > 0); }
  bool HasConductivity() const { return (conductivity_attr.Size() > 0); }
  bool HasLondonDepth() const { return (london_attr.Size() > 0); }

  const auto &GetAttributeToMaterial() const { return attr_mat; }
  const auto &GetBdrAttributeToMaterial() const { return bdr_attr_mat; }

  const auto &GetAttributeGlobalToLocal() const { return loc_attr; }
  const auto &GetBdrAttributeGlobalToLocal() const { return loc_bdr_attr; }

  auto GetAttributeGlobalToLocal(const mfem::Array<int> &attr_list) const
  {
    mfem::Array<int> loc_attr_list(attr_list.Size());
    std::transform(attr_list.begin(), attr_list.end(), loc_attr_list.begin(),
                   [&loc_attr](int attr) { return loc_attr[attr]; }) return loc_attr_list;
  }
  auto GetBdrAttributeGlobalToLocal(const mfem::Array<int> &attr_list) const
  {
    mfem::Array<int> loc_attr_list(attr_list.Size());
    std::transform(attr_list.begin(), attr_list.end(), loc_attr_list.begin(),
                   [&loc_bdr_attr](int attr)
                   { return loc_bdr_attr[attr]; }) return loc_attr_list;
  }

  const auto &GetLocalToSharedFaceMap() const { return local_to_shared; }
};

//
// Material property represented as a piecewise constant coefficient over mesh elements. Can
// be scalar-valued or matrix-valued.
//
class MaterialPropertyCoefficient
{
private:
  // Map attribute to material index (coeff = mat_coeff[attr_mat[attr - 1]], for 1-based
  // attributes).
  mfem::Array<int> attr_mat;

  // Material properry coefficients, ordered by material index.
  mfem::DenseTensor mat_coeff;

public:
  MaterialPropertyCoefficient() {}
  MaterialPropertyCoefficient(const mfem::Array<int> &attr_mat,
                              const mfem::DenseTensor &mat_coeff, double a = 1.0);

  bool empty() const { return mat_coeff.TotalSize() == 0; }

  const auto &GetAttributeToMaterial() const { return attr_mat; }
  const auto &GetMaterialProperties() const { return mat_coeff; }

  void AddCoefficient(const mfem::Array<int> &attr_mat_,
                      const mfem::DenseTensor &mat_coeff_, double a = 1.0);

  template <typename T>
  void AddMaterialProperty(const mfem::Array<int> &attr_list, const T &coeff,
                           double a = 1.0);

  void RestrictCoefficient(const mfem::Array<int> &attr_list);

  void NormalProjectedCoefficient(const mfem::Vector &normal);
};

}  // namespace palace

#endif  // PALACE_MODELS_MATERIAL_OPERATOR_HPP
