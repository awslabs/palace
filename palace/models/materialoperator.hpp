// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <unordered_map>
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
  // Reference to underlying mesh object (not owned).
  const mfem::ParMesh &mesh;

  // Mapping from the local attribute to material index.
  mfem::Array<int> attr_mat;

  // Material properties: relative permeability, relative permittivity, and others (like
  // electrical conductivity and London penetration depth for superconductors.
  mfem::DenseTensor mat_muinv, mat_epsilon, mat_epsilon_imag, mat_epsilon_abs, mat_invz0,
      mat_c0, mat_sigma, mat_invLondon;
  mfem::Array<double> mat_c0_min, mat_c0_max;

  // Domain attributes with nonzero loss tangent, electrical conductivity, London
  // penetration depth.
  mfem::Array<int> losstan_attr, conductivity_attr, london_attr;

  // Attribute mapping for (global, 1-based) domain and boundary attributes to those on this
  // process (still 1-based). For boundaries, the inner map is a mapping from neighboring
  // domain attribute to the resulting local boundary attribute (to discern boundary
  // elements with global boundary attribute which borders more than one domain). Interior
  // boundaries use as neighbor the element with the smaller domain attribute in order to
  // be consistent when the interior boundary element normals are not aligned.
  std::unordered_map<int, int> loc_attr;
  std::unordered_map<int, std::unordered_map<int, int>> loc_bdr_attr;

  void SetUpMaterialProperties(const IoData &iodata, const mfem::ParMesh &mesh);

  const auto AttrToMat(int attr) const
  {
    MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                "Missing local domain attribute for attribute " << attr << "!");
    return attr_mat[loc_attr.at(attr) - 1];
  }

  const auto Wrap(const mfem::DenseTensor &data, int attr) const
  {
    const int k = AttrToMat(attr);
    return mfem::DenseMatrix(const_cast<double *>(data.GetData(k)), data.SizeI(),
                             data.SizeJ());
  }

public:
  MaterialOperator(const IoData &iodata, mfem::ParMesh &mesh);

  int SpaceDimension() const { return mat_muinv.SizeI(); }

  const auto GetInvPermeability(int attr) const { return Wrap(mat_muinv, attr); }
  const auto GetPermittivityReal(int attr) const { return Wrap(mat_epsilon, attr); }
  const auto GetPermittivityImag(int attr) const { return Wrap(mat_epsilon_imag, attr); }
  const auto GetPermittivityAbs(int attr) const { return Wrap(mat_epsilon_abs, attr); }
  const auto GetInvImpedance(int attr) const { return Wrap(mat_invz0, attr); }
  const auto GetLightSpeed(int attr) const { return Wrap(mat_c0, attr); }
  const auto GetConductivity(int attr) const { return Wrap(mat_sigma, attr); }
  const auto GetInvLondonDepth(int attr) const { return Wrap(mat_invLondon, attr); }

  auto GetLightSpeedMin(int attr) const { return mat_c0_min[AttrToMat(attr)]; }
  auto GetLightSpeedMax(int attr) const { return mat_c0_max[AttrToMat(attr)]; }

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
  mfem::Array<int> GetBdrAttributeToMaterial() const;

  const auto &GetAttributeGlobalToLocal() const { return loc_attr; }

  const auto &GetBdrAttributeGlobalToLocal() const { return loc_bdr_attr; }

  template <typename T>
  auto GetAttributeGlobalToLocal(const T &attr_list) const
  {
    // Skip any entries in the input global attribute list which are not on local to this
    // process.
    const auto &loc_attr = GetAttributeGlobalToLocal();
    mfem::Array<int> loc_attr_list;
    for (auto attr : attr_list)
    {
      if (loc_attr.find(attr) != loc_attr.end())
      {
        loc_attr_list.Append(loc_attr.at(attr));
      }
    }
    return loc_attr_list;
  }

  template <typename T>
  auto GetBdrAttributeGlobalToLocal(const T &attr_list) const
  {
    // Skip any entries in the input global boundary attribute list which are not on local
    // to this process.
    const auto &loc_bdr_attr = GetBdrAttributeGlobalToLocal();
    mfem::Array<int> loc_attr_list;
    for (auto attr : attr_list)
    {
      if (loc_bdr_attr.find(attr) != loc_bdr_attr.end())
      {
        const auto &bdr_attr_map = loc_bdr_attr.at(attr);
        for (auto it = bdr_attr_map.begin(); it != bdr_attr_map.end(); ++it)
        {
          loc_attr_list.Append(it->second);
        }
      }
    }
    return loc_attr_list;
  }

  auto GetAttributeGlobalToLocal(const int attr) const
  {
    return GetAttributeGlobalToLocal(std::vector<int>{attr});
  }

  auto GetBdrAttributeGlobalToLocal(const int attr) const
  {
    return GetBdrAttributeGlobalToLocal(std::vector<int>{attr});
  }

  int GetAttributeGlobalToLocal(mfem::ElementTransformation &T) const;

  const auto &GetMesh() const { return mesh; }
};

//
// Material property represented as a piecewise constant coefficient over mesh elements. Can
// be scalar-valued or matrix-valued.
//
class MaterialPropertyCoefficient : public mfem::Coefficient, public mfem::MatrixCoefficient
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Map attribute to material index (coeff = mat_coeff[attr_mat[attr - 1]], for 1-based
  // attributes).
  mfem::Array<int> attr_mat;

  // Material properry coefficients, ordered by material index.
  mfem::DenseTensor mat_coeff;

public:
  MaterialPropertyCoefficient(const MaterialOperator &mat_op)
    : mfem::MatrixCoefficient(0, 0), mat_op(mat_op)
  {
  }
  MaterialPropertyCoefficient(const MaterialOperator &mat_op,
                              const mfem::Array<int> &attr_mat_,
                              const mfem::DenseTensor &mat_coeff_, double a = 1.0);

  bool empty() const { return mat_coeff.TotalSize() == 0; }

  const auto &GetAttributeToMaterial() const { return attr_mat; }
  const auto &GetMaterialProperties() const { return mat_coeff; }

  void AddCoefficient(const mfem::Array<int> &attr_mat_,
                      const mfem::DenseTensor &mat_coeff_, double a = 1.0);

  template <typename T>
  void AddMaterialProperty(const mfem::Array<int> &attr_list, const T &coeff,
                           double a = 1.0);
  template <typename T>
  void AddMaterialProperty(int attr, const T &coeff, double a = 1.0)
  {
    mfem::Array<int> attr_list(1);
    attr_list[0] = attr;
    AddMaterialProperty(attr_list, coeff, a);
  }

  void RestrictCoefficient(const mfem::Array<int> &attr_list);

  void NormalProjectedCoefficient(const mfem::Vector &normal);

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override;
};

}  // namespace palace

#endif  // PALACE_MODELS_MATERIAL_OPERATOR_HPP
