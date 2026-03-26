// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#ifndef PALACE_MODELS_LUMPED_ELEMENT_OPERATOR_HPP
#define PALACE_MODELS_LUMPED_ELEMENT_OPERATOR_HPP
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/lumpedgeometry.hpp"
#include "utils/configfile.hpp"
namespace palace
{
class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;
namespace config
{
struct LumpedElementData;
}  // namespace config
//
// Helper class for passive lumped element boundaries (R/L/C load, no port tracking,
// no excitation).
//
class LumpedElementData
{
public:
  const MaterialOperator &mat_op;
  std::unique_ptr<LumpedGeometry> elem;
  double R, L, C;
  LumpedElementData(const config::LumpedElementData &data,
                    const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  double GetToSquare() const
  {
    return elem->GetGeometryWidth() / elem->GetGeometryLength();
  }
};
//
// A class handling passive lumped element boundaries and their assembly.
//
class LumpedElementOperator
{
private:
  std::map<int, LumpedElementData> elements;
  void SetUpBoundaryProperties(const std::map<int, config::LumpedElementData> &lumpedelement,
                               const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);
public:
  LumpedElementOperator(const std::map<int, config::LumpedElementData> &lumpedelement,
                        const Units &units, const MaterialOperator &mat_op,
                        const mfem::ParMesh &mesh);
  LumpedElementOperator(const IoData &iodata, const MaterialOperator &mat_op,
                        const mfem::ParMesh &mesh);
  const LumpedElementData &GetElement(int idx) const;
  auto begin() const { return elements.begin(); }
  auto end()   const { return elements.end(); }
  auto Size()  const { return elements.size(); }
  mfem::Array<int> GetAttrList()   const;
  mfem::Array<int> GetRsAttrList() const;
  mfem::Array<int> GetLsAttrList() const;
  mfem::Array<int> GetCsAttrList() const;
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients  (double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients     (double coeff, MaterialPropertyCoefficient &fb);
};
}  // namespace palace
#endif  // PALACE_MODELS_LUMPED_ELEMENT_OPERATOR_HPP