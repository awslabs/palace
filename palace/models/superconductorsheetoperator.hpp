// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP
#define PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

//
// A class handling thin-film superconductor sheet boundaries. Each sheet contributes a
// tangential surface term (1/L_ksq) A_t · v_t to the curl-curl operator, where the kinetic
// sheet inductance is L_ksq = mu0 * lambda^2 / d (nondimensionally lambda^2 / d). This is
// the London kinetic-inductance contribution for a superconducting film modeled as a 2D
// sheet.
//
class SuperconductorSheetOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for superconductor sheet attributes: kinetic sheet inductance L_ksq
  // [H/sq] and the per-attribute area scaling for mesh cracking.
  struct SuperconductorSheetData
  {
    double Ls;
    mfem::Array<int> attr_list;
    std::unordered_map<int, double> attr_scaling;
  };
  std::vector<SuperconductorSheetData> boundaries;

  void
  SetUpBoundaryProperties(const std::vector<config::SuperconductorData> &superconductor,
                          const std::unordered_set<int> &cracked_attributes,
                          const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

public:
  SuperconductorSheetOperator(const std::vector<config::SuperconductorData> &superconductor,
                              const std::unordered_set<int> &cracked_attributes,
                              const Units &units, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);
  SuperconductorSheetOperator(const IoData &iodata, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);

  // Returns array of superconductor sheet attributes.
  mfem::Array<int> GetAttrList() const;

  // Returns true if no superconductor sheets are configured.
  bool empty() const { return boundaries.empty(); }

  // Add the tangential surface mass contribution (1/L_ksq) A_t · v_t to the system matrix.
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP
