// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP
#define PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP

#include <memory>
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
// sheet inductance is the finite-thickness London value L_ksq = mu0 * lambda *
// coth(d/lambda) (nondimensionally lambda * coth(d/lambda)). This is the London
// kinetic-inductance contribution for a superconducting film modeled as a 2D sheet; it
// reduces to the thin-film Pearl limit lambda^2/d for d << lambda and saturates at lambda
// for d >> lambda.
//
// A two-sided (two-port) sheet instead models the film as two coupled coincident faces (a
// cracked interface), so the tangential field can jump across it. Each face carries the
// diagonal self-term lambda*coth(d/lambda) and the faces are coupled by the cross-face term
// -1/(mu0*lambda*sinh(d/lambda)) (see BuildTwoPortCoupling), capturing two-sided screening
// for d >~ lambda.
//
class SuperconductorSheetOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for superconductor sheet attributes: kinetic sheet inductance L_ksq
  // [H/sq] and the per-attribute area scaling for mesh cracking. For a two-sided (two-port)
  // sheet, Ls holds the diagonal self-term inductance lambda*tanh(d/lambda) and (lambda_L,
  // d) are retained to build the cross-face coupling -1/(lambda*sinh(d/lambda)).
  struct SuperconductorSheetData
  {
    double Ls;
    double lambda_L = 0.0;
    double thickness = 0.0;
    bool two_sided = false;
    mfem::Array<int> attr_list;
    std::unordered_map<int, double> attr_scaling;
  };
  std::vector<SuperconductorSheetData> boundaries;
  bool has_two_port_ = false;

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

  // Finite-thickness London kinetic sheet inductance L_ksq = lambda * coth(d/lambda) from
  // the penetration depth and film thickness (in consistent units; mu0 absorbed in the
  // nondimensional convention). Reduces to lambda^2/d for d << lambda, saturates at lambda
  // for d >> lambda.
  static double KineticSheetInductance(double lambda_L, double thickness);

  // Returns array of superconductor sheet attributes.
  mfem::Array<int> GetAttrList() const;

  // Returns true if no superconductor sheets are configured.
  bool empty() const { return boundaries.empty(); }

  // Returns true if any sheet is a two-sided (two-port) film needing the cross-face
  // coupling.
  bool HasTwoPort() const { return has_two_port_; }

  // Add the tangential surface mass contribution (1/L_ksq) A_t · v_t to the system matrix.
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb) const;

  // Build the cross-face coupling matrix -1/(mu0*lambda*sinh(d/lambda)) * ∫_Σ A_t^+ · v_t^-
  // for two-sided sheets, on the true dofs of the given ND space (serial). Pairs the two
  // coincident cracked faces by centroid and integrates the cross mass by evaluating each
  // face's Piola basis at shared physical quadrature points. Returns nullptr if none.
  std::unique_ptr<mfem::SparseMatrix>
  BuildTwoPortCoupling(mfem::ParFiniteElementSpace &nd_fespace) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP
