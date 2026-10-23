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

struct CrackSides;
class FiniteElementSpace;
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
// A two-sided (two-port) sheet instead models the film as two coupled coincident faces, so
// the tangential field can jump across it: the solution space is broken across the film,
// and the two versions of the DOFs of each film face are the fields of its two sides. Each
// side carries the diagonal self-term λ·tanh(d/λ) (stiffness coth(d/λ)/(μ₀λ)) and the
// sides are coupled by the cross term −1/(μ₀λ·sinh(d/λ)), capturing two-sided screening
// for d ≳ λ (see AssembleBrokenSheetMatrix).
//
class SuperconductorSheetOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for superconductor sheet attributes: kinetic sheet inductance L_ksq
  // [H/sq] and the per-attribute area scaling for mesh cracking. For a two-sided sheet, Ls
  // holds the diagonal self-term inductance λ·tanh(d/λ), and (λ, d) give the cross term
  // −1/(λ·sinh(d/λ)).
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

  // Returns the boundary attributes of two-sided sheets, across which the solution space is
  // broken (empty without two-sided sheets).
  std::vector<int> GetTwoSidedAttributes() const;

  // Add the tangential surface mass contribution (1/L_ksq) A_t · v_t to the system matrix,
  // for all sheets but the two-sided ones (on a broken space, see
  // AssembleBrokenSheetMatrix instead).
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb) const;

  // Assemble the terms of all sheets as a matrix on the true DOFs of the broken space
  // broken_nd_fespace, a view of nd_fespace broken across GetTwoSidedAttributes() with the
  // given sides. On each two-sided film face with the fields u⁺, u⁻ of its two sides, the
  // bilinear form a (u⁺·v⁺ + u⁻·v⁻) − b (u⁺·v⁻ + u⁻·v⁺), with a = coth(d/λ)/λ and b =
  // 1/(λ·sinh(d/λ)), is assembled as
  //          (a + b) (u⁺·v⁺ + u⁻·v⁻) − b (u⁺ + u⁻)·(v⁺ + v⁻).
  // The first term is integrated by each element over its own film faces, with its own
  // DOFs. The second one only depends on u⁺ + u⁻, the field of nd_fespace whose true DOFs
  // sum the versions read by the two sides, so it is the film boundary mass of nd_fespace
  // between these sums. Neither term needs to identify the two sides of a face. The
  // one-sided sheets (1/L_ksq) u·v are integrated by the element adjacent to each of their
  // boundary elements with its own DOFs, so that a sheet standing on one side of a
  // two-sided film reads the versions of that side.
  std::unique_ptr<mfem::HypreParMatrix>
  AssembleBrokenSheetMatrix(const FiniteElementSpace &nd_fespace,
                            const FiniteElementSpace &broken_nd_fespace,
                            const CrackSides &sides) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SUPERCONDUCTOR_SHEET_OPERATOR_HPP
