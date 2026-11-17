// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_TERMS_HPP
#define PALACE_MODELS_SURFACE_TERMS_HPP

#include <complex>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "models/farfieldboundaryoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

namespace config
{

struct BoundaryData;

}  // namespace config

enum class ProblemType : char;

//
// The surface boundary conditions of one mesh (the 3D mesh, the face mesh of a wave port,
// or a 2D cross-section): absorbing (farfield), finite conductivity, impedance, and
// rational impedance boundaries. Each method adds the terms of all of these boundaries, so
// an assembler which calls it has the terms of every type.
//
// Inside this class, each type is named in each method to which it contributes. To add a
// surface boundary type, edit the member and both constructors, GetAttrList, HasDamping if
// it has a damping term, each Add* method to which it contributes a term, and if a term
// depends on frequency, the list of frequency-dependent terms (sigma_terms and the
// ExtraSystemBdrMassTerm methods, from which the summed ExtraSystemBdrMass methods are
// built). Outside this class, SpaceOperator names each type for its 3D boundary checks
// (CheckBoundaryProperties), and config::BoundaryData::GetCrossSectionEdgeAttributes for
// the boundary edges of 2D cross-sections.
//
class SurfaceTerms
{
private:
  FarfieldBoundaryOperator farfield_op;
  SurfaceConductivityOperator surf_sigma_op;
  SurfaceImpedanceOperator surf_z_op;
  SurfaceRationalImpedanceOperator surf_rz_op;

  // Surface conductivity groups with nonzero conductivity, the first of the
  // frequency-dependent boundary mass terms.
  std::vector<std::size_t> sigma_terms;

  // Add the frequency-dependent boundary mass terms at the (possibly complex) frequency ω,
  // split into real and imaginary parts, to the boundary coefficient of a vector field or
  // of a scalar field (see AddExtraSystemBdrMassCoefficients).
  void AddExtraSystemBdrMassImpl(std::complex<double> omega, bool scalar,
                                 MaterialPropertyCoefficient &fbr,
                                 MaterialPropertyCoefficient &fbi);

public:
  // The boundary attributes are checked on mesh (which also gives the surface normals of
  // the printed summary), and mat_op maps them to the boundary elements of the mesh which
  // is assembled: for a wave port, mesh is the parent 3D mesh and mat_op the material
  // operator of the face mesh of the port. The constructor from IoData also warns when the
  // second-order absorbing boundary is not passive for the Floquet wave vector
  // (FarfieldBoundaryOperator::CheckFloquetPassivity).
  SurfaceTerms(const config::BoundaryData &boundaries, ProblemType problem_type,
               const Units &units, const MaterialOperator &mat_op,
               const mfem::ParMesh &mesh);
  SurfaceTerms(const IoData &iodata, const MaterialOperator &mat_op,
               const mfem::ParMesh &mesh);

  // Returns the attributes of all surface boundaries.
  mfem::Array<int> GetAttrList() const;

  // Whether there are damping terms (resistive impedance or absorbing boundaries).
  bool HasDamping() const;

  // Add the stiffness (inductance), damping (resistance and first-order absorbing
  // boundaries), and mass (capacitance) terms to the boundary coefficient of a vector field
  // (the tangential trace).
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add the damping terms to the boundary coefficient of a scalar field (the out-of-plane
  // field component of a 2D mode solve): as AddDampingBdrCoefficients, but with the scalar
  // inverse impedance of the material next to the absorbing boundaries.
  //
  // Both damping methods add the absorbing terms before the impedance terms. The order
  // changes a coefficient only where MaterialPropertyCoefficient::AddMaterialProperty finds
  // a stored coefficient equal to a new one to a relative 1e-9, which then takes the new
  // value for the boundaries of both.
  void AddScalarDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add the frequency-dependent boundary mass terms (finite conductivity and rational
  // impedance boundaries) at the frequency ω, split into real and imaginary parts, to the
  // boundary coefficient of a vector field (the tangential trace) or of a scalar field (the
  // out-of-plane field component of a 2D mode solve), as the damping methods. These are
  // the frequency-dependent terms of a 2D mode solve. They are the terms of the list
  // below, each its unit coefficient times the real part of its scale at ω in fbr, then
  // times the imaginary part in fbi. At a real frequency, a rational impedance boundary
  // which is not passive warns once.
  void AddExtraSystemBdrMassCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                         MaterialPropertyCoefficient &fbi);
  void AddScalarExtraSystemBdrMassCoefficients(double omega,
                                               MaterialPropertyCoefficient &fbr,
                                               MaterialPropertyCoefficient &fbi);

  // Add all frequency-dependent terms at the frequency ω: the boundary mass terms of
  // AddExtraSystemBdrMassCoefficients (fb) and the second-order absorbing boundary terms
  // (see FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients): boundary curl-curl
  // (dfb), rotated boundary mass (frb), and mixed boundary curl (fpb, dropped when null).
  // The complex-ω overload is the analytic continuation, which reduces to the other one for
  // real ω.
  void AddExtraSystemBdrCoefficients(
      double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
      MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
      MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
      MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi);
  void AddExtraSystemBdrCoefficients(
      std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
      MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
      MaterialPropertyCoefficient &fbi, MaterialPropertyCoefficient &frbr,
      MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
      MaterialPropertyCoefficient *fpbi);

  // The frequency-dependent boundary mass terms one at a time, for the index
  // 0 <= i < GetNumExtraSystemBdrMassTerms(): each a boundary mass with a unit coefficient
  // which is scaled by a function of the (possibly complex) frequency ω. The terms are one
  // for each surface conductivity group with nonzero conductivity, scaled by i ω / Z(ω),
  // then one for each rational impedance boundary, scaled by i ω / Zs(iω). A term gives
  // the unit coefficient of a vector field (the tangential trace, as
  // AddDampingBdrCoefficients) and that of a scalar field (the out-of-plane field
  // component of a 2D mode solve, as AddScalarDampingBdrCoefficients), times coeff; for
  // these types the two are the same.
  int GetNumExtraSystemBdrMassTerms() const;
  void AddExtraSystemBdrMassTermCoefficients(int i, MaterialPropertyCoefficient &fb,
                                             double coeff = 1.0) const;
  void AddScalarExtraSystemBdrMassTermCoefficients(int i, MaterialPropertyCoefficient &fb,
                                                   double coeff = 1.0) const;
  std::complex<double> GetExtraSystemBdrMassTermScale(int i,
                                                      std::complex<double> omega) const;

  // The label of term i in output files and log lines (surfsigma_<group> and
  // rationalz_<boundary>), and the name of its boundary in messages, for example
  // "rational impedance boundary (attribute 5)".
  std::string GetExtraSystemBdrMassTermLabel(int i) const;
  std::string GetExtraSystemBdrMassTermName(int i) const;

  // An exact split of the scale of a term, as a function c(s) of s = iω (the scale at ω is
  // c(iω)): c(s) = P(s) + R(s), with a polynomial P and a strictly proper rational
  // remainder R. GetExtraSystemBdrMassTermSplit returns it for the terms which have one:
  // a rational impedance boundary, for which it is the long division of the Robin
  // coefficient s D(s) / N(s), and not a surface conductivity group. Its functions refer
  // to this SurfaceTerms.
  struct ExactSplit
  {
    // Coefficients of P, highest degree first (empty for the zero polynomial), and the
    // degree of P without leading zero coefficients (-1 for the zero polynomial).
    std::vector<double> polynomial;
    int degree;

    // c(s) and R(s).
    std::function<std::complex<double>(std::complex<double>)> coefficient, remainder;
  };
  std::optional<ExactSplit> GetExtraSystemBdrMassTermSplit(int i) const;

  // Access to the operators of each type, for operations which only one solve has.
  auto &GetFarfieldOp() { return farfield_op; }
  auto &GetRationalImpedanceOp() { return surf_rz_op; }
  const auto &GetFarfieldOp() const { return farfield_op; }
  const auto &GetSurfaceConductivityOp() const { return surf_sigma_op; }
  const auto &GetSurfaceImpedanceOp() const { return surf_z_op; }
  const auto &GetRationalImpedanceOp() const { return surf_rz_op; }
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_TERMS_HPP
