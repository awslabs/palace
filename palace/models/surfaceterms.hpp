// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_TERMS_HPP
#define PALACE_MODELS_SURFACE_TERMS_HPP

#include <complex>
#include <mfem.hpp>
#include "models/farfieldboundaryoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"

namespace palace
{

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
// surface boundary type, edit the member and the constructor, GetAttrList, and each Add*
// method to which it contributes a term. Outside this class, SpaceOperator also names each
// type for its 3D boundary checks (CheckBoundaryProperties).
//
class SurfaceTerms
{
private:
  FarfieldBoundaryOperator farfield_op;
  SurfaceConductivityOperator surf_sigma_op;
  SurfaceImpedanceOperator surf_z_op;
  SurfaceRationalImpedanceOperator surf_rz_op;

  // Add the frequency-dependent boundary mass terms (finite conductivity and rational
  // impedance boundaries) at the real or complex frequency ω.
  template <typename OmegaType>
  void AddExtraSystemBdrMassImpl(OmegaType omega, MaterialPropertyCoefficient &fbr,
                                 MaterialPropertyCoefficient &fbi);

public:
  SurfaceTerms(const config::BoundaryData &boundaries, ProblemType problem_type,
               const Units &units, const MaterialOperator &mat_op,
               const mfem::ParMesh &mesh);

  // Returns the attributes of all surface boundaries.
  mfem::Array<int> GetAttrList() const;

  // Add the stiffness (inductance), damping (resistance and first-order absorbing
  // boundaries), and mass (capacitance) terms to the boundary coefficient of a vector field
  // (the tangential trace).
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add all frequency-dependent terms at the frequency ω: the boundary mass terms of the
  // finite conductivity and rational impedance boundaries (fb) and the second-order
  // absorbing boundary terms (see FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients):
  // boundary curl-curl (dfb), rotated boundary mass (frb), and mixed boundary curl (fpb,
  // dropped when null).
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

  // Access to the operators of each type, for operations which only one solve has.
  auto &GetFarfieldOp() { return farfield_op; }
  auto &GetSurfaceConductivityOp() { return surf_sigma_op; }
  auto &GetSurfaceImpedanceOp() { return surf_z_op; }
  auto &GetRationalImpedanceOp() { return surf_rz_op; }
  const auto &GetFarfieldOp() const { return farfield_op; }
  const auto &GetSurfaceConductivityOp() const { return surf_sigma_op; }
  const auto &GetSurfaceImpedanceOp() const { return surf_z_op; }
  const auto &GetRationalImpedanceOp() const { return surf_rz_op; }
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_TERMS_HPP
