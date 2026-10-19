// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP

#include <complex>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include "models/complexconductivitymodel.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

namespace config
{

struct ConductivityData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling finite conductivity boundaries.
//
class SurfaceConductivityOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Surface properties for finite conductivity boundary attributes: conductor conductivity
  // and London penetration depth, or a dispersive model of σ(ω), permeability, and
  // (optionally) thickness.
  struct ConductivityData
  {
    double sigma, lambda, mu, h;
    std::optional<ComplexConductivityModel> model;
    mfem::Array<int> attr_list;
  };
  std::vector<ConductivityData> boundaries;

  void SetUpBoundaryProperties(const std::vector<config::ConductivityData> &conductivity,
                               ProblemType problem_type, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

  // Per-group complex surface-admittance coefficient i·ω/Z(ω) for the boundary BC,
  // evaluated at (possibly complex) ω. Z is the surface impedance of the conductor with the
  // optional finite-thickness correction. Dispersive models are only evaluated at real ω.
  std::complex<double> EvaluateScalarImpl(std::size_t group_idx,
                                          std::complex<double> omega) const;

public:
  SurfaceConductivityOperator(const std::vector<config::ConductivityData> &conductivity,
                              ProblemType problem_type, const Units &units,
                              const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  SurfaceConductivityOperator(const IoData &iodata, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);

  // Returns array of finite conductivity boundary attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix for a finite conductivity boundary condition.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Complex-ω overload, used by the eigenmode nonlinear solve where the surface-admittance
  // coefficient i·ω/Z(ω) is evaluated at the genuinely complex eigenvalue (ω = -i·λ). For
  // real ω (imag = 0) this produces the same (fbr, fbi) split as the double overload.
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Number of finite-conductivity attribute groups (each with its own σ/μ/h and thus its
  // own frequency-dependent surface admittance i·ω/Z(ω)).
  std::size_t Size() const { return boundaries.size(); }

  // Whether group `g` is active (nonzero conductivity or London depth, or a dispersive
  // model).
  bool IsActive(std::size_t g) const
  {
    return g < boundaries.size() &&
           (boundaries[g].sigma > 0.0 || boundaries[g].lambda > 0.0 || boundaries[g].model);
  }

  // Per-group surface-admittance scalar i·ω/Z_g(ω), evaluated at (possibly complex) ω.
  // Public wrapper so the circuit-synthesis path can sample the scalar for its dispersion
  // fit (the ω-independent boundary mass is factored out via
  // SpaceOperator::GetSurfaceConductivityBoundaryMatrix).
  std::complex<double> EvaluateScalar(std::size_t group_idx,
                                      std::complex<double> omega) const
  {
    return EvaluateScalarImpl(group_idx, omega);
  }

  // Add the ω-independent unit-coefficient boundary mass for group `g` to a material
  // property coefficient. The full per-group contribution to the system matrix at ω is
  // `(i·ω/Z_g(ω))·(this)`. Used to factor out the ω-independent operator for circuit
  // synthesis, mirroring WavePortOperator::AddBoundaryMassBdrCoefficients.
  void AddBoundaryMassBdrCoefficients(std::size_t g, MaterialPropertyCoefficient &fb) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
