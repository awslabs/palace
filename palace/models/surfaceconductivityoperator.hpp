// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

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
  // and permeability, and (optionally) thickness. The attributes are split into those on
  // the exterior of the domain (or on a cracked interior boundary), where the boundary
  // condition applies to one side, and those on an interior boundary (a conducting sheet
  // with two surfaces), where the surface admittance of both surfaces is applied. A
  // boundary with "External" and a finite thickness is split into two groups when it has
  // both kinds of attributes, since the thickness correction for exterior surfaces does not
  // apply to the interior ones. The optional permittivity of the conductor is ε∞ = eps with
  // the frequency-dependent terms.
  struct ConductivityData
  {
    double sigma, eps, mu, h;
    std::vector<config::PermittivityTermData> terms;
    mfem::Array<int> attr_list, ext_attr_list, int_attr_list;
  };
  std::vector<ConductivityData> boundaries;

  void SetUpBoundaryProperties(const std::vector<config::ConductivityData> &conductivity,
                               ProblemType problem_type, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

public:
  SurfaceConductivityOperator(const std::vector<config::ConductivityData> &conductivity,
                              ProblemType problem_type, const Units &units,
                              const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  SurfaceConductivityOperator(const IoData &iodata, const MaterialOperator &mat_op,
                              const mfem::ParMesh &mesh);

  // Returns array of finite conductivity boundary attributes.
  mfem::Array<int> GetAttrList() const;

  // Mesh boundary attributes of attribute group g.
  const mfem::Array<int> &GetAttrList(std::size_t g) const
  {
    return boundaries[g].attr_list;
  }

  // Number of finite-conductivity attribute groups (each with its own σ/μ/h and thus its
  // own frequency-dependent surface admittance i·ω/Z(ω)).
  std::size_t Size() const { return boundaries.size(); }

  // Whether group `g` is active (nonzero conductivity or permittivity).
  bool IsActive(std::size_t g) const
  {
    return g < boundaries.size() && (boundaries[g].sigma > 0.0 || boundaries[g].eps > 0.0);
  }

  // Whether a group has a conductor permittivity term, defined at real frequencies only.
  bool HasRealFrequencyTerms() const;

  // Per-group surface-admittance scalar i·ω/Z_g(ω), evaluated at (possibly complex) ω. Z
  // is the surface impedance of the conductor with the optional finite-thickness
  // correction, analytic in ω (conductor permittivity terms at real ω only).
  std::complex<double> EvaluateScalar(std::size_t group_idx,
                                      std::complex<double> omega) const;

  // Add the ω-independent unit-coefficient boundary mass for group `g`, times coeff, to a
  // material property coefficient. The full per-group contribution to the system matrix at
  // ω is `(i·ω/Z_g(ω))·(this)` for coeff = 1, mirroring
  // WavePortOperator::AddBoundaryMassBdrCoefficients.
  void AddBoundaryMassBdrCoefficients(std::size_t g, MaterialPropertyCoefficient &fb,
                                      double coeff = 1.0) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
