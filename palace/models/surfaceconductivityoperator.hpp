// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

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
public:
  // Surface properties for finite conductivity boundary attributes: conductor conductivity
  // and permeability, and (optionally) thickness. Public so the analytic-continuation
  // path (SurfSigmaFactor) can iterate per-group and obtain the per-group scalar
  // coefficient at complex λ.
  struct ConductivityData
  {
    double sigma, mu, h;
    mfem::Array<int> attr_list;
  };

private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

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

  // Add contributions to system matrix for a finite conductivity boundary condition.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Complex-ω overload (analytic continuation for the complex-frequency cross-section
  // EVP). Forwards to EvaluateScalar with complex ω; for real ω (imag = 0) it produces
  // the same (fbr, fbi) split as the double overload.
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Per-group access for the analytic-continuation path. NumGroups returns the number
  // of distinct (σ, μ, h) attribute groups; AddBoundaryMassBdrCoefficients adds the
  // ω-independent boundary mass for a single group with optional scalar `coeff`;
  // EvaluateScalar returns the per-group complex scalar coefficient i·ω/Z(ω) at
  // complex ω — the analytic-continuation argument is supplied via ω = -i·λ at the
  // call site. Real-ω evaluation also routes through EvaluateScalar (with imag(ω) = 0)
  // so the formula lives in one place.
  std::size_t NumGroups() const { return boundaries.size(); }
  void AddBoundaryMassBdrCoefficients(std::size_t group_idx,
                                      MaterialPropertyCoefficient &fb,
                                      double coeff = 1.0) const;
  std::complex<double> EvaluateScalar(std::size_t group_idx,
                                      std::complex<double> omega) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CONDUCTIVITY_OPERATOR_HPP
