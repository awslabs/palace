// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP

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

enum class ProblemType : char;

//
// A class handling general rational surface impedance boundaries.
//
// The surface impedance per square is a user-provided rational function of frequency,
//
//     Zs(s) = N(s) / D(s),    s = iω,
//
// with N and D real polynomials given as coefficient lists (highest-degree-first). The
// roots of N are the zeros and the roots of D are the poles of Zs, so any passive lumped
// (RLC / network) response can be represented. Unlike the parallel-RLC impedance handled by
// SurfaceImpedanceOperator (which maps onto the constant stiffness/damping/mass matrices), a
// general Zs(ω) is an arbitrary function of frequency and therefore contributes to the
// frequency-dependent "extra" system matrix A2(ω), exactly like the finite-conductivity
// Robin BC (SurfaceConductivityOperator). It is consequently available only for
// frequency-domain problem types (driven, eigenmode, boundary mode).
//
class SurfaceRationalImpedanceOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Frequency scale factor [GHz], used to report physical frequencies in warnings.
  double freq_scale = 1.0;

  // Rational surface impedance per boundary: Zs(s) = N(s)/D(s), s = iω, with N and D real
  // polynomial coefficients stored nondimensionalized and highest-degree-first.
  struct RationalImpedanceData
  {
    std::vector<double> num, den;
    mfem::Array<int> attr_list;
    std::unordered_map<int, double> attr_scaling;
    bool warned_passivity = false;  // One-shot guard for the per-frequency passivity check.
  };
  std::vector<RationalImpedanceData> boundaries;

  void SetUpBoundaryProperties(const std::vector<config::RationalImpedanceData> &impedance,
                               const std::unordered_set<int> &cracked_attributes,
                               ProblemType problem_type, const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

public:
  SurfaceRationalImpedanceOperator(
      const std::vector<config::RationalImpedanceData> &impedance,
      const std::unordered_set<int> &cracked_attributes, ProblemType problem_type,
      const Units &units, const MaterialOperator &mat_op, const mfem::ParMesh &mesh);
  SurfaceRationalImpedanceOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                   const mfem::ParMesh &mesh);

  // Returns array of rational surface impedance attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to the frequency-dependent system matrix A2(ω). The Robin BC term has
  // coefficient iω / Zs(iω) per square, exactly as the parallel-RLC admittance contributes
  // iω·Ys to K + iωC - ω²M.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP
