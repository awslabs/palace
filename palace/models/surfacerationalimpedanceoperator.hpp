// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP

#include <complex>
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
// SurfaceImpedanceOperator (which maps onto the constant stiffness/damping/mass matrices),
// a general Zs(ω) is an arbitrary function of frequency and therefore contributes to the
// frequency-dependent "extra" system matrix A2(ω), exactly like the finite-conductivity
// Robin BC (SurfaceConductivityOperator). It is consequently available only for
// frequency-domain driven and boundary mode problems.
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
    // Long division s·D(s) = P(s)·N(s) + R(s) with deg(R) < deg(N), splitting the Robin
    // coefficient g(s) = s·D/N = P + R/N into a polynomial part P (exactly representable
    // in a polynomial eigenvalue pencil) and a strictly proper remainder R/N (the pole
    // part). Highest-degree-first; an empty vector is the zero polynomial.
    std::vector<double> robin_quotient, robin_remainder;
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

  // Number of configured rational impedance boundaries.
  int GetNumBoundaries() const { return static_cast<int>(boundaries.size()); }

  // Stamp a unit coefficient (including the per-attribute crack scaling) for boundary
  // index idx, so that the assembled boundary mass matrix M_b satisfies: contribution of
  // boundary idx to A2(λ) = g(λ)·M_b with g = EvalRobinCoef(idx, ·). Used by the nonlinear
  // eigensolver seed to fit or freeze the boundary's Robin coefficient.
  void AddUnitBdrCoefficient(int idx, MaterialPropertyCoefficient &fb) const;

  // Evaluate the scalar Robin coefficient g(s) = s·D(s)/N(s) of boundary idx at complex
  // s = iω = λ (nondimensional).
  std::complex<double> EvalRobinCoef(int idx, std::complex<double> s) const;

  // Evaluate the strictly proper (pole) part R(s)/N(s) of the Robin coefficient of
  // boundary idx, from the long division g = P + R/N. The polynomial part P = g - R/N is
  // exactly representable in a polynomial eigenvalue pencil whenever deg(P) <= 2.
  std::complex<double> EvalRobinRemainder(int idx, std::complex<double> s) const;

  // Effective degree of the polynomial part P of the Robin coefficient of boundary idx
  // (-1 for the zero polynomial).
  int GetRobinQuotientDegree(int idx) const;

  // Add contributions to the frequency-dependent system matrix A2(ω). The Robin BC term has
  // coefficient iω / Zs(iω) per square, exactly as the parallel-RLC admittance contributes
  // iω·Ys to K + iωC - ω²M.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Complex-ω overload: evaluates the Robin coefficient s·D(s)/N(s) at genuinely complex
  // s = iω (analytic continuation; rational functions continue trivially off the real-ω
  // axis). Used by the 2D mode assembly when the wave-port cross-section EVP is solved at
  // a complex frequency. For real ω this reduces bit-for-bit to the double overload. The
  // passivity warning applies only on the real-ω axis and is skipped for Im(ω) != 0.
  void AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                     MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RATIONAL_IMPEDANCE_OPERATOR_HPP
