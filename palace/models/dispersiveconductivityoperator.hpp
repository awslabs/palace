// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DISPERSIVE_CONDUCTIVITY_OPERATOR_HPP
#define PALACE_MODELS_DISPERSIVE_CONDUCTIVITY_OPERATOR_HPP

#include <complex>
#include <vector>
#include "models/complexconductivitymodel.hpp"

namespace palace
{

class MaterialOperator;
class MaterialPropertyCoefficient;
class Units;

namespace config
{

struct MaterialData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling the frequency-dependent part of dispersive complex conductivity domains.
//
// A domain with conductivity σ(ω) adds iω σ(ω) M_g to the system matrix, with M_g the mass
// matrix of its group g of attributes. MaterialOperator includes the frequency-independent
// part iω σ_a(ω) = iω σ1(ω_ref) + ω_ref σ2(ω_ref), the two-fluid model matching σ at the
// reference frequency, in the damping and stiffness matrices. This class adds the residual
// δy_g(ω) = iω (σ(ω) - σ_a(ω)) M_g, which vanishes at ω_ref, to A2(ω).
//
class DispersiveConductivityOperator
{
private:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  struct DispersiveGroup
  {
    ComplexConductivityModel model;
    std::vector<int> attributes;
    double omega_ref;
    std::complex<double> sigma_ref;
  };
  std::vector<DispersiveGroup> groups;

  void PrintDomainInfo(const Units &units) const;

public:
  DispersiveConductivityOperator(const std::vector<config::MaterialData> &materials,
                                 ProblemType problem_type, const Units &units,
                                 const MaterialOperator &mat_op);

  // Number of dispersive material groups (the same on all processes).
  std::size_t Size() const { return groups.size(); }
  bool Empty() const { return groups.empty(); }

  // Residual coefficient δy_g(ω) of group g.
  std::complex<double> EvaluateScalar(std::size_t g, double omega) const;

  // Add the residual coefficients of all groups to the real and imaginary parts of the A2
  // domain coefficient.
  void AddExtraSystemCoefficients(double omega, MaterialPropertyCoefficient &fr,
                                  MaterialPropertyCoefficient &fi) const;

  // Add a unit coefficient on the attributes of group g, for the mass matrix M_g.
  void AddUnitCoefficients(std::size_t g, MaterialPropertyCoefficient &f) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DISPERSIVE_CONDUCTIVITY_OPERATOR_HPP
