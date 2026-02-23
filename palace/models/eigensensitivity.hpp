// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_EIGEN_SENSITIVITY_HPP
#define PALACE_MODELS_EIGEN_SENSITIVITY_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class ComplexOperator;
class EigenvalueSolver;
class IoData;
class SpaceOperator;

//
// Computes eigenfrequency sensitivities with respect to piecewise-constant material
// properties (permittivity, inverse permeability, conductivity) for each material domain.
//
// For the polynomial eigenvalue problem F(lambda) = (K + lambda*C + lambda^2*M) u = 0,
// the sensitivity of eigenvalue lambda to a parameter p is:
//
//   dlambda/dp = -u^H (dF/dp) u / u^H (dF/dlambda) u
//
// where dF/dp is the derivative of the operator polynomial with respect to p, and
// dF/dlambda = C + 2*lambda*M is the derivative with respect to the eigenvalue.
//
// For the linear EVP (C = 0), the problem reduces to K u = mu*M u where mu = omega^2,
// and the sensitivity simplifies to:
//
//   dmu/dp = u^H (dK/dp - mu * dM/dp) u / (u^H M u)
//
class EigenSensitivity
{
public:
  // Per-mode, per-domain complex sensitivity of eigenfrequency f = omega / (2*pi).
  // Outer index: mode (0..num_conv-1). Inner index: domain index into domain_attrs.
  struct SensitivityResult
  {
    // The global mesh attributes of the material domains.
    std::vector<int> domain_attrs;

    // df/d(epsilon_r) for each (mode, domain) pair.
    std::vector<std::vector<std::complex<double>>> df_d_epsilon;

    // df/d(mu_inv_r) for each (mode, domain) pair.
    std::vector<std::vector<std::complex<double>>> df_d_muinv;

    // df/d(sigma) for each (mode, domain) pair.
    std::vector<std::vector<std::complex<double>>> df_d_sigma;
  };

  // Compute eigenfrequency sensitivities for all converged eigenpairs. The eigenvectors
  // are assumed to be M-normalized (u^H M u = 1), which Palace enforces via
  // RescaleEigenvectors.
  //
  // Parameters:
  //   iodata      - configuration data (for material domain info)
  //   space_op    - spatial discretization operator (for domain-restricted assembly)
  //   eigen       - eigenvalue solver with converged results
  //   num_conv    - number of converged eigenpairs
  //   has_damping - whether the problem is quadratic (C or A2 present)
  //   C           - assembled damping matrix (may be nullptr for linear EVP)
  //   M           - assembled mass matrix (for denominator computation)
  static SensitivityResult Compute(const IoData &iodata, SpaceOperator &space_op,
                                   const EigenvalueSolver &eigen, int num_conv,
                                   bool has_damping, const ComplexOperator *C,
                                   const ComplexOperator *M);
};

}  // namespace palace

#endif  // PALACE_MODELS_EIGEN_SENSITIVITY_HPP
