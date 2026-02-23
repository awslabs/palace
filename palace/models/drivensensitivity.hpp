// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DRIVEN_SENSITIVITY_HPP
#define PALACE_MODELS_DRIVEN_SENSITIVITY_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>
#include "linalg/ksp.hpp"

namespace palace
{

class ComplexOperator;
class IoData;
class SpaceOperator;

//
// Computes S-parameter sensitivities with respect to piecewise-constant material
// properties (permittivity, inverse permeability, conductivity) for each material domain.
//
// For the frequency-domain system A(omega) u = f, the S-parameter S_ij is a linear
// functional of the solution: S_ij = s_j^H u_i, where s_j is the port mode linear form
// for observation port j.
//
// The adjoint system is:
//   A^H v_j = s_j
//
// Then the sensitivity with respect to material parameter p in domain d is:
//   dS_ij/dp = -v_j^H (dA/dp) u_i
//
// where dA/dp is the domain-restricted bilinear form (reuses P4 infrastructure):
//   dA/d(eps_d)   = -omega^2 * dM_d   (from the mass matrix)
//   dA/d(mu_inv_d) = dK_d              (from the stiffness matrix)
//   dA/d(sigma_d) = i*omega * dC_d     (from the damping matrix)
//
class DrivenSensitivity
{
public:
  // Per-frequency, per-observation-port, per-domain complex sensitivity of S_ij.
  struct SensitivityResult
  {
    // The global mesh attributes of the material domains.
    std::vector<int> domain_attrs;

    // dS/d(epsilon_r) for each (observation_port, domain) pair.
    // Outer index: observation port (matching port_indices order).
    // Inner index: domain index into domain_attrs.
    std::vector<std::vector<std::complex<double>>> dS_d_epsilon;

    // dS/d(mu_inv_r) for each (observation_port, domain) pair.
    std::vector<std::vector<std::complex<double>>> dS_d_muinv;

    // dS/d(sigma) for each (observation_port, domain) pair.
    std::vector<std::vector<std::complex<double>>> dS_d_sigma;

    // Port indices corresponding to the outer dimension of the sensitivity arrays.
    std::vector<int> port_indices;
  };

  // Compute S-parameter sensitivities dS_ij/dp for a given excitation port i at a given
  // frequency omega. The forward solution E (= u_i) must be precomputed.
  //
  // Parameters:
  //   iodata       - configuration data (for material domain info)
  //   space_op     - spatial discretization operator (for domain-restricted assembly)
  //   ksp          - linear solver (reused for the adjoint solve with A^H)
  //   A            - the assembled system matrix at this frequency
  //   P            - the preconditioner matrix at this frequency
  //   E            - forward solution for excitation port i
  //   excitation_idx - index of the excitation port
  //   omega        - angular frequency
  static SensitivityResult Compute(const IoData &iodata, SpaceOperator &space_op,
                                   ComplexKspSolver &ksp, const ComplexOperator &A,
                                   const ComplexOperator &P, const ComplexVector &E,
                                   int excitation_idx, double omega);
};

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SENSITIVITY_HPP
