// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensensitivity.hpp"

#include <cmath>
#include <complex>
#include "linalg/eps.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Compute the complex bilinear form u^H A u for a ComplexOperator A.
std::complex<double> BilinearForm(MPI_Comm comm, const ComplexOperator &A,
                                  const ComplexVector &u, ComplexVector &Au)
{
  A.Mult(u, Au);
  auto dot = linalg::Dot(comm, u, Au);
  return dot;
}

}  // namespace

EigenSensitivity::SensitivityResult
EigenSensitivity::Compute(const IoData &iodata, SpaceOperator &space_op,
                          const EigenvalueSolver &eigen, int num_conv, bool has_damping,
                          const ComplexOperator *C, const ComplexOperator *M)
{
  SensitivityResult result;
  MPI_Comm comm = space_op.GetComm();

  // Collect the list of global mesh attributes for all material domains.
  for (const auto &mat : iodata.domains.materials)
  {
    for (auto attr : mat.attributes)
    {
      result.domain_attrs.push_back(attr);
    }
  }
  const int num_domains = static_cast<int>(result.domain_attrs.size());

  if (num_conv <= 0 || num_domains <= 0)
  {
    return result;
  }

  // Allocate result arrays: [num_modes][num_domains].
  result.df_d_epsilon.resize(num_conv, std::vector<std::complex<double>>(num_domains, 0.0));
  result.df_d_muinv.resize(num_conv, std::vector<std::complex<double>>(num_domains, 0.0));
  result.df_d_sigma.resize(num_conv, std::vector<std::complex<double>>(num_domains, 0.0));

  // Working vectors sized to the local ND true dof count.
  const int n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u(n), Au(n);
  u.UseDevice(true);
  Au.UseDevice(true);

  Mpi::Print("\nComputing eigenvalue sensitivities with respect to material properties\n");

  // Pre-assemble domain-restricted operators outside the mode loop (they depend only on
  // the domain, not the eigenvector). This avoids redundant assembly for each mode.
  // Note: A2 (extra system matrix from wave ports/boundaries) is assumed to be independent
  // of bulk material properties (epsilon, mu, sigma) and is not differentiated here.
  std::vector<std::unique_ptr<ComplexOperator>> dK_vec(num_domains);
  std::vector<std::unique_ptr<ComplexOperator>> dM_vec(num_domains);
  std::vector<std::unique_ptr<ComplexOperator>> dC_vec(num_domains);
  for (int d = 0; d < num_domains; d++)
  {
    const int attr = result.domain_attrs[d];
    dK_vec[d] =
        space_op.GetStiffnessMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
    dM_vec[d] =
        space_op.GetMassMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
    dC_vec[d] =
        space_op.GetDampingMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
  }

  for (int i = 0; i < num_conv; i++)
  {
    // Retrieve eigenvalue and eigenvector.
    std::complex<double> eig_raw = eigen.GetEigenvalue(i);
    eigen.GetEigenvector(i, u);

    // Compute omega from the raw eigenvalue depending on problem type.
    std::complex<double> omega;
    std::complex<double> lambda;  // lambda = i*omega, as used in the polynomial EVP
    if (!has_damping)
    {
      // Linear EVP: eigenvalue is mu = omega^2.
      omega = std::sqrt(eig_raw);
      lambda = 1i * omega;
    }
    else
    {
      // Quadratic EVP: eigenvalue is lambda = i*omega.
      lambda = eig_raw;
      omega = lambda / 1i;
    }

    // Compute the denominator term. For the quadratic EVP:
    //   denom = u^H (C + 2*lambda*M) u
    // For the linear EVP:
    //   denom_mu = u^H M u  (should be 1.0 for M-normalized eigenvectors)
    std::complex<double> denom;
    if (has_damping)
    {
      // Quadratic case: denom = u^H C u + 2*lambda * (u^H M u)
      std::complex<double> uCu = 0.0;
      if (C)
      {
        uCu = BilinearForm(comm, *C, u, Au);
      }
      std::complex<double> uMu = 0.0;
      if (M)
      {
        uMu = BilinearForm(comm, *M, u, Au);
      }
      denom = uCu + 2.0 * lambda * uMu;
    }
    else
    {
      // Linear case: denom = u^H M u (should be ~1.0 for M-normalized vectors).
      if (M)
      {
        denom = BilinearForm(comm, *M, u, Au);
      }
      else
      {
        denom = 1.0;
      }
    }

    // Compute per-domain sensitivities using pre-assembled operators.
    for (int d = 0; d < num_domains; d++)
    {
      std::complex<double> uKdu = 0.0;
      if (dK_vec[d])
      {
        uKdu = BilinearForm(comm, *dK_vec[d], u, Au);
      }

      std::complex<double> uMdu = 0.0;
      if (dM_vec[d])
      {
        uMdu = BilinearForm(comm, *dM_vec[d], u, Au);
      }

      std::complex<double> uCdu = 0.0;
      if (dC_vec[d])
      {
        uCdu = BilinearForm(comm, *dC_vec[d], u, Au);
      }

      if (has_damping)
      {
        // Quadratic EVP: dlambda/dp = -(numerator) / denom
        // For epsilon sensitivity (from M term):
        //   dlambda/d(eps) = -(lambda^2 * uMdu) / denom
        //   df/d(eps) = (-i/(2*pi)) * dlambda/d(eps)
        std::complex<double> dlambda_deps = -(lambda * lambda * uMdu) / denom;
        result.df_d_epsilon[i][d] = (-1i / (2.0 * M_PI)) * dlambda_deps;

        // For mu_inv sensitivity (from K term):
        //   dlambda/d(mu_inv) = -(uKdu) / denom
        //   df/d(mu_inv) = (-i/(2*pi)) * dlambda/d(mu_inv)
        std::complex<double> dlambda_dmuinv = -(uKdu) / denom;
        result.df_d_muinv[i][d] = (-1i / (2.0 * M_PI)) * dlambda_dmuinv;

        // For sigma sensitivity (from C term):
        //   dlambda/d(sigma) = -(lambda * uCdu) / denom
        //   df/d(sigma) = (-i/(2*pi)) * dlambda/d(sigma)
        std::complex<double> dlambda_dsigma = -(lambda * uCdu) / denom;
        result.df_d_sigma[i][d] = (-1i / (2.0 * M_PI)) * dlambda_dsigma;
      }
      else
      {
        // Linear EVP: mu = omega^2, eigenvalue is mu (complex if loss tangent present).
        // dmu/d(eps) = -mu * uMdu / denom  (from -mu*dM/dp term)
        // dmu/d(mu_inv) = uKdu / denom     (from dK/dp term)
        // Then: domega/dp = dmu/dp / (2*omega), df/dp = domega/dp / (2*pi)
        std::complex<double> mu = eig_raw;  // mu = omega^2 (full complex)

        std::complex<double> dmu_deps = -mu * uMdu / denom;
        std::complex<double> domega_deps = dmu_deps / (2.0 * omega);
        result.df_d_epsilon[i][d] = domega_deps / (2.0 * M_PI);

        std::complex<double> dmu_dmuinv = uKdu / denom;
        std::complex<double> domega_dmuinv = dmu_dmuinv / (2.0 * omega);
        result.df_d_muinv[i][d] = domega_dmuinv / (2.0 * M_PI);

        // No sigma contribution in the linear (lossless) case, but compute it anyway
        // in case someone specifies conductivity without it contributing to C.
        result.df_d_sigma[i][d] = 0.0;
      }
    }

    Mpi::Print(" Mode {:d}: sensitivities computed for {:d} domains\n", i + 1,
               num_domains);
  }

  return result;
}

}  // namespace palace
