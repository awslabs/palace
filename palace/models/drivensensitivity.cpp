// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensensitivity.hpp"

#include <complex>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Compute the complex sesquilinear form v^H A u for a ComplexOperator A.
// linalg::Dot(x, y) computes y^H x, so Dot(Au, v) = v^H (Au) = v^H A u.
std::complex<double> SesquilinearForm(MPI_Comm comm, const ComplexOperator &A,
                                      const ComplexVector &v, const ComplexVector &u,
                                      ComplexVector &Au)
{
  A.Mult(u, Au);
  return linalg::Dot(comm, Au, v);
}

// Wrapper operator that redirects Mult to MultHermitianTranspose of the wrapped operator.
// This allows reusing a standard KSP solver for the adjoint system A^H v = g.
class HermitianTransposeOperator : public ComplexOperator
{
private:
  const ComplexOperator &A;

public:
  HermitianTransposeOperator(const ComplexOperator &op)
    : ComplexOperator(op.Width(), op.Height()), A(op)
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    A.MultHermitianTranspose(x, y);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    A.AddMultHermitianTranspose(x, y, a);
  }
};

}  // namespace

DrivenSensitivity::SensitivityResult
DrivenSensitivity::Compute(const IoData &iodata, SpaceOperator &space_op,
                           ComplexKspSolver &ksp, const ComplexOperator &A,
                           const ComplexOperator &P, const ComplexVector &E,
                           int excitation_idx, double omega)
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

  // Collect port indices for which we compute sensitivities.
  for (const auto &[idx, data] : space_op.GetLumpedPortOp())
  {
    result.port_indices.push_back(idx);
  }
  const int num_ports = static_cast<int>(result.port_indices.size());

  if (num_ports <= 0 || num_domains <= 0)
  {
    return result;
  }

  // Allocate result arrays: [num_ports][num_domains].
  result.dS_d_epsilon.resize(num_ports,
                             std::vector<std::complex<double>>(num_domains, 0.0));
  result.dS_d_muinv.resize(num_ports, std::vector<std::complex<double>>(num_domains, 0.0));
  result.dS_d_sigma.resize(num_ports, std::vector<std::complex<double>>(num_domains, 0.0));

  // Working vectors sized to the local ND true dof count.
  const int n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector adjoint_rhs(n), adjoint_sol(n), Au(n);
  adjoint_rhs.UseDevice(true);
  adjoint_sol.UseDevice(true);
  Au.UseDevice(true);

  Mpi::Print("\nComputing S-parameter sensitivities with respect to material properties\n");

  // Set up the adjoint operator A^H for the KSP solver. The preconditioner P is reused
  // from the forward solve without transposition: since it is only an approximation, and
  // A and A^H share the same spectral structure, the forward preconditioner is effective
  // for the adjoint system. This also avoids issues with preconditioner implementations
  // (e.g., direct solvers, AMG) that require assembled matrices rather than operator
  // wrappers.
  HermitianTransposeOperator AH(A);
  ksp.SetOperators(AH, P);

  // Pre-assemble domain-restricted operators outside the port loop (they depend only on
  // the domain, not the port).
  std::vector<std::unique_ptr<ComplexOperator>> dK_vec(num_domains);
  std::vector<std::unique_ptr<ComplexOperator>> dM_vec(num_domains);
  std::vector<std::unique_ptr<ComplexOperator>> dC_vec(num_domains);
  for (int d = 0; d < num_domains; d++)
  {
    const int attr = result.domain_attrs[d];
    dK_vec[d] =
        space_op.GetStiffnessMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
    dM_vec[d] = space_op.GetMassMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
    dC_vec[d] = space_op.GetDampingMatrixDomain<ComplexOperator>(attr, Operator::DIAG_ZERO);
  }

  // For each observation port j, solve the adjoint system and accumulate sensitivities.
  for (int j = 0; j < num_ports; j++)
  {
    const int port_idx = result.port_indices[j];

    // Assemble adjoint RHS: g = s_j, the S-parameter linear form for port j as a
    // true-DOF vector. The S-parameter is S_ij = s_j^T u_i where s_j is a real linear
    // form, so the adjoint RHS is simply s_j (real part only, imaginary part = 0).
    space_op.GetSParameterAdjointRHS(port_idx, adjoint_rhs);

    // Solve adjoint system: A^H v_j = s_j.
    adjoint_sol = 0.0;
    Mpi::Print(" Adjoint solve for observation port {:d}\n", port_idx);
    ksp.Mult(adjoint_rhs, adjoint_sol);

    // Accumulate sensitivities for each material domain.
    //   dS_ij/dp = -v_j^H (dA/dp) u_i
    for (int d = 0; d < num_domains; d++)
    {
      // dA/d(eps_d) = -omega^2 * dM_d
      if (dM_vec[d])
      {
        auto vMu = SesquilinearForm(comm, *dM_vec[d], adjoint_sol, E, Au);
        result.dS_d_epsilon[j][d] = -(-omega * omega) * vMu;
      }

      // dA/d(mu_inv_d) = dK_d
      if (dK_vec[d])
      {
        auto vKu = SesquilinearForm(comm, *dK_vec[d], adjoint_sol, E, Au);
        result.dS_d_muinv[j][d] = -vKu;
      }

      // dA/d(sigma_d) = i*omega * dC_d
      if (dC_vec[d])
      {
        auto vCu = SesquilinearForm(comm, *dC_vec[d], adjoint_sol, E, Au);
        result.dS_d_sigma[j][d] = -(1i * omega) * vCu;
      }
    }

    Mpi::Print(" Port {:d}: sensitivities computed for {:d} domains\n", port_idx,
               num_domains);
  }

  // Restore the original (forward) operators for future solves.
  ksp.SetOperators(A, P);

  return result;
}

}  // namespace palace
