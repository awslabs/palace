// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <cmath>
#include <general/forall.hpp>
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"

namespace palace
{

namespace
{

using mfem::ForallWrap;

void InitializeCoefficients(int order, double max_eig_estimate, mfem::Array<double> &coeff)
{
  // Set up Chebyshev coefficients. See mfem::OperatorChebyshevSmoother or Adams et al.,
  // Parallel multigrid smoothing: polynomial versus Gauss-Seidel, JCP (2003).
  coeff.SetSize(order);
  const double upper_bound = 1.1 * max_eig_estimate;
  const double lower_bound = 0.1 * max_eig_estimate;
  const double theta = 0.5 * (upper_bound + lower_bound);
  const double delta = 0.5 * (upper_bound - lower_bound);
  switch (order - 1)
  {
    case 0:
      {
        coeff[0] = 1.0 / theta;
      }
      break;
    case 1:
      {
        double tmp_0 = 1.0 / (std::pow(delta, 2) - 2.0 * std::pow(theta, 2));
        coeff[0] = -4.0 * theta * tmp_0;
        coeff[1] = 2.0 * tmp_0;
      }
      break;
    case 2:
      {
        double tmp_0 = 3.0 * std::pow(delta, 2);
        double tmp_1 = std::pow(theta, 2);
        double tmp_2 = 1.0 / (-4.0 * std::pow(theta, 3) + theta * tmp_0);
        coeff[0] = tmp_2 * (tmp_0 - 12.0 * tmp_1);
        coeff[1] = 12.0 / (tmp_0 - 4.0 * tmp_1);
        coeff[2] = -4.0 * tmp_2;
      }
      break;
    case 3:
      {
        double tmp_0 = std::pow(delta, 2);
        double tmp_1 = std::pow(theta, 2);
        double tmp_2 = 8.0 * tmp_0;
        double tmp_3 =
            1.0 / (std::pow(delta, 4) + 8.0 * std::pow(theta, 4) - tmp_1 * tmp_2);
        coeff[0] = tmp_3 * (32.0 * std::pow(theta, 3) - 16.0 * theta * tmp_0);
        coeff[1] = tmp_3 * (-48.0 * tmp_1 + tmp_2);
        coeff[2] = 32.0 * theta * tmp_3;
        coeff[3] = -8.0 * tmp_3;
      }
      break;
    case 4:
      {
        double tmp_0 = 5.0 * std::pow(delta, 4);
        double tmp_1 = std::pow(theta, 4);
        double tmp_2 = std::pow(theta, 2);
        double tmp_3 = std::pow(delta, 2);
        double tmp_4 = 60.0 * tmp_3;
        double tmp_5 = 20.0 * tmp_3;
        double tmp_6 =
            1.0 / (16 * std::pow(theta, 5) - std::pow(theta, 3) * tmp_5 + theta * tmp_0);
        double tmp_7 = 160.0 * tmp_2;
        double tmp_8 = 1.0 / (tmp_0 + 16.0 * tmp_1 - tmp_2 * tmp_5);
        coeff[0] = tmp_6 * (tmp_0 + 80 * tmp_1 - tmp_2 * tmp_4);
        coeff[1] = tmp_8 * (tmp_4 - tmp_7);
        coeff[2] = tmp_6 * (-tmp_5 + tmp_7);
        coeff[3] = -80.0 * tmp_8;
        coeff[4] = 16.0 * tmp_6;
      }
      break;
    default:
      MFEM_ABORT("Chebyshev smoother is not implemented for order = " << order);
  }
}

class SymmetricScaledOperator : public mfem::Operator
{
private:
  const mfem::Operator &A;
  const mfem::Vector &d;
  mutable mfem::Vector z;

public:
  SymmetricScaledOperator(const mfem::Operator &op, const mfem::Vector &v)
    : mfem::Operator(op.Height()), A(op), d(v), z(v.Size())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A.Mult(x, z);
    {
      const int N = height;
      const auto *D = d.Read();
      const auto *Z = z.Read();
      auto *Y = y.Write();
      MFEM_FORALL(i, N, { Y[i] = D[i] * Z[i]; });
    }
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    {
      const int N = height;
      const auto *D = d.Read();
      const auto *X = x.Read();
      auto *Z = z.Write();
      MFEM_FORALL(i, N, { Z[i] = D[i] * X[i]; });
    }
    A.Mult(z, y);
  }
};

}  // namespace

ChebyshevSmoother::ChebyshevSmoother(MPI_Comm c, const mfem::Array<int> &tdof_list,
                                     int smooth_it, int poly_order)
  : comm(c), A(nullptr), dbc_tdof_list(tdof_list), pc_it(smooth_it), order(poly_order)
{
}

void ChebyshevSmoother::SetOperator(const mfem::Operator &op)
{
  A = &op;
  height = A->Height();
  width = A->Width();
  const int N = height;
  r.SetSize(N);
  z.SetSize(N);
  dinv.SetSize(N);

  // Configure symmetric diagonal scaling.
  mfem::Vector diag(N);
  A->AssembleDiagonal(diag);
  const auto *D = diag.Read();
  auto *DI = dinv.Write();
  MFEM_FORALL(i, N, {
    MFEM_ASSERT_KERNEL(D[i] != 0.0, "Zero diagonal entry in Chebyshev smoother!");
    DI[i] = 1.0 / D[i];
  });
  const auto *I = dbc_tdof_list.Read();
  MFEM_FORALL(i, dbc_tdof_list.Size(), {
    DI[I[i]] = 1.0;  // Assumes operator DiagonalPolicy::ONE
  });

  // Configure coefficients, using the computed maximum eigenvalue estimate.
  petsc::PetscShellMatrix DinvA(comm, std::make_unique<SymmetricScaledOperator>(*A, dinv));
  double max_eig_estimate = DinvA.Norm2();
  InitializeCoefficients(order, max_eig_estimate, coeff);
}

}  // namespace palace
