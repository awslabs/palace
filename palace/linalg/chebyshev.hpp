// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_CHEBYSHEV_SMOOTHER_HPP
#define PALACE_CHEBYSHEV_SMOOTHER_HPP

#include <mfem.hpp>

namespace palace
{

//
// Matrix-free diagonally-scaled Chebyshev smoothing. This is largely the same as
// mfem::OperatorChebyshevSmoother allows a nonzero initial guess and uses alternative
// methods to estimate the largest eigenvalue.
//
class ChebyshevSmoother : public mfem::Solver
{
private:
  // System matrix (not owned), its communicator, and list of eliminated degrees of freedom.
  MPI_Comm comm;
  const mfem::Operator *A;
  const mfem::Array<int> dbc_tdof_list;

  // Number of smoother iterations and polynomial order.
  const int pc_it, order;

  // Diagonal scaling of the operator.
  mfem::Vector dinv;

  // Array of coefficients for Chebyshev polynomial smoothing.
  mfem::Array<double> coeff;

  // Temporary vectors for smoother application.
  mutable mfem::Vector r, z;

public:
  ChebyshevSmoother(MPI_Comm c, const mfem::Array<int> &tdof_list, int smooth_it,
                    int poly_order);

  void SetOperator(const mfem::Operator &op) override;

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    // y = y + p(A) (x - A y)
    for (int it = 0; it < pc_it; it++)
    {
      if (iterative_mode || it > 0)
      {
        A->Mult(y, r);
        subtract(x, r, r);
      }
      else
      {
        r = x;
        y = 0.0;
      }
      r *= dinv;
      y.Add(coeff[0], r);
      for (int k = 1; k < order; k++)
      {
        z = r;
        A->Mult(z, r);
        r *= dinv;
        y.Add(coeff[k], r);
      }
    }
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    Mult(x, y);  // Assumes operator symmetry
  }
};

}  // namespace palace

#endif  // PALACE_CHEBYSHEV_SMOOTHER_HPP
