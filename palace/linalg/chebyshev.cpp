// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <vector>
#include <mfem/general/forall.hpp>
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"

namespace palace
{

namespace
{

using mfem::ForallWrap;

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

  // Configure symmetric diagonal scaling.
  const int N = height;
  dinv.SetSize(N);
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

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al., Parallel multigrid smoothing:
  // polynomial versus Gauss-Seidel, JCP (2003).
  petsc::PetscShellMatrix DinvA(comm, std::make_unique<SymmetricScaledOperator>(*A, dinv));
  lambda_max = 1.1 * DinvA.Norm2();
}

void ChebyshevSmoother::ArrayMult(const mfem::Array<const mfem::Vector *> &X,
                                  mfem::Array<mfem::Vector *> &Y) const
{
  // Initialize.
  const int nrhs = X.Size();
  mfem::Array<mfem::Vector *> R(nrhs), D(nrhs);
  std::vector<mfem::Vector> rrefs(nrhs), drefs(nrhs);
  if (nrhs * height != r.Size())
  {
    r.SetSize(nrhs * height);
    d.SetSize(nrhs * height);
  }
  for (int j = 0; j < nrhs; j++)
  {
    rrefs[j].MakeRef(r, j * height, height);
    drefs[j].MakeRef(d, j * height, height);
    R[j] = &rrefs[j];
    D[j] = &drefs[j];
  }

  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (iterative_mode || it > 0)
    {
      A->ArrayMult(Y, R);
      for (int j = 0; j < nrhs; j++)
      {
        subtract(*X[j], *R[j], *R[j]);
      }
    }
    else
    {
      for (int j = 0; j < nrhs; j++)
      {
        *R[j] = *X[j];
        *Y[j] = 0.0;
      }
    }

    // 4th-kind Chebyshev smoother
    {
      const auto *DI = dinv.Read();
      for (int j = 0; j < nrhs; j++)
      {
        const auto *RR = R[j]->Read();
        auto *DD = D[j]->ReadWrite();
        MFEM_FORALL(i, height, { DD[i] = 4.0 / (3.0 * lambda_max) * DI[i] * RR[i]; });
      }
    }
    for (int k = 1; k < order; k++)
    {
      A->ArrayAddMult(D, R, -1.0);
      {
        // From Phillips and Fischer or Lottes (with k -> k + 1 shift due to 1-based
        // indexing)
        const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
        const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
        const auto *DI = dinv.Read();
        for (int j = 0; j < nrhs; j++)
        {
          const auto *RR = R[j]->Read();
          auto *YY = Y[j]->ReadWrite();
          auto *DD = D[j]->ReadWrite();
          MFEM_FORALL(i, height, {
            YY[i] += DD[i];
            DD[i] = sd * DD[i] + sr * DI[i] * RR[i];
          });
        }
      }
    }
    for (int j = 0; j < nrhs; j++)
    {
      *Y[j] += *D[j];
    }
  }
}

}  // namespace palace
