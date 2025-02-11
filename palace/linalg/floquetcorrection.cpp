// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "floquetcorrection.hpp"

#include <limits>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

template <typename VecType>
FloquetCorrSolver<VecType>::FloquetCorrSolver(const MaterialOperator &mat_op,
                                              FiniteElementSpace &nd_fespace,
                                              FiniteElementSpace &rt_fespace, double tol,
                                              int max_it, int print)
{
  // Create the mass and cross product operators for Floquet correction.
  {
    constexpr bool skip_zeros = false;
    BilinearForm a(rt_fespace);
    a.AddDomainIntegrator<VectorFEMassIntegrator>();
    std::unique_ptr<Operator> m = a.Assemble(skip_zeros);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      M = std::make_unique<ComplexParOperator>(std::move(m), nullptr, rt_fespace);
    }
    else
    {
      M = std::make_unique<ParOperator>(std::move(m), rt_fespace);
    }
  }

  {
    MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute());
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetFloquetCross(), 1.0);
    constexpr bool skip_zeros = false;
    BilinearForm a(nd_fespace, rt_fespace);
    a.AddDomainIntegrator<VectorFEMassIntegrator>(f);
    std::unique_ptr<Operator> m = a.Assemble(skip_zeros);
    if constexpr (std::is_same<OperType, ComplexOperator>::value)
    {
      Cross = std::make_unique<ComplexParOperator>(std::move(m), nullptr, nd_fespace,
                                                   rt_fespace, false);
    }
    else
    {
      Cross = std::make_unique<ParOperator>(std::move(m), nd_fespace, rt_fespace, false);
    }
  }

  // Setup the linear solver.
  auto pcg = std::make_unique<CgSolver<OperType>>(rt_fespace.GetComm(), print);
  pcg->SetInitialGuess(0);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  auto jac = std::make_unique<JacobiSmoother<OperType>>(rt_fespace.GetComm());
  ksp = std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(jac));
  ksp->SetOperators(*M, *M);

  rhs.SetSize(rt_fespace.GetTrueVSize());
  rhs.UseDevice(true);
}

template <typename VecType>
void FloquetCorrSolver<VecType>::Mult(const VecType &x, VecType &y) const
{
  Cross->Mult(x, rhs);
  ksp->Mult(rhs, y);
}

template <typename VecType>
void FloquetCorrSolver<VecType>::AddMult(const VecType &x, VecType &y, ScalarType a) const
{
  this->Mult(x, rhs);
  rhs *= a;
  y += rhs;
}

template class FloquetCorrSolver<ComplexVector>;

}  // namespace palace
