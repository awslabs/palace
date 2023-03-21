// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorindicators.hpp"

#include "utils/communication.hpp"
#include "utils/quadrules.hpp"

namespace palace
{

void ErrorReductionOperator::operator()(ErrorIndicators &ebar, mfem::Vector &&ind) const
{

  // Compute the global indicator across all processors
  auto comm = Mpi::World();
  double candidate_global_error_indicator = std::accumulate(ind.begin(), ind.end(), 0.0);
  Mpi::GlobalSum(1, &candidate_global_error_indicator, comm);

  ebar.global_error_indicator =
      std::max(ebar.global_error_indicator, candidate_global_error_indicator);

  // update the average local indicator. Using running average update rather
  // than sum and final division to maintain validity at all times.
  auto running_average = [this](const auto &xbar, const auto &x)
  { return (xbar * n + x) / (n + 1); };

  if (n > 0)
  {
    MFEM_VERIFY(ebar.local_error_indicators.Size() == ind.Size(),
                "Local error indicator vectors mismatch.");
    // Combine these error indicators into the current average.
    std::transform(ebar.local_error_indicators.begin(), ebar.local_error_indicators.end(),
                   ind.begin(), ebar.local_error_indicators.begin(), running_average);
  }
  else
  {
    // This is the first sample, just steal the data.
    ebar.local_error_indicators = std::move(ind);
  }

  // Another sample has been added, increment for the running average lambda.
  ++n;
}

// Given a grid function defining a vector solution, compute the error relative
// to a vector coefficient. The default quadrature rule exactly integrates 2p +
// q polynomials, but the quadrature order can be increased or decreased via
// quad_order_increment.
mfem::Vector ComputeElementLpErrors(const mfem::ParGridFunction &sol, double p,
                                    mfem::VectorCoefficient &exsol,
                                    int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();

  mfem::Vector error(fes.GetNE());
  error = 0.0;

  mfem::DenseMatrix vals, exact_vals;
  mfem::Vector loc_errs;

  for (int i = 0; i < fes.GetNE(); i++)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);

    sol.GetVectorValues(T, ir, vals);
    exsol.Eval(exact_vals, T, ir);

    vals -= exact_vals;
    loc_errs.SetSize(vals.Width());

    // compute the lengths of the errors at the integration points thus the
    // vector norm is rotationally invariant
    vals.Norm2(loc_errs);

    for (int j = 0; j < ir.GetNPoints(); j++)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);
      double errj = loc_errs(j);
      if (p < mfem::infinity())
      {
        errj = pow(errj, p);

        error[i] += ip.weight * T.Weight() * errj;
      }
      else
      {
        error[i] = std::max(error[i], errj);
      }
    }
    if (p < mfem::infinity())
    {
      // negative quadrature weights may cause the error to be negative
      if (error[i] < 0.)
      {
        error[i] = -std::pow(-error[i], 1. / p);
      }
      else
      {
        error[i] = std::pow(error[i], 1. / p);
      }
    }
  }
  return error;
}

}  // namespace palace
