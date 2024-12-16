// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_TIME_OPERATOR_HPP
#define PALACE_MODELS_TIME_OPERATOR_HPP

#include <functional>
#include <memory>
#include <mfem.hpp>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

//
// A class handling temporal discretization of the governing equations.
//
class TimeOperator
{
private:
  // Solution vector storage.
  Vector E, B, sol;

  // Time integrator for the first order ODE system.
  std::unique_ptr<mfem::ODESolver> ode;

  // Time-dependent operator for the Edot-E ODE system.
  std::unique_ptr<mfem::TimeDependentOperator> op;

  // Adaptive time-stepping parameters.
  int order;
  mfem::real_t rel_tol, abs_tol;
  bool use_mfem_integrator = false;

public:
  TimeOperator(const IoData &iodata, SpaceOperator &space_op,
               std::function<mfem::real_t(mfem::real_t)> dJ_coef);

  // Access solution vectors for E- and B-fields.
  const Vector &GetE() const { return E; }
  const Vector &GetB() const { return B; }

  // Return the linear solver associated with the implicit or explicit time integrator.
  const KspSolver &GetLinearSolver() const;

  // Initialize time integrators and set 0 initial conditions.
  void Init();

  // Perform time step from t -> t + dt.
  void Step(mfem::real_t &t, mfem::real_t &dt);

  // Print ODE integrator statistics.
  void PrintStats();
};

}  // namespace palace

#endif  // PALACE_MODELS_TIME_OPERATOR_HPP
