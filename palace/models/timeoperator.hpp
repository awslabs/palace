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
  Vector E, dE, En, B;

  // Time integrator for the curl-curl E-field formulation.
  std::unique_ptr<mfem::SecondOrderODESolver> ode;

  // Time-dependent operator for the E-field.
  std::unique_ptr<mfem::SecondOrderTimeDependentOperator> op;

  // Discrete curl for B-field time integration (not owned).
  const Operator *Curl;

public:
  TimeOperator(const IoData &iodata, SpaceOperator &space_op,
               std::function<double(double)> &dJ_coef);

  // Access solution vectors for E- and B-fields.
  const Vector &GetE() const { return E; }
  const Vector &GetEdot() const { return dE; }
  const Vector &GetB() const { return B; }

  // Return the linear solver associated with the implicit or explicit time integrator.
  const KspSolver &GetLinearSolver() const;

  // Return if the time integration scheme explicit or implicit.
  bool isExplicit() const { return op->isExplicit(); }

  // Estimate the maximum stable time step based on the maximum eigenvalue of the
  // undamped system matrix M⁻¹ K.
  double GetMaxTimeStep() const;

  // Initialize time integrators and set 0 initial conditions.
  void Init();

  // Perform time step from t -> t + dt.
  void Step(double &t, double &dt);
};

}  // namespace palace

#endif  // PALACE_MODELS_TIME_OPERATOR_HPP
