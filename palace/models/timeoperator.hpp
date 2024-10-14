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
// Abstract class handling temporal discretization of the governing equations.
//
class TimeOperator 
{
  protected:
    Vector E, B;
  public:
    TimeOperator() {};
    ~TimeOperator() {};

  // Access solution vectors for E- and B-fields.
  virtual const Vector &GetE() const { return E; }
  virtual const Vector &GetB() const { return B; }

  // Return the linear solver associated with the implicit or explicit time integrator.
  virtual const KspSolver &GetLinearSolver() const = 0;

  // Return if the time integration scheme explicit or implicit.
  virtual bool isExplicit() const = 0;

  // Estimate the maximum stable time step based on the maximum eigenvalue of the
  // undamped system matrix M⁻¹ K.
  virtual double GetMaxTimeStep() const = 0;

  // Initialize time integrators and set 0 initial conditions.
  virtual void Init(double &dt) = 0;

  // Perform time step from t -> t + dt.
  virtual void Step(double &t, double &dt) = 0;

  // Print ODE integrator statistics.
  virtual void PrintStats() {};
};

//
// A class handling the second order ODE form of the governing equations.
//
class SecondOrderTimeOperator : public TimeOperator
{
private:
  // Additional vector storage.
  Vector dE, En;

  // Time integrator for the curl-curl E-field formulation.
  std::unique_ptr<mfem::SecondOrderODESolver> ode;

  // Time-dependent operator for the E-field.
  std::unique_ptr<mfem::SecondOrderTimeDependentOperator> op;

  // Discrete curl for B-field time integration (not owned).
  const Operator *Curl;

public:
  SecondOrderTimeOperator(const IoData &iodata, SpaceOperator &space_op,
                          std::function<double(double)> &dJ_coef);

  // Access E-field time derivative
  const Vector &GetEdot() const { return dE; }

  // Return the linear solver associated with the implicit or explicit time integrator.
  const KspSolver &GetLinearSolver() const;

  // Return if the time integration scheme explicit or implicit.
  bool isExplicit() const { return op->isExplicit(); }

  // Estimate the maximum stable time step based on the maximum eigenvalue of the
  // undamped system matrix M⁻¹ K.
  double GetMaxTimeStep() const;

  // Initialize time integrators and set 0 initial conditions.
  void Init(double &dt);

  // Perform time step from t -> t + dt.
  void Step(double &t, double &dt);

};

//
// A class handling the first order ODE system form of the governing equations.
//
class FirstOrderTimeOperator : public TimeOperator
{
private:
  // Solution vector storage.
  Vector sol;

  // Time integrator for the first order E-B ODE system.
  std::unique_ptr<mfem::ODESolver> ode;

  // Time-dependent operator for the E-B system.
  std::unique_ptr<mfem::TimeDependentOperator> op;

  // Adaptive time-stepping parameters.
  int rk_order;
  bool adapt_dt;
  double rel_tol, abs_tol;

  // Discrete curl for B-field time integration (not owned).
  const Operator *Curl;

public:
  FirstOrderTimeOperator(const IoData &iodata, SpaceOperator &space_op,
                         std::function<double(double)> &J_coef);

  // Return the linear solver associated with the implicit or explicit time integrator.
  const KspSolver &GetLinearSolver() const;

  // Return if the time integration scheme explicit or implicit.
  bool isExplicit() const { return op->isExplicit(); }

  // Estimate the maximum stable time step based on the maximum eigenvalue of the
  // undamped system matrix M⁻¹ K.
  double GetMaxTimeStep() const;

  // Initialize time integrators and set 0 initial conditions.
  void Init(double &dt);

  // Perform time step from t -> t + dt.
  void Step(double &t, double &dt);

  // Print ODE integrator statistics.
  void PrintStats();
};

}  // namespace palace

#endif  // PALACE_MODELS_TIME_OPERATOR_HPP
