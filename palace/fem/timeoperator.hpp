// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_TIME_OPERATOR_HPP
#define PALACE_FEM_TIME_OPERATOR_HPP

#include <functional>
#include <memory>
#include <mfem.hpp>

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
  mfem::Vector E, dE, En, B;

  // Time integrator for the curl-curl E-field formulation.
  std::unique_ptr<mfem::SecondOrderODESolver> ode;

  // Time-dependent operator for the E-field.
  std::unique_ptr<mfem::SecondOrderTimeDependentOperator> op;

  // Discrete curl for B-field time integration.
  std::unique_ptr<mfem::Operator> NegCurl;

public:
  TimeOperator(const IoData &iodata, SpaceOperator &spaceop,
               std::function<double(double)> &djcoef);

  // Access solution vectors for E- and B-fields.
  const mfem::Vector &GetE() const { return E; }
  const mfem::Vector &GetEdot() const { return dE; }
  const mfem::Vector &GetB() const { return B; }

  // Is time integration scheme explicit or implicit.
  bool isExplicit() const { return op->isExplicit(); }

  // Return number of linear solves and linear solver iterations performed during time
  // integration.
  int GetTotalKspMult() const;
  int GetTotalKspIter() const;

  // Estimate the maximum stable time step based on the maximum eigenvalue of the
  // undamped system matrix M⁻¹ K.
  double GetMaxTimeStep() const;

  // Initialize time integrators and set 0 initial conditions.
  void Init();

  // Perform time step from t => t + dt.
  void Step(double &t, double &dt);
};

}  // namespace palace

#endif  // PALACE_FEM_TIME_OPERATOR_HPP
