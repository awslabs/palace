// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "utils/configfile.hpp"

mfem::Mesh SingleTetMesh();

namespace palace::test
{

// Restores the global assembly and integration-order settings on destruction. With a
// SolverData argument, first installs its integration orders, as the drivers do.
struct IntegrationSettingsGuard
{
  int pa_order_threshold = BilinearForm::pa_order_threshold;
  int p_trial = fem::DefaultIntegrationOrder::p_trial;
  bool q_order_jac = fem::DefaultIntegrationOrder::q_order_jac;
  int q_order_extra_pk = fem::DefaultIntegrationOrder::q_order_extra_pk;
  int q_order_extra_qk = fem::DefaultIntegrationOrder::q_order_extra_qk;

  IntegrationSettingsGuard() = default;
  explicit IntegrationSettingsGuard(const config::SolverData &solver)
  {
    fem::DefaultIntegrationOrder::p_trial = solver.order;
    fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
    fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;
  }
  IntegrationSettingsGuard(const IntegrationSettingsGuard &) = delete;
  IntegrationSettingsGuard &operator=(const IntegrationSettingsGuard &) = delete;
  ~IntegrationSettingsGuard()
  {
    BilinearForm::pa_order_threshold = pa_order_threshold;
    fem::DefaultIntegrationOrder::p_trial = p_trial;
    fem::DefaultIntegrationOrder::q_order_jac = q_order_jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = q_order_extra_pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = q_order_extra_qk;
  }
};

}  // namespace palace::test

#endif  // PALACE_TEST_HELPERS_HPP
