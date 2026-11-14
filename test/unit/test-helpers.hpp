// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <array>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "utils/configfile.hpp"

mfem::Mesh SingleTetMesh();

// Path to the unit test mesh with the given file name. Most test meshes are stored in
// PALACE_TEST_DATA_DIR/mesh. "star-tri.mesh" and "fichera-tet.mesh" are not stored: they
// are "star-quad.mesh" and "fichera-hex.mesh" after mfem::Mesh::MakeSimplicial, written on
// first use to a temporary file private to the calling process and removed at exit.
std::string TestMeshPath(const std::string &name);

// Sets an environment variable of the calling process, or unsets it when the value is
// null, and restores its previous state on destruction.
class ScopedEnvironmentVariable
{
private:
  std::string name;
  std::optional<std::string> previous;

  void Set(const char *value)
  {
    if (value)
    {
      setenv(name.c_str(), value, 1);
    }
    else
    {
      unsetenv(name.c_str());
    }
  }

public:
  ScopedEnvironmentVariable(const std::string &name, const char *value) : name(name)
  {
    if (const char *old = std::getenv(name.c_str()))
    {
      previous = old;
    }
    Set(value);
  }
  ScopedEnvironmentVariable(const ScopedEnvironmentVariable &) = delete;
  ScopedEnvironmentVariable &operator=(const ScopedEnvironmentVariable &) = delete;
  ~ScopedEnvironmentVariable() { Set(previous ? previous->c_str() : nullptr); }
};

namespace palace::test
{

// Deterministic pseudo-random orthonormal frames (e1, e2, normal) of planes oblique to
// every coordinate axis.
std::vector<std::array<mfem::Vector, 3>> ObliqueFrames(int num_frames);

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
