// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <array>
#include <complex>
#include <cstdlib>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "utils/configfile.hpp"

mfem::Mesh SingleTetMesh();

// A box [0, 2] x [0, 1]^2 of 4 x 2 x 2 hexahedra: the domains x < 1 (attribute 1) and x > 1
// (2), the faces x = 0 (1), x = 2 (2) and the others (3), and the interior boundaries x = 1
// between the domains (4) and x = 3/2 inside the second (5).
mfem::Mesh MakeInteriorBoundaryBox();

// What the calling process prints to stdout while f runs (rethrowing what f throws).
std::string CaptureStdout(const std::function<void()> &f);

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

namespace palace
{
class DrivenSubstructureModel;
}  // namespace palace

namespace palace::test
{

// The standard output of this process while f runs (on the root process, the output of
// Mpi::Print and Mpi::Warning). The standard output is restored also when f throws.
std::string CaptureStdout(const std::function<void()> &f);

// Set a substructuring configuration for an offline run which saves its model at
// model_path (environment.model in the output directory, which it creates), or for an
// online run ("Online") which loads it.
void SetModelRun(nlohmann::json &config, std::string_view mode,
                 const std::string &model_path);

// The record j of the frequencies in the file of a driven substructuring model, as written:
// the lower triangle of S by columns, then g, h and c (see
// DrivenSubstructureModel::PartSizes), read by the calling process.
std::vector<std::complex<double>> ReadDrivenRecord(const DrivenSubstructureModel &model,
                                                   const std::string &path, int j);

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
