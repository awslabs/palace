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

// The standard output of this process while f runs (on the root process, the output of
// Mpi::Print and Mpi::Warning). The standard output is restored also when f throws.
std::string CaptureStdout(const std::function<void()> &f);

// Deterministic pseudo-random orthonormal frames (e1, e2, normal) of planes oblique to
// every coordinate axis.
std::vector<std::array<mfem::Vector, 3>> ObliqueFrames(int num_frames);

// Propagation constant of a mode of the 2D cross-section 0 < x < W, 0 < y < H of a guide in
// vacuum between PEC plates at y = 0 and y = H, with absorbing boundaries of the given
// order at x = 0 and x = W, at the free-space wavenumber k. The fields E = ∇ × (x̂ ψ) with
// ψ = cos(κ(x - W/2)) cos(mπy/H) e^{-ikₙz} and κ² = k² - kₙ² - (mπ/H)² have Eₓ = 0 and, for
// m > 0, an out-of-plane component Eₙ which varies along the absorbing boundaries. On
// these, the field has the tangential wavenumber q, with q² = kₙ² + (mπ/H)² = k² - κ², and
// satisfies the absorbing boundary condition for ∂ₓψ = ∓iγψ at x = W and x = 0: the
// first-order condition for γ = k, and the second-order one, whose term acts on the normal
// component -iωμHₓ = q²ψ of the curl of E, for γ = k - q²/(2k) (the normal wavenumber of
// the condition, see docs/src/reference.md). So κ tan(κW/2) = iγ. For m > 0, the fields E =
// ∇ × ∇ × (x̂ φ) with φ = cos(κ(x - W/2)) sin(mπy/H) e^{-ikₙz} have Hₓ = 0, so that the
// second-order term vanishes on them, and satisfy both conditions for γ = k: they are the
// modes of the first-order condition.
std::complex<double> ExactAbsorbingParallelPlateKn(double W, double H, double k, int m,
                                                   int abc_order);

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
