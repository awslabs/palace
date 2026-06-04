// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using namespace nlohmann;
using namespace Catch::Matchers;

namespace
{

// Mirror of LoadScaleParMesh2 in test-romoperator.cpp.
auto LoadScaleParMesh(IoData &iodata, MPI_Comm world_comm)
{
  std::vector<std::unique_ptr<Mesh>> mesh_;
  std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
  auto smesh = mesh::Load(iodata, world_comm);
  if (iodata.model.Lc <= 0.0)
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, world_comm);
  }
  iodata.NondimensionalizeInputs(smesh);
  mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), world_comm));
  mesh::RefineMesh(iodata, mfem_mesh);
  for (auto &m : mfem_mesh)
  {
    mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
  }
  return mesh_;
}

}  // namespace

// Verify the factorisation invariant
//   A2(ω) v  ==  Σ_p i·kₙ,p(ω) · M^(p)_{μ⁻¹} v   (as complex ND-vector actions),
// where A2(ω) = `GetExtraSystemMatrix(ω)` (the wave-port contribution) and M^(p)_{μ⁻¹}
// is the per-port boundary mass returned by `GetWavePortBoundaryMassMatrix(p)`. The
// wave-port system stamp is i·kₙ·M with the FULL complex propagation constant
// kₙ = β + iα; for a lossy cross-section (loss-tan / conductivity) the real part −α·M
// lands on the real slot and β·M on the imaginary slot, so the identity must use the
// complex kₙ (GetWavePortKnComplex), not just its real part — the cpw_wave_uniform ports
// sit on a slightly lossy substrate (LossTan ~1e-5) so α ≠ 0. The identity holds at every
// ω since both sides assemble the same bilinear form, the only ω-dependent factor being
// the scalar kₙ,p(ω).
TEST_CASE("WavePortOperator-BoundaryMassFactorisation",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto cpw_source_dir = fs::path(PALACE_SOURCE_EXAMPLES_DIR) / "cpw";
  auto config_path = cpw_source_dir / "cpw_wave_uniform.json";

  // Pull the json from the source tree so the mesh reference (relative path) resolves.
  // Override Mesh to absolute so the test does not depend on cwd.
  std::ifstream f(config_path);
  REQUIRE(f.good());
  json setup = json::parse(f, /*cb=*/nullptr, /*allow_exceptions=*/true,
                           /*ignore_comments=*/true);
  auto mesh_rel = setup["Model"]["Mesh"].get<std::string>();
  setup["Model"]["Mesh"] = (cpw_source_dir / mesh_rel).string();
  // Avoid writing any postprocessing output during unit test.
  setup["Problem"]["Output"] = "";

  IoData iodata(setup, /*print=*/false);
  auto mesh_io = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh_io);

  const auto &wp_op = space_op.GetWavePortOp();
  REQUIRE(wp_op.Size() > 0);

  // Per-port boundary masses (ω-independent).
  std::vector<int> port_idxs;
  std::vector<std::unique_ptr<ComplexOperator>> Mwp_p;
  for (const auto &[idx, data] : wp_op)
  {
    auto Mp =
        space_op.GetWavePortBoundaryMassMatrix<ComplexOperator>(idx, Operator::DIAG_ZERO);
    REQUIRE(Mp);  // CPW waveport boundaries should produce non-empty operators.
    port_idxs.push_back(idx);
    Mwp_p.push_back(std::move(Mp));
  }

  // Drive the random-vector check at three frequencies inside the configured sweep band.
  // Use Nondimensionalize since the SpaceOperator works in internal units.
  std::vector<double> omega_GHz = {3.0, 7.0, 14.0};
  std::vector<double> omega_nd;
  omega_nd.reserve(omega_GHz.size());
  for (double f_GHz : omega_GHz)
  {
    omega_nd.push_back(2.0 * M_PI *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f_GHz));
  }

  // Random ND-vector for action comparison.
  ComplexVector v(Mwp_p.front()->Width());
  v.UseDevice(true);
  linalg::SetRandom(comm, v);

  ComplexVector y_lhs(v.Size()), y_rhs(v.Size()), tmp(v.Size());
  y_lhs.UseDevice(true);
  y_rhs.UseDevice(true);
  tmp.UseDevice(true);

  for (double omega : omega_nd)
  {
    // LHS: the full complex A2(ω) acting on v. A2 encodes the wave-port term as
    // i·kₙ·M = (−α·M) + i·(β·M), so its real slot carries the line-attenuation term and
    // its imaginary slot the propagating term; A2·v carries both.
    auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    REQUIRE(A2);
    y_lhs = 0.0;
    A2->Mult(v, y_lhs);

    // RHS: Σ_p i·kₙ,p(ω) · M^(p)_{μ⁻¹} v, using the FULL complex kₙ = β + iα. The M^(p)
    // we built lives in the imaginary slot, so Mwp_p·v already equals i·(M·v). Scaling
    // that by the complex kₙ gives i·kₙ·(M·v), exactly the per-port system stamp — the
    // attenuation (real-slot −α·M) and propagation (imag-slot β·M) terms both reproduced.
    y_rhs = 0.0;
    for (std::size_t i = 0; i < port_idxs.size(); i++)
    {
      std::complex<double> kn =
          space_op.GetWavePortOp().GetWavePortKnComplex(port_idxs[i], omega);
      tmp = 0.0;
      Mwp_p[i]->Mult(v, tmp);
      // tmp = i·(M·v); y_rhs += kn · tmp = i·kₙ·(M·v).
      linalg::AXPY(kn, tmp, y_rhs);
    }

    // Compare actions. Tolerance accounts for finite-precision matrix assembly &
    // reductions across MPI ranks.
    ComplexVector diff(v.Size());
    diff.UseDevice(true);
    diff = y_lhs;
    linalg::AXPY(std::complex<double>(-1.0, 0.0), y_rhs, diff);
    double rel_err =
        linalg::Norml2(comm, diff) / std::max(linalg::Norml2(comm, y_lhs), 1e-300);

    // Tolerance accounts for floating-point matrix assembly and MPI reductions; the
    // identity is algebraic, but the two sides take slightly different code paths
    // (Σ over per-port operator-vector products vs. one fused boundary form), so the
    // sum-of-products of order N ND-DoFs accumulates O(N·ε) noise. Empirically the
    // residual is ~1e-11 across the cpw_wave_uniform test sweep band.
    CAPTURE(omega, rel_err);
    CHECK(rel_err < 1.0e-10);
  }
}
