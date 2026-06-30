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
//   Im{A2(ω) v}  ==  Σ_p k_{n,p}(ω) · M_{μ⁻¹,p} v   (as ND-vector actions),
// where A2(ω) is the imaginary part of `GetExtraSystemMatrix(ω)` (the wave-port
// contribution; cf. waveportoperator.cpp:1080-1088) and M_{μ⁻¹,p} is the new per-port
// boundary mass returned by `GetWavePortBoundaryMassMatrix(p)`. The identity must hold
// at every ω since both sides assemble the same bilinear form, with the only ω-dependent
// factor being the scalar k_{n,p}(ω).
TEST_CASE("WavePortOperator-BoundaryMassFactorisation",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  // The cpw example config and mesh are installed under the test data directory via the
  // regression fixtures (which symlink to examples/cpw and are dereferenced on install).
  auto cpw_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "cpw";
  auto config_path = cpw_dir / "cpw_wave_uniform.json";

  // Override Mesh to absolute so the relative mesh reference resolves regardless of cwd.
  std::ifstream f(config_path);
  REQUIRE(f.good());
  json setup = json::parse(f, /*cb=*/nullptr, /*allow_exceptions=*/true,
                           /*ignore_comments=*/true);
  auto mesh_rel = setup["Model"]["Mesh"].get<std::string>();
  setup["Model"]["Mesh"] = (cpw_dir / mesh_rel).string();
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
    // LHS: imaginary part of A2(ω) acting on v. A2 currently encodes the wave-port term
    // as fbi (imaginary part of the bilinear form), so Im{A2 v} carries the k_n-scaled
    // boundary mass action.
    auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    REQUIRE(A2);
    y_lhs = 0.0;
    A2->Mult(v, y_lhs);

    // RHS: Σ_p k_{n,p}(ω) · M_{μ⁻¹,p} v. The M_p we built lives in the imaginary part,
    // so we want Im{Mwp_p · v} multiplied by the real scalar k_n.
    y_rhs = 0.0;
    for (std::size_t i = 0; i < port_idxs.size(); i++)
    {
      double kn = space_op.GetWavePortOp().GetWavePortKn(port_idxs[i], omega);
      tmp = 0.0;
      Mwp_p[i]->Mult(v, tmp);
      // tmp is a ComplexOperator action on a complex v; Mwp_p is purely imaginary, so its
      // action multiplies the input by i times the real boundary-mass action. Adding the
      // ports together and scaling by kn gives the Σ form above.
      // y_rhs += kn * tmp
      linalg::AXPY(std::complex<double>(kn, 0.0), tmp, y_rhs);
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
