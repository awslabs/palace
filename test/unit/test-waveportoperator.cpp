// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <complex>
#include <fstream>
#include <memory>
#include <numbers>
#include <vector>
#include <Eigen/Dense>
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

json LoadCpwWaveConfig()
{
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
  // Avoid writing any postprocessing output during unit tests.
  setup["Problem"]["Output"] = "";
  return setup;
}

// Dielectric-slab-loaded rectangular guide (see the iris_filter regression fixture). Its
// hybrid/LSM port mode has a frequency-rotating transverse shape, so the modal correction
// is genuinely multi-vector -- the case exercised by ModalCorrectionRotationSubspace below.
json LoadIrisWaveConfig()
{
  auto iris_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "iris_filter";
  std::ifstream f(iris_dir / "driven_wave_synth.json");
  REQUIRE(f.good());
  json setup = json::parse(f, /*cb=*/nullptr, /*allow_exceptions=*/true,
                           /*ignore_comments=*/true);
  auto mesh_rel = setup["Model"]["Mesh"].get<std::string>();
  setup["Model"]["Mesh"] = (iris_dir / mesh_rel).string();
  setup["Problem"]["Output"] = "";
  return setup;
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
  IoData iodata(LoadCpwWaveConfig(), /*print=*/false);
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
    omega_nd.push_back(2.0 * std::numbers::pi *
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

// Verify the applied wave-port operator (the sparse local boundary mass i·k_n·M plus the
// modal correction W = Σ_p (W_full − W_scalar)) is matched to its own mode. Applied to the
// parent-space modal field e_p of the excited port, it must reproduce the injected boundary
// term −iω·s_full,p, which is half the driven excitation vector (RHS₂ = −2iω·s_full,p, cf.
// AddExcitationBdrCoefficients): the enforced operator and the excitation share the same
// full modal n×H response, so the port is reflectionless for its own mode. On the mode the
// mass action i·k_n·M·e = −iω·s_scalar,p cancels the −W_scalar term, leaving the full n×H.
TEST_CASE("WavePortOperator-ModalCorrectionMatchedMode",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(LoadCpwWaveConfig(), /*print=*/false);
  auto mesh_io = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh_io);

  auto &wp_op = space_op.GetWavePortOp();
  REQUIRE(wp_op.Size() > 0);

  const double omega = 2.0 * std::numbers::pi *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(7.0);

  // The full applied wave-port operator: sparse local mass i·k_n·M plus the modal
  // correction (no Floquet ports here). This also triggers the per-port mode/reaction solve
  // at omega.
  auto W = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  REQUIRE(W);

  auto &nd_fespace = space_op.GetNDSpace();
  const auto &mesh = *nd_fespace.Get().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const int n = nd_fespace.GetTrueVSize();

  // Locate the excited port (Excitation == 1 in the cpw config).
  int exc_idx = 0;
  for (const auto &[idx, data] : wp_op)
  {
    if (data.HasExcitation())
    {
      exc_idx = data.excitation;
      break;
    }
  }
  REQUIRE(exc_idx != 0);

  // Reconstruct the excited port's modal field on the parent ND space (its tangential trace
  // over the port boundary; zero elsewhere), from the same mode-field coefficient used for
  // S-projection. sᵀe only samples e on the port boundary, so essential/interior dofs of e
  // do not affect W·e.
  ComplexVector e(n);
  e.UseDevice(true);
  e = 0.0;
  for (const auto &[idx, data] : wp_op)
  {
    if (data.excitation != exc_idx)
    {
      continue;
    }
    mfem::Array<int> marker = mesh::AttrToMarker(bdr_attr_max, data.GetAttrList());
    auto er = data.GetModeFieldCoefficientReal();
    auto ei = data.GetModeFieldCoefficientImag();
    mfem::ParGridFunction e_gf_r(&nd_fespace.Get()), e_gf_i(&nd_fespace.Get());
    e_gf_r = 0.0;
    e_gf_i = 0.0;
    e_gf_r.ProjectBdrCoefficientTangent(*er, marker);
    e_gf_i.ProjectBdrCoefficientTangent(*ei, marker);
    e_gf_r.GetTrueDofs(e.Real());
    e_gf_i.GetTrueDofs(e.Imag());
    break;
  }

  // (i·k_n·M + W)·e should equal −iω·s_full = RHS₂/2 for the excited port.
  ComplexVector We(n), rhs(n), expected(n), diff(n);
  We.UseDevice(true);
  rhs.UseDevice(true);
  expected.UseDevice(true);
  diff.UseDevice(true);
  We = 0.0;
  W->Mult(e, We);

  rhs = 0.0;
  bool nnz = space_op.GetExcitationVector2(exc_idx, omega, rhs);
  REQUIRE(nnz);
  expected = 0.0;
  linalg::AXPY(std::complex<double>(0.5, 0.0), rhs, expected);

  diff = We;
  linalg::AXPY(std::complex<double>(-1.0, 0.0), expected, diff);
  double rel_err =
      linalg::Norml2(comm, diff) / std::max(linalg::Norml2(comm, expected), 1e-300);
  CAPTURE(omega, rel_err);
  // The identity is exact in the FE space up to (i) interpolation of the mode-field
  // coefficient (ProjectBdrCoefficientTangent) versus the linear-form assembly of s, (ii)
  // the submesh-vs-parent reaction consistency (R vs sᵀe), (iii) the cancellation of the
  // sparse local-mass modal action against the −W_scalar term (assembled via different code
  // paths), and (iv) MPI reductions.
  CHECK(rel_err < 1.0e-8);
}

// The complex-ω (eigenmode) correction reconstructs n×H with full complex kₙ; the real-ω
// (driven) correction uses Re{kₙ} and lets the sparse mass carry the Im{kₙ} attenuation. So
// at ω=ω0 they differ by the mode attenuation a≈|Im{kₙ0}|/|Re{kₙ0}|. Check the relative
// action difference is ~a: present (guards the eigenmode path from reverting to Re{kₙ}) and
// bounded (paths track, vanishing as a→0).
TEST_CASE("WavePortOperator-ModalCorrectionComplexAttenuation",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(LoadCpwWaveConfig(), /*print=*/false);
  auto mesh_io = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh_io);

  auto &wp_op = space_op.GetWavePortOp();
  REQUIRE(wp_op.Size() > 0);
  auto &nd_fespace = space_op.GetNDSpace();
  const int n = nd_fespace.GetTrueVSize();
  const double omega = 2.0 * std::numbers::pi *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(7.0);

  // Real-ω operator (triggers Initialize(ω0=omega)), then the complex-ω operator reusing
  // the frozen reference; same (empty) essential-dof list so any difference is purely the
  // kₙ convention.
  mfem::Array<int> dbc;
  auto W_real = wp_op.GetModalCorrectionOperator(omega, nd_fespace, dbc);
  auto W_cplx =
      wp_op.GetModalCorrectionOperator(std::complex<double>(omega, 0.0), nd_fespace, dbc);
  REQUIRE(W_real);
  REQUIRE(W_cplx);

  // Largest per-port attenuation ratio a = |Im{kₙ0}|/|Re{kₙ0}| (kn0 frozen by Initialize).
  double atten = 0.0;
  for (const auto &[idx, data] : wp_op)
  {
    if (data.active && std::abs(data.kn0.real()) > 0.0)
    {
      atten = std::max(atten, std::abs(data.kn0.imag()) / std::abs(data.kn0.real()));
    }
  }
  REQUIRE(atten > 0.0);  // cpw mode is lossy; the two paths must differ.

  // Apply both to a fixed nonzero vector and compare (each is Σ_k g_k (s_kᵀx) s_k).
  ComplexVector x(n), yr(n), yc(n), diff(n);
  for (auto *v : {&x, &yr, &yc, &diff})
  {
    v->UseDevice(true);
  }
  x = 1.0;
  x.Real().Randomize(1);
  x.Imag().Randomize(2);
  yr = 0.0;
  yc = 0.0;
  W_real->Mult(x, yr);
  W_cplx->Mult(x, yc);
  diff = yc;
  linalg::AXPY(std::complex<double>(-1.0, 0.0), yr, diff);
  double rel_err = linalg::Norml2(comm, diff) / std::max(linalg::Norml2(comm, yr), 1e-300);
  CAPTURE(omega, atten, rel_err);
  // The attenuation term largely cancels between W_full and W_scalar, so the net difference
  // is ~0.04a, well above the EVP re-solve floor. Bracket loosely: nonzero guards the
  // full-kₙ term, the upper bound guards against gross divergence.
  CHECK(rel_err > 0.005 * atten);
  CHECK(rel_err < 0.5 * atten);
}

// An inactive wave port is an unloaded boundary, but its unit mass remains part of the
// synthesis coordinate definition when IncludeInSynthesis is true. Verify that exposing
// this unit operator does not accidentally reactivate the physical Robin termination.
TEST_CASE("WavePortOperator-InactiveBoundaryMassForSynthesis",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto setup = LoadCpwWaveConfig();
  // Exercise the per-port doubled-real preconditioner override on the restricted port MPI
  // communicator while leaving the full 3D linear solver on its default real factorization.
  // The global real-valued preconditioner approximation must not silently disable the
  // explicit per-port choice.
  setup["Solver"]["Linear"]["ComplexCoarseSolve"] = false;
  setup["Solver"]["Linear"]["PCMatReal"] = true;
  bool first_port = true;
  for (auto &port : setup["Boundaries"]["WavePort"])
  {
    port["ComplexCoarseSolve"] = true;
    port["Active"] = false;
    port["IncludeInSynthesis"] = true;
    port["Excitation"] = first_port ? 1 : 0;
    first_port = false;
  }

  IoData iodata(setup, /*print=*/false);
  auto mesh_io = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh_io);

  const auto &wp_op = space_op.GetWavePortOp();
  REQUIRE(wp_op.Size() > 0);
  for (const auto &[idx, data] : wp_op)
  {
    CHECK_FALSE(data.active);
    CHECK(data.include_in_synthesis);
    auto Mp =
        space_op.GetWavePortBoundaryMassMatrix<ComplexOperator>(idx, Operator::DIAG_ZERO);
    REQUIRE(Mp);
  }

  const double omega = 2.0 * std::numbers::pi *
                       iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(7.0);
  auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  CHECK_FALSE(A2);
  auto A2_complex = space_op.GetExtraSystemMatrix(std::complex<double>(omega, 1.0e-3),
                                                  Operator::DIAG_ZERO);
  CHECK_FALSE(A2_complex);
}

// The circuit-synthesis modal correction W(ω) = Σ_p g_full,p s_full,p s_full,pᵀ −
// g_scalar,p s_scalar,p s_scalar,pᵀ is folded into a fixed reduced pencil by tracking the
// port modal vectors s(ω) across the band. A homogeneous cross-section has a separable mode
// whose transverse shape is fixed (only k_n(ω) scales), so the sampled s(ω) span a rank-1
// subspace and a center-frozen shape is exact. A transversely inhomogeneous guide (slab)
// supports hybrid modes whose shape rotates with ω, so the samples span rank>=2 and a
// center-frozen single vector is inadequate. Verify both on the slab-loaded iris guide:
// (a) the s_full samples are genuinely rank>=2, (b) projecting onto the band-center shape
// leaves a large residual.
TEST_CASE("WavePortOperator-ModalCorrectionRotationSubspace",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  IoData iodata(LoadIrisWaveConfig(), /*print=*/false);
  auto mesh_io = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh_io);

  auto nd = [&](double f)
  { return iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(f); };
  // Sample from just above the mode cutoff (~6.3 GHz), where the hybrid shape redistributes
  // fastest between slab and air, up through well-separated, to expose the rotation.
  const double f_lo = 6.5, f_hi = 13.0;
  const double w_ref = 2.0 * std::numbers::pi * nd(0.5 * (f_lo + f_hi));

  auto ports = space_op.GetModalCorrectionSynthesisPorts(w_ref);
  REQUIRE(!ports.empty());

  const int msamp = 9;
  bool checked_any = false;
  for (int port_idx : ports)
  {
    // Sample the full modal vector across the band (the band-center sample is the "frozen"
    // reference shape, picked as i_ref below).
    std::vector<std::unique_ptr<ComplexVector>> sf;
    std::vector<double> nf;
    for (int i = 0; i < msamp; i++)
    {
      const double fi = f_lo + (f_hi - f_lo) * i / (msamp - 1);
      auto smp = space_op.SampleModalCorrectionVectors(
          port_idx, std::complex<double>(2.0 * std::numbers::pi * nd(fi), 0.0));
      if (!smp.active)
      {
        continue;
      }
      nf.push_back(linalg::Norml2(comm, *smp.s_full));
      sf.push_back(std::move(smp.s_full));
    }
    const int m = static_cast<int>(sf.size());
    if (m < 2)
    {
      continue;  // inactive/cutoff port modes carry no correction to track
    }
    checked_any = true;

    // (a) Rank of the sampled subspace: the Hermitian Gram of the unit-normalized samples
    // has eigenvalues equal to the squared singular values of the stacked matrix. A
    // rotating mode has a non-negligible second singular value; a separable one is rank-1
    // (σ₂/σ₁ ~ 0).
    Eigen::MatrixXcd G(m, m);
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
      {
        G(i, j) = linalg::Dot(comm, *sf[i], *sf[j]) / (nf[i] * nf[j]);
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(G);
    const auto sv = es.eigenvalues();  // ascending
    const double sigma1 = std::sqrt(std::max(0.0, sv(m - 1)));
    const double sigma2 = std::sqrt(std::max(0.0, sv(m - 2)));
    const double rank_ratio = sigma2 / std::max(sigma1, 1e-300);
    INFO("s_full sampled subspace σ₂/σ₁ = " << rank_ratio);
    CHECK(rank_ratio > 1.0e-3);  // decisively above the rank-1 numerical floor (~1e-13)

    // (b) Center-freeze error: project each sample onto the 1-D span of the band-center
    // shape and measure the relative residual. A frozen (rank-1) model would represent the
    // whole band by this single direction, so a large residual quantifies the rotation it
    // misses.
    const int i_ref = m / 2;
    const std::complex<double> ref_sq =
        linalg::Dot(comm, *sf[i_ref], *sf[i_ref]);  // ‖s_ref‖²
    double max_resid = 0.0;
    for (int i = 0; i < m; i++)
    {
      // Least-squares coefficient c = <s_ref, s_i> / ‖s_ref‖². linalg::Dot(a,b) = <b,a>
      // (conjugate on the first arg), so <s_ref, s_i> = Dot(s_i, s_ref).
      const std::complex<double> proj = linalg::Dot(comm, *sf[i], *sf[i_ref]) / ref_sq;
      ComplexVector d(*sf[i]);
      d.AXPY(-proj, *sf[i_ref]);  // s_i − proj·s_ref
      max_resid = std::max(max_resid, linalg::Norml2(comm, d) / nf[i]);
    }
    INFO("center-freeze max relative residual over band = " << max_resid);
    CHECK(max_resid > 0.1);  // freezing the shape at band center is inadequate
  }
  REQUIRE(checked_any);  // the slab hybrid port mode must be active and tracked
}

// Wave ports on a parallel plate guide along z (width w, gap h in um) with a London slab of
// thickness d above the gap, compared to the exact TM0 mode index from the dispersion
// relation p tanh(p h) = -(q/eps) coth(q d) with a free back face (tanh(q d) for a PEC back
// face), where p^2 = k0^2 (n^2 - 1), q^2 = 1/lambda^2 + p^2, and eps = 1 - 1/(k0 lambda)^2.
// The gap keeps n_eff below the spectral shift of the port mode solve.
TEST_CASE("WavePortOperator London slab", "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  constexpr double lambda = 0.1, w = 1.0, h = 5.0, d = 0.1, l = 2.0, freq_ghz = 500.0;
  for (bool pec_back : {false, true})
  {
    json setup = {
        {"Problem", {{"Type", "Driven"}, {"Output", ""}}},
        {"Model", {{"Mesh", "london_slab.mesh"}, {"L0", 1.0e-6}}},
        {"Domains",
         {{"Materials", json::array({{{"Attributes", {1}}},
                                     {{"Attributes", {2}}, {"LondonDepth", lambda}}})}}},
        {"Boundaries",
         {{"PEC", {{"Attributes", pec_back ? json({1, 2}) : json({1})}}},
          {"WavePort",
           json::array({{{"Index", 1}, {"Attributes", {4}}, {"Excitation", true}},
                        {{"Index", 2}, {"Attributes", {5}}}})}}},
        {"Solver",
         {{"Order", 2},
          {"Driven",
           {{"Samples", json::array({{{"Type", "Point"}, {"Freq", {freq_ghz}}}})}}}}}};
    IoData iodata(setup, /*print=*/false);

    // Hexahedra: attribute 1 for the gap and 2 for the slab; boundary attributes 1 (y = 0),
    // 2 (y = h + d), 3 (x = 0, w), 4 (z = 0), and 5 (z = l).
    constexpr int nx = 2, ny_gap = 8, ny_slab = 8, nz = 2, ny = ny_gap + ny_slab;
    auto smesh = std::make_unique<mfem::Mesh>(
        3, (nx + 1) * (ny + 1) * (nz + 1), nx * ny * nz, 2 * (nx * nz + ny * nz + nx * ny));
    for (int k = 0; k <= nz; k++)
    {
      for (int j = 0; j <= ny; j++)
      {
        const double y = (j <= ny_gap) ? h * j / ny_gap : h + d * (j - ny_gap) / ny_slab;
        for (int i = 0; i <= nx; i++)
        {
          smesh->AddVertex(w * i / nx, y, l * k / nz);
        }
      }
    }
    auto v = [&](int i, int j, int k) { return i + (nx + 1) * (j + (ny + 1) * k); };
    for (int k = 0; k < nz; k++)
    {
      for (int j = 0; j < ny; j++)
      {
        for (int i = 0; i < nx; i++)
        {
          smesh->AddHex(v(i, j, k), v(i + 1, j, k), v(i + 1, j + 1, k), v(i, j + 1, k),
                        v(i, j, k + 1), v(i + 1, j, k + 1), v(i + 1, j + 1, k + 1),
                        v(i, j + 1, k + 1), (j < ny_gap) ? 1 : 2);
        }
      }
    }
    for (int k = 0; k < nz; k++)
    {
      for (int i = 0; i < nx; i++)
      {
        smesh->AddBdrQuad(v(i, 0, k), v(i, 0, k + 1), v(i + 1, 0, k + 1), v(i + 1, 0, k),
                          1);
        smesh->AddBdrQuad(v(i, ny, k), v(i + 1, ny, k), v(i + 1, ny, k + 1),
                          v(i, ny, k + 1), 2);
      }
      for (int j = 0; j < ny; j++)
      {
        smesh->AddBdrQuad(v(0, j, k), v(0, j + 1, k), v(0, j + 1, k + 1), v(0, j, k + 1),
                          3);
        smesh->AddBdrQuad(v(nx, j, k), v(nx, j, k + 1), v(nx, j + 1, k + 1),
                          v(nx, j + 1, k), 3);
      }
    }
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        smesh->AddBdrQuad(v(i, j, 0), v(i, j + 1, 0), v(i + 1, j + 1, 0), v(i + 1, j, 0),
                          4);
        smesh->AddBdrQuad(v(i, j, nz), v(i + 1, j, nz), v(i + 1, j + 1, nz),
                          v(i, j + 1, nz), 5);
      }
    }
    smesh->FinalizeHexMesh(1, 1, true);
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
    iodata.NondimensionalizeInputs(smesh);
    std::vector<std::unique_ptr<Mesh>> mesh_vec;
    mesh_vec.push_back(
        std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *smesh)));
    SpaceOperator space_op(iodata, mesh_vec);

    auto &wp_op = space_op.GetWavePortOp();
    REQUIRE(wp_op.Size() == 2);
    const double omega =
        2.0 * std::numbers::pi *
        iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_ghz);
    wp_op.InitializeModalReference(omega);
    const double k0 =
        2.0 * std::numbers::pi * freq_ghz * 1.0e9 / electromagnetics::c0_ * 1.0e-6;
    const double eps = 1.0 - 1.0 / (k0 * k0 * lambda * lambda);
    auto F = [&](double n)
    {
      const double p = k0 * std::sqrt(n * n - 1.0);
      const double q = std::sqrt(1.0 / (lambda * lambda) + p * p);
      const double g = pec_back ? std::tanh(q * d) : 1.0 / std::tanh(q * d);
      return p * std::tanh(p * h) + q * g / eps;
    };
    double lo = 1.0 + 1.0e-12, hi = 10.0;
    for (int it = 0; it < 200; it++)
    {
      const double mid = 0.5 * (lo + hi);
      ((F(mid) < 0.0) ? lo : hi) = mid;
    }
    const double n_exact = 0.5 * (lo + hi);
    for (const auto &[idx, data] : wp_op)
    {
      const double n_eff = data.kn0.real() / omega;
      CAPTURE(pec_back, idx, n_exact, n_eff, data.kn0.imag() / omega);
      CHECK(std::abs(data.kn0.imag()) < 1.0e-6 * std::abs(data.kn0.real()));
      CHECK_THAT(n_eff, WithinRel(n_exact, 1.0e-7));
    }
  }
}

// The traces of cpw_wave_uniform are interior PEC sheets crossing the wave ports when the
// mesh is not cracked. Refining the elements next to the ports on one side of the trace
// plane makes the port submeshes nonconforming, with master edges on the unrefined side of
// the trace lines which have no boundary elements: the port DoFs on the PEC segments must
// be essential (masters included), and no others.
TEST_CASE("WavePortOperator-InteriorPECOnNonconformingPorts",
          "[waveportoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  json setup = LoadCpwWaveConfig();
  setup["Model"]["CrackInternalBoundaryElements"] = false;
  IoData iodata(setup, /*print=*/false);

  std::vector<std::unique_ptr<Mesh>> mesh_io;
  {
    auto smesh = mesh::Load(iodata, comm);
    if (smesh)
    {
      // Plane of the interior PEC boundary elements (coordinate axis and value).
      const auto &pec = iodata.boundaries.pec.attributes;
      mfem::Vector lo(3), hi(3);
      lo = mfem::infinity();
      hi = -mfem::infinity();
      for (int be = 0; be < smesh->GetNBE(); be++)
      {
        int e1, e2;
        smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
        if (e2 < 0 ||
            std::find(pec.begin(), pec.end(), smesh->GetBdrAttribute(be)) == pec.end())
        {
          continue;
        }
        mfem::Array<int> v;
        smesh->GetBdrElementVertices(be, v);
        for (int i : v)
        {
          for (int d = 0; d < 3; d++)
          {
            lo(d) = std::min(lo(d), smesh->GetVertex(i)[d]);
            hi(d) = std::max(hi(d), smesh->GetVertex(i)[d]);
          }
        }
      }
      int axis = -1;
      for (int d = 0; d < 3; d++)
      {
        if (hi(d) - lo(d) < 1.0e-9 * (1.0 + std::abs(lo(d))))
        {
          axis = d;
        }
      }
      REQUIRE(axis >= 0);

      // Refine the elements next to the wave ports on one side of the trace plane.
      std::vector<int> port_attrs;
      for (const auto &[idx, data] : iodata.boundaries.waveport)
      {
        port_attrs.insert(port_attrs.end(), data.attributes.begin(), data.attributes.end());
      }
      smesh->EnsureNCMesh(true);
      mfem::Array<int> marked;
      mfem::Vector c(3);
      for (int be = 0; be < smesh->GetNBE(); be++)
      {
        if (std::find(port_attrs.begin(), port_attrs.end(), smesh->GetBdrAttribute(be)) ==
            port_attrs.end())
        {
          continue;
        }
        int e1, e2;
        smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
        smesh->GetElementCenter(e1, c);
        if (c(axis) > lo(axis))
        {
          marked.Append(e1);
        }
      }
      marked.Sort();
      marked.Unique();
      smesh->GeneralRefinement(marked, 1, 0);
    }
    if (iodata.model.Lc <= 0.0)
    {
      iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
    }
    iodata.NondimensionalizeInputs(smesh);
    mesh_io.push_back(
        std::make_unique<Mesh>(mesh::Partition(iodata, std::move(smesh), comm)));
  }
  SpaceOperator space_op(iodata, mesh_io);
  const auto &wp_op = space_op.GetWavePortOp();
  REQUIRE(wp_op.Size() > 0);

  // Numbers of checked true DoFs on the Dirichlet segments and off them, and on master
  // edges.
  int counts[3] = {0, 0, 0};
  for (const auto &[idx, data] : wp_op)
  {
    const auto &nd_fes = data.GetNDSpace().Get();
    const auto &h1_fes = data.GetH1Space().Get();
    const auto &port_mesh = *nd_fes.GetParMesh();
    REQUIRE(port_mesh.Nonconforming());
    const int nd_size = nd_fes.GetTrueVSize();
    std::vector<bool> nd_dbc(nd_size, false), h1_dbc(h1_fes.GetTrueVSize(), false);
    for (auto t : data.GetDbcTDofList())
    {
      if (t < nd_size)
      {
        nd_dbc[t] = true;
      }
      else
      {
        h1_dbc[t - nd_size] = true;
      }
    }

    // Segments of the Dirichlet boundary elements of the port submesh (PEC, AuxPEC and the
    // other wave ports, as in the WavePortOperator constructor), gathered from all ranks.
    // Coordinates are taken from the mesh nodes: the vertex coordinates of nonconforming
    // submeshes do not follow their vertex numbering.
    std::vector<int> dbc(iodata.boundaries.pec.attributes.begin(),
                         iodata.boundaries.pec.attributes.end());
    dbc.insert(dbc.end(), iodata.boundaries.auxpec.attributes.begin(),
               iodata.boundaries.auxpec.attributes.end());
    for (const auto &[other_idx, other_data] : iodata.boundaries.waveport)
    {
      if (other_idx != idx && other_data.active)
      {
        dbc.insert(dbc.end(), other_data.attributes.begin(), other_data.attributes.end());
      }
    }
    const int sdim = port_mesh.SpaceDimension();
    std::vector<double> segs;
    for (int be = 0; be < port_mesh.GetNBE(); be++)
    {
      if (std::find(dbc.begin(), dbc.end(), port_mesh.GetBdrAttribute(be)) == dbc.end())
      {
        continue;
      }
      mfem::Array<int> v;
      port_mesh.GetBdrElementVertices(be, v);
      for (int i : v)
      {
        double x[3];
        port_mesh.GetNode(i, x);
        segs.insert(segs.end(), x, x + sdim);
      }
    }
    int nloc = static_cast<int>(segs.size()), nproc = Mpi::Size(comm);
    std::vector<int> sizes(nproc), displs(nproc, 0);
    MPI_Allgather(&nloc, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
    for (int r = 1; r < nproc; r++)
    {
      displs[r] = displs[r - 1] + sizes[r - 1];
    }
    std::vector<double> all_segs(displs.back() + sizes.back());
    MPI_Allgatherv(segs.data(), nloc, MPI_DOUBLE, all_segs.data(), sizes.data(),
                   displs.data(), MPI_DOUBLE, comm);
    auto OnDbc = [&](const double *x)
    {
      for (std::size_t k = 0; k < all_segs.size(); k += 2 * sdim)
      {
        const double *a = &all_segs[k], *b = &all_segs[k + sdim];
        double ab2 = 0.0, ax_ab = 0.0, ax2 = 0.0;
        for (int d = 0; d < sdim; d++)
        {
          ab2 += (b[d] - a[d]) * (b[d] - a[d]);
          ax_ab += (x[d] - a[d]) * (b[d] - a[d]);
          ax2 += (x[d] - a[d]) * (x[d] - a[d]);
        }
        const double t = ax_ab / ab2;
        if (t > -1.0e-9 && t < 1.0 + 1.0e-9 && ax2 - t * t * ab2 < 1.0e-12 * ab2)
        {
          return true;
        }
      }
      return false;
    };
    auto Check = [&](const mfem::ParFiniteElementSpace &fes, const std::vector<bool> &ess,
                     const mfem::Array<int> &dofs, bool on_dbc, bool master)
    {
      for (auto d : dofs)
      {
        const int t = fes.GetLocalTDofNumber((d >= 0) ? d : -1 - d);
        if (t < 0)
        {
          continue;
        }
        CHECK(ess[t] == on_dbc);
        counts[on_dbc ? 0 : 1]++;
        counts[2] += (on_dbc && master);
      }
    };
    const auto &edge_list = port_mesh.ncmesh->GetEdgeList();
    mfem::Array<int> v, dofs;
    mfem::Vector mid(sdim), x0(sdim), x1(sdim);
    for (int e = 0; e < port_mesh.GetNEdges(); e++)
    {
      port_mesh.GetEdgeVertices(e, v);
      port_mesh.GetNode(v[0], x0.GetData());
      port_mesh.GetNode(v[1], x1.GetData());
      add(0.5, x0, 0.5, x1, mid);
      const bool on_dbc = OnDbc(mid.GetData());
      const bool master =
          edge_list.GetMeshIdAndType(e).type == mfem::NCMesh::NCList::MeshIdType::MASTER;
      nd_fes.GetEdgeInteriorDofs(e, dofs);
      Check(nd_fes, nd_dbc, dofs, on_dbc, master);
      h1_fes.GetEdgeInteriorDofs(e, dofs);
      Check(h1_fes, h1_dbc, dofs, on_dbc, master);
    }
    for (int i = 0; i < port_mesh.GetNV(); i++)
    {
      h1_fes.GetVertexDofs(i, dofs);
      port_mesh.GetNode(i, x0.GetData());
      Check(h1_fes, h1_dbc, dofs, OnDbc(x0.GetData()), false);
    }
  }
  Mpi::GlobalSum(3, counts, comm);
  CHECK(counts[0] > 0);
  CHECK(counts[1] > 0);
  CHECK(counts[2] > 0);
}
