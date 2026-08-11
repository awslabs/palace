// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification for the electrostatic operator.
//
// We manufacture a smooth potential on the unit cube [0,1]³ that vanishes on the boundary:
//     V_mms(x,y,z) = sin(πx) sin(πy) sin(πz)
// Since ∇²(sin πx sin πy sin πz) = -3π² V_mms, the volumetric charge that makes V_mms an exact
// solution of  -∇·(ε ∇V) = ρ  (constant ε) is
//     ρ_mms = -ε ∇²V_mms = 3π²ε V_mms.
// We solve -∇·(ε ∇V) = ρ_mms with homogeneous Dirichlet BCs (all six faces grounded, matching
// V_mms = 0 on ∂Ω) and check that the discrete solution converges to V_mms in the L2 norm.

// STUDY NOTE: all `// STUDY:` comments below are a learning layer — REMOVE before committing.

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>                                       // STUDY: the FE library Palace wraps
#include <catch2/catch_test_macros.hpp>                   // STUDY: TEST_CASE, CHECK
#include <catch2/generators/catch_generators_all.hpp>     // STUDY: GENERATE / from_range
#include <catch2/matchers/catch_matchers_floating_point.hpp>  // STUDY: WithinAbs matcher
#include "linalg/ksp.hpp"                                 // STUDY: KspSolver (iterative linear solver)
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"                              // STUDY: palace::Vector (device-aware)
#include "models/laplaceoperator.hpp"                     // STUDY: builds electrostatic stiffness matrix
#include "utils/communication.hpp"                        // STUDY: Mpi::World()
#include "utils/iodata.hpp"                               // STUDY: config object (materials, BCs, solver)
#include "utils/units.hpp"                                // STUDY: nondimensionalization

using namespace palace;
using namespace Catch::Matchers;

namespace
{

constexpr double kEpsilonR = 1.0;  // constant relative permittivity on the domain
// STUDY: ε is a single constant everywhere — that's what lets ρ_mms = -ε∇²V (ε pulls out of the
// STUDY: divergence). Spatially-varying/tensor ε would add a ∇ε·∇V term (a future phase).

// Manufactured solution V_mms = sin(πx) sin(πy) sin(πz).
// STUDY: We DECLARE this to be the exact answer. Chosen because it's 0 on all faces of [0,1]³
// STUDY: (matches grounded BC for free) and is a Laplacian eigenfunction (∇² of it = -3π² times
// STUDY: itself), making the source a trivial multiple of the solution.
double VmmsFunction(const mfem::Vector &x)
{
  return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]) * std::sin(M_PI * x[2]);
}

// Manufactured source ρ_mms = -ε ∇²V_mms = 3π²ε · V_mms (constant ε).
// STUDY: 3 = num dimensions (each sin contributes -π² to the Laplacian). The two minus signs
// STUDY: (-ε and -3π²) cancel → positive source proportional to V_mms.
double RhommsFunction(const mfem::Vector &x)
{
  return 3.0 * M_PI * M_PI * kEpsilonR * VmmsFunction(x);
}

// Solve the manufactured electrostatic problem on an N×N×N unit cube and return the L2 error
// ‖V_h - V_mms‖ against the manufactured solution.
// STUDY: The workhorse — one call = one full solve at resolution n and polynomial order.
// STUDY: The convergence test calls it repeatedly with different (n, order).
double SolveMmsL2Error(int n, int order)
{
  MPI_Comm comm = Mpi::World();  // STUDY: MPI communicator (all ranks); works in serial too.

  // Build the IoData in-memory: constant permittivity on the single domain attribute, all six
  // cube faces grounded (Dirichlet), no mesh file.
  // STUDY: IoData is normally parsed from JSON; here we build it by hand (units-only ctor).
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;                             // STUDY: reference length = 1; mesh coords unchanged.
  iodata.problem.type = ProblemType::ELECTROSTATIC;  // STUDY: selects the physics for MaterialOperator.
  auto &material = iodata.domains.materials.emplace_back();  // STUDY: add one material...
  material.attributes = {1};                         // STUDY: ...on mesh domain attribute 1.
  material.epsilon_r.s = {kEpsilonR, kEpsilonR, kEpsilonR};  // STUDY: .s = isotropic scalar ε.
  iodata.boundaries.pec.attributes = {1, 2, 3, 4, 5, 6};     // STUDY: ground all 6 faces = Dirichlet everywhere.
  iodata.solver.order = order;                       // STUDY: polynomial degree p of the H1 elements.
  iodata.solver.linear.tol = 1.0e-9;                 // STUDY: drive the ITERATIVE-SOLVER (algebraic) error
  // STUDY: below the discretization error we measure. The default 1e-6 leaves an algebraic "error
  // STUDY: floor" that the finest/highest-order discretization error (P=3, N=32 ≈ 6e-8) falls
  // STUDY: below, collapsing the measured rate (the P=3 16→32 failure; also seen as serial≠mpi
  // STUDY: error at the same mesh). 1e-9 sits ~1.5 decades under it — enough margin, not overkill.
  // STUDY: max_it is left at its default (100): PCG converges in 3-9 iters, never near the cap.

  // Structured unit-cube mesh; MakeCartesian3D assigns domain attribute 1 and face attributes
  // 1..6, so h = 1/n.
  // STUDY: mesh GENERATED in-code (not a file) so refining = changing n; h = 1/n exactly.
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      n, n, n, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));  // STUDY: n×n×n hex cells on [0,1]³.
  iodata.NondimensionalizeInputs(serial_mesh);       // STUDY: scale mesh/materials into internal units.
  iodata.CheckConfiguration();  // also initializes DefaultIntegrationOrder used by integrators
  // STUDY: ^ NECESSARY: besides validating, it inits a GLOBAL quadrature-order default the
  // STUDY: integrators read. Skipping it made an earlier version assert during assembly.
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);  // STUDY: partition across ranks.
  std::vector<std::unique_ptr<Mesh>> mesh;           // STUDY: LaplaceOperator wants a vector (MG levels)...
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));  // STUDY: ...wrap ParMesh in palace::Mesh; one level.

  LaplaceOperator laplace_op(iodata, mesh);          // STUDY: builds FE spaces + operator -div(ε∇·).
  auto K = laplace_op.GetStiffnessMatrix();          // STUDY: discrete stiffness K (Dirichlet rows = DIAG_ONE).
  auto &h1_fespace = laplace_op.GetH1Space();        // STUDY: the scalar H1 space for the potential V.

  // Assemble the volumetric source ρ_mms into the RHS. Use mfem::DomainLFIntegrator (a true
  // volume integrator) -- NOT palace::DomainLFIntegrator, which is a boundary integrator.
  // STUDY: a "linear form" is the RHS vector b_i = ∫ ρ_mms·φ_i dV. DomainLFIntegrator does that
  // STUDY: volume integral. GOTCHA: palace::DomainLFIntegrator is aliased to a BOUNDARY
  // STUDY: integrator, so we must fully-qualify mfem::DomainLFIntegrator for a volume source.
  mfem::FunctionCoefficient rho_coef(RhommsFunction);  // STUDY: wraps ρ_mms(x) for mfem to sample.
  mfem::ParLinearForm lf(&h1_fespace.Get());         // STUDY: the RHS linear form on the H1 space.
  lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(rho_coef));  // STUDY: add ∫ρφ (mfem owns the ptr).
  lf.Assemble();                                     // STUDY: compute the integrals.

  Vector RHS(h1_fespace.GetTrueVSize()), V(h1_fespace.GetTrueVSize());
  // STUDY: "true dofs" = the independent unknowns; RHS = b, V = solution of K V = b.
  RHS.UseDevice(true);                               // STUDY: allow GPU storage (no-op on CPU).
  V.UseDevice(true);
  h1_fespace.GetProlongationMatrix()->MultTranspose(lf, RHS);
  // STUDY: lf lives on all local dofs; Pᵀ maps it down to the TRUE dofs → RHS.

  // Homogeneous Dirichlet: zero the RHS at essential (boundary) true dofs. K uses DIAG_ONE on
  // these rows, so the solve pins V = 0 there, matching V_mms = 0 on ∂Ω.
  const mfem::ParMesh &pmesh = h1_fespace.GetParMesh();
  mfem::Array<int> all_bdr(pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0);
  // STUDY: all_bdr = marker over boundary-attribute numbers; size = max attribute present.
  mfem::Array<int> ess_tdofs;                        // STUDY: will hold the boundary true-dof indices.
  all_bdr = 1;                                       // STUDY: mark EVERY boundary attribute essential.
  h1_fespace.Get().GetEssentialTrueDofs(all_bdr, ess_tdofs);  // STUDY: which true dofs are on those faces.
  RHS.SetSubVector(ess_tdofs, 0.0);                  // STUDY: zero RHS there → V=0 on the boundary.

  KspSolver ksp(iodata, laplace_op.GetH1Spaces());   // STUDY: iterative Krylov solver + preconditioner.
  ksp.SetOperators(*K, *K);                          // STUDY: solve with K; build preconditioner from K.
  V = 0.0;                                           // STUDY: initial guess.
  ksp.Mult(RHS, V);                                  // STUDY: THE SOLVE: K V = RHS → V.

  // Recover the grid function and compute the L2 error against V_mms.
  // STUDY: V is raw dof values; to integrate the error we need a field defined everywhere.
  mfem::ParGridFunction V_gf(&h1_fespace.Get());
  V_gf.SetFromTrueDofs(V);                           // STUDY: expand true dofs → full field.

  const int order_quad = order + 2;                  // STUDY: integrate error 2 orders above the FE order,
  // STUDY: so quadrature error doesn't contaminate the measured discretization error.
  const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];  // STUDY: one rule per element geometry.
  for (int i = 0; i < mfem::Geometry::NumGeom; i++)
  {
    irs[i] = &mfem::IntRules.Get(i, order_quad);
  }
  mfem::FunctionCoefficient v_exact(VmmsFunction);   // STUDY: exact V_mms to compare against.
  return V_gf.ComputeL2Error(v_exact, irs);          // STUDY: returns ‖V_h - V_mms‖_L2.
}

// Observed convergence rate between two (mesh size h, L2 error) points: the slope of the error
// in log-log space, log(e0/e1) / log(h0/h1).
// STUDY: theory says error ≈ C·h^rate, so log(error) = logC + rate·log(h) — a line whose slope
// STUDY: IS the rate. This computes that slope between two refinement levels.
double ConvergenceRate(double h0, double e0, double h1, double e1)
{
  return std::log(e0 / e1) / std::log(h0 / h1);
}

// Convergence sweep configuration. Extend these lists to add polynomial orders or refinement
// levels; the test logic below adapts to any number of entries.
// STUDY: these two lists are the ONLY thing to edit to change coverage; the loops below are
// STUDY: generic, so adding e.g. 64 to kResolutions auto-adds a 32→64 rate check.
const std::vector<int> kOrders = {1, 2, 3};
const std::vector<int> kResolutions = {8, 16, 32};  // N per side; h = 1/N (must be ascending)

// Tolerance on the measured convergence rate vs. the theoretical p+1.
constexpr double kRateTol = 0.3;  // STUDY: measured slope must be within ±0.3 of p+1.

}  // namespace

// Single-resolution sanity check: the manufactured problem solves and the L2 error is small.
// STUDY: the cheap "does it even work" smoke test — one solve, one loose threshold.
TEST_CASE("Electrostatic MMS solves to a small L2 error", "[electrostaticmms][Serial][Parallel]")
{
  const double err = SolveMmsL2Error(/*n=*/16, /*order=*/2);
  CHECK(err < 1.0e-3);
}

// Convergence-rate verification: for each polynomial order p, the L2 error must decrease
// monotonically under mesh refinement and every consecutive-level rate must match the
// theoretical p+1 (to within kRateTol). Checking every consecutive pair confirms the rate is
// stable across refinement, not just at the finest level. Order is parametrized so a failure
// reports which order regressed.
TEST_CASE("Electrostatic MMS converges at the optimal rate",
          "[electrostaticmms][Serial][Parallel][Long]")
{
  const int order = GENERATE(from_range(kOrders));
  // STUDY: GENERATE re-runs this whole TEST_CASE once per order (1, then 2, then 3), each a
  // STUDY: separate reported sub-case — so a failure says exactly which order broke.
  CAPTURE(order);  // STUDY: print 'order' if an assertion below fails.

  // Solve on each refinement level and record (h, error).
  std::vector<double> h(kResolutions.size()), err(kResolutions.size());
  for (std::size_t i = 0; i < kResolutions.size(); i++)
  {
    h[i] = 1.0 / kResolutions[i];                    // STUDY: mesh size for level i.
    err[i] = SolveMmsL2Error(kResolutions[i], order);  // STUDY: full solve → L2 error at (N, order).
  }

  // Every consecutive pair: error decreases and the observed rate matches p+1.
  const double expected_rate = order + 1.0;          // STUDY: theoretical L2 rate for degree-p H1.
  for (std::size_t i = 0; i + 1 < kResolutions.size(); i++)  // STUDY: loop adjacent pairs (0-1, 1-2, …).
  {
    CAPTURE(kResolutions[i], kResolutions[i + 1], err[i], err[i + 1]);  // STUDY: print on failure.
    CHECK(err[i + 1] < err[i]);                      // STUDY: finer mesh ⇒ smaller error (monotone).
    const double rate = ConvergenceRate(h[i], err[i], h[i + 1], err[i + 1]);  // STUDY: measured slope.
    CAPTURE(rate, expected_rate);
    CHECK_THAT(rate, WithinAbs(expected_rate, kRateTol));  // STUDY: |measured - (p+1)| ≤ 0.3.
  }
}
