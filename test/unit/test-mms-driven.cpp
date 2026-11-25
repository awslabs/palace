// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification of the driven (frequency-domain)
// operator, following Garcia-Doñoro et al., "Verification Process of Finite-Element Method
// Code for Electromagnetics: Using the Method of Manufactured Solutions," IEEE Antennas
// Propag. Mag. 58(2), 2016. The paper verifies the formulation one term at a time; this
// file covers three of its steps on a common setup, named here (and in the Problem
// enumerators) by what they solve:
//   - the mass projection (the paper's Test 1, eqs. 15–17): the mass term alone;
//   - the PEC problem (Test 3, eqs. 20–22): curl-curl + mass, PEC on every face;
//   - the Neumann-face problem (Test 4, eqs. 23–24): the PEC problem with one natural face.
// The curl-curl term alone (Test 2) is the magnetostatic test (test-mms-magnetostatic.cpp).
// Throughout, the material is lossless with anisotropic μ_r and ε_r, and the angular
// frequency ω is real.
//
// PEC problem. As in the paper, the manufactured field E_mms defines a volume source
// through the strong form of the PDE,
//   J = ∇×(μ⁻¹∇×E_mms) − ω²εE_mms,
// and the code solves
//   ∫ μ⁻¹(∇×E)·(∇×φ) dΩ − ω² ∫ εE·φ dΩ = ∫ J·φ dΩ,   n×E = n×E_mms = 0 on ∂Ω.
// Integrating the curl term by parts recovers the left side for E = E_mms, with a boundary
// term ∮ (n×μ⁻¹∇×E_mms)·φ dS that vanishes because n×φ = 0 on the PEC boundary. The
// right-hand side therefore checks that Palace's discrete operator solves the PDE itself,
// not only that it matches another evaluation of the same bilinear form (the weak-curl
// right-hand side of the magnetostatic test, which is required there because the
// operator is singular). Unlike the magnetostatic case, the operator is nonsingular when
// ω² is not a cavity eigenvalue, so E is unique and both the field error ‖E_h − E_mms‖
// and the curl error ‖∇×E_h − ∇×E_mms‖ (absolute versions of the paper's e_f and e_c) are
// asserted.
//
// Mass projection. Every differential operator is off, which checks the mass term, the
// volume source, and the PEC elimination in isolation. It is the ε-weighted L2 projection
// of E_mms onto ND_p,
//   ∫ εE·φ dΩ = ∫ εE_mms·φ dΩ,   n×E = n×E_mms = 0 on ∂Ω,
// with Palace's driven mass matrix M (the paper's version is the ε = 1 case). There are no
// derivatives, so strong and weak forms coincide, and μ and ω play no part. Only the field
// error is asserted: the L2 projection is the best ε-weighted approximation of E_mms, so
// ‖E_h − E_mms‖ = O(h^p), but nothing controls its curl, which is only reported.
//
// Neumann-face problem. The PEC problem's equation with the x = 1 face natural (PMC)
// instead of PEC. Integrating by parts no longer drops the boundary term there, so the
// right-hand side gains the surface integral −∮_ΓN (n×μ⁻¹∇×E_mms)·φ dS, the inhomogeneous
// Neumann data that the homogeneous natural condition in Palace's weak form must balance.
// Its fields vanish tangentially only on the five PEC faces, so the trace of E_h on Γ_N is
// a genuine unknown. Both errors are asserted, as for the PEC problem.
//
// Backend-specific tests supply E_mms, ∇×E_mms, and ∇×(μ⁻¹∇×E_mms) through MmsCase; this
// file owns the common Palace solve and error checks.

#include "test-mms-driven.hpp"
#include "test-mms-common.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;
using namespace std::complex_literals;

namespace palace::testing::mms::driven
{
namespace
{

// Which of the paper's tests to run on the common setup.
enum class Problem
{
  MASS_PROJECTION,        // M E = b, the ε-weighted L2 projection.
  CURL_CURL_MASS,         // (K − ω²M) E = b, the driven operator, all faces PEC.
  CURL_CURL_MASS_NEUMANN  // As CURL_CURL_MASS, but the x = 1 face is natural (PMC).
};

// The natural face of the Neumann-face problem. MakeCartesian3D numbers the unit-cube
// faces z = 0: 1, y = 0: 2, x = 1: 3, y = 1: 4, x = 0: 5, z = 1: 6, for hexahedra and
// tetrahedra alike, so attribute 3 is x = 1, with outward normal n = +x̂.
constexpr int kNeumannAttribute = 3;

// Optional surface term of the right-hand side: a coefficient integrated against φ over
// the marked boundary faces with the given quadrature rule (see SolveMmsError).
struct BoundarySource
{
  mfem::VectorCoefficient &data;
  mfem::Array<int> &marker;
  const mfem::IntegrationRule &ir;
};

struct MmsError
{
  HYPRE_BigInt primary_dofs;
  double electric;  // ‖E_h − E_mms‖, including the (ideally zero) imaginary part.
  double curl;      // ‖∇×E_h − ∇×E_mms‖, likewise.
};

// Solve the frequency-domain system exactly as DrivenSolver does, with only the
// manufactured right-hand side assembled here:
//   A = K + iω C − ω² M,   P = the production preconditioner matrix for the same ω,
// with K assembled using DIAG_ONE and C, M using DIAG_ZERO, so the PEC rows of A are the
// identity. The Dirichlet data is zero, so zeroing the PEC entries of the right-hand side
// (as SpaceOperator::GetExcitationVector does) imposes n×E = 0. The material is lossless
// and there are no ports, so C and the extra matrix A2(ω) are absent (both asserted
// below, so the solved system cannot silently drift from DrivenSolver's), and the
// right-hand side is real:
//   RHS = R_ndᵀ [ ∫ J·φ dΩ + ∮_Γ g·φ dS ],   R_nd = the ND prolongation (local to true
// dofs), where the surface term with data g is present only for the Neumann-face problem.
// The integrals use the caller's quadrature rules, see SolveMmsError. On a boundary face,
// MFEM's VectorFEDomainLFIntegrator evaluates the physical (tangential) trace of the
// Nédélec basis, so g·φ only sees the tangential part of g.
ComplexVector SolveManufactured(const IoData &iodata, SpaceOperator &space_op, double omega,
                                mfem::VectorCoefficient &source,
                                const mfem::IntegrationRule &ir,
                                const BoundarySource *boundary = nullptr)
{
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(!C);
  REQUIRE(!space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO));
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ksp.SetOperators(*A, *P);

  auto &nd = space_op.GetNDSpace();
  mfem::ParLinearForm b(&nd.Get());
  b.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(source, &ir));
  if (boundary)
  {
    b.AddBoundaryIntegrator(
        new mfem::VectorFEDomainLFIntegrator(boundary->data, &boundary->ir),
        boundary->marker);
  }
  b.Assemble();

  ComplexVector E(nd.GetTrueVSize()), RHS(nd.GetTrueVSize());
  E.UseDevice(true);
  RHS.UseDevice(true);
  E = 0.0;
  RHS = 0.0;
  nd.GetProlongationMatrix()->MultTranspose(b, RHS.Real());
  linalg::SetSubVector(RHS, space_op.GetNDDbcTDofLists().back(), 0.0);

  ksp.Mult(RHS, E);
  // KspSolver::Mult only warns at max_it; see SolveEliminated in test-mms-common.hpp.
  REQUIRE(ksp.GetConverged());
  return E;
}

// Solve the mass projection with SpaceOperator's mass matrix,
//   M = ∫ εE·φ dΩ,   assembled with DIAG_ONE so the PEC rows are the identity,
// and the same right-hand side treatment as SolveManufactured:
//   RHS = R_ndᵀ [ ∫ εE_mms·φ dΩ ],   with the PEC entries zeroed.
// This is the real instantiation GetMassMatrix<Operator>. For a lossless material it
// equals the driven M: both share the real-part assembly, and the complex instantiation
// only adds an imaginary part for lossy materials. Its complex form is exercised by the
// PEC problem, inside K − ω²M.
//
// M is SPD, and on a quasi-uniform mesh its Jacobi-scaled condition number is bounded
// independently of h, so CG with a Jacobi preconditioner (configured by the caller) needs
// a mesh-independent number of iterations. That number still grows with p, steeply for the
// tetrahedral Nédélec basis. Measured: 3–6 iterations on hexahedra at every order; on
// tetrahedra about 20, 65, and 115 for orders 1–3 at the 1e-9 convergence tolerance, and
// about 800 for order 5 at the 1e-13 exactness tolerance (reduction factor ≈0.96 per
// iteration). The solution is returned as a complex vector with zero imaginary part, so
// the error computation is shared with the curl-curl problems.
ComplexVector SolveMassProjection(const IoData &iodata, SpaceOperator &space_op,
                                  mfem::VectorCoefficient &source,
                                  const mfem::IntegrationRule &ir)
{
  auto M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);
  REQUIRE(M);
  KspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ksp.SetOperators(*M, *M);

  auto &nd = space_op.GetNDSpace();
  mfem::ParLinearForm b(&nd.Get());
  b.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(source, &ir));
  b.Assemble();

  Vector X(nd.GetTrueVSize()), RHS(nd.GetTrueVSize());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  nd.GetProlongationMatrix()->MultTranspose(b, RHS);
  linalg::SetSubVector(RHS, space_op.GetNDDbcTDofLists().back(), 0.0);

  ksp.Mult(RHS, X);
  REQUIRE(ksp.GetConverged());

  ComplexVector E(nd.GetTrueVSize());
  E.UseDevice(true);
  E = 0.0;
  E.Real() = X;
  return E;
}

// L2 error of a complex ND (or RT) true-dof vector against a real exact field: the real
// part is compared with the exact field and the imaginary part with zero.
double ComplexL2Error(FiniteElementSpace &fespace, const ComplexVector &x,
                      mfem::VectorCoefficient &exact)
{
  mfem::ParGridFunction x_re(&fespace.Get()), x_im(&fespace.Get());
  x_re.SetFromTrueDofs(x.Real());
  x_im.SetFromTrueDofs(x.Imag());
  mfem::Vector zero(3);
  zero = 0.0;
  mfem::VectorConstantCoefficient zero_coef(zero);
  return std::hypot(x_re.ComputeL2Error(exact), x_im.ComputeL2Error(zero_coef));
}

MmsError ComputeMmsErrors(SpaceOperator &space_op, const ComplexVector &E,
                          mfem::VectorCoefficient &e_exact,
                          mfem::VectorCoefficient &curl_exact)
{
  // ∇×E_h, mapped exactly from the ND space into the RT space by the discrete curl. (The
  // driven solver scales this by −1/(iω) to get B; the curl itself is compared here.)
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector curl_E(Curl.Height());
  curl_E.UseDevice(true);
  curl_E = 0.0;
  Curl.Mult(E.Real(), curl_E.Real());
  Curl.Mult(E.Imag(), curl_E.Imag());

  return {space_op.GlobalTrueVSize(), ComplexL2Error(space_op.GetNDSpace(), E, e_exact),
          ComplexL2Error(space_op.GetRTSpace(), curl_E, curl_exact)};
}

MmsError SolveMmsError(const MmsCase &mms, Problem problem,
                       std::unique_ptr<mfem::Mesh> serial_mesh, int order,
                       double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  // Every face is PEC, so n×E = 0 is imposed on all of ∂Ω, except for the Neumann-face
  // problem, where the x = 1 face is listed as PMC instead. PMC adds no term to Palace's
  // operator (it is the natural condition, the same as leaving the face unassigned), but
  // listing it mirrors a real configuration.
  //
  // The default frequency-domain linear solver is a sparse direct solver when one is built
  // in, so pin an iterative solver to solve the same way on every build and keep the
  // finest meshes affordable. The two curl-curl problems use AMS, GMRES with the shifted
  // AMS preconditioner; at most 29 iterations were measured for either, even for the
  // order-5 exactness solves at 1e-13, serially and on 2 ranks, so the default max_it of
  // 100 leaves margin. The mass projection uses CG with Jacobi on a single level, needing
  // neither the auxiliary-space correction (there is no curl-curl null space to handle)
  // nor multigrid (the condition number of M does not grow with refinement). There the
  // default max_it of 100 is too few for high-order tetrahedra (see SolveMassProjection):
  // order 5 at 1e-13 takes about 800, so 2000 leaves margin, and each CG iteration on M is
  // cheap.
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.mu_r.s = mms.mu_r;
  material.epsilon_r.s = mms.epsilon_r;
  const bool neumann = (problem == Problem::CURL_CURL_MASS_NEUMANN);
  for (int i = 0; i < serial_mesh->bdr_attributes.Size(); i++)
  {
    const int attr = serial_mesh->bdr_attributes[i];
    auto &attributes = (neumann && attr == kNeumannAttribute)
                           ? iodata.boundaries.pmc.attributes
                           : iodata.boundaries.pec.attributes;
    attributes.push_back(attr);
  }
  iodata.solver.order = order;
  iodata.solver.linear.tol = linear_tol;
  switch (problem)
  {
    case Problem::MASS_PROJECTION:
      iodata.solver.linear.type = LinearSolver::JACOBI;
      iodata.solver.linear.krylov_solver = KrylovSolver::CG;
      iodata.solver.linear.mg_max_levels = 1;
      iodata.solver.linear.max_it = 2000;
      break;
    case Problem::CURL_CURL_MASS:
    case Problem::CURL_CURL_MASS_NEUMANN:
      iodata.solver.linear.type = LinearSolver::AMS;
      break;
  }
  // Quadrature for the source integral ∫ J·φ dΩ. Polynomial exactness requires this
  // integral to be exact, or Galerkin orthogonality no longer pins E_h = E_mms. For the
  // polynomial fields of all three tests, J has degree ≤ 2 per variable (total degree
  // ≤ 4), and φ ∈ ND_p has degree ≤ p per variable on hexahedra and ≤ p in total on
  // tetrahedra, so J·φ needs degree p + 2 per variable (Gauss tensor rule) or p + 4 in
  // total (simplex rule). MFEM's default 2p is exact with no slack for order-5 tetrahedra;
  // 2p + 4 covers both element types at every order with margin and also tightens the
  // trigonometric source. (Geometry from the serial mesh, which every rank holds in full,
  // so ranks that own no elements still get a valid rule.)
  const mfem::IntegrationRule &source_ir =
      mfem::IntRules.Get(serial_mesh->GetElementGeometry(0), 2 * order + 4);
  // The same order for the Neumann-face surface integral ∮ (n×H_mms)·φ dS on the faces
  // (squares or triangles). For that problem's polynomial, n×H_mms has degree ≤ 2 per
  // variable (total degree ≤ 3) on x = 1, and the tangential trace of φ has degree ≤ p, so
  // the integrand needs p + 2 per variable or p + 3 in total, well within 2p + 4.
  const mfem::IntegrationRule &boundary_ir =
      mfem::IntRules.Get(serial_mesh->GetBdrElementGeometry(0), 2 * order + 4);

  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  // The manufactured source enters only through the right-hand side, so there is no port
  // or surface current excitation. SpaceOperator(iodata, mesh) rejects a driven problem
  // without one; the component constructor builds the same operator without that check
  // (as in test-libceed.cpp).
  SpaceOperator space_op(iodata.solver, iodata.domains, iodata.boundaries,
                         iodata.problem.type, iodata.units, mesh);
  mfem::VectorFunctionCoefficient e_exact(/*vdim=*/3, mms.electric_field);
  mfem::VectorFunctionCoefficient curl_exact(/*vdim=*/3, mms.curl);
  // The volume source: J = ∇×(μ⁻¹∇×E_mms) − ω²εE_mms for both curl-curl problems (the
  // Neumann-face problem adds only the surface term below), and only the mass part εE_mms
  // (with coefficient 1 instead of −ω²) for the mass projection. The backend supplies the
  // curl-curl term, which carries the μ dependence and the second derivatives. Scaling the
  // mass term component-wise equals εE only because the default (identity) material axes
  // make diag(ε_r) diagonal in the global frame (and likewise for μ in the backends); a
  // rotated tensor would need R diag(ε_r) Rᵀ here.
  const bool mass_only = (problem == Problem::MASS_PROJECTION);
  const double mass_coef = mass_only ? 1.0 : -mms.omega * mms.omega;
  mfem::VectorFunctionCoefficient source(
      /*vdim=*/3,
      [&mms, mass_only, mass_coef](const mfem::Vector &x, mfem::Vector &f)
      {
        if (mass_only)
        {
          f = 0.0;
        }
        else
        {
          mms.curl_curl(x, f);
        }
        mfem::Vector e(3);
        mms.electric_field(x, e);
        for (int d = 0; d < 3; d++)
        {
          f(d) += mass_coef * mms.epsilon_r[d] * e(d);
        }
      });
  // Neumann-face surface data. Integrating the curl term by parts with
  //   ∫ (∇×H)·φ dΩ = ∫ H·(∇×φ) dΩ + ∮ (n×H)·φ dS
  // (from ∇·(H×φ) = φ·(∇×H) − H·(∇×φ)) leaves −∮_ΓN (n×H_mms)·φ dS on the right-hand side
  // with H_mms = μ⁻¹∇×E_mms; the PEC faces drop out because n×φ = 0 there. This is the
  // paper's eq. 23 with Ψ_N = n×(∇×V_mms) (eq. 24). The outward normal of the x = 1 face is
  // +x̂, so −n×H = −(0, −H_z, H_y) = (0, H_z, −H_y). Writing the normal out here, rather
  // than taking it from the boundary element, keeps the sign independent of how the mesh
  // orients its boundary elements. H = μ⁻¹∇×E is diagonal for the same identity-axes reason
  // as the source.
  mfem::VectorFunctionCoefficient neumann_data(
      /*vdim=*/3,
      [&mms](const mfem::Vector &x, mfem::Vector &g)
      {
        mfem::Vector curl(3);
        mms.curl(x, curl);
        g(0) = 0.0;
        g(1) = curl(2) / mms.mu_r[2];
        g(2) = -curl(1) / mms.mu_r[1];
      });
  // The marker selects the natural face; it is only built (and its attribute checked) for
  // the Neumann-face problem, so the others accept any boundary numbering.
  mfem::Array<int> neumann_marker;
  if (neumann)
  {
    REQUIRE(serial_mesh->bdr_attributes.Find(kNeumannAttribute) >= 0);
    neumann_marker.SetSize(serial_mesh->bdr_attributes.Max());
    neumann_marker = 0;
    neumann_marker[kNeumannAttribute - 1] = 1;
  }
  const BoundarySource neumann_source{neumann_data, neumann_marker, boundary_ir};

  const ComplexVector E =
      mass_only ? SolveMassProjection(iodata, space_op, source, source_ir)
                : SolveManufactured(iodata, space_op, mms.omega, source, source_ir,
                                    neumann ? &neumann_source : nullptr);
  return ComputeMmsErrors(space_op, E, e_exact, curl_exact);
}

MmsError
SolveCartesianMmsError(const MmsCase &mms, Problem problem, int n, int order,
                       double linear_tol = 1.0e-9,
                       mfem::Element::Type element_type = mfem::Element::HEXAHEDRON)
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(n, n, n, element_type, 1.0, 1.0, 1.0));
  return SolveMmsError(mms, problem, std::move(serial_mesh), order, linear_tol);
}

const std::vector<int> kOrders = {1, 2, 3};
// The complex GMRES solve with the shifted AMS preconditioner is several times more
// expensive than the magnetostatic CG solve on the same mesh, and order 3 on a 16³ mesh
// alone takes about a minute per element type, so order 3 stops at 8³. It also starts at
// 4³ like the other orders: a 2³ mesh (8 hexahedra, 48 tetrahedra) is too coarse to
// resolve sin πx, so the 2 → 4 rate is pre-asymptotic, and the discrete eigenvalues near
// ω² are least accurate there (see kOmega). Order 3 therefore has a single rate pair.
// The mass projection (a Jacobi-CG solve) is cheaper and has no resonance, so neither
// reason applies to it; it shares these grids so that its report rows line up with the
// curl-curl problems'.
std::vector<int> Resolutions(int order)
{
  return (order < 3) ? std::vector<int>{4, 8, 16} : std::vector<int>{4, 8};
}

void ReportMmsData(const std::string &dataset, const std::string &solution,
                   const std::string &element, int order,
                   const std::string &resolution_kind, int resolution, double h,
                   const MmsError &error)
{
  mms::ReportMmsData(dataset, solution, element, order, resolution_kind, resolution, h,
                     error.primary_dofs, error.electric, error.curl);
}

constexpr double kRateTol = 0.3;
constexpr double kBackendAgreementTol = 1.0e-10;

// Polynomial exactness for any of the problems; the dataset name keeps their rows apart.
void CheckExactness(const MmsCase &polynomial, Problem problem, const std::string &dataset)
{
  // All three operators are nonsingular, so when E_mms lies in the Nédélec space,
  // E_h = E_mms up to the iterative solver tolerance (Galerkin orthogonality with a unique
  // solution; for the mass projection, a projection reproduces anything already in the
  // space), and both the field and curl errors vanish. Unlike the magnetostatic check,
  // E_mms itself must be in ND_p, not only its curl in ∇×ND_p. Both polynomial fields (the
  // bubble field and its Neumann-face variant) have each component of degree ≤ 2 per
  // variable, in the matching Q_{1,2,2} × Q_{2,1,2} × Q_{2,2,1} ⊆ ND_2 pattern, and total
  // degree 4: order 2 suffices on hexahedra, while on tetrahedra ND_p contains the full
  // P_{p-1}, so order 5 is needed. A coarse n = 2 tetrahedral mesh keeps order 5 cheap.
  struct Element
  {
    const char *name;
    mfem::Element::Type type;
    int order, n;
  };
  const std::array<Element, 2> elements = {
      Element{"hexahedron", mfem::Element::HEXAHEDRON, /*order=*/2, /*n=*/4},
      Element{"tetrahedron", mfem::Element::TETRAHEDRON, /*order=*/5, /*n=*/2}};
  for (const auto &[name, type, order, n] : elements)
  {
    DYNAMIC_SECTION(name)
    {
      const auto error = SolveCartesianMmsError(polynomial, problem, n, order,
                                                /*linear_tol=*/1.0e-13, type);
      ReportMmsData(dataset, "polynomial", name, order, "N", n, 1.0 / n, error);
      CAPTURE(error.electric, error.curl);
      CHECK(error.electric < 1.0e-10);
      CHECK(error.curl < 1.0e-10);
    }
  }
}

// Convergence sweep for any of the problems. The field error is always asserted at
// O(h^p); the curl error only when assert_curl is set (the curl-curl problems), and
// otherwise it is reported only (the mass projection).
void CheckConvergence(const MmsCase &smooth, Problem problem, const std::string &dataset,
                      bool assert_curl)
{
  // For order-p (first-kind) Nédélec elements, both the L2 field error and the L2 curl
  // error of the curl-curl solution converge as O(h^p) for a smooth solution (Monk,
  // "Finite Element Methods for Maxwell's Equations," 2003). The paper's eq. 29,
  // e ∝ N^(−2/3) for second-order tetrahedra with N unknowns, is the p = 2 case of this,
  // since N ∝ h^(−3) in 3D. The L2 projection's field error is bounded by the best
  // approximation error, also O(h^p) (ND_p contains the complete polynomials of degree
  // p − 1).
  const std::array<std::pair<const char *, mfem::Element::Type>, 2> elements = {
      std::pair{"hexahedron", mfem::Element::HEXAHEDRON},
      std::pair{"tetrahedron", mfem::Element::TETRAHEDRON}};
  for (const auto &[name, type] : elements)
  {
    for (const int order : kOrders)
    {
      DYNAMIC_SECTION(name << ", order " << order)
      {
        const auto kResolutions = Resolutions(order);
        std::vector<double> h(kResolutions.size());
        std::vector<MmsError> error(kResolutions.size());
        for (std::size_t i = 0; i < kResolutions.size(); i++)
        {
          h[i] = 1.0 / kResolutions[i];
          error[i] = SolveCartesianMmsError(smooth, problem, kResolutions[i], order,
                                            /*linear_tol=*/1.0e-9, type);
          ReportMmsData(dataset, "trig", name, order, "N", kResolutions[i], h[i], error[i]);
        }

        for (std::size_t i = 0; i + 1 < kResolutions.size(); i++)
        {
          CAPTURE(kResolutions[i], kResolutions[i + 1], error[i].electric,
                  error[i + 1].electric, error[i].curl, error[i + 1].curl);

          CHECK(error[i + 1].electric < error[i].electric);
          const double electric_rate =
              ConvergenceRate(h[i], error[i].electric, h[i + 1], error[i + 1].electric);
          CAPTURE(electric_rate);
          CHECK_THAT(electric_rate, WithinAbs(static_cast<double>(order), kRateTol));

          if (assert_curl)
          {
            CHECK(error[i + 1].curl < error[i].curl);
            const double curl_rate =
                ConvergenceRate(h[i], error[i].curl, h[i + 1], error[i + 1].curl);
            CAPTURE(curl_rate);
            CHECK_THAT(curl_rate, WithinAbs(static_cast<double>(order), kRateTol));
          }
        }
      }
    }
  }
}

}  // namespace

void CheckMassProjectionExactness(const MmsCase &polynomial)
{
  CheckExactness(polynomial, Problem::MASS_PROJECTION, "mass-affine-polynomial");
}

void CheckMassProjectionConvergence(const MmsCase &smooth)
{
  CheckConvergence(smooth, Problem::MASS_PROJECTION, "mass-cartesian-smooth",
                   /*assert_curl=*/false);
}

void CheckPolynomialExactness(const MmsCase &polynomial)
{
  CheckExactness(polynomial, Problem::CURL_CURL_MASS, "affine-polynomial");
}

void CheckOptimalConvergence(const MmsCase &smooth)
{
  CheckConvergence(smooth, Problem::CURL_CURL_MASS, "cartesian-smooth",
                   /*assert_curl=*/true);
}

void CheckNeumannExactness(const MmsCase &polynomial)
{
  CheckExactness(polynomial, Problem::CURL_CURL_MASS_NEUMANN, "neumann-affine-polynomial");
}

void CheckNeumannConvergence(const MmsCase &smooth)
{
  CheckConvergence(smooth, Problem::CURL_CURL_MASS_NEUMANN, "neumann-cartesian-smooth",
                   /*assert_curl=*/true);
}

void CheckBackendAgreement(const MmsCase &candidate, const MmsCase &reference)
{
  // Compare a backend's field and curl errors against a reference's on the same problem.
  // This exercises the callbacks, including the strong-form curl-curl term, at every
  // quadrature point of the right-hand side and error integrals. The smooth case keeps the
  // errors far above round-off (and its derivatives are not exact under Enzyme, unlike a
  // low-degree polynomial's), so the relative difference is meaningful.
  REQUIRE(candidate.mu_r == reference.mu_r);
  REQUIRE(candidate.epsilon_r == reference.epsilon_r);
  REQUIRE(candidate.omega == reference.omega);
  constexpr int n = 8, order = 2;
  constexpr double linear_tol = 1.0e-13;
  const auto candidate_error =
      SolveCartesianMmsError(candidate, Problem::CURL_CURL_MASS, n, order, linear_tol);
  const auto reference_error =
      SolveCartesianMmsError(reference, Problem::CURL_CURL_MASS, n, order, linear_tol);
  REQUIRE(candidate_error.primary_dofs == reference_error.primary_dofs);

  const double electric_disagreement =
      std::abs(candidate_error.electric - reference_error.electric) /
      reference_error.electric;
  const double curl_disagreement =
      std::abs(candidate_error.curl - reference_error.curl) / reference_error.curl;
  CAPTURE(candidate_error.electric, reference_error.electric, electric_disagreement,
          candidate_error.curl, reference_error.curl, curl_disagreement);
  CHECK(electric_disagreement < kBackendAgreementTol);
  CHECK(curl_disagreement < kBackendAgreementTol);
}

}  // namespace palace::testing::mms::driven
