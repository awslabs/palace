// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification of the driven (frequency-domain)
// operator, Test 3 of Garcia-Doñoro et al., "Verification Process of Finite-Element Method
// Code for Electromagnetics: Using the Method of Manufactured Solutions," IEEE Antennas
// Propag. Mag. 58(2), 2016. Test 3 adds the mass term to the curl-curl operator of Test 2
// (see test-mms-magnetostatic.cpp) with Dirichlet conditions only. Here the material is
// lossless with anisotropic μ_r and ε_r, and the angular frequency ω is real. As in the
// paper (eq. 20), the manufactured field E_mms defines a volume source through the strong
// form of the PDE,
//   J = ∇×(μ⁻¹∇×E_mms) − ω²εE_mms,
// and the code solves
//   ∫ μ⁻¹(∇×E)·(∇×φ) dΩ − ω² ∫ εE·φ dΩ = ∫ J·φ dΩ,   n×E = n×E_mms = 0 on ∂Ω.
// Integrating the curl term by parts recovers the left side for E = E_mms, with a boundary
// term ∮ (n×μ⁻¹∇×E_mms)·φ dS that vanishes because n×φ = 0 on the PEC boundary. The
// right-hand side therefore checks that Palace's discrete operator solves the PDE itself,
// not only that it matches another evaluation of the same bilinear form (the weak-curl
// right-hand side of the magnetostatic test, which is required there because the
// operator is singular).
//
// Unlike the magnetostatic case, the operator is nonsingular when ω² is not a cavity
// eigenvalue, so E is unique and both the field error ‖E_h − E_mms‖ and the curl error
// ‖∇×E_h − ∇×E_mms‖ (absolute versions of the paper's e_f and e_c) are asserted.
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
//   RHS = R_ndᵀ [ ∫ J·φ dΩ ],   R_nd = the ND prolongation (local to true dofs).
// The source integral uses the caller's quadrature rule ir, see SolveMmsError.
ComplexVector SolveManufactured(const IoData &iodata, SpaceOperator &space_op, double omega,
                                mfem::VectorCoefficient &source,
                                const mfem::IntegrationRule &ir)
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

MmsError SolveMmsError(const MmsCase &mms, std::unique_ptr<mfem::Mesh> serial_mesh,
                       int order, double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  // Every face is PEC, so n×E = 0 is imposed on all of ∂Ω. The default frequency-domain
  // linear solver is a sparse direct solver when one is built in, so pin AMS (GMRES with
  // the shifted AMS preconditioner) to solve the same way on every build and keep the
  // finest meshes affordable.
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::DRIVEN;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.mu_r.s = mms.mu_r;
  material.epsilon_r.s = mms.epsilon_r;
  for (int i = 0; i < serial_mesh->bdr_attributes.Size(); i++)
  {
    iodata.boundaries.pec.attributes.push_back(serial_mesh->bdr_attributes[i]);
  }
  iodata.solver.order = order;
  iodata.solver.linear.type = LinearSolver::AMS;
  iodata.solver.linear.tol = linear_tol;
  // Quadrature for the source integral ∫ J·φ dΩ. Polynomial exactness requires this
  // integral to be exact, or Galerkin orthogonality no longer pins E_h = E_mms. For the
  // bubble field J has degree 2 per variable (total degree 4), and φ ∈ ND_p has degree
  // ≤ p per variable on hexahedra and ≤ p in total on tetrahedra, so J·φ needs degree
  // p + 2 per variable (Gauss tensor rule) or p + 4 in total (simplex rule). MFEM's
  // default 2p is exact with no slack for order-5 tetrahedra; 2p + 4 covers both element
  // types at every order with margin and also tightens the trigonometric source.
  // (Geometry from the serial mesh, which every rank holds in full, so ranks that own no
  // elements still get a valid rule.)
  const mfem::IntegrationRule &source_ir =
      mfem::IntRules.Get(serial_mesh->GetElementGeometry(0), 2 * order + 4);

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
  // The volume source J = ∇×(μ⁻¹∇×E_mms) − ω²εE_mms. The backend supplies the curl-curl
  // term, which carries the μ dependence and the second derivatives. Scaling the mass term
  // component-wise equals εE only because the default (identity) material axes make
  // diag(ε_r) diagonal in the global frame (and likewise for μ in the backends); a rotated
  // tensor would need R diag(ε_r) Rᵀ here.
  mfem::VectorFunctionCoefficient source(
      /*vdim=*/3,
      [&mms](const mfem::Vector &x, mfem::Vector &f)
      {
        mms.curl_curl(x, f);
        mfem::Vector e(3);
        mms.electric_field(x, e);
        for (int d = 0; d < 3; d++)
        {
          f(d) -= mms.omega * mms.omega * mms.epsilon_r[d] * e(d);
        }
      });
  const ComplexVector E = SolveManufactured(iodata, space_op, mms.omega, source, source_ir);
  return ComputeMmsErrors(space_op, E, e_exact, curl_exact);
}

MmsError
SolveCartesianMmsError(const MmsCase &mms, int n, int order, double linear_tol = 1.0e-9,
                       mfem::Element::Type element_type = mfem::Element::HEXAHEDRON)
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(n, n, n, element_type, 1.0, 1.0, 1.0));
  return SolveMmsError(mms, std::move(serial_mesh), order, linear_tol);
}

const std::vector<int> kOrders = {1, 2, 3};
// The complex GMRES solve with the shifted AMS preconditioner is several times more
// expensive than the magnetostatic CG solve on the same mesh, and order 3 on a 16³ mesh
// alone takes about a minute per element type, so order 3 stops at 8³. It also starts at
// 4³ like the other orders: a 2³ mesh (8 hexahedra, 48 tetrahedra) is too coarse to
// resolve sin πx, so the 2 → 4 rate is pre-asymptotic, and the discrete eigenvalues near
// ω² are least accurate there (see kOmega). Order 3 therefore has a single rate pair.
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

}  // namespace

void CheckPolynomialExactness(const MmsCase &polynomial)
{
  // The operator is nonsingular, so when E_mms lies in the Nédélec space, E_h = E_mms up
  // to the iterative solver tolerance (Galerkin orthogonality with a unique solution),
  // and both the field and curl errors vanish. Unlike the magnetostatic check, E_mms
  // itself must be in ND_p, not only its curl in ∇×ND_p. The bubble polynomial has degree
  // 2 per variable (total degree 4): order 2 suffices on hexahedra (ND_2 ⊇ Q_{1,2,2} ×
  // Q_{2,1,2} × Q_{2,2,1}), while on tetrahedra ND_p contains the full P_{p-1}, so order 5
  // is needed. A coarse n = 2 tetrahedral mesh keeps order 5 cheap.
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
      const auto error = SolveCartesianMmsError(polynomial, n, order,
                                                /*linear_tol=*/1.0e-13, type);
      ReportMmsData("affine-polynomial", "polynomial", name, order, "N", n, 1.0 / n, error);
      CAPTURE(error.electric, error.curl);
      CHECK(error.electric < 1.0e-10);
      CHECK(error.curl < 1.0e-10);
    }
  }
}

void CheckOptimalConvergence(const MmsCase &smooth)
{
  // For order-p (first-kind) Nédélec elements, both the L2 field error and the L2 curl
  // error converge as O(h^p) for a smooth solution (Monk, "Finite Element Methods for
  // Maxwell's Equations," 2003). The paper's eq. 29, e ∝ N^(−2/3) for second-order
  // tetrahedra with N unknowns, is the p = 2 case of this, since N ∝ h^(−3) in 3D.
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
          error[i] = SolveCartesianMmsError(smooth, kResolutions[i], order,
                                            /*linear_tol=*/1.0e-9, type);
          ReportMmsData("cartesian-smooth", "trig", name, order, "N", kResolutions[i], h[i],
                        error[i]);
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
  const auto candidate_error = SolveCartesianMmsError(candidate, n, order, linear_tol);
  const auto reference_error = SolveCartesianMmsError(reference, n, order, linear_tol);
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
