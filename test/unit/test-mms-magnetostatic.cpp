// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Method of Manufactured Solutions (MMS) verification of the magnetostatic operator, Test 2
// of Garcia-Doñoro et al., "Verification Process of Finite-Element Method Code for
// Electromagnetics: Using the Method of Manufactured Solutions," IEEE Antennas Propag. Mag.
// 58(2), 2016. The magnetostatic CurlCurlOperator assembles the curl-curl bilinear form
// with no mass term. The paper's eq. 19 is the μ = 1 case; here the relative permeability
// is anisotropic, and the checks run on both hexahedral and tetrahedral meshes:
//   ∫ μ⁻¹ (∇×A)·(∇×φ) dΩ = ∫ μ⁻¹ (∇×A_mms)·(∇×φ) dΩ,   n×A = n×A_mms = 0 on ∂Ω,
// so the test checks the derivatives of the Nédélec basis functions in isolation.
//
// The operator is singular: adding the gradient of any H1 function vanishing on ∂Ω to A
// leaves B = ∇×A unchanged. The vector potential A is therefore only determined up to a
// discrete gradient, and the assertions use the absolute flux density error ‖B_h − B_mms‖
// (an absolute version of the paper's relative curl error e_c). The potential error is
// reported but not asserted; it does not converge under refinement.
//
// Backend-specific tests supply A_mms and B_mms through MmsCase; this file owns the common
// Palace solve and error checks.

#include "test-mms-magnetostatic.hpp"
#include "test-mms-common.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/ksp.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace palace::testing::mms::magnetostatic
{
namespace
{

struct MmsError
{
  HYPRE_BigInt primary_dofs;
  double potential;  // ‖A_h - A_mms‖, defined only up to a discrete gradient: report only.
  double magnetic;   // ‖B_h - B_mms‖ with B = ∇×A, the asserted quantity.
};

// Solve the curl-curl system with the production stiffness matrix, ND/H1 space hierarchies,
// and KspSolver (CG + singular AMS, the magnetostatic defaults), assembling only the
// manufactured right-hand side here. It uses the weak-curl form of the paper's eq. 19,
//   RHS = P^T [ ∫ μ⁻¹(∇×A_mms)·(∇×φ) dΩ ],
// rather than the strong form ∫ J·φ with J = ∇×μ⁻¹∇×A_mms. The weak-curl vector is Cᵀ(...)
// for the discrete curl C, and C G = 0 for the discrete gradient G (the de Rham exact
// sequence), so Gᵀ RHS = 0 exactly: the right-hand side is consistent with the singular
// operator regardless of quadrature error. The strong form's quadrature error leaves a
// small component in the gradient null space, which the singular solve cannot remove. The
// weak-curl form also needs no second derivatives of A_mms. The Dirichlet data is zero, so
// EliminateRHS only zeroes the essential dofs.
Vector SolveManufactured(const IoData &iodata, CurlCurlOperator &curlcurl_op,
                         mfem::VectorCoefficient &h_exact)
{
  auto K = curlcurl_op.GetStiffnessMatrix();
  KspSolver ksp(iodata, curlcurl_op.GetNDSpaces(), &curlcurl_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  auto &nd = curlcurl_op.GetNDSpace();
  mfem::ParLinearForm b(&nd.Get());
  b.AddDomainIntegrator(new mfem::VectorFEDomainLFCurlIntegrator(h_exact));
  b.Assemble();

  Vector X(nd.GetTrueVSize()), RHS(nd.GetTrueVSize());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  nd.GetProlongationMatrix()->MultTranspose(b, RHS);

  SolveEliminated(ksp, *K, X, RHS);
  return X;
}

MmsError ComputeMmsErrors(CurlCurlOperator &curlcurl_op, const Vector &A,
                          mfem::VectorCoefficient &a_exact,
                          mfem::VectorCoefficient &b_exact)
{
  mfem::ParGridFunction A_gf(&curlcurl_op.GetNDSpace().Get());
  A_gf.SetFromTrueDofs(A);

  // B = ∇×A_h, mapped exactly from the ND space into the RT space by the discrete curl.
  const auto &Curl = curlcurl_op.GetCurlMatrix();
  Vector B(Curl.Height());
  B.UseDevice(true);
  B = 0.0;
  Curl.Mult(A, B);
  mfem::ParGridFunction B_gf(&curlcurl_op.GetCurlSpace().Get());
  B_gf.SetFromTrueDofs(B);

  return {curlcurl_op.GlobalTrueVSize(), A_gf.ComputeL2Error(a_exact),
          B_gf.ComputeL2Error(b_exact)};
}

MmsError SolveMmsError(const MmsCase &mms, std::unique_ptr<mfem::Mesh> serial_mesh,
                       int order, double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  // MAGNETOSTATIC selects the curl-curl defaults in CheckConfiguration (CG and
  // singular-operator AMS, since the eliminated operator is only positive semi-definite).
  // Every face is PEC, so n×A = 0 is imposed on all of ∂Ω.
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::MAGNETOSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.mu_r.s = mms.mu_r;
  for (int i = 0; i < serial_mesh->bdr_attributes.Size(); i++)
  {
    iodata.boundaries.pec.attributes.push_back(serial_mesh->bdr_attributes[i]);
  }
  iodata.solver.order = order;
  iodata.solver.linear.tol = linear_tol;

  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  CurlCurlOperator curlcurl_op(iodata, mesh);
  mfem::VectorFunctionCoefficient a_exact(/*vdim=*/3, mms.potential);
  mfem::VectorFunctionCoefficient b_exact(/*vdim=*/3, mms.magnetic_flux_density);
  // The right-hand side data H_mms = μ⁻¹ B_mms. Dividing component-wise is μ⁻¹ only
  // because the material uses the default (identity) material axes, so the tensor
  // diag(mu_r) is diagonal in the global frame. A rotated-tensor case (non-default
  // material_axes) would instead need the full inverse R diag(1/mu_r) Rᵀ here.
  mfem::VectorFunctionCoefficient h_exact(
      /*vdim=*/3,
      [&mms](const mfem::Vector &x, mfem::Vector &f)
      {
        mms.magnetic_flux_density(x, f);
        for (int d = 0; d < 3; d++)
        {
          f(d) /= mms.mu_r[d];
        }
      });
  const Vector A = SolveManufactured(iodata, curlcurl_op, h_exact);
  return ComputeMmsErrors(curlcurl_op, A, a_exact, b_exact);
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
const std::vector<int> kResolutions = {4, 8, 16};

void ReportMmsData(const std::string &dataset, const std::string &solution,
                   const std::string &element, int order,
                   const std::string &resolution_kind, int resolution, double h,
                   const MmsError &error)
{
  mms::ReportMmsData(dataset, solution, element, order, resolution_kind, resolution, h,
                     error.primary_dofs, error.potential, error.magnetic);
}

constexpr double kRateTol = 0.3;
constexpr double kBackendAgreementTol = 1.0e-10;

}  // namespace

void CheckPolynomialExactness(const MmsCase &polynomial)
{
  // The weak-curl Galerkin problem makes B_h = ∇×A_h the μ⁻¹-weighted L2 projection of
  // B_mms onto ∇×ND_p (with n×A_h = 0). B_h is therefore exact, up to the iterative solver
  // tolerance, whenever B_mms lies in ∇×ND_p, even if A_mms itself does not lie in ND_p:
  // the gradient part of A is invisible to B. On a contractible domain the discrete de Rham
  // sequence makes ∇×ND_p the divergence-free fields of RT_{p-1} with B·n = 0 on ∂Ω, which
  // B_mms satisfies, so only the polynomial degree matters:
  //  - hexahedra: B_x ∈ Q_{2,1,1} (similarly for y, z) ⊂ RT_1, so order 2 is exact;
  //  - tetrahedra: B has total degree 3, and RT_{p-1} contains the full P_{p-1}, so order 4
  //    is exact (order 3 gives an error of about 4e-3 on this mesh).
  // Each element uses the lowest exact order, so a mismatch at that order is not hidden by
  // spare polynomial degree. A coarse n = 2 tetrahedral mesh keeps order 4 cheap.
  struct Element
  {
    const char *name;
    mfem::Element::Type type;
    int order, n;
  };
  const std::array<Element, 2> elements = {
      Element{"hexahedron", mfem::Element::HEXAHEDRON, /*order=*/2, /*n=*/4},
      Element{"tetrahedron", mfem::Element::TETRAHEDRON, /*order=*/4, /*n=*/2}};
  for (const auto &[name, type, order, n] : elements)
  {
    DYNAMIC_SECTION(name)
    {
      const auto error = SolveCartesianMmsError(polynomial, n, order,
                                                /*linear_tol=*/1.0e-13, type);
      ReportMmsData("affine-polynomial", "polynomial", name, order, "N", n, 1.0 / n, error);
      CAPTURE(error.potential, error.magnetic);
      CHECK(error.magnetic < 1.0e-10);
    }
  }
}

void CheckOptimalConvergence(const MmsCase &smooth)
{
  // For order-p Nédélec elements the curl lies in the order-(p - 1) Raviart-Thomas space,
  // so the L2 flux density error converges as O(h^p) (see, e.g., Monk, "Finite Element
  // Methods for Maxwell's Equations," 2003). This holds for both element shapes: MFEM's
  // ND_p space contains the full degree-(p - 1) polynomials on tetrahedra (first-kind
  // Nédélec) and Q_{p-1,p,p} × ... on hexahedra, and in both cases ∇×ND_p = RT_{p-1}
  // ∩ ker(div), which approximates the smooth B_mms to O(h^p). MakeCartesian3D splits each
  // hexahedron into six tetrahedra, so h is the same cube edge length for both meshes.
  const std::array<std::pair<const char *, mfem::Element::Type>, 2> elements = {
      std::pair{"hexahedron", mfem::Element::HEXAHEDRON},
      std::pair{"tetrahedron", mfem::Element::TETRAHEDRON}};
  for (const auto &[name, type] : elements)
  {
    for (const int order : kOrders)
    {
      DYNAMIC_SECTION(name << ", order " << order)
      {
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
          CAPTURE(kResolutions[i], kResolutions[i + 1], error[i].magnetic,
                  error[i + 1].magnetic);
          CHECK(error[i + 1].magnetic < error[i].magnetic);
          const double magnetic_rate =
              ConvergenceRate(h[i], error[i].magnetic, h[i + 1], error[i + 1].magnetic);
          CAPTURE(magnetic_rate);
          CHECK_THAT(magnetic_rate, WithinAbs(static_cast<double>(order), kRateTol));
        }
      }
    }
  }
}

void CheckBackendAgreement(const MmsCase &candidate, const MmsCase &reference)
{
  // Compare a backend's flux density error against a reference's on the same problem.
  // Unlike a pointwise comparison, this exercises the callbacks at every quadrature point
  // of the right-hand side and error integrals. The smooth trigonometric case keeps the
  // error far above round-off, so the relative difference is meaningful. The potential
  // error is not compared because A_h is only defined up to a discrete gradient.
  REQUIRE(candidate.mu_r == reference.mu_r);
  constexpr int n = 8, order = 2;
  constexpr double linear_tol = 1.0e-13;
  const auto candidate_error = SolveCartesianMmsError(candidate, n, order, linear_tol);
  const auto reference_error = SolveCartesianMmsError(reference, n, order, linear_tol);
  REQUIRE(candidate_error.primary_dofs == reference_error.primary_dofs);

  const double magnetic_disagreement =
      std::abs(candidate_error.magnetic - reference_error.magnetic) /
      reference_error.magnetic;
  CAPTURE(candidate_error.magnetic, reference_error.magnetic, magnetic_disagreement);
  CHECK(magnetic_disagreement < kBackendAgreementTol);
}

}  // namespace palace::testing::mms::magnetostatic
