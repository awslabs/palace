// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Shared Method of Manufactured Solutions (MMS) verification for the electrostatic
// operator. Backend-specific tests supply V_mms, E_mms, rho_mms, and optional Neumann flux
// callbacks through MmsCase; this file owns the common Palace solve and error checks.

#include "test-mms-electrostatic.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <utility>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;

namespace palace::testing::mms::electrostatic
{
namespace
{

struct MmsError
{
  HYPRE_BigInt primary_dofs;
  double potential;
  double electric;
};

MmsError ComputeMmsErrors(LaplaceOperator &laplace_op, const Vector &V,
                          mfem::Coefficient &v_exact, mfem::VectorCoefficient &e_exact)
{
  mfem::ParGridFunction V_gf(&laplace_op.GetH1Space().Get());
  V_gf.SetFromTrueDofs(V);

  Vector E(laplace_op.GetGradMatrix().Height());
  E = 0.0;
  laplace_op.GetGradMatrix().AddMult(V, E, -1.0);
  mfem::ParGridFunction E_gf(&laplace_op.GetNDSpace().Get());
  E_gf.SetFromTrueDofs(E);

  return {laplace_op.GlobalTrueVSize(), V_gf.ComputeL2Error(v_exact),
          E_gf.ComputeL2Error(e_exact)};
}

// Solve -div(epsilon grad V) = rho with the production stiffness matrix, H1 space
// hierarchy, and KspSolver, assembling only the manufactured data here. The production
// driver has no volumetric source (its RHS comes from unit terminal voltages), so the
// right-hand side mirrors LaplaceOperator::GetExcitationVector with the manufactured
// quantities in place of the terminal excitation:
//   RHS = P^T [ int rho phi dV + int_{Gamma_N} g phi dS ],
//   X   = Dirichlet lift of V_exact on the PEC attributes (or 0 when homogeneous),
// followed by EliminateRHS, which folds -K_ib X_b into RHS and pins the essential dofs.
Vector SolveManufactured(const IoData &iodata, LaplaceOperator &laplace_op,
                         mfem::Coefficient &rho, mfem::Coefficient *dirichlet = nullptr,
                         mfem::Coefficient *neumann = nullptr,
                         const std::vector<int> &neumann_attributes = {})
{
  auto K = laplace_op.GetStiffnessMatrix();
  KspSolver ksp(iodata, laplace_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  auto &h1 = laplace_op.GetH1Space();
  const mfem::ParMesh &mesh = laplace_op.GetMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;

  // Volumetric load. mfem::DomainLFIntegrator is the volume integrator; the similarly
  // named palace::DomainLFIntegrator integrates over boundaries. The Neumann marker must
  // outlive Assemble(), which holds a reference to it.
  mfem::ParLinearForm b(&h1.Get());
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(rho));
  mfem::Array<int> neumann_marker;
  if (neumann)
  {
    neumann_marker = mesh::AttrToMarker(bdr_attr_max, neumann_attributes);
    b.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*neumann), neumann_marker);
  }
  b.Assemble();

  Vector X(h1.GetTrueVSize()), RHS(h1.GetTrueVSize());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  h1.GetProlongationMatrix()->MultTranspose(b, RHS);

  // Non-homogeneous Dirichlet data on every PEC attribute (the essential set when there
  // are no terminals). Values are only correct on master dofs until ParallelProject.
  if (dirichlet)
  {
    const auto dbc_marker =
        mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);
    mfem::ParGridFunction x(&h1.Get());
    x = 0.0;
    x.ProjectBdrCoefficient(*dirichlet, dbc_marker);
    x.ParallelProject(X);
  }

  // Same unwrapping as LaplaceOperator::GetExcitationVector: with multigrid the stiffness
  // matrix is a MultigridOperator whose finest level holds the eliminated ParOperator.
  const auto *mg_K = dynamic_cast<const MultigridOperator *>(K.get());
  const auto *PtAP_K = mg_K ? dynamic_cast<const ParOperator *>(&mg_K->GetFinestOperator())
                            : dynamic_cast<const ParOperator *>(K.get());
  REQUIRE(PtAP_K);
  PtAP_K->EliminateRHS(X, RHS);

  ksp.Mult(RHS, X);
  return X;
}

MmsError SolveMmsError(const MmsCase &mms, std::unique_ptr<mfem::Mesh> serial_mesh,
                       int order, double linear_tol = 1.0e-9)
{
  MPI_Comm comm = Mpi::World();

  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = mms.epsilon_r;
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

  LaplaceOperator laplace_op(iodata, mesh);
  mfem::FunctionCoefficient v_exact(mms.potential);
  mfem::VectorFunctionCoefficient e_exact(/*vdim=*/3, mms.electric_field);
  mfem::FunctionCoefficient rho_coef(mms.charge_density);
  const Vector V = SolveManufactured(iodata, laplace_op, rho_coef,
                                     mms.nonzero_dirichlet ? &v_exact : nullptr);
  return ComputeMmsErrors(laplace_op, V, v_exact, e_exact);
}

MmsError
SolveCartesianMmsError(const MmsCase &mms, int n, int order, double linear_tol = 1.0e-9,
                       mfem::Element::Type element_type = mfem::Element::HEXAHEDRON)
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(n, n, n, element_type, 1.0, 1.0, 1.0));
  return SolveMmsError(mms, std::move(serial_mesh), order, linear_tol);
}

std::unique_ptr<mfem::Mesh> LoadCurvedCylinderMesh(int refinement)
{
  const std::string mesh_path =
      std::string(PALACE_TEST_DATA_DIR) + "/mesh/cylinder-tet-p2.msh";
  std::ifstream mesh_file(mesh_path);
  REQUIRE(mesh_file.good());
  constexpr bool generate_edges = false, refine = true, fix_orientation = true;
  auto serial_mesh =
      std::make_unique<mfem::Mesh>(mesh_file, generate_edges, refine, fix_orientation);
  REQUIRE(serial_mesh->GetNodes());
  REQUIRE(serial_mesh->GetNodes()->FESpace()->GetMaxElementOrder() == 2);
  for (int l = 0; l < refinement; l++)
  {
    serial_mesh->UniformRefinement();
  }
  return serial_mesh;
}

double ConvergenceRate(double h0, double e0, double h1, double e1)
{
  return std::log(e0 / e1) / std::log(h0 / h1);
}

const std::vector<int> kOrders = {1, 2, 3};
const std::vector<int> kResolutions = {4, 8, 16};
const std::vector<int> kCurvedRefinements = {0, 1, 2};

void ReportMmsData(const std::string &dataset, const std::string &solution,
                   const std::string &element, int order,
                   const std::string &resolution_kind, int resolution, double h,
                   const MmsError &error)
{
  if (std::getenv("PALACE_MMS_REPORT"))
  {
    Mpi::Print("MMS_DATA,{},{},{},{},{},{},{:.17g},{},{:.17g},{:.17g}\n", dataset, solution,
               element, order, resolution_kind, resolution, h, error.primary_dofs,
               error.potential, error.electric);
  }
}

constexpr double kRateTol = 0.3;
constexpr double kCurvedRateTol = 0.35;

}  // namespace

void CheckSmallErrors(const MmsCase &homogeneous, const MmsCase &nonhomogeneous)
{
  const auto homogeneous_error = SolveCartesianMmsError(homogeneous, /*n=*/16, /*order=*/2);
  CHECK(homogeneous_error.potential < 1.0e-3);
  CHECK(homogeneous_error.electric < 2.0e-2);

  const auto nonhomogeneous_error =
      SolveCartesianMmsError(nonhomogeneous, /*n=*/16, /*order=*/2);
  CHECK(nonhomogeneous_error.potential < 1.0e-3);
  CHECK(nonhomogeneous_error.electric < 2.0e-2);
}

void CheckPolynomialExactness(const MmsCase &polynomial)
{
  // A degree-2 potential lies exactly in the order-2 H1 space on affine meshes, so only
  // iterative-solver error remains. Exercise both tensor-product and simplex elements.
  const std::array<std::pair<const char *, mfem::Element::Type>, 2> elements = {
      std::pair{"hexahedron", mfem::Element::HEXAHEDRON},
      std::pair{"tetrahedron", mfem::Element::TETRAHEDRON}};
  for (const auto &[name, type] : elements)
  {
    DYNAMIC_SECTION(name)
    {
      const auto error = SolveCartesianMmsError(polynomial, /*n=*/4, /*order=*/2,
                                                /*linear_tol=*/1.0e-13, type);
      ReportMmsData("affine-polynomial", "polynomial", name, /*order=*/2, "N",
                    /*resolution=*/4, /*h=*/0.25, error);
      CHECK(error.potential < 1.0e-10);
      CHECK(error.electric < 1.0e-10);
    }
  }
}

void CheckNeumannBoundary(const MmsCase &polynomial)
{
  // Leave the x = 1 face out of the Dirichlet set and prescribe n·epsilon grad(V) there.
  // The polynomial is FE-exact, so an incorrectly signed or omitted flux produces a large
  // error instead of the near-machine-precision result required below.
  REQUIRE(polynomial.neumann_flux);

  MPI_Comm comm = Mpi::World();
  Units units(1.0, 1.0);
  IoData iodata(units);
  iodata.model.Lc = 1.0;
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  auto &material = iodata.domains.materials.emplace_back();
  material.attributes = {1};
  material.epsilon_r.s = polynomial.epsilon_r;
  iodata.boundaries.pec.attributes = {1, 2, 4, 5, 6};
  iodata.solver.order = 2;
  iodata.solver.linear.tol = 1.0e-13;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0));
  iodata.NondimensionalizeInputs(serial_mesh);
  iodata.CheckConfiguration();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  LaplaceOperator laplace_op(iodata, mesh);
  mfem::FunctionCoefficient v_exact(polynomial.potential);
  mfem::VectorFunctionCoefficient e_exact(/*vdim=*/3, polynomial.electric_field);
  mfem::FunctionCoefficient rho_coef(polynomial.charge_density);
  mfem::FunctionCoefficient neumann_coef(polynomial.neumann_flux);
  const Vector V = SolveManufactured(iodata, laplace_op, rho_coef, &v_exact, &neumann_coef,
                                     /*neumann_attributes=*/{3});

  const auto error = ComputeMmsErrors(laplace_op, V, v_exact, e_exact);
  CHECK(error.potential < 1.0e-10);
  CHECK(error.electric < 1.0e-10);
}

void CheckOptimalConvergence(const MmsCase &nonhomogeneous)
{
  // The non-homogeneous solution exercises the Dirichlet lift in addition to the interior
  // operator. Check every consecutive refinement pair for the expected p+1 potential rate
  // and p electric-field rate.
  for (const int order : kOrders)
  {
    DYNAMIC_SECTION("non-homogeneous, order " << order)
    {
      std::vector<double> h(kResolutions.size());
      std::vector<MmsError> error(kResolutions.size());
      for (std::size_t i = 0; i < kResolutions.size(); i++)
      {
        h[i] = 1.0 / kResolutions[i];
        error[i] = SolveCartesianMmsError(nonhomogeneous, kResolutions[i], order);
        ReportMmsData("cartesian-smooth", "non-homogeneous", "hexahedron", order, "N",
                      kResolutions[i], h[i], error[i]);
      }

      for (std::size_t i = 0; i + 1 < kResolutions.size(); i++)
      {
        CAPTURE(kResolutions[i], kResolutions[i + 1], error[i].potential,
                error[i + 1].potential, error[i].electric, error[i + 1].electric);

        CHECK(error[i + 1].potential < error[i].potential);
        const double potential_rate =
            ConvergenceRate(h[i], error[i].potential, h[i + 1], error[i + 1].potential);
        CAPTURE(potential_rate);
        CHECK_THAT(potential_rate, WithinAbs(order + 1.0, kRateTol));

        CHECK(error[i + 1].electric < error[i].electric);
        const double electric_rate =
            ConvergenceRate(h[i], error[i].electric, h[i + 1], error[i + 1].electric);
        CAPTURE(electric_rate);
        CHECK_THAT(electric_rate, WithinAbs(static_cast<double>(order), kRateTol));
      }
    }
  }
}

void CheckCurvedConvergence(const MmsCase &curved)
{
  // Quadratic isoparametric tetrahedra make the physical-space basis non-polynomial even
  // for the polynomial manufactured potential. Refining the curved map should recover the
  // optimal rates, matching the curved-geometry experiment in the reference paper.
  constexpr int order = 2;
  std::vector<MmsError> error(kCurvedRefinements.size());
  for (std::size_t i = 0; i < kCurvedRefinements.size(); i++)
  {
    error[i] = SolveMmsError(curved, LoadCurvedCylinderMesh(kCurvedRefinements[i]), order);
    ReportMmsData("curved-polynomial", "polynomial", "curved-tetrahedron", order,
                  "refinement", kCurvedRefinements[i], std::pow(0.5, kCurvedRefinements[i]),
                  error[i]);
  }

  for (std::size_t i = 0; i + 1 < kCurvedRefinements.size(); i++)
  {
    CAPTURE(kCurvedRefinements[i], kCurvedRefinements[i + 1], error[i].potential,
            error[i + 1].potential, error[i].electric, error[i + 1].electric);

    CHECK(error[i + 1].potential < error[i].potential);
    const double potential_rate =
        ConvergenceRate(1.0, error[i].potential, 0.5, error[i + 1].potential);
    CAPTURE(potential_rate);
    CHECK_THAT(potential_rate, WithinAbs(order + 1.0, kCurvedRateTol));

    CHECK(error[i + 1].electric < error[i].electric);
    const double electric_rate =
        ConvergenceRate(1.0, error[i].electric, 0.5, error[i + 1].electric);
    CAPTURE(electric_rate);
    CHECK_THAT(electric_rate, WithinAbs(static_cast<double>(order), kCurvedRateTol));
  }
}

}  // namespace palace::testing::mms::electrostatic
