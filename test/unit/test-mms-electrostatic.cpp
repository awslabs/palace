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
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "mms-solvers.hpp"
#include "models/laplaceoperator.hpp"
#include "utils/communication.hpp"
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
  laplace_op.SetRhsSource(rho_coef);
  if (mms.nonzero_dirichlet)
  {
    laplace_op.SetDbcCoefficient(v_exact);
  }

  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  std::vector<Vector> V;
  solver.Solve(V, laplace_op);
  return ComputeMmsErrors(laplace_op, V[0], v_exact, e_exact);
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

void CheckPointwiseAgreement(const MmsCase &enzyme, const MmsCase &analytic,
                             const std::vector<std::array<double, 3>> &points,
                             bool check_neumann)
{
  REQUIRE(enzyme.epsilon_r == analytic.epsilon_r);
  if (check_neumann)
  {
    REQUIRE(enzyme.neumann_flux);
    REQUIRE(analytic.neumann_flux);
  }

  for (const auto &point : points)
  {
    mfem::Vector x(3);
    for (int i = 0; i < 3; i++)
    {
      x[i] = point[i];
    }
    mfem::Vector enzyme_electric(3);
    mfem::Vector analytic_electric(3);
    enzyme.electric_field(x, enzyme_electric);
    analytic.electric_field(x, analytic_electric);

    CHECK_THAT(enzyme.potential(x), WithinAbs(analytic.potential(x), 1.0e-14));
    for (int i = 0; i < 3; i++)
    {
      CHECK_THAT(enzyme_electric[i], WithinAbs(analytic_electric[i], 5.0e-13));
    }
    CHECK_THAT(enzyme.charge_density(x), WithinAbs(analytic.charge_density(x), 5.0e-12));
    if (check_neumann)
    {
      CHECK_THAT(enzyme.neumann_flux(x), WithinAbs(analytic.neumann_flux(x), 5.0e-13));
    }
  }
}

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
  laplace_op.SetRhsSource(rho_coef);
  laplace_op.SetDbcCoefficient(v_exact);
  laplace_op.SetNeumannCoefficient(neumann_coef, {3});

  ExposedElectrostaticSolver solver(iodata, /*root=*/false);
  std::vector<Vector> V;
  solver.Solve(V, laplace_op);

  const auto error = ComputeMmsErrors(laplace_op, V[0], v_exact, e_exact);
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
