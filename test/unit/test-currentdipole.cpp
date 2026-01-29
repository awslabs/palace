// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Verification test for Electrical Current Dipole from:
// https://em.geosci.xyz/content/maxwell1_fundamentals/dipole_sources_in_homogeneous_media/electric_dipole_frequency/analytic_solution.html
//
// For a time-harmonic electrical current dipole in the x-direction (p=Ids·x̂):
//
// Electric field (Eq. 206):
// Ex = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [x² * (-k²r² + 3ikr + 3)/r² + (k²r² - ikr - 1)]
// Ey = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [xy * (-k²r² + 3ikr + 3)/r²]
// Ez = Ids/(4π(σ+iωε)r³) * e^(-ikr) * [xz * (-k²r² + 3ikr + 3)/r²]
//
// Magnetic flux density (B = μ₀H, from Eq. 207):
// Bx = 0
// By = μ₀Ids/(4πr²) * (ikr + 1) * e^(-ikr) * (-z/r)
// Bz = μ₀Ids/(4πr²) * (ikr + 1) * e^(-ikr) * (y/r)
//
// where k² = ω²με - iωμσ, r = √(x² + y² + z²), and p is the dipole moment (A.m).

#include <complex>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mfem/fem/coefficient.hpp>
#include <mfem/linalg/vector.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

// Class for marking elements in the "good" region for error computation
// Excluding center and boundary elements
// Uses centroid-based marking: an element is excluded if its centroid is
// within inner_radius of origin or beyond outer_radius from origin.
// This is more accurate than vertex-based marking which can be overly conservative.
class ElementMarker
{
private:
   mfem::Mesh *mesh;
   int dim;

   // Bounding box of the domain
   mfem::Vector bbox_min, bbox_max;

   // Array marking elements: 1 = include in error computation, 0 = exclude
   mfem::Array<int> elems;

   // Exclusion criteria (fractions of max distance from origin)
   double inner_fraction;
   double outer_fraction;

   // Compute domain bounding box
   void ComputeBoundingBox();

public:
   // Constructor
   ElementMarker(mfem::Mesh *mesh_, double inner_frac = 0.2, double outer_frac = 0.8);

   // Mark elements in the good region
   void MarkElements();

   // Return markers list for elements
   mfem::Array<int> * GetMarkedElements() {return &elems;}
};

ElementMarker::ElementMarker(mfem::Mesh *mesh_, double inner_frac, double outer_frac)
   : mesh(mesh_), inner_fraction(inner_frac), outer_fraction(outer_frac)
{
   dim = mesh->Dimension();
   bbox_min.SetSize(dim);
   bbox_max.SetSize(dim);
   ComputeBoundingBox();
   MarkElements();
}

void ElementMarker::ComputeBoundingBox()
{
   mfem::Vector pmin, pmax;
   mesh->GetBoundingBox(pmin, pmax);
   bbox_min = pmin;
   bbox_max = pmax;
}

void ElementMarker::MarkElements()
{
   int nelem = mesh->GetNE();
   elems.SetSize(nelem);

   // Compute domain characteristic length in each direction
   mfem::Vector domain_size(dim);
   for (int i = 0; i < dim; i++)
   {
      domain_size[i] = bbox_max[i] - bbox_min[i];
   }

   // Compute the maximum distance from origin to any vertex in the mesh
   // For a sphere centered at origin, this gives us the radius
   double max_dist_from_origin = 0.0;
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      double *coords = mesh->GetVertex(i);
      double dist = 0.0;
      for (int d = 0; d < dim; d++)
      {
         dist += coords[d] * coords[d];
      }
      dist = std::sqrt(dist);
      max_dist_from_origin = std::max(max_dist_from_origin, dist);
   }

   // Use maximum distance from origin for spherical exclusion
   double inner_radius = inner_fraction * max_dist_from_origin;
   double outer_radius = outer_fraction * max_dist_from_origin;

   std::cout << "\nElement marking for error computation:";
   std::cout << "\n  Domain bounding box: [" << bbox_min[0] << ", " << bbox_max[0] << "] x ";
   std::cout << "[" << bbox_min[1] << ", " << bbox_max[1] << "] x ";
   std::cout << "[" << bbox_min[2] << ", " << bbox_max[2] << "]";
   std::cout << "\n  Domain size: " << domain_size[0] << " x " << domain_size[1] << " x " << domain_size[2];
   std::cout << "\n  Max distance from origin (sphere radius): " << max_dist_from_origin;
   std::cout << "\n  Inner exclusion fraction: " << inner_fraction << " -> radius: " << inner_radius;
   std::cout << "\n  Outer exclusion fraction: " << outer_fraction << " -> radius: " << outer_radius;

   int included_elements = 0;
   int excluded_center = 0;
   int excluded_boundary = 0;

   // Loop through elements and mark them based on centroid position
   for (int i = 0; i < nelem; i++)
   {
      mfem::Element *el = mesh->GetElement(i);
      mfem::Array<int> vertices;
      el->GetVertices(vertices);
      int nvert = vertices.Size();

      // Compute element centroid
      mfem::Vector centroid(dim);
      centroid = 0.0;
      for (int iv = 0; iv < nvert; iv++)
      {
         int vert_idx = vertices[iv];
         double *coords = mesh->GetVertex(vert_idx);
         for (int d = 0; d < dim; d++)
         {
            centroid[d] += coords[d];
         }
      }
      centroid /= nvert;

      // Calculate distance of centroid from origin
      double dist_from_origin = centroid.Norml2();

      // Mark element based on centroid position
      if (dist_from_origin < inner_radius)
      {
         elems[i] = 0;  // Exclude - too close to center
         excluded_center++;
      }
      else if (dist_from_origin > outer_radius)
      {
         elems[i] = 0;  // Exclude - too close to boundary
         excluded_boundary++;
      }
      else
      {
         elems[i] = 1;  // Include - in good region
         included_elements++;
      }
   }

   std::cout << "\n  Elements included in error computation: " << included_elements;
   std::cout << "\n  Elements excluded (too close to center): " << excluded_center;
   std::cout << "\n  Elements excluded (too close to boundary): " << excluded_boundary;
   std::cout << "\n  Total elements: " << nelem << "\n";
}

namespace palace
{
using namespace Catch;
using namespace electromagnetics;
using namespace Catch::Matchers;
using namespace std::complex_literals;

namespace
{

// Compute the Cartesian components electric field of a time-harmonic electrical
// current dipole aligned in the x-direction (Ids·δ(x)δ(y)δ(z)·x̂) at a given point in space.
// The returned field is non-dimensionalized according to the provided units.
//
// We want the non-dimensional E field to mock a Palace-computed field.
std::array<std::complex<double>, 3>
ComputeCurrentDipoleENonDim(const mfem::Vector &x_nondim, const Units &units, double Ids,
                            double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);
  double x = units.Dimensionalize<Units::ValueType::LENGTH>(x_nondim(0));
  double y = units.Dimensionalize<Units::ValueType::LENGTH>(x_nondim(1));
  double z = units.Dimensionalize<Units::ValueType::LENGTH>(x_nondim(2));

  double omega = 2.0 * M_PI * freq_Hz;
  double k = omega / c0_;
  double kr = k * r;

  std::complex<double> ikr(0, kr);
  std::complex<double> exp_ikr = std::exp(-ikr);
  std::complex<double> factor1 =
      exp_ikr * Ids / (4.0 * M_PI * 1i * omega * epsilon0_ * r * r * r);
  factor1 = units.Nondimensionalize<Units::ValueType::FIELD_E>(factor1);

  std::complex<double> factor2 = (-kr * kr + 3. * ikr + 3.) / (r * r);

  return {{factor1 * (x * x * factor2 + kr * kr - ikr - 1.), factor1 * (x * y * factor2),
           factor1 * (x * z * factor2)}};
}

// Compute the Cartesian components magnetic flux density (B-field) of a time-harmonic electrical
// current dipole aligned in the x-direction (Ids·δ(x)δ(y)δ(z)·x̂) at a given point in space.
// The returned field is non-dimensionalized according to the provided units.
//
// We want the non-dimensional B field to mock a Palace-computed field.
// This is simply B = μ₀ * H, where H is from the reference formulas.
std::array<std::complex<double>, 3>
ComputeCurrentDipoleBNonDim(const mfem::Vector &x_nondim, const Units &units, double Ids,
                            double freq_Hz)
{
  double r_nondim = x_nondim.Norml2();
  if (r_nondim < 1e-12)
    return {0.0, 0.0, 0.0};

  double r = units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);
  double y = units.Dimensionalize<Units::ValueType::LENGTH>(x_nondim(1));
  double z = units.Dimensionalize<Units::ValueType::LENGTH>(x_nondim(2));

  double omega = 2.0 * M_PI * freq_Hz;
  double k = omega / c0_;
  double kr = k * r;

  std::complex<double> ikr(0, kr);
  std::complex<double> exp_ikr = std::exp(-ikr);
  std::complex<double> factor = exp_ikr * Ids * mu0_ * (ikr + 1.0) / (4.0 * M_PI * r * r * r);
  factor = units.Nondimensionalize<Units::ValueType::FIELD_B>(factor);

  return {{0.0, factor * (-z), factor * y}};
}

// Compare the implementation in CurrentDipoleOperator with the analytic solution.
void runCurrentDipoleTest(double freq_Hz, std::unique_ptr<mfem::Mesh> serial_mesh,
                          const std::vector<int> &farfield_attributes,
                          const std::vector<int> &domain_attributes, double L0, double Lc,
                          double inner_frac = 0.2, double outer_frac = 0.8)
{
  constexpr double Ids = 1.;
  int Order = 2;

  Units units(L0, Lc);
  IoData iodata{units};
  iodata.domains.materials.emplace_back().attributes = domain_attributes;
  auto &dipole_config = iodata.domains.current_dipole[1];
  dipole_config.moment = Ids;
  dipole_config.direction = {1, 0, 0};
  dipole_config.center = {0.0, 0.0, 0.0};
  iodata.boundaries.farfield.attributes = farfield_attributes;
  iodata.boundaries.farfield.order = 2;  // TODO: Experiment with order 1
  iodata.problem.type = ProblemType::DRIVEN;
  iodata.solver.linear.type = LinearSolver::AMS;
  iodata.solver.order = Order;
  iodata.solver.linear.max_it = 1000;
  iodata.solver.linear.tol = 1e-10;
  iodata.solver.linear.mg_max_levels = 1;  // Disable multigrid for parallel compatibility
  // Enable ParaView output for visualization
  iodata.solver.driven.save_indices = {0};
  iodata.problem.output_formats.paraview = true;

  iodata.CheckConfiguration();

  auto comm = Mpi::World();
  const int num_ranks = Mpi::Size(comm);
  const int rank = Mpi::Rank(comm);

  if (rank == 0)
  {
    std::cout << "\nRunning current dipole test with " << num_ranks << " MPI ranks\n" << std::endl;

    serial_mesh->PrintInfo();
  }

  // Read parallel mesh.
  const int dim = serial_mesh->Dimension();
  REQUIRE(num_ranks <= serial_mesh->GetNE());
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  std::vector<std::unique_ptr<Mesh>> mesh_vec;
  mesh_vec.push_back(std::make_unique<Mesh>(palace_mesh));
  SpaceOperator space_op(iodata, mesh_vec);

  double omega = 2.0 * M_PI *
                 iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_Hz / 1e9);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  E = 0.0;
  B = 0.0;

  auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get(), A2.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ksp.SetOperators(*A, *P);
  space_op.GetExcitationVector(1, omega, RHS);

  // Solve for E
  ksp.Mult(RHS, E);

  // Compute B = -1/(iω) ∇ x E on the true dofs.
  Curl.Mult(E.Real(), B.Real());
  Curl.Mult(E.Imag(), B.Imag());
  B *= -1.0 / (1i * omega);

  // Recover complete GridFunction from constrained solver solution for error computation.
  GridFunction E_field(space_op.GetNDSpace(), true);
  E_field.Real().SetFromTrueDofs(E.Real());
  E_field.Imag().SetFromTrueDofs(E.Imag());

  GridFunction B_field(space_op.GetRTSpace(), true);
  B_field.Real().SetFromTrueDofs(B.Real());
  B_field.Imag().SetFromTrueDofs(B.Imag());

  // Create analytical solution GridFunction for visualization
  GridFunction E_analytical(space_op.GetNDSpace(), true);
  E_analytical = 0.0;

  GridFunction B_analytical(space_op.GetRTSpace(), true);
  B_analytical = 0.0;

  // Create analytical solution coefficients.
  auto E_exact_real = [=](const mfem::Vector &x, mfem::Vector &E_out)
  {
    auto E_complex = ComputeCurrentDipoleENonDim(x, iodata.units, Ids, freq_Hz);
    E_out(0) = std::real(E_complex[0]);
    E_out(1) = std::real(E_complex[1]);
    E_out(2) = std::real(E_complex[2]);
  };
  auto E_exact_imag = [=](const mfem::Vector &x, mfem::Vector &E_out)
  {
    auto E_complex = ComputeCurrentDipoleENonDim(x, iodata.units, Ids, freq_Hz);
    E_out(0) = std::imag(E_complex[0]);
    E_out(1) = std::imag(E_complex[1]);
    E_out(2) = std::imag(E_complex[2]);
  };
  mfem::VectorFunctionCoefficient E_exact_real_coef(dim, E_exact_real);
  mfem::VectorFunctionCoefficient E_exact_imag_coef(dim, E_exact_imag);

  auto B_exact_real = [=](const mfem::Vector &x, mfem::Vector &B_out)
  {
    auto B_complex = ComputeCurrentDipoleBNonDim(x, iodata.units, Ids, freq_Hz);
    B_out(0) = std::real(B_complex[0]);
    B_out(1) = std::real(B_complex[1]);
    B_out(2) = std::real(B_complex[2]);
  };
  auto B_exact_imag = [=](const mfem::Vector &x, mfem::Vector &B_out)
  {
    auto B_complex = ComputeCurrentDipoleBNonDim(x, iodata.units, Ids, freq_Hz);
    B_out(0) = std::imag(B_complex[0]);
    B_out(1) = std::imag(B_complex[1]);
    B_out(2) = std::imag(B_complex[2]);
  };
  mfem::VectorFunctionCoefficient B_exact_real_coef(dim, B_exact_real);
  mfem::VectorFunctionCoefficient B_exact_imag_coef(dim, B_exact_imag);

  // Project analytical solution onto GridFunction for visualization
  E_analytical.Real().ProjectCoefficient(E_exact_real_coef);
  E_analytical.Imag().ProjectCoefficient(E_exact_imag_coef);

  B_analytical.Real().ProjectCoefficient(B_exact_real_coef);
  B_analytical.Imag().ProjectCoefficient(B_exact_imag_coef);

  // Create element marker to exclude center and boundary regions from error computation
  auto &mesh_ref = palace_mesh.Get();
  ElementMarker element_marker(&mesh_ref, inner_frac, outer_frac);
  mfem::Array<int> *marked_elements = element_marker.GetMarkedElements();

  // Verify that marked_elements size matches local mesh element count
  int local_nelem = mesh_ref.GetNE();
  if (marked_elements->Size() != local_nelem)
  {
    std::cout << "\n*** ERROR: marked_elements size ("
              << marked_elements->Size() << ") != local mesh elements ("
              << local_nelem << ") ***\n";
  }

  // Higher-order integration for accuracy
  int order = 3;
  int order_quad = order + 2;
  const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
  for (int i = 0; i < mfem::Geometry::NumGeom; ++i)
  {
    irs[i] = &(mfem::IntRules.Get(i, order_quad));
  }

  if (rank == 0)
  {
    std::cout << "\n=== Integration Rule Info ===";
    std::cout << "\n  Polynomial order: " << order;
    std::cout << "\n  Quadrature order: " << order_quad;
    std::cout << "\n  Using same integration rules for error and norm computation\n";
  }

  // Compute L2 error only on marked elements (excluding center and boundary regions)
  double E_error_real = E_field.Real().ComputeL2Error(E_exact_real_coef, irs, marked_elements);
  double E_error_imag = E_field.Imag().ComputeL2Error(E_exact_imag_coef, irs, marked_elements);
  double E_total_error = std::sqrt(E_error_real * E_error_real + E_error_imag * E_error_imag);

  // Compute exact solution norm using zero field trick, also only on marked elements
  // When field=0, ComputeL2Error(exact_coef) returns ||0 - exact_coef||_L2 =
  // ||exact_coef||_L2
  GridFunction zero_field_for_norm_computation(space_op.GetNDSpace(), true);
  zero_field_for_norm_computation = 0.0;  // Initialize to zero
  double E_exact_norm_real =
      zero_field_for_norm_computation.Real().ComputeL2Error(E_exact_real_coef, irs, marked_elements);
  double E_exact_norm_imag =
      zero_field_for_norm_computation.Imag().ComputeL2Error(E_exact_imag_coef, irs, marked_elements);
  double E_total_exact_norm =
      std::sqrt(E_exact_norm_real * E_exact_norm_real + E_exact_norm_imag * E_exact_norm_imag);

  // Sanity check: exact norm should be non-zero
  if (rank == 0 && E_total_exact_norm < 1e-14)
  {
    std::cout << "\n*** WARNING: Exact solution norm is nearly zero! ***";
    std::cout << "\n*** This suggests the analytical solution may not be properly defined. ***\n";
  }

  // Also compute errors on the full domain for comparison
  double E_error_real_full = E_field.Real().ComputeL2Error(E_exact_real_coef, irs);
  double E_error_imag_full = E_field.Imag().ComputeL2Error(E_exact_imag_coef, irs);
  double E_total_error_full = std::sqrt(E_error_real_full * E_error_real_full + E_error_imag_full * E_error_imag_full);

  double E_exact_norm_real_full =
      zero_field_for_norm_computation.Real().ComputeL2Error(E_exact_real_coef, irs);
  double E_exact_norm_imag_full =
      zero_field_for_norm_computation.Imag().ComputeL2Error(E_exact_imag_coef, irs);
  double E_total_exact_norm_full =
      std::sqrt(E_exact_norm_real_full * E_exact_norm_real_full + E_exact_norm_imag_full * E_exact_norm_imag_full);

  // Compute relative errors
  double E_relative_error = E_total_error / E_total_exact_norm;
  double E_relative_error_full = E_total_error_full / E_total_exact_norm_full;

  // ========== B-FIELD ERROR COMPUTATION ==========
  // Compute L2 error for B-field only on marked elements
  double B_error_real = B_field.Real().ComputeL2Error(B_exact_real_coef, irs, marked_elements);
  double B_error_imag = B_field.Imag().ComputeL2Error(B_exact_imag_coef, irs, marked_elements);
  double B_total_error = std::sqrt(B_error_real * B_error_real + B_error_imag * B_error_imag);

  // Compute exact B-field norm using zero field trick
  GridFunction zero_B_field_for_norm_computation(space_op.GetRTSpace(), true);
  zero_B_field_for_norm_computation = 0.0;
  double B_exact_norm_real =
      zero_B_field_for_norm_computation.Real().ComputeL2Error(B_exact_real_coef, irs, marked_elements);
  double B_exact_norm_imag =
      zero_B_field_for_norm_computation.Imag().ComputeL2Error(B_exact_imag_coef, irs, marked_elements);
  double B_total_exact_norm =
      std::sqrt(B_exact_norm_real * B_exact_norm_real + B_exact_norm_imag * B_exact_norm_imag);

  // Sanity check: exact B norm should be non-zero
  if (rank == 0 && B_total_exact_norm < 1e-14)
  {
    std::cout << "\n*** WARNING: Exact B solution norm is nearly zero! ***";
    std::cout << "\n*** This suggests the analytical B solution may not be properly defined. ***\n";
  }

  // Also compute B errors on the full domain for comparison
  double B_error_real_full = B_field.Real().ComputeL2Error(B_exact_real_coef, irs);
  double B_error_imag_full = B_field.Imag().ComputeL2Error(B_exact_imag_coef, irs);
  double B_total_error_full = std::sqrt(B_error_real_full * B_error_real_full + B_error_imag_full * B_error_imag_full);

  double B_exact_norm_real_full =
      zero_B_field_for_norm_computation.Real().ComputeL2Error(B_exact_real_coef, irs);
  double B_exact_norm_imag_full =
      zero_B_field_for_norm_computation.Imag().ComputeL2Error(B_exact_imag_coef, irs);
  double B_total_exact_norm_full =
      std::sqrt(B_exact_norm_real_full * B_exact_norm_real_full + B_exact_norm_imag_full * B_exact_norm_imag_full);

  // Compute relative errors for B-field
  double B_relative_error = B_total_error / B_total_exact_norm;
  double B_relative_error_full = B_total_error_full / B_total_exact_norm_full;

  // ======================== VISUALIZATION ========================
  // Write both numerical and analytical solutions to ParaView
  {
    // Compute error field for E
    GridFunction E_error(space_op.GetNDSpace(), true);
    E_error.Real() = E_field.Real();
    E_error.Real() -= E_analytical.Real();
    E_error.Imag() = E_field.Imag();
    E_error.Imag() -= E_analytical.Imag();

    // Compute relative error field (pointwise) for E
    GridFunction E_rel_error(space_op.GetNDSpace(), true);
    E_rel_error = 0.0;

    // For each DOF, compute relative error separately for real and imaginary parts
    auto &E_err_real = E_error.Real();
    auto &E_err_imag = E_error.Imag();
    auto &E_ana_real = E_analytical.Real();
    auto &E_ana_imag = E_analytical.Imag();
    auto &E_rel_real = E_rel_error.Real();
    auto &E_rel_imag = E_rel_error.Imag();

    for (int i = 0; i < E_err_real.Size(); i++)
    {
      // Relative error for real part
      double ana_real_mag = std::abs(E_ana_real(i));
      if (ana_real_mag > 1e-14)
      {
        E_rel_real(i) = E_err_real(i) / ana_real_mag;
      }
      else
      {
        E_rel_real(i) = 0.0;  // Explicitly set to zero when analytical is too small
      }

      // Relative error for imaginary part
      double ana_imag_mag = std::abs(E_ana_imag(i));
      if (ana_imag_mag > 1e-14)
      {
        E_rel_imag(i) = E_err_imag(i) / ana_imag_mag;
      }
      else
      {
        E_rel_imag(i) = 0.0;  // Explicitly set to zero when analytical is too small
      }
    }

    // Compute error field for B
    GridFunction B_error(space_op.GetRTSpace(), true);
    B_error.Real() = B_field.Real();
    B_error.Real() -= B_analytical.Real();
    B_error.Imag() = B_field.Imag();
    B_error.Imag() -= B_analytical.Imag();

    // Compute relative error field (pointwise) for B
    GridFunction B_rel_error(space_op.GetRTSpace(), true);
    B_rel_error = 0.0;

    auto &B_err_real = B_error.Real();
    auto &B_err_imag = B_error.Imag();
    auto &B_ana_real = B_analytical.Real();
    auto &B_ana_imag = B_analytical.Imag();
    auto &B_rel_real = B_rel_error.Real();
    auto &B_rel_imag = B_rel_error.Imag();

    for (int i = 0; i < B_err_real.Size(); i++)
    {
      // Relative error for real part
      double ana_real_mag = std::abs(B_ana_real(i));
      if (ana_real_mag > 1e-14)
      {
        B_rel_real(i) = B_err_real(i) / ana_real_mag;
      }
      else
      {
        B_rel_real(i) = 0.0;  // Explicitly set to zero when analytical is too small
      }

      // Relative error for imaginary part
      double ana_imag_mag = std::abs(B_ana_imag(i));
      if (ana_imag_mag > 1e-14)
      {
        B_rel_imag(i) = B_err_imag(i) / ana_imag_mag;
      }
      else
      {
        B_rel_imag(i) = 0.0;  // Explicitly set to zero when analytical is too small
      }
    }

    auto &mesh = palace_mesh.Get();
    mfem::ParaViewDataCollection paraview_dc("CurrentDipoleComparison", &mesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetDataFormat(mfem::VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(freq_Hz);

    // Register E-field numerical solution
    paraview_dc.RegisterField("E_numerical_real", &E_field.Real());
    paraview_dc.RegisterField("E_numerical_imag", &E_field.Imag());

    // Register E-field analytical solution
    paraview_dc.RegisterField("E_analytical_real", &E_analytical.Real());
    paraview_dc.RegisterField("E_analytical_imag", &E_analytical.Imag());

    // Register E-field relative error (separate for real and imaginary parts)
    paraview_dc.RegisterField("E_relative_error_real", &E_rel_error.Real());
    paraview_dc.RegisterField("E_relative_error_imag", &E_rel_error.Imag());

    // Register B-field numerical solution
    paraview_dc.RegisterField("B_numerical_real", &B_field.Real());
    paraview_dc.RegisterField("B_numerical_imag", &B_field.Imag());

    // Register B-field analytical solution
    paraview_dc.RegisterField("B_analytical_real", &B_analytical.Real());
    paraview_dc.RegisterField("B_analytical_imag", &B_analytical.Imag());

    // Register B-field relative error (separate for real and imaginary parts)
    paraview_dc.RegisterField("B_relative_error_real", &B_rel_error.Real());
    paraview_dc.RegisterField("B_relative_error_imag", &B_rel_error.Imag());

    paraview_dc.Save();

    std::cout << "\n=== ParaView files written to: ParaView/CurrentDipoleComparison/ ===\n";
  }
  // ===============================================================

  // ------------------------Debugging--------------------------------
  std::cout << "\nTesting frequency: " << freq_Hz << " \n";

  // Get the scaling factor to automatically adjust test points for any units configuration
  double scaling_factor = iodata.units.GetMeshLengthRelativeScale();  // Lc_m / L0_m
  double point_scale = 1.0 / scaling_factor;
  std::cout << "\nMesh scaling factor (Lc/L0): " << scaling_factor;
  std::cout << "\nTest points scaled by: " << point_scale;

  // Debug analytical field values at key locations to understand the field pattern
  std::cout << "\n=== Analytical Field Debug at Key Points ===";

  // Define base test points, then scale them to fit within mesh bounds
  std::vector<std::pair<mfem::Vector, std::string>> base_test_points = {
      {mfem::Vector({0.0, 0.0, 0.0}), "Origin"},
      {mfem::Vector({0.02, 0.0, 0.0}), "+0.02 origin on x-axis"},
      {mfem::Vector({-0.1, 0.0, 0.0}), "-0.1 from origin on x-axis"},
      {mfem::Vector({-0.2, 0.3, 0.4}), "{-.2, .3, .4} from origin"},
      {mfem::Vector({0.0, 0.2, 0.0}), "0.2 field on y-axis"},
      {mfem::Vector({0.0, 0.0, -0.3}), "-0.3 field on z-axis"},
      {mfem::Vector({0.4, 0.0, 0.0}), "0.4 field on x-axis"},
      {mfem::Vector({0.0, 0.45, 0.0}), "+0.45 on y-axis"},
      {mfem::Vector({0.0, 0.0, 0.47}), "+0.47 on z-axis"},
      {mfem::Vector({0.495, 0.0, 0.0}), "+0.495 from origin on x-axis"},
      {mfem::Vector({0.0, 0.0, -0.5}), "-0.5 from origin on z-axis"},
  };

  // Scale the test points to fit within mesh bounds for the given units
  std::vector<std::pair<mfem::Vector, std::string>> test_points;
  for (const auto &[base_pt, description] : base_test_points)
  {
    mfem::Vector scaled_pt = base_pt;
    scaled_pt *= point_scale;
    test_points.emplace_back(scaled_pt, description + " (scaled)");
  }

  for (const auto &[pt, description] : test_points)
  {
    double r_nondim = pt.Norml2();
    double r_dim = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(r_nondim);

    // Evaluate nondimensional analytical solution
    auto E_analytical = ComputeCurrentDipoleENonDim(pt, iodata.units, Ids, freq_Hz);

    double analytical_magnitude =
        std::sqrt(std::real(E_analytical[0] * std::conj(E_analytical[0])) +
                  std::real(E_analytical[1] * std::conj(E_analytical[1])) +
                  std::real(E_analytical[2] * std::conj(E_analytical[2])));

    // Try to evaluate Palace FEM solution at this point
    mfem::Vector E_palace_real(dim), E_palace_imag(dim);
    E_palace_real = 0.0;
    E_palace_imag = 0.0;
    bool palace_eval_success = false;

    try
    {
      // Find which element contains this point
      auto &mesh = space_op.GetNDSpace().GetParMesh();
      mfem::Array<int> elem_ids;
      mfem::Array<mfem::IntegrationPoint> ips;

      // Convert Vector to DenseMatrix for FindPoints
      mfem::DenseMatrix point_mat(dim, 1);
      for (int i = 0; i < dim; i++)
      {
        point_mat(i, 0) = pt(i);
      }

      mesh.FindPoints(point_mat, elem_ids, ips);

      if (elem_ids.Size() > 0 && elem_ids[0] >= 0)
      {
        // Evaluate the gridfunction at this point
        E_field.Real().GetVectorValue(elem_ids[0], ips[0], E_palace_real);
        E_field.Imag().GetVectorValue(elem_ids[0], ips[0], E_palace_imag);
        palace_eval_success = true;
      }
    }
    catch (...)
    {
      // Evaluation failed - just continue with zero values
    }

    double palace_magnitude = std::sqrt(
        E_palace_real(0) * E_palace_real(0) + E_palace_real(1) * E_palace_real(1) +
        E_palace_real(2) * E_palace_real(2) + E_palace_imag(0) * E_palace_imag(0) +
        E_palace_imag(1) * E_palace_imag(1) + E_palace_imag(2) * E_palace_imag(2));

    std::cout << "\n"
              << description << " r=" << r_dim << "m:" << " r_nondim=" << r_nondim << "m:";
    std::cout << "\n  PALACE FEM SOLUTION:";
    if (palace_eval_success)
    {
      std::cout << "\n    Palace |E|: " << palace_magnitude;
      std::cout << "\n    Palace Ex: "
                << std::complex<double>(E_palace_real(0), E_palace_imag(0));
      std::cout << "\n    Palace Ey: "
                << std::complex<double>(E_palace_real(1), E_palace_imag(1));
      std::cout << "\n    Palace Ez: "
                << std::complex<double>(E_palace_real(2), E_palace_imag(2));
    }
    else
    {
      std::cout << "\n    Palace evaluation failed at this point";
    }
    std::cout << "\n  ANALYTICAL SOLUTION:";
    std::cout << "\n    Analytical |E|: " << analytical_magnitude;
    std::cout << "\n    Analytical Ex: " << E_analytical[0];
    std::cout << "\n    Analytical Ey: " << E_analytical[1];
    std::cout << "\n    Analytical Ez: " << E_analytical[2];
    if (palace_eval_success && analytical_magnitude > 1e-12)
    {
      std::cout << "\n  COMPARISON:";
      std::cout << "\n    |E| ratio (Palace/Analytical): "
                << palace_magnitude / analytical_magnitude;
    }
  }
  std::cout << "\n============================================\n";

  std::cout << "\n\n=== L2 ERROR ANALYSIS ===";
  std::cout << "Polynomial order " << Order << std::endl;

  std::cout << "\n\n--- E-FIELD: Selective Elements (0.2-0.8 domain) ---";
  std::cout << "\nAbsolute L2 Errors:";
  std::cout << "\n  Real part error: || E_h_Re - E_Re ||_L2 = " << E_error_real;
  std::cout << "\n  Imag part error: || E_h_Im - E_Im ||_L2 = " << E_error_imag;
  std::cout << "\n  Total error:     || E_h - E ||_L2 = " << E_total_error;

  std::cout << "\nExact Solution Norms:";
  std::cout << "\n  Real part norm: || E_Re ||_L2 = " << E_exact_norm_real;
  std::cout << "\n  Imag part norm: || E_Im ||_L2 = " << E_exact_norm_imag;
  std::cout << "\n  Total norm:     || E ||_L2 = " << E_total_exact_norm;

  std::cout << "\nRelative L2 Errors:";
  std::cout << "\n  Real part: || E_h_Re - E_Re || / ||E_Re|| = "
            << E_error_real / E_exact_norm_real;
  std::cout << "\n  Imag part: || E_h_Im - E_Im || / ||E_Im|| = "
            << E_error_imag / E_exact_norm_imag;
  std::cout << "\n  Total:     || E_h - E || / ||E|| = " << E_relative_error;

  std::cout << "\n\n--- E-FIELD: Full Domain (for comparison) ---";
  std::cout << "\nAbsolute L2 Errors:";
  std::cout << "\n  Real part error: || E_h_Re - E_Re ||_L2 = " << E_error_real_full;
  std::cout << "\n  Imag part error: || E_h_Im - E_Im ||_L2 = " << E_error_imag_full;
  std::cout << "\n  Total error:     || E_h - E ||_L2 = " << E_total_error_full;

  std::cout << "\nExact Solution Norms:";
  std::cout << "\n  Real part norm: || E_Re ||_L2 = " << E_exact_norm_real_full;
  std::cout << "\n  Imag part norm: || E_Im ||_L2 = " << E_exact_norm_imag_full;
  std::cout << "\n  Total norm:     || E ||_L2 = " << E_total_exact_norm_full;

  std::cout << "\nRelative L2 Errors:";
  std::cout << "\n  Real part: || E_h_Re - E_Re || / ||E_Re|| = "
            << E_error_real_full / E_exact_norm_real_full;
  std::cout << "\n  Imag part: || E_h_Im - E_Im || / ||E_Im|| = "
            << E_error_imag_full / E_exact_norm_imag_full;
  std::cout << "\n  Total:     || E_h - E || / ||E|| = " << E_relative_error_full;

  std::cout << "\n\n--- E-FIELD: Error Comparison ---";
  std::cout << "\nSelective vs Full Domain Error Ratios:";
  std::cout << "\n  Selective/Full Error Ratio: " << E_total_error / E_total_error_full;
  std::cout << "\n  Selective/Full Relative Error Ratio: " << E_relative_error / E_relative_error_full;

  std::cout << "\n\n--- B-FIELD: Selective Elements (0.2-0.8 domain) ---";
  std::cout << "\nAbsolute L2 Errors:";
  std::cout << "\n  Real part error: || B_h_Re - B_Re ||_L2 = " << B_error_real;
  std::cout << "\n  Imag part error: || B_h_Im - B_Im ||_L2 = " << B_error_imag;
  std::cout << "\n  Total error:     || B_h - B ||_L2 = " << B_total_error;

  std::cout << "\nExact Solution Norms:";
  std::cout << "\n  Real part norm: || B_Re ||_L2 = " << B_exact_norm_real;
  std::cout << "\n  Imag part norm: || B_Im ||_L2 = " << B_exact_norm_imag;
  std::cout << "\n  Total norm:     || B ||_L2 = " << B_total_exact_norm;

  std::cout << "\nRelative L2 Errors:";
  std::cout << "\n  Real part: || B_h_Re - B_Re || / ||B_Re|| = "
            << B_error_real / B_exact_norm_real;
  std::cout << "\n  Imag part: || B_h_Im - B_Im || / ||B_Im|| = "
            << B_error_imag / B_exact_norm_imag;
  std::cout << "\n  Total:     || B_h - B || / ||B|| = " << B_relative_error;

  std::cout << "\n\n--- B-FIELD: Full Domain (for comparison) ---";
  std::cout << "\nAbsolute L2 Errors:";
  std::cout << "\n  Real part error: || B_h_Re - B_Re ||_L2 = " << B_error_real_full;
  std::cout << "\n  Imag part error: || B_h_Im - B_Im ||_L2 = " << B_error_imag_full;
  std::cout << "\n  Total error:     || B_h - B ||_L2 = " << B_total_error_full;

  std::cout << "\nExact Solution Norms:";
  std::cout << "\n  Real part norm: || B_Re ||_L2 = " << B_exact_norm_real_full;
  std::cout << "\n  Imag part norm: || B_Im ||_L2 = " << B_exact_norm_imag_full;
  std::cout << "\n  Total norm:     || B ||_L2 = " << B_total_exact_norm_full;

  std::cout << "\nRelative L2 Errors:";
  std::cout << "\n  Real part: || B_h_Re - B_Re || / ||B_Re|| = "
            << B_error_real_full / B_exact_norm_real_full;
  std::cout << "\n  Imag part: || B_h_Im - B_Im || / ||B_Im|| = "
            << B_error_imag_full / B_exact_norm_imag_full;
  std::cout << "\n  Total:     || B_h - B || / ||B|| = " << B_relative_error_full;

  std::cout << "\n\n--- B-FIELD: Error Comparison ---";
  std::cout << "\nSelective vs Full Domain Error Ratios:";
  std::cout << "\n  Selective/Full Error Ratio: " << B_total_error / B_total_error_full;
  std::cout << "\n  Selective/Full Relative Error Ratio: " << B_relative_error / B_relative_error_full;
  std::cout << "\n\n";
  // -----------------------------------------------------------------------

  CHECK_THAT(E_relative_error, WithinAbs(0.0, 0.1));
  CHECK_THAT(B_relative_error, WithinAbs(0.0, 0.12));
}

TEST_CASE("Electrical Current Dipole in a Cube", "[currentdipole][cube][Serial]")
{
  double freq_Hz = 1e8;
  std::vector<int> attributes = {1, 2, 3, 4, 5, 6};

  // Make mesh for a cube [0, 1] x [0, 1] x [0, 1]
  int resolution = 21;  // TODO: Switch back to 20 after fixing the vertex issue
  // TODO: Try TETRAHEDRON after fixing the vertex issue
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::HEXAHEDRON));
  // Transform to center the origin
  serial_mesh->Transform(
      [](const mfem::Vector &x, mfem::Vector &p)
      {
        p = x;
        p(0) -= 0.5;  // Transform [0,1] -> [-0.5,0.5]
        p(1) -= 0.5;
        p(2) -= 0.5;
      });

  // To make the mesh larger, we can increase the first input of Units, L0, instead of
  // creating a new mesh.
  runCurrentDipoleTest(freq_Hz, std::move(serial_mesh), attributes, {1}, 5.0, 1.0, 0.1, 0.9);
}

TEST_CASE("Electrical Current Dipole in a Sphere", "[currentdipole][sphere][Serial]")
{
  double freq_Hz = 2e8;

  // Load the antenna sphere mesh (n = n_farfield: 5, 7, 10, 13)
  // Corresponding n_circle values: 6, 8, 12, 15
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>("../examples/antenna/mesh/short_n10.msh");

  std::vector<int> domain_attributes = {2};    // Domain volume
  std::vector<int> farfield_attributes = {1};  // Outer boundary (absorbing boundary)
  runCurrentDipoleTest(freq_Hz, std::move(serial_mesh), farfield_attributes,
                       domain_attributes, 1.0, 1.0, 0.1, 0.95);
}

}  // namespace
}  // namespace palace
