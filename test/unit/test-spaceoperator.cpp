// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>
#include <Eigen/Eigenvalues>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "fem/singularfeatures.hpp"
#include "linalg/divfree.hpp"
#include "linalg/vector.hpp"
#include "models/hierarchicalmaxwellestimator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

mfem::Mesh InternalSheetMesh()
{
  mfem::Mesh mesh(3, 6, 4, 10, 3);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(1.0, 1.0, 0.0);
  mesh.AddVertex(0.5, 0.5, 1.0);
  mesh.AddVertex(0.5, 0.5, -1.0);
  mesh.AddTet(0, 1, 2, 4, 1);
  mesh.AddTet(0, 2, 1, 5, 1);
  mesh.AddTet(1, 3, 2, 4, 1);
  mesh.AddTet(1, 2, 3, 5, 1);
  mesh.AddBdrTriangle(0, 1, 2, 7);
  mesh.AddBdrTriangle(1, 3, 2, 7);
  mesh.AddBdrTriangle(0, 1, 4, 1);
  mesh.AddBdrTriangle(2, 0, 4, 1);
  mesh.AddBdrTriangle(1, 3, 4, 1);
  mesh.AddBdrTriangle(3, 2, 4, 1);
  mesh.AddBdrTriangle(1, 0, 5, 1);
  mesh.AddBdrTriangle(0, 2, 5, 1);
  mesh.AddBdrTriangle(3, 1, 5, 1);
  mesh.AddBdrTriangle(2, 3, 5, 1);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

mfem::Mesh InternalLineTipMesh()
{
  mfem::Mesh mesh(2, 9, 8, 9, 2);
  mesh.AddVertex(-1.0, -1.0);
  mesh.AddVertex(0.0, -1.0);
  mesh.AddVertex(1.0, -1.0);
  mesh.AddVertex(-1.0, 0.0);
  mesh.AddVertex(0.0, 0.0);
  mesh.AddVertex(1.0, 0.0);
  mesh.AddVertex(-1.0, 1.0);
  mesh.AddVertex(0.0, 1.0);
  mesh.AddVertex(1.0, 1.0);
  mesh.AddTriangle(0, 1, 4, 1);
  mesh.AddTriangle(0, 4, 3, 1);
  mesh.AddTriangle(3, 4, 7, 1);
  mesh.AddTriangle(3, 7, 6, 1);
  mesh.AddTriangle(1, 2, 5, 1);
  mesh.AddTriangle(1, 5, 4, 1);
  mesh.AddTriangle(4, 5, 8, 1);
  mesh.AddTriangle(4, 8, 7, 1);
  mesh.AddBdrSegment(0, 1, 1);
  mesh.AddBdrSegment(1, 2, 1);
  mesh.AddBdrSegment(2, 5, 1);
  mesh.AddBdrSegment(5, 8, 1);
  mesh.AddBdrSegment(8, 7, 1);
  mesh.AddBdrSegment(7, 6, 1);
  mesh.AddBdrSegment(6, 3, 1);
  mesh.AddBdrSegment(3, 0, 1);
  mesh.AddBdrSegment(3, 4, 7);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

void SetQuadraticGeometry(mfem::Mesh &mesh, bool curved)
{
  const int dimension = mesh.SpaceDimension();
  mesh.SetCurvature(2, false, dimension, mfem::Ordering::byVDIM);
  if (!curved)
  {
    return;
  }
  mfem::VectorFunctionCoefficient geometry(
      dimension,
      [dimension](const mfem::Vector &x, mfem::Vector &value)
      {
        value.SetSize(dimension);
        value = x;
        if (dimension == 2)
        {
          // The selected segment has y = 0 and remains exactly straight.
          value[1] += 0.1 * x[0] * x[1];
        }
        else
        {
          // The selected sheet has z = 0 and its perimeter remains straight.
          value[2] += 0.1 * x[0] * x[2];
        }
      });
  mesh.GetNodes()->ProjectCoefficient(geometry);
}

IoData SingularSpaceData(int dimension, int order = 1, int mg_max_levels = 1,
                         bool heterogeneous = false, bool lumped_port = false,
                         bool open_exterior = false)
{
  config::json materials = {
      {{"Attributes", {1}}, {"Permittivity", 2.3}, {"Permeability", 1.7}}};
  if (heterogeneous)
  {
    materials.push_back(
        {{"Attributes", {2}}, {"Permittivity", 5.1}, {"Permeability", 1.2}});
  }
  config::json boundaries = {
      {"PEC", {{"Attributes", open_exterior ? config::json{7} : config::json{1, 7}}}}};
  if (lumped_port)
  {
    REQUIRE(dimension == 3);
    boundaries = {{"PEC", {{"Attributes", {1}}}},
                  {"LumpedPort",
                   {{{"Index", 1},
                     {"Attributes", {7}},
                     {"R", 50.0},
                     {"L", 0.7},
                     {"C", 0.2},
                     {"Direction", "+X"}}}}};
  }
  config::json config = {
      {"Problem",
       {{"Type",
         dimension == 2 ? "BoundaryMode" : (lumped_port ? "Eigenmode" : "Electrostatic")},
        {"Output", "test_output"}}},
      {"Model", {{"Mesh", "unused.mesh"}}},
      {"Domains", {{"Materials", std::move(materials)}}},
      {"Boundaries", std::move(boundaries)},
      {"Solver",
       {{"Order", 1},
        {"SingularElements",
         {{"Attributes", {7}},
          {"Order", 1},
          {"QuadratureOrder", 8},
          {"AbsTol", 1.0e-3},
          {"RelTol", 1.0e-3},
          {"MaxSubdivisions", 6}}},
        {"Linear",
         {{"Type", "STRUMPACK"}, {"MGMaxLevels", 1}, {"Tol", 1.0e-12}, {"MaxIts", 20}}}}}};
  IoData iodata(config, false);
  iodata.problem.type = ProblemType::EIGENMODE;
  iodata.solver.order = order;
  iodata.solver.linear.mg_max_levels = mg_max_levels;
  return iodata;
}

void FillVector(Vector &vector, double phase)
{
  for (int i = 0; i < vector.Size(); i++)
  {
    vector[i] = std::sin(phase + 0.37 * (i + 1));
  }
}

double RelativeNorm(const Vector &vector, const Vector &reference)
{
  return linalg::Norml2(Mpi::World(), vector) /
         std::max(1.0, linalg::Norml2(Mpi::World(), reference));
}

Eigen::MatrixXd AssembleReducedDenseMatrix(const Operator &matrix,
                                           const mfem::Array<int> &essential)
{
  REQUIRE(matrix.Height() == matrix.Width());
  std::vector<bool> constrained(matrix.Width(), false);
  for (int dof : essential)
  {
    REQUIRE(dof >= 0);
    REQUIRE(dof < matrix.Width());
    constrained[dof] = true;
  }
  std::vector<int> free_dofs;
  for (int dof = 0; dof < matrix.Width(); dof++)
  {
    if (!constrained[dof])
    {
      free_dofs.push_back(dof);
    }
  }
  REQUIRE_FALSE(free_dofs.empty());

  Eigen::MatrixXd dense(free_dofs.size(), free_dofs.size());
  Vector input(matrix.Width()), output(matrix.Height());
  for (std::size_t column = 0; column < free_dofs.size(); column++)
  {
    input = 0.0;
    input[free_dofs[column]] = 1.0;
    matrix.Mult(input, output);
    for (std::size_t row = 0; row < free_dofs.size(); row++)
    {
      dense(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
          output[free_dofs[row]];
    }
  }
  return dense;
}

void CheckLosslessSpectralEmbedding(const Operator &standard_stiffness,
                                    const Operator &standard_mass,
                                    const mfem::Array<int> &standard_essential,
                                    const Operator &combined_stiffness,
                                    const Operator &combined_mass,
                                    const mfem::Array<int> &combined_essential)
{
  const Eigen::MatrixXd standard_K =
      AssembleReducedDenseMatrix(standard_stiffness, standard_essential);
  const Eigen::MatrixXd standard_M =
      AssembleReducedDenseMatrix(standard_mass, standard_essential);
  const Eigen::MatrixXd combined_K =
      AssembleReducedDenseMatrix(combined_stiffness, combined_essential);
  const Eigen::MatrixXd combined_M =
      AssembleReducedDenseMatrix(combined_mass, combined_essential);
  REQUIRE(combined_K.rows() >= standard_K.rows());

  const double standard_scale =
      std::max({1.0, standard_K.cwiseAbs().maxCoeff(), standard_M.cwiseAbs().maxCoeff()});
  const double combined_scale =
      std::max({1.0, combined_K.cwiseAbs().maxCoeff(), combined_M.cwiseAbs().maxCoeff()});
  CHECK((standard_K - standard_K.transpose()).cwiseAbs().maxCoeff() <=
        2.0e-12 * standard_scale);
  CHECK((standard_M - standard_M.transpose()).cwiseAbs().maxCoeff() <=
        2.0e-12 * standard_scale);
  CHECK((combined_K - combined_K.transpose()).cwiseAbs().maxCoeff() <=
        2.0e-12 * combined_scale);
  CHECK((combined_M - combined_M.transpose()).cwiseAbs().maxCoeff() <=
        2.0e-12 * combined_scale);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> standard_mass_solver(standard_M);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> combined_mass_solver(combined_M);
  REQUIRE(standard_mass_solver.info() == Eigen::Success);
  REQUIRE(combined_mass_solver.info() == Eigen::Success);
  CHECK(standard_mass_solver.eigenvalues().minCoeff() >
        256.0 * std::numeric_limits<double>::epsilon() *
            standard_mass_solver.eigenvalues().maxCoeff());
  CHECK(combined_mass_solver.eigenvalues().minCoeff() >
        256.0 * std::numeric_limits<double>::epsilon() *
            combined_mass_solver.eigenvalues().maxCoeff());

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> standard_solver(standard_K,
                                                                            standard_M);
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> combined_solver(combined_K,
                                                                            combined_M);
  REQUIRE(standard_solver.info() == Eigen::Success);
  REQUIRE(combined_solver.info() == Eigen::Success);
  const auto &standard_eigenvalues = standard_solver.eigenvalues();
  const auto &combined_eigenvalues = combined_solver.eigenvalues();
  const Eigen::Index dimension_difference =
      combined_eigenvalues.size() - standard_eigenvalues.size();
  const double spectral_scale = std::max({1.0, standard_eigenvalues.cwiseAbs().maxCoeff(),
                                          combined_eigenvalues.cwiseAbs().maxCoeff()});
  const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() * spectral_scale;
  for (Eigen::Index i = 0; i < standard_eigenvalues.size(); i++)
  {
    CAPTURE(i, dimension_difference, standard_eigenvalues[i], combined_eigenvalues[i],
            combined_eigenvalues[i + dimension_difference], tolerance);
    CHECK(combined_eigenvalues[i] <= standard_eigenvalues[i] + tolerance);
    CHECK(standard_eigenvalues[i] <=
          combined_eigenvalues[i + dimension_difference] + tolerance);
  }
}

void CheckStandardBlockAction(const Operator &combined, const Operator &standard)
{
  REQUIRE(combined.Height() >= standard.Height());
  REQUIRE(combined.Width() >= standard.Width());
  REQUIRE(standard.Height() == standard.Width());

  Vector standard_input(standard.Width()), standard_action(standard.Height());
  Vector combined_input(combined.Width()), combined_action(combined.Height());
  FillVector(standard_input, 0.53);
  combined_input = 0.0;
  linalg::SetSubVector(combined_input, 0, standard_input);

  standard.Mult(standard_input, standard_action);
  combined.Mult(combined_input, combined_action);
  Vector combined_standard_action;
  combined_standard_action.MakeRef(combined_action, 0, standard.Height());
  combined_standard_action -= standard_action;
  CHECK(RelativeNorm(combined_standard_action, standard_action) < 2.0e-12);
}

void CheckStandardBlockAction(const ComplexOperator &combined,
                              const ComplexOperator &standard)
{
  REQUIRE(combined.Height() >= standard.Height());
  REQUIRE(combined.Width() >= standard.Width());
  REQUIRE(standard.Height() == standard.Width());

  ComplexVector standard_input(standard.Width()), standard_action(standard.Height());
  ComplexVector combined_input(combined.Width()), combined_action(combined.Height());
  FillVector(standard_input.Real(), 0.53);
  FillVector(standard_input.Imag(), 0.79);
  combined_input = 0.0;
  linalg::SetSubVector(combined_input.Real(), 0, standard_input.Real());
  linalg::SetSubVector(combined_input.Imag(), 0, standard_input.Imag());

  standard.Mult(standard_input, standard_action);
  combined.Mult(combined_input, combined_action);
  Vector combined_standard_real, combined_standard_imaginary;
  combined_standard_real.MakeRef(combined_action.Real(), 0, standard.Height());
  combined_standard_imaginary.MakeRef(combined_action.Imag(), 0, standard.Height());
  combined_standard_real -= standard_action.Real();
  combined_standard_imaginary -= standard_action.Imag();
  CHECK(RelativeNorm(combined_standard_real, standard_action.Real()) < 2.0e-12);
  CHECK(RelativeNorm(combined_standard_imaginary, standard_action.Imag()) < 2.0e-12);
}

void CheckSymmetricPositive(const Operator &matrix)
{
  Vector x(matrix.Width()), y(matrix.Width()), Ax(matrix.Height()), Ay(matrix.Height());
  FillVector(x, 0.2);
  FillVector(y, 0.9);
  matrix.Mult(x, Ax);
  matrix.Mult(y, Ay);
  const double x_A_y = linalg::Dot(Mpi::World(), x, Ay);
  const double y_A_x = linalg::Dot(Mpi::World(), y, Ax);
  CHECK(std::abs(x_A_y - y_A_x) <=
        2.0e-10 * std::max({1.0, std::abs(x_A_y), std::abs(y_A_x)}));
  CHECK(linalg::Dot(Mpi::World(), x, Ax) > 0.0);
}

void CheckSpaceOperator(mfem::Mesh serial_mesh, bool curved, double loss_tangent = 0.0,
                        int order = 1, int mg_max_levels = 1, bool heterogeneous = false,
                        bool lumped_port = false, bool open_exterior = false,
                        bool hierarchical = false)
{
  REQUIRE(Mpi::Size(Mpi::World()) == 1);
  const int dimension = serial_mesh.Dimension();
  if (heterogeneous)
  {
    for (int element = 0; element < serial_mesh.GetNE(); element++)
    {
      serial_mesh.GetElement(element)->SetAttribute(1 + element % 2);
    }
  }
  SetQuadraticGeometry(serial_mesh, curved);
  const bool tetrahedral = dimension == 3;
  const auto sheet_features =
      tetrahedral ? fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7})
                  : fem::singular::FeatureTopology{};
  const auto line_features =
      tetrahedral ? fem::singular::TriangleFeatureTopology{}
                  : fem::singular::ExtractSerialLineTipFeatures(serial_mesh, {7});

  std::vector<fem::singular::GlobalVertexId> source_vertex_ids(serial_mesh.GetNV());
  std::vector<fem::singular::GlobalVertexId> source_element_ids(serial_mesh.GetNE());
  std::iota(source_vertex_ids.begin(), source_vertex_ids.end(), 0);
  std::iota(source_element_ids.begin(), source_element_ids.end(), 0);
  auto parallel_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  const auto local_sheet_features =
      tetrahedral
          ? fem::singular::DistributeSerialSheetFeatures(
                sheet_features, *parallel_mesh, source_vertex_ids, source_element_ids)
          : fem::singular::FeatureTopology{};
  const auto local_line_features =
      tetrahedral
          ? fem::singular::TriangleFeatureTopology{}
          : fem::singular::DistributeSerialLineTipFeatures(
                line_features, *parallel_mesh, source_vertex_ids, source_element_ids);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));
  auto iodata = SingularSpaceData(dimension, order, mg_max_levels, heterogeneous,
                                  lumped_port, open_exterior);
  for (auto &material : iodata.domains.materials)
  {
    material.tandelta.s = {loss_tangent, loss_tangent, loss_tangent};
  }
  SpaceOperator space_op(iodata, meshes, tetrahedral ? &local_sheet_features : nullptr,
                         tetrahedral ? nullptr : &local_line_features, &source_vertex_ids);
  REQUIRE(space_op.HasSingularEnrichment());
  REQUIRE(space_op.GetSingularMaterialCoefficients().size() ==
          static_cast<std::size_t>(meshes.back()->GetNE()));
  REQUIRE(space_op.GetSingularImagMaterialCoefficients().size() ==
          static_cast<std::size_t>(meshes.back()->GetNE()));
  REQUIRE(space_op.GetSingularAbsMaterialCoefficients().size() ==
          static_cast<std::size_t>(meshes.back()->GetNE()));
  REQUIRE((space_op.GetSingularFeatures() != nullptr) == tetrahedral);
  REQUIRE((space_op.GetTriangleSingularFeatures() != nullptr) == !tetrahedral);
  REQUIRE(space_op.GetNDDbcAttributes().Size() > 0);
  CHECK(space_op.GetSingularAssemblyOptions().quadrature_order ==
        iodata.solver.singular_elements.quadrature_order);

  if (hierarchical)
  {
    HierarchicalMaxwellDomainData hierarchy(space_op);
    CHECK(hierarchy.GetInjection().consistency_error < 1.0e-12);
    CHECK(hierarchy.GetInjection().nonidentity_transformations > 0);
    const std::complex<double> omega(0.7, 0.08);
    const auto complex_operator = hierarchy.BuildComplexPolynomialContributions(omega);
    const auto at_zero = hierarchy.BuildComplexPolynomialContributions(0.0);
    const auto at_one = hierarchy.BuildComplexPolynomialContributions(1.0);
    const auto at_minus_one = hierarchy.BuildComplexPolynomialContributions(-1.0);
    const auto metric = hierarchy.BuildPolynomialMetricContributions(omega);
    const auto metric_one = hierarchy.BuildPolynomialMetricContributions(1.0);
    const auto metric_two = hierarchy.BuildPolynomialMetricContributions(2.0);
    REQUIRE(complex_operator.size() > static_cast<std::size_t>(meshes.back()->GetNE()));
    REQUIRE(metric.size() == complex_operator.size());
    REQUIRE(at_zero.size() == complex_operator.size());
    REQUIRE(at_one.size() == complex_operator.size());
    REQUIRE(at_minus_one.size() == complex_operator.size());
    REQUIRE(metric_one.size() == complex_operator.size());
    REQUIRE(metric_two.size() == complex_operator.size());
    double formula_error = 0.0, metric_error = 0.0, loss_norm = 0.0, damping_norm = 0.0,
           symmetry_error = 0.0;
    const std::complex<double> mass_scale = -omega * omega;
    for (std::size_t element = 0; element < complex_operator.size(); element++)
    {
      const auto &K = at_zero[element].matrix_real;
      for (int row = 0; row < K.Height(); row++)
      {
        for (int column = 0; column < K.Width(); column++)
        {
          // omega=+/-1 separates damping from dielectric loss:
          // A(+1)=K-Mr+i(C-Mi), A(-1)=K-Mr-i(C+Mi).
          const double Mr =
              K(row, column) - 0.5 * (at_one[element].matrix_real(row, column) +
                                      at_minus_one[element].matrix_real(row, column));
          const double C = 0.5 * (at_one[element].matrix_imag(row, column) -
                                  at_minus_one[element].matrix_imag(row, column));
          const double Mi = -0.5 * (at_one[element].matrix_imag(row, column) +
                                    at_minus_one[element].matrix_imag(row, column));
          const double expected_real = K(row, column) - omega.imag() * C +
                                       mass_scale.real() * Mr - mass_scale.imag() * Mi;
          const double expected_imag =
              omega.real() * C + mass_scale.imag() * Mr + mass_scale.real() * Mi;
          formula_error = std::max(
              formula_error,
              std::abs(complex_operator[element].matrix_real(row, column) - expected_real));
          formula_error = std::max(
              formula_error,
              std::abs(complex_operator[element].matrix_imag(row, column) - expected_imag));
          loss_norm = std::max(loss_norm, std::abs(Mi));
          damping_norm = std::max(damping_norm, std::abs(C));
          const double B1 = metric_one[element].matrix(row, column);
          const double B2 = metric_two[element].matrix(row, column);
          const double Mabs = 0.5 * (B2 - 2.0 * B1 + K(row, column));
          const double Cabs = B1 - K(row, column) - Mabs;
          metric_error = std::max(
              metric_error, std::abs(metric[element].matrix(row, column) - K(row, column) -
                                     std::abs(omega) * Cabs - std::norm(omega) * Mabs));
          symmetry_error = std::max(
              symmetry_error, std::abs(complex_operator[element].matrix_real(row, column) -
                                       complex_operator[element].matrix_real(column, row)));
          symmetry_error = std::max(
              symmetry_error, std::abs(complex_operator[element].matrix_imag(row, column) -
                                       complex_operator[element].matrix_imag(column, row)));
        }
      }
    }
    CAPTURE(formula_error, metric_error, loss_norm, damping_norm, symmetry_error);
    CHECK(formula_error < 1.0e-12);
    CHECK(metric_error < 1.0e-12);
    CHECK(loss_norm > 0.0);
    CHECK(damping_norm > 0.0);
    CHECK(symmetry_error < 1.0e-12);
    CHECK_THROWS(hierarchy.BuildPolynomialMetricContributions(0.0));
    const int combined_size =
        hierarchy.GetFineStandardSize() + hierarchy.GetEnrichmentSize();
    mfem::Vector xr(combined_size), xi(combined_size);
    for (int dof = 0; dof < combined_size; dof++)
    {
      xr(dof) = std::sin(0.31 * dof + 0.2);
      xi(dof) = std::cos(0.23 * dof - 0.1);
      if (hierarchy.GetFineEssentialMask()[dof])
      {
        xr(dof) = 0.0;
        xi(dof) = 0.0;
      }
    }
    const auto residual = fem::hierarchical::AssembleComplexResidual(
        combined_size, complex_operator, xr, xi, hierarchy.GetFineEssentialMask());
    const auto lifting = fem::hierarchical::LiftComplexResidualByPatches(
        meshes.back()->Get(), space_op.GetNDSpace().Get(), hierarchy.GetFineNDSpace().Get(),
        hierarchy.GetInjection(), metric, hierarchy.GetFineEssentialMask(),
        hierarchy.GetCoarseEssentialMask(), residual,
        hierarchy.GetElementEnrichmentGuests());
    CHECK(lifting.energy > 0.0);
    CHECK_THAT(std::accumulate(lifting.indicator.begin(), lifting.indicator.end(), 0.0),
               Catch::Matchers::WithinRel(lifting.energy, 1.0e-12));
    CHECK(lifting.real.face_patches > 0);
    CHECK(lifting.imag.face_patches > 0);

    ComplexVector coarse_field(space_op.GetNDTrueVSize());
    for (int dof = 0; dof < coarse_field.Size(); dof++)
    {
      coarse_field.Real()(dof) = std::sin(0.19 * dof + 0.3);
      coarse_field.Imag()(dof) = std::cos(0.29 * dof - 0.2);
    }
    for (int dof : space_op.GetCombinedNDDbcTDofList())
    {
      coarse_field.Real()(dof) = 0.0;
      coarse_field.Imag()(dof) = 0.0;
    }
    const auto parallel_estimate =
        hierarchy.EstimatePolynomialEigenResidual(omega, coarse_field);
    CHECK(parallel_estimate.total_energy > 0.0);
    CHECK_THAT(parallel_estimate.indicator_energy.Sum(),
               Catch::Matchers::WithinRel(parallel_estimate.total_energy, 1.0e-12));
    // The certified entity-patch shape is available on one rank and reproduces the
    // engine's summed real/imaginary lifting energies through the adapter.
    REQUIRE(hierarchy.EntityPatchesAvailable());
    const auto entity_estimate = hierarchy.EstimatePolynomialEigenResidual(
        omega, coarse_field, HierarchicalMaxwellDomainData::PatchShape::ENTITY);
    CHECK(entity_estimate.total_energy > 0.0);
    CHECK_THAT(entity_estimate.indicator_energy.Sum(),
               Catch::Matchers::WithinRel(entity_estimate.total_energy, 1.0e-12));
  }

  auto K = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ONE);
  auto K_zero = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);
  auto M_zero = space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  auto C = space_op.GetDampingMatrix<Operator>(Operator::DIAG_ZERO);
  auto K_bulk = space_op.GetBulkStiffnessMatrix(Operator::DIAG_ZERO);
  auto M_bulk = space_op.GetBulkMassMatrix(Operator::DIAG_ZERO);
  SpaceOperator standard_space_op(iodata.solver, iodata.domains, iodata.boundaries,
                                  iodata.problem.type, iodata.units, meshes);
  auto standard_K = standard_space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ONE);
  auto standard_K_zero =
      standard_space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto standard_M = standard_space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);
  auto standard_M_zero = standard_space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  auto standard_C = standard_space_op.GetDampingMatrix<Operator>(Operator::DIAG_ZERO);
  const auto *hypre_G =
      dynamic_cast<const mfem::HypreParMatrix *>(&space_op.GetGradMatrix());
  REQUIRE(hypre_G);

  CHECK(K->Height() == space_op.GetNDTrueVSize());
  CHECK(M->Height() == space_op.GetNDTrueVSize());
  CHECK(hypre_G->Height() == space_op.GetNDTrueVSize());
  CHECK(hypre_G->Width() == space_op.GetH1TrueVSize());
  CHECK(K->Height() > space_op.GetNDSpace().GetTrueVSize());
  CHECK(hypre_G->Width() > space_op.GetH1Space().GetTrueVSize());
  CHECK(space_op.GlobalTrueVSize() == K->Height());
  CheckSymmetricPositive(*K);
  CheckSymmetricPositive(*M);
  CheckStandardBlockAction(*K, *standard_K);
  CheckStandardBlockAction(*K_zero, *standard_K_zero);
  CheckStandardBlockAction(*M, *standard_M);
  CheckStandardBlockAction(*M_zero, *standard_M_zero);
  REQUIRE((C != nullptr) == lumped_port);
  REQUIRE((standard_C != nullptr) == lumped_port);
  if (C)
  {
    CheckStandardBlockAction(*C, *standard_C);
  }
  if (tetrahedral && !curved && loss_tangent == 0.0 && order == 1 && !heterogeneous &&
      !lumped_port && open_exterior)
  {
    CheckLosslessSpectralEmbedding(*standard_K_zero, *standard_M_zero,
                                   standard_space_op.GetNDDbcTDofLists().back(), *K_zero,
                                   *M_zero, space_op.GetCombinedNDDbcTDofList());
  }

  const auto nd_prolongations = space_op.GetCombinedNDProlongationOperators();
  const auto h1_prolongations = space_op.GetCombinedH1ProlongationOperators();
  const auto gradients = space_op.GetCombinedGradientOperators();
  const auto &nd_spaces = space_op.GetNDSpaces();
  const auto &h1_spaces = space_op.GetH1Spaces();
  const auto &nd_essential = space_op.GetCombinedNDDbcTDofLists();
  const auto &h1_essential = space_op.GetCombinedH1DbcTDofLists();
  const auto &h1_projector_essential = space_op.GetCombinedH1AuxBdrTDofList();
  REQUIRE(nd_spaces.GetNumLevels() == static_cast<std::size_t>(mg_max_levels));
  REQUIRE(h1_spaces.GetNumLevels() == nd_spaces.GetNumLevels());
  REQUIRE(nd_prolongations.size() + 1 == nd_spaces.GetNumLevels());
  REQUIRE(h1_prolongations.size() == nd_prolongations.size());
  REQUIRE(gradients.size() == nd_spaces.GetNumLevels());
  REQUIRE(nd_essential.size() == nd_spaces.GetNumLevels());
  REQUIRE(h1_essential.size() == h1_spaces.GetNumLevels());
  CHECK(std::is_sorted(h1_projector_essential.begin(), h1_projector_essential.end()));
  for (const int dof : h1_projector_essential)
  {
    CHECK(dof >= 0);
    CHECK(dof < gradients.back()->Width());
  }

  for (std::size_t level = 0; level < gradients.size(); level++)
  {
    const int standard_nd_size = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    const int standard_h1_size = h1_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    CHECK(gradients[level]->Height() - standard_nd_size ==
          hypre_G->Height() - space_op.GetNDSpace().GetTrueVSize());
    CHECK(gradients[level]->Width() - standard_h1_size ==
          hypre_G->Width() - space_op.GetH1Space().GetTrueVSize());
    CHECK(std::is_sorted(nd_essential[level].begin(), nd_essential[level].end()));
    CHECK(std::is_sorted(h1_essential[level].begin(), h1_essential[level].end()));
    for (const int dof : nd_essential[level])
    {
      CHECK(dof >= 0);
      CHECK(dof < gradients[level]->Height());
    }
    for (const int dof : h1_essential[level])
    {
      CHECK(dof >= 0);
      CHECK(dof < gradients[level]->Width());
    }
  }

  for (std::size_t level = 0; level < nd_prolongations.size(); level++)
  {
    const auto *nd_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations[level]);
    const auto *h1_prolongation =
        dynamic_cast<const mfem::HypreParMatrix *>(h1_prolongations[level]);
    REQUIRE(nd_prolongation);
    REQUIRE(h1_prolongation);
    CHECK(nd_prolongation->Width() == gradients[level]->Height());
    CHECK(nd_prolongation->Height() == gradients[level + 1]->Height());
    CHECK(h1_prolongation->Width() == gradients[level]->Width());
    CHECK(h1_prolongation->Height() == gradients[level + 1]->Width());

    const int coarse_standard_nd = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
    const int fine_standard_nd = nd_spaces.GetFESpaceAtLevel(level + 1).GetTrueVSize();
    Vector coarse_enrichment(nd_prolongation->Width()),
        fine_enrichment(nd_prolongation->Height());
    coarse_enrichment = 0.0;
    for (int i = coarse_standard_nd; i < coarse_enrichment.Size(); i++)
    {
      coarse_enrichment[i] = 0.13 * (i - coarse_standard_nd + 1);
    }
    nd_prolongation->Mult(coarse_enrichment, fine_enrichment);
    for (int i = 0; i < fine_standard_nd; i++)
    {
      CHECK(fine_enrichment[i] == 0.0);
    }
    for (int i = fine_standard_nd; i < fine_enrichment.Size(); i++)
    {
      CHECK(fine_enrichment[i] ==
            coarse_enrichment[coarse_standard_nd + i - fine_standard_nd]);
    }

    Vector coarse_h1(gradients[level]->Width()),
        prolonged_h1(gradients[level + 1]->Width()),
        gradient_after_prolongation(gradients[level + 1]->Height()),
        coarse_gradient(gradients[level]->Height()),
        prolonged_gradient(gradients[level + 1]->Height());
    FillVector(coarse_h1, 0.83);
    h1_prolongation->Mult(coarse_h1, prolonged_h1);
    gradients[level + 1]->Mult(prolonged_h1, gradient_after_prolongation);
    gradients[level]->Mult(coarse_h1, coarse_gradient);
    nd_prolongation->Mult(coarse_gradient, prolonged_gradient);
    gradient_after_prolongation -= prolonged_gradient;
    CHECK(RelativeNorm(gradient_after_prolongation, prolonged_gradient) < 2.0e-12);
  }

  if (mg_max_levels > 1)
  {
    auto preconditioner = space_op.GetPreconditionerMatrix<ComplexOperator>(
        std::complex<double>(1.0), std::complex<double>(0.0), std::complex<double>(-0.04),
        0.2);
    const auto *hierarchy =
        dynamic_cast<const BaseMultigridOperator<ComplexOperator> *>(preconditioner.get());
    REQUIRE(hierarchy);
    REQUIRE(hierarchy->GetNumLevels() == gradients.size());
    REQUIRE(hierarchy->GetNumAuxiliaryLevels() == gradients.size());
    auto exact_stiffness =
        space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto exact_mass = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto exact_system = space_op.GetSystemMatrix(
        std::complex<double>(1.0), std::complex<double>(0.0), std::complex<double>(-0.04),
        exact_stiffness.get(), static_cast<const ComplexOperator *>(nullptr),
        exact_mass.get());
    ComplexVector probe(exact_system->Width()), exact_action(exact_system->Height()),
        preconditioner_action(exact_system->Height());
    FillVector(probe.Real(), 0.37);
    FillVector(probe.Imag(), 0.79);
    linalg::SetSubVector(probe, nd_essential.back(), 0.0);
    exact_system->Mult(probe, exact_action);
    hierarchy->GetOperatorAtLevel(hierarchy->GetNumLevels() - 1)
        .Mult(probe, preconditioner_action);
    linalg::AXPY(-1.0, exact_action, preconditioner_action);
    CHECK(linalg::Norml2(Mpi::World(), preconditioner_action) <
          2.0e-11 * std::max(1.0, linalg::Norml2(Mpi::World(), exact_action)));
    for (std::size_t level = 0; level < gradients.size(); level++)
    {
      CHECK(hierarchy->GetOperatorAtLevel(level).Height() == gradients[level]->Height());
      CHECK(hierarchy->GetAuxiliaryOperatorAtLevel(level).Height() ==
            gradients[level]->Width());
      CHECK((hierarchy->GetOperatorAtLevel(level).Imag() != nullptr) ==
            (loss_tangent > 0.0));
      CHECK((hierarchy->GetAuxiliaryOperatorAtLevel(level).Imag() != nullptr) ==
            (loss_tangent > 0.0));
    }
    for (std::size_t level = 0; level + 1 < gradients.size(); level++)
    {
      const auto *coarse = hierarchy->GetOperatorAtLevel(level).Real();
      const auto *fine = hierarchy->GetOperatorAtLevel(level + 1).Real();
      const auto *coarse_imaginary = hierarchy->GetOperatorAtLevel(level).Imag();
      const auto *fine_imaginary = hierarchy->GetOperatorAtLevel(level + 1).Imag();
      const auto *prolongation =
          dynamic_cast<const mfem::HypreParMatrix *>(nd_prolongations[level]);
      REQUIRE(coarse);
      REQUIRE(fine);
      REQUIRE(prolongation);
      REQUIRE((coarse_imaginary != nullptr) == (loss_tangent > 0.0));
      REQUIRE((fine_imaginary != nullptr) == (loss_tangent > 0.0));

      const int standard_size = nd_spaces.GetFESpaceAtLevel(level).GetTrueVSize();
      const auto galerkin_error = [&](const Operator &coarse_operator,
                                      const Operator &fine_operator, bool include_standard,
                                      bool include_enrichment)
      {
        Vector coarse_x(coarse_operator.Width()), coarse_action(coarse_operator.Height()),
            fine_x(fine_operator.Width()), fine_action(fine_operator.Height()),
            restricted_action(coarse_operator.Height());
        FillVector(coarse_x, 0.41 + level);
        for (int i = 0; i < coarse_x.Size(); i++)
        {
          if ((i < standard_size && !include_standard) ||
              (i >= standard_size && !include_enrichment))
          {
            coarse_x[i] = 0.0;
          }
        }
        linalg::SetSubVector(coarse_x, nd_essential[level], 0.0);
        coarse_operator.Mult(coarse_x, coarse_action);
        prolongation->Mult(coarse_x, fine_x);
        fine_operator.Mult(fine_x, fine_action);
        prolongation->MultTranspose(fine_action, restricted_action);
        // The recursively projected level is re-eliminated after P^T A P. The
        // V-cycle likewise zeros coarse essential residuals before correction,
        // so compare the Galerkin identity only on the free subspace.
        linalg::SetSubVector(restricted_action, nd_essential[level], 0.0);
        restricted_action -= coarse_action;
        return RelativeNorm(restricted_action, coarse_action);
      };
      const double standard_error = galerkin_error(*coarse, *fine, true, false);
      const double enrichment_error = galerkin_error(*coarse, *fine, false, true);
      const double combined_error = galerkin_error(*coarse, *fine, true, true);
      const double imaginary_error =
          coarse_imaginary ? galerkin_error(*coarse_imaginary, *fine_imaginary, true, true)
                           : 0.0;
      CAPTURE(level, standard_error, enrichment_error, combined_error, imaginary_error);
      CHECK(combined_error < 2.0e-8);
      CHECK(imaginary_error < 2.0e-8);
    }
  }

  auto complex_mass =
      space_op.GetMassMatrix<ComplexOperator>(Operator::DiagonalPolicy::DIAG_ZERO);
  auto complex_stiffness =
      space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DiagonalPolicy::DIAG_ONE);
  auto complex_damping =
      space_op.GetDampingMatrix<ComplexOperator>(Operator::DiagonalPolicy::DIAG_ZERO);
  auto standard_complex_mass =
      standard_space_op.GetMassMatrix<ComplexOperator>(Operator::DiagonalPolicy::DIAG_ZERO);
  auto standard_complex_stiffness = standard_space_op.GetStiffnessMatrix<ComplexOperator>(
      Operator::DiagonalPolicy::DIAG_ONE);
  auto standard_complex_damping = standard_space_op.GetDampingMatrix<ComplexOperator>(
      Operator::DiagonalPolicy::DIAG_ZERO);
  REQUIRE(complex_stiffness);
  REQUIRE(standard_complex_stiffness);
  REQUIRE((complex_damping != nullptr) == lumped_port);
  REQUIRE((standard_complex_damping != nullptr) == lumped_port);
  const auto *complex_mass_real = complex_mass->Real();
  const auto *complex_mass_imag = complex_mass->Imag();
  REQUIRE(complex_mass_real);
  CHECK((complex_mass_imag != nullptr) == (loss_tangent > 0.0));
  REQUIRE(standard_complex_mass->Real());
  CheckStandardBlockAction(*complex_mass_real, *standard_complex_mass->Real());
  CHECK((standard_complex_mass->Imag() != nullptr) == (complex_mass_imag != nullptr));
  if (complex_mass_imag)
  {
    CheckStandardBlockAction(*complex_mass_imag, *standard_complex_mass->Imag());
    Vector probe(complex_mass->Width()), real_action(complex_mass->Height()),
        imaginary_action(complex_mass->Height()), bulk_action(complex_mass->Height());
    FillVector(probe, 0.61);
    complex_mass_real->Mult(probe, real_action);
    complex_mass_imag->Mult(probe, imaginary_action);
    M_bulk->Mult(probe, bulk_action);
    Vector loss_error(imaginary_action);
    loss_error.Add(loss_tangent, bulk_action);
    CHECK(RelativeNorm(loss_error, bulk_action) < 2.0e-11);

    constexpr std::complex<double> coefficient(1.3, -0.4);
    auto system = space_op.GetSystemMatrix(
        std::complex<double>(0.0), std::complex<double>(0.0), coefficient,
        static_cast<const ComplexOperator *>(nullptr),
        static_cast<const ComplexOperator *>(nullptr), complex_mass.get());
    ComplexVector system_probe(system->Width()), system_action(system->Height());
    system_probe = 0.0;
    system_probe.Real() = probe;
    system->Mult(system_probe, system_action);
    Vector expected_real(real_action), expected_imag(real_action);
    expected_real *= coefficient.real();
    expected_real.Add(-coefficient.imag(), imaginary_action);
    expected_imag *= coefficient.imag();
    expected_imag.Add(coefficient.real(), imaginary_action);
    system_action.Real() -= expected_real;
    system_action.Imag() -= expected_imag;
    CHECK(RelativeNorm(system_action.Real(), expected_real) < 2.0e-11);
    CHECK(RelativeNorm(system_action.Imag(), expected_imag) < 2.0e-11);
  }

  constexpr std::complex<double> lambda(-0.021, 0.37);
  auto polynomial = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), lambda,
                                             lambda * lambda, complex_stiffness.get(),
                                             complex_damping.get(), complex_mass.get());
  auto standard_polynomial = standard_space_op.GetSystemMatrix(
      std::complex<double>(1.0, 0.0), lambda, lambda * lambda,
      standard_complex_stiffness.get(), standard_complex_damping.get(),
      standard_complex_mass.get());
  CheckStandardBlockAction(*polynomial, *standard_polynomial);

  Vector h1(hypre_G->Width()), gradient(hypre_G->Height()), curl_gradient(K_zero->Height());
  FillVector(h1, 0.4);
  linalg::SetSubVector(h1, space_op.GetCombinedH1DbcTDofList(), 0.0);
  hypre_G->Mult(h1, gradient);
  K_bulk->Mult(gradient, curl_gradient);
  CHECK(RelativeNorm(curl_gradient, gradient) < 5.0e-8);

  Vector nd(K->Width()), Knd(K->Height());
  FillVector(nd, 1.1);
  K->Mult(nd, Knd);
  for (const int dof : space_op.GetCombinedNDDbcTDofList())
  {
    CHECK(Knd[dof] == nd[dof]);
  }

  // The projector must act on both standard and enrichment coordinates.
  auto scalar_diffusion = space_op.GetBulkScalarDiffusionMatrix();
  Vector scalar(hypre_G->Width()), gradient_mass(hypre_G->Height()),
      algebraic_diffusion(hypre_G->Width()), assembled_diffusion(hypre_G->Width());
  FillVector(scalar, 0.65);
  linalg::SetSubVector(scalar, h1_projector_essential, 0.0);
  hypre_G->Mult(scalar, gradient);
  M_bulk->Mult(gradient, gradient_mass);
  hypre_G->MultTranspose(gradient_mass, algebraic_diffusion);
  scalar_diffusion->Mult(scalar, assembled_diffusion);
  linalg::SetSubVector(algebraic_diffusion, h1_projector_essential, 0.0);
  linalg::SetSubVector(assembled_diffusion, h1_projector_essential, 0.0);
  linalg::AXPY(-1.0, algebraic_diffusion, assembled_diffusion);
  CHECK(linalg::Norml2(Mpi::World(), assembled_diffusion) <=
        2.0e-10 * std::max(1.0, linalg::Norml2(Mpi::World(), algebraic_diffusion)));

  DivFreeSolver<ComplexVector> divfree(iodata, Mpi::World(), *M_bulk, *hypre_G,
                                       std::move(scalar_diffusion), h1_projector_essential);
  ComplexVector electric(M_zero->Width()), projected(M_zero->Width()),
      mass_electric(M_zero->Height());
  FillVector(electric.Real(), 0.25);
  FillVector(electric.Imag(), 0.75);
  linalg::SetSubVector(electric, space_op.GetCombinedNDDbcTDofList(), 0.0);
  projected = electric;
  divfree.Mult(projected);

  ComplexVector weak_divergence(hypre_G->Width());
  M_bulk->Mult(projected.Real(), mass_electric.Real());
  M_bulk->Mult(projected.Imag(), mass_electric.Imag());
  hypre_G->MultTranspose(mass_electric.Real(), weak_divergence.Real());
  hypre_G->MultTranspose(mass_electric.Imag(), weak_divergence.Imag());
  linalg::SetSubVector(weak_divergence, h1_projector_essential, 0.0);
  CHECK(linalg::Norml2(Mpi::World(), weak_divergence) <=
        2.0e-10 * std::max(1.0, linalg::Norml2(Mpi::World(), mass_electric)));

  ComplexVector projected_twice = projected;
  divfree.Mult(projected_twice);
  linalg::AXPY(-1.0, projected, projected_twice);
  CHECK(linalg::Norml2(Mpi::World(), projected_twice) <=
        2.0e-10 * std::max(1.0, linalg::Norml2(Mpi::World(), projected)));
}

}  // namespace

TEST_CASE("Parallel hierarchical Maxwell p-plus-one residual lifting",
          "[spaceoperator][singularelements][hierarchical][Serial][Parallel]")
{
  REQUIRE((Mpi::Size(Mpi::World()) == 1 || Mpi::Size(Mpi::World()) == 2));
  auto serial_mesh = InternalSheetMesh();
  const auto serial_features = fem::singular::ExtractSerialSheetFeatures(serial_mesh, {7});
  std::vector<fem::singular::GlobalVertexId> serial_vertex_ids(serial_mesh.GetNV());
  std::iota(serial_vertex_ids.begin(), serial_vertex_ids.end(), 0);
  const std::array<int, 4> serial_partition{0, 0, 0, 0};
  const std::array<int, 4> parallel_partition{0, 0, 1, 1};
  const int *partition =
      Mpi::Size(Mpi::World()) == 1 ? serial_partition.data() : parallel_partition.data();
  auto parallel_mesh =
      std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh, partition);
  const auto local_vertex_ids =
      fem::singular::MapPartitionedSerialVertexIds(serial_mesh, *parallel_mesh, partition);
  const auto local_features = fem::singular::DistributeSerialSheetFeatures(
      serial_mesh, serial_features, *parallel_mesh, local_vertex_ids);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::move(parallel_mesh)));
  auto iodata = SingularSpaceData(3, 1, 1, false, true);
  for (auto &material : iodata.domains.materials)
  {
    material.tandelta.s = {0.017, 0.017, 0.017};
  }
  SpaceOperator space_op(iodata, meshes, &local_features, nullptr, &local_vertex_ids);
  HierarchicalMaxwellDomainData hierarchy(space_op);

  // A projected smooth physical field is decomposition-independent, unlike raw true-DOF
  // coefficient patterns whose numbering and edge orientations depend on the partition.
  const auto build_field = [](SpaceOperator &space_operator)
  {
    ComplexVector field(space_operator.GetNDTrueVSize());
    field = 0.0;
    mfem::VectorFunctionCoefficient real_field(
        3,
        [](const mfem::Vector &x, mfem::Vector &value)
        {
          value(0) = x(1) * x(2) + 0.3 * x(0);
          value(1) = x(0) - 0.5 * x(2) * x(2);
          value(2) = x(0) * x(1) + 0.25 * x(2);
        });
    mfem::VectorFunctionCoefficient imag_field(
        3,
        [](const mfem::Vector &x, mfem::Vector &value)
        {
          value(0) = x(2) - 0.2 * x(1);
          value(1) = x(0) * x(1);
          value(2) = -x(0) + 0.1 * x(2) * x(2);
        });
    mfem::ParGridFunction projected(&space_operator.GetNDSpace().Get());
    Vector standard_true(space_operator.GetNDSpace().GetTrueVSize());
    projected.ProjectCoefficient(real_field);
    projected.GetTrueDofs(standard_true);
    for (int dof = 0; dof < standard_true.Size(); dof++)
    {
      field.Real()(dof) = standard_true(dof);
    }
    projected.ProjectCoefficient(imag_field);
    projected.GetTrueDofs(standard_true);
    for (int dof = 0; dof < standard_true.Size(); dof++)
    {
      field.Imag()(dof) = standard_true(dof);
    }
    for (int dof : space_operator.GetCombinedNDDbcTDofList())
    {
      field.Real()(dof) = 0.0;
      field.Imag()(dof) = 0.0;
    }
    return field;
  };
  const auto field = build_field(space_op);
  const auto estimate = hierarchy.EstimatePolynomialEigenResidual({0.7, 0.08}, field);
  CAPTURE(estimate.total_energy);
  CHECK(std::isfinite(estimate.total_energy));
  CHECK(estimate.total_energy > 0.0);
  double sum = estimate.indicator_energy.Sum();
  Mpi::GlobalSum(1, &sum, Mpi::World());
  CHECK_THAT(sum, Catch::Matchers::WithinRel(estimate.total_energy, 1.0e-12));
  for (int element = 0; element < estimate.indicator_energy.Size(); element++)
  {
    CHECK(estimate.indicator_energy(element) >= -1.0e-12 * estimate.total_energy);
  }

  if (Mpi::Size(Mpi::World()) == 2)
  {
    // Replicated one-rank oracle: every rank rebuilds the complete problem on its own
    // communicator and the distributed estimate must reproduce it to roundoff, including
    // the per-element localization under the known partition.
    auto replicated_serial_mesh = InternalSheetMesh();
    auto replicated_mesh = std::make_unique<mfem::ParMesh>(
        MPI_COMM_SELF, replicated_serial_mesh, serial_partition.data());
    const auto replicated_vertex_ids = fem::singular::MapPartitionedSerialVertexIds(
        replicated_serial_mesh, *replicated_mesh, serial_partition.data());
    const auto replicated_features = fem::singular::DistributeSerialSheetFeatures(
        replicated_serial_mesh, serial_features, *replicated_mesh, replicated_vertex_ids);
    std::vector<std::unique_ptr<Mesh>> replicated_meshes;
    replicated_meshes.push_back(std::make_unique<Mesh>(std::move(replicated_mesh)));
    auto replicated_iodata = SingularSpaceData(3, 1, 1, false, true);
    for (auto &material : replicated_iodata.domains.materials)
    {
      material.tandelta.s = {0.017, 0.017, 0.017};
    }
    SpaceOperator replicated_space_op(replicated_iodata, replicated_meshes,
                                      &replicated_features, nullptr,
                                      &replicated_vertex_ids);
    HierarchicalMaxwellDomainData replicated_hierarchy(replicated_space_op);
    const auto replicated_field = build_field(replicated_space_op);
    const auto replicated_estimate =
        replicated_hierarchy.EstimatePolynomialEigenResidual({0.7, 0.08}, replicated_field);
    CAPTURE(replicated_estimate.total_energy);
    REQUIRE(replicated_estimate.total_energy > 0.0);
    CHECK_THAT(estimate.total_energy,
               Catch::Matchers::WithinRel(replicated_estimate.total_energy, 1.0e-10));
    // Partition-ordered element mapping: rank 0 owns serial elements {0,1}, rank 1 owns
    // {2,3}.
    REQUIRE(estimate.indicator_energy.Size() == 2);
    const int offset = Mpi::Rank(Mpi::World()) == 0 ? 0 : 2;
    for (int element = 0; element < estimate.indicator_energy.Size(); element++)
    {
      CAPTURE(element, offset);
      CHECK_THAT(estimate.indicator_energy(element),
                 Catch::Matchers::WithinRel(
                     replicated_estimate.indicator_energy(offset + element), 1.0e-9));
    }
  }
}

TEST_CASE("Full-wave singular SpaceOperator preserves Maxwell algebra on high-order maps",
          "[spaceoperator][singularelements][curved][Serial]")
{
  SECTION("2D quadratic affine")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false);
  }
  SECTION("2D genuinely curved")
  {
    CheckSpaceOperator(InternalLineTipMesh(), true);
  }
  SECTION("3D quadratic affine")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017, 1, 1, false, true, false, true);
  }
  SECTION("3D lossless spectral embedding")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.0, 1, 1, false, false, true);
  }
  SECTION("3D genuinely curved")
  {
    CheckSpaceOperator(InternalSheetMesh(), true);
  }
  SECTION("2D isotropic dielectric loss")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false, 0.017);
  }
  SECTION("3D isotropic dielectric loss")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017);
  }
  SECTION("3D resistive and reactive lumped port")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017, 2, 1, false, true);
  }
  SECTION("3D heterogeneous dielectric interface")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017, 1, 1, true);
  }
  SECTION("2D p = 2 combined hierarchy")
  {
    CheckSpaceOperator(InternalLineTipMesh(), false, 0.017, 2, 2);
  }
  SECTION("3D p = 2 combined hierarchy")
  {
    CheckSpaceOperator(InternalSheetMesh(), false, 0.017, 2, 2);
  }
}

}  // namespace palace
