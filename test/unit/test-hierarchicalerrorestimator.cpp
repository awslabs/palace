// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <numeric>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/hierarchicalerrorestimator.hpp"

namespace palace
{
using namespace Catch::Matchers;

TEST_CASE("Hierarchical patch algebra includes domain and boundary contributions",
          "[hierarchicalerrorestimator][Serial]")
{
  using fem::hierarchical::LocalOperatorContribution;
  using fem::hierarchical::SparseColumn;
  LocalOperatorContribution domain, boundary, exterior;
  domain.support_element = boundary.support_element = 0;
  exterior.support_element = 1;
  domain.dofs = boundary.dofs = exterior.dofs = {0, 1};
  domain.matrix.SetSize(2);
  boundary.matrix.SetSize(2);
  exterior.matrix.SetSize(2);
  domain.matrix(0, 0) = 3.0;
  domain.matrix(0, 1) = domain.matrix(1, 0) = -1.0;
  domain.matrix(1, 1) = 2.0;
  boundary.matrix = 0.0;
  boundary.matrix(0, 0) = 0.5;
  boundary.matrix(1, 1) = 1.5;
  exterior.matrix = 0.0;
  exterior.matrix(0, 0) = 100.0;
  exterior.matrix(1, 1) = 100.0;
  domain.rhs.SetSize(2);
  boundary.rhs.SetSize(2);
  exterior.rhs.SetSize(2);
  domain.rhs(0) = 2.0;
  domain.rhs(1) = 1.0;
  boundary.rhs(0) = 0.25;
  boundary.rhs(1) = -0.5;
  exterior.rhs = 0.0;
  const std::vector<LocalOperatorContribution> contributions{domain, boundary, exterior};

  mfem::Vector injected(2);
  injected(0) = 0.2;
  injected(1) = -0.1;
  const auto residual =
      fem::hierarchical::AssembleResidual(2, contributions, injected, {false, false});
  // AssembleResidual includes every contribution, including the exterior element.
  CHECK_THAT(residual(0), WithinAbs(2.25 - 3.5 * 0.2 + 1.0 * -0.1 - 20.0, 1.0e-14));
  CHECK_THAT(residual(1), WithinAbs(0.5 + 1.0 * 0.2 - 3.5 * -0.1 + 10.0, 1.0e-14));

  SparseColumn first{{0}, {1.0}}, second{{0, 1}, {0.5, 1.0}};
  mfem::DenseMatrix restricted;
  mfem::Vector restricted_rhs;
  const long long touched = fem::hierarchical::AssembleRestrictedOperator(
      contributions, {0}, {first, second}, residual, restricted, restricted_rhs);
  // The support filter excludes exterior while retaining both the domain and boundary
  // facet.
  mfem::DenseMatrix expected(2);
  expected(0, 0) = 3.5;
  expected(0, 1) = expected(1, 0) = 0.75;
  expected(1, 1) = 3.375;
  CHECK(touched == 8);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK_THAT(restricted(i, j), WithinAbs(expected(i, j), 1.0e-14));
    }
  }

  mfem::Vector vector(2);
  vector(0) = -0.3;
  vector(1) = 0.7;
  const double energy = fem::hierarchical::Energy({domain, boundary}, vector);
  CHECK_THAT(energy, WithinAbs(2.45, 1.0e-14));
}

TEST_CASE("Hierarchical complex residual retains loss damping and port terms",
          "[hierarchicalerrorestimator][Serial]")
{
  using fem::hierarchical::ComplexLocalOperatorContribution;
  constexpr double omega = 1.7;
  const auto contribution =
      [](const mfem::DenseMatrix &real, const mfem::DenseMatrix &imag, bool source = false)
  {
    ComplexLocalOperatorContribution data;
    data.support_element = 0;
    data.dofs = {0, 1};
    data.matrix_real = real;
    data.matrix_imag = imag;
    data.rhs_real.SetSize(2);
    data.rhs_imag.SetSize(2);
    data.rhs_real = 0.0;
    data.rhs_imag = 0.0;
    if (source)
    {
      data.rhs_real(0) = 1.0;
      data.rhs_real(1) = -0.2;
      data.rhs_imag(0) = 0.4;
      data.rhs_imag(1) = 0.7;
    }
    return data;
  };
  const auto matrix = [](double a00, double a01, double a11)
  {
    mfem::DenseMatrix value(2);
    value(0, 0) = a00;
    value(0, 1) = value(1, 0) = a01;
    value(1, 1) = a11;
    return value;
  };
  const mfem::DenseMatrix zero = matrix(0.0, 0.0, 0.0);
  const mfem::DenseMatrix K = matrix(4.0, -1.0, 3.0);
  const mfem::DenseMatrix C = matrix(0.8, 0.1, 0.5);
  const mfem::DenseMatrix Mr = matrix(2.0, 0.2, 1.0);
  const mfem::DenseMatrix Mi = matrix(0.15, 0.0, 0.08);
  const mfem::DenseMatrix A2r = matrix(0.3, -0.05, 0.2);
  const mfem::DenseMatrix A2i = matrix(0.12, 0.02, 0.18);
  mfem::DenseMatrix omega_C(C), minus_omega2_Mr(Mr), minus_omega2_Mi(Mi);
  omega_C *= omega;
  minus_omega2_Mr *= -omega * omega;
  minus_omega2_Mi *= -omega * omega;
  const std::vector<ComplexLocalOperatorContribution> contributions{
      contribution(K, zero, true),          // K and b
      contribution(zero, omega_C),          // i omega C (resistive port)
      contribution(minus_omega2_Mr, zero),  // -omega^2 Mr
      contribution(zero, minus_omega2_Mi),  // -i omega^2 Mi (dielectric loss)
      contribution(A2r, A2i)};              // general complex A2(omega)

  mfem::Vector xr(2), xi(2);
  xr(0) = 0.25;
  xr(1) = -0.4;
  xi(0) = -0.3;
  xi(1) = 0.2;
  const auto residual =
      fem::hierarchical::AssembleComplexResidual(2, contributions, xr, xi, {false, false});

  // Independent complex arithmetic oracle for
  // A = K + i omega C - omega^2 (Mr + i Mi) + A2(omega).
  mfem::DenseMatrix Ar(K), Ai(omega_C);
  Ar += minus_omega2_Mr;
  Ar += A2r;
  Ai += minus_omega2_Mi;
  Ai += A2i;
  const std::complex<double> b[2]{{1.0, 0.4}, {-0.2, 0.7}};
  const std::complex<double> x[2]{{xr(0), xi(0)}, {xr(1), xi(1)}};
  for (int row = 0; row < 2; row++)
  {
    std::complex<double> action = 0.0;
    for (int column = 0; column < 2; column++)
    {
      action += std::complex<double>(Ar(row, column), Ai(row, column)) * x[column];
    }
    const auto expected = b[row] - action;
    CHECK_THAT(residual.real(row), WithinAbs(expected.real(), 1.0e-14));
    CHECK_THAT(residual.imag(row), WithinAbs(expected.imag(), 1.0e-14));
  }

  // Each non-Hermitian physical family is active: dropping damping, dielectric loss, or
  // A2 changes the exact complex residual instead of silently reducing to K-omega^2 M.
  for (int removed : {1, 3, 4})
  {
    auto ablated = contributions;
    ablated.erase(ablated.begin() + removed);
    const auto residual_ablated =
        fem::hierarchical::AssembleComplexResidual(2, ablated, xr, xi, {false, false});
    const double difference = std::hypot(residual.real(0) - residual_ablated.real(0),
                                         residual.imag(0) - residual_ablated.imag(0)) +
                              std::hypot(residual.real(1) - residual_ablated.real(1),
                                         residual.imag(1) - residual_ablated.imag(1));
    CAPTURE(removed, difference);
    CHECK(difference > 1.0e-2);
  }

  const auto constrained =
      fem::hierarchical::AssembleComplexResidual(2, contributions, xr, xi, {false, true});
  CHECK_THAT(constrained.real(0), WithinAbs(residual.real(0), 1.0e-14));
  CHECK_THAT(constrained.imag(0), WithinAbs(residual.imag(0), 1.0e-14));
  CHECK(constrained.real(1) == 0.0);
  CHECK(constrained.imag(1) == 0.0);
}

TEST_CASE("Hierarchical complex lifting is invariant under global phase",
          "[hierarchicalerrorestimator][Serial]")
{
  auto mesh = mfem::Mesh::MakeCartesian2D(2, 1, mfem::Element::TRIANGLE, true, 2.0, 1.0);
  mfem::ND_FECollection coarse_collection(1, 2), fine_collection(2, 2);
  mfem::FiniteElementSpace coarse_space(&mesh, &coarse_collection);
  mfem::FiniteElementSpace fine_space(&mesh, &fine_collection);
  const auto injection =
      fem::hierarchical::BuildSparsePInjection(mesh, coarse_space, fine_space);

  // Positive H(curl) graph metric in the fine space, scattered in unsigned signed-VDof
  // coordinates exactly as the production element path does.
  std::vector<fem::hierarchical::LocalOperatorContribution> metric(mesh.GetNE());
  mfem::CurlCurlIntegrator curl_curl;
  mfem::VectorFEMassIntegrator mass;
  mfem::Array<int> dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    auto &data = metric[element];
    data.support_element = element;
    mfem::DofTransformation transformation;
    fine_space.GetElementVDofs(element, dofs, transformation);
    for (int dof : dofs)
    {
      data.dofs.push_back(dof >= 0 ? dof : -1 - dof);
    }
    mfem::DenseMatrix local_curl, local_mass;
    auto &T = *mesh.GetElementTransformation(element);
    const auto &fe = *fine_space.GetFE(element);
    curl_curl.AssembleElementMatrix(fe, T, local_curl);
    mass.AssembleElementMatrix(fe, T, local_mass);
    local_curl += local_mass;
    transformation.TransformDual(local_curl);
    data.matrix.SetSize(dofs.Size());
    data.rhs.SetSize(dofs.Size());
    data.rhs = 0.0;
    for (int row = 0; row < dofs.Size(); row++)
    {
      const double row_sign = dofs[row] >= 0 ? 1.0 : -1.0;
      for (int column = 0; column < dofs.Size(); column++)
      {
        const double column_sign = dofs[column] >= 0 ? 1.0 : -1.0;
        data.matrix(row, column) = row_sign * column_sign * local_curl(row, column);
      }
    }
  }

  const int fine_size = fine_space.GetVSize();
  fem::hierarchical::ComplexResidual residual{mfem::Vector(fine_size),
                                              mfem::Vector(fine_size)};
  for (int dof = 0; dof < fine_size; dof++)
  {
    residual.real(dof) = std::sin(0.31 * dof + 0.2);
    residual.imag(dof) = std::cos(0.47 * dof - 0.1);
  }
  const std::vector<bool> fine_essential(fine_size, false);
  const std::vector<bool> coarse_essential(coarse_space.GetVSize(), false);
  const std::vector<std::vector<int>> no_enrichment(mesh.GetNE());
  fem::hierarchical::PatchLiftingOptions options;
  options.sweeps = 3;
  const auto reference = fem::hierarchical::LiftComplexResidualByPatches(
      mesh, coarse_space, fine_space, injection, metric, fine_essential, coarse_essential,
      residual, no_enrichment, options);

  constexpr double phase = 0.63;
  fem::hierarchical::ComplexResidual rotated{mfem::Vector(fine_size),
                                             mfem::Vector(fine_size)};
  for (int dof = 0; dof < fine_size; dof++)
  {
    rotated.real(dof) =
        std::cos(phase) * residual.real(dof) - std::sin(phase) * residual.imag(dof);
    rotated.imag(dof) =
        std::sin(phase) * residual.real(dof) + std::cos(phase) * residual.imag(dof);
  }
  const auto transformed = fem::hierarchical::LiftComplexResidualByPatches(
      mesh, coarse_space, fine_space, injection, metric, fine_essential, coarse_essential,
      rotated, no_enrichment, options);
  REQUIRE(reference.energy > 0.0);
  CHECK_THAT(transformed.energy, WithinRel(reference.energy, 1.0e-12));
  REQUIRE(transformed.indicator.size() == reference.indicator.size());
  for (std::size_t element = 0; element < reference.indicator.size(); element++)
  {
    CHECK_THAT(
        transformed.indicator[element],
        WithinAbs(reference.indicator[element], 1.0e-12 * std::max(reference.energy, 1.0)));
  }
  CHECK_THAT(std::accumulate(reference.indicator.begin(), reference.indicator.end(), 0.0),
             WithinRel(reference.energy, 1.0e-12));
}

}  // namespace palace
