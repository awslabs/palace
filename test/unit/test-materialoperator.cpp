// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <limits>
#include <list>
#include <numbers>
#include <utility>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/operator.hpp"
#include "fixtures.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/hypre.hpp"
#include "linalg/nleps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/domainpostoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/constants.hpp"
#include "utils/geodata.hpp"
#include "utils/permittivitymodel.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch;
using namespace std::complex_literals;

namespace
{

// Two unit cubes of tetrahedra on each process, along x. Attribute 2 is confined to rank
// 0's first cube, separated from the rank interface by attribute 1 (the rest of the mesh):
// with more than one process, the others have no local support of a material on
// attribute 2.
std::vector<std::unique_ptr<Mesh>> MakePartitionedSupportMesh(MPI_Comm comm)
{
  const int size = Mpi::Size(comm);
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      2 * size, 1, 1, mfem::Element::TETRAHEDRON, 2.0 * size, 1.0, 1.0));
  std::vector<int> partitioning(serial_mesh->GetNE());
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    const auto *element = serial_mesh->GetElement(i);
    const auto *vertices = element->GetVertices();
    double center_x = 0.0;
    for (int j = 0; j < element->GetNVertices(); j++)
    {
      center_x += serial_mesh->GetVertex(vertices[j])[0];
    }
    center_x /= element->GetNVertices();
    serial_mesh->SetAttribute(i, center_x < 1.0 ? 2 : 1);
    partitioning[i] = static_cast<int>(center_x / 2.0);
  }
  serial_mesh->SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning.data());
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  return mesh;
}

// Number of elements of the libCEED sub-operators of a partially assembled operator on this
// process.
CeedInt NumCeedOperatorElements(const Operator &op)
{
  const auto *ceed_op = dynamic_cast<const ceed::Operator *>(&op);
  REQUIRE(ceed_op);
  CeedInt num_elem = 0;
  for (std::size_t i = 0; i < ceed_op->Size(); i++)
  {
    Ceed ceed;
    CeedInt num_sub_ops;
    CeedOperator *sub_ops;
    PalaceCeedCallBackend(CeedOperatorGetCeed((*ceed_op)[i], &ceed));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub((*ceed_op)[i], &num_sub_ops));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList((*ceed_op)[i], &sub_ops));
    for (CeedInt j = 0; j < num_sub_ops; j++)
    {
      CeedInt n;
      PalaceCeedCall(ceed, CeedOperatorGetNumElements(sub_ops[j], &n));
      num_elem += n;
    }
  }
  return num_elem;
}

}  // namespace

TEST_CASE("MaterialOperator IsIsotropic", "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  config::PeriodicBoundaryData periodic;

  SECTION("Trivial isotropic material")
  {
    config::MaterialData material;
    material.attributes = {1};
    // Default values should be isotropic (all eigenvalues = 1).

    MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Non-trivial isotropic material")
  {
    config::MaterialData material;
    material.attributes = {1};
    material.mu_r.s[0] = 2.0;
    material.mu_r.s[1] = 2.0;
    material.mu_r.s[2] = 2.0;

    MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC, palace_mesh);
    REQUIRE(mat_op.IsIsotropic(1) == true);
  }

  SECTION("Anisotropic materials")
  {
    config::MaterialData material;
    material.attributes = {1};

    SECTION("Anisotropic permeability")
    {
      material.mu_r.s[0] = 2.0;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic permittivity")
    {
      material.epsilon_r.s[1] = 2.0;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic loss tangent")
    {
      material.tandelta.s[2] = 0.02;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }

    SECTION("Anisotropic conductivity")
    {
      material.sigma.s[0] = 1e6;
      material.sigma.s[2] = 2e6;
      MaterialOperator mat_op({material}, periodic, ProblemType::ELECTROSTATIC,
                              palace_mesh);
      REQUIRE(mat_op.IsIsotropic(1) == false);
    }
  }
}

TEST_CASE("MaterialPropertyCoefficient exact-zero predicate", "[materialoperator][Serial]")
{
  mfem::Array<int> attr_mat(2);
  attr_mat = 0;

  SECTION("Empty storage is mathematically zero")
  {
    MaterialPropertyCoefficient coeff(2);
    CHECK(coeff.empty());
    CHECK(coeff.IsExactlyZero());
  }

  SECTION("Allocated scalar and tensor zeros are exactly zero")
  {
    mfem::DenseTensor scalar(1, 1, 1);
    scalar = 0.0;
    MaterialPropertyCoefficient scalar_coeff(attr_mat, scalar);
    CHECK_FALSE(scalar_coeff.empty());
    CHECK(scalar_coeff.IsExactlyZero());

    mfem::DenseTensor tensor(2, 2, 1);
    tensor = 0.0;
    MaterialPropertyCoefficient tensor_coeff(attr_mat, tensor);
    CHECK_FALSE(tensor_coeff.empty());
    CHECK(tensor_coeff.IsExactlyZero());
  }

  SECTION("Zero-scaled coefficient additions retain exact-zero storage")
  {
    mfem::DenseTensor tensor(2, 2, 1);
    tensor = 0.0;
    tensor(0, 0, 0) = 3.0;
    tensor(1, 1, 0) = -2.0;

    MaterialPropertyCoefficient from_coefficient(2);
    from_coefficient.AddCoefficient(attr_mat, tensor, 0.0);
    CHECK_FALSE(from_coefficient.empty());
    CHECK(from_coefficient.IsExactlyZero());

    MaterialPropertyCoefficient from_property(2);
    from_property.AddMaterialProperty(1, 7.0, 0.0);
    CHECK_FALSE(from_property.empty());
    CHECK(from_property.IsExactlyZero());
  }

  SECTION("Mixed and tiny nonzero values are not exactly zero")
  {
    mfem::DenseTensor mixed(2, 2, 1);
    mixed = 0.0;
    mixed(0, 1, 0) = -4.0;
    MaterialPropertyCoefficient mixed_coeff(attr_mat, mixed);
    CHECK_FALSE(mixed_coeff.IsExactlyZero());

    mfem::DenseTensor tiny(1, 1, 1);
    tiny = std::numeric_limits<double>::min();
    MaterialPropertyCoefficient tiny_coeff(attr_mat, tiny);
    CHECK_FALSE(tiny_coeff.IsExactlyZero());
  }
}

TEST_CASE("MaterialPropertyCoefficient adds a term only to the listed attributes",
          "[materialoperator][Serial]")
{
  using Catch::Matchers::WithinAbs;
  using Catch::Matchers::WithinRel;
  // Attributes with equal properties share one material entry. A later term added to some
  // of the attributes of a shared entry must not reach the others, and a term stamped once
  // per attribute (as the boundary operators do per port element or attribute) must be
  // counted once on each.
  auto Value = [](const MaterialPropertyCoefficient &coeff, int attr)
  {
    const int k = coeff.GetAttributeToMaterial()[attr - 1];
    REQUIRE(k >= 0);
    return coeff.GetMaterialProperties()(0, 0, k);
  };
  auto Attributes = [](std::initializer_list<int> attrs)
  {
    mfem::Array<int> attr_list;
    for (auto attr : attrs)
    {
      attr_list.Append(attr);
    }
    return attr_list;
  };

  SECTION("Per-attribute stamping of a second term")
  {
    MaterialPropertyCoefficient fb(3);
    for (int attr = 1; attr <= 3; attr++)
    {
      fb.AddMaterialProperty(attr, 2.0);
    }
    CHECK(fb.GetMaterialProperties().SizeK() == 1);  // Equal values share one entry
    for (int attr = 1; attr <= 3; attr++)
    {
      fb.AddMaterialProperty(attr, 0.25);
    }
    for (int attr = 1; attr <= 3; attr++)
    {
      CHECK_THAT(Value(fb, attr), WithinRel(2.25, 1.0e-14));
    }
  }

  SECTION("Second term on a subset of a shared entry")
  {
    MaterialPropertyCoefficient fb(3);
    fb.AddMaterialProperty(Attributes({1, 2}), 2.0);
    fb.AddMaterialProperty(3, 2.0);
    CHECK(fb.GetMaterialProperties().SizeK() == 1);
    fb.AddMaterialProperty(3, 0.25);
    CHECK_THAT(Value(fb, 1), WithinRel(2.0, 1.0e-14));
    CHECK_THAT(Value(fb, 2), WithinRel(2.0, 1.0e-14));
    CHECK_THAT(Value(fb, 3), WithinRel(2.25, 1.0e-14));
  }

  SECTION("Sole user of an entry is updated in place")
  {
    MaterialPropertyCoefficient fb(2);
    fb.AddMaterialProperty(Attributes({1, 2}), 1.0);
    fb.AddMaterialProperty(Attributes({1, 2}), 2.0);
    CHECK(fb.GetMaterialProperties().SizeK() == 1);
    CHECK_THAT(Value(fb, 1), WithinRel(3.0, 1.0e-14));
    CHECK_THAT(Value(fb, 2), WithinRel(3.0, 1.0e-14));
  }

  SECTION("Mixed assigned and unassigned attributes are rejected")
  {
    MaterialPropertyCoefficient fb(3);
    fb.AddMaterialProperty(Attributes({1, 2}), 1.0);
    CHECK_THROWS(fb.AddMaterialProperty(Attributes({3, 1}), 0.5));
  }

  SECTION("Tensor properties through AddCoefficient")
  {
    // All attributes share one entry; a second coefficient with a different map adds to
    // attribute 1 alone and to attributes 2 and 3.
    mfem::Array<int> attr_mat_a(3), attr_mat_b(3);
    attr_mat_a = 0;
    attr_mat_b[0] = 0;
    attr_mat_b[1] = 1;
    attr_mat_b[2] = 1;
    mfem::DenseTensor tensor_a(2, 2, 1), tensor_b(2, 2, 2);
    tensor_a = 0.0;
    tensor_b = 0.0;
    tensor_a(0, 0, 0) = tensor_a(1, 1, 0) = 1.0;
    tensor_b(0, 1, 0) = tensor_b(1, 0, 0) = 0.5;
    tensor_b(0, 0, 1) = 2.0;
    MaterialPropertyCoefficient f(3);
    f.AddCoefficient(attr_mat_a, tensor_a);
    f.AddCoefficient(attr_mat_b, tensor_b);
    const auto &attr_mat = f.GetAttributeToMaterial();
    const auto &mat = f.GetMaterialProperties();
    CHECK(attr_mat[0] != attr_mat[1]);
    CHECK(attr_mat[1] == attr_mat[2]);
    CHECK_THAT(mat(0, 0, attr_mat[0]), WithinRel(1.0, 1.0e-14));
    CHECK_THAT(mat(0, 1, attr_mat[0]), WithinRel(0.5, 1.0e-14));
    CHECK_THAT(mat(0, 0, attr_mat[1]), WithinRel(3.0, 1.0e-14));
    CHECK_THAT(mat(0, 1, attr_mat[1]), WithinAbs(0.0, 0.0));
    CHECK_THAT(mat(1, 1, attr_mat[1]), WithinRel(1.0, 1.0e-14));
  }
}

TEST_CASE("MaterialOperator requires materials for retained mesh domains",
          "[materialoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int size = Mpi::Size(comm);

  // Give every rank two contiguous x-directed coarse cells. Attribute 2 is confined to
  // the first cell on rank 0, with an attribute-1 cell separating it from the partition
  // interface. Thus no other rank sees attribute 2, even through a shared-face neighbor.
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      2 * size, 1, 1, mfem::Element::TETRAHEDRON, 2.0 * size, 1.0, 1.0));
  std::vector<int> partitioning(serial_mesh->GetNE());
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    const auto *element = serial_mesh->GetElement(i);
    const auto *vertices = element->GetVertices();
    double center_x = 0.0;
    for (int j = 0; j < element->GetNVertices(); j++)
    {
      center_x += serial_mesh->GetVertex(vertices[j])[0];
    }
    center_x /= element->GetNVertices();
    serial_mesh->SetAttribute(i, center_x < 1.0 ? 2 : 1);
    partitioning[i] = static_cast<int>(center_x / 2.0);
  }
  serial_mesh->SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning.data());
  Mesh palace_mesh(std::move(par_mesh));

  const auto &local_attributes = palace_mesh.GetCeedAttributes();
  const bool has_attr2 = local_attributes.contains(2);
  CHECK(has_attr2 == (Mpi::Rank(comm) == 0));
  CHECK(palace_mesh.GetNE() > 0);

  config::MaterialData material1;
  material1.attributes = {1};
  config::MaterialData material2;
  material2.attributes = {2};
  config::PeriodicBoundaryData periodic;

  SECTION("All retained attributes have materials")
  {
    CHECK_NOTHROW(MaterialOperator({material1, material2}, periodic,
                                   ProblemType::ELECTROSTATIC, palace_mesh));
  }

  SECTION("Missing material is rejected collectively")
  {
    CHECK_THROWS_WITH(
        MaterialOperator({material1}, periodic, ProblemType::ELECTROSTATIC, palace_mesh),
        Catch::Matchers::ContainsSubstring(
            "Mesh domain attribute 2 has no corresponding entry in "
            "config[\"Domains\"][\"Materials\"]!"));
  }

  SECTION("Pole support uses global material indices")
  {
    material2.permittivity_terms.push_back(
        config::PermittivityPoleData{{-2.0, 5.0}, {1.2, -0.4}});
    MaterialOperator mat_op({material1, material2}, periodic, ProblemType::DRIVEN,
                            palace_mesh);
    CHECK(mat_op.HasFrequencyDependentPermittivity());
    REQUIRE(mat_op.NumFrequencyDependentPermittivityMaterials() == 2);
    CHECK(mat_op.HasFrequencyDependentPermittivityA2(1));
    CHECK(mat_op.HasFrequencyDependentPermittivitySupport(1));
    CHECK(mat_op.GetFrequencyDependentPermittivityAttributes(1) == std::vector<int>{2});
    const std::complex<double> s{0.4, 1.1};
    const std::complex<double> pole{-2.0, 5.0}, residue{1.2, -0.4};
    const auto expected =
        residue * s * s / (s - pole) + std::conj(residue) * s * s / (s - std::conj(pole));
    CHECK(std::abs(mat_op.EvaluateFrequencyDependentPermittivityA2(1, s) - expected) <
          1.0e-14 * std::max(1.0, std::abs(expected)));
  }
}

TEST_CASE("MaterialOperator scalar permittivity pole identities",
          "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  config::PeriodicBoundaryData periodic;

  constexpr double eps_inf = 2.5, plasma_frequency = 1.0, collision_frequency = 0.1;
  const config::MaterialData material(
      {{"Attributes", {1}},
       {"Permittivity",
        {{"HighFrequency", eps_inf},
         {"Terms",
          {{{"Type", "Drude"},
            {"PlasmaFrequency", plasma_frequency},
            {"CollisionFrequency", collision_frequency}}}}}}});
  const double wp = 2.0 * std::numbers::pi * 1.0e9 * plasma_frequency;
  const double gamma = 2.0 * std::numbers::pi * 1.0e9 * collision_frequency;
  const double wp2 = wp * wp;

  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);
  const std::complex<double> s{1.25, 4.0};
  const std::complex<double> lhs = eps_inf * s * s + mat_op.GetConductivity(1)(0, 0) * s +
                                   mat_op.EvaluateFrequencyDependentPermittivityA2(0, s);
  const std::complex<double> rhs = eps_inf * s * s + wp2 * s / (s + gamma);
  CHECK(lhs.real() == Approx(rhs.real()));
  CHECK(lhs.imag() == Approx(rhs.imag()));
  CHECK(mat_op.GetConductivity(1)(0, 0) == Approx(wp2 / gamma));
  CHECK(mat_op.HasConductivity());
  CHECK(mat_op.HasFrequencyDependentPermittivityA2());
  CHECK_THROWS_WITH(mat_op.EvaluateFrequencyDependentPermittivityA2(0, {-gamma, 0.0}),
                    Catch::Matchers::ContainsSubstring("singular at s = pole"));

  // This distant-pole value is O(1). The forbidden long-division K/C/A2 split forms
  // O(1e16) terms whose cancellation cannot reproduce this result in binary64.
  config::MaterialData distant;
  distant.attributes = {1};
  distant.permittivity_terms.push_back(
      config::PermittivityPoleData{{-1.0e8, 0.0}, {1.0e8, 0.0}});
  MaterialOperator distant_op({distant}, periodic, ProblemType::DRIVEN, palace_mesh);
  const std::complex<double> distant_expected = -1.0e8 / std::complex<double>(1.0e8, 1.0);
  CHECK(std::abs(distant_op.EvaluateFrequencyDependentPermittivityA2(0, {0.0, 1.0}) -
                 distant_expected) < 2.0e-16);

  auto anisotropic = material;
  anisotropic.epsilon_r.s = {eps_inf, 2.0 * eps_inf, eps_inf};
  CHECK_THROWS(MaterialOperator({anisotropic}, periodic, ProblemType::DRIVEN, palace_mesh));
  auto lossy = material;
  lossy.tandelta = config::SymmetricMatrixData<3>(0.1);
  CHECK_THROWS(MaterialOperator({lossy}, periodic, ProblemType::DRIVEN, palace_mesh));
  auto conducting = material;
  conducting.sigma = config::SymmetricMatrixData<3>(1.0);
  MaterialOperator conducting_op({conducting}, periodic, ProblemType::DRIVEN, palace_mesh);
  CHECK(conducting_op.GetConductivity(1)(0, 0) == Approx(1.0 + wp2 / gamma));
}

TEST_CASE("MaterialOperator named permittivity model evaluation",
          "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  config::PeriodicBoundaryData periodic;

  constexpr double delta_debye = 0.3, tau_ns = 0.04;
  constexpr double delta_lorentz_under = 0.8, f0_under = 6.0, fg_under = 0.2;
  constexpr double delta_lorentz_over = 0.2, f0_over = 1.0, fg_over = 3.0;
  constexpr std::complex<double> pole{-2.0e12, 5.0e12};
  constexpr std::complex<double> residue{1.2e12, -0.4e12};
  constexpr double strength = 0.0575258, f_lower = 9.07157e-5, f_upper = 159.154956;
  const config::MaterialData material({{"Attributes", {1}},
                                       {"Permittivity",
                                        {{"HighFrequency", 2.08},
                                         {"Terms",
                                          {{{"Type", "Debye"},
                                            {"DeltaPermittivity", delta_debye},
                                            {"RelaxationTime", tau_ns}},
                                           {{"Type", "Lorentz"},
                                            {"DeltaPermittivity", delta_lorentz_under},
                                            {"ResonanceFrequency", f0_under},
                                            {"DampingFrequency", fg_under}},
                                           {{"Type", "Lorentz"},
                                            {"DeltaPermittivity", delta_lorentz_over},
                                            {"ResonanceFrequency", f0_over},
                                            {"DampingFrequency", fg_over}},
                                           {{"Type", "PoleResidue"},
                                            {"Pole", {-2.0e12, 5.0e12}},
                                            {"Residue", {1.2e12, -0.4e12}}},
                                           {{"Type", "DjordjevicSarkar"},
                                            {"Strength", strength},
                                            {"LowerFrequency", f_lower},
                                            {"UpperFrequency", f_upper}}}}}}});
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);
  CHECK(mat_op.HasFrequencyDependentPermittivity());
  CHECK(mat_op.HasFrequencyDependentPermittivityA2());

  constexpr double scale = 2.0 * std::numbers::pi * 1.0e9;
  const std::complex<double> s{0.7e9, 2.3e10};
  const auto Lorentz = [s, scale](double delta, double f0, double fg)
  {
    const double w0 = scale * f0, gamma = scale * fg;
    return delta * w0 * w0 / (s * s + gamma * s + w0 * w0);
  };
  const std::complex<double> expected_epsilon =
      delta_debye / (1.0 + s * tau_ns * 1.0e-9) +
      Lorentz(delta_lorentz_under, f0_under, fg_under) +
      Lorentz(delta_lorentz_over, f0_over, fg_over) + residue / (s - pole) +
      std::conj(residue) / (s - std::conj(pole)) +
      strength * std::log((s + scale * f_upper) / (s + scale * f_lower));
  const auto value = mat_op.EvaluateFrequencyDependentPermittivityA2(0, s);
  CHECK(std::abs(value - s * s * expected_epsilon) <
        2.0e-14 * std::max(1.0, std::abs(value)));
}

TEST_CASE("Djordjevic-Sarkar real-axis and analytic evaluation",
          "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  config::PeriodicBoundaryData periodic;

  constexpr double strength = 0.0575258, f_lower = 9.07157e-5, f_upper = 159.154956;
  const config::MaterialData material({{"Attributes", {1}},
                                       {"Permittivity",
                                        {{"HighFrequency", 1.0},
                                         {"Terms",
                                          {{{"Type", "DjordjevicSarkar"},
                                            {"Strength", strength},
                                            {"LowerFrequency", f_lower},
                                            {"UpperFrequency", f_upper}}}}}}});
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  constexpr double f = 2.5;
  const double omega = 2.0 * std::numbers::pi * 1.0e9 * f;
  const std::complex<double> s{0.0, omega};
  const auto epsilon_term = mat_op.EvaluateFrequencyDependentPermittivityA2(0, s) / (s * s);
  const double expected_real =
      0.5 * strength * std::log((f_upper * f_upper + f * f) / (f_lower * f_lower + f * f));
  const double expected_sigma = 2.0 * std::numbers::pi * electromagnetics::epsilon0_ *
                                strength * f * 1.0e9 *
                                (std::atan(f / f_lower) - std::atan(f / f_upper));
  CHECK(epsilon_term.real() == Approx(expected_real).epsilon(2.0e-14));
  CHECK(-omega * electromagnetics::epsilon0_ * epsilon_term.imag() ==
        Approx(expected_sigma).epsilon(2.0e-14));
  CHECK(0.5 * strength == Approx(0.0287629));
  CHECK(2.0 * std::numbers::pi * electromagnetics::epsilon0_ * strength ==
        Approx(3.20031e-12).epsilon(2.0e-6));

  const std::complex<double> s_off_axis{4.0e8, 1.7e10};
  const std::complex<double> expected_off_axis =
      strength * std::log((s_off_axis + 2.0 * std::numbers::pi * 1.0e9 * f_upper) /
                          (s_off_axis + 2.0 * std::numbers::pi * 1.0e9 * f_lower));
  const auto evaluated_off_axis =
      mat_op.EvaluateFrequencyDependentPermittivityA2(0, s_off_axis) /
      (s_off_axis * s_off_axis);
  CHECK(std::abs(evaluated_off_axis - expected_off_axis) <
        1.0e-14 * std::max(1.0, std::abs(expected_off_axis)));
}

TEST_CASE("MaterialOperator anisotropic wave-number bound", "[materialoperator][Serial]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {4.0, 1.0, 1.0};
  material.epsilon_r.s = {1.0, 4.0, 1.0};

  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  // Same-axis products max(mu_i * epsilon_i) are only 4, but a wave polarized with E
  // along y and H along x samples mu_x * epsilon_y = 16.
  CHECK(mat_op.GetMaxMuEpsilon() == Approx(16.0));
}

TEST_CASE("MaterialOperator wave-number bound is global", "[materialoperator][Parallel]")
{
  // The boundary-mode shift target, and with it the mode ranking, derives from this bound,
  // so every rank must see the same value. A rank-local bound made wave-port results
  // depend on the partition (awslabs/palace#565).
  MPI_Comm comm = Mpi::World();
  REQUIRE(Mpi::Size(comm) >= 2);

  // A row of four elements with the dielectric at one end. Rank 1 owns the two elements at
  // the other end, whose shared-face neighbor is vacuum, so the dielectric is not among its
  // local or ghost attributes. Any further ranks own no elements.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 1, 1, mfem::Element::HEXAHEDRON, 4.0, 1.0, 1.0));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    serial_mesh->SetAttribute(i, (i == 0) ? 2 : 1);
  }
  serial_mesh->SetAttributes();
  int partitioning[4] = {0, 0, 1, 1};
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning);
  Mesh palace_mesh(std::move(par_mesh));

  config::MaterialData vacuum, dielectric;
  vacuum.attributes = {1};
  dielectric.attributes = {2};
  dielectric.epsilon_r.s = {9.0, 9.0, 9.0};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({vacuum, dielectric}, periodic, ProblemType::DRIVEN, palace_mesh);
  // Reduce so that rank 0, the only rank that prints, also reports a mismatch elsewhere.
  double bound_min = mat_op.GetMaxMuEpsilon(), bound_max = bound_min;
  Mpi::GlobalMin(1, &bound_min, comm);
  Mpi::GlobalMax(1, &bound_max, comm);
  CHECK(bound_min == Approx(9.0));
  CHECK(bound_max == Approx(9.0));
}

TEST_CASE("MaterialOperator static limit of dispersive materials",
          "[materialoperator][Serial]")
{
  // Electrostatics uses ε(0) = ε∞ + χ(0): Δε for a Debye term, Δε for a Lorentz pair, and
  // Strength log(ω_hi/ω_lo) for Djordjevic-Sarkar. A pole at zero (Drude) has no finite
  // static permittivity and is rejected. Magnetostatics ignores the permittivity.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  config::PeriodicBoundaryData periodic;

  const double eps_inf = 2.0, tau = 0.5, w0 = 3.0, gamma = 0.2;
  config::MaterialData material;
  material.attributes = {1};
  material.epsilon_r.s = {eps_inf, eps_inf, eps_inf};
  material.permittivity_terms.push_back(
      config::PermittivityPoleData{-1.0 / tau, 1.5 / tau});  // Debye, Δε = 1.5
  {
    // Lorentz with Δε = 0.7, as r/(s - p) + conj with p the upper root of s² + γs + ω₀²
    const std::complex<double> p(-0.5 * gamma, std::sqrt(w0 * w0 - 0.25 * gamma * gamma));
    material.permittivity_terms.push_back(
        config::PermittivityPoleData{p, 0.7 * w0 * w0 / (p - std::conj(p))});
  }
  material.permittivity_terms.push_back(config::DjordjevicSarkarData{0.1, 1.0e-3, 1.0e3});
  const double eps0 = eps_inf + 1.5 + 0.7 + 0.1 * std::log(1.0e6);
  CHECK(material.StaticPermittivity() == Approx(eps0));

  MaterialOperator electrostatic({material}, periodic, ProblemType::ELECTROSTATIC,
                                 palace_mesh);
  CHECK_FALSE(electrostatic.HasFrequencyDependentPermittivity());
  CHECK(electrostatic.GetPermittivityReal()(0, 0, 0) == Approx(eps0));
  MaterialOperator magnetostatic({material}, periodic, ProblemType::MAGNETOSTATIC,
                                 palace_mesh);
  CHECK_FALSE(magnetostatic.HasFrequencyDependentPermittivity());

  config::MaterialData drude;
  drude.attributes = {1};
  drude.permittivity_terms = {config::PermittivityPoleData{0.0, 1.0},
                              config::PermittivityPoleData{-0.1, -1.0}};
  CHECK_THROWS_WITH(
      MaterialOperator({drude}, periodic, ProblemType::ELECTROSTATIC, palace_mesh),
      Catch::Matchers::ContainsSubstring("finite static material permittivity"));
}

TEST_CASE("MaterialOperator dispersive wave-number bound", "[materialoperator][Serial]")
{
  // A Debye term χ(s) = Δε/(1 + sτ) raises Re ε above ε∞ to ε∞ + Δε/(1 + ω²τ²), which
  // bounds the wave number of a port mode in the material at ω.
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  const double eps_inf = 2.0, delta = 3.0, tau = 0.5;
  config::MaterialData material;
  material.attributes = {1};
  material.mu_r.s = {2.0, 2.0, 2.0};
  material.epsilon_r.s = {eps_inf, eps_inf, eps_inf};
  material.permittivity_terms.push_back(
      config::PermittivityPoleData{-1.0 / tau, delta / tau});
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  CHECK(mat_op.GetMaxMuEpsilon() == Approx(2.0 * eps_inf));
  for (double omega : {0.1, 1.0, 10.0})
  {
    CAPTURE(omega);
    CHECK(mat_op.GetFrequencyDependentMaxMuEpsilon(omega) ==
          Approx(2.0 * (eps_inf + delta / (1.0 + omega * omega * tau * tau))));
  }
}

TEST_CASE("MaterialOperator utility functions", "[materialoperator][Serial]")
{
  SECTION("IsOrthonormal")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Orthonormal vectors")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 2.0, 3.0};
      REQUIRE(internal::mat::IsOrthonormal(data));
    }

    SECTION("Non-normalized vectors")
    {
      data.v[0] = {2.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsOrthonormal(data));
    }

    SECTION("Non-orthogonal vectors")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.5, 0.866, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsOrthonormal(data));
    }
  }

  SECTION("IsValid")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Valid orthonormal matrix with positive eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsValid(data));
    }

    SECTION("Valid orthonormal matrix with different positive eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 3.0, 4.0};
      REQUIRE(internal::mat::IsValid(data));
    }

    SECTION("Invalid - zero eigenvalue")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 0.0, 1.0};
      REQUIRE(!internal::mat::IsValid(data));
    }
  }

  SECTION("IsIdentity")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Identity matrix")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsIdentity(data));
    }

    SECTION("Non-identity - different eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 1.0, 1.0};
      REQUIRE(!internal::mat::IsIdentity(data));
    }

    SECTION("Identity but rotated basis")
    {
      data.v[0] = {0.0, 1.0, 0.0};
      data.v[1] = {1.0, 0.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 1.0, 1.0};
      REQUIRE(internal::mat::IsIdentity(data));
    }
  }

  SECTION("IsIsotropic")
  {
    config::SymmetricMatrixData<3> data(1.0);

    SECTION("Isotropic material - all eigenvalues equal")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.5, 2.5, 2.5};
      REQUIRE(internal::mat::IsIsotropic(data));
    }

    SECTION("Anisotropic material - different eigenvalues")
    {
      data.v[0] = {1.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {1.0, 2.0, 1.0};
      REQUIRE(!internal::mat::IsIsotropic(data));
    }

    SECTION("Invalid - non-orthonormal but equal eigenvalues")
    {
      data.v[0] = {2.0, 0.0, 0.0};
      data.v[1] = {0.0, 1.0, 0.0};
      data.v[2] = {0.0, 0.0, 1.0};
      data.s = {2.0, 2.0, 2.0};
      REQUIRE(!internal::mat::IsIsotropic(data));
    }
  }
}

TEST_CASE("MaterialOperator ignores ghost-only submesh attributes",
          "[materialoperator][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm);
  REQUIRE(Mpi::Size(comm) >= 2);

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::HEXAHEDRON, 2.0, 1.0, 1.0));
  for (int i = 0; i < serial_mesh->GetNE(); i++)
  {
    serial_mesh->SetAttribute(i, 7);
  }
  serial_mesh->SetAttributes();
  serial_mesh->EnsureNCMesh(true);

  int partitioning[2] = {0, 1};
  auto parent = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning);
  mfem::Array<mfem::Refinement> refinements;
  if (rank == 0)
  {
    refinements.Append(mfem::Refinement(0));
  }
  parent->GeneralRefinement(refinements, -1);

  mfem::Array<int> surface_attributes;
  surface_attributes.Append(1);
  auto submesh = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(*parent, surface_attributes));
  Mesh palace_mesh(std::move(submesh));
  mesh::RemapSubMeshAttributes(static_cast<mfem::ParSubMesh &>(palace_mesh.Get()));
  palace_mesh.RebuildCeedAttributes();

  const bool owns_elements = rank < 2;
  CHECK((palace_mesh.GetNE() > 0) == owns_elements);
  for (int i = 0; i < palace_mesh.GetNE(); i++)
  {
    CHECK(palace_mesh.Get().GetAttribute(i) == 7);
  }
  const auto &local_attributes = palace_mesh.GetCeedAttributes();
  if (owns_elements)
  {
    CHECK(local_attributes.contains(1));
    CHECK(local_attributes.contains(7));
    // The remap must refresh the cached distinct-attribute array.
    CHECK(palace_mesh.Get().attributes.Max() == 7);
  }
  else
  {
    CHECK(local_attributes.empty());
  }

  config::MaterialData material;
  material.attributes = {7};
  config::PeriodicBoundaryData periodic;
  CHECK_NOTHROW(MaterialOperator({material}, periodic, ProblemType::DRIVEN, palace_mesh));
}

}  // namespace palace
