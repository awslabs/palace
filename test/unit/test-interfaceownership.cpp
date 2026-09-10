// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <filesystem>
#include <fstream>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/interfaceownership.hpp"
#include "utils/iodata.hpp"
#include "utils/ownershipquadrature.hpp"
#include "utils/units.hpp"

namespace palace
{
namespace
{
class OwnershipFile
{
public:
  std::string path;
  explicit OwnershipFile(const std::string &rows)
  {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    path = (std::filesystem::temp_directory_path() /
            ("palace-ownership-" + std::to_string(Mpi::Rank(Mpi::World())) + "-" +
             std::to_string(stamp) + ".csv"))
               .string();
    std::ofstream output(path);
    output << "Group,Conductor,Slot,Role,Radius,X0,Y0,Z0,X1,Y1,Z1\n" << rows;
  }
  ~OwnershipFile() { std::filesystem::remove(path); }
};
const std::string two_slots = "1,1,0,Signature,1,0,-2,0,0,2,0\n"
                              "1,1,1,Signature,1,1,-2,0,1,2,0\n";
}  // namespace

TEST_CASE("Interface ownership coordinates and hierarchy",
          "[interfaceownership][Serial][Parallel]")
{
  OwnershipFile file(two_slots);
  InterfaceOwnershipPartition partition(file.path, 1, 1.0), scaled(file.path, 1, 4.0);
  mfem::Vector point(3);
  for (double x : {0.1, 0.49, 0.5, 0.51, 0.9})
  {
    point = 0.0;
    point(0) = x;
    point(2) = 0.2;
    const int expected = x <= 0.5 ? 0 : 1;
    CHECK(partition.SelectSlot(point) == expected);
    point /= 4.0;
    CHECK(scaled.SelectSlot(point) == expected);
  }
  OwnershipFile reversed("1,1,1,Signature,1,1,-2,0,1,2,0\n"
                         "1,1,0,Signature,1,0,-2,0,0,2,0\n");
  point = 0.0;
  point(0) = 0.5;
  CHECK(InterfaceOwnershipPartition(reversed.path, 1, 1.0).SelectSlot(point) == 0);

  OwnershipFile hierarchical("0,1,0,Boundary,1,0,-0.1,0,0,0.1,0\n"
                             "0,2,0,Boundary,1,2,-0.1,0,2,0.1,0\n"
                             "0,1,0,Signature,1,0,-5,0,0,5,0\n"
                             "0,2,1,Signature,1,0.25,-5,0,0.25,5,0\n");
  InterfaceOwnershipPartition sa(hierarchical.path, 0, 1.0);
  point = 0.0;
  point(0) = 0.2;
  point(2) = 10.0;
  CHECK(sa.SelectSlot(point) ==
        0);  // Physical conductor first; z is irrelevant to SA gate.
  point(1) = 4.0;
  CHECK(sa.SelectSlot(point) == 1);  // Beyond 3R, all signatures compete.
  CHECK_THROWS(InterfaceOwnershipPartition(file.path, 99, 1.0));
  CHECK_THROWS(InterfaceOwnershipPartition(file.path, 1, 0.0));
  OwnershipFile malformed(two_slots + "2,1,0,Signature,1,nan,0,0,1,0,0\n");
  CHECK_THROWS(InterfaceOwnershipPartition(malformed.path, 1, 1.0));
  OwnershipFile zero("1,1,0,Signature,1,0,0,0,0,0,0\n");
  CHECK_THROWS(InterfaceOwnershipPartition(zero.path, 1, 1.0));
  OwnershipFile missing("1,1,0,Boundary,1,0,-1,0,0,1,0\n"
                        "1,2,0,Signature,1,1,-1,0,1,1,0\n");
  CHECK_THROWS(InterfaceOwnershipPartition(missing.path, 1, 1.0));

  config::InterfaceDielectricData data;
  Units units(1e-6, 4e-6);
  config::Nondimensionalize(units, data);
  CHECK(data.ownership_coordinate_scale == 4.0);
}

TEST_CASE("Ownership quadrature positive moments", "[interfaceownership][Serial][Parallel]")
{
  for (int order : {2, 12, 20, 32, 48})
  {
    const OwnershipQuadrature quadrature(order);
    const auto &rule = quadrature.Get(mfem::Geometry::TRIANGLE);
    bool positive = true;
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      positive &= rule.IntPoint(q).weight > 0.0;
    }
    REQUIRE(positive);
    for (int a = 0; a <= std::min(6, order); a++)
    {
      for (int b = 0; b <= std::min(6, order) - a; b++)
      {
        double integral = 0.0;
        for (int q = 0; q < rule.GetNPoints(); q++)
        {
          const auto &p = rule.IntPoint(q);
          integral += p.weight * std::pow(p.x, a) * std::pow(p.y, b);
        }
        const double expected =
            std::tgamma(a + 1.) * std::tgamma(b + 1.) / std::tgamma(a + b + 3.);
        CHECK_THAT(integral, Catch::Matchers::WithinAbs(expected, 1e-13));
      }
    }
  }
}

TEST_CASE("Interface ownership exhaustive surface quadrature",
          "[interfaceownership][Serial][Parallel]")
{
  OwnershipFile file(two_slots);
  auto serial = mfem::Mesh::MakeCartesian3D(2, 2, 1, mfem::Element::TETRAHEDRON);
  int top = -1;
  mfem::Vector center(3);
  for (int be = 0; be < serial.GetNBE(); be++)
  {
    auto *T = serial.GetBdrElementTransformation(be);
    T->Transform(mfem::Geometries.GetCenter(T->GetGeometryType()), center);
    if (std::abs(center(2) - 1.0) < 1e-12)
    {
      top = serial.GetBdrAttribute(be);
      break;
    }
  }
  REQUIRE(top > 0);
  Mesh mesh(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
  mfem::H1_FECollection h1(2, 3);
  mfem::ND_FECollection nd(2, 3);
  FiniteElementSpace h1_space(mesh, &h1), nd_space(mesh, &nd);
  fem::DefaultIntegrationOrder::p_trial = 2;
  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator materials({material}, periodic, ProblemType::ELECTROSTATIC, mesh);
  GridFunction field(nd_space, false);
  mfem::Vector value(3);
  value = 0.0;
  value(2) = 1.0;
  mfem::VectorConstantCoefficient coefficient(value);
  field.Real().ProjectCoefficient(coefficient);

  config::InterfaceDielectricData data;
  data.attributes = {top};
  data.type = InterfaceDielectric::MA;
  data.t = 0.2;
  data.epsilon_r = 2.0;
  config::BoundaryPostData plain;
  plain.dielectric.emplace(1, data);
  SurfacePostOperator unmasked(plain, ProblemType::ELECTROSTATIC, materials, h1_space,
                               nd_space);
  const double reference = unmasked.GetInterfaceElectricFieldEnergy(1, field);
  CHECK_THAT(reference, Catch::Matchers::WithinAbs(0.05, 1e-12));

  data.ownership_data_file = file.path;
  data.ownership_group = 1;
  data.ownership_quadrature_order = 12;
  data.ownership_slot = 0;
  config::BoundaryPostData partitioned;
  partitioned.dielectric.emplace(1, data);
  data.ownership_slot = 1;
  partitioned.dielectric.emplace(2, data);
  SurfacePostOperator post(partitioned, ProblemType::ELECTROSTATIC, materials, h1_space,
                           nd_space);
  const double a = post.GetInterfaceElectricFieldEnergy(1, field);
  const double b = post.GetInterfaceElectricFieldEnergy(2, field);
  CHECK_THAT(a + b, Catch::Matchers::WithinAbs(reference, 1e-12));
  CHECK_THAT(a, Catch::Matchers::WithinAbs(reference / 2, 1e-12));
  CHECK_THAT(b, Catch::Matchers::WithinAbs(reference / 2, 1e-12));

  auto matrix_data = partitioned;
  for (auto &[index, entry] : matrix_data.dielectric)
  {
    entry.localize_edge_energy = true;
    entry.save_local_edge_energy = false;
    entry.edge_attributes = {top};
    entry.edge_distances = {2.0};
    entry.edge_frame_normal = std::array<double, 3>{0.0, 0.0, 1.0};
  }
  SurfacePostOperator matrix_post(matrix_data, ProblemType::ELECTROSTATIC, materials,
                                  h1_space, nd_space);
  GridFunction second(nd_space, false);
  value(2) = 2.0;
  mfem::VectorConstantCoefficient second_coefficient(value);
  second.Real().ProjectCoefficient(second_coefficient);
  const auto matrices =
      matrix_post.GetInterfaceElectricFieldEnergyMatrices({&field, &second}, {});
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      const double sum =
          matrices.at(1).at(0).energy_total(i, j) + matrices.at(2).at(0).energy_total(i, j);
      CHECK_THAT(sum, Catch::Matchers::WithinAbs(reference * (i + 1) * (j + 1), 1e-12));
    }
  }
  partitioned.dielectric.erase(2);
  CHECK_THROWS(SurfacePostOperator(partitioned, ProblemType::ELECTROSTATIC, materials,
                                   h1_space, nd_space));
  data.ownership_slot = 0;
  partitioned.dielectric.emplace(2, data);
  CHECK_THROWS(SurfacePostOperator(partitioned, ProblemType::ELECTROSTATIC, materials,
                                   h1_space, nd_space));
  partitioned.dielectric.at(2).ownership_slot = 1;
  partitioned.dielectric.at(2).ownership_quadrature_order = 14;
  CHECK_THROWS(SurfacePostOperator(partitioned, ProblemType::ELECTROSTATIC, materials,
                                   h1_space, nd_space));
}
}  // namespace palace
