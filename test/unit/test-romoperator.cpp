// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "drivers/drivensolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"

using namespace palace;
using namespace nlohmann;
using namespace Catch::Matchers;

class RomOperatorTest : public RomOperator
{
public:
  using RomOperator::RomOperator;
  auto &GetWeightOp() const { return weight_op_W; }
};

auto LoadScaleParMesh2(IoData &iodata, MPI_Comm world_comm)
{
  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::ReadMesh(iodata, world_comm));
    iodata.NondimensionalizeInputs(*mfem_mesh[0]);
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  return mesh_;
}

TEST_CASE("MinimalRationalInterpolation", "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();

  auto fn_tan_shift = [](double z)
  { return std::tan(0.5 * M_PI * (z - std::complex<double>(1., 1.))); };

  // Test scalar case: 2 sample points for 4 x 2 vector
  MinimalRationalInterpolation mri_1(6);

  CHECK(mri_1.GetSamplePoints() == std::vector<double>{});
  CHECK_THROWS(mri_1.FindMaxError(1));

  for (double x_sample : {-2.0, -1.0, 1.0, 2.0})
  {
    auto tan_eval = fn_tan_shift(x_sample) / double(Mpi::Size(comm));
    std::vector<std::complex<double>> tmp = {x_sample * tan_eval, -tan_eval, 5. * tan_eval,
                                             10. * x_sample * tan_eval,
                                             2 * x_sample * x_sample * x_sample * tan_eval};
    ComplexVector c_vec(tmp.size());
    c_vec.Set(tmp.data(), tmp.size(), false);

    mri_1.AddSolutionSample(x_sample, c_vec, comm, Orthogonalization::MGS);
  }

  CHECK(mri_1.GetSamplePoints().size() == 4);
  CHECK(mri_1.GetSamplePoints() == std::vector<double>{-2.0, -1.0, 1.0, 2.0});

  // By symmetry of poles max error should be at zero.
  auto max_err_1 = mri_1.FindMaxError(5);
  REQUIRE(max_err_1.size() == 5);

  // By symmetry highest error should be at zero.
  CHECK_THAT(max_err_1[0], Catch::Matchers::WithinAbsMatcher(0.0, 1e-6));

  // Test that elements of max_error are unique.
  // TODO: get better test for multiple N.
  std::sort(max_err_1.begin(), max_err_1.end());
  CHECK(std::adjacent_find(max_err_1.begin(), max_err_1.end()) == max_err_1.end());

  // TODO: Add more stringent tests of MRI, including estimating poles.
}

// TODO: Do some more basic RomOperator Ctor Checks

TEST_CASE("RomOperator-Synthesis-Port-Cube111", "[romoperator][Serial][Parallel]")
{
  // Work with a simple 1x1x1 Cube

  using VT = palace::Units::ValueType;
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 x 3 different test configuration.
  auto solver_order = GENERATE(1);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DIR) /
                                         "lumpedport_mesh/cube_mesh_1_1_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DIR) /
                                          "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"));

  double L0 = 1.0e-6;
  double Lc = 7.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", PALACE_TEST_DIR}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0},
                                               {"LossTan", 0.0}})})}};

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      {"PEC", json::object({{"Attributes", json::array({})}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = solver_order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "SuperLU"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(std::move(setup_json));

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  const auto &mesh = space_op.GetNDSpace().GetParMesh();

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  auto &weight_op = prom_op.GetWeightOp();
  REQUIRE(weight_op.has_value());

  const auto &W_bulk = weight_op->W_inner_product_weight_bulk;
  const auto &W_port = weight_op->W_inner_product_weight_port;

  CHECK(W_bulk->NumRows() >= 12);
  CHECK(W_bulk->NumCols() >= 12);

  Eigen::MatrixXd print_mat = Eigen::MatrixXd::Zero(12, 12);

  std::cout << "Print port:\n";
  for (int i = 0; i < 12; i++)
  {
    mfem::Vector v;
    v.SetSize(12);
    v.UseDevice(true);

    mfem::Vector w;
    w.UseDevice(true);
    w.SetSize(12);

    w = 0.0;
    w[i] = 1.0;
    W_port->Mult(w, v);

    for (int j = 0; j < 12; j++)
    {
      print_mat(i, j) = std::abs(v(j)) > 1e-8 ? v(j) : 0;
    }
  }

  Eigen::IOFormat HeavyFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
  std::cout << print_mat.format(HeavyFmt) << std::endl;

  std::cout << "Print bulk:\n";
  for (int i = 0; i < 12; i++)
  {
    mfem::Vector v;
    v.SetSize(12);
    v.UseDevice(true);

    mfem::Vector w;
    w.UseDevice(true);
    w.SetSize(12);

    w = 0.0;
    w[i] = 1.0;
    W_bulk->Mult(w, v);

    for (int j = 0; j < 12; j++)
    {
      print_mat(i, j) = std::abs(v(j)) > 1e-8 ? v(j) : 0;
    }
  }

  std::cout << print_mat.format(HeavyFmt) << std::endl;

  std::cout << "Print combo:\n";
  for (int i = 0; i < 12; i++)
  {
    mfem::Vector v;
    v.SetSize(12);
    v.UseDevice(true);

    mfem::Vector w;
    w.UseDevice(true);
    w.SetSize(12);

    w = 0.0;
    w[i] = 1.0;
    weight_op->Mult(w, v);

    for (int j = 0; j < 12; j++)
    {
      print_mat(i, j) = std::abs(v(j)) > 1e-8 ? v(j) : 0;
    }
  }

  std::cout << print_mat.format(HeavyFmt) << std::endl;
}
