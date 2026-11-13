// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <fstream>
#include <numbers>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using namespace nlohmann;
using namespace Catch::Matchers;
using namespace std::complex_literals;

// Assemble operators as Eigen matrices for simpler testing.
template <typename T>
auto toEigenMatrix(const T &op, int n)
{
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
  mfem::Vector w(n), v(n);
  w.UseDevice(true);
  v.UseDevice(true);
  for (int i = 0; i < n; i++)
  {
    w = 0.0;
    w[i] = 1.0;
    op.Mult(w, v);
    for (int j = 0; j < n; j++)
    {
      mat(j, i) = v(j);
    }
  }
  return mat;
}

class DrivenPostOperatorTest : public PostOperator<ProblemType::DRIVEN>
{
public:
  using PostOperator<ProblemType::DRIVEN>::PostOperator;
  const Measurement &GetMeasurement() const { return measurement_cache; }
};

class RomOperatorTest : public RomOperator
{
public:
  using RomOperator::AddAuxBlockDirections;
  using RomOperator::ApplyComplexPolynomialFitCorrections;
  using RomOperator::AugmentedPencil;
  using RomOperator::BuildAugmentedPencil;
  using RomOperator::CalculateNormalizedPROMMatrices;
  using RomOperator::NormalizedMatrices;
  using RomOperator::RomOperator;
  using RomOperator::WavePortAuxBlock;
  using RomOperator::WavePortDispersionFit;
  auto &GetWeightOp() const { return weight_op_W; }
  auto &GetOrthR() const { return orth_R; }
  auto &GetVectors() const { return V; }
  auto &GetKr() const { return Kr; }
  auto &GetCr() const { return Cr; }
  auto &GetMr() const { return Mr; }
  auto &GetAr() const { return Ar; }
};

auto LoadScaleParMesh2(IoData &iodata, MPI_Comm world_comm)
{
  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    auto smesh = mesh::Load(iodata, world_comm);
    int ne = smesh ? smesh->GetNE() : 0;
    Mpi::GlobalMax(1, &ne, world_comm);
    if (Mpi::Size(world_comm) > ne)
    {
      SKIP("Mesh has fewer elements than processes");
    }
    if (iodata.model.Lc <= 0.0)
    {
      iodata.model.Lc = mesh::ComputeReferenceLength(smesh, world_comm);
    }
    iodata.NondimensionalizeInputs(smesh);
    mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), world_comm));
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
  auto fn_tan_shift = [](double z)
  { return std::tan(0.5 * std::numbers::pi * (z - std::complex<double>(1., 1.))); };

  // Test scalar case: 2 sample points for 4 x 2 vector
  MinimalRationalInterpolation mri_1;

  CHECK(mri_1.GetSamplePoints() == std::vector<double>{});
  CHECK_THROWS(mri_1.FindMaxError(1));

  for (double x_sample : {-2.0, -1.0, 1.0, 2.0})
  {
    auto tan_eval = fn_tan_shift(x_sample);

    // The rational vector function we want to evaluate. The vector has common poles all
    // determined by tan. The minimal-rational interpolation method finds coefficients using
    // SVD on the R matrix from a QR decomposition of samples. If the rank of the vector is
    // small or degenerate, this can have trivial solution and fail.
    //
    // Here we pick a 5 component vector where 2 pairs are degenerate (scalar multiples of
    // each other) which is sufficiently large but challenging a problem.
    Eigen::VectorXcd c_vec(5);
    c_vec << -tan_eval, 5. * tan_eval, x_sample * tan_eval, 10. * x_sample * tan_eval,
        2 * x_sample * x_sample * x_sample * tan_eval;

    mri_1.AddSolutionSample(x_sample, c_vec);
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
  std::ranges::sort(max_err_1);
  CHECK(std::adjacent_find(max_err_1.begin(), max_err_1.end()) == max_err_1.end());

  // TODO: Add more stringent tests of MRI, including estimating poles.
}

TEST_CASE("MinimalRationalInterpolation-DistributedSearch",
          "[romoperator][Serial][Parallel]")
{
  // The search split over the processes gives the same point as one scan. The basis grows
  // between the samples.
  MPI_Comm comm = Mpi::World();
  MinimalRationalInterpolation mri;
  int n = 2;
  for (double x : {-2.0, -1.3, 0.4, 1.0, 2.0})
  {
    const auto t = std::tan(0.5 * std::numbers::pi * (x - std::complex<double>(1., 1.)));
    Eigen::VectorXcd c = Eigen::VectorXcd::Zero(++n);
    c.head(3) << -t, 5. * t, 2. * x * x * x * t;
    c(n - 1) += 0.1 * x * t;
    mri.AddSolutionSample(x, c);
  }
  CHECK(mri.FindMaxError(1, comm) == mri.FindMaxError(1));
}

TEST_CASE("MinimalRationalInterpolation-NarrowSpan", "[romoperator][Serial]")
{
  // Samples 2e-11 apart relative to z: delta = span / 1e6 is below the spacing of doubles
  // near z, so the scan must not step by accumulating delta.
  MinimalRationalInterpolation mri;
  const double z0 = 100.0, z1 = 100.0 + 2.0e-9;
  for (int i = 0; i < 2; i++)
  {
    Eigen::VectorXcd c = Eigen::VectorXcd::Zero(2);
    c(i) = 1.0;
    mri.AddSolutionSample(i == 0 ? z0 : z1, c);
  }
  const auto z_max = mri.FindMaxError(1);
  REQUIRE(z_max.size() == 1);
  CHECK(z_max[0] >= z0);
  CHECK(z_max[0] <= z1);
}

// TODO: Do some more basic RomOperator Ctor Checks

// Basic checks of ROM construction in the of synthesis. Checks hybrid domain-boundary
// inner-product weight and port overlap. Works with a simple 1x1x1 Cube. This is a serial
// test as hex mesh only has a single element.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "RomOperator-Synthesis-Port-Cube111",
                 "[romoperator][Serial]")
{
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 different test configuration.
  size_t order = GENERATE(1UL, 2UL, 3UL);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_1_1_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"));

  double L0 = 1.0e-6;

  // Pick a number != 1. This should not matter as it should be be correctly cancelled out
  // in unit conversions, but we want to test this here.
  double Lc = 7.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
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

  double port_ref_R = 50.0;

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      // {"PEC", json::object({{"Attributes", json::array({})}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", port_ref_R},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  // Test hybrid weight operator works as expected.
  const auto &weight_op = prom_op.GetWeightOp();
  REQUIRE(weight_op.has_value());

  const auto &W_bulk = weight_op->W_inner_product_weight_bulk;
  const auto &W_port = weight_op->W_inner_product_weight_port;

  std::size_t nr_tdof_expected = 0;
  std::size_t nr_tdof_port_expected = 0;
  if (mesh_is_hex)
  {
    // Reference vales for nr elements for vector Nédélec elements of the first kind.
    auto nr_tdof_ref_hex = std::vector<std::size_t>{0, 12, 54, 144}.at(order);
    auto nr_tdof_ref_squ = std::vector<std::size_t>{0, 4, 12, 24}.at(order);

    nr_tdof_expected = nr_tdof_ref_hex;
    nr_tdof_port_expected = nr_tdof_ref_squ;
  }
  else
  {
    auto nr_tdof_ref_tet = std::vector<std::size_t>{0, 6, 20, 45}.at(order);
    auto nr_tdof_ref_tri = std::vector<std::size_t>{0, 3, 8, 15}.at(order);

    // Counting tdofs. 3D: There are 6 tets in a cube; remove edges. 2D: There are 2 * 6 =
    // 12 open faces, so 12 double-counted; subtract 6 faces without edges. 1D: add total of
    // 19 edges
    nr_tdof_expected = 6 * (nr_tdof_ref_tet - 6 * order) -
                       12 / 2 * (nr_tdof_ref_tri - 3 * order) + 19 * order;
    nr_tdof_port_expected = 2 * nr_tdof_ref_tri - order;  // 2 triangles, 1 shared edge
  }

  CHECK(W_bulk->NumRows() == nr_tdof_expected);
  CHECK(W_bulk->NumCols() == nr_tdof_expected);

  auto W_port_eigen = toEigenMatrix(*W_port, nr_tdof_expected);
  auto W_bulk_eigen = toEigenMatrix(*W_bulk, nr_tdof_expected);
  auto weight_op_eigen = toEigenMatrix(*weight_op, nr_tdof_expected);

  // Check rows/cols where port matrix is non-zero and ensure that corresponding domain
  // matrix is zero.
  std::size_t nr_rows_count = 0;
  for (int i = 0; i < W_port_eigen.rows(); i++)
  {
    if (W_port_eigen.row(i).norm() > 0.0)
    {
      CHECK(W_bulk_eigen.row(i).norm() == 0.0);
      nr_rows_count++;
    }
    else
    {
      // Check there is some overlap in other rows.
      CHECK(W_bulk_eigen.row(i).norm() > 0.0);
    }
  }
  CHECK(nr_rows_count == nr_tdof_port_expected);

  std::size_t nr_cols_count = 0;
  for (int j = 0; j < W_port_eigen.cols(); j++)
  {

    if (W_port_eigen.col(j).norm() > 0.0)
    {
      CHECK(W_bulk_eigen.col(j).norm() == 0.0);
      nr_cols_count++;
    }
    else
    {
      // Check there is some overlap in other cols.
      CHECK(W_bulk_eigen.col(j).norm() > 0.0);
    }
  }
  CHECK(nr_cols_count == nr_tdof_port_expected);

  // Check element-wise sum is ok.
  for (int i = 0; i < weight_op_eigen.rows(); i++)
  {
    for (int j = 0; j < weight_op_eigen.cols(); j++)
    {
      CHECK_THAT(weight_op_eigen(i, j),
                 WithinRel(W_bulk_eigen(i, j) + W_port_eigen(i, j)) ||
                     WithinAbs(0.0, 1e-18));
    }
  }

  // Now test against port vector. The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e
  // L_e / W_e but with Z_R = 1.0. In this case this = 1.0. See normalization tests of
  // LumpedPort_BasicTests_1ElementPort_Cube321 in test-lumpedportintegration.cpp.
  ComplexVector port_primary_et;
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et);
  Vector port_primary_ht_cn_tmp = port_primary_et.Real();

  W_port->Mult(port_primary_et.Real(), port_primary_ht_cn_tmp);
  auto overlap_port = port_primary_et.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_port, WithinRel(1.0));

  W_bulk->Mult(port_primary_et.Real(), port_primary_ht_cn_tmp);
  auto overlap_bulk = port_primary_et.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_bulk, WithinAbs(0.0, 1e-15));

  auto v_workspace = port_primary_et.Real();

  if (Mpi::Size(world_comm) == 1)
  {
    // Rank local overlap.
    auto overlap_combined_local =
        (*weight_op)
            .InnerProduct(port_primary_et.Real(), port_primary_et.Real(), v_workspace);
    CHECK_THAT(overlap_combined_local, WithinRel(1.0));
  }

  // Global overlap.
  auto overlap_combined = (*weight_op)
                              .InnerProduct(world_comm, port_primary_et.Real(),
                                            port_primary_et.Real(), v_workspace);
  CHECK_THAT(overlap_combined, WithinRel(1.0));

  // Test actually adding port primary vectors to PROM.
  prom_op.AddLumpedPortModesForSynthesis();
  CHECK(prom_op.GetReducedDimension() == 1);
  const auto _norm = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = _norm.L_inv;
  const auto &resistance_R_inv = _norm.R_inv;
  const auto &capacitance_C = _norm.C;
  const auto orth_R = prom_op.GetOrthR();

  CHECK(((inductance_L_inv->rows() == 1) && (inductance_L_inv->cols() == 1)));
  CHECK(((resistance_R_inv->rows() == 1) && (resistance_R_inv->cols() == 1)));
  CHECK(((capacitance_C->rows() == 1) && (capacitance_C->cols() == 1)));
  CHECK(((orth_R.rows() == 1) && (orth_R.cols() == 1)));

  // This should be the same as norm of primary port vector ht_cn with Z_R = 1.0.
  // For single square port this is 1.0
  CHECK_THAT(orth_R(0, 0), WithinRel(1.0));
  CHECK_THAT((*resistance_R_inv)(0, 0).real(), WithinRel(1.0 / port_ref_R));
  // L, C have some component on port 1 from the bulk. So only check they exist.
  CHECK((*inductance_L_inv)(0, 0).real() != 0.0);
  CHECK((*capacitance_C)(0, 0).real() != 0.0);

  CHECK_THAT((*inductance_L_inv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*resistance_R_inv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*capacitance_C)(0, 0).imag(), WithinAbs(0.0, 1e-14));
}

// Basic ROM port construction on a more complex mesh including a composite port and and L,C
// port. Works on a 3x2x1 cubic mesh.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-Synthesis-Port-Cube321",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 different test configuration.
  size_t order = GENERATE(1UL, 2UL, 3UL);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));

  const double L0 = 1.0e-6;
  const double Lc = 7.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  // Reasonable reference values in correct frequency range.
  const double port1_ref_R = 75.0;
  const double port1_ref_C = 5.5e-5;
  const double port1_ref_L = 1.486e-20;

  const double port2_ref_R = 50.0;

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({
           json::object(
               {{"Index", 1},
                {"R", port1_ref_R},
                {"C", port1_ref_C},
                {"L", port1_ref_L},
                {"Excitation", uint(0)},
                {"Elements",
                 json::array(
                     {json::object({{"Attributes", json::array({1})},  // dx,dy=1,1
                                    {"Direction", "+Y"}}),
                      json::object({{"Attributes", json::array({9, 11})},  // dx,dy=1,2
                                    {"Direction", "+Y"}}),
                      json::object({{"Attributes", json::array({2, 6, 10})},  // dx,dy = 3,1
                                    {"Direction", "+X"}})})}}),
           json::object({{"Index", 2},
                         {"R", port2_ref_R},
                         {"Excitation", uint(1)},
                         {"Elements", json::array({json::object(
                                          {{"Attributes", json::array({14})},  // dx,dy=1,1
                                           {"Direction", "+Z"}})})}})          //,
       })}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  // Test hybrid weight operator works as expected.
  const auto &weight_op = prom_op.GetWeightOp();
  REQUIRE(weight_op.has_value());

  const auto &W_bulk = weight_op->W_inner_product_weight_bulk;
  const auto &W_port = weight_op->W_inner_product_weight_port;

  // If serial — concretize as an Eigen matrix on a thread and manually check zeros of row
  // and columns.
  if (Mpi::Size(world_comm) == 1)
  {

    auto W_port_eigen = toEigenMatrix(*W_port, W_port->NumCols());
    auto W_bulk_eigen = toEigenMatrix(*W_bulk, W_port->NumCols());

    // Check rows/cols where port matrix is non-zero and ensure that corresponding domain
    // matrix is zero.
    for (int i = 0; i < W_port_eigen.rows(); i++)
    {
      if (W_port_eigen.row(i).norm() > 0.0)
      {
        CHECK(W_bulk_eigen.row(i).norm() == 0.0);
      }
      else
      {
        // Check there is some overlap in other rows.
        CHECK(W_bulk_eigen.row(i).norm() > 0.0);
      }
    }

    for (int j = 0; j < W_port_eigen.cols(); j++)
    {

      if (W_port_eigen.col(j).norm() > 0.0)
      {
        CHECK(W_bulk_eigen.col(j).norm() == 0.0);
      }
      else
      {
        // Check there is some overlap in other cols.
        CHECK(W_bulk_eigen.col(j).norm() > 0.0);
      }
    }
  }

  // Now test against port vector. The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e
  // L_e / W_e but with Z_R = 1.0. In this case this = 1.0 for Port 2,3 and != 1.0 for Port
  // 1. See normalization tests of LumpedPort_BasicTests_1ElementPort_Cube321 in
  //    test-lumpedportintegration.cpp.

  for (const auto &[port_i, port_norm] :
       std::vector<std::pair<int, double>>{{1, 1.0}, {2, 1.}})
  {
    ComplexVector port_primary_et;
    space_op.GetLumpedPortExcitationVectorPrimaryEt(port_i, port_primary_et);
    Vector port_primary_et_tmp = port_primary_et.Real();

    W_port->Mult(port_primary_et.Real(), port_primary_et_tmp);
    auto overlap_port =
        linalg::Dot(world_comm, port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_port, WithinRel(port_norm));

    W_bulk->Mult(port_primary_et.Real(), port_primary_et_tmp);
    auto overlap_bulk =
        linalg::Dot(world_comm, port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_bulk, WithinAbs(0.0, 1e-15));

    auto overlap_combined = (*weight_op)
                                .InnerProduct(world_comm, port_primary_et.Real(),
                                              port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_combined, WithinRel(port_norm));
  }

  // // Test actually adding port primary vectors to PROM.
  prom_op.AddLumpedPortModesForSynthesis();
  auto rom_dim = prom_op.GetReducedDimension();
  REQUIRE(rom_dim == 2);
  const auto _norm = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = _norm.L_inv;
  const auto &resistance_R_inv = _norm.R_inv;
  const auto &capacitance_C = _norm.C;
  const auto orth_R = prom_op.GetOrthR();

  CHECK(((inductance_L_inv->rows() == rom_dim) && (inductance_L_inv->cols() == rom_dim)));
  CHECK(((resistance_R_inv->rows() == rom_dim) && (resistance_R_inv->cols() == rom_dim)));
  CHECK(((capacitance_C->rows() == rom_dim) && (capacitance_C->cols() == rom_dim)));
  CHECK(((orth_R.rows() == rom_dim) && (orth_R.cols() == rom_dim)));

  // Per-port load matrices are exported separately from the total loaded matrices so
  // downstream circuit assembly can remove internal terminations and keep external ones.
  // For this two-port lumped case they should have the same shape as the total matrices
  // and contain only the diagonal R/L/C loading terms for the corresponding port row.
  REQUIRE(_norm.port_loads.size() == 2);
  auto find_port_load =
      [&](std::string_view label) -> const RomOperatorTest::NormalizedMatrices::PortLoad *
  {
    for (const auto &load : _norm.port_loads)
    {
      if (load.label == label)
      {
        return &load;
      }
    }
    return nullptr;
  };
  const auto *port1_load = find_port_load("port_1_re");
  const auto *port2_load = find_port_load("port_2_re");
  REQUIRE(port1_load != nullptr);
  REQUIRE(port2_load != nullptr);

  auto check_load_shape = [rom_dim](const auto &load)
  {
    REQUIRE(load.L_inv);
    REQUIRE(load.C);
    CHECK(load.L_inv->rows() == rom_dim);
    CHECK(load.L_inv->cols() == rom_dim);
    CHECK(load.C->rows() == rom_dim);
    CHECK(load.C->cols() == rom_dim);
    if (load.R_inv)
    {
      CHECK(load.R_inv->rows() == rom_dim);
      CHECK(load.R_inv->cols() == rom_dim);
    }
  };
  check_load_shape(*port1_load);
  check_load_shape(*port2_load);

  auto check_only_diagonal =
      [rom_dim](const Eigen::MatrixXcd &mat, long row, std::complex<double> expected)
  {
    for (long i = 0; i < rom_dim; i++)
    {
      for (long j = 0; j < rom_dim; j++)
      {
        if (i == row && j == row)
        {
          CHECK_THAT(mat(i, j).real(), WithinRel(expected.real(), 1e-12));
          CHECK_THAT(mat(i, j).imag(), WithinAbs(expected.imag(), 1e-15));
        }
        else
        {
          CHECK_THAT(std::abs(mat(i, j)), WithinAbs(0.0, 1e-15));
        }
      }
    }
  };
  check_only_diagonal(*port1_load->L_inv, 0, 1.0 / port1_ref_L);
  check_only_diagonal(*port1_load->C, 0, port1_ref_C);
  REQUIRE(port1_load->R_inv);
  check_only_diagonal(*port1_load->R_inv, 0, 1.0 / port1_ref_R);

  CHECK(port2_load->L_inv->cwiseAbs().maxCoeff() == 0.0);
  CHECK(port2_load->C->cwiseAbs().maxCoeff() == 0.0);
  REQUIRE(port2_load->R_inv);
  check_only_diagonal(*port2_load->R_inv, 1, 1.0 / port2_ref_R);

  Eigen::MatrixXd inductance_L_inv_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  inductance_L_inv_ref(0, 0) = 1.0 / port1_ref_L;

  Eigen::MatrixXd resistance_R_inv_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  resistance_R_inv_ref(0, 0) = 1.0 / port1_ref_R;
  resistance_R_inv_ref(1, 1) = 1.0 / port2_ref_R;
  // resistance_R_inv_ref(2, 2) = 1.0 / port3_ref_R;

  Eigen::MatrixXd capacitance_C_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  capacitance_C_ref(0, 0) = port1_ref_C;

  auto check_mat_zero_pattern =
      [dim = prom_op.GetReducedDimension()](auto &mat, bool skip_real = false)
  {
    for (long i = 0; i < dim; i++)
    {
      for (long j = 0; j < dim; j++)
      {
        if (i != j && !skip_real)
        {
          CHECK_THAT(std::real(mat(i, j)), WithinAbs(0.0, 1e-15));
        }
        CHECK_THAT(std::imag(mat(i, j)), WithinAbs(0.0, 1e-15));
      }
    }
  };

  check_mat_zero_pattern(orth_R);
  check_mat_zero_pattern(*resistance_R_inv);
  check_mat_zero_pattern(*capacitance_C);
  check_mat_zero_pattern(*inductance_L_inv,
                         true);  // skip: could extras overlap crowded mesh

  // Port vectors are chosen so that after orthogonalization with the hybrid mass matrix
  // their norm is 1.0.
  CHECK_THAT(std::real(orth_R(0, 0)), WithinRel(1.0));
  CHECK_THAT(std::real(orth_R(1, 1)), WithinRel(1.0));

  CHECK_THAT(std::real((*resistance_R_inv)(0, 0)), WithinRel(resistance_R_inv_ref(0, 0)));
  CHECK_THAT(std::real((*resistance_R_inv)(1, 1)), WithinRel(resistance_R_inv_ref(1, 1)));

  // Don't check 11 of due to presence of bulk values.
  // For 00, bulk values are small as compared to port values, put tolerance of ~1e-6.
  CHECK_THAT(std::real((*inductance_L_inv)(0, 0)),
             WithinRel(inductance_L_inv_ref(0, 0), 1e-6));
  CHECK_THAT(std::real((*capacitance_C)(0, 0)), WithinRel(capacitance_C_ref(0, 0), 1e-6));

  // Add some random vectors to mimic solution of system.
  ComplexVector vector_temp;
  vector_temp.Real() = prom_op.GetVectors().at(0);  // Steal-setup
  vector_temp.Imag() = prom_op.GetVectors().at(0);  // Steal-setup

  std::size_t nr_random_vec = 10;

  std::seed_seq seed_gen{2 * Mpi::Rank(world_comm)};
  std::vector<std::uint32_t> seeds(2 * nr_random_vec);
  seed_gen.generate(seeds.begin(), seeds.end());

  for (std::size_t i = 0; i < nr_random_vec; i++)
  {
    vector_temp.Real().Randomize(seeds.at(2 * i));
    vector_temp.Imag().Randomize(seeds.at(2 * i + 1));

    // The returned coefficients reproduce the vector in the updated basis.
    const auto c = prom_op.UpdatePROM(vector_temp, fmt::format("vec_{}", i));
    ComplexVector w(vector_temp.Size());
    linalg::LinearCombination(prom_op.GetVectors(),
                              {c.data(), static_cast<std::size_t>(c.size())}, w);
    w -= vector_temp;
    CHECK(linalg::Norml2(world_comm, w) <=
          1.0e-12 * linalg::Norml2(world_comm, vector_temp));
  }

  // Check orthogonality with s-matrix vectors.
  const auto &rom_vectors = prom_op.GetVectors();

  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    for (std::size_t i = 0; i < rom_vectors.size() - 2; i++)
    {
      GridFunction E(space_op.GetNDSpace(), true);
      E = 0.0;
      E.Real().SetFromTrueDofs(rom_vectors.at(2 + i));
      E.Real().ExchangeFaceNbrData();
      E.Imag().ExchangeFaceNbrData();

      auto S = port_data.GetSParameter(E);
      CHECK_THAT(std::abs(S), WithinAbs(0.0, 1e-14));
    }
  }
}

// Checks failure mode that neighbouring ports have overlap because they share and edge
// degree of freedom. Ports in ROM must be orthogonal for conventional scattering matrix
// interpretation to be meaningful.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "RomOperator-Synthesis-PortOrthogonality",
                 "[romoperator][Serial]")
{
  MPI_Comm world_comm = Mpi::World();

  size_t order = GENERATE(1UL, 2UL);  // Solver order. Should not matter but set non-default

  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_1_1_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"));

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0},
                                               {"LossTan", 0.0}})})}};

  // Port on neighbouring attributes along same direction will have overlap.
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}}),
                                  json::object({{"Index", 2},
                                                {"R", 50.0},
                                                {"Excitation", uint(0)},
                                                {"Attributes", json::array({102})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  std::size_t max_size_per_excitation = 100;

  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  CHECK_THROWS_WITH(prom_op.AddLumpedPortModesForSynthesis(),
                    Catch::Matchers::ContainsSubstring("should have exactly zero overlap"));
}

// Currently all vectors being added to the PROM basis must be different. The logic of
// dropping degenerate (or almost degenerate) vectors is not implemented. This tests checks
// that adding the same vector twice fails.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-UpdatePROM-LinearDependence",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"Excitation", uint(1)},
            {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                    {"Direction", "+Z"}})})}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  std::size_t max_size_per_excitation = 100;

  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  ComplexVector u;
  u.SetSize(space_op.GetNDSpace().GetTrueVSize());
  u.UseDevice(true);
  u.Real().Randomize(42);
  u.Imag() = 0.0;

  prom_op.UpdatePROM(u, "vec_0");
  CHECK_THROWS(prom_op.UpdatePROM(u, "vec_0_duplicate"));
}

// Tests that the per-port "IncludeInSynthesis" flag suppresses basis vector injection
// without changing the boundary condition. We reuse the cube_3_2_1 two-port setup, but
// flip port 1 (non-excited, RLC) to "IncludeInSynthesis": false, leaving port 2 (excited,
// pure R) included. Expectation: AddLumpedPortModesForSynthesis adds exactly one basis
// vector (port 2), the resulting normalized circuit matrices are 1x1, and the diagonal
// resistance entry equals 1 / port2_ref_R.
//
// This also exercises CalculateNormalizedPROMMatrices in the regime where the number of
// port-mode rows (1) is less than the total lumped port count (2): the port-scaling loop
// must iterate over included ports only, not GetLumpedPortOp().Size(), or it writes past
// the end of v_conc.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-Synthesis-IncludeInSynthesis-Skip",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  const double L0 = 1.0e-6;
  const double Lc = 7.0;

  const double port1_ref_R = 75.0;
  const double port1_ref_C = 5.5e-5;
  const double port1_ref_L = 1.486e-20;
  const double port2_ref_R = 50.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  // Port 1 is the passive RLC port flagged out of synthesis; port 2 stays in.
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array(
           {json::object(
                {{"Index", 1},
                 {"R", port1_ref_R},
                 {"C", port1_ref_C},
                 {"L", port1_ref_L},
                 {"Excitation", uint(0)},
                 {"IncludeInSynthesis", false},
                 {"Elements",
                  json::array({json::object(
                                   {{"Attributes", json::array({1})}, {"Direction", "+Y"}}),
                               json::object({{"Attributes", json::array({9, 11})},
                                             {"Direction", "+Y"}}),
                               json::object({{"Attributes", json::array({2, 6, 10})},
                                             {"Direction", "+X"}})})}}),
            json::object(
                {{"Index", 2},
                 {"R", port2_ref_R},
                 {"Excitation", uint(1)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                         {"Direction", "+Z"}})})}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = 1UL;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  // The flag is parsed onto the config struct.
  const auto &cfg_ports = iodata.boundaries.lumpedport;
  REQUIRE(cfg_ports.size() == 2);
  CHECK(cfg_ports.at(1).include_in_synthesis == false);
  CHECK(cfg_ports.at(2).include_in_synthesis == true);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);
  prom_op.AddLumpedPortModesForSynthesis();

  // Only port 2 contributed a basis vector.
  REQUIRE(prom_op.GetReducedDimension() == 1);

  const auto matrices = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = matrices.L_inv;
  const auto &resistance_R_inv = matrices.R_inv;
  const auto &capacitance_C = matrices.C;
  REQUIRE(inductance_L_inv->rows() == 1);
  REQUIRE(resistance_R_inv->rows() == 1);
  REQUIRE(capacitance_C->rows() == 1);

  // The remaining diagonal entry corresponds to port 2 (R = 50 Ohm, no L, no C).
  CHECK_THAT(std::real((*resistance_R_inv)(0, 0)), WithinRel(1.0 / port2_ref_R));
  CHECK_THAT(std::imag((*resistance_R_inv)(0, 0)), WithinAbs(0.0, 1e-15));
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced postprocessing matches full evaluation",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})},
      {"Postprocessing",
       {{"Energy",
         json::array({json::object(
             {{"Index", 1}, {"Attributes", json::array({1, 2, 3, 4, 5, 6})}})})}}}};
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array(
           {json::object(
                {{"Index", 1},
                 {"R", 50.0},
                 {"L", 1.0e-9},
                 {"C", 1.0e-12},
                 {"Excitation", uint(1)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                         {"Direction", "+Z"}})})}}),
            json::object(
                {{"Index", 2},
                 {"R", 50.0},
                 {"Excitation", uint(2)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({1})},
                                                         {"Direction", "+Y"}})})}})})}};
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3},
        {"AdaptiveCircuitSynthesis", true},
        {"MinFreq", 2.0},
        {"MaxFreq", 32.0},
        {"FreqStep", 30.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperator rom_op(iodata, space_op, 4);
  DrivenPostOperatorTest post_op(iodata, space_op);

  REQUIRE(iodata.solver.driven.sample_f.size() == 2);
  const double omega_front = iodata.solver.driven.sample_f.front();
  const double omega_back = iodata.solver.driven.sample_f.back();
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  linalg::SetRandom(comm, sample.Real(), 42);
  linalg::SetRandom(comm, sample.Imag(), 43);
  rom_op.UpdatePROM(sample, "postprocessing_sample");
  REQUIRE(rom_op.GetReducedDimension() > 0);
  post_op.ConfigureReducedPostprocessing(rom_op, iodata.solver.driven.adaptive_tol);
  REQUIRE(post_op.HasReducedPostprocessing());

  // Construct a deterministic complex field directly in the PROM space. This gives both
  // electric and magnetic forms meaningful magnitude while testing the exact reduced
  // quadratic forms independently of PROM solution accuracy.
  Eigen::VectorXcd y(rom_op.GetReducedDimension());
  for (long i = 0; i < y.size(); i++)
  {
    y(i) = std::complex<double>(1.0 / (i + 1.0), -0.25 / (i + 1.0));
  }
  ComplexVector E(space_op.GetNDSpace().GetTrueVSize());
  E.UseDevice(true);
  E = 0.0;
  for (std::size_t i = 0; i < rom_op.GetBasis().size(); i++)
  {
    linalg::AXPY(y(i).real(), rom_op.GetBasis()[i], E.Real());
    linalg::AXPY(y(i).imag(), rom_op.GetBasis()[i], E.Imag());
  }
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y, rom_op);
  REQUIRE(post_op.HasReducedPostprocessing());

  ComplexVector B(space_op.GetCurlMatrix().Height());
  B.UseDevice(true);
  space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
  space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
  B *= -1.0 / (std::complex<double>(0.0, 1.0) * omega_back);
  post_op.MeasureAndPrintAll(1, 1, E, B, omega_back);
  const Measurement full = post_op.GetMeasurement();
  REQUIRE_FALSE(post_op.ReducedPostprocessingNeedsField(1, rom_op));
  post_op.MeasureAndPrintReduced(1, 1, nullptr, omega_back, y, rom_op);
  const Measurement reduced = post_op.GetMeasurement();
  REQUIRE(post_op.HasReducedPostprocessing());

  auto close = [](double x, double y)
  {
    const double scale = std::max(std::abs(x), std::abs(y));
    return scale == 0.0 || std::abs(x - y) <= 1.0e-10 * scale;
  };
  auto check_complex = [](std::complex<double> x, std::complex<double> y)
  {
    const double scale = std::max(std::abs(x), std::abs(y));
    CHECK((scale == 0.0 || std::abs(x - y) <= 1.0e-10 * scale));
  };

  CAPTURE(reduced.domain_E_field_energy_all, full.domain_E_field_energy_all,
          reduced.domain_H_field_energy_all, full.domain_H_field_energy_all);
  REQUIRE(full.domain_E_field_energy_all > 0.0);
  REQUIRE(full.domain_H_field_energy_all > 0.0);
  CHECK(close(reduced.domain_E_field_energy_all, full.domain_E_field_energy_all));
  CHECK(close(reduced.domain_H_field_energy_all, full.domain_H_field_energy_all));
  REQUIRE(full.domain_E_field_energy_i.size() == 1);
  REQUIRE(full.domain_H_field_energy_i.size() == 1);
  REQUIRE(reduced.domain_E_field_energy_i.size() == full.domain_E_field_energy_i.size());
  REQUIRE(reduced.domain_H_field_energy_i.size() == full.domain_H_field_energy_i.size());
  for (std::size_t i = 0; i < full.domain_E_field_energy_i.size(); i++)
  {
    CHECK(reduced.domain_E_field_energy_i[i].idx == full.domain_E_field_energy_i[i].idx);
    CHECK(close(reduced.domain_E_field_energy_i[i].energy,
                full.domain_E_field_energy_i[i].energy));
    CHECK(close(reduced.domain_E_field_energy_i[i].participation_ratio,
                full.domain_E_field_energy_i[i].participation_ratio));
    CHECK(reduced.domain_H_field_energy_i[i].idx == full.domain_H_field_energy_i[i].idx);
    CHECK(close(reduced.domain_H_field_energy_i[i].energy,
                full.domain_H_field_energy_i[i].energy));
    CHECK(close(reduced.domain_H_field_energy_i[i].participation_ratio,
                full.domain_H_field_energy_i[i].participation_ratio));
  }

  const auto &full_port = full.lumped_port_vi.at(1);
  const auto &reduced_port = reduced.lumped_port_vi.at(1);
  REQUIRE(std::abs(full_port.V) > 0.0);
  REQUIRE(std::abs(full_port.I) > 0.0);
  REQUIRE(full_port.inductor_energy > 0.0);
  REQUIRE(full_port.capacitor_energy > 0.0);
  check_complex(reduced_port.V, full_port.V);
  check_complex(reduced_port.I, full_port.I);
  for (std::size_t i = 0; i < full_port.I_RLC.size(); i++)
  {
    check_complex(reduced_port.I_RLC[i], full_port.I_RLC[i]);
  }
  check_complex(reduced_port.S, full_port.S);
  CHECK(close(reduced_port.inductor_energy, full_port.inductor_energy));
  CHECK(close(reduced_port.capacitor_energy, full_port.capacitor_energy));
  CHECK(close(reduced.lumped_port_inductor_energy, full.lumped_port_inductor_energy));
  CHECK(close(reduced.lumped_port_capacitor_energy, full.lumped_port_capacitor_energy));

  // A NaN in the reduced solution reverts to full postprocessing.
  const double tol = iodata.solver.driven.adaptive_tol;
  post_op.ConfigureReducedPostprocessing(rom_op, tol);
  Eigen::VectorXcd y_nan = y;
  y_nan(0) = std::numeric_limits<double>::quiet_NaN();
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y_nan, rom_op);
  CHECK_FALSE(post_op.HasReducedPostprocessing());

  // Each excitation is checked at its first point: after the first passes, a field
  // inconsistent with the reduced solution of the second reverts.
  post_op.ConfigureReducedPostprocessing(rom_op, tol);
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y, rom_op);
  REQUIRE(post_op.HasReducedPostprocessing());
  REQUIRE(post_op.ReducedPostprocessingNeedsField(2, rom_op));
  E *= 1.001;
  post_op.MeasureAndPrintReduced(2, 0, &E, omega_front, y, rom_op);
  CHECK_FALSE(post_op.HasReducedPostprocessing());
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced postprocessing rejects Floquet wave vector",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path = fs::path(PALACE_TEST_DATA_DIR) /
                         "regression/input/floquet_lumped/mesh/floquet_lumped.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path}, {"L0", 1.0e-2}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1, 2})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"Periodic",
       {{"FloquetWaveVector", json::array({0.1, 0.0, 0.0})},
        {"BoundaryPairs",
         json::array({json::object({{"DonorAttributes", json::array({1})},
                                    {"ReceiverAttributes", json::array({2})}}),
                      json::object({{"DonorAttributes", json::array({3})},
                                    {"ReceiverAttributes", json::array({4})}})})}}},
      {"PEC", {{"Attributes", json::array({8})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"Attributes", json::array({7})},
                                                {"Excitation", uint(1)},
                                                {"R", 377.0},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 5.0}, {"MaxFreq", 10.0}, {"FreqStep", 5.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}}}};

  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  REQUIRE(space_op.GetMaterialOp().HasWaveVector());
  RomOperator rom_op(iodata, space_op, 2);
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  sample.Real().Randomize(42);
  sample.Imag() = 0.0;
  rom_op.UpdatePROM(sample, "wave_vector_sample");
  REQUIRE(rom_op.GetReducedDimension() > 0);

  DrivenPostOperatorTest post_op(iodata, space_op);
  post_op.ConfigureReducedPostprocessing(rom_op, iodata.solver.driven.adaptive_tol);
  CHECK_FALSE(post_op.HasReducedPostprocessing());
}

// Check the PROM system assembled online by RomOperator::SolvePROM, from projected
// frequency-dependent terms and port-space wave port pairings, against the Galerkin
// projection of the HDM system matrix and excitation at the same frequency, built from the
// basis of two random samples.
void CheckPROMAgainstHDM(IoData &iodata, const std::vector<int> &excitations)
{
  MPI_Comm comm = Mpi::World();
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists().back();
  const int size = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector sample(size);
  sample.UseDevice(true);
  for (int seed : {42, 44})
  {
    linalg::SetRandom(comm, sample.Real(), seed);
    linalg::SetRandom(comm, sample.Imag(), seed + 1);
    linalg::SetSubVector(sample, dbc_tdof_list, 0.0);
    rom_op.UpdatePROM(sample, fmt::format("sample_{:d}", seed));
  }
  const auto &V = rom_op.GetBasis();
  const auto n = static_cast<Eigen::Index>(V.size());
  REQUIRE(n > 0);

  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  ComplexVector v(size), w(size), rhs(size), b(size);
  for (auto *x : {&v, &w, &rhs, &b})
  {
    x->UseDevice(true);
  }
  auto project = [&](const ComplexVector &x, Eigen::Index j, Eigen::MatrixXcd &X)
  {
    for (Eigen::Index i = 0; i < n; i++)
    {
      X(i, j) = {V[i] * x.Real(), V[i] * x.Imag()};
    }
  };
  for (double omega :
       {iodata.solver.driven.sample_f.front(), iodata.solver.driven.sample_f.back()})
  {
    for (int excitation_idx : excitations)
    {
      CAPTURE(omega, excitation_idx);
      rom_op.SolvePROM(omega);

      // The port modes are cached at this frequency, so the comparison is not affected by
      // the sign ambiguity of a re-solved eigenvector.
      auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                        std::complex<double>(-omega * omega, 0.0), K.get(),
                                        C.get(), M.get(), A2.get());
      b = 0.0;
      if (space_op.GetExcitationVector1(excitation_idx, rhs))
      {
        b.Add(1i * omega, rhs);
      }
      if (space_op.GetExcitationVector2(excitation_idx, omega, rhs))
      {
        b.Add(1.0, rhs);
      }
      Eigen::MatrixXcd Ar(n, n), br(n, 1);
      for (Eigen::Index j = 0; j < n; j++)
      {
        v.Real() = V[j];
        v.Imag() = 0.0;
        A->Mult(v, w);
        project(w, j, Ar);
      }
      project(b, 0, br);
      Mpi::GlobalSum(static_cast<int>(Ar.size()), Ar.data(), comm);
      Mpi::GlobalSum(static_cast<int>(br.size()), br.data(), comm);

      const double error_A = (rom_op.GetAr() - Ar).cwiseAbs().maxCoeff();
      const Eigen::VectorXcd y = Ar.partialPivLu().solve(br.col(0));
      const double error_y = (rom_op.GetReducedSolution(excitation_idx) - y).norm();
      CAPTURE(error_A, Ar.cwiseAbs().maxCoeff(), error_y, y.norm());
      REQUIRE(br.norm() > 0.0);
      CHECK(error_A <= 1.0e-10 * Ar.cwiseAbs().maxCoeff());
      CHECK(error_y <= 1.0e-8 * y.norm());

      // The wave port S-parameters from the pairings, sᵀV y / R, are those of the expanded
      // solution.
      rom_op.ExpandPROMSolution(excitation_idx, v);
      GridFunction E(space_op.GetNDSpace(), true);
      E.Real().SetFromTrueDofs(v.Real());
      E.Imag().SetFromTrueDofs(v.Imag());
      for (const auto &[idx, data] : space_op.GetWavePortOp())
      {
        const auto *sV = rom_op.GetWavePortModePairing(idx);
        REQUIRE(sV);
        const auto S = sV->cwiseProduct(rom_op.GetReducedSolution(excitation_idx)).sum() /
                       data.modal_reaction;
        const auto S_ref = data.GetSParameter(E);
        CAPTURE(idx, S, S_ref);
        CHECK(std::abs(S - S_ref) <= 1.0e-12 * std::max(std::abs(S_ref), 1.0));
      }
    }
  }
}

// Wave ports: the excitation of an inactive port, sharing the excitation of the active
// port or on its own, and the pairings on a nonconformal hexahedral mesh, with hanging
// dofs on the boundary of the port.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator wave ports match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const unsigned int inactive_excitation = GENERATE(1u, 2u);
  const bool hex = GENERATE(false, true);
  CAPTURE(inactive_excitation, hex);

  // 3 cm x 2 cm x 1 cm box: attributes 14 and 18 are non-adjacent 1 cm x 1 cm faces on the
  // y = 2 cm side, used as an active and an inactive (excited) wave port with PEC on every
  // other face (TE10 above cutoff over the 20-30 GHz band).
  json pec_attributes = json::array();
  for (int attr = 1; attr <= 22; attr++)
  {
    if (attr != 14 && attr != 18)
    {
      pec_attributes.push_back(attr);
    }
  }
  json refinement = json::object();
  if (hex)
  {
    // Refine the element below port 1.
    refinement = {
        {"Boxes", json::array({json::object({{"Levels", 1},
                                             {"BoundingBoxMin", {0.0, 1.9, 0.0}},
                                             {"BoundingBoxMax", {0.1, 2.0, 0.1}}})})}};
  }
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh" /
                   (hex ? "cube_mesh_3_2_1_hex.msh" : "cube_mesh_3_2_1_tet.msh")},
      {"L0", 1.0e-2},
      {"Refinement", refinement},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {
      {"PEC", {{"Attributes", pec_attributes}}},
      {"WavePort", json::array({json::object({{"Index", 1},
                                              {"Attributes", json::array({14})},
                                              {"Mode", 1},
                                              {"Offset", 0.0},
                                              {"Excitation", uint(1)}}),
                                json::object({{"Index", 2},
                                              {"Attributes", json::array({18})},
                                              {"Mode", 1},
                                              {"Offset", 0.0},
                                              {"Active", false},
                                              {"Excitation", inactive_excitation}})})}};
  setup_json["Solver"] = {
      {"Order", 2UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 20.0}, {"MaxFreq", 30.0}, {"FreqStep", 10.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};
  IoData iodata(setup_json, false);
  std::vector<int> excitations = {1};
  if (inactive_excitation != 1)
  {
    excitations.push_back(static_cast<int>(inactive_excitation));
  }
  CheckPROMAgainstHDM(iodata, excitations);
}

// Second-order farfield, surface conductivity, and rational surface impedance boundaries,
// with a lumped port excitation.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator boundary terms match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"},
      {"L0", 1.0e-2},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({14})},
                                                {"Direction", "+Z"}})})},
      {"Absorbing", {{"Attributes", json::array({18, 19, 20})}, {"Order", 2}}},
      {"Conductivity", json::array({json::object({{"Attributes", json::array({13, 15})},
                                                  {"Conductivity", 1.0e3}}),
                                    json::object({{"Attributes", json::array({17})},
                                                  {"Conductivity", 5.8e7},
                                                  {"Thickness", 1.0e-6}})})},
      {"RationalImpedance",
       json::array(
           {json::object({{"Attributes", json::array({21, 22})},
                          {"Numerator", json::array({2.0e-6, 0.0})},
                          {"Denominator", json::array({1.0e-19, 6.0e-9, 377.0})}})})}};
  setup_json["Solver"] = {
      {"Order", 2UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 2.0}, {"MaxFreq", 30.0}, {"FreqStep", 28.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};
  IoData iodata(setup_json, false);
  CheckPROMAgainstHDM(iodata, {1});
}

// Floquet ports, with a frequency-dependent (fixed angle of incidence) or fixed Floquet
// wave vector, together with second-order farfield boundaries or a wave port.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator Floquet ports match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const auto [dir, file, excitations] = GENERATE(
      std::make_tuple("floquet_absorbing", "floquet_absorbing_adaptive.json",
                      std::vector<int>{1}),
      std::make_tuple("floquet_absorbing", "floquet_absorbing_scatterer_adaptive.json",
                      std::vector<int>{1}),
      std::make_tuple("floquet_wave", "floquet_wave.json", std::vector<int>{1, 2}));
  CAPTURE(file);
  const auto input_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / dir;
  std::ifstream f(input_dir / file);
  REQUIRE(f.good());
  json setup_json = json::parse(f, nullptr, true, /*ignore_comments=*/true);
  setup_json["Model"]["Mesh"] =
      (input_dir / setup_json["Model"]["Mesh"].get<std::string>()).string();
  setup_json["Problem"]["Output"] = temp_dir;
  setup_json["Problem"]["Verbose"] = 0;
  setup_json["Solver"]["Order"] = 1;
  if (std::string(file) == "floquet_absorbing_scatterer_adaptive.json")
  {
    // Both polarizations of the incident field.
    setup_json["Boundaries"]["FloquetPort"][0]["IncidentPolarization"] = "RHC";
  }
  IoData iodata(setup_json, false);
  CheckPROMAgainstHDM(iodata, excitations);
}

// Several wave ports can observe different modes of the same boundary, for example
// inactive ports extracting the S-parameters of a second mode at the boundaries of the
// active ports. Circuit synthesis must weight such a shared port boundary once in its inner
// product, as for a single port, instead of aborting on the repeated boundary attribute.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-Synthesis-WavePortSharedBoundary",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";
  const auto domain_orthog =
      GENERATE(as<std::string>{}, "Energy", "SpaceOverlap", "FEBasisIdentity");
  CAPTURE(domain_orthog);

  // The box and port faces of the wave port pairing test above: wave ports on the
  // 1 cm x 1 cm faces 14 (active) and 18 (inactive), with PEC on every other face. The
  // shared configuration adds inactive second-mode ports on both faces.
  json pec_attributes = json::array();
  for (int attr = 1; attr <= 22; attr++)
  {
    if (attr != 14 && attr != 18)
    {
      pec_attributes.push_back(attr);
    }
  }
  auto wave_port = [](int index, int attr, int mode, bool active)
  {
    json port = {{"Index", index},
                 {"Attributes", json::array({attr})},
                 {"Mode", mode},
                 {"Offset", 0.0}};
    if (active)
    {
      port["Excitation"] = uint(1);
    }
    else
    {
      port["Active"] = false;
    }
    return port;
  };

  // Inner products of fixed random vectors under the synthesis weight, and the global
  // number of port true dofs it removes from the domain part.
  auto weights = [&](bool shared)
  {
    json ports = json::array({wave_port(1, 14, 1, true), wave_port(2, 18, 1, false)});
    if (shared)
    {
      ports.push_back(wave_port(3, 14, 2, false));
      ports.push_back(wave_port(4, 18, 2, false));
    }
    json setup_json;
    setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
    setup_json["Model"] = {{"Mesh", mesh_path},
                           {"L0", 1.0e-2},
                           {"Refinement", json::object({})},
                           {"CrackInternalBoundaryElements", false}};
    setup_json["Domains"] = {
        {"Materials",
         json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                    {"Permeability", 1.0},
                                    {"Permittivity", 1.0},
                                    {"LossTan", 0.0}})})}};
    setup_json["Boundaries"] = {{"PEC", {{"Attributes", pec_attributes}}},
                                {"WavePort", ports}};
    setup_json["Solver"] = {
        {"Order", 2UL},
        {"Device", "CPU"},
        {"Driven",
         {{"AdaptiveTol", 1.0e-3},
          {"AdaptiveCircuitSynthesis", true},
          {"AdaptiveCircuitSynthesisDomainOrthogonalization", domain_orthog},
          {"MinFreq", 20.0},
          {"MaxFreq", 30.0},
          {"FreqStep", 10.0}}},
        {"Linear",
         {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

    IoData iodata(setup_json, false);
    auto mesh = LoadScaleParMesh2(iodata, comm);
    SpaceOperator space_op(iodata, mesh);
    RomOperatorTest rom_op(iodata, space_op, 4);
    const auto &weight_op = rom_op.GetWeightOp();
    REQUIRE(weight_op.has_value());

    Vector x(space_op.GetNDSpace().GetTrueVSize()), y(x.Size()), work(x.Size());
    x.UseDevice(true);
    y.UseDevice(true);
    work.UseDevice(true);
    linalg::SetRandom(comm, x, 42);
    linalg::SetRandom(comm, y, 43);
    const std::array<double, 3> dots = {weight_op->InnerProduct(comm, x, x, work),
                                        weight_op->InnerProduct(comm, x, y, work),
                                        weight_op->InnerProduct(comm, y, y, work)};
    int n_port_tdofs = weight_op->port_tdof_list.Size();
    Mpi::GlobalSum(1, &n_port_tdofs, comm);
    return std::make_pair(dots, n_port_tdofs);
  };
  const auto [dots_ref, n_port_tdofs_ref] = weights(false);
  const auto [dots, n_port_tdofs] = weights(true);
  REQUIRE(n_port_tdofs_ref > 0);
  CHECK(n_port_tdofs == n_port_tdofs_ref);
  for (std::size_t i = 0; i < dots.size(); i++)
  {
    CHECK_THAT(dots[i], WithinRel(dots_ref[i], 1.0e-12));
  }
}

// Excited ports must always be included in synthesis. The configuration parser
// rejects a port carrying both "Excitation" > 0 and "IncludeInSynthesis": false because
// the excitation vector is unconditionally added to the basis.
TEST_CASE("RomOperator-Synthesis-ExcludedExcitedRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"IncludeInSynthesis", false},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  CHECK_THROWS_WITH((IoData{setup_json, false}),
                    Catch::Matchers::ContainsSubstring("IncludeInSynthesis"));
}

// Same rule for wave ports: an excited wave port carrying "IncludeInSynthesis": false is
// rejected by the configuration parser, because the excitation vector is unconditionally
// added to the synthesis basis (mirrors the lumped-port case above).
TEST_CASE("RomOperator-Synthesis-WavePortExcludedExcitedRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"WavePort", json::array({json::object({{"Index", 1},
                                              {"Mode", 1},
                                              {"Excitation", uint(1)},
                                              {"IncludeInSynthesis", false},
                                              {"Attributes", json::array({100})}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  CHECK_THROWS_WITH((IoData{setup_json, false}),
                    Catch::Matchers::ContainsSubstring("IncludeInSynthesis"));
}

// A Floquet wave vector makes the periodic problem non-reciprocal, so that the synthesized
// circuit matrices would not be symmetric: the configuration parser rejects circuit
// synthesis with a nonzero wave vector, but not periodic boundaries without a phase shift.
TEST_CASE("RomOperator-Synthesis-FloquetWaveVectorRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  const auto wave_vector =
      GENERATE(json::array({0.0, 0.0, 0.0}), json::array({0.0, 0.1, 0.0}));
  setup_json["Boundaries"] = {
      {"Periodic",
       {{"FloquetWaveVector", wave_vector},
        {"BoundaryPairs",
         json::array({json::object({{"DonorAttributes", json::array({1})},
                                    {"ReceiverAttributes", json::array({2})}})})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  INFO("FloquetWaveVector: " << wave_vector.dump());
  if (wave_vector[1] == 0.0)
  {
    CHECK_NOTHROW((IoData{setup_json, false}));
  }
  else
  {
    CHECK_THROWS_WITH((IoData{setup_json, false}),
                      Catch::Matchers::ContainsSubstring("FloquetWaveVector"));
  }
}

// The polynomial dispersion correction helper is the sign bridge between a scalar
// frequency-dependent boundary multiplier and the quadratic pencil K + iωC - ω²M. Lock
// this down directly: evaluating the corrected pencil must reproduce
// (α₀ + α₁ω + α₂ω²) Mp_r for real coefficients (a lossless wave port) and complex ones (a
// lossy wave port, or a scalar-dispersion fit of another BC), also at complex ω (the
// synthesized eigenvalues).
TEST_CASE("RomOperator-PolynomialFitCorrections", "[romoperator][Serial]")
{
  Eigen::MatrixXcd Mp_r(2, 2);
  Mp_r << std::complex<double>(0.0, 2.0), std::complex<double>(0.0, -0.3),
      std::complex<double>(0.0, -0.3), std::complex<double>(0.0, 0.7);

  for (const auto &[alpha0, alpha1, alpha2] :
       {std::array<std::complex<double>, 3>{1.2, -0.4, 0.05},
        std::array<std::complex<double>, 3>{std::complex<double>(0.4, -0.2),
                                            std::complex<double>(-0.1, 0.3),
                                            std::complex<double>(0.02, -0.04)}})
  {
    Eigen::MatrixXcd Kr = Eigen::MatrixXcd::Zero(2, 2);
    Eigen::MatrixXcd Cr = Eigen::MatrixXcd::Zero(2, 2);
    Eigen::MatrixXcd Mr = Eigen::MatrixXcd::Zero(2, 2);
    RomOperatorTest::ApplyComplexPolynomialFitCorrections(alpha0, alpha1, alpha2, Mp_r, Kr,
                                                          Cr, Mr);
    for (std::complex<double> omega :
         {std::complex<double>(0.2, 0.0), std::complex<double>(1.5, 0.0),
          std::complex<double>(3.0, 0.0), std::complex<double>(1.5, 0.4)})
    {
      const Eigen::MatrixXcd got =
          Kr + std::complex<double>(0.0, 1.0) * omega * Cr - omega * omega * Mr;
      const Eigen::MatrixXcd want =
          (alpha0 + alpha1 * omega + alpha2 * omega * omega) * Mp_r;
      CHECK(((got - want).cwiseAbs().maxCoeff() < 1.0e-13));
    }
  }
}

// Modal residue matrices are real symmetric but generally indefinite. Verify the auxiliary
// direction factorization preserves signed eigenvalues rather than replacing M by |M|.
TEST_CASE("RomOperator-AugmentedPencil-IndefiniteDirections", "[romoperator][Serial]")
{
  Eigen::MatrixXd M(2, 2);
  M << 0.0, 1.0, 1.0, 0.0;  // eigenvalues +1 and -1
  const Eigen::MatrixXcd Mp_r =
      std::complex<double>(0.0, 1.0) * M.cast<std::complex<double>>();
  RomOperatorTest::WavePortAuxBlock blk;
  REQUIRE(RomOperatorTest::AddAuxBlockDirections(blk, Mp_r, 1.0e-12));
  REQUIRE(blk.weights.size() == 2);
  REQUIRE(blk.u_dirs.size() == 2);

  Eigen::MatrixXd reconstructed = Eigen::MatrixXd::Zero(2, 2);
  bool has_positive = false, has_negative = false;
  for (std::size_t j = 0; j < blk.weights.size(); j++)
  {
    reconstructed.noalias() += blk.weights[j] * blk.u_dirs[j] * blk.u_dirs[j].transpose();
    has_positive = has_positive || (blk.weights[j] > 0.0);
    has_negative = has_negative || (blk.weights[j] < 0.0);
  }
  CHECK(has_positive);
  CHECK(has_negative);
  CHECK((reconstructed - M).cwiseAbs().maxCoeff() < 1.0e-13);
}

// Modal-correction residue matrices can contain physical signed directions well below the
// fixed wave-port boundary-mass noise floor. The caller-provided tolerance must control
// their factorization; applying the old 3e-3 mass floor here drops the second direction and
// changes the augmented pencil's Schur complement by O(1e-4).
TEST_CASE("RomOperator-AugmentedPencil-ModalResidueRank", "[romoperator][Serial]")
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, 2);
  M(0, 0) = 1.0;
  M(1, 1) = -1.0e-4;
  const Eigen::MatrixXcd Mp_r =
      std::complex<double>(0.0, 1.0) * M.cast<std::complex<double>>();
  RomOperatorTest::WavePortAuxBlock blk;
  REQUIRE(RomOperatorTest::AddAuxBlockDirections(blk, Mp_r, 1.0e-6));
  REQUIRE(blk.weights.size() == 2);

  Eigen::MatrixXd reconstructed = Eigen::MatrixXd::Zero(2, 2);
  for (std::size_t j = 0; j < blk.weights.size(); j++)
  {
    reconstructed.noalias() += blk.weights[j] * blk.u_dirs[j] * blk.u_dirs[j].transpose();
  }
  CHECK((reconstructed - M).cwiseAbs().maxCoeff() < 1.0e-13);
}

// Schur-complement test for the augmented rational realization. This exercises the real
// BuildAugmentedPencil helper with multiple poles and signed symmetric directions, then
// eliminates the aux rows and checks that the effective v-block exactly equals the
// unaugmented rational pencil contribution i Σ r_k/(ω-p_k) M_proj.
TEST_CASE("RomOperator-AugmentedPencil-SchurComplement", "[romoperator][Serial]")
{
  Eigen::MatrixXcd Kr(2, 2), Cr(2, 2), Mr(2, 2);
  Kr << std::complex<double>(1.0, 0.2), std::complex<double>(0.15, -0.05),
      std::complex<double>(0.15, -0.05), std::complex<double>(1.7, -0.1);
  Cr << std::complex<double>(0.05, -0.02), std::complex<double>(0.01, 0.03),
      std::complex<double>(0.01, 0.03), std::complex<double>(-0.04, 0.01);
  Mr << std::complex<double>(0.3, 0.0), std::complex<double>(0.02, -0.01),
      std::complex<double>(0.02, -0.01), std::complex<double>(0.6, 0.02);

  const double theta = 0.37;
  Eigen::VectorXd u0(2), u1(2);
  u0 << std::cos(theta), std::sin(theta);
  u1 << -std::sin(theta), std::cos(theta);
  RomOperatorTest::WavePortAuxBlock blk;
  blk.port_idx = 3;
  blk.weights = {2.0, -0.4};
  blk.u_dirs = {u0, u1};
  blk.poles = {std::complex<double>(1.1, 0.2), std::complex<double>(-0.7, 0.4)};
  blk.residues = {std::complex<double>(0.3, -0.15), std::complex<double>(-0.2, 0.05)};

  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kr, Cr, Mr, {blk}, aux_labels);
  REQUIRE(aug.Kr.rows() == 6);
  REQUIRE(aug.Kr.cols() == 6);
  REQUIRE(aux_labels.size() == 4);
  CHECK(aux_labels.front() == "waveport_3_p0d0");
  CHECK(aux_labels.back() == "waveport_3_p1d1");

  Eigen::MatrixXcd M_proj =
      blk.weights[0] * (u0 * u0.transpose()).cast<std::complex<double>>() +
      blk.weights[1] * (u1 * u1.transpose()).cast<std::complex<double>>();

  for (std::complex<double> omega :
       {std::complex<double>(0.5, 0.1), std::complex<double>(1.6, -0.2),
        std::complex<double>(2.4, 0.0)})
  {
    const Eigen::MatrixXcd A_aug =
        aug.Kr + std::complex<double>(0.0, 1.0) * omega * aug.Cr - omega * omega * aug.Mr;
    const Eigen::MatrixXcd Avv = A_aug.topLeftCorner(2, 2);
    const Eigen::MatrixXcd Ava = A_aug.topRightCorner(2, 4);
    const Eigen::MatrixXcd Aav = A_aug.bottomLeftCorner(4, 2);
    const Eigen::MatrixXcd Aaa = A_aug.bottomRightCorner(4, 4);
    const Eigen::MatrixXcd schur = Avv - Ava * Aaa.fullPivLu().solve(Aav);

    Eigen::MatrixXcd want =
        Kr + std::complex<double>(0.0, 1.0) * omega * Cr - omega * omega * Mr;
    std::complex<double> rational = 0.0;
    for (std::size_t k = 0; k < blk.poles.size(); k++)
    {
      rational += blk.residues[k] / (omega - blk.poles[k]);
    }
    want += std::complex<double>(0.0, 1.0) * rational * M_proj;
    CHECK(((schur - want).cwiseAbs().maxCoeff() < 1.0e-12));
  }
}

// Verify that the regime-2 augmented pencil produced by BuildAugmentedPencil has the
// same eigenvalues as the unaugmented pencil with the rational kₙ-correction applied
// directly. The augmentation is supposed to be exact: eliminating each aux state sₖⱼ
// from the linear pencil reproduces (i·rₖ·σⱼ·uⱼuⱼᵀ)/(ω − pₖ) in the v-row, so the
// eigenvalues of the augmented pencil are exactly the roots of det(K + iωC − ω²M +
// i·Σₖ rₖ/(ω − pₖ)·M_proj) = 0. This test catches sign errors in the coupling formula
// √(−i·rₖ·σⱼ), wrong K_aa/C_aa diagonal entries, and any off-by-one in the aux row
// layout.
TEST_CASE("RomOperator-AugmentedPencil-Eigenvalues", "[romoperator][Serial]")
{
  // 1×1 base pencil. Treat the basis-row entries as the "v" block: K=1, C=0, M=0, and
  // M_proj = 1·u·uᵀ with u=[1], σ=1. Two poles with hand-picked complex residues so the
  // analytical eigenvalue equation is non-trivial in both real and imaginary parts.
  Eigen::MatrixXcd Kr(1, 1), Cr(1, 1), Mr(1, 1);
  Kr(0, 0) = std::complex<double>(1.0, 0.0);
  Cr(0, 0) = std::complex<double>(0.0, 0.0);
  Mr(0, 0) = std::complex<double>(0.0, 0.0);

  RomOperatorTest::WavePortAuxBlock blk;
  blk.port_idx = 1;
  blk.weights = {1.0};
  Eigen::VectorXd u_dir(1);
  u_dir(0) = 1.0;
  blk.u_dirs = {u_dir};
  blk.poles = {std::complex<double>(2.0, 0.1), std::complex<double>(-1.5, 0.05)};
  blk.residues = {std::complex<double>(0.5, 0.2), std::complex<double>(-0.3, 0.1)};

  std::vector<RomOperatorTest::WavePortAuxBlock> aux_blocks = {blk};
  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kr, Cr, Mr, aux_blocks, aux_labels);

  REQUIRE(aug.Kr.rows() == 3);
  REQUIRE(aug.Kr.cols() == 3);
  REQUIRE(aux_labels.size() == 2);

  // Linearize the QEP into a 6×6 generalized eigenvalue problem and solve.
  // (K + iωC − ω²M)·x = 0  ⇔  [iC, K; -I, 0]·[ωx; x] = ω·[M, 0; 0, -I]·[ωx; x].
  const long n = aug.Kr.rows();
  Eigen::MatrixXcd A_lin = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  Eigen::MatrixXcd B_lin = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  A_lin.topLeftCorner(n, n) = std::complex<double>(0.0, 1.0) * aug.Cr;
  A_lin.topRightCorner(n, n) = aug.Kr;
  A_lin.bottomLeftCorner(n, n) = -Eigen::MatrixXcd::Identity(n, n);
  B_lin.topLeftCorner(n, n) = aug.Mr;
  B_lin.bottomRightCorner(n, n) = -Eigen::MatrixXcd::Identity(n, n);

  // Eigen's GeneralizedEigenSolver works on real matrices only, so we use a dense
  // generalized Schur via ComplexEigenSolver on B⁻¹·A. With M=0 in the v-block, B is
  // singular, so we instead build the explicit polynomial whose roots are the
  // augmented eigenvalues (closed form) and verify the augmented pencil annihilates
  // each.
  //
  // The unaugmented effective pencil is (K + iωC − ω²M) + i·Σₖ rₖ/(ω − pₖ)·M_proj
  // applied at the "v" block (size 1 here). Multiplying by Πₖ(ω − pₖ) gives:
  //   p(ω) = (K + iωC − ω²M)·(ω − p₀)·(ω − p₁) + i·M_proj·[r₀·(ω − p₁) + r₁·(ω − p₀)]
  // For this 1×1 case with K=1, C=0, M=0, M_proj=1:
  //   p(ω) = (ω − p₀)·(ω − p₁) + i·[r₀·(ω − p₁) + r₁·(ω − p₀)]
  // a quadratic in ω. The two roots must equal two of the three augmented-pencil
  // eigenvalues (one extra eigenvalue at infinity from the rank-deficient M-block).
  const auto p0 = blk.poles[0];
  const auto p1 = blk.poles[1];
  const auto r0 = blk.residues[0];
  const auto r1 = blk.residues[1];
  // p(ω) = ω² + b·ω + c with
  //   b = -(p₀ + p₁) + i·(r₀ + r₁)
  //   c =  p₀·p₁ - i·(r₀·p₁ + r₁·p₀)
  const std::complex<double> b = -(p0 + p1) + std::complex<double>(0.0, 1.0) * (r0 + r1);
  const std::complex<double> c =
      p0 * p1 - std::complex<double>(0.0, 1.0) * (r0 * p1 + r1 * p0);
  const std::complex<double> disc = std::sqrt(b * b - 4.0 * c);
  const std::complex<double> root_plus = 0.5 * (-b + disc);
  const std::complex<double> root_minus = 0.5 * (-b - disc);

  // Verify each analytical root annihilates the augmented pencil A_aug(ω) = K + iωC −
  // ω²M when treated as a 3×3 complex matrix (det should vanish).
  auto eval_pencil_det = [&](std::complex<double> w)
  {
    Eigen::MatrixXcd A =
        aug.Kr + std::complex<double>(0.0, 1.0) * w * aug.Cr - w * w * aug.Mr;
    return A.determinant();
  };
  CHECK(std::abs(eval_pencil_det(root_plus)) < 1.0e-10);
  CHECK(std::abs(eval_pencil_det(root_minus)) < 1.0e-10);

  // Also verify the unaugmented effective pencil at the same roots gives the same
  // determinant (= 0) — sanity check on the algebraic equivalence.
  auto eval_unaug = [&](std::complex<double> w)
  {
    std::complex<double> kn = r0 / (w - p0) + r1 / (w - p1);
    return Kr(0, 0) + std::complex<double>(0.0, 1.0) * w * Cr(0, 0) - w * w * Mr(0, 0) +
           std::complex<double>(0.0, 1.0) * kn;
  };
  CHECK(std::abs(eval_unaug(root_plus)) < 1.0e-10);
  CHECK(std::abs(eval_unaug(root_minus)) < 1.0e-10);
}

TEST_CASE("RomOperator-SolveReducedSystem", "[romoperator][Serial]")
{
  // LU solves well and ill-conditioned systems, and full-pivoting QR takes over for
  // singular matrices and for LU solutions with a large backward error in any column.
  using Matrix = Eigen::MatrixXcd;
  using VectorC = Eigen::VectorXcd;
  constexpr int n = 12;
  std::srand(1234);
  auto backward_error = [](const Matrix &A, const Matrix &x, const Matrix &b)
  {
    return (b - A * x).lpNorm<Eigen::Infinity>() /
           (A.cwiseAbs().rowwise().sum().maxCoeff() * x.lpNorm<Eigen::Infinity>() +
            b.lpNorm<Eigen::Infinity>());
  };
  auto with_singular_values = [](const Eigen::VectorXd &s)
  {
    Eigen::HouseholderQR<Matrix> qr_u(Matrix::Random(s.size(), s.size())),
        qr_v(Matrix::Random(s.size(), s.size()));
    const Matrix U = qr_u.householderQ(), W = qr_v.householderQ();
    return Matrix(U * s.cast<std::complex<double>>().asDiagonal() * W.adjoint());
  };
  MPI_Comm comm = Mpi::World();
  const VectorC b = VectorC::Random(n);

  SECTION("Conditioning")
  {
    // Condition numbers 1e10 and 1e17 (singular to working precision): the LU solution is
    // backward stable, accurate to cond ε for the first.
    const double s_min = GENERATE(1.0e-10, 1.0e-17);
    CAPTURE(s_min);
    Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(n, 0.0, -3.0)
                            .unaryExpr([](double e) { return std::pow(10.0, e); });
    s(n - 1) = s_min;
    const Matrix A = with_singular_values(s);
    Matrix x = b;
    CHECK_FALSE(SolveReducedSystem(comm, A, x));
    CHECK(backward_error(A, x, b) < 1.0e-14);
    if (s_min > 1.0e-12)
    {
      const VectorC x_ref = A.fullPivHouseholderQr().solve(b);
      CHECK((x - x_ref).norm() < 1.0e-4 * x_ref.norm());
    }
  }

  SECTION("Singular")
  {
    // With a zero column, the LU solution is not finite and the rank-revealing QR gives a
    // solution of norm ~ |b| / 1e-3, the smallest nonzero singular value.
    Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(n, 0.0, -3.0)
                            .unaryExpr([](double e) { return std::pow(10.0, e); });
    Matrix A = with_singular_values(s);
    A.col(3).setZero();
    Matrix x = b;
    CHECK(SolveReducedSystem(comm, A, x));
    CHECK(x.allFinite());
    CHECK(x.norm() < 1.0e4 * b.norm());
  }

  SECTION("Synthesized PROM slice")
  {
    // A 4 x 4 slice (a port row and three sample rows) of the synthesized L⁻¹, R⁻¹, C of a
    // PROM, assembled over its band: LU is used, and agrees with QR.
    Eigen::Matrix4d Linv, Rinv, Cm;
    Linv << 2998545461200.4326, -34515.101667958588, -10465.943783918958,
        3513.9848630862389, -34515.101667958588, 173468600268174.75, 15792824677.864006,
        -336687502.83115691, -10465.943783918958, 15792824677.864006, 180665882176297.28,
        -15959410866.884502, 3513.9848630862389, -336687502.83115691, -15959410866.884502,
        219574225649943.19;
    Rinv << 0.019999999999999997, -2.1956837108787159e-27, 1.6467627831590372e-27,
        2.058453478948796e-28, -2.1956837108787159e-27, 1.0834394534387061e-12,
        8.4753080346694756e-14, -3.3998897726139405e-15, 1.6467627831590372e-27,
        8.4753080346694756e-14, 6.4383392861610713e-13, -6.4336602257112816e-14,
        2.058453478948796e-28, -3.3998897726139405e-15, -6.4336602257112816e-14,
        1.5291552997921576e-12;
    Cm << 2.0238965694568905e-16, -2.8627695194548223e-25, -8.7168790786791186e-26,
        2.9378257887998199e-26, -2.8627695194548223e-25, 1.7708375637436778e-13,
        -7.4949467469211133e-25, 3.0069730146889297e-26, -8.7168790786791186e-26,
        -7.4949467469211133e-25, 1.7708375648713082e-13, 5.5852088951741535e-25,
        2.9378257887998199e-26, 3.0069730146889297e-26, 5.5852088951741535e-25,
        1.7708375640806658e-13;
    for (double f : {2.5e9, 6.0e9, 9.9e9})
    {
      const double w = 2.0 * std::numbers::pi * f;
      Matrix A = Linv.cast<std::complex<double>>() +
                 (1i * w) * Rinv.cast<std::complex<double>>() -
                 (w * w) * Cm.cast<std::complex<double>>();
      A /= A.cwiseAbs().maxCoeff();
      const VectorC b4 = VectorC::Random(4);
      Matrix x = b4;
      CHECK_FALSE(SolveReducedSystem(comm, A, x));
      const VectorC x_ref = A.fullPivHouseholderQr().solve(b4);
      CHECK((x - x_ref).norm() < 1.0e-12 * x_ref.norm());
    }
  }

  SECTION("Pivot growth")
  {
    // Partial pivoting on Wilkinson's matrix (1 on the diagonal and in the last column, -1
    // below the diagonal) does not swap rows and the entries of the last column of U grow
    // as 2^(k-1), so the LU solution has a large backward error. The matrix itself is well
    // conditioned.
    constexpr int m = 64;
    Matrix A = Matrix::Identity(m, m);
    for (int j = 0; j < m; j++)
    {
      for (int i = j + 1; i < m; i++)
      {
        A(i, j) = -1.0;
      }
    }
    A.col(m - 1).setOnes();
    const VectorC bw = VectorC::Random(m);
    VectorC x_lu = A.partialPivLu().solve(bw);
    REQUIRE(!(backward_error(A, x_lu, bw) < 1.0e-12));
    // A zero first column of the right-hand side does not hide the error of the second.
    Matrix x(m, 2);
    x << VectorC::Zero(m), bw;
    CHECK(SolveReducedSystem(comm, A, x));
    CHECK(x.col(0).isZero());
    CHECK(backward_error(A, x.col(1), bw) < 1.0e-14);
  }
}
