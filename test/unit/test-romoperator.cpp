// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
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
#include "drivers/drivensolver.hpp"
#include "fem/mesh.hpp"
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

class RomOperatorTest : public RomOperator
{
public:
  using RomOperator::CalculateNormalizedPROMMatrices;
  using RomOperator::RomOperator;
  auto &GetWeightOp() const { return weight_op_W; }
  auto &GetOrthR() const { return orth_R; }
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

// Basic checks of ROM construction in the of synthesis. Checks hybrid domain-boundary
// inner-product weight and port overlap. This is a serial test as hex mesh only has a
// single element.
TEST_CASE("RomOperator-Synthesis-Port-Cube111", "[romoperator][Serial]")
{
  // Work with a simple 1x1x1 Cube
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 different test configuration.
  size_t order = GENERATE(1ul, 2ul, 3ul);
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

  double port_ref_R = 50.0;

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      {"PEC", json::object({{"Attributes", json::array({})}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", port_ref_R},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(std::move(setup_json));

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

    // Counting tdofs
    // 3D: There are 6 tets in a cube. Remove edges.
    // 2D: There are 2 * 6 = 12 open faces, so 12 double-counted. Subtract 6 faces without
    // edges. interiors. 1D: add total of 19 edges
    nr_tdof_expected = 6 * (nr_tdof_ref_tet - 6 * order) -
                       12 / 2 * (nr_tdof_ref_tri - 3 * order) + 19 * order;
    nr_tdof_port_expected = 2 * nr_tdof_ref_tri - order;  // 2 triangles, 1 shared edge
  }

  CHECK(W_bulk->NumRows() == nr_tdof_expected);
  CHECK(W_bulk->NumCols() == nr_tdof_expected);

  // Assemble operators as Eigen matrices for simpler testing.
  auto toEigenMatrix = [](const auto &op)
  {
    int n = op.NumRows();
    int m = op.NumCols();
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, m);
    mfem::Vector w(n), v(n);
    w.UseDevice(true);
    v.UseDevice(true);
    for (int i = 0; i < m; i++)
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
  };

  auto W_port_eigen = toEigenMatrix(*W_port);
  auto W_bulk_eigen = toEigenMatrix(*W_bulk);
  auto weight_op_eigen = toEigenMatrix(*weight_op);

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

  // Debug Print.
  // if constexpr (false)
  // {
  //   Eigen::IOFormat HeavyFmt(4, 0, ", ", ";\n", "[", "]", "[", "]");
  //   std::cout << W_bulk_eigen.format(HeavyFmt) << "\n";
  //   std::cout << W_port_eigen.format(HeavyFmt) << "\n";
  //   std::cout << weight_op_eigen.format(HeavyFmt) << "\n";
  // }

  // Now test against port vector. The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e
  // L_e / W_e but with Z_R = 1.0. In this case this = 1.0. See normalization tests of
  // LumpedPort_BasicTests_1ElementPort_Cube321 in test-lumpedportintegration.cpp.
  ComplexVector port_primary_ht_cn;
  space_op.GetLumpedPortExcitationVectorPrimary(1, port_primary_ht_cn, true);
  Vector port_primary_ht_cn_tmp = port_primary_ht_cn.Real();

  W_port->Mult(port_primary_ht_cn.Real(), port_primary_ht_cn_tmp);
  auto overlap_port = port_primary_ht_cn.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_port, WithinRel(1.0));

  W_bulk->Mult(port_primary_ht_cn.Real(), port_primary_ht_cn_tmp);
  auto overlap_bulk = port_primary_ht_cn.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_bulk, WithinAbs(0.0, 1e-15));

  if (Mpi::Size(world_comm) == 1)
  {
    // Rank local overlap.
    auto overlap_combined_local =
        weight_op->InnerProduct(port_primary_ht_cn.Real(), port_primary_ht_cn.Real());
    CHECK_THAT(overlap_combined_local, WithinRel(1.0));
  }

  // Global overlap.
  auto overlap_combined = weight_op->InnerProduct(world_comm, port_primary_ht_cn.Real(),
                                                  port_primary_ht_cn.Real());
  CHECK_THAT(overlap_combined, WithinRel(1.0));

  // Test actually adding port primary vectors to PROM.
  prom_op.AddLumpedPortModesForSynthesis(iodata);
  CHECK(prom_op.GetReducedDimension() == 1);
  const auto [m_Linv, m_Rinv, m_C] = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto orth_R = prom_op.GetOrthR();

  CHECK(((m_Linv->rows() == 1) && (m_Linv->cols() == 1)));
  CHECK(((m_Rinv->rows() == 1) && (m_Rinv->cols() == 1)));
  CHECK(((m_C->rows() == 1) && (m_Rinv->cols() == 1)));
  CHECK(((orth_R.rows() == 1) && (orth_R.cols() == 1)));

  // This should be the same as norm of primary port vector ht_cn with Z_R = 1.0.
  // For single square port this is 1.0
  CHECK_THAT(orth_R(0, 0), WithinRel(1.0));
  CHECK_THAT((*m_Rinv)(0, 0).real(), WithinRel(1.0 / port_ref_R));
  // L, C have some component on port 1 from the bulk. So only check they exist.
  CHECK((*m_Linv)(0, 0).real() != 0.0);
  CHECK((*m_C)(0, 0).real() != 0.0);

  CHECK_THAT((*m_Linv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*m_Rinv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*m_C)(0, 0).imag(), WithinAbs(0.0, 1e-14));
}

// double ht_cn_norm_expected = 0.0;
// for (const auto &el : port_1.elems)
// {
//   ht_cn_norm_expected += el->GetGeometryLength() / el->GetGeometryWidth();
// }
// ht_cn_norm_expected /= port_1.elems.size() * port_1.elems.size();

// double port_ref_C = 5.5e-5;
// double port_ref_L = 1.486e-20;
//                                               {"C", port_ref_C},
//                                               {"L", port_ref_L},
