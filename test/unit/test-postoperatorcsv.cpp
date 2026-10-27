// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

// NOLINTBEGIN(cppcoreguidelines-avoid-do-while)

using namespace palace;

// Boiler-plate using inherit from protected idiom for testing non-public members of
// PostOperatorCSV. Use a fixture test cases as member functions for access.

class PostOperatorCSVManualTest : public PostOperatorCSV<ProblemType::DRIVEN>
{
public:
  using PostOperatorCSV<ProblemType::DRIVEN>::PostOperatorCSV;
  bool DefersTableWrites() const { return defer_table_writes; }
  friend class PostOperatorCSVFixture;
};

IoData load_iodata(std::string_view relative_path)
{
  fs::path io_file(PALACE_TEST_DATA_DIR);
  io_file /= relative_path;
  assert(fs::exists(io_file));
  return IoData{io_file.c_str(), false};
}

std::vector<std::unique_ptr<Mesh>> load_mesh(MPI_Comm &world_comm_, IoData &iodata_)
{
  iodata_.model.mesh = TestMeshPath("fichera-tet.mesh");

  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    auto smesh = mesh::Load(iodata_, world_comm_);
    if (iodata_.model.Lc <= 0.0)
    {
      iodata_.model.Lc = mesh::ComputeReferenceLength(smesh, world_comm_);
    }
    iodata_.NondimensionalizeInputs(smesh);
    mfem_mesh.push_back(mesh::Partition(iodata_, std::move(smesh), world_comm_));
    mesh::RefineMesh(iodata_, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  return mesh_;
}

class PostOperatorCSVFixture
{
public:
  MPI_Comm world_comm;
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;
  SpaceOperator space_op;

  PostOperatorCSVFixture(std::string_view relative_path)
    : world_comm(Mpi::World()), iodata(load_iodata(relative_path)),
      mesh(load_mesh(world_comm, iodata)), space_op(iodata, mesh)
  {
  }

  void restart1_fresh_folder()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_test_tmp";
    REQUIRE(!fs::exists(fs::path(iodata.problem.output) /
                        "port-V.csv"));  // Restart is 1 Indexed.

    // No restart, no previous file to load.
    {
      REQUIRE(iodata.solver.driven.restart == 1);  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 0);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});
        CHECK(post_op_csv.HasSingleExIdx());
        CHECK(!post_op_csv.reload_table);  // Default restart
        REQUIRE_NOTHROW(post_op_csv.InitializePortVI(space_op));
        REQUIRE(post_op_csv.port_V.has_value());
        REQUIRE(post_op_csv.port_I.has_value());
      }
    }

    // Finite restart should fail to init table.
    {
      iodata.solver.driven.restart = 3 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 3);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "simulation expected a restart with existing data"));
      }
    }
  }

  void restart1_restart_in_middle()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_c3";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "port-V.csv"));

    // No restart, no previous file to load.
    {
      iodata.solver.driven.restart = 3 + 1;  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 3);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});
        CHECK(post_op_csv.HasSingleExIdx());
        CHECK(post_op_csv.reload_table);

        REQUIRE_NOTHROW(post_op_csv.InitializePortVI(space_op));
        REQUIRE(post_op_csv.port_V.has_value());
        REQUIRE(post_op_csv.port_I.has_value());

        CHECK(post_op_csv.port_V->table.n_rows() == 3);
        REQUIRE(post_op_csv.port_V->table.n_cols() == 8);

        // Validate column names copied from reference table.
        CHECK(post_op_csv.port_V->table[0].name == "idx");
        CHECK(post_op_csv.port_V->table[1].name == "inc1_1");
        CHECK(post_op_csv.port_V->table[2].name == "re1_1");
        CHECK(post_op_csv.port_V->table[3].name == "im1_1");

        // Name lookups go through the index, which must match the restored names.
        CHECK(post_op_csv.port_V->table.has("inc1_1"));
        CHECK(!post_op_csv.port_V->table.has("col_1"));
        CHECK_NOTHROW(post_op_csv.port_V->table["idx"]);
        CHECK(&post_op_csv.port_V->table["re1_1"] == &post_op_csv.port_V->table[2]);

        // Validate properties copied form reference table.
        CHECK(post_op_csv.port_V->table[0].column_group_idx ==
              -1);  // idx is column block -1
        CHECK(post_op_csv.port_V->table[0].min_left_padding == 0);
        CHECK(post_op_csv.port_V->table[0].float_precision == 8);  // set by PrecIndexCol
        CHECK(post_op_csv.port_V->table[0].fmt_sign == "");
        CHECK(post_op_csv.port_V->table[0].print_as_int == false);

        // Rest of columns should all be in column group 1 (matches excitation idx).
        for (std::size_t i = 1; i < post_op_csv.port_V->table.n_cols(); i++)
        {
          CHECK(post_op_csv.port_V->table[i].column_group_idx == 1);
        }
      }
    }
    // Defaults restart overwrites any existing data and should just work.
    {
      iodata.solver.driven.restart = 0 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 0);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});

        CHECK_NOTHROW(post_op_csv.InitializePortVI(space_op));
      }
    }
    // Different restart should fail to init table.
    {
      iodata.solver.driven.restart = 5 + 1;  // Note: 1 <= restart < nr_total_samples
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 5);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart1_restart_with_empty()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_empty";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "port-V.csv"));

    // No restart, no previous file to load.
    {
      REQUIRE(iodata.solver.driven.restart == 1);  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 0);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1});
        CHECK(post_op_csv.HasSingleExIdx());
        CHECK(!post_op_csv.reload_table);  // Default restart

        REQUIRE_NOTHROW(post_op_csv.InitializePortVI(space_op));
        REQUIRE(post_op_csv.port_V.has_value());
        REQUIRE(post_op_csv.port_I.has_value());

        CHECK(post_op_csv.port_V->table.n_rows() == 0);
        REQUIRE(post_op_csv.port_V->table.n_cols() == 8);
      }
    }
    // Finite restart should fail to init table.
    {
      iodata.solver.driven.restart = 3 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart1_mismatch_col_nr()
  {
    iodata.solver.driven.restart = 3;  // non-trivial restart for check to trigger
                                       // Try and load wrong table with incorrect width.
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart2_c03";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "port-V.csv"));

    PostOperatorCSVManualTest post_op_csv{iodata, space_op};
    if (Mpi::Root(Mpi::World()))
    {
      CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                        Catch::Matchers::ContainsSubstring("Mismatched number of columns"));
    }
  }

  void restart1_mismatch_col_headers()
  {
    iodata.solver.driven.restart = 3 + 1;
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_colswap";
    PostOperatorCSVManualTest post_op_csv{iodata, space_op};
    if (Mpi::Root(Mpi::World()))
    {
      CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                        Catch::Matchers::ContainsSubstring("Mismatched column header"));
    }
  }

  void restart1_bad_col_alignment()
  {
    iodata.solver.driven.restart = 3 + 1;
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_badcols";
    PostOperatorCSVManualTest post_op_csv{iodata, space_op};
    if (Mpi::Root(Mpi::World()))
    {
      CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                        Catch::Matchers::ContainsSubstring("Mismatched rows"));
    }
  }

  void restart2_restart_in_middle_ex1()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart2_c03";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "port-V.csv"));

    // No restart, no previous file to load.
    {
      iodata.solver.driven.restart = 3 + 1;  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 3);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1, 2});
        CHECK(!post_op_csv.HasSingleExIdx());
        CHECK(post_op_csv.reload_table);

        REQUIRE_NOTHROW(post_op_csv.InitializePortVI(space_op));
        REQUIRE(post_op_csv.port_V.has_value());
        REQUIRE(post_op_csv.port_I.has_value());

        CHECK(post_op_csv.port_V->table.n_rows() == 3);
        REQUIRE(post_op_csv.port_V->table.n_cols() == 15);

        // Validate column names copied from reference table.
        CHECK(post_op_csv.port_V->table[0].name == "idx");
        CHECK(post_op_csv.port_V->table[1].name == "inc1_1");  // Port 1 hosts Excitation 1
        CHECK(post_op_csv.port_V->table[2].name == "re1_1");
        CHECK(post_op_csv.port_V->table[3].name == "im1_1");
        CHECK(post_op_csv.port_V->table[8].name == "inc2_2");  // Port 2 hosts Excitation 2
        CHECK(post_op_csv.port_V->table[9].name == "re1_2");
        CHECK(post_op_csv.port_V->table[10].name == "im1_2");

        // Validate properties copied form reference table.
        CHECK(post_op_csv.port_V->table[0].column_group_idx ==
              -1);  // idx is column block -1
        CHECK(post_op_csv.port_V->table[0].min_left_padding == 0);
        CHECK(post_op_csv.port_V->table[0].float_precision == 8);  // set by PrecIndexCol
        CHECK(post_op_csv.port_V->table[0].fmt_sign == "");
        CHECK(post_op_csv.port_V->table[0].print_as_int == false);

        for (std::size_t i = 1; i < 8; i++)
        {
          CHECK(post_op_csv.port_V->table[i].column_group_idx == 1);
          CHECK(post_op_csv.port_V->table[i].n_rows() == 3);
        }
        for (std::size_t i = 8; i < post_op_csv.port_V->table.n_cols(); i++)
        {
          CHECK(post_op_csv.port_V->table[i].column_group_idx == 2);
          CHECK(post_op_csv.port_V->table[i].n_rows() == 0);
        }
      }
    }
    // Different restart should fail to init table.
    {
      iodata.solver.driven.restart = 0 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 0);
        CHECK(post_op_csv.ex_idx_i == 0);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1, 2});

        // No throw since 1 restart just overwrite existing table.
        CHECK_NOTHROW(post_op_csv.InitializePortVI(space_op));
      }
    }
    // Different restart should fail to init table.
    {
      iodata.solver.driven.restart = 7 + 1;  // Note: 1 <= restart < nr_total_samples
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 1);
        CHECK(post_op_csv.ex_idx_i == 1);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1, 2});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart2_restart_in_middle_ex2()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart2_c14";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "port-V.csv"));

    // No restart, no previous file to load.
    {
      iodata.solver.driven.restart = 6 + 4 + 1;  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {

        // Check cursor is at zero.
        CHECK(post_op_csv.row_i == 4);
        CHECK(post_op_csv.ex_idx_i == 1);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1, 2});
        CHECK(!post_op_csv.HasSingleExIdx());
        CHECK(post_op_csv.reload_table);

        REQUIRE_NOTHROW(post_op_csv.InitializePortVI(space_op));
        REQUIRE(post_op_csv.port_V.has_value());
        REQUIRE(post_op_csv.port_I.has_value());

        CHECK(post_op_csv.port_V->table.n_rows() == 6);
        REQUIRE(post_op_csv.port_V->table.n_cols() == 15);

        for (std::size_t i = 1; i < 8; i++)
        {
          CHECK(post_op_csv.port_V->table[i].column_group_idx == 1);
          CHECK(post_op_csv.port_V->table[i].n_rows() == 6);
        }
        for (std::size_t i = 8; i < post_op_csv.port_V->table.n_cols(); i++)
        {
          CHECK(post_op_csv.port_V->table[i].column_group_idx == 2);
          CHECK(post_op_csv.port_V->table[i].n_rows() == 4);
        }
      }
    }
    // Restart should not throw since it overwrites table.
    {
      iodata.solver.driven.restart = 0 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_NOTHROW(post_op_csv.InitializePortVI(space_op));
      }
    }
    // Different restart should fail to init table.
    {
      iodata.solver.driven.restart = 3 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }
  void restart3_farfield_reload()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart3_c63";
    REQUIRE(fs::exists(fs::path(iodata.problem.output) / "farfield-rE.csv"));

    SurfacePostOperator surf_post_op(iodata, space_op.GetMaterialOp(),
                                     space_op.GetH1Space(), space_op.GetNDSpace());
    REQUIRE(surf_post_op.farfield.size() == 2);

    // Restart after the 9th of the 12 measurements (excitation 1 complete, excitation 2
    // at 3 of 6 frequencies). The far-field table gains one row block per measurement,
    // so every column group must hold 9 * 2 = 18 rows.
    {
      iodata.solver.driven.restart = 6 + 3 + 1;  // Restart is 1 Indexed.
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};

      if (Mpi::Root(Mpi::World()))
      {
        CHECK(post_op_csv.row_i == 3);
        CHECK(post_op_csv.ex_idx_i == 1);
        CHECK(post_op_csv.nr_expected_measurement_rows == 6);
        CHECK(post_op_csv.ex_idx_v_all == std::vector<std::size_t>{1, 2});
        CHECK(post_op_csv.reload_table);

        REQUIRE_NOTHROW(post_op_csv.InitializeFarFieldE(surf_post_op));
        REQUIRE(post_op_csv.farfield_E.has_value());
        CHECK(post_op_csv.farfield_E->table.n_rows() == 18);
        REQUIRE(post_op_csv.farfield_E->table.n_cols() == 10);

        // Validate column names copied from the reference table.
        CHECK(post_op_csv.farfield_E->table[0].name == "idx");
        CHECK(post_op_csv.farfield_E->table[1].name == "exc");
        CHECK(post_op_csv.farfield_E->table[2].name == "theta");
        CHECK(post_op_csv.farfield_E->table[3].name == "phi");
        CHECK(post_op_csv.farfield_E->table[4].name == "rE0_re");
        CHECK(post_op_csv.farfield_E->table[9].name == "rE2_im");

        // Name lookups go through the index, which must match the restored names.
        CHECK(post_op_csv.farfield_E->table.has("rE2_im"));
        CHECK(!post_op_csv.farfield_E->table.has("col_1"));
        CHECK_NOTHROW(post_op_csv.farfield_E->table["theta"]);
        CHECK(&post_op_csv.farfield_E->table["exc"] == &post_op_csv.farfield_E->table[1]);

        // Validate properties copied from the reference table: the idx column is its
        // own group, all other columns share group 0 and the same 18 rows.
        CHECK(post_op_csv.farfield_E->table[0].column_group_idx == -1);
        CHECK(post_op_csv.farfield_E->table[1].print_as_int == true);
        for (std::size_t i = 1; i < post_op_csv.farfield_E->table.n_cols(); i++)
        {
          CHECK(post_op_csv.farfield_E->table[i].column_group_idx == 0);
          CHECK(post_op_csv.farfield_E->table[i].n_rows() == 18);
        }
      }
    }
    // Default restart overwrites any existing data and should just work.
    {
      iodata.solver.driven.restart = 0 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_NOTHROW(post_op_csv.InitializeFarFieldE(surf_post_op));
      }
    }
    // A restart position not matching the reloaded rows must be rejected, whether it
    // disagrees on the frequency or on the completed excitation.
    {
      iodata.solver.driven.restart = 3 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_THROWS_WITH(post_op_csv.InitializeFarFieldE(surf_post_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
    {
      iodata.solver.driven.restart = 6 + 4 + 1;
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      if (Mpi::Root(Mpi::World()))
      {
        CHECK_THROWS_WITH(post_op_csv.InitializeFarFieldE(surf_post_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart3_farfield_missing_data()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart3_test_tmp";
    REQUIRE(!fs::exists(fs::path(iodata.problem.output) / "farfield-rE.csv"));

    iodata.solver.driven.restart = 3 + 1;  // Restart is 1 Indexed.
    SurfacePostOperator surf_post_op(iodata, space_op.GetMaterialOp(),
                                     space_op.GetH1Space(), space_op.GetNDSpace());
    PostOperatorCSVManualTest post_op_csv{iodata, space_op};
    if (Mpi::Root(Mpi::World()))
    {
      CHECK_THROWS_WITH(post_op_csv.InitializeFarFieldE(surf_post_op),
                        Catch::Matchers::ContainsSubstring(
                            "simulation expected a restart with existing data"));
    }
  }
};

TEST_CASE("PostOperatorCSV_Restart_Helper_ExpectedFilling", "[postoperatorcsv][Serial]")
{
  using vs = std::vector<std::size_t>;
  CHECK(_impl::table_expected_filling(0, 0, 5, 1) == vs{0, 0});
  CHECK(_impl::table_expected_filling(3, 0, 5, 1) == vs{3, 3});

  CHECK(_impl::table_expected_filling(3, 0, 5, 3) == vs{3, 3, 0, 0});
  CHECK(_impl::table_expected_filling(0, 1, 5, 3) == vs{5, 5, 0, 0});
  CHECK(_impl::table_expected_filling(3, 1, 5, 3) == vs{5, 5, 3, 0});
  CHECK(_impl::table_expected_filling(2, 2, 5, 3) == vs{5, 5, 5, 2});
}

// Could also use METHOD_AS_TEST_CASE with fixture. This reuses fixture class, *without*
// inheriting from it, so friendship assignment works.
TEST_CASE("PostOperatorCSV_Restart_OneExcitation", "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture reload_fixture{"postoperatorcsv_restart/restart1.json"};
  SECTION("Driven solver, single excitation: no reload")
  {
    reload_fixture.restart1_fresh_folder();
  }
  SECTION("Driven solver, single excitation: load with restart")
  {
    reload_fixture.restart1_restart_in_middle();
  }
  SECTION("Driven solver, single excitation: load empty table")
  {
    reload_fixture.restart1_restart_with_empty();
  }
  SECTION("Driven solver, single excitation: load mismatch table col nr")
  {
    reload_fixture.restart1_mismatch_col_nr();
  }
  SECTION("Driven solver, single excitation: load mismatch table header")
  {
    reload_fixture.restart1_mismatch_col_headers();
  }
  SECTION("Driven solver, single excitation: bad column data in csv file")
  {
    reload_fixture.restart1_bad_col_alignment();
  }
}
TEST_CASE("PostOperatorCSV_AdaptiveOutputBuffering", "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture fixture{"postoperatorcsv_restart/restart1.json"};
  fixture.iodata.solver.driven.adaptive_tol = 1.0e-3;
  PostOperatorCSVManualTest post_op_csv{fixture.iodata, fixture.space_op};
  if (Mpi::Root(Mpi::World()))
  {
    CHECK(post_op_csv.DefersTableWrites());
  }
}

TEST_CASE("PostOperatorCSV_Restart_TwoExcitation", "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture fixture{"postoperatorcsv_restart/restart2.json"};

  SECTION("Driven solver, two excitations: load with restart in ex 1")
  {
    fixture.restart2_restart_in_middle_ex1();
  }

  SECTION("Driven solver, two excitations: load with restart in ex 2")
  {
    fixture.restart2_restart_in_middle_ex2();
  }
}

TEST_CASE("PostOperatorCSV_Restart_FarField", "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture fixture{"postoperatorcsv_restart/restart3.json"};

  SECTION("Driven solver, far-field table: load with restart")
  {
    fixture.restart3_farfield_reload();
  }

  SECTION("Driven solver, far-field table: no data to restart from")
  {
    fixture.restart3_farfield_missing_data();
  }
}

// NOLINTEND(cppcoreguidelines-avoid-do-while)
