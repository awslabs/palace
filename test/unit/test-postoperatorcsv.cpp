// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
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
  iodata_.model.mesh = fs::path(PALACE_TEST_DATA_DIR) / "mesh/fichera-tet.mesh";

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

  void surface_edge_table(const fs::path &output)
  {
    iodata.problem.output = output;
    const auto csv_path = output / "surface-Q-edge.csv";

    // Write one frequency sample containing two edge diagnostics.
    iodata.solver.driven.restart = 1;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      REQUIRE_NOTHROW(post_op_csv.InitializeSurfaceQEdge(2));
      REQUIRE(post_op_csv.surface_Q_edge.has_value());
      REQUIRE(post_op_csv.surface_Q_edge->table.n_cols() == 8);
      CHECK(post_op_csv.surface_Q_edge->table.n_rows() == 0);
      CHECK(post_op_csv.surface_Q_edge->table[0].header_text == "f (GHz)");
      CHECK(post_op_csv.surface_Q_edge->table[1].header_text == "exc");
      CHECK(post_op_csv.surface_Q_edge->table[2].header_text == "interface");
      CHECK(post_op_csv.surface_Q_edge->table[3].header_text == "R (m)");
      CHECK(post_op_csv.surface_Q_edge->table[4].header_text == "E_out (J)");
      CHECK(post_op_csv.surface_Q_edge->table[5].header_text == "p_out");
      CHECK(post_op_csv.surface_Q_edge->table[6].header_text == "E_ann (J)");
      CHECK(post_op_csv.surface_Q_edge->table[7].header_text == "p_ann");

      post_op_csv.row_idx_v = 2.0;
      post_op_csv.measurement_cache.interface_edge_i = {
          {1, 1.0e-6, 2.0e-18, 3.0e-6, 4.0e-18, 5.0e-6},
          {2, 2.0e-6, 6.0e-18, 7.0e-6, 8.0e-18, 9.0e-6}};
      REQUIRE_NOTHROW(post_op_csv.PrintSurfaceQEdge());
    }

    TableWithCSVFile written(csv_path.string(), true);
    REQUIRE(written.table.n_cols() == 8);
    REQUIRE(written.table.n_rows() == 2);
    CHECK(written.table[0].data == std::vector<double>{2.0, 2.0});
    CHECK(written.table[1].data == std::vector<double>{1.0, 1.0});
    CHECK(written.table[2].data == std::vector<double>{1.0, 2.0});

    // A restart at the next sample accepts exactly the two existing long-form rows.
    iodata.solver.driven.restart = 2;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      REQUIRE_NOTHROW(post_op_csv.InitializeSurfaceQEdge(2));
      REQUIRE(post_op_csv.surface_Q_edge.has_value());
      CHECK(post_op_csv.surface_Q_edge->table.n_rows() == 2);
      CHECK(post_op_csv.surface_Q_edge->table[1].name == "exc");
      CHECK(post_op_csv.surface_Q_edge->table[1].print_as_int);
      CHECK(post_op_csv.surface_Q_edge->table[2].print_as_int);
    }

    // A later restart expects another complete pair of rows.
    iodata.solver.driven.restart = 3;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      CHECK_THROWS_WITH(
          post_op_csv.InitializeSurfaceQEdge(2),
          Catch::Matchers::ContainsSubstring(
              "has 2 rows, but 4 rows were expected at the restart position"));
    }

    // Header changes are rejected even when the row count is compatible.
    written.table[4].header_text = "wrong energy header";
    written.WriteFullTableTrunc();
    iodata.solver.driven.restart = 2;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      CHECK_THROWS_WITH(
          post_op_csv.InitializeSurfaceQEdge(2),
          Catch::Matchers::ContainsSubstring("has an incompatible column header"));
    }
  }

  void surface_local_edge_table(const fs::path &output)
  {
    iodata.problem.output = output;
    const auto csv_path = output / "surface-Q-edge-local.csv";

    iodata.solver.driven.restart = 1;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      REQUIRE_NOTHROW(post_op_csv.InitializeSurfaceQEdgeLocal(2));
      REQUIRE(post_op_csv.surface_Q_edge_local.has_value());
      REQUIRE(post_op_csv.surface_Q_edge_local->table.n_cols() == 56);
      CHECK(post_op_csv.surface_Q_edge_local->table.n_rows() == 0);
      CHECK(post_op_csv.surface_Q_edge_local->table[0].header_text == "f (GHz)");
      CHECK(post_op_csv.surface_Q_edge_local->table[1].header_text == "exc");
      CHECK(post_op_csv.surface_Q_edge_local->table[2].header_text == "interface");
      CHECK(post_op_csv.surface_Q_edge_local->table[3].header_text == "edge");
      CHECK(post_op_csv.surface_Q_edge_local->table[4].header_text == "R (m)");
      CHECK(post_op_csv.surface_Q_edge_local->table[11].header_text == "L_edge (m)");
      CHECK(post_op_csv.surface_Q_edge_local->table[12].header_text == "E_total (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[13].header_text == "p_total");
      CHECK(post_op_csv.surface_Q_edge_local->table[14].header_text == "E_in (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[15].header_text == "p_in");
      CHECK(post_op_csv.surface_Q_edge_local->table[16].header_text == "E_ann (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[17].header_text == "p_ann");
      CHECK(post_op_csv.surface_Q_edge_local->table[18].header_text == "E_bulk_ann (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[19].header_text == "p_bulk_ann");
      CHECK(post_op_csv.surface_Q_edge_local->table[20].header_text == "E_vertex_in (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[21].header_text == "p_vertex_in");
      CHECK(post_op_csv.surface_Q_edge_local->table[22].header_text ==
            "E_bulk_vertex_ann (J)");
      CHECK(post_op_csv.surface_Q_edge_local->table[23].header_text == "p_bulk_vertex_ann");
      CHECK(post_op_csv.surface_Q_edge_local->table[25].header_text == "p_total_n");
      CHECK(post_op_csv.surface_Q_edge_local->table[31].header_text == "p_in_t");
      CHECK(post_op_csv.surface_Q_edge_local->table[33].header_text == "p_ann_n");
      CHECK(post_op_csv.surface_Q_edge_local->table[37].header_text == "p_bulk_top_n_ann");
      CHECK(post_op_csv.surface_Q_edge_local->table[47].header_text ==
            "p_bulk_bottom_l_ann");
      CHECK(post_op_csv.surface_Q_edge_local->table[48].header_text == "automatic");
      CHECK(post_op_csv.surface_Q_edge_local->table[49].header_text == "component");
      CHECK(post_op_csv.surface_Q_edge_local->table[50].header_text == "chain");
      CHECK(post_op_csv.surface_Q_edge_local->table[51].header_text == "v0_type");
      CHECK(post_op_csv.surface_Q_edge_local->table[52].header_text == "v1_type");
      CHECK(post_op_csv.surface_Q_edge_local->table[53].header_text == "process_nx");
      CHECK(post_op_csv.surface_Q_edge_local->table[55].header_text == "process_nz");

      post_op_csv.row_idx_v = 2.0;
      post_op_csv.measurement_cache.interface_local_edge_i = {
          {1,
           1,
           true,
           2,
           3,
           {1, 2},
           {0.0, 0.0, 1.0},
           {1.0e-6, 2.0e-6, 3.0e-6},
           {4.0e-6, 5.0e-6, 6.0e-6},
           7.0e-6,
           8.0e-6,
           9.0e-18,
           1.0e-5,
           {31.0e-18, 33.0e-18},
           {3.2e-5, 3.4e-5},
           11.0e-18,
           1.2e-5,
           {35.0e-18, 37.0e-18},
           {3.6e-5, 3.8e-5},
           13.0e-18,
           1.4e-5,
           {39.0e-18, 41.0e-18},
           {4.0e-5, 4.2e-5},
           14.0e-18,
           1.5e-5,
           15.0e-18,
           1.6e-5,
           16.0e-18,
           1.7e-5,
           {43.0e-18, 45.0e-18, 47.0e-18, 49.0e-18, 51.0e-18, 53.0e-18},
           {4.4e-5, 4.6e-5, 4.8e-5, 5.0e-5, 5.2e-5, 5.4e-5}},
          {2,
           3,
           false,
           -1,
           -1,
           {-1, -1},
           {0.0, 1.0, 0.0},
           {11.0e-6, 12.0e-6, 13.0e-6},
           {14.0e-6, 15.0e-6, 16.0e-6},
           17.0e-6,
           18.0e-6,
           19.0e-18,
           2.0e-5,
           {131.0e-18, 133.0e-18},
           {13.2e-5, 13.4e-5},
           21.0e-18,
           2.2e-5,
           {135.0e-18, 137.0e-18},
           {13.6e-5, 13.8e-5},
           23.0e-18,
           2.4e-5,
           {139.0e-18, 141.0e-18},
           {14.0e-5, 14.2e-5},
           24.0e-18,
           2.5e-5,
           25.0e-18,
           2.6e-5,
           26.0e-18,
           2.7e-5,
           {143.0e-18, 145.0e-18, 147.0e-18, 149.0e-18, 151.0e-18, 153.0e-18},
           {14.4e-5, 14.6e-5, 14.8e-5, 15.0e-5, 15.2e-5, 15.4e-5}}};
      REQUIRE_NOTHROW(post_op_csv.PrintSurfaceQEdgeLocal());
    }

    TableWithCSVFile written(csv_path.string(), true);
    REQUIRE(written.table.n_cols() == 56);
    REQUIRE(written.table.n_rows() == 2);
    CHECK(written.table[0].data == std::vector<double>{2.0, 2.0});
    CHECK(written.table[2].data == std::vector<double>{1.0, 2.0});
    CHECK(written.table[3].data == std::vector<double>{1.0, 3.0});
    CHECK(written.table[11].data == std::vector<double>{7.0e-6, 17.0e-6});
    CHECK(written.table[12].data == std::vector<double>{9.0e-18, 19.0e-18});
    CHECK(written.table[14].data == std::vector<double>{11.0e-18, 21.0e-18});
    CHECK(written.table[18].data == std::vector<double>{15.0e-18, 25.0e-18});
    CHECK(written.table[20].data == std::vector<double>{14.0e-18, 24.0e-18});
    CHECK(written.table[22].data == std::vector<double>{16.0e-18, 26.0e-18});
    CHECK(written.table[24].data == std::vector<double>{31.0e-18, 131.0e-18});
    CHECK(written.table[28].data == std::vector<double>{35.0e-18, 135.0e-18});
    CHECK(written.table[32].data == std::vector<double>{39.0e-18, 139.0e-18});
    CHECK(written.table[36].data == std::vector<double>{43.0e-18, 143.0e-18});
    CHECK(written.table[46].data == std::vector<double>{53.0e-18, 153.0e-18});
    CHECK(written.table[48].data == std::vector<double>{1.0, 0.0});
    CHECK(written.table[49].data == std::vector<double>{2.0, -1.0});
    CHECK(written.table[50].data == std::vector<double>{3.0, -1.0});
    CHECK(written.table[51].data == std::vector<double>{1.0, -1.0});
    CHECK(written.table[52].data == std::vector<double>{2.0, -1.0});
    CHECK(written.table[53].data == std::vector<double>{0.0, 0.0});
    CHECK(written.table[54].data == std::vector<double>{0.0, 1.0});
    CHECK(written.table[55].data == std::vector<double>{1.0, 0.0});

    iodata.solver.driven.restart = 2;
    {
      PostOperatorCSVManualTest post_op_csv{iodata, space_op};
      REQUIRE_NOTHROW(post_op_csv.InitializeSurfaceQEdgeLocal(2));
      REQUIRE(post_op_csv.surface_Q_edge_local.has_value());
      CHECK(post_op_csv.surface_Q_edge_local->table.n_rows() == 2);
      CHECK(post_op_csv.surface_Q_edge_local->table[3].print_as_int);
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

TEST_CASE_METHOD(palace::test::PerRankTempDir, "PostOperatorCSV_SurfaceQEdge",
                 "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture fixture{"postoperatorcsv_restart/restart1.json"};
  fixture.surface_edge_table(temp_dir);
}

TEST_CASE_METHOD(palace::test::PerRankTempDir, "PostOperatorCSV_SurfaceQEdgeLocal",
                 "[postoperatorcsv][Serial]")
{
  PostOperatorCSVFixture fixture{"postoperatorcsv_restart/restart1.json"};
  fixture.surface_local_edge_table(temp_dir);
}

// NOLINTEND(cppcoreguidelines-avoid-do-while)
