#include <complex>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "fem/mesh.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"

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
  fs::path io_file(PALACE_TEST_DIR);
  io_file /= relative_path;
  assert(fs::exists(io_file));
  return {io_file.c_str(), false};
}

std::vector<std::unique_ptr<Mesh>> load_mesh(MPI_Comm &world_comm_, IoData &iodata_)
{
  iodata_.model.mesh = fs::path(PALACE_TEST_DIR) / "mesh/fichera-tet.mesh";

  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::ReadMesh(iodata_, world_comm_));
    iodata_.NondimensionalizeInputs(*mfem_mesh[0]);
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
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_test_tmp";
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "simulation expected a restart with existing data"));
      }
    }
  }

  void restart1_restart_in_middle()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_c3";
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});

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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart1_restart_with_empty()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_empty";
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1});
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
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart2_c03";
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
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_colswap";
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
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_badcols";
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
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart2_c03";
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1, 2});
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1, 2});

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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1, 2});

        CHECK_THROWS_WITH(post_op_csv.InitializePortVI(space_op),
                          Catch::Matchers::ContainsSubstring(
                              "Specified restart position is incompatible with reloaded"));
      }
    }
  }

  void restart2_restart_in_middle_ex2()
  {
    iodata.problem.output =
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart2_c14";
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
        CHECK(post_op_csv.ex_idx_v_all == std::vector<size_t>{1, 2});
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
};

TEST_CASE("PostOperatorCSV_Restart_Helper_ExpectedFilling", "[postoperatorcsv]")
{
  using vs = std::vector<std::size_t>;
  CHECK(_impl::table_expected_filling(0, 0, 5, 1) == vs{0, 0});
  CHECK(_impl::table_expected_filling(3, 0, 5, 1) == vs{3, 3});

  CHECK(_impl::table_expected_filling(3, 0, 5, 3) == vs{3, 3, 0, 0});
  CHECK(_impl::table_expected_filling(0, 1, 5, 3) == vs{5, 5, 0, 0});
  CHECK(_impl::table_expected_filling(3, 1, 5, 3) == vs{5, 5, 3, 0});
  CHECK(_impl::table_expected_filling(2, 2, 5, 3) == vs{5, 5, 5, 2});
}

// Could also use METHOD_AS_TEST_CASE with fixture. This resuses fixture class, *without*
// inheriting from it, so friendship assignment works.
TEST_CASE("PostOperatorCSV_Restart_OneExcitation", "[postoperatorcsv]")
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
TEST_CASE("PostOperatorCSV_Restart_TwoExcitation", "[postoperatorcsv]")
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

// NOLINTEND(cppcoreguidelines-avoid-do-while)
