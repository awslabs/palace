// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <iterator>
#include <sstream>
#include <scn/scan.h>
#include <catch2/catch_test_macros.hpp>
#include "fixtures.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/tablecsv.hpp"

using namespace palace;

// Small tests for parsing using scn library that checks assumptions implemented in table.
TEST_CASE("CheckScnCases", "[tablecsv][Serial]")
{
  {
    auto result = scn::scan<double>("-1.00", "{}");
    CHECK(result->value() == -1.0);
  }
  {
    auto result = scn::scan<double>("+1.00", "{}");
    CHECK(result->value() == 1.0);
  }
  {
    auto result = scn::scan<double>("2", "{}");
    CHECK(result->value() == 2.0);
  }
  {
    auto result = scn::scan<double>("+2.00E-03", "{}");
    CHECK(result->value() == 2.00E-03);
  }
  {
    auto result = scn::scan<double>("-2.00E+03", "{}");
    CHECK(result->value() == -2.00E+03);
  }
}

TEST_CASE("TableCSV", "[tablecsv][Serial]")
{
  Table table{};
  table.col_options.float_precision = 9;
  table.reserve(5, 2);

  // Quick defaults.
  CHECK(table.print_col_separator == ",");
  CHECK(table.print_row_separator == "\n");

  REQUIRE(table.n_rows() == 0);
  REQUIRE(table.n_cols() == 0);

  // Add and access columns.
  {
    auto status_1 = table.insert("col_1");
    REQUIRE(status_1);
    auto &col_1 = table["col_1"];
    CHECK(col_1.header_text == "");

    col_1.header_text = "Header Col 1";
    auto &col_1i = table[0];
    CHECK(col_1i.header_text == "Header Col 1");

    auto status_2 = table.insert("col_2", "Header Col 2");
    REQUIRE(status_2);

    auto &col_2 = table["col_2"];
    CHECK(col_2.header_text == "Header Col 2");
    col_2.data.push_back(2.0);

    CHECK(col_1.data.capacity() == 5);
    CHECK(col_2.data.capacity() == 5);
  }

  REQUIRE(table.n_rows() == 1);
  REQUIRE(table.n_cols() == 2);

  // Invalid access.
  CHECK_THROWS(table["col3"]);
  CHECK_THROWS(table[2]);

  // Check reserved: invalidates references.
  {
    table.reserve(6, 3);
    CHECK(table["col_1"].data.capacity() == 6);
    CHECK(table["col_2"].data.capacity() == 6);
  }

  {
    table.insert(Column("col_3", "Header Col 3"));
    auto &col_3 = table["col_3"];
    CHECK(col_3.header_text == "Header Col 3");
    col_3.data.push_back(3.0);
    col_3 << 6.0;
  }

  std::vector<std::size_t> cols_n_row;
  std::transform(table.cbegin(), table.cend(), std::back_inserter(cols_n_row),
                 [](auto &c) { return c.n_rows(); });

  CHECK(cols_n_row == std::vector<std::size_t>{0, 1, 2});

  CHECK(table.n_cols() == 3);
  CHECK(table.n_rows() == 2);

  // clang-format off
  auto table_str1 = std::string(
    "            Header Col 1,            Header Col 2,            Header Col 3\n"
    "                    NULL,        +2.000000000e+00,        +3.000000000e+00\n"
    "                    NULL,                    NULL,        +6.000000000e+00\n"
  );
  // clang-format on
  CHECK(table.format_table() == table_str1);

  auto &col_2 = table["col_2"];
  col_2.min_left_padding = 5;
  col_2.float_precision = 6;

  // clang-format off
  auto table_str2 = std::string(
    "            Header Col 1,      Header Col 2,            Header Col 3\n"
    "                    NULL,     +2.000000e+00,        +3.000000000e+00\n"
    "                    NULL,              NULL,        +6.000000000e+00\n"
  );
  // clang-format on
  CHECK(table.format_table() == table_str2);

  table["col_2"].fmt_sign = " ";

  col_2.min_left_padding = 0;
  col_2.float_precision = 2;

  // clang-format off
  auto table_str3 = std::string(
    "            Header Col 1,Header Col 2,            Header Col 3\n"
    "                    NULL,    2.00e+00,        +3.000000000e+00\n"
    "                    NULL,        NULL,        +6.000000000e+00\n"
  );
  // clang-format on
  CHECK(table.format_table() == table_str3);

  col_2.fmt_sign.reset();
  col_2.min_left_padding.reset();
  col_2.float_precision.reset();
  CHECK(table.format_table() == table_str1);
}

TEST_CASE("TableCSVParsing1_Basic", "[tablecsv][Serial]")
{
  Table table_expected{};
  {
    table_expected.insert("col_1", "Header Col 1");

    table_expected.insert("col_2", "Header Col 2");
    table_expected["col_2"] << 20.0;

    table_expected.insert("col_3", "Header Col 3");
    table_expected["col_3"] << -3.0 << 6.0;
  }

  auto table_str1 = std::string(
      "            Header Col 1,            Header Col 2,            Header Col 3\n"
      "                    NULL,        +2.000000000e+01,        -3.000000000e+00\n"
      "                    NULL,                    NULL,        +6.000000000e+00\n");

  Table table_parse(table_str1);

  CHECK(table_parse.n_cols() == table_expected.n_cols());
  CHECK(table_parse[0].data == table_expected[0].data);
  CHECK(table_parse[1].data == table_expected[1].data);
  CHECK(table_parse[2].data == table_expected[2].data);

  CHECK(table_parse[0].header_text == table_expected[0].header_text);
  CHECK(table_parse[1].header_text == table_expected[1].header_text);
  CHECK(table_parse[2].header_text == table_expected[2].header_text);
}

TEST_CASE("TableCSVParsing2_NonDefaultSeparators", "[tablecsv][Serial]")
{
  using namespace std::literals;

  Table table_expected{};
  {
    table_expected.insert("col_1", "Header Col 1");

    table_expected.insert("col_2", "Header Col 2");
    table_expected["col_2"] << 20.0;

    table_expected.insert("col_3", "Header Col 3");
    table_expected["col_3"] << -3.0 << 6.0;
  }

  auto table_str1 = std::string(
      "            Header Col 1;            Header Col 2;            Header Col 3\r"
      "                    NULL;        +2.000000000e+01;        -3.000000000e+00\r"
      "                    NULL;                    NULL;        +6.000000000e+00\r");

  Table table_parse(table_str1, ";"sv, "\r"sv);

  CHECK(table_parse.n_cols() == table_expected.n_cols());
  CHECK(table_parse[0].data == table_expected[0].data);
  CHECK(table_parse[1].data == table_expected[1].data);
  CHECK(table_parse[2].data == table_expected[2].data);

  CHECK(table_parse[0].header_text == table_expected[0].header_text);
  CHECK(table_parse[1].header_text == table_expected[1].header_text);
  CHECK(table_parse[2].header_text == table_expected[2].header_text);
}

TEST_CASE("TableCSVParsing3_EmptyCells", "[tablecsv][Serial]")
{
  Table table_expected{};
  {
    table_expected.insert("col_1", "Header Col 1");

    table_expected.insert("col_2", "Header Col 2");
    table_expected["col_2"] << 20.0;

    table_expected.insert("col_3", "Header Col 3");
    table_expected["col_3"] << 3.0;
  }

  auto table_str1 = std::string(
      "            Header Col 1,            Header Col 2,            Header Col 3\n"
      "                        ,        2.000000000e+01,        3.000000000e+00\n"
      "                    NULL,                    NULL,         \n");

  Table table_parse(table_str1);

  CHECK(table_parse.n_cols() == table_expected.n_cols());
  CHECK(table_parse[0].data == table_expected[0].data);
  CHECK(table_parse[1].data == table_expected[1].data);
  CHECK(table_parse[2].data == table_expected[2].data);

  CHECK(table_parse[0].header_text == table_expected[0].header_text);
  CHECK(table_parse[1].header_text == table_expected[1].header_text);
  CHECK(table_parse[2].header_text == table_expected[2].header_text);
}

TEST_CASE("TableCSVParsing_TrimSuffix", "[tablecsv][Serial]")
{
  Table table_expected{};
  {
    table_expected.insert("col_1", "Header Col 1");

    table_expected.insert("col_2", "Header Col 2");
    table_expected["col_2"] << 20.0;

    table_expected.insert("col_3", "Header Col 3");
    table_expected["col_3"] << 3.0;
  }

  auto table_str1 = std::string(
      "            Header Col 1   ,            Header Col 2 ,            Header Col 3 \n  "
      "                         ,        2.000000000e+01  ,        3.000000000e+00\n  "
      "                    NULL  ,                    NULL ,         \n  ");

  Table table_parse(table_str1);

  CHECK(table_parse.n_cols() == table_expected.n_cols());
  CHECK(table_parse[0].data == table_expected[0].data);
  CHECK(table_parse[1].data == table_expected[1].data);
  CHECK(table_parse[2].data == table_expected[2].data);

  CHECK(table_parse[0].header_text == table_expected[0].header_text);
  CHECK(table_parse[1].header_text == table_expected[1].header_text);
  CHECK(table_parse[2].header_text == table_expected[2].header_text);
}

TEST_CASE("TableCSV_LoadFromFile", "[tablecsv][Serial]")
{
  // Make these tests serial to avoid duplicate file access.
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  SECTION("Empty File")
  {
    auto no_file = fs::path(PALACE_TEST_DATA_DIR) /
                   "postoperatorcsv_restart/restart1_all/does-not-exists.csv";
    REQUIRE(!fs::exists(no_file));
    TableWithCSVFile table_w(no_file, true);
    CHECK(table_w.table.empty());
  }

  SECTION("Normal File")
  {
    auto test_file =
        fs::path(PALACE_TEST_DATA_DIR) / "postoperatorcsv_restart/restart1_all/port-V.csv";
    REQUIRE(fs::exists(test_file));

    TableWithCSVFile table_w(test_file, true);

    CHECK(table_w.table.n_rows() == 6);
    CHECK(table_w.table.n_cols() == 8);
    CHECK(table_w.table[0].data == std::vector<double>{2, 8, 14, 20, 26, 32});
    CHECK(table_w.table[1].data == std::vector<double>{1, 1, 1, 1, 1, 1});
  }
}

namespace
{

std::string ReadFileToString(const fs::path &path)
{
  std::ifstream f(path, std::ios_base::in);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

void WriteFileFromString(const fs::path &path, const std::string &content)
{
  std::ofstream f(path, std::ios_base::out | std::ios_base::trunc);
  f << content;
}

}  // namespace

// The incremental write must be byte-identical to a whole-file write of the same table,
// and a reload after every step must round-trip the data, since a restart reads the file
// back at every measurement.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_IncrementalMatchesFullWrite",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto incremental_path = temp_dir / "incremental.csv";
  const auto full_path = temp_dir / "full.csv";

  TableWithCSVFile out(incremental_path);
  out.table.col_options.float_precision = 9;
  out.table.insert("idx", "f (GHz)", -1);
  out.table.insert("v_1", "V1 (V)", 0);

  for (int i = 0; i < 5; i++)
  {
    out.table["idx"] << 0.25 * i;
    out.table["v_1"] << 2.0 * i;
    out.WriteTableIncremental();

    TableWithCSVFile reloaded(incremental_path, true);
    REQUIRE(reloaded.table.n_cols() == out.table.n_cols());
    REQUIRE(reloaded.table.n_rows() == out.table.n_rows());
    for (std::size_t j = 0; j < out.table.n_cols(); j++)
    {
      CHECK(reloaded.table[j].data == out.table[j].data);
    }
  }

  TableWithCSVFile full(full_path);
  full.table = out.table;
  full.WriteFullTableTrunc();
  CHECK(ReadFileToString(incremental_path) == ReadFileToString(full_path));
  CHECK(ReadFileToString(incremental_path) == out.table.format_table());
}

// Rows already on disk are never re-rendered, so the bytes written grow with the step count
// instead of its square. The corrupted byte below survives only if the writer appends.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "TableCSV_IncrementalDoesNotRewriteRowsOnDisk", "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "append.csv";

  TableWithCSVFile out(path);
  out.table.col_options.float_precision = 9;
  out.table.insert("idx", "f (GHz)", -1);
  out.table.insert("v_1", "V1 (V)", 0);
  for (int i = 0; i < 3; i++)
  {
    out.table["idx"] << 0.25 * i;
    out.table["v_1"] << 2.0 * i;
  }
  out.WriteTableIncremental();

  // Corrupt the first byte of the first data row in place, same length.
  auto content = ReadFileToString(path);
  const auto header = out.table.format_header();
  REQUIRE(content.size() > header.size());
  content[header.size()] = '#';
  WriteFileFromString(path, content);

  for (int i = 3; i < 5; i++)
  {
    out.table["idx"] << 0.25 * i;
    out.table["v_1"] << 2.0 * i;
  }
  out.WriteTableIncremental();

  CHECK(ReadFileToString(path) ==
        content + out.table.format_row(3) + out.table.format_row(4));
}

// A multi-excitation sweep fills the table one column group per pass, so rows are
// partially filled mid-sweep. The file must carry that fill state (NULLs included), and
// once a later pass completes the rows, the whole-file fallback must revise them.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_MultiExcitationFillStateRoundTrip",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "sweep.csv";

  TableWithCSVFile out(path);
  Table &table = out.table;
  table.col_options.float_precision = 9;
  table.insert("idx", "f (GHz)", -1);
  table.insert("a_1", "A1 (J)", 0);
  table.insert("c_2", "C2 (J)", 1);
  constexpr std::size_t n_rows = 3;

  // First excitation pass fills the index and group 0 columns; group 1 stays empty.
  for (std::size_t i = 0; i < n_rows; i++)
  {
    table["idx"] << 0.25 * i;
    table["a_1"] << 1.0 * i;
  }
  out.WriteTableIncremental();

  {
    TableWithCSVFile reloaded(path, true);
    REQUIRE(reloaded.table.n_rows() == n_rows);
    CHECK(reloaded.table[0].data == table[0].data);
    CHECK(reloaded.table[1].data == table[1].data);
    CHECK(reloaded.table[2].data.empty());
  }

  // Second excitation pass fills the group 1 columns. The earlier rows were partially
  // filled, so only a whole-file write can now put real values in their group 1 cells.
  for (std::size_t i = 0; i < n_rows; i++)
  {
    table["c_2"] << 2.0 * i;
  }
  out.WriteTableIncremental();

  {
    TableWithCSVFile reloaded(path, true);
    REQUIRE(reloaded.table.n_rows() == n_rows);
    CHECK(reloaded.table[2].data == table[2].data);
  }
  CHECK(ReadFileToString(path) == table.format_table());
}

// A new run at the same path must never append to the previous run's rows.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_NewRunOverwritesPreviousRunOutput",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "previous.csv";

  TableWithCSVFile run_a(path);
  run_a.table.col_options.float_precision = 9;
  run_a.table.insert("idx", "f (GHz)", -1);
  run_a.table.insert("v_1", "V1 (V)", 0);
  for (int i = 0; i < 3; i++)
  {
    run_a.table["idx"] << 0.25 * i;
    run_a.table["v_1"] << 2.0 * i;
  }
  run_a.WriteTableIncremental();
  REQUIRE(fs::exists(path));

  TableWithCSVFile run_b(path);
  run_b.table.col_options.float_precision = 9;
  run_b.table.insert("idx", "f (GHz)", -1);
  run_b.table.insert("v_1", "V1 (V)", 0);
  for (int i = 0; i < 2; i++)
  {
    run_b.table["idx"] << 0.5 * i;
    run_b.table["v_1"] << 3.0 * i;
  }
  run_b.WriteTableIncremental();

  CHECK(ReadFileToString(path) == run_b.table.format_table());
}

// If the file disappears under a running table (deleted externally, or a symlink swap),
// the next write must rebuild it from memory instead of appending into the void.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_WriteRecoversWhenFileVanishes",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "vanishes.csv";

  TableWithCSVFile out(path);
  out.table.col_options.float_precision = 9;
  out.table.insert("idx", "f (GHz)", -1);
  out.table.insert("v_1", "V1 (V)", 0);
  for (int i = 0; i < 2; i++)
  {
    out.table["idx"] << 0.25 * i;
    out.table["v_1"] << 2.0 * i;
  }
  out.WriteTableIncremental();

  fs::remove(path);
  for (int i = 2; i < 4; i++)
  {
    out.table["idx"] << 0.25 * i;
    out.table["v_1"] << 2.0 * i;
  }
  out.WriteTableIncremental();

  CHECK(ReadFileToString(path) == out.table.format_table());
}

// A table with columns but no rows still writes the header, matching a whole-file write,
// and later rows land in the file as usual.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_EmptyTableWritesHeaderOnly",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "empty.csv";

  TableWithCSVFile out(path);
  out.table.col_options.float_precision = 9;
  out.table.insert("idx", "f (GHz)", -1);
  out.table.insert("v_1", "V1 (V)", 0);
  out.WriteTableIncremental();
  CHECK(ReadFileToString(path) == out.table.format_header());

  out.table["idx"] << 0.25;
  out.table["v_1"] << 2.0;
  out.WriteTableIncremental();
  CHECK(ReadFileToString(path) == out.table.format_table());
}
