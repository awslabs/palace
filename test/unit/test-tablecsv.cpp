// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <fmt/format.h>
#include <scn/scan.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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
  std::ranges::transform(table, std::back_inserter(cols_n_row),
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

TEST_CASE("TableCSV nonfinite integer formatting", "[tablecsv][Serial]")
{
  Table table;
  table.insert(Column("status", "status", 0, 0, 0, ""));
  table["status"].print_as_int = true;
  table["status"] << std::numeric_limits<double>::quiet_NaN();
  CHECK(table.format_table() == "status\nnan\n");

  Table floating;
  floating.insert("quality", "quality");
  floating["quality"] << std::numeric_limits<double>::infinity();
  floating["quality"] << -std::numeric_limits<double>::infinity();
  const auto formatted = floating.format_table();
  CHECK(formatted.find("+inf") != std::string::npos);
  CHECK(formatted.find("-inf") != std::string::npos);
  CHECK(formatted.find("nan") == std::string::npos);
}

// Decisions 1829 / 1835 (OI-CSV-SIGNED-NAN-DIRECT-WRITERS): a number of a directly written
// CSV table prints as "{:+.12e}" (finite and +inf / -inf: "+inf" / "-inf"), a NaN of either
// sign as the bare token "nan" (fmt's "+" flag would print "+nan", which pandas reads as
// text); the token is the trimmed cell of a Table column of the same value.
TEST_CASE("TableCSV CsvNumber direct-writer cells", "[tablecsv][Serial]")
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(fmt::format("{:+.12e}", nan) == "+nan");  // the form CsvNumber replaces
  CHECK(CsvNumber(nan) == "nan");
  CHECK(CsvNumber(-nan) == "nan");
  CHECK(CsvNumber(inf) == "+inf");
  CHECK(CsvNumber(-inf) == "-inf");
  CHECK(CsvNumber(0.5) == "+5.000000000000e-01");
  CHECK(CsvNumber(-1.25e-18) == "-1.250000000000e-18");
  CHECK(CsvNumber(0.0) == "+0.000000000000e+00");
  Table table;
  table.insert("value", "value");
  for (const double value : {0.5, -1.25e-18, nan, inf, -inf})
  {
    table["value"] << value;
  }
  for (std::size_t i = 0; i < 5; i++)
  {
    const std::string cell = table.format_row(i);
    const auto begin = cell.find_first_not_of(' ');
    CHECK(cell.substr(begin, cell.size() - 1 - begin) == CsvNumber(table["value"].data[i]));
  }
  for (const double value : {0.5, -1.25e-18, nan, inf, -inf})
  {
    auto result = scn::scan<double>(CsvNumber(value), "{}");
    REQUIRE(result);
    CHECK((std::isnan(value) ? std::isnan(result->value()) : result->value() == value));
  }
}

// Lane L4 (decision 1829 C1): a non-finite numeric cell prints as a BARE token (no padding)
// in float and integer columns, beside text columns; every finite cell, every text cell,
// the header and the padded empty cell ("NULL", out of scope) are printed as before; the
// table reads back with the same values.
TEST_CASE("TableCSV bare non-finite cells", "[tablecsv][Serial]")
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  Table table;
  table.insert("value", "value");
  table.insert(Column::Integer("index", "index"));
  table.insert(Column::Text("status", "coverage status"));
  table["value"] << 0.5 << nan << inf << -inf;
  table["index"] << 3 << nan << inf << -inf;
  table["status"] << std::string("OK") << std::string("nan") << std::string("PARTIAL")
                  << std::string("NO-ENERGY");
  auto Pad = [](const std::string &cell, std::size_t width)
  { return std::string(width - cell.size(), ' ') + cell; };
  const auto index_width = table["index"].col_width(table.col_options);
  CHECK(table.format_header() == Pad("value", 27) + "," + Pad("index", index_width) + "," +
                                     Pad("coverage status", 27) + "\n");
  CHECK(table.format_row(0) == Pad("+5.000000000000e-01", 27) + "," +
                                   Pad("3", index_width) + "," + Pad("OK", 27) + "\n");
  CHECK(table.format_row(1) == "nan,nan," + Pad("nan", 27) + "\n");
  CHECK(table.format_row(2) == "+inf,inf," + Pad("PARTIAL", 27) + "\n");
  CHECK(table.format_row(3) == "-inf,-inf," + Pad("NO-ENERGY", 27) + "\n");

  // A shorter column pads its missing cell with the empty-cell value, as before.
  table["value"] << 1.0;
  CHECK(table.format_row(4) == Pad("+1.000000000000e+00", 27) + "," +
                                   Pad("NULL", index_width) + "," + Pad("NULL", 27) + "\n");

  Table parsed(table.format_table(), std::nullopt, std::nullopt, /*text_columns=*/true);
  REQUIRE(parsed.n_cols() == 3);
  REQUIRE(parsed[0].data.size() == 5);
  CHECK(parsed[0].data[0] == 0.5);
  CHECK(std::isnan(parsed[0].data[1]));
  CHECK(parsed[0].data[2] == inf);
  CHECK(parsed[0].data[3] == -inf);
  CHECK(parsed[0].data[4] == 1.0);
  REQUIRE(parsed[1].data.size() == 4);
  CHECK(parsed[1].data[0] == 3.0);
  CHECK(std::isnan(parsed[1].data[1]));
  CHECK(parsed[1].data[2] == inf);
  CHECK(parsed[1].data[3] == -inf);
  CHECK(parsed[2].text == std::vector<std::string>{"OK", "nan", "PARTIAL", "NO-ENERGY"});
}

// Lane L1 (decision 1634 (3)): a text column (the coverage / regime status columns) prints
// its strings right-aligned at the numeric column width; a numeric column beside it is
// unchanged.
TEST_CASE("TableCSV text column", "[tablecsv][Serial]")
{
  Table table;
  table.insert("value", "value");
  table.insert(Column::Text("status", "coverage status"));
  table["value"] << 0.5;
  table["status"] << std::string("NOT-ACHIEVABLE");
  table["value"] << 1.0;
  table["status"] << std::string("OK");
  CHECK(table.n_rows() == 2);
  // Every column is 27 wide (8 + 12 + 7, the header shorter): right-aligned cells.
  auto Pad = [](const std::string &cell)
  { return std::string(27 - cell.size(), ' ') + cell; };
  CHECK(table.format_table() == Pad("value") + "," + Pad("coverage status") + "\n" +
                                    Pad("+5.000000000000e-01") + "," +
                                    Pad("NOT-ACHIEVABLE") + "\n" +
                                    Pad("+1.000000000000e+00") + "," + Pad("OK") + "\n");
}

// Lane L1: a printed table with text columns reads back (Table(std::string_view) with
// text_columns): a column whose first cell is not a number is a text column; a text cell
// below a number is refused. Without text_columns (every product reader: review m-2) a
// non-numeric cell is refused as before.
TEST_CASE("TableCSV text column round-trips", "[tablecsv][Serial]")
{
  Table table;
  table.insert("value", "value");
  table.insert(Column::Text("type", "coverage Type"));
  table.insert(Column::Text("status", "coverage status"));
  const std::vector<std::string> types = {"SA", "MS", "MA", "Default"};
  const std::vector<std::string> statuses = {"IN", "NO-RANGE", "NOT-TARGET", "NO-ENERGY"};
  for (std::size_t i = 0; i < types.size(); i++)
  {
    table["value"] << 0.25 * static_cast<double>(i);
    table["type"] << types[i];
    table["status"] << statuses[i];
  }
  Table parsed(table.format_table(), std::nullopt, std::nullopt, /*text_columns=*/true);
  REQUIRE(parsed.n_cols() == 3);
  CHECK(!parsed[0].print_as_text);
  CHECK(parsed[1].print_as_text);
  CHECK(parsed[2].print_as_text);
  CHECK(parsed[1].text == types);
  CHECK(parsed[2].text == statuses);
  CHECK(parsed[0].data == std::vector<double>{0.0, 0.25, 0.5, 0.75});
  CHECK(parsed.format_table() == table.format_table());
  CHECK_THROWS(Table(table.format_table()));
  CHECK_THROWS(Table("value\n1.0\nSA\n", std::nullopt, std::nullopt, true));
}

// OI-PATCH-ENERGY-CSV-INDEX-PRECISION (fix round P2, decision 1265): the former index idiom
// Column(name, header, 0, 0, 2, "") printed an integer through the 3-significant-digit
// float form, so an index beyond 999 collapsed ("1.24e+02" for 124 kept 3 digits; 1234 and
// 1235 both read 1.23e+03 on a join). Column::Integer prints the integer exactly and the
// parsed table reads it back exactly.
TEST_CASE("TableCSV integer index column round-trips exactly", "[tablecsv][Serial]")
{
  Table legacy;
  legacy.insert(Column("patch", "patch", 0, 0, 2, ""));
  legacy["patch"] << 1234 << 1235;
  CHECK(legacy.format_table() == "    patch\n 1.23e+03\n 1.24e+03\n");  // the collapse

  Table table;
  table.insert(Column::Integer("patch", "patch"));
  table.insert(Column("energy", "energy (J)"));
  for (const int index : {1, 999, 1000, 1234, 1235, 123456})
  {
    table["patch"] << index;
    table["energy"] << 1.5e-18 * index;
  }
  const std::string text = table.format_table();
  CHECK(text.rfind("patch,", 0) == 0);
  CHECK(text.find("1234,") != std::string::npos);
  CHECK(text.find("1235,") != std::string::npos);
  CHECK(text.find("123456,") != std::string::npos);
  CHECK(text.find("e+03") == std::string::npos);
  Table parsed(text);
  REQUIRE(parsed[0].n_rows() == 6);
  CHECK(parsed[0].data[3] == 1234.0);
  CHECK(parsed[0].data[4] == 1235.0);
  CHECK(parsed[0].data[5] == 123456.0);
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

// The restart path reloads a table from disk and then restores the real column names from a
// reference table. The name index still holds the loader's placeholder names until it is
// rebuilt.
TEST_CASE("TableCSV_RebuildNameIndexAfterRename", "[tablecsv][Serial]")
{
  Table reference{};
  reference.insert("idx", "f (GHz)", -1);
  reference.insert("Ee_1", "E_elec1 (J)", 0);
  reference["idx"] << 1.0;
  reference["Ee_1"] << 2.0;

  const std::string text = reference.format_table();
  Table loaded{text};
  REQUIRE(loaded.n_cols() == reference.n_cols());
  for (std::size_t i = 0; i < loaded.n_cols(); i++)
  {
    loaded[i].name = reference[i].name;
  }

  CHECK(!loaded.has("idx"));
  CHECK_THROWS(loaded["idx"]);

  loaded.RebuildNameIndex();

  CHECK(loaded.has("idx"));
  CHECK(loaded.has("Ee_1"));
  CHECK(!loaded.has("col_0"));
  CHECK(loaded["idx"].data.at(0) == 1.0);
  CHECK(loaded["Ee_1"].data.at(0) == 2.0);

  // insert must still notice a duplicate of the renamed column.
  CHECK(!loaded.insert("idx", "f (GHz)"));
  CHECK(loaded.n_cols() == reference.n_cols());

  // Renaming two columns to the same name breaks the one-name-one-column invariant.
  loaded[1].name = "idx";
  CHECK_THROWS(loaded.RebuildNameIndex());
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

// If the file disappears under a running table (for example, deleted externally), the next
// write must rebuild it from memory instead of appending into the void.
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

// A restarted run reloads the table from disk and keeps appending: the reloaded rows count
// as already written, so they are not rendered again. The corrupted byte below survives
// only if the writer appends.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_IncrementalAppendsAfterReload",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "reload.csv";
  {
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
  }

  auto content = ReadFileToString(path);
  TableWithCSVFile reloaded(path, true);
  REQUIRE(reloaded.table.n_rows() == 3);
  reloaded.table.col_options.float_precision = 9;
  const auto header = reloaded.table.format_header();
  REQUIRE(content.size() > header.size());
  content[header.size()] = '#';
  WriteFileFromString(path, content);

  for (int i = 3; i < 5; i++)
  {
    reloaded.table[0] << 0.25 * i;
    reloaded.table[1] << 2.0 * i;
  }
  reloaded.WriteTableIncremental();
  CHECK(ReadFileToString(path) ==
        content + reloaded.table.format_row(3) + reloaded.table.format_row(4));
}

// A reloaded file cut off before its final row separator must be rewritten whole on the
// next write, not appended to, which would join two rows on one line.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "TableCSV_ReloadWithoutTrailingSeparator",
                 "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const auto path = temp_dir / "truncated.csv";
  {
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
  }
  auto content = ReadFileToString(path);
  REQUIRE(!content.empty());
  REQUIRE(content.back() == '\n');
  content.pop_back();
  WriteFileFromString(path, content);

  TableWithCSVFile reloaded(path, true);
  REQUIRE(reloaded.table.n_rows() == 2);
  reloaded.table.col_options.float_precision = 9;
  reloaded.table[0] << 0.5;
  reloaded.table[1] << 4.0;
  reloaded.WriteTableIncremental();
  CHECK(ReadFileToString(path) == reloaded.table.format_table());
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

// A table object uses ONE incremental writer (SHARED-DECISIONS SD-602): after
// AppendRowsAndClear has started, WriteFullTableTrunc / WriteTableIncremental are refused,
// and AppendRowsAndClear is refused once another writer has put rows on disk.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "TableCSV incremental writers are exclusive per table", "[tablecsv][Serial]")
{
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  auto AddRow = [](TableWithCSVFile &out, double f)
  {
    out.table["idx"] << f;
    out.table["v_1"] << 2.0 * f;
  };
  auto MakeColumns = [](TableWithCSVFile &out)
  {
    out.table.insert("idx", "f (GHz)", -1);
    out.table.insert("v_1", "V1 (V)", 0);
  };
  // The refusals of tablecsv.cpp (SD-602), matched so that an unrelated exception fails.
  using Catch::Matchers::ContainsSubstring;
  const std::string append_only = "is written by AppendRowsAndClear only!";
  const std::string full_or_incremental =
      "is written by WriteFullTableTrunc / WriteTableIncremental!";

  SECTION("AppendRowsAndClear first")
  {
    TableWithCSVFile out(temp_dir / "append-first.csv");
    MakeColumns(out);
    AddRow(out, 0.25);
    REQUIRE_NOTHROW(out.AppendRowsAndClear());
    AddRow(out, 0.5);
    CHECK_THROWS_WITH(out.WriteTableIncremental(), ContainsSubstring(append_only));
    CHECK_THROWS_WITH(out.WriteFullTableTrunc(), ContainsSubstring(append_only));
    REQUIRE_NOTHROW(out.AppendRowsAndClear());
  }

  SECTION("WriteTableIncremental first")
  {
    TableWithCSVFile out(temp_dir / "incremental-first.csv");
    MakeColumns(out);
    AddRow(out, 0.25);
    REQUIRE_NOTHROW(out.WriteTableIncremental());
    AddRow(out, 0.5);
    CHECK_THROWS_WITH(out.AppendRowsAndClear(), ContainsSubstring(full_or_incremental));
    REQUIRE_NOTHROW(out.WriteTableIncremental());
  }

  SECTION("WriteFullTableTrunc first")
  {
    TableWithCSVFile out(temp_dir / "full-first.csv");
    MakeColumns(out);
    AddRow(out, 0.25);
    REQUIRE_NOTHROW(out.WriteFullTableTrunc());
    CHECK_THROWS_WITH(out.AppendRowsAndClear(), ContainsSubstring(full_or_incremental));
  }
}
