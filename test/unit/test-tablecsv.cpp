// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <iterator>
#include <limits>
#include <fmt/format.h>
#include <scn/scan.h>
#include <catch2/catch_test_macros.hpp>
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
