#include <iterator>
#include <fmt/format.h>
#include <scn/scan.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/tablecsv.hpp"

using namespace palace;

// Small tests for parsing using scn library that checks assumptions implemented in table.
TEST_CASE("CheckScnCases", "[tablecsv]")
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

TEST_CASE("TableCSV", "[tablecsv]")
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

  std::vector<size_t> cols_n_row;
  std::transform(table.cbegin(), table.cend(), std::back_inserter(cols_n_row),
                 [](auto &c) { return c.n_rows(); });

  CHECK(cols_n_row == std::vector<size_t>{0, 1, 2});

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

  //   REQUIRE_NOTHROW(boundary_ex_bool.SetUp(*config.find("boundaries_1_pass")));
}

TEST_CASE("TableCSVParsing1_Basic", "[tablecsv]")
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

TEST_CASE("TableCSVParsing2_NonDefaultSeparators", "[tablecsv]")
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

TEST_CASE("TableCSVParsing3_EmptyCells", "[tablecsv]")
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

TEST_CASE("TableCSVParsing_TrimSuffix", "[tablecsv]")
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

TEST_CASE("TableCSV_LoadFromFile", "[tablecsv]")
{
  // Make these tests serial to avoid duplicate file access.
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  SECTION("Empty File")
  {
    auto no_file = fs::path(PALACE_TEST_DIR) /
                   "postoperatorcsv_restart/restart1_all/does-not-exists.csv";
    REQUIRE(!fs::exists(no_file));
    TableWithCSVFile table_w(no_file, true);
    CHECK(table_w.table.empty());
  }

  SECTION("Normal File")
  {
    auto test_file =
        fs::path(PALACE_TEST_DIR) / "postoperatorcsv_restart/restart1_all/port-V.csv";
    REQUIRE(fs::exists(test_file));

    TableWithCSVFile table_w(test_file, true);

    CHECK(table_w.table.n_rows() == 6);
    CHECK(table_w.table.n_cols() == 8);
    CHECK(table_w.table[0].data == std::vector<double>{2, 8, 14, 20, 26, 32});
    CHECK(table_w.table[1].data == std::vector<double>{1, 1, 1, 1, 1, 1});
  }
}
