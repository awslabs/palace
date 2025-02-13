#include <iterator>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "utils/tablecsv.hpp"

using namespace palace;

TEST_CASE("TableCSV", "[tablecsv]")
{
  Table table{};
  table.reserve(5, 2);

  // Quick defaults
  CHECK(table.print_col_separator == ",");
  CHECK(table.print_row_separator == "\n");

  REQUIRE(table.n_rows() == 0);
  REQUIRE(table.n_cols() == 0);

  // Add and access columns
  {
    auto status_1 = table.insert_column("col_1");
    REQUIRE(status_1);
    auto &col_1 = table["col_1"];
    CHECK(col_1.header_text == "");

    col_1.header_text = "Header Col 1";
    auto &col_1i = table[0];
    CHECK(col_1i.header_text == "Header Col 1");

    auto status_2 = table.insert_column("col_2", "Header Col 2");
    REQUIRE(status_2);

    auto &col_2 = table["col_2"];
    CHECK(col_2.header_text == "Header Col 2");
    col_2.data.push_back(2.0);

    CHECK(col_1.data.capacity() == 5);
    CHECK(col_2.data.capacity() == 5);
  }

  REQUIRE(table.n_rows() == 1);
  REQUIRE(table.n_cols() == 2);

  // Invalid access
  CHECK_THROWS(table["col3"]);
  CHECK_THROWS(table[2]);

  // Check reserved: invalidates references
  {
    table.reserve(6, 3);
    CHECK(table["col_1"].data.capacity() == 6);
    CHECK(table["col_2"].data.capacity() == 6);
  }

  {
    table.insert_column(Column("col_3", "Header Col 3"));
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