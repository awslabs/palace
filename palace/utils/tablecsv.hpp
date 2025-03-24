// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TABLECSV_HPP
#define PALACE_UTILS_TABLECSV_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace palace
{

struct ColumnOptions
{
  std::size_t min_left_padding = 8;
  std::size_t float_precision = 9;
  std::string empty_cell_val = {"NULL"};
  std::string fmt_sign = {"+"};
};

class Column
{
  friend class Table;

  // View to default options in table class, will be set when columns are added to Table
  ColumnOptions *defaults = nullptr;

  // Map-like index, to interface via Table class
  std::string name;

public:
  [[nodiscard]] std::size_t col_width() const;

  [[nodiscard]] auto format_header(const std::optional<std::size_t> &width = {}) const;

  [[nodiscard]] auto format_row(std::size_t i, const std::optional<std::size_t> &width = {}) const;

  Column(std::string name, std::string header_text = "",
         std::optional<std::size_t> min_left_padding = {},
         std::optional<std::size_t> float_precision = {},
         std::optional<std::string> empty_cell_val = {},
         std::optional<std::string> fmt_sign = {});

  std::vector<double> data;
  std::string header_text;

  std::optional<std::size_t> min_left_padding;
  std::optional<std::size_t> float_precision;
  std::optional<std::string> empty_cell_val;
  std::optional<std::string> fmt_sign;

  // Quick-fix since leading column of eig is int stored as double (rather then implementing
  // a templated & type erased solution).
  bool print_as_int = false;

  [[nodiscard]] inline std::size_t n_rows() const { return data.size(); }

  // Convenience operator at higher level
  inline auto operator<<(double val)
  {
    data.emplace_back(val);
    return *this;
  }
};

class Table
{
  // Column-wise mini-table for storing data and and printing to csv file for doubles.
  // Future: allow int and other output, allow non-owning memory via span
  std::vector<Column> cols;

  // Cache value to reserve vector space by default
  std::size_t reserve_n_rows = 0;

public:
  // Default column options; can be overwritten column-wise
  ColumnOptions col_options = {};

  // Global printing options
  std::string print_col_separator = ",";
  std::string print_row_separator = "\n";

  // Table properties

  [[nodiscard]] std::size_t n_cols() const { return cols.size(); }
  [[nodiscard]] std::size_t n_rows() const;

  void reserve(std::size_t n_rows, std::size_t n_cols);

  // Insert columns: map like interface
  bool insert(Column &&column);
  template <typename... Args>
  bool insert(Args &&...args)
  {
    return insert(Column(std::forward<Args>(args)...));
  }

  // Access columns via vector position or column name

  inline Column &operator[](std::size_t idx) { return cols.at(idx); }
  inline const Column &operator[](std::size_t idx) const { return (*this)[idx]; }

  Column &operator[](std::string_view name);
  const Column &operator[](std::string_view name) const { return (*this)[name]; }

  inline auto begin() { return cols.begin(); }
  inline auto end() { return cols.end(); }
  inline auto cbegin() const { return cols.begin(); }
  inline auto cend() const { return cols.end(); }

  // Formatting and Printing Options
  template <typename T>
  void append_header(T &buf) const;

  template <typename T>
  void append_row(T &buf, std::size_t row_j) const;

  [[nodiscard]] std::string format_header() const;

  [[nodiscard]] std::string format_row(std::size_t j) const;

  [[nodiscard]] std::string format_table() const;
};

// Wrapper for storing Table to csv file wish row wise updates

class TableWithCSVFile
{
  std::string csv_file_fullpath_;

  // Index to keep track of which row we are currently at the beginning of / printing. Row
  // [-1, 0) is the header, row [0, 1) the first numeric row, etc.
  unsigned long file_append_cursor = -1;

public:
  Table table = {};

  TableWithCSVFile() = default;
  explicit TableWithCSVFile(std::string csv_file_fullpath);

  void WriteFullTableTrunc();

  void AppendHeader();

  void AppendRow();

  [[nodiscard]] auto GetAppendRowCursor() const { return file_append_cursor; }
};

}  // namespace palace

#endif  // PALACE_UTILS_TABLECSV_HPP
