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
  // Options that individual cols can overwrite.
  std::size_t min_left_padding = 8;
  std::size_t float_precision = 12;
  std::string fmt_sign = {"+"};

  // Common options
  std::string empty_cell_val = {"NULL"};
};

class Column
{
  friend class Table;

public:
  // Map-like index, to interface via Table class.
  std::string name;

  [[nodiscard]] std::size_t col_width(const ColumnOptions &defaults = {}) const;

  [[nodiscard]] auto format_header(const ColumnOptions &defaults = {},
                                   const std::optional<std::size_t> &width = {}) const;

  [[nodiscard]] auto format_row(std::size_t i, const ColumnOptions &defaults = {},
                                const std::optional<std::size_t> &width = {}) const;

  Column(std::string name, std::string header_text = "", long column_group_idx = 0,
         std::optional<std::size_t> min_left_padding = {},
         std::optional<std::size_t> float_precision = {},
         std::optional<std::string> fmt_sign = {});

  // Actual Data in Column.
  std::vector<double> data;

  // Pretty text to print in file for column header.
  std::string header_text;

  // Index to group column into blocks, to verify common cursor.
  // We assume that column groups are contiguous in table.
  long column_group_idx = 0;

  // Column-wise printing options that overwrite default.
  std::optional<std::size_t> min_left_padding;
  std::optional<std::size_t> float_precision;
  std::optional<std::string> fmt_sign;

  // Quick-fix since leading column of eig is int stored as double (rather then implementing
  // a templated & type erased solution).
  bool print_as_int = false;

  [[nodiscard]] inline std::size_t n_rows() const { return data.size(); }

  // Convenience operator at higher level.
  inline auto &operator<<(double val)
  {
    data.emplace_back(val);
    return *this;
  }
};

class Table
{
  // Column-wise mini-table for storing data and and printing to csv file for doubles.
  // Future: allow int and other output, allow non-owning memory via span.
  std::vector<Column> cols;

  // Cache value to reserve vector space by default.
  std::size_t reserve_n_rows = 0;

public:
  Table() = default;
  Table(std::string_view table_str,
        std::optional<std::string_view> print_col_separator_ = std::nullopt,
        std::optional<std::string_view> print_row_separator_ = std::nullopt);

  // Default column options; can be overwritten column-wise.
  ColumnOptions col_options = {};

  // Global printing options.
  std::string_view print_col_separator{",", 1};
  std::string_view print_row_separator{"\n", 1};

  // Table properties.

  [[nodiscard]] bool empty() const { return cols.empty(); }
  [[nodiscard]] std::size_t n_cols() const { return cols.size(); }
  [[nodiscard]] std::size_t n_rows() const;

  void reserve(std::size_t n_rows, std::size_t n_cols);

  // Insert columns: map like interface.
  bool insert(Column &&column);
  template <typename... Args>
  bool insert(Args &&...args)
  {
    return insert(Column(std::forward<Args>(args)...));
  }

  // Access columns via vector position or column name.
  inline Column &operator[](std::size_t idx) { return cols.at(idx); }
  Column &operator[](std::string_view name);

  inline auto begin() { return cols.begin(); }
  inline auto end() { return cols.end(); }
  inline auto cbegin() const { return cols.begin(); }
  inline auto cend() const { return cols.end(); }

  // Formatting and Printing Options.
  template <typename T>
  void append_header(T &buf) const;

  template <typename T>
  void append_row(T &buf, std::size_t row_j) const;

  [[nodiscard]] std::string format_header() const;

  [[nodiscard]] std::string format_row(std::size_t j) const;

  [[nodiscard]] std::string format_table() const;
};

// Wrapper for storing Table to csv file.

class TableWithCSVFile
{
  std::string csv_file_fullpath_;

public:
  Table table = {};

  TableWithCSVFile() = default;
  explicit TableWithCSVFile(std::string csv_file_fullpath, bool load_existing_file = false);

  std::string_view get_csv_filepath() const { return {csv_file_fullpath_}; }

  void WriteFullTableTrunc();
};

}  // namespace palace

#endif  // PALACE_UTILS_TABLECSV_HPP
