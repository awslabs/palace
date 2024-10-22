// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_TABLECSV_HPP
#define PALACE_UTILS_TABLECSV_HPP

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace palace::table
{

struct ColumnOptions
{
  size_t min_left_padding = 8;
  size_t float_precision = 9;
  double empty_val = std::nan("1");
  std::string fmt_sign = {"+"};
};

class Column
{
  friend class Table;

  // View to default options in table class, will be set when columns are added to Table
  ColumnOptions *defaults = nullptr;

  // Printing functions are private as they should be called by the Table class only:

  [[nodiscard]] size_t col_width() const
  {
    size_t pad = min_left_padding.value_or(defaults->min_left_padding);
    size_t prec = float_precision.value_or(defaults->float_precision);

    // Normal float in our exponent format needs float_precision + 7 ("+" , leading digit,
    // ".", "e", "+", +2 exponent. Sometimes exponent maybe +3 if very small or large; see
    // std::numeric_limits<double>::max_exponent. So pick +8 to avoid this corner case.
    return std::max(pad + prec + 8, header_text.size());
  }

  [[nodiscard]] auto format_header(const std::optional<size_t> &width = {}) const
  {
    auto w = width.value_or(col_width());
    return fmt::format("{:>{w}s}", header_text, fmt::arg("w", w));
  }

  [[nodiscard]] auto format_row(size_t i, const std::optional<size_t> &width = {}) const
  {
    auto val =
        ((i >= 0) && (i < data.size()) ? data[i] : empty_val.value_or(defaults->empty_val));

    auto sign = fmt_sign.value_or(defaults->fmt_sign);
    auto width_ = width.value_or(col_width());
    auto prec = float_precision.value_or(defaults->float_precision);
    auto fmt_str = fmt::format("{{:>{sign:s}{width}.{prec}e}}", fmt::arg("sign", sign),
                               fmt::arg("width", width_), fmt::arg("prec", prec));
    return fmt::format(fmt::runtime(fmt_str), val);
  }

  // ----

  // Map-like index, to interface via Table class
  std::string name = "";

public:
  Column(std::string name, std::string header_text = "")
    : name(std::move(name)), header_text(std::move(header_text))
  {
  }

  std::vector<double> data = {};
  std::string header_text = "";

  std::optional<size_t> min_left_padding = {};
  std::optional<size_t> float_precision = {};
  std::optional<double> empty_val = {};
  std::optional<std::string> fmt_sign = {};

  [[nodiscard]] size_t n_rows() const { return data.size(); }
};

class Table
{
  // Column-wise mini-table for storing data and and printing to csv file for doubles.
  // Future: allow int and other output, allow non-owning memeory via span
  std::vector<Column> cols = {};

  // Cache value to reserve vector space by default
  size_t reserve_n_rows = 0;

public:
  // Default column options; can be overwritten column-wise
  ColumnOptions col_options = {};

  // Global printing options
  std::string print_col_separator = ",";
  std::string print_row_separator = "\n";

  // Table properties

  [[nodiscard]] size_t n_cols() const { return cols.size(); }
  [[nodiscard]] size_t n_rows() const
  {
    if (n_cols() == 0)
    {
      return 0;
    }
    auto max_col =
        std::max_element(cols.begin(), cols.end(), [](const auto &a, const auto &b)
                         { return a.n_rows() < b.n_rows(); });
    return max_col->n_rows();
  }

  void reserve(size_t n_rows, size_t n_cols)
  {
    reserve_n_rows = n_rows;
    cols.reserve(n_cols);
    for (auto &col : cols)
    {
      col.data.reserve(n_rows);
    }
  }

  // Insert columns: map like interface

  bool insert_column(std::string column_name, std::string column_header = "")
  {
    auto it = std::find_if(cols.begin(), cols.end(),
                           [&column_name](auto &c) { return c.name == column_name; });
    if (it != cols.end())
    {
      return false;
    }
    auto &col = cols.emplace_back(std::move(column_name), std::move(column_header));
    col.defaults = &col_options;
    if (reserve_n_rows > 0)
    {
      col.data.reserve(reserve_n_rows);
    }
    return true;
  }
  bool insert_column(Column &&column)
  {
    auto it = std::find_if(cols.begin(), cols.end(),
                           [&column](auto &c) { return c.name == column.name; });
    if (it != cols.end())
    {
      return false;
    }
    auto &col = cols.emplace_back(std::move(column));
    col.defaults = &col_options;
    if (reserve_n_rows > 0)
    {
      col.data.reserve(reserve_n_rows);
    }
    return true;
  }

  // Access columns via vector position or column name

  Column &at(size_t idx) { return cols.at(idx); }
  const Column &at(size_t idx) const { return cols.at(idx); }

  Column &at(std::string_view name)
  {
    auto it =
        std::find_if(cols.begin(), cols.end(), [&name](auto &c) { return c.name == name; });
    if (it == cols.end())
    {
      throw std::out_of_range(fmt::format("Column {} not found in table", name).c_str());
    }
    return *it;
  }
  const Column &at(std::string_view name) const { return at(name); }

  auto begin() { return cols.begin(); }
  auto end() { return cols.end(); }
  auto cbegin() const { return cols.begin(); }
  auto cend() const { return cols.end(); }

  // Formatting and Printing Options
  // TODO: Improve all the functions below with ranges in C++20

private:
  template <typename T>
  void append_header(T &buf) const
  {
    auto to = [&buf](auto f, auto &&...a)
    { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

    for (int i = 0; i < n_cols(); i++)
    {
      if (i > 0)
      {
        to("{:s}", print_col_separator);
      }
      to("{:s}", cols[i].format_header());
    }
    to("{:s}", print_row_separator);
  }

  template <typename T>
  void append_row(T &buf, size_t row_j) const
  {
    auto to = [&buf](auto f, auto &&...a)
    { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

    for (int i = 0; i < n_cols(); i++)
    {
      if (i > 0)
      {
        to("{:s}", print_col_separator);
      }
      to("{:s}", cols[i].format_row(row_j));
    }
    to("{:s}", print_row_separator);
  }

public:
  [[nodiscard]] std::string format_table()
  {
    fmt::memory_buffer buf{};
    append_header(buf);
    for (int j = 0; j < n_rows(); j++)
    {
      append_row(buf, j);
    }
    return {buf.data(), buf.size()};
  }
};

}  // namespace palace::table

#endif  // PALACE_UTILS_TABLECSV_HPP
