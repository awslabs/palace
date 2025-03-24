// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "tablecsv.hpp"

#include <algorithm>
#include <iterator>
#include <utility>
#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

namespace palace
{

[[nodiscard]] size_t Column::col_width() const
{
  // Quickfix to specify full column width in integer case to match current formatting
  if (print_as_int)
  {
    return std::max(min_left_padding.value_or(defaults->min_left_padding),
                    header_text.size());
  }
  size_t pad = min_left_padding.value_or(defaults->min_left_padding);
  size_t prec = float_precision.value_or(defaults->float_precision);

  // Normal float in our exponent format needs float_precision + 7 ("+" , leading digit,
  // ".", "e", "+", +2 exponent. Sometimes exponent maybe +3 if very small or large; see
  // std::numeric_limits<double>::max_exponent. We pick +7 for consistency, but
  // min_left_padding should be at least 1, which is not currently enforced.
  return std::max(pad + prec + 7, header_text.size());
}

[[nodiscard]] auto Column::format_header(const std::optional<size_t> &width) const
{
  auto w = width.value_or(col_width());
  return fmt::format("{0:>{1}s}", header_text, w);
}

[[nodiscard]] auto Column::format_row(size_t i, const std::optional<size_t> &width) const
{
  auto width_ = width.value_or(col_width());
  // If data available format double
  if ((i >= 0) && (i < data.size()))
  {
    auto val = data[i];
    if (print_as_int)
    {  // Quick-fix to force int printing
      auto fmt_str = fmt::format("{{:>{width}d}}", fmt::arg("width", width_));
      return fmt::format(fmt::runtime(fmt_str), int(val));
    }
    else
    {
      auto sign = fmt_sign.value_or(defaults->fmt_sign);
      auto prec = float_precision.value_or(defaults->float_precision);
      auto fmt_str = fmt::format("{{:>{sign:s}{width}.{prec}e}}", fmt::arg("sign", sign),
                                 fmt::arg("width", width_), fmt::arg("prec", prec));
      return fmt::format(fmt::runtime(fmt_str), val);
    }
  }
  auto empty_cell = empty_cell_val.value_or(defaults->empty_cell_val);
  return fmt::format("{0:>{1}s}", empty_cell, width_);
}

Column::Column(std::string name, std::string header_text,
               std::optional<size_t> min_left_padding,
               std::optional<size_t> float_precision,
               std::optional<std::string> empty_cell_val,
               std::optional<std::string> fmt_sign)
  : name(std::move(name)), header_text(std::move(header_text)),
    min_left_padding(min_left_padding), float_precision(float_precision),
    empty_cell_val(std::move(empty_cell_val)), fmt_sign(std::move(fmt_sign))
{
}

[[nodiscard]] size_t Table::n_rows() const
{
  if (n_cols() == 0)
  {
    return 0;
  }
  auto max_col = std::max_element(cols.begin(), cols.end(), [](const auto &a, const auto &b)
                                  { return a.n_rows() < b.n_rows(); });
  return max_col->n_rows();
}

void Table::reserve(size_t n_rows, size_t n_cols)
{
  reserve_n_rows = n_rows;
  cols.reserve(n_cols);
  for (auto &col : cols)
  {
    col.data.reserve(n_rows);
  }
}

// Insert columns: map like interface
bool Table::insert(Column &&column)
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

Column &Table::operator[](std::string_view name)
{
  auto it =
      std::find_if(cols.begin(), cols.end(), [&name](auto &c) { return c.name == name; });
  if (it == cols.end())
  {
    throw std::out_of_range(fmt::format("Column {} not found in table", name).c_str());
  }
  return *it;
}

// TODO: Improve all the functions below with ranges in C++20
template <typename T>
void Table::append_header(T &buf) const
{
  auto to = [&buf](auto f, auto &&...a)
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  for (size_t i = 0; i < n_cols(); i++)
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
void Table::append_row(T &buf, size_t row_j) const
{
  auto to = [&buf](auto f, auto &&...a)
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  for (size_t i = 0; i < n_cols(); i++)
  {
    if (i > 0)
    {
      to("{:s}", print_col_separator);
    }
    to("{:s}", cols[i].format_row(row_j));
  }
  to("{:s}", print_row_separator);
}

[[nodiscard]] std::string Table::format_header() const
{
  fmt::memory_buffer buf{};
  append_header(buf);
  return {buf.data(), buf.size()};
}

[[nodiscard]] std::string Table::format_row(size_t j) const
{
  fmt::memory_buffer buf{};
  append_row(buf, j);
  return {buf.data(), buf.size()};
}

[[nodiscard]] std::string Table::format_table() const
{
  fmt::memory_buffer buf{};
  append_header(buf);
  for (size_t j = 0; j < n_rows(); j++)
  {
    append_row(buf, j);
  }
  return {buf.data(), buf.size()};
}

// explicit instantiation to avoid fmt inclusion
template void Table::append_header(fmt::memory_buffer &) const;
template void Table::append_row(fmt::memory_buffer &, size_t) const;

TableWithCSVFile::TableWithCSVFile(std::string csv_file_fullpath)
  : csv_file_fullpath_{std::move(csv_file_fullpath)}
{
  // Validate
  auto file_buf = fmt::output_file(
      csv_file_fullpath_, fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);
}

void TableWithCSVFile::WriteFullTableTrunc()
{
  auto file_buf = fmt::output_file(
      csv_file_fullpath_, fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);
  file_buf.print("{}", table.format_table());
  file_append_cursor = table.n_rows();
}

void TableWithCSVFile::AppendHeader()
{
  if (file_append_cursor != -1)
  {
    WriteFullTableTrunc();
  }
  auto file_buf =
      fmt::output_file(csv_file_fullpath_, fmt::file::WRONLY | fmt::file::APPEND);
  file_buf.print("{}", table.format_header());
  file_append_cursor++;
}

void TableWithCSVFile::AppendRow()
{
  if (file_append_cursor < 0)
  {
    WriteFullTableTrunc();
  }
  auto file_buf =
      fmt::output_file(csv_file_fullpath_, fmt::file::WRONLY | fmt::file::APPEND);
  file_buf.print("{}", table.format_row(file_append_cursor));
  file_append_cursor++;
}

}  // namespace palace
