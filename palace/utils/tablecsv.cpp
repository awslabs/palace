// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "tablecsv.hpp"
#include "utils/filesystem.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <utility>
#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/ranges.h>
#include <scn/scan.h>
#include <mfem/general/error.hpp>

namespace palace
{

[[nodiscard]] size_t Column::col_width(const ColumnOptions &defaults) const
{
  // Quickfix to specify full column width in integer case to match current formatting.
  if (print_as_int)
  {
    return std::max(min_left_padding.value_or(defaults.min_left_padding),
                    header_text.size());
  }
  size_t pad = min_left_padding.value_or(defaults.min_left_padding);
  size_t prec = float_precision.value_or(defaults.float_precision);

  // Normal float in our exponent format needs float_precision + 7 ("+" , leading digit,
  // ".", "e", "+", +2 exponent. Sometimes exponent maybe +3 if very small or large; see
  // std::numeric_limits<double>::max_exponent. We pick +7 for consistency, but
  // min_left_padding should be at least 1, which is not currently enforced.
  return std::max(pad + prec + 7, header_text.size());
}

[[nodiscard]] auto Column::format_header(const ColumnOptions &defaults,
                                         const std::optional<size_t> &width) const
{
  auto w = width.value_or(col_width(defaults));
  return fmt::format("{0:>{1}s}", header_text, w);
}

[[nodiscard]] auto Column::format_row(size_t i, const ColumnOptions &defaults,
                                      const std::optional<size_t> &width) const
{
  auto width_ = width.value_or(col_width(defaults));
  // If data available format double.
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
      auto sign = fmt_sign.value_or(defaults.fmt_sign);
      auto prec = float_precision.value_or(defaults.float_precision);
      auto fmt_str = fmt::format("{{:>{sign:s}{width}.{prec}e}}", fmt::arg("sign", sign),
                                 fmt::arg("width", width_), fmt::arg("prec", prec));
      return fmt::format(fmt::runtime(fmt_str), val);
    }
  }
  return fmt::format("{0:>{1}s}", defaults.empty_cell_val, width_);
}

Column::Column(std::string name_, std::string header_text_, long column_group_idx_,
               std::optional<size_t> min_left_padding_,
               std::optional<size_t> float_precision_, std::optional<std::string> fmt_sign_)
  : name(std::move(name_)), header_text(std::move(header_text_)),
    column_group_idx(column_group_idx_), min_left_padding(min_left_padding_),
    float_precision(float_precision_), fmt_sign(std::move(fmt_sign_))
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

// Insert columns: map like interface.
bool Table::insert(Column &&column)
{
  auto it = std::find_if(cols.begin(), cols.end(),
                         [&column](auto &c) { return c.name == column.name; });
  if (it != cols.end())
  {
    return false;
  }
  auto &col = cols.emplace_back(std::move(column));
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

// TODO: Improve all the functions below with ranges in C++20.
template <typename T>
void Table::append_header(T &buf) const
{
  for (size_t i = 0; i < n_cols(); i++)
  {
    if (i > 0)
    {
      fmt::format_to(std::back_inserter(buf), "{:s}", print_col_separator);
    }
    fmt::format_to(std::back_inserter(buf), "{:s}", cols[i].format_header(col_options));
  }
  fmt::format_to(std::back_inserter(buf), "{:s}", print_row_separator);
}

template <typename T>
void Table::append_row(T &buf, size_t row_j) const
{
  for (size_t i = 0; i < n_cols(); i++)
  {
    if (i > 0)
    {
      fmt::format_to(std::back_inserter(buf), "{:s}", print_col_separator);
    }
    fmt::format_to(std::back_inserter(buf), "{:s}", cols[i].format_row(row_j, col_options));
  }
  fmt::format_to(std::back_inserter(buf), "{:s}", print_row_separator);
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

// Helper class to parse split view between delimiters. Not a proper view, but same idea.
class sv_split_r
{
  std::string_view full_view;
  std::string_view delimiters;
  std::size_t row_cursor_start = 0;
  std::size_t row_cursor_end = 0;

public:
  sv_split_r(std::string_view full_view_, std::string_view delimiters_)
    : full_view(full_view_), delimiters(delimiters_)
  {
  }

  std::string_view next()
  {
    row_cursor_end = full_view.find_first_of(delimiters, row_cursor_start);
    auto out = full_view.substr(row_cursor_start, row_cursor_end - row_cursor_start);
    row_cursor_start =
        (row_cursor_end == std::string_view::npos) ? row_cursor_end : row_cursor_end + 1;
    return out;
  }

  bool at_end() const
  {
    return (row_cursor_start == full_view.size()) ||
           (row_cursor_start == std::string_view::npos);
  }
};

std::string_view trim_space(std::string_view str_v)
{
  auto prefix = str_v.find_first_not_of(" \t\f\v");
  str_v.remove_prefix(std::min(prefix, str_v.size()));
  auto suffix = str_v.find_last_not_of(" \t\f\v");  // suffix also counts from 0
  str_v.remove_suffix(std::min(str_v.size() - 1 - suffix, str_v.size()));
  return str_v;
}

Table::Table(std::string_view table_str,
             std::optional<std::string_view> print_col_separator_,
             std::optional<std::string_view> print_row_separator_)
{
  using namespace std::literals;
  if (table_str.empty())
  {
    return;
  }
  try
  {
    print_col_separator = print_col_separator_.value_or(print_col_separator);
    print_row_separator = print_row_separator_.value_or(print_row_separator);

    std::string_view table_full_view{table_str};
    sv_split_r row_split(table_full_view, print_row_separator);
    // Handle first row separately since it defines header & cols.
    {
      sv_split_r entries(row_split.next(), print_col_separator);
      std::size_t col_i = 0;
      while (!entries.at_end())
      {
        std::string entry_trim_str{trim_space(entries.next())};
        this->insert(fmt::format("col_{}", col_i), entry_trim_str);
        col_i++;
      }
    }
    // Loop over other rows
    while (!row_split.at_end())
    {
      sv_split_r entries(row_split.next(), print_col_separator);
      std::size_t col_i = 0;
      while (!entries.at_end())
      {
        auto entry_trim = trim_space(entries.next());
        if ((entry_trim == col_options.empty_cell_val) || (entry_trim.size() == 0))
        {
        }
        else
        {
          auto result = scn::scan<double>(entry_trim, "{}");
          MFEM_VERIFY(result, fmt::format("Could not parse CSV entry \"{}\" as double",
                                          entry_trim));
          (*this)[col_i] << result->value();
        }
        col_i++;
      }
    }
  }
  catch (const std::exception &e)
  {
    MFEM_ABORT("Could not parse CSV table from string!\n  " << e.what());
  }
}

// explicit instantiation to avoid fmt inclusion.
template void Table::append_header(fmt::memory_buffer &) const;
template void Table::append_row(fmt::memory_buffer &, size_t) const;

TableWithCSVFile::TableWithCSVFile(std::string csv_file_fullpath, bool load_existing_file)
  : csv_file_fullpath_{std::move(csv_file_fullpath)}
{
  if (!load_existing_file)
  {
    return;
  }
  if (!fs::exists(csv_file_fullpath_))
  {
    return;
  }

  std::ifstream file_buffer(csv_file_fullpath_, std::ios_base::in);
  if (!file_buffer.good())
  {
    return;
  }
  std::stringstream file_buffer_str;
  file_buffer_str << file_buffer.rdbuf();
  file_buffer.close();
  table = Table(file_buffer_str.str());
}

void TableWithCSVFile::WriteFullTableTrunc()
{
  auto file_buffer = fmt::output_file(
      csv_file_fullpath_, fmt::file::WRONLY | fmt::file::CREATE | fmt::file::TRUNC);
  file_buffer.print("{}", table.format_table());
}

}  // namespace palace
