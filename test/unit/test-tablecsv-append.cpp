// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/tablecsv.hpp"

using namespace palace;

namespace
{

std::string ReadFile(const fs::path &path)
{
  std::ifstream file(path, std::ios_base::in | std::ios_base::binary);
  REQUIRE(file.good());
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

// The per-patch tables' column form (the ES patch table and the Maxwell patch-energy
// table): integer index columns, float columns with the default and a narrower precision,
// nan cells.
void AddColumns(Table &table)
{
  table.insert(Column::Integer("idx", "idx"));
  table.insert(Column::Integer("patch", "patch"));
  table.insert("weight", "weight");
  table.insert("energy", "E_surf (J)", 0, 8, 6);
  table.insert("ratio", "ratio", 0, 0, 3, "");
}

void AddRow(Table &table, int idx, int patch, double weight, double energy, double ratio)
{
  table["idx"] << idx;
  table["patch"] << patch;
  table["weight"] << weight;
  table["energy"] << energy;
  table["ratio"] << ratio;
}

}  // namespace

// Lane fix-response-memory (decisions 1439 / 1446; DESIGN s2.2 / s3 (M2), commit C1): a
// table written in pieces by TableWithCSVFile::AppendRowsAndClear is byte for byte the
// WriteFullTableTrunc of the whole table (the column widths are data-independent).
TEST_CASE("TableCSV AppendRowsAndClear reproduces WriteFullTableTrunc byte for byte",
          "[tablecsv][memfix][Serial]")
{
  // Serial: one writer per file.
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const fs::path dir = fs::temp_directory_path() / "palace-memfix-tablecsv-append";
  fs::create_directories(dir);
  const fs::path whole_path = dir / "whole.csv", pieces_path = dir / "pieces.csv";
  fs::remove(whole_path);
  fs::remove(pieces_path);

  // The reference: every row, one write.
  TableWithCSVFile whole(whole_path.string());
  AddColumns(whole.table);
  // Three "sources" of 3 / 1 / 4 rows (the third holds nan cells and a wide index).
  AddRow(whole.table, 0, 7, 1.0, 1.23456789012345e-18, 0.5);
  AddRow(whole.table, 0, 8, 0.5, 2.5e-19, 0.25);
  AddRow(whole.table, 0, 1024, 0.0, 0.0, 0.0);
  AddRow(whole.table, 1, 7, 1.0, -3.0e-18, 1.0e3);
  AddRow(whole.table, 2, 7, nan, nan, nan);
  AddRow(whole.table, 2, 8, 0.5, 2.5e-19, 0.25);
  AddRow(whole.table, 2, 123456, 1.0, 1.0e-300, 1.0e300);
  AddRow(whole.table, 2, 9, 0.25, 4.0e-19, 0.125);
  whole.WriteFullTableTrunc();
  const std::string reference = ReadFile(whole_path);
  REQUIRE(!reference.empty());

  // The same rows in three pieces (the first piece truncates + writes the header; a
  // stale file at the path first, to show the truncation).
  {
    std::ofstream stale(pieces_path);
    stale << "stale content\n";
  }
  TableWithCSVFile pieces(pieces_path.string());
  AddColumns(pieces.table);
  AddRow(pieces.table, 0, 7, 1.0, 1.23456789012345e-18, 0.5);
  AddRow(pieces.table, 0, 8, 0.5, 2.5e-19, 0.25);
  AddRow(pieces.table, 0, 1024, 0.0, 0.0, 0.0);
  pieces.AppendRowsAndClear();
  CHECK(pieces.appended_rows() == 3);
  CHECK(pieces.table.n_rows() == 0);
  CHECK(pieces.table.n_cols() == 5);
  AddRow(pieces.table, 1, 7, 1.0, -3.0e-18, 1.0e3);
  pieces.AppendRowsAndClear();
  CHECK(pieces.appended_rows() == 4);
  // An empty piece writes nothing (and does not re-truncate).
  pieces.AppendRowsAndClear();
  CHECK(pieces.appended_rows() == 4);
  AddRow(pieces.table, 2, 7, nan, nan, nan);
  AddRow(pieces.table, 2, 8, 0.5, 2.5e-19, 0.25);
  AddRow(pieces.table, 2, 123456, 1.0, 1.0e-300, 1.0e300);
  AddRow(pieces.table, 2, 9, 0.25, 4.0e-19, 0.125);
  pieces.AppendRowsAndClear();
  CHECK(pieces.appended_rows() == 8);
  CHECK(ReadFile(pieces_path) == reference);

  // An empty FIRST piece writes the header only (today's file for a table without rows)
  // and the rows then follow it.
  const fs::path header_first_path = dir / "header-first.csv";
  fs::remove(header_first_path);
  TableWithCSVFile header_first(header_first_path.string());
  AddColumns(header_first.table);
  header_first.AppendRowsAndClear();
  CHECK(header_first.appended_rows() == 0);
  {
    TableWithCSVFile header_only(whole_path.string());
    AddColumns(header_only.table);
    header_only.WriteFullTableTrunc();
    CHECK(ReadFile(header_first_path) == ReadFile(whole_path));
  }
  for (int piece = 0; piece < 2; piece++)
  {
    AddRow(header_first.table, 0, 7, 1.0, 1.23456789012345e-18, 0.5);
    AddRow(header_first.table, 0, 8, 0.5, 2.5e-19, 0.25);
    AddRow(header_first.table, 0, 1024, 0.0, 0.0, 0.0);
    AddRow(header_first.table, 1, 7, 1.0, -3.0e-18, 1.0e3);
    header_first.AppendRowsAndClear();
  }
  for (int piece = 0; piece < 1; piece++)
  {
    AddRow(header_first.table, 2, 7, nan, nan, nan);
    AddRow(header_first.table, 2, 8, 0.5, 2.5e-19, 0.25);
    AddRow(header_first.table, 2, 123456, 1.0, 1.0e-300, 1.0e300);
    AddRow(header_first.table, 2, 9, 0.25, 4.0e-19, 0.125);
    header_first.AppendRowsAndClear();
  }
  // = the reference's header + the first 4 rows twice + the last 4 rows: assembled from the
  // reference's lines.
  std::string expected;
  {
    std::vector<std::string> lines;
    std::stringstream stream(reference);
    for (std::string line; std::getline(stream, line);)
    {
      lines.push_back(line + "\n");
    }
    REQUIRE(lines.size() == 9);
    expected = lines[0];
    for (int piece = 0; piece < 2; piece++)
    {
      for (std::size_t j = 1; j <= 4; j++)
      {
        expected += lines[j];
      }
    }
    for (std::size_t j = 5; j <= 8; j++)
    {
      expected += lines[j];
    }
  }
  CHECK(ReadFile(header_first_path) == expected);
  fs::remove_all(dir);
}
