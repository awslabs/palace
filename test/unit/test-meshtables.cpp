// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>

// The edge (mfem::DSTable) and face (mfem::STable3D) tables of mfem::Mesh number their
// entries in the order they are first pushed: Palace's meshes (and so the finite element
// dof numbering) depend on it. Check the tables against a reference, including tables
// growing past their initial capacity.

TEST_CASE("DSTable numbering and lookups", "[meshtables][Serial]")
{
  std::mt19937 gen(1234);
  for (int nrows : {1, 7, 1000, 50000})
  {
    std::uniform_int_distribution<int> row(0, nrows - 1);
    mfem::DSTable table(nrows);
    std::map<std::pair<int, int>, int> ref;
    const int npush = 20 * nrows + 10;
    for (int i = 0; i < npush; i++)
    {
      int a = row(gen), b = row(gen);
      const auto key = std::minmax(a, b);
      auto it = ref.find(key);
      const int expected = (it == ref.end()) ? static_cast<int>(ref.size()) : it->second;
      if (it == ref.end())
      {
        ref.emplace(key, expected);
      }
      REQUIRE(table.Push(a, b) == expected);
    }
    REQUIRE(table.NumberOfEntries() == static_cast<int>(ref.size()));
    for (const auto &[key, idx] : ref)
    {
      CHECK(table(key.first, key.second) == idx);
      CHECK(table(key.second, key.first) == idx);
    }
    // Entries not in the table, and rows out of range.
    for (int i = 0; i < 1000; i++)
    {
      int a = row(gen), b = row(gen);
      if (!ref.contains(std::minmax(a, b)))
      {
        CHECK(table(a, b) == -1);
      }
    }
    CHECK(table(nrows, 0) == -1);
    // Row iteration: the entries of each row, last pushed first.
    std::map<int, std::vector<std::pair<int, int>>> rows;
    for (const auto &[key, idx] : ref)
    {
      rows[key.first].emplace_back(idx, key.second);
    }
    for (int r = 0; r < nrows; r++)
    {
      std::vector<std::pair<int, int>> expected = rows[r], found;
      std::sort(expected.rbegin(), expected.rend());
      for (mfem::DSTable::RowIterator it(table, r); !it; ++it)
      {
        found.emplace_back(it.Index(), it.Column());
      }
      CHECK(found == expected);
    }
    // Indices changed through the row iterator are seen by the lookups.
    for (mfem::DSTable::RowIterator it(table, 0); !it; ++it)
    {
      it.SetIndex(-2 - it.Column());
    }
    for (mfem::DSTable::RowIterator it(table, 0); !it; ++it)
    {
      CHECK(table(0, it.Column()) == -2 - it.Column());
    }
  }
}

TEST_CASE("STable3D numbering and lookups", "[meshtables][Serial]")
{
  std::mt19937 gen(4321);
  for (int nrows : {4, 9, 1000, 20000})
  {
    std::uniform_int_distribution<int> row(0, nrows - 1);
    mfem::STable3D table(nrows);
    std::map<std::array<int, 3>, int> ref;
    auto Key = [](int a, int b, int c)
    {
      std::array<int, 3> k = {a, b, c};
      std::sort(k.begin(), k.end());
      return k;
    };
    const int npush = 20 * nrows + 10;
    for (int i = 0; i < npush; i++)
    {
      int a = row(gen), b = row(gen), c = row(gen);
      if (a == b || b == c || a == c)
      {
        continue;
      }
      const auto key = Key(a, b, c);
      auto it = ref.find(key);
      const int expected = (it == ref.end()) ? static_cast<int>(ref.size()) : it->second;
      if (it == ref.end())
      {
        ref.emplace(key, expected);
      }
      REQUIRE(table.Push(a, b, c) == expected);
    }
    REQUIRE(table.NumberOfElements() == static_cast<int>(ref.size()));
    for (const auto &[key, idx] : ref)
    {
      CHECK(table(key[2], key[0], key[1]) == idx);
      CHECK(table.Index(key[1], key[2], key[0]) == idx);
    }
    for (int i = 0; i < 1000; i++)
    {
      int a = row(gen), b = row(gen), c = row(gen);
      if (a != b && b != c && a != c && !ref.contains(Key(a, b, c)))
      {
        CHECK(table.Index(a, b, c) == -1);
      }
    }
    CHECK(table.Index(nrows, 0, 1) == -1);
    // Quadrilateral faces are addressed by their three smallest vertices.
    if (nrows >= 9)
    {
      const int n0 = table.NumberOfElements();
      const int q = table.Push4(nrows - 1, 0, 1, 2);
      CHECK(q == table.Push4(2, 1, nrows - 1, 0));
      CHECK(q == table(0, 1, 2, nrows - 1));
      CHECK(q == table.Index(2, 0, 1));
      CHECK(table.NumberOfElements() == n0 + !ref.contains(Key(0, 1, 2)));
    }
  }
}
