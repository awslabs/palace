// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fixtures.hpp"
#include "utils/meshio.hpp"

using namespace palace;

namespace
{

// Pad or truncate to one exactly-8-character Nastran fixed-width field.
std::string Field8(std::string s)
{
  s.resize(8, ' ');
  return s;
}

// Parse the node coordinates out of a Gmsh v2.2 buffer, ASCII or binary depending on how
// the build configured GMSH_BIN.
std::vector<std::array<double, 3>> ParseGmshNodes(const std::string &gmsh)
{
  const bool binary = gmsh.find("2.2 1 8\n") != std::string::npos;
  const auto nodes_pos = gmsh.find("$Nodes\n");
  if (nodes_pos == std::string::npos)
  {
    return {};
  }
  const std::size_t count_pos = nodes_pos + 7;
  const int num_nodes =
      std::stoi(gmsh.substr(count_pos, gmsh.find('\n', count_pos) - count_pos));
  std::vector<std::array<double, 3>> nodes;
  if (!binary)
  {
    std::istringstream in(gmsh.substr(gmsh.find('\n', count_pos) + 1));
    for (int i = 0; i < num_nodes; i++)
    {
      int tag;
      std::array<double, 3> x;
      in >> tag >> x[0] >> x[1] >> x[2];
      nodes.push_back(x);
    }
    return nodes;
  }
  // Binary node records: int tag followed by three doubles.
  std::size_t p = gmsh.find('\n', count_pos) + 1;
  for (int i = 0; i < num_nodes; i++)
  {
    int tag;
    std::array<double, 3> x;
    std::memcpy(&tag, gmsh.data() + p, sizeof(int));
    p += sizeof(int);
    std::memcpy(x.data(), gmsh.data() + p, 3 * sizeof(double));
    p += 3 * sizeof(double);
    nodes.push_back(x);
  }
  return nodes;
}

struct GmshElement
{
  int type, geom;
  std::vector<int> nodes;
};

// Parse the elements out of a Gmsh v2.2 buffer, for the element types the Nastran tests
// use.
std::vector<GmshElement> ParseGmshElements(const std::string &gmsh)
{
  auto NumNodes = [](int type)
  {
    switch (type)
    {
      case 2:  // 3-node triangle
        return 3;
      case 9:  // 6-node triangle
        return 6;
      case 11:  // 10-node tetrahedron
        return 10;
      default:
        FAIL("Unexpected Gmsh element type " << type);
        return 0;
    }
  };
  const bool binary = gmsh.find("2.2 1 8\n") != std::string::npos;
  const auto elems_pos = gmsh.find("$Elements\n");
  if (elems_pos == std::string::npos)
  {
    return {};
  }
  const std::size_t count_pos = elems_pos + 10;
  const int num_elems =
      std::stoi(gmsh.substr(count_pos, gmsh.find('\n', count_pos) - count_pos));
  std::vector<GmshElement> elems;
  if (!binary)
  {
    std::istringstream in(gmsh.substr(gmsh.find('\n', count_pos) + 1));
    for (int i = 0; i < num_elems; i++)
    {
      int tag, num_tags;
      GmshElement e;
      in >> tag >> e.type >> num_tags >> e.geom;
      for (int j = 1; j < num_tags; j++)
      {
        in >> tag;
      }
      e.nodes.resize(NumNodes(e.type));
      for (auto &n : e.nodes)
      {
        in >> n;
      }
      elems.push_back(e);
    }
    return elems;
  }
  // Binary element blocks: header [type, count, number of tags (2)], then per element
  // [tag, physical tag, geometry tag, nodes].
  auto ReadInt = [&gmsh, p = gmsh.find('\n', count_pos) + 1]() mutable
  {
    int i;
    std::memcpy(&i, gmsh.data() + p, sizeof(int));
    p += sizeof(int);
    return i;
  };
  while (static_cast<int>(elems.size()) < num_elems)
  {
    const int type = ReadInt(), count = ReadInt(), num_tags = ReadInt();
    for (int i = 0; i < count; i++)
    {
      GmshElement e{type, 0, std::vector<int>(NumNodes(type))};
      ReadInt();  // Element tag
      e.geom = ReadInt();
      for (int j = 1; j < num_tags; j++)
      {
        ReadInt();
      }
      for (auto &n : e.nodes)
      {
        n = ReadInt();
      }
      elems.push_back(e);
    }
  }
  return elems;
}

// A single second-order tetrahedron and one of its faces, exported by Gmsh 4.15.2 in each
// of its three Nastran formats. Gmsh writes no BEGIN BULK line, and continues the 10-node
// CTETRA with a "+E2" marker in the free and small field formats and with a blank first
// field in the large field format.

// Gmsh free field format (Mesh.BdfFieldFormat = 0).
constexpr const char *gmsh_bdf_free = "$ Created by Gmsh\n"
                                      "GRID,1,0,0.00E+00,0.00E+00,0.00E+00\n"
                                      "GRID,2,0,1.000000,0.00E+00,0.00E+00\n"
                                      "GRID,3,0,0.00E+00,1.000000,0.00E+00\n"
                                      "GRID,4,0,0.00E+00,0.00E+00,1.000000\n"
                                      "GRID,5,0,0.500000,0.00E+00,0.00E+00\n"
                                      "GRID,6,0,0.500000,0.500000,0.00E+00\n"
                                      "GRID,7,0,0.00E+00,0.500000,0.00E+00\n"
                                      "GRID,8,0,0.00E+00,0.00E+00,0.500000\n"
                                      "GRID,9,0,0.500000,0.00E+00,0.500000\n"
                                      "GRID,10,0,0.00E+00,0.500000,0.500000\n"
                                      "CTRIA6,1,1,1,2,3,5,6,7\n"
                                      "CTETRA,2,1,1,3,4,2,7,10,+E2\n"
                                      "+E2,8,5,6,9\n"
                                      "ENDDATA\n";

// Gmsh small field format (Mesh.BdfFieldFormat = 1).
constexpr const char *gmsh_bdf_small =
    "$ Created by Gmsh\n"
    "GRID    1       0       0.00E+000.00E+000.00E+00\n"
    "GRID    2       0       1.0000000.00E+000.00E+00\n"
    "GRID    3       0       0.00E+001.0000000.00E+00\n"
    "GRID    4       0       0.00E+000.00E+001.000000\n"
    "GRID    5       0       0.5000000.00E+000.00E+00\n"
    "GRID    6       0       0.5000000.5000000.00E+00\n"
    "GRID    7       0       0.00E+000.5000000.00E+00\n"
    "GRID    8       0       0.00E+000.00E+000.500000\n"
    "GRID    9       0       0.5000000.00E+000.500000\n"
    "GRID    10      0       0.00E+000.5000000.500000\n"
    "CTRIA6  1       1       1       2       3       5       6       7       \n"
    "CTETRA  2       1       1       3       4       2       7       10      +E2     \n"
    "+E2     8       5       6       9       \n"
    "ENDDATA\n";

// Gmsh large field format (Mesh.BdfFieldFormat = 2).
constexpr const char *gmsh_bdf_large =
    "$ Created by Gmsh\n"
    "GRID*   1               0               0.00000000      0.00000000      \n"
    "*       0.00000000      \n"
    "GRID*   2               0               1.00000000      0.00000000      \n"
    "*       0.00000000      \n"
    "GRID*   3               0               0.00000000      1.00000000      \n"
    "*       0.00000000      \n"
    "GRID*   4               0               0.00000000      0.00000000      \n"
    "*       1.00000000      \n"
    "GRID*   5               0               0.500000000     0.00000000      \n"
    "*       0.00000000      \n"
    "GRID*   6               0               0.500000000     0.500000000     \n"
    "*       0.00000000      \n"
    "GRID*   7               0               0.00000000      0.500000000     \n"
    "*       0.00000000      \n"
    "GRID*   8               0               0.00000000      0.00000000      \n"
    "*       0.500000000     \n"
    "GRID*   9               0               0.500000000     0.00000000      \n"
    "*       0.500000000     \n"
    "GRID*   10              0               0.00000000      0.500000000     \n"
    "*       0.500000000     \n"
    "CTRIA6  1       1       1       2       3       5       6       7       \n"
    "CTETRA  2       1       1       3       4       2       7       10      \n"
    "        8       5       6       9       \n"
    "ENDDATA\n";

}  // namespace

// Nastran fixed-width fields use implicit exponents ("-7.-1" for "-7.E-01"), Fortran 'D'
// exponents, and blank fields (zero by convention). std::stod parses only the valid prefix
// of the implicit forms without an error, so a prefix-only parse silently truncates the
// value; this test guards the full-consumption fix-up in ConvertDoubleNastran.
TEST_CASE_METHOD(palace::test::SharedTempDir, "Nastran special floating point formats",
                 "[meshio][Serial]")
{
  std::ostringstream nas;
  nas << "BEGIN BULK\n";
  // Node 1: implicit exponents for x and y, blank z (zero by Nastran convention).
  nas << Field8("GRID") + Field8("1") + Field8("") + Field8("-7.-1") + Field8("2.3+2") +
             Field8("") + "\n";
  // Node 2: Fortran 'D' exponent.
  nas << Field8("GRID") + Field8("2") + Field8("") + Field8("1.5D+1") + Field8("0.0") +
             Field8("0.0") + "\n";
  // Node 3: ordinary fixed-point fields, with the line truncated after the y field as
  // produced by writers that trim trailing blanks (the absent z field reads as zero).
  nas << Field8("GRID") + Field8("3") + Field8("") + Field8("0.0") + "1.0" + "\n";
  nas << Field8("CTRIA3") + Field8("1") + Field8("1") + Field8("1") + Field8("2") +
             Field8("3") + "\n";
  nas << "ENDDATA\n";

  auto path = temp_dir / "special_float.nas";
  {
    std::ofstream f(path);
    f << nas.str();
  }

  std::stringstream buffer;
  mesh::ConvertMeshNastran(path.string(), buffer);
  auto nodes = ParseGmshNodes(buffer.str());

  REQUIRE(nodes.size() == 3);
  CHECK(nodes[0][0] == Catch::Approx(-0.7));   // "-7.-1" = -7.E-01
  CHECK(nodes[0][1] == Catch::Approx(230.0));  // "2.3+2" = 2.3E+02
  CHECK(nodes[0][2] == Catch::Approx(0.0));    // Blank field
  CHECK(nodes[1][0] == Catch::Approx(15.0));   // "1.5D+1"
  CHECK(nodes[2][1] == Catch::Approx(1.0));    // Short line
  CHECK(nodes[2][2] == Catch::Approx(0.0));    // Absent (unpadded) field
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Nastran invalid input aborts",
                 "[meshio][Serial]")
{
  auto Convert = [this](const std::string &bulk)
  {
    auto path = temp_dir / "invalid.nas";
    {
      std::ofstream f(path);
      f << "BEGIN BULK\n" << bulk;
    }
    std::stringstream buffer;
    mesh::ConvertMeshNastran(path.string(), buffer);
  };
  using Catch::Matchers::ContainsSubstring;
  const std::string tria = "CTRIA3,1,1,1,2,3\nENDDATA\n";
  SECTION("Invalid real")
  {
    CHECK_THROWS_WITH(Convert(Field8("GRID") + Field8("1") + Field8("") + Field8("1.2.3") +
                              Field8("0.0") + Field8("0.0") + "\n" + tria),
                      ContainsSubstring("Invalid number"));
  }
  SECTION("Invalid integer")
  {
    CHECK_THROWS_WITH(Convert("GRID,1,,0.,0.,0.\nCTRIA3,1,1,1,2,x\nENDDATA\n"),
                      ContainsSubstring("Invalid integer"));
  }
  // Palace does not apply coordinate systems, so a GRID point in one (nonzero CP) would be
  // placed wrongly.
  SECTION("GRID coordinate system")
  {
    CHECK_THROWS_WITH(Convert("GRID,1,2,0.,0.,0.\n" + tria),
                      ContainsSubstring("coordinate system"));
    CHECK_THROWS_WITH(Convert(Field8("GRID") + Field8("1") + Field8("2") + Field8("0.") +
                              Field8("0.") + Field8("0.") + "\n" + tria),
                      ContainsSubstring("coordinate system"));
    CHECK_THROWS_WITH(Convert("GRID*   1               2               0.              0.\n"
                              "*       0.\n" +
                              tria),
                      ContainsSubstring("coordinate system"));
    // A GRDSET CP applies to every GRID point with a blank CP field.
    CHECK_THROWS_WITH(Convert(Field8("GRDSET") + Field8("") + Field8("2") +
                              "\nGRID,1,,0.,0.,0.\n" + tria),
                      ContainsSubstring("coordinate system"));
    CHECK_THROWS_WITH(Convert("GRDSET*                 2\nGRID,1,,0.,0.,0.\n" + tria),
                      ContainsSubstring("coordinate system"));
    CHECK_NOTHROW(Convert(Field8("GRDSET") + Field8("") + Field8("0") +
                          "\nGRID,1,,0.,0.,0.\n" + "GRID,2,,1.,0.,0.\nGRID,3,,0.,1.,0.\n" +
                          tria));
  }
  SECTION("Missing ENDDATA")
  {
    CHECK_THROWS_WITH(Convert("GRID,1,,0.,0.,0.\nCTRIA3,1,1,1,2,3\n"),
                      ContainsSubstring("ENDDATA"));
  }
}

// Free field (comma-separated) cards may omit trailing fields, which are then blank. Gmsh
// writes CTRIA3 cards this way.
TEST_CASE_METHOD(palace::test::SharedTempDir, "Nastran free field omitted trailing fields",
                 "[meshio][Serial]")
{
  auto path = temp_dir / "free_field.nas";
  {
    std::ofstream f(path);
    f << "BEGIN BULK\n"
         "GRID,1,,1.,2.\n"   // z omitted
         "GRID,2,,3.,4.,\n"  // z blank
         "GRID,3,,5.,6.,7.\n"
         "CTRIA3,1,1,1,2,3\n"
         "ENDDATA\n";
  }
  std::stringstream buffer;
  mesh::ConvertMeshNastran(path.string(), buffer);
  auto nodes = ParseGmshNodes(buffer.str());
  auto elems = ParseGmshElements(buffer.str());

  REQUIRE(nodes.size() == 3);
  CHECK(nodes[0][2] == 0.0);
  CHECK(nodes[1][2] == 0.0);
  CHECK(nodes[2][2] == 7.0);
  REQUIRE(elems.size() == 1);
  CHECK(elems[0].type == 2);  // 3-node triangle
  CHECK(elems[0].nodes == std::vector<int>{1, 2, 3});
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Nastran meshes exported by Gmsh",
                 "[meshio][Serial]")
{
  const auto [format, bdf] = GENERATE(table<std::string, std::string>(
      {{"free", gmsh_bdf_free}, {"small", gmsh_bdf_small}, {"large", gmsh_bdf_large}}));
  // As Gmsh writes it; with a BEGIN BULK line, to check the cards are parsed the same
  // either way; and without the final newline, which some editors and writers omit.
  const auto variant = GENERATE("as written", "with BEGIN BULK", "no final newline");
  std::string nas = bdf;
  if (std::string(variant) != "as written")
  {
    nas = "BEGIN BULK\n" + nas;
  }
  if (std::string(variant) == "no final newline")
  {
    nas.pop_back();
  }
  CAPTURE(format, variant);

  auto path = temp_dir / "gmsh.bdf";
  {
    std::ofstream f(path);
    f << nas;
  }
  std::stringstream buffer;
  mesh::ConvertMeshNastran(path.string(), buffer);
  auto nodes = ParseGmshNodes(buffer.str());
  auto elems = ParseGmshElements(buffer.str());

  // Gmsh's own MSH 2.2 export of the same mesh.
  const std::vector<std::array<double, 3>> expected_nodes = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.5, 0.0, 0.0},
      {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
  REQUIRE(nodes.size() == expected_nodes.size());
  for (std::size_t i = 0; i < nodes.size(); i++)
  {
    for (int d = 0; d < 3; d++)
    {
      CHECK(nodes[i][d] == Catch::Approx(expected_nodes[i][d]));
    }
  }
  REQUIRE(elems.size() == 2);
  std::ranges::sort(elems, {}, &GmshElement::type);
  CHECK(elems[0].type == 9);  // 6-node triangle
  CHECK(elems[0].geom == 1);
  CHECK(elems[0].nodes == std::vector<int>{1, 2, 3, 5, 6, 7});
  CHECK(elems[1].type == 11);  // 10-node tetrahedron
  CHECK(elems[1].geom == 1);
  CHECK(elems[1].nodes == std::vector<int>{1, 3, 4, 2, 7, 10, 8, 5, 9, 6});
}
