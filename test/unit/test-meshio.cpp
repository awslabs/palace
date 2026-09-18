// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
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
  nas << Field8("GRID") + Field8("2") + Field8("") + Field8("1.0D+0") + Field8("0.0") +
             Field8("0.0") + "\n";
  // Node 3: ordinary fixed-point fields, with the line truncated after the y field as
  // produced by writers that trim trailing blanks (the absent z field reads as zero).
  nas << Field8("GRID") + Field8("3") + Field8("") + Field8("0.0") + "1.0" + "\n";
  nas << Field8("CTRIA3") + Field8("1") + Field8("1") + Field8("1") + Field8("2") +
             Field8("3") + "\n";
  nas << "ENDDATA\n";  // Trailing newline: the reader requires a complete final line.

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
  CHECK(nodes[1][0] == Catch::Approx(1.0));    // "1.0D+0"
  CHECK(nodes[2][1] == Catch::Approx(1.0));    // Short line
  CHECK(nodes[2][2] == Catch::Approx(0.0));    // Absent (unpadded) field
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Nastran invalid number aborts",
                 "[meshio][Serial]")
{
  std::ostringstream nas;
  nas << "BEGIN BULK\n";
  nas << Field8("GRID") + Field8("1") + Field8("") + Field8("1.2.3") + Field8("0.0") +
             Field8("0.0") + "\n";
  nas << "ENDDATA\n";

  auto path = temp_dir / "invalid_float.nas";
  {
    std::ofstream f(path);
    f << nas.str();
  }
  std::stringstream buffer;
  CHECK_THROWS(mesh::ConvertMeshNastran(path.string(), buffer));
}
