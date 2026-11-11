// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fixtures.hpp"
#include "utils/gmshreader.hpp"

using namespace palace;

namespace
{

// Gmsh 2.2 mesh of the unit cube split into tetrahedra (Kuhn subdivision of 2x2x2 cubes),
// with its exterior triangles as boundary elements, plus physical points and lines (not
// part of the MFEM mesh), physical names, non-contiguous node tags, an unused node, and a
// node coordinate written with an exponent and a plus sign.
struct GmshTestMesh
{
  std::vector<std::array<double, 3>> nodes;
  std::vector<int> tags;
  struct Elem
  {
    int type, phys;
    std::vector<int> nodes;
  };
  std::vector<Elem> elems;

  GmshTestMesh()
  {
    const int n = 2;
    auto Id = [n](int i, int j, int k) { return (k * (n + 1) + j) * (n + 1) + i; };
    for (int k = 0; k <= n; k++)
    {
      for (int j = 0; j <= n; j++)
      {
        for (int i = 0; i <= n; i++)
        {
          nodes.push_back({0.5 * i, 0.5 * j, 0.5 * k});
          tags.push_back(3 * Id(i, j, k) + 7);  // Non-contiguous tags
        }
      }
    }
    nodes.push_back({2.0, 2.0, 2.0});  // Unused node
    tags.push_back(1000);
    // Kuhn subdivision: 6 tets per cube along the main diagonal.
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                             {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int k = 0; k < n; k++)
    {
      for (int j = 0; j < n; j++)
      {
        for (int i = 0; i < n; i++)
        {
          for (const auto &p : perms)
          {
            int x[3] = {i, j, k};
            std::vector<int> v = {tags[Id(x[0], x[1], x[2])]};
            for (int d = 0; d < 3; d++)
            {
              x[p[d]]++;
              v.push_back(tags[Id(x[0], x[1], x[2])]);
            }
            elems.push_back({4, 1 + (i + j + k) % 2, v});
          }
        }
      }
    }
    // Exterior triangles: faces of the tets with all vertices on one face of the cube.
    std::map<std::set<int>, int> face_count;
    std::map<std::set<int>, std::vector<int>> face_verts;
    const int tet_faces[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    for (const auto &e : elems)
    {
      for (const auto &f : tet_faces)
      {
        std::vector<int> fv = {e.nodes[f[0]], e.nodes[f[1]], e.nodes[f[2]]};
        std::set<int> key(fv.begin(), fv.end());
        face_count[key]++;
        face_verts[key] = fv;
      }
    }
    std::vector<Elem> bdr;
    for (const auto &[key, count] : face_count)
    {
      if (count == 1)
      {
        bdr.push_back({2, 10 + static_cast<int>(bdr.size() % 3), face_verts[key]});
      }
    }
    // Physical point and line first (other dimensions are ignored by the readers), then the
    // boundary elements, then the tetrahedra.
    std::vector<Elem> all = {{15, 20, {tags[0]}}, {1, 21, {tags[0], tags[1]}}};
    all.insert(all.end(), bdr.begin(), bdr.end());
    all.insert(all.end(), elems.begin(), elems.end());
    elems = all;
  }

  std::string Write(bool binary) const
  {
    std::ostringstream out;
    out << "$MeshFormat\n2.2 " << (binary ? 1 : 0) << " 8\n";
    if (binary)
    {
      const int one = 1;
      out.write(reinterpret_cast<const char *>(&one), sizeof(int));
      out << "\n";
    }
    out << "$EndMeshFormat\n";
    out << "$PhysicalNames\n3\n2 10 \"bdr_a\"\n2 11 \"bdr b\"\n3 1 \"vol_1\"\n"
           "$EndPhysicalNames\n";
    out << "$Nodes\n" << nodes.size() << "\n";
    out.precision(17);
    for (std::size_t i = 0; i < nodes.size(); i++)
    {
      if (binary)
      {
        out.write(reinterpret_cast<const char *>(&tags[i]), sizeof(int));
        out.write(reinterpret_cast<const char *>(nodes[i].data()), 3 * sizeof(double));
      }
      else if (i == 1)
      {
        out << tags[i] << " +5.0e-01 " << nodes[i][1] << " " << nodes[i][2] << "\n";
      }
      else
      {
        out << tags[i] << " " << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2]
            << "\n";
      }
    }
    if (binary)
    {
      out << "\n";
    }
    out << "$EndNodes\n$Elements\n" << elems.size() << "\n";
    for (std::size_t i = 0; i < elems.size(); i++)
    {
      const auto &e = elems[i];
      if (binary)
      {
        const int header[3] = {e.type, 1, 2};
        out.write(reinterpret_cast<const char *>(header), sizeof(header));
        const int num = static_cast<int>(i + 1), elem_tag = 5;
        out.write(reinterpret_cast<const char *>(&num), sizeof(int));
        out.write(reinterpret_cast<const char *>(&e.phys), sizeof(int));
        out.write(reinterpret_cast<const char *>(&elem_tag), sizeof(int));
        out.write(reinterpret_cast<const char *>(e.nodes.data()),
                  e.nodes.size() * sizeof(int));
      }
      else
      {
        out << i + 1 << " " << e.type << " 2 " << e.phys << " 5";
        for (int v : e.nodes)
        {
          out << " " << v;
        }
        out << "\n";
      }
    }
    if (binary)
    {
      out << "\n";
    }
    out << "$EndElements\n";
    return out.str();
  }
};

void WriteFile(const std::string &path, const std::string &data)
{
  std::ofstream fo(path, std::ios::binary);
  fo << data;
}

// The mesh read by MFEM, finalized like the mesh read by the fast reader for comparison.
std::unique_ptr<mfem::Mesh> ReadMfem(const std::string &path)
{
  struct UnfinalizedMesh : public mfem::Mesh
  {
    UnfinalizedMesh(std::istream &input) { Loader(input, 0); }
  };
  mfem::named_ifgzstream fi(path);
  return std::make_unique<mfem::Mesh>(UnfinalizedMesh(fi));
}

void CheckSameMesh(mfem::Mesh &a, mfem::Mesh &b)
{
  REQUIRE(a.Dimension() == b.Dimension());
  REQUIRE(a.SpaceDimension() == b.SpaceDimension());
  REQUIRE(a.GetNV() == b.GetNV());
  REQUIRE(a.GetNE() == b.GetNE());
  REQUIRE(a.GetNBE() == b.GetNBE());
  for (int v = 0; v < a.GetNV(); v++)
  {
    for (int d = 0; d < 3; d++)
    {
      CHECK(a.GetVertex(v)[d] == b.GetVertex(v)[d]);
    }
  }
  auto SameElement = [](const mfem::Element *x, const mfem::Element *y)
  {
    if (x->GetGeometryType() != y->GetGeometryType() ||
        x->GetAttribute() != y->GetAttribute())
    {
      return false;
    }
    return std::equal(x->GetVertices(), x->GetVertices() + x->GetNVertices(),
                      y->GetVertices());
  };
  for (int e = 0; e < a.GetNE(); e++)
  {
    CHECK(SameElement(a.GetElement(e), b.GetElement(e)));
  }
  for (int be = 0; be < a.GetNBE(); be++)
  {
    CHECK(SameElement(a.GetBdrElement(be), b.GetBdrElement(be)));
  }
  CHECK(a.GetNEdges() == b.GetNEdges());
  CHECK(a.GetNumFaces() == b.GetNumFaces());
  for (int be = 0; be < a.GetNBE(); be++)
  {
    CHECK(a.GetBdrElementFaceIndex(be) == b.GetBdrElementFaceIndex(be));
  }
  CHECK(a.bdr_attribute_sets.GetAttributeSetNames() ==
        b.bdr_attribute_sets.GetAttributeSetNames());
  CHECK(a.attribute_sets.GetAttributeSetNames() == b.attribute_sets.GetAttributeSetNames());
  for (const auto &name : a.bdr_attribute_sets.GetAttributeSetNames())
  {
    const auto &sa = a.bdr_attribute_sets.GetAttributeSet(name);
    const auto &sb = b.bdr_attribute_sets.GetAttributeSet(name);
    REQUIRE(sa.Size() == sb.Size());
    for (int i = 0; i < sa.Size(); i++)
    {
      CHECK(sa[i] == sb[i]);
    }
  }
}

}  // namespace

TEST_CASE_METHOD(palace::test::SharedTempDir, "Fast Gmsh 2.2 reader matches MFEM",
                 "[gmshreader][Serial]")
{
  GmshTestMesh gmsh;
  for (bool binary : {false, true})
  {
    DYNAMIC_SECTION((binary ? "binary" : "ASCII"))
    {
      const auto path =
          (temp_dir / (binary ? "cube_binary.msh" : "cube_ascii.msh")).string();
      WriteFile(path, gmsh.Write(binary));

      // With topology: same as mfem::Mesh::Loader.
      auto fast = mesh::ReadGmshMeshFast(path);
      REQUIRE(fast);
      auto ref = ReadMfem(path);
      CheckSameMesh(*fast, *ref);

      // Without topology, the mesh gets the same topology when finalized.
      auto fast_nt = mesh::ReadGmshMeshFast(path, false);
      REQUIRE(fast_nt);
      CHECK(fast_nt->GetNEdges() == 0);
      constexpr bool refine = true, fix_orientation = true;
      fast_nt->Finalize(refine, fix_orientation);
      ref->Finalize(refine, fix_orientation);
      CheckSameMesh(*fast_nt, *ref);
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Fast Gmsh reader falls back to MFEM",
                 "[gmshreader][Serial]")
{
  GmshTestMesh gmsh;
  const std::string ascii = gmsh.Write(false);
  auto Check = [this](const std::string &name, const std::string &data)
  {
    const auto path = (temp_dir / name).string();
    WriteFile(path, data);
    CHECK(!mesh::ReadGmshMeshFast(path));
  };
  SECTION("Format 4.1")
  {
    std::string data = ascii;
    data.replace(data.find("2.2 0 8"), 7, "4.1 0 8");
    Check("v41.msh", data);
  }
  SECTION("Periodic section")
  {
    std::string data = ascii + "$Periodic\n0\n$EndPeriodic\n";
    Check("periodic.msh", data);
  }
  SECTION("High-order element")
  {
    std::string data = ascii;
    const auto pos = data.rfind("\n", data.find("$EndElements") - 2);
    // Turn the last tetrahedron into a 10-node tetrahedron (type 11) with too few nodes:
    // the reader only needs to see the type to fall back.
    const auto line_start = pos + 1;
    const auto type_pos = data.find(' ', line_start) + 1;
    data.replace(type_pos, 1, "11");
    Check("ho.msh", data);
  }
  SECTION("Unknown node tag")
  {
    std::string data = ascii;
    const auto pos = data.rfind(' ', data.find("$EndElements") - 2);
    data.replace(pos + 1, data.find('\n', pos) - pos - 1, "99999");
    Check("tag.msh", data);
  }
  SECTION("Not a Gmsh file")
  {
    Check("mfem.mesh", "MFEM mesh v1.0\n");
  }
}
