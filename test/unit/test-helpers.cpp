// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-helpers.hpp"

#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include "fixtures.hpp"

mfem::Mesh SingleTetMesh()
{
  auto mesh = mfem::Mesh(3, 4, 1, 4);
  mesh.AddVertex(0.0, 0.0, 0.0);
  mesh.AddVertex(1.0, 0.0, 0.0);
  mesh.AddVertex(0.0, 1.0, 0.0);
  mesh.AddVertex(0.0, 0.0, 1.0);
  mesh.AddTet(0, 1, 2, 3);
  mesh.AddBdrTriangle(1, 2, 3, 1);
  mesh.AddBdrTriangle(0, 3, 2, 2);
  mesh.AddBdrTriangle(0, 1, 3, 3);
  mesh.AddBdrTriangle(0, 2, 1, 4);
  mesh.FinalizeTopology();
  mesh.Finalize(true, false);
  return mesh;
}

std::string TestMeshPath(const std::string &name)
{
  const std::string mfem_dir = PALACE_MFEM_DATA_DIR "/";
  static const std::map<std::string, std::string> simplicial_parent = {
      {"star-tri.mesh", "star.mesh"}, {"fichera-tet.mesh", "fichera.mesh"}};
  const auto it = simplicial_parent.find(name);
  if (it == simplicial_parent.end())
  {
    return mfem_dir + name;
  }

  // Each process (CTest job or MPI rank) writes its own copy, so no two processes share a
  // file. MakeSimplicial is deterministic, so all ranks read the same mesh.
  static const palace::test::PerRankTempDir dir;
  const auto path = dir.temp_dir / name;
  if (!std::filesystem::exists(path))
  {
    mfem::Mesh parent(mfem_dir + it->second);
    auto mesh = mfem::Mesh::MakeSimplicial(parent);
    std::ofstream out(path);
    // Full precision, so the file holds exactly the parent's vertex coordinates.
    out.precision(std::numeric_limits<double>::max_digits10);
    mesh.Print(out);
    if (!out)
    {
      throw std::runtime_error("Failed to write " + path.string());
    }
  }
  return path.string();
}
