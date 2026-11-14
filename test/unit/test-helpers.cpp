// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-helpers.hpp"

#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
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
  const std::string stored_dir = PALACE_TEST_DATA_DIR "/mesh/";
  static const std::map<std::string, std::string> simplicial_parent = {
      {"star-tri.mesh", "star-quad.mesh"}, {"fichera-tet.mesh", "fichera-hex.mesh"}};
  const auto it = simplicial_parent.find(name);
  if (it == simplicial_parent.end())
  {
    return stored_dir + name;
  }

  // Each process (CTest job or MPI rank) writes its own copy, so no two processes share a
  // file. MakeSimplicial is deterministic, so all ranks read the same mesh.
  static const palace::test::PerRankTempDir dir;
  const auto path = dir.temp_dir / name;
  if (!std::filesystem::exists(path))
  {
    mfem::Mesh parent(stored_dir + it->second);
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

namespace palace::test
{

std::vector<std::array<mfem::Vector, 3>> ObliqueFrames(int num_frames)
{
  std::mt19937 gen(877);
  std::normal_distribution<double> dist;
  std::vector<std::array<mfem::Vector, 3>> frames;
  while (static_cast<int>(frames.size()) < num_frames)
  {
    mfem::Vector n(3), a(3), e1(3), e2(3);
    for (int d = 0; d < 3; d++)
    {
      n(d) = dist(gen);
      a(d) = dist(gen);
    }
    n /= n.Norml2();
    a.Add(-(a * n), n);
    if (a.Norml2() < 1.0e-3 || std::abs(n(0) * n(1) * n(2)) < 1.0e-2)
    {
      continue;  // Keep the planes oblique to every axis
    }
    e1 = a;
    e1 /= e1.Norml2();
    e2(0) = n(1) * e1(2) - n(2) * e1(1);
    e2(1) = n(2) * e1(0) - n(0) * e1(2);
    e2(2) = n(0) * e1(1) - n(1) * e1(0);
    frames.push_back({e1, e2, n});
  }
  return frames;
}

}  // namespace palace::test
