// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "test-helpers.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <unistd.h>
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

std::string CaptureStdout(const std::function<void()> &f)
{
  std::fflush(stdout);
  std::FILE *file = std::tmpfile();
  const int saved = dup(fileno(stdout));
  dup2(fileno(file), fileno(stdout));
  auto restore = [&]()
  {
    std::fflush(stdout);
    dup2(saved, fileno(stdout));
    close(saved);
  };
  try
  {
    f();
  }
  catch (...)
  {
    restore();
    std::fclose(file);
    throw;
  }
  restore();
  std::rewind(file);
  std::string out;
  for (int c = std::fgetc(file); c != EOF; c = std::fgetc(file))
  {
    out.push_back(static_cast<char>(c));
  }
  std::fclose(file);
  return out;
}

mfem::Mesh MakeInteriorBoundaryBox()
{
  constexpr int nx = 4, n = 2;
  auto vid = [](int i, int j, int k) { return (i * (n + 1) + j) * (n + 1) + k; };
  mfem::Mesh mesh(3, (nx + 1) * (n + 1) * (n + 1), nx * n * n, 2 * n * n * 2 + 4 * nx * n,
                  3);
  for (int i = 0; i <= nx; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      for (int k = 0; k <= n; k++)
      {
        mesh.AddVertex(0.5 * i, 0.5 * j, 0.5 * k);
      }
    }
  }
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < n; j++)
    {
      for (int k = 0; k < n; k++)
      {
        const int v[8] = {vid(i, j, k),
                          vid(i + 1, j, k),
                          vid(i + 1, j + 1, k),
                          vid(i, j + 1, k),
                          vid(i, j, k + 1),
                          vid(i + 1, j, k + 1),
                          vid(i + 1, j + 1, k + 1),
                          vid(i, j + 1, k + 1)};
        mesh.AddHex(v, (i < 2) ? 1 : 2);
      }
    }
  }
  for (int a = 0; a < n; a++)
  {
    for (int b = 0; b < n; b++)
    {
      for (const auto [i, attr] : {std::pair{0, 1}, {2, 4}, {3, 5}, {nx, 2}})
      {
        const int q[4] = {vid(i, a, b), vid(i, a + 1, b), vid(i, a + 1, b + 1),
                          vid(i, a, b + 1)};
        mesh.AddBdrQuad(q, attr);
      }
    }
  }
  for (int i = 0; i < nx; i++)
  {
    for (int a = 0; a < n; a++)
    {
      for (int side = 0; side <= n; side += n)
      {
        const int qy[4] = {vid(i, side, a), vid(i + 1, side, a), vid(i + 1, side, a + 1),
                           vid(i, side, a + 1)};
        const int qz[4] = {vid(i, a, side), vid(i + 1, a, side), vid(i + 1, a + 1, side),
                           vid(i, a + 1, side)};
        mesh.AddBdrQuad(qy, 3);
        mesh.AddBdrQuad(qz, 3);
      }
    }
  }
  mesh.FinalizeHexMesh(1, 0, true);
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

void SetModelRun(nlohmann::json &config, std::string_view mode,
                 const std::string &model_path)
{
  const std::filesystem::path path(model_path);
  if (mode == "Online")
  {
    config["Solver"]["Substructuring"]["ModelFile"] = model_path;
    return;
  }
  if (path.filename() != "environment.model")
  {
    throw std::invalid_argument("An offline run saves its model as environment.model!");
  }
  config["Problem"]["Output"] = path.parent_path().string();
  config["Solver"]["Substructuring"]["SaveModel"] = true;
  std::error_code ec;
  std::filesystem::create_directories(path.parent_path(), ec);
}

std::string CaptureStdout(const std::function<void()> &f)
{
  struct Redirect
  {
    std::FILE *file = std::tmpfile();
    int saved = -1;
    Redirect()
    {
      if (!file)
      {
        throw std::runtime_error("Failed to open a temporary file for the standard output");
      }
      std::fflush(stdout);
      saved = dup(fileno(stdout));
      dup2(fileno(file), fileno(stdout));
    }
    ~Redirect()
    {
      std::fflush(stdout);
      dup2(saved, fileno(stdout));
      close(saved);
      std::fclose(file);
    }
  } redirect;
  f();
  std::fflush(stdout);
  std::rewind(redirect.file);
  std::string out;
  for (int c; (c = std::fgetc(redirect.file)) != EOF;)
  {
    out.push_back(static_cast<char>(c));
  }
  return out;
}

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
