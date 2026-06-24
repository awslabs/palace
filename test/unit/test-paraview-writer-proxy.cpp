// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fixtures.hpp"
#include "utils/ceedparaviewdatacollection.hpp"

namespace
{

namespace fs = std::filesystem;

int GetEnvInt(const char *name, int default_value)
{
  if (const char *value = std::getenv(name))
  {
    return std::atoi(value);
  }
  return default_value;
}

double GetEnvDouble(const char *name, double default_value)
{
  if (const char *value = std::getenv(name))
  {
    return std::atof(value);
  }
  return default_value;
}

fs::path GetOutputRoot(const palace::test::SharedTempDir &tmp)
{
  if (const char *value = std::getenv("PALACE_PROXY_OUTPUT_DIR"))
  {
    return fs::path(value);
  }
  return tmp.temp_dir / "transmon-nc-paraview-proxy";
}

template <typename Op>
double TimeMax(MPI_Comm comm, Op &&op)
{
  MPI_Barrier(comm);
  const auto t0 = std::chrono::steady_clock::now();
  op();
  MPI_Barrier(comm);
  const auto t1 = std::chrono::steady_clock::now();
  const double local = std::chrono::duration<double>(t1 - t0).count();
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, comm);
  return global;
}

std::uintmax_t DirectoryBytes(const fs::path &root)
{
  std::uintmax_t bytes = 0;
  if (!fs::exists(root))
  {
    return bytes;
  }
  for (const auto &entry : fs::recursive_directory_iterator(root))
  {
    if (entry.is_regular_file())
    {
      bytes += entry.file_size();
    }
  }
  return bytes;
}

void FillDeterministic(mfem::ParGridFunction &gf, int seed, int rank)
{
  // Keep this deterministic for a fixed MPI partition so the proxy can be used by
  // autoresearch without measuring random-data noise. The values do not need to be
  // physically meaningful; they just need to force the writer to handle ordinary
  // nonzero scalar/vector grid-function data.
  for (int i = 0; i < gf.Size(); i++)
  {
    const double x = static_cast<double>(i + 1 + 104729 * (rank + 1) + seed);
    gf(i) = std::sin(0.001 * x) + 0.25 * std::cos(0.00017 * x);
  }
}

void WriteGridFunctionFiles(const fs::path &dir, mfem::ParMesh &mesh,
                            mfem::ParGridFunction &scalar, mfem::ParGridFunction &vector)
{
  const int rank = mesh.GetMyRank();
  if (rank == 0)
  {
    fs::create_directories(dir);
  }
  MPI_Barrier(mesh.GetComm());

  auto write_one = [&](const mfem::ParGridFunction &gf, const std::string &name)
  {
    const fs::path path = dir / (name + ".gf." + std::to_string(rank));
    std::ofstream os(path);
    gf.Save(os);
  };

  write_one(scalar, "scalar");
  write_one(vector, "vector");

  // Match Palace's grid-function output enough to include mesh serialization cost.
  mesh.Save(dir / "mesh");
}

}  // namespace

// Hidden writer proxy for the large SingleTransmon AMR ParaView regression. This is not
// a correctness unit test for CI timing; run it explicitly, e.g.
//   palace-unit-tests "[paraview-writer-proxy]" --skip-benchmarks
// or with MPI/GPU through the Palace test environment. It avoids the eigensolve entirely:
// load the SingleTransmon mesh, create a deterministic nonconforming AMR mesh, fill
// high-order grid functions with synthetic data, and write both ParaView VTU and MFEM
// grid-function outputs. The printed METRIC lines are intended for autoresearch harnesses.
TEST_CASE_METHOD(palace::test::SharedTempDir, "SingleTransmon NC ParaView writer proxy",
                 "[.][paraview-writer-proxy][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank = 0, size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  const int order = GetEnvInt("PALACE_PROXY_ORDER", 3);
  const int refine_steps = GetEnvInt("PALACE_PROXY_REFINE_STEPS", 2);
  const double refine_prob = GetEnvDouble("PALACE_PROXY_REFINE_PROB", 0.3);
  const int seed = GetEnvInt("PALACE_PROXY_SEED", 20260624);
  const int compression_level = GetEnvInt("PALACE_PROXY_COMPRESSION", 1);

  REQUIRE(order >= 1);
  REQUIRE(refine_steps >= 0);
  REQUIRE(refine_prob >= 0.0);
  REQUIRE(refine_prob <= 1.0);

  const fs::path mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "examples" / "transmon" / "mesh" / "transmon.msh2";
  REQUIRE(fs::exists(mesh_path));

  fs::path output_root = GetOutputRoot(*this);
  if (rank == 0)
  {
    fs::remove_all(output_root);
    fs::create_directories(output_root);
  }
  MPI_Barrier(comm);

  mfem::Mesh serial_mesh(mesh_path.string().c_str(), 1, 1);
  REQUIRE(serial_mesh.Dimension() == 3);
  serial_mesh.EnsureNCMesh(true);

  // mfem::Mesh::RandomRefinement uses C rand(); seed it identically on every rank before
  // creating the parallel mesh so the serial pre-partitioned mesh is identical everywhere.
  std::srand(seed);
  for (int i = 0; i < refine_steps; i++)
  {
    serial_mesh.RandomRefinement(refine_prob, false, 1);
  }

  mfem::ParMesh mesh(comm, serial_mesh);
  serial_mesh.Clear();
  REQUIRE(mesh.Nonconforming());

  mfem::H1_FECollection fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace scalar_fes(&mesh, &fec);
  mfem::ParFiniteElementSpace vector_fes(&mesh, &fec, mesh.SpaceDimension(),
                                         mfem::Ordering::byVDIM);

  mfem::ParGridFunction scalar(&scalar_fes), vector(&vector_fes);
  FillDeterministic(scalar, seed, rank);
  FillDeterministic(vector, seed + 7919, rank);

  const fs::path paraview_dir = output_root / "paraview";
  const fs::path gridfunction_dir = output_root / "gridfunction";

  palace::CeedParaViewDataCollection pv("transmon_nc_proxy", &mesh);
  pv.SetPrefixPath(paraview_dir.string());
  pv.SetHighOrderOutput(true);
  pv.SetLevelsOfDetail(order);
  pv.SetDataFormat(mfem::VTKFormat::BINARY32);
  pv.SetCompressionLevel(compression_level);
  pv.SetCycle(0);
  pv.SetTime(0.0);
  pv.RegisterField("scalar", &scalar);
  pv.RegisterField("vector", &vector);

  const double paraview_seconds = TimeMax(comm, [&]() { pv.Save(); });
  const double gridfunction_seconds = TimeMax(
      comm, [&]() { WriteGridFunctionFiles(gridfunction_dir, mesh, scalar, vector); });

  std::uintmax_t paraview_bytes = 0;
  std::uintmax_t gridfunction_bytes = 0;
  if (rank == 0)
  {
    paraview_bytes = DirectoryBytes(paraview_dir);
    gridfunction_bytes = DirectoryBytes(gridfunction_dir);
  }

  long long local_elements = mesh.GetNE();
  long long total_elements = 0;
  MPI_Allreduce(&local_elements, &total_elements, 1, MPI_LONG_LONG, MPI_SUM, comm);

  long long local_scalar_dofs = scalar_fes.GetTrueVSize();
  long long local_vector_dofs = vector_fes.GetTrueVSize();
  long long total_scalar_dofs = 0;
  long long total_vector_dofs = 0;
  MPI_Allreduce(&local_scalar_dofs, &total_scalar_dofs, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Allreduce(&local_vector_dofs, &total_vector_dofs, 1, MPI_LONG_LONG, MPI_SUM, comm);

  if (rank == 0)
  {
    std::cout << "METRIC proxy_paraview_seconds=" << paraview_seconds << '\n';
    std::cout << "METRIC proxy_gridfunction_seconds=" << gridfunction_seconds << '\n';
    std::cout << "METRIC proxy_paraview_bytes=" << paraview_bytes << '\n';
    std::cout << "METRIC proxy_gridfunction_bytes=" << gridfunction_bytes << '\n';
    std::cout << "METRIC proxy_mpi_ranks=" << size << '\n';
    std::cout << "METRIC proxy_elements=" << total_elements << '\n';
    std::cout << "METRIC proxy_scalar_true_dofs=" << total_scalar_dofs << '\n';
    std::cout << "METRIC proxy_vector_true_dofs=" << total_vector_dofs << '\n';
    std::cout << "ARTIFACT proxy_output_dir=" << output_root << '\n';
  }

  CHECK(paraview_seconds >= 0.0);
  CHECK(gridfunction_seconds >= 0.0);
  if (rank == 0)
  {
    CHECK(paraview_bytes > 0);
    CHECK(gridfunction_bytes > 0);
  }
}
