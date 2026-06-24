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
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/surfacefunctional.hpp"
#include "models/materialoperator.hpp"
#include "utils/ceedparaviewdatacollection.hpp"
#include "utils/configfile.hpp"

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
                 "[.][paraview-raw-writer-proxy][Serial][Parallel][GPU]")
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


namespace
{

std::unique_ptr<palace::Mesh> MakeProxyPalaceMesh(MPI_Comm comm, const fs::path &mesh_path,
                                                  int refine_steps, double refine_prob,
                                                  int seed)
{
  mfem::Mesh serial_mesh(mesh_path.string().c_str(), 1, 1);
  REQUIRE(serial_mesh.Dimension() == 3);
  serial_mesh.EnsureNCMesh(true);

  std::srand(seed);
  for (int i = 0; i < refine_steps; i++)
  {
    serial_mesh.RandomRefinement(refine_prob, false, 1);
  }

  auto pmesh = std::make_unique<mfem::ParMesh>(comm, serial_mesh);
  serial_mesh.Clear();
  REQUIRE(pmesh->Nonconforming());
  return std::make_unique<palace::Mesh>(std::move(pmesh));
}

std::vector<int> DomainAttributes(const mfem::ParMesh &mesh)
{
  std::vector<int> attrs;
  attrs.reserve(mesh.attributes.Size());
  for (int i = 0; i < mesh.attributes.Size(); i++)
  {
    attrs.push_back(mesh.attributes[i]);
  }
  std::sort(attrs.begin(), attrs.end());
  attrs.erase(std::unique(attrs.begin(), attrs.end()), attrs.end());
  return attrs;
}

void FillDerivedProxyFields(palace::GridFunction &E, palace::GridFunction &B)
{
  mfem::VectorFunctionCoefficient er(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(0.17 * x(1)) + 0.03 * x(2) * x(2);
                                       v(1) = std::cos(0.11 * x(2)) + 0.07 * x(0);
                                       v(2) = 0.02 * x(0) * x(1) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient ei(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 0.05 * x(1) * x(2) - 0.5;
                                       v(1) = std::sin(0.13 * x(0)) - 0.04 * x(2);
                                       v(2) = std::cos(0.19 * x(1)) + 0.01 * x(0) * x(0);
                                     });
  mfem::VectorFunctionCoefficient br(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 0.03 * x(1) - 0.02 * x(2);
                                       v(1) = std::sin(0.07 * x(2)) + 0.5;
                                       v(2) = std::cos(0.09 * x(0)) - 0.01 * x(1) * x(2);
                                     });
  mfem::VectorFunctionCoefficient bi(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::cos(0.05 * x(2)) - 0.2;
                                       v(1) = 0.02 * x(0) * x(2) + 0.1;
                                       v(2) = std::sin(0.03 * x(1)) - 0.04 * x(0);
                                     });
  E.Real().ProjectCoefficient(er);
  E.Imag().ProjectCoefficient(ei);
  B.Real().ProjectCoefficient(br);
  B.Imag().ProjectCoefficient(bi);

  // Force the source fields into device-backed MFEM vectors. The auto2-like path below
  // asks libCEED to evaluate from these device grid functions into device output grid
  // functions, after which the ParaView writer must materialize host VTU arrays.
  E.Real().UseDevice(true);
  E.Imag().UseDevice(true);
  B.Real().UseDevice(true);
  B.Imag().UseDevice(true);
  (void)E.Real().ReadWrite(true);
  (void)E.Imag().ReadWrite(true);
  (void)B.Real().ReadWrite(true);
  (void)B.Imag().ReadWrite(true);
}

void ConfigureVTU(mfem::ParaViewDataCollection &pv, const fs::path &dir, int order,
                  int compression_level)
{
  pv.SetPrefixPath(dir.string());
  pv.SetHighOrderOutput(true);
  pv.SetLevelsOfDetail(order);
  pv.SetDataFormat(mfem::VTKFormat::BINARY32);
  pv.SetCompressionLevel(compression_level);
  pv.SetCycle(0);
  pv.SetTime(0.0);
}

void RegisterSourceFields(mfem::ParaViewDataCollection &pv, palace::GridFunction &E,
                          palace::GridFunction &B)
{
  pv.RegisterField("E_real", &E.Real());
  pv.RegisterField("E_imag", &E.Imag());
  pv.RegisterField("B_real", &B.Real());
  pv.RegisterField("B_imag", &B.Imag());
}

}  // namespace

// Hidden production-like proxy for the branch-specific SingleTransmon AMR ParaView
// regression. Unlike the raw writer microbenchmark above, this exercises the API that
// changed in -auto2: device-resident E/B grid functions are evaluated by libCEED at the
// ParaView/output points, first into materialized high-order grid functions (auto2-like),
// and also into direct VTU point buffers (stream branch). The legacy coefficient mode is
// kept in the same benchmark so a single run can expose whether the eager materialization
// path is the expensive part.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "SingleTransmon NC derived-field ParaView proxy",
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
  const bool include_source_fields =
      GetEnvInt("PALACE_PROXY_INCLUDE_SOURCE_FIELDS", 0) != 0;
  const bool include_direct_buffer = GetEnvInt("PALACE_PROXY_INCLUDE_DIRECT", 0) != 0;

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

  auto mesh = MakeProxyPalaceMesh(comm, mesh_path, refine_steps, refine_prob, seed);
  auto &pmesh = mesh->Get();

  palace::config::MaterialData material;
  material.attributes = DomainAttributes(pmesh);
  REQUIRE(!material.attributes.empty());
  palace::config::PeriodicBoundaryData periodic;
  palace::MaterialOperator mat_op({material}, periodic, palace::ProblemType::DRIVEN, *mesh);

  palace::fem::DefaultIntegrationOrder::p_trial = order;
  palace::fem::DefaultIntegrationOrder::q_order_jac = true;
  palace::fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  palace::fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  palace::FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  palace::FiniteElementSpace rt_fespace(*mesh, &rt_fec);
  palace::GridFunction E(nd_fespace, true), B(rt_fespace, true);
  FillDerivedProxyFields(E, B);

  mfem::L2_FECollection viz_fec(order, pmesh.Dimension());
  mfem::ParFiniteElementSpace viz_scalar(&pmesh, &viz_fec);
  mfem::ParFiniteElementSpace viz_vector(&pmesh, &viz_fec, pmesh.SpaceDimension());

  const double scaling_e = 1.0;
  const double scaling_m = 1.0;

  // Legacy origin/main-style coefficient registration: derived fields are evaluated by
  // MFEM coefficients as ParaView writes the high-order sampling lattice.
  palace::EnergyDensityCoefficient<palace::EnergyDensityType::ELECTRIC> Ue_coeff(
      E, mat_op, scaling_e);
  palace::EnergyDensityCoefficient<palace::EnergyDensityType::MAGNETIC> Um_coeff(
      B, mat_op, scaling_m);
  palace::PoyntingVectorCoefficient S_coeff(E, B, mat_op, scaling_m);

  const fs::path legacy_dir = output_root / "legacy_coeff";
  mfem::ParaViewDataCollection pv_legacy("transmon_nc_legacy_coeff", &pmesh);
  ConfigureVTU(pv_legacy, legacy_dir, order, compression_level);
  if (include_source_fields)
  {
    RegisterSourceFields(pv_legacy, E, B);
  }
  pv_legacy.RegisterCoeffField("U_e", &Ue_coeff);
  pv_legacy.RegisterCoeffField("U_m", &Um_coeff);
  pv_legacy.RegisterVCoeffField("S", &S_coeff);
  const double legacy_seconds = TimeMax(comm, [&]() { pv_legacy.Save(); });

  // -auto2-like path: libCEED evaluates the derived fields from device grid functions
  // into device-backed high-order output grid functions, which are then registered with
  // MFEM ParaView as ordinary fields. This is the path that exploded on large NC AMR.
  palace::DomainFieldEvaluator Ue_eval(palace::DomainFieldEvaluator::Kind::ENERGY_E,
                                       *mesh, mat_op, E.ParFESpace(), nullptr, viz_scalar,
                                       scaling_e);
  palace::DomainFieldEvaluator Um_eval(palace::DomainFieldEvaluator::Kind::ENERGY_M,
                                       *mesh, mat_op, nullptr, B.ParFESpace(), viz_scalar,
                                       scaling_m);
  palace::DomainFieldEvaluator S_eval(palace::DomainFieldEvaluator::Kind::POYNTING,
                                      *mesh, mat_op, E.ParFESpace(), B.ParFESpace(),
                                      viz_vector, scaling_m);
  REQUIRE(Ue_eval.IsValid());
  REQUIRE(Um_eval.IsValid());
  REQUIRE(S_eval.IsValid());

  mfem::ParGridFunction Ue_gf(&viz_scalar), Um_gf(&viz_scalar), S_gf(&viz_vector);
  Ue_gf.UseDevice(true);
  Um_gf.UseDevice(true);
  S_gf.UseDevice(true);
  const double materialize_seconds = TimeMax(comm,
                                             [&]()
                                             {
                                               Ue_eval.Eval(&E, nullptr, Ue_gf);
                                               Um_eval.Eval(nullptr, &B, Um_gf);
                                               S_eval.Eval(&E, &B, S_gf);
                                             });

  const fs::path materialized_dir = output_root / "auto2_materialized";
  mfem::ParaViewDataCollection pv_materialized("transmon_nc_auto2_materialized", &pmesh);
  ConfigureVTU(pv_materialized, materialized_dir, order, compression_level);
  if (include_source_fields)
  {
    RegisterSourceFields(pv_materialized, E, B);
  }
  pv_materialized.RegisterField("U_e", &Ue_gf);
  pv_materialized.RegisterField("U_m", &Um_gf);
  pv_materialized.RegisterField("S", &S_gf);
  const double materialized_save_seconds =
      TimeMax(comm, [&]() { pv_materialized.Save(); });
  const double auto2_like_seconds = materialize_seconds + materialized_save_seconds;

  // Optional stream-branch path: libCEED evaluates the same derived fields into VTU
  // point buffers ordered exactly as the writer consumes them, avoiding the high-order
  // grid-function materialization path. This is disabled by default because the purpose
  // of this proxy is to compare the origin/main coefficient path with the -auto2-like
  // materialized-grid-function path; enabling it is useful when testing the stream fix.
  double direct_buffer_eval_seconds = 0.0;
  double direct_buffer_save_seconds = 0.0;
  double direct_buffer_seconds = 0.0;
  std::uintmax_t direct_buffer_bytes = 0;
  if (include_direct_buffer)
  {
  palace::DomainFieldEvaluator Ue_buf_eval(palace::DomainFieldEvaluator::Kind::ENERGY_E,
                                           *mesh, mat_op, E.ParFESpace(), nullptr, order,
                                           scaling_e);
  palace::DomainFieldEvaluator Um_buf_eval(palace::DomainFieldEvaluator::Kind::ENERGY_M,
                                           *mesh, mat_op, nullptr, B.ParFESpace(), order,
                                           scaling_m);
  palace::DomainFieldEvaluator S_buf_eval(palace::DomainFieldEvaluator::Kind::POYNTING,
                                          *mesh, mat_op, E.ParFESpace(), B.ParFESpace(),
                                          order, scaling_m);
  REQUIRE(Ue_buf_eval.IsValid());
  REQUIRE(Um_buf_eval.IsValid());
  REQUIRE(S_buf_eval.IsValid());

  palace::Vector Ue_buf, Um_buf, S_buf;
  Ue_buf.SetSize(Ue_buf_eval.BufferSize());
  Um_buf.SetSize(Um_buf_eval.BufferSize());
  S_buf.SetSize(S_buf_eval.BufferSize());
  Ue_buf.UseDevice(true);
  Um_buf.UseDevice(true);
  S_buf.UseDevice(true);
  direct_buffer_eval_seconds = TimeMax(comm,
                                           [&]()
                                           {
                                             Ue_buf_eval.EvalBuffer(&E, nullptr, Ue_buf);
                                             Um_buf_eval.EvalBuffer(nullptr, &B, Um_buf);
                                             S_buf_eval.EvalBuffer(&E, &B, S_buf);
                                           });

  const fs::path direct_dir = output_root / "direct_point_buffer";
  palace::CeedParaViewDataCollection pv_direct("transmon_nc_direct_point_buffer", &pmesh);
  ConfigureVTU(pv_direct, direct_dir, order, compression_level);
  if (include_source_fields)
  {
    RegisterSourceFields(pv_direct, E, B);
  }
  pv_direct.RegisterDomainPointField(
      "U_e", Ue_buf, Ue_buf_eval.BufferBases(),
      palace::DomainFieldEvaluator::BufferNumComp(
          palace::DomainFieldEvaluator::Kind::ENERGY_E));
  pv_direct.RegisterDomainPointField(
      "U_m", Um_buf, Um_buf_eval.BufferBases(),
      palace::DomainFieldEvaluator::BufferNumComp(
          palace::DomainFieldEvaluator::Kind::ENERGY_M));
  pv_direct.RegisterDomainPointField(
      "S", S_buf, S_buf_eval.BufferBases(),
      palace::DomainFieldEvaluator::BufferNumComp(
          palace::DomainFieldEvaluator::Kind::POYNTING));
  direct_buffer_save_seconds = TimeMax(comm, [&]() { pv_direct.Save(); });
  direct_buffer_seconds = direct_buffer_eval_seconds + direct_buffer_save_seconds;
  if (rank == 0)
  {
    direct_buffer_bytes = DirectoryBytes(direct_dir);
  }
  }

  std::uintmax_t legacy_bytes = 0;
  std::uintmax_t materialized_bytes = 0;
  if (rank == 0)
  {
    legacy_bytes = DirectoryBytes(legacy_dir);
    materialized_bytes = DirectoryBytes(materialized_dir);
  }

  long long local_elements = pmesh.GetNE();
  long long total_elements = 0;
  MPI_Allreduce(&local_elements, &total_elements, 1, MPI_LONG_LONG, MPI_SUM, comm);
  long long local_e_true_dofs = nd_fespace.Get().GetTrueVSize();
  long long local_b_true_dofs = rt_fespace.Get().GetTrueVSize();
  long long local_scalar_true_dofs = viz_scalar.GetTrueVSize();
  long long local_vector_true_dofs = viz_vector.GetTrueVSize();
  long long total_e_true_dofs = 0, total_b_true_dofs = 0;
  long long total_scalar_true_dofs = 0, total_vector_true_dofs = 0;
  MPI_Allreduce(&local_e_true_dofs, &total_e_true_dofs, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Allreduce(&local_b_true_dofs, &total_b_true_dofs, 1, MPI_LONG_LONG, MPI_SUM, comm);
  MPI_Allreduce(&local_scalar_true_dofs, &total_scalar_true_dofs, 1, MPI_LONG_LONG,
                MPI_SUM, comm);
  MPI_Allreduce(&local_vector_true_dofs, &total_vector_true_dofs, 1, MPI_LONG_LONG,
                MPI_SUM, comm);

  if (rank == 0)
  {
    // Keep proxy_paraview_seconds as the primary autoresearch metric, but now point it at
    // the auto2-like production path rather than the raw writer microbenchmark.
    std::cout << "METRIC proxy_paraview_seconds=" << auto2_like_seconds << '\n';
    std::cout << "METRIC proxy_auto2_like_seconds=" << auto2_like_seconds << '\n';
    std::cout << "METRIC proxy_legacy_coefficient_seconds=" << legacy_seconds << '\n';
    std::cout << "METRIC proxy_materialize_seconds=" << materialize_seconds << '\n';
    std::cout << "METRIC proxy_materialized_save_seconds="
              << materialized_save_seconds << '\n';
    std::cout << "METRIC proxy_direct_buffer_seconds=" << direct_buffer_seconds << '\n';
    std::cout << "METRIC proxy_direct_buffer_eval_seconds="
              << direct_buffer_eval_seconds << '\n';
    std::cout << "METRIC proxy_direct_buffer_save_seconds="
              << direct_buffer_save_seconds << '\n';
    std::cout << "METRIC proxy_materialized_vs_legacy="
              << (auto2_like_seconds / legacy_seconds) << '\n';
    std::cout << "METRIC proxy_direct_vs_materialized="
              << (direct_buffer_seconds / auto2_like_seconds) << '\n';
    std::cout << "METRIC proxy_legacy_bytes=" << legacy_bytes << '\n';
    std::cout << "METRIC proxy_materialized_bytes=" << materialized_bytes << '\n';
    std::cout << "METRIC proxy_direct_buffer_bytes=" << direct_buffer_bytes << '\n';
    std::cout << "METRIC proxy_mpi_ranks=" << size << '\n';
    std::cout << "METRIC proxy_elements=" << total_elements << '\n';
    std::cout << "METRIC proxy_e_true_dofs=" << total_e_true_dofs << '\n';
    std::cout << "METRIC proxy_b_true_dofs=" << total_b_true_dofs << '\n';
    std::cout << "METRIC proxy_scalar_true_dofs=" << total_scalar_true_dofs << '\n';
    std::cout << "METRIC proxy_vector_true_dofs=" << total_vector_true_dofs << '\n';
    std::cout << "ARTIFACT proxy_output_dir=" << output_root << '\n';
  }

  CHECK(legacy_seconds >= 0.0);
  CHECK(materialize_seconds >= 0.0);
  CHECK(materialized_save_seconds >= 0.0);
  CHECK(direct_buffer_eval_seconds >= 0.0);
  CHECK(direct_buffer_save_seconds >= 0.0);
  if (rank == 0)
  {
    CHECK(legacy_bytes > 0);
    CHECK(materialized_bytes > 0);
    CHECK((!include_direct_buffer || direct_buffer_bytes > 0));
  }
}
