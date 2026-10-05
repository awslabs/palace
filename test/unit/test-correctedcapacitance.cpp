// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "driver.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace fs = std::filesystem;
using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

std::string ReadFile(const fs::path &path)
{
  std::ifstream input(path, std::ios::binary);
  REQUIRE(input);
  return std::string(std::istreambuf_iterator<char>(input), {});
}

Table LoadCsv(const fs::path &path)
{
  TableWithCSVFile wrapped(path.string(), /*load_existing_file=*/true);
  return std::move(wrapped.table);
}

// A loaded table names its columns col_<i>: look them up by header text.
const Column &ColumnByHeader(const Table &table, const std::string &header)
{
  for (auto it = table.cbegin(); it != table.cend(); ++it)
  {
    if (it->header_text == header)
    {
      return *it;
    }
  }
  FAIL("No column \"" << header << "\"");
  return *table.cbegin();
}

// A copy of a response-matrix file with its energy columns (from column `first_energy`,
// 0-based) scaled: the fixture matrices (1e-12 J) are 100x the unit-square device's energy
// under its ~19 nondimensional terminal potential and would drive the corrected energy
// negative (the driver's positivity check); at 1e-6 the correction is a few per cent.
void WriteScaledMatrix(const fs::path &source, const fs::path &target,
                       std::size_t first_energy, double scale)
{
  std::ifstream input(source);
  REQUIRE(input);
  std::ofstream output(target);
  std::string line;
  std::getline(input, line);
  output << line << "\n";
  while (std::getline(input, line))
  {
    if (line.empty())
    {
      continue;
    }
    std::vector<std::string> entries;
    std::size_t begin = 0;
    while (begin <= line.size())
    {
      const std::size_t end = line.find(',', begin);
      entries.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
      if (end == std::string::npos)
      {
        break;
      }
      begin = end + 1;
    }
    for (std::size_t i = 0; i < entries.size(); i++)
    {
      if (i >= first_energy)
      {
        output << fmt::format("{:.17e}", std::stod(entries[i]) * scale);
      }
      else
      {
        output << entries[i];
      }
      output << (i + 1 < entries.size() ? "," : "\n");
    }
  }
}

// The electrostatic driver end to end on an in-memory configuration (the unit-test-sized
// solve of the regression harness).
void RunElectrostatic(json config)
{
  IoData iodata(std::move(config), /*print=*/false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  const int omp_threads = utils::ConfigureOmp();
  BlockTimer::Reset();
  palace::Run(iodata, comm, omp_threads, /*git_tag=*/nullptr);
}

// A 2 x 2 matrix from the corrected terminal table's `<name> <variant>[i][j] <unit>`
// columns.
mfem::DenseMatrix CorrectedMatrix(const Table &table, const std::string &name,
                                  const std::string &variant, const std::string &unit)
{
  mfem::DenseMatrix mat(2);
  for (int j = 0; j < 2; j++)
  {
    const Column &col =
        ColumnByHeader(table, fmt::format("{} {}[i][{}] {}", name, variant, j + 1, unit));
    REQUIRE(col.data.size() == 2);
    for (int i = 0; i < 2; i++)
    {
      mat(i, j) = col.data[i];
    }
  }
  return mat;
}

}  // namespace

// Decision 334: the response-corrected capacitance of the fixed-trace DOMAIN defect,
// written next to the raw terminal tables without touching them. The 2D automatic library
// problem (two straight-edge patches on the metal line of MakeAutomatic2DMesh) with the
// line as a second terminal: the raw outputs of the run with ResponseCorrection are
// byte-identical to the run without it; terminal-C-corrected.csv (with the C_m and C⁻¹
// variants and the palace.json record) has the fixed-trace diagonal C_thin x E_ft / E_raw
// of surface-Q-corrected.csv, is symmetric, and its self-consistent column is NaN unless
// every source's corrected solve was accepted (PostprocessOnly: never; an unconverged
// corrected solve: fail closed).
TEST_CASE_METHOD(test::SurfaceResponseFiles, "Electrostatic corrected terminal capacitance",
                 "[electrostaticsolver][surfaceresponseoperator][2d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "corrected-capacitance-2d.mesh";
  const fs::path scaled_library_path =
      temp.temp_dir / "fabrication-process-corrected-capacitance.json";
  {
    mfem::Mesh serial = MakeAutomatic2DMesh();
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream output(mesh_path);
      serial.Print(output);
      constexpr double scale = 1.0e-6;
      std::ifstream library_input(library_path);
      REQUIRE(library_input);
      json library = json::parse(library_input);
      REQUIRE(library["Models"].size() == 1);
      auto &model = library["Models"][0];
      for (const auto &[key, first_energy] :
           {std::pair<const char *, std::size_t>{"FabricatedMatrix", 2},
            std::pair<const char *, std::size_t>{"ThinMatrix", 2},
            std::pair<const char *, std::size_t>{"FabricatedSurfaceMatrix", 5},
            std::pair<const char *, std::size_t>{"ThinSurfaceMatrix", 5}})
      {
        const fs::path source = model[key].get<std::string>();
        const fs::path target =
            temp.temp_dir / ("corrected-capacitance-" + source.filename().string());
        WriteScaledMatrix(source, target, first_energy, scale);
        model[key] = target.string();
      }
      std::ofstream library_output(scaled_library_path);
      library_output << library.dump(2) << "\n";
    }
  }
  Mpi::Barrier(Mpi::World());

  json config = AutomaticConfig2D();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      scaled_library_path.string();
  config["Boundaries"]["Ground"] = {{"Attributes", {1, 3, 4}}};
  config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {2}}},
                                      {{"Index", 2}, {"Attributes", {9, 10}}}};
  // The raw solve converges in 6 PCG iterations (Tol 1e-12); MaxIts bounds the corrected
  // solve of the unconverged run below.
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 12}};
  auto correction = config["Solver"]["Electrostatic"]["ResponseCorrection"];
  const fs::path raw_dir = temp.temp_dir / "corrected-capacitance-raw";
  const fs::path both_dir = temp.temp_dir / "corrected-capacitance-both";
  const fs::path postprocess_dir = temp.temp_dir / "corrected-capacitance-postprocess";
  const fs::path unconverged_dir = temp.temp_dir / "corrected-capacitance-unconverged";

  // The reference: no response correction.
  config["Problem"]["Output"] = raw_dir.string();
  config["Solver"]["Electrostatic"] = json::object();
  RunElectrostatic(config);
  // Postprocessing-only and self-consistent correction (the default CorrectionMode Both),
  // then postprocessing only.
  config["Problem"]["Output"] = both_dir.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"] = correction;
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] = "Both";
  RunElectrostatic(config);
  config["Problem"]["Output"] = postprocess_dir.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] =
      "PostprocessOnly";
  RunElectrostatic(config);
  // Both again with a corrected solve that cannot be accepted: 12 PCG iterations (the
  // recursive residual of this solve falls ~2.5 decades per iteration) cannot reach a
  // relative tolerance of 1e-300, so every source's corrected solve ends unconverged while
  // the raw solve (6 iterations) is unchanged.
  config["Problem"]["Output"] = unconverged_dir.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] = "Both";
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["SolveTol"] = 1.0e-300;
  RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  // Every raw output is byte-identical; the corrected tables exist only with a correction.
  for (const char *file : {"terminal-C.csv", "terminal-Cm.csv", "terminal-Cinv.csv",
                           "terminal-V.csv", "domain-E.csv", "surface-Q.csv"})
  {
    INFO(file);
    REQUIRE(fs::is_regular_file(raw_dir / file));
    CHECK(ReadFile(raw_dir / file) == ReadFile(both_dir / file));
    CHECK(ReadFile(raw_dir / file) == ReadFile(postprocess_dir / file));
    CHECK(ReadFile(raw_dir / file) == ReadFile(unconverged_dir / file));
  }
  for (const char *file : {"terminal-C-corrected.csv", "terminal-Cm-corrected.csv",
                           "terminal-Cinv-corrected.csv"})
  {
    INFO(file);
    CHECK_FALSE(fs::exists(raw_dir / file));
    CHECK(fs::is_regular_file(both_dir / file));
    CHECK(fs::is_regular_file(postprocess_dir / file));
    CHECK(fs::is_regular_file(unconverged_dir / file));
  }

  const Table raw_C = LoadCsv(raw_dir / "terminal-C.csv");
  const Table corrected_C = LoadCsv(both_dir / "terminal-C-corrected.csv");
  const Table corrected_Cm = LoadCsv(both_dir / "terminal-Cm-corrected.csv");
  const Table corrected_Cinv = LoadCsv(both_dir / "terminal-Cinv-corrected.csv");
  const Table energies = LoadCsv(both_dir / "surface-Q-corrected.csv");
  REQUIRE(corrected_C.n_cols() == 5);
  REQUIRE(corrected_C.n_rows() == 2);
  CHECK(corrected_C.cbegin()->header_text == "i");
  CHECK((corrected_C.cbegin() + 1)->header_text == "C fixed-trace[i][1] (F)");
  CHECK((corrected_C.cbegin() + 2)->header_text == "C fixed-trace[i][2] (F)");
  CHECK((corrected_C.cbegin() + 3)->header_text == "C corrected[i][1] (F)");
  CHECK((corrected_C.cbegin() + 4)->header_text == "C corrected[i][2] (F)");
  REQUIRE(energies.n_rows() == 2);
  const Column &E_raw = ColumnByHeader(energies, "E_elec raw (J)");
  const Column &E_ft = ColumnByHeader(energies, "E_elec postprocessed fixed-trace (J)");
  const Column &E_sc = ColumnByHeader(energies, "E_elec corrected (J)");
  auto RawC = [&](int i)
  { return ColumnByHeader(raw_C, fmt::format("C[i][{}] (F)", i + 1)).data[i]; };

  // Fixed trace: C_ft(i, i) = C_thin(i, i) x E_ft / E_raw (the raw bilinear form plus the
  // bilinear domain defect; the energies are printed to 12 digits), symmetric and
  // non-trivially corrected; C_m and C⁻¹ derive from it.
  const auto C_ft = CorrectedMatrix(corrected_C, "C", "fixed-trace", "(F)");
  const auto Cm_ft = CorrectedMatrix(corrected_Cm, "C_m", "fixed-trace", "(F)");
  const auto Cinv_ft = CorrectedMatrix(corrected_Cinv, "C⁻¹", "fixed-trace", "(1/F)");
  for (int i = 0; i < 2; i++)
  {
    const double C_raw = RawC(i);
    CHECK(std::isfinite(C_ft(i, i)));
    CHECK(C_ft(i, i) != C_raw);
    CHECK_THAT(C_ft(i, i), WithinRel(C_raw * E_ft.data[i] / E_raw.data[i], 1.0e-10));
    CHECK_THAT(Cm_ft(i, i), WithinRel(C_ft(i, 0) + C_ft(i, 1), 1.0e-12));
  }
  CHECK(C_ft(0, 1) != 0.0);
  CHECK_THAT(C_ft(0, 1), WithinRel(C_ft(1, 0), 1.0e-12));
  CHECK_THAT(Cm_ft(0, 1), WithinRel(-C_ft(0, 1), 1.0e-12));
  mfem::DenseMatrix identity(2);
  mfem::Mult(C_ft, Cinv_ft, identity);
  CHECK_THAT(identity(0, 0), WithinAbs(1.0, 1.0e-10));
  CHECK_THAT(identity(1, 1), WithinAbs(1.0, 1.0e-10));
  CHECK_THAT(identity(0, 1), WithinAbs(0.0, 1.0e-10));
  CHECK_THAT(identity(1, 0), WithinAbs(0.0, 1.0e-10));

  // Self-consistent: the same form on the corrected fields when every source's corrected
  // solve was accepted (then E_sc is finite for every source), else NaN. The scaled library
  // makes this run's corrected solves accepted (3 PCG iterations to 1e-6).
  const auto C_sc = CorrectedMatrix(corrected_C, "C", "corrected", "(F)");
  const bool accepted = std::isfinite(E_sc.data[0]) && std::isfinite(E_sc.data[1]);
  CHECK(accepted);
  for (int i = 0; i < 2; i++)
  {
    const double C_raw = RawC(i);
    if (accepted)
    {
      CHECK_THAT(C_sc(i, i), WithinRel(C_raw * E_sc.data[i] / E_raw.data[i], 1.0e-10));
    }
    else
    {
      CHECK(std::isnan(C_sc(i, i)));
    }
  }
  CHECK(std::isnan(C_sc(0, 1)) == !accepted);

  // PostprocessOnly: the fixed-trace matrix is the same (the same raw fields and form), the
  // self-consistent column is NaN (no corrected solve).
  const Table postprocess_C = LoadCsv(postprocess_dir / "terminal-C-corrected.csv");
  const auto postprocess_ft = CorrectedMatrix(postprocess_C, "C", "fixed-trace", "(F)");
  const auto postprocess_sc = CorrectedMatrix(postprocess_C, "C", "corrected", "(F)");
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(postprocess_ft(i, j) == C_ft(i, j));
      CHECK(std::isnan(postprocess_sc(i, j)));
    }
  }

  // The palace.json record.
  std::ifstream metadata_input(both_dir / "palace.json");
  REQUIRE(metadata_input);
  const auto metadata = json::parse(metadata_input);
  const auto &record = metadata.at("SurfaceResponse").at("TerminalCapacitance");
  CHECK(record.at("Indices") == json({1, 2}));
  REQUIRE(record.at("FixedTrace").size() == 2);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK_THAT(record.at("FixedTrace")[i][j].get<double>(),
                 WithinRel(C_ft(i, j), 1.0e-11));
    }
  }
  CHECK(record.at("SelfConsistent").is_null() == !accepted);
  std::ifstream postprocess_metadata_input(postprocess_dir / "palace.json");
  REQUIRE(postprocess_metadata_input);
  const auto postprocess_metadata = json::parse(postprocess_metadata_input);
  CHECK(postprocess_metadata.at("SurfaceResponse")
            .at("TerminalCapacitance")
            .at("SelfConsistent")
            .is_null());

  // Unconverged corrected solves (fail closed): the energies report the corrected energy as
  // unavailable, the self-consistent matrix is NaN and the palace.json variant null, while
  // the fixed-trace matrix is that of the same raw fields.
  const Table unconverged_energies = LoadCsv(unconverged_dir / "surface-Q-corrected.csv");
  const Column &unconverged_E_sc =
      ColumnByHeader(unconverged_energies, "E_elec corrected (J)");
  REQUIRE(unconverged_E_sc.data.size() == 2);
  CHECK(std::isnan(unconverged_E_sc.data[0]));
  CHECK(std::isnan(unconverged_E_sc.data[1]));
  const Table unconverged_C = LoadCsv(unconverged_dir / "terminal-C-corrected.csv");
  const auto unconverged_ft = CorrectedMatrix(unconverged_C, "C", "fixed-trace", "(F)");
  const auto unconverged_sc = CorrectedMatrix(unconverged_C, "C", "corrected", "(F)");
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      CHECK(unconverged_ft(i, j) == C_ft(i, j));
      CHECK(std::isnan(unconverged_sc(i, j)));
    }
  }
  std::ifstream unconverged_metadata_input(unconverged_dir / "palace.json");
  REQUIRE(unconverged_metadata_input);
  const auto unconverged_metadata = json::parse(unconverged_metadata_input);
  const auto &unconverged_record =
      unconverged_metadata.at("SurfaceResponse").at("TerminalCapacitance");
  CHECK(unconverged_record.at("FixedTrace") == record.at("FixedTrace"));
  CHECK(unconverged_record.at("SelfConsistent").is_null());
#endif
}

// Decision 352 follow-up (1): the per-patch energies of every applied patch, written on
// request only. The same 2D automatic library problem (two straight-edge patches of one
// model, two terminals) run with and without "PatchEnergy": every output of the run without
// is byte-identical with it (the model-energy, patch, corrected and raw tables); the run
// with it adds surface-response-patch-energy.csv whose rows per (source, evaluation, model)
// sum to the model-energy row (energies, patch count, patch weight) to roundoff, carry the
// patch provenance of surface-response-patches.csv (origin, weight) and the evaluation
// codes 0 / 1 / 2 of the model table; the Maxwell surface response refuses the option.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "Electrostatic surface-response patch energies",
                 "[electrostaticsolver][surfaceresponseoperator][2d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  const fs::path mesh_path = temp.temp_dir / "patch-energy-2d.mesh";
  const fs::path scaled_library_path =
      temp.temp_dir / "fabrication-process-patch-energy.json";
  {
    mfem::Mesh serial = MakeAutomatic2DMesh();
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream output(mesh_path);
      serial.Print(output);
      constexpr double scale = 1.0e-6;
      std::ifstream library_input(library_path);
      REQUIRE(library_input);
      json library = json::parse(library_input);
      REQUIRE(library["Models"].size() == 1);
      auto &model = library["Models"][0];
      for (const auto &[key, first_energy] :
           {std::pair<const char *, std::size_t>{"FabricatedMatrix", 2},
            std::pair<const char *, std::size_t>{"ThinMatrix", 2},
            std::pair<const char *, std::size_t>{"FabricatedSurfaceMatrix", 5},
            std::pair<const char *, std::size_t>{"ThinSurfaceMatrix", 5}})
      {
        const fs::path source = model[key].get<std::string>();
        const fs::path target =
            temp.temp_dir / ("patch-energy-" + source.filename().string());
        WriteScaledMatrix(source, target, first_energy, scale);
        model[key] = target.string();
      }
      std::ofstream library_output(scaled_library_path);
      library_output << library.dump(2) << "\n";
    }
  }
  Mpi::Barrier(Mpi::World());

  json config = AutomaticConfig2D();
  config["Model"]["Mesh"] = mesh_path.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      scaled_library_path.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] = "Both";
  config["Boundaries"]["Ground"] = {{"Attributes", {1, 3, 4}}};
  config["Boundaries"]["Terminal"] = {{{"Index", 1}, {"Attributes", {2}}},
                                      {{"Index", 2}, {"Attributes", {9, 10}}}};
  config["Solver"]["Linear"] = {{"Tol", 1.0e-12}, {"MaxIts", 12}};
  const fs::path without_dir = temp.temp_dir / "patch-energy-without";
  const fs::path with_dir = temp.temp_dir / "patch-energy-with";
  config["Problem"]["Output"] = without_dir.string();
  RunElectrostatic(config);
  config["Problem"]["Output"] = with_dir.string();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["PatchEnergy"] = true;
  RunElectrostatic(config);
  Mpi::Barrier(Mpi::World());

  // The Maxwell surface response has no per-patch export.
  {
    json maxwell = config;
    maxwell["Solver"]["SurfaceResponseCorrection"] =
        maxwell["Solver"]["Electrostatic"]["ResponseCorrection"];
    maxwell["Solver"].erase("Electrostatic");
    CHECK_THROWS_WITH(IoData(maxwell, /*print=*/false), ContainsSubstring("PatchEnergy"));
  }
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  for (const char *file :
       {"terminal-C.csv", "terminal-V.csv", "domain-E.csv", "surface-Q.csv",
        "surface-Q-corrected.csv", "terminal-C-corrected.csv",
        "surface-response-patches.csv", "surface-response-model-energy.csv"})
  {
    INFO(file);
    REQUIRE(fs::is_regular_file(without_dir / file));
    CHECK(ReadFile(without_dir / file) == ReadFile(with_dir / file));
  }
  CHECK_FALSE(fs::exists(without_dir / "surface-response-patch-energy.csv"));
  REQUIRE(fs::is_regular_file(with_dir / "surface-response-patch-energy.csv"));

  const Table models = LoadCsv(with_dir / "surface-response-model-energy.csv");
  const Table patches = LoadCsv(with_dir / "surface-response-patch-energy.csv");
  const Table assignments = LoadCsv(with_dir / "surface-response-patches.csv");
  const Column &model_source = ColumnByHeader(models, "source");
  const Column &model_evaluation = ColumnByHeader(models, "evaluation");
  const Column &model_model = ColumnByHeader(models, "model");
  const Column &model_count = ColumnByHeader(models, "patch count");
  const Column &model_weight = ColumnByHeader(models, "patch weight");
  const Column &model_domain = ColumnByHeader(models, "domain correction (J)");
  const Column &patch_source = ColumnByHeader(patches, "source");
  const Column &patch_evaluation = ColumnByHeader(patches, "evaluation");
  const Column &patch_patch = ColumnByHeader(patches, "patch");
  const Column &patch_model = ColumnByHeader(patches, "model");
  const Column &patch_feature = ColumnByHeader(patches, "feature");
  const Column &patch_weight = ColumnByHeader(patches, "weight");
  const Column &patch_domain = ColumnByHeader(patches, "domain correction (J)");
  const Column &patch_cell_begin = ColumnByHeader(patches, "cell begin (m)");
  const Column &patch_cell_end = ColumnByHeader(patches, "cell end (m)");
  // Two sources x three evaluations x one model; two patches each.
  REQUIRE(model_source.data.size() == 6);
  REQUIRE(patch_source.data.size() == 12);
  REQUIRE(assignments.n_rows() == 2);
  std::vector<std::string> interface_headers;
  for (auto it = models.cbegin(); it != models.cend(); ++it)
  {
    if (it->header_text.rfind("fabricated surface energy[", 0) == 0)
    {
      interface_headers.push_back(it->header_text);
    }
  }
  REQUIRE(!interface_headers.empty());
  auto Sum = [&](const Column &column, int source, int evaluation, int model)
  {
    double sum = 0.0;
    for (std::size_t r = 0; r < column.data.size(); r++)
    {
      if (patch_source.data[r] == source && patch_evaluation.data[r] == evaluation &&
          patch_model.data[r] == model)
      {
        sum += column.data[r];
      }
    }
    return sum;
  };
  auto Count = [&](int source, int evaluation, int model)
  {
    int count = 0;
    for (std::size_t r = 0; r < patch_source.data.size(); r++)
    {
      count += patch_source.data[r] == source && patch_evaluation.data[r] == evaluation &&
               patch_model.data[r] == model;
    }
    return count;
  };
  for (std::size_t r = 0; r < model_source.data.size(); r++)
  {
    const int source = static_cast<int>(model_source.data[r]);
    const int evaluation = static_cast<int>(model_evaluation.data[r]);
    const int model = static_cast<int>(model_model.data[r]);
    INFO("source " << source << " evaluation " << evaluation << " model " << model);
    CHECK(Count(source, evaluation, model) == static_cast<int>(model_count.data[r]));
    CHECK_THAT(Sum(patch_weight, source, evaluation, model),
               WithinRel(model_weight.data[r], 1.0e-12));
    REQUIRE(std::isfinite(model_domain.data[r]));
    CHECK(model_domain.data[r] != 0.0);
    CHECK_THAT(Sum(patch_domain, source, evaluation, model),
               WithinRel(model_domain.data[r], 1.0e-10));
    for (const auto &header : interface_headers)
    {
      const double total = ColumnByHeader(models, header).data[r];
      REQUIRE(std::isfinite(total));
      CHECK(total != 0.0);
      CHECK_THAT(Sum(ColumnByHeader(patches, header), source, evaluation, model),
                 WithinRel(total, 1.0e-10));
    }
  }
  // The patch provenance: 1-based patch indices 1 and 2 (the rows of the patch table), the
  // table's origin and weight, no feature (an explicit 2D library placement), a single
  // cross-section (cell begin = cell end = 0), and no two rows alike.
  const Column &assignment_weight = ColumnByHeader(assignments, "weight");
  for (std::size_t r = 0; r < patch_source.data.size(); r++)
  {
    const int patch = static_cast<int>(patch_patch.data[r]);
    REQUIRE((patch == 1 || patch == 2));
    CHECK(patch_feature.data[r] == -1.0);
    CHECK(patch_weight.data[r] == assignment_weight.data[patch - 1]);
    CHECK(patch_cell_begin.data[r] == 0.0);
    CHECK(patch_cell_end.data[r] == 0.0);
    for (const char *coordinate : {"x", "y", "z"})
    {
      const std::string header = fmt::format("origin {} (m)", coordinate);
      CHECK(ColumnByHeader(patches, header).data[r] ==
            ColumnByHeader(assignments, header).data[patch - 1]);
    }
    for (std::size_t q = 0; q < r; q++)
    {
      const bool same_row = patch_source.data[q] == patch_source.data[r] &&
                            patch_evaluation.data[q] == patch_evaluation.data[r] &&
                            patch_patch.data[q] == patch_patch.data[r];
      CHECK_FALSE(same_row);
    }
  }
  // The two patches differ (the two edges of the line see different fields).
  CHECK(patch_domain.data[0] != patch_domain.data[1]);
#endif
}

}  // namespace palace
