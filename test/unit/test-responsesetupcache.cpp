// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fixtures.hpp"
#include "surfaceresponse-fixtures.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace
{

// Scope of PALACE_RESPONSE_DERIVED_CACHE=0 (the recompute-and-verify switch of the
// derived-matrix cache): set on construction, unset on destruction, so a failing assertion
// inside the scope cannot leak the switch into the following cases of the process.
class DerivedCacheRecomputeEnvGuard
{
public:
  DerivedCacheRecomputeEnvGuard() { setenv("PALACE_RESPONSE_DERIVED_CACHE", "0", 1); }
  DerivedCacheRecomputeEnvGuard(const DerivedCacheRecomputeEnvGuard &) = delete;
  DerivedCacheRecomputeEnvGuard &operator=(const DerivedCacheRecomputeEnvGuard &) = delete;
  ~DerivedCacheRecomputeEnvGuard() { unsetenv("PALACE_RESPONSE_DERIVED_CACHE"); }
};

// The operator applied to one deterministic vector: the bitwise comparator of two
// constructions (the cached derived matrices against the first computation's).
std::vector<double> ApplyToProbe(const SurfaceResponseOperator &response)
{
  Vector x(response.Height()), y(response.Height());
  for (int i = 0; i < x.Size(); i++)
  {
    x(i) = std::sin(0.37 * i + 0.11) + 0.25 * std::cos(1.3 * i);
  }
  y = 0.0;
  response.Mult(x, y);
  y.HostRead();
  return std::vector<double>(y.GetData(), y.GetData() + y.Size());
}

bool BitwiseEqual(const std::vector<double> &a, const std::vector<double> &b)
{
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

}  // namespace

// Decision 734 F1 (the response-setup PROFILE 7.1): BuildDomainResponseMatrices and
// BuildSurfaceResponseMatrices derive a model's unit-scaled, inverted and multiplied
// matrices from inputs that do not change across AMR cycles; their outputs are cached per
// process and keyed on the source files (path / size / time), the blend weights, the
// zero-trace indices, the interface mappings and the unit scale factors. On a private copy
// of the 2D automatic library (its own matrix files, so the shared fixture files are not
// touched): the first construction MISSES (one domain + one surface derivation), a second
// construction HITS both, the recompute-and-verify switch (PALACE_RESPONSE_DERIVED_CACHE=0)
// passes - the recomputation equals the cached matrices bit for bit - and a rewritten
// source file (its write time bumped) INVALIDATES exactly the derivation that read it. The
// operator's action on a probe vector is bitwise equal across all four constructions.
TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator derived-matrix cache across constructions",
                 "[surfaceresponseoperator][setupcache][2d][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // A private library: the fixture's model with every matrix file copied beside it.
  const fs::path private_dir = temp.temp_dir / "derived-cache";
  const fs::path private_fabricated = private_dir / "fabricated.csv";
  const fs::path private_thin = private_dir / "thin.csv";
  const fs::path private_fabricated_surface = private_dir / "fabricated-surface.csv";
  const fs::path private_thin_surface = private_dir / "thin-surface.csv";
  const fs::path private_library = private_dir / "fabrication-process-derived-cache.json";
  if (Mpi::Root(Mpi::World()))
  {
    fs::create_directories(private_dir);
    fs::copy_file(fabricated_path, private_fabricated,
                  fs::copy_options::overwrite_existing);
    fs::copy_file(thin_path, private_thin, fs::copy_options::overwrite_existing);
    fs::copy_file(fabricated_surface_path, private_fabricated_surface,
                  fs::copy_options::overwrite_existing);
    fs::copy_file(thin_surface_path, private_thin_surface,
                  fs::copy_options::overwrite_existing);
    std::ifstream input(library_path);
    REQUIRE(input);
    json library = json::parse(input);
    library["Name"] = "unit-test-process-derived-cache";
    REQUIRE(library["Models"].size() == 1);
    auto &model = library["Models"][0];
    model["FabricatedMatrix"] = private_fabricated.string();
    model["ThinMatrix"] = private_thin.string();
    model["FabricatedSurfaceMatrix"] = private_fabricated_surface.string();
    model["ThinSurfaceMatrix"] = private_thin_surface.string();
    std::ofstream output(private_library);
    output << library.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  json config = AutomaticConfig2D();
  config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
      private_library.string();
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  iodata.boundaries.cracked_attributes.insert(10);
  mfem::Mesh serial = MakeAutomatic2DMesh();
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(
      std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial)));
  LaplaceOperator laplace(iodata, meshes);
  std::shared_ptr<const SurfaceResponseGeometry> geometry;

  // MISS: the first construction derives the model's domain and surface matrices.
  const auto before_first = GetDerivedResponseMatrixCacheCounters();
  SurfaceResponseOperator first(iodata, laplace, &geometry);
  REQUIRE(geometry);
  REQUIRE(first.GetPatchCount() == 2);
  const auto after_first = GetDerivedResponseMatrixCacheCounters();
  CHECK(after_first.misses - before_first.misses == 2);
  CHECK(after_first.hits - before_first.hits == 0);
  CHECK(after_first.entries - before_first.entries == 2);
  const auto first_probe = ApplyToProbe(first);
  REQUIRE(!first_probe.empty());

  // HIT: the second construction (the AMR cycle-2 form: the same geometry passed back)
  // returns the cached derived matrices and acts identically.
  SurfaceResponseOperator second(iodata, laplace, &geometry);
  const auto after_second = GetDerivedResponseMatrixCacheCounters();
  CHECK(after_second.hits - after_first.hits == 2);
  CHECK(after_second.misses - after_first.misses == 0);
  CHECK(after_second.entries == after_first.entries);
  CHECK(second.GetPatchCount() == first.GetPatchCount());
  CHECK(BitwiseEqual(ApplyToProbe(second), first_probe));

  // The recompute-and-verify switch: every hit is recomputed and must equal the cached
  // matrices bit for bit (the proof of the bitwise claim on this library).
  {
    DerivedCacheRecomputeEnvGuard recompute;
    SurfaceResponseOperator verified(iodata, laplace, &geometry);
    const auto after_verified = GetDerivedResponseMatrixCacheCounters();
    CHECK(after_verified.hits - after_second.hits == 2);
    CHECK(after_verified.misses - after_second.misses == 0);
    CHECK(BitwiseEqual(ApplyToProbe(verified), first_probe));
  }

  // INVALIDATION: a rewritten domain source file (the same bytes, a later write time)
  // misses the domain derivation alone; the surface derivation, whose files are untouched,
  // still hits; the result is unchanged.
  if (Mpi::Root(Mpi::World()))
  {
    const auto time = fs::last_write_time(private_fabricated);
    fs::last_write_time(private_fabricated, time + std::chrono::seconds(2));
  }
  Mpi::Barrier(Mpi::World());
  const auto before_invalidated = GetDerivedResponseMatrixCacheCounters();
  SurfaceResponseOperator invalidated(iodata, laplace, &geometry);
  const auto after_invalidated = GetDerivedResponseMatrixCacheCounters();
  CHECK(after_invalidated.misses - before_invalidated.misses == 1);
  CHECK(after_invalidated.hits - before_invalidated.hits == 1);
  CHECK(after_invalidated.entries - before_invalidated.entries == 1);
  CHECK(BitwiseEqual(ApplyToProbe(invalidated), first_probe));
#endif
}

}  // namespace palace
