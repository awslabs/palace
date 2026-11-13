// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The thin-metal sheet energies of SurfacePostOperator on upstream's sheet plumbing (#994:
// sheet attributes, the two-sided interface sum, per-boundary "Crack") and its libCEED
// surface functional (#824), under both crack forms of the integration (SD-000): the legacy
// key (Model.CrackInternalBoundaryElements true: interior metal sheets cracked) and the new
// default (the key unset: PEC sheets uncracked).

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using json = nlohmann::json;
using namespace Catch::Matchers;
namespace fs = std::filesystem;

namespace
{

// The interior plane y = 0.5 of the unit box carries the metal island (x and z in
// [0.25, 0.75]) and, around it, the open dielectric interface; mfem::Mesh::MakeCartesian3D
// numbers the exterior boundary 1 to 6.
constexpr int island_attr = 9;
constexpr int open_attr = 10;
constexpr double plane_y = 0.5;
constexpr double island_min = 0.25, island_max = 0.75;

// The crack forms of the integration (SD-000): Legacy writes
// Model.CrackInternalBoundaryElements true (the bridge; every boundary condition sheet in
// the interior is cracked), Default leaves the key unset (upstream's default: PEC sheets
// are not cracked; impedance-type sheets only with "Crack": true).
enum class CrackPolicy
{
  Legacy,
  Default
};

const char *CrackPolicyName(CrackPolicy policy)
{
  return policy == CrackPolicy::Legacy ? "legacy" : "default";
}

// The unit box with n x n / 2 x n elements whose plane y = 0.5 is covered by interior
// boundary elements (the island and the open plane). The elements below the plane get
// attribute 1 and those above attribute 2.
mfem::Mesh MakeSerialSheetMesh(mfem::Element::Type type, int n = 8)
{
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(n, n / 2, n, type, 1.0, 1.0, 1.0);
  for (int face = 0; face < mesh.GetNumFaces(); face++)
  {
    int e1, e2;
    mesh.GetFaceElements(face, &e1, &e2);
    if (e1 < 0 || e2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    mesh.GetFaceVertices(face, vertices);
    bool on_plane = true, inside = true;
    for (const int v : vertices)
    {
      const double *x = mesh.GetVertex(v);
      on_plane = on_plane && std::abs(x[1] - plane_y) < 1.0e-12;
      inside = inside && x[0] > island_min - 1.0e-12 && x[0] < island_max + 1.0e-12 &&
               x[2] > island_min - 1.0e-12 && x[2] < island_max + 1.0e-12;
    }
    if (on_plane)
    {
      mesh.AddBdrElement(mesh.GetFace(face)->Duplicate(&mesh));
      mesh.SetBdrAttribute(mesh.GetNBE() - 1, inside ? island_attr : open_attr);
    }
  }
  mfem::Vector center(3);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, center);
    mesh.SetAttribute(e, center(1) < plane_y ? 1 : 2);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();
  return mesh;
}

// Write the serial mesh once (root) into the shared temporary directory.
fs::path WriteSheetMesh(const fs::path &dir, const mfem::Mesh &mesh,
                        const std::string &name)
{
  const auto path = dir / name;
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream output(path);
    mesh.Print(output);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

// An eigenmode configuration of the sheet mesh: `boundaries` holds the boundary conditions
// (the metal), the interface dielectric entries are added as its postprocessing, and the
// crack form sets or leaves unset the legacy key.
json SheetConfig(const fs::path &mesh_path, const fs::path &output, CrackPolicy policy,
                 const json &materials, json boundaries, const json &dielectric)
{
  boundaries["Postprocessing"] = {{"Dielectric", dielectric}};
  json config = {
      {"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}, {"Output", output.string()}}},
      {"Model", {{"Mesh", mesh_path.string()}}},
      {"Domains", {{"Materials", materials}}},
      {"Boundaries", boundaries},
      {"Solver", {{"Order", 1}, {"Eigenmode", {{"Target", 1.0}}}}}};
  if (policy == CrackPolicy::Legacy)
  {
    config["Model"]["CrackInternalBoundaryElements"] = true;
  }
  return config;
}

// Substrate (permittivity 4) below the plane, vacuum above it.
json SubstrateVacuumMaterials()
{
  return json::array({{{"Attributes", {1}}, {"Permittivity", 4.0}}, {{"Attributes", {2}}}});
}

// The PEC metal: the exterior boundary and the island.
json PecBoundaries()
{
  return {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6, island_attr}}}}};
}

// An interface dielectric entry on one attribute: with automatic metal edges at the
// matching radii R = 0.1 and 2R (the edge / partitioned / uncovered energies; with
// `localize` also the batched and streamed response evaluators), or, with `region`, inside
// a box enclosing the whole mesh (the same integrand on the coefficient path: an entry with
// a Region never takes the libCEED surface functional).
json DielectricEntry(int index, int attr, const std::string &type, bool region,
                     bool localize = false)
{
  json entry = {{"Index", index},
                {"Attributes", {attr}},
                {"Type", type},
                {"Thickness", 2.0e-3},
                {"Permittivity", 10.0}};
  if (region)
  {
    entry["Region"] = {{"BoxMin", {-10.0, -10.0, -10.0}}, {"BoxMax", {10.0, 10.0, 10.0}}};
  }
  else if (type != "Default")
  {
    entry["AutomaticEdges"] = true;
    entry["EdgeDistances"] = {0.1, 0.2};
    entry["LocalizeEdgeEnergy"] = localize;
  }
  return entry;
}

// The configuration, mesh (loaded with the crack form of the configuration by
// mesh::ReadMesh), materials and first-order spaces of one sheet problem.
struct SheetProblem
{
  IoData iodata;
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<MaterialOperator> mat_op;
  mfem::H1_FECollection h1_fec;
  mfem::ND_FECollection nd_fec;
  std::unique_ptr<FiniteElementSpace> h1_fespace, nd_fespace;

  explicit SheetProblem(const json &config)
    : iodata(config, false), h1_fec(1, 3), nd_fec(1, 3)
  {
    mesh = std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World()));
    mat_op = std::make_unique<MaterialOperator>(iodata, *mesh);
    h1_fespace = std::make_unique<FiniteElementSpace>(*mesh, &h1_fec);
    nd_fespace = std::make_unique<FiniteElementSpace>(*mesh, &nd_fec);
  }

  bool Cracked(int attr) const
  {
    return iodata.boundaries.cracked_attributes.contains(attr);
  }

  // The number of island boundary elements, summed over the ranks (twice the uncracked
  // count on a cracked mesh: one copy per side).
  int GlobalBdrElements(int attr) const
  {
    int count = 0;
    const auto &pmesh = mesh->Get();
    for (int be = 0; be < pmesh.GetNBE(); be++)
    {
      count += (pmesh.GetBdrAttribute(be) == attr);
    }
    Mpi::GlobalSum(1, &count, pmesh.GetComm());
    return count;
  }
};

// A field with a normal component that jumps across the plane (sign and size) and smooth
// tangential components (continuous across the plane, so that a cracked and an uncracked
// mesh project the same field side by side).
void SheetField(const mfem::Vector &x, mfem::Vector &E)
{
  E.SetSize(3);
  const double side = (x(1) > plane_y) ? 1.0 : -0.5;
  E(0) = 0.7 + 0.2 * x(2) + 0.1 * x(1);
  E(1) = side * (1.0 + 0.3 * x(0) - 0.2 * x(2));
  E(2) = -0.4 + 0.1 * x(0) + 0.3 * x(1);
}

std::unique_ptr<GridFunction> ProjectSheetField(FiniteElementSpace &nd_fespace)
{
  auto E = std::make_unique<GridFunction>(nd_fespace, false);
  mfem::VectorFunctionCoefficient f(3, SheetField);
  E->Real().ProjectCoefficient(f);
  Vector e(nd_fespace.GetTrueVSize());
  E->Real().GetTrueDofs(e);
  E->Real().SetFromTrueDofs(e);
  E->Real().ExchangeFaceNbrData();
  return E;
}

// Perimeter portions of the island (on its plane, in mesh coordinates): the edge z = 0.25
// labelled "m0", the edge x = 0.25 labelled "m1" and the edge z = 0.75 labelled "u" (the
// edge x = 0.75 stays unclaimed).
using Portion = SurfacePostOperator::UncoveredPerimeterPortion;
std::vector<Portion> IslandPortions()
{
  return {{{island_min, plane_y, island_min}, {island_max, plane_y, island_min}, "m0", 0},
          {{island_min, plane_y, island_min}, {island_min, plane_y, island_max}, "m1", 1},
          {{island_min, plane_y, island_max}, {island_max, plane_y, island_max}, "u", 2}};
}

double Sum(const std::map<std::string, double> &energies)
{
  double sum = 0.0;
  for (const auto &[label, energy] : energies)
  {
    (void)label;
    sum += energy;
  }
  return sum;
}

}  // namespace

TEST_CASE("SurfacePostOperator raw interface energy equals the sum of its parts",
          "[surfacepostoperator][sheetenergy][Serial][Parallel]")
{
  // P3-28 (SL6-RESOLUTION s3, SD-604): the raw total of a plain entry (upstream's libCEED
  // surface functional with the entry's sum-sides marker) equals the coefficient path of
  // the same entry (inside an enclosing Region), and our parts on the coefficient path
  // partition it: outside(R) = annulus(R) + outside(2R); the partitioned within-2R energies
  // (the portion labels and "unclaimed") = raw - outside(2R); the uncovered energy of a
  // portion = its partitioned share; the weighted portion energies at weight 1 = the
  // partitioned labels; the batched and the streamed response evaluators (one basis
  // field) = raw in total and = raw - outside within each radius. For MA and MS on the
  // island (substrate below, vacuum above) and SA on the open plane around it, on a
  // cracked (legacy key) and an uncracked (default) mesh.
  const auto policy = GENERATE(CrackPolicy::Legacy, CrackPolicy::Default);
  const auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  CAPTURE(CrackPolicyName(policy), static_cast<int>(type));
  test::SharedTempDir temp;
  const auto mesh_path =
      WriteSheetMesh(temp.temp_dir, MakeSerialSheetMesh(type), "sheet.mesh");
  const json dielectric = json::array({DielectricEntry(1, island_attr, "MA", false, true),
                                       DielectricEntry(2, island_attr, "MS", false, true),
                                       DielectricEntry(3, open_attr, "SA", false, true),
                                       DielectricEntry(4, island_attr, "MA", true),
                                       DielectricEntry(5, island_attr, "MS", true),
                                       DielectricEntry(6, open_attr, "SA", true)});
  SheetProblem problem(SheetConfig(mesh_path, temp.temp_dir, policy,
                                   SubstrateVacuumMaterials(), PecBoundaries(),
                                   dielectric));

  // The form is the one named: the island is cracked (one boundary element per side) under
  // the legacy key only; it is a sheet attribute in both forms, the open plane in neither.
  // On the tetrahedral mesh the legacy crack also bisects the two island corners whose
  // triangles join two perimeter vertices (RefineCrackElements: an edge interior to the
  // crack between two seam vertices), so there it has more copies than uncracked faces.
  const int island_faces = (type == mfem::Element::HEXAHEDRON) ? 16 : 32;
  const bool legacy = (policy == CrackPolicy::Legacy);
  CHECK(problem.Cracked(island_attr) == legacy);
  CHECK(!problem.Cracked(open_attr));
  const int island_bdr_elements = problem.GlobalBdrElements(island_attr);
  if (!legacy)
  {
    CHECK(island_bdr_elements == island_faces);
  }
  else if (type == mfem::Element::HEXAHEDRON)
  {
    CHECK(island_bdr_elements == 2 * island_faces);
  }
  else
  {
    CHECK(island_bdr_elements % 2 == 0);
    CHECK(island_bdr_elements > 2 * island_faces);
  }
  const auto sheet_attributes = problem.iodata.boundaries.GetSheetAttributes();
  CHECK(std::ranges::find(sheet_attributes, island_attr) != sheet_attributes.end());
  CHECK(std::ranges::find(sheet_attributes, open_attr) == sheet_attributes.end());

  SurfacePostOperator surf_post_op(problem.iodata, *problem.mat_op, *problem.h1_fespace,
                                   *problem.nd_fespace);
  const auto E = ProjectSheetField(*problem.nd_fespace);
  const auto portions = IslandPortions();
  const std::vector<Portion> uncovered_portions = {portions[2]};
  const auto batched = surf_post_op.GetInterfaceElectricFieldEnergyMatrices({E.get()});
  std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>> streamed;
  for (const auto &samples : surf_post_op.CacheInterfaceResponseSamples())
  {
    std::vector<double> row(samples.RowSize());
    surf_post_op.EvaluateInterfaceResponseRow(samples, *E, nullptr, row.data());
    streamed.emplace(samples.interface_index,
                     surf_post_op.AssembleInterfaceResponseMatrices(
                         samples, row.data(), 1, problem.mesh->GetComm()));
  }
  REQUIRE(batched.size() == 3);
  REQUIRE(streamed.size() == 3);
  for (const int idx : {1, 2, 3})
  {
    INFO("interface " << idx);
    const double raw = surf_post_op.GetInterfaceElectricFieldEnergy(idx, *E);
    const double coefficient_path =
        surf_post_op.GetInterfaceElectricFieldEnergy(idx + 3, *E);
    REQUIRE(raw > 0.0);
    CHECK_THAT(coefficient_path, WithinRel(raw, 1.0e-10));

    // Outside and annulus energies at R and 2R; the outer energy is the one at 2R.
    const auto edges = surf_post_op.GetInterfaceEdgeElectricFieldEnergies(idx, *E);
    REQUIRE(edges.size() == 2);
    CHECK(edges[0].energy_annulus > 0.0);
    CHECK_THAT(edges[0].energy_outside,
               WithinRel(edges[0].energy_annulus + edges[1].energy_outside, 1.0e-10));
    const auto outer = surf_post_op.GetInterfaceOuterElectricFieldEnergy(idx, *E);
    CHECK_THAT(outer.energy_outside, WithinRel(edges[1].energy_outside, 1.0e-12));
    const double within = raw - edges[1].energy_outside;
    REQUIRE(within > 0.0);

    // The exact partition of the within-2R energy over the portions' labels.
    const auto partitioned =
        surf_post_op.GetInterfacePartitionedEdgeEnergies({idx}, *E, nullptr, portions)
            .at(idx);
    REQUIRE(partitioned.size() == 4);
    for (const auto &label : {"m0", "m1", "u", "unclaimed"})
    {
      INFO(label);
      REQUIRE(partitioned.count(label) == 1);
      CHECK(partitioned.at(label) > 0.0);
    }
    CHECK_THAT(Sum(partitioned), WithinRel(within, 1.0e-9));
    CHECK_THAT(edges[1].energy_outside + Sum(partitioned), WithinRel(raw, 1.0e-10));

    // The uncovered energy of a portion is its partitioned share.
    const auto uncovered =
        surf_post_op
            .GetInterfaceUncoveredEdgeEnergies({idx}, *E, nullptr, uncovered_portions)
            .at(idx);
    REQUIRE(uncovered.by_type.size() == 1);
    CHECK_THAT(uncovered.energy, WithinRel(partitioned.at("u"), 1.0e-9));
    CHECK_THAT(uncovered.by_type.at("u"), WithinRel(uncovered.energy, 1.0e-12));

    // The weighted portion energies at weight 1 are the partitioned labels.
    const auto weighted = surf_post_op
                              .GetInterfaceWeightedPortionEdgeEnergies(
                                  {idx}, *E, nullptr, portions, {1.0, 1.0, 1.0})
                              .at(idx);
    REQUIRE(weighted.size() == 3);
    for (const auto &[label, energy] : weighted)
    {
      INFO(label);
      CHECK_THAT(energy, WithinRel(partitioned.at(label), 1.0e-9));
    }

    // The batched Gram and the streamed rows of the response evaluators.
    using ResponseMatrices = std::vector<SurfacePostOperator::InterfaceResponseMatrix>;
    for (const ResponseMatrices *matrices :
         std::array<const ResponseMatrices *, 2>{&batched.at(idx), &streamed.at(idx)})
    {
      REQUIRE(matrices->size() == 2);
      for (std::size_t i = 0; i < 2; i++)
      {
        INFO("radius " << i);
        const auto &matrix = (*matrices)[i];
        REQUIRE(matrix.energy_total.Height() == 1);
        CHECK_THAT(matrix.energy_total(0, 0), WithinRel(raw, 1.0e-10));
        CHECK_THAT(matrix.energy_inside(0, 0),
                   WithinRel(raw - edges[i].energy_outside, 1.0e-9));
      }
    }
  }
}

}  // namespace palace
