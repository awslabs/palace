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
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "surfaceresponse-fixtures.hpp"
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

// A field that is constant on each side of the plane, with the same tangential and a
// different normal component: the lowest-order Nedelec space holds it exactly on any
// tetrahedral or hexahedral mesh whose faces include the plane, cracked or not.
void PiecewiseConstantSheetField(const mfem::Vector &x, mfem::Vector &E)
{
  E.SetSize(3);
  E(0) = 0.6;
  E(1) = (x(1) > plane_y) ? 1.0 : -0.5;
  E(2) = -0.3;
}

using VectorFieldFunction = void (*)(const mfem::Vector &, mfem::Vector &);

std::unique_ptr<GridFunction> ProjectSheetField(FiniteElementSpace &nd_fespace,
                                                VectorFieldFunction field = SheetField)
{
  auto E = std::make_unique<GridFunction>(nd_fespace, false);
  mfem::VectorFunctionCoefficient f(3, field);
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

// The forms of an impedance-type metal sheet (Impedance, RationalImpedance, Conductivity)
// in the integration (USER 1566 (3); SD-000): the legacy key cracks it and splits an
// impedance between the two faces of the crack (twice the impedance on each face); the new
// default cracks it only with the per-boundary "Crack" true (two faces, the given impedance
// on each) and keeps it as one uncracked sheet with "Crack" false.
enum class MetalCrackForm
{
  LegacyKey,
  CrackTrue,
  CrackFalse
};

const char *MetalCrackFormName(MetalCrackForm form)
{
  switch (form)
  {
    case MetalCrackForm::LegacyKey:
      return "legacy-key";
    case MetalCrackForm::CrackTrue:
      return "crack-true";
    case MetalCrackForm::CrackFalse:
      return "crack-false";
  }
  return "";
}

// The boundary law (library "BoundaryCondition" form) of an impedance-type metal in a crack
// form, from its legacy-key law: the Crack-true form describes the same faces with twice
// the impedance (2 Rs, 2 Ls, Cs / 2; a rational impedance's numerator times two), as the
// legacy key splits it; a Conductivity is never split; the Crack-false form (one thin film)
// keeps the numbers.
json CrackFormLaw(const json &legacy_law, MetalCrackForm form)
{
  json law = legacy_law;
  if (form != MetalCrackForm::CrackTrue)
  {
    return law;
  }
  const auto type = law.at("Type").get<std::string>();
  if (type == "Impedance")
  {
    for (const auto *key : {"Rs", "Ls"})
    {
      if (law.contains(key))
      {
        law[key] = 2.0 * law[key].get<double>();
      }
    }
    if (law.contains("Cs"))
    {
      law["Cs"] = 0.5 * law["Cs"].get<double>();
    }
  }
  else if (type == "RationalImpedance")
  {
    for (auto &coefficient : law["Numerator"])
    {
      coefficient = 2.0 * coefficient.get<double>();
    }
  }
  return law;
}

// The placed patches of a response operator (the global list, held by the root) as
// order-free signatures: the model name, the origin and frame axes (to 1e-7 mesh units) and
// the weight (7 digits).
std::multiset<std::string> PatchSignatures(const SurfaceResponseOperator &response)
{
  const auto names = response.GetModelNames();
  std::multiset<std::string> signatures;
  for (const auto &patch : response.GetPatchAssignments())
  {
    std::string signature = names.at(patch.model);
    for (const auto *vector : {&patch.origin, &patch.axis_u, &patch.axis_v, &patch.axis_w})
    {
      for (const double x : *vector)
      {
        signature += "," + std::to_string(std::llround(x * 1.0e7));
      }
    }
    char weight[32];
    std::snprintf(weight, sizeof(weight), ",%.6e", patch.weight);
    signatures.insert(signature + weight);
  }
  return signatures;
}

// A Maxwell surface response problem read from its configuration: the mesh loaded with the
// crack form of the configuration by mesh::ReadMesh, the space operator and the response.
struct MaxwellIslandResponse
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> meshes;
  std::unique_ptr<SpaceOperator> space_op;
  std::unique_ptr<SurfaceResponseOperator> response;

  explicit MaxwellIslandResponse(const json &config) : iodata(config, false)
  {
    meshes.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
    space_op = std::make_unique<SpaceOperator>(iodata, meshes);
    response = std::make_unique<SurfaceResponseOperator>(iodata, *space_op);
  }

  auto ConstantFieldResponse(mfem::VectorCoefficient &field_coefficient) const
  {
    GridFunction field(space_op->GetNDSpace(), true);
    field.Real().ProjectCoefficient(field_coefficient);
    field.Imag() = 0.0;
    return response->GetMaxwellResponse(field, 0.0);
  }
};

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

TEST_CASE("SurfacePostOperator two-qualifying-side sheet energies cracked and uncracked",
          "[surfacepostoperator][sheetenergy][Serial][Parallel]")
{
  // P3-29 (R4 T3 #8, R2 A2): a metal sheet whose two sides both qualify for the interface
  // type (vacuum on both sides for MA, substrate on both sides for MS; any side for a
  // Default entry). Cracked (legacy key), each copy of the sheet has one side; uncracked
  // (default), upstream's sum-sides rule adds the energies of both sides on each face, and
  // that total is the two-sided sum, not the energy of the averaged field. On the
  // hexahedral mesh the two forms have the same elements (the island needs no bisection):
  // the raw totals, the coefficient path and our parts (outside, annulus, partitioned,
  // uncovered) agree for a smooth field. On the tetrahedral mesh the legacy crack bisects
  // two island corners (see the raw = parts case), so the totals are compared there for a
  // field that is constant on each side (exact on both meshes).
  const bool vacuum = GENERATE(true, false);
  const auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  const bool hexahedra = (type == mfem::Element::HEXAHEDRON);
  const std::string metal_type = vacuum ? "MA" : "MS";
  CAPTURE(metal_type, static_cast<int>(type));
  test::SharedTempDir temp;
  const auto mesh_path =
      WriteSheetMesh(temp.temp_dir, MakeSerialSheetMesh(type), "sheet.mesh");
  const json materials =
      vacuum ? json::array({{{"Attributes", {1, 2}}}})
             : json::array({{{"Attributes", {1, 2}}, {"Permittivity", 4.0}}});
  const json dielectric = json::array({DielectricEntry(1, island_attr, metal_type, false),
                                       DielectricEntry(2, island_attr, "Default", false),
                                       DielectricEntry(3, island_attr, metal_type, true),
                                       DielectricEntry(4, island_attr, "Default", true)});
  SheetProblem legacy(SheetConfig(mesh_path, temp.temp_dir, CrackPolicy::Legacy, materials,
                                  PecBoundaries(), dielectric));
  SheetProblem uncracked(SheetConfig(mesh_path, temp.temp_dir, CrackPolicy::Default,
                                     materials, PecBoundaries(), dielectric));
  REQUIRE(legacy.Cracked(island_attr));
  REQUIRE(!uncracked.Cracked(island_attr));
  if (hexahedra)
  {
    REQUIRE(legacy.mesh->Get().GetGlobalNE() == uncracked.mesh->Get().GetGlobalNE());
    REQUIRE(legacy.GlobalBdrElements(island_attr) ==
            2 * uncracked.GlobalBdrElements(island_attr));
  }
  else
  {
    REQUIRE(legacy.mesh->Get().GetGlobalNE() > uncracked.mesh->Get().GetGlobalNE());
  }

  SurfacePostOperator legacy_op(legacy.iodata, *legacy.mat_op, *legacy.h1_fespace,
                                *legacy.nd_fespace);
  SurfacePostOperator uncracked_op(uncracked.iodata, *uncracked.mat_op,
                                   *uncracked.h1_fespace, *uncracked.nd_fespace);
  // The same uncracked entries without the sheet list: the energy of the averaged field.
  SurfacePostOperator averaged_op(
      uncracked.iodata.boundaries.postpro, uncracked.iodata.problem.type, *uncracked.mat_op,
      *uncracked.h1_fespace, *uncracked.nd_fespace, {}, &uncracked.iodata.boundaries);
  const auto field = hexahedra ? SheetField : PiecewiseConstantSheetField;
  const auto E_legacy = ProjectSheetField(*legacy.nd_fespace, field);
  const auto E_uncracked = ProjectSheetField(*uncracked.nd_fespace, field);
  for (const int idx : {1, 2})
  {
    INFO("interface " << idx);
    const double raw_legacy = legacy_op.GetInterfaceElectricFieldEnergy(idx, *E_legacy);
    const double raw_uncracked =
        uncracked_op.GetInterfaceElectricFieldEnergy(idx, *E_uncracked);
    REQUIRE(raw_legacy > 0.0);
    CHECK_THAT(raw_uncracked, WithinRel(raw_legacy, 1.0e-10));
    CHECK_THAT(uncracked_op.GetInterfaceElectricFieldEnergy(idx + 2, *E_uncracked),
               WithinRel(raw_legacy, 1.0e-10));
    CHECK_THAT(legacy_op.GetInterfaceElectricFieldEnergy(idx + 2, *E_legacy),
               WithinRel(raw_legacy, 1.0e-10));
    const double averaged = averaged_op.GetInterfaceElectricFieldEnergy(idx, *E_uncracked);
    CAPTURE(raw_legacy, raw_uncracked, averaged);
    CHECK(std::abs(averaged - raw_uncracked) > 0.1 * raw_uncracked);
  }

  // Our parts of the MA / MS entry, side by side (where the two meshes have the same
  // elements).
  if (!hexahedra)
  {
    return;
  }
  const auto edges_legacy = legacy_op.GetInterfaceEdgeElectricFieldEnergies(1, *E_legacy);
  const auto edges_uncracked =
      uncracked_op.GetInterfaceEdgeElectricFieldEnergies(1, *E_uncracked);
  REQUIRE(edges_legacy.size() == 2);
  REQUIRE(edges_uncracked.size() == 2);
  for (std::size_t i = 0; i < 2; i++)
  {
    INFO("radius " << i);
    CHECK(edges_legacy[i].energy_annulus > 0.0);
    CHECK_THAT(edges_uncracked[i].energy_outside,
               WithinRel(edges_legacy[i].energy_outside, 1.0e-10));
    CHECK_THAT(edges_uncracked[i].energy_annulus,
               WithinRel(edges_legacy[i].energy_annulus, 1.0e-10));
  }
  const auto portions = IslandPortions();
  const auto partitioned_legacy =
      legacy_op.GetInterfacePartitionedEdgeEnergies({1}, *E_legacy, nullptr, portions)
          .at(1);
  const auto partitioned_uncracked =
      uncracked_op.GetInterfacePartitionedEdgeEnergies({1}, *E_uncracked, nullptr, portions)
          .at(1);
  REQUIRE(partitioned_legacy.size() == 4);
  REQUIRE(partitioned_uncracked.size() == partitioned_legacy.size());
  for (const auto &[label, energy] : partitioned_legacy)
  {
    INFO(label);
    CHECK(energy > 0.0);
    REQUIRE(partitioned_uncracked.count(label) == 1);
    CHECK_THAT(partitioned_uncracked.at(label), WithinRel(energy, 1.0e-10));
  }
  const std::vector<Portion> uncovered_portions = {portions[2]};
  CHECK_THAT(
      uncracked_op
          .GetInterfaceUncoveredEdgeEnergies({1}, *E_uncracked, nullptr, uncovered_portions)
          .at(1)
          .energy,
      WithinRel(legacy_op
                    .GetInterfaceUncoveredEdgeEnergies({1}, *E_legacy, nullptr,
                                                       uncovered_portions)
                    .at(1)
                    .energy,
                1.0e-10));
}

TEST_CASE("Batched interface response refuses a two-qualifying-side sheet",
          "[surfacepostoperator][sheetenergy][Serial]")
{
  // SD-401 (d): one amplitude vector per sample cannot hold the energies of two sides, so
  // the batched Gram and the streamed rows FAIL CLOSED on an uncracked sheet whose two
  // sides both qualify (MA with vacuum on both sides), while on the cracked mesh (legacy
  // key; one side per copy) they run and equal the raw energy. Serial: the refusal is
  // rank-local.
  test::SharedTempDir temp;
  const auto mesh_path = WriteSheetMesh(
      temp.temp_dir, MakeSerialSheetMesh(mfem::Element::HEXAHEDRON), "sheet.mesh");
  const json materials = json::array({{{"Attributes", {1, 2}}}});
  const json dielectric = json::array({DielectricEntry(1, island_attr, "MA", false, true)});
  for (const auto policy : {CrackPolicy::Legacy, CrackPolicy::Default})
  {
    INFO(CrackPolicyName(policy));
    SheetProblem problem(SheetConfig(mesh_path, temp.temp_dir, policy, materials,
                                     PecBoundaries(), dielectric));
    SurfacePostOperator surf_post_op(problem.iodata, *problem.mat_op, *problem.h1_fespace,
                                     *problem.nd_fespace);
    const auto E = ProjectSheetField(*problem.nd_fespace);
    const double raw = surf_post_op.GetInterfaceElectricFieldEnergy(1, *E);
    REQUIRE(raw > 0.0);
    const auto samples = surf_post_op.CacheInterfaceResponseSamples();
    REQUIRE(samples.size() == 1);
    std::vector<double> row(samples[0].RowSize());
    if (policy == CrackPolicy::Legacy)
    {
      const auto batched = surf_post_op.GetInterfaceElectricFieldEnergyMatrices({E.get()});
      CHECK_THAT(batched.at(1).at(0).energy_total(0, 0), WithinRel(raw, 1.0e-10));
      surf_post_op.EvaluateInterfaceResponseRow(samples[0], *E, nullptr, row.data());
      const auto streamed = surf_post_op.AssembleInterfaceResponseMatrices(
          samples[0], row.data(), 1, problem.mesh->GetComm());
      CHECK_THAT(streamed.at(0).energy_total(0, 0), WithinRel(raw, 1.0e-10));
    }
    else
    {
      constexpr const char *refusal =
          "Batched interface response fields cannot represent the summed energies of both "
          "sides of an interior sheet";
      CHECK_THROWS_WITH(surf_post_op.GetInterfaceElectricFieldEnergyMatrices({E.get()}),
                        ContainsSubstring(refusal));
      CHECK_THROWS_WITH(
          surf_post_op.EvaluateInterfaceResponseRow(samples[0], *E, nullptr, row.data()),
          ContainsSubstring(refusal));
    }
  }
}

TEST_CASE_METHOD(test::SurfaceResponseFiles,
                 "SurfaceResponseOperator identification on impedance-type metals in every "
                 "crack form",
                 "[surfaceresponseoperator][sheetenergy][maxwell][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // P3-30 (USER 1566 (3): per-boundary cracking of impedance-type metals; R4 T2 / T6): the
  // island metal of the Maxwell island problem as an Impedance, a RationalImpedance or a
  // Conductivity boundary, under the legacy key (cracked, the impedance split), the new
  // default with "Crack" true (cracked) and with "Crack" false (uncracked), loaded by
  // mesh::ReadMesh. The library is the convex PEC library with every model re-keyed to the
  // metal's own law: a test device (decision 1610), as a PEC-built library matches no model
  // on an impedance-type metal (OI-LIBRARY-KEYED-BY-METAL-BOUNDARY-LAW). In every form the
  // identification places the same patches as the PEC-equivalent (the island in Ground,
  // the PEC library) order-free, and the response runs with every edge matched (no refusal:
  // UnmatchedPolicy Error) and the law verified. The corrected values of a field normal to
  // the metal equal the PEC-equivalent's; those of an oblique field (whose tangential trace
  // on the metal enters an impedance-type island, not a PEC one) are the same in every type
  // and form. Whether PEC-built models are valid on impedance-type metals is not tested
  // here (OI-CORRECTION-VALIDATED-ON-PEC-METALS-ONLY).
  constexpr int metal_attr = 9;  // the island of MakeSerialIslandMesh
  const auto mesh_path =
      WriteSheetMesh(temp.temp_dir, MakeSerialIslandMesh(), "island.mesh");
  auto PecIslandConfig = [&](bool legacy_key)
  {
    json config = ConvexMaxwellIslandConfig();
    config["Model"]["Mesh"] = mesh_path.string();
    if (legacy_key)
    {
      config["Model"]["CrackInternalBoundaryElements"] = true;
    }
    else
    {
      config["Model"].erase("CrackInternalBoundaryElements");
    }
    return config;
  };
  MaxwellIslandResponse pec_cracked(PecIslandConfig(true));
  MaxwellIslandResponse pec_uncracked(PecIslandConfig(false));
  REQUIRE(pec_cracked.iodata.boundaries.cracked_attributes.contains(metal_attr));
  REQUIRE(!pec_uncracked.iodata.boundaries.cracked_attributes.contains(metal_attr));
  REQUIRE(pec_cracked.response->GetPatchCount() > 0);
  const auto pec_patches = PatchSignatures(*pec_cracked.response);
  CHECK(PatchSignatures(*pec_uncracked.response) == pec_patches);
  auto field_coefficient = ConstantFieldCoefficient();
  const auto pec_cracked_result = pec_cracked.ConstantFieldResponse(field_coefficient);
  const auto pec_uncracked_result = pec_uncracked.ConstantFieldResponse(field_coefficient);
  mfem::Vector normal_value(3);
  normal_value = 0.0;
  normal_value[1] = -0.4;
  mfem::VectorConstantCoefficient normal_coefficient(normal_value);
  const auto pec_cracked_normal = pec_cracked.ConstantFieldResponse(normal_coefficient);
  const auto pec_uncracked_normal = pec_uncracked.ConstantFieldResponse(normal_coefficient);
  // The response of the oblique field depends on the boundary type through its tangential
  // trace on the metal (essential on a PEC island only); that of a field normal to the
  // metal does not. Both are independent of the crack form.
  for (const auto &[cracked_result, uncracked_result] :
       {std::pair{&pec_cracked_result, &pec_uncracked_result},
        std::pair{&pec_cracked_normal, &pec_uncracked_normal}})
  {
    CHECK_THAT(uncracked_result->domain_correction,
               WithinRel(cracked_result->domain_correction, 1.0e-10));
    CHECK_THAT(uncracked_result->fabricated_surface_energy.at(4),
               WithinRel(cracked_result->fabricated_surface_energy.at(4), 1.0e-10));
  }
  std::optional<double> metal_oblique_correction, metal_oblique_surface_energy;

  for (const json *legacy_law :
       {&impedance_law, &rational_impedance_law, &conductivity_law})
  {
    for (const auto form :
         {MetalCrackForm::LegacyKey, MetalCrackForm::CrackTrue, MetalCrackForm::CrackFalse})
    {
      const json law = CrackFormLaw(*legacy_law, form);
      const auto type = law.at("Type").get<std::string>();
      const std::string label = type + "-" + MetalCrackFormName(form);
      INFO(label << " " << law.dump());
      const auto library_path = temp.temp_dir / ("convex-" + label + ".json");
      if (Mpi::Root(Mpi::World()))
      {
        std::ifstream input(convex_library_3d_path);
        json library = json::parse(input);
        library["Name"] = "unit-test-process-convex-3d-" + label;
        for (auto &model : library["Models"])
        {
          model["BoundaryCondition"] = law;
        }
        std::ofstream output(library_path);
        output << library.dump(2) << "\n";
      }
      Mpi::Barrier(Mpi::World());
      json config = PecIslandConfig(form == MetalCrackForm::LegacyKey);
      config["Boundaries"]["Ground"]["Attributes"] = {1, 2, 3, 4, 5, 6};
      json entry = law;
      entry.erase("Type");
      entry["Attributes"] = json::array({metal_attr});
      if (form != MetalCrackForm::LegacyKey)
      {
        entry["Crack"] = (form == MetalCrackForm::CrackTrue);
      }
      config["Boundaries"][type] = json::array({entry});
      config["Solver"]["SurfaceResponseCorrection"]["Library"] = library_path.string();
      MaxwellIslandResponse metal(config);

      // The mesh form, and the split of the impedance under the legacy key only.
      const bool cracked = (form != MetalCrackForm::CrackFalse);
      CHECK(metal.iodata.boundaries.cracked_attributes.contains(metal_attr) == cracked);
      CHECK(metal.iodata.boundaries.split_impedance_attributes.contains(metal_attr) ==
            (form == MetalCrackForm::LegacyKey));

      // The identification of the PEC-equivalent, order-free; every edge matched and the
      // law verified.
      CHECK(PatchSignatures(*metal.response) == pec_patches);
      const auto result = metal.ConstantFieldResponse(field_coefficient);
      CHECK_THAT(result.matched_length_fraction, WithinAbs(1.0, 1.0e-12));
      CHECK(result.boundary_law_verified);
      CHECK(result.loop_residual < 1.0e-10);

      // The corrected values: for the field normal to the metal those of the
      // PEC-equivalent; for the oblique field the same in every type and form.
      const auto normal = metal.ConstantFieldResponse(normal_coefficient);
      CHECK_THAT(normal.domain_correction,
                 WithinRel(pec_cracked_normal.domain_correction, 1.0e-10));
      CHECK_THAT(normal.fabricated_surface_energy.at(4),
                 WithinRel(pec_cracked_normal.fabricated_surface_energy.at(4), 1.0e-10));
      if (!metal_oblique_correction)
      {
        metal_oblique_correction = result.domain_correction;
        metal_oblique_surface_energy = result.fabricated_surface_energy.at(4);
      }
      CHECK(*metal_oblique_correction > 0.0);
      CHECK_THAT(result.domain_correction, WithinRel(*metal_oblique_correction, 1.0e-10));
      CHECK_THAT(result.fabricated_surface_energy.at(4),
                 WithinRel(*metal_oblique_surface_energy, 1.0e-10));
    }
  }
#endif
}

TEST_CASE_METHOD(
    test::SurfaceResponseFiles,
    "Driven reduced postprocessing is refused with a surface response correction",
    "[postoperator][romoperator][sheetenergy][Serial][Parallel]")
{
#if !defined(MFEM_USE_GSLIB)
  SKIP("SurfaceResponseOperator requires MFEM_USE_GSLIB");
#else
  // P3-32 (SD-505 (1), SD-606 (2), SD-701 (2)): upstream's reduced-coordinate output path
  // of the adaptive driven sweep (PostOperator::ConfigureReducedPostprocessing) is built on
  // the raw PROM's basis; a surface response correction (its gate, corrected energies and
  // corrected field) is measured on the full field only. With the correction (the
  // PostOperator holds a SurfaceResponseOperator) the reduced path stays disabled and the
  // driver's online loop takes MeasureAndPrintAll; without it (and without the interface
  // dielectric postprocessing, which upstream's guard also excludes) upstream's reduced
  // path is configured. The driven island problem: a lumped port on the exterior face z =
  // 0, the island metal in Ground, in both crack forms.
  const bool legacy_key = GENERATE(true, false);
  const bool correction = GENERATE(true, false);
  CAPTURE(legacy_key, correction);
  const auto mesh_path =
      WriteSheetMesh(temp.temp_dir, MakeSerialIslandMesh(), "driven-island.mesh");
  json config = ConvexMaxwellIslandConfig();
  config["Problem"]["Type"] = "Driven";
  config["Problem"]["Verbose"] = 0;
  config["Model"]["Mesh"] = mesh_path.string();
  if (legacy_key)
  {
    config["Model"]["CrackInternalBoundaryElements"] = true;
  }
  else
  {
    config["Model"].erase("CrackInternalBoundaryElements");
  }
  config["Boundaries"]["Ground"]["Attributes"] = {2, 3, 4, 5, 6, 9};
  config["Boundaries"]["LumpedPort"] = {{{"Index", 1},
                                         {"Attributes", {1}},
                                         {"R", 50.0},
                                         {"Direction", "+X"},
                                         {"Excitation", 1}}};
  json solver = {
      {"Order", 1},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 2.0}, {"MaxFreq", 32.0}, {"FreqStep", 30.0}}},
      {"Linear", {{"Type", "Default"}, {"KSPType", "GMRES"}, {"Tol", 1.0e-10}}}};
  if (correction)
  {
    solver["SurfaceResponseCorrection"] = config["Solver"]["SurfaceResponseCorrection"];
  }
  else
  {
    config["Boundaries"].erase("Postprocessing");
  }
  config["Solver"] = solver;
  IoData iodata(config, false);
  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(mesh::ReadMesh(iodata, Mpi::World())));
  CHECK(iodata.boundaries.cracked_attributes.contains(9) == legacy_key);
  SpaceOperator space_op(iodata, meshes);
  RomOperator rom_op(iodata, space_op, 2);
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  linalg::SetRandom(space_op.GetComm(), sample.Real(), 42);
  linalg::SetRandom(space_op.GetComm(), sample.Imag(), 43);
  rom_op.UpdatePROM(sample, "reduced_postprocessing_sample");
  REQUIRE(rom_op.GetReducedDimension() > 0);

  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);
  CHECK((post_op.GetSurfaceResponseOperator() != nullptr) == correction);
  post_op.ConfigureReducedPostprocessing(rom_op);
  CHECK(post_op.HasReducedPostprocessing() == !correction);
#endif
}

}  // namespace palace
