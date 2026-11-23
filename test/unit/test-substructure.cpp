// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/interpolator.hpp"
#include "fem/substructure.hpp"
#include "fem/substructuremodel.hpp"
#include "fixtures.hpp"
#include "linalg/blockcyclic.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"

namespace palace
{

TEST_CASE("Interface true DOFs are identified consistently in parallel",
          "[substructure][Serial][Parallel]")
{
  // The global interface true-DOF count must be independent of the MPI partition. Cube
  // split at x=0.5, H1: 49 interface DOFs at order 1 (7x7 plane), 169 at order 2.
  auto count = [](int order)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::H1_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&mesh, &fec);
    mfem::Array<int> ra(1), ea(1);
    ra[0] = 1;
    ea[0] = 2;
    mfem::Array<int> rm, em;
    return MarkInterfaceTrueDofs(pfes, ra, ea, rm, em);
  };
  CHECK(count(1) == 49);
  CHECK(count(2) == 169);
}

TEST_CASE("Element moments tell sets of elements apart", "[substructure][Serial][Parallel]")
{
  // A box [0, 2] x [-0.3, 0.3]^2, and the same rotated by 90 degrees about the z axis: the
  // same count, volume and distances from the origin, told apart by the centroid moments
  // (of which some vanish, up to rounding). The moments are those of any partition, and a
  // subset differs.
  auto box = [](bool rotate)
  {
    mfem::Mesh serial =
        mfem::Mesh::MakeCartesian3D(4, 2, 2, mfem::Element::HEXAHEDRON, 2.0, 0.6, 0.6);
    serial.Transform(
        [rotate](const mfem::Vector &x, mfem::Vector &y)
        {
          y = x;
          y(1) -= 0.3;
          y(2) -= 0.3;
          if (rotate)
          {
            y(0) = -(x(1) - 0.3);
            y(1) = x(0);
          }
        });
    return serial;
  };
  mfem::Mesh serial = box(false), rotated = box(true);
  mfem::ParMesh mesh(Mpi::World(), serial), mesh_rotated(Mpi::World(), rotated);
  std::vector<int> part(serial.GetNE());
  for (int e = 0; e < serial.GetNE(); e++)
  {
    part[e] = (serial.GetNE() - 1 - e) % Mpi::Size(Mpi::World());  // another partition
  }
  mfem::ParMesh mesh_repartitioned(Mpi::World(), serial, part.data());
  const auto m = AttributeMoments(mesh, {1})[0];
  CHECK(m[0] == 16.0);
  CHECK(std::abs(m[1] - 0.72) <= 1.0e-12);
  CHECK(MomentsMatch(m, AttributeMoments(mesh_repartitioned, {1})[0], 1.0e-12));
  const auto r = AttributeMoments(mesh_rotated, {1})[0];
  CHECK(r[0] == m[0]);
  CHECK(std::abs(r[1] - m[1]) <= 1.0e-12 * m[1]);
  CHECK(std::abs(r[2] - m[2]) <= 1.0e-12 * m[2]);
  CHECK_FALSE(MomentsMatch(m, r, 1.0e-9));
  std::vector<char> half(mesh.GetNE(), 0);
  for (int e = 0; e < mesh.GetNE(); e += 2)
  {
    half[e] = 1;
  }
  CHECK_FALSE(MomentsMatch(m, AttributeMoments(mesh, {1}, half)[0], 1.0e-9));
}

TEST_CASE("Configuration hash of the environment", "[substructure][Serial]")
{
  // The hash of the configuration on an environment (domain attributes 2 and 3, boundary
  // attributes 4, 5, 8 and 9) sees no change off it, of postprocessing, or of the order of
  // entries and lists, but any change on it; without the excitations, neither the indices
  // nor the sources.
  using json = nlohmann::json;
  const json base = json::parse(R"({
    "Problem": {"Type": "Magnetostatic", "Output": "postpro"},
    "Model": {"Mesh": "mesh.msh", "Refinement": {}},
    "Domains": {
      "Materials": [{"Attributes": [1, 2], "Permeability": 2.0},
                    {"Attributes": [3], "Permeability": [5.0, 1.0, 1.0],
                     "MaterialAxes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}],
      "Postprocessing": {"Energy": [{"Index": 1, "Attributes": [2]}]}},
    "Boundaries": {
      "PEC": {"Attributes": [4, 7]},
      "Superconductor": [{"Attributes": [8], "PenetrationDepth": 0.1, "Thickness": 0.01}],
      "SurfaceCurrent": [{"Index": 1, "Elements": [{"Attributes": [5], "Direction": "+X"},
                                                   {"Attributes": [6], "Direction": "-X"}]}],
      "FluxLoop": [{"Index": 2, "FilmAttributes": [8], "HoleAttributes": [9]}],
      "Periodic": {"BoundaryPairs": [
          {"DonorAttributes": [4], "ReceiverAttributes": [9], "Translation": [1, 0, 0]},
          {"DonorAttributes": [6], "ReceiverAttributes": [7], "Translation": [0, 1, 0]}]},
      "Postprocessing": {"SurfaceFlux": [{"Index": 1, "Attributes": [5], "Type": "Magnetic"}]}},
    "Solver": {"Order": 2, "QuadratureOrderExtra": 0,
               "Substructuring": {"SaveModel": true}}})");
  const std::vector<int> attrs = {3, 2}, bdr = {9, 8, 5, 4};
  for (const bool excitations : {false, true})
  {
    CAPTURE(excitations);
    const auto hash = EnvironmentConfigHash(base, attrs, bdr, excitations);
    auto differs = [&](const std::function<void(json &)> &change)
    {
      json config = base;
      change(config);
      return EnvironmentConfigHash(config, attrs, bdr, excitations) != hash;
    };
    auto &mat = base["Domains"]["Materials"];
    // Off the environment, of postprocessing, the order, and the run.
    CHECK_FALSE(differs(
        [&](json &c)
        {
          c["Domains"]["Materials"] = {{{"Attributes", {3}},
                                        {"Permeability", mat[1]["Permeability"]},
                                        {"MaterialAxes", mat[1]["MaterialAxes"]}},
                                       {{"Attributes", {2, 2}}, {"Permeability", 2.0}},
                                       {{"Attributes", {1}}, {"Permeability", 7.0}}};
        }));
    CHECK_FALSE(differs([](json &c) { c["Boundaries"]["PEC"]["Attributes"] = {4}; }));
    CHECK_FALSE(
        differs([](json &c) { c["Boundaries"]["PEC"]["Attributes"] = {7, 4, 1, 4}; }));
    CHECK_FALSE(differs(
        [](json &c)
        { c["Boundaries"]["SurfaceCurrent"][0]["Elements"][1]["Direction"] = "+Y"; }));
    CHECK_FALSE(differs(
        [](json &c)
        {
          c["Domains"]["Postprocessing"]["Energy"][0]["Attributes"] = {3};
          c["Boundaries"]["Postprocessing"]["SurfaceFlux"][0]["Attributes"] = {4};
        }));
    CHECK_FALSE(differs(
        [](json &c)
        {
          c["Problem"]["Output"] = "online";
          c["Model"]["Mesh"] = "region.mesh";
          c["Model"]["Refinement"]["UniformLevels"] = 1;
          c["Solver"]["Substructuring"] = {{"ModelFile", "postpro/environment.model"}};
        }));
    CHECK_FALSE(differs(
        [](json &c)
        { c["Boundaries"]["Periodic"]["BoundaryPairs"][1]["Translation"] = {0, 2, 0}; }));
    // As parsed: a material, a number and a boundary written another way.
    CHECK_FALSE(differs([](json &c) { c["Domains"]["Materials"][0]["Permeability"] = 2; }));
    CHECK_FALSE(differs(
        [](json &c) { c["Domains"]["Materials"][0]["Permeability"] = {2.0, 2.0, 2.0}; }));
    CHECK_FALSE(differs(
        [](json &c)
        {
          c["Domains"]["Materials"][0]["MaterialAxes"] = {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}};
        }));
    CHECK_FALSE(differs(
        [](json &c)
        {
          c["Boundaries"]["Ground"] = c["Boundaries"]["PEC"];
          c["Boundaries"].erase("PEC");
        }));
    // On the environment.
    CHECK(differs(
        [](json &c)
        { c["Boundaries"]["Periodic"]["BoundaryPairs"][0]["Translation"] = {2, 0, 0}; }));
    CHECK(differs([](json &c) { c["Domains"]["Materials"][1]["Permeability"][0] = 6.0; }));
    CHECK(differs([](json &c)
                  { c["Domains"]["Materials"][1]["MaterialAxes"][0] = {0, 1, 0}; }));
    CHECK(differs([](json &c) { c["Domains"]["Materials"][0]["Attributes"] = {1}; }));
    CHECK(differs([](json &c) { c["Domains"]["Materials"][1]["Attributes"] = {2, 3}; }));
    CHECK(differs([](json &c) { c["Boundaries"]["PEC"]["Attributes"] = {7}; }));
    CHECK(differs([](json &c) { c["Boundaries"]["PMC"]["Attributes"] = {9}; }));
    CHECK(differs([](json &c)
                  { c["Boundaries"]["Superconductor"][0]["PenetrationDepth"] = 0.2; }));
    CHECK(differs([](json &c) { c["Solver"]["Order"] = 3; }));
    CHECK(differs([](json &c) { c["Solver"]["QuadratureOrderExtra"] = 1; }));
    CHECK(differs([](json &c) { c["Model"]["CrackInternalBoundaryElements"] = true; }));
    CHECK(differs([](json &c) { c["Model"]["L0"] = 1.0e-5; }));
    // The polarity of a wave port, signal then ground, is ordered.
    auto polarity = [&](int signal, int ground)
    {
      json config = base;
      config["Boundaries"]["WavePort"] = {
          {{"Index", 3}, {"Attributes", {9}}, {"PolarityAttributes", {signal, ground}}}};
      return EnvironmentConfigHash(config, attrs, bdr, excitations);
    };
    CHECK(polarity(4, 5) != polarity(5, 4));
    CHECK(differs([](json &c) { c["Problem"]["Type"] = "Electrostatic"; }));
    // The excitations.
    CHECK(differs([](json &c) { c["Boundaries"]["FluxLoop"][0]["Index"] = 3; }) ==
          excitations);
    CHECK(
        differs([](json &c) { c["Boundaries"]["FluxLoop"][0]["HoleAttributes"] = {5}; }) ==
        excitations);
    CHECK(differs(
              [](json &c)
              {
                c["Boundaries"]["SurfaceCurrent"][0]["Elements"][0]["Direction"] = "+Y";
              }) == excitations);
  }
}

TEST_CASE("Environment records are those of any partition",
          "[substructure][Serial][Parallel]")
{
  // The boundary attributes next to a domain (either side of an interior boundary), the
  // moments of its elements and of the boundary elements within it (not those between the
  // domains), and so the record of an environment, are those of any partition; the record
  // survives its bytes, and other bytes are none.
  mfem::Mesh serial = MakeInteriorBoundaryBox();
  std::vector<int> part(serial.GetNE());
  for (int e = 0; e < serial.GetNE(); e++)
  {
    part[e] = (serial.GetNE() - 1 - e) % Mpi::Size(Mpi::World());
  }
  mfem::ParMesh mesh(Mpi::World(), serial),
      repartitioned(Mpi::World(), serial, part.data());
  CHECK(BoundaryAttributesNextTo(mesh, {2}) == std::vector<int>{2, 3, 4, 5});
  CHECK(BoundaryAttributesNextTo(mesh, {1}) == std::vector<int>{1, 3, 4});
  const auto m = BoundaryAttributeMoments(mesh, {2, 3, 4, 5}, {2});
  REQUIRE(m.size() == 4);
  const double count[4] = {4.0, 16.0, 0.0, 4.0}, area[4] = {1.0, 4.0, 0.0, 1.0};
  for (int a = 0; a < 4; a++)
  {
    CHECK(m[a][0] == count[a]);
    CHECK(std::abs(m[a][1] - area[a]) <= 1.0e-12);
  }
  const auto config = nlohmann::json::parse(
      R"({"Domains": {"Materials": [{"Attributes": [2], "Permittivity": 2.0}]}})");
  std::vector<char> written(mesh.GetNE(), 0), written_rep(repartitioned.GetNE(), 0);
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    written[e] = (c(0) < 1.5);
  }
  for (int e = 0; e < repartitioned.GetNE(); e++)
  {
    repartitioned.GetElementCenter(e, c);
    written_rep[e] = (c(0) < 1.5);
  }
  const auto r = EnvironmentRecord::Make(mesh, {2}, config, written),
             s = EnvironmentRecord::Make(repartitioned, {2}, config, written_rep);
  CHECK(r.bdr_attributes == s.bdr_attributes);
  CHECK(r.config_hash == s.config_hash);
  CHECK(r.excitations_hash == s.excitations_hash);
  for (const auto v : {&EnvironmentRecord::moments, &EnvironmentRecord::region_moments,
                       &EnvironmentRecord::bdr_moments})
  {
    REQUIRE((r.*v).size() == (s.*v).size());
    for (std::size_t a = 0; a < (r.*v).size(); a++)
    {
      CHECK(MomentsMatch((r.*v)[a], (s.*v)[a], 1.0e-12));
    }
  }
  CHECK(r.region_moments[0][0] == 4.0);
  const auto bytes = r.Serialize();
  const auto t = EnvironmentRecord::Deserialize(bytes);
  REQUIRE(t);
  CHECK(t->attributes == r.attributes);
  CHECK(t->bdr_attributes == r.bdr_attributes);
  CHECK(t->moments == r.moments);
  CHECK(t->region_moments == r.region_moments);
  CHECK(t->bdr_moments == r.bdr_moments);
  CHECK(t->config_hash == r.config_hash);
  CHECK(t->excitations_hash == r.excitations_hash);
  CHECK_FALSE(
      EnvironmentRecord::Deserialize(std::vector<char>(bytes.begin(), bytes.end() - 1)));
  auto longer = bytes;
  longer.push_back(0);
  CHECK_FALSE(EnvironmentRecord::Deserialize(longer));
  // The full mesh is told from the mesh of the region, and the record tells whether the
  // excitations of the environment are its own (a flux loop added next to it changes them,
  // not its configuration).
  CHECK_FALSE(r.CheckMesh(mesh, config, 1.0e-9));
  CHECK(r.CheckConfig(config, {2}, false));
  auto excited = config;
  excited["Boundaries"]["FluxLoop"] = {{{"Index", 1},
                                        {"FilmAttributes", {3}},
                                        {"HoleAttributes", {4}},
                                        {"FluxAmounts", {1.0}},
                                        {"Direction", "+Z"}}};
  CHECK_FALSE(r.CheckConfig(excited, {2}, false));
  // A port next to the environment is of its configuration, but whether and how it is
  // excited is of its excitations.
  auto port = config;
  port["Boundaries"]["LumpedPort"] = {
      {{"Index", 1}, {"R", 50.0}, {"Attributes", {5}}, {"Direction", "+Z"}}};
  const auto q = EnvironmentRecord::Make(mesh, {2}, port, written);
  CHECK(q.config_hash != r.config_hash);
  port["Boundaries"]["LumpedPort"][0]["Excitation"] = true;
  port["Boundaries"]["LumpedPort"][0]["ExcitationAmplitude"] = 2.0;
  CHECK_FALSE(q.CheckConfig(port, {2}, false));
}

TEST_CASE("Boundaries and points of the region", "[substructure][Serial][Parallel]")
{
  // Unit cube with the region at x < 1/2 and interior boundaries at x = 1/4 (attribute 7)
  // and x = 1/2 (8), on processes alternating in x: the boundaries at x = 0 (5) and x = 1/4
  // border the region only, the one at x = 1 (3) the environment only, and y = 0 (2) and
  // x = 1/2 both. A probe is in the region if it is found in one of its elements, and one
  // on x = 1/2 is found on one side.
  mfem::Mesh cube = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON);
  mfem::Mesh serial(3, cube.GetNV(), cube.GetNE(), cube.GetNBE() + 32);
  for (int v = 0; v < cube.GetNV(); v++)
  {
    serial.AddVertex(cube.GetVertex(v));
  }
  mfem::Array<int> vtx;
  mfem::Vector c;
  for (int e = 0; e < cube.GetNE(); e++)
  {
    cube.GetElementVertices(e, vtx);
    cube.GetElementCenter(e, c);
    serial.AddHex(vtx.GetData(), (c(0) < 0.5) ? 1 : 2);
  }
  for (int b = 0; b < cube.GetNBE(); b++)
  {
    cube.GetBdrElementVertices(b, vtx);
    serial.AddBdrQuad(vtx.GetData(), cube.GetBdrAttribute(b));
  }
  for (int f = 0; f < cube.GetNFaces(); f++)
  {
    cube.GetFaceVertices(f, vtx);
    double lo = 1.0, hi = 0.0;
    for (int v : vtx)
    {
      lo = std::min(lo, cube.GetVertex(v)[0]);
      hi = std::max(hi, cube.GetVertex(v)[0]);
    }
    if (hi - lo < 1.0e-12 &&
        (std::abs(lo - 0.25) < 1.0e-12 || std::abs(lo - 0.5) < 1.0e-12))
    {
      serial.AddBdrQuad(vtx.GetData(), (lo < 0.4) ? 7 : 8);
    }
  }
  serial.FinalizeHexMesh(1, 0, true);
  std::vector<int> part(serial.GetNE());
  for (int e = 0; e < serial.GetNE(); e++)
  {
    serial.GetElementCenter(e, c);
    part[e] = (static_cast<int>(4.0 * c(0)) + 4 * static_cast<int>(4.0 * c(1))) %
              Mpi::Size(Mpi::World());
  }
  mfem::ParMesh mesh(Mpi::World(), serial, part.data());
  const std::vector<int> region = {1};
  CHECK(BoundaryAttributesInRegion(mesh, region) == std::vector<int>{5, 7});
  CHECK(BoundaryAttributesInRegion(mesh, {2}) == std::vector<int>{3});
  CHECK(BoundaryAttributesInRegion(mesh, {1, 2}) ==
        std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8});
#if defined(MFEM_USE_GSLIB)
  mesh.EnsureNodes();  // as a Palace mesh, for the point search
  std::map<int, config::ProbeData> probes;
  probes[1].center = {0.2, 0.4, 0.6};
  probes[2].center = {0.8, 0.4, 0.6};
  probes[3].center = {2.0, 0.4, 0.6};
  CHECK(fem::ProbesInDomains(probes, mesh, region) == std::vector<int>{1});
  CHECK(fem::ProbesInDomains(probes, mesh, {2}) == std::vector<int>{2});
  // A probe on the face between the domains is in both.
  std::map<int, config::ProbeData> on_gamma;
  on_gamma[1].center = {0.5, 0.4, 0.6};
  CHECK(fem::ProbesInDomains(on_gamma, mesh, region) == std::vector<int>{1});
  CHECK(fem::ProbesInDomains(on_gamma, mesh, {2}) == std::vector<int>{1});
#endif
}

TEST_CASE("Nédélec DOFs at the points of saved ones", "[substructure][Serial]")
{
  // Point-tangent DOFs u(p).t, with signatures (t_b, p_a t_b): a DOF at a saved point
  // matches whatever its tangent (another one, as a face DOF in another basis, or the
  // opposite one), a DOF off the saved points does not.
  auto add = [](std::vector<double> &sig, std::array<double, 3> p, std::array<double, 3> t)
  {
    sig.insert(sig.end(), t.begin(), t.end());
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        sig.push_back(p[a] * t[b]);
      }
    }
  };
  std::vector<double> saved, candidates;
  add(saved, {0.1, 0.2, 0.3}, {1.0, 0.0, 0.0});
  add(saved, {0.5, -0.2, 0.0}, {0.0, 0.6, 0.8});
  add(saved, {0.5, 0.4, 0.0}, {0.0, 0.0, 1.0});
  add(candidates, {0.1, 0.2, 0.3}, {0.0, 1.0, 0.0});
  add(candidates, {0.5, -0.2, 0.0}, {0.0, -0.6, -0.8});
  add(candidates, {0.5, -0.2, 1.0e-6}, {0.0, 0.6, 0.8});
  add(candidates, {0.1 + 1.0e-12, 0.2, 0.3}, {1.0, 0.0, 0.0});
  add(candidates, {0.5, 0.4, 0.0}, {0.6, 0.8, 0.0});
  const auto match = MatchSignaturePoints(candidates, 12, saved, 12);
  REQUIRE(match.size() == 5);
  CHECK(match[0]);
  CHECK(match[1]);
  CHECK_FALSE(match[2]);
  CHECK(match[3]);
  CHECK(match[4]);
}

TEST_CASE("Versions of split Nédélec DOFs are matched by their sides",
          "[substructure][Serial]")
{
  // The versions of a DOF split by a two-sided film share its 12 moments and have side
  // moments of opposite signs; a DOF away from the film has none. An orientation flip
  // negates all 13 columns. Each current version must map to the saved version of its side
  // (of the same 12 moments, which alone would pair the first one with the wrong side).
  auto add = [](std::vector<double> &sig, std::array<double, 3> p, std::array<double, 3> t,
                double side)
  {
    sig.insert(sig.end(), t.begin(), t.end());
    for (int a = 0; a < 3; a++)
    {
      for (int b = 0; b < 3; b++)
      {
        sig.push_back(p[a] * t[b]);
      }
    }
    sig.push_back(side);
  };
  const std::array<double, 3> p = {0.0, 0.3, 0.0}, t = {0.0, 0.1, 0.0}, q = {0.5, 0.2, 0.1},
                              u = {0.1, 0.0, 0.05};
  const double side = 2.0e-3;
  std::vector<double> saved, current;
  add(saved, p, t, -side);
  add(saved, p, t, side);
  add(saved, q, u, 0.0);
  add(current, p, t, -side);
  add(current, q, {-u[0], -u[1], -u[2]}, 0.0);
  add(current, p, {-t[0], -t[1], -t[2]}, -side);
  const SignatureMap map = MatchSignatureBasis(current, saved, 13);
  REQUIRE(map.rows.size() == 3);
  using Entry = std::vector<std::pair<int, double>>;
  CHECK(map.rows[0] == Entry{{0, 1.0}});
  CHECK(map.rows[1] == Entry{{2, -1.0}});
  CHECK(map.rows[2] == Entry{{1, -1.0}});

  // A version whose side moment matches neither saved one is rejected, even though the
  // combinations of the saved versions reproduce it.
  current[12] *= 2.0;
  CHECK_THROWS(MatchSignatureBasis(current, saved, 13));
}

TEST_CASE("H1 DOFs at the points of saved ones", "[substructure][Serial]")
{
  // H1 signatures are the DOF coordinates: a DOF matches a saved one at its point, within
  // the matching tolerance, and not off it.
  const std::vector<double> saved = {0.1, 0.2, 0.3, 0.5, -0.2, 0.0},
                            candidates = {0.5, -0.2, 0.0,  0.1,   0.2 + 1.0e-12,
                                          0.3, 0.5,  -0.2, 1.0e-6};
  const auto match = MatchSignaturePoints(candidates, 3, saved, 3);
  REQUIRE(match.size() == 3);
  CHECK(match[0]);
  CHECK(match[1]);
  CHECK_FALSE(match[2]);
}

TEST_CASE("Saved DOFs are matched in large sets", "[substructure][Serial]")
{
  // Random signatures, permuted, perturbed below the matching tolerance and, for Nédélec
  // DOFs, flipped: the map recovers the permutation and the flips. H1 DOFs (points) are
  // never flipped: points mirrored through the origin do not match.
  std::mt19937 gen(7);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  const int n = 20000;
  for (const int width : {3, 12})
  {
    std::vector<double> saved(static_cast<std::size_t>(n) * width), current(saved.size());
    for (auto &x : saved)
    {
      x = u(gen);
    }
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), gen);
    std::vector<double> sign(n);
    for (int g = 0; g < n; g++)
    {
      sign[g] = (width == 12 && g % 3 == 0) ? -1.0 : 1.0;
      for (int q = 0; q < width; q++)
      {
        current[static_cast<std::size_t>(g) * width + q] =
            sign[g] * saved[static_cast<std::size_t>(perm[g]) * width + q] +
            1.0e-10 * u(gen);
      }
    }
    const SignatureMap map = MatchSignatureBasis(current, saved, width);
    int wrong = 0;
    for (int g = 0; g < n; g++)
    {
      wrong += (map.rows[g] != std::vector<std::pair<int, double>>{{perm[g], sign[g]}});
    }
    CAPTURE(width);
    CHECK(wrong == 0);
    if (width == 3)
    {
      for (auto &x : current)
      {
        x = -x;
      }
      CHECK_THROWS(MatchSignatureBasis(current, saved, width));
    }
  }
}

TEST_CASE("Saved Nédélec DOFs are mapped across meshes and partitions",
          "[substructure][Serial][Parallel]")
{
  // The same tetrahedra in a different element and vertex order (so different face
  // orientations, as a different partition gives): the true-DOF bases of second-order
  // Nédélec elements differ beyond signs on the faces. The signature map must transform a
  // linear form and the mass matrix of one space into those of the other (dual quantities,
  // x_cur = M^-T x_saved).
  const int order = GENERATE(1, 2);
  CAPTURE(order);
  mfem::Mesh base = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::Mesh perm(3, base.GetNV(), base.GetNE());
  const int nv = base.GetNV();
  for (int v = 0; v < nv; v++)
  {
    perm.AddVertex(base.GetVertex((v * 7) % nv));  // 7 and nv = 27 are coprime
  }
  std::vector<int> inv(nv);
  for (int v = 0; v < nv; v++)
  {
    inv[(v * 7) % nv] = v;
  }
  for (int e = base.GetNE() - 1; e >= 0; e--)
  {
    const int *ev = base.GetElement(e)->GetVertices();
    // An even permutation of the vertices keeps the orientation.
    perm.AddTet(inv[ev[1]], inv[ev[2]], inv[ev[0]], inv[ev[3]]);
  }
  perm.FinalizeTetMesh(1, 0, true);
  auto measure = [order](mfem::Mesh &serial, std::vector<double> &sig,
                         std::vector<double> &b, std::vector<double> &M)
  {
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::ND_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace fes(&mesh, &fec);
    const int nt = fes.GetTrueVSize();
    const int n = static_cast<int>(fes.GlobalTrueVSize());
    std::vector<int> index(nt);
    for (int i = 0; i < nt; i++)
    {
      index[i] = static_cast<int>(fes.GetMyTDofOffset()) + i;
    }
    sig = TrueDofSignatures(fes, index, n);
    // A linear form and the mass matrix, replicated (n is small).
    mfem::VectorFunctionCoefficient f(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v.SetSize(3);
                                        v(0) = std::sin(3.0 * x(1)) + x(2);
                                        v(1) = x(0) * x(2) * x(2);
                                        v(2) = std::cos(2.0 * x(0) + x(1));
                                      });
    mfem::ParLinearForm lf(&fes);
    lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f));
    lf.Assemble();
    std::unique_ptr<mfem::HypreParVector> bt(lf.ParallelAssemble());
    mfem::ParBilinearForm a(&fes);
    a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    a.Assemble();
    a.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
    b.assign(n, 0.0);
    M.assign(static_cast<std::size_t>(n) * n, 0.0);
    mfem::Vector e(nt), Ae(nt);
    for (int j = 0; j < n; j++)
    {
      e = 0.0;
      const int jl = j - static_cast<int>(fes.GetMyTDofOffset());
      if (jl >= 0 && jl < nt)
      {
        e(jl) = 1.0;
      }
      A->Mult(e, Ae);
      const double *ae = Ae.HostRead();
      for (int i = 0; i < nt; i++)
      {
        M[static_cast<std::size_t>(j) * n + index[i]] = ae[i];
      }
    }
    const double *bh = bt->HostRead();
    for (int i = 0; i < nt; i++)
    {
      b[index[i]] = bh[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, b.data(), n, MPI_DOUBLE, MPI_SUM, Mpi::World());
    MPI_Allreduce(MPI_IN_PLACE, M.data(), n * n, MPI_DOUBLE, MPI_SUM, Mpi::World());
  };
  std::vector<double> sig_s, b_s, M_s, sig_c, b_c, M_c;
  measure(base, sig_s, b_s, M_s);
  measure(perm, sig_c, b_c, M_c);
  const SignatureMap map = MatchSignatureBasis(sig_c, sig_s, 12);
  int blocks = 0;
  for (const auto &r : map.rows)
  {
    blocks += (r.size() > 1);
  }
  if (order > 1)
  {
    CHECK(blocks > 0);  // face DOFs whose basis is not a signed permutation of the other's
  }
  const auto b = map.DualRows(b_s, 1);
  // M in thirds of its rows (as ranks hold them), each from the saved rows it combines.
  const int n = static_cast<int>(b_s.size()), third = (n + 2) / 3;
  std::vector<double> M;
  for (int g0 = 0; g0 < n; g0 += third)
  {
    const int m = std::min(third, n - g0);
    std::vector<int> read;
    const auto Mg =
        map.DualMatrix<double>(g0, m,
                               [&](int a, double *x)
                               {
                                 read.push_back(a);
                                 std::copy_n(&M_s[static_cast<std::size_t>(a) * n], n, x);
                               });
    CHECK(read == map.SavedRows(g0, m));
    CHECK(read.size() < static_cast<std::size_t>(n));
    M.insert(M.end(), Mg.begin(), Mg.end());
  }
  REQUIRE(M.size() == M_c.size());
  double db = 0.0, mb = 0.0, dM = 0.0, mM = 0.0;
  for (std::size_t i = 0; i < b.size(); i++)
  {
    db = std::max(db, std::abs(b[i] - b_c[i]));
    mb = std::max(mb, std::abs(b_c[i]));
  }
  for (std::size_t q = 0; q < M.size(); q++)
  {
    dM = std::max(dM, std::abs(M[q] - M_c[q]));
    mM = std::max(mM, std::abs(M_c[q]));
  }
  CAPTURE(blocks, db, mb, dM, mM);
  CHECK(db <= 1.0e-12 * mb);
  CHECK(dM <= 1.0e-12 * mM);
}

TEST_CASE("Eliminated rows without a stored diagonal become identity rows",
          "[substructure][Serial][Parallel]")
{
  // Per rank: rows 0 and 3 store no entry (as a device-assembled matrix has for a DOF
  // without element contributions), rows 1 and 2 are coupled, rows 4 and 5 diagonal.
  // Eliminating rows 0, 1 and 3 leaves the identity on them and decouples row 2.
  MPI_Comm comm = Mpi::World();
  const int n = 6, rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  const HYPRE_BigInt off = static_cast<HYPRE_BigInt>(rank) * n,
                     glob = static_cast<HYPRE_BigInt>(size) * n;
  std::vector<HYPRE_BigInt> part;
  if (HYPRE_AssumedPartitionCheck())
  {
    part = {off, off + n};
  }
  else
  {
    for (int r = 0; r <= size; r++)
    {
      part.push_back(static_cast<HYPRE_BigInt>(r) * n);
    }
  }
  std::vector<int> I = {0, 0, 2, 4, 4, 5, 6};
  std::vector<HYPRE_BigInt> J = {off + 1, off + 2, off + 1, off + 2, off + 4, off + 5};
  std::vector<double> V = {2.0, -1.0, -1.0, 2.0, 2.0, 2.0};
  mfem::HypreParMatrix A(comm, n, glob, glob, I.data(), J.data(), V.data(), part.data(),
                         part.data());
  mfem::Array<int> dofs({0, 1, 3});
  std::unique_ptr<mfem::HypreParMatrix> B = EliminateToIdentity(A, dofs);
  mfem::Vector x(n), y(n);
  x = 1.0;
  B->Mult(x, y);
  const double *yy = y.HostRead();
  const std::array<double, 6> expected = {1.0, 1.0, 2.0, 1.0, 2.0, 2.0};
  for (int i = 0; i < n; i++)
  {
    CAPTURE(i);
    CHECK(yy[i] == expected[i]);
  }
}

namespace
{

// The redistribution of the lower triangle of a symmetric n x n matrix S, each process
// holding its share, through a signed permutation of its DOFs (S_cur(g, h) = s_g s_h
// S(p_g, p_h), or S without it), onto the rows of a block of DOFs on each process and onto
// a block cyclic grid of all of them, against S_cur on every process.
template <typename T>
void CheckTriangleRedistribution(int n, bool mapped)
{
  MPI_Comm comm = Mpi::World();
  const int nranks = Mpi::Size(comm), rank = Mpi::Rank(comm);
  auto value = [](int a, int b)
  {
    const int i = std::max(a, b), j = std::min(a, b);
    if constexpr (std::is_same_v<T, double>)
    {
      return std::sin(1.0 + i + 3.1 * j);
    }
    else
    {
      return T(std::sin(1.0 + i + 3.1 * j), std::cos(0.5 * i - 2.0 * j));
    }
  };
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(7);
  std::shuffle(perm.begin(), perm.end(), rng);
  std::vector<int> saved(n);
  std::vector<std::vector<std::pair<int, double>>> to;
  if (mapped)
  {
    to.resize(n);
    for (int a = 0; a < n; a++)
    {
      to[a] = {{perm[a], (a % 3 == 0) ? -1.0 : 1.0}};
      saved[perm[a]] = a;
    }
  }
  auto expected = [&](int g, int h)
  {
    if (!mapped)
    {
      return value(g, h);
    }
    const int a = saved[g], b = saved[h];
    return ((a % 3 == 0) ? -1.0 : 1.0) * ((b % 3 == 0) ? -1.0 : 1.0) * value(a, b);
  };
  auto read = [&](long long e, long long m, T *x)
  {
    long long b = 0;
    while (substructuring::TriangleColumnStart(n, b + 1) <= e)
    {
      b++;
    }
    long long a = b + (e - substructuring::TriangleColumnStart(n, b));
    for (long long k = 0; k < m; k++)
    {
      x[k] = value(static_cast<int>(a), static_cast<int>(b));
      if (++a == n)
      {
        a = ++b;
      }
    }
  };

  // The rows [n r / p, n (r + 1) / p) on process r.
  auto first = [&](int r)
  { return static_cast<int>(static_cast<long long>(n) * r / nranks); };
  std::vector<T> rows(static_cast<std::size_t>(first(rank + 1) - first(rank)) * n, T(0.0));
  substructuring::RedistributeTriangle<T>(
      comm, n, read, to,
      [&](int g, int)
      {
        int q = 0;
        while (first(q + 1) <= g)
        {
          q++;
        }
        return q;
      },
      [&](int g, int h) { return static_cast<std::size_t>(g - first(rank)) * n + h; },
      rows);
  double err = 0.0;
  for (int g = first(rank); g < first(rank + 1); g++)
  {
    for (int h = 0; h < n; h++)
    {
      err = std::max(err, std::abs(rows[static_cast<std::size_t>(g - first(rank)) * n + h] -
                                   expected(g, h)));
    }
  }

  // The grid of all processes, in blocks of 4.
  const BlockCyclicGrid grid(comm, n, nranks, 4);
  std::vector<T> local(grid.LocalSize(), T(0.0));
  substructuring::RedistributeTriangle<T>(
      comm, n, read, to, [&](int g, int h) { return grid.Owner(g, h); },
      [&](int g, int h) { return grid.LocalIndex(g, h); }, local);
  for (int jl = 0; jl < grid.LocalCols(); jl++)
  {
    for (int il = 0; il < grid.LocalRows(); il++)
    {
      const int g = grid.GlobalRow(il), h = grid.GlobalCol(jl);
      err = std::max(err, std::abs(local[grid.LocalIndex(g, h)] - expected(g, h)));
    }
  }
  Mpi::GlobalMax(1, &err, comm);
  CHECK(err == 0.0);
}

}  // namespace

TEST_CASE("Saved lower triangles are redistributed onto rows and process grids",
          "[substructure][Serial][Parallel]")
{
  const int n = GENERATE(1, 2, 37);
  const bool mapped = GENERATE(false, true);
  CAPTURE(n, mapped);
  CheckTriangleRedistribution<double>(n, mapped);
  CheckTriangleRedistribution<std::complex<double>>(n, mapped);
}

TEST_CASE("The parts of a model record pass through rank 0 in chunks",
          "[substructure][Serial][Parallel]")
{
  // Rank q holds 3 q + 2 doubles, but rank 1 none; a chunk of 1 or 4 doubles holds part of
  // those of a rank, the largest one all of them.
  MPI_Comm comm = Mpi::World();
  const int rank = Mpi::Rank(comm), nranks = Mpi::Size(comm);
  const std::size_t chunk = GENERATE(1, 4, 1000);
  auto count = [](int q) { return (q == 1) ? 0 : 3 * q + 2; };
  auto value = [](int q, int i) { return 100.0 * q + i; };
  std::vector<long long> size(nranks);
  for (int q = 0; q < nranks; q++)
  {
    size[q] = count(q);
  }
  std::vector<double> x(count(rank)), expected;
  for (int i = 0; i < count(rank); i++)
  {
    x[i] = value(rank, i);
  }
  std::vector<double> got;
  std::size_t most = 0;
  substructuring::ToRoot<double>(comm, x, size, chunk,
                                 [&](const double *y, std::size_t n)
                                 {
                                   got.insert(got.end(), y, y + n);
                                   most = std::max(most, n);
                                 });
  for (int q = 0; Mpi::Root(comm) && q < nranks; q++)
  {
    for (int i = 0; i < count(q); i++)
    {
      expected.push_back(value(q, i));
    }
  }
  CHECK(got == expected);
  CHECK(most <= chunk);
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Records of a substructuring model file",
                 "[substructure][Serial][Parallel]")
{
  // Records appended to a model file are read back by tag, the last of a tag superseding
  // the earlier ones; a record cut short ends the complete ones, and a file of another
  // container version or none is told apart.
  using substructuring::ModelFile;
  MPI_Comm comm = Mpi::World();
  const std::string path = (temp_dir / "records.model").string();
  const auto A = ModelFile::Tag("AAAA"), B = ModelFile::Tag("BBBB");
  auto append = [&](std::uint32_t tag, std::vector<double> x)
  {
    ModelFile::Append(
        path, comm, tag, static_cast<long long>(sizeof(double) * x.size()),
        [&](std::ostream &f)
        { f.write(reinterpret_cast<const char *>(x.data()), sizeof(double) * x.size()); });
  };
  auto values = [&](const ModelFile &file, std::uint32_t tag)
  {
    const auto *r = file.Find(tag);
    REQUIRE(r);
    const std::vector<char> p = ModelFile::ReadPayload(path, *r, comm);
    std::vector<double> x(p.size() / sizeof(double));
    std::memcpy(x.data(), p.data(), p.size());
    return x;
  };
  ModelFile::Create(path, comm);
  append(A, {1.0, 2.0});
  append(B, {3.0});
  append(A, {4.0, 5.0, 6.0});
  ModelFile file = ModelFile::Read(path, comm);
  CHECK(file.version == static_cast<int>(ModelFile::kVersion));
  CHECK_FALSE(file.incomplete);
  CHECK(file.records.size() == 3);
  CHECK(values(file, A) == std::vector<double>{4.0, 5.0, 6.0});
  CHECK(values(file, B) == std::vector<double>{3.0});
  CHECK(file.Find(ModelFile::Tag("CCCC")) == nullptr);

  // A cut inside the last record leaves the first two.
  if (Mpi::Root(comm))
  {
    std::filesystem::resize_file(path, std::filesystem::file_size(path) - sizeof(double));
  }
  Mpi::Barrier(comm);
  file = ModelFile::Read(path, comm);
  CHECK(file.incomplete);
  CHECK(file.records.size() == 2);
  CHECK(values(file, A) == std::vector<double>{1.0, 2.0});

  // A record of another size than announced stops every rank.
#if defined(MFEM_USE_EXCEPTIONS)
  CHECK_THROWS_WITH(ModelFile::Append(path, comm, A, 8, [](std::ostream &) {}),
                    Catch::Matchers::ContainsSubstring("instead of 8"));
#endif

  // Another container version, a file that is not a model, and no file.
  if (Mpi::Root(comm))
  {
    std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
    const std::uint32_t version = ModelFile::kVersion + 1;
    f.seekp(sizeof(std::uint32_t));
    f.write(reinterpret_cast<const char *>(&version), sizeof(version));
  }
  Mpi::Barrier(comm);
  CHECK(ModelFile::Read(path, comm).version == static_cast<int>(ModelFile::kVersion) + 1);
  if (Mpi::Root(comm))
  {
    std::ofstream(path, std::ios::binary) << "not a model";
  }
  Mpi::Barrier(comm);
  CHECK(ModelFile::Read(path, comm).version == 0);
  if (Mpi::Root(comm))
  {
    std::filesystem::remove(path);
  }
  Mpi::Barrier(comm);
  CHECK(ModelFile::Read(path, comm).version == -1);
}

}  // namespace palace
