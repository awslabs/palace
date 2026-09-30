// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The corner family's REFINED trace basis (corner-basis refinement, USER decision 161 (1),
// 2026-09-30; design doc SURFACE-RESPONSE-IDENTIFICATION.md, Conventions
// CornerTraceBasisRule): the AllRingsFollowMetal layout — every ring of the box, the extra
// ring at MetalThickness + OveretchDepth and the two cap rings included, carries the same
// angle-dependent knot fractions (2 crossings + 5 metal-interior + 9 graded free knots, PEC
// on the two metal rings only, box corners = slaves, cap centres = slaves at the mean of
// the two crossing knots), so the basis has NO events: one interpolation segment per
// convexity, no connectivity angle, hats continuous in the angle. Pinned to the Python
// generator (generate_corner_response.REFINED_RULE, test_generate_corner_response.py). The
// load-time check of every coupon against the rule at its angle, the empty event list and
// the one-segment stencil each fail closed / are exercised here, and the runtime constructs
// the basis of an interpolated corner from the rule (SurfaceMortar lift with the centre
// slaves).

#include "fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "fem/mesh.hpp"
#include "models/cornertracebasis.hpp"
#include "models/laplaceoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace fs = std::filesystem;

using json = nlohmann::json;
using namespace Catch::Matchers;

namespace
{

constexpr double kR = 1.9, kT = 0.1, kOE = 0.05;
constexpr double kDeg = M_PI / 180.0;

ConstructedCornerTraceBasis BuildRefined(double angle, bool convex, double radius = kR,
                                         double t = kT, double oe = kOE)
{
  const auto rule = RefinedCornerTraceBasisRule();
  const auto seed = MakeCornerBoxSeed(radius, t, oe, convex, rule);
  return BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                               angle * kDeg, convex, rule);
}

double TriangleArea(const ConstructedCornerTraceBasis &basis, const std::array<int, 3> &t)
{
  const auto &a = basis.vertices[t[0]].point, &b = basis.vertices[t[1]].point,
             &c = basis.vertices[t[2]].point;
  const std::array<double, 3> ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  const std::array<double, 3> ac = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  const std::array<double, 3> cross = {ab[1] * ac[2] - ab[2] * ac[1],
                                       ab[2] * ac[0] - ab[0] * ac[2],
                                       ab[0] * ac[1] - ab[1] * ac[0]};
  return 0.5 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
}

// The piecewise-linear interpolant of vertex values on the OUTER box faces, sampled at
// (perimeter fraction s, height z) (the bands between the outer rings; the caps are left
// out): barycentric on every outer-face triangle in the (s, z) chart, unwrapped across
// s = 0. NaN where no triangle contains the sample.
std::vector<double> SampleOuterFaces(const ConstructedCornerTraceBasis &basis,
                                     const std::vector<double> &values,
                                     const std::vector<double> &sample_s,
                                     const std::vector<double> &sample_z, double radius)
{
  std::vector<double> s(basis.vertices.size(), std::nan(""));
  std::vector<bool> on_box(basis.vertices.size(), false);
  for (std::size_t v = 0; v < basis.vertices.size(); v++)
  {
    const auto &p = basis.vertices[v].point;
    on_box[v] = std::max(std::abs(p[0]), std::abs(p[1])) >= radius * (1.0 - 1.0e-9);
    if (on_box[v])
    {
      s[v] = SquarePerimeterFraction(radius, p);
    }
  }
  std::vector<double> out(sample_s.size(), std::nan(""));
  for (const auto &t : basis.triangles)
  {
    if (!(on_box[t[0]] && on_box[t[1]] && on_box[t[2]]))
    {
      continue;
    }
    std::array<double, 3> ts = {s[t[0]], s[t[1]], s[t[2]]};
    const std::array<double, 3> tz = {basis.vertices[t[0]].point[2],
                                      basis.vertices[t[1]].point[2],
                                      basis.vertices[t[2]].point[2]};
    const double span =
        *std::max_element(ts.begin(), ts.end()) - *std::min_element(ts.begin(), ts.end());
    if (span > 0.5)
    {
      for (double &f : ts)
      {
        if (f < 0.5)
        {
          f += 1.0;
        }
      }
    }
    const double det =
        (tz[1] - tz[2]) * (ts[0] - ts[2]) + (ts[2] - ts[1]) * (tz[0] - tz[2]);
    if (std::abs(det) < 1.0e-18)
    {
      continue;
    }
    for (std::size_t i = 0; i < sample_s.size(); i++)
    {
      if (!std::isnan(out[i]))
      {
        continue;
      }
      for (const double shift : {0.0, 1.0})
      {
        const double xs = sample_s[i] + shift, z = sample_z[i];
        const double l1 =
            ((tz[1] - tz[2]) * (xs - ts[2]) + (ts[2] - ts[1]) * (z - tz[2])) / det;
        const double l2 =
            ((tz[2] - tz[0]) * (xs - ts[2]) + (ts[0] - ts[2]) * (z - tz[2])) / det;
        const double l3 = 1.0 - l1 - l2;
        if (l1 >= -1.0e-9 && l2 >= -1.0e-9 && l3 >= -1.0e-9)
        {
          out[i] = l1 * values[t[0]] + l2 * values[t[1]] + l3 * values[t[2]];
          break;
        }
      }
    }
  }
  return out;
}

// A smooth trace on the basis: the held-out polynomial at every knot (zero on the PEC
// knots), slaves by their parents.
std::vector<double> SmoothTraceValues(const ConstructedCornerTraceBasis &basis,
                                      double radius)
{
  std::vector<double> values(basis.vertices.size(), 0.0);
  for (std::size_t k = 0; k < basis.knots.size(); k++)
  {
    const auto &p = basis.knots[k];
    const double x = p[0] / radius, y = p[1] / radius, z = p[2] / radius;
    values[k] = basis.zero[k] ? 0.0 : 0.35 + 0.2 * x - 0.15 * y + 0.1 * z + 0.08 * x * y;
  }
  for (std::size_t v = basis.knots.size(); v < basis.vertices.size(); v++)
  {
    const auto &vertex = basis.vertices[v];
    values[v] = vertex.weight_a * values[vertex.parent_a] +
                (1.0 - vertex.weight_a) * values[vertex.parent_b];
  }
  return values;
}

double MaxJumpAcross(const CornerTraceBasisRule &rule, double angle, bool convex,
                     double epsilon)
{
  const auto seed = MakeCornerBoxSeed(kR, kT, kOE, convex, rule);
  auto Sample = [&](double a)
  {
    const auto basis = BuildCornerTraceBasis(
        seed.points, seed.contour_groups, seed.zero_trace_indices, a * kDeg, convex, rule);
    std::vector<double> ss, zs;
    for (int i = 0; i < 400; i++)
    {
      for (const double z :
           {-0.6, -0.3, -0.04, -0.01, 0.02, 0.08, 0.12, 0.2, 0.5, 1.0, 1.5})
      {
        ss.push_back(i / 400.0);
        zs.push_back(z);
      }
    }
    return SampleOuterFaces(basis, SmoothTraceValues(basis, kR), ss, zs, kR);
  };
  const auto below = Sample(angle - epsilon), above = Sample(angle + epsilon);
  double jump = 0.0;
  for (std::size_t i = 0; i < below.size(); i++)
  {
    if (!std::isnan(below[i]) && !std::isnan(above[i]))
    {
      jump = std::max(jump, std::abs(below[i] - above[i]));
    }
  }
  return jump;
}

struct CouponFiles
{
  fs::path points, vertices, triangles;
  std::vector<int> zero_trace_indices;  // 1-based
  std::vector<int> contour_groups;
};

// A unit-test coupon of the refined rule at `angle` (basis points, trace mesh; with
// `perturb_knot` the free knot of that 0-based index is moved 1e-6 R along its ring: a
// coupon that is NOT the rule's).
CouponFiles WriteRefinedCoupon(const fs::path &directory, const std::string &tag,
                               double angle, bool convex, double radius, double t,
                               double oe, std::optional<int> perturb_knot = std::nullopt)
{
  const auto rule = RefinedCornerTraceBasisRule();
  const auto seed = MakeCornerBoxSeed(radius, t, oe, convex, rule);
  const auto basis =
      BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                            angle * kDeg, convex, rule);
  auto points_basis = basis;
  if (perturb_knot)
  {
    auto &p = points_basis.knots[*perturb_knot];
    MFEM_VERIFY(!basis.zero[*perturb_knot], "Perturb a free knot!");
    const double half_width = std::max(std::abs(p[0]), std::abs(p[1]));
    p = SquarePerimeterPoint(half_width, p[2],
                             SquarePerimeterFraction(half_width, p) + 1.0e-6 / 8.0);
  }
  CouponFiles files;
  files.points = directory / ("refined-" + tag + "-points.csv");
  files.vertices = directory / ("refined-" + tag + "-vertices.csv");
  files.triangles = directory / ("refined-" + tag + "-triangles.csv");
  files.contour_groups = seed.contour_groups;
  {
    std::ofstream output(files.points);
    output << std::setprecision(17) << "x,y,z\n";
    for (std::size_t k = 0; k < points_basis.knots.size(); k++)
    {
      output << points_basis.knots[k][0] << "," << points_basis.knots[k][1] << ","
             << points_basis.knots[k][2] << "\n";
      if (basis.zero[k])
      {
        files.zero_trace_indices.push_back(static_cast<int>(k) + 1);
      }
    }
  }
  {
    std::ofstream output(files.vertices);
    output << std::setprecision(17)
           << "vertex,x,y,z,basis,conductor,parent_a,parent_b,weight_a\n";
    for (std::size_t v = 0; v < basis.vertices.size(); v++)
    {
      const auto &vertex = basis.vertices[v];
      output << v + 1 << "," << vertex.point[0] << "," << vertex.point[1] << ","
             << vertex.point[2] << ",";
      if (vertex.basis >= 0)
      {
        output << vertex.basis + 1 << "," << (basis.zero[vertex.basis] ? 1 : 0)
               << ",0,0,0\n";
      }
      else
      {
        output << "0,0," << vertex.parent_a + 1 << "," << vertex.parent_b + 1 << ","
               << vertex.weight_a << "\n";
      }
    }
  }
  {
    std::ofstream output(files.triangles);
    output << "triangle,vertex_i,vertex_j,vertex_k\n";
    for (std::size_t t = 0; t < basis.triangles.size(); t++)
    {
      output << t + 1 << "," << basis.triangles[t][0] + 1 << ","
             << basis.triangles[t][1] + 1 << "," << basis.triangles[t][2] + 1 << "\n";
    }
  }
  return files;
}

json RefinedTraceBasisRecord()
{
  const auto rule = RefinedCornerTraceBasisRule();
  return {{"RingLayout", "AllRingsFollowMetal"},
          {"RingSize", rule.ring_size},
          {"MetalInteriorKnots", rule.metal_interior_knots},
          {"FreeKnots", rule.free_knots},
          {"FreeKnotGrading", rule.free_knot_grading},
          {"ExtraLevelsAboveOverOveretch", rule.extra_levels_above_over_overetch},
          {"Fractions", rule.fractions}};
}

void WriteSyntheticMatrices(const fs::path &domain_path, const fs::path &surface_path,
                            int size, double diagonal, double coupling_scale,
                            double radius_m)
{
  std::ofstream domain(domain_path);
  domain << "basis_i,basis_j,Q_ij (J)\n";
  std::ofstream surface(surface_path);
  surface << "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n";
  for (int i = 0; i < size; i++)
  {
    for (int j = i; j < size; j++)
    {
      const double value =
          (i == j ? diagonal : coupling_scale / (1.0 + std::abs(i - j))) * 1.0e-12;
      domain << i + 1 << "," << j + 1 << "," << value << "\n";
      surface << "1,1," << radius_m << "," << i + 1 << "," << j + 1 << "," << value << ","
              << value << "\n";
    }
  }
}

// A hexahedral box mesh whose metal island (cracked boundary attribute 9 on the plane
// y = 0.5) is an exact star-shaped polygon (test-surfaceresponseoperator.cpp's
// MakePolygonIslandMesh: the radial map of every concentric square of the (x, z) grid onto
// the polygon, blended to the identity between the levels 2 and 3).
std::unique_ptr<mfem::ParMesh>
MakePolygonIslandMesh(const std::vector<std::array<double, 2>> &polygon, double extent,
                      double h)
{
  const int n = static_cast<int>(std::lround(extent / h));
  mfem::Mesh serial =
      mfem::Mesh::MakeCartesian3D(n, 4, n, mfem::Element::HEXAHEDRON, extent, 1.0, extent);
  const double c = 0.5 * extent;
  auto Level = [&](const double *point)
  { return std::max(std::abs(point[0] - c), std::abs(point[2] - c)); };
  for (int face = 0; face < serial.GetNumFaces(); face++)
  {
    int element1, element2;
    serial.GetFaceElements(face, &element1, &element2);
    if (element1 < 0 || element2 < 0)
    {
      continue;
    }
    mfem::Array<int> vertices;
    serial.GetFaceVertices(face, vertices);
    bool on_plane = true;
    double level_max = 0.0;
    for (const int vertex : vertices)
    {
      const double *point = serial.GetVertex(vertex);
      on_plane = on_plane && std::abs(point[1] - 0.5) < 1.0e-12;
      level_max = std::max(level_max, Level(point));
    }
    if (on_plane && level_max <= 1.0 + 1.0e-9)
    {
      serial.AddBdrElement(serial.GetFace(face)->Duplicate(&serial));
      serial.SetBdrAttribute(serial.GetNBE() - 1, 9);
    }
  }
  auto PolygonRadius = [&](double ux, double uz)
  {
    double r = mfem::infinity();
    for (std::size_t k = 0; k < polygon.size(); k++)
    {
      const auto &a = polygon[k], &b = polygon[(k + 1) % polygon.size()];
      const double nx = b[1] - a[1], nz = a[0] - b[0];
      const double denominator = nx * ux + nz * uz;
      if (denominator > 1.0e-14)
      {
        r = std::min(r, (nx * a[0] + nz * a[1]) / denominator);
      }
    }
    return r;
  };
  for (int vertex = 0; vertex < serial.GetNV(); vertex++)
  {
    double *point = serial.GetVertex(vertex);
    const double lx = point[0] - c, lz = point[2] - c;
    const double norm = std::hypot(lx, lz);
    const double level = std::max(std::abs(lx), std::abs(lz));
    if (norm > 1.0e-12 && level < 3.0 - 1.0e-9)
    {
      const double ux = lx / norm, uz = lz / norm;
      const double square = 1.0 / std::max(std::abs(ux), std::abs(uz));
      const double polygon_scale = PolygonRadius(ux, uz) / square;
      const double blend = level <= 2.0 ? 1.0 : (3.0 - level);
      const double scale = 1.0 + blend * (polygon_scale - 1.0);
      point[0] = c + scale * lx;
      point[2] = c + scale * lz;
    }
  }
  serial.FinalizeTopology();
  serial.Finalize();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), serial);
}

}  // namespace

TEST_CASE("CornerRefinedRuleLayout", "[cornerbasisrefinement][Serial][Parallel]")
{
  const auto rule = RefinedCornerTraceBasisRule();
  CHECK(rule.AllRings());
  CHECK(rule.ring_size == 16);
  CHECK(CheckCornerTraceBasisRule(rule).empty());
  CHECK(CheckCornerTraceBasisRule(CornerTraceBasisRule{}).empty());
  {
    auto bad = rule;
    bad.ring_layout = CornerRingLayout::METAL_RINGS_ONLY;
    CHECK_THAT(CheckCornerTraceBasisRule(bad), ContainsSubstring("AllRingsFollowMetal"));
    bad = rule;
    bad.free_knot_grading = {1.0 / 3.0, 2.0 / 3.0, 1.0, 1.5, 2.0};
    CHECK_THAT(CheckCornerTraceBasisRule(bad), ContainsSubstring("more knots"));
    bad = rule;
    bad.free_knot_grading = {2.0 / 3.0, 1.0 / 3.0};
    CHECK_THAT(CheckCornerTraceBasisRule(bad), ContainsSubstring("increasing"));
    bad = rule;
    bad.ring_size = 15;
    CHECK_THAT(CheckCornerTraceBasisRule(bad), ContainsSubstring("RingSize"));
    bad = rule;
    bad.extra_levels_above_over_overetch = {4.0, 1.0};
    CHECK_THAT(CheckCornerTraceBasisRule(bad),
               ContainsSubstring("ExtraLevelsAboveOverOveretch"));
    CornerTraceBasisRule legacy_with_extra;
    legacy_with_extra.extra_levels_above_over_overetch = {1.0};
    CHECK_THAT(CheckCornerTraceBasisRule(legacy_with_extra),
               ContainsSubstring("AllRingsFollowMetal layout option"));
  }
  // The generator's layout at convex 120 and concave 105 degrees (R = 1.9): the perimeter
  // fractions per slot and the zero slots (generate_corner_response.rule_ring_layout with
  // REFINED_RULE; test_generate_corner_response.py pins the same numbers).
  const std::vector<double> convex_120 = {0.905502116982, 0.990696208596, 0.07589030021,
                                          0.161084391824, 0.246278483438, 0.331472575053,
                                          0.416666666667, 0.458333333333, 0.5,
                                          0.553694797275, 0.60738959455,  0.661084391824,
                                          0.714779189099, 0.768473986374, 0.822168783649,
                                          0.863835450315};
  const std::vector<double> concave_105 = {0.5,
                                           0.541666666667,
                                           0.583333333333,
                                           0.602804497065,
                                           0.622275660796,
                                           0.641746824527,
                                           0.661217988258,
                                           0.680689151989,
                                           0.700160315721,
                                           0.741826982387,
                                           0.783493649054,
                                           0.902911374212,
                                           0.022329099369,
                                           0.141746824527,
                                           0.261164549685,
                                           0.380582274842};
  for (const auto &[convex, angle, pins, zero_slots] :
       std::vector<std::tuple<bool, double, std::vector<double>, std::vector<int>>>{
           {true, 120.0, convex_120, {8, 9, 10, 11, 12, 13, 14}},
           {false, 105.0, concave_105, {0, 10, 11, 12, 13, 14, 15}}})
  {
    const auto layout = CornerRuleRingLayout(kR, angle * kDeg, convex, rule, true);
    std::map<int, CornerRingVertex> by_slot;
    int slaves = 0;
    for (const auto &vertex : layout)
    {
      if (vertex.kind == CornerRingVertex::Kind::SLAVE)
      {
        slaves++;
      }
      else
      {
        by_slot[vertex.slot] = vertex;
      }
    }
    REQUIRE(by_slot.size() == 16);
    CHECK(slaves == 4);
    for (int slot = 0; slot < 16; slot++)
    {
      CHECK_THAT(by_slot.at(slot).fraction, WithinAbs(pins[slot], 1.0e-9));
      const bool zero =
          std::find(zero_slots.begin(), zero_slots.end(), slot) != zero_slots.end();
      CHECK((by_slot.at(slot).kind == CornerRingVertex::Kind::ZERO) == zero);
    }
    CHECK(CornerZeroSlots(convex, rule) == zero_slots);
    // Off the metal every knot is free.
    for (const auto &vertex : CornerRuleRingLayout(kR, angle * kDeg, convex, rule, false))
    {
      CHECK(vertex.kind != CornerRingVertex::Kind::ZERO);
    }
  }
  // The constructed basis over the family range: 11 rings (the 7 standard levels, the
  // extra rings at t + oe and t + 4 oe, the two caps) x 16 knots, PEC on the two metal
  // rings only, the graded free knots at R/3 and 2R/3 from each crossing along the
  // perimeter, slaves with a partition of unity, cap centres at the mean of the crossing
  // knots, no degenerate triangle, every knot on its ring's square.
  for (const bool convex : {true, false})
  {
    for (double angle = 75.0; angle <= 180.0 + 1.0e-9; angle += 2.5)
    {
      const auto basis = BuildRefined(angle, convex);
      REQUIRE(basis.knots.size() == 176);
      REQUIRE(basis.contour_groups.size() == 11);
      int zero = 0;
      for (std::size_t k = 0; k < basis.knots.size(); k++)
      {
        zero += basis.zero[k] ? 1 : 0;
      }
      CHECK(zero ==
            14);  // 2 crossings + 5 metal-interior knots on each of the 2 metal rings
      std::set<double> levels;
      for (int r = 0; r < 9; r++)
      {
        levels.insert(std::round(basis.knots[16 * r][2] * 1.0e9) / 1.0e9);
      }
      CHECK(levels.count(std::round((kT + kOE) * 1.0e9) / 1.0e9) == 1);
      CHECK(levels.count(std::round((kT + 4.0 * kOE) * 1.0e9) / 1.0e9) == 1);
      CHECK(levels.count(std::round(-kOE * 1.0e9) / 1.0e9) == 1);
      const auto [first, second] = ArmCrossingFractions(kR, angle * kDeg);
      for (int r = 0; r < 11; r++)
      {
        const double half_width = r < 9 ? kR : kR / 3.0;
        std::vector<double> fractions;
        for (int i = 0; i < 16; i++)
        {
          const auto &p = basis.knots[16 * r + i];
          CHECK_THAT(std::max(std::abs(p[0]), std::abs(p[1])),
                     WithinAbs(half_width, 1.0e-9));
          fractions.push_back(SquarePerimeterFraction(half_width, p));
        }
        // A knot at R/3 and one at 2R/3 along the perimeter from each crossing.
        for (const double crossing : {first, second})
        {
          for (const double g : {1.0 / 3.0, 2.0 / 3.0})
          {
            const bool found =
                std::any_of(fractions.begin(), fractions.end(),
                            [&](double f)
                            {
                              const double d = std::abs(std::fmod(f - crossing + 3.0, 1.0));
                              return std::abs(std::min(d, 1.0 - d) - g / 8.0) < 1.0e-9;
                            });
            CHECK(found);
          }
        }
      }
      int centres = 0;
      for (std::size_t v = basis.knots.size(); v < basis.vertices.size(); v++)
      {
        const auto &vertex = basis.vertices[v];
        CHECK(vertex.basis < 0);
        CHECK(vertex.weight_a >= 0.0);
        CHECK(vertex.weight_a <= 1.0);
        if (std::abs(vertex.point[0]) < 1.0e-12 && std::abs(vertex.point[1]) < 1.0e-12)
        {
          centres++;
          CHECK_THAT(vertex.weight_a, WithinAbs(0.5, 1.0e-15));
          // Parents: the cap ring's crossing knots.
          for (const int parent : {vertex.parent_a, vertex.parent_b})
          {
            const auto &p = basis.knots[parent];
            CHECK_THAT(std::max(std::abs(p[0]), std::abs(p[1])),
                       WithinAbs(kR / 3.0, 1.0e-9));
            const double f = SquarePerimeterFraction(kR / 3.0, p);
            CHECK((std::abs(f - first) < 1.0e-9 ||
                   std::abs(std::fmod(f - second + 2.0, 1.0)) < 1.0e-9 ||
                   std::abs(std::fmod(f - second + 2.0, 1.0) - 1.0) < 1.0e-9));
          }
        }
      }
      CHECK(centres == 2);
      double min_area = mfem::infinity();
      for (const auto &t : basis.triangles)
      {
        min_area = std::min(min_area, TriangleArea(basis, t));
      }
      CHECK(min_area > 1.0e-8);
    }
  }
  // No connectivity angle on this layout.
  {
    const auto seed = MakeCornerBoxSeed(kR, kT, kOE, true, rule);
    CHECK_THROWS_WITH(BuildCornerTraceBasis(seed.points, seed.contour_groups,
                                            seed.zero_trace_indices, 100.0 * kDeg, true,
                                            rule, 112.5 * kDeg),
                      ContainsSubstring("no events"));
  }
}

TEST_CASE("CornerRefinedRuleNoEvents", "[cornerbasisrefinement][Serial][Parallel]")
{
  const auto rule = RefinedCornerTraceBasisRule();
  CHECK(CornerBasisEvents(true, rule).empty());
  CHECK(CornerBasisEvents(false, rule).empty());
  // Hats continuous in the angle across the MetalRingsOnly rule's event angles (its knots
  // passing box corners / side midpoints): the smooth trace's interpolant on the outer
  // faces changes by O(epsilon) only, while the MetalRingsOnly layout jumps O(0.1) at its
  // corner passages (the recorded band-quad flips).
  const CornerTraceBasisRule legacy;
  constexpr double epsilon = 1.0e-3;
  for (const auto &[convex, events] : std::vector<std::pair<bool, std::vector<double>>>{
           {true, {90.0, 135.0, 141.340191745910, 153.434948822922}},
           {false, {90.0, 111.801409486352, 135.0, 158.198590513648}}})
  {
    for (const double angle : events)
    {
      CHECK(MaxJumpAcross(rule, angle, convex, epsilon) < 1.0e-3);
    }
    // The control: the same probe on the recorded rule at its corner passages.
    CHECK(MaxJumpAcross(legacy, 90.0, convex, epsilon) > 1.0e-2);
    CHECK(MaxJumpAcross(legacy, 135.0, convex, epsilon) > 1.0e-2);
    // And the refined rule at a plain angle (the smooth-angle level).
    CHECK(MaxJumpAcross(rule, 127.5, convex, epsilon) < 1.0e-3);
  }
}

TEST_CASE("CornerRefinedRuleStencil", "[cornerbasisrefinement][Serial][Parallel]")
{
  const auto rule = RefinedCornerTraceBasisRule();
  constexpr double tol = 1.0e-2;
  std::vector<CornerFamilyNode> family;
  for (const double angle : {75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0})
  {
    family.push_back({angle, std::nullopt, family.size()});
  }
  for (const bool convex : {true, false})
  {
    CHECK(CheckCornerFamilySegments(family, convex, rule, tol).empty());
    // One segment: cubic on the four nodes around the angle everywhere, one-sided at the
    // range ends, exact at a node, no connectivity angle, refused outside the range.
    auto Angles = [&](const CornerFamilyStencil &stencil)
    {
      std::vector<double> angles;
      double sum = 0.0;
      for (const auto &[index, weight] : stencil.nodes)
      {
        angles.push_back(family[index].angle_degrees);
        sum += weight;
      }
      CHECK_THAT(sum, WithinAbs(1.0, 1.0e-12));
      CHECK(!stencil.connectivity_angle_degrees);
      return angles;
    };
    auto s82 = SelectCornerFamilyStencil(family, 82.5, convex, rule, tol);
    CHECK(s82.rule == "cubic");
    CHECK(Angles(s82) == std::vector<double>{75.0, 90.0, 105.0, 120.0});
    auto s100 = SelectCornerFamilyStencil(family, 100.0, convex, rule, tol);
    CHECK(s100.rule == "cubic");
    CHECK(Angles(s100) == std::vector<double>{75.0, 90.0, 105.0, 120.0});
    auto s127 = SelectCornerFamilyStencil(family, 127.5, convex, rule, tol);
    CHECK(s127.rule == "cubic");
    CHECK(Angles(s127) == std::vector<double>{105.0, 120.0, 135.0, 150.0});
    auto s176 = SelectCornerFamilyStencil(family, 176.0, convex, rule, tol);
    CHECK(s176.rule == "cubic");
    CHECK(Angles(s176) == std::vector<double>{135.0, 150.0, 165.0, 180.0});
    // Across the recorded rule's passages (135, 153.4 / 158.2) without a per-side coupon.
    auto s150 = SelectCornerFamilyStencil(family, 152.0, convex, rule, tol);
    CHECK(s150.rule == "cubic");
    CHECK(Angles(s150) == std::vector<double>{135.0, 150.0, 165.0, 180.0});
    auto exact = SelectCornerFamilyStencil(family, 135.0, convex, rule, tol);
    CHECK(exact.rule == "exact");
    CHECK(exact.nodes.size() == 1);
    CHECK(family[exact.base].angle_degrees == 135.0);
    CHECK_THAT(SelectCornerFamilyStencil(family, 70.0, convex, rule, tol).reason,
               ContainsSubstring("sharper"));
    // A coupon with a connectivity angle: refused (the layout has no events).
    auto stamped = family;
    stamped[2].connectivity_angle_degrees = 112.5;
    CHECK_THROWS_WITH(SelectCornerFamilyStencil(stamped, 100.0, convex, rule, tol),
                      ContainsSubstring("no events"));
    // Two coupons at one angle: a reason at load.
    auto duplicate = family;
    duplicate.push_back({105.0, std::nullopt, duplicate.size()});
    CHECK_THAT(CheckCornerFamilySegments(duplicate, convex, rule, tol),
               ContainsSubstring("two coupons at 105"));
  }
}

TEST_CASE("CornerRefinedRuleLibraryLoad", "[cornerbasisrefinement][Serial][Parallel]")
{
  // The load-time check of every coupon against the rule at its angle (ReadProcessLibrary,
  // fail closed): a family of the refined rule loads and its exact 90-degree corners match
  // the 90 node on a square island; a coupon whose basis points are the rule's at another
  // angle, a coupon carrying a connectivity angle and a coupon without the extra ring are
  // refused before any corner is matched. The runtime constructs the basis of an
  // interpolated corner (the steep house's 112.5-degree corners, cubic on 90 / 105 / 120 /
  // 135) from the rule and evaluates it under the SurfaceMortar lift (the centre slaves
  // enter the mortar mass through MortarVertex::ForEachBasis): every corner patch's
  // fabricated surface energy within a factor two of the exact 135-degree apex patch on
  // the same synthetic matrices.
  constexpr double R = 0.2, t = 0.01, oe = 0.005;
  test::SharedTempDir temp;
  const fs::path isolated_points = temp.temp_dir / "isolated-points.csv";
  const fs::path isolated_domain = temp.temp_dir / "isolated-domain.csv";
  const fs::path isolated_surface = temp.temp_dir / "isolated-surface.csv";
  const fs::path corner_domain = temp.temp_dir / "corner-domain.csv";
  const fs::path corner_surface = temp.temp_dir / "corner-surface.csv";
  std::map<std::string, fs::path> libraries;
  for (const std::string name : {"family", "wrong-angle", "connectivity", "no-extra-ring"})
  {
    libraries[name] = temp.temp_dir / ("library-" + name + ".json");
  }
  if (Mpi::Root(Mpi::World()))
  {
    {
      std::ofstream output(isolated_points);
      output << "x,y,z\n-0.16,-0.12,0.0\n0.16,-0.12,0.0\n0.16,0.12,0.0\n-0.16,0.12,0.0\n";
      std::ofstream domain(isolated_domain);
      domain << "basis_i,basis_j,Q_ij (J)\n";
      std::ofstream surface(isolated_surface);
      surface << "interface,edge,basis_i,basis_j,Q_total_ij (J)\n";
      for (int i = 1; i <= 4; i++)
      {
        for (int j = i; j <= 4; j++)
        {
          const double value = (i == j ? 2.0 : 0.2) * 1.0e-12;
          domain << i << "," << j << "," << value << "\n";
          surface << "1,1," << i << "," << j << "," << value << "\n";
        }
      }
      WriteSyntheticMatrices(corner_domain, corner_surface, 176, 3.0, 0.05, R);
    }
    const json base = {
        {"Version", 3},
        {"TraceLiftVersion", 2},
        {"MatchingRadius", R},
        {"Fabrication",
         {{"MetalThickness", t},
          {"OveretchDepth", oe},
          {"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
        {"Models",
         {{{"Name", "isolated"},
           {"Topology", "IsolatedEdge"},
           {"CouponDepth", R},
           {"FabricatedMatrix", isolated_domain.string()},
           {"ThinMatrix", isolated_domain.string()},
           {"FabricatedSurfaceMatrix", isolated_surface.string()},
           {"ThinSurfaceMatrix", isolated_surface.string()},
           {"BasisPoints", isolated_points.string()},
           {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}}}}};
    auto Model = [&](const std::string &tag, double angle, const CouponFiles &files)
    {
      return json{{"Name", "convex-corner-" + tag},
                  {"Topology", "ConvexCorner"},
                  {"Angle", angle},
                  {"AngleDegrees", angle},
                  {"Convexity", "Convex"},
                  {"AngleTolerance", 1.0e-6},
                  {"CornerRadius", 0.0},
                  {"CornerRadiusTolerance", 0.0},
                  {"FabricatedMatrix", corner_domain.string()},
                  {"ThinMatrix", corner_domain.string()},
                  {"FabricatedSurfaceMatrix", corner_surface.string()},
                  {"ThinSurfaceMatrix", corner_surface.string()},
                  {"BasisPoints", files.points.string()},
                  {"TraceMesh",
                   {{"Vertices", files.vertices.string()},
                    {"Triangles", files.triangles.string()}}},
                  {"ContourGroups", files.contour_groups},
                  {"ZeroTraceIndices", files.zero_trace_indices},
                  {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}},
                  {"TraceBasis", RefinedTraceBasisRecord()}};
    };
    auto Write = [&](const std::string &name, const std::vector<json> &corners)
    {
      json library = base;
      library["Name"] = "unit-test-corner-refined-" + name;
      for (const auto &corner : corners)
      {
        library["Models"].push_back(corner);
      }
      std::ofstream output(libraries.at(name));
      output << library.dump(2) << "\n";
    };
    std::vector<json> family;
    for (const double angle : {90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0})
    {
      const std::string tag = std::to_string(static_cast<int>(angle));
      family.push_back(
          Model(tag, angle, WriteRefinedCoupon(temp.temp_dir, tag, angle, true, R, t, oe)));
    }
    Write("family", family);
    // A free knot of the 90-degree coupon moved 1e-6 R along its ring (the crossings gate
    // still passes: the crossings are PEC knots and no free knot lies on the metal).
    {
      auto wrong = family;
      wrong[0] = Model(
          "90", 90.0,
          WriteRefinedCoupon(temp.temp_dir, "90-wrong", 90.0, true, R, t, oe, 6 * 16 + 0));
      Write("wrong-angle", wrong);
    }
    {
      auto stamped = family;
      stamped[1]["TraceBasis"]["ConnectivityAngleDegrees"] = 112.5;
      Write("connectivity", stamped);
    }
    {
      // A coupon of the refined rule's fractions without the ring at t + 4 oe (the seed of
      // a rule with the t + oe ring only): the level set is not the rule's.
      auto rule = RefinedCornerTraceBasisRule();
      auto one_extra = rule;
      one_extra.extra_levels_above_over_overetch = {1.0};
      const auto seed = MakeCornerBoxSeed(R, t, oe, true, one_extra);
      const auto basis =
          BuildCornerTraceBasis(seed.points, seed.contour_groups, seed.zero_trace_indices,
                                105.0 * kDeg, true, rule);
      CouponFiles files;
      files.points = temp.temp_dir / "refined-105-seven-points.csv";
      files.vertices = temp.temp_dir / "refined-105-seven-vertices.csv";
      files.triangles = temp.temp_dir / "refined-105-seven-triangles.csv";
      files.contour_groups = seed.contour_groups;
      std::ofstream points(files.points), vertices(files.vertices),
          triangles(files.triangles);
      points << std::setprecision(17) << "x,y,z\n";
      for (std::size_t k = 0; k < basis.knots.size(); k++)
      {
        points << basis.knots[k][0] << "," << basis.knots[k][1] << "," << basis.knots[k][2]
               << "\n";
        if (basis.zero[k])
        {
          files.zero_trace_indices.push_back(static_cast<int>(k) + 1);
        }
      }
      vertices << "vertex,x,y,z,basis,conductor,parent_a,parent_b,weight_a\n";
      for (std::size_t v = 0; v < basis.vertices.size(); v++)
      {
        const auto &vertex = basis.vertices[v];
        vertices << std::setprecision(17) << v + 1 << "," << vertex.point[0] << ","
                 << vertex.point[1] << "," << vertex.point[2] << ",";
        if (vertex.basis >= 0)
        {
          vertices << vertex.basis + 1 << "," << (basis.zero[vertex.basis] ? 1 : 0)
                   << ",0,0,0\n";
        }
        else
        {
          vertices << "0,0," << vertex.parent_a + 1 << "," << vertex.parent_b + 1 << ","
                   << vertex.weight_a << "\n";
        }
      }
      triangles << "triangle,vertex_i,vertex_j,vertex_k\n";
      for (std::size_t tri = 0; tri < basis.triangles.size(); tri++)
      {
        triangles << tri + 1 << "," << basis.triangles[tri][0] + 1 << ","
                  << basis.triangles[tri][1] + 1 << "," << basis.triangles[tri][2] + 1
                  << "\n";
      }
      auto missing = family;
      missing[1] = Model("105", 105.0, files);
      Write("no-extra-ring", missing);
    }
  }
  Mpi::Barrier(Mpi::World());

  json config = {
      {"Problem", {{"Type", "Electrostatic"}, {"Output", temp.temp_dir.string()}}},
      {"Model", {{"Mesh", "unused.msh"}}},
      {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
      {"Boundaries",
       {{"Ground", {{"Attributes", {1, 2, 3, 4, 5, 6}}}},
        {"Terminal", {{{"Index", 1}, {"Attributes", {9}}}}},
        {"Postprocessing",
         {{"Dielectric",
           {{{"Index", 4},
             {"Attributes", {9}},
             {"Type", "SA"},
             {"Thickness", 0.002},
             {"Permittivity", 4.0},
             {"AutomaticEdges", true},
             {"EdgeDistances", {R}},
             {"EdgeFrameNormal", {0.0, 1.0, 0.0}}}}}}}}},
      {"Solver",
       {{"Order", 1},
        {"Electrostatic",
         {{"ResponseCorrection",
           {{"Library", ""},
            {"TargetInterfaces", {4}},
            {"UnmatchedPolicy", "Error"},
            {"TraceCoupling", "SurfaceMortar"},
            {"MortarOversampling", 2}}}}}}}};
  // The square island (four exact 90-degree corners) and the steep house (two 112.5-degree
  // corners, a 135-degree apex, two 90-degree corners; vertices on grid rays) in the
  // centred unit frame of MakePolygonIslandMesh.
  const std::vector<std::array<double, 2>> square = {
      {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
  const std::vector<std::array<double, 2>> steep_house = {
      {-1.0, -1.0},
      {1.0, -1.0},
      {1.0, 0.2},
      {0.0, 0.2 + std::tan(22.5 * M_PI / 180.0)},
      {-1.0, 0.2}};
  auto square_mesh = MakePolygonIslandMesh(square, 8.0, 0.1);
  auto steep_mesh = MakePolygonIslandMesh(steep_house, 8.0, 0.1);
  auto ConfigFor = [&](const std::string &name)
  {
    auto result = config;
    result["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        libraries.at(name).string();
    return result;
  };
  auto Preflight = [&](const std::string &name, mfem::ParMesh &mesh)
  {
    IoData iodata(ConfigFor(name), false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto manifest_path = temp.temp_dir / ("requirements-" + name + ".json");
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };
  CHECK_THROWS_WITH(Preflight("wrong-angle", *square_mesh),
                    ContainsSubstring("is not its trace basis rule's coupon") &&
                        ContainsSubstring("basis point"));
  CHECK_THROWS_WITH(Preflight("connectivity", *square_mesh),
                    ContainsSubstring("ConnectivityAngleDegrees") &&
                        ContainsSubstring("no events"));
  CHECK_THROWS_WITH(Preflight("no-extra-ring", *square_mesh),
                    ContainsSubstring("is not its trace basis rule's coupon") &&
                        ContainsSubstring("outer ring levels"));
  {
    const json manifest = Preflight("family", *square_mesh);
    int corners = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner")
      {
        corners++;
        CHECK(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"] == "convex-corner-90");
      }
    }
    CHECK(corners == 4);
  }
  {
    const json manifest = Preflight("family", *steep_mesh);
    std::map<std::string, int> matched;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConvexCorner")
      {
        REQUIRE(feature["Match"]["Status"] == "Matched");
        matched[feature["Match"]["Model"].get<std::string>()]++;
      }
    }
    CHECK(matched["convex-corner-90"] == 2);
    CHECK(matched["convex-corner-135"] == 1);
    int interpolated = 0;
    for (const auto &[name, count] : matched)
    {
      if (name.find("@corner-angle112.5-cubic") != std::string::npos)
      {
        interpolated += count;
      }
    }
    CHECK(interpolated == 2);
    // The runtime: the constructed 112.5-degree basis under the SurfaceMortar lift.
    IoData iodata(ConfigFor("family"), false);
    iodata.boundaries.cracked_attributes.insert(9);
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(*steep_mesh)));
    LaplaceOperator laplace(iodata, meshes);
    SurfaceResponseOperator response(iodata, laplace);
    mfem::ParGridFunction potential(&laplace.GetH1Space().Get());
    mfem::FunctionCoefficient potential_coefficient(
        [](const mfem::Vector &x)
        { return (x[1] - 0.5) * (1.0 + 0.1 * (x[0] - 4.0) - 0.05 * (x[2] - 4.0)); });
    potential.ProjectCoefficient(potential_coefficient);
    Vector potential_true;
    potential.GetTrueDofs(potential_true);
    const auto result = response.GetElectrostaticResponse(potential_true);
    const auto &names = response.GetModelNames();
    std::map<std::string, double> per_patch;
    for (const auto &contribution : result.model_contributions)
    {
      const auto &name = names.at(contribution.model);
      if (name.find("corner") != std::string::npos)
      {
        REQUIRE(contribution.patch_count > 0.0);
        per_patch[name] =
            contribution.fabricated_surface_energy.at(4) / contribution.patch_count;
        CHECK(std::isfinite(per_patch[name]));
        CHECK(per_patch[name] > 0.0);
      }
    }
    REQUIRE(per_patch.size() == 3);  // 90, 135 and the constructed 112.5
    const double apex = per_patch.at("convex-corner-135");
    for (const auto &[name, energy] : per_patch)
    {
      CHECK(energy > 0.5 * apex);
      CHECK(energy < 2.0 * apex);
    }
  }
}

}  // namespace palace
