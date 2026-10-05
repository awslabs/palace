// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// The corner family's REFINED trace basis (corner-basis refinement, USER decision 161 (1),
// 2026-09-30; design doc SURFACE-RESPONSE-IDENTIFICATION.md, Conventions
// CornerTraceBasisRule): the AllRingsFollowMetal layout — every ring of the box, the extra
// rings at MetalThickness + OveretchDepth and MetalThickness + 4 OveretchDepth and the two
// cap rings included, carries the same
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
#include <sstream>
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
CouponFiles
WriteRefinedCoupon(const fs::path &directory, const std::string &tag, double angle,
                   bool convex, double radius, double t, double oe,
                   std::optional<int> perturb_knot = std::nullopt,
                   const CornerTraceBasisRule &rule = RefinedCornerTraceBasisRule())
{
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

// MakePolygonIslandMesh for a STAR-SHAPED (non-convex) polygon: the radial function is the
// ray's hit on the boundary SEGMENT (the half-plane form above is the convex hull of the
// edge lines), so an island with a notch — a concave metal corner whose angle is the
// notch's opening — is exact.
std::unique_ptr<mfem::ParMesh>
MakeStarIslandMesh(const std::vector<std::array<double, 2>> &polygon, double extent,
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
  auto StarRadius = [&](double ux, double uz)
  {
    double r = mfem::infinity();
    for (std::size_t k = 0; k < polygon.size(); k++)
    {
      const auto &a = polygon[k], &b = polygon[(k + 1) % polygon.size()];
      // Ray (0, 0) + s (ux, uz) against the segment a + t (b - a), t in [0, 1].
      const double dx = b[0] - a[0], dz = b[1] - a[1];
      const double determinant = ux * dz - uz * dx;
      if (std::abs(determinant) < 1.0e-14)
      {
        continue;
      }
      const double s = (a[0] * dz - a[1] * dx) / determinant;
      const double t = (a[0] * uz - a[1] * ux) / determinant;
      if (s > 0.0 && t >= -1.0e-12 && t <= 1.0 + 1.0e-12)
      {
        r = std::min(r, s);
      }
    }
    MFEM_VERIFY(std::isfinite(r), "The star-shaped island polygon misses a ray!");
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
      const double polygon_scale = StarRadius(ux, uz) / square;
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

// The unit-square island with a notch of opening `angle` degrees cut from the top side down
// to the apex (0, 1 - 0.3 / tan(angle / 2)) between the grid-ray vertices (+-0.3, 1): one
// CONCAVE metal corner of that angle, two convex corners at the notch mouth.
std::vector<std::array<double, 2>> NotchedIsland(double angle)
{
  const double apex = 1.0 - 0.3 / std::tan(0.5 * angle * kDeg);
  return {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {0.3, 1.0},
          {0.0, apex},  {-0.3, 1.0}, {-1.0, 1.0}};
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

TEST_CASE("CornerRefinedRuleAcuteConcaveGrading",
          "[cornerbasisrefinement][Serial][Parallel]")
{
  // Block (b) family 4 (DESIGN A9; decisions 304 / 318): the graded free-knot distances
  // scale by min(1, FreeArc / Reference), Reference = the concave 60-degree node's free arc
  // 2 - cot 60 deg, so acute concave nodes build with that node's layout proportions (the
  // generator refused them below 56.3 degrees: no inner interval between the two 2R / 3
  // knots) while every node whose free arc reaches the reference is bit-identical.
  const auto rule = RefinedCornerTraceBasisRule();
  CHECK(rule.free_knot_grading_reference_free_arc_over_r ==
        kFreeKnotGradingReferenceFreeArcOverR);
  CHECK(kFreeKnotGradingReferenceFreeArcOverR == 1.4226497308103743);
  CHECK_THAT(kFreeKnotGradingReferenceFreeArcOverR,
             WithinAbs(2.0 - 1.0 / std::sqrt(3.0), 1.0e-15));
  CHECK(CheckCornerTraceBasisRule(rule).empty());
  {
    auto bad = rule;
    bad.free_knot_grading_reference_free_arc_over_r = 4.0 / 3.0;  // = 2 x the outermost
    CHECK_THAT(CheckCornerTraceBasisRule(bad),
               ContainsSubstring("FreeKnotGradingReferenceFreeArcOverR"));
    // The default rule (no grading) ignores the reference.
    CornerTraceBasisRule legacy;
    legacy.free_knot_grading_reference_free_arc_over_r = 0.1;
    CHECK(CheckCornerTraceBasisRule(legacy).empty());
  }
  auto FreeArc = [](double angle, bool convex)
  {
    const auto [first, second] = ArmCrossingFractions(kR, angle * kDeg);
    // The free arc: the sector between the arms for a concave corner, its complement for a
    // convex one (the generator's metal_arc_fractions), over R.
    return 8.0 * (convex ? (first + 1.0) - second : second - first);
  };
  // The concave free arc is (2 - cot theta) R: the 60-degree node's own arc (from its
  // crossing fractions) is the reference within the tolerance -> scale exactly 1; the
  // unscaled grading needs > 4R / 3 (58 ok / 56 fail of the group-A lane).
  CHECK_THAT(FreeArc(60.0, false), WithinAbs(kFreeKnotGradingReferenceFreeArcOverR,
                                             kFreeKnotGradingReferenceToleranceOverR));
  CHECK(FreeKnotGradingScale(FreeArc(60.0, false), rule) == 1.0);
  CHECK(FreeArc(58.0, false) > 4.0 / 3.0);
  CHECK(FreeArc(56.0, false) < 4.0 / 3.0);
  CHECK_THAT(FreeArc(45.0, false), WithinAbs(1.0, 1.0e-12));
  CHECK_THAT(FreeKnotGradingScale(FreeArc(45.0, false), rule),
             WithinAbs(1.0 / kFreeKnotGradingReferenceFreeArcOverR, 1.0e-15));
  CHECK(FreeKnotGradingScale(FreeArc(52.5, false), rule) ==
        FreeArc(52.5, false) / kFreeKnotGradingReferenceFreeArcOverR);
  CHECK(FreeKnotGradingScale(FreeArc(15.0, true), rule) == 1.0);
  CHECK(FreeKnotGradingScale(FreeArc(180.0, true), rule) == 1.0);
  // Pins shared with test_generate_corner_response.AcuteConcaveGradingTest (R = 1.9).
  const std::vector<double> concave_45 = {
      0.5,           0.529288071241, 0.558576142482, 0.559884094988, 0.561192047494,
      0.5625,        0.563807952506, 0.565115905012, 0.566423857518, 0.595711928759,
      0.625,         0.770833333333, 0.916666666667, 0.0625,         0.208333333333,
      0.354166666667};
  const std::vector<double> concave_52p5 = {0.5,
                                            0.536102614993,
                                            0.572205229985,
                                            0.573817507741,
                                            0.575429785496,
                                            0.577042063251,
                                            0.578654341007,
                                            0.580266618762,
                                            0.581878896517,
                                            0.61798151151,
                                            0.654084126503,
                                            0.795070105419,
                                            0.936056084335,
                                            0.077042063251,
                                            0.218028042168,
                                            0.359014021084};
  const std::vector<int> concave_zero_slots = {0, 10, 11, 12, 13, 14, 15};
  for (const auto &[angle, pins, slaves] :
       std::vector<std::tuple<double, std::vector<double>, int>>{{45.0, concave_45, 3},
                                                                 {52.5, concave_52p5, 4}})
  {
    const auto layout = CornerRuleRingLayout(kR, angle * kDeg, false, rule, true);
    std::map<int, CornerRingVertex> by_slot;
    int slave_count = 0;
    for (const auto &vertex : layout)
    {
      if (vertex.kind == CornerRingVertex::Kind::SLAVE)
      {
        slave_count++;
      }
      else
      {
        by_slot[vertex.slot] = vertex;
      }
    }
    REQUIRE(by_slot.size() == 16);
    CHECK(slave_count == slaves);
    for (int slot = 0; slot < 16; slot++)
    {
      CHECK_THAT(by_slot.at(slot).fraction, WithinAbs(pins[slot], 1.0e-9));
      const bool zero = std::find(concave_zero_slots.begin(), concave_zero_slots.end(),
                                  slot) != concave_zero_slots.end();
      CHECK((by_slot.at(slot).kind == CornerRingVertex::Kind::ZERO) == zero);
    }
  }
  // Arc-relative free-knot positions of a scaled node = the 60-degree node's (slot by
  // slot).
  auto Relative = [&](double angle)
  {
    const auto [first, second] = ArmCrossingFractions(kR, angle * kDeg);
    std::map<int, double> relative;
    for (const auto &vertex : CornerRuleRingLayout(kR, angle * kDeg, false, rule, true))
    {
      if (vertex.kind == CornerRingVertex::Kind::FREE)
      {
        const double unwrapped =
            vertex.fraction < first - 1.0e-12 ? vertex.fraction + 1.0 : vertex.fraction;
        relative[vertex.slot] = (unwrapped - first) / (second - first);
      }
    }
    return relative;
  };
  const auto reference = Relative(60.0);
  REQUIRE(reference.size() == 9);
  for (const double angle : {45.0, 46.809687, 51.672903, 52.5, 56.0, 58.0})
  {
    const auto relative = Relative(angle);
    REQUIRE(relative.size() == 9);
    for (const auto &[slot, value] : reference)
    {
      CHECK_THAT(relative.at(slot), WithinAbs(value, 1.0e-12));
    }
  }
  // Bit-identity where the scaling is inactive: the layout under the default reference
  // equals the layout under a reference just above the unscaled minimum (both scale 1, the
  // same floating-point expressions) for every concave node >= 60 degrees and every convex
  // node; it differs at 45 degrees.
  auto unscaled = rule;
  unscaled.free_knot_grading_reference_free_arc_over_r = 4.0 / 3.0 + 1.0e-6;
  REQUIRE(CheckCornerTraceBasisRule(unscaled).empty());
  auto Fractions = [&](double angle, bool convex, const CornerTraceBasisRule &r)
  {
    std::vector<double> fractions;
    for (const auto &vertex : CornerRuleRingLayout(kR, angle * kDeg, convex, r, true))
    {
      fractions.push_back(vertex.fraction);
    }
    return fractions;
  };
  for (const bool convex : {true, false})
  {
    for (double angle = convex ? 15.0 : 60.0; angle <= 180.0 + 1.0e-9; angle += 2.5)
    {
      CHECK(Fractions(angle, convex, rule) == Fractions(angle, convex, unscaled));
    }
    for (const double angle : {70.498378, 63.53, 82.5, 112.5, 142.5, 176.0})
    {
      CHECK(Fractions(angle, convex, rule) == Fractions(angle, convex, unscaled));
    }
  }
  CHECK(Fractions(45.0, false, rule) != Fractions(58.0, false, rule));
  CHECK_THROWS_WITH(Fractions(45.0, false, unscaled), ContainsSubstring("too short"));
  // The constructed basis of the acute nodes and of C3's two keys: 176 knots, 14 PEC, every
  // knot on its ring's square, no degenerate triangle.
  for (const double angle : {45.0, 46.809687, 51.672903, 52.5, 56.0})
  {
    const auto basis = BuildRefined(angle, false);
    REQUIRE(basis.knots.size() == 176);
    int zero = 0;
    for (std::size_t k = 0; k < basis.knots.size(); k++)
    {
      zero += basis.zero[k] ? 1 : 0;
      const double half_width = k < 9 * 16 ? kR : kR / 3.0;
      CHECK_THAT(std::max(std::abs(basis.knots[k][0]), std::abs(basis.knots[k][1])),
                 WithinAbs(half_width, 1.0e-9));
    }
    CHECK(zero == 14);
    double min_area = mfem::infinity();
    for (const auto &t : basis.triangles)
    {
      min_area = std::min(min_area, TriangleArea(basis, t));
    }
    CHECK(min_area > 1.0e-8);
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

TEST_CASE("CornerRefinedRuleAcuteConcaveLibraryLoad",
          "[cornerbasisrefinement][Serial][Parallel]")
{
  // Family 4 at library load and match time (decision 318): a concave family whose
  // 45-degree node records FreeKnotGradingReferenceFreeArcOverR (the scaling is active
  // there) beside 60 / 75 / 90-degree nodes WITHOUT the key (the default reference) is ONE
  // rule: it loads, every node's basis points pass the rule check at its own angle, a
  // 45-degree notch matches the 45 node exactly and a 52.5-degree notch is interpolated
  // (cubic on 45 / 60 / 75 / 90, the constructed basis scaled at 52.5). A node recording
  // another reference is refused ("not built on one TraceBasis rule"); the key on a
  // MetalRingsOnly record is refused at parse. R = 0.05 on the unit island (h = 0.1): the
  // notch arms are 15.7 R long and 12 R apart at the mouth, so the apex is a ConcaveCorner
  // feature (with R = 0.2 the two arms lie within the cluster radius and the notch is one
  // SpatialEdgeCluster).
  constexpr double R = 0.02, t = 0.001, oe = 0.0005;
  test::SharedTempDir temp;
  const fs::path isolated_points = temp.temp_dir / "isolated-points.csv";
  const fs::path isolated_domain = temp.temp_dir / "isolated-domain.csv";
  const fs::path isolated_surface = temp.temp_dir / "isolated-surface.csv";
  const fs::path corner_domain = temp.temp_dir / "corner-domain.csv";
  const fs::path corner_surface = temp.temp_dir / "corner-surface.csv";
  std::map<std::string, fs::path> libraries;
  for (const std::string name : {"family", "other-reference", "legacy-reference"})
  {
    libraries[name] = temp.temp_dir / ("library-acute-" + name + ".json");
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
    auto Model = [&](bool convex, double angle, const CouponFiles &files, json record)
    {
      const std::string topology = convex ? "convex" : "concave";
      std::ostringstream tag;
      tag << angle;
      return json{{"Name", topology + "-corner-" + tag.str()},
                  {"Topology", convex ? "ConvexCorner" : "ConcaveCorner"},
                  {"Angle", angle},
                  {"AngleDegrees", angle},
                  {"Convexity", convex ? "Convex" : "Concave"},
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
                  {"TraceBasis", record}};
    };
    auto Write = [&](const std::string &name, const std::vector<json> &corners)
    {
      json library = base;
      library["Name"] = "unit-test-corner-acute-concave-" + name;
      for (const auto &corner : corners)
      {
        library["Models"].push_back(corner);
      }
      std::ofstream output(libraries.at(name));
      output << library.dump(2) << "\n";
    };
    // The generator writes the reference key only on the node where the scaling is active.
    json active = RefinedTraceBasisRecord();
    active["FreeKnotGradingReferenceFreeArcOverR"] = kFreeKnotGradingReferenceFreeArcOverR;
    std::vector<json> family;
    for (const double angle : {90.0, 105.0, 120.0, 135.0})
    {
      std::ostringstream tag;
      tag << "acute-convex-" << angle;
      family.push_back(Model(
          true, angle, WriteRefinedCoupon(temp.temp_dir, tag.str(), angle, true, R, t, oe),
          RefinedTraceBasisRecord()));
    }
    for (const double angle : {45.0, 60.0, 75.0, 90.0})
    {
      std::ostringstream tag;
      tag << "acute-concave-" << angle;
      family.push_back(
          Model(false, angle,
                WriteRefinedCoupon(temp.temp_dir, tag.str(), angle, false, R, t, oe),
                angle < 60.0 ? active : RefinedTraceBasisRecord()));
    }
    Write("family", family);
    {
      // A 45-degree node built AND recorded with another reference: its own files pass the
      // per-coupon rule check; the family refuses the second rule.
      auto other_rule = RefinedCornerTraceBasisRule();
      other_rule.free_knot_grading_reference_free_arc_over_r = 1.5;
      json other_record = active;
      other_record["FreeKnotGradingReferenceFreeArcOverR"] = 1.5;
      auto other = family;
      other[4] = Model(false, 45.0,
                       WriteRefinedCoupon(temp.temp_dir, "acute-concave-45-other", 45.0,
                                          false, R, t, oe, std::nullopt, other_rule),
                       other_record);
      Write("other-reference", other);
    }
    {
      auto legacy = family;
      legacy[4]["TraceBasis"] = {{"RingSize", 8},
                                 {"MetalInteriorKnots", 1},
                                 {"FreeKnots", 5},
                                 {"Fractions", "PerimeterArcLength"},
                                 {"FreeKnotGradingReferenceFreeArcOverR", 1.5}};
      Write("legacy-reference", legacy);
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
  auto notch_45 = MakeStarIslandMesh(NotchedIsland(45.0), 8.0, 0.1);
  auto notch_52p5 = MakeStarIslandMesh(NotchedIsland(52.5), 8.0, 0.1);
  auto Preflight = [&](const std::string &name, mfem::ParMesh &mesh)
  {
    auto result = config;
    result["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] =
        libraries.at(name).string();
    IoData iodata(result, false);
    iodata.boundaries.cracked_attributes.insert(9);
    const auto manifest_path = temp.temp_dir / ("requirements-acute-" + name + ".json");
    WriteSurfaceResponseRequirements(iodata, mesh, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    return json::parse(input);
  };
  // The exact 45-degree notch resolves by the signature key (no family needed); the
  // interpolated 52.5-degree notch assembles the family and refuses the second rule.
  CHECK_THROWS_WITH(Preflight("other-reference", *notch_52p5),
                    ContainsSubstring("not built on one TraceBasis rule"));
  CHECK_THROWS_WITH(Preflight("legacy-reference", *notch_45),
                    ContainsSubstring("FreeKnotGradingReferenceFreeArcOverR"));
  auto ConcaveMatches = [](const json &manifest)
  {
    std::map<std::string, int> matched;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConcaveCorner")
      {
        REQUIRE(feature["Match"]["Status"] == "Matched");
        matched[feature["Match"]["Model"].get<std::string>()]++;
      }
    }
    return matched;
  };
  {
    const auto matched = ConcaveMatches(Preflight("family", *notch_45));
    REQUIRE(matched.size() == 1);
    CHECK(matched.at("concave-corner-45") == 1);
  }
  {
    const auto matched = ConcaveMatches(Preflight("family", *notch_52p5));
    REQUIRE(matched.size() == 1);
    const std::string name = matched.begin()->first;
    CHECK(matched.begin()->second == 1);
    CHECK(name.find("@corner-angle52.5-cubic") != std::string::npos);
  }
}

TEST_CASE("CornerRefinedRuleAcuteConcaveDeviceMortar",
          "[cornerbasisrefinement][Serial][Parallel]")
{
  // Family 4 round 2 (ii), decision 325: what the DEVICE runtime does with a constructed
  // (interpolated, device-angle) basis whose hats are far below the device mesh resolution.
  // At 48.75 degrees the scaled layout's five inner free knots on the cap rings (|x|, |y|
  // <= R/3 at z = +-R) are 0.0039 R apart (7.4 nm at R = 1.9 um), so a hat's support is a
  // sliver of width 0.0078 R; on this device mesh (hexahedra of h = 0.05 = 2.5 R = 640
  // sliver gaps: the WHOLE matching box is smaller than one element) no sliver hat
  // contains a device DOF. The runtime does not interpolate the
  // device potential at DOF nodes: the SurfaceMortar lift is the L2 projection onto the
  // P1 hats of the coupon's trace triangulation, with the device field (defined everywhere
  // as an FE function) sampled by FindPointsGSLIB at a 2 x 2 Gauss rule on every trace
  // triangle (subdivided by the local mesh size / MortarOversampling), relative to the
  // patch's conductor reference. The rule is exact for a field linear on the box, so every
  // free hat — sliver or not — must read the field's value at its knot (relative to the
  // metal) to rounding, no coefficient is zero or non-finite, and the energies are finite:
  // the device side neither aborts nor zeroes nor silently misreads a sliver hat. (The
  // coupon side is the opposite: the coupon's electrostatic solve prescribes the hat by
  // NODAL interpolation onto its boundary DOFs, so a sliver hat with no DOF in its support
  // is a zero response row, refused by finalize_corner_response.py's positive-definiteness
  // check — the 48.75-degree held-out coupon of decision 325.)
  constexpr double R = 0.02, t = 0.001, oe = 0.0005;
  constexpr double angle_device = 48.75;
  test::SharedTempDir temp;
  const fs::path isolated_points = temp.temp_dir / "isolated-points.csv";
  const fs::path isolated_domain = temp.temp_dir / "isolated-domain.csv";
  const fs::path isolated_surface = temp.temp_dir / "isolated-surface.csv";
  const fs::path corner_domain = temp.temp_dir / "corner-domain.csv";
  const fs::path corner_surface = temp.temp_dir / "corner-surface.csv";
  const fs::path library = temp.temp_dir / "library-acute-device-mortar.json";
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
    json active = RefinedTraceBasisRecord();
    active["FreeKnotGradingReferenceFreeArcOverR"] = kFreeKnotGradingReferenceFreeArcOverR;
    json models = json::array();
    models.push_back({{"Name", "isolated"},
                      {"Topology", "IsolatedEdge"},
                      {"CouponDepth", R},
                      {"FabricatedMatrix", isolated_domain.string()},
                      {"ThinMatrix", isolated_domain.string()},
                      {"FabricatedSurfaceMatrix", isolated_surface.string()},
                      {"ThinSurfaceMatrix", isolated_surface.string()},
                      {"BasisPoints", isolated_points.string()},
                      {"Interfaces", {{{"Type", "SA"}, {"Coupon", 1}}}}});
    auto Model = [&](bool convex, double angle)
    {
      std::ostringstream tag;
      tag << (convex ? "convex" : "concave") << "-" << angle;
      const auto files = WriteRefinedCoupon(temp.temp_dir, "device-mortar-" + tag.str(),
                                            angle, convex, R, t, oe);
      // The generator writes the reference key only where the scaling is active.
      const bool scaled = !convex && angle < 60.0;
      return json{{"Name", tag.str() + "-corner"},
                  {"Topology", convex ? "ConvexCorner" : "ConcaveCorner"},
                  {"Angle", angle},
                  {"AngleDegrees", angle},
                  {"Convexity", convex ? "Convex" : "Concave"},
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
                  {"TraceBasis", scaled ? active : RefinedTraceBasisRecord()}};
    };
    // The notch mouth's two convex corners (114.4 degrees) interpolate on the convex nodes;
    // the 48.75-degree apex interpolates on the family-4 segment 45 / 52.5 / 60 / 75.
    for (const double angle : {90.0, 105.0, 120.0, 135.0})
    {
      models.push_back(Model(true, angle));
    }
    for (const double angle : {45.0, 52.5, 60.0, 75.0, 90.0})
    {
      models.push_back(Model(false, angle));
    }
    const json record = {
        {"Version", 3},
        {"TraceLiftVersion", 2},
        {"Name", "unit-test-corner-acute-concave-device-mortar"},
        {"MatchingRadius", R},
        {"Fabrication",
         {{"MetalThickness", t},
          {"OveretchDepth", oe},
          {"InterfaceLayers", {{"SA", {{"Thickness", 0.002}, {"Permittivity", 4.0}}}}}}},
        {"Models", models}};
    std::ofstream output(library);
    output << record.dump(2) << "\n";
  }
  Mpi::Barrier(Mpi::World());

  const json config = {
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
           {{"Library", library.string()},
            {"TargetInterfaces", {4}},
            // The notch arms at this opening are three SpatialEdgeCluster features the
            // library has no model for (omitted with a warning); the corners are the test.
            {"UnmatchedPolicy", "Warn"},
            {"TraceCoupling", "SurfaceMortar"},
            {"MortarOversampling", 2}}}}}}}};
  auto notch = MakeStarIslandMesh(NotchedIsland(angle_device), 8.0, 0.05);
  IoData iodata(config, false);
  iodata.boundaries.cracked_attributes.insert(9);
  {
    const auto manifest_path = temp.temp_dir / "requirements-acute-device-mortar.json";
    WriteSurfaceResponseRequirements(iodata, *notch, manifest_path.string());
    Mpi::Barrier(Mpi::World());
    std::ifstream input(manifest_path);
    REQUIRE(input);
    const json manifest = json::parse(input);
    int concave = 0;
    for (const auto &feature : manifest["Identification"]["Features"])
    {
      if (feature["Type"] == "ConcaveCorner")
      {
        concave++;
        REQUIRE(feature["Match"]["Status"] == "Matched");
        CHECK(feature["Match"]["Model"].get<std::string>().find(
                  "@corner-angle48.75-cubic") != std::string::npos);
      }
    }
    CHECK(concave == 1);
  }

  // The constructed basis the runtime builds at the device angle: the sliver hats are the
  // five inner free knots of every ring (the 48.75-degree layout, scaled 0.789).
  const auto constructed = BuildRefined(angle_device, false, R, t, oe);
  REQUIRE(constructed.knots.size() == 176);
  std::vector<int> cap_sliver_knots;
  double sliver_gap = mfem::infinity();
  for (int ring = 0; ring < 11; ring++)
  {
    // Concave ring order: crossing1, free1 .. free9, crossing2, metal1 .. metal5; the
    // inner free knots are free3 .. free7 = slots 3 .. 7; rings 9 and 10 are the caps.
    for (int slot = 3; slot <= 7; slot++)
    {
      if (ring >= 9)
      {
        cap_sliver_knots.push_back(16 * ring + slot);
      }
      double gap = 0.0;
      for (int d = 0; d < 3; d++)
      {
        const double delta = constructed.knots[16 * ring + slot][d] -
                             constructed.knots[16 * ring + slot + 1][d];
        gap += delta * delta;
      }
      sliver_gap = std::min(sliver_gap, std::sqrt(gap));
    }
  }
  // The cap rings' inner gap: 0.00392 R; the device element is 0.05 = 2.5 R.
  CHECK_THAT(sliver_gap / R, WithinRel(0.0039178, 1.0e-3));
  CHECK(sliver_gap < 0.05 / 500.0);

  std::vector<std::unique_ptr<Mesh>> meshes;
  meshes.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(*notch)));
  LaplaceOperator laplace(iodata, meshes);
  SurfaceResponseOperator response(iodata, laplace);
  // A potential linear in the height above the metal plane (the island's plane y = 0.5 is
  // the SA interface; the box's local z axis is its normal): zero on the metal, so the
  // conductor reference is zero and every free hat must read +-z_k exactly.
  mfem::ParGridFunction potential(&laplace.GetH1Space().Get());
  mfem::FunctionCoefficient potential_coefficient([](const mfem::Vector &x)
                                                  { return x[1] - 0.5; });
  potential.ProjectCoefficient(potential_coefficient);
  Vector potential_true;
  potential.GetTrueDofs(potential_true);
  const auto &names = response.GetModelNames();
  const auto traces = response.GetSpatialPatchTraces(potential_true);
  int constructed_patches = 0;
  for (const auto &trace : traces)
  {
    const auto &name = names.at(trace.model);
    if (name.find("@corner-angle48.75-cubic") == std::string::npos)
    {
      continue;
    }
    constructed_patches++;
    REQUIRE(trace.contour_size == 176);
    REQUIRE(trace.coefficients.size() >= 176);
    // The sign of the local z axis from the top cap ring (every knot at z = +R).
    const double top = trace.coefficients[16 * 9];
    REQUIRE(std::abs(top) > 0.5 * R);
    const double sign = top > 0.0 ? 1.0 : -1.0;
    for (int k = 0; k < 176; k++)
    {
      REQUIRE(std::isfinite(trace.coefficients[k]));
      if (constructed.zero[k])
      {
        continue;
      }
      // The PEC knots are pinned to the conductor reference (0) while the device field at
      // the metal-top ring z = t reads t: that inconsistency (the thin device vs the
      // coupon's metal thickness, the thin matrices' business) leaks into the free
      // neighbours on the metal rings by a fraction of t (0.22 t measured) and decays away
      // from them; every free hat reads its knot value within t, the hats of the five far
      // rings (|z| >= R / 3, the cap rings included) within 1e-4 R.
      const int ring = k / 16;
      const bool far = ring <= 1 || ring >= 8;
      CHECK_THAT(trace.coefficients[k],
                 WithinAbs(sign * constructed.knots[k][2], far ? 1.0e-4 * R : t));
    }
    for (const int k : cap_sliver_knots)
    {
      REQUIRE(!constructed.zero[k]);
      CHECK_THAT(trace.coefficients[k],
                 WithinAbs(sign * constructed.knots[k][2], 1.0e-4 * R));
    }
  }
  CHECK(constructed_patches == 1);
  const auto result = response.GetElectrostaticResponse(potential_true);
  CHECK(std::isfinite(result.domain_correction));
  CHECK(std::isfinite(result.fabricated_surface_energy.at(4)));
  CHECK(result.fabricated_surface_energy.at(4) > 0.0);
}

TEST_CASE("CornerRefinedRuleLibraryLoad", "[cornerbasisrefinement][Serial][Parallel]")
{
  // The load-time check of every coupon against the rule at its angle (ReadProcessLibrary,
  // fail closed): a family of the refined rule loads and its exact 90-degree corners match
  // the 90 node on a square island; a coupon whose basis points are the rule's at another
  // angle, a coupon carrying a connectivity angle, a coupon without the extra ring and a
  // coupon with a free knot in ZeroTraceIndices are
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
  for (const std::string name :
       {"family", "wrong-angle", "connectivity", "no-extra-ring", "spurious-zero"})
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
      // A free knot of the 120-degree coupon's z = 0 metal ring (ring 3 of the ascending
      // levels) added to ZeroTraceIndices: the crossings gate still passes (the crossings
      // are PEC knots and no FREE knot lies on the metal), the level set is the rule's, and
      // the load-time check refuses the knot's zero flag.
      const auto zero_slots = CornerZeroSlots(true, RefinedCornerTraceBasisRule());
      int free_slot = 0;
      while (std::find(zero_slots.begin(), zero_slots.end(), free_slot) != zero_slots.end())
      {
        free_slot++;
      }
      auto spurious = family;
      auto indices = spurious[2]["ZeroTraceIndices"].get<std::vector<int>>();
      indices.push_back(3 * 16 + free_slot + 1);
      std::sort(indices.begin(), indices.end());
      spurious[2]["ZeroTraceIndices"] = indices;
      Write("spurious-zero", spurious);
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
  CHECK_THROWS_WITH(Preflight("spurious-zero", *square_mesh),
                    ContainsSubstring("is not its trace basis rule's coupon") &&
                        ContainsSubstring("is in ZeroTraceIndices but is a free knot"));
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
