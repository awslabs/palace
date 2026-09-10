// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Read-only mesh audit. No solution space or PDE operator is assembled.
// Usage: audit_coupon_mesh mesh.msh plan-view-boundary.csv fabricated(0|1) output.json
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using Point = std::array<double, 3>;
struct Edge
{
  Point a{}, b{};
};
struct Vertex
{
  Point x{};
  bool physical;
  int index;
};

std::vector<Edge> ReadEdges(const char *path)
{
  std::ifstream f(path);
  if (!f)
  {
    throw std::runtime_error("Cannot read boundary CSV");
  }
  std::string line;
  std::getline(f, line);
  std::map<int, std::vector<Vertex>> loops;
  while (std::getline(f, line))
  {
    if (line.empty())
    {
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    int loop, index, conductor, hole;
    double plane, x, y;
    std::string cls;
    if (!(row >> loop >> index >> conductor >> plane >> hole >> cls >> x >> y))
      throw std::runtime_error("Malformed boundary CSV");
    loops[loop].push_back({{x, y, plane}, cls == "Physical", index});
    if (cls != "Physical" && cls != "Continuation")
      throw std::runtime_error("Unknown edge class");
  }
  std::vector<Edge> result;
  for (auto &[id, vertices] : loops)
  {
    std::sort(vertices.begin(), vertices.end(),
              [](const auto &a, const auto &b) { return a.index < b.index; });
    for (std::size_t i = 0; i < vertices.size(); i++)
      if (vertices[i].physical)
        result.push_back({vertices[i].x, vertices[(i + 1) % vertices.size()].x});
  }
  if (result.empty())
    throw std::runtime_error("No physical edges in boundary");
  return result;
}

double EdgeDistance(const Point &p, const Edge &e)
{
  double dd = 0, t = 0;
  for (int i = 0; i < 2; i++)
  {
    double d = e.b[i] - e.a[i];
    dd += d * d;
    t += (p[i] - e.a[i]) * d;
  }
  t = dd > 0 ? std::clamp(t / dd, 0.0, 1.0) : 0.0;
  double d2 = 0;
  for (int i = 0; i < 2; i++)
  {
    double d = p[i] - e.a[i] - t * (e.b[i] - e.a[i]);
    d2 += d * d;
  }
  return std::sqrt(d2);
}
json Distribution(std::vector<double> a)
{
  if (a.empty())
    return json::object();
  std::sort(a.begin(), a.end());
  json q;
  for (double pct : {0., 1., 10., 50., 90., 95., 99., 99.9, 100.})
  {
    double at = (a.size() - 1) * pct / 100.;
    std::size_t lo = std::floor(at), hi = std::ceil(at);
    q[std::to_string(pct)] = a[lo] + (at - lo) * (a[hi] - a[lo]);
  }
  return q;
}
struct Count
{
  long long n = 0;
  double volume = 0;
  void Add(double v)
  {
    n++;
    volume += v;
  }
};
int Band(double x, const std::vector<double> &bounds)
{
  return std::upper_bound(bounds.begin(), bounds.end(), x) - bounds.begin();
}

int main(int argc, char **argv)
{
  if (argc != 5 && argc != 6)
  {
    std::cerr << "Usage: mesh boundary.csv fabricated(0|1) output.json\n";
    return 2;
  }
  const auto start = std::chrono::steady_clock::now();
  mfem::Mpi::Init(argc, argv);
  try
  {
    auto edges = ReadEdges(argv[2]);
    bool fab = std::stoi(argv[3]) != 0;
    mfem::Mesh mesh(argv[1], 1, 1, true);
    if (mesh.Dimension() != 3 || mesh.GetNE() == 0)
      throw std::runtime_error("Expected 3D mesh");
    double loaded =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::set<double> planes;
    for (const auto &e : edges)
    {
      planes.insert(e.a[2]);
      if (fab)
      {
        planes.insert(e.a[2] - .05);
        planes.insert(e.a[2] + .1);
      }
    }
    // Thresholds and plane offsets are in original mesh coordinates (um for these coupons).
    const std::vector<double> zb = {.01, .05, .1, .2, .5, 1.0},
                              db = {.01, .05, .1, .2, .8, 2.0};
    const std::vector<double> kb = {10, 50, 100, 250, 500, 1000, 2000, 5000};
    std::vector<Count> kh(kb.size() + 1), zh(zb.size() + 1), dh(db.size() + 1);
    std::vector<Count> joint((zb.size() + 1) * (db.size() + 1));
    std::vector<double> kappas, hmin, hmax;
    kappas.reserve(mesh.GetNE());
    hmin.reserve(mesh.GetNE());
    hmax.reserve(mesh.GetNE());
    struct Plan
    {
      double angle, shortest, longest, altitude, kappa;
      Point center;
      long long copies = 0;
    };
    std::map<std::array<double, 6>, Plan> plan;
    std::vector<std::pair<double, json>> worst;
    double volume = 0;
    std::map<int, double> material_volumes, boundary_areas;
    std::map<std::string, long long> geometry_counts;
    Count high, high_far, high_near, bulk;
    double affine_error = 0, vertical_error = 0;
    long long negative_center_jacobians = 0;
    mfem::DenseMatrix J(3), PJ(2), ideal(2);
    mfem::Vector center(3);
    mfem::Array<int> vertices;
    for (int el = 0; el < mesh.GetNE(); el++)
    {
      const auto geometry = mesh.GetElementBaseGeometry(el);
      if (geometry != mfem::Geometry::PRISM && geometry != mfem::Geometry::TETRAHEDRON)
        throw std::runtime_error("Study expects prisms or tetrahedra");
      geometry_counts[mfem::Geometry::Name[geometry]]++;
      mesh.GetElementJacobian(el, J);
      double sm = J.CalcSingularvalue(2), sx = J.CalcSingularvalue(0);
      if (!(sm > 0) || !std::isfinite(sx / sm))
        throw std::runtime_error("Singular element");
      double k = sx / sm;
      kappas.push_back(k);
      hmin.push_back(sm);
      hmax.push_back(sx);
      auto *T = mesh.GetElementTransformation(el);
      T->SetIntPoint(&mfem::Geometries.GetCenter(geometry));
      double det = T->Jacobian().Det();
      if (det <= 0)
        negative_center_jacobians++;
      double v = std::abs(det) * (geometry == mfem::Geometry::PRISM ? .5 : 1.0 / 6.0);
      volume += v;
      material_volumes[mesh.GetAttribute(el)] += v;
      mesh.GetElementCenter(el, center);
      Point p = {center[0], center[1], center[2]};
      double zd = 1e300, ed = 1e300;
      for (double z : planes)
        zd = std::min(zd, std::abs(p[2] - z));
      for (const auto &e : edges)
        ed = std::min(ed, EdgeDistance(p, e));
      kh[Band(k, kb)].Add(v);
      zh[Band(zd, zb)].Add(v);
      dh[Band(ed, db)].Add(v);
      joint[Band(zd, zb) * (db.size() + 1) + Band(ed, db)].Add(v);
      if (k > 1000)
      {
        high.Add(v);
        if (zd > .2)
          high_far.Add(v);
        else
          high_near.Add(v);
      }
      if (zd > .2 && ed > .8)
        bulk.Add(v);
      mesh.GetElementVertices(el, vertices);
      if (geometry == mfem::Geometry::TETRAHEDRON)
      {
        mfem::DenseMatrix corner_j(3);
        const auto *a = mesh.GetVertex(vertices[0]);
        for (int j = 0; j < 3; j++)
          for (int d = 0; d < 3; d++)
            corner_j(d, j) = mesh.GetVertex(vertices[j + 1])[d] - a[d];
        double geometric_volume = std::abs(corner_j.Det()) / 6.0;
        affine_error = std::max(affine_error, std::abs(geometric_volume / v - 1.));
        if (worst.size() < 30 || k > worst.back().first)
        {
          worst.push_back({k,
                           {{"ElementZeroBased", el},
                            {"Attribute", mesh.GetAttribute(el)},
                            {"Kappa", k},
                            {"Center", p},
                            {"SigmaMin", sm},
                            {"SigmaMax", sx},
                            {"DistanceToProcessPlanes", zd},
                            {"PlanDistanceToPhysicalEdge", ed}}});
          std::sort(worst.begin(), worst.end(),
                    [](auto &a, auto &b) { return a.first > b.first; });
          if (worst.size() > 30)
            worst.pop_back();
        }
        continue;
      }
      Point a, b, c;
      std::array<Point, 3> xy;
      for (int j = 0; j < 3; j++)
      {
        const double *q = mesh.GetVertex(vertices[j]);
        const double *r = mesh.GetVertex(vertices[j + 3]);
        xy[j] = {q[0], q[1], 0};
        vertical_error = std::max(vertical_error, std::hypot(q[0] - r[0], q[1] - r[1]));
      }
      a = xy[0];
      b = xy[1];
      c = xy[2];
      std::sort(xy.begin(), xy.end());
      std::array<double, 6> key = {xy[0][0], xy[0][1], xy[1][0],
                                   xy[1][1], xy[2][0], xy[2][1]};
      auto it = plan.find(key);
      if (it == plan.end())
      {
        double ab = std::hypot(b[0] - a[0], b[1] - a[1]),
               ac = std::hypot(c[0] - a[0], c[1] - a[1]),
               bc = std::hypot(c[0] - b[0], c[1] - b[1]);
        double shortest = std::min({ab, ac, bc}), longest = std::max({ab, ac, bc});
        double twicearea =
            std::abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
        double angle = 180;
        for (auto sides :
             {std::array<double, 3>{ab, ac, bc}, std::array<double, 3>{ab, bc, ac},
              std::array<double, 3>{ac, bc, ab}})
          angle = std::min(angle,
                           std::acos(std::clamp((sides[0] * sides[0] + sides[1] * sides[1] -
                                                 sides[2] * sides[2]) /
                                                    (2 * sides[0] * sides[1]),
                                                -1., 1.)) *
                               180 / std::acos(-1.));
        PJ(0, 0) = b[0] - a[0];
        PJ(1, 0) = b[1] - a[1];
        PJ(0, 1) = c[0] - a[0];
        PJ(1, 1) = c[1] - a[1];
        mfem::Geometries.JacToPerfJac(mfem::Geometry::TRIANGLE, PJ, ideal);
        Plan r{angle,
               shortest,
               longest,
               twicearea / longest,
               ideal.CalcSingularvalue(0) / ideal.CalcSingularvalue(1),
               {(a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3, 0},
               0};
        it = plan.emplace(key, r).first;
      }
      it->second.copies++;
      double geometric_volume =
          .5 * std::abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) *
          std::abs(mesh.GetVertex(vertices[3])[2] - mesh.GetVertex(vertices[0])[2]);
      if (v > 0)
        affine_error = std::max(affine_error, std::abs(geometric_volume / v - 1.));
      if (worst.size() < 30 || k > worst.back().first)
      {
        json r = {{"ElementZeroBased", el},
                  {"Attribute", mesh.GetAttribute(el)},
                  {"Kappa", k},
                  {"Center", p},
                  {"SigmaMin", sm},
                  {"SigmaMax", sx},
                  {"PlanKappa", it->second.kappa},
                  {"PlanMinimumAngleDegrees", it->second.angle},
                  {"PlanShortestEdge", it->second.shortest},
                  {"PlanMinimumAltitude", it->second.altitude},
                  {"DistanceToProcessPlanes", zd},
                  {"PlanDistanceToPhysicalEdge", ed}};
        worst.push_back({k, r});
        std::sort(worst.begin(), worst.end(),
                  [](auto &a, auto &b) { return a.first > b.first; });
        if (worst.size() > 30)
          worst.pop_back();
      }
    }
    auto count = [&](Count c)
    {
      return json{{"Elements", c.n},
                  {"ElementFraction", double(c.n) / mesh.GetNE()},
                  {"Volume", c.volume},
                  {"VolumeFraction", c.volume / volume}};
    };
    auto hist = [&](const std::vector<Count> &a, const std::vector<double> &b)
    {
      json r = json::array();
      for (std::size_t i = 0; i < a.size(); i++)
      {
        auto x = count(a[i]);
        x["Lower"] = i ? json(b[i - 1]) : json(0);
        x["Upper"] = i < b.size() ? json(b[i]) : json(nullptr);
        r.push_back(x);
      }
      return r;
    };
    json result = {
        {"Version", 1},
        {"Mesh", argv[1]},
        {"Boundary", argv[2]},
        {"Fabricated", fab},
        {"CoordinateUnits", "input mesh units (um for study inputs)"},
        {"Elements", mesh.GetNE()},
        {"LoadSeconds", loaded},
        {"KappaDefinition", "MFEM GetElementJacobian at reference center, "
                            "perfect-reference normalization, sigma_max/sigma_min"},
        {"VolumeNote", "Center Jacobian times reference volume; validated against "
                       "straight-prism corner volume"},
        {"MaximumVolumeDiscrepancy", affine_error},
        {"MaximumXYDriftAlongSweep", vertical_error},
        {"NonpositiveCenterJacobians", negative_center_jacobians},
        {"KappaPercentiles", Distribution(kappas)},
        {"SigmaMinPercentiles", Distribution(hmin)},
        {"SigmaMaxPercentiles", Distribution(hmax)},
        {"KappaHistogram", hist(kh, kb)},
        {"ProcessPlaneDistanceHistogram", hist(zh, zb)},
        {"PhysicalEdgePlanDistanceHistogram", hist(dh, db)},
        {"KappaAbove1000", count(high)},
        {"KappaAbove1000AndPlaneDistanceAbove0p2", count(high_far)},
        {"BulkPlaneDistanceAbove0p2AndEdgeDistanceAbove0p8", count(bulk)}};
    result["WorstElements"] = json::array();
    for (auto &x : worst)
      result["WorstElements"].push_back(x.second);
    std::vector<double> angles, shortest, altitudes, pk;
    json wp = json::array();
    long long count_short = 0, count_angle = 0;
    std::vector<std::pair<double, json>> worst_plan;
    for (auto &[key, r] : plan)
    {
      angles.push_back(r.angle);
      shortest.push_back(r.shortest);
      altitudes.push_back(r.altitude);
      pk.push_back(r.kappa);
      if (r.shortest < .002)
        count_short++;
      if (r.angle < 1.)
        count_angle++;
      double ed = 1e300;
      for (const auto &e : edges)
        ed = std::min(ed, EdgeDistance(r.center, e));
      json rec = {{"CoordinatesXY", key},
                  {"Center", r.center},
                  {"Kappa", r.kappa},
                  {"MinimumAngleDegrees", r.angle},
                  {"MinimumAltitude", r.altitude},
                  {"ShortestEdge", r.shortest},
                  {"VolumeCopies", r.copies},
                  {"PlanDistanceToPhysicalEdge", ed}};
      if (worst_plan.size() < 30 || r.kappa > worst_plan.back().first)
      {
        worst_plan.push_back({r.kappa, rec});
        std::sort(worst_plan.begin(), worst_plan.end(),
                  [](auto &a, auto &b) { return a.first > b.first; });
        if (worst_plan.size() > 30)
          worst_plan.pop_back();
      }
    }
    result["Plan"] = {{"UniqueTriangles", plan.size()},
                      {"KappaPercentiles", Distribution(pk)},
                      {"MinimumAnglePercentilesDegrees", Distribution(angles)},
                      {"ShortestEdgePercentiles", Distribution(shortest)},
                      {"MinimumAltitudePercentiles", Distribution(altitudes)},
                      {"TrianglesWithEdgeBelow2nm", count_short},
                      {"TrianglesWithAngleBelow1Degree", count_angle}};
    result["Plan"]["WorstTriangles"] = json::array();
    for (auto &x : worst_plan)
      result["Plan"]["WorstTriangles"].push_back(x.second);
    result["DistanceJointHistogram"] = json::array();
    for (std::size_t z = 0; z < zh.size(); z++)
      for (std::size_t d = 0; d < dh.size(); d++)
      {
        auto r = count(joint[z * dh.size() + d]);
        r["PlaneDistanceBin"] = z;
        r["EdgeDistanceBin"] = d;
        result["DistanceJointHistogram"].push_back(r);
      }
    result["TotalSeconds"] =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      auto *T = mesh.GetBdrElementTransformation(be);
      const auto geometry = mesh.GetBdrElementGeometry(be);
      T->SetIntPoint(&mfem::Geometries.GetCenter(geometry));
      boundary_areas[mesh.GetBdrAttribute(be)] +=
          T->Weight() * (geometry == mfem::Geometry::TRIANGLE ? .5 : 1.0);
    }
    result["MaterialVolumes"] = json::object();
    result["BoundaryAreas"] = json::object();
    for (auto [attribute, value] : material_volumes)
      result["MaterialVolumes"][std::to_string(attribute)] = value;
    for (auto [attribute, value] : boundary_areas)
      result["BoundaryAreas"][std::to_string(attribute)] = value;
    result["GeometryCounts"] = geometry_counts;
    result["VolumeNote"] =
        "Center Jacobian volume; compare with corner volume to check affine geometry";
    if (argc == 6)
    {
      if (plan.empty())
        throw std::runtime_error("Cannot export a prism plan from a tetrahedral mesh");
      // Export the existing plan exactly; experiments can change connectivity without
      // rerunning CAD fragmentation or changing vertex coordinates.
      std::map<std::array<double, 2>, int> node_ids;
      for (const auto &[key, data] : plan)
      {
        for (int i = 0; i < 3; i++)
        {
          std::array<double, 2> xy{key[2 * i], key[2 * i + 1]};
          if (!node_ids.count(xy))
            node_ids.emplace(xy, node_ids.size() + 1);
        }
      }
      std::ofstream nodes(std::string(argv[5]) + "-nodes.csv");
      std::ofstream tris(std::string(argv[5]) + "-triangles.csv");
      if (!nodes || !tris)
        throw std::runtime_error("Cannot export plan");
      nodes.precision(17);
      nodes << "node,x,y\n";
      for (const auto &[xy, id] : node_ids)
        nodes << id << ',' << xy[0] << ',' << xy[1] << '\n';
      tris << "triangle,a,b,c\n";
      int index = 0;
      for (const auto &[key, data] : plan)
      {
        std::array<int, 3> ids;
        for (int i = 0; i < 3; i++)
          ids[i] = node_ids.at({key[2 * i], key[2 * i + 1]});
        const double cross =
            (key[2] - key[0]) * (key[5] - key[1]) - (key[3] - key[1]) * (key[4] - key[0]);
        if (cross < 0)
          std::swap(ids[1], ids[2]);
        tris << ++index << ',' << ids[0] << ',' << ids[1] << ',' << ids[2] << '\n';
      }
    }
    std::ofstream out(argv[4]);
    if (!out)
      throw std::runtime_error("Cannot write report");
    out << result.dump(2) << '\n';
    std::cout << argv[4] << " elements=" << mesh.GetNE()
              << " max kappa=" << *std::max_element(kappas.begin(), kappas.end())
              << " seconds=" << result["TotalSeconds"] << std::endl;
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  return 0;
}
