// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Measure actual first surface-element altitude next to physical process edges.
// Matching-boundary edges (attribute 1) and coplanar same-family slot seams are excluded.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using P = std::array<double, 3>;
P sub(P a, P b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
double dot(P a, P b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
P cross(P a, P b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
json quantiles(std::vector<double> a)
{
  json r = json::object();
  if (a.empty())
    return r;
  std::sort(a.begin(), a.end());
  for (double p : {0., 10., 50., 90., 99., 100.})
  {
    double x = (a.size() - 1) * p / 100.;
    auto i = static_cast<std::size_t>(std::floor(x));
    auto j = static_cast<std::size_t>(std::ceil(x));
    r[std::to_string(p)] = a[i] + (x - i) * (a[j] - a[i]);
  }
  return r;
}
struct Record
{
  int attribute;
  P normal;
  double altitude;
};
int main(int argc, char **argv)
{
  if (argc != 3)
    return 2;
  mfem::Mpi::Init(argc, argv);
  mfem::Mesh mesh(argv[1], 1, 1, true);
  std::map<std::pair<int, int>, std::vector<Record>> edges;
  auto point = [&](int i)
  {
    const auto *p = mesh.GetVertex(i);
    return P{p[0], p[1], p[2]};
  };
  mfem::Array<int> vertices;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    mesh.GetBdrElementVertices(be, vertices);
    int n = vertices.Size();
    if (n != 3 && n != 4)
      return 3;
    P a = point(vertices[0]), b = point(vertices[1]), c = point(vertices[2]);
    P normal = cross(sub(b, a), sub(c, a));
    double norm = std::sqrt(dot(normal, normal));
    if (norm == 0)
      return 4;
    for (double &x : normal)
      x /= norm;
    for (int i = 0; i < n; i++)
    {
      int va = vertices[i], vb = vertices[(i + 1) % n];
      P pa = point(va), pb = point(vb), axis = sub(pb, pa);
      double l2 = dot(axis, axis);
      if (l2 == 0)
        return 5;
      double height = 0;
      int count = 0;
      for (int j = 0; j < n; j++)
        if (vertices[j] != va && vertices[j] != vb)
        {
          P d = sub(point(vertices[j]), pa);
          P area_vector = cross(d, axis);
          height += std::sqrt(dot(area_vector, area_vector) / l2);
          count++;
        }
      edges[std::minmax(va, vb)].push_back(
          {mesh.GetBdrAttribute(be), normal, height / count});
    }
  }
  std::map<int, std::vector<double>> heights, lengths;
  std::map<int, double> below_length, total_length;
  long long features = 0;
  for (auto &[edge, records] : edges)
  {
    if (std::any_of(records.begin(), records.end(),
                    [](auto &r) { return r.attribute == 1; }))
      continue;
    bool physical = false;
    for (auto &a : records)
      for (auto &b : records)
        if (a.attribute / 1000 != b.attribute / 1000 ||
            std::abs(dot(a.normal, b.normal)) < std::cos(30. * std::acos(-1.) / 180))
          physical = true;
    if (!physical)
      continue;
    features++;
    double length = std::sqrt(dot(sub(point(edge.second), point(edge.first)),
                                  sub(point(edge.second), point(edge.first))));
    for (auto &r : records)
    {
      int family = r.attribute / 1000;
      heights[family].push_back(r.altitude);
      lengths[family].push_back(length);
      total_length[family] += length;
      if (r.altitude <= .00200001)
        below_length[family] += length;
    }
  }
  json report = {{"Mesh", argv[1]},
                 {"PhysicalMeshEdges", features},
                 {"Units", "input mesh units, um for coupons"},
                 {"Definition",
                  "Adjacent surface triangle altitude, or mean opposing-vertex distance "
                  "for quads. Corner-node geometry; not h/p."},
                 {"Families", json::object()}};
  for (auto &[family, h] : heights)
    report["Families"][std::to_string(family)] = {
        {"Samples", h.size()},
        {"AltitudePercentiles", quantiles(h)},
        {"TangentEdgeLengthPercentiles", quantiles(lengths[family])},
        {"IncidentEdgeLength", total_length[family]},
        {"IncidentEdgeLengthFractionWithAltitudeAtMost2nm",
         below_length[family] / total_length[family]}};
  std::ofstream out(argv[2]);
  if (!out)
    return 6;
  out << report.dump(2) << '\n';
  return out ? 0 : 6;
}
