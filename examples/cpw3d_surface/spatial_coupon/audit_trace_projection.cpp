// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Standalone matching-surface audit: no volume space, assembly, MPI launch or solve.
// Build with -I<scratch> for the header extracted by prepare_trace_projection_audit.py.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "production_trace_coefficient.hpp"

using Point = std::array<double, 3>;
using Sparse = std::vector<std::pair<int, double>>;
using Tri = std::array<Point, 3>;
Point Sub(const Point &a, const Point &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Point Add(const Point &a, const Point &b, double t)
{
  return {a[0] + t * b[0], a[1] + t * b[1], a[2] + t * b[2]};
}
double Dot(const Point &a, const Point &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Point Cross(const Point &a, const Point &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double Norm(const Point &a)
{
  return std::sqrt(Dot(a, a));
}
Point Physical(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  mfem::Vector x(3);
  T.Transform(ip, x);
  return {x[0], x[1], x[2]};
}

struct AffineMap
{
  Point a, b, c;
  double bb, bc, cc, det;
  AffineMap(const Point &a_, const Point &b_, const Point &c_)
    : a(a_), b(Sub(b_, a_)), c(Sub(c_, a_)), bb(Dot(b, b)), bc(Dot(b, c)), cc(Dot(c, c)),
      det(bb * cc - bc * bc)
  {
    MFEM_VERIFY(det > 0.0, "Degenerate affine map");
  }
  mfem::IntegrationPoint Inverse(const Point &x) const
  {
    auto d = Sub(x, a);
    mfem::IntegrationPoint p;
    p.Set2((cc * Dot(b, d) - bc * Dot(c, d)) / det,
           (bb * Dot(c, d) - bc * Dot(b, d)) / det);
    return p;
  }
  Point Forward(const mfem::IntegrationPoint &p) const
  {
    return Add(Add(a, b, p.x), c, p.y);
  }
};

struct TraceBank
{
  std::vector<TracePotentialCoefficient> sources;
  std::vector<Tri> triangles;
  std::vector<std::array<Point, 2>> bounds;
  // The complete input bank is validated externally, then checked again for geometry.
  explicit TraceBank(const char *manifest)
  {
    std::ifstream in(manifest);
    MFEM_VERIFY(in, "Missing source manifest");
    std::string path;
    while (std::getline(in, path))
    {
      sources.emplace_back(path, 3, 1.0, 1.0);
    }
    MFEM_VERIFY(!sources.empty(), "Empty bank");
    for (const auto &t : sources[0].triangles)
    {
      Tri tri{t[0].point, t[1].point, t[2].point};
      triangles.push_back(tri);
      Point lo = tri[0], hi = tri[0];
      for (const auto &p : tri)
      {
        for (int d = 0; d < 3; d++)
        {
          lo[d] = std::min(lo[d], p[d]);
          hi[d] = std::max(hi[d], p[d]);
        }
      }
      bounds.push_back({lo, hi});
    }
    for (const auto &s : sources)
    {
      MFEM_VERIFY(s.curve.empty() && s.triangles.size() == triangles.size(),
                  "Requires consistent surface source bank");
      for (size_t t = 0; t < triangles.size(); t++)
      {
        for (int j = 0; j < 3; j++)
        {
          MFEM_VERIFY(s.triangles[t][j].point == triangles[t][j], "Geometry mismatch");
        }
      }
    }
  }
  Sparse Values(int t, const std::array<double, 3> &w) const
  {
    Sparse result;
    for (size_t s = 0; s < sources.size(); s++)
    {
      const auto &v = sources[s].triangles[t];
      double value = w[0] * v[0].value + w[1] * v[1].value + w[2] * v[2].value;
      if (value != 0.0)
      {
        result.emplace_back(s, value);
      }
    }
    return result;
  }
  Sparse Eval(const Point &point, double &distance) const
  {
    double best = std::numeric_limits<double>::infinity();
    int chosen = -1;
    for (size_t t = 0; t < triangles.size(); t++)
    {
      double lower = 0.0;
      for (int d = 0; d < 3; d++)
      {
        double delta = std::max(
            {bounds[t][0][d] - point[d] - 1e-12, point[d] - bounds[t][1][d] - 1e-12, 0.0});
        lower += delta * delta;
      }
      if (lower > best)
      {
        continue;
      }
      auto result =
          TracePotentialCoefficient::EvaluateTriangle(sources[0].triangles[t], point);
      if (result.first < best)  // Actual production strict comparison / ascending ID tie.
      {
        best = result.first;
        chosen = t;
      }
    }
    MFEM_VERIFY(chosen >= 0, "No nearest source triangle");
    auto tri = sources[0].triangles[chosen];
    std::array<double, 3> w;
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        tri[k].value = (j == k) ? 1.0 : 0.0;
      }
      w[j] = TracePotentialCoefficient::EvaluateTriangle(tri, point).second;
    }
    distance = std::sqrt(best);
    return Values(chosen, w);
  }
};

std::vector<Point> Clip(std::vector<Point> polygon, const Tri &triangle)
{
  Point normal = Cross(Sub(triangle[1], triangle[0]), Sub(triangle[2], triangle[0]));
  for (int j = 0; j < 3 && !polygon.empty(); j++)
  {
    Point a = triangle[j], edge = Sub(triangle[(j + 1) % 3], a);
    auto signed_distance = [&](const Point &p)
    { return Dot(Cross(edge, Sub(p, a)), normal); };
    std::vector<Point> result;
    Point p = polygon.back();
    double dp = signed_distance(p);
    for (const auto &q : polygon)
    {
      double dq = signed_distance(q);
      if ((dp >= 0.0) != (dq >= 0.0))
      {
        result.push_back(Add(p, Sub(q, p), dp / (dp - dq)));
      }
      if (dq >= 0.0)
      {
        result.push_back(q);
      }
      p = q;
      dp = dq;
    }
    polygon = std::move(result);
  }
  return polygon;
}

struct Piece
{
  int source_triangle;
  Tri triangle;
  double twice_area;
};

int Run(const char *mesh_path, const char *manifest, const char *prefix)
{
  TraceBank bank(manifest);
  int ns = static_cast<int>(bank.sources.size());
  mfem::Mesh mesh(mesh_path, 1, 1, false);
  MFEM_VERIFY(mesh.Dimension() == 2 && mesh.SpaceDimension() == 3,
              "Requires matching-only 2D mesh embedded in 3D");
  mfem::H1_FECollection fec(4, 2, mfem::BasisType::GaussLobatto);
  mfem::FiniteElementSpace fes(&mesh, &fec);
  const int nd = fes.GetNDofs();
  mfem::DenseMatrix nodal(ns, nd);
  nodal = 0.0;
  std::vector<int> counts(nd, 0);
  std::vector<Point> nodal_coordinates(nd);
  std::vector<char> ground(nd, false);
  mfem::Array<int> dofs;
  for (int b = 0; b < mesh.GetNBE(); b++)
  {
    int attr = mesh.GetBdrAttribute(b);
    MFEM_VERIFY(attr == 5001 || attr == 6001, "Unverified surface boundary attribute");
    fes.GetBdrElementDofs(b, dofs);
    for (int d : dofs)
    {
      MFEM_VERIFY(d >= 0, "Unexpected signed H1 DOF");
      ground[d] = true;
    }
  }
  double max_distance = 0.0, max_spread = 0.0, max_reference_difference = 0.0;
  double affine_defect = 0.0;
  int reference_checks = 0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    MFEM_VERIFY(mesh.GetAttribute(e) == 1, "Nonmatching face in surface mesh");
    fes.GetElementDofs(e, dofs);
    auto &T = *mesh.GetElementTransformation(e);
    const auto &nodes = fes.GetFE(e)->GetNodes();
    auto geom = mesh.GetElementBaseGeometry(e);
    MFEM_VERIFY(geom == mfem::Geometry::TRIANGLE || geom == mfem::Geometry::SQUARE,
                "Unsupported matching face geometry");
    const auto &corners = *mfem::Geometries.GetVertices(geom);
    AffineMap map(Physical(T, corners.IntPoint(0)), Physical(T, corners.IntPoint(1)),
                  Physical(T, corners.IntPoint(geom == mfem::Geometry::SQUARE ? 3 : 2)));
    for (int j = 0; j < dofs.Size(); j++)
    {
      Point p = Physical(T, nodes.IntPoint(j));
      affine_defect = std::max(affine_defect, Norm(Sub(p, map.Forward(nodes.IntPoint(j)))));
      MFEM_VERIFY(affine_defect < 1e-10,
                  "Nonaffine geometry: split integration is not authorized");
      double distance;
      auto values = bank.Eval(p, distance);
      max_distance = std::max(max_distance, distance);
      int d = dofs[j];
      MFEM_VERIFY(d >= 0, "Unexpected signed H1 DOF");
      std::vector<double> full(ns, 0.0);
      for (auto [s, v] : values)
      {
        full[s] = v;
      }
      if (counts[d])
      {
        for (int s = 0; s < ns; s++)
        {
          max_spread = std::max(max_spread, std::abs(nodal(s, d) / counts[d] - full[s]));
        }
      }
      if (j == 0 && e % std::max(1, mesh.GetNE() / 100) == 0)
      {
        for (int s = 0; s < ns; s++)
        {
          double v = bank.sources[s].Eval(T, nodes.IntPoint(j));
          max_reference_difference =
              std::max(max_reference_difference, std::abs(v - full[s]));
          reference_checks++;
        }
      }
      for (auto [s, v] : values)
      {
        nodal(s, d) += v;
      }
      if (counts[d] == 0)
      {
        nodal_coordinates[d] = p;
      }
      counts[d]++;
    }
  }
  for (int d = 0; d < nd; d++)
  {
    MFEM_VERIFY(counts[d] > 0, "Unused matching DOF");
    for (int s = 0; s < ns; s++)
    {
      nodal(s, d) /= counts[d];
    }
  }
  std::cout << std::setprecision(17) << "Nodal diagnostics: distance=" << max_distance
            << " spread=" << max_spread
            << " production_difference=" << max_reference_difference << std::endl;
  // Approved diagnostic-only guard: Ericson closest-point arithmetic on thin source
  // triangles returns distances up to 2.1e-10 micrometers for covered prism nodes.
  MFEM_VERIFY(max_distance < 1e-8 && max_spread < 1e-9 && max_reference_difference < 1e-12,
              "Source coverage/shared-node/production coefficient discrepancy");
  std::vector<double> ground_max(ns, 0.0), nodal_max(ns, 0.0);
  for (int d = 0; d < nd; d++)
  {
    for (int s = 0; s < ns; s++)
    {
      if (ground[d])
      {
        ground_max[s] = std::max(ground_max[s], std::abs(nodal(s, d)));
      }
      else
      {
        nodal_max[s] = std::max(nodal_max[s], std::abs(nodal(s, d)));
      }
    }
  }
  // Contact-curve diagnostics use the raw production coefficient, not the grounded
  // polynomial. These finite-rule integrals need not resolve every source kink.
  std::array<std::vector<double>, 2> contact_norms{std::vector<double>(ns, 0.0),
                                                   std::vector<double>(ns, 0.0)};
  std::vector<double> contact_split_sq(ns, 0.0), contact_split_max(ns, 0.0);
  for (int b = 0; b < mesh.GetNBE(); b++)
  {
    auto &T = *mesh.GetBdrElementTransformation(b);
    mfem::IntegrationPoint begin, end;
    begin.Set2(0.0, 0.0);
    end.Set2(1.0, 0.0);
    Point a = Physical(T, begin), z = Physical(T, end);
    std::vector<double> cuts{0.0, 1.0};
    for (const auto &tri : bank.triangles)
    {
      AffineMap map(tri[0], tri[1], tri[2]);
      auto normal = Cross(map.b, map.c);
      if (std::abs(Dot(Sub(a, map.a), normal)) > 1e-10 * Norm(normal) ||
          std::abs(Dot(Sub(z, map.a), normal)) > 1e-10 * Norm(normal))
      {
        continue;
      }
      auto u = map.Inverse(a), v = map.Inverse(z);
      std::array<double, 3> first{1 - u.x - u.y, u.x, u.y};
      std::array<double, 3> last{1 - v.x - v.y, v.x, v.y};
      double lo = 0.0, hi = 1.0;
      for (int j = 0; j < 3; j++)
      {
        double delta = last[j] - first[j];
        if (std::abs(delta) < 1e-14)
        {
          if (first[j] < -1e-10)
          {
            lo = 2.0;
          }
        }
        else if (delta > 0.0)
        {
          lo = std::max(lo, -first[j] / delta);
        }
        else
        {
          hi = std::min(hi, -first[j] / delta);
        }
      }
      if (hi >= lo)
      {
        cuts.push_back(std::clamp(lo, 0.0, 1.0));
        cuts.push_back(std::clamp(hi, 0.0, 1.0));
      }
    }
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    auto at = [&](double t)
    {
      double distance;
      std::vector<double> values(ns, 0.0);
      for (auto [s, v] : bank.Eval(Add(a, Sub(z, a), t), distance))
      {
        values[s] = v;
        contact_split_max[s] = std::max(contact_split_max[s], std::abs(v));
      }
      MFEM_VERIFY(distance < 1e-8, "Off-surface split contact sample");
      return values;
    };
    auto previous = at(cuts[0]);
    for (size_t k = 1; k < cuts.size(); k++)
    {
      auto next = at(cuts[k]);
      auto middle = at((cuts[k - 1] + cuts[k]) / 2);
      double length = (cuts[k] - cuts[k - 1]) * Norm(Sub(z, a));
      for (int s = 0; s < ns; s++)
      {
        MFEM_VERIFY(std::abs(middle[s] - (previous[s] + next[s]) / 2) < 1e-8,
                    "Contact split missed a P1 kink");
        contact_split_sq[s] +=
            length *
            (previous[s] * previous[s] + previous[s] * next[s] + next[s] * next[s]) / 3;
      }
      previous = std::move(next);
    }
    for (int level = 0; level < 2; level++)
    {
      const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, level ? 20 : 12);
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &ip = rule.IntPoint(q);
        T.SetIntPoint(&ip);
        double distance;
        for (auto [s, v] : bank.Eval(Physical(T, ip), distance))
        {
          contact_norms[level][s] += ip.weight * T.Weight() * v * v;
        }
        MFEM_VERIFY(distance < 1e-8, "Off-surface ground-contact sample");
      }
    }
  }
  std::ofstream stats(std::string(prefix) + "-nodal-stats.csv");
  stats << std::setprecision(17)
        << "source,raw_min,raw_max,raw_euclidean_norm,raw_nonzero_dofs,grounded_min,"
           "grounded_max,grounded_euclidean_norm,grounded_nonzero_dofs,raw_contact_l2_q12,"
           "raw_contact_l2_q20,raw_contact_l2_split,raw_contact_split_max\n";
  for (int s = 0; s < ns; s++)
  {
    double raw_min = std::numeric_limits<double>::infinity(), raw_max = -raw_min;
    double post_min = raw_min, post_max = raw_max, raw_sq = 0.0, post_sq = 0.0;
    int raw_count = 0, post_count = 0;
    for (int d = 0; d < nd; d++)
    {
      double raw = nodal(s, d), post = ground[d] ? 0.0 : raw;
      raw_min = std::min(raw_min, raw);
      raw_max = std::max(raw_max, raw);
      post_min = std::min(post_min, post);
      post_max = std::max(post_max, post);
      raw_sq += raw * raw;
      post_sq += post * post;
      raw_count += raw != 0.0;
      post_count += post != 0.0;
    }
    stats << s + 1 << ',' << raw_min << ',' << raw_max << ',' << std::sqrt(raw_sq) << ','
          << raw_count << ',' << post_min << ',' << post_max << ',' << std::sqrt(post_sq)
          << ',' << post_count << ',' << std::sqrt(contact_norms[0][s]) << ','
          << std::sqrt(contact_norms[1][s]) << ',' << std::sqrt(contact_split_sq[s]) << ','
          << contact_split_max[s] << '\n';
  }
  // Positive nodal Euclidean Gram matrices: numerical column-rank diagnostics only,
  // not physical L2 Gram matrices or a rigorous rank certificate.
  mfem::DenseMatrix raw_gram(ns), grounded_gram(ns);
  raw_gram = 0.0;
  grounded_gram = 0.0;
  for (int d = 0; d < nd; d++)
  {
    Sparse nonzero;
    for (int s = 0; s < ns; s++)
    {
      if (nodal(s, d) != 0.0)
      {
        nonzero.emplace_back(s, nodal(s, d));
      }
    }
    for (auto [s, a] : nonzero)
    {
      for (auto [t, b] : nonzero)
      {
        raw_gram(s, t) += a * b;
        if (!ground[d])
        {
          grounded_gram(s, t) += a * b;
        }
      }
    }
  }
  for (int post = 0; post < 2; post++)
  {
    std::ofstream gram(std::string(prefix) +
                       (post ? "-grounded-gram.csv" : "-raw-gram.csv"));
    gram << std::setprecision(17);
    const auto &matrix = post ? grounded_gram : raw_gram;
    for (int s = 0; s < ns; s++)
    {
      for (int t = 0; t < ns; t++)
      {
        gram << (t ? "," : "") << matrix(s, t);
      }
      gram << '\n';
    }
  }
  // Preserve every nonzero pre-grounding DOF of the two observed collapsed modes.
  // Their indices are the original complete 135-bank indices, never renumbered.
  if (ns >= 58)
  {
    std::ofstream witness(std::string(prefix) + "-source57-58-dofs.csv");
    witness << std::setprecision(17)
            << "dof_0based,x,y,z,physical_ground,raw57,grounded57,raw58,grounded58\n";
    for (int d = 0; d < nd; d++)
    {
      if (nodal(56, d) != 0.0 || nodal(57, d) != 0.0)
      {
        const auto &x = nodal_coordinates[d];
        witness << d << ',' << x[0] << ',' << x[1] << ',' << x[2] << ','
                << static_cast<int>(ground[d]) << ',' << nodal(56, d) << ','
                << (ground[d] ? 0.0 : nodal(56, d)) << ',' << nodal(57, d) << ','
                << (ground[d] ? 0.0 : nodal(57, d)) << '\n';
      }
    }
  }
  constexpr std::array<int, 3> orders{12, 16, 20};
  struct Norms
  {
    double target = 0, raw_error = 0, error = 0, projected = 0;
  };
  std::vector<std::vector<Norms>> norms(orders.size(), std::vector<Norms>(ns));
  std::vector<int> crossing_count(ns, 0);
  std::vector<bool> example_written(ns, false);
  std::ofstream examples(std::string(prefix) + "-crossings.csv");
  examples << std::setprecision(17)
           << "source,element_0based,source_triangle_a,source_triangle_b,intersection_area_"
              "a,intersection_area_b,gradient_jump,element_centroid_x,element_centroid_y,"
              "element_centroid_z\n";
  double total_area = 0.0, max_area_defect = 0.0, worst_relative_area_defect = 0.0;
  long piece_count = 0;
  int crossing_faces = 0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    fes.GetElementDofs(e, dofs);
    auto &T = *mesh.GetElementTransformation(e);
    const auto *fe = fes.GetFE(e);
    auto geom = mesh.GetElementBaseGeometry(e);
    const auto &corners = *mfem::Geometries.GetVertices(geom);
    std::vector<Point> poly;
    for (int j = 0; j < corners.GetNPoints(); j++)
    {
      poly.push_back(Physical(T, corners.IntPoint(j)));
    }
    AffineMap map(poly[0], poly[1], poly[geom == mfem::Geometry::SQUARE ? 3 : 2]);
    double face_area =
        Norm(Cross(map.b, map.c)) * (geom == mfem::Geometry::TRIANGLE ? 0.5 : 1.0);
    Point lo = poly[0], hi = poly[0], centroid{0, 0, 0};
    for (const auto &p : poly)
    {
      centroid = Add(centroid, p, 1.0 / poly.size());
      for (int d = 0; d < 3; d++)
      {
        lo[d] = std::min(lo[d], p[d]);
        hi[d] = std::max(hi[d], p[d]);
      }
    }
    int axis = 0;
    for (int d = 1; d < 3; d++)
    {
      if (hi[d] - lo[d] < hi[axis] - lo[axis])
      {
        axis = d;
      }
    }
    MFEM_VERIFY(hi[axis] - lo[axis] < 1e-10, "Nonplanar matching face");
    std::vector<Piece> pieces;
    std::map<int, double> intersection_areas;
    for (size_t t = 0; t < bank.triangles.size(); t++)
    {
      const auto &tri = bank.triangles[t];
      if (std::abs(tri[0][axis] - lo[axis]) > 1e-10 ||
          std::abs(tri[1][axis] - lo[axis]) > 1e-10 ||
          std::abs(tri[2][axis] - lo[axis]) > 1e-10)
      {
        continue;
      }
      bool outside = false;
      for (int d = 0; d < 3; d++)
      {
        if (lo[d] > bank.bounds[t][1][d] + 1e-12 || hi[d] < bank.bounds[t][0][d] - 1e-12)
        {
          outside = true;
        }
      }
      if (outside)
      {
        continue;
      }
      auto clipped = Clip(poly, tri);
      for (size_t j = 1; j + 1 < clipped.size(); j++)
      {
        Tri part{clipped[0], clipped[j], clipped[j + 1]};
        double a2 = Norm(Cross(Sub(part[1], part[0]), Sub(part[2], part[0])));
        if (a2 > 1e-25)
        {
          pieces.push_back({static_cast<int>(t), part, a2});
          intersection_areas[t] += a2 / 2;
        }
      }
    }
    double covered = 0.0;
    for (const auto &[t, a] : intersection_areas)
    {
      covered += a;
    }
    max_area_defect = std::max(max_area_defect, std::abs(covered - face_area));
    worst_relative_area_defect =
        std::max(worst_relative_area_defect, std::abs(covered - face_area) / face_area);
    MFEM_VERIFY(std::abs(covered - face_area) < 1e-9 * face_area + 1e-13,
                "Source/FE polygon intersection coverage failure");
    total_area += covered;
    piece_count += pieces.size();
    std::set<int> active;
    for (int d : dofs)
    {
      for (int s = 0; s < ns; s++)
      {
        if (nodal(s, d) != 0.0)
        {
          active.insert(s);
        }
      }
    }
    std::vector<int> interior_triangles;
    for (auto [t, a] : intersection_areas)
    {
      for (int s = 0; s < ns; s++)
      {
        for (const auto &v : bank.sources[s].triangles[t])
        {
          if (v.value != 0.0)
          {
            active.insert(s);
          }
        }
      }
      if (a > std::max(1e-14, face_area * 1e-10))
      {
        interior_triangles.push_back(t);
      }
    }
    if (interior_triangles.size() > 1)
    {
      crossing_faces++;
      int a = interior_triangles[0];
      const auto &ta = bank.triangles[a];
      AffineMap ma(ta[0], ta[1], ta[2]);
      auto gradient =
          [](const AffineMap &m, const std::array<TracePotentialCoefficient::Sample, 3> &v)
      {
        double db = v[1].value - v[0].value, dc = v[2].value - v[0].value;
        Point zero{0, 0, 0};
        return Add(Add(zero, m.b, (m.cc * db - m.bc * dc) / m.det), m.c,
                   (m.bb * dc - m.bc * db) / m.det);
      };
      for (int s : active)
      {
        Point ga = gradient(ma, bank.sources[s].triangles[a]);
        for (size_t j = 1; j < interior_triangles.size(); j++)
        {
          int b = interior_triangles[j];
          const auto &tb = bank.triangles[b];
          AffineMap mb(tb[0], tb[1], tb[2]);
          double jump = Norm(Sub(ga, gradient(mb, bank.sources[s].triangles[b])));
          if (jump > 1e-9)
          {
            crossing_count[s]++;
            if (!example_written[s])
            {
              examples << s + 1 << ',' << e << ',' << a + 1 << ',' << b + 1 << ','
                       << intersection_areas[a] << ',' << intersection_areas[b] << ','
                       << jump << ',' << centroid[0] << ',' << centroid[1] << ','
                       << centroid[2] << '\n';
              example_written[s] = true;
            }
            break;
          }
        }
      }
    }
    mfem::Vector shape(dofs.Size());
    for (const auto &piece : pieces)
    {
      const auto &st = bank.triangles[piece.source_triangle];
      AffineMap source_map(st[0], st[1], st[2]);
      Point ab = Sub(piece.triangle[1], piece.triangle[0]),
            ac = Sub(piece.triangle[2], piece.triangle[0]);
      for (size_t level = 0; level < orders.size(); level++)
      {
        const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, orders[level]);
        for (int q = 0; q < rule.GetNPoints(); q++)
        {
          const auto &ip = rule.IntPoint(q);
          Point x = Add(Add(piece.triangle[0], ab, ip.x), ac, ip.y);
          auto ref = map.Inverse(x), src = source_map.Inverse(x);
          fe->CalcShape(ref, shape);
          double weight = ip.weight * piece.twice_area;
          for (int s : active)
          {
            const auto &v = bank.sources[s].triangles[piece.source_triangle];
            double target =
                (1 - src.x - src.y) * v[0].value + src.x * v[1].value + src.y * v[2].value;
            double raw = 0.0, projected = 0.0;
            for (int j = 0; j < dofs.Size(); j++)
            {
              double value = shape[j] * nodal(s, dofs[j]);
              raw += value;
              if (!ground[dofs[j]])
              {
                projected += value;
              }
            }
            auto &n = norms[level][s];
            n.target += weight * target * target;
            n.raw_error += weight * (raw - target) * (raw - target);
            n.error += weight * (projected - target) * (projected - target);
            n.projected += weight * projected * projected;
          }
        }
      }
    }
  }
  std::ofstream csv(std::string(prefix) + "-norms.csv");
  csv << std::setprecision(17)
      << "source,quadrature_order,target_l2,raw_interpolation_error_l2,grounded_error_l2,"
         "grounded_relative_error,projected_l2,maximum_ground_override,maximum_retained_"
         "nodal_value,kink_crossed_faces\n";
  for (int s = 0; s < ns; s++)
  {
    for (size_t level = 0; level < orders.size(); level++)
    {
      auto n = norms[level][s];
      csv << s + 1 << ',' << orders[level] << ',' << std::sqrt(n.target) << ','
          << std::sqrt(n.raw_error) << ',' << std::sqrt(n.error) << ','
          << (n.target > 0 ? std::sqrt(n.error / n.target) : 0.0) << ','
          << std::sqrt(n.projected) << ',' << ground_max[s] << ',' << nodal_max[s] << ','
          << crossing_count[s] << '\n';
    }
  }
  std::ofstream summary(std::string(prefix) + "-summary.json");
  summary << std::setprecision(17) << "{\n\"Elements\":" << mesh.GetNE()
          << ",\n\"H1TraceDOFs\":" << nd
          << ",\n\"GroundDOFs\":" << std::count(ground.begin(), ground.end(), true)
          << ",\n\"SourceCount\":" << ns << ",\n\"Area\":" << total_area
          << ",\n\"SplitTriangles\":" << piece_count
          << ",\n\"SourceTopologyCrossedFaces\":" << crossing_faces
          << ",\n\"MaximumAffineGeometryDefect\":" << affine_defect
          << ",\n\"MaximumSourceDistanceAtNodes\":" << max_distance
          << ",\n\"MaximumSharedNodeValueSpread\":" << max_spread
          << ",\n\"ProductionReferenceChecks\":" << reference_checks
          << ",\n\"MaximumProductionReferenceDifference\":" << max_reference_difference
          << ",\n\"MaximumFaceCoverageAreaDefect\":" << max_area_defect
          << ",\n\"MaximumRelativeFaceCoverageAreaDefect\":" << worst_relative_area_defect
          << "\n}\n";
  std::cout << "Audited " << ns << " sources, " << mesh.GetNE() << " matching faces, " << nd
            << " H1 trace DOFs; output " << prefix << std::endl;
  return 0;
}

int main(int argc, char **argv)
{
  try
  {
    MFEM_VERIFY(argc == 4, "usage: matching.msh source-files.txt output-prefix");
    return Run(argv[1], argv[2], argv[3]);
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
