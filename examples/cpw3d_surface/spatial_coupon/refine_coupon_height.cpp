// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Experimental spatially graded Z refinement of a coarse-height prism seed.
// Original units must be um; fabrication process is 100 nm metal / 50 nm etch.
// Uses native MFEM nonconforming refinement, never unmatched Gmsh slab interfaces.
#include <algorithm>
#include <array>
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

struct BoundaryPoint
{
  double x, y, z;
  bool physical;
  int index;
};
using Segment = std::array<double, 4>;

int main(int argc, char **argv)
{
  if (argc != 6)
  {
    std::cerr << "seed.msh boundary.csv fabricated(0|1) output.mesh max_elements\n";
    return 2;
  }
  mfem::Mpi::Init(argc, argv);
  try
  {
    const bool fabricated = std::stoi(argv[3]) != 0;
    const int max_elements = std::stoi(argv[5]);
    std::ifstream f(argv[2]);
    if (!f)
      throw std::runtime_error("Cannot open boundary file");
    std::map<int, std::vector<BoundaryPoint>> loops;
    std::string line;
    std::getline(f, line);
    std::set<double> planes;
    while (std::getline(f, line))
    {
      if (line.empty())
        continue;
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      int loop, index, conductor, hole;
      double z, x, y;
      std::string cls;
      if (!(row >> loop >> index >> conductor >> z >> hole >> cls >> x >> y))
        throw std::runtime_error("Invalid boundary row");
      if (cls != "Physical" && cls != "Continuation")
        throw std::runtime_error("Invalid boundary class");
      loops[loop].push_back({x, y, z, cls == "Physical", index});
      planes.insert(z);
      if (fabricated)
      {
        planes.insert(z - .05);
        planes.insert(z + .1);
      }
    }
    std::vector<Segment> edges;
    for (auto &[id, p] : loops)
    {
      std::sort(p.begin(), p.end(), [](auto &a, auto &b) { return a.index < b.index; });
      for (std::size_t i = 0; i < p.size(); i++)
        if (p[i].physical)
        {
          auto &q = p[(i + 1) % p.size()];
          edges.push_back({p[i].x, p[i].y, q.x, q.y});
        }
    }
    if (edges.empty())
      throw std::runtime_error("No physical edges");
    mfem::Mesh mesh(argv[1], 1, 1, true);
    mesh.EnsureNCMesh();
    nlohmann::json report = {{"Input", argv[1]},
                             {"MinimumRequestedHeight", 0.002},
                             {"MaximumRequestedHeight", 0.5},
                             {"DistanceGrowth", 0.4},
                             {"NCLimit", 2},
                             {"SeedElements", mesh.GetNE()},
                             {"Passes", nlohmann::json::array()}};
    mfem::Array<int> vertices;
    bool converged = false;
    for (int pass = 0; pass < 12; pass++)
    {
      mfem::Array<mfem::Refinement> marked;
      for (int el = 0; el < mesh.GetNE(); el++)
      {
        if (mesh.GetElementBaseGeometry(el) != mfem::Geometry::PRISM)
          throw std::runtime_error("Expected prism seed");
        mesh.GetElementVertices(el, vertices);
        double x = 0, y = 0, zmin = 1e300, zmax = -1e300;
        for (int v : vertices)
        {
          const double *p = mesh.GetVertex(v);
          x += p[0] / vertices.Size();
          y += p[1] / vertices.Size();
          zmin = std::min(zmin, p[2]);
          zmax = std::max(zmax, p[2]);
        }
        double radius = 0;
        for (int v : vertices)
        {
          auto *p = mesh.GetVertex(v);
          radius = std::max(radius, std::hypot(p[0] - x, p[1] - y));
        }
        double xy_distance = 1e300;
        for (auto &e : edges)
        {
          double dx = e[2] - e[0], dy = e[3] - e[1], den = dx * dx + dy * dy;
          double t =
              den > 0 ? std::clamp(((x - e[0]) * dx + (y - e[1]) * dy) / den, 0., 1.) : 0.;
          xy_distance =
              std::min(xy_distance, std::hypot(x - e[0] - t * dx, y - e[1] - t * dy));
        }
        // Conservative distance bound: do not miss an edge merely because the
        // element centroid is far from it. Z intervals containing a plane have zero
        // distance.
        xy_distance = std::max(0., xy_distance - radius);
        double z_distance = 1e300;
        for (double plane : planes)
          z_distance = std::min(
              z_distance, plane >= zmin && plane <= zmax
                              ? 0.
                              : std::min(std::abs(plane - zmin), std::abs(plane - zmax)));
        double desired = std::clamp(0.4 * std::hypot(xy_distance, z_distance), 0.002, 0.5);
        if (zmax - zmin > desired * (1. + 1e-10))
          marked.Append(mfem::Refinement(el, mfem::Refinement::Z));
      }
      report["Passes"].push_back(
          {{"Pass", pass}, {"Elements", mesh.GetNE()}, {"Marked", marked.Size()}});
      std::cout << "pass=" << pass << " elements=" << mesh.GetNE()
                << " marked=" << marked.Size() << std::endl;
      if (marked.Size() == 0)
      {
        converged = true;
        break;
      }
      if (static_cast<long long>(mesh.GetNE()) + marked.Size() > max_elements)
        throw std::runtime_error("Predicted mesh budget exceeded");
      mesh.GeneralRefinement(marked, 1, 2);
      if (mesh.GetNE() > max_elements)
        throw std::runtime_error("Mesh budget exceeded after NC closure");
    }
    if (!converged)
      throw std::runtime_error("Refinement did not reach target heights within pass limit");
    // Preserve the planar geometry in quadratic representation, including hanging-node
    // constraints.
    mesh.SetCurvature(2);
    report["Elements"] = mesh.GetNE();
    report["Nonconforming"] = mesh.Nonconforming();
    std::ofstream out(argv[4]);
    if (!out)
      throw std::runtime_error("Cannot write mesh");
    out.precision(17);
    mesh.Print(out);
    out.close();
    std::ofstream metadata(std::string(argv[4]) + ".study.json");
    metadata << report.dump(2) << '\n';
    std::cout << argv[4] << " elements=" << mesh.GetNE() << std::endl;
  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
