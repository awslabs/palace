// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_INTERFACEOWNERSHIP_HPP
#define PALACE_UTILS_INTERFACEOWNERSHIP_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <mfem.hpp>

namespace palace
{

// Opt-in surface-only partition. Signature segments have already-expanded
// endpoints. Boundary segments optionally select the conductor first (the SA
// rule); beyond 3*Radius the signature segments compete globally. One group
// belongs to one process plane. Curves must not be silently replaced by chords.
class InterfaceOwnershipPartition
{
  struct Segment
  {
    int conductor, slot;
    bool boundary;
    std::array<double, 3> first, direction;
    double length_squared;
  };
  std::vector<Segment> segments;
  std::set<int> slots;
  double radius = 0.0;
  double signature_tolerance = 0.0;
  bool hierarchical = false;

  static std::vector<std::string> Split(const std::string &line)
  {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
    {
      const auto first = value.find_first_not_of(" \t\r");
      const auto last = value.find_last_not_of(" \t\r");
      result.push_back(first == std::string::npos ? ""
                                                  : value.substr(first, last - first + 1));
    }
    if (!line.empty() && line.back() == ',')
    {
      result.emplace_back();
    }
    return result;
  }
  static double Number(const std::string &value)
  {
    std::size_t count = 0;
    const double result = std::stod(value, &count);
    MFEM_VERIFY(count == value.size() && std::isfinite(result),
                "Invalid interface ownership number!");
    return result;
  }
  static int Integer(const std::string &value)
  {
    std::size_t count = 0;
    const long result = std::stol(value, &count);
    MFEM_VERIFY(count == value.size() && result >= 0 &&
                    result <= std::numeric_limits<int>::max(),
                "Invalid interface ownership integer!");
    return static_cast<int>(result);
  }
  static double Distance(const mfem::Vector &point, const Segment &segment)
  {
    const int dimension = segment.boundary ? 2 : 3;
    double projection = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      projection += (point(d) - segment.first[d]) * segment.direction[d];
    }
    projection = std::clamp(projection / segment.length_squared, 0.0, 1.0);
    double squared = 0.0;
    for (int d = 0; d < dimension; d++)
    {
      const double delta = point(d) - segment.first[d] - projection * segment.direction[d];
      squared += delta * delta;
    }
    return std::sqrt(squared);
  }

public:
  InterfaceOwnershipPartition(const std::string &filename, int group,
                              double mesh_coordinate_scale)
  {
    MFEM_VERIFY(std::isfinite(mesh_coordinate_scale) && mesh_coordinate_scale > 0.0,
                "Invalid interface ownership coordinate scale!");
    std::ifstream input(filename);
    MFEM_VERIFY(input, "Unable to open interface ownership file \"" << filename << "\"!");
    const std::vector<std::string> header = {
        "Group", "Conductor", "Slot", "Role", "Radius", "X0", "Y0", "Z0", "X1", "Y1", "Z1"};
    std::string line;
    MFEM_VERIFY(std::getline(input, line) && Split(line) == header,
                "Invalid interface ownership CSV header!");
    double plane = 0.0, radius_input = 0.0;
    std::set<int> boundary_conductors, signature_conductors;
    while (std::getline(input, line))
    {
      if (line.find_first_not_of(" \t\r") == std::string::npos)
      {
        continue;
      }
      const auto row = Split(line);
      MFEM_VERIFY(row.size() == header.size(), "Invalid interface ownership CSV row!");
      const int row_group = Integer(row[0]);
      Segment segment;
      segment.conductor = Integer(row[1]);
      segment.slot = Integer(row[2]);
      MFEM_VERIFY(segment.conductor > 0 && (row[3] == "Boundary" || row[3] == "Signature"),
                  "Invalid interface ownership conductor or role!");
      segment.boundary = row[3] == "Boundary";
      const double r = Number(row[4]);
      MFEM_VERIFY(r > 0.0, "Interface ownership radius must be positive!");
      for (int d = 0; d < 3; d++)
      {
        segment.first[d] = Number(row[5 + d]) / mesh_coordinate_scale;
        segment.direction[d] =
            Number(row[8 + d]) / mesh_coordinate_scale - segment.first[d];
      }
      MFEM_VERIFY(std::abs(segment.direction[2]) <= 1e-12 * r / mesh_coordinate_scale,
                  "Interface ownership segments must lie in a process plane!");
      segment.direction[2] = 0.0;
      segment.length_squared = segment.direction[0] * segment.direction[0] +
                               segment.direction[1] * segment.direction[1];
      MFEM_VERIFY(std::isfinite(segment.length_squared) && segment.length_squared > 0.0,
                  "Invalid or zero-length interface ownership segment!");
      if (row_group != group)
      {
        continue;
      }
      if (segments.empty())
      {
        plane = segment.first[2];
        radius_input = r;
      }
      MFEM_VERIFY(r == radius_input && std::abs(segment.first[2] - plane) <=
                                           1e-9 * r / mesh_coordinate_scale,
                  "One interface ownership group must have a common radius and plane!");
      if (segment.boundary)
      {
        hierarchical = true;
        boundary_conductors.insert(segment.conductor);
      }
      else
      {
        slots.insert(segment.slot);
        signature_conductors.insert(segment.conductor);
      }
      segments.push_back(segment);
    }
    MFEM_VERIFY(!input.bad() && !slots.empty(),
                "Empty or unreadable interface ownership group!");
    MFEM_VERIFY(!hierarchical || boundary_conductors == signature_conductors,
                "Interface ownership boundary/signature conductor sets differ!");
    radius = radius_input / mesh_coordinate_scale;
    signature_tolerance = 1e-12 * std::max(radius_input, 1.0) / mesh_coordinate_scale;
    MFEM_VERIFY(std::isfinite(radius) && radius > 0.0 && std::isfinite(signature_tolerance),
                "Interface ownership length scaling overflowed!");
  }

  const std::set<int> &Slots() const { return slots; }

  int SelectSlot(const mfem::Vector &point) const
  {
    MFEM_VERIFY(point.Size() == 3 && std::all_of(point.begin(), point.end(),
                                                 [](double x) { return std::isfinite(x); }),
                "Interface ownership needs a finite three-dimensional point!");
    int conductor = -1;
    if (hierarchical)
    {
      double best = mfem::infinity();
      for (const auto &segment : segments)
      {
        if (segment.boundary)
        {
          best = std::min(best, Distance(point, segment));
        }
      }
      const double tolerance = 1e-9 * radius;
      if (best <= 3 * radius + tolerance)
      {
        conductor = std::numeric_limits<int>::max();
        for (const auto &segment : segments)
        {
          if (segment.boundary && Distance(point, segment) <= best + tolerance)
          {
            conductor = std::min(conductor, segment.conductor);
          }
        }
      }
    }
    double best = mfem::infinity();
    for (const auto &segment : segments)
    {
      if (!segment.boundary && (conductor < 0 || segment.conductor == conductor))
      {
        best = std::min(best, Distance(point, segment));
      }
    }
    int slot = std::numeric_limits<int>::max();
    for (const auto &segment : segments)
    {
      if (!segment.boundary && (conductor < 0 || segment.conductor == conductor) &&
          Distance(point, segment) <= best + signature_tolerance)
      {
        slot = std::min(slot, segment.slot);
      }
    }
    MFEM_VERIFY(slot != std::numeric_limits<int>::max(),
                "Unmatched interface ownership point!");
    return slot;
  }
};

}  // namespace palace

#endif
