// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gridfunctionio.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <unordered_map>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"

namespace palace::fieldio
{

namespace
{

std::string Trim(const std::string &s)
{
  const auto b = s.find_first_not_of(" \t\r\n");
  const auto e = s.find_last_not_of(" \t\r\n");
  return (b == std::string::npos) ? std::string() : s.substr(b, e - b + 1);
}

std::string Unquote(const std::string &s)
{
  std::string t = Trim(s);
  if (t.size() >= 2 && t.front() == '"' && t.back() == '"')
  {
    t = t.substr(1, t.size() - 2);
  }
  return t;
}

FieldData ReadGmsh(std::istream &in, const std::string &filename, const std::string &name)
{
  std::unordered_map<long long, std::array<double, 3>> coords;
  std::unordered_map<long long, double> values;
  std::string line;
  while (std::getline(in, line))
  {
    line = Trim(line);
    if (line == "$Nodes")
    {
      std::getline(in, line);
      const long long n = std::stoll(Trim(line));
      for (long long i = 0; i < n; i++)
      {
        std::getline(in, line);
        std::istringstream ss(line);
        long long tag;
        double x, y, z;
        ss >> tag >> x >> y >> z;
        coords[tag] = {x, y, z};
      }
    }
    else if (line == "$NodeData")
    {
      // String tags (the field name is one of them).
      std::getline(in, line);
      const int num_string = std::stoi(Trim(line));
      std::vector<std::string> string_tags;
      for (int i = 0; i < num_string; i++)
      {
        std::getline(in, line);
        string_tags.push_back(Unquote(line));
      }
      // Real tags.
      std::getline(in, line);
      const int num_real = std::stoi(Trim(line));
      for (int i = 0; i < num_real; i++)
      {
        std::getline(in, line);
      }
      // Integer tags: first is the number of components, second the number of nodes.
      std::getline(in, line);
      const int num_int = std::stoi(Trim(line));
      long long ncomp = 1;
      if (num_int > 0)
      {
        std::getline(in, line);
        std::istringstream ss(line);
        ss >> ncomp;
      }
      const bool match =
          name.empty() || std::find(string_tags.begin(), string_tags.end(), name) !=
                              string_tags.end();
      // Read node entries until $EndNodeData.
      while (std::getline(in, line))
      {
        line = Trim(line);
        if (line == "$EndNodeData")
        {
          break;
        }
        if (line.empty())
        {
          continue;
        }
        std::istringstream ss(line);
        long long tag;
        ss >> tag;
        double value = 0.0;
        for (long long c = 0; c < ncomp; c++)
        {
          double v;
          ss >> v;
          if (c == 0)
          {
            value = v;
          }
        }
        if (match)
        {
          values[tag] = value;
        }
      }
    }
  }
  MFEM_VERIFY(!values.empty(), "Field \"" << name << "\" not found in \"" << filename
                                          << "\"!");
  FieldData data;
  data.points.reserve(values.size());
  data.values.reserve(values.size());
  for (const auto &[tag, value] : values)
  {
    auto it = coords.find(tag);
    MFEM_VERIFY(it != coords.end(),
                "Missing node " << tag << " in field file \"" << filename << "\"!");
    data.points.push_back(it->second);
    data.values.push_back(value);
  }
  return data;
}

// Extract the first (or named) ASCII <DataArray> in the XML substring [begin, end).
std::vector<double> ExtractVtuArray(const std::string &xml, std::size_t begin,
                                    std::size_t end, const std::string &name)
{
  std::size_t pos = begin;
  while (true)
  {
    pos = xml.find("<DataArray", pos);
    if (pos == std::string::npos || pos >= end)
    {
      break;
    }
    const std::size_t tag_end = xml.find('>', pos);
    const std::string tag = xml.substr(pos, tag_end - pos);
    bool match = name.empty();
    const auto npos = tag.find("Name=\"");
    if (npos != std::string::npos)
    {
      const auto nend = tag.find('"', npos + 6);
      if (name == tag.substr(npos + 6, nend - (npos + 6)))
      {
        match = true;
      }
    }
    if (match)
    {
      MFEM_VERIFY(tag.find("format=\"ascii\"") != std::string::npos ||
                      tag.find("format='ascii'") != std::string::npos,
                  "Only ASCII VTU field data is supported!");
      const std::size_t close = xml.find("</DataArray>", tag_end);
      MFEM_VERIFY(close != std::string::npos, "Malformed VTU DataArray!");
      std::istringstream ss(xml.substr(tag_end + 1, close - (tag_end + 1)));
      std::vector<double> vals;
      double v;
      while (ss >> v)
      {
        vals.push_back(v);
      }
      return vals;
    }
    pos = tag_end;
  }
  return {};
}

FieldData ReadVtu(const std::string &filename, const std::string &name)
{
  std::ifstream in(filename);
  MFEM_VERIFY(in.good(), "Unable to open field file \"" << filename << "\"!");
  const std::string xml((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());

  const std::size_t points_begin = xml.find("<Points");
  const std::size_t points_end = xml.find("</Points>");
  MFEM_VERIFY(points_begin != std::string::npos && points_end != std::string::npos,
              "Missing <Points> in VTU file \"" << filename << "\"!");
  const std::vector<double> coords =
      ExtractVtuArray(xml, points_begin, points_end, "");
  MFEM_VERIFY(coords.size() % 3 == 0, "Malformed VTU point coordinates!");

  const std::size_t data_begin = xml.find("<PointData");
  const std::size_t data_end = xml.find("</PointData>");
  MFEM_VERIFY(data_begin != std::string::npos && data_end != std::string::npos,
              "Missing <PointData> in VTU file \"" << filename << "\"!");
  const std::vector<double> values = ExtractVtuArray(xml, data_begin, data_end, name);
  MFEM_VERIFY(!values.empty(), "Field \"" << name << "\" not found in \"" << filename
                                          << "\"!");
  MFEM_VERIFY(values.size() == coords.size() / 3,
              "VTU field \"" << name << "\" size does not match the point count!");

  FieldData data;
  data.points.reserve(values.size());
  data.values = values;
  for (std::size_t i = 0; i < values.size(); i++)
  {
    data.points.push_back({coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]});
  }
  return data;
}

}  // namespace

FieldData ReadFieldFile(const std::string &filename, const std::string &name)
{
  const auto dot = filename.find_last_of('.');
  const std::string ext =
      (dot == std::string::npos) ? std::string() : filename.substr(dot);
  if (ext == ".msh")
  {
    std::ifstream in(filename);
    MFEM_VERIFY(in.good(), "Unable to open field file \"" << filename << "\"!");
    return ReadGmsh(in, filename, name);
  }
  if (ext == ".vtu" || ext == ".vtk")
  {
    return ReadVtu(filename, name);
  }
  MFEM_ABORT("Unsupported field file format \"" << filename
                                                << "\" (expected .msh or .vtu)!");
  return {};
}

FieldGridFunction BuildGridFunction(const Mesh &mesh, const FieldData &data)
{
  FieldGridFunction out;
  out.fec = std::make_unique<mfem::H1_FECollection>(1, mesh.Dimension());
  // The mesh is only used to define the new space, not modified.
  auto &pmesh = const_cast<mfem::ParMesh &>(mesh.Get());
  out.fes = std::make_unique<mfem::ParFiniteElementSpace>(&pmesh, out.fec.get());
  out.gf = std::make_unique<mfem::ParGridFunction>(out.fes.get());

  std::map<std::array<double, 3>, double> lookup;
  for (std::size_t i = 0; i < data.points.size(); i++)
  {
    lookup[data.points[i]] = data.values[i];
  }

  const int sdim = mesh.SpaceDimension();
  mfem::Array<int> dofs;
  for (int v = 0; v < mesh.Get().GetNV(); v++)
  {
    const double *x = mesh.Get().GetVertex(v);
    const std::array<double, 3> key = {x[0], (sdim > 1) ? x[1] : 0.0,
                                       (sdim > 2) ? x[2] : 0.0};
    auto it = lookup.find(key);
    MFEM_VERIFY(it != lookup.end(),
                "Field file does not contain solve mesh vertex " << v << " at (" << key[0]
                                                                 << ", " << key[1] << ", "
                                                                 << key[2] << ")!");
    out.fes->GetVertexDofs(v, dofs);
    MFEM_VERIFY(dofs.Size() == 1,
                "Unexpected vertex dof count in field grid function!");
    (*out.gf)(dofs[0]) = it->second;
  }
  out.gf->ExchangeFaceNbrData();
  return out;
}

}  // namespace palace::fieldio
