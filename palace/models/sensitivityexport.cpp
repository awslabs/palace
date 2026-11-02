// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "sensitivityexport.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <fmt/format.h>
#include "fem/mesh.hpp"
#include "utils/communication.hpp"
#include "utils/enum_string.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Exports of more nodes than this, for the whole mesh, are flagged.
constexpr long long int large_export = 1000000;

// Mark the owned true DOFs of the first component of the mesh nodal space touched by the
// local elements (or boundary elements) for which in_set returns true.
template <typename InSet, typename GetDofs>
Vector MarkNodes(const mfem::ParFiniteElementSpace &nodal_fespace, int num_elem,
                 InSet &&in_set, GetDofs &&get_dofs)
{
  Vector marker_l(nodal_fespace.GetVSize()), marker(nodal_fespace.GetTrueVSize());
  marker_l = 0.0;
  mfem::Array<int> dofs;
  for (int e = 0; e < num_elem; e++)
  {
    if (in_set(e))
    {
      get_dofs(e, dofs);
      for (auto d : dofs)
      {
        marker_l(nodal_fespace.DofToVDof(d, 0)) = 1.0;
      }
    }
  }
  nodal_fespace.GetProlongationMatrix()->MultTranspose(marker_l, marker);
  return marker;
}

// Gather variable-length data from all processes to the root, in rank order.
template <typename T>
std::vector<T> GatherToRoot(const std::vector<T> &data, MPI_Datatype type, MPI_Comm comm)
{
  const int size = static_cast<int>(data.size());
  std::vector<int> counts(Mpi::Size(comm)), displs(Mpi::Size(comm));
  MPI_Gather(&size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  std::vector<T> all;
  if (Mpi::Root(comm))
  {
    for (std::size_t i = 1; i < counts.size(); i++)
    {
      displs[i] = displs[i - 1] + counts[i - 1];
    }
    all.resize(displs.back() + counts.back());
  }
  MPI_Gatherv(data.data(), size, type, all.data(), counts.data(), displs.data(), type, 0,
              comm);
  return all;
}

std::ofstream OpenBinary(const std::filesystem::path &path)
{
  std::ofstream fo(path, std::ios::binary);
  MFEM_VERIFY(fo.is_open(), "Unable to open sensitivity export file " << path << "!");
  return fo;
}

// Write data in little-endian byte order.
template <typename T>
void WriteBinary(std::ofstream &fo, const std::vector<T> &data)
{
  if constexpr (std::endian::native == std::endian::little)
  {
    fo.write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T));
  }
  else
  {
    for (const T v : data)
    {
      auto bytes = std::bit_cast<std::array<char, sizeof(T)>>(v);
      std::reverse(bytes.begin(), bytes.end());
      fo.write(bytes.data(), sizeof(T));
    }
  }
  MFEM_VERIFY(fo.good(), "Failed to write a sensitivity export file!");
}

}  // namespace

SensitivityExport::SensitivityExport(const IoData &iodata, const Mesh &mesh)
  : iodata(iodata), comm(mesh.GetComm()), shape(iodata.sensitivity.shape.has_value()),
    space_dim(mesh.SpaceDimension())
{
  if (!shape)
  {
    return;
  }
  const auto &region = *iodata.sensitivity.shape;
  const auto &pmesh = mesh.Get();
  const auto &nodes = *pmesh.GetNodes();
  const auto *nodal_fespace =
      dynamic_cast<const mfem::ParFiniteElementSpace *>(nodes.FESpace());
  MFEM_VERIFY(nodal_fespace, "Sensitivity export requires a parallel mesh with nodes!");

  // The nodes of the elements with the given attributes.
  Vector in_attributes;
  if (!region.attributes.empty())
  {
    in_attributes = MarkNodes(
        *nodal_fespace, pmesh.GetNE(), [&](int e)
        { return std::ranges::binary_search(region.attributes, pmesh.GetAttribute(e)); },
        [&](int e, mfem::Array<int> &dofs) { nodal_fespace->GetElementDofs(e, dofs); });
  }

  // Owned true DOFs in the region. Constrained (hanging) nodes are not true DOFs.
  const double scale = iodata.units.GetMeshLengthRelativeScale();
  for (int i = 0; i < nodal_fespace->GetNDofs(); i++)
  {
    std::array<int, 3> tdof = {-1, -1, -1};
    std::array<double, 3> x = {0.0, 0.0, 0.0};
    for (int c = 0; c < space_dim; c++)
    {
      const int ldof = nodal_fespace->DofToVDof(i, c);
      tdof[c] = nodal_fespace->GetLocalTDofNumber(ldof);
      x[c] = nodes(ldof);
    }
    if (tdof[0] < 0 || (!region.attributes.empty() && in_attributes(tdof[0]) == 0.0))
    {
      continue;
    }
    if (region.bbmin)
    {
      bool inside = true;
      for (int c = 0; c < space_dim; c++)
      {
        inside = inside && x[c] >= (*region.bbmin)[c] && x[c] <= (*region.bbmax)[c];
      }
      if (!inside)
      {
        continue;
      }
    }
    tdofs.push_back(tdof);
    for (int c = 0; c < 3; c++)
    {
      points.push_back(x[c] * scale);
    }
  }

  // The boundary attributes of each exported node. The list of boundary attributes is
  // global, so every process takes part in each marking.
  if (region.boundary_tags)
  {
    boundary.resize(tdofs.size());
    for (const int attr : pmesh.bdr_attributes)
    {
      const auto on_boundary = MarkNodes(
          *nodal_fespace, pmesh.GetNBE(), [&](int e)
          { return pmesh.GetBdrAttribute(e) == attr; }, [&](int e, mfem::Array<int> &dofs)
          { nodal_fespace->GetBdrElementDofs(e, dofs); });
      for (std::size_t k = 0; k < tdofs.size(); k++)
      {
        if (on_boundary(tdofs[k][0]) > 0.0)
        {
          boundary[k].push_back(attr);
        }
      }
    }
  }

  long long int num_points = tdofs.size();
  Mpi::GlobalSum(1, &num_points, comm);
  if (!region.bbmin && region.attributes.empty() && num_points > large_export)
  {
    Mpi::Warning(comm,
                 "The sensitivity export covers the whole mesh ({:d} nodes): restrict it "
                 "with \"BoundingBoxMin\" and \"BoundingBoxMax\", or \"Attributes\"!\n",
                 num_points);
  }
}

void SensitivityExport::AddObjective(nlohmann::json entry, const Vector *g, double scale)
{
  if (shape)
  {
    MFEM_VERIFY(g,
                "Missing gradient for sensitivity objective " << entry.at("Index") << "!");
    const double s = scale / iodata.units.GetMeshLengthRelativeScale();
    const auto *h_g = g->HostRead();
    auto &gradient = gradients.emplace_back(3 * tdofs.size(), 0.0);
    for (std::size_t k = 0; k < tdofs.size(); k++)
    {
      for (int c = 0; c < space_dim; c++)
      {
        gradient[3 * k + c] = s * h_g[tdofs[k][c]];
      }
    }
    entry["Gradient"] = fmt::format("gradient-{:d}.bin", entry.at("Index").get<int>());
  }
  objectives.push_back(std::move(entry));
}

void SensitivityExport::Write(const std::filesystem::path &dir) const
{
  const auto all_points = GatherToRoot(points, MPI_DOUBLE, comm);
  std::vector<std::vector<double>> all_gradients;
  for (const auto &gradient : gradients)
  {
    all_gradients.push_back(GatherToRoot(gradient, MPI_DOUBLE, comm));
  }
  std::vector<int> counts, attributes;
  for (const auto &attr_list : boundary)
  {
    counts.push_back(static_cast<int>(attr_list.size()));
    attributes.insert(attributes.end(), attr_list.begin(), attr_list.end());
  }
  const auto all_counts = GatherToRoot(counts, MPI_INT, comm);
  const auto all_attributes = GatherToRoot(attributes, MPI_INT, comm);
  if (!Mpi::Root(comm))
  {
    return;
  }

  std::filesystem::create_directories(dir);
  nlohmann::json meta = {{"Format", "PalaceSensitivity"},
                         {"Version", 1},
                         {"ProblemType", ToString(iodata.problem.type)},
                         {"LengthUnit", iodata.model.L0},
                         {"NumPoints", all_points.size() / 3},
                         {"ByteOrder", "LittleEndian"}};
  if (shape)
  {
    auto fp = OpenBinary(dir / "points.bin");
    WriteBinary(fp, all_points);
    meta["Points"] = "points.bin";
    const auto &region = *iodata.sensitivity.shape;
    nlohmann::json json_region = nlohmann::json::object();
    if (region.bbmin)
    {
      const double scale = iodata.units.GetMeshLengthRelativeScale();
      std::array<double, 3> bbmin = *region.bbmin, bbmax = *region.bbmax;
      for (int c = 0; c < 3; c++)
      {
        bbmin[c] *= scale;
        bbmax[c] *= scale;
      }
      json_region["BoundingBoxMin"] = bbmin;
      json_region["BoundingBoxMax"] = bbmax;
    }
    if (!region.attributes.empty())
    {
      json_region["Attributes"] = region.attributes;
    }
    if (!json_region.empty())
    {
      meta["Region"] = json_region;
    }
    if (region.boundary_tags)
    {
      std::vector<std::int64_t> offsets(all_counts.size() + 1, 0);
      for (std::size_t k = 0; k < all_counts.size(); k++)
      {
        offsets[k + 1] = offsets[k] + all_counts[k];
      }
      auto fb = OpenBinary(dir / "boundary.bin");
      WriteBinary(fb, offsets);
      WriteBinary(fb,
                  std::vector<std::int32_t>(all_attributes.begin(), all_attributes.end()));
      meta["Boundary"] = "boundary.bin";
    }
    for (std::size_t i = 0; i < objectives.size(); i++)
    {
      auto fg = OpenBinary(dir / objectives[i].at("Gradient").get<std::string>());
      WriteBinary(fg, all_gradients[i]);
    }
  }
  meta["Objectives"] = objectives;

  std::ofstream fo(dir / "sensitivity.json");
  MFEM_VERIFY(fo.is_open(),
              "Unable to open sensitivity export file " << dir / "sensitivity.json" << "!");
  fo << meta.dump(2) << '\n';
  Mpi::Print(comm, " Wrote the sensitivities of {:d} objective{} to {}\n",
             objectives.size(), (objectives.size() > 1) ? "s" : "", dir.string());
}

}  // namespace palace
