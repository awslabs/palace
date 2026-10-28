// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_ATTRIBUTES_HPP
#define PALACE_MODELS_BOUNDARY_ATTRIBUTES_HPP

#include <set>
#include <string_view>
#include <mfem.hpp>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

// Boundary attributes present in a mesh, for checking attributes from the configuration.
class MeshBoundaryAttributes
{
  mfem::Array<int> marker;

public:
  explicit MeshBoundaryAttributes(const mfem::ParMesh &mesh)
    : marker(mesh::AttrToMarker(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0,
                                mesh.bdr_attributes))
  {
  }

  int Max() const { return marker.Size(); }
  bool Contains(int attr) const
  {
    return attr > 0 && attr <= marker.Size() && marker[attr - 1];
  }
};

// Warn once about configured attributes missing from the mesh and, when check_duplicates is
// true, abort on an attribute claimed twice. Each entry holds its list in .attributes.
template <typename Entries>
void CheckBoundaryAttributes(const MeshBoundaryAttributes &mesh_attrs,
                             const Entries &entries, std::string_view kind,
                             bool check_duplicates = true)
{
  std::set<int> bdr_warn_list;
  mfem::Array<int> marker;
  if (check_duplicates)
  {
    marker.SetSize(mesh_attrs.Max());
    marker = 0;
  }
  for (const auto &entry : entries)
  {
    for (auto attr : entry.attributes)
    {
      if (!mesh_attrs.Contains(attr))
      {
        bdr_warn_list.insert(attr);
        continue;
      }
      if (check_duplicates)
      {
        MFEM_VERIFY(!marker[attr - 1], "Multiple definitions of "
                                           << kind
                                           << " boundary properties for boundary attribute "
                                           << attr << "!");
        marker[attr - 1] = 1;
      }
    }
  }
  if (!bdr_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning("Unknown {} boundary attributes!\nSolver will just ignore them!", kind);
    utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
}

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_ATTRIBUTES_HPP
