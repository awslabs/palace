// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MESH_HPP
#define PALACE_FEM_MESH_HPP

#include <memory>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>

namespace palace
{

//
// Wrapper for MFEM's ParMesh class, with extensions for Palace.
//
class Mesh
{
private:
  // Underlying MFEM object (can also point to a derived class of mfem::ParMesh, such as
  // mfem::ParSubMesh).
  mutable std::unique_ptr<mfem::ParMesh> mesh;

  // Sequence to track mfem::Mesh::sequence and determine if geometry factors need updating.
  mutable long int sequence;

  // Attribute mapping for (global, 1-based) domain and boundary attributes to those on this
  // process (still 1-based). For boundaries, the inner map is a mapping from neighboring
  // domain attribute to the resulting local boundary attribute (to discern boundary
  // elements with global boundary attribute which borders more than one domain). Interior
  // boundaries use as neighbor the element with the smaller domain attribute in order to
  // be consistent when the interior boundary element normals are not aligned.
  mutable std::unordered_map<int, int> loc_attr;
  mutable std::unordered_map<int, std::unordered_map<int, int>> loc_bdr_attr;

  void CheckSequenceRebuild() const
  {
    if (sequence != mesh->GetSequence())
    {
      Rebuild();
      sequence = mesh->GetSequence();
    }
  }
  void Rebuild() const;

public:
  template <typename T>
  Mesh(std::unique_ptr<T> &&mesh) : mesh(std::move(mesh))
  {
    this->mesh->EnsureNodes();
    Rebuild();
    sequence = this->mesh->GetSequence();
  }

  template <typename... T>
  Mesh(T &&...args) : Mesh(std::make_unique<mfem::ParMesh>(std::forward<T>(args)...))
  {
  }

  const auto &Get() const { return *mesh; }
  auto &Get() { return *mesh; }

  operator const mfem::ParMesh &() const { return Get(); }
  operator mfem::ParMesh &() { return Get(); }

  operator const std::unique_ptr<mfem::ParMesh> &() const { return mesh; }
  operator std::unique_ptr<mfem::ParMesh> &() { return mesh; }

  auto Dimension() const { return Get().Dimension(); }
  auto SpaceDimension() const { return Get().SpaceDimension(); }
  auto GetNE() const { return Get().GetNE(); }
  auto GetNBE() const { return Get().GetNBE(); }

  const auto &GetAttributeGlobalToLocal() const
  {
    CheckSequenceRebuild();
    return loc_attr;
  }

  const auto &GetBdrAttributeGlobalToLocal() const
  {
    CheckSequenceRebuild();
    return loc_bdr_attr;
  }

  template <typename T>
  auto GetAttributeGlobalToLocal(const T &attr_list) const
  {
    // Skip any entries in the input global attribute list which are not on local to this
    // process.
    const auto &loc_attr = GetAttributeGlobalToLocal();
    mfem::Array<int> loc_attr_list;
    for (auto attr : attr_list)
    {
      if (loc_attr.find(attr) != loc_attr.end())
      {
        loc_attr_list.Append(loc_attr.at(attr));
      }
    }
    return loc_attr_list;
  }

  template <typename T>
  auto GetBdrAttributeGlobalToLocal(const T &attr_list) const
  {
    // Skip any entries in the input global boundary attribute list which are not on local
    // to this process.
    const auto &loc_bdr_attr = GetBdrAttributeGlobalToLocal();
    mfem::Array<int> loc_attr_list;
    for (auto attr : attr_list)
    {
      if (loc_bdr_attr.find(attr) != loc_bdr_attr.end())
      {
        const auto &bdr_attr_map = loc_bdr_attr.at(attr);
        for (auto it = bdr_attr_map.begin(); it != bdr_attr_map.end(); ++it)
        {
          loc_attr_list.Append(it->second);
        }
      }
    }
    return loc_attr_list;
  }

  auto GetAttributeGlobalToLocal(const int attr) const
  {
    return GetAttributeGlobalToLocal(std::vector<int>{attr});
  }

  auto GetBdrAttributeGlobalToLocal(const int attr) const
  {
    return GetBdrAttributeGlobalToLocal(std::vector<int>{attr});
  }

  int GetAttributeGlobalToLocal(const mfem::ElementTransformation &T) const;

  MPI_Comm GetComm() const { return mesh->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_FEM_MESH_HPP
