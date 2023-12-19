// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"

namespace palace
{

namespace
{

const auto &GetParentMesh(const mfem::ParMesh &mesh)
{
  // Get the parent mesh if the mesh is a boundary submesh (no submesh of submesh
  // capabilities, for now).
  const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(&mesh);
  if (submesh && submesh->GetFrom() == mfem::SubMesh::From::Boundary)
  {
    return *submesh->GetParent();
  }
  return mesh;
}

auto &GetParentMesh(mfem::ParMesh &mesh)
{
  return const_cast<mfem::ParMesh &>(
      GetParentMesh(const_cast<const mfem::ParMesh &>(mesh)));
}

auto BuildAttributeGlobalToLocal(const mfem::ParMesh &mesh)
{
  // Set up sparse map from global domain attributes to local ones on this process.
  // Include ghost elements for all shared faces so we have their material properties
  // stored locally.
  std::unordered_map<int, int> loc_attr;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  int count = 0;
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    const int attr = mesh.GetAttribute(i);
    if (loc_attr.find(attr) == loc_attr.end())
    {
      loc_attr[attr] = ++count;
    }
  }
  for (int i = 0; i < mesh.GetNSharedFaces(); i++)
  {
    mesh.GetSharedFaceTransformations(i, &FET, &T1, &T2);
    int attr = FET.Elem1->Attribute;
    if (loc_attr.find(attr) == loc_attr.end())
    {
      loc_attr[attr] = ++count;
    }
    attr = FET.Elem2->Attribute;
    if (loc_attr.find(attr) == loc_attr.end())
    {
      loc_attr[attr] = ++count;
    }
  }
  return loc_attr;
}

auto GetBdrNeighborAttribute(int i, const mfem::ParMesh &mesh,
                             mfem::FaceElementTransformations &FET,
                             mfem::IsoparametricTransformation &T1,
                             mfem::IsoparametricTransformation &T2)
{
  // For internal boundaries, use the element which corresponds to the domain with lower
  // attribute number (ensures all boundary elements are aligned).
  BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(i, mesh, FET, T1, T2);
  return (FET.Elem2 && FET.Elem2->Attribute < FET.Elem1->Attribute) ? FET.Elem2->Attribute
                                                                    : FET.Elem1->Attribute;
}

auto BuildBdrAttributeGlobalToLocal(const mfem::ParMesh &mesh)
{
  // Set up sparse map from global boundary attributes to local ones on this process. Each
  // original global boundary attribute maps to a key-value pairing of global domain
  // attributes which neighbor the given boundary and local boundary attributes.
  std::unordered_map<int, std::unordered_map<int, int>> loc_bdr_attr;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  int count = 0;
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    const int attr = mesh.GetBdrAttribute(i);
    const int nbr_attr = GetBdrNeighborAttribute(i, mesh, FET, T1, T2);
    auto &bdr_attr_map = loc_bdr_attr[attr];
    if (bdr_attr_map.find(nbr_attr) == bdr_attr_map.end())
    {
      bdr_attr_map[nbr_attr] = ++count;
    }
  }
  return loc_bdr_attr;
}

}  // namespace

void Mesh::Rebuild() const
{
  // Attribute mappings, etc. are always constructed for the parent mesh (use boundary
  // attribute maps for the domain attributes of a boundary submesh, for example).
  auto &parent_mesh = GetParentMesh(*mesh);
  parent_mesh.ExchangeFaceNbrData();
  loc_attr.clear();
  loc_bdr_attr.clear();
  loc_attr = BuildAttributeGlobalToLocal(parent_mesh);
  loc_bdr_attr = BuildBdrAttributeGlobalToLocal(parent_mesh);
}

int Mesh::GetAttributeGlobalToLocal(const mfem::ElementTransformation &T) const
{
  if (T.GetDimension() == T.GetSpaceDim())
  {
    // Domain element.
    auto it = loc_attr.find(T.Attribute);
    MFEM_ASSERT(it != loc_attr.end(), "Invalid domain attribute " << T.Attribute << "!");
    return it->second;
  }
  else
  {
    // Boundary element (or boundary submesh domain).
    auto bdr_attr_map = loc_bdr_attr.find(T.Attribute);
    MFEM_ASSERT(bdr_attr_map != loc_bdr_attr.end(),
                "Invalid domain attribute " << T.Attribute << "!");
    const int nbr_attr = [&]()
    {
      mfem::FaceElementTransformations FET;  // XX TODO: Preallocate these for all elements
      mfem::IsoparametricTransformation T1, T2;
      if (const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(T.mesh))
      {
        MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                    "Unexpected element type in GetAttributeGlobalToLocal!");
        return GetBdrNeighborAttribute(submesh->GetParentElementIDMap()[T.ElementNo],
                                       *submesh->GetParent(), FET, T1, T2);
      }
      else
      {
        MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                    "Unexpected element type in GetAttributeGlobalToLocal!");
        return GetBdrNeighborAttribute(
            T.ElementNo, *static_cast<const mfem::ParMesh *>(T.mesh), FET, T1, T2);
      }
    }();
    auto it = bdr_attr_map->second.find(nbr_attr);
    MFEM_ASSERT(it != bdr_attr_map->second.end(),
                "Invalid domain attribute " << nbr_attr << "!");
    return it->second;
  }
}

}  // namespace palace
