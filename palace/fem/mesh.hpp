// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MESH_HPP
#define PALACE_FEM_MESH_HPP

#include <memory>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace palace
{

namespace ceed
{

//
// Data structure for geometry information stored at quadrature points.
//
struct CeedGeomFactorData
{
  // Dimension of this element topology and space dimension of the underlying mesh.
  int dim, space_dim;

  // Domain or boundary indices from the mesh used to construct Ceed objects with these
  // geometry factors.
  std::vector<int> indices;

  // Mesh geometry factor data: {attr, w * |J|, adj(J)^T / |J|}. Jacobian matrix is
  // space_dim x dim, stored column-major by component.
  CeedVector geom_data;

  // Element restriction for the geometry factor quadrature data.
  CeedElemRestriction geom_data_restr;
};

}  // namespace ceed

//
// Wrapper for MFEM's ParMesh class, with extensions for Palace.
//
class FiniteElementSpace;

class Mesh
{
  friend class FiniteElementSpace;

private:
  // Underlying MFEM object (can also point to a derived class of mfem::ParMesh, such as
  // mfem::ParSubMesh).
  std::unique_ptr<mfem::ParMesh> mesh;

  // Attribute mapping for (global, MFEM, 1-based) domain and boundary attributes to those
  // for libCEED (local to this process, contiguous, also 1-based). For boundaries, the
  // inner map is a mapping from neighboring MFEM domain attribute to the resulting local
  // boundary attribute (to discern boundary elements of a given attribute which border more
  // than one domain). Interior boundaries use as neighbor the element with the smaller
  // domain attribute in order to be consistent when the interior boundary element normals
  // are not aligned.
  std::unordered_map<int, int> loc_attr;
  std::unordered_map<int, std::unordered_map<int, int>> loc_bdr_attr;

  // Mesh data structures for assembling libCEED operators on a (mixed) mesh:
  //   - Mesh element indices for threads and element geometry types.
  //   - Attributes for domain and boundary elements. The attributes are not the same as the
  //     MFEM mesh element attributes, they correspond to the local, contiguous (1-based)
  //     attributes above.
  //   - Geometry factor quadrature point data (w |J| and adj(J)^T / |J|) for domain and
  //     boundary elements.
  mutable ceed::CeedObjectMap<ceed::CeedGeomFactorData> geom_data;

public:
  template <typename... T>
  Mesh(T &&...args) : Mesh(std::make_unique<mfem::ParMesh>(std::forward<T>(args)...))
  {
  }
  template <typename T>
  Mesh(std::unique_ptr<T> &&mesh) : mesh(std::move(mesh))
  {
    this->mesh->EnsureNodes();
    Update();
  }
  ~Mesh() { ResetCeedObjects(); }

  const auto &Get() const { return *mesh; }
  auto &Get() { return *mesh; }


  auto Dimension() const { return Get().Dimension(); }
  auto SpaceDimension() const { return Get().SpaceDimension(); }
  auto GetNE() const { return Get().GetNE(); }
  auto GetNBE() const { return Get().GetNBE(); }

  // Attribute queries.
  auto GetAttribute(int i) const { return Get().GetAttribute(i); }
  auto GetBdrAttribute(int i) const { return Get().GetBdrAttribute(i); }
  int MaxAttribute() const { return Get().attributes.Size() ? Get().attributes.Max() : 0; }
  int MaxBdrAttribute() const
  {
    return Get().bdr_attributes.Size() ? Get().bdr_attributes.Max() : 0;
  }
  auto NumAttributes() const { return Get().attributes.Size(); }
  auto NumBdrAttributes() const { return Get().bdr_attributes.Size(); }
  const auto &Attributes() const { return Get().attributes; }
  const auto &BdrAttributes() const { return Get().bdr_attributes; }

  // Geometry queries.
  auto GetGlobalNE() const { return Get().GetGlobalNE(); }
  auto GetElementGeometry(int i) const { return Get().GetElementGeometry(i); }
  auto GetBdrElementGeometry(int i) const { return Get().GetBdrElementGeometry(i); }
  auto Nonconforming() const { return Get().Nonconforming(); }
  auto Conforming() const { return Get().Conforming(); }
  auto GetNV() const { return Get().GetNV(); }

  // Topology queries.
  auto GetNumFaces() const { return Get().GetNumFaces(); }
  void GetFaceElements(int f, int *e1, int *e2) const { Get().GetFaceElements(f, e1, e2); }

  // High-order nodes.
  auto *GetNodes() const { return Get().GetNodes(); }
  auto *GetNodes() { return Get().GetNodes(); }

  // Element transformations. The overloads writing into a caller-provided T are const
  // because they don't modify the mesh itself.
  auto *GetElementTransformation(int i) { return Get().GetElementTransformation(i); }
  void GetElementTransformation(int i, mfem::IsoparametricTransformation *T) const
  {
    Get().GetElementTransformation(i, T);
  }
  auto *GetBdrElementTransformation(int i) { return Get().GetBdrElementTransformation(i); }
  void GetBdrElementTransformation(int i, mfem::IsoparametricTransformation *T) const
  {
    Get().GetBdrElementTransformation(i, T);
  }
  void GetFaceElementTransformations(int f, mfem::FaceElementTransformations &FET,
                                     mfem::IsoparametricTransformation &T1,
                                     mfem::IsoparametricTransformation &T2)
  {
    Get().GetFaceElementTransformations(f, FET, T1, T2);
  }
  void GetSharedFaceTransformations(int i, mfem::FaceElementTransformations &FET,
                                    mfem::IsoparametricTransformation &T1,
                                    mfem::IsoparametricTransformation &T2)
  {
    Get().GetSharedFaceTransformations(i, FET, T1, T2);
  }

  // Parallel queries.
  auto GetNSharedFaces() const { return Get().GetNSharedFaces(); }
  void ExchangeFaceNbrData() { Get().ExchangeFaceNbrData(); }

  // Refinement.
  void GeneralRefinement(const mfem::Array<int> &m, int nonconf = -1, int nc_limit = 0)
  {
    Get().GeneralRefinement(m, nonconf, nc_limit);
  }
  void GeneralRefinement(const mfem::Array<mfem::Refinement> &m, int nonconf = -1,
                         int nc_limit = 0)
  {
    Get().GeneralRefinement(m, nonconf, nc_limit);
  }

  // Derefinement (encapsulates pncmesh access).
  const mfem::Table &GetDerefinementTable() const;
  void SynchronizeDerefinementData(mfem::Array<double> &elem_error,
                                   const mfem::Table &deref_table) const;

  // Mesh replacement (for rebalance).
  void Reset(std::unique_ptr<mfem::ParMesh> new_mesh);

  const auto &GetCeedAttributes() const { return loc_attr; }
  const auto &GetCeedBdrAttributes() const { return loc_bdr_attr; }

  // Convert a list of global attributes to the corresponding process-local libCEED ones.
  template <typename T>
  auto GetCeedAttributes(const T &attr_list) const
  {
    // Skip any entries in the input global attribute list which are not local to this
    // process.
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

  // Convert a list of global boundary attributes to the corresponding process-local libCEED
  // ones.
  template <typename T>
  auto GetCeedBdrAttributes(const T &attr_list) const
  {
    // Skip any entries in the input global boundary attribute list which are not local
    // to this process.
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

  auto GetCeedAttributes(const int attr) const
  {
    return GetCeedAttributes(std::vector<int>{attr});
  }

  auto GetCeedBdrAttributes(const int attr) const
  {
    return GetCeedBdrAttributes(std::vector<int>{attr});
  }

  auto MaxCeedAttribute() const { return GetCeedAttributes().size(); }
  auto MaxCeedBdrAttribute() const
  {
    std::size_t bdr_attr_max = 0;
    for (const auto &[attr, bdr_attr_map] : GetCeedBdrAttributes())
    {
      bdr_attr_max += bdr_attr_map.size();
    }
    return bdr_attr_max;
  }

  const ceed::GeometryObjectMap<ceed::CeedGeomFactorData> &
  GetCeedGeomFactorData(Ceed ceed) const;

  void ResetCeedObjects();

  void Update();

  MPI_Comm GetComm() const { return mesh->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_FEM_MESH_HPP
