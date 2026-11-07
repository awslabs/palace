// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_MESH_HPP
#define PALACE_FEM_MESH_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>
#include <mfem.hpp>
#include "fem/brokenspace.hpp"
#include "fem/libceed/ceed.hpp"

namespace palace
{

namespace ceed
{

// The element lists which may be used to assemble a libCEED sub-operator.
enum class CeedElementSubset
{
  Full,
  Active
};

//
// Data structure for geometry information stored at quadrature points.
//
struct CeedGeomFactorData
{
  // Dimension of this element topology and space dimension of the underlying mesh.
  int dim, space_dim;

  // Domain or boundary indices from the mesh used to construct Ceed objects with these
  // geometry factors, in the original mesh order. active_indices contains positions in
  // indices for the shared lossy-domain or active-boundary subset, in the same order (empty
  // when no subset is configured, full-sized when every element is active).
  std::vector<int> indices, active_indices;

  // Mesh geometry factor data: {attr, w * |J|, adj(J)^T / |J|}. Jacobian matrix is
  // space_dim x dim, stored column-major by component. For 3D domain elements, this is
  // followed by the physical coordinates x of the quadrature points for the geometry
  // factor data restricted to elements of given attributes (see
  // Mesh::GetCeedGeomFactorData). The vector uses the libCEED backend's strided layout
  // (CEED_STRIDES_BACKEND) and is only accessed through the restrictions below.
  CeedVector geom_data;

  // Element restrictions for all geometry factors and the shared subset. Both address the
  // same geom_data vector. active_geom_data_restr is null when the subset is empty or
  // contains every element, in which case the full strided restriction is used.
  CeedElemRestriction geom_data_restr;
  CeedElemRestriction active_geom_data_restr = nullptr;

  CeedElemRestriction GetGeomDataRestriction(CeedElementSubset subset) const
  {
    if (subset == CeedElementSubset::Active && active_indices.size() < indices.size())
    {
      MFEM_ASSERT(active_geom_data_restr, "Empty active subset for geometry factor data!");
      return active_geom_data_restr;
    }
    return geom_data_restr;
  }
};

}  // namespace ceed

//
// Wrapper for MFEM's ParMesh class, with extensions for Palace.
//
class Mesh
{
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
  bool ceed_from_self = false;  // True after RebuildCeedAttributes()

  // Process-local libCEED attributes in the two shared element sets: lossy domain elements
  // and active boundary elements.
  std::vector<int> lossy_attr, active_bdr_attr;

  // Mesh data structures for assembling libCEED operators on a (mixed) mesh:
  //   - Mesh element indices for threads and element geometry types.
  //   - Attributes for domain and boundary elements. The attributes are not the same as the
  //     MFEM mesh element attributes, they correspond to the local, contiguous (1-based)
  //     attributes above.
  //   - Geometry factor quadrature point data (w |J| and adj(J)^T / |J|) for domain and
  //     boundary elements.
  mutable ceed::CeedObjectMap<ceed::CeedGeomFactorData> geom_data;

  // Interior boundary sides for broken spaces (see fem/brokenspace.hpp), cached per sorted
  // attribute list (error estimator, broken solution space) and valid for the mesh object
  // and sequence they were computed for. On a mesh with hanging entities they are inherited
  // from the sides of the root elements, computed once per list.
  struct CrackSidesEntry
  {
    CrackSides sides;
    const mfem::ParMesh *mesh = nullptr;
    long sequence = -1;
    RootCrackSides root_sides;
    const mfem::ParMesh *root_mesh = nullptr;
    bool root_valid = false;
  };
  mutable std::map<std::vector<int>, CrackSidesEntry> crack_sides;

  // Geometry factor data of the domain elements with given (global) domain attributes,
  // with the physical coordinates of the quadrature points for 3D elements. Built on demand
  // for each Ceed context.
  mutable std::map<
      std::vector<int>,
      std::unordered_map<Ceed,
                         std::optional<ceed::GeometryObjectMap<ceed::CeedGeomFactorData>>>>
      geom_data_attr;

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

  operator const mfem::ParMesh &() const { return Get(); }
  operator mfem::ParMesh &() { return Get(); }

  operator const std::unique_ptr<mfem::ParMesh> &() const { return mesh; }
  operator std::unique_ptr<mfem::ParMesh> &() { return mesh; }

  auto Dimension() const { return Get().Dimension(); }
  auto SpaceDimension() const { return Get().SpaceDimension(); }
  auto GetNE() const { return Get().GetNE(); }
  auto GetNBE() const { return Get().GetNBE(); }

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
      if (loc_attr.contains(attr))
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
      if (loc_bdr_attr.contains(attr))
      {
        const auto &bdr_attr_map = loc_bdr_attr.at(attr);
        for (const auto &[nbr_attr, ceed_attr] : bdr_attr_map)
        {
          loc_attr_list.Append(ceed_attr);
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

  // Return the interior boundary sides of the local elements for the interior boundaries
  // formed by the given boundary attributes (see fem/brokenspace.hpp), computed on first
  // use for each list of attributes. The reference stays valid for the lifetime of the
  // mesh object, and its sides are updated in place. Collective.
  const CrackSides &GetCrackSides(const std::vector<int> &attr_list) const;

  // Set the interior boundary sides from those of the coarse mesh this mesh was refined
  // from (its last operation), so the side labels of the two meshes agree, as the
  // h-prolongation between broken spaces requires. Collective.
  void InheritCrackSides(const std::vector<int> &attr_list, const Mesh &coarse_mesh) const;

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

  // Geometry factor data of the domain elements with the given (global, sorted) domain
  // attributes only, for integrators restricted to these elements (see
  // BilinearForm::AddDomainIntegratorOnAttributes). For 3D elements, it also stores the
  // physical coordinates of the quadrature points, for QFunctions with spatially varying
  // coefficients (such as the PML material tensors). The shared element subsets do not
  // apply to this data (its active subset is empty).
  const ceed::GeometryObjectMap<ceed::CeedGeomFactorData> &
  GetCeedGeomFactorData(Ceed ceed, const std::vector<int> &attr_list) const;

  // Configure the process-local libCEED attributes in the shared lossy-domain and
  // active-boundary element sets. The sets must be configured before geometry data is
  // requested; reconfiguring with the same sets afterwards is a no-op (mesh levels retained
  // across adaptive refinement iterations are reconfigured by every new SpaceOperator).
  void SetCeedActiveAttributes(std::vector<int> lossy_attr_,
                               std::vector<int> active_bdr_attr_);

  const std::vector<int> &GetCeedActiveAttributes(bool use_bdr) const
  {
    return use_bdr ? active_bdr_attr : lossy_attr;
  }

  void ResetCeedObjects();

  void Update();

  // Rebuild CEED attribute maps from the mesh's own elements, ignoring any parent mesh.
  // Use after remapping attributes on a boundary submesh (RemapSubMeshAttributes +
  // RemapSubMeshBdrAttributes) so that CEED data reflects the remapped attributes.
  void RebuildCeedAttributes();

  MPI_Comm GetComm() const { return mesh->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_FEM_MESH_HPP
