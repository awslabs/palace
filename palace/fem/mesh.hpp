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

// XX TODO: Do we need to store qw * |J| separately in each quadrature data? Is it
//          significantly worse if we just use multiple inputs to the QFunction for the
//          different quantities?
// XX TODO: Can we skip adjugate storage and just compute from J on the fly? Or, probably
//          better, can we skip J storage and compute from adj(adj(J)/|J|) = J?
// XX TODO: Rename to CeedMeshData or something like that?

//
// Data structure for geometry information stored at quadrature points. Jacobian matrix is
// dim x space_dim, the adjugate is space_dim x dim, column-major storage by component.
//
struct CeedGeomFactorData_private
{
  // Dimension and space dimension for this element topology.
  int dim, space_dim;

  // Element indices from the mfem::Mesh used to construct Ceed objects with these geometry
  // factors.
  std::vector<int> indices;

  mfem::Vector wdetJ;  // qw * |J|, for H1 conformity with quadrature weights
  mfem::Vector adjJt;  // adj(J)^T / |J|, for H(curl) conformity
  mfem::Vector J;      // J / |J|, for H(div) conformity
  mfem::Vector attr;   // Mesh element attributes

  // Objects for libCEED interface to the quadrature data.
  CeedVector wdetJ_vec, adjJt_vec, J_vec, attr_vec;
  CeedElemRestriction wdetJ_restr, adjJt_restr, J_restr, attr_restr;
  Ceed ceed;

  CeedGeomFactorData_private(Ceed ceed)
    : dim(0), space_dim(0), wdetJ_vec(nullptr), adjJt_vec(nullptr), J_vec(nullptr),
      attr_vec(nullptr), wdetJ_restr(nullptr), adjJt_restr(nullptr), J_restr(nullptr),
      attr_restr(nullptr), ceed(ceed)
  {
  }
  ~CeedGeomFactorData_private();
};

using CeedGeomFactorData = std::unique_ptr<CeedGeomFactorData_private>;

}  // namespace ceed

//
// Wrapper for MFEM's Mesh class, with the addition of data structures for storing mesh
// geometry factors used to construct libCEED operators, as well as a few other useful
// data members not included in MFEM's Mesh or ParMesh types.
//
class Mesh
{
private:
  // Underlying MFEM object.
  std::unique_ptr<mfem::ParMesh> mesh;

  // Sequence to track mfem::Mesh::sequence and determine if geometry factors need updating.
  mutable long int sequence;

  // Attribute mapping for (global, 1-based) domain and boundary attributes to those on this
  // process (still 1-based). An entry of -1 indicates the global attribute is not used on
  // this process.
  mutable std::unordered_map<int, int> loc_attr, loc_bdr_attr;

  // Shared face mapping for boundary coefficients.
  mutable std::unordered_map<int, int> local_to_shared;

  // Mesh data structures for assembling libCEED operators on a (mixed) mesh:
  //   - Mesh element indices for threads and element geometry types.
  //   - Geometry factor quadrature point data (w |J|, adj(J)^T / |J|, J / |J|) for domain
  //     and boundary elements.
  //   - Attributes for domain and boundary elements. The attributes are not the same as the
  //     mesh element attributes, they correspond to the local attributes above converted
  //		 to 0-based indexing.
  mutable ceed::CeedObjectMap<ceed::CeedGeomFactorData> geom_data;

  void CheckSequence() const
  {
    if (sequence != mesh->GetSequence())
    {
      ClearData();
      sequence = mesh->GetSequence();
    }
  }

  std::unordered_map<int, int> &BuildAttributesGlobalToLocal(bool use_bdr = false) const;
  std::unordered_map<int, int> &BuildLocalToSharedFaceMap() const;
  ceed::CeedObjectMap<ceed::CeedGeomFactorData> &BuildCeedGeomFactorData() const;

public:
  template <typename... T>
  Mesh(T &&...args) : Mesh(std::make_unique<mfem::ParMesh>(std::forward<T>(args)...))
  {
  }

  template <typename T>
  Mesh(std::unique_ptr<T> &&mesh) : mesh(std::move(mesh)), sequence(mesh->GetSequence())
  {
    mesh->EnsureNodes();
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
  auto GetNumFaces() const { return Get().GetNumFaces(); }
  auto GetNV() const { return Get().GetNV(); }

  const auto &GetAttributeGlobalToLocal() const
  {
    CheckSequence();
    return !loc_attr.empty() ? loc_attr : BuildAttributesGlobalToLocal();
  }

  const auto &GetBdrAttributeGlobalToLocal() const
  {
    CheckSequence();
    return !loc_bdr_attr.empty() ? loc_bdr_attr : BuildAttributesGlobalToLocal(true);
  }

  const auto &GetLocalToSharedFaceMap() const
  {
    CheckSequence();
    return !local_to_shared.empty() ? local_to_shared : BuildLocalToSharedFaceMap();
  }

  const auto &GetCeedGeomFactorData() const
  {
    CheckSequence();
    return !geom_data.empty() ? geom_data : BuildCeedGeomFactorData();
  }

  void ClearData() const
  {
    loc_attr.clear();
    loc_bdr_attr.clear();
    local_to_shared.clear();
    geom_data.clear();
  }
  void ClearCeedGeomFactorData() const { geom_data.clear(); }

  MPI_Comm GetComm() const { return mesh->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_FEM_MESH_HPP
