// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

#include <algorithm>
#include <functional>
#include <ceed/backend.h>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/libceed/integrator.hpp"

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

auto BuildCeedAttributes(const mfem::ParMesh &mesh)
{
  // Set up sparse map from global domain attributes to local ones on this process.
  // Include ghost elements for all shared faces so we have their material properties
  // stored locally. New attributes for libCEED are contiguous and 1-based.
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
    mesh.GetSharedFaceTransformations(i, FET, T1, T2);
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
  // Handle empty mesh case - ensure loc_attr is never empty
  if (loc_attr.empty())
  {
    loc_attr[1] = 1;  // Add dummy attribute for empty mesh
  }
  return loc_attr;
}

auto BuildCeedBdrAttributes(const mfem::ParMesh &mesh)
{
  // Set up sparse map from global boundary attributes to local ones on this process. Each
  // original global boundary attribute maps to a key-value pairing of global domain
  // attributes which neighbor the given boundary and local boundary attributes. New
  // attributes for libCEED are contiguous and 1-based.
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
  // Handle empty mesh case - boundary attributes can remain empty
  return loc_bdr_attr;
}

auto GetElementIndices(const mfem::ParMesh &mesh, bool use_bdr, int start, int stop)
{
  // Count the number of elements of each type in the local mesh.
  std::unordered_map<mfem::Geometry::Type, int> counts;
  for (int i = start; i < stop; i++)
  {
    const auto geom = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto it = counts.find(geom);
    if (it == counts.end())
    {
      counts[geom] = 1;
    }
    else
    {
      it->second++;
    }
  }

  // Populate the indices arrays for each element geometry.
  std::unordered_map<mfem::Geometry::Type, int> offsets;
  std::unordered_map<mfem::Geometry::Type, std::vector<int>> element_indices;
  for (auto it = counts.begin(); it != counts.end(); ++it)
  {
    offsets[it->first] = 0;
    element_indices[it->first].resize(it->second);
  }
  for (int i = start; i < stop; i++)
  {
    const auto geom = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto &offset = offsets[geom];
    auto &indices = element_indices[geom];
    indices[offset++] = i;
  }

  return element_indices;
}

// Return a function mapping a local domain element index to its libCEED domain attribute.
// For a boundary submesh using its parent mesh's attribute maps, the parent's boundary
// attribute maps are used instead. The transformation objects are used as scratch space by
// the returned function and must outlive it.
std::function<int(int)> GetCeedDomainAttributeMap(
    const mfem::ParMesh &mesh, const std::unordered_map<int, int> &loc_attr,
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr,
    bool ceed_from_self, mfem::FaceElementTransformations &FET,
    mfem::IsoparametricTransformation &T1, mfem::IsoparametricTransformation &T2)
{
  if (!ceed_from_self)
  {
    if (const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(&mesh))
    {
      MFEM_VERIFY(submesh->GetFrom() == mfem::SubMesh::From::Boundary,
                  "Unexpected non-SubMesh object for BuildCeedGeomFactorData with Mesh "
                  "with (dim, space_dim) = ("
                      << mesh.Dimension() << ", " << mesh.SpaceDimension() << ")!");
      return [&mesh, &loc_bdr_attr, submesh, &FET, &T1, &T2](int i)
      {
        // Mesh is a boundary submesh with parent-based CEED data, so we use the
        // boundary attribute mappings from the parent mesh.
        const int attr = mesh.GetAttribute(i);
        const int nbr_attr = GetBdrNeighborAttribute(submesh->GetParentElementIDMap()[i],
                                                     *submesh->GetParent(), FET, T1, T2);
        MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                        loc_bdr_attr.at(attr).find(nbr_attr) != loc_bdr_attr.at(attr).end(),
                    "Missing libCEED boundary attribute for attribute " << attr << "!");
        return loc_bdr_attr.at(attr).at(nbr_attr);
      };
    }
  }
  // Non-submesh mesh, or submesh with self-rebuilt CEED data: use loc_attr directly.
  return [&mesh, &loc_attr](int i)
  {
    const int attr = mesh.GetAttribute(i);
    MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                "Missing libCEED domain attribute for attribute " << attr << "!");
    return loc_attr.at(attr);
  };
}

// Return a function mapping a local boundary element index to its libCEED boundary
// attribute. The transformation objects are used as scratch space by the returned function
// and must outlive it.
std::function<int(int)> GetCeedBdrAttributeMap(
    const mfem::ParMesh &mesh,
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr,
    mfem::FaceElementTransformations &FET, mfem::IsoparametricTransformation &T1,
    mfem::IsoparametricTransformation &T2)
{
  return [&mesh, &loc_bdr_attr, &FET, &T1, &T2](int i)
  {
    const int attr = mesh.GetBdrAttribute(i);
    const int nbr_attr = GetBdrNeighborAttribute(i, mesh, FET, T1, T2);
    MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                    loc_bdr_attr.at(attr).find(nbr_attr) != loc_bdr_attr.at(attr).end(),
                "Missing libCEED boundary attribute for attribute " << attr << "!");
    return loc_bdr_attr.at(attr).at(nbr_attr);
  };
}

auto GetActiveElementIndices(const std::vector<int> &indices,
                             const std::function<int(int)> &GetCeedAttribute,
                             const std::vector<int> &active_attr)
{
  // Store positions into the original geometry-type element list, preserving mesh order.
  std::vector<int> active_indices;
  active_indices.reserve(indices.size());
  for (std::size_t i = 0; i < indices.size(); i++)
  {
    if (std::binary_search(active_attr.begin(), active_attr.end(),
                           GetCeedAttribute(indices[i])))
    {
      active_indices.push_back(static_cast<int>(i));
    }
  }
  return active_indices;
}

auto AssembleGeometryData(Ceed ceed, mfem::Geometry::Type geom, std::vector<int> &indices,
                          std::vector<int> &active_indices,
                          const mfem::GridFunction &mesh_nodes, const Vector &elem_attr)
{
  const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
  const mfem::Mesh &mesh = *mesh_fespace.GetMesh();

  ceed::CeedGeomFactorData data;
  data.dim = mfem::Geometry::Dimension[geom];
  data.space_dim = mesh.SpaceDimension();
  data.indices = std::move(indices);
  data.active_indices = std::move(active_indices);
  const std::size_t num_elem = data.indices.size();
  if (!data.active_indices.empty() && data.active_indices.size() < num_elem)
  {
    // The positions the shared subset leaves out, in the same mesh order, for sub-operators
    // which cover exactly those elements. An empty subset, or one containing every element,
    // needs no complement since the full element list is used in both cases.
    data.complement_indices.reserve(num_elem - data.active_indices.size());
    auto next = data.active_indices.begin();
    for (std::size_t i = 0; i < num_elem; i++)
    {
      if (next != data.active_indices.end() && *next == static_cast<int>(i))
      {
        ++next;
      }
      else
      {
        data.complement_indices.push_back(static_cast<int>(i));
      }
    }
  }

  // Construct mesh node element restriction and basis.
  CeedElemRestriction mesh_restr =
      FiniteElementSpace::BuildCeedElemRestriction(mesh_fespace, ceed, geom, data.indices);
  CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(mesh_fespace, ceed, geom);
  CeedVector mesh_nodes_vec;
  ceed::InitCeedVector(mesh_nodes, ceed, &mesh_nodes_vec);
  CeedInt num_qpts;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &num_qpts));

  // Construct element attribute element restriction and basis.
  CeedElemRestriction attr_restr;
  CeedBasis attr_basis;
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, num_elem, 1, 1, num_elem,
                                                        CEED_STRIDES_BACKEND, &attr_restr));
  {
    // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
    mfem::Vector Bt(num_qpts), Gt(num_qpts), qX(num_qpts), qW(num_qpts);
    Bt = 1.0;
    Gt = 0.0;
    qX = 0.0;
    qW = 0.0;
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_qpts,
                                           Bt.GetData(), Gt.GetData(), qX.GetData(),
                                           qW.GetData(), &attr_basis));
  }
  CeedVector elem_attr_vec;
  ceed::InitCeedVector(elem_attr, ceed, &elem_attr_vec);

  // Allocate storage for geometry factor data (attribute + quadrature weight + Jacobian at
  // each quadrature point). The vector is only ever written and read by libCEED operators
  // on this Ceed context, so it is stored in the backend's own strided layout
  // (CEED_STRIDES_BACKEND, as for the assembled quadrature data). Backends whose E-vector
  // layout matches that layout (the GPU backends) then read the vector directly as the
  // E-vector of every operator instead of gathering a copy of it per operator at setup.
  CeedInt geom_data_size = 2 + data.space_dim * data.dim;
  const CeedSize geom_data_length =
      static_cast<CeedSize>(num_elem) * num_qpts * geom_data_size;
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, geom_data_length, &data.geom_data));
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                           ceed, num_elem, num_qpts, geom_data_size, geom_data_length,
                           CEED_STRIDES_BACKEND, &data.geom_data_restr));
  if (!data.active_indices.empty() && data.active_indices.size() < num_elem)
  {
    // Reuse the full geometry vector without copying any factors. Each offset addresses the
    // first component at one quadrature point of an active element in the full vector; the
    // component stride then selects the remaining components at that point. The vector
    // layout is the backend's, queried from the strided restriction (as [nodes, components,
    // elements] strides) rather than assumed.
    CeedInt layout[3];
    PalaceCeedCall(ceed, CeedElemRestrictionGetLLayout(data.geom_data_restr, layout));
    std::vector<CeedInt> active_offsets(data.active_indices.size() * num_qpts);
    for (std::size_t k = 0; k < data.active_indices.size(); k++)
    {
      const CeedInt elem_offset = data.active_indices[k] * layout[2];
      for (CeedInt q = 0; q < num_qpts; q++)
      {
        active_offsets[k * num_qpts + q] = elem_offset + q * layout[0];
      }
    }
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, data.active_indices.size(), num_qpts, geom_data_size,
                             layout[1], geom_data_length, CEED_MEM_HOST, CEED_COPY_VALUES,
                             active_offsets.data(), &data.active_geom_data_restr));
  }

  // Compute the required geometry factors at quadrature points.
  ceed::AssembleCeedGeometryData(ceed, mesh_restr, mesh_basis, mesh_nodes_vec, attr_restr,
                                 attr_basis, elem_attr_vec, data.geom_data,
                                 data.geom_data_restr);
  PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
  PalaceCeedCall(ceed, CeedVectorDestroy(&elem_attr_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&attr_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&attr_basis));

  return data;
}

auto BuildCeedGeomFactorData(
    const mfem::ParMesh &mesh, const std::unordered_map<int, int> &loc_attr,
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr,
    const std::vector<int> &lossy_attr, const std::vector<int> &active_bdr_attr, Ceed ceed,
    bool ceed_from_self)
{
  // Create a list of the element indices in the mesh corresponding to a given thread and
  // element geometry type and corresponding geometry factor data. libCEED operators will be
  // constructed in parallel over threads, where each thread builds a composite operator
  // with sub-operators for each geometry.
  const std::size_t nt = ceed::internal::NumCeeds();
  auto it = std::find(ceed::internal::GetCeedObjects().begin(),
                      ceed::internal::GetCeedObjects().end(), ceed);
  MFEM_VERIFY(it != ceed::internal::GetCeedObjects().end(),
              "Unable to find matching Ceed context in BuildCeedGeomFactorData!");
  std::size_t i = std::distance(ceed::internal::GetCeedObjects().begin(), it);
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  ceed::GeometryObjectMap<ceed::CeedGeomFactorData> geom_data_map;

  // First domain elements.
  {
    const int num_elem = mesh.GetNE();
    const int stride = (num_elem + nt - 1) / nt;
    const int start = i * stride;
    const int stop = std::min(start + stride, num_elem);
    constexpr bool use_bdr = false;
    auto GetCeedAttribute = GetCeedDomainAttributeMap(mesh, loc_attr, loc_bdr_attr,
                                                      ceed_from_self, FET, T1, T2);
    auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
    for (auto &[geom, indices] : element_indices)
    {
      auto active_indices = GetActiveElementIndices(indices, GetCeedAttribute, lossy_attr);
      Vector elem_attr(indices.size());
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = GetCeedAttribute(indices[k]);
      }
      geom_data_map.emplace(geom, AssembleGeometryData(ceed, geom, indices, active_indices,
                                                       *mesh.GetNodes(), elem_attr));
    }
  }

  // Then boundary elements. For embedded meshes (dim != sdim), boundary element geometry
  // data is only available when the CEED data has been rebuilt from the mesh itself
  // (after attribute remapping on a boundary submesh).
  if (mesh.Dimension() == mesh.SpaceDimension() || ceed_from_self)
  {
    const int nbe = mesh.GetNBE();
    const int stride = (nbe + nt - 1) / nt;
    const int start = i * stride;
    const int stop = std::min(start + stride, nbe);
    constexpr bool use_bdr = true;
    auto GetCeedAttribute = GetCeedBdrAttributeMap(mesh, loc_bdr_attr, FET, T1, T2);
    auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
    for (auto &[geom, indices] : element_indices)
    {
      auto active_indices =
          GetActiveElementIndices(indices, GetCeedAttribute, active_bdr_attr);
      Vector elem_attr(indices.size());
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = GetCeedAttribute(indices[k]);
      }
      geom_data_map.emplace(geom, AssembleGeometryData(ceed, geom, indices, active_indices,
                                                       *mesh.GetNodes(), elem_attr));
    }
  }

  return geom_data_map;
}

}  // namespace

const ceed::GeometryObjectMap<ceed::CeedGeomFactorData> &
Mesh::GetCeedGeomFactorData(Ceed ceed) const
{
  auto it = geom_data.find(ceed);
  MFEM_ASSERT(it != geom_data.end(), "Unknown Ceed context in GetCeedGeomFactorData!");
  auto &geom_data_map = it->second;
  if (geom_data_map.empty() && !loc_attr.empty())
  {
    geom_data_map = BuildCeedGeomFactorData(*mesh, loc_attr, loc_bdr_attr, lossy_attr,
                                            active_bdr_attr, ceed, ceed_from_self);
  }
  return geom_data_map;
}

void Mesh::SetCeedActiveAttributes(std::vector<int> lossy_attr_,
                                   std::vector<int> active_bdr_attr_)
{
  auto Normalize = [](std::vector<int> &attr, std::size_t max_attr)
  {
    std::sort(attr.begin(), attr.end());
    attr.erase(std::unique(attr.begin(), attr.end()), attr.end());
    MFEM_VERIFY(attr.empty() || (attr.front() > 0 && attr.back() <= max_attr),
                "Invalid process-local libCEED element-set attribute!");
  };
  Normalize(lossy_attr_, MaxCeedAttribute());
  Normalize(active_bdr_attr_, MaxCeedBdrAttribute());
  if (lossy_attr_ == lossy_attr && active_bdr_attr_ == active_bdr_attr)
  {
    // Reconfiguring with the same sets is a no-op. This happens for the mesh levels
    // retained across adaptive refinement iterations, whose geometry data stays valid.
    return;
  }
  // The element subsets are baked into the geometry data and into the finite element
  // restrictions built from it, so the sets can only change before any geometry data exists
  // (the sets depend only on the configuration, so this never happens within a run).
  for (const auto &[ceed, geom_data_map] : geom_data)
  {
    MFEM_VERIFY(geom_data_map.empty(),
                "Cannot change the libCEED element sets after geometry data is built!");
  }
  lossy_attr = std::move(lossy_attr_);
  active_bdr_attr = std::move(active_bdr_attr_);
}

void Mesh::ResetCeedObjects()
{
  for (auto &[ceed, geom_data_map] : geom_data)
  {
    for (auto &[key, val] : geom_data_map)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&val.geom_data));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val.geom_data_restr));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val.active_geom_data_restr));
    }
  }
  geom_data.clear();
  for (std::size_t i = 0; i < ceed::internal::GetCeedObjects().size(); i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    geom_data.emplace(ceed, ceed::GeometryObjectMap<ceed::CeedGeomFactorData>());
  }
}

void Mesh::Update()
{
  // Attribute mappings, etc. are always constructed for the parent mesh (use boundary
  // attribute maps for the domain attributes of a boundary submesh, for example).
  ceed_from_self = false;
  auto &parent_mesh = GetParentMesh(*mesh);
  parent_mesh.ExchangeFaceNbrData();
  loc_attr.clear();
  loc_bdr_attr.clear();
  loc_attr = BuildCeedAttributes(parent_mesh);
  loc_bdr_attr = BuildCeedBdrAttributes(parent_mesh);
  lossy_attr.clear();
  active_bdr_attr.clear();
  ResetCeedObjects();
}

void Mesh::RebuildCeedAttributes()
{
  // Rebuild CEED attribute maps from the mesh's own elements (not the parent). This is
  // needed after remapping attributes on a boundary submesh. For ranks with no local
  // elements (e.g., wave port submesh on ranks without port elements), the maps are left
  // empty. Sets the ceed_from_self flag so that BuildCeedGeomFactorData uses loc_attr for
  // domain elements (instead of the ParSubMesh-specific loc_bdr_attr path) and builds
  // boundary element geometry data (even when dim != sdim).
  loc_attr.clear();
  loc_bdr_attr.clear();
  mesh->ExchangeFaceNbrData();  // MPI collective — all ranks must participate
  if (mesh->GetNE() > 0)
  {
    loc_attr = BuildCeedAttributes(*mesh);
    loc_bdr_attr = BuildCeedBdrAttributes(*mesh);
  }
  ceed_from_self = true;
  lossy_attr.clear();
  active_bdr_attr.clear();
  ResetCeedObjects();
}

}  // namespace palace
