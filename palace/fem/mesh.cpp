// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

#include <algorithm>
#include <cmath>
#include "fem/brokenspace.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/libceed/integrator.hpp"
#include "utils/communication.hpp"

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
  // stored locally. New attributes for libCEED are contiguous and 1-based. Ghost elements
  // are read from the face-neighbor data (the shared face transformations crash for
  // nonconforming domain submeshes).
  std::unordered_map<int, int> loc_attr;
  int count = 0;
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    const int attr = mesh.GetAttribute(i);
    if (!loc_attr.contains(attr))
    {
      loc_attr[attr] = ++count;
    }
  }
  for (const mfem::Element *el : mesh.face_nbr_elements)
  {
    const int attr = el->GetAttribute();
    if (!loc_attr.contains(attr))
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
    if (!bdr_attr_map.contains(nbr_attr))
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
  for (const auto &[geom, count] : counts)
  {
    offsets[geom] = 0;
    element_indices[geom].resize(count);
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

auto AssembleGeometryData(Ceed ceed, mfem::Geometry::Type geom, std::vector<int> &indices,
                          const mfem::GridFunction &mesh_nodes, const Vector &elem_attr)
{
  const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
  const mfem::Mesh &mesh = *mesh_fespace.GetMesh();

  ceed::CeedGeomFactorData data;
  data.dim = mfem::Geometry::Dimension[geom];
  data.space_dim = mesh.SpaceDimension();
  data.indices = std::move(indices);
  const std::size_t num_elem = data.indices.size();

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

  // Allocate storage for geometry factor data (stored as attribute + quadrature weight +
  // Jacobian, column-major).
  CeedInt geom_data_size = 2 + data.space_dim * data.dim;
  PalaceCeedCall(
      ceed,
      CeedVectorCreate(ceed, static_cast<CeedSize>(num_elem) * num_qpts * geom_data_size,
                       &data.geom_data));
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                           ceed, num_elem, num_qpts, geom_data_size,
                           static_cast<CeedSize>(num_elem) * num_qpts * geom_data_size,
                           CEED_STRIDES_BACKEND, &data.geom_data_restr));

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
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr, Ceed ceed,
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
    auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
    auto GetCeedAttribute = [&]() -> std::function<int(int)>
    {
      if (!ceed_from_self)
      {
        if (const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(&mesh))
        {
          MFEM_VERIFY(submesh->GetFrom() == mfem::SubMesh::From::Boundary,
                      "Unexpected non-SubMesh object for BuildCeedGeomFactorData with Mesh "
                      "with (dim, space_dim) = ("
                          << mesh.Dimension() << ", " << mesh.SpaceDimension() << ")!");
          return [&, submesh](int i)
          {
            // Mesh is a boundary submesh with parent-based CEED data, so we use the
            // boundary attribute mappings from the parent mesh.
            const int attr = mesh.GetAttribute(i);
            const int nbr_attr = GetBdrNeighborAttribute(
                submesh->GetParentElementIDMap()[i], *submesh->GetParent(), FET, T1, T2);
            MFEM_ASSERT(loc_bdr_attr.contains(attr) &&
                            loc_bdr_attr.at(attr).contains(nbr_attr),
                        "Missing libCEED boundary attribute for attribute " << attr << "!");
            return loc_bdr_attr.at(attr).at(nbr_attr);
          };
        }
      }
      // Non-submesh mesh, or submesh with self-rebuilt CEED data: use loc_attr directly.
      return [&](int i)
      {
        const int attr = mesh.GetAttribute(i);
        MFEM_ASSERT(loc_attr.contains(attr),
                    "Missing libCEED domain attribute for attribute " << attr << "!");
        return loc_attr.at(attr);
      };
    }();
    for (auto &[geom, indices] : element_indices)
    {
      Vector elem_attr(indices.size());
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = GetCeedAttribute(indices[k]);
      }
      geom_data_map.emplace(
          geom, AssembleGeometryData(ceed, geom, indices, *mesh.GetNodes(), elem_attr));
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
    auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
    auto GetCeedAttribute = [&](int i)
    {
      const int attr = mesh.GetBdrAttribute(i);
      const int nbr_attr = GetBdrNeighborAttribute(i, mesh, FET, T1, T2);
      MFEM_ASSERT(loc_bdr_attr.contains(attr) && loc_bdr_attr.at(attr).contains(nbr_attr),
                  "Missing libCEED boundary attribute for attribute " << attr << "!");
      return loc_bdr_attr.at(attr).at(nbr_attr);
    };
    for (auto &[geom, indices] : element_indices)
    {
      Vector elem_attr(indices.size());
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = GetCeedAttribute(indices[k]);
      }
      geom_data_map.emplace(
          geom, AssembleGeometryData(ceed, geom, indices, *mesh.GetNodes(), elem_attr));
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
    geom_data_map =
        BuildCeedGeomFactorData(*mesh, loc_attr, loc_bdr_attr, ceed, ceed_from_self);
  }
  return geom_data_map;
}

void Mesh::ResetCeedObjects()
{
  for (auto &[ceed, geom_data_map] : geom_data)
  {
    for (auto &[key, val] : geom_data_map)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&val.geom_data));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val.geom_data_restr));
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
  ResetCeedObjects();

  // Complete the transfer of interior boundary sides through a rebalance.
  if (crack_gf)
  {
    crack_fespace->Update();
    crack_gf->Update();
    const int ne = mesh->GetNE();
    const auto *h_gf = crack_gf->HostRead();
    for (std::size_t k = 0; k < crack_gf_keys.size(); k++)
    {
      auto &entry = crack_sides[crack_gf_keys[k]];
      if (entry.mesh != mesh.get() || entry.sequence + 1 != mesh->GetSequence() ||
          mesh->GetLastOperation() != mfem::Mesh::REBALANCE)
      {
        continue;
      }
      const auto *h = h_gf + 8 * k * ne;
      auto &sides = entry.sides;
      sides.Reset(ne);
      for (int e = 0; e < ne; e++)
      {
        sides.split[e] = static_cast<std::uint32_t>(std::llround(h[e]));
        for (int w = 0; w < 4; w++)
        {
          sides.version[e][w] =
              static_cast<std::uint32_t>(std::llround(h[(1 + w) * ne + e]));
        }
        sides.carrier[e] = static_cast<std::uint64_t>(std::llround(h[5 * ne + e])) |
                           (static_cast<std::uint64_t>(std::llround(h[6 * ne + e])) << 32);
        sides.side[e] = std::llround(h[7 * ne + e]);
      }
      entry.sequence = mesh->GetSequence();
    }
    crack_gf.reset();
    crack_fespace.reset();
    crack_fec.reset();
    crack_gf_keys.clear();
  }
}

const CrackSides &Mesh::GetCrackSides(const std::vector<int> &attr_list) const
{
  std::vector<int> key(attr_list);
  std::sort(key.begin(), key.end());
  key.erase(std::unique(key.begin(), key.end()), key.end());
  auto &entry = crack_sides[key];
  if (entry.mesh == mesh.get() && entry.sequence == mesh->GetSequence())
  {
    return entry.sides;
  }
  const int bdr_attr_max = mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 0;
  for (auto attr : key)
  {
    if (attr > 0 && attr <= bdr_attr_max)
    {
      marker[attr - 1] = 1;
    }
  }
  if (!mesh::HasHangingEntities(*mesh))
  {
    entry.sides = mesh::ComputeCrackSides(*mesh, marker);
  }
  else if (!mesh::ReconstructCrackSides(*mesh, marker, entry.sides))
  {
    // The sides of a nonconforming mesh refined in this simulation are inherited through
    // the refinements, so this is a nonconforming input mesh whose refinement hierarchy
    // cannot be derefined (Palace itself only refines isotropically, and MFEM does not
    // refine meshes with pyramids and nodes).
    Mpi::Warning(mesh->GetComm(),
                 "The sides of the interior boundaries cannot be identified on this "
                 "nonconforming mesh, which has pyramids, anisotropic refinements, or a "
                 "nonconforming coarsest mesh: error estimation will assume continuous "
                 "fluxes across interior boundaries for the whole simulation, which "
                 "overestimates the error next to them!\n");
    entry.sides.Reset(mesh->GetNE());
  }
  entry.mesh = mesh.get();
  entry.sequence = mesh->GetSequence();
  return entry.sides;
}

void Mesh::InheritCrackSides(const std::vector<int> &attr_list,
                             const Mesh &coarse_mesh) const
{
  // The parents of the coarse-to-fine transformations are checked against the number of
  // coarse elements by mesh::InheritCrackSides.
  MFEM_VERIFY(mesh->GetLastOperation() == mfem::Mesh::REFINE,
              "Interior boundary sides can only be inherited from the mesh this mesh was "
              "refined from!");
  const auto &coarse_sides = coarse_mesh.GetCrackSides(attr_list);
  std::vector<int> key(attr_list);
  std::sort(key.begin(), key.end());
  key.erase(std::unique(key.begin(), key.end()), key.end());
  auto &entry = crack_sides[key];
  entry.sides =
      mesh::InheritCrackSides(*mesh, mesh->GetRefinementTransforms(), coarse_sides);
  entry.mesh = mesh.get();
  entry.sequence = mesh->GetSequence();
}

void Mesh::RefineCrackSides()
{
  // Only a single refinement since the sides were computed can be followed.
  const bool refined =
      mesh->GetLastOperation() == mfem::Mesh::REFINE && mesh->Nonconforming();
  std::vector<CrackSidesEntry *> entries;
  crack_gf_keys.clear();
  for (auto &[key, entry] : crack_sides)
  {
    if (!refined || entry.mesh != mesh.get() || mesh->GetSequence() != entry.sequence + 1)
    {
      entry.mesh = nullptr;
      continue;
    }
    entry.sides =
        mesh::InheritCrackSides(*mesh, mesh->GetRefinementTransforms(), entry.sides);
    entry.sequence = mesh->GetSequence();
    entries.push_back(&entry);
    crack_gf_keys.push_back(key);
  }
  if (entries.empty())
  {
    return;
  }

  // Prepare the transfer of the sides through a possible rebalance, which keeps the mesh
  // object for a nonconforming mesh: the split entities, the versions and carriers (32-bit
  // words) and the side label of each element, all exactly represented as doubles, for
  // each list of interior boundaries in use.
  const int nc = 8 * static_cast<int>(entries.size());
  crack_fec = std::make_unique<mfem::L2_FECollection>(0, mesh->Dimension());
  crack_fespace =
      std::make_unique<mfem::ParFiniteElementSpace>(mesh.get(), crack_fec.get(), nc);
  crack_gf = std::make_unique<mfem::ParGridFunction>(crack_fespace.get());
  const int ne = mesh->GetNE();
  auto *h_gf = crack_gf->HostWrite();
  for (std::size_t k = 0; k < entries.size(); k++)
  {
    const auto &sides = entries[k]->sides;
    auto *h = h_gf + 8 * k * ne;
    for (int e = 0; e < ne; e++)
    {
      h[e] = static_cast<double>(sides.split[e]);
      for (int w = 0; w < 4; w++)
      {
        h[(1 + w) * ne + e] = static_cast<double>(sides.version[e][w]);
      }
      h[5 * ne + e] = static_cast<double>(sides.carrier[e] & 0xFFFFFFFFu);
      h[6 * ne + e] = static_cast<double>(sides.carrier[e] >> 32);
      h[7 * ne + e] = static_cast<double>(sides.side[e]);
    }
  }
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
  ResetCeedObjects();
}

}  // namespace palace
