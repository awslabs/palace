// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

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

template <typename T>
auto AssembleGeometryData(const mfem::GridFunction &mesh_nodes, Ceed ceed,
                          mfem::Geometry::Type geom, std::vector<int> &indices,
                          T GetCeedAttribute)
{
  const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
  const mfem::Mesh &mesh = *mesh_fespace.GetMesh();

  ceed::CeedGeomFactorData data;
  data.dim = mfem::Geometry::Dimension[geom];
  data.space_dim = mesh.SpaceDimension();
  data.indices = std::move(indices);
  const std::size_t num_elem = data.indices.size();

  // Allocate storage for geometry factor data (stored as attribute + quadrature weight +
  // Jacobian, column-major).
  CeedElemRestriction mesh_restr =
      FiniteElementSpace::BuildCeedElemRestriction(mesh_fespace, ceed, geom, data.indices);
  CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(mesh_fespace, ceed, geom);
  CeedInt num_qpts, geom_data_size = 2 + data.space_dim * data.dim;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &num_qpts));
  PalaceCeedCall(
      ceed, CeedVectorCreate(ceed, num_elem * num_qpts * geom_data_size, &data.geom_data));

  // Data for quadrature point i, component j, element k is found at index i * strides[0] +
  // j * strides[1] + k * strides[2].
  CeedMemType mem;
  CeedInt strides[3];
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    // GPU backends have CEED_STRIDES_BACKEND = {1, num_elem * num_qpts, num_qpts}.
    strides[0] = 1;
    strides[1] = num_elem * num_qpts;
    strides[2] = num_qpts;
  }
  else
  {
    // CPU backends have CEED_STRIDES_BACKEND = {1, num_qpts, num_qpts * geom_data_size}.
    strides[0] = 1;
    strides[1] = num_qpts;
    strides[2] = num_qpts * geom_data_size;
  }
  PalaceCeedCall(ceed,
                 CeedElemRestrictionCreateStrided(ceed, num_elem, num_qpts, geom_data_size,
                                                  num_elem * num_qpts * geom_data_size,
                                                  strides, &data.geom_data_restr));

  // Compute the required geometry factors at quadrature points.
  CeedVector mesh_nodes_vec;
  ceed::InitCeedVector(mesh_nodes, ceed, &mesh_nodes_vec);

  ceed::AssembleCeedGeometryData(ceed, mesh_restr, mesh_basis, mesh_nodes_vec,
                                 data.geom_data, data.geom_data_restr);

  PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));

  // Compute element attribute quadrature data. All inputs to a QFunction require the same
  // number of quadrature points, so we store the attribute at each quadrature point. This
  // is the first component of the quadrature data.
  {
    CeedScalar *geom_data_array;
    PalaceCeedCall(ceed,
                   CeedVectorGetArray(data.geom_data, CEED_MEM_HOST, &geom_data_array));
    for (std::size_t k = 0; k < num_elem; k++)
    {
      const auto attr = GetCeedAttribute(data.indices[k]);
      for (CeedInt i = 0; i < num_qpts; i++)
      {
        geom_data_array[i * strides[0] + k * strides[2]] = attr;
      }
    }
    CeedVectorRestoreArray(data.geom_data, &geom_data_array);
  }

  return data;
}

auto BuildCeedGeomFactorData(
    const mfem::ParMesh &mesh, const std::unordered_map<int, int> &loc_attr,
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr, Ceed ceed)
{
  // Create a list of the element indices in the mesh corresponding to a given thread and
  // element geometry type and corresponding geometry factor data. libCEED operators will be
  // constructed in parallel over threads, where each thread builds a composite operator
  // with sub-operators for each geometry.
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
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
      if (const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(&mesh))
      {
        MFEM_VERIFY(submesh->GetFrom() == mfem::SubMesh::From::Boundary,
                    "Unexpected non-SubMesh object for BuildCeedGeomFactorData with Mesh "
                    "with (dim, space_dim) = ("
                        << mesh.Dimension() << ", " << mesh.SpaceDimension() << ")!");
        return [&, submesh](int i)
        {
          // Mesh is actually a boundary submesh, so we use the boundary attribute mappings
          // from the parent mesh.
          const int attr = mesh.GetAttribute(i);
          const int nbr_attr = GetBdrNeighborAttribute(submesh->GetParentElementIDMap()[i],
                                                       *submesh->GetParent(), FET, T1, T2);
          MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                          loc_bdr_attr.at(attr).find(nbr_attr) !=
                              loc_bdr_attr.at(attr).end(),
                      "Missing libCEED boundary attribute for attribute " << attr << "!");
          return loc_bdr_attr.at(attr).at(nbr_attr);
        };
      }
      else
      {
        return [&](int i)
        {
          const int attr = mesh.GetAttribute(i);
          MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                      "Missing libCEED domain attribute for attribute " << attr << "!");
          return loc_attr.at(attr);
        };
      }
    }();
    for (auto &[geom, indices] : element_indices)
    {
      geom_data_map.emplace(geom, AssembleGeometryData(*mesh.GetNodes(), ceed, geom,
                                                       indices, GetCeedAttribute));
    }
  }

  // Then boundary elements (no support for boundary integrators on meshes embedded in
  // higher dimensional space for now).
  if (mesh.Dimension() == mesh.SpaceDimension())
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
      MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                      loc_bdr_attr.at(attr).find(nbr_attr) != loc_bdr_attr.at(attr).end(),
                  "Missing libCEED boundary attribute for attribute " << attr << "!");
      return loc_bdr_attr.at(attr).at(nbr_attr);
    };
    for (auto &[geom, indices] : element_indices)
    {
      geom_data_map.emplace(geom, AssembleGeometryData(*mesh.GetNodes(), ceed, geom,
                                                       indices, GetCeedAttribute));
    }
  }

  return geom_data_map;
}

}  // namespace

const ceed::GeometryObjectMap<ceed::CeedGeomFactorData> &
Mesh::GetCeedGeomFactorData(Ceed ceed) const
{
  MFEM_VERIFY(!loc_attr.empty(),
              "Mesh attribute mappings have not been built for GetCeedGeomFactorData!");
  auto it = geom_data.find(ceed);
  MFEM_ASSERT(it != geom_data.end(), "Unknown Ceed context in GetCeedGeomFactorData!");
  auto &geom_data_map = it->second;
  if (geom_data_map.empty())
  {
    geom_data_map = BuildCeedGeomFactorData(*mesh, loc_attr, loc_bdr_attr, ceed);
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
  auto &parent_mesh = GetParentMesh(*mesh);
  parent_mesh.ExchangeFaceNbrData();
  loc_attr.clear();
  loc_bdr_attr.clear();
  loc_attr = BuildCeedAttributes(parent_mesh);
  loc_bdr_attr = BuildCeedBdrAttributes(parent_mesh);
  ResetCeedObjects();
}

}  // namespace palace
