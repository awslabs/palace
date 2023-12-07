// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/libceed/integrator.hpp"
#include "utils/omp.hpp"

namespace palace
{

namespace ceed
{

namespace
{

CeedGeomFactorData CeedGeomFactorDataCreate(Ceed ceed)
{
  return std::make_unique<CeedGeomFactorData_private>(ceed);
}

}  // namespace

CeedGeomFactorData_private::~CeedGeomFactorData_private()
{
  PalaceCeedCall(ceed, CeedVectorDestroy(&wdetJ_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&adjJt_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&J_vec));
  PalaceCeedCall(ceed, CeedVectorDestroy(&attr_vec));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&wdetJ_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&adjJt_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&J_restr));
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&attr_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&attr_basis));
}

}  // namespace ceed

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

auto BuildLocalToSharedFaceMap(const mfem::ParMesh &mesh)
{
  // Construct shared face mapping for boundary coefficients. The inverse mapping is
  // constructed as part of mfem::ParMesh, but we need this mapping when looping over
  // all mesh faces.
  std::unordered_map<int, int> l2s;
  l2s.reserve(mesh.GetNSharedFaces());
  for (int i = 0; i < mesh.GetNSharedFaces(); i++)
  {
    l2s[mesh.GetSharedFace(i)] = i;
  }
  return l2s;
}

auto BuildAttributeGlobalToLocal(mfem::ParMesh &mesh)
{
  // Set up sparse map from global domain attributes to local ones on this process.
  // Include ghost elements for all shared faces so we have their material properties
  // stored locally.
  std::unordered_map<int, int> loc_attr;
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
    mfem::FaceElementTransformations &FET = *mesh.GetSharedFaceTransformations(i);
    int attr = FET.GetElement1Transformation().Attribute;
    if (loc_attr.find(attr) == loc_attr.end())
    {
      loc_attr[attr] = ++count;
    }
    attr = FET.GetElement2Transformation().Attribute;
    if (loc_attr.find(attr) == loc_attr.end())
    {
      loc_attr[attr] = ++count;
    }
  }
  return loc_attr;
}

auto GetBdrNeighborAttribute(int i, mfem::ParMesh &mesh,
                             const std::unordered_map<int, int> &face_loc_to_shared)
{

  // XX TODO CONST FOR MESH INPUT POSSIBLE?? GetBdrElementTrasformation IS NOT THREAD
  // SAFE...

  // For internal boundaries, use the element which corresponds to the vacuum domain, or
  // at least the one with the higher speed of light.
  mfem::ElementTransformation *T1, *T2;
  BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
      i, mesh, face_loc_to_shared, T1, T2);
  // return (T2 && GetLightSpeedMin(T2->Attribute) > GetLightSpeedMax(T1->Attribute))
  //           ? T2->Attribute
  //           : T1->Attribute;
  return (T2 && T2->Attribute < T1->Attribute) ? T2->Attribute : T1->Attribute;
}

auto BuildBdrAttributeGlobalToLocal(mfem::ParMesh &mesh,
                                    const std::unordered_map<int, int> &face_loc_to_shared)
{
  // Set up sparse map from global boundary attributes to local ones on this process. Each
  // original global boundary attribute maps to a key-value pairing of global domain
  // attributes which neighbor the given boundary and local boundary attributes.
  std::unordered_map<int, std::unordered_map<int, int>> loc_bdr_attr;
  int count = 0;
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    const int attr = mesh.GetBdrAttribute(i);
    const int nbr_attr = GetBdrNeighborAttribute(i, mesh, face_loc_to_shared);
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
    element_indices[it->first] = std::vector<int>(it->second);
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

  // XX TODO: In practice we do not need to compute all of the geometry data for every
  //          simulation type.
  auto data = ceed::CeedGeomFactorDataCreate(ceed);
  data->dim = mfem::Geometry::Dimension[geom];
  data->space_dim = mesh.SpaceDimension();

  // Compute the required geometry factors at quadrature points.
  constexpr auto info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate |
                        ceed::GeomFactorInfo::Jacobian;
  CeedElemRestriction mesh_restr =
      FiniteElementSpace::BuildCeedElemRestriction(mesh_fespace, ceed, geom, indices);
  CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(mesh_fespace, ceed, geom);
  ceed::AssembleCeedGeometryData(info, ceed, mesh_restr, mesh_basis, mesh_nodes, data);

  // XX TODO CLEANUP

  // Compute element attribute quadrature data. All inputs to a QFunction require the same
  // number of quadrature points, so we create a basis to interpolate the single attribute
  // per element to all quadrature points.
  {
    const std::size_t ne = indices.size();
    // CeedInt nqpts;
    // PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
    // data->attr.SetSize(ne * nqpts);
    data->attr.SetSize(ne);
    for (std::size_t i = 0; i < ne; i++)
    {
      // Convert to 0-based indexing for libCEED QFunctions.
      const auto attr = GetCeedAttribute(indices[i]);
      // for (CeedInt q = 0; q < nqpts; q++)
      // {
      //   data->attr[i * nqpts + q] = attr - 0.5;
      // }
      data->attr[i] = attr - 0.5;
    }

    // CeedInt strides[3] = {1, 1, nqpts};
    // PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, 1, ne * nqpts,
    //                                                       strides, &data->attr_restr));

    CeedInt strides[3] = {1, 1, 1};
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, 1, 1, ne, strides,
                                                          &data->attr_restr));

    CeedInt nqpts;
    PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
    mfem::DenseMatrix qX(data->dim, nqpts);
    mfem::Vector qW(nqpts), B(nqpts), G(nqpts * data->dim);
    B = 1.0;
    G = 0.0;
    qX = 0.0;
    qW = 0.0;
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, ceed::GetCeedTopology(geom), 1, 1, nqpts,
                                           B.GetData(), G.GetData(), qX.GetData(),
                                           qW.GetData(), &data->attr_basis));

    ceed::InitCeedVector(data->attr, ceed, &data->attr_vec);
  }

  // Save mesh element indices and clean up.
  data->indices = std::move(indices);
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));

  return data;
}

auto BuildCeedGeomFactorData(
    const mfem::ParMesh &mesh, const std::unordered_map<int, int> &face_loc_to_shared,
    const std::unordered_map<int, int> &loc_attr,
    const std::unordered_map<int, std::unordered_map<int, int>> &loc_bdr_attr, Ceed ceed)
{
  // Create a list of the element indices in the mesh corresponding to a given thread and
  // element geometry type and corresponding geometry factor data. libCEED operators will be
  // constructed in parallel over threads, where each thread builds a composite operator
  // with sub-operators for each geometry.
  std::size_t i;
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  for (i = 0; i < nt; i++)
  {
    if (ceed == ceed::internal::GetCeedObjects()[i])
    {
      break;
    }
  }
  MFEM_VERIFY(i < nt, "Unable to find matching Ceed context in BuildCeedGeomFactorData!");
  ceed::CeedGeomObjectMap<ceed::CeedGeomFactorData> geom_data;

  // First domain elements.
  {
    const int ne = mesh.GetNE();
    const int stride = (ne + nt - 1) / nt;
    const int start = i * stride;
    const int stop = std::min(start + stride, ne);
    constexpr bool use_bdr = false;
    auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
    std::function<int(int)> GetCeedAttribute;
    if (mesh.Dimension() == mesh.SpaceDimension())
    {
      GetCeedAttribute = [&](int i)
      {
        const int attr = mesh.GetAttribute(i);
        MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                    "Missing local domain attribute for attribute " << attr << "!");
        return attr;
      };
    }
    else
    {
      const auto *submesh = dynamic_cast<const mfem::ParSubMesh *>(&mesh);
      MFEM_VERIFY(submesh && submesh->GetFrom() == mfem::SubMesh::From::Boundary,
                  "Unexpected non-SubMesh object for BuildCeedGeomFactorData with Mesh "
                  "with (dim, space_dim) = ("
                      << mesh.Dimension() << ", " << mesh.SpaceDimension() << ")!");
      GetCeedAttribute = [&](int i)
      {
        // Mesh is actually a boundary submesh, so we use the boundary attribute mappings
        // from the parent mesh.
        const int attr = mesh.GetAttribute(i);
        const int nbr_attr = GetBdrNeighborAttribute(
            submesh->GetParentElementIDMap()[i],
            const_cast<mfem::ParMesh &>(*submesh->GetParent()), face_loc_to_shared);
        MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                        loc_bdr_attr.at(attr).find(nbr_attr) != loc_bdr_attr.at(attr).end(),
                    "Missing local boundary attribute for attribute " << attr << "!");
        return loc_bdr_attr.at(attr).at(nbr_attr);
      };
    }
    for (auto &[geom, indices] : element_indices)
    {
      ceed::CeedGeomFactorData data =
          AssembleGeometryData(*mesh.GetNodes(), ceed, geom, indices, GetCeedAttribute);
      geom_data.emplace(geom, std::move(data));
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
      const int nbr_attr =
          GetBdrNeighborAttribute(i, const_cast<mfem::ParMesh &>(mesh), face_loc_to_shared);
      MFEM_ASSERT(loc_bdr_attr.find(attr) != loc_bdr_attr.end() &&
                      loc_bdr_attr.at(attr).find(nbr_attr) != loc_bdr_attr.at(attr).end(),
                  "Missing local boundary attribute for attribute " << attr << "!");
      return loc_bdr_attr.at(attr).at(nbr_attr);
    };
    for (auto &[geom, indices] : element_indices)
    {
      ceed::CeedGeomFactorData data =
          AssembleGeometryData(*mesh.GetNodes(), ceed, geom, indices, GetCeedAttribute);
      geom_data.emplace(geom, std::move(data));
    }
  }

  return geom_data;
}

}  // namespace

void Mesh::Rebuild() const
{
  // Attribute mappings, etc. are always constructed for the parent mesh (use boundary
  // attribute maps for the domain attributes of a boundary submesh, for example).
  auto &parent_mesh = GetParentMesh(*mesh);
  parent_mesh.ExchangeFaceNbrData();
  face_loc_to_shared.clear();
  loc_attr.clear();
  loc_bdr_attr.clear();
  face_loc_to_shared = BuildLocalToSharedFaceMap(parent_mesh);
  loc_attr = BuildAttributeGlobalToLocal(parent_mesh);
  loc_bdr_attr = BuildBdrAttributeGlobalToLocal(parent_mesh, face_loc_to_shared);
  DestroyCeedGeomFactorData();
}

const ceed::CeedGeomObjectMap<ceed::CeedGeomFactorData> &
Mesh::GetCeedGeomFactorData(Ceed ceed) const
{
  // No two threads should ever be calling this simultaneously with the same Ceed context.
  auto it = geom_data.find(ceed);
  if (it == geom_data.end())
  {
    auto val =
        BuildCeedGeomFactorData(*mesh, face_loc_to_shared, loc_attr, loc_bdr_attr, ceed);
    PalacePragmaOmp(critical(InitCeedGeomFactorData))
    {
      it = geom_data.emplace(ceed, std::move(val)).first;
    }
  }
  return it->second;
}

void Mesh::DestroyCeedGeomFactorData() const
{
  geom_data.clear();
}

}  // namespace palace
