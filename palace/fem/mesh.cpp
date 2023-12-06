// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

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

// Count the number of elements of each type in the local mesh.
std::unordered_map<mfem::Geometry::Type, std::vector<int>>
GetElementIndices(const mfem::ParMesh &mesh, bool use_bdr, int start, int stop)
{
  std::unordered_map<mfem::Geometry::Type, int> counts, offsets;
  std::unordered_map<mfem::Geometry::Type, std::vector<int>> element_indices;

  for (int i = start; i < stop; i++)
  {
    const auto key = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto it = counts.find(key);
    if (it == counts.end())
    {
      counts[key] = 1;
    }
    else
    {
      it->second++;
    }
  }

  // Populate the indices arrays for each element geometry.
  for (const auto &[key, val] : counts)
  {
    offsets[key] = 0;
    element_indices[key] = std::vector<int>(val);
  }
  for (int i = start; i < stop; i++)
  {
    const auto key = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto &offset = offsets[key];
    auto &indices = element_indices[key];
    indices[offset++] = i;
  }

  return element_indices;
}

ceed::CeedGeomFactorData AssembleGeometryData(const mfem::GridFunction &mesh_nodes,
                                              const std::unordered_map<int, int> &loc_attr,
                                              Ceed ceed, mfem::Geometry::Type geom,
                                              std::vector<int> &indices)
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

  // Compute element attribute quadrature data. All inputs to a QFunction require the same
  // number of quadrature points, so we create a basis to interpolate the single attribute
  // per element to all quadrature points.
  {
    const std::size_t ne = indices.size();
    const bool use_bdr = (data->dim != mesh.Dimension());
    CeedInt nqpts;
    PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
    data->attr.SetSize(ne * nqpts);
    // data->attr.SetSize(ne);
    for (std::size_t i = 0; i < ne; i++)
    {
      // Convert to 0-based indexing for libCEED QFunctions.
      const int e = indices[i];
      int attr = use_bdr ? mesh.GetBdrAttribute(e) : mesh.GetAttribute(e);
      MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                  "Missing local domain or boundary attribute for attribute " << attr
                                                                              << "!");
      attr = loc_attr.at(attr) - 0.5;
      for (CeedInt q = 0; q < nqpts; q++)
      {
        data->attr[i * nqpts + q] = attr;
      }
      // data->attr[i] = loc_attr.at(attr) - 0.5;
    }

    CeedInt strides[3] = {1, 1, nqpts};
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, 1, ne * nqpts,
                                                          strides, &data->attr_restr));

    // CeedInt strides[3] = {1, 1, 1};
    // PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, 1, 1, ne,
    //                                                       strides, &data->attr_restr));

    // CeedInt nqpts;
    // PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
    // mfem::DenseMatrix qX(data->dim, nqpts);
    // mfem::Vector qW(nqpts), B(nqpts), G(nqpts * data->dim);
    // B = 1.0;
    // G = 0.0;
    // qX = 0.0;
    // qW = 0.0;
    // PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, ceed::GetCeedTopology(geom),
    //                                        1, 1, nqpts, B.GetData(),
    //                                        G.GetData(), qX.GetData(), qW.GetData(),
    //                                        &data->attr_basis));

    ceed::InitCeedVector(data->attr, ceed, &data->attr_vec);
  }

  // Save mesh element indices and clean up.
  data->indices = std::move(indices);
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));

  return data;
}

}  // namespace

void Mesh::CheckSequence() const
{
  PalacePragmaOmp(critical(CheckSequence))
  {
    if (sequence != mesh->GetSequence())
    {
      ClearData();
      sequence = mesh->GetSequence();
    }
  }
}

std::unordered_map<int, int> &Mesh::BuildAttributesGlobalToLocal(bool use_bdr) const
{
  auto &mesh_obj = Get();
  const_cast<mfem::ParMesh &>(mesh_obj).ExchangeFaceNbrData();
  if (!use_bdr)
  {
    // Set up sparse map from global domain attributes to local ones on this process.
    // Include ghost elements for all shared faces so we have their material properties
    // stored locally.
    loc_attr.clear();
    int count = 0;
    for (int i = 0; i < mesh_obj.GetNE(); i++)
    {
      const int attr = mesh_obj.GetAttribute(i);
      auto it = loc_attr.find(attr);
      if (it == loc_attr.end())
      {
        loc_attr[attr] = ++count;
      }
    }
    for (int i = 0; i < mesh_obj.GetNSharedFaces(); i++)
    {
      mfem::FaceElementTransformations &FET =
          *const_cast<mfem::ParMesh &>(mesh_obj).GetSharedFaceTransformations(i);
      int attr = FET.GetElement1Transformation().Attribute;
      auto it = loc_attr.find(attr);
      if (it == loc_attr.end())
      {
        loc_attr[attr] = ++count;
      }
      attr = FET.GetElement2Transformation().Attribute;
      it = loc_attr.find(attr);
      if (it == loc_attr.end())
      {
        loc_attr[attr] = ++count;
      }
    }

    // //XX TODO TESTING
    // for (int i = 0; i < mesh_obj.attributes.Size(); i++)
    // {
    // 	loc_attr[mesh_obj.attributes[i]] = mesh_obj.attributes[i];
    // }

    return loc_attr;
  }
  else
  {
    // Set up sparse map from global boundary attributes to local ones on this process.
    loc_bdr_attr.clear();
    int count = 0;
    for (int i = 0; i < mesh_obj.GetNBE(); i++)
    {
      const int attr = mesh_obj.GetBdrAttribute(i);
      auto it = loc_bdr_attr.find(attr);
      if (it == loc_bdr_attr.end())
      {
        loc_bdr_attr[attr] = ++count;
      }
    }

    // //XX TODO TESTING
    // for (int i = 0; i < mesh_obj.bdr_attributes.Size(); i++)
    // {
    // 	loc_bdr_attr[mesh_obj.bdr_attributes[i]] = mesh_obj.bdr_attributes[i];
    // }

    return loc_bdr_attr;
  }
}

std::unordered_map<int, int> &Mesh::BuildLocalToSharedFaceMap() const
{
  // Construct shared face mapping for boundary coefficients. The inverse mapping is
  // constructed as part of mfem::ParMesh, but we need this mapping when looping over
  // all mesh faces.
  const auto &mesh_obj = Get();
  local_to_shared.clear();
  local_to_shared.reserve(mesh_obj.GetNSharedFaces());
  for (int i = 0; i < mesh_obj.GetNSharedFaces(); i++)
  {
    local_to_shared[mesh_obj.GetSharedFace(i)] = i;
  }
  return local_to_shared;
}

ceed::CeedObjectMap<ceed::CeedGeomFactorData> &Mesh::BuildCeedGeomFactorData() const
{
  // Initialize shared objects outside of the OpenMP region.
  const auto &mesh_obj = Get();
  const mfem::GridFunction &mesh_nodes = *mesh_obj.GetNodes();
  geom_data.clear();
  GetAttributeGlobalToLocal();
  GetBdrAttributeGlobalToLocal();

  // Create a list of the element indices in the mesh corresponding to a given thread and
  // element geometry type. libCEED operators will be constructed in parallel over threads,
  // where each thread builds a composite operator with sub-operators for each geometry.
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];

    // First domain elements.
    {
      const int ne = mesh_obj.GetNE();
      const int stride = (ne + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, ne);
      constexpr bool use_bdr = false;
      auto element_indices = GetElementIndices(mesh_obj, use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data = AssembleGeometryData(
            mesh_nodes, GetAttributeGlobalToLocal(), ceed, geom, indices);
        PalacePragmaOmp(critical(GeomFactorData))
        {
          geom_data.emplace(std::make_pair(ceed, geom), std::move(data));
        }
      }
    }

    // Then boundary elements.
    {
      const int nbe = mesh_obj.GetNBE();
      const int stride = (nbe + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, nbe);
      constexpr bool use_bdr = true;
      auto element_indices = GetElementIndices(mesh_obj, use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data = AssembleGeometryData(
            mesh_nodes, GetBdrAttributeGlobalToLocal(), ceed, geom, indices);
        PalacePragmaOmp(critical(GeomFactorData))
        {
          geom_data.emplace(std::make_pair(ceed, geom), std::move(data));
        }
      }
    }
  }

  return geom_data;
}

}  // namespace palace
