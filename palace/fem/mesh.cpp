// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mesh.hpp"

#include "fem/libceed/integrator.hpp"
#include "utils/omp.hpp"

namespace palace
{

namespace ceed
{

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

CeedGeomFactorData CeedGeomFactorDataCreate(Ceed ceed)
{
  return std::make_unique<CeedGeomFactorData_private>(ceed);
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
  CeedInt nqpts;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
  ceed::AssembleCeedGeometryData(info, ceed, mesh_restr, mesh_basis, mesh_nodes, data);
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));

  // Compute element attribute quadrature data. This should ideally be a single scalar per
  // element but all fields associated with a CeedOperator require the same number of
  // quadrature points.
  {
    const auto ne = indices.size();
    const bool use_bdr = (data->dim != mesh.Dimension());
    data->attr.SetSize(ne * nqpts);
    for (std::size_t j = 0; j < ne; j++)
    {
      const int e = indices[j];
      int attr = use_bdr ? mesh.GetBdrAttribute(e) : mesh.GetAttribute(e);
      MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                  "Missing local domain or boundary attribute for attribute " << attr
                                                                              << "!");
      attr = loc_attr[attr] - 1;  // Convert to 0-based indexing for libCEED QFunctions
      for (CeedInt q = 0; q < nqpts; q++)
      {
        data->attr[j * nqpts + q] = attr;
      }
    }
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, 1, ne * nqpts,
                                                          CEED_STRIDES_BACKEND,
                                                          &data->attr_restr));
    ceed::InitCeedVector(data->attr, ceed, &data->attr_vec);
  }

  // Save mesh element indices.
  data->indices = std::move(indices);

  return data;
}

}  // namespace

std::unordered_map<int, int> &Mesh::BuildAttributesGlobalToLocal(bool use_bdr) const
{
  if (!use_bdr)
  {
    // Set up sparse map from global domain attributes to local ones on this process.
    // Include ghost elements for all shared faces so we have their material properties
    // stored locally.
    loc_attr.clear();
    int count = 0;
    for (int i = 0; i < Mesh().GetNE(); i++)
    {
      const int attr = Mesh().GetAttribute(i);
      auto it = loc_attr.find(attr);
      if (it == loc_attr.end())
      {
        loc_attr[attr] = ++count;
      }
    }
    for (int i = 0; i < Mesh().GetNSharedFaces(); i++)
    {
      const mfem::FaceElementTransformations &FET = *Mesh().GetSharedFaceTransformations(i);
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
    return loc_attr;
  }
  else
  {
    // Set up sparse map from global boundary attributes to local ones on this process.
    loc_bdr_attr.clear();
    int count = 0;
    for (int i = 0; i < Mesh().GetNBE(); i++)
    {
      const int attr = Mesh().GetBdrAttribute(i);
      auto it = loc_bdr_attr.find(attr);
      if (it == loc_bdr_attr.end())
      {
        loc_bdr_attr[attr] = ++count;
      }
    }
    return loc_bdr_attr;
  }
}

std::unordered_map<int, int> &Mesh::BuildLocalToSharedFaceMap() const
{
  // Construct shared face mapping for boundary coefficients. The inverse mapping is
  // constructed as part of mfem::ParMesh, but we need this mapping when looping over
  // all mesh faces.
  local_to_shared.clear();
  local_to_shared.reserve(Mesh().GetNSharedFaces());
  for (int i = 0; i < Mesh().GetNSharedFaces(); i++)
  {
    local_to_shared[Mesh().GetSharedFace(i)] = i;
  }
  return local_to_shared;
}

ceed::CeedObjectMap<ceed::CeedGeomFactorData> &Mesh::BuildCeedGeomFactorData() const
{
  const mfem::GridFunction &mesh_nodes = *Mesh().GetNodes();
  geom_data.clear();

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
      const int ne = Mesh().GetNE();
      const int stride = (ne + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, ne);
      constexpr bool use_bdr = false;
      auto element_indices = GetElementIndices(Mesh(), use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data =
            AssembleGeometryData(mesh_nodes, loc_attr, ceed, geom, indices);
        PalacePragmaOmp(critical(GeomFactorData))
        {
          geom_data.emplace(std::make_pair(ceed, geom), std::move(data));
        }
      }
    }

    // Then boundary elements.
    {
      const int nbe = Mesh().GetNBE();
      const int stride = (nbe + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, nbe);
      constexpr bool use_bdr = true;
      auto element_indices = GetElementIndices(Mesh(), use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data =
            AssembleGeometryData(mesh_nodes, loc_bdr_attr, ceed, geom, indices);
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
