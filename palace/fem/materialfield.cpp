// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialfield.hpp"

#include "fem/integrator.hpp"
#include "fem/mesh.hpp"

namespace palace
{

MaterialField::MaterialField(const mfem::ParGridFunction &gf, double scaling)
  : gf(&gf), scaling(scaling), local_max(scaling * gf.Max()),
    local_min(scaling * gf.Min())
{
}

MaterialField::~MaterialField()
{
  for (auto &[ceed, geom_map] : data)
  {
    for (auto &[geom, d] : geom_map)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&d.data));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&d.data_restr));
    }
  }
}

CeedMaterialFieldData MaterialField::BuildData(Ceed ceed, mfem::Geometry::Type geom,
                                               const std::vector<int> &indices) const
{
  const mfem::Mesh &mesh = *gf->ParFESpace()->GetMesh();
  const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);
  const int num_qpts = ir.GetNPoints();
  const std::size_t num_elem = indices.size();

  // Evaluate the grid function at the reference quadrature points of each element. The
  // grid function is defined on the same mesh, so its element numbering matches the mesh.
  std::vector<double> values(num_elem * num_qpts);
  for (std::size_t e = 0; e < num_elem; e++)
  {
    for (int q = 0; q < num_qpts; q++)
    {
      values[e * num_qpts + q] = scaling * gf->GetValue(indices[e], ir.IntPoint(q));
    }
  }

  CeedMaterialFieldData d;
  const CeedSize size = static_cast<CeedSize>(values.size());
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, size, &d.data));
  PalaceCeedCall(
      ceed, CeedVectorSetArray(d.data, CEED_MEM_HOST, CEED_COPY_VALUES, values.data()));
  PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                           ceed, static_cast<CeedInt>(num_elem), num_qpts, 1, size,
                           CEED_STRIDES_BACKEND, &d.data_restr));
  return d;
}

const CeedMaterialFieldData &MaterialField::GetData(Ceed ceed, mfem::Geometry::Type geom,
                                                    const std::vector<int> &indices) const
{
  auto &geom_map = data[ceed];
  auto it = geom_map.find(geom);
  if (it == geom_map.end())
  {
    it = geom_map.emplace(geom, BuildData(ceed, geom, indices)).first;
  }
  return it->second;
}

void MaterialField::WarmCache(Ceed ceed, const Mesh &mesh) const
{
  // A material field is a domain quantity: only warm the domain element geometries (the
  // boundary geometries in the same map have a lower dimension and different element
  // numbering).
  for (const auto &[geom, geom_data] : mesh.GetCeedGeomFactorData(ceed))
  {
    if (mfem::Geometry::Dimension[geom] == mesh.Dimension())
    {
      GetData(ceed, geom, geom_data.indices);
    }
  }
}

}  // namespace palace
