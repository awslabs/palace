// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialfield.hpp"

#include <array>
#include <map>
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
  for (auto &[key, d] : data)
  {
    Ceed ceed = std::get<0>(key);
    PalaceCeedCall(ceed, CeedVectorDestroy(&d.data));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&d.data_restr));
  }
}

const mfem::ParGridFunction &MaterialField::GetLevelGridFunction(const Mesh &mesh) const
{
  if (gf->ParFESpace()->GetMesh() == &mesh.Get())
  {
    return *gf;
  }

  auto it = level_fields.find(&mesh);
  if (it == level_fields.end())
  {
    // Build a P1 grid function on the coarser mesh by copying the base grid function's
    // nodal values to matching vertices (the coarse mesh of a nested uniform refinement is
    // a subset of the fine mesh).
    LevelGridFunction lf;
    lf.fec = std::make_unique<mfem::H1_FECollection>(1, mesh.Dimension());
    auto &pmesh = const_cast<mfem::ParMesh &>(mesh.Get());
    lf.fes = std::make_unique<mfem::ParFiniteElementSpace>(&pmesh, lf.fec.get());
    lf.gf = std::make_unique<mfem::ParGridFunction>(lf.fes.get());

    const mfem::Mesh &base_mesh = *gf->ParFESpace()->GetMesh();
    const int base_sdim = base_mesh.SpaceDimension();
    std::map<std::array<double, 3>, double> lookup;
    mfem::Array<int> base_dofs;
    for (int v = 0; v < base_mesh.GetNV(); v++)
    {
      const double *x = base_mesh.GetVertex(v);
      gf->FESpace()->GetVertexDofs(v, base_dofs);
      lookup[{x[0], (base_sdim > 1) ? x[1] : 0.0, (base_sdim > 2) ? x[2] : 0.0}] =
          (*gf)(base_dofs[0]);
    }

    const int sdim = mesh.SpaceDimension();
    mfem::Array<int> dofs;
    for (int v = 0; v < mesh.Get().GetNV(); v++)
    {
      const double *x = mesh.Get().GetVertex(v);
      auto lit = lookup.find({x[0], (sdim > 1) ? x[1] : 0.0, (sdim > 2) ? x[2] : 0.0});
      MFEM_VERIFY(lit != lookup.end(),
                  "Coarse mesh vertex not found in the material field grid function!");
      lf.fes->GetVertexDofs(v, dofs);
      (*lf.gf)(dofs[0]) = lit->second;
    }
    lf.gf->ExchangeFaceNbrData();
    it = level_fields.emplace(&mesh, std::move(lf)).first;
  }
  return *it->second.gf;
}

CeedMaterialFieldData MaterialField::BuildData(Ceed ceed, const Mesh &mesh,
                                               mfem::Geometry::Type geom,
                                               const std::vector<int> &indices) const
{
  const mfem::ParGridFunction &field_gf = GetLevelGridFunction(mesh);
  const mfem::Mesh &mfem_mesh = *field_gf.ParFESpace()->GetMesh();
  const int q_order = fem::DefaultIntegrationOrder::Get(mfem_mesh, geom);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);
  const int num_qpts = ir.GetNPoints();
  const std::size_t num_elem = indices.size();

  // Evaluate the grid function at the reference quadrature points of each element.
  std::vector<double> values(num_elem * num_qpts);
  for (std::size_t e = 0; e < num_elem; e++)
  {
    for (int q = 0; q < num_qpts; q++)
    {
      values[e * num_qpts + q] = scaling * field_gf.GetValue(indices[e], ir.IntPoint(q));
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

const CeedMaterialFieldData &MaterialField::GetData(Ceed ceed, const Mesh &mesh,
                                                    mfem::Geometry::Type geom,
                                                    const std::vector<int> &indices) const
{
  const auto key = std::make_tuple(ceed, &mesh, geom);
  auto it = data.find(key);
  if (it == data.end())
  {
    it = data.emplace(key, BuildData(ceed, mesh, geom, indices)).first;
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
      GetData(ceed, mesh, geom, geom_data.indices);
    }
  }
}

}  // namespace palace
