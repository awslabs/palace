// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include <algorithm>
#include "fem/multigrid.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

BoundaryModeOperator::BoundaryModeOperator(const IoData &iodata_,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh,
                                           const SubmeshFrame *frame_in)
  : iodata(iodata_), solver_order(iodata_.solver.order),
    solve_mesh(mesh.back().get()), frame(frame_in)
{
  MFEM_VERIFY(solve_mesh->Dimension() == 2,
              "BoundaryMode solver requires a 2D mesh (waveguide cross-section). When "
              "configured with \"Attributes\", the driver's PreprocessMesh hook extracts "
              "the cross-section before this operator is constructed.");

  // Phase 1: Material operator. If the mesh came from a 3D parent, rotate tensors from
  // the parent 3D frame onto the submesh tangent plane.
  mat_op = std::make_unique<MaterialOperator>(iodata, *solve_mesh);
  if (frame)
  {
    mat_op->RotateMaterialTensors(iodata, frame->e1, frame->e2, frame->normal);
  }

  // Phase 2: FE spaces.
  SetUpFESpaces(mesh);

  // Phase 3: Boundary operators.
  auto &pmesh = solve_mesh->Get();
  surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, *mat_op, pmesh);
  farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, *mat_op, pmesh);
  surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(iodata, *mat_op, pmesh);

  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
             GetNDSpace().GlobalTrueVSize() + GetH1Space().GlobalTrueVSize());
}

void BoundaryModeOperator::SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  const auto &mg = iodata.solver.linear;
  const int dim = solve_mesh->Dimension();

  // Collect Dirichlet boundary attributes.
  {
    const auto &pmesh = solve_mesh->Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    for (auto attr : iodata.boundaries.auxpec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    // When the mesh was extracted from a 3D parent, waveport attributes that landed on
    // the 2D boundary (i.e. other waveports touching this cross-section) behave as PEC
    // for the mode problem. For direct 2D input there is no parent, so no such inherited
    // attributes exist.
    if (frame)
    {
      for (const auto &[idx, data] : iodata.boundaries.waveport)
      {
        const auto &bm_attrs = iodata.solver.boundary_mode.attributes;
        for (auto attr : data.attributes)
        {
          if (std::find(bm_attrs.begin(), bm_attrs.end(), attr) != bm_attrs.end())
          {
            continue;
          }
          if (attr > 0 && attr <= bdr_attr_max)
          {
            dbc_bcs.Append(attr);
          }
        }
      }
    }
    dbc_bcs.Sort();
    dbc_bcs.Unique();
  }

  // FE collections.
  nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_aux_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  // RT collection (estimator flux recovery only; depth follows estimator_mg, mirroring
  // SpaceOperator::rt_fespaces). Order is solver_order - 1 because the recovered flux
  // lives one polynomial degree below the primary ND space.
  const int rt_mg_max_levels = mg.estimator_mg ? mg.mg_max_levels : 1;
  rt_fecs = fem::ConstructFECollections<mfem::RT_FECollection>(
      solver_order - 1, dim, rt_mg_max_levels, mg.mg_coarsening, false);

  nd_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
      mg.mg_max_levels, mesh, nd_fecs, &dbc_bcs, &nd_dbc_tdof_lists);
  h1_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, mesh, h1_fecs, &dbc_bcs, &h1_dbc_tdof_lists);
  h1_aux_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, mesh, h1_aux_fecs, &dbc_bcs, &h1_aux_dbc_tdof_lists);
  rt_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
      rt_mg_max_levels, mesh, rt_fecs);

  // L2 curl space for 2D B-field.
  l2_curl_fec = std::make_unique<mfem::L2_FECollection>(
      solver_order - 1, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL);
  l2_curl_fespace = std::make_unique<FiniteElementSpace>(*solve_mesh, l2_curl_fec.get());
}

}  // namespace palace
