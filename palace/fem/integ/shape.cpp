// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/shape_qf.h"

namespace palace
{

using namespace ceed;

void DiffusionShapeIntegrator::Assemble(Ceed ceed, CeedElemRestriction restr,
                                        CeedBasis basis, CeedVector u, CeedVector w,
                                        CeedElemRestriction mesh_restr,
                                        CeedBasis mesh_basis, CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op) const
{
  CeedQFunctionInfo info;
  CeedInt dim, space_dim, num_comp;
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(basis, &num_comp));
  MFEM_VERIFY(num_comp == 1,
              "DiffusionShapeIntegrator requires a space with a single component!");
  switch (10 * space_dim + dim)
  {
    case 22:
      info.apply_qf = f_shape_diffusion_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_diffusion_22_loc);
      break;
    case 33:
      info.apply_qf = f_shape_diffusion_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_diffusion_33_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dim << ", " << space_dim << ") for DiffusionShapeIntegrator!");
  }
  info.trial_ops = EvalMode::Grad;
  info.test_ops = EvalMode::Grad;

  auto ctx = PopulateCoefficientContext(space_dim, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

void MassShapeIntegrator::Assemble(Ceed ceed, CeedElemRestriction restr, CeedBasis basis,
                                   CeedVector u, CeedVector w,
                                   CeedElemRestriction mesh_restr, CeedBasis mesh_basis,
                                   CeedVector geom_data,
                                   CeedElemRestriction geom_data_restr,
                                   CeedOperator *op) const
{
  CeedQFunctionInfo info;
  CeedInt dim, space_dim, num_comp;
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(basis, &num_comp));
  MFEM_VERIFY(num_comp == 1,
              "MassShapeIntegrator requires a space with a single component!");
  switch (10 * space_dim + dim)
  {
    case 21:
      info.apply_qf = f_shape_mass_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_mass_21_loc);
      break;
    case 22:
      info.apply_qf = f_shape_mass_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_mass_22_loc);
      break;
    case 32:
      info.apply_qf = f_shape_mass_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_mass_32_loc);
      break;
    case 33:
      info.apply_qf = f_shape_mass_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_mass_33_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = (" << dim << ", " << space_dim
                                                         << ") for MassShapeIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Interp;

  auto ctx = PopulateCoefficientContext(1, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

}  // namespace palace
