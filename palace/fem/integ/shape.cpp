// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/integrator.hpp"

#include <ceed/backend.h>
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"

#include "fem/qfunctions/shape_qf.h"

namespace palace
{

using namespace ceed;

namespace
{

// Return the dimensions (10 * space_dim + dim) of the elements of a basis, after checking
// that it is a scalar basis of the given space.
CeedInt GetDimensions(Ceed ceed, CeedBasis basis, CeedElemRestriction geom_data_restr,
                      CeedFESpace fe_space, const char *name)
{
  CeedInt dim, space_dim, num_comp;
  CeedFESpace basis_fe_space;
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &dim));
  PalaceCeedCall(ceed, CeedGeometryDataGetSpaceDimension(geom_data_restr, dim, &space_dim));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(basis, &num_comp));
  PalaceCeedCall(ceed, CeedBasisGetFESpace(basis, &basis_fe_space));
  MFEM_VERIFY(num_comp == 1 && basis_fe_space == fe_space,
              name << " requires an " << ((fe_space == CEED_FE_SPACE_H1) ? "H1" : "H(curl)")
                   << " space with a single component!");
  return 10 * space_dim + dim;
}

}  // namespace

void DiffusionShapeIntegrator::Assemble(Ceed ceed, CeedElemRestriction restr,
                                        CeedBasis basis, CeedVector u, CeedVector w,
                                        CeedElemRestriction mesh_restr,
                                        CeedBasis mesh_basis, CeedVector geom_data,
                                        CeedElemRestriction geom_data_restr,
                                        CeedOperator *op) const
{
  CeedQFunctionInfo info;
  const CeedInt dims = GetDimensions(ceed, basis, geom_data_restr, CEED_FE_SPACE_H1,
                                     "DiffusionShapeIntegrator");
  switch (dims)
  {
    case 22:
      info.apply_qf = f_shape_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_22_loc);
      break;
    case 33:
      info.apply_qf = f_shape_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_33_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dims % 10 << ", " << dims / 10 << ") for DiffusionShapeIntegrator!");
  }
  info.trial_ops = EvalMode::Grad;
  info.test_ops = EvalMode::Grad;

  auto ctx = PopulateCoefficientContext(dims / 10, Q, transpose);
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
  const CeedInt dims =
      GetDimensions(ceed, basis, geom_data_restr, CEED_FE_SPACE_H1, "MassShapeIntegrator");
  switch (dims)
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
      MFEM_ABORT("Invalid value of (dim, space_dim) = (" << dims % 10 << ", " << dims / 10
                                                         << ") for MassShapeIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Interp;

  auto ctx = PopulateCoefficientContext(1, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

void VectorFEMassShapeIntegrator::Assemble(Ceed ceed, CeedElemRestriction restr,
                                           CeedBasis basis, CeedVector u, CeedVector w,
                                           CeedElemRestriction mesh_restr,
                                           CeedBasis mesh_basis, CeedVector geom_data,
                                           CeedElemRestriction geom_data_restr,
                                           CeedOperator *op) const
{
  CeedQFunctionInfo info;
  const CeedInt dims = GetDimensions(ceed, basis, geom_data_restr, CEED_FE_SPACE_HCURL,
                                     "VectorFEMassShapeIntegrator");
  switch (dims)
  {
    case 21:
      info.apply_qf = f_shape_hcurl_21;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_21_loc);
      break;
    case 22:
      info.apply_qf = f_shape_hcurl_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_22_loc);
      break;
    case 32:
      info.apply_qf = f_shape_hcurl_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_32_loc);
      break;
    case 33:
      info.apply_qf = f_shape_hcurl_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_33_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dims % 10 << ", " << dims / 10 << ") for VectorFEMassShapeIntegrator!");
  }
  info.trial_ops = EvalMode::Interp;
  info.test_ops = EvalMode::Interp;

  auto ctx = PopulateCoefficientContext(dims / 10, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

void CurlCurlShapeIntegrator::Assemble(Ceed ceed, CeedElemRestriction restr,
                                       CeedBasis basis, CeedVector u, CeedVector w,
                                       CeedElemRestriction mesh_restr, CeedBasis mesh_basis,
                                       CeedVector geom_data,
                                       CeedElemRestriction geom_data_restr,
                                       CeedOperator *op) const
{
  CeedQFunctionInfo info;
  const CeedInt dims = GetDimensions(ceed, basis, geom_data_restr, CEED_FE_SPACE_HCURL,
                                     "CurlCurlShapeIntegrator");
  switch (dims)
  {
    case 22:
      // Curl in 2D has a single component.
      info.apply_qf = f_shape_l2_22;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_l2_22_loc);
      break;
    case 32:
      info.apply_qf = f_shape_l2_32;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_l2_32_loc);
      break;
    case 33:
      info.apply_qf = f_shape_hdiv_33;
      info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hdiv_33_loc);
      break;
    default:
      MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                 << dims % 10 << ", " << dims / 10 << ") for CurlCurlShapeIntegrator!");
  }
  const bool scalar_curl = (dims % 10 < 3);
  info.trial_ops = EvalMode::Curl;
  info.test_ops = EvalMode::Curl;
  if (scalar_curl)
  {
    info.trial_ops |= EvalMode::Weight;
  }

  auto ctx = PopulateCoefficientContext(scalar_curl ? 1 : 3, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

void VectorFEBoundaryLFShapeIntegrator::Assemble(
    Ceed ceed, CeedElemRestriction restr, CeedBasis basis, CeedVector u, CeedVector w,
    CeedElemRestriction mesh_restr, CeedBasis mesh_basis, CeedVector geom_data,
    CeedElemRestriction geom_data_restr, CeedOperator *op) const
{
  CeedQFunctionInfo info;
  const CeedInt dims = GetDimensions(ceed, basis, geom_data_restr, CEED_FE_SPACE_HCURL,
                                     "VectorFEBoundaryLFShapeIntegrator");
  MFEM_VERIFY(dims == 32, "VectorFEBoundaryLFShapeIntegrator is only available on the "
                          "boundaries of 3D meshes!");
  info.apply_qf = f_shape_hcurl_lf_32;
  info.apply_qf_path = PalaceQFunctionRelativePath(f_shape_hcurl_lf_32_loc);
  info.test_ops = EvalMode::Interp;

  auto ctx = PopulateCoefficientContext(3, Q, transpose);
  AssembleCeedShapeIntegrator(info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar),
                              ceed, u, w, restr, basis, mesh_restr, mesh_basis, geom_data,
                              geom_data_restr, op);
}

}  // namespace palace
