// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include "fem/multigrid.hpp"
#include "linalg/operator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

BoundaryModeOperator::BoundaryModeOperator(const IoData &iodata_,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh,
                                           const MaterialOperator &mat_op_in)
  : iodata(iodata_), solver_order(iodata_.solver.order), solve_mesh(mesh.back().get()),
    mat_op(mat_op_in)
{
  MFEM_VERIFY(solve_mesh->Dimension() == 2,
              "BoundaryMode solver requires a 2D mesh (waveguide cross-section). When "
              "\"Attributes\" is set the driver's Preprocess hook extracts it from "
              "the 3D parent before this operator is constructed.");

  // FE spaces.
  SetUpFESpaces(mesh);

  // Boundary operators.
  auto &pmesh = solve_mesh->Get();
  surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, mat_op, pmesh);
  farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, mat_op, pmesh);
  surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(iodata, mat_op, pmesh);

  // Frequency-independent block matrices. Atn is the ND/H1 gradient coupling,
  // Btn = -Atn^T, Btt is the ND mass.
  std::tie(Atnr, Atni) = mode_assembly::AssembleAtn(GetNDSpace(), GetH1Space(), mat_op);
  Btnr.reset(Atnr->Transpose());
  *Btnr *= -1.0;
  if (Atni)
  {
    Btni.reset(Atni->Transpose());
    *Btni *= -1.0;
  }
  auto [bttr_tmp, btti_tmp] = mode_assembly::AssembleBtt(GetNDSpace(), mat_op);
  Bttr = std::move(bttr_tmp);

  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
             GetNDSpace().GlobalTrueVSize() + GetH1Space().GlobalTrueVSize());
}

BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAtt(std::complex<double> omega, double sigma) const
{
  return mode_assembly::AssembleAtt(GetNDSpace(), mat_op, nullptr, *surf_z_op, *farfield_op,
                                    *surf_sigma_op, omega, sigma);
}

BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAnn(std::complex<double> omega) const
{
  return mode_assembly::AssembleAnn(GetH1Space(), mat_op, nullptr, *surf_z_op, *farfield_op,
                                    *surf_sigma_op, omega);
}

void BoundaryModeOperator::ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn,
                                                ComplexVector &et, ComplexVector &en) const
{
  mode_assembly::ApplyVDBackTransform(e0, kn, GetNDTrueVSize(), GetH1TrueVSize(), et, en);
}

std::complex<double>
BoundaryModeOperator::ComputePoyntingPower(double omega, std::complex<double> kn,
                                           const ComplexVector &et,
                                           const ComplexVector &en) const
{
  if (!Bttr)
  {
    return 0.0;
  }
  auto comm = GetComm();
  std::complex<double> P = 0.5 * std::conj(kn) / omega * linalg::Dot(comm, et, *Bttr, et);
  if (Atnr && en.Size() == GetH1TrueVSize())
  {
    ComplexWrapperOperator Atn(const_cast<mfem::HypreParMatrix *>(Atnr.get()),
                               const_cast<mfem::HypreParMatrix *>(Atni.get()));
    P += std::complex<double>(0.0, 1.0) / (2.0 * omega) * linalg::Dot(comm, en, Atn, et);
  }
  return P;
}

void BoundaryModeOperator::SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  const auto &mg = iodata.solver.linear;
  const int dim = solve_mesh->Dimension();

  // Dirichlet bdr attrs. The solve mesh already carries PEC on every edge that needs it:
  // for direct 2D, the user defined it; for submesh extraction, the driver folded the
  // inherited 3D boundary conditions into 2D attributes before this runs.
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
  // RT collection for estimator flux recovery (order = solver_order - 1 since the
  // recovered flux lives one polynomial degree below the primary ND space).
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
