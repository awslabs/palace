// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include "fem/multigrid.hpp"
#include "linalg/operator.hpp"
#include "models/modeoperatorassembly.hpp"
#include "models/surfaceterms.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/geodata.hpp"

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
  surface_terms = std::make_unique<SurfaceTerms>(iodata, mat_op, solve_mesh->Get());
  CheckBoundaryProperties();

  // Frequency-independent block matrices. Atn is the ND/H1 gradient coupling,
  // Btn = -Atn^T, Btt is the ND mass, all with the in-plane inverse permeability rotated
  // by 90 degrees (mode_assembly::GetTransverseInvPermeability).
  std::tie(Atnr, Atni) =
      mode_assembly::AssembleAtn(GetNDSpace(), GetH1Space(), mat_op, nullptr);
  Btnr.reset(Atnr->Transpose());
  *Btnr *= -1.0;
  if (Atni)
  {
    Btni.reset(Atni->Transpose());
    *Btni *= -1.0;
  }
  auto [bttr_tmp, btti_tmp] = mode_assembly::AssembleBtt(GetNDSpace(), mat_op, nullptr);
  Bttr = std::move(bttr_tmp);

  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
             GetNDSpace().GlobalTrueVSize() + GetH1Space().GlobalTrueVSize());
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
  // P = ½ ∫ (Eₜ × Hₜ*)·ẑ for E = (Eₜ + ẑEₙ) e^{-ikₙz}, whose transverse magnetic field is
  // Hₜ = μ⁻¹ (kₙ ẑ×Eₜ + i ∇ₜEₙ×ẑ) / ω, so that, with (ẑ×a)·μ⁻¹(ẑ×b) = a·νb,
  //   P = (1/2ω) ∫ (kₙ* Eₜ·νEₜ* + i Eₜ·ν∇ₜEₙ*)  (real μ)
  // for the in-plane inverse permeability rotated by 90 degrees, ν, of Btt and Atn
  // (mode_assembly::GetTransverseInvPermeability). With Atn = -(ν∇ₜφₙ, φₜ),
  // etᴴ Atn en = -∫ Eₜ*·ν∇ₜEₙ, the conjugate of minus the second integral. Without the
  // conjugate, the second term is only right when Eₙ is in quadrature with Eₜ, as for a
  // lossless mode: for a lossy mode, it adds a spurious imaginary part to P, which inflates
  // |P| and lowers Z_PV = Re{|V|² / 2P}. The ½ is the time-average weight of the peak
  // phasors.
  constexpr double weight = electromagnetics::TimeAverageWeight(true);
  std::complex<double> P =
      weight * std::conj(kn) / omega * linalg::Dot(comm, et, *Bttr, et);
  if (Atnr && en.Size() == GetH1TrueVSize())
  {
    ComplexWrapperOperator Atn(const_cast<mfem::HypreParMatrix *>(Atnr.get()),
                               const_cast<mfem::HypreParMatrix *>(Atni.get()));
    P -= std::complex<double>(0.0, weight) / omega *
         std::conj(linalg::Dot(comm, en, Atn, et));
  }
  return P;
}

void BoundaryModeOperator::CheckBoundaryProperties() const
{
  // No boundary attribute may be assigned more than one of the boundary conditions of the
  // 2D mode solve (PEC, and the absorbing, finite conductivity, impedance and rational
  // impedance boundaries of SurfaceTerms), as for the 3D operators
  // (SpaceOperator::CheckBoundaryProperties): the terms of each would be added on the
  // boundary. WavePortPEC is not counted: as in the mode solve of a wave port, it may
  // override another boundary condition of the attribute, whose terms then only act on
  // eliminated essential DOFs.
  const auto &pmesh = solve_mesh->Get();
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> pec_bcs;
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr > 0 && attr <= bdr_attr_max)
    {
      pec_bcs.Append(attr);
    }
  }
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, pec_bcs);
  const auto farfield_marker =
      mesh::AttrToMarker(bdr_attr_max, surface_terms->GetFarfieldOp().GetAttrList());
  const auto surf_sigma_marker = mesh::AttrToMarker(
      bdr_attr_max, surface_terms->GetSurfaceConductivityOp().GetAttrList());
  const auto surf_z_marker = mesh::AttrToMarker(
      bdr_attr_max, surface_terms->GetSurfaceImpedanceOp().GetAttrList());
  const auto surf_rz_marker = mesh::AttrToMarker(
      bdr_attr_max, surface_terms->GetRationalImpedanceOp().GetAttrList());
  for (int i = 0; i < bdr_attr_max; i++)
  {
    MFEM_VERIFY(dbc_marker[i] + farfield_marker[i] + surf_sigma_marker[i] +
                        surf_z_marker[i] + surf_rz_marker[i] <=
                    1,
                "Boundary attribute " << i + 1
                                      << " should not be specified with multiple BC!");
  }
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
