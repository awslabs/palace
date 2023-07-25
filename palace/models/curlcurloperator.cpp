// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "curlcurloperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

mfem::Array<int> SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh)
{
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!iodata.boundaries.pec.empty())
  {
    // Check that boundary attributes have been specified correctly.
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    bool first = true;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "PEC boundary attribute tags must be non-negative and correspond to "
      //             "attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr-1],
      //             "Unknown PEC boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
      {
        if (first)
        {
          Mpi::Print("\n");
          first = false;
        }
        Mpi::Warning("Unknown PEC boundary attribute {:d}!\nSolver will just ignore it!\n",
                     attr);
      }
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs, dbc_marker;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  mesh::AttrToMarker(bdr_attr_max, dbc_bcs, dbc_marker);
  return dbc_marker;
}

}  // namespace

CurlCurlOperator::CurlCurlOperator(const IoData &iodata,
                                   const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
  : pa_order_threshold(iodata.solver.pa_order_threshold), skip_zeros(0), print_hdr(true),
    dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())),
    nd_fecs(utils::ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    rt_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    nd_fespaces(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        pa_order_threshold, mesh, nd_fecs, dbc_marker, dbc_tdof_lists)),
    h1_fespaces(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        pa_order_threshold, mesh, h1_fecs)),
    rt_fespace(mesh.back().get(), &rt_fec), mat_op(iodata, *mesh.back()),
    surf_j_op(iodata, GetH1Space())
{
  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_marker.Size() && dbc_marker.Max() > 0)
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrintMarker(dbc_marker);
  }
}

void CurlCurlOperator::CheckBoundaryProperties()
{
  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const auto &surf_j_marker = surf_j_op.GetMarker();
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + surf_j_marker[i] <= 1,
                "Boundary attributes should not be specified with multiple BC!");
  }
}

std::unique_ptr<Operator> CurlCurlOperator::GetStiffnessMatrix()
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1: {:d}, ND: {:d}, RT: {:d}\n Operator assembly level: {}\n",
               GetH1Space().GlobalTrueVSize(), GetNDSpace().GlobalTrueVSize(),
               GetRTSpace().GlobalTrueVSize(),
               GetNDSpace().GetMaxElementOrder() > pa_order_threshold ? "Partial" : "Full");
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  auto K = std::make_unique<MultigridOperator>(nd_fespaces.GetNumLevels());
  for (int l = 0; l < nd_fespaces.GetNumLevels(); l++)
  {
    // Force coarse level operator to be fully assembled always.
    auto &nd_fespace_l = nd_fespaces.GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 nd_fespace_l.GetMaxElementOrder(), nd_fespace_l.GlobalTrueVSize());
    }
    constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
    MaterialPropertyCoefficient<MatType> muinv_func(mat_op);
    auto k = std::make_unique<mfem::SymmetricBilinearForm>(&nd_fespace_l);
    k->AddDomainIntegrator(new mfem::CurlCurlIntegrator(muinv_func));
    auto K_l = std::make_unique<ParOperator>(
        utils::AssembleOperator(std::move(k), true, (l > 0) ? pa_order_threshold : 100,
                                skip_zeros),
        nd_fespace_l);
    if (print_hdr)
    {
      if (const auto *k_spm =
              dynamic_cast<const mfem::SparseMatrix *>(&K_l->LocalOperator()))
      {
        HYPRE_BigInt nnz = k_spm->NumNonZeroElems();
        Mpi::GlobalSum(1, &nnz, nd_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }
  print_hdr = false;
  return K;
}

std::unique_ptr<Operator> CurlCurlOperator::GetCurlMatrix()
{
  // Partial assembly for this operator is only available with libCEED backend.
  auto curl = std::make_unique<mfem::DiscreteLinearOperator>(&GetNDSpace(), &GetRTSpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  return std::make_unique<ParOperator>(
      utils::AssembleOperator(std::move(curl), false, pa_order_threshold - 1), GetNDSpace(),
      GetRTSpace(), true);
}

void CurlCurlOperator::GetExcitationVector(int idx, Vector &RHS)
{
  // Assemble the surface current excitation +J. The SurfaceCurrentOperator assembles -J
  // (meant for time or frequency domain Maxwell discretization, so we multiply by -1 to
  // retrieve +J).
  SumVectorCoefficient fb(GetNDSpace().GetParMesh()->SpaceDimension());
  surf_j_op.AddExcitationBdrCoefficients(idx, fb);
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS = 0.0;
  if (fb.empty())
  {
    return;
  }
  mfem::LinearForm rhs(&GetNDSpace());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(false);
  rhs.Assemble();
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs, RHS, -1.0);
  linalg::SetSubVector(RHS, dbc_tdof_lists.back(), 0.0);
}

}  // namespace palace
