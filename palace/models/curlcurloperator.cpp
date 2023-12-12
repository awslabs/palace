// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "curlcurloperator.hpp"

#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

CurlCurlOperator::CurlCurlOperator(const IoData &iodata,
                                   const std::vector<std::unique_ptr<Mesh>> &mesh)
  : print_hdr(true), dbc_attr(SetUpBoundaryProperties(iodata, *mesh.back())),
    nd_fecs(fem::ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.order, mesh.back()->Get().Dimension(),
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_coarsen_type, false)),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Get().Dimension(),
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_coarsen_type, false)),
    rt_fec(std::make_unique<mfem::RT_FECollection>(iodata.solver.order - 1,
                                                   mesh.back()->Get().Dimension())),
    nd_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        iodata.solver.linear.mg_max_levels, mesh, nd_fecs, &dbc_attr, &dbc_tdof_lists)),
    h1_fespaces(fem::ConstructAuxiliaryFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        nd_fespaces, h1_fecs)),
    rt_fespace(nd_fespaces.GetFinestFESpace(), *mesh.back(), rt_fec.get()),
    mat_op(iodata, *mesh.back()), surf_j_op(iodata, GetH1Space())
{
  // Finalize setup.
  BilinearForm::pa_order_threshold = iodata.solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = iodata.solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = iodata.solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = iodata.solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = iodata.solver.q_order_extra;
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    std::sort(dbc_attr.begin(), dbc_attr.end());
    utils::PrettyPrint(dbc_attr);
  }
}

mfem::Array<int> CurlCurlOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                           const mfem::ParMesh &mesh)
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
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  return dbc_bcs;
}

void CurlCurlOperator::CheckBoundaryProperties()
{
  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto surf_j_marker = mesh::AttrToMarker(bdr_attr_max, surf_j_op.GetAttrList());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + surf_j_marker[i] <= 1,
                "Boundary attributes should not be specified with multiple BC!");
  }
}

namespace
{

void PrintHeader(const mfem::ParFiniteElementSpace &h1_fespace,
                 const mfem::ParFiniteElementSpace &nd_fespace,
                 const mfem::ParFiniteElementSpace &rt_fespace, bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1 (p = {:d}): {:d}, ND (p = {:d}): {:d}, RT (p = {:d}): {:d}\n Operator "
               "assembly level: {}\n",
               h1_fespace.GetMaxElementOrder(), h1_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder(), nd_fespace.GlobalTrueVSize(),
               rt_fespace.GetMaxElementOrder(), rt_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder() > BilinearForm::pa_order_threshold
                   ? "Partial"
                   : "Full");

    const auto &mesh = *nd_fespace.GetParMesh();
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : mesh::CheckElements(mesh).GetGeomTypes())
    {
      const auto *fe = nd_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support ND spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d})\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order);
    }

    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
}

}  // namespace

std::unique_ptr<Operator> CurlCurlOperator::GetStiffnessMatrix()
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);

  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm k(GetNDSpace());
  k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  k.AssembleQuadratureData();
  auto k_vec = k.Assemble(GetNDSpaces(), skip_zeros);
  auto K = std::make_unique<MultigridOperator>(GetNDSpaces().GetNumLevels());
  for (std::size_t l = 0; l < GetNDSpaces().GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = GetNDSpaces().GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 nd_fespace_l.GetMaxElementOrder(), nd_fespace_l.GlobalTrueVSize());
      if (const auto *k_spm = dynamic_cast<const mfem::SparseMatrix *>(k_vec[l].get()))
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
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), nd_fespace_l);
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }

  print_hdr = false;
  return K;
}

void CurlCurlOperator::GetExcitationVector(int idx, Vector &RHS)
{
  // Assemble the surface current excitation +J. The SurfaceCurrentOperator assembles -J
  // (meant for time or frequency domain Maxwell discretization, so we multiply by -1 to
  // retrieve +J).
  SumVectorCoefficient fb(GetNDSpace().SpaceDimension());
  surf_j_op.AddExcitationBdrCoefficients(idx, fb);
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS = 0.0;
  if (fb.empty())
  {
    return;
  }
  mfem::LinearForm rhs(&GetNDSpace().Get());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(false);
  rhs.Assemble();
  GetNDSpace().Get().GetProlongationMatrix()->AddMultTranspose(rhs, RHS, -1.0);
  linalg::SetSubVector(RHS, dbc_tdof_lists.back(), 0.0);
}

}  // namespace palace
