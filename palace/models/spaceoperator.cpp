// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <set>
#include <type_traits>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

using namespace std::complex_literals;

SpaceOperator::SpaceOperator(const IoData &iodata,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : pc_mat_real(iodata.solver.linear.pc_mat_real),
    pc_mat_shifted(iodata.solver.linear.pc_mat_shifted), print_hdr(true),
    print_prec_hdr(true), dbc_attr(SetUpBoundaryProperties(iodata, *mesh.back())),
    nd_fecs(fem::ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        iodata.solver.order - 1, mesh.back()->Dimension(),
        iodata.solver.linear.estimator_mg ? iodata.solver.linear.mg_max_levels : 1,
        iodata.solver.linear.mg_coarsen_type, false)),
    nd_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        iodata.solver.linear.mg_max_levels, mesh, nd_fecs, &dbc_attr, &nd_dbc_tdof_lists)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &h1_dbc_tdof_lists)),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        iodata.solver.linear.estimator_mg ? iodata.solver.linear.mg_max_levels : 1, mesh,
        rt_fecs)),
    mat_op(iodata, *mesh.back()), farfield_op(iodata, mat_op, *mesh.back()),
    surf_sigma_op(iodata, mat_op, *mesh.back()), surf_z_op(iodata, mat_op, *mesh.back()),
    lumped_port_op(iodata, mat_op, *mesh.back()),
    wave_port_op(iodata, mat_op, GetNDSpace(), GetH1Space()),
    surf_j_op(iodata, *mesh.back())
{
  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet PEC BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

mfem::Array<int> SpaceOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                        const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!iodata.boundaries.pec.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "PEC boundary attribute tags must be non-negative and correspond to "
      //             "attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr - 1],
      //             "Unknown PEC boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown PEC boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  return dbc_bcs;
}

void SpaceOperator::CheckBoundaryProperties()
{
  // Mark selected boundary attributes from the mesh as having some Dirichlet, Neumann, or
  // mixed BC applied.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto farfield_marker = mesh::AttrToMarker(bdr_attr_max, farfield_op.GetAttrList());
  const auto surf_sigma_marker =
      mesh::AttrToMarker(bdr_attr_max, surf_sigma_op.GetAttrList());
  const auto surf_z_Rs_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetRsAttrList());
  const auto surf_z_Ls_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetLsAttrList());
  const auto lumped_port_Rs_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetRsAttrList());
  const auto lumped_port_Ls_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetLsAttrList());
  const auto wave_port_marker =
      mesh::AttrToMarker(bdr_attr_max, wave_port_op.GetAttrList());
  mfem::Array<int> aux_bdr_marker(dbc_marker.Size());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    aux_bdr_marker[i] =
        (dbc_marker[i] || farfield_marker[i] || surf_sigma_marker[i] ||
         surf_z_Rs_marker[i] || surf_z_Ls_marker[i] || lumped_port_Rs_marker[i] ||
         lumped_port_Ls_marker[i] || wave_port_marker[i]);
    if (aux_bdr_marker[i])
    {
      aux_bdr_attr.Append(i + 1);
    }
  }
  // aux_bdr_marker = 1;  // Mark all boundaries (including material interfaces
  //                      // added during mesh preprocessing)
  //                      // As tested, this does not eliminate all DC modes!
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    GetH1Spaces().GetFESpaceAtLevel(l).Get().GetEssentialTrueDofs(
        aux_bdr_marker, aux_bdr_tdof_lists.emplace_back());
  }

  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const auto surf_z_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetAttrList());
  const auto lumped_port_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetAttrList());
  const auto surf_j_marker = mesh::AttrToMarker(bdr_attr_max, surf_j_op.GetAttrList());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + farfield_marker[i] + surf_sigma_marker[i] +
                        surf_z_marker[i] + lumped_port_marker[i] + wave_port_marker[i] +
                        surf_j_marker[i] <=
                    1,
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
               (nd_fespace.GetMaxElementOrder() >= BilinearForm::pa_order_threshold)
                   ? "Partial"
                   : "Full");

    const auto &mesh = *nd_fespace.GetParMesh();
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes();
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : geom_types)
    {
      const auto *fe = nd_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support ND spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d}){}\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order,
                 (geom == geom_types.back()) ? "" : ",");
    }
  }
  print_hdr = false;
}

void AddIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *df,
                    const MaterialPropertyCoefficient *f,
                    const MaterialPropertyCoefficient *dfb,
                    const MaterialPropertyCoefficient *fb, bool assemble_q_data = false)
{
  if (df && !df->empty() && f && !f->empty())
  {
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(*df, *f);
  }
  else
  {
    if (df && !df->empty())
    {
      a.AddDomainIntegrator<CurlCurlIntegrator>(*df);
    }
    if (f && !f->empty())
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(*f);
    }
  }
  if (dfb && !dfb->empty() && fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(*dfb, *fb);
  }
  else
  {
    if (dfb && !dfb->empty())
    {
      a.AddBoundaryIntegrator<CurlCurlIntegrator>(*dfb);
    }
    if (fb && !fb->empty())
    {
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(*fb);
    }
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

void AddAuxIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *fb, bool assemble_q_data = false)
{
  if (f && !f->empty())
  {
    a.AddDomainIntegrator<DiffusionIntegrator>(*f);
  }
  if (fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<DiffusionIntegrator>(*fb);
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

auto AssembleOperator(const FiniteElementSpace &fespace,
                      const MaterialPropertyCoefficient *df,
                      const MaterialPropertyCoefficient *f,
                      const MaterialPropertyCoefficient *dfb,
                      const MaterialPropertyCoefficient *fb, bool skip_zeros = false,
                      bool assemble_q_data = false)
{
  BilinearForm a(fespace);
  AddIntegrators(a, df, f, dfb, fb, assemble_q_data);
  return a.Assemble(skip_zeros);
}

auto AssembleOperators(const FiniteElementSpaceHierarchy &fespaces,
                       const MaterialPropertyCoefficient *df,
                       const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *dfb,
                       const MaterialPropertyCoefficient *fb, bool skip_zeros = false,
                       bool assemble_q_data = false, std::size_t l0 = 0)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddIntegrators(a, df, f, dfb, fb, assemble_q_data);
  return a.Assemble(fespaces, skip_zeros, l0);
}

auto AssembleAuxOperators(const FiniteElementSpaceHierarchy &fespaces,
                          const MaterialPropertyCoefficient *f,
                          const MaterialPropertyCoefficient *fb, bool skip_zeros = false,
                          bool assemble_q_data = false, std::size_t l0 = 0)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddAuxIntegrators(a, f, fb, assemble_q_data);
  return a.Assemble(fespaces, skip_zeros, l0);
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(1.0, df, f);
  AddStiffnessBdrCoefficients(1.0, fb);
  int empty = (df.empty() && f.empty() && fb.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto k = AssembleOperator(GetNDSpace(), &df, &f, nullptr, &fb, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto K = std::make_unique<ComplexParOperator>(std::move(k), nullptr, GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
  else
  {
    auto K = std::make_unique<ParOperator>(std::move(k), GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute());
  AddDampingCoefficients(1.0, f);
  AddDampingBdrCoefficients(1.0, fb);
  int empty = (f.empty() && fb.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto c = AssembleOperator(GetNDSpace(), nullptr, &f, nullptr, &fb, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto C = std::make_unique<ComplexParOperator>(std::move(c), nullptr, GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
  else
  {
    auto C = std::make_unique<ParOperator>(std::move(c), GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
}

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute()), fi(mat_op.MaxCeedAttribute()),
      fbr(mat_op.MaxCeedBdrAttribute()), fbi(mat_op.MaxCeedBdrAttribute());
  AddRealMassCoefficients(1.0, fr);
  AddRealMassBdrCoefficients(1.0, fbr);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    AddImagMassCoefficients(1.0, fi);
  }
  int empty[2] = {(fr.empty() && fbr.empty()), (fi.empty() && fbi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(GetNDSpace(), nullptr, &fr, nullptr, &fbr, skip_zeros);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(GetNDSpace(), nullptr, &fi, nullptr, &fbi, skip_zeros);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M =
        std::make_unique<ComplexParOperator>(std::move(mr), std::move(mi), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
  else
  {
    auto M = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
  int empty[2] = {(dfbr.empty() && fbr.empty()), (dfbi.empty() && fbi.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> ar, ai;
  if (!empty[0])
  {
    ar = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbr, &fbr, skip_zeros);
  }
  if (!empty[1])
  {
    ai = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbi, &fbi, skip_zeros);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto A =
        std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
  else
  {
    MFEM_VERIFY(!ai, "Unexpected imaginary part in GetExtraSystemMatrix<Operator>!");
    auto A = std::make_unique<ParOperator>(std::move(ar), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
}

namespace
{

auto BuildParSumOperator(int h, int w, double a0, double a1, double a2,
                         const ParOperator *K, const ParOperator *C, const ParOperator *M,
                         const ParOperator *A2, const FiniteElementSpace &fespace)
{
  auto sum = std::make_unique<SumOperator>(h, w);
  if (K && a0 != 0.0)
  {
    sum->AddOperator(K->LocalOperator(), a0);
  }
  if (C && a1 != 0.0)
  {
    sum->AddOperator(C->LocalOperator(), a1);
  }
  if (M && a2 != 0.0)
  {
    sum->AddOperator(M->LocalOperator(), a2);
  }
  if (A2)
  {
    sum->AddOperator(A2->LocalOperator(), 1.0);
  }
  return std::make_unique<ParOperator>(std::move(sum), fespace);
}

auto BuildParSumOperator(int h, int w, std::complex<double> a0, std::complex<double> a1,
                         std::complex<double> a2, const ComplexParOperator *K,
                         const ComplexParOperator *C, const ComplexParOperator *M,
                         const ComplexParOperator *A2, const FiniteElementSpace &fespace)
{
  // Block 2 x 2 equivalent-real formulation for each term in the sum:
  //                    [ sumr ]  +=  [ ar  -ai ] [ Ar ]
  //                    [ sumi ]      [ ai   ar ] [ Ai ] .
  auto sumr = std::make_unique<SumOperator>(h, w);
  auto sumi = std::make_unique<SumOperator>(h, w);
  if (K)
  {
    if (a0.real() != 0.0)
    {
      if (K->LocalOperator().Real())
      {
        sumr->AddOperator(*K->LocalOperator().Real(), a0.real());
      }
      if (K->LocalOperator().Imag())
      {
        sumi->AddOperator(*K->LocalOperator().Imag(), a0.real());
      }
    }
    if (a0.imag() != 0.0)
    {
      if (K->LocalOperator().Imag())
      {
        sumr->AddOperator(*K->LocalOperator().Imag(), -a0.imag());
      }
      if (K->LocalOperator().Real())
      {
        sumi->AddOperator(*K->LocalOperator().Real(), a0.imag());
      }
    }
  }
  if (C && a1 != 0.0)
  {
    if (a1.real() != 0.0)
    {
      if (C->LocalOperator().Real())
      {
        sumr->AddOperator(*C->LocalOperator().Real(), a1.real());
      }
      if (C->LocalOperator().Imag())
      {
        sumi->AddOperator(*C->LocalOperator().Imag(), a1.real());
      }
    }
    if (a1.imag() != 0.0)
    {
      if (C->LocalOperator().Imag())
      {
        sumr->AddOperator(*C->LocalOperator().Imag(), -a1.imag());
      }
      if (C->LocalOperator().Real())
      {
        sumi->AddOperator(*C->LocalOperator().Real(), a1.imag());
      }
    }
  }
  if (M && a2 != 0.0)
  {
    if (a2.real() != 0.0)
    {
      if (M->LocalOperator().Real())
      {
        sumr->AddOperator(*M->LocalOperator().Real(), a2.real());
      }
      if (M->LocalOperator().Imag())
      {
        sumi->AddOperator(*M->LocalOperator().Imag(), a2.real());
      }
    }
    if (a2.imag() != 0.0)
    {
      if (M->LocalOperator().Imag())
      {
        sumr->AddOperator(*M->LocalOperator().Imag(), -a2.imag());
      }
      if (M->LocalOperator().Real())
      {
        sumi->AddOperator(*M->LocalOperator().Real(), a2.imag());
      }
    }
  }
  if (A2)
  {
    if (A2->LocalOperator().Real())
    {
      sumr->AddOperator(*A2->LocalOperator().Real(), 1.0);
    }
    if (A2->LocalOperator().Imag())
    {
      sumi->AddOperator(*A2->LocalOperator().Imag(), 1.0);
    }
  }
  return std::make_unique<ComplexParOperator>(std::move(sumr), std::move(sumi), fespace);
}

}  // namespace

template <typename OperType, typename ScalarType>
std::unique_ptr<OperType>
SpaceOperator::GetSystemMatrix(ScalarType a0, ScalarType a1, ScalarType a2,
                               const OperType *K, const OperType *C, const OperType *M,
                               const OperType *A2)
{
  using ParOperType =
      typename std::conditional<std::is_same<OperType, ComplexOperator>::value,
                                ComplexParOperator, ParOperator>::type;

  const auto *PtAP_K = (K) ? dynamic_cast<const ParOperType *>(K) : nullptr;
  const auto *PtAP_C = (C) ? dynamic_cast<const ParOperType *>(C) : nullptr;
  const auto *PtAP_M = (M) ? dynamic_cast<const ParOperType *>(M) : nullptr;
  const auto *PtAP_A2 = (A2) ? dynamic_cast<const ParOperType *>(A2) : nullptr;
  MFEM_VERIFY((!K || PtAP_K) && (!C || PtAP_C) && (!M || PtAP_M) && (!A2 || PtAP_A2),
              "SpaceOperator requires ParOperator or ComplexParOperator for system matrix "
              "construction!");

  int height = -1, width = -1;
  if (PtAP_K)
  {
    height = PtAP_K->LocalOperator().Height();
    width = PtAP_K->LocalOperator().Width();
  }
  else if (PtAP_C)
  {
    height = PtAP_C->LocalOperator().Height();
    width = PtAP_C->LocalOperator().Width();
  }
  else if (PtAP_M)
  {
    height = PtAP_M->LocalOperator().Height();
    width = PtAP_M->LocalOperator().Width();
  }
  else if (PtAP_A2)
  {
    height = PtAP_A2->LocalOperator().Height();
    width = PtAP_A2->LocalOperator().Width();
  }
  MFEM_VERIFY(height >= 0 && width >= 0,
              "At least one argument to GetSystemMatrix must not be empty!");

  auto A = BuildParSumOperator(height, width, a0, a1, a2, PtAP_K, PtAP_C, PtAP_M, PtAP_A2,
                               GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), Operator::DiagonalPolicy::DIAG_ONE);
  return A;
}

std::unique_ptr<Operator> SpaceOperator::GetInnerProductMatrix(double a0, double a2,
                                                               const ComplexOperator *K,
                                                               const ComplexOperator *M)
{
  const auto *PtAP_K = (K) ? dynamic_cast<const ComplexParOperator *>(K) : nullptr;
  const auto *PtAP_M = (M) ? dynamic_cast<const ComplexParOperator *>(M) : nullptr;
  MFEM_VERIFY(
      (!K || PtAP_K) && (!M || PtAP_M),
      "SpaceOperator requires ComplexParOperator for inner product matrix construction!");

  int height = -1, width = -1;
  if (PtAP_K)
  {
    height = PtAP_K->LocalOperator().Height();
    width = PtAP_K->LocalOperator().Width();
  }
  else if (PtAP_M)
  {
    height = PtAP_M->LocalOperator().Height();
    width = PtAP_M->LocalOperator().Width();
  }
  MFEM_VERIFY(height >= 0 && width >= 0,
              "At least one argument to GetInnerProductMatrix must not be empty!");

  auto sum = std::make_unique<SumOperator>(height, width);
  if (PtAP_K && a0 != 0.0)
  {
    MFEM_VERIFY(
        PtAP_K->LocalOperator().Real(),
        "Missing real part of stiffness matrix for inner product matrix construction!");
    sum->AddOperator(*PtAP_K->LocalOperator().Real(), a0);
  }
  if (PtAP_M && a2 != 0.0)
  {
    MFEM_VERIFY(PtAP_M->LocalOperator().Real(),
                "Missing real part of mass matrix for inner product matrix construction!");
    sum->AddOperator(*PtAP_M->LocalOperator().Real(), a2);
  }
  return std::make_unique<ParOperator>(std::move(sum), GetNDSpace());
}

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&br, std::unique_ptr<Operator> &&bi,
                           const FiniteElementSpace &fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&br,
                                     std::unique_ptr<Operator> &&bi,
                                     const FiniteElementSpace &fespace)
{
  MFEM_VERIFY(
      !bi,
      "Should not be constructing a real-valued ParOperator with non-zero imaginary part!");
  return std::make_unique<ParOperator>(std::move(br), fespace);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&br,
                                            std::unique_ptr<Operator> &&bi,
                                            const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(br), std::move(bi), fespace);
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetPreconditionerMatrix(double a0, double a1,
                                                                 double a2, double a3)
{
  // XX TODO: Handle complex coeff a0/a1/a2/a3 (like GetSystemMatrix)

  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh. When appropriate, we build the
  // preconditioner on all levels based on the actual complex-valued system matrix. The
  // coarse operator is always fully assembled.
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  MFEM_VERIFY(GetH1Spaces().GetNumLevels() == GetNDSpaces().GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");

  const auto n_levels = GetNDSpaces().GetNumLevels();
  std::vector<std::unique_ptr<Operator>> br_vec(n_levels), bi_vec(n_levels),
      br_aux_vec(n_levels), bi_aux_vec(n_levels);
  constexpr bool skip_zeros = false, assemble_q_data = false;
  if (std::is_same<OperType, ComplexOperator>::value && !pc_mat_real)
  {
    MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()),
        dfi(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
        fi(mat_op.MaxCeedAttribute()), dfbr(mat_op.MaxCeedBdrAttribute()),
        dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
        fbi(mat_op.MaxCeedBdrAttribute());
    AddStiffnessCoefficients(a0, dfr, fr);
    AddStiffnessBdrCoefficients(a0, fbr);
    AddDampingCoefficients(a1, fi);
    AddDampingBdrCoefficients(a1, fbi);
    AddRealMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
    AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
    AddImagMassCoefficients(a2, fi);
    AddExtraSystemBdrCoefficients(a3, dfbr, dfbi, fbr, fbi);
    int empty[2] = {(dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty()),
                    (dfi.empty() && fi.empty() && dfbi.empty() && fbi.empty())};
    Mpi::GlobalMin(2, empty, GetComm());
    if (!empty[0])
    {
      br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, skip_zeros,
                                 assemble_q_data);
      br_aux_vec =
          AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
    }
    if (!empty[1])
    {
      bi_vec = AssembleOperators(GetNDSpaces(), &dfi, &fi, &dfbi, &fbi, skip_zeros,
                                 assemble_q_data);
      bi_aux_vec =
          AssembleAuxOperators(GetH1Spaces(), &fi, &fbi, skip_zeros, assemble_q_data);
    }
  }
  else
  {
    MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()),
        fr(mat_op.MaxCeedAttribute()), dfbr(mat_op.MaxCeedBdrAttribute()),
        fbr(mat_op.MaxCeedBdrAttribute());
    AddStiffnessCoefficients(a0, dfr, fr);
    AddStiffnessBdrCoefficients(a0, fbr);
    AddDampingCoefficients(a1, fr);
    AddDampingBdrCoefficients(a1, fbr);
    AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
    AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
    AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
    int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty());
    Mpi::GlobalMin(1, &empty, GetComm());
    if (!empty)
    {
      br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, skip_zeros,
                                 assemble_q_data);
      br_aux_vec =
          AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
    }
  }

  auto B = std::make_unique<BaseMultigridOperator<OperType>>(n_levels);
  for (bool aux : {false, true})
  {
    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &fespace_l =
          aux ? GetH1Spaces().GetFESpaceAtLevel(l) : GetNDSpaces().GetFESpaceAtLevel(l);
      const auto &dbc_tdof_lists_l = aux ? h1_dbc_tdof_lists[l] : nd_dbc_tdof_lists[l];
      auto &br_l = aux ? br_aux_vec[l] : br_vec[l];
      auto &bi_l = aux ? bi_aux_vec[l] : bi_vec[l];
      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d}{} (p = {:d}): {:d} unknowns", l, aux ? " (auxiliary)" : "",
                   fespace_l.GetMaxElementOrder(), fespace_l.GlobalTrueVSize());
        const auto *b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(br_l.get());
        if (!b_spm)
        {
          b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(bi_l.get());
        }
        if (b_spm)
        {
          HYPRE_BigInt nnz = b_spm->NNZ();
          Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
          Mpi::Print(", {:d} NNZ\n", nnz);
        }
        else
        {
          Mpi::Print("\n");
        }
      }
      auto B_l =
          BuildLevelParOperator<OperType>(std::move(br_l), std::move(bi_l), fespace_l);
      B_l->SetEssentialTrueDofs(dbc_tdof_lists_l, Operator::DiagonalPolicy::DIAG_ONE);
      if (aux)
      {
        B->AddAuxiliaryOperator(std::move(B_l));
      }
      else
      {
        B->AddOperator(std::move(B_l));
      }
    }
  }

  print_prec_hdr = false;
  return B;
}

void SpaceOperator::AddStiffnessCoefficients(double coeff, MaterialPropertyCoefficient &df,
                                             MaterialPropertyCoefficient &f)
{
  // Contribution from material permeability.
  df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvPermeability(), coeff);

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvLondonDepth(), coeff);
  }
}

void SpaceOperator::AddStiffnessBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coeff, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddDampingCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(), coeff);
  }
}

void SpaceOperator::AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coeff, fb);
  surf_z_op.AddDampingBdrCoefficients(coeff, fb);
  lumped_port_op.AddDampingBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), coeff);
}

void SpaceOperator::AddRealMassBdrCoefficients(double coeff,
                                               MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surf_z_op.AddMassBdrCoefficients(coeff, fb);
  lumped_port_op.AddMassBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityImag(), coeff);
  }
}

void SpaceOperator::AddAbsMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityAbs(), coeff);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(double omega,
                                                  MaterialPropertyCoefficient &dfbr,
                                                  MaterialPropertyCoefficient &dfbi,
                                                  MaterialPropertyCoefficient &fbr,
                                                  MaterialPropertyCoefficient &fbi)
{
  // Contribution for second-order farfield boundaries and finite conductivity boundaries.
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  // Contribution for numeric wave ports.
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, Vector &RHS)
{
  // Time domain excitation vector.
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  bool nnz = AddExcitationVector1Internal(excitation_idx, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz;
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, double omega,
                                        ComplexVector &RHS)
{
  // Frequency domain excitation vector: RHS = iω RHS1 + RHS2(ω).
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS.Real());
  RHS *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz1 || nnz2;
}

bool SpaceOperator::GetExcitationVector1(int excitation_idx, ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later).
  RHS1.SetSize(GetNDSpace().GetTrueVSize());
  RHS1.UseDevice(true);
  RHS1 = 0.0;
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS1.Real());
  linalg::SetSubVector(RHS1.Real(), nd_dbc_tdof_lists.back(), 0.0);
  return nnz1;
}

bool SpaceOperator::GetExcitationVector2(int excitation_idx, double omega,
                                         ComplexVector &RHS2)
{
  RHS2.SetSize(GetNDSpace().GetTrueVSize());
  RHS2.UseDevice(true);
  RHS2 = 0.0;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, RHS2);
  linalg::SetSubVector(RHS2, nd_dbc_tdof_lists.back(), 0.0);
  return nnz2;
}

bool SpaceOperator::AddExcitationVector1Internal(int excitation_idx, Vector &RHS1)
{
  // Assemble the time domain excitation -g'(t) J or frequency domain excitation -iω J.
  // The g'(t) or iω factors are not accounted for here, they is accounted for in the time
  // integration or frequency sweep later.
  MFEM_VERIFY(RHS1.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");
  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(excitation_idx, fb);
  surf_j_op.AddExcitationBdrCoefficients(fb);  // No excitation_idx — currently in all
  int empty = (fb.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return false;
  }
  mfem::LinearForm rhs1(&GetNDSpace().Get());
  rhs1.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs1.UseFastAssembly(false);
  rhs1.UseDevice(false);
  rhs1.Assemble();
  rhs1.UseDevice(true);
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1, RHS1);
  return true;
}

bool SpaceOperator::AddExcitationVector2Internal(int excitation_idx, double omega,
                                                 ComplexVector &RHS2)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at the
  // specified frequency.
  MFEM_VERIFY(RHS2.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector2Internal!");
  SumVectorCoefficient fbr(GetMesh().SpaceDimension()), fbi(GetMesh().SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(excitation_idx, omega, fbr, fbi);
  int empty = (fbr.empty() && fbi.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return false;
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Real());
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Imag());
  }
  return true;
}

void SpaceOperator::GetConstantInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 1.0;
  linalg::SetSubVector(v.Real(), nd_dbc_tdof_lists.back(), 0.0);
}

void SpaceOperator::GetRandomInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  linalg::SetRandom(GetNDSpace().GetComm(), v);
  linalg::SetSubVector(v, nd_dbc_tdof_lists.back(), 0.0);
}

template std::unique_ptr<Operator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetSystemMatrix<Operator, double>(double, double, double, const Operator *,
                                                 const Operator *, const Operator *,
                                                 const Operator *);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetSystemMatrix<ComplexOperator, std::complex<double>>(
    std::complex<double>, std::complex<double>, std::complex<double>,
    const ComplexOperator *, const ComplexOperator *, const ComplexOperator *,
    const ComplexOperator *);

template std::unique_ptr<Operator>
SpaceOperator::GetPreconditionerMatrix<Operator>(double, double, double, double);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetPreconditionerMatrix<ComplexOperator>(double, double, double, double);

}  // namespace palace
