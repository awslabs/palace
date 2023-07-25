// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <type_traits>
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

using namespace std::complex_literals;

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

SpaceOperator::SpaceOperator(const IoData &iodata,
                             const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
  : pa_order_threshold(iodata.solver.pa_order_threshold), skip_zeros(0),
    pc_mat_real(iodata.solver.linear.pc_mat_real),
    pc_mat_shifted(iodata.solver.linear.pc_mat_shifted),
    pc_mat_lor(iodata.solver.linear.pc_mat_lor), print_hdr(true), print_prec_hdr(true),
    dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())),
    nd_fecs(utils::ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, pc_mat_lor)),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    rt_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    nd_fespaces(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        pa_order_threshold, mesh, nd_fecs, dbc_marker, nd_dbc_tdof_lists)),
    h1_fespaces(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        pa_order_threshold, mesh, h1_fecs, dbc_marker, h1_dbc_tdof_lists)),
    rt_fespace(mesh.back().get(), &rt_fec), mat_op(iodata, *mesh.back()),
    farfield_op(iodata, mat_op, *mesh.back()), surf_sigma_op(iodata, *mesh.back()),
    surf_z_op(iodata, *mesh.back()), lumped_port_op(iodata, GetH1Space()),
    wave_port_op(iodata, mat_op, GetNDSpace(), GetH1Space()),
    surf_j_op(iodata, GetH1Space())
{
  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_marker.Size() && dbc_marker.Max() > 0)
  {
    Mpi::Print("\nConfiguring Dirichlet PEC BC at attributes:\n");
    utils::PrettyPrintMarker(dbc_marker);
  }
}

void SpaceOperator::CheckBoundaryProperties()
{
  // Mark selected boundary attributes from the mesh as having some Dirichlet, Neumann, or
  // mixed BC applied.
  const auto &farfield_marker = farfield_op.GetMarker();
  const auto &surf_sigma_marker = surf_sigma_op.GetMarker();
  const auto &surf_z_Rs_marker = surf_z_op.GetRsMarker();
  const auto &surf_z_Ls_marker = surf_z_op.GetLsMarker();
  const auto &lumped_port_Rs_marker = lumped_port_op.GetRsMarker();
  const auto &lumped_port_Ls_marker = lumped_port_op.GetLsMarker();
  const auto &wave_port_marker = wave_port_op.GetMarker();
  aux_bdr_marker.SetSize(dbc_marker.Size());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    aux_bdr_marker[i] =
        (dbc_marker[i] || farfield_marker[i] || surf_sigma_marker[i] ||
         surf_z_Rs_marker[i] || surf_z_Ls_marker[i] || lumped_port_Rs_marker[i] ||
         lumped_port_Ls_marker[i] || wave_port_marker[i]);
  }
  // aux_bdr_marker = 1;  // Mark all boundaries (including material interfaces
  //                      // added during mesh preprocessing)
  //                      // As tested, this does not eliminate all DC modes!
  for (int l = 0; l < h1_fespaces.GetNumLevels(); l++)
  {
    h1_fespaces.GetFESpaceAtLevel(l).GetEssentialTrueDofs(
        aux_bdr_marker, aux_bdr_tdof_lists.emplace_back());
  }

  // A final check that no boundary attribute is assigned multiple boundary conditions. The
  // one exception is that a lumped port boundary attribute can be also be assigned some
  // other condition, in which case the fact that it is a port is just used for
  // postprocessing.
  const auto &surf_z_marker = surf_z_op.GetMarker();
  const auto &lumped_port_marker = lumped_port_op.GetMarker();
  const auto &surf_j_marker = surf_j_op.GetMarker();
  bool first = true;
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    if (lumped_port_marker[i])
    {
      if (dbc_marker[i])
      {
        if (first)
        {
          Mpi::Print("\n");
          first = false;
        }
        Mpi::Warning("Lumped port boundary {:d} also marked as PEC!\nBoundary "
                     "condition/excitation will be ignored!\n",
                     i + 1);
      }
    }
    else
    {
      MFEM_VERIFY(dbc_marker[i] + farfield_marker[i] + surf_sigma_marker[i] +
                          surf_z_marker[i] + wave_port_marker[i] + surf_j_marker[i] <=
                      1,
                  "Boundary attributes should not be specified with multiple BC!");
    }
  }
}

namespace
{

void PrintHeader(mfem::ParFiniteElementSpace &h1_fespace,
                 mfem::ParFiniteElementSpace &nd_fespace,
                 mfem::ParFiniteElementSpace &rt_fespace, int pa_order_threshold,
                 bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1: {:d}, ND: {:d}, RT: {:d}\n Operator assembly level: {}\n",
               h1_fespace.GlobalTrueVSize(), nd_fespace.GlobalTrueVSize(),
               rt_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder() > pa_order_threshold ? "Partial" : "Full");
  }
  print_hdr = false;
}

template <typename T1, typename T2, typename T3, typename T4>
std::unique_ptr<Operator> BuildOperator(mfem::ParFiniteElementSpace &fespace, T1 *df, T2 *f,
                                        T3 *dfb, T4 *fb, int pa_order_threshold,
                                        int skip_zeros, bool pc_lor = false)
{
  auto a = std::make_unique<mfem::SymmetricBilinearForm>(&fespace);
  if (df && !df->empty() && f && !f->empty() && mfem::DeviceCanUseCeed())
  {
    a->AddDomainIntegrator(new mfem::CurlCurlMassIntegrator(*df, *f));
  }
  else
  {
    if (df && !df->empty())
    {
      a->AddDomainIntegrator(new mfem::CurlCurlIntegrator(*df));
    }
    if (f && !f->empty())
    {
      a->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(*f));
    }
  }
  if (dfb && !dfb->empty() && fb && !fb->empty() && mfem::DeviceCanUseCeed())
  {
    a->AddBoundaryIntegrator(new mfem::CurlCurlMassIntegrator(*dfb, *fb));
  }
  else
  {
    if (dfb && !dfb->empty())
    {
      a->AddBoundaryIntegrator(new mfem::CurlCurlIntegrator(*dfb));
    }
    if (fb && !fb->empty())
    {
      a->AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(*fb));
    }
  }
  if (pc_lor)
  {
    // After we construct the LOR discretization we deep copy the LOR matrix and the
    // original bilinear form and LOR discretization are no longer needed.
    mfem::Array<int> dummy_dbc_tdof_list;
    mfem::LORDiscretization lor(*a, dummy_dbc_tdof_list);
    return std::make_unique<mfem::SparseMatrix>(lor.GetAssembledMatrix());
  }
  else
  {
    return utils::AssembleOperator(std::move(a), true, pa_order_threshold, skip_zeros);
  }
}

template <typename T1, typename T2>
std::unique_ptr<Operator> BuildAuxOperator(mfem::ParFiniteElementSpace &fespace, T1 *f,
                                           T2 *fb, int pa_order_threshold, int skip_zeros,
                                           bool pc_lor = false)
{
  auto a = std::make_unique<mfem::SymmetricBilinearForm>(&fespace);
  if (f && !f->empty())
  {
    a->AddDomainIntegrator(new mfem::DiffusionIntegrator(*f));
  }
  if (fb && !fb->empty())
  {
    a->AddBoundaryIntegrator(new mfem::DiffusionIntegrator(*fb));
  }
  if (pc_lor)
  {
    // After we construct the LOR discretization we deep copy the LOR matrix and the
    // original bilinear form and LOR discretization are no longer needed.
    mfem::Array<int> dummy_dbc_tdof_list;
    mfem::LORDiscretization lor(*a, dummy_dbc_tdof_list);
    return std::make_unique<mfem::SparseMatrix>(lor.GetAssembledMatrix());
  }
  else
  {
    return utils::AssembleOperator(std::move(a), true, pa_order_threshold, skip_zeros);
  }
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), pa_order_threshold, print_hdr);
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient df(sdim), f(sdim), fb(sdim);
  AddStiffnessCoefficients(1.0, df, f);
  AddStiffnessBdrCoefficients(1.0, fb);
  if (df.empty() && f.empty() && fb.empty())
  {
    return {};
  }

  auto k = BuildOperator(GetNDSpace(), &df, &f, (SumCoefficient *)nullptr, &fb,
                         pa_order_threshold, skip_zeros);
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
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), pa_order_threshold, print_hdr);
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient f(sdim), fb(sdim);
  AddDampingCoefficients(1.0, f);
  AddDampingBdrCoefficients(1.0, fb);
  if (f.empty() && fb.empty())
  {
    return {};
  }

  auto c = BuildOperator(GetNDSpace(), (SumCoefficient *)nullptr, &f,
                         (SumCoefficient *)nullptr, &fb, pa_order_threshold, skip_zeros);
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
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), pa_order_threshold, print_hdr);
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient fr(sdim), fi(sdim), fbr(sdim);
  AddRealMassCoefficients(1.0, fr);
  AddRealMassBdrCoefficients(1.0, fbr);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    AddImagMassCoefficients(1.0, fi);
  }
  if (fr.empty() && fbr.empty() && fi.empty())
  {
    return {};
  }

  std::unique_ptr<Operator> mr, mi;
  if (!fr.empty() || !fbr.empty())
  {
    mr = BuildOperator(GetNDSpace(), (SumCoefficient *)nullptr, &fr,
                       (SumCoefficient *)nullptr, &fbr, pa_order_threshold, skip_zeros);
  }
  if (!fi.empty())
  {
    mi = BuildOperator(GetNDSpace(), (SumCoefficient *)nullptr, &fi,
                       (SumCoefficient *)nullptr, (SumCoefficient *)nullptr,
                       pa_order_threshold, skip_zeros);
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
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), pa_order_threshold, print_hdr);
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient fbr(sdim), fbi(sdim);
  SumCoefficient dfbr, dfbi;
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
  if (dfbr.empty() && fbr.empty() && dfbi.empty() && fbi.empty())
  {
    return {};
  }

  std::unique_ptr<Operator> ar, ai;
  if (!dfbr.empty() || !fbr.empty())
  {
    ar = BuildOperator(GetNDSpace(), (SumCoefficient *)nullptr, (SumCoefficient *)nullptr,
                       &dfbr, &fbr, pa_order_threshold, skip_zeros);
  }
  if (!dfbi.empty() || !fbi.empty())
  {
    ai = BuildOperator(GetNDSpace(), (SumCoefficient *)nullptr, (SumCoefficient *)nullptr,
                       &dfbi, &fbi, pa_order_threshold, skip_zeros);
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
                         const ParOperator *A2, const mfem::ParFiniteElementSpace &fespace)
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
                         const ComplexParOperator *A2,
                         const mfem::ParFiniteElementSpace &fespace)
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
      if (K->LocalOperator().HasReal())
      {
        sumr->AddOperator(*K->LocalOperator().Real(), a0.real());
      }
      if (K->LocalOperator().HasImag())
      {
        sumi->AddOperator(*K->LocalOperator().Imag(), a0.real());
      }
    }
    if (a0.imag() != 0.0)
    {
      if (K->LocalOperator().HasImag())
      {
        sumr->AddOperator(*K->LocalOperator().Imag(), -a0.imag());
      }
      if (K->LocalOperator().HasReal())
      {
        sumi->AddOperator(*K->LocalOperator().Real(), a0.imag());
      }
    }
  }
  if (C && a1 != 0.0)
  {
    if (a1.real() != 0.0)
    {
      if (C->LocalOperator().HasReal())
      {
        sumr->AddOperator(*C->LocalOperator().Real(), a1.real());
      }
      if (C->LocalOperator().HasImag())
      {
        sumi->AddOperator(*C->LocalOperator().Imag(), a1.real());
      }
    }
    if (a1.imag() != 0.0)
    {
      if (C->LocalOperator().HasImag())
      {
        sumr->AddOperator(*C->LocalOperator().Imag(), -a1.imag());
      }
      if (C->LocalOperator().HasReal())
      {
        sumi->AddOperator(*C->LocalOperator().Real(), a1.imag());
      }
    }
  }
  if (M && a2 != 0.0)
  {
    if (a2.real() != 0.0)
    {
      if (M->LocalOperator().HasReal())
      {
        sumr->AddOperator(*M->LocalOperator().Real(), a2.real());
      }
      if (M->LocalOperator().HasImag())
      {
        sumi->AddOperator(*M->LocalOperator().Imag(), a2.real());
      }
    }
    if (a2.imag() != 0.0)
    {
      if (M->LocalOperator().HasImag())
      {
        sumr->AddOperator(*M->LocalOperator().Imag(), -a2.imag());
      }
      if (M->LocalOperator().HasReal())
      {
        sumi->AddOperator(*M->LocalOperator().Real(), a2.imag());
      }
    }
  }
  if (A2)
  {
    if (A2->LocalOperator().HasReal())
    {
      sumr->AddOperator(*A2->LocalOperator().Real(), 1.0);
    }
    if (A2->LocalOperator().HasImag())
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
        PtAP_K->LocalOperator().HasReal(),
        "Missing real part of stiffness matrix for inner product matrix construction!");
    sum->AddOperator(*PtAP_K->LocalOperator().Real(), a0);
  }
  if (PtAP_M && a2 != 0.0)
  {
    MFEM_VERIFY(PtAP_M->LocalOperator().HasReal(),
                "Missing real part of mass matrix for inner product matrix construction!");
    sum->AddOperator(*PtAP_M->LocalOperator().Real(), a2);
  }
  return std::make_unique<ParOperator>(std::move(sum), GetNDSpace());
}

namespace
{

auto BuildLevelOperator(const MultigridOperator &B, std::unique_ptr<Operator> &&br,
                        std::unique_ptr<Operator> &&bi,
                        const mfem::ParFiniteElementSpace &fespace)
{
  return std::make_unique<ParOperator>(std::move(br), fespace);
}

auto BuildLevelOperator(const ComplexMultigridOperator &B, std::unique_ptr<Operator> &&br,
                        std::unique_ptr<Operator> &&bi,
                        const mfem::ParFiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(br), std::move(bi), fespace);
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetPreconditionerMatrix(double a0, double a1,
                                                                 double a2, double a3)
{
  // XX TODO: Test complex PC matrix assembly for l == 0 if coarse solve supports it
  // XX TODO: Handle complex coeff a0/a1/a2/a3 (like GetSystemMatrix)
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == nd_fespaces.GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");
  auto B = std::make_unique<BaseMultigridOperator<OperType>>(nd_fespaces.GetNumLevels());
  for (int s = 0; s < 2; s++)
  {
    auto &fespaces = (s == 0) ? nd_fespaces : h1_fespaces;
    auto &dbc_tdof_lists = (s == 0) ? nd_dbc_tdof_lists : h1_dbc_tdof_lists;
    for (int l = 0; l < fespaces.GetNumLevels(); l++)
    {
      // Force coarse level operator to be fully assembled always.
      auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d}{} (p = {:d}): {:d} unknowns", l,
                   (s == 0) ? "" : " (auxiliary)", fespace_l.GetMaxElementOrder(),
                   fespace_l.GlobalTrueVSize());
      }
      const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
      SumMatrixCoefficient dfr(sdim), fr(sdim), fi(sdim), fbr(sdim), fbi(sdim);
      SumCoefficient dfbr, dfbi;
      if (!std::is_same<OperType, ComplexOperator>::value || pc_mat_real || l == 0)
      {
        // Real-valued system matrix (approximation) for preconditioning.
        AddStiffnessCoefficients(a0, dfr, fr);
        AddStiffnessBdrCoefficients(a0, fbr);
        AddDampingCoefficients(a1, fr);
        AddDampingBdrCoefficients(a1, fbr);
        AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
        AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
        AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
      }
      else
      {
        // Build preconditioner based on the actual complex-valued system matrix.
        AddStiffnessCoefficients(a0, dfr, fr);
        AddStiffnessBdrCoefficients(a0, fbr);
        AddDampingCoefficients(a1, fi);
        AddDampingBdrCoefficients(a1, fbi);
        AddRealMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
        AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
        AddImagMassCoefficients(a2, fi);
        AddExtraSystemBdrCoefficients(a3, dfbr, dfbi, fbr, fbi);
      }

      std::unique_ptr<Operator> br, bi;
      if (!dfr.empty() || !fr.empty() || !dfbr.empty() || !fbr.empty())
      {
        br =
            (s == 0)
                ? BuildOperator(fespace_l, &dfr, &fr, &dfbr, &fbr,
                                (l > 0) ? pa_order_threshold : 100, skip_zeros, pc_mat_lor)
                : BuildAuxOperator(fespace_l, &fr, &fbr, (l > 0) ? pa_order_threshold : 100,
                                   skip_zeros, pc_mat_lor);
      }
      if (!fi.empty() || !dfbi.empty() || !fbi.empty())
      {
        bi =
            (s == 0)
                ? BuildOperator(fespace_l, (SumCoefficient *)nullptr, &fi, &dfbi, &fbi,
                                (l > 0) ? pa_order_threshold : 100, skip_zeros, pc_mat_lor)
                : BuildAuxOperator(fespace_l, &fi, &fbi, (l > 0) ? pa_order_threshold : 100,
                                   skip_zeros, pc_mat_lor);
      }
      if (print_prec_hdr)
      {
        if (const auto *br_spm = dynamic_cast<const mfem::SparseMatrix *>(br.get()))
        {
          HYPRE_BigInt nnz = br_spm->NumNonZeroElems();
          Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
          Mpi::Print(", {:d} NNZ{}\n", nnz, pc_mat_lor ? " (LOR)" : "");
        }
        else
        {
          Mpi::Print("\n");
        }
      }
      auto B_l = BuildLevelOperator(*B, std::move(br), std::move(bi), fespace_l);
      B_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
      if (s == 0)
      {
        B->AddOperator(std::move(B_l));
      }
      else
      {
        B->AddAuxiliaryOperator(std::move(B_l));
      }
    }
  }
  print_prec_hdr = false;
  return B;
}

namespace
{

auto BuildCurl(mfem::ParFiniteElementSpace &nd_fespace,
               mfem::ParFiniteElementSpace &rt_fespace, int pa_order_threshold)
{
  // Partial assembly for this operator is only available with libCEED backend.
  auto curl = std::make_unique<mfem::DiscreteLinearOperator>(&nd_fespace, &rt_fespace);
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  return utils::AssembleOperator(std::move(curl), false, pa_order_threshold - 1);
}

auto BuildGrad(mfem::ParFiniteElementSpace &h1_fespace,
               mfem::ParFiniteElementSpace &nd_fespace, int pa_order_threshold)
{
  auto grad = std::make_unique<mfem::DiscreteLinearOperator>(&h1_fespace, &nd_fespace);
  grad->AddDomainInterpolator(new mfem::GradientInterpolator);
  return utils::AssembleOperator(std::move(grad), true, pa_order_threshold);
}

}  // namespace

template <>
std::unique_ptr<Operator> SpaceOperator::GetCurlMatrix()
{
  return std::make_unique<ParOperator>(
      BuildCurl(GetNDSpace(), GetRTSpace(), pa_order_threshold), GetNDSpace(), GetRTSpace(),
      true);
}

template <>
std::unique_ptr<ComplexOperator> SpaceOperator::GetCurlMatrix()
{
  return std::make_unique<ComplexParOperator>(
      BuildCurl(GetNDSpace(), GetRTSpace(), pa_order_threshold), nullptr, GetNDSpace(),
      GetRTSpace(), true);
}

template <>
std::unique_ptr<Operator> SpaceOperator::GetGradMatrix()
{
  return std::make_unique<ParOperator>(
      BuildGrad(GetH1Space(), GetNDSpace(), pa_order_threshold), GetH1Space(), GetNDSpace(),
      true);
}

template <>
std::unique_ptr<ComplexOperator> SpaceOperator::GetGradMatrix()
{
  return std::make_unique<ComplexParOperator>(
      BuildGrad(GetH1Space(), GetNDSpace(), pa_order_threshold), nullptr, GetH1Space(),
      GetNDSpace(), true);
}

void SpaceOperator::AddStiffnessCoefficients(double coef, SumMatrixCoefficient &df,
                                             SumMatrixCoefficient &f)
{
  constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
  df.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef));

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    constexpr auto MatTypeL = MaterialPropertyType::INV_LONDON_DEPTH;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatTypeL>>(mat_op, coef),
                     mat_op.GetLondonDepthMarker());
  }
}

void SpaceOperator::AddStiffnessBdrCoefficients(double coef, SumMatrixCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coef, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coef, fb);
}

void SpaceOperator::AddDampingCoefficients(double coef, SumMatrixCoefficient &f)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    constexpr auto MatType = MaterialPropertyType::CONDUCTIVITY;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef),
                     mat_op.GetConductivityMarker());
  }
}

void SpaceOperator::AddDampingBdrCoefficients(double coef, SumMatrixCoefficient &fb)
{
  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coef, fb);
  surf_z_op.AddDampingBdrCoefficients(coef, fb);
  lumped_port_op.AddDampingBdrCoefficients(coef, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coef, SumMatrixCoefficient &f)
{
  constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_REAL;
  f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef));
}

void SpaceOperator::AddRealMassBdrCoefficients(double coef, SumMatrixCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surf_z_op.AddMassBdrCoefficients(coef, fb);
  lumped_port_op.AddMassBdrCoefficients(coef, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coef, SumMatrixCoefficient &f)
{
  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_IMAG;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef),
                     mat_op.GetLossTangentMarker());
  }
}

void SpaceOperator::AddAbsMassCoefficients(double coef, SumMatrixCoefficient &f)
{
  constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_ABS;
  f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef));
}

void SpaceOperator::AddExtraSystemBdrCoefficients(double omega, SumCoefficient &dfbr,
                                                  SumCoefficient &dfbi,
                                                  SumMatrixCoefficient &fbr,
                                                  SumMatrixCoefficient &fbi)
{
  // Contribution for second-order farfield boundaries and finite conductivity boundaries.
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  // Contribution for numeric wave ports.
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

bool SpaceOperator::GetExcitationVector(Vector &RHS)
{
  // Time domain excitation vector.
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS = 0.0;
  bool nnz = AddExcitationVector1Internal(RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz;
}

bool SpaceOperator::GetExcitationVector(double omega, ComplexVector &RHS)
{
  // Frequency domain excitation vector: RHS = iω RHS1 + RHS2(ω).
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS = 0.0;
  bool nnz1 = AddExcitationVector1Internal(RHS.Real());
  RHS *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(omega, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz1 || nnz2;
}

bool SpaceOperator::GetExcitationVector1(ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later).
  RHS1.SetSize(GetNDSpace().GetTrueVSize());
  RHS1 = 0.0;
  bool nnz1 = AddExcitationVector1Internal(RHS1.Real());
  linalg::SetSubVector(RHS1.Real(), nd_dbc_tdof_lists.back(), 0.0);
  return nnz1;
}

bool SpaceOperator::GetExcitationVector2(double omega, ComplexVector &RHS2)
{
  RHS2.SetSize(GetNDSpace().GetTrueVSize());
  RHS2 = 0.0;
  bool nnz2 = AddExcitationVector2Internal(omega, RHS2);
  linalg::SetSubVector(RHS2, nd_dbc_tdof_lists.back(), 0.0);
  return nnz2;
}

bool SpaceOperator::AddExcitationVector1Internal(Vector &RHS1)
{
  // Assemble the time domain excitation -g'(t) J or frequency domain excitation -iω J.
  // The g'(t) or iω factors are not accounted for here, they is accounted for in the time
  // integration or frequency sweep later.
  MFEM_VERIFY(RHS1.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");
  SumVectorCoefficient fb(GetNDSpace().GetParMesh()->SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(fb);
  surf_j_op.AddExcitationBdrCoefficients(fb);
  if (fb.empty())
  {
    return false;
  }
  mfem::LinearForm rhs1(&GetNDSpace());
  rhs1.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs1.UseFastAssembly(false);
  rhs1.Assemble();
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1, RHS1);
  return true;
}

bool SpaceOperator::AddExcitationVector2Internal(double omega, ComplexVector &RHS2)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at the
  // specified frequency.
  MFEM_VERIFY(RHS2.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector2Internal!");
  SumVectorCoefficient fbr(GetNDSpace().GetParMesh()->SpaceDimension()),
      fbi(GetNDSpace().GetParMesh()->SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(omega, fbr, fbi);
  if (fbr.empty() && fbi.empty())
  {
    return false;
  }
  mfem::LinearForm rhs2r(&GetNDSpace()), rhs2i(&GetNDSpace());
  rhs2r.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr));
  rhs2i.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
  rhs2r.UseFastAssembly(false);
  rhs2i.UseFastAssembly(false);
  rhs2r.Assemble();
  rhs2i.Assemble();
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2r, RHS2.Real());
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2i, RHS2.Imag());
  return true;
}

void SpaceOperator::GetConstantInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v = 1.0;
  linalg::SetSubVector(v.Real(), nd_dbc_tdof_lists.back(), 0.0);
}

void SpaceOperator::GetRandomInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
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
