// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
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
  int bdr_attr_max = mesh.bdr_attributes.Max();
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

template <typename T1, typename T2, typename T3, typename T4>
auto AddIntegrators(mfem::BilinearForm &a, T1 &df, T2 &f, T3 &dfb, T4 &fb)
{
  if (!df.empty())
  {
    a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(df));
  }
  if (!f.empty())
  {
    a.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(f));
  }
  if (!dfb.empty())
  {
    a.AddBoundaryIntegrator(new mfem::CurlCurlIntegrator(dfb));
  }
  if (!fb.empty())
  {
    a.AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(fb));
  }
}

template <typename T1, typename T2>
auto AddAuxIntegrators(mfem::BilinearForm &a, T1 &f, T2 &fb)
{
  if (!f.empty())
  {
    a.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(f));
  }
  if (!fb.empty())
  {
    a.AddBoundaryIntegrator(new mfem::MixedGradGradIntegrator(fb));
  }
}

}  // namespace

SpaceOperator::SpaceOperator(const IoData &iodata,
                             const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
  : assembly_level(iodata.solver.linear.mat_pa ? mfem::AssemblyLevel::PARTIAL
                                               : mfem::AssemblyLevel::LEGACY),
    skip_zeros(0), pc_mg(iodata.solver.linear.pc_mg),
    pc_lor(iodata.solver.linear.pc_mat_lor),
    pc_shifted(iodata.solver.linear.pc_mat_shifted), print_hdr(true), print_prec_hdr(true),
    dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())),
    nd_fecs(utils::ConstructFECollections<mfem::ND_FECollection>(
        pc_mg, pc_lor, iodata.solver.order, mesh.back()->Dimension())),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        pc_mg, false, iodata.solver.order, mesh.back()->Dimension())),
    rt_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    nd_fespaces(pc_mg ? utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
                            mesh, nd_fecs, &dbc_marker, &nd_dbc_tdof_lists)
                      : utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
                            *mesh.back(), *nd_fecs.back(), &dbc_marker,
                            &nd_dbc_tdof_lists.emplace_back())),
    h1_fespaces(pc_mg ? utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
                            mesh, h1_fecs, &dbc_marker, &h1_dbc_tdof_lists)
                      : utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
                            *mesh.back(), *h1_fecs.back(), &dbc_marker,
                            &h1_dbc_tdof_lists.emplace_back())),
    rt_fespace(mesh.back().get(), &rt_fec), mat_op(iodata, *mesh.back()),
    farfield_op(iodata, mat_op, *mesh.back()), surf_sigma_op(iodata, *mesh.back()),
    surf_z_op(iodata, *mesh.back()), lumped_port_op(iodata, GetH1Space()),
    wave_port_op(iodata, mat_op, GetNDSpace(), GetH1Space()),
    surf_j_op(iodata, GetH1Space())
{
  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_marker.Max() > 0)
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

std::unique_ptr<ParOperator>
SpaceOperator::GetSystemMatrix(SpaceOperator::OperatorType type,
                               Operator::DiagonalPolicy diag_policy)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " ND: {:d}, H1: {:d}, RT: {:d}\n",
               GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
               GetRTSpace().GlobalTrueVSize());
    print_hdr = false;
  }
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient df(sdim), f(sdim), fb(sdim);
  SumCoefficient dfb;
  switch (type)
  {
    case OperatorType::STIFFNESS:
      AddStiffnessCoefficients(1.0, df, f, fb);
      break;
    case OperatorType::DAMPING:
      AddDampingCoefficients(1.0, f, fb);
    case OperatorType::MASS:
      AddRealMassCoefficients(1.0, f, fb);
      break;
    case OperatorType::EXTRA:
    default:
      MFEM_ABORT("Invalid GetSystemMatrix matrix type for HypreParMatrix output!");
  }
  if (df.empty() && f.empty() && dfb.empty() && fb.empty())
  {
    return {};
  }
  auto a = std::make_unique<mfem::SymmetricBilinearForm>(&GetNDSpace());
  AddIntegrators(*a, df, f, dfb, fb);
  a->SetAssemblyLevel(assembly_level);
  a->Assemble(skip_zeros);
  a->Finalize(skip_zeros);
  auto A = std::make_unique<ParOperator>(std::move(a), GetNDSpace(), GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return A;
}

std::unique_ptr<ComplexParOperator>
SpaceOperator::GetComplexSystemMatrix(SpaceOperator::OperatorType type, double omega,
                                      Operator::DiagonalPolicy diag_policy)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " ND: {:d}, H1: {:d}, RT: {:d}\n",
               GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
               GetRTSpace().GlobalTrueVSize());
    print_hdr = false;
  }
  const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient dfr(sdim), dfi(sdim), fr(sdim), fi(sdim), fbr(sdim), fbi(sdim);
  SumCoefficient dfbr, dfbi;
  switch (type)
  {
    case OperatorType::STIFFNESS:
      MFEM_VERIFY(omega == 0.0, "GetComplexSystemMatrix for type OperatorType::STIFFNESS "
                                "does not use omega parameter!");
      AddStiffnessCoefficients(1.0, dfr, fr, fbr);
      break;
    case OperatorType::DAMPING:
      MFEM_VERIFY(omega == 0.0, "GetComplexSystemMatrix for type OperatorType::DAMPING "
                                "does not use omega parameter!");
      AddDampingCoefficients(1.0, fr, fbr);
      break;
    case OperatorType::MASS:
      MFEM_VERIFY(omega == 0.0, "GetComplexSystemMatrix for type OperatorType::MASS does "
                                "not use omega parameter!");
      AddRealMassCoefficients(1.0, fr, fbr);
      AddImagMassCoefficients(1.0, fi, fbi);
      break;
    case OperatorType::EXTRA:
      MFEM_VERIFY(omega > 0.0,
                  "GetComplexSystemMatrix for type OperatorType::EXTRA requires "
                  "use of omega parameter!");
      AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
      break;
  }
  bool has_real = false, has_imag = false;
  std::unique_ptr<mfem::SymmetricBilinearForm> ar, ai;
  if (!dfr.empty() || !fr.empty() || !dfbr.empty() || !fbr.empty())
  {
    has_real = true;
    ar = std::make_unique<mfem::SymmetricBilinearForm>(&GetNDSpace());
    AddIntegrators(*ar, dfr, fr, dfbr, fbr);
    ar->SetAssemblyLevel(assembly_level);
    ar->Assemble(skip_zeros);
    ar->Finalize(skip_zeros);
  }
  if (!dfi.empty() || !fi.empty() || !dfbi.empty() || !fbi.empty())
  {
    has_imag = true;
    ai = std::make_unique<mfem::SymmetricBilinearForm>(&GetNDSpace());
    AddIntegrators(*ai, dfi, fi, dfbi, fbi);
    ai->SetAssemblyLevel(assembly_level);
    ai->Assemble(skip_zeros);
    ai->Finalize(skip_zeros);
  }
  if (!has_real && !has_imag)
  {
    return {};
  }
  auto A = std::make_unique<ComplexParOperator>(
      std::make_unique<ComplexWrapperOperator>(std::move(ar), std::move(ai)), GetNDSpace(),
      GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return A;
}

std::unique_ptr<ParOperator> SpaceOperator::GetSystemMatrix(double a0, double a1, double a2,
                                                            const ParOperator *K,
                                                            const ParOperator *C,
                                                            const ParOperator *M)
{
  int height = -1, width = -1;
  if (K)
  {
    height = K->LocalOperator().Height();
    width = K->LocalOperator().Width();
  }
  else if (C)
  {
    height = C->LocalOperator().Height();
    width = C->LocalOperator().Width();
  }
  else if (M)
  {
    height = M->LocalOperator().Height();
    width = M->LocalOperator().Width();
  }
  MFEM_VERIFY(height >= 0 && width >= 0,
              "At least one argument to GetSystemMatrix must not be empty!");
  auto sum = std::make_unique<SumOperator>(height, width);
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
  auto A = std::make_unique<ParOperator>(std::move(sum), GetNDSpace(), GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), Operator::DiagonalPolicy::DIAG_ONE);
  return A;
}

std::unique_ptr<ComplexParOperator> SpaceOperator::GetComplexSystemMatrix(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2,
    const ComplexParOperator *K, const ComplexParOperator *C, const ComplexParOperator *M,
    const ComplexParOperator *A2)
{
  int height = -1, width = -1;
  if (K)
  {
    height = K->LocalOperator().Height();
    width = K->LocalOperator().Width();
  }
  else if (C)
  {
    height = C->LocalOperator().Height();
    width = C->LocalOperator().Width();
  }
  else if (M)
  {
    height = M->LocalOperator().Height();
    width = M->LocalOperator().Width();
  }
  else if (A2)
  {
    height = A2->LocalOperator().Height();
    width = A2->LocalOperator().Width();
  }
  MFEM_VERIFY(height >= 0 && width >= 0,
              "At least one argument to GetSystemMatrix must not be empty!");
  auto sum = std::make_unique<ComplexSumOperator>(height, width);
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
  auto A = std::make_unique<ComplexParOperator>(std::move(sum), GetNDSpace(), GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), Operator::DiagonalPolicy::DIAG_ONE);
  return A;
}

void SpaceOperator::GetPreconditionerMatrix(double a0, double a1, double a2, double a3,
                                            std::vector<std::unique_ptr<ParOperator>> &B,
                                            std::vector<std::unique_ptr<ParOperator>> &AuxB)
{
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == nd_fespaces.GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");
  for (int s = 0; s < 2; s++)
  {
    auto &B_ = (s == 0) ? B : AuxB;
    auto &fespaces = (s == 0) ? nd_fespaces : h1_fespaces;
    auto &dbc_tdof_lists = (s == 0) ? nd_dbc_tdof_lists : h1_dbc_tdof_lists;
    B_.clear();
    B_.reserve(fespaces.GetNumLevels());
    for (int l = 0; l < fespaces.GetNumLevels(); l++)
    {
      auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
      const int sdim = GetNDSpace().GetParMesh()->SpaceDimension();
      SumMatrixCoefficient df(sdim), f(sdim), fb(sdim);
      SumCoefficient dfb;
      AddStiffnessCoefficients(a0, df, f, fb);
      AddDampingCoefficients(a1, f, fb);
      AddRealMassCoefficients<MaterialPropertyType::PERMITTIVITY_ABS>(
          pc_shifted ? std::abs(a2) : a2, f, fb);
      AddExtraSystemBdrCoefficients(a3, dfb, dfb, fb, fb);
      auto b = std::make_unique<mfem::SymmetricBilinearForm>(&fespace_l);
      if (s == 0)
      {
        AddIntegrators(*b, df, f, dfb, fb);
      }
      else
      {
        // H1 auxiliary space matrix Gᵀ B G.
        AddAuxIntegrators(*b, f, fb);
      }
      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d}{}: {:d} unknowns", l, (s == 0) ? "" : " (auxiliary)",
                   fespace_l.GlobalTrueVSize());
      }
      if (pc_lor)
      {
        // After we construct the LOR discretization we deep copy the LOR matrix and the
        // original bilinear form and LOR discretization are no longer needed.
        mfem::Array<int> dummy_dbc_tdof_list;
        mfem::LORDiscretization lor(*b, dummy_dbc_tdof_list);
        auto b_lor = std::make_unique<mfem::SparseMatrix>(lor.GetAssembledMatrix());
        if (print_prec_hdr)
        {
          HYPRE_BigInt nnz = b_lor->NumNonZeroElems();
          Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
          Mpi::Print(", {:d} NNZ (LOR)\n", nnz);
        }
        B_.push_back(std::make_unique<ParOperator>(std::move(b_lor), fespace_l, fespace_l));
      }
      else
      {
        b->SetAssemblyLevel(assembly_level);
        b->Assemble(skip_zeros);
        b->Finalize(skip_zeros);
        if (print_prec_hdr)
        {
          if (assembly_level == mfem::AssemblyLevel::LEGACY)
          {
            HYPRE_BigInt nnz = b->SpMat().NumNonZeroElems();
            Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
            Mpi::Print(", {:d} NNZ\n", nnz);
          }
          else
          {
            Mpi::Print("\n");
          }
        }
        B_.push_back(std::make_unique<ParOperator>(std::move(b), fespace_l, fespace_l));
      }
      B_.back()->SetEssentialTrueDofs(dbc_tdof_lists[l],
                                      Operator::DiagonalPolicy::DIAG_ONE);
    }
  }
  print_prec_hdr = false;
}

std::unique_ptr<ParOperator> SpaceOperator::GetCurlMatrix()
{
  auto curl = std::make_unique<mfem::DiscreteLinearOperator>(&GetNDSpace(), &GetRTSpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  curl->SetAssemblyLevel(assembly_level);
  curl->Assemble();
  curl->Finalize();
  return std::make_unique<ParOperator>(std::move(curl), GetNDSpace(), GetRTSpace(), true);
}

std::unique_ptr<ComplexParOperator> SpaceOperator::GetComplexCurlMatrix()
{
  auto curl = std::make_unique<mfem::DiscreteLinearOperator>(&GetNDSpace(), &GetRTSpace());
  curl->AddDomainInterpolator(new mfem::CurlInterpolator);
  curl->SetAssemblyLevel(assembly_level);
  curl->Assemble();
  curl->Finalize();
  return std::make_unique<ComplexParOperator>(
      std::make_unique<ComplexWrapperOperator>(std::move(curl), nullptr), GetNDSpace(),
      GetRTSpace(), true);
}

std::unique_ptr<ParOperator> SpaceOperator::GetGradMatrix()
{
  auto grad = std::make_unique<mfem::DiscreteLinearOperator>(&GetH1Space(), &GetNDSpace());
  grad->AddDomainInterpolator(new mfem::GradientInterpolator);
  grad->SetAssemblyLevel(assembly_level);
  grad->Assemble();
  grad->Finalize();
  return std::make_unique<ParOperator>(std::move(grad), GetH1Space(), GetNDSpace(), true);
}

std::unique_ptr<ComplexParOperator> SpaceOperator::GetComplexGradMatrix()
{
  auto grad = std::make_unique<mfem::DiscreteLinearOperator>(&GetH1Space(), &GetNDSpace());
  grad->AddDomainInterpolator(new mfem::GradientInterpolator);
  grad->SetAssemblyLevel(assembly_level);
  grad->Assemble();
  grad->Finalize();
  return std::make_unique<ComplexParOperator>(
      std::make_unique<ComplexWrapperOperator>(std::move(grad), nullptr), GetH1Space(),
      GetNDSpace(), true);
}

void SpaceOperator::AddStiffnessCoefficients(double coef, SumMatrixCoefficient &df,
                                             SumMatrixCoefficient &f,
                                             SumMatrixCoefficient &fb)
{
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
    df.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef));
  }

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_LONDON_DEPTH;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef),
                     mat_op.GetLondonDepthMarker());
  }

  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coef, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coef, fb);
}

void SpaceOperator::AddDampingCoefficients(double coef, SumMatrixCoefficient &f,
                                           SumMatrixCoefficient &fb)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::CONDUCTIVITY;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef),
                     mat_op.GetConductivityMarker());
  }

  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coef, fb);
  surf_z_op.AddDampingBdrCoefficients(coef, fb);
  lumped_port_op.AddDampingBdrCoefficients(coef, fb);
}

template <MaterialPropertyType MatType>
void SpaceOperator::AddRealMassCoefficients(double coef, SumMatrixCoefficient &f,
                                            SumMatrixCoefficient &fb)
{
  f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef));

  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surf_z_op.AddMassBdrCoefficients(coef, fb);
  lumped_port_op.AddMassBdrCoefficients(coef, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coef, SumMatrixCoefficient &f,
                                            SumMatrixCoefficient &fb)
{
  // Contribution for loss tangent: ε => ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::PERMITTIVITY_IMAG;
    f.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(mat_op, coef),
                     mat_op.GetLossTangentMarker());
  }
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
  RHS.SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  return nnz;
}

bool SpaceOperator::GetExcitationVector(double omega, ComplexVector &RHS)
{
  // Frequency domain excitation vector: RHS = iω RHS1 + RHS2(ω).
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS = std::complex<double>(0.0, 0.0);
  bool nnz1 = AddExcitationVector1Internal(RHS.Real());
  RHS *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(omega, RHS);
  RHS.Real().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  RHS.Imag().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  RHS.SyncAlias();
  return nnz1 || nnz2;
}

bool SpaceOperator::GetExcitationVector1(ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later).
  RHS1.SetSize(GetNDSpace().GetTrueVSize());
  RHS1 = std::complex<double>(0.0, 0.0);
  bool nnz1 = AddExcitationVector1Internal(RHS1.Real());
  RHS1.Real().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  RHS1.SyncAlias();
  return nnz1;
}

bool SpaceOperator::GetExcitationVector2(double omega, ComplexVector &RHS2)
{
  RHS2.SetSize(GetNDSpace().GetTrueVSize());
  RHS2 = std::complex<double>(0.0, 0.0);
  bool nnz2 = AddExcitationVector2Internal(omega, RHS2);
  RHS2.Real().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  RHS2.Imag().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  RHS2.SyncAlias();
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
  v = std::complex<double>(1.0, 0.0);
  v.Real().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  v.SyncAlias();
}

void SpaceOperator::GetRandomInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  linalg::SetRandom(GetNDSpace().GetComm(), v);
  v.Real().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  v.Imag().SetSubVector(nd_dbc_tdof_lists.back(), 0.0);
  v.SyncAlias();
}

template void
SpaceOperator::AddRealMassCoefficients<MaterialPropertyType::PERMITTIVITY_REAL>(
    double, SumMatrixCoefficient &, SumMatrixCoefficient &);
template void
SpaceOperator::AddRealMassCoefficients<MaterialPropertyType::PERMITTIVITY_ABS>(
    double, SumMatrixCoefficient &, SumMatrixCoefficient &);

}  // namespace palace
