// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <complex>
#include "linalg/petsc.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"
#include "utils/mfemintegrators.hpp"
#include "utils/mfemoperators.hpp"
#include "utils/multigrid.hpp"
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
auto AddIntegrators(mfem::ParBilinearForm &a, T1 &df, T2 &f, T3 &dfb, T4 &fb)
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
auto AddAuxIntegrators(mfem::ParBilinearForm &a, T1 &f, T2 &fb)
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
  : dbc_marker(SetUpBoundaryProperties(iodata, *mesh.back())), skip_zeros(0),
    pc_gmg(iodata.solver.linear.mat_gmg), pc_lor(iodata.solver.linear.mat_lor),
    pc_shifted(iodata.solver.linear.mat_shifted), print_hdr(true),
    nd_fecs(utils::ConstructFECollections<mfem::ND_FECollection>(
        pc_gmg, pc_lor, iodata.solver.order, mesh.back()->Dimension())),
    h1_fecs(utils::ConstructFECollections<mfem::H1_FECollection>(
        pc_gmg, false, iodata.solver.order, mesh.back()->Dimension())),
    rt_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    nd_fespaces(pc_gmg ? utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
                             mesh, nd_fecs, dbc_marker)
                       : utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
                             *mesh.back(), *nd_fecs.back())),
    h1_fespaces(pc_gmg ? utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
                             mesh, h1_fecs, dbc_marker)
                       : utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
                             *mesh.back(), *h1_fecs.back())),
    rt_fespace(mesh.back().get(), &rt_fec), mat_op(iodata, *mesh.back()),
    farfield_op(iodata, mat_op, *mesh.back()), surf_sigma_op(iodata, *mesh.back()),
    surf_z_op(iodata, *mesh.back()), lumped_port_op(iodata, h1_fespaces.GetFinestFESpace()),
    wave_port_op(iodata, mat_op, nd_fespaces.GetFinestFESpace(),
                 h1_fespaces.GetFinestFESpace()),
    surf_j_op(iodata, h1_fespaces.GetFinestFESpace())
{
  // Finalize setup.
  CheckBoundaryProperties();
  nd_fespaces.GetFinestFESpace().GetEssentialTrueDofs(dbc_marker, dbc_tdof_list);

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

void SpaceOperator::PrintHeader()
{
  if (print_hdr)
  {
    Mpi::Print("\nConfiguring system matrices, number of global unknowns: {:d}\n",
               nd_fespaces.GetFinestFESpace().GlobalTrueVSize());
    print_hdr = false;
  }
}

int SpaceOperator::GetNDof()
{
  int ndof = GetNDSpace().GetTrueVSize();
  Mpi::GlobalSum(1, &ndof, Mpi::World());
  return ndof;
}

std::unique_ptr<petsc::PetscParMatrix>
SpaceOperator::GetSystemMatrixPetsc(SpaceOperator::OperatorType type, double omega,
                                    mfem::Operator::DiagonalPolicy ess_diag, bool print)
{
  // Construct the frequency-dependent complex linear system matrix:
  //                 A = K + iω C - ω² (Mr + i Mi) + A2(ω)
  // or any one of its terms.
  const int sdim = nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient dfr(sdim), dfi(sdim), fr(sdim), fi(sdim), fbr(sdim), fbi(sdim);
  SumCoefficient dfbr, dfbi;
  std::string str;
  switch (type)
  {
    case OperatorType::COMPLETE:
      AddStiffnessCoefficients(1.0, dfr, fr, fbr);
      AddDampingCoefficients(omega, fi, fbi);
      AddRealMassCoefficients(-omega * omega, false, fr, fbr);
      AddImagMassCoefficients(-omega * omega, fi, fbi);
      AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
      str = "A";
      break;
    case OperatorType::STIFFNESS:
      MFEM_VERIFY(omega == 0.0,
                  "GetSystemMatrix for type OperatorType::STIFFNESS does not use omega "
                  "parameter!");
      AddStiffnessCoefficients(1.0, dfr, fr, fbr);
      str = "K";
      break;
    case OperatorType::MASS:
      MFEM_VERIFY(
          omega == 0.0,
          "GetSystemMatrix for type OperatorType::MASS does not use omega parameter!");
      AddRealMassCoefficients(1.0, false, fr, fbr);
      AddImagMassCoefficients(1.0, fi, fbi);
      str = "M";
      break;
    case OperatorType::DAMPING:
      MFEM_VERIFY(
          omega == 0.0,
          "GetSystemMatrix for type OperatorType::DAMPING does not use omega parameter!");
      AddDampingCoefficients(1.0, fr, fbr);
      str = "C";
      break;
    case OperatorType::EXTRA:
      AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);
      str = "A2";
      break;
    default:
      MFEM_ABORT("Invalid GetSystemMatrix matrix type!");
  }
  std::unique_ptr<mfem::HypreParMatrix> hAr, hAi;
  bool has_real = false, has_imag = false;
  if (!dfr.empty() || !fr.empty() || !dfbr.empty() || !fbr.empty())
  {
    has_real = true;
    mfem::ParBilinearForm a(&nd_fespaces.GetFinestFESpace());
    AddIntegrators(a, dfr, fr, dfbr, fbr);
    // a.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    a.Assemble(skip_zeros);
    a.Finalize(skip_zeros);
    hAr.reset(a.ParallelAssemble());
    hAr->EliminateBC(dbc_tdof_list, ess_diag);
  }
  if (!dfi.empty() || !fi.empty() || !dfbi.empty() || !fbi.empty())
  {
    has_imag = true;
    mfem::ParBilinearForm a(&nd_fespaces.GetFinestFESpace());
    AddIntegrators(a, dfi, fi, dfbi, fbi);
    // a.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    a.Assemble(skip_zeros);
    a.Finalize(skip_zeros);
    hAi.reset(a.ParallelAssemble());
    hAi->EliminateBC(dbc_tdof_list, mfem::Operator::DiagonalPolicy::DIAG_ZERO);
  }
  if (!has_real && !has_imag)
  {
    return {};
  }
  auto A = std::make_unique<petsc::PetscShellMatrix>(
      nd_fespaces.GetFinestFESpace().GetComm(), std::move(hAr), std::move(hAi));
  if (!has_imag)
  {
    A->SetRealSymmetric();
  }
  else
  {
    A->SetSymmetric();
  }

  // Print some information.
  PrintHeader();
  if (print)
  {
    if (has_real && has_imag)
    {
      Mpi::Print(" Re{{{}}}: NNZ = {:d}, norm = {:e}\n Im{{{}}}: NNZ = {:d}, norm = {:e}\n",
                 str, A->NNZReal(), A->NormFReal(), str, A->NNZImag(), A->NormFImag());
    }
    else
    {
      Mpi::Print(" {}: NNZ = {:d}, norm = {:e}\n", str,
                 has_real ? A->NNZReal() : A->NNZImag(),
                 has_real ? A->NormFReal() : A->NormFImag());
    }
  }
  return A;
}

std::unique_ptr<mfem::Operator>
SpaceOperator::GetSystemMatrix(SpaceOperator::OperatorType type, double omega,
                               mfem::Operator::DiagonalPolicy ess_diag, bool print)
{
  // Construct the frequency-dependent complex linear system matrix:
  //                 A = K + iω C - ω² (Mr + i Mi) + A2(ω)
  // or any subset of its terms. For output as a HypreParMatrix, only some of
  // the terms are available.
  MFEM_VERIFY(omega == 0.0,
              "GetSystemMatrix for HypreParMatrix does not use omega parameter!");
  const int sdim = nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension();
  SumMatrixCoefficient df(sdim), f(sdim), fb(sdim);
  SumCoefficient dfb;
  std::string str;
  switch (type)
  {
    case OperatorType::STIFFNESS:
      AddStiffnessCoefficients(1.0, df, f, fb);
      str = "K";
      break;
    case OperatorType::MASS:
      AddRealMassCoefficients(1.0, false, f, fb);
      str = "M";
      break;
    case OperatorType::DAMPING:
      AddDampingCoefficients(1.0, f, fb);
      str = "C";
      break;
    case OperatorType::COMPLETE:
    case OperatorType::EXTRA:
    default:
      MFEM_ABORT("Invalid GetSystemMatrix matrix type for HypreParMatrix output!");
  }
  if (df.empty() && f.empty() && fb.empty())
  {
    return {};
  }
  mfem::ParBilinearForm a(&nd_fespaces.GetFinestFESpace());
  AddIntegrators(a, df, f, dfb, fb);
  // a.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  a.Assemble(skip_zeros);
  a.Finalize(skip_zeros);
  std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
  A->EliminateBC(dbc_tdof_list, ess_diag);

  // Print some information.
  PrintHeader();
  if (print)
  {
    Mpi::Print(" {}: NNZ = {:d}, norm = {:e}\n", str, A->NNZ(),
               hypre_ParCSRMatrixFnorm(*A));
  }
  return A;
}

void SpaceOperator::GetPreconditionerInternal(
    const std::function<void(SumMatrixCoefficient &, SumMatrixCoefficient &,
                             SumCoefficient &, SumMatrixCoefficient &)> &AddCoefficients,
    std::vector<std::unique_ptr<mfem::Operator>> &B,
    std::vector<std::unique_ptr<mfem::Operator>> &AuxB, bool print)
{
  // Construct the real, optionally SPD matrix for frequency or time domain preconditioning
  // (Mr > 0, Mi < 0):
  //              B =    K +  ω C + ω² (-/+ Mr - Mi) , or
  //              B = a0 K + a1 C +         Mr .
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == nd_fespaces.GetNumLevels(),
              "Multigrid heirarchy mismatch for auxiliary space preconditioning!");
  for (int s = 0; s < 2; s++)
  {
    auto &B_ = (s == 0) ? B : AuxB;
    B_.clear();
    B_.reserve(nd_fespaces.GetNumLevels());
    for (int l = 0; l < nd_fespaces.GetNumLevels(); l++)
    {
      auto &fespace_l =
          (s == 0) ? nd_fespaces.GetFESpaceAtLevel(l) : h1_fespaces.GetFESpaceAtLevel(l);
      mfem::Array<int> dbc_tdof_list_l;
      if (dbc_marker.Size() > 0)
      {
        fespace_l.GetEssentialTrueDofs(dbc_marker, dbc_tdof_list_l);
      }

      const int sdim = nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension();
      SumMatrixCoefficient df(sdim), f(sdim), fb(sdim);
      SumCoefficient dfb;
      AddCoefficients(df, f, dfb, fb);
      mfem::ParBilinearForm b(&fespace_l);
      if (s == 1)
      {
        // H1 auxiliary space matrix Gᵀ B G.
        AddAuxIntegrators(b, f, fb);
      }
      else
      {
        AddIntegrators(b, df, f, dfb, fb);
      }
      // b.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
      b.Assemble(skip_zeros);
      b.Finalize(skip_zeros);
      std::unique_ptr<mfem::HypreParMatrix> hB;
      if (pc_lor)
      {
        // After we construct the LOR discretization we can extract the LOR matrix and the
        // original bilinear form and LOR discretization are no longer needed.
        mfem::ParLORDiscretization lor(b, dbc_tdof_list_l);
        hB = std::make_unique<mfem::HypreParMatrix>(lor.GetAssembledMatrix());
      }
      else
      {
        hB.reset(b.ParallelAssemble());
      }
      hB->EliminateBC(dbc_tdof_list_l, mfem::Operator::DiagonalPolicy::DIAG_ONE);

      // Print some information.
      PrintHeader();
      if (s == 0 && print)
      {
        std::string str = "";
        if (pc_gmg && pc_lor)
        {
          str = fmt::format(" (Level {:d}, {:d} unknowns, LOR)", l,
                            fespace_l.GlobalTrueVSize());
        }
        else if (pc_gmg)
        {
          str = fmt::format(" (Level {:d}, {:d} unknowns)", l, fespace_l.GlobalTrueVSize());
        }
        else if (pc_lor)
        {
          str = " (LOR)";
        }
        Mpi::Print(" B{}: NNZ = {:d}, norm = {:e}\n", str, hB->NNZ(),
                   hypre_ParCSRMatrixFnorm(*hB));
      }
      B_.push_back(std::move(hB));
    }
  }
}

void SpaceOperator::GetPreconditionerMatrix(
    double omega, std::vector<std::unique_ptr<mfem::Operator>> &B,
    std::vector<std::unique_ptr<mfem::Operator>> &AuxB, bool print)
{
  // Frequency domain preconditioner matrix.
  auto AddCoefficients = [this, omega](SumMatrixCoefficient &df, SumMatrixCoefficient &f,
                                       SumCoefficient &dfb, SumMatrixCoefficient &fb)
  {
    this->AddStiffnessCoefficients(1.0, df, f, fb);
    this->AddDampingCoefficients(omega, f, fb);
    this->AddRealMassCoefficients(pc_shifted ? omega * omega : -omega * omega, true, f, fb);
    this->AddExtraSystemBdrCoefficients(omega, dfb, dfb, fb, fb);
  };
  GetPreconditionerInternal(AddCoefficients, B, AuxB, print);
}

void SpaceOperator::GetPreconditionerMatrix(
    double a0, double a1, std::vector<std::unique_ptr<mfem::Operator>> &B,
    std::vector<std::unique_ptr<mfem::Operator>> &AuxB, bool print)
{
  // Time domain preconditioner matrix.
  auto AddCoefficients = [this, a0, a1](SumMatrixCoefficient &df, SumMatrixCoefficient &f,
                                        SumCoefficient &dfb, SumMatrixCoefficient &fb)
  {
    this->AddStiffnessCoefficients(a0, df, f, fb);
    this->AddDampingCoefficients(a1, f, fb);
    this->AddRealMassCoefficients(1.0, false, f, fb);
  };
  GetPreconditionerInternal(AddCoefficients, B, AuxB, print);
}

std::unique_ptr<mfem::Operator> SpaceOperator::GetNegCurlMatrix()
{
  mfem::ParDiscreteLinearOperator curl(&nd_fespaces.GetFinestFESpace(), &rt_fespace);
  curl.AddDomainInterpolator(new mfem::CurlInterpolator);
  // curl.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  curl.Assemble();
  curl.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> NegCurl(curl.ParallelAssemble());
  *NegCurl *= -1.0;
  return NegCurl;
}

std::unique_ptr<petsc::PetscParMatrix> SpaceOperator::GetNegCurlMatrixPetsc()
{
  return std::make_unique<petsc::PetscShellMatrix>(nd_fespaces.GetFinestFESpace().GetComm(),
                                                   GetNegCurlMatrix());
}

std::unique_ptr<mfem::Operator> SpaceOperator::GetGradMatrix()
{
  mfem::ParDiscreteLinearOperator grad(&h1_fespaces.GetFinestFESpace(),
                                       &nd_fespaces.GetFinestFESpace());
  grad.AddDomainInterpolator(new mfem::GradientInterpolator);
  // grad.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  grad.Assemble();
  grad.Finalize();
  return std::unique_ptr<mfem::HypreParMatrix>(grad.ParallelAssemble());
}

std::unique_ptr<petsc::PetscParMatrix> SpaceOperator::GetGradMatrixPetsc()
{
  return std::make_unique<petsc::PetscShellMatrix>(nd_fespaces.GetFinestFESpace().GetComm(),
                                                   GetGradMatrix());
}

void SpaceOperator::AddStiffnessCoefficients(double coef, SumMatrixCoefficient &df,
                                             SumMatrixCoefficient &f,
                                             SumMatrixCoefficient &fb)
{
  // Contribution for curl-curl term.
  df.AddCoefficient(
      std::make_unique<MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY>>(
          mat_op, coef));

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    f.AddCoefficient(
        std::make_unique<
            MaterialPropertyCoefficient<MaterialPropertyType::INV_LONDON_DEPTH>>(mat_op,
                                                                                 coef),
        mat_op.GetLondonDepthMarker());
  }

  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coef, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coef, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coef, bool abs_coef,
                                            SumMatrixCoefficient &f,
                                            SumMatrixCoefficient &fb)
{
  if (abs_coef)
  {
    f.AddCoefficient(std::make_unique<
                     MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_ABS>>(
        mat_op, coef));
  }
  else
  {
    f.AddCoefficient(std::make_unique<
                     MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL>>(
        mat_op, coef));
  }

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
    f.AddCoefficient(
        std::make_unique<
            MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_IMAG>>(mat_op,
                                                                                  coef),
        mat_op.GetLossTangentMarker());
  }
}

void SpaceOperator::AddDampingCoefficients(double coef, SumMatrixCoefficient &f,
                                           SumMatrixCoefficient &fb)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    f.AddCoefficient(
        std::make_unique<MaterialPropertyCoefficient<MaterialPropertyType::CONDUCTIVITY>>(
            mat_op, coef),
        mat_op.GetConductivityMarker());
  }

  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coef, fb);
  surf_z_op.AddDampingBdrCoefficients(coef, fb);
  lumped_port_op.AddDampingBdrCoefficients(coef, fb);
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

bool SpaceOperator::GetTimeDomainExcitationVector(mfem::Vector &RHS)
{
  return GetExcitationVector1Internal(RHS);
}

bool SpaceOperator::GetFreqDomainExcitationVector(double omega, petsc::PetscParVector &RHS)
{
  mfem::Vector hRHSr, hRHSi;
  bool nnz1 = GetExcitationVector1Internal(hRHSr);
  if (nnz1)
  {
    RHS.SetFromVector(hRHSr);  // Sets into real part
    RHS.Scale(1i * omega);
  }
  else
  {
    RHS.SetZero();
  }
  bool nnz2 = GetExcitationVector2Internal(omega, hRHSr, hRHSi);
  if (nnz2)
  {
    petsc::PetscParVector RHS2(RHS.GetComm(), hRHSr, hRHSi);
    RHS.AXPY(1.0, RHS2);
  }
  return nnz1 || nnz2;
}

bool SpaceOperator::GetFreqDomainExcitationVector1(petsc::PetscParVector &RHS1)
{
  // Assemble the frequency domain excitation term, including only the contributions from
  // lumped ports and surface currents, which is purely imaginary with linear frequency
  // dependence (coefficient iω, it is accounted for later).
  mfem::Vector hRHS1;
  bool nnz = GetExcitationVector1Internal(hRHS1);
  RHS1.SetFromVector(hRHS1);  // Sets into real part
  return nnz;
}

bool SpaceOperator::GetFreqDomainExcitationVector2(double omega,
                                                   petsc::PetscParVector &RHS2)
{
  mfem::Vector hRHS2r, hRHS2i;
  bool nnz = GetExcitationVector2Internal(omega, hRHS2r, hRHS2i);
  RHS2.SetFromVectors(hRHS2r, hRHS2i);
  return nnz;
}

bool SpaceOperator::GetExcitationVector1Internal(mfem::Vector &RHS)
{
  // Assemble the time domain excitation -g'(t) J or -iω J. The g'(t) factor is not
  // accounted for here, it is accounted for in the time integration later. Likewise, the
  // coefficient iω, is accounted for later).
  SumVectorCoefficient fb(nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(fb);
  surf_j_op.AddExcitationBdrCoefficients(fb);
  RHS.SetSize(nd_fespaces.GetFinestFESpace().GetTrueVSize());
  RHS = 0.0;
  if (fb.empty())
  {
    return false;
  }
  mfem::ParLinearForm rhs(&nd_fespaces.GetFinestFESpace());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(true);
  rhs.Assemble();
  rhs.ParallelAssemble(RHS);
  RHS.SetSubVector(dbc_tdof_list, 0.0);
  return true;
}

bool SpaceOperator::GetExcitationVector2Internal(double omega, mfem::Vector &RHSr,
                                                 mfem::Vector &RHSi)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at the
  // specified frequency.
  SumVectorCoefficient fbr(nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension()),
      fbi(nd_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(omega, fbr, fbi);
  RHSr.SetSize(nd_fespaces.GetFinestFESpace().GetTrueVSize());
  RHSi.SetSize(nd_fespaces.GetFinestFESpace().GetTrueVSize());
  RHSr = 0.0;
  RHSi = 0.0;
  if (fbr.empty() && fbi.empty())
  {
    return false;
  }
  mfem::ParLinearForm rhsr(&nd_fespaces.GetFinestFESpace());
  mfem::ParLinearForm rhsi(&nd_fespaces.GetFinestFESpace());
  rhsr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr));
  rhsi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
  rhsr.UseFastAssembly(true);
  rhsi.UseFastAssembly(true);
  rhsr.Assemble();
  rhsi.Assemble();
  rhsr.ParallelAssemble(RHSr);
  rhsi.ParallelAssemble(RHSi);
  RHSr.SetSubVector(dbc_tdof_list, 0.0);
  RHSi.SetSubVector(dbc_tdof_list, 0.0);
  return true;
}

}  // namespace palace
