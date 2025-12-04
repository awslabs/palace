// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bddc.hpp"

#include <petsc.h>
#include "fem/fespace.hpp"
#include "utils/geodata.hpp"
//#if defined(PALACE_WITH_PETSC)

namespace palace
{

BDDCSolver::BDDCSolver(const IoData &iodata, FiniteElementSpace &fespace, int print)
  : mfem::Solver(fespace.Get().GetTrueVSize()), fespace(fespace), comm(fespace.GetComm()), complex_coarse(iodata.solver.linear.complex_coarse_solve)
{
  const char *petscrc_file = "";
  mfem::MFEMInitializePetsc(NULL, NULL, petscrc_file, NULL);
  // anything else??

  if (fespace.GetParMesh().bdr_attributes.Size())
  {
    mfem::Array<int> ess_bdr(fespace.GetParMesh().bdr_attributes.Max());
    // only set PEC as essential?
    ess_bdr = 0;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      std::cout << "Setting ess bdr for PEC attribute: " << attr << "\n";
      ess_bdr[attr] = 1;
    }
    fespace.Get().GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // set everything else as nat?!
    mfem::Array<int> nat_bdr(fespace.GetParMesh().bdr_attributes.Max());
    for (int i = 0; i < ess_bdr.Size(); i++)
    {
      nat_bdr[i] = 1 - ess_bdr[i];
    }
    fespace.Get().GetEssentialTrueDofs(nat_bdr, nat_tdof_list);
  }

  //mfem::Array<int> dbc_marker, dbc_attr;
  //dbc_attr.push_back(-1);
  //int bdr_attr_max = fespace.GetParMesh().bdr_attributes.Size()
  //                         ? fespace.GetParMesh().bdr_attributes.Max()
  //                         : 0;
  //std::cout << "BDDC bdr_attr_max: " << bdr_attr_max << "\n";
  //dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  //std::cout << "BDDC dbc_marker.size: " << dbc_marker.Size() << " min/max: " << dbc_marker.Min() << " " << dbc_marker.Max() << "\n";
  //fespace.Get().GetEssentialTrueDofs(dbc_marker, bdr_tdof_list);
}

void BDDCSolver::SetOperator(const Operator &op)
{
  //const auto *pA = dynamic_cast<const mfem::PetscParMatrix *>(&op); // op -> PetscParMatrix
  //const auto *hA = dynamic_cast<const mfem::HypreParMatrix *>(&op); // op -> HypreParMatrix
  //auto pA = std::make_unique<mfem::PetscParMatrix>(hA, Operator::PETSC_MATIS); // HypreParMatrix -> PetscParMatrix
  //auto *pA = const_cast<mfem::PetscParMatrix*>(dynamic_cast<const mfem::PetscParMatrix*>(&op));
  //auto pA2 = std::make_unique<mfem::PetscParMatrix>(pA, Operator::PETSC_MATIS);

  if (auto *hypre_op = dynamic_cast<const mfem::HypreParMatrix*>(&op)) {
    std::cout << "hypre ensuremulttranspose\n";
    hypre_op->EnsureMultTranspose();  // Ensure fully assembled
  }

//PetscOptionsSetValue(NULL, "-matis_convert_local_nest", "");
//PetscOptionsSetValue(NULL, "-mat_is_convert_local_nest", "");
PetscOptionsSetValue(NULL, "-pc_bddc_check_level", "2");
//PetscOptionsSetValue(NULL, "-mat_is_symmetric", "");
//PetscOptionsSetValue(NULL, "-ksp_monitor", "");
PetscOptionsSetValue(NULL, "-pc_bddc_use_local_mat_graph", "0");
PetscOptionsSetValue(NULL, "-pc_bddc_detect_disconnected", "");
PetscOptionsSetValue(NULL, "-pc_bddc_use_edges", "1"); // ?
PetscOptionsSetValue(NULL, "-pc_bddc_use_faces", "1"); // ?
PetscOptionsSetValue(NULL, "-pc_bddc_corner_selection", "1"); //
PetscOptionsSetValue(NULL, "-pc_bddc_graph_maxcount", "1"); // ?
//PetscOptionsSetValue(NULL, "-pc_bddc_benign_trick", "");
//PetscOptionsSetValue(NULL, "-pc_bddc_nonetflux", "");
PetscOptionsSetValue(NULL, "-pc_bddc_levels", "1"); // ? doesn't help
//PetscOptionsSetValue(NULL, "-pc_bddc_use_deluxe_scaling", "1"); // leads to error?
//PetscOptionsSetValue(NULL, "-pc_bddc_deluxe_zerorows", ""); //? doesn't do much?
//PetscOptionsSetValue(NULL, "-pc_bddc_adaptive_threshold", "2.0"); // leads to error
PetscOptionsSetValue(NULL, "-pc_bddc_neumann_pc_type", "lu");
PetscOptionsSetValue(NULL, "-pc_bddc_neumann_pc_factor_mat_solver_type", "mumps");
PetscOptionsSetValue(NULL, "-pc_bddc_dirichlet_pc_type", "lu");
PetscOptionsSetValue(NULL, "-pc_bddc_dirichlet_pc_factor_mat_solver_type", "mumps");
PetscOptionsSetValue(NULL, "-pc_bddc_coarse_pc_type", "cholesky");
PetscOptionsSetValue(NULL, "-pc_bddc_coarse_pc_factor_mat_solver_type", "mumps");


const auto *pfes = dynamic_cast<const mfem::ParFiniteElementSpace*>(&fespace.Get());
std::cout << " local DOFs: " << pfes->GetVSize() << std::endl;
std::cout << " true DOFs: " << pfes->GetTrueVSize() << std::endl;

// Check for shared DOFs
auto &gc = pfes->GroupComm();
//std::cout << " num groups: " << gc.PrintInfo(mfem::out) << std::endl;
gc.PrintInfo();

  //auto *pA = const_cast<mfem::PetscParMatrix*>(dynamic_cast<const mfem::PetscParMatrix *>(&op)); // op -> PetscParMatrix
  auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATIS); // op -> PetscParMatrix

  //std::cout << "pA height/width: " << pA->Height() << " " << pA->Width() << "\n";
  //std::cout << "ess_tdof_list.Size(): " << ess_tdof_list.Size() << " ess_tdof_list.Max(): " << ess_tdof_list.Max() << "\n";

  mfem::PetscBDDCSolverParams opts;
  opts.SetSpace(&fespace.Get());
  std::cout << "complex_coarse: " << complex_coarse << "\n";
  if (complex_coarse)
  {
    if (ess_tdof_list.Size() > 0)
    {
      // Block matrix has structure [Ar, -Ai; Ai, Ar]
      // Essential DOFs appear in both blocks
      int n = ess_tdof_list.Size();
      int offset = pA->Height() / 2;
      block_ess_tdof_list.SetSize(2 * n);
      for (int i = 0; i < n; i++)
      {
        block_ess_tdof_list[i] = ess_tdof_list[i];           // Block (0,0) and (0,1)
        block_ess_tdof_list[n + i] = ess_tdof_list[i] + offset;  // Block (1,0) and (1,1)
      }
    }
    opts.SetEssBdrDofs(&block_ess_tdof_list);
  }
  else
  {
    opts.SetEssBdrDofs(&ess_tdof_list);
  }
  // opts.SetComputeNetFlux(true); // leads to MFEM abort
  //opts.SetNatBdrDofs(&nat_tdof_list);
  solver = std::make_unique<mfem::PetscBDDCSolver>(*pA, opts);
  std::cout << "bddc.cpp L48 solver->Width: " << solver->Width() << "\n";

  /*
  // try ASM?
  // gmres, restart100, asm, asm_overlap=2, lu. converges but takes many (~250) iterations coarse solve, increasing overlap helps
  // but the performance varies as np is changed...
PetscOptionsSetValue(NULL, "-ksp_type", "gmres"); //preonly does not convegre
PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "100");  // Larger restart
PetscOptionsSetValue(NULL, "-pc_type", "asm");
PetscOptionsSetValue(NULL, "-pc_asm_type", "restrict");//restrict is the default
PetscOptionsSetValue(NULL, "-pc_asm_overlap", "2");  // Increase overlap?
PetscOptionsSetValue(NULL, "-sub_pc_type", "lu");
PetscOptionsSetValue(NULL, "-sub_pc_factor_mat_solver_type", "mumps");
PetscOptionsSetValue(NULL, "-ksp_monitor", "");
  auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATAIJ); // op -> PetscParMatrix
  solver2 = std::make_unique<mfem::PetscLinearSolver>(comm);
  solver2->SetOperator(*pA);
  */
}

}  // namespace palace

//#endif