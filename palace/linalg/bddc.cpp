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
  std::cout << "In BDDC SetOperator\n";
  // See all BDDC options here: https://petsc.org/main/src/ksp/pc/impls/bddc/bddc.c.html
  PetscOptionsSetValue(NULL, "-pc_bddc_check_level", "2"); // int, Verbosity/Checks level, setting to 0 reduces cpu time and memory usage
  //PetscOptionsSetValue(NULL, "-pc_bddc_interface_ext_type", "LUMP"); // enum, Use DIRICHLET or LUMP to extend interface corrections to interior. LUMP leads to error
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_local_mat_graph", "1"); // bool, Use or not adjacency graph of local mat for interface analysis. No noticeable effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_local_mat_graph_square", "1"); // int, Square adjacency graph of local mat for interface analysis. No noticeable effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_graph_maxcount", "1"); // Maximum number of shared subdomains for a connected component
  // 0 or 1 converges well but uses a lot of memory. Not setting or using 2+ leads to much smaller coarse problem but worse convergence
  PetscOptionsSetValue(NULL, "-pc_bddc_corner_selection", "1"); // Activates face-based corner selection. No effect?
  PetscOptionsSetValue(NULL, "-pc_bddc_use_vertices", "1"); // setting to 0 leads to coarse size = 0 if maxcount <2
  PetscOptionsSetValue(NULL, "-pc_bddc_use_edges", "1"); // No effect
  PetscOptionsSetValue(NULL, "-pc_bddc_use_faces", "1"); // Has an effect when graph_maxcount >=2!!cc
  //PetscOptionsSetValue(NULL, "-pc_bddc_vertex_size", "4"); // Connected components smaller or equal to vertex size will be considered as primal vertices
  // 0 -> error, 1 -> same as not setting, 2 -> same as 1, 10 -> same as 1
  // TRY IT AGAIN WITH GRAPH_MAXCOUNT > 1??!!! still no noticeable effect?
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_nnsp", "1"); // No effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_change_of_basis", "1");  // Internal change of basis on local edge nodes. No effect?
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_change_on_faces", "1");  // Internal change of basis on local face nodes. No effect?
  //PetscOptionsSetValue(NULL, "-pc_bddc_switch_static", "1"); // Switch on static condensation ops around the interface preconditioner. Leads to MPI error
  //PetscOptionsSetValue(NULL, "-pc_bddc_coarse_eqs_per_proc", "10000"); // number of cores used for the coarse solve-> ncores = coarse_size/this
  //PetscOptionsSetValue(NULL, "-pc_bddc_coarsening_ratio","4"); // Coarsening ratio, only matters for multilevel
  PetscOptionsSetValue(NULL, "-pc_bddc_levels", "0"); // Maximum number of levels. 0 disables multilevel. Multilevel requires local coordinates
  //PetscOptionsSetValue(NULL, "-pc_bddc_coarse_eqs_limit", "100"); // maximum number of equations on coarsest grid. No effect...?
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_coarse_estimates", "1"); // Use estimated eigenvalues for coarse problem. No effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_deluxe_scaling", "1"); // Use deluxe scaling for BDDC. No effect when maxcount=1, error when maxcount>1
  //PetscOptionsSetValue(NULL, "-pc_bddc_schur_rebuild", "1"); // Whether or not the interface graph for Schur principal minors has to be rebuilt. No effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_schur_layers", "-1"); // Number of dofs' layers for the computation of principal minors (-1 uses all). No effect?
  // but maybe check again after the testing schur  exact?
  //PetscOptionsSetValue(NULL, "-pc_bddc_schur_exact", "1"); // Whether or not to use the exact Schur complement instead of the reduced one. No effect?
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_deluxe_zerorows", "1"); // Zero rows and columns of deluxe operators associated with primal dofs
  //PetscOptionsSetValue(NULL, "-pc_bddc_use_deluxe_singlemat", "1"); // Collapse deluxe operators
  //PetscOptionsSetValue(NULL, "-pc_bddc_adaptive_threshold", "1.0"); // No effect when maxcount=1, error when maxcount>1
  //PetscOptionsSetValue(NULL, "-pc_bddc_adaptive_nmin", "2");
  //PetscOptionsSetValue(NULL, "-pc_bddc_adaptive_nmax", "20");
  PetscOptionsSetValue(NULL, "-pc_bddc_symmetric", "1"); // Symmetric computation of primal basis functions. Setting to 0 uses more memory?
  //PetscOptionsSetValue(NULL, "-pc_bddc_coarse_adj", "1"); // Number of processors where to map the coarse adjacency list. No effect?
  //PetscOptionsSetValue(NULL, "-pc_bddc_detect_disconnected", "1"); // Detects disconnected subdomains. No effect
  //PetscOptionsSetValue(NULL, "-pc_bddc_detect_disconnected_filter", "1"); // Filters out small entries in the local matrix when detecting disconnected subdomains
  //PetscOptionsSetValue(NULL, "-pc_bddc_benign_trick", "0");
  //PetscOptionsSetValue(NULL, "-pc_bddc_benign_change", "0");
  //PetscOptionsSetValue(NULL, "-pc_bddc_benign_compute_correction", "0");
  // PetscOptionsSetValue(NULL, "-pc_bddc_nonetflux", "1"); // Automatic computation of no-net-flux quadrature weights. No effect
  // PetscOptionsSetValue(NULL, "-pc_bddc_eliminate_dirichlet", "1"); // Whether or not we want to eliminate dirichlet dofs during presolve. No effect

  // Linear solvers
  PetscOptionsSetValue(NULL, "-pc_bddc_neumann_pc_type", "lu");
  PetscOptionsSetValue(NULL, "-pc_bddc_neumann_pc_factor_mat_solver_type", "mumps");
  PetscOptionsSetValue(NULL, "-pc_bddc_dirichlet_pc_type", "lu");
  PetscOptionsSetValue(NULL, "-pc_bddc_dirichlet_pc_factor_mat_solver_type", "mumps");
  //PetscOptionsSetValue(NULL, "-pc_bddc_coarse_ksp_type", "preonly");
  PetscOptionsSetValue(NULL, "-pc_bddc_coarse_pc_type", "cholesky");
  PetscOptionsSetValue(NULL, "-pc_bddc_coarse_pc_factor_mat_solver_type", "mumps");

    const mfem::PetscParMatrix *petsc_mat = dynamic_cast<const mfem::PetscParMatrix*>(&op);
    const mfem::HypreParMatrix *hypre_mat = dynamic_cast<const mfem::HypreParMatrix*>(&op);

    if (petsc_mat) {
      std::cout << "it's a PetscParMatrix!\n";
    }
    else if (hypre_mat) {
      std::cout << "it's a HypreParMatrix!\n";
    }
    else {
        std::cout << "it's neither petsc nor hypre!?\n";
    }

  /*
  // Convert HypreParMatrix (op) to PETSc MATAIJ first (not MATIS)
  auto temp = std::make_unique<mfem::PetscParMatrix>(comm, &op, mfem::Operator::PETSC_MATAIJ);
  // Now convert MATAIJ to MATIS
  Mat matAIJ = temp->GetMat();
  Mat matIS;
  MatConvert(matAIJ, MATIS, MAT_INITIAL_MATRIX, &matIS);
  auto pA = std::make_unique<mfem::PetscParMatrix>(matIS, mfem::Operator::PETSC_MATIS);
*/
  auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATIS); // op -> PetscParMatrix. Leads to 0 connected components

  // Set BDDC opts
  mfem::PetscBDDCSolverParams opts;
  opts.SetSpace(&fespace.Get());
  if (complex_coarse)
  {
    MFEM_ABORT("BDDC solver does not support complex coarse matrix!");
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
  //opts.SetNatBdrDofs(&ess_tdof_list);
  solver = std::make_unique<mfem::PetscBDDCSolver>(*pA, opts);


  /*
  // try ASM?
  // gmres, restart100, asm, asm_overlap=2, lu. converges but takes many (~250) iterations coarse solve, increasing overlap helps
  // but the performance varies as np is changed...
PetscOptionsSetValue(NULL, "-ksp_type", "gmres"); //preonly does not convegre
PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "100");  // Larger restart
PetscOptionsSetValue(NULL, "-pc_type", "asm");
PetscOptionsSetValue(NULL, "-pc_asm_type", "restrict");//restrict is the default
//PetscOptionsSetValue(NULL, "-pc_asm_type", "basic");//restrict is the default
PetscOptionsSetValue(NULL, "-pc_asm_overlap", "10");  // Increase overlap?
PetscOptionsSetValue(NULL, "-sub_pc_type", "lu");
PetscOptionsSetValue(NULL, "-sub_pc_factor_mat_solver_type", "mumps");
PetscOptionsSetValue(NULL, "-ksp_monitor", "");
  auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATAIJ); // op -> PetscParMatrix
  solver = std::make_unique<mfem::PetscLinearSolver>(comm);
  solver->SetOperator(*pA);
  */
  // Single transmon on 6 cores.
  // Using overlap 1 does not converge.
  // 2 does (11 outer iterations, ~80 inner gmres its. Total 31 secs
  // 3 leads to same 11 outer its, ~50 inner gmres its. Total 28 secs
  // 7 leads to same 11 outer its, ~17 inner gmres its. Total 25 secs
  // 10 leads to same 11 outer its, ~11 inner gmres its. Total 26 secs
  // Need to try it on higher number of cores. varying overlap
}

}  // namespace palace

//#endif