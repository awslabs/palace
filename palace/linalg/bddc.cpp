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
  : mfem::Solver(), fespace(fespace), comm(fespace.GetComm())
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
    std::cout << "hypre ensuremulttranspose/n";
    hypre_op->EnsureMultTranspose();  // Ensure fully assembled
  }

PetscOptionsSetValue(NULL, "-matis_convert_local_nest", "");
PetscOptionsSetValue(NULL, "-mat_is_convert_local_nest", "");
PetscOptionsSetValue(NULL, "-pc_bddc_check_level", "2");
//PetscOptionsSetValue(NULL, "-mat_is_symmetric", "");
//PetscOptionsSetValue(NULL, "-ksp_monitor", "");
PetscOptionsSetValue(NULL, "-pc_bddc_use_local_mat_graph", "0");
PetscOptionsSetValue(NULL, "-pc_bddc_detect_disconnected", "");
PetscOptionsSetValue(NULL, "-pc_bddc_use_edges", "1"); // ?

  auto pA = std::make_unique<mfem::PetscParMatrix>(comm, &op, Operator::PETSC_MATIS); // op -> PetscParMatrix
  mfem::PetscBDDCSolverParams opts;
  opts.SetSpace(&fespace.Get());
  opts.SetEssBdrDofs(&ess_tdof_list);
  // opts.SetComputeNetFlux(true); // leads to MFEM abort
  //opts.SetNatBdrDofs(&nat_tdof_list);
  solver = std::make_unique<mfem::PetscBDDCSolver>(*pA, opts);
  //std::cout << "bddc.cpp L48 solver->Width: " << solver->Width() << "\n";
}

}  // namespace palace

//#endif