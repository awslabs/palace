// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "curlcurl.hpp"

#include "fem/coefficient.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

CurlCurlMassSolver::CurlCurlMassSolver(const MaterialOperator &mat_op,
                                       const mfem::Array<int> &dbc_marker,
                                       mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
                                       mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                                       double tol, int max_it, int print)
  : mfem::Solver(nd_fespaces.GetFinestFESpace().GetTrueVSize())
{

  // XX TODO NEW ParOperator FRAMEWORK

  constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  constexpr MaterialPropertyType MatTypeEps = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
  MaterialPropertyCoefficient<MatTypeEps> epsilon_func(mat_op);
  MFEM_VERIFY(dbc_marker.Size() ==
                  nd_fespaces.GetFinestFESpace().GetParMesh()->bdr_attributes.Max(),
              "Invalid boundary marker for curl-curl solver!");
  for (int s = 0; s < 2; s++)
  {
    auto &A_ = (s == 0) ? A : AuxA;
    A_.reserve(nd_fespaces.GetNumLevels());
    for (int l = 0; l < nd_fespaces.GetNumLevels(); l++)
    {
      auto &fespace_l =
          (s == 0) ? nd_fespaces.GetFESpaceAtLevel(l) : h1_fespaces.GetFESpaceAtLevel(l);
      mfem::Array<int> dbc_tdof_list_l;
      fespace_l.GetEssentialTrueDofs(dbc_marker, dbc_tdof_list_l);

      mfem::ParBilinearForm a(&fespace_l);
      if (s == 1)
      {
        a.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(epsilon_func));
      }
      else
      {
        a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(muinv_func));
        a.AddDomainIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func));
      }
      // a.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
      a.Assemble();
      a.Finalize();
      mfem::HypreParMatrix *hA = a.ParallelAssemble();
      hA->EliminateBC(dbc_tdof_list_l, mfem::Operator::DiagonalPolicy::DIAG_ONE);
      A_.emplace_back(hA);
    }
  }

  // XX TODO VISIT

  // // The system matrix for the projection is real and SPD. For the coarse-level AMG
  // solve,
  // // we don't use an exact solve on the coarsest level.
  // auto ams = std::make_unique<HypreAmsSolver>(nd_fespaces.GetFESpaceAtLevel(0),
  //                                             &h1_fespaces.GetFESpaceAtLevel(0), 1, 1, 1,
  //                                             false, false, 0);
  // auto gmg = std::make_unique<GeometricMultigridSolver>(std::move(ams), dbc_marker,
  //                                                       nd_fespaces, &h1_fespaces, 1, 1,
  //                                                       2);
  // gmg->SetOperator(A, &AuxA);
  // pc = std::move(gmg);

  ksp = std::make_unique<mfem::CGSolver>(nd_fespaces.GetFinestFESpace().GetComm());
  ksp->SetRelTol(tol);
  ksp->SetMaxIter(max_it);
  ksp->SetPrintLevel(print);
  ksp->SetOperator(*A.back());
  ksp->SetPreconditioner(*pc);

  xr.SetSize(height);
  xi.SetSize(height);
  yr.SetSize(height);
  yi.SetSize(height);
}

}  // namespace palace
