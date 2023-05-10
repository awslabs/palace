// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ams.hpp"

#include "linalg/hypre.hpp"

namespace palace
{

HypreAmsSolver::HypreAmsSolver(mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace *h1_fespace, int cycle_it,
                               int smooth_it, int agg_coarsen, bool vector_interp,
                               bool op_singular, int print_lvl)
  : mfem::HypreSolver(),
    // From the Hypre docs for AMS: cycles 1, 5, 8, 11, 13 are fastest, 7 yields fewest its
    // (MFEM default is 13). 14 is similar to 11/13 but is cheaper in that is uses additive
    // scalar Pi-space corrections.
    cycle_type(vector_interp ? 1 : 14),
    // Control levels of aggressive coarsening based on problem type: SPD/ semi-definite
    // curl-curl operators are easier than indefinite frequency domain problems. When used
    // as the coarse solver of geometric multigrid, always do only a single V-cycle.
    sdim(nd_fespace.GetParMesh()->SpaceDimension()), ams_it(cycle_it),
    ams_smooth_it(smooth_it),
    // Use no aggressive coarsening for frequency domain problems when the preconditioner
    // matrix is not SPD.
    agg_levels(agg_coarsen), ams_singular(op_singular),
    print((print_lvl > 1) ? print_lvl - 1 : 0)
{
  // From MFEM: The AMS preconditioner may sometimes require inverting singular matrices
  // with BoomerAMG, which are handled correctly in hypre's Solve method, but can produce
  // hypre errors in the Setup (specifically in the l1 row norm computation). See the
  // documentation of MFEM's SetErrorMode() for more details.
  error_mode = IGNORE_HYPRE_ERRORS;

  // Set up the AMS solver.
  Initialize();
  ConstructAuxiliaryMatrices(nd_fespace, h1_fespace);
}

HypreAmsSolver::~HypreAmsSolver()
{
  HYPRE_AMSDestroy(ams);
}

void HypreAmsSolver::Initialize()
{
  // Create the Hypre solver object.
  HYPRE_AMSCreate(&ams);
  HYPRE_AMSSetDimension(ams, sdim);
  HYPRE_AMSSetCycleType(ams, cycle_type);

  // Set additional AMS options.
  int coarsen_type = 10;   // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP
  double theta = 0.5;      // AMG strength parameter = 0.25 is 2D optimal (0.5-0.8 for 3D)
  int amg_relax_type = 8;  // 3 = GS, 6 = symm. GS, 8 = l1-symm. GS, 13 = l1-GS,
                           // 18 = l1-Jacobi, 16 = Chebyshev
  int interp_type = 6;     // 6 = Extended+i, 0 = Classical, 13 = FF1
  int Pmax = 4;            // Interpolation width
  int relax_type = 2;      // 2 = l1-SSOR, 4 = trunc. l1-SSOR,
                           // 1 = l1-Jacobi, 16 = Chebyshev
  int relax_sweeps = ams_smooth_it;
  double weight = 1.0;
  double omega = 1.0;

  HYPRE_AMSSetSmoothingOptions(ams, relax_type, relax_sweeps, weight, omega);
  HYPRE_AMSSetAlphaAMGOptions(ams, coarsen_type, agg_levels, amg_relax_type, theta,
                              interp_type, Pmax);
  HYPRE_AMSSetBetaAMGOptions(ams, coarsen_type, agg_levels, amg_relax_type, theta,
                             interp_type, Pmax);

  // int coarse_relax_type = 8;  // Default, l1-symm. GS
  // HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams, coarse_relax_type);
  // HYPRE_AMSSetBetaAMGCoarseRelaxType(ams, coarse_relax_type);

  // Control printing and number of iterations for use as a preconditioner.
  HYPRE_AMSSetPrintLevel(ams, print);
  HYPRE_AMSSetMaxIter(ams, ams_it);
  HYPRE_AMSSetTol(ams, 0.0);
  // HYPRE_AMSSetTol(ams, 1.0e-16);  // Avoid issues with zero RHS

  // Set this option when solving a curl-curl problem with zero mass term.
  if (ams_singular)
  {
    HYPRE_AMSSetBetaPoissonMatrix(ams, nullptr);
  }
}

void HypreAmsSolver::ConstructAuxiliaryMatrices(mfem::ParFiniteElementSpace &nd_fespace,
                                                mfem::ParFiniteElementSpace *h1_fespace)
{
  // Set up the auxiliary spaces for the preconditioner. Mostly the same as MFEM's
  // HypreAMS:Init.
  mfem::FiniteElementCollection *h1_fec = nullptr;
  if (!h1_fespace)
  {

    // XX TODO REUSE FROM INPUT...

    h1_fec = new mfem::H1_FECollection(nd_fespace.GetMaxElementOrder(),
                                       nd_fespace.GetParMesh()->Dimension());
    h1_fespace = new mfem::ParFiniteElementSpace(nd_fespace.GetParMesh(), h1_fec);
  }
  {
    mfem::ParDiscreteLinearOperator grad(h1_fespace, &nd_fespace);
    grad.AddDomainInterpolator(new mfem::GradientInterpolator);
    // grad.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    grad.Assemble();
    grad.Finalize();
    G.reset(grad.ParallelAssemble());
  }
  {
    mfem::ParFiniteElementSpace h1d_fespace(h1_fespace->GetParMesh(), h1_fespace->FEColl(),
                                            sdim, mfem::Ordering::byVDIM);
    mfem::ParDiscreteLinearOperator id_ND(&h1d_fespace, &nd_fespace);
    id_ND.AddDomainInterpolator(new mfem::IdentityInterpolator);
    // id_ND.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
    id_ND.Assemble();
    id_ND.Finalize();
    if (cycle_type < 10)
    {
      Pi.reset(id_ND.ParallelAssemble());
    }
    else
    {
      mfem::Array2D<mfem::HypreParMatrix *> Pi_blocks;
      id_ND.GetParBlocks(Pi_blocks);
      Pix.reset(Pi_blocks(0, 0));
      if (sdim >= 2)
      {
        Piy.reset(Pi_blocks(0, 1));
      }
      if (sdim == 3)
      {
        Piz.reset(Pi_blocks(0, 2));
      }
    }
  }
  if (h1_fec)
  {
    delete h1_fespace;
    delete h1_fec;
  }
}

void HypreAmsSolver::SetOperator(const mfem::Operator &op)
{
  // When the operator changes, we need to rebuild the AMS solver but can use the unchanged
  // auxiliary space matrices.
  bool first = true;
  if (setup_called)
  {
    HYPRE_AMSDestroy(ams);
    Initialize();
    first = false;
  }
  setup_called = 0;

  // From mfem::HypreAMS: Update HypreSolver base class.
  height = op.Height();
  width = op.Width();
  auto *new_A = dynamic_cast<const mfem::HypreParMatrix *>(&op);
  MFEM_VERIFY(new_A, "AMS solver requires a HypreParMatrix operator!");
  A = const_cast<mfem::HypreParMatrix *>(new_A);
  delete X;
  delete B;
  B = X = nullptr;
  auxB.Delete();
  auxB.Reset();
  auxX.Delete();
  auxX.Reset();

  // Eliminate Dirichlet BCs in G, Pi matrices before RAP.
  if (first)
  {
    constexpr double tol = 1.0e-9;
    mfem::Array<int> dbc_rows;
    mfem::Vector diag(A->Height()), rowsums(A->Height());
    A->AssembleDiagonal(diag);
    hypre::hypreParCSRRowSums(*A, rowsums);
    for (int i = 0; i < A->Height(); i++)
    {
      if (std::abs(rowsums(i) - diag(i)) < tol * std::abs(diag(i)))
      {
        dbc_rows.Append(i);
      }
    }
    if (G)
    {
      G->EliminateRows(dbc_rows);
    }
    if (Pi)
    {
      Pi->EliminateRows(dbc_rows);
    }
    if (Pix)
    {
      Pix->EliminateRows(dbc_rows);
    }
    if (Piy)
    {
      Piy->EliminateRows(dbc_rows);
    }
    if (Piz)
    {
      Piz->EliminateRows(dbc_rows);
    }
  }

  // Set the discrete gradient and Nedelec interpolation matrices.
  HYPRE_ParCSRMatrix HY_Pi = (Pi) ? (HYPRE_ParCSRMatrix)*Pi : nullptr;
  HYPRE_ParCSRMatrix HY_Pix = (Pix) ? (HYPRE_ParCSRMatrix)*Pix : nullptr;
  HYPRE_ParCSRMatrix HY_Piy = (Piy) ? (HYPRE_ParCSRMatrix)*Piy : nullptr;
  HYPRE_ParCSRMatrix HY_Piz = (Piz) ? (HYPRE_ParCSRMatrix)*Piz : nullptr;
  HYPRE_AMSSetDiscreteGradient(ams, *G);
  HYPRE_AMSSetInterpolations(ams, HY_Pi, HY_Pix, HY_Piy, HY_Piz);
}

}  // namespace palace
