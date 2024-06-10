// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ams.hpp"

#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/omp.hpp"

namespace palace
{

HypreAmsSolver::HypreAmsSolver(FiniteElementSpace &nd_fespace,
                               FiniteElementSpace &h1_fespace, int cycle_it, int smooth_it,
                               bool vector_interp, bool singular_op, bool agg_coarsen,
                               int print)
  : mfem::HypreSolver(),
    // From the Hypre docs for AMS: cycles 1, 5, 8, 11, 13 are fastest, 7 yields fewest its
    // (MFEM default is 13). 14 is similar to 11/13 but is cheaper in that is uses additive
    // scalar Pi-space corrections.
    cycle_type(vector_interp ? 1 : 14), space_dim(nd_fespace.SpaceDimension()),
    // When used as the coarse solver of geometric multigrid, always do only a single
    // V-cycle.
    ams_it(cycle_it), ams_smooth_it(smooth_it),
    // If we know the operator is singular (no mass matrix, for magnetostatic problems),
    // internally the AMS solver will avoid G-space corrections.
    ams_singular(singular_op),
    // For positive (SPD) operators, we will use aggressive coarsening but not for frequency
    // domain problems when the preconditioner matrix is not SPD.
    agg_coarsen(agg_coarsen), print((print > 1) ? print - 1 : 0)
{
  // From MFEM: The AMS preconditioner may sometimes require inverting singular matrices
  // with BoomerAMG, which are handled correctly in Hypre's Solve method, but can produce
  // Hypre errors in the Setup (specifically in the row l1-norm computation). See the
  // documentation of MFEM's SetErrorMode() for more details.
  error_mode = IGNORE_HYPRE_ERRORS;

  // Set up the AMS solver.
  ConstructAuxiliaryMatrices(nd_fespace, h1_fespace);
  InitializeSolver();
}

HypreAmsSolver::~HypreAmsSolver()
{
  HYPRE_AMSDestroy(ams);
}

void HypreAmsSolver::ConstructAuxiliaryMatrices(FiniteElementSpace &nd_fespace,
                                                FiniteElementSpace &h1_fespace)
{
  // Set up the auxiliary space objects for the preconditioner. Mostly the same as MFEM's
  // HypreAMS:Init. Start with the discrete gradient matrix. We don't skip zeros for the
  // full assembly to accelerate things on GPU and since they shouldn't affect the sparsity
  // pattern of the parallel G^T A G matrix (computed by Hypre).
  const bool skip_zeros_interp = !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  {
    const auto *PtGP =
        dynamic_cast<const ParOperator *>(&nd_fespace.GetDiscreteInterpolator(h1_fespace));
    MFEM_VERIFY(
        PtGP,
        "HypreAmsSolver requires the discrete gradient matrix as a ParOperator operator!");
    G = &PtGP->ParallelAssemble(skip_zeros_interp);
  }

  // Vertex coordinates for the lowest order case, or Nedelec interpolation matrix or
  // matrices for order > 1. Expects that Mesh::SetVerticesFromNodes has been called at some
  // point to avoid calling GridFunction::GetNodalValues here.
  mfem::ParMesh &mesh = h1_fespace.GetParMesh();
  if (h1_fespace.GetMaxElementOrder() == 1)
  {
    mfem::ParGridFunction x_coord(&h1_fespace.Get()), y_coord(&h1_fespace.Get()),
        z_coord(&h1_fespace.Get());
    MFEM_VERIFY(x_coord.Size() == mesh.GetNV(),
                "Unexpected size for vertex coordinates in AMS setup!");
    PalacePragmaOmp(parallel for schedule(static))
    for (int i = 0; i < mesh.GetNV(); i++)
    {
      x_coord(i) = mesh.GetVertex(i)[0];
      if (space_dim > 1)
      {
        y_coord(i) = mesh.GetVertex(i)[1];
      }
      if (space_dim > 2)
      {
        z_coord(i) = mesh.GetVertex(i)[2];
      }
    }
    x.reset(x_coord.ParallelProject());
    x->HypreReadWrite();
    if (space_dim > 1)
    {
      y.reset(y_coord.ParallelProject());
      y->HypreReadWrite();
    }
    if (space_dim > 2)
    {
      z.reset(z_coord.ParallelProject());
      z->HypreReadWrite();
    }
  }
  else
  {
    // Fall back to MFEM legacy assembly for identity interpolator.
    FiniteElementSpace h1d_fespace(h1_fespace.GetMesh(), &h1_fespace.GetFEColl(), space_dim,
                                   mfem::Ordering::byVDIM);
    mfem::DiscreteLinearOperator pi(&h1d_fespace.Get(), &nd_fespace.Get());
    pi.AddDomainInterpolator(new mfem::IdentityInterpolator);
    pi.SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    pi.Assemble(skip_zeros_interp);
    pi.Finalize(skip_zeros_interp);
    ParOperator RAP_Pi(std::make_unique<hypre::HypreCSRMatrix>(pi.SpMat()), h1d_fespace,
                       nd_fespace, true);
    Pi = RAP_Pi.StealParallelAssemble();
    if (cycle_type >= 10)
    {
      // Get blocks of Pi corresponding to each component, and free Pi.
      mfem::Array2D<mfem::HypreParMatrix *> Pi_blocks(1, space_dim);
      Pi->GetBlocks(Pi_blocks, false, true);
      Pix.reset(Pi_blocks(0, 0));
      if (space_dim > 1)
      {
        Piy.reset(Pi_blocks(0, 1));
      }
      if (space_dim > 2)
      {
        Piz.reset(Pi_blocks(0, 2));
      }
      Pi.reset();
    }
  }
}

void HypreAmsSolver::InitializeSolver()
{
  // Create the Hypre solver object.
  HYPRE_AMSCreate(&ams);
  HYPRE_AMSSetDimension(ams, space_dim);
  HYPRE_AMSSetCycleType(ams, cycle_type);

  // Control printing and number of iterations for use as a preconditioner.
  HYPRE_AMSSetPrintLevel(ams, print);
  HYPRE_AMSSetMaxIter(ams, ams_it);
  // HYPRE_AMSSetTol(ams, 1.0e-16);  // Avoid issues with zero RHS

  // Set this option when solving a curl-curl problem with zero mass term.
  if (ams_singular)
  {
    HYPRE_AMSSetBetaPoissonMatrix(ams, nullptr);
  }

  // Set additional AMS options.
  int coarsen_type = 10;                     // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP
  int amg_agg_levels = agg_coarsen ? 1 : 0;  // Number of aggressive coarsening levels
  double theta = 0.5;      // AMG strength parameter = 0.25 is 2D optimal (0.5-0.8 for 3D)
  int amg_relax_type = 8;  // 3 = GS, 6 = symm. GS, 8 = l1-symm. GS, 13 = l1-GS,
                           // 18 = l1-Jacobi, 16 = Chebyshev
  int interp_type = 6;     // 6 = Extended+i, 0 = Classical, 13 = FF1
  int Pmax = 4;            // Interpolation width
  int relax_type = 2;      // 2 = l1-SSOR, 4 = trunc. l1-SSOR, 1 = l1-Jacobi, 16 = Chebyshev
  double weight = 1.0;
  double omega = 1.0;
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    // Modify options for GPU-supported features.
    coarsen_type = 8;
    amg_agg_levels = 0;
    amg_relax_type = 18;
    relax_type = 1;
  }

  HYPRE_AMSSetSmoothingOptions(ams, relax_type, ams_smooth_it, weight, omega);
  HYPRE_AMSSetAlphaAMGOptions(ams, coarsen_type, amg_agg_levels, amg_relax_type, theta,
                              interp_type, Pmax);
  HYPRE_AMSSetBetaAMGOptions(ams, coarsen_type, amg_agg_levels, amg_relax_type, theta,
                             interp_type, Pmax);

  // int coarse_relax_type = 8;  // Default, l1-symm. GS
  // HYPRE_AMSSetAlphaAMGCoarseRelaxType(ams, coarse_relax_type);
  // HYPRE_AMSSetBetaAMGCoarseRelaxType(ams, coarse_relax_type);

  // Set the discrete gradient matrix.
  HYPRE_AMSSetDiscreteGradient(ams, (HYPRE_ParCSRMatrix)*G);

  // Set the mesh vertex coordinates or Nedelec interpolation matrix or matrices.
  HYPRE_ParVector HY_X = (x) ? (HYPRE_ParVector)*x : nullptr;
  HYPRE_ParVector HY_Y = (y) ? (HYPRE_ParVector)*y : nullptr;
  HYPRE_ParVector HY_Z = (z) ? (HYPRE_ParVector)*z : nullptr;
  HYPRE_AMSSetCoordinateVectors(ams, HY_X, HY_Y, HY_Z);

  HYPRE_ParCSRMatrix HY_Pi = (Pi) ? (HYPRE_ParCSRMatrix)*Pi : nullptr;
  HYPRE_ParCSRMatrix HY_Pix = (Pix) ? (HYPRE_ParCSRMatrix)*Pix : nullptr;
  HYPRE_ParCSRMatrix HY_Piy = (Piy) ? (HYPRE_ParCSRMatrix)*Piy : nullptr;
  HYPRE_ParCSRMatrix HY_Piz = (Piz) ? (HYPRE_ParCSRMatrix)*Piz : nullptr;
  HYPRE_AMSSetInterpolations(ams, HY_Pi, HY_Pix, HY_Piy, HY_Piz);
}

void HypreAmsSolver::SetOperator(const Operator &op)
{
  // When the operator changes, we need to rebuild the AMS solver but can use the unchanged
  // auxiliary space matrices.
  if (A)
  {
    HYPRE_AMSDestroy(ams);
    InitializeSolver();
  }
  A = const_cast<mfem::HypreParMatrix *>(dynamic_cast<const mfem::HypreParMatrix *>(&op));
  MFEM_VERIFY(A, "HypreAmsSolver requires a HypreParMatrix operator!");
  height = A->Height();
  width = A->Width();

  // From mfem::HypreAMS: Update HypreSolver base class.
  setup_called = 0;
  delete X;
  delete B;
  B = X = nullptr;
  auxB.Delete();
  auxB.Reset();
  auxX.Delete();
  auxX.Reset();
}

}  // namespace palace
