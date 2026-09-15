// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bddc.hpp"

#if defined(PALACE_WITH_SLEPC)

#include <algorithm>
#include <cmath>
#include <petscksp.h>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "linalg/hypre.hpp"
#include "linalg/petsc.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Extract, for every local (L-vector) dof, the global true dof it maps to and the ±1
// prolongation sign. Returns false if some row is not a single ±1 entry, which is the
// nonconforming (hanging dof) case that MATIS cannot represent.
bool ExtractLocalToGlobal(const mfem::ParFiniteElementSpace &fespace,
                          std::vector<PetscInt> &l2g, std::vector<double> &sign)
{
  const mfem::HypreParMatrix *P =
      const_cast<mfem::ParFiniteElementSpace &>(fespace).Dof_TrueDof_Matrix();
  mfem::SparseMatrix diag, offd;
  HYPRE_BigInt *cmap = nullptr;
  P->GetDiag(diag);
  P->GetOffd(offd, cmap);
  const HYPRE_BigInt col_start = P->ColPart()[0];
  const int n = diag.Height();
  l2g.assign(n, -1);
  sign.assign(n, 0.0);
  for (int i = 0; i < n; i++)
  {
    int count = 0;
    for (int k = diag.GetI()[i]; k < diag.GetI()[i + 1]; k++)
    {
      if (diag.GetData()[k] == 0.0)
      {
        continue;
      }
      l2g[i] = static_cast<PetscInt>(col_start + diag.GetJ()[k]);
      sign[i] = diag.GetData()[k];
      count++;
    }
    for (int k = offd.GetI()[i]; k < offd.GetI()[i + 1]; k++)
    {
      if (offd.GetData()[k] == 0.0)
      {
        continue;
      }
      l2g[i] = static_cast<PetscInt>(cmap[offd.GetJ()[k]]);
      sign[i] = offd.GetData()[k];
      count++;
    }
    if (count != 1 || std::abs(std::abs(sign[i]) - 1.0) > 1.0e-12)
    {
      return false;
    }
  }
  return true;
}

// Convert an MFEM HypreParMatrix to a PETSc AIJ matrix, used for the discrete gradient.
Mat ConvertToPetscAIJ(const mfem::HypreParMatrix &A)
{
  mfem::SparseMatrix diag, offd;
  HYPRE_BigInt *cmap = nullptr;
  const_cast<mfem::HypreParMatrix &>(A).GetDiag(diag);
  const_cast<mfem::HypreParMatrix &>(A).GetOffd(offd, cmap);
  const HYPRE_BigInt row_start = A.RowPart()[0], col_start = A.ColPart()[0];

  Mat B;
  PalacePetscCall(MatCreate(A.GetComm(), &B));
  PalacePetscCall(MatSetType(B, MATAIJ));
  PalacePetscCall(MatSetSizes(B, A.Height(), A.Width(), A.M(), A.N()));
  std::vector<PetscInt> d_nnz(diag.Height()), o_nnz(diag.Height());
  for (int i = 0; i < diag.Height(); i++)
  {
    d_nnz[i] = diag.GetI()[i + 1] - diag.GetI()[i];
    o_nnz[i] = offd.GetI()[i + 1] - offd.GetI()[i];
  }
  PalacePetscCall(MatMPIAIJSetPreallocation(B, 0, d_nnz.data(), 0, o_nnz.data()));
  PalacePetscCall(MatSeqAIJSetPreallocation(B, 0, d_nnz.data()));
  for (int i = 0; i < diag.Height(); i++)
  {
    const PetscInt gi = static_cast<PetscInt>(row_start + i);
    for (int k = diag.GetI()[i]; k < diag.GetI()[i + 1]; k++)
    {
      const PetscInt gj = static_cast<PetscInt>(col_start + diag.GetJ()[k]);
      const PetscScalar v = diag.GetData()[k];
      PalacePetscCall(MatSetValues(B, 1, &gi, 1, &gj, &v, INSERT_VALUES));
    }
    for (int k = offd.GetI()[i]; k < offd.GetI()[i + 1]; k++)
    {
      const PetscInt gj = static_cast<PetscInt>(cmap[offd.GetJ()[k]]);
      const PetscScalar v = offd.GetData()[k];
      PalacePetscCall(MatSetValues(B, 1, &gi, 1, &gj, &v, INSERT_VALUES));
    }
  }
  PalacePetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
  PalacePetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  return B;
}

}  // namespace

struct BddcSolver::PetscData
{
  Mat A = nullptr;
  Mat A_loc = nullptr;
  ISLocalToGlobalMapping l2g_map = nullptr;
  PC pc = nullptr;
  Vec x = nullptr, y = nullptr;

  ~PetscData()
  {
    if (pc)
    {
      PCDestroy(&pc);
    }
    if (x)
    {
      VecDestroy(&x);
    }
    if (y)
    {
      VecDestroy(&y);
    }
    if (A)
    {
      MatDestroy(&A);
    }
    if (A_loc)
    {
      MatDestroy(&A_loc);
    }
    if (l2g_map)
    {
      ISLocalToGlobalMappingDestroy(&l2g_map);
    }
  }
};

BddcSolver::BddcSolver(const FiniteElementSpace &fespace,
                       const FiniteElementSpace *aux_fespace, int print)
  : Solver<Operator>(), fespace(fespace), aux_fespace(aux_fespace), print(print)
{
}

BddcSolver::~BddcSolver() = default;

void BddcSolver::SetOperator(const Operator &op)
{
  // BDDC needs the unassembled subdomain matrix, which is the local operator behind
  // Palace's parallel RAP operator.
  const auto *PtAP = dynamic_cast<const ParOperator *>(&op);
  MFEM_VERIFY(PtAP, "BDDC requires a ParOperator to access the subdomain matrix!");
  MFEM_VERIFY(!mfem::Device::Allows(mfem::Backend::DEVICE_MASK),
              "The BDDC preconditioner does not support GPU execution!");

  auto data = std::make_unique<PetscData>();
  MPI_Comm comm = fespace.GetComm();

  std::vector<PetscInt> l2g;
  std::vector<double> sign;
  int ok = ExtractLocalToGlobal(fespace.Get(), l2g, sign) ? 1 : 0;
  Mpi::GlobalMin(1, &ok, comm);
  MFEM_VERIFY(ok, "BDDC requires a conforming mesh: the prolongation operator must have a "
                  "single ±1 entry per local dof!");

  // BDDC's subdomain solves need explicit matrix entries, so assemble matrix-free
  // (libCEED) operators here.
  const Operator &A_loc_op = PtAP->LocalOperator();
  const auto *csr = dynamic_cast<const hypre::HypreCSRMatrix *>(&A_loc_op);
  std::unique_ptr<hypre::HypreCSRMatrix> csr_data;
  if (!csr)
  {
    const auto *ceed_op = dynamic_cast<const ceed::Operator *>(&A_loc_op);
    MFEM_VERIFY(ceed_op,
                "BDDC requires the local operator as a sparse or libCEED operator!");
    csr_data = BilinearForm::FullAssemble(*ceed_op, false, false);
    csr = csr_data.get();
  }

  // Mark essential (Dirichlet) dofs on local dofs. They are eliminated from the local
  // matrix below so the assembled MATIS matches the operator's DIAG_ONE elimination.
  const int n_loc = csr->Height();
  MFEM_VERIFY(static_cast<std::size_t>(n_loc) == l2g.size(),
              "Mismatch between the local operator size and the local dof count!");
  std::vector<char> is_dbc(n_loc, 0);
  if (const mfem::Array<int> *dbc_tdofs = PtAP->GetEssentialTrueDofs())
  {
    Vector marker_t(fespace.GetTrueVSize()), marker_l(n_loc);
    marker_t = 0.0;
    linalg::SetSubVector(marker_t, *dbc_tdofs, 1.0);
    fespace.GetProlongationMatrix()->Mult(marker_t, marker_l);
    const auto *m = marker_l.HostRead();
    for (int i = 0; i < n_loc; i++)
    {
      is_dbc[i] = (std::abs(m[i]) > 0.5);
    }
  }

  // Multiplicity of each local dof, so a unit diagonal can be split across the subdomains
  // sharing a constrained dof. |P|ᵀ1 counts per true dof; P maps the count back.
  std::vector<double> mult_loc(n_loc, 1.0);
  {
    Vector ones_l(n_loc), mult_t(fespace.GetTrueVSize()), back_l(n_loc);
    ones_l = 1.0;
    const auto *P =
        const_cast<mfem::ParFiniteElementSpace &>(fespace.Get()).Dof_TrueDof_Matrix();
    P->AbsMultTranspose(1.0, ones_l, 0.0, mult_t);
    P->Mult(mult_t, back_l);
    const auto *b = back_l.HostRead();
    for (int i = 0; i < n_loc; i++)
    {
      mult_loc[i] = std::max(1.0, std::abs(b[i]));
    }
  }

  // Local matrix with the prolongation signs absorbed: Ã = S A S.
  PalacePetscCall(MatCreate(PETSC_COMM_SELF, &data->A_loc));
  PalacePetscCall(MatSetType(data->A_loc, MATSEQAIJ));
  PalacePetscCall(MatSetSizes(data->A_loc, n_loc, n_loc, n_loc, n_loc));
  {
    std::vector<PetscInt> nnz(n_loc);
    for (int i = 0; i < n_loc; i++)
    {
      nnz[i] = static_cast<PetscInt>(csr->GetI()[i + 1] - csr->GetI()[i]);
    }
    PalacePetscCall(MatSeqAIJSetPreallocation(data->A_loc, 0, nnz.data()));
  }
  for (int i = 0; i < n_loc; i++)
  {
    const PetscInt row = i;
    if (is_dbc[i])
    {
      // Unit diagonal, split across the subdomains sharing this dof.
      const PetscScalar v = 1.0 / mult_loc[i];
      PalacePetscCall(MatSetValues(data->A_loc, 1, &row, 1, &row, &v, INSERT_VALUES));
      continue;
    }
    for (auto k = csr->GetI()[i]; k < csr->GetI()[i + 1]; k++)
    {
      const PetscInt j = static_cast<PetscInt>(csr->GetJ()[k]);
      if (is_dbc[j])
      {
        continue;  // drop coupling to constrained dofs
      }
      const PetscScalar v = csr->GetData()[k] * sign[i] * sign[j];
      PalacePetscCall(MatSetValues(data->A_loc, 1, &row, 1, &j, &v, INSERT_VALUES));
    }
  }
  PalacePetscCall(MatAssemblyBegin(data->A_loc, MAT_FINAL_ASSEMBLY));
  PalacePetscCall(MatAssemblyEnd(data->A_loc, MAT_FINAL_ASSEMBLY));

  // Assemble the MATIS operator over the subdomain partition.
  PalacePetscCall(ISLocalToGlobalMappingCreate(comm, 1, static_cast<PetscInt>(l2g.size()),
                                               l2g.data(), PETSC_COPY_VALUES,
                                               &data->l2g_map));
  const PetscInt n_owned = fespace.GetTrueVSize();
  const PetscInt n_global = fespace.GlobalTrueVSize();
  PalacePetscCall(MatCreateIS(comm, 1, n_owned, n_owned, n_global, n_global, data->l2g_map,
                              data->l2g_map, &data->A));
  PalacePetscCall(MatISSetLocalMat(data->A, data->A_loc));
  PalacePetscCall(MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY));
  PalacePetscCall(MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY));
  PalacePetscCall(MatSetOption(data->A, MAT_SYMMETRIC, PETSC_TRUE));

  PalacePetscCall(PCCreate(comm, &data->pc));
  PalacePetscCall(PCSetType(data->pc, PCBDDC));
  PalacePetscCall(PCSetOperators(data->pc, data->A, data->A));
  // Constraints are already in the local matrix; do not declare them to BDDC as well.

  // H(curl): the discrete gradient lets BDDC handle the curl-curl kernel.
  if (aux_fespace)
  {
    const Operator &G = fespace.GetDiscreteInterpolator(*aux_fespace);
    const auto *G_par = dynamic_cast<const ParOperator *>(&G);
    MFEM_VERIFY(G_par, "BDDC requires the discrete gradient as a ParOperator!");
    Mat G_petsc = ConvertToPetscAIJ(G_par->ParallelAssemble());
    const PetscInt order = fespace.Get().FEColl()->GetOrder();
    PalacePetscCall(
        PCBDDCSetDiscreteGradient(data->pc, G_petsc, order, 0, PETSC_TRUE, PETSC_FALSE));
    PalacePetscCall(MatDestroy(&G_petsc));
  }

  PalacePetscCall(PCSetFromOptions(data->pc));
  PalacePetscCall(PCSetUp(data->pc));
  if (print > 1)
  {
    PalacePetscCall(PCView(data->pc, PETSC_VIEWER_STDOUT_(comm)));
  }
  PalacePetscCall(MatCreateVecs(data->A, &data->x, &data->y));

  petsc = std::move(data);
  height = op.Height();
  width = op.Width();
}

void BddcSolver::Mult(const Vector &x, Vector &y) const
{
  MFEM_VERIFY(petsc, "BddcSolver::SetOperator must be called before Mult!");
  MFEM_ASSERT(x.Size() == height && y.Size() == width,
              "Incompatible dimensions for BddcSolver::Mult!");
  {
    PetscScalar *px;
    PalacePetscCall(VecGetArrayWrite(petsc->x, &px));
    const auto *hx = x.HostRead();
    for (int i = 0; i < x.Size(); i++)
    {
      px[i] = hx[i];
    }
    PalacePetscCall(VecRestoreArrayWrite(petsc->x, &px));
  }
  PalacePetscCall(PCApply(petsc->pc, petsc->x, petsc->y));
  {
    // PETSc is built with complex scalars for SLEPc; discard the zero imaginary part.
    const PetscScalar *py;
    PalacePetscCall(VecGetArrayRead(petsc->y, &py));
    auto *hy = y.HostWrite();
    for (int i = 0; i < y.Size(); i++)
    {
      hy[i] = PetscRealPart(py[i]);
    }
    PalacePetscCall(VecRestoreArrayRead(petsc->y, &py));
  }
}

}  // namespace palace

#endif
