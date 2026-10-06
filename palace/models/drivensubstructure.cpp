// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensubstructure.hpp"

#include <algorithm>
#include <utility>
#include "linalg/mumpsschur.hpp"
#include "linalg/rap.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

#if !defined(MFEM_USE_MUMPS)
class MumpsSchurSolver
{
};
#endif

namespace
{

// Local true DOFs touched by the elements with the given attributes.
std::vector<char> MarkTrueDofs(const mfem::ParFiniteElementSpace &fes,
                               const std::vector<int> &attrs)
{
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  mfem::Vector m(fes.GetVSize()), mt(fes.GetTrueVSize());
  m = 0.0;
  mfem::Array<int> vdofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (std::ranges::find(attrs, mesh.GetAttribute(e)) != attrs.end())
    {
      fes.GetElementVDofs(e, vdofs);
      for (int d : vdofs)
      {
        m(d >= 0 ? d : -1 - d) = 1.0;
      }
    }
  }
  fes.GetProlongationMatrix()->MultTranspose(m, mt);
  std::vector<char> mark(mt.Size());
  for (int i = 0; i < mt.Size(); i++)
  {
    mark[i] = (mt(i) > 0.0);
  }
  return mark;
}

// The assembled real (imag = false) or imaginary part of an operator, if any.
const mfem::HypreParMatrix *Part(const ComplexOperator *A, bool imag)
{
  const auto *part =
      A ? dynamic_cast<const ParOperator *>(imag ? A->Imag() : A->Real()) : nullptr;
  return part ? &part->ParallelAssemble() : nullptr;
}

// sum_k a_k X_k over the given (coefficient, matrix) terms (null matrices skipped).
std::unique_ptr<mfem::HypreParMatrix>
Sum(const std::vector<std::pair<double, const mfem::HypreParMatrix *>> &terms)
{
  std::unique_ptr<mfem::HypreParMatrix> sum;
  for (const auto &[a, X] : terms)
  {
    if (!X)
    {
      continue;
    }
    if (!sum)
    {
      sum = std::make_unique<mfem::HypreParMatrix>(*X);
      *sum *= a;
    }
    else
    {
      sum.reset(mfem::Add(1.0, *sum, a, *X));
    }
  }
  return sum;
}

}  // namespace

DrivenSubstructure::DrivenSubstructure(SpaceOperator &space_op,
                                       const std::vector<int> &region_attributes,
                                       const std::vector<int> &environment_attributes)
  : space_op(space_op), env_attrs(environment_attributes)
{
  const auto &fes = space_op.GetNDSpace().Get();
  const auto &mesh = *fes.GetParMesh();
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int a = mesh.GetAttribute(e);
    MFEM_VERIFY(std::ranges::find(region_attributes, a) != region_attributes.end() ||
                    std::ranges::find(env_attrs, a) != env_attrs.end(),
                "Domain attribute " << a
                                    << " is neither in the region nor in the "
                                       "environment!");
  }

  // Interface DOFs: on region and environment elements; environment interior: on
  // environment elements only (both without the Dirichlet DOFs).
  const std::vector<char> rm = MarkTrueDofs(fes, region_attributes),
                          em = MarkTrueDofs(fes, env_attrs);
  const int nt = fes.GetTrueVSize();
  std::vector<char> dbc(nt, 0);
  for (int d : space_op.GetNDDbcTDofLists().back())
  {
    dbc[d] = 1;
  }
  is_env_int.assign(nt, 0);
  is_gamma.assign(nt, 0);
  for (int i = 0; i < nt; i++)
  {
    is_gamma[i] = !dbc[i] && rm[i] && em[i];
    is_env_int[i] = !dbc[i] && em[i] && !rm[i];
    if (!is_gamma[i] && !is_env_int[i])
    {
      other.Append(i);
    }
  }

  // Interface order: by owner rank, then local index.
  std::vector<HYPRE_BigInt> mine;
  const HYPRE_BigInt tstart = fes.GetMyTDofOffset();
  for (int i = 0; i < nt; i++)
  {
    if (is_gamma[i])
    {
      mine.push_back(tstart + i);
    }
  }
  MPI_Comm comm = fes.GetComm();
  const int nranks = Mpi::Size(comm);
  int nloc = static_cast<int>(mine.size());
  std::vector<int> cnt(nranks), disp(nranks, 0);
  MPI_Allgather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
  for (int r = 1; r < nranks; r++)
  {
    disp[r] = disp[r - 1] + cnt[r - 1];
  }
  gamma_tdofs.resize(disp[nranks - 1] + cnt[nranks - 1]);
  MPI_Allgatherv(mine.data(), nloc, HYPRE_MPI_BIG_INT, gamma_tdofs.data(), cnt.data(),
                 disp.data(), HYPRE_MPI_BIG_INT, comm);

  // The frequency-independent environment operators.
  SpaceOperator::AssemblyRestriction environment(space_op, env_attrs);
  K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
}

DrivenSubstructure::~DrivenSubstructure() = default;

std::unique_ptr<mfem::HypreParMatrix> DrivenSubstructure::BlockOperator(double omega)
{
  // A_E(ω) = K + iω C - ω² M + A2(ω): real part Kr - ω Ci - ω² Mr + A2r, imaginary part
  // Ki + ω Cr - ω² Mi + A2i.
  std::unique_ptr<ComplexOperator> A2;
  {
    SpaceOperator::AssemblyRestriction environment(space_op, env_attrs);
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  }
  const double w2 = omega * omega;
  auto Ar = Sum({{1.0, Part(K.get(), false)},
                 {-omega, Part(C.get(), true)},
                 {-w2, Part(M.get(), false)},
                 {1.0, Part(A2.get(), false)}});
  auto Ai = Sum({{1.0, Part(K.get(), true)},
                 {omega, Part(C.get(), false)},
                 {-w2, Part(M.get(), true)},
                 {1.0, Part(A2.get(), true)}});
  MFEM_VERIFY(Ar, "Missing real part of the environment operator!");
  if (!Ai)
  {
    // A lossless environment: a zero imaginary part with the pattern of the real one.
    Ai = std::make_unique<mfem::HypreParMatrix>(*Ar);
    *Ai = 0.0;
  }
  // Pin the DOFs outside the environment interior and Γ, which leaves the Schur
  // complement unchanged.
  Ar->EliminateBC(other, Operator::DIAG_ONE);
  Ai->EliminateBC(other, Operator::DIAG_ZERO);
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  mfem::Array2D<double> coeffs(2, 2);
  blocks(0, 0) = Ar.get();
  blocks(0, 1) = Ai.get();
  blocks(1, 0) = Ai.get();
  blocks(1, 1) = Ar.get();
  coeffs(0, 0) = coeffs(0, 1) = coeffs(1, 0) = 1.0;
  coeffs(1, 1) = -1.0;
  return std::unique_ptr<mfem::HypreParMatrix>(
      mfem::HypreParMatrixFromBlocks(blocks, &coeffs));
}

void DrivenSubstructure::Condense(double omega)
{
#if !defined(MFEM_USE_MUMPS)
  MFEM_ABORT("Driven substructuring requires MUMPS!");
#else
  auto B = BlockOperator(omega);
  const int nG = InterfaceSize();
  if (!schur)
  {
    // Schur variables of the real form: the real, then the imaginary interface parts. Its
    // local rows are the local real rows, then the local imaginary rows.
    const auto &fes = space_op.GetNDSpace().Get();
    const HYPRE_BigInt tstart = fes.GetMyTDofOffset(), bstart = B->GetRowStarts()[0];
    const int nt = fes.GetTrueVSize();
    std::vector<HYPRE_BigInt> vars(2 * nG);
    for (int g = 0; g < nG; g++)
    {
      const HYPRE_BigInt i = gamma_tdofs[g] - tstart;  // local only on its owner
      vars[g] = vars[nG + g] = -1;
      if (i >= 0 && i < nt)
      {
        vars[g] = bstart + i;
        vars[nG + g] = bstart + nt + i;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, vars.data(), 2 * nG, HYPRE_MPI_BIG_INT, MPI_MAX,
                  fes.GetComm());
    schur = std::make_unique<MumpsSchurSolver>(*B, vars, 0.0, false, true);
  }
  else
  {
    schur->Refactor(*B);
  }

  // The Schur complement of [[Ar, Ai], [Ai, -Ar]] is [[Sr, Si], [Si, -Sr]].
  S.clear();
  if (Mpi::Root(space_op.GetComm()))
  {
    const auto &SB = schur->Schur();
    const std::size_t ld = 2 * static_cast<std::size_t>(nG);
    S.resize(static_cast<std::size_t>(nG) * nG);
    for (int b = 0; b < nG; b++)
    {
      for (int a = 0; a < nG; a++)
      {
        S[static_cast<std::size_t>(b) * nG + a] = {SB[b * ld + a], SB[(nG + b) * ld + a]};
      }
    }
  }
#endif
}

}  // namespace palace
