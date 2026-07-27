// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "blockprecond.hpp"

#include <mfem.hpp>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

void CopySubVector(const Vector &src, Vector &dst, int src_offset, int dst_offset, int size)
{
  const bool use_dev = src.UseDevice() || dst.UseDevice();
  const auto *sd = src.Read(use_dev);
  auto *dd = dst.Write(use_dev);
  mfem::forall_switch(use_dev, size, [=] MFEM_HOST_DEVICE(int i)
                      { dd[dst_offset + i] = sd[src_offset + i]; });
}

}  // namespace

template <>
void BlockDiagonalPreconditioner<Operator>::ExtractSubVector(const Vector &src, Vector &dst,
                                                             int offset, int size)
{
  CopySubVector(src, dst, offset, 0, size);
}

template <>
void BlockDiagonalPreconditioner<Operator>::InsertSubVector(const Vector &src, Vector &dst,
                                                            int offset, int size)
{
  CopySubVector(src, dst, 0, offset, size);
}

template <>
void BlockDiagonalPreconditioner<ComplexOperator>::ExtractSubVector(
    const ComplexVector &src, ComplexVector &dst, int offset, int size)
{
  CopySubVector(src.Real(), dst.Real(), offset, 0, size);
  CopySubVector(src.Imag(), dst.Imag(), offset, 0, size);
}

template <>
void BlockDiagonalPreconditioner<ComplexOperator>::InsertSubVector(const ComplexVector &src,
                                                                   ComplexVector &dst,
                                                                   int offset, int size)
{
  CopySubVector(src.Real(), dst.Real(), 0, offset, size);
  CopySubVector(src.Imag(), dst.Imag(), 0, offset, size);
}

template <typename OperType>
BlockDiagonalPreconditioner<OperType>::BlockDiagonalPreconditioner(
    int block0_size, std::unique_ptr<Solver<OperType>> &&pc0,
    std::unique_ptr<Solver<OperType>> &&pc1)
  : Solver<OperType>(), block0_size(block0_size), pc0(std::move(pc0)), pc1(std::move(pc1))
{
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::SetBlockOperators(const OperType &op0,
                                                              const OperType &op1)
{
  block0_size = op0.Height();
  pc0->SetOperator(op0);
  pc1->SetOperator(op1);
  this->height = op0.Height() + op1.Height();
  this->width = this->height;
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::SetOperator(const OperType &op)
{
  // No-op: use SetBlockOperators to set each block's operator independently.
  this->height = op.Height();
  this->width = op.Width();
}

template <typename OperType>
void BlockDiagonalPreconditioner<OperType>::Mult(const VecType &x, VecType &y) const
{
  const int n0 = block0_size;
  const int n1 = x.Size() - n0;

  // Lazy allocation of temporary vectors.
  if (x0.Size() != n0)
  {
    x0.SetSize(n0);
    y0.SetSize(n0);
    x0.UseDevice(true);
    y0.UseDevice(true);
  }
  if (x1.Size() != n1)
  {
    x1.SetSize(n1);
    y1.SetSize(n1);
    x1.UseDevice(true);
    y1.UseDevice(true);
    if (L10)
    {
      t1.SetSize(n1);
      t1.UseDevice(true);
    }
  }

  // Block 0: z0 = P0^{-1} r0
  ExtractSubVector(x, x0, 0, n0);
  pc0->Mult(x0, y0);
  InsertSubVector(y0, y, 0, n0);

  // Block 1: z1 = P1^{-1} (r1 - L10 z0)  [or z1 = P1^{-1} r1 if L10 is null]
  ExtractSubVector(x, x1, n0, n1);
  if (L10)
  {
    L10->Mult(y0, t1);
    linalg::AXPBY(1.0, x1, -1.0, t1);
    pc1->Mult(t1, y1);
  }
  else
  {
    pc1->Mult(x1, y1);
  }
  InsertSubVector(y1, y, n0, n1);
}

template class BlockDiagonalPreconditioner<Operator>;
template class BlockDiagonalPreconditioner<ComplexOperator>;

AdditivePatchSolver::AdditivePatchSolver(
    int size, std::vector<mfem::Array<int>> &&patch_dofs,
    std::vector<std::unique_ptr<mfem::Solver>> &&patch_solvers)
  : mfem::Solver(size)
{
  MFEM_VERIFY(size > 0 && !patch_dofs.empty() && patch_dofs.size() == patch_solvers.size(),
              "Additive patch correction requires compatible nonempty patches!");
  patches.reserve(patch_dofs.size());
  for (std::size_t i = 0; i < patch_dofs.size(); i++)
  {
    int previous = -1;
    bool valid = patch_solvers[i] != nullptr;
    for (const int dof : patch_dofs[i])
    {
      valid = valid && dof > previous && dof < size;
      previous = dof;
    }
    MFEM_VERIFY(valid, "Additive patch correction received invalid true DOFs!");
    patches.push_back({std::move(patch_dofs[i]), std::move(patch_solvers[i]), {}, {}});
  }
}

void AdditivePatchSolver::Restrict(const Vector &source, Vector &destination,
                                   const mfem::Array<int> &indices)
{
  destination.SetSize(indices.Size());
  const auto *source_data = source.HostRead();
  auto *destination_data = destination.HostWrite();
  for (int i = 0; i < indices.Size(); i++)
  {
    destination_data[i] = source_data[indices[i]];
  }
}

void AdditivePatchSolver::AddProlongation(const Vector &source, Vector &destination,
                                          const mfem::Array<int> &indices)
{
  MFEM_VERIFY(source.Size() == indices.Size(),
              "Additive patch correction has an inconsistent local size!");
  const auto *source_data = source.HostRead();
  auto *destination_data = destination.HostReadWrite();
  for (int i = 0; i < indices.Size(); i++)
  {
    destination_data[indices[i]] += source_data[i];
  }
}

void AdditivePatchSolver::SetPatchOperators(
    const std::vector<const Operator *> &patch_operators)
{
  MFEM_VERIFY(patch_operators.size() == patches.size(),
              "Additive patch correction received inconsistent operators!");
  for (std::size_t i = 0; i < patches.size(); i++)
  {
    MFEM_VERIFY(patch_operators[i] &&
                    patch_operators[i]->Height() == patches[i].dofs.Size() &&
                    patch_operators[i]->Width() == patches[i].dofs.Size(),
                "Additive patch correction received an invalid principal operator!");
    patches[i].solver->SetOperator(*patch_operators[i]);
  }
}

void AdditivePatchSolver::SetOperator(const Operator &full_operator)
{
  MFEM_VERIFY(full_operator.Height() == height && full_operator.Width() == width,
              "Additive patch correction requires a compatible square operator!");
}

void AdditivePatchSolver::Mult(const Vector &x, Vector &y) const
{
  MFEM_VERIFY(x.Size() == width,
              "Additive patch correction received an inconsistent vector!");
  y.SetSize(height);
  y = 0.0;
  for (auto &patch : patches)
  {
    Restrict(x, patch.rhs, patch.dofs);
    patch.correction.SetSize(patch.rhs.Size());
    patch.solver->Mult(patch.rhs, patch.correction);
    AddProlongation(patch.correction, y, patch.dofs);
  }
}

SymmetricPatchSubspacePreconditioner::SymmetricPatchSubspacePreconditioner(
    int standard_size, std::unique_ptr<mfem::Solver> &&standard_pc,
    std::unique_ptr<mfem::Solver> &&patch_pc)
  : Solver<Operator>(), standard_size(standard_size), standard_pc(std::move(standard_pc)),
    patch_pc(std::move(patch_pc))
{
  MFEM_VERIFY(this->standard_pc && this->patch_pc,
              "Symmetric patch correction requires two subspace solvers!");
}

void SymmetricPatchSubspacePreconditioner::SetSubspaceOperators(
    const Operator &full_operator, const Operator &standard_operator)
{
  bool valid = full_operator.Height() == full_operator.Width() &&
               standard_operator.Height() == standard_operator.Width() &&
               standard_operator.Height() == standard_size && standard_size >= 0 &&
               standard_size <= full_operator.Height();
  MFEM_VERIFY(valid, "Invalid operators or true DOFs for symmetric patch correction!");

  standard_pc->SetOperator(standard_operator);
  patch_pc->SetOperator(full_operator);
  SetOperator(full_operator);
}

void SymmetricPatchSubspacePreconditioner::SetOperator(const Operator &full_operator)
{
  MFEM_VERIFY(full_operator.Height() == full_operator.Width() &&
                  standard_size <= full_operator.Height(),
              "Symmetric patch correction requires a compatible square operator!");
  op = &full_operator;
  height = full_operator.Height();
  width = full_operator.Width();
}

void SymmetricPatchSubspacePreconditioner::AddStandardCorrection(const Vector &source,
                                                                 Vector &destination,
                                                                 double coefficient) const
{
  MFEM_VERIFY(source.Size() == width && destination.Size() == height,
              "Standard correction received an inconsistent vector size!");
  standard_rhs.SetSize(standard_size);
  const auto *source_data = source.HostRead();
  auto *rhs_data = standard_rhs.HostWrite();
  for (int i = 0; i < standard_size; i++)
  {
    rhs_data[i] = source_data[i];
  }
  standard_correction.SetSize(standard_size);
  standard_pc->Mult(standard_rhs, standard_correction);
  const auto *correction_data = standard_correction.HostRead();
  auto *destination_data = destination.HostReadWrite();
  for (int i = 0; i < standard_size; i++)
  {
    destination_data[i] += coefficient * correction_data[i];
  }
}

void SymmetricPatchSubspacePreconditioner::UpdateResidual(const Vector &rhs,
                                                          const Vector &solution) const
{
  action.SetSize(height);
  residual.SetSize(height);
  op->Mult(solution, action);
  linalg::AXPBY(1.0, rhs, -1.0, action);
  residual = action;
}

void SymmetricPatchSubspacePreconditioner::Mult(const Vector &x, Vector &y) const
{
  MFEM_VERIFY(op && x.Size() == width,
              "Symmetric patch correction is not configured for this vector!");
  y.SetSize(height);
  y = 0.0;

  AddStandardCorrection(x, y);
  UpdateResidual(x, y);

  patch_correction.SetSize(residual.Size());
  patch_pc->Mult(residual, patch_correction);
  linalg::AXPBY(1.0, patch_correction, 1.0, y);

  action.SetSize(height);
  op->Mult(patch_correction, action);
  AddStandardCorrection(action, y, -1.0);
}

}  // namespace palace
