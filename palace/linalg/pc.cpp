// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "pc.hpp"

#include <petsc.h>
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/mumps.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "utils/iodata.hpp"

namespace palace
{

std::unique_ptr<mfem::Solver>
ConfigurePreconditioner(const IoData &iodata, const mfem::Array<int> &dbc_marker,
                        mfem::ParFiniteElementSpaceHierarchy &fespaces,
                        mfem::ParFiniteElementSpaceHierarchy *aux_fespaces)
{
  config::LinearSolverData::Type type = iodata.solver.linear.type;
  if (type == config::LinearSolverData::Type::DEFAULT)
  {
    if (iodata.problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        (iodata.problem.type == config::ProblemData::Type::TRANSIENT &&
         iodata.solver.transient.type == config::TransientSolverData::Type::CENTRAL_DIFF))
    {
      type = config::LinearSolverData::Type::BOOMER_AMG;
    }
    else if (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
             iodata.problem.type == config::ProblemData::Type::TRANSIENT)
    {
      type = config::LinearSolverData::Type::AMS;
    }
    else
    {
      // Prefer sparse direct solver for frequency domain problems if available.
#if defined(MFEM_USE_SUPERLU)
      type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      type = config::LinearSolverData::Type::MUMPS;
#else
      type = config::LinearSolverData::Type::AMS;
#endif
    }
  }
  int print = iodata.problem.verbose - 1;
  MPI_Comm comm = fespaces.GetFESpaceAtLevel(0).GetComm();
  std::unique_ptr<mfem::Solver> pc;
  switch (type)
  {
    case config::LinearSolverData::Type::AMS:
      // Can either be the coarse solve for geometric multigrid or the solver at the finest
      // space (in which case fespaces.GetNumLevels() == 1).
      pc = std::make_unique<HypreAmsSolver>(
          iodata, fespaces.GetFESpaceAtLevel(0),
          aux_fespaces ? &aux_fespaces->GetFESpaceAtLevel(0) : nullptr, print);
      break;
    case config::LinearSolverData::Type::BOOMER_AMG:
      pc = std::make_unique<BoomerAmgSolver>(iodata, print);
      break;
    case config::LinearSolverData::Type::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
      pc = std::make_unique<SuperLUSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::STRUMPACK:
#if defined(MFEM_USE_STRUMPACK)
      pc = std::make_unique<StrumpackSolver>(comm, iodata, print);
      break;
#endif
    case config::LinearSolverData::Type::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK) &&                                   \
    (STRUMPACK_VERSION_MAJOR >= 6 && STRUMPACK_VERSION_MINOR >= 3 && \
     STRUMPACK_VERSION_PATCH > 1)
      pc = std::make_unique<StrumpackMixedPrecisionSolver>(comm, iodata, print);
#else
      MFEM_ABORT("Solver was not built with STRUMPACK support or uses STRUMPACK older than "
                 "6.3.1 which does not include mixed-precision support, please choose a "
                 "different solver!");
#endif
      break;
    case config::LinearSolverData::Type::MUMPS:
#if defined(MFEM_USE_MUMPS)
      pc = std::make_unique<MumpsSolver>(comm, iodata, print);
#else
      MFEM_ABORT(
          "Solver was not built with MUMPS support, please choose a different solver!");
#endif
      break;
    default:
      MFEM_ABORT("Unexpected type for KspPreconditioner configuration!");
      break;
  }
  if (iodata.solver.linear.mat_gmg)
  {
    // This will construct the multigrid hierarchy using pc as the coarse solver
    // (ownership of pc is transferred to the GeometricMultigridSolver). When a special
    // auxiliary space smoother for pre-/post-smoothing is not desired, h1_fespace is just
    // a nullptr.
    return std::make_unique<GeometricMultigridSolver>(iodata, std::move(pc), dbc_marker,
                                                      fespaces, aux_fespaces);
  }
  else
  {
    return pc;
  }
}

void KspPreconditioner::Init(int n)
{
  // Set up temporary vector storage.
#if defined(PETSC_USE_COMPLEX)
  if (x_.Size() == 2 * n && y_.Size() == 2 * n)
  {
    return;
  }
  x_.SetSize(2 * n);
  y_.SetSize(2 * n);
#else
  if (x_.Size() == n && y_.Size() == n)
  {
    return;
  }
  x_.SetSize(n);
  y_.SetSize(n);
#endif
}

void KspPreconditioner::SetOperator(const mfem::Operator &op)
{
  pc_->SetOperator(op);
  Init(op.Height());
}

void KspPreconditioner::SetOperator(
    const std::vector<std::unique_ptr<mfem::Operator>> &ops,
    const std::vector<std::unique_ptr<mfem::Operator>> *aux_ops)
{
  auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pc_.get());
  if (gmg)
  {
    gmg->SetOperator(ops, aux_ops);
    Init(ops.back()->Height());
  }
  else
  {
    SetOperator(*ops.back());
  }
}

void KspPreconditioner::Mult(const petsc::PetscParVector &x, petsc::PetscParVector &y) const
{
#if defined(PETSC_USE_COMPLEX)
  mfem::Vector xr_, xi_, yr_, yi_;
  mfem::Array<const mfem::Vector *> X(2);
  mfem::Array<mfem::Vector *> Y(2);
  xr_.MakeRef(x_, 0, x_.Size() / 2);
  xi_.MakeRef(x_, x_.Size() / 2, x_.Size() / 2);
  yr_.MakeRef(y_, 0, y_.Size() / 2);
  yi_.MakeRef(y_, y_.Size() / 2, y_.Size() / 2);
  X[0] = &xr_;
  X[1] = &xi_;
  Y[0] = &yr_;
  Y[1] = &yi_;
  // yr_ = 0.0;
  // yi_ = 0.0;
  x.GetToVectors(xr_, xi_);
  pc_->ArrayMult(X, Y);
  y.SetFromVectors(yr_, yi_);
#else
  // y_ = 0.0;
  x.GetToVector(x_);
  pc_->Mult(x_, y_);
  y.SetFromVector(y_);
#endif
}

PetscErrorCode KspPreconditioner::PCSetUp(PC pc)
{
  // The preconditioner operators are set up outside of the linear solve by the user, so
  // this method does nothing.
  PetscFunctionBeginUser;
  PetscFunctionReturn(0);
}

PetscErrorCode KspPreconditioner::PCApply(PC pc, Vec x, Vec y)
{
  // Apply the preconditioner. If PETSc is compiled with complex number support, the real
  // preconditioner applied in block diagonal form.
  KspPreconditioner *op;
  petsc::PetscParVector xx(x, true), yy(y, true);
  PetscFunctionBeginUser;

  PetscCall(PCShellGetContext(pc, (void **)&op));
  MFEM_VERIFY(op, "Invalid PETSc shell PC context!");
  op->Mult(xx, yy);
  PetscFunctionReturn(0);
}

PetscErrorCode KspPreconditioner::PCDestroy(PC pc)
{
  // Ownership of the preconditioner context is not inherited by the shell preconditioner,
  // so this does nothing.
  PetscFunctionBeginUser;
  PetscFunctionReturn(0);
}

}  // namespace palace
