// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_KSP_SOLVER_HPP
#define PALACE_LINALG_KSP_SOLVER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class ComplexParOperator;
class ComplexVector;
class IoData;

class KspSolver : public mfem::Solver
{
protected:
  // The actual solver and preconditioner objects.
  std::unique_ptr<mfem::IterativeSolver> ksp_;
  std::unique_ptr<mfem::Solver> pc_;

private:
  // Counters for number of calls to Mult method for linear solves, and cumulative number
  // of iterations.
  mutable int ksp_mult, ksp_mult_it;

protected:
  KspSolver() : ksp_(nullptr), pc_(nullptr), ksp_mult(0), ksp_mult_it(0) {}

public:
  KspSolver(const IoData &iodata, mfem::ParFiniteElementSpaceHierarchy &fespaces,
            mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr);
  KspSolver(std::unique_ptr<mfem::IterativeSolver> &&ksp,
            std::unique_ptr<mfem::Solver> &&pc)
    : ksp_(std::move(ksp)), pc_(std::move(pc)), ksp_mult(0), ksp_mult_it(0)
  {
  }

  int NumTotalMult() const { return ksp_mult; }
  int NumTotalMultIter() const { return ksp_mult_it; }

  void SetOperator(const Operator &op) override
  {
    MFEM_ABORT("SetOperator with a single operator is not implemented for KspSolver, you "
               "must specify the preconditioner operator as well!");
  }
  virtual void
  SetOperator(const Operator &op, const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
              const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr);

  void Mult(const Vector &x, Vector &y) const override;
};

class ComplexKspSolver : public KspSolver
{
public:
  ComplexKspSolver(const IoData &iodata, mfem::ParFiniteElementSpaceHierarchy &fespaces,
                   mfem::ParFiniteElementSpaceHierarchy *aux_fespaces = nullptr);

  using KspSolver::SetOperator;
  void SetOperator(
      const Operator &op, const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
      const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr) override
  {
    MFEM_ABORT("SetOperator with a real-valued operator is not implemented for "
               "ComplexKspSolver, use the complex-valued signature instead!");
  }
  void SetOperator(const ComplexOperator &op,
                   const std::vector<std::unique_ptr<ParOperator>> &pc_ops,
                   const std::vector<std::unique_ptr<ParOperator>> *pc_aux_ops = nullptr);

  void Mult(const Vector &x, Vector &y) const override
  {
    MFEM_ABORT("Mult with a real-valued vector is not implemented for "
               "ComplexKspSolver, use the complex-valued signature instead!");
  }
  void Mult(const ComplexVector &x, ComplexVector &y) const;
};

// XX TODO REMOVE

// class IoData;
// class KspPreconditioner;

// namespace petsc
// {

// class PetscParMatrix;
// class PetscParVector;

// }  // namespace petsc

// //
// // A wrapper of PETSc's KSP class for solving linear systems.
// //
// class KspSolver
// {
// public:
//   enum class Type
//   {
//     CG,
//     CGSYM,
//     FCG,
//     MINRES,
//     GMRES,
//     FGMRES,
//     BCGS,
//     BCGSL,
//     FBCGS,
//     QMRCGS,
//     TFQMR,
//     CHOLESKY,
//     LU
//   };

// private:
//   // The actual PETSc object.
//   KSP ksp;

//   // Boolean to handle SetFromOptions calls.
//   mutable bool clcustom;

//   // Control print level for debugging.
//   int print;

//   // Print PETSc options database prior to solve.
//   bool print_opts;

//   // Check for final residual if not converged. Defaults to true.
//   bool check_final;

//   // Counter for number of calls to Mult method for a linear solve.
//   mutable PetscInt solve;

//   // Set up debugging output and configure the solver based on user specified parameters.
//   void Configure(const IoData &iodata);
//   void ConfigureVerbose(int print, const std::string &prefix);

//   // Customize object with command line options set.
//   void Customize() const;

// public:
//   // Calls PETSc's KSPCreate.
//   KspSolver(MPI_Comm comm, const IoData &iodata, const std::string &prefix =
//   std::string()); KspSolver(MPI_Comm comm, int print_lvl, const std::string &prefix =
//   std::string());

//   // Calls PETSc's KSPDestroy.
//   ~KspSolver();

//   // Sets the solver type.
//   void SetType(Type type, bool piped = false);

//   // Set solver tolerance.
//   void SetTol(PetscReal tol);

//   // Set solver tolerance.
//   void SetAbsTol(PetscReal tol);

//   // Set maximum number of iterations.
//   void SetMaxIter(PetscInt maxits);

//   // Set options specific to GMRES and FGMRES solvers.
//   void SetGMRESOptions(PetscInt maxsize, bool mgs, bool cgs2);

//   // Sets the tab level for KSP output.
//   void SetTabLevel(PetscInt l);

//   // Set flag to print PETSc options database at start of solve.
//   void SetPrintOptions(bool opts) { print_opts = opts; }

//   // Set flag to check final residual if unconverged.
//   void SetCheckFinal(bool check) { check_final = check; }

//   // Set an initial vector for the solution subspace.
//   void SetNonzeroInitialGuess(bool guess);

//   // Sets the MVP and preconditioner matrix.
//   void SetOperator(const petsc::PetscParMatrix &A, bool copy_prefix = true);

//   // Configures a shell preconditioner based on the given preconditioner object.
//   void SetPreconditioner(const KspPreconditioner &op);

//   // Application of the solver.
//   void Mult(const petsc::PetscParVector &b, petsc::PetscParVector &x) const;

//   // Call KSPReset, for example if the operator dimension has changed.
//   void Reset();

//   // Get number of solver calls.
//   PetscInt GetTotalNumMult() const;

//   // Get number of solver iterations.
//   PetscInt GetNumIter() const;
//   PetscInt GetTotalNumIter() const;

//   // Get the associated MPI communicator.
//   MPI_Comm GetComm() const;

//   // Conversion function to PETSc's KSP type.
//   operator KSP() const { return ksp; }

//   // Typecasting to PETSc object.
//   operator PetscObject() const { return reinterpret_cast<PetscObject>(ksp); }

//   // Simple static linear solve methods. The sym variable defines the matrix type: 0 for
//   // general, 1 for SPD, 2 for symmetric indefinite (definitions from MUMPS).
//   static void SolveJacobi(const petsc::PetscParMatrix &A, const petsc::PetscParVector &b,
//                           petsc::PetscParVector &x, PetscInt sym, double PetscReal
//                           = 1.0e-9, PetscInt max_it = 5000);
//   static void SolveDirect(const petsc::PetscParMatrix &A, const petsc::PetscParVector &b,
//                           petsc::PetscParVector &x, PetscInt sym);
// };

}  // namespace palace

#endif  // PALACE_LINALG_KSP_SOLVER_HPP
