// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "timeoperator.hpp"

#include <vector>
#include "linalg/gmg.hpp"
#include "linalg/pc.hpp"
#include "linalg/petsc.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

class CurlCurlSystemOperator : public mfem::SecondOrderTimeDependentOperator
{
private:
  // MPI communicator for the parallel operators.
  MPI_Comm comm;

  // System matrices and excitation RHS.
  std::unique_ptr<mfem::Operator> K, M, C;
  mfem::Vector NegJ;

  // Reference to essential boundary true degrees of freedom from SpaceOperator (not owned).
  const mfem::Array<int> dbc_tdof_list;

  // Time dependence of current pulse for excitation: -J'(t) = -g'(t) J. This function
  // returns g'(t).
  std::function<double(double)> &dJcoef;

  // Internal objects for solution of linear systems during time stepping.
  mutable double a0_, a1_;
  mutable mfem::Vector RHS;
  mutable std::vector<std::unique_ptr<mfem::Operator>> P, AuxP;
  std::function<void(double, double, std::vector<std::unique_ptr<mfem::Operator>> &,
                     std::vector<std::unique_ptr<mfem::Operator>> &)>
      GetPreconditionerMatrix;

  // Linear system solvers and settings for implicit time integration.
  std::unique_ptr<mfem::IterativeSolver> kspM, kspA;
  std::unique_ptr<mfem::Solver> pcM, pcA;
  mutable int kspM_mult, kspA_mult, kspM_it, kspA_it;

  void FormRHS(const mfem::Vector &u, const mfem::Vector &du, mfem::Vector &rhs) const
  {
    // Multiply: rhs = -(K u + C du) - g'(t) J.
    rhs = 0.0;
    K->AddMult(u, rhs, -1.0);
    if (C)
    {
      C->AddMult(du, rhs, -1.0);
    }
    rhs.Add(dJcoef(t), NegJ);
  }

public:
  CurlCurlSystemOperator(const IoData &iodata, SpaceOperator &spaceop,
                         std::function<double(double)> &djcoef, double t0,
                         mfem::TimeDependentOperator::Type type)
    : mfem::SecondOrderTimeDependentOperator(spaceop.GetNDSpace().GetTrueVSize(), t0, type),
      comm(spaceop.GetNDSpace().GetComm()), dbc_tdof_list(spaceop.GetDbcTDofList()),
      dJcoef(djcoef)
  {
    // Construct the system matrices defining the linear operator. PEC boundaries are
    // handled simply by setting diagonal entries of the mass matrix for the corresponding
    // dofs. Because the Dirichlet BC is always homogenous, no special elimination is
    // required on the RHS. Diagonal entries are set in M (so M is non-singular).
    K = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::STIFFNESS,
                                mfem::Operator::DIAG_ZERO);
    M = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::MASS,
                                mfem::Operator::DIAG_ONE);
    C = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::DAMPING,
                                mfem::Operator::DIAG_ZERO);

    // Set up RHS vector for the current source term: -g'(t) J, where g(t) handles the time
    // dependence.
    spaceop.GetTimeDomainExcitationVector(NegJ);
    RHS.SetSize(NegJ.Size());

    // Set up linear solvers (SetOperator will be called later on at first time step).
    mfem::IterativeSolver::PrintLevel print =
        mfem::IterativeSolver::PrintLevel().Warnings().Errors();
    if (iodata.problem.verbose > 0)
    {
      print.Summary();
      if (iodata.problem.verbose > 1)
      {
        print.Iterations();
        if (iodata.problem.verbose > 2)
        {
          print.All();
        }
      }
    }
    {
      // PCG with a simple smoother preconditioner for mass matrix systems.
      mfem::Vector diag(M->Height());
      M->AssembleDiagonal(diag);
      pcM = std::make_unique<mfem::OperatorJacobiSmoother>(diag, spaceop.GetDbcTDofList());

      auto pcg = std::make_unique<mfem::CGSolver>(comm);
      pcg->iterative_mode = iodata.solver.linear.ksp_initial_guess;
      pcg->SetRelTol(iodata.solver.linear.tol);
      pcg->SetMaxIter(iodata.solver.linear.max_it);
      pcg->SetPrintLevel(0);
      pcg->SetOperator(*M);
      pcg->SetPreconditioner(*pcM);
      kspM = std::move(pcg);
    }
    {
      // For explicit schemes, recommended to just use cheaper preconditioners. Otherwise,
      // use AMS or a direct solver. The system matrix is formed as a sequence of matrix
      // vector products, and is only assembled for preconditioning.
      pcA = ConfigurePreconditioner(iodata, spaceop.GetDbcMarker(), spaceop.GetNDSpaces(),
                                    &spaceop.GetH1Spaces());

      auto pcg = std::make_unique<mfem::CGSolver>(comm);
      pcg->iterative_mode = iodata.solver.linear.ksp_initial_guess;
      pcg->SetRelTol(iodata.solver.linear.tol);
      pcg->SetMaxIter(iodata.solver.linear.max_it);
      pcg->SetPrintLevel(print);
      pcg->SetOperator(*this);
      pcg->SetPreconditioner(*pcA);
      kspA = std::move(pcg);
      if (iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::DEFAULT &&
          iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::CG)
      {
        Mpi::Warning("Transient problem type always uses CG as the Krylov solver!\n");
      }

      // The assembled matrix for preconditioning is constructed as a function of the
      // coefficients defined by the time integrator.
      GetPreconditionerMatrix = [&](double a0, double a1,
                                    std::vector<std::unique_ptr<mfem::Operator>> &B,
                                    std::vector<std::unique_ptr<mfem::Operator>> &AuxB)
      { spaceop.GetPreconditionerMatrix(a0, a1, B, AuxB, true); };
    }
    kspM_mult = kspA_mult = kspM_it = kspA_it = 0;
  }

  MPI_Comm GetComm() const { return comm; }
  const mfem::Operator &GetK() const { return *K; }
  const mfem::Operator &GetM() const { return *M; }
  const mfem::Operator &GetC() const { return *C; }
  const mfem::Array<int> &GetDbcTDofList() const { return dbc_tdof_list; }

  int GetNumMult() const { return kspM_mult; }
  int GetNumMultIter() const { return kspM_it; }
  int GetNumImplicitSolve() const { return kspA_mult; }
  int GetNumImplicitSolveIter() const { return kspA_it; }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    // Multiply: y = (a0 K + a1 C + M) x.
    M->Mult(x, y);
    K->AddMult(x, y, a0_);
    if (C)
    {
      C->AddMult(x, y, a1_);
    }
  }

  void Mult(const mfem::Vector &u, const mfem::Vector &du, mfem::Vector &ddu) const override
  {
    // Solve: M ddu = -(K u + C du) - g'(t) J.
    Mpi::Print("\n");
    if (kspM_mult == 0)
    {
      // Operators have already been set in constructor.
      ddu = 0.0;
    }
    FormRHS(u, du, RHS);
    kspM->Mult(RHS, ddu);
    if (!kspM->GetConverged())
    {
      Mpi::Warning("Linear solver did not converge in {:d} iterations!\n",
                   kspM->GetNumIterations());
    }
    kspM_mult++;
    kspM_it += kspM->GetNumIterations();
  }

  void ImplicitSolve(const double a0, const double a1, const mfem::Vector &u,
                     const mfem::Vector &du, mfem::Vector &k) override
  {
    // Solve: (a0 K + a1 C + M) k = -(K u + C du) - g'(t) J, where a0 may be 0 in the
    // explicit case. At first iteration, construct the solver. Also don't print a newline
    // if already done by the mass matrix solve at the first iteration.
    if (kspA_mult > 0)
    {
      Mpi::Print("\n");
    }
    if (kspA_mult == 0 || a0 != a0_ || a1 != a1_)
    {
      // Configure the matrix (matrices) from which the preconditioner will be constructed.
      GetPreconditionerMatrix(a0, a1, P, AuxP);
      auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pcA.get());
      if (gmg)
      {
        gmg->SetOperator(P, &AuxP);
      }
      else
      {
        pcA->SetOperator(*P.back());
      }
      a0_ = a0;
      a1_ = a1;
      k = 0.0;
    }
    FormRHS(u, du, RHS);
    kspA->Mult(RHS, k);
    if (!kspA->GetConverged())
    {
      Mpi::Warning("Linear solver did not converge in {:d} iterations!\n",
                   kspA->GetNumIterations());
    }
    kspA_mult++;
    kspA_it += kspA->GetNumIterations();
  }
};

class SymmetricProductOperator : public mfem::Operator
{
private:
  const mfem::Operator &A, &B;
  mutable mfem::Vector z;

public:
  SymmetricProductOperator(const mfem::Operator &opA, const mfem::Operator &opB)
    : mfem::Operator(opA.Height(), opB.Width()), A(opA), B(opB), z(opB.Height())
  {
  }

  void Mult(const mfem::Vector &x, mfem::Vector &y) const override
  {
    B.Mult(x, z);
    A.Mult(z, y);
  }

  void MultTranspose(const mfem::Vector &x, mfem::Vector &y) const override
  {
    A.Mult(x, z);
    B.Mult(z, y);
  }
};

}  // namespace

TimeOperator::TimeOperator(const IoData &iodata, SpaceOperator &spaceop,
                           std::function<double(double)> &djcoef)
{
  // Construct discrete curl matrix for B-field time integration.
  NegCurl = spaceop.GetNegCurlMatrix();

  // Allocate space for solution vectors.
  E.SetSize(NegCurl->Width());
  dE.SetSize(NegCurl->Width());
  En.SetSize(NegCurl->Width());
  B.SetSize(NegCurl->Height());

  // Create ODE solver for 2nd-order IVP.
  mfem::TimeDependentOperator::Type type;
  switch (iodata.solver.transient.type)
  {
    case config::TransientSolverData::Type::GEN_ALPHA:
    case config::TransientSolverData::Type::DEFAULT:
      {
        constexpr double rho_inf = 1.0;
        ode = std::make_unique<mfem::GeneralizedAlpha2Solver>(rho_inf);
        type = mfem::TimeDependentOperator::IMPLICIT;
      }
      break;
    case config::TransientSolverData::Type::NEWMARK:
      {
        constexpr double beta = 0.25, gamma = 0.5;
        ode = std::make_unique<mfem::NewmarkSolver>(beta, gamma);
        type = mfem::TimeDependentOperator::IMPLICIT;
      }
      break;
    case config::TransientSolverData::Type::CENTRAL_DIFF:
      {
        ode = std::make_unique<mfem::CentralDifferenceSolver>();
        type = mfem::TimeDependentOperator::EXPLICIT;
      }
      break;
    default:
      MFEM_ABORT("Invalid transient solver type!");
      type = mfem::TimeDependentOperator::EXPLICIT;  // For compiler warning
      break;
  }

  // Set up time-dependent operator for 2nd-order curl-curl equation for E.
  op = std::make_unique<CurlCurlSystemOperator>(iodata, spaceop, djcoef, 0.0, type);
}

int TimeOperator::GetTotalKspMult() const
{
  const auto &curlcurl = dynamic_cast<CurlCurlSystemOperator &>(*op);
  return curlcurl.GetNumMult() + curlcurl.GetNumImplicitSolve();
}

int TimeOperator::GetTotalKspIter() const
{
  const auto &curlcurl = dynamic_cast<CurlCurlSystemOperator &>(*op);
  return curlcurl.GetNumMultIter() + curlcurl.GetNumImplicitSolveIter();
}

double TimeOperator::GetMaxTimeStep() const
{
  const auto &curlcurl = dynamic_cast<CurlCurlSystemOperator &>(*op);
  const mfem::Operator &M = curlcurl.GetM();
  const mfem::Operator &K = curlcurl.GetK();

  // Solver for M⁻¹.
  constexpr double lin_tol = 1.0e-9;
  constexpr int max_lin_it = 500;
  mfem::CGSolver pcg(curlcurl.GetComm());
  pcg.SetRelTol(lin_tol);
  pcg.SetMaxIter(max_lin_it);
  pcg.SetPrintLevel(0);
  pcg.SetOperator(M);

  mfem::Vector diag(M.Height());
  M.AssembleDiagonal(diag);
  mfem::OperatorJacobiSmoother prec(diag, curlcurl.GetDbcTDofList());
  pcg.SetPreconditioner(prec);

  // Power iteration to estimate largest eigenvalue of undamped system matrix M⁻¹ K.
  petsc::PetscShellMatrix MinvK(curlcurl.GetComm(),
                                std::make_unique<SymmetricProductOperator>(pcg, K));
  double lam = MinvK.Norm2();
  MFEM_VERIFY(lam > 0.0, "Error during power iteration, λ = " << lam << "!");
  return 2.0 / std::sqrt(lam);
}

void TimeOperator::Init()
{
  // Always use zero initial conditions.
  E = 0.0;
  dE = 0.0;
  B = 0.0;
  ode->Init(*op);
}

void TimeOperator::Step(double &t, double &dt)
{
  // Single time step for E-field.
  En = E;
  ode->Step(E, dE, t, dt);

  // Trapezoidal integration for B-field: dB/dt = -∇ x E.
  En.Add(1.0, E);
  NegCurl->AddMult(En, B, 0.5 * dt);
}

}  // namespace palace
