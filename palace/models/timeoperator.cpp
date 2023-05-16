// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "timeoperator.hpp"

#include <vector>
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

class TimeDependentCurlCurlOperator : public mfem::SecondOrderTimeDependentOperator
{
public:
  // System matrices and excitation RHS.
  std::unique_ptr<ParOperator> K, M, C;
  Vector NegJ;

  // Time dependence of current pulse for excitation: -J'(t) = -g'(t) J. This function
  // returns g'(t).
  std::function<double(double)> &dJcoef;

  // Internal objects for solution of linear systems during time stepping.
  double a0_, a1_;
  std::unique_ptr<KspSolver> kspM, kspA;
  std::unique_ptr<ParOperator> A;
  std::vector<std::unique_ptr<ParOperator>> B, AuxB;
  mutable Vector RHS;

  // Bindings to SpaceOperator functions to get the system matrix and preconditioner, and
  // construct the linear solver.
  std::function<std::unique_ptr<KspSolver>(double a0, double a1)> ConfigureLinearSolver;

public:
  TimeDependentCurlCurlOperator(const IoData &iodata, SpaceOperator &spaceop,
                                std::function<double(double)> &djcoef, double t0,
                                mfem::TimeDependentOperator::Type type)
    : mfem::SecondOrderTimeDependentOperator(spaceop.GetNDSpace().GetTrueVSize(), t0, type),
      dJcoef(djcoef)
  {
    // Construct the system matrices defining the linear operator. PEC boundaries are
    // handled simply by setting diagonal entries of the mass matrix for the corresponding
    // dofs. Because the Dirichlet BC is always homogenous, no special elimination is
    // required on the RHS. Diagonal entries are set in M (so M is non-singular).
    K = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::STIFFNESS,
                                Operator::DIAG_ZERO);
    M = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::MASS, Operator::DIAG_ONE);
    C = spaceop.GetSystemMatrix(SpaceOperator::OperatorType::DAMPING, Operator::DIAG_ZERO);

    // Set up RHS vector for the current source term: -g'(t) J, where g(t) handles the time
    // dependence.
    spaceop.GetExcitationVector(NegJ);
    RHS.SetSize(NegJ.Size());

    // Set up linear solvers.
    {
      // PCG with a simple Jacobi preconditioner for mass matrix systems.
      auto pcg = std::make_unique<mfem::CGSolver>(M->GetComm());
      pcg->iterative_mode = iodata.solver.linear.ksp_initial_guess;
      pcg->SetRelTol(iodata.solver.linear.tol);
      pcg->SetMaxIter(iodata.solver.linear.max_it);
      pcg->SetPrintLevel(0);
      kspM =
          std::make_unique<KspSolver>(std::move(pcg), std::make_unique<JacobiSmoother>());
      kspM->SetOperator(*M, *M);
    }
    {
      // For explicit schemes, recommended to just use cheaper preconditioners. Otherwise,
      // use AMS or a direct solver. The system matrix is formed as a sequence of matrix
      // vector products, and is only assembled for preconditioning.
      ConfigureLinearSolver = [this, &iodata,
                               &spaceop](double a0, double a1) -> std::unique_ptr<KspSolver>
      {
        // Configure the system matrix and also the matrix (matrices) from which the
        // preconditioner will be constructed.
        A = spaceop.GetSystemMatrix(a0, a1, 1.0, K.get(), C.get(), M.get());
        spaceop.GetPreconditionerMatrix(a0, a1, 1.0, 0.0, B, AuxB);

        // Configure the solver.
        auto ksp = std::make_unique<KspSolver>(iodata, spaceop.GetNDSpaces(),
                                               &spaceop.GetH1Spaces());
        ksp->SetOperator(*A, B, &AuxB);
        return ksp;
      };
    }
  }

  void FormRHS(const Vector &u, const Vector &du, Vector &rhs) const
  {
    // Multiply: rhs = -(K u + C du) - g'(t) J.
    K->Mult(u, rhs);
    if (C)
    {
      C->AddMult(du, rhs, 1.0);
    }
    add(-1.0, rhs, dJcoef(t), NegJ, rhs);
  }

  void Mult(const Vector &u, const Vector &du, Vector &ddu) const override
  {
    // Solve: M ddu = -(K u + C du) - g'(t) J.
    Mpi::Print("\n");
    if (kspM->NumTotalMult() == 0)
    {
      // Operators have already been set in constructor.
      ddu = 0.0;
    }
    FormRHS(u, du, RHS);
    kspM->Mult(RHS, ddu);
  }

  void ImplicitSolve(const double a0, const double a1, const Vector &u, const Vector &du,
                     Vector &k) override
  {
    // Solve: (a0 K + a1 C + M) k = -(K u + C du) - g'(t) J, where a0 may be 0 in the
    // explicit case. At first iteration, construct the solver. Also don't print a newline
    // if already done by the mass matrix solve at the first iteration.
    if (kspA && kspA->NumTotalMult() > 0)
    {
      Mpi::Print("\n");
    }
    if (!kspA || a0 != a0_ || a1 != a1_)
    {
      // Configure the linear solver, including the system matrix and also the matrix
      // (matrices) from which the preconditioner will be constructed.
      kspA = ConfigureLinearSolver(a0, a1);
      a0_ = a0;
      a1_ = a1;
      k = 0.0;
    }
    FormRHS(u, du, RHS);
    kspA->Mult(RHS, k);
  }
};

}  // namespace

TimeOperator::TimeOperator(const IoData &iodata, SpaceOperator &spaceop,
                           std::function<double(double)> &djcoef)
{
  // Construct discrete curl matrix for B-field time integration.
  Curl = spaceop.GetCurlMatrix();

  // Allocate space for solution vectors.
  E.SetSize(Curl->Width());
  dE.SetSize(Curl->Width());
  En.SetSize(Curl->Width());
  B.SetSize(Curl->Height());

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
  op = std::make_unique<TimeDependentCurlCurlOperator>(iodata, spaceop, djcoef, 0.0, type);
}

const KspSolver &TimeOperator::GetLinearSolver() const
{
  const auto &curlcurl = dynamic_cast<const TimeDependentCurlCurlOperator &>(*op);
  MFEM_VERIFY(curlcurl.kspA,
              "No linear solver for time-depdendent operator has been constructed!\n");
  return *curlcurl.kspA;
}

double TimeOperator::GetMaxTimeStep() const
{
  const auto &curlcurl = dynamic_cast<const TimeDependentCurlCurlOperator &>(*op);
  const ParOperator &M = *curlcurl.M;
  const ParOperator &K = *curlcurl.K;

  // Solver for M⁻¹.
  constexpr double lin_tol = 1.0e-9;
  constexpr int max_lin_it = 500;
  mfem::CGSolver pcg(M.GetComm());
  pcg.SetRelTol(lin_tol);
  pcg.SetMaxIter(max_lin_it);
  pcg.SetPrintLevel(0);
  pcg.SetOperator(M);

  JacobiSmoother jac;
  jac.SetOperator(M);
  pcg.SetPreconditioner(jac);

  // Power iteration to estimate largest eigenvalue of undamped system matrix M⁻¹ K.
  SymmetricProductOperator op(pcg, K);
  double lam = linalg::SpectralNorm(M.GetComm(), op, false);
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
  En += E;
  Curl->AddMult(En, B, -0.5 * dt);
}

}  // namespace palace
