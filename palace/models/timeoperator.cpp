// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "timeoperator.hpp"

#include <vector>
#include "linalg/gmg.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/pc.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

class TimeDependentCurlCurlOperator : public mfem::SecondOrderTimeDependentOperator
{
private:
  // MPI communicator for the parallel operators.
  MPI_Comm comm;

  // System matrices and excitation RHS.
  std::unique_ptr<ParOperator> K, M, C;
  Vector NegJ;

  // Time dependence of current pulse for excitation: -J'(t) = -g'(t) J. This function
  // returns g'(t).
  std::function<double(double)> &dJcoef;

  // Internal objects for solution of linear systems during time stepping.
  double a0_, a1_;
  std::unique_ptr<ParOperator> A;
  std::vector<std::unique_ptr<ParOperator>> B, AuxB;
  mutable Vector RHS;

  // XX TODO REMOVE
  //  std::function<std::unique_ptr<ParOperator>(double, double)> GetSystemMatrix;
  //  std::function<void(double, double, std::vector<std::unique_ptr<ParOperator>> &,
  //                     std::vector<std::unique_ptr<ParOperator>> &)>
  //    GetPreconditionerMatrix;
  //  std::function<std::unique_ptr<mfem::IterativeSolver>(double, double)> GetSystemMatrix;

  // Linear system solvers and settings for implicit time integration.
  std::unique_ptr<mfem::IterativeSolver> kspM, kspA;
  std::unique_ptr<mfem::Solver> pcM, pcA;
  mutable int kspM_mult, kspA_mult, kspM_it, kspA_it;

  // Bindings to SpaceOperator functions to get the system matrix and preconditioner, and
  // construct the linear solver.
  std::function<std::unique_ptr<mfem::IterativeSolver>(double a0, double a1)>
      ConfigureLinearSolver;

  void FormRHS(const Vector &u, const Vector &du, Vector &rhs) const
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
  TimeDependentCurlCurlOperator(const IoData &iodata, SpaceOperator &spaceop,
                                std::function<double(double)> &djcoef, double t0,
                                mfem::TimeDependentOperator::Type type)
    : mfem::SecondOrderTimeDependentOperator(spaceop.GetNDSpace().GetTrueVSize(), t0, type),
      comm(spaceop.GetNDSpace().GetComm()), dJcoef(djcoef)
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
      Vector diag(M->Height());
      M->AssembleDiagonal(diag);

      // XX TODO: Should not need DBC TDOF LIST as the diagonal is already 1 upon
      // assembly... (see ParOperator)
      //          Maybe avoid MFEM's JAcobi smoother and write our own like in Chebyshev??
      // pcM = std::make_unique<mfem::OperatorJacobiSmoother>(diag,
      // spaceop.GetDbcTDofList());
      pcM = std::make_unique<JacobiSmoother>(diag);

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

      // // XX TODO ADDRESS, WITH BCS, ETC.....
      // pcA = ConfigurePreconditioner(iodata, spaceop.GetDbcMarker(),
      // spaceop.GetNDSpaces(),
      //                               &spaceop.GetH1Spaces());

      // XX TODO TEST IF THE BELOW WORKS?

      // auto pcg = std::make_unique<mfem::CGSolver>(comm);
      // pcg->iterative_mode = iodata.solver.linear.ksp_initial_guess;
      // pcg->SetRelTol(iodata.solver.linear.tol);
      // pcg->SetMaxIter(iodata.solver.linear.max_it);
      // pcg->SetPrintLevel(print);
      // pcg->SetOperator(*this);
      // pcg->SetPreconditioner(*pcA);
      // kspA = std::move(pcg);

      // XX TODO REMOVE
      //  GetSystemMatrix = [this, &spaceop](double a0, double a1) ->
      //  std::unique_ptr<ParOperator>
      //  {
      //    return spaceop.GetSystemMatrix(a0, a1, 1.0, this->K.get(), this->C.get(),
      //    this->M.get());
      //  }
      //  GetPreconditionerMatrix = [&spaceop](double a0, double a1,
      //                                std::vector<std::unique_ptr<ParOperator>> &B,
      //                                std::vector<std::unique_ptr<ParOperator>> &AuxB)
      //  { spaceop.GetPreconditionerMatrix(a0, a1, 1.0, 0.0, B, AuxB); };
      //  ConfigureLinearSolver = [=](std::unique_ptr<ParOperator> &A,
      //  std::unique_ptr<mfem::Solver> &pc) -> std::unique_ptr<mfem::IterativeSolver>
      //  {
      //    auto pcg = std::make_unique<mfem::CGSolver>(comm);
      //    pcg->iterative_mode = iodata.solver.linear.ksp_initial_guess;
      //    pcg->SetRelTol(iodata.solver.linear.tol);
      //    pcg->SetMaxIter(iodata.solver.linear.max_it);
      //    pcg->SetPrintLevel(print);
      //    pcg->SetOperator(*A);
      //    pcg->SetPreconditioner(*pc);
      //  }

      // The time domain system matrix is A = a0 K + a1 C + M, which constructed using the
      // assembled K, C, and M matrices and the coefficients a0 and a1 defined by the time
      // integrator.
      if (iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::DEFAULT &&
          iodata.solver.linear.ksp_type != config::LinearSolverData::KspType::CG)
      {
        Mpi::Warning("Transient problem type always uses CG as the Krylov solver!\n");
      }
      bool iterative_mode = iodata.solver.linear.ksp_initial_guess;
      double tol = iodata.solver.linear.tol;
      int max_it = iodata.solver.linear.max_it;
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
      ConfigureLinearSolver = [this, &spaceop, iterative_mode, tol, max_it,
                               print](double a0,
                                      double a1) -> std::unique_ptr<mfem::IterativeSolver>
      {
        // Configure the system matrix and also the matrix (matrices) from which the
        // preconditioner will be constructed.
        this->A = spaceop.GetSystemMatrix(a0, a1, 1.0, this->K.get(), this->C.get(),
                                          this->M.get());
        spaceop.GetPreconditionerMatrix(a0, a1, 1.0, 0.0, this->B, this->AuxB);

        // Configure the preconditioner.
        auto *gmg = dynamic_cast<GeometricMultigridSolver *>(this->pcA.get());

        // XX TODO WIP
        //  if (gmg)
        //  {
        //    gmg->SetOperator(this->B, &this->AuxB);
        //  }
        //  else
        //  {
        //    this->pcA->SetOperator(*this->B.back());
        //  }

        // Construct and return the linear solver.
        auto pcg = std::make_unique<mfem::CGSolver>(this->comm);
        pcg->iterative_mode = iterative_mode;
        pcg->SetRelTol(tol);
        pcg->SetMaxIter(max_it);
        pcg->SetPrintLevel(print);
        pcg->SetOperator(*this->A);
        pcg->SetPreconditioner(*this->pcA);
        return pcg;
      };
    }
    kspM_mult = kspA_mult = kspM_it = kspA_it = 0;
  }

  MPI_Comm GetComm() const { return comm; }
  const ParOperator &GetK() const { return *K; }
  const ParOperator &GetM() const { return *M; }
  const ParOperator &GetC() const { return *C; }

  int GetNumMult() const { return kspM_mult; }
  int GetNumMultIter() const { return kspM_it; }
  int GetNumImplicitSolve() const { return kspA_mult; }
  int GetNumImplicitSolveIter() const { return kspA_it; }

  void Mult(const Vector &u, const Vector &du, Vector &ddu) const override
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

  void ImplicitSolve(const double a0, const double a1, const Vector &u, const Vector &du,
                     Vector &k) override
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
      // Configure the linear solver, including the system matrix and also the matrix
      // (matrices) from which the preconditioner will be constructed.
      kspA = ConfigureLinearSolver(a0, a1);

      // XX TODO WORKING: REMOVE THE BELOW IF THIS WORKS...

      // A = GetSystemMatrix(a0, a1);
      // GetPreconditionerMatrix(a0, a1, P, AuxP);
      // auto *gmg = dynamic_cast<GeometricMultigridSolver *>(pcA.get());
      // if (gmg)
      // {
      //   gmg->SetOperator(P, &AuxP);
      // }
      // else
      // {
      //   pcA->SetOperator(*P.back());
      // }

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

int TimeOperator::GetTotalKspMult() const
{
  const auto &curlcurl = dynamic_cast<TimeDependentCurlCurlOperator &>(*op);
  return curlcurl.GetNumMult() + curlcurl.GetNumImplicitSolve();
}

int TimeOperator::GetTotalKspIter() const
{
  const auto &curlcurl = dynamic_cast<TimeDependentCurlCurlOperator &>(*op);
  return curlcurl.GetNumMultIter() + curlcurl.GetNumImplicitSolveIter();
}

double TimeOperator::GetMaxTimeStep() const
{
  const auto &curlcurl = dynamic_cast<TimeDependentCurlCurlOperator &>(*op);
  const ParOperator &M = curlcurl.GetM();
  const ParOperator &K = curlcurl.GetK();

  // Solver for M⁻¹.
  constexpr double lin_tol = 1.0e-9;
  constexpr int max_lin_it = 500;
  mfem::CGSolver pcg(curlcurl.GetComm());
  pcg.SetRelTol(lin_tol);
  pcg.SetMaxIter(max_lin_it);
  pcg.SetPrintLevel(0);
  pcg.SetOperator(M);

  Vector diag(M.Height());
  M.AssembleDiagonal(diag);
  JacobiSmoother prec(diag);
  pcg.SetPreconditioner(prec);

  // Power iteration to estimate largest eigenvalue of undamped system matrix M⁻¹ K.
  SymmetricProductOperator op(pcg, K);
  double lam = linalg::SpectralNorm(curlcurl.GetComm(), op, false);
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
  Curl->AddMult(En, B, -0.5 * dt);
}

}  // namespace palace
