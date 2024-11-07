// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "timeoperator.hpp"

#include <limits>
#include <vector>
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/solver.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

class TimeDependentFirstOrderOperator : public mfem::TimeDependentOperator
{
public:
  // MPI communicator.
  MPI_Comm comm;

  // System matrices and excitation RHS.
  std::unique_ptr<Operator> K, M, C;
  Vector NegJ;

  // Time dependence of current pulse for excitation: -J'(t) = -g'(t) J. This function
  // returns g'(t).
  std::function<double(double)> &dJ_coef;

  // Internal objects for solution of linear systems during time stepping.
  double dt_, saved_gamma;
  std::unique_ptr<KspSolver> kspM, kspA;
  std::unique_ptr<Operator> A, B;
  mutable Vector RHS;
  int size_E;

  // Bindings to SpaceOperator functions to get the system matrix and preconditioner, and
  // construct the linear solver.
  std::function<void(double dt)> ConfigureLinearSolver;

public:
  TimeDependentFirstOrderOperator(const IoData &iodata, SpaceOperator &space_op,
                                  std::function<double(double)> &dJ_coef, double t0,
                                  mfem::TimeDependentOperator::Type type)
    : mfem::TimeDependentOperator(2 * space_op.GetNDSpace().GetTrueVSize(), t0, type),
      comm(space_op.GetComm()), dJ_coef(dJ_coef)
  {
    // Get dimensions of E and Edot vectors.
    size_E = space_op.GetNDSpace().GetTrueVSize();

    // Construct the system matrices defining the linear operator. PEC boundaries are
    // handled simply by setting diagonal entries of the mass matrix for the corresponding
    // dofs. Because the Dirichlet BC is always homogeneous, no special elimination is
    // required on the RHS. Diagonal entries are set in M (so M is non-singular).
    K = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
    C = space_op.GetDampingMatrix<Operator>(Operator::DIAG_ZERO);
    M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);

    // Set up RHS vector for the current source term: -g'(t) J, where g(t) handles the time
    // dependence.
    space_op.GetExcitationVector(NegJ);
    RHS.SetSize(2 * size_E);
    RHS.UseDevice(true);

    // Set up linear solvers.
    {
      const int print = iodata.problem.verbose;
      auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
      pcg->SetInitialGuess(0);
      pcg->SetRelTol(iodata.solver.linear.tol);
      pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
      pcg->SetMaxIter(iodata.solver.linear.max_it);
      auto jac = std::make_unique<JacobiSmoother<Operator>>(comm);
      kspM = std::make_unique<KspSolver>(std::move(pcg), std::move(jac));
      kspM->SetOperators(*M, *M);
    }
    {
      // For explicit schemes, recommended to just use cheaper preconditioners. Otherwise,
      // use AMS or a direct solver. The system matrix is formed as a sequence of matrix
      // vector products, and is only assembled for preconditioning.
      ConfigureLinearSolver = [this, &iodata, &space_op](double dt)
      {
        // Configure the system matrix and also the matrix (matrices) from which the
        // preconditioner will be constructed.
        A = space_op.GetSystemMatrix(dt * dt, dt, 1.0, K.get(), C.get(), M.get());
        B = space_op.GetPreconditionerMatrix<Operator>(dt * dt, dt, 1.0, 0.0);

        // Configure the solver.
        if (!kspA)
        {
          kspA = std::make_unique<KspSolver>(iodata, space_op.GetNDSpaces(),
                                             &space_op.GetH1Spaces());
        }
        kspA->SetOperators(*A, *B);
      };
    }
  }

  // Form the RHS for the first-order ODE system
  void FormRHS(const Vector &u, Vector &rhs) const
  {
    Vector u1(u.GetData() + 0, size_E);
    Vector u2(u.GetData() + size_E, size_E);
    Vector rhs1(rhs.GetData() + 0, size_E);
    Vector rhs2(rhs.GetData() + size_E, size_E);
    // u1 = Edot, u2 = E
    // rhs1 = -(K * u2 + C * u1) - J(t)
    // rhs2 = u1
    K->Mult(u2, rhs1);
    if (C)
    {
      C->AddMult(u1, rhs1, 1.0);
    }
    linalg::AXPBYPCZ(-1.0, rhs1, dJ_coef(t), NegJ, 0.0, rhs1);

    rhs2 = u1;
  }

  // Solve M du = rhs
  // |M 0| |du1| = |-(K * u2 + C * u1) - J(t) |
  // |0 I| |du2|   | u1                       |
  void Mult(const Vector &u, Vector &du) const override
  {
    if (kspM->NumTotalMult() == 0)
    {
      // Operators have already been set in constructor.
      du = 0.0;
    }
    FormRHS(u, RHS);
    Vector du1(du.GetData() + 0, size_E);
    Vector du2(du.GetData() + size_E, size_E);
    Vector rhs1(RHS.GetData() + 0, size_E);
    Vector rhs2(RHS.GetData() + size_E, size_E);
    kspM->Mult(rhs1, du1);
    du2 = rhs2;
  }

  void ImplicitSolve(double dt, const Vector &u, Vector &k) override
  {
    // Solve: M k = f(u + dt k, t)
    // Use block elimination to avoid solving a 2n x 2n linear system
    if (!kspA || dt != dt_)
    {
      // Configure the linear solver, including the system matrix and also the matrix
      // (matrices) from which the preconditioner will be constructed.
      ConfigureLinearSolver(dt);
      dt_ = dt;
      k = 0.0;
    }
    Mpi::Print("\n");
    FormRHS(u, RHS);
    Vector k1(k.GetData() + 0, size_E);
    Vector k2(k.GetData() + size_E, size_E);
    Vector rhs1(RHS.GetData() + 0, size_E);
    Vector rhs2(RHS.GetData() + size_E, size_E);
    // A k1 = rhs1 - dt K rhs2
    K->AddMult(rhs2, rhs1, -dt);
    kspA->Mult(rhs1, k1);

    // k2 = rhs2 + dt k1
    linalg::AXPBYPCZ(1.0, rhs2, dt, k1, 0.0, k2);
  }

  void ExplicitMult(const Vector &u, Vector &v) const override { Mult(u, v); }

  // Setup A = M - gamma J = M + gamma C + gamma^2 K
  int SUNImplicitSetup(const Vector &y, const Vector &fy, int jok, int *jcur, double gamma) override
  {
    // Update Jacobian matrix
    if (!kspA || gamma != saved_gamma)
    {
      ConfigureLinearSolver(gamma);
    }

    // Indicate Jacobian was updated
    *jcur = 1;

    // Save gamma for use in solve
    saved_gamma = gamma;

    return 0;
  }

  // Solve (Mass - dt Jacobian) x = Mass b
  int SUNImplicitSolve(const Vector &b, Vector &x, double tol) override
  {
    Vector b1(b.GetData() + 0, size_E);
    Vector b2(b.GetData() + size_E, size_E);
    Vector x1(x.GetData() + 0, size_E);
    Vector x2(x.GetData() + size_E, size_E);
    Vector rhs(RHS.GetData() + 0, size_E);

    // A x1 = M b1 - dt K b2
    M->Mult(b1, rhs);
    K->AddMult(b2, rhs, -saved_gamma);
    kspA->Mult(rhs, x1);

    // x2 = b2 + dt x1
    linalg::AXPBYPCZ(1.0, b2, saved_gamma, x1, 0.0, x2);

    return 0;
  }
};

}  // namespace

TimeOperator::TimeOperator(const IoData &iodata, SpaceOperator &space_op,
                           std::function<double(double)> &dJ_coef)
  : rel_tol(iodata.solver.transient.rel_tol), abs_tol(iodata.solver.transient.abs_tol),
    order(iodata.solver.transient.order)
{
  // Construct discrete curl matrix for B-field time integration.
  Curl = &space_op.GetCurlMatrix();

  // Get sizes.
  int size_E = space_op.GetNDSpace().GetTrueVSize();
  int size_B = space_op.GetRTSpace().GetTrueVSize();

  // Allocate space for solution vectors.
  sol.SetSize(2 * size_E);
  E.SetSize(size_E);
  En.SetSize(size_E);
  B.SetSize(size_B);
  sol.UseDevice(true);
  E.UseDevice(true);
  En.UseDevice(true);
  B.UseDevice(true);

  // Sol = [Edot, E]
  E.MakeRef(sol, size_E);

  // Create ODE solver for 1st-order IVP.
  mfem::TimeDependentOperator::Type type = mfem::TimeDependentOperator::IMPLICIT;
  op = std::make_unique<TimeDependentFirstOrderOperator>(iodata, space_op, dJ_coef,
                                                         0.0, type);
  switch (iodata.solver.transient.type)
  {
    case config::TransientSolverData::Type::GEN_ALPHA:
    case config::TransientSolverData::Type::DEFAULT:
      {
        constexpr double rho_inf = 1.0;
        use_mfem_integrator = true;
        ode = std::make_unique<mfem::GeneralizedAlphaSolver>(rho_inf);
      }
      break;
    case config::TransientSolverData::Type::RUNGE_KUTTA:
      {
        constexpr int gamma_opt = 2;
        use_mfem_integrator = true;
        ode = std::make_unique<mfem::SDIRK23Solver>(gamma_opt);
      }
      break;
    case config::TransientSolverData::Type::ARKODE:
      {
#if defined(MFEM_USE_SUNDIALS)
        // SUNDIALS ARKODE solver.
        std::unique_ptr<mfem::ARKStepSolver> arkode;
        arkode = std::make_unique<mfem::ARKStepSolver>(space_op.GetComm(),
                                                       mfem::ARKStepSolver::IMPLICIT);
        // Initialize ARKODE.
        arkode->Init(*op);
        // Use implicit setup/solve defined in SUNImplicit*.
        arkode->UseMFEMLinearSolver();
        // Implicit solve is linear and J is not time-dependent.
        ARKStepSetLinear(arkode->GetMem(), 0);
        // Relative and absolute tolerances.
        arkode->SetSStolerances(rel_tol, abs_tol);
        // Set the order of the RK scheme.
        ARKStepSetOrder(arkode->GetMem(), order);
        // Set the ODE solver to ARKODE.
        ode = std::move(arkode);
#else
        MFEM_ABORT("Solver was not built with SUNDIALS support, please choose a "
                   "different transient solver type!");
#endif
      }
      break;
    case config::TransientSolverData::Type::CVODE:
      {
#if defined(MFEM_USE_SUNDIALS)
        // SUNDIALS CVODE solver.
        std::unique_ptr<mfem::CVODESolver> cvode;
        cvode = std::make_unique<mfem::CVODESolver>(space_op.GetComm(), CV_BDF);
        // Initialize CVODE.
        cvode->Init(*op);
        // Relative and absolute tolerances for time step control.
        cvode->SetSStolerances(rel_tol, abs_tol);
        // Use implicit setup/solve defined in SUNImplicit*.
        cvode->UseMFEMLinearSolver();
        // Set the max order of the multistep scheme.
        // CV_BDF can go up to 5, but >= 3 is not unconditionally stable.
        cvode->SetMaxOrder(order);
        // Set the max number of steps allowed in one CVODE step() call.
        cvode->SetMaxNSteps(10000);
        // Set the ODE solver to CVODE.
        ode = std::move(cvode);
#else
        MFEM_ABORT("Solver was not built with SUNDIALS support, please choose a "
                   "different transient solver type!");
#endif
      }
      break;
  }
}

const KspSolver &TimeOperator::GetLinearSolver() const
{
  const auto &firstOrder = dynamic_cast<const TimeDependentFirstOrderOperator &>(*op);
  MFEM_VERIFY(firstOrder.kspA,
              "No linear solver for time-dependent operator has been constructed!\n");
  return *firstOrder.kspA;
}

double TimeOperator::GetMaxTimeStep() const
{
  const auto &firstOrder = dynamic_cast<const TimeDependentFirstOrderOperator &>(*op);
  MPI_Comm comm = firstOrder.comm;
  const Operator &M = *firstOrder.M;
  const Operator &K = *firstOrder.K;

  // Solver for M⁻¹.
  constexpr double lin_tol = 1.0e-9;
  constexpr int max_lin_it = 10000;
  CgSolver<Operator> pcg(comm, 0);
  pcg.SetRelTol(lin_tol);
  pcg.SetMaxIter(max_lin_it);
  pcg.SetOperator(M);
  JacobiSmoother<Operator> jac(comm);
  jac.SetOperator(M);
  pcg.SetPreconditioner(jac);

  // Power iteration to estimate largest eigenvalue of undamped system matrix M⁻¹ K (can use
  // Hermitian eigenvalue solver as M, K are SPD).
  ProductOperator op(pcg, K);
  double lam = linalg::SpectralNorm(comm, op, true);
  MFEM_VERIFY(lam > 0.0, "Error during power iteration, λ = " << lam << "!");
  return 2.0 / std::sqrt(lam);
}

void TimeOperator::Init()
{
  // Always use zero initial conditions.
  sol = 0.0;
  B = 0.0;
  if (use_mfem_integrator)
  {
    ode->Init(*op);
  }
}

void TimeOperator::Step(double &t, double &dt)
{
  En = E;
  double dt_input = dt;
  ode->Step(sol, t, dt);
  // Ensure user-specified dt does not change.
  dt = dt_input;

  // Trapezoidal integration for B-field: dB/dt = -∇ x E.
  En += E;
  Curl->AddMult(En, B, -0.5 * dt);
}

void TimeOperator::PrintStats()
{
#if defined(MFEM_USE_SUNDIALS)
  if (mfem::ARKStepSolver *arkode = dynamic_cast<mfem::ARKStepSolver *>(ode.get()))
  {
    long int expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, netfails;
    ARKStepGetTimestepperStats(arkode->GetMem(), &expsteps, &accsteps, &step_attempts,
                               &nfe_evals, &nfi_evals, &nlinsetups, &netfails);

    long int nniters;
    ARKStepGetNumNonlinSolvIters(arkode->GetMem(), &nniters);

    Mpi::Print("\nARKODE time-stepper statistics\n");
    Mpi::Print(" Stability-limited steps: {:d}\n", expsteps);
    Mpi::Print(" Accuracy-limited steps: {:d}\n", accsteps);
    Mpi::Print(" Calls to explicit RHS function: {:d}\n", nfe_evals);
    Mpi::Print(" Calls to implicit RHS function: {:d}\n", nfi_evals);
    Mpi::Print(" Calls to linear solver setup function: {:d}\n", nlinsetups);
    Mpi::Print(" Calls to linear solver solve function: {:d}\n", nniters);
    Mpi::Print(" Number of error test failures: {:d}\n", netfails);
  }
  else if (mfem::CVODESolver *cvode = dynamic_cast<mfem::CVODESolver *>(ode.get()))
  {
    long int nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;

    // Get integrator stats.
    CVodeGetIntegratorStats(cvode->GetMem(), &nsteps, &nfevals, &nlinsetups, &netfails,
                            &qlast, &qcur, &hinused, &hlast, &hcur, &tcur);
    Mpi::Print("\n CVODE time-stepper statistics\n");
    Mpi::Print(" Number of steps: {:d}\n", nsteps);
    Mpi::Print(" Calls to RHS function: {:d}\n", nfevals);
    Mpi::Print(" Calls to linear solver setup function: {:d}\n", nlinsetups);
    Mpi::Print(" Number of error test failures: {:d}\n", netfails);
    Mpi::Print("\n");
  }
#endif
}

}  // namespace palace
