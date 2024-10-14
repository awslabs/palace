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

class TimeDependentEBSystemOperator: public mfem::TimeDependentOperator
{
public:
  // MPI communicator.
  MPI_Comm comm;

  // System matrices and excitation RHS.
  std::unique_ptr<Operator> K, M, C;
  Vector NegJ;

  // Time dependence of current pulse for excitation: -J(t) = -g(t) J. This function
  // returns g(t).
  std::function<double(double)> &J_coef;

  // Internal objects for solution of linear systems during time stepping.
  double dt_;
  std::unique_ptr<KspSolver> kspM, kspA;
  std::unique_ptr<Operator> A, B;
  mutable Vector RHS, rhsE;
  int size_E, size_B;
  double saved_gamma; 

  // Weak curl operator.
  std::unique_ptr<Operator> weakCurl;

  // Discrete curl operator.
  const Operator *Curl;

  // Bindings to SpaceOperator functions to get the system matrix and preconditioner, and
  // construct the linear solver.
  std::function<void(double dt)> ConfigureLinearSolver;

public:
  TimeDependentEBSystemOperator(const IoData &iodata, SpaceOperator &space_op,
                                std::function<double(double)> &J_coef, double t0, 
                                mfem::TimeDependentOperator::Type type)
    : mfem::TimeDependentOperator(space_op.GetNDSpace().GetTrueVSize()+space_op.GetRTSpace().GetTrueVSize(),
                                  t0, type),
      comm(space_op.GetComm()), J_coef(J_coef)
  {
    // Get dimensions of E and B vectors
    size_E = space_op.GetNDSpace().GetTrueVSize();
    size_B = space_op.GetRTSpace().GetTrueVSize();

    // Construct the system matrices defining the linear operator. PEC boundaries are
    // handled simply by setting diagonal entries of the mass matrix for the corresponding
    // dofs. Because the Dirichlet BC is always homogeneous, no special elimination is
    // required on the RHS. Diagonal entries are set in M (so M is non-singular).
    K = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
    C = space_op.GetDampingMatrix<Operator>(Operator::DIAG_ZERO);
    M = space_op.GetMassMatrix<Operator>(Operator::DIAG_ONE);

    // Set up RHS vector for the current source term: -g(t) J, where g(t) handles the time
    // dependence.
    space_op.GetExcitationVector(NegJ);
    RHS.SetSize(size_E+size_B);
    RHS.UseDevice(true);
    rhsE.SetSize(size_E);
    rhsE.UseDevice(true);

    // Discrete curl and weak curl.
    weakCurl = space_op.GetWeakCurlMatrix<Operator>();
    Curl = &space_op.GetCurlMatrix();
    
    // Set up linear solvers.
    {
      auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
      pcg->SetInitialGuess(iodata.solver.linear.initial_guess);
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
        A = space_op.GetSystemMatrix(dt*dt, dt, 1.0, K.get(), C.get(), M.get());
        B = space_op.GetPreconditionerMatrix<Operator>(dt*dt, dt, 1.0, 0.0);

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

  // Form the RHS for the explicit formulation:
  // rhsE = -C*E + 1/mu curl B - J(t)
  // rhsB = -curl E
  void FormRHS(const Vector &u, Vector &rhs) const
  {
    Vector uE(u.GetData()     +      0, size_E);
    Vector uB(u.GetData()     + size_E, size_B);
    Vector rhsE(rhs.GetData() +      0, size_E);
    Vector rhsB(rhs.GetData() + size_E, size_B);
    weakCurl->Mult(uB, rhsE);
    if (C)
    {
      C->AddMult(uE, rhsE, -1.0);
    }
    linalg::AXPBYPCZ(1.0, rhsE, J_coef(t), NegJ, 0.0, rhsE);
    Curl->Mult(uE, rhsB);
    rhsB *= -1.0;
  }
  
  // Form the RHS of the E-field for the implicit formulation
  //(M + dt*C + dt^2*K)kE = -C*E - dt*K*E + 1/mu curl B - J(t)
  //                   kB = -curl (E+dt*kE)
  void FormRHSImplicit(const Vector &u, const double dt, Vector &rhs) const
  {
    Vector uE(u.GetData() +      0, size_E);
    Vector uB(u.GetData() + size_E, size_B);
    weakCurl->Mult(uB, rhs);
    if (C)
    {
      C->AddMult(uE, rhs, -1.0);
    }
    K->AddMult(uE, rhs, -dt);
    linalg::AXPBYPCZ(1.0, rhs, J_coef(t), NegJ, 0.0, rhs);
  }

  // Solve M du = rhs
  void Mult(const Vector &u, Vector &du) const override
  {
    FormRHS(u, RHS);
    Vector duE(du.GetData()   +      0, size_E);
    Vector duB(du.GetData()   + size_E, size_B);
    Vector rhsE(RHS.GetData() +      0, size_E);
    Vector rhsB(RHS.GetData() + size_E, size_B);
    kspM->Mult(rhsE, duE);
    duB = rhsB;
  }
         
  void ImplicitSolve(double dt, const Vector &u, Vector &k) override
  {
    // Solve for k = [kE,kB]
    // After substituting kB expression into kE equation, kE eqn is independent of kB
    // and kE can be solved using the second-order curlcurl operator 
    // (M + dt*C + dt^2*K)kE = -C*E - dt*K*E + 1/mu curl B - J(t)
    //                    kB = -curl (E+dt*kE)
    if (!kspA || dt != dt_)
    {
      // Configure the linear solver, including the system matrix and also the matrix
      // (matrices) from which the preconditioner will be constructed.
      ConfigureLinearSolver(dt);
      dt_ = dt;
      k = 0.0;
    }
    Mpi::Print("\n");
    // Solve A kE = rhsE
    FormRHSImplicit(u, dt, rhsE);
    Vector kE(k.GetData() +      0, size_E);
    Vector kB(k.GetData() + size_E, size_B);
    Vector uE(u.GetData() +      0, size_E);
    Vector uB(u.GetData() + size_E, size_B);
    kspA->Mult(rhsE, kE);

    // kB = -curl (E+dt*kE)
    linalg::AXPBYPCZ(1.0, uE, dt, kE, 0.0, rhsE);
    Curl->Mult(rhsE, kB);
    kB *= -1.0;    
  }

  void ExplicitMult(const Vector &u, Vector &v) const override
  {
    FormRHS(u, v);
  }

  // Setup A = M - gamma J = M + gamma C + gamma^2 K 
  int SUNImplicitSetup(const Vector &y, const Vector &fy,
                       int jok, int *jcur, double gamma)
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

  // Solve (Mass - dt Jacobian) x = b 
  // In the present ODE system:
  // | M 0 | - dt | -C     1/mu curl | |xE| = |bE|
  // | 0 I |      | -curl  0         | |xB|   |bB|
  // (M + dt C)xE - dt/mu curl xB = bE 
  //             -dt curl xE + xB = bB 
  // Substitute xB = bB - dt curl xE  into xE equation
  // -> (M + dt C + dt^2 K) xE = bE + dt/mu curl bB
  int SUNImplicitSolve(const Vector &b, Vector &x, double tol)
  {
    Vector bE(b.GetData() +      0, size_E);
    Vector bB(b.GetData() + size_E, size_B);
    Vector xE(x.GetData() +      0, size_E);
    Vector xB(x.GetData() + size_E, size_B);

    // Solve A xE = bE + dt/mu curl bB
    weakCurl->Mult(bB, rhsE);
    rhsE *= saved_gamma;
    rhsE += bE;
    kspA->Mult(rhsE, xE);

    // xB = bB - dt curl xE
    Curl->Mult(xE, xB);
    xB *= -saved_gamma;
    xB += bB;

    return 0;
  }

  int SUNMassSetup() 
  { 
    // Already set M in the constructor.
    return 0;
  }

  int SUNMassSolve(const Vector &b, Vector &x, double tol)
  {
    Vector bE(b.GetData() +      0, size_E);
    Vector bB(b.GetData() + size_E, size_B);
    Vector xE(x.GetData() +      0, size_E);
    Vector xB(x.GetData() + size_E, size_B);
    kspM->Mult(bE, xE);
    xB = bB;
  
    return 0;
  }

  int SUNMassMult(const Vector &x, Vector &v)
  {
    Vector vE(v.GetData() +      0, size_E);
    Vector vB(v.GetData() + size_E, size_B);
    Vector xE(x.GetData() +      0, size_E);
    Vector xB(x.GetData() + size_E, size_B);
    M->Mult(xE, vE);
    vB = xB;
   
    return 0;
  }
};

class TimeDependentCurlCurlOperator : public mfem::SecondOrderTimeDependentOperator
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
  double a0_, a1_;
  std::unique_ptr<KspSolver> kspM, kspA;
  std::unique_ptr<Operator> A, B;
  mutable Vector RHS;

  // Bindings to SpaceOperator functions to get the system matrix and preconditioner, and
  // construct the linear solver.
  std::function<void(double a0, double a1)> ConfigureLinearSolver;

public:
  TimeDependentCurlCurlOperator(const IoData &iodata, SpaceOperator &space_op,
                                std::function<double(double)> &dJ_coef, double t0,
                                mfem::TimeDependentOperator::Type type)
    : mfem::SecondOrderTimeDependentOperator(space_op.GetNDSpace().GetTrueVSize(), t0,
                                             type),
      comm(space_op.GetComm()), dJ_coef(dJ_coef)
  {
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
    RHS.SetSize(NegJ.Size());
    RHS.UseDevice(true);

    // Set up linear solvers.
    {
      auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
      pcg->SetInitialGuess(iodata.solver.linear.initial_guess);
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
      ConfigureLinearSolver = [this, &iodata, &space_op](double a0, double a1)
      {
        // Configure the system matrix and also the matrix (matrices) from which the
        // preconditioner will be constructed.
        A = space_op.GetSystemMatrix(a0, a1, 1.0, K.get(), C.get(), M.get());
        B = space_op.GetPreconditionerMatrix<Operator>(a0, a1, 1.0, 0.0);

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

  void FormRHS(const Vector &u, const Vector &du, Vector &rhs) const
  {
    // Multiply: rhs = -(K u + C du) - g'(t) J.
    K->Mult(u, rhs);
    if (C)
    {
      C->AddMult(du, rhs, 1.0);
    }
    linalg::AXPBYPCZ(-1.0, rhs, dJ_coef(t), NegJ, 0.0, rhs);
  }

  void Mult(const Vector &u, const Vector &du, Vector &ddu) const override
  {
    // Solve: M ddu = -(K u + C du) - g'(t) J.
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
    if (!kspA || a0 != a0_ || a1 != a1_)
    {
      // Configure the linear solver, including the system matrix and also the matrix
      // (matrices) from which the preconditioner will be constructed.
      ConfigureLinearSolver(a0, a1);
      a0_ = a0;
      a1_ = a1;
      k = 0.0;
    }
    Mpi::Print("\n");
    FormRHS(u, du, RHS);
    kspA->Mult(RHS, k);
  }
};

}  // namespace

SecondOrderTimeOperator::SecondOrderTimeOperator(const IoData &iodata, 
                                                 SpaceOperator &space_op,
                                                 std::function<double(double)> &dJ_coef)
{
  // Construct discrete curl matrix for B-field time integration.
  Curl = &space_op.GetCurlMatrix();

  // Allocate space for solution vectors.
  E.SetSize(Curl->Width());
  dE.SetSize(Curl->Width());
  En.SetSize(Curl->Width());
  B.SetSize(Curl->Height());
  E.UseDevice(true);
  dE.UseDevice(true);
  En.UseDevice(true);
  B.UseDevice(true);

  // Create ODE solver for 2nd-order IVP.
  mfem::TimeDependentOperator::Type type = mfem::TimeDependentOperator::EXPLICIT;
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
  }

  // Set up time-dependent operator for 2nd-order curl-curl equation for E.
  op =
      std::make_unique<TimeDependentCurlCurlOperator>(iodata, space_op, dJ_coef, 0.0, type);
}

const KspSolver &SecondOrderTimeOperator::GetLinearSolver() const
{
  const auto &curlcurl = dynamic_cast<const TimeDependentCurlCurlOperator &>(*op);
  MFEM_VERIFY(curlcurl.kspA,
              "No linear solver for time-depdendent operator has been constructed!\n");
  return *curlcurl.kspA;
}

double SecondOrderTimeOperator::GetMaxTimeStep() const
{
  const auto &curlcurl = dynamic_cast<const TimeDependentCurlCurlOperator &>(*op);
  MPI_Comm comm = curlcurl.comm;
  const Operator &M = *curlcurl.M;
  const Operator &K = *curlcurl.K;

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

void SecondOrderTimeOperator::Init(double &dt)
{
  // Always use zero initial conditions.
  E = 0.0;
  dE = 0.0;
  B = 0.0;
  ode->Init(*op);
}

void SecondOrderTimeOperator::Step(double &t, double &dt)
{
  // Single time step for E-field.
  En = E;
  ode->Step(E, dE, t, dt);

  // Trapezoidal integration for B-field: dB/dt = -∇ x E.
  En += E;
  Curl->AddMult(En, B, -0.5 * dt);
}

FirstOrderTimeOperator::FirstOrderTimeOperator(const IoData &iodata, SpaceOperator &space_op,
                             std::function<double(double)> &J_coef)
{
  // Get sizes.
  int size_E = space_op.GetNDSpace().GetTrueVSize();
  int size_B = space_op.GetRTSpace().GetTrueVSize();
  
  // Allocate space for solution vectors.
  sol.SetSize(size_E+size_B);
  E.SetSize(size_E);
  B.SetSize(size_B);
  sol.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);

  // Sol = [E, B]
  E.MakeRef(sol, 0);
  B.MakeRef(sol, size_E);

  // SUNDIALS adaptive time-stepping parameters.
  adapt_dt = iodata.solver.transient.adaptive_dt;
  rel_tol = iodata.solver.transient.rel_tol;
  abs_tol = iodata.solver.transient.abs_tol;
  rk_order = iodata.solver.transient.rk_order;
  // Minimum RK order is 2.
  rk_order = std::max(rk_order, 2);

  // Create ODE solver for 1st-order IVP.
  mfem::TimeDependentOperator::Type type = mfem::TimeDependentOperator::EXPLICIT;
  switch (iodata.solver.transient.type)
  {
    case config::TransientSolverData::Type::IMPLICIT_RK:
    case config::TransientSolverData::Type::DEFAULT:
      {
#if defined(MFEM_USE_SUNDIALS)
        // SUNDIALS ARKode solver.
        std::unique_ptr<mfem::ARKStepSolver> arkode;
        arkode = std::make_unique<mfem::ARKStepSolver>(space_op.GetComm(), 
                                                       mfem::ARKStepSolver::IMPLICIT);
        type = mfem::TimeDependentOperator::IMPLICIT;
        // Operator for first-order ODE system.
        op = std::make_unique<TimeDependentEBSystemOperator>(iodata, space_op, J_coef, 0.0, type);
        // Initialize ARKode.
        arkode->Init(*op);
        // IRK settings
        // Use implicit setup/solve defined in SUNImplicit*
        arkode->UseMFEMLinearSolver();
        // Use mass matrix operations defined in SUNMass*
        arkode->UseMFEMMassLinearSolver(0);
        // Implicit solve is linear and J is not time-dependent.
        ARKStepSetLinear(arkode->GetMem(), 0);
        // Maximum IRK order is 5.
        rk_order = std::min(rk_order, 5);
        // Relative and absolute tolerances.
        arkode->SetSStolerances(rel_tol, abs_tol);  
        // Set the order of the RK scheme.
        ARKStepSetOrder(arkode->GetMem(), rk_order); 
        // Set the ODE solver to ARKode.
        ode = std::move(arkode);
#endif
      }
      break;
    case config::TransientSolverData::Type::EXPLICIT_RK:
      {
#if defined(MFEM_USE_SUNDIALS)
        // SUNDIALS ARKode solver.
        std::unique_ptr<mfem::ARKStepSolver> arkode;
        arkode = std::make_unique<mfem::ARKStepSolver>(space_op.GetComm(), 
                                                       mfem::ARKStepSolver::EXPLICIT);
        type = mfem::TimeDependentOperator::EXPLICIT;
        // Operator for first-order ODE system.
        op = std::make_unique<TimeDependentEBSystemOperator>(iodata, space_op, J_coef, 0.0, type);
        // Initialize ARKode.
        arkode->Init(*op);
        // Maximum ERK order is 8.
        rk_order = std::min(rk_order, 8);
        // Relative and absolute tolerances.
        arkode->SetSStolerances(rel_tol, abs_tol);  
        // Set the order of the RK scheme.
        ARKStepSetOrder(arkode->GetMem(), rk_order); 
        // Set the ODE solver to ARKode.
        ode = std::move(arkode);
#endif
      }
      break;
  }
}

const KspSolver &FirstOrderTimeOperator::GetLinearSolver() const
{
  const auto &ebsystem = dynamic_cast<const TimeDependentEBSystemOperator &>(*op);
  if (isExplicit()) 
  {
    MFEM_VERIFY(ebsystem.kspM,
                "No linear solver for time-dependent operator has been constructed!\n");
    return *ebsystem.kspM;
  }
  else 
  {
    MFEM_VERIFY(ebsystem.kspA,
                "No linear solver for time-dependent operator has been constructed!\n");
    return *ebsystem.kspA;
  }
}

double FirstOrderTimeOperator::GetMaxTimeStep() const
{
  const auto &ebsystem = dynamic_cast<const TimeDependentEBSystemOperator &>(*op);
  MPI_Comm comm = ebsystem.comm;
  const Operator &M = *ebsystem.M;
  const Operator &K = *ebsystem.K;

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

void FirstOrderTimeOperator::Init(double &dt)
{
  // Always use zero initial conditions.
  sol = 0.0;
  ode->Init(*op);
#if defined(MFEM_USE_SUNDIALS)
  if (mfem::ARKStepSolver* arkode = dynamic_cast<mfem::ARKStepSolver*>(ode.get()))
  {
    // Setting a max internal time step can sometimes be beneficial.
    // We could add this as an optional user input, but what should the default be?
    //arkode->SetMaxStep(dt); 
    if(!adapt_dt) 
    {
      // Disable adaptive time stepping.
      arkode->SetFixedStep(dt);
    }
  }
#endif
}

void FirstOrderTimeOperator::Step(double &t, double &dt)
{
  double dt_input = dt;
  ode->Step(sol, t, dt);
  // Ensure user-specified dt does not change.
  dt = dt_input; 
}

void FirstOrderTimeOperator::PrintStats()
{
#if defined(MFEM_USE_SUNDIALS)
  if (mfem::ARKStepSolver* arkode = dynamic_cast<mfem::ARKStepSolver*>(ode.get()))
  {
    long int expsteps, accsteps, step_attempts, nfe_evals, nfi_evals, nlinsetups, netfails;
    ARKStepGetTimestepperStats(arkode->GetMem(), &expsteps, &accsteps,
                               &step_attempts, &nfe_evals, &nfi_evals,
                               &nlinsetups, &netfails);
    long int nniters;
    ARKStepGetNumNonlinSolvIters(arkode->GetMem(), &nniters); 

    Mpi::Print("\nARKode time-stepper statistics\n");
    Mpi::Print(" Stability-limited steps: {:d}\n", expsteps);
    Mpi::Print(" Accuracy-limited steps: {:d}\n", accsteps);
    Mpi::Print(" Calls to explicit RHS function: {:d}\n", nfe_evals);
    Mpi::Print(" Calls to implicit RHS function: {:d}\n", nfi_evals);
    Mpi::Print(" Calls to linear solver setup function: {:d}\n", nlinsetups);
    Mpi::Print(" Calls to linear solver solve function: {:d}\n", nniters);
    Mpi::Print(" Number of error test failures: {:d}\n", netfails);
  } 
#endif
}

}  // namespace palace
