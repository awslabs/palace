// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeanalysissolver.hpp"

#include <complex>

#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/mumps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

std::pair<ErrorIndicator, long long int>
ModeAnalysisSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Mode analysis requires a 2D mesh (waveguide cross-section).
  MFEM_VERIFY(mesh.back()->Dimension() == 2,
              "ModeAnalysis solver requires a 2D mesh (waveguide cross-section)!");

  const auto &ma_data = iodata.solver.mode_analysis;
  const double freq_GHz = ma_data.freq;
  const int num_modes = ma_data.n;
  const double tol = ma_data.tol;

  // Nondimensionalize the operating frequency.
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz (omega = "
             "{:.6e})\n",
             freq_GHz, omega);

  // Construct FE spaces: ND for tangential E, H1 for normal (out-of-plane) E component.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                         mesh.back()->Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                         mesh.back()->Dimension());
  FiniteElementSpace nd_fespace(*mesh.back(), nd_fec.get());
  FiniteElementSpace h1_fespace(*mesh.back(), h1_fec.get());

  // Apply PEC boundary conditions (essential BCs).
  const auto &pmesh = mesh.back()->Get();
  int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> dbc_marker =
      palace::mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.pec.attributes);

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  const int nd_size = nd_fespace.GetTrueVSize();
  const int h1_size = h1_fespace.GetTrueVSize();
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             nd_fespace.GlobalTrueVSize(), h1_fespace.GlobalTrueVSize(),
             nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize());

  // Construct material operator.
  MaterialOperator mat_op(iodata, *mesh.back());

  // Compute shift parameter: sigma = -omega^2 * max(eps_r * mu_r) with safety factor.
  double c_min = mat_op.GetLightSpeedMax().Min();
  Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light in ModeAnalysisSolver!");
  const double mu_eps_max = 1.0 / (c_min * c_min) * 1.1;
  const double sigma = -omega * omega * mu_eps_max;
  Mpi::Print(" Shift: sigma = {:.6e}\n", sigma);

  // Assemble block system matrices.
  // Att = (mu^-1 curl et, curl ft) - omega^2 (eps et, ft) - sigma (mu^-1 et, ft)
  MaterialPropertyCoefficient muinv_cc(mat_op.GetAttributeToMaterial(),
                                       mat_op.GetCurlCurlInvPermeability());
  MaterialPropertyCoefficient eps_shifted(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetPermittivityReal(), -omega * omega);
  eps_shifted.AddCoefficient(mat_op.GetAttributeToMaterial(),
                             mat_op.GetInvPermeability(), -sigma);
  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc, eps_shifted);
  auto Att = ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  // Atn = -(mu^-1 grad en, ft)
  MaterialPropertyCoefficient muinv_neg(mat_op.GetAttributeToMaterial(),
                                        mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_neg);
  auto Atn =
      ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
          .StealParallelAssemble();

  // Ant = -(eps et, grad fn)
  MaterialPropertyCoefficient eps_pos(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetPermittivityReal(), 1.0);
  BilinearForm ant(nd_fespace, h1_fespace);
  ant.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(eps_pos);
  auto Ant =
      ParOperator(ant.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
          .StealParallelAssemble();

  // Ann = -(eps_zz en, fn) — uses the z-z (out-of-plane) permittivity component (1×1).
  // For a 2D mesh, the normal (propagation) direction is z-hat. The scalar permittivity
  // is stored as mat_epsilon_scalar from the (2,2) entry of the full 3×3 tensor.
  MaterialPropertyCoefficient eps_neg(mat_op.GetAttributeToMaterial(),
                                      mat_op.GetPermittivityScalar(), -1.0);
  BilinearForm ann(h1_fespace);
  ann.AddDomainIntegrator<MassIntegrator>(eps_neg);
  auto Ann = ParOperator(ann.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  // Btt = (mu^-1 et, ft)
  MaterialPropertyCoefficient muinv_pos(mat_op.GetAttributeToMaterial(),
                                        mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_pos);
  auto Btt = ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  // Build 2x2 block A matrix and eliminate essential DOFs.
  std::unique_ptr<mfem::HypreParMatrix> opAr;
  {
    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    blocks(0, 0) = Att.get();
    blocks(0, 1) = Atn.get();
    blocks(1, 0) = Ant.get();
    blocks(1, 1) = Ann.get();
    opAr.reset(mfem::HypreParMatrixFromBlocks(blocks));

    mfem::Array<int> dbc_tdof_list;
    dbc_tdof_list.Append(nd_dbc_tdof_list);
    for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
    {
      dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
    }
    opAr->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  }

  // Build 2x2 block B matrix (Btt in (0,0), zero elsewhere).
  std::unique_ptr<mfem::HypreParMatrix> opBr;
  {
    mfem::SparseMatrix zero_sp(h1_size, h1_size);
    zero_sp.Finalize();
    mfem::HypreParMatrix zero_h1(h1_fespace.GetComm(), h1_fespace.GlobalTrueVSize(),
                                  h1_fespace.GlobalTrueVSize(),
                                  h1_fespace.Get().GetTrueDofOffsets(),
                                  h1_fespace.Get().GetTrueDofOffsets(), &zero_sp);

    mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
    blocks(0, 0) = Btt.get();
    blocks(0, 1) = nullptr;
    blocks(1, 0) = nullptr;
    blocks(1, 1) = &zero_h1;
    opBr.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Wrap as complex operators (real-valued system, no imaginary part for lossless).
  auto opA = std::make_unique<ComplexWrapperOperator>(std::move(opAr), nullptr);
  auto opB = std::make_unique<ComplexWrapperOperator>(std::move(opBr), nullptr);

  // Configure eigenvalue solver.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving for {:d} propagation mode(s)...\n", num_modes);

  std::unique_ptr<EigenvalueSolver> eigen;
#if defined(PALACE_WITH_SLEPC)
  {
    auto slepc = std::make_unique<slepc::SlepcEPSSolver>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc->SetProblemType(
        slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetNumModes(num_modes, std::max(2 * num_modes + 1, 20));
    slepc->SetTol(tol);
    slepc->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    eigen = std::move(slepc);
  }
#elif defined(PALACE_WITH_ARPACK)
  {
    auto arpack = std::make_unique<arpack::ArpackEPSSolver>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    arpack->SetNumModes(num_modes, std::max(2 * num_modes + 1, 20));
    arpack->SetTol(tol);
    arpack->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    eigen = std::move(arpack);
  }
#else
  MFEM_ABORT("ModeAnalysis solver requires SLEPc or ARPACK!");
#endif

  // Set up the KSP linear solver for the shift-and-invert spectral transformation.
  // Uses GMRES with a sparse direct preconditioner (same pattern as the wave port solver).
  std::unique_ptr<ComplexKspSolver> ksp;
  {
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(
        mesh.back()->GetComm(), iodata.problem.verbose);
    gmres->SetInitialGuess(false);
    gmres->SetRelTol(iodata.solver.linear.tol);
    gmres->SetMaxIter(iodata.solver.linear.max_it);
    gmres->SetRestartDim(iodata.solver.linear.max_it);

    // Choose a sparse direct solver as preconditioner.
#if defined(MFEM_USE_SUPERLU)
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        std::make_unique<SuperLUSolver>(mesh.back()->GetComm(),
                                        SymbolicFactorization::DEFAULT, false, true,
                                        iodata.problem.verbose - 1));
#elif defined(MFEM_USE_STRUMPACK)
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        std::make_unique<StrumpackSolver>(mesh.back()->GetComm(),
                                           SymbolicFactorization::DEFAULT,
                                           SparseCompression::NONE, 0.0, 0, 0, true,
                                           iodata.problem.verbose - 1));
#elif defined(MFEM_USE_MUMPS)
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        std::make_unique<MumpsSolver>(mesh.back()->GetComm(),
                                      mfem::MUMPSSolver::UNSYMMETRIC,
                                      SymbolicFactorization::DEFAULT, 0.0, true,
                                      iodata.problem.verbose - 1));
#else
    MFEM_ABORT("ModeAnalysis solver requires SuperLU_DIST, STRUMPACK, or MUMPS!");
    std::unique_ptr<MfemWrapperSolver<ComplexOperator>> pc;
#endif
    pc->SetSaveAssembled(false);
    pc->SetDropSmallEntries(false);
    ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
  }

  // Precondition with the real part of A.
  ComplexWrapperOperator opP(opA->Real(), nullptr);
  ksp->SetOperators(*opA, opP);
  eigen->SetLinearSolver(*ksp);

  // Solve the generalized eigenvalue problem B x = lambda A x (shift-and-invert).
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);
  int num_conv = eigen->Solve();
  Mpi::Print(" Found {:d} converged eigenvalue(s)\n\n", num_conv);

  // Extract and report mode properties.
  BlockTimer bt2(Timer::POSTPRO);
  const double kc = 1.0 / iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);

  Mpi::Print(" {:>5s}  {:>22s}  {:>22s}  {:>16s}  {:>16s}\n", "Mode", "Re{kn} (1/m)",
             "Im{kn} (1/m)", "Re{n_eff}", "Im{n_eff}");
  Mpi::Print(" {}\n", std::string(88, '-'));

  for (int i = 0; i < std::min(num_conv, num_modes); i++)
  {
    std::complex<double> lambda = eigen->GetEigenvalue(i);

    // The extracted eigenvalue is lambda = 1 / (-kn^2 - sigma).
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);

    // Effective index: n_eff = kn / k0, where k0 = omega/c is the free-space wavenumber.
    // In Palace's nondimensional units: k0 = omega (since c = 1/sqrt(mu0*eps0) = 1 nondim).
    std::complex<double> n_eff = kn / omega;

    // Dimensionalize kn.
    std::complex<double> kn_dim = kn * kc;

    Mpi::Print(" {:5d}  {:>+22.12e}  {:>+22.12e}  {:>+16.8e}  {:>+16.8e}\n", i + 1,
               kn_dim.real(), kn_dim.imag(), n_eff.real(), n_eff.imag());
  }
  Mpi::Print("\n");

  // TODO: Save mode field patterns to Paraview (Phase 1b).
  // TODO: Compute Z0 from V/I integrals, then L and C per unit length (Phase 2).

  ErrorIndicator indicator;
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
