// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_SOLVER_HPP
#define PALACE_SCHEMA_TYPES_SOLVER_HPP

// Mirrors palace::config::SolverData and its arm types. Descriptions match
// PR 716's scripts/schema/config/solver.json.
//
// Phase 1 notes:
//
//   - `Samples` is a TaggedUnion<"Type", Point, Linear, Log>; PR 716 uses
//     `oneOf`. We opt into oneOf via `palace::schema::utils::schema_composition` at the bottom
//     of this header.
//
//   - Enums like `AdaptiveCircuitSynthesisDomainOrthogonalization` are new
//     in PR 716. We add new enum types here rather than reusing Palace's
//     labels.hpp enums (Phase 2 will consolidate).
//
//   - `InitialGuess` in PR 716 is a boolean; Palace's current runtime stores
//     an int sentinel (-1 = default, 0/1 otherwise). We follow PR 716's
//     boolean wire format for Phase 1.

#include <optional>
#include <string>
#include <vector>

#include <rfl.hpp>

#include "schema/utils/annotations.hpp"
#include "common.hpp"

namespace palace::schema {

// New enums introduced by PR 716 (not present in labels.hpp).

enum class AdaptiveCircuitSynthesisDomainOrthogonalization {
    Energy,
    FEBasisIdentity,
    SpaceOverlap,
};

// --- Samples (tagged union arms) -------------------------------------------

struct SamplesPoint {
    PALACE_SCHEMA_TAG(Type, "Point", "Explicit list of frequency sample points.");

    PALACE_SCHEMA_DESC(Freq, "Explicit frequencies to sample, GHz.",
             std::vector<double>) = {};

    PALACE_SCHEMA_DESC(SaveStep,
             "Save fields every N steps within this sample. `0` disables saving.",
             std::optional<int>) = 0;

    PALACE_SCHEMA_DESC(AddToPROM,
             "Force inclusion of these points in the PROM for adaptive sweep "
             "(primarily a debugging tool).",
             std::optional<bool>) = false;
};

struct SamplesLinear {
    PALACE_SCHEMA_TAG(Type, "Linear", "Linearly-spaced frequency samples.");

    PALACE_SCHEMA_DESC(MinFreq, "Lower bound, GHz.", palace::schema::utils::Min<double, 0.0>) = 0.0;

    PALACE_SCHEMA_DESC(MaxFreq, "Upper bound, GHz.", palace::schema::utils::Min<double, 0.0>) = 0.0;

    PALACE_SCHEMA_DESC(FreqStep,
             "Step size, GHz. Mutually exclusive with `\"NSample\"`.",
             std::optional<double>) = std::nullopt;

    PALACE_SCHEMA_DESC(NSample,
             "Number of samples. Mutually exclusive with `\"FreqStep\"`.",
             std::optional<int>) = std::nullopt;

    PALACE_SCHEMA_DESC(SaveStep, "Save fields every N steps. `0` disables saving.",
             std::optional<int>) = 0;

    PALACE_SCHEMA_DESC(AddToPROM, "Force inclusion in the PROM for adaptive sweep.",
             std::optional<bool>) = false;
};

struct SamplesLog {
    PALACE_SCHEMA_TAG(Type, "Log", "Logarithmically-spaced frequency samples.");

    PALACE_SCHEMA_DESC(MinFreq, "Lower bound, GHz.", palace::schema::utils::XMin<double, 0.0>) = 1.0;

    PALACE_SCHEMA_DESC(MaxFreq, "Upper bound, GHz.", palace::schema::utils::XMin<double, 0.0>) = 1.0;

    PALACE_SCHEMA_DESC(NSample, "Number of samples.", palace::schema::utils::Min<int, 1>) = 1;

    PALACE_SCHEMA_DESC(SaveStep, "Save fields every N steps. `0` disables saving.",
             std::optional<int>) = 0;

    PALACE_SCHEMA_DESC(AddToPROM, "Force inclusion in the PROM for adaptive sweep.",
             std::optional<bool>) = false;
};

using Sample =
    rfl::TaggedUnion<"Type", SamplesPoint, SamplesLinear, SamplesLog>;

// --- DrivenSolverData ------------------------------------------------------

struct DrivenSolverData {
    // x-palace-deprecated
    PALACE_SCHEMA_DESC(MinFreq,
             "Lower bound of the frequency sweep interval, GHz. Deprecated: use "
             "[`Linear Samples`](@ref config-solver-driven-samples-linear) "
             "interface instead.",
             std::optional<double>) = std::nullopt;

    // x-palace-deprecated
    PALACE_SCHEMA_DESC(MaxFreq,
             "Upper bound of the frequency sweep interval, GHz. Deprecated: use "
             "[`Linear Samples`](@ref config-solver-driven-samples-linear) "
             "interface instead.",
             std::optional<double>) = std::nullopt;

    // x-palace-deprecated
    PALACE_SCHEMA_DESC(FreqStep,
             "Frequency step size for the frequency sweep, GHz. Deprecated: use "
             "[`Linear Samples`](@ref config-solver-driven-samples-linear) "
             "interface instead.",
             std::optional<double>) = std::nullopt;

    // x-palace-deprecated
    PALACE_SCHEMA_DESC(SaveStep,
             "Controls how often, in number of frequency steps, to save computed "
             "fields to disk for [visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization). Files are saved in the "
             "`paraview/` (and/or `gridfunction/`) directory under "
             "[`/Problem/Output`](@ref config-problem-output). Deprecated: use "
             "[`Linear Samples`](@ref config-solver-driven-samples-linear) "
             "interface instead.",
             std::optional<int>) = 0;

    PALACE_SCHEMA_DESC(Samples,
             "Array of frequency sample specifications. Combined with "
             "`\"MinFreq\"`/`\"MaxFreq\"`/`\"FreqStep\"` to form a sorted, "
             "unique set of samples.",
             std::vector<Sample>) = {};

    PALACE_SCHEMA_DESC(Save,
             "Additional frequencies at which to save computed fields to disk "
             "for [visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization), in addition to those "
             "specified by `\"SaveStep\"`. Files are saved in the `paraview/` "
             "(and/or `gridfunction/`) directory under [`/Problem/Output`]"
             "(@ref config-problem-output).",
             std::vector<double>) = {};

    PALACE_SCHEMA_DESC(Restart,
             "1-based sample index from which to restart a partial frequency "
             "sweep (i.e. `\"Restart\": x` starts from the *x*-th sample of the "
             "combined sample set). Not valid for adaptive sweep.",
             palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(AdaptiveTol,
             "Relative error convergence tolerance for adaptive fast frequency "
             "sweep. A value of `0` disables adaptive sweep and the full-order "
             "model is solved at each frequency step. A positive value ensures "
             "the reduced-order model is reliable relative to the full-order "
             "model in the frequency band of interest.",
             palace::schema::utils::Min<double, 0.0>) = 0.0;

    PALACE_SCHEMA_DESC(AdaptiveMaxSamples,
             "Maximum number of frequency samples used to construct the "
             "reduced-order model for adaptive sweep, if the tolerance "
             "(`\"AdaptiveTol\"`) is not met first. In simulations with "
             "multiple excitations, this is the maximum per excitation.",
             palace::schema::utils::XMin<int, 0>) = 20;

    PALACE_SCHEMA_DESC(AdaptiveConvergenceMemory,
             "Number of consecutive samples satisfying the error tolerance "
             "required to declare convergence of the adaptive sampling "
             "algorithm.",
             palace::schema::utils::XMin<int, 0>) = 2;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(AdaptiveGSOrthogonalization,
             "Gram-Schmidt variant for orthogonalizing the adaptive "
             "reduced-order model basis. Same options as "
             "[`/Linear/GSOrthogonalization`](@ref "
             "config-solver-linear-gsorthogonalization).",
             Orthogonalization) = Orthogonalization::CGS2;

    PALACE_SCHEMA_DESC(AdaptiveCircuitSynthesis,
             "Use the adaptive reduced-order model to print synthesized "
             "circuit-like matrices (L⁻¹, R⁻¹, C). Requires adaptive sweep to "
             "be enabled, all `LumpedPort` fields to be orthogonal, and only "
             "LRC-type frequency dependence (no `WavePort`, `Conductivity`, or "
             "second-order farfield BCs).",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(AdaptiveCircuitSynthesisDomainOrthogonalization,
             "Weight matrix type for domain orthogonalization when building "
             "synthesized circuit matrices.",
             AdaptiveCircuitSynthesisDomainOrthogonalization) =
        AdaptiveCircuitSynthesisDomainOrthogonalization::Energy;
};

// --- EigenSolverData -------------------------------------------------------

struct EigenSolverData {
    PALACE_SCHEMA_DESC(Target,
             "(Nonzero) frequency target above which to search for "
             "eigenvalues, GHz.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0;

    PALACE_SCHEMA_DESC(Tol,
             "Relative convergence tolerance for the eigenvalue solver.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-6;

    PALACE_SCHEMA_DESC(MaxIts,
             "Maximum number of iterations for the iterative eigenvalue solver. "
             "A value less than 1 uses the solver default.",
             palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(MaxSize,
             "Maximum subspace dimension for the eigenvalue solver. A value "
             "less than 1 uses the solver default.",
             palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(N, "Number of eigenvalues to compute.", palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(Save,
             "Number of computed field modes to save to disk for "
             "[visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization). Files are saved in "
             "the `paraview/` (and/or `gridfunction/`) directory under "
             "[`/Problem/Output`](@ref config-problem-output).",
             int) = 0;

    PALACE_SCHEMA_DESC(Type, "Specifies the eigenvalue solver backend.",
             EigenSolverBackend) = EigenSolverBackend::Default;

    PALACE_SCHEMA_DESC(NonlinearType,
             "Specifies the nonlinear eigenvalue solver for problems with "
             "frequency-dependent boundary conditions.",
             NonlinearEigenSolver) = NonlinearEigenSolver::Hybrid;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(TargetUpper,
             "Upper frequency bound for the eigenvalue search, GHz. Only used "
             "for nonlinear problems. A value less than or equal to zero uses "
             "`3 × Target` automatically. An inaccurate upper bound can "
             "negatively affect convergence of the nonlinear eigensolver.",
             double) = -1.0;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(PEPLinear,
             "Linearize the polynomial eigenvalue problem before solving.",
             bool) = true;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(Scaling, "Enable scaling of the eigenvalue problem.",
             bool) = true;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(StartVector,
             "Use a random start vector for the eigensolver.", bool) = true;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(StartVectorConstant,
             "Use a constant start vector instead of a random one.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MassOrthogonal,
             "Orthogonalize eigenvectors with respect to the mass matrix.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(RefineNonlinear,
             "Refine nonlinear eigenvalue solutions after the initial solve.",
             bool) = true;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(LinearTol,
             "Tolerance for the inner linear solve within the nonlinear "
             "eigensolver.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-3;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(PreconditionerLag,
             "Number of eigensolver iterations between preconditioner updates.",
             palace::schema::utils::Min<int, 0>) = 10;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(PreconditionerLagTol,
             "Residual tolerance below which preconditioner updates are skipped.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-4;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MaxRestart, "Maximum number of restarts for the eigensolver.",
             palace::schema::utils::Min<int, 0>) = 2;
};

// --- TransientSolverData ---------------------------------------------------

struct TransientSolverData {
    PALACE_SCHEMA_DESC(Type, "Time integration scheme for the second-order ODE system.",
             TimeSteppingScheme) = TimeSteppingScheme::Default;

    PALACE_SCHEMA_DESC(Excitation, "Controls the time dependence of the source excitation.",
             Excitation) = Excitation::Sinusoidal;

    PALACE_SCHEMA_DESC(ExcitationFreq,
             "Center frequency for harmonic source excitations, GHz. Only "
             "relevant when `\"Excitation\"` is `\"Sinusoidal\"`, "
             "`\"Gaussian\"`, `\"DifferentiatedGaussian\"`, or "
             "`\"ModulatedGaussian\"`.",
             double) = 0.0;

    PALACE_SCHEMA_DESC(ExcitationWidth,
             "Pulse width for Gaussian-type source excitations, ns. Only "
             "relevant when `\"Excitation\"` is `\"Gaussian\"`, "
             "`\"DifferentiatedGaussian\"`, `\"ModulatedGaussian\"`, or "
             "`\"SmoothStep\"`.",
             double) = 0.0;

    PALACE_SCHEMA_DESC(MaxTime, "End of simulation time interval, ns.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0;

    PALACE_SCHEMA_DESC(TimeStep, "Uniform time step size, ns.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0e-2;

    PALACE_SCHEMA_DESC(SaveStep,
             "Controls how often, in number of time steps, to save computed "
             "fields to disk for [visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization). Files are saved in "
             "the `paraview/` (and/or `gridfunction/`) directory under "
             "[`/Problem/Output`](@ref config-problem-output).",
             int) = 0;

    PALACE_SCHEMA_DESC(Order,
             "Order of adaptive Runge-Kutta integrators or maximum multistep "
             "method order, must be within `[2, 5]`. Only relevant when "
             "`\"Type\"` is `\"ARKODE\"` or `\"CVODE\"`.",
             palace::schema::utils::Closed<int, 2, 5>) = 2;

    PALACE_SCHEMA_DESC(RelTol,
             "Relative tolerance for adaptive time-stepping. Only relevant when "
             "`\"Type\"` is `\"ARKODE\"` or `\"CVODE\"`.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0e-4;

    PALACE_SCHEMA_DESC(AbsTol,
             "Absolute tolerance for adaptive time-stepping. Only relevant when "
             "`\"Type\"` is `\"ARKODE\"` or `\"CVODE\"`.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0e-9;
};

// --- Electrostatic / Magnetostatic -----------------------------------------

struct ElectrostaticSolverData {
    PALACE_SCHEMA_DESC(Save,
             "Number of computed electric field solutions to save to disk for "
             "[visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization), ordered by the "
             "entries in the computed capacitance matrix. Files are saved in "
             "the `paraview/` (and/or `gridfunction/`) directory under "
             "[`/Problem/Output`](@ref config-problem-output).",
             palace::schema::utils::Min<int, 0>) = 0;
};

struct MagnetostaticSolverData {
    PALACE_SCHEMA_DESC(Save,
             "Number of computed magnetic field solutions to save to disk for "
             "[visualization with ParaView]"
             "(../guide/postprocessing.md#Visualization), ordered by the "
             "entries in the computed inductance matrix. Files are saved in "
             "the `paraview/` (and/or `gridfunction/`) directory under "
             "[`/Problem/Output`](@ref config-problem-output).",
             palace::schema::utils::Min<int, 0>) = 0;
};

// --- LinearSolverData ------------------------------------------------------

struct LinearSolverData {
    PALACE_SCHEMA_DESC(Type,
             "Specifies the solver used for preconditioning the linear system.",
             LinearSolver) = LinearSolver::Default;

    PALACE_SCHEMA_DESC(KSPType, "Specifies the iterative Krylov subspace solver.",
             KrylovSolver) = KrylovSolver::Default;

    PALACE_SCHEMA_DESC(Tol,
             "Relative residual convergence tolerance for the iterative linear "
             "solver.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-6;

    PALACE_SCHEMA_DESC(MaxIts,
             "Maximum number of iterations for the iterative linear solver.",
             palace::schema::utils::XMin<int, 0>) = 100;

    PALACE_SCHEMA_DESC(MaxSize,
             "Maximum Krylov space size for GMRES/FGMRES. A value less than 1 "
             "defaults to `\"MaxIts\"`.",
             palace::schema::utils::XMin<int, 0>) = 1;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(InitialGuess,
             "Use the previous solution as an initial guess for the iterative "
             "solver.",
             bool) = true;

    PALACE_SCHEMA_DESC(MGMaxLevels,
             "When greater than 1, enables [geometric multigrid preconditioning]"
             "(https://en.wikipedia.org/wiki/Multigrid_method) using p- and "
             "h-multigrid coarsening. The solver specified by `\"Type\"` is "
             "used on the coarsest level; fine levels use Chebyshev smoothing.",
             palace::schema::utils::Min<int, 1>) = 100;

    PALACE_SCHEMA_DESC(MGCoarsenType,
             "Coarsening strategy for constructing p-multigrid levels.",
             MultigridCoarsening) = MultigridCoarsening::Logarithmic;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MGUseMesh, "Use the mesh hierarchy for h-multigrid coarsening.",
             bool) = true;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MGAuxiliarySmoother,
             "Use an auxiliary space smoother on multigrid levels.",
             bool) = true;

    PALACE_SCHEMA_DESC(MGCycleIts,
             "Number of V-cycle iterations per preconditioner application for "
             "multigrid preconditioners (when geometric multigrid is enabled "
             "or `\"Type\"` is `\"AMS\"` or `\"BoomerAMG\"`). A value less "
             "than 1 defaults to 2 for AMS frequency domain problems, 1 "
             "otherwise.",
             palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(MGSmoothIts,
             "Number of pre- and post-smooth iterations for multigrid "
             "preconditioners, when the geometric multigrid preconditioner is "
             "enabled.",
             palace::schema::utils::XMin<int, 0>) = 1;

    PALACE_SCHEMA_DESC(MGSmoothOrder,
             "Polynomial smoothing order for geometric multigrid. A value less "
             "than 1 defaults to `max(2 × solution order, 4)`.",
             palace::schema::utils::XMin<int, 0>) = 1;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MGSmoothEigScaleMax,
             "Maximum eigenvalue scale for the Chebyshev smoother.",
             palace::schema::utils::XMin<double, 0.0>) = 1.0;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MGSmoothEigScaleMin,
             "Minimum eigenvalue scale for the Chebyshev smoother.",
             palace::schema::utils::Min<double, 0.0>) = 0.0;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(MGSmoothChebyshev4th,
             "Use 4th-kind Chebyshev polynomial smoother.", bool) = true;

    PALACE_SCHEMA_DESC(PCMatReal,
             "When `true`, builds the preconditioner for frequency domain "
             "problems using a real-valued approximation of the system matrix. "
             "Always used on the coarsest multigrid level regardless of this "
             "setting.",
             bool) = false;

    PALACE_SCHEMA_DESC(PCMatShifted,
             "When `true`, builds the preconditioner using a positive-definite "
             "approximation by flipping the sign of the mass matrix "
             "contribution. Can improve performance at high frequencies "
             "relative to the lowest nonzero eigenfrequencies of the model.",
             bool) = false;

    PALACE_SCHEMA_DESC(ComplexCoarseSolve,
             "When `true`, uses the true complex-valued system matrix for the "
             "coarse-level solver. When `false`, uses the real-valued "
             "approximation.",
             bool) = false;

    PALACE_SCHEMA_DESC(DropSmallEntries,
             "When `true`, drops entries smaller than double-precision machine "
             "epsilon from the system matrix used in the sparse direct solver.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(ReorderingReuse,
             "Reuse the matrix reordering from the previous solve.",
             bool) = true;

    PALACE_SCHEMA_DESC(PCSide,
             "Preconditioning side. Not all options are available for all "
             "iterative solver choices.",
             PreconditionerSide) = PreconditionerSide::Default;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(ColumnOrdering,
             "Column ordering algorithm for sparse direct solvers: "
             "`\"METIS\"`, `\"ParMETIS\"`, `\"Scotch\"`, `\"PTScotch\"`, "
             "`\"PORD\"`, `\"AMD\"`, `\"RCM\"`, `\"Default\"`.",
             SymbolicFactorization) = SymbolicFactorization::Default;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(STRUMPACKCompressionType,
             "STRUMPACK compression type: `\"None\"`, `\"BLR\"`, `\"HSS\"`, "
             "`\"HODLR\"`, `\"ZFP\"`, `\"BLR-HODLR\"`, `\"ZFP-BLR-HODLR\"`.",
             SparseCompression) = SparseCompression::None;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(STRUMPACKCompressionTol,
             "Tolerance for STRUMPACK lossy compression.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-3;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(STRUMPACKLossyPrecision,
             "Precision bits for STRUMPACK ZFP lossy compression.",
             palace::schema::utils::Min<int, 0>) = 16;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(STRUMPACKButterflyLevels,
             "Number of butterfly levels for STRUMPACK HODLR compression.",
             palace::schema::utils::Min<int, 0>) = 1;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(SuperLU3DCommunicator,
             "Use a 3D process grid communicator for SuperLU_DIST.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(AMSVectorInterpolation, "Use vector interpolation in AMS.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(AMSSingularOperator,
             "Indicate to AMS that the operator is singular (e.g. "
             "magnetostatics with no essential boundary conditions).",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(AMGAggressiveCoarsening,
             "Enable aggressive coarsening in BoomerAMG.", bool) = false;

    PALACE_SCHEMA_DESC(AMSMaxIts,
             "Number of AMS iterations per preconditioner application when "
             "geometric multigrid is enabled (`\"MGMaxLevels\"` > 1). A value "
             "less than 1 defaults to the solution order ([`/Solver/Order`]"
             "(@ref config-solver-order)).",
             palace::schema::utils::Min<int, 0>) = 0;

    PALACE_SCHEMA_DESC(DivFreeTol,
             "Relative tolerance for divergence-free cleaning used in the "
             "eigenmode simulation type. Ignored when a nonzero Floquet wave "
             "vector is specified or when a nonzero [LondonDepth]"
             "(@ref config-domains-materials-londondepth) is used.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-12;

    PALACE_SCHEMA_DESC(DivFreeMaxIts,
             "Maximum number of iterations for divergence-free cleaning. "
             "Ignored under the same conditions as `\"DivFreeTol\"`.",
             palace::schema::utils::Min<int, 0>) = 1000;

    PALACE_SCHEMA_DESC(EstimatorTol,
             "Relative tolerance for the flux projection solve used in the "
             "error estimate calculation.",
             palace::schema::utils::Min<double, 0.0>) = 1.0e-6;

    PALACE_SCHEMA_DESC(EstimatorMaxIts,
             "Maximum number of iterations for the flux projection solve used "
             "in the error estimate calculation.",
             palace::schema::utils::Min<int, 0>) = 10000;

    PALACE_SCHEMA_DESC(EstimatorMG,
             "When `true`, uses a multigrid preconditioner with AMG coarse "
             "solve for the error estimator linear solver instead of Jacobi.",
             bool) = false;

    PALACE_SCHEMA_DESC(GSOrthogonalization,
             "Gram-Schmidt variant used to orthogonalize vectors in Krylov "
             "subspace methods and other parts of the solver.",
             Orthogonalization) = Orthogonalization::MGS;
};

// --- Top-level SolverData --------------------------------------------------

struct SolverData {
    PALACE_SCHEMA_DESC(Order,
             "Finite element order (degree). Arbitrary high-order spaces are "
             "supported.",
             palace::schema::utils::Min<int, 1>) = 1;

    PALACE_SCHEMA_DESC(PartialAssemblyOrder,
             "Order at which to switch from full assembly of finite element "
             "operators to [partial assembly]"
             "(https://mfem.org/howto/assembly_levels/). Setting to `1` fully "
             "activates partial assembly on all levels; a large value (greater "
             "than [Order](@ref config-solver-order)) results in fully "
             "assembled sparse matrix operators.",
             palace::schema::utils::Min<int, 1>) = 1;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(QuadratureOrderJacobian,
             "Use the Jacobian-based quadrature order instead of the default.",
             bool) = false;

    // x-palace-advanced
    PALACE_SCHEMA_DESC(QuadratureOrderExtra,
             "Extra quadrature order added on top of the default.",
             int) = 0;

    PALACE_SCHEMA_DESC(Device,
             "The runtime device configuration passed to [MFEM]"
             "(https://mfem.org/howto/assembly_levels/) to activate different "
             "computation backends. When *Palace* is built with OpenMP support "
             "(`PALACE_WITH_OPENMP=ON`), `omp` is automatically added to the "
             "MFEM device list.",
             Device) = Device::CPU;

    PALACE_SCHEMA_DESC(Backend,
             "Specifies the [libCEED backend]"
             "(https://libceed.org/en/latest/gettingstarted/#backends) to use. "
             "If empty, a suitable default is selected based on [Device]"
             "(@ref config-solver-device).",
             std::string) = "";

    PALACE_SCHEMA_DESC(Driven,
             "Configuration for the frequency domain driven solver. Only "
             "relevant when [`/Problem/Type`](@ref config-problem-type) is "
             "`\"Driven\"`.",
             DrivenSolverData) = {};

    PALACE_SCHEMA_DESC(Eigenmode,
             "Configuration for the eigenvalue solver. Only relevant when "
             "[`/Problem/Type`](@ref config-problem-type) is `\"Eigenmode\"`.",
             EigenSolverData) = {};

    PALACE_SCHEMA_DESC(Electrostatic,
             "Configuration for the electrostatic solver. Only relevant when "
             "[`/Problem/Type`](@ref config-problem-type) is "
             "`\"Electrostatic\"`.",
             ElectrostaticSolverData) = {};

    PALACE_SCHEMA_DESC(Magnetostatic,
             "Configuration for the magnetostatic solver. Only relevant when "
             "[`/Problem/Type`](@ref config-problem-type) is "
             "`\"Magnetostatic\"`.",
             MagnetostaticSolverData) = {};

    PALACE_SCHEMA_DESC(Transient,
             "Configuration for the time domain driven solver. Only relevant "
             "when [`/Problem/Type`](@ref config-problem-type) is "
             "`\"Transient\"`. Simulations always start from rest at *t* = 0.",
             TransientSolverData) = {};

    PALACE_SCHEMA_DESC(Linear,
             "Configuration for the linear solver used by all simulation types.",
             LinearSolverData) = {};
};

}  // namespace palace::schema

// Opt into `oneOf` for the Samples union so the emitted schema matches PR
// 716's shape (reflect-cpp's default is `anyOf`).
template <>
struct palace::schema::utils::schema_composition<palace::schema::Sample> {
    static constexpr auto value = palace::schema::utils::Compose::OneOf;
};

#endif  // PALACE_SCHEMA_TYPES_SOLVER_HPP
