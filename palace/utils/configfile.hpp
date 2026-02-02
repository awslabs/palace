// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CONFIG_FILE_HPP
#define PALACE_UTILS_CONFIG_FILE_HPP

#include <array>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "labels.hpp"

namespace palace::config
{

using json = nlohmann::json;

//
// Data structures for storing configuration file data.
//
namespace internal
{

// An ElementData consists of a list of attributes making up a single element of a
// potentially multielement boundary, and a direction and/or a normal defining the incident
// field. These are used for lumped ports, terminals, surface currents, and other boundary
// postprocessing objects.
struct ElementData
{
  // Vector defining the direction for this port. In a Cartesian system, "X", "Y", and "Z"
  // map to (1,0,0), (0,1,0), and (0,0,1), respectively.
  std::array<double, 3> direction{{0.0, 0.0, 0.0}};

  CoordinateSystem coordinate_system = CoordinateSystem::CARTESIAN;

  // List of boundary attributes for this element.
  std::vector<int> attributes = {};
};

}  // namespace internal

// Problem & Model Config.

struct OutputFormatsData
{
public:
  // Enable Paraview output format.
  bool paraview = true;

  // Enable MFEM GLVis grid function output format.
  bool gridfunction = false;
};

struct ProblemData
{
public:
  // Simulation type.
  ProblemType type = ProblemType::DRIVEN;

  // Level of printing.
  int verbose = 1;

  // Output path for storing results.
  std::string output = "";

  // Output formats configuration.
  OutputFormatsData output_formats = {};

  ProblemData() = default;
  ProblemData(const json &problem);
};

struct BoxRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Region bounding box limits [m].
  std::array<double, 3> bbmin{{0.0, 0.0, 0.0}}, bbmax{{0.0, 0.0, 0.0}};
};

struct SphereRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Sphere radius [m].
  double r = 0.0;

  // Sphere center [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};
};

struct RefinementData
{
public:
  // Non-dimensional tolerance used to specify convergence of adaptive mesh refinement.
  double tol = 1.0e-2;

  // Maximum number of iterations to perform during adaptive mesh refinement.
  int max_it = 0;

  // If a refinement results in a greater number of DOFs than this value, no future
  // refinement will be allowed.
  int max_size = 0;

  // Whether or not to perform nonconformal adaptation.
  bool nonconformal = true;

  // Maximum difference in nonconformal refinements between two adjacent elements. Zero
  // implies there is no constraint on local nonconformity.
  int max_nc_levels = 1;

  // Dörfler update fraction. The set of marked elements is the minimum set that contains
  // update_fraction of the total error.
  double update_fraction = 0.7;

  // Maximum allowable ratio of number of elements across processors before rebalancing is
  // performed.
  double maximum_imbalance = 1.1;

  // Whether to save off results of each adaptation iteration as a subfolder within the post
  // processing directory.
  bool save_adapt_iterations = true;

  // Whether to write a (serial) mesh to file after mesh modification during AMR.
  bool save_adapt_mesh = false;

  // Parallel uniform mesh refinement levels.
  int uniform_ref_levels = 0;

  // Serial uniform mesh refinement levels.
  int ser_uniform_ref_levels = 0;

private:
  // Refinement data for mesh regions.
  std::vector<BoxRefinementData> box_list = {};
  std::vector<SphereRefinementData> sphere_list = {};

public:
  auto &GetBox(int i) { return box_list[i]; }
  const auto &GetBoxes() const { return box_list; }
  auto &GetBoxes() { return box_list; }

  auto &GetSphere(int i) { return sphere_list[i]; }
  const auto &GetSpheres() const { return sphere_list; }
  auto &GetSpheres() { return sphere_list; }

  RefinementData() = default;
  RefinementData(const json &refinement);
};

struct ModelData
{
public:
  // Mesh file.
  std::string mesh = "";

  // Mesh length unit and optional characteristic length scale for nondimensionalization
  // [m].
  double L0 = 1.0e-6;
  double Lc = -1.0;

  // Remove high-order curvature information from the mesh.
  bool remove_curvature = false;

  // Convert mesh to simplex elements.
  bool make_simplex = false;

  // Convert mesh to hexahedral elements (using tet-to-hex algorithm).
  bool make_hex = false;

  // Reorder elements based on spatial location after loading the serial mesh, which can
  // potentially increase memory coherency.
  bool reorder_elements = false;

  // Remove elements (along with any associated unattached boundary elements) from the mesh
  // which do not have any material properties specified.
  bool clean_unused_elements = true;

  // Split, or "crack", boundary elements lying on internal boundaries to decouple the
  // elements on either side.
  bool crack_bdr_elements = true;

  // When required, refine elements neighboring a split or crack in order to enable the
  // decoupling.
  bool refine_crack_elements = true;

  // Factor for displacing duplicated interior boundary elements, usually added just for
  // visualization.
  double crack_displ_factor = 1.0e-12;

  // Add new boundary elements for faces are on the computational domain boundary or which
  // have attached elements on either side with different domain attributes.
  bool add_bdr_elements = true;

  // Export mesh after pre-processing but before cracking.
  bool export_prerefined_mesh = false;

  // Call MFEM's ReorientTetMesh as a check of mesh orientation after partitioning.
  bool reorient_tet_mesh = false;

  // Partitioning file (if specified, does not compute a new partitioning).
  std::string partitioning = "";

  // Object controlling mesh refinement.
  RefinementData refinement = {};

  ModelData() = default;
  ModelData(const json &model);
};

// Domain Config.

// Store symmetric matrix data as set of outer products: Σᵢ sᵢ * vᵢ *  vᵢᵀ.
template <std::size_t N>
struct SymmetricMatrixData
{
public:
  std::array<double, N> s;
  std::array<std::array<double, N>, N> v;

  SymmetricMatrixData(double diag)
  {
    s.fill(diag);
    std::size_t i = 0;
    for (auto &x : v)
    {
      x.fill(0.0);
      x[i++] = 1.0;
    }
  }
};

struct MaterialData
{
public:
  // Relative permeability.
  SymmetricMatrixData<3> mu_r = 1.0;

  // Relative permittivity.
  SymmetricMatrixData<3> epsilon_r = 1.0;

  // Loss tangent.
  SymmetricMatrixData<3> tandelta = 0.0;

  // Conductivity [S/m].
  SymmetricMatrixData<3> sigma = 0.0;

  // London penetration depth [m].
  double lambda_L = 0.0;

  // List of domain attributes for this material.
  std::vector<int> attributes = {};

  MaterialData() = default;
  MaterialData(const json &domain);
};

struct DomainEnergyData
{
public:
  // List of domain attributes for this domain postprocessing index.
  std::vector<int> attributes = {};

  DomainEnergyData() = default;
  DomainEnergyData(const json &domain);
};

struct ProbeData
{
public:
  // Physical space coordinates for the probe location [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};

  ProbeData() = default;
  ProbeData(const json &probe);
};

struct DomainPostData
{
public:
  // List of all postprocessing domain attributes.
  std::vector<int> attributes = {};

  // Domain postprocessing objects.
  std::map<int, DomainEnergyData> energy = {};
  std::map<int, ProbeData> probe = {};

  DomainPostData() = default;
  DomainPostData(const json &postpro);
};

struct CurrentDipoleData
{
public:
  // Current dipole direction (normalized unit vector).
  std::array<double, 3> direction{{0.0, 0.0, 0.0}};

  // Current dipole moment magnitude [A·m].
  double moment = 0.0;

  // Current dipole center position.
  std::array<double, 3> center{{0.0, 0.0, 0.0}};

  CurrentDipoleData() = default;
  CurrentDipoleData(const json &source);
};

struct DomainData
{
public:
  // List of all domain attributes (excluding postprocessing).
  std::vector<int> attributes = {};

  // Domain objects.
  std::vector<MaterialData> materials = {};
  std::map<int, CurrentDipoleData> current_dipole = {};
  DomainPostData postpro = {};

  DomainData() = default;
  DomainData(const json &domains);
};

// Boundary Configuration.

struct PecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  PecBoundaryData() = default;
  PecBoundaryData(const json &pec);
};

struct PmcBoundaryData
{
public:
  // List of boundary attributes with PMC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  PmcBoundaryData() = default;
  PmcBoundaryData(const json &pmc);
};

struct WavePortPecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions for wave ports.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  WavePortPecBoundaryData() = default;
  WavePortPecBoundaryData(const json &auxpec);
};

struct FarfieldBoundaryData
{
public:
  // Approximation order for farfield ABC.
  int order = 1;

  // List of boundary attributes with farfield absorbing boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  FarfieldBoundaryData() = default;
  FarfieldBoundaryData(const json &absorbing);
};

struct ConductivityData
{
public:
  // Electrical conductivity of the conductor [S/m].
  double sigma = 0.0;

  // Conductor relative permeability.
  double mu_r = 1.0;

  // Optional conductor thickness [m].
  double h = 0.0;

  // Optional flag for an external boundary surface, relevant for the thickness correction.
  bool external = false;

  // List of boundary attributes for this surface conductivity boundary condition.
  std::vector<int> attributes = {};

  ConductivityData() = default;
  ConductivityData(const json &boundary);
};

struct ImpedanceData
{
public:
  // Boundary surface resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double Rs = 0.0;
  double Ls = 0.0;
  double Cs = 0.0;

  // List of boundary attributes for this impedance boundary condition.
  std::vector<int> attributes = {};

  ImpedanceData() = default;
  ImpedanceData(const json &boundary);
};

struct LumpedPortData
{
public:
  // Port circuit resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double R = 0.0;
  double L = 0.0;
  double C = 0.0;

  // Port surface resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double Rs = 0.0;
  double Ls = 0.0;
  double Cs = 0.0;

  // Input excitation for driven & transient solver:
  // - Wave/Lumped ports with same index are excited together.
  // - 1-based index if excited; 0 if not excited.
  int excitation = 0;

  // Flag for boundary damping term in driven and transient simulations.
  bool active = true;

  // For each lumped port index, each element contains a list of attributes making up a
  // single element of a potentially multielement lumped port.
  std::vector<internal::ElementData> elements = {};

  LumpedPortData() = default;
  LumpedPortData(const json &port);
};

struct TerminalData
{
public:
  // List of boundary attributes for this terminal.
  std::vector<int> attributes = {};

  TerminalData() = default;
  TerminalData(const json &terminal);
};

struct PeriodicData
{
public:
  // Vector defining the affine transformation matrix for this periodic boundary condition.
  std::array<double, 16> affine_transform = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // List of boundary donor attributes for this periodic boundary condition.
  std::vector<int> donor_attributes = {};

  // List of boundary receiver attributes for this periodic boundary condition.
  std::vector<int> receiver_attributes = {};
};

struct PeriodicBoundaryData
{
public:
  // Vector of periodic boundary pairs.
  std::vector<PeriodicData> boundary_pairs = {};

  // Floquet/Bloch wavevector specifying the phase delay in the X/Y/Z directions.
  std::array<double, 3> wave_vector = {0.0, 0.0, 0.0};

  PeriodicBoundaryData() = default;
  PeriodicBoundaryData(const json &periodic);
};

struct WavePortData
{
public:
  // Mode index for the numeric wave port.
  int mode_idx = 1;

  // Port offset for de-embedding [m].
  double d_offset = 0.0;

  // Eigenvalue solver type for boundary mode calculation.
  EigenSolverBackend eigen_solver = EigenSolverBackend::DEFAULT;

  // Input excitation for driven & transient solver:
  // - Wave/Lumped ports with same index are excited together.
  // - 1-based index if excited; 0 if not excited.
  int excitation = 0;

  // Flag for boundary damping term in driven and transient simulations.
  bool active = true;

  // List of boundary attributes for this wave port.
  std::vector<int> attributes = {};

  // Maximum number of iterations in linear solver.
  int ksp_max_its = 45;

  // Tolerance for linear solver.
  double ksp_tol = 1e-8;

  // Tolerance for eigenvalue solver.
  double eig_tol = 1e-6;

  // Print level for linear and eigenvalue solvers.
  int verbose = 0;

  WavePortData() = default;
  WavePortData(const json &port);
};

struct SurfaceCurrentData
{
public:
  // For each surface current source index, each element contains a list of attributes
  // making up a single element of a potentially multielement current source.
  std::vector<internal::ElementData> elements = {};

  SurfaceCurrentData() = default;
  SurfaceCurrentData(const json &source);
};

struct SurfaceFluxData
{
public:
  // Surface flux type.
  SurfaceFlux type = SurfaceFlux::ELECTRIC;

  // Flag for whether or not to consider the boundary as an infinitely thin two-sided
  // boundary for postprocessing.
  bool two_sided = false;

  // Coordinates of a point away from which to compute the outward flux (for orienting the
  // surface normal) [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};

  // Flag which indicates whether or not the center point was specified.
  bool no_center = true;

  // List of boundary attributes for this surface flux postprocessing index.
  std::vector<int> attributes = {};

  SurfaceFluxData() = default;
  SurfaceFluxData(const json &flux);
};

struct InterfaceDielectricData
{
public:
  // Type of interface dielectric for computing electric field energy participation ratios.
  InterfaceDielectric type = InterfaceDielectric::DEFAULT;

  // Dielectric interface thickness [m].
  double t = 0.0;

  // Relative permittivity.
  double epsilon_r = 0.0;

  // Loss tangent.
  double tandelta = 0.0;

  // List of boundary attributes for this interface dielectric postprocessing index.
  std::vector<int> attributes = {};

  InterfaceDielectricData() = default;
  InterfaceDielectricData(const json &dielectric);
};

struct FarFieldPostData
{
public:
  // List of boundary attributes to use for the surface integral.
  std::vector<int> attributes = {};

  // List of (theta, phi) where the wave-zone fields should be evaluated.
  // Units are radians.
  std::vector<std::pair<double, double>> thetaphis = {};

  FarFieldPostData() = default;
  FarFieldPostData(const json &farfield);

  bool empty() const { return thetaphis.empty(); };
};

struct BoundaryPostData
{
public:
  // List of all postprocessing boundary attributes.
  std::vector<int> attributes = {};

  // Boundary postprocessing objects.
  std::map<int, SurfaceFluxData> flux = {};
  std::map<int, InterfaceDielectricData> dielectric = {};
  FarFieldPostData farfield = {};

  BoundaryPostData() = default;
  BoundaryPostData(const json &postpro);
};

struct BoundaryData
{
public:
  // List of all boundary attributes (excluding postprocessing).
  std::vector<int> attributes = {};

  // List of all boundary attributes affected by mesh cracking.
  std::unordered_set<int> cracked_attributes = {};

  // Boundary objects.
  PecBoundaryData pec = {};
  PmcBoundaryData pmc = {};
  WavePortPecBoundaryData auxpec = {};
  FarfieldBoundaryData farfield = {};
  std::vector<ConductivityData> conductivity = {};
  std::vector<ImpedanceData> impedance = {};
  std::map<int, LumpedPortData> lumpedport = {};
  std::map<int, TerminalData> terminal = {};
  std::map<int, WavePortData> waveport = {};
  std::map<int, SurfaceCurrentData> current = {};
  PeriodicBoundaryData periodic = {};
  BoundaryPostData postpro = {};

  BoundaryData() = default;
  BoundaryData(const json &boundaries);
};

// Solver Configuration.

struct DrivenSolverData
{
public:
  // Explicit frequency samples [GHz].
  std::vector<double> sample_f = {};

  // Indices of frequency samples to explicitly add to the PROM.
  std::vector<std::size_t> prom_indices;

  // Indices of frequency samples on which to save fields to disk.
  std::vector<std::size_t> save_indices;

  // Restart iteration for a partial sweep. 1-based indexing. So 1 <= restart <= nr_freq *
  // nr_excitations.
  std::size_t restart = 1;

  // Error tolerance for enabling adaptive frequency sweep.
  double adaptive_tol = 0.0;

  // Maximum number of frequency samples for adaptive frequency sweep.
  std::size_t adaptive_max_size = 20;

  // Memory required for adaptive sampling convergence.
  std::size_t adaptive_memory = 2;

  // Gram-Schmidt orthogonalization used in PROM construction, separate from linear solver
  // orthogonalization option. Default to CGS2 for higher quality.
  Orthogonalization adaptive_solver_gs_orthog_type = Orthogonalization::CGS2;

  // Return circuit matrices from port with port excitation vectors.
  bool adaptive_circuit_synthesis = false;

  // Domain orthogonalization type for circuit synthesis weight matrix.
  DomainOrthogonalizationWeight adaptive_circuit_synthesis_domain_orthog =
      DomainOrthogonalizationWeight::ENERGY;

  DrivenSolverData() = default;
  DrivenSolverData(const json &driven);
};

struct EigenSolverData
{
public:
  // Target for shift-and-invert spectral transformation [GHz].
  double target = 0.0;

  // Eigenvalue solver relative tolerance.
  double tol = 1.0e-6;

  // Maximum iterations for eigenvalue solver.
  int max_it = -1;

  // Eigenvalue solver subspace dimension or maximum dimension before restart.
  int max_size = -1;

  // Desired number of eigenmodes.
  int n = 1;

  // Number of modes to write to disk.
  int n_post = 0;

  // Use operator scaling in order to increase numerical robustness.
  bool scale = true;

  // Compute and set a starting vector for the eigenvalue solver.
  bool init_v0 = true;
  bool init_v0_const = false;

  // Orthogonalize basis vectors using a mass matrix inner product, instead of generating
  // using a standard ℓ² (Euclidean) norm.
  bool mass_orthog = false;

  // Eigenvalue solver type.
  EigenSolverBackend type = EigenSolverBackend::DEFAULT;

  // For SLEPc eigenvalue solver, use linearized formulation for quadratic eigenvalue
  // problems.
  bool pep_linear = true;

  // Nonlinear eigenvalue solver type.
  NonlinearEigenSolver nonlinear_type = NonlinearEigenSolver::HYBRID;

  // For nonlinear problems, refine the linearized solution with a nonlinear eigensolver.
  bool refine_nonlinear = true;

  // For nonlinear problems using the hybrid approach, relative tolerance of the linear
  // eigenvalue solver used to generate the initial guess.
  double linear_tol = 1e-3;

  // Upper end of the target range for nonlinear eigenvalue solver [GHz]. A value <0
  // will use the default (3 * target).
  double target_upper = -1;

  // Update frequency of the preconditioner in the quasi-Newton nonlinear eigenvalue solver.
  int preconditioner_lag = 10;

  // Relative tolerance below which the preconditioner is not updated, regardless of the
  // lag.
  double preconditioner_lag_tol = 1e-4;

  // Maximum number of failed attempts with a given initial guess in the quasi-Newton
  // nonlinear eigenvalue solver.
  int max_restart = 2;

  EigenSolverData() = default;
  EigenSolverData(const json &eigenmode);
};

struct ElectrostaticSolverData
{
public:
  // Number of fields to write to disk.
  int n_post = 0;

  ElectrostaticSolverData() = default;
  ElectrostaticSolverData(const json &electrostatic);
};

struct MagnetostaticSolverData
{
public:
  // Number of fields to write to disk.
  int n_post = 0;

  MagnetostaticSolverData() = default;
  MagnetostaticSolverData(const json &magnetostatic);
};

struct TransientSolverData
{
public:
  // Time integration scheme type.
  TimeSteppingScheme type = TimeSteppingScheme::DEFAULT;

  // Excitation type for port excitation.
  Excitation excitation = Excitation::SINUSOIDAL;

  // Excitation parameters: frequency [GHz] and pulse width [ns].
  double pulse_f = 0.0;
  double pulse_tau = 0.0;

  // Upper bound of time interval [ns].
  double max_t = 1.0;

  // Step size for time stepping [ns].
  double delta_t = 1.0e-2;

  // Step increment for saving fields to disk.
  int delta_post = 0;

  // RK scheme order for SUNDIALS ARKODE integrators.
  // Max order for SUNDIALS CVODE integrator.
  // Not used for generalized α and Runge-Kutta integrators.
  int order = 2;

  // Adaptive time-stepping tolerances for CVODE and ARKODE.
  double rel_tol = 1e-4;
  double abs_tol = 1e-9;

  TransientSolverData() = default;
  TransientSolverData(const json &transient);
};

struct LinearSolverData
{
public:
  // Solver type.
  LinearSolver type = LinearSolver::DEFAULT;

  // Krylov solver type.
  KrylovSolver krylov_solver = KrylovSolver::DEFAULT;

  // Iterative solver relative tolerance.
  double tol = 1.0e-6;

  // Maximum iterations for iterative solver.
  int max_it = 100;

  // Maximum Krylov space dimension for GMRES/FGMRES iterative solvers.
  int max_size = -1;

  // Reuse previous solution as initial guess for Krylov solvers.
  int initial_guess = -1;

  // Maximum number of levels for geometric multigrid (set to 1 to disable multigrid).
  int mg_max_levels = 100;

  // Type of coarsening for p-multigrid.
  MultigridCoarsening mg_coarsening = MultigridCoarsening::LOGARITHMIC;

  // Controls whether or not to include in the geometric multigrid hierarchy the mesh levels
  // from uniform refinement.
  bool mg_use_mesh = true;

  // Number of iterations for preconditioners which support it. For multigrid, this is the
  // number of V-cycles per Krylov solver iteration.
  int mg_cycle_it = -1;

  // Use auxiliary space smoothers on geometric multigrid levels.
  int mg_smooth_aux = -1;

  // Number of pre-/post-smoothing iterations at each geometric or algebraic multigrid
  // level.
  int mg_smooth_it = 1;

  // Order of polynomial smoothing for geometric multigrid.
  int mg_smooth_order = -1;

  // Safety factors for eigenvalue estimates associated with Chebyshev smoothing for
  // geometric multigrid.
  double mg_smooth_sf_max = 1.0;
  double mg_smooth_sf_min = 0.0;

  // Smooth based on 4th-kind Chebyshev polynomials for geometric multigrid, otherwise
  // use standard 1st-kind polynomials.
  bool mg_smooth_cheby_4th = true;

  // For frequency domain applications, precondition linear systems with a real-valued
  // approximation to the system matrix.
  bool pc_mat_real = false;

  // For frequency domain applications, precondition linear systems with a shifted matrix
  // (makes the preconditioner matrix SPD).
  int pc_mat_shifted = -1;

  // Matrix symmetry type for sparse direct solvers (computed from problem type and
  // boundary conditions).
  MatrixSymmetry pc_mat_sym = MatrixSymmetry::UNSYMMETRIC;

  // For frequency domain applications, use the complex-valued system matrix in the sparse
  // direct solver.
  bool complex_coarse_solve = false;

  // Drop small entries (< numerical ε) in the system matrix used in the sparse direct
  // solver.
  bool drop_small_entries = false;

  // Reuse the sparsity pattern (reordering) for repeated factorizations.
  bool reorder_reuse = true;

  // Choose left or right preconditioning.
  PreconditionerSide pc_side = PreconditionerSide::DEFAULT;

  // Specify details for the column ordering method in the symbolic factorization for sparse
  // direct solvers.
  SymbolicFactorization sym_factorization = SymbolicFactorization::DEFAULT;

  // Low-rank and butterfly compression parameters for sparse direct solvers which support
  // it (mainly STRUMPACK).
  SparseCompression strumpack_compression_type = SparseCompression::NONE;

  double strumpack_lr_tol = 1.0e-3;
  int strumpack_lossy_precision = 16;
  int strumpack_butterfly_l = 1;

  // Option to enable 3D process grid for SuperLU_DIST solver.
  bool superlu_3d = false;

  // Option to use vector or scalar Pi-space corrections for the AMS preconditioner.
  bool ams_vector_interp = false;

  // Option to tell the AMS solver that the operator is singular, like for magnetostatic
  // problems.
  int ams_singular_op = -1;

  // Option to use aggressive coarsening for Hypre AMG solves (with BoomerAMG or AMS).
  // Typically use this when the operator is positive definite.
  int amg_agg_coarsen = -1;

  // Maximum number of iterations of the AMS solver.
  int ams_max_it = -1;

  // Relative tolerance for solving linear systems in divergence-free projector.
  double divfree_tol = 1.0e-12;

  // Maximum number of iterations for solving linear systems in divergence-free projector.
  int divfree_max_it = 1000;

  // Relative tolerance for solving linear systems in the error estimator.
  double estimator_tol = 1.0e-6;

  // Maximum number of iterations for solving linear systems in the error estimator.
  int estimator_max_it = 10000;

  // Use geometric multigrid + AMG for error estimator linear solver preconditioner (instead
  // of just Jacobi).
  bool estimator_mg = false;

  // Enable different variants of Gram-Schmidt orthogonalization for GMRES/FGMRES iterative
  // solvers and SLEPc eigenvalue solver.
  Orthogonalization gs_orthog = Orthogonalization::MGS;

  LinearSolverData() = default;
  LinearSolverData(const json &linear);
};

struct SolverData
{
public:
  // Approximation order.
  int order = 1;

  // Order above which to use partial assembly instead of full assembly.
  int pa_order_threshold = 1;

  // Include the order of det(J) in the order of accuracy for quadrature rule selection.
  bool q_order_jac = false;

  // Additional quadrature order of accuracy (in addition to 2p or 2p + order(|J|)) for
  // quadrature rule selection.
  int q_order_extra = 0;

  // Device used to configure MFEM.
  Device device = Device::CPU;

  // Backend for libCEED (https://libceed.org/en/latest/gettingstarted/#backends).
  std::string ceed_backend = "";

  // Solver objects.
  DrivenSolverData driven = {};
  EigenSolverData eigenmode = {};
  ElectrostaticSolverData electrostatic = {};
  MagnetostaticSolverData magnetostatic = {};
  TransientSolverData transient = {};
  LinearSolverData linear = {};

  SolverData() = default;
  SolverData(const json &solver);
};

// Calculate the number of steps from [start, end) in increments of delta. Will only include
// end if it is a multiple of delta beyond start.
int GetNumSteps(double start, double end, double delta);

// Parse a string as a direction vector, returning the direction and coordinate system.
std::pair<std::array<double, 3>, CoordinateSystem>
ParseStringAsDirection(std::string str, bool required = true);

// Validation functions for cross-field checks. Return empty string if valid, otherwise
// return all error messages concatenated.
std::optional<std::string> Validate(const BoundaryData &boundaries);

}  // namespace palace::config

// Forward declare Units for Nondimensionalize functions.
namespace palace
{
class Units;
}  // namespace palace

namespace palace::config
{

// Nondimensionalization functions. Each function nondimensionalizes the fields of the
// given struct using the provided Units.
void Nondimensionalize(const Units &units, RefinementData &data);
void Nondimensionalize(const Units &units, MaterialData &data);
void Nondimensionalize(const Units &units, ProbeData &data);
void Nondimensionalize(const Units &units, CurrentDipoleData &data);
void Nondimensionalize(const Units &units, ConductivityData &data);
void Nondimensionalize(const Units &units, ImpedanceData &data);
void Nondimensionalize(const Units &units, LumpedPortData &data);
void Nondimensionalize(const Units &units, PeriodicBoundaryData &data);
void Nondimensionalize(const Units &units, WavePortData &data);
void Nondimensionalize(const Units &units, SurfaceFluxData &data);
void Nondimensionalize(const Units &units, InterfaceDielectricData &data);
void Nondimensionalize(const Units &units, EigenSolverData &data);
void Nondimensionalize(const Units &units, DrivenSolverData &data);
void Nondimensionalize(const Units &units, TransientSolverData &data);

}  // namespace palace::config

#endif  // PALACE_UTILS_CONFIG_FILE_HPP
