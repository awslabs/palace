// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_CONFIGFILE_HPP
#define PALACE_CONFIGFILE_HPP

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json_fwd.hpp>

namespace palace::config
{

using json = nlohmann::json;

//
// Data structures for storing configuration file data.
//

namespace internal
{

template <typename DataType>
struct DataVector
{
protected:
  std::vector<DataType> vecdata = {};

public:
  [[nodiscard]] auto size() const { return vecdata.size(); }
  [[nodiscard]] auto empty() const { return vecdata.empty(); }
  [[nodiscard]] auto begin() const { return vecdata.begin(); }
  [[nodiscard]] auto end() const { return vecdata.end(); }
  [[nodiscard]] auto begin() { return vecdata.begin(); }
  [[nodiscard]] auto end() { return vecdata.end(); }
  [[nodiscard]] auto front() const { return vecdata.front(); }
  [[nodiscard]] auto back() const { return vecdata.back(); }
  [[nodiscard]] auto front() { return vecdata.front(); }
  [[nodiscard]] auto back() { return vecdata.back(); }
};

template <typename DataType>
struct DataMap
{
protected:
  // Map keys are the object indices for postprocessing.
  std::map<int, DataType> mapdata = {};

public:
  [[nodiscard]] auto &operator[](int i) { return mapdata[i]; }
  [[nodiscard]] auto &at(int i) { return mapdata.at(i); }
  [[nodiscard]] auto size() const { return mapdata.size(); }
  [[nodiscard]] auto empty() const { return mapdata.empty(); }
  [[nodiscard]] auto begin() const { return mapdata.begin(); }
  [[nodiscard]] auto end() const { return mapdata.end(); }
  [[nodiscard]] auto begin() { return mapdata.begin(); }
  [[nodiscard]] auto end() { return mapdata.end(); }
};

}  // namespace internal

struct ProblemData
{
public:
  // Simulation type.
  enum class Type
  {
    DRIVEN,
    EIGENMODE,
    ELECTROSTATIC,
    MAGNETOSTATIC,
    TRANSIENT,
    INVALID = -1
  };
  Type type = Type::INVALID;

  // Level of printing.
  int verbose = 1;

  // Output path for storing results.
  std::string output = "";

  void SetUp(json &config);
};

struct BoxRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Region bounding box limits [m].
  std::vector<double> bbmin = {}, bbmax = {};
};

struct SphereRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Sphere radius [m].
  double r = 0.0;

  // Sphere center [m].
  std::vector<double> center = {};
};

// Stores data specifying the adaptive mesh refinement algorithm.
struct AdaptiveRefinementData
{
  // Non-dimensional tolerance used to specify convergence of the AMR.
  double tolerance = 1e-3;
  // Maximum number of iterations to perform during the AMR.
  int max_its = 0;
  // Minimum number of iterations to perform during the AMR.
  int min_its = 0;
  // Dorfler update fraction. The set of marked elements is the minimum set
  // that contains update_fraction of the total error.
  double update_fraction = 0.25;
  // Whether or not to perform coarsening during the AMR. Incompatible with
  // `construct_geometric_multigrid`.
  double coarsening_fraction = 0.0;
  // Whether to construct a geometric multigrid using the generated sequence of
  // meshes.
  bool construct_geometric_multigrid = true;
  // Maximum difference in non-conformal refinements between two adjacent
  // elements. Default = 0 implies there is no constraint on local non-conformity.
  int max_nc_levels = 0;
  // If a refinement results in a greater number of DOFs than this value, no
  // future refinement will be allowed unless coarsening is allowed to occur.
  int dof_limit = 0;
  // Used in transient AMR. AMR is triggered when the error estimate rises
  // above tolerance, and then proceeds until the error estimate is below
  // on_update_tolerance_ratio * tolerance. This ensures that a mesh
  // generalizes into the future.
  double on_update_tolerance_ratio = 0.5;
  // Frequency with which to store the post processing results for a given
  // adaptation, e.g. save_step = 3 means save every third adaptation.
  int save_step = 0;
};

struct RefinementData
{
public:
  // Parallel uniform mesh refinement levels.
  int uniform_ref_levels = 0;
  // Whether to perform adaptive mesh refinement;
  bool use_adaptive_mesh_refinement = false;
  // Structure containing adaptivity specification data.
  // Only used if `use_adaptive_mesh_refinement==true`.
  AdaptiveRefinementData adaptation;

private:
  // Refinement data for mesh regions.
  std::vector<BoxRefinementData> boxlist = {};
  std::vector<SphereRefinementData> spherelist = {};

public:
  auto &GetBox(int i) { return boxlist[i]; }
  const auto &GetBoxes() const { return boxlist; }
  auto &GetBoxes() { return boxlist; }

  auto &GetSphere(int i) { return spherelist[i]; }
  const auto &GetSpheres() const { return spherelist; }
  auto &GetSpheres() { return spherelist; }

  void SetUp(json &model);
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

  // Partitioning file (if specified, does not compute a new partitioning).
  std::string partition = "";

  // Call MFEM's ReorientTetMesh as a check of mesh orientation after partitioning.
  bool reorient_tet = false;

  // Object controlling mesh refinement.
  RefinementData refinement = {};

  void SetUp(json &config);
};

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
};

struct MaterialDomainData : public internal::DataVector<MaterialData>
{
public:
  void SetUp(json &domains);
};

struct DomainDielectricData
{
public:
  // List of domain attributes for this domain dielectric postprocessing index.
  std::vector<int> attributes = {};
};

struct DomainDielectricPostData : public internal::DataMap<DomainDielectricData>
{
public:
  void SetUp(json &postpro);
};

struct ProbeData
{
public:
  // Physical space coordinates for the probe location [m].
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct ProbePostData : public internal::DataMap<ProbeData>
{
public:
  void SetUp(json &postpro);
};

struct DomainPostData
{
public:
  // Set of all postprocessing domain attributes.
  std::set<int> attributes;

  // Domain postprocessing objects.
  DomainDielectricPostData dielectric;
  ProbePostData probe;

  void SetUp(json &domains);
};

struct DomainData
{
public:
  // Set of all domain attributes.
  std::set<int> attributes = {};

  // Domain objects.
  MaterialDomainData materials = {};
  DomainPostData postpro = {};

  void SetUp(json &config);
};

struct PecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct PmcBoundaryData
{
public:
  // List of boundary attributes with PMC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct WavePortPecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions for wave ports.
  std::vector<int> attributes;

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct FarfieldBoundaryData
{
public:
  // Approximation order for farfield ABC.
  int order = 1;

  // List of boundary attributes with farfield absortbing boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
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
};

struct ConductivityBoundaryData : public internal::DataVector<ConductivityData>
{
public:
  void SetUp(json &boundaries);
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
};

struct ImpedanceBoundaryData : public internal::DataVector<ImpedanceData>
{
public:
  void SetUp(json &boundaries);
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

  // Flag for source term in driven and transient simulations.
  bool excitation = false;

  // For each lumped port index, each Node contains a list of attributes making up a single
  // element of a potentially multielement port.
  struct Node
  {
    // String defining source excitation field direction. Options are "X", "Y", "Z", or "R",
    // preceeded with a "+" or "-" for the direction. The first three options are for
    // uniform lumped ports and the last is for a coaxial lumped port (radial excitation).
    std::string direction = "";

    // List of boundary attributes for this lumped port element.
    std::vector<int> attributes = {};
  };
  std::vector<Node> nodes = {};
};

struct LumpedPortBoundaryData : public internal::DataMap<LumpedPortData>
{
public:
  void SetUp(json &boundaries);
};

struct WavePortData
{
public:
  // Flag for source term in driven and transient simulations.
  bool excitation = false;

  // Mode index for the numeric wave port.
  int mode_idx = 1;

  // Port offset for de-embedding [m].
  double d_offset = 0.0;

  // List of boundary attributes for this wave port.
  std::vector<int> attributes;
};

struct WavePortBoundaryData : public internal::DataMap<WavePortData>
{
public:
  void SetUp(json &boundaries);
};

struct SurfaceCurrentData
{
public:
  // For each surface current source index, each Node contains a list of attributes making
  // up a single element of a potentially multielement current source.
  struct Node
  {
    // String defining surface current source excitation direction. Options are "X", "Y",
    // "Z", or "R", preceeded with a "+" or "-" for the direction. The first three options
    // are for uniform sources and the last is for a coaxial source (radial excitation).
    std::string direction = "";

    // List of boundary attributes for this surface current element.
    std::vector<int> attributes = {};
  };
  std::vector<Node> nodes = {};
};

struct SurfaceCurrentBoundaryData : public internal::DataMap<SurfaceCurrentData>
{
public:
  void SetUp(json &boundaries);
};

struct CapacitanceData
{
public:
  // List of boundary attributes for this capacitance postprocessing index.
  std::vector<int> attributes;
};

struct CapacitancePostData : public internal::DataMap<CapacitanceData>
{
public:
  void SetUp(json &postpro);
};

struct InductanceData
{
public:
  // String defining global direction with which to orient the computed flux (influences the
  // computed sign). Options are "X", "Y", or "Z", preceeded with a "+" or "-" for the
  // direction.
  std::string direction = "";

  // List of boundary attributes for this inductance postprocessing index.
  std::vector<int> attributes = {};
};

struct InductancePostData : public internal::DataMap<InductanceData>
{
public:
  void SetUp(json &postpro);
};

struct InterfaceDielectricData
{
public:
  // Dielectric interface thickness [m].
  double ts = 0.0;

  // Loss tangent.
  double tandelta = 0.0;

  // Relative permittivity.
  double epsilon_r = 0.0;

  // Optional relative permittivity for metal-substrate, metal-air, and substrate-air
  // layers.
  double epsilon_r_ma = 0.0;
  double epsilon_r_ms = 0.0;
  double epsilon_r_sa = 0.0;

  // For each dielectric postprocessing index, each Node contains a list of attributes
  // sharing the same side value.
  struct Node
  {
    // String defining surface side for interior surfaces. Options are "X", "Y", or "Z",
    // preceeded with a "+" or "-" for the direction.
    std::string side = "";

    // List of domain or boundary attributes for this interface dielectric postprocessing
    // index.
    std::vector<int> attributes = {};
  };
  std::vector<Node> nodes = {};
};

struct InterfaceDielectricPostData : public internal::DataMap<InterfaceDielectricData>
{
public:
  void SetUp(json &postpro);
};

struct BoundaryPostData
{
public:
  // Set of all postprocessing boundary attributes.
  std::set<int> attributes;

  // Boundary postprocessing objects.
  CapacitancePostData capacitance = {};
  InductancePostData inductance = {};
  InterfaceDielectricPostData dielectric = {};

  void SetUp(json &boundaries);
};

struct BoundaryData
{
public:
  // Set of all boundary attributes.
  std::set<int> attributes = {};

  // Boundary objects.
  PecBoundaryData pec = {};
  PmcBoundaryData pmc = {};
  WavePortPecBoundaryData auxpec = {};
  FarfieldBoundaryData farfield = {};
  ConductivityBoundaryData conductivity = {};
  ImpedanceBoundaryData impedance = {};
  LumpedPortBoundaryData lumpedport = {};
  WavePortBoundaryData waveport = {};
  SurfaceCurrentBoundaryData current = {};
  BoundaryPostData postpro = {};

  void SetUp(json &config);
};

struct DrivenSolverData
{
public:
  // Lower bound of frequency sweep [GHz].
  double min_f = 0.0;

  // Upper bound of frequency sweep [GHz].
  double max_f = 0.0;

  // Step size for frequency sweep [GHz].
  double delta_f = 0.0;

  // Step increment for saving fields to disk.
  int delta_post = 0;

  // Only perform postprocessing on port boundaries, skipping domain interior.
  bool only_port_post = false;

  // Error tolerance for enabling adaptive frequency sweep.
  double adaptive_tol = 0.0;

  // Maximum number of frequency samples for adaptive frequency sweep.
  int adaptive_nmax = 0;

  // Number of candidate points for error metric calculation in adaptive
  // frequency sweep.
  int adaptive_ncand = 0;

  // Use error metric based on an a posteriori residual error estimate. Otherwise just use
  // the 2-norm of the HDM residual.
  bool adaptive_metric_aposteriori = false;

  // Restart iteration for a partial sweep.
  int rst = 1;

  void SetUp(json &solver);
};

struct EigenSolverData
{
public:
  // Target for shift-and-invert spectral transformation [GHz].
  double target = 0.0;

  // Eigenvalue solver relative tolerance.
  double tol = 1e-6;

  // Maximum iterations for eigenvalue solver.
  int max_it = -1;

  // Eigensolver subspace dimension or maximum dimension before restart.
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
  enum class Type
  {
    ARPACK,
    SLEPC,
    FEAST,
    DEFAULT,
    INVALID = -1
  };
  Type type = Type::DEFAULT;

  // For SLEPc eigenvalue solver, use linearized formulation for quadratic eigenvalue
  // problems.
  bool pep_linear = true;

  // Number of integration points used for the FEAST eigenvalue solver contour.
  int feast_contour_np = 4;

  // Parameters for the FEAST eigenvalue solver contour.
  double feast_contour_ub = 0.0;
  double feast_contour_ar = 1.0;

  // Use more than just the standard single moment for FEAST subspace construction.
  int feast_moments = 1;

  void SetUp(json &solver);
};

struct ElectrostaticSolverData
{
public:
  // Number of fields to write to disk.
  int n_post = 0;

  void SetUp(json &solver);
};

struct MagnetostaticSolverData
{
public:
  // Number of fields to write to disk.
  int n_post = 0;

  void SetUp(json &solver);
};

struct TransientSolverData
{
public:
  // Time integration scheme type.
  enum class Type
  {
    GEN_ALPHA,
    NEWMARK,
    CENTRAL_DIFF,
    DEFAULT,
    INVALID = -1
  };
  Type type = Type::DEFAULT;

  // Excitation type for port excitation.
  enum class ExcitationType
  {
    SINUSOIDAL,
    GAUSSIAN,
    DIFF_GAUSSIAN,
    MOD_GAUSSIAN,
    RAMP_STEP,
    SMOOTH_STEP,
    INVALID = -1
  };
  ExcitationType excitation = ExcitationType::INVALID;

  // Excitation parameters: frequency [GHz] and pulse width [ns].
  double pulse_f = 0.0;
  double pulse_tau = 0.0;

  // Upper bound of time interval [ns].
  double max_t = 1.0;

  // Step size for time stepping [ns].
  double delta_t = 1e-2;

  // Step increment for saving fields to disk.
  int delta_post = 0;

  // Only perform postprocessing on port boundaries, skipping domain interior.
  bool only_port_post = false;

  void SetUp(json &solver);
};

struct LinearSolverData
{
public:
  // Solver type.
  enum class Type
  {
    AMS,
    BOOMER_AMG,
    MUMPS,
    SUPERLU,
    STRUMPACK,
    STRUMPACK_MP,
    DEFAULT,
    INVALID = -1
  };
  Type type = Type::DEFAULT;

  // Krylov solver type.
  enum class KspType
  {
    CG,
    CGSYM,
    FCG,
    MINRES,
    GMRES,
    FGMRES,
    BCGS,
    BCGSL,
    FBCGS,
    QMRCGS,
    TFQMR,
    DEFAULT,
    INVALID = -1
  };
  KspType ksp_type = KspType::DEFAULT;

  // Iterative solver relative tolerance.
  double tol = 1e-6;

  // Maximum iterations for iterative solver.
  int max_it = 100;

  // Maximum Krylov space dimension for GMRES/FGMRES iterative solvers.
  int max_size = -1;

  // Enable modified Gram-Schmidt orthogonalization instead of classical for GMRES/FGMRES
  // Krylov solvers and SLEPc eigenvalue solver.
  bool orthog_mgs = false;
  bool orthog_cgs2 = false;

  // Reuse previous solution as initial guess for Krylov solvers.
  int ksp_initial_guess = -1;

  // Enable pipelined Krylov solver variants to reduce blocking communications.
  bool ksp_piped = false;

  // Enable hp-geometric multigrid coarsening, using the solver specified by the type member
  // at the coarsest level.
  bool mat_gmg = true;

  // Enable low-order refined (LOR) preconditioner construction. Only available for meshes
  // based on tensor elements.
  bool mat_lor = false;

  // For frequency domain applications, precondition linear systems with a shifted matrix
  // (makes the preconditoner matrix SPD).
  int mat_shifted = -1;

  // Number of iterations for preconditioners which support it. For multigrid, this is the
  // number of V-cycles per Krylov solver iteration.
  int mg_cycle_it = 1;

  // Number of pre-/post-smoothing iterations at each geometric or algebraic multigrid
  // level.
  int mg_smooth_it = 1;

  // Order of polynomial smoothing for geometric multigrid.
  int mg_smooth_order = 4;

  // Choose left or right preconditioning.
  enum class SideType
  {
    RIGHT,
    LEFT,
    DEFAULT,
    INVALID = -1
  };
  SideType pc_side_type = SideType::DEFAULT;

  // Choose left or right preconditioning.
  enum class SymFactType
  {
    METIS,
    PARMETIS,
    DEFAULT,
    INVALID = -1
  };
  SymFactType sym_fact_type = SymFactType::DEFAULT;

  // Low-rank and butterfly compression parameters for sparse direct solvers which support
  // it.
  enum class CompressionType
  {
    NONE,
    BLR,
    HSS,
    HODLR,
    ZFP,
    BLR_HODLR,
    ZFP_BLR_HODLR,
    INVALID = -1
  };
  CompressionType strumpack_compression_type = CompressionType::NONE;
  double strumpack_lr_tol = 1.0e-3;
  int strumpack_lossy_precision = 16;
  int strumpack_butterfly_l = 1;

  // Option to enable 3D process grid for SuperLU_DIST solver.
  bool superlu_3d = false;

  // Option to use vector or scalar Pi-space corrections for the AMS preconditioner.
  bool ams_vector = false;

  // Relative tolerance for solving linear systems in divergence-free projector.
  double divfree_tol = 1.0e-12;

  // Maximum number of iterations for solving linear systems in divergence-free projector.
  int divfree_max_it = 100;

  void SetUp(json &solver);
};

struct SolverData
{
public:
  // Approximation order.
  int order = 1;

  // Solver objects.
  DrivenSolverData driven = {};
  EigenSolverData eigenmode = {};
  ElectrostaticSolverData electrostatic = {};
  MagnetostaticSolverData magnetostatic = {};
  TransientSolverData transient = {};
  LinearSolverData linear = {};

  void SetUp(json &config);
};

}  // namespace palace::config

#endif  // PALACE_CONFIGFILE_HPP
