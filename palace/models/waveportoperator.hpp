// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
#define PALACE_MODELS_WAVE_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <unordered_map>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
class SumVectorCoefficient;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;
class Units;

namespace config
{

struct BoundaryData;
struct SolverData;
struct WavePortData;

}  // namespace config

enum class ProblemType : char;

//
// Helper class for wave port boundaries in a model.
//
template <typename OperType>
class BlockDiagonalPreconditioner;
class FiniteElementSpace;
class FiniteElementSpaceHierarchy;

// (ModeEigenSolverConfig removed — solver parameters passed directly to constructor.)

//
// Configuration for p-multigrid preconditioning of the boundary mode eigenvalue problem.
// When provided, the ModeEigenSolver uses a block-diagonal GMG preconditioner
// (ND p-multigrid + H1 p-multigrid) instead of a sparse direct solver.
//
struct ModeEigenSolverMultigridConfig
{
  // ND and H1 FE space hierarchies for the diagonal blocks (not owned).
  FiniteElementSpaceHierarchy *nd_fespaces = nullptr;
  FiniteElementSpaceHierarchy *h1_fespaces = nullptr;

  // H1 auxiliary space hierarchy for Hiptmair distributive relaxation in the ND block.
  FiniteElementSpaceHierarchy *h1_aux_fespaces = nullptr;

  // Per-level essential BC true DOF lists for ND, H1, and H1 auxiliary blocks (not owned).
  std::vector<mfem::Array<int>> *nd_dbc_tdof_lists = nullptr;
  std::vector<mfem::Array<int>> *h1_dbc_tdof_lists = nullptr;
  std::vector<mfem::Array<int>> *h1_aux_dbc_tdof_lists = nullptr;
};

//
// Linear eigenvalue solver for 2D boundary mode computation using a direct linearization
// of the transverse curl-curl (Equation 1) and normal curl-curl (Equation 2) equations
// with the Vardapetyan-Demkowicz substitution e_n_tilde = i*kn*E_n, e_t = E_t.
//
// Assembles and solves the generalized eigenvalue problem:
//
//   [Att  Atn] [et]           [Btt  Btn] [et]
//   [ 0   Ann] [en] = lambda  [ 0    0 ] [en]
//
// where:
//   Att = mu_cc^{-1}(curl_t u, curl_t v) - omega^2(eps u, v) + BC-t
//   Atn = -(mu^{-1} grad_t u, v)
//   Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n (H1 stiffness + mass +
//          impedance boundary)
//   Btt = (mu^{-1} u, v)
//   Btn = -Atn^T (negative transpose coupling from Equation 2)
//
// The A matrix is upper block-triangular (Ant = 0). The B matrix has Btn coupling from
// Equation 2 and Bnn = 0. This is a standard GEP with shift-and-invert.
//
// The eigenvector contains [e_t; e_n_tilde] where e_n_tilde = i*kn*E_n. Recovery of
// E_n = e_n_tilde / (i*kn) is performed in the driver, not here.
//
class ModeEigenSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  // The constructor assembles frequency-independent matrices (Atn, Btn = Atn^T, Btt)
  // and configures linear and eigenvalue solvers. Frequency-dependent matrices (Att, Ann)
  // are assembled at solve time via AssembleFrequencyDependent(). Matrix assembly uses
  // the FE space communicator (all processes). If solver_comm != MPI_COMM_NULL, the
  // linear and eigenvalue solvers are configured on that communicator.
  // Constructor. Material properties and attribute mapping come from mat_op directly.
  // For 3D wave port submeshes, normal is the outward surface normal; for 2D domain
  // meshes, normal is nullptr. Boundary operators are optional (nullptr = none).
  ModeEigenSolver(const MaterialOperator &mat_op, const mfem::Vector *normal,
                  SurfaceImpedanceOperator *surf_z_op,
                  FarfieldBoundaryOperator *farfield_op,
                  SurfaceConductivityOperator *surf_sigma_op,
                  const FiniteElementSpace &nd_fespace,
                  const FiniteElementSpace &h1_fespace,
                  const mfem::Array<int> &dbc_tdof_list, int num_modes, int num_vec,
                  double eig_tol, EigenvalueSolver::WhichType which_eig,
                  const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
                  int verbose, MPI_Comm solver_comm = MPI_COMM_NULL,
                  const ModeEigenSolverMultigridConfig *mg_config = nullptr);

  ~ModeEigenSolver();

  // Assemble frequency-dependent matrices and solve the eigenvalue problem. The shift
  // sigma = -kn_target^2 is applied. An optional initial space vector can be provided
  // for eigenvalue solver warm-starting. When has_solver is true (default), the calling
  // process participates in the eigenvalue solve; when false (wave port non-port process),
  // only assembly is performed.
  SolveResult Solve(double omega, double sigma, bool has_solver = true,
                    const ComplexVector *initial_space = nullptr);

  // Access converged eigenvalues and eigenvectors.
  std::complex<double> GetEigenvalue(int i) const;
  void GetEigenvector(int i, ComplexVector &x) const;

  // Get the eigenpair error for the i-th converged eigenvalue.
  double GetError(int i, EigenvalueSolver::ErrorType type) const;

  // Access the assembled Btt matrix (needed for impedance postprocessing).
  const mfem::HypreParMatrix *GetBtt() const { return Bttr.get(); }

  // Get the true vector sizes for the ND and H1 FE spaces.
  int GetNDTrueVSize() const { return nd_size; }
  int GetH1TrueVSize() const { return h1_size; }

  // Access the linear solver (for metadata reporting). Returns nullptr if this process
  // does not have a solver configured (non-port process in wave port mode).
  const ComplexKspSolver *GetLinearSolver() const { return ksp.get(); }

private:
  // Solver parameters.
  int num_modes, num_vec;
  double eig_tol;
  EigenvalueSolver::WhichType which_eig;
  const config::LinearSolverData &linear;
  EigenSolverBackend eigen_backend;
  int verbose;

  // Material operator and boundary operators (not owned).
  const MaterialOperator &mat_op;
  const mfem::Vector *normal;
  SurfaceImpedanceOperator *surf_z_op;
  FarfieldBoundaryOperator *farfield_op;
  SurfaceConductivityOperator *surf_sigma_op;

  // References to FE spaces (not owned).
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Essential boundary condition true DOF list for the combined block system.
  mfem::Array<int> dbc_tdof_list;

  // Cached FE space sizes.
  int nd_size, h1_size;

  // Frequency-independent assembled matrices.
  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  // Btn: negative transpose of Atn: (mu^{-1} u, grad_t v).
  // Bttr: positive copy of Btt for external use (postprocessing).
  std::unique_ptr<mfem::HypreParMatrix> Atnr, Atni;
  std::unique_ptr<mfem::HypreParMatrix> Btnr, Btni;
  std::unique_ptr<mfem::HypreParMatrix> Bttr;
  std::unique_ptr<ComplexOperator> opB;

  // Frequency-dependent block A operator (rebuilt each Solve).
  std::unique_ptr<ComplexOperator> opA;

  // Eigenvalue and linear solvers (null on processes without solver_comm).
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Permutation that maps external mode index to eigensolver index, sorted by ascending
  // Re{kn}. This ensures consistent mode ordering regardless of eigensolver backend.
  std::vector<int> mode_perm;

  // Assemble frequency-dependent Att and Ann, then build block A (MPI collective on FE
  // space comm).
  void AssembleFrequencyDependent(double omega, double sigma);

  // Private helper methods for bilinear form assembly.
  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

  // Att: curl-curl + mass + BC-t.
  ComplexHypreParMatrix AssembleAtt(double omega, double sigma) const;

  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  ComplexHypreParMatrix AssembleAtn() const;

  // Ann: H1 stiffness + mass + BC-n.
  // Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n impedance.
  ComplexHypreParMatrix AssembleAnn(double omega) const;

  // Btt: ND mass (mu^{-1} u, v) (positive).
  ComplexHypreParMatrix AssembleBtt() const;

  // Build the 2x2 block A matrix. The (1,0) block is -sigma * Btn from the
  // shift-and-invert transformation (nullptr when sigma = 0).
  ComplexHypreParMatrix
  BuildSystemMatrixA(const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
                     const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
                     const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
                     const mfem::HypreParMatrix *shifted_Btnr = nullptr) const;

  // Build the 2x2 block B matrix: [Btt, 0; Btn, 0].
  ComplexHypreParMatrix BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                           const mfem::HypreParMatrix *Btti,
                                           const mfem::HypreParMatrix *Btnr,
                                           const mfem::HypreParMatrix *Btni,
                                           const mfem::HypreParMatrix *Dnn) const;

  // Optional multigrid configuration (not owned, may be nullptr for sparse direct path).
  const ModeEigenSolverMultigridConfig *mg_config;

  // Non-owning pointer to the block preconditioner (for setting operators in Solve).
  BlockDiagonalPreconditioner<ComplexOperator> *block_pc_ptr = nullptr;

  // Multigrid preconditioner operators (owned, must outlive the GMG solver application).
  std::unique_ptr<ComplexMultigridOperator> att_mg_op, ann_mg_op;

  // Shifted off-diagonal operator -sigma*Btn for block lower-triangular preconditioning.
  std::unique_ptr<ComplexOperator> shifted_Btn_op;

  // Assemble preconditioner operators at all multigrid levels for the Att (ND) block.
  // Returns a ComplexMultigridOperator with primary (ND) and auxiliary (H1) operators.
  std::unique_ptr<ComplexMultigridOperator> AssembleAttPreconditioner(double omega,
                                                                      double sigma) const;

  // Assemble preconditioner operators at all multigrid levels for the Ann (H1) block.
  std::unique_ptr<ComplexMultigridOperator> AssembleAnnPreconditioner(double omega) const;

  // Set up the linear solver (GMRES + sparse direct preconditioner).
  void SetUpLinearSolver(MPI_Comm comm);

  // Set up the linear solver with p-multigrid block-diagonal preconditioning.
  void SetUpMultigridLinearSolver(MPI_Comm comm);

  // Set up the eigenvalue solver (SLEPc or ARPACK).
  void SetUpEigenSolver(MPI_Comm comm);
};

//
// Wave port data for a single port boundary.
//
class WavePortData
{
public:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Wave port properties.
  int mode_idx;
  double d_offset;
  int excitation;
  bool active;
  std::complex<double> kn0;
  double omega0;
  mfem::Vector port_normal;

private:
  // List of all boundary attributes making up this port boundary.
  mfem::Array<int> attr_list;

  // SubMesh data structures to define finite element spaces and grid functions on the
  // SubMesh corresponding to this port boundary.
  std::unique_ptr<Mesh> port_mesh;
  std::unique_ptr<mfem::FiniteElementCollection> port_nd_fec, port_h1_fec;
  std::unique_ptr<FiniteElementSpace> port_nd_fespace, port_h1_fespace;
  std::unique_ptr<mfem::ParTransferMap> port_nd_transfer, port_h1_transfer;
  std::unordered_map<int, int> submesh_parent_elems;
  mfem::Array<int> port_dbc_tdof_list;
  double mu_eps_max;

  // Submesh-specific material operator and boundary condition operators. Constructed on
  // the remapped submesh with rebuilt CEED data for correct attribute mapping.
  std::unique_ptr<MaterialOperator> port_mat_op;
  std::unique_ptr<SurfaceImpedanceOperator> port_surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> port_farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> port_surf_sigma_op;

  // Boundary mode eigenvalue problem solver.
  std::unique_ptr<ModeEigenSolver> mode_solver;
  ComplexVector v0, e0;

  // Communicator for processes which have elements for this port.
  MPI_Comm port_comm = MPI_COMM_NULL;
  int port_root = 0;

  // Grid functions storing the last computed electric field mode on the port, and stored
  // objects for computing functions of the port modes for use as an excitation or in
  // postprocessing.
  std::unique_ptr<GridFunction> port_E0t, port_E0n, port_S0t, port_E;
  std::unique_ptr<mfem::LinearForm> port_sr, port_si;

  // Voltage path for line integral (optional, for impedance postprocessing).
  std::vector<mfem::Vector> voltage_path;
  int voltage_n_samples = 100;
  bool has_voltage_coords = false;

  // Whether this is a 2D direct port (no submesh).
  bool is_2d_direct = false;

public:
  // 3D submesh constructor (existing): extracts submesh from parent mesh.
  WavePortData(const config::WavePortData &data, const IoData &iodata,
               const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
               mfem::ParFiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_attr);

  // 2D direct constructor: uses provided mesh/spaces directly. No submesh extraction,
  // no parent DOF transfer. Used by BoundaryModeOperator.
  WavePortData(const MaterialOperator &mat_op, const mfem::Vector *normal,
               SurfaceImpedanceOperator *surf_z_op, FarfieldBoundaryOperator *farfield_op,
               SurfaceConductivityOperator *surf_sigma_op, FiniteElementSpace &nd_fespace,
               FiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_tdof_list,
               int num_modes, int num_vec, double eig_tol,
               EigenvalueSolver::WhichType which_eig,
               const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
               int verbose, MPI_Comm solver_comm = MPI_COMM_NULL,
               const ModeEigenSolverMultigridConfig *mg_config = nullptr);

  ~WavePortData();

  [[nodiscard]] constexpr bool HasExcitation() const { return excitation != 0; }
  [[nodiscard]] bool HasVoltageCoords() const { return has_voltage_coords; }

  const auto &GetAttrList() const { return attr_list; }

  void Initialize(double omega);

  // Access the eigenvalue solver directly (for 2D direct mode used by
  // BoundaryModeOperator).
  ModeEigenSolver &GetEigenSolver() { return *mode_solver; }
  const ModeEigenSolver &GetEigenSolver() const { return *mode_solver; }

  HYPRE_BigInt GlobalTrueNDSize() const { return port_nd_fespace->GlobalTrueVSize(); }
  HYPRE_BigInt GlobalTrueH1Size() const { return port_h1_fespace->GlobalTrueVSize(); }

  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientReal() const;
  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientImag() const;

  std::unique_ptr<mfem::VectorCoefficient>
  GetModeFieldCoefficientReal(double scaling = 1.0) const;
  std::unique_ptr<mfem::VectorCoefficient>
  GetModeFieldCoefficientImag(double scaling = 1.0) const;

  // Characteristic impedance Z = |V|^2 / (2P) from the port mode voltage and unit power.
  // Requires voltage coordinates to be configured.
  std::complex<double> GetCharacteristicImpedance() const;

  double GetExcitationPower() const;

  // Excitation voltage from the normalized port mode field.
  std::complex<double> GetExcitationVoltage() const;

  std::complex<double> GetPower(GridFunction &E, GridFunction &B) const;
  std::complex<double> GetSParameter(GridFunction &E) const;

  // Voltage line integral V = integral of E . dl on the port, using the 3D E field.
  // Requires voltage coordinates to be configured.
  std::complex<double> GetVoltage(GridFunction &E) const;
};

//
// A class handling wave port boundaries and their postprocessing.
//
class WavePortOperator
{
private:
  // Mapping from port index to data structure containing port information.
  std::map<int, WavePortData> ports;

  // Flag which forces no printing during WavePortData::Print().
  bool suppress_output;
  double fc, kc;

  void SetUpBoundaryProperties(const IoData &iodata, const MaterialOperator &mat_op,
                               mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh);

  // Compute boundary modes for all wave port boundaries at the specified frequency.
  void Initialize(double omega);

public:
  WavePortOperator(const config::BoundaryData &boundaries, const config::SolverData &solver,
                   ProblemType problem_type, const Units &units,
                   const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace);
  WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                   mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace);

  // Access data structures for the wave port with the given index.
  const WavePortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Enable or suppress all outputs (log printing and fields to disk).
  void SetSuppressOutput(bool suppress) { suppress_output = suppress; }

  // Returns array of wave port attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix from wave ports.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries.
  void AddExcitationBdrCoefficients(int excitation_idx, double omega,
                                    SumVectorCoefficient &fbr, SumVectorCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
