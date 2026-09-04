// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MODE_EIGEN_SOLVER_HPP
#define PALACE_MODELS_MODE_EIGEN_SOLVER_HPP

#include <complex>
#include <memory>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/modeoperatorassembly.hpp"
#include "models/waveportreducedmodel.hpp"
#include "utils/labels.hpp"

namespace palace
{

class BoundaryModeOperator;
class FarfieldBoundaryOperator;
class FiniteElementSpace;
class FiniteElementSpaceHierarchy;
class MaterialOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;
class SurfaceRationalImpedanceOperator;

template <typename OperType>
class BlockDiagonalPreconditioner;

namespace config
{
struct LinearSolverData;
}  // namespace config

//
// Linear eigenvalue solver for the 2D boundary-mode GEP. Uses the Vardapetyan–Demkowicz
// substitution e_n_tilde = i*kn*E_n (with e_t = E_t) to linearize the coupled transverse
// and normal curl-curl equations into
//
//   [Att  Atn] [et]           [Btt  Btn] [et]
//   [ 0   Ann] [en] = lambda  [ 0    0 ] [en]
//
// solved via shift-and-invert on sigma = -kn_target^2. Frequency-independent matrices
// (Atn, Btn = -Atn^T, Btt) come either from a BoundaryModeOperator (2D domain path) or
// local assembly (3D wave-port submesh path). Frequency-dependent Att/Ann blocks are
// evaluated from the shared mode_assembly component bank. The solver builds the block
// system, drives KSP + EPS, and returns eigenpairs; wave-port projection is delegated to
// WavePortReducedModel.

class ModeEigenSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  using ReducedModelStats = WavePortReducedModel::Stats;

  // Bare FE space constructor (WavePort). No multigrid — Att/Ann are evaluated from the
  // shared parametric component bank. `solver_comm` restricts solver setup to
  // port ranks; pass MPI_COMM_NULL to fall back to the FE space comm. For 3D wave port
  // submeshes `normal` is the outward surface normal.
  ModeEigenSolver(const MaterialOperator &mat_op, const mfem::Vector *normal,
                  SurfaceImpedanceOperator &surf_z_op,
                  FarfieldBoundaryOperator &farfield_op,
                  SurfaceConductivityOperator &surf_sigma_op,
                  SurfaceRationalImpedanceOperator &surf_rz_op,
                  const FiniteElementSpace &nd_fespace,
                  const FiniteElementSpace &h1_fespace,
                  const mfem::Array<int> &dbc_tdof_list, int num_modes, int num_vec,
                  double eig_tol, EigenvalueSolver::WhichType which_eig,
                  const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
                  int verbose, MPI_Comm solver_comm = MPI_COMM_NULL);

  // BoundaryModeOperator constructor (BoundaryMode 2D domain path). Aliases frequency-
  // independent matrices from the operator, delegates frequency-dependent assembly to
  // it, and uses its hierarchies for p-multigrid preconditioning when GetNumLevels() > 1.
  // The driver bakes any 3D-parent frame into iodata before BMO construction, so this
  // path is frame-agnostic (normal unused on the BMO side).
  ModeEigenSolver(BoundaryModeOperator &bmo, const mfem::Array<int> &dbc_tdof_list,
                  int num_modes, int num_vec, double eig_tol,
                  EigenvalueSolver::WhichType which_eig,
                  const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
                  int verbose);

  ~ModeEigenSolver() = default;

  // Solve the shifted GEP with shift sigma = -kn_target^2. Ranks configured without a
  // solver (WavePort non-port ranks with empty FE spaces) contribute to assembly only
  // and return num_converged = 0. omega may be complex (eigenmode nonlinear solve, where
  // the cross-section EVP is solved at the genuinely complex eigenvalue ω = -i·λ); the
  // shift sigma stays real (a pure algebraic centering). For real omega (imag = 0) the
  // assembly and solve reduce bit-for-bit to the real-frequency form.
  SolveResult Solve(std::complex<double> omega, double sigma,
                    const ComplexVector *initial_space = nullptr);

  // Collect exact real-frequency eigenvectors for a per-port reduced model. Reduced
  // evaluation remains disabled until EnableReducedModel is called, so adaptive offline
  // HDM samples always use the exact eigensolver while seeding the basis.
  void ConfigureReducedModelTraining(std::size_t max_basis_size);

  // Enable guarded Rayleigh-Ritz evaluation for subsequent real-frequency solves. The
  // reduced tolerance is derived conservatively from the adaptive driven tolerance and
  // the configured eigensolver tolerance. Complex-frequency solves always remain exact.
  void EnableReducedModel(double adaptive_tol);

  const ReducedModelStats &GetReducedModelStats() const
  {
    return reduced_model->GetStats();
  }
  std::size_t GetReducedBasisSize() const { return reduced_model->GetBasisSize(); }
  double GetReducedTolerance() const { return reduced_model->GetTolerance(); }

  std::complex<double> GetEigenvalue(int i) const;
  void GetEigenvector(int i, ComplexVector &x) const;
  double GetError(int i, EigenvalueSolver::ErrorType type) const;

  // Propagation constant kn for mode i, recovered from the shift of the most recent Solve.
  std::complex<double> GetPropagationConstant(int i) const;

  // Heuristic classifier: |Im kn| < 0.1*|Re kn| and Re kn > 0.
  static bool IsPropagating(std::complex<double> kn);

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
  SurfaceImpedanceOperator &surf_z_op;
  FarfieldBoundaryOperator &farfield_op;
  SurfaceConductivityOperator &surf_sigma_op;
  SurfaceRationalImpedanceOperator &surf_rz_op;

  // References to FE spaces (not owned). Finest level of bmo's hierarchies when set.
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Owning BoundaryModeOperator for the 2D domain path (null for WavePort). Used for
  // accessing FE hierarchies + per-level DBC lists during p-multigrid preconditioning.
  BoundaryModeOperator *bmo = nullptr;

  // Complete frequency-parametric components shared by exact and reduced solves.
  std::unique_ptr<mode_assembly::ModeOperatorModel> mode_op_model;

  // Essential boundary condition true DOF list for the combined block system.
  mfem::Array<int> dbc_tdof_list;

  // Cached FE space sizes.
  int nd_size, h1_size;

  // Parameters from the latest exact assembly (also used for one-time ROM validation).
  double sigma_cached = 0.0;
  std::complex<double> last_assembled_omega = 0.0;
  double last_assembled_sigma = 0.0;

  // Frequency-independent block matrices. On the BMO path these alias BMO-owned matrices;
  // on the bare ctor path they point at the owned_* members below.
  const mfem::HypreParMatrix *Atnr = nullptr;
  const mfem::HypreParMatrix *Atni = nullptr;
  const mfem::HypreParMatrix *Btnr = nullptr;
  const mfem::HypreParMatrix *Btni = nullptr;
  const mfem::HypreParMatrix *Bttr = nullptr;

  // Owning storage for the bare ctor (WavePort) path. BMO path leaves these null.
  std::unique_ptr<mfem::HypreParMatrix> owned_Atnr, owned_Atni;
  std::unique_ptr<mfem::HypreParMatrix> owned_Btnr, owned_Btni;
  std::unique_ptr<mfem::HypreParMatrix> owned_Bttr;

  std::unique_ptr<ComplexOperator> opB;

  // Frequency-dependent block A operator (rebuilt each Solve).
  std::unique_ptr<ComplexOperator> opA;

  // Eigenvalue and linear solvers (null on processes without solver_comm).
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Communicator used by the eigensolver. For wave ports this is the sub-communicator of
  // ranks owning port unknowns, not the parent FE-space communicator.
  MPI_Comm solver_comm = MPI_COMM_NULL;

  // Wave-port-only projection/training state is isolated from the exact eigensolver.
  std::unique_ptr<WavePortReducedModel> reduced_model;

  // Best portable single-vector initial space for subsequent exact solves.
  ComplexVector warm_start;

  // Permutation that maps external mode index to eigensolver index, sorted by ascending
  // Re{kn}. This ensures consistent mode ordering regardless of eigensolver backend.
  std::vector<int> mode_perm;

  // Evaluate frequency-dependent Att and Ann from shared components, then build block A
  // (MPI collective on the FE-space communicator). Omega may be complex.
  void AssembleFrequencyDependent(std::complex<double> omega, double sigma);

  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

  // Build the 2x2 block B matrix: [Btt, 0; Btn, 0].
  ComplexHypreParMatrix BuildSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                           const mfem::HypreParMatrix *Btti,
                                           const mfem::HypreParMatrix *Btnr,
                                           const mfem::HypreParMatrix *Btni,
                                           const mfem::HypreParMatrix *Dnn) const;

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

  // Shared ctor body: frequency-independent matrix assembly and linear + eigenvalue
  // solver setup. `solver_comm` is the WavePort sub-communicator when not null; it's
  // derived from the FE space comm for the bare-ctor default path and for the
  // hierarchy ctor.
  void Init(MPI_Comm solver_comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_MODE_EIGEN_SOLVER_HPP
