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

template <typename OperType>
class BlockDiagonalPreconditioner;

namespace config
{
struct LinearSolverData;
}  // namespace config

//
// Linear eigenvalue solver for 2D boundary mode analysis. Uses the Vardapetyan–Demkowicz
// substitution e_n_tilde = i*kn*E_n (with e_t = E_t) to linearize the coupled transverse
// and normal curl-curl equations into the shifted GEP
//
//   [Att  Atn] [et]           [Btt  Btn] [et]
//   [ 0   Ann] [en] = lambda  [ 0    0 ] [en]
//
// where Att is the frequency-dependent ND block (curl-curl + mass + BC-t), Atn is the
// ND/H1 gradient coupling, Ann is the H1 block (stiffness + mass + BC-n), Btt is the ND
// mass, and Btn = -Atn^T. Solved via shift-and-invert on sigma = -kn_target^2.
//
// Callers get raw eigenvectors [et; e_n_tilde] via GetEigenvector and propagation
// constants via GetPropagationConstant. ApplyVDBackTransform recovers physical en = e_n_tilde
// / (i*kn); ComputePoyntingPower provides the shared Poynting integral. Each driver
// handles its own phase / power normalization on top of these primitives.
//
class ModeEigenSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  // Bare FE space constructor, used by WavePort (no p-multigrid). Matrix assembly runs
  // on the FE space communicator; the solver runs on `solver_comm` when non-null (port
  // ranks only) or on the FE space communicator otherwise. For 3D wave port submeshes
  // `normal` is the outward surface normal; boundary operators may each be null.
  ModeEigenSolver(const MaterialOperator &mat_op, const mfem::Vector *normal,
                  SurfaceImpedanceOperator *surf_z_op,
                  FarfieldBoundaryOperator *farfield_op,
                  SurfaceConductivityOperator *surf_sigma_op,
                  const FiniteElementSpace &nd_fespace,
                  const FiniteElementSpace &h1_fespace,
                  const mfem::Array<int> &dbc_tdof_list, int num_modes, int num_vec,
                  double eig_tol, EigenvalueSolver::WhichType which_eig,
                  const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
                  int verbose, MPI_Comm solver_comm = MPI_COMM_NULL);

  // Hierarchy constructor, used by BoundaryMode. p-multigrid preconditioning is enabled
  // automatically when the hierarchies have > 1 level; otherwise sparse-direct. For a 2D
  // domain mesh pass normal = nullptr.
  ModeEigenSolver(const MaterialOperator &mat_op, const mfem::Vector *normal,
                  SurfaceImpedanceOperator *surf_z_op,
                  FarfieldBoundaryOperator *farfield_op,
                  SurfaceConductivityOperator *surf_sigma_op,
                  FiniteElementSpaceHierarchy &nd_fespaces,
                  FiniteElementSpaceHierarchy &h1_fespaces,
                  FiniteElementSpaceHierarchy &h1_aux_fespaces,
                  const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
                  const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists,
                  const std::vector<mfem::Array<int>> &h1_aux_dbc_tdof_lists,
                  const mfem::Array<int> &dbc_tdof_list, int num_modes, int num_vec,
                  double eig_tol, EigenvalueSolver::WhichType which_eig,
                  const config::LinearSolverData &linear, EigenSolverBackend eigen_backend,
                  int verbose);

  ~ModeEigenSolver() = default;

  // Solve the shifted GEP with shift sigma = -kn_target^2. Ranks configured without a
  // solver (WavePort non-port ranks with empty FE spaces) contribute to assembly only
  // and return num_converged = 0.
  SolveResult Solve(double omega, double sigma,
                    const ComplexVector *initial_space = nullptr);

  // Access converged eigenvalues and eigenvectors.
  std::complex<double> GetEigenvalue(int i) const;
  void GetEigenvector(int i, ComplexVector &x) const;

  // Get the eigenpair error for the i-th converged eigenvalue.
  double GetError(int i, EigenvalueSolver::ErrorType type) const;

  // Propagation constant kn for mode i, recovered from the shift of the most recent Solve.
  std::complex<double> GetPropagationConstant(int i) const;

  // Alias the ND and H1 halves of a pre-loaded eigenvector e0 = [e_t_tilde; e_n_tilde] as
  // et / en, and apply the Vardapetyan–Demkowicz back-transform en := ẽn / (i·kn) so en
  // holds the physical En. Used by both drivers after GetEigenvector (and any desired
  // phase normalization) to produce the physical mode fields; each driver applies its
  // own power normalization flavor afterward.
  void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, ComplexVector &et,
                            ComplexVector &en) const;

  // Poynting power P = (1/2) conj(kn)/omega · etᴴ Btt et + i/(2·omega) · etᴴ Atn En for
  // physical (et, En). Used by BoundaryMode for power normalization and impedance
  // postprocessing; WavePort uses a different normalization and does not call this.
  std::complex<double> ComputePoyntingPower(double omega, std::complex<double> kn,
                                            const ComplexVector &et,
                                            const ComplexVector &en) const;

  // Heuristic classifier: |Im kn| < 0.1·|Re kn| and Re kn > 0.
  static bool IsPropagating(std::complex<double> kn);

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

  // Shift sigma = -kn_target^2 from the last Solve (used by GetPropagationConstant).
  double sigma_cached = 0.0;

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

  // p-multigrid context. Non-null only on the hierarchy ctor path (BoundaryMode); the
  // bare ctor (WavePort) leaves them null and the solver uses sparse-direct. The finest
  // level of nd_fespaces / h1_fespaces is what nd_fespace / h1_fespace above reference.
  FiniteElementSpaceHierarchy *nd_fespaces = nullptr;
  FiniteElementSpaceHierarchy *h1_fespaces = nullptr;
  FiniteElementSpaceHierarchy *h1_aux_fespaces = nullptr;
  const std::vector<mfem::Array<int>> *nd_dbc_tdof_lists = nullptr;
  const std::vector<mfem::Array<int>> *h1_dbc_tdof_lists = nullptr;
  const std::vector<mfem::Array<int>> *h1_aux_dbc_tdof_lists = nullptr;

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
