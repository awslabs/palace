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
// Linear eigenvalue solver for the 2D boundary-mode GEP. Uses the Vardapetyan–Demkowicz
// substitution e_n_tilde = i*kn*E_n (with e_t = E_t) to linearize the coupled transverse
// and normal curl-curl equations into
//
//   [Att  Atn] [et]           [Btt  Btn] [et]
//   [ 0   Ann] [en] = lambda  [ 0    0 ] [en]
//
// solved via shift-and-invert on sigma = -kn_target^2. Frequency-independent matrices
// (Atn, Btn = -Atn^T, Btt) come either from a BoundaryModeOperator (2D domain path) or
// from local assembly through mode_assembly:: free functions (3D wave port submesh path,
// which has no BMO). The solver itself just builds the block system, drives KSP + EPS,
// and returns eigenpairs; physical-mode reconstruction (VD back-transform, Poynting
// power) lives on BoundaryModeOperator for the 2D path and inline in WavePortData for
// the 3D path.

// Frequency-independent and frequency-dependent assembly of the 2D boundary-mode GEP
// blocks, shared by BoundaryModeOperator (2D domain path) and ModeEigenSolver's bare
// ctor (3D wave port submesh path).
namespace mode_assembly
{

using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                         std::unique_ptr<mfem::HypreParMatrix>>;

// Atn = -(mu^{-1} grad_t u, v). ND / H1 gradient coupling, real-only.
ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op);

// Btt = (mu^{-1} u, v). ND mass, real-only (positive).
ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op);

// Att = mu_cc^{-1} curl-curl  -  omega^2 eps mass  -  sigma (mu^{-1} mass) + BC-t
// (impedance / absorbing / conductivity). Frequency- and shift-dependent. `omega` may
// be complex (analytic continuation of the cross-section modal problem onto the complex
// frequency plane); for real omega the result is bit-for-bit identical to the legacy
// real-omega assembly. The shift `sigma` stays real (pure algebraic spectral shift).
ComplexHypreParMatrix AssembleAtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal,
                                  SurfaceImpedanceOperator &surf_z_op,
                                  FarfieldBoundaryOperator &farfield_op,
                                  SurfaceConductivityOperator &surf_sigma_op,
                                  std::complex<double> omega, double sigma);

// Ann = -(mu^{-1} grad u, grad v) + omega^2 (eps u, v) + BC-n. Frequency-dependent.
// farfield_op and surf_sigma_op contribute impedance / loss terms on the H1 block.
// `omega` may be complex (see AssembleAtt); real omega reproduces the legacy assembly.
ComplexHypreParMatrix
AssembleAnn(const FiniteElementSpace &h1_fespace, const MaterialOperator &mat_op,
            const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
            FarfieldBoundaryOperator &farfield_op,
            SurfaceConductivityOperator &surf_sigma_op, std::complex<double> omega);

// Alias the ND and H1 halves of a pre-loaded eigenvector e0 = [e_t_tilde; e_n_tilde] as
// et / en, and apply the Vardapetyan–Demkowicz back-transform en := ẽn / (i·kn) so en
// holds the physical En. Pure scalar op, no MPI.
void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en);

}  // namespace mode_assembly

class ModeEigenSolver
{
public:
  struct SolveResult
  {
    int num_converged;
    double sigma;
  };

  // Bare FE space constructor (WavePort). No multigrid — Att/Ann are re-assembled at
  // each Solve via mode_assembly free functions. `solver_comm` restricts solver setup to
  // port ranks; pass MPI_COMM_NULL to fall back to the FE space comm. For 3D wave port
  // submeshes `normal` is the outward surface normal.
  ModeEigenSolver(const MaterialOperator &mat_op, const mfem::Vector *normal,
                  SurfaceImpedanceOperator &surf_z_op,
                  FarfieldBoundaryOperator &farfield_op,
                  SurfaceConductivityOperator &surf_sigma_op,
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
  // and return num_converged = 0. `omega` may be complex to compute the cross-section
  // mode at a complex frequency (the BC stamping i*kn(omega)*M then carries the exact
  // complex kn); for real omega the assembly and result are unchanged. `sigma` is the
  // real algebraic spectral shift and is independent of Im(omega).
  SolveResult Solve(std::complex<double> omega, double sigma,
                    const ComplexVector *initial_space = nullptr);

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
  const mfem::Vector *normal;
  SurfaceImpedanceOperator &surf_z_op;
  FarfieldBoundaryOperator &farfield_op;
  SurfaceConductivityOperator &surf_sigma_op;

  // References to FE spaces (not owned). Finest level of bmo's hierarchies when set.
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Owning BoundaryModeOperator for the 2D domain path (null for WavePort). Used for
  // frequency-dependent Att/Ann assembly and for accessing FE hierarchies + per-level
  // DBC lists during p-multigrid preconditioner assembly.
  BoundaryModeOperator *bmo = nullptr;

  // Essential boundary condition true DOF list for the combined block system.
  mfem::Array<int> dbc_tdof_list;

  // Cached FE space sizes.
  int nd_size, h1_size;

  // Shift sigma = -kn_target^2 from the last Solve (used by GetPropagationConstant).
  double sigma_cached = 0.0;

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

  // Permutation that maps external mode index to eigensolver index, sorted by ascending
  // Re{kn}. This ensures consistent mode ordering regardless of eigensolver backend.
  std::vector<int> mode_perm;

  // Assemble frequency-dependent Att and Ann, then build block A (MPI collective on FE
  // space comm). Uses bmo->AssembleAtt/Ann on the BMO path; falls through to
  // mode_assembly:: free functions on the bare path.
  void AssembleFrequencyDependent(std::complex<double> omega, double sigma);

  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

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
