// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_MODE_EIGEN_SOLVER_HPP
#define PALACE_LINALG_MODE_EIGEN_SOLVER_HPP

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
// E_n = e_n_tilde / (i*kn) is performed by the caller (see GetPhysicalMode on the
// owning operator).
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
                  BoundaryModeOperator *bmo = nullptr);

  ~ModeEigenSolver() = default;

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

  // Propagation constant kn for mode i, recovered from the shift of the most recent Solve.
  std::complex<double> GetPropagationConstant(int i) const;

  // Load eigenvector i, back-transform the H1 block from the VD variable ẽn = i·kn·En to
  // physical En = ẽn/(i·kn), and power-normalize the combined vector to unit Poynting
  // power. Returns the propagation constant kn. et and en are aliased views into e0.
  std::complex<double> GetPhysicalMode(int i, double omega, ComplexVector &e0,
                                       ComplexVector &et, ComplexVector &en) const;

  // Poynting power P = (1/2) conj(kn)/omega · etᴴ·Btt·et + i/(2·omega) · etᴴ·Atn·En
  // for physical (et, En).
  std::complex<double> ComputePoyntingPower(double omega, std::complex<double> kn,
                                            const ComplexVector &et,
                                            const ComplexVector &en) const;

  // Heuristic classifier: |Im kn| < 0.1·|Re kn| and Re kn > 0.
  static bool IsPropagating(std::complex<double> kn);

  // Access the assembled Btt and Atn matrices (needed for power normalization).
  const mfem::HypreParMatrix *GetBtt() const { return Bttr.get(); }
  const mfem::HypreParMatrix *GetAtnr() const { return Atnr.get(); }
  const mfem::HypreParMatrix *GetAtni() const { return Atni.get(); }

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

  // Optional boundary-mode operator reference; when non-null, its FE space hierarchies
  // and per-level DBC tdof lists drive p-multigrid preconditioning. When null (3D wave
  // port submesh path), the sparse-direct preconditioner is used.
  BoundaryModeOperator *bmo = nullptr;

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

}  // namespace palace

#endif  // PALACE_LINALG_MODE_EIGEN_SOLVER_HPP
