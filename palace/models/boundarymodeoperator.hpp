// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
#define PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP

#include <complex>
#include <memory>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class FiniteElementSpace;
class MaterialOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;

//
// Configuration for the BoundaryModeOperator, parameterizing the differences between the
// 2D mode analysis (ModeAnalysisSolver) and 3D wave port (WavePortOperator) eigenvalue
// problems.
//
struct BoundaryModeOperatorConfig
{
  // Material property mappings. For mode analysis, these come from
  // GetAttributeToMaterial() (volume attrs). For wave ports, from
  // GetBdrAttributeToMaterial() (boundary attrs).
  const mfem::Array<int> *attr_to_material = nullptr;

  // Inverse permeability tensor. Used for Atn, Ant, and Btt bilinear forms.
  const mfem::DenseTensor *inv_permeability = nullptr;

  // Curl-curl inverse permeability. In 2D mode analysis this is the scalar mu_inv
  // (GetCurlCurlInvPermeability); for 3D wave ports it is the full tensor
  // (GetInvPermeability), with subsequent normal projection.
  const mfem::DenseTensor *curlcurl_inv_permeability = nullptr;

  // Real part of permittivity tensor. Used for Att, Ant bilinear forms.
  const mfem::DenseTensor *permittivity_real = nullptr;

  // Scalar permittivity for the Ann mass matrix. In 2D mode analysis this is
  // GetPermittivityScalar(); for 3D wave ports it is GetPermittivityReal() with
  // subsequent normal projection.
  const mfem::DenseTensor *permittivity_scalar = nullptr;

  // Imaginary part of permittivity (loss tangent contribution). Set to nullptr when
  // there is no loss tangent.
  const mfem::DenseTensor *permittivity_imag = nullptr;
  // Scalar (1x1) imaginary permittivity for H1 mass in 2D. Set to nullptr for 3D.
  const mfem::DenseTensor *permittivity_imag_scalar = nullptr;

  // Surface normal vector for 3D wave port boundaries. Set to nullptr for 2D domain
  // meshes (mode analysis), which do not need normal projection.
  const mfem::Vector *normal = nullptr;

  // Whether the material operator reports a nonzero loss tangent.
  bool has_loss_tangent = false;

  // Domain conductivity tensor. Used if non-null for damping contribution to Att.
  const mfem::DenseTensor *conductivity = nullptr;
  bool has_conductivity = false;

  // London penetration depth (1/lambda_L^2 tensor). Adds a mass term to the stiffness
  // operator: (1/lambda_L^2)(E, F). Used for superconductor modeling.
  const mfem::DenseTensor *inv_london_depth = nullptr;
  // Scalar (1x1) London depth for the H1 Ann mass in 2D. Set to nullptr for 3D (which
  // uses normal projection of inv_london_depth instead).
  const mfem::DenseTensor *inv_london_depth_scalar = nullptr;
  bool has_london_depth = false;

  // Reference to MaterialOperator for ceed attribute access. Required when boundary
  // operators are specified.
  const MaterialOperator *mat_op = nullptr;

  // Boundary operators for Robin BC contributions to Att. All optional (nullptr = none).
  // Non-const because the Add*Coefficients methods are not const.
  SurfaceImpedanceOperator *surf_z_op = nullptr;
  FarfieldBoundaryOperator *farfield_op = nullptr;
  SurfaceConductivityOperator *surf_sigma_op = nullptr;

  // Eigenvalue solver configuration.
  int num_modes = 1;
  int num_vec = -1;
  double eig_tol = 1.0e-6;
  EigenvalueSolver::WhichType which_eig = EigenvalueSolver::WhichType::LARGEST_REAL;

  // Linear solver configuration.
  const config::LinearSolverData *linear = nullptr;
  EigenSolverBackend eigen_backend = EigenSolverBackend::DEFAULT;

  // Verbosity level for solvers.
  int verbose = 0;
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
//   Btt = -(mu^{-1} u, v)
//   Btn = -(mu^{-1} u, grad v) = Atn^T (transpose coupling from Equation 2)
//
// The A matrix is upper block-triangular (Ant = 0). The B matrix has Btn coupling from
// Equation 2 and Bnn = 0. This is a standard GEP with shift-and-invert.
//
// The eigenvector contains [e_t; e_n_tilde] where e_n_tilde = i*kn*E_n. Recovery of
// E_n = e_n_tilde / (i*kn) is performed in the driver, not here.
//
class BoundaryModeOperator
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
  BoundaryModeOperator(const BoundaryModeOperatorConfig &config,
                       const FiniteElementSpace &nd_fespace,
                       const FiniteElementSpace &h1_fespace,
                       const mfem::Array<int> &dbc_tdof_list,
                       MPI_Comm solver_comm = MPI_COMM_NULL);

  ~BoundaryModeOperator();

  // Assemble frequency-dependent matrices and solve the eigenvalue problem. The shift
  // sigma = -kn_target^2 is applied. An optional initial space vector can be provided
  // for eigenvalue solver warm-starting.
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
  // Configuration (stored for Solve-time assembly).
  BoundaryModeOperatorConfig config;

  // References to FE spaces (not owned).
  const FiniteElementSpace &nd_fespace;
  const FiniteElementSpace &h1_fespace;

  // Essential boundary condition true DOF list for the combined block system.
  mfem::Array<int> dbc_tdof_list;

  // Cached FE space sizes.
  int nd_size, h1_size;

  // Frequency-independent assembled matrices.
  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  // Btn: transpose of Atn, coupling from Equation 2: -(mu^{-1} u, grad_t v).
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

  // Btt: ND mass -(mu^{-1} u, v).
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

  // Set up the linear solver (GMRES + sparse direct preconditioner).
  void SetUpLinearSolver(MPI_Comm comm);

  // Set up the eigenvalue solver (SLEPc or ARPACK).
  void SetUpEigenSolver(MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
