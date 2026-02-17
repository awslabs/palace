// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_ROM_OPERATOR_HPP
#define PALACE_MODELS_ROM_OPERATOR_HPP

#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/units.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

// Helper classs to define inner product used in the PROM when doing circuit synthesis. The
// weight matrix is the simple domain mass matrix of the space with the tdof of the port
// zeroed out (using the diag policy materials) and replaced by the boundary port matrix.
//
// This 'hybrid' approach is different from the conventional mass matrix $M_{ab} = (v,
// w)_{\Omega} + \sum_i (v, w)_{\Gamma_i}$ where the tdof of the port contributed both to
// the domain $\Omega$ overlap and the port overlap $\Gamma_i$. Our method basically the
// tdof on the ports fully with $\Gamma_i$. This ensures correct orthogonality of PROM basis
// with the boundary form used when calculating the system response matrix (S-Matrix).
//
// The ParSumOperator defined in linalg/rap.hpp can only have a single diag policy so we
// can't use that approach here.
struct HybridBulkBoundaryOperator
{
  // Need this since SetEssentialTrueDofs makes pointer to these values, even though it
  // looks like it is passing by references.
  mfem::Array<int> port_tdof_list;

  std::unique_ptr<ParOperator> W_inner_product_weight_bulk = {};
  std::unique_ptr<ParOperator> W_inner_product_weight_port = {};

  HybridBulkBoundaryOperator(mfem::Array<int> &&port_tdof_list_,
                             std::unique_ptr<ParOperator> &&W_inner_product_weight_bulk_,
                             std::unique_ptr<ParOperator> &&W_inner_product_weight_port_)
    : port_tdof_list{std::move(port_tdof_list_)},
      W_inner_product_weight_bulk{std::move(W_inner_product_weight_bulk_)},
      W_inner_product_weight_port{std::move(W_inner_product_weight_port_)}
  {
    // Don't need to zero of PEC, since vectors entering inner product should already have
    // the PEC contribution removed.
    W_inner_product_weight_bulk->SetEssentialTrueDofs(port_tdof_list, Operator::DIAG_ZERO);
    MFEM_VERIFY((W_inner_product_weight_bulk->NumCols() ==
                 W_inner_product_weight_port->NumCols()) &&
                    (W_inner_product_weight_bulk->NumRows() ==
                     W_inner_product_weight_port->NumRows()),
                "Mismatch sizes of inner product matrices!")
  }

  auto NumCols() const { return W_inner_product_weight_bulk->NumCols(); }

  auto NumRows() const { return W_inner_product_weight_bulk->NumRows(); }

  void Mult(const Vector &x, Vector &y) const
  {
    W_inner_product_weight_bulk->Mult(x, y);
    W_inner_product_weight_port->AddMult(x, y);
  }

  double InnerProduct(const Vector &x, const Vector &y, Vector &v_workspace) const
  {
    v_workspace.SetSize(x.Size());
    v_workspace.UseDevice(x.UseDevice());
    Mult(x, v_workspace);
    return linalg::LocalDot(y, v_workspace);
  }

  double InnerProduct(MPI_Comm comm, const Vector &x, const Vector &y,
                      Vector &v_workspace) const
  {
    auto dot = InnerProduct(x, y, v_workspace);
    Mpi::GlobalSum(1, &dot, comm);
    return dot;
  }
};

// Class for handling minimal-rational interpolation of solutions in frequency space. Used
// as an error indicator and for selecting the next frequency sample points in PROM
// construction. Each excitation gets a separate MRI, so sample frequencies are not shared.
class MinimalRationalInterpolation
{
private:
  // (Complex-valued) upper-triangular matrix R from orthogonalization of the HDM samples.
  // Minimal rational interpolant (MRI) defined by the vector q of interpolation weights and
  // support points z is used as an error indicator.
  std::vector<ComplexVector> Q;
  std::size_t dim_Q = 0;
  Eigen::MatrixXcd R;
  Eigen::VectorXcd q;
  std::vector<double> z;

public:
  MinimalRationalInterpolation(std::size_t max_size);
  void AddSolutionSample(double omega, const ComplexVector &u, MPI_Comm comm,
                         Orthogonalization orthog_type);
  std::vector<double> FindMaxError(std::size_t N) const;

  const auto &GetSamplePoints() const { return z; }
};

//
// A class handling projection-based reduced order model (PROM) construction and use for
// adaptive fast frequency sweeps.
//
class RomOperator
{
protected:
  // Reference to HDM discretization (not owned).
  // TODO(C++20): Use std::reference_wrapper with incomplete types.
  SpaceOperator &space_op;

  // Used for constructing & reuse of RHS1.
  int excitation_idx_cache = 0;

  // HDM system matrices and excitation:
  // - System matrix is: A(ω) = K + iω C - ω² M + A2(ω).
  // - Excitation / drive: = iω RHS1 + RHS2(ω).
  // - Vector r is internal vector workspace of size RHS
  // - The non-quadratic in ω operators A2(ω) and RHS2(ω) are built on fly in SolveHDM.
  // - Need to recompute RHS1 when excitation index changes (cf excitation_idx_cache).
  std::unique_ptr<ComplexOperator> K, M, C, A2;
  ComplexVector RHS1, RHS2, r;

  // System properties: will be set when calling SetExcitationIndex & SolveHDM.
  bool has_A2 = true;
  bool has_RHS1 = true;
  bool has_RHS2 = true;

  // HDM linear system solver and preconditioner.
  std::unique_ptr<ComplexKspSolver> ksp;

  // PROM matrices and vectors. Projected matrices are Mr = Vᴴ M V where V is the reduced
  // order basis defined below.
  Eigen::MatrixXcd Kr, Mr, Cr;  // Extend during UpdatePROM as modes are added
  Eigen::VectorXcd RHS1r;       // Need to recompute drive vector on excitation change.

  // Frequency dependant PROM matrix Ar and RHSr are assembled and used only during
  // SolvePROM. Define them here so memory allocation can be reused in "online" evaluation.
  Eigen::MatrixXcd Ar;
  Eigen::VectorXcd RHSr;

  // PROM reduced-order basis (real-valued).
  std::vector<Vector> V;
  Orthogonalization orthog_type;

  // Weight operator for PROM basis if doing synthesis, in order to have correct
  // orthogonality to port vectors and converge with mesh / order. Default empty, which
  // means the identity weight in the finite element basis.
  std::optional<HybridBulkBoundaryOperator> weight_op_W = {};

  // Label to distinguish port modes from solution projection and to print PROM matrices.
  std::vector<std::string> v_node_label;

  // Upper-triangular matrix R from orthogonalization procedure U = VR. Here U the HDM
  // fields added by `UpdatePROM`.
  Eigen::MatrixXd orth_R;

  // MRIs: one for each excitation index. Only used to pick new frequency sample point.
  std::map<int, MinimalRationalInterpolation> mri;

  // Helper function that normalizes PROM matrices so that they correspond to proper
  // admittance matrices on the ports. Also does unit conversion to physical (input) units.
  // Define so that 1.0 on port i corresponds to full (un-normalized solution), so you can
  // use Linv, Rinv, C directly.
  std::tuple<std::unique_ptr<Eigen::MatrixXcd>, std::unique_ptr<Eigen::MatrixXcd>,
             std::unique_ptr<Eigen::MatrixXcd>>
  CalculateNormalizedPROMMatrices(const Units &units) const;

public:
  RomOperator(const IoData &iodata, SpaceOperator &space_op,
              std::size_t max_size_per_excitation);

  // Return the HDM linear solver.
  const ComplexKspSolver &GetLinearSolver() const { return *ksp; }

  // Return PROM dimension.
  auto GetReducedDimension() const { return V.size(); }

  // Return set of sampled parameter points for basis construction.
  const auto &GetSamplePoints(int excitation_idx) const
  {
    return mri.at(excitation_idx).GetSamplePoints();
  }

  // Set excitation index to build corresponding RHS vector (linear in frequency part).
  void SetExcitationIndex(int excitation_idx);

  // Assemble and solve the HDM at the specified frequency.
  void SolveHDM(int excitation_idx, double omega, ComplexVector &u);

  // Add port fields to PROM for circuit synthesis to connect to the outside world. Requires
  // fields that are orthogonal to boundary overlap and primary fields, which requires
  // custom weight matrix in orthogonality — see weight_op_W.
  void AddLumpedPortModesForSynthesis();

  // Add field configuration to the reduced-order basis and update the PROM. Requires a name
  // "node_label". This will be printed in the header of the csv files when printing PROM
  // matrices. It is needed to distinguish port and solution field configuration as well as
  // to reconstruct if field configuration are pure real, imaginary or complex.
  void UpdatePROM(const ComplexVector &u, std::string_view node_label);

  // Add solution u to the minimal-rational interpolation for error estimation. MRI are
  // separated by excitation index.
  void UpdateMRI(int excitation_idx, double omega, const ComplexVector &u);

  // Assemble and solve the PROM at the specified frequency, expanding the solution back
  // into the high-dimensional space.
  void SolvePROM(int excitation_idx, double omega, ComplexVector &u);

  // Compute the location(s) of the maximum error in the range of the previously sampled
  // parameter points.
  std::vector<double> FindMaxError(int excitation_idx, std::size_t N = 1) const
  {
    return mri.at(excitation_idx).FindMaxError(N);
  }

  // Compute eigenvalue estimates for the current PROM system.
  std::vector<std::complex<double>> ComputeEigenvalueEstimates() const;

  // Print PROM matrices to file include in input (SI) units.
  void PrintPROMMatrices(const Units &units, const fs::path &post_dir) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_ROM_OPERATOR_HPP
