// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
#define PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP

#include <memory>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;
class HypreParMatrix;

}  // namespace mfem

namespace palace
{

//
// Geometric multigrid preconditioner using a given coarse solver for the provided
// hierarchy of finite element spaces. Optionally can be configured to use auxiliary space
// smoothing at each level, and one-sided V-cycles.
//
template <OperatorType OperType>
class GeometricMultigridSolver : public Solver<OperType>
{
  using VecType = typename Solver<OperType>::VecType;

private:
  // Number of V-cycles per preconditioner application.
  const int pc_it;

  // Prolongation operators (not owned).
  std::vector<const Operator *> P;

  // System matrices at each multigrid level (not owned).
  std::vector<const OperType *> A;
  std::vector<const mfem::Array<int> *> dbc_tdof_lists;

  // Smoothers for each level. Coarse-level solver is B[0].
  mutable std::vector<std::unique_ptr<Solver<OperType>>> B;

  // Temporary vectors for preconditioner application. The type of these is dictated by the
  // MFEM Operator interface for multiple RHS.
  mutable std::vector<VecType> X, Y, R;

  // Enable timer contribution for Timer::KSP_COARSE_SOLVE.
  bool use_timer;

  // Polynomial order of the pre- and post-smoothers of symmetric V-cycles.
  const int smooth_order;

  // One-sided V-cycles, which pre-smooth but do not post-smooth, at twice the smoothing
  // order (Phillips and Fischer, 2022): the preconditioner is not symmetric, so only for
  // GMRES and FGMRES.
  bool one_sided;

  // Optional subdomain correction at the finest level, Y += S A_S⁻¹ Sᵀ (X - A Y), with S
  // the selection of the subdomain true DOFs of the multigrid operator and A_S = Sᵀ A S the
  // corresponding principal submatrix (complex-valued operators only).
  std::unique_ptr<Solver<OperType>> subdomain_solver;
  std::unique_ptr<mfem::HypreParMatrix> subdomain_S;
  mutable VecType subdomain_r, subdomain_x;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

  // Set up and apply the subdomain correction.
  void SetUpSubdomainCorrection(const BaseMultigridOperator<OperType> &mg_op);
  void SubdomainCorrection(int l) const;

public:
  GeometricMultigridSolver(MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G, int cycle_it,
                           int smooth_it, int cheby_order, double cheby_sf_max,
                           double cheby_sf_min, bool cheby_4th_kind,
                           bool one_sided = false);

  ~GeometricMultigridSolver() override;

  // Set the solver for the subdomain correction, configured with the subdomain of the
  // multigrid operator at the next call to SetOperator.
  void SetSubdomainSolver(std::unique_ptr<Solver<OperType>> &&solver)
  {
    subdomain_solver = std::move(solver);
  }

  void SetOperator(const OperType &op) override;

  void Mult(const VecType &x, VecType &y) const override;

  void EnableTimer() { use_timer = true; }

  auto GetNumLevels() const { return A.size(); }

  // Switch between one-sided and symmetric V-cycles.
  void SetOneSided(bool one_sided_cycle);
  bool IsOneSided() const { return one_sided; }
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
