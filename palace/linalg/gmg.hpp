// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
#define PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/sliverpatch.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "utils/iodata.hpp"

namespace mfem
{

template <typename T>
class Array;
class ParMesh;

}  // namespace mfem

namespace palace
{

class FiniteElementSpace;

// Choice of the coarse solve of a complex operator by a multigrid solver (see
// GeometricMultigridSolver::CheckCoarseSolve), for later multigrid solvers of the same
// problem on the same processes: whether the complex-valued coarse solve was chosen, the
// global size of the coarse level, and on the node of the process the memory of the setup
// of the real-valued approximation and the memory in use at the check (bytes).
struct CoarseSolveChoice
{
  bool complex = false;
  long long coarse_size = 0;
  double real_memory = 0.0, check_memory = 0.0;
};

//
// Geometric multigrid preconditioner using a given coarse solver for the provided
// hierarchy of finite element spaces. Optionally can be configured to use auxiliary space
// smoothing at each level.
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

  // Enable timer contribution for Timer::KSP_COARSE_SOLVE, unless the coarse solver times
  // itself (setup vs. apply).
  bool use_timer, coarse_timer;

  // Whether the coarse solver has been set up for an operator.
  bool coarse_ready = false;

  // Exact solves on the patches of sliver tetrahedra at each level above the coarsest (see
  // SliverPatchSolver and SetSliverPatches), or nullptr, with the minimum dihedral angle
  // (degrees, 0 to disable) which defines slivers, the extent of their patches, the largest
  // patch size, and the sliver elements of the mesh they were found on.
  std::vector<std::unique_ptr<SliverPatchSolver<OperType>>> patches;
  std::vector<const FiniteElementSpace *> patch_spaces;
  // Levels whose patches were checked to cover at most max_patch_fraction of the true dofs
  // (when they have more than max_patch_dofs_small true dofs), or are disabled.
  std::vector<bool> patches_checked;
  static constexpr double max_patch_fraction = 0.01;
  static constexpr long long max_patch_dofs_small = 10000;

  // With automatic patch layers (negative sliver_layers), a level keeps no patch solves
  // when the dense work of its vertex-star patches, max over the processes of Σₚ nₚ² /
  // (true dofs of the process × order of the Chebyshev smoother), exceeds
  // max_patch_work_ratio; the layers (-1 for none) and work ratio of each level.
  static constexpr double max_patch_work_ratio = 12.0;
  std::vector<int> patch_layers;
  std::vector<double> patch_work;
  double sliver_angle = 0.0;
  int sliver_layers = 1, sliver_max_size = 2000;
  const mfem::ParMesh *sliver_mesh = nullptr;
  std::vector<int> slivers;
  long long sliver_count = 0;

  // Performance tally ids per level: pre-/post-smoothing, residual, restriction,
  // prolongation, sliver patch solves (index 0 holds the coarse solve).
  struct LevelTallies
  {
    int pre = -1, post = -1, res = -1, restr = -1, prol = -1, patch = -1;
  };
  std::vector<LevelTallies> perf;

  // Set up the sliver patch solvers for the level operators.
  void SetSliverPatchOperators();

  // The Chebyshev smoothers may raise their order at setup, from the spectrum of their own
  // level operator (see ChebyshevSmoother). The raised orders are kept only if a short
  // GMRES probe with the whole cycle does not predict that they need more than 25% more
  // iterations than the base orders (see CompareWithGmresProbe): the coarse-grid correction
  // can be inaccurate (for example the real-valued coarse solve of a complex operator), and
  // then the higher orders can slow down the outer Krylov solver. The decision is reused
  // for the later operators of the hierarchy (for example those of the nonlinear eigenvalue
  // iterations) whose levels raise the same orders.
  std::vector<int> probe_orders;
  bool probe_keep = true, print_probe = true;
  void CheckSmootherOrders();
  std::vector<std::size_t> GetRaisedSmootherLevels(std::vector<int> &orders,
                                                   std::vector<int> &base_orders) const;
  void SetSmootherOrders(const std::vector<std::size_t> &levels,
                         const std::vector<int> &orders);

  // The coarse solver of a complex operator can factorize the real-valued approximation
  // Ar + Ai of the coarse operator, which is a poor approximation for some problems (for
  // example with Floquet periodic or absorbing boundaries): the multigrid cycle is then
  // far from a contraction. With a factory for the coarse solver of the complex-valued
  // matrix (see SetComplexCoarseSolver), the first setup with the approximation checks it
  // (CheckCoarseSolve):
  //  1. A short GMRES probe of the coarse operator preconditioned by the coarse solve, from
  //     a fixed random start vector, finds the approximation poor if the residual is above
  //     coarse_probe_max_res after coarse_probe_max_its iterations (it stops once below).
  //  2. If the complex-valued factorization fits in the available memory (see
  //     CoarseMemoryFits), it is set up next to the approximation, and short GMRES probes
  //     of the whole cycle compare the approximation with it (see CompareWithGmresProbe),
  //     at the current Chebyshev orders and, if they were raised, with the approximation at
  //     the base orders. The complex-valued coarse solve is kept if it is predicted (from
  //     the ratio of the logarithms of the residuals) to need coarse_complex_min_ratio
  //     times fewer iterations than the better of the two: a more accurate coarse solve
  //     does not help when the cycle is limited by its smoothing, and it costs memory and
  //     time. A failure of the complex-valued setup aborts the run, after a hint printed
  //     before it (see TryComplexCoarseSolve).
  // With the complex-valued coarse solve, the auxiliary space smoothers which follow the
  // order of the primary space use complex_coarse_aux_order (as with an explicit one). The
  // choice is kept for the later operators of the hierarchy, and can be given to later
  // multigrid solvers of the same problem (see SetCoarseSolveChoice).
  static constexpr int coarse_probe_max_its = 10;
  static constexpr double coarse_probe_max_res = 1.0e-4, coarse_complex_min_ratio = 1.5;
  // Memory model of the complex-valued coarse solver relative to the memory R used by the
  // setup of the approximation on a node (with its factorization): over its first and
  // later factorizations (for example those of a nonlinear eigenvalue solver), it uses at
  // most complex_coarse_memory_ratio R plus coarse_memory_floor above the memory before its
  // setup (the ratio is given per solver by SetComplexCoarseSolver), where R is at least
  // coarse_memory_per_nonzero per nonzero of the real part. The allocations after
  // the setup (eigenvalue solver and Krylov bases, error estimation, postprocessing) are
  // reserved as coarse_memory_reserve times the memory of the node after the setup. A
  // previous choice is reused for coarse levels at most coarse_choice_max_growth times
  // larger or smaller than the one it was made for.
  static constexpr double coarse_memory_floor = 256.0 * 1024 * 1024,
                          coarse_memory_per_nonzero = 100.0, coarse_memory_reserve = 0.6,
                          coarse_choice_max_growth = 1.5;
  std::function<std::unique_ptr<Solver<OperType>>()> complex_coarse_factory;
  int complex_coarse_aux_order = 0;
  double complex_coarse_memory_ratio = 0.0;
  bool coarse_checked = false, complex_coarse = false, tried_complex_coarse = false;
  std::optional<CoarseSolveChoice> coarse_choice, prev_coarse_choice;
  std::vector<bool> aux_order_set;
  void CheckCoarseSolve(long coarse_memory);
  int CoarseMemoryFits(MPI_Comm comm, double real_memory, const CoarseSolveChoice *checked,
                       bool same_size, double *estimate, double *available,
                       double *node_memory = nullptr) const;
  void TryComplexCoarseSolve(bool trial, double estimate, double available);
  void SetAuxiliarySmootherOrders(int order);

  // Solver role for the performance report.
  std::string role;

  // Inclusive tally of the coarse solve for nested solvers (-1 for none).
  int coarse_tally;

  // Internal function to perform a single V-cycle iteration.
  void VCycle(int l, bool initial_guess) const;

public:
  GeometricMultigridSolver(
      MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
      const std::vector<const Operator *> &P, const std::vector<const Operator *> *G,
      int cycle_it, int smooth_it, int cheby_order, double cheby_sf_max,
      double cheby_sf_min, bool cheby_4th_kind, int cheby_max_order = 0,
      int cheby_order_aux = -1,
      MultigridAuxiliaryResidual aux_residual = MultigridAuxiliaryResidual::EXACT);
  GeometricMultigridSolver(const IoData &iodata, MPI_Comm comm,
                           std::unique_ptr<Solver<OperType>> &&coarse_solver,
                           const std::vector<const Operator *> &P,
                           const std::vector<const Operator *> *G = nullptr)
    : GeometricMultigridSolver(
          comm, std::move(coarse_solver), P, G, iodata.solver.linear.mg_cycle_it,
          iodata.solver.linear.mg_smooth_it, iodata.solver.linear.mg_smooth_order,
          iodata.solver.linear.mg_smooth_sf_max, iodata.solver.linear.mg_smooth_sf_min,
          iodata.solver.linear.mg_smooth_cheby_4th,
          iodata.solver.linear.mg_smooth_order_max,
          iodata.solver.linear.mg_smooth_order_aux, iodata.solver.linear.mg_aux_residual)
  {
  }

  void SetOperator(const OperType &op) override { SetOperator(op, false); }

  // Set the operators of all levels. With reuse_coarse, the coarse solver keeps the
  // operator (and factorization) of a previous call, if any, as an approximate coarse
  // solve for the new coarse operator.
  void SetOperator(const OperType &op, bool reuse_coarse);

  void Mult(const VecType &x, VecType &y) const override;

  void SetPerfTallies(const std::string &prefix, Timer::Index parent) override;

  void EnableTimer();

  void SetCoarseTally(int id) { coarse_tally = id; }

  // Use the coarse solver made by factory, which factorizes the complex-valued coarse
  // matrix, if a check at setup finds the coarse solver given at construction, of the
  // real-valued approximation of a complex operator, to be a poor approximation which
  // slows down the multigrid cycle (see CheckCoarseSolve). The auxiliary space smoothers
  // which follow the order of the primary space then use aux_order, if positive. The memory
  // of the complex-valued coarse solver is estimated as memory_ratio times that of the
  // setup of the approximation (see CoarseMemoryFits). Call before SetOperator.
  void SetComplexCoarseSolver(std::function<std::unique_ptr<Solver<OperType>>()> factory,
                              int aux_order, double memory_ratio)
  {
    complex_coarse_factory = std::move(factory);
    complex_coarse_aux_order = aux_order;
    complex_coarse_memory_ratio = memory_ratio;
  }

  // Whether the coarse solver made by the factory is in use, and whether the check tried
  // it (it did not if the approximation passed the coarse probe, if the complex-valued
  // coarse solve was estimated not to fit in memory or the available memory is unknown, or
  // if a previous choice was reused).
  bool UsesComplexCoarseSolver() const { return complex_coarse; }
  bool TriedComplexCoarseSolver() const { return tried_complex_coarse; }

  // The choice of the coarse solve made by the check (none before it), and a choice made
  // for an earlier operator of the same problem (for example by the multigrid solver of
  // an earlier shift or mesh), used without the check at the next SetOperator if the
  // coarse level has a similar size. Call SetCoarseSolveChoice before SetOperator.
  const std::optional<CoarseSolveChoice> &GetCoarseSolveChoice() const
  {
    return coarse_choice;
  }
  void SetCoarseSolveChoice(const CoarseSolveChoice &choice)
  {
    prev_coarse_choice = choice;
  }

  // Change the residuals between the primary and auxiliary space smoothers of the
  // distributive relaxation smoothers (see DistRelaxationSmoother).
  void SetAuxiliaryResidual(MultigridAuxiliaryResidual residual);

  // Check the residual updates of the distributive relaxation smoothers on the auxiliary
  // space too (the default) or only globally (see DistRelaxationSmoother).
  void SetAuxiliaryGradientCheck(bool check);

  // Add exact solves on the patches of the tetrahedra with a minimum dihedral angle below
  // min_angle (degrees; 0 disables them) to every level above the coarsest (see
  // SliverPatchSolver): after the pre-smoothing, the correction δ = Σₚ Rₚᵀ Aₚₚ⁻¹ Rₚ r of
  // the level residual r is added to the iterate. For complex operators, the residual is
  // then updated to r - A δ before its restriction (the patch solves precede the
  // coarse-grid correction); for real operators, which may be preconditioned for CG, the
  // same residual is restricted, so that the V-cycle stays symmetric. The patches hold the
  // entities incident to the vertices of each sliver (layers 1) or its own entities
  // (layers 0); with negative layers, the vertex stars unless they are expensive compared
  // with the smoothing of the level (see max_patch_work_ratio). Call before SetOperator.
  void SetSliverPatches(double min_angle, int layers = -1, int max_size = 2000);
};

}  // namespace palace

#endif  // PALACE_LINALG_GEOMETRIC_MULTIGRID_HPP
