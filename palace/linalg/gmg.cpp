// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gmg.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <new>
#include <string>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/fespace.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/gmresprobe.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/enum_string.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

namespace palace
{

template <OperatorType OperType>
GeometricMultigridSolver<OperType>::GeometricMultigridSolver(
    MPI_Comm comm, std::unique_ptr<Solver<OperType>> &&coarse_solver,
    const std::vector<const Operator *> &P, const std::vector<const Operator *> *G,
    int cycle_it, int smooth_it, int cheby_order, double cheby_sf_max, double cheby_sf_min,
    bool cheby_4th_kind, int cheby_max_order, int cheby_order_aux,
    MultigridAuxiliaryResidual aux_residual)
  : Solver<OperType>(), pc_it(cycle_it), P(P.begin(), P.end()), A(P.size() + 1),
    dbc_tdof_lists(P.size()), B(P.size() + 1), X(P.size() + 1), Y(P.size() + 1),
    R(P.size() + 1), use_timer(false), coarse_timer(false), patches(P.size() + 1),
    perf(P.size() + 1), role(perf::CurrentRole()), coarse_tally(-1)
{
  // Configure levels of geometric coarsening. Multigrid vectors will be configured at first
  // call to Mult. The multigrid operator size is set based on the finest space dimension.
  const auto n_levels = P.size() + 1;
  MFEM_VERIFY(n_levels > 0,
              "Empty finite element space hierarchy during multigrid solver setup!");
  MFEM_VERIFY(!G || G->size() == n_levels,
              "Invalid input for distributive relaxation smoother auxiliary space transfer "
              "operators (mismatch in number of levels)!");

  // Use the supplied level 0 (coarse) solver.
  B[0] = std::move(coarse_solver);

  // Configure level smoothers. Use distributive relaxation smoothing if an auxiliary
  // finite element space was provided.
  for (std::size_t l = 1; l < n_levels; l++)
  {
    if (G)
    {
      const int cheby_smooth_it = 1;
      B[l] = std::make_unique<DistRelaxationSmoother<OperType>>(
          comm, *(*G)[l], smooth_it, cheby_smooth_it, cheby_order, cheby_sf_max,
          cheby_sf_min, cheby_4th_kind, cheby_max_order, cheby_order_aux, aux_residual);
    }
    else
    {
      const int cheby_smooth_it = smooth_it;
      if (cheby_4th_kind)
      {
        B[l] = std::make_unique<ChebyshevSmoother<OperType>>(
            comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_max_order);
      }
      else
      {
        B[l] = std::make_unique<ChebyshevSmoother1stKind<OperType>>(
            comm, cheby_smooth_it, cheby_order, cheby_sf_max, cheby_sf_min);
      }
    }
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::EnableTimer()
{
  use_timer = true;
  if (auto *wrapper = dynamic_cast<MfemWrapperSolver<OperType> *>(B[0].get()))
  {
    wrapper->EnableTimer(true, Timer::KSP_COARSE_SOLVE);
    coarse_timer = true;
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetOperator(const OperType &op, bool reuse_coarse)
{
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;

  const auto *mg_op = dynamic_cast<const BaseMultigridOperator<OperType> *>(&op);
  MFEM_VERIFY(mg_op, "GeometricMultigridSolver requires a MultigridOperator or "
                     "ComplexMultigridOperator argument provided to SetOperator!");

  const auto n_levels = A.size();
  MFEM_VERIFY(
      mg_op->GetNumLevels() == n_levels &&
          (!mg_op->HasAuxiliaryOperators() || mg_op->GetNumAuxiliaryLevels() == n_levels),
      "Invalid number of levels for operators in multigrid solver setup!");
  long coarse_memory = -1;  // Memory of the coarse solver setup, if done here
  bool reuse_complex = false;
  double reuse_estimate = 0.0, reuse_available = 0.0;
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // A choice made for an earlier operator of the problem with a coarse level of similar
    // size is used without the check: the complex-valued coarse solver is then set up
    // directly, without the approximation, if its estimated memory fits.
    if (complex_coarse_factory && !coarse_checked && prev_coarse_choice && n_levels > 1)
    {
      const auto &Ac =
          static_cast<const ComplexParOperator &>(mg_op->GetOperatorAtLevel(0));
      long long n = Ac.Height();
      Mpi::GlobalSum(1, &n, Ac.GetComm());
      const double growth = double(n) / std::max(prev_coarse_choice->coarse_size, 1LL);
      if (growth <= coarse_choice_max_growth && growth >= 1.0 / coarse_choice_max_growth)
      {
        coarse_checked = true;
        coarse_choice = prev_coarse_choice;
        const int fits = prev_coarse_choice->complex
                             ? CoarseMemoryFits(Ac.GetComm(),
                                                prev_coarse_choice->real_memory *
                                                    std::pow(growth, 4.0 / 3.0),
                                                &*prev_coarse_choice,
                                                n == prev_coarse_choice->coarse_size,
                                                &reuse_estimate, &reuse_available)
                             : 0;
        reuse_complex = (fits > 0);
        Mpi::Print(Ac.GetComm(), " Coarse solve chosen for an earlier operator: {}\n",
                   reuse_complex                  ? "complex-valued"
                   : !prev_coarse_choice->complex ? "real-valued approximation"
                   : (fits < 0)
                       ? "real-valued approximation (the available memory is unknown)"
                       : "real-valued approximation (the complex-valued one is estimated "
                         "not to fit in memory)");
        perf::SetInfo(perf::LinearSolverPath(role, "CoarseSolveProbe"),
                      {{"Reused", true}, {"Complex", reuse_complex}});
      }
    }
  }
  for (std::size_t l = 0; l < n_levels; l++)
  {
    A[l] = &mg_op->GetOperatorAtLevel(l);
    MFEM_VERIFY(
        A[l]->Width() == A[l]->Height() &&
            (n_levels == 1 ||
             (A[l]->Height() == ((l < n_levels - 1) ? P[l]->Width() : P[l - 1]->Height()))),
        "Invalid operator sizes for GeometricMultigridSolver!");

    const auto *PtAP_l = dynamic_cast<const ParOperType *>(&mg_op->GetOperatorAtLevel(l));
    MFEM_VERIFY(
        PtAP_l,
        "GeometricMultigridSolver requires ParOperator or ComplexParOperator operators!");
    if (l < n_levels - 1)
    {
      dbc_tdof_lists[l] = PtAP_l->GetEssentialTrueDofs();
    }

    auto *dist_smoother = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get());
    if (dist_smoother)
    {
      MFEM_VERIFY(mg_op->HasAuxiliaryOperators(),
                  "Distributive relaxation smoother relies on both primary space and "
                  "auxiliary space operators for multigrid smoothing!");
      dist_smoother->SetOperators(mg_op->GetOperatorAtLevel(l),
                                  mg_op->GetAuxiliaryOperatorAtLevel(l));
    }
    else if ((l > 0 || !reuse_coarse || !coarse_ready) && !(l == 0 && reuse_complex))
    {
      if (l == 0)
      {
        if (complex_coarse_factory && !coarse_checked)
        {
          // The memory of the setup is measured for the check: the free memory of the heap
          // would hide it.
          memory_reporting::TrimHeap();
        }
        coarse_memory = memory_reporting::GetCurrentMemory();
      }
      B[l]->SetOperator(mg_op->GetOperatorAtLevel(l));
      coarse_ready = coarse_ready || (l == 0);
    }

    X[l].SetSize(A[l]->Height());
    Y[l].SetSize(A[l]->Height());
    R[l].SetSize(A[l]->Height());
    X[l].UseDevice(true);
    Y[l].UseDevice(true);
    R[l].UseDevice(true);
  }

  SetSliverPatchOperators();
  if (reuse_complex)
  {
    B[0].reset();
    TryComplexCoarseSolve(false, reuse_estimate, reuse_available);
    SetAuxiliarySmootherOrders(complex_coarse_aux_order);
    complex_coarse = true;
    coarse_ready = true;
  }
  if (coarse_memory >= 0)
  {
    CheckCoarseSolve(coarse_memory);
  }
  CheckSmootherOrders();

  // Resolved Chebyshev smoother parameters per level (level 0 is the coarse solver).
  auto Cheby = [](const Solver<OperType> *B) -> nlohmann::json
  {
    if (const auto *c = dynamic_cast<const ChebyshevSmoother<OperType> *>(B))
    {
      return {{"Order", c->GetOrder()},
              {"LambdaMax", c->GetMaxEigenvalue()},
              {"SinglePrecisionDirections", c->UsesSinglePrecision()}};
    }
    if (const auto *c = dynamic_cast<const ChebyshevSmoother1stKind<OperType> *>(B))
    {
      return {{"Order", c->GetOrder()},
              {"LambdaMax", c->GetMaxEigenvalue()},
              {"SinglePrecisionDirections", c->UsesSinglePrecision()}};
    }
    return nullptr;
  };
  auto levels = nlohmann::json::array();
  for (std::size_t l = 1; l < n_levels; l++)
  {
    if (const auto *d = dynamic_cast<const DistRelaxationSmoother<OperType> *>(B[l].get()))
    {
      nlohmann::json level = {{"Primary", Cheby(d->GetSmoother())},
                              {"Auxiliary", Cheby(d->GetAuxiliarySmoother())},
                              {"Residual", ToString(d->GetAuxiliaryResidual())}};
      if (d->GetGradientError() >= 0.0)
      {
        level["GradientError"] = d->GetGradientError();
      }
      if (d->GetAuxiliaryGradientError() >= 0.0)
      {
        level["AuxiliaryGradientError"] = d->GetAuxiliaryGradientError();
      }
      if (d->GetExactAuxiliaryInput())
      {
        level["ExactAuxiliaryInput"] = true;
      }
      levels.push_back(std::move(level));
    }
    else
    {
      levels.push_back(Cheby(B[l].get()));
    }
  }
  if (!levels.empty())
  {
    perf::SetInfo(perf::LinearSolverPath(role, "Chebyshev"), levels);
  }

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::CheckSmootherOrders()
{
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;
  const auto &op = *A.back();
  MPI_Comm comm = static_cast<const ParOperType &>(op).GetComm();

  std::vector<int> orders, base_orders;
  const auto raised = GetRaisedSmootherLevels(orders, base_orders);
  if (raised.empty())
  {
    return;
  }
  auto SetOrders = [&](bool base)
  { SetSmootherOrders(raised, base ? base_orders : orders); };

  if (orders == probe_orders)
  {
    SetOrders(!probe_keep);
    return;
  }

  // The probe of the raised (candidate) against the base (reference) orders, from a random
  // start vector seeded by process (see CompareWithGmresProbe).
  const int n = op.Height();
  VecType v0(n), t(n);
  v0.UseDevice(true);
  t.UseDevice(true);
  linalg::SetRandom(comm, v0);
  linalg::Normalize(comm, v0);
  const auto probe = CompareWithGmresProbe<VecType>(
      comm, v0, [&](bool base) { SetOrders(base); },
      [&](const VecType &x, VecType &y)
      {
        op.Mult(x, t);
        Mult(t, y);
      });
  const bool keep = !probe.use_reference;
  probe_orders = orders;
  probe_keep = keep;
  const int its = probe.its;
  const double res_raised = probe.res_candidate, res_base = probe.res_reference,
               ratio = probe.ratio;
  if (print_probe)
  {
    Mpi::Print(
        comm,
        " Raised Chebyshev smoothing order {}: GMRES probe residual {:.3e} after {:d} "
        "iterations, {:.3e} with the base order (predicted iteration ratio {:.2f})\n",
        keep ? "kept" : "not used", res_raised, its, res_base, ratio);
    print_probe = false;
  }
  perf::SetInfo(perf::LinearSolverPath(role, "ChebyshevOrderProbe"),
                {{"Iterations", its},
                 {"RaisedResidual", res_raised},
                 {"BaseResidual", res_base},
                 {"PredictedIterationRatio",
                  std::isfinite(ratio) ? nlohmann::json(ratio) : nlohmann::json()},
                 {"Raised", keep}});
}

template <OperatorType OperType>
std::vector<std::size_t> GeometricMultigridSolver<OperType>::GetRaisedSmootherLevels(
    std::vector<int> &orders, std::vector<int> &base_orders) const
{
  // The levels whose smoother order was raised at setup, with the orders of the Chebyshev
  // smoothers of the primary space per level (0 for none), through the distributive
  // relaxation smoother when there is one (whose auxiliary space smoother follows). The
  // orders are the maximum over the processes, so that all of them probe and use the same
  // orders.
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;
  MPI_Comm comm = static_cast<const ParOperType *>(A.back())->GetComm();
  std::vector<std::size_t> raised;
  orders.assign(A.size(), 0);
  base_orders.assign(A.size(), 0);
  for (std::size_t l = 1; l < A.size(); l++)
  {
    const Solver<OperType> *smoother = B[l].get();
    if (const auto *d = dynamic_cast<const DistRelaxationSmoother<OperType> *>(smoother))
    {
      smoother = d->GetSmoother();
    }
    if (const auto *c = dynamic_cast<const ChebyshevSmoother<OperType> *>(smoother))
    {
      orders[l] = c->GetOrder();
      base_orders[l] = c->GetBaseOrder();
    }
  }
  Mpi::GlobalMax(static_cast<int>(orders.size()), orders.data(), comm);
  for (std::size_t l = 1; l < A.size(); l++)
  {
    if (orders[l] > base_orders[l])
    {
      raised.push_back(l);
    }
  }
  return raised;
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetSmootherOrders(
    const std::vector<std::size_t> &levels, const std::vector<int> &orders)
{
  for (auto l : levels)
  {
    if (auto *d = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get()))
    {
      d->SetSmootherOrder(orders[l]);
    }
    else
    {
      static_cast<ChebyshevSmoother<OperType> &>(*B[l]).SetOrder(orders[l]);
    }
  }
}

template <OperatorType OperType>
int GeometricMultigridSolver<OperType>::CoarseMemoryFits(MPI_Comm comm, double real_memory,
                                                         const CoarseSolveChoice *checked,
                                                         bool same_size, double *estimate,
                                                         double *available,
                                                         double *node_memory) const
{
  // Whether the complex-valued coarse solver fits in the available memory of every node
  // (1), not (0), or the available memory is unknown (-1), collective, from the memory
  // real_memory of the setup of the approximation on the node: held now (checked is
  // nullptr), or estimated from a check (for the multigrid solver of a later operator,
  // before any coarse setup). It must fit next to the held approximation, and after the
  // setup (replacing it) with the reserve for the later allocations. The reserve is reduced
  // by the memory allocated since the check for an operator of the same size (for example
  // of a repeated eigenvalue solve, whose basis exists). Returns the largest estimate and
  // the smallest available memory over the nodes.
  memory_reporting::TrimHeap();
  const auto &node = memory_reporting::GetNodeComms(comm).node;
  long data[2] = {memory_reporting::GetCurrentMemory(),
                  memory_reporting::GetAvailableMemory()};
  MPI_Allreduce(MPI_IN_PLACE, &data[0], 1, MPI_LONG, MPI_SUM, node);
  MPI_Allreduce(MPI_IN_PLACE, &data[1], 1, MPI_LONG, MPI_MIN, node);
  const double used = data[0], avail = data[1];
  if (node_memory)
  {
    *node_memory = used;
  }
  double reserve;
  if (!checked)
  {
    reserve = coarse_memory_reserve * used;
  }
  else if (same_size)
  {
    const double later =
        std::max(used - (checked->check_memory - checked->real_memory), 0.0);
    reserve = std::max(coarse_memory_reserve * checked->check_memory - later, 0.0);
  }
  else
  {
    reserve = coarse_memory_reserve * (used + real_memory);
  }
  const double complex_memory =
      complex_coarse_memory_ratio * real_memory + coarse_memory_floor;
  int status[2] = {(avail >= 0.0 && complex_coarse_memory_ratio > 0.0 &&
                    complex_memory <= avail &&
                    complex_memory - (checked ? 0.0 : real_memory) + reserve <= avail),
                   (avail >= 0.0)};
  MPI_Allreduce(MPI_IN_PLACE, status, 2, MPI_INT, MPI_MIN, comm);
  double values[2] = {complex_memory, -avail};
  MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_DOUBLE, MPI_MAX, comm);
  if (estimate)
  {
    *estimate = values[0];
  }
  if (available)
  {
    *available = -values[1];
  }
  return !status[1] ? -1 : status[0];
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::TryComplexCoarseSolve(bool trial, double estimate,
                                                               double available)
{
  // Replace the coarse solver by the complex-valued one, set it up, and factorize it with
  // one solve (lazy solvers factorize at the first solve). A failure cannot be recovered
  // from on all processes: the constructors and the setup of the sparse direct solvers
  // communicate, and an exception (for example out of memory) may be raised on some
  // processes only while the others wait in a collective operation. A failure which is not
  // an exception (an abort, or the process killed when out of memory) does not return
  // either: print the hint first.
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;
  MPI_Comm comm = static_cast<const ParOperType *>(A[0])->GetComm();
  Mpi::Print(comm,
             " Setting up the complex-valued coarse solve: estimated {} per node ({} "
             "available); set config[\"Solver\"][\"Linear\"][\"ComplexCoarseSolve\"] "
             "to false if this runs out of memory\n",
             memory_reporting::FormatBytes(estimate),
             memory_reporting::FormatBytes(available));
  std::fflush(stdout);
  try
  {
    {
      perf::ScopedRole scoped_role(role);
      B[0] = complex_coarse_factory();
    }
    auto &wrapper = static_cast<MfemWrapperSolver<OperType> &>(*B[0]);
    if (trial)
    {
      // The counters and memory of the performance report wait for the outcome; the
      // factorization is timed as an explicit one.
      wrapper.DeferSetupRecord();
    }
    if (coarse_timer)
    {
      perf::ScopedRole scoped_role(role);
      wrapper.EnableTimer(true, Timer::KSP_COARSE_SOLVE);
    }
    B[0]->SetOperator(*A[0]);
    VecType x(A[0]->Height()), y(A[0]->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    linalg::SetRandom(comm, x);
    B[0]->SetInitialGuess(false);
    B[0]->Mult(x, y);
  }
  catch (const std::exception &e)
  {
    std::cerr << fmt::format(
                     "Error: the setup of the complex-valued coarse solve chosen at setup "
                     "failed ({}){}. Set config[\"Solver\"][\"Linear\"]"
                     "[\"ComplexCoarseSolve\"] to false to use the real-valued "
                     "approximation.",
                     (*e.what()) ? e.what() : "unknown error",
                     dynamic_cast<const std::bad_alloc *>(&e) ? ", out of memory" : "")
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::CheckCoarseSolve(long coarse_memory)
{
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // Only the real-valued approximation of an operator with real and imaginary parts is
    // inexact (see MfemWrapperSolver).
    const auto &Ac = static_cast<const ComplexParOperator &>(*A[0]);
    if (!complex_coarse_factory || coarse_checked || A.size() < 2 || !Ac.Real() ||
        !Ac.Imag())
    {
      return;
    }
    coarse_checked = true;
    MPI_Comm comm = Ac.GetComm();
    CoarseSolveChoice choice;
    choice.coarse_size = Ac.Height();
    Mpi::GlobalSum(1, &choice.coarse_size, comm);

    // GMRES probe of the coarse operator preconditioned by the approximation, from a random
    // unit start vector seeded by process, until the residual is below coarse_probe_max_res
    // (a run which fails on non-finite values is poor). The stopping tests are reduced over
    // the processes, so that all of them take the same decision.
    B[0]->SetInitialGuess(false);
    bool poor = true;
    int its = 0;
    double res = 1.0;
    {
      const int n = Ac.Height();
      VecType v0(n), t(n);
      v0.UseDevice(true);
      t.UseDevice(true);
      linalg::SetRandom(comm, v0);
      linalg::Normalize(comm, v0);
      ProbeGmres<VecType> run(comm, v0, coarse_probe_max_its);
      while (true)
      {
        run.Step(
            [&](const VecType &x, VecType &y)
            {
              Ac.Mult(x, t);
              B[0]->Mult(t, y);
            });
        bool flags[3] = {run.Stopped() && !run.Failed(), run.Failed(),
                         run.Residual() <= coarse_probe_max_res};
        Mpi::GlobalOr(3, flags, comm);
        run.SetFlags(flags[0], flags[1]);
        poor = flags[1] || !flags[2];
        if (flags[2] || run.Stopped() || run.Iterations() >= coarse_probe_max_its)
        {
          break;
        }
      }
      its = run.Iterations();
      res = run.Failed() ? std::nan("") : run.Residual();
    }
    int fits = 0;
    double estimate = 0.0, available = 0.0;
    if (poor)
    {
      // The memory of the setup of the approximation on the node, with the factorization
      // (here after the first solve of the probe for solvers which factorize lazily), as
      // measured by the solver (or the memory growth since the start of the coarse setup),
      // and at least coarse_memory_per_nonzero per nonzero of the real part.
      const auto &node = memory_reporting::GetNodeComms(comm).node;
      const auto *wrapper = dynamic_cast<const MfemWrapperSolver<OperType> *>(B[0].get());
      double real = (wrapper && wrapper->GetSetupMemory() >= 0.0)
                        ? wrapper->GetSetupMemory()
                        : memory_reporting::GetCurrentMemory() - coarse_memory;
      if (const auto *Ar_loc =
              dynamic_cast<const hypre::HypreCSRMatrix *>(Ac.LocalOperator().Real()))
      {
        real = std::max(real, coarse_memory_per_nonzero * Ar_loc->NNZ());
      }
      real = std::max(real, 0.0);
      MPI_Allreduce(MPI_IN_PLACE, &real, 1, MPI_DOUBLE, MPI_SUM, node);
      choice.real_memory = real;
      fits = CoarseMemoryFits(comm, real, nullptr, false, &estimate, &available,
                              &choice.check_memory);
    }
    std::vector<GmresProbeResult> cycle_probes;
    if (poor && fits > 0)
    {
      // Set up the complex-valued coarse solver next to the approximation, and compare
      // GMRES probes of the whole cycle with the approximation (candidate) and with the
      // complex-valued coarse solve (reference, see CompareWithGmresProbe) at the current
      // Chebyshev orders, as CheckSmootherOrders probes them: the complex-valued coarse
      // solve is used if it is predicted to need coarse_complex_min_ratio times fewer
      // iterations, and, if the orders were raised, also fewer than the approximation at
      // the base orders, which CheckSmootherOrders may fall back to.
      std::vector<int> orders, base_orders;
      const auto raised = GetRaisedSmootherLevels(orders, base_orders);
      auto other = std::move(B[0]);
      tried_complex_coarse = true;
      TryComplexCoarseSolve(true, estimate, available);
      std::swap(B[0], other);
      bool complex_selected = false, base_selected = false;
      auto Select = [&](bool complex, bool base)
      {
        if (complex != complex_selected)
        {
          std::swap(B[0], other);
          SetAuxiliarySmootherOrders(complex ? complex_coarse_aux_order : -1);
          complex_selected = complex;
        }
        if (base != base_selected)
        {
          SetSmootherOrders(raised, base ? base_orders : orders);
          base_selected = base;
        }
      };
      const auto &op = static_cast<const ComplexParOperator &>(*A.back());
      MPI_Comm comm_op = op.GetComm();
      const int n = op.Height();
      VecType v0(n), t(n);
      v0.UseDevice(true);
      t.UseDevice(true);
      linalg::SetRandom(comm_op, v0);
      linalg::Normalize(comm_op, v0);
      GmresProbeOptions opts;
      opts.max_ratio = coarse_complex_min_ratio;
      for (bool base : {false, true})
      {
        if (base && raised.empty())
        {
          break;
        }
        cycle_probes.push_back(CompareWithGmresProbe<VecType>(
            comm_op, v0, [&](bool complex) { Select(complex, base && !complex); },
            [&](const VecType &x, VecType &y)
            {
              op.Mult(x, t);
              Mult(t, y);
            },
            opts));
        if (!cycle_probes.back().use_reference)
        {
          break;
        }
      }
      complex_coarse = cycle_probes.back().use_reference;
      Select(complex_coarse, false);
      other.reset();
      if (complex_coarse)
      {
        static_cast<MfemWrapperSolver<OperType> &>(*B[0]).RecordDeferredSetup();
        probe_orders.clear();
      }
      memory_reporting::TrimHeap();
    }
    choice.complex = complex_coarse;
    coarse_choice = choice;
    std::string cycle;
    auto cycle_info = nlohmann::json::array();
    for (std::size_t k = 0; k < cycle_probes.size(); k++)
    {
      const auto &probe = cycle_probes[k];
      cycle += fmt::format(
          "{} {:.3e} after {:d} iterations{}, {:.3e} with the complex-valued "
          "coarse solve (predicted iteration ratio {:.2f})",
          (k == 0) ? ", multigrid GMRES probe residual" : ";", probe.res_candidate,
          probe.its, (k == 0) ? "" : " at the base Chebyshev order", probe.res_reference,
          probe.ratio);
      cycle_info.push_back(
          {{"BaseOrder", k > 0},
           {"Iterations", probe.its},
           {"Residual", probe.res_candidate},
           {"ComplexResidual", probe.res_reference},
           {"PredictedIterationRatio",
            std::isfinite(probe.ratio) ? nlohmann::json(probe.ratio) : nlohmann::json()}});
    }
    Mpi::Print(
        comm,
        " Real-valued coarse solve approximation: GMRES probe residual {:.3e} after "
        "{:d} iterations{}{}\n",
        res, its,
        !poor         ? ""
        : (fits < 0)  ? ", the available memory is unknown"
        : (fits == 0) ? fmt::format(", the complex-valued coarse solve is estimated not to "
                                    "fit in memory ({} per node, {} available)",
                                    memory_reporting::FormatBytes(estimate),
                                    memory_reporting::FormatBytes(available))
                      : cycle,
        cycle_probes.empty() ? ""
                             : (complex_coarse ? ": using the complex-valued one"
                                               : ": not using the complex-valued one"));
    perf::SetInfo(
        perf::LinearSolverPath(role, "CoarseSolveProbe"),
        {{"Iterations", its},
         {"Residual", std::isfinite(res) ? nlohmann::json(res) : nlohmann::json()},
         {"CycleProbes", cycle_info},
         {"Complex", complex_coarse}});
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetAuxiliarySmootherOrders(int order)
{
  // Give the auxiliary space smoothers which follow the order of the primary space their
  // own order, or (order <= 0) let them follow it again.
  if (complex_coarse_aux_order <= 0)
  {
    return;
  }
  aux_order_set.resize(A.size(), false);
  for (std::size_t l = 1; l < A.size(); l++)
  {
    auto *d = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get());
    if (!d || !dynamic_cast<const ChebyshevSmoother<OperType> *>(d->GetAuxiliarySmoother()))
    {
      continue;
    }
    if (order > 0 && d->GetAuxiliarySmootherOrder() <= 0)
    {
      d->SetAuxiliarySmootherOrder(order);
      aux_order_set[l] = true;
    }
    else if (order <= 0 && aux_order_set[l])
    {
      d->SetAuxiliarySmootherOrder(order);
      aux_order_set[l] = false;
    }
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetSliverPatches(double min_angle, int layers,
                                                          int max_size)
{
  MFEM_VERIFY(min_angle >= 0.0 && layers <= 1 && max_size > 0,
              "Invalid sliver patch parameters for GeometricMultigridSolver!");
  sliver_angle = min_angle;
  sliver_layers = layers;
  sliver_max_size = max_size;
  sliver_mesh = nullptr;
  for (auto &patch : patches)
  {
    patch.reset();
  }
  patches_checked.clear();
  patch_spaces.clear();
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetSliverPatchOperators()
{
  using ParOperType = std::conditional_t<std::is_same_v<OperType, ComplexOperator>,
                                         ComplexParOperator, ParOperator>;
  if (sliver_angle <= 0.0 || A.size() < 2)
  {
    return;
  }
  const auto t0 = std::chrono::steady_clock::now();
  auto Space = [this](std::size_t l) -> const FiniteElementSpace &
  { return static_cast<const ParOperType *>(A[l])->TrialFiniteElementSpace(); };
  MPI_Comm comm = Space(A.size() - 1).GetComm();
  patch_spaces.resize(A.size(), nullptr);
  patches_checked.resize(A.size(), false);
  bool print = false;
  for (std::size_t l = 1; l < A.size(); l++)
  {
    const auto &fespace = Space(l);
    if (patch_spaces[l] != &fespace)
    {
      // Find the slivers of the mesh once (the levels of a p-multigrid hierarchy share it).
      const auto &mesh = fespace.GetParMesh();
      if (&mesh != sliver_mesh)
      {
        SliverStatistics stats;
        slivers = FindSliverElements(mesh, sliver_angle, &stats);
        sliver_mesh = &mesh;
        sliver_count = stats.num_slivers;
        if (sliver_count > 0)
        {
          Mpi::Print(comm,
                     " Sliver patch solves: {:d} tetrahedra with a minimum dihedral angle "
                     "below {:g} degrees (smallest {:.3g})\n",
                     sliver_count, sliver_angle, stats.min_angle);
          if (!SliverPatchesSupported())
          {
            Mpi::Warning(comm, "Sliver patch solves are disabled with device and OpenMP "
                               "backends!\n");
          }
          print = true;
        }
      }
      patches[l].reset();
      patches_checked[l] = false;
      patch_layers.resize(A.size(), 0);
      patch_work.resize(A.size(), 0.0);
      if (sliver_count > 0)
      {
        auto MakePatches = [&](int layers)
        {
          patch_layers[l] = layers;
          return std::make_unique<SliverPatchSolver<OperType>>(
              fespace, slivers, layers,
              static_cast<const ParOperType *>(A[l])->GetEssentialTrueDofs(),
              sliver_max_size, std::is_same_v<OperType, ComplexOperator>);
        };
        // The dense work of an application, Σₚ nₚ² on the process of the patches, relative
        // to the smoothing work of the level on that process, its true dofs times the order
        // of the Chebyshev smoother: when the vertex-star patches are expensive (large
        // merged clusters of slivers), the level keeps no patch solves.
        const Solver<OperType> *smoother = B[l].get();
        if (const auto *d =
                dynamic_cast<const DistRelaxationSmoother<OperType> *>(smoother))
        {
          smoother = d->GetSmoother();
        }
        int order = 1;
        if (const auto *c = dynamic_cast<const ChebyshevSmoother<OperType> *>(smoother))
        {
          order = c->GetOrder();
        }
        else if (const auto *c =
                     dynamic_cast<const ChebyshevSmoother1stKind<OperType> *>(smoother))
        {
          order = c->GetOrder();
        }
        patches[l] = MakePatches((sliver_layers < 0) ? 1 : sliver_layers);
        patch_work[l] =
            patches[l]->PatchWork() / (std::max(fespace.GetTrueVSize(), 1) * double(order));
        Mpi::GlobalMax(1, &patch_work[l], comm);
        if (sliver_layers < 0 && patch_work[l] > max_patch_work_ratio)
        {
          patches[l].reset();
          patch_layers[l] = -1;
        }
      }
      patch_spaces[l] = &fespace;
    }
    if (patches[l] && patch_spaces[l] == &fespace && !patches_checked[l])
    {
      // Slivers are rare: patches covering a large fraction of a large level are not
      // slivers (for example, with too large a threshold for the mesh), and many large
      // patches are expensive and can be unstable for indefinite operators.
      long long data[2] = {patches[l]->NumTrueDofs(), 0};
      Mpi::GlobalSum(1, data, comm);
      data[1] = fespace.GlobalTrueVSize();
      patches_checked[l] = true;
      if (data[0] > max_patch_fraction * data[1] && data[0] > max_patch_dofs_small)
      {
        Mpi::Warning(
            comm,
            "Sliver patch solves disabled on multigrid level {:d}: the patches of the "
            "tetrahedra with a minimum dihedral angle below {:g} degrees cover {:.2f}% "
            "of the true dofs (more than {:g}%)!\n",
            l, sliver_angle, 100.0 * data[0] / data[1], 100.0 * max_patch_fraction);
        patches[l].reset();
      }
    }
    if (patches[l])
    {
      patches[l]->SetOperator(*A[l]);
    }
  }
  if (print)
  {
    // Summary of the first setup for the slivers of a mesh (global).
    double setup =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Mpi::GlobalMax(1, &setup, comm);
    auto info = nlohmann::json::array();
    for (std::size_t l = 1; l < A.size(); l++)
    {
      if (!patches[l] && patch_layers[l] < 0)
      {
        Mpi::Print(comm,
                   "  Level {:d} (p = {:d}): no patch solves (patch work ratio {:.3g}, "
                   "more than {:g})\n",
                   l, Space(l).GetMaxElementOrder(), patch_work[l], max_patch_work_ratio);
        info.push_back({{"Level", l}, {"Patches", 0}, {"WorkRatio", patch_work[l]}});
        continue;
      }
      long long data[7] = {0, 0, 0, 0, 0, 0, 0};
      int max_size = 0;
      if (patches[l])
      {
        data[0] = patches[l]->NumPatches();
        data[1] = patches[l]->NumTrueDofs();
        data[2] = patches[l]->NumSlivers();
        data[3] = patches[l]->NumRemoteTrueDofs();
        data[4] = (patches[l]->NumPatches() > 0 && !patches[l]->IsActive());
        data[5] = patches[l]->NumExcludedLDofs()[0];
        data[6] = patches[l]->NumExcludedLDofs()[1];
        max_size = patches[l]->MaxPatchSize();
      }
      Mpi::GlobalSum(7, data, comm);
      Mpi::GlobalMax(1, &max_size, comm);
      Mpi::Print(comm,
                 "  Level {:d} (p = {:d}): {:d} patches with {:d} true dofs (largest {:d}, "
                 "{:d} on other processes) for {:d} slivers{}\n",
                 l, Space(l).GetMaxElementOrder(), data[0], data[1], max_size, data[3],
                 data[2],
                 data[4] ? fmt::format(", inactive on {:d} processes", data[4]) : "");
      if (patches[l])
      {
        Mpi::Print(comm, "   Patch layers: {:d} (patch work ratio {:.3g})\n",
                   patch_layers[l], patch_work[l]);
      }
      if (data[5] > 0 || data[6] > 0)
      {
        Mpi::Print(comm, "   L-dofs left out: {:d} constrained, {:d} essential\n", data[5],
                   data[6]);
      }
      info.push_back({{"Level", l},
                      {"Patches", data[0]},
                      {"TrueDofs", data[1]},
                      {"MaxSize", max_size},
                      {"Slivers", data[2]},
                      {"RemoteTrueDofs", data[3]},
                      {"Inactive", data[4]},
                      {"Layers", patch_layers[l]},
                      {"WorkRatio", patch_work[l]}});
    }
    Mpi::Print(comm, "  Patch setup time: {:.3f} s (max)\n", setup);
    perf::SetInfo(perf::LinearSolverPath(role, "SliverPatches"), info);
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetPerfTallies(const std::string &prefix,
                                                        Timer::Index parent)
{
  // The coarse solve of the timed multigrid is charged to its own exclusive row.
  perf[0].pre = perf::RegisterTally(prefix + "/MG/Coarse",
                                    (use_timer && parent == Timer::KSP_PRECONDITIONER)
                                        ? Timer::KSP_COARSE_SOLVE
                                        : parent);
  for (std::size_t l = 1; l < B.size(); l++)
  {
    const auto lp = fmt::format("{}/MG/L{}/", prefix, l);
    perf[l].pre = perf::RegisterTally(lp + "PreSmooth", parent);
    perf[l].post = perf::RegisterTally(lp + "PostSmooth", parent);
    perf[l].res = perf::RegisterTally(lp + "Residual", parent);
    perf[l].restr = perf::RegisterTally(lp + "Restriction", parent);
    perf[l].prol = perf::RegisterTally(lp + "Prolongation", parent);
    perf[l].patch = perf::RegisterTally(lp + "SliverPatch", parent);
    B[l]->SetPerfTallies(lp + "Smoother", parent);
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::Mult(const VecType &x, VecType &y) const
{
  // Initialize.
  const auto n_levels = A.size();
  MFEM_ASSERT(!this->initial_guess,
              "Geometric multigrid solver does not use initial guess!");
  MFEM_ASSERT(n_levels > 1 || pc_it == 1,
              "Single-level geometric multigrid will not work with multiple iterations!");

  // Apply V-cycle. The initial guess for y is zero'd at the first pre-smooth iteration.
  X.back() = x;
  for (int it = 0; it < pc_it; it++)
  {
    VCycle(n_levels - 1, (it > 0));
  }
  y = Y.back();
}

namespace
{

inline void RealMult(const Operator &op, const Vector &x, Vector &y)
{
  op.Mult(x, y);
}

inline void RealMult(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.Mult(x.Real(), y.Real());
  op.Mult(x.Imag(), y.Imag());
}

inline void RealMultTranspose(const Operator &op, const Vector &x, Vector &y)
{
  op.MultTranspose(x, y);
}

inline void RealMultTranspose(const Operator &op, const ComplexVector &x, ComplexVector &y)
{
  op.MultTranspose(x.Real(), y.Real());
  op.MultTranspose(x.Imag(), y.Imag());
}

}  // namespace

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::VCycle(int l, bool initial_guess) const
{
  // Pre-smooth, with zero initial guess (Y = 0 set inside). This is the coarse solve at
  // level 0. Important to note that the smoothers must respect the initial guess flag
  // correctly (given X, Y, compute Y <- Y + B (X - A Y)) .
  B[l]->SetInitialGuess(initial_guess);
  const auto &t = perf[l];
  if (l == 0)
  {
    BlockTimer bt(Timer::KSP_COARSE_SOLVE, use_timer && !coarse_timer);
    perf::ScopedTally st(t.pre);
    perf::ScopedTally st_coarse(coarse_tally);
    B[l]->Mult(X[l], Y[l]);
    return;
  }
  {
    perf::ScopedTally st(t.pre);
    B[l]->Mult2(X[l], Y[l], R[l]);
  }

  // Compute residual.
  {
    perf::ScopedTally st(t.res);
    A[l]->Mult(Y[l], R[l]);
    linalg::AXPBY(1.0, X[l], -1.0, R[l]);
  }

  // Exact sliver patch solves, before the coarse-grid correction (complex operators, with
  // the residual updated) or for the same residual (real operators, for symmetry).
  if (patches[l] && patches[l]->IsActive())
  {
    perf::ScopedTally st(t.patch);
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      patches[l]->AddMultUpdate(R[l], Y[l]);
    }
    else
    {
      patches[l]->AddMult(R[l], Y[l]);
    }
  }

  // Coarse grid correction.
  {
    perf::ScopedTally st(t.restr);
    RealMultTranspose(*P[l - 1], R[l], X[l - 1]);
    if (dbc_tdof_lists[l - 1])
    {
      linalg::SetSubVector(X[l - 1], *dbc_tdof_lists[l - 1], 0.0);
    }
  }
  VCycle(l - 1, false);

  // Prolongate and add.
  {
    perf::ScopedTally st(t.prol);
    RealMult(*P[l - 1], Y[l - 1], R[l]);
    Y[l] += R[l];
  }

  // Post-smooth, with nonzero initial guess.
  perf::ScopedTally st(t.post);
  B[l]->SetInitialGuess(true);
  B[l]->MultTranspose2(X[l], Y[l], R[l]);
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetAuxiliaryResidual(
    MultigridAuxiliaryResidual residual)
{
  for (std::size_t l = 1; l < B.size(); l++)
  {
    if (auto *d = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get()))
    {
      d->SetAuxiliaryResidual(residual);
    }
  }
}

template <OperatorType OperType>
void GeometricMultigridSolver<OperType>::SetAuxiliaryGradientCheck(bool check)
{
  for (std::size_t l = 1; l < B.size(); l++)
  {
    if (auto *d = dynamic_cast<DistRelaxationSmoother<OperType> *>(B[l].get()))
    {
      d->SetAuxiliaryGradientCheck(check);
    }
  }
}

template class GeometricMultigridSolver<Operator>;
template class GeometricMultigridSolver<ComplexOperator>;

}  // namespace palace
