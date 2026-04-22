// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iodata.hpp"

#include <sstream>
#include <string>
#include <nlohmann/json.hpp>
#include "utils/enum_string.hpp"

namespace palace
{

namespace
{

template <typename E>
std::string EnumString(E e)
{
  std::ostringstream oss;
  oss << e;
  return oss.str();
}

// MatrixSymmetry has no PALACE_JSON_SERIALIZE_ENUM (it is always computed, never
// deserialized from JSON), so provide a local stringifier.
std::string MatrixSymmetryString(MatrixSymmetry s)
{
  switch (s)
  {
    case MatrixSymmetry::SPD:
      return "SPD";
    case MatrixSymmetry::SYMMETRIC:
      return "Symmetric";
    case MatrixSymmetry::UNSYMMETRIC:
      return "Unsymmetric";
  }
  return "Unsymmetric";
}

using json = nlohmann::json;

// Helpers for each sub-structure. Each helper writes every field of the struct back into
// the JSON object. Fields left at struct defaults are still written so that the
// `-resolved.json` is a complete record of the run, independent of which user-facing
// defaults happen to match struct defaults today.

void ConcretizeProblem(const config::ProblemData &problem, json &j_problem)
{
  j_problem["Type"] = EnumString(problem.type);
  j_problem["Verbose"] = problem.verbose;
  j_problem["Output"] = problem.output;
  auto &j_out = j_problem["OutputFormats"];
  if (!j_out.is_object())
  {
    j_out = json::object();
  }
  j_out["Paraview"] = problem.output_formats.paraview;
  j_out["GridFunction"] = problem.output_formats.gridfunction;
}

void ConcretizeLinear(const config::LinearSolverData &linear, json &j_linear)
{
  j_linear["Type"] = EnumString(linear.type);
  j_linear["KSPType"] = EnumString(linear.krylov_solver);
  j_linear["Tol"] = linear.tol;
  j_linear["MaxIts"] = linear.max_it;
  j_linear["MaxSize"] = linear.max_size;
  // Sentinel (-1) → concrete 0/1 int fields are declared as boolean in solver.json; emit
  // as bool so the resolved config round-trips through schema validation cleanly.
  j_linear["InitialGuess"] = static_cast<bool>(linear.initial_guess);
  j_linear["MGMaxLevels"] = linear.mg_max_levels;
  j_linear["MGCoarsenType"] = EnumString(linear.mg_coarsening);
  j_linear["MGUseMesh"] = linear.mg_use_mesh;
  j_linear["MGCycleIts"] = linear.mg_cycle_it;
  j_linear["MGSmoothIts"] = linear.mg_smooth_it;
  j_linear["MGSmoothOrder"] = linear.mg_smooth_order;
  j_linear["MGSmoothEigScaleMax"] = linear.mg_smooth_sf_max;
  j_linear["MGSmoothEigScaleMin"] = linear.mg_smooth_sf_min;
  j_linear["MGSmoothChebyshev4th"] = linear.mg_smooth_cheby_4th;
  j_linear["MGAuxiliarySmoother"] = static_cast<bool>(linear.mg_smooth_aux);
  j_linear["PCMatReal"] = linear.pc_mat_real;
  j_linear["PCMatShifted"] = static_cast<bool>(linear.pc_mat_shifted);
  j_linear["ComplexCoarseSolve"] = linear.complex_coarse_solve;
  j_linear["DropSmallEntries"] = linear.drop_small_entries;
  j_linear["ReorderingReuse"] = linear.reorder_reuse;
  j_linear["PCMatSymmetry"] = MatrixSymmetryString(linear.pc_mat_sym);
  j_linear["PCSide"] = EnumString(linear.pc_side);
  j_linear["ColumnOrdering"] = EnumString(linear.sym_factorization);
  j_linear["STRUMPACKCompressionType"] = EnumString(linear.strumpack_compression_type);
  j_linear["STRUMPACKCompressionTol"] = linear.strumpack_lr_tol;
  j_linear["STRUMPACKLossyPrecision"] = linear.strumpack_lossy_precision;
  j_linear["STRUMPACKButterflyLevels"] = linear.strumpack_butterfly_l;
  j_linear["SuperLU3DCommunicator"] = linear.superlu_3d;
  j_linear["AMSVectorInterpolation"] = linear.ams_vector_interp;
  j_linear["AMSSingularOperator"] = static_cast<bool>(linear.ams_singular_op);
  j_linear["AMGAggressiveCoarsening"] = static_cast<bool>(linear.amg_agg_coarsen);
  j_linear["AMSMaxIts"] = linear.ams_max_it;
  j_linear["DivFreeTol"] = linear.divfree_tol;
  j_linear["DivFreeMaxIts"] = linear.divfree_max_it;
  j_linear["EstimatorTol"] = linear.estimator_tol;
  j_linear["EstimatorMaxIts"] = linear.estimator_max_it;
  j_linear["EstimatorMG"] = linear.estimator_mg;
  j_linear["GSOrthogonalization"] = EnumString(linear.gs_orthog);
}

void ConcretizeEigenmode(const config::EigenSolverData &eigenmode, json &j_eigen)
{
  j_eigen["Target"] = eigenmode.target;
  j_eigen["Tol"] = eigenmode.tol;
  j_eigen["MaxIts"] = eigenmode.max_it;
  j_eigen["MaxSize"] = eigenmode.max_size;
  j_eigen["N"] = eigenmode.n;
  j_eigen["Save"] = eigenmode.n_post;
  j_eigen["Type"] = EnumString(eigenmode.type);
  j_eigen["PEPLinear"] = eigenmode.pep_linear;
  j_eigen["Scaling"] = eigenmode.scale;
  j_eigen["StartVector"] = eigenmode.init_v0;
  j_eigen["StartVectorConstant"] = eigenmode.init_v0_const;
  j_eigen["MassOrthogonal"] = eigenmode.mass_orthog;
  j_eigen["NonlinearType"] = EnumString(eigenmode.nonlinear_type);
  j_eigen["RefineNonlinear"] = eigenmode.refine_nonlinear;
  j_eigen["LinearTol"] = eigenmode.linear_tol;
  j_eigen["TargetUpper"] = eigenmode.target_upper;
  j_eigen["PreconditionerLag"] = eigenmode.preconditioner_lag;
  j_eigen["PreconditionerLagTol"] = eigenmode.preconditioner_lag_tol;
  j_eigen["MaxRestart"] = eigenmode.max_restart;
}

void ConcretizeTransient(const config::TransientSolverData &transient, json &j_transient)
{
  j_transient["Type"] = EnumString(transient.type);
  j_transient["Excitation"] = EnumString(transient.excitation);
  j_transient["ExcitationFreq"] = transient.pulse_f;
  j_transient["ExcitationWidth"] = transient.pulse_tau;
  j_transient["MaxTime"] = transient.max_t;
  j_transient["TimeStep"] = transient.delta_t;
  j_transient["SaveStep"] = transient.delta_post;
  j_transient["Order"] = transient.order;
  j_transient["RelTol"] = transient.rel_tol;
  j_transient["AbsTol"] = transient.abs_tol;
}

void ConcretizeDriven(const config::DrivenSolverData &driven, json &j_driven)
{
  j_driven["Restart"] = driven.restart;
  j_driven["AdaptiveTol"] = driven.adaptive_tol;
  j_driven["AdaptiveMaxSamples"] = driven.adaptive_max_size;
  j_driven["AdaptiveConvergenceMemory"] = driven.adaptive_memory;
  j_driven["AdaptiveGSOrthogonalization"] =
      EnumString(driven.adaptive_solver_gs_orthog_type);
  j_driven["AdaptiveCircuitSynthesis"] = driven.adaptive_circuit_synthesis;
  j_driven["AdaptiveCircuitSynthesisDomainOrthogonalization"] =
      EnumString(driven.adaptive_circuit_synthesis_domain_orthog);
  // Note: `Samples`, `Save` (save_indices), and the back-compat MinFreq/MaxFreq/FreqStep
  // block are left as the user wrote them — they are the primary user input, not defaults.
}

}  // namespace

void IoData::ConcretizeDefaults(const IoData &iodata, json &config)
{
  // Write resolved IoData values back into the JSON config. Produces a self-describing
  // record of every Palace decision — every sentinel replaced with the concrete value
  // CheckConfiguration resolved — so the `-resolved.json` is enough to re-run the
  // simulation without needing to consult any defaults. Called before nondimensionalization
  // so values stay in SI units.

  // Problem section.
  if (!config.contains("Problem"))
  {
    config["Problem"] = json::object();
  }
  ConcretizeProblem(iodata.problem, config["Problem"]);

  // Solver section and its sub-structures.
  if (!config.contains("Solver"))
  {
    config["Solver"] = json::object();
  }
  auto &j_solver = config["Solver"];

  j_solver["Order"] = iodata.solver.order;
  j_solver["PartialAssemblyOrder"] = iodata.solver.pa_order_threshold;
  j_solver["QuadratureOrderJacobian"] = iodata.solver.q_order_jac;
  j_solver["QuadratureOrderExtra"] = iodata.solver.q_order_extra;
  j_solver["Device"] = EnumString(iodata.solver.device);
  j_solver["Backend"] = iodata.solver.ceed_backend;

  if (!j_solver.contains("Linear"))
  {
    j_solver["Linear"] = json::object();
  }
  ConcretizeLinear(iodata.solver.linear, j_solver["Linear"]);

  // Problem-type-specific solver sections.
  switch (iodata.problem.type)
  {
    case ProblemType::EIGENMODE:
      if (!j_solver.contains("Eigenmode"))
      {
        j_solver["Eigenmode"] = json::object();
      }
      ConcretizeEigenmode(iodata.solver.eigenmode, j_solver["Eigenmode"]);
      break;
    case ProblemType::TRANSIENT:
      if (!j_solver.contains("Transient"))
      {
        j_solver["Transient"] = json::object();
      }
      ConcretizeTransient(iodata.solver.transient, j_solver["Transient"]);
      break;
    case ProblemType::DRIVEN:
      if (!j_solver.contains("Driven"))
      {
        j_solver["Driven"] = json::object();
      }
      ConcretizeDriven(iodata.solver.driven, j_solver["Driven"]);
      break;
    case ProblemType::ELECTROSTATIC:
      if (!j_solver.contains("Electrostatic"))
      {
        j_solver["Electrostatic"] = json::object();
      }
      j_solver["Electrostatic"]["Save"] = iodata.solver.electrostatic.n_post;
      break;
    case ProblemType::MAGNETOSTATIC:
      if (!j_solver.contains("Magnetostatic"))
      {
        j_solver["Magnetostatic"] = json::object();
      }
      j_solver["Magnetostatic"]["Save"] = iodata.solver.magnetostatic.n_post;
      break;
  }

  // Boundaries.WavePort — resolved eigen solver backend for each port.
  if (config.contains("Boundaries"))
  {
    auto &j_boundaries = config["Boundaries"];
    if (j_boundaries.contains("WavePort"))
    {
      for (auto &j_port : j_boundaries["WavePort"])
      {
        int idx = j_port.at("Index").get<int>();
        auto it = iodata.boundaries.waveport.find(idx);
        if (it != iodata.boundaries.waveport.end())
        {
          j_port["SolverType"] = EnumString(it->second.eigen_solver);
        }
      }
    }
  }
}

}  // namespace palace
