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

using json = nlohmann::json;

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

// Small helper: write `value` under `key` only if the user left it blank or wrote the
// explicit sentinel string "Default". Every other user entry passes through untouched.
// The sentinel handling is load-bearing — CheckConfiguration resolves enum DEFAULT to
// a concrete backend in-memory, and we have to propagate that resolution back into the
// JSON so the written config contains no defaults.
template <typename T>
void FillIfMissing(json &j, const char *key, const T &value)
{
  auto it = j.find(key);
  if (it == j.end())
  {
    j[key] = value;
    return;
  }
  if (it->is_string() && it->template get<std::string>() == "Default")
  {
    *it = value;
  }
}

// Helpers for each sub-structure. Each fills in any missing keys using the resolved
// IoData value; user-provided entries are left untouched.

void ConcretizeProblem(const config::ProblemData &problem, json &j_problem)
{
  FillIfMissing(j_problem, "Type", EnumString(problem.type));
  FillIfMissing(j_problem, "Verbose", problem.verbose);
  FillIfMissing(j_problem, "Output", problem.output);
  if (!j_problem.contains("OutputFormats"))
  {
    j_problem["OutputFormats"] = json::object();
  }
  auto &j_out = j_problem["OutputFormats"];
  FillIfMissing(j_out, "Paraview", problem.output_formats.paraview);
  FillIfMissing(j_out, "GridFunction", problem.output_formats.gridfunction);
}

void ConcretizeLinear(const config::LinearSolverData &linear, json &j_linear)
{
  FillIfMissing(j_linear, "Type", EnumString(linear.type));
  FillIfMissing(j_linear, "KSPType", EnumString(linear.krylov_solver));
  FillIfMissing(j_linear, "Tol", linear.tol);
  FillIfMissing(j_linear, "MaxIts", linear.max_it);
  FillIfMissing(j_linear, "MaxSize", linear.max_size);
  // Sentinel (-1) → concrete 0/1 int fields are declared as boolean in solver.json; cast
  // so the filled-in value matches the schema type.
  FillIfMissing(j_linear, "InitialGuess", static_cast<bool>(linear.initial_guess));
  FillIfMissing(j_linear, "MGMaxLevels", linear.mg_max_levels);
  FillIfMissing(j_linear, "MGCoarsenType", EnumString(linear.mg_coarsening));
  FillIfMissing(j_linear, "MGUseMesh", linear.mg_use_mesh);
  FillIfMissing(j_linear, "MGCycleIts", linear.mg_cycle_it);
  FillIfMissing(j_linear, "MGSmoothIts", linear.mg_smooth_it);
  FillIfMissing(j_linear, "MGSmoothOrder", linear.mg_smooth_order);
  FillIfMissing(j_linear, "MGSmoothEigScaleMax", linear.mg_smooth_sf_max);
  FillIfMissing(j_linear, "MGSmoothEigScaleMin", linear.mg_smooth_sf_min);
  FillIfMissing(j_linear, "MGSmoothChebyshev4th", linear.mg_smooth_cheby_4th);
  FillIfMissing(j_linear, "MGAuxiliarySmoother", static_cast<bool>(linear.mg_smooth_aux));
  FillIfMissing(j_linear, "PCMatReal", linear.pc_mat_real);
  FillIfMissing(j_linear, "PCMatShifted", static_cast<bool>(linear.pc_mat_shifted));
  FillIfMissing(j_linear, "ComplexCoarseSolve", linear.complex_coarse_solve);
  FillIfMissing(j_linear, "DropSmallEntries", linear.drop_small_entries);
  FillIfMissing(j_linear, "ReorderingReuse", linear.reorder_reuse);
  FillIfMissing(j_linear, "PCMatSymmetry", MatrixSymmetryString(linear.pc_mat_sym));
  FillIfMissing(j_linear, "PCSide", EnumString(linear.pc_side));
  FillIfMissing(j_linear, "ColumnOrdering", EnumString(linear.sym_factorization));
  FillIfMissing(j_linear, "STRUMPACKCompressionType",
                EnumString(linear.strumpack_compression_type));
  FillIfMissing(j_linear, "STRUMPACKCompressionTol", linear.strumpack_lr_tol);
  FillIfMissing(j_linear, "STRUMPACKLossyPrecision", linear.strumpack_lossy_precision);
  FillIfMissing(j_linear, "STRUMPACKButterflyLevels", linear.strumpack_butterfly_l);
  FillIfMissing(j_linear, "SuperLU3DCommunicator", linear.superlu_3d);
  FillIfMissing(j_linear, "AMSVectorInterpolation", linear.ams_vector_interp);
  FillIfMissing(j_linear, "AMSSingularOperator", static_cast<bool>(linear.ams_singular_op));
  FillIfMissing(j_linear, "AMGAggressiveCoarsening",
                static_cast<bool>(linear.amg_agg_coarsen));
  FillIfMissing(j_linear, "AMSMaxIts", linear.ams_max_it);
  FillIfMissing(j_linear, "DivFreeTol", linear.divfree_tol);
  FillIfMissing(j_linear, "DivFreeMaxIts", linear.divfree_max_it);
  FillIfMissing(j_linear, "EstimatorTol", linear.estimator_tol);
  FillIfMissing(j_linear, "EstimatorMaxIts", linear.estimator_max_it);
  FillIfMissing(j_linear, "EstimatorMG", linear.estimator_mg);
  FillIfMissing(j_linear, "GSOrthogonalization", EnumString(linear.gs_orthog));
}

void ConcretizeEigenmode(const config::EigenSolverData &eigenmode, json &j_eigen)
{
  FillIfMissing(j_eigen, "Target", eigenmode.target);
  FillIfMissing(j_eigen, "Tol", eigenmode.tol);
  FillIfMissing(j_eigen, "MaxIts", eigenmode.max_it);
  FillIfMissing(j_eigen, "MaxSize", eigenmode.max_size);
  FillIfMissing(j_eigen, "N", eigenmode.n);
  FillIfMissing(j_eigen, "Save", eigenmode.n_post);
  FillIfMissing(j_eigen, "Type", EnumString(eigenmode.type));
  FillIfMissing(j_eigen, "PEPLinear", eigenmode.pep_linear);
  FillIfMissing(j_eigen, "Scaling", eigenmode.scale);
  FillIfMissing(j_eigen, "StartVector", eigenmode.init_v0);
  FillIfMissing(j_eigen, "StartVectorConstant", eigenmode.init_v0_const);
  FillIfMissing(j_eigen, "MassOrthogonal", eigenmode.mass_orthog);
  FillIfMissing(j_eigen, "NonlinearType", EnumString(eigenmode.nonlinear_type));
  FillIfMissing(j_eigen, "RefineNonlinear", eigenmode.refine_nonlinear);
  FillIfMissing(j_eigen, "LinearTol", eigenmode.linear_tol);
  FillIfMissing(j_eigen, "TargetUpper", eigenmode.target_upper);
  FillIfMissing(j_eigen, "PreconditionerLag", eigenmode.preconditioner_lag);
  FillIfMissing(j_eigen, "PreconditionerLagTol", eigenmode.preconditioner_lag_tol);
  FillIfMissing(j_eigen, "MaxRestart", eigenmode.max_restart);
}

void ConcretizeTransient(const config::TransientSolverData &transient, json &j_transient)
{
  FillIfMissing(j_transient, "Type", EnumString(transient.type));
  FillIfMissing(j_transient, "Excitation", EnumString(transient.excitation));
  FillIfMissing(j_transient, "ExcitationFreq", transient.pulse_f);
  FillIfMissing(j_transient, "ExcitationWidth", transient.pulse_tau);
  FillIfMissing(j_transient, "MaxTime", transient.max_t);
  FillIfMissing(j_transient, "TimeStep", transient.delta_t);
  FillIfMissing(j_transient, "SaveStep", transient.delta_post);
  FillIfMissing(j_transient, "Order", transient.order);
  FillIfMissing(j_transient, "RelTol", transient.rel_tol);
  FillIfMissing(j_transient, "AbsTol", transient.abs_tol);
}

void ConcretizeDriven(const config::DrivenSolverData &driven, json &j_driven)
{
  FillIfMissing(j_driven, "Restart", driven.restart);
  FillIfMissing(j_driven, "AdaptiveTol", driven.adaptive_tol);
  FillIfMissing(j_driven, "AdaptiveMaxSamples", driven.adaptive_max_size);
  FillIfMissing(j_driven, "AdaptiveConvergenceMemory", driven.adaptive_memory);
  FillIfMissing(j_driven, "AdaptiveGSOrthogonalization",
                EnumString(driven.adaptive_solver_gs_orthog_type));
  FillIfMissing(j_driven, "AdaptiveCircuitSynthesis", driven.adaptive_circuit_synthesis);
  FillIfMissing(j_driven, "AdaptiveCircuitSynthesisDomainOrthogonalization",
                EnumString(driven.adaptive_circuit_synthesis_domain_orthog));
}

}  // namespace

json IoData::ConcretizeDefaults(const IoData &iodata, json config)
{
  // Walk the resolved IoData and fill in any keys the user left out, so the returned
  // config is enough to re-run the simulation without consulting any Palace defaults.
  // User-provided entries are preserved byte-for-byte.

  if (!config.contains("Problem"))
  {
    config["Problem"] = json::object();
  }
  ConcretizeProblem(iodata.problem, config["Problem"]);

  if (!config.contains("Solver"))
  {
    config["Solver"] = json::object();
  }
  auto &j_solver = config["Solver"];

  FillIfMissing(j_solver, "Order", iodata.solver.order);
  FillIfMissing(j_solver, "PartialAssemblyOrder", iodata.solver.pa_order_threshold);
  FillIfMissing(j_solver, "QuadratureOrderJacobian", iodata.solver.q_order_jac);
  FillIfMissing(j_solver, "QuadratureOrderExtra", iodata.solver.q_order_extra);
  FillIfMissing(j_solver, "Device", EnumString(iodata.solver.device));
  FillIfMissing(j_solver, "Backend", iodata.solver.ceed_backend);

  if (!j_solver.contains("Linear"))
  {
    j_solver["Linear"] = json::object();
  }
  ConcretizeLinear(iodata.solver.linear, j_solver["Linear"]);

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
      FillIfMissing(j_solver["Electrostatic"], "Save", iodata.solver.electrostatic.n_post);
      break;
    case ProblemType::MAGNETOSTATIC:
      if (!j_solver.contains("Magnetostatic"))
      {
        j_solver["Magnetostatic"] = json::object();
      }
      FillIfMissing(j_solver["Magnetostatic"], "Save", iodata.solver.magnetostatic.n_post);
      break;
  }

  // Boundaries.WavePort — fill in the resolved eigen solver backend for each port that
  // did not specify SolverType.
  if (config.contains("Boundaries"))
  {
    auto &j_boundaries = config["Boundaries"];
    if (j_boundaries.contains("WavePort"))
    {
      for (auto &j_port : j_boundaries["WavePort"])
      {
        if (j_port.contains("SolverType"))
        {
          continue;
        }
        int idx = j_port.at("Index").get<int>();
        auto it = iodata.boundaries.waveport.find(idx);
        if (it != iodata.boundaries.waveport.end())
        {
          j_port["SolverType"] = EnumString(it->second.eigen_solver);
        }
      }
    }
  }

  return config;
}

}  // namespace palace
