// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iodata.hpp"

#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
#include "utils/enum_string.hpp"

namespace palace
{

namespace
{

using json = nlohmann::json;
using Entry = std::pair<std::string, json>;

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

// Write `value` under `key` only when the user left it blank or wrote the explicit
// sentinel string "Default". Every other user entry passes through untouched.
// Sentinel handling is load-bearing — CheckConfiguration resolves enum DEFAULT to a
// concrete backend in-memory, and we propagate that resolution back into the JSON so
// the written config contains no defaults.
void Concretize(json &j, const std::string &key, const json &value)
{
  auto it = j.find(key);
  if (it == j.end())
  {
    j[key] = value;
    return;
  }
  if (it->is_string() && it->get<std::string>() == "Default")
  {
    *it = value;
  }
}

void ApplyEntries(json &j, const std::vector<Entry> &entries)
{
  for (const auto &[key, value] : entries)
  {
    Concretize(j, key, value);
  }
}

void ConcretizeProblem(const config::ProblemData &problem, json &j_problem)
{
  ApplyEntries(j_problem, {{"Type", EnumString(problem.type)},
                           {"Verbose", problem.verbose},
                           {"Output", problem.output}});
  if (!j_problem.contains("OutputFormats"))
  {
    j_problem["OutputFormats"] = json::object();
  }
  ApplyEntries(j_problem["OutputFormats"],
               {{"Paraview", problem.output_formats.paraview},
                {"GridFunction", problem.output_formats.gridfunction}});
}

void ConcretizeLinear(const config::LinearSolverData &linear, json &j_linear)
{
  // Sentinel (-1) → concrete 0/1 int fields are declared as boolean in solver.json; cast
  // so the filled-in value matches the schema type.
  ApplyEntries(j_linear,
               {{"Type", EnumString(linear.type)},
                {"KSPType", EnumString(linear.krylov_solver)},
                {"Tol", linear.tol},
                {"MaxIts", linear.max_it},
                {"MaxSize", linear.max_size},
                {"InitialGuess", static_cast<bool>(linear.initial_guess)},
                {"MGMaxLevels", linear.mg_max_levels},
                {"MGCoarsenType", EnumString(linear.mg_coarsening)},
                {"MGUseMesh", linear.mg_use_mesh},
                {"MGCycleIts", linear.mg_cycle_it},
                {"MGSmoothIts", linear.mg_smooth_it},
                {"MGSmoothOrder", linear.mg_smooth_order},
                {"MGSmoothEigScaleMax", linear.mg_smooth_sf_max},
                {"MGSmoothEigScaleMin", linear.mg_smooth_sf_min},
                {"MGSmoothChebyshev4th", linear.mg_smooth_cheby_4th},
                {"MGAuxiliarySmoother", static_cast<bool>(linear.mg_smooth_aux)},
                {"PCMatReal", linear.pc_mat_real},
                {"PCMatShifted", static_cast<bool>(linear.pc_mat_shifted)},
                {"ComplexCoarseSolve", linear.complex_coarse_solve},
                {"DropSmallEntries", linear.drop_small_entries},
                {"ReorderingReuse", linear.reorder_reuse},
                {"PCMatSymmetry", MatrixSymmetryString(linear.pc_mat_sym)},
                {"PCSide", EnumString(linear.pc_side)},
                {"ColumnOrdering", EnumString(linear.sym_factorization)},
                {"STRUMPACKCompressionType", EnumString(linear.strumpack_compression_type)},
                {"STRUMPACKCompressionTol", linear.strumpack_lr_tol},
                {"STRUMPACKLossyPrecision", linear.strumpack_lossy_precision},
                {"STRUMPACKButterflyLevels", linear.strumpack_butterfly_l},
                {"SuperLU3DCommunicator", linear.superlu_3d},
                {"AMSVectorInterpolation", linear.ams_vector_interp},
                {"AMSSingularOperator", static_cast<bool>(linear.ams_singular_op)},
                {"AMGAggressiveCoarsening", static_cast<bool>(linear.amg_agg_coarsen)},
                {"AMSMaxIts", linear.ams_max_it},
                {"DivFreeTol", linear.divfree_tol},
                {"DivFreeMaxIts", linear.divfree_max_it},
                {"EstimatorTol", linear.estimator_tol},
                {"EstimatorMaxIts", linear.estimator_max_it},
                {"EstimatorMG", linear.estimator_mg},
                {"GSOrthogonalization", EnumString(linear.gs_orthog)}});
}

void ConcretizeEigenmode(const config::EigenSolverData &eigenmode, json &j_eigen)
{
  ApplyEntries(j_eigen, {{"Target", eigenmode.target},
                         {"Tol", eigenmode.tol},
                         {"MaxIts", eigenmode.max_it},
                         {"MaxSize", eigenmode.max_size},
                         {"N", eigenmode.n},
                         {"Save", eigenmode.n_post},
                         {"Type", EnumString(eigenmode.type)},
                         {"PEPLinear", eigenmode.pep_linear},
                         {"Scaling", eigenmode.scale},
                         {"StartVector", eigenmode.init_v0},
                         {"StartVectorConstant", eigenmode.init_v0_const},
                         {"MassOrthogonal", eigenmode.mass_orthog},
                         {"NonlinearType", EnumString(eigenmode.nonlinear_type)},
                         {"RefineNonlinear", eigenmode.refine_nonlinear},
                         {"LinearTol", eigenmode.linear_tol},
                         {"TargetUpper", eigenmode.target_upper},
                         {"PreconditionerLag", eigenmode.preconditioner_lag},
                         {"PreconditionerLagTol", eigenmode.preconditioner_lag_tol},
                         {"MaxRestart", eigenmode.max_restart}});
}

void ConcretizeTransient(const config::TransientSolverData &transient, json &j_transient)
{
  ApplyEntries(j_transient, {{"Type", EnumString(transient.type)},
                             {"Excitation", EnumString(transient.excitation)},
                             {"ExcitationFreq", transient.pulse_f},
                             {"ExcitationWidth", transient.pulse_tau},
                             {"MaxTime", transient.max_t},
                             {"TimeStep", transient.delta_t},
                             {"SaveStep", transient.delta_post},
                             {"Order", transient.order},
                             {"RelTol", transient.rel_tol},
                             {"AbsTol", transient.abs_tol}});
}

void ConcretizeDriven(const config::DrivenSolverData &driven, json &j_driven)
{
  ApplyEntries(j_driven, {{"Restart", driven.restart},
                          {"AdaptiveTol", driven.adaptive_tol},
                          {"AdaptiveMaxSamples", driven.adaptive_max_size},
                          {"AdaptiveConvergenceMemory", driven.adaptive_memory},
                          {"AdaptiveGSOrthogonalization",
                           EnumString(driven.adaptive_solver_gs_orthog_type)},
                          {"AdaptiveCircuitSynthesis", driven.adaptive_circuit_synthesis},
                          {"AdaptiveCircuitSynthesisDomainOrthogonalization",
                           EnumString(driven.adaptive_circuit_synthesis_domain_orthog)}});
}

}  // namespace

json IoData::ConcretizeDefaults(const IoData &iodata, json config)
{
  // Walk the resolved IoData and fill in any keys the user left out, so the returned
  // config is enough to re-run the simulation without consulting any Palace defaults.
  // User-provided entries are preserved byte-for-byte; user-written "Default" sentinels
  // are replaced with the concrete value from the IoData.

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

  ApplyEntries(j_solver, {{"Order", iodata.solver.order},
                          {"PartialAssemblyOrder", iodata.solver.pa_order_threshold},
                          {"QuadratureOrderJacobian", iodata.solver.q_order_jac},
                          {"QuadratureOrderExtra", iodata.solver.q_order_extra},
                          {"Device", EnumString(iodata.solver.device)},
                          {"Backend", iodata.solver.ceed_backend}});

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
      Concretize(j_solver["Electrostatic"], "Save", iodata.solver.electrostatic.n_post);
      break;
    case ProblemType::MAGNETOSTATIC:
      if (!j_solver.contains("Magnetostatic"))
      {
        j_solver["Magnetostatic"] = json::object();
      }
      Concretize(j_solver["Magnetostatic"], "Save", iodata.solver.magnetostatic.n_post);
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
        int idx = j_port.at("Index").get<int>();
        auto it = iodata.boundaries.waveport.find(idx);
        if (it != iodata.boundaries.waveport.end())
        {
          Concretize(j_port, "SolverType", EnumString(it->second.eigen_solver));
        }
      }
    }
  }

  return config;
}

}  // namespace palace
