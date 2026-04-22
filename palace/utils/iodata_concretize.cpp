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

}  // namespace

using json = nlohmann::json;

void IoData::ConcretizeDefaults(const IoData &iodata, json &config)
{
  // Write resolved IoData values back into the JSON config. This produces a self-describing
  // config with all sentinels replaced by the concrete values that CheckConfiguration()
  // resolved — a complete record of every Palace decision for a given run. Called BEFORE
  // nondimensionalization so all values are in SI units.

  // Solver section — always present after CheckConfiguration (even if omitted by the user,
  // defaults are resolved).
  if (!config.contains("Solver"))
  {
    config["Solver"] = json::object();
  }
  auto &j_solver = config["Solver"];

  // Top-level Solver fields.
  j_solver["Order"] = iodata.solver.order;

  // Solver.Linear — the bulk of resolved sentinels live here.
  if (!j_solver.contains("Linear"))
  {
    j_solver["Linear"] = json::object();
  }
  auto &j_linear = j_solver["Linear"];
  const auto &linear = iodata.solver.linear;

  j_linear["Type"] = EnumString(linear.type);
  j_linear["KSPType"] = EnumString(linear.krylov_solver);
  j_linear["Tol"] = linear.tol;
  j_linear["MaxIts"] = linear.max_it;
  j_linear["MaxSize"] = linear.max_size;
  // Sentinel (-1) → concrete 0/1 int fields are declared as boolean in solver.json; emit
  // as bool so the resolved config round-trips through schema validation cleanly.
  j_linear["InitialGuess"] = static_cast<bool>(linear.initial_guess);
  j_linear["PCMatShifted"] = static_cast<bool>(linear.pc_mat_shifted);
  j_linear["MGAuxiliarySmoother"] = static_cast<bool>(linear.mg_smooth_aux);
  j_linear["MGSmoothOrder"] = linear.mg_smooth_order;
  j_linear["AMSSingularOperator"] = static_cast<bool>(linear.ams_singular_op);
  j_linear["AMGAggressiveCoarsening"] = static_cast<bool>(linear.amg_agg_coarsen);
  j_linear["AMSMaxIts"] = linear.ams_max_it;
  j_linear["MGCycleIts"] = linear.mg_cycle_it;
  j_linear["PCMatSymmetry"] = MatrixSymmetryString(linear.pc_mat_sym);
  j_linear["ReorderingReuse"] = linear.reorder_reuse;
  j_linear["PCSide"] = EnumString(linear.pc_side);
  j_linear["ColumnOrdering"] = EnumString(linear.sym_factorization);
  j_linear["GSOrthogonalization"] = EnumString(linear.gs_orthog);

  // Solver.Eigenmode — resolved eigen solver backend and associated defaults.
  if (iodata.problem.type == ProblemType::EIGENMODE)
  {
    if (!j_solver.contains("Eigenmode"))
    {
      j_solver["Eigenmode"] = json::object();
    }
    auto &j_eigen = j_solver["Eigenmode"];
    const auto &eigenmode = iodata.solver.eigenmode;
    j_eigen["Type"] = EnumString(eigenmode.type);
    j_eigen["Tol"] = eigenmode.tol;
    j_eigen["MaxIts"] = eigenmode.max_it;
    j_eigen["MaxSize"] = eigenmode.max_size;
    j_eigen["N"] = eigenmode.n;
    j_eigen["Save"] = eigenmode.n_post;
    j_eigen["PEPLinear"] = eigenmode.pep_linear;
    j_eigen["Scaling"] = eigenmode.scale;
    j_eigen["StartVector"] = eigenmode.init_v0;
    j_eigen["StartVectorConstant"] = eigenmode.init_v0_const;
    j_eigen["MassOrthogonal"] = eigenmode.mass_orthog;
  }

  // Solver.Transient — DEFAULT already resolved to the concrete scheme in
  // CheckConfiguration, so EnumString produces the concrete name.
  if (iodata.problem.type == ProblemType::TRANSIENT)
  {
    if (!j_solver.contains("Transient"))
    {
      j_solver["Transient"] = json::object();
    }
    auto &j_transient = j_solver["Transient"];
    const auto &transient = iodata.solver.transient;
    j_transient["Type"] = EnumString(transient.type);
    j_transient["Order"] = transient.order;
    j_transient["RelTol"] = transient.rel_tol;
    j_transient["AbsTol"] = transient.abs_tol;
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
