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

// Fill the key if absent, or if the user wrote the explicit "Default" sentinel.
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
                             {"SaveStep", transient.delta_post}});
  // Order and RelTol/AbsTol are only meaningful for the adaptive CVODE/ARKODE schemes.
  // GeneralizedAlpha and RungeKutta warn at parse time if these keys are present, so we
  // must only emit them for schemes that actually consume them.
  const bool adaptive = (transient.type == TimeSteppingScheme::CVODE ||
                         transient.type == TimeSteppingScheme::ARKODE);
  if (adaptive)
  {
    ApplyEntries(j_transient, {{"Order", transient.order},
                               {"RelTol", transient.rel_tol},
                               {"AbsTol", transient.abs_tol}});
  }
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

void ConcretizeBoundaryMode(const config::BoundaryModeSolverData &bm, json &j_bm)
{
  ApplyEntries(j_bm, {{"Freq", bm.freq},
                      {"N", bm.n},
                      {"Save", bm.n_post},
                      {"Target", bm.target},
                      {"Tol", bm.tol},
                      {"MaxSize", bm.max_size},
                      {"Type", EnumString(bm.type)}});
}

void ConcretizeModel(const config::ModelData &model, json &j_model)
{
  // Note: Lc is deliberately NOT emitted. It defaults to a sentinel (-1) that
  // NondimensionalizeInputs() resolves from the mesh bounding box at runtime, after
  // ConcretizeDefaults has already produced the sidecar. Emitting the sentinel would
  // fail schema validation (Lc > 0); emitting a concrete value is impossible here
  // because the mesh has not been loaded yet. If the user wrote Lc explicitly, the
  // original key already survives untouched.
  ApplyEntries(j_model, {{"L0", model.L0},
                         {"RemoveCurvature", model.remove_curvature},
                         {"MakeSimplex", model.make_simplex},
                         {"MakeHexahedral", model.make_hex},
                         {"ReorderElements", model.reorder_elements},
                         {"CleanUnusedElements", model.clean_unused_elements},
                         {"CrackInternalBoundaryElements", model.crack_bdr_elements},
                         {"RefineCrackElements", model.refine_crack_elements},
                         {"CrackDisplacementFactor", model.crack_displ_factor},
                         {"AddInterfaceBoundaryElements", model.add_bdr_elements},
                         {"ExportPrerefinedMesh", model.export_prerefined_mesh},
                         {"ReorientTetMesh", model.reorient_tet_mesh},
                         {"Partitioning", model.partitioning}});
  if (!j_model.contains("Refinement"))
  {
    j_model["Refinement"] = json::object();
  }
  const auto &ref = model.refinement;
  ApplyEntries(j_model["Refinement"],
               {{"Tol", ref.tol},
                {"MaxIts", ref.max_it},
                {"MaxSize", ref.max_size},
                {"Nonconformal", ref.nonconformal},
                {"MaxNCLevels", ref.max_nc_levels},
                {"UpdateFraction", ref.update_fraction},
                {"MaximumImbalance", ref.maximum_imbalance},
                {"SaveAdaptIterations", ref.save_adapt_iterations},
                {"SaveAdaptMesh", ref.save_adapt_mesh},
                {"UniformLevels", ref.uniform_ref_levels},
                {"SerialUniformLevels", ref.ser_uniform_ref_levels}});
}

void ConcretizeDomains(const config::DomainData &domains, json &j_domains)
{
  // Materials: positional match to the C++ vector. Emit the scalar (isotropic) baseline
  // when a physical property key is absent, preserving any user-written scalar or
  // tensor form. SymmetricMatrixData<3> with s={v,v,v} and identity v is exactly what
  // the parser reconstructs from a single number.
  if (j_domains.contains("Materials"))
  {
    auto &j_mats = j_domains["Materials"];
    const std::size_t n = std::min(j_mats.size(), domains.materials.size());
    for (std::size_t i = 0; i < n; ++i)
    {
      const auto &m = domains.materials[i];
      ApplyEntries(j_mats[i], {{"Permeability", m.mu_r.s[0]},
                               {"Permittivity", m.epsilon_r.s[0]},
                               {"LossTan", m.tandelta.s[0]},
                               {"Conductivity", m.sigma.s[0]},
                               {"LondonDepth", m.lambda_L}});
    }
  }
}

void ConcretizeBoundaries(const config::BoundaryData &boundaries, json &j_boundaries)
{
  // Absorbing (farfield) boundary. Only touch it if the user declared it.
  if (j_boundaries.contains("Absorbing"))
  {
    Concretize(j_boundaries["Absorbing"], "Order", boundaries.farfield.order);
  }

  // Conductivity: JSON array of objects, C++ vector (positional). Match by index.
  if (j_boundaries.contains("Conductivity"))
  {
    auto &j_cond = j_boundaries["Conductivity"];
    const std::size_t n = std::min(j_cond.size(), boundaries.conductivity.size());
    for (std::size_t i = 0; i < n; ++i)
    {
      const auto &c = boundaries.conductivity[i];
      ApplyEntries(
          j_cond[i],
          {{"Permeability", c.mu_r}, {"Thickness", c.h}, {"External", c.external}});
    }
  }

  // Impedance: JSON array of objects, C++ vector (positional).
  if (j_boundaries.contains("Impedance"))
  {
    auto &j_imp = j_boundaries["Impedance"];
    const std::size_t n = std::min(j_imp.size(), boundaries.impedance.size());
    for (std::size_t i = 0; i < n; ++i)
    {
      const auto &imp = boundaries.impedance[i];
      ApplyEntries(j_imp[i], {{"Rs", imp.Rs}, {"Ls", imp.Ls}, {"Cs", imp.Cs}});
    }
  }

  // LumpedPort: JSON array with Index, C++ map keyed by Index.
  if (j_boundaries.contains("LumpedPort"))
  {
    for (auto &j_port : j_boundaries["LumpedPort"])
    {
      auto idx_it = j_port.find("Index");
      if (idx_it == j_port.end())
      {
        continue;
      }
      auto it = boundaries.lumpedport.find(idx_it->get<int>());
      if (it == boundaries.lumpedport.end())
      {
        continue;
      }
      const auto &lp = it->second;
      ApplyEntries(j_port, {{"R", lp.R},
                            {"L", lp.L},
                            {"C", lp.C},
                            {"Rs", lp.Rs},
                            {"Ls", lp.Ls},
                            {"Cs", lp.Cs},
                            {"Excitation", lp.excitation},
                            {"Active", lp.active}});
    }
  }

  // WavePort: JSON array with Index, C++ map keyed by Index.
  if (j_boundaries.contains("WavePort"))
  {
    for (auto &j_port : j_boundaries["WavePort"])
    {
      auto idx_it = j_port.find("Index");
      if (idx_it == j_port.end())
      {
        continue;
      }
      auto it = boundaries.waveport.find(idx_it->get<int>());
      if (it == boundaries.waveport.end())
      {
        continue;
      }
      const auto &wp = it->second;
      ApplyEntries(j_port, {{"Mode", wp.mode_idx},
                            {"Offset", wp.d_offset},
                            {"SolverType", EnumString(wp.eigen_solver)},
                            {"Excitation", wp.excitation},
                            {"Active", wp.active},
                            {"MaxIts", wp.ksp_max_its},
                            {"KSPTol", wp.ksp_tol},
                            {"EigenTol", wp.eig_tol},
                            {"MaxSize", wp.max_size},
                            {"Verbose", wp.verbose},
                            {"NSamples", wp.n_samples}});
    }
  }

  // Periodic: single object with a Floquet wave vector.
  if (j_boundaries.contains("Periodic"))
  {
    Concretize(
        j_boundaries["Periodic"], "FloquetWaveVector",
        json::array({boundaries.periodic.wave_vector[0], boundaries.periodic.wave_vector[1],
                     boundaries.periodic.wave_vector[2]}));
  }

  // Postprocessing: per-index maps for flux/dielectric/impedance/voltage.
  if (j_boundaries.contains("Postprocessing"))
  {
    auto &j_pp = j_boundaries["Postprocessing"];
    const auto &pp = boundaries.postpro;

    auto ApplyIndexed = [&](const char *key, auto &&per_entry)
    {
      auto sect_it = j_pp.find(key);
      if (sect_it == j_pp.end())
      {
        return;
      }
      for (auto &j_entry : *sect_it)
      {
        auto idx_it = j_entry.find("Index");
        if (idx_it == j_entry.end())
        {
          continue;
        }
        per_entry(idx_it->get<int>(), j_entry);
      }
    };

    ApplyIndexed("SurfaceFlux",
                 [&](int idx, json &j_entry)
                 {
                   auto it = pp.flux.find(idx);
                   if (it != pp.flux.end())
                   {
                     ApplyEntries(j_entry, {{"Type", EnumString(it->second.type)},
                                            {"TwoSided", it->second.two_sided}});
                   }
                 });

    ApplyIndexed("Dielectric",
                 [&](int idx, json &j_entry)
                 {
                   auto it = pp.dielectric.find(idx);
                   if (it != pp.dielectric.end())
                   {
                     ApplyEntries(j_entry, {{"Type", EnumString(it->second.type)},
                                            {"LossTan", it->second.tandelta}});
                   }
                 });

    ApplyIndexed("Impedance",
                 [&](int idx, json &j_entry)
                 {
                   auto it = pp.impedance.find(idx);
                   if (it != pp.impedance.end())
                   {
                     Concretize(j_entry, "NSamples", it->second.n_samples);
                   }
                 });

    ApplyIndexed("Voltage",
                 [&](int idx, json &j_entry)
                 {
                   auto it = pp.voltage.find(idx);
                   if (it != pp.voltage.end())
                   {
                     Concretize(j_entry, "NSamples", it->second.n_samples);
                   }
                 });
  }
}

}  // namespace

json IoData::ConcretizeDefaults(const IoData &iodata, json config)
{
  if (!config.contains("Problem"))
  {
    config["Problem"] = json::object();
  }
  ConcretizeProblem(iodata.problem, config["Problem"]);

  if (!config.contains("Model"))
  {
    config["Model"] = json::object();
  }
  ConcretizeModel(iodata.model, config["Model"]);

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
    case ProblemType::BOUNDARYMODE:
      if (!j_solver.contains("BoundaryMode"))
      {
        j_solver["BoundaryMode"] = json::object();
      }
      ConcretizeBoundaryMode(iodata.solver.boundary_mode, j_solver["BoundaryMode"]);
      break;
  }

  if (config.contains("Domains"))
  {
    ConcretizeDomains(iodata.domains, config["Domains"]);
  }

  if (config.contains("Boundaries"))
  {
    ConcretizeBoundaries(iodata.boundaries, config["Boundaries"]);
  }

  return config;
}

}  // namespace palace
