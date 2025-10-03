// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iodata.hpp"

#include <charconv>
#include <fstream>
#include <functional>
#include <iostream>
#include <regex>
#include <sstream>
#include <stack>
#include <string>
#include <string_view>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

namespace palace
{

std::stringstream PreprocessFile(const char *filename)
{
  // Read configuration file into memory and return as a stringstream.
  std::string file;
  {
    std::ifstream fi(filename);
    std::stringstream buf;
    if (!fi.is_open())
    {
      MFEM_ABORT("Unable to open configuration file \"" << filename << "\"!");
    }
    buf << fi.rdbuf();
    fi.close();
    file = buf.str();
  }

  // Strip C and C++ style comments (//, /* */) using regex. Correctly handles comments
  // within strings and escaped comment markers (see tinyurl.com/2s3n8dkr). An alternative
  // for the middle line is: R"(|(\/\*([^*]|[\r\n]|(\*+([^*\/]|[\r\n])))*\*+\/))", but this
  // seems to sometimes lead to issues with std::regex_replace for long files.
  {
    std::regex rgx(R"((([\"'])(?:(?=(\\?))\3.)*?\2))"
                   R"(|(\/\*(.|[\r\n])*?\*\/))"
                   R"(|(\/\/.*))");
    file = std::regex_replace(file, rgx, "$1");
  }

  // Also strip whitespace.
  {
    std::regex rgx(R"((([\"'])(?:(?=(\\?))\3.)*?\2))"
                   R"(|(\s+))");
    file = std::regex_replace(file, rgx, "$1");
  }

  // Also strip erroneous trailing commas.
  {
    std::regex rgx(R"((([\"'])(?:(?=(\\?))\3.)*?\2))"
                   R"(|,+(?=\s*?[\}\]]))");
    file = std::regex_replace(file, rgx, "$1");
  }

  // Perform integer range expansion for arrays ([a - b, c] = [a-b,c] =
  // [a,a+1,...,b-1,b,c]). The whole file is now one line and arrays have no spaces after
  // whitespace stripping.
  std::stringstream output;
  auto RangeExpand = [](std::string_view str) -> std::string
  {
    // Handle the given string which is only numeric with possible hyphens.
    if (str.empty())
    {
      return "";
    }
    int num;
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.length(), num);
    MFEM_VERIFY(
        ec == std::errc(),
        "Invalid integer conversion in range expansion"
            << (ec == std::errc::result_out_of_range ? " (integer out of range)!" : "!"));
    if (ptr == str.data() + str.length())
    {
      return std::string(str);
    }
    // Range specified, expand the bounds.
    int num2;
    auto [ptr2, ec2] = std::from_chars(ptr + 1, str.data() + str.length(), num2);
    MFEM_VERIFY(
        ec2 == std::errc(),
        "Invalid integer conversion in range expansion"
            << (ec2 == std::errc::result_out_of_range ? " (integer out of range)!" : "!"));
    std::string rng;
    while (num < num2)
    {
      rng += std::to_string(num++) + ",";
    }
    rng += std::to_string(num);
    return rng;
  };
  {
    const std::string range_vals = "-0123456789,";
    auto start = file.begin();
    bool inside = false;
    for (auto it = start; it != file.end(); ++it)
    {
      if (inside)
      {
        if (*it == ']')
        {
          // Apply integer range expansion (as needed) to the array, which includes only
          // digits, commas, and '-'. Exclude the outer square brackets.
          std::string_view str(file.data() + (start - file.cbegin() + 1), it - start - 1);
          std::size_t s = 0, pos;
          output << '[';
          while ((pos = str.find(',', s)) != std::string::npos)
          {
            output << RangeExpand(str.substr(s, pos - s)) << ',';
            s = pos + 1;
          }
          output << RangeExpand(str.substr(s)) << ']';
          start = it + 1;
          inside = false;
        }
        else if (*it == '[')
        {
          output << std::string(start, it);
          start = it;
        }
        else if (range_vals.find(*it) == std::string::npos)
        {
          output << std::string(start, it);
          start = it;
          inside = false;
        }
      }
      else if (*it == '[')
      {
        output << std::string(start, it);
        start = it;
        inside = true;
      }
    }
    output << std::string(start, file.end());
  }
  return output;
}

using json = nlohmann::json;

IoData::IoData(const char *filename, bool print) : units(1.0, 1.0), init(false)
{
  // Open configuration file and preprocess: strip whitespace, comments, and expand integer
  // ranges.
  std::stringstream buffer = PreprocessFile(filename);

  // Parse the configuration file. Use a callback function to detect and throw errors for
  // duplicate keys.
  json config;
  std::stack<std::set<json>> parse_stack;
  json::parser_callback_t check_duplicate_keys =
      [&](int, json::parse_event_t event, json &parsed)
  {
    switch (event)
    {
      case json::parse_event_t::object_start:
        parse_stack.push(std::set<json>());
        break;
      case json::parse_event_t::object_end:
        parse_stack.pop();
        break;
      case json::parse_event_t::key:
        {
          const auto result = parse_stack.top().insert(parsed);
          if (!result.second)
          {
            MFEM_ABORT("Error parsing configuration file!\nDuplicate key "
                       << parsed << " was already seen in this object!");
            return false;
          }
        }
        break;
      default:
        break;
    }
    return true;
  };
  try
  {
    config = json::parse(buffer, check_duplicate_keys);
  }
  catch (json::parse_error &e)
  {
    MFEM_ABORT("Error parsing configuration file!\n  " << e.what());
  }
  if (print)
  {
    Mpi::Print("\n{}\n", config.dump(2));
  }

  // Set up configuration option data structures.
  problem.SetUp(config);
  model.SetUp(config);
  domains.SetUp(config);
  boundaries.SetUp(config);
  solver.SetUp(config);

  // Cleanup and error checking.
  config.erase("Problem");
  config.erase("Model");
  config.erase("Domains");
  config.erase("Boundaries");
  config.erase("Solver");
  MFEM_VERIFY(config.empty(), "Found an unsupported configuration file section!\n"
                                  << config.dump(2));

  // Check compatibility of configuration file and problem type.
  CheckConfiguration();
}

void IoData::CheckConfiguration()
{
  // Check that the provided domain and boundary objects are all supported by the requested
  // problem type.
  if (problem.type == ProblemType::DRIVEN)
  {
    // No unsupported domain or boundary objects for frequency domain driven simulations.
  }
  else if (problem.type == ProblemType::EIGENMODE)
  {
    // No unsupported domain or boundary objects for frequency domain driven simulations.
  }
  else if (problem.type == ProblemType::ELECTROSTATIC)
  {
    if (!boundaries.farfield.empty())
    {
      Mpi::Warning(
          "Electrostatic problem type does not support absorbing boundary conditions!\n");
    }
    if (!boundaries.conductivity.empty())
    {
      Mpi::Warning("Electrostatic problem type does not support surface conductivity "
                   "boundary conditions!\n");
    }
    if (!boundaries.impedance.empty())
    {
      Mpi::Warning("Electrostatic problem type does not support surface impedance boundary "
                   "conditions!\n");
    }
    if (!boundaries.auxpec.empty() || !boundaries.waveport.empty())
    {
      Mpi::Warning(
          "Electrostatic problem type does not support wave port boundary conditions!\n");
    }
    if (!boundaries.current.empty())
    {
      Mpi::Warning(
          "Electrostatic problem type does not support surface current excitation!\n");
    }
  }
  else if (problem.type == ProblemType::MAGNETOSTATIC)
  {
    if (!boundaries.farfield.empty())
    {
      Mpi::Warning(
          "Magnetostatic problem type does not support absorbing boundary conditions!\n");
    }
    if (!boundaries.conductivity.empty())
    {
      Mpi::Warning("Magnetostatic problem type does not support surface conductivity "
                   "boundary conditions!\n");
    }
    if (!boundaries.impedance.empty())
    {
      Mpi::Warning("Magnetostatic problem type does not support surface impedance boundary "
                   "conditions!\n");
    }
    if (!boundaries.lumpedport.empty())
    {
      Mpi::Warning(
          "Magnetostatic problem type does not support lumped port boundary conditions!\n");
    }
    if (!boundaries.auxpec.empty() || !boundaries.waveport.empty())
    {
      Mpi::Warning(
          "Magnetostatic problem type does not support wave port boundary conditions!\n");
    }
    if (!boundaries.postpro.dielectric.empty())
    {
      Mpi::Warning("Magnetostatic problem type does not support boundary interface "
                   "dielectric loss postprocessing!\n");
    }
  }
  else if (problem.type == ProblemType::TRANSIENT)
  {
    if (!boundaries.conductivity.empty())
    {
      Mpi::Warning("Transient problem type does not support surface conductivity boundary "
                   "conditions!\n");
    }
    if (!boundaries.auxpec.empty() || !boundaries.waveport.empty())
    {
      Mpi::Warning(
          "Transient problem type does not support wave port boundary conditions!\n");
    }
    if (!boundaries.farfield.empty() && boundaries.farfield.order > 1)
    {
      Mpi::Warning("Transient problem type does not support absorbing boundary conditions "
                   "with order > 1!\n");
    }
  }

  // Resolve default values in configuration file.
  if (solver.linear.type == LinearSolver::DEFAULT)
  {
    if (problem.type == ProblemType::ELECTROSTATIC)
    {
      solver.linear.type = LinearSolver::BOOMER_AMG;
    }
    else if (problem.type == ProblemType::MAGNETOSTATIC ||
             problem.type == ProblemType::TRANSIENT)
    {
      solver.linear.type = LinearSolver::AMS;
    }
    else
    {
      // Prefer sparse direct solver for frequency domain problems if available.
#if defined(MFEM_USE_SUPERLU)
      solver.linear.type = LinearSolver::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      solver.linear.type = LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      solver.linear.type = LinearSolver::MUMPS;
#else
      solver.linear.type = LinearSolver::AMS;
#endif
    }
  }
  if (solver.linear.krylov_solver == KrylovSolver::DEFAULT)
  {
    // Problems with SPD operators use CG by default, else GMRES.
    if (problem.type == ProblemType::ELECTROSTATIC ||
        problem.type == ProblemType::MAGNETOSTATIC ||
        problem.type == ProblemType::TRANSIENT)
    {
      solver.linear.krylov_solver = KrylovSolver::CG;
    }
    else
    {
      solver.linear.krylov_solver = KrylovSolver::GMRES;
    }
  }
  if (solver.linear.max_size < 0)
  {
    solver.linear.max_size = solver.linear.max_it;
  }
  if (solver.linear.initial_guess < 0)
  {
    if ((problem.type == ProblemType::DRIVEN && solver.driven.adaptive_tol <= 0.0) ||
        problem.type == ProblemType::TRANSIENT ||
        problem.type == ProblemType::ELECTROSTATIC ||
        problem.type == ProblemType::MAGNETOSTATIC)
    {
      // Default true only driven simulations without adaptive frequency sweep, transient
      // simulations, electrostatics, or magnetostatics.
      solver.linear.initial_guess = 1;
    }
    else
    {
      solver.linear.initial_guess = 0;
    }
  }
  if (solver.linear.pc_mat_shifted < 0)
  {
    if (problem.type == ProblemType::DRIVEN && solver.linear.type == LinearSolver::AMS)
    {
      // Default true only driven simulations using AMS (false for most cases).
      solver.linear.pc_mat_shifted = 1;
    }
    else
    {
      solver.linear.pc_mat_shifted = 0;
    }
  }
  if (solver.linear.mg_smooth_aux < 0)
  {
    if (problem.type == ProblemType::ELECTROSTATIC ||
        problem.type == ProblemType::MAGNETOSTATIC)
    {
      // Disable auxiliary space smoothing using distributive relaxation by default for
      // problems which don't need it.
      solver.linear.mg_smooth_aux = 0;
    }
    else
    {
      solver.linear.mg_smooth_aux = 1;
    }
  }
  if (solver.linear.mg_smooth_order < 0)
  {
    solver.linear.mg_smooth_order = std::max(2 * solver.order, 4);
  }
  if (solver.linear.ams_singular_op < 0)
  {
    solver.linear.ams_singular_op = (problem.type == ProblemType::MAGNETOSTATIC);
  }
  if (solver.linear.amg_agg_coarsen < 0)
  {
    solver.linear.amg_agg_coarsen = (problem.type == ProblemType::ELECTROSTATIC ||
                                     problem.type == ProblemType::MAGNETOSTATIC ||
                                     problem.type == ProblemType::TRANSIENT);
  }
  if (solver.linear.reorder_reuse && solver.linear.drop_small_entries &&
      solver.linear.complex_coarse_solve && (problem.type == ProblemType::EIGENMODE) &&
      (!boundaries.waveport.empty() || !boundaries.conductivity.empty() ||
       (!boundaries.farfield.empty() && boundaries.farfield.order > 1)))
  {
    // Do not reuse the sparsity pattern for nonlinear eigenmode simulations with complex
    // coarse preconditioners when dropping small entries. In those cases, the sparsity
    // pattern of the first preconditioner (purely real coefficients) will be different from
    // subsequent preconditioners with complex coefficients.
    solver.linear.reorder_reuse = false;
  }
  // Configure settings for quadrature rules and partial assembly.
  BilinearForm::pa_order_threshold = solver.pa_order_threshold;
  fem::DefaultIntegrationOrder::p_trial = solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = solver.q_order_extra;
}

namespace
{

template <std::size_t N>
constexpr config::SymmetricMatrixData<N> &operator/=(config::SymmetricMatrixData<N> &data,
                                                     double s)
{
  for (auto &x : data.s)
  {
    x /= s;
  }
  return data;
}

}  // namespace

void IoData::NondimensionalizeInputs(mfem::ParMesh &mesh)
{
  // Nondimensionalization of the equations is based on a given length Lc in [m], typically
  // the largest domain dimension. Configuration file lengths and the mesh coordinates are
  // provided with units of model.L0 x [m].
  MFEM_VERIFY(!init, "NondimensionalizeInputs should only be called once!");
  init = true;

  // Calculate the reference length and time. A user specified model.Lc is in mesh length
  // units.
  if (model.Lc <= 0.0)
  {
    mfem::Vector bbmin, bbmax;
    mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
    bbmax -= bbmin;
    model.Lc = *std::max_element(bbmax.begin(), bbmax.end());
  }
  // Define units now mesh length set. Note: In model field Lc is measured in units of L0.
  units = Units(model.L0, model.Lc * model.L0);

  // Mesh refinement parameters.
  auto DivideLengthScale = [Lc0 = units.GetMeshLengthRelativeScale()](double val)
  { return val / Lc0; };
  for (auto &box : model.refinement.GetBoxes())
  {
    std::transform(box.bbmin.begin(), box.bbmin.end(), box.bbmin.begin(),
                   DivideLengthScale);
    std::transform(box.bbmax.begin(), box.bbmax.end(), box.bbmax.begin(),
                   DivideLengthScale);
  }
  for (auto &sphere : model.refinement.GetSpheres())
  {
    sphere.r /= units.GetMeshLengthRelativeScale();
    std::transform(sphere.center.begin(), sphere.center.end(), sphere.center.begin(),
                   DivideLengthScale);
  }

  // Materials: conductivity and London penetration depth.
  for (auto &data : domains.materials)
  {
    data.sigma /= units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>();
    data.lambda_L /= units.GetMeshLengthRelativeScale();
  }

  // Probe location coordinates.
  for (auto &[idx, data] : domains.postpro.probe)
  {
    std::transform(data.center.begin(), data.center.end(), data.center.begin(),
                   DivideLengthScale);
  }

  // Finite conductivity boundaries.
  for (auto &data : boundaries.conductivity)
  {
    data.sigma /= units.GetScaleFactor<Units::ValueType::CONDUCTIVITY>();
    data.h /= units.GetMeshLengthRelativeScale();
  }

  // Impedance boundaries and lumped ports.
  for (auto &data : boundaries.impedance)
  {
    data.Rs /= units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    data.Ls /= units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
    data.Cs /= units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  }
  for (auto &[idx, data] : boundaries.lumpedport)
  {
    data.R /= units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    data.L /= units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
    data.C /= units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
    data.Rs /= units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
    data.Ls /= units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
    data.Cs /= units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  }

  // Floquet periodic boundaries.
  for (auto &k : boundaries.periodic.wave_vector)
  {
    k *= units.GetMeshLengthRelativeScale();
  }

  // Wave port offset distance.
  for (auto &[idx, data] : boundaries.waveport)
  {
    data.d_offset /= units.GetMeshLengthRelativeScale();
  }

  // Center coordinates for surface flux.
  for (auto &[idx, data] : boundaries.postpro.flux)
  {
    std::transform(data.center.begin(), data.center.end(), data.center.begin(),
                   DivideLengthScale);
  }

  // Dielectric interface thickness.
  for (auto &[idx, data] : boundaries.postpro.dielectric)
  {
    data.t /= units.GetMeshLengthRelativeScale();
  }

  // For eigenmode simulations:
  solver.eigenmode.target /= units.GetScaleFactor<Units::ValueType::FREQUENCY>();
  solver.eigenmode.target_upper /= units.GetScaleFactor<Units::ValueType::FREQUENCY>();

  // For driven simulations:
  for (auto &f : solver.driven.sample_f)
    f /= units.GetScaleFactor<Units::ValueType::FREQUENCY>();

  // For transient simulations:
  solver.transient.pulse_f /= units.GetScaleFactor<Units::ValueType::FREQUENCY>();
  solver.transient.pulse_tau /= units.GetScaleFactor<Units::ValueType::TIME>();
  solver.transient.max_t /= units.GetScaleFactor<Units::ValueType::TIME>();
  solver.transient.delta_t /= units.GetScaleFactor<Units::ValueType::TIME>();

  // Scale mesh vertices for correct nondimensionalization.
  mesh::NondimensionalizeMesh(mesh, units.GetMeshLengthRelativeScale());

  // Print some information.
  Mpi::Print(mesh.GetComm(),
             "\nCharacteristic length and time scales:\n L₀ = {:.3e} m, t₀ = {:.3e} ns\n",
             units.GetScaleFactor<Units::ValueType::LENGTH>(),
             units.GetScaleFactor<Units::ValueType::TIME>());
}

}  // namespace palace
