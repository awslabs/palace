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
#include "utils/constants.hpp"
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

IoData::IoData(const char *filename, bool print) : Lc(1.0), tc(1.0), init(false)
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
  if (problem.type == config::ProblemData::Type::DRIVEN)
  {
    // No unsupported domain or boundary objects for frequency domain driven simulations.
  }
  else if (problem.type == config::ProblemData::Type::EIGENMODE)
  {
    if (!boundaries.conductivity.empty())
    {
      Mpi::Warning("Eigenmode problem type does not support surface conductivity boundary "
                   "conditions!\n");
    }
    if (!boundaries.auxpec.empty() || !boundaries.waveport.empty())
    {
      Mpi::Warning(
          "Eigenmode problem type does not support wave port boundary conditions!\n");
    }
    if (!boundaries.farfield.empty() && boundaries.farfield.order > 1)
    {
      Mpi::Warning("Eigenmode problem type does not support absorbing boundary conditions "
                   "with order > 1!\n");
    }
  }
  else if (problem.type == config::ProblemData::Type::ELECTROSTATIC)
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
  else if (problem.type == config::ProblemData::Type::MAGNETOSTATIC)
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
  else if (problem.type == config::ProblemData::Type::TRANSIENT)
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
  }

  // Resolve default values in configuration file.
  if (solver.linear.type == config::LinearSolverData::Type::DEFAULT)
  {
    if (problem.type == config::ProblemData::Type::ELECTROSTATIC)
    {
      solver.linear.type = config::LinearSolverData::Type::BOOMER_AMG;
    }
    else if (problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
             problem.type == config::ProblemData::Type::TRANSIENT)
    {
      solver.linear.type = config::LinearSolverData::Type::AMS;
    }
    else
    {
      // Prefer sparse direct solver for frequency domain problems if available.
#if defined(MFEM_USE_SUPERLU)
      solver.linear.type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      solver.linear.type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      solver.linear.type = config::LinearSolverData::Type::MUMPS;
#else
      solver.linear.type = config::LinearSolverData::Type::AMS;
#endif
    }
  }
  if (solver.linear.ksp_type == config::LinearSolverData::KspType::DEFAULT)
  {
    // Problems with SPD operators use CG by default, else GMRES.
    if (problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
        problem.type == config::ProblemData::Type::TRANSIENT)
    {
      solver.linear.ksp_type = config::LinearSolverData::KspType::CG;
    }
    else
    {
      solver.linear.ksp_type = config::LinearSolverData::KspType::GMRES;
    }
  }
  if (solver.linear.max_size < 0)
  {
    solver.linear.max_size = solver.linear.max_it;
  }
  if (solver.linear.initial_guess < 0)
  {
    if ((problem.type == config::ProblemData::Type::DRIVEN &&
         solver.driven.adaptive_tol <= 0.0) ||
        problem.type == config::ProblemData::Type::TRANSIENT ||
        problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        problem.type == config::ProblemData::Type::MAGNETOSTATIC)
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
    if (problem.type == config::ProblemData::Type::DRIVEN &&
        solver.linear.type == config::LinearSolverData::Type::AMS)
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
    if (problem.type == config::ProblemData::Type::ELECTROSTATIC ||
        problem.type == config::ProblemData::Type::MAGNETOSTATIC)
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
    solver.linear.ams_singular_op =
        (problem.type == config::ProblemData::Type::MAGNETOSTATIC);
  }
  if (solver.linear.amg_agg_coarsen < 0)
  {
    solver.linear.amg_agg_coarsen =
        (problem.type == config::ProblemData::Type::ELECTROSTATIC ||
         problem.type == config::ProblemData::Type::MAGNETOSTATIC ||
         problem.type == config::ProblemData::Type::TRANSIENT);
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
  Lc = model.Lc * model.L0;                 // [m]
  tc = 1.0e9 * Lc / electromagnetics::c0_;  // [ns]

  // Mesh refinement parameters.
  auto DivideLengthScale = [Lc0 = GetMeshLengthScale()](double val) { return val / Lc0; };
  for (auto &box : model.refinement.GetBoxes())
  {
    std::transform(box.bbmin.begin(), box.bbmin.end(), box.bbmin.begin(),
                   DivideLengthScale);
    std::transform(box.bbmax.begin(), box.bbmax.end(), box.bbmax.begin(),
                   DivideLengthScale);
  }
  for (auto &sphere : model.refinement.GetSpheres())
  {
    sphere.r /= GetMeshLengthScale();
    std::transform(sphere.center.begin(), sphere.center.end(), sphere.center.begin(),
                   DivideLengthScale);
  }

  // Materials: conductivity and London penetration depth.
  for (auto &data : domains.materials)
  {
    data.sigma /= 1.0 / (electromagnetics::Z0_ * Lc);
    data.lambda_L /= GetMeshLengthScale();
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
    data.sigma /= 1.0 / (electromagnetics::Z0_ * Lc);
    data.h /= GetMeshLengthScale();
  }

  // Impedance boundaries and lumped ports.
  for (auto &data : boundaries.impedance)
  {
    data.Rs /= electromagnetics::Z0_;
    data.Ls /= electromagnetics::mu0_ * Lc;
    data.Cs /= electromagnetics::epsilon0_ * Lc;
  }
  for (auto &[idx, data] : boundaries.lumpedport)
  {
    data.R /= electromagnetics::Z0_;
    data.L /= electromagnetics::mu0_ * Lc;
    data.C /= electromagnetics::epsilon0_ * Lc;
    data.Rs /= electromagnetics::Z0_;
    data.Ls /= electromagnetics::mu0_ * Lc;
    data.Cs /= electromagnetics::epsilon0_ * Lc;
  }

  // Wave port offset distance.
  for (auto &[idx, data] : boundaries.waveport)
  {
    data.d_offset /= GetMeshLengthScale();
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
    data.t /= GetMeshLengthScale();
  }

  // For eigenmode simulations:
  solver.eigenmode.target *= 2.0 * M_PI * tc;

  // For driven simulations:
  solver.driven.min_f *= 2.0 * M_PI * tc;
  solver.driven.max_f *= 2.0 * M_PI * tc;
  solver.driven.delta_f *= 2.0 * M_PI * tc;

  // For transient simulations:
  solver.transient.pulse_f *= 2.0 * M_PI * tc;
  solver.transient.pulse_tau /= tc;
  solver.transient.max_t /= tc;
  solver.transient.delta_t /= tc;

  // Scale mesh vertices for correct nondimensionalization.
  mesh::NondimensionalizeMesh(mesh, GetMeshLengthScale());

  // Print some information.
  Mpi::Print(mesh.GetComm(),
             "\nCharacteristic length and time scales:\n L₀ = {:.3e} m, t₀ = {:.3e} ns\n",
             Lc, tc);
}

template <typename T>
T IoData::DimensionalizeValue(IoData::ValueType type, T v) const
{
  // Characteristic reference magnetic field strength Hc² = 1 / (Zc * Lc²) A/m (with Ec =
  // Hc Zc). Yields Pc = Hc² Zc Lc² = 1 W.
  const T Hc = 1.0 / std::sqrt(electromagnetics::Z0_ * Lc * Lc);  // [A/m]
  T sf = 1.0;
  switch (type)
  {
    case ValueType::TIME:
      sf = tc;  // [ns]
      break;
    case ValueType::FREQUENCY:
      sf = 1.0 / (2.0 * M_PI * tc);  // [GHz/rad]
      break;
    case ValueType::LENGTH:
      sf = Lc;  // [m]
      break;
    case ValueType::IMPEDANCE:
      sf = electromagnetics::Z0_;  // [Ω]
      break;
    case ValueType::INDUCTANCE:
      sf = electromagnetics::mu0_ * Lc;  // [H]
      break;
    case ValueType::CAPACITANCE:
      sf = electromagnetics::epsilon0_ * Lc;  // [F]
      break;
    case ValueType::CONDUCTIVITY:
      sf = 1.0 / (electromagnetics::Z0_ * Lc);  // [S/m]
      break;
    case ValueType::VOLTAGE:
      sf = Hc * electromagnetics::Z0_ * Lc;  // [V]
      break;
    case ValueType::CURRENT:
      sf = Hc * Lc;  // [A]
      break;
    case ValueType::POWER:
      sf = Hc * Hc * electromagnetics::Z0_ * Lc * Lc;  // [W]
      break;
    case ValueType::ENERGY:
      sf = Hc * Hc * electromagnetics::Z0_ * Lc * Lc * tc;  // [J]
      break;
    case ValueType::FIELD_E:
      sf = Hc * electromagnetics::Z0_;  // [V/m]
      break;
    case ValueType::FIELD_D:
      sf = electromagnetics::epsilon0_ * Hc * electromagnetics::Z0_;  // [C/m²]
      break;
    case ValueType::FIELD_H:
      sf = Hc;  // [A/m]
      break;
    case ValueType::FIELD_B:
      sf = electromagnetics::mu0_ * Hc;  // [Wb/m²]
      break;
  }
  return v * sf;
}

template double IoData::DimensionalizeValue(IoData::ValueType, double) const;

}  // namespace palace
