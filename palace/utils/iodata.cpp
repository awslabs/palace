// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "iodata.hpp"

#include <fstream>
#include <functional>
#include <iostream>
#include <regex>
#include <stack>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "utils/communication.hpp"
#include "utils/constants.hpp"

namespace palace
{

using json = nlohmann::json;

IoData::IoData(const char *filename, bool print) : Lc(1.0), tc(1.0), init(false)
{
  // Open configuration file and preprocess: strip whitespace, comments, and expand integer
  // ranges.
  std::stringstream buffer;
  PreprocessFile(filename, buffer);

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

void IoData::PreprocessFile(const char *filename, std::stringstream &buffer)
{
  // Read configuration file into memory.
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
  // within strings and escaped comment markers (see tinyurl.com/2s3n8dkr).
  {
    std::regex rgx(R"((([\"'])(?:(?=(\\?))\3.)*?\2))"
                   R"(|(\/\*([^*]|[\r\n]|(\*+([^*\/]|[\r\n])))*\*+\/))"
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
  auto RangeExpand = [](const std::string &str, std::size_t pos) -> std::string
  {
    // Handle the substring str.substr(0, pos), which is only numeric with possible hyphens.
    MFEM_VERIFY(pos != std::string::npos, "Invalid string size in range expansion!");
    std::string rng;
    std::size_t size, size2;
    int number = std::stoi(str, &size);
    MFEM_VERIFY(size <= pos, "Unexpected stoi result in range expansion!");
    if (size < pos)
    {
      // Range specified, expand the bounds.
      MFEM_VERIFY(str[size] == '-', "Invalid character encountered in range expansion!");
      int number2 = std::stoi(str.substr(size + 1), &size2);
      MFEM_VERIFY(size + size2 + 1 == pos, "Unexpected stoi result in range expansion!");
      MFEM_VERIFY(number < number2, "Invalid range bounds in range expansion!");
      while (number < number2)
      {
        rng += std::to_string(number++) + ",";
      }
      rng += std::to_string(number);
    }
    else
    {
      // Just push back the number.
      rng += str.substr(0, size);
    }
    if (pos != str.length())
    {
      rng += ",";
    }
    return rng;
  };
  {
    buffer.str(std::string(""));  // Clear the output buffer
    std::regex rgx(R"(\[(-?[0-9][\-\,0-9]*[0-9])\])");
    auto it = file.cbegin();
    const auto end = file.cend();
    for (std::smatch match; std::regex_search(it, end, match, rgx); it = match[0].second)
    {
      // Apply integer range expansion (as needed) to the first capture group. The match
      // includes only digits, commas, and '-'.
      std::string str = match[1].str(), range;
      std::size_t pos = 0;
      MFEM_VERIFY(str.find_first_not_of(",-0123456789") == std::string::npos,
                  "Range expansion expects only integer values!");
      buffer << match.prefix() << '[';
      while ((pos = str.find(',')) != std::string::npos)
      {
        buffer << RangeExpand(str, pos);
        str.erase(0, pos + 1);
      }
      buffer << RangeExpand(str, str.length()) << ']';
    }
    buffer << std::string(it, end);
  }
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
    if (!boundaries.postpro.inductance.empty())
    {
      Mpi::Warning("Electrostatic problem type does not support boundary inductance "
                   "postprocessing!\n");
    }
  }
  else if (problem.type == config::ProblemData::Type::MAGNETOSTATIC)
  {
    if (!domains.postpro.dielectric.empty())
    {
      Mpi::Warning("Magnetostatic problem type does not support domain bulk dielectric "
                   "loss postprocessing!\n");
    }
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
    if (!boundaries.postpro.capacitance.empty())
    {
      Mpi::Warning("Magnetostatic problem type does not support boundary capacitance "
                   "postprocessing!\n");
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
      // simulations, or electrostatic or magnetostatics.
      solver.linear.initial_guess = 1;
    }
    else
    {
      solver.linear.initial_guess = 0;
    }
  }
  if (solver.linear.pc_mat_shifted < 0)
  {
    solver.linear.pc_mat_shifted = 0;  // Default false for most cases
    if (problem.type == config::ProblemData::Type::DRIVEN)
    {
#if defined(MFEM_USE_SUPERLU) || defined(MFEM_USE_STRUMPACK) || defined(MFEM_USE_MUMPS)
      if (solver.linear.type == config::LinearSolverData::Type::AMS)
#else
      if (solver.linear.type == config::LinearSolverData::Type::AMS ||
          solver.linear.type == config::LinearSolverData::Type::DEFAULT)
#endif
      {
        // Default true only driven simulations using AMS.
        solver.linear.pc_mat_shifted = 1;
      }
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
  // Nondimensionalization of the equations is based on a given length Lc in[m], typically
  // the largest domain dimension. Configuration file lengths and the mesh coordinates are
  // provided with units of model.L0 x [m].
  MFEM_VERIFY(!init, "NondimensionalizeInputs should only be called once!");
  init = true;

  // Calculate the reference length and time.
  if (model.Lc > 0.0)
  {
    // User specified Lc in mesh length units.
    Lc = model.Lc * model.L0;  // [m]
  }
  else
  {
    mfem::Vector bbmin, bbmax;
    mesh.GetBoundingBox(bbmin, bbmax);
    bbmin *= model.L0;  // [m]
    bbmax *= model.L0;
    Lc = -mfem::infinity();
    for (int d = 0; d < mesh.SpaceDimension(); d++)
    {
      double l = bbmax(d) - bbmin(d);
      if (Lc < l)
      {
        Lc = l;  // [m]
      }
    }
  }
  tc = 1.0e9 * Lc / electromagnetics::c0_;  // [ns]

  // Mesh refinement parameters.
  auto Divides = [this](double val) { return val / (Lc / model.L0); };
  for (auto &box : model.refinement.GetBoxes())
  {
    std::transform(box.bbmin.begin(), box.bbmin.end(), box.bbmin.begin(), Divides);
    std::transform(box.bbmax.begin(), box.bbmax.end(), box.bbmax.begin(), Divides);
  }
  for (auto &sphere : model.refinement.GetSpheres())
  {
    sphere.r /= Lc / model.L0;
    std::transform(sphere.center.begin(), sphere.center.end(), sphere.center.begin(),
                   Divides);
  }

  // Materials: conductivity and London penetration depth.
  for (auto &data : domains.materials)
  {
    data.sigma /= 1.0 / (electromagnetics::Z0_ * Lc);
    data.lambda_L /= Lc / model.L0;
  }

  // Probe location coordinates.
  for (auto &[idx, data] : domains.postpro.probe)
  {
    data.x /= Lc / model.L0;
    data.y /= Lc / model.L0;
    data.z /= Lc / model.L0;
  }

  // Finite conductivity boundaries.
  for (auto &data : boundaries.conductivity)
  {
    data.sigma /= 1.0 / (electromagnetics::Z0_ * Lc);
    data.h /= Lc / model.L0;
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
    data.d_offset /= Lc / model.L0;
  }

  // Dielectric interface thickness.
  for (auto &[idx, data] : boundaries.postpro.dielectric)
  {
    data.ts /= Lc / model.L0;
  }

  // For eigenmode simulations:
  solver.eigenmode.target *= 2.0 * M_PI * tc;
  solver.eigenmode.feast_contour_ub *= 2.0 * M_PI * tc;

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
  for (int i = 0; i < mesh.GetNV(); i++)
  {
    double *v = mesh.GetVertex(i);
    std::transform(v, v + mesh.SpaceDimension(), v, Divides);
  }
  if (mesh.GetNodes())
  {
    *mesh.GetNodes() /= Lc / model.L0;
  }

  // Print some information.
  Mpi::Print(mesh.GetComm(),
             "\nCharacteristic length and time scales:\n L₀ = {:.3e} m, t₀ = {:.3e} ns\n",
             Lc, tc);
}

double IoData::DimensionalizeValue(IoData::ValueType type, double v) const
{
  // XX TODO: Add more for fields, currents, voltages, energies
  double sf;
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
    default:
      MFEM_ABORT("Unsupported value type for dimensionalization!");
      sf = 1.0;  // For compiler warning
      break;
  }
  return v * sf;
}

}  // namespace palace
