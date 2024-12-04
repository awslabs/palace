// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "configfile.hpp"

#include <algorithm>
#include <iterator>
#include <string_view>
#include <fmt/format.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "models/portexcitationhelper.hpp"

// This is similar to NLOHMANN_JSON_SERIALIZE_ENUM, but results in an error if an enum
// value corresponding to the string cannot be found. Also adds an overload for stream
// printing enum values.
#define PALACE_JSON_SERIALIZE_ENUM(ENUM_TYPE, ...)                                  \
  template <typename BasicJsonType>                                                 \
  inline void to_json(BasicJsonType &j, const ENUM_TYPE &e)                         \
  {                                                                                 \
    static_assert(std::is_enum<ENUM_TYPE>::value, #ENUM_TYPE " must be an enum!");  \
    static const std::pair<ENUM_TYPE, BasicJsonType> m[] = __VA_ARGS__;             \
    auto it = std::find_if(std::begin(m), std::end(m),                              \
                           [e](const std::pair<ENUM_TYPE, BasicJsonType> &ej_pair)  \
                           { return ej_pair.first == e; });                         \
    MFEM_VERIFY(it != std::end(m),                                                  \
                "Invalid value for " << #ENUM_TYPE " given when parsing to JSON!"); \
    j = it->second;                                                                 \
  }                                                                                 \
  template <typename BasicJsonType>                                                 \
  inline void from_json(const BasicJsonType &j, ENUM_TYPE &e)                       \
  {                                                                                 \
    static_assert(std::is_enum<ENUM_TYPE>::value, #ENUM_TYPE " must be an enum!");  \
    static const std::pair<ENUM_TYPE, BasicJsonType> m[] = __VA_ARGS__;             \
    auto it = std::find_if(std::begin(m), std::end(m),                              \
                           [j](const std::pair<ENUM_TYPE, BasicJsonType> &ej_pair)  \
                           { return ej_pair.second == j; });                        \
    MFEM_VERIFY(it != std::end(m),                                                  \
                "Invalid value (" << j << ") for "                                  \
                                  << #ENUM_TYPE                                     \
                    " given in the configuration file when parsing from JSON!");    \
    e = it->first;                                                                  \
  }                                                                                 \
  std::ostream &operator<<(std::ostream &os, const ENUM_TYPE &e)                    \
  {                                                                                 \
    static const std::pair<ENUM_TYPE, const char *> m[] = __VA_ARGS__;              \
    os << std::find_if(std::begin(m), std::end(m),                                  \
                       [e](const std::pair<ENUM_TYPE, const char *> &ej_pair)       \
                       { return ej_pair.first == e; })                              \
              ->second;                                                             \
    return os;                                                                      \
  }

namespace palace::config
{

using json = nlohmann::json;

namespace internal
{

// Helper for converting string keys to enum for ElementData::CoordinateSystem.
PALACE_JSON_SERIALIZE_ENUM(ElementData::CoordinateSystem,
                           {{ElementData::CoordinateSystem::CARTESIAN, "Cartesian"},
                            {ElementData::CoordinateSystem::CYLINDRICAL, "Cylindrical"}})

}  // namespace internal

namespace
{

int AtIndex(json::iterator &port_it, std::string_view errmsg_parent)
{
  MFEM_VERIFY(
      port_it->find("Index") != port_it->end(),
      fmt::format("Missing {} \"Index\" in the configuration file!", errmsg_parent));
  auto index = port_it->at("Index").get<int>();
  MFEM_VERIFY(index > 0, fmt::format("The {} \"Index\" should be an integer > 0; got {}",
                                     errmsg_parent, index));
  return index;
}

template <std::size_t N>
void ParseSymmetricMatrixData(json &mat, const std::string &name,
                              SymmetricMatrixData<N> &data)
{
  auto it = mat.find(name);
  if (it != mat.end() && it->is_array())
  {
    // Attempt to parse as an array.
    data.s = it->get<std::array<double, N>>();
  }
  else
  {
    // Fall back to scalar parsing with default.
    double s = mat.value(name, data.s[0]);
    data.s.fill(s);
  }
  data.v = mat.value("MaterialAxes", data.v);
}

// Helper function for extracting element data from the configuration file, either from a
// provided keyword argument of from a specified vector. In extracting the direction various
// checks are performed for validity of the input combinations.
void ParseElementData(json &elem, const std::string &name, bool required,
                      internal::ElementData &data)
{
  data.attributes = elem.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(data.attributes.begin(), data.attributes.end());
  auto it = elem.find(name);
  if (it != elem.end() && it->is_array())
  {
    // Attempt to parse as an array.
    data.direction = it->get<std::array<double, 3>>();
    data.coordinate_system = elem.value("CoordinateSystem", data.coordinate_system);
  }
  else
  {
    // Fall back to parsing as a string (value is optional).
    MFEM_VERIFY(elem.find("CoordinateSystem") == elem.end(),
                "Cannot specify \"CoordinateSystem\" when specifying a direction or side "
                "using a string in the configuration file!");
    std::string direction;
    direction = elem.value(name, direction);
    for (auto &c : direction)
    {
      c = std::tolower(c);
    }
    const auto xpos = direction.find("x");
    const auto ypos = direction.find("y");
    const auto zpos = direction.find("z");
    const auto rpos = direction.find("r");
    const bool xfound = xpos != std::string::npos;
    const bool yfound = ypos != std::string::npos;
    const bool zfound = zpos != std::string::npos;
    const bool rfound = rpos != std::string::npos;
    if (xfound)
    {
      MFEM_VERIFY(direction.length() == 1 || direction[xpos - 1] == '-' ||
                      direction[xpos - 1] == '+',
                  "Missing required sign specification on \"X\" for \""
                      << name << "\" in the configuration file!");
      MFEM_VERIFY(!yfound && !zfound && !rfound,
                  "\"X\" cannot be combined with \"Y\", \"Z\", or \"R\" for \""
                      << name << "\" in the configuration file!");
      data.direction[0] =
          (direction.length() == 1 || direction[xpos - 1] == '+') ? 1.0 : -1.0;
      data.coordinate_system = internal::ElementData::CoordinateSystem::CARTESIAN;
    }
    if (yfound)
    {
      MFEM_VERIFY(direction.length() == 1 || direction[ypos - 1] == '-' ||
                      direction[ypos - 1] == '+',
                  "Missing required sign specification on \"Y\" for \""
                      << name << "\" in the configuration file!");
      MFEM_VERIFY(!xfound && !zfound && !rfound,
                  "\"Y\" cannot be combined with \"X\", \"Z\", or \"R\" for \""
                      << name << "\" in the configuration file!");
      data.direction[1] =
          direction.length() == 1 || direction[ypos - 1] == '+' ? 1.0 : -1.0;
      data.coordinate_system = internal::ElementData::CoordinateSystem::CARTESIAN;
    }
    if (zfound)
    {
      MFEM_VERIFY(direction.length() == 1 || direction[zpos - 1] == '-' ||
                      direction[zpos - 1] == '+',
                  "Missing required sign specification on \"Z\" for \""
                      << name << "\" in the configuration file!");
      MFEM_VERIFY(!xfound && !yfound && !rfound,
                  "\"Z\" cannot be combined with \"X\", \"Y\", or \"R\" for \""
                      << name << "\" in the configuration file!");
      data.direction[2] =
          direction.length() == 1 || direction[zpos - 1] == '+' ? 1.0 : -1.0;
      data.coordinate_system = internal::ElementData::CoordinateSystem::CARTESIAN;
    }
    if (rfound)
    {
      MFEM_VERIFY(direction.length() == 1 || direction[rpos - 1] == '-' ||
                      direction[rpos - 1] == '+',
                  "Missing required sign specification on \"R\" for \""
                      << name << "\" in the configuration file!");
      MFEM_VERIFY(!xfound && !yfound && !zfound,
                  "\"R\" cannot be combined with \"X\", \"Y\", or \"Z\" for \""
                      << name << "\" in the configuration file!");
      data.direction[0] =
          direction.length() == 1 || direction[rpos - 1] == '+' ? 1.0 : -1.0;
      data.direction[1] = 0.0;
      data.direction[2] = 0.0;
      data.coordinate_system = internal::ElementData::CoordinateSystem::CYLINDRICAL;
    }
  }
  MFEM_VERIFY(data.coordinate_system !=
                      internal::ElementData::CoordinateSystem::CYLINDRICAL ||
                  (data.direction[1] == 0.0 && data.direction[2] == 0.0),
              "Parsing azimuthal and longitudinal directions for cylindrical coordinate "
              "system directions from the configuration file is not currently supported!");
  MFEM_VERIFY(
      !required || data.direction[0] != 0.0 || data.direction[1] != 0.0 ||
          data.direction[2] != 0.0,
      "Missing \"" << name
                   << "\" for an object which requires it in the configuration file!");
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &data)
{
  bool first = true;
  for (const auto &x : data)
  {
    if (!first)
    {
      os << ' ';
    }
    os << x;
    first = false;
  }
  return os;
}

template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &data)
{
  bool first = true;
  for (const auto &x : data)
  {
    if (!first)
    {
      os << ' ';
    }
    os << x;
    first = false;
  }
  return os;
}

template <std::size_t N>
std::ostream &operator<<(std::ostream &os, const SymmetricMatrixData<N> &data)
{
  os << "s: " << data.s;
  int j = 0;
  for (const auto &x : data.v)
  {
    os << ", v" << j++ << ": " << x;
  }
  return os;
}

constexpr bool JSON_DEBUG = false;

}  // namespace

// Helper for converting string keys to enum for ProblemData::Type.
PALACE_JSON_SERIALIZE_ENUM(ProblemData::Type,
                           {{ProblemData::Type::DRIVEN, "Driven"},
                            {ProblemData::Type::EIGENMODE, "Eigenmode"},
                            {ProblemData::Type::ELECTROSTATIC, "Electrostatic"},
                            {ProblemData::Type::MAGNETOSTATIC, "Magnetostatic"},
                            {ProblemData::Type::TRANSIENT, "Transient"}})

void ProblemData::SetUp(json &config)
{
  auto problem = config.find("Problem");
  MFEM_VERIFY(problem != config.end(),
              "\"Problem\" must be specified in the configuration file!");
  MFEM_VERIFY(problem->find("Type") != problem->end(),
              "Missing config[\"Problem\"][\"Type\"] in the configuration file!");
  type = problem->at("Type");  // Required
  verbose = problem->value("Verbose", verbose);
  output = problem->value("Output", output);

  // Check for provided solver configuration data (not required for electrostatics or
  // magnetostatics since defaults can be used for every option).
  auto solver = config.find("Solver");
  if (type == ProblemData::Type::DRIVEN)
  {
    MFEM_VERIFY(solver->find("Driven") != solver->end(),
                "config[\"Problem\"][\"Type\"] == \"Driven\" should be accompanied by a "
                "config[\"Solver\"][\"Driven\"] configuration!");
  }
  else if (type == ProblemData::Type::EIGENMODE)
  {
    MFEM_VERIFY(solver->find("Eigenmode") != solver->end(),
                "config[\"Problem\"][\"Type\"] == \"Eigenmode\" should be accompanied by a "
                "config[\"Solver\"][\"Eigenmode\"] configuration!");
  }
  else if (type == ProblemData::Type::ELECTROSTATIC)
  {
    // MFEM_VERIFY(
    //     solver->find("Electrostatic") != solver->end(),
    //     "config[\"Problem\"][\"Type\"] == \"Electrostatic\" should be accompanied by a "
    //     "config[\"Solver\"][\"Electrostatic\"] configuration!");
  }
  else if (type == ProblemData::Type::MAGNETOSTATIC)
  {
    // MFEM_VERIFY(
    //     solver->find("Magnetostatic") != solver->end(),
    //     "config[\"Problem\"][\"Type\"] == \"Magnetostatic\" should be accompanied by a "
    //     "config[\"Solver\"][\"Magnetostatic\"] configuration!");
  }
  else if (type == ProblemData::Type::TRANSIENT)
  {
    MFEM_VERIFY(solver->find("Transient") != solver->end(),
                "config[\"Problem\"][\"Type\"] == \"Transient\" should be accompanied by a "
                "config[\"Solver\"][\"Transient\"] configuration!");
  }

  // Cleanup
  problem->erase("Type");
  problem->erase("Verbose");
  problem->erase("Output");
  MFEM_VERIFY(problem->empty(),
              "Found an unsupported configuration file keyword under \"Problem\"!\n"
                  << problem->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Type: " << type << '\n';
    std::cout << "Verbose: " << verbose << '\n';
    std::cout << "Output: " << output << '\n';
  }
}

void RefinementData::SetUp(json &model)
{
  auto refinement = model.find("Refinement");
  if (refinement == model.end())
  {
    return;
  }

  // Options for AMR.
  tol = refinement->value("Tol", tol);
  max_it = refinement->value("MaxIts", max_it);
  max_size = refinement->value("MaxSize", max_size);
  nonconformal = refinement->value("Nonconformal", nonconformal);
  max_nc_levels = refinement->value("MaxNCLevels", max_nc_levels);
  update_fraction = refinement->value("UpdateFraction", update_fraction);
  maximum_imbalance = refinement->value("MaximumImbalance", maximum_imbalance);
  save_adapt_iterations = refinement->value("SaveAdaptIterations", save_adapt_iterations);
  save_adapt_mesh = refinement->value("SaveAdaptMesh", save_adapt_mesh);
  MFEM_VERIFY(tol > 0.0, "config[\"Refinement\"][\"Tol\"] must be strictly positive!");
  MFEM_VERIFY(max_it >= 0, "config[\"Refinement\"][\"MaxIts\"] must be non-negative!");
  MFEM_VERIFY(max_size >= 0, "config[\"Refinement\"][\"MaxSize\"] must be non-negative!");
  MFEM_VERIFY(max_nc_levels >= 0,
              "config[\"Refinement\"][\"MaxNCLevels\"] must be non-negative!");
  MFEM_VERIFY(update_fraction > 0 && update_fraction < 1,
              "config[\"Refinement\"][\"UpdateFraction\" must be in (0,1)!");
  MFEM_VERIFY(
      maximum_imbalance >= 1,
      "config[\"Refinement\"][\"MaximumImbalance\"] must be greater than or equal to 1!");

  // Options for a priori refinement.
  uniform_ref_levels = refinement->value("UniformLevels", uniform_ref_levels);
  ser_uniform_ref_levels = refinement->value("SerialUniformLevels", ser_uniform_ref_levels);
  MFEM_VERIFY(uniform_ref_levels >= 0 && ser_uniform_ref_levels >= 0,
              "Number of uniform mesh refinement levels must be non-negative!");
  auto boxes = refinement->find("Boxes");
  if (boxes != refinement->end())
  {
    MFEM_VERIFY(boxes->is_array(), "config[\"Refinement\"][\"Boxes\"] should specify an "
                                   "array in the configuration file!");
    for (auto it = boxes->begin(); it != boxes->end(); ++it)
    {
      MFEM_VERIFY(
          it->find("Levels") != it->end(),
          "Missing \"Boxes\" refinement region \"Levels\" in the configuration file!");
      auto bbmin = it->find("BoundingBoxMin");
      auto bbmax = it->find("BoundingBoxMax");
      MFEM_VERIFY(bbmin != it->end() && bbmax != it->end(),
                  "Missing \"Boxes\" refinement region \"BoundingBoxMin/Max\" in the "
                  "configuration file!");
      MFEM_VERIFY(bbmin->is_array() && bbmin->is_array(),
                  "config[\"Refinement\"][\"Boxes\"][\"BoundingBoxMin/Max\"] should "
                  "specify an array in the configuration file!");
      BoxRefinementData &data = box_list.emplace_back();
      data.ref_levels = it->at("Levels");                // Required
      data.bbmin = bbmin->get<std::array<double, 3>>();  // Required
      data.bbmax = bbmax->get<std::array<double, 3>>();  // Required

      // Cleanup
      it->erase("Levels");
      it->erase("BoundingBoxMin");
      it->erase("BoundingBoxMax");
      MFEM_VERIFY(it->empty(), "Found an unsupported configuration file keyword under "
                               "config[\"Refinement\"][\"Boxes\"]!\n"
                                   << it->dump(2));

      // Debug
      if constexpr (JSON_DEBUG)
      {
        std::cout << "Levels: " << data.ref_levels << '\n';
        std::cout << "BoundingBoxMin: " << data.bbmin << '\n';
        std::cout << "BoundingBoxMax: " << data.bbmax << '\n';
      }
    }
  }
  auto spheres = refinement->find("Spheres");
  if (spheres != refinement->end())
  {
    MFEM_VERIFY(spheres->is_array(), "config[\"Refinement\"][\"Spheres\"] should specify "
                                     "an array in the configuration file!");
    for (auto it = spheres->begin(); it != spheres->end(); ++it)
    {
      MFEM_VERIFY(
          it->find("Levels") != it->end(),
          "Missing \"Spheres\" refinement region \"Levels\" in the configuration file!");
      auto ctr = it->find("Center");
      MFEM_VERIFY(ctr != it->end() && it->find("Radius") != it->end(),
                  "Missing \"Spheres\" refinement region \"Center\" or \"Radius\" in "
                  "configuration file!");
      MFEM_VERIFY(ctr->is_array(),
                  "config[\"Refinement\"][\"Spheres\"][\"Center\"] should specify "
                  "an array in the configuration file!");
      SphereRefinementData &data = sphere_list.emplace_back();
      data.ref_levels = it->at("Levels");               // Required
      data.r = it->at("Radius");                        // Required
      data.center = ctr->get<std::array<double, 3>>();  // Required

      // Cleanup
      it->erase("Levels");
      it->erase("Radius");
      it->erase("Center");
      MFEM_VERIFY(it->empty(), "Found an unsupported configuration file keyword under "
                               "config[\"Refinement\"][\"Spheres\"]!\n"
                                   << it->dump(2));

      // Debug
      if constexpr (JSON_DEBUG)
      {
        std::cout << "Levels: " << data.ref_levels << '\n';
        std::cout << "Radius: " << data.r << '\n';
        std::cout << "Center: " << data.center << '\n';
      }
    }
  }

  // Cleanup
  refinement->erase("Tol");
  refinement->erase("MaxIts");
  refinement->erase("MaxSize");
  refinement->erase("Nonconformal");
  refinement->erase("MaxNCLevels");
  refinement->erase("UpdateFraction");
  refinement->erase("MaximumImbalance");
  refinement->erase("SaveAdaptIterations");
  refinement->erase("SaveAdaptMesh");
  refinement->erase("UniformLevels");
  refinement->erase("SerialUniformLevels");
  refinement->erase("Boxes");
  refinement->erase("Spheres");
  MFEM_VERIFY(refinement->empty(),
              "Found an unsupported configuration file keyword under \"Refinement\"!\n"
                  << refinement->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Tol: " << tol << '\n';
    std::cout << "MaxIts: " << max_it << '\n';
    std::cout << "MaxSize: " << max_size << '\n';
    std::cout << "Nonconformal: " << nonconformal << '\n';
    std::cout << "MaxNCLevels: " << max_nc_levels << '\n';
    std::cout << "UpdateFraction: " << update_fraction << '\n';
    std::cout << "MaximumImbalance: " << maximum_imbalance << '\n';
    std::cout << "SaveAdaptIterations: " << save_adapt_iterations << '\n';
    std::cout << "SaveAdaptMesh: " << save_adapt_mesh << '\n';
    std::cout << "UniformLevels: " << uniform_ref_levels << '\n';
    std::cout << "SerialUniformLevels: " << ser_uniform_ref_levels << '\n';
  }
}

void ModelData::SetUp(json &config)
{
  auto model = config.find("Model");
  MFEM_VERIFY(model != config.end(),
              "\"Model\" must be specified in the configuration file!");
  MFEM_VERIFY(model->find("Mesh") != model->end(),
              "Missing config[\"Model\"][\"Mesh\"] file in the configuration file!");
  mesh = model->at("Mesh");  // Required
  L0 = model->value("L0", L0);
  Lc = model->value("Lc", Lc);
  remove_curvature = model->value("RemoveCurvature", remove_curvature);
  make_simplex = model->value("MakeSimplex", make_simplex);
  make_hex = model->value("MakeHexahedral", make_hex);
  reorder_elements = model->value("ReorderElements", reorder_elements);
  clean_unused_elements = model->value("CleanUnusedElements", clean_unused_elements);
  crack_bdr_elements = model->value("CrackInternalBoundaryElements", crack_bdr_elements);
  refine_crack_elements = model->value("RefineCrackElements", refine_crack_elements);
  crack_displ_factor = model->value("CrackDisplacementFactor", crack_displ_factor);
  add_bdr_elements = model->value("AddInterfaceBoundaryElements", add_bdr_elements);
  reorient_tet_mesh = model->value("ReorientTetMesh", reorient_tet_mesh);
  partitioning = model->value("Partitioning", partitioning);
  refinement.SetUp(*model);

  // Cleanup
  model->erase("Mesh");
  model->erase("L0");
  model->erase("Lc");
  model->erase("RemoveCurvature");
  model->erase("MakeSimplex");
  model->erase("MakeHexahedral");
  model->erase("ReorderElements");
  model->erase("CleanUnusedElements");
  model->erase("CrackInternalBoundaryElements");
  model->erase("RefineCrackElements");
  model->erase("CrackDisplacementFactor");
  model->erase("AddInterfaceBoundaryElements");
  model->erase("ReorientTetMesh");
  model->erase("Partitioning");
  model->erase("Refinement");
  MFEM_VERIFY(model->empty(),
              "Found an unsupported configuration file keyword under \"Model\"!\n"
                  << model->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Mesh: " << mesh << '\n';
    std::cout << "L0: " << L0 << '\n';
    std::cout << "Lc: " << Lc << '\n';
    std::cout << "RemoveCurvature: " << remove_curvature << '\n';
    std::cout << "MakeSimplex: " << make_simplex << '\n';
    std::cout << "MakeHexahedral: " << make_hex << '\n';
    std::cout << "ReorderElements: " << reorder_elements << '\n';
    std::cout << "CleanUnusedElements: " << clean_unused_elements << '\n';
    std::cout << "CrackInternalBoundaryElements: " << crack_bdr_elements << '\n';
    std::cout << "RefineCrackElements: " << refine_crack_elements << '\n';
    std::cout << "CrackDisplacementFactor: " << crack_displ_factor << '\n';
    std::cout << "AddInterfaceBoundaryElements: " << add_bdr_elements << '\n';
    std::cout << "ReorientTetMesh: " << reorient_tet_mesh << '\n';
    std::cout << "Partitioning: " << partitioning << '\n';
  }
}

void DomainMaterialData::SetUp(json &domains)
{
  auto materials = domains.find("Materials");
  MFEM_VERIFY(materials != domains.end() && materials->is_array(),
              "\"Materials\" must be specified as an array in the configuration file!");
  for (auto it = materials->begin(); it != materials->end(); ++it)
  {
    MFEM_VERIFY(
        it->find("Attributes") != it->end(),
        "Missing \"Attributes\" list for \"Materials\" domain in the configuration file!");
    MaterialData &data = vecdata.emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    ParseSymmetricMatrixData(*it, "Permeability", data.mu_r);
    ParseSymmetricMatrixData(*it, "Permittivity", data.epsilon_r);
    ParseSymmetricMatrixData(*it, "LossTan", data.tandelta);
    ParseSymmetricMatrixData(*it, "Conductivity", data.sigma);
    data.lambda_L = it->value("LondonDepth", data.lambda_L);

    // Cleanup
    it->erase("Attributes");
    it->erase("Permeability");
    it->erase("Permittivity");
    it->erase("LossTan");
    it->erase("Conductivity");
    it->erase("MaterialAxes");
    it->erase("LondonDepth");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Materials\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Permeability: " << data.mu_r << '\n';
      std::cout << "Permittivity: " << data.epsilon_r << '\n';
      std::cout << "LossTan: " << data.tandelta << '\n';
      std::cout << "Conductivity: " << data.sigma << '\n';
      std::cout << "LondonDepth: " << data.lambda_L << '\n';
    }
  }
}

void DomainEnergyPostData::SetUp(json &postpro)
{
  auto energy = postpro.find("Energy");
  if (energy == postpro.end())
  {
    return;
  }
  MFEM_VERIFY(energy->is_array(),
              "\"Energy\" should specify an array in the configuration file!");
  for (auto it = energy->begin(); it != energy->end(); ++it)
  {
    auto index = AtIndex(it, "\"Energy\" domain");
    MFEM_VERIFY(
        it->find("Attributes") != it->end(),
        "Missing \"Attributes\" list for \"Energy\" domain in the configuration file!");
    auto ret = mapdata.insert(std::make_pair(index, DomainEnergyData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"Energy\" domains "
                            "in the configuration file!");
    auto &data = ret.first->second;
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());

    // Cleanup
    it->erase("Index");
    it->erase("Attributes");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Energy\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "Attributes: " << data.attributes << '\n';
    }
  }
}

void ProbePostData::SetUp(json &postpro)
{
  auto probe = postpro.find("Probe");
  if (probe == postpro.end())
  {
    return;
  }
  MFEM_VERIFY(probe->is_array(),
              "\"Probe\" should specify an array in the configuration file!");
  for (auto it = probe->begin(); it != probe->end(); ++it)
  {
    auto index = AtIndex(it, "\"Probe\" point");
    auto ctr = it->find("Center");
    MFEM_VERIFY(ctr != it->end() && ctr->is_array(),
                "Missing \"Probe\" point \"Center\" or \"Center\" should specify an array "
                "in the configuration file!");
    auto ret = mapdata.insert(std::make_pair(index, ProbeData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"Probe\" points in "
                            "the configuration file!");
    auto &data = ret.first->second;
    data.center = ctr->get<std::array<double, 3>>();  // Required

    // Cleanup
    it->erase("Index");
    it->erase("Center");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Probe\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "Center: " << data.center << '\n';
    }
  }
}

void DomainPostData::SetUp(json &domains)
{
  auto postpro = domains.find("Postprocessing");
  if (postpro == domains.end())
  {
    return;
  }
  energy.SetUp(*postpro);
  probe.SetUp(*postpro);

  // Store all unique postprocessing domain attributes.
  for (const auto &[idx, data] : energy)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();

  // Cleanup
  postpro->erase("Energy");
  postpro->erase("Probe");
  MFEM_VERIFY(postpro->empty(),
              "Found an unsupported configuration file keyword under \"Postprocessing\"!\n"
                  << postpro->dump(2));
}

void DomainData::SetUp(json &config)
{
  auto domains = config.find("Domains");
  MFEM_VERIFY(domains != config.end(),
              "\"Domains\" must be specified in the configuration file!");
  materials.SetUp(*domains);
  postpro.SetUp(*domains);

  // Store all unique domain attributes.
  for (const auto &data : materials)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();
  for (const auto &attr : postpro.attributes)
  {
    MFEM_VERIFY(std::lower_bound(attributes.begin(), attributes.end(), attr) !=
                    attributes.end(),
                "Domain postprocessing can only be enabled on domains which have a "
                "corresponding \"Materials\" entry!");
  }

  // Cleanup
  domains->erase("Materials");
  domains->erase("Postprocessing");
  MFEM_VERIFY(domains->empty(),
              "Found an unsupported configuration file keyword under \"Domains\"!\n"
                  << domains->dump(2));
}

void PecBoundaryData::SetUp(json &boundaries)
{
  auto pec = boundaries.find("PEC");
  auto ground = boundaries.find("Ground");
  if (pec == boundaries.end() && ground == boundaries.end())
  {
    return;
  }
  if (pec == boundaries.end())
  {
    pec = ground;
  }
  else if (ground == boundaries.end())  // Do nothing
  {
  }
  else
  {
    MFEM_ABORT(
        "Configuration file should not specify both \"PEC\" and \"Ground\" boundaries!");
  }
  MFEM_VERIFY(
      pec->find("Attributes") != pec->end(),
      "Missing \"Attributes\" list for \"PEC\" boundary in the configuration file!");
  attributes = pec->at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());

  // Cleanup
  pec->erase("Attributes");
  MFEM_VERIFY(pec->empty(),
              "Found an unsupported configuration file keyword under \"PEC\"!\n"
                  << pec->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "PEC:" << attributes << '\n';
  }
}

void PmcBoundaryData::SetUp(json &boundaries)
{
  auto pmc = boundaries.find("PMC");
  auto zeroq = boundaries.find("ZeroCharge");
  if (pmc == boundaries.end() && zeroq == boundaries.end())
  {
    return;
  }
  if (pmc == boundaries.end())
  {
    pmc = zeroq;
  }
  else if (zeroq == boundaries.end())  // Do nothing
  {
  }
  else
  {
    MFEM_ABORT("Configuration file should not specify both \"PMC\" and \"ZeroCharge\" "
               "boundaries!");
  }
  MFEM_VERIFY(
      pmc->find("Attributes") != pmc->end(),
      "Missing \"Attributes\" list for \"PMC\" boundary in the configuration file!");
  attributes = pmc->at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());

  // Cleanup
  pmc->erase("Attributes");
  MFEM_VERIFY(pmc->empty(),
              "Found an unsupported configuration file keyword under \"PMC\"!\n"
                  << pmc->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "PMC:" << attributes << '\n';
  }
}

void WavePortPecBoundaryData::SetUp(json &boundaries)
{
  auto pec = boundaries.find("WavePortPEC");
  if (pec == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(pec->find("Attributes") != pec->end(),
              "Missing \"Attributes\" list for \"WavePortPEC\" boundary in the "
              "configuration file!");
  attributes = pec->at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());

  // Cleanup
  pec->erase("Attributes");
  MFEM_VERIFY(pec->empty(),
              "Found an unsupported configuration file keyword under \"WavePortPEC\"!\n"
                  << pec->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "WavePortPEC:" << attributes << '\n';
  }
}

void FarfieldBoundaryData::SetUp(json &boundaries)
{
  auto absorbing = boundaries.find("Absorbing");
  if (absorbing == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(
      absorbing->find("Attributes") != absorbing->end(),
      "Missing \"Attributes\" list for \"Absorbing\" boundary in the configuration file!");
  attributes = absorbing->at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  order = absorbing->value("Order", order);
  MFEM_VERIFY(order == 1 || order == 2,
              "Supported values for config[\"Absorbing\"][\"Order\"] are 1 or 2!");

  // Cleanup
  absorbing->erase("Attributes");
  absorbing->erase("Order");
  MFEM_VERIFY(absorbing->empty(),
              "Found an unsupported configuration file keyword under \"Absorbing\"!\n"
                  << absorbing->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Absorbing:" << attributes << '\n';
    std::cout << "Order: " << order << '\n';
  }
}

void ConductivityBoundaryData::SetUp(json &boundaries)
{
  auto conductivity = boundaries.find("Conductivity");
  if (conductivity == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(conductivity->is_array(),
              "\"Conductivity\" should specify an array in the configuration file!");
  for (auto it = conductivity->begin(); it != conductivity->end(); ++it)
  {
    MFEM_VERIFY(it->find("Attributes") != it->end(),
                "Missing \"Attributes\" list for \"Conductivity\" boundary in the "
                "configuration file!");
    MFEM_VERIFY(
        it->find("Conductivity") != it->end(),
        "Missing \"Conductivity\" boundary \"Conductivity\" in the configuration file!");
    ConductivityData &data = vecdata.emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.sigma = it->at("Conductivity");  // Required
    data.mu_r = it->value("Permeability", data.mu_r);
    data.h = it->value("Thickness", data.h);
    data.external = it->value("External", data.external);

    // Cleanup
    it->erase("Attributes");
    it->erase("Conductivity");
    it->erase("Permeability");
    it->erase("Thickness");
    it->erase("External");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Conductivity\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Conductivity: " << data.sigma << '\n';
      std::cout << "Permeability: " << data.mu_r << '\n';
      std::cout << "Thickness: " << data.h << '\n';
      std::cout << "External: " << data.external << '\n';
    }
  }
}

void ImpedanceBoundaryData::SetUp(json &boundaries)
{
  auto impedance = boundaries.find("Impedance");
  if (impedance == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(impedance->is_array(),
              "\"Impedance\" should specify an array in the configuration file!");
  for (auto it = impedance->begin(); it != impedance->end(); ++it)
  {
    MFEM_VERIFY(it->find("Attributes") != it->end(),
                "Missing \"Attributes\" list for \"Impedance\" boundary in the "
                "configuration file!");
    ImpedanceData &data = vecdata.emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.Rs = it->value("Rs", data.Rs);
    data.Ls = it->value("Ls", data.Ls);
    data.Cs = it->value("Cs", data.Cs);

    // Cleanup
    it->erase("Attributes");
    it->erase("Rs");
    it->erase("Ls");
    it->erase("Cs");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Impedance\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Rs: " << data.Rs << '\n';
      std::cout << "Ls: " << data.Ls << '\n';
      std::cout << "Cs: " << data.Cs << '\n';
    }
  }
}

int ParsePortExcitation(json::iterator &port_it, TriBool &excitation_used_bool_input,
                        int default_excitation = 0)
{
  auto it_excitation = port_it->find("Excitation");
  if (it_excitation == port_it->end())
  {
    // Keep default; don't set input flag
    return default_excitation;
  }
  // For error printing
  int port_idx = port_it->at("Index");

  if (it_excitation->is_boolean())
  {
    if (excitation_used_bool_input == TriBool::Uninitalized)
    {
      excitation_used_bool_input = TriBool::True;
    }
    MFEM_VERIFY(excitation_used_bool_input == TriBool::True,
                fmt::format("\"Excitation\" on port index {:d} specified with a bool after "
                            "using integers; consitantly use either one or the other!",
                            port_idx));
    return int(it_excitation->get<bool>());  // 0 false; 1 true
  }
  else if (it_excitation->is_number_integer())
  {
    if (excitation_used_bool_input == TriBool::Uninitalized)
    {
      excitation_used_bool_input = TriBool::False;
    }
    MFEM_VERIFY(
        excitation_used_bool_input == TriBool::False,
        fmt::format("\"Excitation\" on port index {:d} specified with an "
                    "integer after using bools; consitantly use either one or the other!",
                    port_idx));

    auto excitation_parse = it_excitation->get<int>();
    MFEM_VERIFY(excitation_parse >= 0,
                fmt::format("\"Excitation\" on port index {:d} should be an "
                            "positive integer (excited) or zero (not excited); got {:d}",
                            port_idx, excitation_parse));
    return excitation_parse;
  }
  else
  {
    MFEM_ABORT(fmt::format("\"Excitation\" on port index {:d} could not be parsed "
                           "as a bool or non-negative integer; got {}",
                           port_idx, it_excitation->dump(2)));
  }
}

LumpedPortBoundaryData::SetUpReturnInfo LumpedPortBoundaryData::SetUp(json &boundaries)
{
  SetUpReturnInfo out = {};
  auto port = boundaries.find("LumpedPort");
  auto terminal = boundaries.find("Terminal");
  if (port == boundaries.end() && terminal == boundaries.end())
  {
    return out;
  }

  if (port == boundaries.end())
  {
    port = terminal;
    out.has_terminal_spec = TriBool::True;
  }
  else if (terminal == boundaries.end())  // Do nothing
  {
    out.has_terminal_spec = TriBool::False;
  }
  else
  {
    MFEM_ABORT("Configuration file should not specify both \"LumpedPort\" and \"Terminal\" "
               "boundaries!");
  }
  std::string err_label =
      (out.has_terminal_spec == TriBool::True) ? "\"Terminal\"" : "\"LumpedPort\"";
  std::string err_label_b = fmt::format("{} boundary", err_label);
  MFEM_VERIFY(
      port->is_array(),
      fmt::format("{} should specify an array in the configuration file!", err_label));
  for (auto it = port->begin(); it != port->end(); ++it)
  {
    auto index = AtIndex(it, err_label);
    auto ret = mapdata.insert(std::make_pair(index, LumpedPortData()));
    MFEM_VERIFY(ret.second, fmt::format("Repeated \"Index\" found when processing {} "
                                        "boundaries in the configuration file!",
                                        err_label));
    auto &data = ret.first->second;
    data.R = it->value("R", data.R);
    data.L = it->value("L", data.L);
    data.C = it->value("C", data.C);
    data.Rs = it->value("Rs", data.Rs);
    data.Ls = it->value("Ls", data.Ls);
    data.Cs = it->value("Cs", data.Cs);

    data.excitation =
        ParsePortExcitation(it, out.excitation_input_is_bool, data.excitation);
    data.active = it->value("Active", data.active);

    if (it->find("Attributes") != it->end())
    {
      MFEM_VERIFY(
          it->find("Elements") == it->end(),
          fmt::format(
              "Cannot specify both top-level \"Attributes\" list and \"Elements\" for "
              "{} in the configuration file!",
              err_label_b));
      auto &elem = data.elements.emplace_back();
      ParseElementData(*it, "Direction", terminal == boundaries.end(), elem);
    }
    else
    {
      auto elements = it->find("Elements");
      MFEM_VERIFY(elements != it->end(),
                  fmt::format("Missing top-level \"Attributes\" list or \"Elements\" for "
                              "{} in the configuration file!",
                              err_label_b));
      for (auto elem_it = elements->begin(); elem_it != elements->end(); ++elem_it)
      {
        MFEM_VERIFY(elem_it->find("Attributes") != elem_it->end(),
                    fmt::format("Missing \"Attributes\" list for {} element in "
                                "the configuration file!",
                                err_label_b));
        auto &elem = data.elements.emplace_back();
        ParseElementData(*elem_it, "Direction", terminal == boundaries.end(), elem);

        // Cleanup
        elem_it->erase("Attributes");
        elem_it->erase("Direction");
        elem_it->erase("CoordinateSystem");
        MFEM_VERIFY(elem_it->empty(), fmt::format("Found an unsupported configuration file "
                                                  "keyword under {} element!\n{}",
                                                  err_label_b, elem_it->dump(2)));
      }
    }

    // Cleanup
    it->erase("Index");
    it->erase("R");
    it->erase("L");
    it->erase("C");
    it->erase("Rs");
    it->erase("Ls");
    it->erase("Cs");
    it->erase("Excitation");
    it->erase("Active");
    it->erase("Attributes");
    it->erase("Direction");
    it->erase("CoordinateSystem");
    it->erase("Elements");
    MFEM_VERIFY(it->empty(),
                fmt::format("Found an unsupported configuration file keyword under {}!\n{}",
                            err_label, it->dump(2)));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "R: " << data.R << '\n';
      std::cout << "L: " << data.L << '\n';
      std::cout << "C: " << data.C << '\n';
      std::cout << "Rs: " << data.Rs << '\n';
      std::cout << "Ls: " << data.Ls << '\n';
      std::cout << "Cs: " << data.Cs << '\n';
      std::cout << "Excitation: " << data.excitation << '\n';
      std::cout << "Active: " << data.active << '\n';
      for (const auto &elem : data.elements)
      {
        std::cout << "Attributes: " << elem.attributes << '\n';
        std::cout << "Direction: " << elem.direction << '\n';
        std::cout << "CoordinateSystem: " << elem.coordinate_system << '\n';
      }
    }
  }
  return out;
}

void PeriodicBoundaryData::SetUp(json &boundaries)
{
  auto periodic = boundaries.find("Periodic");
  if (periodic == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(periodic->is_array(),
              "\"Periodic\" should specify an array in the configuration file!");
  for (auto it = periodic->begin(); it != periodic->end(); ++it)
  {
    MFEM_VERIFY(it->find("DonorAttributes") != it->end(),
                "Missing \"DonorAttributes\" list for \"Periodic\" boundary in the "
                "configuration file!");
    MFEM_VERIFY(it->find("ReceiverAttributes") != it->end(),
                "Missing \"ReceiverAttributes\" list for \"Periodic\" boundary in the "
                "configuration file!");
    MFEM_VERIFY(it->find("Translation") != it->end(),
                "Missing \"Translation\" vector for \"Periodic\" boundary in the "
                "configuration file!");
    PeriodicData &data = vecdata.emplace_back();
    data.donor_attributes = it->at("DonorAttributes").get<std::vector<int>>();  // Required
    data.receiver_attributes =
        it->at("ReceiverAttributes").get<std::vector<int>>();               // Required
    data.translation = it->at("Translation").get<std::array<double, 3>>();  // Required

    // Cleanup
    it->erase("DonorAttributes");
    it->erase("ReceiverAttributes");
    it->erase("Translation");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Periodic\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "DonorAttributes: " << data.donor_attributes << '\n';
      std::cout << "ReceiverAttributes: " << data.receiver_attributes << '\n';
      std::cout << "Translation: " << data.translation << '\n';
    }
  }
}

// Helper for converting string keys to enum for WavePortData::EigenSolverType.
PALACE_JSON_SERIALIZE_ENUM(WavePortData::EigenSolverType,
                           {{WavePortData::EigenSolverType::DEFAULT, "Default"},
                            {WavePortData::EigenSolverType::SLEPC, "SLEPc"},
                            {WavePortData::EigenSolverType::ARPACK, "ARPACK"}})

WavePortBoundaryData::SetUpReturnInfo WavePortBoundaryData::SetUp(json &boundaries)
{
  SetUpReturnInfo out{};
  auto port = boundaries.find("WavePort");
  if (port == boundaries.end())
  {
    return out;
  }
  MFEM_VERIFY(port->is_array(),
              "\"WavePort\" should specify an array in the configuration file!");
  for (auto it = port->begin(); it != port->end(); ++it)
  {
    MFEM_VERIFY(
        it->find("Attributes") != it->end(),
        "Missing \"Attributes\" list for \"WavePort\" boundary in the configuration file!");
    auto index = AtIndex(it, "\"WavePort\" boundary");
    auto ret = mapdata.insert(std::make_pair(index, WavePortData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"WavePort\" "
                            "boundaries in the configuration file!");
    auto &data = ret.first->second;
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.mode_idx = it->value("Mode", data.mode_idx);
    MFEM_VERIFY(data.mode_idx > 0,
                "\"WavePort\" boundary \"Mode\" must be positive (1-based)!");
    data.d_offset = it->value("Offset", data.d_offset);
    data.eigen_type = it->value("SolverType", data.eigen_type);

    data.excitation =
        ParsePortExcitation(it, out.excitation_input_is_bool, data.excitation);
    data.active = it->value("Active", data.active);
    data.ksp_max_its = it->value("MaxIts", data.ksp_max_its);
    data.ksp_tol = it->value("KSPTol", data.ksp_tol);
    data.eig_tol = it->value("EigenTol", data.eig_tol);
    data.verbose = it->value("Verbose", data.verbose);

    // Cleanup
    it->erase("Index");
    it->erase("Attributes");
    it->erase("Mode");
    it->erase("Offset");
    it->erase("SolverType");
    it->erase("Excitation");
    it->erase("Active");
    it->erase("MaxIts");
    it->erase("KSPTol");
    it->erase("EigenTol");
    it->erase("Verbose");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"WavePort\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Mode: " << data.mode_idx << '\n';
      std::cout << "Offset: " << data.d_offset << '\n';
      std::cout << "SolverType: " << data.eigen_type << '\n';
      std::cout << "Excitation: " << data.excitation << '\n';
      std::cout << "Active: " << data.active << '\n';
      std::cout << "MaxIts: " << data.ksp_max_its << '\n';
      std::cout << "KSPTol: " << data.ksp_tol << '\n';
      std::cout << "EigenTol: " << data.eig_tol << '\n';
      std::cout << "Verbose: " << data.verbose << '\n';
    }
  }
  return out;
}

void SurfaceCurrentBoundaryData::SetUp(json &boundaries)
{
  auto source = boundaries.find("SurfaceCurrent");
  if (source == boundaries.end())
  {
    return;
  }
  MFEM_VERIFY(source->is_array(),
              "\"SurfaceCurrent\" should specify an array in the configuration file!");
  for (auto it = source->begin(); it != source->end(); ++it)
  {
    auto index = AtIndex(it, "\"SurfaceCurrent\" source");
    auto ret = mapdata.insert(std::make_pair(index, SurfaceCurrentData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"SurfaceCurrent\" "
                            "boundaries in the configuration file!");
    auto &data = ret.first->second;
    if (it->find("Attributes") != it->end())
    {
      MFEM_VERIFY(it->find("Elements") == it->end(),
                  "Cannot specify both top-level \"Attributes\" list and \"Elements\" for "
                  "\"SurfaceCurrent\" boundary in the configuration file!");
      auto &elem = data.elements.emplace_back();
      ParseElementData(*it, "Direction", true, elem);
    }
    else
    {
      auto elements = it->find("Elements");
      MFEM_VERIFY(
          elements != it->end(),
          "Missing top-level \"Attributes\" list or \"Elements\" for \"SurfaceCurrent\" "
          "boundary in the configuration file!");
      for (auto elem_it = elements->begin(); elem_it != elements->end(); ++elem_it)
      {
        MFEM_VERIFY(
            elem_it->find("Attributes") != elem_it->end(),
            "Missing \"Attributes\" list for \"SurfaceCurrent\" boundary element in "
            "configuration file!");
        auto &elem = data.elements.emplace_back();
        ParseElementData(*elem_it, "Direction", true, elem);

        // Cleanup
        elem_it->erase("Attributes");
        elem_it->erase("Direction");
        elem_it->erase("CoordinateSystem");
        MFEM_VERIFY(elem_it->empty(), "Found an unsupported configuration file keyword "
                                      "under \"SurfaceCurrent\" boundary element!\n"
                                          << elem_it->dump(2));
      }
    }

    // Cleanup
    it->erase("Index");
    it->erase("Attributes");
    it->erase("Direction");
    it->erase("Elements");
    it->erase("CoordinateSystem");
    MFEM_VERIFY(
        it->empty(),
        "Found an unsupported configuration file keyword under \"SurfaceCurrent\"!\n"
            << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      for (const auto &elem : data.elements)
      {
        std::cout << "Attributes: " << elem.attributes << '\n';
        std::cout << "Direction: " << elem.direction << '\n';
        std::cout << "CoordinateSystem: " << elem.coordinate_system << '\n';
      }
    }
  }
}

// Helper for converting string keys to enum for SurfaceFluxPostData::Type.
PALACE_JSON_SERIALIZE_ENUM(SurfaceFluxData::Type,
                           {{SurfaceFluxData::Type::ELECTRIC, "Electric"},
                            {SurfaceFluxData::Type::MAGNETIC, "Magnetic"},
                            {SurfaceFluxData::Type::POWER, "Power"}})

void SurfaceFluxPostData::SetUp(json &postpro)
{
  auto flux = postpro.find("SurfaceFlux");
  if (flux == postpro.end())
  {
    return;
  }
  MFEM_VERIFY(flux->is_array(),
              "\"SurfaceFlux\" should specify an array in the configuration file!");
  for (auto it = flux->begin(); it != flux->end(); ++it)
  {
    auto index = AtIndex(it, "\"SurfaceFlux\" boundary");
    MFEM_VERIFY(it->find("Attributes") != it->end() && it->find("Type") != it->end(),
                "Missing \"Attributes\" list or \"Type\" for \"SurfaceFlux\" boundary "
                "in the configuration file!");
    auto ret = mapdata.insert(std::make_pair(index, SurfaceFluxData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"SurfaceFlux\" "
                            "boundaries in the configuration file!");
    auto &data = ret.first->second;
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.type = it->at("Type");  // Required
    data.two_sided = it->value("TwoSided", data.two_sided);
    auto ctr = it->find("Center");
    if (ctr != it->end())
    {
      MFEM_VERIFY(ctr->is_array(),
                  "\"Center\" should specify an array in the configuration file!");
      data.center = ctr->get<std::array<double, 3>>();
      data.no_center = false;
    }

    // Cleanup
    it->erase("Index");
    it->erase("Attributes");
    it->erase("Type");
    it->erase("TwoSided");
    it->erase("Center");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"SurfaceFlux\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Type: " << data.type << '\n';
      std::cout << "TwoSided: " << data.two_sided << '\n';
      std::cout << "Center: " << data.center << '\n';
    }
  }
}

// Helper for converting string keys to enum for InterfaceDielectricData::Type.
PALACE_JSON_SERIALIZE_ENUM(InterfaceDielectricData::Type,
                           {{InterfaceDielectricData::Type::DEFAULT, "Default"},
                            {InterfaceDielectricData::Type::MA, "MA"},
                            {InterfaceDielectricData::Type::MS, "MS"},
                            {InterfaceDielectricData::Type::SA, "SA"}})

void InterfaceDielectricPostData::SetUp(json &postpro)
{
  auto dielectric = postpro.find("Dielectric");
  if (dielectric == postpro.end())
  {
    return;
  }
  MFEM_VERIFY(dielectric->is_array(),
              "\"Dielectric\" should specify an array in the configuration file!");
  for (auto it = dielectric->begin(); it != dielectric->end(); ++it)
  {
    auto index = AtIndex(it, "\"Dielectric\" boundary");
    MFEM_VERIFY(it->find("Attributes") != it->end() && it->find("Thickness") != it->end() &&
                    it->find("Permittivity") != it->end(),
                "Missing \"Dielectric\" boundary \"Attributes\" list, \"Thickness\", or "
                "\"Permittivity\" in the configuration file!");
    auto ret = mapdata.insert(std::make_pair(index, InterfaceDielectricData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"Dielectric\" "
                            "boundaries in the configuration file!");
    auto &data = ret.first->second;
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.type = it->value("Type", data.type);
    data.t = it->at("Thickness");             // Required
    data.epsilon_r = it->at("Permittivity");  // Required
    data.tandelta = it->value("LossTan", data.tandelta);

    // Cleanup
    it->erase("Index");
    it->erase("Attributes");
    it->erase("Type");
    it->erase("Thickness");
    it->erase("Permittivity");
    it->erase("LossTan");
    MFEM_VERIFY(it->empty(),
                "Found an unsupported configuration file keyword under \"Dielectric\"!\n"
                    << it->dump(2));

    // Debug
    if constexpr (JSON_DEBUG)
    {
      std::cout << "Index: " << ret.first->first << '\n';
      std::cout << "Attributes: " << data.attributes << '\n';
      std::cout << "Type: " << data.type << '\n';
      std::cout << "Thickness: " << data.t << '\n';
      std::cout << "Permittivity: " << data.epsilon_r << '\n';
      std::cout << "LossTan: " << data.tandelta << '\n';
    }
  }
}
void BoundaryPostData::SetUp(json &boundaries)
{
  auto postpro = boundaries.find("Postprocessing");
  if (postpro == boundaries.end())
  {
    return;
  }
  flux.SetUp(*postpro);
  dielectric.SetUp(*postpro);

  // Store all unique postprocessing boundary attributes.
  for (const auto &[idx, data] : flux)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : dielectric)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();

  // Cleanup
  postpro->erase("SurfaceFlux");
  postpro->erase("Dielectric");
  MFEM_VERIFY(postpro->empty(),
              "Found an unsupported configuration file keyword under \"Postprocessing\"!\n"
                  << postpro->dump(2));
}

void BoundaryData::SetUp(json &config)
{
  auto boundaries = config.find("Boundaries");
  MFEM_VERIFY(boundaries != config.end(),
              "\"Boundaries\" must be specified in the configuration file!");
  pec.SetUp(*boundaries);
  pmc.SetUp(*boundaries);
  auxpec.SetUp(*boundaries);
  farfield.SetUp(*boundaries);
  conductivity.SetUp(*boundaries);
  impedance.SetUp(*boundaries);
  auto lumpedport_setup_info = lumpedport.SetUp(*boundaries);
  periodic.SetUp(*boundaries);
  auto waveport_setup_info = waveport.SetUp(*boundaries);
  current.SetUp(*boundaries);
  postpro.SetUp(*boundaries);

  // Ensure consistent excitation specifier
  MFEM_VERIFY((lumpedport_setup_info.excitation_input_is_bool == TriBool::Uninitalized ||
               waveport_setup_info.excitation_input_is_bool == TriBool::Uninitalized ||
               lumpedport_setup_info.excitation_input_is_bool ==
                   waveport_setup_info.excitation_input_is_bool),
              "\"Excitation\" on lumped and wave ports should be specified using "
              "either integers or using bools, but not both.")

  // Ensure unique indexing of lumpedport, waveport, current
  {
    std::vector<std::pair<int, std::string_view>> index_map;
    std::string lumpedport_str = (lumpedport_setup_info.has_terminal_spec == TriBool::True)
                                     ? "Terminal"
                                     : "LumpedPort";
    std::string waveport_str = "WavePort";
    std::string current_str = "SurfaceCurrent";

    for (const auto &data : lumpedport)
    {
      index_map.emplace_back(data.first, std::string_view(lumpedport_str));
    }
    for (const auto &data : waveport)
    {
      index_map.emplace_back(data.first, std::string_view(waveport_str));
    }
    for (const auto &data : current)
    {
      index_map.emplace_back(data.first, std::string_view(current_str));
    }
    if (!index_map.empty())
    {
      std::sort(index_map.begin(), index_map.end(),
                [](auto &a, auto &b) { return a.first < b.first; });

      // Check all duplicates / triplicates for pretty-printing warning: Two lag iterator
      // Future: Can make this cleaner with C++20 std::views::chunk_by and filter.
      fmt::memory_buffer buf{};
      auto to = [&buf](auto f, auto &&...a)  // mini-lambda for cleaner code
      { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

      for (auto it = index_map.begin(); it != index_map.end(); it++)
      {
        auto it_1 = it;
        while ((++it_1 != index_map.end()) && (it_1->first == it->first))
        {
        }
        if (std::distance(it, it_1) > 1)
        {
          to(" Index {:d} specified multiple times:", it->first);
          for (; std::distance(it, it_1) >= 1; it++)
          {
            to(" {},", it->second);
          }
          buf.resize(buf.size() - 1);  // remove trailing ","
          to("\n");
        }
      }
      MFEM_VERIFY(buf.size() == 0,
                  fmt::format("Port index must be unique between \"LumpedPort\"s, "
                              "\"Terminal\"s, \"WavePort\"s, \"SurfaceCurrent\"s:\n{}",
                              fmt::to_string(buf)));
    }
  }

  // Typical usecase: If the excitation is specified as a bool and there is only a single
  // index (wave-port or lumped port) then set excitation index to be the port index, not 1.
  if (lumpedport_setup_info.excitation_input_is_bool == TriBool::True ||
      waveport_setup_info.excitation_input_is_bool == TriBool::True)
  {
    PortExcitationHelper excitation_view(lumpedport, waveport, current);
    if (excitation_view.Size() == 1)
    {
      auto ex = excitation_view.excitations.begin();
      auto n_lumped = ex->second.lumped_port.size();
      auto n_wave = ex->second.wave_port.size();
      auto n_current = ex->second.current_port.size();
      if (n_current == 0 && n_lumped == 1 && n_wave == 0)
      {
        int port_idx = ex->second.lumped_port.at(0);
        lumpedport.at(port_idx).excitation = port_idx;
      }
      else if (n_current == 0 && n_lumped == 0 && n_wave == 1)
      {
        int port_idx = ex->second.wave_port.at(0);
        waveport.at(port_idx).excitation = port_idx;
      }
    }
  }

  // Store all unique boundary attributes.
  attributes.insert(attributes.end(), pec.attributes.begin(), pec.attributes.end());
  attributes.insert(attributes.end(), pmc.attributes.begin(), pmc.attributes.end());
  attributes.insert(attributes.end(), auxpec.attributes.begin(), auxpec.attributes.end());
  attributes.insert(attributes.end(), farfield.attributes.begin(),
                    farfield.attributes.end());
  for (const auto &data : conductivity)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &data : impedance)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : lumpedport)
  {
    for (const auto &elem : data.elements)
    {
      attributes.insert(attributes.end(), elem.attributes.begin(), elem.attributes.end());
    }
  }
  for (const auto &[idx, data] : waveport)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : current)
  {
    for (const auto &elem : data.elements)
    {
      attributes.insert(attributes.end(), elem.attributes.begin(), elem.attributes.end());
    }
  }
  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();

  // Cleanup
  boundaries->erase("PEC");
  boundaries->erase("PMC");
  boundaries->erase("WavePortPEC");
  boundaries->erase("Absorbing");
  boundaries->erase("Conductivity");
  boundaries->erase("Impedance");
  boundaries->erase("LumpedPort");
  boundaries->erase("Periodic");
  boundaries->erase("WavePort");
  boundaries->erase("SurfaceCurrent");
  boundaries->erase("Ground");
  boundaries->erase("ZeroCharge");
  boundaries->erase("Terminal");
  boundaries->erase("Postprocessing");
  MFEM_VERIFY(boundaries->empty(),
              "Found an unsupported configuration file keyword under \"Boundaries\"!\n"
                  << boundaries->dump(2));
}

void DrivenSolverData::SetUp(json &solver)
{
  auto driven = solver.find("Driven");
  if (driven == solver.end())
  {
    return;
  }
  MFEM_VERIFY(driven->find("MinFreq") != driven->end() &&
                  driven->find("MaxFreq") != driven->end() &&
                  driven->find("FreqStep") != driven->end(),
              "Missing \"Driven\" solver \"MinFreq\", \"MaxFreq\", or \"FreqStep\" in "
              "configuration file!");
  min_f = driven->at("MinFreq");     // Required
  max_f = driven->at("MaxFreq");     // Required
  delta_f = driven->at("FreqStep");  // Required
  delta_post = driven->value("SaveStep", delta_post);
  rst = driven->value("Restart", rst);
  adaptive_tol = driven->value("AdaptiveTol", adaptive_tol);
  adaptive_max_size = driven->value("AdaptiveMaxSamples", adaptive_max_size);
  adaptive_memory = driven->value("AdaptiveConvergenceMemory", adaptive_memory);

  // Cleanup
  driven->erase("MinFreq");
  driven->erase("MaxFreq");
  driven->erase("FreqStep");
  driven->erase("SaveStep");
  driven->erase("Restart");
  driven->erase("AdaptiveTol");
  driven->erase("AdaptiveMaxSamples");
  driven->erase("AdaptiveConvergenceMemory");
  MFEM_VERIFY(driven->empty(),
              "Found an unsupported configuration file keyword under \"Driven\"!\n"
                  << driven->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "MinFreq: " << min_f << '\n';
    std::cout << "MaxFreq: " << max_f << '\n';
    std::cout << "FreqStep: " << delta_f << '\n';
    std::cout << "SaveStep: " << delta_post << '\n';
    std::cout << "Restart: " << rst << '\n';
    std::cout << "AdaptiveTol: " << adaptive_tol << '\n';
    std::cout << "AdaptiveMaxSamples: " << adaptive_max_size << '\n';
    std::cout << "AdaptiveConvergenceMemory: " << adaptive_memory << '\n';
  }
}

void EigenSolverData::SetUp(json &solver)
{
  auto eigenmode = solver.find("Eigenmode");
  if (eigenmode == solver.end())
  {
    return;
  }
  MFEM_VERIFY(eigenmode->find("Target") != eigenmode->end() ||
                  solver.find("Driven") != solver.end(),
              "Missing \"Eigenmode\" solver \"Target\" in the configuration file!");
  target = eigenmode->value("Target", target);  // Required (only for eigenmode simulations)
  tol = eigenmode->value("Tol", tol);
  max_it = eigenmode->value("MaxIts", max_it);
  max_size = eigenmode->value("MaxSize", max_size);
  n = eigenmode->value("N", n);
  n_post = eigenmode->value("Save", n_post);
  type = eigenmode->value("Type", type);
  pep_linear = eigenmode->value("PEPLinear", pep_linear);
  scale = eigenmode->value("Scaling", scale);
  init_v0 = eigenmode->value("StartVector", init_v0);
  init_v0_const = eigenmode->value("StartVectorConstant", init_v0_const);
  mass_orthog = eigenmode->value("MassOrthogonal", mass_orthog);

  // Cleanup
  eigenmode->erase("Target");
  eigenmode->erase("Tol");
  eigenmode->erase("MaxIts");
  eigenmode->erase("MaxSize");
  eigenmode->erase("N");
  eigenmode->erase("Save");
  eigenmode->erase("Type");
  eigenmode->erase("PEPLinear");
  eigenmode->erase("Scaling");
  eigenmode->erase("StartVector");
  eigenmode->erase("StartVectorConstant");
  eigenmode->erase("MassOrthogonal");
  MFEM_VERIFY(eigenmode->empty(),
              "Found an unsupported configuration file keyword under \"Eigenmode\"!\n"
                  << eigenmode->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Target: " << target << '\n';
    std::cout << "Tol: " << tol << '\n';
    std::cout << "MaxIts: " << max_it << '\n';
    std::cout << "MaxSize: " << max_size << '\n';
    std::cout << "N: " << n << '\n';
    std::cout << "Save: " << n_post << '\n';
    std::cout << "Type: " << type << '\n';
    std::cout << "PEPLinear: " << pep_linear << '\n';
    std::cout << "Scaling: " << scale << '\n';
    std::cout << "StartVector: " << init_v0 << '\n';
    std::cout << "StartVectorConstant: " << init_v0_const << '\n';
    std::cout << "MassOrthogonal: " << mass_orthog << '\n';
  }
}

void ElectrostaticSolverData::SetUp(json &solver)
{
  auto electrostatic = solver.find("Electrostatic");
  if (electrostatic == solver.end())
  {
    return;
  }
  n_post = electrostatic->value("Save", n_post);

  // Cleanup
  electrostatic->erase("Save");
  MFEM_VERIFY(electrostatic->empty(),
              "Found an unsupported configuration file keyword under \"Electrostatic\"!\n"
                  << electrostatic->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Save: " << n_post << '\n';
  }
}

void MagnetostaticSolverData::SetUp(json &solver)
{
  auto magnetostatic = solver.find("Magnetostatic");
  if (magnetostatic == solver.end())
  {
    return;
  }
  n_post = magnetostatic->value("Save", n_post);

  // Cleanup
  magnetostatic->erase("Save");
  MFEM_VERIFY(magnetostatic->empty(),
              "Found an unsupported configuration file keyword under \"Magnetostatic\"!\n"
                  << magnetostatic->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Save: " << n_post << '\n';
  }
}

// Helper for converting string keys to enum for TransientSolverData::Type and
// TransientSolverData::ExcitationType.
PALACE_JSON_SERIALIZE_ENUM(TransientSolverData::Type,
                           {{TransientSolverData::Type::DEFAULT, "Default"},
                            {TransientSolverData::Type::GEN_ALPHA, "GeneralizedAlpha"},
                            {TransientSolverData::Type::RUNGE_KUTTA, "RungeKutta"},
                            {TransientSolverData::Type::CVODE, "CVODE"},
                            {TransientSolverData::Type::ARKODE, "ARKODE"}})
PALACE_JSON_SERIALIZE_ENUM(
    TransientSolverData::ExcitationType,
    {{TransientSolverData::ExcitationType::SINUSOIDAL, "Sinusoidal"},
     {TransientSolverData::ExcitationType::GAUSSIAN, "Gaussian"},
     {TransientSolverData::ExcitationType::DIFF_GAUSSIAN, "DifferentiatedGaussian"},
     {TransientSolverData::ExcitationType::MOD_GAUSSIAN, "ModulatedGaussian"},
     {TransientSolverData::ExcitationType::RAMP_STEP, "Ramp"},
     {TransientSolverData::ExcitationType::SMOOTH_STEP, "SmoothStep"}})

void TransientSolverData::SetUp(json &solver)
{
  auto transient = solver.find("Transient");
  if (transient == solver.end())
  {
    return;
  }
  MFEM_VERIFY(
      transient->find("Excitation") != transient->end(),
      "Missing \"Transient\" solver \"Excitation\" type in the configuration file!");
  MFEM_VERIFY(transient->find("MaxTime") != transient->end() &&
                  transient->find("TimeStep") != transient->end(),
              "Missing \"Transient\" solver \"MaxTime\" or \"TimeStep\" in the "
              "configuration file!");
  type = transient->value("Type", type);
  excitation = transient->at("Excitation");  // Required
  pulse_f = transient->value("ExcitationFreq", pulse_f);
  pulse_tau = transient->value("ExcitationWidth", pulse_tau);
  max_t = transient->at("MaxTime");     // Required
  delta_t = transient->at("TimeStep");  // Required
  delta_post = transient->value("SaveStep", delta_post);
  order = transient->value("Order", order);
  rel_tol = transient->value("RelTol", rel_tol);
  abs_tol = transient->value("AbsTol", abs_tol);

  if (type == Type::GEN_ALPHA || type == Type::RUNGE_KUTTA)
  {
    if (transient->contains("Order"))
    {
      MFEM_WARNING("GeneralizedAlpha and RungeKutta transient solvers do not use "
                   "config[\"Transient\"][\"Order\"]!");
    }
    if (transient->contains("RelTol") || transient->contains("AbsTol"))
    {
      MFEM_WARNING(
          "GeneralizedAlpha and RungeKutta transient solvers do not use\n"
          "config[\"Transient\"][\"RelTol\"] and config[\"Transient\"][\"AbsTol\"]!");
    }
  }
  else
  {
    MFEM_VERIFY(rel_tol > 0,
                "config[\"Transient\"][\"RelTol\"] must be strictly positive!");
    MFEM_VERIFY(abs_tol > 0,
                "config[\"Transient\"][\"AbsTol\"] must be strictly positive!");
    MFEM_VERIFY(order >= 2 && order <= 5,
                "config[\"Transient\"][\"Order\"] must be between 2 and 5!");
  }

  // Cleanup
  transient->erase("Type");
  transient->erase("Excitation");
  transient->erase("ExcitationFreq");
  transient->erase("ExcitationWidth");
  transient->erase("MaxTime");
  transient->erase("TimeStep");
  transient->erase("SaveStep");
  transient->erase("Order");
  transient->erase("RelTol");
  transient->erase("AbsTol");
  MFEM_VERIFY(transient->empty(),
              "Found an unsupported configuration file keyword under \"Transient\"!\n"
                  << transient->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Type: " << type << '\n';
    std::cout << "Excitation: " << excitation << '\n';
    std::cout << "ExcitationFreq: " << pulse_f << '\n';
    std::cout << "ExcitationWidth: " << pulse_tau << '\n';
    std::cout << "MaxTime: " << max_t << '\n';
    std::cout << "TimeStep: " << delta_t << '\n';
    std::cout << "SaveStep: " << delta_post << '\n';
    std::cout << "Order: " << order << '\n';
    std::cout << "RelTol: " << rel_tol << '\n';
    std::cout << "AbsTol: " << abs_tol << '\n';
  }
}

// Helpers for converting string keys to enum for LinearSolverData::Type,
// LinearSolverData::KspType, LinearSolverData::SideType,
// LinearSolverData::MultigridCoarsenType, LinearSolverData::SymFactType,
// LinearSolverData::CompressionType, and LinearSolverData::OrthogType.
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::Type,
                           {{LinearSolverData::Type::DEFAULT, "Default"},
                            {LinearSolverData::Type::AMS, "AMS"},
                            {LinearSolverData::Type::BOOMER_AMG, "BoomerAMG"},
                            {LinearSolverData::Type::MUMPS, "MUMPS"},
                            {LinearSolverData::Type::SUPERLU, "SuperLU"},
                            {LinearSolverData::Type::STRUMPACK, "STRUMPACK"},
                            {LinearSolverData::Type::STRUMPACK_MP, "STRUMPACK-MP"},
                            {LinearSolverData::Type::JACOBI, "Jacobi"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::KspType,
                           {{LinearSolverData::KspType::DEFAULT, "Default"},
                            {LinearSolverData::KspType::CG, "CG"},
                            {LinearSolverData::KspType::MINRES, "MINRES"},
                            {LinearSolverData::KspType::GMRES, "GMRES"},
                            {LinearSolverData::KspType::FGMRES, "FGMRES"},
                            {LinearSolverData::KspType::BICGSTAB, "BiCGSTAB"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::SideType,
                           {{LinearSolverData::SideType::DEFAULT, "Default"},
                            {LinearSolverData::SideType::RIGHT, "Right"},
                            {LinearSolverData::SideType::LEFT, "Left"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::MultigridCoarsenType,
                           {{LinearSolverData::MultigridCoarsenType::LINEAR, "Linear"},
                            {LinearSolverData::MultigridCoarsenType::LOGARITHMIC,
                             "Logarithmic"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::SymFactType,
                           {{LinearSolverData::SymFactType::DEFAULT, "Default"},
                            {LinearSolverData::SymFactType::METIS, "METIS"},
                            {LinearSolverData::SymFactType::PARMETIS, "ParMETIS"},
                            {LinearSolverData::SymFactType::SCOTCH, "Scotch"},
                            {LinearSolverData::SymFactType::PTSCOTCH, "PTScotch"},
                            {LinearSolverData::SymFactType::PORD, "PORD"},
                            {LinearSolverData::SymFactType::AMD, "AMD"},
                            {LinearSolverData::SymFactType::RCM, "RCM"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::CompressionType,
                           {{LinearSolverData::CompressionType::NONE, "None"},
                            {LinearSolverData::CompressionType::BLR, "BLR"},
                            {LinearSolverData::CompressionType::HSS, "HSS"},
                            {LinearSolverData::CompressionType::HODLR, "HODLR"},
                            {LinearSolverData::CompressionType::ZFP, "ZFP"},
                            {LinearSolverData::CompressionType::BLR_HODLR, "BLR-HODLR"},
                            {LinearSolverData::CompressionType::ZFP_BLR_HODLR,
                             "ZFP-BLR-HODLR"}})
PALACE_JSON_SERIALIZE_ENUM(LinearSolverData::OrthogType,
                           {{LinearSolverData::OrthogType::MGS, "MGS"},
                            {LinearSolverData::OrthogType::CGS, "CGS"},
                            {LinearSolverData::OrthogType::CGS2, "CGS2"}})

void LinearSolverData::SetUp(json &solver)
{
  auto linear = solver.find("Linear");
  if (linear == solver.end())
  {
    return;
  }
  type = linear->value("Type", type);
  ksp_type = linear->value("KSPType", ksp_type);
  tol = linear->value("Tol", tol);
  max_it = linear->value("MaxIts", max_it);
  max_size = linear->value("MaxSize", max_size);
  initial_guess = linear->value("InitialGuess", initial_guess);

  // Options related to multigrid.
  mg_max_levels = linear->value("MGMaxLevels", mg_max_levels);
  mg_coarsen_type = linear->value("MGCoarsenType", mg_coarsen_type);
  mg_use_mesh = linear->value("MGUseMesh", mg_use_mesh);
  mg_cycle_it = linear->value("MGCycleIts", mg_cycle_it);
  mg_smooth_aux = linear->value("MGAuxiliarySmoother", mg_smooth_aux);
  mg_smooth_it = linear->value("MGSmoothIts", mg_smooth_it);
  mg_smooth_order = linear->value("MGSmoothOrder", mg_smooth_order);
  mg_smooth_sf_max = linear->value("MGSmoothEigScaleMax", mg_smooth_sf_max);
  mg_smooth_sf_min = linear->value("MGSmoothEigScaleMin", mg_smooth_sf_min);
  mg_smooth_cheby_4th = linear->value("MGSmoothChebyshev4th", mg_smooth_cheby_4th);

  // Preconditioner-specific options.
  pc_mat_real = linear->value("PCMatReal", pc_mat_real);
  pc_mat_shifted = linear->value("PCMatShifted", pc_mat_shifted);
  pc_side_type = linear->value("PCSide", pc_side_type);
  sym_fact_type = linear->value("ColumnOrdering", sym_fact_type);
  strumpack_compression_type =
      linear->value("STRUMPACKCompressionType", strumpack_compression_type);
  strumpack_lr_tol = linear->value("STRUMPACKCompressionTol", strumpack_lr_tol);
  strumpack_lossy_precision =
      linear->value("STRUMPACKLossyPrecision", strumpack_lossy_precision);
  strumpack_butterfly_l = linear->value("STRUMPACKButterflyLevels", strumpack_butterfly_l);
  superlu_3d = linear->value("SuperLU3DCommunicator", superlu_3d);
  ams_vector_interp = linear->value("AMSVectorInterpolation", ams_vector_interp);
  ams_singular_op = linear->value("AMSSingularOperator", ams_singular_op);
  amg_agg_coarsen = linear->value("AMGAggressiveCoarsening", amg_agg_coarsen);

  // Other linear solver options.
  divfree_tol = linear->value("DivFreeTol", divfree_tol);
  divfree_max_it = linear->value("DivFreeMaxIts", divfree_max_it);
  estimator_tol = linear->value("EstimatorTol", estimator_tol);
  estimator_max_it = linear->value("EstimatorMaxIts", estimator_max_it);
  estimator_mg = linear->value("EstimatorMG", estimator_mg);
  gs_orthog_type = linear->value("GSOrthogonalization", gs_orthog_type);

  // Cleanup
  linear->erase("Type");
  linear->erase("KSPType");
  linear->erase("Tol");
  linear->erase("MaxIts");
  linear->erase("MaxSize");
  linear->erase("InitialGuess");

  linear->erase("MGMaxLevels");
  linear->erase("MGCoarsenType");
  linear->erase("MGUseMesh");
  linear->erase("MGCycleIts");
  linear->erase("MGAuxiliarySmoother");
  linear->erase("MGSmoothIts");
  linear->erase("MGSmoothOrder");
  linear->erase("MGSmoothEigScaleMax");
  linear->erase("MGSmoothEigScaleMin");
  linear->erase("MGSmoothChebyshev4th");

  linear->erase("PCMatReal");
  linear->erase("PCMatShifted");
  linear->erase("PCSide");
  linear->erase("ColumnOrdering");
  linear->erase("STRUMPACKCompressionType");
  linear->erase("STRUMPACKCompressionTol");
  linear->erase("STRUMPACKLossyPrecision");
  linear->erase("STRUMPACKButterflyLevels");
  linear->erase("SuperLU3DCommunicator");
  linear->erase("AMSVectorInterpolation");
  linear->erase("AMSSingularOperator");
  linear->erase("AMGAggressiveCoarsening");

  linear->erase("DivFreeTol");
  linear->erase("DivFreeMaxIts");
  linear->erase("EstimatorTol");
  linear->erase("EstimatorMaxIts");
  linear->erase("EstimatorMG");
  linear->erase("GSOrthogonalization");
  MFEM_VERIFY(linear->empty(),
              "Found an unsupported configuration file keyword under \"Linear\"!\n"
                  << linear->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Type: " << type << '\n';
    std::cout << "KSPType: " << ksp_type << '\n';
    std::cout << "Tol: " << tol << '\n';
    std::cout << "MaxIts: " << max_it << '\n';
    std::cout << "MaxSize: " << max_size << '\n';
    std::cout << "InitialGuess: " << initial_guess << '\n';

    std::cout << "MGMaxLevels: " << mg_max_levels << '\n';
    std::cout << "MGCoarsenType: " << mg_coarsen_type << '\n';
    std::cout << "MGUseMesh: " << mg_use_mesh << '\n';
    std::cout << "MGCycleIts: " << mg_cycle_it << '\n';
    std::cout << "MGAuxiliarySmoother: " << mg_smooth_aux << '\n';
    std::cout << "MGSmoothIts: " << mg_smooth_it << '\n';
    std::cout << "MGSmoothOrder: " << mg_smooth_order << '\n';
    std::cout << "MGSmoothEigScaleMax: " << mg_smooth_sf_max << '\n';
    std::cout << "MGSmoothEigScaleMin: " << mg_smooth_sf_min << '\n';
    std::cout << "MGSmoothChebyshev4th: " << mg_smooth_cheby_4th << '\n';

    std::cout << "PCMatReal: " << pc_mat_real << '\n';
    std::cout << "PCMatShifted: " << pc_mat_shifted << '\n';
    std::cout << "PCSide: " << pc_side_type << '\n';
    std::cout << "ColumnOrdering: " << sym_fact_type << '\n';
    std::cout << "STRUMPACKCompressionType: " << strumpack_compression_type << '\n';
    std::cout << "STRUMPACKCompressionTol: " << strumpack_lr_tol << '\n';
    std::cout << "STRUMPACKLossyPrecision: " << strumpack_lossy_precision << '\n';
    std::cout << "STRUMPACKButterflyLevels: " << strumpack_butterfly_l << '\n';
    std::cout << "SuperLU3DCommunicator: " << superlu_3d << '\n';
    std::cout << "AMSVectorInterpolation: " << ams_vector_interp << '\n';
    std::cout << "AMSSingularOperator: " << ams_singular_op << '\n';
    std::cout << "AMGAggressiveCoarsening: " << amg_agg_coarsen << '\n';

    std::cout << "DivFreeTol: " << divfree_tol << '\n';
    std::cout << "DivFreeMaxIts: " << divfree_max_it << '\n';
    std::cout << "EstimatorTol: " << estimator_tol << '\n';
    std::cout << "EstimatorMaxIts: " << estimator_max_it << '\n';
    std::cout << "EstimatorMG: " << estimator_mg << '\n';
    std::cout << "GSOrthogonalization: " << gs_orthog_type << '\n';
  }
}

// Helpers for converting string keys to enum for SolverData::Device.
PALACE_JSON_SERIALIZE_ENUM(SolverData::Device, {{SolverData::Device::CPU, "CPU"},
                                                {SolverData::Device::GPU, "GPU"},
                                                {SolverData::Device::DEBUG, "Debug"}})

void SolverData::SetUp(json &config)
{
  auto solver = config.find("Solver");
  if (solver == config.end())
  {
    return;
  }
  order = solver->value("Order", order);
  pa_order_threshold = solver->value("PartialAssemblyOrder", pa_order_threshold);
  q_order_jac = solver->value("QuadratureOrderJacobian", q_order_jac);
  q_order_extra = solver->value("QuadratureOrderExtra", q_order_extra);
  device = solver->value("Device", device);
  ceed_backend = solver->value("Backend", ceed_backend);

  driven.SetUp(*solver);
  eigenmode.SetUp(*solver);
  electrostatic.SetUp(*solver);
  magnetostatic.SetUp(*solver);
  transient.SetUp(*solver);
  linear.SetUp(*solver);

  // Cleanup
  solver->erase("Order");
  solver->erase("PartialAssemblyOrder");
  solver->erase("QuadratureOrderJacobian");
  solver->erase("QuadratureOrderExtra");
  solver->erase("Device");
  solver->erase("Backend");

  solver->erase("Driven");
  solver->erase("Eigenmode");
  solver->erase("Electrostatic");
  solver->erase("Magnetostatic");
  solver->erase("Transient");
  solver->erase("Linear");
  MFEM_VERIFY(solver->empty(),
              "Found an unsupported configuration file keyword under \"Solver\"!\n"
                  << solver->dump(2));

  // Debug
  if constexpr (JSON_DEBUG)
  {
    std::cout << "Order: " << order << '\n';
    std::cout << "PartialAssemblyOrder: " << pa_order_threshold << '\n';
    std::cout << "QuadratureOrderJacobian: " << q_order_jac << '\n';
    std::cout << "QuadratureOrderExtra: " << q_order_extra << '\n';
    std::cout << "Device: " << device << '\n';
    std::cout << "Backend: " << ceed_backend << '\n';
  }
}

}  // namespace palace::config
