// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "configfile.hpp"

#include <algorithm>
#include <iterator>
#include <string_view>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

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

using json = nlohmann::json;
namespace palace
{
// Helpers for converting enums specified in labels.hpp. Must be done in palace scope rather
// than palace::config scope to ensure argument-dependent-lookup succeeds in json.

// Helper for converting string keys to enum for CoordinateSystem.
PALACE_JSON_SERIALIZE_ENUM(CoordinateSystem,
                           {{CoordinateSystem::CARTESIAN, "Cartesian"},
                            {CoordinateSystem::CYLINDRICAL, "Cylindrical"}})

// Helper for converting string keys to enum for ProblemType.
PALACE_JSON_SERIALIZE_ENUM(ProblemType,
                           {{ProblemType::DRIVEN, "Driven"},
                            {ProblemType::EIGENMODE, "Eigenmode"},
                            {ProblemType::ELECTROSTATIC, "Electrostatic"},
                            {ProblemType::MAGNETOSTATIC, "Magnetostatic"},
                            {ProblemType::TRANSIENT, "Transient"},
                            {ProblemType::MODEANALYSIS, "ModeAnalysis"}})

// Helper for converting string keys to enum for EigenSolverBackend.
PALACE_JSON_SERIALIZE_ENUM(EigenSolverBackend, {{EigenSolverBackend::DEFAULT, "Default"},
                                                {EigenSolverBackend::SLEPC, "SLEPc"},
                                                {EigenSolverBackend::ARPACK, "ARPACK"}})

// Helper for converting string keys to enum for EigenSolverBackend.
PALACE_JSON_SERIALIZE_ENUM(NonlinearEigenSolver, {{NonlinearEigenSolver::HYBRID, "Hybrid"},
                                                  {NonlinearEigenSolver::SLP, "SLP"}})

// Helper for converting string keys to enum for SurfaceFlux.
PALACE_JSON_SERIALIZE_ENUM(SurfaceFlux, {{SurfaceFlux::ELECTRIC, "Electric"},
                                         {SurfaceFlux::MAGNETIC, "Magnetic"},
                                         {SurfaceFlux::POWER, "Power"}})

// Helper for converting string keys to enum for InterfaceDielectric.
PALACE_JSON_SERIALIZE_ENUM(InterfaceDielectric, {{InterfaceDielectric::DEFAULT, "Default"},
                                                 {InterfaceDielectric::MA, "MA"},
                                                 {InterfaceDielectric::MS, "MS"},
                                                 {InterfaceDielectric::SA, "SA"}})

// Helper for converting string keys to enum for FrequencySampling.
PALACE_JSON_SERIALIZE_ENUM(FrequencySampling, {{FrequencySampling::DEFAULT, "Default"},
                                               {FrequencySampling::LINEAR, "Linear"},
                                               {FrequencySampling::LOG, "Log"},
                                               {FrequencySampling::POINT, "Point"}})

// Helper for converting string keys to enum for TimeSteppingScheme and Excitation.
PALACE_JSON_SERIALIZE_ENUM(TimeSteppingScheme,
                           {{TimeSteppingScheme::DEFAULT, "Default"},
                            {TimeSteppingScheme::GEN_ALPHA, "GeneralizedAlpha"},
                            {TimeSteppingScheme::RUNGE_KUTTA, "RungeKutta"},
                            {TimeSteppingScheme::CVODE, "CVODE"},
                            {TimeSteppingScheme::ARKODE, "ARKODE"}})
PALACE_JSON_SERIALIZE_ENUM(Excitation,
                           {{Excitation::SINUSOIDAL, "Sinusoidal"},
                            {Excitation::GAUSSIAN, "Gaussian"},
                            {Excitation::DIFF_GAUSSIAN, "DifferentiatedGaussian"},
                            {Excitation::MOD_GAUSSIAN, "ModulatedGaussian"},
                            {Excitation::RAMP_STEP, "Ramp"},
                            {Excitation::SMOOTH_STEP, "SmoothStep"}})

// Helper for converting string keys to enum for LinearSolver, KrylovSolver, and
// MultigridCoarsening
PALACE_JSON_SERIALIZE_ENUM(LinearSolver, {{LinearSolver::DEFAULT, "Default"},
                                          {LinearSolver::AMS, "AMS"},
                                          {LinearSolver::BOOMER_AMG, "BoomerAMG"},
                                          {LinearSolver::MUMPS, "MUMPS"},
                                          {LinearSolver::SUPERLU, "SuperLU"},
                                          {LinearSolver::STRUMPACK, "STRUMPACK"},
                                          {LinearSolver::STRUMPACK_MP, "STRUMPACK-MP"},
                                          {LinearSolver::JACOBI, "Jacobi"}})
PALACE_JSON_SERIALIZE_ENUM(KrylovSolver, {{KrylovSolver::DEFAULT, "Default"},
                                          {KrylovSolver::CG, "CG"},
                                          {KrylovSolver::MINRES, "MINRES"},
                                          {KrylovSolver::GMRES, "GMRES"},
                                          {KrylovSolver::FGMRES, "FGMRES"},
                                          {KrylovSolver::BICGSTAB, "BiCGSTAB"}})
PALACE_JSON_SERIALIZE_ENUM(MultigridCoarsening,
                           {{MultigridCoarsening::LINEAR, "Linear"},
                            {MultigridCoarsening::LOGARITHMIC, "Logarithmic"}})

// Helpers for converting string keys to enum for PreconditionerSide, SymbolicFactorization,
// SparseCompression, and Orthogonalization.
PALACE_JSON_SERIALIZE_ENUM(PreconditionerSide, {{PreconditionerSide::DEFAULT, "Default"},
                                                {PreconditionerSide::RIGHT, "Right"},
                                                {PreconditionerSide::LEFT, "Left"}})
PALACE_JSON_SERIALIZE_ENUM(SymbolicFactorization,
                           {{SymbolicFactorization::DEFAULT, "Default"},
                            {SymbolicFactorization::METIS, "METIS"},
                            {SymbolicFactorization::PARMETIS, "ParMETIS"},
                            {SymbolicFactorization::SCOTCH, "Scotch"},
                            {SymbolicFactorization::PTSCOTCH, "PTScotch"},
                            {SymbolicFactorization::PORD, "PORD"},
                            {SymbolicFactorization::AMD, "AMD"},
                            {SymbolicFactorization::RCM, "RCM"}})
PALACE_JSON_SERIALIZE_ENUM(SparseCompression,
                           {{SparseCompression::NONE, "None"},
                            {SparseCompression::BLR, "BLR"},
                            {SparseCompression::HSS, "HSS"},
                            {SparseCompression::HODLR, "HODLR"},
                            {SparseCompression::ZFP, "ZFP"},
                            {SparseCompression::BLR_HODLR, "BLR-HODLR"},
                            {SparseCompression::ZFP_BLR_HODLR, "ZFP-BLR-HODLR"}})
PALACE_JSON_SERIALIZE_ENUM(Orthogonalization, {{Orthogonalization::MGS, "MGS"},
                                               {Orthogonalization::CGS, "CGS"},
                                               {Orthogonalization::CGS2, "CGS2"}})

PALACE_JSON_SERIALIZE_ENUM(DomainOrthogonalizationWeight,
                           {{DomainOrthogonalizationWeight::ENERGY, "Energy"},
                            {DomainOrthogonalizationWeight::FE_BASIS_IDENTITY,
                             "FEBasisIdentity"},
                            {DomainOrthogonalizationWeight::SPACE_OVERLAP, "SpaceOverlap"}})

// Helpers for converting string keys to enum for Device.
PALACE_JSON_SERIALIZE_ENUM(Device, {{Device::CPU, "CPU"},
                                    {Device::GPU, "GPU"},
                                    {Device::DEBUG, "Debug"}})
}  // namespace palace

namespace palace::config
{

namespace
{

int AtIndex(json::const_iterator port_it, std::string_view errmsg_parent)
{
  return port_it->at("Index").get<int>();
}

template <std::size_t N>
void ParseSymmetricMatrixData(const json &mat, const std::string &name,
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
void ParseElementData(const json &elem, bool required, internal::ElementData &data)
{
  data.attributes = elem.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(data.attributes.begin(), data.attributes.end());
  auto it = elem.find("Direction");
  if (it != elem.end() && it->is_array())
  {
    data.direction = it->get<std::array<double, 3>>();
    data.coordinate_system = elem.value("CoordinateSystem", data.coordinate_system);
  }
  else
  {
    // String direction - CoordinateSystem is implicit in the string value.
    MFEM_VERIFY(elem.find("CoordinateSystem") == elem.end(),
                "Cannot specify \"CoordinateSystem\" with string \"Direction\"!");
    std::tie(data.direction, data.coordinate_system) =
        ParseStringAsDirection(elem.value("Direction", ""), required);
  }
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &data)
{
  os << fmt::format("{}", fmt::join(data, " "));
  return os;
}

template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &data)
{
  os << fmt::format("{}", fmt::join(data, " "));
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

}  // namespace

// Parse optional JSON field into type T, returns default-constructed T if missing.
template <typename T>
T ParseOptional(const json &j, const std::string &key)
{
  if (auto it = j.find(key); it != j.end())
  {
    return T(*it);
  }
  return T{};
}

// Parse optional JSON array field into std::vector<T>, returns empty vector if missing.
template <typename T>
std::vector<T> ParseOptionalVector(const json &j, const std::string &key)
{
  std::vector<T> result;
  if (auto it = j.find(key); it != j.end())
  {
    for (const auto &elem : *it)
    {
      result.emplace_back(elem);
    }
  }
  return result;
}

// Parse optional JSON array into std::map<int, T> keyed by "Index", checking for
// duplicates.
template <typename T>
std::map<int, T> ParseOptionalMap(const json &j, const std::string &key,
                                  const std::string &type_name)
{
  std::map<int, T> result;
  if (auto it = j.find(key); it != j.end())
  {
    for (auto elem = it->begin(); elem != it->end(); ++elem)
    {
      auto index = AtIndex(elem, type_name);
      auto [iter, inserted] = result.try_emplace(index, *elem);
      MFEM_VERIFY(inserted, "Repeated \"Index\" found when processing "
                                << type_name << " in the configuration file!");
    }
  }
  return result;
}

ProblemData::ProblemData(const json &problem)
{
  type = problem.at("Type");  // Required
  verbose = problem.value("Verbose", verbose);
  output = problem.value("Output", output);

  // Parse output formats.
  auto output_formats_it = problem.find("OutputFormats");
  if (output_formats_it != problem.end())
  {
    output_formats.paraview = output_formats_it->value("Paraview", output_formats.paraview);
    output_formats.gridfunction =
        output_formats_it->value("GridFunction", output_formats.gridfunction);
  }
}

RefinementData::RefinementData(const json &refinement)
{
  // Options for AMR.
  tol = refinement.value("Tol", tol);
  max_it = refinement.value("MaxIts", max_it);
  max_size = refinement.value("MaxSize", max_size);
  nonconformal = refinement.value("Nonconformal", nonconformal);
  max_nc_levels = refinement.value("MaxNCLevels", max_nc_levels);
  update_fraction = refinement.value("UpdateFraction", update_fraction);
  maximum_imbalance = refinement.value("MaximumImbalance", maximum_imbalance);
  save_adapt_iterations = refinement.value("SaveAdaptIterations", save_adapt_iterations);
  save_adapt_mesh = refinement.value("SaveAdaptMesh", save_adapt_mesh);

  // Options for a priori refinement.
  uniform_ref_levels = refinement.value("UniformLevels", uniform_ref_levels);
  ser_uniform_ref_levels = refinement.value("SerialUniformLevels", ser_uniform_ref_levels);
  auto boxes = refinement.find("Boxes");
  if (boxes != refinement.end())
  {
    for (auto it = boxes->begin(); it != boxes->end(); ++it)
    {
      auto bbmin = it->find("BoundingBoxMin");
      auto bbmax = it->find("BoundingBoxMax");
      BoxRefinementData &data = box_list.emplace_back();
      data.ref_levels = it->at("Levels");                // Required
      data.bbmin = bbmin->get<std::array<double, 3>>();  // Required
      data.bbmax = bbmax->get<std::array<double, 3>>();  // Required
    }
  }
  auto spheres = refinement.find("Spheres");
  if (spheres != refinement.end())
  {
    for (auto it = spheres->begin(); it != spheres->end(); ++it)
    {
      auto ctr = it->find("Center");
      SphereRefinementData &data = sphere_list.emplace_back();
      data.ref_levels = it->at("Levels");               // Required
      data.r = it->at("Radius");                        // Required
      data.center = ctr->get<std::array<double, 3>>();  // Required
    }
  }
}

ModelData::ModelData(const json &model)
{
  mesh = model.at("Mesh");  // Required
  L0 = model.value("L0", L0);
  Lc = model.value("Lc", Lc);
  remove_curvature = model.value("RemoveCurvature", remove_curvature);
  make_simplex = model.value("MakeSimplex", make_simplex);
  make_hex = model.value("MakeHexahedral", make_hex);
  reorder_elements = model.value("ReorderElements", reorder_elements);
  clean_unused_elements = model.value("CleanUnusedElements", clean_unused_elements);
  crack_bdr_elements = model.value("CrackInternalBoundaryElements", crack_bdr_elements);
  refine_crack_elements = model.value("RefineCrackElements", refine_crack_elements);
  crack_displ_factor = model.value("CrackDisplacementFactor", crack_displ_factor);
  add_bdr_elements = model.value("AddInterfaceBoundaryElements", add_bdr_elements);
  export_prerefined_mesh = model.value("ExportPrerefinedMesh", export_prerefined_mesh);
  reorient_tet_mesh = model.value("ReorientTetMesh", reorient_tet_mesh);
  partitioning = model.value("Partitioning", partitioning);
  refinement = ParseOptional<RefinementData>(model, "Refinement");
}

MaterialData::MaterialData(const json &domain)
{
  attributes = domain.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  ParseSymmetricMatrixData(domain, "Permeability", mu_r);
  ParseSymmetricMatrixData(domain, "Permittivity", epsilon_r);
  ParseSymmetricMatrixData(domain, "LossTan", tandelta);
  ParseSymmetricMatrixData(domain, "Conductivity", sigma);
  lambda_L = domain.value("LondonDepth", lambda_L);
}

DomainEnergyData::DomainEnergyData(const json &domain)
{
  attributes = domain.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
}

ProbeData::ProbeData(const json &probe)
{
  center = probe.at("Center").get<std::array<double, 3>>();  // Required
}

DomainPostData::DomainPostData(const json &postpro)
{
  energy = ParseOptionalMap<DomainEnergyData>(postpro, "Energy", "\"Energy\" domain");
  probe = ParseOptionalMap<ProbeData>(postpro, "Probe", "\"Probe\" point");

  // Store all unique postprocessing domain attributes.
  for (const auto &[idx, data] : energy)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();
}

CurrentDipoleData::CurrentDipoleData(const json &source)
{
  auto dir = source.find("Direction");
  if (dir->is_array())
  {
    direction = dir->get<std::array<double, 3>>();
    double norm = direction[0] * direction[0] + direction[1] * direction[1] +
                  direction[2] * direction[2];
    for (auto &x : direction)
    {
      x /= norm;
    }
  }
  else
  {
    auto direction_and_coord = ParseStringAsDirection(dir->get<std::string>());
    MFEM_VERIFY(direction_and_coord.second == CoordinateSystem::CARTESIAN,
                "\"R\" is not a valid \"Direction\" for \"CurrentDipole\"!");
    direction = direction_and_coord.first;
  }
  center = source.at("Center").get<std::array<double, 3>>();  // Required
  moment = source.at("Moment");                               // Required
}

DomainData::DomainData(const json &domains)
{
  for (const auto &d : *domains.find("Materials"))
  {
    materials.emplace_back(d);
  }
  current_dipole = ParseOptionalMap<CurrentDipoleData>(domains, "CurrentDipole",
                                                       "\"CurrentDipole\" source");
  postpro = ParseOptional<DomainPostData>(domains, "Postprocessing");

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
}

PecBoundaryData::PecBoundaryData(const json &pec)
{
  attributes = pec.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
}

PmcBoundaryData::PmcBoundaryData(const json &pmc)
{
  attributes = pmc.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
}

WavePortPecBoundaryData::WavePortPecBoundaryData(const json &auxpec)
{
  attributes = auxpec.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
}

FarfieldBoundaryData::FarfieldBoundaryData(const json &absorbing)
{
  attributes = absorbing.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  order = absorbing.value("Order", order);
}

ConductivityData::ConductivityData(const json &boundary)
{
  attributes = boundary.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  sigma = boundary.at("Conductivity");  // Required
  mu_r = boundary.value("Permeability", mu_r);
  h = boundary.value("Thickness", h);
  external = boundary.value("External", external);
}

ImpedanceData::ImpedanceData(const json &boundary)
{
  attributes = boundary.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  Rs = boundary.value("Rs", Rs);
  Ls = boundary.value("Ls", Ls);
  Cs = boundary.value("Cs", Cs);
}

int ParsePortExcitation(const json &port, int index)
{
  auto it = port.find("Excitation");
  if (it == port.end())
  {
    return 0;  // Not excited
  }
  else if (it->is_boolean())
  {
    return int(it->get<bool>());  // 0 false; 1 true
  }
  else if (it->is_number_unsigned())
  {
    return it->get<int>();
  }
  else
  {
    MFEM_ABORT(fmt::format("\"Excitation\" on port index {:d} could not be parsed "
                           "as a bool or unsigned (non-negative) integer; got {}",
                           index, it->dump(2)));
  }
}

LumpedPortData::LumpedPortData(const json &port)
{
  int index = port.at("Index");  // Required
  R = port.value("R", R);
  L = port.value("L", L);
  C = port.value("C", C);
  Rs = port.value("Rs", Rs);
  Ls = port.value("Ls", Ls);
  Cs = port.value("Cs", Cs);

  excitation = ParsePortExcitation(port, index);
  active = port.value("Active", active);
  if (port.find("Attributes") != port.end())
  {
    MFEM_VERIFY(port.find("Elements") == port.end(),
                "Cannot specify both top-level \"Attributes\" list and \"Elements\" for "
                "\"LumpedPort\" in the configuration file!");
    auto &elem = elements.emplace_back();
    ParseElementData(port, true, elem);
  }
  else
  {
    auto elems = port.find("Elements");
    MFEM_VERIFY(elems != port.end(),
                "Missing top-level \"Attributes\" list or \"Elements\" for "
                "\"LumpedPort\" in the configuration file!");
    for (const auto &e : *elems)
    {
      auto &elem = elements.emplace_back();
      ParseElementData(e, true, elem);
    }
  }
}

TerminalData::TerminalData(const json &terminal)
{
  attributes = terminal.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
}

PeriodicBoundaryData::PeriodicBoundaryData(const json &periodic)
{
  auto floquet = periodic.find("FloquetWaveVector");
  if (floquet != periodic.end())
  {
    MFEM_VERIFY(floquet->is_array(),
                "\"FloquetWaveVector\" should specify an array in the configuration file!");
    wave_vector = floquet->get<std::array<double, 3>>();
  }

  auto pairs = periodic.find("BoundaryPairs");
  MFEM_VERIFY(pairs->is_array(),
              "\"BoundaryPairs\" should specify an array in the configuration file!");
  for (auto it = pairs->begin(); it != pairs->end(); ++it)
  {
    MFEM_VERIFY(it->find("DonorAttributes") != it->end(),
                "Missing \"DonorAttributes\" list for \"Periodic\" boundary in the "
                "configuration file!");
    MFEM_VERIFY(it->find("ReceiverAttributes") != it->end(),
                "Missing \"ReceiverAttributes\" list for \"Periodic\" boundary in the "
                "configuration file!");

    PeriodicData &data = boundary_pairs.emplace_back();
    data.donor_attributes = it->at("DonorAttributes").get<std::vector<int>>();  // Required
    data.receiver_attributes =
        it->at("ReceiverAttributes").get<std::vector<int>>();  // Required
    auto translation = it->find("Translation");
    if (translation != it->end())
    {
      MFEM_VERIFY(translation->is_array(),
                  "\"Translation\" should specify an array in the configuration file!");
      std::array<double, 3> translation_array = translation->get<std::array<double, 3>>();
      for (int i = 0; i < 3; i++)
      {
        data.affine_transform[i * 4 + i] = 1.0;
        data.affine_transform[i * 4 + 3] = translation_array[i];
      }
      data.affine_transform[3 * 4 + 3] = 1.0;
    }
    auto transformation = it->find("AffineTransformation");
    if (transformation != it->end())
    {
      MFEM_VERIFY(
          transformation->is_array(),
          "\"AffineTransformation\" should specify an array in the configuration file!");
      data.affine_transform = transformation->get<std::array<double, 16>>();
    }
  }
}

WavePortData::WavePortData(const json &port)
{
  int index = port.at("Index");                                // Required
  attributes = port.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  mode_idx = port.value("Mode", mode_idx);
  d_offset = port.value("Offset", d_offset);
  eigen_solver = port.value("SolverType", eigen_solver);
  excitation = ParsePortExcitation(port, index);
  active = port.value("Active", active);
  ksp_max_its = port.value("MaxIts", ksp_max_its);
  ksp_tol = port.value("KSPTol", ksp_tol);
  eig_tol = port.value("EigenTol", eig_tol);
  max_size = port.value("MaxSize", max_size);
  verbose = port.value("Verbose", verbose);
  if (auto it = port.find("VoltagePath"); it != port.end())
  {
    for (const auto &pt : *it)
    {
      voltage_path.push_back(pt.get<std::vector<double>>());
    }
  }
  integration_order = port.value("IntegrationOrder", integration_order);
}

SurfaceCurrentData::SurfaceCurrentData(const json &source)
{
  if (source.find("Attributes") != source.end())
  {
    MFEM_VERIFY(source.find("Elements") == source.end(),
                "Cannot specify both top-level \"Attributes\" list and \"Elements\" for "
                "\"SurfaceCurrent\" boundary in the configuration file!");
    auto &elem = elements.emplace_back();
    ParseElementData(source, true, elem);
  }
  else
  {
    auto elems = source.find("Elements");
    MFEM_VERIFY(
        elems != source.end(),
        "Missing top-level \"Attributes\" list or \"Elements\" for \"SurfaceCurrent\" "
        "boundary in the configuration file!");
    for (const auto &e : *elems)
    {
      auto &elem = elements.emplace_back();
      ParseElementData(e, true, elem);
    }
  }
}

SurfaceFluxData::SurfaceFluxData(const json &flux)
{
  attributes = flux.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  type = flux.at("Type");  // Required
  two_sided = flux.value("TwoSided", two_sided);
  auto ctr = flux.find("Center");
  if (ctr != flux.end())
  {
    center = ctr->get<std::array<double, 3>>();
    no_center = false;
  }
}

InterfaceDielectricData::InterfaceDielectricData(const json &dielectric)
{
  attributes = dielectric.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());
  type = dielectric.value("Type", type);
  t = dielectric.at("Thickness");             // Required
  epsilon_r = dielectric.at("Permittivity");  // Required
  tandelta = dielectric.value("LossTan", tandelta);
}

ModeImpedanceData::ModeImpedanceData(const json &imp)
{
  if (auto it = imp.find("VoltageAttributes"); it != imp.end())
  {
    voltage_attributes = it->get<std::vector<int>>();
    std::sort(voltage_attributes.begin(), voltage_attributes.end());
  }
  if (auto it = imp.find("CurrentAttributes"); it != imp.end())
  {
    current_attributes = it->get<std::vector<int>>();
    std::sort(current_attributes.begin(), current_attributes.end());
  }
  if (auto it = imp.find("VoltagePath"); it != imp.end())
  {
    for (const auto &pt : *it)
    {
      voltage_path.push_back(pt.get<std::vector<double>>());
    }
  }
  if (auto it = imp.find("CurrentPath"); it != imp.end())
  {
    for (const auto &pt : *it)
    {
      current_path.push_back(pt.get<std::vector<double>>());
    }
  }
  integration_order = imp.value("IntegrationOrder", integration_order);
  MFEM_VERIFY(!voltage_attributes.empty() || voltage_path.size() >= 2,
              "Impedance boundary requires either \"VoltageAttributes\" or "
              "\"VoltagePath\" in the configuration file!");
}

ModeVoltageData::ModeVoltageData(const json &volt)
{
  if (auto it = volt.find("VoltageAttributes"); it != volt.end())
  {
    voltage_attributes = it->get<std::vector<int>>();
    std::sort(voltage_attributes.begin(), voltage_attributes.end());
  }
  if (auto it = volt.find("VoltagePath"); it != volt.end())
  {
    for (const auto &pt : *it)
    {
      voltage_path.push_back(pt.get<std::vector<double>>());
    }
  }
  integration_order = volt.value("IntegrationOrder", integration_order);
  MFEM_VERIFY(!voltage_attributes.empty() || voltage_path.size() >= 2,
              "Voltage boundary requires either \"VoltageAttributes\" or "
              "\"VoltagePath\" in the configuration file!");
}

FarFieldPostData::FarFieldPostData(const json &farfield)
{
  attributes = farfield.at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());

  // Generate NSample points with the following properties:
  // - If NSample >= 2, the generated points are precisely NSample, otherwise NSample = 2.
  // - The poles, the equator, and the XZ plane are always included.
  // - The points are almost uniformly on a sphere, with a small bias due to satisfying the
  //   previous condition.
  // - The points are on rings of constant theta.

  auto nsample_json = farfield.find("NSample");
  int nsample = 0;
  if (nsample_json != farfield.end())
  {
    nsample = nsample_json->get<int>();
    if (nsample > 0)
    {
      // Always include poles.
      thetaphis.emplace_back(0.0, 0.0);   // North pole.
      thetaphis.emplace_back(M_PI, 0.0);  // South pole.

      if (nsample > 2)
      {
        int remaining = nsample - 2;

        // Distribute all remaining points across rings with number weighted by the
        // local circumference.
        int n_theta = std::max(1, static_cast<int>(std::sqrt(remaining)));
        n_theta = std::min(n_theta, remaining);  // Can't have more rings than points.

        std::vector<int> points_per_level(n_theta);
        std::vector<double> sin_theta_values(n_theta);
        double total_sin_theta = 0.0;

        // Calculate sin(theta) for each ring and total (sin(theta) is proportional to the
        // circumference).
        for (int i = 0; i < n_theta; ++i)
        {
          double theta = std::acos(1.0 - 2.0 * (i + 1) / (n_theta + 1.0));
          sin_theta_values[i] = std::sin(theta);
          total_sin_theta += sin_theta_values[i];
        }

        // Distribute points proportional to sin(theta).
        int assigned_points = 0;
        for (int i = 0; i < n_theta - 1; ++i)
        {
          points_per_level[i] =
              static_cast<int>(remaining * sin_theta_values[i] / total_sin_theta + 0.5);
          assigned_points += points_per_level[i];
        }
        // Assign remaining points to last ring to ensure exact total.
        points_per_level[n_theta - 1] = remaining - assigned_points;

        for (int i = 1; i <= n_theta; ++i)
        {
          // Ensure equator and XZ plane inclusion.
          bool is_equator = (i == (n_theta + 1) / 2);
          double theta = is_equator ? M_PI / 2 : std::acos(1.0 - 2.0 * i / (n_theta + 1.0));
          int points_in_level = points_per_level[i - 1];

          for (int j = 0; j < points_in_level; ++j)
          {
            double phi = 2.0 * M_PI * j / points_in_level;

            // Force XZ plane points (phi = 0 or π).
            if (j == 0)
            {
              phi = 0.0;
            }
            else if (j == points_in_level / 2)
            {
              phi = M_PI;
            }

            thetaphis.emplace_back(theta, phi);
          }
        }
      }

      if (nsample > 2)
        // Cast to avoid compiler warnings about types.
        MFEM_ASSERT(static_cast<int>(thetaphis.size()) == nsample,
                    "Sampled number of points is not NSample!");
    }
  }

  auto thetaphis_json = farfield.find("ThetaPhis");
  if (thetaphis_json != farfield.end())
  {
    MFEM_VERIFY(thetaphis_json->is_array(),
                "\"ThetaPhis\" should specify an array in the configuration file!");

    // JSON does not support the notion of pair, so we read the theta and phis as vectors
    // of vectors, and then cast them to vectors of pairs.
    //
    // Convert to radians in the process.
    auto vec_of_vec = thetaphis_json->get<std::vector<std::vector<double>>>();
    for (const auto &vec : vec_of_vec)
    {
      thetaphis.emplace_back(vec[0] * M_PI / 180, vec[1] * M_PI / 180);
    }
  }

  // Remove duplicate entries with numerical tolerance.
  constexpr double tol = 1e-6;
  std::sort(thetaphis.begin(), thetaphis.end());
  auto it = std::unique(thetaphis.begin(), thetaphis.end(),
                        [tol](const auto &a, const auto &b)
                        {
                          // At poles (theta ≈ 0 or π), phi is irrelevant.
                          if ((std::abs(a.first) < tol || std::abs(a.first - M_PI) < tol) &&
                              (std::abs(b.first) < tol || std::abs(b.first - M_PI) < tol))
                          {
                            return std::abs(a.first - b.first) < tol;
                          }

                          // Check direct match.
                          if (std::abs(a.first - b.first) < tol)
                          {
                            double phi_diff = std::abs(a.second - b.second);
                            return phi_diff < tol || std::abs(phi_diff - 2.0 * M_PI) < tol;
                          }

                          // Check theta periodicity: (θ, φ) ≡ (π-θ, φ+π).
                          if (std::abs(a.first - (M_PI - b.first)) < tol)
                          {
                            double phi_diff = std::abs(a.second - (b.second + M_PI));
                            if (phi_diff > M_PI)
                              phi_diff = 2.0 * M_PI - phi_diff;
                            return phi_diff < tol;
                          }

                          return false;
                        });
  thetaphis.erase(it, thetaphis.end());

  if (thetaphis.empty())
  {
    MFEM_WARNING("No target points specified under farfield \"FarField\"!\n");
  }
}
BoundaryPostData::BoundaryPostData(const json &postpro)
{
  flux =
      ParseOptionalMap<SurfaceFluxData>(postpro, "SurfaceFlux", "\"SurfaceFlux\" boundary");
  dielectric = ParseOptionalMap<InterfaceDielectricData>(postpro, "Dielectric",
                                                         "\"Dielectric\" boundary");
  impedance = ParseOptionalMap<ModeImpedanceData>(postpro, "Impedance",
                                                   "\"Impedance\" boundary");
  voltage = ParseOptionalMap<ModeVoltageData>(postpro, "Voltage",
                                               "\"Voltage\" boundary");
  farfield = ParseOptional<FarFieldPostData>(postpro, "FarField");

  // Store all unique postprocessing boundary attributes.
  for (const auto &[idx, data] : flux)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : dielectric)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : impedance)
  {
    attributes.insert(attributes.end(), data.voltage_attributes.begin(),
                      data.voltage_attributes.end());
    attributes.insert(attributes.end(), data.current_attributes.begin(),
                      data.current_attributes.end());
  }

  attributes.insert(attributes.end(), farfield.attributes.begin(),
                    farfield.attributes.end());

  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();
}

BoundaryData::BoundaryData(const json &boundaries)
{
  // PEC can be specified as "PEC" or "Ground".
  auto pec_it = boundaries.find("PEC");
  auto ground_it = boundaries.find("Ground");
  MFEM_VERIFY(
      pec_it == boundaries.end() || ground_it == boundaries.end(),
      "Configuration file should not specify both \"PEC\" and \"Ground\" boundaries!");
  if (pec_it != boundaries.end())
  {
    pec = PecBoundaryData(*pec_it);
  }
  else if (ground_it != boundaries.end())
  {
    pec = PecBoundaryData(*ground_it);
  }

  // PMC can be specified as "PMC" or "ZeroCharge".
  auto pmc_it = boundaries.find("PMC");
  auto zeroq_it = boundaries.find("ZeroCharge");
  MFEM_VERIFY(pmc_it == boundaries.end() || zeroq_it == boundaries.end(),
              "Configuration file should not specify both \"PMC\" and \"ZeroCharge\" "
              "boundaries!");
  if (pmc_it != boundaries.end())
  {
    pmc = PmcBoundaryData(*pmc_it);
  }
  else if (zeroq_it != boundaries.end())
  {
    pmc = PmcBoundaryData(*zeroq_it);
  }

  auxpec = ParseOptional<WavePortPecBoundaryData>(boundaries, "WavePortPEC");
  farfield = ParseOptional<FarfieldBoundaryData>(boundaries, "Absorbing");
  conductivity = ParseOptionalVector<ConductivityData>(boundaries, "Conductivity");
  impedance = ParseOptionalVector<ImpedanceData>(boundaries, "Impedance");
  lumpedport = ParseOptionalMap<LumpedPortData>(boundaries, "LumpedPort", "\"LumpedPort\"");
  terminal = ParseOptionalMap<TerminalData>(boundaries, "Terminal", "\"Terminal\"");
  periodic = ParseOptional<PeriodicBoundaryData>(boundaries, "Periodic");
  waveport = ParseOptionalMap<WavePortData>(boundaries, "WavePort", "\"WavePort\"");
  current = ParseOptionalMap<SurfaceCurrentData>(boundaries, "SurfaceCurrent",
                                                 "\"SurfaceCurrent\"");
  postpro = ParseOptional<BoundaryPostData>(boundaries, "Postprocessing");

  // Ensure unique indexing of lumpedport, waveport, current.
  {
    std::map<int, std::string> index_map;
    std::map<int, std::vector<int>> excitation_map;
    const std::string lumpedport_str = "\"LumpedPort\"";
    const std::string waveport_str = "WavePort";
    const std::string current_str = "SurfaceCurrent";

    for (const auto &data : lumpedport)
    {
      auto result = index_map.insert({data.first, lumpedport_str});
      MFEM_VERIFY(result.second, "Duplicate \"Index\": " << data.first << " in "
                                                         << index_map[data.first] << "!");
      excitation_map[data.second.excitation].emplace_back(data.first);
    }
    for (const auto &data : waveport)
    {
      auto result = index_map.insert({data.first, waveport_str});
      MFEM_VERIFY(result.second, "Duplicate \"Index\": " << data.first << " in "
                                                         << index_map[data.first] << "!");
      excitation_map[data.second.excitation].emplace_back(data.first);
    }
    for (const auto &data : current)
    {
      auto result = index_map.insert({data.first, current_str});
      MFEM_VERIFY(result.second, "Duplicate \"Index\": " << data.first << " in "
                                                         << index_map[data.first] << "!");
    }
    // Typical usecase: If each excitation is simple, S-parameters will be calculated.
    //    If there were multiple excitations specified, check their indices match the
    //    port indices. If there was only one, assign it.
    excitation_map.erase(0);  // zeroth index is unexcited.
    bool calc_s_params = std::all_of(excitation_map.begin(), excitation_map.end(),
                                     [](const auto &x) { return x.second.size() == 1; });
    if (calc_s_params && !excitation_map.empty())
    {
      // If there's one excitation, needs to be 1 (set with bool) or the port index.
      const auto &ext1 = *excitation_map.begin();
      MFEM_VERIFY(
          (excitation_map.size() == 1 &&
           (ext1.first == 1 || ext1.second[0] == ext1.first)) ||
              std::all_of(excitation_map.begin(), excitation_map.end(),
                          [](const auto &x) { return x.first == x.second[0]; }),
          "\"Excitation\" must match \"Index\" for single ports to avoid ambiguity!");

      for (auto &[port_idx, lp] : lumpedport)
      {
        if (lp.excitation == 1)
        {
          lp.excitation = port_idx;
        }
      }
      for (auto &[port_idx, wp] : waveport)
      {
        if (wp.excitation == 1)
        {
          wp.excitation = port_idx;
        }
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
}

std::vector<double> ConstructLinearRange(double start, double end, double delta)
{
  auto n_step = GetNumSteps(start, end, delta);
  std::vector<double> f(n_step);
  std::iota(f.begin(), f.end(), 0);
  std::for_each(f.begin(), f.end(), [=](double &x) { x = start + x * delta; });
  return f;
}
std::vector<double> ConstructLinearRange(double start, double end, int n_sample)
{
  std::vector<double> f(n_sample);
  for (int i = 0; i < n_sample; i++)
  {
    f[i] = start + (double(i) / (n_sample - 1)) * (end - start);
  }
  return f;
}
std::vector<double> ConstructLogRange(double start, double end, int n_sample)
{
  MFEM_VERIFY(start > 0 && end > 0,
              "\"Type\": \"Log\" only valid for non-zero start and end!");
  std::vector<double> f(n_sample);
  double log_start = std::log10(start);
  double log_end = std::log10(end);
  for (int i = 0; i < n_sample; i++)
  {
    double log_val = log_start + (double(i) / (n_sample - 1)) * (log_end - log_start);
    f[i] = std::pow(10.0, log_val);
  }
  return f;
}

// Helper to find entry closest to x in vec, up to tol. If no match, returns end.
auto FindNearestValue(const std::vector<double> &vec, double x, double tol)
{
  // Find the first element not less than x.
  auto it = std::lower_bound(vec.begin(), vec.end(), x);
  // Check if we found an exact match or a close enough value.
  if (it != vec.end() && std::abs(*it - x) <= tol)
  {
    return it;
  }
  // If we're not at the beginning, check the previous element too.
  if (it != vec.begin())
  {
    auto prev = std::prev(it);
    if (std::abs(*prev - x) <= tol)
    {
      return prev;
    }
  }
  // No value within tol found
  return vec.end();
}

DrivenSolverData::DrivenSolverData(const json &driven)
{
  restart = driven.value("Restart", restart);
  adaptive_tol = driven.value("AdaptiveTol", adaptive_tol);
  adaptive_max_size = driven.value("AdaptiveMaxSamples", adaptive_max_size);
  adaptive_memory = driven.value("AdaptiveConvergenceMemory", adaptive_memory);
  adaptive_solver_gs_orthog_type =
      driven.value("AdaptiveGSOrthogonalization", adaptive_solver_gs_orthog_type);
  adaptive_circuit_synthesis =
      driven.value("AdaptiveCircuitSynthesis", adaptive_circuit_synthesis);
  adaptive_circuit_synthesis_domain_orthog =
      driven.value("AdaptiveCircuitSynthesisDomainOrthogonalization",
                   adaptive_circuit_synthesis_domain_orthog);

  MFEM_VERIFY(!(restart != 1 && adaptive_tol > 0.0),
              "\"Restart\" is incompatible with adaptive frequency sweep!");

  std::vector<double> save_f, prom_f;  // samples to be saved to paraview and added to prom
  // Backwards compatible top level interface.
  if (driven.find("MinFreq") != driven.end() && driven.find("MaxFreq") != driven.end() &&
      driven.find("FreqStep") != driven.end())
  {
    double min_f = driven.at("MinFreq");     // Required
    double max_f = driven.at("MaxFreq");     // Required
    double delta_f = driven.at("FreqStep");  // Required
    sample_f = ConstructLinearRange(min_f, max_f, delta_f);
    if (int save_step = driven.value("SaveStep", 0); save_step > 0)
    {
      for (std::size_t n = 0; n < sample_f.size(); n += save_step)
      {
        save_f.emplace_back(sample_f[n]);
      }
    }
  }
  if (auto freq_samples = driven.find("Samples"); freq_samples != driven.end())
  {
    for (auto &r : *freq_samples)
    {
      auto type = r.value("Type", r.find("Freq") != r.end() ? FrequencySampling::POINT
                                                            : FrequencySampling::DEFAULT);
      auto f = [&]()
      {
        switch (type)
        {
          case FrequencySampling::LINEAR:
            {
              auto min_f = r.at("MinFreq");
              auto max_f = r.at("MaxFreq");
              auto delta_f = r.value("FreqStep", 0.0);
              auto n_sample = r.value("NSample", 0);
              MFEM_VERIFY((delta_f > 0) ^ (n_sample > 0),
                          "Only one of \"FreqStep\" or \"NSample\" can be specified for "
                          "\"Type\": \"Linear\"!");
              if (delta_f > 0)
              {
                return ConstructLinearRange(min_f, max_f, delta_f);
              }
              if (n_sample > 0)
              {
                return ConstructLinearRange(min_f, max_f, n_sample);
              }
            }
          case FrequencySampling::LOG:
            {
              auto min_f = r.at("MinFreq");
              auto max_f = r.at("MaxFreq");
              auto n_sample = r.at("NSample");
              return ConstructLogRange(min_f, max_f, n_sample);
            }
          case FrequencySampling::POINT:
            return r.at("Freq").get<std::vector<double>>();
        }
        return std::vector<double>{};
      }();
      sample_f.insert(sample_f.end(), f.begin(), f.end());

      if (auto save_step = r.value("SaveStep", 0); save_step > 0)
      {
        for (std::size_t n = 0; n < f.size(); n += save_step)
        {
          save_f.emplace_back(f[n]);
        }
      }
      if (auto prom_sample = r.value("AddToPROM", false); prom_sample)
      {
        if (adaptive_tol == 0)
        {
          MFEM_WARNING("Ignoring \"AddToPROM\" for non-adaptive simulation!");
        }
        prom_f.insert(prom_f.end(), f.begin(), f.end());
      }
    }
  }

  // Deduplicate all samples, and find indices of save and prom samples.
  constexpr double delta_eps = 1.0e-9;  // Precision in frequency comparisons (Hz)
  auto equal_f = [=](auto x, auto y) { return std::abs(x - y) < delta_eps; };
  auto deduplicate = [&equal_f](auto &f)
  {
    std::sort(f.begin(), f.end());
    f.erase(std::unique(f.begin(), f.end(), equal_f), f.end());
  };

  // Enforce explicit saves exactly match the sample frequencies.
  deduplicate(sample_f);
  auto explicit_save_f = driven.value("Save", std::vector<double>());
  for (auto &f : explicit_save_f)
  {
    auto it = FindNearestValue(sample_f, f, delta_eps);
    MFEM_VERIFY(it != sample_f.end(),
                "Entry " << f << " in \"Save\" must be an explicitly sampled frequency!");
    f = *it;
  }
  save_f.insert(save_f.end(), explicit_save_f.begin(), explicit_save_f.end());
  deduplicate(save_f);
  deduplicate(prom_f);

  // Given the matched ordering, and values are assigned by copying, can do a
  // paired-iterator scan.
  for (auto it_sample = sample_f.begin(), it_save = save_f.begin(); it_save != save_f.end();
       ++it_save)
  {
    while (*it_sample != *it_save)  // safe because save samples is a subset of samples
    {
      ++it_sample;
      MFEM_VERIFY(it_sample != sample_f.end(),
                  "Save frequency " << *it_save << " not found in sample frequencies!");
    }
    save_indices.emplace_back(std::distance(sample_f.begin(), it_sample));
  }
  // PROM sampling always begins with the minimum and maximum frequencies. Exclude them from
  // extra samples. Can use equality comparison given no floating point operations have been
  // done.
  prom_indices = {0, sample_f.size() - 1};
  if (prom_f.size() > 0 && prom_f.back() == sample_f.back())
  {
    prom_f.pop_back();
  }
  if (prom_f.size() > 0 && prom_f.front() == sample_f.front())
  {
    prom_f.erase(prom_f.begin(), std::next(prom_f.begin()));
  }
  for (auto it_sample = sample_f.begin(), it_prom = prom_f.begin(); it_prom != prom_f.end();
       ++it_prom)
  {
    while (*it_sample != *it_prom)  // safe because prom samples is a subset of samples
    {
      ++it_sample;
      MFEM_VERIFY(it_sample != sample_f.end(), "PROM sample frequency "
                                                   << *it_prom
                                                   << " not found in sample frequencies!");
    }
    prom_indices.emplace_back(std::distance(sample_f.begin(), it_sample));
  }

  MFEM_VERIFY(!sample_f.empty(), "No sample frequency samples specified in \"Driven\"!");
}

EigenSolverData::EigenSolverData(const json &eigenmode)
{
  target = eigenmode.value("Target", target);
  tol = eigenmode.value("Tol", tol);
  max_it = eigenmode.value("MaxIts", max_it);
  max_size = eigenmode.value("MaxSize", max_size);
  n = eigenmode.value("N", n);
  n_post = eigenmode.value("Save", n_post);
  type = eigenmode.value("Type", type);
  pep_linear = eigenmode.value("PEPLinear", pep_linear);
  scale = eigenmode.value("Scaling", scale);
  init_v0 = eigenmode.value("StartVector", init_v0);
  init_v0_const = eigenmode.value("StartVectorConstant", init_v0_const);
  mass_orthog = eigenmode.value("MassOrthogonal", mass_orthog);
  nonlinear_type = eigenmode.value("NonlinearType", nonlinear_type);
  refine_nonlinear = eigenmode.value("RefineNonlinear", refine_nonlinear);
  linear_tol = eigenmode.value("LinearTol", linear_tol);
  target_upper = eigenmode.value("TargetUpper", target_upper);
  preconditioner_lag = eigenmode.value("PreconditionerLag", preconditioner_lag);
  preconditioner_lag_tol = eigenmode.value("PreconditionerLagTol", preconditioner_lag_tol);
  max_restart = eigenmode.value("MaxRestart", max_restart);

  target_upper = (target_upper < 0) ? 3 * target : target_upper;  // default = 3 * target
  MFEM_VERIFY(target_upper > target, "config[\"Eigenmode\"][\"TargetUpper\"] must be "
                                     "greater than config[\"Eigenmode\"][\"Target\"]!");
}

ElectrostaticSolverData::ElectrostaticSolverData(const json &electrostatic)
{
  n_post = electrostatic.value("Save", n_post);
}

MagnetostaticSolverData::MagnetostaticSolverData(const json &magnetostatic)
{
  n_post = magnetostatic.value("Save", n_post);
}

TransientSolverData::TransientSolverData(const json &transient)
{
  type = transient.value("Type", type);
  excitation = transient.at("Excitation");  // Required
  pulse_f = transient.value("ExcitationFreq", pulse_f);
  pulse_tau = transient.value("ExcitationWidth", pulse_tau);
  max_t = transient.at("MaxTime");     // Required
  delta_t = transient.at("TimeStep");  // Required
  delta_post = transient.value("SaveStep", delta_post);
  order = transient.value("Order", order);
  rel_tol = transient.value("RelTol", rel_tol);
  abs_tol = transient.value("AbsTol", abs_tol);

  if (type == TimeSteppingScheme::GEN_ALPHA || type == TimeSteppingScheme::RUNGE_KUTTA)
  {
    if (transient.contains("Order"))
    {
      MFEM_WARNING("GeneralizedAlpha and RungeKutta transient solvers do not use "
                   "config[\"Transient\"][\"Order\"]!");
    }
    if (transient.contains("RelTol") || transient.contains("AbsTol"))
    {
      MFEM_WARNING(
          "GeneralizedAlpha and RungeKutta transient solvers do not use\n"
          "config[\"Transient\"][\"RelTol\"] and config[\"Transient\"][\"AbsTol\"]!");
    }
  }
}

ModeAnalysisSolverData::ModeAnalysisSolverData(const json &ma)
{
  freq = ma.at("Freq");  // Required
  n = ma.value("N", n);
  n_post = ma.value("Save", n_post);
  target = ma.value("Target", target);
  tol = ma.value("Tol", tol);
  max_size = ma.value("MaxSize", max_size);
  type = ma.value("Type", type);
  if (auto it = ma.find("Attributes"); it != ma.end())
  {
    attributes = it->get<std::vector<int>>();
    std::sort(attributes.begin(), attributes.end());
  }
}

LinearSolverData::LinearSolverData(const json &linear)
{
  type = linear.value("Type", type);
  krylov_solver = linear.value("KSPType", krylov_solver);
  tol = linear.value("Tol", tol);
  max_it = linear.value("MaxIts", max_it);
  max_size = linear.value("MaxSize", max_size);
  initial_guess = linear.value("InitialGuess", initial_guess);

  // Options related to multigrid.
  mg_max_levels = linear.value("MGMaxLevels", mg_max_levels);
  mg_coarsening = linear.value("MGCoarsenType", mg_coarsening);
  mg_use_mesh = linear.value("MGUseMesh", mg_use_mesh);
  mg_cycle_it = linear.value("MGCycleIts", mg_cycle_it);
  mg_smooth_aux = linear.value("MGAuxiliarySmoother", mg_smooth_aux);
  mg_smooth_it = linear.value("MGSmoothIts", mg_smooth_it);
  mg_smooth_order = linear.value("MGSmoothOrder", mg_smooth_order);
  mg_smooth_sf_max = linear.value("MGSmoothEigScaleMax", mg_smooth_sf_max);
  mg_smooth_sf_min = linear.value("MGSmoothEigScaleMin", mg_smooth_sf_min);
  mg_smooth_cheby_4th = linear.value("MGSmoothChebyshev4th", mg_smooth_cheby_4th);

  // Preconditioner-specific options.
  pc_mat_real = linear.value("PCMatReal", pc_mat_real);
  pc_mat_shifted = linear.value("PCMatShifted", pc_mat_shifted);
  complex_coarse_solve = linear.value("ComplexCoarseSolve", complex_coarse_solve);
  drop_small_entries = linear.value("DropSmallEntries", drop_small_entries);
  reorder_reuse = linear.value("ReorderingReuse", reorder_reuse);
  pc_side = linear.value("PCSide", pc_side);
  sym_factorization = linear.value("ColumnOrdering", sym_factorization);
  strumpack_compression_type =
      linear.value("STRUMPACKCompressionType", strumpack_compression_type);
  strumpack_lr_tol = linear.value("STRUMPACKCompressionTol", strumpack_lr_tol);
  strumpack_lossy_precision =
      linear.value("STRUMPACKLossyPrecision", strumpack_lossy_precision);
  strumpack_butterfly_l = linear.value("STRUMPACKButterflyLevels", strumpack_butterfly_l);
  superlu_3d = linear.value("SuperLU3DCommunicator", superlu_3d);
  ams_vector_interp = linear.value("AMSVectorInterpolation", ams_vector_interp);
  ams_singular_op = linear.value("AMSSingularOperator", ams_singular_op);
  amg_agg_coarsen = linear.value("AMGAggressiveCoarsening", amg_agg_coarsen);
  ams_max_it = linear.value("AMSMaxIts", ams_max_it);

  // Other linear solver options.
  divfree_tol = linear.value("DivFreeTol", divfree_tol);
  divfree_max_it = linear.value("DivFreeMaxIts", divfree_max_it);
  estimator_tol = linear.value("EstimatorTol", estimator_tol);
  estimator_max_it = linear.value("EstimatorMaxIts", estimator_max_it);
  estimator_mg = linear.value("EstimatorMG", estimator_mg);
  gs_orthog = linear.value("GSOrthogonalization", gs_orthog);
}

SolverData::SolverData(const json &solver)
{
  order = solver.value("Order", order);
  pa_order_threshold = solver.value("PartialAssemblyOrder", pa_order_threshold);
  q_order_jac = solver.value("QuadratureOrderJacobian", q_order_jac);
  q_order_extra = solver.value("QuadratureOrderExtra", q_order_extra);
  device = solver.value("Device", device);
  ceed_backend = solver.value("Backend", ceed_backend);

  driven = ParseOptional<DrivenSolverData>(solver, "Driven");
  eigenmode = ParseOptional<EigenSolverData>(solver, "Eigenmode");
  electrostatic = ParseOptional<ElectrostaticSolverData>(solver, "Electrostatic");
  magnetostatic = ParseOptional<MagnetostaticSolverData>(solver, "Magnetostatic");
  transient = ParseOptional<TransientSolverData>(solver, "Transient");
  mode_analysis = ParseOptional<ModeAnalysisSolverData>(solver, "ModeAnalysis");
  linear = ParseOptional<LinearSolverData>(solver, "Linear");
}

int GetNumSteps(double start, double end, double delta)
{
  if (end < start)
  {
    return 1;
  }
  constexpr double delta_eps = 1.0e-9;  // 9 digits of precision comparing endpoint
  double dn = std::abs(end - start) / std::abs(delta);
  int n_step = 1 + static_cast<int>(dn);
  double dfinal = start + n_step * delta;
  return n_step + ((delta < 0.0 && dfinal - end > -delta_eps * end) ||
                   (delta > 0.0 && dfinal - end < delta_eps * end));
}

std::pair<std::array<double, 3>, CoordinateSystem> ParseStringAsDirection(std::string str,
                                                                          bool required)
{
  if (str.empty())
  {
    MFEM_VERIFY(!required, "Missing required \"Direction\" in the configuration file!");
    return {std::array{0.0, 0.0, 0.0}, CoordinateSystem::CARTESIAN};
  }
  const bool is_positive = str.length() == 1 || str[0] == '+';
  const char axis = std::tolower(str.back());
  switch (axis)
  {
    case 'x':
      return {std::array{is_positive ? 1.0 : -1.0, 0.0, 0.0}, CoordinateSystem::CARTESIAN};
    case 'y':
      return {std::array{0.0, is_positive ? 1.0 : -1.0, 0.0}, CoordinateSystem::CARTESIAN};
    case 'z':
      return {std::array{0.0, 0.0, is_positive ? 1.0 : -1.0}, CoordinateSystem::CARTESIAN};
    case 'r':
      return {std::array{is_positive ? 1.0 : -1.0, 0.0, 0.0},
              CoordinateSystem::CYLINDRICAL};
    default:
      return {std::array{0.0, 0.0, 0.0}, CoordinateSystem::CARTESIAN};
  }
}

}  // namespace palace::config
