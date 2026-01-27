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
PALACE_JSON_SERIALIZE_ENUM(ProblemType, {{ProblemType::DRIVEN, "Driven"},
                                         {ProblemType::EIGENMODE, "Eigenmode"},
                                         {ProblemType::ELECTROSTATIC, "Electrostatic"},
                                         {ProblemType::MAGNETOSTATIC, "Magnetostatic"},
                                         {ProblemType::TRANSIENT, "Transient"}})

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

void ProblemData::SetUp(const json &config)
{
  auto problem = config.find("Problem");
  MFEM_VERIFY(problem != config.end(),
              "\"Problem\" must be specified in the configuration file!");
  MFEM_VERIFY(problem->find("Type") != problem->end(),
              "Missing config[\"Problem\"][\"Type\"] in the configuration file!");
  type = problem->at("Type");  // Required
  verbose = problem->value("Verbose", verbose);
  output = problem->value("Output", output);

  // Parse output formats.
  auto output_formats_it = problem->find("OutputFormats");
  if (output_formats_it != problem->end())
  {
    output_formats.paraview = output_formats_it->value("Paraview", output_formats.paraview);
    output_formats.gridfunction =
        output_formats_it->value("GridFunction", output_formats.gridfunction);
  }
}

void RefinementData::SetUp(const json &model)
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

  // Options for a priori refinement.
  uniform_ref_levels = refinement->value("UniformLevels", uniform_ref_levels);
  ser_uniform_ref_levels = refinement->value("SerialUniformLevels", ser_uniform_ref_levels);
  auto boxes = refinement->find("Boxes");
  if (boxes != refinement->end())
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
  auto spheres = refinement->find("Spheres");
  if (spheres != refinement->end())
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

void ModelData::SetUp(const json &config)
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
  export_prerefined_mesh = model->value("ExportPrerefinedMesh", export_prerefined_mesh);
  reorient_tet_mesh = model->value("ReorientTetMesh", reorient_tet_mesh);
  partitioning = model->value("Partitioning", partitioning);
  refinement.SetUp(*model);
}

void DomainMaterialData::SetUp(const json &domains)
{
  auto materials = domains.find("Materials");
  if (materials == domains.end())
  {
    return;
  }
  for (auto it = materials->begin(); it != materials->end(); ++it)
  {
    MaterialData &data = emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    ParseSymmetricMatrixData(*it, "Permeability", data.mu_r);
    ParseSymmetricMatrixData(*it, "Permittivity", data.epsilon_r);
    ParseSymmetricMatrixData(*it, "LossTan", data.tandelta);
    ParseSymmetricMatrixData(*it, "Conductivity", data.sigma);
    data.lambda_L = it->value("LondonDepth", data.lambda_L);
  }
}

void DomainEnergyPostData::SetUp(const json &postpro)
{
  auto energy = postpro.find("Energy");
  if (energy == postpro.end())
  {
    return;
  }
  for (auto it = energy->begin(); it != energy->end(); ++it)
  {
    auto index = AtIndex(it, "\"Energy\" domain");
    auto ret = mapdata.insert(std::make_pair(index, DomainEnergyData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"Energy\" domains "
                            "in the configuration file!");
    auto &data = ret.first->second;
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
  }
}

void ProbePostData::SetUp(const json &postpro)
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
  }
}

void DomainPostData::SetUp(const json &domains)
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
}

void CurrentDipoleSourceData::SetUp(const json &domains)
{
  auto current_dipole = domains.find("CurrentDipole");
  if (current_dipole == domains.end())
  {
    return;
  }
  MFEM_VERIFY(current_dipole->is_array(),
              "\"CurrentDipole\" should specify an array in the configuration file!");

  for (auto it = current_dipole->begin(); it != current_dipole->end(); ++it)
  {
    auto index = AtIndex(it, "\"CurrentDipole\" source");
    auto ret = mapdata.insert(std::make_pair(index, CurrentDipoleData()));
    MFEM_VERIFY(ret.second, "Repeated \"Index\" found when processing \"CurrentDipole\" "
                            "sources in the configuration file!");
    auto &data = ret.first->second;

    MFEM_VERIFY(
        it->find("Direction") != it->end(),
        "Missing \"CurrentDipole\" source \"Direction\" in the configuration file!");
    MFEM_VERIFY(it->find("Center") != it->end(),
                "Missing \"CurrentDipole\" source \"Center\" in the configuration file!");
    MFEM_VERIFY(
        it->find("Moment") != it->end(),
        "Missing \"CurrentDipole\" source \"Moment\" magnitude in the configuration file!");
    auto direction = it->find("Direction");
    auto center = it->find("Center");
    MFEM_VERIFY(center->is_array(),
                "\"CurrentDipole\" source \"Center\" should specify an array "
                "in the configuration file!");

    if (direction->is_array())
    {
      // Attempt to parse as an array.
      data.direction = direction->get<std::array<double, 3>>();
      double norm = data.direction[0] * data.direction[0] +
                    data.direction[1] * data.direction[1] +
                    data.direction[2] * data.direction[2];
      for (auto &x : data.direction)
        x /= norm;
    }
    else
    {
      auto direction_and_coord =
          ParseStringAsDirection(direction->get<std::string>());  // Required
      MFEM_VERIFY(direction_and_coord.second == CoordinateSystem::CARTESIAN,
                  "\"R\" is not a valid \"Direction\" for \"CurrentDipole\"!");
      data.direction = direction_and_coord.first;
    }
    data.center = center->get<std::array<double, 3>>();  // Required
    data.moment = it->at("Moment");                      // Required
  }
}

void DomainData::SetUp(const json &config)
{
  auto domains = config.find("Domains");
  MFEM_VERIFY(domains != config.end(),
              "\"Domains\" must be specified in the configuration file!");
  materials.SetUp(*domains);
  current_dipole.SetUp(*domains);
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
}

void PecBoundaryData::SetUp(const json &boundaries)
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
}

void PmcBoundaryData::SetUp(const json &boundaries)
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
}

void WavePortPecBoundaryData::SetUp(const json &boundaries)
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
}

void FarfieldBoundaryData::SetUp(const json &boundaries)
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
}

void ConductivityBoundaryData::SetUp(const json &boundaries)
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
    ConductivityData &data = emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.sigma = it->at("Conductivity");  // Required
    data.mu_r = it->value("Permeability", data.mu_r);
    data.h = it->value("Thickness", data.h);
    data.external = it->value("External", data.external);
  }
}

void ImpedanceBoundaryData::SetUp(const json &boundaries)
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
    ImpedanceData &data = emplace_back();
    data.attributes = it->at("Attributes").get<std::vector<int>>();  // Required
    std::sort(data.attributes.begin(), data.attributes.end());
    data.Rs = it->value("Rs", data.Rs);
    data.Ls = it->value("Ls", data.Ls);
    data.Cs = it->value("Cs", data.Cs);
  }
}

int ParsePortExcitation(json::const_iterator port_it, int default_excitation)
{
  auto it_excitation = port_it->find("Excitation");
  if (it_excitation == port_it->end())
  {
    // Keep default; don't set input flag.
    return default_excitation;
  }
  else if (it_excitation->is_boolean())
  {
    return int(it_excitation->get<bool>());  // 0 false; 1 true
  }
  else if (it_excitation->is_number_unsigned())
  {
    return it_excitation->get<int>();
  }
  else
  {
    MFEM_ABORT(fmt::format("\"Excitation\" on port index {:d} could not be parsed "
                           "as a bool or unsigned (non-negative) integer; got {}",
                           int(port_it->at("Index")), it_excitation->dump(2)));
  }
}

void LumpedPortBoundaryData::SetUp(const json &boundaries)
{
  auto port = boundaries.find("LumpedPort");
  auto terminal = boundaries.find("Terminal");
  if (port == boundaries.end() && terminal == boundaries.end())
  {
    return;
  }

  if (port == boundaries.end())
  {
    port = terminal;
  }
  else if (terminal == boundaries.end())  // Do nothing
  {
  }
  else
  {
    MFEM_ABORT("Configuration file should not specify both \"LumpedPort\" and \"Terminal\" "
               "boundaries!");
  }

  std::string label = (terminal != boundaries.end()) ? "\"Terminal\"" : "\"LumpedPort\"";
  MFEM_VERIFY(port->is_array(),
              label << " should specify an array in the configuration file!");
  for (auto it = port->begin(); it != port->end(); ++it)
  {
    auto index = AtIndex(it, label);
    auto ret = mapdata.insert(std::make_pair(index, LumpedPortData()));
    MFEM_VERIFY(ret.second, fmt::format("Repeated \"Index\" found when processing {} "
                                        "boundaries in the configuration file!",
                                        label));
    auto &data = ret.first->second;
    data.R = it->value("R", data.R);
    data.L = it->value("L", data.L);
    data.C = it->value("C", data.C);
    data.Rs = it->value("Rs", data.Rs);
    data.Ls = it->value("Ls", data.Ls);
    data.Cs = it->value("Cs", data.Cs);

    data.excitation = ParsePortExcitation(it, data.excitation);
    data.active = it->value("Active", data.active);
    if (it->find("Attributes") != it->end())
    {
      MFEM_VERIFY(
          it->find("Elements") == it->end(),
          fmt::format(
              "Cannot specify both top-level \"Attributes\" list and \"Elements\" for "
              "{} in the configuration file!",
              label));
      auto &elem = data.elements.emplace_back();
      ParseElementData(*it, terminal == boundaries.end(), elem);
    }
    else
    {
      auto elements = it->find("Elements");
      MFEM_VERIFY(elements != it->end(),
                  fmt::format("Missing top-level \"Attributes\" list or \"Elements\" for "
                              "{} in the configuration file!",
                              label));
      for (auto elem_it = elements->begin(); elem_it != elements->end(); ++elem_it)
      {
        MFEM_VERIFY(elem_it->find("Attributes") != elem_it->end(),
                    fmt::format("Missing \"Attributes\" list for {} element in "
                                "the configuration file!",
                                label));
        auto &elem = data.elements.emplace_back();
        ParseElementData(*elem_it, terminal == boundaries.end(), elem);
      }
    }
  }
}

void PeriodicBoundaryData::SetUp(const json &boundaries)
{
  auto periodic = boundaries.find("Periodic");
  if (periodic == boundaries.end())
  {
    return;
  }
  auto floquet = periodic->find("FloquetWaveVector");
  if (floquet != periodic->end())
  {
    MFEM_VERIFY(floquet->is_array(),
                "\"FloquetWaveVector\" should specify an array in the configuration file!");
    wave_vector = floquet->get<std::array<double, 3>>();
  }

  auto pairs = periodic->find("BoundaryPairs");
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

void WavePortBoundaryData::SetUp(const json &boundaries)
{
  auto port = boundaries.find("WavePort");
  if (port == boundaries.end())
  {
    return;
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
    data.d_offset = it->value("Offset", data.d_offset);
    data.eigen_solver = it->value("SolverType", data.eigen_solver);

    data.excitation = ParsePortExcitation(it, data.excitation);
    data.active = it->value("Active", data.active);
    data.ksp_max_its = it->value("MaxIts", data.ksp_max_its);
    data.ksp_tol = it->value("KSPTol", data.ksp_tol);
    data.eig_tol = it->value("EigenTol", data.eig_tol);
    data.verbose = it->value("Verbose", data.verbose);
  }
}

void SurfaceCurrentBoundaryData::SetUp(const json &boundaries)
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
      ParseElementData(*it, true, elem);
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
        ParseElementData(*elem_it, true, elem);
      }
    }
  }
}

void SurfaceFluxPostData::SetUp(const json &postpro)
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
  }
}

void InterfaceDielectricPostData::SetUp(const json &postpro)
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
  }
}
void FarFieldPostData::SetUp(const json &postpro)
{
  auto farfield = postpro.find("FarField");
  if (farfield == postpro.end())
  {
    return;
  }

  MFEM_VERIFY(farfield->find("Attributes") != farfield->end(),
              "Missing \"Attributes\" list for \"FarField\" postprocessing in the "
              "configuration file!");

  attributes = farfield->at("Attributes").get<std::vector<int>>();  // Required
  std::sort(attributes.begin(), attributes.end());

  // Generate NSample points with the following properties:
  // - If NSample >= 2, the generated points are precisely NSample, otherwise NSample = 2.
  // - The poles, the equator, and the XZ plane are always included.
  // - The points are almost uniformly on a sphere, with a small bias due to satisfying the
  //   previous condition.
  // - The points are on rings of constant theta.

  auto nsample_json = farfield->find("NSample");
  int nsample = 0;
  if (nsample_json != farfield->end())
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

  auto thetaphis_json = farfield->find("ThetaPhis");
  if (thetaphis_json != farfield->end())
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
void BoundaryPostData::SetUp(const json &boundaries)
{
  auto postpro = boundaries.find("Postprocessing");
  if (postpro == boundaries.end())
  {
    return;
  }
  flux.SetUp(*postpro);
  dielectric.SetUp(*postpro);
  farfield.SetUp(*postpro);

  // Store all unique postprocessing boundary attributes.
  for (const auto &[idx, data] : flux)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }
  for (const auto &[idx, data] : dielectric)
  {
    attributes.insert(attributes.end(), data.attributes.begin(), data.attributes.end());
  }

  attributes.insert(attributes.end(), farfield.attributes.begin(),
                    farfield.attributes.end());

  std::sort(attributes.begin(), attributes.end());
  attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
  attributes.shrink_to_fit();
}

void BoundaryData::SetUp(const json &config)
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
  lumpedport.SetUp(*boundaries);
  periodic.SetUp(*boundaries);
  waveport.SetUp(*boundaries);
  current.SetUp(*boundaries);
  postpro.SetUp(*boundaries);

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

void DrivenSolverData::SetUp(const json &solver)
{
  auto driven = solver.find("Driven");
  if (driven == solver.end())
  {
    return;
  }

  restart = driven->value("Restart", restart);
  adaptive_tol = driven->value("AdaptiveTol", adaptive_tol);
  adaptive_max_size = driven->value("AdaptiveMaxSamples", adaptive_max_size);
  adaptive_memory = driven->value("AdaptiveConvergenceMemory", adaptive_memory);

  MFEM_VERIFY(!(restart != 1 && adaptive_tol > 0.0),
              "\"Restart\" is incompatible with adaptive frequency sweep!");

  std::vector<double> save_f, prom_f;  // samples to be saved to paraview and added to prom
  // Backwards compatible top level interface.
  if (driven->find("MinFreq") != driven->end() &&
      driven->find("MaxFreq") != driven->end() && driven->find("FreqStep") != driven->end())
  {
    double min_f = driven->at("MinFreq");     // Required
    double max_f = driven->at("MaxFreq");     // Required
    double delta_f = driven->at("FreqStep");  // Required
    sample_f = ConstructLinearRange(min_f, max_f, delta_f);
    if (int save_step = driven->value("SaveStep", 0); save_step > 0)
    {
      for (std::size_t n = 0; n < sample_f.size(); n += save_step)
      {
        save_f.emplace_back(sample_f[n]);
      }
    }
  }
  if (auto freq_samples = driven->find("Samples"); freq_samples != driven->end())
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
  auto explicit_save_f = driven->value("Save", std::vector<double>());
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

void EigenSolverData::SetUp(const json &solver)
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
  nonlinear_type = eigenmode->value("NonlinearType", nonlinear_type);
  refine_nonlinear = eigenmode->value("RefineNonlinear", refine_nonlinear);
  linear_tol = eigenmode->value("LinearTol", linear_tol);
  target_upper = eigenmode->value("TargetUpper", target_upper);
  preconditioner_lag = eigenmode->value("PreconditionerLag", preconditioner_lag);
  preconditioner_lag_tol = eigenmode->value("PreconditionerLagTol", preconditioner_lag_tol);
  max_restart = eigenmode->value("MaxRestart", max_restart);

  target_upper = (target_upper < 0) ? 3 * target : target_upper;  // default = 3 * target
  MFEM_VERIFY(target_upper > target, "config[\"Eigenmode\"][\"TargetUpper\"] must be "
                                     "greater than config[\"Eigenmode\"][\"Target\"]!");
}

void ElectrostaticSolverData::SetUp(const json &solver)
{
  auto electrostatic = solver.find("Electrostatic");
  if (electrostatic == solver.end())
  {
    return;
  }
  n_post = electrostatic->value("Save", n_post);
}

void MagnetostaticSolverData::SetUp(const json &solver)
{
  auto magnetostatic = solver.find("Magnetostatic");
  if (magnetostatic == solver.end())
  {
    return;
  }
  n_post = magnetostatic->value("Save", n_post);
}

void TransientSolverData::SetUp(const json &solver)
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

  if (type == TimeSteppingScheme::GEN_ALPHA || type == TimeSteppingScheme::RUNGE_KUTTA)
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
}

void LinearSolverData::SetUp(const json &solver)
{
  auto linear = solver.find("Linear");
  if (linear == solver.end())
  {
    return;
  }
  type = linear->value("Type", type);
  krylov_solver = linear->value("KSPType", krylov_solver);
  tol = linear->value("Tol", tol);
  max_it = linear->value("MaxIts", max_it);
  MFEM_VERIFY(max_it > 0,
              "config[\"Solver\"][\"Linear\"][\"MaxIts\"] must be strictly positive!");
  max_size = linear->value("MaxSize", max_size);
  initial_guess = linear->value("InitialGuess", initial_guess);

  // Options related to multigrid.
  mg_max_levels = linear->value("MGMaxLevels", mg_max_levels);
  mg_coarsening = linear->value("MGCoarsenType", mg_coarsening);
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
  complex_coarse_solve = linear->value("ComplexCoarseSolve", complex_coarse_solve);
  drop_small_entries = linear->value("DropSmallEntries", drop_small_entries);
  reorder_reuse = linear->value("ReorderingReuse", reorder_reuse);
  pc_side = linear->value("PCSide", pc_side);
  sym_factorization = linear->value("ColumnOrdering", sym_factorization);
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
  gs_orthog = linear->value("GSOrthogonalization", gs_orthog);
}

void SolverData::SetUp(const json &config)
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
