// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/enum_string.hpp"

#include <algorithm>
#include <iterator>
#include <utility>
#include <mfem.hpp>

// The single mapping table per enum lives here. ToString/FromString use linear find_if in
// each direction, matching the previous PALACE_JSON_SERIALIZE_ENUM behavior: where DEFAULT
// is an alias (same underlying value as a concrete entry), it is listed last so ToString
// returns the concrete name (first match) while FromString still accepts "Default".

namespace palace
{

#define PALACE_ENUM_STRING_DEFINE(ENUM_TYPE, ...)                               \
  std::string ToString(ENUM_TYPE e)                                             \
  {                                                                             \
    static const std::pair<ENUM_TYPE, const char *> m[] = __VA_ARGS__;          \
    auto it = std::find_if(std::begin(m), std::end(m),                          \
                           [e](const std::pair<ENUM_TYPE, const char *> &p)     \
                           { return p.first == e; });                           \
    MFEM_VERIFY(it != std::end(m),                                              \
                "Invalid value for " #ENUM_TYPE " when converting to string!"); \
    return it->second;                                                          \
  }                                                                             \
  void FromString(std::string_view s, ENUM_TYPE &e)                             \
  {                                                                             \
    static const std::pair<ENUM_TYPE, const char *> m[] = __VA_ARGS__;          \
    auto it = std::find_if(std::begin(m), std::end(m),                          \
                           [s](const std::pair<ENUM_TYPE, const char *> &p)     \
                           { return s == p.second; });                          \
    MFEM_VERIFY(it != std::end(m), "Invalid value ("                            \
                                       << s << ") for " #ENUM_TYPE              \
                                       << " given in the configuration file!"); \
    e = it->first;                                                              \
  }

PALACE_ENUM_STRING_DEFINE(CoordinateSystem,
                          {{CoordinateSystem::CARTESIAN, "Cartesian"},
                           {CoordinateSystem::CYLINDRICAL, "Cylindrical"}})

PALACE_ENUM_STRING_DEFINE(ProblemType, {{ProblemType::DRIVEN, "Driven"},
                                        {ProblemType::EIGENMODE, "Eigenmode"},
                                        {ProblemType::ELECTROSTATIC, "Electrostatic"},
                                        {ProblemType::MAGNETOSTATIC, "Magnetostatic"},
                                        {ProblemType::TRANSIENT, "Transient"},
                                        {ProblemType::BOUNDARYMODE, "BoundaryMode"}})

// The list order matters: concrete entries must come before DEFAULT (the alias) so
// ToString emits the concrete name ("SLEPc" or "ARPACK") rather than "Default"; FromString
// on either string still maps to the same underlying value as the alias.
PALACE_ENUM_STRING_DEFINE(EigenSolverBackend, {{EigenSolverBackend::SLEPC, "SLEPc"},
                                               {EigenSolverBackend::ARPACK, "ARPACK"},
                                               {EigenSolverBackend::DEFAULT, "Default"}})

PALACE_ENUM_STRING_DEFINE(NonlinearEigenSolver, {{NonlinearEigenSolver::HYBRID, "Hybrid"},
                                                 {NonlinearEigenSolver::SLP, "SLP"}})

PALACE_ENUM_STRING_DEFINE(SurfaceFlux, {{SurfaceFlux::ELECTRIC, "Electric"},
                                        {SurfaceFlux::MAGNETIC, "Magnetic"},
                                        {SurfaceFlux::POWER, "Power"}})

PALACE_ENUM_STRING_DEFINE(InterfaceDielectric, {{InterfaceDielectric::DEFAULT, "Default"},
                                                {InterfaceDielectric::MA, "MA"},
                                                {InterfaceDielectric::MS, "MS"},
                                                {InterfaceDielectric::SA, "SA"}})

PALACE_ENUM_STRING_DEFINE(FrequencySampling, {{FrequencySampling::DEFAULT, "Default"},
                                              {FrequencySampling::LINEAR, "Linear"},
                                              {FrequencySampling::LOG, "Log"},
                                              {FrequencySampling::POINT, "Point"}})

// DEFAULT is an alias for GEN_ALPHA; placing it last ensures ToString emits
// "GeneralizedAlpha" rather than "Default" while still accepting "Default" as input.
PALACE_ENUM_STRING_DEFINE(TimeSteppingScheme,
                          {{TimeSteppingScheme::GEN_ALPHA, "GeneralizedAlpha"},
                           {TimeSteppingScheme::RUNGE_KUTTA, "RungeKutta"},
                           {TimeSteppingScheme::CVODE, "CVODE"},
                           {TimeSteppingScheme::ARKODE, "ARKODE"},
                           {TimeSteppingScheme::DEFAULT, "Default"}})

PALACE_ENUM_STRING_DEFINE(Excitation,
                          {{Excitation::SINUSOIDAL, "Sinusoidal"},
                           {Excitation::GAUSSIAN, "Gaussian"},
                           {Excitation::DIFF_GAUSSIAN, "DifferentiatedGaussian"},
                           {Excitation::MOD_GAUSSIAN, "ModulatedGaussian"},
                           {Excitation::RAMP_STEP, "Ramp"},
                           {Excitation::SMOOTH_STEP, "SmoothStep"}})

PALACE_ENUM_STRING_DEFINE(LinearSolver, {{LinearSolver::DEFAULT, "Default"},
                                         {LinearSolver::AMS, "AMS"},
                                         {LinearSolver::BOOMER_AMG, "BoomerAMG"},
                                         {LinearSolver::MUMPS, "MUMPS"},
                                         {LinearSolver::SUPERLU, "SuperLU"},
                                         {LinearSolver::STRUMPACK, "STRUMPACK"},
                                         {LinearSolver::STRUMPACK_MP, "STRUMPACK-MP"},
                                         {LinearSolver::JACOBI, "Jacobi"}})

PALACE_ENUM_STRING_DEFINE(KrylovSolver, {{KrylovSolver::DEFAULT, "Default"},
                                         {KrylovSolver::CG, "CG"},
                                         {KrylovSolver::MINRES, "MINRES"},
                                         {KrylovSolver::GMRES, "GMRES"},
                                         {KrylovSolver::FGMRES, "FGMRES"},
                                         {KrylovSolver::BICGSTAB, "BiCGSTAB"}})

PALACE_ENUM_STRING_DEFINE(MultigridCoarsening,
                          {{MultigridCoarsening::LINEAR, "Linear"},
                           {MultigridCoarsening::LOGARITHMIC, "Logarithmic"}})

PALACE_ENUM_STRING_DEFINE(PreconditionerSide, {{PreconditionerSide::DEFAULT, "Default"},
                                               {PreconditionerSide::RIGHT, "Right"},
                                               {PreconditionerSide::LEFT, "Left"}})

PALACE_ENUM_STRING_DEFINE(SymbolicFactorization,
                          {{SymbolicFactorization::DEFAULT, "Default"},
                           {SymbolicFactorization::METIS, "METIS"},
                           {SymbolicFactorization::PARMETIS, "ParMETIS"},
                           {SymbolicFactorization::SCOTCH, "Scotch"},
                           {SymbolicFactorization::PTSCOTCH, "PTScotch"},
                           {SymbolicFactorization::PORD, "PORD"},
                           {SymbolicFactorization::AMD, "AMD"},
                           {SymbolicFactorization::RCM, "RCM"}})

PALACE_ENUM_STRING_DEFINE(SparseCompression,
                          {{SparseCompression::NONE, "None"},
                           {SparseCompression::BLR, "BLR"},
                           {SparseCompression::HSS, "HSS"},
                           {SparseCompression::HODLR, "HODLR"},
                           {SparseCompression::ZFP, "ZFP"},
                           {SparseCompression::BLR_HODLR, "BLR-HODLR"},
                           {SparseCompression::ZFP_BLR_HODLR, "ZFP-BLR-HODLR"}})

PALACE_ENUM_STRING_DEFINE(Orthogonalization, {{Orthogonalization::MGS, "MGS"},
                                              {Orthogonalization::CGS, "CGS"},
                                              {Orthogonalization::CGS2, "CGS2"}})

PALACE_ENUM_STRING_DEFINE(DomainOrthogonalizationWeight,
                          {{DomainOrthogonalizationWeight::ENERGY, "Energy"},
                           {DomainOrthogonalizationWeight::FE_BASIS_IDENTITY,
                            "FEBasisIdentity"},
                           {DomainOrthogonalizationWeight::SPACE_OVERLAP, "SpaceOverlap"}})

PALACE_ENUM_STRING_DEFINE(Device, {{Device::CPU, "CPU"},
                                   {Device::GPU, "GPU"},
                                   {Device::DEBUG, "Debug"}})

PALACE_ENUM_STRING_DEFINE(PMLStretchFormulation,
                          {{PMLStretchFormulation::FIXED, "Fixed"},
                           {PMLStretchFormulation::CFS, "CFS"},
                           {PMLStretchFormulation::FREQUENCY_DEPENDENT,
                            "FrequencyDependent"}})

PALACE_ENUM_STRING_DEFINE(PMLCoordinateType,
                          {{PMLCoordinateType::CARTESIAN, "Cartesian"}})

#undef PALACE_ENUM_STRING_DEFINE

}  // namespace palace
