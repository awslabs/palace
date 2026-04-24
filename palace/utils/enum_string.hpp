// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_ENUM_STRING_HPP
#define PALACE_UTILS_ENUM_STRING_HPP

#include <string>
#include <string_view>
#include "labels.hpp"

// Single source of truth for converting the configuration enums in labels.hpp to and from
// their string representations. The bodies (and the one mapping table per enum) live in
// enum_string.cpp; the JSON to_json/from_json shims in configfile.cpp and the resolved
// config writer both delegate to ToString/FromString so there is exactly one mapping.

namespace palace
{

#define PALACE_ENUM_STRING_DECL(ENUM_TYPE) \
  std::string ToString(ENUM_TYPE e);       \
  void FromString(std::string_view s, ENUM_TYPE &e);

PALACE_ENUM_STRING_DECL(CoordinateSystem)
PALACE_ENUM_STRING_DECL(ProblemType)
PALACE_ENUM_STRING_DECL(EigenSolverBackend)
PALACE_ENUM_STRING_DECL(NonlinearEigenSolver)
PALACE_ENUM_STRING_DECL(SurfaceFlux)
PALACE_ENUM_STRING_DECL(InterfaceDielectric)
PALACE_ENUM_STRING_DECL(FrequencySampling)
PALACE_ENUM_STRING_DECL(TimeSteppingScheme)
PALACE_ENUM_STRING_DECL(Excitation)
PALACE_ENUM_STRING_DECL(LinearSolver)
PALACE_ENUM_STRING_DECL(KrylovSolver)
PALACE_ENUM_STRING_DECL(MultigridCoarsening)
PALACE_ENUM_STRING_DECL(PreconditionerSide)
PALACE_ENUM_STRING_DECL(SymbolicFactorization)
PALACE_ENUM_STRING_DECL(SparseCompression)
PALACE_ENUM_STRING_DECL(Orthogonalization)
PALACE_ENUM_STRING_DECL(DomainOrthogonalizationWeight)
PALACE_ENUM_STRING_DECL(Device)
PALACE_ENUM_STRING_DECL(PMLStretchFormulation)
PALACE_ENUM_STRING_DECL(PMLCoordinateType)

#undef PALACE_ENUM_STRING_DECL

}  // namespace palace

#endif  // PALACE_UTILS_ENUM_STRING_HPP
