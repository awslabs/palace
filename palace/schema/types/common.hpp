// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_COMMON_HPP
#define PALACE_SCHEMA_TYPES_COMMON_HPP

// Shared enum types and type aliases used across every Palace configuration
// section. Enumerator names are the JSON wire format — reflect-cpp emits the
// enumerator name as the JSON value, so the casing here is load-bearing.
// Palace values containing hyphens (e.g. "STRUMPACK-MP", "BLR-HODLR") are not
// legal C++ identifiers and become underscored here; the Phase 1.5 enum-
// description pass maps them back to the hyphenated wire value.

#include <vector>

#include <rfl.hpp>

#include "schema/utils/annotations.hpp"

namespace palace::schema {

// Shared alias for "list of integer attributes." Palace's hand-written schema
// defines `#/$defs/Attributes` with a `minItems: 1` constraint and $refs it
// from every boundary / material entry. reflect-cpp's Validator<T, Minimum>
// only composes with scalar numeric types, so minItems is enforced at the
// application layer via MFEM_VERIFY (same as Palace does today).
using AttributeList = std::vector<int>;

// --- Problem section -------------------------------------------------------

enum class ProblemType {
    Eigenmode,
    Driven,
    Transient,
    Electrostatic,
    Magnetostatic,
};

// --- Boundary / shared ------------------------------------------------------

enum class CoordinateSystem {
    Cartesian,
    Cylindrical,
};

enum class EigenSolverBackend {
    Default,
    SLEPc,
    ARPACK,
};

enum class NonlinearEigenSolver {
    Hybrid,
    SLP,
};

enum class SurfaceFlux {
    Electric,
    Magnetic,
    Power,
};

enum class InterfaceDielectric {
    Default,
    MA,
    MS,
    SA,
};

// --- Solver section --------------------------------------------------------

enum class TimeSteppingScheme {
    Default,
    GeneralizedAlpha,
    RungeKutta,
    CVODE,
    ARKODE,
};

enum class Excitation {
    Sinusoidal,
    Gaussian,
    DifferentiatedGaussian,
    ModulatedGaussian,
    Ramp,
    SmoothStep,
};

enum class LinearSolver {
    Default,
    AMS,
    BoomerAMG,
    MUMPS,
    SuperLU,
    STRUMPACK,
    // Palace JSON spells this "STRUMPACK-MP"; hyphen is not a legal C++
    // identifier character. Phase 1.5 remaps the underscored form to the
    // hyphenated wire value during schema emission.
    STRUMPACK_MP,
    Jacobi,
};

enum class KrylovSolver {
    Default,
    CG,
    MINRES,
    GMRES,
    FGMRES,
    BiCGSTAB,
};

enum class MultigridCoarsening {
    Linear,
    Logarithmic,
};

enum class PreconditionerSide {
    Default,
    Right,
    Left,
};

enum class SymbolicFactorization {
    Default,
    METIS,
    ParMETIS,
    Scotch,
    PTScotch,
    PORD,
    AMD,
    RCM,
};

enum class SparseCompression {
    None,
    BLR,
    HSS,
    HODLR,
    ZFP,
    // See STRUMPACK_MP above — hyphen → underscore.
    BLR_HODLR,
    ZFP_BLR_HODLR,
};

enum class Orthogonalization {
    MGS,
    CGS,
    CGS2,
};

enum class Device {
    CPU,
    GPU,
    Debug,
};

}  // namespace palace::schema

#endif  // PALACE_SCHEMA_TYPES_COMMON_HPP
