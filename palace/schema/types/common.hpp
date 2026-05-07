// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_COMMON_HPP
#define PALACE_SCHEMA_TYPES_COMMON_HPP

// Shared enum types and type aliases used across every Palace configuration
// section. Enumerator names are the JSON wire format — reflect-cpp emits the
// enumerator name as the JSON value, so the casing here is load-bearing.
// Palace values containing hyphens (e.g. "STRUMPACK-MP", "BLR-HODLR") are not
// legal C++ identifiers and become underscored here; the Phase 2 rename pass
// will map them back to the hyphenated wire value.
//
// Enums that carry PR-716 per-value descriptions are declared via
// `PALACE_SCHEMA_ENUM` — the macro emits the enum body *and* the
// `enum_descriptions<E>` specialization from a single list. Enums without
// descriptions stay as plain `enum class`.

#include <array>
#include <vector>

#include <rfl.hpp>

#include "schema/utils/annotations.hpp"

namespace palace::schema
{

// Shared alias for "list of integer attributes." Palace's hand-written schema
// defines `#/$defs/Attributes` with a `minItems: 1` constraint and $refs it
// from every boundary / material entry. reflect-cpp's Validator<T, Minimum>
// only composes with scalar numeric types, so minItems is enforced at the
// application layer via MFEM_VERIFY (same as Palace does today).
using AttributeList = std::vector<int>;

// Shared alias for the "direction" field used by ports and sources
// (LumpedPort, Element, SurfaceCurrent). PR 716 accepts either an axis
// keyword (Cartesian X/Y/Z or cylindrical R, with optional sign and case)
// or an explicit 3-element numeric vector. We emit that as
// `anyOf: [{string, enum}, {array, 3 numbers}]` via `rfl::Variant`.
//
// `PortDirectionLabel` is the `rfl::Literal` enumerating the 24 allowed
// axis keywords; reflect-cpp renders it inline as `{"type": "string",
// "enum": [...]}`. Exposed as a nested type so in-class initializers can
// construct a PortDirection from a C++ string literal (the variant itself
// is not constructible from `const char*`).
using PortDirectionLabel = rfl::Literal<"R", "X", "Y", "Z",      //
                                        "+R", "+X", "+Y", "+Z",  //
                                        "-R", "-X", "-Y", "-Z",  //
                                        "r", "x", "y", "z",      //
                                        "+r", "+x", "+y", "+z",  //
                                        "-r", "-x", "-y", "-z">;
using PortDirection = rfl::Variant<PortDirectionLabel, std::array<double, 3>>;

// Cartesian-only counterpart used by `CurrentDipole`. PR 716's
// `domains.json` constrains the dipole direction to the 18 Cartesian axis
// keywords (no `R`/`r` cylindrical variants), since a Dirac current
// source has no notion of a cylindrical frame.
using DipoleDirectionLabel = rfl::Literal<"X", "Y", "Z",      //
                                          "+X", "+Y", "+Z",  //
                                          "-X", "-Y", "-Z",  //
                                          "x", "y", "z",      //
                                          "+x", "+y", "+z",  //
                                          "-x", "-y", "-z">;
using DipoleDirection = rfl::Variant<DipoleDirectionLabel, std::array<double, 3>>;

}  // namespace palace::schema

// Per-arm aliases for the `PortDirection` and `DipoleDirection` variants.
// The variant itself stays inline as `anyOf: [...]`, but each arm is
// rewritten to `{"$ref": "#/$defs/<alias>"}` so the emitted JSON Schema
// matches PR 716's hand-authored layout:
//
//   "Direction": { "anyOf": [
//       { "$ref": "#/$defs/PortDirection" },
//       { "$ref": "#/$defs/Vector3" } ] }
//
// PortDirectionLabel / DipoleDirectionLabel are `rfl::Literal<...>` enums
// rendered inline by reflect-cpp; the arm-alias pass hoists their bodies
// into `$defs/PortDirection` and `$defs/DipoleDirection` respectively.
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::PortDirectionLabel>
{
  static constexpr std::string_view value = "PortDirection";
};

template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::DipoleDirectionLabel>
{
  static constexpr std::string_view value = "DipoleDirection";
};

// Hoist every fixed 3-element `double` array into the shared
// `$defs/Vector3` entry. Used by Center / Translation / FloquetWaveVector
// / etc. as plain field types, and by the `Direction` variants as their
// numeric-vector arm.
template <>
struct palace::schema::utils::schema_alias_name<::std::array<double, 3>>
{
  static constexpr std::string_view value = "Vector3";
};

namespace palace::schema
{

// --- Problem section -------------------------------------------------------

PALACE_SCHEMA_ENUM(ProblemType,
                   (Eigenmode, "Perform an undamped or damped eigenfrequency analysis."),
                   (Driven, "Perform a frequency-domain driven simulation."),
                   (Transient, "Perform a time-domain excitation response simulation."),
                   (Electrostatic, "Perform an electrostatic analysis to compute the "
                                   "capacitance matrix for a set of voltage terminals."),
                   (Magnetostatic, "Perform a magnetostatic analysis to compute the "
                                   "inductance matrix for a set of current sources."));

// --- Boundary / shared ------------------------------------------------------

PALACE_SCHEMA_ENUM(CoordinateSystem, (Cartesian, "Standard Cartesian coordinates."),
                   (Cylindrical,
                    "Cylindrical coordinates (enables `+R`/`-R` directions)."));

PALACE_SCHEMA_ENUM(EigenSolverBackend,
                   (Default, "Use the default solver (currently SLEPc Krylov-Schur)."),
                   (SLEPc, "Krylov-Schur solver from SLEPc."),
                   (ARPACK, "ARPACK eigensolver."));

PALACE_SCHEMA_ENUM(NonlinearEigenSolver,
                   (Hybrid,
                    "Hybrid algorithm: solve a polynomial (quadratic) approximation first, "
                    "then refine with a quasi-Newton nonlinear eigensolver."),
                   (SLP, "SLEPc Successive Linear Problem (SLP) nonlinear eigensolver."));

PALACE_SCHEMA_ENUM(SurfaceFluxType, (Electric, "Integrate the electric flux density."),
                   (Magnetic, "Integrate the magnetic flux density."),
                   (Power, "Integrate the Poynting vector (energy flux)."));

PALACE_SCHEMA_ENUM(InterfaceDielectric,
                   (Default, "Use the full electric field evaluated at the boundary."),
                   (MA, "Use metal-air interface boundary conditions."),
                   (MS, "Use metal-substrate interface boundary conditions."),
                   (SA, "Use substrate-air interface boundary conditions."));

// --- Solver section --------------------------------------------------------

PALACE_SCHEMA_ENUM(
    TimeSteppingScheme, (Default, "Use the default `\"GeneralizedAlpha\"` scheme."),
    (GeneralizedAlpha, "Second-order implicit generalized-α method with ``\\rho_\\infty = "
                       "1``. Unconditionally stable."),
    (RungeKutta,
     "Two-stage singly diagonal implicit Runge-Kutta (SDIRK). Second-order, L-stable."),
    (CVODE,
     "SUNDIALS CVODE implicit multistep method with adaptive time-stepping. Requires "
     "SUNDIALS support (see [installation options](../install.md#Configuration-options))."),
    (ARKODE,
     "SUNDIALS ARKode implicit Runge-Kutta with adaptive time-stepping. Requires SUNDIALS "
     "support (see [installation options](../install.md#Configuration-options))."));

PALACE_SCHEMA_ENUM(
    Excitation, (Sinusoidal, "Sinusoidal excitation at a user-specified frequency."),
    (Gaussian, "Gaussian pulse with a user-specified width (defines the bandwidth)."),
    (DifferentiatedGaussian, "Differentiated Gaussian pulse with a user-specified width."),
    (ModulatedGaussian,
     "Modulated Gaussian pulse at a center frequency and width, with no DC component."),
    (Ramp, "Differentiable unit step function to model the ramp up to a DC signal."),
    (SmoothStep, "Smooth many-times differentiable unit step over a specified width."));

// Palace JSON spells `STRUMPACK_MP` as "STRUMPACK-MP" on the wire; that
// rename is a Phase 2 concern. For now, the enumerator name is the wire
// value.
PALACE_SCHEMA_ENUM(
    LinearSolver,
    (Default,
     "Use `\"AMS\"` for curl-curl and time-domain problems; a sparse direct solver (if "
     "available) for frequency domain; `\"BoomerAMG\"` for electrostatics."),
    (AMS, "Hypre's [Auxiliary-space Maxwell Solver "
          "(AMS)](https://hypre.readthedocs.io/en/latest/solvers-ams.html), an algebraic "
          "multigrid (AMG)-based preconditioner for curl-curl operators."),
    (BoomerAMG,
     "Hypre's [BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) "
     "algebraic multigrid solver."),
    (MUMPS, "[MUMPS](http://mumps.enseeiht.fr/) sparse direct solver in real double "
            "precision. Requires MUMPS support (see [installation "
            "options](../install.md#Configuration-options))."),
    (SuperLU, "[SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist) sparse direct "
              "solver in real double precision. For frequency domain problems uses a real "
              "approximation to the complex system matrix. Requires SuperLU_DIST support "
              "(see [installation options](../install.md#Configuration-options))."),
    (STRUMPACK,
     "[STRUMPACK](https://portal.nersc.gov/project/sparse/strumpack) sparse direct solver "
     "in real double precision. Not compatible with magnetostatics (singular curl-curl "
     "operator); use `\"AMS\"` instead. Requires STRUMPACK support (see [installation "
     "options](../install.md#Configuration-options))."),
    (STRUMPACK_MP, ""),
    (Jacobi, "Diagonal Jacobi preconditioner (not recommended in general)."));

// MINRES and BiCGSTAB have no PR-716 description — empty strings opt out.
PALACE_SCHEMA_ENUM(
    KrylovSolver,
    (Default, "Use `\"GMRES\"` for frequency domain problems; `\"CG\"` for real symmetric "
              "positive-definite problems (transient, electrostatic, magnetostatic)."),
    (CG, "Preconditioned conjugate gradient."), (MINRES, ""), (GMRES, "GMRES."),
    (FGMRES, "Flexible GMRES."), (BiCGSTAB, ""));

PALACE_SCHEMA_ENUM(MultigridCoarsening, (Linear, "Linear coarsening."),
                   (Logarithmic, "Logarithmic coarsening."));

PALACE_SCHEMA_ENUM(PreconditionerSide, (Default, "Solver-default side."),
                   (Right, "Right preconditioning."), (Left, "Left preconditioning."));

enum class SymbolicFactorization
{
  Default,
  METIS,
  ParMETIS,
  Scotch,
  PTScotch,
  PORD,
  AMD,
  RCM,
};

// See STRUMPACK_MP: `BLR_HODLR` / `ZFP_BLR_HODLR` are underscored C++
// enumerators for wire values `"BLR-HODLR"` / `"ZFP-BLR-HODLR"`.
enum class SparseCompression
{
  None,
  BLR,
  HSS,
  HODLR,
  ZFP,
  BLR_HODLR,
  ZFP_BLR_HODLR,
};

PALACE_SCHEMA_ENUM(Orthogonalization, (MGS, "Modified Gram-Schmidt."),
                   (CGS, "Classical Gram-Schmidt."),
                   (CGS2, "Two-step classical Gram-Schmidt with reorthogonalization."));

PALACE_SCHEMA_ENUM(Device, (CPU, "Run on CPU."),
                   (GPU,
                    "Run on GPU via CUDA (`MFEM_USE_CUDA=ON`) or HIP (`MFEM_USE_HIP=ON`)."),
                   (Debug, "MFEM debug device, useful for diagnosing GPU-related issues."));

}  // namespace palace::schema

// Hoist every `CoordinateSystem` field into a shared `$defs/CoordinateSystem`
// entry. Runs after `inject_enum_descriptions`, so the canonical body is
// the fully-decorated `oneOf` with per-value `const`/`description` pairs.
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::CoordinateSystem>
{
  static constexpr std::string_view value = "CoordinateSystem";
};

#endif  // PALACE_SCHEMA_TYPES_COMMON_HPP
