// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SCHEMA_TYPES_BOUNDARIES_HPP
#define PALACE_SCHEMA_TYPES_BOUNDARIES_HPP

// Mirrors palace::config::Boundary and its children. Descriptions match
// PR 716's scripts/schema/config/boundaries.json.
//
// Phase 1 deviations from PR 716:
//
//   - PEC/Ground and PMC/ZeroCharge are exposed as distinct optional sections;
//     Palace merges them at runtime into a single canonical block. The schema
//     `allOf + not` mutual-exclusion between the two pairs is enforced in
//     Phase 1.5's emit-binary fragment (same as the `if`/`then` rules).
//
//   - LumpedPort.Excitation / WavePort.Excitation: PR 716 accepts
//     `boolean | integer>=0` via `oneOf`; reflect-cpp cannot emit this mixed
//     shape from a single C++ type. We model as `int` (Palace stores it as
//     an int internally: 0/false = inactive, >0 = excitation group index).
//
//   - Materials.*.oneOf scalar-or-3-array pattern: see domains.hpp.

#include <array>
#include <string>
#include <vector>

#include <rfl.hpp>

#include "common.hpp"
#include "schema/utils/annotations.hpp"

namespace palace::schema
{

struct Element
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC_REQUIRED(
      Direction,
      "Excitation direction. Axis-aligned Cartesian directions can be "
      "specified using keywords: `\"+X\"`, `\"-X\"`, `\"+Y\"`, `\"-Y\"`, "
      "`\"+Z\"`, `\"-Z\"`. Coaxial directions use `\"+R\"`, `\"-R\"`. "
      "Alternatively, specify a normalized 3-element array, e.g. `[0.0, "
      "1.0, 0.0]`. The coordinate system is determined by "
      "`\"CoordinateSystem\"`.",
      PortDirection) = {};

  PALACE_SCHEMA_DESC(CoordinateSystem,
                     "Coordinate system for this element's `\"Direction\"` vector.",
                     CoordinateSystem) = CoordinateSystem::Cartesian;
};

// --- Dirichlet-like attribute-only blocks ----------------------------------

struct PEC
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

struct PMC
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

struct Ground
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

struct ZeroCharge
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

struct WavePortPEC
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

// --- Absorbing / Conductivity / Impedance ----------------------------------

struct FarField
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC_ADVANCED(
      Order,
      "Specify a first- or second-order approximation for the absorbing "
      "boundary condition. Second-order is only available for frequency "
      "domain driven simulations.",
      palace::schema::utils::Closed<int, 1, 2>) = 1;
};

struct Conductivity
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Conductivity,
                              "Electrical conductivity for this boundary, S/m.",
                              palace::schema::utils::Min<double, 0.0>) = 0.0;

  PALACE_SCHEMA_DESC(Permeability, "Relative permeability for this boundary.",
                     palace::schema::utils::Min<double, 0.0>) = 1.0;

  PALACE_SCHEMA_DESC(Thickness,
                     "Optional conductor thickness in mesh length units. Activates a "
                     "finite-thickness boundary condition for metal.",
                     palace::schema::utils::Min<double, 0.0>) = 0.0;

  PALACE_SCHEMA_DESC(External,
                     "Whether this boundary is on the exterior of the computational "
                     "domain. Relevant for the thickness correction.",
                     bool) = false;
};

struct Impedance
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC(Rs, "Surface resistance for this impedance boundary, Ω/sq.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Ls, "Surface inductance for this impedance boundary, H/sq.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Cs, "Surface capacitance for this impedance boundary, F/sq.",
                     double) = 0.0;
};

// --- LumpedPort / WavePort / Terminal / SurfaceCurrent ---------------------

// Common circuit/surface parameters and excitation knobs shared by the
// `Attributes` and `Elements` forms of `LumpedPort`. Folded inline at
// each variant arm via inheritance — reflect-cpp aggregates the base
// fields into the derived struct's emitted properties.
struct LumpedPortCommon
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index,
      "Index of this lumped port, used in postprocessing output files. "
      "Must be unique across all port and source types.",
      palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(R,
                     "Circuit resistance, Ω. Use with `\"L\"` and `\"C\"`; do not mix "
                     "with surface parameters `\"Rs\"`, `\"Ls\"`, `\"Cs\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(L,
                     "Circuit inductance, H. Use with `\"R\"` and `\"C\"`; do not mix "
                     "with surface parameters `\"Rs\"`, `\"Ls\"`, `\"Cs\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(C,
                     "Circuit capacitance, F. Use with `\"R\"` and `\"L\"`; do not mix "
                     "with surface parameters `\"Rs\"`, `\"Ls\"`, `\"Cs\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Rs,
                     "Surface resistance, Ω/sq. Use with `\"Ls\"` and `\"Cs\"`; do not "
                     "mix with circuit parameters `\"R\"`, `\"L\"`, `\"C\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Ls,
                     "Surface inductance, H/sq. Use with `\"Rs\"` and `\"Cs\"`; do not "
                     "mix with circuit parameters `\"R\"`, `\"L\"`, `\"C\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Cs,
                     "Surface capacitance, F/sq. Use with `\"Rs\"` and `\"Ls\"`; do not "
                     "mix with circuit parameters `\"R\"`, `\"L\"`, `\"C\"`.",
                     double) = 0.0;

  PALACE_SCHEMA_DESC(Excitation,
                     "Turns on or off port excitation for driven or transient "
                     "simulations. Can be specified as a boolean or as a non-negative "
                     "integer (excitation group index). See the [boundary conditions "
                     "guide](../guide/boundaries.md#Lumped-and-wave-port-excitation) "
                     "for details.",
                     palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC(Active,
                     "Turns on or off the damping boundary condition for this port for "
                     "driven or transient simulations.",
                     bool) = true;
};

// Single-element form: a port covers a contiguous attribute set with a
// single direction. The `LumpedPortCommon` fields are inlined via
// `rfl::Flatten` so the emitted schema lists `Index`/`R`/`L`/... at the
// top level rather than nesting them under a `Common` key. Required
// keys: `Index` (from the common block) and `Attributes`.
struct LumpedPortAttributes
{
  rfl::Flatten<LumpedPortCommon> common = {};

  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes,
      "Integer array of mesh boundary attributes for this lumped port "
      "boundary.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC(Direction,
                     "Excitation direction keyword or 3-array (see PortDirection schema).",
                     PortDirection) = PortDirectionLabel("+X");

  PALACE_SCHEMA_DESC(CoordinateSystem,
                     "Coordinate system used to express the `\"Direction\"` vector. If "
                     "a keyword argument is used for `\"Direction\"` this value is "
                     "ignored.",
                     CoordinateSystem) = CoordinateSystem::Cartesian;
};

// Multi-element form: the port spans multiple disjoint surfaces, each
// with its own direction. Required keys: `Index` (from the common
// block) and `Elements`.
struct LumpedPortElements
{
  rfl::Flatten<LumpedPortCommon> common = {};

  PALACE_SCHEMA_DESC_REQUIRED(
      Elements,
      "Sub-elements for a multielement lumped port. Each element provides "
      "its own attributes / direction / coordinate system. Elements add "
      "in parallel.",
      std::vector<Element>) = {};
};

using LumpedPort = rfl::Variant<LumpedPortAttributes, LumpedPortElements>;

struct Terminal
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index,
      "Index of this terminal, used in postprocessing output files and to "
      "index the computed capacitance matrix.",
      palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};
};

struct WavePort
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index,
      "Index of this wave port, used in postprocessing output files. Must "
      "be unique across all port and source types.",
      palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC(Mode,
                     "Mode index (1-based) for the characteristic port mode of this wave "
                     "port, ranked in order of decreasing wave number.",
                     palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(Offset,
                     "Offset distance used for S-parameter de-embedding for this wave "
                     "port, specified in mesh length units.",
                     palace::schema::utils::Min<double, 0.0>) = 0.0;

  PALACE_SCHEMA_DESC(SolverType,
                     "Eigenvalue solver for computing the boundary mode. Accepts the "
                     "same options as [`/Solver/Eigenmode/Type`](@ref "
                     "config-solver-eigenmode-type).",
                     EigenSolverBackend) = EigenSolverBackend::Default;

  PALACE_SCHEMA_DESC(Excitation,
                     "Turns on or off port excitation for driven simulations. Can be "
                     "specified as a boolean or as a non-negative integer (excitation "
                     "group index). See the [boundary conditions "
                     "guide](../guide/boundaries.md#Lumped-and-wave-port-excitation) "
                     "for details.",
                     palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC(Active,
                     "Turns on or off the damping boundary condition for this port for "
                     "driven simulations.",
                     bool) = true;

  PALACE_SCHEMA_DESC_ADVANCED(
      MaxIts,
      "Maximum number of iterations for the GMRES solver used in the "
      "wave port boundary mode analysis.",
      palace::schema::utils::XMin<int, 0>) = 45;

  PALACE_SCHEMA_DESC_ADVANCED(
      KSPTol,
      "Tolerance for the linear solver used in the wave port boundary "
      "mode analysis.",
      palace::schema::utils::XMin<double, 0.0>) = 1.0e-8;

  PALACE_SCHEMA_DESC_ADVANCED(EigenTol,
                              "Tolerance for the eigenvalue solver used in the wave port "
                              "boundary mode analysis.",
                              palace::schema::utils::XMin<double, 0.0>) = 1.0e-6;

  PALACE_SCHEMA_DESC_ADVANCED(Verbose,
                              "Verbosity level for the wave port linear and eigensolvers.",
                              palace::schema::utils::Min<int, 0>) = 0;
};

// Surface current sources mirror lumped ports' single-vs-multi-element
// shape: the variant has one arm with a flat `Attributes` list and
// another with an `Elements` array.
struct SurfaceCurrentCommon
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Index,
      "Index of this surface current source, used in postprocessing "
      "output files. Must be unique across all port and source types.",
      palace::schema::utils::XMin<int, 0>) = 1;
};

struct SurfaceCurrentAttributes
{
  rfl::Flatten<SurfaceCurrentCommon> common = {};

  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes,
      "Integer array of mesh boundary attributes for this surface current "
      "boundary.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC(Direction,
                     "Excitation direction keyword or 3-array (see PortDirection schema).",
                     PortDirection) = PortDirectionLabel("+X");

  PALACE_SCHEMA_DESC(CoordinateSystem,
                     "Coordinate system for the `\"Direction\"` vector. Same options as "
                     "[`/LumpedPort/CoordinateSystem`](@ref "
                     "config-boundaries-lumpedport-coordinatesystem).",
                     CoordinateSystem) = CoordinateSystem::Cartesian;
};

struct SurfaceCurrentElements
{
  rfl::Flatten<SurfaceCurrentCommon> common = {};

  PALACE_SCHEMA_DESC_REQUIRED(
      Elements,
      "Sub-elements for a multielement surface current source. Each "
      "element provides its own attributes / direction / coordinate "
      "system. Elements add in parallel to give the same total current "
      "as a single-element source.",
      std::vector<Element>) = {};
};

using SurfaceCurrent = rfl::Variant<SurfaceCurrentAttributes, SurfaceCurrentElements>;

// --- Periodic --------------------------------------------------------------

struct BoundaryPair
{
  PALACE_SCHEMA_DESC_REQUIRED(
      DonorAttributes,
      "Integer array of donor mesh boundary attributes for a periodic "
      "boundary pair.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC_REQUIRED(ReceiverAttributes,
                              "Integer array of receiver mesh boundary attributes for a "
                              "periodic boundary pair.",
                              AttributeList) = {};

  PALACE_SCHEMA_DESC(Translation,
                     "3-element translation vector `[dx, dy, dz]` from donor to "
                     "receiver boundary, in mesh length units. If neither "
                     "`\"Translation\"` nor `\"AffineTransformation\"` are specified, "
                     "the transformation is detected automatically.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC(AffineTransformation,
                     "16-element row-major 4×4 affine transformation matrix from donor "
                     "to receiver boundary, in mesh length units. If neither "
                     "`\"Translation\"` nor `\"AffineTransformation\"` are specified, "
                     "the transformation is detected automatically.",
                     std::array<double, 16>) = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
  };
};

struct Periodic
{
  PALACE_SCHEMA_DESC(FloquetWaveVector,
                     "3-element Floquet wave vector `[kx, ky, kz]` defining the phase "
                     "delay between periodic boundaries, in radians per mesh length "
                     "unit.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };

  PALACE_SCHEMA_DESC_REQUIRED(
      BoundaryPairs,
      "Array of donor–receiver boundary pairs defining the periodic "
      "mapping.",
      std::vector<BoundaryPair>) = {};
};

// --- Boundary postprocessing -----------------------------------------------

struct SurfaceFlux
{
  PALACE_SCHEMA_DESC_REQUIRED(Index,
                              "Index of this surface flux postprocessing boundary, used in "
                              "output files.",
                              palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC_REQUIRED(Type, "Type of surface flux to integrate over the boundary.",
                              SurfaceFluxType) = SurfaceFluxType::Electric;

  PALACE_SCHEMA_DESC(TwoSided,
                     "For internal boundary surfaces: when `false`, the flux on both "
                     "sides is averaged; when `true`, it is summed with opposite "
                     "normal direction.",
                     bool) = false;

  PALACE_SCHEMA_DESC(Center,
                     "Point used to determine the outward normal orientation, in mesh "
                     "length units. Only used when `\"TwoSided\"` is `false`. If not "
                     "specified, the point will be computed as the centroid of the "
                     "axis-aligned bounding box for all elements making up the "
                     "postprocessing boundary.",
                     std::array<double, 3>) = {
    {0.0, 0.0, 0.0}
  };
};

struct Dielectric
{
  PALACE_SCHEMA_DESC(Index, "Index of this dielectric interface, used in output files.",
                     palace::schema::utils::XMin<int, 0>) = 1;

  PALACE_SCHEMA_DESC(Attributes,
                     "Integer array of mesh boundary attributes this object applies to.",
                     AttributeList) = {};

  PALACE_SCHEMA_DESC(Type,
                     "Interface type used to determine the boundary conditions for "
                     "computing the EPR. See also [theory "
                     "reference](../reference.md#Bulk-and-interface-dielectric-loss).",
                     InterfaceDielectric) = InterfaceDielectric::Default;

  PALACE_SCHEMA_DESC_REQUIRED(
      Thickness, "Thickness of this dielectric interface, in mesh length units.",
      palace::schema::utils::Min<double, 0.0>) = 0.0;

  PALACE_SCHEMA_DESC_REQUIRED(
      Permittivity,
      "Relative permittivity of this dielectric interface layer. This "
      "should be the interface layer permittivity for the specific "
      "\"Type\" of interface specified.",
      palace::schema::utils::Min<double, 0.0>) = 0.0;

  PALACE_SCHEMA_DESC(LossTan, "Loss tangent of this dielectric interface.",
                     palace::schema::utils::Min<double, 0.0>) = 0.0;
};

struct FarFieldPostprocessing
{
  PALACE_SCHEMA_DESC_REQUIRED(
      Attributes, "Integer array of mesh boundary attributes this object applies to.",
      AttributeList) = {};

  PALACE_SCHEMA_DESC(NSample,
                     "Number of uniformly-spaced points used to discretize the "
                     "far-field sphere.",
                     palace::schema::utils::Min<int, 0>) = 0;

  PALACE_SCHEMA_DESC(ThetaPhis,
                     "Additional specific (θ, φ) angle pairs in degrees at which to "
                     "evaluate the far field. θ ∈ [0°, 180°] is the polar angle, "
                     "φ ∈ [0°, 360°] is the azimuthal angle.",
                     std::vector<std::array<double, 2>>) = {};
};

struct BoundaryPostprocessing
{
  PALACE_SCHEMA_DESC(SurfaceFlux,
                     "Array of surface flux postprocessing boundaries. Results are "
                     "written to `surface-F.csv` in the output directory.",
                     std::vector<SurfaceFlux>) = {};

  PALACE_SCHEMA_DESC(Dielectric,
                     "Array of interface dielectric loss postprocessing boundaries. "
                     "Computes energy participation ratios (EPR) and quality factors "
                     "for dielectric interfaces. See also the [reference "
                     "documentation](../reference.md#Bulk-and-interface-dielectric-loss"
                     "). Results are written to `surface-Q.csv` in the output "
                     "directory.",
                     std::vector<Dielectric>) = {};

  PALACE_SCHEMA_DESC(FarField,
                     "Far-field electric field extraction. The boundary attributes "
                     "must enclose the system and be on an external boundary.",
                     FarFieldPostprocessing) = {};
};

// --- Top-level Boundary ------------------------------------------------

struct Boundaries
{
  PALACE_SCHEMA_DESC(PEC,
                     "Perfect electric conductor (PEC) boundary condition: enforces "
                     "zero tangential electric field. This is a homogeneous Dirichlet "
                     "condition for frequency/time domain and magnetostatic "
                     "formulations.",
                     PEC) = {};

  PALACE_SCHEMA_DESC(PMC,
                     "Perfect magnetic conductor (PMC) boundary condition: enforces "
                     "zero tangential magnetic field. This is the natural "
                     "(homogeneous Neumann) boundary condition; it also imposes "
                     "symmetry of the electric field across the surface.",
                     PMC) = {};

  PALACE_SCHEMA_DESC(Ground,
                     "Zero-voltage (ground) boundary condition for electrostatic "
                     "simulations. Mutually exclusive with [PEC](@ref "
                     "config-boundaries-pec).",
                     Ground) = {};

  PALACE_SCHEMA_DESC(ZeroCharge,
                     "Zero surface charge (homogeneous Neumann) boundary condition "
                     "for electrostatic simulations. Also imposes symmetry of the "
                     "electric field across the surface. Mutually exclusive with "
                     "[PMC](@ref config-boundaries-pmc).",
                     ZeroCharge) = {};

  PALACE_SCHEMA_DESC(WavePortPEC,
                     "Additional PEC boundary conditions for the 2D eigensolve used "
                     "in wave port mode analysis, along with those already specified "
                     "under [PEC](@ref config-boundaries-pec) and [Conductivity](@ref "
                     "config-boundaries-conductivity). Only relevant when [WavePort]"
                     "(@ref config-boundaries-waveport) boundaries are present.",
                     WavePortPEC) = {};

  PALACE_SCHEMA_DESC(Absorbing,
                     "Farfield absorbing (scattering) boundary conditions. These are "
                     "artificial boundary conditions applied at farfield boundaries "
                     "to minimize reflections.",
                     FarField) = {};

  PALACE_SCHEMA_DESC(Conductivity,
                     "Array of finite conductivity surface impedance boundaries. "
                     "Models the effect of a boundary with non-infinite conductivity "
                     "for conductors with thickness much larger than the skin depth. "
                     "Only available for frequency domain driven and eigenmode "
                     "simulations.",
                     std::vector<Conductivity>) = {};

  PALACE_SCHEMA_DESC(Impedance,
                     "Array of surface impedance boundary conditions. The surface "
                     "impedance relates the tangential electric and magnetic fields "
                     "using the parallel combination of the specified resistance, "
                     "inductance, and capacitance per square.",
                     std::vector<Impedance>) = {};

  PALACE_SCHEMA_DESC(LumpedPort,
                     "Array of lumped port boundary conditions. Lumped ports can be "
                     "specified on boundaries internal to the computational domain.",
                     std::vector<LumpedPort>) = {};

  PALACE_SCHEMA_DESC(Terminal,
                     "Array of terminal boundaries for electrostatic simulations. "
                     "Capacitance matrix entries are extracted for each terminal.",
                     std::vector<Terminal>) = {};

  PALACE_SCHEMA_DESC(WavePort,
                     "Array of numeric wave port boundary conditions. Wave ports can "
                     "only be specified on the true boundary of the computational "
                     "domain (they must be \"one-sided\"). A 2D boundary mode "
                     "eigenproblem is solved on each wave port to compute the port "
                     "mode shape. Only available for frequency domain driven and "
                     "eigenmode simulations.",
                     std::vector<WavePort>) = {};

  PALACE_SCHEMA_DESC(SurfaceCurrent,
                     "Array of surface current source boundaries. Prescribes a unit "
                     "source surface current excitation on the given boundary to "
                     "excite a driven, transient, or magnetostatic simulation. For "
                     "magnetostatic simulations, inductance matrix entries are "
                     "extracted for each surface current boundary.",
                     std::vector<SurfaceCurrent>) = {};

  PALACE_SCHEMA_DESC(Periodic,
                     "Periodic boundary conditions for surfaces whose meshes are "
                     "identical after translation and/or rotation. Floquet periodic "
                     "boundary conditions with a phase shift are also supported.",
                     Periodic) = {};

  PALACE_SCHEMA_DESC(Postprocessing, "Configuration for boundary postprocessing.",
                     BoundaryPostprocessing) = {};
};

}  // namespace palace::schema

// --- Variant arm aliases for LumpedPort / SurfaceCurrent ------------------
//
// Each array element of `LumpedPort` / `SurfaceCurrent` is a variant
// over the `Attributes` and `Elements` shapes. The variant is emitted
// inline as `oneOf: [{<arm0 body>}, {<arm1 body>}]`; these alias
// specializations promote each arm body into its own `$defs` entry and
// rewrite the inline arm to a `$ref`, matching PR-716's hand-authored
// layout.
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::LumpedPortAttributes>
{
  static constexpr std::string_view value = "LumpedPortAttributes";
};
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::LumpedPortElements>
{
  static constexpr std::string_view value = "LumpedPortElements";
};
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::SurfaceCurrentAttributes>
{
  static constexpr std::string_view value = "SurfaceCurrentAttributes";
};
template <>
struct palace::schema::utils::schema_alias_name<::palace::schema::SurfaceCurrentElements>
{
  static constexpr std::string_view value = "SurfaceCurrentElements";
};

// Force the variant arm composition keyword to `oneOf`. reflect-cpp's
// default for `rfl::Variant` is `anyOf`; PR-716 expresses these
// mutually-exclusive shapes via `oneOf` (exactly one arm matches).
template <>
struct palace::schema::utils::schema_composition<::palace::schema::LumpedPort>
{
  static constexpr auto value = palace::schema::utils::Compose::OneOf;
};
template <>
struct palace::schema::utils::schema_composition<::palace::schema::SurfaceCurrent>
{
  static constexpr auto value = palace::schema::utils::Compose::OneOf;
};

#endif  // PALACE_SCHEMA_TYPES_BOUNDARIES_HPP
