// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PML_HPP
#define PALACE_MODELS_PML_HPP

#include <array>
#include <complex>
#include <vector>
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace::pml
{

// A compiled, nondimensionalized per-region PML profile. Built from config::PMLData at
// setup time. Holds the geometric extents of the PML layer (inner interface + thickness
// along each signed face) and the profile parameters (order, σ_max, κ_max, α_max).
//
// Geometry convention: direction_signs and thickness are indexed per signed face in the
// order {−x, +x, −y, +y, −z, +z}. For each active face:
//   - sign != 0 ⇒ absorption is active in that direction
//   - thickness[i] > 0 ⇒ PML layer depth
//   - interface_coord[axis][0] = inner interface on the negative-x side (min of PML
//     region becomes outer, max is inner)
//   - interface_coord[axis][1] = inner interface on the positive-x side
//
// The "inner interface" is the PML/physical boundary; the "outer boundary" is the
// truncation face where a PEC (or PMC) condition closes the domain.
struct Profile
{
  // Bulk properties (nondimensional) of the background material filling this PML region.
  // These are multiplied by the Jacobian combination to produce the effective μ̃⁻¹, ε̃.
  double epsilon_r = 1.0;
  double mu_r = 1.0;

  // Grading order n in σ(x) = σ_max · ((x − x₀) / d)^n.
  int order = 3;

  // Per-axis peak σ, κ, α (ready to use; sentinels already resolved).
  std::array<double, 3> sigma_max{{0.0, 0.0, 0.0}};
  std::array<double, 3> kappa_max{{1.0, 1.0, 1.0}};
  std::array<double, 3> alpha_max{{0.0, 0.0, 0.0}};

  // Per-axis PML thickness from the inner interface. One value per axis: we assume the
  // negative and positive faces of the same axis use the same thickness (typical usage).
  // If different thicknesses per face are needed, revisit.
  std::array<double, 3> thickness{{0.0, 0.0, 0.0}};

  // Inner-interface coordinate per signed face. interface_coord[axis][0] is the x
  // coordinate of the inner edge of the negative-side PML on that axis;
  // interface_coord[axis][1] is the inner edge of the positive-side PML. Elements
  // outside these limits are in the PML along that axis.
  std::array<std::array<double, 2>, 3> interface_coord{};

  // Which signed faces are active (1 = active, 0 = inactive), layout identical to
  // config::PMLData::direction_signs.
  std::array<int, 6> direction_active{{0, 0, 0, 0, 0, 0}};

  PMLStretchFormulation formulation = PMLStretchFormulation::FIXED;

  // Reference ω₀ for FIXED and CFS (nondimensional). Must be set (> 0) before use.
  double reference_frequency = 0.0;

  bool allow_refinement = false;
};

// Compute the per-axis depth d_i(x) into the PML layer from the inner interface, along
// each active direction. Returns 0 for axes where the element is not in the PML.
// Inputs x are in the same frame as `profile.interface_coord` (nondimensional).
std::array<double, 3> ComputeDepth(const Profile &profile, const std::array<double, 3> &x);

// Compute the local σ, κ, α at position x using the polynomial grading profile.
// σ_i(x) = σ_max_i · (d_i(x) / thickness_i)^n; similarly for α_i. κ_i(x) goes from 1 at
// the inner interface to κ_max_i at the outer, using the same polynomial shape.
struct LocalStretchParams
{
  std::array<double, 3> sigma{{0.0, 0.0, 0.0}};
  std::array<double, 3> kappa{{1.0, 1.0, 1.0}};
  std::array<double, 3> alpha{{0.0, 0.0, 0.0}};
};
LocalStretchParams ComputeLocalStretchParams(const Profile &profile,
                                             const std::array<double, 3> &x);

// Compute the diagonal stretch s_i(ω) for the three formulations (Palace's e^{+iωt}
// time convention):
//   FIXED / FREQUENCY_DEPENDENT: s_i = 1 − i σ_i/ω
//   CFS                        : s_i = κ_i + σ_i/(α_i + iω)   (ε₀ absorbed into nondim)
// The callers decide which ω to pass in (profile.reference_frequency for FIXED/CFS,
// the live solve frequency for FREQUENCY_DEPENDENT).
std::array<std::complex<double>, 3> ComputeStretch(const LocalStretchParams &local,
                                                   double omega,
                                                   PMLStretchFormulation formulation);

// The complex anisotropic μ̃⁻¹ and ε̃ tensors needed by Palace's real/imag-paired
// integrator path, for a Cartesian (diagonal) PML. Each tensor is stored as the three
// real-part diagonal entries plus the three imaginary-part diagonal entries:
//
//   μ̃⁻¹_ii = mu_inv_re[i] + i · mu_inv_im[i]
//   ε̃_ii   = eps_re[i]    + i · eps_im[i]
struct StretchTensors
{
  std::array<double, 3> mu_inv_re{{1.0, 1.0, 1.0}};
  std::array<double, 3> mu_inv_im{{0.0, 0.0, 0.0}};
  std::array<double, 3> eps_re{{1.0, 1.0, 1.0}};
  std::array<double, 3> eps_im{{0.0, 0.0, 0.0}};
};

// Top-level helper: given a profile, a sample position x, and a solve frequency ω,
// return the four tensors. This is the single source of truth for PML stretch math;
// every caller in the solver pipeline goes through this function.
StretchTensors ComputeStretchTensors(const Profile &profile, const std::array<double, 3> &x,
                                     double omega);

// Resolve "auto" σ_max sentinels (-1) into concrete values using the standard formula
// σ_max_i = −(n + 1) · ln(R) / (2 · d_i · √(εᵣ · μᵣ)). Called when building a Profile
// from a config::PMLData. Leaves positive user-specified values untouched.
void ResolveSigmaMaxDefaults(const config::PMLData &data, double mu_r, double epsilon_r,
                             std::array<double, 3> &sigma_max_out);

// Build a Profile from a parsed config::PMLData. Assumes config data is already
// nondimensionalized. For FREQUENCY_DEPENDENT, reference_frequency is left at zero
// and refilled each time ComputeStretchTensors is called with the live ω.
Profile BuildProfile(const config::PMLData &data, double mu_r, double epsilon_r);

// Representative sample point inside a PML attribute's meshed slab, given the global
// mesh bounding box. For axis-aligned PML layers this is the midpoint of the active
// face's slab along each direction.
std::array<double, 3> ComputeSlabCentroid(const config::PMLData &data,
                                          const std::array<double, 3> &bbmin,
                                          const std::array<double, 3> &bbmax);

// Auto-detect which faces of a PML region are active and how thick the layer is, from
// the per-attribute bounding box and the global mesh bounding box. A face is
// considered "active" if the PML bbox reaches the global bbox on that side (within
// `tol` relative to the global extent) while being strictly inside the global bbox on
// the opposite side. The thickness along each axis is the PML bbox extent in that
// direction. Returns (direction_signs, thickness) in PMLData's convention.
struct SlabGeometry
{
  std::array<int, 6> direction_signs{{0, 0, 0, 0, 0, 0}};
  std::array<double, 6> thickness{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
};
SlabGeometry DetectSlabGeometry(const std::array<double, 3> &attr_min,
                                const std::array<double, 3> &attr_max,
                                const std::array<double, 3> &global_min,
                                const std::array<double, 3> &global_max,
                                double rel_tol = 1.0e-6);

}  // namespace palace::pml

#endif  // PALACE_MODELS_PML_HPP
