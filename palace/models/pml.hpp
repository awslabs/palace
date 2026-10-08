// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PML_HPP
#define PALACE_MODELS_PML_HPP

#include <array>
#include <complex>
#include <string>
#include <vector>
#include "utils/configfile.hpp"

// Forward-declare the union used by libCEED QFunction contexts (defined in
// fem/qfunctions/coeff/coeff_qf.h).
union CeedIntScalar;

namespace palace::pml
{

//
// Cartesian perfectly matched layer (PML) model. Within the layer, the coordinates normal
// to each active face of the PML box are stretched by the complex factor (Palace's e^{+iωt}
// time convention)
//                       s(x, ω) = κ(x) + σ(x) / (α(x) + iω) ,
// with σ, κ - 1, α graded polynomially from zero at the inner interface to their maximum
// values at the outer boundary. The stretch is applied through the equivalent anisotropic
// material tensors μ̃⁻¹ = S μ⁻¹ S / det(S) and ε̃ = det(S) S⁻¹ ε S⁻¹, S = diag(s_x, s_y,
// s_z), evaluated at each quadrature point (see fem/qfunctions/coeff/pml_qf.h). Static
// profiles use a fixed real reference frequency ω₀ in the stretch, while
// frequency-dependent profiles use the (possibly complex) solve frequency.
//
// Faces are indexed {-x, +x, -y, +y, -z, +z}: face f is on axis f / 2 and on the positive
// side if f % 2 == 1. All quantities are nondimensional.
//

// Coordinate of the inner interface (the boundary with the physical region) and layer
// thickness of each face. Faces with zero thickness are inactive.
struct LayerGeometry
{
  std::array<double, 6> inner{};
  std::array<double, 6> thickness{};
};

// A compiled PML profile: the stretch parameters of one PML material together with its
// background material tensors (3 x 3, column-major).
struct Profile
{
  LayerGeometry geometry;
  std::array<double, 6> sigma_max{};
  std::array<double, 3> kappa_max{{1.0, 1.0, 1.0}};
  std::array<double, 3> alpha_max{};
  int order = 3;
  bool frequency_dependent = false;
  double reference_frequency = 0.0;
  std::array<double, 9> mu_inv{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
  std::array<double, 9> epsilon_real{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
  std::array<double, 9> epsilon_imag{};
};

// Detect the PML layer geometry from the bounding box of the physical (non-PML) region,
// whose faces are the inner PML interfaces, and the bounding box of the whole mesh, whose
// faces are the outer PML boundaries. A face is active if the two boxes differ on that face
// by more than rel_tol times the extent of the mesh.
LayerGeometry DetectLayerGeometry(const std::array<double, 3> &inner_min,
                                  const std::array<double, 3> &inner_max,
                                  const std::array<double, 3> &outer_min,
                                  const std::array<double, 3> &outer_max,
                                  double rel_tol = 1.0e-6);

// Layer geometry from explicitly configured directions and thicknesses, measured inward
// from the faces of the bounding box of the whole mesh.
LayerGeometry ConfiguredLayerGeometry(const config::PMLData &data,
                                      const std::array<double, 3> &outer_min,
                                      const std::array<double, 3> &outer_max);

// Smallest refractive index of a background material over its principal directions.
double RefractiveIndex(const std::array<double, 9> &mu_inv,
                       const std::array<double, 9> &epsilon_real);

// Peak conductivity σ_max of each active face, either the configured value for the face's
// axis, or, if not specified, σ_max = -(n + 1) ln(R) / (2 d n_r) for a target
// normal-incidence reflection coefficient R, layer thickness d, grading order n, and
// refractive index n_r of each face. For the stretch to be the same in all materials on
// a face, n_r is the smallest refractive index of the PML materials on that face, so that
// the reflection target is met in each of them.
std::array<double, 6> ResolveSigmaMax(const config::PMLData &data,
                                      const LayerGeometry &geometry,
                                      const std::array<double, 6> &n_r);

// Build a profile from a (nondimensionalized) configuration, layer geometry, refractive
// indices of the faces for the default σ_max, and background material tensors μ⁻¹ and
// ε = ε' + i ε'' (3 x 3, column-major).
Profile BuildProfile(const config::PMLData &data, const LayerGeometry &geometry,
                     const std::array<double, 6> &n_r, const std::array<double, 9> &mu_inv,
                     const std::array<double, 9> &epsilon_real,
                     const std::array<double, 9> &epsilon_imag);

// Check that two profiles define the same stretch on the faces they share (with tolerance
// tol for coordinates). Returns the names of inconsistent configuration parameters, or an
// empty string.
std::string CheckStretchConsistency(const Profile &p, const Profile &q, double tol);

// Fractional depth d / t ∈ [0, 1] into the layer along each axis at the point x (zero for
// axes along which x is not in the layer).
std::array<double, 3> ComputeDepthFraction(const Profile &profile,
                                           const std::array<double, 3> &x);

// Output part of the complex PML tensor terms assembled by a PML integrator.
enum class TensorPart : int
{
  REAL = 0,
  IMAG = 1,
  ABS = 2  // Re{c} |T| with entrywise magnitude, for real-valued approximations
};

// Per-integrator data of the PML QFunction context: the integrator assembles
// part(c_muinv μ̃⁻¹) and/or part(c_eps ε̃) terms, with the stretch of frequency-dependent
// profiles evaluated at omega.
struct ContextHeader
{
  TensorPart part = TensorPart::REAL;
  std::complex<double> c_muinv = 0.0, c_eps = 0.0, omega = 0.0;
  std::array<double, 9> wave_vector_cross{};  // [k ×], column-major
};

// Pack the QFunction context (layout in fem/qfunctions/coeff/pml_qf.h), mapping each
// libCEED attribute to its profile index in attr_to_profile (-1 for non-PML attributes).
// Only profiles whose frequency dependence matches frequency_dependent are included.
// Returns an empty context if no local attribute maps to an included profile.
std::vector<CeedIntScalar> PackContext(const ContextHeader &header,
                                       const std::vector<int> &attr_to_profile,
                                       const std::vector<Profile> &profiles,
                                       bool frequency_dependent);

}  // namespace palace::pml

#endif  // PALACE_MODELS_PML_HPP
