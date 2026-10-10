// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PML_HPP
#define PALACE_MODELS_PML_HPP

#include <array>
#include <complex>
#include <vector>
#include <mfem.hpp>
#include "utils/configfile.hpp"

// Forward-declare the union used by libCEED QFunction contexts (defined in
// fem/qfunctions/coeff/coeff_qf.h).
union CeedIntScalar;

namespace palace
{

class Mesh;

namespace pml
{

//
// Cartesian perfectly matched layer (PML) model. Within the layer, the coordinates normal
// to each active face of the PML box are stretched by the complex factor (Palace's e^{+iωt}
// time convention)
//                       s(x, ω) = κ(x) + σ(x) / (α(x) + iω) ,
// with σ, κ - 1, α graded polynomially from zero at the inner interface to their maximum
// values at the outer boundary. The stretch is applied through the equivalent anisotropic
// material tensors μ̃⁻¹ = S μ⁻¹ S / det(S) and ε̃ = det(S) S⁻¹ ε S⁻¹, S = diag(s_x, s_y,
// s_z), of the background material properties μ⁻¹ and ε of each PML region, evaluated at
// each quadrature point (see fem/qfunctions/coeff/pml_qf.h). A static stretch uses a fixed
// real reference frequency ω₀, while a frequency-dependent stretch uses the (possibly
// complex) solve frequency.
//
// The PML regions are configured in blocks, each with its own stretch: the stretch of a
// block is the same function of position in all of its PML regions, so that the PML is a
// coordinate transformation, also across the material interfaces inside of the layer (for
// example, a substrate crossing the PML). The stretch must be continuous across the
// interfaces between the PML regions of different blocks, and equal to one at the
// interfaces with the physical (non-PML) region (see CheckStretchContinuity).
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

// The coordinate stretch of the PML.
struct Stretch
{
  LayerGeometry geometry;
  std::array<double, 6> sigma_max{};
  std::array<double, 3> kappa_max{{1.0, 1.0, 1.0}};
  std::array<double, 3> alpha_max{};
  int order = 3;
  bool frequency_dependent = false;
  double reference_frequency = 0.0;
};

// Background material tensors of a PML region (3 x 3, column-major).
struct Background
{
  std::array<double, 9> mu_inv{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
  std::array<double, 9> epsilon_real{{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
  std::array<double, 9> epsilon_imag{};
};

// Detect the PML layer geometry from the bounding box of the physical (non-PML) region,
// whose faces are the inner PML interfaces, and the bounding box of the PML regions of a
// block, whose faces are the outer PML boundaries. A face is active if the two boxes differ
// on that face by more than rel_tol times the extent of the outer box.
LayerGeometry DetectLayerGeometry(const std::array<double, 3> &inner_min,
                                  const std::array<double, 3> &inner_max,
                                  const std::array<double, 3> &outer_min,
                                  const std::array<double, 3> &outer_max,
                                  double rel_tol = 1.0e-6);

// Layer geometry from explicitly configured directions and thicknesses, measured inward
// from the faces of the bounding box of the PML regions of a block.
LayerGeometry ConfiguredLayerGeometry(const config::PMLData &data,
                                      const std::array<double, 3> &outer_min,
                                      const std::array<double, 3> &outer_max);

// Smallest refractive index of a background material over its principal directions, for
// the 3 x 3 (column-major) tensors μ⁻¹ and Re{ε}.
double RefractiveIndex(const std::array<double, 9> &mu_inv,
                       const std::array<double, 9> &epsilon_real);

// Stretch from a (nondimensionalized) configuration and layer geometry. The peak
// conductivity σ_max of each active face is either the configured value for the face's
// axis, or, if not specified, σ_max = -(n + 1) ln(R) / (2 d n_r) for a target
// normal-incidence reflection coefficient R, layer thickness d, grading order n, and
// refractive index n_r: the smallest refractive index of the background materials of the
// PML regions of the block, so that each of them meets the reflection target.
Stretch BuildStretch(const config::PMLData &data, const LayerGeometry &geometry,
                     double n_r);

// Output part of the complex PML tensor terms assembled by a PML integrator.
enum class TensorPart : int
{
  REAL = 0,
  IMAG = 1,
  ABS = 2  // Re{c} |T| with entrywise magnitude, for real-valued approximations
};

// Per-integrator data of the PML QFunction context: the integrator assembles
// part(c_muinv μ̃⁻¹) and/or part(c_eps ε̃) terms, with the frequency-dependent stretch
// evaluated at omega (ignored for a static stretch).
struct ContextHeader
{
  TensorPart part = TensorPart::REAL;
  std::complex<double> c_muinv = 0.0, c_eps = 0.0, omega = 0.0;
  std::array<double, 9> wave_vector_cross{};  // [k ×], column-major
};

// Stretch factors s_a(x, ω) at the point x, with the frequency-dependent stretch evaluated
// at omega (ignored for a static stretch).
std::array<std::complex<double>, 3> EvaluateStretch(const Stretch &stretch,
                                                    const std::array<double, 3> &x,
                                                    std::complex<double> omega);

// Pack the QFunction context of a PML integrator (layout in fem/qfunctions/coeff/pml_qf.h),
// with attr_background the background index of each (1-based) libCEED attribute (-1 for
// attributes outside of the PML regions). Returns an empty context if no attribute is in a
// PML region.
std::vector<CeedIntScalar> PackContext(const ContextHeader &header, const Stretch &stretch,
                                       const std::vector<int> &attr_background,
                                       const std::vector<Background> &backgrounds);

//
// A block of PML regions of a 3D mesh, for frequency domain problems: the stretch and the
// background material of each PML attribute (the material properties of its material).
//
class Layer
{
private:
  Stretch stretch;

  // Domain attributes of the PML regions (global mesh attributes, sorted).
  std::vector<int> attributes;

  // Background material index of each local libCEED attribute (-1 for attributes outside
  // of the PML regions), and background materials.
  std::vector<int> attr_background;
  std::vector<Background> backgrounds;

public:
  // Set up the PML regions of the block configured by data, with pml_attributes the
  // attributes of the PML regions of all blocks (the remaining attributes form the physical
  // region). The background materials are the materials with the given properties
  // (indexed by material, with attr_mat the map from libCEED attribute to material index).
  // Checks that each element of the PML regions of the block is in the layer.
  Layer(const config::PMLData &data, const std::vector<int> &pml_attributes,
        const std::vector<config::MaterialData> &materials, const Mesh &mesh,
        const mfem::Array<int> &attr_mat, const mfem::DenseTensor &mu_inv,
        const mfem::DenseTensor &epsilon_real, const mfem::DenseTensor &epsilon_imag);

  const Stretch &GetStretch() const { return stretch; }
  bool IsFrequencyDependent() const { return stretch.frequency_dependent; }
  const std::vector<int> &GetAttributes() const { return attributes; }

  // Whether the local libCEED attribute is in a PML region of the block.
  bool IsPMLCeedAttribute(int ceed_attr) const
  {
    return ceed_attr > 0 && ceed_attr <= static_cast<int>(attr_background.size()) &&
           attr_background[ceed_attr - 1] >= 0;
  }

  // Pack the QFunction context of a PML integrator (see pml::PackContext). Returns an empty
  // context if no local attribute is in a PML region.
  std::vector<CeedIntScalar> PackContext(const ContextHeader &header) const;
};

// Check that the stretch is continuous across the interfaces between the PML regions of
// different blocks, and equal to one at the interfaces between the PML regions and the
// physical region: at points of each mesh face on these interfaces, the stretch factors of
// the two sides must be the same functions of the frequency. Otherwise, the
// PML is not a coordinate transformation and the interfaces reflect.
void CheckStretchContinuity(const std::vector<Layer> &layers, const Mesh &mesh);

}  // namespace pml

}  // namespace palace

#endif  // PALACE_MODELS_PML_HPP
