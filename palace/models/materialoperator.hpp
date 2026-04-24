// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <array>
#include <vector>
#include <mfem.hpp>
#include "fem/mesh.hpp"
#include "models/pml.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class IoData;

//
// A class handling material attributes.
//
class MaterialOperator
{
private:
  // Reference to underlying mesh object (not owned).
  const Mesh &mesh;

  // Mapping from the local libCEED attribute to material index.
  mfem::Array<int> attr_mat;

  // Material properties: relative permeability, relative permittivity, and others (like
  // electrical conductivity, London penetration depth for superconductors and Floquet wave
  // vector).
  mfem::DenseTensor mat_muinv, mat_epsilon, mat_epsilon_imag, mat_epsilon_abs, mat_invz0,
      mat_c0, mat_sigma, mat_invLondon, mat_kxTmuinv, mat_muinvkx, mat_kxTmuinvkx, mat_kx;
  mfem::DenseMatrix wave_vector_cross;
  mfem::Array<double> mat_c0_min, mat_c0_max;

  // PML material tensors. For each libCEED-attribute-indexed position (matching attr_mat
  // size), these store the complex anisotropic μ̃⁻¹ and ε̃ split into real and imaginary
  // parts. Entries are zero for non-PML attributes so they can be assembled
  // unconditionally; SpaceOperator gates with HasPML().
  //
  // Two parallel tensor sets:
  //   *_static_{re,im} — populated only for FIXED and CFS PML attributes. Baked into
  //     K and M at setup. Static during a solve (constant across frequencies).
  //   *_freq_{re,im}   — populated only for FREQUENCY_DEPENDENT PML attributes. Refilled
  //     per solve frequency by RebuildPMLTensors and added into GetExtraSystemMatrix(ω).
  mfem::DenseTensor mat_muinv_pml_static_re, mat_muinv_pml_static_im,
      mat_epsilon_pml_static_re, mat_epsilon_pml_static_im;
  mfem::DenseTensor mat_muinv_pml_freq_re, mat_muinv_pml_freq_im, mat_epsilon_pml_freq_re,
      mat_epsilon_pml_freq_im;

  // Per-PML-attribute profile (one per user-declared PML material block). Indexed by
  // slot in pml_profiles; pml_attr_to_profile maps libCEED attribute → slot, or -1 for
  // non-PML attributes.
  std::vector<pml::Profile> pml_profiles;
  mfem::Array<int> pml_attr_to_profile;

  // Centroid of each PML attribute's meshed region, in nondimensional coordinates. Used
  // as the sample point for ComputeStretchTensors. One entry per pml_profiles slot.
  std::vector<std::array<double, 3>> pml_centroid;

  bool has_pml_attr = false;
  bool has_pml_freq_dependent_attr = false;

  // Are materials isotropic? True when all the material properties are effectively
  // scalar-valued (ie, true scalars or vectors with identical entries). Also true when a
  // material is isotropic, the intersection is true when all are isotropic.
  mfem::Array<bool> attr_is_isotropic;

  // Flag for global domain attributes with nonzero loss tangent, electrical conductivity,
  // London penetration depth, or Floquet wave vector.
  bool has_losstan_attr, has_conductivity_attr, has_london_attr, has_wave_attr;

  void SetUpMaterialProperties(const std::vector<config::MaterialData> &materials,
                               const config::PeriodicBoundaryData &periodic,
                               ProblemType problem_type, const mfem::ParMesh &mesh);
  void SetUpFloquetWaveVector(const config::PeriodicBoundaryData &periodic,
                              ProblemType problem_type, const mfem::ParMesh &mesh);

  // Map from an attribute (specified on a mesh) to a material index (location in the
  // property vector).
  auto AttrToMat(int attr) const
  {
    const auto &loc_attr = mesh.GetCeedAttributes();
    MFEM_ASSERT(loc_attr.find(attr) != loc_attr.end(),
                "Missing libCEED domain attribute for attribute " << attr << "!");
    return attr_mat[loc_attr.at(attr) - 1];
  }

  auto Wrap(const mfem::DenseTensor &data, int attr) const
  {
    const int k = AttrToMat(attr);
    return mfem::DenseMatrix(const_cast<double *>(data.GetData(k)), data.SizeI(),
                             data.SizeJ());
  }

public:
  MaterialOperator(const std::vector<config::MaterialData> &materials,
                   const config::PeriodicBoundaryData &periodic, ProblemType problem_type,
                   const Mesh &mesh);
  MaterialOperator(const IoData &iodata, const Mesh &mesh);

  int SpaceDimension() const { return mat_muinv.SizeI(); }

  auto GetInvPermeability(int attr) const { return Wrap(mat_muinv, attr); }
  auto GetPermittivityReal(int attr) const { return Wrap(mat_epsilon, attr); }
  auto GetPermittivityImag(int attr) const { return Wrap(mat_epsilon_imag, attr); }
  auto GetPermittivityAbs(int attr) const { return Wrap(mat_epsilon_abs, attr); }
  auto GetInvImpedance(int attr) const { return Wrap(mat_invz0, attr); }
  auto GetLightSpeed(int attr) const { return Wrap(mat_c0, attr); }
  auto GetConductivity(int attr) const { return Wrap(mat_sigma, attr); }
  auto GetInvLondonDepth(int attr) const { return Wrap(mat_invLondon, attr); }
  auto GetFloquetCurl(int attr) const { return Wrap(mat_muinvkx, attr); }
  auto GetFloquetMass(int attr) const { return Wrap(mat_kxTmuinvkx, attr); }
  auto GetFloquetCross(int attr) const { return Wrap(mat_kx, attr); }

  auto GetLightSpeedMin(int attr) const { return mat_c0_min[AttrToMat(attr)]; }
  auto GetLightSpeedMax(int attr) const { return mat_c0_max[AttrToMat(attr)]; }

  bool IsIsotropic(int attr) const { return attr_is_isotropic[AttrToMat(attr)]; }

  const auto &GetInvPermeability() const { return mat_muinv; }
  const auto &GetPermittivityReal() const { return mat_epsilon; }
  const auto &GetPermittivityImag() const { return mat_epsilon_imag; }
  const auto &GetPermittivityAbs() const { return mat_epsilon_abs; }
  const auto &GetInvImpedance() const { return mat_invz0; }
  const auto &GetLightSpeed() const { return mat_c0; }
  const auto &GetConductivity() const { return mat_sigma; }
  const auto &GetInvLondonDepth() const { return mat_invLondon; }
  const auto &GetFloquetCurl() const { return mat_muinvkx; }
  const auto &GetFloquetMass() const { return mat_kxTmuinvkx; }
  const auto &GetFloquetCross() const { return mat_kx; }

  const auto &GetLightSpeedMin() const { return mat_c0_min; }
  const auto &GetLightSpeedMax() const { return mat_c0_max; }

  bool HasLossTangent() const { return has_losstan_attr; }
  bool HasConductivity() const { return has_conductivity_attr; }
  bool HasLondonDepth() const { return has_london_attr; }
  bool HasWaveVector() const { return has_wave_attr; }
  bool HasPML() const { return has_pml_attr; }
  bool HasFrequencyDependentPML() const { return has_pml_freq_dependent_attr; }

  // PML tensor accessors (per libCEED attribute). Only meaningful when HasPML() is true.
  //
  // "Static" accessors hold contributions from FIXED/CFS PML attributes — added once to
  // K and M at setup. "Freq" accessors hold contributions from FREQUENCY_DEPENDENT PML
  // attributes — added to the extra system matrix A2(ω) at each solve frequency.
  const auto &GetInvPermeabilityPMLStaticReal() const { return mat_muinv_pml_static_re; }
  const auto &GetInvPermeabilityPMLStaticImag() const { return mat_muinv_pml_static_im; }
  const auto &GetPermittivityPMLStaticReal() const { return mat_epsilon_pml_static_re; }
  const auto &GetPermittivityPMLStaticImag() const { return mat_epsilon_pml_static_im; }
  const auto &GetInvPermeabilityPMLFreqReal() const { return mat_muinv_pml_freq_re; }
  const auto &GetInvPermeabilityPMLFreqImag() const { return mat_muinv_pml_freq_im; }
  const auto &GetPermittivityPMLFreqReal() const { return mat_epsilon_pml_freq_re; }
  const auto &GetPermittivityPMLFreqImag() const { return mat_epsilon_pml_freq_im; }

  // Populate the eight PML tensors from pml_profiles. Two modes:
  //
  //   RebuildPMLTensors(ω, /*only_freq_dependent=*/false):
  //     Rebuild both static and freq tensors. FIXED/CFS go into *_static_* at their
  //     profile's reference_frequency; FREQUENCY_DEPENDENT goes into *_freq_* at the
  //     passed ω (or zero if ω ≤ 0). Called once at setup.
  //
  //   RebuildPMLTensors(ω, /*only_freq_dependent=*/true):
  //     Rebuild ONLY *_freq_* at the passed ω. *_static_* is left untouched. Called by
  //     GetExtraSystemMatrix(ω) during a frequency sweep or NEP iteration.
  void RebuildPMLTensors(double omega, bool only_freq_dependent = false);

  const auto &GetAttributeToMaterial() const { return attr_mat; }
  mfem::Array<int> GetBdrAttributeToMaterial() const;

  template <typename T>
  auto GetCeedAttributes(const T &attr_list) const
  {
    return mesh.GetCeedAttributes(attr_list);
  }
  template <typename T>
  auto GetCeedBdrAttributes(const T &attr_list) const
  {
    return mesh.GetCeedBdrAttributes(attr_list);
  }

  auto MaxCeedAttribute() const { return mesh.MaxCeedAttribute(); }
  auto MaxCeedBdrAttribute() const { return mesh.MaxCeedBdrAttribute(); }

  const auto &GetMesh() const { return mesh; }
};

//
// Material property represented as a piecewise constant coefficient over domain or boundary
// mesh elements. Can be scalar-valued or matrix-valued. This should probably always operate
// at the level of libCEED attribute numbers (contiguous, 1-based) for consistency.
//
class MaterialPropertyCoefficient
{
private:
  // Map attribute to material index (coeff = mat_coeff[attr_mat[attr - 1]], for 1-based
  // attributes).
  mfem::Array<int> attr_mat;

  // Material property coefficients, ordered by material index.
  mfem::DenseTensor mat_coeff;

public:
  MaterialPropertyCoefficient(int attr_max);
  MaterialPropertyCoefficient(const mfem::Array<int> &attr_mat_,
                              const mfem::DenseTensor &mat_coeff_, double a = 1.0);

  bool empty() const { return mat_coeff.TotalSize() == 0; }

  const auto &GetAttributeToMaterial() const { return attr_mat; }
  const auto &GetMaterialProperties() const { return mat_coeff; }

  void AddCoefficient(const mfem::Array<int> &attr_mat_,
                      const mfem::DenseTensor &mat_coeff_, double a = 1.0);

  template <typename T>
  void AddMaterialProperty(const mfem::Array<int> &attr_list, const T &coeff,
                           double a = 1.0);
  template <typename T>
  void AddMaterialProperty(int attr, const T &coeff, double a = 1.0)
  {
    mfem::Array<int> attr_list(1);
    attr_list[0] = attr;
    AddMaterialProperty(attr_list, coeff, a);
  }

  MaterialPropertyCoefficient &operator*=(double a);

  void RestrictCoefficient(const mfem::Array<int> &attr_list);

  void NormalProjectedCoefficient(const mfem::Vector &normal);
};

}  // namespace palace

namespace palace::internal::mat
{

template <std::size_t N>
bool IsOrthonormal(const config::SymmetricMatrixData<N> &data);

template <std::size_t N>
bool IsValid(const config::SymmetricMatrixData<N> &data);

template <std::size_t N>
bool IsIsotropic(const config::SymmetricMatrixData<N> &data);

template <std::size_t N>
bool IsIdentity(const config::SymmetricMatrixData<N> &data);

}  // namespace palace::internal::mat

#endif  // PALACE_MODELS_MATERIAL_OPERATOR_HPP
