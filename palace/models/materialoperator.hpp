// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MATERIAL_OPERATOR_HPP
#define PALACE_MODELS_MATERIAL_OPERATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>
#include "fem/mesh.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class IoData;

//
// A class handling material attributes.
//
class MaterialOperator
{
public:
  // One nonzero pole-residue term r/(s - p) of a material susceptibility.
  struct PermittivityPoleTerm
  {
    std::complex<double> pole, residue;
  };

private:
  // Reference to underlying mesh object (not owned).
  const Mesh &mesh;

  // Mapping from the local libCEED attribute to material index.
  mfem::Array<int> attr_mat;

  // Material properties: relative permeability, relative permittivity, and others (like
  // electrical conductivity, London penetration depth for superconductors and Floquet wave
  // vector).
  mfem::DenseTensor mat_muinv, mat_muinv_scalar, mat_epsilon, mat_epsilon_scalar,
      mat_epsilon_imag, mat_epsilon_imag_scalar, mat_epsilon_abs, mat_invz0, mat_c0,
      mat_sigma, mat_sigma_scalar, mat_invLondon, mat_invLondon_scalar, mat_kxTmuinv,
      mat_muinvkx, mat_kxTmuinvkx, mat_kx;
  mfem::DenseMatrix wave_vector_cross;
  mfem::Vector wave_vector;        // BZ-wrapped k_F (fixed) or k₀ = k_F/ω (freq-scaled).
  mfem::Vector wave_vector_bz;     // BZ-wrapped k_F (always k_F, never k₀). For BZ offset.
  double floquet_omega_ref = 0.0;  // Nondimensional; when > 0, k_F scales with frequency.
  mfem::Array<double> mat_c0_min, mat_c0_max, mat_mu_eps_max;

  // Are materials isotropic? True when all the material properties are effectively
  // scalar-valued (ie, true scalars or vectors with identical entries). Also true when a
  // material is isotropic, the intersection is true when all are isotropic.
  mfem::Array<bool> attr_is_isotropic;

  // Flag for global domain attributes with nonzero loss tangent, electrical conductivity,
  // London penetration depth, Floquet wave vector, or frequency-dependent permittivity.
  bool has_losstan_attr, has_conductivity_attr, has_london_attr, has_wave_attr,
      has_frequency_dependent_permittivity;

  // Frequency-dependent permittivity state is indexed by the global position in
  // Domains.Materials, unlike attr_mat which is compacted independently on each rank.
  // Exact zero poles are folded into mat_sigma during setup; all remaining terms share a
  // cached spatial mass operator for each material.
  std::vector<std::vector<int>> frequency_dependent_permittivity_attributes;
  std::vector<std::vector<PermittivityPoleTerm>> permittivity_pole_terms;
  std::vector<std::vector<config::DjordjevicSarkarData>> djordjevic_sarkar_terms;
  // The terms as configured, for the stored energy of their models.
  std::vector<std::vector<config::PermittivityTermData>> permittivity_terms;
  // Global support of the nonlinear terms (A2) of each material, and of its configured
  // model (any frequency-dependent material, including one with only exact zero poles).
  std::vector<bool> frequency_dependent_permittivity_support,
      frequency_dependent_permittivity_model_support;
  // Global material index (position in Domains.Materials) of each rank-local material.
  mfem::Array<int> local_mat_to_material;
  // Largest relative permeability and the high-frequency permittivity of each
  // frequency-dependent material, for GetFrequencyDependentMaxMuEpsilon.
  std::vector<double> frequency_dependent_permittivity_mu_max,
      frequency_dependent_permittivity_eps_inf;
  // The user conductivity of each frequency-dependent material, if it is isotropic.
  std::vector<double> frequency_dependent_permittivity_sigma;
  std::vector<bool> frequency_dependent_permittivity_sigma_isotropic;

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
    MFEM_ASSERT(loc_attr.contains(attr),
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
  // Returns the scalar z-z (out-of-plane) component of the inverse permeability for a
  // given attribute. In 2D, this is stored in mat_muinv_scalar; in 3D, it is the (2,2)
  // entry of the full inverse permeability tensor.
  double GetInvPermeabilityZZ(int attr) const
  {
    if (mat_muinv.SizeI() == 2)
    {
      return Wrap(mat_muinv_scalar, attr)(0, 0);
    }
    return Wrap(mat_muinv, attr)(2, 2);
  }
  auto GetPermittivityReal(int attr) const { return Wrap(mat_epsilon, attr); }
  auto GetPermittivityImag(int attr) const { return Wrap(mat_epsilon_imag, attr); }
  auto GetPermittivityAbs(int attr) const { return Wrap(mat_epsilon_abs, attr); }
  // Scalar (1x1) permittivity for the out-of-plane component in 2D mode analysis.
  const auto &GetPermittivityScalar() const { return mat_epsilon_scalar; }
  // Scalar (1x1) imaginary permittivity for the out-of-plane component in 2D.
  const auto &GetPermittivityImagScalar() const { return mat_epsilon_imag_scalar; }
  // Scalar (1x1) conductivity for the out-of-plane component in 2D mode analysis.
  const auto &GetConductivityScalar() const { return mat_sigma_scalar; }
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

  // Returns scalar (1x1) inverse permeability for 2D curl-curl (curl is scalar in 2D),
  // or the full matrix for 3D.
  const auto &GetCurlCurlInvPermeability() const
  {
    return (mat_muinv.SizeI() == 2) ? mat_muinv_scalar : mat_muinv;
  }
  const auto &GetInvPermeability() const { return mat_muinv; }
  const auto &GetPermittivityReal() const { return mat_epsilon; }
  const auto &GetPermittivityImag() const { return mat_epsilon_imag; }
  const auto &GetPermittivityAbs() const { return mat_epsilon_abs; }
  const auto &GetInvImpedance() const { return mat_invz0; }
  const auto &GetLightSpeed() const { return mat_c0; }
  const auto &GetConductivity() const { return mat_sigma; }
  const auto &GetInvLondonDepth() const { return mat_invLondon; }
  // Scalar (1x1) out-of-plane London depth for the H1 Ann mass in 2D mode analysis.
  const auto &GetInvLondonDepthScalar() const { return mat_invLondon_scalar; }
  const auto &GetFloquetCurl() const { return mat_muinvkx; }
  const auto &GetFloquetMass() const { return mat_kxTmuinvkx; }
  const auto &GetFloquetCross() const { return mat_kx; }

  const auto &GetLightSpeedMin() const { return mat_c0_min; }
  const auto &GetLightSpeedMax() const { return mat_c0_max; }

  // Conservative upper bound on (k_n / omega)^2 for bulk waves, computed as
  // lambda_max(mu_r) * lambda_max(epsilon_r) and maximized over this mesh's materials.
  // Keeping the tensor extrema separate is required for anisotropic media because the
  // electric and magnetic polarizations generally lie on different material axes.
  double GetMaxMuEpsilon() const;

  // max(mu_r) * Re(epsilon_r(omega)) over the globally supported frequency-dependent
  // materials, or zero without any. GetMaxMuEpsilon uses the high-frequency permittivity,
  // which a dispersive material exceeds below its resonances. Not collective.
  double GetFrequencyDependentMaxMuEpsilon(double omega) const;

  // Split the libCEED boundary attributes of the given mesh boundary attributes by the
  // frequency-dependent material (any configured model) of their neighboring domain.
  // Returns those next to no such material, and the index and attributes of each material
  // that neighbors the boundaries on any process (the attributes may be empty on this one).
  struct FrequencyDependentBdrSplit
  {
    mfem::Array<int> other;
    std::vector<std::pair<std::size_t, mfem::Array<int>>> materials;
  };
  FrequencyDependentBdrSplit SplitCeedBdrAttributesByFrequencyDependentPermittivity(
      const mfem::Array<int> &bdr_attrs) const;
  double GetFrequencyDependentPermittivityHighFrequency(std::size_t material_idx) const
  {
    return frequency_dependent_permittivity_eps_inf.at(material_idx);
  }

  bool HasLossTangent() const { return has_losstan_attr; }
  bool HasConductivity() const { return has_conductivity_attr; }
  bool HasFrequencyDependentPermittivity() const
  {
    return has_frequency_dependent_permittivity;
  }
  std::size_t NumFrequencyDependentPermittivityMaterials() const
  {
    return frequency_dependent_permittivity_attributes.size();
  }
  bool HasFrequencyDependentPermittivityA2(std::size_t material_idx) const
  {
    return !permittivity_pole_terms.at(material_idx).empty() ||
           !djordjevic_sarkar_terms.at(material_idx).empty();
  }
  // Whether the nonlinear terms (A2) of the material have mesh support on any process. This
  // is the predicate of the operators that hold A2, not of the complete permittivity.
  bool HasFrequencyDependentPermittivitySupport(std::size_t material_idx) const
  {
    return frequency_dependent_permittivity_support.at(material_idx);
  }
  // Whether the material has a configured frequency-dependent model with mesh support on
  // any process, including a model with only exact zero poles (no A2). This is the
  // predicate of every consumer of the complete permittivity ε(s) or the stored energy.
  bool HasFrequencyDependentPermittivityModel(std::size_t material_idx) const
  {
    return frequency_dependent_permittivity_model_support.at(material_idx);
  }
  const std::vector<int> &
  GetFrequencyDependentPermittivityAttributes(std::size_t material_idx) const
  {
    return frequency_dependent_permittivity_attributes.at(material_idx);
  }
  // Evaluate the complete frequency-dependent contribution to the system coefficient.
  // The high-frequency permittivity and exact zero poles are handled by the ordinary mass
  // and conductivity paths, respectively.
  std::complex<double>
  EvaluateFrequencyDependentPermittivityA2(std::size_t material_idx,
                                           std::complex<double> s) const;
  // The two parts of EvaluateFrequencyDependentPermittivityA2: the pole terms, with complex
  // poles expanded into conjugate pairs and each contributing r·s²/(s − p), and the
  // remaining non-rational Djordjevic-Sarkar contribution.
  const std::vector<PermittivityPoleTerm> &
  GetPermittivityPoleTerms(std::size_t material_idx) const
  {
    return permittivity_pole_terms.at(material_idx);
  }
  std::complex<double> EvaluateDjordjevicSarkarPermittivityA2(std::size_t material_idx,
                                                              std::complex<double> s) const;
  // The complete relative permittivity ε∞ + χ(s) of a frequency-dependent material, at
  // s ≠ 0, from its configured terms, including exact zero poles. Errors at a singularity.
  std::complex<double> EvaluateFrequencyDependentPermittivity(std::size_t material_idx,
                                                              std::complex<double> s) const;
  // The wave permittivity ε(s) + σ/s of a frequency-dependent material, with its user
  // Conductivity σ (a free current, not part of D or the stored energy), for the
  // propagation constant and impedance of an absorbing boundary or Floquet port. Unlike the
  // conductivity matrix, σ excludes the zero poles, which ε(s) already contains. Errors for
  // an anisotropic σ.
  std::complex<double>
  EvaluateFrequencyDependentWavePermittivity(std::size_t material_idx,
                                             std::complex<double> s) const;
  // ε∞ + StoredEnergySusceptibility at real ω > 0: the constant permittivity with the same
  // time-averaged stored electric energy as the model, possibly ≤ 0 for a general
  // PoleResidue model. Errors at a singularity.
  double GetFrequencyDependentStoredEnergyPermittivity(std::size_t material_idx,
                                                       double omega) const;

  // Material tensors at the real angular frequency ω > 0, in the layout of
  // GetPermittivityReal() (indexed by the rank-local material, for
  // MaterialPropertyCoefficient). Constant materials keep their entries, Re{ε}, Im{ε} =
  // −Re{ε} tan δ, and Re{ε} for the stored energy; a frequency-dependent material has
  // Re ε(iω), Im ε(iω), or GetFrequencyDependentStoredEnergyPermittivity. These are values,
  // not state: the getters above always return the high-frequency permittivity ε∞.
  mfem::DenseTensor EvaluatePermittivityReal(double omega) const;
  mfem::DenseTensor EvaluatePermittivityImag(double omega) const;
  mfem::DenseTensor EvaluateStoredEnergyPermittivity(double omega) const;

  // The light speed (μ ε)^(-1/2) at the real angular frequency ω > 0, in the layout of
  // GetLightSpeed(), with Re ε(iω) for the given frequency-dependent materials (global
  // indices, the same on every process): the wave speed of a lossless medium with this
  // permittivity. Errors unless Re ε(iω) > 0 for each, on every process. Other
  // frequency-dependent materials keep ε∞ and are not evaluated.
  mfem::DenseTensor EvaluateLightSpeed(double omega,
                                       const std::vector<std::size_t> &materials) const;

  // Whether the domain attribute is vacuum for interface postprocessing: a light speed of
  // one and no frequency-dependent permittivity, independent of the frequency (a
  // frequency-dependent material with ε∞ = 1 is not vacuum).
  bool IsVacuum(int attr) const;

  // A copy of the material tensor data (in the layout of GetPermittivityReal()), with
  // f(m, M) applied to the entry M of each rank-local material with the configured
  // frequency-dependent model m.
  template <typename F>
  mfem::DenseTensor MapFrequencyDependentMaterials(const mfem::DenseTensor &data,
                                                   F &&f) const
  {
    mfem::DenseTensor out(data);
    for (int k = 0; k < out.SizeK(); k++)
    {
      const int m = GetLocalFrequencyDependentPermittivityModelMaterial(k);
      if (m >= 0)
      {
        f(static_cast<std::size_t>(m), out(k));
      }
    }
    return out;
  }

  // The frequency-dependent material with a configured model on the given domain attribute
  // (unlike GetFrequencyDependentPermittivityMaterial, also one with only zero poles), or
  // -1 if there is none. The same for a rank-local material index.
  int GetFrequencyDependentPermittivityModelMaterial(int attr) const;
  int GetLocalFrequencyDependentPermittivityModelMaterial(int local_mat) const;

  // The rank-local material index of a domain attribute, for the evaluated tensors above.
  int GetLocalMaterial(int attr) const { return AttrToMat(attr); }
  // The entry of a material tensor (in the layout of GetPermittivityReal()) on a domain
  // attribute.
  auto GetMaterialEntry(const mfem::DenseTensor &data, int attr) const
  {
    return Wrap(data, attr);
  }
  // The time-averaged stored electric energy of the polarization of a frequency-dependent
  // material at real ω > 0, as a susceptibility (permittivity::StoredEnergySusceptibility):
  // the stored energy is that of a constant permittivity ε∞ plus this value. It includes
  // exact zero poles, which store none.
  double GetFrequencyDependentStoredEnergySusceptibility(std::size_t material_idx,
                                                         double omega) const;
  // Re χ(iω) of the nonlinear terms of a frequency-dependent material, the field energy at
  // its permittivity ε(iω) less that at ε∞ (a zero pole contributes none).
  double GetFrequencyDependentSusceptibilityReal(std::size_t material_idx,
                                                 double omega) const;
  // The frequency-dependent material with a nonlinear term on the given domain attribute,
  // or -1 if there is none.
  int GetFrequencyDependentPermittivityMaterial(int attr) const;
  bool HasDjordjevicSarkarPermittivity(std::size_t material_idx) const
  {
    return !djordjevic_sarkar_terms.at(material_idx).empty();
  }
  // True when a globally supported material has a nonlinear permittivity term. This is
  // deliberately frequency independent: valid additive contributions can cancel at an
  // interpolation or ROM probe frequency without ceasing to be nonlinear.
  bool HasFrequencyDependentPermittivityA2() const;
  // Distinct poles (in s = iω) of the nonlinear permittivity terms of all globally
  // supported materials, i.e. the poles of their contribution to A2(s).
  std::vector<std::complex<double>> GetFrequencyDependentPermittivityPoles() const;
  bool HasLondonDepth() const { return has_london_attr; }
  bool HasWaveVector() const { return has_wave_attr; }
  const mfem::Vector &GetWaveVector() const { return wave_vector; }
  const mfem::Vector &GetWaveVectorBZ() const { return wave_vector_bz; }
  bool HasFloquetFrequencyScaling() const { return floquet_omega_ref > 0.0; }
  double GetFloquetOmegaRef() const { return floquet_omega_ref; }

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

  // Mathematically exact zero: empty storage is zero, otherwise every stored scalar must
  // compare exactly equal to zero. This deliberately has no tolerance.
  bool IsExactlyZero() const;

  // Append the (libCEED, 1-based) attributes whose material property tensor has at least
  // one entry which is not exactly zero. Attributes with no assigned material, or an
  // exactly zero tensor, contribute exactly zero and are not appended. The output is not
  // sorted or made unique by this call.
  void AddNonzeroAttributes(std::vector<int> &attr_list) const;

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

// Rotate each 3x3 config-file symmetric tensor (mu_r, epsilon_r, tandelta, sigma) in
// `materials` into a local frame whose in-plane axes are (e1, e2) and out-of-plane axis
// is `normal`, zeroing cross terms between the two. After this call, a 2D MaterialOperator
// constructed from the mutated iodata yields the physically correct in-plane 2x2 block
// (leading submatrix = R^T M R) and scalar out-of-plane components ((2,2) = n^T M n) via
// the default 2D construction path — no additional frame information is required
// downstream. `lambda_L` is a scalar length, unaffected.
void RotateMaterialDefinitions(std::vector<config::MaterialData> &materials,
                               const mfem::Vector &e1, const mfem::Vector &e2,
                               const mfem::Vector &normal);

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
