// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialoperator.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include "fem/libceed/ceed.hpp"  // for <ceed.h> before coeff_qf.h
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "linalg/densematrix.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace internal::mat
{

template <std::size_t N>
bool IsOrthonormal(const config::SymmetricMatrixData<N> &data)
{

  // All the vectors are normalized.
  constexpr auto tol = 1.0e-6;
  auto UnitNorm = [&](const std::array<double, N> &x)
  {
    double s = -1.0;
    for (const auto &i : x)
    {
      s += std::pow(i, 2);
    }
    return std::abs(s) < tol;
  };
  bool valid = std::all_of(data.v.begin(), data.v.end(), UnitNorm);

  // All the vectors are orthogonal.
  for (std::size_t i1 = 0; i1 < N; i1++)
  {
    const auto &v1 = data.v.at(i1);
    for (std::size_t i2 = i1 + 1; i2 < N; i2++)
    {
      const auto &v2 = data.v.at(i2);
      double s = 0.0;
      for (std::size_t j = 0; j < N; j++)
      {
        s += v1[j] * v2[j];
      }
      valid &= std::abs(s) < tol;
    }
  }
  return valid;
}

template <std::size_t N>
bool IsValid(const config::SymmetricMatrixData<N> &data)
{
  return IsOrthonormal(data) && std::all_of(data.s.begin(), data.s.end(),
                                            [](auto d) { return std::abs(d) > 0.0; });
}

template <std::size_t N>
bool IsIsotropic(const config::SymmetricMatrixData<N> &data)
{
  return IsOrthonormal(data) &&
         std::all_of(data.s.begin(), data.s.end(), [&](auto d) { return d == data.s[0]; });
}

template <std::size_t N>
bool IsIdentity(const config::SymmetricMatrixData<N> &data)
{
  return IsOrthonormal(data) &&
         std::all_of(data.s.begin(), data.s.end(), [](auto d) { return d == 1.0; });
}

template <std::size_t N>
mfem::DenseMatrix ToDenseMatrix(const config::SymmetricMatrixData<N> &data)
{
  mfem::DenseMatrix M(N, N);
  mfem::Vector V(N);
  for (std::size_t i = 0; i < N; i++)
  {
    for (std::size_t j = 0; j < N; j++)
    {
      V(j) = data.v[i][j];
    }
    AddMult_a_VVt(data.s[i], V, M);
  }
  return M;
}

// Extract the leading sdim x sdim submatrix from a SymmetricMatrixData<N> matrix.
template <std::size_t N>
mfem::DenseMatrix ToDenseMatrixTruncated(const config::SymmetricMatrixData<N> &data,
                                         int sdim)
{
  auto M = ToDenseMatrix(data);
  if (sdim >= static_cast<int>(N))
  {
    return M;
  }
  mfem::DenseMatrix Msub(sdim, sdim);
  for (int i = 0; i < sdim; i++)
  {
    for (int j = 0; j < sdim; j++)
    {
      Msub(i, j) = M(i, j);
    }
  }
  return Msub;
}

}  // namespace internal::mat

namespace
{

// Union-bounding-box pre-pass over all PML materials that use auto-detected geometry,
// so N stacked slabs (multiple attributes, possibly across multiple materials) share
// one logical profile geometry: same inner interface, same total thickness, same
// σ_max. Without this, each slab's σ(z) resets to 0 at its inner face and auto-σ_max
// is computed from the per-slab (smaller) thickness — which gives substantially worse
// absorption than an unsplit layer of the same total thickness.
pml::SlabGeometry
DetectUnifiedPMLSlabGeometry(const std::vector<config::MaterialData> &materials,
                             const mfem::ParMesh &mesh, int attr_max, int attr_max_local,
                             const std::array<double, 3> &bbmin,
                             const std::array<double, 3> &bbmax)
{
  mfem::Array<int> union_marker(attr_max);
  union_marker = 0;
  bool any_autodetect = false;
  for (const auto &data : materials)
  {
    if (!data.pml.has_value() || !data.pml->autodetect_geometry)
    {
      continue;
    }
    any_autodetect = true;
    for (auto attr : data.attributes)
    {
      if (attr > 0 && attr <= attr_max_local)
      {
        union_marker[attr - 1] = 1;
      }
    }
  }
  if (!any_autodetect)
  {
    return {};
  }
  // GetAxisAlignedBoundingBox does a global MPI reduction, so every rank must call it.
  mfem::Vector u_min_mfem, u_max_mfem;
  mesh::GetAxisAlignedBoundingBox(mesh, union_marker, /*bdr=*/false, u_min_mfem,
                                  u_max_mfem);
  const std::array<double, 3> u_min{{u_min_mfem[0], u_min_mfem[1], u_min_mfem[2]}};
  const std::array<double, 3> u_max{{u_max_mfem[0], u_max_mfem[1], u_max_mfem[2]}};

  // If none of the auto-detect PML attributes are present anywhere on this mesh, the
  // union bounding box is empty (min > max on every axis: GetAxisAlignedBoundingBox
  // seeds min/max with +inf/-inf and no element updates them). This happens when the
  // MaterialOperator is built on a mesh that does not contain the PML region — e.g. the
  // wave-port boundary submesh, which carries the full material list but only the port-
  // face attributes. Return an empty geometry quietly; the per-material loop below maps
  // nothing on such a mesh, so the PML simply does not participate there.
  const bool union_empty =
      (u_min[0] > u_max[0]) && (u_min[1] > u_max[1]) && (u_min[2] > u_max[2]);
  if (union_empty)
  {
    return {};
  }

  auto geom = pml::DetectSlabGeometry(u_min, u_max, bbmin, bbmax);
  const bool any_active =
      std::any_of(geom.direction_signs.begin(), geom.direction_signs.end(),
                  [](int s) { return s != 0; });
  MFEM_VERIFY(any_active,
              "Auto-detected PML geometry found no PML attributes touching the global "
              "mesh bounding box. Ensure the PML region sits on the outer boundary of "
              "the mesh, or specify \"Direction\" and \"Thickness\" explicitly in each "
              "PML material.");
  return geom;
}

}  // namespace

MaterialOperator::MaterialOperator(const std::vector<config::MaterialData> &materials,
                                   const config::PeriodicBoundaryData &periodic,
                                   ProblemType problem_type, const Mesh &mesh)
  : mesh(mesh)
{
  SetUpMaterialProperties(materials, periodic, problem_type, mesh);
}

MaterialOperator::MaterialOperator(const IoData &iodata, const Mesh &mesh)
  : MaterialOperator(iodata.domains.materials, iodata.boundaries.periodic,
                     iodata.problem.type, mesh)
{
}

void MaterialOperator::SetUpMaterialProperties(
    const std::vector<config::MaterialData> &materials,
    const config::PeriodicBoundaryData &periodic, ProblemType problem_type,
    const mfem::ParMesh &mesh)
{
  // Check material attributes. Only verify positivity here — the per-material local
  // presence check happens below via mat_marker (some ranks may not have all attributes
  // in their partition).
  MFEM_VERIFY(!materials.empty(), "Materials must be non-empty!");
  for (const auto &data : materials)
  {
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(attr > 0, "Material attribute tags must be positive!");
    }
  }

  // Set up material properties of the different domain regions, represented with element-
  // wise constant matrix-valued coefficients for the relative permeability, permittivity,
  // and other material properties.
  const auto &loc_attr = this->mesh.GetCeedAttributes();
  mfem::Array<int> mat_marker(materials.size());
  mat_marker = 0;
  int nmats = 0;
  for (std::size_t i = 0; i < materials.size(); i++)
  {
    const auto &data = materials[i];
    for (auto attr : data.attributes)
    {
      if (loc_attr.find(attr) != loc_attr.end())
      {
        mat_marker[i] = 1;
        nmats++;
        break;
      }
    }
  }
  attr_mat.SetSize(loc_attr.size());
  attr_mat = -1;

  attr_is_isotropic.SetSize(nmats);

  const int sdim = mesh.SpaceDimension();
  mat_muinv.SetSize(sdim, sdim, nmats);
  if (sdim == 2)
  {
    mat_muinv_scalar.SetSize(1, 1, nmats);
    mat_epsilon_scalar.SetSize(1, 1, nmats);
    mat_epsilon_imag_scalar.SetSize(1, 1, nmats);
    mat_invLondon_scalar.SetSize(1, 1, nmats);
  }
  mat_epsilon.SetSize(sdim, sdim, nmats);
  mat_epsilon_imag.SetSize(sdim, sdim, nmats);
  mat_epsilon_abs.SetSize(sdim, sdim, nmats);
  mat_invz0.SetSize(sdim, sdim, nmats);
  mat_c0.SetSize(sdim, sdim, nmats);
  mat_sigma.SetSize(sdim, sdim, nmats);
  mat_invLondon.SetSize(sdim, sdim, nmats);
  mat_c0_min.SetSize(nmats);
  mat_c0_max.SetSize(nmats);
  mat_muinvkx.SetSize(sdim, sdim, nmats);
  mat_kxTmuinvkx.SetSize(sdim, sdim, nmats);
  mat_kx.SetSize(sdim, sdim, nmats);
  has_losstan_attr = has_conductivity_attr = has_london_attr = has_wave_attr = false;

  // PML profile registry: a per-libCEED-attribute index into pml_profiles (or -1 for
  // non-PML attributes). The PML QFunction reads from the packed context (pml_ctx)
  // built at the end of this function; it looks up the profile index from the element
  // attribute and evaluates the stretch tensor at each quadrature point from the
  // cached physical coordinate.
  pml_attr_to_profile.assign(attr_mat.Size(), -1);
  pml_profiles.clear();
  pml_ctx.clear();
  has_pml_attr = false;
  has_pml_freq_dependent_attr = false;

  // Set up Floquet wave vector for periodic meshes with phase-delay constraints.
  SetUpFloquetWaveVector(periodic, problem_type, mesh);

  int count = 0;
  for (std::size_t i = 0; i < materials.size(); i++)
  {
    if (!mat_marker[i])
    {
      continue;
    }
    const auto &data = materials[i];
    if (problem_type == ProblemType::ELECTROSTATIC)
    {
      MFEM_VERIFY(internal::mat::IsValid(data.epsilon_r),
                  "Material has no valid permittivity defined!");
      if (!internal::mat::IsIdentity(data.mu_r) || internal::mat::IsValid(data.sigma) ||
          std::abs(data.lambda_L) > 0.0)
      {
        Mpi::Warning(
            "Electrostatic problem type does not account for material permeability,\n"
            "electrical conductivity, or London depth!\n");
      }
    }
    else if (problem_type == ProblemType::MAGNETOSTATIC)
    {
      MFEM_VERIFY(internal::mat::IsValid(data.mu_r),
                  "Material has no valid permeability defined!");
      if (!internal::mat::IsIdentity(data.epsilon_r) ||
          internal::mat::IsValid(data.tandelta) || internal::mat::IsValid(data.sigma) ||
          std::abs(data.lambda_L) > 0.0)
      {
        Mpi::Warning(
            "Magnetostatic problem type does not account for material permittivity,\n"
            "loss tangent, electrical conductivity, or London depth!\n");
      }
    }
    else
    {
      MFEM_VERIFY(internal::mat::IsValid(data.mu_r) &&
                      internal::mat::IsValid(data.epsilon_r),
                  "Material has no valid permeability or no valid permittivity defined!");
      if (problem_type == ProblemType::TRANSIENT)
      {
        MFEM_VERIFY(!internal::mat::IsValid(data.tandelta),
                    "Transient problem type does not support material loss tangent, use "
                    "electrical conductivity instead!");
      }
      else
      {
        MFEM_VERIFY(
            !(internal::mat::IsValid(data.tandelta) && internal::mat::IsValid(data.sigma)),
            "Material loss model should probably use only one of loss tangent or "
            "electrical conductivity!");
      }
    }

    attr_is_isotropic[count] = internal::mat::IsIsotropic(data.mu_r) &&
                               internal::mat::IsIsotropic(data.epsilon_r) &&
                               internal::mat::IsIsotropic(data.tandelta) &&
                               internal::mat::IsIsotropic(data.sigma);

    // Map all attributes to this material property index.
    for (auto attr : data.attributes)
    {
      auto it = loc_attr.find(attr);
      if (it != loc_attr.end())
      {
        MFEM_VERIFY(
            attr_mat[it->second - 1] < 0,
            "Detected multiple definitions of material properties for domain attribute "
                << attr << "!");
        attr_mat[it->second - 1] = count;
      }
    }

    // Compute the inverse of the input permeability matrix. Use truncated versions
    // of the config 3x3 matrices for 2D meshes.
    mfem::DenseMatrix mat_mu = internal::mat::ToDenseMatrixTruncated(data.mu_r, sdim);
    mfem::DenseMatrixInverse(mat_mu, true).GetInverseMatrix(mat_muinv(count));
    if (sdim == 2)
    {
      // In 2D, compute out-of-plane scalar components using n^T M_3x3 n projection.
      // For axis-aligned meshes n = [0, 0, 1], giving the (2, 2) entry; this remains
      // correct for non-axis-aligned cross-sections once the driver has pre-rotated
      // the config materials into the local frame via RotateMaterialDefinitions.
      const mfem::Vector normal({0.0, 0.0, 1.0});
      auto ProjectNormal = [&normal](const mfem::DenseMatrix &M3) -> double
      { return M3.InnerProduct(normal, normal); };

      mfem::DenseMatrix mat_mu_3d = internal::mat::ToDenseMatrix(data.mu_r);
      mfem::DenseMatrix mat_muinv_3d(3, 3);
      mfem::DenseMatrixInverse(mat_mu_3d, true).GetInverseMatrix(mat_muinv_3d);
      mat_muinv_scalar(count)(0, 0) = ProjectNormal(mat_muinv_3d);

      mfem::DenseMatrix mat_eps_3d = internal::mat::ToDenseMatrix(data.epsilon_r);
      mat_epsilon_scalar(count)(0, 0) = ProjectNormal(mat_eps_3d);

      mfem::DenseMatrix mat_td_3d = internal::mat::ToDenseMatrix(data.tandelta);
      mfem::DenseMatrix epstd_3d(3, 3);
      Mult(mat_eps_3d, mat_td_3d, epstd_3d);
      mat_epsilon_imag_scalar(count)(0, 0) = -ProjectNormal(epstd_3d);
    }

    // Material permittivity: Re{ε} = ε, Im{ε} = -ε * tan(δ)
    mfem::DenseMatrix T(sdim, sdim);
    mat_epsilon(count).Set(1.0,
                           internal::mat::ToDenseMatrixTruncated(data.epsilon_r, sdim));
    Mult(mat_epsilon(count), internal::mat::ToDenseMatrixTruncated(data.tandelta, sdim), T);
    T *= -1.0;
    mat_epsilon_imag(count).Set(1.0, T);
    if (mat_epsilon_imag(count).MaxMaxNorm() > 0.0)
    {
      has_losstan_attr = true;
    }

    // ε * √(I + tan(δ) * tan(δ)ᵀ)
    MultAAt(internal::mat::ToDenseMatrixTruncated(data.tandelta, sdim), T);
    for (int d = 0; d < T.Height(); d++)
    {
      T(d, d) += 1.0;
    }
    Mult(mat_epsilon(count), linalg::MatrixSqrt(T), mat_epsilon_abs(count));

    // √(μ⁻¹ ε)
    Mult(mat_muinv(count), mat_epsilon(count), mat_invz0(count));
    mat_invz0(count).Set(1.0, linalg::MatrixSqrt(mat_invz0(count)));

    // √((μ ε)⁻¹)
    Mult(mat_mu, mat_epsilon(count), T);
    mat_c0(count).Set(1.0, linalg::MatrixPow(T, -0.5));
    mat_c0_min[count] = linalg::SingularValueMin(mat_c0(count));
    mat_c0_max[count] = linalg::SingularValueMax(mat_c0(count));

    // Electrical conductivity, σ
    mat_sigma(count).Set(1.0, internal::mat::ToDenseMatrixTruncated(data.sigma, sdim));
    if (mat_sigma(count).MaxMaxNorm() > 0.0)
    {
      has_conductivity_attr = true;
    }

    // λ⁻² * μ⁻¹
    {
      double invL2 = std::abs(data.lambda_L) > 0.0 ? std::pow(data.lambda_L, -2.0) : 0.0;
      mat_invLondon(count).Set(1.0, mat_muinv(count));
      mat_invLondon(count) *= invL2;
      if (sdim == 2)
      {
        // Scalar out-of-plane London depth: λ⁻² * n^T μ⁻¹ n.
        mat_invLondon_scalar(count)(0, 0) = mat_muinv_scalar(count)(0, 0) * invL2;
      }
      if (mat_invLondon(count).MaxMaxNorm() > 0.0)
      {
        has_london_attr = true;
      }
    }

    // μ⁻¹ [k x]
    Mult(mat_muinv(count), wave_vector_cross, mat_muinvkx(count));

    // [k x]^T μ⁻¹ [k x]
    T.Transpose(wave_vector_cross);
    Mult(T, mat_muinvkx(count), mat_kxTmuinvkx(count));

    // [k x]
    mat_kx(count).Set(1.0, wave_vector_cross);

    count++;
  }
  bool has_attr[4] = {has_losstan_attr, has_conductivity_attr, has_london_attr,
                      has_wave_attr};
  Mpi::GlobalOr(4, has_attr, mesh.GetComm());
  has_losstan_attr = has_attr[0];
  has_conductivity_attr = has_attr[1];
  has_london_attr = has_attr[2];
  has_wave_attr = has_attr[3];

  // ---- PML setup (second pass) --------------------------------------------------
  // For each material with a PML block, build a pml::Profile and zero out the bulk
  // material tensors on PML attributes so the standard coefficient path contributes
  // nothing. PML contributions are applied per-quadrature-point in the PML QFunction
  // (see fem/qfunctions/coeff/pml_qf.h) driven by a packed context built below.
  mfem::Vector bbmin_mfem, bbmax_mfem;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin_mfem, bbmax_mfem);
  const std::array<double, 3> bbmin{{bbmin_mfem[0], bbmin_mfem[1], bbmin_mfem[2]}};
  const std::array<double, 3> bbmax{{bbmax_mfem[0], bbmax_mfem[1], bbmax_mfem[2]}};

  // Per-rank attribute-max for marker sizing; global max across ranks lets us size the
  // per-attribute bounding-box marker uniformly so every rank makes the same collective
  // calls below even if a PML material's elements happen to all live on other ranks.
  int attr_max_local = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
  int attr_max = attr_max_local;
  Mpi::GlobalMax(1, &attr_max, mesh.GetComm());

  const auto unified_geom =
      DetectUnifiedPMLSlabGeometry(materials, mesh, attr_max, attr_max_local, bbmin, bbmax);

  for (std::size_t i = 0; i < materials.size(); i++)
  {
    if (!materials[i].pml.has_value())
    {
      continue;
    }
    const auto &data = materials[i];
    // Local copy that may have direction_signs / thickness filled in by auto-detection.
    auto pml_cfg = *data.pml;

    // Auto-detected geometry uses the unified slab geometry computed in the pre-pass
    // above, so all PML materials that participate in auto-detection share the same
    // thickness and inner interface. This makes N-slab and 1-slab layouts equivalent.
    if (pml_cfg.autodetect_geometry)
    {
      pml_cfg.direction_signs = unified_geom.direction_signs;
      pml_cfg.thickness = unified_geom.thickness;
    }

    // Build the Profile. Use the isotropic average of the anisotropic ε_r / μ_r for the
    // reflection-target σ_max auto-computation (a PML's reflection coefficient depends
    // on the scalar wave impedance, which the averaged values approximate).
    const auto mu_avg = (data.mu_r.s[0] + data.mu_r.s[1] + data.mu_r.s[2]) / 3.0;
    const auto eps_avg =
        (data.epsilon_r.s[0] + data.epsilon_r.s[1] + data.epsilon_r.s[2]) / 3.0;
    auto profile = pml::BuildProfile(pml_cfg, mu_avg, eps_avg);

    // Interface coordinates from the mesh bounding box. For each active face, the
    // inner interface is at bbmin[axis] + thickness_neg (for −face) or
    // bbmax[axis] − thickness_pos (for +face). The PML QFunction uses these to compute
    // the depth into the absorbing layer from the physical QP coordinate.
    for (int axis = 0; axis < 3; axis++)
    {
      profile.interface_coord[axis][0] = bbmin[axis] + pml_cfg.thickness[2 * axis + 0];
      profile.interface_coord[axis][1] = bbmax[axis] - pml_cfg.thickness[2 * axis + 1];
    }

    const int profile_idx = static_cast<int>(pml_profiles.size());
    pml_profiles.push_back(profile);

    // Map each of this material's libCEED attributes to this profile, and zero out the
    // bulk mat_muinv / mat_epsilon entries for those attributes so the standard
    // coefficient path contributes nothing. Only local attributes (those present in
    // this rank's CEED attribute map) are mapped.
    for (auto attr : data.attributes)
    {
      auto it = loc_attr.find(attr);
      if (it == loc_attr.end())
      {
        continue;
      }
      const int ceed_attr = it->second;  // 1-based
      pml_attr_to_profile[ceed_attr - 1] = profile_idx;

      // Zero out the bulk material tensors for this attribute so the bulk coefficient
      // path contributes nothing (PML tensors are assembled separately via per-QP
      // integrators).
      const int mat_idx = attr_mat[ceed_attr - 1];
      if (mat_idx >= 0)
      {
        mat_muinv(mat_idx) = 0.0;
        mat_epsilon(mat_idx) = 0.0;
        mat_epsilon_imag(mat_idx) = 0.0;
        mat_epsilon_abs(mat_idx) = 0.0;
      }
    }

    has_pml_attr = true;
    if (pml_cfg.frequency_dependent)
    {
      has_pml_freq_dependent_attr = true;
    }
  }

  bool has_pml_buf[2] = {has_pml_attr, has_pml_freq_dependent_attr};
  Mpi::GlobalOr(2, has_pml_buf, mesh.GetComm());
  has_pml_attr = has_pml_buf[0];
  has_pml_freq_dependent_attr = has_pml_buf[1];

  // Pack the PML QFunction context once at setup. Per-region ω fields for
  // FIXED/CFS profiles are already set (to their reference_frequency). FREQUENCY_DEPENDENT
  // regions have ω = 0 and will be refreshed by SpaceOperator::GetExtraSystemMatrix(ω).
  if (has_pml_attr)
  {
    pml_ctx = pml::PackProfileContextAll(pml_attr_to_profile, pml_profiles);
  }

  // Print a summary of PML regions on the root rank.
  if (has_pml_attr && Mpi::Root(mesh.GetComm()) && !pml_profiles.empty())
  {
    Mpi::Print("\nConfigured PML regions ({} profile{} on rank 0):\n", pml_profiles.size(),
               pml_profiles.size() == 1 ? "" : "s");
    constexpr std::array<const char *, 6> face_name = {"-x", "+x", "-y", "+y", "-z", "+z"};
    for (std::size_t p = 0; p < pml_profiles.size(); p++)
    {
      const auto &prof = pml_profiles[p];
      std::string faces;
      for (int f = 0; f < 6; f++)
      {
        if (prof.direction_active[f])
        {
          if (!faces.empty())
            faces += ",";
          faces += face_name[f];
        }
      }
      Mpi::Print(" Profile {}: faces [{}], thickness [{:.3e}, {:.3e}, {:.3e}]\n", p, faces,
                 prof.thickness[0], prof.thickness[1], prof.thickness[2]);
    }
  }
}

void MaterialOperator::RefreshPMLContextFrequency(double omega) const
{
  if (!has_pml_attr || pml_ctx.empty())
  {
    return;
  }
  pml::RefreshPMLContextFrequency(pml_ctx.data(), GetPMLNumAttributes(),
                                  GetPMLNumProfiles(), omega);
}

void MaterialOperator::SetUpFloquetWaveVector(const config::PeriodicBoundaryData &periodic,
                                              ProblemType problem_type,
                                              const mfem::ParMesh &mesh)
{
  const int sdim = mesh.SpaceDimension();
  const double tol = std::numeric_limits<double>::epsilon();

  // Get Floquet wave vector (stored as member for access by FloquetPortData).
  wave_vector.SetSize(sdim);
  wave_vector = 0.0;
  MFEM_VERIFY(static_cast<int>(periodic.wave_vector.size()) >= sdim,
              "Floquet wave vector size must be at least the spatial dimension.");
  std::copy_n(periodic.wave_vector.begin(), sdim, wave_vector.GetData());
  has_wave_attr = (wave_vector.Norml2() > tol);

  MFEM_VERIFY(!has_wave_attr || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE,
              "Quasi-periodic Floquet boundary conditions are only available for "
              " frequency domain driven or eigenmode simulations!");
  MFEM_VERIFY(!has_wave_attr || sdim == 3,
              "Quasi-periodic Floquet periodic boundary conditions are only available "
              " in 3D!");

  // Store nondimensional reference angular frequency for frequency-dependent k_F scaling.
  // Already nondimensionalized in IoData::NondimensionalizeStandard.
  floquet_omega_ref = periodic.floquet_reference_freq;
  MFEM_VERIFY(floquet_omega_ref == 0.0 || problem_type == ProblemType::DRIVEN,
              "FloquetReferenceFrequency (frequency-dependent k_F) is only supported for "
              "driven simulations!");

  // Get mesh dimensions in x/y/z coordinates.
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  bbmax -= bbmin;

  // BZ wrapping: wrap Floquet wave vector to the first Brillouin zone [-π/L, π/L].
  // DISABLED when frequency scaling is active: the BZ offset is frequency-independent
  // but kF scales with ω, creating a mismatch at frequencies away from the reference.
  // Without wrapping, both volume and Floquet port use the unwrapped kF consistently.
  if (floquet_omega_ref == 0.0)
  {
    for (int i = 0; i < sdim; i++)
    {
      double half_bz = M_PI / bbmax[i];
      if (wave_vector[i] > half_bz || wave_vector[i] < -half_bz)
      {
        wave_vector[i] = std::remainder(wave_vector[i], 2.0 * half_bz);
      }
    }
  }

  // Save BZ-wrapped k_F, then convert to k₀ = k_F/ω for frequency-independent tensors.
  wave_vector_bz = wave_vector;
  if (floquet_omega_ref > 0.0)
  {
    wave_vector *= 1.0 / floquet_omega_ref;
  }

  // Matrix representation of cross product with wave vector (or k₀ when scaling active).
  // [k x] = | 0  -k3  k2|
  //         | k3  0  -k1|
  //         |-k2  k1  0 |
  // For 2D, the cross product is not defined; set to sdim x sdim zeros so that
  // subsequent Mult operations have matching dimensions.
  wave_vector_cross.SetSize(sdim);
  wave_vector_cross = 0.0;
  if (sdim == 3)
  {
    wave_vector_cross(0, 1) = -wave_vector[2];
    wave_vector_cross(0, 2) = wave_vector[1];
    wave_vector_cross(1, 0) = wave_vector[2];
    wave_vector_cross(1, 2) = -wave_vector[0];
    wave_vector_cross(2, 0) = -wave_vector[1];
    wave_vector_cross(2, 1) = wave_vector[0];
  }
}

mfem::Array<int> MaterialOperator::GetBdrAttributeToMaterial() const
{
  // Construct map from all (contiguous) local libCEED boundary attributes to the material
  // index in the neighboring element.
  mfem::Array<int> bdr_attr_mat(mesh.MaxCeedBdrAttribute());
  bdr_attr_mat = -1;
  for (const auto &[attr, bdr_attr_map] : mesh.GetCeedBdrAttributes())
  {
    for (auto it = bdr_attr_map.begin(); it != bdr_attr_map.end(); ++it)
    {
      MFEM_ASSERT(it->second > 0 && it->second <= bdr_attr_mat.Size(),
                  "Invalid libCEED boundary attribute " << it->second << "!");
      bdr_attr_mat[it->second - 1] = AttrToMat(it->first);
    }
  }
  return bdr_attr_mat;
}

void RotateMaterialDefinitions(std::vector<config::MaterialData> &materials,
                               const mfem::Vector &e1, const mfem::Vector &e2,
                               const mfem::Vector &normal)
{
  // Build the 3x2 rotation R = [e1 | e2] mapping local 2D → global 3D.
  mfem::DenseMatrix R(3, 2);
  for (int d = 0; d < 3; d++)
  {
    R(d, 0) = e1(d);
    R(d, 1) = e2(d);
  }

  // Rotate a 3x3 symmetric M into the local frame and encode the result as a
  // SymmetricMatrixData<3>. The rotated tensor is block-diagonal in the local frame:
  //   leading 2x2 = R^T M R   (in-plane)
  //   (2, 2)     = n^T M n    (out-of-plane)
  //   cross terms 0.
  // We represent the 2x2 block via its 2x2 eigendecomposition (two unit eigenvectors
  // lifted to 3D with a trailing 0) and the scalar out-of-plane via (0, 0, 1).
  auto Rotated = [&](const mfem::DenseMatrix &M3)
  {
    mfem::DenseMatrix temp(3, 2), A(2, 2);
    Mult(M3, R, temp);
    MultAtB(R, temp, A);
    const double a = A(0, 0), b = 0.5 * (A(0, 1) + A(1, 0)), d = A(1, 1);
    const double tr = a + d;
    const double disc = std::sqrt(std::max(0.25 * (a - d) * (a - d) + b * b, 0.0));
    double s1 = 0.5 * tr + disc, s2 = 0.5 * tr - disc;
    std::array<double, 2> w1, w2;
    const double scale = std::abs(a) + std::abs(d) + 1.0;
    if (std::abs(b) > 1.0e-14 * scale)
    {
      // Eigenvector for s1: proportional to (b, s1 - a) (or (s1 - d, b)).
      double nx = s1 - d, ny = b;
      const double nrm = std::sqrt(nx * nx + ny * ny);
      w1 = {nx / nrm, ny / nrm};
      w2 = {-w1[1], w1[0]};
    }
    else
    {
      // Already diagonal; eigenvalues are a and d themselves.
      s1 = a;
      s2 = d;
      w1 = {1.0, 0.0};
      w2 = {0.0, 1.0};
    }
    const double sn = M3.InnerProduct(normal, normal);
    config::SymmetricMatrixData<3> out(0.0);
    out.s = {s1, s2, sn};
    out.v = {{{w1[0], w1[1], 0.0}, {w2[0], w2[1], 0.0}, {0.0, 0.0, 1.0}}};
    return out;
  };

  for (auto &data : materials)
  {
    const mfem::DenseMatrix mu_3d = internal::mat::ToDenseMatrix(data.mu_r);
    const mfem::DenseMatrix eps_3d = internal::mat::ToDenseMatrix(data.epsilon_r);
    const mfem::DenseMatrix tandelta_3d = internal::mat::ToDenseMatrix(data.tandelta);
    const mfem::DenseMatrix sigma_3d = internal::mat::ToDenseMatrix(data.sigma);

    data.mu_r = Rotated(mu_3d);
    data.epsilon_r = Rotated(eps_3d);
    data.tandelta = Rotated(tandelta_3d);
    data.sigma = Rotated(sigma_3d);
    // `lambda_L` is a scalar length; unchanged by rotation.
  }
}

MaterialPropertyCoefficient::MaterialPropertyCoefficient(int attr_max)
{
  attr_mat.SetSize(attr_max);
  attr_mat = -1;
}

MaterialPropertyCoefficient::MaterialPropertyCoefficient(
    const mfem::Array<int> &attr_mat_, const mfem::DenseTensor &mat_coeff_, double a)
  : attr_mat(attr_mat_), mat_coeff(mat_coeff_)
{
  *this *= a;
}

namespace
{

void UpdateProperty(mfem::DenseTensor &mat_coeff, int k, double coeff, double a)
{
  // Constant diagonal coefficient.
  if (mat_coeff.SizeI() == 0 && mat_coeff.SizeJ() == 0)
  {
    // Initialize the coefficient material properties.
    MFEM_VERIFY(k == 0 && mat_coeff.SizeK() == 1,
                "Unexpected initial size for MaterialPropertyCoefficient!");
    mat_coeff.SetSize(1, 1, mat_coeff.SizeK());
    mat_coeff(0, 0, k) = a * coeff;
  }
  else
  {
    MFEM_VERIFY(mat_coeff.SizeI() == mat_coeff.SizeJ(),
                "Invalid dimensions for MaterialPropertyCoefficient update!");
    for (int i = 0; i < mat_coeff.SizeI(); i++)
    {
      mat_coeff(i, i, k) += a * coeff;
    }
  }
}

void UpdateProperty(mfem::DenseTensor &mat_coeff, int k, const mfem::DenseMatrix &coeff,
                    double a)
{
  if (mat_coeff.SizeI() == 0 && mat_coeff.SizeJ() == 0)
  {
    // Initialize the coefficient material properties.
    MFEM_VERIFY(k == 0 && mat_coeff.SizeK() == 1,
                "Unexpected initial size for MaterialPropertyCoefficient!");
    mat_coeff.SetSize(coeff.Height(), coeff.Width(), mat_coeff.SizeK());
    mat_coeff(k).Set(a, coeff);
  }
  else if (coeff.Height() == mat_coeff.SizeI() && coeff.Width() == mat_coeff.SizeJ())
  {
    // Add as full matrix.
    mat_coeff(k).Add(a, coeff);
  }
  else if (coeff.Height() == 1 && coeff.Width() == 1)
  {
    // Add as diagonal.
    UpdateProperty(mat_coeff, k, coeff(0, 0), a);
  }
  else if (mat_coeff.SizeI() == 1 && mat_coeff.SizeJ() == 1)
  {
    // Convert to matrix coefficient and previous data add as diagonal.
    mfem::DenseTensor mat_coeff_scalar(mat_coeff);
    mat_coeff.SetSize(coeff.Height(), coeff.Width(), mat_coeff_scalar.SizeK());
    mat_coeff = 0.0;
    for (int l = 0; l < mat_coeff.SizeK(); l++)
    {
      UpdateProperty(mat_coeff, l, mat_coeff_scalar(0, 0, l), 1.0);
    }
    mat_coeff(k).Add(a, coeff);
  }
  else
  {
    MFEM_ABORT("Invalid dimensions when updating material property at index " << k << "!");
  }
}

bool Equals(const mfem::DenseMatrix &mat_coeff, double coeff, double a)
{
  MFEM_VERIFY(mat_coeff.Height() == mat_coeff.Width(),
              "Invalid dimensions for MaterialPropertyCoefficient update!");
  constexpr double tol = 1.0e-9;
  for (int i = 0; i < mat_coeff.Height(); i++)
  {
    if (std::abs(mat_coeff(i, i) - a * coeff) >= tol * std::abs(mat_coeff(i, i)))
    {
      return false;
    }
    for (int j = 0; j < mat_coeff.Width(); j++)
    {
      if (j != i && std::abs(mat_coeff(i, j)) > 0.0)
      {
        return false;
      }
    }
  }
  return true;
}

bool Equals(const mfem::DenseMatrix &mat_coeff, const mfem::DenseMatrix &coeff, double a)
{
  if (coeff.Height() == 1 && coeff.Width() == 1)
  {
    return Equals(mat_coeff, coeff(0, 0), a);
  }
  else
  {
    constexpr double tol = 1.0e-9;
    mfem::DenseMatrix T(mat_coeff);
    T.Add(-a, coeff);
    return (T.MaxMaxNorm() < tol * mat_coeff.MaxMaxNorm());
  }
}

}  // namespace

void MaterialPropertyCoefficient::AddCoefficient(const mfem::Array<int> &attr_mat_,
                                                 const mfem::DenseTensor &mat_coeff_,
                                                 double a)
{
  if (empty())
  {
    MFEM_VERIFY(attr_mat_.Size() == attr_mat.Size(),
                "Invalid resize of attribute to material property map in "
                "MaterialPropertyCoefficient::AddCoefficient!");
    attr_mat = attr_mat_;
    mat_coeff = mat_coeff_;
    *this *= a;
  }
  else if (attr_mat_ == attr_mat)
  {
    MFEM_VERIFY(mat_coeff_.SizeK() == mat_coeff.SizeK(),
                "Invalid dimensions for MaterialPropertyCoefficient::AddCoefficient!");
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      UpdateProperty(mat_coeff, k, mat_coeff_(k), a);
    }
  }
  else
  {
    for (int k = 0; k < mat_coeff_.SizeK(); k++)
    {
      // Get list of all attributes which use this material property.
      mfem::Array<int> attr_list;
      attr_list.Reserve(attr_mat_.Size());
      for (int i = 0; i < attr_mat_.Size(); i++)
      {
        if (attr_mat_[i] == k)
        {
          attr_list.Append(i + 1);
        }
      }

      // Add or update the material property.
      AddMaterialProperty(attr_list, mat_coeff_(k), a);
    }
  }
}

template <typename T>
void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &attr_list,
                                                      const T &coeff, double a)
{
  // Preprocess the attribute list. If any of the given attributes already have material
  // properties assigned, then they all need to point to the same material and it is
  // updated in place. Otherwise a new material is added for these attributes.
  if (attr_list.Size() == 0)
  {
    // No attributes, nothing to add.
    return;
  }

  int mat_idx = -1;
  for (auto attr : attr_list)
  {
    MFEM_VERIFY(attr <= attr_mat.Size(),
                "Out of bounds access for attribute "
                    << attr << " in MaterialPropertyCoefficient::AddMaterialProperty!");
    if (mat_idx < 0)
    {
      mat_idx = attr_mat[attr - 1];
    }
    else
    {
      MFEM_VERIFY(mat_idx == attr_mat[attr - 1],
                  "All attributes for MaterialPropertyCoefficient::AddMaterialProperty "
                  "must correspond to the same "
                  "existing material if it exists!");
    }
  }

  if (mat_idx < 0)
  {
    // Check if we can reuse an existing material.
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      if (Equals(mat_coeff(k), coeff, a))
      {
        mat_idx = k;
        break;
      }
    }
    if (mat_idx < 0)
    {
      // Append a new material and assign the attributes to it.
      const mfem::DenseTensor mat_coeff_backup(mat_coeff);
      mat_coeff.SetSize(mat_coeff_backup.SizeI(), mat_coeff_backup.SizeJ(),
                        mat_coeff_backup.SizeK() + 1);
      for (int k = 0; k < mat_coeff_backup.SizeK(); k++)
      {
        mat_coeff(k).Set(1.0, mat_coeff_backup(k));
      }
      mat_idx = mat_coeff.SizeK() - 1;
    }
    mat_coeff(mat_idx) = 0.0;  // Zero out so we can add

    // Assign all attributes to this new material.
    for (auto attr : attr_list)
    {
      attr_mat[attr - 1] = mat_idx;
    }
  }
  UpdateProperty(mat_coeff, mat_idx, coeff, a);
}

MaterialPropertyCoefficient &MaterialPropertyCoefficient::operator*=(double a)
{
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    mat_coeff(k) *= a;
  }
  return *this;
}

void MaterialPropertyCoefficient::RestrictCoefficient(const mfem::Array<int> &attr_list)
{
  // Create a new material property coefficient with materials corresponding to only the
  // unique ones in the given attribute list.
  const mfem::Array<int> attr_mat_orig(attr_mat);
  const mfem::DenseTensor mat_coeff_orig(mat_coeff);
  attr_mat = -1;
  mat_coeff.SetSize(mat_coeff_orig.SizeI(), mat_coeff_orig.SizeJ(), 0);
  for (auto attr : attr_list)
  {
    if (attr_mat[attr - 1] >= 0)
    {
      // Attribute has already been processed.
      continue;
    }

    // Find all attributes in restricted list of attributes which map to this material index
    // and process them together.
    const int orig_mat_idx = attr_mat_orig[attr - 1];
    const int new_mat_idx = mat_coeff.SizeK();
    for (auto attr2 : attr_list)
    {
      if (attr_mat_orig[attr2 - 1] == orig_mat_idx)
      {
        attr_mat[attr2 - 1] = new_mat_idx;
      }
    }

    // Append the new material property.
    const mfem::DenseTensor mat_coeff_backup(mat_coeff);
    mat_coeff.SetSize(mat_coeff_backup.SizeI(), mat_coeff_backup.SizeJ(),
                      mat_coeff_backup.SizeK() + 1);
    for (int k = 0; k < mat_coeff_backup.SizeK(); k++)
    {
      mat_coeff(k).Set(1.0, mat_coeff_backup(k));
    }
    mat_coeff(new_mat_idx).Set(1.0, mat_coeff_orig(orig_mat_idx));
  }
}

void MaterialPropertyCoefficient::NormalProjectedCoefficient(const mfem::Vector &normal)
{
  mfem::DenseTensor mat_coeff_backup(mat_coeff);
  mat_coeff.SetSize(1, 1, mat_coeff_backup.SizeK());
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    mat_coeff(k) = mat_coeff_backup(k).InnerProduct(normal, normal);
  }
}

template void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &,
                                                               const mfem::DenseMatrix &,
                                                               double);
template void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &,
                                                               const double &, double);

// Explicit template instantiations for internal::mat functions.
template bool internal::mat::IsOrthonormal(const config::SymmetricMatrixData<3> &);
template bool internal::mat::IsValid(const config::SymmetricMatrixData<3> &);
template bool internal::mat::IsIsotropic(const config::SymmetricMatrixData<3> &);
template bool internal::mat::IsIdentity(const config::SymmetricMatrixData<3> &);

}  // namespace palace
