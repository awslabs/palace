// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "pml.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "fem/qfunctions/coeff/pml_qf.h"
#include "linalg/densematrix.hpp"

namespace palace::pml
{

LayerGeometry DetectLayerGeometry(const std::array<double, 3> &inner_min,
                                  const std::array<double, 3> &inner_max,
                                  const std::array<double, 3> &outer_min,
                                  const std::array<double, 3> &outer_max, double rel_tol)
{
  double extent = 0.0;
  for (int a = 0; a < 3; a++)
  {
    extent = std::max(extent, outer_max[a] - outer_min[a]);
  }
  const double tol = rel_tol * extent;
  LayerGeometry geometry;
  for (int a = 0; a < 3; a++)
  {
    const double t_neg = inner_min[a] - outer_min[a], t_pos = outer_max[a] - inner_max[a];
    geometry.inner[2 * a] = inner_min[a];
    geometry.inner[2 * a + 1] = inner_max[a];
    geometry.thickness[2 * a] = (t_neg > tol) ? t_neg : 0.0;
    geometry.thickness[2 * a + 1] = (t_pos > tol) ? t_pos : 0.0;
  }
  return geometry;
}

LayerGeometry ConfiguredLayerGeometry(const config::PMLData &data,
                                      const std::array<double, 3> &outer_min,
                                      const std::array<double, 3> &outer_max)
{
  LayerGeometry geometry;
  for (int f = 0; f < 6; f++)
  {
    const int a = f / 2;
    const double t =
        (data.direction_signs[f] != 0) ? std::max(data.thickness[f], 0.0) : 0.0;
    geometry.thickness[f] = t;
    geometry.inner[f] = (f % 2 == 0) ? outer_min[a] + t : outer_max[a] - t;
  }
  return geometry;
}

double RefractiveIndex(const std::array<double, 9> &mu_inv,
                       const std::array<double, 9> &epsilon_real)
{
  // Smallest refractive index n = √(λ_min(ε) λ_min(μ)) over the principal directions of the
  // (symmetric positive definite) background tensors, with λ_min(μ) = 1 / λ_max(μ⁻¹). The
  // eigenvalues are the singular values.
  auto ToDenseMatrix = [](const std::array<double, 9> &A)
  {
    mfem::DenseMatrix M(3, 3);
    std::copy(A.begin(), A.end(), M.Data());
    return M;
  };
  const double eps_min = linalg::SingularValueMin(ToDenseMatrix(epsilon_real));
  const double mu_inv_max = linalg::SingularValueMax(ToDenseMatrix(mu_inv));
  MFEM_VERIFY(eps_min > 0.0 && mu_inv_max > 0.0,
              "Invalid background material properties for PML!");
  return std::sqrt(eps_min / mu_inv_max);
}

std::array<double, 6> ResolveSigmaMax(const config::PMLData &data,
                                      const LayerGeometry &geometry,
                                      const std::array<double, 6> &n_r)
{
  std::array<double, 6> sigma_max{};
  for (int f = 0; f < 6; f++)
  {
    const double t = geometry.thickness[f];
    if (t <= 0.0)
    {
      continue;
    }
    sigma_max[f] =
        (data.sigma_max[f / 2] >= 0.0)
            ? data.sigma_max[f / 2]
            : -(data.order + 1) * std::log(data.reflection_target) / (2.0 * t * n_r[f]);
  }
  return sigma_max;
}

Profile BuildProfile(const config::PMLData &data, const LayerGeometry &geometry,
                     const std::array<double, 6> &n_r, const std::array<double, 9> &mu_inv,
                     const std::array<double, 9> &epsilon_real,
                     const std::array<double, 9> &epsilon_imag)
{
  Profile profile;
  profile.geometry = geometry;
  profile.sigma_max = ResolveSigmaMax(data, geometry, n_r);
  profile.kappa_max = data.kappa_max;
  profile.alpha_max = data.alpha_max;
  profile.order = data.order;
  profile.frequency_dependent = data.frequency_dependent;
  profile.reference_frequency = data.frequency_dependent ? 0.0 : data.reference_frequency;
  profile.mu_inv = mu_inv;
  profile.epsilon_real = epsilon_real;
  profile.epsilon_imag = epsilon_imag;
  return profile;
}

std::string CheckStretchConsistency(const Profile &p, const Profile &q, double tol)
{
  // The stretch factors of all profiles must be the same functions of position on the faces
  // they share: the PML is only reflectionless if it is a coordinate transformation, also
  // across material interfaces inside the layer.
  auto Differ = [](double a, double b, double tol)
  { return std::abs(a - b) > tol * std::max({std::abs(a), std::abs(b), 1.0e-300}); };
  if (p.frequency_dependent != q.frequency_dependent ||
      Differ(p.reference_frequency, q.reference_frequency, 1.0e-12))
  {
    return "\"FrequencyDependent\" and \"ReferenceFrequency\"";
  }
  for (int f = 0; f < 6; f++)
  {
    if (p.geometry.thickness[f] <= 0.0 || q.geometry.thickness[f] <= 0.0)
    {
      continue;
    }
    const int a = f / 2;
    if (std::abs(p.geometry.inner[f] - q.geometry.inner[f]) > tol ||
        std::abs(p.geometry.thickness[f] - q.geometry.thickness[f]) > tol)
    {
      return "\"Thickness\"";
    }
    if (p.order != q.order)
    {
      return "\"Order\"";
    }
    if (Differ(p.sigma_max[f], q.sigma_max[f], 1.0e-12))
    {
      return "\"SigmaMax\" and \"ReflectionTarget\"";
    }
    if (Differ(p.kappa_max[a], q.kappa_max[a], 1.0e-12) ||
        Differ(p.alpha_max[a], q.alpha_max[a], 1.0e-12))
    {
      return "\"KappaMax\" and \"AlphaMax\"";
    }
  }
  return {};
}

std::array<double, 3> ComputeDepthFraction(const Profile &profile,
                                           const std::array<double, 3> &x)
{
  const auto &g = profile.geometry;
  std::array<double, 3> r{};
  for (int a = 0; a < 3; a++)
  {
    if (g.thickness[2 * a] > 0.0 && x[a] < g.inner[2 * a])
    {
      r[a] = (g.inner[2 * a] - x[a]) / g.thickness[2 * a];
    }
    else if (g.thickness[2 * a + 1] > 0.0 && x[a] > g.inner[2 * a + 1])
    {
      r[a] = (x[a] - g.inner[2 * a + 1]) / g.thickness[2 * a + 1];
    }
    r[a] = std::min(r[a], 1.0);
  }
  return r;
}

std::vector<CeedIntScalar> PackContext(const ContextHeader &header,
                                       const std::vector<int> &attr_to_profile,
                                       const std::vector<Profile> &profiles,
                                       bool frequency_dependent)
{
  // Compact the profile list to the included profiles.
  std::vector<int> profile_map(profiles.size(), -1);
  int num_profiles = 0;
  for (std::size_t k = 0; k < profiles.size(); k++)
  {
    if (profiles[k].frequency_dependent == frequency_dependent)
    {
      profile_map[k] = num_profiles++;
    }
  }
  bool any = false;
  for (auto k : attr_to_profile)
  {
    any = any || (k >= 0 && profile_map[k] >= 0);
  }
  if (!any)
  {
    return {};
  }

  const int num_attr = static_cast<int>(attr_to_profile.size());
  std::vector<CeedIntScalar> ctx(PALACE_PML_HEADER_SIZE + num_attr +
                                 PALACE_PML_PROFILE_SIZE * num_profiles);
  ctx[0].first = num_attr;
  ctx[1].first = static_cast<int>(header.part);
  ctx[2].second = header.c_muinv.real();
  ctx[3].second = header.c_muinv.imag();
  ctx[4].second = header.c_eps.real();
  ctx[5].second = header.c_eps.imag();
  ctx[6].second = header.omega.real();
  ctx[7].second = header.omega.imag();
  for (int i = 0; i < 9; i++)
  {
    ctx[8 + i].second = header.wave_vector_cross[i];
  }
  for (int i = 0; i < num_attr; i++)
  {
    const int k = attr_to_profile[i];
    ctx[PALACE_PML_HEADER_SIZE + i].first = (k >= 0) ? profile_map[k] : -1;
  }
  for (std::size_t k = 0; k < profiles.size(); k++)
  {
    if (profile_map[k] < 0)
    {
      continue;
    }
    const auto &p = profiles[k];
    CeedIntScalar *out = ctx.data() + PALACE_PML_HEADER_SIZE + num_attr +
                         PALACE_PML_PROFILE_SIZE * profile_map[k];
    out[0].first = p.frequency_dependent ? 1 : 0;
    out[1].first = p.order;
    out[2].second = p.reference_frequency;
    for (int f = 0; f < 6; f++)
    {
      out[3 + f].second = p.geometry.inner[f];
      out[9 + f].second = p.geometry.thickness[f];
      out[15 + f].second = p.sigma_max[f];
    }
    for (int a = 0; a < 3; a++)
    {
      out[21 + a].second = p.kappa_max[a];
      out[24 + a].second = p.alpha_max[a];
    }
    for (int i = 0; i < 9; i++)
    {
      out[27 + i].second = p.mu_inv[i];
      out[36 + i].second = p.epsilon_real[i];
      out[45 + i].second = p.epsilon_imag[i];
    }
  }
  return ctx;
}

}  // namespace palace::pml
