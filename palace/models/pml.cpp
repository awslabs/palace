// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "pml.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "fem/qfunctions/coeff/pml_qf.h"
#include "linalg/densematrix.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"

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
    const double t = data.directions[f] ? std::max(data.thickness[f], 0.0) : 0.0;
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

Stretch BuildStretch(const config::PMLData &data, const LayerGeometry &geometry, double n_r)
{
  Stretch stretch;
  stretch.geometry = geometry;
  for (int f = 0; f < 6; f++)
  {
    const double t = geometry.thickness[f];
    if (t <= 0.0)
    {
      continue;
    }
    stretch.sigma_max[f] =
        (data.sigma_max[f / 2] >= 0.0)
            ? data.sigma_max[f / 2]
            : -(data.order + 1) * std::log(data.reflection_target) / (2.0 * t * n_r);
  }
  stretch.kappa_max = data.kappa_max;
  stretch.alpha_max = data.alpha_max;
  stretch.order = data.order;
  stretch.frequency_dependent = data.frequency_dependent;
  stretch.reference_frequency = data.frequency_dependent ? 0.0 : data.reference_frequency;
  return stretch;
}

namespace
{

std::string FaceName(int f)
{
  return fmt::format("{}{}", (f % 2) ? "+" : "-", "xyz"[f / 2]);
}

}  // namespace

Layer::Layer(const config::PMLData &data,
             const std::vector<config::MaterialData> &materials, const Mesh &mesh,
             const mfem::Array<int> &attr_mat, const mfem::DenseTensor &mu_inv,
             const mfem::DenseTensor &epsilon_real, const mfem::DenseTensor &epsilon_imag)
  : attributes(data.attributes)
{
  const mfem::ParMesh &pmesh = mesh.Get();
  std::ranges::sort(attributes);
  attributes.erase(std::ranges::unique(attributes).begin(), attributes.end());
  MFEM_VERIFY(!attributes.empty(), "\"PML.Attributes\" must not be empty!");
  MFEM_VERIFY(data.frequency_dependent || data.reference_frequency > 0.0,
              "Static PML requires a positive reference frequency!");
  auto IsPMLAttribute = [this](int attr)
  { return std::ranges::binary_search(attributes, attr); };

  // Smallest refractive index of the materials of the PML regions, for the default σ_max
  // (the same on all processes).
  double n_r = std::numeric_limits<double>::infinity();
  for (const auto &material : materials)
  {
    if (std::ranges::none_of(material.attributes, IsPMLAttribute))
    {
      continue;
    }
    MFEM_VERIFY(!internal::mat::IsValid(material.sigma) && material.lambda_L == 0.0,
                "PML regions do not support materials with electrical conductivity or "
                "London penetration depth (attributes "
                    << fmt::format("{}", fmt::join(material.attributes, ", ")) << ")!");
    mfem::DenseMatrix muinv(3, 3);
    mfem::DenseMatrixInverse(internal::mat::ToDenseMatrix(material.mu_r), true)
        .GetInverseMatrix(muinv);
    const auto eps = internal::mat::ToDenseMatrix(material.epsilon_r);
    std::array<double, 9> muinv_data, eps_data;
    std::copy_n(muinv.Data(), 9, muinv_data.begin());
    std::copy_n(eps.Data(), 9, eps_data.begin());
    n_r = std::min(n_r, RefractiveIndex(muinv_data, eps_data));
  }
  MFEM_VERIFY(std::isfinite(n_r), "No material found for the PML regions!");

  // Bounding boxes of the whole mesh and of the non-PML (physical) region, whose faces are
  // the outer PML boundaries and the inner PML interfaces, respectively. These are global
  // reductions, called on all ranks.
  int attr_max = pmesh.attributes.Size() ? pmesh.attributes.Max() : 0;
  Mpi::GlobalMax(1, &attr_max, pmesh.GetComm());
  mfem::Array<int> physical_marker(attr_max);
  physical_marker = 1;
  for (auto attr : attributes)
  {
    if (attr <= attr_max)
    {
      physical_marker[attr - 1] = 0;
    }
  }
  mfem::Vector bbmin, bbmax, phys_bbmin, phys_bbmax;
  mesh::GetAxisAlignedBoundingBox(pmesh, bbmin, bbmax);
  mesh::GetAxisAlignedBoundingBox(pmesh, physical_marker, false, phys_bbmin, phys_bbmax);
  MFEM_VERIFY(phys_bbmin(0) <= phys_bbmax(0),
              "The PML regions require a physical (non-PML) region in the mesh!");
  const std::array<double, 3> outer_min{{bbmin(0), bbmin(1), bbmin(2)}},
      outer_max{{bbmax(0), bbmax(1), bbmax(2)}},
      inner_min{{phys_bbmin(0), phys_bbmin(1), phys_bbmin(2)}},
      inner_max{{phys_bbmax(0), phys_bbmax(1), phys_bbmax(2)}};
  double extent = 0.0;
  for (int a = 0; a < 3; a++)
  {
    extent = std::max(extent, outer_max[a] - outer_min[a]);
  }
  const double tol = 1.0e-6 * extent;

  // Layer geometry and stretch.
  const auto geometry =
      data.autodetect_geometry
          ? DetectLayerGeometry(inner_min, inner_max, outer_min, outer_max)
          : ConfiguredLayerGeometry(data, outer_min, outer_max);
  MFEM_VERIFY(std::ranges::any_of(geometry.thickness, [](double t) { return t > 0.0; }),
              "No active PML faces found: "
                  << (data.autodetect_geometry
                          ? "the PML regions must lie outside of the bounding box of the "
                            "non-PML regions of the mesh (or specify \"PML.Direction\" "
                            "and \"PML.Thickness\")!"
                          : "\"PML.Thickness\" must be positive for at least one "
                            "\"PML.Direction\"!"));
  stretch = BuildStretch(data, geometry, n_r);

  // The non-PML regions must lie inside of the layer.
  for (int f = 0; f < 6; f++)
  {
    const bool outside = geometry.thickness[f] > 0.0 &&
                         ((f % 2 == 0) ? phys_bbmin(f / 2) < geometry.inner[f] - tol
                                       : phys_bbmax(f / 2) > geometry.inner[f] + tol);
    MFEM_VERIFY(!outside, "Non-PML regions of the mesh extend into the PML layer on the "
                              << FaceName(f)
                              << " face: their attributes must be PML regions as well!");
  }

  // Every element of a PML region must be in the layer (an element of the physical region
  // would have no material properties, and no stretch would be applied in a region beyond
  // an inactive face).
  {
    int outside = 0;
    for (int e = 0; e < pmesh.GetNE() && !outside; e++)
    {
      if (!IsPMLAttribute(pmesh.GetAttribute(e)))
      {
        continue;
      }
      mfem::IsoparametricTransformation T;
      mfem::Vector center;
      pmesh.GetElementTransformation(e, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), center);
      bool in_layer = false;
      for (int f = 0; f < 6; f++)
      {
        in_layer = in_layer || (geometry.thickness[f] > 0.0 &&
                                ((f % 2 == 0) ? center(f / 2) < geometry.inner[f]
                                              : center(f / 2) > geometry.inner[f]));
      }
      outside = in_layer ? 0 : pmesh.GetAttribute(e);
    }
    Mpi::GlobalMax(1, &outside, pmesh.GetComm());
    MFEM_VERIFY(!outside, "PML attribute "
                              << outside
                              << " has elements outside of the PML layer, where no PML "
                                 "stretch is applied (check \"PML.Direction\" and "
                                 "\"PML.Thickness\")!");
  }

  // Background material of each local libCEED attribute in a PML region.
  MFEM_VERIFY(mu_inv.SizeI() == 3 && epsilon_real.SizeI() == 3 && epsilon_imag.SizeI() == 3,
              "PML regions require 3D material tensors!");
  const auto &loc_attr = mesh.GetCeedAttributes();
  attr_background.assign(attr_mat.Size(), -1);
  std::unordered_map<int, int> mat_background;
  for (auto attr : attributes)
  {
    const auto it = loc_attr.find(attr);
    if (it == loc_attr.end())
    {
      continue;
    }
    const int mat = attr_mat[it->second - 1];
    MFEM_VERIFY(mat >= 0, "Missing material properties for PML attribute " << attr << "!");
    auto [bg, inserted] =
        mat_background.try_emplace(mat, static_cast<int>(backgrounds.size()));
    if (inserted)
    {
      auto &b = backgrounds.emplace_back();
      std::copy_n(mu_inv(mat).Data(), 9, b.mu_inv.begin());
      std::copy_n(epsilon_real(mat).Data(), 9, b.epsilon_real.begin());
      std::copy_n(epsilon_imag(mat).Data(), 9, b.epsilon_imag.begin());
    }
    attr_background[it->second - 1] = bg->second;
  }
}

std::vector<CeedIntScalar> Layer::PackContext(const ContextHeader &header) const
{
  return pml::PackContext(header, stretch, attr_background, backgrounds);
}

std::vector<CeedIntScalar> PackContext(const ContextHeader &header, const Stretch &stretch,
                                       const std::vector<int> &attr_background,
                                       const std::vector<Background> &backgrounds)
{
  if (std::ranges::none_of(attr_background, [](int k) { return k >= 0; }))
  {
    return {};
  }
  const int num_attr = static_cast<int>(attr_background.size());
  std::vector<CeedIntScalar> ctx(PALACE_PML_HEADER_SIZE + PALACE_PML_STRETCH_SIZE +
                                 num_attr +
                                 PALACE_PML_BACKGROUND_SIZE * backgrounds.size());
  const std::complex<double> omega =
      stretch.frequency_dependent ? header.omega : stretch.reference_frequency;
  ctx[0].first = num_attr;
  ctx[1].first = static_cast<int>(header.part);
  ctx[2].second = header.c_muinv.real();
  ctx[3].second = header.c_muinv.imag();
  ctx[4].second = header.c_eps.real();
  ctx[5].second = header.c_eps.imag();
  ctx[6].second = omega.real();
  ctx[7].second = omega.imag();
  for (int i = 0; i < 9; i++)
  {
    ctx[8 + i].second = header.wave_vector_cross[i];
  }
  CeedIntScalar *p = ctx.data() + PALACE_PML_HEADER_SIZE;
  p[0].first = stretch.order;
  for (int f = 0; f < 6; f++)
  {
    p[1 + f].second = stretch.geometry.inner[f];
    p[7 + f].second = stretch.geometry.thickness[f];
    p[13 + f].second = stretch.sigma_max[f];
  }
  for (int a = 0; a < 3; a++)
  {
    p[19 + a].second = stretch.kappa_max[a];
    p[22 + a].second = stretch.alpha_max[a];
  }
  CeedIntScalar *attr_map = p + PALACE_PML_STRETCH_SIZE;
  for (int i = 0; i < num_attr; i++)
  {
    MFEM_ASSERT(attr_background[i] < static_cast<int>(backgrounds.size()),
                "Invalid PML background material index!");
    attr_map[i].first = attr_background[i];
  }
  CeedIntScalar *b = attr_map + num_attr;
  for (const auto &background : backgrounds)
  {
    for (int i = 0; i < 9; i++)
    {
      b[i].second = background.mu_inv[i];
      b[9 + i].second = background.epsilon_real[i];
      b[18 + i].second = background.epsilon_imag[i];
    }
    b += PALACE_PML_BACKGROUND_SIZE;
  }
  return ctx;
}

}  // namespace palace::pml
