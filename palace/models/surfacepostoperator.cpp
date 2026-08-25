// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <set>
#include <tuple>
#include <utility>
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/strattonchu.hpp"
#include "utils/communication.hpp"
#include "utils/edgedistance.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/metaledge.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

template <typename T>
mfem::Array<int> SetUpBoundaryProperties(const T &data,
                                         const mfem::Array<int> &bdr_attr_marker)
{
  mfem::Array<int> attr_list;
  attr_list.Reserve(static_cast<int>(data.attributes.size()));
  std::set<int> bdr_warn_list;
  for (auto attr : data.attributes)
  {
    // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
    //             "Boundary postprocessing attribute tags must be non-negative and "
    //             "correspond to attributes in the mesh!");
    // MFEM_VERIFY(bdr_attr_marker[attr - 1],
    //             "Unknown boundary postprocessing attribute " << attr << "!");
    if (attr <= 0 || attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
    {
      bdr_warn_list.insert(attr);
    }
    else
    {
      attr_list.Append(attr);
    }
  }
  if (!bdr_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning(
        "Unknown boundary postprocessing attributes!\nSolver will just ignore them!");
    utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  return attr_list;
}

double EdgeDistanceOutsideWeight(double distance, double radius, double smoothing)
{
  if (!std::isfinite(radius))
  {
    return 0.0;
  }
  if (smoothing == 0.0)
  {
    return distance >= radius ? 1.0 : 0.0;
  }

  const double lower = radius * (1.0 - smoothing);
  const double upper = radius * (1.0 + smoothing);
  if (distance <= lower)
  {
    return 0.0;
  }
  if (distance >= upper)
  {
    return 1.0;
  }
  const double x = (distance - lower) / (upper - lower);
  return x * x * (3.0 - 2.0 * x);
}

double EdgeDistanceWindowWeight(double distance, double distance_min, double distance_max,
                                double smoothing)
{
  return std::max(0.0, EdgeDistanceOutsideWeight(distance, distance_min, smoothing) -
                           EdgeDistanceOutsideWeight(distance, distance_max, smoothing));
}

class EdgeDistanceCoefficient : public mfem::Coefficient
{
private:
  mfem::Coefficient &coefficient;
  const EdgeDistanceTree &edge_distance_tree;
  const double distance_min;
  const double distance_max;
  const double smoothing;

public:
  EdgeDistanceCoefficient(mfem::Coefficient &coefficient,
                          const EdgeDistanceTree &edge_distance_tree, double distance_min,
                          double distance_max, double smoothing)
    : coefficient(coefficient), edge_distance_tree(edge_distance_tree),
      distance_min(distance_min), distance_max(distance_max), smoothing(smoothing)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector point(T.GetSpaceDim());
    T.Transform(ip, point);
    const double distance = std::sqrt(edge_distance_tree.DistanceSquared(point));
    return EdgeDistanceWindowWeight(distance, distance_min, distance_max, smoothing) *
           coefficient.Eval(T, ip);
  }
};

}  // namespace

SurfacePostOperator::SurfaceFluxData::SurfaceFluxData(
    const config::SurfaceFluxData &data, const mfem::ParMesh &mesh,
    const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);

  // Store the type of flux.
  switch (data.type)
  {
    case SurfaceFlux::ELECTRIC:
      type = SurfaceFlux::ELECTRIC;
      break;
    case SurfaceFlux::MAGNETIC:
      type = SurfaceFlux::MAGNETIC;
      break;
    case SurfaceFlux::POWER:
      type = SurfaceFlux::POWER;
      break;
  }

  // Store information about the global direction for orientation. Note the true boundary
  // normal is used in calculating the flux, this is just used to determine the sign.
  two_sided = data.two_sided;
  if (!two_sided)
  {
    center.SetSize(mesh.SpaceDimension());
    if (data.no_center)
    {
      // Compute the center as the bounding box centroid for all boundary elements making up
      // this postprocessing boundary.
      mfem::Vector bbmin, bbmax;
      mesh::GetAxisAlignedBoundingBox(
          mesh, mesh::AttrToMarker(bdr_attr_marker.Size(), attr_list), true, bbmin, bbmax);
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        center(d) = 0.5 * (bbmin(d) + bbmax(d));
      }
    }
    else
    {
      std::copy(data.center.begin(), data.center.end(), center.begin());
    }
  }
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::SurfaceFluxData::GetCoefficient(const mfem::ParGridFunction *E,
                                                     const mfem::ParGridFunction *B,
                                                     const MaterialOperator &mat_op) const
{
  switch (type)
  {
    case SurfaceFlux::ELECTRIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>>(
          attr_list, E, nullptr, mat_op, two_sided, center);
    case SurfaceFlux::MAGNETIC:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>>>(
          attr_list, nullptr, B, mat_op, two_sided, center);
    case SurfaceFlux::POWER:
      return std::make_unique<
          RestrictedCoefficient<BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>>>(
          attr_list, E, B, mat_op, two_sided, center);
  }
  return {};
}

SurfacePostOperator::InterfaceDielectricData::InterfaceDielectricData(
    const config::InterfaceDielectricData &data, const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);
  edge_distances = data.edge_distances;
  edge_distance_smoothing = data.edge_distance_smoothing;
  localize_edge_energy = data.localize_edge_energy;
  save_local_edge_energy = data.save_local_edge_energy;
  edge_frame_normal = data.edge_frame_normal.value_or(std::array<double, 3>{});
  flux_recovery = data.flux_recovery;

  // Calculate surface dielectric loss according to the formulas from J. Wenner et al.,
  // Surface loss simulations of superconducting coplanar waveguide resonators, Appl. Phys.
  // Lett. (2011). If only a general layer permittivity is specified and not any special
  // metal-air (MA), metal-substrate (MS), or substrate-air (SA) permittivity, compute the
  // numerator of the participation ratio according to the regular formula
  //                       p * E_elec = 1/2 t Re{∫ (ε E)ᴴ E_m dS} .
  switch (data.type)
  {
    case InterfaceDielectric::DEFAULT:
      type = InterfaceDielectric::DEFAULT;
      break;
    case InterfaceDielectric::MA:
      type = InterfaceDielectric::MA;
      break;
    case InterfaceDielectric::MS:
      type = InterfaceDielectric::MS;
      break;
    case InterfaceDielectric::SA:
      type = InterfaceDielectric::SA;
      break;
  }
  t = data.t;
  epsilon = data.epsilon_r;
  tandelta = data.tandelta;
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    const GridFunction &E, const GridFunction *D, const MaterialOperator &mat_op,
    InterfaceDielectricComponent component) const
{
  const GridFunction *recovered_flux = flux_recovery ? D : nullptr;
  MFEM_VERIFY(!flux_recovery || recovered_flux,
              "Interface dielectric flux recovery was requested but no recovered electric "
              "flux was supplied!");
  switch (type)
  {
    case InterfaceDielectric::DEFAULT:
      return std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>>(
          attr_list, E, mat_op, t, epsilon, recovered_flux, component);
    case InterfaceDielectric::MA:
      return std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>>(
          attr_list, E, mat_op, t, epsilon, recovered_flux, component);
    case InterfaceDielectric::MS:
      return std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>>(
          attr_list, E, mat_op, t, epsilon, recovered_flux, component);
    case InterfaceDielectric::SA:
      return std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>>(
          attr_list, E, mat_op, t, epsilon, recovered_flux, component);
  }
  return {};  // For compiler warning
}

SurfacePostOperator::FarFieldData::FarFieldData(const config::FarFieldPostData &data,
                                                const mfem::ParMesh &mesh,
                                                const mfem::Array<int> &bdr_attr_marker)
  : thetaphis(data.thetaphis)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);
}

SurfacePostOperator::SurfacePostOperator(
    const config::BoundaryPostData &postpro, ProblemType problem_type,
    const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &h1_fespace,
    mfem::ParFiniteElementSpace &nd_fespace,
    const std::unordered_set<int> *cracked_attributes,
    const config::BoundaryData *boundaries,
    std::shared_ptr<const SurfacePostGeometry> *automatic_geometry)
  : mat_op(mat_op), h1_fespace(h1_fespace), nd_fespace(nd_fespace)
{
  BlockTimer setup_timer(Timer::CONSTRUCT_SURFACE_POST);

  // Check that boundary attributes have been specified correctly.
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!postpro.flux.empty() || !postpro.dielectric.empty() || !postpro.farfield.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
  }

  // Surface flux postprocessing.
  for (const auto &[idx, data] : postpro.flux)
  {
    MFEM_VERIFY(problem_type != ProblemType::ELECTROSTATIC ||
                    data.type == SurfaceFlux::ELECTRIC,
                "Magnetic field or power surface flux postprocessing are not available "
                "for electrostatic problems!");
    MFEM_VERIFY(problem_type != ProblemType::MAGNETOSTATIC ||
                    data.type == SurfaceFlux::MAGNETIC,
                "Electric field or power surface flux postprocessing are not available "
                "for magnetostatic problems!");
    flux_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh(), bdr_attr_marker);
  }

  // Interface dielectric postprocessing.
  MFEM_VERIFY(postpro.dielectric.empty() || problem_type != ProblemType::MAGNETOSTATIC,
              "Interface dielectric loss postprocessing is not available for "
              "magnetostatic problems!");
  using EdgeDistanceTreeKey =
      std::tuple<std::vector<int>, std::vector<int>, std::optional<std::array<double, 3>>>;
  std::map<EdgeDistanceTreeKey, std::shared_ptr<const EdgeDistanceTree>>
      edge_distance_trees;
  using AutomaticEdgeDistanceTreeKey =
      std::pair<std::vector<std::size_t>, std::optional<std::array<double, 3>>>;
  std::map<AutomaticEdgeDistanceTreeKey, std::shared_ptr<const EdgeDistanceTree>>
      automatic_edge_distance_trees;
  const SurfacePostGeometry *cached_automatic_geometry =
      automatic_geometry && *automatic_geometry ? automatic_geometry->get() : nullptr;
  std::shared_ptr<SurfacePostGeometry> new_automatic_geometry;
  if (automatic_geometry && !*automatic_geometry)
  {
    new_automatic_geometry = std::make_shared<SurfacePostGeometry>();
  }
  MetalEdgeGeometry metal_edges;
  const bool use_automatic_edges = std::any_of(
      postpro.dielectric.begin(), postpro.dielectric.end(),
      [cached_automatic_geometry](const auto &entry)
      {
        const auto &[idx, data] = entry;
        if (!data.automatic_edges)
        {
          return false;
        }
        return data.localize_edge_energy || !cached_automatic_geometry ||
               cached_automatic_geometry->automatic_edge_distance_trees.find(idx) ==
                   cached_automatic_geometry->automatic_edge_distance_trees.end();
      });
  if (use_automatic_edges)
  {
    MFEM_VERIFY(boundaries,
                "Automatic metal edge extraction requires complete boundary data!");
    metal_edges = ExtractMetalEdgeGeometry(mesh, *boundaries);
  }
  for (const auto &[idx, data] : postpro.dielectric)
  {
    MFEM_VERIFY(
        !data.flux_recovery || !cracked_attributes ||
            std::none_of(
                data.attributes.begin(), data.attributes.end(),
                [cracked_attributes](int attr)
                { return cracked_attributes->find(attr) != cracked_attributes->end(); }),
        "Interface dielectric \"FluxRecovery\" is not supported on cracked internal "
        "boundaries: the volume L2 projection does not provide a controlled normal trace "
        "on a zero-thickness PEC surface!");
    auto it = eps_surfs.try_emplace(idx, data, bdr_attr_marker).first;
    if (data.edge_distances.empty())
    {
      continue;
    }

    if (data.automatic_edges)
    {
      if (!data.localize_edge_energy && cached_automatic_geometry)
      {
        const auto cached_tree =
            cached_automatic_geometry->automatic_edge_distance_trees.find(idx);
        if (cached_tree != cached_automatic_geometry->automatic_edge_distance_trees.end())
        {
          it->second.edge_distance_tree = cached_tree->second;
          continue;
        }
      }

      auto segment_indices =
          GetInterfaceMetalEdgeSegmentIndices(metal_edges, idx, data.type);
      ExcludeMetalEdgeSegmentIndices(mesh, metal_edges, data.edge_exclude_attributes,
                                     segment_indices);
      const AutomaticEdgeDistanceTreeKey tree_key{segment_indices, data.edge_frame_normal};
      auto tree_it = automatic_edge_distance_trees.find(tree_key);
      if (tree_it == automatic_edge_distance_trees.end())
      {
        std::vector<std::array<double, 3>> process_normals;
        if (data.localize_edge_energy)
        {
          process_normals = BuildMetalEdgeProcessNormals(
              mesh, metal_edges, segment_indices, [this](int attribute)
              { return this->mat_op.GetLightSpeedMax(attribute); }, data.edge_frame_normal);
        }
        tree_it =
            automatic_edge_distance_trees
                .try_emplace(tree_key, BuildEdgeDistanceTree(metal_edges, segment_indices,
                                                             std::move(process_normals)))
                .first;
      }
      it->second.edge_distance_tree = tree_it->second;
      if (!data.localize_edge_energy && new_automatic_geometry)
      {
        new_automatic_geometry->automatic_edge_distance_trees.emplace(idx, tree_it->second);
      }
    }
    else
    {
      const EdgeDistanceTreeKey tree_key{data.edge_attributes, data.edge_exclude_attributes,
                                         data.edge_frame_normal};
      auto tree_it = edge_distance_trees.find(tree_key);
      if (tree_it == edge_distance_trees.end())
      {
        tree_it =
            edge_distance_trees
                .try_emplace(tree_key, BuildEdgeDistanceTree(mesh, data.edge_attributes,
                                                             data.edge_exclude_attributes,
                                                             data.edge_frame_normal))
                .first;
      }
      it->second.edge_distance_tree = tree_it->second;
    }
  }
  if (automatic_geometry && !*automatic_geometry && new_automatic_geometry &&
      !new_automatic_geometry->automatic_edge_distance_trees.empty())
  {
    *automatic_geometry = std::move(new_automatic_geometry);
  }

  // FarField postprocessing.
  MFEM_VERIFY(postpro.farfield.empty() || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE,
              "Far-field extraction is only available for driven and eigenmode problems!");

  // Check that we don't have anisotropic materials.
  if (!postpro.farfield.empty())
  {
    const auto &mesh = *nd_fespace.GetParMesh();
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker =
        mesh::AttrToMarker(bdr_attr_max, postpro.farfield.attributes);

    std::set<int> domain_attrs;

    for (int i = 0; i < mesh.GetNBE(); i++)
    {
      if (bdr_attr_marker[mesh.GetBdrAttribute(i) - 1])
      {
        int elem_id, _face_id;
        mesh.GetBdrElementAdjacentElement(i, elem_id, _face_id);
        if (elem_id >= 0)
        {
          domain_attrs.insert(mesh.GetAttribute(elem_id));
        }
      }
    }

    for (int attr : domain_attrs)
    {
      MFEM_VERIFY(mat_op.IsIsotropic(attr),
                  "FarField requires isotropic materials, but attribute " +
                      std::to_string(attr) + " is not.");
    }
  }

  farfield = FarFieldData(postpro.farfield, *nd_fespace.GetParMesh(), bdr_attr_marker);
}

SurfacePostOperator::SurfacePostOperator(
    const IoData &iodata, const MaterialOperator &mat_op,
    mfem::ParFiniteElementSpace &h1_fespace, mfem::ParFiniteElementSpace &nd_fespace,
    std::shared_ptr<const SurfacePostGeometry> *automatic_geometry)
  : SurfacePostOperator(iodata.boundaries.postpro, iodata.problem.type, mat_op, h1_fespace,
                        nd_fespace, &iodata.boundaries.cracked_attributes,
                        &iodata.boundaries, automatic_geometry)
{
}

std::complex<double> SurfacePostOperator::GetSurfaceFlux(int idx, const GridFunction *E,
                                                         const GridFunction *B) const
{
  // For complex-valued fields, output the separate real and imaginary parts for the time-
  // harmonic quantity. For power flux (Poynting vector), output only the stationary real
  // part and not the part which has double the frequency.
  auto it = flux_surfs.find(idx);
  MFEM_VERIFY(it != flux_surfs.end(),
              "Unknown surface flux postprocessing index requested!");
  const bool has_imag = (E) ? E->HasImag() : B->HasImag();
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);
  auto f =
      it->second.GetCoefficient(E ? &E->Real() : nullptr, B ? &B->Real() : nullptr, mat_op);
  std::complex<double> dot(GetLocalSurfaceIntegral(*f, attr_marker), 0.0);
  if (has_imag)
  {
    f = it->second.GetCoefficient(E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr,
                                  mat_op);
    double doti = GetLocalSurfaceIntegral(*f, attr_marker);
    if (it->second.type == SurfaceFlux::POWER)
    {
      dot += doti;
    }
    else
    {
      dot.imag(doti);
    }
  }
  Mpi::GlobalSum(1, &dot, (E) ? E->GetComm() : B->GetComm());
  return dot;
}

double SurfacePostOperator::GetInterfaceLossTangent(int idx) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  return it->second.tandelta;
}

double SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx, const GridFunction &E,
                                                            const GridFunction *D) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);
  auto f = it->second.GetCoefficient(E, D, mat_op);
  double dot = GetLocalSurfaceIntegral(*f, attr_marker);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::vector<SurfacePostOperator::InterfaceEdgeEnergy>
SurfacePostOperator::GetInterfaceEdgeElectricFieldEnergies(int idx, const GridFunction &E,
                                                           const GridFunction *D) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &data = it->second;
  if (data.edge_distances.empty())
  {
    return {};
  }

  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
  auto coefficient = data.GetCoefficient(E, D, mat_op);

  std::vector<InterfaceEdgeEnergy> energies;
  energies.reserve(data.edge_distances.size());
  std::vector<double> local_energy(2 * data.edge_distances.size());
  for (std::size_t i = 0; i < data.edge_distances.size(); i++)
  {
    const double distance = data.edge_distances[i];
    EdgeDistanceCoefficient outside(*coefficient, *data.edge_distance_tree, distance,
                                    mfem::infinity(), data.edge_distance_smoothing);
    EdgeDistanceCoefficient annulus(*coefficient, *data.edge_distance_tree, distance,
                                    2.0 * distance, data.edge_distance_smoothing);
    local_energy[2 * i] = GetLocalSurfaceIntegral(outside, attr_marker);
    local_energy[2 * i + 1] = GetLocalSurfaceIntegral(annulus, attr_marker);
  }
  Mpi::GlobalSum(static_cast<int>(local_energy.size()), local_energy.data(), E.GetComm());
  for (std::size_t i = 0; i < data.edge_distances.size(); i++)
  {
    energies.push_back(
        {data.edge_distances[i], local_energy[2 * i], local_energy[2 * i + 1]});
  }
  return energies;
}

SurfacePostOperator::InterfaceEdgeEnergy
SurfacePostOperator::GetInterfaceOuterElectricFieldEnergy(int idx, const GridFunction &E,
                                                          const GridFunction *D) const
{
  auto energies = GetInterfaceOuterElectricFieldEnergies({idx}, E, D);
  MFEM_ASSERT(energies.size() == 1 && energies.find(idx) != energies.end(),
              "Missing response-corrected interface energy!");
  return energies.at(idx);
}

std::map<int, SurfacePostOperator::InterfaceEdgeEnergy>
SurfacePostOperator::GetInterfaceOuterElectricFieldEnergies(const std::set<int> &indices,
                                                            const GridFunction &E,
                                                            const GridFunction *D) const
{
  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  std::map<int, InterfaceEdgeEnergy> energies;
  std::vector<double> local_energies;
  local_energies.reserve(indices.size());
  for (const int idx : indices)
  {
    auto it = eps_surfs.find(idx);
    MFEM_VERIFY(it != eps_surfs.end(),
                "Unknown interface dielectric postprocessing index requested!");
    const auto &data = it->second;
    MFEM_VERIFY(!data.edge_distances.empty(),
                "Response-corrected interface requires a matching edge distance!");

    const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
    auto coefficient = data.GetCoefficient(E, D, mat_op);
    const double distance = data.edge_distances.back();
    EdgeDistanceCoefficient outside(*coefficient, *data.edge_distance_tree, distance,
                                    mfem::infinity(), data.edge_distance_smoothing);
    local_energies.push_back(GetLocalSurfaceIntegral(outside, attr_marker));
    energies.emplace(idx, InterfaceEdgeEnergy{distance, 0.0, 0.0});
  }
  Mpi::GlobalSum(static_cast<int>(local_energies.size()), local_energies.data(),
                 E.GetComm());
  std::size_t i = 0;
  for (auto &[idx, energy] : energies)
  {
    (void)idx;
    energy.energy_outside = local_energies[i++];
  }
  return energies;
}

std::size_t SurfacePostOperator::GetNInterfaceEdgeEntries() const
{
  std::size_t size = 0;
  for (const auto &[idx, data] : eps_surfs)
  {
    size += data.edge_distances.size();
  }
  return size;
}

std::vector<SurfacePostOperator::InterfaceLocalEdgeEnergy>
SurfacePostOperator::GetInterfaceLocalEdgeElectricFieldEnergies(int idx,
                                                                const GridFunction &E,
                                                                const GridFunction *D,
                                                                bool include_volume) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &data = it->second;
  if (!data.localize_edge_energy)
  {
    return {};
  }

  const auto &mesh = *h1_fespace.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);
  const bool polarized = data.type != InterfaceDielectric::DEFAULT;
  auto coefficient = polarized ? nullptr : data.GetCoefficient(E, D, mat_op);
  auto normal_coefficient =
      polarized ? data.GetCoefficient(E, D, mat_op, InterfaceDielectricComponent::NORMAL)
                : nullptr;
  auto tangential_coefficient =
      polarized
          ? data.GetCoefficient(E, D, mat_op, InterfaceDielectricComponent::TANGENTIAL)
          : nullptr;
  const std::size_t edge_count = data.edge_distance_tree->Size();
  const std::size_t radius_count = data.edge_distances.size();
  std::vector<double> local_total_energy(edge_count, 0.0);
  std::vector<double> local_total_polarized_energy(2 * edge_count, 0.0);
  std::vector<double> local_energy(2 * radius_count * edge_count, 0.0);
  std::vector<double> local_vertex_energy(radius_count * edge_count, 0.0);
  std::vector<double> local_polarized_energy(4 * radius_count * edge_count, 0.0);
  const double max_distance =
      2.0 * data.edge_distances.back() * (1.0 + data.edge_distance_smoothing);

  mfem::Vector point(mesh.SpaceDimension());
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
    const auto *fe = h1_fespace.GetBE(be);
    const auto &ir =
        mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      T->Transform(ip, point);
      const auto nearest = data.edge_distance_tree->Nearest(point);
      const double distance = std::sqrt(nearest.distance_squared);
      const double vertex_distance =
          data.edge_distance_tree->DistanceAlongEdgeToNonregularVertex(point,
                                                                       nearest.segment);
      const double integration_weight = ip.weight * T->Weight();
      const double normal_energy =
          polarized ? integration_weight * normal_coefficient->Eval(*T, ip) : 0.0;
      const double tangential_energy =
          polarized ? integration_weight * tangential_coefficient->Eval(*T, ip) : 0.0;
      const double energy = polarized ? normal_energy + tangential_energy
                                      : integration_weight * coefficient->Eval(*T, ip);
      local_total_energy[nearest.segment] += energy;
      local_total_polarized_energy[2 * nearest.segment] += normal_energy;
      local_total_polarized_energy[2 * nearest.segment + 1] += tangential_energy;
      if (distance >= max_distance)
      {
        continue;
      }

      for (std::size_t radius = 0; radius < radius_count; radius++)
      {
        const std::size_t entry = radius * edge_count + nearest.segment;
        const double matching_distance = data.edge_distances[radius];
        const double outside_weight = EdgeDistanceOutsideWeight(
            distance, matching_distance, data.edge_distance_smoothing);
        const double annulus_weight =
            EdgeDistanceWindowWeight(distance, matching_distance, 2.0 * matching_distance,
                                     data.edge_distance_smoothing);
        const double vertex_weight =
            1.0 - EdgeDistanceOutsideWeight(vertex_distance, matching_distance,
                                            data.edge_distance_smoothing);
        local_energy[2 * entry] += (1.0 - outside_weight) * energy;
        local_energy[2 * entry + 1] += annulus_weight * energy;
        local_vertex_energy[entry] += (1.0 - outside_weight) * vertex_weight * energy;
        local_polarized_energy[4 * entry] += (1.0 - outside_weight) * normal_energy;
        local_polarized_energy[4 * entry + 1] += annulus_weight * normal_energy;
        local_polarized_energy[4 * entry + 2] += (1.0 - outside_weight) * tangential_energy;
        local_polarized_energy[4 * entry + 3] += annulus_weight * tangential_energy;
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(local_total_energy.size()), local_total_energy.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_total_polarized_energy.size()),
                 local_total_polarized_energy.data(), E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_energy.size()), local_energy.data(), E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_vertex_energy.size()), local_vertex_energy.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(local_polarized_energy.size()),
                 local_polarized_energy.data(), E.GetComm());
  const LocalVolumeEdgeEnergyCache *local_volume_energy = nullptr;
  if (include_volume)
  {
    local_volume_energy = &GetLocalVolumeEdgeElectricFieldEnergies(data, E);
  }

  std::vector<InterfaceLocalEdgeEnergy> energies;
  energies.reserve(local_energy.size());
  for (std::size_t radius = 0; radius < radius_count; radius++)
  {
    for (std::size_t edge = 0; edge < edge_count; edge++)
    {
      const auto &segment = data.edge_distance_tree->GetSegment(edge);
      double length_squared = 0.0;
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        const double delta = segment.p1[d] - segment.p0[d];
        length_squared += delta * delta;
      }
      const std::size_t entry = radius * edge_count + edge;
      const auto frame =
          BuildEdgeFrame(segment,
                         data.edge_distance_tree->HasProcessNormals()
                             ? data.edge_distance_tree->GetProcessNormal(edge)
                             : data.edge_frame_normal,
                         mesh.SpaceDimension());
      int component = -1;
      int chain = -1;
      std::array<int, 2> vertex_types{-1, -1};
      if (data.edge_distance_tree->HasMetadata())
      {
        const auto &metadata = data.edge_distance_tree->GetMetadata(edge);
        component = metadata.component;
        chain = metadata.chain;
        vertex_types = metadata.vertex_types;
      }
      InterfaceLocalEdgeEnergy energy{
          static_cast<int>(edge + 1),
          data.edge_distance_tree->HasMetadata(),
          component,
          chain,
          vertex_types,
          frame.normal,
          segment.p0,
          segment.p1,
          std::sqrt(length_squared),
          data.edge_distances[radius],
          local_total_energy[edge],
          local_energy[2 * entry],
          local_energy[2 * entry + 1],
          local_vertex_energy[entry],
          {local_total_polarized_energy[2 * edge],
           local_total_polarized_energy[2 * edge + 1]},
          {local_polarized_energy[4 * entry], local_polarized_energy[4 * entry + 2]},
          {local_polarized_energy[4 * entry + 1], local_polarized_energy[4 * entry + 3]},
          local_volume_energy ? local_volume_energy->energies[entry] : 0.0,
          local_volume_energy ? local_volume_energy->vertex_energies[entry] : 0.0,
          {}};
      if (local_volume_energy)
      {
        for (std::size_t component = 0; component < 6; component++)
        {
          energy.energy_volume_annulus_polarized[component] =
              local_volume_energy->polarized_energies[component][entry];
        }
      }
      energies.push_back(std::move(energy));
    }
  }
  return energies;
}

template <InterfaceDielectric Type>
std::vector<SurfacePostOperator::InterfaceResponseMatrix>
SurfacePostOperator::GetInterfaceElectricFieldEnergyMatricesImpl(
    const InterfaceDielectricData &data, const std::vector<const GridFunction *> &E,
    const std::vector<const GridFunction *> &D) const
{
  MFEM_VERIFY(!E.empty() && (D.empty() || D.size() == E.size()),
              "Invalid batched interface response fields!");
  MFEM_VERIFY(!data.flux_recovery || !D.empty(),
              "Batched interface response requires recovered electric flux!");
  const auto &mesh = *h1_fespace.GetParMesh();
  const int sdim = mesh.SpaceDimension();
  const int basis_size = static_cast<int>(E.size());
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto attr_marker = mesh::AttrToMarker(bdr_attr_max, data.attr_list);

  std::vector<std::unique_ptr<InterfaceDielectricCoefficient<Type>>> evaluators;
  evaluators.reserve(E.size());
  for (std::size_t i = 0; i < E.size(); i++)
  {
    evaluators.push_back(std::make_unique<InterfaceDielectricCoefficient<Type>>(
        *E[i], mat_op, data.t, data.epsilon, data.flux_recovery ? D[i] : nullptr));
  }

  std::vector<double> weights, distances, amplitudes;
  mfem::Vector point(sdim), field;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attr <= 0 || attr > attr_marker.Size() || !attr_marker[attr - 1])
    {
      continue;
    }
    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(be);
    const auto *fe = h1_fespace.GetBE(be);
    const auto &ir =
        mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T->SetIntPoint(&ip);
      T->Transform(ip, point);
      const auto nearest = data.edge_distance_tree->Nearest(point);
      distances.push_back(std::sqrt(nearest.distance_squared));
      weights.push_back(ip.weight * T->Weight());
      for (const auto &evaluator : evaluators)
      {
        evaluator->EvalEnergyField(*T, ip, field);
        MFEM_ASSERT(field.Size() == sdim + 1,
                    "Invalid batched interface energy field size!");
        amplitudes.insert(amplitudes.end(), field.begin(), field.end());
      }
    }
  }

  const int sample_count = static_cast<int>(weights.size());
  auto Assemble = [&](const std::vector<double> *inside_weights, mfem::DenseMatrix &normal,
                      mfem::DenseMatrix &tangential)
  {
    normal.SetSize(basis_size);
    normal = 0.0;
    tangential.SetSize(basis_size);
    tangential = 0.0;
    if (sample_count == 0)
    {
      return;
    }
    mfem::DenseMatrix normal_samples(basis_size, sample_count);
    mfem::DenseMatrix tangential_samples(basis_size, sdim * sample_count);
    for (int sample = 0; sample < sample_count; sample++)
    {
      const double window = inside_weights ? (*inside_weights)[sample] : 1.0;
      const double scale = std::sqrt(0.5 * weights[sample] * window);
      for (int basis = 0; basis < basis_size; basis++)
      {
        const std::size_t offset =
            (static_cast<std::size_t>(sample) * basis_size + basis) * (sdim + 1);
        normal_samples(basis, sample) = scale * amplitudes[offset];
        for (int d = 0; d < sdim; d++)
        {
          tangential_samples(basis, sdim * sample + d) = scale * amplitudes[offset + d + 1];
        }
      }
    }
    mfem::AddMult_a_AAt(1.0, normal_samples, normal);
    mfem::AddMult_a_AAt(1.0, tangential_samples, tangential);
  };
  auto Reduce = [&](mfem::DenseMatrix &matrix)
  {
    Mpi::GlobalSum(matrix.Height() * matrix.Width(), matrix.GetData(),
                   E.front()->GetComm());
  };

  mfem::DenseMatrix total_normal, total_tangential;
  Assemble(nullptr, total_normal, total_tangential);
  Reduce(total_normal);
  Reduce(total_tangential);

  std::vector<InterfaceResponseMatrix> result;
  result.reserve(data.edge_distances.size());
  std::vector<double> inside_weights(sample_count);
  for (const double radius : data.edge_distances)
  {
    for (int sample = 0; sample < sample_count; sample++)
    {
      inside_weights[sample] =
          1.0 - EdgeDistanceOutsideWeight(distances[sample], radius,
                                          data.edge_distance_smoothing);
    }
    InterfaceResponseMatrix entry;
    entry.distance = radius;
    entry.energy_total_normal = total_normal;
    entry.energy_total_tangential = total_tangential;
    Assemble(&inside_weights, entry.energy_inside_normal, entry.energy_inside_tangential);
    Reduce(entry.energy_inside_normal);
    Reduce(entry.energy_inside_tangential);
    entry.energy_total = entry.energy_total_normal;
    entry.energy_total += entry.energy_total_tangential;
    entry.energy_inside = entry.energy_inside_normal;
    entry.energy_inside += entry.energy_inside_tangential;
    result.push_back(std::move(entry));
  }
  return result;
}

std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>>
SurfacePostOperator::GetInterfaceElectricFieldEnergyMatrices(
    const std::vector<const GridFunction *> &E,
    const std::vector<const GridFunction *> &D) const
{
  std::map<int, std::vector<InterfaceResponseMatrix>> result;
  for (const auto &[idx, data] : eps_surfs)
  {
    if (!data.localize_edge_energy)
    {
      continue;
    }
    switch (data.type)
    {
      case InterfaceDielectric::DEFAULT:
        result.emplace(
            idx, GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::DEFAULT>(
                     data, E, D));
        break;
      case InterfaceDielectric::MA:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::MA>(
                           data, E, D));
        break;
      case InterfaceDielectric::MS:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::MS>(
                           data, E, D));
        break;
      case InterfaceDielectric::SA:
        result.emplace(idx,
                       GetInterfaceElectricFieldEnergyMatricesImpl<InterfaceDielectric::SA>(
                           data, E, D));
        break;
    }
  }
  return result;
}

const SurfacePostOperator::LocalVolumeEdgeEnergyCache &
SurfacePostOperator::GetLocalVolumeEdgeElectricFieldEnergies(
    const InterfaceDielectricData &data, const GridFunction &E) const
{
  for (const auto &entry : local_volume_edge_energy_cache)
  {
    if (entry.edge_distance_tree == data.edge_distance_tree.get() &&
        entry.edge_distances == data.edge_distances &&
        entry.edge_distance_smoothing == data.edge_distance_smoothing &&
        entry.edge_frame_normal == data.edge_frame_normal)
    {
      return entry;
    }
  }

  const auto &mesh = *E.ParFESpace()->GetParMesh();
  const std::size_t edge_count = data.edge_distance_tree->Size();
  const std::size_t radius_count = data.edge_distances.size();
  auto &entry = local_volume_edge_energy_cache.emplace_back(
      LocalVolumeEdgeEnergyCache{data.edge_distance_tree.get(),
                                 data.edge_distances,
                                 data.edge_distance_smoothing,
                                 data.edge_frame_normal,
                                 std::vector<double>(radius_count * edge_count, 0.0),
                                 std::vector<double>(radius_count * edge_count, 0.0),
                                 {}});
  for (auto &component : entry.polarized_energies)
  {
    component.assign(radius_count * edge_count, 0.0);
  }
  const double max_distance =
      2.0 * data.edge_distances.back() * (1.0 + data.edge_distance_smoothing);

  std::vector<EdgeFrame> frames;
  frames.reserve(edge_count);
  for (std::size_t edge = 0; edge < edge_count; edge++)
  {
    const auto &normal = data.edge_distance_tree->HasProcessNormals()
                             ? data.edge_distance_tree->GetProcessNormal(edge)
                             : data.edge_frame_normal;
    frames.push_back(BuildEdgeFrame(data.edge_distance_tree->GetSegment(edge), normal,
                                    mesh.SpaceDimension()));
  }

  mfem::Vector point(mesh.SpaceDimension());
  mfem::Vector field(mesh.SpaceDimension()), displacement(mesh.SpaceDimension());
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mfem::IsoparametricTransformation T;
    mesh.GetElementTransformation(e, &T);
    const auto &ir =
        mfem::IntRules.Get(T.GetGeometryType(), fem::DefaultIntegrationOrder::Get(T));
    for (int q = 0; q < ir.GetNPoints(); q++)
    {
      const auto &ip = ir.IntPoint(q);
      T.SetIntPoint(&ip);
      T.Transform(ip, point);
      const auto nearest = data.edge_distance_tree->Nearest(point);
      const double distance = std::sqrt(nearest.distance_squared);
      if (distance >= max_distance)
      {
        continue;
      }

      const auto &frame = frames[nearest.segment];
      const auto &segment = data.edge_distance_tree->GetSegment(nearest.segment);
      const double vertex_distance =
          data.edge_distance_tree->DistanceAlongEdgeToNonregularVertex(point,
                                                                       nearest.segment);
      std::array<double, 6> component_energy{};
      auto AddFieldEnergy = [&](const mfem::ParGridFunction &grid_function)
      {
        grid_function.GetVectorValue(T, ip, field);
        mat_op.GetPermittivityReal(T.Attribute).Mult(field, displacement);
        const auto polarized =
            GetPolarizedEdgeEnergyDensity(point, segment, frame, field, displacement);
        for (std::size_t component = 0; component < 6; component++)
        {
          component_energy[component] += polarized[component];
        }
      };
      AddFieldEnergy(E.Real());
      if (E.HasImag())
      {
        AddFieldEnergy(E.Imag());
      }
      const double integration_weight = ip.weight * T.Weight();
      for (std::size_t radius = 0; radius < radius_count; radius++)
      {
        const double matching_distance = data.edge_distances[radius];
        const double annulus_weight =
            EdgeDistanceWindowWeight(distance, matching_distance, 2.0 * matching_distance,
                                     data.edge_distance_smoothing);
        const double vertex_weight =
            1.0 - EdgeDistanceOutsideWeight(vertex_distance, matching_distance,
                                            data.edge_distance_smoothing);
        const std::size_t result = radius * edge_count + nearest.segment;
        for (std::size_t component = 0; component < 6; component++)
        {
          const double energy =
              annulus_weight * integration_weight * component_energy[component];
          entry.polarized_energies[component][result] += energy;
          entry.energies[result] += energy;
          entry.vertex_energies[result] += vertex_weight * energy;
        }
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(entry.energies.size()), entry.energies.data(),
                 E.GetComm());
  Mpi::GlobalSum(static_cast<int>(entry.vertex_energies.size()),
                 entry.vertex_energies.data(), E.GetComm());
  for (auto &component : entry.polarized_energies)
  {
    Mpi::GlobalSum(static_cast<int>(component.size()), component.data(), E.GetComm());
  }
  return entry;
}

std::size_t SurfacePostOperator::GetNInterfaceLocalEdgeEntries() const
{
  std::size_t size = 0;
  for (const auto &[idx, data] : eps_surfs)
  {
    if (data.localize_edge_energy && data.save_local_edge_energy)
    {
      size += data.edge_distances.size() * data.edge_distance_tree->Size();
    }
  }
  return size;
}

bool SurfacePostOperator::NeedsFluxRecovery() const
{
  return std::any_of(eps_surfs.begin(), eps_surfs.end(),
                     [](const auto &entry) { return entry.second.flux_recovery; });
}

void SurfacePostOperator::ResetInterfaceLocalEdgeEnergyCache() const
{
  local_volume_edge_energy_cache.clear();
}

double
SurfacePostOperator::GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                             const mfem::Array<int> &attr_marker) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  mfem::LinearForm s(&h1_fespace);
  s.AddBoundaryIntegrator(new BoundaryLFIntegrator(f),
                          const_cast<mfem::Array<int> &>(attr_marker));
  s.UseFastAssembly(false);
  s.UseDevice(false);
  s.Assemble();
  s.UseDevice(true);
  return linalg::LocalSum(s);
}

std::vector<std::array<std::complex<double>, 3>> SurfacePostOperator::GetFarFieldrE(
    const std::vector<std::pair<double, double>> &theta_phi_pairs, const GridFunction &E,
    const GridFunction &B, double omega_re, double omega_im) const
{
  if (theta_phi_pairs.empty())
    return {};
  MFEM_VERIFY(nd_fespace.GetParMesh()->SpaceDimension() == 3,
              "Far-field computation is only available for 3D simulations!");
  BlockTimer bt0(Timer::POSTPRO_FARFIELD);

  // Compute target unit vectors from the given theta and phis.
  std::vector<std::array<double, 3>> r_naughts;
  r_naughts.reserve(theta_phi_pairs.size());

  r_naughts.reserve(theta_phi_pairs.size());
  for (const auto &[theta, phi] : theta_phi_pairs)
  {
    r_naughts.emplace_back(std::array<double, 3>{
        std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)});
  }

  const auto &mesh = *nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, farfield.attr_list);

  // Integrate. Each MPI process computes its contribution and we will reduce
  // everything at the end. We make them std::vector<std::array<double, 3>>
  // because we want a very simple memory layout so that we can reduce
  // everything with two MPI calls.
  std::vector<std::array<double, 3>> integrals_r(theta_phi_pairs.size());
  std::vector<std::array<double, 3>> integrals_i(theta_phi_pairs.size());

  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    if (!attr_marker[mesh.GetBdrAttribute(i) - 1])
      continue;

    auto *T = const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(i);
    const auto *fe = nd_fespace.GetBE(i);
    const auto *ir =
        &mfem::IntRules.Get(fe->GetGeomType(), fem::DefaultIntegrationOrder::Get(*T));

    AddStrattonChuIntegrandAtElement(E, B, mat_op, omega_re, omega_im, r_naughts, *T, *ir,
                                     integrals_r, integrals_i);
  }

  double *data_r_ptr = integrals_r.data()->data();
  double *data_i_ptr = integrals_i.data()->data();
  std::size_t total_elements = integrals_r.size() * 3;
  Mpi::GlobalSum(total_elements, data_i_ptr, E.GetComm());
  Mpi::GlobalSum(total_elements, data_r_ptr, E.GetComm());

  // Finally, we apply cross product to reduced integrals and package the result
  // in a neatly accessible vector of arrays of complex numbers.
  std::vector<std::array<std::complex<double>, 3>> result(theta_phi_pairs.size());
  StaticVector<3> tmp_r, tmp_i;
  for (std::size_t k = 0; k < theta_phi_pairs.size(); k++)
  {
    linalg::Cross3(r_naughts[k], integrals_r[k], tmp_r);
    linalg::Cross3(r_naughts[k], integrals_i[k], tmp_i);
    for (std::size_t d = 0; d < 3; d++)
    {
      result[k][d] = std::complex<double>{tmp_r[d], tmp_i[d]};
    }
  }
  return result;
}

}  // namespace palace
