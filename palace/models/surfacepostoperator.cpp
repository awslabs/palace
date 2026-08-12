// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacepostoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/singulargeometry.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/strattonchu.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"

namespace palace
{

// Globally replicated, read-only geometry for a hard physical cutoff around the
// topological boundary of a selected interface surface. The primitives are points in 2D
// and finite straight segments in 3D. A balanced AABB tree keeps point queries logarithmic
// for the large sheet perimeters encountered in device-scale meshes.
class SurfaceEdgeCutoff
{
public:
  using Point = std::array<double, 3>;

private:
  struct Primitive
  {
    Point first, second, minimum, maximum, center;
  };
  struct Node
  {
    Point minimum, maximum;
    int begin = 0, end = 0, left = -1, right = -1;
  };

  int space_dimension;
  std::vector<Primitive> primitives;
  std::vector<int> order;
  std::vector<Node> nodes;

  static double BoxDistanceSquared(const Point &point, const Node &node,
                                   int space_dimension)
  {
    double distance_squared = 0.0;
    for (int d = 0; d < space_dimension; d++)
    {
      const double distance =
          point[d] < node.minimum[d]
              ? node.minimum[d] - point[d]
              : (point[d] > node.maximum[d] ? point[d] - node.maximum[d] : 0.0);
      distance_squared += distance * distance;
    }
    return distance_squared;
  }

  static double PrimitiveDistanceSquared(const Point &point, const Primitive &primitive,
                                         int space_dimension)
  {
    double length_squared = 0.0, projection = 0.0;
    for (int d = 0; d < space_dimension; d++)
    {
      const double tangent = primitive.second[d] - primitive.first[d];
      length_squared += tangent * tangent;
      projection += (point[d] - primitive.first[d]) * tangent;
    }
    const double coordinate =
        length_squared > 0.0 ? std::clamp(projection / length_squared, 0.0, 1.0) : 0.0;
    double distance_squared = 0.0;
    for (int d = 0; d < space_dimension; d++)
    {
      const double closest =
          primitive.first[d] + coordinate * (primitive.second[d] - primitive.first[d]);
      const double distance = point[d] - closest;
      distance_squared += distance * distance;
    }
    return distance_squared;
  }

  int BuildNode(int begin, int end)
  {
    Node node;
    node.begin = begin;
    node.end = end;
    node.minimum.fill(std::numeric_limits<double>::infinity());
    node.maximum.fill(-std::numeric_limits<double>::infinity());
    Point center_minimum, center_maximum;
    center_minimum.fill(std::numeric_limits<double>::infinity());
    center_maximum.fill(-std::numeric_limits<double>::infinity());
    for (int i = begin; i < end; i++)
    {
      const auto &primitive = primitives[order[i]];
      for (int d = 0; d < space_dimension; d++)
      {
        node.minimum[d] = std::min(node.minimum[d], primitive.minimum[d]);
        node.maximum[d] = std::max(node.maximum[d], primitive.maximum[d]);
        center_minimum[d] = std::min(center_minimum[d], primitive.center[d]);
        center_maximum[d] = std::max(center_maximum[d], primitive.center[d]);
      }
    }
    const int node_index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    constexpr int LeafSize = 8;
    if (end - begin <= LeafSize)
    {
      return node_index;
    }
    int axis = 0;
    for (int d = 1; d < space_dimension; d++)
    {
      if (center_maximum[d] - center_minimum[d] >
          center_maximum[axis] - center_minimum[axis])
      {
        axis = d;
      }
    }
    const int middle = begin + (end - begin) / 2;
    std::nth_element(
        order.begin() + begin, order.begin() + middle, order.begin() + end,
        [this, axis](int left, int right)
        { return primitives[left].center[axis] < primitives[right].center[axis]; });
    const int left = BuildNode(begin, middle);
    const int right = BuildNode(middle, end);
    nodes[node_index].left = left;
    nodes[node_index].right = right;
    return node_index;
  }

public:
  SurfaceEdgeCutoff(int space_dimension,
                    const std::vector<std::pair<Point, Point>> &segments)
    : space_dimension(space_dimension)
  {
    primitives.reserve(segments.size());
    for (const auto &[first, second] : segments)
    {
      Primitive primitive{first, second, first, first, first};
      for (int d = 0; d < space_dimension; d++)
      {
        primitive.minimum[d] = std::min(first[d], second[d]);
        primitive.maximum[d] = std::max(first[d], second[d]);
        primitive.center[d] = 0.5 * (first[d] + second[d]);
      }
      primitives.push_back(primitive);
    }
    order.resize(primitives.size());
    std::iota(order.begin(), order.end(), 0);
    nodes.reserve(2 * primitives.size());
    if (!primitives.empty())
    {
      BuildNode(0, static_cast<int>(primitives.size()));
    }
  }

  bool Empty() const { return primitives.empty(); }
  std::size_t Size() const { return primitives.size(); }

  bool Contains(const mfem::Vector &physical_point, double cutoff) const
  {
    if (primitives.empty() || !(cutoff > 0.0))
    {
      return false;
    }
    Point point{0.0, 0.0, 0.0};
    std::copy_n(physical_point.begin(), space_dimension, point.begin());
    const double cutoff_squared = cutoff * cutoff;
    std::array<int, 64> stack{};
    int stack_size = 1;
    stack[0] = 0;
    while (stack_size > 0)
    {
      const int node_index = stack[--stack_size];
      const auto &node = nodes[node_index];
      if (BoxDistanceSquared(point, node, space_dimension) >= cutoff_squared)
      {
        continue;
      }
      if (node.left < 0)
      {
        for (int i = node.begin; i < node.end; i++)
        {
          if (PrimitiveDistanceSquared(point, primitives[order[i]], space_dimension) <
              cutoff_squared)
          {
            return true;
          }
        }
      }
      else
      {
        MFEM_ASSERT(stack_size <= static_cast<int>(stack.size()) - 2,
                    "Standard surface EdgeCutoff BVH query stack overflow!");
        stack[stack_size++] = node.left;
        stack[stack_size++] = node.right;
      }
    }
    return false;
  }
};

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

using CutoffPoint = SurfaceEdgeCutoff::Point;

struct CutoffFaceKey
{
  int size;
  std::array<CutoffPoint, 4> points;

  bool operator<(const CutoffFaceKey &other) const
  {
    return std::tie(size, points) < std::tie(other.size, other.points);
  }
};

struct CutoffFace
{
  int size;
  std::array<CutoffPoint, 4> points;
  bool selected;
};

std::shared_ptr<const SurfaceEdgeCutoff>
BuildSurfaceEdgeCutoff(const mfem::ParMesh &mesh, const std::vector<int> &attributes)
{
  constexpr int MaximumVertices = 4;
  constexpr int RecordWidth = 2 + 3 * MaximumVertices;
  const int space_dimension = mesh.SpaceDimension();
  MFEM_VERIFY(space_dimension == 2 || space_dimension == 3,
              "Standard surface EdgeCutoff requires a two- or three-dimensional mesh!");
  if (space_dimension == 3 && mesh.Nonconforming())
  {
    MFEM_VERIFY(mesh.ncmesh,
                "Standard surface EdgeCutoff found an invalid nonconforming mesh!");
    const auto &edges = mesh.ncmesh->GetEdgeList();
    const auto &faces = mesh.ncmesh->GetFaceList();
    MFEM_VERIFY(edges.masters.Size() == 0 && edges.slaves.Size() == 0 &&
                    faces.masters.Size() == 0 && faces.slaves.Size() == 0,
                "Standard surface EdgeCutoff currently requires a conforming initial "
                "mesh; hanging boundary seams are not yet canonicalized!");
  }
  const std::set<int> selected_attributes(attributes.begin(), attributes.end());

  std::vector<double> local_data;
  local_data.reserve(RecordWidth * mesh.GetNBE());
  for (int boundary = 0; boundary < mesh.GetNBE(); boundary++)
  {
    const auto &element = *mesh.GetBdrElement(boundary);
    const bool selected = selected_attributes.count(mesh.GetBdrAttribute(boundary)) != 0;
    if (selected)
    {
      auto *transformation =
          const_cast<mfem::ParMesh &>(mesh).GetBdrElementTransformation(boundary);
      MFEM_VERIFY(transformation &&
                      fem::singular::IsAffineElementTransformation(*transformation),
                  "Standard surface EdgeCutoff currently requires affine selected "
                  "interface faces with straight boundary features!");
    }
    const int vertices = element.GetNVertices();
    const int expected_vertices = space_dimension == 2 ? 2 : vertices;
    MFEM_VERIFY(vertices == expected_vertices && vertices >= 2 &&
                    vertices <= MaximumVertices,
                "Standard surface EdgeCutoff requires segment boundaries in 2D and "
                "triangle or quadrilateral boundaries in 3D!");
    local_data.push_back(static_cast<double>(vertices));
    local_data.push_back(selected ? 1.0 : 0.0);
    for (int vertex = 0; vertex < MaximumVertices; vertex++)
    {
      for (int d = 0; d < 3; d++)
      {
        const double coordinate = vertex < vertices && d < space_dimension
                                      ? mesh.GetVertex(element.GetVertices()[vertex])[d]
                                      : 0.0;
        MFEM_VERIFY(std::isfinite(coordinate),
                    "Standard surface EdgeCutoff found a nonfinite mesh vertex!");
        local_data.push_back(coordinate == 0.0 ? 0.0 : coordinate);
      }
    }
  }
  MFEM_VERIFY(local_data.size() <=
                  static_cast<std::size_t>(std::numeric_limits<int>::max()),
              "Too many local boundary faces for standard surface EdgeCutoff!");
  const int local_count = static_cast<int>(local_data.size());
  std::vector<int> counts(Mpi::Size(mesh.GetComm()));
  Mpi::Allgather(1, &local_count, counts.data(), mesh.GetComm());
  std::vector<int> offsets(counts.size());
  int global_count = 0;
  for (std::size_t rank = 0; rank < counts.size(); rank++)
  {
    MFEM_VERIFY(counts[rank] >= 0 && counts[rank] % RecordWidth == 0 &&
                    counts[rank] <= std::numeric_limits<int>::max() - global_count,
                "Invalid distributed boundary-face count for standard surface "
                "EdgeCutoff!");
    offsets[rank] = global_count;
    global_count += counts[rank];
  }
  std::vector<double> global_data(global_count);
  Mpi::Allgatherv(local_count, local_data.data(), global_data.data(), counts.data(),
                  offsets.data(), mesh.GetComm());

  std::map<CutoffFaceKey, CutoffFace> faces;
  for (int offset = 0; offset < global_count; offset += RecordWidth)
  {
    const int vertices = static_cast<int>(global_data[offset]);
    MFEM_VERIFY(vertices >= 2 && vertices <= MaximumVertices &&
                    (space_dimension == 3 || vertices == 2),
                "A gathered standard surface EdgeCutoff face is invalid!");
    CutoffFace face{vertices, {}, global_data[offset + 1] != 0.0};
    int cursor = offset + 2;
    for (int vertex = 0; vertex < MaximumVertices; vertex++)
    {
      for (int d = 0; d < 3; d++)
      {
        face.points[vertex][d] = global_data[cursor++];
      }
    }
    CutoffFaceKey key{vertices, face.points};
    std::sort(key.points.begin(), key.points.begin() + vertices);
    const auto [it, inserted] = faces.emplace(key, face);
    MFEM_VERIFY(inserted || it->second.selected == face.selected,
                "A duplicated standard surface EdgeCutoff face has inconsistent "
                "attributes!");
  }

  std::vector<std::pair<CutoffPoint, CutoffPoint>> primitives;
  if (space_dimension == 2)
  {
    std::map<CutoffPoint, std::array<int, 2>> incidence;
    for (const auto &[key, face] : faces)
    {
      for (int vertex = 0; vertex < face.size; vertex++)
      {
        auto &counts = incidence[face.points[vertex]];
        counts[0]++;
        counts[1] += face.selected ? 1 : 0;
      }
    }
    for (const auto &[point, counts] : incidence)
    {
      if (counts[1] > 0 && (counts[1] < counts[0] || counts[0] == 1))
      {
        primitives.emplace_back(point, point);
      }
    }
  }
  else
  {
    using EdgeKey = std::array<CutoffPoint, 2>;
    std::map<EdgeKey, std::array<int, 2>> incidence;
    for (const auto &[key, face] : faces)
    {
      for (int edge = 0; edge < face.size; edge++)
      {
        EdgeKey key{face.points[edge], face.points[(edge + 1) % face.size]};
        std::sort(key.begin(), key.end());
        auto &counts = incidence[key];
        counts[0]++;
        counts[1] += face.selected ? 1 : 0;
      }
    }
    for (const auto &[edge, counts] : incidence)
    {
      if (counts[1] > 0 && (counts[1] < counts[0] || counts[0] == 1))
      {
        primitives.emplace_back(edge[0], edge[1]);
      }
    }
  }
  return std::make_shared<SurfaceEdgeCutoff>(space_dimension, primitives);
}

class EdgeCutoffCoefficient : public mfem::Coefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coefficient;
  std::shared_ptr<const SurfaceEdgeCutoff> geometry;
  double cutoff;
  int space_dimension;

public:
  EdgeCutoffCoefficient(std::unique_ptr<mfem::Coefficient> coefficient,
                        std::shared_ptr<const SurfaceEdgeCutoff> geometry, double cutoff,
                        int space_dimension)
    : coefficient(std::move(coefficient)), geometry(std::move(geometry)), cutoff(cutoff),
      space_dimension(space_dimension)
  {
  }

  void SetTime(double time) override
  {
    mfem::Coefficient::SetTime(time);
    coefficient->SetTime(time);
  }

  double Eval(mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &point) override
  {
    double physical_data[3];
    mfem::Vector physical_point(physical_data, space_dimension);
    transformation.Transform(point, physical_point);
    return geometry->Contains(physical_point, cutoff)
               ? 0.0
               : coefficient->Eval(transformation, point);
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
    const config::InterfaceDielectricData &data, const mfem::ParMesh &mesh,
    const mfem::Array<int> &bdr_attr_marker)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);

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
  edge_cutoff = data.edge_cutoff;
  cutoff_attributes = data.attributes;
  MFEM_VERIFY(std::isfinite(edge_cutoff) && edge_cutoff >= 0.0,
              "Standard surface postprocessing received an invalid EdgeCutoff!");
}

std::unique_ptr<mfem::Coefficient>
SurfacePostOperator::InterfaceDielectricData::GetCoefficient(
    const GridFunction &E, const MaterialOperator &mat_op) const
{
  std::unique_ptr<mfem::Coefficient> coefficient;
  switch (type)
  {
    case InterfaceDielectric::DEFAULT:
      coefficient = std::make_unique<RestrictedCoefficient<
          InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>>>(
          attr_list, E, mat_op, t, epsilon);
      break;
    case InterfaceDielectric::MA:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MA>>>(
          attr_list, E, mat_op, t, epsilon);
      break;
    case InterfaceDielectric::MS:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::MS>>>(
          attr_list, E, mat_op, t, epsilon);
      break;
    case InterfaceDielectric::SA:
      coefficient = std::make_unique<
          RestrictedCoefficient<InterfaceDielectricCoefficient<InterfaceDielectric::SA>>>(
          attr_list, E, mat_op, t, epsilon);
      break;
  }
  if (edge_cutoff > 0.0 && cutoff_geometry && !cutoff_geometry->Empty())
  {
    coefficient = std::make_unique<EdgeCutoffCoefficient>(
        std::move(coefficient), cutoff_geometry, edge_cutoff,
        E.ParFESpace()->GetParMesh()->SpaceDimension());
  }
  return coefficient;
}

SurfacePostOperator::FarFieldData::FarFieldData(const config::FarFieldPostData &data,
                                                const mfem::ParMesh &mesh,
                                                const mfem::Array<int> &bdr_attr_marker)
  : thetaphis(data.thetaphis)
{
  // Store boundary attributes for this postprocessing boundary.
  attr_list = SetUpBoundaryProperties(data, bdr_attr_marker);
}

SurfacePostOperator::SurfacePostOperator(const config::BoundaryPostData &postpro,
                                         ProblemType problem_type,
                                         const MaterialOperator &mat_op,
                                         mfem::ParFiniteElementSpace &h1_fespace,
                                         mfem::ParFiniteElementSpace &nd_fespace)
  : mat_op(mat_op), h1_fespace(h1_fespace), nd_fespace(nd_fespace)
{
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
  for (const auto &[idx, data] : postpro.dielectric)
  {
    eps_surfs.try_emplace(idx, data, *h1_fespace.GetParMesh(), bdr_attr_marker);
  }
  // Geometry is independent of dielectric type and cutoff distance, so measurements using
  // the same selected interface attributes share one globally gathered BVH.
  std::map<std::vector<int>, std::shared_ptr<const SurfaceEdgeCutoff>> cutoff_geometries;
  for (auto &[idx, data] : eps_surfs)
  {
    if (!(data.edge_cutoff > 0.0))
    {
      continue;
    }
    auto attributes = data.cutoff_attributes;
    std::sort(attributes.begin(), attributes.end());
    attributes.erase(std::unique(attributes.begin(), attributes.end()), attributes.end());
    const auto [it, inserted] = cutoff_geometries.emplace(attributes, nullptr);
    if (inserted)
    {
      it->second = BuildSurfaceEdgeCutoff(mesh, attributes);
      if (it->second->Empty())
      {
        Mpi::Warning(mesh.GetComm(),
                     "Standard surface EdgeCutoff found no boundary features for "
                     "interface attributes {}! The cutoff will have no effect.\n",
                     fmt::join(attributes, ", "));
      }
      else
      {
        Mpi::Print(mesh.GetComm(),
                   " Standard surface EdgeCutoff: {} globally gathered boundary "
                   "features for interface attributes {}\n",
                   it->second->Size(), fmt::join(attributes, ", "));
      }
    }
    data.cutoff_geometry = it->second;
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

SurfacePostOperator::SurfacePostOperator(const IoData &iodata,
                                         const MaterialOperator &mat_op,
                                         mfem::ParFiniteElementSpace &h1_fespace,
                                         mfem::ParFiniteElementSpace &nd_fespace)
  : SurfacePostOperator(iodata.boundaries.postpro, iodata.problem.type, mat_op, h1_fespace,
                        nd_fespace)
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

double SurfacePostOperator::GetInterfaceElectricFieldEnergy(int idx,
                                                            const GridFunction &E) const
{
  auto it = eps_surfs.find(idx);
  MFEM_VERIFY(it != eps_surfs.end(),
              "Unknown interface dielectric postprocessing index requested!");
  const auto &mesh = *h1_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, it->second.attr_list);
  auto f = it->second.GetCoefficient(E, mat_op);
  // The cutoff mask is discontinuous inside cut elements. A high fixed minimum rule avoids
  // the severe quadrature-node thresholding of the legacy default-order boundary rule.
  constexpr int EdgeCutoffQuadratureOrder = 24;
  double dot = GetLocalSurfaceIntegral(
      *f, attr_marker, it->second.edge_cutoff > 0.0 ? EdgeCutoffQuadratureOrder : 0);
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

double SurfacePostOperator::GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                                    const mfem::Array<int> &attr_marker,
                                                    int minimum_quadrature_order) const
{
  // Integrate the coefficient over the boundary attributes making up this surface index.
  mfem::LinearForm s(&h1_fespace);
  s.AddBoundaryIntegrator(new BoundaryLFIntegrator(f, minimum_quadrature_order),
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
