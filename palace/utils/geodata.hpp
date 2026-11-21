// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class Units;

namespace config
{

struct ConductivityData;

}  // namespace config

namespace mesh
{

//
// Functions for mesh related functionality.
//

// Load a serial mesh from disk and perform all serial-stage preparation: AMR compat
// checks, cleanup, simplex/hex conversion, element reordering, serial uniform refinement,
// region-based (box/sphere) refinement, boundary cracking, and finalization. Returns a
// null pointer on ranks that do not hold a copy of the serial mesh. Called by main.cpp
// before Preprocess hooks on the solver mutate the serial mesh (e.g. BoundaryMode
// submesh extraction). With `parallel_nc`, a nonconforming mesh which Palace saved is read
// in parallel instead (see ReadNCMeshCoarse): every process gets its coarse mesh, which
// Partition refines, and iodata keeps the rest of the file for Partition.
std::unique_ptr<mfem::Mesh> Load(IoData &iodata, MPI_Comm comm, bool parallel_nc = true);

// Interior conductivity boundary attributes without "Crack" (all of their elements are
// interior), which are uncracked sheets by default. Serial mesh, no communication.
std::vector<int> GetDefaultConductivitySheetAttributes(
    const mfem::Mesh &mesh, const std::vector<config::ConductivityData> &conductivity);

// Partition and distribute a serial mesh prepared by Load, producing a parallel mesh.
// `smesh` is non-null only on loading ranks (see Load's contract).
std::unique_ptr<mfem::ParMesh> Partition(IoData &iodata, std::unique_ptr<mfem::Mesh> smesh,
                                         MPI_Comm comm);

// Convenience wrapper: Load followed by Partition with no Preprocess hook.
std::unique_ptr<mfem::ParMesh> ReadMesh(IoData &iodata, MPI_Comm comm);

// Maximum axis-aligned bbox extent of the pre-partition serial mesh, reduced over
// `comm` so every rank sees the same value. Returns 0 when no rank holds a mesh.
double ComputeReferenceLength(const std::unique_ptr<mfem::Mesh> &mesh, MPI_Comm comm);

// Refine the provided mesh according to the data in the input file (parallel uniform
// refinement only; box / sphere region refinement happens in Load on the serial mesh).
// If levels of refinement are requested, the refined meshes are stored in order of
// increased refinement. Ownership of the initial coarse mesh is inherited by the fine
// meshes and it should not be deleted. The fine mesh hierarchy is owned by the user.
void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

// Print parallel mesh statistics (element/vertex/edge/face counts, element size h, and
// shape regularity kappa). When `full` is true, also prints the (refinement-invariant)
// mesh curvature order and bounding box; pass false to print only the statistics that
// change under refinement (used after each adaptive mesh refinement iteration). `prefix`
// is prepended to the stats block (e.g. "Coarse "/"Refined ").
void PrintMeshInfo(mfem::ParMesh &mesh, const IoData &iodata, bool full = true,
                   const std::string &prefix = "");

// Dimensionalize a mesh for use in exporting a mesh. Scales vertices and nodes by L.
void DimensionalizeMesh(mfem::Mesh &mesh, double L);

// Nondimensionalize a mesh for use in the solver. Scales vertices and nodes by 1/L.
void NondimensionalizeMesh(mfem::Mesh &mesh, double L);
void Nondimensionalize(const Units &units, mfem::Mesh &mesh);

// Struct containing flags for the (global) mesh element types.
struct ElementTypeInfo
{
  bool has_simplices;
  bool has_hexahedra;
  bool has_prisms;
  bool has_pyramids;
  std::vector<mfem::Geometry::Type> GetGeomTypes(int dim = 3) const;
};

// Simplified helper for describing the element types in a (Par)Mesh.
ElementTypeInfo CheckElements(const mfem::Mesh &mesh);

// Check if a tetrahedral (Par)Mesh is ready for local refinement.
bool CheckRefinementFlags(const mfem::Mesh &mesh);

// Helper function to convert a set of attribute numbers to a marker array. The marker array
// will be of size max_attr and it will contain only zeroes and ones. Ones indicate which
// attribute numbers are present in the list array. In the special case when list has a
// single entry equal to -1 the marker array will contain all ones.
void AttrToMarker(int max_attr, const int *attr_list, int attr_list_size,
                  mfem::Array<int> &marker, bool skip_invalid = false);

template <typename T>
inline void AttrToMarker(int max_attr, const T &attr_list, mfem::Array<int> &marker,
                         bool skip_invalid = false)
{
  const auto size = std::distance(attr_list.begin(), attr_list.end());
  AttrToMarker(max_attr, (size > 0) ? &attr_list[0] : nullptr, size, marker, skip_invalid);
}

template <typename T>
inline mfem::Array<int> AttrToMarker(int max_attr, const T &attr_list,
                                     bool skip_invalid = false)
{
  mfem::Array<int> marker;
  AttrToMarker(max_attr, attr_list, marker, skip_invalid);
  return marker;
}

// Maximum boundary attribute over all local boundary elements. This can differ from
// mesh.bdr_attributes.Max() for internally-added boundary elements.
int GetMaxBdrAttribute(const mfem::ParMesh &mesh);

// Helper function to convert boundary attribute numbers to a marker array sized by the
// maximum actual boundary element attribute.
template <typename T>
inline mfem::Array<int> BdrAttrToMarker(const mfem::ParMesh &mesh, const T &attr_list,
                                        bool skip_invalid = false)
{
  return AttrToMarker(GetMaxBdrAttribute(mesh), attr_list, skip_invalid);
}

// Helper function to construct the axis-aligned bounding box for all elements with the
// given attribute.
void GetAxisAlignedBoundingBox(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                               bool bdr, mfem::Vector &min, mfem::Vector &max);

inline void GetAxisAlignedBoundingBox(const mfem::ParMesh &mesh, int attr, bool bdr,
                                      mfem::Vector &min, mfem::Vector &max)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetAxisAlignedBoundingBox(mesh, marker, bdr, min, max);
}

inline void GetAxisAlignedBoundingBox(const mfem::ParMesh &mesh, mfem::Vector &min,
                                      mfem::Vector &max)
{
  mfem::Array<int> marker(mesh.attributes.Max());
  marker = 1;
  GetAxisAlignedBoundingBox(mesh, marker, false, min, max);
}

// Struct describing a bounding box in terms of its center and the vectors from the center
// to its faces, as returned by GetBoundingBall for a bounding ball. Supports both 2D and
// 3D: in 2D, center has 2 entries and axes is 2x2; in 3D, center has 3 entries and axes is
// 3x3.
struct BoundingBox
{
  // The central point of the bounding box (size 2 or 3).
  mfem::Vector center;

  // Vectors from center to the midpoint of each face, stored as columns of a dense matrix.
  // In 3D this is 3x3, in 2D this is 2x2. Column i is the i-th axis vector.
  mfem::DenseMatrix axes;

  // Whether or not this bounding box is two dimensional (i.e. planar in the highest
  // dimension). In 2D meshes this is always true. In 3D meshes this indicates a planar
  // surface.
  bool planar;

  // Return the spatial dimension (number of axes).
  int Dim() const { return axes.Width(); }

  // Compute the lengths along each axis.
  mfem::Vector Lengths() const;
};

// Given a mesh and a marker, compute the bounding circle/sphere of the marked elements. In
// this case the axes of the bounding box object are arbitrary, but the Lengths function
// returns the ball diameter. This function implements Welzl's algorithm.
BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr);

inline BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetBoundingBall(mesh, marker, bdr);
}

// Helper functions for the geometry of the marked boundary elements, sampled at their
// vertices, or on a curved mesh at the points of 4 * order subdivisions of each edge.
// Directions need not be normalized.

// Compute the interval spanned along a direction: the smallest and largest values of
// x · dir / |dir| over the points of the elements.
std::pair<double, double> GetProjectedInterval(const mfem::ParMesh &mesh,
                                               const mfem::Array<int> &marker,
                                               const mfem::Vector &dir);

// Compute the cross-section by the plane (a line in 2D) of points x with
// x · dir / |dir| = offset. Returns its measure: the length of the cut through the surface
// in 3D, or the number of points where the boundary crosses the line in 2D. The centroid
// of the cross-section is returned in centroid, and is zero when the cross-section is
// empty. Points with x · dir / |dir| >= offset are above the plane, and points within
// round-off of the plane are on it, so a plane through an edge shared by two elements
// cuts only one of them there.
double GetCrossSection(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                       const mfem::Vector &dir, double offset, mfem::Vector &centroid);

// Compute the area (length in 2D), centroid, and covariance (the second moments about the
// centroid, per unit area) of the marked boundary elements, integrated on their geometry.
double GetSurfaceMoments(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                         mfem::Vector &centroid, mfem::DenseMatrix &covariance);

// Helper function for computing the closest distance of a marked group to a given point,
// by brute force searching over the entire point set. Optionally compute the furthest
// distance instead of the closest.
double GetDistanceFromPoint(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr, const mfem::Vector &origin, bool max = false);

inline double GetDistanceFromPoint(const mfem::ParMesh &mesh, int attr, bool bdr,
                                   const mfem::Vector &dir, bool max = false)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetDistanceFromPoint(mesh, marker, bdr, dir, max);
}

// Helper function to compute the average surface normal for all elements with the given
// attributes.
mfem::Vector GetSurfaceNormal(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                              bool average = true);

inline mfem::Vector GetSurfaceNormal(const mfem::ParMesh &mesh, int attr,
                                     bool average = true)
{
  const bool bdr = (mesh.Dimension() == mesh.SpaceDimension());
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetSurfaceNormal(mesh, marker, average);
}

inline mfem::Vector GetSurfaceNormal(const mfem::ParMesh &mesh, bool average = true)
{
  const bool bdr = (mesh.Dimension() == mesh.SpaceDimension());
  const auto &attributes = bdr ? mesh.bdr_attributes : mesh.attributes;
  return GetSurfaceNormal(mesh, AttrToMarker(attributes.Max(), attributes), average);
}

// Serial overload of GetSurfaceNormal for use during pre-partition extraction.
mfem::Vector GetSurfaceNormal(const mfem::Mesh &mesh, const mfem::Array<int> &marker,
                              bool average = true);

inline mfem::Vector GetSurfaceNormal(const mfem::Mesh &mesh, bool average = true)
{
  const bool bdr = (mesh.Dimension() == mesh.SpaceDimension());
  const auto &attributes = bdr ? mesh.bdr_attributes : mesh.attributes;
  return GetSurfaceNormal(mesh, AttrToMarker(attributes.Max(), attributes), average);
}

// Submesh post-extraction helpers. Templated over SubMeshT (mfem::SubMesh for the serial
// BoundaryMode path, mfem::ParSubMesh for the parallel WavePort path). MPI reductions
// degenerate to no-ops on the serial (MPI_COMM_SELF) instantiation.

// Remap domain element attributes to match the parent's neighboring 3D element's
// attribute (so material definitions apply directly to the submesh).
template <class SubMeshT>
void RemapSubMeshAttributes(SubMeshT &submesh);

// Trace each submesh boundary edge back to its parent boundary face and assign that
// face's attribute. For edges shared by multiple parent faces, the non-surface face
// wins. The parallel instantiation resolves cross-rank contributions via MPI_Allgather.
template <class SubMeshT>
void RemapSubMeshBdrAttributes(SubMeshT &submesh, const mfem::Array<int> &surface_attrs);

// Add boundary elements at edges where the selected surface meets a parent boundary
// face whose attribute is in internal_bdr_attrs. CreateFromBoundary only emits boundary
// elements at the geometric boundary of the selected region, so these internal
// intersections have to be added explicitly.
template <class SubMeshT>
void AddSubMeshInternalBoundaryElements(SubMeshT &submesh,
                                        const mfem::Array<int> &surface_attrs,
                                        const std::vector<int> &internal_bdr_attrs);

// Project a planar 2D submesh (3D ambient coords) to true 2D coordinates. Replaces node
// coordinates with their projection onto the tangent plane so SpaceDimension() == 2 and
// all 2D infrastructure works as for a native 2D mesh. Returns the surface normal (3D)
// and fills `centroid`, `e1`, `e2` with the tangent frame used for the projection
// (needed by the caller for projecting additional 3D coordinates into the same 2D
// frame). Serial only — extraction runs before partitioning.
mfem::Vector ProjectSubmeshTo2D(mfem::Mesh &submesh, mfem::Vector &centroid,
                                mfem::Vector &e1, mfem::Vector &e2);

// Copy a serial 2D nonconforming mesh to a new mesh with the same elements, boundary
// elements and hanging vertices, whose nonconforming structure is rebuilt from them (the
// leaf elements become the root elements). MFEM keeps the boundary of a nonconforming mesh
// in its mfem::NCMesh, which is used to print or distribute it, and does not update it for
// boundary elements added to the mfem::Mesh or for their attributes. The nonconforming
// mesh of a SubMesh of a nonconforming mesh has no boundary at all.
std::unique_ptr<mfem::Mesh> RebuildNonconformingMesh(const mfem::Mesh &mesh);

// Project a 3D point to 2D local coordinates using a previously computed tangent frame.
inline mfem::Vector Project3Dto2D(const mfem::Vector &p3d, const mfem::Vector &centroid,
                                  const mfem::Vector &e1, const mfem::Vector &e2)
{
  mfem::Vector p2d(2);
  p2d(0) = (p3d(0) - centroid(0)) * e1(0) + (p3d(1) - centroid(1)) * e1(1) +
           (p3d(2) - centroid(2)) * e1(2);
  p2d(1) = (p3d(0) - centroid(0)) * e2(0) + (p3d(1) - centroid(1)) * e2(1) +
           (p3d(2) - centroid(2)) * e2(2);
  return p2d;
}

// Helper functions to compute the volume or area for all domain or boundary elements with
// the given attributes.
double GetSurfaceArea(const mfem::ParMesh &mesh, const mfem::Array<int> &marker);

inline double GetSurfaceArea(const mfem::ParMesh &mesh, int attr)
{
  mfem::Array<int> marker(mesh.bdr_attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetSurfaceArea(mesh, marker);
}

double GetVolume(const mfem::ParMesh &mesh, const mfem::Array<int> &marker);

inline double GetVolume(const mfem::ParMesh &mesh, int attr)
{
  mfem::Array<int> marker(mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetVolume(mesh, marker);
}

// Distribute a serial mesh from the root processor across all MPI ranks using METIS
// partitioning and the MeshPartitioner-based distribution pipeline. The serial mesh need
// only be valid on the root rank (non-root ranks may pass an empty unique_ptr). The serial
// mesh is consumed (released) during distribution.
std::unique_ptr<mfem::ParMesh> DistributeSerialMesh(MPI_Comm comm,
                                                    std::unique_ptr<mfem::Mesh> &smesh);

// Order-independent (conforming) topological entity counts for a mesh, broken down per
// geometry. Populated only on the root rank when an adapted mesh is written (see
// RebalanceMesh and CompleteMeshEntityCounts).
struct MeshEntityCounts
{
  int dim = 0;
  long long true_vertices = 0;
  long long true_edges = 0;
  // TRUE codim-1 faces by geometry (3D only; empty in 2D). Keys are TRIANGLE, SQUARE.
  std::map<mfem::Geometry::Type, long long> true_faces;
  // Raw cells by geometry (2D: TRIANGLE, SQUARE; 3D: TETRAHEDRON, CUBE, PRISM, PYRAMID).
  std::map<mfem::Geometry::Type, long long> cells;
  std::vector<int> domain_attributes, boundary_attributes;
  bool valid = false;  // set true only when populated (root, save happened)
};

// What RebalanceMesh did, the same on all processes: the imbalance ratio (largest over
// smallest number of elements per process) before and after, and whether a rebalance was
// applied (moved; not every element need change process). A conformal mesh for which no
// partition comes back keeps its distribution: its ratio after is the ratio before.
struct RebalanceResult
{
  double ratio_before = 1.0, ratio_after = 1.0;
  bool moved = false;
  // Why no partition came back, when a conformal mesh keeps its distribution.
  std::string reason;
};

// The line reporting a rebalance which the tolerance `tol` called for, without its
// indentation: what moved and the new ratio, or that the distribution was kept.
std::string RebalanceMessage(const RebalanceResult &result, double tol);

// Helper function responsible for rebalancing the mesh, and optionally writing meshes from
// the intermediate stages to disk. Returns what was done (see RebalanceResult). When
// `out_counts` is non-null and the adapted mesh is written, geometry-resolved counts are
// filled on the root rank. Nonconforming vertex and edge counts must then be completed
// collectively with CompleteMeshEntityCounts.
RebalanceResult RebalanceMesh(const IoData &iodata, std::unique_ptr<mfem::ParMesh> &mesh,
                              MeshEntityCounts *out_counts = nullptr);

// Complete true vertex and edge counts for a nonconforming mesh using distributed
// lowest-order H1 and ND spaces. This is collective over the mesh communicator and must be
// called on every rank after RebalanceMesh has populated out_counts for the final mesh.
void CompleteMeshEntityCounts(mfem::ParMesh &mesh, MeshEntityCounts &counts);

// Helper for creating a hexahedral mesh from a tetrahedral mesh.
mfem::Mesh MeshTetToHex(const mfem::Mesh &orig_mesh);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
