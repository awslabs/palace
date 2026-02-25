// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;

namespace mesh
{

//
// Functions for mesh related functionality.
//

// Read and partition a serial mesh from file, returning a pointer to the new parallel mesh
// object, which should be destroyed by the user.
std::unique_ptr<mfem::ParMesh> ReadMesh(IoData &iodata, MPI_Comm comm);

// Refine the provided mesh according to the data in the input file. If levels of refinement
// are requested, the refined meshes are stored in order of increased refinement. Ownership
// of the initial coarse mesh is inherited by the fine meshes and it should not be deleted.
// The fine mesh hierarchy is owned by the user.
void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

// Dimensionalize a mesh for use in exporting a mesh. Scales vertices and nodes by L.
void DimensionalizeMesh(mfem::Mesh &mesh, double L);

// Nondimensionalize a mesh for use in the solver. Scales vertices and nodes by 1/L.
void NondimensionalizeMesh(mfem::Mesh &mesh, double L);

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

// Struct describing a bounding box in terms of the center and face normals. The normals
// specify the direction from the center of the box. Supports both 2D and 3D: in 2D,
// center has 2 entries and axes is 2x2; in 3D, center has 3 entries and axes is 3x3.
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

  // Compute the area of the bounding box spanned by the first two normals.
  double Area() const;

  // Compute the volume of the 3D bounding box. Returns zero if planar or 2D.
  double Volume() const;

  // Compute the normalized axes of the bounding box, returned as columns of a dense
  // matrix.
  mfem::DenseMatrix Normals() const;

  // Compute the lengths along each axis.
  mfem::Vector Lengths() const;

  // Compute the deviations in degrees of a vector from each of the axis directions. Angles
  // are returned in the interval [0, 180].
  mfem::Vector Deviations(const mfem::Vector &direction) const;
};

// Helper functions for computing bounding boxes from a mesh and markers. These do not need
// to be axis-aligned. Note: This function only returns a minimum oriented bounding box for
// points whose convex hull exactly forms a rectangle or rectangular prism, implementing a
// vastly simplified version of QuickHull for this case. For other shapes, the result is
// less predictable, and may not even form a bounding box of the sampled point cloud.
BoundingBox GetBoundingBox(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                           bool bdr);

inline BoundingBox GetBoundingBox(const mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetBoundingBox(mesh, marker, bdr);
}

// Given a mesh and a marker, compute the bounding circle/sphere of the marked elements. In
// this case the normals of the bounding box object are arbitrary, and the Area and Volume
// members should not be used, but the Lengths function returns the ball diameter. This
// function implements Welzl's algorithm.
BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr);

inline BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetBoundingBall(mesh, marker, bdr);
}

// Helper function for computing the direction aligned length of a marked group.
double GetProjectedLength(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                          bool bdr, const mfem::Vector &dir);

inline double GetProjectedLength(const mfem::ParMesh &mesh, int attr, bool bdr,
                                 const mfem::Vector &dir)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetProjectedLength(mesh, marker, bdr, dir);
}

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

// Remap domain element attributes of a boundary ParSubMesh from parent boundary face
// attributes to the neighboring domain element attributes in the parent mesh. After this
// call, each submesh element has the attribute of its adjacent domain element in the parent
// mesh, matching the material definitions in the config. This enables creating a
// MaterialOperator directly on the submesh.
void RemapSubMeshAttributes(mfem::ParSubMesh &submesh);

// Remap boundary element attributes of a boundary ParSubMesh. By default, MFEM assigns
// all submesh boundary elements the same attribute. This function traces each submesh
// boundary edge back to the parent 3D mesh to find which parent boundary face contains it,
// and assigns that face's attribute. After this call, the submesh boundary attributes match
// the parent mesh's boundary conditions (PEC, absorbing, wave port, etc.).
void RemapSubMeshBdrAttributes(mfem::ParSubMesh &submesh,
                               const mfem::Array<int> &surface_attrs);

// Project a planar 2D submesh (with 3D ambient coordinates from ParSubMesh) to true 2D
// coordinates. Computes the surface normal and tangent frame from the mesh, then replaces
// each node coordinate with its projection onto the tangent plane. After this call, the
// mesh has SpaceDimension() == 2 and all downstream 2D infrastructure (FE spaces, GSLIB)
// works as for a native 2D mesh. Returns the surface normal (3D) for use in material
// tensor projection (BoundaryModeOperator). The optional centroid and tangent vectors
// (e1, e2) are output parameters for transforming additional 3D coordinates (e.g.,
// voltage/current path points) to the same 2D frame.
mfem::Vector ProjectSubmeshTo2D(mfem::ParMesh &submesh, mfem::Vector *centroid = nullptr,
                                mfem::Vector *e1 = nullptr, mfem::Vector *e2 = nullptr);

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

// Helper function responsible for rebalancing the mesh, and optionally writing meshes from
// the intermediate stages to disk. Returns the imbalance ratio before rebalancing.
double RebalanceMesh(const IoData &iodata, std::unique_ptr<mfem::ParMesh> &mesh);

// Helper for creating a hexahedral mesh from a tetrahedral mesh.
mfem::Mesh MeshTetToHex(const mfem::Mesh &orig_mesh);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
