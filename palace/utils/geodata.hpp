// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <array>
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
std::unique_ptr<mfem::ParMesh> ReadMesh(const IoData &iodata, MPI_Comm comm);

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
  std::vector<mfem::Geometry::Type> GetGeomTypes() const;
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
// specify the direction from the center of the box.
struct BoundingBox
{
  // The central point of the bounding box.
  std::array<double, 3> center;

  // Vectors from center to the midpoint of each face.
  std::array<std::array<double, 3>, 3> axes;

  // Whether or not this bounding box is two dimensional.
  bool planar;

  // Compute the area of the bounding box spanned by the first two normals.
  double Area() const;

  // Compute the volume of the 3D bounding box. Returns zero if planar.
  double Volume() const;

  // Compute the normalized axes of the bounding box.
  std::array<std::array<double, 3>, 3> Normals() const;

  // Compute the lengths along each axis.
  std::array<double, 3> Lengths() const;

  // Compute the deviations in degrees of a vector from each of the axis directions. Angles
  // are returned in the interval [0, 180].
  std::array<double, 3> Deviations(const std::array<double, 3> &direction) const;
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
                          bool bdr, const std::array<double, 3> &dir);

inline double GetProjectedLength(const mfem::ParMesh &mesh, int attr, bool bdr,
                                 const std::array<double, 3> &dir)
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
                            bool bdr, const std::array<double, 3> &origin,
                            bool max = false);

inline double GetDistanceFromPoint(const mfem::ParMesh &mesh, int attr, bool bdr,
                                   const std::array<double, 3> &dir, bool max = false)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetDistanceFromPoint(mesh, marker, bdr, dir);
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

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
