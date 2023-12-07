// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include <mfem.hpp>

// XX TODO CONST CORRECTNESS FOR MESH IF GETELEMENTTRANSFORMATION IS CONST AND THREAD SAFE
//         VERSION IS USED

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
std::unique_ptr<mfem::ParMesh> ReadMesh(MPI_Comm comm, const IoData &iodata, bool reorder,
                                        bool clean, bool add_bdr, bool unassembled);

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

// Helper function to convert a set of attribute numbers to a marker array. The marker array
// will be of size max_attr and it will contain only zeroes and ones. Ones indicate which
// attribute numbers are present in the list array. In the special case when list has a
// single entry equal to -1 the marker array will contain all ones.
template <typename T>
void AttrToMarker(int max_attr, const T &attr_list, mfem::Array<int> &marker);
template <typename T>
mfem::Array<int> AttrToMarker(int max_attr, const T &attr_list)
{
  mfem::Array<int> marker;
  AttrToMarker(max_attr, attr_list, marker);
  return marker;
}

// Helper function to construct the bounding box for all elements with the given attribute.
void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                               bool bdr, mfem::Vector &min, mfem::Vector &max);
void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                               mfem::Vector &max);

// Struct describing a bounding box in terms of the center and face normals. The normals
// specify the direction from the center of the box.
struct BoundingBox
{
  // The central point of the bounding box.
  std::array<double, 3> center;

  // Vectors from center to the midpoint of each face.
  std::array<std::array<double, 3>, 3> normals;

  // Whether or not this bounding box is two dimensional.
  bool planar;

  // Compute the area of the bounding box spanned by the first two normals.
  double Area() const;

  // Compute the volume of a 3D bounding box. Returns zero if planar.
  double Volume() const;

  // Compute the lengths of each axis.
  std::array<double, 3> Lengths() const;

  // Compute the deviation in degrees of a vector from each of the normal directions.
  std::array<double, 3> Deviation(const std::array<double, 3> &direction) const;
};

// Struct describing a bounding ball in terms of a center and radius. If a ball is two
// dimensional, additionally provides a normal to the plane.
struct BoundingBall
{
  // The centroid of the ball.
  std::array<double, 3> center;

  // The radius of the ball from the center.
  double radius;

  // If the ball is two dimensional, the normal defining the planar surface. Zero magnitude
  // if a sphere.
  std::array<double, 3> planar_normal;

  // Whether or not this bounding ball is two dimensional.
  bool planar;

  // Compute the area of the bounding box spanned by the first two normals.
  double Area() const { return M_PI * std::pow(radius, 2.0); }

  // Compute the volume of a 3D bounding box. Returns zero if planar.
  double Volume() const { return planar ? 0.0 : (4 * M_PI / 3) * std::pow(radius, 3.0); }
};

// Helper functions for computing bounding boxes from a mesh and markers.
BoundingBox GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr);
BoundingBox GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr);

// Helper function for computing the direction aligned length of a marked group.
double GetProjectedLength(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                          const std::array<double, 3> &dir);
double GetProjectedLength(mfem::ParMesh &mesh, int attr, bool bdr,
                          const std::array<double, 3> &dir);

// Given a mesh and a marker, compute the diameter of a bounding circle/sphere, assuming
// that the extrema points are in the marked group.
BoundingBall GetBoundingBall(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr);
BoundingBall GetBoundingBall(mfem::ParMesh &mesh, int attr, bool bdr);

// Helper function to compute the average surface normal for all elements with the given
// attribute.
mfem::Vector GetSurfaceNormal(mfem::ParMesh &mesh, const mfem::Array<int> &attr_list,
                              bool average = true);
mfem::Vector GetSurfaceNormal(mfem::ParMesh &mesh, int attr, bool average = true);
mfem::Vector GetSurfaceNormal(mfem::ParMesh &mesh, bool average = true);

// Helper function responsible for rebalancing the mesh, and optionally writing meshes from
// the intermediate stages to disk. Returns the imbalance ratio before rebalancing.
double RebalanceMesh(const IoData &iodata, double tol, mfem::ParMesh &mesh);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
