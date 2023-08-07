// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <array>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <mpi.h>

namespace mfem
{

template <typename T>
class Array;
class ParMesh;
class Vector;

}  // namespace mfem

namespace palace
{

class IoData;
class Timer;

namespace mesh
{

//
// Functions for mesh related functionality.
//

// Read and partition a serial mesh from file, returning a pointer to the new parallel mesh
// object, which should be destroyed by the user.
std::unique_ptr<mfem::ParMesh> ReadMesh(MPI_Comm comm, const IoData &iodata, bool reorder,
                                        bool clean, bool add_bdr, bool unassembled,
                                        Timer &timer);

// Refine the provided mesh according to the data in the input file. If levels of refinement
// are requested, the refined meshes are stored in order of increased refinement. Ownership
// of the initial coarse mesh is inherited by the fine meshes and it should not be deleted.
// The fine mesh hierarchy is owned by the user.
void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

// Helper function to convert a set of attribute numbers to a marker array. The marker array
// will be of size max_attr and it will contain only zeroes and ones. Ones indicate which
// attribute numbers are present in the attrs array. In the special case when attrs has a
// single entry equal to -1 the marker array will contain all ones.
void AttrToMarker(int max_attr, const mfem::Array<int> &attrs, mfem::Array<int> &marker);
void AttrToMarker(int max_attr, const std::vector<int> &attrs, mfem::Array<int> &marker);

// Helper function to construct the bounding box for all elements with the given attribute.
void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                               mfem::Vector &max);
void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                               bool bdr, mfem::Vector &min, mfem::Vector &max);

// Struct describing a bounding box in terms of the center and face normals. The
// normals specify the direction from the center of the box.
struct BoundingBox
{
  // The central point of the bounding box.
  Eigen::Vector3d center;
  // Vectors from center to the midpoint of each face.
  std::array<Eigen::Vector3d, 3> normals;
  // Whether or not this bounding box is two dimensional.
  bool planar;
  // Compute the area of the bounding box spanned by the first two normals.
  double Area() const;
  // Compute the volume of a 3D bounding box. Returns zero if planar.
  double Volume() const;
};

// Computes a bounding box for a point cloud, assuming a distributed point cloud
// of vertices.
BoundingBox BoundingBoxFromPointCloud(MPI_Comm comm, std::vector<Eigen::Vector3d> vertices);
// Helper functions for computing bounding boxes from a mesh and markers.
BoundingBox GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr);
BoundingBox GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr);

// Struct describing a bounding ball in terms of a center and radius. If a ball
// is two dimensional, additionally provides a normal to the plane.
struct BoundingBall
{
  // The centroid of the ball.
  Eigen::Vector3d center;
  // The radius of the ball from the center.
  double radius;
  // If the ball is two dimensional, the normal defining the planar surface.
  // Zero magnitude if a sphere.
  Eigen::Vector3d planar_normal;
};

// Helper for computing a diameter from a point cloud. The cloud is assumed to
// contain a set of points on a circle.
BoundingBall BoundingBallFromPointCloud(MPI_Comm comm,
                                        std::vector<Eigen::Vector3d> vertices);
// Given a mesh and a marker, compute the diameter of a bounding circle/sphere,
// assuming that the extrema points are in the marked group.
BoundingBall GetBoundingBall(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr);
BoundingBall GetBoundingBall(mfem::ParMesh &mesh, int attr, bool bdr);

// Helper function to compute the average surface normal for all elements with the given
// attribute.
void GetSurfaceNormal(mfem::ParMesh &mesh, int attr, mfem::Vector &normal);
void GetSurfaceNormal(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                      mfem::Vector &normal);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
