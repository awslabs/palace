// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_IMPL_HPP
#define PALACE_UTILS_GEODATA_IMPL_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "utils/iodata.hpp"

namespace palace
{

struct BoundingBox;

namespace mesh
{

//
// Implementations Functions for mesh related functionality. Isolated to avoid Eigen
// propagation.
//

// For the public interface, we can use a BoundingBox as a generalization of a BoundingBall.
// Internally, however, it's nice to work with a specific ball data type.
struct BoundingBall
{
  Eigen::Vector3d origin;
  double radius;
  bool planar;
};

// Index along a Hilbert curve of a point with dim coordinates of b bits each (Skilling,
// "Programming the Hilbert curve", AIP Conf. Proc. 707, 2004).
std::uint64_t HilbertIndex(std::array<std::uint32_t, 3> X, int dim, int b);

// Rank of each part of a partition with one part per rank, which groups the parts sharing
// faces on the same node. pairs lists (p, q, w) for parts p != q sharing faces of weight
// w > 0, and rank_node the node of each rank, as a label in [0, size) such as its lowest
// rank. For each of the orders of the parts, the parts are assigned to the nodes in
// contiguous blocks of the order by ascending label, swapped in pairs between nodes while
// this reduces the weight between nodes, and given the ranks of their node in the order.
// Returns the result with the least weight between nodes, the first on ties.
std::vector<int> NumberPartsOnNodes(const std::vector<std::vector<int>> &orders,
                                    const std::vector<std::array<int, 3>> &pairs,
                                    const std::vector<int> &rank_node);

// Helper for collecting a point cloud from a mesh, used in calculating bounding boxes and
// bounding balls. Returns the dominant rank, for which the vertices argument will be
// filled, while all other ranks will have an empty vector. Vertices are de-duplicated to a
// certain floating point precision.
int CollectPointCloudOnRoot(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr, std::vector<Eigen::Vector3d> &vertices);

// Calculates a bounding box from a point cloud, result is broadcast across all processes.
BoundingBox BoundingBoxFromPointCloud(MPI_Comm comm,
                                      const std::vector<Eigen::Vector3d> &vertices,
                                      int dominant_rank);

// Compute the distance from a point orthogonal to the list of normal axes, relative to
// the given origin.
auto PerpendicularDistance(const std::initializer_list<Eigen::Vector3d> &normals,
                           const Eigen::Vector3d &origin, const Eigen::Vector3d &v);

// Use 4 points to define a sphere in 3D. If the points are coplanar, 3 of them are used to
// define a circle which is interpreted as the equator of the sphere. We assume the points
// are unique and not collinear.
BoundingBall SphereFromPoints(const std::vector<std::size_t> &indices,
                              const std::vector<Eigen::Vector3d> &vertices);

// Implementation of the recursive Welzl algorithm kernel.
BoundingBall Welzl(std::vector<std::size_t> P, std::vector<std::size_t> R,
                   const std::vector<Eigen::Vector3d> &vertices);

// Calculates a bounding ball from a point cloud using Welzl's algorithm, result is
// broadcast across all processes. We don't operate on the convex hull, since the number of
// points should be small enough that operating on the full set should be OK. If only three
// points are provided, the bounding circle is computed (likewise for if the points are
// coplanar).
BoundingBox BoundingBallFromPointCloud(MPI_Comm comm,
                                       const std::vector<Eigen::Vector3d> &vertices,
                                       int dominant_rank);

// Compute a normal vector from an element transformation, optionally ensure aligned
// (| normal ⋅ align | > 0)
void Normal(mfem::ElementTransformation &T, mfem::Vector &normal,
            const mfem::Vector *const align);

// Determine the vertex mapping between donor and receiver boundary attributes.
// Uses the translation vector or affine transformation matrix specified in the
// configuration file. If not provided, attempts to automatically detect the
// affine transformation between donor and receiver boundary vertices.
std::vector<int>
DeterminePeriodicVertexMapping(std::unique_ptr<mfem::Mesh> &mesh,
                               const struct palace::config::PeriodicData &data,
                               const double tol = 1e-8);

// Conformally bisect the selected mesh edges by splitting each incident tetrahedron and
// boundary triangle. Used by periodic crack refinement.
int LocalEdgeSplit(std::unique_ptr<mfem::Mesh> &mesh,
                   const std::vector<std::pair<int, int>> &edges);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_IMPL_HPP
