// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GEODATA_HPP
#define PALACE_UTILS_GEODATA_HPP

#include <memory>
#include <vector>
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
void GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                    mfem::Vector &max);
void GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                    mfem::Vector &min, mfem::Vector &max);

// Helper function to compute the average surface normal for all elements with the given
// attribute.
void GetSurfaceNormal(mfem::ParMesh &mesh, int attr, mfem::Vector &normal);
void GetSurfaceNormal(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                      mfem::Vector &normal);

// Rebalance a conformal mesh across processor ranks, using the MeshPartitioner.
// Gathers the mesh onto the root rank before scattering the partitioned mesh.
void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &mesh);

}  // namespace mesh

}  // namespace palace

#endif  // PALACE_UTILS_GEODATA_HPP
