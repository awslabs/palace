// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_NCMESH_DIST_HPP
#define PALACE_UTILS_NCMESH_DIST_HPP

#include <memory>
#include <mfem.hpp>

namespace palace::mesh
{

// Distribute a conforming serial mesh as a nonconforming mfem::ParMesh. The serial mesh
// and the element partitioning are held by the first process of each shared-memory node
// (smesh is nullptr elsewhere). The result is identical to every process constructing
// mfem::ParMesh(comm, smesh, partitioning) from its own copy of the serial mesh after
// EnsureNCMesh, but the other processes only receive the element connectivity (to build
// the root-element NCMesh, which mfem::ParNCMesh replicates on every process) and the
// curved mesh nodes of their own elements, instead of the full serial mesh. With
// local_reorder, the root elements are numbered along a Hilbert curve through their
// centroids instead of in the serial order, with the same partition, for the locality of
// the local elements of every process (see SetRootElementOrigin).
std::unique_ptr<mfem::ParMesh> DistributeNCMesh(MPI_Comm comm,
                                                std::unique_ptr<mfem::Mesh> &smesh,
                                                const int *partitioning,
                                                bool local_reorder = false);

// Distribute a nonconforming serial mesh with a refinement tree (for example, the adapted
// mesh saved by a previous adaptive simulation), held by the first process of each
// shared-memory node, as a nonconforming mfem::ParMesh with the given partitioning of its
// elements. Every process builds the NCMesh of the root elements, the refinement is
// replayed level by level in parallel, the curved mesh nodes of each element are copied
// from the serial mesh, and the mesh is rebalanced to the partitioning. The result has the
// same elements, nodes, attributes and partition as mfem::ParMesh(comm, smesh,
// partitioning) on every process, but no process other than the node roots holds the
// serial mesh. Returns nullptr, with smesh untouched, on every process if the refinement
// tree cannot be replayed: anisotropic refinement, root elements with hanging vertices, or
// boundary elements which do not lie on the faces of the root elements.
std::unique_ptr<mfem::ParMesh> DistributeNCMeshTree(MPI_Comm comm,
                                                    std::unique_ptr<mfem::Mesh> &smesh,
                                                    const int *partitioning);

}  // namespace palace::mesh

#endif  // PALACE_UTILS_NCMESH_DIST_HPP
